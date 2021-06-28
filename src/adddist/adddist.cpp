/***************************************************************************
 *   Copyright (C) 2013-2020 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#include "liblib/mybase.h"

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>

#include <vector>
#include <numeric>
#include <algorithm>

#include "liblib/mygetopt.h"
#include "liblib/msg.h"
#include "incconf/localconfig.h"
#include "libpro/srcpro/MOptions.h"
#include "libpro/srcpro/datapro.h"
#include "libpro/srcpro/SEGProfile.h"
#include "libpro/srcpro/PMTransModel.h"
#include "libpro/srcpro/PMProfileModel.h"
#include "libHDP/HDPbase.h"
#include "libpro/srcaln/MSA.h"
#include "adddist.h"


// =========================================================================
// declarations of functions:
void GetDstColumns(mystring, std::vector<int>&);
void ReadDistances(const char* filename, const int prolen, 
            const std::vector<int>& colvec, float prbthrhld, 
            std::vector<std::vector<float>>& distogram);
void GetPositionPairHelper(int,  const mystring&, const char*&, int&, 
        const std::vector<int>&, int[2]);
void GetPairDistanceHelper(int, const mystring&, const char*&, int&, 
        const std::vector<int>&, float, float*);
void ProcessDistancesForPosition(int position, 
        std::vector<float>&, std::vector<float>&,
        std::vector<int>&, std::vector<int>&,
        std::vector<std::vector<std::pair<float,int>>>&,
        std::vector<std::vector<float>>&);
// =========================================================================


int main( int argc, char *argv[] )
{
    int             c;
    char*           p;
    float           fval;
    //string values of options
    mystring        myoptarg;
    mystring        dstfile;
    mystring        profile;
    mystring        output;
    mystring        optfile;
    mystring        dstcolist;//list of distance columns
    mystring        prbthshld;//probability threshold
    std::vector<int> vdcolumns(1, 3);//vector of unique values of distance columns (def=3)
    float           prbthshldval = 0.0f;
    bool            suppress = true;//suppress warnings

    SetArguments( &argc, &argv );
    SetProgramName( argv[0], version );

    if( argc <= 1 ) {
        fprintf( stdout, "%s", usage(argv[0],makeinst,version,verdate).c_str());
        return EXIT_SUCCESS;
    }

    static struct myoption long_options[] = {
        {"i", my_required_argument, 'i'},
        {"j", my_required_argument, 'j'},
        {"o", my_required_argument, 'o'},
        {"p", my_required_argument, 'p'},
        {"v", my_no_argument,       'v'},
        {"h", my_no_argument,       'h'},
        {"dst", my_required_argument, 'L'},
        {"prb", my_required_argument, 'P'},
        { NULL, my_n_targflags, 0 }
    };

    try {
        try {
            MyGetopt mygetopt( long_options, (const char**)argv, argc );
            while(( c = mygetopt.GetNextOption( &myoptarg )) >= 0 ) {
                switch( c ) {
                    case ':':   fprintf( stdout, "Argument missing. Please try option -h for help.%s", NL );
                                return EXIT_FAILURE;
                    case '?':   fprintf( stdout, "Unrecognized option. Please try option -h for help.%s", NL );
                                return EXIT_FAILURE;
                    case '!':   fprintf( stdout, "Ill-formed option. Please try option -h for help.%s", NL );
                                return EXIT_FAILURE;
                    case 'h':   fprintf( stdout, "%s", usage(argv[0],makeinst,version,verdate).c_str());
                                return EXIT_SUCCESS;
                    case 'i':   dstfile = myoptarg; break;
                    case 'j':   profile = myoptarg; break;
                    case 'o':   output = myoptarg; break;
                    case 'p':   optfile = myoptarg; break;
                    case 'v':   suppress = false; break;
                    //
                    case 'L':   dstcolist = myoptarg; break;
                    case 'P':   prbthshld = myoptarg; break;
                    default:    break;
                }
            }
        } catch( myexception const& ex ) {
            error( dynamic_cast<myruntime_error const&>(ex).pretty_format().c_str());
            return EXIT_FAILURE;
        }
    } catch ( ... ) {
        error("Unknown exception caught.");
        return EXIT_FAILURE;
    }

    SetVerboseMode( !suppress );

    if( dstfile.empty()) {
        error( "File of distances is not specified." );
        return EXIT_FAILURE;
    }
    if( profile.empty()) {
        error( "Input profile is not specified." );
        return EXIT_FAILURE;
    }
    if( output.empty()) {
        error( "Output file is not specified." );
        return EXIT_FAILURE;
    }


    TRY
        if( !dstcolist.empty())
            GetDstColumns(dstcolist, vdcolumns);
        if(vdcolumns.size() < 1) {
            error( "No columns of distances given." );
            return EXIT_FAILURE;
        }
        //fprintf(stderr,"Columns:\n");
        //std::for_each(vdcolumns.begin(), vdcolumns.end(), [](const int& n){fprintf(stderr,"%d\n",n);});
    CATCH_ERROR_RETURN(;);


    if( !prbthshld.empty()) {
        fval = strtof( prbthshld.c_str(), &p );
        if( errno || *p ) {
            error( "Probability threshold is invalid." );
            return EXIT_FAILURE;
        }
        if( fval < 0.0f || 1.0f <= fval ) {
            error( "Probability threshold is outside the interval [0,1)." );
            return EXIT_FAILURE;
        }
        prbthshldval = fval;
    }


    mystring        insoptfile = GetFullOptionsFilename();
    mystring        altinsoptfile =
                        mystring( my_dirname( argv[0])) + DIRSEP +
                        UPDIR + DIRSEP +
                        GetParamDirectory() + DIRSEP +
                        GetOptionsFilename();

    if( !file_exists( insoptfile.c_str()))
        insoptfile = altinsoptfile;

    if( optfile.empty())
        optfile = insoptfile;


    TRY
        MOptions::Read( optfile.c_str());
    CATCH_ERROR_RETURN(;);



    int     ret = EXIT_SUCCESS;
    FILE*   fp = NULL;

    TRY
        message( "Reading input profile...");

        pmodel::PMProfileModel  pmpm;
        pmodel::PMTransModel    pmtm;

        if(( fp = fopen(profile.c_str(), "r")) ==  NULL )
            throw MYRUNTIME_ERROR("Main: Failed to open profile.");

        TextReadProfileCOMER(fp, pmpm, pmtm);

        fclose(fp);
        fp = NULL;

        message( "Reading distances...");
        std::vector<std::vector<float>> distogram(pmpm.GetSize());

        ReadDistances(dstfile.c_str(), pmpm.GetSize(), 
            vdcolumns, prbthshldval, distogram);

//         c = 0;
//         for(const std::vector<float>& dn: distogram) {
//             if(c++, dn.size() < 1)
//                 continue;
//             fprintf(stdout,"%3d:",c);
//             std::for_each(dn.begin(), dn.end(), [](const float& d){fprintf(stdout," %.0f",d);});
//             fprintf(stdout,"\n");
//         }

        pmpm.SetDistogram(std::move(distogram));

        message( "Writing output profile...");

        fp = fopen( output.c_str(), "w" );
        if( fp == NULL )
            throw MYRUNTIME_ERROR("Failed to open file for writing.");

        TextWriteProfileCondensed(fp, pmpm, pmtm);

        message( "Done.");

    CATCH_ERROR_RETURN(if(fp) fclose(fp));

    return ret;
}

// -------------------------------------------------------------------------
// GetDstColumns: get column numbers corresponding to inter-residue 
// distances;
// colist, string of column numbers separated by commas;
// colvec, sorted vector of unique column numbers
//
void GetDstColumns(mystring colist, std::vector<int>& colvec)
{
    MYMSG( "Main::GetDstColumns", 4 );
    char* p;
    size_t pos;
    colvec.clear();
    for( pos = colist.find(','); ; pos = colist.find(',')) {
        mystring strval = colist.substr( 0, pos );
        errno = 0;
        int c = (int)strtol( strval.c_str(), &p, 10 );
        if( errno || *p || strval.empty())
            throw MYRUNTIME_ERROR("Invalid column specified by option dst.");
        if( c < 3 || 100 < c )
            throw MYRUNTIME_ERROR(
            "A column specified by command-line option dst is outside the interval [3,100].");
        colvec.push_back(c-1);//insert zero-based column indices
        if( pos == mystring::npos )
            break;
        colist = colist.substr(pos+1);
    }
    std::sort(colvec.begin(), colvec.end());
    auto last = std::unique(colvec.begin(), colvec.end());
    colvec.erase(last, colvec.end());
}

// -------------------------------------------------------------------------
// ReadDistances: read distance information from file;
// filename, name of file of distances (and optionally probabilities);
// prolen, profile length;
// colvec, distance columns to consider;
// prbthrhld, distance probability threshold;
// distogram, output distogram
//
void ReadDistances(const char* filename, const int prolen, 
            const std::vector<int>& colvec, float prbthrhld, 
            std::vector<std::vector<float>>& distogram)
{
    MYMSG( "Main::ReadDistances", 4 );
    mystring preamb = "Main::ReadDistances: ";
    myruntime_error mre;
    const int optdstnmax = MAX_N_IR_DISTANCE_VALUES;//MOptions::GetDSTNMAX();
    std::vector<float> vtmpdsts;//temporary vector of distances
    std::vector<float> vdsts;//vector of distances
    std::vector<int> vposs;//vector of positions
    std::vector<int> vndxs;//vector of indices
    char msgbuf[BUF_MAX];
    mystring buffer;
    const char* p;
    int emsg;
    int rpair[2] = {0,0}, prevpos = 1;//positions in file are 1-based
    //temporary lower part distogram and indices:
    std::vector<std::vector<std::pair<float,int>>> tmpltridgram(prolen);

    if(colvec.size() < 1)
        return;

    if(!filename)
        throw MYRUNTIME_ERROR(preamb + "Empty filename.");

    FILE* fp = fopen( filename, "r" );

    if(!fp)
        throw MYRUNTIME_ERROR(preamb + "Failed to open file: " + filename);

    try {
        vtmpdsts.reserve(optdstnmax*2);

        for(int ll = 0; !feof(fp); ll++)
        {
            if(( emsg = skip_comments(fp, buffer)) != 0 )
                throw MYRUNTIME_ERROR(preamb + TranslateReadError(emsg));

            if( buffer.empty())
                continue;

            int c = 0;
            float dst = -1.f;

            p = buffer.c_str();

            GetPositionPairHelper(ll, buffer, p, c, colvec, rpair);
            GetPairDistanceHelper(ll, buffer, p, c, colvec, prbthrhld, &dst);

            if( rpair[0] < prevpos || rpair[1] <= rpair[0]) {
                sprintf(msgbuf, "line %d, comments excluded.", ll);
                throw MYRUNTIME_ERROR(preamb + "Invalid distance file format: "
                "The positions (1-2 columns) should be non-decreasing: " + msgbuf);
            }

            if( prevpos < rpair[0]) {
                //the next position has begun;
                //position under consideration is prevpos
                if( prolen < rpair[0] || prolen < rpair[1]) {
                    sprintf(msgbuf, " %d at line %d (comments excluded).", prolen, ll);
                    throw MYRUNTIME_ERROR(preamb + "Invalid distance file format: "
                    "Position(s) are greater than the profile length " + msgbuf);
                }

                ProcessDistancesForPosition(
                        prevpos-1,//0-based 
                        vtmpdsts,
                        vdsts, vposs, vndxs,
                        tmpltridgram,
                        distogram);

                vdsts.clear();
                vposs.clear();
                vtmpdsts.clear();

                prevpos = rpair[0];
            }

            if( dst >= 0.f ) {
                vdsts.push_back(dst);
                vposs.push_back(rpair[1]-1);//0-based
                //
                tmpltridgram[rpair[1]-1].push_back({dst,rpair[0]-1});//-1 to make 0-based
            }
        }

        //ensure all information is included
        ProcessDistancesForPosition(
                prevpos-1,//0-based 
                vtmpdsts,
                vdsts, vposs, vndxs,
                tmpltridgram,
                distogram);

    } catch( myruntime_error const& ex ) {
        mre = ex;
    } catch( myexception const& ex ) {
        mre = ex;
    } catch( ... ) {
        mre = myruntime_error("Unknown exception caught.");
    }

    fclose(fp);

    if( mre.isset())
        throw mre;
}

// -------------------------------------------------------------------------
// GetPositionPairHelper: get a pair of residue indices from one file line;
// linenumber, file line number;
// buffer, file line read;
// p, current processing position buffer;
// c, starting column number in the read line for processing;
// colvec, distance columns to consider;
// rpair, array of an output pair of residue indices;
//
void GetPositionPairHelper(
        int linenumber, 
        const mystring& buffer, const char*& p, int& c, 
        const std::vector<int>& /*colvec*/, 
        int rpair[2])
{
    MYMSG( "Main::GetPositionPairHelper", 4 );
    mystring preamb = "Main::GetPositionPairHelper: ";
    char msgbuf[BUF_MAX];
    size_t rbts;
    int emsg;
    int intval;

    if( !p)
        return;

    //read a pair of residue indices
    for( ; c < 2; c++ )
    {
        if( buffer.length() <= (size_t)(p-buffer.c_str())) {
            sprintf(msgbuf, "line %d, comments excluded.", linenumber);
            throw MYRUNTIME_ERROR(preamb + "Invalid distance file format: "
            "Missing a pair of residue indices in the first two columns: " + msgbuf);
        }

        if(( emsg = read_integer(p, buffer.length() - (size_t)(p-buffer.c_str()), &intval, &rbts )) != 0) {
            sprintf(msgbuf, "(line %d, comments excluded): ", linenumber);
            throw MYRUNTIME_ERROR(preamb + 
            "Invalid residue indices in the first two columns " + 
            msgbuf + TranslateReadError(emsg));
        }

        p += rbts;

        if(intval < 1 || 10000 < intval) {
            sprintf(msgbuf, "line %d, comments excluded.", linenumber);
            throw MYRUNTIME_ERROR(preamb + "Residue index out of bounds: " + msgbuf);
        }

        rpair[c] = intval;
    }
}

// -------------------------------------------------------------------------
// GetPairDistanceHelper: process distance information from one file line;
// linenumber, file line number;
// buffer, file line read;
// p, current processing position buffer;
// c, starting column number in the read line for processing;
// colvec, distance columns to consider;
// prbthrhld, distance probability threshold;
// distance, output calculated distance
//
void GetPairDistanceHelper(
        int linenumber, 
        const mystring& buffer, const char*& p, int& c, 
        const std::vector<int>& colvec, float prbthrhld, 
        float* distance)
{
    MYMSG( "Main::GetPairDistanceHelper", 4 );
    mystring preamb = "Main::GetPairDistanceHelper: ";
    char msgbuf[BUF_MAX];
    const char* pp;
    float dst, fval;
    int ndx, cve;
    int ndsts = 0;
    float dstsum = 0.0f;
    size_t rbts;
    int emsg;

    int coldstlast = colvec.back();
    int colast = (prbthrhld > 0.0f)? coldstlast+1: coldstlast;

    if( !distance || !p)
        return;

    *distance = -1.f;

    //read distances and, optionally, their probabilities
    for( ndx = 0; c <= colast && ndx < (int)colvec.size(); c++ )
    {
        if( buffer.length() <= (size_t)(p-buffer.c_str())) {
            sprintf(msgbuf, "No.%d: line %d, comments excluded.", ndx+1, linenumber);
            throw MYRUNTIME_ERROR(preamb + "Invalid distance file format: "
            "Missing a value for distance " + msgbuf);
        }

        for( ; *p==' ' || *p=='\t'; p++ );
        for( pp = p; pp < buffer.c_str() + buffer.length() &&
            *pp && *pp!=' ' && *pp!='\t' && *pp!='\n' && *pp!='\r'; pp++ );

        cve = colvec[ndx];

        if( c < cve ) {//c is zero-based index
            p = pp;
            continue;
        }

        ndx++;//move to the next column of interest in advance

        if( *p == '-' ) {
            //ignore fields starting with `-'
            p = pp;
            continue;
        }

        if(( emsg = read_float(p, buffer.length() - (size_t)(p-buffer.c_str()), &fval, &rbts )) != 0) {
            sprintf(msgbuf, "No.%d (line %d, comments excluded): ", ndx, linenumber);
            throw MYRUNTIME_ERROR(preamb + 
            "Invalid value for distance " + msgbuf + TranslateReadError(emsg));
        }

        p += rbts;

        if(fval < 0.f || 64.f < fval) {
            sprintf(msgbuf, "No.%d: line %d, comments excluded.", ndx, linenumber);
            throw MYRUNTIME_ERROR(preamb + "Invalid value for distance " + msgbuf);
        }

        dst = fval;

        if(prbthrhld > 0.0f) {
            //distances only of sufficient confidence are considered;
            if(++c > colast) {
                sprintf(msgbuf, "No.%d: line %d, comments excluded.", ndx, linenumber);
                throw MYRUNTIME_ERROR(preamb + "Invalid distance file format: "
                "Missing a probability value for distance " + msgbuf);
            }

            //read probability for the distance just read
            if(( emsg = read_float(p, buffer.length() - (size_t)(p-buffer.c_str()), &fval, &rbts )) != 0) {
                sprintf(msgbuf, "No.%d (line %d, comments excluded): ", ndx, linenumber);
                throw MYRUNTIME_ERROR(preamb + 
                "Invalid probability value for distance " + msgbuf + TranslateReadError(emsg));
            }

            p += rbts;

            if(fval < 0.f || 1.f < fval) {
                sprintf(msgbuf, "No.%d: line %d, comments excluded.", ndx, linenumber);
                throw MYRUNTIME_ERROR(preamb + "Invalid probability value for distance " + msgbuf);
            }

            if(fval < prbthrhld)
                //distance to be omitted due to low confidence
                continue;
        }

        dstsum += dst;
        ndsts++;
    }

    if(ndx < (int)colvec.size()) {
        sprintf(msgbuf, "%zu (line %d, comments excluded).", colvec.size(), linenumber);
        throw MYRUNTIME_ERROR(preamb + "Missing distance values; expected number: " + msgbuf);
    }

    if(ndsts)
        if(!prbthrhld || (prbthrhld && ndsts==(int)colvec.size()))
            *distance = dstsum / ndsts;
}

// -------------------------------------------------------------------------
// ProcessDistancesForPosition: process distance information for one profile 
// position and add compiled distances at the corresponding position of the 
// distogram;
// position, profile position under consideration;
// vtmpdsts, temporary vector of distances [passed to avoid allocation];
// vdsts, vector of distances [passed to avoid allocation];
// vposs, vector of positions [passed to avoid allocation];
// vndxs, vector of indices [passed to avoid allocation];
// tmpltridgram, temporary lower part distogram and indices, which is 
// contantly updated;
// distogram, output distogram updated at <position>;
//
void ProcessDistancesForPosition(
        int position, 
        std::vector<float>& vtmpdsts,
        std::vector<float>& vdsts,
        std::vector<int>& vposs,
        std::vector<int>& vndxs,
        std::vector<std::vector<std::pair<float,int>>>& tmpltridgram,
        std::vector<std::vector<float>>& distogram)
{
    MYMSG( "Main::GetPairDistanceHelper", 4 );
    mystring preamb = "Main::GetPairDistanceHelper: ";
    const static int optdstnmax = MAX_N_IR_DISTANCE_VALUES; //MOptions::GetDSTNMAX();
    const static int optdstsegm = MOptions::GetDSTSEGM();
    const static int optdstnpos = MOptions::GetDSTNPOS();
    const static float optdstsep6 = MOptions::GetDSTSEP6();
    const static float optdstgapl = MOptions::GetDSTGAPL();
    const static int optdstfree = MOptions::GetDSTFREE();
    float dstpen;//distance penalty value

    //copy/move all (column) elements saved up to this position so far;
    //this columnn (and row, major position) will be done
    for(const std::pair<float,int>& dn: tmpltridgram[position]) {
        vdsts.push_back(dn.first);
        vposs.push_back(dn.second);
    }

    //initialize indices
    vndxs.resize(vdsts.size());
    std::iota(vndxs.begin(), vndxs.end(), 0);

    //sort by how many positions the distances are apart from this 
    // position under consideration
    std::sort(vndxs.begin(), vndxs.end(),
        [position, &vposs](size_t n1, size_t n2) {
            return abs(position-vposs[n1]) < abs(position-vposs[n2]);
        });

    //take only a limited number of values
    if(optdstnmax < (int)vndxs.size())
        vndxs.resize(optdstnmax);

    //sort selected distances by position
    std::sort(vndxs.begin(), vndxs.end(),
        [&vposs](size_t n1, size_t n2) {return vposs[n1] < vposs[n2];}
    );

    //fill in gaps between distance values with penalty distances in the 
    // sorted vector of distances vdsts
    if(vndxs.size() && 
       optdstfree < abs(position - vposs[vndxs[0]]))
        vtmpdsts.push_back(vdsts[vndxs[0]]);
    int maxmainsep = 1;
    for(size_t i = 1; i < vndxs.size(); i++) {
        int posep = vposs[vndxs[i]] - vposs[vndxs[i-1]];
        int mainsepp = abs(position - vposs[vndxs[i-1]]);
        int mainsep = abs(position - vposs[vndxs[i]]);
        //get the farthest distance in sequence from position
        int mainsp2 = SLC_MAX(mainsepp, mainsep);
        if(mainsep <= optdstfree) {
            //do not penalize positions falling into the band of the main diagonal;
            if(posep > 1 && maxmainsep < mainsp2)// && !(posep == 2 && vposs[vndxs[i]] - position == 1))
                maxmainsep = mainsp2;
            continue;
        }
        //if( posep >= optdstnpos || maxposep > optdstnpos)//prev.
        //if there's been a gap between the positions including the 
        // region which is ignored, then include penalties:
        if((posep > 1 && mainsp2 >= optdstnpos) || maxmainsep >= optdstnpos)
            vtmpdsts.push_back(dstpen = optdstgapl);
        else if(posep > 1 || maxmainsep > 1)
            vtmpdsts.push_back(dstpen = optdstsep6);
        //NOTE: make sure that the number of missing values greater than DSTSEGM
        // are balanced with distance penalties of the corresponding length
        for(int k = 2; k < posep && k < optdstsegm; k++)
            vtmpdsts.push_back(dstpen);
        vtmpdsts.push_back(vdsts[vndxs[i]]);
        maxmainsep = 1;//reset max value, it won't be needed
    }

    //take only a limited number of values
    if(optdstnmax < (int)vtmpdsts.size()) {
        if(position + optdstnmax > (int)distogram.size())
            //abandon the first entries if position is no more than half optdstnmax 
            // positions apart from the end
            vtmpdsts.erase(vtmpdsts.begin(), vtmpdsts.begin() + vtmpdsts.size()-optdstnmax);
        else
            vtmpdsts.resize(optdstnmax);
    }

    //assign compiled distances to the distogram's corresponding position
    distogram[position] = vtmpdsts;
}
