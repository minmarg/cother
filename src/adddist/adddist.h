/***************************************************************************
 *   Copyright (C) 2013-2020 by Mindaugas Margelevicius                    *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __adddist_h__
#define __adddist_h__

// Version history:
// 1.01   Initial version

static const char*  version = "1.01";
static const char*  verdate = "";

static const char*  makeinst = "\n\
<>\n\
\n\
Add pairwise distance information to corresponding positions of the profile.\n\
(C)2013-2020 Mindaugas Margelevicius, Institute of Biotechnology, Vilnius University\n\
\n\
\n\
Usage:\n\
<> -i <distances> -j <profile> -o <output> [-p <options>] [<other_options>]\n\
\n\
Profile construction options:\n\
\n\
-i <distances_file>         Input file of distances, where the first two \n\
                            columns give residue (one-based) indices followed\n\
                            by distance information (distance(s), optional\n\
                            probabilities). The file is required to provide\n\
                            the ordered upper triangle of the distance matrix.\n\
-j <profile_file>           Input COMER profile.\n\
-o <output_file>            Output COTHER profile with added distances.\n\
-p <options_file>           Input file of options.\n\
                        By default, the options file in the installation\n\
                            directory is used.\n\
\n\
Distance specification options:\n\
--dst=<column_list>         Comma-separated list of columns representing\n\
                            distances to be averaged.\n\
                        Default=3\n\
--prb=<probability>         Probability lower-bound threshold at which to\n\
                            consider predicted distances. Probabilities are\n\
                            assumed to be given next to distances (i.e.,\n\
                            <column_list>+1 for each entry in the list).\n\
                            If 0, probabilities are ignored and may be absent.\n\
                        Default=0\n\
\n\
Other options:\n\
-v                          Verbose mode.\n\
-h                          This text.\n\
\n\
\n\
Examples:\n\
<> -i mydist.txt -j myprofile.pro -o my_output_cother_profile\n\
<> -i mydist.txt -j myprofile.pro -o my_output_cother_profile --dst=3,5,7 --prb=0.2\n\
\n\
";

#endif//__adddist_h__
