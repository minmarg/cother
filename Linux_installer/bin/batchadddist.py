#!/usr/bin/env python
##(C)2019-2021 Mindaugas Margelevicius, Institute of Biotechnology, Vilnius University
import sys, os, re
import fnmatch
import multiprocessing
from joblib import Parallel, delayed
from optparse import OptionParser

description = "Include pairwise distance information in each profile "+ \
    "present in the given directory."
proext = '.pro'
tproext = '.tpro'

bindir = os.path.dirname(os.path.abspath(__file__))
adddist = os.path.join(bindir, 'adddist')



def ParseArguments():
    """Parse command-line options.
    """
    parser = OptionParser(description=description)

    parser.add_option("--pro", dest="prodir",
                      help="Input directory of COMER profiles (with extension "+
                      "%s)."%(proext), 
                      metavar="DIR")
    parser.add_option("--ddd", dest="dstdir",
                      help="Input directory of distances files (see adddist). "+
                      "The prefix of the filenames should match the filenames of "+
                      "profiles in the directory specified by option --pro.", 
                      metavar="DIR")
    parser.add_option("--sfx", dest="suffix", type=str, default='_pred.prb',
                      help="Suffix for distances files (with extension). Distance "+
                      "files will be read from directory specified by option --dst "+
                      "with filenames constructed from the filenames of profiles "+
                      "appended by this suffix (default='%default').",
                      metavar="SUFFIX")
    parser.add_option("-o", "--outdir", dest="outdir",
                      help="Output directory of resulting cother profiles (with "+
                      "extension %s)."%(tproext),
                      metavar="DIR")
    parser.add_option("-p", "--optfile", dest="optfile",
                      help="Input file of options. By default, the options file "+
                      "in the installation directory is used. "+
                      "(Option passed to adddist.)",
                      metavar="FILE")
    parser.add_option("--dst", dest="dstlist", type=str, default='3',
                      help="Comma-separated list of columns in distances files "+
                      "representing distances to be averaged (default='%default'). "+
                      "(Option passed to adddist.)",
                      metavar="LIST")
    parser.add_option("--prb", dest="prob", type=float, default=0,
                      help="Probability lower-bound threshold at which to consider "+
                      "predicted distances. Probabilities are assumed to be given "+
                      "next to distances (see adddist) (default='%default'). "+
                      "(Option passed to adddist.)",
                      metavar="PROB")
    parser.add_option("--cpus", dest="cpus", type=int, default=multiprocessing.cpu_count()//2,
                      help="Number of CPUs to use; default=%default",
                      metavar="CPUs")

    (options, args) = parser.parse_args()

    if not options.prodir:
        sys.exit("ERROR: Input directory of profiles is not specified.\n")
    if not os.path.isdir(options.prodir):
        sys.exit("ERROR: Invalid input directory of profiles: "+str(options.prodir)+"\n")

    if not options.dstdir:
        sys.exit("ERROR: Input directory of distances files is not specified.\n")
    if not os.path.isdir(options.dstdir):
        sys.exit("ERROR: Invalid input directory of distances files: "+str(options.dstdir)+"\n")

    if options.optfile and not os.path.isfile(options.optfile):
        sys.exit("ERROR: Invalid options file specified: "+str(options.optfile)+"\n")

    if options.prob < 0 or 1 < options.prob:
        sys.exit("ERROR: Invalid probability specified: "+str(options.prob)+"\n")

    if options.cpus < 1:
        sys.exit("ERROR: Invalid number of CPUs specified: "+str(options.cpus)+"\n")

    if not options.outdir:
        sys.exit("ERROR: Output directory of cother profiles is not specified.\n")
    if not os.path.exists(options.outdir):
        os.mkdir(options.outdir)
    elif os.path.isfile(options.outdir):
        sys.exit("ERROR: Output directory is a file: %s"%(options.outdir))

    return options



def getFilelist(directory, ext):
    """Read `directory' for files with extension `ext', sort them by 
    size and return the resulting list.
    """
    files = []

    with os.scandir(directory) as entries:
        for entry in entries:
            if entry.is_file() and fnmatch.fnmatch(entry, '*' + ext):
                files.append(entry.name.rsplit(ext)[0]) ##add filename without extension

    for i in range(len(files)):
        files[i] = (files[i], os.path.getsize(os.path.join(directory, files[i] + ext)))

    #sort files by size
    files.sort(key=lambda name: name[1], reverse=True)

    for i in range(len(files)):
        files[i] = files[i][0]

    return files



def processProfile(profilename, options):
    """Process one profile by adding distance information to it and 
    writing the resulting cother profile to file. 
    """
    dstfile = os.path.join(options.dstdir, profilename + options.suffix)

    if not os.path.isfile(dstfile):
        sys.stderr.write('WARNING: Distances file not found: %s\n'%(dstfile))
        return

    profile = os.path.join(options.prodir, profilename + proext)
    outfile = os.path.join(options.outdir, profilename + tproext)

    if os.path.isfile(outfile):
        return

    optfileopt = ''
    if options.optfile: optfileopt = '-p ' + options.optfile

    cmdline = '%s -v -i %s -j %s -o %s %s --dst=%s --prb=%s'%(
        adddist, dstfile, profile, outfile, optfileopt,
        options.dstlist, options.prob)

    #sys.stderr.write('%s\n'%(cmdline))
    os.system(cmdline)



if __name__ == "__main__":

    options = ParseArguments()

    if not os.path.isfile(adddist):
        sys.exit("ERROR: Program adddist not found: "+adddist+"\n")

    profiles = getFilelist(options.prodir, proext)

    # parallel execution
    n_cores = options.cpus
    if n_cores < 1: n_cores = 1

    Parallel(n_jobs=n_cores, prefer="threads")(delayed(processProfile)(
        profiles[i], options) for i in range(len(profiles)))

