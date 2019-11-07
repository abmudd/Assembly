#!/usr/bin/env python

# Import modules
# ============================================================

import argparse, os, subprocess, sys


# Parse arguments
# ============================================================

parser = argparse.ArgumentParser(description='This script downloads organelle sequences from species close to the sequenced species using NCBI accessions, creates the respective NOVOPlasty config files, and runs NOVOPlasty.')
parser.add_argument("-e", "--efetch", metavar='STR', help="full path to efetch [efetch]", type=str, default='efetch')
parser.add_argument("-n", "--novoplasty", metavar='STR', help="full path to NOVOPlasty.pl [NOVOPlasty.pl]", type=str, default='seqkit')
parser.add_argument("-s", "--seqkit", metavar='STR', help="full path to seqkit [seqkit]", type=str, default='seqkit')
parser.add_argument("-v", "--version", action='version', version='%(prog)s 1.0')
parser.add_argument("-x", "--tenx", help="flag if input data is 10X [False]", action='store_true')
required = parser.add_argument_group('required arguments')
required.add_argument("-a", "--accessions", metavar='STR', help="list of NCBI accessions", type=str, required=True)
required.add_argument("-f", "--freads", metavar='STR', help="full path to forward reads", type=str, required=True)
required.add_argument("-i", "--isize", metavar='INT', help="insert size", type=str, required=True)
required.add_argument("-l", "--rlen", metavar='INT', help="read length", type=str, required=True)
required.add_argument("-r", "--rreads", metavar='STR', help="full path to reverse reads", type=str, required=True)
args = parser.parse_args()


# Function to check executables
# ============================================================

def which(program):
    def is_exe(fpath):
        return os.path.exists(fpath) and os.access(fpath, os.X_OK)
    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return os.path.abspath(program)
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file
    return None

def check_exe_return(function):
    if not which(function):
        sys.stderr.write('['+os.path.basename(sys.argv[0])+']: '+function+' not found in specified location or not executable\n')
        sys.exit(127)
    else:
        return which(function)


# Set variables
# ============================================================

efetch = check_exe_return(args.efetch)
novoplasty = check_exe_return(args.novoplasty)


# Convert 10X if applicable
# ============================================================

if args.tenx:
    seqkit = check_exe_return(args.seqkit)
    subprocess.check_call(seqkit+' subseq -j 8 -r 17:200 -o cut16.fastq.gz '+args.freads,shell=True)
    freads = os.getcwd()+'/cut16.fastq.gz'
else:
    freads = args.freads


# Write out config file
# ============================================================

batch = open('NOVOPlasty.batch', 'a')

for line in open(str(args.accessions), 'r'):
    item = line.rstrip()

    # Download accession from NCBI
    fasta = open(item+'.fa', 'w')
    subprocess.check_call(efetch+' -db nucleotide -id '+str(item)+' -format fasta',stdout=fasta,shell=True)
    fasta.flush()
    fasta.close()

    # Write out NOVOPlasty config file
    config = open(item+'.config', 'w')
    config.write('Project name         = '+item+'\n'
                 +'Insert size          = '+args.isize+'\n'
                 +'Insert size auto     = yes\n'
                 +'Read Length          = '+args.rlen+'\n'
                 +'Type                 = mito\n'
                 +'Genome Range         = 120000-200000\n'
                 +'K-mer                = 61\n'
                 +'Insert Range         = 1.6\n'
                 +'Insert Range strict  = 1.2\n'
                 +'Single/Paired        = PE\n'
                 +'Max memory           = 120\n'
                 +'Extended log         = 0\n'
                 +'Save assembled reads = no\n'
                 +'Combined reads       = \n'
                 +'Forward reads        = '+freads+'\n'
                 +'Reverse reads        = '+args.rreads+'\n'
                 +'Seed Input           = '+item+'.fa\n'
                 +'Reference            = '+item+'.fa\n')
    config.flush()
    config.close()

    # Write out batch file
    batch.write('perl '+novoplasty+' -c '+item+'.config\n')

batch.close()
subprocess.call('sh NOVOPlasty.batch', shell=True)
