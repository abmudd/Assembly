#!/usr/bin/env python

import argparse, datetime, subprocess, sys, os

# Function to check executables
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
        sys.stderr.write('ERROR: '+function+' not found in specified location or not executable\n')
        sys.exit(127)
    else:
        return which(function)

# Parse arguments
parser = argparse.ArgumentParser(description='Pipeline to run cactus with two species and filter the output. Cactus jobs are submitted to a SLURM cluster with sbatch using the defined memory, processor, and time flags, whereas filter jobs are run locally. All cactus variables must either be in the appropriate paths before running or have the -l flag set to load cactus as a module.')
parser.add_argument("-a", "--hal", metavar='STR', help="full path to hal directory [~/tools/hal]", type=str, default='~/tools/hal')
parser.add_argument("-c", "--chain", metavar='STR', help="chain file to liftover reference assembly [none]", type=str, default='stdin')
parser.add_argument("-f", "--filter", metavar='STR', help="netFilter synteny level: chimp for human/chimp, mouse for human/mouse, or none for no synteny netFilter [none]", type=str, default='none')
parser.add_argument("-g", "--gat", metavar='STR', help="full path to GenomeAlignmentTools directory [~/tools/GenomeAlignmentTools]", type=str, default='~/tools/GenomeAlignmentTools')
parser.add_argument("-k", "--kent", metavar='STR', help="full path to kent directory [~/tools/kent]", type=str, default='~/tools/kent')
parser.add_argument("-l", "--load", help="load cactus and samtools modules", action='store_true')
parser.add_argument("-m", "--memory", metavar='INT', help="RAM memory for job submission in GB [110]", type=int, default='110')
parser.add_argument("-n", "--nodedir", metavar='STR', help="SLURM node temp directory [$TMPDIR]", type=str, default='$TMPDIR')
parser.add_argument("-o", "--output", metavar='STR', help="output prefix [out]", type=str, default='out')
parser.add_argument("-p", "--processor", metavar='INT', help="processors for job submission [64]", type=int, default='64')
parser.add_argument("-r", "--colinrun", metavar='INT', help="cutoff size for runs of colinearity [1000]", type=int, default='1000')
parser.add_argument("-s", "--samtools", metavar='STR', help="path to samtools [samtools]", type=str, default='samtools')
parser.add_argument("-t", "--time", metavar='INT', help="hours for job submission [168]", type=int, default='168')
parser.add_argument("-v", "--version", action='version', version='%(prog)s 1.0')
required = parser.add_argument_group('required arguments')
required.add_argument("qry_name", help="name of query species in tree", type=str)
required.add_argument("qry_fasta", help="query fasta file", type=str)
required.add_argument("ref_name", help="name of reference species in tree", type=str)
required.add_argument("ref_fasta", help="reference fasta file", type=str)
required.add_argument("tree", help="phylogenetic tree [string or file]", type=str)
args = parser.parse_args()

# Output date, time, and command line
log = open(args.output+'.log','a')
log.write(str(datetime.datetime.now()).split('.')[0]+' - '+' '.join(sys.argv)+'\n')
log.close()

# Set output directory
output = args.output
if not os.path.exists(os.path.dirname(args.output)):
    output = os.getcwd()+'/'+args.output

# Cactus function
filter_flag = False
cactus_in = open(output+'.cactus.in','w')
if os.path.isfile(args.tree):
    for line in open(args.tree,'r'):
        cactus_in.write(line.rstrip()+'\n')
else:
    cactus_in.write(args.tree+'\n')
cactus_in.write(args.ref_name+' '+args.ref_fasta+'\n'+args.qry_name+' '+args.qry_fasta+'\n')
cactus_in.close()
if os.path.isfile(output+'.hal'): # Run filter function locally
    filter_flag = True
elif os.path.exists(os.path.dirname(output+'/jobStore')): # Restart cactus run
    cactus_sh = open(output+'.cactus.sh','w')
    cactus_sh.write('#!/bin/sh\n#SBATCH --cpus-per-task='+str(args.processor)+'\n#SBATCH --mem='
                    +str(args.memory)+'G\n#SBATCH --time='+str(args.time)+':0:0\ncd '+os.getcwd()+'\n')
    if args.load:
        cactus_sh.write('module purge\n')
        cactus_sh.write('unset PYTHONPATH\n')
        cactus_sh.write('module load cactus\n')
    cactus_sh.write('count=$RANDOM\n')
    cactus_sh.write('mkdir -p '+output+'/workdir\n')
    cactus_sh.write('cactus_args="--maxCores '+str(int(args.processor*0.9))+' --binariesMode local --workDir'
                    +' '+args.nodedir+'/workdir_$count '+args.nodedir+'/jobStore_$count '
                    +os.path.realpath(output+'.cactus.in')+' '+output+'.hal"\n')
    cactus_sh.write('mv '+output+'/workdir '+args.nodedir+'/workdir_$count &\n')
    cactus_sh.write('mv '+output+'/jobStore '+args.nodedir+'/jobStore_$count &\n')
    cactus_sh.write('wait\n')
    cactus_sh.write('timeout '+str(int(args.time-1))+'h cactus --restart $cactus_args &>>'+output+'.log\n')
    cactus_sh.write('mv '+args.nodedir+'/workdir_$count '+output+'/workdir &\n')
    cactus_sh.write('mv '+args.nodedir+'/jobStore_$count '+output+'/jobStore &\n')
    cactus_sh.write('wait\n')
    cactus_sh.close()
    subprocess.call('sbatch '+output+'.cactus.sh',shell=True) # Submit job to SLURM cluster
else: # Begin cactus run
    cactus_sh = open(output+'.cactus.sh','w')
    cactus_sh.write('#!/bin/sh\n#SBATCH --cpus-per-task='+str(args.processor)+'\n#SBATCH --mem='
                    +str(args.memory)+'G\n#SBATCH --time='+str(args.time)+':0:0\ncd '+os.getcwd()+'\n')
    if args.load:
        cactus_sh.write('module purge\n')
        cactus_sh.write('unset PYTHONPATH\n')
        cactus_sh.write('module load cactus\n')
    cactus_sh.write('count=$RANDOM\n')
    cactus_sh.write('mkdir -p '+output+'/workdir\n')
    cactus_sh.write('cactus_args="--maxCores '+str(int(args.processor*0.9))+' --binariesMode local --workDir'
                    +' '+args.nodedir+'/workdir_$count '+args.nodedir+'/jobStore_$count '
                    +os.path.realpath(output+'.cactus.in')+' '+output+'.hal"\n')
    cactus_sh.write('mv '+output+'/workdir '+args.nodedir+'/workdir_$count\n')
    cactus_sh.write('timeout '+str(int(args.time-1))+'h cactus $cactus_args &>>'+output+'.log\n')
    cactus_sh.write('mv '+args.nodedir+'/workdir_$count '+output+'/workdir &\n')
    cactus_sh.write('mv '+args.nodedir+'/jobStore_$count '+output+'/jobStore &\n')
    cactus_sh.write('wait\n')
    cactus_sh.close()
    subprocess.call('sbatch '+output+'.cactus.sh',shell=True) # Submit job to SLURM cluster

# Filter function
if not filter_flag:
    exit()
halLiftover = check_exe_return(os.path.expanduser(args.hal)+'/bin/halLiftover')
pslMap = check_exe_return(os.path.expanduser(args.kent)+'/src/pslMap')
axtChain = check_exe_return(os.path.expanduser(args.kent)+'/src/axtChain')
faToTwoBit = check_exe_return(os.path.expanduser(args.kent)+'/src/faToTwoBit')
chainPreNet = check_exe_return(os.path.expanduser(args.kent)+'/src/chainPreNet')
chainCleaner = check_exe_return(os.path.expanduser(args.gat)+'/bin/chainCleaner')
chainNet = check_exe_return(os.path.expanduser(args.kent)+'/src/chainNet')
netSyntenic = check_exe_return(os.path.expanduser(args.kent)+'/src/netSyntenic')
netFilter = check_exe_return(os.path.expanduser(args.kent)+'/src/netFilter')
netToAxt = check_exe_return(os.path.expanduser(args.kent)+'/src/netToAxt')
axtSort = check_exe_return(os.path.expanduser(args.kent)+'/src/axtSort')
axtToMaf = check_exe_return(os.path.expanduser(args.kent)+'/src/axtToMaf')
maf2stats = check_exe_return(os.path.dirname(os.path.expanduser(sys.argv[0]))+'/bin/maf2stats.py')
stats2colinearity = check_exe_return(os.path.dirname(os.path.expanduser(sys.argv[0]))
                                     +'/bin/stats2colinearity.py')
maf2indels = check_exe_return(os.path.dirname(os.path.expanduser(sys.argv[0]))+'/bin/maf2indels.py')
mafsInRegion = check_exe_return(os.path.expanduser(args.kent)+'/src/mafsInRegion')
if args.load:
    samtools = 'samtools'
else:
    samtools = check_exe_return(args.samtools)
filter_sh = open(output+'.filter.sh','w')
filter_sh.write('#!/bin/bash\nexport PATH='+os.path.expanduser(args.gat)+'/src/:'
                +os.path.expanduser(args.kent)+'/src:$PATH\n')
if args.load:
    filter_sh.write('module load samtools\n')
    filter_sh.write('if [[ ! -x "$(command -v samtools)" ]]; then\n')
    filter_sh.write('\techo "ERROR: samtools not found in specified location or not executable."\n')
    filter_sh.write('\texit 127\n')
    filter_sh.write('fi\n')
if not os.path.isfile(args.qry_fasta+'.fai'):
    filter_sh.write(samtools+' faidx '+args.qry_fasta+'\n')
if not os.path.isfile(args.ref_fasta+'.fai'):
    filter_sh.write(samtools+' faidx '+args.ref_fasta+'\n')
filter_sh.write('awk \'OFS="\\t" {print $1,"0",$2}\' '+args.qry_fasta+'.fai | '+halLiftover+' --outPSL '
                +output+'.hal '+args.qry_name+' stdin '+args.ref_name+' '+output+'.psl &\n')
filter_sh.write('sed "s/>/>'+args.qry_name+'./" '+args.qry_fasta+' >'+args.qry_name+'.sed.fasta &\n')
filter_sh.write('sed "s/>/>'+args.ref_name+'./" '+args.ref_fasta+' >'+args.ref_name+'.sed.fasta &\nwait\n')
filter_sh.write('cut -f1-2 '+args.ref_fasta+'.fai | awk \'{print "chain 0 "$1" "$2" + 0 "$2" "$1" "$2" + 0 '
                +'"$2" "NR"\\n"$2"\\n"}\' | '+pslMap+' -chainMapFile '+output+'.psl '+args.chain+' stdout | '
                +'awk \'OFS="\\t" {print $1,$2,$3,$4,$5,$6,$7,$8,$9,"'+args.qry_name+'."$10,$11,$12,$13,"'
                +args.ref_name+'."$14,$15,$16,$17,$18,$19,$20,$21}\' | sort -k10,10 -k12,12n -k13,13n >'
                +output+'.sign.psl\n')
filter_sh.write(axtChain+' -psl -faQ -faT -linearGap=loose '+output+'.sign.psl '+args.ref_name+'.sed.fasta '
                +args.qry_name+'.sed.fasta '+output+'.chain &>>'+output+'.log &\n')
filter_sh.write(faToTwoBit+' '+args.qry_name+'.sed.fasta '+args.qry_name+'.sed.2bit &\n')
filter_sh.write(faToTwoBit+' '+args.ref_name+'.sed.fasta '+args.ref_name+'.sed.2bit &\n')
filter_sh.write(samtools+' faidx '+args.qry_name+'.sed.fasta &\n')
filter_sh.write(samtools+' faidx '+args.ref_name+'.sed.fasta &\nwait\n')
filter_sh.write(chainPreNet+' '+output+'.chain '+args.ref_name+'.sed.fasta.fai '+args.qry_name
                +'.sed.fasta.fai '+output+'.prenet.chain &>>'+output+'.log\n')
filter_sh.write(chainCleaner+' '+output+'.prenet.chain '+args.ref_name+'.sed.2bit '+args.qry_name+'.sed.2bit'
                +' '+output+'.prenet.clean.chain '+output+'.prenet.clean.bed -tSizes='+args.ref_name+'.sed.'
                +'fasta.fai -qSizes='+args.qry_name+'.sed.fasta.fai -linearGap=loose &>>'+output+'.log\n')
filter_sh.write(chainNet+' -minSpace=1 -minFill=1 '+output+'.prenet.clean.chain '+args.ref_name+'.sed.fasta.'
                +'fai '+args.qry_name+'.sed.fasta.fai '+output+'.prenet.clean.'+args.ref_name+'.net '+output
                +'.prenet.clean.'+args.qry_name+'.net &>>'+output+'.log\n')
filter_sh.write(netSyntenic+' '+output+'.prenet.clean.'+args.ref_name+'.net '+output+'.prenet.clean.'
                +args.ref_name+'.syn.net &>>'+output+'.log\n')
if args.filter.lower() == 'chimp':
    filter_sh.write(netFilter+' -chimpSyn '+output+'.prenet.clean.'+args.ref_name+'.syn.net >'+output
                    +'.prenet.clean.'+args.ref_name+'.syn.filter.net\n')
elif args.filter.lower() == 'mouse':
    filter_sh.write(netFilter+' -syn '+output+'.prenet.clean.'+args.ref_name+'.syn.net >'+output
                    +'.prenet.clean.'+args.ref_name+'.syn.filter.net\n')
if args.filter.lower() == 'chimp' or args.filter.lower() == 'mouse':
    filter_sh.write(netToAxt+' '+output+'.prenet.clean.'+args.ref_name+'.syn.filter.net '+output
                    +'.prenet.clean.chain '+args.ref_name+'.sed.2bit '+args.qry_name+'.sed.2bit '+output
                    +'.filter.axt &>>'+output+'.log\n')
else:
    filter_sh.write(netToAxt+' '+output+'.prenet.clean.'+args.ref_name+'.syn.net '+output
                    +'.prenet.clean.chain '+args.ref_name+'.sed.2bit '+args.qry_name+'.sed.2bit '+output
                    +'.filter.axt &>>'+output+'.log\n')    
filter_sh.write(axtSort+' '+output+'.filter.axt '+output+'.filter.sort.axt &>>'+output+'.log\n')
filter_sh.write(axtToMaf+' '+output+'.filter.sort.axt '+args.ref_name+'.sed.fasta.fai '+args.qry_name
                +'.sed.fasta.fai '+output+'.filter.sort.maf &>>'+output+'.log\n')
filter_sh.write(maf2stats+' '+output+'.filter.sort.maf '+args.ref_name+' '+args.qry_name+' >'+output
                +'.filter.sort.stats\n')
filter_sh.write('tail -n +2 '+output+'.filter.sort.stats | sort -k5,5 -k7,7n | awk \'{print $0 "\\t" NR}\' '
                +'| sort -k1,1 -k3,3n | '+stats2colinearity+' - '+args.ref_name+' '+args.qry_name+' > '
                +output+'.filter.sort.colinearity\n')
filter_sh.write(maf2indels+' '+output+'.filter.sort.maf >'+output+'.filter.sort.indels\n')
filter_sh.write('awk \'{if ($9>='+str(args.colinrun)+') {print $1 "\\t" $3 "\\t" $4}}\' '+output
                +'.filter.sort.colinearity | cut -f2- -d\'.\' | '+mafsInRegion+' stdin '+output
                +'.filter.sort.colin'+str(args.colinrun)+'.maf '+output+'.filter.sort.maf\n')
filter_sh.close()
subprocess.call('sh '+output+'.filter.sh',shell=True)
exit()
