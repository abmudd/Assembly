#!/usr/bin/env python

# Import modules
# ============================================================

import argparse, os, pysam, sys


# Parse arguments
# ============================================================

parser = argparse.ArgumentParser(description='This script splits up single- and multi-exon Trinity transcripts from the input bam. Single-exon transcripts should be later filtered with TransDecoder. Multi-exon transcripts are tossed if the first and/or last exon is less than 60 bp, if an intron is less than 60 bp or greater than 300,000 bp, or if the aligned transcript size is less than 250 bp.')
parser.add_argument('-o', '--output', metavar='STR', help='output file prefix [output]', type=str, default='output')
#parser.add_argument('-t', '--transdecoder', metavar='STR', help='path to TransDecoder.LongOrfs [TransDecoder.LongOrfs]', type=str, default='TransDecoder.LongOrfs')
parser.add_argument("-v", "--version", action='version', version='%(prog)s 1.0')
required = parser.add_argument_group('required arguments')
required.add_argument("bam", help="input bam of Trinity transcripts aligned to genome", type=str)
args = parser.parse_args()


# Function to check executables
# ============================================================

#def which(program):
#    def is_exe(fpath):
#        return os.path.exists(fpath) and os.access(fpath, os.X_OK)
#    fpath, fname = os.path.split(program)
#    if fpath:
#        if is_exe(program):
#            return os.path.abspath(program)
#    else:
#        for path in os.environ["PATH"].split(os.pathsep):
#            exe_file = os.path.join(path, program)
#            if is_exe(exe_file):
#                return exe_file
#    return None

#def check_exe_return(function):
#    if not which(function):
#        sys.stderr.write('ERROR: '+function+' not found in specified location or not executable\n')
#        sys.exit(127)
#    else:
#        return which(function)


# Check executables
# ============================================================

#TransDecoder_LongOrfs = check_exe_return(args.transdecoder)
#grep = check_exe_return('grep')
#cut = check_exe_return('cut')


# Set variables
# ============================================================

seen_single = [] # Check if single exon seq has already been written out
seen_multi = [] # Check if multi exon seq has already been written out
out_single = open(args.output+'.single_exon.fasta','w') # Open single exon output file
out_multi = open(args.output+'.multi_exon.fasta','w') # Open multi exon output file


###############################################################################
# Run
###############################################################################

for read in pysam.AlignmentFile(args.bam,'r'):

    if read.is_secondary or read.is_supplementary: # Ignore secondary and supplementary alignments
        continue

    else:
        # Extract intron start and end
        pos = read.pos
        cigar = read.cigar
        introns = []

        cigar_type = [i[0] for i in cigar]
        cigar_len = [i[1] for i in cigar]

        for i in [i for i, l in enumerate(cigar_type) if l == 3]:
            size = cigar_len[i]
            start = pos
            for j in range(len(cigar_type[:i])):
                if cigar_type[j] in [0, 2, 3]:
                    start += cigar_len[j]
            end = start + size - 1

            introns.append([start, end])

        # Identify start and end locations
        new_start = read.pos
        new_end = read.positions[-1]

        # Output single exon transcripts
        if not introns:
            if read.seq not in seen_single:
                out_single.write('>'+read.query_name+'\n'+read.seq+'\n')
                seen_single.append(read.seq)

        # Toss two exon transcripts where each exon is <60 bp
        elif len(introns) == 1 and introns[0][0] - (read.pos + 1) < 60 and new_end - introns[-1][1] < 60:
            continue

        # Filter remaining multi-exon transcripts
        else:
            keep = True
            new_pos = False

            if introns[0][0] - (read.pos + 1) < 60: # Set new start if first exon is <60 bp
                keep = False
#                new_start = introns[0][1] + 1
#                del introns[0]
#                new_pos = True

            if introns and new_end - introns[-1][1] < 60: # Set new end if last exon is <60 bp
                keep = False
#                new_end = introns[-1][0] - 1
#                del introns[-1]
#                new_pos = True

#            if new_pos:
#                temp_seq = ''
#                for i in range(0,len(read.seq)): # Trim read using new start / end locations of alignments
#                    if read.aligned_pairs[i][1] and read.aligned_pairs[i][1] >= new_start and read.aligned_pairs[i][1] <= new_end:
#                        temp_seq += read.seq[i]
#            else:
#                temp_seq = read.seq
            if read.alen < 250: # Toss transcript if aligned sequence is <250 bp
                keep = False

            for i in range(0,len(introns)): # Toss transcripts if intron(s) is <60 bp or >300000 bp
                if introns[i][1] - introns[i][0] < 60 or introns[i][1] - introns[i][0] > 300000:
                    keep = False

#            if len(temp_seq) < 250: # Toss transcript if aligned sequence is <250 bp
#                keep = False ##

            if keep and read.seq not in seen_multi: # Write out transcripts
                out_multi.write('>'+read.query_name+'\n'+read.seq+'\n')
                seen_multi.append(read.seq)


# Close output files
# ============================================================

out_single.close()
out_multi.close()

#subprocess.call(TransDecoder_LongOrfs+' -S -m 80 -t '+args.output+'.single_exon.fasta',shell=True)
#subprocess.call('paste - - <'+args.output+'.single_exon.fasta.transdecoder_dir/longest_orfs.cds | grep complete | sed \'s/\t/\n/g\' | cut -f3 -d\':\' | paste - - | awk \'{print ">" $0}\' | sed \'s/\t/\n/g\' >>'+args.output+'.final.fasta',shell=True)
