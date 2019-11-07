#!/usr/bin/env python

import itertools, os, pysam, sys

if len(sys.argv) != 3 or sys.argv[1] in ["-h","-help","--help"]:
    sys.stderr.write("Usage: "+os.path.basename(sys.argv[0])+" <read1.bam> <read2.bam>\n")
    sys.stderr.write("This script parses bam files and output the SSPACE tab file.\n")
    sys.stderr.write("Note: Bam files should be sorted by name.\n")
    sys.exit(1)

seen = {}
with pysam.AlignmentFile(sys.argv[1], 'rb') as bam1, pysam.AlignmentFile(sys.argv[2], 'rb') as bam2:
    for read1, read2 in itertools.izip(bam1,bam2):
        if read1.is_unmapped or read2.is_unmapped: 
            continue
        elif read1.reference_name == read2.reference_name: 
            continue
        elif read1.mapping_quality <= 30 or read2.mapping_quality <= 30:
            continue
        elif read1.is_duplicate or read2.is_duplicate: 
            continue
        elif read1.is_qcfail or read2.is_qcfail: 
            continue
        elif read1.is_secondary or read2.is_secondary: 
            continue
        elif read1.is_supplementary or read2.is_supplementary: 
            continue
        elif not read1.has_tag("NM") or not read2.has_tag("NM"):
            continue
        elif read1.get_tag("NM") >= 5 or read2.get_tag("NM") >= 5:
            continue
        elif float(read1.query_alignment_length)/read1.query_length <= 0.9 or float(read2.query_alignment_length)/read2.query_length <= 0.9:
            continue
        else:
            query = read1.reference_name+':'+str(read1.reference_start)+':'+str(read1.reference_end)+':'+read2.reference_name+':'+str(read2.reference_start)+':'+str(read2.reference_end)
            query_rev = read2.reference_name+':'+str(read2.reference_start)+':'+str(read2.reference_end)+':'+read1.reference_name+':'+str(read1.reference_start)+':'+str(read1.reference_end)
            if query not in seen and query_rev not in seen:
                sys.stdout.write(read1.reference_name+'\t'+str(read1.reference_start)+'\t'+str(read1.reference_end)+'\t'+read2.reference_name+'\t'+str(read2.reference_start)+'\t'+str(read2.reference_end)+'\n')
                seen[query] = True
