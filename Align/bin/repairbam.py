#!/usr/bin/env python

# This script is derived from a script written by Jessen Bredeson <jessenbredeson@berkeley.edu>.

import argparse, collections, numpy, os, pysam, sys

def _by_reference(record):
    return (record.reference_id, record.reference_start)

def repair(queries, orientation, min_insert, max_insert):
    pidx = -1
    primary = [None, None]
    for query in queries:
        query.flag |= 0x1  # paired
        if query.is_secondary or query.is_supplementary:
            if pidx >= 0:
                query.flag |= 0x40 << pidx
        else:
            pidx += 1
            primary[pidx] = query
            query.flag |= 0x40 << pidx

    if primary[0] and primary[1]:
        primary = sorted(primary, key=_by_reference)
        if primary[0].is_unmapped or primary[1].is_unmapped:
            return
        if primary[0].reference_id != primary[1].reference_id:
            return
        len_insert = abs(primary[1].reference_end - primary[0].reference_start)
        if orientation >= 0 and \
           orientation != (int(primary[0].is_reverse) | int(primary[1].is_reverse) << 1):
            return
        if min_insert > len_insert:
            return
        if max_insert > min_insert and max_insert < len_insert:
            return

        for query in queries:
            query.flag |= 0x2  # proper

        return len_insert


parser = argparse.ArgumentParser(description='This script sets proper pairs based on provided insert size and oreintation.')
parser.add_argument("-p", "--pair", metavar='STR', help="Proper pair orientation: {FF,FR,RF,RR,all} [all]", type=str, default="ALL", required=True)
required = parser.add_argument_group('required arguments')
required.add_argument("-i", "--iave", metavar='INT', help="Average of insert size for proper pairing", type=int, required=True)
required.add_argument("-s", "--istd", metavar='INT', help="Std dev of insert size for proper pairing", type=int, required=True)
required.add_argument("in_bam", help="Input bam", type=str)
required.add_argument("in_config", help="Input config", type=str)
required.add_argument("out_prefix", help="Output prefix", type=str)
args = parser.parse_args()

proper_orientation = args.pair.upper()
max_isize = args.iave + (4*args.istd)
if (args.iave - (4*args.istd)) > 0: 
    min_isize = (args.iave - (4*args.istd))
else:
    min_isize = 0

if proper_orientation not in {'FF','FR','RF','RR','ALL'}:
    sys.stderr.write("Invalid enumerative value: %s\n" % proper_orientation)
    sys.exit(1)
    
proper_orientation = -1 if proper_orientation == 'ALL' \
    else (int(proper_orientation[0] == 'R') | int(proper_orientation[1] == 'R') << 1)

ibam = pysam.AlignmentFile(args.in_bam)
iconfig = open(args.in_config, 'r')
obam = pysam.AlignmentFile(args.out_prefix+'.bam', mode='wb', template=ibam)
oconfig = open(args.out_prefix+'.config', 'w')
oins = open(args.out_prefix+'.insert', 'w')

queries = collections.deque()
num_primary = 0
prev_query_name = None
all_len_ins = []
for curr in ibam:
    len_ins = False
    if curr.query_name != prev_query_name:
        if num_primary == 2:
            len_ins = repair(queries, orientation=proper_orientation, min_insert=min_isize, max_insert=max_isize)

        for query in queries:
            obam.write(query)

        if len_ins:
            all_len_ins.append(len_ins)
            oins.write(str(len_ins)+'\n')

        queries = collections.deque()
        num_primary = 0
        
    if not curr.is_secondary and not curr.is_supplementary:
        num_primary += 1

    prev_query_name = curr.query_name
    queries.append(curr)

if num_primary == 2:
    repair(queries, orientation=proper_orientation, min_insert=min_isize, max_insert=max_isize)
    
for query in queries:
    obam.write(query)

if len_ins:
    all_len_ins.append(len_ins)
    oins.write(str(len_ins)+'\n')

for item in iconfig:
    line = item.rstrip().split('\t')
    if line[0] in args.out_prefix:
        if len(all_len_ins) > 0:
            line[4] = str(int(numpy.average(all_len_ins)))
            line[5] = str(int(numpy.std(all_len_ins)))
            oconfig.write('\t'.join(line)+'\n')
            break
        else:
            oconfig.write(item)

ibam.close()
iconfig.close()
obam.close()
oconfig.close()
oins.close()
