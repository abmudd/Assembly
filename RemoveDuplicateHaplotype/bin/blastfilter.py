#!/usr/bin/env python

import os, pysam, sys

if len(sys.argv) != 3 or sys.argv[1] in ["-h","-help","--help"]:
    sys.stderr.write("Usage: %s <halfdepth.blast.out> <in.fa>\n" % os.path.basename(sys.argv[0]))
    sys.stderr.write("This script filters half-depth contigs and outputs the contigs to remove.\n")
    sys.exit(1)

fasta = pysam.FastaFile(sys.argv[2])
length = dict(zip(fasta.references, fasta.lengths))

blast_bit = {}
blast_alen = {}
blast_pid = {}
blast_convert = {}

def pop_all(value):
    blast_bit.pop(value)
    blast_alen.pop(value)
    blast_pid.pop(value)
    blast_convert.pop(value)

def set_all(value1, value2, bit, alen, pid):
    blast_bit[value1] = bit
    blast_bit[value2] = bit
    blast_alen[value1] = alen
    blast_alen[value2] = alen
    blast_pid[value1] = pid
    blast_pid[value2] = pid
    blast_convert[value1] = value2
    blast_convert[value2] = value1

def pop_and_set(value1, value2, value3, value4, bit, alen, pid):
    pop_all(value1)
    pop_all(value2)
    if value3 != value2 and value3 != value1:
        pop_all(value3)
    if value4 != value3 and value4 != value2 and value4 != value1:
        pop_all(value4)
    set_all(value1, value2, bit, alen, pid)

def partial_pop_and_set(value1, value2, value3, bit, alen, pid):
    pop_all(value1)
    pop_all(value3)
    set_all(value1, value2, bit, alen, pid)

file = sys.stdin
if sys.argv[1] != '-':
    file = open(sys.argv[1], 'r')

for blast in file:
    blast = blast.rstrip().split('\t')
    if blast[0] in blast_convert and blast[1] in blast_convert:
        if blast[11] > blast_bit[blast[0]] and blast[11] > blast_bit[blast[1]]:
            pop_and_set(blast[0], blast[1], blast_convert[blast[0]], blast_convert[blast[1]], blast[11], blast[3], blast[2])
        elif blast[11] == blast_bit[blast[0]] and blast[11] == blast_bit[blast[1]]:
            if blast[3] > blast_alen[blast[0]] and blast[3] > blast_alen[blast[1]]:
                pop_and_set(blast[0], blast[1], blast_convert[blast[0]], blast_convert[blast[1]], blast[11], blast[3], blast[2])
            elif blast[3] == blast_alen[blast[0]] and blast[3] == blast_alen[blast[1]]:
                if blast[2] > blast_pid[blast[0]] and blast[2] > blast_pid[blast[1]]:
                    pop_and_set(blast[0], blast[1], blast_convert[blast[0]], blast_convert[blast[1]], blast[11], blast[3], blast[2])
    elif blast[0] in blast_convert and blast[1] not in blast_convert:
        if blast[11] > blast_bit[blast[0]]:
            partial_pop_and_set(blast[0], blast[1], blast_convert[blast[0]], blast[11], blast[3], blast[2])
        elif blast[11] == blast_bit[blast[0]]:
            if blast[3] > blast_alen[blast[0]]:
                partial_pop_and_set(blast[0], blast[1], blast_convert[blast[0]], blast[11], blast[3], blast[2])
            elif blast[3] == blast_alen[blast[0]]:
                if blast[2] > blast_pid[blast[0]]:
                    partial_pop_and_set(blast[0], blast[1], blast_convert[blast[0]], blast[11], blast[3], blast[2])
    elif blast[0] not in blast_convert and blast[1] in blast_convert:
        if blast[11] > blast_bit[blast[1]]:
            partial_pop_and_set(blast[1], blast[0], blast_convert[blast[1]], blast[11], blast[3], blast[2])
        elif blast[11] == blast_bit[blast[1]]:
            if blast[3] > blast_alen[blast[1]]:
                partial_pop_and_set(blast[1], blast[0], blast_convert[blast[1]], blast[11], blast[3], blast[2])
            elif blast[3] == blast_alen[blast[1]]:
                if blast[2] > blast_pid[blast[1]]:
                    partial_pop_and_set(blast[1], blast[0], blast_convert[blast[1]], blast[11], blast[3], blast[2])
    else:
        set_all(blast[0], blast[1], blast[11], blast[3], blast[2])

seen = {}

for item in blast_convert:
    if item not in seen:
        length_i1 = length[item]
        length_i2 = length[blast_convert[item]]
        if length_i1 < length_i2:
            sys.stdout.write(item+'\n')
        else:
            sys.stdout.write(blast_convert[item]+'\n')
        seen[blast_convert[item]] = True
