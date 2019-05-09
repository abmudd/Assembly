#!/usr/bin/env perl

# AUTHOR = Jessen Bredeson, Department of Molecular and Cell Biology, University of California, Berkeley.
$CONTACT = 'Jessen Bredeson <jessenbredeson@berkeley.edu>';

use strict;
use warnings;

use Getopt::Std;
use IO::Compress::Gzip;
use IO::Uncompress::Gunzip;

my %opts;
getopts('l:',\%opts) && @ARGV == 2 || die "Usage: single2pair.pl [-l INT] <in.fastq> <out.prefix>\n";

my $minlen = $opts{'l'} > 0 ? $opts{'l'} : 1;

my $GENERIC_HEADER = qr/^\@([^\s\/]+)/;
my $PAIRED_HEADER  = qr/^\@([^\s\/]+)#?[ATCG]*\/([12])/;
my $CASAVA_HEADER  = qr/^\@([^\s\/]+)\s([12]):[YN]:\d+:[ATCG]*/;
my $N  = 'N' x $minlen;
my $Q  = '#' x $minlen;

my $in = IO::Uncompress::Gunzip->new($ARGV[0]) || die "Cannot open $ARGV[0] for reading: $!\n";
my $o1 = IO::Compress::Gzip->new("$ARGV[1]_1.fastq.gz") || die "Cannot open $ARGV[1]_1.fastq.gz for writing: $!\n";
my $o2 = IO::Compress::Gzip->new("$ARGV[1]_2.fastq.gz") || die "Cannot open $ARGV[1]_2.fastq.gz for writing: $!\n";
while (my $h = <$in>) {
    my $s = <$in>;
    my $p = <$in>;
    my $q = <$in>;
    
    if ($h =~ $CASAVA_HEADER || $h =~ $PAIRED_HEADER || $h =~ $GENERIC_HEADER) {
	if ($2 == 1) { 
	    print($o1 "\@$1/1\n",$s,"+\n",$q); 
	    print($o2 "\@$1/2\n","$N\n","+\n","$Q\n");
	} else {
	    print($o1 "\@$1/1\n","$N\n","+\n","$Q\n");
	    print($o2 "\@$1/2\n",$s,"+\n",$q);
	}
    } else {
	die "Unsupported header format: $h\n";
    }
}

