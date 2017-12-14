#this is for extracting valid reads (e.g. lev distance mismatch <= 3) out from the alignment results;
#
use strict;
use warnings;

my $algn_file = $ARGV[0];


open (IN,"<$algn_file");
my @lines = <IN>;
chomp(@lines);
close(IN);

my $thr = 3;

foreach my $line(@lines) {
	if ($line =~ /^\>([AGCTN]+)\_(\d+)\_VS_(.+)\:\t(\d+)\t(\d+)\t(\d+)$/) {
		my $bc = $1;
		my $read_num = $2;
		my $st = $4;
		my $ed = $5;
		my $ms = $6;
                if ($ms <= $thr && $ed > $st) {
			print "$bc\_$read_num\n";
		#	print "$line\n";
		}
		#else {
		#	print "$bc\_$read_num\n";
		#}
	}
}

