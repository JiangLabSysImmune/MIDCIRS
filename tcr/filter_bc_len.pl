#based on the primer alignment results, keep only barcode of certain length
#
use strict;
use warnings;

my $algn_file = $ARGV[0];
my $rd_len = $ARGV[1];

open (IN,"<$algn_file");
my @lines = <IN>;
chomp(@lines);
close(IN);

foreach my $line(@lines) {
	if ($line =~ /^\>([AGCTN]+)\_(\d+)\_VS_(.+)\:\t(\d+)\t(\d+)\t(\d+)$/) {
		my $bc = $1;
		my $read_num = $2;
		my $st = $4;
		my $ed = $5;
		my $ms = $6;
		if ($rd_len == 150) {
			if ($ed == $rd_len-12) {
				print "$bc\_$read_num\n";
			}
		}
		else {
			print "$bc\_$read_num\n";
		}
	}
}
