#for extracting reads of FASTA file into barcode groups
#
use strict;
use warnings;

my $fs_file = $ARGV[0];

my %fs_info;
my $cur_name;
my $cur_bc;
open (FS,"<$fs_file");
while (<FS>) {
	chomp;
	my $line = $_;
	if ($line =~ /^\>([AGCTN]+)\_\d+/) {
		$cur_name = $line;
		#$cur_bc = $1;
	}
	elsif ($line =~ /^[AGCTN]/) {
		#my $line2 = $cur_name."\n".$line;
		#$fs_info{$line2} = $cur_bc;
		$fs_info{$cur_name} = $line;
	}
}
close(FS);

#sort into groups
my %rec;
foreach my $readname(sort keys %fs_info) {
	print "$readname\n";
	if ($readname =~ /^\>([AGCTN]+)\_\d+/) {
		my $bc = $1;
		open (OUT,">>./bc_groups/$bc.fasta");
		print OUT "$readname\n$fs_info{$readname}\n";
		close (OUT);
	}
}
