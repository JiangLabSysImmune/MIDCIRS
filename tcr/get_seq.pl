#for extracting reads seqs based on list
#
use strict;
use warnings;

my $orgfile = $ARGV[0];
my $lsfile = $ARGV[1];

my %or_info;
open (OR,"<$orgfile");
my $cur_name;
while(<OR>) {
	chomp;
	if ($_ =~ /\>(.+)/) {
		$cur_name = $1;
	}
	else{
		$or_info{$cur_name} = $_;
	}
}
close(OR);

open (LS,"<$lsfile");
while(<LS>) {
	chomp;
	my $name = $_;
	print ">$name\n$or_info{$name}\n";
}
close(LS);
