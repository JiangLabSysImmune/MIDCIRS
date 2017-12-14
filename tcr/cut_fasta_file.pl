#for cutting fasta file into certain files
#
use strict;
use warnings;

my $fafile = $ARGV[0];
my $cutnum = $ARGV[1];
my $outfile = $ARGV[2];


print "counting seq number\n";
my $seq_num = 0;
open (FA,"<$fafile");
while(<FA>) {
	chomp;
	if ($_ =~ /\>/) {
		$seq_num++;
	}
}
print "seq num:$seq_num\n";
close(FA);

my $sub_num = int($seq_num/$cutnum)+1;

open (FA,"<$fafile");
my $count = 0;
my $temp = "";
my $filecount = 0;
while(<FA>) {
	#chomp;
	my $line = $_;
	if ($line =~ /^\>/) {
		$count++;
	}
	
	if ($count > $sub_num) {
		$filecount++;
		my $outname = $outfile."_".$filecount.".fasta";
		open (OUT,">$outname");
		print OUT "$temp";
		print "$outname\n";
		close (OUT);
		$count = 0;
		$temp = "";
	}

	$temp = $temp.$line;
	if (eof(FA)) {
                $filecount++;
                my $outname = $outfile."_".$filecount.".fasta";
                open (OUT,">$outname");
                print OUT "$temp";
                print "$outname\n";
                close (OUT);
                $count = 0;
                $temp = "";
	}

}
close(FA);

	
	
