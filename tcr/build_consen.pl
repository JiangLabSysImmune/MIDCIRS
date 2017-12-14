#for building the consensus based on sub-clustering info and Quality score info
#this is just for temporary use, when building the pipeline, will use JAVA 
#
use strict;
use warnings;

#input three files: 1). list file; 2).sequence file; 3).quality score file
my $lsfile = $ARGV[0];

my $sqfile = $ARGV[1];
my $scfile = $ARGV[2];

my %sq_info;
my $cur_name0;
open (SQ,"<$sqfile");
while (<SQ>) {
	chomp;
	my $line = $_;
	if ($line =~ /\>([AGCTN]+\_\d+)/) {
		$cur_name0 = $1;
	}
	else {
		$sq_info{$cur_name0} = $line;
	}
}
close(SQ);

my %sc_info;
my $cur_name;
open (SC,"<$scfile");
while (<SC>) {
	chomp;
	my $line = $_;
	if ($line =~ /^\|([AGCTN]+\_\d+)/) {
		$cur_name = $1;
	}
	else {
		$sc_info{$cur_name} = $line;
	}
}
close(SC);

open (LS,"<$lsfile");
while (<LS>) {
	chomp;
	my $line = $_;
	if ($line =~ /([AGCTN]+)\./) {
		my $bc = $1;

		my %gp_info;
		my $gpfile = "./levGroup/".$bc."_group.res";
		open (GP,"<$gpfile");
		while (<GP>) {
			chomp;
			my $line2 = $_;
			if ($line2) {
				my @cts = split(/\t/,$line2);
				$gp_info{$cts[0]} = $cts[1];
			}
		}
		close(GP);

		#find max sub-clust index
		my $max = 0;
		my $tmp;
		foreach my $name(sort keys %gp_info) {
			my $id = $gp_info{$name};
			$tmp = $name;
			if ($id > $max) {
				$max = $id;
			}
		}
		
		my $size = keys %gp_info;
		if ($size == 1) {
			my $cons_name = $bc."_"."0"."_"."1";
			my $cons_seq = $sq_info{$tmp};
			print ">$cons_name\n$cons_seq\n";
		}
		elsif ($size > 1) {
			for (my $i = 1;$i <= $max;$i++) {
				my %gp_seq;
				my %gp_sco;
				foreach my $name(sort keys %gp_info) {
					if ($gp_info{$name} == $i) {
						$gp_seq{$name} = $sq_info{$name};
						$gp_sco{$name} = $sc_info{$name};
					}
				}
				my $len = keys %gp_seq;
				my $cons_name = $bc."_".$i."_".$len;
				my $cons_seq = &consens(\%gp_seq,\%gp_sco);
				print ">$cons_name\n$cons_seq\n";
			}
		}
		#delete all small files within bc_groups,disRe,levRCM,levGroup
		my $seqfile = "./bc_groups/".$bc.".fasta";
		my $grpfile = "./levGroup/".$bc."_group.res";

		unlink $seqfile;
		unlink $grpfile;

		#deleting end
	}
}
close(LS);

#deleting directories
#system ("rm -rf ./bc_groups/");
#system ("rm -rf ./disRe/");
#system ("rm -rf ./levRCM/");
#system ("rm -rf ./levGroup/");
		






sub consens {
	my %gp_seq = %{$_[0]};
	my %gp_sco = %{$_[0]};
	
	my $sc_all = 100;

	my $cons = "";
	if (keys %gp_seq == 1) {
		my @key = keys %gp_seq;
		$cons = $gp_seq{$key[0]};
	}
	elsif (keys %gp_seq > 1) {
		my @key = keys %gp_seq;
		my $tmp = $gp_seq{$key[0]};
		for (my $i = 0;$i < length($tmp);$i++) {
			my %sub;
			$sub{"A"} = 0;
			$sub{"G"} = 0;
			$sub{"C"} = 0;
			$sub{"T"} = 0;
			$sub{"N"} = 0;

			foreach my $name(keys %gp_seq) {
				my $sub_seq = substr($gp_seq{$name},$i,1);
				my $sub_sco = substr($gp_sco{$name},$i,1);
				my $wgh = ord($sub_sco)/$sc_all;
				
				$sub{$sub_seq} = $sub{$sub_seq}+$wgh;
			}
			my $sub_cons = &find_max(\%sub);
			$cons = $cons.$sub_cons;
		}
	}
	return $cons;

}

sub find_max{
	#print "mmmmmm".$_[0];
	my %sub = %{$_[0]};

	#foreach my $kk(sort keys %sub) {
	#	print "$kk\t$sub{$kk}\n";
	#}


	my $max_char = "A";
	my $max = $sub{$max_char};

	foreach my $kk(keys %sub) {
		if ($sub{$kk} > $max) {
			$max = $sub{$kk};
			$max_char = $kk;
		}
	}
	return $max_char;
}


