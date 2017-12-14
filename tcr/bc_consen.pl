#for buiding consensus directly based on barcode information instead of doing subclustering
#-singleton barcode: _0_1
#-barcode with read == 2: require two to be the same (_1_2), otherwise _0_2
#
#perl bc_cons.pl ./Jill52/all.list ./Jill52/Jill52.150.score ./Jill52/bc_groups
use strict;
use warnings;

#input list file, quality score file

my $lsfile = $ARGV[0];
my $scfile = $ARGV[1];
my $dir = $ARGV[2];

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
	my $line0 = $_;
	my @cts = split(/\t/,$line0);

	my $line = $cts[0];
	if ($line =~ /([AGCTN]+)\./) {
		my $bc = $1;
		
		my %gp_seq;
		my %gp_sco;

		my $gpfile = "$dir/$line"; #same as sequence file
		open (GP,"<$gpfile");
		my $name2;
		while (<GP>) {
                        chomp;
                        my $line2 = $_;
			if ($line2 =~ /^\>(.+)/) {
				$name2 = $1;
			}
			elsif ($line2 =~ /^[AGCTN]/) {
				$gp_seq{$name2} = $line2;
				$gp_sco{$name2} = $sc_info{$name2};
			}
                }
                close(GP);
		
		my $size = keys %gp_seq;
		if ($size == 1) {
			my $cons_name = $bc."_0_1";
			my @kk = keys %gp_seq;
			my $cons_seq = $gp_seq{$kk[0]};
			print ">$cons_name\n$cons_seq\n";
		}
		elsif ($size == 2) {
			my @kk = keys %gp_seq;
			my $nn1 = $kk[0];
			my $nn2 = $kk[1];
			if ($gp_seq{$nn1} ne $gp_seq{$nn2}) {
				my $cons_name1 = $bc."_0_1";
				my $cons_seq1 = $gp_seq{$nn1};
				my $cons_name2 = $bc."_0_2";
				my $cons_seq2 = $gp_seq{$nn2};
				print ">$cons_name1\n$cons_seq1\n>$cons_name2\n$cons_seq2\n";
			}
			else {
				my $cons_name = $bc."_1_2";
				my $cons_seq = $gp_seq{$nn1};
				print ">$cons_name\n$cons_seq\n";
			}
		}
		else {
			my $len = keys %gp_seq;
			my $cons_name = $bc."_1"."_".$len;
			my $cons_seq = &consens(\%gp_seq,\%gp_sco);
			print ">$cons_name\n$cons_seq\n";
		}

	}

}

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
        my %sub = %{$_[0]};

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





