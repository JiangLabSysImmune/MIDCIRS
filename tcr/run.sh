javac *.java

ind=$1;

rd_len=250;   #250 for MiSeq
#FASTQ to FASTA
perl trans_raw.pl $ind

#group into barcodes
perl group_into_barcode.pl $ind

perl cut_seq_len.pl $ind.250.fasta 150 >$ind.150.fasta
perl cut_sc_len.pl $ind.250.score 150 >$ind.150.score

#align to the primer set
##name your constant primers into cons.primer.fasta
#java -XX:+UseSerialGC alignHam $ind.150.fasta cons.primer.fasta >all.algn.primer.res
perl filter_bc_len.pl all.algn.primer.res $rd_len >primer.valid.list
perl get_seq.pl $ind.150.fasta primer.valid.list >$ind.150.primer.fasta

#align to J and constant region
java -XX:+UseSerialGC alignHam $ind.150.primer.fasta JC.combine >algn.res
perl get_valid_list.pl algn.res >jc3.valid.list
perl get_seq.pl $ind.150.fasta jc3.valid.list >jc3.valid.fasta

rm -rf bc_groups
rm -rf levGroup

mkdir bc_groups
mkdir levGroup

#group into barcode groups
perl get_bcgroups.pl jc3.valid.fasta
rm -f all.list
cd bc_groups
ls >../all.list
cd ..

#do SUB-clustering
java allClust all.list;
#build consensus sequences
perl build_consen.pl all.list $ind.150.fasta $ind.150.score >consen.txt

#rm -rf ./bc_groups/
#rm -rf ./levGroup

#later stage-use IMGT/MIGEC to identif V-D-J etc...

