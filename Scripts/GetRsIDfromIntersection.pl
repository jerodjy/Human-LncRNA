use strict;

#This script uses a list of .bed positions to obtain their respective rsID from a file with the same first 3 columns and a 4th rsID column

#INPUT FILES
#Interest SNP list. -intersec format:
#Chromosome	Initial_position	Final_position
#RsID list. -snps format:
#Chromosome	Initial_position	Final_position	rsID

#OUTPUT FILES
#RS_intersection_list file:
#RsID
#RS_notfound:
#RsID

my $i=0;
my @line=();
my (%snp_list,%test_list)=();
my %opts =();

use Getopt::Long;
GetOptions(\%opts,'-intersec=s','-snps=s');

open(intrsc,"$opts{intersec}");
open(snps,"$opts{snps}");
open(OUTPUT,">/LUSTRE/usuario/selene/JE_ROTACION/Data/RS_intersection_list");
open(OUTPUT2,">/LUSTRE/usuario/selene/JE_ROTACION/Data/RS_notfound");

while(<intrsc>){
	chomp($_);
	@line=split("\t",$_);
	$snp_list{"$line[0]\t$line[1]\t$line[2]"}=1;
	$i++;
}

while(<snps>){
	chomp($_);
	@line=split("\t",$_);
	$test_list{"$line[0]\t$line[1]\t$line[2]"}=$line[3];
	if(exists($snp_list{"$line[0]\t$line[1]\t$line[2]"})){
		print OUTPUT "$line[3]\n";
	}else{
		print "RsID of 1000 genomes SNPs not present in the intersection: $_\n";
	}
}

open(intrsc,"$opts{intersec}");
while(<intrsc>){
	chomp($_);
        @line=split("\t",$_);
        if(!exists($test_list{"$line[0]\t$line[1]\t$line[2]"})){
                print OUTPUT2 "List of intersection rsID not found in 1000 genomes: $test_list{'$line[0]\t$line[1]\t$line[2]'}\n";
        }
}

print "Number of read SNPs from intersection file: $i";

close(OUTPUT);
close(snps);
close(intrsc);
close(OUTPUT2);
