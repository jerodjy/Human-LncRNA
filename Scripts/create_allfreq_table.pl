use strict;
#use Data:Dumper;

#This script takes a concatenated table (by column) of --freq results from plink. It uses a threshold to simplify the data, assuming fixation when crossing that threshold of an allele
#or marking an NA otherwise.

#INPUT FILES
#Concatenated table of --freq results from plink:
#Pop1_Chromosome	Pop1_RsID	Pop1_A1	Pop1_A2	Pop1_MAF	Pop1_NCHROBS Pop2_Chromosome	Pop2_RsID	Pop2_A1	Pop2_A2	...
#Populations list:
#Population

#OUTPUT FILES
#File with all alleles from all populations assuming fixation with a given threshold
#RsID	Pop1_allele	Pop2_allele	Pop3_allele	...
#File with only the RsIDs with more than 1 fixated allele and with variation among these alleles of a certain RsID.
#RsID   Pop1_allele     Pop2_allele     Pop3_allele     ...

my ($i,$cut,$temp,$maf,$temp2,$useful)=0;
my (@pop_array,@line,@NA)=();
my %frq_hash=();
my %opts=();

use Getopt::Long;
GetOptions(\%opts,'-table=s','-threshold=f','-poplist=s');

open(tab,"$opts{table}");
open(pops,"$opts{poplist}");
open(OUT,">/LUSTRE/usuario/selene/JE_ROTACION/Data/Frequencies/Filtered_table/All_pops_lncRNA_frq_filtered.txt");
open(FILTER,">/LUSTRE/usuario/selene/JE_ROTACION/Data/Frequencies/Filtered_table/All_pops_lncRNA_frq_filtered_informative90.txt");
$cut=$opts{threshold};
print $cut;

while(<pops>){
	chomp($_);
	push(@pop_array,$_);
}

$temp=join("\t",@pop_array);
print OUT "RsID\t$temp\n";
print FILTER "RsID\t$temp\n";

while(<tab>){
	chomp($_);
	$_=~s/\s+/\t/g;
	@line=split("\t",$_);
	$frq_hash{$line[1]}[0]="NA";
	for($i=0;$i<$#line;$i=$i+6){ #Check this numbers
		if($i<6){
			$maf=$line[2];
			if($line[4]<$cut){
				$frq_hash{$line[1]}[0]=$line[2];
			}else{
                               	$frq_hash{$line[1]}[0]="NA";	#Change for NA	
			}
		}elsif($line[$i+4]<$cut){
                        push(@{$frq_hash{$line[1]}},$line[$i+2]);
                }else{
                        push(@{$frq_hash{$line[1]}},"NA"); #Change for NA
		}
	}

	$temp=join("\t",@{$frq_hash{$line[1]}});
	print OUT "$line[1]\t$temp\n";


	@NA=$temp=~/NA/g;	
	print (scalar @NA);
	print " $#{$frq_hash{$line[1]}} ";
	if((scalar @NA)<$#{$frq_hash{$line[1]}}){	#First I make sure there are more than 2 fixed alleles to evaluate afterwards if all fixed alleles are the same (which are not of interest to us)
		for($i=0;$i<=$#{$frq_hash{$line[1]}};$i++){	#Then all fixed alleles are compared to discard all RsIDs without variation
			if($frq_hash{$line[1]}[$i] ne "NA"){	#This works
				if($temp2 ne "A" && $temp2 ne "T" && $temp2 ne "G" && $temp2 ne "C"){			#Here used to be the problem
					$temp2=$frq_hash{$line[1]}[$i];
 	                             	print "$frq_hash{$line[1]}[$i]p-$temp2";
				}elsif($temp2 ne $frq_hash{$line[1]}[$i]){
  	                            	print "$frq_hash{$line[1]}[$i]u-$temp2";	
					$useful=1;
				}else{
					print "$frq_hash{$line[1]}[$i]-$temp2 ";
				}
			}
		}	
                        print "\n";	
	}

	if($useful==1){        
		print FILTER "$line[1]\t$temp\n";
	}
	$useful=0;
	$temp2=0;

}

close(tab);
close(OUT);
close(pops);
close(FILTER);
