#!/usr/bin/perl

#use strict;
#use warnings;

#This script obtains the relative position of SNPs within a gene. For example, a SNP is registered with its position within the genome
#nevertheless, for future uses, this script translates these positions to the relative position of the spliced transcript where it is located.
#In such case a SNP could have a 10,000 position in a chromosome but it could be located in the first position (0) of the gene or at a lower position relative to the gene start because of upstream introns.

#INPUT FILES
#List of genes with exons in bed12 format merged with the intersected SNPs (this is the output of overlapSelect from uscs genome browser):
#Chromosome	Starting_position	Ending_position	Name	Score	Strand	Thick_start	Thick_end	ItemRGB	Block_count	Block_sizes	Block_starts	SNP_Chromosome	SNP_Starting_position	SNP_Ending_position	SNP_RsID

#For more information visit: https://genome.ucsc.edu/FAQ/FAQformat.html#format1
#(Sorted_LncRNAs_ExonBlockCoords_Complete_Human_SecondDataset_StrandSpecific_clean_bed12_rsID.bed) #Check this name

#OUTPUT FILES
#Same information as the input: gene with exons in bed12 format, intersected SNP, plus the SNP's relative position to the spliced RNA.
#Gene_Chromosome	Gene_Starting_position	Gene_Ending_position	Gene_Name	Gene_Score	Gene_Strand	Gene_Thick_start	Gene_Thick_end	Gene_ItemRGB	Gene_Block_count	Gene_Block_sizes	Gene_Block_starts	SNP_Chromosome	SNP_Starting_position	SNP_Ending_position	SNP_RsID	Relative_position


#perl /LUSTRE/usuario/selene/JE_ROTACION/Data/Scripts/SNP_rel_position.pl -intersect /LUSTRE/usuario/selene/JE_ROTACION/Data/Sorted_LncRNAs_ExonBlockCoords_Complete_Human_SecondDataset_StrandSpecific_clean_bed12_rsID.bed
#output: /LUSTRE/usuario/selene/JE_ROTACION/Data/SNPs_relative_position/SNPs_rel_pos.txt

my ($i,$for_reverse,$position,$sum)=0;
my (@exon_start,@exon_length,@line)=();
my %opts=();

use Getopt::Long;
GetOptions(\%opts,'-intersect=s');

open(IN,"$opts{intersect}");
open(OUT,">/LUSTRE/usuario/selene/JE_ROTACION/Data/SNPs_relative_position/SNPs_rel_pos.txt");

while(<IN>){
	chomp($_);
	@line=split("\t",$_);
	if($line[5] eq "+"){	#For the forward strand genes
		
		(@exon_length,@exon_start)=();
	
		$line[10]=~s/,$//;
		$line[11]=~s/,$//;		
		@exon_length=split(",",$line[10]);
		@exon_start=split(",",$line[11]);
		
		$i=0;
		$sum=0;

		while((($line[13]-$line[1])>($exon_start[$i]+$exon_length[$i])) && $i<$#exon_start){	#This makes sure the SNP is within the considered exons. $i will increase until the SNP is equal or smaller than the current exon.	
		
			$sum+=$exon_start[$i+1]-($exon_start[$i]+$exon_length[$i]);	#Intergenic_size = next_exon_start - (prev_exon_start + prev_exon_length)
			$i++;	
			
#			if(($line[13]-$line[1])>($exon_start[$i]+$exon_length[$i])&&($line[13]-$line[1])<$exon_start[$i+1]){	#There is a case of an intronic SNP... this doesn't works: chr1	11444532	11445106	LocusForward_60499	0	+	11444532	11445106	0	2	147,78,	0,496,	chr1	11445097	11445098	rs537081094	... CHECK IT OUT
#				print OUT "PROBLEM:INTRONIC SNP ";
#			}
		}
			
		$position=$line[13]-$line[1]-$sum;	# $sum is intended to sum all the intergenic positions to substract it to the SNP position relative to the start of the gene (gene_start-SNP_start). It will consider all exons before the current one.

		#The relative positions are based on the starting positions of the SNP. A SNP located in the very beginning of the gene will have
		#a 0 position.

		print OUT "$_\t$position\n";
			
	}else{			#For the reverse strand genes. It is exactly the same as the previous code, but it also sums all exon lengths
				#to substract the current position. In that way it will be obtained a relative position starting at the ending position. reverse_relative_position = total_exon_lengths - forward_relative_position	
		(@exon_length,@exon_start)=();
	
		$line[10]=~s/,$//;
		$line[11]=~s/,$//;		
		@exon_length=split(",",$line[10]);
		@exon_start=split(",",$line[11]);
		
		$i=0;
		$sum=0;

		while((($line[13]-$line[1])>($exon_start[$i]+$exon_length[$i])) && $i<$#exon_start){		

			$sum+=$exon_start[$i+1]-($exon_start[$i]+$exon_length[$i]);	#Intergenic_size = next_exon_start - (prev_exon_start + prev_exon_length)
			print "($line[13]-$line[1]) > ($exon_start[$i]+$exon_length[$i])";			

			$i++;	
		}
		
		$for_reverse=0;
		foreach $i (@exon_length){
			$for_reverse+=$i;
		}
			
		$position=$for_reverse-($line[13]-$line[1]-$sum);	#The obtained position for the +strand will be substracted to the total spliced RNA ($for_reverse).

		print OUT "$_\t$position\n";
	
	}
	
	

}








