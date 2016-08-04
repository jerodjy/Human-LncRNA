#!/usr/bin/perl -X

if (!($ARGV[0])) {

        print "
exons2Bed.pl input.bed [filter]

Takes a 6 field bed with common identifiers and assembles them into a 12 field bed
file.
The optional filter field can be used to truncate accession IDs following a given
string.
Scores from common identifiers are summed and put into the score of the output file.

";
        exit;
}

if ($ARGV[1]) { $filter = $ARGV[1]; }

open IN,$ARGV[0] or die "Can't find $ARGV[0]\n";

while (<IN>) {
        $_ =~ s/\n|\r//g;
        my @fields = split("\t",$_);
        if (scalar @fields == 4) {
                ($chr,$start,$stop,$strnd) = @fields;
                ($id,$scre) = ();
        }
        else {
                ($chr,$start,$stop,$id,$scre,$strnd) = @fields;
        }
        if ($filter) { $id =~ s/$filter.*//g; }
        $lines->{$id}->{$start} = join "\t",($chr,$start,$stop,$id,$score,$strnd);
        $chrom{$id} = $chr;
        if ($strnd eq "") { $strand{$id} = "+"; }
        else { $strand{$id} = $strnd; }
        $score{$id} += $scre;
}
close IN;

my $reserved = "0";

foreach $id (sort keys %{$lines}) {
        my ($blockStarts,$blockSizes,$chrStart,$chrStop,$blockCount) = ();
        my @list = keys %{$lines->{$id}};
        foreach $start (sort { $a <=> $b } @list) {
                ($chr,$chrBlockStart,$chrBlockStop,$id,$scre,$strnd) =
split("\t",$lines->{$id}->{$start});
                if (!($chrStart)) { $chrStart = $chrBlockStart; }
                my $blockStart = $chrBlockStart - $chrStart;
                $blockStarts = $blockStarts."$blockStart,";
                my $blockSize = $chrBlockStop - $chrBlockStart;
                $blockSizes = $blockSizes."$blockSize,";
                $blockCount++;
        }
        $chrStop = $chrBlockStop;

        print join
"\t",($chrom{$id},$chrStart,$chrStop,$id,$score{$id},$strand{$id},$chrStart,$chrStop,$reserved,$blockCount,$blockSizes,$blockStarts);
        print "\n";
}
