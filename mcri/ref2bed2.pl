#!/usr/bin/perl -w

# estimate the length of the Introns/intron/transcript from GENCODE GTF file
# Contact: Shicheng Guo
# Version 1.3
# Update: 2017-04-03
# promoter: [Tss-2K,Tss+2K]
# Enhancer: [Tss+2k,Tss+5K]
# mm9: 797657 mm9.refGene.bed
# mm9.chrom.sizes, download from http://hgdownload.cse.ucsc.edu/goldenPath/mm9/bigZips/mm9.chrom.sizes

use strict;
use warnings;
die "Usage: perl $0 refGene.hg19.bed.txt hg19.chrom.sizes > refGene.bed\n" if scalar(@ARGV<2);

my $refGene=shift @ARGV;
my $ChromSize=shift @ARGV;

open F,$ChromSize;
my %chromesize;
while(<F>){
chomp;
my ($chr,$size)=split /\s+/;
$chromesize{$chr}=$size;
}

open IN, $refGene or die "Can't open $refGene\n";
my %refGene;
while(<IN>){
#    next if /"random"|"hap"|"chrUN"/;
    my ($bin,$NM,$chr,$strand,$txStart, $txEnd, $cdsStart, $cdsEnd, $exonCount, $exonStarts, $exonEnds,undef,$genesymbol,undef) = split /\t/;
    my @exonStarts = split(",", $exonStarts);
    my @exonEnds = split(",", $exonEnds);
    if($strand eq "+"){
    # gene in positive strand	
    my @enhancer=($txStart-4000,$txStart-2001);
    my @promter=($txStart-2000,$txStart+2000);
	my @UTR5=($txStart,$cdsStart);
	my @UTR3=($cdsEnd,$txEnd);
	next if $enhancer[0] <=0;
    next if $UTR3[1] >= $chromesize{$chr};
    print "$chr\t$enhancer[0]\t$enhancer[1]\t$strand\t$NM\t$genesymbol\tEnhancer\tEnhancer\n";
    print "$chr\t$promter[0]\t$promter[1]\t$strand\t$NM\t$genesymbol\tPromoter\tPromoter\n";
    print "$chr\t$UTR5[0]\t$UTR5[1]\t$strand\t$NM\t$genesymbol\tUTR5\tUTR5\n";
    my $exonrank;
    foreach my $i(0..($#exonStarts-1)){
    $exonrank++;
    print "$chr\t$exonStarts[$i]\t$exonEnds[$i]\t$strand\t$NM\t$genesymbol\tExon\tExon$exonrank\n";
    print "$chr\t$exonEnds[$i]\t$exonStarts[$i+1]\t$strand\t$NM\t$genesymbol\tIntron\tIntron$exonrank\n";
    }
    $exonrank++;
    print "$chr\t$exonStarts[$#exonStarts]\t$exonEnds[$#exonStarts]\t+\t$NM\t$genesymbol\tExon\tExon$exonrank\n";
    print "$chr\t$UTR3[0]\t$UTR3[1]\t$strand\t$NM\t$genesymbol\tUTR3\tUTR3\n";
    }elsif($strand eq '-'){
    # gene in negative strand
    my @enhancer=($txEnd+2001,$txEnd+4000);
    my @promter=($txEnd-2000,$txEnd+2000);
	my @UTR5=($cdsEnd,$txEnd);
	my @UTR3=($txStart,$cdsStart);
        next if $enhancer[1] >= $chromesize{$chr};
    next if $UTR3[0]<=0;
    print "$chr\t$enhancer[0]\t$enhancer[1]\t$strand\t$NM\t$genesymbol\tEnhancer\tEnhancer\n";
    print "$chr\t$promter[0]\t$promter[1]\t$strand\t$NM\t$genesymbol\tPromoter\tPromoter\n";
    print "$chr\t$UTR5[0]\t$UTR5[1]\t$strand\t$NM\t$genesymbol\tUTR5\tUTR5\n";
    my $exonrank;
    foreach my $i(reverse((1..$#exonStarts))){
    $exonrank++;
    print "$chr\t$exonStarts[$i]\t$exonEnds[$i]\t$strand\t$NM\t$genesymbol\tExon\tExon$exonrank\n";
    print "$chr\t$exonEnds[$i-1]\t$exonStarts[$i]\t$strand\t$NM\t$genesymbol\tIntron\tIntron$exonrank\n";
    }
    $exonrank++;
    print "$chr\t$exonStarts[0]\t$exonEnds[0]\t$strand\t$NM\t$genesymbol\tExon\tExon$exonrank\n";
    print "$chr\t$UTR3[0]\t$UTR3[1]\t$strand\t$NM\t$genesymbol\tUTR3\tUTR3\n";
    }else{
    die "Please check the input file reference.txt, the thrid column should be strand: + or -\n";
    }        
}

        