my %db;
my $rs=shift @ARGV;
open DB,"ROI.refGene.vcf.txt";
my %rsimpute;
while(<DB>){
chomp;
my @line=split/\s+/;
if(/IMPUTED/){
$rsimpute{$line[2]}="imputed";
}else{
$rsimpute{$line[2]}="typed";
}
}
open DB1,"ROI.assoc.logistic";
while(<DB1>){
chomp;
next if $_ !~ /ADD/i;
my @line=split/\s+/;
$db{$line[2]}="$line[2]\t$line[3]\t$line[12]\t$rsimpute{$line[2]}";
}
close DB1;

my %snp;
open DB2, "$rs.ld";
print "SNP\tPOS\tPVAL\tTYPE\tRSQR\n";
while(<DB2>){
chomp;
next if /CHR_A/;
my @line=split /\s+/,$_;
next if ! defined $db{$line[6]};
next if defined $snp{$line[6]};
print "$db{$line[6]}\t$line[7]\n";
$snp{$line[6]}=$line[6];
}
