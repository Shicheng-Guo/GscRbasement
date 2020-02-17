open DB,"/gpfs/home/guosa/hpc/rheumatology/RA/RA500/RA2020-B1.fam";
my %db;
while(<DB>){
my @line=split/\s+/;
$db{$line[1]}=$_;
}
close DB;
open F, shift @ARGV;
while(<F>){
my @line=split /\s+/,$_;
if(defined $db{$line[1]}){
print $db{$line[1]};
}else{
print "$_";
}
}