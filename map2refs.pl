use warnings;
use strict;

use Getopt::Std;

my %opts;
getopts('r:q:m:l:', \%opts);

my $lib      = $opts{l};
my $readdir  = $opts{q};
my $mapdir   = $opts{m};
my $refs     = $opts{r};

my $bam     = mapFiles($readdir,$lib,$refs,$mapdir);

sub mapFiles {
    my ($readdir,$lib,$refs,$mapdir) = @_;
    my $file1 = $readdir . $lib . '_1_final.fastq.gz';
    my $file2 = $file1; $file2 =~ s/_1_/_2_/;
    my $fileu = $file1; $fileu =~ s/_1_/_u_/;

    my $sam      = $mapdir . $lib . ".sam";
    my $bam      = $mapdir . $lib . ".bam";
    my $finalout = $mapdir . $lib . ".sorted";

    my $call1 = system("bowtie2 -5 5 -3 5 -p 1 -x $refs -U $file1,$file2,$fileu -S $sam");
    my $call2 = system("samtools view -bS $sam > $bam");
    my $call3 = system("samtools sort $bam $finalout");
    $finalout .= ".bam";
    return($finalout);
}


