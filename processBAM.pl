use warnings;
use strict;
use Getopt::Std;

my %opts;
getopts('l:r:m:p:g:a:', \%opts);

my $lib      = $opts{l};
my $mapdir   = $opts{m};
my $refs     = $opts{r};
my $seqdict  = $opts{p};
my $gatk     = $opts{g};
my $addrg    = $opts{a};



my $bam     = $mapdir . $lib . ".sorted.bam";
my $bamrg   = rg($lib, $bam, $addrg);
my $ibamrg  = indbam($bamrg);
my $bamproc = processBam($ibamrg, $refs, $seqdict, $gatk);


### sub routines ###


sub rg {
    my ($lib, $bam, $addrg)   = @_;
    my $bamrg         = $bam;
       $bamrg         =~ s/.bam/.rg.bam/;
    my ($lane, $index) = split(/_/,$lib);

    my $AddOrReplCall = "java -Xmx8g -jar $addrg";
    print "\nread group file name: $bamrg\n";
    system("$AddOrReplCall INPUT=$bam OUTPUT=$bamrg RGID=$lib RGLB=$lib RGPU=$lane RGPL=illumina RGSM=$lib");
    return($bamrg);
}

sub indbam {
    my $bamrg   = $_[0];
    my $ibamrg = $bamrg;
    $ibamrg =~ s/\.bam/\.ind\.bam/;
    system("mv $bamrg $ibamrg");
    system("samtools index $ibamrg");
    return($ibamrg);
}



sub processBam {
    my ($bam, $ref, $seqdict, $gatk) = @_;

    my $bamint = $bam;
       $bamint =~ s/.bam/.intervals/g;

    my $bamproc = $bam;
       $bamproc =~ s/.sorted.rg.ind.bam/.proc.bam/g;

    my $seqdictCall = "java -Xmx8g -jar $seqdict";
    my $gatkCall    = "java -Xmx8g -jar $gatk";
    my $dict  = $ref;
       $dict  =~ s/\.fasta/\.fasta.fai/;

    my $rtc   = "RealignerTargetCreator";
    my $ir    = "IndelRealigner";

   # unless (-e $dict) { system("$seqdictCall R=$ref O=$dict"); }
    system("$gatkCall -R $ref -T $rtc  -I $bam -o $bamint");
    system("$gatkCall -R $ref -T $ir -I $bam --targetIntervals $bamint -o $bamproc");

    return($bamproc);
}



__END__;

