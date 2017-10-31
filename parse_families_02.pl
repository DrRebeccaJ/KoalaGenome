#!/net/isi-software/server/bin/perl
use strict;
use warnings;
use lib '/tgac/software/testing/perl/5.22.1/x86_64/lib/site_perl/5.22.0/';
use Statistics::Descriptive;
use POSIX;
use lib '/tgac/software/testing/perl/5.22.1/x86_64/lib/site_perl/5.22.0/Scalar-List-Utils-1.45/lib';
use List::Util 'shuffle';
use List::Util 'min';
use List::Util 'max';

my $project=$ARGV[0];
#CREATE_LOG(\@ARGV);

my ($input, $outpath)=@ARGV[1 .. $#ARGV];

#exit;

my $count_1_up=0;
my $count_1_down=0;
my $count_2_up=0;
my $count_2_down=0;
my $count_3_up=0;
my $count_3_down=0;
my $count_4_up=0;
my $count_4_down=0;

my @groups=("KoPC", "Shar", "Meug", "Hsap", "Mmus", "Cfam", "Mdom", "Oana", "Ggal");
my %ortho_groups;
open(OUT, ">$outpath/families_parsed.out");
open(OUT1, ">$outpath/koala_expansions_contractions.out");
open(OUT2, ">$outpath/marsupials_expansions_contractions.out");
open(OUT3, ">$outpath/diprotodontia_expansions_contractions.out");
open(OUT4, ">$outpath/autralidelphia_expansions_contractions.out");
print OUT join ("\t", "GROUP", "PATTERN", "nb_KoPC", "nb_Shar", "nb_Meug", "nb_Hsap", "nb_Mmus", "nb_Cfam", "nb_Mdom", "nb_Oana", "nb_Ggal", "KoPC", "Shar", "Meug", "Hsap", "Mmus", "Cfam", "Mdom", "Oana", "Ggal"),"\n";
my %done;
open(IN, "$input")||die"IN $input\n";
while(<IN>){
	chomp;
	my @split=split /\t/, $_;
	for my $n (1 .. $#split){
		my @var=split /\|/, $split[$n];
		push @{$ortho_groups{$var[0]}}, $var[1];
	}
	my @results;
	my @counts;
	my @genes;
	for my $i (0 .. $#groups){
		if($ortho_groups{$groups[$i]}){
			if(scalar @{$ortho_groups{$groups[$i]}} == 1){
				push @results, 1;
				push @genes, @{$ortho_groups{$groups[$i]}};
			}
			else{
				push @results, "m";
				push @genes, join (";", @{$ortho_groups{$groups[$i]}});
			}
			push @counts, scalar @{$ortho_groups{$groups[$i]}};
		}
		else{
			push @results, 0;
			push @counts, 0;
			push @genes, "NA";
		}
	}
	my $min_1=0;
	my $min_2=0;
	my $max_1=0;
	my $max_2=0;
	my @temp1=("KoPC", "Shar", "Meug", "Mdom");
	my @temp2=("Hsap", "Cfam", "Mmus");
	for my $j (0 ..  $#temp1){
		my $k=0;
		if($ortho_groups{$temp1[$j]}){
			if($j==0){
				$min_1 = scalar @{$ortho_groups{$temp1[$j]}};
				$max_1 = scalar @{$ortho_groups{$temp1[$j]}};
			}
			else{
				if(scalar @{$ortho_groups{$temp1[$j]}} > $max_1){
					$max_1 = scalar @{$ortho_groups{$temp1[$j]}};
				}
				if(scalar @{$ortho_groups{$temp1[$j]}} < $min_1){
					$min_1 = scalar @{$ortho_groups{$temp1[$j]}};
				}
			}
		}
		else{
			$min_1=0;
		}
	}
	for my $j (0 ..  $#temp2){
		my $k=0;
		if($ortho_groups{$temp2[$j]}){
			if($j==0){
				$min_2 = scalar @{$ortho_groups{$temp2[$j]}};
				$max_2 = scalar @{$ortho_groups{$temp2[$j]}};
			}
			else{
				if(scalar @{$ortho_groups{$temp2[$j]}} > $max_2){
					$max_2 = scalar @{$ortho_groups{$temp2[$j]}};
				}
				if($min_2 > scalar @{$ortho_groups{$temp2[$j]}}){
					$min_2 = scalar @{$ortho_groups{$temp2[$j]}};
				}
			}
		}
		else{
			$min_2=0;
		}
	}
	if($min_1 > $max_2 || $max_1 < $min_2){
		$done{$split[0]}=1;
		if($min_1 > $max_2){
			$count_1_up++;
			print OUT2 join ("\t", "up", $split[0], join ("", @results), @counts, @genes), "\n";
		}
		elsif($max_1 < $min_2){
			$count_1_down++;
			print OUT2 join ("\t", "down", $split[0], join ("", @results), @counts, @genes), "\n";
		}
	}
	if(!$done{$split[0]}){
		$min_1=0;
		$min_2=0;
		$max_1=0;
		$max_2=0;
		@temp1=();
		@temp2=();
		@temp1=("KoPC", "Shar", "Meug");
		@temp2=("Hsap", "Cfam", "Mmus", "Mdom");
		for my $j (0 ..  $#temp1){
			my $k=0;
			if($ortho_groups{$temp1[$j]}){
				if($j==0){
					$min_1 = scalar @{$ortho_groups{$temp1[$j]}};
					$max_1 = scalar @{$ortho_groups{$temp1[$j]}};
				}
				else{
					if($min_1 > scalar @{$ortho_groups{$temp1[$j]}}){
						$min_1 = scalar @{$ortho_groups{$temp1[$j]}};
					}
					if($max_1 < scalar @{$ortho_groups{$temp1[$j]}}){
						$max_1 = scalar @{$ortho_groups{$temp1[$j]}};
					}
				}
			}
			else{
				$min_1=0;
			}
		}
		for my $j (0 ..  $#temp2){
			my $k=0;
			if($ortho_groups{$temp2[$j]}){
				if($j==0){
					$min_2 = scalar @{$ortho_groups{$temp2[$j]}};
					$max_2 = scalar @{$ortho_groups{$temp2[$j]}};
				}
				else{
					if($min_2 > scalar @{$ortho_groups{$temp2[$j]}}){
						$min_2 = scalar @{$ortho_groups{$temp2[$j]}};
					}
					if($max_2 < scalar @{$ortho_groups{$temp2[$j]}}){
						$max_2 = scalar @{$ortho_groups{$temp2[$j]}};
					}
				}
			}
			else{
				$min_2=0;
			}
		}
		if($min_1 > $max_2 || $max_1 < $min_2){
			$done{$split[0]}=1;
			if($min_1 > $max_2){
				$count_2_up++;
				print OUT4 join ("\t", "up", $split[0], join ("", @results), @counts, @genes), "\n";
			}
			elsif($max_1 < $min_2){
				$count_2_down++;
				print OUT4 join ("\t", "down", $split[0], join ("", @results), @counts, @genes), "\n";
			}
		}
	}
	if(!$done{$split[0]}){
		$min_1=0;
		$min_2=0;
		$max_1=0;
		$max_2=0;
		@temp1=();
		@temp2=();
		@temp1=("KoPC", "Meug");
		@temp2=("Hsap", "Cfam", "Mmus", "Shar", "Mdom");
		for my $j (0 ..  $#temp1){
			my $k=0;
			if($ortho_groups{$temp1[$j]}){
				if($j==0){
					$min_1 = scalar @{$ortho_groups{$temp1[$j]}};
					$max_1 = scalar @{$ortho_groups{$temp1[$j]}};
				}
				else{
					if($min_1 > scalar @{$ortho_groups{$temp1[$j]}}){
						$min_1 = scalar @{$ortho_groups{$temp1[$j]}};
					}
					if($max_1 < scalar @{$ortho_groups{$temp1[$j]}}){
						$max_1 = scalar @{$ortho_groups{$temp1[$j]}};
					}
				}
			}
			else{
				$min_1=0;
			}
		}
		for my $j (0 ..  $#temp2){
			my $k=0;
			if($ortho_groups{$temp2[$j]}){
				if($j==0){
					$min_2 = scalar @{$ortho_groups{$temp2[$j]}};
					$max_2 = scalar @{$ortho_groups{$temp2[$j]}};
				}
				else{
					if($min_2 > scalar @{$ortho_groups{$temp2[$j]}}){
						$min_2 = scalar @{$ortho_groups{$temp2[$j]}};
					}
					if($max_2 < scalar @{$ortho_groups{$temp2[$j]}}){
						$max_2 = scalar @{$ortho_groups{$temp2[$j]}};
					}
				}
			}	
			else{
				$min_2=0;
			}
		}
		if($min_1 > $max_2 || $max_1 < $min_2){
			$done{$split[0]}=1;
			if($min_1 > $max_2){
				$count_3_up++;
				print OUT3 join ("\t", "up", $split[0], join ("", @results), @counts, @genes), "\n";
			}
			elsif($max_1 < $min_2){
				$count_3_down++;
				print OUT3 join ("\t", "down", $split[0], join ("", @results), @counts, @genes), "\n";
			}
		}
	}
	if(!$done{$split[0]}){
		$min_1=0;
		$min_2=0;
		$max_1=0;
		$max_2=0;
		@temp1=();
		@temp2=();
		if($ortho_groups{"KoPC"}){
			$min_2=scalar @{$ortho_groups{"KoPC"}}	;
			$max_2=scalar @{$ortho_groups{"KoPC"}}  ;
		}
		@temp2=("Shar", "Meug", "Mdom");
		for my $j (0 ..  $#temp2){
			my $k=0;
			if($ortho_groups{$temp2[$j]}){
				if($j==0){
					$min_1 = scalar @{$ortho_groups{$temp2[$j]}};
				}
				else{
					if($min_1 > scalar @{$ortho_groups{$temp2[$j]}}){
						$min_1 = scalar @{$ortho_groups{$temp2[$j]}};
					}
				}
				if($max_1 < scalar @{$ortho_groups{$temp2[$j]}}){
					$max_1 = scalar @{$ortho_groups{$temp2[$j]}};
				}
			}	
			else{
				$min_1=0;
			}
		}
		if(($min_2 > $max_1)||($max_2 < $min_1)){
			if($min_2 > $max_1){
				$count_4_up++;
				print OUT1 join ("\t", "up", $split[0], join ("", @results), @counts, @genes), "\n";	
			}
			elsif($max_2 < $min_1){
				$count_4_down++;
				print OUT1 join ("\t", "down", $split[0], join ("", @results), @counts, @genes), "\n";
			}
		}
		@temp2=();
	}
	print OUT join ("\t", $split[0], join ("", @results), @counts, @genes), "\n";
	@results=();
	@counts=();
	@genes=();
	%ortho_groups=();
}
close IN;
close OUT;
close OUT1;
close OUT2;
close OUT3;
close OUT4;

print "NUMBER OF GENES\n";
print join ("\t", "Expansions/contractions in Marsupials:", $count_1_up, $count_1_down),"\n";
print join ("\t", "Expansions/contractions in Australidelphia:", $count_2_up, $count_2_down),"\n";
print join ("\t", "Expansions/contractions in Diprotodontia:", $count_3_up, $count_3_down),"\n";
print join ("\t", "Expansions/contractions in Koala:", $count_4_up, $count_4_down),"\n";

#########################################################
#########################################################
sub CREATE_LOG{
	my ($array)=(@_);
	my $date=localtime();
	if(-e "$project"){
		open ("IN", "$project")||die"IN $project\n";
		my @project_log=<IN>;
		close IN;
		open(LOG, ">$project");
		for my $n (0 .. $#project_log){
			print LOG $project_log[$n];
		}
	}
	else{
		open(LOG, ">$project");
	}
	print LOG "\n\n#########################################################\n#########################################################\n
	$date\n
#########################################################\n#########################################################\n
	\n";
	print LOG "$0\n";
	for my $n (1 .. $#{$array}){
		print LOG join ("\t", $array->[$n]),"\n";
	}
	close LOG;
}

