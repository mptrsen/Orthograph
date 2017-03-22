#!/usr/bin/env perl

###############################################################################################
#Calculate statistics of Orthograph results for one or more species of the analysis
###############################################################################################
#
#Calculate no. of hits (no. assigned orthologous protein sequences for each species)
#absolute total no. of aa sites (including ambiguous sites) for all the hits of each species, 
#number of ambiguous sites (X and *), N50 statistic for protein hits,
#average, maximum, minimum and median length of protein hits
#
#The script must be executed in the directory where the result folders (1 or more species)
#of Orthograph are placed -> eg /home/user/Documents/Orthograph_results. Results are printed
#in a tab delimited (excel format) .txt file in the same directory
#
#Usage: calculate_statistics_orthograph_results.pl species_directory_1 [species_directory_2 ..]
#or calculate_statistics_orthograph_results.pl list_with_names_of_species_directories.txt
#Copyright 2017 A.Vasilikopoulos
###############################################################################################

use strict;
use warnings;
use Cwd;

die "Usage: $0 species_directory_1 [species_directory_2 ..]\nAlternatively: $0 list_of_taxa.txt\n" unless (scalar @ARGV);

my $pwd = cwd();

#If one argument is given and is a .txt file (list of names)
if (scalar @ARGV == 1 and -f $ARGV[0]) {
	
	my $list=$ARGV[0];	
	open (my $fh_list, '<', $list)
		or die"Could not open file \"$list\":$!\n";
	
	my @taxa = <$fh_list>;
	chomp @taxa;
	
	my $out = "statistics_orthograph.txt";
	open (my $fh_out, '>', $out)  
    or die "Could not open file \"$out\":$!\n"; 
    
	#print labels        
	print {$fh_out} "Species_name\tno.hits\ttotal_no_aa\tno_X\tno_stop\tN50\taverage_length\tmedian_length\tmax_length\tmin_length\n";
	
	foreach (@taxa) {
	
		chdir "$pwd";
	
		die "$_ is either not present or not a directory" unless -d $_;
	
		#return statistics for dir and print
		my @statistics = &calculate_statistics_orthograph_results($_);
	
		printf {$fh_out} "%s\t%d\t%d\t%d\t%d\t%d\t%d\t%.1f\t%d\t%d\n", $_,
			$statistics[0],
			$statistics[1],
			$statistics[2],
			$statistics[3],
			$statistics[4],
			$statistics[5],
			$statistics[6],
			$statistics[7],
			$statistics[8];
	}	
}

#If multiple directories are given as arguments
else {
	
	my $out = "statistics_orthograph.txt";
	open (my $fh_out, '>', $out)  
		or die "Could not open file \"$out\":$!\n"; 
  
	#print labels          
	print {$fh_out} "Species_name\tno.hits\ttotal_no_aa\tno_X\tno_stop\tN50\taverage_length\tmedian_length\tmax_length\tmin_length\n";
		 
	foreach (@ARGV) {
	
		chdir "$pwd";
	
		die "$_ is either not present or not a directory" unless -d $_;
	
		#return statistics for dir and print
		my @statistics = &calculate_statistics_orthograph_results($_);
	
		printf {$fh_out} "%s\t%d\t%d\t%d\t%d\t%d\t%d\t%.1f\t%d\t%d\n",
			$_,
			$statistics[0],
			$statistics[1],
			$statistics[2],
			$statistics[3],
			$statistics[4],
			$statistics[5],
			$statistics[6],
			$statistics[7],
			$statistics[8];
	}                    
}

#############################################################################################################################################

sub calculate_statistics_orthograph_results {
	
	my $input_dir = shift @_;
	
	#generate statistics variables
	my @lengths;
	my $no_hits;
	my $total_aa_sites=0;
	my $total_Xs=0;
	my $total_stop=0;
	my $N50=0;
	my $mean=0;
	my $median_length=0;
	my $max_length=0;
	my $min_length=0;
	
	#Change to the directory where the COG fasta files with the aminoacid sequeces are placed
	my $path_aa="$input_dir/aa";	
	if (-d $path_aa) {chdir("$path_aa");}	   
	else{ die"Could not access directory $path_aa:$!\n";}
		
	#Count no. of hits for the species
	$no_hits = qx(ls -1 | wc -l);	
	chomp $no_hits;
    
	#Calculate total aa sites, total X sites and total stop codons
    my @COG_names = glob "*.fa";
    if (@COG_names){
	
	     #loop through the output aa files for each species
	     foreach my $COG (@COG_names) {
		
		    open (my $fh_COG, '<', $COG) 
			or die"Could not open file \"$COG\":$!\n";
		
		     while (my $line = <$fh_COG>) {			
			   chomp $line;
						
			   if ($line =~ m/>.+\[translate\([1-3]\)\].+/) {				 
				   $line = <$fh_COG>;	
				   chomp $line;				 
			 
				   #Count no of ambiguous sites (* and X)
				   my @aa = split ("", $line);				
				   foreach my $aa (@aa){
					   if    ($aa eq "*") { ++$total_stop }
					   elsif ($aa eq "X") { ++$total_Xs   }
				   }
				   #calculate total aa length
				   my $length = length $line;			 
				   $total_aa_sites += $length;   
				   push (@lengths, $length);								 
			   }
		     }		
	    	close $fh_COG;
	      }	
	
	     #calculate average length of protein hits
	     $mean = $total_aa_sites/$no_hits;

	     #Calculate N50 for protein hit lengths
	     my @lengths_sorted = sort { $b <=> $a } @lengths;    
	     my $N50_threshold = $total_aa_sites/2;   
	     my $sum_check = 0;
    
	     for (my $i = 0; $i < @lengths_sorted; ++$i){
		     $sum_check += $lengths_sorted[$i];		
		     if ($sum_check >= $N50_threshold) {
			     $N50 = $lengths_sorted[$i];last;
		     }
	     }

	     #Calculate median length of protein hits
	     if ($no_hits%2==0) {
		     my $first = $lengths_sorted[scalar @lengths_sorted/2];		
		     my $second = $lengths_sorted[(scalar @lengths_sorted/2)-1];		
		     $median_length = ($first+$second)/2;
	     }
	     else {
		     $median_length = $lengths_sorted[((scalar @lengths_sorted)-1)/2];
	     }
	     #calculate max, min lengths
	     $max_length = shift @lengths_sorted;
	     $min_length = pop @lengths_sorted;
    }
	my @statistics = ($no_hits, $total_aa_sites, $total_Xs, $total_stop, $N50, $mean, $median_length, $max_length, $min_length);
	return @statistics;
}
