#!/usr/bin/env perl

# Copyright (C) 2014, Oliver Niehuis
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# Please contact: o.niehuis@zfmk.de
#
# The program summarizes the ouput of orthograph, Status: 2014-08-24
# The program work only with genuine Orthograph headers
# The program will fail summaring Orthograph output with HaMStR headers
# usage: perl summarize_orthograph_results.pl -id PATH -sd PATH
#	-i input directory: output directories of orthograph: aa, nt, log
#	-o -> provide output/summarize directory, must be outside out of the input directory.
# 	-> provide the absolute path for both!

use strict;
use warnings;
use File::Spec::Functions;
use File::Basename;
use Getopt::Long;
use Tie::File;
use Carp;

# Get options from command line
my %options;
GetOptions( \%options, 'i=s', 'o=s', 'help', 'h' );

# Check user input
if (    $options{ 'help' }
     or $options{ 'h' } 
     or !exists $options{ 'i' }
     or !exists $options{ 'o' } ) {

    print "USAGE: summarize_orthograph_results.pl -i INPUTDIRECTORY -o SUMMARYDIRECTORY\n";
    exit;

}
if ( $options{ 'i' } eq $options{ 'o' } ) {

    print "INPUTDIRECTORY and SUMMARYDIRECTORY must be different\n";
    exit;

}

# Names for summary directories
my $aa_sum_dir = 'aa_summarized';
my $nt_sum_dir = 'nt_summarized';

# Create directories with summary results
&check_summary_directory( catdir( $options{ 'o' }, $aa_sum_dir ) );
&check_summary_directory( catdir( $options{ 'o' }, $nt_sum_dir ) );

# Read all directories and files in input directory into an array
my @dir = &read_dir( $options{ 'i' }, { 'hide' => 1 } );

# Process files of each directory
foreach my $dir ( @dir ) {

	print "Processing directory: $dir\n";
	
	# Get all files from the species directory
	my @items = &read_dir( catfile( $options{ 'i' }, $dir ), { 'hide' => 1 } );

    # Names of subfolders for amino acids and nucleotides
    my $aa = 'aa';
    my $nt = 'nt';

	# Get all file names except system files
	my @aa_file_names = &read_dir( catdir( $options{ 'i' }, $dir, $aa ), {'hide' => 1} );
	my @nt_file_names = &read_dir( catdir( $options{ 'i' }, $dir, $nt ), {'hide' => 1} );

    # Process all aa files of current taxon
	foreach my $aa_file_name ( @aa_file_names ) {

	    # Name of the summary file
	    my $summary_filename_aa;
	    $aa_file_name =~ m/^(.+)(\.aa\.fa)$/i;
	    $summary_filename_aa = $1.'.summarized'.$2;

	    # Full path for inputfile and outputfile without coreorthologs
	    my $input_file  = catfile( $options{ 'i' }, $dir, $aa, $aa_file_name );
	    my $output_file = catfile( $options{ 'o' }, $aa_sum_dir, $summary_filename_aa );

	    # File does not exist
	    if ( ! -e "$output_file" ) {
		    # Save each line of the inputfile in the outpufile
		    &move_fasta_file( $input_file, $output_file );
	    }

	    # File does already exist
	    else {
	    
	    	&append_non_reference_sequences_from_to( $input_file, $output_file );	
	    
	    }
    }
    
    # Process all nt files of current taxon
    foreach my $nt_file_name ( @nt_file_names ) {

		# Name of summary file
		$nt_file_name =~ m/^(.*)(\.nt\.fa)$/i;
		my $summary_filename_nt = $1.'.summarized'.$2;
			
		# Full path for input and outputfile
		my $input_file  = catfile( $options{ 'i' }, $dir, $nt, $nt_file_name );
		my $output_file = catfile( $options{ 'o' }, $nt_sum_dir, $summary_filename_nt );

		# File does not exist
	    if ( ! -e "$output_file" ) {
		    # Save each line of the inputfile in the outpufile
		    &move_fasta_file( $input_file, $output_file );
	    }

	    # File does already exist
	    else {   
	    	
	    	&append_non_reference_sequences_from_to( $input_file, $output_file );
	    }
	}
}

##########################################################################################

sub append_non_reference_sequences_from_to {
	
	# Unpack @_
	my ( $input_file, $output_file ) = @_;
	
	open( my $input_file_fh, '<', $input_file )
		or die "Could not open file\"$input_file\": $!\n";
	
	open ( my $output_file_fh, '>>', $output_file )
		or die "Couldn't open file \" $output_file\": $!\n";

	
	my $is_reference_taxon = 0;
	while( my $line = <$input_file_fh> ) {
		chomp $line;
		if( $line =~ m/^>/ ) {
			my @columns = split( '\|', $line );
			if( $columns[1] eq '.' or scalar @columns == 3 ) {
				$is_reference_taxon = 1;
				next;
			}
			else{
				$is_reference_taxon = 0;
				print {$output_file_fh} $line, "\n";
			}
		}
		else{
				print {$output_file_fh} $line, "\n" if !$is_reference_taxon;
		
		}
	}
	close $input_file_fh or die "Could not close file\"$input_file\": $!\n";
	close $output_file_fh or die "Couldn't close file \" $output_file\": $!\n";
}


##########################################################################################

sub check_summary_directory {
	
	# Unpack @_
	my $sd = shift @_;
	
	if ( -e $sd ) {
	
		# If exists, check whether or not it contains files
		my @sum_dir = &read_dir( $sd, {'hide' => 1 } );
	
		# If it contains files, ask whether they should be deleted or exit program
		if ( @sum_dir )  {
		    die "Directory \"$sd\" is not empty! Please delete the directory's content ".
		        "or choose a different directory.\n";
		}
	}
	
	# If directory does not exist, create it
	else {
		mkdir $sd;
	}
}

##########################################################################################

sub read_dir {
		
	# Unpack @_
	my ( $path, $arg_ref ) = @_; 

	# Defaults
	my %DEFAULT_OF = ( 'hide' => 1,  # 1: hide system files; 0: don't 
	);

	# Check provided arguments
	die 'Missing or superfluous arguments'        if @_ < 1  || @_ > 2;
	die 'Option(s) not passed via anonymous hash' if @_ == 2 && ref $arg_ref ne 'HASH';

	foreach my $provided_options ( keys %{ $arg_ref } ) {
		die 'Unknown option(s)' if !exists $DEFAULT_OF{ $provided_options }; 
	}

	# Set defaults
	#          If option given...            use option             else default
	my $hide = exists $arg_ref->{'hide'}  ?  $arg_ref->{'hide'}  :  $DEFAULT_OF{'hide'};

	# Open directory handle
	opendir ( my $dir, $path ) or 
		die "Couldn't find path \"$path\": $!";

	# Read file names
	my @files = readdir( $dir ) or
		die "Couldn't read directory \"$path\": $!";

	# Close directory handle
	closedir ( $dir ) or
		die "Couldn't close directory \"$path\": $!";

	# Filter hidden system files out
	if ( $hide ) {
		@files = grep {! /^\./ } @files;
	}

	# Return file names
	return @files;
}

##########################################################################################

sub move_fasta_file {

	# Unpack @_
	my $inputfile  = shift @_;
	my $outputfile = shift @_;
	
	# Open filehandles for both files
	open ( my $FH , '<', $inputfile  ) or die "Couldn't open file \"$inputfile\": $!\n" ;
	open ( my $FH2, '>', $outputfile ) or die "Couldn't open file \"$outputfile\": $!\n";
	
	# Iterate through each line of the file and save it in the outputfile
	my $first_line = 1;
	while ( my $line = <$FH> ) {
        print {$FH2} $line;
    }
	
	# Close filehandles
	close $FH;
	close $FH2;

}
