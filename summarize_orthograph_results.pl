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
# Version: 2014-08-28

use strict;
use warnings;
use File::Spec::Functions;
use File::Basename;
use Getopt::Long;
use Tie::File;
use Carp;


my $header_sector_separator = '&&';

# Get options from command line
my %options;
GetOptions( \%options, 'i=s', 'o=s' , 'c', 'm', 'd=s', 'help', 'h' );

my $usage =
	"\nUSAGE: summarize_orthograph_results.pl -i INPUTFOLDER -o OUTPUTFOLDER\n\n"
   ."Mandatory parameters: -i FOLDER with Orthograph result folders\n"
    ."                      -o existing FOLDER for summary files to be stored\n"
   ."\nOptional parameters:  -c Make corresponding AA and NT headers identical\n"
   ."                      -d FILENAME of file with a list of species you want\n"
   ."                         to exclude. Sequences of these species will not\n"
   ."                         be in the output. Species names have to be in\n"
   ."                         separate lines.\n"
   ."                      -m Mask stop symbols (*) with X\n\n";
;

# Check user input
if (    $options{ 'help' }
     or $options{ 'h' } 
     or !exists $options{ 'i' }
     or !exists $options{ 'o' }
    ) {

    print $usage;
    exit;

}

if ( $options{ 'i' } eq $options{ 'o' } ) {

    print "INPUTFOLDER and OUTPUTFOLDER must be different!\n";
    exit;

}

if( ! -e $options{ 'i' } ) {
	print "INPUTFOLDER \"$options{ 'i' }\" does not exist!\n";
	exit;
}


if( ! -e $options{ 'o' } ) {
	print "OUTPUTFOLDER \"$options{ 'o' }\" does not exist!\n";
	exit;
}


my %is_to_discard;
if ( exists $options{ 'd' } ) {

    print $usage if !defined $options{ 'd' }; 

	open( my $fh, '<', $options{ 'd' } )
		or die "Could not open file \"$options{ 'd' }\": $!\n";

	while( my $line = <$fh> ) {
		chomp $line;
		next if $line =~ m/^\s*$/;
		$is_to_discard{$line} = 1;
	}
}

my $corresponding = exists $options{ 'c' } ? 1 : 0;
my $mask_stop_symbols = exists $options{ 'm' } ? 1 : 0;

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
	print "Processing amino acids\n";

	foreach my $aa_file_name ( @aa_file_names ) {

	    # Name of the summary file
	    $aa_file_name =~ m/^(.+)\.fa$/i;
	    my $summary_filename_aa = $1.'.summarized.fa';

	    # Full path for inputfile and outputfile without coreorthologs
	    my $input_file  = catfile( $options{ 'i' }, $dir, $aa, $aa_file_name );
	    my $output_file = catfile( $options{ 'o' }, $aa_sum_dir, $summary_filename_aa );

	    &summarize( $input_file, $output_file, $header_sector_separator, $corresponding,
	                \%is_to_discard, $mask_stop_symbols, 
	               )
	    ;	
    }
        
    # Process all nt files of current taxon
	print "Processing nucleotides\n";
	
    foreach my $nt_file_name ( @nt_file_names ) {

		# Name of summary file
		$nt_file_name =~ m/^(.+)\.fa$/i;
		my $summary_filename_nt = $1.'.summarized.fa';
			
		# Full path for input and outputfile
		my $input_file  = catfile( $options{ 'i' }, $dir, $nt, $nt_file_name );
		my $output_file = catfile( $options{ 'o' }, $nt_sum_dir, $summary_filename_nt );

	    &summarize( $input_file, $output_file, $header_sector_separator, $corresponding,
	                \%is_to_discard, $mask_stop_symbols, 
	               )
	    ;	
	}	
}

##########################################################################################

sub summarize {
	
	# Unpack @_
	my ( $input_file, $output_file, $header_sector_separator, $corresponding,
	     $is_to_discard_REF, $mask_stop_symbols                                 ) = @_;
		
	open( my $input_file_fh, '<', $input_file )
		or die "Could not open file\"$input_file\": $!\n";
	
	my $file_exists = ( -e "$output_file" ) ?  1  :  0;

	open( my $output_file_fh, '>>', $output_file )
		or die "Couldn't open file \" $output_file\": $!\n";

	my $print_sequences = 1;
	
	while( my $line = <$input_file_fh> ) {
		chomp $line;
		
		next if $line =~ m/^\s+$/;
		
		# Headers
		if( $line =~ m/^>/ ) {
			my @header_sectors       = split( $header_sector_separator, $line );
			my @columns_first_sector = split( '\|', $header_sectors[0] );
			
			# Subsequent file of a given ortholog group
			if( $file_exists ) {

				# Orthograph header of reference taxon 
				if(     scalar @columns_first_sector == 6 
				    and $columns_first_sector[5] eq '.' ) {
				    
					$print_sequences = 0;
				
				}
				# HaMStR header of reference taxon 
				elsif ( scalar @columns_first_sector == 3 ) {
				
					$print_sequences = 0;

				}				
				# Header of non-reference taxon
				else {
					
					$line = &erase_range_and_translation_information(
					        $header_sector_separator, @header_sectors )
						if $corresponding;
					print { $output_file_fh } $line, "\n";
					$print_sequences = 1;
					
				}
			}
			# First file of a given ortholog group
			else {
				# Orthograph header of reference taxon 
				if(     scalar @columns_first_sector == 6
				    and $columns_first_sector[5] eq '.' ) {
				    
					my $taxon = $columns_first_sector[1];
					
					if( exists ${$is_to_discard_REF}{$taxon} ) {
						$print_sequences = 0;
						next;
					}
					
					$line = &erase_range_and_translation_information(
					        $header_sector_separator, @header_sectors )
						if $corresponding;
					print { $output_file_fh } $line, "\n";
					$print_sequences = 1;
					
				}
				# HaMStR header of reference taxon 
				elsif ( scalar @columns_first_sector == 3 ) {
					my $taxon = $columns_first_sector[1];
					if( exists ${$is_to_discard_REF}{$taxon} ) {
						$print_sequences = 0;
						next;
					}
					print { $output_file_fh } $line, "\n";
					$print_sequences = 1;
				}
				# Header of non-reference taxon
				else {
					$line = &erase_range_and_translation_information(
					        $header_sector_separator, @header_sectors )
						if $corresponding;
					print { $output_file_fh } $line, "\n";
					$print_sequences = 1;
				}
			}
		}
		
		# Sequences
		else{
			$line =~ tr/*/X/ if $mask_stop_symbols;
			print { $output_file_fh } $line, "\n" if $print_sequences;
		
		}
	}
	close $input_file_fh  or die "Could not close file\"$input_file\": $!\n";
	close $output_file_fh or die "Couldn't close file \" $output_file\": $!\n";
}


sub erase_range_and_translation_information {

	my ( $header_sector_separator, @header_sectors ) = @_;

	my @new_header_sectors;
	my $first = 1;
	while( my $header_sector = shift @header_sectors ) {
		my @columns = split( '\|', $header_sector );
		if ( $first ) {
			( $columns[3], $columns[4] ) = ('.', '.');
			$first = 0;
		}
		else {
			( $columns[1], $columns[2] ) = ('.', '.');
		}
		push( @new_header_sectors, join( '|',  @columns ) );
	}
	return join( $header_sector_separator, @new_header_sectors ); 

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