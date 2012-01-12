#!/usr/bin/perl
#--------------------------------------------------
# This file is part of Forage.
# Copyright 2011 Malte Petersen <mptrsen@uni-bonn.de>
# 
# Forage is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
# 
# Forage is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License along with
# Forage. If not, see http://www.gnu.org/licenses/.
#-------------------------------------------------- 

use strict;   # make me write good code
use warnings; # cry if something seems odd
use Data::Dumper;
use Config;		# allows checking for system configuration
use Getopt::Long;
use Carp;		# alternative warn and die
use File::Path qw(make_path);	# mkdir with parent dirs; this also uses File::Spec
use File::Basename;	# parsing path names
use IO::File;	# object-oriented access to files
use IO::Dir;	# object-oriented access to dirs
use Tie::File;
(my $libdir = $0) =~ s/forage\.pl$//; 
use lib qw($libdir);
#use Forage::Item;
#use Genetic::Codes;

#--------------------------------------------------
# Only use threads if the system supports it.
# The whole threads system is totally not implemented yet,
#	do not attempt to use it!
#-------------------------------------------------- 

my $use_threads = 0;
if ($use_threads == 1 and $Config{'useithreads'}) {
	print "Using threads.\n";
	use Forage::Threaded;
}
elsif ($use_threads == 1 and !$Config{'useithreads'}) {
	die "Fatal: Cannot use threads: Your version of Perl was not compiled with threading support. Not using threads.\n";
}
else {
	print "Not using threads. ";
	use Forage::Unthreaded;
	print "Loaded Forage::Unthreaded.\n";
	$use_threads = 0;
}

#--------------------------------------------------
# # Variable initialisation
#-------------------------------------------------- 
my $version = 0.00001;#{{{

#--------------------------------------------------
# # Programs
#-------------------------------------------------- 
my $translateprog = '/home/mpetersen/src/exonerate-2.2.0/src/util/fastatranslate';
my $indexprog = 'fastaindex';
my $hmmsearchprog = 'hmmsearch';

#--------------------------------------------------
# # Other variables
#-------------------------------------------------- 
my $estfile = '';
my $protfile = '';
my $indexfile = '';
my $backup_ext = '.bak';
my $eval_threshold;
my $hmmdir = '';
my $hmmfile = '';
my $hmmfullout = 0;
my $hmmoutdir = $hmmsearchprog;
my $hmmoutopt;
my $hmmresultfileref;
my $outdir;
my $score_threshold;
my $verbose = 0;
my @eval_option = ();
my @hmmfiles;
my @hmmsearchcmd;
my @score_option = ();
my @seqobjs;
my $hitcount;
my $i;
my $header = <<EOF;
Forage: Find Orthologs using Reciprocity Among Genes and ESTs
Copyright 2011 Malte Petersen <mptrsen\@uni-bonn.de>
Version $version

EOF

#--------------------------------------------------
# # Get command line options
#-------------------------------------------------- 
GetOptions(	'v'         => \$verbose,
			'threads'         => \$use_threads,		# make using threads optional
			'estfile=s'       => \$estfile,
			'E=s'             => \$estfile,
			'eval=s'          => \$eval_threshold,
			'score=s'         => \$score_threshold,
			'H=s'             => \$hmmfile,
			'hmmdir=s'        => \$hmmdir,
			'hmmsearchprog=s'	=> \$hmmsearchprog,
			'hmmfullout'      => \$hmmfullout,
);#}}}

#--------------------------------------------------
# # Input error checking, reporting etc
#-------------------------------------------------- 
&intro;

#--------------------------------------------------
# # create list of HMM files
#-------------------------------------------------- 
@hmmfiles = &hmmlist;

print "Using HMM dir $hmmdir with ", scalar @hmmfiles, " HMMs\n" 
	if $hmmdir;
print "Using HMM file $hmmfile.\n" 
	if $hmmfile;

print "e-Value cutoff: $eval_threshold.\n" 
	if $eval_threshold;

print "Score cutoff: $score_threshold.\n"
	if $score_threshold;

#--------------------------------------------------
# # translate the ESTs to protein, index 
#-------------------------------------------------- 
print "Translating $estfile in all six reading frames...\t";
$protfile = &translate_est(File::Spec->catfile($estfile));
print "\n";
#--------------------------------------------------
# print "Indexing $estfile for fast access...\t";
# $indexfile = &index_est(File::Spec->catfile($protfile));
# print "\n";
#-------------------------------------------------- 

#--------------------------------------------------
# # hmmsearch the protfile using all HMMs
#-------------------------------------------------- 
print "Hmmsearching the protein file using all HMMs in $hmmdir...\n";
$i = 0;
$hitcount = 0;

#--------------------------------------------------
# # Setup the Forage module
#-------------------------------------------------- 

# verbose output; this is a class method
Forage::Unthreaded->verbose(1) if $verbose;

# the output directory
Forage::Unthreaded->hmmoutdir($hmmoutdir);

# the hmmsearch program
# Forage::Unthreaded->hmmsearchprog($hmmsearchprog);

# the hmmsearch command
Forage::Unthreaded->hmmsearchcmd(\@hmmsearchcmd);

# whether or not we want full output
Forage::Unthreaded->hmmfullout(0);

#--------------------------------------------------
# # Do the pHMM search
#-------------------------------------------------- 

foreach my $hmmfile (@hmmfiles) {
	++$i;
	# create new hmmobject with a hmm file, should have all the necessary info for doing hmmsearch
	my $hmmobj = Forage::Unthreaded->new($hmmfile);	
	# now do the hmmsearch on the protfile
	$hmmobj->hmmsearch($protfile);
	# count the hmmsearch hits
	unless ($hmmobj->hmmhitcount()) {	# do not care further with HMM files that did not return any result
		printf "%4d hits detected for %s\n", $hmmobj->hmmhitcount, basename($hmmobj->hmmfile) if $verbose;
		next;
	}
	printf "%4d hits detected for %s\n", $hmmobj->hmmhitcount, basename($hmmobj->hmmfile) if $verbose;
	++$hitcount;
	
	#--------------------------------------------------
	# #TODO next:
	#-------------------------------------------------- 
	# find the hit seqs in the EST file and re-search them against the core ortholog db
	# using either blast or hmmsearch [hmmsearch: don't we need profiles for that?]

	#--------------------------------------------------
	# # TODO then:
	#-------------------------------------------------- 
	# for the re-hits, gather nuc seq and compile everything that Karen wants output :)
}

printf "%d HMMs hit something.   %d HMM files processed. \n", $hitcount, $i;
print "Done!\n";
exit;


# parse the protfile
#--------------------------------------------------
# &parse_protfile($protfile, @hmmfiles);
#-------------------------------------------------- 

###################################################
# # Functions follow
###################################################

# Sub: intro
# Checks input, file/dir presence, etc.
# Returns True if everything is OK.
sub intro {#{{{
	die "Fatal: At least two arguments required: EST file (-E) and HMM file/dir (-H or -hmmdir)!\n"
		unless ($estfile and ($hmmdir or $hmmfile));

	die "Fatal: Can't use both e-value and score thresholds!\n"
		if ($eval_threshold and $score_threshold);

	$outdir = File::Spec->catdir('out_'.basename($estfile, '.fa'));
	$hmmoutdir = File::Spec->catdir($outdir, $hmmsearchprog);

	# build hmmsearch command line
	# full output only if desired, table otherwise
	if ($eval_threshold) {
		@eval_option = ('-E', $eval_threshold);
	}
	if ($score_threshold) {
		@score_option = ('-T', $score_threshold);
	}
	if ($hmmfullout) {
		@hmmsearchcmd = ($hmmsearchprog, @eval_option, @score_option, '-o');
	}
	else {
		@hmmsearchcmd = ($hmmsearchprog, @eval_option, @score_option, '--tblout');
	}

	print $header;

	print "EST file " and &doesexist($estfile);

	print "HMMsearch output dir " and &createdir($hmmoutdir);
}#}}}

# Sub: hmmlist
# Expects: nothing (?)
# Returns: array hmmfiles
sub hmmlist {#{{{
	if ($hmmdir) {
		my $dir = IO::Dir->new(File::Spec->catdir($hmmdir));
		while (my $file = $dir->read) {
			push(@hmmfiles, File::Spec->catfile($hmmdir, $file)) if ($file =~ /\.hmm$/);
		}
		$dir->close;
		return sort @hmmfiles;
	}
	else {
		print "single hmm\n"
			if $verbose;
		push(@hmmfiles, $hmmfile);
		return @hmmfiles;
	}
}#}}}

# Sub: translate_est
# Translate a nucleotide fasta file to protein in all six reading frames
# Expects: scalar string filename
# Returns: scalar string filename (protfile)
sub translate_est {#{{{
	my ($infile) = shift;
	(my $outfile = $infile) =~ s/(\.fa$)/_prot$1/;
	if (-e $outfile) {
		print "$outfile exists, using this one";
		return $outfile;
	}
	&backup_old_output_files($outfile);
	my $translateline = $translateprog . " " . $infile . ">$outfile";
	die "Fatal: Could not translate $infile: $!\n"
		if system($translateline);

	return $outfile;
}#}}}

# Sub: index_est
# index a fasta file for fast fetching of sequences using fastafetch
# Expects: scalar string filename
# Returns: scalar string filename (indexfile)
sub index_est {#{{{
	my ($infile) = shift;
	(my $indexfile = $infile) =~ s/\.fa$/.idx/;
	if (-e $indexfile) {
		print "$indexfile exists, using this one.\n";
		return $indexfile;
	}
	&backup_old_output_files($indexfile);
	my @indexline = ($indexprog, '--index', $indexfile, '--fasta', $infile);
	die "Fatal: Could not index $infile: $!\n"
		if system(@indexline);
}#}}}

# Sub: gethmmscores
# Parse the hmmsearch result, populate the data structure
# Expects: scalar reference to hmmresult file
# Returns: reference to result object
sub gethmmscores {#{{{
	my $hmmresultref = shift;
}#}}}

# sub: clean_up_old_output_files
# input: reference to list of relevant contigs
sub backup_old_output_files {#{{{
	my ($outfile) = shift @_;
	if (-e $outfile) {
		rename( $outfile, $outfile.$backup_ext ) or die "wah: could not rename $outfile during backup: $!\n";
		print "backed up old file $outfile to ", $outfile.$backup_ext, "\n";
	}
}#}}}

# Sub: helpmessage
# Prints the help message
sub helpmessage {#{{{
	my $helpmessage = <<EOF;
'estfile=s' 
'E=s'				
'H=s'				
'hmmdir=s'	
EOF
	print $helpmessage;
}#}}}

# Sub: createdir
# Creates a directory with parent directories if it does not already exist
# Expects: scalar string dirname
# Returns: True if successful
sub createdir {#{{{
	my $dir = shift;
	unless (-d $dir) {
		if (-e $dir and not -d $dir) {
			die "Warning: $dir exists, but is not a directory! This will most likely lead to problems later.\n";
		}
		else {
			print "does not exist, creating $dir... ";
			make_path $dir, { verbose => 0 } or die "Fatal: Could not create $dir: $!\n";
			print "OK\n";
			return 1;
		}
	}
	print "$dir exists.\n";
	return 1;
}#}}}

# Sub: doesexist
# Checks whether an item exists on disk
# Expects: scalar string filename/dirname
# Returns: True if exists, False otherwise
sub doesexist {#{{{
	my $item = shift;
	unless (-e $item) {
		warn "$item does not exist!\n";
		return 0;
	}
	print "$item exists.\n";
	return 1;
}#}}}
