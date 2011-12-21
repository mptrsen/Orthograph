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
use File::Path qw(make_path);	# mkdir with parent dirs
use File::Basename;	# parsing path names
use File::Spec;
use Tie::File;
(my $libdir = $0) =~ s/forage\.pl$//; 
use lib qw($libdir);
#use Forage::Item;
#use Genetic::Codes;

#--------------------------------------------------
# only use threads if the system supports it
# the whole threads system is totally not implemented yet,
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
my $translatecmd = 'fastatranslate';
my $hmmsearchprog = 'hmmsearch';

my $estfile = '';
my $eval_threshold;
my $hmmdir = '';
my $hmmervars;
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
GetOptions(	'v'					=> \$verbose,
			'threads'			=> \$use_threads,		# make using threads optional
			'estfile=s'			=> \$estfile,
			'E=s'				=> \$estfile,
			'eval=s'			=> \$eval_threshold,
			'score=s'			=> \$score_threshold,
			'H=s'				=> \$hmmfile,
			'hmmdir=s'			=> \$hmmdir,
			'hmmsearchprog=s'	=> \$hmmsearchprog,
			'hmmfullout'		=> \$hmmfullout,
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
# # translate the ESTs to protein
#-------------------------------------------------- 
print "Translating EST file to protein... ";
my $protfile = &translate_est(File::Spec->catfile($estfile));
tie(my @targets, 'Tie::File', $estfile) or die "Fatal: Could not tie to $estfile\: $!\n";

#--------------------------------------------------
# # hmmsearch the protfile using all HMMs
#-------------------------------------------------- 
print "Hmmsearching the protein file using all HMMs in $hmmdir...\n";
$i = 0;
$hitcount = 0;

#--------------------------------------------------
# # Create HMMer variable package
#-------------------------------------------------- 
$hmmervars = {
	'hmmfile'				=> \$hmmfile,
	'protfile'			=> \$protfile,
	'outdir'				=> \$outdir,
	'hmmfullout'		=> \$hmmfullout,
	'hmmsearchprog' => \$hmmsearchprog,
	'hmmsearchcmd'	=> \@hmmsearchcmd,
	'verbose'				=> \$verbose,
};

foreach my $hmmfile (@hmmfiles) {
	++$i;
	$hmmervars->{'hmmfile'} = \$hmmfile;
	$hmmervars->{'protfile'} = \$protfile;
	# create new hmmobject, has all the necessary info for doing hmmsearch
	my $hmmobj = Forage::Unthreaded->new($hmmervars);	
	# now do the hmmsearch
	$hmmobj->hmmsearch($hmmervars);
	# count the hmmsearch hits
	if ($hmmobj->hitcount == 0) {
		print "No hits detected\n";
		next;
	}
	print 'Hits detected: ' . $hmmobj->hitcount . "\n" if $verbose;
	foreach ($hmmobj->hits) {
		print $_ . "\n";
	}
	++$hitcount;
	
	#--------------------------------------------------
	# # TODO next:
	#-------------------------------------------------- 
	# find the hit seqs in the EST file and re-search them against the core ortholog db
	# using either blast or hmmsearch [hmmsearch: don't we need profiles for that?]
	# use a pre-prepared Tie::File array, faster for large files
	# TODO do this the object-oriented way :P
	print "found back in EST file: \n" if grep( $hmmobj->{'hmmhits'}, @targets );

	#--------------------------------------------------
	# # TODO then:
	#-------------------------------------------------- 
	# for the re-hits, gather nuc seq and compile everything that Karen wants output :)
}

printf "%d HMM files processed.\n", $i;
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

	$outdir = File::Spec->catdir('out_'.basename($estfile));
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
		my $dir = File::Spec->catdir($hmmdir);
		opendir(my $dir_handle, $dir) or die "Fatal: Could not open HMM dir: $!\n";
		#opendir(my $dir, $hmmdir) or die "Fatal: Could not open HMM dir: $!\n";
		while (my $file = readdir $dir_handle) {
			push(@hmmfiles, File::Spec->catfile($dir, $file)) if ($file =~ /\.hmm$/);
		}
		closedir($dir_handle);
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
	my ($infile) = @_;
	(my $outfile = $infile) =~ s/(\.fa$)/_prot$1/;
	if (-e $outfile) {
		print "$outfile exists, using this one.\n";
		return $outfile;
	}
	&backup_old_output_files($outfile);
	my $translateline = $translatecmd . " " . $infile . ">$outfile";
	print "translating $infile, please wait...\t";
	die "Fatal: Could not translate $infile: $!\n"
		if system($translateline);
	print "done!\n";
	return $outfile;
}#}}}

# Sub: gethmmscores
# Parse the hmmsearch result, populate the data structure
# Expects: scalar reference to hmmresult file
# Returns: reference to result object
sub gethmmscores {#{{{
	my $hmmresultref = shift;
	print "Tieing hmmresult file $hmmresultref to an array\n"
		if $verbose;
	tie(my @file, 'Tie::File', $hmmresultref) 
		or die "Fatal: Could not open HMM result file: $!\n";
	untie @file 
		or die "Warning: Could not close HMM result file: $!\n";
	print "$hmmresultref successfully untied\n"
		if $verbose;
}#}}}

# sub: clean_up_old_output_files
# input: reference to list of relevant contigs
sub backup_old_output_files {#{{{
	my ($outfile) = shift @_;
	my $backup_ext = '.bak';
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
