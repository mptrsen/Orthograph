#!/usr/bin/perl
#--------------------------------------------------
# This file is part of Orthograph.
# Copyright 2014 Malte Petersen <mptrsen@uni-bonn.de>
# 
# Orthograph is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
# 
# Orthograph is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License along with
# Orthograph. If not, see http://www.gnu.org/licenses/.
#-------------------------------------------------- 
use strict;
use warnings;
use autodie;
use File::Temp;
use File::Basename;
use File::Spec;
use Getopt::Long;
use Data::Dumper;
use Cwd;

my $skipall         = 0;
my $outdir          = cwd();
my $exonerate       = 'exonerate';
my $score_threshold = 30;
my $num_threads     = 1;
my $help            = 0;

print "Called: $0 @ARGV\n";

my $usage  = <<"END_OF_USAGE";
Usage: $0 [OPTIONS] PEPTIDEFILE CDSFILE
Generates 100% corresponding amino acid and nucleotide sequences for OGS
Options:
  --help, -h          print this help
  --skipall           skip all sequences that could not be found in both peptide and cds file
  --outdir PATH       set output directory to PATH (default: '/tmp')
  --exonerate PATH    set path to Exonerate executable (default: 'exonerate')
  --score-threshold N set score threshold for the alignments (default: 30)
  --num-threads N     set number of CPU threads to use [requires the Perl module Parallel::ForkManager] (default: 1)
END_OF_USAGE

GetOptions(
	'skipall'           => \$skipall,
	'outdir=s'          => \$outdir,
	'exonerate=s'       => \$exonerate,
	'score-threshold=i' => \$score_threshold,
	'num-threads=i'     => \$num_threads,
	'help|h'            => \$help,
) or print $usage and exit 1;

# see whether Parallel::ForkManager is installed
my $forkmanager = eval {
	require Parallel::ForkManager;
	Parallel::ForkManager->new( $num_threads );
};
unless ($forkmanager) {
	if ($num_threads != 1) {
		print "Error: Requested to multi-thread but Parallel::ForkManager not usable\n" and exit 1;	
	}
}

# need exactly two input files
scalar @ARGV == 2 or die $usage;

# read input files and check they have equal numbers of sequences
my $ogs         = slurp_fasta($ARGV[0]) or die "Fatal error, exiting\n";
my $transcripts = slurp_fasta($ARGV[1]) or die "Fatal error, exiting\n";
unless (scalar keys %$ogs == scalar keys %$transcripts) {
	warn "Warning: Unequal number of sequences!\n";
}

# initialize more variables
my @strange              = ();
my @doublestrange        = ();
my $nupd                 = 0;
my $nunchgd              = 0;
my $n                    = 0;
my $nseqs                = scalar keys %$ogs;
my $c                    = 0;
my $new_ogs_file         = File::Spec->catfile($outdir, 'corresp-' . basename($ARGV[0]));
my $new_transcripts_file = File::Spec->catfile($outdir, 'corresp-' . basename($ARGV[1]));

# open output files
open my $new_ogs, '>', $new_ogs_file;
open my $new_transcripts, '>', $new_transcripts_file;

# go through the peptide sequences sorted by header
foreach my $hdr (sort {$a cmp $b} keys %$ogs) {

	++$c;

	# fork this
	if ($forkmanager) { $forkmanager->start and next }

	# if the headers aren't identical, the sequences can't be correlated. 
	# in this case, they may be skipped at the user's discretion (--skipall).
	# otherwise, the program will exit.
	unless (exists $transcripts->{$hdr}) {
		if ($skipall) {
			print "Not found in transcriptome: '$hdr', skipping\n";
			next;
		}
		else {
			print "Not found in transcriptome: '$hdr', exiting\n" and exit(1);
		}
	}

	# generate input files 
	my $aafn    = fastaify($hdr, $ogs->{$hdr});
	my $ntfn    = fastaify($hdr, $transcripts->{$hdr});
	my $outfile = File::Temp->new();

	# some settings for exonerate
	my $exonerate_model = 'protein2dna';
	my $exonerate_ryo   = '>nucleotide\n%tcs>aminoacid\n%qas';

	my $reportline = "$hdr ($c of $nseqs) ";
	print $reportline . "aligning...\n";

	# run exonerate
	my @command = QW( "$exonerate
		--bestn 1
		--score $score_threshold
		--ryo '$exonerate_ryo'
		--model $exonerate_model
		--querytype protein
		--targettype dna
		--verbose 0
		--showalignment no
		--showvulgar no
		--query $aafn
		--target $ntfn
		> $outfile"
	);
	system("@command") and die "Fatal: Exonerate failed with errcode $?: $!\n";

	# get the alignment results
	my $res = slurp_fasta($outfile);
	unless ($res) {
		print $reportline . "strange sequences, skipping\n";
		push @doublestrange, $hdr;
		next;
	}

	# no alignment found, something is wrong, skip this pair
	if (!$res->{'nucleotide'}) {
		print $reportline . "no alignment found, skipping\n";
		push @strange, $hdr;
		next;
	}
	# the sequence has been updated
	elsif ($res->{'nucleotide'} ne $transcripts->{$hdr}) {
		printf $new_ogs ">%s\n%s\n", $hdr, $res->{'aminoacid'};
		printf $new_transcripts ">%s\n%s\n", $hdr, $res->{'nucleotide'};
		# count
		++$nupd;
		$reportline .= "updated: ";
	}
	# the sequence is the same
	else {
		printf $new_ogs ">%s\n%s\n", $hdr, $res->{'aminoacid'};
		printf $new_transcripts ">%s\n%s\n", $hdr, $res->{'nucleotide'};
		# count
		++$nunchgd;
		$reportline .= "unchanged: ";
	}

	# check whether lengths correlate
	if (length($res->{'nucleotide'}) == length($res->{'aminoacid'}) * 3) {
		print $reportline . "lengths equal\n";
	}
	else {
		print $reportline . "lengths differ\n";
		push @strange, $hdr;
	}

	# set output buffer back to normal
	local $| = 1;

	# delete tempfiles
	undef $aafn;
	undef $ntfn;
	undef $outfile;

	# catch fork
	if ($forkmanager) { $forkmanager->finish }
}

# wait for all children
if ($forkmanager) { $forkmanager->wait_all_children }

print "Done\n";
printf "Updated %d of %d sequences, wrote %d sequence%s to %s and %s\n",
	$nupd,
	$c,
	$nupd + $nunchgd,
	$n == 1 ? '' : 's',
	$new_ogs_file,
	$new_transcripts_file,
;

if (@strange) {
	print "The following sequences are not of equal length or no alignment could be found. You should double-check these:\n";
	while (my $hdr = shift @strange) {
		print $hdr, "\n";
	}
	print "The most common reason are ambiguity characters (X) in the amino acid sequence.\nExonerate does not insert a corresponding gap for those in the nucleotide sequence.\nIt is recommended to remove all X from the amino acid sequences.\n";
}

if (@doublestrange) {
	print "\nThe following sequences turned up a very strange alignment. You should _really_ double-check these:\n\n";
	while (my $hdr = shift @doublestrange) {
		print $hdr, "\n";
	}
	print "\nThe most common reason is a subsequence missing from either amino acid or nucleotide sequence so that Exonerate finds more than one alignment.\nIt is _recommended_ to fix these sequences.\n";
}

print "All done.\n";

exit;

# generates a temporary fasta file
sub fastaify {
	my $fh = File::Temp->new();
	printf $fh ">%s\n%s\n", $_[0], $_[1];
	close $fh;
	return $fh;
}

# loads a Fasta file into a hashref
# arguments: scalar string path to file
# returns: hashref (header => sequence)
sub slurp_fasta {
	my $infile = shift;
	my $sequences = {};
	my $infh = Seqload::Fasta->open($infile);
	while (my ($h, $s) = $infh->next_seq()) {
		# make sure the header is unique
		print "Warning: Non-unique header '$h': " and return 0 if $sequences->{$h};
		# ok you got it
		$sequences->{$h} = $s;
	}
	undef $infh;
	return $sequences;
}

# return an array from a whitespace-separated string
sub QW {
	my $s = shift;
	$s =~ s/\n/ /g;
	return split /\s+/, $s;
}

# Object-oriented access to Fasta files
package Seqload::Fasta;
use strict;
use warnings;
use Carp;

# Constructor. Returns a sequence database object.
sub open {
  my ($class, $filename) = @_;
  open (my $fh, '<', $filename)
    or confess "Fatal: Could not open $filename\: $!\n";
  my $self = {
    'filename' => $filename,
    'fh'       => $fh
  };
  bless($self, $class);
  return $self;
}

# Returns the next sequence as an array (hdr, seq). 
# Useful for looping through a seq database.
sub next_seq {
  my $self = shift;
  my $fh = $self->{'fh'};
	# this is the trick that makes this work
  local $/ = "\n>"; # change the line separator
  return unless defined(my $item = readline($fh));  # read the line(s)
  chomp $item;
  
  if ($. == 1 and $item !~ /^>/) {  # first line is not a header
    croak "Fatal: " . $self->{'filename'} . "is not a FASTA file: Missing descriptor line\n";
  }

	# remove the '>'
  $item =~ s/^>//;

	# split to a maximum of two items (header, sequence)
  my ($hdr, $seq) = split(/\n/, $item, 2);
	$hdr =~ s/\s+$//;	# remove all trailing whitespace
  $seq =~ s/>//g if defined $seq;
  $seq =~ s/\s+//g if defined $seq; # remove all whitespace, including newlines

  return($hdr, $seq);
}

# Closes the file and undefs the database object.
sub close {
  my $self = shift;
  my $fh = $self->{'fh'};
  my $filename = $self->{'filename'};
  close($fh) or carp("Warning: Could not close $filename\: $!\n");
  undef($self);
}

# Destructor. This is called when you undef() an object
sub DESTROY {
  my $self = shift;
  $self->close;
}

