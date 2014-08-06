#!/usr/bin/perl
use strict;
use warnings;
use autodie;
use File::Temp;
use File::Basename;
use File::Spec;
use Getopt::Long;

my $skipfirst = 0;
my $skipall   = 0;
my $outdir    = '/tmp';
my $exonerate = 'exonerate';
my $score_threshold = 30;

GetOptions(
	'skipfirst'         => \$skipfirst,
	'skipall'           => \$skipall,
	'outdir=s'          => \$outdir,
	'exonerate=s'       => \$exonerate,
	'score-threshold=d' => \$score_threshold,
);

my $usage = "Usage: $0 OGSFILE TRANSCRIPTOMEFILE\n";

scalar @ARGV == 2 or die $usage;

print "Call: $0 @ARGV\n";

my $ogs = slurp_fasta($ARGV[0]);

my $transcripts = slurp_fasta($ARGV[1]);

unless (scalar keys %$ogs == scalar keys %$transcripts) {
	print "Unequal number of sequences!\n";
}

my %seen    = ();
my $nupd    = 0;
my $nunchgd = 0;
my $n       = 0;
my $nseqs   = scalar keys %$ogs;
my $c       = 0;

my $new_ogs_file = File::Spec->catfile($outdir, 'corresp-' . basename($ARGV[0]));
my $new_transcripts_file = File::Spec->catfile($outdir, 'corresp-' . basename($ARGV[1]));
open my $new_ogs, '>', $new_ogs_file;
open my $new_transcripts, '>', $new_transcripts_file;

foreach my $hdr (sort {$a cmp $b} keys %$ogs) {

	++$c;

	# if the headers aren't identical, the sequences can't be correlated. 
	# in this case, they may be skipped at the user's discretion (--skipfirst or
	# --skipall options). otherwise, the program will exit.
	unless (exists $transcripts->{$hdr}) {
		if ($skipall) {
			print "Not found in transcriptome: '$hdr', skipping\n";
			next;
		}
		elsif ($skipfirst) {
			print "Not found in transcriptome: '$hdr', skipping\n";
			$skipfirst = 0;
			next;
		}
		else {
			print "Not found in transcriptome: '$hdr', exiting\n" and exit(1);
		}
	}

	# print to files 
	my $aafn    = fastaify($hdr, $ogs->{$hdr});
	my $ntfn    = fastaify($hdr, $transcripts->{$hdr});
	my $outfile = File::Temp->new();

	# some settings for exonerate
	my $exonerate_model = 'protein2genome';
	my $exonerate_ryo   = '>ca\n%tcs>qa\n%qas';

	# set output buffer to flush immediately
	local $| = 0;
	print "Checking $hdr ($c of $nseqs)... ";
	
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

	system("@command") and die $?;

	# get the alignment results
	my $res = slurp_fasta($outfile);
	# no alignment found, something is wrong, skip this pair
	if (!$res->{'ca'}) {
		print "done, no alignment found, skipping\n";
		next;
	}
	# the sequence has been updated
	elsif ($res->{'ca'} ne $transcripts->{$hdr}) {
		printf $new_ogs ">%s\n%s\n", $hdr, $res->{'qa'};
		printf $new_transcripts ">%s\n%s\n", $hdr, $res->{'ca'};
		# count
		++$nupd;
		++$n;
		print "done, updated\n";
	}
	# the sequence is the same
	else {
		printf $new_ogs ">%s\n%s\n", $hdr, $res->{'qa'};
		printf $new_transcripts ">%s\n%s\n", $hdr, $res->{'ca'};
		# count
		++$nunchgd;
		++$n;
		print "done, unchanged\n";
	}

	# set output buffer back to normal
	local $| = 1;

	# delete tempfile
	undef $outfile;
}

printf "Done, updated %d of %d sequences, wrote %d sequences to %s and %s\n",
	$nupd,
	$n,
	$nupd + $nunchgd,
	$new_ogs_file,
	$new_transcripts_file,
;

exit;

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
		#--------------------------------------------------
		# # remove possible taxon shorthand
		# $h =~ s/^[A-Z]{5} //;
		# # we only need the id field
		# my @fields = split /\s+/, $h;
		# $h = $fields[0];
		# # remove possible -R*/-P* suffixes
		# $h =~ s/-[RPT][A-H]$//;
		# # remove pipes
		# $h =~ s/\|/_/g;
		#-------------------------------------------------- 
		# make sure the header is unique
		die "Non-unique header: $h\n" if $sequences->{$h};
		# ok you got it
		$sequences->{$h} = $s;
	}
	undef $infh;
	return $sequences;
}

sub QW {
	my $s = shift;
	$s =~ s/\n/ /g;
	return split /\s+/, $s;
}

package Seqload::Fasta;
use strict;
use warnings;
use Carp;
require Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw( fasta2csv check_if_fasta );

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

# Convert a fasta file to a csv file the easy way
# Usage: Seqload::Fasta::fasta2csv($fastafile, $csvfile);
sub fasta2csv {
  my $fastafile = shift;
  my $csvfile = shift;

  my $fastafh = Seqload::Fasta->open($fastafile);
  CORE::open(my $outfh, '>', $csvfile)
    or confess "Fatal: Could not open $csvfile\: $!\n";
  while (my ($hdr, $seq) = $fastafh->next_seq) {
		$hdr =~ s/,/_/g;	# remove commas from header, they mess up a csv file
    print $outfh $hdr . ',' . $seq . "\n"
			or confess "Fatal: Could not write to $csvfile\: $!\n";
  }
  CORE::close $outfh;
  $fastafh->close;

  return 1;
}

# validates a fasta file by looking at the FIRST (header, sequence) pair
# arguments: scalar string path to file
# returns: true on validation, false otherwise
sub check_if_fasta {
	my $infile = shift;
	my $infh = Seqload::Fasta->open($infile);
	my ($h, $s) = $infh->next_seq() or return 0;
	return 1;
}

