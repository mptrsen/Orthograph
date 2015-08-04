#--------------------------------------------------
# This file is part of Orthograph.
# Copyright 2014 Malte Petersen <mptrsen@uni-bonn.de>
# 
# Orthograph is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
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
package Wrapper::Exonerate;

use strict;
use warnings;
use autodie;
use File::Basename; # basename of files
use File::Temp;	# temporary files
use IO::File; # object-oriented access to files
use Carp; # extended dying functions
use Data::Dumper;

use Seqload::Fasta;	# object-oriented access to fasta files

my $verbose          = 0;
my $debug            = 0;
my $exhaustive       = 0;
my $genetic_code     = 1;
my $alignment_model  = 'protein2genome';
my $outdir           = File::Spec->catdir('.');
my $searchprog       = 'exonerate';
my $translateprog    = 'fastatranslate';
my $evalue_threshold = 10;
my $score_threshold  = 10;
my @searchcmd;

sub new {
	my ($class, $query, $target) = @_;

	my $self = {
		'query'      => { 'header' => 'query',  'sequence' => $query  },
		'target'     => { 'header' => 'target', 'sequence' => $target },
		'resultfile' => undef,
		'hitcount'   => 0,
	};
	bless $self, $class;
	return $self;
}

=head1 Class methods

=head2 verbose

Sets verbose output on (TRUE) or off (FALSE). Default is FALSE.

=cut

sub verbose {
  my $class = shift;
  if (ref $class) { confess("Class method called as object method") }
  unless (scalar @_ == 1) { confess("Usage: Wrapper::Exonerate->verbose(1|0)") }
  $verbose = shift;
}

=head2 debug

Sets debug output on (TRUE) or off (FALSE). Default is FALSE.

=cut

sub debug {
  my $class = shift;
  if (ref $class) { confess("Class method called as object method") }
  unless (scalar @_ == 1) { confess("Usage: Wrapper::Exonerate->debug(1|0)") }
  $debug = shift;
}

=head2 exhaustive

Sets exhaustive search. Default is FALSE.

=cut

sub exhaustive {
  my $class = shift;
  if (ref $class) { confess("Class method called as object method") }
  unless (scalar @_ == 1) { confess("Usage: Wrapper::Exonerate->exhaustive(1|0)") }
  $exhaustive = shift;
}

=head2 outdir

Sets output directory for the Exonerate output files. Expects a reference to
scalar F<pathname>. Defaults to F<.>.

=cut

sub outdir {
  my $class = shift;
  if (ref $class) { confess("Class method called as object method") }
  unless (scalar @_ == 1) { confess("Usage: Wrapper::Exonerate->outdir(OUTDIR)") }
  $outdir = shift;
}

=head2 searchprog

Sets the Exonerate program. Expects a string. Defaults to 'F<exonerate>'.

=cut

sub searchprog {
  my $class = shift;
  if (ref $class) { confess("Class method called as object method") }
  unless (scalar @_ == 1) { confess("Usage: Wrapper::Exonerate->searchprog(COMMAND)") }
  $searchprog = shift;
}

=head2 translateprog

Sets the translation program. Expects a string. Defaults to 'F<exonerate>'.

NOTE: the program must be able to be called like this:

PROGRAM -F 1 FASTAFILE

where the -F option specifies the reading frame. It must provide output on
STDOUT and in Fasta format.

=cut

sub translateprog {
  my $class = shift;
  if (ref $class) { confess("Class method called as object method") }
  unless (scalar @_ == 1) { confess("Usage: Wrapper::Exonerate->translateprog(COMMAND)") }
  $translateprog = shift;
}

=head2 genetic_code

Sets or returns the genetic code to use for the exonerate search. Defaults to
1 (standard genetic code). 

=cut

sub genetic_code {
	my $class = shift;
	if (ref $class) { confess("Class method used as object method\n") }
	if (scalar(@_) == 0) { return $genetic_code }
	if (scalar(@_) >  1) { confess("Usage: Wrapper::Exonerate->score_threshold(N)\n") }
	$genetic_code = shift(@_);
	unless ($genetic_code =~ /^[0-9]+$/) { confess("Invalid argument (must be integer): $genetic_code\n") }
}

=head2 alignment_model

Sets or returns the alignment model. Defaults to 'protein2genome' (alignment of a protein sequence to genomic DNA).

=cut

sub alignment_model {
	my $class = shift;
	if (ref $class) { confess("Class method used as object method\n") }
	if (scalar(@_) == 0) { return $alignment_model }
	if (scalar(@_) >  1) { confess("Usage: Wrapper::Exonerate->alignment_model(MODEL)\n") }
	$alignment_model = shift(@_);
	unless ($alignment_model =~ /^protein2(genome|dna)$/) { confess("Invalid argument: $alignment_model\n") }
}

=head2 score_threshold

Sets or returns the score threshold to use for the exonerate search. Defaults to
0 (disabled). 

=cut

sub score_threshold {
	my $class = shift;
	if (ref $class) { confess("Class method used as object method\n") }
	if (scalar(@_) == 0) { return $score_threshold }
	if (scalar(@_) >  1) { confess("Usage: Wrapper::Exonerate->score_threshold(N)\n") }
	$score_threshold = shift(@_);
	unless ($score_threshold =~ /^[0-9]+$/) { confess("Invalid argument (must be integer): $score_threshold\n") }
	$evalue_threshold = 0;
}


=head1 Object methods

=head2 resultfile

Sets or returns the path to the result file.

=cut

sub resultfile {
	my $self       = shift;
	if    (scalar @_ == 0) { return $self->{'resultfile'} }
	elsif (scalar @_ > 1 ) { confess 'Usage: Wrapper::Exonerate->query_sequence($sequence)', "\n" }
	$self->{'resultfile'} = shift;
	return 1;
}

sub delete_resultfile {
	my $self = shift;
	if ($self->resultfile()) {
		unlink $self->resultfile() or warn "Could not delete result file " . $self->resultfile() . "\n";
		return 1;
	}
	return undef;
}

=head2 search

Searches query against target sequence(s) and stores the result in the object (actually, only the output file location is stored, but this is parsed once you access the result).

=cut

sub search {
	my $self = shift;
	# some exonerate options
	my $exonerate_model = $exhaustive ? $alignment_model . ':bestfit' : $alignment_model;
	my $exhaustive = $exhaustive ? '--exhaustive yes' : '';

	my ($queryfile, $targetfile);
	unless ($self->query_file()) {
		$queryfile = &fastaify($self->{'query'}->{'header'}, $self->{'query'}->{'sequence'});
	}
	unless ($self->target_file()) {
		$targetfile = &fastaify($self->{'target'}->{'header'}, $self->{'target'}->{'sequence'});
	}

	$self->query_file($queryfile);
	$self->target_file($targetfile);
  my $outfile = File::Spec->catfile($outdir, $self->{'query'}->{'header'} . '-' . $self->{'target'}->{'header'} . '.exonerateout');

	# roll your own output for exonerate
	# %tcb/%tce : target coding region beginning/end
	# %tcs      : target coding sequence
	# %qab/%qae : query alignment region beginning/end 
	# %qas      : query alignment sequence
	my $exonerate_ryo = '>cdna %tcb %tce\n%tcs>aa %qab %qae\n%qas';

	# the complete command line
	my $exonerate_cmd = qq($searchprog --bestn 1 --score $score_threshold --ryo '$exonerate_ryo' --geneticcode $genetic_code --model $exonerate_model --querytype protein --targettype dna --verbose 0 --showalignment no --showvulgar no $exhaustive --query $queryfile --target $targetfile > $outfile);
	print "$exonerate_cmd\n" if $debug;

	# run the beast now
	system($exonerate_cmd) and confess "Error running exonerate: $!\n";
	$self->resultfile($outfile);
	return 1;
}

sub translated_cdna {
	my $self = shift;
	unless ($self->{'cdna_translated'}) { translate_cdna($self) }
	return $self->{'cdna_translated'};
}

sub translate_cdna {
	my $self = shift;
	my $outfile = File::Spec->catfile($outdir, 'translatethis.fa');
	my $translatefile = fastaify('cdna', $self->{'cdna_sequence'});
	my $translate_cmd = qq($translateprog --geneticcode $genetic_code -F 1 $translatefile);
	if ($debug) {
		print 'Translating CDNA...', "\n";
		print $translate_cmd, "\n";
	}
	my $translated_cdna_fasta = [ `$translate_cmd` ] or croak "Fatal: Couldn't translate CDNA sequence using command '$translate_cmd'\n";
	shift @$translated_cdna_fasta;
	chomp @$translated_cdna_fasta;
	$self->{'cdna_translated'} = join '', @$translated_cdna_fasta;
	return 1;
}

sub cdna_start {
	my $self = shift;
	# add 1 since exonerate uses 0-based coordinates
	return $self->{'cdna_start'} + 1;
}

sub cdna_end {
	my $self = shift;
	return $self->{'cdna_end'};
}

sub aa_start {
	my $self = shift;
	# add 1 since exonerate uses 0-based coordinates
	return $self->{'aa_start'} + 1;
}

sub aa_end {
	my $self = shift;
	return $self->{'aa_end'};
}

sub result {
	my $self = shift;
	if ($self->{'result'}) { return $self->{'result'} }
	else {
		my $fh = IO::File->new($self->resultfile());
		$self->{'result'} = [ <$fh> ];
		$fh->close;
		chomp @{$self->{'result'}};
		return $self->{'result'};
	}
}

sub aa_sequence {
	my $self = shift;
	unless ($self->{'aa_sequence'}) { $self->parse_result() }
	unless ($self->{'aa_sequence'}) { return undef }
	$self->{'aa_sequence'} =~ s/\s//g;
	return $self->{'aa_sequence'};
}

sub cdna_sequence {
	my $self = shift;
	unless ($self->{'cdna_sequence'}) { $self->parse_result() }
	unless ($self->{'cdna_sequence'}) { return undef }
	$self->{'cdna_sequence'} =~ s/\s//g;
	return $self->{'cdna_sequence'};
}

sub parse_result {
	my $self = shift;
	my $header;
	# if there was no result
	if (-z $self->resultfile()) { return undef }

	# otherwise, continue
	my $fh = Seqload::Fasta->open($self->resultfile());
	# watch for the order of sequences, they must correspond to the order in
	# the --ryo option in the exonerate call
	# get the cdna sequence
	($header, $self->{'cdna_sequence'}) = $fh->next_seq();
	# get the cdna coordinates 
	my @fields = split ' ', $header;
	$self->{'cdna_start'} = $fields[-2];
	$self->{'cdna_end'}   = $fields[-1];
	# get the aa sequence
	($header, $self->{'aa_sequence'})   = $fh->next_seq();
	# get the aa coordinates
	@fields = split ' ', $header;
	$self->{'aa_start'} = $fields[-2];
	$self->{'aa_end'}   = $fields[-1];
	undef $fh;
}


=head2 query

Sets or returns the query sequence.

=cut

sub query {
	my $self = shift;
	if    (scalar @_ == 0) { return $self->{'query'} }
	elsif (scalar @_ > 1 ) { confess 'Usage: Wrapper::Exonerate->query_sequence($sequence)', "\n" }
	$self->{'query'} = shift @_;
	return 1;
}

=head2 target

Sets or returns the target sequence.

=cut

sub target {
	my $self = shift;
	if    (scalar @_ == 0) { return $self->{'target'} }
	elsif (scalar @_ > 1 ) { confess 'Usage: Wrapper::Exonerate->target_sequence($sequence)', "\n" }
	$self->{'target'} = shift @_;
	return 1;
}

sub query_file {
	my $self = shift;
	if (scalar @_ == 0) { return $self->{'queryfile'} }
	else { $self->{'queryfile'} = shift @_ }
}

sub target_file {
	my $self = shift;
	if (scalar @_ == 0) { return $self->{'targetfile'} }
	else { $self->{'targetfile'} = shift @_ }
}

sub fastaify {
	my $header = shift;
	my $sequence = shift;
	my $fh = File::Temp->new(UNLINK=>1);
	printf { $fh } ">%s\n%s\n", $header, $sequence;
	close $fh;
	if ($debug > 1) {
		printf "Wrote this sequence to Fasta file '%s':\n>%s\n%s\n", $fh, $header, $sequence;
	}
	return $fh;
}

'this line intentionally left true'
