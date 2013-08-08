#--------------------------------------------------
# This file is part of Orthograph.
# Copyright 2012 Malte Petersen <mptrsen@uni-bonn.de>
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
use File::Temp;
use IO::File; # object-oriented access to files
use Carp; # extended dying functions
use Data::Dumper;

my $verbose    = 0;
my $debug      = 1;
my $exhaustive = 0;
my $outdir     = File::Spec->catdir('.');
my $searchprog = 'exonerate';
my @searchcmd;
my $evalue_threshold = 10;
my $score_threshold = 10;

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
	if    (scalar @_ == 0) { return $self->{'query'}->{'sequence'} }
	elsif (scalar @_ > 1 ) { confess 'Usage: Wrapper::Exonerate->query_sequence($sequence)', "\n" }
	$self->{'resultfile'} = shift;
	return 1;
}

=head2 search

Searches query against target sequence(s) and stores the result in the object (actually, only the output file location is stored, but this is parsed once you access the result).

=cut

sub search {
	my $self = shift;
	# some exonerate options
	my $exonerate_model = $exhaustive ? 'protein2genome:bestfit' : 'protein2genome';
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
	#my $exonerate_ryo = "Score: %s\n%V\n>%qi_%ti_[%tcb:%tce]_cdna\n%tcs//\n>%qi[%qab:%qae]_query\n%qas//\n>%ti[%tab:%tae]_target\n%tas//\n";
	# just the target coding sequence (tcs)
	my $exonerate_ryo = "%tcs";

	# the complete command line
	my $exonerate_cmd = qq($searchprog --bestn 1 --score $score_threshold --ryo '$exonerate_ryo' --model $exonerate_model --verbose 0 --showalignment no --showvulgar no $exhaustive $queryfile $targetfile > $outfile);
	print "$exonerate_cmd\n" if $debug;

	# run the beast now
	system($exonerate_cmd) and confess "Error running exonerate: $!\n";
	$self->{'resultfile'} = $outfile;
	return 1;
}

sub cdna_start {
	my $self = shift;
	return $self->{'cdna_start'};
}

sub cdna_end {
	my $self = shift;
	return $self->{'cdna_end'};
}

sub result {
	my $self = shift;
	if ($self->{'result'}) { return $self->{'result'} }
	else {
		my $fh = IO::File->new($self->{'resultfile'});
		$self->{'result'} = [ <$fh> ];
		$fh->close;
		chomp @{$self->{'result'}};
		return $self->{'result'};
	}
}

sub cdna_sequence {
	my $self = shift;
	$self->{'cdna_sequence'} = join '', @{ $self->result };
	$self->{'cdna_sequence'} =~ s/\s//g;
	return $self->{'cdna_sequence'};
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
	$fh->unlink_on_destroy(0) if $debug;
	printf { $fh } ">%s\n%s\n", $header, $sequence;
	close $fh;
	print "Wrote '$header' to Fasta file '$fh'\n" if $debug;
	return $fh;
}

