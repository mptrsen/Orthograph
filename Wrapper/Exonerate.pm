package Wrapper::Exonerate;

use strict;
use warnings;
use File::Basename; # basename of files
use IO::File; # object-oriented access to files
use Carp; # extended dying functions
use Data::Dumper;

my $verbose = 0;
my $debug = 0;
my $outdir = File::Spec->catdir('.');
my $searchprog = 'exonerate';
my @searchcmd;
my $evalue_threshold = 10;
my $score_threshold = 10;

sub new {
	my ($class, $query, $target) = @_;

	my $self = {
		'query'      => $query,
		'target'     => $target,
		'resultfile' => '',
		'hitcount'   => 0,
	}
	bless $self, $class;
	return $self;
}

=head1 Class methods

=head2 verbose

Sets verbose output on (TRUE) or off (FALSE). Default is FALSE.

=cut

sub verbose {#{{{
  my $class = shift;
  if (ref $class) { confess("Class method called as object method") }
  unless (scalar @_ == 1) { confess("Usage: Wrapper::Exonerate->verbose(1|0)") }
  $verbose = shift;
}#}}}

=head2 debug

Sets debug output on (TRUE) or off (FALSE). Default is FALSE.

=cut

sub debug {#{{{
  my $class = shift;
  if (ref $class) { confess("Class method called as object method") }
  unless (scalar @_ == 1) { confess("Usage: Wrapper::Exonerate->debug(1|0)") }
  $debug = shift;
}#}}}

=head2 outdir

Sets output directory for the hmmsearch output files. Expects a reference to
scalar F<pathname>. Defaults to F<.>.

=cut

sub outdir {#{{{
  my $class = shift;
  if (ref $class) { confess("Class method called as object method") }
  unless (scalar @_ == 1) { confess("Usage: Wrapper::Exonerate->outdir(OUTDIR)") }
  $outdir = shift;
}#}}}

=head2 searchprog

Sets the HMMsearch program. Expects a string. Defaults to 'F<hmmsearch>'.

=cut

sub searchprog {#{{{
  my $class = shift;
  if (ref $class) { confess("Class method called as object method") }
  unless (scalar @_ == 1) { confess("Usage: Wrapper::Exonerate->searchprog(COMMAND)") }
  $searchprog = shift;
}#}}}

=head2 score_threshold

Sets or returns the score threshold to use for the blastp search. Defaults to
0 (disabled). Note that e-value and score thresholds are mutually exclusive; if
you set one, this automatically unsets the other.

=cut

sub score_threshold {#{{{
	my $class = shift;
	if (ref $class) { confess("Class method used as object method\n") }
	if (scalar(@_) == 0) { return $score_threshold }
	if (scalar(@_) >  1) { confess("Usage: Wrapper::Exonerate->score_threshold(N)\n") }
	$score_threshold = shift(@_);
	unless ($score_threshold =~ /^[0-9]+$/) { confess("Invalid argument (must be integer): $score_threshold\n") }
	$evalue_threshold = 0;
}#}}}

sub resultfile {
	my $self       = shift;
	my $resultfile = shift;
	$self->{'resultfile'} = $resultfile;
	return 1;
}

=head2 search

Searches query against target sequence(s) and stores the result in the object (actually, only the output file location is stored, but this is parsed once you access the result).

=cut

sub search {
	my $self = shift;
	# some exonerate options
	my $exonerate_model = $exhaustive ? 'protein2genome:bestfit' : 'protein2genome';
	my $exonerate_exhaustive = $exhaustive ? '--exhaustive yes' : '';

	# roll your own output for exonerate
	my $exonerate_ryo = "Score: %s\n%V\n>%qi_%ti_[%tcb:%tce]_cdna\n%tcs//\n>%qi[%qab:%qae]_query\n%qas//\n>%ti[%tab:%tae]_target\n%tas//\n";

	# the complete command line
	my $exonerate_cmd = qq( $searchprog --score $self->score_threshold --ryo '$exonerate_ryo' --model $exonerate_model --verbose 0 --showalignment no --showvulgar no $exonerate_exhaustive $self->query_sequence $self->target_sequence 2> /dev/null );
}

=head2 query_sequence

Sets or returns the query sequence.

=cut

sub query_sequence {
	my $self = shift;
	if    (scalar @_ == 0) { return $self->{'query'}->{'sequence'} }
	elsif (scalar @_ > 1 ) { confess 'Usage: Wrapper::Exonerate->query_sequence($sequence)', "\n" }
	$self->{'query'}->{'sequence'} = shift @_;
}

=head2 target_sequence

Sets or returns the target sequence.

=cut

sub target_sequence {
	my $self = shift;
	if    (scalar @_ == 0) { return $self->{'target'}->{'sequence'} }
	elsif (scalar @_ > 1 ) { confess 'Usage: Wrapper::Exonerate->target_sequence($sequence)', "\n" }
	$self->{'target'}->{'sequence'} = shift @_;
}
