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
package Wrapper::Hmmsearch;
use strict;
use warnings;
use autodie;
use File::Basename; # basename of files
use IO::File; # object-oriented access to files
use Carp; # extended dying functions
use Data::Dumper;
my $verbose          = 0;
my $debug            = 0;
my $outdir           = File::Spec->catdir('.');
my $searchprog       = 'hmmsearch';
my @searchcmd;
my $evalue_threshold = 10;
my $score_threshold  = 10;
my $num_threads      = 1;

sub new {
  my ($class, $hmmfile) = @_;

  my $self = {
    'hmmfile'    => $hmmfile,
    'resultfile' => '',
    'hitcount'   => 0,
  };

  bless ($self, $class);
  return $self;
}

=head1 Class methods

=head2 verbose

Sets verbose output on (TRUE) or off (FALSE). Default is FALSE.

=cut

sub verbose {#{{{
  my $class = shift;
  if (ref $class) { confess("Class method called as object method") }
  unless (scalar @_ == 1) { confess("Usage: Wrapper::Hmmsearch->verbose(1|0)") }
  $verbose = shift;
}#}}}

=head2 debug

Sets debug output on (TRUE) or off (FALSE). Default is FALSE.

=cut

sub debug {#{{{
  my $class = shift;
  if (ref $class) { confess("Class method called as object method") }
  unless (scalar @_ == 1) { confess("Usage: Wrapper::Hmmsearch->debug(1|0)") }
  $debug = shift;
}#}}}

=head2 outdir

Sets output directory for the hmmsearch output files. Expects a reference to
scalar F<pathname>. Defaults to F<.>.

=cut

sub outdir {#{{{
  my $class = shift;
  if (ref $class) { confess("Class method called as object method") }
  unless (scalar @_ == 1) { confess("Usage: Wrapper::Hmmsearch->outdir(OUTDIR)") }
  $outdir = shift;
}#}}}

=head2 searchprog

Sets the HMMsearch program. Expects a string. Defaults to 'F<hmmsearch>'.

=cut

sub searchprog {#{{{
  my $class = shift;
  if (ref $class) { confess("Class method called as object method") }
  unless (scalar @_ == 1) { confess("Usage: Wrapper::Hmmsearch->searchprog(COMMAND)") }
  $searchprog = shift;
}#}}}

=head2 evalue_threshold

Sets or returns the e-value threshold to use for the HMM search. Defaults to
10. Note that e-value and score thresholds are mutually exclusive; if you set
one, this automatically unsets the other.

=cut

sub evalue_threshold {#{{{
	my $class = shift;
	if (ref $class) { confess("Class method used as object method\n") }
	if (scalar(@_) == 0) { return $score_threshold };
	if (scalar(@_) >  1) { confess("Usage: Wrapper::Hmmsearch->evalue_threshold(N)\n") }
	$evalue_threshold = shift(@_);
	# check whether this is a valid number
	unless ($evalue_threshold =~ /^[0-9]+(\.[0-9]+)?([eE][-+]?[0-9]+)?$/) { confess("Invalid argument (must be integer, float or exponential): $evalue_threshold\n") }
	$score_threshold = 0;
}#}}}

=head2 score_threshold

Sets or returns the score threshold to use for the HMM search. Defaults to
0 (disabled). Note that e-value and score thresholds are mutually exclusive; if
you set one, this automatically unsets the other.

=cut

sub score_threshold {#{{{
	my $class = shift;
	if (ref $class) { confess("Class method used as object method\n") }
	if (scalar(@_) == 0) { return $score_threshold }
	if (scalar(@_) >  1) { confess("Usage: Wrapper::Hmmsearch->score_threshold(N)\n") }
	$score_threshold = shift(@_);
	unless ($score_threshold =~ /^[0-9]+$/) { confess("Invalid argument (must be integer): $score_threshold\n") }
	$evalue_threshold = 0;
}#}}}

=head2 num_threads

Sets or returns the number of CPU threads to use for multithreaded searching. Defaults to 1.

=cut

sub num_threads {
	my $class = shift;
	if (ref $class) { confess("Class method used as object method\n") }
	if (scalar(@_) == 0) { return $num_threads }
	if (scalar(@_) >  1) { confess("Usage: Wrapper::Hmmsearch->num_threads(N)\n") }
	$num_threads = shift(@_);
	unless ($num_threads =~ /^[0-9]+$/) { confess("Invalid argument (must be integer): $num_threads\n") }
}

=head1 Object methods

=head2 search(FILE)

Searches the FILE for matches using the hidden Markov model.

Arguments: scalar string filename 

=cut

sub search {#{{{
  my $self = shift;
  unless (scalar @_ == 1) { confess("Usage: OBJECT->search(FILE)") }
  my $protfile  = shift;
  my $hmmfile   = $self->hmmfile;
  # full output if desired, table only otherwise; reflects in outfile extension
  my $outfile = File::Spec->catfile($outdir, basename($hmmfile).'.domtbl');
  # return right away if this search has been conducted before
  if (-e $outfile) {
		print STDERR "HMMsearch output file exists in '$outfile'\n" if $debug;
    $self->resultfile($outfile);
    return $self;
  }
  else {
		print STDERR "HMMsearch output file does not exist in '$outfile', conducting new search\n" if $debug;
		my $threshold_option = $evalue_threshold ? qq(-E $evalue_threshold) : qq(-T $score_threshold);
    my $hmmsearchline = qq($searchprog --domtblout $outfile $threshold_option --cpu $num_threads $hmmfile $protfile);
    print STDERR "\n$hmmsearchline\n\n" if $debug;
    # do the search
    my $result = [ `$hmmsearchline` ];
    confess("Fatal: hmmsearch failed on '$protfile' with HMM '$hmmfile': $!")
      unless (scalar @$result);
		$self->{'resultfile'} = $outfile;
		return $self;
  }
}#}}}


=head2 hitcount

Returns the number of hits from a hmmsearch result file.

=cut

sub hitcount {#{{{
  my $self = shift;
  if ($self->{'hitcount'}) { 
    return $self->{'hitcount'};
  }
	$self->{'hitcount'} = scalar(@{$self->result}); 
	return $self->{'hitcount'};
}#}}}

=head2 result

Returns the hmmsearch result as it is in the result file, sans the first 3 lines of the table (comments)

=cut

sub result {#{{{
  my $self = shift;
  if ($self->{'result'}) {
    return $self->{'result'};
  }
  my $fh = IO::File->new($self->resultfile())
		or croak("Fatal: Could not open resultfile");
  $self->{'result'} = [ <$fh> ];
  $fh->close;
	# remove comment lines
	@{$self->{'result'}} = grep { /^[^#]/ } @{$self->{'result'}};
  return $self->{'result'};
}#}}}

=head2 hmmhits_arrayref

Returns an array reference to a list of lists, e.g., like so:

  $hmmhits->[$i][0..3]  # of line $i, fields 1, 3, 5, 6 of the hmmsearch table output

=cut
sub hits_arrayref {#{{{
  my $self = shift;
  if ($self->{'hits'}) {
    return $self->{'hits'};
  }
  $self->{'hits'} = [ ];
  foreach (@{$self->result}) {
    # maximum of 19 columns, the last one may contain whitespace
    my @line = split(/\s+/);  
    push(@{$self->{'hits'}}, {
      'target'     => $line[0],		# target ID
      'query'      => $line[3],		# query ID
      'evalue'     => $line[12],	# e-value of the best domain
      'score'      => $line[13],	# score of the best domain
      'hmm_start'  => $line[15],	# beginning of domain
      'hmm_end'    => $line[16],	# end of domain
			'ali_start'  => $line[17],	# beginning of domain
			'ali_end'    => $line[18],	# end of domain
      'env_start'  => $line[19],	# beginning of domain
      'env_end'    => $line[20],	# end of domain
    });
  }
  # this is an array reference
  return $self->{'hits'};
}#}}}

=head2 hmmname

Sets or returns the name of the HMM that was used (may differ from the HMM filename).

=cut

sub hmmname {#{{{
  my $self = shift;
  if ($self->{'hmmname'}) {
    return $self->{'hmmname'};
  }
	open my $fh, '<', $self->hmmfile();
	# the second line contains the HMM name
	<$fh>;
	my $line = <$fh>;
	close $fh;
	chomp $line;
	my @fields = split /\s+/, $line;
  $self->{'hmmname'} = $fields[1];
  return $self->{'hmmname'};
}#}}}

=head2 hmmfile

Sets or returns the HMM filename.

=cut

sub hmmfile {#{{{
  my $self = shift;
  if (scalar @_ == 1) {
    $self->{'hmmfile'} = shift;
    return 1;
  }
  return $self->{'hmmfile'};
}#}}}

# protein file (normally: EST input file)
sub protfile {#{{{
  my $self = shift;
  return $self->{'protfile'};
}#}}}

=head2 resultfile

Sets or returns the HMMsearch result filename.

=cut

sub resultfile {#{{{
  my $self = shift;
  if (scalar @_ == 1) {
    $self->{'resultfile'} = shift; 
    return 1;
  }
  return $self->{'resultfile'};
}#}}}

1;

=head1 AUTHOR

Written by Malte Petersen <mptrsen@uni-bonn.de>.

=head1 COPYRIGHT

Copyright (c) 2012 by Malte Petersen. All rights reserved.

This program is free software; you may redistribute and/or modify it under the
same terms as Orthograph itself.

=cut
