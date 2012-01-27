#--------------------------------------------------
# This file is part of Forage.
# Copyright 2012 Malte Petersen <mptrsen@uni-bonn.de>
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
package Forage::Hmmsearch;
use strict;
use warnings;
use File::Basename; # basename of files
use IO::File; # object-oriented access to files
use Carp; # extended dying functions
use Data::Dumper;
my $verbose = 0;
my $debug = 0;
my $hmmoutdir = File::Spec->catdir('.');
my $hmmsearchprog = 'hmmsearch';
my @hmmsearchcmd;
my $hmmfullout = 0;
my $evalue_threshold = 10;
my $score_threshold = 0;
1;

sub new {
  my ($class, $hmmfile) = @_;

  my $self = {
    'hmmfile'       => $hmmfile,
    'hmmresultfile' => '',
    'hmmhitcount'   => 0,
  };

  bless ($self, $class);
  return $self;
}

#--------------------------------------------------
# # Class methods
#-------------------------------------------------- 

=head2 verbose

Sets verbose output on (TRUE) or off (FALSE). Default is FALSE.

=cut

sub verbose {#{{{
  my $class = shift;
  if (ref $class) { confess "Class method called as object method" }
  unless (scalar @_ == 1) { confess "Usage: Forage::Unthreaded->verbose(1|0)" }
  $verbose = shift;
}#}}}

=head2 debug

Sets debug output on (TRUE) or off (FALSE). Default is FALSE.

=cut

sub debug {#{{{
  my $class = shift;
  if (ref $class) { confess "Class method called as object method" }
  unless (scalar @_ == 1) { confess "Usage: Forage::Unthreaded->debug(1|0)" }
  $debug = shift;
}#}}}

=head2 outdir

Sets output directory for the hmmsearch output files. Expects a reference to scalar F<pathname>. Defaults to F<.>.

=cut

sub hmmoutdir {#{{{
  my $class = shift;
  if (ref $class) { confess "Class method called as object method" }
  unless (scalar @_ == 1) { confess "Usage: Forage::Unthreaded->hmmoutdir(OUTDIR)" }
  $hmmoutdir = shift;
}#}}}

=head2 hmmsearchprog

Sets the HMMsearch program. Expects a string. Defaults to 'F<hmmsearch>'.

=cut

sub hmmsearchprog {#{{{
  my $class = shift;
  if (ref $class) { confess "Class method called as object method" }
  unless (scalar @_ == 1) { confess "Usage: Forage::Unthreaded->hmmsearchcmd(COMMAND)" }
  $hmmsearchprog = shift;
}#}}}

=head2 hmmfullout

Sets the option for full hmmsearch output (parsing of which is not yet implemented). Can be TRUE or FALSE. Defaults to FALSE.

=cut

sub hmmfullout {#{{{
  my $class = shift;
  if (ref $class) { confess "Class method called as object method" }
  unless (scalar @_ == 1) { confess "Usage: Forage::Unthreaded->hmmfullout(1|0)" }
  $hmmfullout = shift;
  if ($hmmfullout) { 
    @hmmsearchcmd = qq($hmmsearchprog -E $evalue_threshold -T $score_threshold);
  }
  else { @hmmsearchcmd = qq($hmmsearchprog -E $evalue_threshold -T $score_threshold) }
}#}}}

=head3 evalue_threshold

Sets or returns the e-value threshold to use for the blastp search. Defaults to 10.

=cut

sub evalue_threshold {#{{{
	my $class = shift;
	if (ref $class) { confess "Class method used as object method\n" }
	unless (@_ == 1) { confess "Usage: Forage::Blast->evalue_threshold(N)\n" }
	$evalue_threshold = shift;
}#}}}

=head3 score_threshold

Sets or returns the e-value threshold to use for the blastp search. Defaults to 0.

=cut

sub score_threshold {#{{{
	my $class = shift;
	if (ref $class) { confess "Class method used as object method\n" }
	unless (@_ == 1) { confess "Usage: Forage::Blast->score_threshold(N)\n" }
	$score_threshold = shift;
}#}}}

#--------------------------------------------------
# # Object methods
#-------------------------------------------------- 


=head2 hmmsearch

HMM-searches a sequence file using F<hmmfile>, leaving F<hmmoutfile> for later processing of the results. 

Expects: Reference to sequence object, scalar string filename to F<hmmfile>. 

=cut

sub hmmsearch {#{{{
  my $self = shift;
  unless (scalar @_ == 1) { confess "Usage: OBJECT->hmmsearch(FILE)" }
  my $protfile  = shift;
  my $hmmfile   = $self->hmmfile;
  # full output if desired, table only otherwise; reflects in outfile extension
  my $hmmoutfile = $hmmfullout ? 
    File::Spec->catfile($hmmoutdir, basename($hmmfile).'.out') : 
    File::Spec->catfile($hmmoutdir, basename($hmmfile).'.tbl');
  # return right away if this search has been conducted before
  if (-e $hmmoutfile) {
    $self->hmmresultfile($hmmoutfile);
    return $self;
  }
  else {
    my @hmmsearchline = $hmmfullout ?
      qq($hmmsearchprog -o $hmmoutfile $hmmfile $protfile) :
      qq($hmmsearchprog --tblout $hmmoutfile $hmmfile $protfile);
    print "\n@hmmsearchline\n\n"
      if $debug;
    # do the search
    my $hmmresult = [ `@hmmsearchline` ];
    croak "Fatal: hmmsearch failed on $protfile with HMM $hmmfile: $!\n" 
      unless (scalar @$hmmresult);
    # only save those results that actually match something
    unless (grep( /No hits detected/, @$hmmresult )) {
      $self->{'hmmresultfile'} = $hmmoutfile;
      return $self;
    }
    # empty result
    else {
      $self->{'hmmresultfile'} = $hmmoutfile;
      return $self;
    }
  }
}#}}}


=head2 hmmhitcount

Returns the number of hits from a hmmsearch result file.

=cut

sub hmmhitcount {#{{{
  my $self = shift;
  if ($self->{'hmmhitcount'}) { 
    return $self->{'hmmhitcount'};
  }
  unless ($hmmfullout) {
    $self->{'hmmhitcount'} = scalar(@{$self->hmmresult}); 
    return $self->{'hmmhitcount'};
  }
  # dunno what do with hmmfullout yet... TODO implement!
}#}}}

=head2 hmmresult

Returns the hmmsearch result as it is in the result file, sans the first 3 lines of the table (comments)

=cut

sub hmmresult {#{{{
  my $self = shift;
  if ($self->{'hmmresult'}) {
    return $self->{'hmmresult'};
  }
  my $fh = IO::File->new($self->hmmresultfile());
  $self->{'hmmresult'} = [ <$fh> ];
  $fh->close;
  splice(@{$self->{'hmmresult'}}, 0, 3);
  return $self->{'hmmresult'};
}#}}}

=head2 hmmhits_arrayref

Returns an array reference to a list of lists, e.g., like so:
  $hmmhits->[$i][0..3]  # of line $i, fields 1, 3, 5, 6 of the hmmsearch table output

=cut
sub hmmhits_arrayref {#{{{
  my $self = shift;
  if ($self->{'hmmhits'}) {
    return $self->{'hmmhits'};
  }
  $self->{'hmmhits'} = [ ];
  foreach (@{$self->hmmresult}) {
    # maximum of 19 columns, the last one may contain whitespace
    my @line = split(/\s+/);  
    push(@{$self->{'hmmhits'}}, {
      'target' => $line[0],
      'query'  => $line[2],
      'evalue' => $line[4],
      'score'  => $line[5],
    });
  }
  # this is an array reference
  return $self->{'hmmhits'};
}#}}}

=head2 hmmname

Sets or returns the name of the HMM that was used (may differ from the HMM filename).

=cut
sub hmmname {#{{{
  my $self = shift;
  if ($self->{'hmmname'}) {
    return $self->{'hmmname'};
  }
	my @line = split(/\s+/, ${$self->hmmresult}[0]);
  $self->{'hmmname'} = $line[2];
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

=head2 hmmresultfile

Sets or returns the HMMsearch result filename.

=cut

sub hmmresultfile {#{{{
  my $self = shift;
  if (scalar @_ == 1) {
    $self->{'hmmresultfile'} = shift; 
    return 1;
  }
  return $self->{'hmmresultfile'};
}#}}}

