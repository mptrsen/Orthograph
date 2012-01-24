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
package Forage::Hmmsearch;
use strict;
use warnings;
use File::Basename; # basename of files
use IO::File; # object-oriented access to files
use Carp; # extended dying functions
my $verbose = 0;
my $hmmoutdir = File::Spec->catdir('.');
my $hmmsearchprog = 'hmmsearch';
my $hmmsearchcmd;
my $hmmfullout = 0;
1;

sub new {
  my ($class, $hmmfile, $protfile) = @_;

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

# Sub: verbose
# sets $verbose
sub verbose {#{{{
  my $class = shift;
  if (ref $class) { confess "Class method called as object method" }
  unless (scalar @_ == 1) { confess "Usage: Forage::Unthreaded->verbose(1|0)" }
  $verbose = shift;
}#}}}

# Sub: outdir
# sets output dir $outdir
# Expects: reference to scalar path
sub hmmoutdir {#{{{
  my $class = shift;
  if (ref $class) { confess "Class method called as object method" }
  unless (scalar @_ == 1) { confess "Usage: Forage::Unthreaded->hmmoutdir(OUTDIR)" }
  $hmmoutdir = shift;
}#}}}

# Sub: hmmsearchcmd
# sets hmmsearch command
# Expects: reference to array
sub hmmsearchcmd {#{{{
  my $class = shift;
  if (ref $class) { confess "Class method called as object method" }
  unless (scalar @_ == 1) { confess "Usage: Forage::Unthreaded->hmmsearchcmd(COMMAND)" }
  $hmmsearchcmd = shift;
}#}}}

# Sub: hmmfullout
# sets option for full hmmsearch output
# Expects: TRUE or FALSE
sub hmmfullout {#{{{
  my $class = shift;
  if (ref $class) { confess "Class method called as object method" }
  unless (scalar @_ == 1) { confess "Usage: Forage::Unthreaded->hmmfullout(1|0)" }
  $hmmfullout = shift;
}#}}}

#--------------------------------------------------
# # Object methods
#-------------------------------------------------- 

# Sub: hmmsearch
# HMM search a sequence using HMM, leaving an outfile for later processing
# Expects: reference to sequence object, scalar string filename to HMM
# Returns: scalar reference to hmmoutfile
sub hmmsearch {#{{{
  my $self = shift;
  unless (scalar @_ == 1) { confess "Usage: OBJECT->hmmsearch(FILE)" }
  my $protfile  = shift;
  my $hmmfile   = $self->hmmfile;
  # full output if desired, table only otherwise; reflects in outfile extension
  my $hmmoutfile = $hmmfullout ? 
    File::Spec->catfile($hmmoutdir, basename($hmmfile).'.out') : 
    File::Spec->catfile($hmmoutdir, basename($hmmfile).'.tbl');
  # e-value and score options if desired
  if (-e $hmmoutfile) {
    $self->hmmresultfile($hmmoutfile);
    return $self;
  }
  else {
    my @hmmsearchline = (@$hmmsearchcmd, $hmmoutfile, $hmmfile, $protfile);
    print join " ", @hmmsearchline, "\n"
      if $verbose;
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

# sub: hmmhitcount
# extracts the number of hits from a hmm search result file
# returns: int number of hmm hits
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

# sub: hmmresult
# returns: the hmmsearch result as it is in the result file, sans the first 3 lines of the table (comments)
sub hmmresult {#{{{
  my $self = shift;
  if ($self->{'hmmresult'}) {
    return $self->{'hmmresult'};
  }
  my $fh = IO::File->new($self->{'hmmresultfile'});
  $self->{'hmmresult'} = [ <$fh> ];
  $fh->close;
  splice(@{$self->{'hmmresult'}}, 0, 3);
  return $self->{'hmmresult'};
}#}}}

# sub: hmmhits_arrayref
# returns: array reference to list of list
# $hmmhits->[$i][0..3] (of line $i, fields 1, 3, 5, 6 of hmmsearch table output)
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
      'eval'   => $line[4],
      'score'  => $line[5],
    });
  }
  # this is an array reference
  return $self->{'hmmhits'};
}#}}}

# sub: hmmname
# returns: name of the HMM that was used (may differ from the HMM file)
sub hmmname {#{{{
  my $self = shift;
  if ($self->{'hmmname'}) {
    return $self->{'hmmname'};
  }
  return ${$self->hmmresult}[1];
}#}}}

# hmm file used for searching
sub hmmfile {#{{{
  my $self = shift;
  return $self->{'hmmfile'};
}#}}}

# protein file (normally: EST input file)
sub protfile {#{{{
  my $self = shift;
  return $self->{'protfile'};
}#}}}

# hmmsearch result file (table or fullout)
sub hmmresultfile {#{{{
  my $self = shift;
  if (scalar @_ == 1) {
    $self->{'hmmresultfile'} = shift; 
    return 1;
  }
  return ${$self->{'hmmresultfile'}};
}#}}}

