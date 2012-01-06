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
package Forage::Unthreaded;
use strict;
use warnings;
use File::Basename;	# basename of files
use IO::File;	# object-oriented access to files
use Carp;	# extended dying functions
use Data::Dumper;
my $verbose = 0;
my $hmmoutdir = File::Spec->catdir('.');
my $hmmsearchprog = 'hmmsearch';
my $hmmsearchcmd;
my $hmmfullout = 0;
1;

sub new {
	my ($class, $hmmfile, $protfile) = @_;

	my $self = {
		'hmmfile'				=> $hmmfile,
		'hmmresultfile'	=> '',
		'hmmhits'				=> 0,
	};

	bless ($self, $class);
	return $self;
}

#--------------------------------------------------
# # Class methods
#-------------------------------------------------- 

# Sub: verbose
# sets $verbose
sub verbose {
	my $class = shift;
	if (ref $class) { confess "Class method called as object method" }
	unless (scalar @_ == 1) { confess "Usage: Forage::Unthreaded->verbose(1|0)" }
	$verbose = shift;
}

# Sub: outdir
# sets output dir $outdir
# Expects: reference to scalar path
sub hmmoutdir {
	my $class = shift;
	if (ref $class) { confess "Class method called as object method" }
	unless (scalar @_ == 1) { confess "Usage: Forage::Unthreaded->hmmoutdir(OUTDIR)" }
	$hmmoutdir = shift;
}

# Sub: hmmsearchcmd
# sets hmmsearch command
# Expects: reference to array
sub hmmsearchcmd {
	my $class = shift;
	if (ref $class) { confess "Class method called as object method" }
	unless (scalar @_ == 1) { confess "Usage: Forage::Unthreaded->hmmsearchcmd(COMMAND)" }
	$hmmsearchcmd = shift;
}

sub hmmfullout {
	my $class = shift;
	if (ref $class) { confess "Class method called as object method" }
	unless (scalar @_ == 1) { confess "Usage: Forage::Unthreaded->hmmfullout(1|0)" }
	$hmmfullout = shift;
}

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
	my $protfile	= shift;
	my $hmmfile		= $self->hmmfile;
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
		confess "Fatal: hmmsearch failed on $protfile with HMM $hmmfile: $!\n" 
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
sub hmmhitcount {
	my $self = shift;
	if ($self->{'hmmhits'}) { 
		return $self->{'hmmhits'};
	}
	if ($self->{'hmmresult'}) {
		unless ($hmmfullout) {
			$self->{'hmmhits'} = scalar(@{$self->{'hmmresult'}}) - 3;	# -3 because the first 3 lines of the table are comments
			return $self->{'hmmhits'};
		}
		croak "THIS NEEDS TO BE IMPLEMENTED, CAN'T WORK WITH FULL HMMSEARCH OUTPUT YET\n";
	}
	unless ($hmmfullout) {
		$self->{'hmmhits'} = scalar(@{$self->hmmresult}) - 3;	# -3 because the first 3 lines of the table are comments
		return $self->{'hmmhits'};
	}
}

# sub: hmmresult
# returns: the hmmsearch result as it is in the result file, sans the first 3 lines of the table (comments)
sub hmmresult {
	my $self = shift;
	if ($self->{'hmmresult'}) {
		return $self->{'hmmresult'};
	}
	my $fh = IO::File->new();
	$fh->open($self->{'hmmresultfile'})
		or croak 'Fatal: Could not open hmmsearch result file ', $self->{'hmmresultfile'}, ': ', $!, "\n";
	$self->{'hmmresult'} = [ <$fh> ];
	$fh->close;
	return $self->{'hmmresult'};
}

# hmm file used for searching
sub hmmfile {
	my $self = shift;
	return $self->{'hmmfile'};
}

# protein file (normally: EST input file)
sub protfile {
	my $self = shift;
	return $self->{'protfile'};
}

# hmmsearch result file (table or fullout)
sub hmmresultfile {
	my $self = shift;
	if (scalar @_ == 1) {
		$self->{'hmmresultfile'} = shift; 
		return 1;
	}
	return ${$self->{'hmmresultfile'}};
}

