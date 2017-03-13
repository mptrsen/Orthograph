#--------------------------------------------------
# This file is part of Orthograph.
# Copyright 2013 Malte Petersen <mptrsen@uni-bonn.de>
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

=head1 NAME 

B<Wrapper::Swipe>

=head1 DESCRIPTION

The B<Wrapper::Swipe> module provides an object-oriented interface to the
SWIPE sequence search algorithm. It does not conform to the Guidelines of
the BioPerl package and as such does not return or handle Bioperl's Bio::Seq or
Bio::Search objects. Instead, it is designed to be a lot simpler and more
straightforward to use for applications that only want to use the SWIPE
package without the immense overhead of the entire Bioperl backend.

=head1 SYNOPSIS

  use Wrapper::Swipe;

  # set up the swipr program
  Wrapper::Swipe->searchprog('swipe');

  # create a new Swipe object
  my $searchobj = Wrapper::Swipe->new('/path/to/blast/database');

  # do the swipe search
  $searchobj->search($infile, $outfile);

  # get results 
  $searchobj->hits_arrayref()

=cut

package Wrapper::Swipe;
use strict;
use warnings;
use File::Basename; # basename of files
use IO::File; # object-oriented access to files
use Carp; # extended dying functions
use Data::Dumper;
my $verbose          = 0;
my $debug            = 0;
my $outdir      = File::Spec->catdir('.');
my $searchprog       = 'swipe';
my $makeblastdbprog  = 'makeblastdb';
my $evalue_threshold = 10;
my $score_threshold  = 0;
my $max_hits         = 100;
my $cmd;
my $num_threads      = 1;

sub new {
	my ($class, $db) = @_;

	my $self = {
		'db'       => $db,
		'query'    => undef,
		'hitcount' => undef,
	};

	bless($self, $class);
	return $self;
}

=head1 CLASS METHODS

=head3 verbose()

Sets $verbose. Defaults to 0.

=cut

sub verbose {
	my $class = shift;
	if (ref $class) { confess("Class method used as object method\n") }
	unless (scalar @_ == 1) { confess("Usage: Wrapper::Swipe->verbose(1)\n") }
	$verbose = shift;
}

=head3 debug()

Sets $debug. Defaults to 0.

=cut

sub debug {
	my $class = shift;
	if (ref $class) { confess("Class method used as object method\n") }
	unless (scalar @_ == 1) { confess("Usage: Wrapper::Swipe->debug(1)\n") }
	$debug = shift;
}

=head3 searchprog()

Sets the swipe program. Defaults to B<swipe>.

=cut

sub searchprog {
	my $class = shift;
	if (ref $class) { confess("Class method used as object method\n") }
	unless (scalar @_ == 1) { confess("Usage: Wrapper::Swipe->searchprog(COMMAND)\n") }
	$searchprog = shift;
}

=head3 set_makeblastdb()

Sets the makeblastdb program. Defaults to B<makeblastdb>.

=cut

sub set_makeblastdb {
	my $class = shift;
	if (ref $class) { confess("Class method used as object method\n") }
	unless (scalar @_ == 1) { confess("Usage: Wrapper::Swipe->set_makeblastdb(COMMAND)\n") }
	$makeblastdbprog = shift;
}


=head3 outdir()

Sets the output directory. Defaults to 'F<.>'.

=cut

sub outdir {
	my $class = shift;
	if (ref $class) { confess("Class method used as object method\n") }
	unless (@_ == 1) { confess("Usage: Wrapper::Swipe->outdir(FILENAME)\n") }
	my $outdir = shift;
}

=head3 evalue_threshold()

Sets or returns the e-value threshold to use for the swipe search. Defaults to 10.

=cut

sub evalue_threshold {
	my $class = shift;
	if (ref $class) { confess("Class method used as object method\n") }
	if (scalar(@_) == 0) { return $evalue_threshold };
	if (scalar(@_) >  1) { confess("Usage: Wrapper::Swipe->evalue_threshold(N)\n") }
	$evalue_threshold = shift(@_);
	unless ($evalue_threshold =~ /^[0-9]+(\.[0-9]+)?([eE][-+]?[0-9]+)?$/) { confess("Invalid argument (must be integer, float or exponential): $evalue_threshold\n") }
}

=head3 score_threshold()

Sets or returns the score threshold to use for the swipe search. Defaults to 10.

=cut

sub score_threshold {
	my $class = shift;
	if (ref $class) { confess("Class method used as object method\n") }
	if (scalar(@_) == 0) { return $score_threshold };
	if (scalar(@_) >  1) { confess("Usage: Wrapper::Swipe->score_threshold(N)\n") }
	$score_threshold = shift(@_);
	unless ($score_threshold =~ /^[0-9]+(\.[0-9]+)?([eE][-+]?[0-9]+)?$/) { confess("Invalid argument (must be integer, float or exponential): $score_threshold\n") }
}

=head2 num_threads

Sets or returns the number of CPU threads to use for multithreaded searching. Defaults to 1.

=cut

sub num_threads {
	my $class = shift;
	if (ref $class) { confess("Class method used as object method\n") }
	if (scalar(@_) == 0) { return $num_threads }
	if (scalar(@_) >  1) { confess("Usage: Wrapper::Swipe->num_threads(N)\n") }
	$num_threads = shift(@_);
	unless ($num_threads =~ /^[0-9]+$/) { confess("Invalid argument (must be integer): $num_threads\n") }
}

=head3 max_hits()

Sets the maximum number of hits to be returned. Defaults to 100.

=cut

sub max_hits {
	my $class = shift;
	if (ref $class) { confess("Class method used as object method\n") }
	unless (@_ == 1) { confess("Usage: Wrapper::Swipe->max_hits(N)\n") }
	$max_hits = shift;
	unless ($max_hits =~ /^[0-9]+$/) { confess("Invalid argument (must be integer): $max_hits\n") }
}

=head1 OBJECT METHODS

=head3 db()

Returns the BLAST database that was selected upon creating a new object

=cut

sub db {
	my $self = shift;
	unless ($self->{'db'}) { confess("I do not have a BLAST db\n") }
	return $self->{'db'};
}

=head3 search()

performs the SWIPE search with the command line set via searchprog(), using
the F<queryfile> on the BLAST DB previously set via db() and with respect to
the e-value set via evalue_threshold().

=cut

sub search {
	my $self = shift;
	unless (scalar @_ == 2) { confess("Usage: Wrapper::Swipe->search(FILE, OUTFILE)\n") }
	my $queryfile = shift;
	my $outfile = shift;
	my $db = $self->db;
	# return right away if this search has been conducted before
	if (-e $outfile) {
		print STDERR "SWIPE output file exists in '$outfile'\n" if $debug;
		$self->resultfile($outfile);
		return $self;
	}
	else {
		print STDERR "SWIPE output file does not exist in '$outfile', conducting new search\n" if $debug;
		# use outfmt 9 for comment lines
		# the columns and their order is different from blast
		my $cmd = qq($searchprog --outfmt 9 --evalue $evalue_threshold --min_score $score_threshold --num_threads $num_threads --db $db --query "$queryfile" --out "$outfile");

		# do the search or die
		print STDERR "\n$cmd\n\n"
			if $debug;
		croak "Fatal: SWIPE search failed: $!\n"
			if system("$cmd");

		# store the resultfile path
		$self->{'resultfile'} = $outfile;
		return $self;
	}
}

=head3 resultfile()

Sets or returns the SWIPE result filename as a path.

=cut

sub resultfile {
	my $self = shift;
	if (scalar @_ == 1) {
		$self->{'resultfile'} = shift;
		return 1;
	}
	return $self->{'resultfile'};
}

=head3 hitcount()

Returns the number of SWIPE hits

=cut

sub hitcount {
	my $self = shift;
	if ($self->{'hitcount'}) {
		return $self->{'hitcount'};
	}
	$self->{'hitcount'} = scalar @{$self->hits_arrayref};
	return $self->{'hitcount'};
}

=head3 result()

Returns the SWIPE result as an array of strings, just as it is in the output file. 

=cut

sub result {
	my $self = shift;
	if ($self->{'result'}) {
		return $self->{'result'};
	}
	my $fh = IO::File->new($self->resultfile);
	$self->{'result'} = [ <$fh> ];
	$fh->close;
	chomp @{$self->{'result'}};
	# remove the four comment lines
	splice @{$self->{'result'}}, 0, 4;
	return $self->{'result'};
}

=head3 hits_arrayref()

Returns the SWIPE result as an array of arrays.

=cut

sub hits_arrayref {
	my $self = shift;
	if ($self->{'hits'}) {
		return $self->{'hits'};
	}
	$self->{'hits'} = [ ];
	my $query = '';
	foreach (@{$self->result}) {
		my @line = split(/\s+/);
		$line[1] =~ s/^lcl\|//;  # remove the "lcl|" prefix from the id
		push(@{$self->{'hits'}}, {
			# the columns are: query, subject, perc identity, alignment length, mismatches, gaps, q.start, q.end, s.start, s.end, e-value, bit score
			'query'  => $line[0],
			'target' => $line[1],
			'score'  => $line[-1],
			'evalue' => $line[-2],
			'start'  => $line[-4],
			'end'    => $line[-3],
		});
	}
	return $self->{'hits'};
}

=head1 INDEPENDENT FUNCTIONS

=head3 makeblastdb()

Create a BLAST database from a (fasta) file.

=cut

sub makeblastdb {
	croak('Usage: makeblastdb($infile, $outfile, $title)' . "\n") unless scalar(@_) == 3;
	my $infile = shift or croak('Usage: makeblastdb($infile, $outfile, $title)' . "\n");
	my $outfile = shift or croak('Usage: makeblastdb($infile, $outfile, $title)' . "\n");
	my $title = shift or croak('Usage: makeblastdb($infile, $outfile, $title)' . "\n");
	my @cmd = qw($makeblastdbprog -in $infile -out $outfile -input_type fasta -title $title);
}

# return true
1;

=head1 AUTHOR

Written by Malte Petersen <mptrsen@uni-bonn.de>.

=head1 COPYRIGHT

Copyright (c) 2012 by Malte Petersen. All rights reserved.

This program is free software; you may redistribute and/or modify it under the
same terms as Orthograph itself.

=cut

