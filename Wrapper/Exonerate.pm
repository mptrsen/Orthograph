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

my $verbose    = 0;
my $debug      = 0;
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
	my $exonerate_ryo = '>cdna %tcb %tce\n%tcs>aa %qab %qae\n%qas';

	# the complete command line
	my $exonerate_cmd = qq($searchprog --bestn 1 --score $score_threshold --ryo '$exonerate_ryo' --model $exonerate_model --querytype protein --targettype dna --verbose 0 --showalignment no --showvulgar no $exhaustive --query $queryfile --target $targetfile > $outfile);
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

sub aa_sequence {
	my $self = shift;
	unless ($self->{'aa_sequence'}) { $self->get_orf() }
	unless ($self->{'aa_sequence'}) { return undef }
	$self->{'aa_sequence'} =~ s/\s//g;
	return $self->{'aa_sequence'};
}

sub cdna_sequence {
	my $self = shift;
	unless ($self->{'cdna_sequence'}) { $self->get_orf() }
	unless ($self->{'cdna_sequence'}) { return undef }
	$self->{'cdna_sequence'} =~ s/\s//g;
	return $self->{'cdna_sequence'};
}

sub get_orf {
	my $self = shift;
	# if there was no result
	if (-z $self->{'resultfile'}) { return undef }

	my $cdna_seq = '';
	my $aa_seq   = '';

	# otherwise, continue
	my $orf_list = slurp_orfs_from_fasta($self);

	# start is the beginning of the first orf,
	# end the end of the last
	$self->{'cdna_start'} = $orf_list->[0]->{'cdna_start'};
	$self->{'cdna_end'} = $orf_list->[-1]->{'cdna_end'};

	for (my $i = 0; $i < scalar @$orf_list; $i++) {
		$cdna_seq .= $orf_list->[$i]->{'cdna_seq'};
		$aa_seq .= $orf_list->[$i]->{'aa_seq'};
		# don't do this for the last orf
		if ($i < scalar(@$orf_list) - 1) {
			# the cdna sequence
			my $num_missing = $orf_list->[$i+1]->{'cdna_start'} - $orf_list->[$i]->{'cdna_end'} - 1;
			my $indel = lc(substr($self->{'target'}->{'sequence'}, $orf_list->[$i]->{'cdna_end'}, $num_missing));
			# append gap characters until codon filled
			# no need to do that for the aa sequence
			while (length($indel) % 3 != 0) { $indel .= '-' }
			# append indel sequence
			$cdna_seq .= $indel;

			# the aa sequence accordingly
			$num_missing = $orf_list->[$i+1]->{'aa_start'} - $orf_list->[$i]->{'aa_end'} - 1;
			$indel = lc(substr($self->{'query'}->{'sequence'}, $orf_list->[$i]->{'aa_end'}, $num_missing));
			# append indel
			$aa_seq .= $indel;
		}
	}
	
	$self->{'cdna_sequence'} = $cdna_seq;
	return $self;
	
}


sub slurp_orfs_from_fasta {
	my $self = shift;
	my $list = [ ];
	my $fh = Seqload::Fasta->open($self->{'resultfile'});
	# watch for the order of sequences, they must correspond to the order in
	# the --ryo option in the exonerate call
	while (my ($h_cdna, $s_cdna) = $fh->next_seq()) {
		# there will always be two sequences per alignment, so fetching the next one is ok
		my ($h_aa, $s_aa) = $fh->next_seq();
		# but test anyway
		croak "Fatal: no sequence found for ORF (aa)\n" unless $h_aa and $s_aa;
		my @fields = split ' ', $h_cdna;
		my $cdna_start = $fields[1];
		my $cdna_end   = $fields[2];
		@fields = split ' ', $h_aa;
		my $aa_start   = $fields[1];
		my $aa_end     = $fields[2];
		push @$list, {
			'cdna_start' => $cdna_start,
			'cdna_end'   => $cdna_end,
			'aa_start'   => $aa_start,
			'aa_end'     => $aa_end,
			'cdna_seq'   => $s_cdna,
			'aa_seq'     => $s_aa,
		};
	}
	undef $fh;
	return $list;
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
	if ($debug) {
		printf "Wrote this sequence to Fasta file '%s':\n>%s\n%s\n", $fh, $header, $sequence;
	}
	return $fh;
}

