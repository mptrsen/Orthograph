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
package Orthograph::Ortholog;

use strict;
use warnings;
use Wrapper::Exonerate;

=head2 new

Constructor. Returns blessed hashref.

=cut

sub new {
	my ($class, $hit) = @_;
	my $self = $hit;
	return bless ($self, $class);
}

=head2 hmm_target

Sets or gets the hmm target. If called with an argument, sets the hmm target and returns 1, otherwise, returns the hmm target as a scalar.

=cut

sub hmm_target {
	my $self = shift @_;
	if (scalar @_ > 0) {
		$self->{'hmm_target'} = shift @_;
		return 1;
	}
	return $self->{'hmm_target'};
}

=head2 hmm_name

Sets or gets the hmm name, i.e., the hmm id. If called with an argument, sets the hmm name and returns 1, otherwise, returns the hmm name.

=cut

sub hmm_name {
	my $self = shift @_;
	if (scalar @_ > 0) {
		$self->{'hmm_name'} = shift @_;
		return 1;
	}
	return $self->{'hmm_name'};
}

=head2 hmm_start

Sets or gets the start of the hmm target sequence on hmm level. If called with an argument, sets the start and returns 1, otherwise, returns the start.

=cut

sub hmm_start {
	my $self = shift @_;
	if (scalar @_ > 0) {
		$self->{'hmm_start'} = shift @_;
		return 1;
	}
	return $self->{'hmm_start'};
}

=head2 hmm_end

Sets or gets the end of the hmm target sequence on hmm level. If called with an argument, sets the end and returns 1, otherwise, returns the end.

=cut

sub hmm_end {
	my $self = shift @_;
	if (scalar @_ > 0) {
		$self->{'hmm_end'} = shift @_;
		return 1;
	}
	return $self->{'hmm_end'};
}

=head2 env_start

Sets or gets the start of the hmm target sequence on transcript level. If called with an argument, sets the start and returns 1, otherwise, returns the start.

=cut

sub env_start {
	my $self = shift @_;
	if (scalar @_ > 0) {
		$self->{'env_start'} = shift @_;
		return 1;
	}
	return $self->{'env_start'};
}

=head2 env_end

Sets or gets the end of the hmm target sequence on transcript level. If called with an argument, sets the end and returns 1, otherwise, returns the end.

=cut

sub env_end {
	my $self = shift @_;
	if (scalar @_ > 0) {
		$self->{'env_end'} = shift @_;
		return 1;
	}
	return $self->{'env_end'};
}

=head2 blast_target

Sets or gets the blast target sequence identifier. If called with an argument, sets the target and returns 1, otherwise, returns the target.

=cut

sub blast_target {
	my $self = shift @_;
	if (scalar @_ > 0) {
		$self->{'blast_target'} = shift @_;
		return 1;
	}
	return $self->{'blast_target'};
}

=head2 blast_query

Sets or gets the blast query sequence identifier. If called with an argument, sets the query and returns 1, otherwise, returns the query.

=cut

sub blast_query {
	my $self = shift @_;
	if (scalar @_ > 0) {
		$self->{'blast_query'} = shift @_;
		return 1;
	}
	return $self->{'blast_query'};
}

=head2 as_hashref

Returns the object as a hash reference. The hash has all the keys that have been defined so far via the constructor or the methods, i.e., its structure is as follows:

  {
    'hmm_name' => 'EOG50000K',
    'hmm_target' => '7ee8467b316db5c07b2c76f865376035',
    'hmm_start => 356,
    'hmm_end' => 1282,
    'blast_query' => '7ee8467b316db5c07b2c76f865376035',
    'blast_target' => 248,
    ...
	}

=cut

sub as_hashref {
	my $self = shift;
	my $r = {};
	while (my ($k, $v)) {
		$r->{$k} = $v;
	}
	return $r;
}

=head2 real_header

Sets or returns the real header of the orthologous transcript sequence. To be removed.

=cut

sub real_header {
	my $self = shift;
	if (scalar @_ > 0) {
		$self->{'real_header'} = shift;
		return 1;
	}
	return $self->{'real_header'};
}


=head crop_to_hmm

Crops the ORF sequence including coordinates to the HMM region of the alignment.

=cut

sub crop_to_orf {
	my $self = shift;
	unless ($self->orf_aa_start() and $self->orf_aa_end()) {
		return 0;
	}
	else {
		
	}
}

'This line intentionally left false';
