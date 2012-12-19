package Orthograph::Sequence;

use strict;
use warnings;

use Data::Dumper;

sub new {
	my ($class, $hashref) = @_;
	my $self = {
		'header'    => $hashref->{'header'},
		'sequence'  => $hashref->{'sequence'},
		'digest'    => $hashref->{'digest'},
		'blast_hit' => $hashref->{'blast_hit'},
		'start'     => $hashref->{'start'},
		'end'       => $hashref->{'end'},
		'reftaxon'  => $hashref->{'reftaxon'},
	}
	bless $self, $class;
	return $self;
}

=head1 Object Methods

=head2 header

Get or set the header. Any arguments beyond the first will be ignored.

=cut

sub header {
	my $self = shift;
	if (scalar @_ > 0) {
		my $thing = shift @_;
		unless (ref $thing) {
			$self->{'header'} = $thing;
			return 1;
		}
	}
	else { return $self->{'header'} }
}

=head2 sequence

Get or set the sequence. Any arguments beyond the first will be ignored.

=cut

sub sequence {
	my $self = shift;
	if (scalar @_ > 0) {
		my $thing = shift @_;
		unless (ref $thing) {
			$self->{'sequence'} = $thing;
			return 1;
		}
	}
	else { return $self->{'sequence'} }
}

=head2 hmmhit

Get or set the hmmhit. Any arguments beyond the first will be ignored.

=cut

sub hmmhit {
	my $self = shift;
	if (scalar @_ > 0) {
		my $thing = shift @_;
		unless (ref $thing) {
			$self->{'hmmhit'} = $thing;
			return 1;
		}
	}
	else { return $self->{'hmmhit'} }
}

=head2 blast_hit

Get or set the blast_hit. Any arguments beyond the first will be ignored.

=cut

sub blast_hit {
	my $self = shift;
	if (scalar @_ > 0) {
		my $thing = shift @_;
		unless (ref $thing) {
			$self->{'blast_hit'} = $thing;
			return 1;
		}
	}
	else { return $self->{'blast_hit'} }
}

=head2 start

Get or set the start. Any arguments beyond the first will be ignored.

=cut

sub start {
	my $self = shift;
	if (scalar @_ > 0) {
		my $thing = shift @_;
		unless (ref $thing) {
			$self->{'start'} = $thing;
			return 1;
		}
	}
	else { return $self->{'start'} }
}

=head2 end

Get or set the end. Any arguments beyond the first will be ignored.

=cut

sub end {
	my $self = shift;
	if (scalar @_ > 0) {
		my $thing = shift @_;
		unless (ref $thing) {
			$self->{'end'} = $thing;
			return 1;
		}
	}
	else { return $self->{'end'} }
}

=head2 reftaxon

Get or set the reftaxon. Any arguments beyond the first will be ignored.

=cut

sub reftaxon {
	my $self = shift;
	if (scalar @_ > 0) {
		my $thing = shift @_;
		unless (ref $thing) {
			$self->{'reftaxon'} = $thing;
			return 1;
		}
	}
	else { return $self->{'reftaxon'} }
}
