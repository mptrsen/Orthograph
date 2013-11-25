package Orthograph::Ortholog;

use strict;
use warnings;

sub new {
	my ($class, $hmmhit) = @_;
	my $self = $hmmhit;
	return bless ($self, $class);
}

sub blast_target {
	my $self = shift @_;
	if (scalar @_ > 0) {
		$self->{'blast_target'} = shift @_;
		return 1;
	}
	return $self->{'blast_target'};
}

sub blast_query {
	my $self = shift @_;
	if (scalar @_ > 0) {
		$self->{'blast_query'} = shift @_;
		return 1;
	}
	return $self->{'blast_query'};
}

sub as_hashref {
	my $self = shift;
	my $r = {};
	while (my ($k, $v)) {
		$r->{$k} = $v;
	}
	return $r;
}

'This line intentionally left false';
