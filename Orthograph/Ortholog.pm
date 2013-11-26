package Orthograph::Ortholog;

use strict;
use warnings;

sub new {
	my ($class, $hmmhit) = @_;
	my $self = $hmmhit;
	return bless ($self, $class);
}

sub hmm_target {
	my $self = shift @_;
	if (scalar @_ > 0) {
		$self->{'hmm_target'} = shift @_;
		return 1;
	}
	return $self->{'hmm_target'};
}

sub hmm_name {
	my $self = shift @_;
	if (scalar @_ > 0) {
		$self->{'hmm_name'} = shift @_;
		return 1;
	}
	return $self->{'hmm_name'};
}

sub hmm_start {
	my $self = shift @_;
	if (scalar @_ > 0) {
		$self->{'hmm_start'} = shift @_;
		return 1;
	}
	return $self->{'hmm_start'};
}

sub hmm_end {
	my $self = shift @_;
	if (scalar @_ > 0) {
		$self->{'hmm_end'} = shift @_;
		return 1;
	}
	return $self->{'hmm_end'};
}

sub env_start {
	my $self = shift @_;
	if (scalar @_ > 0) {
		$self->{'env_start'} = shift @_;
		return 1;
	}
	return $self->{'env_start'};
}

sub env_end {
	my $self = shift @_;
	if (scalar @_ > 0) {
		$self->{'env_end'} = shift @_;
		return 1;
	}
	return $self->{'env_end'};
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
