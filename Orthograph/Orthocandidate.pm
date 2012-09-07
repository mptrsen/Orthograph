#--------------------------------------------------
# This file is part of Orthograph.
# Copyright 2012 Malte Petersen <mptrsen@uni-bonn.de>
# 
# Orthograph is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
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

B<Orthograph::Orthocandidate>

=head1 DESCRIPTION

The B<Orthograph::Orthocandidate> module provides a simple, object-oriented
interface to a orthology candidate database, however managed. It provides
methods for adding and removing candidate objects, looping through the list
backwards and forwards as well as dumping the entire list.

Note that this is only a front-end. The backend may use whatever different
means of storage. I plan to write a MySQL database interface so that the RAM
shall not be cluttered during large analyses.

=head1 SYNOPSIS

  use Orthograph::Orthocandidate;

  # create a new candidate object
	my $candidate = Orthograph::Orthocandidate->new($ortholog_id, $digest);

	# add the newly created object to the candidate list
	Orthograph::Orthocandidate->add($candidate);

	# output all objects in the list
	while (my $o = Orthograph::Orthocandidate->next()) {
		printf("%s\t%s\n", $o->orthoid(), $o->digest());
	}

=cut

package Orthograph::Orthocandidate;
use strict;
use warnings;
use IO::File;
use DBI;
use DBD::mysql;
use Carp;
use Data::Dumper;
my $verbose   = 0;
my $debug     = 0;
my $count     = 0;
my $list      = [ ];
my $index     = 0;
my $num_items = 0;

sub new {
	my $class   = shift;
	unless (scalar(@_) == 2) { confess('Usage: Orthograph::Orthocandidate->new($orthoid, $digest)') }
	my ($orthoid, $digest) = splice(@_, 0, 2);
	my $self    = {
		'orthoid' => $orthoid,
		'digest'  => $digest,
	};
	bless($self, $class);
	#Orthograph::Orthocandidates->add($self);	# TODO ??
	return $self;
}

=head1 OBJECT METHODS

=head3 orthoid()

Returns the ortholog id of the current candidate object.

=cut

sub orthoid {
	my $self = shift;
	return $$self{'orthoid'};
}

=head3 digest()

Returns the SHA1 digest of the current candidate object.

=cut

sub digest {
	my $self = shift;
	return $$self{'digest'};
}

=head1 CLASS METHODS

=head3 getindex()

Returns the index of the current candidate object.

=cut

sub getindex {
	my $class = shift;
	if (ref($class)) { confess("Class method used as object method") }
	return $index;
}

=head3 list()

Returns all existing objects in an array of references.

=cut

sub all {
	my $class = shift;
	if (ref($class)) { confess("Class method used as object method") }
	return @{$list};
}

=head3 add(OBJECT)

Adds OBJECT to the list of candidate objects.

=cut

sub add {
	my $class = shift;
	if (ref($class)) { confess("Class method used as object method") }
	unless (scalar(@_) == 1) { confess("Usage: Orthograph::Orthocandidate->add(SCALAR)") }
	my $cand = shift;
	unless (ref($cand)) { confess("Argument to Orthograph::Orthocandidate->add() must be hash reference. RTFM") }
	$$cand{'index'} = $index++;
	push(@{$list}, $cand);
	$num_items++;
	return 1;
}

=head3 cut()

Removes the current candidate object from the list and returns it.

=cut

sub cut {
	my $class = shift;
	if (ref($class)) { confess("Class method used as object method") }
	my $item = splice(@{$list}, $index, 1);
	--$num_items;
	if ($index > $num_items - 1) { $index = $num_items - 1 }
	return $item;
}

=head3 next()

Returns the next candidate object in the list.

=cut

sub next {
	my $class = shift;
	if (ref($class)) { confess("Class method used as object method") }
	if ($index > $num_items - 1) {
		$index = $num_items - 1;
		return 0;
	}
	else {
		return $$list[$index++];
	}
}

=head3 prev()

Returns the previous candidate object in the list.

=cut

sub prev {
	my $class = shift;
	if (ref($class)) { confess("Class method used as object method") }
	if ($index < 0) {
		$index = 0;
		return 0;
	}
	else {
		return $$list[$index--];
	}
}

=head3 rewind()

Sets the list index to 0 (the first element). Useful before looping through the list again.

=cut

sub rewind {
	my $class = shift;
	if (ref($class)) { confess("Class method used as object method") }
	$index = 0;
	return 1;
}

1;

=head1 AUTHOR

Written by Malte Petersen <mptrsen@uni-bonn.de>.

=head1 COPYRIGHT

Copyright (c) 2012 by Malte Petersen. All rights reserved.

This program is free software; you may redistribute and/or modify it under the
same terms as Orthograph itself.

=cut
