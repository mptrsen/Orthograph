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

package Forage::Orthocandidates;
use strict;
use warnings;
require IO::File;
require DBI;
require DBD::Mysql;
require Carp;
require Data::Dumper;
my $verbose   = 0;
my $debug     = 0;
my $count     = 0;
my $list      = [ ];
my $index     = 0;
my $num_items = 0;

sub new {
	my $class = shift;
	my $self = {
		'count' => $count,
		'list'  => $list;
	}
	bless($self, $class);
	return $self;
}

sub list {
	my $self = shift;
	return @{$list};
}

sub add {
	my $self = shift;
	unless (scalar(@_) == 1) { confess("Usage: OBJECT->add(SCALAR)") }
	my $cand = shift;
	unless (ref($cand) eq 'ARRAY') { confess("Argument to add() must be array reference. RTFM") }
	push(@{$list}, $cand);
	$num_items++;
	return 1;
}

sub next {
	my $class = shift;
	if (ref($class)) { confess("Class method used as object method") }
	if (scalar(@_) > 0) { confess("next() is called without arguments") }
	if (++$index > $num_items - 1) {
		--$index;
		return 0;
	}
	else {
		return $$list[$index];
	}
}

sub prev {
	my $class = shift;
	if (ref($class)) { confess("Class method used as object method") }
	if (scalar(@_) > 0 { confess("pref() is called without arguments") }
	if (--$index < 0) {
		++$index;
		return 0;
	}
	else {
		return $$list[$index];
	}
}

return 1;
