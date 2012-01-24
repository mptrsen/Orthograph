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
package Forage::Blast;
use strict;
use warnings;
use File::Basename; # basename of files
use IO::File; # object-oriented access to files
use Carp; # extended dying functions
my $verbose = 0;
my $blastoutdir = File::Spec->catdir('.');
my $blastprog = 'blastp';

sub new {
	my ($class, $db, $queryfile) = @_;

	my $self = {
		'db'       => $db,
		'query'    => $queryfile,
		'hitcount' => undef;
	};

	bless($self, $class);
	return $self;
}

sub blastp {
	my $self = shift;
	unless (scalar @_ == 1) { confess "Usage: Forage::Blast->blastp(FILE)\n" }
	my $queryfile = shift;
	# TODO complete the blast search
}
1;
