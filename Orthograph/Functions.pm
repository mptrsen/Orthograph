#--------------------------------------------------
# This file is part of Orthograph.
# Copyright 2014 Malte Petersen <mptrsen@uni-bonn.de>
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
package Orthograph::Functions;

=head1 Functions common to all Orthograph tools

=cut

use strict;
use warnings;
use autodie;

=head2 file2arrayref

Reads a file, line by line. Removes linebreaks and returns an arrayref.

=cut

sub file2arrayref {
	my $f = shift;
	open my $fh, '<', $f or die "Fatal: Could not open file '$f': $!\n";
	my $l = [ <$fh> ];
	chomp @$l;
	return $l;
}

sub touch {
	my $now = time;
	my $file = shift @_;
	utime($now, $now, $file)
	|| open(my $fn, '>>', $file)
	|| croak("Couldn't touch file $file: $!\n");
	return 1;
}

1;
