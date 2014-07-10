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
use File::Basename;
use File::Path qw( make_path );	# this also uses File::Spec

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

=head2 touch

Touches a file (updates its access times, creates an empty file if it doesn't exist)

=cut

sub touch {
	my $now = time;
	my $file = shift @_;
	my $dir = dirname($file);
	unless (-e $file) { 
		unless (-e $dir) { make_path($dir) }
		open my $fh, '>>', $file or die "Fatal: Couldn't open file '$file': $!\n";
		close $fh or warn "Warning: Couldn't close file '$file': $!\n";
		return 1;
	}
	utime $now, $now, $file or die $! or croak("Fatal: Couldn't touch file $file: $!\n");
	return 1;
}

1;
