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
use Carp;
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
	utime $now, $now, $file or croak("Fatal: Couldn't touch file $file: $!\n");
	return 1;
}

=head2 progress_bar

Prints a self-overwriting, wget-style progress bar.
This is not written to the log file so it doesn't get cluttered.
Arguments: scalar int so-far, scalar int total, scalar int width, scalar char "what-to-use-as-char"

=cut

sub progress_bar {#{{{
	my ($got, $total, $width, $char) = @_;
	$width ||= 25;
	$char ||= '=';
	my $num_width = length($total);
	local $| = 1;
	printf("|%-${width}s| Progress: %${num_width}s of %s (%.2f%%)\r",	
		$char x (($width-1)*$got/$total) . '>',
		$got,
		$total,
		100 * $got / $total
	);
	local $| = 0;
	return 1;
}#}}}

# Sub: file_is_empty
# tests whether a file is empty (i.e., contains nothing or only empty lines)
# Expects: scalar string path to file
# Returns: True if file is empty, false otherwise
sub file_is_empty {
	my $file = shift;
	my $fh = IO::File->new($file);
	while (<$fh>) {
		/^\s*$/ ? next : return 0;	# skip empty lines
	}
	return 1;
}

# Sub: makedir
# Creates a directory with parent directories if it does not already exist
# Expects: scalar string dirname
# Returns: True if successful
sub makedir {#{{{
  my $dir = shift;
  if (-e $dir and not -d $dir) {
    print "Fatal: $dir exists, but is not a directory! This will most likely lead to problems later. Exiting.\n" and exit(1);
  }
  elsif (!make_path($dir, { verbose => 0 })) {
    warn "Warning: Could not create $dir: $!\n";
    return 0;
  }
  return 1;
}#}}}

# Sub: cleardir
# Empties a directory of non-dotfiles
# Arguments: scalar string dirname
sub cleardir {#{{{
	my $dir = shift;
	opendir(my $dirh, $dir)
		or croak "Fatal: Couldn't open dir $dir: $!";
	foreach (readdir($dirh)) {
		next if $_ =~ /^\.\.?$/;
		unlink(File::Spec->catfile($dir, $_))
			or warn "Warning: Could not delete file " . File::Spec->catfile($dir, $_) . ": $!";
	}
	closedir($dirh)
		or warn "Fatal: Couldn't close dir $dir: $!";
	return 1;
}#}}}

1;
