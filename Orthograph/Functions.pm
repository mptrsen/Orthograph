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
use File::Which qw( which );    # test executability of programs in PATH

=head2 file2arrayref

Reads a file, line by line. Removes linebreaks and returns an arrayref.

=cut

sub file2arrayref {
	my $f = shift;
	open my $fh, '<', $f or die "Fatal: Could not open file '$f': $!\n";
	my $f_content = do { local $/; <$fh> }; # slurp entire file
	$f_content =~ s/\r\n?/\n/g; # convert to unix-style line breaks
	my $l = [ split /\n/, $f_content ];
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
		or croak "Fatal: Couldn't open dir '$dir': $!";
	foreach (readdir($dirh)) {
		next if $_ =~ /^\.\.?$/;
		unlink(File::Spec->catfile($dir, $_))
			or warn "Warning: Could not delete file '" . File::Spec->catfile($dir, $_) . "': $!\n";
	}
	closedir($dirh)
		or warn "Warning: Couldn't close dir '$dir': $!";
	return 1;
}#}}}

# Sub: usage
# Returns a usage message based on the config hash from the Orthograph::Config
# module
# Arguments: hashref
# Returns: scalar string (usage message)
sub print_usage {
	my $config = shift;
	print "Usage: $0 [OPTIONS] INPUTFILE\n";
	print "Options:\n";
	foreach my $opt (sort { $a cmp $b } keys %$config) {
		next if $opt =~ /db_/;
		print "\t--$opt\n";
	}
	print "See the documentation for a description.\n";
	return 1;
}

# test whether a (program) file exists and is executable
sub program_exists {
	my $path = shift;
	# remove optional arguments if they exist. this requires that arguments start
	# with a - and are separated by a space from the command
	$path =~ s/\s+-.+//; 
	if (which($path))              { return 1     }
	elsif (-f $path and -x $path)  { return 1     }
	else                           { return undef }
}

# test whether the program versions are OK
sub fastatranslate_version_ok {
	my $program = shift;
	my $ret = [ `$program --version` ];
	$$ret[0] =~ /version ([0-9]+)\.([0-9]+)/;
	my $version = "$1.$2";
	if ("$1.$2" >= 2.2) { return ($version, 1) }
	return ($version, 0);
}

sub alignment_program_version_ok {
	my $program = shift;
	$program =~ s/\s.+$//;
	my $ret = [ `$program --version 2>&1` ];
	$$ret[0] =~ /v([0-9.]+)/;
	my $version = $1;
	#if ($version eq '7.023b') { return $version, 1 }
	my @nums = split /\./, $version;
	if ($nums[0] > 7 or ($nums[0] == 7 and $nums[1] >= 273)) { return ($nums[0] . (defined $nums[1] ? '.' . $nums[1] : ''), 1)  };
	return ($version, 0);
}

sub hmmbuild_version_ok {
	my $program = shift;
	my $ret = [ `$program -h` ];
	$$ret[1] =~ /HMMER ([0-9.]+)/;
	my $version = $1;
	my @nums = split /\./, $version;
	if ($nums[0] > 3 or ($nums[0] >= 3 and $nums[1] >= 1)) { return ($version, 1) };
	return ($version, 0);
}

sub hmmsearch_version_ok {
	my $program = shift;
	my $ret = [ `$program -h` ];
	$$ret[1] =~ /HMMER ([0-9.]+)/;
	my $version = $1;
	my @nums = split /\./, $version;
	if ($nums[0] > 3 or ($nums[0] >= 3 and $nums[1] >= 1)) { return ($version, 1) };
	return ($version, 0);
}

sub makeblastdb_version_ok {
	my $program = shift;
	my $ret = [ `$program -version` ];
	$$ret[0] =~ /makeblastdb: ([0-9.]+)/;
	my $version = $1;
	my @nums = split /\./, $version;
	if ($nums[0] > 2 or ($nums[0] == 2 and $nums[1] > 2) or ($nums[0] == 2 and $nums[1] == 2 and $nums[2] >= 28)) { return ($version, 1) }
	return ($version, 0);
}

sub blastp_version_ok {
	my $program = shift;
	my $ret = [ `$program -version` ];
	$$ret[0] =~ /blastp: ([0-9.]+)/;
	my $version = $1;
	my @nums = split /\./, $version;
	if ($nums[0] > 2 or ($nums[0] == 2 and $nums[1] > 2) or ($nums[0] == 2 and $nums[1] == 2 and $nums[2] >= 28)) { return ($version, 1) }
	return ($version, 0);
}

sub swipe_version_ok {
	my $program = shift;
	my $ret = [ `$program -h` ];
	$$ret[0] =~ /SWIPE ([0-9.]+)/;
	my $version = $1;
	my @nums = split /\./, $version;
	if ($nums[0] > 2 or ($nums[0] == 2 and $nums[1] > 0) or ($nums[0] == 2 and $nums[1] == 0 and $nums[2] >= 12)) { return ($version, 1) }
	return ($version, 0);
}

sub exonerate_version_ok {
	my $program = shift;
	my $ret = [ `$program --version` ];
	$$ret[0] =~ /version ([0-9.]+)/;
	my $version = $1;
	my @nums = split /\./, $version;
	if ($nums[0] > 2 or ($nums[0] == 2 and $nums[1] > 2) or ($nums[0] == 2 and $nums[1] == 2 and $nums[2] >= 0)) { return ($version, 1) }
	return ($version, 0);
}

sub test_dependencies {

	my ($translate_program, $alignment_program, $hmmbuild_program, $makeblastdb_program, $hmmsearch_program, $blast_program, $exonerate_program, $swipe_program) = @_;

	my $version = 0;
	my $ok = 0;

	# test whether the programs exist where specified
	program_exists($translate_program)   or die "Fatal: Fastatranslate not executable at '$translate_program'. Verify path and/or permissions.\n";
	program_exists($exonerate_program)   or die "Fatal: Exonerate not executable at '$exonerate_program'. Verify path and/or permissions.\n";
	program_exists($alignment_program)   or die "Fatal: Alignment program not executable at '$alignment_program'. Verify path and/or permissions.\n";
	program_exists($hmmbuild_program)    or die "Fatal: HMMbuild not executable at '$hmmbuild_program'. Verify path and/or permissions.\n";
	program_exists($makeblastdb_program) or die "Fatal: Makeblastdb not executable at '$hmmbuild_program'. Verify path and/or permissions.\n";
	program_exists($hmmsearch_program)   or die "Fatal: HMMsearch not executable at '$hmmsearch_program'. Verify path and/or permissions.\n";
	program_exists($blast_program)       or die "Fatal: BLASTP not executable at '$blast_program'. Verify path and/or permissions.\n";
	# swipe is optional
	if ($swipe_program) {
		program_exists($swipe_program)       or die "Fatal: SWIPE not executable at '$swipe_program'. Verify path and/or permissions.\n";
	}

	# test whether the versions are correct
	($version, $ok) = fastatranslate_version_ok($translate_program);
	$ok or die "Fatal: fastatranslate failed version check. Requires at least version 2.2.0. You have version '$version'\n";
	print "OK: '$translate_program' version $version\n";
	($version, $ok) = exonerate_version_ok($exonerate_program);
	$ok or die "Fatal: exonerate failed version check. Requires at least version 2.2.0. You have version '$version'\n";
	print "OK: '$exonerate_program' version $version\n";
	($version, $ok) = alignment_program_version_ok($alignment_program);
	$ok or die "Fatal: mafft failed version check. Requires at least version 7.273. You have version '$version'\n";
	print "OK: '$alignment_program' version $version\n";
	($version, $ok) = hmmbuild_version_ok($hmmbuild_program);
	$ok or die "Fatal: hmmbuild failed version check. Requires at least version 3.1b1. You have version '$version'\n";
	print "OK: '$hmmbuild_program' version $version\n";
	($version, $ok) = hmmsearch_version_ok($hmmsearch_program);
	$ok or die "Fatal: hmmsearch failed version check. Requires at least version 3.1b1. You have version '$version'\n";
	print "OK: '$hmmsearch_program' version $version\n";
	($version, $ok) = makeblastdb_version_ok($makeblastdb_program);
	$ok or die "Fatal: makeblastdb failed version check. Requires at least version 2.2.28. You have version '$version'\n";
	print "OK: '$makeblastdb_program' version $version \n";
	($version, $ok) = blastp_version_ok($blast_program);
	$ok or die "Fatal: blastp failed version check. Requires at least version 2.2.28. You have version '$version'\n";
	print "OK: '$blast_program' version $version \n";
	if ($swipe_program) {
		($version, $ok) = swipe_version_ok($swipe_program);
		$ok or die "Fatal: swipe failed version check. Requires at least version 2.0.12. You have version '$version'\n";
		print "OK: '$swipe_program' version $version \n";
	}

	return 1;

}


1;
