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
package Orthograph::Functions;

=head1 Functions common to all Orthograph tools

=cut

use IO::File;

=head2 get_configfile

mini argument parser to get the config file name

=cut

# mini argument parser for the configfile
sub get_configfile {
	my $configfile = shift(@_);
	for (my $i = 0; $i < scalar @ARGV; ++$i) {
		if ($ARGV[$i] =~ /-c\b/) {
			if ($ARGV[$i+1] !~ /^-/) {
				$configfile = $ARGV[$i+1];
				if (!-e $configfile) { die "Fatal: Config file '$configfile' not found!\n" }
			}
			# the file name starts with a hyphen, may be a stray option, so warn the
			# user and don't use this name
			else { warn "Warning: Config file name '$ARGV[$i+1]' starts with a hyphen (-), could be a stray option. Use './$ARGV[$i+1]' if you mean it. Falling back to '$configfile' for now.\n" }
		}
	}
	return $configfile;
}

=head2 Sub: parse_config

Parse a simple, ini-style config file where keys are separated from values by '='. 
Sections are not supported. 
outputdir = /home/foo/bar

=cut 

sub parse_config {#{{{
	my $file = shift;
	my $conf = { };
	my $fh = IO::File->new($file) or die "Fatal: Could not open config file '$file'\: $!\n";

	while (my $line = $fh->getline()) {
		next if $line =~ /^\s*$/; # skip empty lines
		next if $line =~ /^\s*#/; # skip comment lines starting with '#'
		if ($line !~ /^\s*\w+\s*=\s*[\/]?\w+/) {
			die "Fatal: Invalid format in line $. of config file $file:\n$line\n"
		}
		
		# split by '=' producing a maximum of two items
		my ($key, $val) = split('=', $line, 2);

		foreach ($key, $val) {
		  s/\s+$//; # remove all trailing whitespace
		  s/^\s+//; # remove all leading whitespace
		}

		die "Fatal: Configuration option '$key' defined twice in line $. of config file '$file'\n"
		  if defined $conf->{$key};
		$conf->{$key} = $val;
	}
	close($fh);
	return $conf;
}#}}}


1;
