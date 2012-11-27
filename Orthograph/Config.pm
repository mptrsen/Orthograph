package Orthograph::Config;
use strict;
use warnings;

use File::Spec;             
use FindBin;                # locate the dir of this script during compile time
use lib $FindBin::Bin;      # $Bin is the directory of the original script
require Exporter;
our @ISA = qw( Exporter );
our @EXPORT_OK = qw( $config );

my $program_name = 'Orthograph';
our $config = &getconfig; 

exit(1) unless $config;


sub getconfig {
	my $configfile = File::Spec->catfile($FindBin::Bin, lc($program_name) . '.conf');
	# in case the user tells us to use a different one with -c
	$configfile = &get_configfile($configfile);

	# parse if exists
	if (-e $configfile) {
		print "Parsing config file '$configfile'.\n";
		$config = &parse_config($configfile);
	}#}}}
	else {
		warn "Fatal: Config file '$configfile' not found!\n";
		return 0;
	}

	return $config;
}

=head2 get_configfile

mini argument parser to get the config file name

=cut

# mini argument parser for the configfile
sub get_configfile {
	my $configfile = shift(@_);
	# every argument
	for (my $i = 0; $i < scalar @ARGV; ++$i) {
		# is this '-c'?
		if ($ARGV[$i] =~ /-c\b/) {
			# does the next one begin with a hyphen?
			if ($ARGV[$i+1] !~ /^-/) {
				$configfile = $ARGV[$i+1];
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
