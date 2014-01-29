#--------------------------------------------------
# This file is part of Orthograph.
# Copyright 2013 Malte Petersen <mptrsen@uni-bonn.de>
# 
# Orthograph is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
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
package Orthograph::Config;
use strict;
use warnings;

use Data::Dumper;
use File::Basename;
use File::Spec;             
use FindBin;                # locate the dir of this script during compile time
use Getopt::Long;           # parse command line arguments
use lib $FindBin::Bin;      # $Bin is the directory of the original script
require Exporter;
our @ISA = qw( Exporter );
our @EXPORT_OK = qw( $config );

my $program_name = 'Orthograph';
my $configfile = File::Spec->catfile($FindBin::Bin, lc($program_name) . '.conf');
our $config = &getconfig; 

#--------------------------------------------------
# # Get command line options. These may override variables set in the config file.
#-------------------------------------------------- 
#{{{
GetOptions( $config,
	'continue',
	'create',
	'delete-ogs=s',
	'delete-set=s',
	'destroy',
	'evalue-bin-size=i',
	'list-ests|le',
	'list-ogs|lo',
	'list-sets|ls',
	'list-taxa|lt',
	'load-ogs-nucleotide=s',
	'load-ogs-peptide=s',
	'prepare',
  'aaoutdir=s',
  'alignment-program=s',
	'backup',
  'backup-extension',
  'blast-evalue-threshold',
  'blast-evalue-threshold=f',
  'blast-max-hits=i',
  'blast-score-threshold=i',
  'blastoutdir=s',
  'blastp-output-dir=s',
  'clear-database!',
  'clear-files!',
  'configfile|c=s',
  'database-engine=s',
  'db-prefix=s',
  'debug|d',
  'estfile',
  'estfile|E=s',
	'header-separator=s',
  'hmmsearch-evalue-threshold',
  'hmmsearch-evalue-threshold=f',
  'hmmsearch-output-dir',
  'hmmsearch-score=i',
  'hmmsearch-score-threshold',
  'hmmsearchprog=s',
  'logfile|log=s',
	'make-set',
  'max-blast-searches=i',
  'mysql-database=s',
  'mysql-password=s',
  'mysql-server=s',
  'mysql-username=s',
  'db_table_aaseqs',
  'db_table_blast',
  'db_table_blastdbs',
  'db_table_ests',
  'db_table_hmmsearch',
  'db_table_log_evalues',
  'db_table_orthologs',
  'db_table_sequence_pairs',
  'db_table_sequence_types',
  'db_table_set_details',
  'db_table_taxa',
	'num-threads=i',
  'ortholog-set=s',
  'output-directory=s',
	'overwrite|o',
  'preparedb',
  'quiet',
  'quiet|q',
  'reference-taxa=s',
  'reference-taxon=s',
  'sets-dir=s',
  'soft-threshold=i',
  'species-name=s',
	'sqlite-database=s',
  'substitute-u-with=s',
  'verbose|v',
) or print "Fatal: I don't know what you want me to do. Terminating.\n" and exit(1);#}}}

#--------------------------------------------------
# # These variables can be set in the config file
#-------------------------------------------------- 
#{{{



# MySQL settings
$config->{'database-engine'}            //= 'mysql';
$config->{'mysql-database'}             //= 'orthograph';
$config->{'mysql-password'}             //= 'root';
$config->{'mysql-server'}               //= 'localhost';
$config->{'mysql-username'}             //= 'root';
$config->{'mysql-prefix'}               //= 'orthograph';
$config->{'mysql-timeout'}              //= 600;

# database tables
$config->{'db_table_aaseqs'}         //= 'aaseqs';
$config->{'db_table_blast'}          //= 'blast';
$config->{'db_table_blastdbs'}       //= 'blastdbs';
$config->{'db_table_ests'}           //= 'ests';
$config->{'db_table_hmmsearch'}      //= 'hmmsearch';
$config->{'db_table_log_evalues'}    //= 'log_evalues';
$config->{'db_table_scores'}         //= 'scores';
$config->{'db_table_ntseqs'}         //= 'ntseqs';
$config->{'db_table_ogs'}            //= 'ogs';
$config->{'db_table_orthologs'}      //= 'orthologs';
$config->{'db_table_sequence_pairs'} //= 'sequence_pairs';
$config->{'db_table_sequence_types'} //= 'sequence_types';
$config->{'db_table_set_details'}    //= 'set_details';
$config->{'db_table_temp'}           //= 'temp';
$config->{'db_table_taxa'}           //= 'taxa';
$config->{'db_table_users'}          //= 'users';

# make sure there is exactly one underscore at the end of the prefix
(my $db_prefix = $config->{'db-prefix'}) =~ s/_*$/_/;

# temporary hash to prepend the prefix
my %C = map { $_ => $db_prefix . $config->{$_} } grep { $_ =~ /^db_table_/ } keys %$config;

# merge the modified entries into the config hash
$config->{$_} = $C{$_} foreach keys %C;

# free memory... well, it's bound to go out of scope anyway
undef %C;

# more variables

$config->{'aaoutdir'}                   //= 'aa';
$config->{'alignment-program'}          //= 'mafft-linsi --anysymbol';
$config->{'backup'}                     //= 1;
$config->{'backup-extension'}           //= '.bak';
$config->{'blast-evalue-threshold'}     //= 10;
$config->{'blast-max-hits'}             //= 100;
$config->{'blast-program'}              //= 'blastp';
$config->{'blast-score-threshold'}      //= 10;
$config->{'blastoutdir'}                //= basename($config->{'blast-program'});
$config->{'clear-files'}                //= 0;
$config->{'clear-database'}             //= 1;
$config->{'configfile'}                 //= $configfile;
$config->{'continue'}                   //= 0;
$config->{'create'}                     //= 0;
$config->{'debug'}                      //= 0;
$config->{'delete-ogs'}                 //= '';
$config->{'delete-set'}                 //= '';
$config->{'destroy'}                    //= 0;
$config->{'estfile'}                    //= '';
$config->{'evalue-bin-size'}            //= 500;
$config->{'header-separator'}           //= '|';
$config->{'hmmbuild-program'}           //= 'hmmbuild';
$config->{'hmmsearch-evalue-threshold'} //= defined $config->{'hmmsearch-score-threshold'} ? undef : 10;
$config->{'hmmsearch-program'}          //= 'hmmsearch';
$config->{'hmmsearch-program'}          //= 'hmmsearch';
$config->{'hmmsearch-score-threshold'}  //= defined $config->{'hmmsearch-evalue-threshold'} ? undef : 10;
$config->{'hmmsearchoutdir'}            //= basename($config->{'hmmsearch-program'});
$config->{'load-ogs-nucleotide'}        //= '';
$config->{'load-ogs-peptide'}           //= '';
$config->{'logfile'}                    //= '';
$config->{'make-set'}                   //= 0;
$config->{'makeblastdb-program'}        //= 'makeblastdb';
$config->{'max-blast-searches'}         //= 100;
$config->{'ntoutdir'}                   //= 'nt';
$config->{'num-threads'}                //= 1;
$config->{'ortholog-set'}               //= '';
$config->{'output-directory'}           //= '';
$config->{'overwrite'}                  //= 0;
$config->{'prepare'}                    //= 0;  
$config->{'quiet'}                      //= 0;  # I like my quiet
$config->{'reference-taxa'}             //= '';
$config->{'sets-dir'}                   //= 'sets';
$config->{'soft-threshold'}             //= 5;
$config->{'species-name'}               //= '';
$config->{'sqlite-program'}             //= '/usr/bin/sqlite3';
# substitution character for selenocysteine, which normally leads to blast freaking out
$config->{'substitute-u-with'}          //= 'X';
$config->{'translate-program'}          //= 'fastatranslate';
$config->{'verbose'}                    //= 0;
#}}}

# compound options
if ($config->{'continue'}) {
	$config->{'clear-files'}    = 0;
	$config->{'clear-database'} = 0;
}

# mutually exclusive options
if ($config->{'database-engine'} ne 'mysql' and $config->{'database-engine'} ne 'sqlite') {
	print STDERR "Fatal: Database engine not set correctly! Must be 'mysql' or 'sqlite'.\n";
	exit 1;
}

if ($config->{'database-engine'} eq 'sqlite' and not defined $config->{'sqlite-database'}) {
	print STDERR "Fatal: SQLite database not specified\n";
	exit 1;
}


# un-quote the header separator
$config->{'header-separator'} =~ s/^('|")//;
$config->{'header-separator'} =~ s/('|")$//;

# if something went wrong
die "Fatal: Error parsing config" unless $config;

=head2 getconfig 

get the config filename from the command line and parse it.

Returns a hashref that contains all config variables from both the config file and the command line.

=cut

sub getconfig {
	# in case the user tells us to use a different one with -c
	$configfile = get_configfile($configfile);

	# parse if exists
	if (-e $configfile) {
		print "Parsing config file '$configfile'.\n";
		$config = parse_config($configfile);
	}#}}}
	else {
		print STDERR "Fatal: Config file '$configfile' not found!\n";
		exit 1;
	}

	return $config;
}

=head2 get_configfile([CONFIGFILE])

mini argument parser to get the config file name

Returns the config filename as provided with the option B<-c> on the command line.

As an optional argument, this function accepts a default config filename which
will be returned if there was no -c on the command line. Note: unless a default
is provided, this function will return undef.

=cut

# mini argument parser for the configfile
sub get_configfile {
	my $configfile = shift(@_);
	# every argument
	for (my $i = 0; $i < scalar @ARGV; ++$i) {
		# is this '-c'?
		if ($ARGV[$i] =~ /(-c|-configfile)\b/) {
			unless ($ARGV[$i+1]) { print "Config file name missing (you specified -c but didn't provide a name). Exiting.\n" and exit(1)}
			# does the next one begin with a hyphen?
			if ($ARGV[$i+1] !~ /^-/) {
				$configfile = $ARGV[$i+1];
				splice @ARGV, $i, 2;
			}
			# the file name starts with a hyphen, may be a stray option, so warn the
			# user and don't use this name
			else { warn "Warning: Config file name '$ARGV[$i+1]' starts with a hyphen (-), could be a stray option. Use './$ARGV[$i+1]' if you mean it. Falling back to '$configfile' for now.\n" }
		}
	}
	return $configfile;
}

=head2 parse_config

Parse a simple, ini-style config file where keys are separated from values by '='. 
Sections are not supported. 

Returns a hashref.

Config file example: 

  outputdir = /home/foo/bar

=cut 

sub parse_config {#{{{
	my $file = shift;
	my $conf = { };
	my $fh = IO::File->new($file) or print "Fatal: Could not open config file '$file'\: $!\n" and exit(1);

	while (my $line = $fh->getline()) {
		next if $line =~ /^\s*$/; # skip empty lines
		next if $line =~ /^\s*#/; # skip comment lines starting with '#'
		if ($line !~ /^\s*\S+\s*=\s*[\/]?\S+/) {
			print "Fatal: Invalid format in line $. of config file $file:\n$line\n" and exit(1);
		}
		
		# split by '=' producing a maximum of two items; the value may contain whitespace
		my ($key, $val) = split('=', $line, 2);

		foreach ($key, $val) {
		  s/\s+$//; # remove all trailing whitespace
		  s/^\s+//; # remove all leading whitespace
		}

		print "Fatal: Configuration option '$key' defined twice in line $. of config file '$file'\n" and exit(1)
		  if defined $conf->{$key};
		$conf->{$key} = $val;
	}
	close($fh);
	return $conf;
}#}}}


1;
