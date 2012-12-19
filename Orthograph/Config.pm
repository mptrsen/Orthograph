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
	'delete-ogs',
	'delete-set',
	'destroy',
	'evalue-bin-size=i',
	'list-ests|le',
	'list-ogs|lo',
	'list-sets|ls',
	'list-taxa|lt',
	'load-ogs-nucleotide=s',
	'load-ogs-peptide=s',
	'prepare',
  'aaoutdir',
  'alignment-program',
	'backup',
  'backup-extension',
  'blast-evalue-threshold',
  'blast-evalue-threshold=f',
  'blast-max-hits=i',
  'blast-score-threshold',
  'blastoutdir',
  'blastp-output-dir',
  'clear-database!',
  'clear-files!',
  'configfile|c=s',
  'debug',
  'debug|d',
  'estfile',
  'estfile|E=s',
  'hmmsearch-evalue-threshold',
  'hmmsearch-evalue-threshold=f',
  'hmmsearch-output-dir',
  'hmmsearch-score=i',
  'hmmsearch-score-threshold',
  'hmmsearchprog=s',
  'logfile',
  'logfile|log=s',
  'max-blast-searches',
  'mysql-database',
  'mysql-password',
  'mysql-prefix',
  'mysql-server',
  'mysql_table_aaseqs',
  'mysql_table_blast',
  'mysql_table_blastdbs',
  'mysql_table_ests',
  'mysql_table_hmmsearch',
  'mysql_table_log_evalues',
  'mysql_table_orthologs',
  'mysql_table_sequence_pairs',
  'mysql_table_sequence_types',
  'mysql_table_set_details',
  'mysql_table_taxa',
  'mysql-username',
  'ortholog-set',
  'output-directory',
  'preparedb',
  'quiet',
  'quiet|q',
  'reference-taxa=s',
  'reference-taxon=s',
  'sets-dir=s',
  'soft-threshold=i',
  'species-name=s',
  'substitute-u-with=s',
  'verbose|v',
) or print "Fatal: I don't know what you want me to do. Terminating.\n" and exit(1);#}}}

#--------------------------------------------------
# # These variables can be set in the config file
#-------------------------------------------------- 
#{{{



# MySQL settings
defined $config->{'mysql-database'}              or $config->{'mysql-database'}               = 'orthograph';
defined $config->{'mysql-password'}              or $config->{'mysql-password'}               = 'root';
defined $config->{'mysql-server'}                or $config->{'mysql-server'}                 = 'localhost';
defined $config->{'mysql-username'}              or $config->{'mysql-username'}               = 'root';
defined $config->{'mysql-prefix'}                or $config->{'mysql-prefix'}                 = 'orthograph';

# MySQL tables
defined $config->{'mysql_table_aaseqs'}          or $config->{'mysql_table_aaseqs'}         = 'aaseqs';
defined $config->{'mysql_table_blast'}           or $config->{'mysql_table_blast'}          = 'blast';
defined $config->{'mysql_table_blastdbs'}        or $config->{'mysql_table_blastdbs'}       = 'blastdbs';
defined $config->{'mysql_table_ests'}            or $config->{'mysql_table_ests'}           = 'ests';
defined $config->{'mysql_table_hmmsearch'}       or $config->{'mysql_table_hmmsearch'}      = 'hmmsearch';
defined $config->{'mysql_table_log_evalues'}     or $config->{'mysql_table_log_evalues'}    = 'log_evalues';
defined $config->{'mysql_table_ntseqs'}          or $config->{'mysql_table_ntseqs'}         = 'ntseqs';
defined $config->{'mysql_table_ogs'}             or $config->{'mysql_table_ogs'}            = 'ogs';
defined $config->{'mysql_table_orthologs'}       or $config->{'mysql_table_orthologs'}      = 'orthologs';
defined $config->{'mysql_table_sequence_pairs'}  or $config->{'mysql_table_sequence_pairs'} = 'sequence_pairs';
defined $config->{'mysql_table_sequence_types'}  or $config->{'mysql_table_sequence_types'} = 'sequence_types';
defined $config->{'mysql_table_set_details'}     or $config->{'mysql_table_set_details'}    = 'set_details';
defined $config->{'mysql_table_temp'}            or $config->{'mysql_table_temp'}           = 'temp';
defined $config->{'mysql_table_taxa'}            or $config->{'mysql_table_taxa'}           = 'taxa';
defined $config->{'mysql_table_users'}           or $config->{'mysql_table_users'}          = 'users';

# make sure there is exactly one underscore at the end of the prefix
(my $mysql_prefix = $config->{'mysql-prefix'}) =~ s/_*$/_/;

# temporary hash to prepend the prefix
my %C = map { $_ => $mysql_prefix . $config->{$_} } grep { $_ =~ /^mysql_table_/ } keys %$config;

# merge the modified entries into the config hash
$config->{$_} = $C{$_} foreach keys %C;

# free memory... well, it's bound to go out of scope anyway
undef %C;

# more variables

defined $config->{'aaoutdir'}                    or $config->{'aaoutdir'}                   = 'aa';
defined $config->{'alignment-program'}           or $config->{'alignment-program'}          = 'mafft-linsi --anysymbol';
defined $config->{'backup'}                      or $config->{'backup'}                     = 1;
defined $config->{'backup-extension'}            or $config->{'backup-extension'}           = '.bak';
defined $config->{'blast-evalue-threshold'}      or $config->{'blast-evalue-threshold'}     = 10;
defined $config->{'blast-max-hits'}              or $config->{'blast-max-hits'}             = 100;
defined $config->{'blast-program'}               or $config->{'blast-program'}              = 'blastp';
defined $config->{'blast-score-threshold'}       or $config->{'blast-score-threshold'}      = 10;
defined $config->{'blastoutdir'}                 or $config->{'blastoutdir'}                = basename($config->{'blast-program'});
defined $config->{'clear-files'}                 or $config->{'clear-files'}                = 0;
defined $config->{'clear-database'}              or $config->{'clear-database'}             = 1;
defined $config->{'configfile'}                  or $config->{'configfile'}                 = $configfile;
defined $config->{'continue'}                    or $config->{'continue'}                   = 0;
defined $config->{'create'}                      or $config->{'create'}                     = 0;
defined $config->{'debug'}                       or $config->{'debug'}                      = 0;
defined $config->{'delete-ogs'}                  or $config->{'delete-ogs'}                 = '';
defined $config->{'delete-set'}                  or $config->{'delete-set'}                 = '';
defined $config->{'destroy'}                     or $config->{'destroy'}                    = 0;
defined $config->{'estfile'}                     or $config->{'estfile'}                    = '';
defined $config->{'evalue-bin-size'}             or $config->{'evalue-bin-size'}            = 500;
defined $config->{'hmmbuild-program'}            or $config->{'hmmbuild-program'}           = 'hmmbuild';
defined $config->{'hmmsearch-evalue-threshold'}  or $config->{'hmmsearch-evalue-threshold'} = defined $config->{'hmmsearch-score-threshold'} ? undef : 10;
defined $config->{'hmmsearch-program'}           or $config->{'hmmsearch-program'}          = 'hmmsearch';
defined $config->{'hmmsearch-program'}           or $config->{'hmmsearch-program'}          = 'hmmsearch';
defined $config->{'hmmsearch-score-threshold'}   or $config->{'hmmsearch-score-threshold'}  = defined $config->{'hmmsearch-evalue-threshold'} ? undef : 10;
defined $config->{'hmmsearchoutdir'}             or $config->{'hmmsearchoutdir'}            = basename($config->{'hmmsearch-program'});
defined $config->{'load-ogs-nucleotide'}         or $config->{'load-ogs-nucleotide'}        = '';
defined $config->{'load-ogs-peptide'}            or $config->{'load-ogs-peptide'}           = '';
defined $config->{'logfile'}                     or $config->{'logfile'}                    = '';
defined $config->{'makeblastdb-program'}         or $config->{'makeblastdb-program'}        = 'makeblastdb';
defined $config->{'max-blast-searches'}          or $config->{'max-blast-searches'}         = 100;
defined $config->{'ortholog-set'}                or $config->{'ortholog-set'}               = '';
defined $config->{'output-directory'}            or $config->{'output-directory'}           = '';
defined $config->{'prepare'}                     or $config->{'prepare'}                    = 0;  
defined $config->{'quiet'}                       or $config->{'quiet'}                      = 0;  # I like my quiet
defined $config->{'reference-taxa'}              or $config->{'reference-taxa'}             = '';
defined $config->{'sets-dir'}                    or $config->{'sets-dir'}                   = 'sets';
defined $config->{'soft-threshold'}              or $config->{'soft-threshold'}             = 5;
defined $config->{'species-name'}                or $config->{'species-name'}               = '';
# substitution character for selenocysteine, which normally leads to blast freaking out
defined $config->{'substitute-u-with'}           or $config->{'substitute-u-with'}          = 'X';
defined $config->{'translate-program'}           or $config->{'translate-program'}          = 'fastatranslate';
defined $config->{'verbose'}                     or $config->{'verbose'}                    = 0;
#}}}

# compound options
if ($config->{'continue'}) {
	$config->{'clear-files'}    = 0;
	$config->{'clear-database'} = 0;
}

# if something went wrong
die unless $config;

=head2 getconfig 

get the config filename from the command line and parse it.

Returns a hashref that contains all config variables from both the config file and the command line.

=cut

sub getconfig {
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

=head2 get_configfile([CONFIGFILE])

mini argument parser to get the config file name

Returns the config filename as provided on the command line.

As an optional argument, this function accepts a default config filename which
will be returned if there was no -c on the command line. Note: unless a default
is provided, this function will return B<undef>.

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
