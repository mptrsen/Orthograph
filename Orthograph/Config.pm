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
our $config = &getconfig; 

#--------------------------------------------------
# # Get command line options. These may override variables set in the config file.
#-------------------------------------------------- 
#{{{
GetOptions( $config,
  'aaoutdir',
  'alignment_program',
  'backup_extension',
  'blast_evalue_threshold',
  'blast_max_hits',
  'blast_score_threshold',
  'blastoutdir',
  'blastp_output_dir',
  'clear_data',
  'debug',
  'estfile',
  'hmmfile',
  'hmmsearch_evalue_threshold',
  'hmmsearch_output_dir',
  'hmmsearch_score_threshold',
	'list-ests|le',
	'list-ogs|lo',
	'list-sets|ls',
	'list-taxa|lt',
  'logfile',
  'max_blast_searches',
  'mysql_dbname',
  'mysql_dbpassword',
  'mysql_dbserver',
  'mysql_dbuser',
  'mysql_prefix',
  'mysql_table_aaseqs',
  'mysql_table_blast',
  'mysql_table_blastdbs',
  'mysql_table_ests',
  'mysql_table_hmmsearch',
  'mysql_table_orthologs',
  'mysql_table_sequence_pairs',
  'mysql_table_sequence_types',
  'mysql_table_set_details',
  'mysql_table_taxa',
  'ortholog_set',
  'output_directory',
  'quiet',
  'reference_taxa=s',
  'reference_taxon=s',
  'sets_dir=s',
  'soft_threshold=i',
  'substitute_u_with=s',
  'blast_evalue_threshold=f',
  'blast_max_hits=i',
  'c=s',
  'debug|d',
  'estfile|E=s',
  'hmmsearch_evalue_threshold=f',
  'hmmsearch_score=i',
  'hmmsearchprog=s',
  'list_species|l',
  'logfile|log=s',
  'preparedb',
  'quiet',
  'species_name=s',
  'verbose|v',
) or die("Fatal: I don't know what you want me to do. Terminating.\n");#}}}

#--------------------------------------------------
# # These variables can be set in the config file
#-------------------------------------------------- 
#{{{



# MySQL settings
defined $config->{'mysql_database'}             or $config->{'mysql_database'}               = 'orthograph';
defined $config->{'mysql_password'}             or $config->{'mysql_password'}               = 'root';
defined $config->{'mysql_server'}               or $config->{'mysql_server'}                 = 'localhost';
defined $config->{'mysql_username'}             or $config->{'mysql_username'}               = 'root';
defined $config->{'mysql_prefix'}               or $config->{'mysql_prefix'}                 = 'orthograph';

# MySQL tables
defined $config->{'mysql_table_aaseqs'}         or $config->{'mysql_table_aaseqs'}         = 'aaseqs';
defined $config->{'mysql_table_blast'}          or $config->{'mysql_table_blast'}          = 'blast';
defined $config->{'mysql_table_blastdbs'}       or $config->{'mysql_table_blastdbs'}       = 'blastdbs';
defined $config->{'mysql_table_ests'}           or $config->{'mysql_table_ests'}           = 'ests';
defined $config->{'mysql_table_hmmsearch'}      or $config->{'mysql_table_hmmsearch'}      = 'hmmsearch';
defined $config->{'mysql_table_ntseqs'}         or $config->{'mysql_table_ntseqs'}         = 'ntseqs';
defined $config->{'mysql_table_ogs'}            or $config->{'mysql_table_ogs'}            = 'ogs';
defined $config->{'mysql_table_orthologs'}      or $config->{'mysql_table_orthologs'}      = 'orthologs';
defined $config->{'mysql_table_sequence_pairs'} or $config->{'mysql_table_sequence_pairs'} = 'sequence_pairs';
defined $config->{'mysql_table_sequence_types'} or $config->{'mysql_table_sequence_types'} = 'sequence_types';
defined $config->{'mysql_table_set_details'}    or $config->{'mysql_table_set_details'}    = 'set_details';
defined $config->{'mysql_table_taxa'}           or $config->{'mysql_table_taxa'}           = 'taxa';
defined $config->{'mysql_table_users'}          or $config->{'mysql_table_users'}          = 'users';

# make sure there is exactly one underscore at the end of the prefix
(my $mysql_prefix = $config->{'mysql_prefix'}) =~ s/_*$/_/;

# temporary hash to prepend the prefix
my %C = map { $_ => $mysql_prefix . $config->{$_} } grep { $_ =~ /^mysql_table_/ } keys %$config;

# merge the modified entries into the config hash
$config->{$_} = $C{$_} foreach keys %C;

# free memory... well, it's bound to go out of scope anyway
undef %C;

# more variables

defined $config->{'aaoutdir'}                   or $config->{'aaoutdir'}                   = 'aa';
defined $config->{'alignment_program'}          or $config->{'alignment_program'}          = 'alignment';
defined $config->{'backup_extension'}           or $config->{'backup_extension'}           = '.bak';
defined $config->{'blast_evalue_threshold'}     or $config->{'blast_evalue_threshold'}     = 10;
defined $config->{'blast_max_hits'}             or $config->{'blast_max_hits'}             = 10;
defined $config->{'blast_program'}              or $config->{'blast_program'}              = 'blastp';
defined $config->{'blast_score_threshold'}      or $config->{'blast_score_threshold'}      = 10;
defined $config->{'blastoutdir'}                or $config->{'blastoutdir'}                = basename($config->{'blast_program'});
defined $config->{'clear_data'}                 or $config->{'clear_data'}                 = 1;
defined $config->{'debug'}                      or $config->{'debug'}                      = 0;
defined $config->{'estfile'}                    or $config->{'estfile'}                    = '';
defined $config->{'hmmbuild_program'}           or $config->{'hmmbuild_program'}           = 'hmmbuild';
defined $config->{'hmmfile'}                    or $config->{'hmmfile'}                    = '';
defined $config->{'hmmsearch_evalue_threshold'} or $config->{'hmmsearch_evalue_threshold'} = undef;
defined $config->{'hmmsearch_program'}          or $config->{'hmmsearch_program'}          = 'hmmsearch';
defined $config->{'hmmsearchoutdir'}            or $config->{'hmmsearchoutdir'}            = basename($config->{'hmmsearch_program'});
defined $config->{'makeblastdb_program'}        or $config->{'makeblastdb_program'}        = 'makeblastdb';
defined $config->{'translate_program'}          or $config->{'translate_program'}          = 'fastatranslate';
defined $config->{'hmmsearch_score_threshold'}  or $config->{'hmmsearch_score_threshold'}  = $config->{'hmmsearch_evalue_threshold'} ? undef : 10;
defined $config->{'logfile'}                    or $config->{'logfile'}                    = '';
defined $config->{'max_blast_searches'}         or $config->{'max_blast_searches'}         = 100;
defined $config->{'ortholog_set'}               or $config->{'ortholog_set'}               = '';
defined $config->{'output_directory'}           or $config->{'output_directory'}           = '';
defined $config->{'quiet'}                      or $config->{'quiet'}                      = 0;  # I like my quiet
defined $config->{'reference_taxa'}             or $config->{'reference_taxa'}             = '';
defined $config->{'reference_taxon'}            or $config->{'reference_taxon'}            = '';
defined $config->{'sets_dir'}                   or $config->{'sets_dir'}                   = 'sets';
defined $config->{'species_name'}               or $config->{'species_name'}               = '';
defined $config->{'soft_threshold'}             or $config->{'soft_threshold'}             = 5;
# substitution character for selenocysteine, which normally leads to blast freaking out
defined $config->{'substitute_u_with'}          or $config->{'substitute_u_with'}          = 'X';
defined $config->{'verbose'}                    or $config->{'verbose'}                    = 0;
#}}}

# if something went wrong
die unless $config;

=head2 getconfig 

get the config filename from the command line and parse it.

Returns a hashref that contains all config variables from both the config file and the command line.

=cut

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
		if ($ARGV[$i] =~ /-c\b/) {
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
	my $fh = IO::File->new($file) or die "Fatal: Could not open config file '$file'\: $!\n";

	while (my $line = $fh->getline()) {
		next if $line =~ /^\s*$/; # skip empty lines
		next if $line =~ /^\s*#/; # skip comment lines starting with '#'
		if ($line !~ /^\s*\w+\s*=\s*[\/]?\w+/) {
			die "Fatal: Invalid format in line $. of config file $file:\n$line\n"
		}
		
		# split by '=' producing a maximum of two items; the value may contain whitespace
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
