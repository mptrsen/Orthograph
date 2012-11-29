package Orthograph::Config;
use strict;
use warnings;

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
  'mysql_table_set_details',
  'mysql_table_taxa',
  'ortholog_set',
  'output_directory',
  'quiet',
  'reference_taxon=s@',
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
defined $config->{'mysql_dbname'}               or $config->{'mysql_dbname'}               = 'orthograph';
defined $config->{'mysql_dbpassword'}           or $config->{'mysql_dbpassword'}           = 'root';
defined $config->{'mysql_dbserver'}             or $config->{'mysql_dbserver'}             = 'localhost';
defined $config->{'mysql_dbuser'}               or $config->{'mysql_dbuser'}               = 'root';
defined $config->{'mysql_prefix'}               or $config->{'mysql_prefix'}               = 'orthograph';

# MySQL tables
defined $config->{'mysql_table_aaseqs'}         or $config->{'mysql_table_aaseqs'}         = 'aaseqs';
defined $config->{'mysql_table_blast'}          or $config->{'mysql_table_blast'}          = 'blast';
defined $config->{'mysql_table_blastdbs'}       or $config->{'mysql_table_blastdbs'}       = 'blastdbs';
defined $config->{'mysql_table_ests'}           or $config->{'mysql_table_ests'}           = 'ests';
defined $config->{'mysql_table_hmmsearch'}      or $config->{'mysql_table_hmmsearch'}      = 'hmmsearch';
defined $config->{'mysql_table_orthologs'}      or $config->{'mysql_table_orthologs'}      = 'orthologs';
defined $config->{'mysql_table_sequence_pairs'} or $config->{'mysql_table_sequence_pairs'} = 'sequence_pairs';
defined $config->{'mysql_table_set_details'}    or $config->{'mysql_table_set_details'}    = 'set_details';
defined $config->{'mysql_table_taxa'}           or $config->{'mysql_table_taxa'}           = 'taxa';

# more variables
defined $config->{'aaoutdir'}                   or $config->{'aaoutdir'}                   = 'aa';
defined $config->{'backup_extension'}           or $config->{'backup_extension'}           = '.bak';
defined $config->{'blast_evalue_threshold'}     or $config->{'blast_evalue_threshold'}     = 10;
defined $config->{'blast_max_hits'}             or $config->{'blast_max_hits'}             = 10;
defined $config->{'blast_score_threshold'}      or $config->{'blast_score_threshold'}      = 10;
defined $config->{'blastoutdir'}                or $config->{'blastoutdir'}                = 'blastp';
defined $config->{'clear_data'}                 or $config->{'clear_data'}                 = 1;
defined $config->{'debug'}                      or $config->{'debug'}                      = 0;
defined $config->{'estfile'}                    or $config->{'estfile'}                    = '';
defined $config->{'hmmfile'}                    or $config->{'hmmfile'}                    = '';
defined $config->{'hmmsearch_evalue_threshold'} or $config->{'hmmsearch_evalue_threshold'} = undef;
defined $config->{'hmmsearch_output_dir'}       or $config->{'hmmsearch_output_dir'}       = basename($config->{'hmmsearch_program'});
# mutually exclusive options
defined $config->{'hmmsearch_score_threshold'}  or $config->{'hmmsearch_score_threshold'}  = $config->{'hmmsearch_evalue_threshold'} ? undef : 10;
defined $config->{'logfile'}                    or $config->{'logfile'}                    = '';
defined $config->{'max_blast_searches'}         or $config->{'max_blast_searches'}         = 20;

defined $config->{'ortholog_set'}               or $config->{'ortholog_set'}               = '';
defined $config->{'output_directory'}           or $config->{'output_directory'}           = '';
defined $config->{'quiet'}                      or $config->{'quiet'}                      = 0;  # I like my quiet
defined $config->{'reference_taxa'}             or $config->{'reference_taxa'}             = [ ];
defined $config->{'reference_taxon'}            or $config->{'reference_taxon'}            = [ ];
defined $config->{'sets_dir'}                   or $config->{'sets_dir'}                   = 'sets';
defined $config->{'species_name'}               or $config->{'species_name'}               = '';
# substitution character for selenocysteine, which normally leads to blast freaking out
defined $config->{'substitute_u_with'}          or $config->{'substitute_u_with'}          = 'X';
defined $config->{'verbose'}                    or $config->{'verbose'}                    = 0;
#}}}

# make sure there is exactly one underscore at the end of the prefix
(my $mysql_prefix = $config->{'mysql_prefix'}) =~ s/_*$/_/;

# temporary hash to prepend the prefix
my %C = map { $_ => $mysql_prefix . $config->{$_} } grep { $_ =~ /^mysql_table_/ } keys %$config;

# merge the modified entries into the config hash
$config->{$_} = $C{$_} foreach keys %C;

# free memory... well, it's bound to go out of scope anyway
undef %C;

# special option: reference_taxa is a LIST
# if provided in the config file, it's a string that must be transformed into an arrayref
if (defined $config->{'reference_taxa'}) {
	$config->{'reference_taxa'} = [ split(/\s*,\s*/, $config->{'reference_taxa'}) ] 
}
# if provided on the command line as multiple --reference_taxon, it is an arrayref already
if (defined $config->{'reference_taxon'}) {
	push(@{$config->{'reference_taxa'}}, @{$config->{'reference_taxon'}});
}

printf("%26s => %s\n", $_, $config->{$_}) foreach sort keys %$config;
printf "<%s>", $_ foreach @{$config->{'reference_taxa'}};
exit;
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
