package Wrapper::Mysql;
use strict;
use warnings;
use Carp;
use FindBin;        # locate the dir of this script during compile time
use lib $FindBin::Bin;                 # $Bin is the directory of the original script
use Orthograph::Config;                # provides config in exported hashref $config

#--------------------------------------------------
# # Parse config file
#-------------------------------------------------- 
my $config = Orthograph::Config->getconfig();

#--------------------------------------------------
# # These variables can be set in the config file
#-------------------------------------------------- 
#{{{
my $debug          = $config->{'debug'}                ? $config->{'debug'}                : undef;
my $mysql_dbname   = $config->{'mysql_dbname'}         ? $config->{'mysql_dbname'}         : 'orthograph';
my $mysql_dbpwd    = $config->{'mysql_dbpassword'}     ? $config->{'mysql_dbpassword'}     : 'root';
my $mysql_dbserver = $config->{'mysql_dbserver'}       ? $config->{'mysql_dbserver'}       : 'localhost';
my $mysql_dbuser   = $config->{'mysql_dbuser'}          ? $config->{'mysql_dbuser'}         : 'root';
my $mysql_table_prefix = $config->{'mysql_table_prefix'} ? $config->{'mysql_table_prefix'} : 'orthograph';

# make sure there is exactly one underscore at the end of the prefix
$mysql_table_prefix =~ s/_*$/_/;

my $mysql_table_blastdbs = $config->{'mysql_table_blastdbs'} ?
	$mysql_table_prefix . $config->{'mysql_table_blastdbs'} :
	$mysql_table_prefix . 'blastdbs';
my $mysql_table_ests = $config->{'mysql_table_ests'} ?
	$mysql_table_prefix . $config->{'mysql_table_ests'} :
	$mysql_table_prefix . 'ests';
my $mysql_table_hmmsearch = $config->{'mysql_table_hmmsearch'} ?
	$mysql_table_prefix . $config->{'mysql_table_hmmsearch'} :
	$mysql_table_prefix . 'hmmsearch';
my $mysql_table_set_details = $config->{'mysql_table_set_details'} ?
	$mysql_table_prefix . $config->{'mysql_table_set_details'} :
	$mysql_table_prefix . 'set_details';
my $mysql_table_aaseqs = $config->{'mysql_table_aaseqs'} ?
	$mysql_table_prefix . $config->{'mysql_table_aaseqs'} :
	$mysql_table_prefix . 'aaseqs';
my $mysql_table_ntseqs = $config->{'mysql_table_ntseqs'} ?
	$mysql_table_prefix . $config->{'mysql_table_ntseqs'} :
	$mysql_table_prefix . 'ntseqs';
my $mysql_table_ogs = $config->{'mysql_table_ogs'} ?
	$mysql_table_prefix . $config->{'mysql_table_ogs'} :
	$mysql_table_prefix . 'ogs';
my $mysql_table_orthologs = $config->{'mysql_table_orthologs'} ?
	$mysql_table_prefix . $config->{'mysql_table_orthologs'} :
	$mysql_table_prefix . 'orthologs';
my $mysql_table_seqpairs       = $config->{'mysql_table_sequence_pairs'} ?
	$mysql_table_prefix . $config->{'mysql_table_sequence_pairs'} :
	$mysql_table_prefix . 'sequence_pairs';
my $mysql_table_seqtypes       = $config->{'mysql_table_sequence_types'} ?
	$mysql_table_prefix . $config->{'mysql_table_sequence_types'} :
	$mysql_table_prefix . 'sequence_types';
my $mysql_table_taxa       = $config->{'mysql_table_taxa'} ?
	$mysql_table_prefix . $config->{'mysql_table_taxa'} :
	$mysql_table_prefix . 'taxa';
my $mysql_table_temp       = $config->{'mysql_table_temp'} ?
	$mysql_table_prefix . $config->{'mysql_table_temp'} :
	$mysql_table_prefix . 'temp';
my $mysql_table_users       = $config->{'mysql_table_users'} ?
	$mysql_table_prefix . $config->{'mysql_table_users'} :
	$mysql_table_prefix . 'users';

# table names. too lazy to change all of them, so i'll just reuse the old hash structure
my %t            = (
	'aaseqs'       => $mysql_table_aaseqs,
	'blastdbs'     => $mysql_table_blastdbs,
	'ntseqs'       => $mysql_table_ntseqs,
	'ogs'          => $mysql_table_ogs,
	'orthologs'    => $mysql_table_orthologs,
	'seqpairs'     => $mysql_table_seqpairs,
	'seqtypes'     => $mysql_table_seqtypes,
	'set_details'  => $mysql_table_set_details,
	'taxa'         => $mysql_table_taxa,
	'temp'         => $mysql_table_temp,
	'users'        => $mysql_table_users,
);

# Sub: mysql_dbh
# Returns a MySQL database handle
sub mysql_dbh {#{{{
	return DBI->connect("DBI:mysql:$mysql_dbname:$mysql_dbserver", $mysql_dbuser, $mysql_dbpwd);
}#}}}

# Sub: mysql_get
# Get from the database the result of a SQL query
# Arguments: QUERY as a string literal
# Returns: Reference to array of arrays (result lines->fields)
sub mysql_get {#{{{
	my $query = shift;
	unless ($query) { croak "Usage: mysql_get(QUERY)\n" }
  # prepare anonymous array
	my $results = [ ];
  # connect and fetch stuff
	my $dbh = &mysql_dbh;
	my $sql = $dbh->prepare($query);
	$sql->execute();
	while (my @result = $sql->fetchrow_array() ) {
		push(@$results, \@result);
	}
	$sql->finish();
	$dbh->disconnect; # disconnect ASAP
	return $results;
}#}}}


sub get_taxa_in_all_sets {
	my %setlist = ();
	my $query = "SELECT DISTINCT $t{'set_details'}.name, $t{'taxa'}.name
		FROM $t{'seqpairs'}
		INNER JOIN $t{'taxa'}
			ON $t{'seqpairs'}.taxid = $t{'taxa'}.id
		INNER JOIN $t{'orthologs'}
			ON $t{'orthologs'}.sequence_pair = $t{'seqpairs'}.id 
		INNER JOIN $t{'set_details'}
			ON $t{'orthologs'}.setid = $t{'set_details'}.id"
	;
	my $data = &mysql_get($query);
	foreach my $row (@$data) {
		$setlist{$$row[0]} .= ' ' . $$row[1];
	}
	return %setlist;
}

sub get_taxa_in_set {
	my $setname = shift @_;
	unless ($setname) { croak("Usage: get_taxa_in_set(SETNAME)") }
	my @reftaxa;
	my $query = "SELECT DISTINCT $t{'set_details'}.name, $t{'taxa'}.name
		FROM $t{'seqpairs'}
		INNER JOIN $t{'taxa'}
			ON $t{'seqpairs'}.taxid = $t{'taxa'}.id
		INNER JOIN $t{'orthologs'}
			ON $t{'orthologs'}.sequence_pair = $t{'seqpairs'}.id 
		INNER JOIN $t{'set_details'}
			ON $t{'orthologs'}.setid = $t{'set_details'}.id
		WHERE $t{'set_details'} = '$setname'"
	;
	my $data = &mysql_get($query);
	foreach my $row (@$data) {
		push(@reftaxa, $$row[1]);
	}
	return @reftaxa;
}


1;
