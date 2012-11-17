package Wrapper::Mysql;
use strict;
use warnings;
use Carp;
use Exporter;
use FindBin;        # locate the dir of this script during compile time
use lib $FindBin::Bin;                 # $Bin is the directory of the original script
use Orthograph::Config;                # configuration parser getconfig()
use Data::Dumper;

my $config = $Orthograph::Config::config;  # copy config
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

my $mysql_table_blast = $config->{'mysql_table_blast'} ?
	$mysql_table_prefix . $config->{'mysql_table_blast'} :
	$mysql_table_prefix . 'blast';
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
	$sql->execute() or die;
	while (my @result = $sql->fetchrow_array() ) {
		push(@$results, \@result);
	}
	$sql->finish();
	$dbh->disconnect; # disconnect ASAP
	return $results;
}#}}}


=head2 get_taxa_in_all_sets

Get a list of sets associated with the included taxa.

Arguments: None

Returns: hash of scalars

=cut
sub get_taxa_in_all_sets {
	my %setlist = ();
	my $query = "SELECT DISTINCT $mysql_table_set_details.name, $mysql_table_taxa.name
		FROM $mysql_table_seqpairs
		INNER JOIN $mysql_table_taxa
			ON $mysql_table_seqpairs.taxid = $mysql_table_taxa.id
		INNER JOIN $mysql_table_orthologs
			ON $mysql_table_orthologs.sequence_pair = $mysql_table_seqpairs.id 
		INNER JOIN $mysql_table_set_details
			ON $mysql_table_orthologs.setid = $mysql_table_set_details.id"
	;
	my $data = &mysql_get($query) or croak();
	foreach my $row (@$data) {
		$setlist{$$row[0]} .= ' ' . $$row[1];
	}
	return %setlist;
}

=head2 get_taxa_in_set(SETNAME)

Returns a list of taxon names for a named set.

Arguments: scalar string SETNAME

Returns: array of scalars

=cut
sub get_taxa_in_set {
	my $setname = shift @_;
	unless ($setname) { croak("Usage: get_taxa_in_set(SETNAME)") }
	my @reftaxa;
	my $query = "SELECT DISTINCT $mysql_table_set_details.name, $mysql_table_taxa.name
		FROM $mysql_table_seqpairs
		INNER JOIN $mysql_table_taxa
			ON $mysql_table_seqpairs.taxid = $mysql_table_taxa.id
		INNER JOIN $mysql_table_orthologs
			ON $mysql_table_orthologs.sequence_pair = $mysql_table_seqpairs.id 
		INNER JOIN $mysql_table_set_details
			ON $mysql_table_orthologs.setid = $mysql_table_set_details.id
		WHERE $mysql_table_set_details = '$setname'"
	;
	my $data = &mysql_get($query);
	foreach my $row (@$data) {
		push(@reftaxa, $$row[1]);
	}
	return @reftaxa;
}

=head2 get_taxid_for_species(SPECIESNAME)

Returns the taxid for a named species.

Arguments: scalar string SPECIESNAME

Returns: scalar int TAXID

=cut
sub get_taxid_for_species {
	my $species_name = shift(@_);
	unless ($species_name) { croak("Usage: get_taxid_for_species(SPECIESNAME)") }
	my $query = "SELECT id FROM $mysql_table_taxa WHERE core = 0 AND longname = '$species_name'";
	my $result = &mysql_get($query);
	if ($result) { return $$result[0][0] }
	return 0;
}

=head2 get_set_id(SETNAME)

get the set id for a named set. 

Arguments: scalar string SETNAME

Returns: scalar int SETID

=cut
sub get_set_id {
	my $setname = shift(@_);
	unless ($setname) { croak("Usage: get_set_id(SETNAME)") }
	my $query = "SELECT id FROM $mysql_table_set_details WHERE name = '$setname'";
	my $result = &mysql_get($query);
	if ( scalar(@$result) > 1 ) { 
		warn("Warning: Multiple sets of the same name!\n");
		return $$result[0][0];
	}
	return $$result[0][0];
}

# get a orthoid => list_of_aaseq_ids relationship from the db
sub get_orthologs_for_set_hashref {
	my $setid = shift(@_);
	unless ($setid) { croak("Usage: get_orthologs_for_set(SETID)") }
	my $query = "SELECT $mysql_table_orthologs.ortholog_gene_id, $mysql_table_aaseqs.id 
		FROM $mysql_table_orthologs 
		INNER JOIN $mysql_table_seqpairs 
			ON $mysql_table_orthologs.sequence_pair = $mysql_table_seqpairs.id
		INNER JOIN $mysql_table_aaseqs
			ON $mysql_table_seqpairs.aa_seq = $mysql_table_aaseqs.id
		INNER JOIN $mysql_table_set_details 
			ON $mysql_table_orthologs.setid = $mysql_table_set_details.id
		WHERE $mysql_table_set_details.id = $setid";
	my $data = &mysql_get($query);
	my $result = { };
	foreach my $line (@$data) {
		push( @{$$result{$$line[0]}}, $$line[1] );
	}
	return $result;
}

=head2 get_hitlist_hashref(SPECIESID, SETID)


Get the results in the form:

  evalue => {
    orthoid => [
      reciprocal_hit,
      reciprocal_hit,
      reciprocal_hit,
      etc.
    ]
		orthoid2 => [
      reciprocal_hit,
      reciprocal_hit,
		]
  }

Arguments: scalar int SPECIESID, scalar int SETID

Returns: hashref of hashrefs of arrays of hashes - lol

=cut
sub get_hitlist_hashref {
	my $specid = shift(@_) or croak("Usage: get_hitlist_for(SPECIESID, SETID)");
	my $setid  = shift(@_) or croak("Usage: get_hitlist_for(SPECIESID, SETID)");
	my $query = "SELECT 
		$mysql_table_hmmsearch.evalue,
		$mysql_table_orthologs.ortholog_gene_id, 
		$mysql_table_hmmsearch.target,
		$mysql_table_blast.target,
		$mysql_table_blast.evalue,
		$mysql_table_taxa.name
		FROM $mysql_table_hmmsearch
		INNER JOIN $mysql_table_orthologs
			ON $mysql_table_hmmsearch.query = $mysql_table_orthologs.ortholog_gene_id
		INNER JOIN $mysql_table_blast
			ON $mysql_table_hmmsearch.target = $mysql_table_blast.query
		INNER JOIN $mysql_table_aaseqs
			ON $mysql_table_blast.target = $mysql_table_aaseqs.id
		INNER JOIN  $mysql_table_taxa
			ON $mysql_table_aaseqs.taxid = $mysql_table_taxa.id
		INNER JOIN $mysql_table_set_details
			ON $mysql_table_orthologs.setid = $mysql_table_set_details.id
		ORDER BY $mysql_table_hmmsearch.evalue DESC
		";
	my $data = &mysql_get($query);
	my $result = { };
	foreach my $line ( @$data ) {
		$$result{$$line[0]} = $$line[1];
	}
	return $result;
}
	
1;
