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
	return DBI->connect("DBI:mysql:$mysql_dbname:$mysql_dbserver;mysql_local_infile=1", $mysql_dbuser, $mysql_dbpwd);
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
	print "$query;\n" 
		if $debug;
	my $dbh = &mysql_dbh;
	my $sql = $dbh->prepare($query);
	$sql->execute() or croak;
	while (my @result = $sql->fetchrow_array() ) {
		push(@$results, \@result);
	}
	$sql->finish();
	$dbh->disconnect; # disconnect ASAP
	return $results;
}#}}}

# Sub: get_ortholog_sets
# Get list of ortholog sets from the database
# Arguments: none
# Returns: hash reference of set names => description
sub get_ortholog_sets {#{{{
	my %sets = ();
	my $query = "SELECT * FROM $mysql_table_set_details";
	my $data = &Wrapper::Mysql::mysql_get($query);
	foreach my $item (@$data) {
		$sets{$$item[1]} = $$item[2];
	}
	return(\%sets);
}#}}}

#TODO merge with get_ortholog_sets() into one function that accepts a query
# Sub: list_ogs
# Get list of OGS in the database
# Arguments: none
# Returns: array reference (list of OGS)
sub get_list_of_ogs {#{{{
	my %ogslist = ();
	# TODO rewrite this part using parametrized queries to protect from SQL injections?
	my $query = "SELECT DISTINCT $mysql_table_taxa.name , $mysql_table_ogs.version
		FROM $mysql_table_aaseqs
		INNER JOIN $mysql_table_seqpairs
			ON $mysql_table_aaseqs.id  = $mysql_table_seqpairs.aa_seq
		INNER JOIN $mysql_table_taxa
			ON $mysql_table_seqpairs.taxid = $mysql_table_taxa.id
		INNER JOIN $mysql_table_ogs
			ON $mysql_table_taxa.id = $mysql_table_ogs.taxid"
	;
	my $data = &Wrapper::Mysql::mysql_get($query);
	foreach my $item (@$data) {
		$ogslist{$$item[0]} = $$item[1];
	}
	return(\%ogslist);
}#}}}


=head2 get_taxa_in_all_sets

Get a list of sets associated with the included taxa.

Arguments: None

Returns: hash of scalars

=cut
sub get_taxa_in_all_sets {
	my %setlist = ();
	# TODO rewrite this part using parametrized queries to protect from SQL injections?
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
	my $set_id = shift @_;
	unless ($set_id) { croak("Usage: get_taxa_in_set(SETNAME)") }
	my @reftaxa;
	# TODO rewrite this part using parametrized queries to protect from SQL injections?
	my $query = "SELECT DISTINCT $mysql_table_set_details.name, $mysql_table_taxa.name
		FROM $mysql_table_seqpairs
		INNER JOIN $mysql_table_taxa
			ON $mysql_table_seqpairs.taxid = $mysql_table_taxa.id
		INNER JOIN $mysql_table_orthologs
			ON $mysql_table_orthologs.sequence_pair = $mysql_table_seqpairs.id 
		INNER JOIN $mysql_table_set_details
			ON $mysql_table_orthologs.setid = $mysql_table_set_details.id
		WHERE $mysql_table_set_details.id = '$set_id'"
	;
	my $data = &mysql_get($query);
	foreach my $row (@$data) {
		push(@reftaxa, $$row[1]);
	}
	return @reftaxa;
}

sub get_number_of_ests_for_specid {
	my $specid = shift @_ or croak "Usage: get_number_of_ests_for_specid(SPECID)";

	# TODO rewrite this part using parametrized queries to protect from SQL injections?
	my $result = &mysql_get("SELECT COUNT(*) FROM $mysql_table_ests WHERE taxid = '$specid'");

	return $$result[0][0];
}

sub get_taxids_in_set {
	my $setid = shift @_ or croak "Usage: get_taxids_in_set(SETID)";

	# get list of taxon names for this set
	my @taxa = &get_taxa_in_set($setid);

	# make a string for the next query
	my $taxa_string = join ',', map { $_ = "'$_'" } @taxa;

	# get the taxids for them
	# TODO rewrite this part using parametrized queries to protect from SQL injections?
	my $taxids = &mysql_get("SELECT id FROM $mysql_table_taxa WHERE name IN ($taxa_string)");

	return $taxids;
}

sub get_number_of_ests_for_set {
	my $setid = shift @_ or croak "Usage: get_taxids_in_set(SETID)";

	# get list of taxids for this set
	my $taxids = &get_taxids_in_set($setid);

	# make a fucking string out of these fucking fucks
	my $taxids_string = join ',', map { $$_[0] = "'$$_[0]'" } @$taxids;

	# get the number of aaseqs for those taxids
	# TODO rewrite this part using parametrized queries to protect from SQL injections?
	my $aaseqs = &mysql_get("SELECT COUNT(*) FROM  $mysql_table_aaseqs WHERE $mysql_table_aaseqs.taxid IN ($taxids_string)");

	return $$aaseqs[0][0];
}

=head2 get_aaseqs_for_set(SETID)

Returns a hashref for all aa sequences of taxa in a set, for creation of a BLAST database.
Don't forget to undef this ref (or let it go out of scope), as it is potentially very large!

Arguments: scalar int TAXID

Returns: hashref { ID => SEQ }

=cut
sub get_aaseqs_for_set {
	my $setid = shift @_ or croak "Usage: get_aaseqs_for_set(SETID)";

	# get the taxids for them
	my $taxids = &get_taxids_in_set($setid);

	# make a fucking string out of these fucking fucks
	my $taxids_string = join ',', map { $$_[0] = "'$$_[0]'" } @$taxids;

	# get the aaseqs for those taxids
	# this is a potentially very large collection, i hope that's fine with you
	# TODO rewrite this part using parametrized queries to protect from SQL injections?
	my $aaseqs = &mysql_get("SELECT $mysql_table_aaseqs.id, $mysql_table_aaseqs.sequence FROM  $mysql_table_aaseqs WHERE $mysql_table_aaseqs.taxid IN ($taxids_string)");

	$aaseqs = { map { $$_[0] => $$_[1] } @$aaseqs };
	return $aaseqs;
}


=head2 get_taxid_for_species(SPECIESNAME)

Returns the taxid for a named species.

Arguments: scalar string SPECIESNAME

Returns: scalar int TAXID

=cut
sub get_taxid_for_species {
	my $species_name = shift(@_);
	unless ($species_name) { croak("Usage: get_taxid_for_species(SPECIESNAME)") }
	# TODO rewrite this part using parametrized queries to protect from SQL injections?
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
	# TODO rewrite this part using parametrized queries to protect from SQL injections?
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
	# TODO rewrite this part using parametrized queries to protect from SQL injections?
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
      blast_hit => {
				blasteval => E,
				taxname_of_hit => S,
			},
      etc.
    ]
		orthoid2 => [
      blast_hit,
      blast_hit,
		]
  }

Arguments: scalar int SPECIESID, scalar int SETID

Returns: hashref of hashrefs of arrayrefs of hashrefs - lol

=cut
sub get_hitlist_hashref {
	my $specid = shift(@_) or croak("Usage: get_hitlist_for(SPECIESID, SETID)");
	my $setid  = shift(@_) or croak("Usage: get_hitlist_for(SPECIESID, SETID)");
	# TODO rewrite this part using parametrized queries to protect from SQL injections?
	my $query = "SELECT DISTINCT
		$mysql_table_hmmsearch.evalue,
		$mysql_table_orthologs.ortholog_gene_id, 
		$mysql_table_hmmsearch.target,
		$mysql_table_hmmsearch.start,
		$mysql_table_hmmsearch.end,
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
		WHERE $mysql_table_set_details.id = $setid
		AND $mysql_table_hmmsearch.taxid  = $specid
		";
	my $data = &mysql_get($query);
	my $result = { };
	foreach my $line ( @$data ) {
		push( @{ $result->{$$line[0]}->{$$line[1]} }, {
			'hmmhit'       => $$line[2],
			'start'        => $$line[3],
			'end'          => $$line[4],
			'blast_hit'    => $$line[5],
			'blast_evalue' => $$line[6],
			'species_name' => $$line[7],
		});
	}
	return $result;
}

sub get_hit_transcripts {
	my $specid = shift(@_) or croak("Usage: get_hitlist_for(SPECIESID, SETID)");
	my $setid  = shift(@_) or croak("Usage: get_hitlist_for(SPECIESID, SETID)");
	# TODO rewrite this part using parametrized queries to protect from SQL injections?
	my $query = "SELECT DISTINCT
		$mysql_table_hmmsearch.target
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
		WHERE $mysql_table_set_details.id = $setid
		AND $mysql_table_hmmsearch.taxid  = $specid
	";
	my $data = &mysql_get($query);
	my @result;
	push(@result, ${shift(@$data)}[0]) while @$data;
	return @result;
}
	
1;
