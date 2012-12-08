=head1 NAME 

Wrapper::Mysql

=head1 SYNOPSIS

  use Wrapper::Mysql;

  my $dbh = Wrapper::Mysql::mysql_dbh();
  do_stuff_with_dbh();
  undef $dbh;

  Wrapper::Mysql::do($query);

  my $result = Wrapper::Mysql::get($query);

=head1 DESCRIPTION

Wrapper module that provides MySQL functions. Deals with specific queries to
handle the complex database structure in the Orthograph pipeline.

=cut

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

# MySQL settings
my $mysql_dbname               = $config->{'mysql_database'};
my $mysql_dbpwd                = $config->{'mysql_password'};
my $mysql_dbserver             = $config->{'mysql_server'};
my $mysql_dbuser               = $config->{'mysql_username'};

my $mysql_table_blast          = $config->{'mysql_table_blast'};
my $mysql_table_blastdbs       = $config->{'mysql_table_blastdbs'};
my $mysql_table_ests           = $config->{'mysql_table_ests'};
my $mysql_table_ogs            = $config->{'mysql_table_ogs'};
my $mysql_table_hmmsearch      = $config->{'mysql_table_hmmsearch'};
my $mysql_table_log_evalues    = $config->{'mysql_table_log_evalues'};
my $mysql_table_set_details    = $config->{'mysql_table_set_details'};
my $mysql_table_aaseqs         = $config->{'mysql_table_aaseqs'};
my $mysql_table_seqpairs       = $config->{'mysql_table_sequence_pairs'};
my $mysql_table_taxa           = $config->{'mysql_table_taxa'};
my $mysql_table_orthologs      = $config->{'mysql_table_orthologs'};
my $mysql_col_aaseq            = 'aa_seq';
my $mysql_col_digest           = 'digest';
my $mysql_col_end              = 'end';
my $mysql_col_evalue           = 'evalue';
my $mysql_col_header           = 'header';
my $mysql_col_id               = 'id';
my $mysql_col_log_evalue       = 'log_evalue';
my $mysql_col_name             = 'name';
my $mysql_col_orthoid          = 'ortholog_gene_id';
my $mysql_col_query            = 'query';
my $mysql_col_setid            = 'setid';
my $mysql_col_sequence         = 'sequence';
my $mysql_col_start            = 'start';
my $mysql_col_target           = 'target';
my $mysql_col_taxid            = 'taxid';
my $outdir                     = $config->{'output_directory'};
my $orthoset                   = $config->{'ortholog_set'};
my $quiet                      = $config->{'quiet'};
my $reftaxa                    = $config->{'reference_taxa'};
# substitution character for selenocysteine, which normally leads to blast freaking out
my $u_subst                    = $config->{'substitute_u_with'};
my $sets_dir                   = $config->{'sets_dir'};
my $species_name               = $config->{'species_name'};
my $verbose                    = $config->{'verbose'};
#}}}

=head1 FUNCTIONS

=head2 mysql_dbh()

Get a database handle

Arguments: -

Returns: Database handle

=cut

sub mysql_dbh {#{{{
	return DBI->connect("DBI:mysql:$mysql_dbname:$mysql_dbserver;mysql_local_infile=1", $mysql_dbuser, $mysql_dbpwd);
}#}}}

=head2 mysql_get($query)

Get from the database the result of a SQL query

Expects: QUERY as a string literal

Returns: Reference to array of arrays (result lines->fields)

=cut

sub mysql_get {#{{{
	my $query = shift;
	unless ($query) { croak "Usage: mysql_get(QUERY)\n" }
  # prepare anonymous array
	my $results = [ ];
  # connect and fetch stuff
	my $dbh = &mysql_dbh();
	my $sql = $dbh->prepare($query);
	$sql->execute() or return 0;
	while (my @result = $sql->fetchrow_array() ) {
		push(@$results, \@result);
	}
	$sql->finish();
	$dbh->disconnect; # disconnect ASAP
	return $results;
}#}}}

=head2 mysql_do($query)

Connect to a database, execute a single $query (for repetitive queries, you
better do that by hand for performance reasons).

Expects: scalar string SQL query. 

Returns 1 on result, dies otherwise.

=cut

sub mysql_do {#{{{
	my $query = shift;
	unless ($query) { croak "Usage: mysql_do(QUERY)\n" }
	my @fields = @_;
	my $dbh = &mysql_dbh();
	my $sql = $dbh->prepare($query);
	$sql->execute(@fields) or die;
	$dbh->disconnect();
	return 1;
}#}}}

=head2 get_ortholog_sets()

Get list of ortholog sets from the database

Arguments: none

Returns: hash reference of set names => description

=cut

sub get_ortholog_sets {#{{{
	my %sets = ();
	my $query = "SELECT * FROM $mysql_table_set_details";
	my $data = &Wrapper::Mysql::mysql_get($query);
	foreach my $item (@$data) {
		$sets{$$item[1]} = $$item[2];
	}
	return(\%sets);
}#}}}

=head2 list_ogs

Get list of OGS in the database

Arguments: none

Returns: array reference (list of OGS)

=cut

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


=head2 get_ortholog_groups_for_set($setid)

Returns a hashref of hashrefs to create an ortholog set from. Each key in the hashref (the ortholog group ID) is a hashref of sequence_ID => sequence.

=cut

sub get_ortholog_groups_for_set {
	my $setid = shift @_ or croak "Usage: Wrapper::Mysql::get_ortholog_groups_for_set(SETID)";
	my $data = {};
	my $query = "SELECT o.ortholog_gene_id, a.id, a.sequence
		FROM $mysql_table_orthologs         AS o
    INNER JOIN $mysql_table_seqpairs    AS p
    ON o.sequence_pair = p.id
    INNER JOIN $mysql_table_aaseqs      AS a
    ON a.id = p.aa_seq
    INNER JOIN $mysql_table_set_details AS d
    ON d.id = o.setid
    WHERE d.id = ?";

	my $dbh = &mysql_dbh();
	my $sql = $dbh->prepare($query);
	$sql->execute( $setid ) or die;
	while (my @row = $sql->fetchrow_array()) {
		# load the whole set into memory, i don't give a frak
		$$data{$row[0]}{$row[1]} = $row[2];
	}
	$dbh->disconnect();	# disc asap

	return $data;
}

sub get_transcripts {
	my $specid = shift or croak "Usage: Wrapper::Mysql::get_transcripts(SPECIESID)";
	# TODO rewrite this part using parametrized queries to protect from SQL injections?
	my $query = "SELECT digest, sequence
		FROM $mysql_table_ests
		WHERE taxid = ?";
	my $dbh = Wrapper::Mysql::mysql_dbh();
	my $sth = $dbh->prepare($query);
	$sth->execute($specid);
	my $data = $sth->fetchall_arrayref();
	$sth->finish();
	$dbh->disconnect();
	return $data;
}

# Sub: get_hmmresults
# Get hmmsearch results from the database.
# Arguments: scalar string hmmsearch query
# Returns: reference to array of arrays
sub get_hmmresults {#{{{
	my $hmmquery = shift or croak "Usage: Wrapper::Mysql::get_hmmresults(HMMQUERY)";
	my $query_get_sequences = "SELECT $mysql_table_ests.digest,
		  $mysql_table_ests.sequence,
		  $mysql_table_hmmsearch.start,
		  $mysql_table_hmmsearch.end
		FROM $mysql_table_ests 
		INNER JOIN $mysql_table_hmmsearch
		ON $mysql_table_hmmsearch.target = $mysql_table_ests.digest
		WHERE $mysql_table_hmmsearch.query = ?";

	# get the sequences from the database (as array->array reference)
	my $dbh = &mysql_dbh();
	my $sth = $dbh->prepare($query_get_sequences);
	$sth->execute($hmmquery);
	my $results = $sth->fetchall_arrayref();
	$sth->finish();
	$dbh->disconnect();
	return $results;
}#}}}

# Sub: reciprocal_match
# Test whether this ortholog id led to a reciprocal match in the database
# Arguments: ortholog id, est digest
# Returns: 1 on match, 0 otherwise
sub reciprocal_match {
}

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

=head2 set_exists(SETNAME)

Tests whether a named set exists in the database. Returns 1 on success (the set
exists), 0 otherwise.

=cut

sub set_exists {
	my $set = shift;
	my $dbh = &mysql_dbh();
	my $sth = $dbh->prepare("SELECT * FROM $mysql_table_set_details WHERE $mysql_col_name = ? LIMIT 1");
	$sth->execute($set);
	my $result = $sth->fetchrow_arrayref;
	if ( $$result[0] ) { return 1 }
	return 0;
}

=head2 insert_taxon_into_table(TAXON_NAME)

Inserts a (non-core) taxon into the database. The taxon shorthand will be NULL
and the 'core' switch will be 0.

Returns the newly generated taxon ID.

=cut

sub insert_taxon_into_table {
	my $species_name = shift(@_);
	unless ($species_name) { croak("Usage: Wrapper::Mysql::insert_taxon_into_table(SPECIESNAME)") }
	my $query = "INSERT IGNORE INTO $mysql_table_taxa (longname, core) VALUES (?, ?)";
	my $dbh = &mysql_dbh();
	my $sth = $dbh->prepare($query);
	$sth->execute( $species_name, 0 )
		or return 0;
	$dbh->disconnect();
	return &get_taxid_for_species($species_name);
}

sub create_log_evalues_view {
	unless (scalar @_ == 1) { croak 'Usage: Wrapper::Mysql::create_log_evalues_view($species_id)' }
	my $taxid = shift;
	my $query_create_log_evalues = "CREATE OR REPLACE VIEW $mysql_table_log_evalues AS
	  SELECT $mysql_table_hmmsearch.$mysql_col_log_evalue AS $mysql_col_log_evalue,
	    COUNT($mysql_table_hmmsearch.$mysql_col_log_evalue) AS `count`
	  FROM $mysql_table_hmmsearch
	  WHERE $mysql_table_hmmsearch.$mysql_col_taxid = ?
	  GROUP BY $mysql_table_hmmsearch.$mysql_col_log_evalue
	  ORDER BY $mysql_table_hmmsearch.$mysql_col_log_evalue";
	my $dbh = &mysql_dbh();
	my $sth = $dbh->prepare($query_create_log_evalues);
	$sth->execute( $taxid ) or return 0;
	$dbh->disconnect();
	return 1;
}
	

# get a orthoid => list_of_aaseq_ids relationship from the db
sub get_orthologs_for_set_hashref {
	my $setid = shift(@_);
	unless ($setid) { croak("Usage: get_orthologs_for_set(SETID)") }
	my $query = "SELECT DISTINCT
		$mysql_table_orthologs.ortholog_gene_id,
		$mysql_table_aaseqs.id 
		FROM $mysql_table_orthologs 
		INNER JOIN $mysql_table_seqpairs 
			ON $mysql_table_orthologs.sequence_pair = $mysql_table_seqpairs.id
		INNER JOIN $mysql_table_aaseqs
			ON $mysql_table_seqpairs.aa_seq = $mysql_table_aaseqs.id
		INNER JOIN $mysql_table_set_details 
			ON $mysql_table_orthologs.setid = $mysql_table_set_details.id
		WHERE $mysql_table_set_details.id = ?";
	my $dbh = &mysql_dbh();
	my $sth = $dbh->prepare($query);
	$sth->execute( $setid );
	my $result = { };
	while (my $line = $sth->fetchrow_arrayref()) {
		push( @{$$result{$$line[0]}}, $$line[1] );
	}
	$sth->finish;
	$dbh->disconnect;
	return $result;
}

=head2 get_hitlist_hashref(SPECIESID, SETID)

Get the results in the form:

  evalue => {
    orthoid => [
      blast_hit => {
        blasteval => e-value,
        taxname_of_hit => string,
      },
      etc.
    ],
    orthoid2 => [
      blast_hit,
      blast_hit,
    ]
  }

Arguments: scalar int SPECIESID, scalar int SETID

Returns: hashref of hashrefs of arrayrefs of hashrefs - lol

=cut

sub get_hitlist_hashref {
	scalar @_ == 4 or croak("Usage: get_hitlist_for(SPECIESID, SETID, LIMIT, OFFSET)");
	my ($specid, $setid, $limit, $offset) = @_;
	my $query = "SELECT DISTINCT
		$mysql_table_hmmsearch.evalue,
		$mysql_table_orthologs.ortholog_gene_id, 
		$mysql_table_hmmsearch.target,
		$mysql_table_ests.header,
		$mysql_table_ests.sequence,
		$mysql_table_hmmsearch.start,
		$mysql_table_hmmsearch.end,
		$mysql_table_blast.target,
		$mysql_table_blast.evalue,
		$mysql_table_taxa.name
		FROM $mysql_table_hmmsearch
		INNER JOIN $mysql_table_ests
			ON $mysql_table_hmmsearch.target = $mysql_table_ests.digest
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
		WHERE $mysql_table_set_details.id = ?
		AND $mysql_table_hmmsearch.taxid  = ?
		ORDER BY $mysql_table_hmmsearch.log_evalue ASC
		LIMIT $limit 
		OFFSET $offset
		";
	print "fetching:\n$query\n";
	my $dbh = &mysql_dbh();
	my $sth = $dbh->prepare($query);
	$sth->execute( $setid, $specid ) or croak;
	my $result = { };
	while (my $line = $sth->fetchrow_arrayref()) {
		my $start = $$line[5] - 1;
		my $length = $$line[6] - $start;
		# first key is the hmmsearch evalue, second key is the orthoid
		push( @{ $result->{$$line[0]}->{$$line[1]} }, {
			'hmmhit'       => $$line[2],
			'header'       => $$line[3],
			'sequence'     => substr($$line[4], $start, $length),
			'start'        => $$line[5],
			'end'          => $$line[6],
			'blast_hit'    => $$line[7],
			'blast_evalue' => $$line[8],
			'species_name' => $$line[9],
		});
	}
	$sth->finish();
	$dbh->disconnect();
	scalar keys %$result > 0 ? return $result : return 0;
}

=head2 get_logevalue_count()

Returns a hashref as $hash->{$log_evalue} = number_of_occurences (an int)

=cut

sub get_logevalue_count {
	my $query_get_logevalues = "SELECT $mysql_col_log_evalue, count FROM $mysql_table_log_evalues";
	my $dbh = &mysql_dbh();
	my $sth = $dbh->prepare($query_get_logevalues);
	$sth->execute();
	my $d = $sth->fetchall_arrayref();
	$sth->finish();
	$dbh->disconnect();
	my $num_of_logevalues = { };
	foreach my $row (@$d) {
		$num_of_logevalues->{$$row[0]} = $$row[1];
	}
	return $num_of_logevalues;
}

=head2 get_results_for_logevalue_range($setid, $taxonid, $min, $max)

Fetch results from the database that are BETWEEN $min and $max in terms of log-evalue.

Returns a [ rows->[ (columns) ] ] arrayref.

=cut

sub get_results_for_logevalue_range {
	my ($setid, $taxid, $min, $max) = @_;
	# complex parametrized query
	my $query = "SELECT DISTINCT $mysql_table_hmmsearch.$mysql_col_evalue,
			$mysql_table_orthologs.$mysql_col_orthoid,
			$mysql_table_hmmsearch.$mysql_col_target,
			$mysql_table_ests.$mysql_col_header,
			$mysql_table_ests.$mysql_col_sequence,
			$mysql_table_hmmsearch.$mysql_col_start,
			$mysql_table_hmmsearch.$mysql_col_end,
			$mysql_table_blast.$mysql_col_target,
			$mysql_table_blast.$mysql_col_evalue,
			$mysql_table_taxa.$mysql_col_name
		FROM $mysql_table_log_evalues
		LEFT JOIN $mysql_table_hmmsearch
			ON $mysql_table_log_evalues.$mysql_col_log_evalue = $mysql_table_hmmsearch.$mysql_col_log_evalue
		LEFT JOIN $mysql_table_ests
			ON $mysql_table_hmmsearch.$mysql_col_target = $mysql_table_ests.$mysql_col_digest
		LEFT JOIN $mysql_table_orthologs
			ON $mysql_table_hmmsearch.$mysql_col_query = $mysql_table_orthologs.$mysql_col_orthoid
		LEFT JOIN $mysql_table_blast
			ON $mysql_table_hmmsearch.$mysql_col_target = $mysql_table_blast.$mysql_col_query
		LEFT JOIN $mysql_table_aaseqs
			ON $mysql_table_blast.$mysql_col_target = $mysql_table_aaseqs.$mysql_col_id
		LEFT JOIN $mysql_table_taxa
			ON $mysql_table_aaseqs.$mysql_col_taxid = $mysql_table_taxa.$mysql_col_id
		LEFT JOIN $mysql_table_set_details
			ON $mysql_table_orthologs.$mysql_col_setid = $mysql_table_set_details.$mysql_col_id
		WHERE $mysql_table_hmmsearch.$mysql_col_log_evalue IS NOT NULL
			AND $mysql_table_ests.$mysql_col_digest          IS NOT NULL
			AND $mysql_table_orthologs.$mysql_col_orthoid    IS NOT NULL
			AND $mysql_table_blast.$mysql_col_query          IS NOT NULL
			AND $mysql_table_aaseqs.$mysql_col_id            IS NOT NULL
			AND $mysql_table_taxa.$mysql_col_id              IS NOT NULL
			AND $mysql_table_set_details.$mysql_col_id       IS NOT NULL
			AND $mysql_table_set_details.$mysql_col_id       = ?
			AND $mysql_table_hmmsearch.$mysql_col_taxid      = ?
			AND $mysql_table_hmmsearch.$mysql_col_log_evalue BETWEEN ? AND ?";
	my $dbh = &mysql_dbh();
	my $sth = $dbh->prepare($query);
	$sth->execute( $setid, $taxid, $min, $max );
	my $result = { };
	while (my $line = $sth->fetchrow_arrayref()) {
		my $start = $$line[5] - 1;
		my $length = $$line[6] - $start;
		# first key is the hmmsearch evalue, second key is the orthoid
		push( @{ $result->{$$line[0]}->{$$line[1]} }, {
			'hmmhit'       => $$line[2],
			'header'       => $$line[3],
			'sequence'     => substr($$line[4], $start, $length),
			'start'        => $$line[5],
			'end'          => $$line[6],
			'blast_hit'    => $$line[7],
			'blast_evalue' => $$line[8],
			'species_name' => $$line[9],
		});
	}
	$sth->finish();
	$dbh->disconnect();
	scalar keys %$result > 0 ? return $result : return 0;
}

=head2 get_results_for_logevalue($setid, $taxonid, $logevalue)

Fetch results from the database that have $logevalue.

Returns a [ rows->[ (columns) ] ] arrayref.

=cut

sub get_results_for_logevalue {
	my $setid   = shift;
	my $taxid   = shift;
	my $logeval = shift;
	my $query = "SELECT DISTINCT $mysql_table_hmmsearch.$mysql_col_evalue,
			$mysql_table_orthologs.$mysql_col_orthoid,
			$mysql_table_hmmsearch.$mysql_col_target,
			$mysql_table_ests.$mysql_col_header,
			$mysql_table_ests.$mysql_col_sequence,
			$mysql_table_hmmsearch.$mysql_col_start,
			$mysql_table_hmmsearch.$mysql_col_end,
			$mysql_table_blast.$mysql_col_target,
			$mysql_table_blast.$mysql_col_evalue,
			$mysql_table_taxa.$mysql_col_name
		FROM $mysql_table_log_evalues
		LEFT JOIN $mysql_table_hmmsearch
			ON $mysql_table_log_evalues.$mysql_col_log_evalue = $mysql_table_hmmsearch.$mysql_col_log_evalue
		LEFT JOIN $mysql_table_ests
			ON $mysql_table_hmmsearch.$mysql_col_target = $mysql_table_ests.$mysql_col_digest
		LEFT JOIN $mysql_table_orthologs
			ON $mysql_table_hmmsearch.$mysql_col_query = $mysql_table_orthologs.$mysql_col_orthoid
		LEFT JOIN $mysql_table_blast
			ON $mysql_table_hmmsearch.$mysql_col_target = $mysql_table_blast.$mysql_col_query
		LEFT JOIN $mysql_table_aaseqs
			ON $mysql_table_blast.$mysql_col_target = $mysql_table_aaseqs.$mysql_col_id
		LEFT JOIN $mysql_table_taxa
			ON $mysql_table_aaseqs.$mysql_col_taxid = $mysql_table_taxa.$mysql_col_id
		LEFT JOIN $mysql_table_set_details
			ON $mysql_table_orthologs.$mysql_col_setid = $mysql_table_set_details.$mysql_col_id
		WHERE $mysql_table_hmmsearch.$mysql_col_log_evalue IS NOT NULL
			AND $mysql_table_ests.$mysql_col_digest          IS NOT NULL
			AND $mysql_table_orthologs.$mysql_col_orthoid    IS NOT NULL
			AND $mysql_table_blast.$mysql_col_query          IS NOT NULL
			AND $mysql_table_aaseqs.$mysql_col_id            IS NOT NULL
			AND $mysql_table_taxa.$mysql_col_id              IS NOT NULL
			AND $mysql_table_set_details.$mysql_col_id       IS NOT NULL
			AND $mysql_table_set_details.$mysql_col_id       = ?
			AND $mysql_table_hmmsearch.$mysql_col_taxid      = ?
			AND $mysql_table_hmmsearch.$mysql_col_log_evalue = ?";
	my $dbh = &mysql_dbh();
	my $sth = $dbh->prepare($query);
	$sth->execute( $setid, $taxid, $logeval );
	my $result = { };
	while (my $line = $sth->fetchrow_arrayref()) {
		my $start = $$line[5] - 1;
		my $length = $$line[6] - $start;
		# first key is the hmmsearch evalue, second key is the orthoid
		push( @{ $result->{$$line[0]}->{$$line[1]} }, {
			'hmmhit'       => $$line[2],
			'header'       => $$line[3],
			'sequence'     => substr($$line[4], $start, $length),
			'start'        => $$line[5],
			'end'          => $$line[6],
			'blast_hit'    => $$line[7],
			'blast_evalue' => $$line[8],
			'species_name' => $$line[9],
		});
	}
	$sth->finish();
	$dbh->disconnect();
	scalar keys %$result > 0 ? return $result : return 0;
}

=head2 get_hit_transcripts($species_id, $set_id)

Returns a list of transcript digests that were hit during the HMM search for
this $species_id and this $set_id. 

=cut

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
