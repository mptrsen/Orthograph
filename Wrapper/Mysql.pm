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
=head1 NAME 

Wrapper::Mysql

=head1 SYNOPSIS

  use Wrapper::Mysql;

  my $dbh = Wrapper::Mysql::get_dbh();
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
use DBI;            # database interface
use DBD::mysql;     # MySQL database driver

my $config = $Orthograph::Config::config;  # copy config

# MySQL settings
my $mysql_dbname               = $config->{'mysql-database'};
my $mysql_dbpwd                = $config->{'mysql-password'};
my $mysql_dbserver             = $config->{'mysql-server'};
my $mysql_dbuser               = $config->{'mysql-username'};
my $mysql_timeout              = $config->{'mysql-timeout'};
my $sleep_for                  = 10;

my $mysql_table_aaseqs         = $config->{'db_table_aaseqs'};
my $mysql_table_blast          = $config->{'db_table_blast'};
my $mysql_table_blastdbs       = $config->{'db_table_blastdbs'};
my $mysql_table_ests           = $config->{'db_table_ests'};
my $mysql_table_hmmsearch      = $config->{'db_table_hmmsearch'};
my $mysql_table_log_evalues    = $config->{'db_table_log_evalues'};
my $mysql_table_scores         = $config->{'db_table_scores'};
my $mysql_table_ntseqs         = $config->{'db_table_ntseqs'};
my $mysql_table_ogs            = $config->{'db_table_ogs'};
my $mysql_table_orthologs      = $config->{'db_table_orthologs'};
my $mysql_table_seqpairs       = $config->{'db_table_sequence_pairs'};
my $mysql_table_set_details    = $config->{'db_table_set_details'};
my $mysql_table_taxa           = $config->{'db_table_taxa'};
my $mysql_col_aaseq            = 'aa_seq';
my $mysql_col_digest           = 'digest';
my $mysql_col_end              = 'end';
my $mysql_col_env_end          = 'env_end';
my $mysql_col_env_start        = 'env_start';
my $mysql_col_evalue           = 'evalue';
my $mysql_col_hmm_end          = 'hmm_end';
my $mysql_col_hmm_start        = 'hmm_start';
my $mysql_col_header           = 'header';
my $mysql_col_id               = 'id';
my $mysql_col_log_evalue       = 'log_evalue';
my $mysql_col_score            = 'score';
my $mysql_col_name             = 'name';
my $mysql_col_ntseq            = 'nt_seq';
my $mysql_col_orthoid          = 'ortholog_gene_id';
my $mysql_col_query            = 'query';
my $mysql_col_setid            = 'setid';
my $mysql_col_sequence         = 'sequence';
my $mysql_col_seqpair          = 'sequence_pair';
my $mysql_col_start            = 'start';
my $mysql_col_target           = 'target';
my $mysql_col_taxid            = 'taxid';
my $outdir                     = $config->{'output-directory'};
my $orthoset                   = $config->{'ortholog-set'};
my $quiet                      = $config->{'quiet'};
my $reftaxa                    = $config->{'reference-taxa'};
# substitution character for selenocysteine, which normally leads to blast freaking out
my $u_subst                    = $config->{'substitute-u-with'};
my $sets_dir                   = $config->{'sets-dir'};
my $species_name               = $config->{'species-name'};
my $g_species_id               = undef;	# global variable
my $verbose                    = $config->{'verbose'};
my $debug                      = $config->{'debug'};
#}}}

# Check whether all information was provided in the configuration
defined $mysql_dbname   or fail_and_exit('MySQL database name not specified');
defined $mysql_dbuser   or fail_and_exit('MySQL database username not specified');
defined $mysql_dbpwd    or fail_and_exit('MySQL database password not specified');
defined $mysql_dbserver or fail_and_exit('MySQL database server not specified');


=head1 FUNCTIONS

=cut

sub fail_and_exit {
	my $msg = shift @_;
	print STDERR 'Fatal: ' . $msg . "\n";
	exit 1;
}

=head2 get_dbh()

Get a database handle

Arguments: -

Returns: Database handle

=cut

sub get_dbh {#{{{
	my $dbh = undef;
	my $slept = 0;

	until ($dbh = DBI->connect("DBI:mysql:$mysql_dbname:$mysql_dbserver;mysql_local_infile=1", $mysql_dbuser, $mysql_dbpwd)) {
		if ($slept >= $mysql_timeout) { 
			carp "Warning: Connection retry timeout exceeded\n" and return undef;
		}
		carp "Warning: Connection failed, retrying in $sleep_for seconds\n";
		sleep $sleep_for;
		$slept += $sleep_for;
	}

	if ($dbh) { return $dbh }
	return undef;
}#}}}

=head2 mysql_get($query)

Get from the database the result of a SQL query

Expects: QUERY as a string literal

Returns: Reference to array of arrays (result lines->fields)

=cut

sub mysql_get {#{{{
	my $query = shift;
	unless ($query) { croak "Usage: mysql_get(QUERY, ARGS)\n" }
	my @args = @_;
  # prepare anonymous array
	my $results = [ ];
  # connect and fetch stuff
	my $dbh = get_dbh()
		or return undef;
	my $sth = $dbh->prepare($query);
	$sth = execute($sth, $mysql_timeout, @args);
	while (my @result = $sth->fetchrow_array() ) {
		push(@$results, \@result);
	}
	$sth->finish();
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
	my $dbh = get_dbh()
		or return undef;
	my $sth = $dbh->prepare($query);
	$sth = execute($sth, $mysql_timeout, @fields);
	$dbh->disconnect();
	return 1;
}#}}}

# Sub: check
# Check whether the result of a query is present or not, return appropriate
# Arguments: Scalar string QUERY
# Returns: 1 or 0 depending on presence of result
sub check {#{{{
	my $query = shift;
	unless ($query) { croak "Usage: check(QUERY)\n"; }
	my @results;
	my $dbh = get_dbh();
	my $sql = $dbh->prepare($query);
	$sql->execute();
	if ($sql->fetchrow_array()) {
		return 1;
	}
	return 0;
}#}}}

sub drop_tables {
	my %t = @_;
	print 'DROPing tables: ', join(", ", values(%t)), "\n" if $verbose;
	my $dbh = get_dbh() or fail_and_exit("Couldn't get database connection");
	foreach my $table (keys(%t)) {
		$dbh->do("DROP TABLE IF EXISTS $t{$table}") or die "Could not execute drop query: $!\n";
	}
	$dbh->disconnect();
}

sub create_tables {
	my %t = shift;
	# the queries for the individual tables
	my %create_table = (#{{{
		# table: blastdbs
		'blastdbs' => "CREATE TABLE `$t{'blastdbs'}` (
			`id`           INT UNSIGNED NOT NULL AUTO_INCREMENT PRIMARY KEY,
			`setid`        INT UNSIGNED DEFAULT NULL, UNIQUE(setid),
			`blastdb_path` VARCHAR(255) DEFAULT NULL)",
		
		# table: ogs
		'ogs' => "CREATE TABLE `$t{'ogs'}` (
			`id`           INT UNSIGNED NOT NULL AUTO_INCREMENT PRIMARY KEY,
			`type`         INT(1),
			`taxid`        INT UNSIGNED NOT NULL, UNIQUE(taxid),
			`version`      VARCHAR(255))",
		
		# table: ortholog_set
		'ortholog_set' => "CREATE TABLE `$t{'orthologs'}` (
			`id`               INT UNSIGNED NOT NULL AUTO_INCREMENT PRIMARY KEY,
			`setid`            INT UNSIGNED NOT NULL,
			`ortholog_gene_id` VARCHAR(10)  NOT NULL,
			`sequence_pair`    INT UNSIGNED NOT NULL,
			UNIQUE INDEX (setid, ortholog_gene_id, sequence_pair))",

		# table: sequence_pairs
		'sequence_pairs' => "CREATE TABLE `$t{'seqpairs'}` (
			`id`           BIGINT UNSIGNED NOT NULL AUTO_INCREMENT PRIMARY KEY,
			`taxid`        INT    UNSIGNED,
			`ogs_id`       INT    UNSIGNED,
			`aa_seq`       INT    UNSIGNED, UNIQUE(aa_seq),
			`nt_seq`       INT    UNSIGNED, UNIQUE(nt_seq), 
			`date`         INT    UNSIGNED,
			`user`         INT    UNSIGNED)",

		# table: sequences_aa
		'aa_sequences' => "CREATE TABLE `$t{'aaseqs'}` (
			`id`           BIGINT UNSIGNED NOT NULL AUTO_INCREMENT PRIMARY KEY,
			`taxid`        INT             NOT NULL, INDEX(taxid),
			`header`       VARCHAR(512),             INDEX(header), UNIQUE(header),
			`sequence`     MEDIUMBLOB,
			`user`         INT UNSIGNED,
			`date`         INT UNSIGNED)",

		# table: sequences_nt
		'nt_sequences' => "CREATE TABLE `$t{'ntseqs'}` (
			`id`           BIGINT UNSIGNED NOT NULL AUTO_INCREMENT PRIMARY KEY,
			`taxid`        INT             NOT NULL, INDEX(taxid),
			`header`       VARCHAR(512),             INDEX(header), UNIQUE(header),
			`sequence`     MEDIUMBLOB,
			`user`         INT UNSIGNED,
			`date`         INT UNSIGNED)",

		# table: set_details
		'set_details' => "CREATE TABLE `$t{'set_details'}` (
			`id`           INT UNSIGNED NOT NULL AUTO_INCREMENT PRIMARY KEY,
			`name`         VARCHAR(255), UNIQUE(name),
			`description`  BLOB)",

		# table: taxa
		'taxa' => "CREATE TABLE `$t{'taxa'}` (
			`id`           INT UNSIGNED NOT NULL AUTO_INCREMENT PRIMARY KEY,
			`name`         VARCHAR(20),  UNIQUE(name),
			`longname`     VARCHAR(255), 
			`core`         TINYINT UNSIGNED NOT NULL)",
		
		# table: users
		'users' => "CREATE TABLE `$t{'users'}` (
			`id`           INT UNSIGNED NOT NULL AUTO_INCREMENT PRIMARY KEY,
			`name`         VARCHAR(255), UNIQUE(name))",
		# table: seqtypes
		'seqtypes' => "CREATE TABLE `$t{'seqtypes'}` (
			`id`           INT UNSIGNED NOT NULL AUTO_INCREMENT PRIMARY KEY,
			`type`         CHAR(3),     UNIQUE(type))",
	);#}}}

	# to start off with nt and aa sequence types
	my $insert_seqtypes = "INSERT IGNORE INTO $t{'seqtypes'} (type) VALUES ('nt'),('aa')";

	my $dbh = get_dbh();
	foreach (values %create_table) {
		print $_, ";\n" if $verbose;
		$dbh->do($_) or die "Could not exec query: $!\n";
	}
	$dbh->do($insert_seqtypes);	# start off with 'nt' and 'aa' seqtypes
	$dbh->disconnect;
}

sub create_temp_table {
	my $temptable = shift @_;
	my $dbh = get_dbh();
	my $create_temp_table_query = "CREATE TABLE $temptable (
			`name`     VARCHAR(255), INDEX(name),
			`longname` VARCHAR(255),
			`orthoset` VARCHAR(255), INDEX(orthoset),
			`orthoid`  VARCHAR(255), INDEX(orthoid),
			`blastdb`  VARCHAR(255),
			`header`   VARCHAR(512), INDEX(header),
			`sequence` MEDIUMBLOB,
			`description` VARCHAR(255))";
	$dbh->do("DROP TABLE IF EXISTS $temptable") or die "Fatal: Could not DROP TABLE $temptable\n";
	$dbh->do($create_temp_table_query) or die "Fatal: Could not CREATE TABLE $temptable\n";
}

sub load_csv_into_temptable {
	my $csvfile   = shift @_;
	my $temptable = shift @_;
	my $loadquery = "LOAD DATA LOCAL INFILE '$csvfile' 
		INTO TABLE $temptable FIELDS TERMINATED BY ',' LINES TERMINATED BY '\n' (
			name,
			longname,
			orthoset,
			orthoid,
			blastdb,
			header,
			sequence,
			description)";
	my $dbh = get_dbh();
	$dbh->do($loadquery) or die "Fatal: Could not LOAD DATA into temporary table $temptable\n";
	$dbh->disconnect;
}

sub fill_tables_from_temp_table {
	my $t = shift @_;
	my $temptable = shift @_;
	my @queries = (
		# user name
		"INSERT IGNORE INTO $t->{'users'} (name) VALUES ('$mysql_dbuser')",
		# taxa (name, longname)
		"INSERT IGNORE INTO $t->{'taxa'} (name, longname, core) 
			SELECT DISTINCT $temptable.name, $temptable.longname, 1 
			FROM $temptable",
		# set name + description
		"INSERT IGNORE INTO $t->{'set_details'} (name, description)
			SELECT DISTINCT $temptable.orthoset, $temptable.description 
			FROM $temptable LIMIT 1",
		# blast databases
		"INSERT IGNORE INTO $t->{'blastdbs'} (setid, blastdb_path) 
			SELECT DISTINCT $t->{'set_details'}.id, $temptable.blastdb 
			FROM $temptable
			LEFT JOIN $t->{'set_details'} 
				ON $t->{'set_details'}.name = $temptable.orthoset",
		# pep sequences
		"INSERT IGNORE INTO $t->{'aaseqs'} (taxid, header, sequence, user, date) 
			SELECT $t->{'taxa'}.id, $temptable.header, $temptable.sequence, $t->{'users'}.id, UNIX_TIMESTAMP()
			FROM $temptable
				LEFT JOIN $t->{'taxa'} 
			ON $temptable.name  = $t->{'taxa'}.name
				INNER JOIN $t->{'users'}
			ON $t->{'users'}.name = '$mysql_dbuser'",
		# delete everything where header or sequence is NULL or empty
		"DELETE FROM $t->{'aaseqs'}
			WHERE $t->{'aaseqs'}.header IS NULL
			OR $t->{'aaseqs'}.sequence IS NULL
			OR $t->{'aaseqs'}.header = ''
			OR $t->{'aaseqs'}.sequence = ''",
		# sequence pairs (pep-nuc)
		"INSERT IGNORE INTO $t->{'seqpairs'} (taxid, ogs_id, aa_seq, nt_seq, date, user)
			SELECT $t->{'taxa'}.id, $t->{'ogs'}.id, $t->{'aaseqs'}.id, $t->{'ntseqs'}.id, UNIX_TIMESTAMP(), $t->{'users'}.id
			FROM $t->{'taxa'}
			INNER JOIN $t->{'aaseqs'}
				ON $t->{'aaseqs'}.taxid = $t->{'taxa'}.id
			LEFT JOIN $t->{'ogs'}
				ON $t->{'taxa'}.id = $t->{'ogs'}.taxid
			LEFT JOIN $t->{'ntseqs'}
				ON $t->{'aaseqs'}.header = $t->{'ntseqs'}.header
			INNER JOIN $t->{'users'}
				ON $t->{'users'}.name = '$mysql_dbuser'",
		# orthologous groups
		"INSERT IGNORE INTO $t->{'orthologs'} (setid, ortholog_gene_id, sequence_pair) 
			SELECT $t->{'set_details'}.id, $temptable.orthoid, $t->{'seqpairs'}.id 
			FROM $t->{'aaseqs'} 
			INNER JOIN $temptable 
				ON $t->{'aaseqs'}.header = $temptable.header 
			INNER JOIN $t->{'seqpairs'} 
				ON $t->{'seqpairs'}.aa_seq = $t->{'aaseqs'}.id 
			INNER JOIN $t->{'set_details'} 
				ON $t->{'set_details'}.name = $temptable.orthoset",
	);

	my $dbh = get_dbh();
	my $nrows;
	foreach (@queries) {
		print $_ . ";\n";
		$nrows = $dbh->do($_) or die();
		($nrows > 0) ? printf("Query OK, %d rows affected\n", $nrows) : print "Query OK\n";
	}
	$dbh->disconnect;
	return $nrows;
}


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

	my $dbh = get_dbh()
		or return undef;
	my $sth = $dbh->prepare($query);
	$sth = execute($sth, $mysql_timeout, $setid);
	while (my @row = $sth->fetchrow_array()) {
		# load the whole set into memory, i don't give a frak
		$$data{$row[0]}{$row[1]} = $row[2];
	}
	$dbh->disconnect();	# disc asap

	return $data;
}

sub get_transcripts {
	my $specid = shift or croak "Usage: Wrapper::Mysql::get_transcripts(SPECIESID, TYPE)";
	my $type = shift or croak "Usage: Wrapper::Mysql::get_transcripts(SPECIESID, TYPE)";
	my $query = "SELECT digest, sequence
		FROM $mysql_table_ests
		WHERE taxid = ?
		AND type = ?";
	my $dbh = get_dbh()
		or return undef;
	my $sth = $dbh->prepare($query);
	$sth = execute($sth, $mysql_timeout, $specid, $type);
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
	my ($hmmquery, $taxid) = @_ or croak "Usage: Wrapper::Mysql::get_hmmresults(HMMQUERY)";
	# disable query cache for this one
	my $query_get_sequences = "SELECT SQL_NO_CACHE $mysql_table_ests.digest,
		  $mysql_table_ests.sequence,
		  $mysql_table_hmmsearch.env_start,
		  $mysql_table_hmmsearch.env_end
		FROM $mysql_table_ests 
		INNER JOIN $mysql_table_hmmsearch
		ON $mysql_table_hmmsearch.target = $mysql_table_ests.digest
		WHERE $mysql_table_hmmsearch.query = ?
		AND $mysql_table_hmmsearch.taxid = ?";

	# get the sequences from the database (as array->array reference)
	my $dbh = get_dbh()
		or return undef;
	my $sth = $dbh->prepare($query_get_sequences);
	do {
		$sth->execute($hmmquery, $taxid);
	} while ($sth->err);
	my $results = $sth->fetchall_arrayref();
	$sth->finish();
	$dbh->disconnect();
	return $results;
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

=head2 get_number_of_ests_for_specid(ID)

Returns the number of EST sequences (transcripts) for a given species id.

Argument: scalar int ID

Returns: scalar int 

=cut

sub get_number_of_ests_for_specid {
	my $specid = shift @_ or croak "Usage: get_number_of_ests_for_specid(SPECID)";

	# TODO rewrite this part using parametrized queries to protect from SQL injections?
	my $result = &mysql_get("SELECT COUNT(*) FROM $mysql_table_ests");

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
	print "getting taxid for $species_name\n";
	unless ($species_name) { croak("Usage: get_taxid_for_species(SPECIESNAME)") }
	# TODO rewrite this part using parametrized queries to protect from SQL injections?
	my $query = "SELECT id FROM $mysql_table_taxa WHERE core = 0 AND longname = '$species_name'";
	my $result = mysql_get($query);
	if ($result) { 
		$g_species_id = $$result[0][0];
		return $$result[0][0];
	}
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
	my $dbh = get_dbh()
		or return undef;
	my $sth = $dbh->prepare("SELECT * FROM $mysql_table_set_details WHERE $mysql_col_name = ? LIMIT 1");
	$sth = execute($sth, $mysql_timeout, $set);
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
	if (my $taxid = &get_taxid_for_species($species_name)) { return $taxid }
	my $query = "INSERT IGNORE INTO $mysql_table_taxa (longname, core) VALUES (?, ?)";
	my $dbh = get_dbh()
		or return undef;
	my $sth = $dbh->prepare($query);
	$sth = execute($sth, $mysql_timeout, $species_name, 0);
	$dbh->disconnect();

	$g_species_id = &get_taxid_for_species($species_name) or croak;
	return $g_species_id;
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
	my $dbh = get_dbh()
		or return undef;
	my $sth = $dbh->prepare($query_create_log_evalues);
	$sth = execute($sth, $mysql_timeout, $taxid);
	$dbh->disconnect();
	return 1;
}
	
sub create_scores_view {
	unless (scalar @_ == 1) { croak 'Usage: Wrapper::Mysql::create_scores_view($species_id)' }
	my $taxid = shift;
	my $query_create_scores_view = "CREATE OR REPLACE VIEW $mysql_table_scores AS
	  SELECT $mysql_table_hmmsearch.$mysql_col_score AS $mysql_col_score,
	    COUNT($mysql_table_hmmsearch.$mysql_col_score) AS `count`
	  FROM $mysql_table_hmmsearch
	  WHERE $mysql_table_hmmsearch.$mysql_col_taxid = ?
	  GROUP BY $mysql_table_hmmsearch.$mysql_col_score
	  ORDER BY $mysql_table_hmmsearch.$mysql_col_score DESC";
	my $dbh = get_dbh()
		or return undef;
	my $sth = $dbh->prepare($query_create_scores_view);
	$sth = execute($sth, $mysql_timeout, $taxid);
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
	my $dbh = get_dbh()
		or return undef;
	my $sth = $dbh->prepare($query);
	$sth = execute($sth, $mysql_timeout, $setid);
	my $result = { };
	while (my $line = $sth->fetchrow_arrayref()) {
		push( @{$$result{$$line[0]}}, $$line[1] );
	}
	$sth->finish;
	$dbh->disconnect;
	return $result;
}

=head2 get_ortholog_group($orthoid)

Get a specific ortholog group, i.e. aa headers and sequences.

=cut

sub get_ortholog_group {
	my $setid   = shift;
	my $orthoid = shift;
	my $query = "SELECT 
		$mysql_table_aaseqs.$mysql_col_header, $mysql_table_aaseqs.$mysql_col_sequence
		FROM $mysql_table_aaseqs
		INNER JOIN $mysql_table_seqpairs
			ON $mysql_table_aaseqs.$mysql_col_id = $mysql_table_seqpairs.$mysql_col_aaseq
		INNER JOIN $mysql_table_orthologs
			ON $mysql_table_seqpairs.$mysql_col_id = $mysql_table_orthologs.$mysql_col_seqpair
		AND   $mysql_table_orthologs.$mysql_col_setid = ?
		AND   $mysql_table_orthologs.$mysql_col_orthoid = ?";
	my $dbh = get_dbh()
		or return undef;
	my $sth = $dbh->prepare($query);
	$sth = execute($sth, $mysql_timeout, $setid, $orthoid);
	my $data = $sth->fetchall_arrayref();
	return $data;
}

sub get_ortholog_group_nucleotide {
	my $setid   = shift;
	my $orthoid = shift;
	my $query = "SELECT 
		$mysql_table_ntseqs.$mysql_col_header, $mysql_table_ntseqs.$mysql_col_sequence
		FROM $mysql_table_ntseqs
		INNER JOIN $mysql_table_seqpairs
			ON $mysql_table_ntseqs.$mysql_col_id = $mysql_table_seqpairs.$mysql_col_ntseq
		INNER JOIN $mysql_table_orthologs
			ON $mysql_table_seqpairs.$mysql_col_id = $mysql_table_orthologs.$mysql_col_seqpair
		AND   $mysql_table_orthologs.$mysql_col_setid = ?
		AND   $mysql_table_orthologs.$mysql_col_orthoid = ?";
	my $dbh = get_dbh()
		or return undef;
	my $sth = $dbh->prepare($query);
	$sth = execute($sth, $mysql_timeout, $setid, $orthoid);
	my $data = $sth->fetchall_arrayref();
	return $data;
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
		$mysql_table_hmmsearch.hmm_start,
		$mysql_table_hmmsearch.hmm_end,
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
	my $dbh = get_dbh()
		or return undef;
	my $sth = $dbh->prepare($query);
	$sth = execute($sth, $mysql_timeout, $setid, $specid);
	my $result = { };
	while (my $line = $sth->fetchrow_arrayref()) {
		my $start = $$line[5] - 1;
		my $length = $$line[6] - $start;
		# first key is the hmmsearch evalue, second key is the orthoid
		push( @{ $result->{$$line[0]}->{$$line[1]} }, {
			'hmmhit'       => $$line[2],
			'header'       => $$line[3],
			'sequence'     => substr($$line[4], $start, $length),
			'hmm_start'        => $$line[5],
			'hmm_end'          => $$line[6],
			'blast_hit'    => $$line[7],
			'blast_evalue' => $$line[8],
			'reftaxon'     => $$line[9],
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
	my $dbh = get_dbh()
		or return undef;
	my $sth = $dbh->prepare($query_get_logevalues);
	$sth = execute($sth, $mysql_timeout);
	my $d = $sth->fetchall_arrayref();
	$sth->finish();
	$dbh->disconnect();
	my $num_of_logevalues = { };
	foreach my $row (@$d) {
		$num_of_logevalues->{$$row[0]} = $$row[1];
	}
	return $num_of_logevalues;
}

sub get_scores_count {
	my $query_get_scores = "SELECT $mysql_col_score, count FROM $mysql_table_scores";
	my $dbh = get_dbh()
		or return undef;
	my $sth = $dbh->prepare($query_get_scores);
	$sth = execute($sth, $mysql_timeout);
	my $d = $sth->fetchall_arrayref();
	$sth->finish();
	$dbh->disconnect();
	my $num_of_scores = { };
	foreach my $row (@$d) {
		$num_of_scores->{$$row[0]} = $$row[1];
	}
	return $num_of_scores;
}

sub execute {
	my $sth = shift or croak "Usage: execute(STH, TIMEOUT, ARGS)\n";
	my $timeout = shift or croak "Usage: execute(STH, TIMEOUT, ARGS)\n";
	my @args = @_;
	my $slept = 0;
	until ($sth->execute(@args)) {
		carp "Warning: execution failed, retrying in $sleep_for seconds...\n";
		if ($slept > $timeout) { croak "Fatal: execution ultimately failed, failing this transaction\n" }
		sleep $sleep_for;
		$slept += $sleep_for;
	}
	return $sth;
}

=head2 get_results_for_logevalue($setid, $taxonid, $min [, $max])

Fetch results from the database that have e-value $min or are BETWEEN $min AND $max.

The function intelligently does the correct query depending on the number of arguments.

Returns a [ rows->[ (columns) ] ] arrayref.

=cut

sub get_results_for_logevalue {
	my $setid   = shift;
	my $taxid   = shift;
	my $min     = shift;
	my $max     = shift;
	# generic query
	my $query = "SELECT DISTINCT $mysql_table_hmmsearch.$mysql_col_evalue,
			$mysql_table_orthologs.$mysql_col_orthoid,
			$mysql_table_hmmsearch.$mysql_col_target,
			$mysql_table_ests.$mysql_col_header,
			$mysql_table_ests.$mysql_col_sequence,
			$mysql_table_hmmsearch.$mysql_col_hmm_start,
			$mysql_table_hmmsearch.$mysql_col_hmm_end,
			$mysql_table_hmmsearch.$mysql_col_env_start,
			$mysql_table_hmmsearch.$mysql_col_env_end,
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
			AND $mysql_table_hmmsearch.$mysql_col_taxid      = ?";

	# modify the generic query
	# e-value range
	if ($max) { $query .= "\n			AND $mysql_table_hmmsearch.$mysql_col_log_evalue BETWEEN ? AND ?" }
	# single e-value
	else      { $query .= "\n			AND $mysql_table_hmmsearch.$mysql_col_log_evalue = ?" }

	# good for debugging
	print $query . "\n" if $debug;

	my $dbh = get_dbh()
		or return undef;
	my $sth = $dbh->prepare($query);

	# e-value range
	if ($max) {
		$sth = execute($sth, $mysql_timeout, $setid, $taxid, $min, $max);
	}
	# single e-value
	else      {
		$sth = execute($sth, $mysql_timeout, $setid, $taxid, $min);
	} 

	# will hold the result
	my $result = { };

	while (my $line = $sth->fetchrow_arrayref()) {
		my $start = $$line[7] - 1;
		my $length = $$line[8] - $start;
		# first key is the hmmsearch evalue, second key is the orthoid
		push( @{ $result->{$$line[0]}->{$$line[1]} }, {
			'hmmhit'       => $$line[2],
			'header'       => $$line[3],
			'sequence'     => substr($$line[4], $start, $length),
			'hmm_start'    => $$line[5],
			'hmm_end'      => $$line[6],
			'env_start'    => $$line[7],
			'env_end'      => $$line[8],
			'blast_hit'    => $$line[9],
			'blast_evalue' => $$line[10],
			'reftaxon'     => $$line[11],
		});
	}
	$sth->finish();
	$dbh->disconnect();
	scalar keys %$result > 0 ? return $result : return undef;
}

sub get_results_for_score {
	my $setid   = shift;
	my $taxid   = shift;
	my $min     = shift;
	my $max     = shift;
	# generic query
	my $query = "SELECT DISTINCT $mysql_table_hmmsearch.$mysql_col_score,
			$mysql_table_orthologs.$mysql_col_orthoid,
			$mysql_table_hmmsearch.$mysql_col_target,
			$mysql_table_ests.$mysql_col_header,
			$mysql_table_ests.$mysql_col_sequence,
			$mysql_table_hmmsearch.$mysql_col_hmm_start,
			$mysql_table_hmmsearch.$mysql_col_hmm_end,
			$mysql_table_hmmsearch.$mysql_col_env_start,
			$mysql_table_hmmsearch.$mysql_col_env_end,
			$mysql_table_blast.$mysql_col_target,
			$mysql_table_blast.$mysql_col_evalue,
			$mysql_table_taxa.$mysql_col_name
		FROM $mysql_table_scores
		LEFT JOIN $mysql_table_hmmsearch
			ON $mysql_table_scores.$mysql_col_score = $mysql_table_hmmsearch.$mysql_col_score
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
		WHERE $mysql_table_hmmsearch.$mysql_col_score      IS NOT NULL
			AND $mysql_table_ests.$mysql_col_digest          IS NOT NULL
			AND $mysql_table_orthologs.$mysql_col_orthoid    IS NOT NULL
			AND $mysql_table_blast.$mysql_col_query          IS NOT NULL
			AND $mysql_table_aaseqs.$mysql_col_id            IS NOT NULL
			AND $mysql_table_taxa.$mysql_col_id              IS NOT NULL
			AND $mysql_table_set_details.$mysql_col_id       IS NOT NULL
			AND $mysql_table_set_details.$mysql_col_id       = ?
			AND $mysql_table_hmmsearch.$mysql_col_taxid      = ?";

	# modify the generic query
	# score range
	if ($max) { $query .= "\n			AND $mysql_table_hmmsearch.$mysql_col_score BETWEEN ? AND ?" }
	# single score
	else      { $query .= "\n			AND $mysql_table_hmmsearch.$mysql_col_score = ?" }

	# good for debugging
	print $query . "\n" if $debug;

	my $dbh = get_dbh()
		or return undef;
	my $sth = $dbh->prepare($query);

	# score range
	if ($max) {
		$sth = execute($sth, $mysql_timeout, $setid, $taxid, $min, $max);
	}
	# single score
	else      {
		$sth = execute($sth, $mysql_timeout, $setid, $taxid, $min);
	} 

	# will hold the result
	my $result = { };

	while (my $line = $sth->fetchrow_arrayref()) {
		my $start = $$line[7] - 1;
		my $length = $$line[8] - $start;
		# first key is the hmmsearch score, second key is the orthoid
		push( @{ $result->{$$line[0]}->{$$line[1]} }, {
			'hmmhit'       => $$line[2],
			'header'       => $$line[3],
			'sequence'     => substr($$line[4], $start, $length),
			'hmm_start'    => $$line[5],
			'hmm_end'      => $$line[6],
			'env_start'    => $$line[7],
			'env_end'      => $$line[8],
			'blast_hit'    => $$line[9],
			'blast_evalue' => $$line[10],
			'reftaxon'     => $$line[11],
		});
	}
	$sth->finish();
	$dbh->disconnect();
	scalar keys %$result > 0 ? return $result : return undef;
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
	
=head2 get_reference_sequence(scalar int ID)

Fetches the amino acid sequence for ID from the database. Returns a string.

=cut

sub get_reference_sequence {
	my $id = shift @_ or croak "Usage: get_reference_sequence(ID)\n";
	my $query = "SELECT $mysql_col_sequence 
		FROM $mysql_table_aaseqs
		WHERE $mysql_col_id = '$id'";
	my $result = &mysql_get($query);
	return $result->[0]->[0];
}

sub get_transcript_for {
	my $digest = shift @_ or croak "Usage: get_transcript_for(ID)\n";
	my $query  = "SELECT $mysql_col_sequence
		FROM $mysql_table_ests
		WHERE $mysql_col_digest = ?";
	my $result = &mysql_get($query, $digest);
	return $result->[0]->[0];
}

sub get_nucleotide_transcript_for {
	my $digest = shift @_ or croak "Usage: get_transcript_for(ID)\n";
	my $query  = "SELECT $mysql_col_header
		FROM $mysql_table_ests
		WHERE $mysql_col_digest = ?";
	my $result = &mysql_get($query, $digest);
	# remove the revcomp/translate portion
	print "translated header: <$result->[0]->[0]>\n" if $debug;
	(my $original_header = $result->[0]->[0]) =~ s/ ?(\[revcomp]:)?\[translate\(\d\)\]$//;
	print "original header: <$original_header>\n" if $debug;
	$query = "SELECT $mysql_col_sequence
		FROM $mysql_table_ests
		WHERE $mysql_col_header = ?";
	$result = &mysql_get($query, $original_header);
	return $result->[0]->[0];
}

=head2 get_nuc_for_pep(scalar int ID)

Fetches the nucleotide sequence for a given amino acid sequence with id ID from the database. Returns a string.

=cut

sub get_nuc_for_pep {
	my $pepid = shift @_ or croak "Usage: get_nuc_for_pep(PEPTIDE_ID)\n";
	my $query = "SELECT $mysql_table_seqpairs.$mysql_col_ntseq 
		FROM $mysql_table_seqpairs
		WHERE $mysql_table_seqpairs.$mysql_col_aaseq = ?";
	print $query, "\n", $pepid, "\n";
	my $dbh = get_dbh()
		or return undef;
	my $sth = $dbh->prepare($query);
	$sth = execute($sth, $mysql_timeout, $pepid);
	my $data = $sth->fetchall_arrayref();
	print Dumper($data); exit;
}

=head2 get_real_table_names(int ID, string EST_TABLE, string HMMSEARCH_TABLE, string BLAST_TABLE)

Renames the table names according to ID. Returns a list of the three table names.

=cut

sub get_real_table_names {
	my $specid = shift @_;
	my $real_table_ests      = $mysql_table_ests      . '_' . $specid;
	my $real_table_hmmsearch = $mysql_table_hmmsearch . '_' . $specid;
	my $real_table_blast     = $mysql_table_blast     . '_' . $specid;
	$mysql_table_ests        = $real_table_ests;
	$mysql_table_hmmsearch   = $real_table_hmmsearch;
	$mysql_table_blast       = $real_table_blast;
	return ($real_table_ests, $real_table_hmmsearch, $real_table_blast);
}

=head2 get_scores_list

Returns list of scores as present in the scores view

=cut

sub get_scores_list {
	my $q = "SELECT `score` FROM $mysql_table_scores ORDER BY `$mysql_table_scores`.`$mysql_col_score` DESC";
	return map { $_->[0] } @{mysql_get($q)};
}

=head2 get_hmmresult_for_score(SCORE)

Gets a list of hmmsearch hits for a given score

Arguments: scalar float 

Returns: arrayref of arrayrefs

  [
   [
    query,
    target,
    log_evalue,
    env_start,
    env_end,
    hmm_start,
    hmm_end
   ],
   [
    ...
   ]
  ]

=cut

sub get_hmmresult_for_score {
	my $score = shift;
	my $q_score_row = "SELECT 
		$mysql_table_hmmsearch.$mysql_col_query,
		$mysql_table_hmmsearch.$mysql_col_target,
		$mysql_table_hmmsearch.$mysql_col_score,
		$mysql_table_hmmsearch.$mysql_col_log_evalue,
		$mysql_table_hmmsearch.$mysql_col_env_start,
		$mysql_table_hmmsearch.$mysql_col_env_end,
		$mysql_table_hmmsearch.$mysql_col_hmm_start,
		$mysql_table_hmmsearch.$mysql_col_hmm_end
		FROM $mysql_table_hmmsearch
		WHERE $mysql_table_hmmsearch.$mysql_col_score = ?
		ORDER BY $mysql_table_hmmsearch.$mysql_col_log_evalue";
	my $d = mysql_get($q_score_row, $score);
	my $r = [];
	foreach (@$d) {
		push @$r, {
			'query'      => $_->[0],
			'target'     => $_->[1],
			'score'      => $_->[2],
			'log_evalue' => $_->[3],
			'env_start'  => $_->[4],
			'env_end'    => $_->[5],
			'hmm_start'  => $_->[6],
			'hmm_end'    => $_->[7],
		}
	}
	return $r;
}

sub get_blastresult_for_digest {
	my $digest = shift;
	my $q_blastresult = "SELECT
		$mysql_table_blast.$mysql_col_query,
		$mysql_table_blast.$mysql_col_target,
		$mysql_table_blast.$mysql_col_score,
		$mysql_table_blast.$mysql_col_log_evalue,
		$mysql_table_blast.$mysql_col_start,
		$mysql_table_blast.$mysql_col_end
		FROM $mysql_table_blast
		WHERE $mysql_table_blast.$mysql_col_query = ?
		ORDER BY $mysql_table_blast.$mysql_col_score";
	my $d = mysql_get($q_blastresult, $digest);
	my $r = [];
	foreach (@$d) {
		push @$r, {
			'query'      => $_->[0],
			'target'     => $_->[1],
			'score'      => $_->[2],
			'log_evalue' => $_->[3],
			'start'      => $_->[4],
			'end'        => $_->[5],
		}
	}
	return $r;
}

sub get_real_header {
	my $digest = shift;
	my $q = "SELECT $mysql_table_ests.$mysql_col_header
		FROM $mysql_table_ests
		WHERE $mysql_table_ests.$mysql_col_digest = ?
		LIMIT 1";
	print $q, "\n" if $debug;
	my $d = mysql_get($q, $digest);
	return $d->[0]->[0];
}

1;
