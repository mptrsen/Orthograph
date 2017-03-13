#--------------------------------------------------
# This file is part of Orthograph.
# Copyright 2014 Malte Petersen <mptrsen@uni-bonn.de>
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

Wrapper::Sqlite

=head1 SYNOPSIS

  use Wrapper::Sqlite;

  my $dbh = Wrapper::Sqlite::db_dbh();
  do_stuff_with_dbh();
  undef $dbh;

  Wrapper::Sqlite::do($query);

  my $result = Wrapper::Sqlite::get($query);

=head1 DESCRIPTION

Wrapper module that provides db functions. Deals with specific queries to
handle the complex database structure in the Orthograph pipeline.

=cut

package Wrapper::Sqlite;
use strict;
use warnings;
use Carp;
use Exporter;
use FindBin;        # locate the dir of this script during compile time
use lib $FindBin::RealBin;             # $RealBin is the directory of the original script
use Orthograph::Config;                # configuration parser getconfig()
use Orthograph::Functions;
use Data::Dumper;
use DBI;
use DBD::SQLite;

my $config = $Orthograph::Config::config;  # copy config

# db settings
my $database                = $config->{'sqlite-database'};
my $db_timeout              = 600;
my $sqlite                  = $config->{'sqlite-program'};
my $sleep_for               = 1;

my $db_attached             = 'species_database';
my $db_table_aaseqs         = $config->{'db_table_aaseqs'};
my $db_table_blast          = $config->{'db_table_blast'};
my $db_table_blastdbs       = $config->{'db_table_blastdbs'};
my $db_table_ests           = $config->{'db_table_ests'};
my $db_table_hmmsearch      = $config->{'db_table_hmmsearch'};
my $db_table_log_evalues    = $config->{'db_table_log_evalues'};
my $db_table_scores         = $config->{'db_table_scores'};
my $db_table_ntseqs         = $config->{'db_table_ntseqs'};
my $db_table_ogs            = $config->{'db_table_ogs'};
my $db_table_orthologs      = $config->{'db_table_orthologs'};
my $db_table_seqpairs       = $config->{'db_table_sequence_pairs'};
my $db_table_seqtypes       = $config->{'db_table_sequence_types'};
my $db_table_set_details    = $config->{'db_table_set_details'};
my $db_table_species_info   = $config->{'db_table_species_info'};
my $db_table_taxa           = $config->{'db_table_taxa'};
my $db_table_temp           = $config->{'db_table_temp'};
my $db_col_aaseq            = 'aa_seq';
my $db_col_ali_end          = 'ali_end';
my $db_col_ali_start        = 'ali_start';
my $db_col_blastdb          = 'blastdb';
my $db_col_blastdb_path     = 'blastdb_path';
my $db_col_core             = 'core';
my $db_col_date             = 'date';
my $db_col_digest           = 'digest';
my $db_col_description      = 'description';
my $db_col_end              = 'end';
my $db_col_env_end          = 'env_end';
my $db_col_env_start        = 'env_start';
my $db_col_evalue           = 'evalue';
my $db_col_hmm_end          = 'hmm_end';
my $db_col_hmm_start        = 'hmm_start';
my $db_col_hmmsearch_id     = 'hmmsearch_id';
my $db_col_header           = 'header';
my $db_col_id               = 'id';
my $db_col_log_evalue       = 'log_evalue';
my $db_col_longname         = 'longname';
my $db_col_score            = 'score';
my $db_col_name             = 'name';
my $db_col_ntseq            = 'nt_seq';
my $db_col_ogsid            = 'ogs_id';
my $db_col_ogsversion       = 'version';
my $db_col_orthoid          = 'ortholog_gene_id';
my $db_col_orthoset         = 'orthoset';
my $db_col_query            = 'query';
my $db_col_rebuild          = 'rebuild';
my $db_col_setid            = 'setid';
my $db_col_sequence         = 'sequence';
my $db_col_seqpair          = 'sequence_pair';
my $db_col_start            = 'start';
my $db_col_target           = 'target';
my $db_col_taxid            = 'taxid';
my $db_col_type             = 'type';
my $db_col_version          = 'version';
my $outdir                  = $config->{'output-directory'};
my $orthoset                = $config->{'ortholog-set'};
my $quiet                   = $config->{'quiet'};
my $reftaxa                 = $config->{'reference-taxa'};
# substitution character for selenocysteine, which normally leads to blast freaking out
my $u_subst                 = $config->{'substitute-u-with'};
my $sets_dir                = $config->{'sets-dir'};
# species name; remove all whitespace
(my $species_name = $config->{'species-name'}) =~ s/\s/_/g;
my $g_species_id            = undef;	# global variable
my $verbose                 = $config->{'verbose'};
my $debug                   = $config->{'debug'};
my $stdout = *STDOUT;
my $stderr = *STDERR;
my $attached_db_file        = File::Spec->catfile($config->{'output-directory'}, $species_name . '.sqlite');
my $query_attach_file       = "ATTACH DATABASE '$attached_db_file' as '$db_attached'";

# report that this module is loaded
print "Using SQLite database file '$database'\n" if $verbose;
print "Using SQLite database file '$attached_db_file' as species-specific database\n" if $verbose;
unless (-f $attached_db_file) { Orthograph::Functions::touch($attached_db_file) }

# test whether the sqlite binary exists where specified
Orthograph::Functions::program_exists($sqlite) or print "Fatal: SQLite program not executable where specified at '$sqlite'. Verify path and/or permissions\n" and exit(1);

=head1 FUNCTIONS

=head2 pass_stderr

Reassign STDERR to a different filehandle

=cut

sub pass_stderr {
	$stderr = shift;
}


=head2 pass_stdout

Reassign STDOUT to a different filehandle

=cut

sub pass_stdout {
	$stdout = shift;
}


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
	my $args = shift @_;
	my $dbh = undef;
	my $slept = 0;

	until ($dbh = DBI->connect("DBI:SQLite:$database")) {
		if ($slept >= $db_timeout) { 
			carp "Warning: Connection retry timeout exceeded\n" and return undef;
		}
		carp "Warning: Connection failed, retrying in $sleep_for seconds\n";
		sleep $sleep_for;
		$slept += $sleep_for;
	}

	if ($dbh) {
		$dbh->sqlite_busy_timeout($db_timeout * 1000);
		if ($debug > 1) { print $query_attach_file, "\n" }
		$dbh->do($query_attach_file) or die "Fatal: Could not ATTACH DATABASE: $DBI::errstr";
		return $dbh;
	}
	fail_and_exit("Could not get a database connection: $DBI::errstr");
}#}}}

=head2 db_get($query)

Get from the database the result of a SQL query

Expects: QUERY as a string literal

Returns: Reference to array of arrays (result lines->fields)

=cut

sub db_get {#{{{
	my $query = shift;
	unless ($query) { croak "Usage: db_get(QUERY, ARGS)\n" }
	my @args = @_;
  # prepare anonymous array
	my $results = [ ];

	print Dumper( $query, @args) if $debug > 1;

  # connect and fetch stuff
	my $dbh = get_dbh()
		or return undef;
	my $sth = $dbh->prepare($query) or croak;
	$sth = execute($sth, $db_timeout, @args);
	while (my @result = $sth->fetchrow_array() ) {
		push(@$results, \@result);
	}
	$sth->finish();
	$dbh->disconnect; # disconnect ASAP
	return $results;
}#}}}

=head2 db_do($query)

Connect to a database, execute a single $query (for repetitive queries, you
better do that by hand for performance reasons).

Expects: scalar string SQL query. 

Returns 1 on result, dies otherwise.

=cut

sub db_do {#{{{
	my $query = shift;
	unless ($query) { croak "Usage: db_do(QUERY, ARG, ARG, ...)\n" }
	my @args = @_;
	my $dbh = get_dbh()
		or return undef;
	print $query, "\n" if $debug;
	my $sth = $dbh->prepare($query) or return 0;
	$sth = execute($sth, $db_timeout, @args);
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
	my $sth = $dbh->prepare($query) or fail_and_exit("Could not prepare query: $DBI::errstr");
	$sth->execute(@_);
	if ($sth->fetchrow_array()) {
		return 1;
	}
	return 0;
}#}}}

sub execute {
	my $sth = shift or croak "Usage: execute(STH, TIMEOUT, ARGS)\n";
	my $timeout = shift or croak "Usage: execute(STH, TIMEOUT, ARGS)\n";
	my @args = @_;
	my $slept = 0;
	until ($sth->execute(@args)) {
		warn "Warning: execution failed ($DBI::errstr), retrying in $sleep_for seconds...\n";
		if ($slept > $timeout) { croak "Fatal: execution ultimately failed, failing this transaction\n" }
		sleep $sleep_for;
		$slept += $sleep_for;
	}
	return $sth;
}

sub attached_db_file {
	return $attached_db_file;
}

sub drop_tables {
	my $t = shift @_;
	print 'DROPing tables: ', join(", ", values(%$t)), "\n" if $verbose;
	my $dbh = get_dbh() or fail_and_exit("Couldn't get database connection");
	foreach my $table (keys(%$t)) {
		$dbh->do("DROP TABLE IF EXISTS $t->{$table}") or die "Could not execute drop query: $!\n";
	}
	$dbh->disconnect();
}

sub create_tables {
	my $t = shift @_;
	# the queries for the individual tables
	my %create_table = (#{{{
		# table: blastdbs
		'blastdbs' => "CREATE TABLE `$db_table_blastdbs` (
			`$db_col_id`           INTEGER PRIMARY KEY,
			`$db_col_setid`        INTEGER UNSIGNED DEFAULT NULL UNIQUE,
			`$db_col_rebuild`      INT(1))",
		
		# table: ogs
		'ogs' => "CREATE TABLE `$db_table_ogs` (
			`$db_col_id`           INTEGER PRIMARY KEY,
			`$db_col_type`         INT(1),
			`$db_col_taxid`        INTEGER UNSIGNED NOT NULL,
			`$db_col_version`      TEXT(255),
			UNIQUE ($db_col_type, $db_col_taxid, $db_col_version))",
		
		# table: ortholog_set
		'ortholog_set' => "CREATE TABLE `$db_table_orthologs` (
			`$db_col_id`               INTEGER PRIMARY KEY,
			`$db_col_setid`            INTEGER UNSIGNED NOT NULL,
			`$db_col_orthoid`          TEXT(10) NOT NULL,
			`$db_col_seqpair`          INTEGER UNSIGNED NOT NULL,
			UNIQUE ($db_col_setid, $db_col_orthoid, $db_col_seqpair))",

		# table: sequence_pairs
		'sequence_pairs' => "CREATE TABLE `$db_table_seqpairs` (
			`$db_col_id`           INTEGER PRIMARY KEY,
			`$db_col_taxid`        INTEGER UNSIGNED,
			`$db_col_ogsid`        INTEGER UNSIGNED,
			`$db_col_aaseq`        INTEGER UNSIGNED UNIQUE DEFAULT NULL,
			`$db_col_ntseq`        INTEGER UNSIGNED UNIQUE DEFAULT NULL, 
			`$db_col_date`         INTEGER UNSIGNED DEFAULT CURRENT_TIMESTAMP)",

		# table: sequences_aa
		'aa_sequences' => "CREATE TABLE `$db_table_aaseqs` (
			`$db_col_id`           INTEGER PRIMARY KEY,
			`$db_col_taxid`        INTEGER     NOT NULL, 
			`$db_col_ogsid`        INTEGER     NOT NULL,
			`$db_col_header`       TEXT(512)   UNIQUE,
			`$db_col_sequence`     MEDIUMBLOB,
			`$db_col_date`         INTEGER     UNSIGNED DEFAULT CURRENT_TIMESTAMP)",

		# table: sequences_nt
		'nt_sequences' => "CREATE TABLE `$db_table_ntseqs` (
			`$db_col_id`           INTEGER     PRIMARY KEY,
			`$db_col_taxid`        INTEGER     NOT NULL, 
			`$db_col_ogsid`        INTEGER     NOT NULL,
			`$db_col_header`       TEXT(512)   UNIQUE,
			`$db_col_sequence`     MEDIUMBLOB,
			`$db_col_date`         INTEGER     UNSIGNED DEFAULT CURRENT_TIMESTAMP)",

		# table: set_details
		'set_details' => "CREATE TABLE `$db_table_set_details` (
			`$db_col_id`           INTEGER PRIMARY KEY,
			`$db_col_name`         TEXT(255) UNIQUE,
			`$db_col_description`  BLOB)",

		# table: taxa
		'taxa' => "CREATE TABLE `$db_table_taxa` (
			`$db_col_id`           INTEGER PRIMARY KEY,
			`$db_col_name`         TEXT(20)  UNIQUE,
			`$db_col_core`         TINYINTEGER UNSIGNED NOT NULL)",

		# table: seqtypes
		'seqtypes' => "CREATE TABLE `$db_table_seqtypes` (
			`$db_col_id`           INTEGER PRIMARY KEY,
			`$db_col_type`         TEXT(3)     UNIQUE)",
	);#}}}

	my @indices = (
		# indices for sequences_aa
"CREATE INDEX IF NOT EXISTS $t->{'aaseqs'}_$db_col_taxid  ON $t->{'aaseqs'} ($db_col_taxid)",
"CREATE INDEX IF NOT EXISTS $t->{'ntseqs'}_$db_col_taxid  ON $t->{'ntseqs'} ($db_col_taxid)",
"CREATE INDEX IF NOT EXISTS $t->{'aaseqs'}_$db_col_header  ON $t->{'aaseqs'} ($db_col_header)",
"CREATE INDEX IF NOT EXISTS $t->{'ntseqs'}_$db_col_header  ON $t->{'ntseqs'} ($db_col_header)",
"CREATE INDEX IF NOT EXISTS $t->{'seqpairs'}_$db_col_aaseq  ON $t->{'seqpairs'} ($db_col_aaseq)",
"CREATE INDEX IF NOT EXISTS $t->{'seqpairs'}_$db_col_ntseq  ON $t->{'seqpairs'} ($db_col_ntseq)",
	);

	# useful pragmas for performance?
	my @pragmas = (
		"PRAGMA main.page_size=4096",
		"PRAGMA main.cache_size=10000",
		"PRAGMA main.synchronous=NORMAL",
		"PRAGMA main.journal_mode=WAL",
	);
	# to start off with nt and aa sequence types
	my $insert_seqtypes = "INSERT OR IGNORE INTO $t->{'seqtypes'} (type) VALUES ('nt'),('aa')";

	my $dbh = get_dbh();
	foreach (values %create_table, @indices) {
		print $_, ";\n" if $debug;
		$dbh->do($_) or die "Could not exec query: $DBI::errstr\n";
	}
	$dbh->do($insert_seqtypes);	# start off with 'nt' and 'aa' seqtypes
	$dbh->disconnect;
}


sub create_temp_table {
	my $temptable = shift @_;
	my $create_temp_table_query = "CREATE TABLE $temptable (
			`$db_col_name`        TEXT(255),
			`$db_col_longname`    TEXT(255),
			`$db_col_orthoset`    TEXT(255),
			`$db_col_orthoid`     TEXT(255),
			`$db_col_blastdb`     TEXT(255),
			`$db_col_header`      TEXT(512),
			`$db_col_sequence`    MEDIUMBLOB,
			`$db_col_description` TEXT(255))";
	my $create_temp_indices_query = "BEGIN;
CREATE INDEX IF NOT EXISTS ${temptable}_name ON $temptable ($db_col_name);
CREATE INDEX IF NOT EXISTS ${temptable}_orthoset ON $temptable ($db_col_orthoset);
CREATE INDEX IF NOT EXISTS ${temptable}_orthoid ON $temptable ($db_col_orthoid);
CREATE INDEX IF NOT EXISTS ${temptable}_header ON $temptable ($db_col_header);
COMMIT;";
	my $dbh = get_dbh();
	$dbh->do("DROP TABLE IF EXISTS $temptable") or die "Fatal: Could not DROP TABLE $temptable\n";
	$dbh->do($create_temp_table_query) or die "Fatal: Could not CREATE TABLE $temptable\n";
	$dbh->do($create_temp_indices_query) or die "Fatal: Could not create indices on temporary table\n";
	$dbh->disconnect;
}

sub load_csv_into_temptable {
	my $csvfile   = shift @_;
	my $temptable = shift @_;
	my @loadqueries = (
		".separator ,",
		".mode csv",
		".import $csvfile $temptable",
		".mode list",
	);
	foreach (@loadqueries) {
		if ($debug > 1) {
			print $_, "\n";
			print "execute? "; 
			<STDIN>;
		}
		system qq{$sqlite -separator "," $database "$_"} and die "Fatal: Could not import CSV file '$csvfile' into temporary table $temptable\n";
	}
}

sub fill_tables_from_temp_table {
	my $t = shift @_;
	my $temptable = shift @_;
	my @queries = (
		# taxa (name, longname)
		"INSERT OR IGNORE INTO $t->{'taxa'} ($db_col_name, $db_col_longname, $db_col_core) 
			SELECT DISTINCT $temptable.name, $temptable.$db_col_longname, 1 
			FROM $temptable",
		# set name + description
		"INSERT OR IGNORE INTO $t->{'set_details'} ($db_col_name, $db_col_description)
			SELECT DISTINCT $temptable.orthoset, $temptable.$db_col_description 
			FROM $temptable LIMIT 1",
		# blast databases
		"INSERT OR IGNORE INTO $t->{'blastdbs'} ($db_col_setid, blastdb_path) 
			SELECT DISTINCT $t->{'set_details'}.id, $temptable.$db_col_blastdb 
			FROM $temptable
			LEFT JOIN $t->{'set_details'} 
				ON $t->{'set_details'}.name = $temptable.$db_col_orthoset",
		# pep sequences
		"INSERT OR IGNORE INTO $t->{'aaseqs'} ($db_col_taxid, $db_col_header, $db_col_sequence, $db_col_date) 
			SELECT $t->{'taxa'}.id, $temptable.$db_col_header, $temptable.$db_col_sequence, CURRENT_TIMESTAMP
			FROM $temptable
				LEFT JOIN $t->{'taxa'} 
			ON $temptable.$db_col_name  = $t->{'taxa'}.$db_col_name",
		# delete everything where header or sequence is NULL or empty
		"DELETE FROM $t->{'aaseqs'}
			WHERE $t->{'aaseqs'}.$db_col_header IS NULL
			OR $t->{'aaseqs'}.$db_col_sequence IS NULL
			OR $t->{'aaseqs'}.$db_col_header = ''
			OR $t->{'aaseqs'}.$db_col_sequence = ''",
		# sequence pairs (pep-nuc)
		"INSERT OR IGNORE INTO $t->{'seqpairs'} ($db_col_taxid, $db_col_ogsid, $db_col_aaseq, $db_col_ntseq, $db_col_date)
			SELECT $t->{'taxa'}.id, $t->{'ogs'}.id, $t->{'aaseqs'}.id, $t->{'ntseqs'}.id, CURRENT_TIMESTAMP
			FROM $t->{'taxa'}
			INNER JOIN $t->{'aaseqs'}
				ON $t->{'aaseqs'}.taxid = $t->{'taxa'}.$db_col_id
			LEFT JOIN $t->{'ogs'}
				ON $t->{'taxa'}.id = $t->{'ogs'}.$db_col_taxid
			LEFT JOIN $t->{'ntseqs'}
				ON $t->{'aaseqs'}.header = $t->{'ntseqs'}.$db_col_header",
		# orthologous groups
		"INSERT OR IGNORE INTO $t->{'orthologs'} ($db_col_setid, $db_col_orthoid, $db_col_seqpair) 
			SELECT $t->{'set_details'}.id, $temptable.$db_col_orthoid, $t->{'seqpairs'}.$db_col_id 
			FROM $t->{'aaseqs'} 
			INNER JOIN $temptable 
				ON $t->{'aaseqs'}.$db_col_header = $temptable.$db_col_header 
			INNER JOIN $t->{'seqpairs'} 
				ON $t->{'seqpairs'}.$db_col_aaseq = $t->{'aaseqs'}.$db_col_id 
			INNER JOIN $t->{'set_details'} 
				ON $t->{'set_details'}.$db_col_name = $temptable.$db_col_orthoset",
	);

	my $dbh = get_dbh();
	my $nrows;
	foreach (@queries) {
		print $_ . ";\n" if $debug;
		$nrows = $dbh->do($_) or fail_and_exit("Query failed: $_");
		if ($debug) {
			($nrows > 0) ? printf("Query OK, %d rows affected\n", $nrows) : print "Query OK\n";
		}
	}
	$dbh->disconnect;
	return $nrows;
}

sub get_number_of_cogs_for_set {
	my $setn = shift @_;
	my $q = "
		SELECT COUNT(DISTINCT $db_table_orthologs.ortholog_gene_id)
		FROM $db_table_orthologs
		INNER JOIN $db_table_set_details
			ON $db_table_orthologs.setid = $db_table_set_details.id
		WHERE $db_table_set_details.$db_col_id = ?";
	my $r = db_get($q, $setn);
	return $$r[0][0];
}

=head2 get_ortholog_sets()

Get list of ortholog sets from the database

Arguments: none

Returns: hash reference of set names => description

=cut

sub get_ortholog_sets {#{{{
	my %sets = ();
	my $query = "
		SELECT
			$db_table_set_details.$db_col_name,
			$db_table_set_details.$db_col_description,
			COUNT($db_table_orthologs.$db_col_id),
			COUNT(distinct $db_table_orthologs.$db_col_orthoid)
		FROM $db_table_set_details
		INNER JOIN $db_table_orthologs
			ON $db_table_set_details.$db_col_id = $db_table_orthologs.$db_col_setid
		GROUP BY $db_table_set_details.$db_col_id";
	my $data = db_get($query);
	return($data);
}#}}}


=head2 get_ortholog_groups_for_set($setid)

Returns a hashref of hashrefs to create an ortholog set from. Each key in the hashref (the ortholog group ID) is a hashref of sequence_ID => sequence.

=cut

sub get_ortholog_groups_for_set {
	my $setid = shift @_ or croak "Usage: Wrapper::Sqlite::get_ortholog_groups_for_set(SETID)";
	my $data = {};
	my $query = "SELECT o.ortholog_gene_id, a.header, a.sequence
		FROM $db_table_orthologs         AS o
    INNER JOIN $db_table_seqpairs    AS p
    ON o.sequence_pair = p.id
    INNER JOIN $db_table_aaseqs      AS a
    ON a.id = p.aa_seq
    INNER JOIN $db_table_set_details AS s
    ON s.id = o.setid
    WHERE s.id = ?";

	my $dbh = get_dbh()
		or return undef;
	my $sth = $dbh->prepare($query);
	$sth = execute($sth, $db_timeout, $setid);
	while (my @row = $sth->fetchrow_array()) {
		# load the whole set into memory, i don't give a frak
		$$data{$row[0]}{$row[1]} = $row[2];
	}
	$dbh->disconnect();	# disc asap

	return $data;
}

sub preparedb {
	my $query_create_ests = "CREATE TABLE $db_attached.$db_table_ests ( 
		'$db_col_id'        INTEGER NOT NULL PRIMARY KEY,
		'$db_col_digest'    TEXT(32)     NOT NULL,           
		'$db_col_taxid'     UNSIGNED INTEGER NOT NULL,       
		'$db_col_type'      UNSIGNED TINYINT(4) NOT NULL,
		'$db_col_date'      UNSIGNED INT,
		'$db_col_header'    TEXT      NOT NULL,       
		'$db_col_sequence'  MEDIUMBLOB DEFAULT NULL
		)";

	my $query_create_hmmsearch = "CREATE TABLE $db_attached.$db_table_hmmsearch (
		'$db_col_id'         INTEGER NOT NULL PRIMARY KEY,
		'$db_col_taxid'      UNSIGNED INTEGER NOT NULL,       
		'$db_col_query'      TEXT(255) NOT NULL,       
		'$db_col_target'     TEXT(32)     NOT NULL,       
		'$db_col_score'      DOUBLE       NOT NULL,
		'$db_col_evalue'     TEXT(8)      NOT NULL,
		'$db_col_log_evalue' DOUBLE       NOT NULL DEFAULT '-999',
		'$db_col_env_start'  UNSIGNED INTEGER NOT NULL,
		'$db_col_env_end'    UNSIGNED INTEGER NOT NULL,
		'$db_col_ali_start'  UNSIGNED INTEGER NOT NULL,
		'$db_col_ali_end'    UNSIGNED INTEGER NOT NULL,
		'$db_col_hmm_start'  UNSIGNED INTEGER NOT NULL,
		'$db_col_hmm_end'    UNSIGNED INTEGER NOT NULL
		)";

	my $query_create_blast = "CREATE TABLE $db_attached.$db_table_blast (
		'$db_col_id'            INTEGER NOT NULL PRIMARY KEY,
		'$db_col_taxid'         UNSIGNED INTEGER NOT NULL,       
		'$db_col_query'         TEXT(32)     NOT NULL,       
		'$db_col_target'        UNSIGNED INTEGER NOT NULL,       
		'$db_col_score'         DOUBLE       NOT NULL,
		'$db_col_evalue'        TEXT(8)      NOT NULL,
		'$db_col_log_evalue'    DOUBLE       NOT NULL DEFAULT '-999',
		'$db_col_start'         UNSIGNED INTEGER NOT NULL,
		'$db_col_end'           UNSIGNED INTEGER NOT NULL,
		'$db_col_hmmsearch_id'  UNSIGNED INTEGER NOT NULL
		)";

	my $query_create_species_info = "CREATE TABLE $db_attached.$db_table_species_info (
		'$db_col_id'           INTEGER PRIMARY KEY,
		'$db_col_name'         TEXT(255)
		)",

	# the CREATE INDEX statement does not like the attached db prefix 
	# generate table names without it
	(my $simple_table_ests = $db_table_ests)           =~ s/$db_attached\.//;
	(my $simple_table_hmmsearch = $db_table_hmmsearch) =~ s/$db_attached\.//;
	(my $simple_table_blast = $db_table_blast)         =~ s/$db_attached\.//;
	my @query_create_indices = (
		"CREATE INDEX IF NOT EXISTS $db_attached.${db_table_ests}_header ON $simple_table_ests ($db_col_header)",
		"CREATE INDEX IF NOT EXISTS $db_attached.${db_table_ests}_digest ON $simple_table_ests ($db_col_digest)",
		"CREATE INDEX IF NOT EXISTS $db_attached.${db_table_ests}_taxid ON $simple_table_ests ($db_col_taxid)",
		"CREATE INDEX IF NOT EXISTS $db_attached.${db_table_ests}_header ON $simple_table_ests ($db_col_header)",
		"CREATE INDEX IF NOT EXISTS $db_attached.${db_table_hmmsearch}_taxid ON $simple_table_hmmsearch ($db_col_taxid)",
		"CREATE INDEX IF NOT EXISTS $db_attached.${db_table_hmmsearch}_query ON $simple_table_hmmsearch ($db_col_query)",
		"CREATE INDEX IF NOT EXISTS $db_attached.${db_table_hmmsearch}_target ON $simple_table_hmmsearch ($db_col_target)",
		"CREATE INDEX IF NOT EXISTS $db_attached.${db_table_hmmsearch}_evalue ON $simple_table_hmmsearch ($db_col_log_evalue)",
		"CREATE INDEX IF NOT EXISTS $db_attached.${db_table_hmmsearch}_score ON $simple_table_hmmsearch ($db_col_score)",
		"CREATE INDEX IF NOT EXISTS $db_attached.${db_table_blast}_taxid ON $simple_table_blast ($db_col_taxid)",
		"CREATE INDEX IF NOT EXISTS $db_attached.${db_table_blast}_query ON $simple_table_blast ($db_col_query)",
		"CREATE INDEX IF NOT EXISTS $db_attached.${db_table_blast}_target ON $simple_table_blast ($db_col_target)",
		"CREATE INDEX IF NOT EXISTS $db_attached.${db_table_blast}_evalue ON $simple_table_blast ($db_col_log_evalue)",
		"CREATE INDEX IF NOT EXISTS $db_attached.${db_table_blast}_hmmsearch_id ON $simple_table_blast ($db_col_hmmsearch_id)",
	);

	# drop all tables (just delete the database file)
	unlink $attached_db_file or fail_and_exit("Could not delete database file '$attached_db_file': $!");

	# open connection
	my $dbh = get_dbh()
		or croak "Fatal: Could not connect to database: $DBI::errstr\n" and exit 1;

	# create all tables
	foreach my $query ($query_create_ests, $query_create_hmmsearch, $query_create_blast, $query_create_species_info, @query_create_indices) {
		print "$query;\n" if $debug;
		my $sql = $dbh->prepare($query);
		$sql->execute()
		  or croak "Fatal: Could not execute SQL query: $DBI::errstr\n" and exit(1);
	}

	# disconnect
	$dbh->disconnect();
}

sub get_transcripts_for_species {
	my $specid = shift or croak "Usage: Wrapper::Sqlite::get_transcripts(SPECIESID, TYPE)";
	my $type = shift or croak "Usage: Wrapper::Sqlite::get_transcripts(SPECIESID, TYPE)";
	my $query = "SELECT digest, sequence
		FROM $db_table_ests
		WHERE taxid = ?
		AND type = ?";
	my $dbh = get_dbh()
		or return undef;
	my $sth = $dbh->prepare($query);
	$sth = execute($sth, $db_timeout, $specid, $type);
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
	my ($hmmquery, $taxid) = @_ or croak "Usage: Wrapper::Sqlite::get_hmmresults(HMMQUERY)";
	# disable query cache for this one
	my $query_get_sequences = "SELECT $db_table_ests.digest,
		  $db_table_ests.sequence,
		  $db_table_hmmsearch.env_start,
		  $db_table_hmmsearch.env_end,
			$db_table_hmmsearch.id
		FROM $db_table_ests 
		INNER JOIN $db_table_hmmsearch
		ON $db_table_hmmsearch.target = $db_table_ests.digest
		WHERE $db_table_hmmsearch.query = ?
		AND $db_table_hmmsearch.taxid = ?";

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
	# TODO rewrite this part using parametrized queries to protect from SQL injections?
	my $query = "SELECT DISTINCT $db_table_set_details.name, $db_table_taxa.name
		FROM $db_table_seqpairs
		INNER JOIN $db_table_taxa
			ON $db_table_seqpairs.taxid = $db_table_taxa.id
		INNER JOIN $db_table_orthologs
			ON $db_table_orthologs.sequence_pair = $db_table_seqpairs.id 
		INNER JOIN $db_table_set_details
			ON $db_table_orthologs.setid = $db_table_set_details.id"
	;
	my $data = db_get($query) or croak();
	return $data;
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
	my $query = "SELECT DISTINCT $db_table_set_details.name, $db_table_taxa.name
		FROM $db_table_seqpairs
		INNER JOIN $db_table_taxa
			ON $db_table_seqpairs.taxid = $db_table_taxa.id
		INNER JOIN $db_table_orthologs
			ON $db_table_orthologs.sequence_pair = $db_table_seqpairs.id 
		INNER JOIN $db_table_set_details
			ON $db_table_orthologs.setid = $db_table_set_details.id
		WHERE $db_table_set_details.id = '$set_id'"
	;
	my $data = db_get($query);
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
	my $result = db_get("SELECT COUNT(*) FROM $db_table_ests WHERE $db_col_type = 2");

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
	my $taxids = db_get("SELECT id FROM $db_table_taxa WHERE name IN ($taxa_string)");

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
	my $aaseqs = db_get("SELECT COUNT(*) FROM  $db_table_aaseqs WHERE $db_table_aaseqs.taxid IN ($taxids_string)");

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
	my $aaseqs = db_get("SELECT $db_table_aaseqs.id, $db_table_aaseqs.sequence FROM  $db_table_aaseqs WHERE $db_table_aaseqs.taxid IN ($taxids_string)");

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
	my $query = "SELECT $db_table_species_info.$db_col_id FROM $db_table_species_info WHERE $db_col_name = '$species_name'";
	my $result = db_get($query);
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
	my $query = "SELECT $db_col_id FROM $db_table_set_details WHERE name = ?";
	my $result = db_get($query,  $setname);
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
	my $sth = $dbh->prepare("SELECT * FROM $db_table_set_details WHERE $db_col_name = ? LIMIT 1");
	$sth = execute($sth, $db_timeout, $set);
	my $result = $sth->fetchrow_arrayref;
	if ( $$result[0] ) { return 1 }
	return 0;
}

=head2 insert_taxon_into_database(TAXON_NAME)

Inserts a core taxon into the database.

Returns the newly generated taxon ID.

=cut

sub insert_taxon_into_database {
	my $name = shift;
	my $core = shift;
	db_do("INSERT OR IGNORE INTO $db_table_taxa ($db_col_name, $db_col_core) VALUES (?, ?)", $name, $core) or croak;
	my $res = db_get("SELECT $db_col_id FROM $db_table_taxa WHERE $db_col_name = ? AND $db_col_core = ?", $name, $core);
	return $res->[0]->[0];
}

=head2 insert_species_info(TAXON_NAME)

Inserts a (non-core) taxon into the database. The taxon shorthand will be NULL
and the 'core' switch will be 0.

Returns the newly generated taxon ID.

=cut

sub insert_species_info {
	my $species_name = shift(@_);
	unless ($species_name) { croak("Usage: Wrapper::Mysql::insert_species_info(SPECIESNAME)") }
	if (my $taxid = get_taxid_for_species($species_name)) { return $taxid }
	my $query = "INSERT OR IGNORE INTO $db_table_species_info ($db_col_name) VALUES (?)";
	db_do($query, $species_name);

	$g_species_id = get_taxid_for_species($species_name) or croak;
	return $g_species_id;
}

sub insert_ogs_info_into_database {
	my $type       = shift;
	my $taxid      = shift;
	my $ogsversion = shift;
	db_do("INSERT OR IGNORE INTO $db_table_ogs ($db_col_type, $db_col_taxid, $db_col_version) VALUES (?, ?, ?)", $type, $taxid, $ogsversion) or croak;
	my $res = db_get("SELECT $db_col_id FROM $db_table_ogs WHERE $db_col_taxid = ? AND $db_col_version = ? AND $db_col_type = ?", $taxid, $ogsversion, $type);
	return $res->[0]->[0];
}

sub upload_ogs_sequences {
	my ($inf, $hdrs, $taxid, $type, $ogsid) = @_;

	# determine correct table and columns
	my $seqtable      = $type > 1 ? $db_table_aaseqs : $db_table_ntseqs;
	my $otherseqtable = $type > 1 ? $db_table_ntseqs : $db_table_aaseqs;
	my $seqcol        = $type > 1 ? $db_col_aaseq    : $db_col_ntseq;
	my $otherseqcol   = $type > 1 ? $db_col_ntseq    : $db_col_aaseq;

	load_csv_into_temptable($inf, $db_table_temp);

	my $query_insert_sequences = "
		INSERT OR IGNORE INTO $seqtable ($db_col_taxid, $db_col_ogsid, $db_col_header, $db_col_sequence, $db_col_date)
		SELECT $db_table_taxa.$db_col_id, $db_table_temp.$db_col_ogsid, $db_table_temp.$db_col_header, $db_table_temp.$db_col_sequence, CURRENT_TIMESTAMP
		FROM $db_table_temp 
		LEFT JOIN $db_table_taxa 
			ON $db_table_temp.$db_col_taxid = $db_table_taxa.$db_col_id
	";

	my $query_insert_seqpairs = "
		INSERT OR IGNORE INTO $db_table_seqpairs ($db_col_taxid, $db_col_ogsid, $seqcol, $otherseqcol, $db_col_date)
		SELECT $db_table_taxa.$db_col_id, $seqtable.$db_col_ogsid, $seqtable.$db_col_id, $otherseqtable.$db_col_id, CURRENT_TIMESTAMP
		FROM $seqtable
		INNER JOIN $db_table_taxa
			ON $seqtable.$db_col_taxid = $db_table_taxa.$db_col_id
		LEFT JOIN $otherseqtable
			ON $seqtable.$db_col_header = $otherseqtable.$db_col_header
		WHERE $seqtable.$db_col_ogsid = $ogsid
	";

	my $dbh = get_dbh();
	# if this is amino acid data, just blindly load it into the resp. table
	if ($type > 1) {
		$dbh->do($query_insert_sequences) or fail_and_exit("OGS loading failed: $DBI::errstr");
		$dbh->do($query_insert_seqpairs)  or fail_and_exit("OGS loading failed: $DBI::errstr");
	}
	# otherwise, this is nucleotide data, need to check for each sequence
	# whether the corresponding aa seq exists
	else {
		upload_sequences_individually($hdrs, $seqtable, $otherseqtable, $seqcol, $otherseqcol)
	}
	$dbh->disconnect();
	return 1;
}

sub upload_sequences_individually {
	my $hdrs          = shift;
	my $seqtable      = shift;
	my $otherseqtable = shift;
	my $seqcol        = shift;
	my $otherseqcol   = shift;
	my $dbh = get_dbh();
	my $sth_select = $dbh->prepare("SELECT $db_col_id FROM $otherseqtable WHERE $db_col_header = ?");
	my $query_insert = "
		INSERT OR IGNORE INTO $seqtable ($db_col_taxid, $db_col_ogsid, $db_col_header, $db_col_sequence, $db_col_date)
		SELECT $db_table_taxa.$db_col_id, $db_table_temp.$db_col_ogsid, $db_table_temp.$db_col_header, $db_table_temp.$db_col_sequence, CURRENT_TIMESTAMP
		FROM $db_table_temp 
		LEFT JOIN $db_table_taxa 
			ON $db_table_temp.$db_col_taxid = $db_table_taxa.$db_col_id
		WHERE $db_col_header = ?
	";
	my $sth_insert = $dbh->prepare($query_insert);
	my $query_update = "
		UPDATE $db_table_seqpairs
		SET $seqcol = (SELECT $db_col_id FROM $seqtable WHERE $db_col_header = ? LIMIT 1)
		WHERE $otherseqcol = ?
	";
	my $sth_update = $dbh->prepare($query_update);

	my $res = [ ];
	my $n   = scalar @$hdrs;
	my $c   = 0;
	foreach my $hdr (@$hdrs) {
		# check whether the hdr exists for an aa seq
		$sth_select->execute($hdr);
		$res = $sth_select->fetchall_arrayref();
		if (scalar @$res == 1) {
			# ok, upload
			$sth_insert->execute($hdr);
			# and update seqpairs
			$sth_update->execute($hdr, $$res[0][0]);
			$c++;
			Orthograph::Functions::progress_bar($c, $n, 25, '-');
		}
		else {
			fail_and_exit("Corresponding amino acid sequence not found for nucleotide sequence with ID '$hdr'")
		}
	}
	return 1;
}

sub insert_orthologs {
	my $query = "
		INSERT INTO $db_table_orthologs ($db_col_setid, $db_col_orthoid, $db_col_seqpair)
		SELECT
			$db_table_temp.$db_col_setid,
			$db_table_temp.$db_col_orthoid,
			$db_table_seqpairs.$db_col_id
		FROM $db_table_temp
		INNER JOIN $db_table_aaseqs
			ON $db_table_temp.$db_col_header = $db_table_aaseqs.$db_col_header
		INNER JOIN $db_table_seqpairs 
			ON $db_table_aaseqs.$db_col_id = $db_table_seqpairs.$db_col_aaseq
	";
	db_do($query) or return 0;
	return 1;
}

sub get_number_of_orthologs_for_set {
	my $setid = shift;
	my $query = "SELECT COUNT($db_col_id) FROM $db_table_orthologs WHERE $db_col_setid = ?";
	my $res = db_get($query, $setid);
	return $res->[0]->[0];
}

sub load_ests_from_file {
	my $csvfile = shift;
	my $list = shift;
	my $specid = shift;

	# load data from csv file into database
	# create temporary table first

	my $q_drop_temp   = "DROP TABLE IF EXISTS $db_table_temp";
	my $q_create_temp = "CREATE TABLE $db_table_temp (
	  '$db_col_digest' TEXT,
		'$db_col_taxid'  INT,
		'$db_col_type'   INT,
		'$db_col_date'   INT,
		'$db_col_header' TEXT,
		'$db_col_sequence' TEXT)
	";

	# load data into temptable
	# needs to work directly on the species database since apparently
	# sqlite cannot use .import on attached databases. if anyone has an idea 
	# or a solution on how to mitigate this, please let me know!
	my @loadqueries = (
		".mode csv",
		$q_drop_temp,
		$q_create_temp,
		".import $csvfile $db_table_temp",
		".mode list",
	);
	foreach (@loadqueries) {
		my $cmd = qq{$sqlite -separator "," "$attached_db_file" "$_"};
		if ($debug) {
			print "$cmd\n";
			print "execute? "; 
			<STDIN>;
		}
		system("$cmd") and die "Fatal: Could not import CSV file '$csvfile' into temporary table $db_table_temp\n";
	}

	# transfer data from temptable into main table
	my $q_transfer = "INSERT INTO $db_attached.$db_table_ests (
		$db_col_digest,
  	$db_col_taxid,
  	$db_col_type,
  	$db_col_date,
  	$db_col_header,
  	$db_col_sequence)
  	SELECT 
    $db_col_digest,
  	$db_col_taxid,
  	$db_col_type,
  	$db_col_date,
  	$db_col_header,
  	$db_col_sequence
		FROM $db_attached.$db_table_temp
	";
	my $dbh = get_dbh();
	print $q_transfer, "\n" if $debug > 1;
	my $sth = $dbh->prepare($q_transfer);
	my $num_ests = $sth->execute();
	$sth->finish();
	$dbh->do("DROP TABLE $db_table_temp");
	$dbh->disconnect;
	if (defined($DBI::errstr)) { print "$DBI::errstr\n" and exit(1) }
	return $num_ests;
}


sub create_log_evalues_view {
	unless (scalar @_ == 1) { croak 'Usage: Wrapper::Sqlite::create_log_evalues_view($species_id)' }
	my $taxid = shift;
	my $query_drop_log_evalues = "DROP VIEW IF EXISTS $db_table_log_evalues";
	my $query_create_log_evalues = "CREATE VIEW $db_table_log_evalues AS
	  SELECT $db_table_hmmsearch.$db_col_log_evalue AS $db_col_log_evalue,
	    COUNT($db_table_hmmsearch.$db_col_log_evalue) AS `count`
	  FROM $db_table_hmmsearch
	  WHERE $db_table_hmmsearch.$db_col_taxid = $taxid
	  GROUP BY $db_table_hmmsearch.$db_col_log_evalue
	  ORDER BY $db_table_hmmsearch.$db_col_log_evalue";
	my $dbh = get_dbh()
		or return undef;
	$dbh->do($query_drop_log_evalues);
	my $sth = $dbh->prepare($query_create_log_evalues);
	$sth = execute($sth, $db_timeout);
	$dbh->disconnect();
	return 1;
}

sub create_scores_view {
	unless (scalar @_ == 1) { croak 'Usage: Wrapper::Sqlite::create_scores_view($species_id)' }
	my $taxid = shift;
	my $query_drop_scores_view = "DROP VIEW IF EXISTS $db_table_scores";
	my $query_create_scores_view = "CREATE VIEW IF NOT EXISTS $db_table_scores AS
	  SELECT $db_table_hmmsearch.$db_col_score AS $db_col_score,
	    COUNT($db_table_hmmsearch.$db_col_score) AS `count`
	  FROM $db_table_hmmsearch
	  WHERE $db_table_hmmsearch.$db_col_taxid = $taxid
	  GROUP BY $db_table_hmmsearch.$db_col_score
	  ORDER BY $db_table_hmmsearch.$db_col_score DESC";
	my $dbh = get_dbh()
		or return undef;
	$dbh->do($query_drop_scores_view);
	$dbh->do($query_create_scores_view);
	$dbh->disconnect();
	return 1;
}
	


# get a orthoid => list_of_aaseq_ids relationship from the db
sub get_orthologs_for_set_hashref {
	my $setid = shift(@_);
	unless ($setid) { croak("Usage: get_orthologs_for_set(SETID)") }
	my $query = "SELECT DISTINCT
		$db_table_orthologs.ortholog_gene_id,
		$db_table_aaseqs.id 
		FROM $db_table_orthologs 
		INNER JOIN $db_table_seqpairs 
			ON $db_table_orthologs.sequence_pair = $db_table_seqpairs.id
		INNER JOIN $db_table_aaseqs
			ON $db_table_seqpairs.aa_seq = $db_table_aaseqs.id
		INNER JOIN $db_table_set_details 
			ON $db_table_orthologs.setid = $db_table_set_details.id
		WHERE $db_table_set_details.id = ?";
	my $dbh = get_dbh()
		or return undef;
	my $sth = $dbh->prepare($query);
	$sth = execute($sth, $db_timeout, $setid);
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
			$db_table_taxa.$db_col_name,
			$db_table_aaseqs.$db_col_header,
			$db_table_aaseqs.$db_col_sequence
		FROM $db_table_aaseqs
		INNER JOIN $db_table_seqpairs
			ON $db_table_aaseqs.$db_col_id = $db_table_seqpairs.$db_col_aaseq
		INNER JOIN $db_table_orthologs
			ON $db_table_seqpairs.$db_col_id = $db_table_orthologs.$db_col_seqpair
		INNER JOIN $db_table_taxa
			ON $db_table_aaseqs.$db_col_taxid = $db_table_taxa.$db_col_id
		AND   $db_table_orthologs.$db_col_setid = ?
		AND   $db_table_orthologs.$db_col_orthoid = ?
		ORDER BY $db_table_taxa.$db_col_name";
	my $dbh = get_dbh()
		or return undef;
	my $sth = $dbh->prepare($query);
	$sth = execute($sth, $db_timeout, $setid, $orthoid);
	my $data = $sth->fetchall_arrayref();
	return $data;
}

sub get_ortholog_group_nucleotide {
	my $setid   = shift;
	my $orthoid = shift;
	my $query = "SELECT 
			$db_table_taxa.$db_col_name,
			$db_table_ntseqs.$db_col_header,
			$db_table_ntseqs.$db_col_sequence
		FROM $db_table_ntseqs
		INNER JOIN $db_table_seqpairs
			ON $db_table_ntseqs.$db_col_id = $db_table_seqpairs.$db_col_ntseq
		INNER JOIN $db_table_orthologs
			ON $db_table_seqpairs.$db_col_id = $db_table_orthologs.$db_col_seqpair
		INNER JOIN $db_table_taxa
			ON $db_table_ntseqs.$db_col_taxid = $db_table_taxa.$db_col_id
		AND   $db_table_orthologs.$db_col_setid = ?
		AND   $db_table_orthologs.$db_col_orthoid = ?
		ORDER BY $db_table_taxa.$db_col_name";
	my $dbh = get_dbh()
		or return undef;
	my $sth = $dbh->prepare($query);
	$sth = execute($sth, $db_timeout, $setid, $orthoid);
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
		$db_table_hmmsearch.evalue,
		$db_table_orthologs.ortholog_gene_id, 
		$db_table_hmmsearch.target,
		$db_table_ests.header,
		$db_table_ests.sequence,
		$db_table_hmmsearch.hmm_start,
		$db_table_hmmsearch.hmm_end,
		$db_table_blast.target,
		$db_table_blast.evalue,
		$db_table_taxa.name
		FROM $db_table_hmmsearch
		INNER JOIN $db_table_ests
			ON $db_table_hmmsearch.target = $db_table_ests.digest
		INNER JOIN $db_table_orthologs
			ON $db_table_hmmsearch.query = $db_table_orthologs.ortholog_gene_id
		INNER JOIN $db_table_blast
			ON $db_table_hmmsearch.target = $db_table_blast.query
		INNER JOIN $db_table_aaseqs
			ON $db_table_blast.target = $db_table_aaseqs.id
		INNER JOIN  $db_table_taxa
			ON $db_table_aaseqs.taxid = $db_table_taxa.id
		INNER JOIN $db_table_set_details
			ON $db_table_orthologs.setid = $db_table_set_details.id
		WHERE $db_table_set_details.id = ?
		AND $db_table_hmmsearch.taxid  = ?
		ORDER BY $db_table_hmmsearch.log_evalue ASC
		LIMIT $limit 
		OFFSET $offset
		";
	print "fetching:\n$query\n";
	my $dbh = get_dbh()
		or return undef;
	my $sth = $dbh->prepare($query);
	$sth = execute($sth, $db_timeout, $setid, $specid);
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
	my $query_get_logevalues = "SELECT $db_col_log_evalue, count FROM $db_table_log_evalues";
	my $dbh = get_dbh()
		or return undef;
	my $sth = $dbh->prepare($query_get_logevalues);
	$sth = execute($sth, $db_timeout);
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
	my $query_get_scores = "SELECT $db_col_score, count FROM $db_table_scores";
	my $dbh = get_dbh()
		or return undef;
	my $sth = $dbh->prepare($query_get_scores);
	$sth = execute($sth, $db_timeout);
	my $d = $sth->fetchall_arrayref();
	$sth->finish();
	$dbh->disconnect();
	my $num_of_scores = { };
	foreach my $row (@$d) {
		$num_of_scores->{$$row[0]} = $$row[1];
	}
	return $num_of_scores;
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
	my $query = "SELECT DISTINCT $db_table_hmmsearch.$db_col_evalue,
			$db_table_orthologs.$db_col_orthoid,
			$db_table_hmmsearch.$db_col_target,
			$db_table_ests.$db_col_header,
			$db_table_ests.$db_col_sequence,
			$db_table_hmmsearch.$db_col_hmm_start,
			$db_table_hmmsearch.$db_col_hmm_end,
			$db_table_hmmsearch.$db_col_env_start,
			$db_table_hmmsearch.$db_col_env_end,
			$db_table_blast.$db_col_target,
			$db_table_blast.$db_col_evalue,
			$db_table_taxa.$db_col_name
		FROM $db_table_log_evalues
		LEFT JOIN $db_table_hmmsearch
			ON $db_table_log_evalues.$db_col_log_evalue = $db_table_hmmsearch.$db_col_log_evalue
		LEFT JOIN $db_table_ests
			ON $db_table_hmmsearch.$db_col_target = $db_table_ests.$db_col_digest
		LEFT JOIN $db_table_orthologs
			ON $db_table_hmmsearch.$db_col_query = $db_table_orthologs.$db_col_orthoid
		LEFT JOIN $db_table_blast
			ON $db_table_hmmsearch.$db_col_target = $db_table_blast.$db_col_query
		LEFT JOIN $db_table_aaseqs
			ON $db_table_blast.$db_col_target = $db_table_aaseqs.$db_col_id
		LEFT JOIN $db_table_taxa
			ON $db_table_aaseqs.$db_col_taxid = $db_table_taxa.$db_col_id
		LEFT JOIN $db_table_set_details
			ON $db_table_orthologs.$db_col_setid = $db_table_set_details.$db_col_id
		WHERE $db_table_hmmsearch.$db_col_log_evalue IS NOT NULL
			AND $db_table_ests.$db_col_digest          IS NOT NULL
			AND $db_table_orthologs.$db_col_orthoid    IS NOT NULL
			AND $db_table_blast.$db_col_query          IS NOT NULL
			AND $db_table_aaseqs.$db_col_id            IS NOT NULL
			AND $db_table_taxa.$db_col_id              IS NOT NULL
			AND $db_table_set_details.$db_col_id       IS NOT NULL
			AND $db_table_set_details.$db_col_id       = ?
			AND $db_table_hmmsearch.$db_col_taxid      = ?";

	# modify the generic query
	# e-value range
	if ($max) { $query .= "\n			AND $db_table_hmmsearch.$db_col_log_evalue BETWEEN ? AND ?" }
	# single e-value
	else      { $query .= "\n			AND $db_table_hmmsearch.$db_col_log_evalue = ?" }

	# good for debugging
	print $query . "\n" if $debug > 1;

	my $dbh = get_dbh()
		or return undef;
	my $sth = $dbh->prepare($query);

	# e-value range
	if ($max) {
		$sth = execute($sth, $db_timeout, $setid, $taxid, $min, $max);
	}
	# single e-value
	else      {
		$sth = execute($sth, $db_timeout, $setid, $taxid, $min);
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

sub get_hmmresults_for_single_score {
	my $setid   = shift;
	my $score   = shift;
	# generic query
	my $query = "SELECT DISTINCT
			$db_table_orthologs.$db_col_orthoid,
			$db_table_hmmsearch.$db_col_id,
			$db_table_hmmsearch.$db_col_target,
			$db_table_hmmsearch.$db_col_evalue,
			$db_table_hmmsearch.$db_col_hmm_start,
			$db_table_hmmsearch.$db_col_hmm_end,
			$db_table_hmmsearch.$db_col_ali_start,
			$db_table_hmmsearch.$db_col_ali_end,
			$db_table_hmmsearch.$db_col_env_start,
			$db_table_hmmsearch.$db_col_env_end,
			$db_table_ests.$db_col_header
		FROM $db_table_hmmsearch
		LEFT JOIN $db_table_ests
			ON $db_table_hmmsearch.$db_col_target = $db_table_ests.$db_col_digest
		LEFT JOIN $db_table_orthologs
			ON $db_table_hmmsearch.$db_col_query = $db_table_orthologs.$db_col_orthoid
		LEFT JOIN $db_table_species_info
			ON $db_table_hmmsearch.$db_col_taxid = $db_table_species_info.$db_col_id
		LEFT JOIN $db_table_set_details
			ON $db_table_orthologs.$db_col_setid = $db_table_set_details.$db_col_id
		WHERE $db_table_ests.$db_col_digest          IS NOT NULL
			AND $db_table_orthologs.$db_col_orthoid    IS NOT NULL
			AND $db_table_species_info.$db_col_id              IS NOT NULL
			AND $db_table_set_details.$db_col_id       IS NOT NULL
			AND $db_table_set_details.$db_col_id       = ?
			AND $db_table_hmmsearch.$db_col_score      = ?
	";

	# good for debugging
	print $query . "\n" if $debug > 1;

	my $dbh = get_dbh()
		or return undef;
	my $sth = $dbh->prepare($query);

	# single score
	$sth = execute($sth, $db_timeout, $setid, $score);

	# will hold the result
	my $result = [ ];

	while (my $line = $sth->fetchrow_arrayref()) {
		# first key is the hmmsearch score, second key is the orthoid
		push( @$result, {
			'orthoid'      => $$line[0],
			'hmmsearch_id' => $$line[1],
			'hmmhit'       => $$line[2],
			'hmm_evalue'   => $$line[3],
			'hmm_start'    => $$line[4],
			'hmm_end'      => $$line[5],
			'ali_start'    => $$line[6],
			'ali_end'      => $$line[7],
			'env_start'    => $$line[8],
			'env_end'      => $$line[9],
			'header'       => $$line[10],
		});
	}
	$sth->finish();
	$dbh->disconnect();
	return $result;
}

sub get_blastresults_for_hmmsearch_id {
	my $setid          = shift;
	my $hmmsearch_id   = shift;
	# generic query
	my $query = "SELECT DISTINCT
			$db_table_blast.$db_col_target,
			$db_table_blast.$db_col_score,
			$db_table_blast.$db_col_evalue,
			$db_table_blast.$db_col_start,
			$db_table_blast.$db_col_end,
			$db_table_hmmsearch.$db_col_target,
			$db_table_hmmsearch.$db_col_evalue,
			$db_table_hmmsearch.$db_col_env_start,
			$db_table_hmmsearch.$db_col_env_end,
			$db_table_hmmsearch.$db_col_ali_start,
			$db_table_hmmsearch.$db_col_ali_end,
			$db_table_hmmsearch.$db_col_hmm_start,
			$db_table_hmmsearch.$db_col_hmm_end,
			$db_table_ests.$db_col_header
		FROM $db_table_hmmsearch
		LEFT JOIN $db_table_blast
			ON $db_table_hmmsearch.$db_col_id = $db_table_blast.$db_col_hmmsearch_id
		LEFT JOIN $db_table_ests
			ON $db_table_ests.$db_col_digest = $db_table_hmmsearch.$db_col_target
		WHERE $db_table_hmmsearch.$db_col_id         IS NOT NULL
			AND $db_table_blast.$db_col_hmmsearch_id   IS NOT NULL
			AND $db_table_hmmsearch.$db_col_target     IS NOT NULL
			AND $db_table_hmmsearch.$db_col_id         = ?
		ORDER BY $db_table_blast.$db_col_score DESC
	";

	# good for debugging
	if ($debug > 1) {
		print $query . "\n";
		print "Executing this query with $hmmsearch_id\n";
	}

	my $dbh = get_dbh()
		or return undef;
	my $sth = $dbh->prepare($query);

	# single score
	$sth = execute($sth, $db_timeout, $hmmsearch_id);

	# will hold the result
	my $result = [ ];

	while (my $line = $sth->fetchrow_arrayref()) {
		# first key is the hmmsearch score, second key is the orthoid
		push( @$result, {
			'blast_hit'    => $$line[0],
			'blast_score'  => $$line[1],
			'blast_evalue' => $$line[2],
			'blast_start'  => $$line[3],
			'blast_end'    => $$line[4],
			'hmmhit'       => $$line[5],
			'hmm_evalue'   => $$line[6],
			'env_start'    => $$line[7],
			'env_end'      => $$line[8],
			'ali_start'    => $$line[9],
			'ali_end'      => $$line[10],
			'hmm_start'    => $$line[11],
			'hmm_end'      => $$line[12],
			'header'       => $$line[13],
		});
	}
	$sth->finish();
	$dbh->disconnect();
	return $result;
}


sub get_results_for_single_score {
	my $setid   = shift;
	my $score   = shift;
	# generic query
	my $query = "SELECT DISTINCT
			$db_table_orthologs.$db_col_orthoid,
			$db_table_hmmsearch.$db_col_target,
			$db_table_hmmsearch.$db_col_evalue,
			$db_table_hmmsearch.$db_col_hmm_start,
			$db_table_hmmsearch.$db_col_hmm_end,
			$db_table_hmmsearch.$db_col_env_start,
			$db_table_hmmsearch.$db_col_env_end,
			$db_table_blast.$db_col_target,
			$db_table_blast.$db_col_score,
			$db_table_blast.$db_col_evalue,
			$db_table_blast.$db_col_start,
			$db_table_blast.$db_col_end,
			$db_table_ests.$db_col_header
		FROM $db_table_hmmsearch
		LEFT JOIN $db_table_ests
			ON $db_table_hmmsearch.$db_col_target = $db_table_ests.$db_col_digest
		LEFT JOIN $db_table_orthologs
			ON $db_table_hmmsearch.$db_col_query = $db_table_orthologs.$db_col_orthoid
		LEFT JOIN $db_table_blast
			ON $db_table_hmmsearch.$db_col_target = $db_table_blast.$db_col_query
		LEFT JOIN $db_table_aaseqs
			ON $db_table_blast.$db_col_target = $db_table_aaseqs.$db_col_id
		LEFT JOIN $db_table_taxa
			ON $db_table_aaseqs.$db_col_taxid = $db_table_taxa.$db_col_id
		LEFT JOIN $db_table_set_details
			ON $db_table_orthologs.$db_col_setid = $db_table_set_details.$db_col_id
		WHERE $db_table_ests.$db_col_digest          IS NOT NULL
			AND $db_table_orthologs.$db_col_orthoid    IS NOT NULL
			AND $db_table_blast.$db_col_query          IS NOT NULL
			AND $db_table_aaseqs.$db_col_id            IS NOT NULL
			AND $db_table_taxa.$db_col_id              IS NOT NULL
			AND $db_table_set_details.$db_col_id       IS NOT NULL
			AND $db_table_set_details.$db_col_id       = ?
			AND $db_table_hmmsearch.$db_col_score      = ?
		ORDER BY $db_table_blast.$db_col_score DESC
	";

	# good for debugging
	print $query . "\n" if $debug > 1;

	my $dbh = get_dbh()
		or return undef;
	my $sth = $dbh->prepare($query);

	# single score
	$sth = execute($sth, $db_timeout, $setid, $score);

	# will hold the result
	my $result = [ ];

	while (my $line = $sth->fetchrow_arrayref()) {
		# first key is the hmmsearch score, second key is the orthoid
		push( @$result, {
			'orthoid'      => $$line[0],
			'hmmhit'       => $$line[1],
			'hmm_evalue'   => $$line[2],
			'hmm_start'    => $$line[3],
			'hmm_end'      => $$line[4],
			'env_start'    => $$line[5],
			'env_end'      => $$line[6],
			'blast_hit'    => $$line[7],
			'blast_score'  => $$line[8],
			'blast_evalue' => $$line[9],
			'blast_start'  => $$line[10],
			'blast_end'    => $$line[11],
			'header'       => $$line[12],
		});
	}
	$sth->finish();
	$dbh->disconnect();
	return $result;
}

sub get_results_for_score {
	my $setid   = shift;
	my $taxid   = shift;
	my $min     = shift;
	my $max     = shift;
	# generic query
	my $query = "SELECT DISTINCT $db_table_hmmsearch.$db_col_score,
			$db_table_orthologs.$db_col_orthoid,
			$db_table_hmmsearch.$db_col_target,
			$db_table_ests.$db_col_header,
			$db_table_ests.$db_col_sequence,
			$db_table_hmmsearch.$db_col_hmm_start,
			$db_table_hmmsearch.$db_col_hmm_end,
			$db_table_hmmsearch.$db_col_env_start,
			$db_table_hmmsearch.$db_col_env_end,
			$db_table_blast.$db_col_target,
			$db_table_blast.$db_col_evalue,
			$db_table_taxa.$db_col_name
		FROM $db_table_scores
		LEFT JOIN $db_table_hmmsearch
			ON $db_table_scores.$db_col_score = $db_table_hmmsearch.$db_col_score
		LEFT JOIN $db_table_ests
			ON $db_table_hmmsearch.$db_col_target = $db_table_ests.$db_col_digest
		LEFT JOIN $db_table_orthologs
			ON $db_table_hmmsearch.$db_col_query = $db_table_orthologs.$db_col_orthoid
		LEFT JOIN $db_table_blast
			ON $db_table_hmmsearch.$db_col_target = $db_table_blast.$db_col_query
		LEFT JOIN $db_table_aaseqs
			ON $db_table_blast.$db_col_target = $db_table_aaseqs.$db_col_id
		LEFT JOIN $db_table_taxa
			ON $db_table_aaseqs.$db_col_taxid = $db_table_taxa.$db_col_id
		LEFT JOIN $db_table_set_details
			ON $db_table_orthologs.$db_col_setid = $db_table_set_details.$db_col_id
		WHERE $db_table_hmmsearch.$db_col_score      IS NOT NULL
			AND $db_table_ests.$db_col_digest          IS NOT NULL
			AND $db_table_orthologs.$db_col_orthoid    IS NOT NULL
			AND $db_table_blast.$db_col_query          IS NOT NULL
			AND $db_table_aaseqs.$db_col_id            IS NOT NULL
			AND $db_table_taxa.$db_col_id              IS NOT NULL
			AND $db_table_set_details.$db_col_id       IS NOT NULL
			AND $db_table_set_details.$db_col_id       = ?
			AND $db_table_hmmsearch.$db_col_taxid      = ?";

	# modify the generic query
	# score range
	if ($max) { $query .= "\n			AND $db_table_hmmsearch.$db_col_score BETWEEN ? AND ?" }
	# single score
	else      { $query .= "\n			AND $db_table_hmmsearch.$db_col_score = ?" }

	# good for debugging
	print $query . "\n" if $debug > 1;

	my $dbh = get_dbh()
		or return undef;
	my $sth = $dbh->prepare($query);

	# score range
	if ($max) {
		$sth = execute($sth, $db_timeout, $setid, $taxid, $min, $max);
	}
	# single score
	else      {
		$sth = execute($sth, $db_timeout, $setid, $taxid, $min);
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
		$db_table_hmmsearch.target
		FROM $db_table_hmmsearch
		INNER JOIN $db_table_orthologs
			ON $db_table_hmmsearch.query = $db_table_orthologs.ortholog_gene_id
		INNER JOIN $db_table_blast
			ON $db_table_hmmsearch.target = $db_table_blast.query
		INNER JOIN $db_table_aaseqs
			ON $db_table_blast.target = $db_table_aaseqs.id
		INNER JOIN  $db_table_taxa
			ON $db_table_aaseqs.taxid = $db_table_taxa.id
		INNER JOIN $db_table_set_details
			ON $db_table_orthologs.setid = $db_table_set_details.id
		WHERE $db_table_set_details.id = $setid
		AND $db_table_hmmsearch.taxid  = $specid
	";
	my $data = db_get($query);
	my @result;
	push(@result, ${shift(@$data)}[0]) while @$data;
	return @result;
}
	
=head2 get_reference_sequence(scalar int ID)

Fetches the amino acid sequence for ID from the database. Returns a string.

=cut

sub get_reference_sequence {
	my $id = shift @_ or croak "Usage: get_reference_sequence(ID)\n";
	my $query = "SELECT
			$db_table_aaseqs.$db_col_sequence, 
			$db_table_taxa.$db_col_name
		FROM $db_table_aaseqs
		INNER JOIN $db_table_taxa
			ON $db_table_aaseqs.$db_col_taxid = $db_table_taxa.id
		WHERE $db_table_aaseqs.$db_col_id = '$id'";
	my $result = db_get($query);
	return ($result->[0]->[0], $result->[0]->[1]);
}

sub get_transcript_for {
	my $digest = shift @_ or croak "Usage: get_transcript_for(ID)\n";
	my $query  = "SELECT $db_col_sequence
		FROM $db_table_ests
		WHERE $db_col_digest = ?";
	my $result = db_get($query, $digest);
	return $result->[0]->[0];
}

sub get_nucleotide_transcript_for {
	my $original_header = shift @_ or croak "Usage: get_transcript_for(ID)\n";
	my $query = "SELECT $db_col_sequence
		FROM $db_table_ests
		WHERE $db_col_header = ?";
	my $result = db_get($query, $original_header);
	return ($result->[0]->[0]);
}

=head2 get_nuc_for_pep(scalar int ID)

Fetches the nucleotide sequence for a given amino acid sequence with id ID from the database. Returns a string.

=cut

sub get_nuc_for_pep {
	my $pepid = shift @_ or croak "Usage: get_nuc_for_pep(PEPTIDE_ID)\n";
	my $query = "SELECT $db_table_seqpairs.$db_col_ntseq 
		FROM $db_table_seqpairs
		WHERE $db_table_seqpairs.$db_col_aaseq = ?";
	print $query, "\n", $pepid, "\n";
	my $dbh = get_dbh()
		or return undef;
	my $sth = $dbh->prepare($query);
	$sth = execute($sth, $db_timeout, $pepid);
	my $data = $sth->fetchall_arrayref();
}

=head2 get_real_table_names(int ID, string EST_TABLE, string HMMSEARCH_TABLE, string BLAST_TABLE)

Renames the table names according to ID. Returns a list of the three table names.

=cut

sub get_real_table_names {
	#--------------------------------------------------
	# my $specid = shift @_;
	# my $real_table_ests      = $db_attached . '.' . $db_table_ests         . '_' . $specid;
	# my $real_table_hmmsearch = $db_attached . '.' . $db_table_hmmsearch    . '_' . $specid;
	# my $real_table_blast     = $db_attached . '.' . $db_table_blast        . '_' . $specid;
	# my $real_table_temp      = $db_attached . '.' . $db_table_temp         . '_' . $specid;
	# my $real_table_species   = $db_attached . '.' . $db_table_species_info . '_' . $specid;
	# $db_table_ests           = $real_table_ests;
	# $db_table_hmmsearch      = $real_table_hmmsearch;
	# $db_table_blast          = $real_table_blast;
	# $db_table_temp           = $real_table_temp;
	# $db_table_species_info   = $real_table_temp;
	#-------------------------------------------------- 
	return ($db_table_ests, $db_table_hmmsearch, $db_table_blast, $db_table_temp, $db_table_species_info);
}


=head2 get_scores_list

Returns list of scores as present in the scores view

=cut

sub get_scores_list {
	my $specid = shift;
	my $setid = shift;
	my $threshold = shift;
	my $r = db_get("SELECT score FROM $db_table_hmmsearch
		WHERE score >= $threshold
		GROUP BY score
		ORDER BY score DESC");
	# flatten multidimensional array
	$r = [ map { @$_ } @$r ];
	return $r;
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
		$db_table_hmmsearch.$db_col_query,
		$db_table_hmmsearch.$db_col_target,
		$db_table_hmmsearch.$db_col_score,
		$db_table_hmmsearch.$db_col_log_evalue,
		$db_table_hmmsearch.$db_col_env_start,
		$db_table_hmmsearch.$db_col_env_end,
		$db_table_hmmsearch.$db_col_hmm_start,
		$db_table_hmmsearch.$db_col_hmm_end
		FROM $db_table_hmmsearch
		WHERE $db_table_hmmsearch.$db_col_score = ?
		ORDER BY $db_table_hmmsearch.$db_col_log_evalue";
	my $d = db_get($q_score_row, $score);
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
		$db_table_blast.$db_col_query,
		$db_table_blast.$db_col_target,
		$db_table_blast.$db_col_score,
		$db_table_blast.$db_col_log_evalue,
		$db_table_blast.$db_col_start,
		$db_table_blast.$db_col_end
		FROM $db_table_blast
		WHERE $db_table_blast.$db_col_query = ?
		ORDER BY $db_table_blast.$db_col_score";
	my $d = db_get($q_blastresult, $digest);
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
	my $q = "SELECT $db_table_ests.$db_col_header
		FROM $db_table_ests
		WHERE $db_table_ests.$db_col_digest = ?
		LIMIT 1";
	print $q, "\n" if $debug > 1;
	my $d = db_get($q, $digest);
	return $d->[0]->[0];
}

sub insert_results_into_blast_table {
	my $hits = shift;
	my $species_id = shift;
	my $hmmsearch_id = shift;
	my $hitcount = 0;

	my $query_insert_result = "INSERT OR IGNORE INTO $db_table_blast (
		`$db_col_taxid`,
		`$db_col_query`,
		`$db_col_target`,
		`$db_col_score`,
		`$db_col_evalue`,
		`$db_col_log_evalue`,
		`$db_col_start`,
		`$db_col_end`,
		`$db_col_hmmsearch_id`
		) VALUES (
		?,
		?,
		?,
		?,
		?,
		?,
		?,
		?,
		?
	)";

	my $dbh = get_dbh()
		or print "Fatal: Could not connect to database: $DBI::errstr\n" and exit 1;
	$dbh->do("BEGIN");
	my $sql = $dbh->prepare($query_insert_result);

	# this is a reference to an array of hashes
	foreach my $hit (@$hits) {
		$sql->execute(
			$species_id,
			$hit->{'query'},  # query (HMM)
			$hit->{'target'}, # target (header)
			$hit->{'score'},  # score
			$hit->{'evalue'},  # evalue
			$hit->{'evalue'} != 0 ? log($hit->{'evalue'}) : -999,  # natural logarithm only if not 0
			$hit->{'start'},
			$hit->{'end'},
			$hmmsearch_id,
		) or print "Fatal: Could not push to database!\n" and exit(1);
		++$hitcount;
	}
	$dbh->do("COMMIT");
	$dbh->disconnect;
	return $hitcount;
}

sub insert_results_into_hmmsearch_table {
	my $hits = shift;
	my $species_id = shift;
	my $hitcount = 0;

	# SQL query for pushing HMMsearch results to the db
	my $query_insert_result = "INSERT OR IGNORE INTO $db_table_hmmsearch (
		`$db_col_taxid`,
		`$db_col_query`,
		`$db_col_target`,
		`$db_col_score`,
		`$db_col_evalue`,
		`$db_col_log_evalue`,
		`$db_col_hmm_start`,
		`$db_col_hmm_end`,
		`$db_col_ali_start`,
		`$db_col_ali_end`,
		`$db_col_env_start`,
		`$db_col_env_end`
		) VALUES (
		?,
		?,
		?,
		?,
		?,
		?,
		?,
		?,
		?,
		?,
		?,
		?)";

	my $dbh = get_dbh();
	$dbh->do("BEGIN");
	my $sql = $dbh->prepare($query_insert_result);

	# this is a reference to an array of hashes
	foreach my $hit (@$hits) {
		my $affected_rows = $sql->execute(
			$species_id,
			$hit->{'query'},  # query (HMM)
			$hit->{'target'}, # target (header)
			$hit->{'score'},  # score
			$hit->{'evalue'}, # evalue
			$hit->{'evalue'} != 0 ? log($hit->{'evalue'}) : -999,	# natural logarithm only if evalue not 0
			$hit->{'hmm_start'},  # start of hit domain on the HMM
			$hit->{'hmm_end'},    # end of hit domain on the HMM
			$hit->{'ali_start'},  # start of hit domain on the target seq
			$hit->{'ali_end'},    # end of hit domain on the target seq
			$hit->{'env_start'},  # start of hit domain on the target seq (envelope)
			$hit->{'env_end'},    # end of hit domain on the target seq (envelope)
		) or print "Fatal: Could not push to database!\n" and exit(1);
		if ($affected_rows > 0) { ++$hitcount }
	}
	$dbh->do("COMMIT");
	$dbh->disconnect;
	return $hitcount;
}


sub get_orthologs {
	my ($orthoid, $estdigest) = @_;
	# what do we want from the database?
	# TODO rewrite this part using parametrized queries to protect from SQL injections?
	my $query = "SELECT DISTINCT
			$db_table_orthologs.$db_col_orthoid AS orthogroup,
			$db_table_taxa.$db_col_name	        AS name,
			$db_table_ests.$db_col_digest       AS EST_digest,
      $db_table_aaseqs.$db_col_header     AS AA_header,
      $db_table_aaseqs.$db_col_sequence   AS AA_seq,
      $db_table_ests.$db_col_header       AS EST_hdr,
      $db_table_ests.$db_col_sequence     AS EST_seq
    FROM $db_table_aaseqs
		INNER JOIN $db_table_taxa
			ON $db_table_aaseqs.$db_col_taxid = $db_table_taxa.$db_col_id
    INNER JOIN $db_table_blast 
			ON $db_table_aaseqs.$db_col_id = $db_table_blast.$db_col_target 
    INNER JOIN $db_table_hmmsearch 
			ON $db_table_blast.$db_col_query = $db_table_hmmsearch.$db_col_target  
    INNER JOIN $db_table_ests 
			ON $db_table_hmmsearch.$db_col_target = $db_table_ests.$db_col_digest 
    INNER JOIN $db_table_orthologs 
			ON $db_table_hmmsearch.$db_col_query = $db_table_orthologs.$db_col_orthoid
		INNER JOIN $db_table_seqpairs
			ON $db_table_orthologs.$db_col_seqpair = $db_table_seqpairs.$db_col_id 
		WHERE $db_table_orthologs.$db_col_orthoid = ?
			AND $db_table_ests.$db_col_digest       = ?
			AND $db_table_seqpairs.$db_col_aaseq    = $db_table_aaseqs.$db_col_id
		ORDER BY $db_table_hmmsearch.$db_col_evalue, $db_table_blast.$db_col_evalue";

	# open connection and do the transaction
	my $dbh = get_dbh()
		or print $stderr "Fatal: Could not connect to database: $DBI::errstr\n" and exit 1;
	my $sql = $dbh->prepare($query);
	$sql->execute( $orthoid, $estdigest );
	my $result = $sql->fetchall_arrayref();
	$dbh->disconnect();
	return $result;
}

sub get_taxon_shorthands {
	my $q = "SELECT * FROM $db_table_taxa WHERE `core` = '1'";
	my $dbh = get_dbh();
	my $sth = $dbh->prepare($q);
	$sth->execute();
	my $res = $sth->fetchall_arrayref();
	$dbh->disconnect;
	return $res;
}

sub create_temptable_for_ogs_data {
	# create the temporary table
	my $c = "CREATE TABLE $db_table_temp (
	`$db_col_taxid`    INTEGER      NOT NULL,
	`$db_col_type`     INTEGER      NOT NULL,
	`$db_col_ogsid`    INTEGER      NOT NULL,
	`$db_col_header`   TEXT(255) NOT NULL,
	`$db_col_sequence` TEXT)";
	my $i = "CREATE INDEX temp_header ON $db_table_temp (header)";
	my $dbh = get_dbh();
	print $stdout $c, "\n" if $debug > 1;
	$dbh->do("DROP TABLE IF EXISTS $db_table_temp");
	$dbh->do($c);
	$dbh->do($i);
	$dbh->disconnect();
}

sub import_ogs_into_database {
	my ($f, $hdrs, $seqtable, $otherseqtable, $seqcol, $otherseqcol, $type, $taxid, $ogsversion) = @_;

	load_csv_into_temptable($f, $db_table_temp);

		# insert data into main table. IGNORE is important to avoid duplicates without throwing errors.
	my $query_insert_sequences = "
		INSERT OR IGNORE INTO $seqtable (taxid, header, sequence)
		SELECT $db_table_taxa.id, $db_table_temp.header, $db_table_temp.sequence 
		FROM $db_table_temp 
		LEFT JOIN $db_table_taxa 
			ON $db_table_temp.taxid = $db_table_taxa.id
	";

		# insert sequence pairs relationships. 
		# this is an ugly but required hack because sqlite does not support INSERT
		# OR UPDATE.  we need the original ids for the orthologous groups, though,
		# so REPLACE isn't an option. 
		# I offer a bottle of champagne as well as my eternal gratitude to anyone
		# who can find a cleaner solution.  contact me if you have one!
	my $query_insert_pair =	"
		INSERT OR IGNORE INTO $db_table_seqpairs (
			$db_col_taxid,
			$seqcol,
			$otherseqcol,
			$db_col_date
		)
		SELECT 
			$db_table_taxa.$db_col_id,
			$db_table_aaseqs.$db_col_id,
			$db_table_ntseqs.$db_col_id,
			CURRENT_TIMESTAMP
		FROM $seqtable
		LEFT JOIN $otherseqtable
			ON $seqtable.$db_col_header = $otherseqtable.$db_col_header
		LEFT JOIN $db_table_taxa
			ON $seqtable.$db_col_taxid = $db_table_taxa.$db_col_id
		WHERE $db_table_taxa.$db_col_id = ?
		AND $seqtable.$db_col_header = ?
		UNION 
		SELECT 
			$db_table_taxa.$db_col_id,
			$seqtable.$db_col_id,
			$otherseqtable.$db_col_id,
			CURRENT_TIMESTAMP
		FROM $seqtable
		LEFT JOIN $otherseqtable
			ON $seqtable.$db_col_header = $otherseqtable.$db_col_header
		LEFT JOIN $db_table_taxa
			ON $seqtable.$db_col_taxid = $db_table_taxa.$db_col_id
		WHERE $db_table_taxa.$db_col_id = ?
		AND $otherseqtable.$db_col_header = ?
	";

	my $query_get_pair_ids = "
		SELECT
			$seqtable.$db_col_id,
			$otherseqtable.$db_col_id
		FROM
			$seqtable
		LEFT JOIN $otherseqtable
			ON $seqtable.$db_col_header = $otherseqtable.$db_col_header
		WHERE $seqtable.$db_col_header = ?
			AND $seqtable.$db_col_taxid = ?
	";

	my $query_get_pair_id = "
		SELECT 
			$db_table_seqpairs.$db_col_id,
			$seqtable.$db_col_id,
			$otherseqtable.$db_col_id
		FROM $db_table_seqpairs
		LEFT JOIN $seqtable
			ON $db_table_seqpairs.$seqcol = $seqtable.$db_col_id
		LEFT JOIN $otherseqtable
			ON $db_table_seqpairs.$otherseqcol = $otherseqtable.$db_col_id 
		WHERE ($seqtable.$db_col_header = ?
			OR $otherseqtable.$db_col_header = ?)
			AND $db_table_seqpairs.$db_col_taxid = ?;
	";

	my $query_select_pair = "
		SELECT
			$db_table_seqpairs.$db_col_id
		FROM $db_table_seqpairs
		WHERE $db_table_seqpairs.$seqcol = ?
			OR $db_table_seqpairs.$otherseqcol = ?
	";

	my $query_update_pair = "
		UPDATE OR IGNORE $db_table_seqpairs 
		SET 
			$db_col_taxid = ?,
			$db_col_ogsid = ?,
			$seqcol = ?,
			$otherseqcol = ?,
			$db_col_date = CURRENT_TIMESTAMP
		WHERE
			$db_col_id = ?
	";

	my $query_get_seq_id = "
		SELECT $seqtable.$db_col_id FROM $seqtable WHERE $db_col_header = ? AND $db_col_taxid = ?
	";

	my $query_get_otherseq_id = "
		SELECT $otherseqtable.$db_col_id FROM $otherseqtable WHERE $db_col_header = ? AND $db_col_taxid = ?
	";

	my $query_insert_seqpair = "
		INSERT OR IGNORE INTO $db_table_seqpairs (
			$db_col_taxid,
			$seqcol,
			$otherseqcol,
			$db_col_date
		)
		VALUES (
			?,
			?,
			?,
			CURRENT_TIMESTAMP
		)
	";

	# if # affected > 0 : ok
	# else: 

	my $query_get_seqpair_id = "
		SELECT
			$db_table_seqpairs.$db_col_id
		FROM $db_table_seqpairs
		WHERE $db_table_seqpairs.$seqcol = ?
			OR $db_table_seqpairs.$otherseqcol = ?
	";

	my $query_update_seqpair = "
		UPDATE OR IGNORE $db_table_seqpairs 
		SET 
			$db_col_taxid = ?,
			$db_col_ogsid = ?,
			$seqcol = ?,
			$otherseqcol = ?,
			$db_col_date = CURRENT_TIMESTAMP
		WHERE
			$db_col_id = ?
	";

	my $dbh = get_dbh();
	$dbh->do($query_insert_sequences) or fail_and_exit("OGS loading failed: $DBI::errstr");
	# update OGS table
	my $query_insert_ogs = "INSERT OR IGNORE INTO $db_table_ogs (`type`, `taxid`, `version`) VALUES ('$type', '$taxid', '$ogsversion')";
	if ($debug) {
		print $query_insert_ogs, "\n";
		printf "Execute this with <%s>, <%s>, and <%s>? ", $type, $taxid, $ogsversion;
		quit_on_q();
	}
	$dbh->do($query_insert_ogs) or fail_and_exit("Could not update OGS table: no entry added");
	my $ogsid = $dbh->selectall_arrayref("SELECT $db_col_id FROM $db_table_ogs WHERE $db_col_taxid = $taxid AND $db_col_ogsversion = '$ogsversion' AND $db_col_type = $type");
	$ogsid = $$ogsid[0][0];
	print "Got OGS ID $ogsid for taxon ID $taxid\n" if $debug;

	my $sth_ins  = $dbh->prepare( $query_insert_pair  ) or die;
	my $sth_sel  = $dbh->prepare( $query_get_pair_ids ) or die;
	my $sth_selp = $dbh->prepare( $query_select_pair  ) or die;
	my $sth_upd  = $dbh->prepare( $query_update_pair  ) or die;

	my $sth_get_seq_id      = $dbh->prepare( $query_get_seq_id      ) or die;
	my $sth_get_otherseq_id = $dbh->prepare( $query_get_otherseq_id ) or die;
	my $sth_insert_seqpair  = $dbh->prepare( $query_insert_seqpair  ) or die;
	my $sth_get_seqpair_id  = $dbh->prepare( $query_get_seqpair_id  ) or die;
	my $sth_update_seqpair  = $dbh->prepare( $query_update_seqpair  ) or die;


	# for each header
	#		get the corresponding id from the other seq table, if there is one
	#		try to insert both ids 
	#		test if a row was updated
	#		if not:
	#		update the id
	my $c = 0;
	my $n = scalar @$hdrs;
	foreach my $hdr (@$hdrs) {
		# report progress
		$c++;
		Orthograph::Functions::progress_bar($c, $n, 25, '-');

		# reset variables
		my ($seqid, $otherseqid);

		# get the seq id
		if ($debug) {
			printf "Query: %s\nExecute with <%s> and <%s>?",
				$query_get_seq_id,
				$hdr,
				$taxid
			;
			quit_on_q();
		}
		$sth_get_seq_id->execute($hdr, $taxid);
		my $ids = $sth_get_seq_id->fetchall_arrayref();
		if (scalar @$ids > 1) { fail_and_exit("Found more than one record with ID '$hdr'! Database corrupted?") }
		if (scalar @$ids == 0) {
			$seqid = undef;
		}
		else {
			$seqid = $$ids[0][0];
		}
		if ($debug) {
			printf "got id: <%s>\n", $seqid ? $seqid : 'NULL';
		}

		# get the otherseq id
		if ($debug) {
			printf "Query: %s\nExecute with <%s> and <%s>?",
				$query_get_otherseq_id,
				$hdr,
				$taxid
			;
			quit_on_q();
		}
		$sth_get_otherseq_id->execute($hdr, $taxid);
		$ids = $sth_get_otherseq_id->fetchall_arrayref();
		if (scalar @$ids > 1) { fail_and_exit("Found more than one record with ID '$hdr'! Database corrupted?") }
		if (scalar @$ids == 0) {
			$otherseqid = undef;
		}
		else {
			$otherseqid = $$ids[0][0];
		}
		if ($debug) {
			printf "got id: <%s>\n", defined $otherseqid ? $otherseqid : 'NULL';
		}

		# try to insert new sequence pair
		if ($debug) {
			printf "Query: %s\nExecute with <%s>, <%s> and <%s>?",
			$query_insert_seqpair, 
			$taxid,
			$seqid      ? $seqid      : 'NULL',
			$otherseqid ? $otherseqid : 'NULL'
			;
			quit_on_q();
		}
		$sth_insert_seqpair->execute($taxid, $seqid, $otherseqid);

		# check if any rows were affected (i.e., a new sequence pair was inserted)

		# no rows affected, a pair already present
		if ($sth_insert_seqpair->rows() == 0) { 
			print "no rows affected, pair already present\n" if $debug;

			# get extant pair id
			if ($debug) {
				printf "Query: %s\nExecute with <%s> and <%s>?",
				$query_get_seqpair_id, 
				$seqid      ? $seqid      : '',
				$otherseqid ? $otherseqid : ''
				;
				quit_on_q();
			}
			$sth_get_seqpair_id->execute($seqid, $otherseqid);
			$ids = $sth_get_seqpair_id->fetchall_arrayref();		

			if (scalar @$ids > 1) { fail_and_exit("Found more than one record with ID '$hdr'! Database corrupted?") }
			elsif (scalar @$ids == 0) {
				fail_and_exit("No extant sequence pair with $seqcol ID $seqid or $otherseqcol $otherseqid found! Database corrupted?");
			}

			# there is exactly one sequence pair with these ids present
			# update it with the new values
			else {
				my $seqpairid = $$ids[0][0];
				if ($debug) {
					printf "got id: <%s>\n", $seqpairid;
					printf "Query: %s\nExecute with <%s>, <%s>, <%s>, <%s> and <%s>?",
					$query_update_seqpair, 
					$taxid,
					$ogsid,
					$seqid      ? $seqid      : '',
					$otherseqid ? $otherseqid : '',
					$seqpairid
					;
					quit_on_q();
				}
				$sth_update_seqpair->execute($taxid, $ogsid, $seqid, $otherseqid, $seqpairid);
				if ($sth_update_seqpair->rows() == 0) {
					warn "Warning: No sequence pair updated for '%s' and IDs %d (%s) %d (%s) %d (%s), %d (seqpair)\n",
						$hdr,
						$seqid,
						$seqcol,
						$otherseqid,
						$otherseqcol,
						$seqpairid;
				}
				else {
					print "\nUpdated sequence pair $seqpairid for $hdr \n" if $verbose;
				}

			}

		}
		elsif ($sth_insert_seqpair->rows() == 1) {
			print "\nInserted new sequence pair for $hdr \n" if $verbose;
		}
		else {
			fail_and_exit('Something went wrong here...');
		}



		# skip the rest as this is old code
		next;

	
		if ($debug) {
			print $sth_ins->{Statement};
			printf "Execute this with <%s>, <%s>, <%s> and <%s>? ", $taxid, $hdr, $taxid, $hdr;
			quit_on_q();
		}
		$sth_ins->execute($taxid, $hdr, $taxid, $hdr);
		# no rows affected, nothing has been inserted
		if ($sth_ins->rows() == 0) {
			if ($debug) {
				print "no rows affected, sequence pair already exists. attempting update...\n";
				print $sth_sel->{Statement};
				printf "Execute this with <%s>, and <%s>? ", $hdr, $taxid;
				quit_on_q();
			}
			# determine the nt and aa sequence ids
			$sth_sel->execute($hdr, $taxid);
			my $ids = $sth_sel->fetchall_arrayref();
			if (scalar @$ids > 1) { fail_and_exit("Found more than one record with ID '$hdr'! Database corrupted?") }
			elsif (scalar @$ids == 0) { fail_and_exit("Could not find amino acid or nucleotide sequence with ID '$hdr'! Make sure the IDs correspond.") }
			if ($debug) {
				print "got these ids: \n";
				printf "<%s> ", defined $_ ? $_ : 'NULL' foreach (@{$$ids[0]});
				print "\n";
			}
			# get the sequence pair id
			if ($debug) {
				print $sth_selp->{Statement};
				printf "Execute this with <%s> and <%s>? ", @{$$ids[0]};
				quit_on_q();
			}
			$sth_selp->execute(@{$$ids[0]});
			my $seqpairid = $dbh->selectcol_arrayref($sth_selp);
			if ($debug) {
				print "got sequence pair ID $$seqpairid[0]\n";
				print $sth_upd->{Statement};
				printf "Execute this with <%s>, <%s>, <%s>, <%s> and <%s>? ", $taxid, $ogsid, $$ids[0][0], $$ids[0][1], $$seqpairid[0];
				quit_on_q();
			}
			$sth_upd->execute($taxid, $ogsid, $$ids[0][0], $$ids[0][1], $$seqpairid[0]);
			if ($sth_upd->rows() == 0) {
				warn "Warning: Sequence already present in database: '$hdr'\n";
				warn sprintf "These IDs messed up: %d (%s) %d (%s) %d (%s)\n",
					$$ids[0][0],
					$seqcol,
					$$ids[0][1],
					$otherseqcol,
					$$seqpairid[0],
					'sequence pair',
				;
			}
		}
		
	}
	return 1;
}

sub quit_on_q {
	my $response = <STDIN>;
	chomp $response;
	if ($response =~ /^q$/) { exit }
}


sub get_sequence_count_for_taxon {
	my $taxid = shift;
	my $q = "SELECT COUNT(*) FROM $db_table_aaseqs WHERE $db_table_aaseqs.$db_col_taxid = '$taxid'";
	my $r = db_get($q);
	return $$r[0][0];
}

sub get_list_of_taxa {
	my $q = "SELECT name, $db_col_longname FROM $db_table_taxa WHERE $db_col_core = '1'";
	my $r = db_get($q);
	return $r;
}

sub delete_sequences_with_headers {
	my $headers = shift;
	my $count = 0;
	my $q = "DELETE FROM $db_table_aaseqs WHERE header IN (SELECT HEADER from $db_table_aaseqs WHERE header = ? LIMIT 1)";
	my $dbh = get_dbh();
	my $sth = $dbh->prepare($q);
	while (shift @$headers) {
		print $sth->{Statement}, "\n" if $debug > 1;
		$sth->execute($_);
		$count += $sth->rows();
	}
	$dbh->disconnect();
	return $count;
}

sub delete_taxon {
	my $ti = shift;
	my $q = "DELETE FROM $db_table_taxa WHERE id = ?";
	db_do($q) or return 0;
	return 1;
}

sub clear_db {
	my $n = 0;
	my $dbh = get_dbh()
		or print $stderr "Fatal: Could not connect to database: $DBI::errstr\n" and exit 1;
	foreach my $table ($db_table_ests, $db_table_hmmsearch, $db_table_blast) {
		my $query = "DROP TABLE $table";
		my $sth = $dbh->prepare($query);
		$n += $sth->execute();
	}
	# disconnect ASAP
	$dbh->disconnect();
	return $n;
}

sub db_structure_present {
	if (! -e $database) { return 0 }
	my $r = get_list_of_tables();
	return grep /$db_table_orthologs/, @$r;
}

sub get_list_of_tables {
	my $q = '.tables';
	my $cmd = "$sqlite $database '$q'";
	my $r = [ `$cmd` ];
	if ($?) { die "Fatal: Could not get list of tables\n" }
	return $r;
}

sub get_reftaxon_id {
	my $shorthand = shift;
	my $result = db_get("SELECT $db_col_id FROM $db_table_taxa WHERE $db_col_name = ?", $shorthand);
	return $$result[0][0];
}

sub get_reftaxon_name {
	my $id = shift;
	my $result = db_get("SELECT $db_table_taxa.$db_col_name FROM $db_table_taxa INNER JOIN $db_table_aaseqs ON $db_table_taxa.$db_col_id = $db_table_aaseqs.$db_col_taxid WHERE $db_table_aaseqs.$db_col_id = ?", $id);
	return $$result[0][0];
}

sub species_tables_present {
	if (!-e $database) { return 0 }
	my $r = get_list_of_tables();
	unless (grep /$db_table_blast/, @$r) { return 0 }
	else    { return 1 }
}

sub get_list_of_ogs {
	my $query = "
		SELECT
			$db_table_ogs.$db_col_id,
			$db_table_taxa.$db_col_name,
			$db_table_ogs.$db_col_version,
			$db_table_ogs.$db_col_type,
			COUNT($db_table_aaseqs.$db_col_id)
		FROM $db_table_taxa
		INNER JOIN $db_table_ogs
			ON $db_table_taxa.$db_col_id = $db_table_ogs.$db_col_taxid
		INNER JOIN $db_table_aaseqs
			ON $db_table_ogs.$db_col_id = $db_table_aaseqs.$db_col_ogsid
		GROUP BY $db_table_aaseqs.$db_col_taxid
		UNION
		SELECT
			$db_table_ogs.$db_col_id,
			$db_table_taxa.$db_col_name,
			$db_table_ogs.$db_col_version,
			$db_table_ogs.$db_col_type,
			COUNT($db_table_ntseqs.$db_col_id)
		FROM $db_table_taxa
		INNER JOIN $db_table_ogs
			ON $db_table_taxa.$db_col_id = $db_table_ogs.$db_col_taxid
		INNER JOIN $db_table_ntseqs
			ON $db_table_ogs.$db_col_id = $db_table_ntseqs.$db_col_ogsid
		GROUP BY $db_table_ntseqs.$db_col_taxid
		ORDER BY $db_table_ogs.$db_col_id
	";
	return db_get($query);
}

sub insert_new_set {
	my $name = shift;
	my $descript = shift;
	my $query = "
		INSERT OR IGNORE INTO $db_table_set_details ($db_col_name, $db_col_description)
		VALUES (?, ?)
	";
	db_do($query, $name, $descript) or return 0;
	my $id = db_get("SELECT $db_col_id FROM $db_table_set_details WHERE $db_col_name = ?", $name);
	return $$id[0][0];
}

sub load_set_into_temptable {
	my $csvfile = shift;

	my $q_drop_temp = "DROP TABLE IF EXISTS $db_table_temp";
	my $q_create_temp = "CREATE TABLE $db_table_temp (
	  '$db_col_orthoid' INT,
		'$db_col_header'  INT,
		'$db_col_ogsid'   INT,
		'$db_col_setid'   INT
	)
	";

	# load data into temptable
	# needs to work directly on the species database since apparently
	# sqlite cannot use .import on attached databases. if anyone has an idea 
	# or a solution on how to mitigate this, please let me know!
	my @loadqueries = (
		".mode csv",
		$q_drop_temp,
		$q_create_temp,
		".import $csvfile $db_table_temp",
		".mode list",
	);
	foreach (@loadqueries) {
		my $cmd = qq{$sqlite -separator "," "$database" "$_"};
		if ($debug) {
			print "$cmd\n";
		}
		system("$cmd") and fail_and_exit("Could not import CSV file '$csvfile' into temporary table $db_table_temp\n");
	}
	return 1;
}

sub delete_set {
	my $setid = shift;
	db_do("DELETE FROM $db_table_orthologs WHERE $db_col_setid = ?", $setid) or croak;
	db_do("DELETE FROM $db_table_blastdbs WHERE $db_col_setid = ?", $setid) or croak;
	return db_do("DELETE FROM $db_table_set_details WHERE $db_col_id = ?", $setid);
}

sub insert_blastdb {
	my $setid = shift;
	return db_do("INSERT INTO $db_table_blastdbs ($db_col_setid, $db_col_rebuild) VALUES (?, ?)", $setid, 1);
}

sub set_blastdb_to_rebuild {
	my $setid = shift;
	my $flag  = shift;
	return db_do("UPDATE $db_table_blastdbs SET $db_col_rebuild = ? WHERE $db_col_setid = ?", $flag, $setid);
}

sub blastdb_needs_rebuilding {
	my $setid = shift;
	return check("SELECT * from $db_table_blastdbs where setid = ? and rebuild = 1", $setid);
}

1;
