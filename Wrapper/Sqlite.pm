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
use lib $FindBin::Bin;                 # $Bin is the directory of the original script
use Orthograph::Config;                # configuration parser getconfig()
use Data::Dumper;
use DBI;
use DBD::SQLite;

my $config = $Orthograph::Config::config;  # copy config

# db settings
my $database                = $config->{'sqlite-database'};
my $db_timeout              = 600;
my $sqlite                  = $config->{'sqlite-program'};
my $sleep_for               = 1;
my $db_dbuser               = $config->{'username'} || $ENV{"LOGNAME"} || $ENV{"USER"} || getpwuid $<;

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
my $db_table_taxa           = $config->{'db_table_taxa'};
my $db_table_temp           = $config->{'db_table_temp'};
my $db_col_aaseq            = 'aa_seq';
my $db_col_date             = 'date';
my $db_col_digest           = 'digest';
my $db_col_end              = 'end';
my $db_col_env_end          = 'env_end';
my $db_col_env_start        = 'env_start';
my $db_col_evalue           = 'evalue';
my $db_col_hmm_end          = 'hmm_end';
my $db_col_hmm_start        = 'hmm_start';
my $db_col_header           = 'header';
my $db_col_id               = 'id';
my $db_col_log_evalue       = 'log_evalue';
my $db_col_score            = 'score';
my $db_col_name             = 'name';
my $db_col_ntseq            = 'nt_seq';
my $db_col_ogsid            = 'ogs_id';
my $db_col_ogsversion       = 'ogs_version';
my $db_col_orthoid          = 'ortholog_gene_id';
my $db_col_query            = 'query';
my $db_col_setid            = 'setid';
my $db_col_sequence         = 'sequence';
my $db_col_seqpair          = 'sequence_pair';
my $db_col_start            = 'start';
my $db_col_target           = 'target';
my $db_col_taxid            = 'taxid';
my $db_col_type             = 'type';
my $outdir                  = $config->{'output-directory'};
my $orthoset                = $config->{'ortholog-set'};
my $quiet                   = $config->{'quiet'};
my $reftaxa                 = $config->{'reference-taxa'};
# substitution character for selenocysteine, which normally leads to blast freaking out
my $u_subst                 = $config->{'substitute-u-with'};
my $sets_dir                = $config->{'sets-dir'};
my $species_name            = $config->{'species-name'};
my $g_species_id            = undef;	# global variable
my $verbose                 = $config->{'verbose'};
my $debug                   = $config->{'debug'};
my $stdout = *STDOUT;
my $stderr = *STDERR;
#}}}



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
		return $dbh;
	}
	return undef;
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
  # connect and fetch stuff
	my $dbh = get_dbh()
		or return undef;
	my $sth = $dbh->prepare($query);
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
	unless ($query) { croak "Usage: db_do(QUERY)\n" }
	my @fields = @_;
	my $dbh = get_dbh()
		or return undef;
	my $sth = $dbh->prepare($query);
	$sth = execute($sth, $db_timeout, @fields);
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
		'blastdbs' => "CREATE TABLE `$t->{'blastdbs'}` (
			`id`           INTEGER PRIMARY KEY,
			`setid`        INTEGER UNSIGNED DEFAULT NULL UNIQUE,
			`blastdb_path` TEXT(255) DEFAULT NULL)",
		
		# table: ogs
		'ogs' => "CREATE TABLE `$t->{'ogs'}` (
			`id`           INTEGER PRIMARY KEY,
			`type`         INT(1),
			`taxid`        INTEGER UNSIGNED NOT NULL UNIQUE,
			`version`      TEXT(255))",
		
		# table: ortholog_set
		'ortholog_set' => "CREATE TABLE `$t->{'orthologs'}` (
			`id`               INTEGER PRIMARY KEY,
			`setid`            INTEGER UNSIGNED NOT NULL,
			`ortholog_gene_id` TEXT(10)  NOT NULL,
			`sequence_pair`    INTEGER UNSIGNED NOT NULL,
			UNIQUE (setid, ortholog_gene_id, sequence_pair))",

		# table: sequence_pairs
		'sequence_pairs' => "CREATE TABLE `$t->{'seqpairs'}` (
			`id`           INTEGER PRIMARY KEY,
			`taxid`        INTEGER    UNSIGNED,
			`ogs_id`       INTEGER    UNSIGNED,
			`aa_seq`       INTEGER    UNSIGNED UNIQUE,
			`nt_seq`       INTEGER    UNSIGNED UNIQUE, 
			`date`         INTEGER    UNSIGNED DEFAULT CURRENT_TIMESTAMP,
			`user`         INTEGER    UNSIGNED)",

		# table: sequences_aa
		'aa_sequences' => "CREATE TABLE `$t->{'aaseqs'}` (
			`id`           INTEGER PRIMARY KEY,
			`taxid`        INTEGER             NOT NULL, 
			`header`       TEXT(512)    UNIQUE,
			`sequence`     MEDIUMBLOB,
			`user`         INTEGER UNSIGNED,
			`date`         INTEGER UNSIGNED DEFAULT CURRENT_TIMESTAMP)",

		# table: sequences_nt
		'nt_sequences' => "CREATE TABLE `$t->{'ntseqs'}` (
			`id`           INTEGER PRIMARY KEY,
			`taxid`        INTEGER             NOT NULL, 
			`header`       TEXT(512)    UNIQUE,
			`sequence`     MEDIUMBLOB,
			`user`         INTEGER UNSIGNED,
			`date`         INTEGER UNSIGNED DEFAULT CURRENT_TIMESTAMP)",

		# table: set_details
		'set_details' => "CREATE TABLE `$t->{'set_details'}` (
			`id`           INTEGER PRIMARY KEY,
			`name`         TEXT(255) UNIQUE,
			`description`  BLOB)",

		# table: taxa
		'taxa' => "CREATE TABLE `$t->{'taxa'}` (
			`id`           INTEGER PRIMARY KEY,
			`name`         TEXT(20)  UNIQUE,
			`longname`     TEXT(255), 
			`core`         TINYINTEGER UNSIGNED NOT NULL)",
		
		# table: users
		'users' => "CREATE TABLE `$t->{'users'}` (
			`id`           INTEGER PRIMARY KEY,
			`name`         TEXT(255) UNIQUE)",
		# table: seqtypes
		'seqtypes' => "CREATE TABLE `$t->{'seqtypes'}` (
			`id`           INTEGER PRIMARY KEY,
			`type`         TEXT(3)     UNIQUE)",
	);#}}}

	my @indices = (
		# indices for sequences_aa
"CREATE INDEX IF NOT EXISTS $t->{'aaseqs'}_taxid  ON $t->{'aaseqs'} (taxid)",
"CREATE INDEX IF NOT EXISTS $t->{'ntseqs'}_taxid  ON $t->{'ntseqs'} (taxid)",
"CREATE INDEX IF NOT EXISTS $t->{'aaseqs'}_header  ON $t->{'aaseqs'} (header)",
"CREATE INDEX IF NOT EXISTS $t->{'ntseqs'}_header  ON $t->{'ntseqs'} (header)",
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
	foreach (values %create_table) {
		print $_, ";\n" if $verbose;
		$dbh->do($_) or die "Could not exec query: $!\n";
	}
	foreach (@indices) {
		print $_, ";\n" if $verbose;
		$dbh->do($_) or die "Could not exec query: $!\n";
	}
	$dbh->do($insert_seqtypes);	# start off with 'nt' and 'aa' seqtypes
	$dbh->disconnect;
}


sub create_temp_table {
	my $temptable = shift @_;
	my $create_temp_table_query = "CREATE TABLE $temptable (
			`name`        TEXT(255),
			`longname`    TEXT(255),
			`orthoset`    TEXT(255),
			`orthoid`     TEXT(255),
			`blastdb`     TEXT(255),
			`header`      TEXT(512),
			`sequence`    MEDIUMBLOB,
			`description` TEXT(255))";
	my $create_temp_indices_query = "BEGIN;
CREATE INDEX IF NOT EXISTS ${temptable}_name ON $temptable (name);
CREATE INDEX IF NOT EXISTS ${temptable}_orthoset ON $temptable (orthoset);
CREATE INDEX IF NOT EXISTS ${temptable}_orthoid ON $temptable (orthoid);
CREATE INDEX IF NOT EXISTS ${temptable}_header ON $temptable (header);
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
		print $_, "\n" if $debug;
		system qq{$sqlite -separator "," $database "$_"} and die "Fatal: Could not import CSV file into temporary table $temptable\n";
	}
}

sub fill_tables_from_temp_table {
	my $t = shift @_;
	my $temptable = shift @_;
	my @queries = (
		# user name
		"INSERT OR IGNORE INTO $t->{'users'} (name) VALUES ('$db_dbuser')",
		# taxa (name, longname)
		"INSERT OR IGNORE INTO $t->{'taxa'} (name, longname, core) 
			SELECT DISTINCT $temptable.name, $temptable.longname, 1 
			FROM $temptable",
		# set name + description
		"INSERT OR IGNORE INTO $t->{'set_details'} (name, description)
			SELECT DISTINCT $temptable.orthoset, $temptable.description 
			FROM $temptable LIMIT 1",
		# blast databases
		"INSERT OR IGNORE INTO $t->{'blastdbs'} (setid, blastdb_path) 
			SELECT DISTINCT $t->{'set_details'}.id, $temptable.blastdb 
			FROM $temptable
			LEFT JOIN $t->{'set_details'} 
				ON $t->{'set_details'}.name = $temptable.orthoset",
		# pep sequences
		"INSERT OR IGNORE INTO $t->{'aaseqs'} (taxid, header, sequence, user, date) 
			SELECT $t->{'taxa'}.id, $temptable.header, $temptable.sequence, $t->{'users'}.id, CURRENT_TIMESTAMP
			FROM $temptable
				LEFT JOIN $t->{'taxa'} 
			ON $temptable.name  = $t->{'taxa'}.name
				INNER JOIN $t->{'users'}
			ON $t->{'users'}.name = '$db_dbuser'",
		# delete everything where header or sequence is NULL or empty
		"DELETE FROM $t->{'aaseqs'}
			WHERE $t->{'aaseqs'}.header IS NULL
			OR $t->{'aaseqs'}.sequence IS NULL
			OR $t->{'aaseqs'}.header = ''
			OR $t->{'aaseqs'}.sequence = ''",
		# sequence pairs (pep-nuc)
		"INSERT OR IGNORE INTO $t->{'seqpairs'} (taxid, ogs_id, aa_seq, nt_seq, date, user)
			SELECT $t->{'taxa'}.id, $t->{'ogs'}.id, $t->{'aaseqs'}.id, $t->{'ntseqs'}.id, CURRENT_TIMESTAMP, $t->{'users'}.id
			FROM $t->{'taxa'}
			INNER JOIN $t->{'aaseqs'}
				ON $t->{'aaseqs'}.taxid = $t->{'taxa'}.id
			LEFT JOIN $t->{'ogs'}
				ON $t->{'taxa'}.id = $t->{'ogs'}.taxid
			LEFT JOIN $t->{'ntseqs'}
				ON $t->{'aaseqs'}.header = $t->{'ntseqs'}.header
			INNER JOIN $t->{'users'}
				ON $t->{'users'}.name = '$db_dbuser'",
		# orthologous groups
		"INSERT OR IGNORE INTO $t->{'orthologs'} (setid, ortholog_gene_id, sequence_pair) 
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
		print $_ . ";\n" if $verbose;
		$nrows = $dbh->do($_) or fail_and_exit("Query failed: $_");
		if ($verbose) {
			($nrows > 0) ? printf("Query OK, %d rows affected\n", $nrows) : print "Query OK\n";
		}
	}
	$dbh->disconnect;
	return $nrows;
}

sub get_number_of_cogs_for_set {
	my $setn = shift @_;
	my $q = "SELECT COUNT(DISTINCT $db_table_orthologs.ortholog_gene_id) FROM $db_table_orthologs INNER JOIN $db_table_set_details ON $db_table_orthologs.setid = $db_table_set_details.id WHERE $db_table_set_details.name = ?";
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
	my $query = "SELECT * FROM $db_table_set_details";
	my $data = db_get($query);
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
	my $query = "SELECT DISTINCT $db_table_taxa.name , $db_table_ogs.version
		FROM $db_table_aaseqs
		INNER JOIN $db_table_seqpairs
			ON $db_table_aaseqs.id  = $db_table_seqpairs.aa_seq
		INNER JOIN $db_table_taxa
			ON $db_table_seqpairs.taxid = $db_table_taxa.id
		INNER JOIN $db_table_ogs
			ON $db_table_taxa.id = $db_table_ogs.taxid"
	;
	my $data = db_get($query);
	foreach my $item (@$data) {
		$ogslist{$$item[0]} = $$item[1];
	}
	return(\%ogslist);
}#}}}


=head2 get_ortholog_groups_for_set($setid)

Returns a hashref of hashrefs to create an ortholog set from. Each key in the hashref (the ortholog group ID) is a hashref of sequence_ID => sequence.

=cut

sub get_ortholog_groups_for_set {
	my $setid = shift @_ or croak "Usage: Wrapper::Sqlite::get_ortholog_groups_for_set(SETID)";
	my $data = {};
	my $query = "SELECT o.ortholog_gene_id, a.id, a.sequence
		FROM $db_table_orthologs         AS o
    INNER JOIN $db_table_seqpairs    AS p
    ON o.sequence_pair = p.id
    INNER JOIN $db_table_aaseqs      AS a
    ON a.id = p.aa_seq
    INNER JOIN $db_table_set_details AS d
    ON d.id = o.setid
    WHERE d.id = ?";

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
	my $query_create_ests = "CREATE TABLE $db_table_ests ( 
		`$db_col_id`        INTEGER NOT NULL PRIMARY KEY,
		`$db_col_digest`    TEXT(32)     NOT NULL,           
		`$db_col_taxid`     UNSIGNED INTEGER NOT NULL,       
		`$db_col_type`      UNSIGNED TINYINT(4) NOT NULL,
		`$db_col_date`      UNSIGNED INT,
		`$db_col_header`    TEXT      NOT NULL,       
		`$db_col_sequence`  MEDIUMBLOB DEFAULT NULL
		)";

	my $query_create_hmmsearch = "CREATE TABLE $db_table_hmmsearch (
		`$db_col_id`         INTEGER NOT NULL PRIMARY KEY,
		`$db_col_taxid`      UNSIGNED INTEGER NOT NULL,       
		`$db_col_query`      TEXT(255) NOT NULL,       
		`$db_col_target`     TEXT(32)     NOT NULL,       
		`$db_col_score`      DOUBLE       NOT NULL,
		`$db_col_evalue`     TEXT(8)      NOT NULL,
		`$db_col_log_evalue` DOUBLE       NOT NULL DEFAULT '-999',
		`$db_col_env_start`  UNSIGNED INTEGER NOT NULL,
		`$db_col_env_end`    UNSIGNED INTEGER NOT NULL,
		`$db_col_hmm_start`  UNSIGNED INTEGER NOT NULL,
		`$db_col_hmm_end`    UNSIGNED INTEGER NOT NULL
		)";

	my $query_create_blast = "CREATE TABLE $db_table_blast (
		`$db_col_id`            INTEGER NOT NULL PRIMARY KEY,
		`$db_col_taxid`         UNSIGNED INTEGER NOT NULL,       
		`$db_col_query`         TEXT(32)     NOT NULL,       
		`$db_col_target`        UNSIGNED INTEGER NOT NULL,       
		`$db_col_score`         DOUBLE       NOT NULL,
		`$db_col_evalue`        TEXT(8)      NOT NULL,
		`$db_col_log_evalue`    DOUBLE       NOT NULL DEFAULT '-999',
		`$db_col_start`         UNSIGNED INTEGER NOT NULL,
		`$db_col_end`           UNSIGNED INTEGER NOT NULL
		)";

	my @query_create_indices = (
"CREATE INDEX IF NOT EXISTS ${db_table_ests}_header ON $db_table_ests (header)",
"CREATE INDEX IF NOT EXISTS ${db_table_ests}_digest ON $db_table_ests ($db_col_digest)",
"CREATE INDEX IF NOT EXISTS ${db_table_ests}_taxid ON $db_table_ests ($db_col_taxid)",
"CREATE INDEX IF NOT EXISTS ${db_table_ests}_header ON $db_table_ests ($db_col_header)",
"CREATE INDEX IF NOT EXISTS ${db_table_hmmsearch}_taxid ON $db_table_hmmsearch ($db_col_taxid)",
"CREATE INDEX IF NOT EXISTS ${db_table_hmmsearch}_query ON $db_table_hmmsearch ($db_col_query)",
"CREATE INDEX IF NOT EXISTS ${db_table_hmmsearch}_target ON $db_table_hmmsearch ($db_col_target)",
"CREATE INDEX IF NOT EXISTS ${db_table_hmmsearch}_evalue ON $db_table_hmmsearch ($db_col_log_evalue)",
"CREATE INDEX IF NOT EXISTS ${db_table_hmmsearch}_score ON $db_table_hmmsearch ($db_col_score)",
"CREATE INDEX IF NOT EXISTS ${db_table_blast}_taxid ON $db_table_blast ($db_col_taxid)",
"CREATE INDEX IF NOT EXISTS ${db_table_blast}_query ON $db_table_blast ($db_col_query)",
"CREATE INDEX IF NOT EXISTS ${db_table_blast}_target ON $db_table_blast ($db_col_target)",
"CREATE INDEX IF NOT EXISTS ${db_table_blast}_evalue ON $db_table_blast ($db_col_log_evalue)"
	);

	# open connection
	my $dbh = get_dbh()
		or croak "Fatal: Could not connect to database: $DBI::errstr\n" and exit 1;

	# drop all tables
	foreach ($db_table_ests, $db_table_hmmsearch, $db_table_blast) {
		my $query_drop = "DROP TABLE IF EXISTS $_";
		print "$query_drop;\n" if $verbose;
		my $sql = $dbh->prepare($query_drop);
		$sql->execute()
		  or croak "Fatal: Could not execute SQL query: $DBI::errstr\n" and exit(1);
	}

	# create all tables
	foreach my $query ($query_create_ests, $query_create_hmmsearch, $query_create_blast, @query_create_indices) {
		print "$query;\n" if $verbose;
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
		  $db_table_hmmsearch.env_end
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
	my %setlist = ();
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
	my $data = &db_get($query) or croak();
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
	my $result = &db_get("SELECT COUNT(*) FROM $db_table_ests");

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
	my $taxids = &db_get("SELECT id FROM $db_table_taxa WHERE name IN ($taxa_string)");

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
	my $aaseqs = &db_get("SELECT COUNT(*) FROM  $db_table_aaseqs WHERE $db_table_aaseqs.taxid IN ($taxids_string)");

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
	my $aaseqs = &db_get("SELECT $db_table_aaseqs.id, $db_table_aaseqs.sequence FROM  $db_table_aaseqs WHERE $db_table_aaseqs.taxid IN ($taxids_string)");

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
	my $query = "SELECT $db_table_taxa.id FROM $db_table_taxa WHERE core = 0 AND longname = '$species_name'";
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
	my $query = "SELECT id FROM $db_table_set_details WHERE name = '$setname'";
	my $result = &db_get($query);
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

=head2 insert_taxon_into_table(TAXON_NAME)

Inserts a (non-core) taxon into the database. The taxon shorthand will be NULL
and the 'core' switch will be 0.

Returns the newly generated taxon ID.

=cut

sub insert_taxon_into_table {
	my $species_name = shift(@_);
	unless ($species_name) { croak("Usage: Wrapper::Sqlite::insert_taxon_into_table(SPECIESNAME)") }
	if (my $taxid = get_taxid_for_species($species_name)) { return $taxid }
	my $query = "INSERT OR IGNORE INTO $db_table_taxa (longname, core) VALUES (?, ?)";
	my $dbh = get_dbh()
		or return undef;
	my $sth = $dbh->prepare($query);
	$sth = execute($sth, $db_timeout, $species_name, 0);
	$dbh->disconnect();

	$g_species_id = &get_taxid_for_species($species_name) or croak;
	return $g_species_id;
}


sub load_ests_from_file {
	my $f = shift;
	my $list = shift;

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
	my $dbh = get_dbh();
	foreach ($q_drop_temp, $q_create_temp) {
		print $_, "\n" if $debug;
		$dbh->do($_);
	}
	$dbh->disconnect;
	load_csv_into_temptable($f, $db_table_temp);

	# transfer data from temptable into main table
	my $q_transfer = "INSERT INTO $db_table_ests (
	  '$db_col_digest',
		'$db_col_taxid',
		'$db_col_type',
		'$db_col_date',
		'$db_col_header',
		'$db_col_sequence')
		SELECT 
	  $db_table_temp.'$db_col_digest',
		$db_table_temp.'$db_col_taxid',
		$db_table_temp.'$db_col_type',
		$db_table_temp.'$db_col_date',
		$db_table_temp.'$db_col_header',
		$db_table_temp.'$db_col_sequence'
		FROM $db_table_temp
	";
	$dbh = get_dbh();
	print $q_transfer, "\n" if $debug;
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
		$db_table_aaseqs.$db_col_header, $db_table_aaseqs.$db_col_sequence
		FROM $db_table_aaseqs
		INNER JOIN $db_table_seqpairs
			ON $db_table_aaseqs.$db_col_id = $db_table_seqpairs.$db_col_aaseq
		INNER JOIN $db_table_orthologs
			ON $db_table_seqpairs.$db_col_id = $db_table_orthologs.$db_col_seqpair
		AND   $db_table_orthologs.$db_col_setid = ?
		AND   $db_table_orthologs.$db_col_orthoid = ?";
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
		$db_table_ntseqs.$db_col_header, $db_table_ntseqs.$db_col_sequence
		FROM $db_table_ntseqs
		INNER JOIN $db_table_seqpairs
			ON $db_table_ntseqs.$db_col_id = $db_table_seqpairs.$db_col_ntseq
		INNER JOIN $db_table_orthologs
			ON $db_table_seqpairs.$db_col_id = $db_table_orthologs.$db_col_seqpair
		AND   $db_table_orthologs.$db_col_setid = ?
		AND   $db_table_orthologs.$db_col_orthoid = ?";
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
	print $query . "\n" if $debug;

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
	print $query . "\n" if $debug;

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
	print $query . "\n" if $debug;

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
	my $data = &db_get($query);
	my @result;
	push(@result, ${shift(@$data)}[0]) while @$data;
	return @result;
}
	
=head2 get_reference_sequence(scalar int ID)

Fetches the amino acid sequence for ID from the database. Returns a string.

=cut

sub get_reference_sequence {
	my $id = shift @_ or croak "Usage: get_reference_sequence(ID)\n";
	my $query = "SELECT $db_col_sequence 
		FROM $db_table_aaseqs
		WHERE $db_col_id = '$id'";
	my $result = &db_get($query);
	return $result->[0]->[0];
}

sub get_transcript_for {
	my $digest = shift @_ or croak "Usage: get_transcript_for(ID)\n";
	my $query  = "SELECT $db_col_sequence
		FROM $db_table_ests
		WHERE $db_col_digest = ?";
	my $result = &db_get($query, $digest);
	return $result->[0]->[0];
}

sub get_nucleotide_transcript_for {
	my $digest = shift @_ or croak "Usage: get_transcript_for(ID)\n";
	my $query  = "SELECT $db_col_header
		FROM $db_table_ests
		WHERE $db_col_digest = ?";
	my $result = &db_get($query, $digest);
	# remove the revcomp/translate portion
	print "translated header: <$result->[0]->[0]>\n" if $debug;
	(my $original_header = $result->[0]->[0]) =~ s/ ?(\[revcomp]:)?\[translate\(\d\)\]$//;
	print "original header: <$original_header>\n" if $debug;
	$query = "SELECT $db_col_sequence
		FROM $db_table_ests
		WHERE $db_col_header = ?";
	$result = &db_get($query, $original_header);
	return $result->[0]->[0];
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
	print Dumper($data); exit;
}

=head2 get_real_table_names(int ID, string EST_TABLE, string HMMSEARCH_TABLE, string BLAST_TABLE)

Renames the table names according to ID. Returns a list of the three table names.

=cut

sub get_real_table_names {
	my $specid = shift @_;
	my $real_table_ests      = $db_table_ests      . '_' . $specid;
	my $real_table_hmmsearch = $db_table_hmmsearch . '_' . $specid;
	my $real_table_blast     = $db_table_blast     . '_' . $specid;
	$db_table_ests        = $real_table_ests;
	$db_table_hmmsearch   = $real_table_hmmsearch;
	$db_table_blast       = $real_table_blast;
	return ($db_table_ests, $db_table_hmmsearch, $db_table_blast);
}


=head2 get_scores_list

Returns list of scores as present in the scores view

=cut

sub get_scores_list {
	my $specid = shift;
	my $setid = shift;
	my $r = db_get("SELECT score FROM $db_table_hmmsearch
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
	print $q, "\n" if $debug;
	my $d = db_get($q, $digest);
	return $d->[0]->[0];
}

sub insert_results_into_blast_table {
	my $hits = shift;
	my $species_id = shift;
	my $hitcount = 0;

	my $query_insert_result = "INSERT OR IGNORE INTO $db_table_blast (
		`$db_col_taxid`,
		`$db_col_query`,
		`$db_col_target`,
		`$db_col_score`,
		`$db_col_evalue`,
		`$db_col_log_evalue`,
		`$db_col_start`,
		`$db_col_end`
		) VALUES (
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
			$hit->{'end'},
			$hit->{'start'},
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
			$hit->{'env_start'},  # start of hit domain on the target seq
			$hit->{'env_end'},    # end of hit domain on the target seq
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
	my $q = "CREATE TABLE $db_table_temp (
	`taxid`    TEXT(5)      NOT NULL,
	`header`   TEXT(255) NOT NULL,
	`sequence` TEXT)";
	my $dbh = get_dbh();
	print $stdout $q, "\n" if $debug;
	$dbh->do("DROP TABLE IF EXISTS $db_table_temp");
	$dbh->do($q);
	$dbh->disconnect();
}

sub import_ogs_into_database {
	my ($f, $hdrs, $seqtable, $otherseqtable, $seqcol, $otherseqcol, $type, $taxon, $ogsversion) = @_;

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

	my $query_select_pair = "
		SELECT 
			$db_table_seqpairs.$db_col_id,
			$seqtable.$db_col_id,
			$otherseqtable.$db_col_id
		FROM $db_table_seqpairs
		LEFT JOIN $seqtable
			ON $db_table_seqpairs.$seqcol = $seqtable.$db_col_id
		LEFT JOIN $otherseqtable
			ON $seqtable.$db_col_header = $otherseqtable.$db_col_header
		LEFT JOIN $db_table_taxa
			ON $db_table_seqpairs.$db_col_taxid = $db_table_taxa.$db_col_id
		WHERE $db_table_taxa.$db_col_id = ?
		AND $seqtable.$db_col_header = ?
		UNION 
		SELECT 
			$db_table_seqpairs.$db_col_id,
			$seqtable.$db_col_id,
			$otherseqtable.$db_col_id
		FROM $db_table_seqpairs
		LEFT JOIN $seqtable
			ON $db_table_seqpairs.$seqcol = $seqtable.$db_col_id
		LEFT JOIN $otherseqtable
			ON $seqtable.$db_col_header = $otherseqtable.$db_col_header
		LEFT JOIN $db_table_taxa
			ON $db_table_seqpairs.$db_col_taxid = $db_table_taxa.$db_col_id
		WHERE $db_table_taxa.$db_col_id = ?
		AND $otherseqtable.$db_col_header = ?
	";

	my $query_update_pair = "
		UPDATE $db_table_seqpairs 
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
	my $sth_ins = $dbh->prepare($query_insert_pair);
	my $sth_sel = $dbh->prepare($query_select_pair);
	my $sth_upd = $dbh->prepare($query_update_pair);
	

	foreach my $hdr (@$hdrs) {
	# for each header
	#		get the corresponding id from the other seq table, if there is one
	#		try to insert both ids 
	#		test if a row was updated
	#		if not:
	#		update the id
	
		print $sth_ins->{Statement} if $debug;
		printf "Execute this with %s, %s, %s and %s?\n", $taxon, $hdr, $taxon, $hdr;
		<STDIN>;
		$sth_ins->execute($taxon, $hdr, $taxon, $hdr);
		# no rows affected, nothing has been inserted
		if ($sth_ins->rows() == 0) {
			# determine the seqpairs id
			print $sth_sel->{Statement} if $debug;
			$sth_sel->execute($taxon, $hdr, $taxon, $hdr);
			if ($sth_sel->rows() > 1) { croak "Fatal: SELECT statement returned more than one row!\n" }
			elsif ($sth_sel->rows() == 0) { croak "Fatal: SELECT statement returned zero rows!\n" }
			my $ids = $sth_sel->fetchall_arrayref();
			if ($debug) {
				print "got these ids: \n";
				print Dumper $ids;
				print $sth_upd->{Statement};
				printf "Execute this with %s, %s, %s, %s and %s?\n", $taxon, $ogsversion, $$ids[0][1], $$ids[0][2], $$ids[0][0];
				<STDIN>;
			}
			$sth_upd->execute($taxon, $ogsversion, $$ids[0][1], $$ids[0][2], $$ids[0][0]);
			if ($sth_upd->rows() == 0) { croak "Fatal: UPDATE didn't affect anything (no rows updated)!\n" }
		}
		
	}


		# update OGS table
		"INSERT OR IGNORE INTO $db_table_ogs (`type`, `taxid`, `version`) VALUES ('$type', '$taxon', '$ogsversion')";

}

sub get_sequence_count_for_taxon {
	my $taxon = shift;
	my $q = "SELECT COUNT(*) FROM $db_table_aaseqs WHERE $db_table_aaseqs.taxid = '$taxon'";
	my $r = db_get($q);
	return $$r[0][0];
}

sub get_list_of_taxa {
	my $q = "SELECT name, longname FROM $db_table_taxa WHERE core = '1'";
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
	if (!-e $database) { return 0 }
	my $r = get_list_of_tables();
	if ($r) { return 1 }
	else    { return 0 }
}

sub get_list_of_tables {
	my $q = '.tables';
	my $cmd = "$sqlite $database '$q'";
	my $r = [ `$cmd` ];
	if ($!) { die "Fatal: Could not get list of tables: $!\n" }
	unless (grep /$db_table_orthologs/, @$r) { return 0 }
	return 1;
}

1;
