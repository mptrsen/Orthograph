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
use lib $FindBin::RealBin;             # $RealBin is the directory of the original script
use Orthograph::Config;                # configuration parser getconfig()
use Data::Dumper;
use DBI;            # database interface
use DBD::mysql;     # MySQL database driver

my $config = $Orthograph::Config::config;  # copy config

# MySQL settings
my $db_dbname               = $config->{'mysql-database'};
my $db_dbpwd                = $config->{'mysql-password'};
my $db_dbserver             = $config->{'mysql-server'};
my $db_dbuser               = $config->{'mysql-username'};
my $db_timeout              = $config->{'mysql-timeout'};
my $sleep_for               = 10;

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
my $db_col_core             = 'core';
my $db_col_date             = 'date';
my $db_col_description      = 'description';
my $db_col_digest           = 'digest';
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
my $db_col_orthoid          = 'ortholog_gene_id';
my $db_col_query            = 'query';
my $db_col_setid            = 'setid';
my $db_col_sequence         = 'sequence';
my $db_col_seqpair          = 'sequence_pair';
my $db_col_start            = 'start';
my $db_col_target           = 'target';
my $db_col_taxid            = 'taxid';
my $db_col_type             = 'type';
my $db_col_version          = 'version';
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
my $stderr = *STDERR;
my $stdout = *STDOUT;
#}}}

# Check whether all information was provided in the configuration
defined $db_dbname   or fail_and_exit('MySQL database name not specified');
defined $db_dbuser   or fail_and_exit('MySQL database username not specified');
defined $db_dbpwd    or fail_and_exit('MySQL database password not specified');
defined $db_dbserver or fail_and_exit('MySQL database server not specified');

# report that this module is loaded
print "Using MySQL database '$db_dbname' on $db_dbserver\n" unless $quiet;

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

	until ($dbh = DBI->connect("DBI:mysql:$db_dbname:$db_dbserver;db_local_infile=1", $db_dbuser, $db_dbpwd)) {
		if ($slept >= $db_timeout) { 
			carp "Warning: Connection retry timeout exceeded\n" and return undef;
		}
		carp "Warning: Connection failed, retrying in $sleep_for seconds\n";
		sleep $sleep_for;
		$slept += $sleep_for;
	}

	if ($dbh) { return $dbh }
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
			`$db_col_id`           INT UNSIGNED NOT NULL AUTO_INCREMENT PRIMARY KEY,
			`$db_col_setid`        INT UNSIGNED DEFAULT NULL, UNIQUE(setid))",

		# table: ogs
		'ogs' => "CREATE TABLE `$t->{'ogs'}` (
			`$db_col_id`           INT UNSIGNED NOT NULL AUTO_INCREMENT PRIMARY KEY,
			`$db_col_type`         INT(1),
			`$db_col_taxid`        INT UNSIGNED NOT NULL,
			`$db_col_version`      VARCHAR(255),
			UNIQUE(taxid, version))",

		# table: ortholog_set
		'ortholog_set' => "CREATE TABLE `$t->{'orthologs'}` (
			`$db_col_id`               INT UNSIGNED NOT NULL AUTO_INCREMENT PRIMARY KEY,
			`$db_col_setid`            INT UNSIGNED NOT NULL,
			`$db_col_orthoid`          VARCHAR(10)  NOT NULL,
			`$db_col_seqpair`          INT UNSIGNED NOT NULL,
			UNIQUE INDEX (setid, ortholog_gene_id, sequence_pair))",

		# table: sequence_pairs
		'sequence_pairs' => "CREATE TABLE `$t->{'seqpairs'}` (
			`$db_col_id`           BIGINT UNSIGNED NOT NULL AUTO_INCREMENT PRIMARY KEY,
			`$db_col_taxid`        INT    UNSIGNED,
			`$db_col_ogsid`       INT    UNSIGNED,
			`$db_col_aaseq`       INT    UNSIGNED, UNIQUE(aa_seq),
			`$db_col_ntseq`       INT    UNSIGNED, UNIQUE(nt_seq), 
			`$db_col_date`         INT    UNSIGNED)",

		# table: sequences_aa
		'aa_sequences' => "CREATE TABLE `$t->{'aaseqs'}` (
			`$db_col_id`           BIGINT UNSIGNED NOT NULL AUTO_INCREMENT PRIMARY KEY,
			`$db_col_taxid`        INT             NOT NULL, INDEX(taxid),
			`$db_col_ogsid`        INT             NOT NULL,
			`$db_col_header`       VARCHAR(4096),            INDEX(header(24)), UNIQUE(header(24)),
			`$db_col_sequence`     MEDIUMBLOB,
			`$db_col_date`         INT UNSIGNED)",

		# table: sequences_nt
		'nt_sequences' => "CREATE TABLE `$t->{'ntseqs'}` (
			`$db_col_id`           BIGINT UNSIGNED NOT NULL AUTO_INCREMENT PRIMARY KEY,
			`$db_col_taxid`        INT             NOT NULL, INDEX(taxid),
			`$db_col_ogsid`        INT             NOT NULL,
			`$db_col_header`       VARCHAR(4096),            INDEX(header(24)), UNIQUE(header(24)),
			`$db_col_sequence`     MEDIUMBLOB,
			`$db_col_date`         INT UNSIGNED)",

		# table: set_details
		'set_details' => "CREATE TABLE `$t->{'set_details'}` (
			`$db_col_id`           INT UNSIGNED NOT NULL AUTO_INCREMENT PRIMARY KEY,
			`$db_col_name`         VARCHAR(255), UNIQUE(name),
			`$db_col_description`  BLOB)",

		# table: taxa
		'taxa' => "CREATE TABLE `$t->{'taxa'}` (
			`$db_col_id`           INT UNSIGNED NOT NULL AUTO_INCREMENT PRIMARY KEY,
			`$db_col_name`         VARCHAR(255),  UNIQUE(name),
			`$db_col_longname`     VARCHAR(255), 
			`$db_col_core`         TINYINT UNSIGNED NOT NULL)",
		
		# table: seqtypes
		'seqtypes' => "CREATE TABLE `$t->{'seqtypes'}` (
			`$db_col_id`           INT UNSIGNED NOT NULL AUTO_INCREMENT PRIMARY KEY,
			`$db_col_type`         CHAR(3),     UNIQUE(type))",
	);#}}}

	# to start off with nt and aa sequence types
	my $insert_seqtypes = "INSERT IGNORE INTO $t->{'seqtypes'} (type) VALUES ('nt'),('aa')";

	my $dbh = get_dbh();
	foreach (values %create_table) {
		print $_, ";\n" if $debug;
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
			`header`   VARCHAR(512), INDEX(header(24)),
			`sequence` MEDIUMBLOB,
			`description` VARCHAR(255))";
	$dbh->do("DROP TABLE IF EXISTS $temptable") or die "Fatal: Could not DROP TABLE $temptable\n";
	$dbh->do($create_temp_table_query) or die "Fatal: Could not CREATE TABLE $temptable\n";
}

sub load_csv_into_temptable {
	my $csvfile   = shift @_;
	my $temptable = shift @_;
	my $loadquery = "LOAD DATA LOCAL INFILE '$csvfile' 
		INTO TABLE $temptable FIELDS TERMINATED BY ',' LINES TERMINATED BY '\n'
		(
			$db_col_taxid,
			$db_col_type,
			$db_col_ogsid,
			$db_col_header,
			$db_col_sequence
		)";
	my $dbh = get_dbh();
	$dbh->do($loadquery) or die "Fatal: Could not LOAD DATA into temporary table $temptable\n";
	$dbh->disconnect;
}

sub fill_tables_from_temp_table {
	my $t = shift @_;
	my $temptable = shift @_;
	my @queries = (
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
		"INSERT IGNORE INTO $t->{'aaseqs'} (taxid, header, sequence, date) 
			SELECT $t->{'taxa'}.id, $temptable.header, $temptable.sequence, UNIX_TIMESTAMP()
			FROM $temptable
				LEFT JOIN $t->{'taxa'} 
			ON $temptable.name  = $t->{'taxa'}.name",
		# delete everything where header or sequence is NULL or empty
		"DELETE FROM $t->{'aaseqs'}
			WHERE $t->{'aaseqs'}.header IS NULL
			OR $t->{'aaseqs'}.sequence IS NULL
			OR $t->{'aaseqs'}.header = ''
			OR $t->{'aaseqs'}.sequence = ''",
		# sequence pairs (pep-nuc)
		"INSERT IGNORE INTO $t->{'seqpairs'} (taxid, ogs_id, aa_seq, nt_seq, date)
			SELECT $t->{'taxa'}.id, $t->{'ogs'}.id, $t->{'aaseqs'}.id, $t->{'ntseqs'}.id, UNIX_TIMESTAMP()
			FROM $t->{'taxa'}
			INNER JOIN $t->{'aaseqs'}
				ON $t->{'aaseqs'}.taxid = $t->{'taxa'}.id
			LEFT JOIN $t->{'ogs'}
				ON $t->{'taxa'}.id = $t->{'ogs'}.taxid
			LEFT JOIN $t->{'ntseqs'}
				ON $t->{'aaseqs'}.header = $t->{'ntseqs'}.header",
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
	my $step = 1;
	foreach (@queries) {
		print $_ . ";\n" if $debug;
		$nrows = $dbh->do($_) or die();
		print "Step $step completed\n" if $verbose;
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

=head2 get_list_of_ogs

Get list of OGS in the database

Arguments: none

Returns: array reference (list of OGS)

=cut

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


sub get_list_of_taxa {
	my $q = "SELECT `name`, `longname` FROM $db_table_taxa WHERE `core` = '1'";
	my $r = db_get($q);
	return $r;
}

sub get_sequence_count_for_taxon {
	my $taxid = shift;
	my $q = "SELECT COUNT(*) FROM $db_table_aaseqs WHERE $db_table_aaseqs.$db_col_taxid = ?";
	my $r = db_get($q, $taxid);
	return $$r[0][0];
}


=head2 get_ortholog_groups_for_set($setid)

Returns a hashref of hashrefs to create an ortholog set from. Each key in the hashref (the ortholog group ID) is a hashref of sequence_ID => sequence.

=cut

sub get_ortholog_groups_for_set {
	my $setid = shift @_ or croak "Usage: Wrapper::Mysql::get_ortholog_groups_for_set(SETID)";
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
	my $query_create_ests = "CREATE TABLE $db_table_ests ( 
		`$db_col_id`        BIGINT UNSIGNED     NOT NULL AUTO_INCREMENT,
		`$db_col_digest`    CHAR(32)            NOT NULL,           
		`$db_col_taxid`     INT UNSIGNED        NOT NULL,       
		`$db_col_type`      TINYINT(1) UNSIGNED NOT NULL,
		`$db_col_date`      INT UNSIGNED,
		`$db_col_header`    VARCHAR(4096)       NOT NULL,       
		`$db_col_sequence`  MEDIUMBLOB DEFAULT NULL,
		PRIMARY KEY (`$db_col_id`),
		INDEX (`$db_col_digest`(4)),
		INDEX (`$db_col_taxid`),
		INDEX (`$db_col_header`(10))
	) ENGINE=MYISAM";

	my $query_create_hmmsearch = "CREATE TABLE $db_table_hmmsearch (
		`$db_col_id`         BIGINT UNSIGNED NOT NULL AUTO_INCREMENT,
		`$db_col_taxid`      INT UNSIGNED NOT NULL,       
		`$db_col_query`      VARCHAR(255) NOT NULL,       
		`$db_col_target`     CHAR(32)     NOT NULL,       
		`$db_col_score`      DOUBLE       NOT NULL,
		`$db_col_evalue`     CHAR(8)      NOT NULL,
		`$db_col_log_evalue` DOUBLE       NOT NULL DEFAULT '-999',
		`$db_col_env_start`  INT UNSIGNED NOT NULL,
		`$db_col_env_end`    INT UNSIGNED NOT NULL,
		`$db_col_ali_start`  INT UNSIGNED NOT NULL,
		`$db_col_ali_end`    INT UNSIGNED NOT NULL,
		`$db_col_hmm_start`  INT UNSIGNED NOT NULL,
		`$db_col_hmm_end`    INT UNSIGNED NOT NULL,
		PRIMARY KEY (`$db_col_id`),
		INDEX (`$db_col_taxid`),
		INDEX (`$db_col_query`),
		INDEX (`$db_col_target`(4)),
		INDEX (`$db_col_log_evalue`),
		INDEX (`$db_col_score`)
	)";

	my $query_create_blast = "CREATE TABLE $db_table_blast (
		`$db_col_id`            BIGINT UNSIGNED NOT NULL AUTO_INCREMENT, 
		`$db_col_taxid`         INT UNSIGNED NOT NULL,       
		`$db_col_query`         CHAR(32)     NOT NULL,       
		`$db_col_target`        INT UNSIGNED NOT NULL,       
		`$db_col_score`         DOUBLE       NOT NULL,
		`$db_col_evalue`        CHAR(8)      NOT NULL,
		`$db_col_log_evalue`    DOUBLE       NOT NULL DEFAULT '-999',
		`$db_col_start`         INT UNSIGNED NOT NULL,
		`$db_col_end`           INT UNSIGNED NOT NULL,
		`$db_col_hmmsearch_id`  INT UNSIGNED NOT NULL,
		PRIMARY KEY (`$db_col_id`),
		INDEX (`$db_col_taxid`),
		INDEX (`$db_col_query`(4)),
		INDEX (`$db_col_target`),
		INDEX (`$db_col_log_evalue`),
		INDEX (`$db_col_hmmsearch_id`)
	)";

	# open connection
	my $dbh = get_dbh()
		or croak "Fatal: Could not connect to database: $DBI::errstr\n" and exit 1;

	# drop all tables
	foreach ($db_table_ests, $db_table_hmmsearch, $db_table_blast) {
		my $query_drop = "DROP TABLE IF EXISTS $_";
		print "$query_drop;\n" if $debug;
		my $sql = $dbh->prepare($query_drop);
		$sql->execute()
		  or croak "Fatal: Could not execute SQL query: $DBI::errstr\n" and exit(1);
	}

	# create all tables
	foreach my $query ($query_create_ests, $query_create_hmmsearch, $query_create_blast) {
		print "$query;\n" if $debug;
		my $sql = $dbh->prepare($query);
		$sql->execute()
		  or croak "Fatal: Could not execute SQL query: $DBI::errstr\n" and exit(1);
	}

	# disconnect
	$dbh->disconnect();
}


sub get_transcripts_for_species {
	my $specid = shift or croak "Usage: Wrapper::Mysql::get_transcripts(SPECIESID, TYPE)";
	my $type = shift or croak "Usage: Wrapper::Mysql::get_transcripts(SPECIESID, TYPE)";
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
	my ($hmmquery, $taxid) = @_ or croak "Usage: Wrapper::Mysql::get_hmmresults(HMMQUERY)";
	# disable query cache for this one
	my $query_get_sequences = "SELECT SQL_NO_CACHE $db_table_ests.$db_col_digest,
		  $db_table_ests.$db_col_sequence,
		  $db_table_hmmsearch.$db_col_env_start,
		  $db_table_hmmsearch.$db_col_env_end,
			$db_table_hmmsearch.$db_col_id
		FROM $db_table_ests 
		INNER JOIN $db_table_hmmsearch
		ON $db_table_hmmsearch.$db_col_target = $db_table_ests.$db_col_digest
		WHERE $db_table_hmmsearch.$db_col_query = ?
		AND $db_table_hmmsearch.$db_col_taxid = ?";

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
	my $data = &db_get($query) or croak();
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
	my $data = &db_get($query);
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
	print "getting taxid for $species_name\n";
	unless ($species_name) { croak("Usage: get_taxid_for_species(SPECIESNAME)") }
	# TODO rewrite this part using parametrized queries to protect from SQL injections?
	my $query = "SELECT id FROM $db_table_taxa WHERE core = 0 AND longname = '$species_name'";
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

=head2 insert_taxon_into_database(TAXON_NAME)

Inserts a core taxon into the database.

Returns the newly generated taxon ID.

=cut

sub insert_taxon_into_database {
	my $name = shift;
	my $core = shift;
	db_do("INSERT IGNORE INTO $db_table_taxa ($db_col_name, $db_col_core) VALUES (?, ?)", $name, $core) or croak;
	my $res = db_get("SELECT $db_col_id FROM $db_table_taxa WHERE $db_col_name = ? AND $db_col_core = ?", $name, $core);
	return $res->[0]->[0];
}

sub insert_ogs_info_into_database {
	my $type       = shift;
	my $taxid      = shift;
	my $ogsversion = shift;
	db_do("INSERT IGNORE INTO $db_table_ogs ($db_col_type, $db_col_taxid, $db_col_version) VALUES ($type, $taxid, $ogsversion)") or croak;
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
		INSERT IGNORE INTO $seqtable ($db_col_taxid, $db_col_ogsid, $db_col_header, $db_col_sequence, $db_col_date)
		SELECT $db_table_taxa.$db_col_id, $db_table_temp.$db_col_ogsid, $db_table_temp.$db_col_header, $db_table_temp.$db_col_sequence, CURRENT_TIMESTAMP
		FROM $db_table_temp 
		LEFT JOIN $db_table_taxa 
			ON $db_table_temp.$db_col_taxid = $db_table_taxa.$db_col_id
	";

	my $query_insert_seqpairs = "
		INSERT IGNORE INTO $db_table_seqpairs ($db_col_taxid, $db_col_ogsid, $seqcol, $otherseqcol, $db_col_date)
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


=head2 insert_taxon_into_table(TAXON_NAME)

Inserts a (non-core) taxon into the database. The taxon shorthand will be NULL
and the 'core' switch will be 0.

Returns the newly generated taxon ID.

=cut

sub insert_taxon_into_table {
	my $species_name = shift(@_);
	unless ($species_name) { croak("Usage: Wrapper::Mysql::insert_taxon_into_table(SPECIESNAME)") }
	if (my $taxid = &get_taxid_for_species($species_name)) { return $taxid }
	my $query = "INSERT IGNORE INTO $db_table_taxa (longname, core) VALUES (?, ?)";
	my $dbh = get_dbh()
		or return undef;
	my $sth = $dbh->prepare($query);
	$sth = execute($sth, $db_timeout, $species_name, 0);
	$dbh->disconnect();

	$g_species_id = &get_taxid_for_species($species_name) or croak;
	return $g_species_id;
}

sub create_log_evalues_view {
	unless (scalar @_ == 1) { croak 'Usage: Wrapper::Mysql::create_log_evalues_view($species_id)' }
	my $taxid = shift;
	my $query_create_log_evalues = "CREATE OR REPLACE VIEW $db_table_log_evalues AS
	  SELECT $db_table_hmmsearch.$db_col_log_evalue AS $db_col_log_evalue,
	    COUNT($db_table_hmmsearch.$db_col_log_evalue) AS `count`
	  FROM $db_table_hmmsearch
	  WHERE $db_table_hmmsearch.$db_col_taxid = ?
	  GROUP BY $db_table_hmmsearch.$db_col_log_evalue
	  ORDER BY $db_table_hmmsearch.$db_col_log_evalue";
	my $dbh = get_dbh()
		or return undef;
	my $sth = $dbh->prepare($query_create_log_evalues);
	$sth = execute($sth, $db_timeout, $taxid);
	$dbh->disconnect();
	return 1;
}
	
sub create_scores_view {
	unless (scalar @_ == 1) { croak 'Usage: Wrapper::Mysql::create_scores_view($species_id)' }
	my $taxid = shift;
	my $query_create_scores_view = "CREATE OR REPLACE VIEW $db_table_scores AS
	  SELECT $db_table_hmmsearch.$db_col_score AS $db_col_score,
	    COUNT($db_table_hmmsearch.$db_col_score) AS `count`
	  FROM $db_table_hmmsearch
	  WHERE $db_table_hmmsearch.$db_col_taxid = ?
	  GROUP BY $db_table_hmmsearch.$db_col_score
	  ORDER BY $db_table_hmmsearch.$db_col_score DESC";
	my $dbh = get_dbh()
		or return undef;
	my $sth = $dbh->prepare($query_create_scores_view);
	$sth = execute($sth, $db_timeout, $taxid);
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
			AND $db_table_species_info.$db_col_id      IS NOT NULL
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
			'env_start'    => $$line[6],
			'env_end'      => $$line[7],
			'ali_start'    => $$line[8],
			'ali_end'      => $$line[9],
			'hmm_start'    => $$line[10],
			'hmm_end'      => $$line[11],
			'header'       => $$line[12],
		});
	}
	$sth->finish();
	$dbh->disconnect();
	return $result;
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
	my $result = &db_get($query, $digest);
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
	my $specid = shift @_;
	my $real_table_ests      = $db_table_ests      . '_' . $specid;
	my $real_table_hmmsearch = $db_table_hmmsearch . '_' . $specid;
	my $real_table_blast     = $db_table_blast     . '_' . $specid;
	$db_table_ests        = $real_table_ests;
	$db_table_hmmsearch   = $real_table_hmmsearch;
	$db_table_blast       = $real_table_blast;
	return ($real_table_ests, $real_table_hmmsearch, $real_table_blast);
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
	print $q, "\n" if $debug;
	my $d = db_get($q, $digest);
	return $d->[0]->[0];
}

sub get_number_of_orthologs_for_set {
	my $setid = shift;
	my $query = "SELECT COUNT($db_col_id) FROM $db_table_orthologs WHERE $db_col_setid = ?";
	my $res = db_get($query, $setid);
	return $res->[0]->[0];
}

sub load_ests_from_file {
	my $f = shift;
	my $list = shift;

	# load data from csv file into database
	my $query_disable_keys = "ALTER TABLE $db_table_ests DISABLE KEYS";
	my $query_enable_keys = "ALTER TABLE $db_table_ests ENABLE KEYS";
	my $query = "LOAD DATA LOCAL INFILE '$f' INTO TABLE $db_table_ests FIELDS TERMINATED BY ',' ($list)";

	# open connection and do the transaction
	my $dbh = Wrapper::Mysql::get_dbh()
		or print "Fatal: Could not connect to database: $DBI::errstr\n" and exit 1;

	# flush tables, then disable indexes before loading the data. 
	# afterwards re-index the table
	print "Disabling indices on $db_table_ests...\n" if $verbose;
	$dbh->do($query_disable_keys);
	my $sth = $dbh->prepare($query);
	print "Uploading data...\n" if $verbose;
	my $num_ests = $sth->execute() or die "Fatal: MySQL transaction failed\n";
	print "Re-indexing $db_table_ests...\n" if $verbose;
	$dbh->do($query_enable_keys) or die "Fatal: MySQL transaction failed\n";
	# disconnect ASAP and die if errors
	$dbh->disconnect;
	if (defined($DBI::errstr)) { print "$DBI::errstr\n" and exit(1) }
	return $num_ests;
}

sub insert_results_into_blast_table {
	my $hits = shift;
	my $species_id = shift;
	my $hmmsearch_id = shift;
	my $hitcount = 0;

	my $query_insert_result = "INSERT IGNORE INTO $db_table_blast (
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
			$hit->{'end'},
			$hit->{'start'},
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
	my $query_insert_result = "INSERT IGNORE INTO $db_table_hmmsearch(
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
	$dbh->do("START TRANSACTION");
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
		) or fail_and_exit('Could not push to database!');
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
	`$db_col_header`   TEXT(255)    NOT NULL,
	`$db_col_sequence` TEXT)";
	my $i = "CREATE INDEX temp_header ON $db_table_temp (header(8))";
	my $dbh = get_dbh();
	print $stdout $c, "\n" if $debug > 1;
	$dbh->do("DROP TABLE IF EXISTS $db_table_temp");
	$dbh->do($c);
	$dbh->do($i);
	$dbh->disconnect();
}

sub import_ogs_into_database {
	my ($tmpfh, $hdrs, $seqtable, $otherseqtable, $seqcol, $otherseqcol, $type, $taxon, $ogsversion) = @_;
	my @q = (
		# load data into temp table
		"LOAD DATA LOCAL INFILE '$tmpfh' 
		INTO TABLE $db_table_temp 
		FIELDS TERMINATED BY ',' 
		(taxid, header, sequence)",

		# insert data into main table. IGNORE is important to avoid duplicates without throwing errors.
		"INSERT IGNORE INTO $seqtable (taxid, header, sequence)
		SELECT $db_table_taxa.id, $db_table_temp.header, $db_table_temp.sequence 
		FROM $db_table_temp 
		LEFT JOIN $db_table_taxa 
		ON $db_table_temp.taxid = $db_table_taxa.id",

		# insert sequence pairs relationships
		"INSERT INTO $db_table_seqpairs (taxid, ogs_id, $otherseqcol, $seqcol, date) 
		SELECT $db_table_taxa.id, $db_table_ogs.id, $otherseqtable.id, $seqtable.id, UNIX_TIMESTAMP()
		FROM $db_table_taxa
		RIGHT JOIN $seqtable
		ON $seqtable.taxid = $db_table_taxa.id
		LEFT JOIN $db_table_ogs
		ON $db_table_taxa.id = $db_table_ogs.taxid
		LEFT JOIN $otherseqtable
		ON $otherseqtable.header = $seqtable.header
		WHERE $db_table_taxa.id = '$taxon'
		ON DUPLICATE KEY UPDATE $db_table_seqpairs.$seqcol = $seqtable.id,
		$db_table_seqpairs.$otherseqcol = $otherseqtable.id",

		# update OGS table
		"INSERT IGNORE INTO $db_table_ogs (`type`, `taxid`, `version`) VALUES ('$type', '$taxon', '$ogsversion')"
	);

	my $dbh = get_dbh();
	my $ret = undef;
	foreach my $query (@q) {
		if ($debug) {
			print $query, "\n";
		}
		$ret = $dbh->do($query);
		unless (defined $ret) {
			return 0;
		}
	}
	return 1;
}

sub delete_sequences_with_headers {
	my $headers = shift;
	my $count = 0;
	my $q = "DELETE FROM $db_table_aaseqs WHERE header IN (SELECT HEADER from $db_table_aaseqs WHERE header = ? LIMIT 1)";
	my $dbh = get_dbh();
	my $sth = $dbh->prepare($q);
	while (shift @$headers) {
		$sth->execute($q, $_);
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
	my $q = "SELECT * FROM $db_table_seqtypes";
	my $r = db_get($q);
	if ($r) { return 1 }
	else    { return 0 }
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

sub delete_set {
	my $setid = shift;
	return db_do("DELETE FROM $db_table_set_details WHERE $db_col_id = ?", $setid);
}

1;
