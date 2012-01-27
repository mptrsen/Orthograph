#!/usr/bin/perl
# Documentation before the code#{{{
=head1 NAME 

preparedb - prepare a MySQL database for running Forage

=head1 SYNOPSIS

preparedb FASTAFILE

=head1 DESCRIPTION

preparedb prepares a db (duh) for running Forage. Basically, it assumes a
certain database structure, reads a fasta file and loads its content into the
database. Due to its usage of the MySQL LOAD INFILE command, it is pretty fast.

=cut#}}}

use strict;
use warnings;
use Getopt::Long;
use File::Temp;     # temporary files
use Seqload::Fasta qw(fasta2csv); # object-oriented access to fasta files, fasta2csv converter
use DBI;
use DBD::mysql;
use Seqload::Fasta;


my $dbname   = 'forage';
my $dbserver = '127.0.0.1';
my $table    = 'core_orthologs';
my $user     = undef;
my $pwd      = undef;

=head1 OPTIONS#{{{

=head2 -D database

Database name to use. Defaults to 'forage'.

=head2 -h host

Database server. Defaults to '127.0.0.1', which is equivalent to 'localhost'.

=head2 -t table

Table name to use. Defaults to 'core_orthologs'.

=head2 -u username

Username for the connection. No default; this must be set.

=head2 -p password

Password for the connection. No default; this must be set.

=cut#}}}

GetOptions(
	'D=s' => \$dbname,
	'h=s' => \$dbserver,
	't=s' => \$table,
	'u=s' => \$user,
	'p=s' => \$pwd,
	);

unless (defined $user and defined $pwd) {
	die "Fatal: You must specify username (-u) and password (-p).\n";
}

my $infile = shift @ARGV or die "Fatal: Argument INFILE missing\n";

my $tmpfile = File::Temp->new();

my $faobj = Seqload::Fasta->open($infile);
while (my ($hdr, $seq) = $faobj->next_seq()) {
	my @line = split(/\|/, $hdr);
	printf $tmpfile "%s,%s,%s,%s,%s\n", 
		$line[1], 
		$line[0], 
		$line[1] . '_prot',
		$line[2],
		$seq;
}
$faobj->close();

my $loadquery = "LOAD DATA LOCAL INFILE '$tmpfile' IGNORE INTO TABLE $table
	FIELDS TERMINATED BY ',' (
		taxon,
		hmm,
		blastdb,
		hdr,
		seq)";
my $dbh = DBI->connect("DBI:mysql:$dbname:$dbserver", $user, $pwd);
my $sql = $dbh->prepare($loadquery);
$sql->execute();
$dbh->disconnect;

exit; 
