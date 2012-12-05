#!/usr/bin/perl
use strict;
use warnings;
use DBI;
use DBD::mysql;
use Benchmark;
use Time::HiRes qw( time );
use Data::Dumper;
use List::Util qw( min max sum );
use autodie;

my $dbname   = 'orthograph';
my $dbserver = 'localhost';
my $username = 'mpetersen';
my $pwd      = 'mpetersen';
my $max_bin_size = 500;
my $limit    = 2000;
my $offset   = 0;
my $from     = -999;
my $plus     = 10;
my $until    = 100;
my $t0       = scalar time;
my $t1;
my $start    = scalar time;

my $query_get_logevalues = "SELECT log_evalue, count FROM o_logevalues";
my $dbh      =	DBI->connect("DBI:mysql:$dbname:$dbserver", $username, $pwd);
my $sql = $dbh->prepare($query_get_logevalues);
$sql->execute();
my $d = $sql->fetchall_arrayref();
$dbh->disconnect();
my %num_of_logevalues = ();

foreach my $row (@$d) {
	$num_of_logevalues{$$row[0]} = $$row[1];
}

my @logevalues = sort { $a <=> $b } keys %num_of_logevalues;

my @bin = ( );
my $sum = 0;
my $query = "select o_hmmsearch.evalue, o_orthologs.ortholog_gene_id, o_hmmsearch.start, o_hmmsearch.end, o_blast.target, o_blast.evalue, o_taxa.name from o_logevalues left join o_hmmsearch on o_logevalues.log_evalue = o_hmmsearch.log_evalue left join o_ests on o_hmmsearch.target = o_ests.digest left join o_orthologs on o_hmmsearch.query = o_orthologs.ortholog_gene_id left join o_blast on o_hmmsearch.target = o_blast.query left join o_aaseqs on o_blast.target = o_aaseqs.id left join o_taxa on o_aaseqs.taxid = o_taxa.id  left join o_set_details on o_orthologs.setid = o_set_details.id where o_set_details.id='1' and o_hmmsearch.taxid = '45' and o_hmmsearch.log_evalue between ? and ?";
my $query_single = "select o_hmmsearch.evalue, o_orthologs.ortholog_gene_id, o_hmmsearch.start, o_hmmsearch.end, o_blast.target, o_blast.evalue, o_taxa.name from o_logevalues left join o_hmmsearch on o_logevalues.log_evalue = o_hmmsearch.log_evalue left join o_ests on o_hmmsearch.target = o_ests.digest left join o_orthologs on o_hmmsearch.query = o_orthologs.ortholog_gene_id left join o_blast on o_hmmsearch.target = o_blast.query left join o_aaseqs on o_blast.target = o_aaseqs.id left join o_taxa on o_aaseqs.taxid = o_taxa.id  left join o_set_details on o_orthologs.setid = o_set_details.id where o_set_details.id='1' and o_hmmsearch.taxid = '45' and o_hmmsearch.log_evalue = ?";
$dbh =	DBI->connect("DBI:mysql:$dbname:$dbserver", $username, $pwd);
$sql = $dbh->prepare($query);
for my $i (0..$#logevalues) {
	# if the number of rows would not exceed the bin size
	if ($sum + $num_of_logevalues{$logevalues[$i]} < $max_bin_size or $i >= scalar(@logevalues) - 1) {
		push @bin, $logevalues[$i];
		$sum += $num_of_logevalues{$logevalues[$i]};
		printf "got %f: %d x\n", $logevalues[$i], $num_of_logevalues{$logevalues[$i]};
	}
	elsif (scalar @bin == 0) {
		# this one exceeds the bin size alone, fetch
		push @bin, $logevalues[$i];
		print "overload on single bin. captured ", sum( @num_of_logevalues{@bin} ), " so far\n";
		my $val = $logevalues[$i];
		printf "fetching single %f...\n", $logevalues[$i];
		$t0 = scalar time;
		my $nrows = &fetch_single( $logevalues[$i] );
		$t1 = scalar time;
		printf "%f: fetched %10d rows in %.2f secs\n", $logevalues[$i], $nrows, $t1-$t0;
		@bin = ();
		$sum = 0;

	}
	# otherwise, fetch and clear
	else {
		print "overload. captured ", sum( @num_of_logevalues{@bin} ), " so far\n";
		my $min = min(@bin);
		my $max = max(@bin);
		printf "fetching from %f to %f...\n", $min, $max;
		$t0 = scalar time;
		my $nrows = &fetch( $min, $max );
		$t1 = scalar time;
		printf "%f to %f: fetched %10d rows in %.2f secs\n", $min, $max, $nrows, $t1-$t0;
		@bin = ();
		$sum = 0;
		push @bin, $logevalues[$i];
		$sum += $num_of_logevalues{$logevalues[$i]};
		printf "got %f: %d x\n", $logevalues[$i], $num_of_logevalues{$logevalues[$i]};
	}
}
print "captured ", sum( @bin ), "\n";
print "did ", scalar @logevalues, "\n";

exit;

sub fetch {
	my $min = shift;
	my $max = shift;
	$sql->execute( $min, $max );
	my $rows = $sql->fetchall_arrayref();
	my $nrows = scalar @$rows;
	return $nrows;
}

sub fetch_single {
	my $logeval = shift;
	$sql = $dbh->prepare($query_single);
	$sql->execute( $logeval );
	my $rows = $sql->fetchall_arrayref();
	my $nrows = scalar @$rows;
	$sql = $dbh->prepare($query);
	return $nrows;
}

# bins of equal size, won't work
for (my $i = $from; $i < $until; $i += $plus) {

	$t0 = scalar time;
	$sql->execute( $i, $i + $plus );
	my $rows = $sql->fetchall_arrayref();
	my $nrows = scalar @$rows;
	$t1 = scalar time;
	printf "%4d to %4d: fetched %10d rows in %.2f secs\n", $i, $i + $plus, $nrows, $t1-$t0;
}

exit;

$dbh      =	DBI->connect("DBI:mysql:$dbname:$dbserver", $username, $pwd);
$sql = $dbh->prepare($query);

# single evalues, won't work
foreach my $logeval (@logevalues) {
	$sql->execute($logeval);
	my $rows = $sql->fetchall_arrayref();
	print Dumper($rows); exit;
}

# limit and offset, won't work
open my $fh, '>', 'benchmark.csv';
while (1) {
	#$t0 = Benchmark->new();
	$t0 = scalar time;
	$sql->execute( $limit, $offset ) or die();
	my $rows = $sql->fetchall_arrayref();
	$t1 = scalar time;
	#$t1 = Benchmark->new();
	#printf "block %10d took %s\n", $offset, timestr(timediff($t1, $t0));
	printf "block %10d took %5.2f secs, %5.2f total\n", $offset, $t1 - $t0, $t1 - $start;
	printf $fh "%d,%5.2f,%5.2f\n", $offset, $t1-$t0, $t1-$start;
	last unless defined $$rows[0][0];
	$offset += $limit;
	
}
$dbh->disconnect();
close $fh;
