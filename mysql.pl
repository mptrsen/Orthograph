#!/usr/bin/perl
use strict;
use warnings;
use DBI;
use DBD::mysql;
use Data::Dumper;

my $count    = 0;
my $data     = { };
my $db       = 'orthograph';
my $dbpwd    = 'malty';
my $dbserver = 'localhost';
my $dbuser   = 'malty';
my $fuzzy    = 0;
my $out      = { };
my $sort_by  = 'blasteval';
my $strict   = 1;
my @reftaxa  = qw(AAEGY CFLOR ISCAP);

#--------------------------------------------------
# this query selects all orthology candidates for which
# the reciprocal search also hit something good
#-------------------------------------------------- 
my $query = "
SELECT DISTINCT
	orthograph_orthologs.ortholog_gene_id AS orthoid,
	orthograph_hmmsearch.target           AS transcript,
	orthograph_hmmsearch.evalue           AS hmm_eval,
	orthograph_blast.target               AS blast_target,
	orthograph_blast.evalue               AS blast_eval,
	orthograph_taxa.name                  AS reftax
FROM orthograph_aaseqs
INNER JOIN orthograph_taxa
	ON orthograph_aaseqs.taxid            = orthograph_taxa.id
INNER JOIN orthograph_blast
	ON orthograph_aaseqs.id               = orthograph_blast.target
INNER JOIN orthograph_hmmsearch
	ON orthograph_blast.query             = orthograph_hmmsearch.target
INNER JOIN orthograph_ests
	ON orthograph_hmmsearch.target        = orthograph_ests.digest
INNER JOIN orthograph_orthologs
	ON orthograph_hmmsearch.query         = orthograph_orthologs.ortholog_gene_id
INNER JOIN orthograph_sequence_pairs
	ON orthograph_orthologs.sequence_pair = orthograph_sequence_pairs.id
INNER JOIN orthograph_set_details
	ON orthograph_set_details.name        = 'alltaxa'
WHERE orthograph_ests.spec              = 'Mengenilla'
ORDER BY orthograph_hmmsearch.evalue
";

# get data, store in hash->array->hash
my $dbh = DBI->connect("dbi:mysql:$db:$dbserver", $dbuser, $dbpwd);
my $sql = $dbh->prepare($query);
$sql->execute();
while (my $line = $sql->fetchrow_arrayref()) {
	push @{$$data{$$line[0]}}, {
		'digest'      => $$line[1],
		'hmmeval'     => $$line[2],
		'blasttarget' => $$line[3],
		'blasteval'   => $$line[4],
		'reftaxon'    => $$line[5]
	};
	$count++;
}
$dbh->disconnect;

# sort by blast evalue or hmmsearch evalue or the number of your mom's chest hairs
foreach my $eog (keys %$data) {
	@{$$data{$eog}} = sort {
		$$a{$sort_by} <=> $$b{$sort_by};
	} @{$$data{$eog}};
}

# output each eog with the correct transcript
# each transcript shall be assigned to only one ortholog group
foreach my $eog (keys %$data) {
	# for all hits, see if they are in the reftaxa list
	REFTAXON:
	for my $i (0..$#{$$data{$eog}}) {
		# is this reftaxon in our list?
		if ( grep /$$data{$eog}[$i]{'reftaxon'}/ , @reftaxa ) {
			# ok it's there
			# take the best (hmm|blast) evalue
			# using a hash makes sure each ortholog group gets only one transcript
			if ( defined $$out{$$data{$eog}[$i]{'digest'}} ) {
				if ( $$data{$eog}[$i]{$sort_by} < $$out{$$data{$eog}[$i]{'digest'}}{$sort_by} ) {
					$$out{$$data{$eog}[$i]{'digest'}} = {
						'eog'       => $eog,
						'hmmeval'   => $$data{$eog}[$i]{'hmmeval'},
						'blasteval' => $$data{$eog}[$i]{'blasteval'},
						'reftaxon'  => $$data{$eog}[$i]{'reftaxon'},
					};
				}
				# is this a strict search, i.e., do we need to hit ALL reftaxa?
				unless ($strict) { last; }
				else { next REFTAXON }
			}
			else {
				$$out{$$data{$eog}[$i]{'digest'}} = {
					'eog'       => $eog,
					'hmmeval'   => $$data{$eog}[$i]{'hmmeval'},
					'blasteval' => $$data{$eog}[$i]{'blasteval'},
					'reftaxon'  => $$data{$eog}[$i]{'reftaxon'},
				};
			}
		}
	}
}

foreach my $digest (keys %$out) {
	printf "%s => %s (eval %1.1e)\n", $$out{$digest}{'eog'}, $digest, $$out{$digest}{'blasteval'};
}

#print Dumper($out);
printf "%d hits %d 1:1 hits %d total hits\n", scalar keys %$data, scalar keys %$out, $count;

