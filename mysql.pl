#!/usr/bin/perl
use strict;
use warnings;
use DBI;
use DBD::mysql;
use Data::Dumper;

my $count    = 0;
my $data     = { };
my $db       = 'orthograph';
my $dbuser   = 'mpetersen';
my $dbpwd    = 'mpetersen';
my $dbserver = 'localhost';
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
	o_orthologs.ortholog_gene_id AS orthoid,
	o_hmmsearch.target           AS transcript,
	o_hmmsearch.evalue           AS hmm_eval,
	o_blast.target               AS blast_target,
	o_blast.evalue               AS blast_eval,
	o_taxa.name                  AS reftax
FROM o_aaseqs
INNER JOIN o_taxa
	ON o_aaseqs.taxid            = o_taxa.id
INNER JOIN o_blast
	ON o_aaseqs.id               = o_blast.target
INNER JOIN o_hmmsearch
	ON o_blast.query             = o_hmmsearch.target
INNER JOIN o_ests
	ON o_hmmsearch.target        = o_ests.digest
INNER JOIN o_orthologs
	ON o_hmmsearch.query         = o_orthologs.ortholog_gene_id
INNER JOIN o_sequence_pairs
	ON o_orthologs.sequence_pair = o_sequence_pairs.id
INNER JOIN o_set_details
	ON o_set_details.name        = 'notmany'
WHERE o_ests.spec              = 'Mengenilla'
ORDER BY o_hmmsearch.evalue
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

