#!/usr/bin/perl
use strict;
use warnings;
use DBI;
use DBD::mysql;
use Data::Dumper;

my $count    = 0;
my $data     = { };
my $db       = 'orthograph';
my $dbpwd    = 'mpetersen';
my $dbserver = 'localhost';
my $dbuser   = 'mpetersen';
my $out = { };
my $sort_by  = 'blasteval';
my @reftaxa  = ('AAEGY');

#--------------------------------------------------
# get query
#-------------------------------------------------- 
my $query = "
SELECT DISTINCT
	o_hmmsearch.target           AS transcript,
	o_orthologs.ortholog_gene_id AS orthoid,
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
	ON o_set_details.name        = 'mosq'
WHERE o_ests.spec              = 'Mengenilla'
ORDER BY o_hmmsearch.evalue
";

# get data, store in hash->array->hash
my $dbh = DBI->connect("dbi:mysql:$db:$dbserver", $dbuser, $dbpwd);
my $sql = $dbh->prepare($query);
$sql->execute();
while (my $line = $sql->fetchrow_arrayref()) {
	push @{$$data{$$line[0]}}, {
		'orthoid'     => $$line[1],
		'hmmeval'     => $$line[2],
		'blasttarget' => $$line[3],
		'blasteval'   => $$line[4],
		'reftaxon'    => $$line[5]
	};
	$count++;
}
$dbh->disconnect;

# sort by blast evalue or hmmsearch evalue or the number of your mom's chest hairs
foreach my $est (keys %$data) {
	@{$$data{$est}} = sort {
		$$a{$sort_by} <=> $$b{$sort_by};
	} @{$$data{$est}};
}

# output each eog with the correct transcript
# each transcript shall be assigned to only one ortholog group
foreach my $est (keys %$data) {
	# for all hits, see if they are in the reftaxa list
	for my $i (0..$#{$$data{$est}}) {
		# is this reftaxon in our list?
		if ( grep /$$data{$est}[$i]{'reftaxon'}/ , @reftaxa ) {
			# ok it's there
			# take the largest (hmm|blast) evalue
			# using a hash makes sure each ortholog group gets only one transcript
			if ( defined $$out{$$data{$est}[$i]{'orthoid'}} ) {
				if ( $$data{$est}[$i]{$sort_by} < $$out{$$data{$est}[$i]{'orthoid'}}{$sort_by} ) {
					$$out{$$data{$est}[$i]{'orthoid'}} = {
						'est'       => $est,
						'hmmeval'   => $$data{$est}[$i]{'hmmeval'},
						'blasteval' => $$data{$est}[$i]{'blasteval'},
						'reftaxon'  => $$data{$est}[$i]{'reftaxon'},
					};
				}
			}
			else {
				$$out{$$data{$est}[$i]{'orthoid'}} = {
					'est'       => $est,
					'hmmeval'   => $$data{$est}[$i]{'hmmeval'},
					'blasteval' => $$data{$est}[$i]{'blasteval'},
					'reftaxon'  => $$data{$est}[$i]{'reftaxon'},
				};
			}
		}
	}
}

print Dumper($out);
printf "%d hits %d total\n", scalar keys %$data, $count;
printf "%d eogs %d total\n", scalar keys %$out, $count;

