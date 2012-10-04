#!/usr/bin/perl
use strict;
use warnings;
use DBI;
use DBD::mysql;
use Data::Dumper;

my $db       = 'orthograph';
my $dbserver = 'localhost';
my $dbuser   = 'malty';
my $dbpwd    = 'malty';
my $data     = { };
my $count    = 0;
my $sort_by  = 'blasteval';
my @reftaxa  = ('AAEGY');

#--------------------------------------------------
# get query
#-------------------------------------------------- 
my $query = "
SELECT DISTINCT
	orthograph_hmmsearch.target           AS transcript,
	orthograph_orthologs.ortholog_gene_id AS orthoid,
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
	ON orthograph_set_details.name        = 'test'
WHERE orthograph_ests.spec              = 'Mengenilla moldrzyki'
ORDER BY orthograph_hmmsearch.evalue
";

# get data, store in hash->array->hash
my $dbh = DBI->connect("dbi:mysql:$db:$dbserver", $dbuser, $dbpwd);
my $sql = $dbh->prepare($query);
$sql->execute();
while (my $line = $sql->fetchrow_arrayref()) {
	#print join("\t", @$line), "\n";
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

printf "%d hits %d total\n", scalar keys %$data, $count;

my $out = { };
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
			if ( ( defined $$out{$$data{$est}[$i]{'orthoid'}} ) and ( $$data{$est}[$i]{$sort_by} < $$out{$$data{$est}[$i]{'orthoid'}}{$sort_by} ) ) {
				$$out{$$data{$est}[$i]{'orthoid'}} = {
					'est'       => $est,
					'hmmeval'   => $$data{$est}[$i]{'hmmeval'},
					'blasteval' => $$data{$est}[$i]{'blasteval'},
					'reftaxon'  => $$data{$est}[$i]{'reftaxon'},
				};
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
