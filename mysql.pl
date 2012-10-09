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
my $transcripts = { };
my $table       = [ ];
my $hashtable   = { };

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
	AND o_aaseqs.id              = o_sequence_pairs.aa_seq
INNER JOIN o_set_details
	ON o_set_details.name        = 'notmany'
WHERE o_ests.spec              = 'Mengenilla'
AND o_hmmsearch.evalue         < 1e-05
AND o_blast.evalue             < 1e-05
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
	$$hashtable{$$line[1]}{$$line[0]} = $$line[2];
	$count++;
}
$dbh->disconnect;

#--------------------------------------------------
# my $i=0;
# foreach (keys %$hashtable) { print ++$i, "\t" }
# print "\n";
# foreach my $k (keys %$hashtable) {
# 	printf "%s\t", substr $k, 0, 4;
# 	foreach my $k2 (keys %{$$hashtable{$k}}) {
# 		print $$hashtable{$k}{$k2}, "\t";
# 	}
# 	print "\n";
# }
# 
# exit;
# 
#-------------------------------------------------- 
# sort by blast evalue or hmmsearch evalue or the number of your mom's chest hairs
foreach my $eog (keys %$data) {
	@{$$data{$eog}} = sort {
		$$a{$sort_by} <=> $$b{$sort_by};
	} @{$$data{$eog}};
}

# reverse the hash so that we get a transcript->orthoid assignment
foreach my $eog (keys %$data) {
	for my $i (0..$#{$$data{$eog}}) {
		push @{$$transcripts{$$data{$eog}[$i]{'digest'}}}, {
			'orthoid'     => $eog,
			'hmmeval'     => $$data{$eog}[$i]{'hmmeval'},
			'blasttarget' => $$data{$eog}[$i]{'blasttarget'},
			'blasteval'   => $$data{$eog}[$i]{'blasteval'},
			'reftaxon'    => $$data{$eog}[$i]{'reftaxon'},
		};
	}
}

my @keys_data = keys %$data;
my @keys_transcripts = keys %$transcripts;

# make a HTML table!
print "<html><head><title>Table</title></head><body>\n";
print "<table>\n";
print "<tr><th>&nbsp;</th>\n";	# first cell is empty 
printf "<th align='left'>%s</th>\n", $_ foreach @keys_data;	# header
print "</tr>\n";
for (my $x = 0; $x < scalar @keys_transcripts; ++$x) {
	print "<tr>\n";
	print "<td><b>", $keys_transcripts[$x], "</b></td>\n";
	for (my $y = 0; $y < scalar @keys_data; ++$y) {

		# get the matching transcript from the list
		foreach my $hit (@{$$data{$keys_data[$y]}}) {
			if ($$hit{'digest'} eq $keys_transcripts[$x]) {
				print "<td>", $$hit{'hmmeval'}, "</td>\n";
				last;
			}
			else {
				print "<td>NULL</td>\n";
			}
			#printf "%s\n", defined $$_{'digest'} ? $$_{'digest'} : 'NULL';
			#print grep $$_{'digest'} =~ /$keys_transcripts[$x]/, @{$$data{$keys_data[$y]}};
		}

		#--------------------------------------------------
		# if (defined $$transcripts{$keys_transcripts[$x]}) {
		# 	$$table[$x][$y] = $$data{$keys_data[$y]}[$x]{'hmmeval'};
		# 	print "<td>", defined $$transcripts{$keys_transcripts[$x]}[$y]{'hmmeval'} ? $$transcripts{$keys_transcripts[$x]}[$y]{'hmmeval'} : '<b>NULL</b>', "</td>\n";
		# }
		# else {
		# 	$$table[$x][$y] = '<b>NULL</b>';
		# 	print "<td>NULL</td>\n";
		# }
		#-------------------------------------------------- 
	}
	print "</tr>\n";
}
print "</table>\n";
print "</body></html>\n";
warn Dumper($data);
exit;

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

#--------------------------------------------------
# foreach my $digest (keys %$out) {
# 	printf "%s => %s (eval %1.1e)\n", $$out{$digest}{'eog'}, $digest, $$out{$digest}{'blasteval'};
# }
# 
# #print Dumper($out);
# printf "%d hits %d 1:1 hits %d total records\n", scalar keys %$data, scalar keys %$out, $count;
#-------------------------------------------------- 

printf "%s => %s\n", $_, $$data{$_}[0]{'digest'} foreach keys %$data;
