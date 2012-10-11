#!/usr/bin/perl
use strict;
use warnings;
use DBI;
use DBD::mysql;
use Data::Dumper;

my $count       = 0;
my $data        = { };
my $db          = 'orthograph';
my $dbuser      = $ENV{LOGNAME};
my $dbpwd       = $ENV{LOGNAME};
my $dbserver    = 'localhost';
my $setname     = 'notmany';
my $speciesname = 'Mengenilla';
my $out         = { };
my $sort_by     = 'blasteval';
my @reftaxa     = qw(AMELL ACEPH DPULE);
my $transcripts = { };
my $table       = [ ];
my $strict      = 1;

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
	ON o_set_details.name        = '$setname'
WHERE o_ests.spec              = '$speciesname'
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
	$count++;
}
$dbh->disconnect;

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
my $fh = IO::File->new("$speciesname.html", 'w');
print $fh "<html><head><title>$speciesname</title></head><body>\n";
print $fh "<table>\n";
print $fh "<tr><th>&nbsp;</th>\n";	# first cell is empty 
printf $fh "<th align='left'>%s</th>\n", $_ foreach @keys_data;	# header
print $fh "</tr>\n";
# columns
for (my $x = 0; $x < scalar @keys_transcripts; ++$x) {
	print $fh "<tr>\n";
	print $fh "<td><b>", $keys_transcripts[$x], "</b></td>\n";
	# rows
	for (my $y = 0; $y < scalar @keys_data; ++$y) {
		# get the matching transcript from the list
		my $flag = 0;
		foreach my $hit (@{$$transcripts{$keys_transcripts[$x]}}) {
			if ($$hit{'orthoid'} eq $keys_data[$y]) {
				print $fh "<td>", $$hit{'hmmeval'}, "</td>\n" if $flag == 0;
				$flag = 1;
				# connect these in the result table
				push @{$$table[$x][$y]}, {
					'hmmeval'     => $$hit{'hmmeval'},
					'blasteval'   => $$hit{'blasteval'},
					'blasttarget' => $$hit{'blasttarget'},
					'reftaxon'    => $$hit{'reftaxon'},
				};
			}
		}
		if ($flag == 0) {
			# this is a non-match
			print $fh "<td>NULL</td>\n";
			$$table[$x][$y] = 0;
		}
		else {
			$flag = 0;
			@{$$table[$x][$y]} = sort {$$a{'hmmeval'} <=> $$b{'hmmeval'}} @{$$table[$x][$y]};
		}
	}
	print $fh "</tr>\n";
}
print $fh "</table>\n";
print $fh "</body></html>\n";
# close file
undef $fh;
#warn Dumper($table);

# take only the reference taxa
my $num_reftaxa = scalar @reftaxa;
for my $y (0 .. $#keys_transcripts) {
	TR:
	for my $x (0 .. $#keys_data) {
		# skip non-matches and already-removed entries
		next if not defined $$table[$x][$y] or $$table[$x][$y] == 0;


		# local list of reftaxa for this match
		my @this_reftaxa;
		push @this_reftaxa, $$_{'reftaxon'} foreach (@{$$table[$x][$y]});
		
		# intersection of @reftaxa and @this_reftaxa
		my %union = my %isect = ();
		foreach my $e (@reftaxa) { $union{$e} = 1 }
		foreach my $e (@this_reftaxa) { $isect{$e} = 1 if ($union{$e}) }
		

		# This orthoid-transcript combination matches all reference taxa!
		if ($strict and scalar keys %isect == $num_reftaxa) { 
			print "strict for $keys_data[$y] and $keys_transcripts[$x]\n";
			# but it may still be redundant... fuck.
			# remove it from the table to hopefully eliminate redundancy
			splice @$table, $x, 1;
			last;
		}
		elsif (!$strict and scalar keys %isect) {
			print "fuzzy for  $keys_data[$y] and $keys_transcripts[$x]\n";
			# but it may still be redundant... fuck.
			# remove it from the table to hopefully eliminate redundancy
			splice @$table, $x, 1;
			last;
		}

	}
}
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

printf "%s => %s\n", $_, $$data{$_}[0]{'digest'} foreach keys %$data;
