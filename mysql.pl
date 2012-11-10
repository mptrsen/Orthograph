#!/usr/bin/env perl
use strict;
use warnings;
use DBI;
use DBD::mysql;
use Data::Dumper;
use IO::File;

my $count       = 0;
my $data_by_orthoid        = { };
my $db          = 'orthograph';
my $dbuser      = $ENV{LOGNAME};
my $dbpwd       = $ENV{LOGNAME};
my $dbserver    = 'localhost';
my $setname     = 'mosquitoes';
my $speciesname = 'Mengenilla';
my $out         = { };
my $sort_by     = 'blasteval';
my @reftaxa     = qw(CQUIN AAEGY );
my $data_by_transcript = { };
my $data_by_evalue = { };
my $table       = [ ];
my $strict      = 0;

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
	AND orthograph_aaseqs.id              = orthograph_sequence_pairs.aa_seq
INNER JOIN orthograph_set_details
	ON orthograph_set_details.name        = '$setname'
WHERE orthograph_ests.spec              = '$speciesname'
AND orthograph_hmmsearch.evalue         < 1e-05
AND orthograph_blast.evalue             < 1e-05
ORDER BY orthograph_hmmsearch.evalue
";

# get data, store in hash->array->hash
my $dbh = DBI->connect("dbi:mysql:$db:$dbserver", $dbuser, $dbpwd);
my $sql = $dbh->prepare($query);
$sql->execute();
while (my $line = $sql->fetchrow_arrayref()) {
	# the orthoid shall be the key
	push @{$$data_by_orthoid{$$line[0]}}, {	# 0: orthoid
		'digest'      => $$line[1],	# 1: SHA256 digest
		'hmmeval'     => sprintf("%e", $$line[2]),	# 2: hmmsearch evalue
		'blasttarget' => $$line[3],	# 3: blast target id	
		'blasteval'   => sprintf("%e", $$line[4]),	# 4: blast evalue
		'reftaxon'    => $$line[5]	# 5: reference taxon
	};
	# the same structure, but keyed by transcript digest
	push @{$$data_by_transcript{$$line[1]}}, {
		'orthoid'     => $$line[0],
		'hmmeval'     => sprintf("%e", $$line[2]),	# 2: hmmsearch evalue
		'blasttarget' => $$line[3],
		'blasteval'   => sprintf("%e", $$line[4]),
		'reftaxon'    => $$line[5]
	};
	# the same structure, but keyed by e-value
	push @{$$data_by_evalue{sprintf("%e", $$line[2])}}, {
		'orthoid'     => $$line[0],
		'digest'      => $$line[1],
		'blasttarget' => $$line[3],
		'blasteval'   => sprintf("%e", $$line[4]),
		'reftaxon'    => $$line[5]
	};
	$count++;
}
$dbh->disconnect;


# sort by blast evalue or hmmsearch evalue or the number of your mom's chest hairs
foreach my $eog (keys %$data_by_orthoid) {
	@{$$data_by_orthoid{$eog}} = sort {
		$$a{$sort_by} <=> $$b{$sort_by}
	} @{$$data_by_orthoid{$eog}};
}
foreach my $eval (keys %$data_by_evalue) {
	@{$$data_by_evalue{$eval}} = sort {
		$$a{$sort_by} <=> $$b{$sort_by}
	} @{$$data_by_evalue{$eval}};
}
foreach my $transcript (keys %$data_by_transcript) {
	@{$$data_by_transcript{$transcript}} = sort {
		$$a{$sort_by} <=> $$b{$sort_by}
	} @{$$data_by_transcript{$transcript}};
}

#--------------------------------------------------
# # a normal table for screen output
# # header
# printf "%-8s %-9s %-32s %-8s\n",
# 	'#hmmeval',
# 	'orthoid',
# 	'est_digest',
# 	'blasteval';
# # data
# foreach my $hmmeval (sort {$a <=> $b} keys %$data_by_evalue) {
# 	for my $i (0 .. scalar(@{$$data_by_evalue{$hmmeval}}) - 1) {
# 		printf "%-8.1e %8s %32s %-8.1e %5s\n",
# 			$hmmeval,
# 			$$data_by_evalue{$hmmeval}[$i]{'orthoid'},
# 			$$data_by_evalue{$hmmeval}[$i]{'digest'},
# 			$$data_by_evalue{$hmmeval}[$i]{'blasteval'},
# 			$$data_by_evalue{$hmmeval}[$i]{'reftaxon'},
# 		;
# 	}
# }
#-------------------------------------------------- 

# keys to the hashes
my @keys_orthoids    = keys %$data_by_orthoid    ;
my @keys_transcript = keys %$data_by_transcript;
my @keys_evalues     = keys %$data_by_evalue     ;


# make a HTML table!
my $fh = IO::File->new("$speciesname.html", 'w');
print $fh "<html><head><title>$speciesname</title></head><body>\n";
print $fh "<table>\n";
print $fh "<tr><th>&nbsp;</th>\n";	# first cell is empty 
printf $fh "<th align='left'>%s</th>\n", $_ foreach @keys_orthoids;	# header
print $fh "</tr>\n";
# columns
for (my $x = 0; $x < scalar @keys_transcript; ++$x) {
	print $fh "<tr>\n";
	print $fh "<td><b>", $keys_transcript[$x], "</b></td>";
	# rows
	for (my $y = 0; $y < scalar @keys_orthoids; ++$y) {
		# get the matching transcript from the list
		my $flag = 0;
		foreach my $hit (@{$$data_by_transcript{$keys_transcript[$x]}}) {
			if ($$hit{'orthoid'} eq $keys_orthoids[$y]) {
				print $fh "<td>", $$hit{'hmmeval'}, "</td>" if $flag == 0;
				# connect these in the result table
				push @{$$table[$x][$y]}, {
					'hmmeval'     => $$hit{'hmmeval'},
					'blasteval'   => $$hit{'blasteval'},
					'blasttarget' => $$hit{'blasttarget'},
					'reftaxon'    => $$hit{'reftaxon'},
				};
				$flag = 1;
			}
		}
		if ($flag == 0) {
			# this is a non-match
			print $fh "<td>NULL</td>";
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


# TODO walk through the e-values and examine each hit pair

foreach my $heval (sort {$a <=> $b} keys %$data_by_evalue) {
	# each (hmmsearch) e-value is associated with a list of hits
	foreach my $i (0 .. scalar( @{$$data_by_evalue{$heval}} ) - 1) {
		# only if both orthoid and transcript have not been assigned before
		if (
			grep( /$$data_by_evalue{$heval}[$i]{'digest'}/, @keys_transcript )
			and grep( /$$data_by_evalue{$heval}[$i]{'orthoid'}/, @keys_orthoids )
		)
		{
			# is the reftaxon of this hit part of the reftaxa list?
			if ( grep /$$data_by_evalue{$heval}[$i]{'reftaxon'}/, @reftaxa ) {
				printf "%.1e : %s and %s (%s: %.1e, hit #%d)\n",
					$heval,
					$$data_by_evalue{$heval}[$i]{'orthoid'},
					$$data_by_evalue{$heval}[$i]{'digest'},
					$$data_by_evalue{$heval}[$i]{'reftaxon'},
					$$data_by_evalue{$heval}[$i]{'blasteval'},
					$i,
				;
				for (0 .. scalar(@keys_orthoids) -1 ) {
					if ($keys_orthoids[$_] eq $$data_by_evalue{$heval}[$i]{'orthoid'}) {
						splice(@keys_orthoids, $_, 1) and last;
					}
				}
				for (0 .. scalar(@keys_transcript) -1 ) {
					if ($keys_transcript[$_] eq $$data_by_evalue{$heval}[$i]{'digest'}) {
						splice(@keys_transcript, $_, 1) and last;
					}
				}
				# make sure they only get assigned once
				#delete $$data_by_transcript{ $$data_by_evalue{$heval}[$i]{'digest'} };
				#delete $$data_by_orthoid{ $$data_by_evalue{$heval}[$i]{'orthoid'} };
			}
		}
	}
}

exit;

# take only the reference taxa
my $num_reftaxa = scalar @reftaxa;
TR:
for my $x (0 .. $#keys_transcript) {
	OG:
	for my $y (0 .. $#keys_orthoids) {
		# skip non-matches and already-removed entries
		next if not defined $$table[$x][$y] or $$table[$x][$y] == 0;

		# local list of reftaxa for this match
		my @this_reftaxa;
		push @this_reftaxa, $$_{'reftaxon'} foreach (@{$$table[$x][$y]});
		
		# intersection of @reftaxa and @this_reftaxa
		my %union = my %isect = ();
		foreach my $el (@reftaxa) { $union{$el} = 1 }
		foreach my $el (@this_reftaxa) { $isect{$el} = 1 if ($union{$el}) }
		

		# This orthoid-transcript combination matches all reference taxa!
		if (scalar keys %isect == $num_reftaxa) { 
			# strict match
			print "+\t$keys_orthoids[$y] and $keys_transcript[$x]: ";
			printf "%s ", $_ foreach keys %isect;
			printf "with %.2e\n", $$table[$x][$y][0]{'hmmeval'};
			# but it may still be redundant... fuck.
			# remove it from the table to hopefully eliminate redundancy
			splice @{$$table[$_]}, $y, 1 foreach (0 .. $#{$table});
			splice @$table, $x, 1;
			last;
		}
		# this one matches at least one reference taxon
		elsif (scalar keys %isect) {
			# fuzzy match
			print "=\t$keys_orthoids[$y] and $keys_transcript[$x]: ";
			printf "%s ", $_ foreach keys %isect;
			printf "with %.2e\n", $$table[$x][$y][0]{'hmmeval'};
			# but it may still be redundant... fuck.
			# remove it from the table to hopefully eliminate redundancy
			splice @{$$table[$_]}, $y, 1 foreach (0 .. $#{$table});
			splice @$table, $x, 1;
			last;
		}
		# this one matches none
		else {
			print "-\t$keys_orthoids[$y] and $keys_transcript[$x]: ";
			printf "%s ", $_ foreach @this_reftaxa;
			print "\n";
		}

	}
}
print "Reftaxa: @reftaxa\n";
exit;

# output each eog with the correct transcript
# each transcript shall be assigned to only one ortholog group
foreach my $eog (keys %$data_by_orthoid) {
	# for all hits, see if they are in the reftaxa list
	REFTAXON:
	for my $i (0..$#{$$data_by_orthoid{$eog}}) {
		# is this reftaxon in our list?
		if ( grep /$$data_by_orthoid{$eog}[$i]{'reftaxon'}/ , @reftaxa ) {
			# ok it's there
			# take the best (hmm|blast) evalue
			# using a hash makes sure each ortholog group gets only one transcript
			if ( defined $$out{$$data_by_orthoid{$eog}[$i]{'digest'}} ) {
				if ( $$data_by_orthoid{$eog}[$i]{$sort_by} < $$out{$$data_by_orthoid{$eog}[$i]{'digest'}}{$sort_by} ) {
					$$out{$$data_by_orthoid{$eog}[$i]{'digest'}} = {
						'eog'       => $eog,
						'hmmeval'   => $$data_by_orthoid{$eog}[$i]{'hmmeval'},
						'blasteval' => $$data_by_orthoid{$eog}[$i]{'blasteval'},
						'reftaxon'  => $$data_by_orthoid{$eog}[$i]{'reftaxon'},
					};
				}
				# is this a strict search, i.e., do we need to hit ALL reftaxa?
				unless ($strict) { last; }
				else { next REFTAXON }
			}
			else {
				$$out{$$data_by_orthoid{$eog}[$i]{'digest'}} = {
					'eog'       => $eog,
					'hmmeval'   => $$data_by_orthoid{$eog}[$i]{'hmmeval'},
					'blasteval' => $$data_by_orthoid{$eog}[$i]{'blasteval'},
					'reftaxon'  => $$data_by_orthoid{$eog}[$i]{'reftaxon'},
				};
			}
		}
	}
}

printf "%s => %s\n", $_, $$data_by_orthoid{$_}[0]{'digest'} foreach keys %$data_by_orthoid;
