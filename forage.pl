#!/usr/bin/perl
#--------------------------------------------------
# FORAGE: Find Orthologs using Reciprocity Among Genes and ESTs
# Copyright 2011 Malte Petersen <mptrsen@uni-bonn.de>
# 
# Forage is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
# 
# Forage is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License along with
# Forage. If not, see http://www.gnu.org/licenses/.
#-------------------------------------------------- 

use strict;   # make me write good code
use warnings; # cry if something seems odd
use Getopt::Long;
use Bio::Perl;
use Bio::SeqIO;
use Bio::PrimarySeq;

# Forage: Find Orthologs using Reciprocity Among Genes and ESTs
my $version = 0.00001;
my $estfile = '';
my @seqobjs;

GetOptions(	'estfile=s' => \$estfile,
						'E=s'				=> \$estfile,
);

my $header = <<EOF;
Forage: Find Orthologs using Reciprocity Among Genes and ESTs
Copyright 2011 Malte Petersen <mptrsen\@uni-bonn.de>
Version $version
Using EST file $estfile
EOF
print $header;
&translate($estfile);

# Sub: translate
# Translate a sequence into all six reading frames, print out everything
# Expects: scalar string filename
# Returns: -
sub translate {
	my ($infile) = @_;
	my $estfileobj = Bio::SeqIO->new( '-file' => $infile, '-format' => 'fasta');
	while (my $seqobj = $estfileobj->next_seq) {
		my @sixframesobj = Bio::SeqUtils->translate_6frames($seqobj);
		print '>', $seqobj->display_id, "\n";
		print $seqobj->seq, "\n";
		for(my $i=0; $i < @sixframesobj; ++$i) {
			my $protseq = $sixframesobj[$i];
			print '>', $protseq->display_id, "|frame_", $i+1, "\n";
			print $protseq->seq, "\n";
		}
	}
}
