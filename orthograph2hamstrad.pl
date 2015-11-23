#!/usr/bin/perl
use strict;
use warnings;
use autodie;
use File::Spec::Functions;

my $usage = <<"__EOT__";

$0 - convert Orthograph results into HaMStRad format.

USAGE:

  $0 INPUTDIRECTORY OUTPUTDIRECTORY

Writes all Fasta files in INPUTDIRECTORY/aa and INPUTDIRECTORY/nt to
OUTPUTDIRECTORY/aa and OUTPUTDIRECTORY/nt, respectively. OUTPUTDIRECTORY
must exist, aa and nt subdirectory will be created.
Fasta headers will be corresponding.

__EOT__


my $indir  = shift @ARGV or print $usage and exit;
my $outdir = shift @ARGV or print $usage and exit;

# input and output directories
my $aaind = catdir($indir,  'aa');
my $ntind = catdir($indir,  'nt');
my $ntod  = catdir($outdir, 'nt');
my $aaod  = catdir($outdir, 'aa');

# see whether we have summarized data or not, change input directory accordingly
if (! -d $aaind) {
	if (-d catdir($indir, 'aa_summarized') and -d  catdir($indir, 'nt_summarized')) {
		$aaind = catdir($indir,  'aa_summarized');
		$ntind = catdir($indir,  'nt_summarized');
		$ntod  = catdir($outdir, 'nt_summarized');
		$aaod  = catdir($outdir, 'aa_summarized');
	}
	else {
		die "Fatal: could not find aa and nt subdirectories. Is this an Orthograph output directory?\n"
	}
}

# create nt output dir unless it exists
create_dir($ntod);

# create aa output dir unless it exists
create_dir($aaod);

# directory handles
my $aadh;
opendir $aadh, $aaind;

while (my $f = readdir($aadh)) {

	# skip everything but aa fasta files
	next unless $f =~ /\.fa$/;
	(my $basename = $f) =~ s/\.aa(\.summarized)?(\.mafft)?\.fa$//;

	my $summarized = defined $1 ? $1 : '';
	my $mafft      = defined $2 ? $2 : '';
	
	my $aafile = catfile($aaind, $f);
	my $ntfile = catfile($ntind, $basename . '.nt' . $summarized . $mafft . '.fa');

	# the corresponding nt file must exist
	if (!-e catfile($ntfile)) {
		warn "Warning: $ntfile does not exist, skipping...\n";
		next;
	}

	# read them both in
	my $aadata = fasta2arrayref($aafile);
	my $ntdata = fasta2arrayref($ntfile);
	
	# open output files
	open my $ntofh, '>', catfile($ntod, $basename . $summarized . $mafft . '.nt.fa');
	open my $aaofh, '>', catfile($aaod, $basename . $summarized . $mafft . '.aa.fa');

	foreach my $item (@$aadata) {

		# get corresponding nt sequence
		my $ntitem = find_sequence_for_taxon($ntdata, $item->{'tax'});

		# print aa headers and aa sequences to aa output file
		printf $aaofh ">%s\n%s\n", magic_hamstr_format($item), $item->{'seq'};

		# print aa headers, but nt sequences to nt output file
		printf $ntofh ">%s\n%s\n", magic_hamstr_format($item), $ntitem->{'seq'};

	}

	close $ntofh;
	close $aaofh;

}

exit;

# create a directory or leave it alone if it exists
sub create_dir {
	my $dir = shift;
	if (-e $dir) {
		if (! -d $dir) {
			die "Fatal: output dir $dir exists, but is not a directory!\n";
		}
	}
	else {
		mkdir $dir;
	}
	return 1;
}


# find corresponding item in a list of items that belongs to the taxon of
# another item
sub find_sequence_for_taxon {
	my $data = shift;
	my $tax = shift;
	foreach my $item (@$data) {
		if (defined $item->{'tax'} and $item->{'tax'} eq $tax) { return $item }
	}
	# not found
	warn "Warning: in $data->[0]->{'cog'}: no nucleotide sequence for $tax, leaving it empty\n";
	return { 'seq' => '' };
}

# slurp a fasta file into fancy arrayref structure
# arguments: filename
# returns: arrayref of hashrefs, each having:
#		'cog' => cog id
#		'tax' => taxon name
#		'hdr' => header string
#		'seq' => sequence string
#		'ref' => boolean, is a reftaxon or not
sub fasta2arrayref {
	my $f = shift;
	my $seqs = [ ];
	my $fh = Seqload::Fasta->open($f);
	while (my ($h, $s) = $fh->next_seq()) {
		my @sectors = split /&&/, $h;
		my @first_sector = split /\|/, $sectors[0];
		my $is_reftaxon = $first_sector[-1] eq '.' ? 1 : 0;
		push @$seqs, { 'cog' => $first_sector[0], 'tax' => $first_sector[1], 'hdr' => $h, 'seq' => $s, 'ref' => $is_reftaxon };
	}
	$fh->close();
	return $seqs;
}

# convert a string from orthograph to hamstr format
# arguments: hashref containing the following keys:
#		'hdr' => header string
#		'reftaxon' => 0 or 1
# returns: header string in hamstr format
sub magic_hamstr_format {
	my $args = shift;
	my $h = $args->{'hdr'};
	my $is_reftaxon = $args->{'ref'};
	my @sectors = split /&&/, $h;
	my @first_sector = split /\|/, $sectors[0];
	my $hamstr_hdr = '';
	if ($is_reftaxon) {
		$hamstr_hdr = join('|', $first_sector[0], $first_sector[1], $first_sector[2]);
	}
	else {
		$hamstr_hdr = sprintf "%s|%s|%s|%s-%d",
			$first_sector[0],                      # cog id
			$first_sector[-1],                     # reftaxon
			$first_sector[1],                      # taxon
			$first_sector[2],                      # seq id
			$first_sector[3] ne '.' ? abs( eval(  $first_sector[3] ) ) + 1 : 0,  # length
		;
		if (scalar @sectors > 1) {
			for (my $i = 1; $i < scalar @sectors; ++$i) {
				my @sectorfields = split /\|/, $sectors[$i];
				$hamstr_hdr .= sprintf "PP%s-%d",
					$sectorfields[0],                     # seq id
					$sectorfields[1] ne '.' ? abs( eval(  $sectorfields[1] ) ) + 1 : 0, # length
				;
			}
		}
	}
	return $hamstr_hdr;
}

# Module for object-oriented Fasta parsing
package Seqload::Fasta;
use strict;
use warnings;
use Carp;

# Constructor. Returns a sequence database object.
sub open {
  my ($class, $filename) = @_;
  open (my $fh, '<', $filename)
    or confess "Fatal: Could not open $filename\: $!\n";
  my $self = {
    'filename' => $filename,
    'fh'       => $fh
  };
  bless($self, $class);
  return $self;
}

# Returns the next sequence as an array (hdr, seq). 
# Useful for looping through a seq database.
sub next_seq {
  my $self = shift;
  my $fh = $self->{'fh'};
	# this is the trick that makes this work
  local $/ = "\n>"; # change the line separator
  return unless defined(my $item = readline($fh));  # read the line(s)
  chomp $item;
  
  if ($. == 1 and $item !~ /^>/) {  # first line is not a header
    croak "Fatal: " . $self->{'filename'} . " is not a FASTA file: Missing descriptor line\n";
  }

	# remove the '>'
  $item =~ s/^>//;

	# split to a maximum of two items (header, sequence)
  my ($hdr, $seq) = split(/\n/, $item, 2);
	$hdr =~ s/\s+$//;	# remove all trailing whitespace
  $seq =~ s/>//g if defined $seq;
  $seq =~ s/\s+//g if defined $seq; # remove all whitespace, including newlines

  return($hdr, $seq);
}

# Closes the file and undefs the database object.
sub close {
  my $self = shift;
  my $fh = $self->{'fh'};
  my $filename = $self->{'filename'};
  close($fh) or carp("Warning: Could not close $filename\: $!\n");
  undef($self);
}

