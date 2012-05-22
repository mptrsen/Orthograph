# Documentation before the code#{{{
=head1 NAME

Seqload::Fasta

=head1 DESCRIPTION

A library for handling FASTA sequences in an object-oriented fashion. 

=head1 SYNOPSIS

  use Seqload::Fasta qw(fasta2csv);
  
  # open the file, return fasta file object
  my $file = Seqload::Fasta->open($filename);
  
  # loop through the sequences
  while (my ($hdr, $seq) = $file->next_seq) {
    print $hdr . "\n" . $seq . "\n";
  }

  # close the file object
  $file->close;

  # convert a fasta file to a csv file
  fasta2csv($fastafile, $csvfile);


=head1 METHODS

=head2 open(FILENAME)

Opens a fasta file. Returns a sequence database object.

=head2 next_seq

Returns the next sequence in a sequence database object as an array (HEADER,
SEQUENCE). Note that the '>' character is truncated from the header.

=head2 close

Closes the file.

=head1 FUNCTIONS

=head2 fasta2csv($fastafile, $csvfile)

Converts a fasta file into a csv file where each line consists of
'HEADER,SEQUENCE'. Manages opening, parsing and closing of the files, no
additional file handles necessary..

=cut#}}}

package Seqload::Fasta;
use Carp;
require Exporter;
@ISA = qw(Exporter);
@EXPORT_OK = qw(fasta2csv check_if_fasta);

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
  local $/ = "\n>"; # change the line separator
  return unless $item = readline($fh);  # read the line(s)
  chomp $item;
  
  if ($. == 1 and $item !~ /^>/) {  # first line is not a header
    croak "Fatal: " . $self->{'filename'} . "is not a FASTA file: Missing descriptor line\n";
  }

  $item =~ s/^>//;

  my ($hdr, $seq) = split(/\n/, $item, 2);
  $seq =~ s/>//g if defined $seq;
  $seq =~ s/\s+//g if defined $seq; # remove all whitespace, including newlines

  return($hdr, $seq);
}

# Destructor. Closes the file and undefs the database object.
sub close {
  my $self = shift;
  my $fh = $self->{'fh'};
  my $filename = $self->{'filename'};
  close($fh) or croak "Fatal: Could not close $filename\: $!\n";
  undef($self);
}

# I dunno if this is required but I guess this is called when you undef() an object
sub DESTROY {
  my $self = shift;
  $self->close;
}

# Convert a fasta file to a csv file the easy way
# Usage: Seqload::Fasta::fasta2csv($fastafile, $csvfile);
sub fasta2csv {
  my $fastafile = shift;
  my $csvfile = shift;
  my $fastafh = Seqload::Fasta->open($fastafile);

  open(my $outfh, '>', $csvfile)
    or confess "Fatal: Could not open $outfile\: $!\n";
  while ((my $hdr, $seq) = $fastafh->next_seq) {
    print $outfh $hdr . ',' . $seq . "\n";
  }
  close $outfh;

  $fastafh->close;
  return 1;
}


# validates a fasta file by looking at the FIRST (header, sequence) pair
# arguments: scalar string path to file
# returns: true on validation, false otherwise
sub check_if_fasta {
	my $infile = shift;
	my $infh = Seqload::Fasta->open($infile);
	my ($h, $s) = $infh->next_seq() or return 0;
	return 1;
}
# return true
1;
