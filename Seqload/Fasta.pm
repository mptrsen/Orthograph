# Documentation before the code#{{{
=head1 NAME

Seqload::Fasta

=head1 DESCRIPTION

A library for handling FASTA sequences in an object-oriented fashion. 

=head1 SYNOPSIS

  use Seqload::Fasta;
  
  my $file = Seqload::Fasta->new($filename);
  
  while (my ($hdr, $seq) = $file->next) {
    print $hdr . "\n" . $seq . "\n";
  }

	$file->close;

=head1 METHODS

=head2 open(FILENAME)

Opens a fasta file. Returns a sequence database object.

=head2 next

Returns the next sequence in a sequence database object as an array (HEADER, SEQUENCE).

=head2 close

Closes the file.

=cut#}}}

package Seqload::Fasta;
use Carp;

# Constructor. Returns a sequence database object.
sub open {
	my ($class, $filename) = @_;
	open (my $fh, '<', $filename)
		or croak "Fatal: Could not open $filename\: $!\n";
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
	local $/ = "\n>";	# change the line separator
	return unless $item = readline($fh);	# read the line(s)
	chomp $item;
	
	if ($. == 1 and $item !~ /^>/) {	# first line is not a header
		croak "Fatal: " . $self->{'filename'} . "is not a FASTA file: Missing descriptor line\n";
	}

	$item =~ s/^>//;

	my ($hdr, $seq) = split(/\n/, $item, 2);
	$seq =~ s/>//g if defined $seq;
	$seq =~ s/\s+//g if defined $seq;	# remove all whitespace, including newlines

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

# return true
1;
