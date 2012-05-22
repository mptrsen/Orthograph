package Forage::Orthoset;
use strict;
use warnings;
use Carp;
use IO::File;
use IO::Dir;
use File::Basename;
use Seqload::Fasta;
use Data::Dumper;

my $verbose = 0;
my $debug = 0;
my $setdir = '';
my $hmmlist = [ ];
my $alignment_program = 'mafft --localpair --maxiterate 1000';
my $hmmbuild_program = 'hmmbuild';
my $stockholm_header_width = 50;

sub new {
	my $class = shift(@_);
	my $setname = shift(@_);

	my $self = {
		'setname' => $setname,
		'setdir'  => $setdir,
		'hmmlist' => $hmmlist;
	}
	bless($self, $class);
	return($self);
}

=head1 Class methods

=head2 setdir()

Sets the ortholog set directory.
Argument: Scalar string pathname
Returns: -
=cut

sub setdir {
	my $self = shift(@_);
	if (ref($self)) { confess('Class method called as object method') }
	unless (scalar(@_) == 1) { confess('Usage: Forage::Orthoset->setdir($dirname)') }
	$setdir = shift(@_);
	unless (-e $setdir) { confess("Fatal: set dir $setdir does not exist") }
}

=head2 verbose()

Sets or returns verbosity status
Argument: Scalar true/false or nothing
Returns: -
=cut

sub verbose {
	my $self = shift(@_);
	if (ref($self)) { confess('Class method called as object method') }
	if (scalar(@_) == 0) { return $verbose }
	unless (scalar(@_) == 1) { confess('Usage: Forage::Orthoset->verbose(1|0)') }
	$verbose = shift(@_);
}

=head2 debug()

Sets or returns debug output
Argument: Scalar true/false or nothing
Returns: -
=cut

sub debug {
	my $self = shift(@_);
	if (ref($self)) { confess('Class method called as object method') }
	if (scalar(@_) == 0) { return $debug }
	unless (scalar(@_) == 1) { confess('Usage: Forage::Orthoset->debug(1|0)') }
	$debug = shift(@_);
}

=head2 stockholm_header_width()

Sets or returns the stockholm format header width.
Arguments: Scalar integer or nothing
Returns: Scalar integer stockholm format header width or false if not provided an int

=cut

sub stockholm_header_width {
	my $self = shift(@_);
	if (ref($self)) { confess('Class method called as object method') }
	if (scalar(@_) == 0) { return($stockholm_header_width) }
	unless (scalar(@_) == 1) { confess('Usage: Forage::Orthoset->stockholm_header_width($integer)') }
	$stockholm_header_width = shift(@_);
	unless ($stockholm_header_width =~ /^\d+$/) { confess('Forage::Orthoset->stockholm_header_width: Not an integer') }
}

=head1 Object methods

=head2 hmmlist

Returns a list of HMM files in the set directory.

=cut

sub hmmlist {
	my $self = shift(@_);
	my $dirh = IO::Dir->new($self->setdir);
	while (my $file = $dirh->read()) {
		next unless $file =~ /\.hmm$/;
		push(@{$self->{'hmmlist'}}, File::Spec->catfile($self->setdir(), $file));
	}
	return $self->{'hmmlist'};
}

=head2 make_hmm

Generates a HMM from a fasta file

=cut

sub make_hmm {
	my $self = shift(@_);
	my $fafile = shift(@_);
	my $alnfile = File::Spec->catfile($fafile . '.aln');
	my $stockholmfile = File::Spec->catfile($alnfile . '.stockh');
	my $hmmfile = File::Spec->catfile($setdir, basename($fafile, '.fa') . '.hmm');
	# TODO align
	my @align_cmd = qq($alignment_program $fafile);	# prints to stdout, must be catched
	my $alnfh = IO::File->new($alnfile, 'w');
	print $alnfh `@align_cmd`
		or croak("Fatal: Could not write to alignment file $alnfile\: $!");
	$alnf->close();
	# TODO convert to stockholm
	$alnfh = Seqload::Fasta->open($alnfile);
	my $stockhfh = IO::File->new($stockholmfile, 'w');
	print $stockhfh "# STOCKHOLM 1.0\n";	# stockholm header
	while (my ($h, $s) = $alnfh->next_seq()) {
		$h =~ s/\s+/_/g;	# replace all whitespace with underscore
		printf($stockhfh "%-50s %s\n", $h, $s)
			or croak("Fatal: Could not write to file $stockholmfile\: $!\n");
	}
	print $stockhfh "//";	# stockholm footer
	$stockhfh->close();

	# TODO hmmbuild
	my @hmmbuild_cmd = qq($hmmbuild_program $hmmfile $alnfile);
	system(@hmmbuild_cmd) 
		and croak("Fatal: $hmmbuild failed: $!\n");
}

return 1;
