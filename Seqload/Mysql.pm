# Documentation before the code#{{{
=head1 NAME

Seqload::Mysql

=head1 DESCRIPTION

A library for handling simple sequences from a MySQL database in an object-oriented fashion.

=head1 SYNOPSIS

  use Seqload::Mysql;

  # open a database connection
  # $table, $hdr_col, $seq_col must be defined separately
	Seqload::Mysql->set_table($table);
	Seqload::Mysql->set_hdr_col($hdr_col);
	Seqload::Mysql->set_seq_col($seq_col);
  my $sequences = Seqload::Mysql->open($db, $dbserver, $dbuser, $dbpwd);
 
  # loop through the results
  while (my ($hdr, $seq) = $sequences->next) {
  	print '>' . $hdr . "\n" . $seq . "\n";
  }
  
  # close the connection
  $sequences->close;

=head1 METHODS

=head2 open(DB, DBSERVER, DBUSER, DBPWD)

Establishes a database connection. Returns a sequence database object.

TABLE, HDR_COL and SEQ_COL must be defined via Seqload::Mysql->set_table(), Seqload::Mysql->set_hdr_col() and Seqload::Mysql->set_seq_col(). 

=head2 next

Returns the next sequence in a sequence database object as an array (HEADER, SEQUENCE).

=head2 close

Closes the database connection.

=head1 CLASS METHODS

=head2 set_table

Sets the table name to use for Seqload::Mysql->open().

=head2 set_hdr_col

Sets the header column to use for Seqload::Mysql->open().

=head2 set_seq_col

Sets the sequence column to use for Seqload::Mysql->open().

=cut#}}}

package Seqload::Mysql;

use Carp;
use DBI;
use DBD::mysql;
use Data::Dumper;

my $table;
my $id_col;
my $date_col;
my $hdr_col;
my $seq_col;

#--------------------------------------------------
# # constructor
#-------------------------------------------------- 
sub open {
	confess "USAGE: Seqload::Mysql->connect(DB, DBSERVER, DBUSER, DBPWD)\n" 
		unless scalar @_ == 5;
	my ($class, $db, $dbserver, $dbuser, $dbpwd) = @_;

	if (!defined $table) { confess "TABLE not defined" }
	elsif (!defined $hdr_col) { confess "HDR_COL not defined\n" }
	elsif (!defined $seq_col) { confess "SEQ_COL not defined\n" }

	my $dbh = DBI->connect("dbi:mysql:$db:$dbserver", $dbuser, $dbpwd)
		or croak "Fatal: Could not connect to database: $!\n";
	
	my $query = "SELECT $hdr_col, $seq_col FROM $table";
	my $sql = $dbh->prepare($query);
	$sql->execute();

	my $self = { 
		'sql' => $sql, 
		'dbh' => $dbh
	};
	bless($self, $class);
	return $self;
}

#--------------------------------------------------
# class methods
#-------------------------------------------------- 

sub set_table {
	my $self = shift;
	if (ref $self) { confess "Class method called as object method\n" }
	unless (scalar @_ == 1) { confess "USAGE: Seqload::Mysql->table(TABLENAME)\n" }
	$table = shift;
}

sub set_hdr_col {
	my $self = shift;
	if (ref $self) { confess "Class method called as object method\n" }
	unless (scalar @_ == 1) { confess "USAGE: Seqload::Mysql->table(HEADER_COL)\n" }
	$hdr_col = shift;
}

sub set_seq_col {
	my $self = shift;
	if (ref $self) { confess "Class method called as object method\n" }
	unless (scalar @_ == 1) { confess "USAGE: Seqload::Mysql->table(SEQUENCE_COL)\n" }
	$seq_col = shift;
}

#--------------------------------------------------
# object methods
#-------------------------------------------------- 

sub next_seq {
	my $self = shift;
	my $sql = $$self{'sql'};
	if (my ($hdr, $seq) = $sql->fetchrow_array) {
		return ($hdr, $seq);
	}
	return;
}

sub close {
	my $self = shift;
	$$self{'dbh'}->disconnect;
	undef $self;
}

sub DESTROY {
	my $self = shift;
	$self->disconnect;
}

# return true
1;
