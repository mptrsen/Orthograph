package IO::Tee;

require 5.004;
use strict;
use Carp;
use Symbol;
use IO::Handle;
use IO::File;
use vars qw($VERSION @ISA);
$VERSION = '0.64';
@ISA = 'IO::Handle';

# Constructor -- bless array reference into our class

sub new
{
    my $class = shift;
    my $self = gensym;
    @{*$self} = map {
        ! ref($_) ? IO::File->new($_)
        : ref($_) eq 'ARRAY' ? IO::File->new(@$_)
        : ref($_) eq 'GLOB' ? bless $_, 'IO::Handle'
        : $_ or return undef } @_;
    bless $self, $class;
    tie *$self, $class, $self;
    return $self;
}

# Return a list of all associated handles

sub handles
{
    @{*{$_[0]}};
}

# Proxy routines for various IO::Handle and IO::File operations

sub _method_return_success
{
    my $method = (caller(1))[3];
    $method =~ s/.*:://;

    my $self = shift;
    my $ret = 1;
    foreach my $fh (@{*$self}) { undef $ret unless $fh->$method(@_) }
    return $ret;
}

sub close        { _method_return_success(@_) }
sub truncate     { _method_return_success(@_) }
sub write        { _method_return_success(@_) }
sub syswrite     { _method_return_success(@_) }
sub format_write { _method_return_success(@_) }
sub fcntl        { _method_return_success(@_) }
sub ioctl        { _method_return_success(@_) }
sub flush        { _method_return_success(@_) }
sub clearerr     { _method_return_success(@_) }
sub seek         { _method_return_success(@_) }

sub formline
{
    my $self = shift;
    my $picture = shift;
    local($^A) = $^A;
    local($\) = "";
    formline($picture, @_);

    my $ret = 1;
    foreach my $fh (@{*$self}) { undef $ret unless print $fh $^A }
    return $ret;
}

sub _state_modify
{
    my $method = (caller(1))[3];
    $method =~ s/.*:://;
    croak "$method values cannot be retrieved collectively" if @_ <= 1;

    my $self = shift;
    if (ref $self)
    {
        foreach my $fh (@{*$self}) { $fh->$method(@_) }
    }
    else
    {
        IO::Handle->$method(@_);
    }
    # Note that we do not return any "previous value" here
}

sub autoflush                    { _state_modify(@_) }
sub output_field_separator       { _state_modify(@_) }
sub output_record_separator      { _state_modify(@_) }
sub format_page_number           { _state_modify(@_) }
sub format_lines_per_page        { _state_modify(@_) }
sub format_lines_left            { _state_modify(@_) }
sub format_name                  { _state_modify(@_) }
sub format_top_name              { _state_modify(@_) }
sub format_line_break_characters { _state_modify(@_) }
sub format_formfeed              { _state_modify(@_) }

sub input_record_separator
{
    my $self = shift;
    my $ret = (ref $self ? ${*$self}[0] : 'IO::Handle')
        ->input_record_separator(@_);
    $ret; # This works around an apparent bug in Perl 5.004_04
}

sub input_line_number
{
    my $self = shift;
    my $ret = ${*$self}[0]->input_line_number(@_);
    $ret; # This works around an apparent bug in Perl 5.004_04
}

# File handle tying interface

sub TIEHANDLE
{
    my ($class, $self) = @_;
    return bless *$self{ARRAY}, $class;
}

sub PRINT
{
    my $self = shift;
    my $ret = 1;
    foreach my $fh (@$self) { undef $ret unless print $fh @_ }
    return $ret;
}

sub PRINTF
{
    my $self = shift;
    my $fmt = shift;
    my $ret = 1;
    foreach my $fh (@$self) { undef $ret unless printf $fh $fmt, @_ }
    return $ret;
}

sub _multiplex_input
{
    my ($self, $input) = @_;
    my $ret = 1;
    if (length $input)
    {
        for (my $i = 1; $i < @$self; ++$i)
        {
            undef $ret unless print {$self->[$i]} $input;
        }
    }
    $ret;
}

sub READ
{
    my $self = shift;
    my $bytes = $self->[0]->read(@_);
    $bytes and $self->_multiplex_input(substr($_[0], $_[2], $bytes));
    $bytes;
}

sub READLINE
{
    my $self = shift;
    my $infh = $self->[0];
    if (wantarray)
    {
        my @data;
        my $data;
        while (defined($data = <$infh>) and length($data))
        {
            push @data, $data;
            $self->_multiplex_input($data);
        }
        @data;
    }
    else
    {
        my $data = <$infh>;
        defined $data and $self->_multiplex_input($data);
        $data;
    }
}

sub GETC
{
    my $self = shift;
    my $data = getc($self->[0]);
    defined $data and $self->_multiplex_input($data);
    $data;
}

sub sysread
{
    my $self = shift;
    my $bytes = ${*$self}[0]->sysread(@_);
    $bytes and (\@{*$self})->
        _multiplex_input(substr($_[0], $_[2] || 0, $bytes));
    $bytes;
}

sub EOF
{
    my $self = shift;
    return $self->[0]->eof;
}

1;
__END__

=head1 NAME

IO::Tee - Multiplex output to multiple output handles

=head1 SYNOPSIS

    use IO::Tee;

    $tee = IO::Tee->new($handle1, $handle2);
    print $tee "foo", "bar";
    my $input = <$tee>;

=head1 DESCRIPTION

C<IO::Tee> objects can be used to multiplex input and output in two
different ways.  The first way is to multiplex output to zero or more
output handles.  The C<IO::Tee> constructor, given a list of output
handles, returns a tied handle that can be written to.  When written
to (using print or printf), the C<IO::Tee> object multiplexes the
output to the list of handles originally passed to the constructor.
As a shortcut, you can also directly pass a string or an array
reference to the constructor, in which case C<IO::File::new> is called
for you with the specified argument or arguments.

The second way is to multiplex input from one input handle to zero or
more output handles as it is being read.  The C<IO::Tee> constructor,
given an input handle followed by a list of output handles, returns a
tied handle that can be read from as well as written to.  When written
to, the C<IO::Tee> object multiplexes the output to all handles passed
to the constructor, as described in the previous paragraph.  When read
from, the C<IO::Tee> object reads from the input handle given as the
first argument to the C<IO::Tee> constructor, then writes any data
read to the output handles given as the remaining arguments to the
constructor.

The C<IO::Tee> class supports certain C<IO::Handle> and C<IO::File>
methods related to input and output.  In particular, the following
methods will iterate themselves over all handles associated with the
C<IO::Tee> object, and return TRUE indicating success if and only if
all associated handles returned TRUE indicating success:

=over 4

=item close

=item truncate

=item write

=item syswrite

=item format_write

=item formline

=item fcntl

=item ioctl

=item flush

=item clearerr

=item seek

=back

The following methods perform input multiplexing as described above:

=over 4

=item read

=item sysread

=item readline

=item getc

=item gets

=item eof

=item getline

=item getlines

=back

The following methods can be used to set (but not retrieve) the
current values of output-related state variables on all associated
handles:

=over 4

=item autoflush

=item output_field_separator

=item output_record_separator

=item format_page_number

=item format_lines_per_page

=item format_lines_left

=item format_name

=item format_top_name

=item format_line_break_characters

=item format_formfeed

=back

The following methods are directly passed on to the input handle given
as the first argument to the C<IO::Tee> constructor:

=over 4

=item input_record_separator

=item input_line_number

=back

Note that the return value of input multiplexing methods (such as
C<print>) is always the return value of the input action, not the
return value of subsequent output actions.  In particular, no error is
indicated by the return value if the input action itself succeeds but
subsequent output multiplexing fails.

=head1 EXAMPLE

    use IO::Tee;
    use IO::File;

    my $tee = new IO::Tee(\*STDOUT,
        new IO::File(">tt1.out"), ">tt2.out");

    print join(' ', $tee->handles), "\n";

    for (1..10) { print $tee $_, "\n" }
    for (1..10) { $tee->print($_, "\n") }
    $tee->flush;

    $tee = new IO::Tee('</etc/passwd', \*STDOUT);
    my @lines = <$tee>;
    print scalar(@lines);

=head1 AUTHOR

Chung-chieh Shan, ken@digitas.harvard.edu

=head1 COPYRIGHT

Copyright (c) 1998-2001 Chung-chieh Shan.  All rights reserved.
This program is free software; you can redistribute it and/or
modify it under the same terms as Perl itself.

=head1 SEE ALSO

L<perlfunc>, L<IO::Handle>, L<IO::File>.

=cut
