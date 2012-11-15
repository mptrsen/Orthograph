package Orthograph::Config;
use strict;
use warnings;
require Exporter;
use Orthograph::Functions;

our @ISA = qw( Exporter );
our @EXPORT_OK = qw( $config $configfile );

my $program_name = 'Orthograph';
my $configfile = lc($program_name) . '.conf';#{{{
our $config;
$configfile = Orthograph::Functions::get_configfile($configfile);
if (-e $configfile) {
  $config = Orthograph::Functions::parse_config($configfile);
}#}}}

1;
