#!/usr/bin/perl
#--------------------------------------------------
# This file is part of Forage.
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
use Data::Dumper;
use Config;   # allows checking for system configuration
use Getopt::Long;
use Carp;   # alternative warn and die
use File::Path qw(make_path); # mkdir with parent dirs; this also uses File::Spec
use File::Basename; # parsing path names
use File::Temp; # temporary files
use IO::File; # object-oriented access to files
use IO::Dir;  # object-oriented access to dirs
<<<<<<< HEAD
=======
use Tie::File;
use Forage::Hmmsearch;
>>>>>>> ffbb25968d9e0cdc1dac57195419dd04e634aa20
use Seqload::Fasta qw(fasta2csv); # object-oriented access to fasta files, fasta2csv converter
use Seqload::Mysql; # object-oriented access to fasta-style MySQL databases
(my $libdir = $0) =~ s/forage\.pl$//; 
use lib qw($libdir);

my $version = 0.00002;
print "Forage: Find Orthologs using Reciprocity Among Genes and ESTs\n";
print "Copyright 2012 Malte Petersen <mptrsen\@uni-bonn.de>\n";
print "Version $version\n\n";

#--------------------------------------------------
# Only use threads if the system supports it.
# The whole threads system is totally not implemented yet,
# do not attempt to use it!
#-------------------------------------------------- 
my $use_threads = 0;#{{{
if ($use_threads == 1 and $Config{'useithreads'}) {
  print "Using threads.\n";
  use Forage::Threaded;
}
elsif ($use_threads == 1 and !$Config{'useithreads'}) {
  die "Fatal: Cannot use threads: Your version of Perl was not compiled with threading support. Not using threads.\n";
}
else {
<<<<<<< HEAD
  print "Not using threads. ";
  use Forage::Unthreaded;
  print "Loaded Forage::Unthreaded.\n";
  $use_threads = 0;
=======
	print "Not using threads.\n";
	$use_threads = 0;
>>>>>>> ffbb25968d9e0cdc1dac57195419dd04e634aa20
}#}}}

#--------------------------------------------------
# # Variable initialisation
#-------------------------------------------------- 
my $config;                       # will hold the configuration from config file

#--------------------------------------------------
# # Parse config file
#-------------------------------------------------- 
(my $configfile = $0) =~ s/(\.pl)?$/.conf/;#{{{
# mini argument parser for the configfile
for (my $i = 0; $i < scalar @ARGV; ++$i) {
<<<<<<< HEAD
  if ($ARGV[$i] =~ /-c/) {
    if ($ARGV[$i+1] !~ /^-/) {
      $configfile = $ARGV[$i+1];
    }
    else { warn "Warning: Config file name '$ARGV[$i+1]' not accepted (use './$ARGV[$i+1]' if you mean it). Falling back to '$configfile'\n" }
  }
=======
	if ($ARGV[$i] =~ /-c/) {
		if ($ARGV[$i+1] !~ /^-/) {
		  $configfile = $ARGV[$i+1];
		}
		else { warn "Warning: Config file name '$ARGV[$i+1]' not accepted (use './$ARGV[$i+1]' if you mean it). Falling back to '$configfile'\n" }
	}
>>>>>>> ffbb25968d9e0cdc1dac57195419dd04e634aa20
}
if (-e $configfile) {
  print "Parsing config file '$configfile'.\n";
  $config = &parse_config($configfile);
}#}}}

#--------------------------------------------------
# # Programs in the order of their use
#-------------------------------------------------- 
my $translateprog = $config->{'translate_program'} ? $config->{'translate_program'} : 'fastatranslate';#{{{
my $hmmsearchprog = $config->{'hmmsearch_program'} ? $config->{'hmmsearch_program'} : 'hmmsearch';
my $blastprog     = $config->{'blast_program'}     ? $config->{'blast_program'}     : 'blastp';#}}}

#--------------------------------------------------
# # These variables can be set in the config file
#-------------------------------------------------- 
my $blast_max_hits = $config->{'blast_max_hits'}               ? $config->{'blast_max_hits'}              : 3;#{{{
my $backup_ext     = $config->{'backup_extension'}      ? $config->{'backup_extension'}     : '.bak';
my $blastoutdir    = $config->{'blastoutdir'}           ? $config->{'blastoutdir'}          : 'blastp';
my $estfile        = $config->{'estfile'}               ? $config->{'estfile'}              : '';
my $eval_threshold = $config->{'eval_threshold'}        ? $config->{'eval_threshold'}       : undef;
my $hmmdir         = $config->{'hmmdir'}                ? $config->{'hmmdir'}               : '';
my $hmmfile        = $config->{'hmmfile'}               ? $config->{'hmmfile'}              : '';
my $hmmoutdir      = $config->{'hmmsearch_output_dir'}  ? $config->{'hmmsearch_output_dir'} : basename($hmmsearchprog);
my $mysql_dbname   = $config->{'mysql_dbname'}          ? $config->{'mysql_dbname'}         : 'forage';
my $mysql_dbpwd    = $config->{'mysql_dbpassword'}      ? $config->{'mysql_dbpassword'}     : 'root';
my $mysql_dbserver = $config->{'mysql_dbserver'}        ? $config->{'mysql_dbserver'}       : 'localhost';
my $mysql_dbuser   = $config->{'mysql_dbuser'}          ? $config->{'mysql_dbuser'}         : 'root';
my $mysql_table_blast = $config->{'mysql_table_blast'}  ? $config->{'mysql_table_blast'}    : 'blast';
my $mysql_table_core_orthologs = $config->{'mysql_table_core_orthologs'}    ? $config->{'mysql_table_core_orthologs'}     : 'core_orthologs';
my $mysql_table_ests = $config->{'mysql_table_ests'}    ? $config->{'mysql_table_ests'}     : 'ests';
my $mysql_table_hmmsearch = $config->{'mysql_table_hmmsearch'} ? $config->{'mysql_table_hmmsearch'} : 'hmmsearch';
my $outdir         = $config->{'output_directory'}      ? $config->{'output_directory'}     : undef;
my $quiet          = $config->{'quiet'}                 ? $config->{'quiet'}                : undef;  # I like my quiet
my $score_threshold = $config->{'score_threshold'}      ? $config->{'score_threshold'}      : undef;
my $species_name   = $config->{'species_name'}          ? $config->{'species_name'}         : undef;
my $verbose        = $config->{'verbose'}               ? $config->{'verbose'}              : undef;

#}}}

#--------------------------------------------------
# # More variables
#-------------------------------------------------- 
my $count               = 0;#{{{
my $hmmresultfileref;
my $hmmfullout          = 0;
my $hitcount;
my $hmmhitcount;
my $mysql_dbi           = "dbi\:mysql\:$mysql_dbname\:$mysql_dbserver";
my $mysql_col_date      = 'date';
my $mysql_col_eval      = 'eval';
my $mysql_col_hdr       = 'hdr';
my $mysql_col_hmm       = 'hmm';
my $mysql_col_hmmhit    = 'hmmhit';
my $mysql_col_id        = 'id';
my $mysql_col_query     = 'query';
my $mysql_col_score     = 'score';
my $mysql_col_seq       = 'seq';
my $mysql_col_spec      = 'spec';
my $mysql_col_target    = 'target';
my $mysql_col_taxon     = 'taxon';
my $num_ests;
my $preparedb;
my $protfile            = '';
my $timestamp           = time();
my @eval_option;
my @hmmfiles;
my @hmmsearchcmd;
my @score_option;
my @seqobjs;
my $i;#}}}

#--------------------------------------------------
# # Get command line options. These may override variables set via the config file.
#-------------------------------------------------- 
GetOptions( 'v'     => \$verbose,#{{{
<<<<<<< HEAD
  'c'               => \$configfile,
  'threads'         => \$use_threads,   # make using threads optional
  'estfile=s'       => \$estfile,
  'E=s'             => \$estfile,
  'eval=s'          => \$eval_threshold,
  'score=s'         => \$score_threshold,
  'H=s'             => \$hmmfile,
  'hmmdir=s'        => \$hmmdir,
  'hmmsearchprog=s' => \$hmmsearchprog,
  'hmmfullout'      => \$hmmfullout,
  'preparedb'       => \$preparedb,
  'quiet'           => \$quiet,
  'taxon=s'         => \$species_name,
=======
	'c'               => \$configfile,
	'threads'         => \$use_threads,   # make using threads optional
	'estfile=s'       => \$estfile,
	'E=s'             => \$estfile,
	'eval=s'          => \$eval_threshold,
	'score=s'         => \$score_threshold,
	'H=s'             => \$hmmfile,
	'hmmdir=s'        => \$hmmdir,
	'hmmsearchprog=s' => \$hmmsearchprog,
	'hmmfullout'      => \$hmmfullout,
	'preparedb'       => \$preparedb,
	'quiet'           => \$quiet,
	'taxon=s'         => \$species_name,
>>>>>>> ffbb25968d9e0cdc1dac57195419dd04e634aa20
);#}}}

#--------------------------------------------------
# # Special case: Prepare the MySQL database by dropping and recreating all
# # tables if requested, then exit
#-------------------------------------------------- 
if ($preparedb) {#{{{
  print "Setting MySQL database $mysql_dbname to a clean slate...\n";
  &preparedb;
  print "OK; MySQL database now ready to run Forage.\n";
  exit;
}#}}}

#--------------------------------------------------
# # Normal run. Input error checking, reporting etc
#-------------------------------------------------- 
&intro;

#--------------------------------------------------
# # create list of HMM files
#-------------------------------------------------- 
@hmmfiles = &hmmlist;

unless ($quiet) {#{{{
  print "Using HMM dir $hmmdir with ", scalar @hmmfiles, " HMMs\n" 
    if $hmmdir;
  print "Using HMM file $hmmfile.\n" 
    if $hmmfile;
  print "e-Value cutoff: $eval_threshold.\n" 
    if $eval_threshold;
  print "Score cutoff: $score_threshold.\n"
    if $score_threshold;
}#}}}

#--------------------------------------------------
# # translate the ESTs to protein, feed that shit to the database
#-------------------------------------------------- 
$| = 1;
$protfile = &translate_est(File::Spec->catfile($estfile));#{{{

# Clear database of data from the same species
print "Clearing database of previous results from '$species_name'... " unless $quiet;
&clear_db;
print "done.\n" unless $quiet;

print "Storing translated sequences to MySQL database '$mysql_dbname' on $mysql_dbserver... " unless $quiet;

# Create temporary csv file for high-speed reading into database
my $tmpfh = File::Temp->new('UNLINK' => 1);
fasta2csv($protfile, $tmpfh) or die "Could not fasta2csv $protfile into $tmpfh\: $!\n";

# load data from csv file into database
my $query = "LOAD DATA LOCAL INFILE '$tmpfh' INTO TABLE $mysql_table_ests FIELDS TERMINATED BY ',' ($mysql_col_hdr, $mysql_col_seq)";
# add date and spec
my $query_update = "UPDATE $mysql_table_ests SET 
<<<<<<< HEAD
  $mysql_col_date='" . time() . "', 
  $mysql_col_spec='$species_name' WHERE $mysql_col_date IS NULL";
=======
	$mysql_col_date='" . time() . "', 
	$mysql_col_spec='$species_name' WHERE $mysql_col_date IS NULL";
>>>>>>> ffbb25968d9e0cdc1dac57195419dd04e634aa20

# open connection
my $dbh = DBI->connect($mysql_dbi, $mysql_dbuser, $mysql_dbpwd);
$num_ests = $dbh->do($query);
$dbh->do($query_update);
# disconnect ASAP
$dbh->disconnect;#}}}

# report
print " done.\n" unless $quiet;
$| = 0;

printf "%d sequences stored to database '%s' on %s.\n",
  $num_ests,
  $mysql_dbname,
  $mysql_dbserver unless $quiet;

#--------------------------------------------------
# # Setup the Forage module
#-------------------------------------------------- 

# These are all class methods#{{{

# verbose output; this is a class method
Forage::Hmmsearch->verbose(1) if $verbose;

# the output directory
Forage::Hmmsearch->hmmoutdir($hmmoutdir);

# the hmmsearch program
# Forage::Hmmsearch->hmmsearchprog($hmmsearchprog);

# the hmmsearch command
Forage::Hmmsearch->hmmsearchcmd(\@hmmsearchcmd);

# whether or not we want full output
Forage::Hmmsearch->hmmfullout(0);#}}}

#--------------------------------------------------
# # HMMsearch the protfile using all HMMs
#-------------------------------------------------- 
print "Hmmsearching the protein file using all HMMs in $hmmdir\:\n" unless $quiet;
$i = 0;
$hitcount = 0;

# SQL query for pushing HMMsearch results to the db
my $query_insert_hmmresult = "INSERT INTO $mysql_table_hmmsearch (
  $mysql_col_spec,
  $mysql_col_hmmhit,
  $mysql_col_hmm,
  $mysql_col_score,
  $mysql_col_eval) VALUES (
  ?,
  ?,
  ?,
  ?,
  ?)";

#--------------------------------------------------
# # Do the HMM search - this may be pipelined in the future
#-------------------------------------------------- 

HMMFILE:
foreach my $hmmfile (@hmmfiles) {#{{{
<<<<<<< HEAD
  ++$i;
  # create new hmmobject with a hmm file, should have all the necessary info for doing hmmsearch
  my $hmmobj = Forage::Unthreaded->new($hmmfile); 
  # now do the hmmsearch on the protfile
  $hmmobj->hmmsearch($protfile);
  # count the hmmsearch hits
  unless ($hmmobj->hmmhitcount()) { # do not care further with HMM files that did not return any result
    printf "%4d hits detected for %s.\n", 0, basename($hmmobj->hmmfile) unless $quiet;
    next;
  }
  printf "%4d hits detected for %s\n", $hmmobj->hmmhitcount, basename($hmmobj->hmmfile) unless $quiet;

  # push results to database
  $dbh = DBI->connect($mysql_dbi, $mysql_dbuser, $mysql_dbpwd);
  my $sql = $dbh->prepare($query_insert_hmmresult);

  # this is an array reference
  foreach my $hmmhit (@{$hmmobj->hmmhits_arrayref}) {
    $sql->execute(
      $species_name,
      $hmmhit->[0], # target (header)
      $hmmhit->[1], # query (HMM)
      $hmmhit->[2], # score
      $hmmhit->[3]  # evalue
    ) or die "Fatal: Could not push to database!\n";
    ++$hmmhitcount;
  }
  $dbh->disconnect;
  print "       ... pushed to database.\n" if $verbose;;
  ++$hitcount;
  

  #--------------------------------------------------
  # #TODO next:
  #-------------------------------------------------- 
  # Get the hits from the database, store to a tempfile, blast against the core ortholog db
  # using either blast or hmmsearch [hmmsearch: don't we need profiles for that?]
  my $query = "SELECT * 
    FROM $mysql_table_hmmsearch 
    WHERE $mysql_col_hmm = " . $hmmobj->hmmhits_arrayref->[0][1] . " 
    ORDER BY $mysql_col_eval";

  $dbh = DBI->connect($mysql_dbi, $mysql_dbuser, $mysql_dbpwd);
  $sql = $dbh->prepare($query);
  $sql->execute();
  printf "     %-60s %-7s%-5s\n", 'Target', 'Score', 'eval';
  while (my @result = $sql->fetchrow_array()) {
    printf "     %-60s %-7.2f%-4.1e\n", $result[3], $result[5], $result[4];
  }
  $dbh->disconnect;


  #--------------------------------------------------
  # # TODO then:
  #-------------------------------------------------- 
  # for the re-hits, gather nuc seq and compile everything that Karen wants output :)
=======
	++$i;
	# create new hmmobject with a hmm file, should have all the necessary info for doing hmmsearch
	my $hmmobj = Forage::Hmmsearch->new($hmmfile); 
	# now do the hmmsearch on the protfile
	$hmmobj->hmmsearch($protfile);
	# count the hmmsearch hits
	unless ($hmmobj->hmmhitcount()) { # do not care further with HMM files that did not return any result
		printf "%4d hits detected for %s.\n", 0, basename($hmmobj->hmmfile) unless $quiet;
		next;
	}
	printf "%4d hits detected for %s:\n", $hmmobj->hmmhitcount, basename($hmmobj->hmmfile) unless $quiet;
	# print list of hits if verbose
	if ($verbose) {
		printf "     %s\n", $_->{'target'} foreach (@{$hmmobj->hmmhits_arrayref});
	}

	#--------------------------------------------------
	# # push results to database
	#-------------------------------------------------- 
	$dbh = DBI->connect($mysql_dbi, $mysql_dbuser, $mysql_dbpwd);
	my $sql = $dbh->prepare($query_insert_hmmresult);

	# this is a reference to an array of hashes
	foreach my $hmmhit (@{$hmmobj->hmmhits_arrayref}) {
		$sql->execute(
		  $species_name,
		  $hmmhit->{'target'}, # target (header)
		  $hmmhit->{'query'},  # query (HMM)
		  $hmmhit->{'score'},  # score
		  $hmmhit->{'eval'}    # evalue
		) or die "Fatal: Could not push to database!\n";
		++$hmmhitcount;
	}
	$dbh->disconnect;
	print "       ... pushed to database.\n" if $verbose;;
	++$hitcount;
	

	#--------------------------------------------------
	# # the reciprocal search
	#-------------------------------------------------- 

	# setup SQL query; use the first array item since they all share the query (HMM) ID 
	my $query_get_sequences = "SELECT $mysql_col_hdr, $mysql_col_seq
		FROM $mysql_table_ests INNER JOIN $mysql_table_hmmsearch
		ON $mysql_table_hmmsearch.$mysql_col_hmmhit = $mysql_table_ests.$mysql_col_hdr 
		WHERE $mysql_table_hmmsearch.$mysql_col_hmm = " . $hmmobj->hmmhits_arrayref->[0]{'query'};

  # get the sequences from the database (as array->array reference)
  $count = 0;
  foreach my $result (@{&mysql_get($query_get_sequences)}) {
    ++$count;
    # run blastp on it: setup a temporary file that will not be preserved
    my $tmpfh = File::Temp->new('UNLINK' => 1);
    # fasta header and sequence
    print $tmpfh '>' . $hmmobj->hmmhits_arrayref->[0]{'query'} . ':' . $$result[0] . "\n"; 
    print $tmpfh $$result[1] . "\n";

    my $blastoutfile = File::Spec->catfile($blastoutdir, $hmmobj->hmmhits_arrayref->[0]{'query'} . '_' . $count . '.out');
    my $blastp_command = "blastp -db /home/mpetersen/hamstr/blast_dir/dappu_1391/dappu_1391_prot -query $tmpfh -outfmt 7 -max_target_seqs 3 -out $blastoutfile";
    print "      running $blastp_command\n" if $verbose;
    if (system($blastp_command)) { die "blastp for result $count failed: $!\n" }
  }


	#--------------------------------------------------
	# # TODO then:
	#-------------------------------------------------- 
	# for the re-hits, gather nuc seq and compile everything that Karen wants output :)
>>>>>>> ffbb25968d9e0cdc1dac57195419dd04e634aa20
}#}}}

# report, end the program
printf "%d HMMs hit something.   %d HMM files processed. \n", $hitcount, $i unless $quiet;
printf "Forage analysis complete. Searched %d EST sequences using %d HMMs and obtained %d hits.\n", $num_ests, $i, $hmmhitcount;
exit;


###################################################
# # Functions follow
###################################################

# Sub: mysql_get
# Get from the database the result of a SQL query
# Expects: QUERY as a string literal
# Returns: Reference to array of arrays (result lines->fields)
sub mysql_get {#{{{
	my $query = shift;
	unless ($query) {
		croak "Usage: mysql_get(QUERY)\n";
	}
  # prepare anonymous array
	my $results = [ ];
  # connect and fetch stuff
	my $dbh = DBI->connect($mysql_dbi, $mysql_dbuser, $mysql_dbpwd);
	my $sql = $dbh->prepare($query);
	$sql->execute();
	while (my @result = $sql->fetchrow_array() ) {
		push(@$results, \@result);
	}
	$sql->finish();
	$dbh->disconnect; # disconnect ASAP
	return $results;
}#}}}

# Sub: parse_config
# Parse a simple, ini-style config file where keys are separated from values by '='.
# E.g.
# outputdir = /home/foo/bar
# translateprog=/usr/bin/translate
sub parse_config {#{{{
<<<<<<< HEAD
  my $file = shift;
  my $conf = { };
  open(my $fh, '<', $file) or die "Fatal: Could not open config file $file\: $!\n";

  while (my $line = <$fh>) {
    next if $line =~ /^\s*$/; # skip empty lines
    next if $line =~ /^\s*#/; # skip comment lines starting with '#'
    
    # split by '=' producing a maximum of two items
    my ($key, $val) = split('=', $line, 2);

    foreach ($key, $val) {
      s/\s+$//; # remove all trailing whitespace
      s/^\s+//; # remove all leading whitespace
    }

    die "Fatal: Configuration option '$key' defined twice in line $. of config file $file\n"
      if defined $conf->{$key};
    $conf->{$key} = $val;
  }
  close($fh);
  return $conf;
=======
	my $file = shift;
	my $conf = { };
	open(my $fh, '<', $file) or die "Fatal: Could not open config file $file\: $!\n";

	while (my $line = <$fh>) {
		next if $line =~ /^\s*$/; # skip empty lines
		next if $line =~ /^\s*#/; # skip comment lines starting with '#'
		
		# split by '=' producing a maximum of two items
		my ($key, $val) = split('=', $line, 2);

		foreach ($key, $val) {
		  s/\s+$//; # remove all trailing whitespace
		  s/^\s+//; # remove all leading whitespace
		}

		die "Fatal: Configuration option '$key' defined twice in line $. of config file $file\n"
		  if defined $conf->{$key};
		$conf->{$key} = $val;
	}
	close($fh);
	return $conf;
>>>>>>> ffbb25968d9e0cdc1dac57195419dd04e634aa20
}#}}}

# Sub: intro
# Checks input, file/dir presence, etc.
# Returns True if everything is OK.
sub intro {#{{{
<<<<<<< HEAD
  die "Fatal: At least two arguments required: EST file (-E) and HMM file/dir (-H or -hmmdir)!\n"
    unless ($estfile and ($hmmdir or $hmmfile));
  
  die "Fatal: Species name needed (-taxon NAME)!\n"
    unless ($species_name);

  # mutually exclusive options
  die "Fatal: Can't operate in both verbose and quiet mode\n"
    if ($verbose and $quiet);

  # mutually exclusive options
  die "Fatal: Can't use both e-value and score thresholds\n"
    if ($eval_threshold and $score_threshold);

  # construct output directory paths
  # outdir may be defined in the config file
  $outdir = defined($outdir) ? $outdir : $species_name;
  $hmmoutdir = defined($hmmoutdir) ? File::Spec->catdir($outdir, $hmmoutdir) : File::Spec->catdir($outdir, basename($hmmsearchprog));

  # build hmmsearch command line
  # full output only if desired, table otherwise
  if ($eval_threshold) {
    @eval_option = ('-E', $eval_threshold);
  }
  if ($score_threshold) {
    @score_option = ('-T', $score_threshold);
  }
  if ($hmmfullout) {
    @hmmsearchcmd = ($hmmsearchprog, @eval_option, @score_option, '-o');
  }
  else {
    @hmmsearchcmd = ($hmmsearchprog, @eval_option, @score_option, '--tblout');
  }

  # the HMM directory
  if (-d $hmmdir) {
    print "HMM dir $hmmdir exists.\n" unless $quiet;
  }
  else {
    die "Fatal: HMM dir $hmmdir does not exist or is not a directory!\n";
  }

  # the EST file
  if (-e $estfile) {
    print "EST file $estfile exists.\n" unless $quiet;
  }
  else {
    die "Fatal: EST file $estfile does not exist!\n";
  }

  # the HMMsearch output directory
  if (-d $hmmoutdir) {
    print "HMMsearch output dir $hmmoutdir exists.\n" unless $quiet;
  }
  else {
    $| = 1;
    print "HMMsearch output dir $hmmoutdir does not exist, creating... " unless $quiet;
    &createdir($hmmoutdir) and print "done.\n";
=======
	die "Fatal: At least two arguments required: EST file (-E) and HMM file/dir (-H or -hmmdir)!\n"
		unless ($estfile and ($hmmdir or $hmmfile));
	
	die "Fatal: Species name needed (-taxon NAME)!\n"
		unless ($species_name);

	# mutually exclusive options
	die "Fatal: Can't operate in both verbose and quiet mode\n"
		if ($verbose and $quiet);

	# mutually exclusive options
	die "Fatal: Can't use both e-value and score thresholds\n"
		if ($eval_threshold and $score_threshold);

	# construct output directory paths
	# outdir may be defined in the config file
	$outdir = defined($outdir) ? $outdir : $species_name;
	$hmmoutdir = defined($hmmoutdir) ? File::Spec->catdir($outdir, $hmmoutdir) : File::Spec->catdir($outdir, basename($hmmsearchprog));
  $blastoutdir = defined($blastoutdir) ? File::Spec->catdir($outdir, $blastoutdir) : File::Spec->catdir($outdir, basename($blastprog));

	# build hmmsearch command line
	# full output only if desired, table otherwise
	if ($eval_threshold) {
		@eval_option = ('-E', $eval_threshold);
	}
	if ($score_threshold) {
		@score_option = ('-T', $score_threshold);
	}
	if ($hmmfullout) {
		@hmmsearchcmd = ($hmmsearchprog, @eval_option, @score_option, '-o');
	}
	else {
		@hmmsearchcmd = ($hmmsearchprog, @eval_option, @score_option, '--tblout');
	}

	# the HMM directory
	if (-d $hmmdir) {
		print "HMM dir $hmmdir exists.\n" unless $quiet;
	}
	else {
		die "Fatal: HMM dir $hmmdir does not exist or is not a directory!\n";
	}

	# the EST file
	if (-e $estfile) {
		print "EST file $estfile exists.\n" unless $quiet;
	}
	else {
		die "Fatal: EST file $estfile does not exist!\n";
	}

	# the HMMsearch output directory
	if (-d $hmmoutdir) {
		print "HMMsearch output dir $hmmoutdir exists.\n" unless $quiet;
	}
	else {
		$| = 1;
		print "HMMsearch output dir $hmmoutdir does not exist, creating... " unless $quiet;
		&createdir($hmmoutdir) and print "done.\n";
		$| = 0;
	}

  # the BLASTP output directory
  if (-d $blastoutdir) {
    print "BLASTP output dir $blastoutdir exists.\n" unless $quiet;
  }
  else {
    $| = 1;
    print "BLASTP output dir $blastoutdir does not exist, creating... " unless $quiet;
    &createdir($blastoutdir) and print "done.\n";
>>>>>>> ffbb25968d9e0cdc1dac57195419dd04e634aa20
    $| = 0;
  }
}#}}}

# Sub: hmmlist
# Expects: nothing (?)
# Returns: array hmmfiles
sub hmmlist {#{{{
<<<<<<< HEAD
  if ($hmmdir) {
    my $dir = IO::Dir->new(File::Spec->catdir($hmmdir));
    while (my $file = $dir->read) {
      push(@hmmfiles, File::Spec->catfile($hmmdir, $file)) if ($file =~ /\.hmm$/);
    }
    $dir->close;
    return sort @hmmfiles;
  }
  else {
    print "single hmm\n"
      if $verbose;
    push(@hmmfiles, $hmmfile);
    return @hmmfiles;
  }
=======
	if ($hmmdir) {
		my $dir = IO::Dir->new(File::Spec->catdir($hmmdir));
		while (my $file = $dir->read) {
		  push(@hmmfiles, File::Spec->catfile($hmmdir, $file)) if ($file =~ /\.hmm$/);
		}
		$dir->close;
		return sort @hmmfiles;
	}
	else {
		print "single hmm\n"
		  if $verbose;
		push(@hmmfiles, $hmmfile);
		return @hmmfiles;
	}
>>>>>>> ffbb25968d9e0cdc1dac57195419dd04e634aa20
}#}}}

# Sub: translate_est
# Translate a nucleotide fasta file to protein in all six reading frames
# Expects: scalar string filename
# Returns: scalar string filename (protfile)
sub translate_est {#{{{
  my ($infile) = shift;
  (my $outfile = $infile) =~ s/(\.fa$)/_prot$1/;
  $| = 1;
  print "Translating $estfile in all six reading frames... " unless $quiet;
  if (-e $outfile) {
    print "$outfile exists, using this one.\n" unless $quiet;
    return $outfile;
  }
  my $translateline = qq('$translateprog' $infile > $outfile);
  die "Fatal: Could not translate $infile: $!\n"
    if system($translateline);

  print "done.\n";
  $| = 0;
  return $outfile;
}#}}}

# Sub: gethmmscores
# Parse the hmmsearch result, populate the data structure
# Expects: scalar reference to hmmresult file
# Returns: reference to result object
sub gethmmscores {#{{{
  my $hmmresultref = shift;
}#}}}

# sub: clean_up_old_output_files
# input: reference to list of relevant contigs
sub backup_old_output_files {#{{{
  my ($outfile) = shift @_;
  if (-e $outfile) {
    rename($outfile, File::Spec->catfile($outfile . '.' . $backup_ext)) or die "wah: could not rename $outfile during backup: $!\n";
    print "backed up old file $outfile to ", $outfile.$backup_ext, "\n";
  }
}#}}}

# Sub: helpmessage
# Prints the help message
sub helpmessage {#{{{
  my $helpmessage = <<EOF;
estfile ESTFILE
E ESTFILE 
H HMMFILE
hmmdir HMMDIR
EOF
  print $helpmessage;
}#}}}

# Sub: createdir
# Creates a directory with parent directories if it does not already exist
# Expects: scalar string dirname
# Returns: True if successful
sub createdir {#{{{
  my $dir = shift;
  if (-e $dir and not -d $dir) {
    die "Fatal: $dir exists, but is not a directory! This will most likely lead to problems later.\n";
  }
  else {
    make_path $dir, { verbose => 0 } or die "Fatal: Could not create $dir: $!\n";
    return 1;
  }
  return 1;
}#}}}

# Sub: preparedb
# Generate a clean database, deleting all existing tables and starting from scratch
# Returns: True on success
sub preparedb {#{{{
<<<<<<< HEAD
  my $query_create_ests = "CREATE TABLE $mysql_table_ests ( 
    $mysql_col_id   INT NOT NULL PRIMARY KEY AUTO_INCREMENT,
    $mysql_col_spec VARCHAR(255) NOT NULL,
    $mysql_col_date INT(11) UNSIGNED,
    $mysql_col_hdr  VARCHAR(255) NOT NULL, INDEX ($mysql_col_hdr),
    $mysql_col_seq  MEDIUMBLOB DEFAULT NULL)";  # BLOB data is stored independently of the row data and does not fall into the 65535 B limit.

  my $query_create_hmmsearch = "CREATE TABLE $mysql_table_hmmsearch (
    $mysql_col_id        INT NOT NULL PRIMARY KEY AUTO_INCREMENT,
    $mysql_col_spec      VARCHAR(255) NOT NULL,
    $mysql_col_hmm       VARCHAR(255) NOT NULL,
    $mysql_col_hmmhit    VARCHAR(255) NOT NULL, INDEX ($mysql_col_hmmhit),
    $mysql_col_score     FLOAT NOT NULL,
    $mysql_col_eval      FLOAT NOT NULL)";
  
  my $query_create_blast = "CREATE TABLE $mysql_table_blast (
    $mysql_col_id     INT NOT NULL PRIMARY KEY AUTO_INCREMENT,
    $mysql_col_spec   VARCHAR(255) NOT NULL,
    $mysql_col_query  VARCHAR(255) NOT NULL, INDEX ($mysql_col_query),
    $mysql_col_target VARCHAR(255) NOT NULL, INDEX ($mysql_col_target),
    $mysql_col_score  FLOAT NOT NULL,
    $mysql_col_eval   FLOAT NOT NULL)";
  
  my $query_create_core_orthologs = "CREATE TABLE $mysql_table_core_orthologs (
    $mysql_col_id     INT NOT NULL PRIMARY KEY AUTO_INCREMENT,
    $mysql_col_spec   VARCHAR(255) NOT NULL,
    $mysql_col_taxon  VARCHAR(25) NOT NULL,
    $mysql_col_hmm    VARCHAR(255) NOT NULL,
    $mysql_col_hdr    VARCHAR(255) NOT NULL, INDEX ($mysql_col_hdr),
    $mysql_col_seq    MEDIUMBLOB DEFAULT NULL)";  # BLOB data is stored independently of the row data and does not fall into the 65535 B limit.

  # open connection
  my $dbh = DBI->connect($mysql_dbi, $mysql_dbuser, $mysql_dbpwd);

  # drop all tables
  foreach ($mysql_table_ests, $mysql_table_hmmsearch, $mysql_table_blast, $mysql_table_core_orthologs) {
    my $query_drop = "DROP TABLE IF EXISTS $_";
    print "$query_drop\n" if $verbose;
    my $sql = $dbh->prepare($query_drop);
    $sql->execute()
      or die "Could not execute SQL query: $!\n";
  }

  # create all tables
  foreach my $query ($query_create_ests, $query_create_hmmsearch, $query_create_blast, $query_create_core_orthologs) {
    printf "$query\n" if $verbose;
    my $sql = $dbh->prepare($query);
    $sql->execute()
      or die "Could not execute SQL query: $!\n";
  }

  # disconnect
  $dbh->disconnect();

  return 1;
=======
	my $query_create_ests = "CREATE TABLE $mysql_table_ests ( 
		$mysql_col_id   INT NOT NULL PRIMARY KEY AUTO_INCREMENT,
		$mysql_col_spec VARCHAR(255) NOT NULL,
		$mysql_col_date INT(11) UNSIGNED,
		$mysql_col_hdr  VARCHAR(255) NOT NULL, INDEX ($mysql_col_hdr),
		$mysql_col_seq  MEDIUMBLOB DEFAULT NULL)";  # BLOB data is stored independently of the row data and does not fall into the 65535 B limit.

	my $query_create_hmmsearch = "CREATE TABLE $mysql_table_hmmsearch (
		$mysql_col_id        INT NOT NULL PRIMARY KEY AUTO_INCREMENT,
		$mysql_col_spec      VARCHAR(255) NOT NULL,
		$mysql_col_hmm       VARCHAR(255) NOT NULL,
		$mysql_col_hmmhit    VARCHAR(255) NOT NULL, INDEX ($mysql_col_hmmhit),
		$mysql_col_score     FLOAT NOT NULL,
		$mysql_col_eval      FLOAT NOT NULL)";
	
	my $query_create_blast = "CREATE TABLE $mysql_table_blast (
		$mysql_col_id     INT NOT NULL PRIMARY KEY AUTO_INCREMENT,
		$mysql_col_spec   VARCHAR(255) NOT NULL,
		$mysql_col_query  VARCHAR(255) NOT NULL, INDEX ($mysql_col_query),
		$mysql_col_target VARCHAR(255) NOT NULL, INDEX ($mysql_col_target),
		$mysql_col_score  FLOAT NOT NULL,
		$mysql_col_eval   FLOAT NOT NULL)";
	
	my $query_create_core_orthologs = "CREATE TABLE $mysql_table_core_orthologs (
		$mysql_col_id     INT NOT NULL PRIMARY KEY AUTO_INCREMENT,
		$mysql_col_spec   VARCHAR(255) NOT NULL,
		$mysql_col_taxon  VARCHAR(25) NOT NULL,
		$mysql_col_hmm    VARCHAR(255) NOT NULL,
		$mysql_col_hdr    VARCHAR(255) NOT NULL, INDEX ($mysql_col_hdr),
		$mysql_col_seq    MEDIUMBLOB DEFAULT NULL)";  # BLOB data is stored independently of the row data and does not fall into the 65535 B limit.

	# open connection
	my $dbh = DBI->connect($mysql_dbi, $mysql_dbuser, $mysql_dbpwd);

	# drop all tables
	foreach ($mysql_table_ests, $mysql_table_hmmsearch, $mysql_table_blast, $mysql_table_core_orthologs) {
		my $query_drop = "DROP TABLE IF EXISTS $_";
		print "$query_drop\n" if $verbose;
		my $sql = $dbh->prepare($query_drop);
		$sql->execute()
		  or die "Could not execute SQL query: $!\n";
	}

	# create all tables
	foreach my $query ($query_create_ests, $query_create_hmmsearch, $query_create_blast, $query_create_core_orthologs) {
		printf "$query\n" if $verbose;
		my $sql = $dbh->prepare($query);
		$sql->execute()
		  or die "Could not execute SQL query: $!\n";
	}

	# disconnect
	$dbh->disconnect();

	return 1;
>>>>>>> ffbb25968d9e0cdc1dac57195419dd04e634aa20
} #}}}

# Sub: clear_db
# clears the database of previous results from the same species 
sub clear_db {#{{{
<<<<<<< HEAD
  # clear previous results from the same species
  my $query_clear_ests           = "DELETE FROM $mysql_table_ests 
                                    WHERE $mysql_col_spec = '$species_name'";
  my $query_clear_hmmsearch      = "DELETE FROM $mysql_table_hmmsearch 
                                    WHERE $mysql_col_spec = '$species_name'";
  my $query_clear_blast          = "DELETE FROM $mysql_table_blast 
                                    WHERE $mysql_col_spec = '$species_name'";
  my $query_clear_core_orthologs = "DELETE FROM $mysql_table_core_orthologs 
                                    WHERE $mysql_col_spec = '$species_name'";

  # open connection
  my $dbh = DBI->connect($mysql_dbi, $mysql_dbuser, $mysql_dbpwd);

  my $sql = $dbh->prepare($query_clear_ests);
  $sql->execute() or die "$!\n";
  $sql = $dbh->prepare($query_clear_hmmsearch);
  $sql->execute() or die "$!\n";
  $sql = $dbh->prepare($query_clear_blast);
  $sql->execute() or die "$!\n";
  $sql = $dbh->prepare($query_clear_core_orthologs);
  $sql->execute() or die "$!\n";

  # disconnect ASAP
  $dbh->disconnect;
  return 1;
=======
	# clear previous results from the same species
	my $query_clear_ests           = "DELETE FROM $mysql_table_ests 
		                                WHERE $mysql_col_spec = '$species_name'";
	my $query_clear_hmmsearch      = "DELETE FROM $mysql_table_hmmsearch 
		                                WHERE $mysql_col_spec = '$species_name'";
	my $query_clear_blast          = "DELETE FROM $mysql_table_blast 
		                                WHERE $mysql_col_spec = '$species_name'";
	my $query_clear_core_orthologs = "DELETE FROM $mysql_table_core_orthologs 
		                                WHERE $mysql_col_spec = '$species_name'";

	# open connection
	my $dbh = DBI->connect($mysql_dbi, $mysql_dbuser, $mysql_dbpwd);

	my $sql = $dbh->prepare($query_clear_ests);
	$sql->execute() or die "$!\n";
	$sql = $dbh->prepare($query_clear_hmmsearch);
	$sql->execute() or die "$!\n";
	$sql = $dbh->prepare($query_clear_blast);
	$sql->execute() or die "$!\n";
	$sql = $dbh->prepare($query_clear_core_orthologs);
	$sql->execute() or die "$!\n";

	# disconnect ASAP
	$dbh->disconnect;
	return 1;
>>>>>>> ffbb25968d9e0cdc1dac57195419dd04e634aa20
}#}}}

# Documentation#{{{
=head1 NAME

Forage

=head1 DESCRIPTION

F<F>ind F<O>rthologs using F<R>eciprocity F<A>mong F<G>enes and F<E>STs

=head1 SYNOPSIS

forage.pl [OPTIONS]

=head1 OPTIONS

=head2 -c CONFIGFILE

Use CONFIGFILE instead of forage.conf. Any options on the command line override options set in the config file.

=head1 CONFIG FILE

The configuration file allows to set all options available on the command line so that the user is spared having to use a very long command every time.

The config file has to be in ini-style: 

	# this is a comment
	translateprog  = fastatranslate
	hmmdir         = /home/malty/thesis/forage/hmms
	estfile        = /home/malty/data/cleaned/Andrena_vaga.fa
	mysql_dbname   = forage
	mysql_dbserver = localhost
	mysql_dbuser   = root
	mysql_dbpwd    = root
	mysql_table    = ests

etc. Empty lines and comments are ignored, keys and values have to be separated by an equal sign. 

=head2 AVAILABLE OPTIONS:

=head2 estfile

The fasta file containing the EST sequence database.

=head2 eval_threshold

e-Value threshold for the HMM search. Must be a number in scientific notation like 1e-05 or something.

=head2 hmmdir

The directory containing the HMM files that are used to search the EST database. See also F<hmmfile>.

=head2 hmmfile

The HMM file to use if you want to run Forage on a single HMM.

=head2 hmmsearch_output_dir

The directory name where the hmmsearch output files will be placed.

=head2 hmmsearchprog

The hmmsearch program. On most systems, this is 'hmmsearch', which is the default. Change this if your binary has a different name.

=head2 mysql_dbname

MySQL: Database name. Ask your administrator if you don't know.

=head2 mysql_dbpassword

MySQL: Database password. Ask your administrator if you don't know.

=head2 mysql_dbserver

MySQL: Database server. Normally 'localhost'. Ask your administrator if you don't know.

=head2 mysql_dbuser

MySQL: Database user name. Ask your administrator if you don't know.

=head2 mysql_table_ests

MySQL: Table name for the EST sequences. Ask your administrator if you don't know.

=head2 mysql_table_hmmsearch

MySQL: Table name for the HMMsearch results. Ask your administrator if you don't know.

=head2 mysql_table_blast

MySQL: Table name for the BLAST results. Ask your administrator if you don't know.

=head2 mysql_table_core_orthologs

MySQL: Table name for the core ortholog sequences. Ask your administrator if you don't know.

=head2 output_directory

Output directory to use by Forage. It will be created if it does not exist.

=head2 quiet

Quiet operation: Forage will retain most of the status messages.

=head2 species_name

Name of the analyzed species. Will occur in the database as well as in all output files in appropriate places.

=head2 score_threshold

Score threshold for the HMM search. Must be a number.

=head2 translate_program

The tool to use for translation of the EST sequence file. The program must accept a fasta file name as input and provide its output on STDOUT, which is then being redirected into a file that Forage uses.

=head2 verbose

Be verbose (output more information about what is going on). 

=head1 AUTHOR

Written by Malte Petersen <mptrsen@uni-bonn.de>

=head1 COPYRIGHT

Copyright (C) 2012 Malte Petersen 

This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

=cut
#}}}
