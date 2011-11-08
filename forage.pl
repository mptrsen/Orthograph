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
use Getopt::Long;
use File::Path qw(mkpath);	# mkdir with parent dirs
use Tie::File;
use lib './';
use Forage::Item;
use Genetic::Codes;

#--------------------------------------------------
# # Variable initialisation
#-------------------------------------------------- 
my $version = 0.00001;#{{{
my $verbose = 0;
my $outdir;
my $estfile = '';
my $hmmfile = '';
my $hmmdir = '';
my @hmmfiles;
my $hmmsearchprog = 'hmmsearch';
my @hmmsearchcmd;
my $hmmfullout = 0;
my $hmmoutopt;
my $hmmresultfileref;
my $eval_threshold;
my @eval_option = ();
my $score_threshold;
my @score_option = ();
my @seqobjs;
my $translatecmd = 'fastatranslate';
my $i;
my $header = <<EOF;
Forage: Find Orthologs using Reciprocity Among Genes and ESTs
Copyright 2011 Malte Petersen <mptrsen\@uni-bonn.de>
Version $version

EOF

GetOptions(	'v'					=> \$verbose,
						'estfile=s' => \$estfile,
						'E=s'				=> \$estfile,
						'eval=s'		=> \$eval_threshold,
						'score=s'		=> \$score_threshold,
						'H=s'				=> \$hmmfile,
						'hmmdir=s'	=> \$hmmdir,
						'hmmsearchprog=s'	=> \$hmmsearchprog,
						'hmmfullout'	=> \$hmmfullout,
);#}}}

#--------------------------------------------------
# # Input error checking, reporting etc
#-------------------------------------------------- 
&intro;

#--------------------------------------------------
# # create list of HMM files
#-------------------------------------------------- 
@hmmfiles = &hmmlist;

print "Using HMM dir $hmmdir with ", scalar @hmmfiles, " HMMs\n" 
	if $hmmdir;
print "Using HMM file $hmmfile.\n" 
	if $hmmfile;

print "e-Value cutoff: $eval_threshold.\n" 
	if $eval_threshold;

print "Score cutoff: $score_threshold.\n"
	if $score_threshold;

#--------------------------------------------------
# # translate the ESTs to protein
#-------------------------------------------------- 
print "Translating EST file to protein... ";
my $protfile = &fastatranslate_est($estfile);

#--------------------------------------------------
# # hmmsearch the protfile using all HMMs
#-------------------------------------------------- 
print "Hmmsearching the protein file using all HMMs... ";
foreach my $hmmfile (@hmmfiles) {
	++$i;
	my $hmmresultref = &hmmsearch($hmmfile, $protfile, $outdir);
	#print $hmmresultref;
	#&gethmmscores($hmmresultfileref);
}
print "$i HMM files processed.\n";


# parse the protfile
#--------------------------------------------------
# &parse_protfile($protfile, @hmmfiles);
#-------------------------------------------------- 

###################################################
# # Functions follow
###################################################

# Sub: intro
# Checks input etc.
# Returns True if everything is OK.
sub intro {#{{{
	die "Fatal: At least two arguments required: EST file (-E) and HMM file/dir (-H or -hmmdir)!\n"
		unless ($estfile and ($hmmdir or $hmmfile));

	die "Fatal: Can't use both e-value and score thresholds!\n"
		if ($eval_threshold and $score_threshold);

	$outdir = "out_$estfile";


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

	print $header;

	print "EST file " and &doesexist($estfile);

	print "HMMsearch output dir " and &createdir("$outdir/$hmmdir");
}#}}}

# Sub: hmmlist
# Expects: nothing (?)
# Returns: array hmmfiles
sub hmmlist {#{{{
	if ($hmmdir) {
		$hmmdir =~ s/\/+$//;
		opendir(my $dir, $hmmdir) or die "Fatal: Could not open HMM dir: $!\n";
		while (my $file = readdir $dir) {
			push(@hmmfiles, "$hmmdir/$file") if ($file =~ /\.hmm$/);
		}
		closedir($dir);
		return sort @hmmfiles;
	}
	else {
		print "single hmm\n"
			if $verbose;
		push(@hmmfiles, $hmmfile);
		return @hmmfiles;
	}
}#}}}

# Sub: fastatranslate_est
# Translate a nucleotide fasta file to protein in all six reading frames
# Expects: scalar string filename
# Returns: scalar string filename (protfile)
sub fastatranslate_est {#{{{
	my ($infile) = @_;
	(my $outfile = $infile) =~ s/(\.fa$)/_prot$1/;
	if (-e $outfile) {
		print "$outfile exists, using this one.\n";
		return $outfile;
	}
	&backup_old_output_files($outfile);
	my $translateline = $translatecmd . " " . $infile . ">$outfile";
	print "Translating $infile, please wait...\t";
	die "Fatal: Could not translate $infile: $!"
		if system($translateline);
	print "done!\n";
	return $outfile;
}#}}}

# Sub: hmmsearch
# HMM search a sequence using HMM, leaving an outfile for later processing
# Expects: reference to sequence object, scalar string filename to HMM
# Returns: scalar reference to hmmoutfile
sub hmmsearch {#{{{
	my ($hmmfile, $protfile, $outdir) = @_;
	# full output if desired, table only otherwise; reflects in outfile extension
	my $hmmoutfile = $hmmfullout ? $hmmfile.'.out' : $hmmfile.'.tbl';
	# e-value and score options if desired
	if (-e "$outdir/$hmmoutfile") {
		print "hmmsearch result file $hmmoutfile already exists, skipping this HMM\n"
			if $verbose;
		return $hmmoutfile;
	}
	else {
		print join " ", (@hmmsearchcmd, "$outdir/$hmmoutfile", $hmmfile, $protfile, "\n")
			if $verbose;
		die "Fatal: hmmsearch failed on $protfile with HMM $hmmfile\n" 
			if system(@hmmsearchcmd, "$outdir/$hmmoutfile", $hmmfile, $protfile);
		return $hmmoutfile;
	}
}#}}}

# Sub: gethmmscores
# Parse the hmmsearch result, populate the data structure
# Expects: scalar reference to hmmresult file
# Returns: reference to result object
sub gethmmscores {#{{{
	my $hmmresultref = shift;
	print "Tieing hmmresult file $hmmresultref to an array\n"
		if $verbose;
	tie(my @file, 'Tie::File', $hmmresultref) 
		or die "Fatal: Could not open HMM result file: $!\n";
	untie @file 
		or die "Warning: Could not close HMM result file: $!\n";
	print "$hmmresultref successfully untied\n"
		if $verbose;
}#}}}

# sub: clean_up_old_output_files
# input: reference to list of relevant contigs
sub backup_old_output_files {#{{{
	my ($outfile) = shift @_;
	my $backup_ext = '.bak';
	if (-e $outfile) {
		rename( $outfile, $outfile.$backup_ext ) or die "wah: could not rename $outfile during backup: $!\n";
		print "backed up old file $outfile to ", $outfile.$backup_ext, "\n";
	}
}#}}}

# Sub: helpmessage
# Prints the help message
sub helpmessage {#{{{
	my $helpmessage = <<EOF;
'estfile=s' 
'E=s'				
'H=s'				
'hmmdir=s'	
EOF
	print $helpmessage;
}#}}}

# Sub: createdir
# Creates a directory with parent directories if it does not already exist
# Expects: scalar string dirname
# Returns: True if successful
sub createdir {#{{{
	my $dir = shift;
	unless (-d $dir) {
		if (-e $dir and not -d $dir) {
			die "Warning: $dir exists, but is not a directory! This will most likely lead to problems later.\n";
		}
		else {
			print "Creating $dir... ";
			File::Path->make_path($dir) or die "Fatal: Could not create $dir: $!\n";
			print "OK\n";
			return 1;
		}
	}
	print "$dir exists.\n";
	return 1;
}#}}}

# Sub: doesexist
# Checks whether an item exists on disk
# Expects: scalar string filename/dirname
# Returns: True if exists, False otherwise
sub doesexist {#{{{
	my $item = shift;
	unless (-e $item) {
		warn "$item does not exist!\n";
		return 0;
	}
	print "$item exists.\n";
	return 1;
}#}}}
