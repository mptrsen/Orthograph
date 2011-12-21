package Forage::Unthreaded;
use strict;
use warnings;
use File::Basename;
use File::Spec;
use Data::Dumper;
1;

sub new {
	my ($class, $hmmervars) = @_;

	
	my $self = {
		'hmmfile'				=> $hmmervars->{'hmmfile'},
		'protfile'			=> $hmmervars->{'protfile'},
		'hmmresultfile'	=> '',
		'hmmresult'			=> '',
		'hmmhits'				=> 0,
	};
	bless ($self, $class);
	return $self;
}

# Sub: hmmsearch
# HMM search a sequence using HMM, leaving an outfile for later processing
# Expects: reference to sequence object, scalar string filename to HMM
# Returns: scalar reference to hmmoutfile
sub hmmsearch {#{{{
	my ($self, $hmmervars) = @_;
	my $hmmfile		= $self->{'hmmfile'};
	my $protfile	= $self->{'protfile'};
	my $outdir		= $hmmervars->{'outdir'};
	my $hmmfullout	= $hmmervars->{'hmmfullout'};
	my $hmmsearchcmd = $hmmervars->{'hmmsearchcmd'};
	my $hmmsearchprog	= $hmmervars->{'hmmsearchprog'};
	my $verbose			= $hmmervars->{'verbose'};
	# full output if desired, table only otherwise; reflects in outfile extension
	my $hmmoutfile = $hmmfullout ? File::Spec->catfile($$outdir, $$hmmsearchprog, basename($$hmmfile).'.out') : File::Spec->catfile($$outdir, $$hmmsearchprog, basename($$hmmfile).'.tbl');
	# e-value and score options if desired
	if (-e $hmmoutfile) {
		$self->{'hmmresultfile'} = $hmmoutfile;
		print "hmmsearch result file $hmmoutfile already exists, using this one\n"
			if $$verbose;
		open(my $fh, '<', $self->{'hmmresultfile'})
			or die 'Fatal: Could not open hmmsearch result file ', $self->{'hmmresultfile'}, ': ', $!, "\n";
		$self->{'hmmresult'} = [ <$fh> ];
		close $fh;
		splice( @{$self->{'hmmresult'}}, 0, 3 );	# remove first 3 lines since they are comments anyway
		return $self;
	}
	else {
		my @hmmsearchline = (@$hmmsearchcmd, $hmmoutfile, $$hmmfile, $$protfile);
		print join " ", @hmmsearchline, "\n"
			if $$verbose;
		my $hmmresult = [ `@hmmsearchline` ];
		die "Fatal: hmmsearch failed on $$protfile with HMM $$hmmfile: $!\n" 
			unless (scalar @$hmmresult);
		# only save those results that actually match something
		unless (grep( /No hits detected/, @$hmmresult )) {
			$self->{'hmmresult'} = $hmmresult;
			splice( @{$self->{'hmmresult'}}, 0, 3 );	# remove first 3 lines since they are comments anyway
			$self->{'hmmresultfile'} = $hmmoutfile;
			return $self;
		}
		else {
			$self->{'hmmresult'} = 0;
			$self->{'hmmresultfile'} = $hmmoutfile;
			return $self;
		}
	}
}#}}}

sub hitcount {
	my ($self, $hmmvars) = @_;
	return scalar(@{$self->{'hmmresult'}}) ;	
}

sub hits {
	my ($self, $hmmvars) = @_;
	open(my $fh, '<', $self->{'hmmresultfile'})
		or die 'Fatal: Could not open hmmsearch result file ', $self->{'hmmresultfile'}, ': ', $!, "\n";
	$self->{'hmmresult'} = [ <$fh> ];
	close $fh;
	unless ($hmmvars->{'hmmfullout'}) {
		$self->{'hmmhits'} = scalar(@{$self->{'hmmresult'}}) ;	# -3 because the first 3 lines of the table are comments
	}
}
