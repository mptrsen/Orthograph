#!/usr/bin/perl
use strict;   # make me write good code
use warnings; # cry if something seems odd

# Forage: Find Orthologs using Reciprocity Among Genes and ESTs

print "FORAGE: Find Orthologs using Reciprocity Among Genes and ESTs\n";
print "(c) 2011 Malte Petersen\n";
print "version 0.00000001\n";
print "\n";
print "What HaMStR does:\n";
print "\n";
print "Step 0: hmmbuild the HMMs from the core orthologs\n";
print "Step 1: translate the EST library to protein sequences somehow\n";
print "Step 1: hmmsearch the ProtEST library using the HMM built in Step 0\n";
print "Step 2: if obtained hits, re-blast them against a reference proteome\n";
print "Step 3: if obtained hits, translate the EST to protein using genewise (?)\n";
exit;
