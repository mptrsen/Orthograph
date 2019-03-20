ORTHOGRAPH: Orthology prediction using a Graph-based, Reciprocal Approach with Profile Hidden Markov models
===========================================================================================================

CITATION
========

Orthograph was published in [BMC Bioinformatics
(2017)](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-017-1529-8).
Please cite this paper in your academic publications. Thank you.

DOCUMENTATION
=============

A manual page containing a description of all options is in the `doc`
directory. It can be read using `man -l doc/orthograph.man` or directly with
`man orthograph` if installed properly (ask your system administrator).

This is a quickstart guide to help you start off. 

SYSTEM REQUIREMENTS
===================

Orthograph requires the following software packages to be installed on the
system. The whole package was tested and works with these versions, and newer
versions should be no problem. Downwards compatibility (i.e., older versions)
has not been tested. Should be fine, too, but be wary with MAFFT, its developers
change or drop features frequently across versions. Also, there has been
a report about an older SQLite version not working, so make sure you
have a reasonably recent one.

You need a database backend, either SQLite (the default) or MySQL. If you don't
know which one to pick, use SQLite. It's easier to install and use, but doesn't
have the capabilities for a centralized server-client setup. In grid computing
environments, this is an advantage since no network connections are required,
but the database is a flat file on the hard drive.

Package        | Version   | Download from
-------------- | --------- | -------------------------------------------------------------
Perl           | 5.14      | <http://www.perl.org>
SQLite         | 3.8.2     | <http://sqlite.org/download.html>
MySQL          | 5.6.17    | <http://dev.mysql.com/downloads/mysql/>
MAFFT          | 7.273     | <http://mafft.cbrc.jp/alignment/software/>
HMMer          | 3.1b1     | <http://hmmer.janelia.org/software/>
NCBI BLAST+    | 2.2.28+   | <ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/>
Exonerate      | 2.2.0     | <https://www.ebi.ac.uk/about/vertebrate-genomics/software/exonerate>

Make sure you get the correct versions for your operating system (most packages
are provided as 32 bit and 64 bit versions).

Don't worry about your Perl version if you run a reasonably recent Linux or
UNIX system, such as Mac OS X, but see [here][1] and [there][2] if you are
trying to install the required packages on a Windows system under Cygwin (or
natively, for that matter). You might run into problems getting the MySQL
driver, HMMer3 and Exonerate to work. Really, don't use Windows. 

Make sure you have the database interface module DBI for your Perl version as
well as the proper database driver for your database engine of choice
(DBD::mysql or DBD::SQLite)

If you get an error like "Can't locate Some/Module.pm in @INC", you need to
update your Perl installation; see REQUIRED PERL MODULES for more information.


TEST DATA
=========

To verify that your dependencies are correctly installed and work, you can test
Orthograph right away using the example config file orthograph.conf and the test
dataset in the test_data directory. It includes part of a transcriptome assembly
from Misof et al. (2014) (<http://www.ncbi.nlm.nih.gov/Traces/wgs/?val=GAUJ02>).
Test your dependencies using:

	$ ./orthograph-analyzer --test-deps

Review the path settings in the config file and adjust them accordingly if you
installed some dependencies locally.

Then run Orthograph on the supplied test file (it is already specified as the
input file in the example config):

	$ ./orthograph-analyzer

If all goes smooth, this should run without problems.


DATA REQUIREMENTS
=================

Orthograph is an orthology prediction pipeline that attempts to map transcript
sequences to known orthologous groups. It is not meant to infer orthology de
novo among genomes or transcriptomes, but relies on the orthology information in
pre-defined sets. Because of this requirement, it is essential that your input
data conforms to some data standards as listed below.


a) The clusters of orthologous groups (COGs)
--------------------------------------------

The information about which genes belong to a COG consists of a tab-delimited
file that contains at least one line per COG. If you get your COGs from
[OrthoDB 7][3] (recommended), this is the case. The table must have three
columns:

	COG_ID	sequence_ID	taxon_name

In an OrthoDB 7 file, each line has a number of tab-delimited fields like this:

	EOG7M10DZ	AECH19093	AECH19093-PA	Acromyrmex echinatior	Panamanian leafcutter ant	AECHI	1282	IPR006121,IPR008250,IPR005834,IPR023214 

The 'EOG7M10DZ' field is the COG ID. The line must also have an unambiguous
sequence ID. For OrthoDB version 7, it is always in the third field, i.e.,
'AECH19093-PA'. The taxon name is in the fourth field and the sixth field
contains the taxon shorthand. An OrthoDB 7 table can be re-formatted using
`cut`, picking out the columns 1, 3, and 4:

	$ cut -f1,3,4 ORTHODBFILE > ORTHOGRAPH_INPUT_FILE

You may need to filter the OrthoDB 7 file to contain only those taxa you want.
Probably the easiest way is to do this using `grep`:

- create a file containing the taxon shorthands you want. The shorthands must
	be unambiguously identical to the ones found in the OrthoDB file, and must
	be in a single line each, with no extra whitespace anywhere in the file (no
	trailing empty lines!).

- use `grep -f LISTFILE ORTHODBFILE > OUTPUTFILE` to filter the ORTHODBFILE
	according to the shorthand list in LISTFILE. The output will be written to
	OUTPUTFILE. 

- double-check the resulting file. Is it empty? Does it contain the right
	number of lines (number of taxa * number of COGs)?

If you want to create these files yourself, make sure that the COG ID is in the
first tab-separated field, the sequence ID in the second field, and the taxon
shorthand in the third field. 

For newer versions of OrthoDB, the table columns may be different. In that
case, the `cut` and filtering command may need to be adjusted to select
different columns.  Also make sure the table does not have any column headers.


b) The official gene set (OGS)
------------------------------

You need the OGS for each of the taxa that are present in the ortholog set,
i.e., contributed sequences to the COGs. The sequences for a single taxon must
be in Fasta format, and all (!) sequences of this taxon that occur in the COGs
must be in the OGS database as well. Their sequence identifier must be
*unambiguously identical* to the sequence ID in the file from OrthoDB. For
example, this would be the related sequence for Acromyrmex echinatior from the
OGS:

	>AECH19093-PA

Note that these IDs must not contain any commas (',').

The OGS must be present as amino acid sequences, and may be present as
nucleotide sequences. If the latter is missing, no nucleotide sequences can be
output for the reference taxa. The sequence identifiers of corresponding
nucleotide and amino acid sequences must be identical.

The helper script make-ogs-corresponding.pl is provided that uses Exonerate to
generate 100% corresponding amino acid and nucleotide sequences for the OGS. It
only requires that the identifiers are unambiguously identical. You can call it
like this:

	perl make-ogs-corresponding.pl [OPTIONS] PEPTIDEFILE CDSFILE

Provide the OGS on amino acid level as PEPTIDEFILE and on nucleotide level as
CDSFILE. A list of options is printed when you call the script without arguments
or with the -h option.


QUICKSTART
==========

For a fast start into analysis using Orthograph, follow these steps. They
assume that all required software packages are installed and work. 


## 1. Some preliminary remarks

Throughout this manual, the following marks a shell command:

	$ some_command

This means you should enter some_command at the shell prompt. The `$` sigil is
not part of the command.

All Orthograph tools adhere to the `-c` option that allows you to specify the
path to a config file, e.g.:

	$ orthograph-analyzer -c different_species.conf

This is useful if you use different SQLite or MySQL databases, for example, or
if you want to generate config files for dozens of species.


## 2. Create a config file. 

Make a copy of the supplied config file orthograph.conf.example and rename this
copy to orthograph.conf. Edit it. The file is commented to help you start off.
Most options are self-explanatory. 

For now, it is sufficient to supply the database settings (for SQLite: path to
the database file; for MySQL: database, username, password). If you don't know
how to create a MySQL user and a database, see the appendix (the very bottom of
this file).

Note: The config file is read by all three Orthograph tools. Especially the
mandatory part of the file is important as it contains settings for the database
interaction. Without these settings none of the programs will work.


## 3. Create the tables.

Make `orthograph_manager` executable and use it:

	$ chmod +x orthograph-manager
	$ ./orthograph-manager --create

The `-create` switch tells the script to setup the database structure for
Orthograph.  It will create a number of tables, all of which start with the
prefix `orthograph_` (or whatever prefix you chose in the config file). 

At this point, if you use the SQLite backend, the SQLite database file is
created at the location you specified in the config file.


## 4. Load OGS sequences into the database :

	$ ./orthograph-manager --load-ogs-peptide FASTAFILE

You need to load data for every taxon separately (this may change in the
future), and you should load data for as many taxa in the orthologs file as
possible (preferably all). You will be asked for information about the file you
are loading.  These sequences are used for the reciprocal BLAST databases. If
you don't load any peptide data, Orthograph will not run properly.

Some of the protein IDs (in the first header field of the OGS file) must
correspond to the third field in the OrthoDB file. This is important, as these
are regarded as part of ortholog groups. Repeat this step for as many taxa in
your ortholog set as you have peptide sequences for.

If you made a mistake during input of the information, you may cancel the
process at any time using `Ctrl+C`. To check what OGS are present in the
database, call `orthograph-manager -lo`.

There are options to specify OGS taxon and version so you can easily automate
the process:

	--ogs-version VERSION
	--ogs-taxon-name TAXON

Both version and name can be an arbitrary string (but enclose it in quotes if it
contains spaces!). 


## 5. Optional: Load nucleotide sequence data into the database:

	$ ./orthograph-manager --load-ogs-nucleotide FASTAFILE

The nucleotide sequences are expected to be the coding sequences for the
peptide sequences you loaded earlier. The headers have to correspond EXACTLY,
or they cannot be correlated and nucleotide output is not possible. For
example, these two headers are not identical, even though they look similar:

	>AAEL1307285-PA Q0C738 Bystin IPR007955
	>AAEL1307285-RA Q0C738 Bystin IPR007955

Make sure the nucleotide sequences can be correlated to the proteome sequences
by giving identical names to corresponding sequences.

The same options for specifying OGS taxon and version number (see step 5) are
available for nucleotide data, as well.

I know. Nucleotide OGS loading is slow. I'm working on improving that
performance. For now, be happy that you only need to do this once.


## 6. Upload your ortholog set into the database:

	$ ./orthograph-manager FILE

You will be asked for information about the ortholog set you are loading. If you
made a mistake during input of the information, you may cancel the process at
any time using `Ctrl+C`. I'm sorry if you mistyped the last taxon name and have
to start over :)

The file must have a three-column tab-delimited format as described
in the DATA REQUIREMENTS section. If you downloaded your ortholog set from
OrthoDB 7, you can use the instructions above to format the table accordingly.

Dependent on your OrthoDB query, you may need to filter your file so that it
contains only the taxa you want in the ortholog set. 

Make a note of your ortholog set name. You will need it later.

To check what sets have been uploaded so far, call `orthograph-manager -ls`. 


## 7. Create the required database structure for the Orthograph search program.

Make `orthograph-analyzer` executable if it isn't already and run it with the
-prepare option:

	$ chmod +x orthograph-analyzer
	$ ./orthograph-analyzer --prepare

At this point, if you use the SQLite backend, a species-specific SQLite database
file is created in the output directory you specified in the config file.


## 8. Complete your config file if you didn't fill out everything before.

Change the setting of `ortholog-set` to the name of your custom ortholog set.
Also supply the rest of the mandatory settings (`species-name`, `input-file`,
and `output-directory` (recommended)) and consider changing some optional ones
that will affect the searches.


## 9. Start Orthograph! 

	$ ./orthograph-analyzer


## 10. Get some coffee, sit back and watch (or do something else).

Orthograph generates the profile hidden Markov models (pHMMs) for your ortholog
set, if they don't exist. This may take a long time depending on the size of
your set, but rest assured, once the pHMMs exist, subsequent analyses will be
faster by that amount of time. 

After that, your input file is translated into all six possible reading frames.
The translated file is placed in your output directory. This library is searched
using all pHMMs. Candidate orthologs are verified with a reciprocal search
against the BLAST database of all proteomes in your ortholog set. The results
are cached in the database.

The HMMER and BLAST output tables are placed in the 'hmmsearch' and 'blast'
subdirectories in the output directory. The 'aa' and 'nt' directories are
created but left empty. A log file called 'orthograph-analyzer-DATE.log' is
created in the 'log' directory (or wherever you specified).


## 11. Start the reporter

Make `orthograph-reporter` executable if it isn't already and run it:

	$ chmod +x orthograph-reporter
	$ ./orthograph-reporter

`orthograph-reporter` fetches the search results from the database and assigns
ortholog relations by clustering best reciprocal hits. After it is finished,
the output directory contains five subdirectories:

a) hmmsearch. This contains all the HMMsearch output tables. 

b) blast. This contains all the BLASTP output tables.

c) aa. This contains the actual result of your analysis. The hit sequences are
output along with the core-ortholog sequences from your ortholog set, grouped
by ortholog groups. See the appendix for information on Fasta header structure.

d) nt. This contains the actual results on nucleotide level.

e) log. This contains log files, such as the entire standard output for both
`orthograph-analyzer` and `orthograph-reporter`:


The main output directory contains two report tables after running
orthograph-reporter:

- best-reciprocal-hits.txt: a tabular listing of all sequences (sections) that
  fulfilled the BRH criterion, _irrespective_ of whether they overlap with
  other sequences. It is structured as follows (tab-separated):
  - the COG ID
  - transcript ID
  - start on the transcript
  - end on the transcript
  - HMM score
  - HMM e-value

- non-overlapping-best-reciprocal-hits.txt: a tabular listing of all sequences
  (sections) that fulfilled the BRH criterion but do _not_ overlap with other
  hits on transcript level. They may overlap on HMM level, though (paralogous
  copies perhaps). Its structure is identical to the one of
  best-reciprocal-hits.txt (see above).

- summary.txt: a tabular listing of all transcripts that were eventually mapped
  to each COG. These do not overlap, fulfill all other criteria (length, number
  of mismatches, etc.) and there are aa and nt output files for them.


## 12. Optional: Summarize output

For further analyses, it may be required to summarize the sequences from
multiple Orthograph output directories. This is accomplished using
`summarize_orthograph_results.pl` by Oliver Niehuis, which can be run like
this:

	$ perl summarize_orthograph_results.pl -i INPUT_DIRECTORY -o OUTPUT_DIRECTORY [OPTIONS]

The program generates an amino acid resp. nucleotide sequence file in the
OUTPUT_DIRECTORY for each COG containing the sequences from all taxa, including
reference taxa. The INPUT_DIRECTORY must be a directory containing Orthograph
output directories. For further information, read the usage instructions
that are printed when calling the program with the `-h` flag.


## 13. Optional: Convert to HaMStRad format

For downstream analyses that depend on HaMStRad output, it is necessary to
convert the Fasta header format accordingly. There is a converter script to do
that; it is run like this:

	$ perl orthograph2hamstrad.pl INPUT_DIRECTORY OUTPUT_DIRECTORY

The program converts all files in INPUT_DIRECTORY/aa and INPUT_DIRECTORY/nt and
places them in OUTPUT_DIRECTORY/aa and OUTPUT_DIRECTORY/nt, respectively.
OUTPUT_DIRECTORY must exist, the aa and nt subdirectories will be created.



APPENDIX
========

In the output Fasta files, the headers are structured like this:

	>COG_ID|Taxon_name|Sequence_ID|Start-End|Reading_frame|Reference_taxon

For example:

	>EOG7VQW3C|Nasonia_vitripennis|NV15365-RA|4-604|[translate(1)]|Tribolium castaneum

This would mean: For the ortholog group EOG7VQW3C and the analyzed taxon with
the name Nasonia_vitripennis, the transcript with the ID NV15365-RA, translated
in reading frame 1 from coordinates 4 to 604 was assigned, and the reference
taxon (the one with the best reciprocal hit) was Tribolium castaneum.

If multiple transcripts were concatenated, the last four fields are repeated
for each transcript, separated by '&&':

	>EOG7096NB|Nasonia_vitripennis|NV50397-RA|128-553|[translate(1)]|Nasonia vitripennis&&NV50397-RA|578-1309|[translate(1)]|Nasonia vitripennis

The separators can be specified via the options "header-separator" and
"concatenation-header-separator", respectively.


MySQL performance
-----------------

Orthograph makes extensive use of MySQL query caching. For optimal performance,
make sure that query caching is enabled in your MySQL instance. Also, if you
plan to run multiple orthograph instances on the same database, you should
definitely use table files instead of database files (`innodb_file_per_table`)
and increase the lock wait timeout (`innodb_lock_wait_timeout`) to between 200
and 600 to make sure that long queries (which may occur) do not block other
orthograph processes.


Required Perl modules
---------------------
Orthograph uses only modules that are present in any standard Perl distribution. 
Additional modules such as IO::Tee and File::Which are shipped with Orthograph.
If you get an error like "Can't locate Some/Thing.pm in @INC" it means that
Perl can't find one of the required modules. Update your Perl installation if
this is the case and make sure it contains these modules:

-  autodie
-  strict
-  warnings
-  Archive::Tar
-  Benchmark 
-  Carp
-  Config
-  DBD::mysql
-  DBD::SQLite
-  DBI
-  Data::Dumper
-  Digest::SHA
-  File::Basename
-  File::Path
-  File::Spec
-  File::Temp
-  FindBin
-  Getopt::Long
-  IO::Dir
-  IO::File
-  List::Util 
-  Time::HiRes

Part of my Diploma thesis at the ZFMK/zmb, Bonn, Germany 
(c) 2011-2012 Malte Petersen <mptrsen@uni-bonn.de>

[1]: <http://cpansearch.perl.org/src/JWIED/DBD-mysql-2.1028/INSTALL.html#special%20systems>
[2]: <http://forums.mysql.com/read.php?51,389833,389833>
[3]: <http://orthodb.org/>


## Set up the database. 

Skip this step if you already have a MySQL account and a database (and know
what your credentials are).

	mysql> SOME STATEMENT;

This means you should enter SOME STATEMENT at the MySQL prompt. The `mysql>` is
not part of the statement, but the semicolon is.

There must exist a database for your usage and a database user must have the
usual permissions on it. If you are the administrator of your local machine and
installed MySQL yourself, you have those privileges. However, you should never
do ordinary work as root (the administrator), so let's create a working user.
Login to MySQL as root:

	$ mysql -u root -p
	Enter password:

Enter the MySQL administrator password. You will be greeted by MySQL and then
given a new prompt. Issue the following statement:

	mysql> CREATE USER 'USERNAME'@'localhost' IDENTIFIED BY 'PASSWORD';

Substitute USERNAME and PASSWORD with your username and password of choice. The
single quotes are important. Next, you will create a database and grant the
working user the necessary permissions on it:

	mysql> CREATE DATABASE orthograph;
	mysql> GRANT ALL ON orthograph.* TO 'USERNAME'@'localhost';

You may use a different database name, of course. Make a note of your new user
name, the password and the database name. You will need it later.


Common problems
---------------

"No sequences found. Something went wrong. Check your input file. Exiting." 

Really, check your input file. Is it empty? Is it a valid fasta file? If it
looks ok, try deleting the translated version so Orthograph is forced to
recreate it. 


Options
-------

Orthograph accepts a number of options, both on the command line and in the
config file. In fact, you can customize just about anything. All options are
documented in the Orthograph manual.
