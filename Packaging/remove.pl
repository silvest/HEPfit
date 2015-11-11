#!/usr/bin/perl
use strict;
use warnings;
no warnings qw(once);
use File::Copy;

if (@ARGV < 1) {
    print "Usage: remove.pl <file>\n";
    exit 1;
}
my $FILE = $ARGV[0];

# create a backup file
my $backup = "$FILE" . ".bak";
move($FILE, $backup);

open(INFILE,"<$backup") || die "cannot open the input file!";
open(OUTFILE,">$FILE") || die "cannot open the output file!";

# output 	
my $IsREMOVED = "NO";
while(<INFILE>) {
    if (/BEGIN: REMOVE FROM THE PACKAGE/) {
	$IsREMOVED = "YES";
    } 
    if ($IsREMOVED eq "NO") {
	print OUTFILE $_;
    }
    if (/END: REMOVE FROM THE PACKAGE/) {
	$IsREMOVED = "NO";
    } 
}

close(INFILE);
close(OUTFILE);
	
