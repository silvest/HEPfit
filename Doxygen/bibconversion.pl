#/
# Copyright (C) 2014 SusyFit Collaboration
# All rights reserved.
#
# For the licensing terms see doc/COPYING.
#/

#!/usr/bin/perl -w
# -*- perl -*-
# vim: syntax=perl
eval 'exec /usr/bin/perl -w -S $0 ${1+"$@"}'
if 0;

$version = '@VERSION@';

use Getopt::Std;
use open IO => ':crlf';
use Term::ANSIColor;
use Config;
use strict;
use warnings;

my @tmpfiles = ();
my @bibfiles = ();
my @doxygen_list = ();
my $OS = $Config{osname};
my @keyword_list = ("QCD", "Model");

sub usage {
    my $program = $0;
    $program =~ s+^.*/++;
    print STDERR <<_EOF_;
    
        USAGE: perl $program <list of bibtex files>
    
_EOF_

    exit(1);
}

sub header {
    print color 'bold blue on_white';
    print STDOUT <<_EOF_;


    Copyright (C) 2014 SusyFit Collaboration
    All rights reserved.
    
    For the licensing terms see doc/COPYING.
_EOF_
    print color 'reset';
    print "\n\n";
}

&usage if ($#ARGV < 0);
&header;

unlink ($ARGV[-1]);

foreach (@ARGV) {
    if (/\.bib$/ && \$_ != \$ARGV[-1]) {
        my $bibfile = $_;
        print "\t--> $bibfile added to list\n" if (push(@bibfiles,$bibfile));
    } elsif (/\.bib$/ && (-e $ARGV[-1])){
        print colored ['Red'], "\n\n\tERROR:";
        print " Output file exists. Please choose another name or move the existing file.\n\n";
        exit(1)
    } elsif (!(/\.bib$/)) {
        print colored ['Red'], "\n\n\tERROR:";
        print " Please list only .bib files in the correct BibTex format.\n";
        &usage;
    } else {
        open BIBOUT, "+>", pop @ARGV;
        open BIBOUT_TEMP, ">", "temp.bib";
    }
}

foreach (@bibfiles) {
    open BIBIN , "<", $_;
    while (<BIBIN>){
        s/\%\n//g;
        my $z = 0;
        $_ =~ s/\"\{/\"\[/;
        $_ =~ s/\}\"/\]\"/;
        $_ =~ s/\\/\[BS\]/g;
        $_ =~ s/\$/(++$z % 2) != 0 ? "MJ\[" : "\$"/eg;
        $_ =~ s/\$/(++$z % 2) == 0 ? "\]MJ" : "\$"/eg;
        if (/MJ\[/ || /\]MJ/){
            $_ =~ s/\{/\\\{/g;
            $_ =~ s/\}/\\\}/g;
        }
        $_ =~ s/\"\[/\"\{/;
        $_ =~ s/\]\"/\}\"/;
        $_ =~ s/\n//g if (/year/ && /\=/);
        print BIBOUT_TEMP $_;
    }
    close BIBIN;
}
close BIBOUT_TEMP;
open BIBOUT_TEMP, "<temp.bib";
push(@tmpfiles,"temp.bib");

while (<BIBOUT_TEMP>){
    if (/year/){
        $_ =~ s/\s+eprint\s+/\[arXiv\:/g;
        $_ =~ s/(\",)/ /;
        $_ =~ s/\s+note/\",\n      note/g;
        $_ =~ s/\s+pages/\",\n      pages/g;
        $_ =~ s/\s+SLACcitation/\",\n      SLACcitation/g;
        $_ =~ s/\s+reportNumber/\",\n      reportNumber/g;
    }
    $_ =~ s/(arXiv.*)(\",)/$1\]$2/ if (!(/archivePrefix/));
    $_ =~ s/\:= \"/\:/g;
    foreach my $key(@keyword_list){
        $_ =~ s/\b$key\b/\%$key/g;
    }
    print BIBOUT $_;
}
    
close BIBOUT;
close BIBOUT_TEMP;
system("cp SusyFit.bib temp.bib");
system("sed -e \'s\/\\\$\/\\\]MJ/g\' temp.bib > SusyFit.bib"); # A Fistful of dollars get left over no matter what!!
            
print "\n";
print colored ['Green'], "\n\tStarting bibtex conversion...\n\n";


if (`which doxygen`){
    @doxygen_list = `which doxygen`;
}
if ($OS eq 'darwin'){
    my $doxyapp = "/Applications/Doxygen\.app/Contents/Resources/doxygen";
    push(@doxygen_list, $doxyapp) if (-e $doxyapp);
}
if (@doxygen_list == 0){
    print colored ['Red'], "\n\tERROR:";
    print " No Doxygen installation found.\n\n";
    exit(1)
}
my $doxygen_asked = 1;
my $doxygen_used = "";
if ($#doxygen_list > 0){
    my $doxygen_size = @doxygen_list;
#   @doxygen_index = (1..$doxygen_size);
    my $doxygen_index = 0;
    print colored ['Green'], "\tFound more than 1 doxygen installation: \n";
    print "\t[0] NO Doxygen run.\n";
    foreach (@doxygen_list){
        $doxygen_index += 1;
        print "\t[$doxygen_index] $_";
    }
    print colored ['Green'], "\n\n\tPlease enter the number for the doxygen executable (DEFAULT: 1).\n\t";
    $doxygen_asked = <STDIN>;
    $doxygen_asked = 1 if ($doxygen_asked eq "\n");
    while (($doxygen_asked > $doxygen_size) || ($doxygen_asked < 0)){
        print colored ['Red'], "\n\tERROR:";
        print " illegal option: $doxygen_asked";
        print "\tPlease enter a number between 0 and $doxygen_size.\n\t";
        $doxygen_asked = <STDIN>;
    }
}
if ($doxygen_asked != 0){
    chomp($doxygen_used = $doxygen_list[$doxygen_asked - 1]);
    print "\tUsing doxygen executable at $doxygen_used.\n";

    print colored ['Green'], "\n\tDo you want a fresh Doxygen build? (yes or no) (DEFAULT: no)\n";
    print "\tThis will remove the html directory in the current directory\n\t";
    chomp(my $delete_asked = <STDIN>);
    $delete_asked = "no" if ($delete_asked eq "");
    while (($delete_asked ne "yes") && ($delete_asked ne "no")){
        print colored ['Red'], "\n\tERROR:";
        print " illegal option: $delete_asked";
        print "\tPlease enter yes or no.\n\t";
        $delete_asked = <STDIN>;
    }
    if ($delete_asked eq "yes"){
        if (system(`rm -r html/*`)){
            print "\n\tDeleted the html directory.\n";
        } else {
            print "\n\tFailed to delete the html directory.\n";
            exit(1)
        }
    }
    print colored ['Green'], "\n\tRunning Doxygen...\n";
    system($doxygen_used);
} else {
    print "\tNO DOXYGEN RUN.\n\n";
}

#print colored ['Green'], "\n\tPatching index.html...\n";
#
#chomp(my $index_path = `find ./html -name "index.html"`);
#if (!(-e $index_path)){
#    print colored ['Red'], "\n\tERROR:";
#    print " index.html file not found recursively in the html directory.\n";
#    exit(1)
#} else {
#    print "\tindex.html found at: $index_path\n";
#    open INDEX, "<", $index_path;
#    open INDEXOUT, ">index.html";
#}
#
#chomp(my $citelist_path = `find ./html -name "citelist.html"`);
#if (!(-e $citelist_path)){
#    print colored ['Red'], "\n\tERROR:";
#    print " citelist.html file not found recursively in the html directory.\n";
#    exit(1)
#} else {
#print "\tcitelist.html found at: $citelist_path\n";
#}
#$citelist_path =~ s/\.\/html\///;

#while (<INDEX>){
#    #$_ =~ s/citelist.html/$citelist_path/g;
#    print INDEXOUT $_;
#}
#close INDEX;
#close INDEXOUT;
#
#system("mv index.html html/index.html");

#if ($doxygen_asked != 0){
#    print colored ['Green'], "\n\tRe-Running Doxygen to generate main page...\n";
#    system($doxygen_used);
#}

#exec("open html\/index.html") if ($OS eq 'darwin');

unlink(@tmpfiles);

exit(0);