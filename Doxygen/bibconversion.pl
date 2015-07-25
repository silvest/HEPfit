#/
# Copyright (C) 2014 HEPfit Collaboration
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
use File::Copy;

my @tmpfiles = ();
my @bibfiles = ();
my @doxygen_list = ();
my $OS = $Config{osname};
my @keyword_list = ("QCD", "Model");
my $end_index = 0;
my $bibfile_out = "";
my $doxyfile_path = "";
my $doxyfile_specified = "";

sub usage {
    my $program = $0;
    $program =~ s+^.*/++;
    print STDERR <<_EOF_;
    
        USAGE:
        perl $program <list of input bibtex files (*.bib)> [-o[-of]] [<output bibtex file (<output>.bib)>] [-dox Doxyfile]
    
        The option specifications are positional.
        The above format must be followed if the options are used.
        
        -o      Specifies the output file name. (MUST end in .bib)
        -of     Deletes the output file if it exists.
        -dox    Specifies the Doxygen configuration file from the commandline
            
        Default Output File HEPfit.bib can be specified with:
        perl $program <list of input bibtex files (*.bib)> [-of]
            
_EOF_
            
    exit(1);
}

sub header {
    print color 'bold blue on_white';
    print STDOUT <<_EOF_;
    
    
    Copyright (C) 2014 HEPfit Collaboration
    All rights reserved.
    
    For the licensing terms see doc/COPYING.
_EOF_
    print color 'reset';
    print "\n\n";
}

&usage if ($#ARGV < 0);
&header;

if ($ARGV[-1] eq "-dox"){
    print colored ['Red'], "\n\n\tERROR:";
    print " Please choose a name for the Doxygen configuration file.\n\n";
    &usage;
    exit(1)
}

if ($ARGV[-2] eq "-dox"){
    $doxyfile_path = pop @ARGV;
    $doxyfile_specified = pop @ARGV;
}

if (($ARGV[-1] eq "-of") || ($ARGV[-2] eq "-of")){
    unlink ("HEPfit.bib");
    unlink ($ARGV[-1]);
}

if ($ARGV[-1] eq "-o"){
    print colored ['Red'], "\n\n\tERROR:";
    print " Please choose a name for the output file.\n\n";
    &usage;
    exit(1)
}

if ($ARGV[-2] ne "-o" && $ARGV[-2] ne "-of"){
    print colored ['Yellow'], "\n\n\tWARNING:";
    print " Output file not specified or the optional -o option not set.\n";
    if (-e "HEPfit.bib"){
        print colored ['Red'], "\n\n\tERROR:";
        print " HEPfit.bib exists. Please choose a name for the output file or move HEPfit.bib.\n\n";
        &usage;
        exit(1)
    } else {
        print colored ['Yellow'], "\tWARNING:";
        print " Output bib file being set to HEPfit.bib.\n\n";
        $bibfile_out = "HEPfit.bib";
        open BIBOUT, "+>", $bibfile_out;
        $end_index = 0;
    }
} elsif ($ARGV[-2] eq "-o" && (-e $ARGV[-1])){
    print colored ['Red'], "\n\n\tERROR:";
    print " Output file exists. Please choose another name or move the existing file.\n\n";
    exit(1)
} else {
    $bibfile_out = pop @ARGV;
    if (index($bibfile_out, ".bib") != -1){
        open BIBOUT, "+>", $bibfile_out;
        $end_index = 1;
    } else {
        print colored ['Red'], "\n\n\tERROR:";
        print " Please list only .bib files in the correct BibTex format.\n";
        &usage;
        exit(1)
    }
}

$end_index = 1 if ($ARGV[-1] eq "-of");

for(my $i=0; $i < (@ARGV - $end_index); $i++) {
    if (index($ARGV[$i], ".bib") != -1) {
        my $bibfile = $ARGV[$i];
        print "\t--> $bibfile added to list\n" if (push(@bibfiles,$bibfile));
    } else {
        print colored ['Red'], "\n\n\tERROR:";
        print " Please list only .bib files in the correct BibTex format.\n";
        &usage;
        exit(1)
    }
}

open BIBOUT_TEMP, ">", "temp.bib";

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
        $_ =~ s/\"\[/\"\{/ if (!(/note/));
            $_ =~ s/\]\"/\}\"/ if (!(/note/));
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
system("cp $bibfile_out temp.bib");
system("sed -e \'s\/\\\$\/\\\]MJ/g\' temp.bib > $bibfile_out"); # A Fistful of dollars get left over no matter what!!

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
    unlink(@tmpfiles);
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
    
    if ($doxyfile_specified ne "-dox"){
        print colored ['Green'], "\n\tDo you want to use the default Doxyfile (DEFAULT: yes)\n\t";
        chomp(my $doxyfile_asked = <STDIN>);
        $doxyfile_asked = "yes" if ($doxyfile_asked eq "");
        while (($doxyfile_asked ne "yes") && ($doxyfile_asked ne "no")){
            print colored ['Red'], "\n\tERROR:";
            print " illegal option: $doxyfile_asked";
            print "\tPlease enter yes or no.\n\t";
            $doxyfile_asked = <STDIN>;
        }
        $doxyfile_path = "Doxyfile";
        if ($doxyfile_asked eq "no"){
            print colored ['Green'], "\n\tPath to the configuration file to be used by Doxygen?\n\t";
            chomp($doxyfile_path = <STDIN>);
            while (!(-e $doxyfile_path)){
                print colored ['Red'], "\n\tERROR:";
                print " Specified Doxygen configuration file does not exist.\n";
                print "\tPlease enter the correct path.\n\t";
                chomp($doxyfile_path = <STDIN>);
            }
        }
    }
    print "\tUsing doxygen configuration file $doxyfile_path.\n";
    
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
            unlink(@tmpfiles);
            exit(1)
        }
    } else {
        print "\tDoxygen building on existing documentation.\n";
    }
    print colored ['Green'], "\n\tRunning Doxygen...\n";
    system(join(" ", $doxygen_used, $doxyfile_path));
} else {
    print "\tNO DOXYGEN RUN.\n";
}

print colored ['Green'], "\n\tPatching citelist.html...\n";

chomp(my $citelist_path = `find ./html -name "citelist.html"`);
if (!(-e $citelist_path)){
    print colored ['Red'], "\n\tERROR:";
    print " citelist.html file not found recursively in the html directory.\n";
    unlink(@tmpfiles);
    exit(1)
} else {
    print "\tcitelist.html found at: $citelist_path\n";
}

if(!(copy($citelist_path,"citelist_temp.html"))){
    print colored ['Red'], "\n\tERROR:";
    print " citelist.html file not copied to the current directory.\n";
    unlink(@tmpfiles);
    exit(1)
}
open CITELISTIN, "<citelist_temp.html";
open CITELISTOUT, ">citelist.html";
push(@tmpfiles,"citelist_temp.html");
while (<CITELISTIN>){
    if (/target\=blank/){
        print colored ['Red'], "\n\tERROR:";
        print " patched version of citelist.html exists. Please rerun Doxygen.\n";
        unlink(@tmpfiles);
        exit(1)
    }
    $_ =~ s/([A-Z])\./$1\.\&\#160/g;
    $_ =~ s/\&\#160\&\#160/\&\#160/g;
    $_ =~ s/MJ\[/\\\(/g;
    $_ =~ s/\]MJ/\\\)/g;
    $_ =~ s/\[BS\]/\\/g;
    $_ =~ s/(arXiv)\:([0-9]+\.[0-9]+)/\<a href=\"http\:\/\/inspirehep\.net\/search\?ln\=en\&p\=$1\%3A$2\" style\=\"color\:#00009C\" target\=blank\>$1\:$2<\/a>/;
    $_ =~ s/(arXiv:)(hep-.+)\/([0-9]+)/\<a href=\"http\:\/\/inspirehep\.net\/search\?ln\=en\&p\=$2\%2F$3\" style\=\"color\:#00009C\" target\=blank>$1$2\/$3<\/a>/;
    print CITELISTOUT $_;
}

close CITELISTIN;
close CITELISTOUT;
move("citelist.html", $citelist_path);

print colored ['Green'], "\n\tPatching html/index.html...\n";
open(INFILE,"<MainPage.md") || die "cannot open MainPage.md!";
my $MainPageHeading = <INFILE>;
close(INFILE);
open(INFILE,"<html/index.html") || die "cannot open html/index.html!";
open(OUTFILE,">index_new.html") || die "cannot open index_new.html!";
my $MainPageHeadingORG = "HEPfit Documentation";
while(<INFILE>) {
    $_ =~ s/$MainPageHeadingORG/$MainPageHeading/;
    print OUTFILE $_;
}
print "\thtml/index.html has been modified\n";
move("index_new.html", "html/index.html");
close(INFILE);
close(OUTFILE);

if(!(copy("images/Model_inherit_graph.svg", "html/"))){
    print colored ['Red'], "\n\tERROR:";
    print " images/Model_inherit_graph.svg file not found.\n";
    unlink(@tmpfiles);
    exit(1)
}

print colored ['Green'], "\n\tPatching _page_models.html...\n";
chomp(my $page_model_path = `find ./html -name "_page_models.html"`);
if (!(-e $page_model_path)){
    print colored ['Red'], "\n\tERROR:";
    print " _page_models.html file not found recursively in the html directory.\n";
    unlink(@tmpfiles);
    exit(1)
} else {
    print "\t_page_models.html found at: $page_model_path\n";
}
if(!(copy($page_model_path,"_page_models_temp.html"))){
    print colored ['Red'], "\n\tERROR:";
    print " _page_models.html file not copied to the current directory.\n";
    unlink(@tmpfiles);
    exit(1)
}
open PAGEMODELIN, "<_page_models_temp.html";
open PAGEMODELOUT, ">_page_models.html";
push(@tmpfiles,"_page_models_temp.html");
while (<PAGEMODELIN>){
    if (/Model\_inherit\_graph\.svg/){
        print colored ['Red'], "\n\tERROR:";
        print " patched version of _page_models_temp.html exists. Please rerun Doxygen.\n";
        unlink(@tmpfiles);
        exit(1)
    }
    $_ =~ s/\<p\>MODEL\_GRAPH\_INHERITE\_SVG\<\/p\>/\<div class\=\"center\"\>\<iframe scrolling\=\"no\" frameborder\=\"0\" src\=\"\.\.\/\.\.\/Model\_inherit\_graph\.svg\" width\=\"544\" height\=\"396\"\>This browser is not able to show SVG\: try Firefox, Chrome, Safari, or Opera instead\.\<\/iframe\>\<\/div\>/;
    print PAGEMODELOUT $_;
}

close PAGEMODELIN;
close PAGEMODELOUT;
move("_page_models.html", $page_model_path);

unlink(@tmpfiles);
print colored ['Green'], "\n\tDONE!\n\n";




exit(0);
