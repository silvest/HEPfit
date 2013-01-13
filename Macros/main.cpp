/* 
 * Copyright (C) 2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <streambuf>
#include <algorithm>
#include <TApplication.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TString.h>
#include <TObject.h>
#include <TF1.h>
#include <TLegend.h>
#include "BaseMacros.h"
#include "SFH1D.h"
#include "Pull.h"
#include "SFH2D.h"

using namespace std;

/**
 *  @addtogroup Macros
 *  A project for drawing histograms and outputting results. 
 *  @{
 */

//-- set probability ranges --
//const double prob68 = 0.6827;
//const double prob95 = 0.9545;
const double prob68 = 0.68;
const double prob95 = 0.95;

/**
 * @brief 
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
int main(int argc, char** argv) 
{
    
    if (argc < 3) {
        cout << "#################################################################" << endl;
        cout << " Usage:                                                          " << endl;
        cout << "   1-D histogram:                                                " << endl;
        cout << "     > macros rootfile plotname --oneDim <optional parameters>   " << endl;
        cout << "   Compatibility plot:                                           " << endl;
        cout << "     > macros rootfile plotname --compat <optional parameters>   " << endl;
        cout << "   2-D histogram:                                                " << endl;
        cout << "     > macros rootfile plotname --twoDim <optional parameters>   " << endl;
        cout << "       'rootfile' (with extension)                               " << endl;
        cout << "       'plotname' (e.g. sin2b)                                   " << endl;
        cout << "                                                                 " << endl;
        cout << "   Check the contents of a rootfile:                             " << endl;
        cout << "     > macros rootfile info                                      " << endl;
        cout << "                                                                 " << endl;
        cout << " Optional parameters for 1-D histograms:                         " << endl;
        cout << "   --oneDim         -> 1-D histogram (mandatory)                 " << endl;
        cout << "   --orig           -> superimpose the original histogram        " << endl;
        cout << "   --leftLegend     -> put the legend in the left area           " << endl;
        cout << "   --outputTxt      -> output results to a text file             " << endl;
        cout << "   -output=filename -> the base name of output eps and text files (without extension)" << endl;
        cout << "                       (default: plotname)                       " << endl;
        cout << "   -maxDigits=n     -> set max digits in the axis labels         " << endl;
        cout << "                       (default: n=8)                            " << endl;
        cout << "   -precision=n     -> precision of values in the std output     " << endl;
        cout << "                       (default: n=6)                            " << endl;
        cout << "   -xlab=namex      -> x label                                   " << endl;
        cout << "   -ylab=namey      -> y label                                   " << endl;
        cout << "   -addtext=text    -> attach additional information             " << endl;
        cout << "   -addtextAt=\"[x,y]\" -> position of the text                  " << endl;
        cout << "   -smooth=ntime    -> iterative smoothing with TH1::Smooth()    " << endl;
        cout << "                       (default: ntime=0)                        " << endl;
        cout << "   -moreBins=nbin   -> increase the number of bins               " << endl;
        cout << "                       (default: nbin=100)                       " << endl;
        cout << "   -col68=index     -> color index of the 68% interval           " << endl;
        cout << "                       (default: index=1393)                     " << endl;
        cout << "   -col95=index     -> color index of the 95% interval           " << endl;
        cout << "                       (default: index=1392)                     " << endl;
        cout << "   -fill=index      -> index of the fill area style              " << endl;
        cout << "                       (default: index=1001)                     " << endl;
        cout << "   -leg=legend      -> legend for the histogram                  " << endl;
        cout << "                       (default: no legend)                      " << endl;
        cout << "   -legScale=scale  -> scale factor for the legend               " << endl;
        cout << "                       (default: scale=1.0)                      " << endl;
        cout << "   *** superimpose the second histogram ***                      " << endl;
        cout << "   -plot2=name      -> name of the histogram                     " << endl;
        cout << "   -rootfile2=filename2 -> rootfile name for plot2 (with extension)" << endl;
        cout << "                           (default: same as rootfile)           " << endl;
        cout << "   -smooth2=ntime   -> iterative smoothing for plot2             " << endl;
        cout << "                       (default: ntime=0)                        " << endl;
        cout << "   -moreBins2=nbin  -> increase the number of bins for plot2     " << endl;
        cout << "                       (default: nbin=100)                       " << endl;
        cout << "   -col682=index    -> color index of the 68% interval for plot2 " << endl;
        cout << "                       (default: index=2)                        " << endl;
        cout << "   -col952=index    -> color index of the 95% interval for plot2 " << endl;
        cout << "                       (default: index=5)                        " << endl;
        cout << "   -fill2=index     -> index of the fill area style              " << endl;
        cout << "                       (default: index=1001)                     " << endl;
        cout << "   -leg2=legend     -> legend for plot2                          " << endl;
        cout << "                       (default: no legend)                      " << endl;
        cout << "   *** superimpose a Gaussian function ***                       " << endl;
        cout << "   -priorMean=mean  -> mean value for the Gaussian function      " << endl;
        cout << "   -priorSigma=sig  -> standard deviation for the Gaussian function" << endl;
        cout << "                                                                 " << endl;
        cout << "   -leg3=legend     -> legend for the Gaussian function          " << endl;
        cout << "                       (default: no legend)                      " << endl;
        cout << " Optional parameters for compatibility plots:                    " << endl;
        cout << "   --compat           -> Compatibility plot (mandatory)          " << endl;
        cout << "   --outputTxt      -> output results to a text file             " << endl;
        cout << "   -output=filename -> the base name of output eps and text files (without extension)" << endl;
        cout << "                       (default: plotname)                       " << endl;
        cout << "   -maxDigits=n     -> set max digits in the axis labels         " << endl;
        cout << "                       (default: n=8)                            " << endl;
        cout << "   -precision=n     -> precision of values in the std output     " << endl;
        cout << "                       (default: n=6)                            " << endl;
        cout << "   -addtext=text    -> attach additional information             " << endl;
        cout << "   -addtextAt=\"[x,y]\" -> position of the text                  " << endl;
        cout << "   -range=\"[xmin,xmax]x[ymin,ymax]\" -> define the graph range  " << endl;
        cout << "   -bins=\"[xbins]x[ybins]\"          -> define the graph binning" << endl;
        cout << "                                         (default: 100x20)       " << endl;
        cout << "   -xlab=namex      -> x label                                   " << endl;
        cout << "   -ylab=namey      -> y label                                   " << endl;
        cout << "   *** put a measured point ***                                  " << endl;
        cout << "   -val=xval        -> xval is the measured point                " << endl;
        cout << "   -err=xerr        -> xerr is the error of xval                 " << endl;
        cout << "                                                                 " << endl;
        cout << " Optional parameters for 2-D histograms:                         " << endl;
        cout << "   --twoDim         -> 2-D histogram (mandatory)                 " << endl;
        cout << "   --drawlines      -> draw contour lines                        " << endl;
        cout << "   --outputTxt      -> output results to a text file             " << endl;
        cout << "   --only951        -> draw only the 95% contour                 " << endl;
        cout << "   -output=filename -> the base name of output eps and text files (without extension)" << endl;
        cout << "                       (default: plotname)                       " << endl;
        cout << "   -maxDigits=n     -> set max digits in the axis labels         " << endl;
        cout << "                       (default: n=8)                            " << endl;
        cout << "   -precision=n     -> precision of values in the std output     " << endl;
        cout << "                       (default: n=6)                            " << endl;
        cout << "   -addtext=text    -> attach additional information             " << endl;
        cout << "   -addtextAt=\"[x,y]\" -> position of the text                  " << endl;
        cout << "   -range=\"[xmin,xmax]x[ymin,ymax]\" -> define the graph range  " << endl;
        cout << "   -xlab=namex      -> x label                                   " << endl;
        cout << "   -ylab=namey      -> y label                                   " << endl;
        cout << "   -smooth=ntime    -> iterative smoothing with TH2::Smooth()    " << endl;
        cout << "                       (default: ntime=0)                        " << endl;
        cout << "   -col68=index     -> color index of the 68% interval           " << endl;
        cout << "                       (default: index=1393)                     " << endl;
        cout << "   -col95=index     -> color index of the 95% interval           " << endl;
        cout << "                       (default: index=1392)                     " << endl;
        cout << "   -line=index      -> index of the line style                   " << endl;
        cout << "                       (default: index=1)                        " << endl;
        cout << "   -fill=index      -> index of the fill area style              " << endl;
        cout << "                       (default: index=1001)                     " << endl;
        cout << "   -leg=legend      -> legend for the histogram                  " << endl;
        cout << "                       (default: no legend)                      " << endl;
        cout << "   -legScale=scale  -> scale factor for the legend               " << endl;
        cout << "                       (default: scale=1.0)                      " << endl;
        cout << "   *** superimpose the second histogram ***                      " << endl;
        cout << "   --only952        -> draw only the 95% contour                 " << endl;
        cout << "   -plot2=name      -> name of the histogram                     " << endl;
        cout << "   -rootfile2=filename -> rootfile name for plot2 (with extension)" << endl;
        cout << "                          (default: same as rootfile)            " << endl;
        cout << "   -smooth2=ntime   -> iterative smoothing for plot2             " << endl;
        cout << "                       (default: ntime=0)                        " << endl;
        cout << "   -col682=index    -> color index of the 68% interval for plot2 " << endl;
        cout << "                       (default: index=2)                        " << endl;
        cout << "   -col952=index    -> color index of the 95% interval for plot2 " << endl;
        cout << "                       (default: index=5)                        " << endl;
        cout << "   -line2=index     -> index of the line style                   " << endl;
        cout << "                       (default: index=1)                        " << endl;
        cout << "   -fill2=index     -> index of the fill area style              " << endl;
        cout << "                       (default: index=1001)                     " << endl;
        cout << "   -leg2=legend     -> legend for plot2                          " << endl;
        cout << "                       (default: no legend)                      " << endl;
        cout << "   *** superimpose the third histogram ***                       " << endl;
        cout << "   --only953        -> draw only the 95% contour                 " << endl;
        cout << "   -plot3=name      -> name of the histogram                     " << endl;
        cout << "   -rootfile3=filename -> rootfile name for plot3 (with extension)" << endl;
        cout << "                          (default: same as rootfile)            " << endl;
        cout << "   -smooth3=ntime   -> iterative smoothing for plot3             " << endl;
        cout << "                       (default: ntime=0)                        " << endl;
        cout << "   -col683=index    -> color index of the 68% interval for plot3 " << endl;
        cout << "                       (default: index=2)                        " << endl;
        cout << "   -col953=index    -> color index of the 95% interval for plot3 " << endl;
        cout << "                       (default: index=5)                        " << endl;
        cout << "   -line3=index     -> index of the line style                   " << endl;
        cout << "                       (default: index=1)                        " << endl;
        cout << "   -fill3=index     -> index of the fill area style              " << endl;
        cout << "                       (default: index=1001)                     " << endl;
        cout << "   -leg3=legend     -> legend for plot3                          " << endl;
        cout << "                       (default: no legend)                      " << endl;
        cout << "   *** superimpose the fourth histogram ***                      " << endl;
        cout << "   --only954        -> draw only the 95% contour                 " << endl;
        cout << "   -plot4=name      -> name of the histogram                     " << endl;
        cout << "   -rootfile4=filename -> rootfile name for plot4 (with extension)" << endl;
        cout << "                          (default: same as rootfile)            " << endl;
        cout << "   -smooth4=ntime   -> iterative smoothing for plot4             " << endl;
        cout << "                       (default: ntime=0)                        " << endl;
        cout << "   -col684=index    -> color index of the 68% interval for plot4 " << endl;
        cout << "                       (default: index=2)                        " << endl;
        cout << "   -col954=index    -> color index of the 95% interval for plot4 " << endl;
        cout << "                       (default: index=5)                        " << endl;
        cout << "   -line4=index     -> index of the line style                   " << endl;
        cout << "                       (default: index=1)                        " << endl;
        cout << "   -fill4=index     -> index of the fill area style              " << endl;
        cout << "                       (default: index=1001)                     " << endl;
        cout << "   -leg4=legend     -> legend for plot4                          " << endl;
        cout << "                       (default: no legend)                      " << endl;
        cout << "   *** put a measured point ***                                  " << endl;
        cout << "   -Xval=xval2      -> xval2 is the x value of the point         " << endl;
        cout << "   -Xerr=xerr2      -> xerr2 is the error of xval2               " << endl;
        cout << "   -Yval=yval2      -> yval2 is the y value of the point         " << endl;
        cout << "   -Yerr=yerr2      -> yerr2 is the error of yval2               " << endl;
        cout << "   -colP=index      -> color index of the measured point         " << endl;
        cout << "                       (default: index=1)                        " << endl;
        cout << "   -legP=legend     -> legend for the point                      " << endl;
        cout << "                       (default: no legend)                      " << endl;
        cout << "#################################################################" << endl;
        return 0;
    }

    // parameters which can be changed by command-line arguments
    bool bOneDim = false, bCompat = false, bTwoDim = false;
    bool bOrig = false, bOutputTxt = false, bContLines = false, bLeftLegend = false;
    int maxDig = 8, prec = 6;
    int nx = 100, ny = 20;
    double xval = -999.0, xerr = 0.0, x_low = 0.0, x_up = 0.0, y_low = 0.0, y_up = 0.0;
    double prior_mean = 0.0, prior_sigma = 0.0;
    double legend_scale = 1.0;
    TString addtext = "";
    double addtext_x = 0.0, addtext_y = 0.0;
    TString xlab = "", ylab = "";
    
    double xval2 = -999.0, xerr2 = 0.0, yval2 = -999.0, yerr2 = 0.0;
    TString legP="";
    int colP = 1;
    
    bool bOnly95 = false;
    string plotname = "", filename = "";
    TString leg="";
    int smooth = 0, col68 = 1393, col95 = 1392, fillStyle = 1001, lineStyle = 1, newNbins = 100 ;
    //
    bool b2ndplot = false, bOnly952 = false;
    string plotname2 = "", filename2 = "";
    TString leg2 = "";
    int smooth2 = 0, col682 = 2, col952 = 5, fillStyle2 = 1001, lineStyle2 = 1, newNbins2 = 100;
    //
    bool b3rdplot = false, bOnly953 = false;
    string plotname3 = "", filename3 = "";
    TString leg3 = "";
    int smooth3 = 0, col683 = 2, col953 = 5, fillStyle3 = 1001, lineStyle3 = 1;
    //
    bool b4thplot = false, bOnly954 = false;
    string plotname4 = "", filename4 = "";
    TString leg4 = "";
    int smooth4 = 0, col684 = 2, col954 = 5, fillStyle4 = 1001, lineStyle4 = 1;

    // root file
    filename = argv[1];
    TFile *datafile = new TFile(filename.c_str());
    if (argc == 3 && strncmp(argv[2], "info", 4) == 0) {
        datafile->ls();
        return 0;
    }
    filename2 = filename;
    filename3 = filename;
    filename4 = filename;
    
    // plot name
    plotname = argv[2];
    TObject* tobj;
    tobj = datafile->Get(plotname.c_str());
    if (tobj == NULL) {
        cout << "Error: plot \"" << plotname << "\" does not exist..." << endl;
        return 1;
    }
    
    // options
    string outputFileName = plotname;
    for (int i = 3; i < argc; i++) {
        char str[100];
        if (strncmp(argv[i], "--oneDim", 8) == 0) bOneDim = true;
        if (strncmp(argv[i], "--compat", 8) == 0) bCompat = true;
        if (strncmp(argv[i], "--twoDim", 8) == 0) bTwoDim = true;    
        if (strncmp(argv[i], "--orig", 6) == 0) bOrig = true;
        if (strncmp(argv[i], "--outputTxt", 11) == 0) bOutputTxt = true;        
        if (strncmp(argv[i], "--drawlines", 11) == 0) bContLines = true;        
        if (strncmp(argv[i], "--leftLegend", 12) == 0) bLeftLegend = true;        
        if (strncmp(argv[i], "--only951", 9) == 0) bOnly95 = true;
        if (strncmp(argv[i], "--only952", 9) == 0) bOnly952 = true; 
        if (strncmp(argv[i], "--only953", 9) == 0) bOnly953 = true;
        if (strncmp(argv[i], "--only954", 9) == 0) bOnly954 = true;
        
        if (strncmp(argv[i], "-output=", 8) == 0) {
            sscanf(argv[i], "-output=%s", str);
            outputFileName = str;
        }

        if (strncmp(argv[i], "-xlab=", 6) == 0) {
            sscanf(argv[i], "-xlab=%s", str);
            xlab = str;
        }
        
        if (strncmp(argv[i], "-ylab=", 6) == 0) {
            sscanf(argv[i], "-ylab=%s", str);
            ylab = str;
        }        
        
        if (strncmp(argv[i], "-addtext=", 9) == 0) {
            sscanf(argv[i], "-addtext=%s", str);
            addtext = str;
        }
        if (strncmp(argv[i], "-addtextAt=", 11) == 0)
            sscanf(argv[i], "-addtextAt=[%lf,%lf]", &addtext_x, &addtext_y);
        
        if (strncmp(argv[i], "-maxDigits=", 11) == 0) 
            sscanf(argv[i], "-maxDigits=%d", &maxDig);

        if (strncmp(argv[i], "-precision=", 11) == 0) 
            sscanf(argv[i], "-precision=%d", &prec);
        
        if (strncmp(argv[i], "-priorMean=", 11) == 0) 
            sscanf(argv[i], "-priorMean=%lf", &prior_mean);
        if (strncmp(argv[i], "-priorSigma=", 12) == 0) 
            sscanf(argv[i], "-priorSigma=%lf", &prior_sigma);

        if (strncmp(argv[i], "-range=", 7) == 0) {
            TString stmp(argv[i] + 7);
            sscanf(stmp.Data(), "[%lf,%lf]x[%lf,%lf]", &x_low, &x_up, &y_low, &y_up);
        }

        if (strncmp(argv[i], "-bins=", 6) == 0) {
            TString stmp(argv[i] + 6);
            sscanf(stmp.Data(), "[%d]x[%d]", &nx, &ny);
        }

        if (strncmp(argv[i], "-val=", 5) == 0) 
            sscanf(argv[i], "-val=%lf", &xval);
        if (strncmp(argv[i], "-err=", 5) == 0) 
            sscanf(argv[i], "-err=%lf", &xerr);
        
        if (strncmp(argv[i], "-Xval=", 6) == 0) 
            sscanf(argv[i], "-Xval=%lf", &xval2);
        if (strncmp(argv[i], "-Xerr=", 6) == 0) 
            sscanf(argv[i], "-Xerr=%lf", &xerr2);
        if (strncmp(argv[i], "-Yval=", 6) == 0) 
            sscanf(argv[i], "-Yval=%lf", &yval2);
        if (strncmp(argv[i], "-Yerr=", 6) == 0) 
            sscanf(argv[i], "-Yerr=%lf", &yerr2);
        if (strncmp(argv[i], "-colP=", 6) == 0) 
            sscanf(argv[i], "-colP=%d", &colP);

        if (strncmp(argv[i], "-moreBins=", 10) == 0)
            sscanf(argv[i], "-moreBins=%d", &newNbins);
        if (strncmp(argv[i], "-moreBins2=", 11) == 0)
            sscanf(argv[i], "-moreBins2=%d", &newNbins2);
        
        if (strncmp(argv[i], "-smooth=", 8) == 0)
            sscanf(argv[i], "-smooth=%d", &smooth);
        if (strncmp(argv[i], "-smooth2=", 9) == 0) 
            sscanf(argv[i], "-smooth2=%d", &smooth2);
        if (strncmp(argv[i], "-smooth3=", 9) == 0) 
            sscanf(argv[i], "-smooth3=%d", &smooth3);
        if (strncmp(argv[i], "-smooth4=", 9) == 0) 
            sscanf(argv[i], "-smooth4=%d", &smooth4);

        if (strncmp(argv[i], "-col68=", 7) == 0) 
            sscanf(argv[i], "-col68=%d", &col68);
        if (strncmp(argv[i], "-col95=", 7) == 0) 
            sscanf(argv[i], "-col95=%d", &col95);
        if (strncmp(argv[i], "-col682=", 8) == 0) 
            sscanf(argv[i], "-col682=%d", &col682);
        if (strncmp(argv[i], "-col952=", 8) == 0) 
            sscanf(argv[i], "-col952=%d", &col952);
        if (strncmp(argv[i], "-col683=", 8) == 0) 
            sscanf(argv[i], "-col683=%d", &col683);
        if (strncmp(argv[i], "-col953=", 8) == 0) 
            sscanf(argv[i], "-col953=%d", &col953);
        if (strncmp(argv[i], "-col684=", 8) == 0) 
            sscanf(argv[i], "-col684=%d", &col684);
        if (strncmp(argv[i], "-col954=", 8) == 0) 
            sscanf(argv[i], "-col954=%d", &col954);

        if (strncmp(argv[i], "-line=", 6) == 0) 
            sscanf(argv[i], "-line=%d", &lineStyle);
        if (strncmp(argv[i], "-line2=", 7) == 0) 
            sscanf(argv[i], "-line2=%d", &lineStyle2);
        if (strncmp(argv[i], "-line3=", 7) == 0) 
            sscanf(argv[i], "-line3=%d", &lineStyle3);
        if (strncmp(argv[i], "-line4=", 7) == 0) 
            sscanf(argv[i], "-line4=%d", &lineStyle4);
        
        if (strncmp(argv[i], "-fill=", 6) == 0) 
            sscanf(argv[i], "-fill=%d", &fillStyle);
        if (strncmp(argv[i], "-fill2=", 7) == 0) 
            sscanf(argv[i], "-fill2=%d", &fillStyle2);
        if (strncmp(argv[i], "-fill3=", 7) == 0) 
            sscanf(argv[i], "-fill3=%d", &fillStyle3);
        if (strncmp(argv[i], "-fill4=", 7) == 0) 
            sscanf(argv[i], "-fill4=%d", &fillStyle4);
        
        if (strncmp(argv[i], "-plot2=", 7) == 0) {
            sscanf(argv[i], "-plot2=%s", str);
            plotname2 = str;
            b2ndplot = true;
        }
        if (strncmp(argv[i], "-plot3=", 7) == 0) {
            sscanf(argv[i], "-plot3=%s", str);
            plotname3 = str;
            b3rdplot = true;
        }
        if (strncmp(argv[i], "-plot4=", 7) == 0) {
            sscanf(argv[i], "-plot4=%s", str);
            plotname4 = str;
            b4thplot = true;
        }
        
        if (strncmp(argv[i], "-rootfile2=", 11) == 0) {
            sscanf(argv[i], "-rootfile2=%s", str);
            filename2 = str;
        }        
        if (strncmp(argv[i], "-rootfile3=", 11) == 0) {
            sscanf(argv[i], "-rootfile3=%s", str);
            filename3 = str;
        }        
        if (strncmp(argv[i], "-rootfile4=", 11) == 0) {
            sscanf(argv[i], "-rootfile4=%s", str);
            filename4 = str;
        }        

        if (strncmp(argv[i], "-legScale=", 10) == 0) 
            sscanf(argv[i], "-legScale=%lf", &legend_scale);
        if (strncmp(argv[i], "-leg=", 5) == 0) {
            sscanf(argv[i], "-leg=%s", str);
            leg = str;
        }
        if (strncmp(argv[i], "-leg2=", 6) == 0) {
            sscanf(argv[i], "-leg2=%s", str);
            leg2 = str;
        }
        if (strncmp(argv[i], "-leg3=", 6) == 0) {
            sscanf(argv[i], "-leg3=%s", str);
            leg3 = str;
        }
        if (strncmp(argv[i], "-leg4=", 6) == 0) {
            sscanf(argv[i], "-leg4=%s", str);
            leg4 = str;
        }
        if (strncmp(argv[i], "-legP=", 6) == 0) {
            sscanf(argv[i], "-legP=%s", str);
            legP = str;
        }
    }        

    if (!(bOneDim | bCompat | bTwoDim)) {
        cout << "Error: A mandatory option --oneDim, --compat or --twoDim is missing!" << endl;
        return 1;
    }
    if ( (bOneDim&bCompat) | (bOneDim&bTwoDim) | (bTwoDim&bCompat) ) {
        cout << "Error: Use only one of --oneDim, --compat and --twoDim!" << endl;
        return 1;
    }

    // output files
    string epsFileName = outputFileName + ".eps";
    string txtFileName = outputFileName + ".txt";

    // 2rd plot
    TFile *datafile2;
    TObject* tobj2;
    if (b2ndplot) {
        datafile2 = new TFile(filename2.c_str());
        tobj2 = datafile2->Get(plotname2.c_str());
        if (tobj2 == NULL && b2ndplot) {
            cout << "Error: plot \"" << plotname2 << "\" does not exist..." << endl;
            return 1;
        }
    }
    
    // 3rd plot
    TFile *datafile3;
    TObject* tobj3;
    if (b3rdplot) {
        datafile3 = new TFile(filename3.c_str());
        tobj3 = datafile3->Get(plotname3.c_str());
        if (tobj3 == NULL && b3rdplot) {
            cout << "Error: plot \"" << plotname3 << "\" does not exist..." << endl;
            return 1;
        }
    }

    // 4th plot
    TFile *datafile4;
    TObject* tobj4;
    if (b4thplot) {
        datafile4 = new TFile(filename4.c_str());
        tobj4 = datafile4->Get(plotname4.c_str());
        if (tobj4 == NULL && b4thplot) {
            cout << "Error: plot \"" << plotname4 << "\" does not exist..." << endl;
            return 1;
        }
    }

    //----------------------------------------------------------------------

    ofstream* fout;
    streambuf* buf;
    if (bOutputTxt) {
        fout = new ofstream();
        fout->open(txtFileName.c_str(), ios::out);
        if (!fout->is_open())
            cout << "cannot open " << txtFileName << endl;
        buf = fout->rdbuf();        
    } else
        buf = cout.rdbuf();
    ostream os(buf);
    os.precision(prec);

    BaseMacros::DefineNewColours();
    BaseMacros myMacros;
    
    TCanvas TC("TC", "TC", 3);
    TC.SetLeftMargin(0.20);
    // TC.SetRightMargin(0.20);
    TC.SetBottomMargin(0.15);
    TC.SetFillColor(0);
    TC.SetBorderMode(0);
    TC.SetFrameFillColor(0);
    TC.SetHighLightColor(0);

    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    //gStyle->SetStripDecimals(false);   

    TLegend *legend;
    int num_leg = 0;
    if ( leg.CompareTo("")!= 0 ) num_leg++;
    if ( leg2.CompareTo("")!= 0 ) num_leg++;
    if ( leg3.CompareTo("")!= 0 ) num_leg++;
    if ( leg4.CompareTo("")!= 0 ) num_leg++;
    if ( legP.CompareTo("")!= 0 ) num_leg++;
    double legend_ymin;
    if (num_leg == 1) legend_ymin = 0.80;
    else if (num_leg == 2) legend_ymin = 0.74;
    else if (num_leg == 3) legend_ymin = 0.68;
    else if (num_leg == 4) legend_ymin = 0.62;
    else if (num_leg == 5) legend_ymin = 0.56;
    else legend_ymin = 0.47;
    legend_ymin = (0.88 - 0.02)*(1.0 - legend_scale) + legend_ymin*legend_scale;
    if (!bLeftLegend)
        legend = new TLegend(0.63,legend_ymin,0.75,0.88);
    else
        //legend = new TLegend(0.23,legend_ymin,0.48,0.88);
        legend = new TLegend(0.23,legend_ymin,0.35,0.88);
    legend->SetFillColor(0);
    legend->SetBorderSize(0);
    legend->SetTextFont(42);
    legend->SetTextSize(0.043*legend_scale);
    legend->SetMargin(0.7);

    //----------------------------------------------------------------------

    if (bOneDim) {
        // 1-D histogram
        TH1D* hist = (TH1D*) tobj->Clone();   
        os << hist->GetXaxis()->GetTitle() << " in " << plotname << endl;

        // output the 1-D histogram
        SFH1D SFHisto1D(*hist, prob68, prob95);
        SFHisto1D.smoothHist(smooth);
        SFHisto1D.increaseNbins(newNbins);
        SFHisto1D.Draw(xlab, ylab, col68, col95, fillStyle, maxDig, bOrig, false); 
    
        // rescale (for mHl)
        //SFHisto1D.getNewHist()->Scale(10.0);
        //SFHisto1D.getNewHist68()->Scale(10.0);
        //SFHisto1D.getNewHist95()->Scale(10.0);
        //SFHisto1D.getNewHist()->GetXaxis()->SetRange(400,1400);
        //leg += " [x10]";
            
        // superimpose another 1-D histogram (e.g. for a posterior)
        TH1D* plot2_pt = NULL;
        if (b2ndplot) {
            TH1D* hist2 = (TH1D*) tobj2->Clone();
            SFH1D SFHisto1D2(*hist2, prob68, prob95);
            SFHisto1D2.smoothHist(smooth2);
            SFHisto1D2.increaseNbins(newNbins2);
            SFHisto1D2.Draw("", "", col682, col952, fillStyle2, maxDig, bOrig, true);    
            plot2_pt = SFHisto1D2.getNewHist68();

            // rescale the y axis
            double ymax_new = max( SFHisto1D.getNewHist()->GetMaximum(), 
                                   1.1*SFHisto1D2.getNewHist()->GetMaximum() );
            SFHisto1D.getNewHist()->SetMaximum(ymax_new);
        }
            
        // superimpose a Gaussian (prior) function
        TF1* prior = NULL;
        if (prior_sigma != 0.0) {
            prior = new TF1("prior",
                    "1./sqrt(2.*TMath::Pi())/[1]* exp(- (x-[0])*(x-[0])/2./[1]/[1])",
                    SFHisto1D.getNewHist()->GetXaxis()->GetXmin(),
                    SFHisto1D.getNewHist()->GetXaxis()->GetXmax());
            prior->SetParameter(0, prior_mean);
            prior->SetParameter(1, prior_sigma);    
            prior->SetLineStyle(2);
            prior->SetLineWidth(4);
            prior->SetNpx(1000);
            prior->Draw("SAME");
            
            // rescale the y axis
            double ymax_new = max( SFHisto1D.getNewHist()->GetMaximum(),
                                   1.1*prior->GetMaximum() );
            SFHisto1D.getNewHist()->SetMaximum(ymax_new);
        }

        // legends: Change the order if needed. 
        if (prior != NULL) legend->AddEntry(prior, myMacros.ConvertTitle(leg3), "L");
        if (plot2_pt != NULL) legend->AddEntry(plot2_pt, myMacros.ConvertTitle(leg2), "F");
        legend->AddEntry(SFHisto1D.getNewHist68(), myMacros.ConvertTitle(leg), "F");
        
        // output results to os
        SFHisto1D.OutputResults(os, smooth, true);

        // output results for the original histogram before smoothing 
        // nor increasing the number of bins for comparison
        SFH1D orig1D(*hist, prob68, prob95);
        os << endl << "[Original histogram]" << endl;
        orig1D.OutputResults(os, 0, false);
        
    } else if (bCompat) {
        // 1-D histogram
        TH1D* hist = (TH1D*) tobj->Clone();
        os << hist->GetXaxis()->GetTitle() << " in " << plotname << endl;
        os << "  Num of bins: " << nx << " x " << ny << endl;
        
        // output the compatibility plot
        Pull CompatPlot(*hist, nx, ny, x_low, x_up, y_low, y_up);
        CompatPlot.Draw(xlab, ylab, xval, xerr, maxDig);
        
    } else if (bTwoDim) {    
        // 2-D histogram
        TH2D* hist = (TH2D*) tobj->Clone();
        os << hist->GetXaxis()->GetTitle() << " vs " << hist->GetYaxis()->GetTitle() 
           << " in " << plotname << endl;
        os << "  smooth: " << smooth << " time(s)" << endl;
        
        // output the 2-D histogram
        SFH2D SFHisto2D(*hist, os, prob68, prob95, x_low, x_up, y_low, y_up);
        SFHisto2D.smoothHist(smooth);
        SFHisto2D.draw(xlab, ylab, col68, col95, lineStyle, fillStyle, maxDig, 
                       bContLines, bOnly95, false);
        TGraph* contour_pt = SFHisto2D.getContour();
        
        // superimpose the 2nd histogram
        TGraph* contour2_pt = NULL;
        if (b2ndplot) {
            os << "[Second graph]" << endl;
            os << "  smooth: " << smooth2 << " time(s)" << endl;
            TH2D* hist2 = (TH2D*) tobj2->Clone();
            SFH2D SFHisto2D2(*hist2, os, prob68, prob95);
            SFHisto2D2.smoothHist(smooth2);
            SFHisto2D2.draw("", "", col682, col952, lineStyle2, fillStyle2, maxDig, 
                            bContLines, bOnly952, true);
            contour2_pt = SFHisto2D2.getContour();
        }

        // superimpose the 3rd histogram
        TGraph* contour3_pt = NULL;
        if (b3rdplot) {
            os << "[Third graph]" << endl;
            os << "  smooth: " << smooth3 << " time(s)" << endl;
            TH2D* hist3 = (TH2D*) tobj3->Clone();
            SFH2D SFHisto2D3(*hist3, os, prob68, prob95);
            SFHisto2D3.smoothHist(smooth3);
            SFHisto2D3.draw("", "", col683, col953, lineStyle3, fillStyle3, maxDig, 
                            bContLines, bOnly953, true);
            contour3_pt = SFHisto2D3.getContour();
        }
        
        // superimpose the 4th histogram
        TGraph* contour4_pt = NULL;
        SFH2D* SFHisto2D4 = NULL;
        if (b4thplot) {
            os << "[Fourth graph]" << endl;
            os << "  smooth: " << smooth4 << " time(s)" << endl;
            TH2D* hist4 = (TH2D*) tobj4->Clone();
            SFHisto2D4 = new SFH2D(*hist4, os, prob68, prob95);
            SFHisto2D4->smoothHist(smooth4);
            SFHisto2D4->draw("", "", col684, col954, lineStyle4, fillStyle4, maxDig, 
                             bContLines, bOnly954, true);
            contour4_pt = SFHisto2D4->getContour();
        }
        
        // draw a given point with error bars
        TGraphErrors* g1 = NULL;
        if (xval2 != -999.0 && yval2 != -999.0) {
            TLine *lx = new TLine();
            lx->SetLineWidth(3);
            lx->SetLineColor(colP);
            double zero = 0, err;
            err = xerr2;
            g1 = new TGraphErrors(1, &xval2, &yval2, &err, &zero);
            g1->SetLineWidth(3);
            g1->SetLineStyle(1);
            g1->SetLineColor(colP);
            g1->SetMarkerStyle(20);
            g1->SetMarkerSize(1);
            g1->SetMarkerColor(colP);
            g1->Draw("P");
            
            double xmin = hist->GetXaxis()->GetXmin();
            double xmax = hist->GetXaxis()->GetXmax();
            double ymin = hist->GetYaxis()->GetXmin();
            double ymax = hist->GetYaxis()->GetXmax();
            
            err = xerr2;
            double min_x = std::max(xval2 - err, xmin);
            double max_x = std::min(xval2 + err, xmax);
            lx->DrawLine(min_x, yval2, max_x, yval2);
            
            TLine *ly = new TLine();
            ly->SetLineWidth(3);
            ly->SetLineColor(colP);
            err = yerr2;
            TGraphErrors* g2 = new TGraphErrors(1, &xval2, &yval2, &zero, &err);
            g2->SetLineWidth(3);
            g2->SetLineStyle(1);
            g2->SetLineColor(colP);
            g2->SetMarkerStyle(20);
            g2->SetMarkerSize(1);
            g2->SetMarkerColor(colP);
            g2->Draw("P");
            
            double min_y = std::max(yval2 - err, ymin);
            double max_y = std::min(yval2 + err, ymax);
            ly->DrawLine(xval2, min_y, xval2, max_y);
        }
        
        // legends: Change the order if needed. 
        string leg_opt;
        if (contour4_pt != NULL && leg4.CompareTo("")!= 0) {
            if (fillStyle4!=0) leg_opt = "F";
            else leg_opt = "L";
            legend->AddEntry(contour4_pt, myMacros.ConvertTitle(leg4), leg_opt.c_str());
        }
        if (contour3_pt != NULL && leg3.CompareTo("")!= 0) {
            if (fillStyle3!=0) leg_opt = "F";
            else leg_opt = "L";
            legend->AddEntry(contour3_pt, myMacros.ConvertTitle(leg3), leg_opt.c_str()); 
        }
        if (contour2_pt != NULL && leg2.CompareTo("")!= 0) {
            if (fillStyle2!=0) leg_opt = "F";
            else leg_opt = "L";
            legend->AddEntry(contour2_pt, myMacros.ConvertTitle(leg2), leg_opt.c_str());
        }
        if (leg.CompareTo("")!= 0) {
            if (fillStyle!=0) leg_opt = "F";
            else leg_opt = "L";
            legend->AddEntry(contour_pt, myMacros.ConvertTitle(leg), leg_opt.c_str());
        }
        if (g1 != NULL && legP.CompareTo("")!= 0) 
            legend->AddEntry(g1, myMacros.ConvertTitle(legP), "LP");
    } 

    // additional text
    if ( addtext.CompareTo("")!= 0 ) {
        TLatex* TL = new TLatex();
        TL->SetTextSize(0.043);
        TL->SetTextFont(42);
        TL->SetTextAlign(11);
        TL->SetNDC(1);
        TL->DrawLatex(addtext_x, addtext_y, myMacros.ConvertTitle(addtext));
    }
       
    // draw the legend if not empty
    if (num_leg>0) legend->Draw("");
    
    gPad->RedrawAxis();
    TC.Print(epsFileName.c_str());
    
    if (bOutputTxt) {
        fout->close();          
        cout << txtFileName << " has been created." << endl;
    }

    return 0;
}

/** 
 * @}
 */


