/* 
 * Copyright (C) 2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include <cmath>
#include <vector>
#include <functional>
#include <TPad.h>
#include <TCanvas.h>
#include <TString.h>
#include <TGaxis.h>
#include <TLine.h>
#include <TGraphErrors.h>
#include <TROOT.h>
#include "SFH2D.h"

SFH2D::SFH2D(TH2D& hist, std::ostream& os_in, 
             const double prob68_in, const double prob95_in) 
: os(os_in), prob68(prob68_in), prob95(prob95_in), origHist(hist), 
        origName(hist.GetName())
{
    std::string NewName = origName + "_new";
    newHist = (TH2D*) hist.Clone(NewName.c_str()); 
    double sum = newHist->Integral();
    newHist->Scale(1.0/sum);
}


void SFH2D::smoothHist(const int smooth)
{
    if (smooth!=0) { 
        for (int k = 0; k < smooth; k++)
            newHist->Smooth();
    }
}


void SFH2D::draw(const TString xlab, const TString ylab, 
                 const int col68, const int col95, const int maxDigits, 
                 const double xval2, const double xerr2,
                 const double yval2, const double yerr2,            
                 const bool bLine, const bool superImpose) 
{
    newHist->SetTitle("");
    newHist->SetStats(0);
    newHist->GetXaxis()->SetTitleSize(0.075);
    newHist->GetYaxis()->SetTitleSize(0.075);
    newHist->GetXaxis()->SetTitleOffset(0.85);
    newHist->GetYaxis()->SetTitleOffset(0.90);    
    newHist->GetXaxis()->SetNdivisions(505);
    newHist->GetYaxis()->SetNdivisions(505);
    newHist->GetXaxis()->SetLabelSize(0.043);
    newHist->GetYaxis()->SetLabelSize(0.043);
    newHist->GetXaxis()->SetLabelOffset(0.013);
    newHist->GetYaxis()->SetLabelOffset(0.013);
    newHist->SetLabelFont(42,"X");
    newHist->SetLabelFont(42,"Y");
    newHist->SetTitleFont(42,"X");
    newHist->SetTitleFont(42,"Y");
    //newHist->SetLabelFont(62,"X");
    //newHist->SetLabelFont(62,"Y");
    //newHist->SetTitleFont(62,"X");
    //newHist->SetTitleFont(62,"Y");
    ((TGaxis*) newHist->GetXaxis())->SetMaxDigits(maxDigits);
    ((TGaxis*) newHist->GetYaxis())->SetMaxDigits(maxDigits);

    // Titles of the axes 
    if (xlab.CompareTo("")==0) {
        TString Xtitle = ConvertTitle(newHist->GetXaxis()->GetTitle());
        newHist->GetXaxis()->SetTitle(Xtitle);
    } else
        newHist->GetXaxis()->SetTitle(xlab);
    if (ylab.CompareTo("")==0) {
        TString Ytitle = ConvertTitle(newHist->GetYaxis()->GetTitle());
        newHist->GetYaxis()->SetTitle(Ytitle);     
    } else 
        newHist->GetYaxis()->SetTitle(ylab);    
    
    if (!superImpose) {
        TH2D* null2D = (TH2D*) newHist->Clone("null2D");
        null2D->Reset();
        null2D->Draw();
    }

    drawFromGraph(0, "AREA", col95);
    drawFromGraph(1, "AREA", col68);

    // draw the 95% contour line. 
    if (bLine)
        drawFromGraph(0, "CONT", col68);

    // draw error bars
    if (xval2 != -999.0 && yval2 != -999.0 && !superImpose) {
        TLine *lx = new TLine();
        lx->SetLineWidth(3);
        double zero = 0, err;
        err = xerr2;
        TGraphErrors *g1 = new TGraphErrors(1, &xval2, &yval2, &err, &zero);
        g1->SetLineWidth(3);
        g1->SetLineStyle(1);
        g1->SetMarkerStyle(20);
        g1->SetMarkerSize(1);
        g1->Draw("P");

        double xmin = newHist->GetXaxis()->GetXmin();
        double xmax = newHist->GetXaxis()->GetXmax();
        double ymin = newHist->GetYaxis()->GetXmin();
        double ymax = newHist->GetYaxis()->GetXmax();

        err = xerr2;
        double min_x = std::max(xval2 - err, xmin);
        double max_x = std::min(xval2 + err, xmax);
        lx->DrawLine(min_x, yval2, max_x, yval2);

        TLine *ly = new TLine();
        ly->SetLineWidth(3);
        err = yerr2;
        TGraphErrors* g2 = new TGraphErrors(1, &xval2, &yval2, &zero, &err);
        g2->SetLineWidth(3);
        g2->SetLineStyle(1);
        g2->SetMarkerStyle(20);
        g2->SetMarkerSize(1);
        g2->Draw("P");

        double min_y = std::max(yval2 - err, ymin);
        double max_y = std::min(yval2 + err, ymax);
        ly->DrawLine(xval2, min_y, xval2, max_y);
    }
    
    gPad->RedrawAxis();
}



//////////////////////////////////////////////////////////////////////// 

double SFH2D::getLevel(const double Prob) const 
{
    int nbinx = newHist->GetNbinsX();
    int nbiny = newHist->GetNbinsY();
    
    // normalize the histogram
    double sum = newHist->Integral();
    newHist->Scale(1.0/sum);

    std::vector<double> orderedContents;
    for (int ix = 1; ix <= nbinx; ix++) 
        for (int iy = 1; iy <= nbiny; iy++) 
            orderedContents.push_back(newHist->GetBinContent(ix, iy));
    
    // sort into descending order
    std::sort(orderedContents.begin(), orderedContents.end(), std::greater<double>());

    // get the level where the probability is equal to or a bit larger than Prob.
    double level = 0.0, area = 0.0;
    for (int i = 0; i < (int) orderedContents.size(); i++) {
        area += orderedContents.at(i);
        level = orderedContents.at(i);
        //std::cout << i << " " << level << " " << area << std::endl; // for debug
        if (area >= Prob) break;
    }
    
    return level;
}


TObjArray* SFH2D::getContours() const
{
    // set levels
    double levels[2];
    levels[0] = getLevel(prob95);
    levels[1] = getLevel(prob68);
    
    // get contours
    TCanvas c1("c1", "contours", 3);
    TH2D* tmpHist = (TH2D*) newHist->Clone("tmp");
    tmpHist->SetContour(2, levels);
    tmpHist->Draw("contlist");
    c1.Update(); // Needed to force the plotting and retrieve the contours in TGraphs
    TObjArray *contours = (TObjArray*) gROOT->GetListOfSpecials()->FindObject("contours");

    // output
    os << "Contours:" << std::endl;
    if (contours == NULL)
        os << "  No contours were extracted!" << std::endl;
    os << "  The number of the levels: " << contours->GetSize() << std::endl;
    for(int i = 0; i < contours->GetSize(); i++){
        TList* contLevel = (TList*)contours->At(i);
        os << "  Contour " << i << " for the level " << levels[i]
           << " has " << contLevel->GetSize() << " graph(s)" << std::endl;
    }

    return contours;
}


void SFH2D::drawFromGraph(const int ind, const std::string DrawOpts, 
                          const int col) const 
{
    TList* contour = (TList*) getContours()->At(ind);
    
    TAxis* tmp = (TAxis*) newHist->GetXaxis();
    double xmin = tmp->GetXmin();
    double xmax = tmp->GetXmax();
    tmp = (TAxis*) newHist->GetYaxis();
    double ymin = tmp->GetXmin();
    double ymax = tmp->GetXmax();
    
    double epsp = 1.0;
    double epsm = 1.0 - 1.0e-6;
    double x01 = newHist->GetBinContent(newHist->GetXaxis()->FindBin(xmin * epsp), 
                                        newHist->GetYaxis()->FindBin(ymax * epsm));
    double x11 = newHist->GetBinContent(newHist->GetXaxis()->FindBin(xmax * epsm), 
                                        newHist->GetYaxis()->FindBin(ymax * epsm));
    double x00 = newHist->GetBinContent(newHist->GetXaxis()->FindBin(xmin * epsp), 
                                        newHist->GetYaxis()->FindBin(ymin * epsp));
    double x10 = newHist->GetBinContent(newHist->GetXaxis()->FindBin(xmax * epsm), 
                                        newHist->GetYaxis()->FindBin(ymin * epsp));

    TGraph* curv = (TGraph*) contour->First();
    
    for (int i = 0; i < contour->GetSize(); i++) {

        bool l01 = false;
        bool l11 = false;
        bool l00 = false;
        bool l10 = false;
        double binwx = newHist->GetXaxis()->GetBinWidth(1);
        double binwy = newHist->GetYaxis()->GetBinWidth(1);

        double val00 = -1., val01 = -1., val10 = -1., val11 = -1.;
        double val;
        for (int j = 0; j < curv->GetN(); j++) {
            double xx, yy;
            curv->GetPoint(j, xx, yy);
            if (xx >= xmax - binwx / 2.) curv->SetPoint(j, xmax, yy);
            if (xx <= xmin + binwx / 2.) curv->SetPoint(j, xmin, yy);
            if (yy >= ymax - binwy / 2.) curv->SetPoint(j, xx, ymax);
            if (yy <= ymin + binwy / 2.) curv->SetPoint(j, xx, ymin);
            val = newHist->GetBinContent(newHist->GetXaxis()->FindBin(xx), 
                                         newHist->GetYaxis()->FindBin(yy));
            if (val < x01) {
                l01 = true;
                val01 = val;
            }
            if (val < x00) {
                l00 = true;
                val00 = val;
            }
            if (val < x10) {
                l10 = true;
                val10 = val;
            }
            if (val < x11) {
                l11 = true;
                val11 = val;
            }
        }

        //cout << "--------------" << endl;
        //cout << "x00  " << x00 << endl;
        //cout << "x01  " << x01 << endl;
        //cout << "x10  " << x10 << endl;
        //cout << "x11  " << x11 << endl;
        //cout << "val  " << val << endl;
        //cout << "--------------" << endl;
        
        std::vector<SFH2D_Point> vp;
        double xl, yl;
        curv->GetPoint(curv->GetN() - 1, xl, yl);
        if (l01) {
            SFH2D_Point p(xmin, ymax);
            p.R(sqrt(pow(xmin - xl, 2) + pow(ymax - yl, 2)));
            vp.push_back(p);
        }
        if (l11) {
            SFH2D_Point p(xmax, ymax);
            p.R(sqrt(pow(xmax - xl, 2) + pow(ymax - yl, 2)));
            vp.push_back(p);
        }
        if (l10) {
            SFH2D_Point p(xmax, ymin);
            p.R(sqrt(pow(xmax - xl, 2) + pow(ymin - yl, 2)));
            vp.push_back(p);
        }
        if (l00) {
            SFH2D_Point p(xmin, ymin);
            p.R(sqrt(pow(xmin - xl, 2) + pow(ymin - yl, 2)));
            vp.push_back(p);
        }
        std::sort(vp.begin(), vp.end());
        
        if (DrawOpts == "AREA") {
            for (unsigned int i = 0; i < vp.size(); i++) {
                curv->SetPoint(curv->GetN(), vp[i].m_x, vp[i].m_y);
            }
        }

        TGraph* tmpTGraph = CloseTGraph(curv);
        //TGraph* tmpTGraph = (TGraph*) curv->Clone(); // test
        if (DrawOpts.compare("AREA") == 0) {
            tmpTGraph->SetLineWidth(2);
            tmpTGraph->SetLineColor(col);
            tmpTGraph->SetFillColor(col);
            tmpTGraph->Draw("F");
        }

        // necessary to ensure that the contour lines of plot1 survive 
        // after superimposing plot2 
        TGraph* curv2 = (TGraph*) curv->Clone(); 
        
        if (DrawOpts.compare("CONT") == 0) {
            curv2->SetLineWidth(2);
            curv2->SetLineColor(col);
            curv2->Draw();
        }
        curv = (TGraph*) contour->After(curv); // Get Next graph
    }

    return;
}


TGraph* SFH2D::CloseTGraph(TGraph* inputgraph) const 
{

    TAxis* tmp = (TAxis*) newHist->GetXaxis();
    double xmin = tmp->GetXmin();
    double xmax = tmp->GetXmax();
    double binx = tmp->GetNbins();

    tmp = (TAxis*) newHist->GetYaxis();
    double ymin = tmp->GetXmin();
    double etamin = ymin;
    double ymax = tmp->GetXmax();    
    double biny = tmp->GetNbins();
    
    double x_i, x_j, y_i, y_j;
    inputgraph->GetPoint(0, x_i, y_i);
    inputgraph->GetPoint(inputgraph->GetN() - 1, x_j, y_j);

    double xleft_i(0), xright_i(0);
    double ybottom_i(0), ytop_i(0);
    double xleft_j(0), xright_j(0);
    double ybottom_j(0), ytop_j(0);

    double deltax = (xmax - xmin) / binx;
    double deltay = (ymax - etamin) / biny;

    if (fabs(x_i - xmin) < deltax) xleft_i = 1.;
    if (fabs(x_i - xmax) < deltax) xright_i = 1.;
    if (fabs(y_i - etamin) < deltay) ybottom_i = 1.;
    if (fabs(y_i - ymax) < deltay) ytop_i = 1.;

    if (fabs(x_j - xmin) < deltax) xleft_j = 1.;
    if (fabs(x_j - xmax) < deltax) xright_j = 1.;
    if (fabs(y_j - etamin) < deltay) ybottom_j = 1.;
    if (fabs(y_j - ymax) < deltay) ytop_j = 1.;

    double xnew[inputgraph->GetN() + 3];
    double ynew[inputgraph->GetN() + 3];

    for (int i = 0; i < inputgraph->GetN(); i++) {
        inputgraph->GetPoint(i, xnew[i], ynew[i]);
    }

    if (xleft_i == 1. && ybottom_j == 1.) {
        // we go from bottom to left
        // passing through the left-bottom corner
        xnew[inputgraph->GetN()] = x_j;
        ynew[inputgraph->GetN()] = etamin;
        xnew[inputgraph->GetN() + 1] = xmin;
        ynew[inputgraph->GetN() + 1] = etamin;
        xnew[inputgraph->GetN() + 2] = xmin;
        ynew[inputgraph->GetN() + 2] = y_i;
    } else if (xleft_j == 1. && ybottom_i == 1.) {
        // we go from left to bottom
        // passing through the left-bottom corner
        xnew[inputgraph->GetN()] = xmin;
        ynew[inputgraph->GetN()] = y_j;
        xnew[inputgraph->GetN() + 1] = xmin;
        ynew[inputgraph->GetN() + 1] = etamin;
        xnew[inputgraph->GetN() + 2] = x_i;
        ynew[inputgraph->GetN() + 2] = etamin;
    } else if (xleft_j == 1. && ytop_i == 1.) {
        // we go from left to top
        // passing through the left-top corner
        xnew[inputgraph->GetN()] = xmin;
        ynew[inputgraph->GetN()] = y_j;
        xnew[inputgraph->GetN() + 1] = xmin;
        ynew[inputgraph->GetN() + 1] = ymax;
        xnew[inputgraph->GetN() + 2] = x_i;
        ynew[inputgraph->GetN() + 2] = ymax;
    } else if (xleft_i == 1. && ytop_j == 1.) {
        // we go from left to top
        // passing through the left-top corner
        //xnew[inputgraph->GetN()]   = x_j;    ynew[inputgraph->GetN()]   = ymax;
        //xnew[inputgraph->GetN()+1] = xmin;   ynew[inputgraph->GetN()+1] = ymax;
        //xnew[inputgraph->GetN()+2] = xmin;   ynew[inputgraph->GetN()+2] = y_i;
    } else if (xright_i == 1. && ytop_j == 1.) {
        // we go from right to top
        // passing through the tight-top corner
        xnew[inputgraph->GetN()] = x_j;
        ynew[inputgraph->GetN()] = ymax;
        xnew[inputgraph->GetN() + 1] = xmax;
        ynew[inputgraph->GetN() + 1] = ymax;
        xnew[inputgraph->GetN() + 2] = xmax;
        ynew[inputgraph->GetN() + 2] = y_i;
    } else if (xright_j == 1. && ytop_i == 1.) {
        // we go from top to right
        // passing through the tight-top corner
        xnew[inputgraph->GetN()] = xmax;
        ynew[inputgraph->GetN()] = y_j;
        xnew[inputgraph->GetN() + 1] = xmax;
        ynew[inputgraph->GetN() + 1] = ymax;
        xnew[inputgraph->GetN() + 2] = x_i;
        ynew[inputgraph->GetN() + 2] = ymax;
    } else if (xright_i == 1. && ybottom_j == 1.) {
        // we go from right to bottom
        // passing through the right-bottom corner
        xnew[inputgraph->GetN()] = x_j;
        ynew[inputgraph->GetN()] = etamin;
        xnew[inputgraph->GetN() + 1] = xmax;
        ynew[inputgraph->GetN() + 1] = etamin;
        xnew[inputgraph->GetN() + 2] = xmax;
        ynew[inputgraph->GetN() + 2] = y_i;
    } else if (xright_j == 1. && ybottom_i == 1.) {
        // we go from bottom to right
        // passing through the right-bottom corner
        xnew[inputgraph->GetN()] = xmax;
        ynew[inputgraph->GetN()] = y_j;
        xnew[inputgraph->GetN() + 1] = xmax;
        ynew[inputgraph->GetN() + 1] = etamin;
        xnew[inputgraph->GetN() + 2] = x_i;
        ynew[inputgraph->GetN() + 2] = etamin;
    } else if (xleft_i == 1. && ybottom_j == 1.) {
        // we go from left to bottom
        // passing through the left-bottom corner
        xnew[inputgraph->GetN()] = x_j;
        ynew[inputgraph->GetN()] = etamin;
        xnew[inputgraph->GetN() + 1] = xmin;
        ynew[inputgraph->GetN() + 1] = etamin;
        xnew[inputgraph->GetN() + 2] = xmin;
        ynew[inputgraph->GetN() + 2] = y_i;
    } else if (xleft_j == 1. && ybottom_i == 1.) {
        // we go from bottom to left
        // passing through the left-bottom corner
        xnew[inputgraph->GetN()] = xmin;
        ynew[inputgraph->GetN()] = y_j;
        xnew[inputgraph->GetN() + 1] = xmin;
        ynew[inputgraph->GetN() + 1] = etamin;
        xnew[inputgraph->GetN() + 2] = x_i;
        ynew[inputgraph->GetN() + 2] = etamin;
        //    NP C_Bs vs phi_Bs
    } else if (ybottom_i == 1. && ytop_j == 1.) { //from top to bottom
        xnew[inputgraph->GetN()] = x_j;
        ynew[inputgraph->GetN()] = ymax;
        xnew[inputgraph->GetN() + 1] = 1.0;
        ynew[inputgraph->GetN() + 1] = ymax;
        xnew[inputgraph->GetN() + 2] = 1.0;
        ynew[inputgraph->GetN() + 2] = etamin;
    } else if (ybottom_j == 1. && ytop_i == 1.) { //from bottom to top
        xnew[inputgraph->GetN()] = x_j;
        ynew[inputgraph->GetN()] = etamin;
        xnew[inputgraph->GetN() + 1] = xmin;
        ynew[inputgraph->GetN() + 1] = etamin;
        xnew[inputgraph->GetN() + 2] = xmin;
        ynew[inputgraph->GetN() + 2] = ymax;
        //xnew[inputgraph->GetN()]   = x_j;    ynew[inputgraph->GetN()]   = etamin;
        //xnew[inputgraph->GetN()+1] = 1.0;    ynew[inputgraph->GetN()+1] = etamin;
        //xnew[inputgraph->GetN()+2] = 1.0;    ynew[inputgraph->GetN()+2] = ymax;

        // achilleplot Bd
    } else if (xleft_i == 1. && xright_j == 1.) { //from left to right
        xnew[inputgraph->GetN()] = xmax;
        ynew[inputgraph->GetN()] = y_j;
        xnew[inputgraph->GetN() + 1] = xmax;
        ynew[inputgraph->GetN() + 1] = etamin;
        xnew[inputgraph->GetN() + 2] = xmin;
        ynew[inputgraph->GetN() + 2] = etamin;
    } else if (xleft_j == 1. && xright_i == 1.) { //from right to left
        xnew[inputgraph->GetN()] = xmin;
        ynew[inputgraph->GetN()] = y_j;
        xnew[inputgraph->GetN() + 1] = xmin;
        ynew[inputgraph->GetN() + 1] = etamin;
        xnew[inputgraph->GetN() + 2] = xmax;
        ynew[inputgraph->GetN() + 2] = etamin;
    } else {
        // nominal
        xnew[inputgraph->GetN()] = x_j;
        ynew[inputgraph->GetN()] = y_j;
        xnew[inputgraph->GetN() + 1] = x_j;
        ynew[inputgraph->GetN() + 1] = y_j;
        xnew[inputgraph->GetN() + 2] = x_j;
        ynew[inputgraph->GetN() + 2] = y_j;
    }

    TGraph* newTGraph = new TGraph(inputgraph->GetN() + 3, xnew, ynew);
    return newTGraph;
}



