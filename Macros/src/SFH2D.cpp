/* 
 * Copyright (C) 2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include <vector>
#include <functional>
#include <TPad.h>
#include <TCanvas.h>
#include <TString.h>
#include <TGaxis.h>
#include <TLine.h>
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
    
    myCurv = NULL;
    g1 = NULL;
    
    xLow = 0.0;
    xUp = 0.0;
    yLow = 0.0;
    yUp = 0.0;    
}


void SFH2D::SetRange(const double x_low, const double x_up, 
                     const double y_low, const double y_up) {
    xLow = x_low;
    xUp = x_up;
    yLow = y_low;
    yUp = y_up;
}


void SFH2D::smoothHist(const int smooth)
{
    if (smooth!=0) { 
        for (int k = 0; k < smooth; k++)
            newHist->Smooth();
    }
}


void SFH2D::draw(const TString xlab, const TString ylab, 
                 const int col68, const int col95, const int lineStyle, 
                 const int fillStyle, const int maxDigits, 
                 const double xval2, const double xerr2,
                 const double yval2, const double yerr2,            
                 const bool bLine, const bool bOnly95, const bool superImpose) 
{
    // draw the axes 
    if (!superImpose) {
        TH2D* null2D;
        if (xLow==0.0 && xUp ==0.0 && yLow==0.0 && yUp ==0.0)
            null2D = (TH2D*) newHist->Clone("null2D");
        else {
            null2D = (TH2D*) newHist->Clone("null2D");
            null2D->Reset("M");
            null2D = new TH2D("null2D","null2D", 100, xLow, xUp, 100, yLow, yUp); 
            null2D->SetXTitle(newHist->GetXaxis()->GetTitle());
            null2D->SetYTitle(newHist->GetYaxis()->GetTitle());
        }
        null2D->SetTitle("");
        null2D->SetStats(0);
        null2D->GetXaxis()->SetTitleSize(0.075);
        null2D->GetYaxis()->SetTitleSize(0.075);
        null2D->GetXaxis()->SetTitleOffset(0.85);
        null2D->GetYaxis()->SetTitleOffset(0.90);    
        null2D->GetXaxis()->SetNdivisions(505);
        null2D->GetYaxis()->SetNdivisions(505);
        null2D->GetXaxis()->SetLabelSize(0.043);
        null2D->GetYaxis()->SetLabelSize(0.043);
        null2D->GetXaxis()->SetLabelOffset(0.013);
        null2D->GetYaxis()->SetLabelOffset(0.013);
        null2D->SetLabelFont(42,"X");
        null2D->SetLabelFont(42,"Y");
        null2D->SetTitleFont(42,"X");
        null2D->SetTitleFont(42,"Y");
        //null2D->SetLabelFont(62,"X");
        //null2D->SetLabelFont(62,"Y");
        //null2D->SetTitleFont(62,"X");
        //null2D->SetTitleFont(62,"Y");
        ((TGaxis*) null2D->GetXaxis())->SetMaxDigits(maxDigits);
        ((TGaxis*) null2D->GetYaxis())->SetMaxDigits(maxDigits);
        
        // Titles of the axes 
        if (xlab.CompareTo("")==0) {
            TString Xtitle = ConvertTitle(null2D->GetXaxis()->GetTitle());
            null2D->GetXaxis()->SetTitle(Xtitle);
        } else
            null2D->GetXaxis()->SetTitle(xlab);
        if (ylab.CompareTo("")==0) {
            TString Ytitle = ConvertTitle(null2D->GetYaxis()->GetTitle());
            null2D->GetYaxis()->SetTitle(Ytitle);     
        } else 
            null2D->GetYaxis()->SetTitle(ylab);         
        
        null2D->Reset();
        null2D->Draw();
    }

    drawFromGraph(0, "AREA", col95, lineStyle, fillStyle); // 95%
    if (!bOnly95) drawFromGraph(1, "AREA", col68, lineStyle, fillStyle); // 68%
        
    // draw the contour lines. 
    if (bLine) {
        drawFromGraph(0, "CONT", col68, lineStyle, fillStyle); // 95%
        if (!bOnly95) drawFromGraph(1, "CONT", col68, lineStyle, fillStyle); // 68%
    }
    
    // draw a given point with error bars
    if (xval2 != -999.0 && yval2 != -999.0 && !superImpose) {
        TLine *lx = new TLine();
        lx->SetLineWidth(3);
        double zero = 0, err;
        err = xerr2;
        g1 = new TGraphErrors(1, &xval2, &yval2, &err, &zero);
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

    // get the level where the probability corresponds to a given Prob.
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
    os << " *Contours:" << std::endl;
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
                          const int col, const int lineStyle, const int fillStyle) 
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
    TGraph* curv_new[contour->GetSize()];
    
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
        
        if (DrawOpts == "AREA")
            for (unsigned int i = 0; i < vp.size(); i++)
                curv->SetPoint(curv->GetN(), vp[i].m_x, vp[i].m_y);

        if (contour->GetSize()==2)
            curv_new[i] = (TGraph*) curv->Clone();
        else {
            curv_new[i] = CloseTGraph(curv);
            curv_new[i]->SetLineWidth(2);
            curv_new[i]->SetLineColor(col);            
            curv_new[i]->SetLineStyle(lineStyle);
            curv_new[i]->SetFillColor(col);
            curv_new[i]->SetFillStyle(fillStyle);
            if (DrawOpts.compare("AREA") == 0 && fillStyle!=0)
                curv_new[i]->Draw("F");
            else if (DrawOpts.compare("CONT") == 0)
                curv_new[i]->Draw();
        }
        myCurv = curv_new[i];
        curv = (TGraph*) contour->After(curv); // Get Next graph
    }

    if (contour->GetSize()==2) {
        TGraph* curv_combined = CloseTwoTGraphs(curv_new[0], curv_new[1]); 
        curv_combined->SetLineWidth(2);
        curv_combined->SetLineColor(col);
        curv_combined->SetLineStyle(lineStyle);
        curv_combined->SetFillColor(col);
        curv_combined->SetFillStyle(fillStyle);
        if (DrawOpts.compare("AREA") == 0 && fillStyle!=0) {
            curv_combined->Draw("F");
        } else if (DrawOpts.compare("CONT") == 0)
            curv_combined->Draw();
        myCurv = curv_combined;
    }
    
    return;
}


TGraph* SFH2D::CloseTGraph(TGraph* inputgraph) const 
{
    TAxis* tmp = (TAxis*) newHist->GetXaxis();
    double xmin = tmp->GetXmin();
    double xmax = tmp->GetXmax();
    double deltax = tmp->GetBinWidth(1);

    tmp = (TAxis*) newHist->GetYaxis();
    double ymin = tmp->GetXmin();
    double ymax = tmp->GetXmax();    
    double biny = tmp->GetNbins();
    double deltay = tmp->GetBinWidth(1);

    // get the end points of the contour line
    double x_i, x_j, y_i, y_j;
    inputgraph->GetPoint(0, x_i, y_i);
    inputgraph->GetPoint(inputgraph->GetN() - 1, x_j, y_j);
    //std::cout << "  (x_i,y_i)=(" << x_i << "," << y_i 
    //          << "), (x_j,y_j)=(" << x_j << "," << y_j << ")" << std::endl;

    // check if the point (x_i,y_i) is on the boundary
    bool xleft_i = false, xright_i = false, ybottom_i = false, ytop_i = false;
    if (fabs(x_i - xmin) < deltax) xleft_i = true;
    if (fabs(x_i - xmax) < deltax) xright_i = true;
    if (fabs(y_i - ymin) < deltay) ybottom_i = true;
    if (fabs(y_i - ymax) < deltay) ytop_i = true;

    // check if the point (x_j,y_j) is on the boundary
    bool xleft_j = false, xright_j = false, ybottom_j = false, ytop_j = false;
    if (fabs(x_j - xmin) < deltax) xleft_j = true;
    if (fabs(x_j - xmax) < deltax) xright_j = true;
    if (fabs(y_j - ymin) < deltay) ybottom_j = true;
    if (fabs(y_j - ymax) < deltay) ytop_j = true;
    //std::cout << "  i[LRBT]: " << xleft_i << " " << xright_i << " " << ybottom_i << " " << ytop_i
    //          << "  j[LRBT]: " << xleft_j << " " << xright_j << " " << ybottom_j << " " << ytop_j << std::endl;
    
    // get all the points of the contour line
    double xnew[inputgraph->GetN() + 3];
    double ynew[inputgraph->GetN() + 3];
    for (int i = 0; i < inputgraph->GetN(); i++)
        inputgraph->GetPoint(i, xnew[i], ynew[i]);

    // set additional 3 points to connect the end points of the contour line
    int n_start = inputgraph->GetN();
    int n_middle = inputgraph->GetN() + 1;
    int n_last = inputgraph->GetN() + 2;
    if (xleft_i && ybottom_j) {
        // we go from bottom to left
        // passing through the left-bottom corner
        xnew[n_start] = x_j;   ynew[n_start] = ymin;
        xnew[n_middle] = xmin; ynew[n_middle] = ymin;
        xnew[n_last] = xmin;   ynew[n_last] = y_i;
    } else if (xleft_j && ybottom_i) {
        // we go from left to bottom
        // passing through the left-bottom corner
        xnew[n_start] = xmin;  ynew[n_start] = y_j;
        xnew[n_middle] = xmin; ynew[n_middle] = ymin;
        xnew[n_last] = x_i;    ynew[n_last] = ymin;
    } else if (xleft_j && ytop_i) {
        // we go from left to top
        // passing through the left-top corner
        xnew[n_start] = xmin;  ynew[n_start] = y_j;
        xnew[n_middle] = xmin; ynew[n_middle] = ymax;
        xnew[n_last] = x_i;    ynew[n_last] = ymax;
    } else if (xleft_i && ytop_j) {
        // we go from left to top
        // passing through the left-top corner
        xnew[n_start]   = x_j; ynew[n_start]   = ymax;
        xnew[n_middle] = xmin; ynew[n_middle] = ymax;
        xnew[n_last] = xmin;   ynew[n_last] = y_i;
    } else if (xright_i && ytop_j) {
        // we go from right to top
        // passing through the tight-top corner
        xnew[n_start] = x_j;   ynew[n_start] = ymax;
        xnew[n_middle] = xmax; ynew[n_middle] = ymax;
        xnew[n_last] = xmax;   ynew[n_last] = y_i;
    } else if (xright_j && ytop_i) {
        // we go from top to right
        // passing through the tight-top corner
        xnew[n_start] = xmax;  ynew[n_start] = y_j;
        xnew[n_middle] = xmax; ynew[n_middle] = ymax;
        xnew[n_last] = x_i;    ynew[n_last] = ymax;
    } else if (xright_i && ybottom_j) {
        // we go from right to bottom
        // passing through the right-bottom corner
        xnew[n_start] = x_j;   ynew[n_start] = ymin;
        xnew[n_middle] = xmax; ynew[n_middle] = ymin;
        xnew[n_last] = xmax;   ynew[n_last] = y_i;
    } else if (xright_j && ybottom_i) {
        // we go from bottom to right
        // passing through the right-bottom corner
        xnew[n_start] = xmax;  ynew[n_start] = y_j;
        xnew[n_middle] = xmax; ynew[n_middle] = ymin;
        xnew[n_last] = x_i;    ynew[n_last] = ymin;
    } else if (ybottom_i && ytop_j) { //from top to bottom
        xnew[n_start] = x_j;  ynew[n_start] = ymax;
        xnew[n_middle] = xmax; ynew[n_middle] = ymax; // the top-right corner
        xnew[n_last] = xmax;   ynew[n_last] = ymin;   // the bottom-right corner
    } else if (ybottom_j && ytop_i) { //from bottom to top
        xnew[n_start] = x_j;   ynew[n_start] = ymin;
        xnew[n_middle] = xmin; ynew[n_middle] = ymin; // the bottom-left corner
        xnew[n_last] = xmin;   ynew[n_last] = ymax;   // the top-left corner
    } else if (xleft_i && xright_j) { //from left to right
        xnew[n_start] = xmax;  ynew[n_start] = y_j;
        xnew[n_middle] = xmax; ynew[n_middle] = ymin; // the bottom-right corner
        xnew[n_last] = xmin;   ynew[n_last] = ymin;   // the bottom-left corner
    } else if (xleft_j && xright_i) { //from right to left
        xnew[n_start] = xmin;  ynew[n_start] = y_j;
        xnew[n_middle] = xmin; ynew[n_middle] = ymax; // the top-left corner
        xnew[n_last] = xmax;   ynew[n_last] = ymax;   // the top-right corner
    } else {
        // nominal
        xnew[n_start] = x_j;  ynew[n_start] = y_j;
        xnew[n_middle] = x_j; ynew[n_middle] = y_j;
        xnew[n_last] = x_j;   ynew[n_last] = y_j;
    }

    //std::cout << "  Additional points: ";
    //for (int i = inputgraph->GetN(); i < inputgraph->GetN() + 3; i++)
    //    std::cout << "(" << xnew[i] << "," << ynew[i] << ") ";
    //std::cout << std::endl;
    
    // add 3 additional points to TGraph to connect the end point of the contour
    TGraph* newTGraph = new TGraph(inputgraph->GetN() + 3, xnew, ynew);
    return newTGraph;
}


TGraph* SFH2D::CloseTwoTGraphs(TGraph* inputgraph1, TGraph* inputgraph2) const
{
    double xmin = newHist->GetXaxis()->GetXmin();
    double xmax = newHist->GetXaxis()->GetXmax();
    double ymin = newHist->GetYaxis()->GetXmin();
    double ymax = newHist->GetYaxis()->GetXmax();    

    // the end points
    std::vector<SFH2D_Point> p_xmin, p_xmax, p_ymin, p_ymax;
    
    // get all the points of the contour lines    
    int n_all = inputgraph1->GetN() + inputgraph2->GetN();
    std::vector<SFH2D_Point> vp_org;
    for (int i = 0; i < n_all; i++) {
        double xtmp, ytmp;
        if (i < inputgraph1->GetN() ) inputgraph1->GetPoint(i, xtmp, ytmp);
        else inputgraph2->GetPoint(i - inputgraph1->GetN(), xtmp, ytmp);
        SFH2D_Point p(xtmp, ytmp);
        vp_org.push_back(p);

        // store the index of an end point
        if (xtmp==xmin) p_xmin.push_back(p);
        if (xtmp==xmax) p_xmax.push_back(p);
        if (ytmp==ymin) p_ymin.push_back(p);
        if (ytmp==ymax) p_ymax.push_back(p);
    }

    // add NP points in the middle of an end point of the first contour and 
    // that of the second contour
    int NP = 20;
    if (p_xmin.size()==2) {
        double xtmp, ytmp;
        for (int i=0; i<NP; i++) {
            xtmp = p_xmin.at(1).m_x + (p_xmin.at(0).m_x - p_xmin.at(1).m_x)/(double)NP*(double)i;
            ytmp = p_xmin.at(1).m_y + (p_xmin.at(0).m_y - p_xmin.at(1).m_y)/(double)NP*(double)i;
            SFH2D_Point p(xtmp, ytmp);
            vp_org.push_back(p);
        }
    }
    if (p_xmax.size()==2) {
        double xtmp, ytmp;
        for (int i=0; i<NP; i++) {
            xtmp = p_xmax.at(1).m_x + (p_xmax.at(0).m_x - p_xmax.at(1).m_x)/(double)NP*(double)i;
            ytmp = p_xmax.at(1).m_y + (p_xmax.at(0).m_y - p_xmax.at(1).m_y)/(double)NP*(double)i;
            SFH2D_Point p(xtmp, ytmp);
            vp_org.push_back(p);
        }
    }
    if (p_ymin.size()==2) {
        double xtmp, ytmp;
        for (int i=0; i<NP; i++) {
            xtmp = p_ymin.at(1).m_x + (p_ymin.at(0).m_x - p_ymin.at(1).m_x)/(double)NP*(double)i;
            ytmp = p_ymin.at(1).m_y + (p_ymin.at(0).m_y - p_ymin.at(1).m_y)/(double)NP*(double)i;
            SFH2D_Point p(xtmp, ytmp);
            vp_org.push_back(p);
        }
    }
    if (p_ymax.size()==2) {
        double xtmp, ytmp;
        for (int i=0; i<NP; i++) {
            xtmp = p_ymax.at(1).m_x + (p_ymax.at(0).m_x - p_ymax.at(1).m_x)/(double)NP*(double)i;
            ytmp = p_ymax.at(1).m_y + (p_ymax.at(0).m_y - p_ymax.at(1).m_y)/(double)NP*(double)i;
            SFH2D_Point p(xtmp, ytmp);
            vp_org.push_back(p);
        }
    }
    
    std::vector<SFH2D_Point> vp_new;
    std::vector<SFH2D_Point>::iterator it, it_minimal;
    int ind, ind_minimal;
    double dist, dist_tmp;
    double dist_max = sqrt( pow(xmax - xmin, 2.0) + pow(ymax - ymin, 2.0) );

    // set the first point
    it = vp_org.begin();
    ind = std::distance(vp_org.begin(), it);
    vp_new.push_back(vp_org.at(ind));
    vp_org.erase(it);

    while (vp_org.size()>0) {
        dist = dist_max;
        for (it = vp_org.begin(); it < vp_org.end(); it++) {
            ind = std::distance(vp_org.begin(), it);
            dist_tmp = vp_new.back().distance(vp_org.at(ind));
            if (dist_tmp < dist) {
                it_minimal = it;
                ind_minimal = ind;
                dist = dist_tmp;
            }
        }
        vp_new.push_back(vp_org.at(ind_minimal));
        vp_org.erase(it_minimal);
        //std::cout << vp_org.size() << " " << vp_new.size() << " " << ind_minimal 
        //          << " " << dist << std::endl; //debug
    }

    // Debug
    //for (int i = 0; i < vp_new.size(); i++) {
    //    int ind;
    //    if(i==0) ind = 0;
    //    else ind = i-1;        
    //    std::cout << vp_new.at(i).m_x << " " << vp_new.at(i).m_y << " "
    //              << vp_new.at(i).distance(vp_new.at(ind)) << std::endl;
    //}
    
    TGraph* newTGraph = new TGraph(n_all);
    for (int i = 0; i < vp_new.size(); i++)
        newTGraph->SetPoint(i, vp_new.at(i).m_x, vp_new.at(i).m_y);
    
    return newTGraph;
}




