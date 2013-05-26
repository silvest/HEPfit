/* 
 * Copyright (C) 2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include <iostream>
#include <cmath>
#include <TH2D.h>
#include <TStyle.h>
#include <TObjArray.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include "BaseMacros.h"

#ifndef SFH2D_H
#define	SFH2D_H

/**
 * @class SFH2D
 * @ingroup Macros
 * @brief A class for handling a 2-D histogram.
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to draw the smallest intervals containing the 68% 
 * and 95% probability ranges. 
 */
class SFH2D : public BaseMacros  {
public:

    /**
     * A constructor.
     * @param[in] hist A 2-D ROOT histogram 
     * @param[in] os_in The standard output stream. 
     * @param[in] prob68_in The probability of the 68% interval. 
     * @param[in] prob95_in The probability of the 95% interval. 
     * @param[in] x_low
     * @param[in] x_up
     * @param[in] y_low
     * @param[in] y_up
     * @param[in] NumNewPoints_in
     */
    SFH2D(TH2D& hist, std::ostream& os_in, 
          const double prob68_in=0.68, const double prob95_in=0.95,
          const double x_low=0.0, const double x_up=0.0, 
          const double y_low=0.0, const double y_up=0.0,
          const int NumNewPoints_in=20);

    /**
     * @return A modified 2-D histogram. 
     */
    TH2D* getNewHist() const 
    {
        return newHist;
    }
    
    /**
     * @return The last contour drawn by drawFromGraph(). 
     */    
    TGraph* getContour() const 
    {
        return myCurv;
    }
    
    /**
     * Smooth the bin contents of th histogram by using TH2::Smooth(). 
     * @param smooth
     */
    void smoothHist(const int smooth=0);

    /**
     * Draw the modified histogram. 
     * @param[in] xlab The label of the x axis. 
     * @param[in] ylab The label of the y axis. 
     * @param[in] lineWidth The line width. 
     * @param[in] lineColor The index of the line color. 
     * @param[in] col68 The color index for the 68% interval. 
     * @param[in] col95 The color index for the 95% interval. 
     * @param[in] lineStyle The index of the line style. 
     * @param[in] lineStyle68 The index of the line style for the 68% contour line. 
     * @param[in] fillStyle The index of the fill area style. 
     * @param[in] maxDigits The maximum digits of axis labels.
     * @param[in] bLine  
     * @param[in] bOnly95
     * @param[in] superImpose
     * @param[in] YTitleOffset 
     */
    void Draw(const TString xlab, const TString ylab, 
              const int lineWidth, const int lineColor,  
              const int col68, const int col95, const int lineStyle, 
              const int lineStyle68,
              const int fillStyle, const int maxDigits, 
              const bool bLine=true, const bool bOnly95=false, 
              const bool superImpose=false, const double YTitleOffset=1.5);

        
private:
    
    /**
     * The standard output stream. 
     */
    std::ostream& os;

    /**
     * Probability of the 68% interval.
     */
    const double prob68;
    
    /**
     * Probability of the 95% interval. 
     */
    const double prob95;
    
    /**
     * The original 2-D histogram.
     */
    TH2D& origHist;

    /**
     * The name of the original histogram.
     */
    const std::string origName;

    /**
     * A modified histogram.  
     */
    TH2D* newHist;
    
    TGraph* myCurv;
    
    double xLow;
    double xUp;
    double yLow;
    double yUp;    

    int NumNewPoints;

    /**
     * @param[in] Prob A probability.
     * @return The contour level at a given probability. 
     */    
    double getLevel(const double Prob) const;

    /**
     * Get 68% and 95% contours. 
     */
    TObjArray* getContours() const;

    /**
     * 
     * @param[in] ind
     * @param[in] DrawOpts
     * @param[in] lineWidth The line width. 
     * @param[in] lineColor The index of the line color. 
     * @param[in] col
     * @param[in] lineStyle
     * @param[in] fillStyle
     */
    void drawFromGraph(const int ind, const std::string DrawOpts, 
                       const int lineWidth, const int lineColor, 
                       const int col, const int lineStyle, const int fillStyle); 
    
    /**
     * @param[in] inputgraph
     * @return 
     */
    TGraph* CloseTGraph(TGraph* inputgraph) const;

    /**
     * Connect an end point of a contour with an end point of the other contour 
     * in the case where if both points are on the same boundary of the plot. 
     * @param[in] cont_ind
     * @param[in] inputgraph1
     * @param[in] inputgraph2
     * @return 
     */
    TGraph* CloseTwoTGraphs(const int cont_ind, TGraph* inputgraph1, 
                            TGraph* inputgraph2) const;
    
};



/**
 * @class SFH2D_Point
 * @ingroup Macros
 * @brief A supplementary class to SFH2D. 
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class SFH2D_Point {
public:

    SFH2D_Point(double x, double y) 
    {
        m_x = x;
        m_y = y;
    }

    void R(double r) 
    {
        m_r = r;
    }

    bool operator<(const SFH2D_Point& b) const 
    {
        return m_r < b.m_r;
    }
    
    double distance(const SFH2D_Point& b) const
    {
        return sqrt( pow(this->m_x - b.m_x, 2.0) + pow(this->m_y - b.m_y, 2.0) );
    }
    
    double m_x;
    double m_y;
    double m_r;
};



#endif	/* SFH2D_H */

