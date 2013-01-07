/* 
 * Copyright (C) 2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include <iostream>
#include <TH2D.h>
#include <TStyle.h>
#include <TObjArray.h>
#include <TGraph.h>
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
     */
    SFH2D(TH2D& hist, std::ostream& os_in, 
          const double prob68_in=0.68, const double prob95_in=0.95);

    /**
     * @return A modified 2-D histogram. 
     */
    TH2D* getNewHist() const 
    {
        return newHist;
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
     * @param[in] col68 The color index for the 68% interval. 
     * @param[in] col95 The color index for the 95% interval. 
     * @param[in] maxDigits The maximum digits of axis labels. 
     * @param xval2
     * @param xerr2
     * @param yval2
     * @param yerr2
     * @param[in] bLine  
     * @param[in] superImpose
     */
    void draw(const TString xlab, const TString ylab, 
              const int col68, const int col95, const int maxDigits, 
              const double xval2 = -999.0, const double xerr2 = 0.0,
              const double yval2 = -999.0, const double yerr2 = 0.0,            
              const bool bLine=true, const bool superImpose=false);

        
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
     * @param[in] col
     */
    void drawFromGraph(const int ind, const std::string DrawOpts, const int col) const;    
    
    /**
     * @param[in] inputgraph
     * @return 
     */
    TGraph* CloseTGraph(TGraph* inputgraph) const;
    
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
    
    double m_x;
    double m_y;
    double m_r;
};

#endif	/* SFH2D_H */

