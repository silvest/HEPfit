/* 
 * Copyright (C) 2012 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "Observable2D.h"

Observable2D::Observable2D(const std::string name_i, const std::string thname_i,
        const std::string thname2_i, const std::string label_i, 
        const std::string label2_i, const bool tMCMC_i, const double min_i,
        const double max_i, const double min2_i, const double max2_i, 
        ThObservable * tho_i, ThObservable * tho2_i) : Observable (name_i, thname_i,
        label_i, tMCMC_i, min_i, max_i, tho_i){
    thname2 = thname2_i;
    label2 = label2_i;
    min2 = min2_i;
    max2 = max2_i;
    tho2 = tho2_i;
}

Observable2D::Observable2D(const Observable& o1d) :  Observable (o1d){
    thname2 = "";
    label2 = "";
    min2 = 0.;
    max2 = 0.;
    tho2 = NULL;   
}

Observable2D::Observable2D(const Observable2D& orig) :  Observable (orig.name, orig.thname,
        orig.label, orig.tMCMC, orig.min, orig.max, orig.tho){
   
    distr = orig.distr; 
    filename = orig.filename; 
    histoname = orig.histoname;
    ave = orig.ave; 
    errg = orig.errg; 
    errf = orig.errf;
    
    thname2 = orig.thname2;
    label2 = orig.label2;
    min2 = orig.min2;
    max2 = orig.max2;
    tho2 = orig.tho2;   
}

Observable2D::~Observable2D() {
}

double Observable2D::getTheoryValue2(){
    return tho2->getThValue();
}
