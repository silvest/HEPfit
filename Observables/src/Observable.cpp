/* 
 * Copyright (C) 2012 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "Observable.h"

Observable::Observable (const std::string name_i, const std::string thname_i,
        const std::string label_i, const bool tMCMC_i, const double min_i,
        const double max_i, ThObservable * tho_i) {
    name = name_i;
    thname = thname_i;
    label = label_i;
    min = min_i;
    max = max_i;
    tMCMC = tMCMC_i;
    tho = tho_i;
    distr = "";
    filename = "";
    histoname = "";
    ave = 0.;
    errg = 0.;
    errf = 0.;
}

Observable::Observable(const Observable& orig) {
    name = orig.name;
    thname = orig.thname;
    label = orig.label;
    min = orig.min;
    max = orig.max;
    tMCMC = orig.tMCMC;
    tho = orig.tho;
    distr = orig.distr; 
    filename = orig.filename; 
    histoname = orig.histoname;
    ave = orig.ave; 
    errg = orig.errg; 
    errf = orig.errf;
}

Observable::~Observable() {
}

std::ostream& operator<<(std::ostream& output, const Observable& o)
  {
    output << "Observable name, tMCMC, min, max, distribution, distribution parameters" << std::endl;
    output << o.name << " " << o.tMCMC << " " << o.min << " " << o.max << " " 
            << o.distr << " " << o.filename << " " << o.histoname << " " << o.ave << 
            " " << o.errg << " " << o.errf << std::endl;
    return output;
  }

double Observable::getTheoryValue(){
    return tho->getThValue();
}

