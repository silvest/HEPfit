/* 
 * File:   Observable.cpp
 * Author: silvest
 * 
 * Created on February 22, 2011, 11:45 AM
 */

#include "Observable.h"

void Observable::Init(const std::string distr_i, const std::string filename_i,
        const std::string histoname_i, const double ave_i, const double errg_i, const double errf_i) {
    distr = distr_i;
    filename = filename_i;
    histoname = histoname_i;
    ave = ave_i;
    errg = errg_i;
    errf = errf_i;
}

Observable::Observable (const std::string name_i, const bool tMCMC_i, 
        const double min_i, const double max_i, ThObservable * tho_i) {
    name = name_i;
    min = min_i;
    max = max_i;
    tMCMC = tMCMC_i;
    tho = tho_i;
}

Observable::Observable(const Observable& orig) {
    name = orig.name;
    min = orig.min;
    max = orig.max;
    tMCMC = orig.tMCMC;
    tho = orig.tho;
    Init(orig.distr, orig.filename, orig.histoname,orig.ave, orig.errg, orig.errf);
}

Observable::~Observable() {
}

void Observable::Set(const std::string distr_i) {
    if((distr_i.compare("noweight")!=0) || tMCMC) {
        std::cout << "Wrong Observable set called: Set(" <<
                 name <<", " <<  distr_i << ")"<<std::endl;
    }
    Init(distr_i, "", "", 0., 0., 0.);
}

void Observable::Set(const std::string distr_i, const std::string filename_i, const std::string histoname_i) {
    if(distr_i.compare("file")!=0) {
        std::cout << "Wrong Observable set called: Set(" << 
                name <<", "  << distr_i << ", " << filename_i << ", " << histoname_i << ")"<<std::endl;
    }
    Init(distr_i, filename_i, histoname_i, 0., 0., 0.);
}

void Observable::Set(const std::string distr_i, const double ave_i,
    const double errg_i, const double errf_i) {
    if(distr_i.compare("weight")!=0) {
        std::cout << "Wrong Observable set called: Set(" <<
                name <<", "  << distr_i
                << ", " << ave_i <<", " << errg_i <<", " << errf_i <<", "  << ")"<<std::endl;
    }
    Init(distr_i, "", "", ave_i, errg_i, errf_i);
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

