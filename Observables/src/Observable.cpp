/* 
 * File:   Observable.cpp
 * Author: silvest
 * 
 * Created on February 22, 2011, 11:45 AM
 */

#include "Observable.h"

void Observable::Init(const bool tMCMC_i, const double min_i,
    const double max_i, const std::string distr_i, const std::string filename_i,
        const double ave_i, const double errg_i, const double errf_i) {
    min = min_i;
    max = max_i;
    tMCMC = tMCMC_i;
    distr = distr_i;
    filename = filename_i;
    ave = ave_i;
    errg = errg_i;
    errf = errf_i;
}

Observable::Observable (const std::string name_i) {
    name = name_i;
}


Observable::Observable(const Observable& orig) {
	name = orig.name;
    Init(orig.tMCMC, orig.min, orig.max, orig.distr, orig.filename,
            orig.ave, orig.errg, orig.errf);
}

Observable::~Observable() {
}

void Observable::Set(const bool tMCMC_i, const double min_i,
    const double max_i, const std::string distr_i) {
    if((distr_i.compare("noweight")!=0) || tMCMC_i) {
        std::cout << "Wrong Observable set called: Set(" <<
                 name << tMCMC_i <<", " << min_i <<", " << max_i <<", " << distr_i << ")"<<std::endl;
    }
    Init(tMCMC_i, min_i, max_i, distr_i, "", 0., 0., 0.);
}

void Observable::Set(const bool tMCMC_i, const double min_i,
    const double max_i, const std::string distr_i, const std::string filename_i) {
    if(distr_i.compare("file")!=0) {
        std::cout << "Wrong Observable set called: Set(" << 
                name <<", " << tMCMC_i <<", " << min_i <<", " << max_i <<", " << distr_i << ", " << filename_i << ")"<<std::endl;
    }
    Init(tMCMC_i, min_i, max_i, distr_i, filename_i, 0., 0., 0.);
}

void Observable::Set(const bool tMCMC_i, const double min_i,
    const double max_i, const std::string distr_i, const double ave_i,
    const double errg_i, const double errf_i) {
    if(distr_i.compare("weight")!=0) {
        std::cout << "Wrong Observable set called: Set(" <<
                name <<", " << tMCMC_i <<", " << min_i <<", " << max_i <<", " << distr_i
                << ", " << ave_i <<", " << errg_i <<", " << errf_i <<", "  << ")"<<std::endl;
    }
    Init(tMCMC_i, min_i, max_i, distr_i, "", ave_i, errg_i, errf_i);
}

std::ostream& operator<<(std::ostream& output, const Observable& o)
  {
    output << "Observable name, tMCMC, min, max, distribution, distribution parameters" << std::endl;
    output << o.name << " " << o.tMCMC << " " << o.min << " " << o.max << " " 
            << o.distr << " " << o.filename << " " << o.ave << 
            " " << o.errg << " " << o.errf << std::endl;
    return output;
  }
