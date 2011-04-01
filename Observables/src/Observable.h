/* 
 * File:   Observable.h
 * Author: silvest
 *
 * Created on February 22, 2011, 11:45 AM
 */

#ifndef OBSERVABLE_H
#define	OBSERVABLE_H

#include <string>
#include <iostream>
#include "ThObservable.h"

class Observable {
public:
    Observable(const std::string name_i, const bool tMCMC_i, 
        const double min_i, const double max_i, ThObservable * tho_i);
    Observable(const Observable& orig);
    double getTheoryValue();
    virtual ~Observable();
    void Set(const std::string distr);
    void Set(const std::string distr, const std::string filename, const std::string histoname);
    void Set(const std::string distr, const double ave,
            const double errg, const double errf);
    std::string name, distr, filename, histoname;
    double ave, errg, errf, min, max;
    bool tMCMC;
    friend std::ostream& operator<<(std::ostream& output, const Observable& o);
private:
    void Init(const std::string distr, const std::string filename, 
            const std::string histoname, const double ave, 
            const double errg, const double errf);
    ThObservable * tho;
};

#endif	/* OBSERVABLE_H */

