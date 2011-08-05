/* 
 * File:   StandardModelMatching.h
 * Author: silvest
 *
 * Created on June 9, 2011, 2:16 PM
 */

#ifndef STANDARDMODELMATCHING_H
#define	STANDARDMODELMATCHING_H

#include "StandardModel.h"

#define LEPS 1.e-5 // tolerance in the limit of S(x,y) to S(x)

class StandardModelMatching {
public:
    StandardModelMatching(const StandardModel& SM_i);
    virtual const std::vector<WilsonCoefficient>& CMdbd2(const StandardModel& SM_i);
    virtual const std::vector<WilsonCoefficient>& CMdbs2(const StandardModel& SM_i);
protected:
    std::vector<WilsonCoefficient> vmcdb, vmcds;
private:
    StandardModel SM;
    double S0(double, double) const;
    double S0(double) const;
    double S0p(double x) const;
    double S11(double x) const;
    double S18(double x) const;
    double S1(double x) const;
    WilsonCoefficient mcdbd2, mcdbs2;
};

#endif	/* STANDARDMODELMATCHING_H */

