/*
 * Copyright (C) 2017 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef THDMWSIGNALSTRENGTHS_H
#define THDMWSIGNALSTRENGTHS_H

#include "NPbase.h"
#include "THDMW.h"
#include "THDMWcache.h"

class THDMWSignalStrengths : public NPbase {
public:
    /**
     * @brief Constructor.
     */
    THDMWSignalStrengths(const StandardModel& SM_i);

    virtual double muggH(const double sqrt_s) const;
    virtual double muVBF(const double sqrt_s) const;
    virtual double mueeWBF(const double sqrt_s) const;
    virtual double muWH(const double sqrt_s) const;
    virtual double muZH(const double sqrt_s) const;
    virtual double mueeZH(const double sqrt_s) const;
    virtual double muVH(const double sqrt_s) const;
    virtual double muVBFpVH(const double sqrt_s) const;
    virtual double muttH(const double sqrt_s) const;
//    virtual double muggHpttH(const double sqrt_s) const;
//    virtual double mueettH(const double sqrt_s) const;
//    virtual double Gammagg() const;
//    virtual double GammaWW() const;
//    virtual double GammaZZ() const;
//    virtual double GammaZga() const;
//    virtual double Gammagaga() const;
//    virtual double Gammamumu() const;
//    virtual double Gammatautau() const;
//    virtual double Gammacc() const;
//    virtual double Gammabb() const;
//    virtual double GammaTotal() const;
    virtual double BrHggRatio() const;
    virtual double BrHWWRatio() const;
    virtual double BrHZZRatio() const;
    virtual double BrHZgaRatio() const;
    virtual double BrHgagaRatio() const;
    virtual double BrHmumuRatio() const;
    virtual double BrHtautauRatio() const;
    virtual double BrHccRatio() const;
    virtual double BrHbbRatio() const;
    virtual double computeGammaTotalRatio() const;

private:
    const THDMW& myTHDMW;
};

#endif /* THDMWSIGNALSTRENGTHS_H */

