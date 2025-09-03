/*
 * Copyright (C) 2017 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef THDMWSIGNALSTRENGTHS_H
#define THDMWSIGNALSTRENGTHS_H

#include "NPbase.h"

class THDMW;
class THDMWcache;

class THDMWSignalStrengths : public NPbase {
public:
    /**
     * @brief Constructor.
     */
    THDMWSignalStrengths(const StandardModel& SM_i);

    virtual const double muggH(const double sqrt_s) const;
    virtual const double muVBF(const double sqrt_s) const;
    virtual const double mueeWBF(const double sqrt_s, const double Pol_em, const double Pol_ep) const;
    virtual const double muWH(const double sqrt_s) const;
    virtual const double muZH(const double sqrt_s) const;
    virtual const double mueeZH(const double sqrt_s, const double Pol_em, const double Pol_ep) const;
    virtual const double muVH(const double sqrt_s) const;
    virtual const double muVBFpVH(const double sqrt_s) const;
    virtual const double muttH(const double sqrt_s) const;
//    virtual const double muggHpttH(const double sqrt_s) const;
//    virtual const double mueettH(const double sqrt_s) const;
//    virtual const double Gammagg() const;
//    virtual const double GammaWW() const;
//    virtual const double GammaZZ() const;
//    virtual const double GammaZga() const;
//    virtual const double Gammagaga() const;
//    virtual const double Gammamumu() const;
//    virtual const double Gammatautau() const;
//    virtual const double Gammacc() const;
//    virtual const double Gammabb() const;
//    virtual const double GammaTotal() const;
    virtual const double BrHggRatio() const;
    virtual const double BrHWWRatio() const;
    virtual const double BrHZZRatio() const;
    virtual const double BrHZgaRatio() const;
    virtual const double BrHgagaRatio() const;
    virtual const double BrHmumuRatio() const;
    virtual const double BrHtautauRatio() const;
    virtual const double BrHccRatio() const;
    virtual const double BrHbbRatio() const;
    virtual const double computeGammaTotalRatio() const;

private:
    const THDMW& myTHDMW;
};

#endif /* THDMWSIGNALSTRENGTHS_H */

