/* 
 * Copyright (C) 2017 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "THDMWSignalStrengths.h"
#include "THDMW.h"
#include "THDMWcache.h"


THDMWSignalStrengths::THDMWSignalStrengths(const StandardModel& SM_i)
: myTHDMW(static_cast<const THDMW&> (SM_i))
{}

const double THDMWSignalStrengths::muggH(const double sqrt_s) const
{
    return myTHDMW.getMyTHDMWCache()->rh_gg;
}

const double THDMWSignalStrengths::muVBF(const double sqrt_s) const
{
    return myTHDMW.getMyTHDMWCache()->rh_VV;
}

const double THDMWSignalStrengths::mueeWBF(const double sqrt_s) const
{
    return myTHDMW.getMyTHDMWCache()->rh_VV;
}

const double THDMWSignalStrengths::muWH(const double sqrt_s) const
{
    return myTHDMW.getMyTHDMWCache()->rh_VV;
}

const double THDMWSignalStrengths::muZH(const double sqrt_s) const
{
    return myTHDMW.getMyTHDMWCache()->rh_VV;
}

const double THDMWSignalStrengths::mueeZH(const double sqrt_s) const
{
    return myTHDMW.getMyTHDMWCache()->rh_VV;
}

const double THDMWSignalStrengths::muVH(const double sqrt_s) const
{
    return myTHDMW.getMyTHDMWCache()->rh_VV;
}

const double THDMWSignalStrengths::muVBFpVH(const double sqrt_s) const
{
    return myTHDMW.getMyTHDMWCache()->rh_VV;
}

const double THDMWSignalStrengths::muttH(const double sqrt_s) const
{
    return myTHDMW.getMyTHDMWCache()->rh_QuQu;
}

const double THDMWSignalStrengths::BrHggRatio() const
{
    return myTHDMW.getMyTHDMWCache()->rh_gg / computeGammaTotalRatio();
}

const double THDMWSignalStrengths::BrHWWRatio() const
{
    return myTHDMW.getMyTHDMWCache()->rh_VV / computeGammaTotalRatio();
}

const double THDMWSignalStrengths::BrHZZRatio() const
{
    return myTHDMW.getMyTHDMWCache()->rh_VV / computeGammaTotalRatio();
}

const double THDMWSignalStrengths::BrHZgaRatio() const
{
    return myTHDMW.getMyTHDMWCache()->rh_Zga / computeGammaTotalRatio();
}

const double THDMWSignalStrengths::BrHgagaRatio() const
{
    return myTHDMW.getMyTHDMWCache()->rh_gaga / computeGammaTotalRatio();
}

const double THDMWSignalStrengths::BrHmumuRatio() const
{
    return myTHDMW.getMyTHDMWCache()->rh_ll / computeGammaTotalRatio();
}

const double THDMWSignalStrengths::BrHtautauRatio() const
{
    return myTHDMW.getMyTHDMWCache()->rh_ll / computeGammaTotalRatio();
}

const double THDMWSignalStrengths::BrHccRatio() const
{
    return myTHDMW.getMyTHDMWCache()->rh_QuQu / computeGammaTotalRatio();
}

const double THDMWSignalStrengths::BrHbbRatio() const
{
    return myTHDMW.getMyTHDMWCache()->rh_QdQd / computeGammaTotalRatio();
}

const double THDMWSignalStrengths::computeGammaTotalRatio() const
{
    return myTHDMW.getMyTHDMWCache()->sumModBRs;
}
