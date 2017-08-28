/* 
 * Copyright (C) 2017 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "THDMWSignalStrengths.h"

THDMWSignalStrengths::THDMWSignalStrengths(const StandardModel& SM_i)
: myTHDMW(static_cast<const THDMW&> (SM_i))
{}

double THDMWSignalStrengths::muggH(const double sqrt_s) const
{
    return myTHDMW.getMyTHDMWCache()->rh_gg;
}

double THDMWSignalStrengths::muVBF(const double sqrt_s) const
{
    return myTHDMW.getMyTHDMWCache()->rh_VV;
}

double THDMWSignalStrengths::mueeWBF(const double sqrt_s) const
{
    return myTHDMW.getMyTHDMWCache()->rh_VV;
}

double THDMWSignalStrengths::muWH(const double sqrt_s) const
{
    return myTHDMW.getMyTHDMWCache()->rh_VV;
}

double THDMWSignalStrengths::muZH(const double sqrt_s) const
{
    return myTHDMW.getMyTHDMWCache()->rh_VV;
}

double THDMWSignalStrengths::mueeZH(const double sqrt_s) const
{
    return myTHDMW.getMyTHDMWCache()->rh_VV;
}

double THDMWSignalStrengths::muVH(const double sqrt_s) const
{
    return myTHDMW.getMyTHDMWCache()->rh_VV;
}

double THDMWSignalStrengths::muVBFpVH(const double sqrt_s) const
{
    return myTHDMW.getMyTHDMWCache()->rh_VV;
}

double THDMWSignalStrengths::muttH(const double sqrt_s) const
{
    return myTHDMW.getMyTHDMWCache()->rh_QuQu;
}

double THDMWSignalStrengths::BrHggRatio() const
{
    return myTHDMW.getMyTHDMWCache()->rh_gg / computeGammaTotalRatio();
}

double THDMWSignalStrengths::BrHWWRatio() const
{
    return myTHDMW.getMyTHDMWCache()->rh_VV / computeGammaTotalRatio();
}

double THDMWSignalStrengths::BrHZZRatio() const
{
    return myTHDMW.getMyTHDMWCache()->rh_VV / computeGammaTotalRatio();
}

double THDMWSignalStrengths::BrHZgaRatio() const
{
    return myTHDMW.getMyTHDMWCache()->rh_Zga / computeGammaTotalRatio();
}

double THDMWSignalStrengths::BrHgagaRatio() const
{
    return myTHDMW.getMyTHDMWCache()->rh_gaga / computeGammaTotalRatio();
}

double THDMWSignalStrengths::BrHmumuRatio() const
{
    return myTHDMW.getMyTHDMWCache()->rh_ll / computeGammaTotalRatio();
}

double THDMWSignalStrengths::BrHtautauRatio() const
{
    return myTHDMW.getMyTHDMWCache()->rh_ll / computeGammaTotalRatio();
}

double THDMWSignalStrengths::BrHccRatio() const
{
    return myTHDMW.getMyTHDMWCache()->rh_QuQu / computeGammaTotalRatio();
}

double THDMWSignalStrengths::BrHbbRatio() const
{
    return myTHDMW.getMyTHDMWCache()->rh_QdQd / computeGammaTotalRatio();
}

double THDMWSignalStrengths::computeGammaTotalRatio() const
{
    return myTHDMW.getMyTHDMWCache()->sumModBRs;
}
