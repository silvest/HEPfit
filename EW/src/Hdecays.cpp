/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "Hdecays.h"
#include "StandardModel.h"

double Htobb::computeThValue()
{   
    return SM.GammaHtobb();
}

double Htocc::computeThValue()
{   
    return SM.GammaHtocc();
}

double Htoss::computeThValue()
{   
    return SM.GammaHtoss();
}

double Htotautau::computeThValue()
{   
    return SM.GammaHtotautau();
}

double Htomumu::computeThValue()
{   
    return SM.GammaHtomumu();
}

double HtoWW::computeThValue()
{   
    return SM.GammaHtoWWstar();
}

double HtoZZ::computeThValue()
{   
    return SM.GammaHtoZZstar();
}

double Htogaga::computeThValue()
{   
    return SM.GammaHtogaga();
}

double HtoZga::computeThValue()
{   
    return SM.GammaHtoZga();
}

double Htogg::computeThValue()
{   
    return SM.GammaHtogg();
}

double Hwidth::computeThValue()
{   
    return SM.GammaHTot();
}