/*
 * Copyright (C) 2015 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "GeneralTHDMSTU.h"
#include "GeneralTHDM.h"
#include "GeneralTHDMcache.h"

GeneralTHDMSTU::GeneralTHDMSTU(const StandardModel& SM_i)
: ThObservable(SM_i), myGTHDM(static_cast<const GeneralTHDM*> (&SM_i))

{
// No needed anymore
//    mycache = new GeneralTHDMcache(SM_i);
};

double GeneralTHDMSTU::computeThValue()
{
    return 0.0;
}
/////////////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////////////////
GTHDMDeltaS::GTHDMDeltaS(const StandardModel& SM_i)
: GeneralTHDMSTU(SM_i)
{}

double GTHDMDeltaS::computeThValue()
{
    return myGTHDM->GTHDMDeltaS();
}

GTHDMDeltaT::GTHDMDeltaT(const StandardModel& SM_i)
: GeneralTHDMSTU(SM_i)
{}

double GTHDMDeltaT::computeThValue()
{
    return myGTHDM->GTHDMDeltaT();
}



GTHDMDeltaU::GTHDMDeltaU(const StandardModel& SM_i)
: GeneralTHDMSTU(SM_i)
{}

double GTHDMDeltaU::computeThValue()
{
    return myGTHDM->GTHDMDeltaU();
}