/*
 * Copyright (C) 2014 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "OtherThObservables.h"
#include "NPbase.h"


//-----  Collider observables: LHC dilepton events  ----------

/* -------------------------------------*/

NevLHCee13::NevLHCee13(const StandardModel& SM_i, const int i_bin_i)
: ThObservable(SM_i), i_bin(i_bin_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("NevLHCee13 called with a class whose parent is not NPbase");
}


NevLHCee13::~NevLHCee13()
{}

double NevLHCee13::computeThValue()
{    
    return (myNPbase->NevLHCppee13(i_bin));
}

/* -------------------------------------*/

NevLHCmumu13::NevLHCmumu13(const StandardModel& SM_i, const int i_bin_i)
: ThObservable(SM_i), i_bin(i_bin_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("NevLHCmumu13 called with a class whose parent is not NPbase");
}


NevLHCmumu13::~NevLHCmumu13()
{}

double NevLHCmumu13::computeThValue()
{    
    return (myNPbase->NevLHCppmumu13(i_bin));
}

/* -------------------------------------*/

NevLHCtautau13::NevLHCtautau13(const StandardModel& SM_i, const int i_bin_i)
: ThObservable(SM_i), i_bin(i_bin_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("NevLHCtautau13 called with a class whose parent is not NPbase");
}


NevLHCtautau13::~NevLHCtautau13()
{}

double NevLHCtautau13::computeThValue()
{    
    return (myNPbase->NevLHCpptautau13(i_bin));
}

/* -------------------------------------*/


//-----  Collider observables: LHC mono-lepton events  ----------

NevLHCenu13::NevLHCenu13(const StandardModel& SM_i, const int i_bin_i)
: ThObservable(SM_i), i_bin(i_bin_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("NevLHCenu13 called with a class whose parent is not NPbase");
}


NevLHCenu13::~NevLHCenu13()
{}

double NevLHCenu13::computeThValue()
{    
    return (myNPbase->NevLHCppenu13(i_bin));
}

/* -------------------------------------*/

NevLHCmunu13::NevLHCmunu13(const StandardModel& SM_i, const int i_bin_i)
: ThObservable(SM_i), i_bin(i_bin_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("NevLHCmunu13 called with a class whose parent is not NPbase");
}


NevLHCmunu13::~NevLHCmunu13()
{}

double NevLHCmunu13::computeThValue()
{    
    return (myNPbase->NevLHCppmunu13(i_bin));
}

/* -------------------------------------*/

NevLHCtaunu13::NevLHCtaunu13(const StandardModel& SM_i, const int i_bin_i)
: ThObservable(SM_i), i_bin(i_bin_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("NevLHCtaunu13 called with a class whose parent is not NPbase");
}


NevLHCtaunu13::~NevLHCtaunu13()
{}

double NevLHCtaunu13::computeThValue()
{    
    return (myNPbase->NevLHCpptaunu13(i_bin));
}

/* -------------------------------------*/