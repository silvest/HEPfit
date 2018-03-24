/*
 * Copyright (C) 2014 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "NP_couplings.h"


//-----  Zff couplings observables  ----------

/* -------------------------------------*/

deltagZveveL::deltagZveveL(const StandardModel& SM_i):

        ThObservable(SM_i), 
        myNPbase(static_cast<const NPbase*> (&SM_i))
{
}


deltagZveveL::~deltagZveveL()
{}

double deltagZveveL::computeThValue()
{
    double dgV = myNPbase->deltaGV_f(SM.getLeptons(StandardModel::NEUTRINO_1));
    double dgA = myNPbase->deltaGA_f(SM.getLeptons(StandardModel::NEUTRINO_1));
    double gVSM = SM.gV_f(SM.getLeptons(StandardModel::NEUTRINO_1)).real();
    double gASM = SM.gA_f(SM.getLeptons(StandardModel::NEUTRINO_1)).real();
    
    return (dgV + dgA)/(gVSM + gASM);
}

/* -------------------------------------*/

deltagZvmuvmuL::deltagZvmuvmuL(const StandardModel& SM_i):

        ThObservable(SM_i), 
        myNPbase(static_cast<const NPbase*> (&SM_i))
{
}


deltagZvmuvmuL::~deltagZvmuvmuL()
{}

double deltagZvmuvmuL::computeThValue()
{
    double dgV = myNPbase->deltaGV_f(SM.getLeptons(StandardModel::NEUTRINO_2));
    double dgA = myNPbase->deltaGA_f(SM.getLeptons(StandardModel::NEUTRINO_2));
    double gVSM = SM.gV_f(SM.getLeptons(StandardModel::NEUTRINO_2)).real();
    double gASM = SM.gA_f(SM.getLeptons(StandardModel::NEUTRINO_2)).real();
    
    return (dgV + dgA)/(gVSM + gASM);
}

/* -------------------------------------*/

deltagZvtavtaL::deltagZvtavtaL(const StandardModel& SM_i):

        ThObservable(SM_i), 
        myNPbase(static_cast<const NPbase*> (&SM_i))
{
}


deltagZvtavtaL::~deltagZvtavtaL()
{}

double deltagZvtavtaL::computeThValue()
{
    double dgV = myNPbase->deltaGV_f(SM.getLeptons(StandardModel::NEUTRINO_3));
    double dgA = myNPbase->deltaGA_f(SM.getLeptons(StandardModel::NEUTRINO_3));
    double gVSM = SM.gV_f(SM.getLeptons(StandardModel::NEUTRINO_3)).real();
    double gASM = SM.gA_f(SM.getLeptons(StandardModel::NEUTRINO_3)).real();
    
    return (dgV + dgA)/(gVSM + gASM);
}


/* -------------------------------------*/

deltagZeeL::deltagZeeL(const StandardModel& SM_i):

        ThObservable(SM_i), 
        myNPbase(static_cast<const NPbase*> (&SM_i))
{
}


deltagZeeL::~deltagZeeL()
{}

double deltagZeeL::computeThValue()
{
    double dgV = myNPbase->deltaGV_f(SM.getLeptons(StandardModel::ELECTRON));
    double dgA = myNPbase->deltaGA_f(SM.getLeptons(StandardModel::ELECTRON));
    double gVSM = SM.gV_f(SM.getLeptons(StandardModel::ELECTRON)).real();
    double gASM = SM.gA_f(SM.getLeptons(StandardModel::ELECTRON)).real();
    
    return (dgV + dgA)/(gVSM + gASM);
}

/* -------------------------------------*/

deltagZeeR::deltagZeeR(const StandardModel& SM_i):

        ThObservable(SM_i), 
        myNPbase(static_cast<const NPbase*> (&SM_i))
{
}


deltagZeeR::~deltagZeeR()
{}

double deltagZeeR::computeThValue()
{
    double dgV = myNPbase->deltaGV_f(SM.getLeptons(StandardModel::ELECTRON));
    double dgA = myNPbase->deltaGA_f(SM.getLeptons(StandardModel::ELECTRON));
    double gVSM = SM.gV_f(SM.getLeptons(StandardModel::ELECTRON)).real();
    double gASM = SM.gA_f(SM.getLeptons(StandardModel::ELECTRON)).real();

    return (dgV - dgA)/(gVSM - gASM);
}

/* -------------------------------------*/

deltagZmumuL::deltagZmumuL(const StandardModel& SM_i):

        ThObservable(SM_i), 
        myNPbase(static_cast<const NPbase*> (&SM_i))
{
}


deltagZmumuL::~deltagZmumuL()
{}

double deltagZmumuL::computeThValue()
{
    double dgV = myNPbase->deltaGV_f(SM.getLeptons(StandardModel::MU));
    double dgA = myNPbase->deltaGA_f(SM.getLeptons(StandardModel::MU));
    double gVSM = SM.gV_f(SM.getLeptons(StandardModel::MU)).real();
    double gASM = SM.gA_f(SM.getLeptons(StandardModel::MU)).real();
    
    return (dgV + dgA)/(gVSM + gASM);
}

/* -------------------------------------*/

deltagZmumuR::deltagZmumuR(const StandardModel& SM_i):

        ThObservable(SM_i), 
        myNPbase(static_cast<const NPbase*> (&SM_i))
{
}


deltagZmumuR::~deltagZmumuR()
{}

double deltagZmumuR::computeThValue()
{
    double dgV = myNPbase->deltaGV_f(SM.getLeptons(StandardModel::MU));
    double dgA = myNPbase->deltaGA_f(SM.getLeptons(StandardModel::MU));
    double gVSM = SM.gV_f(SM.getLeptons(StandardModel::MU)).real();
    double gASM = SM.gA_f(SM.getLeptons(StandardModel::MU)).real();

    return (dgV - dgA)/(gVSM - gASM);
}

/* -------------------------------------*/

deltagZtataL::deltagZtataL(const StandardModel& SM_i):

        ThObservable(SM_i), 
        myNPbase(static_cast<const NPbase*> (&SM_i))
{
}


deltagZtataL::~deltagZtataL()
{}

double deltagZtataL::computeThValue()
{
    double dgV = myNPbase->deltaGV_f(SM.getLeptons(StandardModel::TAU));
    double dgA = myNPbase->deltaGA_f(SM.getLeptons(StandardModel::TAU));
    double gVSM = SM.gV_f(SM.getLeptons(StandardModel::TAU)).real();
    double gASM = SM.gA_f(SM.getLeptons(StandardModel::TAU)).real();
    
    return (dgV + dgA)/(gVSM + gASM);
}

/* -------------------------------------*/

deltagZtataR::deltagZtataR(const StandardModel& SM_i):

        ThObservable(SM_i), 
        myNPbase(static_cast<const NPbase*> (&SM_i))
{
}


deltagZtataR::~deltagZtataR()
{}

double deltagZtataR::computeThValue()
{
    double dgV = myNPbase->deltaGV_f(SM.getLeptons(StandardModel::TAU));
    double dgA = myNPbase->deltaGA_f(SM.getLeptons(StandardModel::TAU));
    double gVSM = SM.gV_f(SM.getLeptons(StandardModel::TAU)).real();
    double gASM = SM.gA_f(SM.getLeptons(StandardModel::TAU)).real();

    return (dgV - dgA)/(gVSM - gASM);
}


/* -------------------------------------*/

deltagZuuL::deltagZuuL(const StandardModel& SM_i):

        ThObservable(SM_i), 
        myNPbase(static_cast<const NPbase*> (&SM_i))
{
}


deltagZuuL::~deltagZuuL()
{}

double deltagZuuL::computeThValue()
{
    double dgV = myNPbase->deltaGV_f(SM.getQuarks(StandardModel::UP));
    double dgA = myNPbase->deltaGA_f(SM.getQuarks(StandardModel::UP));
    double gVSM = SM.gV_f(SM.getQuarks(StandardModel::UP)).real();
    double gASM = SM.gA_f(SM.getQuarks(StandardModel::UP)).real();
    
    return (dgV + dgA)/(gVSM + gASM);
}

/* -------------------------------------*/

deltagZuuR::deltagZuuR(const StandardModel& SM_i):

        ThObservable(SM_i), 
        myNPbase(static_cast<const NPbase*> (&SM_i))
{
}


deltagZuuR::~deltagZuuR()
{}

double deltagZuuR::computeThValue()
{
    double dgV = myNPbase->deltaGV_f(SM.getQuarks(StandardModel::UP));
    double dgA = myNPbase->deltaGA_f(SM.getQuarks(StandardModel::UP));
    double gVSM = SM.gV_f(SM.getQuarks(StandardModel::UP)).real();
    double gASM = SM.gA_f(SM.getQuarks(StandardModel::UP)).real();

    return (dgV - dgA)/(gVSM - gASM);
}

/* -------------------------------------*/

deltagZccL::deltagZccL(const StandardModel& SM_i):

        ThObservable(SM_i), 
        myNPbase(static_cast<const NPbase*> (&SM_i))
{
}


deltagZccL::~deltagZccL()
{}

double deltagZccL::computeThValue()
{
    double dgV = myNPbase->deltaGV_f(SM.getQuarks(StandardModel::CHARM));
    double dgA = myNPbase->deltaGA_f(SM.getQuarks(StandardModel::CHARM));
    double gVSM = SM.gV_f(SM.getQuarks(StandardModel::CHARM)).real();
    double gASM = SM.gA_f(SM.getQuarks(StandardModel::CHARM)).real();
    
    return (dgV + dgA)/(gVSM + gASM);
}

/* -------------------------------------*/

deltagZccR::deltagZccR(const StandardModel& SM_i):

        ThObservable(SM_i), 
        myNPbase(static_cast<const NPbase*> (&SM_i))
{
}


deltagZccR::~deltagZccR()
{}

double deltagZccR::computeThValue()
{
    double dgV = myNPbase->deltaGV_f(SM.getQuarks(StandardModel::CHARM));
    double dgA = myNPbase->deltaGA_f(SM.getQuarks(StandardModel::CHARM));
    double gVSM = SM.gV_f(SM.getQuarks(StandardModel::CHARM)).real();
    double gASM = SM.gA_f(SM.getQuarks(StandardModel::CHARM)).real();

    return (dgV - dgA)/(gVSM - gASM);
}


/* -------------------------------------*/

deltagZttL::deltagZttL(const StandardModel& SM_i):

        ThObservable(SM_i), 
        myNPbase(static_cast<const NPbase*> (&SM_i))
{
}


deltagZttL::~deltagZttL()
{}

double deltagZttL::computeThValue()
{
//    Ztt eff. couplings not available in StandardModel class. Compare with tree level 
//    Corrections to Ztt eff. couplings are 0 by default in NPBase, unless overrriden. 
    double dgV = myNPbase->deltaGV_f(SM.getQuarks(StandardModel::TOP));
    double dgA = myNPbase->deltaGA_f(SM.getQuarks(StandardModel::TOP));
    double gSM = (SM.getQuarks(StandardModel::TOP)).getIsospin() 
    - ((SM.getQuarks(StandardModel::TOP)).getCharge())*(SM.sW2());

    return 0.5*(dgV + dgA)/gSM;
}

/* -------------------------------------*/

deltagZttR::deltagZttR(const StandardModel& SM_i):

        ThObservable(SM_i), 
        myNPbase(static_cast<const NPbase*> (&SM_i))
{
}


deltagZttR::~deltagZttR()
{}

double deltagZttR::computeThValue()
{
//    Ztt eff. couplings not available in StandardModel class. Compare with tree level 
    double dgV = myNPbase->deltaGV_f(SM.getQuarks(StandardModel::TOP));
    double dgA = myNPbase->deltaGA_f(SM.getQuarks(StandardModel::TOP));
    double gSM = - ((SM.getQuarks(StandardModel::TOP)).getCharge())*(SM.sW2());

    return 0.5*(dgV - dgA)/gSM;
}

/* -------------------------------------*/

deltagZttV::deltagZttV(const StandardModel& SM_i):

        ThObservable(SM_i), 
        myNPbase(static_cast<const NPbase*> (&SM_i))
{
}


deltagZttV::~deltagZttV()
{}

double deltagZttV::computeThValue()
{
//    Ztt eff. couplings not available in StandardModel class. Compare with tree level 
//    Corrections to Ztt eff. couplings are 0 by default in NPBase, unless overrriden. 
    double dgV = myNPbase->deltaGV_f(SM.getQuarks(StandardModel::TOP));
    double gSM = ((SM.getQuarks(StandardModel::TOP)).getIsospin()) * (1.0 - 4.0*fabs(SM.getQuarks(StandardModel::TOP).getCharge())*(SM.sW2()));

    return dgV/gSM;
}


/* -------------------------------------*/

deltagZttA::deltagZttA(const StandardModel& SM_i):

        ThObservable(SM_i), 
        myNPbase(static_cast<const NPbase*> (&SM_i))
{
}


deltagZttA::~deltagZttA()
{}

double deltagZttA::computeThValue()
{
//    Ztt eff. couplings not available in StandardModel class. Compare with tree level 
//    Corrections to Ztt eff. couplings are 0 by default in NPBase, unless overrriden. 
    double dgA = myNPbase->deltaGA_f(SM.getQuarks(StandardModel::TOP));
    double gSM = (SM.getQuarks(StandardModel::TOP)).getIsospin();

    return dgA/gSM;
}

/* -------------------------------------*/

deltagZddL::deltagZddL(const StandardModel& SM_i):

        ThObservable(SM_i), 
        myNPbase(static_cast<const NPbase*> (&SM_i))
{
}


deltagZddL::~deltagZddL()
{}

double deltagZddL::computeThValue()
{
    double dgV = myNPbase->deltaGV_f(SM.getQuarks(StandardModel::DOWN));
    double dgA = myNPbase->deltaGA_f(SM.getQuarks(StandardModel::DOWN));
    double gVSM = SM.gV_f(SM.getQuarks(StandardModel::DOWN)).real();
    double gASM = SM.gA_f(SM.getQuarks(StandardModel::DOWN)).real();
    
    return (dgV + dgA)/(gVSM + gASM);
}

/* -------------------------------------*/

deltagZddR::deltagZddR(const StandardModel& SM_i):

        ThObservable(SM_i), 
        myNPbase(static_cast<const NPbase*> (&SM_i))
{
}


deltagZddR::~deltagZddR()
{}

double deltagZddR::computeThValue()
{
    double dgV = myNPbase->deltaGV_f(SM.getQuarks(StandardModel::DOWN));
    double dgA = myNPbase->deltaGA_f(SM.getQuarks(StandardModel::DOWN));
    double gVSM = SM.gV_f(SM.getQuarks(StandardModel::DOWN)).real();
    double gASM = SM.gA_f(SM.getQuarks(StandardModel::DOWN)).real();

    return (dgV - dgA)/(gVSM - gASM);
}

/* -------------------------------------*/

deltagZssL::deltagZssL(const StandardModel& SM_i):

        ThObservable(SM_i), 
        myNPbase(static_cast<const NPbase*> (&SM_i))
{
}


deltagZssL::~deltagZssL()
{}

double deltagZssL::computeThValue()
{
    double dgV = myNPbase->deltaGV_f(SM.getQuarks(StandardModel::STRANGE));
    double dgA = myNPbase->deltaGA_f(SM.getQuarks(StandardModel::STRANGE));
    double gVSM = SM.gV_f(SM.getQuarks(StandardModel::STRANGE)).real();
    double gASM = SM.gA_f(SM.getQuarks(StandardModel::STRANGE)).real();
    
    return (dgV + dgA)/(gVSM + gASM);
}

/* -------------------------------------*/

deltagZssR::deltagZssR(const StandardModel& SM_i):

        ThObservable(SM_i), 
        myNPbase(static_cast<const NPbase*> (&SM_i))
{
}


deltagZssR::~deltagZssR()
{}

double deltagZssR::computeThValue()
{
    double dgV = myNPbase->deltaGV_f(SM.getQuarks(StandardModel::STRANGE));
    double dgA = myNPbase->deltaGA_f(SM.getQuarks(StandardModel::STRANGE));
    double gVSM = SM.gV_f(SM.getQuarks(StandardModel::STRANGE)).real();
    double gASM = SM.gA_f(SM.getQuarks(StandardModel::STRANGE)).real();

    return (dgV - dgA)/(gVSM - gASM);
}

/* -------------------------------------*/

deltagZbbL::deltagZbbL(const StandardModel& SM_i):

        ThObservable(SM_i), 
        myNPbase(static_cast<const NPbase*> (&SM_i))
{
}


deltagZbbL::~deltagZbbL()
{}

double deltagZbbL::computeThValue()
{
    double dgV = myNPbase->deltaGV_f(SM.getQuarks(StandardModel::BOTTOM));
    double dgA = myNPbase->deltaGA_f(SM.getQuarks(StandardModel::BOTTOM));
    double gVSM = SM.gV_f(SM.getQuarks(StandardModel::BOTTOM)).real();
    double gASM = SM.gA_f(SM.getQuarks(StandardModel::BOTTOM)).real();
    
    return (dgV + dgA)/(gVSM + gASM);
}

/* -------------------------------------*/

deltagZbbR::deltagZbbR(const StandardModel& SM_i):

        ThObservable(SM_i), 
        myNPbase(static_cast<const NPbase*> (&SM_i))
{
}


deltagZbbR::~deltagZbbR()
{}

double deltagZbbR::computeThValue()
{
    double dgV = myNPbase->deltaGV_f(SM.getQuarks(StandardModel::BOTTOM));
    double dgA = myNPbase->deltaGA_f(SM.getQuarks(StandardModel::BOTTOM));
    double gVSM = SM.gV_f(SM.getQuarks(StandardModel::BOTTOM)).real();
    double gASM = SM.gA_f(SM.getQuarks(StandardModel::BOTTOM)).real();

    return (dgV - dgA)/(gVSM - gASM);
}

/* -------------------------------------*/


//-----  Wff couplings observables  ----------

/* -------------------------------------*/

deltaUWeve::deltaUWeve(const StandardModel& SM_i):

        ThObservable(SM_i), 
        myNPbase(static_cast<const NPbase*> (&SM_i))
{
}


deltaUWeve::~deltaUWeve()
{}

double deltaUWeve::computeThValue()
{
    double dU = myNPbase->deltaGL_Wff(SM.getLeptons(StandardModel::NEUTRINO_1), SM.getLeptons(StandardModel::ELECTRON)).real();
    double gSM = 1.;
    
    return dU/gSM;
}

/* -------------------------------------*/

deltaUWmuvmu::deltaUWmuvmu(const StandardModel& SM_i):

        ThObservable(SM_i), 
        myNPbase(static_cast<const NPbase*> (&SM_i))
{
}


deltaUWmuvmu::~deltaUWmuvmu()
{}

double deltaUWmuvmu::computeThValue()
{
    double dU = myNPbase->deltaGL_Wff(SM.getLeptons(StandardModel::NEUTRINO_2), SM.getLeptons(StandardModel::MU)).real();
    double gSM = 1.;
    
    return dU/gSM;
}

/* -------------------------------------*/

deltaUWtavta::deltaUWtavta(const StandardModel& SM_i):

        ThObservable(SM_i), 
        myNPbase(static_cast<const NPbase*> (&SM_i))
{
}


deltaUWtavta::~deltaUWtavta()
{}

double deltaUWtavta::computeThValue()
{
    double dU = myNPbase->deltaGL_Wff(SM.getLeptons(StandardModel::NEUTRINO_3), SM.getLeptons(StandardModel::TAU)).real();
    double gSM = 1.;
    
    return dU/gSM;
}

/* -------------------------------------*/

deltaVudL::deltaVudL(const StandardModel& SM_i):

        ThObservable(SM_i), 
        myNPbase(static_cast<const NPbase*> (&SM_i))
{
}


deltaVudL::~deltaVudL()
{}

double deltaVudL::computeThValue()
{
    double dV = myNPbase->deltaGL_Wff(SM.getQuarks(StandardModel::UP), SM.getQuarks(StandardModel::DOWN)).real();
    double gSM = 1.;
    
    return dV/gSM;
}

/* -------------------------------------*/

deltaVudR::deltaVudR(const StandardModel& SM_i):

        ThObservable(SM_i), 
        myNPbase(static_cast<const NPbase*> (&SM_i))
{
}


deltaVudR::~deltaVudR()
{}

double deltaVudR::computeThValue()
{
    double dV = myNPbase->deltaGR_Wff(SM.getQuarks(StandardModel::UP), SM.getQuarks(StandardModel::DOWN)).real();
    
    return dV;
}

/* -------------------------------------*/


deltaVcsL::deltaVcsL(const StandardModel& SM_i):

        ThObservable(SM_i), 
        myNPbase(static_cast<const NPbase*> (&SM_i))
{
}


deltaVcsL::~deltaVcsL()
{}

double deltaVcsL::computeThValue()
{
    double dV = myNPbase->deltaGL_Wff(SM.getQuarks(StandardModel::CHARM), SM.getQuarks(StandardModel::STRANGE)).real();
    double gSM = 1.;
    
    return dV/gSM;
}

/* -------------------------------------*/

deltaVcsR::deltaVcsR(const StandardModel& SM_i):

        ThObservable(SM_i), 
        myNPbase(static_cast<const NPbase*> (&SM_i))
{
}


deltaVcsR::~deltaVcsR()
{}

double deltaVcsR::computeThValue()
{
    double dV = myNPbase->deltaGR_Wff(SM.getQuarks(StandardModel::CHARM), SM.getQuarks(StandardModel::STRANGE)).real();
    
    return dV;
}

/* -------------------------------------*/


deltaVtbL::deltaVtbL(const StandardModel& SM_i):

        ThObservable(SM_i), 
        myNPbase(static_cast<const NPbase*> (&SM_i))
{
}


deltaVtbL::~deltaVtbL()
{}

double deltaVtbL::computeThValue()
{
    double dV = myNPbase->deltaGL_Wff(SM.getQuarks(StandardModel::TOP), SM.getQuarks(StandardModel::BOTTOM)).real();
    double gSM = 1.;
    
    return dV/gSM;
}

/* -------------------------------------*/

deltaVtbR::deltaVtbR(const StandardModel& SM_i):

        ThObservable(SM_i), 
        myNPbase(static_cast<const NPbase*> (&SM_i))
{
}


deltaVtbR::~deltaVtbR()
{}

double deltaVtbR::computeThValue()
{
    double dV = myNPbase->deltaGR_Wff(SM.getQuarks(StandardModel::TOP), SM.getQuarks(StandardModel::BOTTOM)).real();
    
    return dV;
}

/* -------------------------------------*/


//-----  Hff couplings observables  ----------

/* -------------------------------------*/

deltagHee::deltagHee(const StandardModel& SM_i):

        ThObservable(SM_i), 
        myNPbase(static_cast<const NPbase*> (&SM_i))
{
}


deltagHee::~deltagHee()
{}

double deltagHee::computeThValue()
{
    double dg = myNPbase->deltaG_hff(SM.getLeptons(StandardModel::ELECTRON)).real();
    double gSM = -(SM.getLeptons(StandardModel::ELECTRON)).getMass() / (SM.v());
    
    return dg/gSM;
}

/* -------------------------------------*/

deltagHmumu::deltagHmumu(const StandardModel& SM_i):

        ThObservable(SM_i), 
        myNPbase(static_cast<const NPbase*> (&SM_i))
{
}


deltagHmumu::~deltagHmumu()
{}

double deltagHmumu::computeThValue()
{
    double dg = myNPbase->deltaG_hff(SM.getLeptons(StandardModel::MU)).real();
    double gSM = -(SM.getLeptons(StandardModel::MU)).getMass() / (SM.v());
    
    return dg/gSM;
}

/* -------------------------------------*/

deltagHtata::deltagHtata(const StandardModel& SM_i):

        ThObservable(SM_i), 
        myNPbase(static_cast<const NPbase*> (&SM_i))
{
}


deltagHtata::~deltagHtata()
{}

double deltagHtata::computeThValue()
{
    double dg = myNPbase->deltaG_hff(SM.getLeptons(StandardModel::TAU)).real();
    double gSM = -(SM.getLeptons(StandardModel::TAU)).getMass() / (SM.v());
    
    return dg/gSM;
}

/* -------------------------------------*/

deltagHuu::deltagHuu(const StandardModel& SM_i):

        ThObservable(SM_i), 
        myNPbase(static_cast<const NPbase*> (&SM_i))
{
}


deltagHuu::~deltagHuu()
{}

double deltagHuu::computeThValue()
{
    double dg = myNPbase->deltaG_hff(SM.getQuarks(StandardModel::UP)).real();
    double gSM = -(SM.getQuarks(StandardModel::UP)).getMass() / (SM.v());
    
    return dg/gSM;
}

/* -------------------------------------*/

deltagHcc::deltagHcc(const StandardModel& SM_i):

        ThObservable(SM_i), 
        myNPbase(static_cast<const NPbase*> (&SM_i))
{
}


deltagHcc::~deltagHcc()
{}

double deltagHcc::computeThValue()
{
    double dg = myNPbase->deltaG_hff(SM.getQuarks(StandardModel::CHARM)).real();
    double gSM = -(SM.getQuarks(StandardModel::CHARM)).getMass() / (SM.v());
    
    return dg/gSM;
}

/* -------------------------------------*/

deltagHtt::deltagHtt(const StandardModel& SM_i):

        ThObservable(SM_i), 
        myNPbase(static_cast<const NPbase*> (&SM_i))
{
}


deltagHtt::~deltagHtt()
{}

double deltagHtt::computeThValue()
{
    double dg = myNPbase->deltaG_hff(SM.getQuarks(StandardModel::TOP)).real();
    double gSM = -(SM.getQuarks(StandardModel::TOP)).getMass() / (SM.v());
    
    return dg/gSM;
}

/* -------------------------------------*/

deltagHdd::deltagHdd(const StandardModel& SM_i):

        ThObservable(SM_i), 
        myNPbase(static_cast<const NPbase*> (&SM_i))
{
}


deltagHdd::~deltagHdd()
{}

double deltagHdd::computeThValue()
{
    double dg = myNPbase->deltaG_hff(SM.getQuarks(StandardModel::DOWN)).real();
    double gSM = -(SM.getQuarks(StandardModel::DOWN)).getMass() / (SM.v());
    
    return dg/gSM;
}

/* -------------------------------------*/

deltagHss::deltagHss(const StandardModel& SM_i):

        ThObservable(SM_i), 
        myNPbase(static_cast<const NPbase*> (&SM_i))
{
}


deltagHss::~deltagHss()
{}

double deltagHss::computeThValue()
{
    double dg = myNPbase->deltaG_hff(SM.getQuarks(StandardModel::STRANGE)).real();
    double gSM = -(SM.getQuarks(StandardModel::STRANGE)).getMass() / (SM.v());
    
    return dg/gSM;
}

/* -------------------------------------*/

deltagHbb::deltagHbb(const StandardModel& SM_i):

        ThObservable(SM_i), 
        myNPbase(static_cast<const NPbase*> (&SM_i))
{
}


deltagHbb::~deltagHbb()
{}

double deltagHbb::computeThValue()
{
    double dg = myNPbase->deltaG_hff(SM.getQuarks(StandardModel::BOTTOM)).real();
    double gSM = -(SM.getQuarks(StandardModel::BOTTOM)).getMass() / (SM.v());
    
    return dg/gSM;
}

/* -------------------------------------*/

//-----  HGG couplings observables  ----------

/* -------------------------------------*/

deltagHGG::deltagHGG(const StandardModel& SM_i):

        ThObservable(SM_i), 
        myNPbase(static_cast<const NPbase*> (&SM_i))
{
}


deltagHGG::~deltagHGG()
{}

double deltagHGG::computeThValue()
{
    double dgRatio = myNPbase->deltaG_hggRatio();
    
    return dgRatio;
}

/* -------------------------------------*/

//-----  HZZ couplings observables  ----------

/* -------------------------------------*/

deltagHZZ::deltagHZZ(const StandardModel& SM_i):

        ThObservable(SM_i), 
        myNPbase(static_cast<const NPbase*> (&SM_i))
{
}


deltagHZZ::~deltagHZZ()
{}

double deltagHZZ::computeThValue()
{
    double dg = myNPbase->deltaG3_hZZ();
    double gSM = (SM.getMz()) * (SM.getMz()) / (SM.v());
    
    return dg/gSM;
}

/* -------------------------------------*/

gHZZ1::gHZZ1(const StandardModel& SM_i):

        ThObservable(SM_i), 
        myNPbase(static_cast<const NPbase*> (&SM_i))
{
}


gHZZ1::~gHZZ1()
{}

double gHZZ1::computeThValue()
{
    double gNP = myNPbase->deltaG1_hZZ();
    
    return gNP;
}

/* -------------------------------------*/

gHZZ2::gHZZ2(const StandardModel& SM_i):

        ThObservable(SM_i), 
        myNPbase(static_cast<const NPbase*> (&SM_i))
{
}


gHZZ2::~gHZZ2()
{}

double gHZZ2::computeThValue()
{
    double gNP = myNPbase->deltaG2_hZZ();
    
    return gNP;
}

/* -------------------------------------*/

//-----  HAA couplings observables  ----------

/* -------------------------------------*/

deltagHAA::deltagHAA(const StandardModel& SM_i):

        ThObservable(SM_i), 
        myNPbase(static_cast<const NPbase*> (&SM_i))
{
}


deltagHAA::~deltagHAA()
{}

double deltagHAA::computeThValue()
{
    double dgRatio = myNPbase->deltaG_hAARatio();
    
    return dgRatio;
}

/* -------------------------------------*/

//-----  HZA couplings observables  ----------

/* -------------------------------------*/

deltagHZA::deltagHZA(const StandardModel& SM_i):

        ThObservable(SM_i), 
        myNPbase(static_cast<const NPbase*> (&SM_i))
{
}


deltagHZA::~deltagHZA()
{}

double deltagHZA::computeThValue()
{
    double dgRatio = myNPbase->deltaG1_hZARatio();
    
    return dgRatio;
}

/* -------------------------------------*/

gHZA2::gHZA2(const StandardModel& SM_i):

        ThObservable(SM_i), 
        myNPbase(static_cast<const NPbase*> (&SM_i))
{
}


gHZA2::~gHZA2()
{}

double gHZA2::computeThValue()
{
    double gNP = myNPbase->deltaG2_hZA();
    
    return gNP;
}

/* -------------------------------------*/

//-----  HWW couplings observables  ----------
/* -------------------------------------*/

deltagHWW::deltagHWW(const StandardModel& SM_i):

        ThObservable(SM_i), 
        myNPbase(static_cast<const NPbase*> (&SM_i))
{
}


deltagHWW::~deltagHWW()
{}

double deltagHWW::computeThValue()
{
    double dg = myNPbase->deltaG3_hWW();
    double gSM = 2.0 * (SM.Mw_tree())* (SM.Mw_tree()) / (SM.v());
    
    return dg/gSM;
}

/* -------------------------------------*/

gHWW1::gHWW1(const StandardModel& SM_i):

        ThObservable(SM_i), 
        myNPbase(static_cast<const NPbase*> (&SM_i))
{
}


gHWW1::~gHWW1()
{}

double gHWW1::computeThValue()
{
    double gNP = myNPbase->deltaG1_hWW();
    
    return gNP;
}

/* -------------------------------------*/

gHWW2::gHWW2(const StandardModel& SM_i):

        ThObservable(SM_i), 
        myNPbase(static_cast<const NPbase*> (&SM_i))
{
}


gHWW2::~gHWW2()
{}

double gHWW2::computeThValue()
{
    double gNP = myNPbase->deltaG2_hWW();
    
    return gNP;
}

/* -------------------------------------*/

//-----  HHH couplings observables  ----------

/* -------------------------------------*/

deltalHHH::deltalHHH(const StandardModel& SM_i):

        ThObservable(SM_i), 
        myNPbase(static_cast<const NPbase*> (&SM_i))
{
}


deltalHHH::~deltalHHH()
{}

double deltalHHH::computeThValue()
{
    double dg = 0.;
    double gSM = 1.;
    
    return dg/gSM;
}

/* -------------------------------------*/

//-----  VVV couplings observables  ----------

// See aTGC in EW

