/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "Lepton_Decays.h"
#include "StandardModel.h"

double GammaMuon::computeThValue()
{    
    return SM.Gamma_muon();
}

double GammaTautoMuon::computeThValue()
{    
    return SM.Gamma_tau_l_nunu(SM.getLeptons(SM.MU));
}

double GammaTautoElectron::computeThValue()
{    
    return SM.Gamma_tau_l_nunu(SM.getLeptons(SM.ELECTRON));
}


// LFU Tests in Tau decays

double gmugeTauLFU::computeThValue()
{    
    return SM.TauLFU_gmuge();
}

double gtaugmuTauLFU::computeThValue()
{    
    return SM.TauLFU_gtaugmu();
}

double gtaugeTauLFU::computeThValue()
{    
    return SM.TauLFU_gtauge();
}
