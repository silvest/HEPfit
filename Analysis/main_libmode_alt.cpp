/* 
 * Copyright (C) 2014 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include <iostream>
#include <ComputeObservables.h>

int main(int argc, char** argv) 
{
    try {
        std::map<std::string, double> DObs; /* Map of observables */
        std::map<std::string, double> DPars_IN; /* Map of initialization parameters */
        std::map<std::string, double> DPars; /* Map of parameters being varied. */
        std::map<std::string, std::string> DFlags; /* Map of flags */
        
        /* Initialize  the parameters of the model */
        DPars_IN["GF"] = 1.1663787e-5;
        DPars_IN["ale"] = 7.2973525698e-3;
        DPars_IN["AlsMz"] = 0.1184;
        DPars_IN["dAle5Mz"] = 0.02750;
        DPars_IN["Mz"] = 91.1875;
        DPars_IN["mtop"] = 173.2;
        DPars_IN["mHl"] = 125.6;
        DPars_IN["delMw"] = 0.;
        DPars_IN["delSin2th_l"] = 0.;
        DPars_IN["delGammaZ"] = 0.;
        DPars_IN["mup"] = 0.0023;
        DPars_IN["mdown"] = 0.0048;
        DPars_IN["mstrange"] = 0.095;
        DPars_IN["mcharm"] = 1.275;
        DPars_IN["mbottom"] = 4.18;
        DPars_IN["muc"] = 1.275;
        DPars_IN["mub"] = 4.18;
        DPars_IN["mut"] = 164.;
        DPars_IN["mneutrino_1"] = 0.;
        DPars_IN["mneutrino_2"] = 0.;
        DPars_IN["mneutrino_3"] = 0.;
        DPars_IN["melectron"] = 5.109989e-4;
        DPars_IN["mmu"] = 0.10565837;
        DPars_IN["mtau"] = 1.77682;
        DPars_IN["MBd"] = 0.;
        DPars_IN["tBd"] = 0.;
        DPars_IN["MBs"] = 0.;
        DPars_IN["tBs"] = 0.;
        DPars_IN["MBp"] = 0.;
        DPars_IN["MK0"] = 0.;
        DPars_IN["MKp"] = 0.;
        DPars_IN["FK"] = 0.;
        DPars_IN["FBs"] = 0.;
        DPars_IN["FBsoFBd"] = 0.;
        DPars_IN["BBsoBBd"] = 0.;
        DPars_IN["BBs1"] = 0.;
        DPars_IN["BBs2"] = 0.;
        DPars_IN["BBs3"] = 0.;
        DPars_IN["BBs4"] = 0.;
        DPars_IN["BBs5"] = 0.;
        DPars_IN["BBsscale"] = 0.;
        DPars_IN["BBsscheme"] = 0.;
        DPars_IN["lambda"] = 0.;
        DPars_IN["A"] = 0.;
        DPars_IN["rhob"] = 0.;
        DPars_IN["etab"] = 0.;
        DPars_IN["muw"] = 0.;
        DPars_IN["phiEpsK"] = 0.;
        DPars_IN["KbarEpsK"] = 0.;
        DPars_IN["DeltaMK"] = 0.;
        DPars_IN["Dmk"] = 0.;
        DPars_IN["SM_M12D"] = 0.;
        DPars_IN["MD"] = 0.;
        DPars_IN["FD"] = 0.;
        DPars_IN["BD1"] = 0.;
        DPars_IN["BD2"] = 0.;
        DPars_IN["BD3"] = 0.;
        DPars_IN["BD4"] = 0.;
        DPars_IN["BD5"] = 0.;
        DPars_IN["BDscale"] = 0.;
        DPars_IN["BDscheme"] = 0.;
        DPars_IN["BK1"] = 0.;
        DPars_IN["BK2"] = 0.;
        DPars_IN["BK3"] = 0.;
        DPars_IN["BK4"] = 0.;
        DPars_IN["BK5"] = 0.;
        DPars_IN["BKscale"] = 0.;
        DPars_IN["BKscheme"] = 0.;
        DPars_IN["EpsK"] = 0.;
        DPars_IN["BK(1/2)1"] = 0.;
        DPars_IN["BK(1/2)2"] = 0.;
        DPars_IN["BK(1/2)3"] = 0.;
        DPars_IN["BK(1/2)4"] = 0.;
        DPars_IN["BK(1/2)5"] = 0.;
        DPars_IN["BK(1/2)6"] = 0.;
        DPars_IN["BK(1/2)7"] = 0.;
        DPars_IN["BK(1/2)8"] = 0.;
        DPars_IN["BK(1/2)9"] = 0.;
        DPars_IN["BK(1/2)10"] = 0.;
        DPars_IN["BKd_scale"] = 0.;
        DPars_IN["BKd_scheme"] = 0.;
        DPars_IN["BK(3/2)1"] = 0.;
        DPars_IN["BK(3/2)2"] = 0.;
        DPars_IN["BK(3/2)3"] = 0.;
        DPars_IN["BK(3/2)4"] = 0.;
        DPars_IN["BK(3/2)5"] = 0.;
        DPars_IN["BK(3/2)6"] = 0.;
        DPars_IN["BK(3/2)7"] = 0.;
        DPars_IN["BK(3/2)8"] = 0.;
        DPars_IN["BK(3/2)9"] = 0.;
        DPars_IN["BK(3/2)10"] = 0.;
        DPars_IN["ReA2_Kd"] = 0.;
        DPars_IN["ReA0_Kd"] = 0.;
        DPars_IN["Omega_eta_etap"] = 0.;
        DPars_IN["Br_Kp_P0enu"] = 0.;
        DPars_IN["Br_Kp_munu"] = 0.;
        DPars_IN["Br_B_Xcenu"] = 0.;
        DPars_IN["DeltaP_cu"] = 0.;
        DPars_IN["IB_Kl"] = 0.;
        DPars_IN["IB_Kp"] = 0.;
        DPars_IN["tKl"] = 0.;
        DPars_IN["tKp"] = 0.;
        
        /* Initialize the Observables to be returned */
        DObs["Mw"] = 0;
        DObs["GammaW"] = 0.;
        DObs["GammaZ"] = 0.;
        DObs["Mz"] = 0.;
        ComputeObservables CO("StandardModel", DPars_IN, DObs);
        CO.setFlags(DFlags);
        
        for (int i = 0; i < 1000; i++) {
        /* Pass values to the parameters being varied */
            DPars["Mz"] = 91.1875 + 0.0001 * i;
            DPars["AlsMz"] = 0.1184 + 0.000001 * i;
            
            DObs = CO.compute(DPars);
        
            std::cout << "\nParameters[" << i + 1 << "]:"<< std::endl;
            for (std::map<std::string, double>::iterator it = DPars.begin(); it != DPars.end(); it++){
                std::cout << it->first << " = " << it->second << std::endl;
            }
            std::cout << "\nObservables:[" << i + 1 << "]:" << std::endl;
            for (std::map<std::string, double>::iterator it = DObs.begin(); it != DObs.end(); it++){
                std::cout << it->first << " = " << it->second << std::endl;
            }
        }
        
        return EXIT_SUCCESS;
    } catch (const std::runtime_error& e) {
        std::cerr << e.what() << std::endl;
        return EXIT_FAILURE;
    }
}
