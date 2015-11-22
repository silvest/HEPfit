/* 
 * Copyright (C) 2014 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef INPUTPARAMETERS_H
#define	INPUTPARAMETERS_H

/**
 * @class InputParameters
 * @ingroup ComputeObservables
 * @brief A class for defining the default values of the mandatory parameters of
 * the model being used on the library mode.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 */
class InputParameters {
public:
    
    /**
     * @brief The default constructor.
     */
    InputParameters()
    {};
    
    /**
     * @brief The default destructor.
     */
    virtual ~InputParameters()
    {};
    
    /**
     * @brief A method that returns the map of the mandatory model parameters given the model name. 
     * @param[in] ModelName the name of the model being used
     * @return the map of the mandatory parameters 
     */
    std::map<std::string, double> getInputParameters(std::string& ModelName)
    {
        if (ModelName.compare("StandardModel") == 0) return(StandardModel());
        else if (ModelName.compare("NPSTU") == 0) return(NPSTU());
        else if (ModelName.compare("NPSTUZbbbarLR") == 0) return (NPSTUZbbbarLR());
        else if (ModelName.compare("NPEpsilons") == 0) return(NPEpsilons());
        else if (ModelName.compare("NPEpsilons_pureNP") == 0) return(NPEpsilons_pureNP());
        else if (ModelName.compare("NPHiggs") == 0) return(NPHiggs());
        else if (ModelName.compare("NPZbbbar") == 0) return(NPZbbbarVA());
        else if (ModelName.compare("NPZbbbarLR") == 0) return(NPZbbbarLR());
        else if (ModelName.compare("NPZbbbarLinearized") == 0) return(NPZbbbarVA());
        else if (ModelName.compare("NPZbbbarLinearizedLR") == 0) return(NPZbbbarLR());
        else if (ModelName.compare("NPEffectiveBS") == 0) return(NPEffectiveBS());
        else if (ModelName.compare("NPEffectiveBS_LFU") == 0) return(NPEffectiveBS());
        else if (ModelName.compare("NPEffectiveBS_QFU") == 0) return(NPEffectiveBS());
        else if (ModelName.compare("NPEffectiveBS_LFU_QFU") == 0) return(NPEffectiveBS());
        else if (ModelName.compare("NPEffectiveGIMR") == 0) return(NPEffectiveGIMR());
        //else if (ModelName.compare("NPEffectiveGIMR_LFU") == 0) return(NPEffectiveGIMR());
        //else if (ModelName.compare("NPEffectiveGIMR_QFU") == 0) return(NPEffectiveGIMR());
        else if (ModelName.compare("NPEffectiveGIMR_LFU_QFU") == 0) return(NPEffectiveGIMR());
        else throw std::runtime_error("\nERROR: Incorrect model name passed to InputParameters():  " + ModelName + "\n");
    };
    
private:
    
    /**
     * @brief A method that generates the map of the mandatory model parameters for StandardModel.
     * @return the map of the mandatory parameters
     *
     * @attention The parameters are initialized to the central values of the 
     * experimental data used in @cite Ciuchini:2013pca. 
     */
    std::map<std::string, double> StandardModel()
    {
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
        DPars_IN["MKstar"] = 0.;
        DPars_IN["tKstar"] = 0.;
        DPars_IN["Mphi"] = 0.;
        DPars_IN["tphi"] = 0.;
        DPars_IN["MBp"] = 0.;
        DPars_IN["tBp"] = 0.;
        DPars_IN["MK0"] = 0.;
        DPars_IN["MKp"] = 0.;
        DPars_IN["FK"] = 0.;
        DPars_IN["FBs"] = 0.;
        DPars_IN["FKstar"] = 0.;
        DPars_IN["FKstarp"] = 0.;
        DPars_IN["Fphi"] = 0.;
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
        DPars_IN["reh_0"] = 0.;
        DPars_IN["reh_p"] = 0.;
        DPars_IN["reh_m"] = 0.;
        DPars_IN["imh_0"] = 0.;
        DPars_IN["imh_p"] = 0.;
        DPars_IN["imh_m"] = 0.;
        DPars_IN["reh_0_1"] = 0.;
        DPars_IN["reh_p_1"] = 0.;
        DPars_IN["reh_m_1"] = 0.;
        DPars_IN["imh_0_1"] = 0.;
        DPars_IN["imh_p_1"] = 0.;
        DPars_IN["imh_m_1"] = 0.;
        DPars_IN["reh_0_2"] = 0.;
        DPars_IN["reh_p_2"] = 0.;
        DPars_IN["reh_m_2"] = 0.;
        DPars_IN["imh_0_2"] = 0.;
        DPars_IN["imh_p_2"] = 0.;
        DPars_IN["imh_m_2"] = 0.;
        DPars_IN["reh_0_MP"] = 0.;
        DPars_IN["imh_0_MP"] = 0.;
        DPars_IN["reh_0_1_MP"] = 0.;
        DPars_IN["imh_0_1_MP"] = 0.;
        DPars_IN["a_0V"] = 0.;
        DPars_IN["a_1V"] = 0.;
        DPars_IN["a_2V"] = 0.;
        DPars_IN["MRV"] = 0.;
        DPars_IN["a_0A0"] = 0.;
        DPars_IN["a_1A0"] = 0.;
        DPars_IN["a_2A0"] = 0.;
        DPars_IN["MRA0"] = 0.;
        DPars_IN["a_0A1"] = 0.;
        DPars_IN["a_1A1"] = 0.;
        DPars_IN["a_2A1"] = 0.;
        DPars_IN["MRA1"] = 0.;
        DPars_IN["a_0A12"] = 0.;
        DPars_IN["a_1A12"] = 0.;
        DPars_IN["a_2A12"] = 0.;
        DPars_IN["MRA12"] = 0.;
        DPars_IN["a_0T1"] = 0.;
        DPars_IN["a_1T1"] = 0.;
        DPars_IN["a_2T1"] = 0.;
        DPars_IN["MRT1"] = 0.;
        DPars_IN["a_0T2"] = 0.;
        DPars_IN["a_1T2"] = 0.;
        DPars_IN["a_2T2"] = 0.;
        DPars_IN["MRT2"] = 0.;
        DPars_IN["a_0T23"] = 0.;
        DPars_IN["a_1T23"] = 0.;
        DPars_IN["a_2T23"] = 0.;
        DPars_IN["MRT23"] = 0.;
        DPars_IN["a_0Vphi"] = 0.;
        DPars_IN["a_1Vphi"] = 0.;
        DPars_IN["a_2Vphi"] = 0.;
        DPars_IN["MRVphi"] = 0.;
        DPars_IN["a_0A0phi"] = 0.;
        DPars_IN["a_1A0phi"] = 0.;
        DPars_IN["a_2A0phi"] = 0.;
        DPars_IN["MRA0phi"] = 0.;
        DPars_IN["a_0A1phi"] = 0.;
        DPars_IN["a_1A1phi"] = 0.;
        DPars_IN["a_2A1phi"] = 0.;
        DPars_IN["MRA1phi"] = 0.;
        DPars_IN["a_0A12phi"] = 0.;
        DPars_IN["a_1A12phi"] = 0.;
        DPars_IN["a_2A12phi"] = 0.;
        DPars_IN["MRA12phi"] = 0.;
        DPars_IN["a_0T1phi"] = 0.;
        DPars_IN["a_1T1phi"] = 0.;
        DPars_IN["a_2T1phi"] = 0.;
        DPars_IN["MRT1phi"] = 0.;
        DPars_IN["a_0T2phi"] = 0.;
        DPars_IN["a_1T2phi"] = 0.;
        DPars_IN["a_2T2phi"] = 0.;
        DPars_IN["MRT2phi"] = 0.;
        DPars_IN["a_0T23phi"] = 0.;
        DPars_IN["a_1T23phi"] = 0.;
        DPars_IN["a_2T23phi"] = 0.;
        DPars_IN["MRT23phi"] = 0.;
        DPars_IN["r_1_fplus"] = 0.;
        DPars_IN["r_2_fplus"] = 0.;
        DPars_IN["m_fit2_fplus"] = 0.;
        DPars_IN["r_1_fT"] = 0.;
        DPars_IN["r_2_fT"] = 0.;
        DPars_IN["m_fit2_fT"] = 0.;
        DPars_IN["r_2_f0"] = 0.;
        DPars_IN["m_fit2_f0"] = 0.;
        DPars_IN["bsgamma_E0"] = 0.;
        DPars_IN["BLNPcorr"] = 0.;
        DPars_IN["Gambino_mukin"] = 0.;
        DPars_IN["Gambino_BRsem"] = 0.;
        DPars_IN["Gambino_Mbkin"] = 0.;
        DPars_IN["Gambino_Mcatmuc"] = 0.;
        DPars_IN["Gambino_mupi2"] = 0.;
        DPars_IN["Gambino_rhoD3"] = 0.;
        DPars_IN["Gambino_muG2"] = 0.;
        DPars_IN["Gambino_rhoLS3"] = 0.;
        DPars_IN["lambdaB"] = 0.3;
        DPars_IN["alpha1kst"] = 0.;
        DPars_IN["alpha2kst"] = 0.;
        DPars_IN["DGs_Gs"] = 0.;
        
        return (DPars_IN);
    };
    
    /**
     * @brief A method that generates the map of the mandatory model parameters for NPEpsilons.
     * @return the map of the mandatory parameters
     *
     * @attention The values of the epsilon parameters are initialized to their 
     * SM values @cite Ciuchini:2013pca.
     */
    std::map<std::string, double> NPEpsilons()
    {
        DPars_IN = StandardModel();
        
        DPars_IN["epsilon_1"] = 0.00521;
        DPars_IN["epsilon_2"] = -0.00737;
        DPars_IN["epsilon_3"] = 0.00528;
        DPars_IN["epsilon_b"] = -0.00694;
        
        return (DPars_IN);
    };
    
    /**
     * @brief A method that generates the map of the mandatory model parameters for NPEpsilons_pureNP.
     * @return the map of the mandatory parameters 
     */
    std::map<std::string, double> NPEpsilons_pureNP()
    {
        DPars_IN = StandardModel();
        
        DPars_IN["delEps_1"] = 0.;
        DPars_IN["delEps_2"] = 0.;
        DPars_IN["delEps_3"] = 0.;
        DPars_IN["delEps_b"] = 0.;
        
        return (DPars_IN);
    };
    
    /**
     * @brief A method that generates the map of the mandatory model parameters for NPSTU.
     * @return the map of the mandatory parameters 
     */
    std::map<std::string, double> NPSTU()
    {
        DPars_IN = StandardModel();
        
        DPars_IN["obliqueS"] = 0.;
        DPars_IN["obliqueT"] = 0.;
        DPars_IN["obliqueU"] = 0.;
        
        return (DPars_IN);
    };
    
    /**
     * @brief A method that generates the map of the mandatory model parameters for NPHiggs.
     * @return the map of the mandatory parameters 
     *
     * @attention The cutoff scale LambdaNP is initialize to 0, i.e., it is taken
     * to be @f$\Lambda = 4\pi v/\sqrt{|1-a^2|}@f$ as explained in the description
     * of @ref NPHiggsParameters "NPHiggs" class.
     */
    std::map<std::string, double> NPHiggs()
    {
        DPars_IN = StandardModel();
                
        DPars_IN["a"] = 1.;
        DPars_IN["b"] = 1.;
        DPars_IN["c_u"] = 1.;
        DPars_IN["c_d"] = 1.;
        DPars_IN["c_e"] = 1.;
        DPars_IN["d_3"] = 1.;
        DPars_IN["d_4"] = 1.;
        DPars_IN["LambdaNP"] = 0.;
        
        return (DPars_IN);
    };

    /**
     * @brief A method that generates the map of the mandatory model parameters for NPSTUZbbbarLR.
     * @return the map of the mandatory parameters
     */ 
    std::map<std::string, double> NPSTUZbbbarLR()
    {
        DPars_IN = NPSTU();

        DPars_IN["deltaGLb"] = 0.;
        DPars_IN["deltaGRb"] = 0.;

        return (DPars_IN);
    };

    /**
     * @brief A method that generates the map of the mandatory model parameters for NPZbbbar and NPZbbbarLinearized.
     * @return the map of the mandatory parameters
     */
    std::map<std::string, double> NPZbbbarVA()
    {
        DPars_IN = StandardModel();
        
        DPars_IN["deltaGVb"] = 0.;
        DPars_IN["deltaGAb"] = 0.;
        
        return (DPars_IN);
    };

    /**
     * @brief A method that generates the map of the mandatory model parameters for NPZbbbar and NPZbbbarLinearized.
     * @return the map of the mandatory parameters
     */
    std::map<std::string, double> NPZbbbarLR()
    {
        DPars_IN = StandardModel();

        DPars_IN["deltaGLb"] = 0.;
        DPars_IN["deltaGRb"] = 0.;

        return (DPars_IN);
    };

    /**
     * @brief A method that generates the map of the mandatory model parameters for NPEffectiveBS.
     * @return the map of the mandatory parameters 
     *
     * @attention The NP scale Lambda_NP is initialize to 1 TeV.
     */
    std::map<std::string, double> NPEffectiveBS()
    {
        DPars_IN = StandardModel();
                        
        DPars_IN["cWB_NP"] = 0.;
        DPars_IN["cH_NP"] = 0.;
        DPars_IN["cLL_NP"] = 0.;
        DPars_IN["cHLp_NP"] = 0.;
        DPars_IN["cHQp_NP"] = 0.;
        DPars_IN["cHL_NP"] = 0.;
        DPars_IN["cHQ_NP"] = 0.;
        DPars_IN["cHE_NP"] = 0.;
        DPars_IN["cHU_NP"] = 0.;
        DPars_IN["cHD_NP"] = 0.;
        DPars_IN["Lambda_NP"] = 1000.;
        
        return (DPars_IN);
    };
    
    /**
     * @brief A method that generates the map of the mandatory model parameters for NPEffectiveGIMR.
     * @return the map of the mandatory parameters
     *
     * @attention The NP scale Lambda_NP is initialize to 1 TeV.
     */
    std::map<std::string, double> NPEffectiveGIMR()
    {
        DPars_IN = StandardModel();
        
        DPars_IN["cWB_NP"] = 0.;
        DPars_IN["cH_NP"] = 0.;
        DPars_IN["cLL_NP"] = 0.;
        DPars_IN["cHLp_NP"] = 0.;
        DPars_IN["cHQ1p_NP"] = 0.;
        DPars_IN["cHQ2p_NP"] = 0.;
        DPars_IN["cHQ3p_NP"] = 0.;
        DPars_IN["cHL_NP"] = 0.;
        DPars_IN["cHQ1_NP"] = 0.;
        DPars_IN["cHQ2_NP"] = 0.;
        DPars_IN["cHQ3_NP"] = 0.;
        DPars_IN["cHE_NP"] = 0.;
        DPars_IN["cHU1_NP"] = 0.;
        DPars_IN["cHU2_NP"] = 0.;
        DPars_IN["cHU3_NP"] = 0.;
        DPars_IN["cHD1_NP"] = 0.;
        DPars_IN["cHD2_NP"] = 0.;
        DPars_IN["cHD3_NP"] = 0.;
        DPars_IN["Lambda_NP"] = 1000.;
        
        return (DPars_IN);
    };
    
    std::map<std::string, double> DPars_IN;///< A map for the list of mandatory parameters in the model being used.
};

#endif	/* INPUTPARAMETERS_H */
    