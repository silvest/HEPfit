/* 
 * Copyright (C) 2014 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef INPUTPARAMETERS_H
#define	INPUTPARAMETERS_H

#include <map>

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
        else if (ModelName.compare("NPZbbbar") == 0) return(NPZbbbarVA());
        else if (ModelName.compare("NPZbbbarLR") == 0) return(NPZbbbarLR());
        else if (ModelName.compare("NPZbbbarLinearized") == 0) return(NPZbbbarVA());
        else if (ModelName.compare("NPZbbbarLinearizedLR") == 0) return(NPZbbbarLR());
        else if (ModelName.compare("NPSMEFT6dtopquark") == 0) return(NPSMEFT6dtopquark());
        else if (ModelName.compare("HiggsKvKf") == 0) return(NPHiggsKvKf());
        else if (ModelName.compare("HiggsKvKfgen") == 0) return(NPHiggsKvKfgen());
        else if (ModelName.compare("HiggsKvgenKfgen") == 0) return(NPHiggsKvgenKfgen());
        else if (ModelName.compare("HiggsKigen") == 0) return(NPHiggsKigen());
        else if (ModelName.compare("HiggsChiral") == 0) return(NPHiggsChiral());
        else if (ModelName.compare("NPDF2") == 0) return(NPDF2());
        else if (ModelName.compare("CMFV") == 0) return(NPCMFV());
        else if (ModelName.compare("FlavourWilsonCoefficient") == 0) return(NPFlavourWilsonCoefficient());
        else if (ModelName.compare("FlavourWilsonCoefficient_DF2") == 0) return(NPFlavourWilsonCoefficient_DF2());
        else if (ModelName.compare("RealWeakEFTLFV") == 0) return(NPRealWeakEFTLFV());
        else if (ModelName.compare("RealWeakEFTCC") == 0) return(NPRealWeakEFTCC());
        else if (ModelName.compare("RealWeakEFTCCPM") == 0) return(NPRealWeakEFTCCPM());
        else if (ModelName.compare("THDM") == 0) return(NPTHDM());
        else if (ModelName.compare("GeorgiMachacek") == 0) return(NPGeorgiMachacek());
        else if (ModelName.compare("GeneralTHDM") == 0) return(NPGeneralTHDM());
        else if (ModelName.compare("THDMW") == 0) return(NPTHDMW());
        else if (ModelName.compare("GeneralSUSY") == 0) return(NPGeneralSUSY());
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
        DPars_IN["delR0b"] = 0.0;
        DPars_IN["delSin2th_q"] = 0.;
        DPars_IN["delSin2th_b"] = 0.;
        DPars_IN["delsigma0H"] = 0.;
        DPars_IN["delR0l"] = 0.;
        DPars_IN["delR0c"] = 0.;
        DPars_IN["lambda"] = 0.22506;
        DPars_IN["A"] = 0.811;
        DPars_IN["rhob"] = 0.124;
        DPars_IN["etab"] = 0.356;
        DPars_IN["muw"] = 0.;
        
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
     * @brief A method that generates the map of the mandatory model parameters for NPSMEFT6dtopquark
     * @return the map of the mandatory parameters
     */
    std::map<std::string, double> NPSMEFT6dtopquark()
    {
        DPars_IN = StandardModel();

        DPars_IN["C_phit"] = 0.;
        DPars_IN["C_phiQ3"] = 0.;
        DPars_IN["C_phiQ1"] = 0.;
        DPars_IN["C_tW"] = 0.;
        DPars_IN["C_tB"] = 0.;
        DPars_IN["C_tphi"] = 0.;
        DPars_IN["C_phib"] = 0.;
        DPars_IN["C_bW"] = 0.;
        DPars_IN["C_bB"] = 0.;
        DPars_IN["C_phitb"] = 0.;
        DPars_IN["C_ed"] = 0.;
        DPars_IN["C_eq"] = 0.;
        DPars_IN["C_ld"] = 0.;
        DPars_IN["C_lqP"] = 0.;
        DPars_IN["C_eu"] = 0.;
        DPars_IN["C_lu"] = 0.;
        DPars_IN["C_lqM"] = 0.;

        return (DPars_IN);
    };
   
    /**
     * @brief A method that generates the map of the mandatory model parameters for HiggsKvKf
     * @return the map of the mandatory parameters
     */
    std::map<std::string, double> NPHiggsKvKf()
    {
        DPars_IN = StandardModel();

        DPars_IN["Kv"] = 0.;
        DPars_IN["Kf"] = 0.;
        DPars_IN["BrHinv"] = 0.;

        return (DPars_IN);
    };
    
    /**
     * @brief A method that generates the map of the mandatory model parameters for HiggsKvKf
     * @return the map of the mandatory parameters
     */
    std::map<std::string, double> NPHiggsKvKfgen()
    {
        DPars_IN = StandardModel();

        DPars_IN["Kv"] = 0.;
        DPars_IN["Ku"] = 0.;
        DPars_IN["Kd"] = 0.;
        DPars_IN["Kl"] = 0.;
        DPars_IN["BrHinv"] = 0.;

        return (DPars_IN);
    };
    
    /**
     * @brief A method that generates the map of the mandatory model parameters for HiggsKvgenKfgen
     * @return the map of the mandatory parameters
     */
    std::map<std::string, double> NPHiggsKvgenKfgen()
    {
        DPars_IN = StandardModel();

        DPars_IN["KW"] = 0.;
        DPars_IN["KZ"] = 0.;
        DPars_IN["Ku"] = 0.;
        DPars_IN["Kd"] = 0.;
        DPars_IN["Kl"] = 0.;
        DPars_IN["BrHinv"] = 0.;

        return (DPars_IN);
    };
    
    /**
     * @brief A method that generates the map of the mandatory model parameters for HiggsKigen
     * @return the map of the mandatory parameters
     */
    std::map<std::string, double> NPHiggsKigen()
    {
        DPars_IN = StandardModel();
        
        DPars_IN["Kw"] = 0.;
        DPars_IN["Kz"] = 0.;
        DPars_IN["Kg"] = 0.;
        DPars_IN["Kga"] = 0.;
        DPars_IN["Kzga"] = 0.;
        DPars_IN["Ku"] = 0.;
        DPars_IN["Kc"] = 0.;
        DPars_IN["Kt"] = 0.;
        DPars_IN["Kd"] = 0.;
        DPars_IN["Ks"] = 0.;
        DPars_IN["Kb"] = 0.;
        DPars_IN["Ke"] = 0.;
        DPars_IN["Kmu"] = 0.;
        DPars_IN["Ktau"] = 0.;
        DPars_IN["KH"] = 0.;
        DPars_IN["BrHinv"] = 0.;
        DPars_IN["BrHexo"] = 0.;
        DPars_IN["eggFint"] = 0.;
        DPars_IN["eggFpar"] = 0.;
        DPars_IN["ettHint"] = 0.;
        DPars_IN["ettHpar"] = 0.;
        DPars_IN["eVBFint"] = 0.;
        DPars_IN["eVBFpar"] = 0.;
        DPars_IN["eWHint"] = 0.;
        DPars_IN["eWHpar"] = 0.;
        DPars_IN["eZHint"] = 0.;
        DPars_IN["eZHpar"] = 0.;
        DPars_IN["eeeWBFint"] = 0.;
        DPars_IN["eeeWBFpar"] = 0.;
        DPars_IN["eeeZHint"] = 0.;
        DPars_IN["eeeZHpar"] = 0.;
        DPars_IN["eeettHint"] = 0.;
        DPars_IN["eeettHpar"] = 0.;
        DPars_IN["eepWBFint"] = 0.;
        DPars_IN["eepWBFpar"] = 0.;
        DPars_IN["eepZBFint"] = 0.;
        DPars_IN["eepZBFpar"] = 0.;
        DPars_IN["eHggint"] = 0.;
        DPars_IN["eHggpar"] = 0.;
        DPars_IN["eHWWint"] = 0.;
        DPars_IN["eHWWpar"] = 0.;
        DPars_IN["eHZZint"] = 0.;
        DPars_IN["eHZZpar"] = 0.;
        DPars_IN["eHZgaint"] = 0.;
        DPars_IN["eHZgapar"] = 0.;
        DPars_IN["eHgagaint"] = 0.;
        DPars_IN["eHgagapar"] = 0.;
        DPars_IN["eHmumuint"] = 0.;
        DPars_IN["eHmumupar"] = 0.;
        DPars_IN["eHtautauint"] = 0.;
        DPars_IN["eHtautaupar"] = 0.;
        DPars_IN["eHccint"] = 0.;
        DPars_IN["eHccpar"] = 0.;
        DPars_IN["eHbbint"] = 0.;
        DPars_IN["eHbbpar"] = 0.;
        DPars_IN["eggFHgaga"] = 0.;
        DPars_IN["eggFHZga"] = 0.;
        DPars_IN["eggFHZZ"] = 0.;
        DPars_IN["eggFHWW"] = 0.;
        DPars_IN["eggFHtautau"] = 0.;
        DPars_IN["eggFHbb"] = 0.;
        DPars_IN["eggFHmumu"] = 0.;
        DPars_IN["eVBFHgaga"] = 0.;
        DPars_IN["eVBFHZga"] = 0.;
        DPars_IN["eVBFHZZ"] = 0.;
        DPars_IN["eVBFHWW"] = 0.;
        DPars_IN["eVBFHtautau"] = 0.;
        DPars_IN["eVBFHbb"] = 0.;
        DPars_IN["eVBFHmumu"] = 0.;
        DPars_IN["eWHgaga"] = 0.;
        DPars_IN["eWHZga"] = 0.;
        DPars_IN["eWHZZ"] = 0.;
        DPars_IN["eWHWW"] = 0.;
        DPars_IN["eWHtautau"] = 0.;
        DPars_IN["eWHbb"] = 0.;
        DPars_IN["eWHmumu"] = 0.;
        DPars_IN["eZHgaga"] = 0.;
        DPars_IN["eZHZga"] = 0.;
        DPars_IN["eZHZZ"] = 0.;
        DPars_IN["eZHWW"] = 0.;
        DPars_IN["eZHtautau"] = 0.;
        DPars_IN["eZHbb"] = 0.;
        DPars_IN["eZHmumu"] = 0.;
        DPars_IN["ettHgaga"] = 0.;
        DPars_IN["ettHZga"] = 0.;
        DPars_IN["ettHZZ"] = 0.;
        DPars_IN["ettHWW"] = 0.;
        DPars_IN["ettHtautau"] = 0.;
        DPars_IN["ettHbb"] = 0.;
        DPars_IN["ettHmumu"] = 0.;
        DPars_IN["eVBFHinv"] = 0.;
        DPars_IN["eVHinv"] = 0.;

        return (DPars_IN);
    };
    
    /**
     * @brief A method that generates the map of the mandatory model parameters for HiggsChiral
     * @return the map of the mandatory parameters
     */
    std::map<std::string, double> NPHiggsChiral()
    {
        DPars_IN = StandardModel();

        DPars_IN["cv"] = 0.;
        DPars_IN["ct"] = 0.;
        DPars_IN["cb"] = 0.;
        DPars_IN["cc"] = 0.;
        DPars_IN["ctau"] = 0.;
        DPars_IN["cmu"] = 0.;
        DPars_IN["cg"] = 0.;
        DPars_IN["cga"] = 0.;
        DPars_IN["cZga"] = 0.;
        DPars_IN["obsZgaLimitATLAS13"] = 0.;
        DPars_IN["obsZgaLimitCMS13"] = 0.;
        DPars_IN["obsZgaLimitATLAS"] = 0.;
        DPars_IN["obsZgaLimitCMS"] = 0.;
        DPars_IN["expZgaLimitATLAS13"] = 0.;
        DPars_IN["expZgaLimitCMS13"] = 0.;
        DPars_IN["expZgaLimitATLAS"] = 0.;
        DPars_IN["expZgaLimitCMS"] = 0.;
        
        return (DPars_IN);
    };
    
    /**
     * @brief A method that generates the map of the mandatory model parameters for NPDF2
     * @return the map of the mandatory parameters
     */
    std::map<std::string, double> NPDF2()
    {
        DPars_IN = StandardModel();

        DPars_IN["CepsK"] = 0.;
        DPars_IN["CDMK"] = 0.;
        DPars_IN["CBd"] = 0.;
        DPars_IN["PhiBd"] = 0.;
        DPars_IN["CBs"] = 0.;
        DPars_IN["PhiBs"] = 0.;
        
        return (DPars_IN);
    };
    
    /**
     * @brief A method that generates the map of the mandatory model parameters for CMFV
     * @return the map of the mandatory parameters
     */
    std::map<std::string, double> NPCMFV()
    {
        DPars_IN = StandardModel();

        DPars_IN["Ftt"] = 0.;
        
        return (DPars_IN);
    };
    
    /**
     * @brief A method that generates the map of the mandatory model 
     * parameters for FlavourWilsonCoefficient
     * @return the map of the mandatory parameters
     */
    std::map<std::string, double> NPFlavourWilsonCoefficient()
    {
        DPars_IN = StandardModel();

        DPars_IN["reDC1"] = 0.;
        DPars_IN["reDC2"] = 0.;
        DPars_IN["reDC3"] = 0.;
        DPars_IN["reDC4"] = 0.;
        DPars_IN["reDC5"] = 0.;
        DPars_IN["reDC6"] = 0.;
        DPars_IN["reDC7"] = 0.;
        DPars_IN["reDC8"] = 0.;
        DPars_IN["reDC9"] = 0.;
        DPars_IN["reDC10"] = 0.;
        DPars_IN["reDC7g"] = 0.;
        DPars_IN["reDC8g"] = 0.;
        DPars_IN["imDC1"] = 0.;
        DPars_IN["imDC2"] = 0.;
        DPars_IN["imDC3"] = 0.;
        DPars_IN["imDC4"] = 0.;
        DPars_IN["imDC5"] = 0.;
        DPars_IN["imDC6"] = 0.;
        DPars_IN["imDC7"] = 0.;
        DPars_IN["imDC8"] = 0.;
        DPars_IN["imDC9"] = 0.;
        DPars_IN["imDC10"] = 0.;
        DPars_IN["imDC7g"] = 0.;
        DPars_IN["imDC8g"] = 0.;
        DPars_IN["reDC7p"] = 0.;
        DPars_IN["reDC8p"] = 0.;
        DPars_IN["reDC9p"] = 0.;
        DPars_IN["reDC10p"] = 0.;
        DPars_IN["reDC7gp"] = 0.;
        DPars_IN["reDC8gp"] = 0.;
        DPars_IN["imDC7p"] = 0.;
        DPars_IN["imDC8p"] = 0.;
        DPars_IN["imDC9p"] = 0.;
        DPars_IN["imDC10p"] = 0.;
        DPars_IN["imDC7gp"] = 0.;
        DPars_IN["imDC8gp"] = 0.;
        DPars_IN["WCscale"] = 0.;
        
        return (DPars_IN);
    };
    
    /**
     * @brief A method that generates the map of the mandatory model 
     * parameters for FlavourWilsonCoefficient_DF2
     * @return the map of the mandatory parameters
     */
    std::map<std::string, double> NPFlavourWilsonCoefficient_DF2()
    {
        DPars_IN = StandardModel();

        DPars_IN["reC1_s"] = 0.;
        DPars_IN["reC2_s"] = 0.;
        DPars_IN["reC3_s"] = 0.;
        DPars_IN["reC4_s"] = 0.;
        DPars_IN["reC5_s"] = 0.;
        DPars_IN["imC1_s"] = 0.;
        DPars_IN["imC2_s"] = 0.;
        DPars_IN["imC3_s"] = 0.;
        DPars_IN["imC4_s"] = 0.;
        DPars_IN["imC5_s"] = 0.;
        DPars_IN["WCscale_s"] = 0.;
        DPars_IN["reC1_c"] = 0.;
        DPars_IN["reC2_c"] = 0.;
        DPars_IN["reC3_c"] = 0.;
        DPars_IN["reC4_c"] = 0.;
        DPars_IN["reC5_c"] = 0.;
        DPars_IN["imC1_c"] = 0.;
        DPars_IN["imC2_c"] = 0.;
        DPars_IN["imC3_c"] = 0.;
        DPars_IN["imC4_c"] = 0.;
        DPars_IN["imC5_c"] = 0.;
        DPars_IN["WCscale_c"] = 0.;
        DPars_IN["reC1_bd"] = 0.;
        DPars_IN["reC2_bd"] = 0.;
        DPars_IN["reC3_bd"] = 0.;
        DPars_IN["reC4_bd"] = 0.;
        DPars_IN["reC5_bd"] = 0.;
        DPars_IN["imC1_bd"] = 0.;
        DPars_IN["imC2_bd"] = 0.;
        DPars_IN["imC3_bd"] = 0.;
        DPars_IN["imC4_bd"] = 0.;
        DPars_IN["imC5_bd"] = 0.;
        DPars_IN["WCscale_bd"] = 0.;
        DPars_IN["reC1_bs"] = 0.;
        DPars_IN["reC2_bs"] = 0.;
        DPars_IN["reC3_bs"] = 0.;
        DPars_IN["reC4_bs"] = 0.;
        DPars_IN["reC5_bs"] = 0.;
        DPars_IN["imC1_bs"] = 0.;
        DPars_IN["imC2_bs"] = 0.;
        DPars_IN["imC3_bs"] = 0.;
        DPars_IN["imC4_bs"] = 0.;
        DPars_IN["imC5_bs"] = 0.;
        DPars_IN["WCscale_bs"] = 0.;
        
        return (DPars_IN);
    };
    
    /**
     * @brief A method that generates the map of the mandatory model parameters for RealWeakEFTLFV
     * @return the map of the mandatory parameters
     */
    std::map<std::string, double> NPRealWeakEFTLFV()
    {
        DPars_IN = StandardModel();

        DPars_IN["C7"] = 0.;
        DPars_IN["C7p"] = 0.;
        DPars_IN["C8"] = 0.;
        DPars_IN["C8p"] = 0.;
        DPars_IN["C9_11"] = 0.;
        DPars_IN["C9p_11"] = 0.;
        DPars_IN["C10_11"] = 0.;
        DPars_IN["C10p_11"] = 0.;
        DPars_IN["CS_11"] = 0.;
        DPars_IN["CSp_11"] = 0.;
        DPars_IN["CP_11"] = 0.;
        DPars_IN["CPp_11"] = 0.;
        DPars_IN["C9_22"] = 0.;
        DPars_IN["C9p_22"] = 0.;
        DPars_IN["C10_22"] = 0.;
        DPars_IN["C10p_22"] = 0.;
        DPars_IN["CS_22"] = 0.;
        DPars_IN["CSp_22"] = 0.;
        DPars_IN["CP_22"] = 0.;
        DPars_IN["CPp_22"] = 0.;
        DPars_IN["WCscale"] = 0.;
        
        return (DPars_IN);
    };
    
    /**
     * @brief A method that generates the map of the mandatory model parameters for RealWeakEFTCC
     * @return the map of the mandatory parameters
     */
    std::map<std::string, double> NPRealWeakEFTCC()
    {
        DPars_IN = StandardModel();

        DPars_IN["CS"] = 0.;
        DPars_IN["CP"] = 0.;
        DPars_IN["CV"] = 0.;
        DPars_IN["CA"] = 0.;
        DPars_IN["CT"] = 0.;
        
        return (DPars_IN);
    };
    
    /**
     * @brief A method that generates the map of the mandatory model parameters for RealWeakEFTCCPM
     * @return the map of the mandatory parameters
     */
    std::map<std::string, double> NPRealWeakEFTCCPM()
    {
        DPars_IN = StandardModel();

        DPars_IN["CSL"] = 0.;
        DPars_IN["CSR"] = 0.;
        DPars_IN["CVL"] = 0.;
        DPars_IN["CVR"] = 0.;
        DPars_IN["CT"] = 0.;
        
        return (DPars_IN);
    };
    
    /**
     * @brief A method that generates the map of the mandatory model parameters for the THDM.
     * @return the map of the mandatory parameters
     */
    std::map<std::string, double> NPTHDM()
    {
        DPars_IN = StandardModel();

        DPars_IN["logtb"] = 0.;
        DPars_IN["bma"] = 0.;
        DPars_IN["mHh1"] = 0.;
        DPars_IN["mA1"] = 0.;
        DPars_IN["mHp1"] = 0.;
        DPars_IN["mHh2"] = 0.;
        DPars_IN["mA2"] = 0.;
        DPars_IN["mHp2"] = 0.;
        DPars_IN["m12_2"] = 0.;
        DPars_IN["BDtaunu_SM"] = 0.;
        DPars_IN["BDtaunu_A"] = 0.;
        DPars_IN["BDtaunu_B"] = 0.;
        DPars_IN["BDstartaunu_SM"] = 0.;
        DPars_IN["BDstartaunu_A"] = 0.;
        DPars_IN["BDstartaunu_B"] = 0.;
        DPars_IN["bsgamma_theoryerror"] = 0.;
        DPars_IN["Q_THDM"] = 0.;
        DPars_IN["Rpeps"] = 0.;
        DPars_IN["NLOuniscale"] = 0.;

        return (DPars_IN);
    };
            
    /**
     * @brief A method that generates the map of the mandatory model parameters for the GeorgiMachacek.
     * @return the map of the mandatory parameters
     */
    std::map<std::string, double> NPGeorgiMachacek()
    {
        DPars_IN = StandardModel();

        DPars_IN["vDelta"] = 0.;
        DPars_IN["alpha"] = 0.;
        DPars_IN["mHh"] = 0.;
        DPars_IN["mA"] = 0.;
        DPars_IN["mH5"] = 0.;
        DPars_IN["Mu1"] = 0.;
        DPars_IN["Mu2"] = 0.;
        DPars_IN["Q_GM"] = 0.;

        return (DPars_IN);
    };
     
    /**
     * @brief A method that generates the map of the mandatory model parameters for the GeneralTHDM.
     * @return the map of the mandatory parameters
     */
    std::map<std::string, double> NPGeneralTHDM()
    {
        DPars_IN = StandardModel();

        DPars_IN["logtb"] = 0.;
        DPars_IN["mHp1"] = 0.;
        DPars_IN["mH21"] = 0.;
        DPars_IN["mH31"] = 0.;
        DPars_IN["mHp2"] = 0.;
        DPars_IN["mH2sq"] = 0.;
        DPars_IN["mH3sq"] = 0.;
        DPars_IN["alpha1"] = 0.;
        DPars_IN["alpha2"] = 0.;
        DPars_IN["alpha3"] = 0.;
        DPars_IN["Relambda5"] = 0.;
        DPars_IN["Imlambda5"] = 0.;
        DPars_IN["Relambda6"] = 0.;
        DPars_IN["Relambda7"] = 0.;

        DPars_IN["yu1R_GTHDM"] = 0.;
        DPars_IN["yd1R_GTHDM"] = 0.;
        DPars_IN["yl1R_GTHDM"] = 0.;



        DPars_IN["Nu_11r"] = 0.;
        DPars_IN["Nu_11i"] = 0.;
        DPars_IN["Nu_12r"] = 0.;
        DPars_IN["Nu_12i"] = 0.;
        DPars_IN["Nu_13r"] = 0.;
        DPars_IN["Nu_13i"] = 0.;
        DPars_IN["Nu_21r"] = 0.;
        DPars_IN["Nu_21i"] = 0.;
        DPars_IN["Nu_22r"] = 0.;
        DPars_IN["Nu_22i"] = 0.;
        DPars_IN["Nu_23r"] = 0.;
        DPars_IN["Nu_23i"] = 0.;
        DPars_IN["Nu_31r"] = 0.;
        DPars_IN["Nu_31i"] = 0.;
        DPars_IN["Nu_32r"] = 0.;
        DPars_IN["Nu_32i"] = 0.;
        DPars_IN["Nu_33r"] = 0.;
        DPars_IN["Nu_33i"] = 0.;
        DPars_IN["Nd_11r"] = 0.;
        DPars_IN["Nd_11i"] = 0.;
        DPars_IN["Nd_12r"] = 0.;
        DPars_IN["Nd_12i"] = 0.;
        DPars_IN["Nd_13r"] = 0.;
        DPars_IN["Nd_13i"] = 0.;
        DPars_IN["Nd_21r"] = 0.;
        DPars_IN["Nd_21i"] = 0.;
        DPars_IN["Nd_22r"] = 0.;
        DPars_IN["Nd_22i"] = 0.;
        DPars_IN["Nd_23r"] = 0.;
        DPars_IN["Nd_23i"] = 0.;
        DPars_IN["Nd_31r"] = 0.;
        DPars_IN["Nd_31i"] = 0.;
        DPars_IN["Nd_32r"] = 0.;
        DPars_IN["Nd_32i"] = 0.;
        DPars_IN["Nd_33r"] = 0.;
        DPars_IN["Nd_33i"] = 0.;
        DPars_IN["Nl_11r"] = 0.;
        DPars_IN["Nl_11i"] = 0.;
        DPars_IN["Nl_12r"] = 0.;
        DPars_IN["Nl_12i"] = 0.;
        DPars_IN["Nl_13r"] = 0.;
        DPars_IN["Nl_13i"] = 0.;
        DPars_IN["Nl_21r"] = 0.;
        DPars_IN["Nl_21i"] = 0.;
        DPars_IN["Nl_22r"] = 0.;
        DPars_IN["Nl_22i"] = 0.;
        DPars_IN["Nl_23r"] = 0.;
        DPars_IN["Nl_23i"] = 0.;
        DPars_IN["Nl_31r"] = 0.;
        DPars_IN["Nl_31i"] = 0.;
        DPars_IN["Nl_32r"] = 0.;
        DPars_IN["Nl_32i"] = 0.;
        DPars_IN["Nl_33r"] = 0.;
        DPars_IN["Nl_33i"] = 0.;
        DPars_IN["Q_GTHDM"] = 0.;
        DPars_IN["RpepsGTHDM"] = 0.;
        DPars_IN["NLOuniscaleGTHDM"] = 0.;

        return (DPars_IN);
    };
    
    /**
     * @brief A method that generates the map of the mandatory model parameters for the THDMW.
     * @return the map of the mandatory parameters
     */
    std::map<std::string, double> NPTHDMW()
    {
        DPars_IN = StandardModel();

        DPars_IN["THDMW_logtb"] = 0.;
        DPars_IN["THDMW_bma"] = 0.;
        DPars_IN["THDMW_lambda1"] = 0.;
        DPars_IN["THDMW_lambda2"] = 0.;
        DPars_IN["THDMW_lambda3"] = 0.;
        DPars_IN["THDMW_lambda4"] = 0.;
        DPars_IN["THDMW_lambda5"] = 0.;
        DPars_IN["THDMW_mS2"] = 0.;
        DPars_IN["THDMW_mu1"] = 0.;
        DPars_IN["THDMW_mu2"] = 0.;
        DPars_IN["THDMW_mu3"] = 0.;
        DPars_IN["THDMW_mu4"] = 0.;
        DPars_IN["THDMW_mu5"] = 0.;
        DPars_IN["THDMW_mu6"] = 0.;
        DPars_IN["THDMW_nu1"] = 0.;
        DPars_IN["THDMW_nu2"] = 0.;
        DPars_IN["THDMW_nu3"] = 0.;
        DPars_IN["THDMW_nu4"] = 0.;
        DPars_IN["THDMW_nu5"] = 0.;
        DPars_IN["THDMW_omega1"] = 0.;
        DPars_IN["THDMW_omega2"] = 0.;
        DPars_IN["THDMW_omega3"] = 0.;
        DPars_IN["THDMW_omega4"] = 0.;
        DPars_IN["THDMW_omega5"] = 0.;
        DPars_IN["THDMW_kappa1"] = 0.;
        DPars_IN["THDMW_kappa2"] = 0.;
        DPars_IN["THDMW_kappa3"] = 0.;
        DPars_IN["THDMW_etaU"] = 0.;
        DPars_IN["THDMW_etaD"] = 0.;
        DPars_IN["THDMW_rho_b"] = 0.;
        DPars_IN["THDMW_S_b"] = 0.;\
        DPars_IN["THDMW_kappa3"] = 0.;
        DPars_IN["Q_THDMW"] = 0.;
        DPars_IN["RpepsTHDMW"] = 0.;
        DPars_IN["NLOuniscaleTHDMW"] = 0.;

        return (DPars_IN);
    };
    
    /**
     * @brief A method that generates the map of the mandatory model parameters for the SUSY.
     * @return the map of the mandatory parameters
     */
    std::map<std::string, double> NPSUSY()
    {
        DPars_IN = StandardModel();

        DPars_IN["m1r"] = 0.;
        DPars_IN["m1i"] = 0.;
        DPars_IN["m2r"] = 0.;
        DPars_IN["m2i"] = 0.;
        DPars_IN["m3"] = 0.;\
        DPars_IN["muHr"] = 0.;
        DPars_IN["muHi"] = 0.;
        DPars_IN["mHptree"] = 0.;
        DPars_IN["tanb"] = 0.;
        DPars_IN["Q_SUSY"] = 0.;

        return (DPars_IN);
    };
    
    /**
     * @brief A method that generates the map of the mandatory model parameters for the GeneralSUSY.
     * @return the map of the mandatory parameters
     */
    std::map<std::string, double> NPGeneralSUSY()
    {
        DPars_IN = StandardModel();
        
        DPars_IN.insert(NPSUSY().begin(), NPSUSY().end());

        DPars_IN["msQhat2_11r"] = 0.;
        DPars_IN["msQhat2_12r"] = 0.;
        DPars_IN["msQhat2_12i"] = 0.;
        DPars_IN["msQhat2_13r"] = 0.;
        DPars_IN["msQhat2_13i"] = 0.;
        DPars_IN["msQhat2_22r"] = 0.;
        DPars_IN["msQhat2_23r"] = 0.;
        DPars_IN["msQhat2_23i"] = 0.;
        DPars_IN["msQhat2_33r"] = 0.;
        DPars_IN["msUhat2_11r"] = 0.;
        DPars_IN["msUhat2_12r"] = 0.;
        DPars_IN["msUhat2_12i"] = 0.;
        DPars_IN["msUhat2_13r"] = 0.;
        DPars_IN["msUhat2_13i"] = 0.;
        DPars_IN["msUhat2_22r"] = 0.;
        DPars_IN["msUhat2_23r"] = 0.;
        DPars_IN["msUhat2_23i"] = 0.;
        DPars_IN["msUhat2_33r"] = 0.;
        DPars_IN["msDhat2_11r"] = 0.;
        DPars_IN["msDhat2_12r"] = 0.;
        DPars_IN["msDhat2_12i"] = 0.;
        DPars_IN["msDhat2_13r"] = 0.;
        DPars_IN["msDhat2_13i"] = 0.;
        DPars_IN["msDhat2_22r"] = 0.;
        DPars_IN["msDhat2_23r"] = 0.;
        DPars_IN["msDhat2_23i"] = 0.;
        DPars_IN["msDhat2_33r"] = 0.;
        DPars_IN["msLhat2_11r"] = 0.;
        DPars_IN["msLhat2_12r"] = 0.;
        DPars_IN["msLhat2_12i"] = 0.;
        DPars_IN["msLhat2_13r"] = 0.;
        DPars_IN["msLhat2_13i"] = 0.;
        DPars_IN["msLhat2_22r"] = 0.;
        DPars_IN["msLhat2_23r"] = 0.;
        DPars_IN["msLhat2_23i"] = 0.;
        DPars_IN["msLhat2_33r"] = 0.;
        DPars_IN["msEhat2_11r"] = 0.;
        DPars_IN["msEhat2_12r"] = 0.;
        DPars_IN["msEhat2_12i"] = 0.;
        DPars_IN["msEhat2_13r"] = 0.;
        DPars_IN["msEhat2_13i"] = 0.;
        DPars_IN["msEhat2_22r"] = 0.;
        DPars_IN["msEhat2_23r"] = 0.;
        DPars_IN["msEhat2_23i"] = 0.;
        DPars_IN["msEhat2_33r"] = 0.;
        DPars_IN["msNhat2_11r"] = 0.;
        DPars_IN["msNhat2_12r"] = 0.;
        DPars_IN["msNhat2_12i"] = 0.;
        DPars_IN["msNhat2_13r"] = 0.;
        DPars_IN["msNhat2_13i"] = 0.;
        DPars_IN["msNhat2_22r"] = 0.;
        DPars_IN["msNhat2_23r"] = 0.;
        DPars_IN["msNhat2_23i"] = 0.;
        DPars_IN["msNhat2_33r"] = 0.;
        DPars_IN["TUhat_11r"] = 0.;
        DPars_IN["TUhat_12r"] = 0.;
        DPars_IN["TUhat_13r"] = 0.;
        DPars_IN["TUhat_21r"] = 0.;
        DPars_IN["TUhat_22r"] = 0.;
        DPars_IN["TUhat_23r"] = 0.;
        DPars_IN["TUhat_31r"] = 0.;
        DPars_IN["TUhat_32r"] = 0.;
        DPars_IN["TUhat_33r"] = 0.;
        DPars_IN["TUhat_11i"] = 0.;
        DPars_IN["TUhat_12i"] = 0.;
        DPars_IN["TUhat_13i"] = 0.;
        DPars_IN["TUhat_21i"] = 0.;
        DPars_IN["TUhat_22i"] = 0.;
        DPars_IN["TUhat_23i"] = 0.;
        DPars_IN["TUhat_31i"] = 0.;
        DPars_IN["TUhat_32i"] = 0.;
        DPars_IN["TUhat_33i"] = 0.;
        DPars_IN["TDhat_11r"] = 0.;
        DPars_IN["TDhat_12r"] = 0.;
        DPars_IN["TDhat_13r"] = 0.;
        DPars_IN["TDhat_21r"] = 0.;
        DPars_IN["TDhat_22r"] = 0.;
        DPars_IN["TDhat_23r"] = 0.;
        DPars_IN["TDhat_31r"] = 0.;
        DPars_IN["TDhat_32r"] = 0.;
        DPars_IN["TDhat_33r"] = 0.;
        DPars_IN["TDhat_11i"] = 0.;
        DPars_IN["TDhat_12i"] = 0.;
        DPars_IN["TDhat_13i"] = 0.;
        DPars_IN["TDhat_21i"] = 0.;
        DPars_IN["TDhat_22i"] = 0.;
        DPars_IN["TDhat_23i"] = 0.;
        DPars_IN["TDhat_31i"] = 0.;
        DPars_IN["TDhat_32i"] = 0.;
        DPars_IN["TDhat_33i"] = 0.;
        DPars_IN["TEhat_11r"] = 0.;
        DPars_IN["TEhat_12r"] = 0.;
        DPars_IN["TEhat_13r"] = 0.;
        DPars_IN["TEhat_21r"] = 0.;
        DPars_IN["TEhat_22r"] = 0.;
        DPars_IN["TEhat_23r"] = 0.;
        DPars_IN["TEhat_31r"] = 0.;
        DPars_IN["TEhat_32r"] = 0.;
        DPars_IN["TEhat_33r"] = 0.;
        DPars_IN["TEhat_11i"] = 0.;
        DPars_IN["TEhat_12i"] = 0.;
        DPars_IN["TEhat_13i"] = 0.;
        DPars_IN["TEhat_21i"] = 0.;
        DPars_IN["TEhat_22i"] = 0.;
        DPars_IN["TEhat_23i"] = 0.;
        DPars_IN["TEhat_31i"] = 0.;
        DPars_IN["TEhat_32i"] = 0.;
        DPars_IN["TEhat_33i"] = 0.;
        DPars_IN["TNhat_11r"] = 0.;
        DPars_IN["TNhat_12r"] = 0.;
        DPars_IN["TNhat_13r"] = 0.;
        DPars_IN["TNhat_21r"] = 0.;
        DPars_IN["TNhat_22r"] = 0.;
        DPars_IN["TNhat_23r"] = 0.;
        DPars_IN["TNhat_31r"] = 0.;
        DPars_IN["TNhat_32r"] = 0.;
        DPars_IN["TNhat_33r"] = 0.;
        DPars_IN["TNhat_11i"] = 0.;
        DPars_IN["TNhat_12i"] = 0.;
        DPars_IN["TNhat_13i"] = 0.;
        DPars_IN["TNhat_21i"] = 0.;
        DPars_IN["TNhat_22i"] = 0.;
        DPars_IN["TNhat_23i"] = 0.;
        DPars_IN["TNhat_31i"] = 0.;
        DPars_IN["TNhat_32i"] = 0.;
        DPars_IN["TNhat_33i"] = 0.;

        return (DPars_IN);
    };
    
    std::map<std::string, double> DPars_IN;///< A map for the list of mandatory parameters in the model being used.
};

#endif	/* INPUTPARAMETERS_H */
    