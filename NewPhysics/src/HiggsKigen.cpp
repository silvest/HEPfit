/* 
 * Copyright (C) 2014 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "HiggsKigen.h"

const std::string HiggsKigen::HKvKfgenvars[NHKvKfgenvars] = {
    "Kw", "Kz", "Kg", "Kga", "Kzga", "Ku", "Kc", "Kt", "Kd", "Ks", "Kb", "Ke", "Kmu", "Ktau", "KH", "BrHinv", "BrHexo",
    "eggFint","eggFpar","ettHint","ettHpar",
    "eVBFint","eVBFpar","eWHint","eWHpar","eZHint","eZHpar",
    "eeeWBFint","eeeWBFpar","eeeZHint","eeeZHpar","eeettHint","eeettHpar",
    "eepWBFint","eepWBFpar","eepZBFint","eepZBFpar",
    "eHggint","eHggpar","eHWWint","eHWWpar","eHZZint","eHZZpar","eHZgaint","eHZgapar",
    "eHgagaint","eHgagapar","eHmumuint","eHmumupar","eHtautauint","eHtautaupar",
    "eHccint","eHccpar","eHbbint","eHbbpar",
    "eggFHgaga","eggFHZga","eggFHZZ","eggFHWW","eggFHtautau","eggFHbb","eggFHmumu",   
    "eVBFHgaga","eVBFHZga","eVBFHZZ","eVBFHWW","eVBFHtautau","eVBFHbb","eVBFHmumu",   
    "eWHgaga","eWHZga","eWHZZ","eWHWW","eWHtautau","eWHbb","eWHmumu",    
    "eZHgaga","eZHZga","eZHZZ","eZHWW","eZHtautau","eZHbb","eZHmumu",
    "ettHgaga","ettHZga","ettHZZ","ettHWW","ettHtautau","ettHbb","ettHmumu",
    "eVBFHinv","eVHinv"
};

HiggsKigen::HiggsKigen()
: NPbase()
{
    FlagKiLoop = true;
    FlagCustodial = false;
    FlagUniversalKf = false;
    FlagUniversalK = false;
    FlagUseKH = false;
    
    ModelParamMap.insert(std::make_pair("Kw", std::cref(Kw)));
    ModelParamMap.insert(std::make_pair("Kz", std::cref(Kz)));
    ModelParamMap.insert(std::make_pair("Kg", std::cref(Kg)));
    ModelParamMap.insert(std::make_pair("Kga", std::cref(Kga)));
    ModelParamMap.insert(std::make_pair("Kzga", std::cref(Kzga)));
    ModelParamMap.insert(std::make_pair("Ku", std::cref(Ku)));
    ModelParamMap.insert(std::make_pair("Kc", std::cref(Kc)));
    ModelParamMap.insert(std::make_pair("Kt", std::cref(Kt)));
    ModelParamMap.insert(std::make_pair("Kd", std::cref(Kd)));
    ModelParamMap.insert(std::make_pair("Ks", std::cref(Ks)));
    ModelParamMap.insert(std::make_pair("Kb", std::cref(Kb)));
    ModelParamMap.insert(std::make_pair("Ke", std::cref(Ke)));
    ModelParamMap.insert(std::make_pair("Kmu", std::cref(Kmu)));
    ModelParamMap.insert(std::make_pair("Ktau", std::cref(Ktau)));
    ModelParamMap.insert(std::make_pair("KH", std::cref(KH)));
    ModelParamMap.insert(std::make_pair("BrHinv", std::cref(BrHinv)));
    ModelParamMap.insert(std::make_pair("BrHexo", std::cref(BrHexo)));
    ModelParamMap.insert(std::make_pair("eggFint", std::cref(eggFint)));
    ModelParamMap.insert(std::make_pair("eggFpar", std::cref(eggFpar)));
    ModelParamMap.insert(std::make_pair("ettHint", std::cref(ettHint)));
    ModelParamMap.insert(std::make_pair("ettHpar", std::cref(ettHpar)));
    ModelParamMap.insert(std::make_pair("eVBFint", std::cref(eVBFint)));
    ModelParamMap.insert(std::make_pair("eVBFpar", std::cref(eVBFpar)));
    ModelParamMap.insert(std::make_pair("eWHint", std::cref(eWHint)));
    ModelParamMap.insert(std::make_pair("eWHpar", std::cref(eWHpar)));
    ModelParamMap.insert(std::make_pair("eZHint", std::cref(eZHint)));
    ModelParamMap.insert(std::make_pair("eZHpar", std::cref(eZHpar)));
    ModelParamMap.insert(std::make_pair("eeeWBFint", std::cref(eeeWBFint)));
    ModelParamMap.insert(std::make_pair("eeeWBFpar", std::cref(eeeWBFpar)));
    ModelParamMap.insert(std::make_pair("eeeZHint", std::cref(eeeZHint)));
    ModelParamMap.insert(std::make_pair("eeeZHpar", std::cref(eeeZHpar)));
    ModelParamMap.insert(std::make_pair("eeettHint", std::cref(eeettHint)));
    ModelParamMap.insert(std::make_pair("eeettHpar", std::cref(eeettHpar)));
    ModelParamMap.insert(std::make_pair("eepWBFint", std::cref(eepWBFint)));
    ModelParamMap.insert(std::make_pair("eepWBFpar", std::cref(eepWBFpar)));
    ModelParamMap.insert(std::make_pair("eepZBFint", std::cref(eepZBFint)));
    ModelParamMap.insert(std::make_pair("eepZBFpar", std::cref(eepZBFpar)));
    ModelParamMap.insert(std::make_pair("eHggint", std::cref(eHggint)));
    ModelParamMap.insert(std::make_pair("eHggpar", std::cref(eHggpar)));
    ModelParamMap.insert(std::make_pair("eHWWint", std::cref(eHWWint)));
    ModelParamMap.insert(std::make_pair("eHWWpar", std::cref(eHWWpar)));
    ModelParamMap.insert(std::make_pair("eHZZint", std::cref(eHZZint)));
    ModelParamMap.insert(std::make_pair("eHZZpar", std::cref(eHZZpar)));
    ModelParamMap.insert(std::make_pair("eHZgaint", std::cref(eHZgaint)));
    ModelParamMap.insert(std::make_pair("eHZgapar", std::cref(eHZgapar)));
    ModelParamMap.insert(std::make_pair("eHgagaint", std::cref(eHgagaint)));
    ModelParamMap.insert(std::make_pair("eHgagapar", std::cref(eHgagapar)));
    ModelParamMap.insert(std::make_pair("eHmumuint", std::cref(eHmumuint)));
    ModelParamMap.insert(std::make_pair("eHmumupar", std::cref(eHmumupar)));
    ModelParamMap.insert(std::make_pair("eHtautauint", std::cref(eHtautauint)));
    ModelParamMap.insert(std::make_pair("eHtautaupar", std::cref(eHtautaupar)));
    ModelParamMap.insert(std::make_pair("eHccint", std::cref(eHccint)));
    ModelParamMap.insert(std::make_pair("eHccpar", std::cref(eHccpar)));
    ModelParamMap.insert(std::make_pair("eHbbint", std::cref(eHbbint)));
    ModelParamMap.insert(std::make_pair("eHbbpar", std::cref(eHbbpar)));    
    ModelParamMap.insert(std::make_pair("eggFHgaga", std::cref(eggFHgaga)));
    ModelParamMap.insert(std::make_pair("eggFHZga", std::cref(eggFHZga)));
    ModelParamMap.insert(std::make_pair("eggFHZZ", std::cref(eggFHZZ)));
    ModelParamMap.insert(std::make_pair("eggFHWW", std::cref(eggFHWW)));
    ModelParamMap.insert(std::make_pair("eggFHtautau", std::cref(eggFHtautau)));
    ModelParamMap.insert(std::make_pair("eggFHbb", std::cref(eggFHbb)));
    ModelParamMap.insert(std::make_pair("eggFHmumu", std::cref(eggFHmumu)));   
    ModelParamMap.insert(std::make_pair("eVBFHgaga", std::cref(eVBFHgaga)));
    ModelParamMap.insert(std::make_pair("eVBFHZga", std::cref(eVBFHZga)));
    ModelParamMap.insert(std::make_pair("eVBFHZZ", std::cref(eVBFHZZ)));
    ModelParamMap.insert(std::make_pair("eVBFHWW", std::cref(eVBFHWW)));
    ModelParamMap.insert(std::make_pair("eVBFHtautau", std::cref(eVBFHtautau)));
    ModelParamMap.insert(std::make_pair("eVBFHbb", std::cref(eVBFHbb)));
    ModelParamMap.insert(std::make_pair("eVBFHmumu", std::cref(eVBFHmumu)));   
    ModelParamMap.insert(std::make_pair("eWHgaga", std::cref(eWHgaga)));
    ModelParamMap.insert(std::make_pair("eWHZga", std::cref(eWHZga)));
    ModelParamMap.insert(std::make_pair("eWHZZ", std::cref(eWHZZ)));
    ModelParamMap.insert(std::make_pair("eWHWW", std::cref(eWHWW)));
    ModelParamMap.insert(std::make_pair("eWHtautau", std::cref(eWHtautau)));
    ModelParamMap.insert(std::make_pair("eWHbb", std::cref(eWHbb)));
    ModelParamMap.insert(std::make_pair("eWHmumu", std::cref(eWHmumu)));    
    ModelParamMap.insert(std::make_pair("eZHgaga", std::cref(eZHgaga)));
    ModelParamMap.insert(std::make_pair("eZHZga", std::cref(eZHZga)));
    ModelParamMap.insert(std::make_pair("eZHZZ", std::cref(eZHZZ)));
    ModelParamMap.insert(std::make_pair("eZHWW", std::cref(eZHWW)));
    ModelParamMap.insert(std::make_pair("eZHtautau", std::cref(eZHtautau)));
    ModelParamMap.insert(std::make_pair("eZHbb", std::cref(eZHbb)));
    ModelParamMap.insert(std::make_pair("eZHmumu", std::cref(eZHmumu)));
    ModelParamMap.insert(std::make_pair("ettHgaga", std::cref(ettHgaga)));
    ModelParamMap.insert(std::make_pair("ettHZga", std::cref(ettHZga)));
    ModelParamMap.insert(std::make_pair("ettHZZ", std::cref(ettHZZ)));
    ModelParamMap.insert(std::make_pair("ettHWW", std::cref(ettHWW)));
    ModelParamMap.insert(std::make_pair("ettHtautau", std::cref(ettHtautau)));
    ModelParamMap.insert(std::make_pair("ettHbb", std::cref(ettHbb)));
    ModelParamMap.insert(std::make_pair("ettHmumu", std::cref(ettHmumu)));
    ModelParamMap.insert(std::make_pair("eVBFHinv", std::cref(eVBFHinv)));
    ModelParamMap.insert(std::make_pair("eVHinv", std::cref(eVHinv)));
}

bool HiggsKigen::PostUpdate()
{
    if (!NPbase::PostUpdate()) return (false);

//  Check first the flags that control the values of the parameters 
    
    if (FlagCustodial) {
//  Assign to all Kz the value of Kw
        Kz = Kw;
    }
    
    if (FlagUniversalKf) {
//  Assign to all Kf the value of Kt
        Ku = Kt;
        Kc = Kt;
        Kd = Kt;
        Ks = Kt;
        Kb = Kt;
        Ke = Kt;
        Kmu = Kt;
        Ktau = Kt;
    }
    
    if (FlagUniversalK) {
//  Assign to all K the value of Kt 
        Kz = Kt;
        Kw = Kt;
        
        Ku = Kt;
        Kc = Kt;
        Kd = Kt;
        Ks = Kt;
        Kb = Kt;
        Ke = Kt;
        Kmu = Kt;
        Ktau = Kt;

//  If the following are given by modifications of SM loops only they will be replaced below        
        Kg = Kt;
        Kga = Kt;
        Kzga = Kt;
    }
    
//  Calculation of some quantities repeteadly used in the code
    
//  Then the flag to add the values of the loops to the cache  
    
    if (!FlagKiLoop) {
        Kg = computeKgLoop();
        Kga = computeKgagaLoop();
        Kzga = computeKZgaLoop();
    }
    
//  Ratio of the total Higgs width with respect to SM and Exotic BR
    if (FlagUseKH) {
        GammaHTotR = KH*KH;  
        GammaHTotSMR = computeGammaTotalSMRatio();
        
        GammaHTotnoER = KH*KH;  
        GammaHTotSMnoER = computeGammaTotalSMnoErrorRatio();
        
        //  Compute BrHexo from KH
        BrHexo = 1.0 - BrHinv - GammaHTotSMR/GammaHTotR;     
        
        //  BRHinv is positive by definition and the bound is controlled by the corresponding parameter. 
        //  The exotic must be in [0,1]. Demand that 0<BrHexo<1.
        //  Otherwise, exclude the point from likelihood by assigning a nan to all parameters
        if ( (BrHexo > 1.0) || (BrHexo < 0.0) ) {

            GammaHTotR = std::numeric_limits<double>::quiet_NaN();

            Kw = std::numeric_limits<double>::quiet_NaN();            
            Kz = std::numeric_limits<double>::quiet_NaN();
            Kg = std::numeric_limits<double>::quiet_NaN();
            Kga = std::numeric_limits<double>::quiet_NaN();
            Kzga = std::numeric_limits<double>::quiet_NaN();        
            Ku = std::numeric_limits<double>::quiet_NaN();
            Kc = std::numeric_limits<double>::quiet_NaN();
            Kt = std::numeric_limits<double>::quiet_NaN();
            Kd = std::numeric_limits<double>::quiet_NaN();
            Ks = std::numeric_limits<double>::quiet_NaN();
            Kb = std::numeric_limits<double>::quiet_NaN();
            Ke = std::numeric_limits<double>::quiet_NaN();
            Kmu = std::numeric_limits<double>::quiet_NaN();
            Ktau = std::numeric_limits<double>::quiet_NaN();  
            KH = std::numeric_limits<double>::quiet_NaN();
            
            BrHinv = std::numeric_limits<double>::quiet_NaN();
            BrHexo = std::numeric_limits<double>::quiet_NaN();
        }        
        
    } else {
        GammaHTotR = computeGammaTotalRatio();
        GammaHTotSMR = computeGammaTotalSMRatio();
        
        GammaHTotnoER = computeGammaTotalnoErrorRatio();
        GammaHTotSMnoER = computeGammaTotalSMnoErrorRatio();
    }
        
    // The total theory error in the H width
    eHwidth = ( GammaHTotR - GammaHTotnoER ) / GammaHTotnoER;

    return (true);
}

void HiggsKigen::setParameter(const std::string name, const double& value)
{
    if (name.compare("Kw") == 0) {
        Kw = value;
    } else if (name.compare("Kz") == 0) {
        Kz = value;
    } else if (name.compare("Kg") == 0) {
        Kg = value;
    } else if (name.compare("Kga") == 0) {
        Kga = value;
    } else if (name.compare("Kzga") == 0) {
        Kzga = value;
    } else if (name.compare("Ku") == 0) {
        Ku = value;
    } else if (name.compare("Kc") == 0) {
        Kc = value;
    } else if (name.compare("Kt") == 0) {
        Kt = value;
    } else if (name.compare("Kd") == 0) {
        Kd = value;
    } else if (name.compare("Ks") == 0) {
        Ks = value;
    } else if (name.compare("Kb") == 0) {
        Kb = value;
    } else if (name.compare("Ke") == 0) {
        Ke = value;
    } else if (name.compare("Kmu") == 0) {
        Kmu = value;
    } else if (name.compare("Ktau") == 0) {
        Ktau = value;
    } else if (name.compare("KH") == 0) {
        KH = value;
    } else if (name.compare("BrHinv") == 0) {
//  Always positive
        BrHinv = fabs(value);
    } else if (name.compare("BrHexo") == 0) {
//  Always positive
        BrHexo = fabs(value);
    } else if (name.compare("eggFint") == 0) {
        eggFint = value;
    } else if (name.compare("eggFpar") == 0) {
        eggFpar = value;
    } else if (name.compare("ettHint") == 0) {
        ettHint = value;
    } else if (name.compare("ettHpar") == 0) {
        ettHpar = value;
    } else if (name.compare("eVBFint") == 0) {
        eVBFint = value;
    } else if (name.compare("eVBFpar") == 0) {
        eVBFpar = value;
    } else if (name.compare("eWHint") == 0) {
        eWHint = value;
    } else if (name.compare("eWHpar") == 0) {
        eWHpar = value;
    } else if (name.compare("eZHint") == 0) {
        eZHint = value;
    } else if (name.compare("eZHpar") == 0) {
        eZHpar = value;
    } else if (name.compare("eeeWBFint") == 0) {
        eeeWBFint = value;
    } else if (name.compare("eeeWBFpar") == 0) {
        eeeWBFpar = value;
    } else if (name.compare("eeeZHint") == 0) {
        eeeZHint = value;
    } else if (name.compare("eeeZHpar") == 0) {
        eeeZHpar = value;
    } else if (name.compare("eeettHint") == 0) {
        eeettHint = value;
    } else if (name.compare("eeettHpar") == 0) {
        eeettHpar = value;
    } else if (name.compare("eepWBFint") == 0) {
        eepWBFint = value;
    } else if (name.compare("eepWBFpar") == 0) {
        eepWBFpar = value;
    } else if (name.compare("eepZBFint") == 0) {
        eepZBFint = value;
    } else if (name.compare("eepZBFpar") == 0) {
        eepZBFpar = value;
    } else if (name.compare("eHggint") == 0) {
        eHggint = value;
    } else if (name.compare("eHggpar") == 0) {
        eHggpar = value;
    } else if (name.compare("eHWWint") == 0) {
        eHWWint = value;
    } else if (name.compare("eHWWpar") == 0) {
        eHWWpar = value;
    } else if (name.compare("eHZZint") == 0) {
        eHZZint = value;
    } else if (name.compare("eHZZpar") == 0) {
        eHZZpar = value;
    } else if (name.compare("eHZgaint") == 0) {
        eHZgaint = value;
    } else if (name.compare("eHZgapar") == 0) {
        eHZgapar = value;
    } else if (name.compare("eHgagaint") == 0) {
        eHgagaint = value;
    } else if (name.compare("eHgagapar") == 0) {
        eHgagapar = value;
    } else if (name.compare("eHmumuint") == 0) {
        eHmumuint = value;
    } else if (name.compare("eHmumupar") == 0) {
        eHmumupar = value;
    } else if (name.compare("eHtautauint") == 0) {
        eHtautauint = value;
    } else if (name.compare("eHtautaupar") == 0) {
        eHtautaupar = value;
    } else if (name.compare("eHccint") == 0) {
        eHccint = value;
    } else if (name.compare("eHccpar") == 0) {
        eHccpar = value;
    } else if (name.compare("eHbbint") == 0) {
        eHbbint = value;
    } else if (name.compare("eHbbpar") == 0) {
        eHbbpar = value;   
    } else if (name.compare("eggFHgaga") == 0) {
        eggFHgaga = value;
    } else if (name.compare("eggFHZga") == 0) {
        eggFHZga = value;
    } else if (name.compare("eggFHZZ") == 0) {
        eggFHZZ = value;
    } else if (name.compare("eggFHWW") == 0) {
        eggFHWW = value;
    } else if (name.compare("eggFHtautau") == 0) {
        eggFHtautau = value;
    } else if (name.compare("eggFHbb") == 0) {
        eggFHbb = value;
    } else if (name.compare("eggFHmumu") == 0) {
        eggFHmumu = value;     
    } else if (name.compare("eVBFHgaga") == 0) {
        eVBFHgaga = value;
    } else if (name.compare("eVBFHZga") == 0) {
        eVBFHZga = value;
    } else if (name.compare("eVBFHZZ") == 0) {
        eVBFHZZ = value;
    } else if (name.compare("eVBFHWW") == 0) {
        eVBFHWW = value;
    } else if (name.compare("eVBFHtautau") == 0) {
        eVBFHtautau = value;
    } else if (name.compare("eVBFHbb") == 0) {
        eVBFHbb = value;
    } else if (name.compare("eVBFHmumu") == 0) {
        eVBFHmumu = value;  
    } else if (name.compare("eWHgaga") == 0) {
        eWHgaga = value;
    } else if (name.compare("eWHZga") == 0) {
        eWHZga = value;
    } else if (name.compare("eWHZZ") == 0) {
        eWHZZ = value;
    } else if (name.compare("eWHWW") == 0) {
        eWHWW = value;
    } else if (name.compare("eWHtautau") == 0) {
        eWHtautau = value;
    } else if (name.compare("eWHbb") == 0) {
        eWHbb = value;
    } else if (name.compare("eWHmumu") == 0) {
        eWHmumu = value; 
    } else if (name.compare("eZHgaga") == 0) {
        eZHgaga = value;
    } else if (name.compare("eZHZga") == 0) {
        eZHZga = value;
    } else if (name.compare("eZHZZ") == 0) {
        eZHZZ = value;
    } else if (name.compare("eZHWW") == 0) {
        eZHWW = value;
    } else if (name.compare("eZHtautau") == 0) {
        eZHtautau = value;
    } else if (name.compare("eZHbb") == 0) {
        eZHbb = value;
    } else if (name.compare("eZHmumu") == 0) {
        eZHmumu = value; 
    } else if (name.compare("ettHgaga") == 0) {
        ettHgaga = value;
    } else if (name.compare("ettHZga") == 0) {
        ettHZga = value;
    } else if (name.compare("ettHZZ") == 0) {
        ettHZZ = value;
    } else if (name.compare("ettHWW") == 0) {
        ettHWW = value;
    } else if (name.compare("ettHtautau") == 0) {
        ettHtautau = value;
    } else if (name.compare("ettHbb") == 0) {
        ettHbb = value;
    } else if (name.compare("ettHmumu") == 0) {
        ettHmumu = value;
    } else if (name.compare("eVBFHinv") == 0) {
        eVBFHinv = value;
    } else if (name.compare("eVHinv") == 0) {
        eVHinv = value;
    } else
        NPbase::setParameter(name, value);
}

bool HiggsKigen::CheckParameters(const std::map<std::string, double>& DPars)
{
    for (int i = 0; i < NHKvKfgenvars; i++) {
        if (DPars.find(HKvKfgenvars[i]) == DPars.end()) {
            std::cout << "ERROR: missing mandatory HiggsKigen parameter " << HKvKfgenvars[i] << std::endl;
            raiseMissingModelParameterCount();
            addMissingModelParameter(HKvKfgenvars[i]);
        }
    }
    return (NPbase::CheckParameters(DPars));
}

bool HiggsKigen::setFlag(const std::string name, const bool value)
{
    bool res = false;
    if (name.compare("KiLoop") == 0) {
        FlagKiLoop = value;
        res = true;
    } else if (name.compare("Custodial") == 0) {
        FlagCustodial = value;
        res = true;
    } else if (name.compare("UniversalKf") == 0) {
        FlagUniversalKf = value;
        res = true;
    } else if (name.compare("UniversalK") == 0) {
        FlagUniversalK = value;
        res = true;
    } else if (name.compare("UseKH") == 0) {
        FlagUseKH = value;
        res = true;
    } else
        res = NPbase::setFlag(name, value);

    return (res);
}

////////////////////////////////////////////////////////////////////////

double HiggsKigen::obliqueS() const
{
    if (FlagCustodial) {
        double Lambda, Kv;
    
        Kv=Kw;
    
        if (fabs(1.0 - Kv * Kv) < pow(10.0, -32.0))
            Lambda = pow(10.0, 19.0);
        else
            Lambda = 4.0 * M_PI * v() / sqrt(fabs(1.0 - Kv * Kv));

        return ( 1.0 / 12.0 / M_PI * (1.0 - Kv * Kv) * log(Lambda * Lambda / mHl / mHl));
    
    } else {
        return 0.0;
    }
}

double HiggsKigen::obliqueT() const
{
    if (FlagCustodial) {
        double Lambda, Kv;
        double cW2_SM = trueSM.cW2();
        
        Kv=Kw;
        
        if (fabs(1.0 - Kv * Kv) < pow(10.0, -32.0))
            Lambda = pow(10.0, 19.0);
        else
            Lambda = 4.0 * M_PI * v() / sqrt(fabs(1.0 - Kv * Kv));

        return ( -3.0 / 16.0 / M_PI / cW2_SM * (1.0 - Kv * Kv) * log(Lambda * Lambda / mHl / mHl));
    
    } else {
        return 0.0;
    }
}

double HiggsKigen::obliqueU() const
{
    return 0.0;
}

////////////////////////////////////////////////////////////////////////

double HiggsKigen::muggH(const double sqrt_s) const
{
    double mu = 1.0;
    
    if (FlagKiLoop) {
        
        // H->gg width and ggH cross section assumed to scale in the same way ~Kg^2
        mu = (computeKg() * computeKg());

    } else {
        
        // Do not use the expression for Kg, which is computed from the width.
        // Use instead the following
        if (sqrt_s == 1.96){
        // Exact expression not available for Tevatron. Assume Kt^2 scaling.
            mu = (computeKt() * computeKt());
            
        } else {
            
            double sigmatt_SM = trueSM.computeSigmaggH_tt(sqrt_s);
            double sigmabb_SM = trueSM.computeSigmaggH_bb(sqrt_s);
            double sigmatb_SM = trueSM.computeSigmaggH_tb(sqrt_s);
            
            mu = (computeKt() * computeKt() * sigmatt_SM
                + computeKb() * computeKb() * sigmabb_SM
                + computeKt() * computeKb() * sigmatb_SM)
                / (sigmatt_SM + sigmabb_SM + sigmatb_SM);
        }
    }
    
    return (mu*(1.0 + eggFint + eggFpar));    
}

double HiggsKigen::muVBF(const double sqrt_s) const
{
    double mu = 1.0;
    double sigmaWF_SM = trueSM.computeSigmaWF(sqrt_s);
    double sigmaZF_SM = trueSM.computeSigmaZF(sqrt_s);
    double sigmaZWF_SM = trueSM.computeSigmaZWF(sqrt_s);
    mu = (computeKW() * computeKW() * sigmaWF_SM
            + computeKZ() * computeKZ() * sigmaZF_SM
            + computeKW() * computeKZ() * sigmaZWF_SM)
            / (sigmaWF_SM + sigmaZF_SM + sigmaZWF_SM);
    return (mu*(1.0 + eVBFint + eVBFpar));
}

double HiggsKigen::muVBFgamma(const double sqrt_s) const
{
    // Assume photon comes from the internal W lines -> no ZBF
    double mu = 1.0;    
    mu = (computeKW() * computeKW());
    return (mu*(1.0 + eVBFint + eVBFpar));
}

double HiggsKigen::mueeWBF(const double sqrt_s) const
{
    double mu = 1.0;
    mu = (computeKW() * computeKW());
    return (mu*(1.0 + eeeWBFint + eeeWBFpar));
}

double HiggsKigen::mueeWBFPol(const double sqrt_s, const double Pol_em, const double Pol_ep) const
{
    double mu = 1.0;
    mu = (computeKW() * computeKW());
    return (mu*(1.0 + eeeWBFint + eeeWBFpar));
}

double HiggsKigen::mueeHvv(const double sqrt_s) const
{
//  In this model we neglect the contributions from HZ->vv with off-shell Z
//  Dependence is exact if Custodial symmetry option is on
    double mu = 1.0;
    mu = (computeKW() * computeKW());
    return (mu*(1.0 + eeeWBFint + eeeWBFpar));
}

double HiggsKigen::mueeHvvPol(const double sqrt_s, const double Pol_em, const double Pol_ep) const
{
//  In this model we neglect the contributions from HZ->vv with off-shell Z
//  Dependence is exact if Custodial symmetry option is on
    double mu = 1.0;
    mu = (computeKW() * computeKW());
    return (mu*(1.0 + eeeWBFint + eeeWBFpar));
}

double HiggsKigen::mueeZBF(const double sqrt_s) const
{
    double mu = 1.0;
    mu = (computeKZ() * computeKZ());
    return (mu*(1.0 + eeeWBFint + eeeWBFpar));
}

double HiggsKigen::mueeZBFPol(const double sqrt_s, const double Pol_em, const double Pol_ep) const
{
    double mu = 1.0;
    mu = (computeKZ() * computeKZ());
    return (mu*(1.0 + eeeWBFint + eeeWBFpar));
}

double HiggsKigen::muepWBF(const double sqrt_s) const
{
    double mu = 1.0;
    mu = (computeKW() * computeKW());    
    return (mu*(1.0 + eepWBFint + eepWBFpar));
}

double HiggsKigen::muepZBF(const double sqrt_s) const
{
    double mu = 1.0;
    mu = (computeKZ() * computeKZ());
    return (mu*(1.0 + eepZBFint + eepZBFpar));
}

double HiggsKigen::muWH(const double sqrt_s) const
{
    double mu = 1.0;
    mu = (computeKW() * computeKW());
    return (mu*(1.0 + eWHint + eWHpar));
}

double HiggsKigen::muWHpT250(const double sqrt_s) const
{
    return muWH(sqrt_s);
}

double HiggsKigen::muZH(const double sqrt_s) const
{
    double mu = 1.0;
    mu = (computeKZ() * computeKZ());
    return (mu*(1.0 + eZHint + eZHpar));
}

double HiggsKigen::muZHpT250(const double sqrt_s) const
{
    return muZH(sqrt_s);
}

double HiggsKigen::mueeZH(const double sqrt_s) const
{
    double mu = 1.0;
    mu = (computeKZ() * computeKZ());
    return (mu*(1.0 + eeeZHint + eeeZHpar));
}

double HiggsKigen::mueeZllH(const double sqrt_s) const
{
    double mu = 1.0;
    mu = (computeKZ() * computeKZ());
    return (mu*(1.0 + eeeZHint + eeeZHpar));
}

double HiggsKigen::mueeZqqH(const double sqrt_s) const
{
    double mu = 1.0;
    mu = (computeKZ() * computeKZ());
    return (mu*(1.0 + eeeZHint + eeeZHpar));
}

double HiggsKigen::mueeZHPol(const double sqrt_s, const double Pol_em, const double Pol_ep) const
{
    double mu = 1.0;
    mu = (computeKZ() * computeKZ());
    return (mu*(1.0 + eeeZHint + eeeZHpar));
}

double HiggsKigen::mueeZllHPol(const double sqrt_s, const double Pol_em, const double Pol_ep) const
{
    double mu = 1.0;
    mu = (computeKZ() * computeKZ());
    return (mu*(1.0 + eeeZHint + eeeZHpar));
}

double HiggsKigen::mueeZqqHPol(const double sqrt_s, const double Pol_em, const double Pol_ep) const
{
    double mu = 1.0;
    mu = (computeKZ() * computeKZ());
    return (mu*(1.0 + eeeZHint + eeeZHpar));
}

double HiggsKigen::muVH(const double sqrt_s) const
{
    double sigmaWH_SM = trueSM.computeSigmaWH(sqrt_s);
    double sigmaZH_SM = trueSM.computeSigmaZH(sqrt_s);
    return ((computeKW() * computeKW()*(1.0 + eWHint + eWHpar) * sigmaWH_SM
            + computeKZ() * computeKZ()*(1.0 + eZHint + eZHpar) * sigmaZH_SM)
            / (sigmaWH_SM + sigmaZH_SM));
}

double HiggsKigen::muVHpT250(const double sqrt_s) const
{
    //Use MG SM values
    double sigmaWH_SM = 0.26944e-01;
    double sigmaZH_SM = 0.14600e-01;
    return ((computeKW() * computeKW()*(1.0 + eWHint + eWHpar) * sigmaWH_SM
            + computeKZ() * computeKZ()*(1.0 + eZHint + eZHpar) * sigmaZH_SM)
            / (sigmaWH_SM + sigmaZH_SM));
}

double HiggsKigen::muVBFpVH(const double sqrt_s) const
{
    double sigmaWH_SM = trueSM.computeSigmaWH(sqrt_s);
    double sigmaZH_SM = trueSM.computeSigmaZH(sqrt_s);
    double sigmaWF_SM = trueSM.computeSigmaWF(sqrt_s);
    double sigmaZF_SM = trueSM.computeSigmaZF(sqrt_s);
    double sigmaZWF_SM = trueSM.computeSigmaZWF(sqrt_s);    
    double sigmaVBF_SM = sigmaWF_SM + sigmaZF_SM + sigmaZWF_SM;
    
    double sigmaWH = muWH(sqrt_s) * sigmaWH_SM;
    double sigmaZH = muZH(sqrt_s) * sigmaZH_SM;
    double sigmaVBF = muVBF(sqrt_s) * sigmaVBF_SM;

    return ((sigmaWH + sigmaZH + sigmaVBF) / (sigmaWH_SM + sigmaZH_SM + sigmaVBF_SM));
}

double HiggsKigen::muttH(const double sqrt_s) const
{
    double mu = 1.0;
    mu = (computeKt() * computeKt());
    return (mu*(1.0 + ettHint + ettHpar));
}

double HiggsKigen::mutHq(const double sqrt_s) const
{
    double mu = 1.0;
    double Kt = computeKt();
    double Kw = computeKW();
    
//  From ATLAS parametrization of the K-framework (at 14 TeV?)   
    mu = 2.63 * Kt * Kt + 3.58 * Kw * Kw - 5.21 * Kt * Kw;
    
    return mu;
}

double HiggsKigen::muggHpttH(const double sqrt_s) const
{
    double sigmaggH_SM = trueSM.computeSigmaggH(sqrt_s);
    double sigmattH_SM = trueSM.computeSigmattH(sqrt_s);

    double sigmaggH = muggH(sqrt_s) * sigmaggH_SM;
    double sigmattH = muttH(sqrt_s) * sigmattH_SM;

    return ((sigmaggH + sigmattH) / (sigmaggH_SM + sigmattH_SM));
} 

double HiggsKigen::mueettH(const double sqrt_s) const
{
    double mu = 1.0;
    mu = (computeKt() * computeKt());
    return (mu*(1.0 + eeettHint + eeettHpar));
}

double HiggsKigen::mummH(const double sqrt_s) const
{
    double mu = 1.0;
//    mu = (computeKmu() * computeKmu());
    
    if (sqrt_s == 0.125) {
        
//  Peak production cross section mu mu -> H -> X = 4 pi/mH^2 * BR(H->mu mu) * BR(H-> X)
//  Use mu mu -> H = 4 pi/mH^2 * BR(H->mu mu), so the xs BR formulae still applies
        mu = BrHmumuRatio();
        
    } else
        throw std::runtime_error("Bad argument in HiggsKigen::mummH()");  
    
    return mu;
}

double HiggsKigen::mummHNWA(const double sqrt_s) const
{
    double mu = 1.0;
    mu = (computeKmu() * computeKmu());
    
    return mu;
}

double HiggsKigen::mueettHPol(const double sqrt_s, const double Pol_em, const double Pol_ep) const
{
    double mu = 1.0;
    mu = (computeKt() * computeKt());
    return (mu*(1.0 + eeettHint + eeettHpar));
}

double HiggsKigen::BrHggRatio() const
{
    double Br = 1.0;
    Br = (computeKg() * computeKg() / GammaHTotR);
    return (Br * (1.0 + eHggint + eHggpar));
}

double HiggsKigen::BrHWWRatio() const
{
    double Br = 1.0;
    Br = (computeKW() * computeKW() / GammaHTotR);
    return (Br * (1.0 + eHWWint + eHWWpar));
}

double HiggsKigen::BrHWW2l2vRatio() const
{
    return BrHWWRatio();
}

double HiggsKigen::BrHZZRatio() const
{
    double Br = 1.0;
    Br = (computeKZ() * computeKZ() / GammaHTotR);
    return (Br * (1.0 + eHZZint + eHZZpar));
}

double HiggsKigen::BrHZZ4lRatio() const
{
    return BrHZZRatio();
}

double HiggsKigen::BrHZZ4eRatio() const
{
    return BrHZZRatio();
}

double HiggsKigen::BrHZZ2e2muRatio() const
{
    return BrHZZRatio();
}

double HiggsKigen::BrHZZ4muRatio() const
{
    return BrHZZRatio();
}

double HiggsKigen::BrHZgaRatio() const
{
    double Br = 1.0;
    Br = (computeKZga() * computeKZga() / GammaHTotR);
    return (Br * (1.0 + eHZgaint + eHZgapar));
}

double HiggsKigen::BrHZgallRatio() const
{
    return BrHZgaRatio();
}

double HiggsKigen::BrHZgaeeRatio() const
{
    return BrHZgaRatio();
}

double HiggsKigen::BrHZgamumuRatio() const
{
    return BrHZgaRatio();
}

double HiggsKigen::BrHgagaRatio() const
{
    double Br = 1.0;
    Br = (computeKgaga() * computeKgaga() / GammaHTotR);
    return (Br * (1.0 + eHgagaint + eHgagapar));
}

double HiggsKigen::BrHmumuRatio() const
{
    double Br = 1.0;
    Br = (computeKmu() * computeKmu() / GammaHTotR);
    return (Br * (1.0 + eHmumuint + eHmumupar));
}

double HiggsKigen::BrHtautauRatio() const
{
    double Br = 1.0;
    Br = (computeKtau() * computeKtau() / GammaHTotR);
    return (Br * (1.0 + eHtautauint + eHtautaupar));
}

double HiggsKigen::BrHccRatio() const
{
    double Br = 1.0;
    Br = (computeKc() * computeKc() / GammaHTotR);
    return (Br * (1.0 + eHccint + eHccpar));
}

double HiggsKigen::BrHbbRatio() const
{
    double Br = 1.0;
    Br = (computeKb() * computeKb() / GammaHTotR);
    return (Br * (1.0 + eHbbint + eHbbpar));
}

double HiggsKigen::muttHZbbboost(const double sqrt_s) const
{    
    return computeKt() * computeKt() * computeKb() * computeKb() / GammaHTotR;
    
}


double HiggsKigen::muggHgaga(const double sqrt_s) const
{
    return muggH(sqrt_s)*BrHgagaRatio();
}
    
double HiggsKigen::muggHgagaInt(const double sqrt_s) const
{
    double muNWA, GammaRatio, fki;
    double kt,kb,kc,ks,ku,kd,kta,km,ke;
    double kw;
    double kg,kgamma;
    
    muNWA = muggH(sqrt_s)*BrHgagaRatio();
    
    if (!FlagKiLoop) {
    
        kt = computeKt();
        kb = computeKb();
        kc = computeKc();
        ks = computeKs();
        ku = computeKu();
        kd = computeKd();
        kta = computeKtau();
        km = computeKmu();
        ke = computeKe();
    
        kw = computeKW();
    
        kg = 0.0;
        kgamma = 0.0;
    
        GammaRatio = computeGammaTotalRatio();
    
    
        fki = 0.000802422 *kb*kb + 0.000312884 *kb*kc + 0.0000182107 *kc*kc +
            5.94769e-9 *kb*kd + 9.62554e-10 *kc*kd + 6.2785e-15 *kd*kd + 
            5.53251e-10 *kb*ke + 3.51863e-11 *kc*ke + 1.09243e-15 *kd*ke - 
            0.00905016 *kb*kg - 0.00190706 *kc*kg - 5.97591e-9 *kd*kg - 
            6.72288e-10 *ke*kg - 0.0271505 *kb*kgamma - 0.00143029 *kc*kgamma - 
            1.79277e-8 *kd*kgamma - 0.174392 *kg*kgamma + 8.97565e-6 *kb*km + 
            6.21013e-7 *kc*km + 2.37616e-11 *kd*km - 0.0000460022 *kg*km + 
            4.93348e-6 *kb*ks + 8.51176e-7 *kc*ks + 1.29558e-11 *kd*ks + 
            1.16267e-12 *ke*ks - 0.0000123381 *kg*ks - 0.0000370143 *kgamma*ks + 
            2.22544e-8 *km*ks + 6.08665e-9 *ks*ks - 0.0467672 *kb*kt - 
            0.00394193 *kc*kt - 3.08808e-8 *kd*kt - 6.94817e-10 *ke*kt - 
            0.240315 *kg*kt - 0.180236 *kgamma*kt - 0.0000475437 *km*kt - 
            0.0000637578 *ks*kt - 0.248368 *kt*kt + 0.00100168 *kb*kta + 
            0.0000759092 *kc*kta + 3.44671e-9 *kd*kta - 0.00975386 *kg*kta + 
            2.93009e-6 *ks*kta - 0.0100807 *kt*kta + 5.30126e-8 *kb*ku + 
            5.54256e-9 *kc*ku + 1.15815e-13 *kd*ku + 4.05206e-15 *ke*ku - 
            1.03323e-7 *kg*ku - 7.74926e-8 *kgamma*ku + 
            8.62762e-11 *km*ku + 1.17664e-10 *ks*ku - 2.13572e-7 *kt*ku + 
            1.2332e-8 *kta*ku + 3.40922e-13 *ku*ku + 0.169912 *kb*kw + 
            0.00895098 *kc*kw + 1.12194e-7 *kd*kw + 1.09137 *kg*kw + 
            0.000231641 *ks*kw + 1.12795 *kt*kw + 4.8496e-7 *ku*kw;
    
        return (muNWA - 0.022 *GammaRatio * fki)/0.978;
    } else {
        return muNWA;
    }
}
    
double HiggsKigen::muVBFHgaga(const double sqrt_s) const
{
    return muVBF(sqrt_s)*BrHgagaRatio();
}

double HiggsKigen::muZHgaga(const double sqrt_s) const
{
    return muZH(sqrt_s)*BrHgagaRatio();
}

double HiggsKigen::muWHgaga(const double sqrt_s) const
{
    return muWH(sqrt_s)*BrHgagaRatio();
}

double HiggsKigen::muVHgaga(const double sqrt_s) const
{
    return muVH(sqrt_s)*BrHgagaRatio();
}

double HiggsKigen::muttHgaga(const double sqrt_s) const
{
    return muttH(sqrt_s)*BrHgagaRatio();
}

double HiggsKigen::muggHZga(const double sqrt_s) const
{
    return muggH(sqrt_s)*BrHZgaRatio();
}

double HiggsKigen::muVBFHZga(const double sqrt_s) const
{
    return muVBF(sqrt_s)*BrHZgaRatio();
}

double HiggsKigen::muZHZga(const double sqrt_s) const
{
    return muZH(sqrt_s)*BrHZgaRatio();
}

double HiggsKigen::muWHZga(const double sqrt_s) const
{
    return muWH(sqrt_s)*BrHZgaRatio();
}

double HiggsKigen::muVHZga(const double sqrt_s) const
{
    return muVH(sqrt_s)*BrHZgaRatio();
}

double HiggsKigen::muttHZga(const double sqrt_s) const
{
    return muttH(sqrt_s)*BrHZgaRatio();
}

double HiggsKigen::muggHZZ(const double sqrt_s) const
{
    return muggH(sqrt_s)*BrHZZRatio();
}

double HiggsKigen::muVBFHZZ(const double sqrt_s) const
{
    return muVBF(sqrt_s)*BrHZZRatio();
}

double HiggsKigen::muZHZZ(const double sqrt_s) const
{
    return muZH(sqrt_s)*BrHZZRatio();
}

double HiggsKigen::muWHZZ(const double sqrt_s) const
{
    return muWH(sqrt_s)*BrHZZRatio();
}

double HiggsKigen::muVHZZ(const double sqrt_s) const
{
    return muVH(sqrt_s)*BrHZZRatio();
}

double HiggsKigen::muttHZZ(const double sqrt_s) const
{
    return muttH(sqrt_s)*BrHZZRatio();
}

double HiggsKigen::muggHZZ4l(const double sqrt_s) const
{
    return muggH(sqrt_s)*BrHZZ4lRatio();
}

double HiggsKigen::muVBFHZZ4l(const double sqrt_s) const
{
    return muVBF(sqrt_s)*BrHZZ4lRatio();
}

double HiggsKigen::muZHZZ4l(const double sqrt_s) const
{
    return muZH(sqrt_s)*BrHZZ4lRatio();
}

double HiggsKigen::muWHZZ4l(const double sqrt_s) const
{
    return muWH(sqrt_s)*BrHZZ4lRatio();
}

double HiggsKigen::muVHZZ4l(const double sqrt_s) const
{
    return muVH(sqrt_s)*BrHZZ4lRatio();
}

double HiggsKigen::muttHZZ4l(const double sqrt_s) const
{
    return muttH(sqrt_s)*BrHZZ4lRatio();
}

double HiggsKigen::muggHWW(const double sqrt_s) const
{
    return muggH(sqrt_s)*BrHWWRatio();
}

double HiggsKigen::muVBFHWW(const double sqrt_s) const
{
    return muVBF(sqrt_s)*BrHWWRatio();
}

double HiggsKigen::muZHWW(const double sqrt_s) const
{
    return muZH(sqrt_s)*BrHWWRatio();
}

double HiggsKigen::muWHWW(const double sqrt_s) const
{
    return muWH(sqrt_s)*BrHWWRatio();
}

double HiggsKigen::muVHWW(const double sqrt_s) const
{
    return muVH(sqrt_s)*BrHWWRatio();
}

double HiggsKigen::muttHWW(const double sqrt_s) const
{
    return muttH(sqrt_s)*BrHWWRatio();
}

double HiggsKigen::muggHWW2l2v(const double sqrt_s) const
{
    return muggH(sqrt_s)*BrHWW2l2vRatio();
}

double HiggsKigen::muVBFHWW2l2v(const double sqrt_s) const
{
    return muVBF(sqrt_s)*BrHWW2l2vRatio();
}

double HiggsKigen::muZHWW2l2v(const double sqrt_s) const
{
    return muZH(sqrt_s)*BrHWW2l2vRatio();
}

double HiggsKigen::muWHWW2l2v(const double sqrt_s) const
{
    return muWH(sqrt_s)*BrHWW2l2vRatio();
}

double HiggsKigen::muVHWW2l2v(const double sqrt_s) const
{
    return muVH(sqrt_s)*BrHWW2l2vRatio();
}

double HiggsKigen::muttHWW2l2v(const double sqrt_s) const
{
    return muttH(sqrt_s)*BrHWW2l2vRatio();
}

double HiggsKigen::muggHmumu(const double sqrt_s) const
{
    return muggH(sqrt_s)*BrHmumuRatio();
}

double HiggsKigen::muVBFHmumu(const double sqrt_s) const
{
    return muVBF(sqrt_s)*BrHmumuRatio();
}

double HiggsKigen::muZHmumu(const double sqrt_s) const
{
    return muZH(sqrt_s)*BrHmumuRatio();
}

double HiggsKigen::muWHmumu(const double sqrt_s) const
{
    return muWH(sqrt_s)*BrHmumuRatio();
}

double HiggsKigen::muVHmumu(const double sqrt_s) const
{
    return muVH(sqrt_s)*BrHmumuRatio();
}

double HiggsKigen::muttHmumu(const double sqrt_s) const
{
    return muttH(sqrt_s)*BrHmumuRatio();
}

double HiggsKigen::muggHtautau(const double sqrt_s) const
{
    return muggH(sqrt_s)*BrHtautauRatio();
}

double HiggsKigen::muVBFHtautau(const double sqrt_s) const
{
    return muVBF(sqrt_s)*BrHtautauRatio();
}

double HiggsKigen::muZHtautau(const double sqrt_s) const
{
    return muZH(sqrt_s)*BrHtautauRatio();
}

double HiggsKigen::muWHtautau(const double sqrt_s) const
{
    return muWH(sqrt_s)*BrHtautauRatio();
}

double HiggsKigen::muVHtautau(const double sqrt_s) const
{
    return muVH(sqrt_s)*BrHtautauRatio();
}

double HiggsKigen::muttHtautau(const double sqrt_s) const
{
    return muttH(sqrt_s)*BrHtautauRatio();
}

double HiggsKigen::muggHbb(const double sqrt_s) const
{
    return muggH(sqrt_s)*BrHbbRatio();
}

double HiggsKigen::muVBFHbb(const double sqrt_s) const
{
    return muVBF(sqrt_s)*BrHbbRatio();
}

double HiggsKigen::muZHbb(const double sqrt_s) const
{
    return muZH(sqrt_s)*BrHbbRatio();
}

double HiggsKigen::muWHbb(const double sqrt_s) const
{
    return muWH(sqrt_s)*BrHbbRatio();
}

double HiggsKigen::muVHbb(const double sqrt_s) const
{
    return muVH(sqrt_s)*BrHbbRatio();
}

double HiggsKigen::muttHbb(const double sqrt_s) const
{
    return muttH(sqrt_s)*BrHbbRatio();
}


////////////////////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------------------
//-- Special Hadron collider signal strengths with separate full TH unc U(prod x decay) ---
//-----------------------------------------------------------------------------------------
//////////////////////////////////////////////////////////////////////////////////////////// 

double HiggsKigen::muTHUggHgaga(const double sqrt_s) const
{
    return ( muggH(sqrt_s)*BrHgagaRatio() * (1.0 + eggFHgaga ) * (1.0 + eHwidth)/(1.0 + eggFint + eggFpar)/(1.0 + eHgagaint + eHgagapar) );
}
   
double HiggsKigen::muTHUVBFHgaga(const double sqrt_s) const
{
    return ( muVBF(sqrt_s)*BrHgagaRatio() * (1.0 + eVBFHgaga ) * (1.0 + eHwidth)/(1.0 + eVBFint + eVBFpar)/(1.0 + eHgagaint + eHgagapar) );
}

double HiggsKigen::muTHUZHgaga(const double sqrt_s) const
{
    return ( muZH(sqrt_s)*BrHgagaRatio() * (1.0 + eZHgaga ) * (1.0 + eHwidth)/(1.0 + eZHint + eZHpar)/(1.0 + eHgagaint + eHgagapar) );
}

double HiggsKigen::muTHUWHgaga(const double sqrt_s) const
{
    return ( muWH(sqrt_s)*BrHgagaRatio() * (1.0 + eWHgaga ) * (1.0 + eHwidth)/(1.0 + eWHint + eWHpar)/(1.0 + eHgagaint + eHgagapar) );
}

double HiggsKigen::muTHUVHgaga(const double sqrt_s) const
{    
    //  Theory uncertainty in VH production, from the WH and ZH ones
    double sigmaWH_SM = trueSM.computeSigmaWH(sqrt_s);
    double sigmaZH_SM = trueSM.computeSigmaZH(sqrt_s);    
    double eVHtot,eVHgaga;
    
    eVHtot = ((eWHint + eWHpar) * sigmaWH_SM + (eZHint + eZHpar) * sigmaZH_SM) / (sigmaWH_SM + sigmaZH_SM);

    eVHgaga = (eWHgaga * sigmaWH_SM + eZHgaga * sigmaZH_SM) / (sigmaWH_SM + sigmaZH_SM);
    
    return ( muVH(sqrt_s)*BrHgagaRatio() * (1.0 + eVHgaga ) * (1.0 + eHwidth)/(1.0 + eVHtot)/(1.0 + eHgagaint + eHgagapar) );
}

double HiggsKigen::muTHUttHgaga(const double sqrt_s) const
{
    return ( muttH(sqrt_s)*BrHgagaRatio() * (1.0 + ettHgaga ) * (1.0 + eHwidth)/(1.0 + eeettHint + eeettHpar)/(1.0 + eHgagaint + eHgagapar) );
}

double HiggsKigen::muTHUggHZga(const double sqrt_s) const
{
    return ( muggH(sqrt_s)*BrHZgaRatio() * (1.0 + eggFHZga ) * (1.0 + eHwidth)/(1.0 + eggFint + eggFpar)/(1.0 + eHZgaint + eHZgapar) );
}

double HiggsKigen::muTHUVBFHZga(const double sqrt_s) const
{
    return ( muVBF(sqrt_s)*BrHZgaRatio() * (1.0 + eVBFHZga ) * (1.0 + eHwidth)/(1.0 + eVBFint + eVBFpar)/(1.0 + eHZgaint + eHZgapar) );
}

double HiggsKigen::muTHUZHZga(const double sqrt_s) const
{
    return ( muZH(sqrt_s)*BrHZgaRatio() * (1.0 + eZHZga ) * (1.0 + eHwidth)/(1.0 + eZHint + eZHpar)/(1.0 + eHZgaint + eHZgapar) );
}

double HiggsKigen::muTHUWHZga(const double sqrt_s) const
{
    return ( muWH(sqrt_s)*BrHZgaRatio() * (1.0 + eWHZga ) * (1.0 + eHwidth)/(1.0 + eWHint + eWHpar)/(1.0 + eHZgaint + eHZgapar) );
}

double HiggsKigen::muTHUVHZga(const double sqrt_s) const
{
    //  Theory uncertainty in VH production, from the WH and ZH ones
    double sigmaWH_SM = trueSM.computeSigmaWH(sqrt_s);
    double sigmaZH_SM = trueSM.computeSigmaZH(sqrt_s);    
    double eVHtot,eVHZga;
    
    eVHtot = ((eWHint + eWHpar) * sigmaWH_SM + (eZHint + eZHpar) * sigmaZH_SM) / (sigmaWH_SM + sigmaZH_SM);

    eVHZga = (eWHZga * sigmaWH_SM + eZHZga * sigmaZH_SM) / (sigmaWH_SM + sigmaZH_SM);
    
    return ( muVH(sqrt_s)*BrHZgaRatio() * (1.0 + eVHZga ) * (1.0 + eHwidth)/(1.0 + eVHtot)/(1.0 + eHZgaint + eHZgapar) );
}

double HiggsKigen::muTHUttHZga(const double sqrt_s) const
{
    return ( muttH(sqrt_s)*BrHZgaRatio() * (1.0 + ettHZga ) * (1.0 + eHwidth)/(1.0 + eeettHint + eeettHpar)/(1.0 + eHZgaint + eHZgapar) );
}

double HiggsKigen::muTHUggHZZ(const double sqrt_s) const
{
    return ( muggH(sqrt_s)*BrHZZRatio() * (1.0 + eggFHZZ ) * (1.0 + eHwidth)/(1.0 + eggFint + eggFpar)/(1.0 + eHZZint + eHZZpar) );
}

double HiggsKigen::muTHUVBFHZZ(const double sqrt_s) const
{
    return ( muVBF(sqrt_s)*BrHZZRatio() * (1.0 + eVBFHZZ ) * (1.0 + eHwidth)/(1.0 + eVBFint + eVBFpar)/(1.0 + eHZZint + eHZZpar) );
}

double HiggsKigen::muTHUZHZZ(const double sqrt_s) const
{
    return ( muZH(sqrt_s)*BrHZZRatio() * (1.0 + eZHZZ ) * (1.0 + eHwidth)/(1.0 + eZHint + eZHpar)/(1.0 + eHZZint + eHZZpar) );
}

double HiggsKigen::muTHUWHZZ(const double sqrt_s) const
{
    return ( muWH(sqrt_s)*BrHZZRatio() * (1.0 + eWHZZ ) * (1.0 + eHwidth)/(1.0 + eWHint + eWHpar)/(1.0 + eHZZint + eHZZpar) );
}

double HiggsKigen::muTHUVHZZ(const double sqrt_s) const
{
    //  Theory uncertainty in VH production, from the WH and ZH ones
    double sigmaWH_SM = trueSM.computeSigmaWH(sqrt_s);
    double sigmaZH_SM = trueSM.computeSigmaZH(sqrt_s);    
    double eVHtot,eVHZZ;
    
    eVHtot = ((eWHint + eWHpar) * sigmaWH_SM + (eZHint + eZHpar) * sigmaZH_SM) / (sigmaWH_SM + sigmaZH_SM);

    eVHZZ = (eWHZZ * sigmaWH_SM + eZHZZ * sigmaZH_SM) / (sigmaWH_SM + sigmaZH_SM);
    
    return ( muVH(sqrt_s)*BrHZZRatio() * (1.0 + eVHZZ ) * (1.0 + eHwidth)/(1.0 + eVHtot)/(1.0 + eHZZint + eHZZpar) );
}

double HiggsKigen::muTHUttHZZ(const double sqrt_s) const
{
    return ( muttH(sqrt_s)*BrHZZRatio() * (1.0 + ettHZZ ) * (1.0 + eHwidth)/(1.0 + eeettHint + eeettHpar)/(1.0 + eHZZint + eHZZpar) );
}

double HiggsKigen::muTHUggHZZ4l(const double sqrt_s) const
{
    return ( muggH(sqrt_s)*BrHZZ4lRatio() * (1.0 + eggFHZZ ) * (1.0 + eHwidth)/(1.0 + eggFint + eggFpar)/(1.0 + eHZZint + eHZZpar) );
}

double HiggsKigen::muTHUVBFHZZ4l(const double sqrt_s) const
{
    return ( muVBF(sqrt_s)*BrHZZ4lRatio() * (1.0 + eVBFHZZ ) * (1.0 + eHwidth)/(1.0 + eVBFint + eVBFpar)/(1.0 + eHZZint + eHZZpar) );
}

double HiggsKigen::muTHUZHZZ4l(const double sqrt_s) const
{
    return ( muZH(sqrt_s)*BrHZZ4lRatio() * (1.0 + eZHZZ ) * (1.0 + eHwidth)/(1.0 + eZHint + eZHpar)/(1.0 + eHZZint + eHZZpar) );
}

double HiggsKigen::muTHUWHZZ4l(const double sqrt_s) const
{
    return ( muWH(sqrt_s)*BrHZZ4lRatio() * (1.0 + eWHZZ ) * (1.0 + eHwidth)/(1.0 + eWHint + eWHpar)/(1.0 + eHZZint + eHZZpar) );
}

double HiggsKigen::muTHUVHZZ4l(const double sqrt_s) const
{
    //  Theory uncertainty in VH production, from the WH and ZH ones
    double sigmaWH_SM = trueSM.computeSigmaWH(sqrt_s);
    double sigmaZH_SM = trueSM.computeSigmaZH(sqrt_s);    
    double eVHtot,eVHZZ;
    
    eVHtot = ((eWHint + eWHpar) * sigmaWH_SM + (eZHint + eZHpar) * sigmaZH_SM) / (sigmaWH_SM + sigmaZH_SM);

    eVHZZ = (eWHZZ * sigmaWH_SM + eZHZZ * sigmaZH_SM) / (sigmaWH_SM + sigmaZH_SM);
    
    return ( muVH(sqrt_s)*BrHZZ4lRatio() * (1.0 + eVHZZ ) * (1.0 + eHwidth)/(1.0 + eVHtot)/(1.0 + eHZZint + eHZZpar) );
}

double HiggsKigen::muTHUttHZZ4l(const double sqrt_s) const
{
    return ( muttH(sqrt_s)*BrHZZ4lRatio() * (1.0 + ettHZZ ) * (1.0 + eHwidth)/(1.0 + eeettHint + eeettHpar)/(1.0 + eHZZint + eHZZpar) );
}

double HiggsKigen::muTHUggHWW(const double sqrt_s) const
{
    return ( muggH(sqrt_s)*BrHWWRatio() * (1.0 + eggFHWW ) * (1.0 + eHwidth)/(1.0 + eggFint + eggFpar)/(1.0 + eHWWint + eHWWpar) );
}

double HiggsKigen::muTHUVBFHWW(const double sqrt_s) const
{
    return ( muVBF(sqrt_s)*BrHWWRatio() * (1.0 + eVBFHWW ) * (1.0 + eHwidth)/(1.0 + eVBFint + eVBFpar)/(1.0 + eHWWint + eHWWpar) );
}

double HiggsKigen::muTHUZHWW(const double sqrt_s) const
{
    return ( muZH(sqrt_s)*BrHWWRatio() * (1.0 + eZHWW ) * (1.0 + eHwidth)/(1.0 + eZHint + eZHpar)/(1.0 + eHWWint + eHWWpar) );
}

double HiggsKigen::muTHUWHWW(const double sqrt_s) const
{
    return ( muWH(sqrt_s)*BrHWWRatio() * (1.0 + eWHWW ) * (1.0 + eHwidth)/(1.0 + eWHint + eWHpar)/(1.0 + eHWWint + eHWWpar) );
}

double HiggsKigen::muTHUVHWW(const double sqrt_s) const
{
    //  Theory uncertainty in VH production, from the WH and ZH ones
    double sigmaWH_SM = trueSM.computeSigmaWH(sqrt_s);
    double sigmaZH_SM = trueSM.computeSigmaZH(sqrt_s);    
    double eVHtot,eVHWW;
    
    eVHtot = ((eWHint + eWHpar) * sigmaWH_SM + (eZHint + eZHpar) * sigmaZH_SM) / (sigmaWH_SM + sigmaZH_SM);

    eVHWW = (eWHWW * sigmaWH_SM + eZHWW * sigmaZH_SM) / (sigmaWH_SM + sigmaZH_SM);
    
    return ( muVH(sqrt_s)*BrHWWRatio() * (1.0 + eVHWW ) * (1.0 + eHwidth)/(1.0 + eVHtot)/(1.0 + eHWWint + eHWWpar) );
}

double HiggsKigen::muTHUttHWW(const double sqrt_s) const
{
    return ( muttH(sqrt_s)*BrHWWRatio() * (1.0 + ettHWW ) * (1.0 + eHwidth)/(1.0 + eeettHint + eeettHpar)/(1.0 + eHWWint + eHWWpar) );
}

double HiggsKigen::muTHUggHWW2l2v(const double sqrt_s) const
{
    return ( muggH(sqrt_s)*BrHWW2l2vRatio() * (1.0 + eggFHWW ) * (1.0 + eHwidth)/(1.0 + eggFint + eggFpar)/(1.0 + eHWWint + eHWWpar) );
}

double HiggsKigen::muTHUVBFHWW2l2v(const double sqrt_s) const
{
    return ( muVBF(sqrt_s)*BrHWW2l2vRatio() * (1.0 + eVBFHWW ) * (1.0 + eHwidth)/(1.0 + eVBFint + eVBFpar)/(1.0 + eHWWint + eHWWpar) );
}

double HiggsKigen::muTHUZHWW2l2v(const double sqrt_s) const
{
    return ( muZH(sqrt_s)*BrHWW2l2vRatio() * (1.0 + eZHWW ) * (1.0 + eHwidth)/(1.0 + eZHint + eZHpar)/(1.0 + eHWWint + eHWWpar) );
}

double HiggsKigen::muTHUWHWW2l2v(const double sqrt_s) const
{
    return ( muWH(sqrt_s)*BrHWW2l2vRatio() * (1.0 + eWHWW ) * (1.0 + eHwidth)/(1.0 + eWHint + eWHpar)/(1.0 + eHWWint + eHWWpar) );
}

double HiggsKigen::muTHUVHWW2l2v(const double sqrt_s) const
{
        //  Theory uncertainty in VH production, from the WH and ZH ones
    double sigmaWH_SM = trueSM.computeSigmaWH(sqrt_s);
    double sigmaZH_SM = trueSM.computeSigmaZH(sqrt_s);    
    double eVHtot,eVHWW;
    
    eVHtot = ((eWHint + eWHpar) * sigmaWH_SM + (eZHint + eZHpar) * sigmaZH_SM) / (sigmaWH_SM + sigmaZH_SM);

    eVHWW = (eWHWW * sigmaWH_SM + eZHWW * sigmaZH_SM) / (sigmaWH_SM + sigmaZH_SM);
    
    return ( muVH(sqrt_s)*BrHWW2l2vRatio() * (1.0 + eVHWW ) * (1.0 + eHwidth)/(1.0 + eVHtot)/(1.0 + eHWWint + eHWWpar) );
}

double HiggsKigen::muTHUttHWW2l2v(const double sqrt_s) const
{
    return ( muttH(sqrt_s)*BrHWW2l2vRatio() * (1.0 + ettHWW ) * (1.0 + eHwidth)/(1.0 + eeettHint + eeettHpar)/(1.0 + eHWWint + eHWWpar) );
}

double HiggsKigen::muTHUggHmumu(const double sqrt_s) const
{
    return ( muggH(sqrt_s)*BrHmumuRatio() * (1.0 + eggFHmumu ) * (1.0 + eHwidth)/(1.0 + eggFint + eggFpar)/(1.0 + eHmumuint + eHmumupar) );
}

double HiggsKigen::muTHUVBFHmumu(const double sqrt_s) const
{
    return ( muVBF(sqrt_s)*BrHmumuRatio() * (1.0 + eVBFHmumu ) * (1.0 + eHwidth)/(1.0 + eVBFint + eVBFpar)/(1.0 + eHmumuint + eHmumupar) );
}

double HiggsKigen::muTHUZHmumu(const double sqrt_s) const
{
    return ( muZH(sqrt_s)*BrHmumuRatio() * (1.0 + eZHmumu ) * (1.0 + eHwidth)/(1.0 + eZHint + eZHpar)/(1.0 + eHmumuint + eHmumupar) );
}

double HiggsKigen::muTHUWHmumu(const double sqrt_s) const
{
    return ( muWH(sqrt_s)*BrHmumuRatio() * (1.0 + eWHmumu ) * (1.0 + eHwidth)/(1.0 + eWHint + eWHpar)/(1.0 + eHmumuint + eHmumupar) );
}

double HiggsKigen::muTHUVHmumu(const double sqrt_s) const
{
    //  Theory uncertainty in VH production, from the WH and ZH ones
    double sigmaWH_SM = trueSM.computeSigmaWH(sqrt_s);
    double sigmaZH_SM = trueSM.computeSigmaZH(sqrt_s);    
    double eVHtot,eVHmumu;
    
    eVHtot = ((eWHint + eWHpar) * sigmaWH_SM + (eZHint + eZHpar) * sigmaZH_SM) / (sigmaWH_SM + sigmaZH_SM);

    eVHmumu = (eWHmumu * sigmaWH_SM + eZHmumu * sigmaZH_SM) / (sigmaWH_SM + sigmaZH_SM);
    
    return ( muVH(sqrt_s)*BrHmumuRatio() * (1.0 + eVHmumu ) * (1.0 + eHwidth)/(1.0 + eVHtot)/(1.0 + eHmumuint + eHmumupar) );
}

double HiggsKigen::muTHUttHmumu(const double sqrt_s) const
{
    return ( muttH(sqrt_s)*BrHmumuRatio() * (1.0 + ettHmumu ) * (1.0 + eHwidth)/(1.0 + eeettHint + eeettHpar)/(1.0 + eHmumuint + eHmumupar) );
}

double HiggsKigen::muTHUggHtautau(const double sqrt_s) const
{
    return ( muggH(sqrt_s)*BrHtautauRatio() * (1.0 + eggFHtautau ) * (1.0 + eHwidth)/(1.0 + eggFint + eggFpar)/(1.0 + eHtautauint + eHtautaupar) );
}

double HiggsKigen::muTHUVBFHtautau(const double sqrt_s) const
{
    return ( muVBF(sqrt_s)*BrHtautauRatio() * (1.0 + eVBFHtautau ) * (1.0 + eHwidth)/(1.0 + eVBFint + eVBFpar)/(1.0 + eHtautauint + eHtautaupar) );
}

double HiggsKigen::muTHUZHtautau(const double sqrt_s) const
{
    return ( muZH(sqrt_s)*BrHtautauRatio() * (1.0 + eZHtautau ) * (1.0 + eHwidth)/(1.0 + eZHint + eZHpar)/(1.0 + eHtautauint + eHtautaupar) );
}

double HiggsKigen::muTHUWHtautau(const double sqrt_s) const
{
    return ( muWH(sqrt_s)*BrHtautauRatio() * (1.0 + eWHtautau ) * (1.0 + eHwidth)/(1.0 + eWHint + eWHpar)/(1.0 + eHtautauint + eHtautaupar) );
}

double HiggsKigen::muTHUVHtautau(const double sqrt_s) const
{
    //  Theory uncertainty in VH production, from the WH and ZH ones
    double sigmaWH_SM = trueSM.computeSigmaWH(sqrt_s);
    double sigmaZH_SM = trueSM.computeSigmaZH(sqrt_s);    
    double eVHtot,eVHtautau;
    
    eVHtot = ((eWHint + eWHpar) * sigmaWH_SM + (eZHint + eZHpar) * sigmaZH_SM) / (sigmaWH_SM + sigmaZH_SM);

    eVHtautau = (eWHtautau * sigmaWH_SM + eZHtautau * sigmaZH_SM) / (sigmaWH_SM + sigmaZH_SM);
    
    return ( muVH(sqrt_s)*BrHtautauRatio() * (1.0 + eVHtautau ) * (1.0 + eHwidth)/(1.0 + eVHtot)/(1.0 + eHtautauint + eHtautaupar) );
}

double HiggsKigen::muTHUttHtautau(const double sqrt_s) const
{
    return ( muttH(sqrt_s)*BrHtautauRatio() * (1.0 + ettHtautau ) * (1.0 + eHwidth)/(1.0 + eeettHint + eeettHpar)/(1.0 + eHtautauint + eHtautaupar) );
}

double HiggsKigen::muTHUggHbb(const double sqrt_s) const
{
    return ( muggH(sqrt_s)*BrHbbRatio() * (1.0 + eggFHbb ) * (1.0 + eHwidth)/(1.0 + eggFint + eggFpar)/(1.0 + eHbbint + eHbbpar) );
}

double HiggsKigen::muTHUVBFHbb(const double sqrt_s) const
{
    return ( muVBF(sqrt_s)*BrHbbRatio() * (1.0 + eVBFHbb ) * (1.0 + eHwidth)/(1.0 + eVBFint + eVBFpar)/(1.0 + eHbbint + eHbbpar) );
}

double HiggsKigen::muTHUZHbb(const double sqrt_s) const
{
    return ( muZH(sqrt_s)*BrHbbRatio() * (1.0 + eZHbb ) * (1.0 + eHwidth)/(1.0 + eZHint + eZHpar)/(1.0 + eHbbint + eHbbpar) );
}

double HiggsKigen::muTHUWHbb(const double sqrt_s) const
{
    return ( muWH(sqrt_s)*BrHbbRatio() * (1.0 + eWHbb ) * (1.0 + eHwidth)/(1.0 + eWHint + eWHpar)/(1.0 + eHbbint + eHbbpar) );
}

double HiggsKigen::muTHUVHbb(const double sqrt_s) const
{
    //  Theory uncertainty in VH production, from the WH and ZH ones
    double sigmaWH_SM = trueSM.computeSigmaWH(sqrt_s);
    double sigmaZH_SM = trueSM.computeSigmaZH(sqrt_s);    
    double eVHtot,eVHbb;
    
    eVHtot = ((eWHint + eWHpar) * sigmaWH_SM + (eZHint + eZHpar) * sigmaZH_SM) / (sigmaWH_SM + sigmaZH_SM);

    eVHbb = (eWHbb * sigmaWH_SM + eZHbb * sigmaZH_SM) / (sigmaWH_SM + sigmaZH_SM);
    
    return ( muVH(sqrt_s)*BrHbbRatio() * (1.0 + eVHbb ) * (1.0 + eHwidth)/(1.0 + eVHtot)/(1.0 + eHbbint + eHbbpar) );
}

double HiggsKigen::muTHUttHbb(const double sqrt_s) const
{
    return ( muttH(sqrt_s)*BrHbbRatio() * (1.0 + ettHbb ) * (1.0 + eHwidth)/(1.0 + eeettHint + eeettHpar)/(1.0 + eHbbint + eHbbpar) );
}

double HiggsKigen::muTHUVBFBRinv(const double sqrt_s) const
{
    return ( muVBF(sqrt_s)*Br_H_inv() * (1.0 + eVBFHinv ) /(1.0 + eVBFint + eVBFpar) );
}

double HiggsKigen::muTHUVBFHinv(const double sqrt_s) const
{
    return ( muVBF(sqrt_s)*BrHtoinvRatio() * (1.0 + eVBFHinv ) /(1.0 + eVBFint + eVBFpar) );
}

double HiggsKigen::muTHUVHBRinv(const double sqrt_s) const
{
    //  Theory uncertainty in VH production, from the WH and ZH ones
    double sigmaWH_SM = trueSM.computeSigmaWH(sqrt_s);
    double sigmaZH_SM = trueSM.computeSigmaZH(sqrt_s);    
    double eVHtot;
    
    eVHtot = ((eWHint + eWHpar) * sigmaWH_SM + (eZHint + eZHpar) * sigmaZH_SM) / (sigmaWH_SM + sigmaZH_SM);
    
    return ( muVH(sqrt_s)*Br_H_inv() * (1.0 + eVHinv )/(1.0 + eVHtot) );
}

double HiggsKigen::muTHUVHinv(const double sqrt_s) const
{
    //  Theory uncertainty in VH production, from the WH and ZH ones
    double sigmaWH_SM = trueSM.computeSigmaWH(sqrt_s);
    double sigmaZH_SM = trueSM.computeSigmaZH(sqrt_s);    
    double eVHtot;
    
    eVHtot = ((eWHint + eWHpar) * sigmaWH_SM + (eZHint + eZHpar) * sigmaZH_SM) / (sigmaWH_SM + sigmaZH_SM);
    
    return ( muVH(sqrt_s)*BrHtoinvRatio() * (1.0 + eVHinv )/(1.0 + eVHtot) );
}

double HiggsKigen::muTHUggHZZ4mu(const double sqrt_s) const
{
    return ( muggH(sqrt_s)*BrHZZ4muRatio() * (1.0 + eggFHZZ ) * (1.0 + eHwidth)/(1.0 + eggFint + eggFpar)/(1.0 + eHZZint + eHZZpar) );
}

double HiggsKigen::muTHUggHZgamumu(const double sqrt_s) const
{
    return ( muggH(sqrt_s)*BrHZgamumuRatio() * (1.0 + eggFHZga ) * (1.0 + eHwidth)/(1.0 + eggFint + eggFpar)/(1.0 + eHZgaint + eHZgapar) );
}

////////////////////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------------------
//-- Decays ---
//-----------------------------------------------------------------------------------------
//////////////////////////////////////////////////////////////////////////////////////////// 

double HiggsKigen::computeGammaTotalRatio() const
{
    double GammaRatioSMmodes;

    GammaRatioSMmodes = computeGammaTotalSMRatio();
        
    return (GammaRatioSMmodes/(1.0 - BrHinv - BrHexo));

}

double HiggsKigen::computeGammaTotalnoErrorRatio() const
{
    double GammaRatioSMmodes;

    GammaRatioSMmodes = computeGammaTotalSMnoErrorRatio();
        
    return (GammaRatioSMmodes/(1.0 - BrHinv - BrHexo));

}

double HiggsKigen::computeGammaTotalSMRatio() const
{
    double GammaRatioSMmodes;

    GammaRatioSMmodes = ((computeKg() * computeKg() * (1.0 + eHggint + eHggpar) * trueSM.computeBrHtogg()
            + computeKW() * computeKW() * (1.0 + eHWWint + eHWWpar) * trueSM.computeBrHtoWW()
            + computeKZ() * computeKZ() * (1.0 + eHZZint + eHZZpar) * trueSM.computeBrHtoZZ()
            + computeKZga() * computeKZga() * (1.0 + eHZgaint + eHZgapar) * trueSM.computeBrHtoZga()
            + computeKgaga() * computeKgaga() * (1.0 + eHgagaint + eHgagapar) * trueSM.computeBrHtogaga()
            + computeKmu() * computeKmu() * (1.0 + eHmumuint + eHmumupar) * trueSM.computeBrHtomumu()
            + computeKtau() * computeKtau() * (1.0 + eHtautauint + eHtautaupar) * trueSM.computeBrHtotautau()
            + computeKc() * computeKc() * (1.0 + eHccint + eHccpar) * trueSM.computeBrHtocc()
            + computeKb() * computeKb() * (1.0 + eHbbint + eHbbpar) * trueSM.computeBrHtobb())
            / (trueSM.computeBrHtogg()
            + trueSM.computeBrHtoWW()
            + trueSM.computeBrHtoZZ()
            + trueSM.computeBrHtoZga()
            + trueSM.computeBrHtogaga()
            + trueSM.computeBrHtomumu()
            + trueSM.computeBrHtotautau()
            + trueSM.computeBrHtocc()
            + trueSM.computeBrHtobb()));
        
    return GammaRatioSMmodes;

}

double HiggsKigen::computeGammaTotalSMnoErrorRatio() const
{
    double GammaRatioSMmodes;

    GammaRatioSMmodes = ((computeKg() * computeKg() * trueSM.computeBrHtogg()
            + computeKW() * computeKW() * trueSM.computeBrHtoWW()
            + computeKZ() * computeKZ() * trueSM.computeBrHtoZZ()
            + computeKZga() * computeKZga() * trueSM.computeBrHtoZga()
            + computeKgaga() * computeKgaga() * trueSM.computeBrHtogaga()
            + computeKmu() * computeKmu() * trueSM.computeBrHtomumu()
            + computeKtau() * computeKtau() * trueSM.computeBrHtotautau()
            + computeKc() * computeKc() * trueSM.computeBrHtocc()
            + computeKb() * computeKb() * trueSM.computeBrHtobb())
            / (trueSM.computeBrHtogg()
            + trueSM.computeBrHtoWW()
            + trueSM.computeBrHtoZZ()
            + trueSM.computeBrHtoZga()
            + trueSM.computeBrHtogaga()
            + trueSM.computeBrHtomumu()
            + trueSM.computeBrHtotautau()
            + trueSM.computeBrHtocc()
            + trueSM.computeBrHtobb()));
        
    return GammaRatioSMmodes;

}

double HiggsKigen::Br_H_exo() const
{
    return BrHexo;    
};


double HiggsKigen::Br_H_inv() const
{
    return ( (trueSM.computeBrHto4v())*BrHZZRatio() + BrHinv );
}

double HiggsKigen::Br_H_inv_NP() const
{
    return BrHinv;
}


double HiggsKigen::BrHvisRatio() const
{    
        return (1.0 - Br_H_inv());   
}

double HiggsKigen::BrHtoinvRatio() const
{    
        return (Br_H_inv()/(trueSM.computeBrHto4v()));   
}

////////////////////////////////////////////////////////////////////////

double HiggsKigen::computeKg() const
{
    return Kg;
}

double HiggsKigen::computeKW() const
{
    return Kw;
}

double HiggsKigen::computeKZ() const
{
    return Kz;
}

double HiggsKigen::computeKZga() const
{
    return Kzga;
}

double HiggsKigen::computeKgaga() const
{
    return Kga;
}

double HiggsKigen::computeKe() const
{
    return Ke;
}

double HiggsKigen::computeKmu() const
{
    return Kmu;
}

double HiggsKigen::computeKtau() const
{
    return Ktau;
}

double HiggsKigen::computeKu() const
{
    return Ku;
}

double HiggsKigen::computeKc() const
{
    return Kc;
}

double HiggsKigen::computeKt() const
{
    return Kt;
}

double HiggsKigen::computeKd() const
{
    return Kd;
}

double HiggsKigen::computeKs() const
{
    return Ks;
}

double HiggsKigen::computeKb() const
{
    return Kb;
}

//  USE THESE TO DEFINE KAPPAg,ETC

double HiggsKigen::computeKgLoop() const
{
    double Mt=trueSM.getQuarks(QCD::TOP).getMass();
    double Mb=trueSM.getQuarks(QCD::BOTTOM).getMass();
    double Mc=trueSM.getQuarks(QCD::CHARM).getMass();
    double TAUt=4.0*Mt*Mt/(mHl*mHl);
    double TAUb=4.0*Mb*Mb/(mHl*mHl);
    double TAUc=4.0*Mc*Mc/(mHl*mHl);
    double KgEff;

    KgEff = ( ( Kt * 0.5*TAUt*(1.0+(1.0-TAUt)*f_func(TAUt))
                     +Kb * 0.5*TAUb*(1.0+(1.0-TAUb)*f_func(TAUb)) 
                     +Kc * 0.5*TAUc*(1.0+(1.0-TAUc)*f_func(TAUc)) ) * (1.0+11.0*AlsMz/(4.0*M_PI))).abs2();

    KgEff = KgEff / ( ( 0.5*TAUt*(1.0+(1.0-TAUt)*f_func(TAUt))
                     + 0.5*TAUb*(1.0+(1.0-TAUb)*f_func(TAUb))
                     + 0.5*TAUc*(1.0+(1.0-TAUc)*f_func(TAUc)) ) * (1.0+11.0*AlsMz/(4.0*M_PI))).abs2();

    return (sqrt(KgEff));
}

double HiggsKigen::computeKZgaLoop() const
{
    double Mt=trueSM.getQuarks(QCD::TOP).getMass();
    double Mb=trueSM.getQuarks(QCD::BOTTOM).getMass();
    double Mc=trueSM.getQuarks(QCD::CHARM).getMass();
    double Mtau=trueSM.getLeptons(StandardModel::TAU).getMass();
    double Mmu=trueSM.getLeptons(StandardModel::MU).getMass();
    double MW=trueSM.Mw_tree();
    double cW2=trueSM.c02();
    double sW2=1.0-cW2;
    double TAUt=4.0*Mt*Mt/(mHl*mHl);
    double TAUb=4.0*Mb*Mb/(mHl*mHl);
    double TAUc=4.0*Mc*Mc/(mHl*mHl);
    double TAUtau=4.0*Mtau*Mtau/(mHl*mHl);
    double TAUmu=4.0*Mmu*Mmu/(mHl*mHl);
    double TAUw=4.0*MW*MW/(mHl*mHl);
    double LAMt=4.0*Mt*Mt/(Mz*Mz);
    double LAMb=4.0*Mb*Mb/(Mz*Mz);
    double LAMc=4.0*Mc*Mc/(Mz*Mz);
    double LAMtau=4.0*Mtau*Mtau/(Mz*Mz);
    double LAMmu=4.0*Mmu*Mmu/(Mz*Mz);
    double LAMw=4.0*MW*MW/(Mz*Mz);
    double KZgaEff;

    KZgaEff = ((-Kt * 4.0*(0.5-4.0/3.0*sW2)*(Int1(TAUt,LAMt)-Int2(TAUt,LAMt)) * (1.0-AlsMz/M_PI)
                      +Kb * 2.0*(-0.5+2.0/3.0*sW2)*(Int1(TAUb,LAMb)-Int2(TAUb,LAMb))
                      -Kc * 4.0*(0.5-4.0/3.0*sW2)*(Int1(TAUc,LAMc)-Int2(TAUc,LAMc))            
                      +Ktau * 2.0*(-0.5+2.0*sW2)*(Int1(TAUtau,LAMtau)-Int2(TAUtau,LAMtau))
                      +Kmu * 2.0*(-0.5+2.0*sW2)*(Int1(TAUmu,LAMmu)-Int2(TAUmu,LAMmu)) )/sqrt(sW2*cW2)
                     -Kw * sqrt(cW2/sW2)*(4.0*(3.0-sW2/cW2)*Int2(TAUw,LAMw)
                            +((1.0+2.0/TAUw)*sW2/cW2-(5.0+2.0/TAUw))*Int1(TAUw,LAMw))).abs2();
    
    KZgaEff = KZgaEff / ((-4.0*(0.5-4.0/3.0*sW2)*(Int1(TAUt,LAMt)-Int2(TAUt,LAMt)) * (1.0-AlsMz/M_PI)
                      + 2.0*(-0.5+2.0/3.0*sW2)*(Int1(TAUb,LAMb)-Int2(TAUb,LAMb))
                      - 4.0*(0.5-4.0/3.0*sW2)*(Int1(TAUc,LAMc)-Int2(TAUc,LAMc))
                      + 2.0*(-0.5+2.0*sW2)*(Int1(TAUtau,LAMtau)-Int2(TAUtau,LAMtau))
                      + 2.0*(-0.5+2.0*sW2)*(Int1(TAUmu,LAMmu)-Int2(TAUmu,LAMmu)) )/sqrt(sW2*cW2)
                      - sqrt(cW2/sW2)*(4.0*(3.0-sW2/cW2)*Int2(TAUw,LAMw)
                            +((1.0+2.0/TAUw)*sW2/cW2-(5.0+2.0/TAUw))*Int1(TAUw,LAMw))).abs2();
    
    return (sqrt(KZgaEff));
}

double HiggsKigen::computeKgagaLoop() const
{
    double Mt=trueSM.getQuarks(QCD::TOP).getMass();
    double Mb=trueSM.getQuarks(QCD::BOTTOM).getMass();
    double Mc=trueSM.getQuarks(QCD::CHARM).getMass();
    double Mtau=trueSM.getLeptons(StandardModel::TAU).getMass();
    double Mmu=trueSM.getLeptons(StandardModel::MU).getMass();
    double MW=trueSM.Mw_tree();
    double TAUt=4.0*Mt*Mt/(mHl*mHl);
    double TAUb=4.0*Mb*Mb/(mHl*mHl);
    double TAUc=4.0*Mc*Mc/(mHl*mHl);
    double TAUtau=4.0*Mtau*Mtau/(mHl*mHl);
    double TAUmu=4.0*Mmu*Mmu/(mHl*mHl);
    double TAUw=4.0*MW*MW/(mHl*mHl);
    double KgagaEff;

    KgagaEff = ( Kt * (8./3.)*TAUt*(1.+(1.-TAUt)*f_func(TAUt)) * (1.0-AlsMz/M_PI)
                     +Kb * (2./3.)*TAUb*(1.+(1.-TAUb)*f_func(TAUb))
                     +Kc * (8./3.)*TAUc*(1.+(1.-TAUc)*f_func(TAUc))
                     +Ktau * 2.0*TAUtau*(1.+(1.-TAUtau)*f_func(TAUtau))
                     +Kmu * 2.0*TAUmu*(1.+(1.-TAUmu)*f_func(TAUmu))
                     -Kw * (2.0+3.0*TAUw+3.0*TAUw*(2.0-TAUw)*f_func(TAUw)) ).abs2();
    
    KgagaEff = KgagaEff / ( (8./3.)*TAUt*(1.+(1.-TAUt)*f_func(TAUt)) * (1.0-AlsMz/M_PI)
                     + (2./3.)*TAUb*(1.+(1.-TAUb)*f_func(TAUb))
                     + (8./3.)*TAUc*(1.+(1.-TAUc)*f_func(TAUc))
                     + 2.0*TAUtau*(1.+(1.-TAUtau)*f_func(TAUtau))
                     + 2.0*TAUmu*(1.+(1.-TAUmu)*f_func(TAUmu))
                     - (2.0+3.0*TAUw+3.0*TAUw*(2.0-TAUw)*f_func(TAUw)) ).abs2();
    
    return (sqrt(KgagaEff));
}

////////////////////////////////////////////////////////////////////////

double HiggsKigen::deltaG_hgg() const
{
    // Not needed/implemented yet
    return 0.0;
}

double HiggsKigen::deltaG_hggRatio() const{
    return (computeKg() - 1.0);
}

double HiggsKigen::deltaG1_hWW() const{
    return 0.0;
}

double HiggsKigen::deltaG2_hWW() const{
    return 0.0;
}

double HiggsKigen::deltaG3_hWW() const{
    
    double gSM = 2.0 * (trueSM.Mw_tree())* (trueSM.Mw_tree()) / (trueSM.v());
    
    return ( gSM * (computeKW() - 1.0 ) );
}

double HiggsKigen::deltaG1_hZZ() const{
    return 0.0;
}

double HiggsKigen::deltaG2_hZZ() const{
    return 0.0;
}

double HiggsKigen::deltaG3_hZZ() const{
    
    double gSM = (trueSM.getMz()) * (trueSM.getMz()) / (trueSM.v());
    
    return ( gSM * (computeKZ() - 1.0 ) );
}

double HiggsKigen::deltaG1_hZA() const{
    // Not needed/implemented yet
    return 0.0;
}

double HiggsKigen::deltaG1_hZARatio() const{
    return (computeKZga() - 1.0);
}

double HiggsKigen::deltaG2_hZA() const{
    return 0.0;
}

double HiggsKigen::deltaG_hAA() const{
    // Not needed/implemented yet
    return 0.0;
}

double HiggsKigen::deltaG_hAARatio() const{
    return (computeKgaga() - 1.0);
}

gslpp::complex HiggsKigen::deltaG_hff(const Particle p) const{
    
    double gSM = -(p.getMass()) / (trueSM.v());
    
    if ( p.is("ELECTRON") ) {
        return ( gSM * (computeKe() - 1.0 ) );
    } else if ( p.is("MU") ) {
         return ( gSM * (computeKmu() - 1.0 ) );       
    } else if ( p.is("TAU") ) {
        return ( gSM * (computeKtau() - 1.0 ) );        
    } else if ( p.is("UP") ) {
        return ( gSM * (computeKu() - 1.0 ) );        
    } else if ( p.is("CHARM") ) {
        return ( gSM * (computeKc() - 1.0 ) );        
    } else if ( p.is("TOP") ) {
        return ( gSM * (computeKt() - 1.0 ) );        
    } else if ( p.is("DOWN") ) {
        return ( gSM * (computeKd() - 1.0 ) );        
    } else if ( p.is("STRANGE") ) {
        return ( gSM * (computeKs() - 1.0 ) );        
    } else if ( p.is("BOTTOM") ) {
        return ( gSM * (computeKb() - 1.0 ) );        
    } else {
        return 0.0;
    }
}
 
///////////////////////////////////////////////////////////////////////////////

double HiggsKigen::kappamueff() const
{
      return computeKmu();
}

double HiggsKigen::kappataueff() const
{
      return computeKtau();
}

double HiggsKigen::kappaceff() const
{
      return computeKc();
}

double HiggsKigen::kappabeff() const
{
      return computeKb();
}

double HiggsKigen::kappaGeff() const
{
      return computeKg();
}

double HiggsKigen::kappaZeff() const
{
      return computeKZ();
}

double HiggsKigen::kappaWeff() const
{
      return computeKW();
}

double HiggsKigen::kappaAeff() const
{
      return computeKgaga();
}

double HiggsKigen::kappaZAeff() const
{
      return computeKZga();
}

//////////////////////////////////////////////////////////////////////

gslpp::complex HiggsKigen::f_func(const double x) const{
    if(x<1) {
    gslpp::complex z = -gslpp::complex::i()*M_PI;
    return -pow(log((1+sqrt(1-x))/(1-sqrt(1-x)))+z,2)/4.0;
    }
    else {
        return pow(asin(sqrt(1.0/x)),2);
    }
}

gslpp::complex HiggsKigen::g_func(const double x) const{
    if(x<1) {
    gslpp::complex z = -gslpp::complex::i()*M_PI;
    gslpp::complex gs1 = sqrt(1.0-x)*(log((1.0+sqrt(1.0-x))/(1.0-sqrt(1.0-x)))+z)/2.0;
    return gs1;
    }
    else {
        gslpp::complex gg1 = sqrt(x-1.0)*asin(sqrt(1.0/x));
        return gg1;
    }
}

gslpp::complex HiggsKigen::Int1(const double tau, const double lambda) const{
    return tau*lambda/(tau-lambda)/2.0+tau*tau*lambda*lambda/((tau-lambda)
           *(tau-lambda))/2.0*(f_func(tau)-f_func(lambda))+tau*tau*lambda/((tau-lambda)
           *(tau-lambda))*(g_func(tau)-g_func(lambda));
}

gslpp::complex HiggsKigen::Int2(const double tau, const double lambda) const{
    return -tau*lambda/(tau-lambda)/2.0*(f_func(tau)-f_func(lambda));
}


