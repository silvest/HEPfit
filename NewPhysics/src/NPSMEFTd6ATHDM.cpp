/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/cppFiles/file.cc to edit this template
 */


#include "NPSMEFTd6ATHDM.h"

const std::string NPSMEFTd6ATHDM::NPSMEFTd6ATHDMVars[NNPSMEFTd6ATHDMVars] = {
    "sqY2", "Z2", "Z3", "Z4", "Z5", "Z6", "Z7", "etau", "etad", "etae"
};




NPSMEFTd6ATHDM::NPSMEFTd6ATHDM()
: NPSMEFTd6MFV()
{
    
    setModelName("NPSMEFTd6ATHDM");
    
   
    ModelParamMap.insert(std::make_pair("sqY2", std::cref(sqY2)));
    //ModelParamMap.insert(std::make_pair("Lambda_NP", std::cref(sqY2)));
    ModelParamMap.insert(std::make_pair("sqY2", std::cref(Lambda_NP)));
    ModelParamMap.insert(std::make_pair("Z2", std::cref(Z2)));
    ModelParamMap.insert(std::make_pair("Z3", std::cref(Z3)));
    ModelParamMap.insert(std::make_pair("Z4", std::cref(Z4)));
    ModelParamMap.insert(std::make_pair("Z5", std::cref(Z5)));
    ModelParamMap.insert(std::make_pair("Z6", std::cref(Z6)));
    ModelParamMap.insert(std::make_pair("Z7", std::cref(Z7)));
    ModelParamMap.insert(std::make_pair("etau", std::cref(etau)));
    ModelParamMap.insert(std::make_pair("etad", std::cref(etad)));
    ModelParamMap.insert(std::make_pair("etae", std::cref(etae)));
    
    
}



void NPSMEFTd6ATHDM::setParameter(const std::string name, const double& value)
{

    if (name.compare("sqY2") == 0) {

        sqY2 = value;
        Lambda_NP = value;

    } else if (name.compare("Z2") == 0) {
        Z2 = value;

    } else if (name.compare("Z3") == 0) {
        Z3 = value;

    }  else if (name.compare("Z4") == 0) {
        
        Z4 = value;
                       
    }  else if (name.compare("Z5") == 0) {
        
        Z5 = value;
                       
    }  else if (name.compare("Z6") == 0) {
        
        Z6 = value;
                       
    } else if (name.compare("Z7") == 0) {
        
        Z7 = value;
        
    } else if (name.compare("etau") == 0) {
        
        etau = value;
        
    } else if (name.compare("etad") == 0) {
        etad = value;
        
    } else if (name.compare("etae") == 0){
        etae = value;
    } else
        NPSMEFTd6MFV::setParameter(name, value);
}



void NPSMEFTd6ATHDM::setNPSMEFTd6MFVParameters()
{
    
    //std::cout<<"Lambda_NP = "<<Lambda_NP<<std::endl;
    
    //double Y1 = -getSMEFTCoeffEW("mh2") //check better the sign which depends on the definition in RGEsolver
    //double Z1 = 2*getSMEFTCoeffEW("lambda")
    // We do not need to include these two parameters since they are SM parameters. If we had to use them we'd just include
    // these relations to extract them at the NP scale. However, the object getSMEFTCoeffEW takes the values from 
    // the object SMEFTEvolEW whose values are initialised and run up in NPSMEFTd6General::GenerateSMInitialConditions();
    // WE WOULD NEED TO CALL THIS FUNCTION IN OUR postupdate as is done in the postupdate of NPSMEFTd6MFV if we'd like to 
    // access these parameters at this stage. There is no need of doing so since they do not enter in the tree-level matching.

    

    CH_LNP = Z6*Z6;
    //OuH[mif1_, mif2_] -> (Z6*lamU2bar[mif2, mif1])/MH2^2
    CuH_0_LNP = Z6*etau;
    CdH_0_LNP = Z6*etad;
    CeH_0_LNP = Z6*etae;

    
    // alphaOqu1[mif1_, mif2_, mif3_, mif4_] -> -1/6*(lamU2[mif3, mif2]*lamU2bar[mif4, mif1])/MH2^2
    // The implemented term in the MFV model is Cqu1_y_LNP*YucL(mif1, mif4)*YuL(mif3, mif2) 
    // Check conventions carefully
    Cqu1_y_LNP = -1/6*etau*etau;
    Cqu8_y_LNP = -1*etau*etau;
    
    Cqd1_y_LNP = -1/6*etad*etad;
    Cqd8_y_LNP = -1*etad*etad;
    
    //alphaOquqd1[mif1_, mif2_, mif3_, mif4_] -> (lamD2bar[mif4, mif3]*lamU2bar[mif2, mif1])/MH2^2
    // The implemented term in the MFV model is Cquqd1_00_LNP*YdcL(mif3, mif4)*YuL(mif1, mif2) 
    Cquqd1_00_LNP = etad*etau;
    
    //alphaOle[mif1_, mif2_, mif3_, mif4_] -> -((lamL2[mif3, mif2] lamL2bar[mif4, mif1])/(2 MH2^2))
    Cle_y_LNP = (-0.5*etae*etae);
    
    //alphaOledq[mif1_, mif2_, mif3_, mif4_] -> (lamD2[mif3, mif4] lamL2bar[mif2, mif1])/MH2^2
    Cledq_00_LNP = (1.*etad*etae);
    
    //alphaOlequ1[mif1_, mif2_, mif3_, mif4_] -> -((lamL2bar[mif2, mif1] lamU2bar[mif4, mif3])/MH2^2)
    Clequ1_00_LNP = (-1.*etau*etae); 
    
    
}

bool NPSMEFTd6ATHDM::PostUpdate() {
    
	
    setNPSMEFTd6MFVParameters();
    
    if (!NPSMEFTd6MFV::PostUpdate()) return (false);
    
    return (true);    
}


