/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "NPSMEFTd6U3.h"

const std::string NPSMEFTd6U3::NPSMEFTd6U3Vars[NNPSMEFTd6U3Vars] = {
    "CG_LNP", "CW_LNP", "CHG_LNP", "CHW_LNP", "CHB_LNP",
    "CHWB_LNP", "CHD_LNP", "CHbox_LNP", "CH_LNP",
    "CHl1_LNP","CHl3_LNP","CHe_LNP","CHq1_LNP","CHq3_LNP","CHu_LNP","CHd_LNP","CHud_LNP",
    "CeH_LNP","CuH_LNP","CdH_LNP","CuG_LNP","CuW_LNP","CuB_LNP","CdG_LNP","CdW_LNP","CdB_LNP","CeW_LNP","CeB_LNP",
    "Cll_LNP","Clq1_LNP","Clq3_LNP","Cee_LNP","Ceu_LNP","Ced_LNP","Cle_LNP","Clu_LNP","Cld_LNP","Cqe_LNP",
    "Cledq_LNP", "Cqq1_LNP", "Cqq3_LNP", "Cuu_LNP", "Cdd_LNP", "Cud1_LNP", "Cud8_LNP", "Cqu1_LNP", "Cqu8_LNP", "Cqd1_LNP", "Cqd8_LNP",
    "Cquqd1_LNP", "Cquqd8_LNP", "Clequ1_LNP", "Clequ3_LNP"
};

NPSMEFTd6U3::NPSMEFTd6U3()
: NPSMEFTd6General()
{
    
   
    ModelParamMap.insert(std::make_pair("CG_LNP", std::cref(CG_LNP)));
    ModelParamMap.insert(std::make_pair("CW_LNP", std::cref(CW_LNP)));
    ModelParamMap.insert(std::make_pair("CHG_LNP", std::cref(CHG_LNP)));
    ModelParamMap.insert(std::make_pair("CHW_LNP", std::cref(CHW_LNP)));
    ModelParamMap.insert(std::make_pair("CHB_LNP", std::cref(CHB_LNP)));
    ModelParamMap.insert(std::make_pair("CHWB_LNP", std::cref(CHWB_LNP)));
    ModelParamMap.insert(std::make_pair("CHD_LNP", std::cref(CHD_LNP)));
    ModelParamMap.insert(std::make_pair("CHbox_LNP", std::cref(CHbox_LNP)));
    ModelParamMap.insert(std::make_pair("CH_LNP", std::cref(CH_LNP)));
    
    ModelParamMap.insert(std::make_pair("CHl1_LNP", std::cref(CHl1_LNP)));
    ModelParamMap.insert(std::make_pair("CHl3_LNP", std::cref(CHl3_LNP)));
    ModelParamMap.insert(std::make_pair("CHe_LNP", std::cref(CHe_LNP)));
    ModelParamMap.insert(std::make_pair("CHq1_LNP", std::cref(CHq1_LNP)));
    ModelParamMap.insert(std::make_pair("CHq3_LNP", std::cref(CHq3_LNP)));
    ModelParamMap.insert(std::make_pair("CHu_LNP", std::cref(CHu_LNP)));
    ModelParamMap.insert(std::make_pair("CHd_LNP", std::cref(CHd_LNP)));
    ModelParamMap.insert(std::make_pair("CHud_LNP", std::cref(CHud_LNP)));
    
    
    ModelParamMap.insert(std::make_pair("CeH_LNP", std::cref(CeH_LNP)));
    ModelParamMap.insert(std::make_pair("CuH_LNP", std::cref(CuH_LNP)));
    ModelParamMap.insert(std::make_pair("CdH_LNP", std::cref(CdH_LNP)));
    ModelParamMap.insert(std::make_pair("CuG_LNP", std::cref(CuG_LNP)));
    ModelParamMap.insert(std::make_pair("CuW_LNP", std::cref(CuW_LNP)));
    ModelParamMap.insert(std::make_pair("CuB_LNP", std::cref(CuB_LNP)));
    ModelParamMap.insert(std::make_pair("CdG_LNP", std::cref(CdG_LNP)));
    ModelParamMap.insert(std::make_pair("CdW_LNP", std::cref(CdW_LNP)));
    ModelParamMap.insert(std::make_pair("CdB_LNP", std::cref(CdB_LNP)));
    ModelParamMap.insert(std::make_pair("CeW_LNP", std::cref(CeW_LNP)));
    ModelParamMap.insert(std::make_pair("CeB_LNP", std::cref(CeB_LNP)));
    
    ModelParamMap.insert(std::make_pair("Cll_LNP", std::cref(Cll_LNP)));
    ModelParamMap.insert(std::make_pair("Clq1_LNP", std::cref(Clq1_LNP)));
    ModelParamMap.insert(std::make_pair("Clq3_LNP", std::cref(Clq3_LNP)));
    ModelParamMap.insert(std::make_pair("Cee_LNP", std::cref(Cee_LNP)));
    ModelParamMap.insert(std::make_pair("Ceu_LNP", std::cref(Ceu_LNP)));
    ModelParamMap.insert(std::make_pair("Ced_LNP", std::cref(Ced_LNP)));
    ModelParamMap.insert(std::make_pair("Cle_LNP", std::cref(Cle_LNP)));
    ModelParamMap.insert(std::make_pair("Clu_LNP", std::cref(Clu_LNP)));
    ModelParamMap.insert(std::make_pair("Cld_LNP", std::cref(Cld_LNP)));
    ModelParamMap.insert(std::make_pair("Cqe_LNP", std::cref(Cqe_LNP)));
    
    
    ModelParamMap.insert(std::make_pair("Cledq_LNP", std::cref(Cledq_LNP)));
    ModelParamMap.insert(std::make_pair("Cqq1_LNP", std::cref(Cqq1_LNP)));
    ModelParamMap.insert(std::make_pair("Cqq3_LNP", std::cref(Cqq3_LNP)));
    ModelParamMap.insert(std::make_pair("Cuu_LNP", std::cref(Cuu_LNP)));
    ModelParamMap.insert(std::make_pair("Cdd_LNP", std::cref(Cdd_LNP)));
    ModelParamMap.insert(std::make_pair("Cud1_LNP", std::cref(Cud1_LNP)));
    ModelParamMap.insert(std::make_pair("Cud8_LNP", std::cref(Cud8_LNP)));
    ModelParamMap.insert(std::make_pair("Cqu1_LNP", std::cref(Cqu1_LNP)));
    ModelParamMap.insert(std::make_pair("Cqu8_LNP", std::cref(Cqu8_LNP)));
    ModelParamMap.insert(std::make_pair("Cqd1_LNP", std::cref(Cqd1_LNP)));
    ModelParamMap.insert(std::make_pair("Cqd8_LNP", std::cref(Cqd8_LNP)));
    
    
    ModelParamMap.insert(std::make_pair("Cquqd1_LNP", std::cref(Cquqd1_LNP)));
    ModelParamMap.insert(std::make_pair("Cquqd8_LNP", std::cref(Cquqd8_LNP)));
    ModelParamMap.insert(std::make_pair("Clequ1_LNP", std::cref(Clequ1_LNP)));
    ModelParamMap.insert(std::make_pair("Clequ3_LNP", std::cref(Clequ3_LNP)));
    
}




void NPSMEFTd6U3::setParameter(const std::string name, const double& value)
{


    //Decide if we want to set the bosic WC here or in the General class
    if (name.compare("CG_LNP") == 0) {
        CG_LNP = value;
        SMEFTEvol.SetCoefficient("CG", value  * LambdaNPm2);
    } else if (name.compare("CW_LNP") == 0) {
        CW_LNP = value;
        SMEFTEvol.SetCoefficient("CW", value * LambdaNPm2);
    } else if (name.compare("CHG_LNP") == 0) {
        CHG_LNP = value;
        SMEFTEvol.SetCoefficient("CHG", value * LambdaNPm2);
    } else if (name.compare("CHW_LNP") == 0) {
        CHW_LNP = value;
        SMEFTEvol.SetCoefficient("CHW", value * LambdaNPm2);
    } else if (name.compare("CHB_LNP") == 0) {
        CHB_LNP = value;
        SMEFTEvol.SetCoefficient("CHB", value * LambdaNPm2);
    } else if (name.compare("CHWB_LNP") == 0) {
        CHWB_LNP = value;
        SMEFTEvol.SetCoefficient("CHWB", value * LambdaNPm2);
    } else if (name.compare("CHD_LNP") == 0) {
        CHD_LNP = value;
        SMEFTEvol.SetCoefficient("CHD", value * LambdaNPm2);
    } else if (name.compare("CHbox_LNP") == 0) {
        CHbox_LNP = value;
        SMEFTEvol.SetCoefficient("CHbox", value * LambdaNPm2);
    } else if (name.compare("CH_LNP") == 0) {
        CH_LNP = value;
        SMEFTEvol.SetCoefficient("CH", value * LambdaNPm2);
        
        
    } else if (name.compare("CHl1_LNP") == 0) {
        
        //std::cout<<"\033[1;33m   LambdaNPm2 =  \033[0m "<< LambdaNPm2 <<std::endl;
                
        CHl1_LNP = value;
        
        CHl1_11r_LNP = value;
        SMEFTEvol.SetCoefficient("CHl1R", value * LambdaNPm2, 0, 0);
          
        CHl1_22r_LNP = value;
        SMEFTEvol.SetCoefficient("CHl1R", value * LambdaNPm2, 1, 1);
            
        CHl1_33r_LNP = value;
        SMEFTEvol.SetCoefficient("CHl1R", value * LambdaNPm2, 2, 2);
            
    } else if (name.compare("CHl3_LNP") == 0) {
        
        CHl3_LNP = value;
        
        CHl3_11r_LNP = value;
        SMEFTEvol.SetCoefficient("CHl3R", value * LambdaNPm2, 0, 0);
          
        CHl3_22r_LNP = value;
        SMEFTEvol.SetCoefficient("CHl3R", value * LambdaNPm2, 1, 1);
            
        CHl3_33r_LNP = value;
        SMEFTEvol.SetCoefficient("CHl3R", value * LambdaNPm2, 2, 2);
            
    } else if (name.compare("CHe_LNP") == 0) {
        
        CHe_LNP = value;
        
        CHe_11r_LNP = value;
        SMEFTEvol.SetCoefficient("CHeR", value * LambdaNPm2, 0, 0);
          
        CHe_22r_LNP = value;
        SMEFTEvol.SetCoefficient("CHeR", value * LambdaNPm2, 1, 1);
            
        CHe_33r_LNP = value;
        SMEFTEvol.SetCoefficient("CHeR", value * LambdaNPm2, 2, 2);
            
    } else if (name.compare("CHq1_LNP") == 0) {
        
        CHq1_LNP = value;
        
        CHq1_11r_LNP = value;
        SMEFTEvol.SetCoefficient("CHq1R", value * LambdaNPm2, 0, 0);
          
        CHq1_22r_LNP = value;
        SMEFTEvol.SetCoefficient("CHq1R", value * LambdaNPm2, 1, 1);
            
        CHq1_33r_LNP = value;
        SMEFTEvol.SetCoefficient("CHq1R", value * LambdaNPm2, 2, 2);
            
    } else if (name.compare("CHq3_LNP") == 0) {
        
        CHq3_LNP = value;
        
        CHq3_11r_LNP = value;
        SMEFTEvol.SetCoefficient("CHq3R", value * LambdaNPm2, 0, 0);
          
        CHq3_22r_LNP = value;
        SMEFTEvol.SetCoefficient("CHq3R", value * LambdaNPm2, 1, 1);
            
        CHq3_33r_LNP = value;
        SMEFTEvol.SetCoefficient("CHq3R", value * LambdaNPm2, 2, 2);
            
    } else if (name.compare("CHu_LNP") == 0) {
        
        CHu_LNP = value;
                
        CHu_11r_LNP = value;
        SMEFTEvol.SetCoefficient("CHuR", value * LambdaNPm2, 0, 0);
          
        CHu_22r_LNP = value;
        SMEFTEvol.SetCoefficient("CHuR", value * LambdaNPm2, 1, 1);
            
        CHu_33r_LNP = value;
        SMEFTEvol.SetCoefficient("CHuR", value * LambdaNPm2, 2, 2);
            
    } else if (name.compare("CHd_LNP") == 0) {
        
        CHd_LNP = value;
                
        CHd_11r_LNP = value;
        SMEFTEvol.SetCoefficient("CHdR", value * LambdaNPm2, 0, 0);
          
        CHd_22r_LNP = value;
        SMEFTEvol.SetCoefficient("CHdR", value * LambdaNPm2, 1, 1);
            
        CHd_33r_LNP = value;
        SMEFTEvol.SetCoefficient("CHdR", value * LambdaNPm2, 2, 2);
            
    } else if (name.compare("CHud_LNP") == 0) {
        
        CHud_LNP = value;
        
        CHud_11r_LNP = value;
        SMEFTEvol.SetCoefficient("CHudR", value * LambdaNPm2, 0, 0);
          
        CHud_22r_LNP = value;
        SMEFTEvol.SetCoefficient("CHudR", value * LambdaNPm2, 1, 1);
            
        CHud_33r_LNP = value;
        SMEFTEvol.SetCoefficient("CHudR", value * LambdaNPm2, 2, 2);
            
    } else if (name.compare("CeH_LNP") == 0) {
        //// Think a bit more if we need to add the Yukawas here!
        //// Basically we're relaxing the U(3) symmetry...
        
        CeH_LNP = value;
        
        CeH_11r_LNP = value;
        SMEFTEvol.SetCoefficient("CeHR", value * LambdaNPm2, 0, 0);
          
        CeH_22r_LNP = value;
        SMEFTEvol.SetCoefficient("CeHR", value * LambdaNPm2, 1, 1);
            
        CeH_33r_LNP = value;
        SMEFTEvol.SetCoefficient("CeHR", value * LambdaNPm2, 2, 2);
            
    } else if (name.compare("CuH_LNP") == 0) {
        //// Think a bit more if we need to add the Yukawas here!
        
        CuH_LNP = value;
        
        CuH_11r_LNP = value;
        SMEFTEvol.SetCoefficient("CuHR", value * LambdaNPm2, 0, 0);
          
        CuH_22r_LNP = value;
        SMEFTEvol.SetCoefficient("CuHR", value * LambdaNPm2, 1, 1);
            
        CuH_33r_LNP = value;
        SMEFTEvol.SetCoefficient("CuHR", value * LambdaNPm2, 2, 2);
            
    } else if (name.compare("CdH_LNP") == 0) {
        //// Think a bit more if we need to add the Yukawas here!
        
        CdH_LNP = value;
        
        CdH_11r_LNP = value;
        SMEFTEvol.SetCoefficient("CdHR", value * LambdaNPm2, 0, 0);
          
        CdH_22r_LNP = value;
        SMEFTEvol.SetCoefficient("CdHR", value * LambdaNPm2, 1, 1);
            
        CdH_33r_LNP = value;
        SMEFTEvol.SetCoefficient("CdHR", value * LambdaNPm2, 2, 2);
            
    } else if (name.compare("CuG_LNP") == 0) {
        
        CuG_LNP = value;
        
        CuG_11r_LNP = value;
        SMEFTEvol.SetCoefficient("CuGR", value * LambdaNPm2, 0, 0);
          
        CuG_22r_LNP = value;
        SMEFTEvol.SetCoefficient("CuGR", value * LambdaNPm2, 1, 1);
            
        CuG_33r_LNP = value;
        SMEFTEvol.SetCoefficient("CuGR", value * LambdaNPm2, 2, 2);
            
    } else if (name.compare("CuW_LNP") == 0) {
        
        CuW_LNP = value;
        
        CuW_11r_LNP = value;
        SMEFTEvol.SetCoefficient("CuWR", value * LambdaNPm2, 0, 0);
          
        CuW_22r_LNP = value;
        SMEFTEvol.SetCoefficient("CuWR", value * LambdaNPm2, 1, 1);
            
        CuW_33r_LNP = value;
        SMEFTEvol.SetCoefficient("CuWR", value * LambdaNPm2, 2, 2);
            
    } else if (name.compare("CuB_LNP") == 0) {
        
        CuB_LNP = value;
        
        CuB_11r_LNP = value;
        SMEFTEvol.SetCoefficient("CuBR", value * LambdaNPm2, 0, 0);
          
        CuB_22r_LNP = value;
        SMEFTEvol.SetCoefficient("CuBR", value * LambdaNPm2, 1, 1);
            
        CuB_33r_LNP = value;
        SMEFTEvol.SetCoefficient("CuBR", value * LambdaNPm2, 2, 2);
            
    } else if (name.compare("CdG_LNP") == 0) {
        
        CdG_LNP = value;
                
        CdG_11r_LNP = value;
        SMEFTEvol.SetCoefficient("CdGR", value * LambdaNPm2, 0, 0);
          
        CdG_22r_LNP = value;
        SMEFTEvol.SetCoefficient("CdGR", value * LambdaNPm2, 1, 1);
            
        CdG_33r_LNP = value;
        SMEFTEvol.SetCoefficient("CdGR", value * LambdaNPm2, 2, 2);
            
    } else if (name.compare("CdW_LNP") == 0) {
        
        CdW_LNP = value;
        
        CdW_11r_LNP = value;
        SMEFTEvol.SetCoefficient("CdWR", value * LambdaNPm2, 0, 0);
          
        CdW_22r_LNP = value;
        SMEFTEvol.SetCoefficient("CdWR", value * LambdaNPm2, 1, 1);
            
        CdW_33r_LNP = value;
        SMEFTEvol.SetCoefficient("CdWR", value * LambdaNPm2, 2, 2);
            
    } else if (name.compare("CdB_LNP") == 0) {
        
        CdB_LNP = value;
        
        CdB_11r_LNP = value;
        SMEFTEvol.SetCoefficient("CdBR", value * LambdaNPm2, 0, 0);
          
        CdB_22r_LNP = value;
        SMEFTEvol.SetCoefficient("CdBR", value * LambdaNPm2, 1, 1);
            
        CdB_33r_LNP = value;
        SMEFTEvol.SetCoefficient("CdBR", value * LambdaNPm2, 2, 2);
            
    } else if (name.compare("CeW_LNP") == 0) {
        
        CeW_LNP = value;
        
        CeW_11r_LNP = value;
        SMEFTEvol.SetCoefficient("CeWR", value * LambdaNPm2, 0, 0);
          
        CeW_22r_LNP = value;
        SMEFTEvol.SetCoefficient("CeWR", value * LambdaNPm2, 1, 1);
            
        CeW_33r_LNP = value;
        SMEFTEvol.SetCoefficient("CeWR", value * LambdaNPm2, 2, 2);
            
    } else if (name.compare("CeB_LNP") == 0) {
        
        CeB_LNP = value;
        
        CeB_11r_LNP = value;
        SMEFTEvol.SetCoefficient("CeBR", value * LambdaNPm2, 0, 0);
          
        CeB_22r_LNP = value;
        SMEFTEvol.SetCoefficient("CeBR", value * LambdaNPm2, 1, 1);
            
        CeB_33r_LNP = value;
        SMEFTEvol.SetCoefficient("CeBR", value * LambdaNPm2, 2, 2);
            
    } else if (name.compare("Cll_LNP") == 0) {
        
        Cll_LNP = value; 
        
        Cll_1111r_LNP = value;
        SMEFTEvol.SetCoefficient("CllR", value * LambdaNPm2, 0, 0, 0, 0);
        
        Cll_1221r_LNP = value;
        SMEFTEvol.SetCoefficient("CllR", value * LambdaNPm2, 0, 1, 1, 0);
           
        Cll_1122r_LNP = value;
        SMEFTEvol.SetCoefficient("CllR", value * LambdaNPm2, 0, 0, 1, 1);
                
        Cll_1133r_LNP = value;
        SMEFTEvol.SetCoefficient("CllR", value * LambdaNPm2, 0, 0, 2, 2);
        
        Cll_1331r_LNP = value;
        SMEFTEvol.SetCoefficient("CllR", value * LambdaNPm2, 0, 2, 2, 0);
        
        Cll_2222r_LNP = value;
        SMEFTEvol.SetCoefficient("CllR", value * LambdaNPm2, 1, 1, 1, 1);
        
        Cll_2233r_LNP = value;
        SMEFTEvol.SetCoefficient("CllR", value * LambdaNPm2, 1, 1, 2, 2);
            
        Cll_2332r_LNP = value;
        SMEFTEvol.SetCoefficient("CllR", value * LambdaNPm2, 1, 2, 2, 1);
        

        Cll_3333r_LNP = value;
        SMEFTEvol.SetCoefficient("CllR", value * LambdaNPm2, 2, 2, 2, 2);
        
            
    } else if (name.compare("Clq1_LNP") == 0) {
        
        Clq1_LNP = value;
        
        Clq1_1111r_LNP = value;
        SMEFTEvol.SetCoefficient("Clq1R", value * LambdaNPm2, 0, 0, 0, 0);
        
        Clq1_1122r_LNP = value;
        SMEFTEvol.SetCoefficient("Clq1R", value * LambdaNPm2, 0, 0, 1, 1);
            
        Clq1_2211r_LNP = value;
        SMEFTEvol.SetCoefficient("Clq1R", value * LambdaNPm2, 1, 1, 0, 0);
            
        Clq1_1133r_LNP = value;
        SMEFTEvol.SetCoefficient("Clq1R", value * LambdaNPm2, 0, 0, 2, 2);
            
        Clq1_2222r_LNP = value;
        SMEFTEvol.SetCoefficient("Clq1R", value * LambdaNPm2, 1, 1, 1, 1);
            
        Clq1_2233r_LNP = value;
        SMEFTEvol.SetCoefficient("Clq1R", value * LambdaNPm2, 1, 1, 2, 2);
            
        Clq1_3311r_LNP = value;
        SMEFTEvol.SetCoefficient("Clq1R", value * LambdaNPm2, 2, 2, 0, 0);
            
        Clq1_3322r_LNP = value;
        SMEFTEvol.SetCoefficient("Clq1R", value * LambdaNPm2, 2, 2, 1, 1);
            
        Clq1_3333r_LNP = value;
        SMEFTEvol.SetCoefficient("Clq1R", value * LambdaNPm2, 2, 2, 2, 2);
        
            
    } else if (name.compare("Clq3_LNP") == 0) {
        
        Clq3_LNP = value;
        
        Clq3_1111r_LNP = value;
        SMEFTEvol.SetCoefficient("Clq3R", value * LambdaNPm2, 0, 0, 0, 0);
        
        Clq3_1122r_LNP = value;
        SMEFTEvol.SetCoefficient("Clq3R", value * LambdaNPm2, 0, 0, 1, 1);
            
        Clq3_2211r_LNP = value;
        SMEFTEvol.SetCoefficient("Clq3R", value * LambdaNPm2, 1, 1, 0, 0);
            
        Clq3_1133r_LNP = value;
        SMEFTEvol.SetCoefficient("Clq3R", value * LambdaNPm2, 0, 0, 2, 2);
            
        Clq3_2222r_LNP = value;
        SMEFTEvol.SetCoefficient("Clq3R", value * LambdaNPm2, 1, 1, 1, 1);
            
        Clq3_2233r_LNP = value;
        SMEFTEvol.SetCoefficient("Clq3R", value * LambdaNPm2, 1, 1, 2, 2);
            
        Clq3_3311r_LNP = value;
        SMEFTEvol.SetCoefficient("Clq3R", value * LambdaNPm2, 2, 2, 0, 0);
            
        Clq3_3322r_LNP = value;
        SMEFTEvol.SetCoefficient("Clq3R", value * LambdaNPm2, 2, 2, 1, 1);
            
        Clq3_3333r_LNP = value;
        SMEFTEvol.SetCoefficient("Clq3R", value * LambdaNPm2, 2, 2, 2, 2);
        
            
    } else if (name.compare("Cee_LNP") == 0) {
        
        Cee_LNP = value;
        
        Cee_1111r_LNP = value;
        SMEFTEvol.SetCoefficient("CeeR", value * LambdaNPm2, 0, 0, 0, 0);
        
        Cee_1122r_LNP = value;
        SMEFTEvol.SetCoefficient("CeeR", value * LambdaNPm2, 0, 0, 1, 1);
            
        Cee_1133r_LNP = value;
        SMEFTEvol.SetCoefficient("CeeR", value * LambdaNPm2, 0, 0, 2, 2);
            
        Cee_2222r_LNP = value;
        SMEFTEvol.SetCoefficient("CeeR", value * LambdaNPm2, 1, 1, 1, 1);
            
        Cee_2233r_LNP = value;
        SMEFTEvol.SetCoefficient("CeeR", value * LambdaNPm2, 1, 1, 2, 2);
            
        Cee_3333r_LNP = value;
        SMEFTEvol.SetCoefficient("CeeR", value * LambdaNPm2, 2, 2, 2, 2);
    
    } else if (name.compare("Ceu_LNP") == 0) {
        
        Ceu_LNP = value;
        
        Ceu_1111r_LNP = value;
        SMEFTEvol.SetCoefficient("CeuR", value * LambdaNPm2, 0, 0, 0, 0);
        
        Ceu_1122r_LNP = value;
        SMEFTEvol.SetCoefficient("CeuR", value * LambdaNPm2, 0, 0, 1, 1);
            
        Ceu_1133r_LNP = value;
        SMEFTEvol.SetCoefficient("CeuR", value * LambdaNPm2, 0, 0, 2, 2);
            
        Ceu_2211r_LNP = value;
        SMEFTEvol.SetCoefficient("CeuR", value * LambdaNPm2, 1, 1, 0, 0);
        
        Ceu_2222r_LNP = value;
        SMEFTEvol.SetCoefficient("CeuR", value * LambdaNPm2, 1, 1, 1, 1);
            
        Ceu_2233r_LNP = value;
        SMEFTEvol.SetCoefficient("CeuR", value * LambdaNPm2, 1, 1, 2, 2);
            
        Ceu_3311r_LNP = value;
        SMEFTEvol.SetCoefficient("CeuR", value * LambdaNPm2, 2, 2, 0, 0);
            
        Ceu_3322r_LNP = value;
        SMEFTEvol.SetCoefficient("CeuR", value * LambdaNPm2, 2, 2, 1, 1);
            
        Ceu_3333r_LNP = value;
        SMEFTEvol.SetCoefficient("CeuR", value * LambdaNPm2, 2, 2, 2, 2);
                       
    } else if (name.compare("Ced_LNP") == 0) {
        
        Ced_LNP = value;
        
        Ced_1111r_LNP = value;
        SMEFTEvol.SetCoefficient("CedR", value * LambdaNPm2, 0, 0, 0, 0);
        
        Ced_1122r_LNP = value;
        SMEFTEvol.SetCoefficient("CedR", value * LambdaNPm2, 0, 0, 1, 1);
            
        Ced_1133r_LNP = value;
        SMEFTEvol.SetCoefficient("CedR", value * LambdaNPm2, 0, 0, 2, 2);
            
        Ced_2211r_LNP = value;
        SMEFTEvol.SetCoefficient("CedR", value * LambdaNPm2, 1, 1, 0, 0);
        
        Ced_2222r_LNP = value;
        SMEFTEvol.SetCoefficient("CedR", value * LambdaNPm2, 1, 1, 1, 1);
            
        Ced_2233r_LNP = value;
        SMEFTEvol.SetCoefficient("CedR", value * LambdaNPm2, 1, 1, 2, 2);
            
        Ced_3311r_LNP = value;
        SMEFTEvol.SetCoefficient("CedR", value * LambdaNPm2, 2, 2, 0, 0);
            
        Ced_3322r_LNP = value;
        SMEFTEvol.SetCoefficient("CedR", value * LambdaNPm2, 2, 2, 1, 1);
            
        Ced_3333r_LNP = value;
        SMEFTEvol.SetCoefficient("CedR", value * LambdaNPm2, 2, 2, 2, 2);
                       
    }  else if (name.compare("Cle_LNP") == 0) {
        
        Cle_LNP = value;
        
        Cle_1111r_LNP = value;
        SMEFTEvol.SetCoefficient("CleR", value * LambdaNPm2, 0, 0, 0, 0);
        
        Cle_1122r_LNP = value;
        SMEFTEvol.SetCoefficient("CleR", value * LambdaNPm2, 0, 0, 1, 1);
            
        Cle_1133r_LNP = value;
        SMEFTEvol.SetCoefficient("CleR", value * LambdaNPm2, 0, 0, 2, 2);
            
        Cle_2211r_LNP = value;
        SMEFTEvol.SetCoefficient("CleR", value * LambdaNPm2, 1, 1, 0, 0);
        
        Cle_2222r_LNP = value;
        SMEFTEvol.SetCoefficient("CleR", value * LambdaNPm2, 1, 1, 1, 1);
            
        Cle_2233r_LNP = value;
        SMEFTEvol.SetCoefficient("CleR", value * LambdaNPm2, 1, 1, 2, 2);
            
        Cle_3311r_LNP = value;
        SMEFTEvol.SetCoefficient("CleR", value * LambdaNPm2, 2, 2, 0, 0);
            
        Cle_3322r_LNP = value;
        SMEFTEvol.SetCoefficient("CleR", value * LambdaNPm2, 2, 2, 1, 1);
            
        Cle_3333r_LNP = value;
        SMEFTEvol.SetCoefficient("CleR", value * LambdaNPm2, 2, 2, 2, 2);
                       
    }  else if (name.compare("Clu_LNP") == 0) {
        
        Clu_LNP = value;
        
        Clu_1111r_LNP = value;
        SMEFTEvol.SetCoefficient("CluR", value * LambdaNPm2, 0, 0, 0, 0);
        
        Clu_1122r_LNP = value;
        SMEFTEvol.SetCoefficient("CluR", value * LambdaNPm2, 0, 0, 1, 1);
            
        Clu_1133r_LNP = value;
        SMEFTEvol.SetCoefficient("CluR", value * LambdaNPm2, 0, 0, 2, 2);
            
        Clu_2211r_LNP = value;
        SMEFTEvol.SetCoefficient("CluR", value * LambdaNPm2, 1, 1, 0, 0);
        
        Clu_2222r_LNP = value;
        SMEFTEvol.SetCoefficient("CluR", value * LambdaNPm2, 1, 1, 1, 1);
            
        Clu_2233r_LNP = value;
        SMEFTEvol.SetCoefficient("CluR", value * LambdaNPm2, 1, 1, 2, 2);
            
        Clu_3311r_LNP = value;
        SMEFTEvol.SetCoefficient("CluR", value * LambdaNPm2, 2, 2, 0, 0);
            
        Clu_3322r_LNP = value;
        SMEFTEvol.SetCoefficient("CluR", value * LambdaNPm2, 2, 2, 1, 1);
            
        Clu_3333r_LNP = value;
        SMEFTEvol.SetCoefficient("CluR", value * LambdaNPm2, 2, 2, 2, 2);
                       
    }  else if (name.compare("Cld_LNP") == 0) {
        
        Cld_LNP = value;
        
        Cld_1111r_LNP = value;
        SMEFTEvol.SetCoefficient("CldR", value * LambdaNPm2, 0, 0, 0, 0);
        
        Cld_1122r_LNP = value;
        SMEFTEvol.SetCoefficient("CldR", value * LambdaNPm2, 0, 0, 1, 1);
            
        Cld_1133r_LNP = value;
        SMEFTEvol.SetCoefficient("CldR", value * LambdaNPm2, 0, 0, 2, 2);
            
        Cld_2211r_LNP = value;
        SMEFTEvol.SetCoefficient("CldR", value * LambdaNPm2, 1, 1, 0, 0);
        
        Cld_2222r_LNP = value;
        SMEFTEvol.SetCoefficient("CldR", value * LambdaNPm2, 1, 1, 1, 1);
            
        Cld_2233r_LNP = value;
        SMEFTEvol.SetCoefficient("CldR", value * LambdaNPm2, 1, 1, 2, 2);
            
        Cld_3311r_LNP = value;
        SMEFTEvol.SetCoefficient("CldR", value * LambdaNPm2, 2, 2, 0, 0);
            
        Cld_3322r_LNP = value;
        SMEFTEvol.SetCoefficient("CldR", value * LambdaNPm2, 2, 2, 1, 1);
            
        Cld_3333r_LNP = value;
        SMEFTEvol.SetCoefficient("CldR", value * LambdaNPm2, 2, 2, 2, 2);
                       
    }  else if (name.compare("Cqe_LNP") == 0) {
        
        Cqe_LNP = value;
        
        Cqe_1111r_LNP = value;
        SMEFTEvol.SetCoefficient("CqeR", value * LambdaNPm2, 0, 0, 0, 0);
        
        Cqe_1122r_LNP = value;
        SMEFTEvol.SetCoefficient("CqeR", value * LambdaNPm2, 0, 0, 1, 1);
            
        Cqe_1133r_LNP = value;
        SMEFTEvol.SetCoefficient("CqeR", value * LambdaNPm2, 0, 0, 2, 2);
            
        Cqe_2211r_LNP = value;
        SMEFTEvol.SetCoefficient("CqeR", value * LambdaNPm2, 1, 1, 0, 0);
        
        Cqe_2222r_LNP = value;
        SMEFTEvol.SetCoefficient("CqeR", value * LambdaNPm2, 1, 1, 1, 1);
            
        Cqe_2233r_LNP = value;
        SMEFTEvol.SetCoefficient("CqeR", value * LambdaNPm2, 1, 1, 2, 2);
            
        Cqe_3311r_LNP = value;
        SMEFTEvol.SetCoefficient("CqeR", value * LambdaNPm2, 2, 2, 0, 0);
            
        Cqe_3322r_LNP = value;
        SMEFTEvol.SetCoefficient("CqeR", value * LambdaNPm2, 2, 2, 1, 1);
            
        Cqe_3333r_LNP = value;
        SMEFTEvol.SetCoefficient("CqeR", value * LambdaNPm2, 2, 2, 2, 2);
                       
    } else if (name.compare("Cledq_LNP") == 0) {
        
        Cledq_LNP = value;
        
        Cledq_1111r_LNP = value;
        SMEFTEvol.SetCoefficient("CledqR", value * LambdaNPm2, 0, 0, 0, 0);
        
        Cledq_1122r_LNP = value;
        SMEFTEvol.SetCoefficient("CledqR", value * LambdaNPm2, 0, 0, 1, 1);
            
        Cledq_1133r_LNP = value;
        SMEFTEvol.SetCoefficient("CledqR", value * LambdaNPm2, 0, 0, 2, 2);
            
        Cledq_2211r_LNP = value;
        SMEFTEvol.SetCoefficient("CledqR", value * LambdaNPm2, 1, 1, 0, 0);
        
        Cledq_2222r_LNP = value;
        SMEFTEvol.SetCoefficient("CledqR", value * LambdaNPm2, 1, 1, 1, 1);
            
        Cledq_2233r_LNP = value;
        SMEFTEvol.SetCoefficient("CledqR", value * LambdaNPm2, 1, 1, 2, 2);
            
        Cledq_3311r_LNP = value;
        SMEFTEvol.SetCoefficient("CledqR", value * LambdaNPm2, 2, 2, 0, 0);
            
        Cledq_3322r_LNP = value;
        SMEFTEvol.SetCoefficient("CledqR", value * LambdaNPm2, 2, 2, 1, 1);
            
        Cledq_3333r_LNP = value;
        SMEFTEvol.SetCoefficient("CledqR", value * LambdaNPm2, 2, 2, 2, 2);
                       
    } else if (name.compare("Cqq1_LNP") == 0) {
        
        Cqq1_LNP = value;
        
        Cqq1_1111r_LNP = value;
        SMEFTEvol.SetCoefficient("Cqq1R", value * LambdaNPm2, 0, 0, 0, 0);
        
        Cqq1_1221r_LNP = value;
        SMEFTEvol.SetCoefficient("Cqq1R", value * LambdaNPm2, 0, 1, 1, 0);
           
        Cqq1_1122r_LNP = value;
        SMEFTEvol.SetCoefficient("Cqq1R", value * LambdaNPm2, 0, 0, 1, 1);
                
        Cqq1_1133r_LNP = value;
        SMEFTEvol.SetCoefficient("Cqq1R", value * LambdaNPm2, 0, 0, 2, 2);
        
        Cqq1_1331r_LNP = value;
        SMEFTEvol.SetCoefficient("Cqq1R", value * LambdaNPm2, 0, 2, 2, 0);
        
        Cqq1_2222r_LNP = value;
        SMEFTEvol.SetCoefficient("Cqq1R", value * LambdaNPm2, 1, 1, 1, 1);
        
        Cqq1_2233r_LNP = value;
        SMEFTEvol.SetCoefficient("Cqq1R", value * LambdaNPm2, 1, 1, 2, 2);
            
        Cqq1_2332r_LNP = value;
        SMEFTEvol.SetCoefficient("Cqq1R", value * LambdaNPm2, 1, 2, 2, 1);
        
        Cqq1_3333r_LNP = value;
        SMEFTEvol.SetCoefficient("Cqq1R", value * LambdaNPm2, 2, 2, 2, 2);

    } else if (name.compare("Cqq3_LNP") == 0) {
        
        Cqq3_LNP = value;
        
        Cqq3_1111r_LNP = value;
        SMEFTEvol.SetCoefficient("Cqq3R", value * LambdaNPm2, 0, 0, 0, 0);
        
        Cqq3_1221r_LNP = value;
        SMEFTEvol.SetCoefficient("Cqq3R", value * LambdaNPm2, 0, 1, 1, 0);
           
        Cqq3_1122r_LNP = value;
        SMEFTEvol.SetCoefficient("Cqq3R", value * LambdaNPm2, 0, 0, 1, 1);
                
        Cqq3_1133r_LNP = value;
        SMEFTEvol.SetCoefficient("Cqq3R", value * LambdaNPm2, 0, 0, 2, 2);
        
        Cqq3_1331r_LNP = value;
        SMEFTEvol.SetCoefficient("Cqq3R", value * LambdaNPm2, 0, 2, 2, 0);
        
        Cqq3_2222r_LNP = value;
        SMEFTEvol.SetCoefficient("Cqq3R", value * LambdaNPm2, 1, 1, 1, 1);
        
        Cqq3_2233r_LNP = value;
        SMEFTEvol.SetCoefficient("Cqq3R", value * LambdaNPm2, 1, 1, 2, 2);
            
        Cqq3_2332r_LNP = value;
        SMEFTEvol.SetCoefficient("Cqq3R", value * LambdaNPm2, 1, 2, 2, 1);
        
        Cqq3_3333r_LNP = value;
        SMEFTEvol.SetCoefficient("Cqq3R", value * LambdaNPm2, 2, 2, 2, 2);

    } else if (name.compare("Cuu_LNP") == 0) {
        
        Cuu_LNP = value;
        
        Cuu_1111r_LNP = value;
        SMEFTEvol.SetCoefficient("CuuR", value * LambdaNPm2, 0, 0, 0, 0);
        
        Cuu_1221r_LNP = value;
        SMEFTEvol.SetCoefficient("CuuR", value * LambdaNPm2, 0, 1, 1, 0);
           
        Cuu_1122r_LNP = value;
        SMEFTEvol.SetCoefficient("CuuR", value * LambdaNPm2, 0, 0, 1, 1);
                
        Cuu_1133r_LNP = value;
        SMEFTEvol.SetCoefficient("CuuR", value * LambdaNPm2, 0, 0, 2, 2);
        
        Cuu_1331r_LNP = value;
        SMEFTEvol.SetCoefficient("CuuR", value * LambdaNPm2, 0, 2, 2, 0);
        
        Cuu_2222r_LNP = value;
        SMEFTEvol.SetCoefficient("CuuR", value * LambdaNPm2, 1, 1, 1, 1);
        
        Cuu_2233r_LNP = value;
        SMEFTEvol.SetCoefficient("CuuR", value * LambdaNPm2, 1, 1, 2, 2);
            
        Cuu_2332r_LNP = value;
        SMEFTEvol.SetCoefficient("CuuR", value * LambdaNPm2, 1, 2, 2, 1);
        
        Cuu_3333r_LNP = value;
        SMEFTEvol.SetCoefficient("CuuR", value * LambdaNPm2, 2, 2, 2, 2);

    } else if (name.compare("Cdd_LNP") == 0) {
        
        Cdd_LNP = value;
        
        Cdd_1111r_LNP = value;
        SMEFTEvol.SetCoefficient("CddR", value * LambdaNPm2, 0, 0, 0, 0);
        
        Cdd_1221r_LNP = value;
        SMEFTEvol.SetCoefficient("CddR", value * LambdaNPm2, 0, 1, 1, 0);
           
        Cdd_1122r_LNP = value;
        SMEFTEvol.SetCoefficient("CddR", value * LambdaNPm2, 0, 0, 1, 1);
                
        Cdd_1133r_LNP = value;
        SMEFTEvol.SetCoefficient("CddR", value * LambdaNPm2, 0, 0, 2, 2);
        
        Cdd_1331r_LNP = value;
        SMEFTEvol.SetCoefficient("CddR", value * LambdaNPm2, 0, 2, 2, 0);
        
        Cdd_2222r_LNP = value;
        SMEFTEvol.SetCoefficient("CddR", value * LambdaNPm2, 1, 1, 1, 1);
        
        Cdd_2233r_LNP = value;
        SMEFTEvol.SetCoefficient("CddR", value * LambdaNPm2, 1, 1, 2, 2);
            
        Cdd_2332r_LNP = value;
        SMEFTEvol.SetCoefficient("CddR", value * LambdaNPm2, 1, 2, 2, 1);
        
        Cdd_3333r_LNP = value;
        SMEFTEvol.SetCoefficient("CddR", value * LambdaNPm2, 2, 2, 2, 2);

    } else if (name.compare("Cud1_LNP") == 0) {
        
        Cud1_LNP = value;
        
        Cud1_1111r_LNP = value;
        SMEFTEvol.SetCoefficient("Cud1R", value * LambdaNPm2, 0, 0, 0, 0);
        
        Cud1_1122r_LNP = value;
        SMEFTEvol.SetCoefficient("Cud1R", value * LambdaNPm2, 0, 0, 1, 1);
            
        Cud1_1133r_LNP = value;
        SMEFTEvol.SetCoefficient("Cud1R", value * LambdaNPm2, 0, 0, 2, 2);
            
        Cud1_2211r_LNP = value;
        SMEFTEvol.SetCoefficient("Cud1R", value * LambdaNPm2, 1, 1, 0, 0);
        
        Cud1_2222r_LNP = value;
        SMEFTEvol.SetCoefficient("Cud1R", value * LambdaNPm2, 1, 1, 1, 1);
            
        Cud1_2233r_LNP = value;
        SMEFTEvol.SetCoefficient("Cud1R", value * LambdaNPm2, 1, 1, 2, 2);
            
        Cud1_3311r_LNP = value;
        SMEFTEvol.SetCoefficient("Cud1R", value * LambdaNPm2, 2, 2, 0, 0);
            
        Cud1_3322r_LNP = value;
        SMEFTEvol.SetCoefficient("Cud1R", value * LambdaNPm2, 2, 2, 1, 1);
            
        Cud1_3333r_LNP = value;
        SMEFTEvol.SetCoefficient("Cud1R", value * LambdaNPm2, 2, 2, 2, 2);
                       
    } else if (name.compare("Cud8_LNP") == 0) {
        
        Cud8_LNP = value;
        
        Cud8_1111r_LNP = value;
        SMEFTEvol.SetCoefficient("Cud8R", value * LambdaNPm2, 0, 0, 0, 0);
        
        Cud8_1122r_LNP = value;
        SMEFTEvol.SetCoefficient("Cud8R", value * LambdaNPm2, 0, 0, 1, 1);
            
        Cud8_1133r_LNP = value;
        SMEFTEvol.SetCoefficient("Cud8R", value * LambdaNPm2, 0, 0, 2, 2);
            
        Cud8_2211r_LNP = value;
        SMEFTEvol.SetCoefficient("Cud8R", value * LambdaNPm2, 1, 1, 0, 0);
        
        Cud8_2222r_LNP = value;
        SMEFTEvol.SetCoefficient("Cud8R", value * LambdaNPm2, 1, 1, 1, 1);
            
        Cud8_2233r_LNP = value;
        SMEFTEvol.SetCoefficient("Cud8R", value * LambdaNPm2, 1, 1, 2, 2);
            
        Cud8_3311r_LNP = value;
        SMEFTEvol.SetCoefficient("Cud8R", value * LambdaNPm2, 2, 2, 0, 0);
            
        Cud8_3322r_LNP = value;
        SMEFTEvol.SetCoefficient("Cud8R", value * LambdaNPm2, 2, 2, 1, 1);
            
        Cud8_3333r_LNP = value;
        SMEFTEvol.SetCoefficient("Cud8R", value * LambdaNPm2, 2, 2, 2, 2);
                       
    } else if (name.compare("Cqu1_LNP") == 0) {
        
        Cqu1_LNP = value;
        
        Cqu1_1111r_LNP = value;
        SMEFTEvol.SetCoefficient("Cqu1R", value * LambdaNPm2, 0, 0, 0, 0);
        
        Cqu1_1122r_LNP = value;
        SMEFTEvol.SetCoefficient("Cqu1R", value * LambdaNPm2, 0, 0, 1, 1);
            
        Cqu1_1133r_LNP = value;
        SMEFTEvol.SetCoefficient("Cqu1R", value * LambdaNPm2, 0, 0, 2, 2);
            
        Cqu1_2211r_LNP = value;
        SMEFTEvol.SetCoefficient("Cqu1R", value * LambdaNPm2, 1, 1, 0, 0);
        
        Cqu1_2222r_LNP = value;
        SMEFTEvol.SetCoefficient("Cqu1R", value * LambdaNPm2, 1, 1, 1, 1);
            
        Cqu1_2233r_LNP = value;
        SMEFTEvol.SetCoefficient("Cqu1R", value * LambdaNPm2, 1, 1, 2, 2);
            
        Cqu1_3311r_LNP = value;
        SMEFTEvol.SetCoefficient("Cqu1R", value * LambdaNPm2, 2, 2, 0, 0);
            
        Cqu1_3322r_LNP = value;
        SMEFTEvol.SetCoefficient("Cqu1R", value * LambdaNPm2, 2, 2, 1, 1);
            
        Cqu1_3333r_LNP = value;
        SMEFTEvol.SetCoefficient("Cqu1R", value * LambdaNPm2, 2, 2, 2, 2);
                       
    } else if (name.compare("Cqu8_LNP") == 0) {
        
        Cqu8_LNP = value;
        
        Cqu8_1111r_LNP = value;
        SMEFTEvol.SetCoefficient("Cqu8R", value * LambdaNPm2, 0, 0, 0, 0);
        
        Cqu8_1122r_LNP = value;
        SMEFTEvol.SetCoefficient("Cqu8R", value * LambdaNPm2, 0, 0, 1, 1);
            
        Cqu8_1133r_LNP = value;
        SMEFTEvol.SetCoefficient("Cqu8R", value * LambdaNPm2, 0, 0, 2, 2);
            
        Cqu8_2211r_LNP = value;
        SMEFTEvol.SetCoefficient("Cqu8R", value * LambdaNPm2, 1, 1, 0, 0);
        
        Cqu8_2222r_LNP = value;
        SMEFTEvol.SetCoefficient("Cqu8R", value * LambdaNPm2, 1, 1, 1, 1);
            
        Cqu8_2233r_LNP = value;
        SMEFTEvol.SetCoefficient("Cqu8R", value * LambdaNPm2, 1, 1, 2, 2);
            
        Cqu8_3311r_LNP = value;
        SMEFTEvol.SetCoefficient("Cqu8R", value * LambdaNPm2, 2, 2, 0, 0);
            
        Cqu8_3322r_LNP = value;
        SMEFTEvol.SetCoefficient("Cqu8R", value * LambdaNPm2, 2, 2, 1, 1);
            
        Cqu8_3333r_LNP = value;
        SMEFTEvol.SetCoefficient("Cqu8R", value * LambdaNPm2, 2, 2, 2, 2);
                       
    } else if (name.compare("Cqd1_LNP") == 0) {
        
        Cqd1_LNP = value;
        
        Cqd1_1111r_LNP = value;
        SMEFTEvol.SetCoefficient("Cqd1R", value * LambdaNPm2, 0, 0, 0, 0);
        
        Cqd1_1122r_LNP = value;
        SMEFTEvol.SetCoefficient("Cqd1R", value * LambdaNPm2, 0, 0, 1, 1);
            
        Cqd1_1133r_LNP = value;
        SMEFTEvol.SetCoefficient("Cqd1R", value * LambdaNPm2, 0, 0, 2, 2);
            
        Cqd1_2211r_LNP = value;
        SMEFTEvol.SetCoefficient("Cqd1R", value * LambdaNPm2, 1, 1, 0, 0);
        
        Cqd1_2222r_LNP = value;
        SMEFTEvol.SetCoefficient("Cqd1R", value * LambdaNPm2, 1, 1, 1, 1);
            
        Cqd1_2233r_LNP = value;
        SMEFTEvol.SetCoefficient("Cqd1R", value * LambdaNPm2, 1, 1, 2, 2);
            
        Cqd1_3311r_LNP = value;
        SMEFTEvol.SetCoefficient("Cqd1R", value * LambdaNPm2, 2, 2, 0, 0);
            
        Cqd1_3322r_LNP = value;
        SMEFTEvol.SetCoefficient("Cqd1R", value * LambdaNPm2, 2, 2, 1, 1);
            
        Cqd1_3333r_LNP = value;
        SMEFTEvol.SetCoefficient("Cqd1R", value * LambdaNPm2, 2, 2, 2, 2);
                       
    } else if (name.compare("Cqd8_LNP") == 0) {
        
        Cqd8_LNP = value;
        
        Cqd8_1111r_LNP = value;
        SMEFTEvol.SetCoefficient("Cqd8R", value * LambdaNPm2, 0, 0, 0, 0);
        
        Cqd8_1122r_LNP = value;
        SMEFTEvol.SetCoefficient("Cqd8R", value * LambdaNPm2, 0, 0, 1, 1);
            
        Cqd8_1133r_LNP = value;
        SMEFTEvol.SetCoefficient("Cqd8R", value * LambdaNPm2, 0, 0, 2, 2);
            
        Cqd8_2211r_LNP = value;
        SMEFTEvol.SetCoefficient("Cqd8R", value * LambdaNPm2, 1, 1, 0, 0);
        
        Cqd8_2222r_LNP = value;
        SMEFTEvol.SetCoefficient("Cqd8R", value * LambdaNPm2, 1, 1, 1, 1);
            
        Cqd8_2233r_LNP = value;
        SMEFTEvol.SetCoefficient("Cqd8R", value * LambdaNPm2, 1, 1, 2, 2);
            
        Cqd8_3311r_LNP = value;
        SMEFTEvol.SetCoefficient("Cqd8R", value * LambdaNPm2, 2, 2, 0, 0);
            
        Cqd8_3322r_LNP = value;
        SMEFTEvol.SetCoefficient("Cqd8R", value * LambdaNPm2, 2, 2, 1, 1);
            
        Cqd8_3333r_LNP = value;
        SMEFTEvol.SetCoefficient("Cqd8R", value * LambdaNPm2, 2, 2, 2, 2);
                       
    } else if (name.compare("Cquqd1_LNP") == 0) {
        
        Cquqd1_LNP = value;
        
        Cquqd1_1111r_LNP = value;
        SMEFTEvol.SetCoefficient("Cquqd1R", value * LambdaNPm2, 0, 0, 0, 0);
        
        Cquqd1_1122r_LNP = value;
        SMEFTEvol.SetCoefficient("Cquqd1R", value * LambdaNPm2, 0, 0, 1, 1);
            
        Cquqd1_1133r_LNP = value;
        SMEFTEvol.SetCoefficient("Cquqd1R", value * LambdaNPm2, 0, 0, 2, 2);
            
        Cquqd1_2211r_LNP = value;
        SMEFTEvol.SetCoefficient("Cquqd1R", value * LambdaNPm2, 1, 1, 0, 0);
        
        Cquqd1_2222r_LNP = value;
        SMEFTEvol.SetCoefficient("Cquqd1R", value * LambdaNPm2, 1, 1, 1, 1);
            
        Cquqd1_2233r_LNP = value;
        SMEFTEvol.SetCoefficient("Cquqd1R", value * LambdaNPm2, 1, 1, 2, 2);
            
        Cquqd1_3311r_LNP = value;
        SMEFTEvol.SetCoefficient("Cquqd1R", value * LambdaNPm2, 2, 2, 0, 0);
            
        Cquqd1_3322r_LNP = value;
        SMEFTEvol.SetCoefficient("Cquqd1R", value * LambdaNPm2, 2, 2, 1, 1);
            
        Cquqd1_3333r_LNP = value;
        SMEFTEvol.SetCoefficient("Cquqd1R", value * LambdaNPm2, 2, 2, 2, 2);
                       
    } else if (name.compare("Cquqd8_LNP") == 0) {
        
        Cquqd8_LNP = value;
        
        Cquqd8_1111r_LNP = value;
        SMEFTEvol.SetCoefficient("Cquqd8R", value * LambdaNPm2, 0, 0, 0, 0);
        
        Cquqd8_1122r_LNP = value;
        SMEFTEvol.SetCoefficient("Cquqd8R", value * LambdaNPm2, 0, 0, 1, 1);
            
        Cquqd8_1133r_LNP = value;
        SMEFTEvol.SetCoefficient("Cquqd8R", value * LambdaNPm2, 0, 0, 2, 2);
            
        Cquqd8_2211r_LNP = value;
        SMEFTEvol.SetCoefficient("Cquqd8R", value * LambdaNPm2, 1, 1, 0, 0);
        
        Cquqd8_2222r_LNP = value;
        SMEFTEvol.SetCoefficient("Cquqd8R", value * LambdaNPm2, 1, 1, 1, 1);
            
        Cquqd8_2233r_LNP = value;
        SMEFTEvol.SetCoefficient("Cquqd8R", value * LambdaNPm2, 1, 1, 2, 2);
            
        Cquqd8_3311r_LNP = value;
        SMEFTEvol.SetCoefficient("Cquqd8R", value * LambdaNPm2, 2, 2, 0, 0);
            
        Cquqd8_3322r_LNP = value;
        SMEFTEvol.SetCoefficient("Cquqd8R", value * LambdaNPm2, 2, 2, 1, 1);
            
        Cquqd8_3333r_LNP = value;
        SMEFTEvol.SetCoefficient("Cquqd8R", value * LambdaNPm2, 2, 2, 2, 2);
                       
    } else if (name.compare("Clequ1_LNP") == 0) {
        
        Clequ1_LNP = value;
        
        Clequ1_1111r_LNP = value;
        SMEFTEvol.SetCoefficient("Clequ1R", value * LambdaNPm2, 0, 0, 0, 0);
        
        Clequ1_1122r_LNP = value;
        SMEFTEvol.SetCoefficient("Clequ1R", value * LambdaNPm2, 0, 0, 1, 1);
            
        Clequ1_1133r_LNP = value;
        SMEFTEvol.SetCoefficient("Clequ1R", value * LambdaNPm2, 0, 0, 2, 2);
            
        Clequ1_2211r_LNP = value;
        SMEFTEvol.SetCoefficient("Clequ1R", value * LambdaNPm2, 1, 1, 0, 0);
        
        Clequ1_2222r_LNP = value;
        SMEFTEvol.SetCoefficient("Clequ1R", value * LambdaNPm2, 1, 1, 1, 1);
            
        Clequ1_2233r_LNP = value;
        SMEFTEvol.SetCoefficient("Clequ1R", value * LambdaNPm2, 1, 1, 2, 2);
            
        Clequ1_3311r_LNP = value;
        SMEFTEvol.SetCoefficient("Clequ1R", value * LambdaNPm2, 2, 2, 0, 0);
            
        Clequ1_3322r_LNP = value;
        SMEFTEvol.SetCoefficient("Clequ1R", value * LambdaNPm2, 2, 2, 1, 1);
            
        Clequ1_3333r_LNP = value;
        SMEFTEvol.SetCoefficient("Clequ1R", value * LambdaNPm2, 2, 2, 2, 2);
                       
    }  else if (name.compare("Clequ3_LNP") == 0) {
        
        Clequ3_LNP = value;
        
        Clequ3_1111r_LNP = value;
        SMEFTEvol.SetCoefficient("Clequ3R", value * LambdaNPm2, 0, 0, 0, 0);
        
        Clequ3_1122r_LNP = value;
        SMEFTEvol.SetCoefficient("Clequ3R", value * LambdaNPm2, 0, 0, 1, 1);
            
        Clequ3_1133r_LNP = value;
        SMEFTEvol.SetCoefficient("Clequ3R", value * LambdaNPm2, 0, 0, 2, 2);
            
        Clequ3_2211r_LNP = value;
        SMEFTEvol.SetCoefficient("Clequ3R", value * LambdaNPm2, 1, 1, 0, 0);
        
        Clequ3_2222r_LNP = value;
        SMEFTEvol.SetCoefficient("Clequ3R", value * LambdaNPm2, 1, 1, 1, 1);
            
        Clequ3_2233r_LNP = value;
        SMEFTEvol.SetCoefficient("Clequ3R", value * LambdaNPm2, 1, 1, 2, 2);
            
        Clequ3_3311r_LNP = value;
        SMEFTEvol.SetCoefficient("Clequ3R", value * LambdaNPm2, 2, 2, 0, 0);
            
        Clequ3_3322r_LNP = value;
        SMEFTEvol.SetCoefficient("Clequ3R", value * LambdaNPm2, 2, 2, 1, 1);
            
        Clequ3_3333r_LNP = value;
        SMEFTEvol.SetCoefficient("Clequ3R", value * LambdaNPm2, 2, 2, 2, 2);
                       
     
    } else if (name.compare("Lambda_NP") == 0)
        Lambda_NP = value;
    //Let's go directly for the NPbase setParameter, we should change it to the NPSMEFTd6General removing reading the unnecessary parameters
    else
        NPbase::setParameter(name, value);


}