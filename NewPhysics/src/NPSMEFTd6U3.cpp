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
    "CHl1_LNP","CHl3_LNP","CHe_LNP","CHq1_LNP","CHq3_LNP","CHu_LNP","CHd_LNP",
    "Cll_aabb_LNP","Cll_abba_LNP","Clq1_LNP","Clq3_LNP","Cee_LNP","Ceu_LNP","Ced_LNP","Cle_LNP","Clu_LNP","Cld_LNP","Cqe_LNP",
    "Cqq1_aabb_LNP","Cqq1_abba_LNP", "Cqq3_aabb_LNP", "Cqq3_abba_LNP", "Cuu_aabb_LNP", "Cuu_abba_LNP", "Cdd_aabb_LNP", "Cdd_abba_LNP",
    "Cud1_LNP", "Cud8_LNP", "Cqu1_LNP", "Cqu8_LNP", "Cqd1_LNP", "Cqd8_LNP", "Lambda_NP"
};

NPSMEFTd6U3::NPSMEFTd6U3()
: NPSMEFTd6General()
{
    setModelName("NPSMEFTd6U3");
   
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

    
    
    ModelParamMap.insert(std::make_pair("Cll_aabb_LNP", std::cref(Cll_aabb_LNP)));
    ModelParamMap.insert(std::make_pair("Cll_abba_LNP", std::cref(Cll_abba_LNP)));
    ModelParamMap.insert(std::make_pair("Clq1_LNP", std::cref(Clq1_LNP)));
    ModelParamMap.insert(std::make_pair("Clq3_LNP", std::cref(Clq3_LNP)));
    ModelParamMap.insert(std::make_pair("Cee_LNP", std::cref(Cee_LNP)));
    ModelParamMap.insert(std::make_pair("Ceu_LNP", std::cref(Ceu_LNP)));
    ModelParamMap.insert(std::make_pair("Ced_LNP", std::cref(Ced_LNP)));
    ModelParamMap.insert(std::make_pair("Cle_LNP", std::cref(Cle_LNP)));
    ModelParamMap.insert(std::make_pair("Clu_LNP", std::cref(Clu_LNP)));
    ModelParamMap.insert(std::make_pair("Cld_LNP", std::cref(Cld_LNP)));
    ModelParamMap.insert(std::make_pair("Cqe_LNP", std::cref(Cqe_LNP)));
    
    
    ModelParamMap.insert(std::make_pair("Cqq1_aabb_LNP", std::cref(Cqq1_aabb_LNP)));
    ModelParamMap.insert(std::make_pair("Cqq1_abba_LNP", std::cref(Cqq1_abba_LNP)));
    ModelParamMap.insert(std::make_pair("Cqq3_aabb_LNP", std::cref(Cqq3_aabb_LNP)));
    ModelParamMap.insert(std::make_pair("Cqq3_abba_LNP", std::cref(Cqq3_abba_LNP)));
    ModelParamMap.insert(std::make_pair("Cuu_aabb_LNP", std::cref(Cuu_aabb_LNP)));
    ModelParamMap.insert(std::make_pair("Cuu_abba_LNP", std::cref(Cuu_abba_LNP)));
    ModelParamMap.insert(std::make_pair("Cdd_aabb_LNP", std::cref(Cdd_aabb_LNP)));
    ModelParamMap.insert(std::make_pair("Cdd_abba_LNP", std::cref(Cdd_abba_LNP)));

    ModelParamMap.insert(std::make_pair("Cud1_LNP", std::cref(Cud1_LNP)));
    ModelParamMap.insert(std::make_pair("Cud8_LNP", std::cref(Cud8_LNP)));
    ModelParamMap.insert(std::make_pair("Cqu1_LNP", std::cref(Cqu1_LNP)));
    ModelParamMap.insert(std::make_pair("Cqu8_LNP", std::cref(Cqu8_LNP)));
    ModelParamMap.insert(std::make_pair("Cqd1_LNP", std::cref(Cqd1_LNP)));
    ModelParamMap.insert(std::make_pair("Cqd8_LNP", std::cref(Cqd8_LNP)));
    
    
}




void NPSMEFTd6U3::setParameter(const std::string name, const double& value)
{



    if (name.compare("CG_LNP") == 0) {
        CG_LNP = value;

    } else if (name.compare("CW_LNP") == 0) {
        CW_LNP = value;

    } else if (name.compare("CHG_LNP") == 0) {
        CHG_LNP = value;

    } else if (name.compare("CHW_LNP") == 0) {
        CHW_LNP = value;

    } else if (name.compare("CHB_LNP") == 0) {
        CHB_LNP = value;

    } else if (name.compare("CHWB_LNP") == 0) {
        CHWB_LNP = value;

    } else if (name.compare("CHD_LNP") == 0) {
        CHD_LNP = value;

    } else if (name.compare("CHbox_LNP") == 0) {
        CHbox_LNP = value;

    } else if (name.compare("CH_LNP") == 0) {
        CH_LNP = value;

    } else if (name.compare("CHl1_LNP") == 0) {
        
        //std::cout<<"\033[1;33m   LambdaNPm2 =  \033[0m "<< LambdaNPm2 <<std::endl;
                
        CHl1_LNP = value;

    } else if (name.compare("CHl3_LNP") == 0) {
        
        CHl3_LNP = value;

    } else if (name.compare("CHe_LNP") == 0) {
        
        CHe_LNP = value;

    } else if (name.compare("CHq1_LNP") == 0) {
        
        CHq1_LNP = value;

    } else if (name.compare("CHq3_LNP") == 0) {
        
        CHq3_LNP = value;

    } else if (name.compare("CHu_LNP") == 0) {
        
        CHu_LNP = value;

    } else if (name.compare("CHd_LNP") == 0) {
        
        CHd_LNP = value;
  
    } else if (name.compare("Cll_aabb_LNP") == 0) {
        
        Cll_aabb_LNP = value; 
        
    } else if (name.compare("Cll_abba_LNP") == 0) {
        
        Cll_abba_LNP = value;   
        
    } else if (name.compare("Clq1_LNP") == 0) {
        
        Clq1_LNP = value;
 
    } else if (name.compare("Clq3_LNP") == 0) {
        
        Clq3_LNP = value;
     
    } else if (name.compare("Cee_LNP") == 0) {
        
        Cee_LNP = value;
    
    } else if (name.compare("Ceu_LNP") == 0) {
        
        Ceu_LNP = value;
                       
    } else if (name.compare("Ced_LNP") == 0) {
        
        Ced_LNP = value;
                       
    }  else if (name.compare("Cle_LNP") == 0) {
        
        Cle_LNP = value;
                       
    }  else if (name.compare("Clu_LNP") == 0) {
        
        Clu_LNP = value;
                       
    }  else if (name.compare("Cld_LNP") == 0) {
        
        Cld_LNP = value;
                       
    }  else if (name.compare("Cqe_LNP") == 0) {
        
        Cqe_LNP = value;
                       
    } else if (name.compare("Cqq1_aabb_LNP") == 0) {
        
        Cqq1_aabb_LNP = value;
        
    } else if (name.compare("Cqq1_abba_LNP") == 0) {
        
        Cqq1_abba_LNP = value;
        

    } else if (name.compare("Cqq3_aabb_LNP") == 0) {
        
        Cqq3_aabb_LNP = value;
        
    } else if (name.compare("Cqq3_abba_LNP") == 0) {
        
        Cqq3_abba_LNP = value;
        

    } else if (name.compare("Cuu_aabb_LNP") == 0) {
        
        Cuu_aabb_LNP = value;
        
    } else if (name.compare("Cuu_abba_LNP") == 0) {
        
        Cuu_abba_LNP = value;
        

    } else if (name.compare("Cdd_aabb_LNP") == 0) {
        
        Cdd_aabb_LNP = value;
        
    } else if (name.compare("Cdd_abba_LNP") == 0) {
        
        Cdd_abba_LNP = value;
        

    } else if (name.compare("Cud1_LNP") == 0) {
        
        Cud1_LNP = value;
                       
    } else if (name.compare("Cud8_LNP") == 0) {
        
        Cud8_LNP = value;
                       
    } else if (name.compare("Cqu1_LNP") == 0) {
        
        Cqu1_LNP = value;
        
                       
    } else if (name.compare("Cqu8_LNP") == 0) {
        
        Cqu8_LNP = value;
        
                       
    } else if (name.compare("Cqd1_LNP") == 0) {
        
        Cqd1_LNP = value;
                       
    } else if (name.compare("Cqd8_LNP") == 0) {
        
        Cqd8_LNP = value;
                              
     
    } else if (name.compare("Lambda_NP") == 0)
        Lambda_NP = value;
    else
        NPSMEFTd6General::setParameter(name, value);


}



void NPSMEFTd6U3::setNPSMEFTd6GeneralParameters()
{

    
    CHl1_11r_LNP = CHl1_LNP;
    CHl1_22r_LNP = CHl1_LNP;
    CHl1_33r_LNP = CHl1_LNP;
    
    CHl3_11r_LNP = CHl3_LNP;
    CHl3_22r_LNP = CHl3_LNP;
    CHl3_33r_LNP = CHl3_LNP;
    
    CHe_11r_LNP = CHe_LNP;
    CHe_22r_LNP = CHe_LNP;
    CHe_33r_LNP = CHe_LNP;
    
    CHq1_11r_LNP = CHq1_LNP;
    CHq1_22r_LNP = CHq1_LNP;
    CHq1_33r_LNP = CHq1_LNP;
    
    CHq3_11r_LNP = CHq3_LNP;
    CHq3_22r_LNP = CHq3_LNP;
    CHq3_33r_LNP = CHq3_LNP;
    
    CHu_11r_LNP = CHu_LNP;
    CHu_22r_LNP = CHu_LNP;
    CHu_33r_LNP = CHu_LNP;
    
    CHd_11r_LNP = CHd_LNP;
    CHd_22r_LNP = CHd_LNP;
    CHd_33r_LNP = CHd_LNP;
    
    
    
    
    
    
    Cll_1111r_LNP = Cll_aabb_LNP+Cll_abba_LNP;
    Cll_2222r_LNP = Cll_aabb_LNP+Cll_abba_LNP;
    Cll_3333r_LNP = Cll_aabb_LNP+Cll_abba_LNP;

    Cll_1122r_LNP = Cll_aabb_LNP;
    Cll_1133r_LNP = Cll_aabb_LNP;
    Cll_2233r_LNP = Cll_aabb_LNP;
    
    Cll_1221r_LNP = Cll_abba_LNP;
    Cll_1331r_LNP = Cll_abba_LNP;
    Cll_2332r_LNP = Cll_abba_LNP;
    //////////////////////////////
    
    
    
    Clq1_1111r_LNP = Clq1_LNP;
    Clq1_1122r_LNP = Clq1_LNP;
    Clq1_1133r_LNP = Clq1_LNP;
    Clq1_2211r_LNP = Clq1_LNP;
    Clq1_2222r_LNP = Clq1_LNP;
    Clq1_2233r_LNP = Clq1_LNP;
    Clq1_3311r_LNP = Clq1_LNP;
    Clq1_3322r_LNP = Clq1_LNP;
    Clq1_3333r_LNP = Clq1_LNP;
    
    
    
    Clq3_1111r_LNP = Clq3_LNP;
    Clq3_1122r_LNP = Clq3_LNP;
    Clq3_1133r_LNP = Clq3_LNP;
    Clq3_2211r_LNP = Clq3_LNP;
    Clq3_2222r_LNP = Clq3_LNP;
    Clq3_2233r_LNP = Clq3_LNP;
    Clq3_3311r_LNP = Clq3_LNP;
    Clq3_3322r_LNP = Clq3_LNP;
    Clq3_3333r_LNP = Clq3_LNP;
    
    
    
    Cee_1111r_LNP = Cee_LNP;
    Cee_1122r_LNP = Cee_LNP;
    Cee_1133r_LNP = Cee_LNP;
    Cee_2222r_LNP = Cee_LNP;
    Cee_2233r_LNP = Cee_LNP;
    Cee_3333r_LNP = Cee_LNP;
    
    
    
    Ceu_1111r_LNP = Ceu_LNP;
    Ceu_1122r_LNP = Ceu_LNP;
    Ceu_1133r_LNP = Ceu_LNP;
    Ceu_2211r_LNP = Ceu_LNP;
    Ceu_2222r_LNP = Ceu_LNP;
    Ceu_2233r_LNP = Ceu_LNP;
    Ceu_3311r_LNP = Ceu_LNP;
    Ceu_3322r_LNP = Ceu_LNP;
    Ceu_3333r_LNP = Ceu_LNP;

    
    
    Ced_1111r_LNP = Ced_LNP;
    Ced_1122r_LNP = Ced_LNP;
    Ced_1133r_LNP = Ced_LNP;
    Ced_2211r_LNP = Ced_LNP;
    Ced_2222r_LNP = Ced_LNP;
    Ced_2233r_LNP = Ced_LNP;
    Ced_3311r_LNP = Ced_LNP;
    Ced_3322r_LNP = Ced_LNP;
    Ced_3333r_LNP = Ced_LNP;

    
    
    Cle_1111r_LNP = Cle_LNP;
    Cle_1122r_LNP = Cle_LNP;
    Cle_1133r_LNP = Cle_LNP;
    Cle_2211r_LNP = Cle_LNP;
    Cle_2222r_LNP = Cle_LNP;
    Cle_2233r_LNP = Cle_LNP;
    Cle_3311r_LNP = Cle_LNP;
    Cle_3322r_LNP = Cle_LNP;
    Cle_3333r_LNP = Cle_LNP; 

    
    
    Clu_1111r_LNP = Clu_LNP;
    Clu_1122r_LNP = Clu_LNP;
    Clu_1133r_LNP = Clu_LNP;
    Clu_2211r_LNP = Clu_LNP;
    Clu_2222r_LNP = Clu_LNP;
    Clu_2233r_LNP = Clu_LNP;
    Clu_3311r_LNP = Clu_LNP;
    Clu_3322r_LNP = Clu_LNP;
    Clu_3333r_LNP = Clu_LNP; 
    

    
    
    Cld_1111r_LNP = Cld_LNP;
    Cld_1122r_LNP = Cld_LNP;
    Cld_1133r_LNP = Cld_LNP;
    Cld_2211r_LNP = Cld_LNP;
    Cld_2222r_LNP = Cld_LNP;
    Cld_2233r_LNP = Cld_LNP;
    Cld_3311r_LNP = Cld_LNP;
    Cld_3322r_LNP = Cld_LNP;
    Cld_3333r_LNP = Cld_LNP; 
    

    
    
    Cqe_1111r_LNP = Cqe_LNP;
    Cqe_1122r_LNP = Cqe_LNP;
    Cqe_1133r_LNP = Cqe_LNP;
    Cqe_2211r_LNP = Cqe_LNP;
    Cqe_2222r_LNP = Cqe_LNP;
    Cqe_2233r_LNP = Cqe_LNP;
    Cqe_3311r_LNP = Cqe_LNP;
    Cqe_3322r_LNP = Cqe_LNP;
    Cqe_3333r_LNP = Cqe_LNP; 
    
    
    
    
    
    
    
    
    Cqq1_1111r_LNP = Cqq1_aabb_LNP+Cqq1_abba_LNP;
    Cqq1_2222r_LNP = Cqq1_aabb_LNP+Cqq1_abba_LNP;
    Cqq1_3333r_LNP = Cqq1_aabb_LNP+Cqq1_abba_LNP;

    Cqq1_1122r_LNP = Cqq1_aabb_LNP;
    Cqq1_1133r_LNP = Cqq1_aabb_LNP;
    Cqq1_2233r_LNP = Cqq1_aabb_LNP;
    
    Cqq1_1221r_LNP = Cqq1_abba_LNP;
    Cqq1_1331r_LNP = Cqq1_abba_LNP;
    Cqq1_2332r_LNP = Cqq1_abba_LNP;
    
    
    
    
    
    Cqq3_1111r_LNP = Cqq3_aabb_LNP+Cqq3_abba_LNP;
    Cqq3_2222r_LNP = Cqq3_aabb_LNP+Cqq3_abba_LNP;
    Cqq3_3333r_LNP = Cqq3_aabb_LNP+Cqq3_abba_LNP;

    Cqq3_1122r_LNP = Cqq3_aabb_LNP;
    Cqq3_1133r_LNP = Cqq3_aabb_LNP;
    Cqq3_2233r_LNP = Cqq3_aabb_LNP;
    
    Cqq3_1221r_LNP = Cqq3_abba_LNP;
    Cqq3_1331r_LNP = Cqq3_abba_LNP;
    Cqq3_2332r_LNP = Cqq3_abba_LNP;
   
    
    
    
    
    Cuu_1111r_LNP = Cuu_aabb_LNP+Cuu_abba_LNP;
    Cuu_2222r_LNP = Cuu_aabb_LNP+Cuu_abba_LNP;
    Cuu_3333r_LNP = Cuu_aabb_LNP+Cuu_abba_LNP;

    Cuu_1122r_LNP = Cuu_aabb_LNP;
    Cuu_1133r_LNP = Cuu_aabb_LNP;
    Cuu_2233r_LNP = Cuu_aabb_LNP;
    
    Cuu_1221r_LNP = Cuu_abba_LNP;
    Cuu_1331r_LNP = Cuu_abba_LNP;
    Cuu_2332r_LNP = Cuu_abba_LNP;
    
    
    
    
    
    
    Cdd_1111r_LNP = Cdd_aabb_LNP+Cdd_abba_LNP;
    Cdd_2222r_LNP = Cdd_aabb_LNP+Cdd_abba_LNP;
    Cdd_3333r_LNP = Cdd_aabb_LNP+Cdd_abba_LNP;

    Cdd_1122r_LNP = Cdd_aabb_LNP;
    Cdd_1133r_LNP = Cdd_aabb_LNP;
    Cdd_2233r_LNP = Cdd_aabb_LNP;
    
    Cdd_1221r_LNP = Cdd_abba_LNP;
    Cdd_1331r_LNP = Cdd_abba_LNP;
    Cdd_2332r_LNP = Cdd_abba_LNP;
    //////////////////////////////
    
    
    
    
    Cud1_1111r_LNP = Cud1_LNP;
    Cud1_1122r_LNP = Cud1_LNP;
    Cud1_1133r_LNP = Cud1_LNP;
    Cud1_2211r_LNP = Cud1_LNP;
    Cud1_2222r_LNP = Cud1_LNP;
    Cud1_2233r_LNP = Cud1_LNP;
    Cud1_3311r_LNP = Cud1_LNP;
    Cud1_3322r_LNP = Cud1_LNP;
    Cud1_3333r_LNP = Cud1_LNP; 
    
    
    
    
    Cud8_1111r_LNP = Cud8_LNP;
    Cud8_1122r_LNP = Cud8_LNP;
    Cud8_1133r_LNP = Cud8_LNP;
    Cud8_2211r_LNP = Cud8_LNP;
    Cud8_2222r_LNP = Cud8_LNP;
    Cud8_2233r_LNP = Cud8_LNP;
    Cud8_3311r_LNP = Cud8_LNP;
    Cud8_3322r_LNP = Cud8_LNP;
    Cud8_3333r_LNP = Cud8_LNP; 
    
    
    
    
    Cqu1_1111r_LNP = Cqu1_LNP;
    Cqu1_1122r_LNP = Cqu1_LNP;
    Cqu1_1133r_LNP = Cqu1_LNP;
    Cqu1_2211r_LNP = Cqu1_LNP;
    Cqu1_2222r_LNP = Cqu1_LNP;
    Cqu1_2233r_LNP = Cqu1_LNP;
    Cqu1_3311r_LNP = Cqu1_LNP;
    Cqu1_3322r_LNP = Cqu1_LNP;
    Cqu1_3333r_LNP = Cqu1_LNP; 
    
    
    
    
    Cqu8_1111r_LNP = Cqu8_LNP;
    Cqu8_1122r_LNP = Cqu8_LNP;
    Cqu8_1133r_LNP = Cqu8_LNP;
    Cqu8_2211r_LNP = Cqu8_LNP;
    Cqu8_2222r_LNP = Cqu8_LNP;
    Cqu8_2233r_LNP = Cqu8_LNP;
    Cqu8_3311r_LNP = Cqu8_LNP;
    Cqu8_3322r_LNP = Cqu8_LNP;
    Cqu8_3333r_LNP = Cqu8_LNP; 
    
    
    
    
    
    
    
    Cqd1_1111r_LNP = Cqd1_LNP;
    Cqd1_1122r_LNP = Cqd1_LNP;
    Cqd1_1133r_LNP = Cqd1_LNP;
    Cqd1_2211r_LNP = Cqd1_LNP;
    Cqd1_2222r_LNP = Cqd1_LNP;
    Cqd1_2233r_LNP = Cqd1_LNP;
    Cqd1_3311r_LNP = Cqd1_LNP;
    Cqd1_3322r_LNP = Cqd1_LNP;
    Cqd1_3333r_LNP = Cqd1_LNP; 
    
    
    
    
    Cqd8_1111r_LNP = Cqd8_LNP;
    Cqd8_1122r_LNP = Cqd8_LNP;
    Cqd8_1133r_LNP = Cqd8_LNP;
    Cqd8_2211r_LNP = Cqd8_LNP;
    Cqd8_2222r_LNP = Cqd8_LNP;
    Cqd8_2233r_LNP = Cqd8_LNP;
    Cqd8_3311r_LNP = Cqd8_LNP;
    Cqd8_3322r_LNP = Cqd8_LNP;
    Cqd8_3333r_LNP = Cqd8_LNP; 
    
}


bool NPSMEFTd6U3::PostUpdate()
{
    NPSMEFTd6General::GenerateSMInitialConditions();
    
    setNPSMEFTd6GeneralParameters();
    
    if (!NPSMEFTd6General::PostUpdate()) return (false);
    
    return (true);
        
}
