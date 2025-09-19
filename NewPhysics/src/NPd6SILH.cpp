/* Copyright (C) 2025 HEPfit Collaboration
* 
* For the licensing terms see doc/COPYING.
*
* Created on Tue 29 Jul 2025 16:05:01
*
*/

#include "NPd6SILH.h"


std::string NPd6SILH::NPd6SILHVars[NNPd6SILHVars] = {
    "cH_LNP", "cT_LNP", "c6_LNP", "cB_LNP", "cW_LNP", "c2B_LNP", "c2W_LNP", "c2G_LNP", "c3W_LNP", "c3G_LNP", "cHW_LNP", "cHB_LNP", "cgam_LNP", "cg_LNP", "cHq1_LNP", "cHq3_LNP", "cHt_LNP", "ctD_LNP", "cqD1_LNP", "cqD3_LNP", "cqq1_LNP", "cqq3_LNP", "cqt1_LNP", "cqt8_LNP", "ctt_LNP", "ctG_LNP", "ctB_LNP", "ctW_LNP", "cu_LNP", "cd_LNP", "ce_LNP", "Lambda_NP" 
};


NPd6SILH::NPd6SILH()
: NPSMEFTd6General(), YuUV(3,3,0.), YuUVhc(3,3,0.), YdUV(3,3,0.), YdUVhc(3,3,0.), YeUV(3,3,0.), YeUVhc(3,3,0.) {
    setModelName("NPd6SILH");

    FlagRGEci = true;

    ModelParamMap.insert(std::make_pair("cH_LNP", std::cref(cH_LNP)));
    ModelParamMap.insert(std::make_pair("cT_LNP", std::cref(cT_LNP)));
    ModelParamMap.insert(std::make_pair("c6_LNP", std::cref(c6_LNP)));
    ModelParamMap.insert(std::make_pair("cB_LNP", std::cref(cB_LNP)));
    ModelParamMap.insert(std::make_pair("cW_LNP", std::cref(cW_LNP)));
    ModelParamMap.insert(std::make_pair("c2B_LNP", std::cref(c2B_LNP)));
    ModelParamMap.insert(std::make_pair("c2W_LNP", std::cref(c2W_LNP)));
    ModelParamMap.insert(std::make_pair("c2G_LNP", std::cref(c2G_LNP)));
    ModelParamMap.insert(std::make_pair("c3W_LNP", std::cref(c3W_LNP)));
    ModelParamMap.insert(std::make_pair("c3G_LNP", std::cref(c3G_LNP)));
    ModelParamMap.insert(std::make_pair("cHW_LNP", std::cref(cHW_LNP)));
    ModelParamMap.insert(std::make_pair("cHB_LNP", std::cref(cHB_LNP)));
    ModelParamMap.insert(std::make_pair("cgam_LNP", std::cref(cgam_LNP)));
    ModelParamMap.insert(std::make_pair("cg_LNP", std::cref(cg_LNP)));
    ModelParamMap.insert(std::make_pair("cHq1_LNP", std::cref(cHq1_LNP)));
    ModelParamMap.insert(std::make_pair("cHq3_LNP", std::cref(cHq3_LNP)));
    ModelParamMap.insert(std::make_pair("cHt_LNP", std::cref(cHt_LNP)));
    ModelParamMap.insert(std::make_pair("ctD_LNP", std::cref(ctD_LNP)));
    ModelParamMap.insert(std::make_pair("cqD1_LNP", std::cref(cqD1_LNP)));
    ModelParamMap.insert(std::make_pair("cqD3_LNP", std::cref(cqD3_LNP)));
    ModelParamMap.insert(std::make_pair("cqq1_LNP", std::cref(cqq1_LNP)));
    ModelParamMap.insert(std::make_pair("cqq3_LNP", std::cref(cqq3_LNP)));
    ModelParamMap.insert(std::make_pair("cqt1_LNP", std::cref(cqt1_LNP)));
    ModelParamMap.insert(std::make_pair("cqt8_LNP", std::cref(cqt8_LNP)));
    ModelParamMap.insert(std::make_pair("ctt_LNP", std::cref(ctt_LNP)));
    ModelParamMap.insert(std::make_pair("ctG_LNP", std::cref(ctG_LNP)));
    ModelParamMap.insert(std::make_pair("ctB_LNP", std::cref(ctB_LNP)));
    ModelParamMap.insert(std::make_pair("ctW_LNP", std::cref(ctW_LNP)));
    ModelParamMap.insert(std::make_pair("cu_LNP", std::cref(cu_LNP)));
    ModelParamMap.insert(std::make_pair("cd_LNP", std::cref(cd_LNP)));
    ModelParamMap.insert(std::make_pair("ce_LNP", std::cref(ce_LNP)));

}


void NPd6SILH::setParameter(const std::string name, const double& value)
{

    if(name.compare("cH_LNP") == 0 ) { 

        cH_LNP = value;

    } else if(name.compare("cT_LNP")==0) {

        cT_LNP = value;

    } else if(name.compare("c6_LNP")==0) {

        c6_LNP = value;

    } else if(name.compare("cB_LNP")==0) {

        cB_LNP = value;

    } else if(name.compare("cW_LNP")==0) {

        cW_LNP = value;

    } else if(name.compare("c2B_LNP")==0) {

        c2B_LNP = value;

    } else if(name.compare("c2W_LNP")==0) {

        c2W_LNP = value;

    } else if(name.compare("c2G_LNP")==0) {

        c2G_LNP = value;

    } else if(name.compare("c3W_LNP")==0) {

        c3W_LNP = value;

    } else if(name.compare("c3G_LNP")==0) {

        c3G_LNP = value;

    } else if(name.compare("cHW_LNP")==0) {

        cHW_LNP = value;

    } else if(name.compare("cHB_LNP")==0) {

        cHB_LNP = value;

    } else if(name.compare("cgam_LNP")==0) {

        cgam_LNP = value;

    } else if(name.compare("cg_LNP")==0) {

        cg_LNP = value;

    } else if(name.compare("cHq1_LNP")==0) {

        cHq1_LNP = value;

    } else if(name.compare("cHq3_LNP")==0) {

        cHq3_LNP = value;

    } else if(name.compare("cHt_LNP")==0) {

        cHt_LNP = value;

    } else if(name.compare("ctD_LNP")==0) {

        ctD_LNP = value;

    } else if(name.compare("cqD1_LNP")==0) {

        cqD1_LNP = value;

    } else if(name.compare("cqD3_LNP")==0) {

        cqD3_LNP = value;

    } else if(name.compare("cqq1_LNP")==0) {

        cqq1_LNP = value;

    } else if(name.compare("cqq3_LNP")==0) {

        cqq3_LNP = value;

    } else if(name.compare("cqt1_LNP")==0) {

        cqt1_LNP = value;

    } else if(name.compare("cqt8_LNP")==0) {

        cqt8_LNP = value;

    } else if(name.compare("ctt_LNP")==0) {

        ctt_LNP = value;

    } else if(name.compare("ctG_LNP")==0) {

        ctG_LNP = value;

    } else if(name.compare("ctB_LNP")==0) {

        ctB_LNP = value;

    } else if(name.compare("ctW_LNP")==0) {

        ctW_LNP = value;

    } else if(name.compare("cu_LNP")==0) {

        cu_LNP = value;

    } else if(name.compare("cd_LNP")==0) {

        cd_LNP = value;

    } else if(name.compare("ce_LNP")==0) {

        ce_LNP = value;

    } else if(name.compare("Lambda_NP")==0) {

        Lambda_NP = value;

    } else {
        NPSMEFTd6General::setParameter(name, value);
    }
}


bool NPd6SILH::setFlag(const std::string name, const bool value) {
    bool res = false;
    if (name.compare("RGEci") == 0) {
        FlagRGEci = value;
        //We need to fix FlagRGEci also in the NPSMEFTd6General
        res = NPSMEFTd6General::setFlag(name, value);

    } else
        res = NPSMEFTd6General::setFlag(name, value);

    return (res);
}


void NPd6SILH::setNPSMEFTd6GeneralParameters()
{    
    // SM parameters at the UV scale
    g1UV = getSMEFTCoeffEW("g1");
    g2UV = getSMEFTCoeffEW("g2");
    g3UV = getSMEFTCoeffEW("g3");
    lambdaHUV = getSMEFTCoeffEW("lambda");

    g1UV2 = g1UV * g1UV;
    g2UV2 = g2UV * g2UV;
    g3UV2 = g3UV * g3UV;

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            YuUV.assignre(i, j, getSMEFTCoeffEW("YuR", i, j) );
            YuUV.assignim(i, j, getSMEFTCoeffEW("YuI", i, j) );
            YdUV.assignre(i, j, getSMEFTCoeffEW("YdR", i, j) );
            YdUV.assignim(i, j, getSMEFTCoeffEW("YdI", i, j) );
            YeUV.assignre(i, j, getSMEFTCoeffEW("YeR", i, j) );
            YeUV.assignim(i, j, getSMEFTCoeffEW("YeI", i, j) );
        }
    }

    YuUVhc = YuUV.hconjugate();
    YdUVhc = YdUV.hconjugate();
    YeUVhc = YeUV.hconjugate();
    
    
    // Matching

    CG_LNP = c3G_LNP;
    CW_LNP = c3W_LNP;
    CH_LNP = -c6_LNP - 2*cHW_LNP*g2UV*lambdaHUV - 4*cW_LNP*g2UV*lambdaHUV - 2*c2W_LNP*g2UV2*lambdaHUV;
    CHbox_LNP = -cH_LNP - cT_LNP - (cB_LNP*g1UV)/2. - (cHB_LNP*g1UV)/4. - (c2B_LNP*g1UV2)/4. - (3*cHW_LNP*g2UV)/4. - (3*cW_LNP*g2UV)/2. - (3*c2W_LNP*g2UV2)/4.;
    CHD_LNP = -4*cT_LNP - 2*cB_LNP*g1UV - cHB_LNP*g1UV - c2B_LNP*g1UV2;
    CHG_LNP = cg_LNP;
    CHB_LNP = cgam_LNP + (cHB_LNP*g1UV)/4.;
    CHW_LNP = (cHW_LNP*g2UV)/4.;
    CHWB_LNP = (cHW_LNP*g1UV)/4. + (cHB_LNP*g2UV)/4.;
    CuH_11r_LNP = ( cu_LNP*YuUVhc(0,0) - (cHW_LNP*g2UV*YuUVhc(0,0))/2. - cW_LNP*g2UV*YuUVhc(0,0) - (c2W_LNP*g2UV2*YuUVhc(0,0))/2. ).real();
    CuH_11i_LNP = ( cu_LNP*YuUVhc(0,0) - (cHW_LNP*g2UV*YuUVhc(0,0))/2. - cW_LNP*g2UV*YuUVhc(0,0) - (c2W_LNP*g2UV2*YuUVhc(0,0))/2. ).imag();
    CuH_12r_LNP = ( cu_LNP*YuUVhc(0,1) - (cHW_LNP*g2UV*YuUVhc(0,1))/2. - cW_LNP*g2UV*YuUVhc(0,1) - (c2W_LNP*g2UV2*YuUVhc(0,1))/2. ).real();
    CuH_12i_LNP = ( cu_LNP*YuUVhc(0,1) - (cHW_LNP*g2UV*YuUVhc(0,1))/2. - cW_LNP*g2UV*YuUVhc(0,1) - (c2W_LNP*g2UV2*YuUVhc(0,1))/2. ).imag();
    CuH_13r_LNP = ( cu_LNP*YuUVhc(0,2) - (cHW_LNP*g2UV*YuUVhc(0,2))/2. - cW_LNP*g2UV*YuUVhc(0,2) - (c2W_LNP*g2UV2*YuUVhc(0,2))/2. ).real();
    CuH_13i_LNP = ( cu_LNP*YuUVhc(0,2) - (cHW_LNP*g2UV*YuUVhc(0,2))/2. - cW_LNP*g2UV*YuUVhc(0,2) - (c2W_LNP*g2UV2*YuUVhc(0,2))/2. ).imag();
    CuH_21r_LNP = ( cu_LNP*YuUVhc(1,0) - (cHW_LNP*g2UV*YuUVhc(1,0))/2. - cW_LNP*g2UV*YuUVhc(1,0) - (c2W_LNP*g2UV2*YuUVhc(1,0))/2. ).real();
    CuH_21i_LNP = ( cu_LNP*YuUVhc(1,0) - (cHW_LNP*g2UV*YuUVhc(1,0))/2. - cW_LNP*g2UV*YuUVhc(1,0) - (c2W_LNP*g2UV2*YuUVhc(1,0))/2. ).imag();
    CuH_22r_LNP = ( cu_LNP*YuUVhc(1,1) - (cHW_LNP*g2UV*YuUVhc(1,1))/2. - cW_LNP*g2UV*YuUVhc(1,1) - (c2W_LNP*g2UV2*YuUVhc(1,1))/2. ).real();
    CuH_22i_LNP = ( cu_LNP*YuUVhc(1,1) - (cHW_LNP*g2UV*YuUVhc(1,1))/2. - cW_LNP*g2UV*YuUVhc(1,1) - (c2W_LNP*g2UV2*YuUVhc(1,1))/2. ).imag();
    CuH_23r_LNP = ( cu_LNP*YuUVhc(1,2) - (cHW_LNP*g2UV*YuUVhc(1,2))/2. - cW_LNP*g2UV*YuUVhc(1,2) - (c2W_LNP*g2UV2*YuUVhc(1,2))/2. ).real();
    CuH_23i_LNP = ( cu_LNP*YuUVhc(1,2) - (cHW_LNP*g2UV*YuUVhc(1,2))/2. - cW_LNP*g2UV*YuUVhc(1,2) - (c2W_LNP*g2UV2*YuUVhc(1,2))/2. ).imag();
    CuH_31r_LNP = ( cu_LNP*YuUVhc(2,0) - (cHW_LNP*g2UV*YuUVhc(2,0))/2. - cW_LNP*g2UV*YuUVhc(2,0) - (c2W_LNP*g2UV2*YuUVhc(2,0))/2. ).real();
    CuH_31i_LNP = ( cu_LNP*YuUVhc(2,0) - (cHW_LNP*g2UV*YuUVhc(2,0))/2. - cW_LNP*g2UV*YuUVhc(2,0) - (c2W_LNP*g2UV2*YuUVhc(2,0))/2. ).imag();
    CuH_32r_LNP = ( cu_LNP*YuUVhc(2,1) - (cHW_LNP*g2UV*YuUVhc(2,1))/2. - cW_LNP*g2UV*YuUVhc(2,1) - (c2W_LNP*g2UV2*YuUVhc(2,1))/2. ).real();
    CuH_32i_LNP = ( cu_LNP*YuUVhc(2,1) - (cHW_LNP*g2UV*YuUVhc(2,1))/2. - cW_LNP*g2UV*YuUVhc(2,1) - (c2W_LNP*g2UV2*YuUVhc(2,1))/2. ).imag();
    CuH_33r_LNP = ( cu_LNP*YuUVhc(2,2) - (cHW_LNP*g2UV*YuUVhc(2,2))/2. - cW_LNP*g2UV*YuUVhc(2,2) - (c2W_LNP*g2UV2*YuUVhc(2,2))/2. ).real();
    CuH_33i_LNP = ( cu_LNP*YuUVhc(2,2) - (cHW_LNP*g2UV*YuUVhc(2,2))/2. - cW_LNP*g2UV*YuUVhc(2,2) - (c2W_LNP*g2UV2*YuUVhc(2,2))/2. ).imag();
    CdH_11r_LNP = ( cd_LNP*YdUVhc(0,0) - (cHW_LNP*g2UV*YdUVhc(0,0))/2. - cW_LNP*g2UV*YdUVhc(0,0) - (c2W_LNP*g2UV2*YdUVhc(0,0))/2. ).real();
    CdH_11i_LNP = ( cd_LNP*YdUVhc(0,0) - (cHW_LNP*g2UV*YdUVhc(0,0))/2. - cW_LNP*g2UV*YdUVhc(0,0) - (c2W_LNP*g2UV2*YdUVhc(0,0))/2. ).imag();
    CdH_12r_LNP = ( cd_LNP*YdUVhc(0,1) - (cHW_LNP*g2UV*YdUVhc(0,1))/2. - cW_LNP*g2UV*YdUVhc(0,1) - (c2W_LNP*g2UV2*YdUVhc(0,1))/2. ).real();
    CdH_12i_LNP = ( cd_LNP*YdUVhc(0,1) - (cHW_LNP*g2UV*YdUVhc(0,1))/2. - cW_LNP*g2UV*YdUVhc(0,1) - (c2W_LNP*g2UV2*YdUVhc(0,1))/2. ).imag();
    CdH_13r_LNP = ( cd_LNP*YdUVhc(0,2) - (cHW_LNP*g2UV*YdUVhc(0,2))/2. - cW_LNP*g2UV*YdUVhc(0,2) - (c2W_LNP*g2UV2*YdUVhc(0,2))/2. ).real();
    CdH_13i_LNP = ( cd_LNP*YdUVhc(0,2) - (cHW_LNP*g2UV*YdUVhc(0,2))/2. - cW_LNP*g2UV*YdUVhc(0,2) - (c2W_LNP*g2UV2*YdUVhc(0,2))/2. ).imag();
    CdH_21r_LNP = ( cd_LNP*YdUVhc(1,0) - (cHW_LNP*g2UV*YdUVhc(1,0))/2. - cW_LNP*g2UV*YdUVhc(1,0) - (c2W_LNP*g2UV2*YdUVhc(1,0))/2. ).real();
    CdH_21i_LNP = ( cd_LNP*YdUVhc(1,0) - (cHW_LNP*g2UV*YdUVhc(1,0))/2. - cW_LNP*g2UV*YdUVhc(1,0) - (c2W_LNP*g2UV2*YdUVhc(1,0))/2. ).imag();
    CdH_22r_LNP = ( cd_LNP*YdUVhc(1,1) - (cHW_LNP*g2UV*YdUVhc(1,1))/2. - cW_LNP*g2UV*YdUVhc(1,1) - (c2W_LNP*g2UV2*YdUVhc(1,1))/2. ).real();
    CdH_22i_LNP = ( cd_LNP*YdUVhc(1,1) - (cHW_LNP*g2UV*YdUVhc(1,1))/2. - cW_LNP*g2UV*YdUVhc(1,1) - (c2W_LNP*g2UV2*YdUVhc(1,1))/2. ).imag();
    CdH_23r_LNP = ( cd_LNP*YdUVhc(1,2) - (cHW_LNP*g2UV*YdUVhc(1,2))/2. - cW_LNP*g2UV*YdUVhc(1,2) - (c2W_LNP*g2UV2*YdUVhc(1,2))/2. ).real();
    CdH_23i_LNP = ( cd_LNP*YdUVhc(1,2) - (cHW_LNP*g2UV*YdUVhc(1,2))/2. - cW_LNP*g2UV*YdUVhc(1,2) - (c2W_LNP*g2UV2*YdUVhc(1,2))/2. ).imag();
    CdH_31r_LNP = ( cd_LNP*YdUVhc(2,0) - (cHW_LNP*g2UV*YdUVhc(2,0))/2. - cW_LNP*g2UV*YdUVhc(2,0) - (c2W_LNP*g2UV2*YdUVhc(2,0))/2. ).real();
    CdH_31i_LNP = ( cd_LNP*YdUVhc(2,0) - (cHW_LNP*g2UV*YdUVhc(2,0))/2. - cW_LNP*g2UV*YdUVhc(2,0) - (c2W_LNP*g2UV2*YdUVhc(2,0))/2. ).imag();
    CdH_32r_LNP = ( cd_LNP*YdUVhc(2,1) - (cHW_LNP*g2UV*YdUVhc(2,1))/2. - cW_LNP*g2UV*YdUVhc(2,1) - (c2W_LNP*g2UV2*YdUVhc(2,1))/2. ).real();
    CdH_32i_LNP = ( cd_LNP*YdUVhc(2,1) - (cHW_LNP*g2UV*YdUVhc(2,1))/2. - cW_LNP*g2UV*YdUVhc(2,1) - (c2W_LNP*g2UV2*YdUVhc(2,1))/2. ).imag();
    CdH_33r_LNP = ( cd_LNP*YdUVhc(2,2) - (cHW_LNP*g2UV*YdUVhc(2,2))/2. - cW_LNP*g2UV*YdUVhc(2,2) - (c2W_LNP*g2UV2*YdUVhc(2,2))/2. ).real();
    CdH_33i_LNP = ( cd_LNP*YdUVhc(2,2) - (cHW_LNP*g2UV*YdUVhc(2,2))/2. - cW_LNP*g2UV*YdUVhc(2,2) - (c2W_LNP*g2UV2*YdUVhc(2,2))/2. ).imag();
    CeH_11r_LNP = ( ce_LNP*YeUVhc(0,0) - (cHW_LNP*g2UV*YeUVhc(0,0))/2. - cW_LNP*g2UV*YeUVhc(0,0) - (c2W_LNP*g2UV2*YeUVhc(0,0))/2. ).real();
    CeH_11i_LNP = ( ce_LNP*YeUVhc(0,0) - (cHW_LNP*g2UV*YeUVhc(0,0))/2. - cW_LNP*g2UV*YeUVhc(0,0) - (c2W_LNP*g2UV2*YeUVhc(0,0))/2. ).imag();
    CeH_12r_LNP = ( ce_LNP*YeUVhc(0,1) - (cHW_LNP*g2UV*YeUVhc(0,1))/2. - cW_LNP*g2UV*YeUVhc(0,1) - (c2W_LNP*g2UV2*YeUVhc(0,1))/2. ).real();
    CeH_12i_LNP = ( ce_LNP*YeUVhc(0,1) - (cHW_LNP*g2UV*YeUVhc(0,1))/2. - cW_LNP*g2UV*YeUVhc(0,1) - (c2W_LNP*g2UV2*YeUVhc(0,1))/2. ).imag();
    CeH_13r_LNP = ( ce_LNP*YeUVhc(0,2) - (cHW_LNP*g2UV*YeUVhc(0,2))/2. - cW_LNP*g2UV*YeUVhc(0,2) - (c2W_LNP*g2UV2*YeUVhc(0,2))/2. ).real();
    CeH_13i_LNP = ( ce_LNP*YeUVhc(0,2) - (cHW_LNP*g2UV*YeUVhc(0,2))/2. - cW_LNP*g2UV*YeUVhc(0,2) - (c2W_LNP*g2UV2*YeUVhc(0,2))/2. ).imag();
    CeH_21r_LNP = ( ce_LNP*YeUVhc(1,0) - (cHW_LNP*g2UV*YeUVhc(1,0))/2. - cW_LNP*g2UV*YeUVhc(1,0) - (c2W_LNP*g2UV2*YeUVhc(1,0))/2. ).real();
    CeH_21i_LNP = ( ce_LNP*YeUVhc(1,0) - (cHW_LNP*g2UV*YeUVhc(1,0))/2. - cW_LNP*g2UV*YeUVhc(1,0) - (c2W_LNP*g2UV2*YeUVhc(1,0))/2. ).imag();
    CeH_22r_LNP = ( ce_LNP*YeUVhc(1,1) - (cHW_LNP*g2UV*YeUVhc(1,1))/2. - cW_LNP*g2UV*YeUVhc(1,1) - (c2W_LNP*g2UV2*YeUVhc(1,1))/2. ).real();
    CeH_22i_LNP = ( ce_LNP*YeUVhc(1,1) - (cHW_LNP*g2UV*YeUVhc(1,1))/2. - cW_LNP*g2UV*YeUVhc(1,1) - (c2W_LNP*g2UV2*YeUVhc(1,1))/2. ).imag();
    CeH_23r_LNP = ( ce_LNP*YeUVhc(1,2) - (cHW_LNP*g2UV*YeUVhc(1,2))/2. - cW_LNP*g2UV*YeUVhc(1,2) - (c2W_LNP*g2UV2*YeUVhc(1,2))/2. ).real();
    CeH_23i_LNP = ( ce_LNP*YeUVhc(1,2) - (cHW_LNP*g2UV*YeUVhc(1,2))/2. - cW_LNP*g2UV*YeUVhc(1,2) - (c2W_LNP*g2UV2*YeUVhc(1,2))/2. ).imag();
    CeH_31r_LNP = ( ce_LNP*YeUVhc(2,0) - (cHW_LNP*g2UV*YeUVhc(2,0))/2. - cW_LNP*g2UV*YeUVhc(2,0) - (c2W_LNP*g2UV2*YeUVhc(2,0))/2. ).real();
    CeH_31i_LNP = ( ce_LNP*YeUVhc(2,0) - (cHW_LNP*g2UV*YeUVhc(2,0))/2. - cW_LNP*g2UV*YeUVhc(2,0) - (c2W_LNP*g2UV2*YeUVhc(2,0))/2. ).imag();
    CeH_32r_LNP = ( ce_LNP*YeUVhc(2,1) - (cHW_LNP*g2UV*YeUVhc(2,1))/2. - cW_LNP*g2UV*YeUVhc(2,1) - (c2W_LNP*g2UV2*YeUVhc(2,1))/2. ).real();
    CeH_32i_LNP = ( ce_LNP*YeUVhc(2,1) - (cHW_LNP*g2UV*YeUVhc(2,1))/2. - cW_LNP*g2UV*YeUVhc(2,1) - (c2W_LNP*g2UV2*YeUVhc(2,1))/2. ).imag();
    CeH_33r_LNP = ( ce_LNP*YeUVhc(2,2) - (cHW_LNP*g2UV*YeUVhc(2,2))/2. - cW_LNP*g2UV*YeUVhc(2,2) - (c2W_LNP*g2UV2*YeUVhc(2,2))/2. ).real();
    CeH_33i_LNP = ( ce_LNP*YeUVhc(2,2) - (cHW_LNP*g2UV*YeUVhc(2,2))/2. - cW_LNP*g2UV*YeUVhc(2,2) - (c2W_LNP*g2UV2*YeUVhc(2,2))/2. ).imag();
    CuG_33r_LNP = ctG_LNP;
    CuG_33i_LNP = 0.0;
    CuW_33r_LNP = ctW_LNP;
    CuW_33i_LNP = 0.0;
    CuB_33r_LNP = ctB_LNP;
    CuB_33i_LNP = 0.0;
    CHl1_11r_LNP = (cB_LNP*g1UV)/2. + (cHB_LNP*g1UV)/4. + (c2B_LNP*g1UV2)/2.;
    CHl1_22r_LNP = (cB_LNP*g1UV)/2. + (cHB_LNP*g1UV)/4. + (c2B_LNP*g1UV2)/2.;
    CHl1_33r_LNP = (cB_LNP*g1UV)/2. + (cHB_LNP*g1UV)/4. + (c2B_LNP*g1UV2)/2.;
    CHl3_11r_LNP = -0.25*(cHW_LNP*g2UV) - (cW_LNP*g2UV)/2. - (c2W_LNP*g2UV2)/2.;
    CHl3_22r_LNP = -0.25*(cHW_LNP*g2UV) - (cW_LNP*g2UV)/2. - (c2W_LNP*g2UV2)/2.;
    CHl3_33r_LNP = -0.25*(cHW_LNP*g2UV) - (cW_LNP*g2UV)/2. - (c2W_LNP*g2UV2)/2.;
    CHe_11r_LNP = cB_LNP*g1UV + (cHB_LNP*g1UV)/2. + c2B_LNP*g1UV2;
    CHe_22r_LNP = cB_LNP*g1UV + (cHB_LNP*g1UV)/2. + c2B_LNP*g1UV2;
    CHe_33r_LNP = cB_LNP*g1UV + (cHB_LNP*g1UV)/2. + c2B_LNP*g1UV2;
    CHq1_11r_LNP = -0.16666666666666666*(cB_LNP*g1UV) - (cHB_LNP*g1UV)/12. - (c2B_LNP*g1UV2)/6.;
    CHq1_22r_LNP = -0.16666666666666666*(cB_LNP*g1UV) - (cHB_LNP*g1UV)/12. - (c2B_LNP*g1UV2)/6.;
    CHq1_33r_LNP = cHq1_LNP - (cB_LNP*g1UV)/6. - (cHB_LNP*g1UV)/12. - (cqD1_LNP*g1UV)/2. - (c2B_LNP*g1UV2)/6.;
    CHq3_11r_LNP = -0.25*(cHW_LNP*g2UV) - (cW_LNP*g2UV)/2. - (c2W_LNP*g2UV2)/2.;
    CHq3_22r_LNP = -0.25*(cHW_LNP*g2UV) - (cW_LNP*g2UV)/2. - (c2W_LNP*g2UV2)/2.;
    CHq3_33r_LNP = cHq3_LNP - (cHW_LNP*g2UV)/4. - (cqD3_LNP*g2UV)/2. - (cW_LNP*g2UV)/2. - (c2W_LNP*g2UV2)/2.;
    CHu_11r_LNP = (-2*cB_LNP*g1UV)/3. - (cHB_LNP*g1UV)/3. - (2*c2B_LNP*g1UV2)/3.;
    CHu_22r_LNP = (-2*cB_LNP*g1UV)/3. - (cHB_LNP*g1UV)/3. - (2*c2B_LNP*g1UV2)/3.;
    CHu_33r_LNP = cHt_LNP - (2*cB_LNP*g1UV)/3. - (cHB_LNP*g1UV)/3. - (ctD_LNP*g1UV)/2. - (2*c2B_LNP*g1UV2)/3.;
    CHd_11r_LNP = (cB_LNP*g1UV)/3. + (cHB_LNP*g1UV)/6. + (c2B_LNP*g1UV2)/3.;
    CHd_22r_LNP = (cB_LNP*g1UV)/3. + (cHB_LNP*g1UV)/6. + (c2B_LNP*g1UV2)/3.;
    CHd_33r_LNP = (cB_LNP*g1UV)/3. + (cHB_LNP*g1UV)/6. + (c2B_LNP*g1UV2)/3.;
    Cll_1111r_LNP = -0.25*(c2B_LNP*g1UV2) - (c2W_LNP*g2UV2)/4.;
    Cll_1122r_LNP = -0.25*(c2B_LNP*g1UV2) + (c2W_LNP*g2UV2)/4.;
    Cll_1133r_LNP = -0.25*(c2B_LNP*g1UV2) + (c2W_LNP*g2UV2)/4.;
    Cll_1221r_LNP = -0.5*(c2W_LNP*g2UV2);
    Cll_1331r_LNP = -0.5*(c2W_LNP*g2UV2);
    Cll_2222r_LNP = -0.25*(c2B_LNP*g1UV2) - (c2W_LNP*g2UV2)/4.;
    Cll_2233r_LNP = -0.25*(c2B_LNP*g1UV2) + (c2W_LNP*g2UV2)/4.;
    Cll_2332r_LNP = -0.5*(c2W_LNP*g2UV2);
    Cll_3333r_LNP = -0.25*(c2B_LNP*g1UV2) - (c2W_LNP*g2UV2)/4.;
    Cqq1_1111r_LNP = -0.027777777777777776*(c2B_LNP*g1UV2) - (c2G_LNP*g3UV2)/12.;
    Cqq1_1122r_LNP = -0.027777777777777776*(c2B_LNP*g1UV2) + (c2G_LNP*g3UV2)/6.;
    Cqq1_1133r_LNP = -0.16666666666666666*(cqD1_LNP*g1UV) - (c2B_LNP*g1UV2)/36. + (c2G_LNP*g3UV2)/6.;
    Cqq1_1221r_LNP = -0.25*(c2G_LNP*g3UV2);
    Cqq1_1331r_LNP = -0.25*(c2G_LNP*g3UV2);
    Cqq1_2222r_LNP = -0.027777777777777776*(c2B_LNP*g1UV2) - (c2G_LNP*g3UV2)/12.;
    Cqq1_2233r_LNP = -0.16666666666666666*(cqD1_LNP*g1UV) - (c2B_LNP*g1UV2)/36. + (c2G_LNP*g3UV2)/6.;
    Cqq1_2332r_LNP = -0.25*(c2G_LNP*g3UV2);
    Cqq1_3333r_LNP = cqq1_LNP - (cqD1_LNP*g1UV)/6. - (c2B_LNP*g1UV2)/36. - (c2G_LNP*g3UV2)/12.;
    Cqq3_1111r_LNP = -(c2W_LNP*g2UV2) - (c2G_LNP*g3UV2)/4.;
    Cqq3_1122r_LNP = -(c2W_LNP*g2UV2);
    Cqq3_1133r_LNP = -0.5*(cqD3_LNP*g2UV) - c2W_LNP*g2UV2;
    Cqq3_1221r_LNP = -0.25*(c2G_LNP*g3UV2);
    Cqq3_1331r_LNP = -0.25*(c2G_LNP*g3UV2);
    Cqq3_2222r_LNP = -(c2W_LNP*g2UV2) - (c2G_LNP*g3UV2)/4.;
    Cqq3_2233r_LNP = -0.5*(cqD3_LNP*g2UV) - c2W_LNP*g2UV2;
    Cqq3_2332r_LNP = -0.25*(c2G_LNP*g3UV2);
    Cqq3_3333r_LNP = cqq3_LNP - (cqD3_LNP*g2UV)/2. - c2W_LNP*g2UV2 - (c2G_LNP*g3UV2)/4.;
    Clq1_1111r_LNP = (c2B_LNP*g1UV2)/6.;
    Clq1_1122r_LNP = (c2B_LNP*g1UV2)/6.;
    Clq1_1133r_LNP = (c2B_LNP*g1UV2)/6.;
    Clq1_2211r_LNP = (c2B_LNP*g1UV2)/6.;
    Clq1_2222r_LNP = (c2B_LNP*g1UV2)/6.;
    Clq1_2233r_LNP = (c2B_LNP*g1UV2)/6.;
    Clq1_3311r_LNP = (c2B_LNP*g1UV2)/6.;
    Clq1_3322r_LNP = (c2B_LNP*g1UV2)/6.;
    Clq1_3333r_LNP = (c2B_LNP*g1UV2)/6.;
    Clq3_1111r_LNP = -0.5*(c2W_LNP*g2UV2);
    Clq3_1122r_LNP = -0.5*(c2W_LNP*g2UV2);
    Clq3_1133r_LNP = -0.5*(c2W_LNP*g2UV2);
    Clq3_2211r_LNP = -0.5*(c2W_LNP*g2UV2);
    Clq3_2222r_LNP = -0.5*(c2W_LNP*g2UV2);
    Clq3_2233r_LNP = -0.5*(c2W_LNP*g2UV2);
    Clq3_3311r_LNP = -0.5*(c2W_LNP*g2UV2);
    Clq3_3322r_LNP = -0.5*(c2W_LNP*g2UV2);
    Clq3_3333r_LNP = -0.5*(c2W_LNP*g2UV2);
    Cee_1111r_LNP = -(c2B_LNP*g1UV2);
    Cee_1122r_LNP = -0.5*(c2B_LNP*g1UV2);
    Cee_1133r_LNP = -0.5*(c2B_LNP*g1UV2);
    Cee_2222r_LNP = -(c2B_LNP*g1UV2);
    Cee_2233r_LNP = -0.5*(c2B_LNP*g1UV2);
    Cee_3333r_LNP = -(c2B_LNP*g1UV2);
    Cuu_1111r_LNP = (-4*c2B_LNP*g1UV2)/9. - (c2G_LNP*g3UV2)/3.;
    Cuu_1122r_LNP = (-4*c2B_LNP*g1UV2)/9. + (c2G_LNP*g3UV2)/6.;
    Cuu_1133r_LNP = (-2*ctD_LNP*g1UV)/3. - (4*c2B_LNP*g1UV2)/9. + (c2G_LNP*g3UV2)/6.;
    Cuu_1221r_LNP = -0.5*(c2G_LNP*g3UV2);
    Cuu_1331r_LNP = -0.5*(c2G_LNP*g3UV2);
    Cuu_2222r_LNP = (-4*c2B_LNP*g1UV2)/9. - (c2G_LNP*g3UV2)/3.;
    Cuu_2233r_LNP = (-2*ctD_LNP*g1UV)/3. - (4*c2B_LNP*g1UV2)/9. + (c2G_LNP*g3UV2)/6.;
    Cuu_2332r_LNP = -0.5*(c2G_LNP*g3UV2);
    Cuu_3333r_LNP = ctt_LNP - (2*ctD_LNP*g1UV)/3. - (4*c2B_LNP*g1UV2)/9. - (c2G_LNP*g3UV2)/3.;
    Cdd_1111r_LNP = -0.1111111111111111*(c2B_LNP*g1UV2) - (c2G_LNP*g3UV2)/3.;
    Cdd_1122r_LNP = -0.1111111111111111*(c2B_LNP*g1UV2) + (c2G_LNP*g3UV2)/6.;
    Cdd_1133r_LNP = -0.1111111111111111*(c2B_LNP*g1UV2) + (c2G_LNP*g3UV2)/6.;
    Cdd_1221r_LNP = -0.5*(c2G_LNP*g3UV2);
    Cdd_1331r_LNP = -0.5*(c2G_LNP*g3UV2);
    Cdd_2222r_LNP = -0.1111111111111111*(c2B_LNP*g1UV2) - (c2G_LNP*g3UV2)/3.;
    Cdd_2233r_LNP = -0.1111111111111111*(c2B_LNP*g1UV2) + (c2G_LNP*g3UV2)/6.;
    Cdd_2332r_LNP = -0.5*(c2G_LNP*g3UV2);
    Cdd_3333r_LNP = -0.1111111111111111*(c2B_LNP*g1UV2) - (c2G_LNP*g3UV2)/3.;
    Ceu_1111r_LNP = (4*c2B_LNP*g1UV2)/3.;
    Ceu_1122r_LNP = (4*c2B_LNP*g1UV2)/3.;
    Ceu_1133r_LNP = ctD_LNP*g1UV + (4*c2B_LNP*g1UV2)/3.;
    Ceu_2211r_LNP = (4*c2B_LNP*g1UV2)/3.;
    Ceu_2222r_LNP = (4*c2B_LNP*g1UV2)/3.;
    Ceu_2233r_LNP = ctD_LNP*g1UV + (4*c2B_LNP*g1UV2)/3.;
    Ceu_3311r_LNP = (4*c2B_LNP*g1UV2)/3.;
    Ceu_3322r_LNP = (4*c2B_LNP*g1UV2)/3.;
    Ceu_3333r_LNP = ctD_LNP*g1UV + (4*c2B_LNP*g1UV2)/3.;
    Ced_1111r_LNP = (-2*c2B_LNP*g1UV2)/3.;
    Ced_1122r_LNP = (-2*c2B_LNP*g1UV2)/3.;
    Ced_1133r_LNP = (-2*c2B_LNP*g1UV2)/3.;
    Ced_2211r_LNP = (-2*c2B_LNP*g1UV2)/3.;
    Ced_2222r_LNP = (-2*c2B_LNP*g1UV2)/3.;
    Ced_2233r_LNP = (-2*c2B_LNP*g1UV2)/3.;
    Ced_3311r_LNP = (-2*c2B_LNP*g1UV2)/3.;
    Ced_3322r_LNP = (-2*c2B_LNP*g1UV2)/3.;
    Ced_3333r_LNP = (-2*c2B_LNP*g1UV2)/3.;
    Cud1_1111r_LNP = (4*c2B_LNP*g1UV2)/9.;
    Cud1_1122r_LNP = (4*c2B_LNP*g1UV2)/9.;
    Cud1_1133r_LNP = (4*c2B_LNP*g1UV2)/9.;
    Cud1_2211r_LNP = (4*c2B_LNP*g1UV2)/9.;
    Cud1_2222r_LNP = (4*c2B_LNP*g1UV2)/9.;
    Cud1_2233r_LNP = (4*c2B_LNP*g1UV2)/9.;
    Cud1_3311r_LNP = (ctD_LNP*g1UV)/3. + (4*c2B_LNP*g1UV2)/9.;
    Cud1_3322r_LNP = (ctD_LNP*g1UV)/3. + (4*c2B_LNP*g1UV2)/9.;
    Cud1_3333r_LNP = (ctD_LNP*g1UV)/3. + (4*c2B_LNP*g1UV2)/9.;
    Cud8_1111r_LNP = -2*c2G_LNP*g3UV2;
    Cud8_1122r_LNP = -2*c2G_LNP*g3UV2;
    Cud8_1133r_LNP = -2*c2G_LNP*g3UV2;
    Cud8_2211r_LNP = -2*c2G_LNP*g3UV2;
    Cud8_2222r_LNP = -2*c2G_LNP*g3UV2;
    Cud8_2233r_LNP = -2*c2G_LNP*g3UV2;
    Cud8_3311r_LNP = -2*c2G_LNP*g3UV2;
    Cud8_3322r_LNP = -2*c2G_LNP*g3UV2;
    Cud8_3333r_LNP = -2*c2G_LNP*g3UV2;
    Cle_1111r_LNP = -(c2B_LNP*g1UV2);
    Cle_1122r_LNP = -(c2B_LNP*g1UV2);
    Cle_1133r_LNP = -(c2B_LNP*g1UV2);
    Cle_2211r_LNP = -(c2B_LNP*g1UV2);
    Cle_2222r_LNP = -(c2B_LNP*g1UV2);
    Cle_2233r_LNP = -(c2B_LNP*g1UV2);
    Cle_3311r_LNP = -(c2B_LNP*g1UV2);
    Cle_3322r_LNP = -(c2B_LNP*g1UV2);
    Cle_3333r_LNP = -(c2B_LNP*g1UV2);
    Clu_1111r_LNP = (2*c2B_LNP*g1UV2)/3.;
    Clu_1122r_LNP = (2*c2B_LNP*g1UV2)/3.;
    Clu_1133r_LNP = (ctD_LNP*g1UV)/2. + (2*c2B_LNP*g1UV2)/3.;
    Clu_2211r_LNP = (2*c2B_LNP*g1UV2)/3.;
    Clu_2222r_LNP = (2*c2B_LNP*g1UV2)/3.;
    Clu_2233r_LNP = (ctD_LNP*g1UV)/2. + (2*c2B_LNP*g1UV2)/3.;
    Clu_3311r_LNP = (2*c2B_LNP*g1UV2)/3.;
    Clu_3322r_LNP = (2*c2B_LNP*g1UV2)/3.;
    Clu_3333r_LNP = (ctD_LNP*g1UV)/2. + (2*c2B_LNP*g1UV2)/3.;
    Cld_1111r_LNP = -0.3333333333333333*(c2B_LNP*g1UV2);
    Cld_1122r_LNP = -0.3333333333333333*(c2B_LNP*g1UV2);
    Cld_1133r_LNP = -0.3333333333333333*(c2B_LNP*g1UV2);
    Cld_2211r_LNP = -0.3333333333333333*(c2B_LNP*g1UV2);
    Cld_2222r_LNP = -0.3333333333333333*(c2B_LNP*g1UV2);
    Cld_2233r_LNP = -0.3333333333333333*(c2B_LNP*g1UV2);
    Cld_3311r_LNP = -0.3333333333333333*(c2B_LNP*g1UV2);
    Cld_3322r_LNP = -0.3333333333333333*(c2B_LNP*g1UV2);
    Cld_3333r_LNP = -0.3333333333333333*(c2B_LNP*g1UV2);
    Cqe_1111r_LNP = (c2B_LNP*g1UV2)/3.;
    Cqe_1122r_LNP = (c2B_LNP*g1UV2)/3.;
    Cqe_1133r_LNP = (c2B_LNP*g1UV2)/3.;
    Cqe_2211r_LNP = (c2B_LNP*g1UV2)/3.;
    Cqe_2222r_LNP = (c2B_LNP*g1UV2)/3.;
    Cqe_2233r_LNP = (c2B_LNP*g1UV2)/3.;
    Cqe_3311r_LNP = cqD1_LNP*g1UV + (c2B_LNP*g1UV2)/3.;
    Cqe_3322r_LNP = cqD1_LNP*g1UV + (c2B_LNP*g1UV2)/3.;
    Cqe_3333r_LNP = cqD1_LNP*g1UV + (c2B_LNP*g1UV2)/3.;
    Cqu1_1111r_LNP = (-2*c2B_LNP*g1UV2)/9.;
    Cqu1_1122r_LNP = (-2*c2B_LNP*g1UV2)/9.;
    Cqu1_1133r_LNP = -0.16666666666666666*(ctD_LNP*g1UV) - (2*c2B_LNP*g1UV2)/9.;
    Cqu1_2211r_LNP = (-2*c2B_LNP*g1UV2)/9.;
    Cqu1_2222r_LNP = (-2*c2B_LNP*g1UV2)/9.;
    Cqu1_2233r_LNP = -0.16666666666666666*(ctD_LNP*g1UV) - (2*c2B_LNP*g1UV2)/9.;
    Cqu1_3311r_LNP = (-2*cqD1_LNP*g1UV)/3. - (2*c2B_LNP*g1UV2)/9.;
    Cqu1_3322r_LNP = (-2*cqD1_LNP*g1UV)/3. - (2*c2B_LNP*g1UV2)/9.;
    Cqu1_3333r_LNP = cqt1_LNP - (2*cqD1_LNP*g1UV)/3. - (ctD_LNP*g1UV)/6. - (2*c2B_LNP*g1UV2)/9.;
    Cqu8_1111r_LNP = -2*c2G_LNP*g3UV2;
    Cqu8_1122r_LNP = -2*c2G_LNP*g3UV2;
    Cqu8_1133r_LNP = -2*c2G_LNP*g3UV2;
    Cqu8_2211r_LNP = -2*c2G_LNP*g3UV2;
    Cqu8_2222r_LNP = -2*c2G_LNP*g3UV2;
    Cqu8_2233r_LNP = -2*c2G_LNP*g3UV2;
    Cqu8_3311r_LNP = -2*c2G_LNP*g3UV2;
    Cqu8_3322r_LNP = -2*c2G_LNP*g3UV2;
    Cqu8_3333r_LNP = cqt8_LNP - 2*c2G_LNP*g3UV2;
    Cqd1_1111r_LNP = (c2B_LNP*g1UV2)/9.;
    Cqd1_1122r_LNP = (c2B_LNP*g1UV2)/9.;
    Cqd1_1133r_LNP = (c2B_LNP*g1UV2)/9.;
    Cqd1_2211r_LNP = (c2B_LNP*g1UV2)/9.;
    Cqd1_2222r_LNP = (c2B_LNP*g1UV2)/9.;
    Cqd1_2233r_LNP = (c2B_LNP*g1UV2)/9.;
    Cqd1_3311r_LNP = (cqD1_LNP*g1UV)/3. + (c2B_LNP*g1UV2)/9.;
    Cqd1_3322r_LNP = (cqD1_LNP*g1UV)/3. + (c2B_LNP*g1UV2)/9.;
    Cqd1_3333r_LNP = (cqD1_LNP*g1UV)/3. + (c2B_LNP*g1UV2)/9.;
    Cqd8_1111r_LNP = -2*c2G_LNP*g3UV2;
    Cqd8_1122r_LNP = -2*c2G_LNP*g3UV2;
    Cqd8_1133r_LNP = -2*c2G_LNP*g3UV2;
    Cqd8_2211r_LNP = -2*c2G_LNP*g3UV2;
    Cqd8_2222r_LNP = -2*c2G_LNP*g3UV2;
    Cqd8_2233r_LNP = -2*c2G_LNP*g3UV2;
    Cqd8_3311r_LNP = -2*c2G_LNP*g3UV2;
    Cqd8_3322r_LNP = -2*c2G_LNP*g3UV2;
    Cqd8_3333r_LNP = -2*c2G_LNP*g3UV2;
}


bool NPd6SILH::PostUpdate()
{    
    //  1) Post-update operations involving SM parameters only 

    v2LambdaNP2 = v() * v() / Lambda_NP/Lambda_NP;
    
//  Obtain the SM parameters at the UV scale
    GenerateSMInitialConditions();

//  Perform the matching
    setNPSMEFTd6GeneralParameters();

    if (!NPSMEFTd6General::PostUpdate()) return (false);

    return (true);

}

const double NPd6SILH::obliqueS() const
{
    return 0.0;
}

const double NPd6SILH::obliqueT() const
{
    return 0.0;
}

const double NPd6SILH::obliqueW() const
{
    return ( g2UV2 * c2W_LNP * v2LambdaNP2 / 2.0);
}

const double NPd6SILH::obliqueY() const
{
    return ( g2UV2 *  c2B_LNP * v2LambdaNP2 / 2.0);
}


const double NPd6SILH::AuxObs_NP30() const
{
    // To be used for some temporary observable

    // WY analysis at 13 TeV for HL-LHC 3/ab
    double Wpar, Ypar, Wpar2, Ypar2;
    double Chi2NC13, Chi2CC13, Chi2Tot;

    Wpar = 10000.0 * obliqueW();
    Ypar = 10000.0 * obliqueY();

    Wpar2 = Wpar*Wpar;
    Ypar2 = Ypar*Ypar;

    Chi2CC13 = Wpar2 * (18.365037149441695 + 2.422904241798858 * Wpar + 0.12120594308623695 * Wpar2);

    Chi2NC13 = 0.032772034538390675 * Wpar2 * Wpar2 + 2.815243944990361 * Ypar2 - 0.36522061776278516 * Ypar2 * Ypar
            + 0.017375258924241194 * Ypar2 * Ypar2 + Wpar2 * Wpar * (-0.7059117582389635 + 0.006816297425306027 * Ypar)
            + Wpar * Ypar * (7.988302197022343 + Ypar * (-0.5450119819316416 + 0.0050292149953719766 * Ypar))
            + Wpar2 * (5.68581760491364 + Ypar * (-0.5794111075840261 + 0.048026245835369625 * Ypar));

    Chi2Tot = Chi2CC13 + Chi2NC13;

    // To be used as Gaussian observable with mean=0, var=1 I must return the sqrt.
    return sqrt(Chi2Tot);
}