/* Copyright (C) 2025 HEPfit Collaboration
* 
* For the licensing terms see doc/COPYING.
*
* Created on Sun 27 Jul 2025 17:40:34
*
*/

#include "NPd6SILH.h"


std::string NPd6SILH::NPd6SILHVars[NNPd6SILHVars] = {
    "cH", "cT", "c6", "cB", "cW", "c2B", "c2W", "c2G", "c3W", "c3G", "cHW", "cHB", "cgam", "cg", "ctD", "cqD1", "cqD3", "cqq1", "cqq3", "cqt1", "cqt8", "ctt", "ctG", "ctB", "ctW", "cu", "cd", "ce", "Lambda_NP" 
};


NPd6SILH::NPd6SILH()
: NPSMEFTd6General(), YuUV(3,3,0.), YuUVhc(3,3,0.), YdUV(3,3,0.), YdUVhc(3,3,0.), YeUV(3,3,0.), YeUVhc(3,3,0.) {
    setModelName("NPd6SILH");

    FlagRGEci = true;

    ModelParamMap.insert(std::make_pair("cH", std::cref(cH)));
    ModelParamMap.insert(std::make_pair("cT", std::cref(cT)));
    ModelParamMap.insert(std::make_pair("c6", std::cref(c6)));
    ModelParamMap.insert(std::make_pair("cB", std::cref(cB)));
    ModelParamMap.insert(std::make_pair("cW", std::cref(cW)));
    ModelParamMap.insert(std::make_pair("c2B", std::cref(c2B)));
    ModelParamMap.insert(std::make_pair("c2W", std::cref(c2W)));
    ModelParamMap.insert(std::make_pair("c2G", std::cref(c2G)));
    ModelParamMap.insert(std::make_pair("c3W", std::cref(c3W)));
    ModelParamMap.insert(std::make_pair("c3G", std::cref(c3G)));
    ModelParamMap.insert(std::make_pair("cHW", std::cref(cHW)));
    ModelParamMap.insert(std::make_pair("cHB", std::cref(cHB)));
    ModelParamMap.insert(std::make_pair("cgam", std::cref(cgam)));
    ModelParamMap.insert(std::make_pair("cg", std::cref(cg)));
    ModelParamMap.insert(std::make_pair("ctD", std::cref(ctD)));
    ModelParamMap.insert(std::make_pair("cqD1", std::cref(cqD1)));
    ModelParamMap.insert(std::make_pair("cqD3", std::cref(cqD3)));
    ModelParamMap.insert(std::make_pair("cqq1", std::cref(cqq1)));
    ModelParamMap.insert(std::make_pair("cqq3", std::cref(cqq3)));
    ModelParamMap.insert(std::make_pair("cqt1", std::cref(cqt1)));
    ModelParamMap.insert(std::make_pair("cqt8", std::cref(cqt8)));
    ModelParamMap.insert(std::make_pair("ctt", std::cref(ctt)));
    ModelParamMap.insert(std::make_pair("ctG", std::cref(ctG)));
    ModelParamMap.insert(std::make_pair("ctB", std::cref(ctB)));
    ModelParamMap.insert(std::make_pair("ctW", std::cref(ctW)));
    ModelParamMap.insert(std::make_pair("cu", std::cref(cu)));
    ModelParamMap.insert(std::make_pair("cd", std::cref(cd)));
    ModelParamMap.insert(std::make_pair("ce", std::cref(ce)));

}


void NPd6SILH::setParameter(const std::string name, const double& value)
{

    if(name.compare("cH") == 0 ) { 

        cH = value;

    } else if(name.compare("cT")==0) {

        cT = value;

    } else if(name.compare("c6")==0) {

        c6 = value;

    } else if(name.compare("cB")==0) {

        cB = value;

    } else if(name.compare("cW")==0) {

        cW = value;

    } else if(name.compare("c2B")==0) {

        c2B = value;

    } else if(name.compare("c2W")==0) {

        c2W = value;

    } else if(name.compare("c2G")==0) {

        c2G = value;

    } else if(name.compare("c3W")==0) {

        c3W = value;

    } else if(name.compare("c3G")==0) {

        c3G = value;

    } else if(name.compare("cHW")==0) {

        cHW = value;

    } else if(name.compare("cHB")==0) {

        cHB = value;

    } else if(name.compare("cgam")==0) {

        cgam = value;

    } else if(name.compare("cg")==0) {

        cg = value;

    } else if(name.compare("ctD")==0) {

        ctD = value;

    } else if(name.compare("cqD1")==0) {

        cqD1 = value;

    } else if(name.compare("cqD3")==0) {

        cqD3 = value;

    } else if(name.compare("cqq1")==0) {

        cqq1 = value;

    } else if(name.compare("cqq3")==0) {

        cqq3 = value;

    } else if(name.compare("cqt1")==0) {

        cqt1 = value;

    } else if(name.compare("cqt8")==0) {

        cqt8 = value;

    } else if(name.compare("ctt")==0) {

        ctt = value;

    } else if(name.compare("ctG")==0) {

        ctG = value;

    } else if(name.compare("ctB")==0) {

        ctB = value;

    } else if(name.compare("ctW")==0) {

        ctW = value;

    } else if(name.compare("cu")==0) {

        cu = value;

    } else if(name.compare("cd")==0) {

        cd = value;

    } else if(name.compare("ce")==0) {

        ce = value;

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

    CG_LNP = c3G;
    CW_LNP = c3W;
    CH_LNP = -c6 - 2*cHW*g2UV*lambdaHUV - 4*cW*g2UV*lambdaHUV - 2*c2W*g2UV2*lambdaHUV;
    CHbox_LNP = -cH - cT - (cB*g1UV)/2. - (cHB*g1UV)/4. - (c2B*g1UV2)/4. - (3*cHW*g2UV)/4. - (3*cW*g2UV)/2. - (3*c2W*g2UV2)/4.;
    CHD_LNP = -4*cT - 2*cB*g1UV - cHB*g1UV - c2B*g1UV2;
    CHG_LNP = cg;
    CHB_LNP = cgam + (cHB*g1UV)/4.;
    CHW_LNP = (cHW*g2UV)/4.;
    CHWB_LNP = (cHW*g1UV)/4. + (cHB*g2UV)/4.;
    CuH_11r_LNP = ( cu*YuUVhc(0,0) - (cHW*g2UV*YuUVhc(0,0))/2. - cW*g2UV*YuUVhc(0,0) - (c2W*g2UV2*YuUVhc(0,0))/2. ).real();
    CuH_11i_LNP = ( cu*YuUVhc(0,0) - (cHW*g2UV*YuUVhc(0,0))/2. - cW*g2UV*YuUVhc(0,0) - (c2W*g2UV2*YuUVhc(0,0))/2. ).imag();
    CuH_12r_LNP = ( cu*YuUVhc(0,1) - (cHW*g2UV*YuUVhc(0,1))/2. - cW*g2UV*YuUVhc(0,1) - (c2W*g2UV2*YuUVhc(0,1))/2. ).real();
    CuH_12i_LNP = ( cu*YuUVhc(0,1) - (cHW*g2UV*YuUVhc(0,1))/2. - cW*g2UV*YuUVhc(0,1) - (c2W*g2UV2*YuUVhc(0,1))/2. ).imag();
    CuH_13r_LNP = ( cu*YuUVhc(0,2) - (cHW*g2UV*YuUVhc(0,2))/2. - cW*g2UV*YuUVhc(0,2) - (c2W*g2UV2*YuUVhc(0,2))/2. ).real();
    CuH_13i_LNP = ( cu*YuUVhc(0,2) - (cHW*g2UV*YuUVhc(0,2))/2. - cW*g2UV*YuUVhc(0,2) - (c2W*g2UV2*YuUVhc(0,2))/2. ).imag();
    CuH_21r_LNP = ( cu*YuUVhc(1,0) - (cHW*g2UV*YuUVhc(1,0))/2. - cW*g2UV*YuUVhc(1,0) - (c2W*g2UV2*YuUVhc(1,0))/2. ).real();
    CuH_21i_LNP = ( cu*YuUVhc(1,0) - (cHW*g2UV*YuUVhc(1,0))/2. - cW*g2UV*YuUVhc(1,0) - (c2W*g2UV2*YuUVhc(1,0))/2. ).imag();
    CuH_22r_LNP = ( cu*YuUVhc(1,1) - (cHW*g2UV*YuUVhc(1,1))/2. - cW*g2UV*YuUVhc(1,1) - (c2W*g2UV2*YuUVhc(1,1))/2. ).real();
    CuH_22i_LNP = ( cu*YuUVhc(1,1) - (cHW*g2UV*YuUVhc(1,1))/2. - cW*g2UV*YuUVhc(1,1) - (c2W*g2UV2*YuUVhc(1,1))/2. ).imag();
    CuH_23r_LNP = ( cu*YuUVhc(1,2) - (cHW*g2UV*YuUVhc(1,2))/2. - cW*g2UV*YuUVhc(1,2) - (c2W*g2UV2*YuUVhc(1,2))/2. ).real();
    CuH_23i_LNP = ( cu*YuUVhc(1,2) - (cHW*g2UV*YuUVhc(1,2))/2. - cW*g2UV*YuUVhc(1,2) - (c2W*g2UV2*YuUVhc(1,2))/2. ).imag();
    CuH_31r_LNP = ( cu*YuUVhc(2,0) - (cHW*g2UV*YuUVhc(2,0))/2. - cW*g2UV*YuUVhc(2,0) - (c2W*g2UV2*YuUVhc(2,0))/2. ).real();
    CuH_31i_LNP = ( cu*YuUVhc(2,0) - (cHW*g2UV*YuUVhc(2,0))/2. - cW*g2UV*YuUVhc(2,0) - (c2W*g2UV2*YuUVhc(2,0))/2. ).imag();
    CuH_32r_LNP = ( cu*YuUVhc(2,1) - (cHW*g2UV*YuUVhc(2,1))/2. - cW*g2UV*YuUVhc(2,1) - (c2W*g2UV2*YuUVhc(2,1))/2. ).real();
    CuH_32i_LNP = ( cu*YuUVhc(2,1) - (cHW*g2UV*YuUVhc(2,1))/2. - cW*g2UV*YuUVhc(2,1) - (c2W*g2UV2*YuUVhc(2,1))/2. ).imag();
    CuH_33r_LNP = ( cu*YuUVhc(2,2) - (cHW*g2UV*YuUVhc(2,2))/2. - cW*g2UV*YuUVhc(2,2) - (c2W*g2UV2*YuUVhc(2,2))/2. ).real();
    CuH_33i_LNP = ( cu*YuUVhc(2,2) - (cHW*g2UV*YuUVhc(2,2))/2. - cW*g2UV*YuUVhc(2,2) - (c2W*g2UV2*YuUVhc(2,2))/2. ).imag();
    CdH_11r_LNP = ( cd*YdUVhc(0,0) - (cHW*g2UV*YdUVhc(0,0))/2. - cW*g2UV*YdUVhc(0,0) - (c2W*g2UV2*YdUVhc(0,0))/2. ).real();
    CdH_11i_LNP = ( cd*YdUVhc(0,0) - (cHW*g2UV*YdUVhc(0,0))/2. - cW*g2UV*YdUVhc(0,0) - (c2W*g2UV2*YdUVhc(0,0))/2. ).imag();
    CdH_12r_LNP = ( cd*YdUVhc(0,1) - (cHW*g2UV*YdUVhc(0,1))/2. - cW*g2UV*YdUVhc(0,1) - (c2W*g2UV2*YdUVhc(0,1))/2. ).real();
    CdH_12i_LNP = ( cd*YdUVhc(0,1) - (cHW*g2UV*YdUVhc(0,1))/2. - cW*g2UV*YdUVhc(0,1) - (c2W*g2UV2*YdUVhc(0,1))/2. ).imag();
    CdH_13r_LNP = ( cd*YdUVhc(0,2) - (cHW*g2UV*YdUVhc(0,2))/2. - cW*g2UV*YdUVhc(0,2) - (c2W*g2UV2*YdUVhc(0,2))/2. ).real();
    CdH_13i_LNP = ( cd*YdUVhc(0,2) - (cHW*g2UV*YdUVhc(0,2))/2. - cW*g2UV*YdUVhc(0,2) - (c2W*g2UV2*YdUVhc(0,2))/2. ).imag();
    CdH_21r_LNP = ( cd*YdUVhc(1,0) - (cHW*g2UV*YdUVhc(1,0))/2. - cW*g2UV*YdUVhc(1,0) - (c2W*g2UV2*YdUVhc(1,0))/2. ).real();
    CdH_21i_LNP = ( cd*YdUVhc(1,0) - (cHW*g2UV*YdUVhc(1,0))/2. - cW*g2UV*YdUVhc(1,0) - (c2W*g2UV2*YdUVhc(1,0))/2. ).imag();
    CdH_22r_LNP = ( cd*YdUVhc(1,1) - (cHW*g2UV*YdUVhc(1,1))/2. - cW*g2UV*YdUVhc(1,1) - (c2W*g2UV2*YdUVhc(1,1))/2. ).real();
    CdH_22i_LNP = ( cd*YdUVhc(1,1) - (cHW*g2UV*YdUVhc(1,1))/2. - cW*g2UV*YdUVhc(1,1) - (c2W*g2UV2*YdUVhc(1,1))/2. ).imag();
    CdH_23r_LNP = ( cd*YdUVhc(1,2) - (cHW*g2UV*YdUVhc(1,2))/2. - cW*g2UV*YdUVhc(1,2) - (c2W*g2UV2*YdUVhc(1,2))/2. ).real();
    CdH_23i_LNP = ( cd*YdUVhc(1,2) - (cHW*g2UV*YdUVhc(1,2))/2. - cW*g2UV*YdUVhc(1,2) - (c2W*g2UV2*YdUVhc(1,2))/2. ).imag();
    CdH_31r_LNP = ( cd*YdUVhc(2,0) - (cHW*g2UV*YdUVhc(2,0))/2. - cW*g2UV*YdUVhc(2,0) - (c2W*g2UV2*YdUVhc(2,0))/2. ).real();
    CdH_31i_LNP = ( cd*YdUVhc(2,0) - (cHW*g2UV*YdUVhc(2,0))/2. - cW*g2UV*YdUVhc(2,0) - (c2W*g2UV2*YdUVhc(2,0))/2. ).imag();
    CdH_32r_LNP = ( cd*YdUVhc(2,1) - (cHW*g2UV*YdUVhc(2,1))/2. - cW*g2UV*YdUVhc(2,1) - (c2W*g2UV2*YdUVhc(2,1))/2. ).real();
    CdH_32i_LNP = ( cd*YdUVhc(2,1) - (cHW*g2UV*YdUVhc(2,1))/2. - cW*g2UV*YdUVhc(2,1) - (c2W*g2UV2*YdUVhc(2,1))/2. ).imag();
    CdH_33r_LNP = ( cd*YdUVhc(2,2) - (cHW*g2UV*YdUVhc(2,2))/2. - cW*g2UV*YdUVhc(2,2) - (c2W*g2UV2*YdUVhc(2,2))/2. ).real();
    CdH_33i_LNP = ( cd*YdUVhc(2,2) - (cHW*g2UV*YdUVhc(2,2))/2. - cW*g2UV*YdUVhc(2,2) - (c2W*g2UV2*YdUVhc(2,2))/2. ).imag();
    CeH_11r_LNP = ( ce*YeUVhc(0,0) - (cHW*g2UV*YeUVhc(0,0))/2. - cW*g2UV*YeUVhc(0,0) - (c2W*g2UV2*YeUVhc(0,0))/2. ).real();
    CeH_11i_LNP = ( ce*YeUVhc(0,0) - (cHW*g2UV*YeUVhc(0,0))/2. - cW*g2UV*YeUVhc(0,0) - (c2W*g2UV2*YeUVhc(0,0))/2. ).imag();
    CeH_12r_LNP = ( ce*YeUVhc(0,1) - (cHW*g2UV*YeUVhc(0,1))/2. - cW*g2UV*YeUVhc(0,1) - (c2W*g2UV2*YeUVhc(0,1))/2. ).real();
    CeH_12i_LNP = ( ce*YeUVhc(0,1) - (cHW*g2UV*YeUVhc(0,1))/2. - cW*g2UV*YeUVhc(0,1) - (c2W*g2UV2*YeUVhc(0,1))/2. ).imag();
    CeH_13r_LNP = ( ce*YeUVhc(0,2) - (cHW*g2UV*YeUVhc(0,2))/2. - cW*g2UV*YeUVhc(0,2) - (c2W*g2UV2*YeUVhc(0,2))/2. ).real();
    CeH_13i_LNP = ( ce*YeUVhc(0,2) - (cHW*g2UV*YeUVhc(0,2))/2. - cW*g2UV*YeUVhc(0,2) - (c2W*g2UV2*YeUVhc(0,2))/2. ).imag();
    CeH_21r_LNP = ( ce*YeUVhc(1,0) - (cHW*g2UV*YeUVhc(1,0))/2. - cW*g2UV*YeUVhc(1,0) - (c2W*g2UV2*YeUVhc(1,0))/2. ).real();
    CeH_21i_LNP = ( ce*YeUVhc(1,0) - (cHW*g2UV*YeUVhc(1,0))/2. - cW*g2UV*YeUVhc(1,0) - (c2W*g2UV2*YeUVhc(1,0))/2. ).imag();
    CeH_22r_LNP = ( ce*YeUVhc(1,1) - (cHW*g2UV*YeUVhc(1,1))/2. - cW*g2UV*YeUVhc(1,1) - (c2W*g2UV2*YeUVhc(1,1))/2. ).real();
    CeH_22i_LNP = ( ce*YeUVhc(1,1) - (cHW*g2UV*YeUVhc(1,1))/2. - cW*g2UV*YeUVhc(1,1) - (c2W*g2UV2*YeUVhc(1,1))/2. ).imag();
    CeH_23r_LNP = ( ce*YeUVhc(1,2) - (cHW*g2UV*YeUVhc(1,2))/2. - cW*g2UV*YeUVhc(1,2) - (c2W*g2UV2*YeUVhc(1,2))/2. ).real();
    CeH_23i_LNP = ( ce*YeUVhc(1,2) - (cHW*g2UV*YeUVhc(1,2))/2. - cW*g2UV*YeUVhc(1,2) - (c2W*g2UV2*YeUVhc(1,2))/2. ).imag();
    CeH_31r_LNP = ( ce*YeUVhc(2,0) - (cHW*g2UV*YeUVhc(2,0))/2. - cW*g2UV*YeUVhc(2,0) - (c2W*g2UV2*YeUVhc(2,0))/2. ).real();
    CeH_31i_LNP = ( ce*YeUVhc(2,0) - (cHW*g2UV*YeUVhc(2,0))/2. - cW*g2UV*YeUVhc(2,0) - (c2W*g2UV2*YeUVhc(2,0))/2. ).imag();
    CeH_32r_LNP = ( ce*YeUVhc(2,1) - (cHW*g2UV*YeUVhc(2,1))/2. - cW*g2UV*YeUVhc(2,1) - (c2W*g2UV2*YeUVhc(2,1))/2. ).real();
    CeH_32i_LNP = ( ce*YeUVhc(2,1) - (cHW*g2UV*YeUVhc(2,1))/2. - cW*g2UV*YeUVhc(2,1) - (c2W*g2UV2*YeUVhc(2,1))/2. ).imag();
    CeH_33r_LNP = ( ce*YeUVhc(2,2) - (cHW*g2UV*YeUVhc(2,2))/2. - cW*g2UV*YeUVhc(2,2) - (c2W*g2UV2*YeUVhc(2,2))/2. ).real();
    CeH_33i_LNP = ( ce*YeUVhc(2,2) - (cHW*g2UV*YeUVhc(2,2))/2. - cW*g2UV*YeUVhc(2,2) - (c2W*g2UV2*YeUVhc(2,2))/2. ).imag();
    CuG_33r_LNP = ctG;
    CuG_33i_LNP = 0.0;
    CuW_33r_LNP = ctW;
    CuW_33i_LNP = 0.0;
    CuB_33r_LNP = ctB;
    CuB_33i_LNP = 0.0;
    CHl1_11r_LNP = (cB*g1UV)/2. + (cHB*g1UV)/4. + (c2B*g1UV2)/2.;
    CHl1_22r_LNP = (cB*g1UV)/2. + (cHB*g1UV)/4. + (c2B*g1UV2)/2.;
    CHl1_33r_LNP = (cB*g1UV)/2. + (cHB*g1UV)/4. + (c2B*g1UV2)/2.;
    CHl3_11r_LNP = -0.25*(cHW*g2UV) - (cW*g2UV)/2. - (c2W*g2UV2)/2.;
    CHl3_22r_LNP = -0.25*(cHW*g2UV) - (cW*g2UV)/2. - (c2W*g2UV2)/2.;
    CHl3_33r_LNP = -0.25*(cHW*g2UV) - (cW*g2UV)/2. - (c2W*g2UV2)/2.;
    CHe_11r_LNP = cB*g1UV + (cHB*g1UV)/2. + c2B*g1UV2;
    CHe_22r_LNP = cB*g1UV + (cHB*g1UV)/2. + c2B*g1UV2;
    CHe_33r_LNP = cB*g1UV + (cHB*g1UV)/2. + c2B*g1UV2;
    CHq1_11r_LNP = -0.16666666666666666*(cB*g1UV) - (cHB*g1UV)/12. - (c2B*g1UV2)/6.;
    CHq1_22r_LNP = -0.16666666666666666*(cB*g1UV) - (cHB*g1UV)/12. - (c2B*g1UV2)/6.;
    CHq1_33r_LNP = -0.16666666666666666*(cB*g1UV) - (cHB*g1UV)/12. - (cqD1*g1UV)/2. - (c2B*g1UV2)/6.;
    CHq3_11r_LNP = -0.25*(cHW*g2UV) - (cW*g2UV)/2. - (c2W*g2UV2)/2.;
    CHq3_22r_LNP = -0.25*(cHW*g2UV) - (cW*g2UV)/2. - (c2W*g2UV2)/2.;
    CHq3_33r_LNP = -0.25*(cHW*g2UV) - (cqD3*g2UV)/2. - (cW*g2UV)/2. - (c2W*g2UV2)/2.;
    CHu_11r_LNP = (-2*cB*g1UV)/3. - (cHB*g1UV)/3. - (2*c2B*g1UV2)/3.;
    CHu_22r_LNP = (-2*cB*g1UV)/3. - (cHB*g1UV)/3. - (2*c2B*g1UV2)/3.;
    CHu_33r_LNP = (-2*cB*g1UV)/3. - (cHB*g1UV)/3. - (ctD*g1UV)/2. - (2*c2B*g1UV2)/3.;
    CHd_11r_LNP = (cB*g1UV)/3. + (cHB*g1UV)/6. + (c2B*g1UV2)/3.;
    CHd_22r_LNP = (cB*g1UV)/3. + (cHB*g1UV)/6. + (c2B*g1UV2)/3.;
    CHd_33r_LNP = (cB*g1UV)/3. + (cHB*g1UV)/6. + (c2B*g1UV2)/3.;
    Cll_1111r_LNP = -0.25*(c2B*g1UV2) - (c2W*g2UV2)/4.;
    Cll_1122r_LNP = -0.25*(c2B*g1UV2) + (c2W*g2UV2)/4.;
    Cll_1133r_LNP = -0.25*(c2B*g1UV2) + (c2W*g2UV2)/4.;
    Cll_1221r_LNP = -0.5*(c2W*g2UV2);
    Cll_1331r_LNP = -0.5*(c2W*g2UV2);
    Cll_2222r_LNP = -0.25*(c2B*g1UV2) - (c2W*g2UV2)/4.;
    Cll_2233r_LNP = -0.25*(c2B*g1UV2) + (c2W*g2UV2)/4.;
    Cll_2332r_LNP = -0.5*(c2W*g2UV2);
    Cll_3333r_LNP = -0.25*(c2B*g1UV2) - (c2W*g2UV2)/4.;
    Cqq1_1111r_LNP = -0.027777777777777776*(c2B*g1UV2) - (c2G*g3UV2)/12.;
    Cqq1_1122r_LNP = -0.027777777777777776*(c2B*g1UV2) + (c2G*g3UV2)/6.;
    Cqq1_1133r_LNP = -0.16666666666666666*(cqD1*g1UV) - (c2B*g1UV2)/36. + (c2G*g3UV2)/6.;
    Cqq1_1221r_LNP = -0.25*(c2G*g3UV2);
    Cqq1_1331r_LNP = -0.25*(c2G*g3UV2);
    Cqq1_2222r_LNP = -0.027777777777777776*(c2B*g1UV2) - (c2G*g3UV2)/12.;
    Cqq1_2233r_LNP = -0.16666666666666666*(cqD1*g1UV) - (c2B*g1UV2)/36. + (c2G*g3UV2)/6.;
    Cqq1_2332r_LNP = -0.25*(c2G*g3UV2);
    Cqq1_3333r_LNP = cqq1 - (cqD1*g1UV)/6. - (c2B*g1UV2)/36. - (c2G*g3UV2)/12.;
    Cqq3_1111r_LNP = -(c2W*g2UV2) - (c2G*g3UV2)/4.;
    Cqq3_1122r_LNP = -(c2W*g2UV2);
    Cqq3_1133r_LNP = -0.5*(cqD3*g2UV) - c2W*g2UV2;
    Cqq3_1221r_LNP = -0.25*(c2G*g3UV2);
    Cqq3_1331r_LNP = -0.25*(c2G*g3UV2);
    Cqq3_2222r_LNP = -(c2W*g2UV2) - (c2G*g3UV2)/4.;
    Cqq3_2233r_LNP = -0.5*(cqD3*g2UV) - c2W*g2UV2;
    Cqq3_2332r_LNP = -0.25*(c2G*g3UV2);
    Cqq3_3333r_LNP = cqq3 - (cqD3*g2UV)/2. - c2W*g2UV2 - (c2G*g3UV2)/4.;
    Clq1_1111r_LNP = (c2B*g1UV2)/6.;
    Clq1_1122r_LNP = (c2B*g1UV2)/6.;
    Clq1_1133r_LNP = (c2B*g1UV2)/6.;
    Clq1_2211r_LNP = (c2B*g1UV2)/6.;
    Clq1_2222r_LNP = (c2B*g1UV2)/6.;
    Clq1_2233r_LNP = (c2B*g1UV2)/6.;
    Clq1_3311r_LNP = (c2B*g1UV2)/6.;
    Clq1_3322r_LNP = (c2B*g1UV2)/6.;
    Clq1_3333r_LNP = (c2B*g1UV2)/6.;
    Clq3_1111r_LNP = -0.5*(c2W*g2UV2);
    Clq3_1122r_LNP = -0.5*(c2W*g2UV2);
    Clq3_1133r_LNP = -0.5*(c2W*g2UV2);
    Clq3_2211r_LNP = -0.5*(c2W*g2UV2);
    Clq3_2222r_LNP = -0.5*(c2W*g2UV2);
    Clq3_2233r_LNP = -0.5*(c2W*g2UV2);
    Clq3_3311r_LNP = -0.5*(c2W*g2UV2);
    Clq3_3322r_LNP = -0.5*(c2W*g2UV2);
    Clq3_3333r_LNP = -0.5*(c2W*g2UV2);
    Cee_1111r_LNP = -(c2B*g1UV2);
    Cee_1122r_LNP = -(c2B*g1UV2);
    Cee_1133r_LNP = -(c2B*g1UV2);
    Cee_2222r_LNP = -(c2B*g1UV2);
    Cee_2233r_LNP = -(c2B*g1UV2);
    Cee_3333r_LNP = -(c2B*g1UV2);
    Cuu_1111r_LNP = (-4*c2B*g1UV2)/9. - (c2G*g3UV2)/3.;
    Cuu_1122r_LNP = (-4*c2B*g1UV2)/9. + (c2G*g3UV2)/6.;
    Cuu_1133r_LNP = (-2*ctD*g1UV)/3. - (4*c2B*g1UV2)/9. + (c2G*g3UV2)/6.;
    Cuu_1221r_LNP = -0.5*(c2G*g3UV2);
    Cuu_1331r_LNP = -0.5*(c2G*g3UV2);
    Cuu_2222r_LNP = (-4*c2B*g1UV2)/9. - (c2G*g3UV2)/3.;
    Cuu_2233r_LNP = (-2*ctD*g1UV)/3. - (4*c2B*g1UV2)/9. + (c2G*g3UV2)/6.;
    Cuu_2332r_LNP = -0.5*(c2G*g3UV2);
    Cuu_3333r_LNP = ctt - (2*ctD*g1UV)/3. - (4*c2B*g1UV2)/9. - (c2G*g3UV2)/3.;
    Cdd_1111r_LNP = -0.1111111111111111*(c2B*g1UV2) - (c2G*g3UV2)/3.;
    Cdd_1122r_LNP = -0.1111111111111111*(c2B*g1UV2) + (c2G*g3UV2)/6.;
    Cdd_1133r_LNP = -0.1111111111111111*(c2B*g1UV2) + (c2G*g3UV2)/6.;
    Cdd_1221r_LNP = -0.5*(c2G*g3UV2);
    Cdd_1331r_LNP = -0.5*(c2G*g3UV2);
    Cdd_2222r_LNP = -0.1111111111111111*(c2B*g1UV2) - (c2G*g3UV2)/3.;
    Cdd_2233r_LNP = -0.1111111111111111*(c2B*g1UV2) + (c2G*g3UV2)/6.;
    Cdd_2332r_LNP = -0.5*(c2G*g3UV2);
    Cdd_3333r_LNP = -0.1111111111111111*(c2B*g1UV2) - (c2G*g3UV2)/3.;
    Ceu_1111r_LNP = (4*c2B*g1UV2)/3.;
    Ceu_1122r_LNP = (4*c2B*g1UV2)/3.;
    Ceu_1133r_LNP = ctD*g1UV + (4*c2B*g1UV2)/3.;
    Ceu_2211r_LNP = (4*c2B*g1UV2)/3.;
    Ceu_2222r_LNP = (4*c2B*g1UV2)/3.;
    Ceu_2233r_LNP = ctD*g1UV + (4*c2B*g1UV2)/3.;
    Ceu_3311r_LNP = (4*c2B*g1UV2)/3.;
    Ceu_3322r_LNP = (4*c2B*g1UV2)/3.;
    Ceu_3333r_LNP = ctD*g1UV + (4*c2B*g1UV2)/3.;
    Ced_1111r_LNP = (-2*c2B*g1UV2)/3.;
    Ced_1122r_LNP = (-2*c2B*g1UV2)/3.;
    Ced_1133r_LNP = (-2*c2B*g1UV2)/3.;
    Ced_2211r_LNP = (-2*c2B*g1UV2)/3.;
    Ced_2222r_LNP = (-2*c2B*g1UV2)/3.;
    Ced_2233r_LNP = (-2*c2B*g1UV2)/3.;
    Ced_3311r_LNP = (-2*c2B*g1UV2)/3.;
    Ced_3322r_LNP = (-2*c2B*g1UV2)/3.;
    Ced_3333r_LNP = (-2*c2B*g1UV2)/3.;
    Cud1_1111r_LNP = (4*c2B*g1UV2)/9.;
    Cud1_1122r_LNP = (4*c2B*g1UV2)/9.;
    Cud1_1133r_LNP = (4*c2B*g1UV2)/9.;
    Cud1_2211r_LNP = (4*c2B*g1UV2)/9.;
    Cud1_2222r_LNP = (4*c2B*g1UV2)/9.;
    Cud1_2233r_LNP = (4*c2B*g1UV2)/9.;
    Cud1_3311r_LNP = (ctD*g1UV)/3. + (4*c2B*g1UV2)/9.;
    Cud1_3322r_LNP = (ctD*g1UV)/3. + (4*c2B*g1UV2)/9.;
    Cud1_3333r_LNP = (ctD*g1UV)/3. + (4*c2B*g1UV2)/9.;
    Cud8_1111r_LNP = -2*c2G*g3UV2;
    Cud8_1122r_LNP = -2*c2G*g3UV2;
    Cud8_1133r_LNP = -2*c2G*g3UV2;
    Cud8_2211r_LNP = -2*c2G*g3UV2;
    Cud8_2222r_LNP = -2*c2G*g3UV2;
    Cud8_2233r_LNP = -2*c2G*g3UV2;
    Cud8_3311r_LNP = -2*c2G*g3UV2;
    Cud8_3322r_LNP = -2*c2G*g3UV2;
    Cud8_3333r_LNP = -2*c2G*g3UV2;
    Cle_1111r_LNP = -(c2B*g1UV2);
    Cle_1122r_LNP = -(c2B*g1UV2);
    Cle_1133r_LNP = -(c2B*g1UV2);
    Cle_2211r_LNP = -(c2B*g1UV2);
    Cle_2222r_LNP = -(c2B*g1UV2);
    Cle_2233r_LNP = -(c2B*g1UV2);
    Cle_3311r_LNP = -(c2B*g1UV2);
    Cle_3322r_LNP = -(c2B*g1UV2);
    Cle_3333r_LNP = -(c2B*g1UV2);
    Clu_1111r_LNP = (2*c2B*g1UV2)/3.;
    Clu_1122r_LNP = (2*c2B*g1UV2)/3.;
    Clu_1133r_LNP = (ctD*g1UV)/2. + (2*c2B*g1UV2)/3.;
    Clu_2211r_LNP = (2*c2B*g1UV2)/3.;
    Clu_2222r_LNP = (2*c2B*g1UV2)/3.;
    Clu_2233r_LNP = (ctD*g1UV)/2. + (2*c2B*g1UV2)/3.;
    Clu_3311r_LNP = (2*c2B*g1UV2)/3.;
    Clu_3322r_LNP = (2*c2B*g1UV2)/3.;
    Clu_3333r_LNP = (ctD*g1UV)/2. + (2*c2B*g1UV2)/3.;
    Cld_1111r_LNP = -0.3333333333333333*(c2B*g1UV2);
    Cld_1122r_LNP = -0.3333333333333333*(c2B*g1UV2);
    Cld_1133r_LNP = -0.3333333333333333*(c2B*g1UV2);
    Cld_2211r_LNP = -0.3333333333333333*(c2B*g1UV2);
    Cld_2222r_LNP = -0.3333333333333333*(c2B*g1UV2);
    Cld_2233r_LNP = -0.3333333333333333*(c2B*g1UV2);
    Cld_3311r_LNP = -0.3333333333333333*(c2B*g1UV2);
    Cld_3322r_LNP = -0.3333333333333333*(c2B*g1UV2);
    Cld_3333r_LNP = -0.3333333333333333*(c2B*g1UV2);
    Cqe_1111r_LNP = (c2B*g1UV2)/3.;
    Cqe_1122r_LNP = (c2B*g1UV2)/3.;
    Cqe_1133r_LNP = (c2B*g1UV2)/3.;
    Cqe_2211r_LNP = (c2B*g1UV2)/3.;
    Cqe_2222r_LNP = (c2B*g1UV2)/3.;
    Cqe_2233r_LNP = (c2B*g1UV2)/3.;
    Cqe_3311r_LNP = cqD1*g1UV + (c2B*g1UV2)/3.;
    Cqe_3322r_LNP = cqD1*g1UV + (c2B*g1UV2)/3.;
    Cqe_3333r_LNP = cqD1*g1UV + (c2B*g1UV2)/3.;
    Cqu1_1111r_LNP = (-2*c2B*g1UV2)/9.;
    Cqu1_1122r_LNP = (-2*c2B*g1UV2)/9.;
    Cqu1_1133r_LNP = -0.16666666666666666*(ctD*g1UV) - (2*c2B*g1UV2)/9.;
    Cqu1_2211r_LNP = (-2*c2B*g1UV2)/9.;
    Cqu1_2222r_LNP = (-2*c2B*g1UV2)/9.;
    Cqu1_2233r_LNP = -0.16666666666666666*(ctD*g1UV) - (2*c2B*g1UV2)/9.;
    Cqu1_3311r_LNP = (-2*cqD1*g1UV)/3. - (2*c2B*g1UV2)/9.;
    Cqu1_3322r_LNP = (-2*cqD1*g1UV)/3. - (2*c2B*g1UV2)/9.;
    Cqu1_3333r_LNP = cqt1 - (2*cqD1*g1UV)/3. - (ctD*g1UV)/6. - (2*c2B*g1UV2)/9.;
    Cqu8_1111r_LNP = -2*c2G*g3UV2;
    Cqu8_1122r_LNP = -2*c2G*g3UV2;
    Cqu8_1133r_LNP = -2*c2G*g3UV2;
    Cqu8_2211r_LNP = -2*c2G*g3UV2;
    Cqu8_2222r_LNP = -2*c2G*g3UV2;
    Cqu8_2233r_LNP = -2*c2G*g3UV2;
    Cqu8_3311r_LNP = -2*c2G*g3UV2;
    Cqu8_3322r_LNP = -2*c2G*g3UV2;
    Cqu8_3333r_LNP = cqt8 - 2*c2G*g3UV2;
    Cqd1_1111r_LNP = (c2B*g1UV2)/9.;
    Cqd1_1122r_LNP = (c2B*g1UV2)/9.;
    Cqd1_1133r_LNP = (c2B*g1UV2)/9.;
    Cqd1_2211r_LNP = (c2B*g1UV2)/9.;
    Cqd1_2222r_LNP = (c2B*g1UV2)/9.;
    Cqd1_2233r_LNP = (c2B*g1UV2)/9.;
    Cqd1_3311r_LNP = (cqD1*g1UV)/3. + (c2B*g1UV2)/9.;
    Cqd1_3322r_LNP = (cqD1*g1UV)/3. + (c2B*g1UV2)/9.;
    Cqd1_3333r_LNP = (cqD1*g1UV)/3. + (c2B*g1UV2)/9.;
    Cqd8_1111r_LNP = -2*c2G*g3UV2;
    Cqd8_1122r_LNP = -2*c2G*g3UV2;
    Cqd8_1133r_LNP = -2*c2G*g3UV2;
    Cqd8_2211r_LNP = -2*c2G*g3UV2;
    Cqd8_2222r_LNP = -2*c2G*g3UV2;
    Cqd8_2233r_LNP = -2*c2G*g3UV2;
    Cqd8_3311r_LNP = -2*c2G*g3UV2;
    Cqd8_3322r_LNP = -2*c2G*g3UV2;
    Cqd8_3333r_LNP = -2*c2G*g3UV2;
}


bool NPd6SILH::PostUpdate()
{

// Obtain the values of the SM parameters at the UV scale, to do the matching
    ChangeToEvolutorsBasisPureSM();
    double Mu_LEW[3] = {mu_LEW, mc_LEW, mt_LEW};
    double Md_LEW[3] = {md_LEW, ms_LEW, mb_LEW};
    double Me_LEW[3] = {me_LEW, mmu_LEW, mtau_LEW};

    if (FlagRGEci) {
        // SM initial conditions at the UV scale. Use RGEsolver SMEFTEvolEW
        SMEFTEvolEW.GenerateSMInitialConditions(muw, Lambda_NP, SMEFTBasisFlag, "Numeric",
            g1_LEW, g2_LEW, g3_LEW, lambdaH_LEW, mH2_LEW,
            Mu_LEW, Md_LEW, Me_LEW, s12CKM_LEW, s13CKM_LEW, s23CKM_LEW, dCKM_LEW);

    } else {
        // SM initial conditions at the UV scale. Use RGEsolver SMEFTEvolEW
        // Skip RGE by setting the two scales at Lambda_NP for the EFT and to muw for the SM pars
        SMEFTEvolEW.GenerateSMInitialConditions(muw, muw, SMEFTBasisFlag, "Numeric",
            g1_LEW, g2_LEW, g3_LEW, lambdaH_LEW, mH2_LEW,
            Mu_LEW, Md_LEW, Me_LEW, s12CKM_LEW, s13CKM_LEW, s23CKM_LEW, dCKM_LEW);
    }

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


    setNPSMEFTd6GeneralParameters();

    if (!NPSMEFTd6General::PostUpdate()) return (false);

    return (true);

}
