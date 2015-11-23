/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include <math.h>
#include "MFV.h"

const std::string MFV::MFVvars[NMFVvars] = {
    "a1", "a2", "a3", "a4r", "a4i", "a5r", "a5i", "a6", "a7", "a8r", "a8i",
    "x1", "x2",
    "y1", "y2r", "y2i", "y3", "y4r", "y4i", "y5r", "y5i", "y6", "y7",
    "w1r", "w1i", "w2r", "w2i", "w3r", "w3i", "w4r", "w4i", "w5r", "w5i"
};

MFV::MFV()
: SUSY(), X()
{
}

bool MFV::InitializeModel()
{
    setModelInitialized(SUSY::InitializeModel());
    return (IsModelInitialized());
}

bool MFV::Init(const std::map<std::string, double>& DPars)
{
    return(SUSY::Init(DPars));
}

bool MFV::PreUpdate()
{    
    if(!SUSY::PreUpdate()) return (false);
    return (true);
}

bool MFV::Update(const std::map<std::string, double>& DPars)
{
    if(!PreUpdate()) return (false);
    
    UpdateError = false;
    
    for (std::map<std::string, double>::const_iterator it = DPars.begin(); it != DPars.end(); it++)
        setParameter(it->first, it->second);
    
    if (UpdateError) return (false);
    
    if(!PostUpdate()) return (false);
    
    return (true);
}

bool MFV::PostUpdate()
{
    if (!SUSY::PostUpdate()) return (false);
    return (true);
}

void MFV::setParameter(const std::string name, const double& value)
{
    if(name.compare("a1") == 0)
        a1 = value;
    else if(name.compare("a2") == 0)
        a2 = value;
    else if(name.compare("a3") == 0)
        a3 = value;
    else if(name.compare("a4r") == 0)
        a4.real() = value;
    else if(name.compare("a4i") == 0)
        a4.imag() = value;
    else if(name.compare("a5r") == 0)
        a5.real() = value;
    else if(name.compare("a5i") == 0)
        a5.imag() = value;
    else if(name.compare("a6") == 0)
        a6 = value;
    else if(name.compare("a7") == 0)
        a7 = value;
    else if(name.compare("a8r") == 0)
        a8.real() = value;
    else if(name.compare("a8i") == 0)
        a8.imag() = value;
    else if(name.compare("x1") == 0)
        x1 = value;
    else if(name.compare("x2") == 0)
        x2 = value;
    else if(name.compare("y1") == 0)
        y1 = value;
    else if(name.compare("y2r") == 0)
        y2.real() = value;
    else if(name.compare("y2i") == 0)
        y2.imag() = value;
    else if(name.compare("y3") == 0)
        y3 = value;
    else if(name.compare("y4r") == 0)
        y4.real() = value;
    else if(name.compare("y4i") == 0)
        y4.imag() = value;
    else if(name.compare("y5r") == 0)
        y5.real() = value;
    else if(name.compare("y5i") == 0)
        y5.imag() = value;
    else if(name.compare("y6") == 0)
        y6 = value;
    else if(name.compare("y7") == 0)
        y7 = value;
    else if(name.compare("w1r") == 0)
        w1.real() = value;
    else if(name.compare("w1i") == 0)
        w1.imag() = value;
    else if(name.compare("w2r") == 0)
        w2.real() = value;
    else if(name.compare("w2i") == 0)
        w2.imag() = value;
    else if(name.compare("w3r") == 0)
        w3.real() = value;
    else if(name.compare("w3i") == 0)
        w3.imag() = value;
    else if(name.compare("w4r") == 0)
        w4.real() = value;
    else if(name.compare("w4i") == 0)
        w4.imag() = value;
    else if(name.compare("w5r") == 0)
        w5.real() = value;
    else if(name.compare("w5i") == 0)
        w5.imag() = value;
    else
        SUSY::setParameter(name, value);
}

bool MFV::CheckParameters(const std::map<std::string, double>& DPars)
{
    for (int i = 0; i < NMFVvars; i++) {
        if (DPars.find(MFVvars[i]) == DPars.end()) {
            std::cout << "missing mandatory MFV parameter " << MFVvars[i] << std::endl;
            return false;
        }
    }
    return(SUSY::CheckParameters(DPars));
}

void MFV::SetSoftTerms()
{
    // Colangelo's expressions in Colangelo's basis
    X.Update(myCKM);
    msQhat2 = gslpp::matrix<gslpp::complex>::Id(3) * a1 + X.GetX13() * x1 + X.GetX1() * y1 +
               X.GetX5() * y2 + X.GetX9() * y2.conjugate();
    msUhat2 = gslpp::matrix<gslpp::complex>::Id(3) * a2 + X.GetX1() * x2;
    msDhat2 = gslpp::matrix<gslpp::complex>::Id(3) * a3 + X.GetX1() * y3 + X.GetX3() * w1 +
               X.GetX4() * w1.conjugate();
    TUhat = X.GetX5() * a4 + X.GetX1() * y4 + X.GetX6() * w2;
    TDhat = X.GetX1() * a5 + X.GetX5() * y5 + X.GetX2() * w3 + X.GetX4() * w4;
    msLhat2 = gslpp::matrix<gslpp::complex>::Id(3) * a6 + X.GetX1() * y6;
    msEhat2 = gslpp::matrix<gslpp::complex>::Id(3) * a7 + X.GetX1() * y7;
    TEhat = X.GetX1() * a8 + X.GetX2() * w5;
    
    // rotation to the SCKM basis according to SLHA notation
    
    gslpp::matrix<gslpp::complex> ckm(3,3,0.);
    myCKM.getCKM(ckm);
  
    TUhat = sqrt(2.) * TUhat * ckm.hconjugate();
    TDhat = sqrt(2.) * TDhat;
    TEhat = sqrt(2.) * TEhat;

    // msNhat2 and TNhat remain 0.
}

