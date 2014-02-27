/* 
 * Copyright (C) 2012 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include <stdexcept>
#include <EWSM.h>
#include "NPSTUVWXY.h"


const std::string NPSTUVWXY::STUVWXYvars[NSTUVWXYvars] 
= {"obliqueShat", "obliqueThat", "obliqueUhat", 
   "obliqueV", "obliqueW", "obliqueX", "obliqueY"};


NPSTUVWXY::NPSTUVWXY()
: NPbase()
{
}


bool NPSTUVWXY::InitializeModel()
{
    setModelInitialized(NPbase::InitializeModel());
    return (IsModelInitialized());
}


bool NPSTUVWXY::Init(const std::map<std::string, double>& DPars) 
{
    Update(DPars);
    return(CheckParameters(DPars)); 
}


bool NPSTUVWXY::Update(const std::map<std::string,double>& DPars)
{
    for (std::map<std::string, double>::const_iterator it = DPars.begin(); it != DPars.end(); it++)
        setParameter(it->first, it->second);
    if(!NPbase::Update(DPars)) return (false);
    return (true);
}

    
void NPSTUVWXY::setParameter(const std::string name, const double& value)
{
    if (name.compare("obliqueShat") == 0)
        myObliqueShat = value;
    else if (name.compare("obliqueThat") == 0)
        myObliqueThat = value;
    else if (name.compare("obliqueUhat") == 0)
        myObliqueUhat = value;    
    else if (name.compare("obliqueV") == 0)
        myObliqueV = value;    
    else if (name.compare("obliqueW") == 0)
        myObliqueW = value;    
    else if (name.compare("obliqueX") == 0)
        myObliqueX = value;    
    else if (name.compare("obliqueY") == 0)
        myObliqueY = value;    
    else
        NPbase::setParameter(name, value);
}


bool NPSTUVWXY::CheckParameters(const std::map<std::string, double>& DPars)
{
    for (int i = 0; i < NSTUVWXYvars; i++) {
        if (DPars.find(STUVWXYvars[i]) == DPars.end()) {
            std::cout << "ERROR: Missing mandatory NPSTUVWXY parameter "
                      << STUVWXYvars[i] << std::endl;
            return false;
        }
    }
    return(NPbase::CheckParameters(DPars));
}


bool NPSTUVWXY::setFlag(const std::string name, const bool value)
{
    bool res = false;
    res = NPbase::setFlag(name,value);
    return(res);
}


bool NPSTUVWXY::CheckFlags() const
{
    return(NPbase::CheckFlags());
}


////////////////////////////////////////////////////////////////////////

double NPSTUVWXY::epsilon1() const
{
    double c2 = myEWSM->cW2_SM();
    double s2 = myEWSM->sW2_SM();
    double eps1 = myEWSM->epsilon1_SM();
    eps1 += myObliqueThat - myObliqueW + 2.0*sqrt(s2)/sqrt(c2)*myObliqueX
            - s2/c2*myObliqueY;
    return eps1;
}


double NPSTUVWXY::epsilon2() const
{
    double c2 = myEWSM->cW2_SM();
    double s2 = myEWSM->sW2_SM();
    double eps2 = myEWSM->epsilon2_SM();
    eps2 += myObliqueUhat - myObliqueV - myObliqueW + 2.0*sqrt(s2)/sqrt(c2)*myObliqueX;
    return eps2;
}


double NPSTUVWXY::epsilon3() const
{
    double c2 = myEWSM->cW2_SM();
    double s2 = myEWSM->sW2_SM();
    double eps3 = myEWSM->epsilon3_SM();
    eps3 += myObliqueShat  - myObliqueW + myObliqueX/sqrt(s2)/sqrt(c2) - myObliqueY;
    return eps3;
}


double NPSTUVWXY::epsilonb() const
{
    double epsb = myEWSM->epsilonb_SM();
    return epsb;
}


////////////////////////////////////////////////////////////////////////     

double NPSTUVWXY::obliqueS() const
{
    double sW2_SM = myEWSM->sW2_SM();
    double sW_SM = sqrt(sW2_SM);
    double cW_SM = myEWSM->cW2_SM();
    return ( ( myObliqueShat - myObliqueW + myObliqueX/(sW_SM*cW_SM) - myObliqueY )
            * 4.0*sW2_SM/alphaMz() );
}


double NPSTUVWXY::obliqueT() const
{
    double sW2_SM = myEWSM->sW2_SM();
    double sW_SM = sqrt(sW2_SM);
    double cW2_SM = myEWSM->cW2_SM();
    double cW_SM = sqrt(cW2_SM);
    return ( ( myObliqueThat - myObliqueW + 2.0*sW_SM/cW_SM*myObliqueX
            - sW2_SM/cW2_SM*myObliqueY )/alphaMz() );
}


double NPSTUVWXY::obliqueU() const
{
    double sW2_SM = myEWSM->sW2_SM();
    double sW_SM = sqrt(sW2_SM);
    double cW_SM = sqrt(myEWSM->cW2_SM());
    return ( ( - myObliqueUhat + myObliqueV + myObliqueW
            - 2.0*sW_SM/cW_SM*myObliqueX )*4.0*sW2_SM/alphaMz() );
}


////////////////////////////////////////////////////////////////////////

double NPSTUVWXY::GammaW() const
{
    double Gamma_W = myEWSM->GammaW_SM();

    double Wbar = (obliqueV() - obliqueW())/alphaMz();

    double alpha = StandardModel::alphaMz();
    double c2 = myEWSM->cW2_SM();
    double s2 = myEWSM->sW2_SM();

    Gamma_W *= 1.0 - 3.0*alpha/4.0/(c2-s2)
               *( obliqueS() - 2.0*c2*obliqueT()
                  - (c2-s2)*obliqueU()/2.0/s2 - 2.0*(c2 - s2)*Wbar )
               - (1.0 + c2)/2.0/(c2-s2)*DeltaGF();

    return Gamma_W;
}


