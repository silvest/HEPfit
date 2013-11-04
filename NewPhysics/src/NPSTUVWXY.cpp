/* 
 * Copyright (C) 2012-2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include <stdexcept>
#include "NPSTUVWXY.h"


const std::string NPSTUVWXY::STUVWXYvars[NSTUVWXYvars] 
= {"obliqueShat", "obliqueThat", "obliqueUhat", 
   "obliqueV", "obliqueW", "obliqueX", "obliqueY"};


NPSTUVWXY::NPSTUVWXY()
: NPZbbbar() 
{
}


bool NPSTUVWXY::Update(const std::map<std::string,double>& DPars)
{
    for (std::map<std::string, double>::const_iterator it = DPars.begin(); it != DPars.end(); it++)
        setParameters(it->first, it->second);
    if(!NPZbbbar::Update(DPars)) return (false);

    return (true);
}


bool NPSTUVWXY::Init(const std::map<std::string, double>& DPars) 
{
    Update(DPars);
    return(CheckParameters(DPars)); 
}


bool NPSTUVWXY::CheckParameters(const std::map<std::string, double>& DPars) 
{
    for (int i = 0; i < NSTUVWXYvars; i++) {
        if (DPars.find(STUVWXYvars[i]) == DPars.end()) {
            std::cout << "ERROR: Missing mandatory NPSTUVWXY parameter" 
                      << STUVWXYvars[i] << std::endl;
            return false;
        }
    }
    return(NPZbbbar::CheckParameters(DPars));
}

    
void NPSTUVWXY::setParameters(const std::string name, const double& value) 
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
        NPZbbbar::setParameters(name, value);       
}


bool NPSTUVWXY::InitializeModel() 
{
    setModelInitialized(NPZbbbar::InitializeModel());
    return (IsModelInitialized());
}


void NPSTUVWXY::SetEWSMflags(EWSM& myEWSM) 
{
    StandardModel::SetEWSMflags(myEWSM);
}


bool NPSTUVWXY::SetFlag(const std::string name, const bool& value)
{
    bool res = false;
    if (name.compare("EWABC") == 0)
        throw std::runtime_error("ERROR: Flag EWABC is not applicable to NPSTUVWXY"); 
    else if (name.compare("EWABC2") == 0)
        throw std::runtime_error("ERROR: Flag EWABC2 is not applicable to NPSTUVWXY"); 
    else
        res = NPZbbbar::SetFlag(name,value);

    return(res);
}


////////////////////////////////////////////////////////////////////////

double NPSTUVWXY::epsilon1() const
{
    double c2 = StandardModel::cW2();
    double s2 = StandardModel::sW2();
    double eps1 = epsilon1_SM();
    eps1 += myObliqueThat - myObliqueW + 2.0*sqrt(s2)/sqrt(c2)*myObliqueX
            - s2/c2*myObliqueY;
    return eps1;
}


double NPSTUVWXY::epsilon2() const
{
    double c2 = StandardModel::cW2();
    double s2 = StandardModel::sW2();
    double eps2 = epsilon2_SM();
    eps2 += myObliqueUhat - myObliqueV - myObliqueW + 2.0*sqrt(s2)/sqrt(c2)*myObliqueX;
    return eps2;
}


double NPSTUVWXY::epsilon3() const
{
    double c2 = StandardModel::cW2();
    double s2 = StandardModel::sW2();
    double eps3 = epsilon3_SM();
    eps3 += myObliqueShat  - myObliqueW + myObliqueX/sqrt(s2)/sqrt(c2) - myObliqueY;
    return eps3;
}


double NPSTUVWXY::epsilonb() const
{
    double epsb = epsilonb_SM();
    return epsb;
}


////////////////////////////////////////////////////////////////////////     

double NPSTUVWXY::obliqueS() const
{
    double sW2_SM = StandardModel::sW2();
    double sW_SM = sqrt(sW2_SM);
    double cW_SM = sqrt(StandardModel::cW2());
    return ( ( myObliqueShat - myObliqueW + myObliqueX/(sW_SM*cW_SM) - myObliqueY )
            * 4.0*sW2_SM/alphaMz() );
}


double NPSTUVWXY::obliqueT() const
{
    double sW2_SM = StandardModel::sW2();
    double sW_SM = sqrt(sW2_SM);
    double cW2_SM = StandardModel::cW2();
    double cW_SM = sqrt(cW2_SM);
    return ( ( myObliqueThat - myObliqueW + 2.0*sW_SM/cW_SM*myObliqueX
            - sW2_SM/cW2_SM*myObliqueY )/alphaMz() );
}


double NPSTUVWXY::obliqueU() const
{
    double sW2_SM = StandardModel::sW2();
    double sW_SM = sqrt(sW2_SM);
    double cW_SM = sqrt(StandardModel::cW2());
    return ( ( - myObliqueUhat + myObliqueV + myObliqueW
            - 2.0*sW_SM/cW_SM*myObliqueX )*4.0*sW2_SM/alphaMz() );
}


////////////////////////////////////////////////////////////////////////

double NPSTUVWXY::Mw() const
{
    double myMw = StandardModel::Mw();

    if (IsFlagEWBURGESS()) {
        myMw *= 1.0 - 0.00723/2.0*obliqueS() + 0.0111/2.0*obliqueT() + 0.00849/2.0*obliqueU();
        return myMw;
    }

    if (!IsFlagNotLinearizedNP() ) {
        double alpha = StandardModel::alphaMz();
        double c2 = StandardModel::cW2();
        double s2 = StandardModel::sW2();

        myMw *= 1.0 - alpha/4.0/(c2-s2)
                *( obliqueS() - 2.0*c2*obliqueT() - (c2-s2)*obliqueU()/2.0/s2 );
    } else
        if (obliqueS()!=0.0 || obliqueT()!=0.0 || obliqueU()!=0.0)
            throw std::runtime_error("NPSTUVWXY::Mw(): The oblique corrections STU cannot be used with flag NotLinearizedNP=1");

    return myMw;
}


double NPSTUVWXY::cW2() const
{
    return ( Mw()*Mw()/Mz/Mz );
}


double NPSTUVWXY::sW2() const
{
    return ( 1.0 - cW2() );
}


double NPSTUVWXY::GammaW() const
{
    double Gamma_W = StandardModel::GammaW();

    double Wbar = (obliqueV() - obliqueW())/alphaMz();

    if (IsFlagEWBURGESS()) {
        Gamma_W *= 1.0 - 0.00723*obliqueS() + 0.0111*obliqueT()
                   + 0.00849*obliqueU() + 0.00781*Wbar;
        return Gamma_W;
    }

    if (!IsFlagNotLinearizedNP() ) {
        double alpha = StandardModel::alphaMz();
        double c2 = StandardModel::cW2();
        double s2 = StandardModel::sW2();

        Gamma_W *= 1.0 - 3.0*alpha/4.0/(c2-s2)
                   *( obliqueS() - 2.0*c2*obliqueT()
                      - (c2-s2)*obliqueU()/2.0/s2 - 2.0*(c2 - s2)*Wbar );
        } else
            if (obliqueS()!=0.0 || obliqueT()!=0.0 || obliqueU()!=0.0)
                throw std::runtime_error("NPSTUVWXY::GammaW(): The oblique corrections STU cannot be used with flag NotLinearizedNP=1");

    return Gamma_W;
}






