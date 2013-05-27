/* 
 * Copyright (C) 2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include <stdexcept>
#include "NPEffective.h"


const std::string NPEffective::NPEffectiveVars[NNPEffectiveVars] 
= {"cWB_NP", "cH_NP", "cLL_NP", "cHLp_NP", "cHQp_NP", 
   "cHL_NP", "cHQ_NP", "cHE_NP", "cHU_NP", "cHD_NP", "Lambda_NP"};


NPEffective::NPEffective() 
: StandardModel() 
{
}


bool NPEffective::Update(const std::map<std::string,double>& DPars) 
{
    for (std::map<std::string, double>::const_iterator it = DPars.begin(); it != DPars.end(); it++)
        SetParameter(it->first, it->second);
    if(!StandardModel::Update(DPars)) return (false);

    return (true);
}


bool NPEffective::Init(const std::map<std::string, double>& DPars) 
{
    Update(DPars);
    return(CheckParameters(DPars)); 
}


bool NPEffective::CheckParameters(const std::map<std::string, double>& DPars) 
{
    for (int i = 0; i < NNPEffectiveVars; i++) {
        if (DPars.find(NPEffectiveVars[i]) == DPars.end()) {
            std::cout << "ERROR: Missing mandatory NPEffective parameter" 
                      << NPEffectiveVars[i] << std::endl;
            return false;
        }
    }
    return(StandardModel::CheckParameters(DPars));
}

    
void NPEffective::SetParameter(const std::string name, const double& value) 
{
    if (name.compare("cWB_NP") == 0)
        cWB = value;
    else if (name.compare("cH_NP") == 0)
        cH = value;
    else if (name.compare("cLL_NP") == 0)
        cLL = value;
    else if (name.compare("cHLp_NP") == 0)
        cHLp = value;
    else if (name.compare("cHQp_NP") == 0)
        cHQp = value;
    else if (name.compare("cHL_NP") == 0)
        cHL = value;
    else if (name.compare("cHQ_NP") == 0)
        cHQ = value;
    else if (name.compare("cHE_NP") == 0)
        cHE = value;
    else if (name.compare("cHU_NP") == 0)
        cHU = value;
    else if (name.compare("cHD_NP") == 0)
        cHD = value;
    else if (name.compare("Lambda_NP") == 0)
        LambdaNP = value;
    else
        StandardModel::SetParameter(name, value);       
}


bool NPEffective::InitializeModel() 
{
    SetModelInitialized(StandardModel::InitializeModel());
    return (IsModelInitialized());
}


void NPEffective::SetEWSMflags(EWSM& myEWSM) 
{
    StandardModel::SetEWSMflags(myEWSM);
}


bool NPEffective::SetFlag(const std::string name, const bool& value) 
{
    bool res = false;
    
    if (name.compare("EWABC") == 0) 
        throw std::runtime_error("ERROR: Flag EWABC is not applicable to NPEffective"); 
    else if (name.compare("EWABC2") == 0)
        throw std::runtime_error("ERROR: Flag EWABC2 is not applicable to NPEffective"); 
    else if (name.compare("epsilon1SM") == 0) 
        throw std::runtime_error("ERROR: Flag epsilon1SM is not applicable to NPEffective"); 
    else if (name.compare("epsilon2SM") == 0) 
        throw std::runtime_error("ERROR: Flag epsilon2SM is not applicable to NPEffective"); 
    else if (name.compare("epsilon3SM") == 0) 
        throw std::runtime_error("ERROR: Flag epsilon3SM is not applicable to NPEffective"); 
    else if (name.compare("epsilonbSM") == 0) 
        throw std::runtime_error("ERROR: Flag epsilonbSM is not applicable to NPEffective"); 
    else if (name.compare("NotLinearizedNP") == 0)
        throw std::runtime_error("ERROR: Flag NotLinearizedNP is not applicable to NPEffective");
    else
        res = StandardModel::SetFlag(name,value);

    return(res);
}


////////////////////////////////////////////////////////////////////////

double NPEffective::v() const
{
    double ratio = StandardModel::v()*StandardModel::v()/LambdaNP/LambdaNP;
    double DelGF = - (cLL - 2.0*cHLp)*ratio;

    return ( sqrt( (1.0 + DelGF)/sqrt(2.0)/GF ) );
}


double NPEffective::Mw_tree() const
{
    double GF0 = GF*(1.0 - DeltaGF());
    double tmp = 4.0*M_PI*ale/sqrt(2.0)/GF0/Mz/Mz;
    return ( Mz/sqrt(2.0) * sqrt(1.0 + sqrt(1.0 - tmp)) );
}


double NPEffective::DeltaGF() const
{
    double ratio = v()*v()/LambdaNP/LambdaNP;

    return ( - (cLL - 2.0*cHLp)*ratio );
}


////////////////////////////////////////////////////////////////////////     

double NPEffective::obliqueS() const 
{
    double ratio = v()*v()/LambdaNP/LambdaNP;
    double sW_SM = sqrt(StandardModel::sW2()); /* This has to be the SM value. */
    double cW_SM = sqrt(StandardModel::cW2()); /* This has to be the SM value. */

    return ( 4.0*sW_SM*cW_SM*cWB/alphaMz()*ratio );
}


double NPEffective::obliqueT() const
{
    double ratio = v()*v()/LambdaNP/LambdaNP;

    return ( - cH/2.0/alphaMz()*ratio );    
}


double NPEffective::obliqueU() const
{
    return 0.0;
}


////////////////////////////////////////////////////////////////////////     

double NPEffective::deltaGLl(StandardModel::lepton l) const
{
    double ratio = v()*v()/LambdaNP/LambdaNP;
    double sW2_SM = StandardModel::sW2(); /* This has to be the SM value. */
    double cW2_SM = StandardModel::cW2(); /* This has to be the SM value. */
    double gVf_SM = StandardModel::gVl(l).real(); /* This has to be the SM value. */
    double gAf_SM = StandardModel::gAl(l).real(); /* This has to be the SM value. */
    switch (l) {
        case StandardModel::NEUTRINO_1:
        case StandardModel::NEUTRINO_2:
        case StandardModel::NEUTRINO_3:
            return ( (cHLp - cHL)/2.0*ratio - (gVf_SM + gAf_SM)/4.0*DeltaGF() );
        case StandardModel::ELECTRON:
        case StandardModel::MU:
        case StandardModel::TAU:
            return ( - (cHLp + cHL)/2.0*ratio - (gVf_SM + gAf_SM)/4.0*DeltaGF()
                     + (gVf_SM - gAf_SM)*cW2_SM/2.0/(cW2_SM - sW2_SM)*DeltaGF() );
        default:
            throw std::runtime_error("Error in NPEffective::deltaGLl()");        
    }   
}


double NPEffective::deltaGLq(StandardModel::quark q) const
{
    double ratio = v()*v()/LambdaNP/LambdaNP;
    double sW2_SM = StandardModel::sW2(); /* This has to be the SM value. */
    double cW2_SM = StandardModel::cW2(); /* This has to be the SM value. */
    double gVf_SM = StandardModel::gVq(q).real(); /* This has to be the SM value. */
    double gAf_SM = StandardModel::gAq(q).real(); /* This has to be the SM value. */
    switch (q) {
        case StandardModel::UP:
        case StandardModel::CHARM:
           return ( (cHQp - cHQ)/2.0*ratio - (gVf_SM + gAf_SM)/4.0*DeltaGF()
                     + (gVf_SM - gAf_SM)*cW2_SM/2.0/(cW2_SM - sW2_SM)*DeltaGF() );
        case StandardModel::TOP:
            return 0.0;
        case StandardModel::DOWN:
        case StandardModel::STRANGE:
        case StandardModel::BOTTOM:
            return ( - (cHQp + cHQ)/2.0*ratio - (gVf_SM + gAf_SM)/4.0*DeltaGF()
                     + (gVf_SM - gAf_SM)*cW2_SM/2.0/(cW2_SM - sW2_SM)*DeltaGF() );
        default:
            throw std::runtime_error("Error in NPEffective::deltaGLq()");        
    }
}


double NPEffective::deltaGRl(StandardModel::lepton l) const
{
    double ratio = v()*v()/LambdaNP/LambdaNP;
    double sW2_SM = StandardModel::sW2(); /* This has to be the SM value. */
    double cW2_SM = StandardModel::cW2(); /* This has to be the SM value. */
    double gVf_SM = StandardModel::gVl(l).real(); /* This has to be the SM value. */
    double gAf_SM = StandardModel::gAl(l).real(); /* This has to be the SM value. */
    switch (l) {
        case StandardModel::NEUTRINO_1:
        case StandardModel::NEUTRINO_2:
        case StandardModel::NEUTRINO_3:
            return 0.0;
        case StandardModel::ELECTRON:
        case StandardModel::MU:
        case StandardModel::TAU:
            return ( - cHE/2.0*ratio 
                     + (gVf_SM - gAf_SM)/4.0/(cW2_SM - sW2_SM)*DeltaGF() );
        default:
            throw std::runtime_error("Error in NPEffective::deltaGRl()");        
    } 
}


double NPEffective::deltaGRq(StandardModel::quark q) const
{
    double ratio = v()*v()/LambdaNP/LambdaNP;
    double sW2_SM = StandardModel::sW2(); /* This has to be the SM value. */
    double cW2_SM = StandardModel::cW2(); /* This has to be the SM value. */
    double gVf_SM = StandardModel::gVq(q).real(); /* This has to be the SM value. */
    double gAf_SM = StandardModel::gAq(q).real(); /* This has to be the SM value. */
    switch (q) {
        case StandardModel::UP:
        case StandardModel::CHARM:
            return ( - cHU/2.0*ratio 
                     + (gVf_SM - gAf_SM)/4.0/(cW2_SM - sW2_SM)*DeltaGF() );
        case StandardModel::TOP:
            return 0.0;
        case StandardModel::DOWN:
        case StandardModel::STRANGE:
        case StandardModel::BOTTOM:
            return ( - cHD/2.0*ratio 
                     + (gVf_SM - gAf_SM)/4.0/(cW2_SM - sW2_SM)*DeltaGF() );
        default:
            throw std::runtime_error("Error in NPEffective::deltaGRq()");        
    }   
}


double NPEffective::deltaGVl(StandardModel::lepton l) const
{
    return ( deltaGLl(l) + deltaGRl(l) + StandardModel::deltaGVl(l) );
}


double NPEffective::deltaGVq(StandardModel::quark q) const
{
    return ( deltaGLq(q) + deltaGRq(q) + StandardModel::deltaGVq(q) );
}


double NPEffective::deltaGAl(StandardModel::lepton l) const
{
    return ( deltaGLl(l) - deltaGRl(l) + StandardModel::deltaGAl(l) );
}


double NPEffective::deltaGAq(StandardModel::quark q) const
{
    return ( deltaGLq(q) - deltaGRq(q) + StandardModel::deltaGAq(q) );
}


////////////////////////////////////////////////////////////////////////

double NPEffective::epsilon1() const
{
    throw std::runtime_error("ERROR: NPEffective::epsilon1() is not implemented");
}


double NPEffective::epsilon2() const
{
    throw std::runtime_error("ERROR: NPEffective::epsilon2() is not implemented");
}


double NPEffective::epsilon3() const
{
    throw std::runtime_error("ERROR: NPEffective::epsilon3() is not implemented");
}


double NPEffective::epsilonb() const
{
    throw std::runtime_error("ERROR: NPEffective::epsilonb() is not implemented");
}


////////////////////////////////////////////////////////////////////////

double NPEffective::Mw() const
{
    double myMw = StandardModel::Mw();

    if (!IsFlagNotLinearizedNP() ) {
        double alpha = StandardModel::alphaMz();
        double c2 = StandardModel::cW2();
        double s2 = StandardModel::sW2();

        myMw *= 1.0 - alpha/4.0/(c2-s2)
                *( obliqueS() - 2.0*c2*obliqueT() - (c2-s2)*obliqueU()/2.0/s2 )
                - s2/2.0/(c2-s2)*DeltaGF();
    } else
        if (obliqueS()!=0.0 || obliqueT()!=0.0 || obliqueU()!=0.0)
            throw std::runtime_error("NPEffective::Mw(): The oblique corrections STU cannot be used with flag NotLinearizedNP=1");

    return myMw;
}


double NPEffective::cW2() const
{
    return ( Mw()*Mw()/Mz/Mz );
}


double NPEffective::sW2() const
{
    return ( 1.0 - cW2() );
}


double NPEffective::GammaW() const
{
    double Gamma_W = StandardModel::GammaW();

    if (!IsFlagNotLinearizedNP() ) {
        double alpha = StandardModel::alphaMz();
        double c2 = StandardModel::cW2();
        double s2 = StandardModel::sW2();

        Gamma_W *= 1.0 - 3.0*alpha/4.0/(c2-s2)
                   *( obliqueS() - 2.0*c2*obliqueT()
                      - (c2-s2)*obliqueU()/2.0/s2 )
                    - 3.0*s2/2.0/(c2-s2)*DeltaGF();
        } else
            if (obliqueS()!=0.0 || obliqueT()!=0.0 || obliqueU()!=0.0)
                throw std::runtime_error("NPEffective::GammaW(): The oblique corrections STU cannot be used with flag NotLinearizedNP=1");

    return Gamma_W;
}






