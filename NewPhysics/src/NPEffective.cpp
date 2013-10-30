/* 
 * Copyright (C) 2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include <stdexcept>
#include "NPEffective.h"


NPEffective::NPEffective() 
: StandardModel() 
{
}


bool NPEffective::Update(const std::map<std::string,double>& DPars) 
{
    for (std::map<std::string, double>::const_iterator it = DPars.begin(); it != DPars.end(); it++)
        parseParameters(it->first, it->second);
    if(!StandardModel::Update(DPars)) return (false);

    return (true);
}


bool NPEffective::Init(const std::map<std::string, double>& DPars) 
{
    Update(DPars);
    return(CheckParameters(DPars)); 
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
    //return ( sqrt( (1.0 - (cL1L2 - cHL1p - cHL2p)/sqrt(2.0)/GF/LambdaNP/LambdaNP)
    //               /sqrt(2.0)/GF ) );

    /* use the tree-level relation */
    return StandardModel::v();
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

    return ( - (cL1L2 - cHL1p - cHL2p)*ratio );
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
            return ( (cHL1p - cHL1)/2.0*ratio - (gVf_SM + gAf_SM)/4.0*DeltaGF() );
        case StandardModel::NEUTRINO_2:
            return ( (cHL2p - cHL2)/2.0*ratio - (gVf_SM + gAf_SM)/4.0*DeltaGF() );
        case StandardModel::NEUTRINO_3:
            return ( (cHL3p - cHL3)/2.0*ratio - (gVf_SM + gAf_SM)/4.0*DeltaGF() );
        case StandardModel::ELECTRON:
            return ( - (cHL1p + cHL1)/2.0*ratio - (gVf_SM + gAf_SM)/4.0*DeltaGF()
                     + (gVf_SM - gAf_SM)*cW2_SM/2.0/(cW2_SM - sW2_SM)*DeltaGF() );
        case StandardModel::MU:
            return ( - (cHL2p + cHL2)/2.0*ratio - (gVf_SM + gAf_SM)/4.0*DeltaGF()
                     + (gVf_SM - gAf_SM)*cW2_SM/2.0/(cW2_SM - sW2_SM)*DeltaGF() );
        case StandardModel::TAU:
            return ( - (cHL3p + cHL3)/2.0*ratio - (gVf_SM + gAf_SM)/4.0*DeltaGF()
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
           return ( (cHQ1p - cHQ1)/2.0*ratio - (gVf_SM + gAf_SM)/4.0*DeltaGF()
                     + (gVf_SM - gAf_SM)*cW2_SM/2.0/(cW2_SM - sW2_SM)*DeltaGF() );
        case StandardModel::CHARM:
           return ( (cHQ2p - cHQ2)/2.0*ratio - (gVf_SM + gAf_SM)/4.0*DeltaGF()
                     + (gVf_SM - gAf_SM)*cW2_SM/2.0/(cW2_SM - sW2_SM)*DeltaGF() );
        case StandardModel::TOP:
            return 0.0;
        case StandardModel::DOWN:
            return ( - (cHQ1p + cHQ1)/2.0*ratio - (gVf_SM + gAf_SM)/4.0*DeltaGF()
                     + (gVf_SM - gAf_SM)*cW2_SM/2.0/(cW2_SM - sW2_SM)*DeltaGF() );
        case StandardModel::STRANGE:
            return ( - (cHQ2p + cHQ2)/2.0*ratio - (gVf_SM + gAf_SM)/4.0*DeltaGF()
                     + (gVf_SM - gAf_SM)*cW2_SM/2.0/(cW2_SM - sW2_SM)*DeltaGF() );
        case StandardModel::BOTTOM:
            return ( - (cHQ3p + cHQ3)/2.0*ratio - (gVf_SM + gAf_SM)/4.0*DeltaGF()
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
            return ( - cHE1/2.0*ratio
                     + (gVf_SM - gAf_SM)/4.0/(cW2_SM - sW2_SM)*DeltaGF() );
        case StandardModel::MU:
            return ( - cHE2/2.0*ratio
                     + (gVf_SM - gAf_SM)/4.0/(cW2_SM - sW2_SM)*DeltaGF() );
        case StandardModel::TAU:
            return ( - cHE3/2.0*ratio
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
            return ( - cHU1/2.0*ratio
                     + (gVf_SM - gAf_SM)/4.0/(cW2_SM - sW2_SM)*DeltaGF() );
        case StandardModel::CHARM:
            return ( - cHU2/2.0*ratio
                     + (gVf_SM - gAf_SM)/4.0/(cW2_SM - sW2_SM)*DeltaGF() );
        case StandardModel::TOP:
            return 0.0;
        case StandardModel::DOWN:
            return ( - cHD1/2.0*ratio
                     + (gVf_SM - gAf_SM)/4.0/(cW2_SM - sW2_SM)*DeltaGF() );
        case StandardModel::STRANGE:
            return ( - cHD2/2.0*ratio
                     + (gVf_SM - gAf_SM)/4.0/(cW2_SM - sW2_SM)*DeltaGF() );
        case StandardModel::BOTTOM:
            return ( - cHD3/2.0*ratio
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
        double ratio = v()*v()/LambdaNP/LambdaNP;
        
        Gamma_W *= 1.0 - 3.0*alpha/4.0/(c2-s2)
                         *( obliqueS() - 2.0*c2*obliqueT() - (c2-s2)*obliqueU()/2.0/s2 )
                   //- 3.0*s2/2.0/(c2-s2)*DeltaGF()
                   - (1.0 + c2)/2.0/(c2-s2)*DeltaGF()
                   + (cHL1p + cHL2p + cHL3p + cHQ1p + cHQ2p)*ratio;
    } else
        if (obliqueS()!=0.0 || obliqueT()!=0.0 || obliqueU()!=0.0)
            throw std::runtime_error("NPEffective::GammaW(): The oblique corrections STU cannot be used with flag NotLinearizedNP=1");

    return Gamma_W;
}






