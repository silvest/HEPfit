/* 
 * Copyright (C) 2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include <stdexcept>
#include "NPEffective.h"


NPEffective::NPEffective() 
: NPbase()
{
}

bool NPEffective::PostUpdate() {
    bool NPbaseup = NPbase::PostUpdate();
    trueNPbase = *this;
    return (NPbaseup);
}

////////////////////////////////////////////////////////////////////////

double NPEffective::getCoeff(const std::string name) const
{
    if (name.compare("cWB") == 0)
        return cWB;
    else if (name.compare("cH") == 0)
        return cH;
    else if (name.compare("cL1L1") == 0)
        return cL1L1;
    else if (name.compare("cL1L2") == 0)
        return cL1L2;
    else if (name.compare("cL1L3") == 0)
        return cL1L3;
    else if (name.compare("cL2L2") == 0)
        return cL2L2;
    else if (name.compare("cL2L3") == 0)
        return cL2L3;
    else if (name.compare("cL3L3") == 0)
        return cL3L3;
    else if (name.compare("cHL1p") == 0)
        return cHL1p;
    else if (name.compare("cHL2p") == 0)
        return cHL2p;
    else if (name.compare("cHL3p") == 0)
        return cHL3p;
    else if (name.compare("cHQ1p") == 0)
        return cHQ1p;
    else if (name.compare("cHQ2p") == 0)
        return cHQ2p;
    else if (name.compare("cHQ3p") == 0)
        return cHQ3p;
    else if (name.compare("cHL1") == 0)
        return cHL1;
    else if (name.compare("cHL2") == 0)
        return cHL2;
    else if (name.compare("cHL3") == 0)
        return cHL3;
    else if (name.compare("cHQ1") == 0)
        return cHQ1;
    else if (name.compare("cHQ2") == 0)
        return cHQ2;
    else if (name.compare("cHQ3") == 0)
        return cHQ3;
    else if (name.compare("cHE1") == 0)
        return cHE1;
    else if (name.compare("cHE2") == 0)
        return cHE2;
    else if (name.compare("cHE3") == 0)
        return cHE3;
    else if (name.compare("cHU1") == 0)
        return cHU1;
    else if (name.compare("cHU2") == 0)
        return cHU2;
    else if (name.compare("cHU3") == 0)
        return cHU3;
    else if (name.compare("cHD1") == 0)
        return cHD1;
    else if (name.compare("cHD2") == 0)
        return cHD2;
    else if (name.compare("cHD3") == 0)
        return cHD3;
    else
        throw std::runtime_error("NPEffective::getCoeff(): " + name + " is not available");
}


////////////////////////////////////////////////////////////////////////

double NPEffective::v() const
{
    //return ( sqrt( (1.0 - (cL1L2 - cHL1p - cHL2p)/sqrt(2.0)/GF/LambdaNP/LambdaNP)
    //               /sqrt(2.0)/GF ) );

    /* use the tree-level relation */
    return trueSM.v();
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

double NPEffective::GammaW() const
{
    double Gamma_W = trueSM.GammaW();

    double alpha = alphaMz();
    double c2 = trueSM.cW2();
    double s2 = trueSM.sW2();
    double ratio = v()*v()/LambdaNP/LambdaNP;

    Gamma_W *= 1.0 - 3.0*alpha/4.0/(c2-s2)
                     *( obliqueS() - 2.0*c2*obliqueT() - (c2-s2)*obliqueU()/2.0/s2 )
               //- 3.0*s2/2.0/(c2-s2)*DeltaGF()
               - (1.0 + c2)/2.0/(c2-s2)*DeltaGF()
               + (cHL1p + cHL2p + cHL3p + cHQ1p + cHQ2p)*ratio;

    return Gamma_W;
}


////////////////////////////////////////////////////////////////////////     

double NPEffective::obliqueS() const 
{
    double ratio = v()*v()/LambdaNP/LambdaNP;
    double sW_SM = sqrt(trueSM.sW2()); /* This has to be the SM value. */
    double cW_SM = sqrt(trueSM.cW2()); /* This has to be the SM value. */

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

double NPEffective::deltaGV_f(const Particle p) const
{
    return ( deltaGL_f_tmp(p) + deltaGR_f_tmp(p) + trueNPbase.deltaGV_f(p) );
}


double NPEffective::deltaGA_f(const Particle p) const
{
    return ( deltaGL_f_tmp(p) - deltaGR_f_tmp(p) + trueNPbase.deltaGA_f(p) );
}

////////////////////////////////////////////////////////////////////////

double NPEffective::deltaGL_f_tmp(const Particle p) const
{
    double ratio = v()*v()/LambdaNP/LambdaNP;
    if(p.is("NEUTRINO_1"))
            return ( (cHL1p - cHL1)/2.0*ratio );
    if(p.is("NEUTRINO_2"))
            return ( (cHL2p - cHL2)/2.0*ratio );
    if(p.is("NEUTRINO_3"))
            return ( (cHL3p - cHL3)/2.0*ratio );
    if(p.is("ELECTRON"))
            return ( - (cHL1p + cHL1)/2.0*ratio );
    if(p.is("MU"))
            return ( - (cHL2p + cHL2)/2.0*ratio );
    if(p.is("TAU"))
            return ( - (cHL3p + cHL3)/2.0*ratio );
    if(p.is("UP"))
           return ( (cHQ1p - cHQ1)/2.0*ratio );
    if(p.is("CHARM"))
           return ( (cHQ2p - cHQ2)/2.0*ratio );
    if(p.is("TOP"))
            return 0.0;
    if(p.is("DOWN"))
            return ( - (cHQ1p + cHQ1)/2.0*ratio );
    if(p.is("STRANGE"))
            return ( - (cHQ2p + cHQ2)/2.0*ratio );
    if(p.is("BOTTOM"))
            return ( - (cHQ3p + cHQ3)/2.0*ratio );
    throw std::runtime_error("Error in NPEffective::deltaGL_f()");
}

double NPEffective::deltaGR_f_tmp(const Particle p) const
{
    double ratio = v()*v()/LambdaNP/LambdaNP;
    if(p.is("NEUTRINO_1"))
            return ( 0. );
    if(p.is("NEUTRINO_2"))
            return ( 0. );
    if(p.is("NEUTRINO_3"))
            return ( 0. );
    if(p.is("ELECTRON"))
            return ( - cHE1/2.0*ratio );
    if(p.is("MU"))
            return ( - cHE2/2.0*ratio );
    if(p.is("TAU"))
            return ( - cHE3/2.0*ratio );
    if(p.is("UP"))
           return ( - cHU1/2.0*ratio );
    if(p.is("CHARM"))
           return ( - cHU2/2.0*ratio );
    if(p.is("TOP"))
            return 0.0;
    if(p.is("DOWN"))
            return ( - cHD1/2.0*ratio );
    if(p.is("STRANGE"))
            return ( - cHD2/2.0*ratio );
    if(p.is("BOTTOM"))
            return ( - cHD3/2.0*ratio );
    throw std::runtime_error("Error in NPEffective::deltaGL_f()");
}


