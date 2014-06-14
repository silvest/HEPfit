/*
 * Copyright (C) 2014 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef NPEFFECTIVEGIMR_H
#define	NPEFFECTIVEGIMR_H

#include <string.h>
#include <stdexcept>
#include <gslpp.h>
#include "NPbase.h"

class NPEffectiveGIMR : public NPbase {
public:

    static const int NNPEffectiveGIMRVars = 89;  /* Only a subset of the coefficients! */
    static const std::string NPEffectiveGIMRVars[NNPEffectiveGIMRVars];

    NPEffectiveGIMR();

    virtual bool PostUpdate();

    /**
     * @brief @copybrief Model::CheckParameters()
     * @copydetails Model::CheckParameters()
     */
    virtual bool CheckParameters(const std::map<std::string, double>& DPars);

    
    ////////////////////////////////////////////////////////////////////////

    virtual double DeltaGF() const;
    virtual double obliqueS() const;
    virtual double obliqueT() const;
    virtual double obliqueU() const;
    virtual double Mw() const;
    virtual double GammaW() const;
    virtual double deltaGV_f(const Particle p) const;
    virtual double deltaGA_f(const Particle p) const;
    double deltaGL_f(const Particle p) const;
    double deltaGR_f(const Particle p) const;

    ////////////////////////////////////////////////////////////////////////
protected:

    /**
     * @brief @copybrief Model::setParameter()
     * @copydetails Model::setParameter()
     */
    virtual void setParameter(const std::string name, const double& value);

    double CW; ///< The dimension-6 operator coefficient \f$C_{W}\f$.
    double CHG;
    double CHW;
    double CHB;
    double CHWB;
    double CHD;
    double CHbox;
    double CH;
    double CHL1_11, CHL1_12, CHL1_13, CHL1_22, CHL1_23, CHL1_33;
    double CHL3_11, CHL3_12, CHL3_13, CHL3_22, CHL3_23, CHL3_33;
    double CHe_11, CHe_12, CHe_13, CHe_22, CHe_23, CHe_33;
    double CHQ1_11, CHQ1_12, CHQ1_13, CHQ1_22, CHQ1_23, CHQ1_33;
    double CHQ3_11, CHQ3_12, CHQ3_13, CHQ3_22, CHQ3_23, CHQ3_33;
    double CHu_11, CHu_12, CHu_13, CHu_22, CHu_23, CHu_33;
    double CHd_11, CHd_12, CHd_13, CHd_22, CHd_23, CHd_33;
    double CeH_11r, CeH_12r, CeH_13r, CeH_22r, CeH_23r, CeH_33r;
    double CeH_11i, CeH_12i, CeH_13i, CeH_22i, CeH_23i, CeH_33i;
    double CuH_11r, CuH_12r, CuH_13r, CuH_22r, CuH_23r, CuH_33r;
    double CuH_11i, CuH_12i, CuH_13i, CuH_22i, CuH_23i, CuH_33i;
    double CdH_11r, CdH_12r, CdH_13r, CdH_22r, CdH_23r, CdH_33r;
    double CdH_11i, CdH_12i, CdH_13i, CdH_22i, CdH_23i, CdH_33i;
    double CLL_1221, CLL_2112;
    double LambdaNP;

    double v2_over_LambdaNP2;
    double cW_tree, sW_tree;
    double cW2_tree, sW2_tree;


    ////////////////////////////////////////////////////////////////////////
private:

};

#endif	/* NPEFFECTIVEGIMR_H */

