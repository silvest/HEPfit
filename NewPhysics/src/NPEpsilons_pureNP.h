/*
 * Copyright (C) 2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef NPEPSILONS_PURENP_H
#define	NPEPSILONS_PURENP_H

#include "NPbase.h"

/**
 * @class NPEpsilons_pureNP
 * @brief A class for new physics in the form of contributions to the \f$\varepsilon_{1,2,3,b}\f$ parameters.
 * Only new physics contributions are parameterized.
 * @ingroup NewPhysics
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details
 */
class NPEpsilons_pureNP : public NPbase {
public:
    static const int NEPSILONpureNPvars = 4;
    static const std::string EPSILONpureNPvars[NEPSILONpureNPvars];

    NPEpsilons_pureNP();

    virtual std::string ModelName() const
    {
        return "NPEpsilons_pureNP";
    }

    virtual bool InitializeModel();
    virtual void setEWSMflags(EWSM& myEWSM);

    virtual bool Init(const std::map<std::string, double>& DPars);
    virtual bool Update(const std::map<std::string, double>& DPars);
    virtual bool CheckParameters(const std::map<std::string, double>& DPars);

    virtual bool setFlag(const std::string, const bool&);
    virtual bool CheckFlags() const;


    ////////////////////////////////////////////////////////////////////////
    
    /**
     * @return the pure new physics contribution to \f$\varepsilon_1\f$
     */
    virtual double epsilon1() const;

    /**
     * @return the pure new physics contribution to \f$\varepsilon_2\f$
     */
    virtual double epsilon2() const;
    
    /**
     * @return the pure new physics contribution to \f$\varepsilon_3\f$
     */
    virtual double epsilon3() const;
    
    /**
     * @return the pure new physics contribution to \f$\varepsilon_b\f$
     */
    virtual double epsilonb() const;


    ////////////////////////////////////////////////////////////////////////

    /**
     * @return the \f$W\f$-boson mass in GeV
     */
    virtual double Mw() const;

    /**
     * @return the (square of the) cosine of the weak angle in the On-mass-shell renormalization scheme,
     *  \f$\cos^2{\theta_W}=\frac{M_W^2}{M_Z^2}\f$
     */
    virtual double cW2() const;

    /**
     * @return the (square of the) sine of the weak angle in the On-mass-shell renormalization scheme,
     *  \f$\sin^2{\theta_W}=1-\frac{M_W^2}{M_Z^2}\f$
     */
    virtual double sW2() const;

    /**
     * @return the total width of the \f$W\f$ boson in GeV [NOT IMPLEMENTED YET]
     */
    virtual double GammaW() const;


    ////////////////////////////////////////////////////////////////////////

    virtual double deltaGVl(StandardModel::lepton l) const;

    virtual double deltaGVq(StandardModel::quark q) const;

    virtual double deltaGAl(StandardModel::lepton l) const;

    virtual double deltaGAq(StandardModel::quark q) const;


    ////////////////////////////////////////////////////////////////////////
protected:
    double deltaEps_1, deltaEps_2, deltaEps_3, deltaEps_b;
    virtual void setParameter(const std::string name, const double& value);

};

#endif	/* NPEPSILONS_PURENP_H */

