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
 * @brief A class for new physics with the @f$\varepsilon@f$ parameters, where
 * only NP contributions are parameterized.
 * @ingroup NewPhysics
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details
 */
class NPEpsilons_pureNP : public NPbase {
public:
    static const int NEPSILONpureNPvars = 4;
    static const std::string EPSILONpureNPvars[NEPSILONpureNPvars];
    //static const int NEPSILONpureNPflags = 0;
    //static const std::string EPSILONpureNPflags[NEPSILONpureNPflags];

    NPEpsilons_pureNP();

    virtual std::string ModelName() const
    {
        return "NPEpsilons_pureNP";
    }

    virtual bool Update(const std::map<std::string, double>& DPars);
    virtual bool Init(const std::map<std::string, double>& DPars);
    virtual bool CheckParameters(const std::map<std::string, double>& DPars);

    virtual bool InitializeModel();
    virtual void setEWSMflags(EWSM& myEWSM);

    virtual bool setFlag(const std::string, const bool&);
    virtual bool CheckFlags() const;


    ////////////////////////////////////////////////////////////////////////

    virtual double epsilon1() const;

    virtual double epsilon2() const;

    virtual double epsilon3() const;

    virtual double epsilonb() const;


    ////////////////////////////////////////////////////////////////////////

    /**
     * @return the W boson mass
     */
    virtual double Mw() const;

    /**
     * @return Mw^2/Mz^2
     */
    virtual double cW2() const;

    /**
     * @return 1-Mw^2/Mz^2
     */
    virtual double sW2() const;

    virtual double deltaGVl(StandardModel::lepton l) const;

    virtual double deltaGVq(StandardModel::quark q) const;

    virtual double deltaGAl(StandardModel::lepton l) const;

    virtual double deltaGAq(StandardModel::quark q) const;

    /**
     * @return the total width of the W boson
     */
    virtual double GammaW() const;


    ////////////////////////////////////////////////////////////////////////

protected:
    double deltaEps_1, deltaEps_2, deltaEps_3, deltaEps_b;
    virtual void setParameter(const std::string name, const double& value);


    ////////////////////////////////////////////////////////////////////////

private:

    ////////////////////////////////////////////////////////////////////////

private:

};

#endif	/* NPEPSILONS_PURENP_H */

