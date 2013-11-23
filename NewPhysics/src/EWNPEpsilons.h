/*
 * Copyright (C) 2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef EWNPEPSILONS_H
#define	EWNPEPSILONS_H

#include <EWSM.h>
#include <EWepsilons.h>

/**
 * @class EWNPEpsilons
 * @brief A class for ...
 * @ingroup NewPhysics
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details
 */
class EWNPEpsilons : public EWSM {
public:

    EWNPEpsilons(const StandardModel& SM_i);


    ////////////////////////////////////////////////////////////////////////

    double Mw_NPEpsilons() const;


    ////////////////////////////////////////////////////////////////////////

    /**
     * @param[in] l Name of a lepton.
     * @return The effective coupling @f$\rho_Z^l@f$.
     */
    virtual complex rhoZ_l(const StandardModel::lepton l) const;

    /**
     * @param[in] q Name of a quark.
     * @return The effective coupling neutral-current interactions @f$\rho_Z^q@f$.
     */
    virtual complex rhoZ_q(const StandardModel::quark q) const;

    /**
     * @param[in] l Name of a lepton.
     * @return The effective coupling neutral-current interactions @f$\kappa_Z^l@f$.
     */
    virtual complex kappaZ_l(const StandardModel::lepton l) const;

    /**
     * @param[in] q Name of a quark.
     * @return The effective coupling neutral-current interactions @f$\kappa_Z^q@f$.
     */
    virtual complex kappaZ_q(const StandardModel::quark q) const;

    /**
     * @param[in] l Name of a lepton.
     * @return The effective vector coupling for neutral-current interactions @f$g_V^l@f$.
     */
    virtual complex gVl(const StandardModel::lepton l) const;

    /**
     * @param[in] q Name of a quark.
     * @return The effective vector coupling for neutral-current interactions @f$g_V^q@f$.
     */
    virtual complex gVq(const StandardModel::quark q) const;

    /**
     * @param[in] l Name of a lepton.
     * @return The effective axial-vector coupling for neutral-current interactions @f$g_A^l@f$.
     */
    virtual complex gAl(const StandardModel::lepton l) const;

    /**
     * @param[in] q Name of a quark.
     * @return The effective axial-vector coupling for neutral-current interactions @f$g_A^q@f$.
     */
    virtual complex gAq(const StandardModel::quark q) const;


    ////////////////////////////////////////////////////////////////////////
private:
    const StandardModel& SM;
    EWepsilons* myEWepsilons;
    
};

#endif	/* EWNPEPSILONS_H */

