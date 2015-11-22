/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef EWSMONELOOPEW_HV_H
#define	EWSMONELOOPEW_HV_H

#include <PVfunctions.h>
#include "StandardModel.h"


/**
 * @class EWSMOneLoopEW_HV
 * @ingroup StandardModel
 * @brief A class for @f$O(\alpha)@f$ one-loop corrections to the %EW
 * precision observables in the 't Hooft-Feynman gauge.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class EWSMOneLoopEW_HV {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    EWSMOneLoopEW_HV(const StandardModel& SM_i);


    ////////////////////////////////////////////////////////////////////////    

    /**
     * @param[in] l name of lepton
     * @return mass of lepton
     */
    double ml(const StandardModel::lepton l) const
    {
        return SM.getLeptons(l).getMass();
    }

    /**
     * @param[in] q name of quark
     * @param[in] mu renormalization scale
     * @param[in] order (=LO, NLO, NNLO, FULLNLO, FULLNNLO)
     * @return the MSbar mass of u, d, s, c, b or the pole mass of t
     */
    double mq(const QCD::quark q, const double mu, const orders order = FULLNLO) const
    {
        switch (q) {
            case QCD::UP:
            case QCD::DOWN:
            case QCD::STRANGE:
                return SM.Mrun(mu, SM.getQuarks(q).getMass_scale(), SM.getQuarks(q).getMass(), order);
            case QCD::CHARM:
            case QCD::BOTTOM:
                return SM.Mrun(mu, SM.getQuarks(q).getMass(), order);
            case QCD::TOP:
                return SM.getMtpole(); // the pole mass
            default:
                throw std::runtime_error("Error in EWSMOneLoopEW_HV::mq()");
        }
    }

    /**
     * @param[in] l name of lepton
     * @param[in] Mw the W-boson mass
     * @return the tree-level vector coupling for Z->l lbar
     *
     * @attention depends on sW2
     */
    double vl(const StandardModel::lepton l, const double Mw) const
    {
        double sW2 = 1.0 - Mw * Mw / SM.getMz() / SM.getMz();
        return ( al(l) - 2.0 * SM.getLeptons(l).getCharge() * sW2);
    }

    /**
     * @param[in] q name of quark
     * @param[in] Mw the W-boson mass
     * @return the tree-level vector coupling for Z->q qbar
     * 
     * @attention depends on sW2
     */
    double vq(const QCD::quark q, const double Mw) const
    {
        double sW2 = 1.0 - Mw * Mw / SM.getMz() / SM.getMz();
        return ( aq(q) - 2.0 * SM.getQuarks(q).getCharge() * sW2);
    }

    /**
     * @param[in] l name of lepton
     * @return the tree-level axial-vector coupling for Z->l lbar
     */
    double al(const StandardModel::lepton l) const
    {
        return ( SM.getLeptons(l).getIsospin());
    }

    /**
     * @param[in] q name of quark
     * @return the tree-level axial-vector coupling for Z->q qbar
     */
    double aq(const QCD::quark q) const
    {
        return ( SM.getQuarks(q).getIsospin());
    }


    ////////////////////////////////////////////////////////////////////////    

    /**
     * @param[in] mu renormalization scale
     * @param[in] s momentum-squared
     * @param[in] Mw the W-boson mass
     * @return bosonic contribution to the self-energy function of W boson
     */
    gslpp::complex SigmaWW_bos(const double mu, const double s, const double Mw) const;

    /**
     * @param[in] mu renormalization scale
     * @param[in] muForMq renormalization scale for the running quark mass
     * @param[in] s momentum-squared
     * @return fermionic contribution to the self-energy function of W boson
     */
    gslpp::complex SigmaWW_fer(const double mu, const double muForMq, const double s) const;

    /**
     * @param[in] mu renormalization scale
     * @param[in] s momentum-squared
     * @param[in] Mw the W-boson mass
     * @return bosonic contribution to the self-energy function of Z boson
     */
    gslpp::complex SigmaZZ_bos(const double mu, const double s, const double Mw) const;

    /**
     * @param[in] mu renormalization scale
     * @param[in] muForMq renormalization scale for the running quark mass
     * @param[in] s momentum-squared
     * @param[in] Mw the W-boson mass
     * @return fermionic contribution to the self-energy function of Z boson
     */
    gslpp::complex SigmaZZ_fer(const double mu, const double muForMq, const double s,
            const double Mw) const;

    /**
     * @param[in] mu renormalization scale
     * @param[in] s momentum-squared
     * @param[in] Mw the W-boson mass
     * @return bosonic contribution to the self-energy function of photon 
     */
    gslpp::complex SigmaGammaGamma_bos(const double mu, const double s, const double Mw) const;

    /**
     * @param[in] mu renormalization scale
     * @param[in] muForMq renormalization scale for the running quark mass
     * @param[in] s momentum-squared
     * @return fermionic contribution to the self-energy function of photon 
     */
    gslpp::complex SigmaGammaGamma_fer(const double mu, const double muForMq, const double s) const;

    /**
     * @param[in] mu renormalization scale
     * @param[in] s momentum-squared
     * @param[in] Mw the W-boson mass
     * @return bosonic contribution to the self-energy function of photon 
     */
    gslpp::complex PiGammaGamma_bos(const double mu, const double s, const double Mw) const;

    /**
     * @param[in] mu renormalization scale
     * @param[in] s momentum-squared
     * @param[in] l name of lepton
     * @return contribution to the self-energy function of photon from lepton l
     */
    gslpp::complex PiGammaGamma_fer_l(const double mu, const double s, const StandardModel::lepton l) const;

    /**
     * @param[in] mu renormalization scale
     * @param[in] muForMq renormalization scale for the running quark mass
     * @param[in] s momentum-squared
     * @param[in] q name of quark
     * @return contribution to the self-energy function of photon from quark q
     */
    gslpp::complex PiGammaGamma_fer_q(const double mu, const double muForMq,
            const double s, const QCD::quark q) const;

    /**
     * @param[in] mu renormalization scale
     * @param[in] muForMq renormalization scale for the running quark mass
     * @param[in] s momentum-squared
     * @return fermionic contribution to the self-energy function of photon 
     */
    gslpp::complex PiGammaGamma_fer(const double mu, const double muForMq, const double s) const;

    /**
     * @param[in] mu renormalization scale
     * @param[in] s momentum-squared
     * @param[in] Mw the W-boson mass
     * @return bosonic contribution to the self-energy function of the Z-gamma mixing
     */
    gslpp::complex SigmaZgamma_bos(const double mu, const double s, const double Mw) const;

    /**
     * @param[in] mu renormalization scale
     * @param[in] muForMq renormalization scale for the running quark mass
     * @param[in] s momentum-squared
     * @param[in] Mw the W-boson mass
     * @return fermionic contribution to the self-energy function of the Z-gamma mixing
     */
    gslpp::complex SigmaZgamma_fer(const double mu, const double muForMq,
            const double s, const double Mw) const;


    ////////////////////////////////////////////////////////////////////////      

    /**
     * @param[in] s momentum-squared
     * @param[in] m1 mass 
     * @param[in] m2 mass
     * @return one-loop function F(s, m1, m2)
     */
    gslpp::complex F_Hollik(const double s, const double m1, const double m2) const;

    /**
     * @param[in] muIR renormalization scale for a possible IR divergence
     * @param[in] s momentum-squared
     * @param[in] m1 mass 
     * @param[in] m2 mass
     * @return dF(s, m1, m2)/ds
     */
    gslpp::complex Fprime_Hollik(const double muIR, const double s, const double m1, const double m2) const;

    /**
     * @param[in] mu renormalization scale
     * @param[in] s momentum-squared
     * @param[in] Mw the W-boson mass
     * @return bosonic contribution to the self-energy function of W boson
     */
    gslpp::complex SigmaWW_bos_Hollik(const double mu, const double s, const double Mw) const;

    /**
     * @param[in] mu renormalization scale
     * @param[in] s momentum-squared
     * @param[in] Mw the W-boson mass
     * @return bosonic contribution to the self-energy function of Z boson
     */
    gslpp::complex SigmaZZ_bos_Hollik(const double mu, const double s, const double Mw) const;

    /**
     * @param[in] mu renormalization scale
     * @param[in] s momentum-squared
     * @param[in] Mw the W-boson mass
     * @return bosonic contribution to the self-energy function of photon
     */
    gslpp::complex SigmaGammaGamma_bos_Hollik(const double mu, const double s, const double Mw) const;

    /**
     * @param[in] mu renormalization scale
     * @param[in] s momentum-squared
     * @param[in] Mw the W-boson mass
     * @return bosonic contribution to the self-energy function of photon
     */
    gslpp::complex PiGammaGamma_bos_Hollik(const double mu, const double s, const double Mw) const;

    /**
     * @param[in] mu renormalization scale
     * @param[in] s momentum-squared
     * @param[in] Mw the W-boson mass
     * @return bosonic contribution to the self-energy function of the Z-gamma mixing
     */
    gslpp::complex SigmaZgamma_bos_Hollik(const double mu, const double s, const double Mw) const;


    ////////////////////////////////////////////////////////////////////////      

private:
    const StandardModel& SM; ///< A reference to an object of type StandardModel.
    const PVfunctions PV; ///< An object of type PVfunctions.

};

#endif	/* EWSMONELOOPEW_HV_H */

