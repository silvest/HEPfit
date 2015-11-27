/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef SUSY_H
#define	SUSY_H

#include <gslpp.h>
#include <StandardModel.h>
#include "SUSYMatching.h"

class SUSYSpectrum; // forward reference to Spectrum class

/**
 * @class SUSY
 * @ingroup SUSY
 * @brief A base class for SUSY models. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details  
 */
class SUSY: public StandardModel {
public:

    static const int NSUSYvars = 10;
    static const std::string SUSYvars[NSUSYvars];
    
    /**
     * @brief Friend classes of SUSY class.
     */
    friend class SUSYSpectrum;
    //friend class FeynHiggsWrapper;

    /**
     * @brief A SUSY constructor.
     */
    SUSY();

    /**
     * @brief A SUSY destructor. 
     */
    ~SUSY();
    ///////////////////////////////////////////////////////////////////////////
    // Initialization

    virtual bool InitializeModel();

    virtual SUSYMatching* getMyMatching() const
    {
        return mySUSYMatching;
    }

    
    ///////////////////////////////////////////////////////////////////////////
    // Parameters 

    virtual bool Init(const std::map<std::string, double>& DPars);

    virtual bool PreUpdate();

    virtual bool Update(const std::map<std::string, double>& DPars);

    virtual bool PostUpdate();

    virtual bool CheckParameters(const std::map<std::string, double>& DPars);

    
    ///////////////////////////////////////////////////////////////////////////
    // Flags

    virtual bool setFlag(const std::string, const bool);

    ///////////////////////////////////////////////////////////////////////////
    // functions for the input parameters of SUSY model

    /**
     * @brief Gets the bino mass.
     * @return The bino mass.
     */
    gslpp::complex getM1() const
    {
        return m1;
    }

    /**
     * @brief Gets the wino mass.
     * @return The wino mass.
     */
    gslpp::complex getM2() const
    {
        return m2;
    }

    /**
     * @brief Gets the gluino mass.
     * @return The gluino mass.
     */
    double getM3() const
    {
        return m3;
    }

    /**
     * @brief Gets the charged-Higgs mass at tree-level.
     * @return The charged-Higgs mass at tree-level.
     */
    double getMHptree() const
    {
        return mHptree;
    }

    /**
     * @brief Gets the @f$\mu@f$ parameter in the superpotential.
     * @return The @f$\mu@f$ parameter in the superpotential.
     */
    gslpp::complex getMuH() const
    {
        return muH;
    }

    /**
     * @brief Gets @f$\tan\beta@f$.
     * @return @f$\tan\beta@f$.
     */
    double getTanb() const
    {
        return tanb;
    }

    /**
     * @brief Gets the scale of the input parameters.
     * @return The scale of the input parameters.
     */
    double getQ_SUSY() const
    {
        return Q_SUSY;
    }


    ///////////////////////////////////////////////////////////////////////////
    // functions for the parameters in the Higgs sector.

    /**
     * @return The Higgs VEV associated with the down-type quarks.
     */
    double v1() const;

    /**
     * @return The Higgs VEV associated with the up-type quarks.
     */
    double v2() const;
    
    /**
     * @brief Gets @f$\sin\beta@f$.
     * @return @f$\sin\beta@f$.
     */
    double getSinb() const
    {
        return sinb;
    }

    /**
     * @brief Gets @f$\cos\beta@f$.
     * @return @f$\cos\beta@f$.
     */
    double getCosb() const
    {
        return cosb;
    }

    /**
     * @brief Gets the sine of the effective mixing angle for the CP-even neutral Higgs bosons. 
     * @return
     */
    gslpp::complex getSaeff() const
    {
        return saeff;
    }

    /**
     * @brief Gets the light Higgs mass.
     * @return The light Higgs mass.
     */
    virtual double getMHl() const
    {
        return mh[0];
    }

    /**
     * @brief Gets the heavy Higgs mass.
     * @return The heavy Higgs mass.
     */
    double getMHh() const
    {
        return mh[1];
    }

    /**
     * @brief Gets the pseudo-scalar Higgs mass.
     * @return The pseudo-scalar Higgs mass. 
     */
    double getMHa() const
    {
        return mh[2];
    }

    /**
     * @brief Gets the charged Higgs mass.
     * @return The charged Higgs mass.
     */
    double getMHp() const
    {
        return mh[3];
    }


    ///////////////////////////////////////////////////////////////////////////
    // functions for the parameters in the gaugino sector.

    /**
     * @brief Gets the gluino mass obtained from FHGetPara().
     * @return The gluino mass obtained from FHGetPara().
     */
    double getMGl() const;

    
    ///////////////////////////////////////////////////////////////////////////
    // functions for the parameters in the chargino sector.

    /**
     * @brief Gets the chargino masses.
     * @return The chargino masses.
     */
    gslpp::vector<double> getMch() const
    {
        return mch;
    }

    /**
     * @brief Gets the rotation matrix for negative charginos.
     * @return The rotation matrix for negative charginos.
     */
    gslpp::matrix<gslpp::complex> getU() const
    {
        return U;
    }

    /**
     * @brief Gets the rotation matrix for positive charginos.
     * @return The rotation matrix for positive charginos.
     */
    gslpp::matrix<gslpp::complex> getV() const
    {
        return V;
    }


    ///////////////////////////////////////////////////////////////////////////
    // functions for the parameters in the neutralino sector.
    
    /**
     * @brief Gets the neutralino masses.
     * @return The neutralino masses. 
     */
    gslpp::vector<double> getMneu() const
    {
        return mneu;
    }

    /**
     * @brief Gets the rotation matrix for neutralinos.
     * @return The rotation matrix for neutralinos.
     */
    gslpp::matrix<gslpp::complex> getN() const
    {
        return N;
    }


    ///////////////////////////////////////////////////////////////////////////
    // functions for the parameters in the squark sector.

    gslpp::matrix<gslpp::complex> getMsQhat2() const
    {
        return msQhat2;
    }

    gslpp::matrix<gslpp::complex> getMsUhat2() const
    {
        return msUhat2;
    }

    gslpp::matrix<gslpp::complex> getMsDhat2() const
    {
        return msDhat2;
    }

    /**
     * @brief Gets the trilinear-coupling matrix for up-type squarks.
     * @return The trilinear-coupling matrix for up-type squarks.
     */
    gslpp::matrix<gslpp::complex> getTUhat() const
    {
        return TUhat;
    }

    /**
     * @brief Gets the trilinear-coupling matrix for down-type squarks.
     * @return The trilinear-coupling matrix for down-type squarks.
     */
    gslpp::matrix<gslpp::complex> getTDhat() const
    {
        return TDhat;
    }

    /**
     * @brief Gets the up-type squark mass squared.
     * @return The up-type squark mass squared.
     */
    gslpp::vector<double> getMsu2() const
    {
        return m_su2;
    }

    /**
     * @brief Gets the down-type squark mass squared.
     * @return The down-type squark mass squared. 
     */
    gslpp::vector<double> getMsd2() const
    {
        return m_sd2;
    }

    /**
     * @brief Gets the down-type squark mass squared with the
     * @f$\Delta_b@f$ corrections in the off-diagonal entries. 
     * @return The down-type squark mass squared.
     */
    gslpp::vector<double> getMsdresum2() const
    {
        return m_sdresum2;
    }

    /**
     * @brief Gets the rotation matrix for up-type squarks.
     * @return The rotation matrix for up-type squarks.
     */
    gslpp::matrix<gslpp::complex> getRu() const
    {
        return Ru;
    }

    /**
     * @brief Gets the rotation matrix for down-type squarks.
     * @return The rotation matrix for down-type squarks.
     */
    gslpp::matrix<gslpp::complex> getRd() const
    {
        return Rd;
    }

    /**
     * @brief Gets the rotation matrix for down-type squarks with the
     * @f$\Delta_b@f$ corrections in the off-diagonal entries. 
     * @return The rotation matrix for down-type squarks.
     */
    gslpp::matrix<gslpp::complex> getRdresum() const
    {
        return Rdresum;
    }

    ///////////////////////////////////////////////////////////////////////////
    // functions for the parameters in the slepton sector.

    gslpp::matrix<gslpp::complex> getMsLhat2() const
    {
        return msLhat2;
    }

    gslpp::matrix<gslpp::complex> getMsNhat2() const
    {
        return msNhat2;
    }

    gslpp::matrix<gslpp::complex> getMsEhat2() const
    {
        return msEhat2;
    }

    /**
     * @brief Gets the trilinear-coupling matrix for sneutrinos.
     * @return The trilinear-coupling matrix for sneutrinos.
     */
    gslpp::matrix<gslpp::complex> getTNhat() const
    {
        return TNhat;
    }

    /**
     * @brief Gets the trilinear-coupling matrix for charged sleptons.
     * @return The trilinear-coupling matrix for charged sleptons.
     */
    gslpp::matrix<gslpp::complex> getTEhat() const
    {
        return TEhat;
    }

    /**
     * @brief Gets the sneutrino mass squared.
     * @return The sneutrino mass squared.
     */
    gslpp::vector<double> getMsn2() const
    {
        return m_sn2;
    }

    /**
     * @brief Gets the charged slepton mass squared.
     * @return The charged slepton mass squared.
     */
    gslpp::vector<double> getMse2() const
    {
        return m_se2;
    }

    /**
     * @brief Gets the rotation matrix for sneutrinos.
     * @return The rotation matrix for sneutrinos.
     */
    gslpp::matrix<gslpp::complex> getRn() const
    {
        return Rn;
    }

    /**
     * @brief Gets the rotation matrix for charged sleptons.
     * @return The rotation matrix for charged sleptons.
     */
    gslpp::matrix<gslpp::complex> getRl() const
    {
        return Rl;
    }

    
    ///////////////////////////////////////////////////////////////////////////
    // functions for SM fermions

    double Mq_Q(const quark q) const
    {
        switch (q) {
            case UP:
            case CHARM:
            case TOP:
                return mu_Q[(int)(q - UP)/2];
            case DOWN:
            case STRANGE:
            case BOTTOM:
                return md_Q[((int)(q - DOWN))/2];
            default:
                throw std::runtime_error("SUSY::Mq_Q(): Error!");
        }
    }

    double Ml_Q(const lepton l) const
    {
        switch (l) {
            case ELECTRON:
            case MU:
            case TAU:
                return me_Q[(int)(l - ELECTRON)/2];
            case NEUTRINO_1:
            case NEUTRINO_2:
            case NEUTRINO_3:
                return mn_Q[((int)(l - NEUTRINO_1))/2];
            default:
                throw std::runtime_error("SUSY::Ml_Q(): Error!");
        }
    }

protected:
    virtual void setParameter(const std::string name , const double& value);
    virtual void SetTanb(const double tanb);
    virtual void computeYukawas();

    SUSYSpectrum* mySUSYSpectrum;

    // model parameters at scale Q
    gslpp::complex m1, m2, muH;
    double m3, mHptree, tanb, Q_SUSY;

    // sinb and cosb computed with setTanb()
    double sinb, cosb;

    // soft-breaking parameters associated with squark and slepton mass terms and
    // trilinear couplings in the SCKM basis, which will be set with SetSoftTerms()
    gslpp::matrix<gslpp::complex> msQhat2, msUhat2, msDhat2, msLhat2, msNhat2, msEhat2;
    gslpp::matrix<gslpp::complex> TUhat, TDhat, TNhat, TEhat;

    // soft-breaking parameters computed with FeynHiggs
    double mHp, mh[4];
    gslpp::complex saeff;
    gslpp::vector<double> mch, mneu;
    gslpp::vector<double> m_su2, m_sd2, m_sdresum2, m_sn2, m_se2;

    // rotation matrices
    gslpp::matrix<gslpp::complex> U, V, N, Ru, Rd, Rdresum, Rn, Rl;

    // quark and lepton masses at scale Q, computed in setYukawas()
    double mu_Q[3], md_Q[3], me_Q[3], mn_Q[3];
    
    ///////////////////////////////////////////////////////////////////////////
private:    
    bool flag_h, flag_g, flag_ch, flag_ne;
    SUSYMatching* mySUSYMatching;

};

#endif	/* SUSY_H */

