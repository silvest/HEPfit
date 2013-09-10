/* 
 * Copyright (C) 2012-2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef SUSY_H
#define	SUSY_H

#include <gslpp.h>
#include <StandardModel.h>
#include <EWSM.h>
#include "SUSYMatching.h"

using namespace gslpp;

class FeynHiggsWrapper; // forward reference to FeynHiggsWrapper class
class SUSYSpectrum; // forward reference to Spectrum class
class EWSUSY; // forward reference to EWSUSY class

/**
 * @addtogroup SUSY
 * @brief A project for a basis of SUSY models. 
 * @{
 */

/**
 * @class SUSY
 * @brief A base class for SUSY models. 
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details  
 */
class SUSY: public StandardModel {
public:

    static const int NSUSYvars = 10;
    static const std::string SUSYvars[NSUSYvars];

    static const int NSUSYFlags = 4;
    static const std::string SUSYFlags[NSUSYFlags];
    
    /**
     * @brief Friend classes of SUSY class.
     */
    friend class SUSYSpectrum;
    friend class FeynHiggsWrapper;

    /**
     * @brief A SUSY constructor. 
     */
    SUSY();
    
    virtual std::string ModelName() const
    {
        return "SUSY";
    }

    
    ///////////////////////////////////////////////////////////////////////////
    // Initialization and Matching

    virtual bool InitializeModel();

    virtual void SetEWSMflags(EWSM& myEWSM);

    virtual SUSYMatching* GetMyMatching() const
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

    virtual bool SetFlag(const std::string, const bool&);

    bool IsFlag_h() const
    {
        return flag_h;
    }

    bool IsFlag_g() const
    {
        return flag_g;
    }

    bool IsFlag_ch() const
    {
        return flag_ch;
    }

    bool IsFlag_ne() const
    {
        return flag_ne;
    }


    ///////////////////////////////////////////////////////////////////////////

    /**
     * @brief
     * @return
     */
    EWSUSY* getMyEWSUSY() const
    {
        return myEWSUSY;
    }

    /**
     * @brief
     * @return
     */
    FeynHiggsWrapper* getMyFH() const
    {
        return myFH;
    }


    ///////////////////////////////////////////////////////////////////////////
    // functions for the input parameters of SUSY model

    /**
     * @brief Gets the bino mass.
     * @return The bino mass.
     */
    complex getM1() const
    {
        return m1;
    }

    /**
     * @brief Gets the wino mass.
     * @return The wino mass.
     */
    complex getM2() const
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
    complex getMuH() const
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
    double getQ() const
    {
        return Q;
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
    complex getSaeff() const
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
    vector<double> getMch() const
    {
        return mch;
    }

    /**
     * @brief Gets the rotation matrix for negative charginos.
     * @return The rotation matrix for negative charginos.
     */
    matrix<complex> getU() const
    {
        return U;
    }

    /**
     * @brief Gets the rotation matrix for positive charginos.
     * @return The rotation matrix for positive charginos.
     */
    matrix<complex> getV() const
    {
        return V;
    }


    ///////////////////////////////////////////////////////////////////////////
    // functions for the parameters in the neutralino sector.
    
    /**
     * @brief Gets the neutralino masses.
     * @return The neutralino masses. 
     */
    vector<double> getMneu() const
    {
        return mneu;
    }

    /**
     * @brief Gets the rotation matrix for neutralinos.
     * @return The rotation matrix for neutralinos.
     */
    matrix<complex> getN() const
    {
        return N;
    }


    ///////////////////////////////////////////////////////////////////////////
    // functions for the parameters in the squark sector.

    matrix<complex> getMsQhat2() const
    {
        return msQhat2;
    }

    matrix<complex> getMsUhat2() const
    {
        return msUhat2;
    }

    matrix<complex> getMsDhat2() const
    {
        return msDhat2;
    }

    /**
     * @brief Gets the trilinear-coupling matrix for up-type squarks.
     * @return The trilinear-coupling matrix for up-type squarks.
     */
    matrix<complex> getTUhat() const
    {
        return TUhat;
    }

    /**
     * @brief Gets the trilinear-coupling matrix for down-type squarks.
     * @return The trilinear-coupling matrix for down-type squarks.
     */
    matrix<complex> getTDhat() const
    {
        return TDhat;
    }

    /**
     * @brief Gets the up-type squark mass squared.
     * @return The up-type squark mass squared.
     */
    vector<double> getMsu2() const
    {
        return m_su2;
    }

    /**
     * @brief Gets the down-type squark mass squared.
     * @return The down-type squark mass squared. 
     */
    vector<double> getMsd2() const
    {
        return m_sd2;
    }

    /**
     * @brief Gets the rotation matrix for up-type squarks.
     * @return The rotation matrix for up-type squarks.
     */
    matrix<complex> getRu() const
    {
        return Ru;
    }

    /**
     * @brief Gets the rotation matrix for down-type squarks.
     * @return The rotation matrix for down-type squarks.
     */
    matrix<complex> getRd() const
    {
        return Rd;
    }


    ///////////////////////////////////////////////////////////////////////////
    // functions for the parameters in the slepton sector.

    matrix<complex> getMsLhat2() const
    {
        return msLhat2;
    }

    matrix<complex> getMsNhat2() const
    {
        return msNhat2;
    }

    matrix<complex> getMsEhat2() const
    {
        return msEhat2;
    }

    /**
     * @brief Gets the trilinear-coupling matrix for sneutrinos.
     * @return The trilinear-coupling matrix for sneutrinos.
     */
    matrix<complex> getTNhat() const
    {
        return TNhat;
    }

    /**
     * @brief Gets the trilinear-coupling matrix for charged sleptons.
     * @return The trilinear-coupling matrix for charged sleptons.
     */
    matrix<complex> getTEhat() const
    {
        return TEhat;
    }

    /**
     * @brief Gets the sneutrino mass squared.
     * @return The sneutrino mass squared.
     */
    vector<double> getMsn2() const
    {
        return m_sn2;
    }

    /**
     * @brief Gets the charged slepton mass squared.
     * @return The charged slepton mass squared.
     */
    vector<double> getMse2() const
    {
        return m_se2;
    }

    /**
     * @brief Gets the rotation matrix for sneutrinos.
     * @return The rotation matrix for sneutrinos.
     */
    matrix<complex> getRn() const
    {
        return Rn;
    }

    /**
     * @brief Gets the rotation matrix for charged sleptons.
     * @return The rotation matrix for charged sleptons.
     */
    matrix<complex> getRl() const
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


    ///////////////////////////////////////////////////////////////////////////
    // EW precision observables

    /**
     * @brief The W boson mass.
     * @return @f$M_W@f$.
     */
    double Mw() const;

    /**
     * @brief The W boson mass in the @f$\Delta\rho@f$ approximation.
     * @return @f$M_W@f$ in the @f$\Delta\rho@f$ approximation.
     */
    double Mw_dRho() const;

    /**
     * @brief
     * @return
     */
    double cW2() const;
    
    /**
     * @brief
     * @return
     */
    double sW2() const;


    ///////////////////////////////////////////////////////////////////////////

protected:
    virtual void SetParameter(const std::string, const double&);
    virtual void SetTanb(const double tanb);
    virtual void SetYukawas();
    virtual void SetSoftTerms();

    FeynHiggsWrapper* myFH;
    EWSUSY* myEWSUSY;

    // model parameters at scale Q
    complex m1, m2, muH;
    double m3, mHptree, tanb, Q;

    // sinb and cosb computed with setTanb()
    double sinb, cosb;

    // soft-breaking parameters associated with squark and slepton mass terms and
    // trilinear couplings in the SCKM basis, which will be set with SetSoftTerms()
    matrix<complex> msQhat2, msUhat2, msDhat2, msLhat2, msNhat2, msEhat2;
    matrix<complex> TUhat, TDhat, TNhat, TEhat;

    // soft-breaking parameters computed with FeynHiggs
    double mHp, mh[4];
    complex saeff;
    vector<double> mch, mneu;
    vector<double> m_su2, m_sd2, m_sn2, m_se2;

    // rotation matrices
    matrix<complex> U, V, N, Ru, Rd, Rn, Rl;

    // quark and lepton masses at scale Q, computed in setYukawas()
    double mu_Q[3], md_Q[3], me_Q[3], mn_Q[3];

    
    ///////////////////////////////////////////////////////////////////////////
private:    
    bool flag_h, flag_g, flag_ch, flag_ne;
    SUSYMatching* mySUSYMatching;

};

/** 
 * @}
 */

#endif	/* SUSY_H */

