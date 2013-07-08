/* 
 * Copyright (C) 2012-2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef SUSY_H
#define	SUSY_H

#include <StandardModel.h>
#include <EWSM.h>
#include "SUSYMatching.h"

using namespace gslpp;

class FeynHiggs;
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
     * @brief A friend classes of SUSY class.
     */
    friend class SUSYSpectrum;
    friend class FeynHiggs;

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
    // functions for the input parameters of SUSY model

    /**
     * @brief Gets the bino mass.
     * @return The bino mass.
     */
    double getM1() const
    {
        return m3;
    }

    /**
     * @brief Gets the wino mass.
     * @return The wino mass.
     */
    double getM2() const
    {
        return m3;
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
     * @brief Gets the light Higgs mass.
     * @return The light Higgs mass.
     */
    virtual double getMHl() const
    {
        return mh[0];
    }

    /**
     * @brief Gets the Heavy Higgs mass.
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
     * @brief Gets the chargino mass.
     * @return The chargino mass.
     */
    vector<double> getMch() const
    {
        return Mch;
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
     * @brief Gets the neutralino mass.
     * @return The neutralino mass
     */
    vector<double> getMneu() const
    {
        return Mneu;
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

    /**
     * @brief Gets the trilinear-coupling matrix for up-type squarks.
     * @return The trilinear-coupling matrix for up-type squarks.
     */
    matrix<complex> getTU() const
    {
        return TU;
    }

    /**
     * @brief Gets the trilinear-coupling matrix for down-type squarks.
     * @return The trilinear-coupling matrix for down-type squarks.
     */
    matrix<complex> getTD() const
    {
        return TD;
    }

    /**
     * @brief Gets the up-type squark mass squared.
     * @return The up-type squark mass squared.
     */
    vector<double> getMsu2() const
    {
        return Msu2;
    }

    /**
     * @brief Gets the down-type squark mass squared.
     * @return The down-type squark mass squared. 
     */
    vector<double> getMsd2() const
    {
        return Msd2;
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

    /**
     * @brief Gets the trilinear-coupling matrix for sneutrinos.
     * @return The trilinear-coupling matrix for sneutrinos.
     */
    matrix<complex> getTN() const
    {
        return TN;
    }

    /**
     * @brief Gets the trilinear-coupling matrix for charged sleptons.
     * @return The trilinear-coupling matrix for charged sleptons.
     */
    matrix<complex> getTE() const
    {
        return TE;
    }

    /**
     * @brief Gets the sneutrino mass squared.
     * @return The sneutrino mass squared.
     */
    vector<double> getMsn2() const
    {
        return Msn2;
    }

    /**
     * @brief Gets the charged slepton mass squared.
     * @return The charged slepton mass squared.
     */
    vector<double> getMsl2() const
    {
        return Msl2;
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
    // EW precision observables

    /**
     * @brief Computes the W boson mass. 
     * @return The W boson mass. 
     */
    double Mw() const;

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

    /**
     * @brief
     * @param[in] INDF fermion index [0-9] (see EWphysics::flavour_st_to_int())
     * @return The ratio of the effective vector coupling constants @f$g_Z^f=g_V^f/g_A^f@f$ for INDF.
     */
    complex gZf(const int INDF) const; // gZf = gVf/gAf

    /**
     * @brief
     * @param[in] INDF fermion index [0-9] (see EWphysics::flavour_st_to_int())
     * @return The weak form factor for INDF.
     */
    complex rhoZf(const int INDF) const;

    /**
     * @brief
     * @return The radiative-correction factor @f$\Delta r@f$. 
     */
    double Delta_r() const;


    ///////////////////////////////////////////////////////////////////////////

protected:
    virtual void SetParameter(const std::string, const double&);
    virtual void SetSoftTerms();
    void setTanb(const double tanb);
    void setYukawas();
    FeynHiggs* myFH;
    EWSUSY* myEWSUSY;

    // model parameters at scale Q
    complex m1, m2, muH;
    double m3, mHptree, tanb, Q;

    // sinb and cosb computed with setTanb()
    double sinb, cosb;

    // squark and slepton mass matrices and trilinear couplings in the SCKM basis,
    // which are set with SetSoftTerms()
    matrix<complex> MsQ2, MsU2, MsD2, MsL2, MsN2, MsE2, TU, TD, TN, TE;

    // other soft-breaking parameters computed with FeynHiggs
    double mHp, mh[4];
    complex saeff;
    vector<double> Mch, Mneu;
    vector<double> Msu2, Msd2, Msn2, Msl2;

    // rotation matrices
    matrix<complex> U, V, N, Ru, Rd, Rn, Rl;

    // quark masses at scale Q, computed in setYukawas()
    double mu_Q[3], md_Q[3];

    
    ///////////////////////////////////////////////////////////////////////////
private:    
    bool flag_h, flag_g, flag_ch, flag_ne;
    SUSYMatching* mySUSYMatching;
    int updateFlag;

};

/** 
 * @}
 */

#endif	/* SUSY_H */

