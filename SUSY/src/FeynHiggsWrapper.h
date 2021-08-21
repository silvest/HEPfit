/*
 * Copyright (C) 2013 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */


#ifndef FEYNHIGGSWRAPPER_H
#define	FEYNHIGGSWRAPPER_H

#include "SUSY.h"
#if FEYNHIGGS
#include <CFeynHiggs.h>
#include <gslpp.h>

/**
 * @class FeynHiggsWrapper
 * @ingroup SUSY
 * @brief A wrapper class for FeynHiggs library.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details
 */
class FeynHiggsWrapper {
public:

    /**
     * @brief A FeynHiggsWrapper constructor.
     * @param[in] SUSY_in An object of SUSY class.
     */
    FeynHiggsWrapper(SUSY& SUSY_in);

    
    ///////////////////////////////////////////////////////////////////////////

    /**
     * @brief Sets the FeynHiggs parameters. 
     * @return Zero if successful.
     */
    bool SetFeynHiggsPars();

    /**
     * @brief Computes the Higgs masses and mixings with FeynHiggs.
     * @return Zero if successful. 
     */
    bool CalcHiggsSpectrum();

    /**
     * @brief Computes the mass spectra of squarks, sleptons, charginos and neutralinos with FeynHiggs.
     * @return Zero if successful. 
     */
    bool CalcSpectrum();

    /**
     * @brief Sets the FeynHiggs input parameters with an SLHA file, used for
     * test and debug. 
     */
    void SetFeynHiggsParsSLHA(const char *filename) const;

    /**
     * @brief Writes FeynHiggs outputs to an SLHA file for test and debug.
     * @param[in] filename 
     */
    void OutputSLHA(const char *filename) const;


    ///////////////////////////////////////////////////////////////////////////

    /**
     * @brief The W-boson mass used as an input to FHSetSMPara(). 
     * @return The W-boson mass used as an input to FHSetSMPara(). 
     */
    double getMw_FHinput() const
    {
        return Mw_FHinput;
    }

    /**
     * @brief Gets the correction to the bottom Yukawa coupling @f$\Delta_b@f$. 
     * @return The correction to the bottom Yukawa coupling @f$\Delta_b@f$
     */
    gslpp::complex getFHDeltab() const
    {
        return FHDeltab;
    }

    /**
     * @brief Gets the gluino mass.
     * @return The gluino mass. 
     */
    double getMGl() const
    {
        return FHMGl;
    }

    /**
     * @brief Gets the tree-level Higgs masses.
     * @return The tree-level Higgs masses. 
     */
    double getFHMHtree(const int i) const
    {
        return FHMHtree[i];
    }

    /**
     * @brief Gets the tree-level @f$\sin(\alpha)@f$.
     * @return The tree-level @f$\sin(\alpha)@f$. 
     */
    double getFHSAtree() const
    {
        return FHSAtree;
    }
    
    /**
     * @brief Gets @f$\alpha_s(m_{top})@f$.
     * @return @f$\alpha_s(m_{top})@f$. 
     */
    double getFHAlfasMT() const
    {
        return AlfasMT;
    }
    

    /**
     * @brief Gets the muon anomalous magnetic moment.
     * @return The muon @f$g-2@f$.
     */
    double getFHgm2() 
    {
        if (computeConstraints) CalcConstraints();
        return FHgm2;
    }

    /**
     * @brief Gets @f$\Delta\r@f$.
     * @return @f$\Delta\r@f$.
     */
    double getFHdeltar()
    {
        if (computeConstraints) CalcConstraints();
        return FHdeltar;
    }

    /**
     * @brief Gets @f$\Delta\rho@f$.
     * @return @f$\Delta\rho@f$.
     */
    double getFHdeltarho()
    {
        if (computeConstraints) CalcConstraints();
        return FHdeltarho;
    }

    /**
     * @brief Gets the W-boson mass in the MSSM.
     * @return The W-boson mass in the MSSM.
     */
    double getFHMWMSSM()
    {
        if (computeConstraints) CalcConstraints();
        return FHMWMSSM;
    }

    /**
     * @brief Gets the W-boson mass in the SM.
     * @return The W-boson mass in the SM.
     */
    double getFHMWSM()
    {
        if (computeConstraints) CalcConstraints();
        return FHMWSM;
    }

    /**
     * @brief Gets the effective weak mixing angle in the MSSM.
     * @return The effective weak mixing angle in the MSSM.
     */
    double getFHSW2MSSM()
    {
        if (computeConstraints) CalcConstraints();
        return FHSW2MSSM;
    }

    /**
     * @brief Gets the effective weak mixing angle in the SM.
     * @return The effective weak mixing angle in the SM.
     */
    double getFHSW2SM()
    {
        if (computeConstraints) CalcConstraints();
        return FHSW2SM;
    }

    /**
     * @brief Gets the electron EDM.
     * @return The electron EDM.
     */
    double getFHedmeTh()
    {
        if (computeConstraints) CalcConstraints();
        return FHedmeTh;
    }

    /**
     * @brief Gets the neutron EDM. 
     * @return The neutron EDM.
     */
    double getFHedmn()
    {
        if (computeConstraints) CalcConstraints();
        return FHedmn;
    }

    /**
     * @brief Gets the mercury EDM.
     * @return The mercury EDM.
     */
    double getFHedmHg()
    {
        if (computeConstraints) CalcConstraints();
        return FHedmHg;
    }
    
    /**
     * @brief Gets the branching ratio for @f$B\to X_s\gamma@f$ in the MSSM.
     * @return @f$Br(B\to X_s\gamma)@f$ in the MSSM.
     */
    double getFHbsgMSSM()
    {
        if (computeFlavour) CalcFlavour();
        return FHbsgMSSM;
    }

    /**
     * @brief Gets the branching ratio for @f$B\to X_s\gamma@f$ in the SM.
     * @return @f$Br(B\to X_s\gamma)@f$ in the SM.
     */
    double getFHbsgSM()
    {
        if (computeFlavour) CalcFlavour();
        return FHbsgSM;
    }

    /**
     * @brief Gets @f$\Delta M_s@f$ in the MSSM.
     * @return @f$\Delta M_s@f$ in the MSSM.
     */
    double getFHdeltaMsMSSM()
    {
        if (computeFlavour) CalcFlavour();
        return FHdeltaMsMSSM;
    }

    /**
     * @brief Gets @f$\Delta M_s@f$ in the SM. 
     * @return @f$\Delta M_s@f$ in the SM.
     */
    double getFHdeltaMsSM()
    {
        if (computeFlavour) CalcFlavour();
        return FHdeltaMsSM;
    }

    /**
     * @brief Gets the branching ratio for @f$B_s\to\mu^+\mu^-@f$ in the MSSM.
     * @return @f$Br(B_s\to\mu^+\mu^-)@f$ in the MSSM.
     */
    double getFHbsmumuMSSM()
    {
        if (computeFlavour) CalcFlavour();
        return FHbsmumuMSSM;
    }

    /**
     * @brief Gets the branching ratio for @f$B_s\to\mu^+\mu^-@f$ in the SM.
     * @return @f$Br(B_s\to\mu^+\mu^-)@f$ in the SM.
     */
    double getFHbsmumuSM()
    {
        if (computeFlavour) CalcFlavour();
        return FHbsmumuSM;
    }

    /**
     * Sort sfermion masses in increasing order.
     * @param[in,out] m_sf2 A vector of sfermion mass squared.
     * @param[in,out] Rf The corresponding rotation matrix.
     */
    void SortSfermionMasses(gslpp::vector<double>& m_sf2, gslpp::matrix<gslpp::complex>& Rf) const;


    ///////////////////////////////////////////////////////////////////////////
private:
    SUSY& mySUSY;

    double Mw_FHinput; /* The W-boson mass used as an input to FHSetSMPara() */

    bool NMFVu, NMFVd, NMFVe; /* true if off-diagonal entries exist in the sfermion mass matrix */
    //bool NMFVnu;

    // see CalcSpectrum()
    gslpp::complex FHDeltab; /* the correction to the bottom Yukawa coupling */

    // see FHGetPara()
    double FHMGl, FHMHtree[4], FHSAtree, AlfasMT; 

    // see CalcConstraints()
    double FHgm2; /* muon g-2 */
    double FHdeltar;
    double FHdeltarho; /* Delta rho */
    double FHMWMSSM, FHMWSM; /* the W mass */
    double FHSW2MSSM, FHSW2SM; /* the effective weak mixing angle */
    double FHedmeTh, FHedmn, FHedmHg; /* the EDMs of electron, neutron and mercury */

    // see CalcFlavour()
    double FHbsgMSSM, FHbsgSM; /* B -> X_s gamma */
    double FHdeltaMsMSSM, FHdeltaMsSM; /* Delta Ms */
    double FHbsmumuMSSM, FHbsmumuSM; /* Bs -> mu+ mu- */

    bool computeHiggsCouplings, computeHiggsProd, computeConstraints, computeFlavour;

    /**
     * @brief Computes the Higgs couplings, decay widths and branching ratios with FeynHiggs.
     * @return Zero if successful.
     */
    bool CalcHiggsCouplings();

    /**
     * @brief Computes Higgs production cross-sections with FeynHiggs.
     * @return Zero if successful.
     */
    bool CalcHiggsProd(const double&);

    /**
     * @brief Calculates electroweak precision observables.
     * @return Zero if successful.
     */
    bool CalcConstraints();

    /**
     * @brief Calculates flavour observables.
     * @return Zero if successful.
     */
    bool CalcFlavour();

};

#endif	/* FEYNHIGGSWRAPPER_H */
#endif
