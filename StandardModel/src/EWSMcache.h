/* 
 * Copyright (C) 2012-2014 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef EWSMCACHE_H
#define	EWSMCACHE_H

#include <cmath>
#include <PVfunctions.h>
#include <Polylogarithms.h>
#include <ClausenFunctions.h>
#include "StandardModel.h"

using namespace gslpp;

/**
 * @class EWSMcache
 * @ingroup StandardModel
 * @brief A class for cache variables used in computing radiative corrections
 * to the %EW precision observables.
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class EWSMcache {

public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    EWSMcache(const StandardModel& SM_i);

    
    ////////////////////////////////////////////////////////////////////////

    void setBDebug(bool bDebug)
    {
        this->bDebug = bDebug;
    }

    
    ////////////////////////////////////////////////////////////////////////     
    
    /**
     * @return a reference to the StandardModel object
     */
    const StandardModel& getSM() const 
    {
        return SM;
    }
    
    /**
     * @return the object of PVfunctions class
     */
    const PVfunctions getPV() const 
    {
        return PV;
    }

    /**
     * @return the object of Polylogarithms class
     */
    const Polylogarithms getPolyLog() const 
    {
        return PolyLog;
    }

    /**
     * @return the object of ClausenFunctions class
     */
    const ClausenFunctions getClausen() const 
    {
        return Clausen;
    }

    
    ////////////////////////////////////////////////////////////////////////         
    // get functions for constants
    
    /**
     * @return the zeta function zeta(2)
     */
    double GetZeta2() const 
    {
        return zeta2;
    }

    /**
     * @return the zeta function zeta(3)
     */
    double GetZeta3() const 
    {
        return zeta3;
    }
    
    /**
     * @return the zeta function zeta(4)
     */
    double GetZeta4() const 
    {
        return zeta4;
    }        
    
    /**
     * @return the zeta function zeta(5)
     */
    double GetZeta5() const 
    {
        return zeta5;
    }
    
    double GetS2() const 
    {
        return S2;
    }
 
    double GetD3() const 
    {
        return D3;
    }   
    
    double GetB4() const 
    {
        return B4;
    }    

    /**
     * @return the logarithmic function log(2)
     */
    double GetLog2() const 
    {
        return log2;
    }    

    
    //////////////////////////////////////////////////////////////////////// 

    /**
     * @param[in] l name of lepton
     * @return mass of lepton
     */
    double ml(const StandardModel::lepton l) const;

    /**
     * @param[in] l name of lepton
     * @return mass squared of lepton
     */
    double ml2(const StandardModel::lepton l) const
    {
        return ( ml(l)*ml(l) );
    }

    /**
     * @param[in] q name of quark
     * @param[in] mu renormalization scale
     * @param[in] order (=LO, NLO, NNLO, FULLNLO, FULLNNLO[defalut])
     * @return the MSbar mass of u, d, s, c, b or the pole mass of t
     */
    double mq(const StandardModel::quark q, const double mu, 
              const orders order=FULLNNLO) const;

    /**
     * @param[in] q name of quark
     * @param[in] mu renormalization scale
     * @param[in] order (=LO, NLO, NNLO, FULLNLO, FULLNNLO[defalut])
     * @return the MSbar mass squared of u, d, s, c, b or the pole mass squared of t
     */
    double mq2(const StandardModel::quark q, const double mu,
               const orders order=FULLNNLO) const
    {
        return ( mq(q, mu, order)*mq(q, mu, order) );
    }
    
    /**
     * @return The top-quark mass.
     */
    double Mt() const 
    {
        return SM.getMtpole();
    }
    
    /**
     * @return The top-quark mass squared.
     */
    double Mt2() const 
    {
        return ( SM.getMtpole()*SM.getMtpole() ); 
    }

    /**
     * @return alpha_s(M_z^2)
     */
    double alsMz() const 
    {
        return SM.getAlsMz();
    }
    
    /**
     * @return G_F
     */
    double GF() const
    {
        return SM.getGF();
    }
    
    /**
     * @return electromagnetic coupling at q^2=0
     */
    double ale() const
    {
        return SM.getAle();
    }

    /**
     * @return five-flavour hadronic correction to alpha at Mz^2
     */
    double dAle5Mz() const 
    {
        return SM.getDAle5Mz();
    }

    /**
     * @return The Z-boson mass.
     */
    double Mz() const 
    {
        return SM.getMz();
    }

    /**
     * @return The Z-boson mass squared.
     */
    double Mz2() const
    {
        return ( SM.getMz()*SM.getMz() );
    }

    /**
     * @return The higgs mass.
     */
    double mh() const 
    {
        return SM.getMHl();
    }
    
    /**
     * @return The higgs mass squared.
     */
    double mh2() const
    {
        return ( SM.getMHl()*SM.getMHl() );
    }

    /**
     * @param[in] Mw_i The W-boson mass. 
     * @return The W-boson mass. 
     */
    double Mw(const double Mw_i) const
    {
        return Mw_i;
    }
    
    /**
     * @param[in] Mw_i The W-boson mass.
     * @return The W-boson mass squared.
     */
    double Mw2(const double Mw_i) const
    {
        return ( Mw_i*Mw_i );
    }

    /**
     * @param[in] Mw_i the W-boson mass
     * @return c_W^2
     */
    double cW2(const double Mw_i) const 
    {
        return ( Mw(Mw_i)*Mw(Mw_i)/Mz()/Mz() );
    }   
    
    /**
     * @param[in] Mw_i the W-boson mass
     * @return s_W^2
     */
    double sW2(const double Mw_i) const 
    {
        return ( 1.0 - cW2(Mw_i) );
    } 
    
    /**
     * @param[in] l name of lepton
     * @return electric charge of a lepton "l"
     */
    double Ql(const StandardModel::lepton l) const
    {
        return SM.getLeptons(l).getCharge();
    }    

    /**
     * @param[in] q name of quark
     * @return electric charge of a quark "q"
     */
    double Qq(const StandardModel::quark q) const
    {
        return SM.getQuarks(q).getCharge();
    }

    /**
     * @param[in] l name of lepton
     * @param[in] Mw_i the W-boson mass
     * @return the tree-level vector coupling for Z->l lbar
     * @attention depends on sW2
     */
    double vl(const StandardModel::lepton l, const double Mw_i) const 
    {
        return ( al(l) - 2.0*Ql(l)*sW2(Mw_i) );
    }

    /**
     * @param[in] q name of quark
     * @param[in] Mw_i the W-boson mass
     * @return the tree-level vector coupling for Z->q qbar
     * @attention depends on sW2
     */
    double vq(const StandardModel::quark q, const double Mw_i) const
    {
        return ( aq(q) - 2.0*Qq(q)*sW2(Mw_i) );
    }
    
    /**
     * @param[in] l name of lepton
     * @return the tree-level axial-vector coupling for Z->l lbar
     */
    double al(const StandardModel::lepton l) const
    {
        return ( SM.getLeptons(l).getIsospin() );
    }

    /**
     * @param[in] q name of quark
     * @return the tree-level axial-vector coupling for Z->q qbar
     */
    double aq(const StandardModel::quark q) const
    {
        return ( SM.getQuarks(q).getIsospin() );
    }
    
    /**
     * @param[in] l name of lepton
     * @param[in] Mw_i the W-boson mass
     * @return |v_l+a_l| 
     * @attention depends on sW2
     */
    double sigmal(const StandardModel::lepton l, const double Mw_i) const
    {
        return ( 1.0 - 2.0*fabs(Ql(l))*sW2(Mw_i) );
    }

    /**
     * @param[in] q name of quark
     * @param[in] Mw_i the W-boson mass
     * @return |v_q+a_q| 
     * @attention depends on sW2
     */
    double sigmaq(const StandardModel::quark q, const double Mw_i) const 
    {
        return ( 1.0 - 2.0*fabs(Qq(q))*sW2(Mw_i) );
    }  
    
    /**
     * @param[in] l name of lepton
     * @param[in] Mw_i the W-boson mass
     * @return v_l-a_l 
     * @attention depends on sW2
     */    
    double deltal(const StandardModel::lepton l, const double Mw_i) const 
    {
        return ( - 2.0*Ql(l)*sW2(Mw_i) );   
    } 

    /**
     * @param[in] q name of quark
     * @param[in] Mw_i the W-boson mass
     * @return v_q-a_q 
     * @attention depends on sW2
     */    
    double deltaq(const StandardModel::quark q, const double Mw_i) const 
    {
        return ( - 2.0*Qq(q)*sW2(Mw_i) );   
    }  
    
    /**
     * @param[in] Mw_i the W-boson mass
     * @return the conversion factor from alpha to GF
     */
    double f_AlphaToGF(const double Mw_i) const 
    {
        return ( sqrt(2.0)*GF()*pow(Mz(),2.0)*sW2(Mw_i)*cW2(Mw_i)/M_PI/ale() );
    }
   
    /**
     * @return X_t with G_F
     */
    double Xt_GF() const 
    {
        return ( GF()*Mt()*Mt()/8.0/sqrt(2.0)/M_PI/M_PI );
    }
    
    /**
     * @param[in] Mw_i the W-boson mass
     * @return X_t with alpha(0)
     */
    double Xt_alpha(const double Mw_i) const 
    {
        return ( Xt_GF()/f_AlphaToGF(Mw_i) );
    }
    
    /**
     * @return alpha_s(M_t^2)
     */
    double alsMt() const 
    {
        if (bDebug)
            return ( 0.1074432788 );// for debug
        else
            return ( SM.Als(Mt(),FULLNNLO) );
    }

    /**
     * @brief for two-loop QCD corrections
     * @return 
     */
    double Phi_QCD2() const 
    {
        double r_QCD2 = Mz()*Mz()/4.0/Mt()/Mt();
        return ( asin(sqrt(r_QCD2)) );
    }

    /**
     * @brief for two-loop QCD corrections
     * @return 
     */
    double gamma_QCD2() const 
    {
        double r_QCD2 = Mz()*Mz()/4.0/Mt()/Mt();
        return ( log(2.0*sqrt(r_QCD2)) );
    }

    /**
     * @brief for two-loop QCD corrections
     * @return 
     */
    double h_QCD2() const 
    {
        double r_QCD2 = Mz()*Mz()/4.0/Mt()/Mt();
        return ( log(2.0*sqrt(1.0-r_QCD2)) );
    }
    
    /**
     * @brief for two-loop QCD corrections
     * @return 
     */
    double logV1primeAndA1prime() const 
    {
        gsl_complex OneMinusE2Iphi = gsl_complex_rect(1.0-cos(2.0*Phi_QCD2()), 
                                     -sin(2.0*Phi_QCD2()));
        gsl_complex OneMinusE4Iphi = gsl_complex_rect(1.0-cos(4.0*Phi_QCD2()), 
                                     -sin(4.0*Phi_QCD2()));
        return ( GSL_REAL(gsl_complex_log(OneMinusE2Iphi))
                - 2.0*GSL_REAL(gsl_complex_log(OneMinusE4Iphi)) );
    }
    
    /**
     * @brief for two-loop QCD corrections
     * @return 
     */
    double Cl2_2Phi() const 
    {
        double Phi= asin(Mz()/2.0/Mt());
        return ( Clausen.Cl2(2.0*Phi) );
    }

    /**
     * @brief for two-loop QCD corrections
     * @return 
     */
    double Cl2_4Phi() const 
    {
        double Phi= asin(Mz()/2.0/Mt());
        return ( Clausen.Cl2(4.0*Phi) );
    }

    /**
     * @brief for two-loop QCD corrections
     * @return 
     */
    double Cl3_2Phi() const 
    {
        double Phi= asin(Mz()/2.0/Mt());
        return ( Clausen.Cl3(2.0*Phi) );
    }

    /**
     * @brief for two-loop QCD corrections
     * @return 
     */
    double Cl3_4Phi() const 
    {
        double Phi= asin(Mz()/2.0/Mt());
        return ( Clausen.Cl3(4.0*Phi) );
    }

    
    ////////////////////////////////////////////////////////////////////////  
    // Functions for the caches
    
    double logMZtoME() const;
    double logMZtoMMU() const;
    double logMZtoMTAU() const;
    double logMZtoMTOP() const;
    double logMTOPtoMH() const;
    double log_cW2(const double Mw_i) const;
    
    double Li2_MW2toMTOP2(const double Mw_i) const;
    double Li3_MW2toMTOP2(const double Mw_i) const;
    double Li3_for_F1(const double Mw_i) const;
    
    double A0_Mz2_Mw2(const double Mw_i) const;
    double A0_Mz2_mh2() const;
    double A0_Mw2_Mz2(const double Mw_i) const;
    double A0_Mw2_mh2(const double Mw_i) const;
    double A0_Mz2_Mz2() const;
    double A0_Mw2_Mw2(const double Mw_i) const;
    
    complex B0_Mz2_Mw2_mh2_Mw2(const double Mw_i) const;
    complex B0_Mz2_0_mh2_Mw2(const double Mw_i) const;
    complex B0_Mw2_Mz2_Mt2_Mt2(const double Mw_i) const;
    complex B0_Mz2_Mz2_Mw2_Mw2(const double Mw_i) const;
    complex B0_Mz2_Mz2_mh2_Mz2() const;
    complex B0_Mz2_Mw2_Mz2_Mw2(const double Mw_i) const;
    complex B0_Mz2_Mw2_0_Mw2(const double Mw_i) const;
    complex B0_Mz2_0_Mz2_Mw2(const double Mw_i) const;
    complex B0_Mz2_0_0_Mw2(const double Mw_i) const;
    complex B0_Mw2_Mz2_Mw2_Mw2(const double Mw_i) const;
    complex B0_Mw2_Mw2_Mz2_Mw2(const double Mw_i) const;
    complex B0_Mw2_Mw2_mh2_Mw2(const double Mw_i) const;
    complex B0_Mw2_Mw2_0_Mw2(const double Mw_i) const;
    complex B0_Mz2_Mz2_ml2_ml2(const StandardModel::lepton l) const;
    complex B0_Mz2_Mz2_mq2_mq2(const StandardModel::quark q) const;
    
    complex B0p_Mz2_0_mh2_Mw2(const double Mw_i) const;
    complex B0p_Mz2_Mz2_mh2_Mz2() const;
    complex B0p_Mz2_0_Mz2_Mw2(const double Mw_i) const;
    complex B0p_Mz2_Mz2_Mw2_Mw2(const double Mw_i) const;
    complex B0p_Mw2_Mw2_Mz2_Mw2(const double Mw_i) const;
    complex B0p_Mw2_Mw2_mh2_Mw2(const double Mw_i) const;
    complex B0p_Mw2_Mw2_0_Mw2(const double Mw_i) const;
    complex B0p_Mz2_Mz2_ml2_ml2(const StandardModel::lepton l) const;
    complex B0p_Mz2_Mz2_mq2_mq2(const StandardModel::quark q) const;
    
    complex B1_Mz2_0_ml2_mlprime2(const int gen) const;
    complex B1_Mz2_0_mq2_mqprime2(const int gen) const;
    complex B1_Mz2_0_mlprime2_ml2(const int gen) const;
    complex B1_Mz2_0_mqprime2_mq2(const int gen) const;
    complex B1_Mz2_Mw2_ml2_mlprime2(const int gen, const double Mw_i) const;
    complex B1_Mz2_Mw2_mq2_mqprime2(const int gen, const double Mw_i) const;
    complex B1_Mz2_Mw2_mlprime2_ml2(const int gen, const double Mw_i) const;
    complex B1_Mz2_Mw2_mqprime2_mq2(const int gen, const double Mw_i) const;
    
    complex B1p_Mw2_Mw2_ml2_mlprime2(const int gen, const double Mw_i) const;
    complex B1p_Mw2_Mw2_mq2_mqprime2(const int gen, const double Mw_i) const;
    complex B1p_Mw2_Mw2_mlprime2_ml2(const int gen, const double Mw_i) const;
    complex B1p_Mw2_Mw2_mqprime2_mq2(const int gen, const double Mw_i) const;
    
    complex Bf_Mz2_Mz2_ml2_ml2(const StandardModel::lepton l) const;
    complex Bf_Mz2_Mz2_mq2_mq2(const StandardModel::quark q) const;
    complex Bf_Mz2_0_ml2_ml2(const StandardModel::lepton l) const;
    complex Bf_Mz2_0_mq2_mq2(const StandardModel::quark q) const;
    complex Bf_Mz2_Mw2_mlprime2_ml2(const int gen, const double Mw_i) const;
    complex Bf_Mz2_Mw2_mqprime2_mq2(const int gen, const double Mw_i) const;
    complex Bf_Mz2_0_mlprime2_ml2(const int gen) const;
    complex Bf_Mz2_0_mqprime2_mq2(const int gen) const;
    complex Bf_Mw2_Mw2_mlprime2_ml2(const int gen, const double Mw_i) const;
    complex Bf_Mw2_Mw2_mqprime2_mq2(const int gen, const double Mw_i) const;
    
    complex Bfp_Mz2_Mz2_ml2_ml2(const StandardModel::lepton l) const;
    complex Bfp_Mz2_Mz2_mq2_mq2(const StandardModel::quark q) const;
    complex Bfp_Mw2_Mw2_mlprime2_ml2(const int gen, const double Mw_i) const;
    complex Bfp_Mw2_Mw2_mqprime2_mq2(const int gen, const double Mw_i) const;
    
    complex C0_Mz2_Mw2_Mt2_Mw2(const double Mw_i) const;
    complex C0_Mz2_Mt2_Mw2_Mt2(const double Mw_i) const;
    complex C0_Mz2_0_Mw2_0(const double Mw_i) const;
    complex C0_Mz2_Mw2_0_Mw2(const double Mw_i) const;
    complex C0_Mw2_Mw2_0_Mz2(const double Mw_i) const;
    complex C0_Mw2_0_Mz2_0(const double Mw_i) const;
    complex C0_Mz2_0_Mz2_0() const;
    
    
    //////////////////////////////////////////////////////////////////////// 

private:
    bool bDebug;// for debug
    bool bUseCacheEWSMcache;// true for caching
    
    const StandardModel& SM;///< A reference to an object of type StandardModel.
    const PVfunctions PV;///< An object of type PVfunctions.
    const ClausenFunctions Clausen;///< An object of type ClausenFunctions.
    const Polylogarithms PolyLog;///< An object of type Polylogarithms.
    
    // Constants 
    double zeta2, zeta3, zeta4, zeta5;
    double S2, D3, B4;    
    double log2;

    
    ////////////////////////////////////////////////////////////////////////     
    // Caches for double variables
    // The last element is the cache for a double variable. 
    
    /* Logarithms */
    mutable double logMZtoME_cache[3];
    mutable double logMZtoMMU_cache[3];
    mutable double logMZtoMTAU_cache[3];
    mutable double logMZtoMTOP_cache[3]; 
    mutable double logMTOPtoMH_cache[3]; 
    mutable double log_cW2_cache[3]; 
    
    /* Dilogarithm and Trilogarithm for two-loop QCD corrections */
    mutable double Li2_MW2toMTOP2_cache[3]; 
    mutable double Li3_MW2toMTOP2_cache[3]; 
    mutable double Li3_for_F1_cache[3]; 
    
    /* One-loop functions */
    mutable double A0_Mz2_Mw2_cache[3];
    mutable double A0_Mz2_mh2_cache[3];
    mutable double A0_Mw2_Mz2_cache[3];
    mutable double A0_Mw2_mh2_cache[3];
    mutable double A0_Mz2_Mz2_cache[2];
    mutable double A0_Mw2_Mw2_cache[2];

    ////////////////////////////////////////////////////////////////////////     
    // Caches for complex variables
    // The last two elements are the caches for a complex variable.
    
    mutable double B0_Mz2_Mw2_mh2_Mw2_cache[5];
    mutable double B0_Mz2_0_mh2_Mw2_cache[5];
    mutable double B0_Mw2_Mz2_Mt2_Mt2_cache[5];
    mutable double B0_Mz2_Mz2_Mw2_Mw2_cache[4];
    mutable double B0_Mz2_Mz2_mh2_Mz2_cache[4];
    mutable double B0_Mz2_Mw2_Mz2_Mw2_cache[4];
    mutable double B0_Mz2_Mw2_0_Mw2_cache[4];
    mutable double B0_Mz2_0_Mz2_Mw2_cache[4];
    mutable double B0_Mz2_0_0_Mw2_cache[4];
    mutable double B0_Mw2_Mz2_Mw2_Mw2_cache[4];
    mutable double B0_Mw2_Mw2_Mz2_Mw2_cache[4];
    mutable double B0_Mw2_Mw2_mh2_Mw2_cache[4];
    mutable double B0_Mw2_Mw2_0_Mw2_cache[3];
    mutable double B0_Mz2_Mz2_ml2_ml2_cache[6][4];
    mutable double B0_Mz2_Mz2_mq2_mq2_cache[6][4];
    
    mutable double B0p_Mz2_0_mh2_Mw2_cache[5];
    mutable double B0p_Mz2_Mz2_mh2_Mz2_cache[4];
    mutable double B0p_Mz2_0_Mz2_Mw2_cache[4];
    mutable double B0p_Mz2_Mz2_Mw2_Mw2_cache[4];
    mutable double B0p_Mw2_Mw2_Mz2_Mw2_cache[4];
    mutable double B0p_Mw2_Mw2_mh2_Mw2_cache[4];
    mutable double B0p_Mw2_Mw2_0_Mw2_cache[3];
    mutable double B0p_Mz2_Mz2_ml2_ml2_cache[6][4];
    mutable double B0p_Mz2_Mz2_mq2_mq2_cache[6][4];
    
    mutable double B1_Mz2_0_ml2_mlprime2_cache[3][5];
    mutable double B1_Mz2_0_mq2_mqprime2_cache[3][5];
    mutable double B1_Mz2_0_mlprime2_ml2_cache[3][5];
    mutable double B1_Mz2_0_mqprime2_mq2_cache[3][5];
    mutable double B1_Mz2_Mw2_ml2_mlprime2_cache[3][6];
    mutable double B1_Mz2_Mw2_mq2_mqprime2_cache[3][6];
    mutable double B1_Mz2_Mw2_mlprime2_ml2_cache[3][6];
    mutable double B1_Mz2_Mw2_mqprime2_mq2_cache[3][6];
    
    mutable double B1p_Mw2_Mw2_ml2_mlprime2_cache[3][5];
    mutable double B1p_Mw2_Mw2_mq2_mqprime2_cache[3][5];
    mutable double B1p_Mw2_Mw2_mlprime2_ml2_cache[3][5];
    mutable double B1p_Mw2_Mw2_mqprime2_mq2_cache[3][5];
        
    mutable double Bf_Mz2_Mz2_ml2_ml2_cache[6][4];
    mutable double Bf_Mz2_Mz2_mq2_mq2_cache[6][4];
    mutable double Bf_Mz2_0_ml2_ml2_cache[6][4];
    mutable double Bf_Mz2_0_mq2_mq2_cache[6][4];
    mutable double Bf_Mz2_Mw2_mlprime2_ml2_cache[3][6];
    mutable double Bf_Mz2_Mw2_mqprime2_mq2_cache[3][6];
    mutable double Bf_Mz2_0_mlprime2_ml2_cache[3][5];
    mutable double Bf_Mz2_0_mqprime2_mq2_cache[3][5];
    mutable double Bf_Mw2_Mw2_mlprime2_ml2_cache[3][5];
    mutable double Bf_Mw2_Mw2_mqprime2_mq2_cache[3][5];

    mutable double Bfp_Mz2_Mz2_ml2_ml2_cache[6][4];
    mutable double Bfp_Mz2_Mz2_mq2_mq2_cache[6][4];
    mutable double Bfp_Mw2_Mw2_mlprime2_ml2_cache[3][5];
    mutable double Bfp_Mw2_Mw2_mqprime2_mq2_cache[3][5];
    
    mutable double C0_Mz2_Mw2_Mt2_Mw2_cache[5];
    mutable double C0_Mz2_Mt2_Mw2_Mt2_cache[5];
    mutable double C0_Mz2_0_Mw2_0_cache[4];
    mutable double C0_Mz2_Mw2_0_Mw2_cache[4];
    mutable double C0_Mw2_Mw2_0_Mz2_cache[4];
    mutable double C0_Mw2_0_Mz2_0_cache[4];
    mutable double C0_Mz2_0_Mz2_0_cache[3];

    
    ////////////////////////////////////////////////////////////////////////     

    bool CacheCheck(const double cache[], 
                    const int NumPar, const double params[]) const 
    {
        if (!bUseCacheEWSMcache) return false;
        bool bCache = true;
        for(int i=0; i<NumPar; ++i)
            bCache &= (params[i] == cache[i]);
        return bCache;
    }

    void newCacheForDouble(double cache[], const int NumPar, 
                           const double params[], const double newResult) const 
    {
        if (!bUseCacheEWSMcache) return;
        for(int i=0; i<NumPar; ++i)
            cache[i] = params[i];
        cache[NumPar] = newResult;
    }
    
    void newCacheForComplex(double cache[], const int NumPar, 
                            const double params[], const complex newResult) const 
    {
        if (!bUseCacheEWSMcache) return;
        for(int i=0; i<NumPar; ++i)
            cache[i] = params[i];
        cache[NumPar] = newResult.real();
        cache[NumPar+1] = newResult.imag();
    }
    
};

#endif	/* EWSMCACHE_H */

