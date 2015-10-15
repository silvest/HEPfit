/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef SUSYMATCHING_H
#define	SUSYMATCHING_H

#include <gslpp.h>
#include <complex>
#include <PVfunctions.h>
#include "StandardModelMatching.h"

#define LEPS 1.e-5     // tolerance in the limit of S(x,y) to S(x) 
#define SUSYLEPS 5.e-3 // tolerance in the limits of D0[x,y,z,t] and D2[x,y,z,t]
#define SUSYLEPS2 1.e-10 // tolerance in the limits of B0[x,y]
#define SUSYLEPS3 5.e-3 // tolerance in the limits of D0N

class SUSY;

/**
 * @class SUSYMatching
 * @ingroup SUSY
 * @brief A class for the matching in the MSSM. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details  
 */
class SUSYMatching : public StandardModelMatching {
public:
    SUSYMatching(const SUSY & SUSY_i);
    
    /**
     * 
     * @brief \f$ \Delta B = 2 \f$, \f$ B_{d} \f$ 
     * @return return the vector of SUSY Wilson coefficients
     */
    virtual  std::vector<WilsonCoefficient>& CMdbd2();
    
    /**
     * 
     * @brief \f$ \Delta B = 2 \f$, \f$ B_{s} \f$ 
     * @return return the vector of SUSY Wilson coefficients
     */
    virtual  std::vector<WilsonCoefficient>& CMdbs2();
    
//    /**
//     * 
//     * @brief \f$ \Delta C = 2 \f$,
//     * @return return the vector of SUSY Wilson coefficients
//     */
    virtual  std::vector<WilsonCoefficient>& CMdd2();
    
    /**
     * 
     * @brief \f$ \Delta S = 2 \f$ 
     * @return return the vector of SUSY Wilson coefficients
     */
    virtual  std::vector<WilsonCoefficient>& CMdk2();
    
//    /** 
//     * 
//     * @brief operator basis: 
//     * @brief current current; qcd penguins; magnetic and chromomagnetic penguins; semileptonic  
//     * @param a, a=0 deltaS=0 deltaC=0;  a=1 deltaS=1 deltaC=0;
//     * @return Wilson coefficients Buras base for non-leptonic B decays 
//     */
//    virtual  std::vector<WilsonCoefficient>& CMbnlep(const int& a);
//    
//    /**
//     * 
//     * @brief operator basis: - current current opertors  
//     * @param a, a=0 deltaS=0 deltaC=0;  a=1 1,0 ;  a=2 0,1 ; a=3 1,1 
//     * @return Wilson coefficients, Buras basis, for non-leptonic B decays 
//     */
//    virtual  std::vector<WilsonCoefficient>& CMbnlepCC(const int& a);
    
    /** 
     * 
     * @brief operator basis: current current; qcd penguins; 
     * magnetic and chromomagnetic penguins; semileptonic 
     * @return Wilson coefficients, Misiak basis, for \f$ B \rightarrow X_{s} \gamma, l^{+} l{-} \f$
     */
    virtual  std::vector<WilsonCoefficient>& CMbsg();

    /**
     *
     * @brief li -> lj gamma
     * @return
     */
    virtual std::vector<WilsonCoefficient>& CMDLij(int li_lj);

    /**
     *
     * @brief li -> lj lj lj
     * @return
     */
    virtual std::vector<WilsonCoefficient>& CMDLi3j(int li_lj);

//   /** Calculates amplitudes for m->(3)e (1), t->(3)m (2) and t->(3)e (3)**/
    virtual gslpp::vector<gslpp::complex> AFunctions(int n);

//   /** Calculates amplitudes for m->3e (1), t->3m (2) and t->3e (3)**/
    virtual gslpp::vector<gslpp::complex> BFunctions(int n);

//   /** Calculates amplitudes for m->3e (1), t->3m (2) and t->3e (3)**/
    virtual gslpp::vector<gslpp::complex> FFunctions(int n);

   /** Calculates C7 and C7' for m->e (1), t->m (2) and t->e (3)**/
    virtual gslpp::vector<gslpp::complex> C7_Lepton(int n);

   /** Calculates C9 and C9' for m->eee (1), t->mmm (2) and t->eee (3)**/
    virtual gslpp::vector<gslpp::complex> C9_Lepton(int n);

   /** Calculates C10 and C10' for m->eee (1), t->mmm (2) and t->eee (3)**/
    virtual gslpp::vector<gslpp::complex> C10_Lepton(int n);

   /** Calculates CS and CS' for m->eee (1), t->mmm (2) and t->eee (3)**/
    virtual gslpp::vector<gslpp::complex> CS_Lepton(int n);

    
//    /**
//     * 
//     * @brief current-current oerators, Misiak basis
//     * @return Wilson coefficients for \f$ D^{0} \rightarrow \pi \pi , K K \f$
//     */
//    virtual  std::vector<WilsonCoefficient>& CMd1();
//    
//    /**
//     * 
//     * @brief current-current oerators, Buras basis
//     * @return Wilson coefficients for \f$ D^{0} \rightarrow \pi \pi , K K \f$
//     */
//    virtual  std::vector<WilsonCoefficient>& CMd1Buras();
    
    /**
     *
     * @brief Updates to new SUSY parameter sets.
     * @return
     */
    
    void updateSUSYParameters();

   
    /** Large tan beta correction to CKM, down-quark mass **/
    
    void Comp_DeltaMd();
    void Comp_mySUSY_CKM();
    
    /** Quark Masses run to the SUSY scale Q **/
    void Comp_mySUSYMQ();
    
    /** Calcolous of Charged Higgs - squrks down - squarks up vertex **/
    
    void Comp_VUDHH();
    
    /** Calcolous of down quark - down squark - neutralino vertex **/
    
    void Comp_VdDNL(int flag);
    void Comp_VdDNR(int flag);
    
    /** Calcolous of down quark - up squark - chargino vertex **/

    void Comp_VdUCL();
    void Comp_VdUCR(int flag);     
    
    
    /** Calcolous of up quark - up squark - neutralino vertex **/
    
    void Comp_VuUN();
    
    /** Calcolous of up quark - down squark - chargino vertex **/
    
    void Comp_VuDCL();
    void Comp_VuDCR();
    
    
    /** Tan beta corrections to down quark masses **/
    
    void Comp_Eps_J();
    void Comp_Lambda0EpsY();
    
    /** Calcolous of charged Higgs - quark vertex **/
    
    void Comp_PHLR();
    void Comp_PHRL();
    
    void Comp_DeltaDL();

    
    
private:
    const SUSY & mySUSY;
    const PVfunctions PV;

    WilsonCoefficient mcdbd2, mcdbd2Hp, mcdbd2gg, mcdbd2ChiChi, mcdbd2Chi0Chi0, mcdbd2Chi0g,
                      mcdbd2HpT, mcdbd2ggT, mcdbd2ChiChiT, mcdbd2Chi0Chi0T, mcdbd2Chi0gT;
    WilsonCoefficient mcdbs2, mcdbs2Hp, mcdbs2gg, mcdbs2ChiChi, mcdbs2Chi0Chi0, mcdbs2Chi0g,
                      mcdbs2HpT, mcdbs2ggT, mcdbs2ChiChiT, mcdbs2Chi0Chi0T, mcdbs2Chi0gT;
    WilsonCoefficient mcdk2, mcdk2Hp, mcdk2gg, mcdk2ChiChi, mcdk2Chi0Chi0, mcdk2Chi0g,
                      mcdk2HpT, mcdk2ggT, mcdk2ChiChiT, mcdk2Chi0Chi0T, mcdk2Chi0gT;
    WilsonCoefficient mcdd2, mcdd2Hp, mcdd2gg, mcdd2ChiChi, mcdd2Chi0Chi0, mcdd2Chi0g,
                      mcdd2HpT, mcdd2ggT, mcdd2ChiChiT, mcdd2Chi0Chi0T, mcdd2Chi0gT;
    WilsonCoefficient mcDLij,mcDLi3j;
    
    WilsonCoefficient mcbsg, mcbnlep, mcbnlepCC, mcd1, mcd1Buras;
    
    std::vector<WilsonCoefficient> vmdbd2;                           
    std::vector<WilsonCoefficient> vmdbs2;
    std::vector<WilsonCoefficient> vmdk2;
    std::vector<WilsonCoefficient> vmdd2;
    std::vector<WilsonCoefficient> vmcbsg;
    std::vector<WilsonCoefficient> vmDLij;
    std::vector<WilsonCoefficient> vmDLi3j;
    
    gslpp::matrix<gslpp::complex> myCKM;
    gslpp::matrix<gslpp::complex> myRu;
    gslpp::matrix<gslpp::complex> myRd;
    gslpp::matrix<gslpp::complex> myRl;
    gslpp::matrix<gslpp::complex> myRn;
//    gslpp::matrix<double> mym_sn_sq; //Why?
//    gslpp::matrix<double> mym_se_sq; //Why?
    gslpp::vector<double> mym_sn_sq;
    gslpp::vector<double> mym_se_sq;
    gslpp::vector<double> MChi0;
    gslpp::vector<double> MChi;
    gslpp::matrix<gslpp::complex> myN;
    gslpp::matrix<gslpp::complex> myV;
    gslpp::matrix<gslpp::complex> myU;

    gslpp::matrix<gslpp::complex> myNreal;
    gslpp::matrix<double> Lepty;
    gslpp::matrix<double> Leptz;
    gslpp::matrix<double> Leptfa1;
    gslpp::matrix<double> Leptfa2;
    gslpp::matrix<double> Leptf1;
    gslpp::matrix<double> Leptf2;
    gslpp::matrix<double> Leptf3;
    gslpp::matrix<double> Leptf4;
    gslpp::matrix<gslpp::complex> CRlE;
    gslpp::matrix<gslpp::complex> CRlMU;
    gslpp::matrix<gslpp::complex> CRlTAU;
    gslpp::matrix<gslpp::complex> CLlE;
    gslpp::matrix<gslpp::complex> CLlMU;
    gslpp::matrix<gslpp::complex> CLlTAU;
    gslpp::matrix<gslpp::complex> NRlE;
    gslpp::matrix<gslpp::complex> NRlMU;
    gslpp::matrix<gslpp::complex> NRlTAU;
    gslpp::matrix<gslpp::complex> NLlE;
    gslpp::matrix<gslpp::complex> NLlMU;
    gslpp::matrix<gslpp::complex> NLlTAU;
    gslpp::matrix<gslpp::complex> AmpA1LN;
    gslpp::matrix<gslpp::complex> AmpA1RN;
    gslpp::matrix<gslpp::complex> AmpA1LC;
    gslpp::matrix<gslpp::complex> AmpA1RC;
    gslpp::matrix<gslpp::complex> AmpTauA1LN;
    gslpp::matrix<gslpp::complex> AmpTauA1RN;
    gslpp::matrix<gslpp::complex> AmpTauA1LC;
    gslpp::matrix<gslpp::complex> AmpTauA1RC;
    gslpp::matrix<gslpp::complex> AmpTEA1LN;
    gslpp::matrix<gslpp::complex> AmpTEA1RN;
    gslpp::matrix<gslpp::complex> AmpTEA1LC;
    gslpp::matrix<gslpp::complex> AmpTEA1RC;
    gslpp::matrix<gslpp::complex> AmpALN;
    gslpp::matrix<gslpp::complex> AmpARN;
    gslpp::matrix<gslpp::complex> AmpALC;
    gslpp::matrix<gslpp::complex> AmpARC;
    gslpp::matrix<gslpp::complex> AmpTauALN;
    gslpp::matrix<gslpp::complex> AmpTauARN;
    gslpp::matrix<gslpp::complex> AmpTauALC;
    gslpp::matrix<gslpp::complex> AmpTauARC;
    gslpp::matrix<gslpp::complex> AmpTEALN;
    gslpp::matrix<gslpp::complex> AmpTEARN;
    gslpp::matrix<gslpp::complex> AmpTEALC;
    gslpp::matrix<gslpp::complex> AmpTEARC;
    std::complex<double> Leptfzn[6][4][4];
    std::complex<double> Leptgzn[6][4][4];
    std::complex<double> Leptfzc[3][2][2];
    std::complex<double> Leptgzc[3][2][2];

    double Q_S;
    double mu2R;
    double Als;
    double Mg;
    double tanb;
    double sinb;
    double cosb;
    double v;
    double v1;
    double v2;
    double gW;
    
    
    /** Passarino Veltman functions, Dk with k = 0,2 , and their limits **/
    
    double Dk(double x, double y, double z, double t, int k);

    double D0N(double x, double y, double z, double t);
    double D2LL0(double a, double b);
    double DL0(double a, double b, double c, int k);
    double DL(double a, double b, double c, int k);
    double DLL(double a, double b, int k);
    double DLLp(double a, double b, int k);
    double DLLL(double a, int k);
    
    
    /** Calcolous of Charged Higgs contributions **/
    gslpp::vector<gslpp::complex> CdF2dHp(int b, int q, int Dmixingflag);
    
    /** Calcolous Gluinos Box **/
    gslpp::vector<gslpp::complex> CdF2dgg(int b, int q, int Dmixingflag);
    
    /** Calcolous Neutralino - Gluino Box **/
    gslpp::vector<gslpp::complex> CdF2dChi0g(int b, int q, int Dmixingflag);
    
    /** Calcolous Charginos contribution to Wilson Coefficents **/
    gslpp::vector<gslpp::complex> CdF2dChiChi(int b, int q, int Dmixingflag);
    
    /** Calcolous Neutralino contribution to Wilson Coefficents **/
    gslpp::vector<gslpp::complex> CdF2dChi0Chi0(int b, int q, int Dmixingflag);

    gslpp::vector<gslpp::complex> AFunctions();
    gslpp::vector<gslpp::complex> BFunctions();
    gslpp::vector<gslpp::complex> FFunctions();
    gslpp::vector<gslpp::complex> C7_Lepton();
    gslpp::vector<gslpp::complex> C9_Lepton();
    gslpp::vector<gslpp::complex> C10_Lepton();
    gslpp::vector<gslpp::complex> CS_Lepton();
    
    
    /** Feynmann rule for the Chargino - down quark - up squarks vertex with tan beta correction  **/
    
    gslpp::complex VdUCL_cache[3][6][2];
    gslpp::complex VdUCR_cache[3][6][2][2];
    
    gslpp::complex VdUCL(int b, int k, int j);
    gslpp::complex VdUCR(int b, int k, int j, int flag);
    
    /** Feynmann rule for the Neutralino - down quark - down squarks vertex with tan beta correction  **/
    
    gslpp::complex VdDNL_cache[3][6][4][2];
    gslpp::complex VdDNR_cache[3][6][4][2];
    
    gslpp::complex VdDNL(int b, int k, int j, int flag);
    gslpp::complex VdDNR(int b, int k, int j, int flag);
    
    /** Feynmann rule for the Chargino - up quark - down squarks vertex with tan beta correction  **/
    
    gslpp::complex VuDCL_cache[3][6][2];
    gslpp::complex VuDCR_cache[3][6][2];
    
    gslpp::complex VuDCL(int b, int k, int j);
    gslpp::complex VuDCR(int b, int k, int j);
    
    /** Vertices Neutralino (N) - up quark (u) - up squarks (U) from Buras arXiv:hep-ph/0210145v2 
        in SLHA convention usefull in D - Dbar mixing **/
    
    gslpp::complex VuUNL_cache[3][6][4];
    gslpp::complex VuUNR_cache[3][6][4];
    
    gslpp::complex VuUN(int b, int k, int j, const std::string);
  
  
    /** Rules to include in the same formulas the D - D bar mixing  **/
    gslpp::complex VdUCL(int b, int k, int j, int Dmixingflag);
    gslpp::complex VdUCR(int b, int k, int j, int flag, int Dmixingflag);
        
    gslpp::complex VdDNL(int b, int k, int j, int flag, int Dmixingflag);
    gslpp::complex VdDNR(int b, int k, int j, int flag, int Dmixingflag);
    
    
    
    
 
    double CL(double a, double b, int k);
    double CLL(double a, int k);
    double Ck(double x, double y, double z,int k);
    double BL(double a, int k);
    double Bk(double x, double y, int k);
    
    gslpp::complex DeltaMd(int J, int I);
    
    gslpp::vector<gslpp::complex> Eps_JCache;
    gslpp::complex Eps_J(int J);
    
    gslpp::matrix<gslpp::complex> Lambda0EpsYCache;
    
    gslpp::complex Lambda0EpsY(int J, int I);
    
    gslpp::matrix<gslpp::complex> DeltaDL_Cache;
    
    
    gslpp::complex DeltaDL(int J, int I);
    gslpp::complex DeltaDR(int J, int I);
    
    
    gslpp::complex xdS(int S);
    gslpp::complex xuS(int S);
    gslpp::complex XRLS(int J, int I, int S);
    gslpp::complex XLRS(int J, int I, int S);
    
    
   
    gslpp::matrix<gslpp::complex> PHLRCache;
    
    
    gslpp::complex PHLR(int j, int i);
    gslpp::complex DeltaFHL(int j, int i);
    
    gslpp::matrix<gslpp::complex> PHRLCache;
    
    gslpp::complex PHRL(int j, int i);
    
    
    gslpp::complex VUDHH(int i, int j);
    gslpp::complex PGLR(int j, int i);
    gslpp::complex PGRL(int j, int i);
    gslpp::complex PLRk(int j, int i, int k);
    gslpp::complex PRLk(int j, int i, int k);
    
    gslpp::complex PRLk(int j, int i, int k, int Dmixingflag);
    gslpp::complex PLRk(int j, int i, int k, int Dmixingflag);
    
    gslpp::matrix<gslpp::complex> mySUSY_CKM();
    
    gslpp::complex EpsPrime(int J, int I);
    gslpp::vector<gslpp::complex> CalcC7(int b, int q);
    double F7k(double x, int k);
   
    
    gslpp::matrix<gslpp::complex> myCKM_cache;
    gslpp::matrix<gslpp::complex> VUDHH_cache;
    gslpp::matrix<gslpp::complex> DeltaMd_cache;
    gslpp::vector<double> mySUSYMQ;
    
};

#endif	/* SUSYMATCHING_H */
