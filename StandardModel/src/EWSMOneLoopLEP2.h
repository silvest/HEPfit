/* 
 * File:   EWSMOneloopLEP2.h
 * Author: giovannigrilli
 */

#ifndef EWSMONELOOPLEP2_H
#define	EWSMONELOOPLEP2_H

#include "EWSMOneLoopEW.h"
using namespace gslpp;


class EWSMOneLoopLEP2{
public:
    
    /**
     * @brief EWSMOneLoopLEP2 constructor
     * @param[in] cache_i reference to an EWSMcache object
     */
    EWSMOneLoopLEP2(const EWSMcache& cache_i,const StandardModel& SM_i);
    
    /**
     * @brief gamma propagator
     * @param[in] mu renormalization scale
     * @param[in] s momentum squared
     * @param[in] Mw_i the W-boson mass
     * @return Chi_gamma
     */
    complex Chi_gamma(const double mu, const double s, const double Mw_i,
                      const double W, const double X, const double Y) const;
    
    /**
     * @brief Z propagator
     * @param[in] mu renormalization scale
     * @param[in] s momentum squared
     * @param[in] Mw_i the W-boson mass
     * @return Chi_Z
     */
    complex Chi_Z(const double mu, const double s, const double Mw_i,
                  const double W, const double X, const double Y) const;
    
    /**
     * @brief gamma-Z propagator
     * @param[in] mu renormalization scale
     * @param[in] s momentum squared
     * @param[in] Mw_i the W-boson mass
     * @return Chi_gammaZ
     */
    complex Chi_gammaZ(const double mu, const double s, const double Mw_i,
                       const double W, const double X, const double Y) const;
    
    /**
     * @brief couplings constant for quarks to the Z gauge boson
     * @param[in] q quark 
     * @param[in] rho  the elicity
     * @param[in] Mw_i the W-boson mass
     * @return Chi_gammaZ
     */
    double g_rhofq(const QCD::quark q, const double rho, const double Mw_i) const;
    
    /**
     * @brief couplings constant for lepton to the Z gauge boson
     * @param[in] l lepton
     * @param[in] rho  the elicity
     * @param[in] Mw_i the W-boson mass
     * @return Chi_gammaZ
     */
    double g_rhofl(const StandardModel::lepton l, const double rho, const double Mw_i) const;
    
    /**
     * @brief couplings constant for electron to the Z gauge boson
     * @param[in] rho  the elicity
     * @param[in] Mw_i the W-boson mass
     * @return g_rhoe
     */
    double g_rhoe(const double rho, const double Mw_i) const;
    
    
    
    /**
     * @brief Mandelstam variable t
     * @param[in] mf fermion mass
     * @param[in] s momentum squared
     * @param[in] theta scattering angle
     * @return t
     */
    double t(const double mf, const double s,const double theta) const;
    
    /**
     * @brief Mandelstam variable u
     * @param[in] mf fermion mass
     * @param[in] s momentum squared
     * @param[in] theta scattering angle
     * @return u
     */
    double u(const double mf, const double s,const double theta) const;
    
    /**
     * @brief coefficient for the e+e- to leptons dressed Born amplitude
     * @param[in] l lepton
     * @param[in] rho  the elicity
     * @param[in] s  the momentum squared
     * @param[in] Mw_i  the W mass
     * @return Chi_gammaZ
     */
    complex Al(const double mu,const StandardModel::lepton l, const double rho, const double k, const double s, const double Mw_i,
                      const double W, const double X, const double Y) const;
    
    /**
     * @brief coefficient for the e+e- to quarks dressed Born amplitude
     * @param[in] mu renormalization scale
     * @param[in] q quark
     * @param[in] rho  the elicity
     * @param[in] s  the momentum squared
     * @param[in] Mw_i  the W mass
     * @return Chi_gammaZ
     */
    complex Aq(const double mu,const QCD::quark q, const double rho, const double k, const double s, const double Mw_i,
                      const double W, const double X, const double Y) const;
    
    /**
     * @brief coefficient for the e+e- to leptons self energies corrections
     * @param[in] mu renormalization scale
     * @param[in] l lepton
     * @param[in] rho  the elicity
     * @param[in] k  the elicity
     * @param[in] s  the momentum squared
     * @param[in] Mw_i  the W mass
     * @return Chi_gammaZ
     */
    complex Bl(const double mu,const StandardModel::lepton l, const double rho, const double k, const double s, const double Mw_i,
                      const double W, const double X, const double Y) const;
    
    /**
     * @brief coefficient for the e+e- to quarks self energies corrections
     * @param[in] mu renormalization scale
     * @param[in] q quark
     * @param[in] rho  the elicity
     * @param[in] k  the elicity
     * @param[in] s  the momentum squared
     * @param[in] Mw_i  the W mass
     * @return Chi_gammaZ
     */
    complex Bq(const double mu,const QCD::quark q, const double rho, const double k, const double s, const double Mw_i,
                      const double W, const double X, const double Y) const;
    
    /**
     * @brief invariant function Lambda2
     * @param[in] m the mass 
     * @param[in] s  the momentum squared
     * @return Lambda2
     */
    complex Lambda2(const double m, const double s) const;
    
    /**
     * @brief invariant function Lambda3
     * @param[in] m the mass
     * @param[in] s  the momentum squared
     * @return Lambda3
     */
    complex Lambda3(const double m, const double s) const;
    
    
    /**
     * @brief factor for the contribution M_eegamma
     * @param[in] s  the momentum squared
     * @param[in] Mw_i  the W mass
     * @return FLgammal
     */
    complex FLgammal(const double s, const double Mw_i) const;
    
    /**
     * @brief factor for the contribution M_eeZ
     * @param[in] s  the momentum squared
     * @param[in] Mw_i  the W mass
     * @return FLZl
     */
    complex FLZl(const double s, const double Mw_i) const;
    
    
    /**
     * @brief coefficient for weak vertex eegamma correction 
     * @param[in] mu renormalization scale
     * @param[in] l lepton
     * @param[in] rho  the elicity
     * @param[in] k  the elicity
     * @param[in] s  the momentum squared
     * @param[in] Mw_i  the W mass
     * @return Cl
     */
    complex Cl(const double mu,const StandardModel::lepton l, const double rho, const double k, const double s, const double Mw_i,
                      const double W, const double X, const double Y) const;
    /**
     * @brief coefficient for weak vertex eegamma correction 
     * @param[in] mu renormalization scale
     * @param[in] q lepton
     * @param[in] rho  the elicity
     * @param[in] k  the elicity
     * @param[in] s  the momentum squared
     * @param[in] Mw_i  the W mass
     * @return Cq
     */
    complex Cq(const double mu,const QCD::quark q, const double rho, const double k, const double s, const double Mw_i,
                      const double W, const double X, const double Y) const;
    /**
     * @brief coefficient for weak vertex eegamma correction 
     * @param[in] mu renormalization scale
     * @param[in] l lepton
     * @param[in] rho  the elicity
     * @param[in] k  the elicity
     * @param[in] s  the momentum squared
     * @param[in] Mw_i  the W mass
     * @return Dl_rho
     */
    complex Dl_rho(const double mu,const StandardModel::lepton l, const double rho, const double k, const double s, const double Mw_i,
                      const double W, const double X, const double Y) const;
    /**
     * @brief coefficient for weak vertex eegamma correction 
     * @param[in] mu renormalization scale
     * @param[in] l lepton
     * @param[in] rho  the elicity
     * @param[in] k  the elicity
     * @param[in] s  the momentum squared
     * @param[in] Mw_i  the W mass
     * @return Dq_rho
     */
    complex Dq_rho(const double mu,const QCD::quark q, const double rho, const double k, const double s, const double Mw_i,
                      const double W, const double X, const double Y) const;
    
    
    /**
     * @brief vector couplings constant for e
     * @param[in] Mw_i the W-boson mass
     * @return ve
     */
    double ve(const double Mw_i) const;
    
    /**
     * @brief axial vector couplings constant for e
     * @param[in] Mw_i the W-boson mass
     * @return ae
     */
    double ae(const double Mw_i) const;
    
    /**
     * @brief vector couplings constant for q
     * @param[in] q quark
     * @param[in] Mw_i the W-boson mass
     * @return vq
     */
    double vq(const QCD::quark q, const double Mw_i) const;
    
    /**
     * @brief axial vector couplings constant for q
     * @param[in] q quark
     * @param[in] Mw_i the W-boson mass
     * @return aq
     */
    double aq(const QCD::quark q,const double Mw_i) const;
    
    /**
     * @brief vector couplings constant for q
     * @param[in] l lepton
     * @param[in] Mw_i the W-boson mass
     * @return vl
     */
    double vl(const StandardModel::lepton l, const double Mw_i) const;
    
    /**
     * @brief axial vector couplings constant for q
     * @param[in] l lepton
     * @param[in] Mw_i the W-boson mass
     * @return al
     */
    double al(const StandardModel::lepton l,const double Mw_i) const;
    
    /**
     * @brief coefficients C for the final vertex corrections
     * @param[in] mu the renormalization scale
     * @param[in] l lepton
     * @param[in] m1 mass
     * @param[in] m2 mass
     * @param[in] m3 mass
     * @param[in] s the momentum squared
     * @return C1plus_l
     */
    complex C1plus_l(const double mu, const StandardModel::lepton l,
                    const double m1, const double m2, const double m3,
                    const double s) const;
    
    /**
     * @brief coefficients C for the final vertex corrections
     * @param[in] mu the renormalization scale
     * @param[in] q quark
     * @param[in] m1 mass
     * @param[in] m2 mass
     * @param[in] m3 mass
     * @param[in] s the momentum squared
     * @return C1plus_q
     */
    complex C1plus_q(const double mu, const QCD::quark q,
                    const double m1, const double m2, const double m3,
                    const double s) const;
    
    
    /**
     * @brief coefficients C for the final vertex corrections
     * @param[in] mu the renormalization scale
     * @param[in] l lepton
     * @param[in] m1 mass
     * @param[in] m2 mass
     * @param[in] m3 mass
     * @param[in] s the momentum squared
     * @return C1minus_l
     */
    complex C1minus_l(const double mu, const StandardModel::lepton l,
                    const double m1, const double m2, const double m3,
                    const double s) const;
    
    /**
     * @brief coefficients C for the final vertex corrections
     * @param[in] mu the renormalization scale
     * @param[in] q quark
     * @param[in] m1 mass
     * @param[in] m2 mass
     * @param[in] m3 mass
     * @param[in] s the momentum squared
     * @return C1minus_q
     */
    complex C1minus_q(const double mu, const QCD::quark q,
                    const double m1, const double m2, const double m3,
                    const double s) const;
    
    /**
     * @brief coefficients C for the final vertex corrections
     * @param[in] mu the renormalization scale
     * @param[in] l lepton
     * @param[in] m1 mass
     * @param[in] m2 mass
     * @param[in] m3 mass
     * @param[in] s the momentum squared
     * @return C20_l
     */
    complex C20_l(const double mu, const StandardModel::lepton l,
                    const double m1, const double m2, const double m3,
                    const double s) const;
    
    /**
     * @brief coefficients C for the final vertex corrections
     * @param[in] mu the renormalization scale
     * @param[in] q quark
     * @param[in] m1 mass
     * @param[in] m2 mass
     * @param[in] m3 mass
     * @param[in] s the momentum squared
     * @return C20_q
     */
    complex C20_q(const double mu, const QCD::quark q,
                    const double m1, const double m2, const double m3,
                    const double s) const;
    
    /**
     * @brief coefficients C for the final vertex corrections
     * @param[in] mu the renormalization scale
     * @param[in] l lepton
     * @param[in] m1 mass
     * @param[in] m2 mass
     * @param[in] m3 mass
     * @param[in] s the momentum squared
     * @return C2plus_l
     */
    complex C2plus_l(const double mu, const StandardModel::lepton l,
                    const double m1, const double m2, const double m3,
                    const double s) const;
    
    /**
     * @brief coefficients C for the final vertex corrections
     * @param[in] mu the renormalization scale
     * @param[in] q quark
     * @param[in] m1 mass
     * @param[in] m2 mass
     * @param[in] m3 mass
     * @param[in] s the momentum squared
     * @return C2plus_q
     */
    complex C2plus_q(const double mu, const QCD::quark q,
                    const double m1, const double m2, const double m3,
                    const double s) const;
    /**
     * @brief coefficients C for the final vertex corrections
     * @param[in] mu the renormalization scale
     * @param[in] l lepton
     * @param[in] m1 mass
     * @param[in] m2 mass
     * @param[in] m3 mass
     * @param[in] s the momentum squared
     * @return C2minus_l
     */
    complex C2minus_l(const double mu, const StandardModel::lepton l,
                    const double m1, const double m2, const double m3,
                    const double s) const;
    
    /**
     * @brief coefficients C for the final vertex corrections
     * @param[in] mu the renormalization scale
     * @param[in] q quark
     * @param[in] m1 mass
     * @param[in] m2 mass
     * @param[in] m3 mass
     * @param[in] s the momentum squared
     * @return C2minus_q
     */
    complex C2minus_q(const double mu, const QCD::quark q,
                    const double m1, const double m2, const double m3,
                    const double s) const;
    
    
    /**
     * @brief coefficients C for the final vertex corrections
     * @param[in] mu the renormalization scale
     * @param[in] l lepton
     * @param[in] m1 mass
     * @param[in] m2 mass
     * @param[in] m3 mass
     * @param[in] s the momentum squared
     * @return C2plusminus_l
     */
    complex C2plusminus_l(const double mu, const StandardModel::lepton l,
                    const double m1, const double m2, const double m3,
                    const double s) const;
    
    /**
     * @brief coefficients C for the final vertex corrections
     * @param[in] mu the renormalization scale
     * @param[in] q quark
     * @param[in] m1 mass
     * @param[in] m2 mass
     * @param[in] m3 mass
     * @param[in] s the momentum squared
     * @return C2plusminus_q
     */
    complex C2plusminus_q(const double mu, const QCD::quark q,
                    const double m1, const double m2, const double m3,
                    const double s) const;
    
    /**
     * @brief abbreviation muf for the final vertex corrections
     * @param[in] Mw_i the W mass
     * @param[in] mf fermion mass
     * @return muf_q
     */
    double muf(const double Mw_i, const double mf) const;
    
    /**
     * @brief abbreviation muf for the final vertex corrections
     * @param[in] Mw_i the W mass
     * @param[in] l lepton
     * @return muf_l
     */
    //double muf_l(const double Mw_i, const StandardModel::lepton l) const;
    
    /**
     * @brief abbreviation alphaf_q for the final vertex corrections
     * @param[in] Mw_i the W mass
     * @param[in] q quark
     * @return alphaf_q
     */ 
    double alphaf_q(const double Mw_i,const QCD::quark q) const;
    
    /**
     * @brief abbreviation alphaf_q for the final vertex corrections
     * @param[in] Mw_i the W mass
     * @param[in] q quark
     * @return alphaf_q
     */ 
    double alphaf_qprime(const double Mw_i,const QCD::quark q) const;
    
    
    
    /**
     * @brief quark of the internal line for the final vertex corrections
     * @param[in] q quark
     * @return qprime
     */ 
    QCD::quark qprime(const QCD::quark q) const;
    
    
    /**
     * @brief abbreviation alphaf_l for the final vertex corrections
     * @param[in] Mw_i the W mass
     * @param[in] l lepton
     * @return alphaf_l
     */ 
    double alphaf_l(const double Mw_i,const StandardModel::lepton l) const;
    
    /**
     * @brief weak renormalization constant C1f for leptons
     * @param[in] l lepton
     * @param[in] Mw_i the W mass
     * @param[in] mu the renormalization scale
     * @return C1f_l
     */ 
    double C1f_l(const StandardModel::lepton l,const double Mw_i, const double mu) const;
    
    /**
     * @brief weak renormalization constant C1f for quark
     * @param[in] q quark
     * @param[in] Mw_i the W mass
     * @param[in] mu the renormalization scale
     * @return C1f_q
     */ 
    double C1f_q(const QCD::quark q,const double Mw_i, const double mu) const;
    
    /**
     * @brief mass of the internal quark line
     * @param[in] q quark
     * @return mqprime
     */ 
    //double mqprime(const QCD::quark q) const;
    
    /**
     * @brief mass of the internal quark line
     * @param[in] q quark
     * @return Qqprime
     */ 
    double Qqprime(const QCD::quark q) const;
    
    /**
     * @brief mass of the internal quark line
     * @param[in] q quark
     * @return I3qprime
     */ 
    double I3qprime(const QCD::quark q) const;
    
    /**
     * @brief weak renormalization constant C1f for leptons
     * @param[in] l lepton
     * @param[in] Mw_i the W mass
     * @param[in] mu the renormalization scale
     * @return C2f_l
     */ 
    double C2f_l(const StandardModel::lepton l,const double Mw_i, const double mu) const;
    
    /**
     * @brief weak renormalization constant C1f for quark
     * @param[in] q quark
     * @param[in] Mw_i the W mass
     * @param[in] mu the renormalization scale
     * @return C2f_q
     */ 
    double C2f_q(const QCD::quark q,const double Mw_i, const double mu) const;
    
    /**
     * @brief contribtion to vector form factor CV
     * @param[in] l lepton
     * @param[in] Mw_i the W mass
     * @param[in] mu the renormalization scale
     * @return CVgamma_l
     */ 
    double CVgammal(const StandardModel::lepton l, const double Mw_i, const double mu) const;
    
    /**
     * @brief contribtion to vector form factor CV
     * @param[in] l lepton
     * @param[in] Mw_i the W mass
     * @param[in] mu the renormalization scale
     * @return CVgamma_q
     */ 
    double CVgammaq(const QCD::quark q, const double Mw_i, const double mu) const;
    
    /**
     * @brief contribtion to axial vector form factor CA
     * @param[in] l lepton
     * @param[in] Mw_i the W mass
     * @param[in] mu the renormalization scale
     * @return CAgamma_l
     */ 
    double CAgammal(const StandardModel::lepton l, const double Mw_i, const double mu) const;
    
    /**
     * @brief contribtion to axial vector form factor CA
     * @param[in] l lepton
     * @param[in] Mw_i the W mass
     * @param[in] mu the renormalization scale
     * @return CAgamma_q
     */ 
    double CAgammaq(const QCD::quark q, const double Mw_i, const double mu) const;
    
    
    
    
    /**
     * @brief contribtion to vector form factor CV
     * @param[in] l lepton
     * @param[in] Mw_i the W mass
     * @param[in] mu the renormalization scale
     * @return CVZ_l
     */ 
    double CVZl(const StandardModel::lepton l, const double Mw_i, const double mu) const;
    
    /**
     * @brief contribtion to vector form factor CV
     * @param[in] l lepton
     * @param[in] Mw_i the W mass
     * @param[in] mu the renormalization scale
     * @return CVZ_q
     */ 
    double CVZq(const QCD::quark q, const double Mw_i, const double mu) const;
    
    /**
     * @brief contribtion to axial vector form factor CA
     * @param[in] l lepton
     * @param[in] Mw_i the W mass
     * @param[in] mu the renormalization scale
     * @return CAgamma_l
     */ 
    double CAZl(const StandardModel::lepton l, const double Mw_i, const double mu) const;
    
    /**
     * @brief contribtion to axial vector form factor CA
     * @param[in] q quark
     * @param[in] Mw_i the W mass
     * @param[in] mu the renormalization scale
     * @return CAgamma_q
     */ 
    double CAZq(const QCD::quark q, const double Mw_i, const double mu) const;
    
    
    
/////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////   VECTOR FORM FACTOR      /////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////
    
    
    
    /**
     * @brief class I vertex correction: vector form factor 'a' for lepton
     * @param[in] l lepton
     * @param[in] Mw_i the W mass
     * @param[in] mu the renormalization scale
     * @param[in] s the momentum squared
     * @return FI_Val
     */ 
    complex FI_Val(const StandardModel::lepton l, const double Mw_i, const double mu,
                  const double s) const;
    
    /**
     * @brief class I vertex correction: vector form factor 'a' for quark
     * @param[in] q quark
     * @param[in] Mw_i the W mass
     * @param[in] mu the renormalization scale
     * @param[in] s the momentum squared
     * @return FI_Vaq
     */ 
    complex FI_Vaq(const QCD::quark q, const double Mw_i, const double mu,
                  const double s) const;
    
    /**
     * @brief class I vertex correction: vector form factor 'b' for lepton
     * @param[in] l lepton
     * @param[in] Mw_i the W mass
     * @param[in] mu the renormalization scale
     * @param[in] s the momentum squared
     * @return FI_Val
     */ 
    complex FI_Vbl(const StandardModel::lepton l, const double Mw_i, const double mu,
                  const double s) const;
    
    /**
     * @brief class I vertex correction: vector form factor 'b' for quark
     * @param[in] q quark
     * @param[in] Mw_i the W mass
     * @param[in] mu the renormalization scale
     * @param[in] s the momentum squared
     * @return FI_Vaq
     */ 
    complex FI_Vbq(const QCD::quark q, const double Mw_i, const double mu,
                  const double s) const;
    
    /**
     * @brief class II vertex correction: vector form factor 'c' for lepton
     * @param[in] l lepton
     * @param[in] Mw_i the W mass
     * @param[in] mu the renormalization scale
     * @param[in] s the momentum squared
     * @return FI_Val
     */ 
    complex FII_Vcl(const StandardModel::lepton l, const double Mw_i, const double mu,
                  const double s) const;
    
    /**
     * @brief class II vertex correction: vector form factor 'c' for quark
     * @param[in] q quark
     * @param[in] Mw_i the W mass
     * @param[in] mu the renormalization scale
     * @param[in] s the momentum squared
     * @return FI_Vaq
     */ 
    complex FII_Vcq(const QCD::quark q, const double Mw_i, const double mu,
                  const double s) const;
    
    /**
     * @brief class III vertex correction: vector form factor 'd' for lepton
     * @param[in] l lepton
     * @param[in] Mw_i the W mass
     * @param[in] mu the renormalization scale
     * @param[in] s the momentum squared
     * @return FI_Val
     */ 
    complex FIII_Vdl(const StandardModel::lepton l, const double Mw_i, const double mu,
                  const double s) const;
    
    /**
     * @brief class II vertex correction: vector form factor 'd' for quark
     * @param[in] q quark
     * @param[in] Mw_i the W mass
     * @param[in] mu the renormalization scale
     * @param[in] s the momentum squared
     * @return FI_Vaq
     */ 
    complex FIII_Vdq(const QCD::quark q, const double Mw_i, const double mu,
                  const double s) const;
    
    /**
     * @brief class III vertex correction: vector form factor 'e' for lepton
     * @param[in] l lepton
     * @param[in] Mw_i the W mass
     * @param[in] mu the renormalization scale
     * @param[in] s the momentum squared
     * @return FIII_Vel
     */ 
    complex FIII_Vel(const StandardModel::lepton l, const double Mw_i, const double mu,
                  const double s) const;
    
    /**
     * @brief class II vertex correction: vector form factor 'd' for quark
     * @param[in] q quark
     * @param[in] Mw_i the W mass
     * @param[in] mu the renormalization scale
     * @param[in] s the momentum squared
     * @return FI_Vaq
     */ 
    complex FIII_Veq(const QCD::quark q, const double Mw_i, const double mu,
                  const double s) const;
    
    
    /**
     * @brief class III vertex correction: vector form factor 'd' for lepton
     * @param[in] l lepton
     * @param[in] Mw_i the W mass
     * @param[in] mu the renormalization scale
     * @param[in] s the momentum squared
     * @return FI_Val
     */ 
    complex FIII_Vfl(const StandardModel::lepton l, const double Mw_i, const double mu,
                  const double s) const;
    
    /**
     * @brief class II vertex correction: vector form factor 'd' for quark
     * @param[in] q quark
     * @param[in] Mw_i the W mass
     * @param[in] mu the renormalization scale
     * @param[in] s the momentum squared
     * @return FI_Vaq
     */ 
    complex FIII_Vfq(const QCD::quark q, const double Mw_i, const double mu,
                  const double s) const;
    
    /**
     * @brief class IV vertex correction: vector form factor 'g' for lepton
     * @param[in] l lepton
     * @param[in] Mw_i the W mass
     * @param[in] mu the renormalization scale
     * @param[in] s the momentum squared
     * @return FIV_Vgl
     */ 
    complex FIV_Vgl(const StandardModel::lepton l, const double Mw_i, const double mu,
                  const double s) const;
    
    /**
     * @brief class IV vertex correction: vector form factor 'g' for quark
     * @param[in] q quark
     * @param[in] Mw_i the W mass
     * @param[in] mu the renormalization scale
     * @param[in] s the momentum squared
     * @return FIV_Vgq
     */ 
    complex FIV_Vgq(const QCD::quark q, const double Mw_i, const double mu,
                  const double s) const;
    
    /**
     * @brief class V vertex correction: vector form factor 'h' for lepton
     * @param[in] l lepton
     * @param[in] Mw_i the W mass
     * @param[in] mu the renormalization scale
     * @param[in] s the momentum squared
     * @return FV_Vhl
     */ 
    complex FV_Vhl(const StandardModel::lepton l, const double Mw_i, const double mu,
                  const double s) const;
    
    /**
     * @brief class V vertex correction: vector form factor 'h' for quark
     * @param[in] q quark
     * @param[in] Mw_i the W mass
     * @param[in] mu the renormalization scale
     * @param[in] s the momentum squared
     * @return FV_Vhq
     */ 
    complex FV_Vhq(const QCD::quark q, const double Mw_i, const double mu,
                  const double s) const;
    
    /**
     * @brief class VI vertex correction: vector form factor 'i' for lepton
     * @param[in] l lepton
     * @param[in] Mw_i the W mass
     * @param[in] mu the renormalization scale
     * @param[in] s the momentum squared
     * @return FVI_Vil
     */ 
    complex FVI_Vil(const StandardModel::lepton l, const double Mw_i, const double mu,
                  const double s) const;
    
    /**
     * @brief class VI vertex correction: vector form factor 'i' for quark
     * @param[in] q quark
     * @param[in] Mw_i the W mass
     * @param[in] mu the renormalization scale
     * @param[in] s the momentum squared
     * @return FVI_Viq
     */ 
    complex FVI_Viq(const QCD::quark q, const double Mw_i, const double mu,
                  const double s) const;
    
    
    /**
     * @brief class I vertex correction: vector form factor 'j' for lepton
     * @param[in] l lepton
     * @param[in] Mw_i the W mass
     * @param[in] mu the renormalization scale
     * @param[in] s the momentum squared
     * @return FI_Vjl
     */ 
    complex FI_Vjl(const StandardModel::lepton l, const double Mw_i, const double mu,
                  const double s) const;
    
    /**
     * @brief class I vertex correction: vector form factor 'a' for quark
     * @param[in] q quark
     * @param[in] Mw_i the W mass
     * @param[in] mu the renormalization scale
     * @param[in] s the momentum squared
     * @return FI_Vaq
     */ 
    complex FI_Vjq(const QCD::quark q, const double Mw_i, const double mu,
                  const double s) const;
    
    /**
     * @brief class I vertex correction: vector form factor 'k' for lepton
     * @param[in] l lepton
     * @param[in] Mw_i the W mass
     * @param[in] mu the renormalization scale
     * @param[in] s the momentum squared
     * @return FI_Vkl
     */ 
    complex FI_Vkl(const StandardModel::lepton l, const double Mw_i, const double mu,
                  const double s) const;
    
    /**
     * @brief class I vertex correction: vector form factor 'k' for quark
     * @param[in] q quark
     * @param[in] Mw_i the W mass
     * @param[in] mu the renormalization scale
     * @param[in] s the momentum squared
     * @return FI_Vkq
     */ 
    complex FI_Vkq(const QCD::quark q, const double Mw_i, const double mu,
                  const double s) const;
    
    /**
     * @brief class II vertex correction: vector form factor 'l' for lepton
     * @param[in] l lepton
     * @param[in] Mw_i the W mass
     * @param[in] mu the renormalization scale
     * @param[in] s the momentum squared
     * @return FI_Val
     */ 
    complex FII_Vll(const StandardModel::lepton l, const double Mw_i, const double mu,
                  const double s) const;
    
    /**
     * @brief class II vertex correction: vector form factor 'l' for quark
     * @param[in] q quark
     * @param[in] Mw_i the W mass
     * @param[in] mu the renormalization scale
     * @param[in] s the momentum squared
     * @return FI_Vlq
     */ 
    complex FII_Vlq(const QCD::quark q, const double Mw_i, const double mu,
                  const double s) const;
    
    /**
     * @brief class III vertex correction: vector form factor 'm' for lepton
     * @param[in] l lepton
     * @param[in] Mw_i the W mass
     * @param[in] mu the renormalization scale
     * @param[in] s the momentum squared
     * @return FI_Vml
     */ 
    complex FIII_Vml(const StandardModel::lepton l, const double Mw_i, const double mu,
                  const double s) const;
    
    /**
     * @brief class III vertex correction: vector form factor 'm' for quark
     * @param[in] q quark
     * @param[in] Mw_i the W mass
     * @param[in] mu the renormalization scale
     * @param[in] s the momentum squared
     * @return FI_Vmq
     */ 
    complex FIII_Vmq(const QCD::quark q, const double Mw_i, const double mu,
                  const double s) const;
    
    /**
     * @brief class III vertex correction: vector form factor 'n' for lepton
     * @param[in] l lepton
     * @param[in] Mw_i the W mass
     * @param[in] mu the renormalization scale
     * @param[in] s the momentum squared
     * @return FIII_Vnl
     */ 
    complex FIII_Vnl(const StandardModel::lepton l, const double Mw_i, const double mu,
                  const double s) const;
    
    /**
     * @brief class II vertex correction: vector form factor 'n' for quark
     * @param[in] q quark
     * @param[in] Mw_i the W mass
     * @param[in] mu the renormalization scale
     * @param[in] s the momentum squared
     * @return FI_Vnq
     */ 
    complex FIII_Vnq(const QCD::quark q, const double Mw_i, const double mu,
                  const double s) const;
    
    
    /**
     * @brief class III vertex correction: vector form factor 'o' for lepton
     * @param[in] l lepton
     * @param[in] Mw_i the W mass
     * @param[in] mu the renormalization scale
     * @param[in] s the momentum squared
     * @return FI_Vol
     */ 
    complex FIII_Vol(const StandardModel::lepton l, const double Mw_i, const double mu,
                  const double s) const;
    
    /**
     * @brief class II vertex correction: vector form factor 'o' for quark
     * @param[in] q quark
     * @param[in] Mw_i the W mass
     * @param[in] mu the renormalization scale
     * @param[in] s the momentum squared
     * @return FI_Voq
     */ 
    complex FIII_Voq(const QCD::quark q, const double Mw_i, const double mu,
                  const double s) const;
    
    
    /**
     * @brief class IV vertex correction: vector form factor 'p' for lepton
     * @param[in] l lepton
     * @param[in] Mw_i the W mass
     * @param[in] mu the renormalization scale
     * @param[in] s the momentum squared
     * @return FIV_Vpl
     */ 
    complex FIV_Vpl(const StandardModel::lepton l, const double Mw_i, const double mu,
                  const double s) const;
    
    /**
     * @brief class IV vertex correction: vector form factor 'p' for quark
     * @param[in] q quark
     * @param[in] Mw_i the W mass
     * @param[in] mu the renormalization scale
     * @param[in] s the momentum squared
     * @return FIV_Vpq
     */ 
    complex FIV_Vpq(const QCD::quark q, const double Mw_i, const double mu,
                  const double s) const;
    
    /**
     * @brief class V vertex correction: vector form factor 's' for lepton
     * @param[in] l lepton
     * @param[in] Mw_i the W mass
     * @param[in] mu the renormalization scale
     * @param[in] s the momentum squared
     * @return FV_Vsl
     */ 
    complex FV_Vsl(const StandardModel::lepton l, const double Mw_i, const double mu,
                  const double s) const;
    
    /**
     * @brief class V vertex correction: vector form factor 's' for quark
     * @param[in] q quark
     * @param[in] Mw_i the W mass
     * @param[in] mu the renormalization scale
     * @param[in] s the momentum squared
     * @return FV_Vsq
     */ 
    complex FV_Vsq(const QCD::quark q, const double Mw_i, const double mu,
                  const double s) const;
    
    /**
     * @brief class V vertex correction: vector form factor 't' for lepton
     * @param[in] l lepton
     * @param[in] Mw_i the W mass
     * @param[in] mu the renormalization scale
     * @param[in] s the momentum squared
     * @return FV_Vtl
     */ 
    complex FV_Vtl(const StandardModel::lepton l, const double Mw_i, const double mu,
                  const double s) const;
    
    /**
     * @brief class V vertex correction: vector form factor 't' for quark
     * @param[in] q quark
     * @param[in] Mw_i the W mass
     * @param[in] mu the renormalization scale
     * @param[in] s the momentum squared
     * @return FV_Vtq
     */ 
    complex FV_Vtq(const QCD::quark q, const double Mw_i, const double mu,
                  const double s) const;
    
    /**
     * @brief class VI vertex correction: vector form factor 'u' for lepton
     * @param[in] l lepton
     * @param[in] Mw_i the W mass
     * @param[in] mu the renormalization scale
     * @param[in] s the momentum squared
     * @return FVI_Vul
     */ 
    complex FVI_Vul(const StandardModel::lepton l, const double Mw_i, const double mu,
                  const double s) const;
    
    /**
     * @brief class VI vertex correction: vector form factor 'u' for quark
     * @param[in] q quark
     * @param[in] Mw_i the W mass
     * @param[in] mu the renormalization scale
     * @param[in] s the momentum squared
     * @return FVI_Vuq
     */ 
    complex FVI_Vuq(const QCD::quark q, const double Mw_i, const double mu,
                  const double s) const;
    
    /**
     * @brief class VI vertex correction: vector form factor 'v' for lepton
     * @param[in] l lepton
     * @param[in] Mw_i the W mass
     * @param[in] mu the renormalization scale
     * @param[in] s the momentum squared
     * @return FVI_Vvl
     */ 
    complex FVI_Vvl(const StandardModel::lepton l, const double Mw_i, const double mu,
                  const double s) const;
    
    /**
     * @brief class VI vertex correction: vector form factor 'v' for quark
     * @param[in] q quark
     * @param[in] Mw_i the W mass
     * @param[in] mu the renormalization scale
     * @param[in] s the momentum squared
     * @return FVI_Vvq
     */ 
    complex FVI_Vvq(const QCD::quark q, const double Mw_i, const double mu,
                  const double s) const;
    
    /**
     * @brief weak vector gammaff form factor 
     * @param[in] l lepton
     * @param[in] Mw_i the W mass
     * @param[in] mu the renormalization scale
     * @param[in] s the momentum squared
     * @return FVgammal_weak
     */ 
    complex FVgammal_weak(const StandardModel::lepton l, const double Mw_i, const double mu,
                  const double s) const;
    
    /**
     * @brief weak vector gammaff form factor 
     * @param[in] q quark
     * @param[in] Mw_i the W mass
     * @param[in] mu the renormalization scale
     * @param[in] s the momentum squared
     * @return FVgammaq_weak
     */ 
    complex FVgammaq_weak(const QCD::quark q, const double Mw_i, const double mu,
                  const double s) const;
    
    /**
     * @brief weak vector gammaff form factor 
     * @param[in] l lepton
     * @param[in] Mw_i the W mass
     * @param[in] mu the renormalization scale
     * @param[in] s the momentum squared
     * @return FVZl_weak
     */ 
    complex FVZl_weak(const StandardModel::lepton l, const double Mw_i, const double mu,
                  const double s) const;
    
    /**
     * @brief weak vector gammaff form factor 
     * @param[in] l lepton
     * @param[in] Mw_i the W mass
     * @param[in] mu the renormalization scale
     * @param[in] s the momentum squared
     * @return FVZq_weak
     */ 
    complex FVZq_weak(const QCD::quark q, const double Mw_i, const double mu,
                  const double s) const;
    
    
/////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////   MAGNETIC FORM FACTOR      /////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////
     
    
    
    /**
     * @brief class I vertex correction: vector form factor 'a' for lepton
     * @param[in] l lepton
     * @param[in] Mw_i the W mass
     * @param[in] mu the renormalization scale
     * @param[in] s the momentum squared
     * @return FI_Val
     */ 
    complex FI_Mal(const StandardModel::lepton l, const double Mw_i, const double mu,
                  const double s) const;
    
    /**
     * @brief class I vertex correction: vector form factor 'a' for quark
     * @param[in] q quark
     * @param[in] Mw_i the W mass
     * @param[in] mu the renormalization scale
     * @param[in] s the momentum squared
     * @return FI_Vaq
     */ 
    complex FI_Maq(const QCD::quark q, const double Mw_i, const double mu,
                  const double s) const;
    
    /**
     * @brief class I vertex correction: vector form factor 'b' for lepton
     * @param[in] l lepton
     * @param[in] Mw_i the W mass
     * @param[in] mu the renormalization scale
     * @param[in] s the momentum squared
     * @return FI_Val
     */ 
    complex FI_Mbl(const StandardModel::lepton l, const double Mw_i, const double mu,
                  const double s) const;
    
    /**
     * @brief class I vertex correction: vector form factor 'b' for quark
     * @param[in] q quark
     * @param[in] Mw_i the W mass
     * @param[in] mu the renormalization scale
     * @param[in] s the momentum squared
     * @return FI_Vaq
     */ 
    complex FI_Mbq(const QCD::quark q, const double Mw_i, const double mu,
                  const double s) const;
    
    /**
     * @brief class II vertex correction: vector form factor 'c' for lepton
     * @param[in] l lepton
     * @param[in] Mw_i the W mass
     * @param[in] mu the renormalization scale
     * @param[in] s the momentum squared
     * @return FI_Val
     */ 
    complex FII_Mcl(const StandardModel::lepton l, const double Mw_i, const double mu,
                  const double s) const;
    
    /**
     * @brief class II vertex correction: vector form factor 'c' for quark
     * @param[in] q quark
     * @param[in] Mw_i the W mass
     * @param[in] mu the renormalization scale
     * @param[in] s the momentum squared
     * @return FI_Vaq
     */ 
    complex FII_Mcq(const QCD::quark q, const double Mw_i, const double mu,
                  const double s) const;
    
    /**
     * @brief class III vertex correction: vector form factor 'd' for lepton
     * @param[in] l lepton
     * @param[in] Mw_i the W mass
     * @param[in] mu the renormalization scale
     * @param[in] s the momentum squared
     * @return FI_Val
     */ 
    complex FIII_Mdl(const StandardModel::lepton l, const double Mw_i, const double mu,
                  const double s) const;
    
    /**
     * @brief class II vertex correction: vector form factor 'd' for quark
     * @param[in] q quark
     * @param[in] Mw_i the W mass
     * @param[in] mu the renormalization scale
     * @param[in] s the momentum squared
     * @return FI_Vaq
     */ 
    complex FIII_Mdq(const QCD::quark q, const double Mw_i, const double mu,
                  const double s) const;
    
    /**
     * @brief class III vertex correction: vector form factor 'e' for lepton
     * @param[in] l lepton
     * @param[in] Mw_i the W mass
     * @param[in] mu the renormalization scale
     * @param[in] s the momentum squared
     * @return FIII_Vel
     */ 
    complex FIII_Mel(const StandardModel::lepton l, const double Mw_i, const double mu,
                  const double s) const;
    
    /**
     * @brief class II vertex correction: vector form factor 'd' for quark
     * @param[in] q quark
     * @param[in] Mw_i the W mass
     * @param[in] mu the renormalization scale
     * @param[in] s the momentum squared
     * @return FI_Vaq
     */ 
    complex FIII_Meq(const QCD::quark q, const double Mw_i, const double mu,
                  const double s) const;
    
    
    /**
     * @brief class III vertex correction: vector form factor 'd' for lepton
     * @param[in] l lepton
     * @param[in] Mw_i the W mass
     * @param[in] mu the renormalization scale
     * @param[in] s the momentum squared
     * @return FI_Val
     */ 
    complex FIII_Mfl(const StandardModel::lepton l, const double Mw_i, const double mu,
                  const double s) const;
    
    /**
     * @brief class II vertex correction: vector form factor 'd' for quark
     * @param[in] q quark
     * @param[in] Mw_i the W mass
     * @param[in] mu the renormalization scale
     * @param[in] s the momentum squared
     * @return FI_Vaq
     */ 
    complex FIII_Mfq(const QCD::quark q, const double Mw_i, const double mu,
                  const double s) const;
    
    /**
     * @brief class IV vertex correction: vector form factor 'g' for lepton
     * @param[in] l lepton
     * @param[in] Mw_i the W mass
     * @param[in] mu the renormalization scale
     * @param[in] s the momentum squared
     * @return FIV_Vgl
     */ 
    complex FIV_Mgl(const StandardModel::lepton l, const double Mw_i, const double mu,
                  const double s) const;
    
    /**
     * @brief class IV vertex correction: vector form factor 'g' for quark
     * @param[in] q quark
     * @param[in] Mw_i the W mass
     * @param[in] mu the renormalization scale
     * @param[in] s the momentum squared
     * @return FIV_Vgq
     */ 
    complex FIV_Mgq(const QCD::quark q, const double Mw_i, const double mu,
                  const double s) const;
    
    /**
     * @brief class V vertex correction: vector form factor 'h' for lepton
     * @param[in] l lepton
     * @param[in] Mw_i the W mass
     * @param[in] mu the renormalization scale
     * @param[in] s the momentum squared
     * @return FV_Vhl
     */ 
    complex FV_Mhl(const StandardModel::lepton l, const double Mw_i, const double mu,
                  const double s) const;
    
    /**
     * @brief class V vertex correction: vector form factor 'h' for quark
     * @param[in] q quark
     * @param[in] Mw_i the W mass
     * @param[in] mu the renormalization scale
     * @param[in] s the momentum squared
     * @return FV_Vhq
     */ 
    complex FV_Mhq(const QCD::quark q, const double Mw_i, const double mu,
                  const double s) const;
    
    /**
     * @brief class VI vertex correction: vector form factor 'i' for lepton
     * @param[in] l lepton
     * @param[in] Mw_i the W mass
     * @param[in] mu the renormalization scale
     * @param[in] s the momentum squared
     * @return FVI_Vil
     */ 
    complex FVI_Mil(const StandardModel::lepton l, const double Mw_i, const double mu,
                  const double s) const;
    
    /**
     * @brief class VI vertex correction: vector form factor 'i' for quark
     * @param[in] q quark
     * @param[in] Mw_i the W mass
     * @param[in] mu the renormalization scale
     * @param[in] s the momentum squared
     * @return FVI_Viq
     */ 
    complex FVI_Miq(const QCD::quark q, const double Mw_i, const double mu,
                  const double s) const;
    
    
    /**
     * @brief class I vertex correction: vector form factor 'j' for lepton
     * @param[in] l lepton
     * @param[in] Mw_i the W mass
     * @param[in] mu the renormalization scale
     * @param[in] s the momentum squared
     * @return FI_Vjl
     */ 
    complex FI_Mjl(const StandardModel::lepton l, const double Mw_i, const double mu,
                  const double s) const;
    
    /**
     * @brief class I vertex correction: vector form factor 'a' for quark
     * @param[in] q quark
     * @param[in] Mw_i the W mass
     * @param[in] mu the renormalization scale
     * @param[in] s the momentum squared
     * @return FI_Vaq
     */ 
    complex FI_Mjq(const QCD::quark q, const double Mw_i, const double mu,
                  const double s) const;
    
    /**
     * @brief class I vertex correction: vector form factor 'k' for lepton
     * @param[in] l lepton
     * @param[in] Mw_i the W mass
     * @param[in] mu the renormalization scale
     * @param[in] s the momentum squared
     * @return FI_Vkl
     */ 
    complex FI_Mkl(const StandardModel::lepton l, const double Mw_i, const double mu,
                  const double s) const;
    
    /**
     * @brief class I vertex correction: vector form factor 'k' for quark
     * @param[in] q quark
     * @param[in] Mw_i the W mass
     * @param[in] mu the renormalization scale
     * @param[in] s the momentum squared
     * @return FI_Vkq
     */ 
    complex FI_Mkq(const QCD::quark q, const double Mw_i, const double mu,
                  const double s) const;
    
    /**
     * @brief class II vertex correction: vector form factor 'l' for lepton
     * @param[in] l lepton
     * @param[in] Mw_i the W mass
     * @param[in] mu the renormalization scale
     * @param[in] s the momentum squared
     * @return FI_Val
     */ 
    complex FII_Mll(const StandardModel::lepton l, const double Mw_i, const double mu,
                  const double s) const;
    
    /**
     * @brief class II vertex correction: vector form factor 'l' for quark
     * @param[in] q quark
     * @param[in] Mw_i the W mass
     * @param[in] mu the renormalization scale
     * @param[in] s the momentum squared
     * @return FI_Vlq
     */ 
    complex FII_Mlq(const QCD::quark q, const double Mw_i, const double mu,
                  const double s) const;
    
    /**
     * @brief class III vertex correction: vector form factor 'm' for lepton
     * @param[in] l lepton
     * @param[in] Mw_i the W mass
     * @param[in] mu the renormalization scale
     * @param[in] s the momentum squared
     * @return FI_Vml
     */ 
    complex FIII_Mml(const StandardModel::lepton l, const double Mw_i, const double mu,
                  const double s) const;
    
    /**
     * @brief class III vertex correction: vector form factor 'm' for quark
     * @param[in] q quark
     * @param[in] Mw_i the W mass
     * @param[in] mu the renormalization scale
     * @param[in] s the momentum squared
     * @return FI_Vmq
     */ 
    complex FIII_Mmq(const QCD::quark q, const double Mw_i, const double mu,
                  const double s) const;
    
    /**
     * @brief class III vertex correction: vector form factor 'n' for lepton
     * @param[in] l lepton
     * @param[in] Mw_i the W mass
     * @param[in] mu the renormalization scale
     * @param[in] s the momentum squared
     * @return FIII_Vnl
     */ 
    complex FIII_Mnl(const StandardModel::lepton l, const double Mw_i, const double mu,
                  const double s) const;
    
    /**
     * @brief class II vertex correction: vector form factor 'n' for quark
     * @param[in] q quark
     * @param[in] Mw_i the W mass
     * @param[in] mu the renormalization scale
     * @param[in] s the momentum squared
     * @return FI_Vnq
     */ 
    complex FIII_Mnq(const QCD::quark q, const double Mw_i, const double mu,
                  const double s) const;
    
    
    /**
     * @brief class III vertex correction: vector form factor 'o' for lepton
     * @param[in] l lepton
     * @param[in] Mw_i the W mass
     * @param[in] mu the renormalization scale
     * @param[in] s the momentum squared
     * @return FI_Vol
     */ 
    complex FIII_Mol(const StandardModel::lepton l, const double Mw_i, const double mu,
                  const double s) const;
    
    /**
     * @brief class II vertex correction: vector form factor 'o' for quark
     * @param[in] q quark
     * @param[in] Mw_i the W mass
     * @param[in] mu the renormalization scale
     * @param[in] s the momentum squared
     * @return FI_Voq
     */ 
    complex FIII_Moq(const QCD::quark q, const double Mw_i, const double mu,
                  const double s) const;
    
    
    /**
     * @brief class IV vertex correction: vector form factor 'p' for lepton
     * @param[in] l lepton
     * @param[in] Mw_i the W mass
     * @param[in] mu the renormalization scale
     * @param[in] s the momentum squared
     * @return FIV_Vpl
     */ 
    complex FIV_Mpl(const StandardModel::lepton l, const double Mw_i, const double mu,
                  const double s) const;
    
    /**
     * @brief class IV vertex correction: vector form factor 'p' for quark
     * @param[in] q quark
     * @param[in] Mw_i the W mass
     * @param[in] mu the renormalization scale
     * @param[in] s the momentum squared
     * @return FIV_Vpq
     */ 
    complex FIV_Mpq(const QCD::quark q, const double Mw_i, const double mu,
                  const double s) const;
    
    /**
     * @brief class V vertex correction: vector form factor 's' for lepton
     * @param[in] l lepton
     * @param[in] Mw_i the W mass
     * @param[in] mu the renormalization scale
     * @param[in] s the momentum squared
     * @return FV_Vsl
     */ 
    complex FV_Msl(const StandardModel::lepton l, const double Mw_i, const double mu,
                  const double s) const;
    
    /**
     * @brief class V vertex correction: vector form factor 's' for quark
     * @param[in] q quark
     * @param[in] Mw_i the W mass
     * @param[in] mu the renormalization scale
     * @param[in] s the momentum squared
     * @return FV_Vsq
     */ 
    complex FV_Msq(const QCD::quark q, const double Mw_i, const double mu,
                  const double s) const;
    
    /**
     * @brief class V vertex correction: vector form factor 't' for lepton
     * @param[in] l lepton
     * @param[in] Mw_i the W mass
     * @param[in] mu the renormalization scale
     * @param[in] s the momentum squared
     * @return FV_Vtl
     */ 
    complex FV_Mtl(const StandardModel::lepton l, const double Mw_i, const double mu,
                  const double s) const;
    
    /**
     * @brief class V vertex correction: vector form factor 't' for quark
     * @param[in] q quark
     * @param[in] Mw_i the W mass
     * @param[in] mu the renormalization scale
     * @param[in] s the momentum squared
     * @return FV_Vtq
     */ 
    complex FV_Mtq(const QCD::quark q, const double Mw_i, const double mu,
                  const double s) const;
    
    /**
     * @brief class VI vertex correction: vector form factor 'u' for lepton
     * @param[in] l lepton
     * @param[in] Mw_i the W mass
     * @param[in] mu the renormalization scale
     * @param[in] s the momentum squared
     * @return FVI_Vul
     */ 
    complex FVI_Mul(const StandardModel::lepton l, const double Mw_i, const double mu,
                  const double s) const;
    
    /**
     * @brief class VI vertex correction: vector form factor 'u' for quark
     * @param[in] q quark
     * @param[in] Mw_i the W mass
     * @param[in] mu the renormalization scale
     * @param[in] s the momentum squared
     * @return FVI_Vuq
     */ 
    complex FVI_Muq(const QCD::quark q, const double Mw_i, const double mu,
                  const double s) const;
    
    /**
     * @brief class VI vertex correction: vector form factor 'v' for lepton
     * @param[in] l lepton
     * @param[in] Mw_i the W mass
     * @param[in] mu the renormalization scale
     * @param[in] s the momentum squared
     * @return FVI_Vvl
     */ 
    complex FVI_Mvl(const StandardModel::lepton l, const double Mw_i, const double mu,
                  const double s) const;
    
    /**
     * @brief class VI vertex correction: vector form factor 'v' for quark
     * @param[in] q quark
     * @param[in] Mw_i the W mass
     * @param[in] mu the renormalization scale
     * @param[in] s the momentum squared
     * @return FVI_Vvq
     */ 
    complex FVI_Mvq(const QCD::quark q, const double Mw_i, const double mu,
                  const double s) const;
    
    /**
     * @brief weak vector gammaff form factor 
     * @param[in] l lepton
     * @param[in] Mw_i the W mass
     * @param[in] mu the renormalization scale
     * @param[in] s the momentum squared
     * @return FVgammal_weak
     */ 
    complex FMgammal_weak(const StandardModel::lepton l, const double Mw_i, const double mu,
                  const double s) const;
    /**
     * @brief weak vector gammaff form factor 
     * @param[in] q quark
     * @param[in] Mw_i the W mass
     * @param[in] mu the renormalization scale
     * @param[in] s the momentum squared
     * @return FVgammaq_weak
     */ 
    complex FMgammaq_weak(const QCD::quark q, const double Mw_i, const double mu,
                  const double s) const;
    
    /**
     * @brief weak vector gammaff form factor 
     * @param[in] l lepton
     * @param[in] Mw_i the W mass
     * @param[in] mu the renormalization scale
     * @param[in] s the momentum squared
     * @return FVZl_weak
     */ 
    complex FMZl_weak(const StandardModel::lepton l, const double Mw_i, const double mu,
                  const double s) const;
    
    /**
     * @brief weak vector gammaff form factor 
     * @param[in] l lepton
     * @param[in] Mw_i the W mass
     * @param[in] mu the renormalization scale
     * @param[in] s the momentum squared
     * @return FVZq_weak
     */ 
    complex FMZq_weak(const QCD::quark q, const double Mw_i, const double mu,
                  const double s) const;
    
    
/////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////   AXIAL FORM FACTOR      /////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////
    
    
    
    /**
     * @brief class I vertex correction: axial form factor 'a' for lepton
     * @param[in] l lepton
     * @param[in] Mw_i the W mass
     * @param[in] mu the renormalization scale
     * @param[in] s the momentum squared
     * @return FI_Val
     */ 
    complex FI_Aal(const StandardModel::lepton l, const double Mw_i, const double mu,
                  const double s) const;
    
    /**
     * @brief class I vertex correction: vector form factor 'a' for quark
     * @param[in] q quark
     * @param[in] Mw_i the W mass
     * @param[in] mu the renormalization scale
     * @param[in] s the momentum squared
     * @return FI_Vaq
     */ 
    complex FI_Aaq(const QCD::quark q, const double Mw_i, const double mu,
                  const double s) const;
    
    /**
     * @brief class I vertex correction: vector form factor 'b' for lepton
     * @param[in] l lepton
     * @param[in] Mw_i the W mass
     * @param[in] mu the renormalization scale
     * @param[in] s the momentum squared
     * @return FI_Val
     */ 
    complex FI_Abl(const StandardModel::lepton l, const double Mw_i, const double mu,
                  const double s) const;
    
    /**
     * @brief class I vertex correction: vector form factor 'b' for quark
     * @param[in] q quark
     * @param[in] Mw_i the W mass
     * @param[in] mu the renormalization scale
     * @param[in] s the momentum squared
     * @return FI_Vaq
     */ 
    complex FI_Abq(const QCD::quark q, const double Mw_i, const double mu,
                  const double s) const;
    
    /**
     * @brief class II vertex correction: vector form factor 'c' for lepton
     * @param[in] l lepton
     * @param[in] Mw_i the W mass
     * @param[in] mu the renormalization scale
     * @param[in] s the momentum squared
     * @return FI_Val
     */ 
    complex FII_Acl(const StandardModel::lepton l, const double Mw_i, const double mu,
                  const double s) const;
    
    /**
     * @brief class II vertex correction: vector form factor 'c' for quark
     * @param[in] q quark
     * @param[in] Mw_i the W mass
     * @param[in] mu the renormalization scale
     * @param[in] s the momentum squared
     * @return FI_Vaq
     */ 
    complex FII_Acq(const QCD::quark q, const double Mw_i, const double mu,
                  const double s) const;
    
    /**
     * @brief class IV vertex correction: vector form factor 'g' for lepton
     * @param[in] l lepton
     * @param[in] Mw_i the W mass
     * @param[in] mu the renormalization scale
     * @param[in] s the momentum squared
     * @return FIV_Vgl
     */ 
    complex FIII_Afl(const StandardModel::lepton l, const double Mw_i, const double mu,
                  const double s) const;
    
    /**
     * @brief class IV vertex correction: vector form factor 'g' for quark
     * @param[in] q quark
     * @param[in] Mw_i the W mass
     * @param[in] mu the renormalization scale
     * @param[in] s the momentum squared
     * @return FIV_Vgq
     */ 
    complex FIII_Afq(const QCD::quark q, const double Mw_i, const double mu,
                  const double s) const;
    
    /**
     * @brief class IV vertex correction: vector form factor 'g' for lepton
     * @param[in] l lepton
     * @param[in] Mw_i the W mass
     * @param[in] mu the renormalization scale
     * @param[in] s the momentum squared
     * @return FIV_Vgl
     */ 
    complex FIV_Agl(const StandardModel::lepton l, const double Mw_i, const double mu,
                  const double s) const;
    
    /**
     * @brief class IV vertex correction: vector form factor 'g' for quark
     * @param[in] q quark
     * @param[in] Mw_i the W mass
     * @param[in] mu the renormalization scale
     * @param[in] s the momentum squared
     * @return FIV_Vgq
     */ 
    complex FIV_Agq(const QCD::quark q, const double Mw_i, const double mu,
                  const double s) const;
    
    /**
     * @brief class V vertex correction: vector form factor 'h' for lepton
     * @param[in] l lepton
     * @param[in] Mw_i the W mass
     * @param[in] mu the renormalization scale
     * @param[in] s the momentum squared
     * @return FV_Vhl
     */ 
    complex FV_Ahl(const StandardModel::lepton l, const double Mw_i, const double mu,
                  const double s) const;
    
    /**
     * @brief class V vertex correction: vector form factor 'h' for quark
     * @param[in] q quark
     * @param[in] Mw_i the W mass
     * @param[in] mu the renormalization scale
     * @param[in] s the momentum squared
     * @return FV_Vhq
     */ 
    complex FV_Ahq(const QCD::quark q, const double Mw_i, const double mu,
                  const double s) const;
    
    /**
     * @brief class VI vertex correction: vector form factor 'i' for lepton
     * @param[in] l lepton
     * @param[in] Mw_i the W mass
     * @param[in] mu the renormalization scale
     * @param[in] s the momentum squared
     * @return FVI_Vil
     */ 
    complex FVI_Ail(const StandardModel::lepton l, const double Mw_i, const double mu,
                  const double s) const;
    
    /**
     * @brief class VI vertex correction: vector form factor 'i' for quark
     * @param[in] q quark
     * @param[in] Mw_i the W mass
     * @param[in] mu the renormalization scale
     * @param[in] s the momentum squared
     * @return FVI_Viq
     */ 
    complex FVI_Aiq(const QCD::quark q, const double Mw_i, const double mu,
                  const double s) const;
    
    
    /**
     * @brief class I vertex correction: vector form factor 'j' for lepton
     * @param[in] l lepton
     * @param[in] Mw_i the W mass
     * @param[in] mu the renormalization scale
     * @param[in] s the momentum squared
     * @return FI_Vjl
     */ 
    complex FI_Ajl(const StandardModel::lepton l, const double Mw_i, const double mu,
                  const double s) const;
    
    /**
     * @brief class I vertex correction: vector form factor 'a' for quark
     * @param[in] q quark
     * @param[in] Mw_i the W mass
     * @param[in] mu the renormalization scale
     * @param[in] s the momentum squared
     * @return FI_Vaq
     */ 
    complex FI_Ajq(const QCD::quark q, const double Mw_i, const double mu,
                  const double s) const;
    
    /**
     * @brief class I vertex correction: vector form factor 'k' for lepton
     * @param[in] l lepton
     * @param[in] Mw_i the W mass
     * @param[in] mu the renormalization scale
     * @param[in] s the momentum squared
     * @return FI_Vkl
     */ 
    complex FI_Akl(const StandardModel::lepton l, const double Mw_i, const double mu,
                  const double s) const;
    
    /**
     * @brief class I vertex correction: vector form factor 'k' for quark
     * @param[in] q quark
     * @param[in] Mw_i the W mass
     * @param[in] mu the renormalization scale
     * @param[in] s the momentum squared
     * @return FI_Vkq
     */ 
    complex FI_Akq(const QCD::quark q, const double Mw_i, const double mu,
                  const double s) const;
    
    /**
     * @brief class II vertex correction: vector form factor 'l' for lepton
     * @param[in] l lepton
     * @param[in] Mw_i the W mass
     * @param[in] mu the renormalization scale
     * @param[in] s the momentum squared
     * @return FI_Val
     */ 
    complex FII_All(const StandardModel::lepton l, const double Mw_i, const double mu,
                  const double s) const;
    
    /**
     * @brief class II vertex correction: vector form factor 'l' for quark
     * @param[in] q quark
     * @param[in] Mw_i the W mass
     * @param[in] mu the renormalization scale
     * @param[in] s the momentum squared
     * @return FI_Vlq
     */ 
    complex FII_Alq(const QCD::quark q, const double Mw_i, const double mu,
                  const double s) const;
    
    /**
     * @brief class III vertex correction: vector form factor 'm' for lepton
     * @param[in] l lepton
     * @param[in] Mw_i the W mass
     * @param[in] mu the renormalization scale
     * @param[in] s the momentum squared
     * @return FI_Vml
     */ 
    complex FIII_Aml(const StandardModel::lepton l, const double Mw_i, const double mu,
                  const double s) const;
    
    /**
     * @brief class III vertex correction: vector form factor 'm' for quark
     * @param[in] q quark
     * @param[in] Mw_i the W mass
     * @param[in] mu the renormalization scale
     * @param[in] s the momentum squared
     * @return FI_Vmq
     */ 
    complex FIII_Amq(const QCD::quark q, const double Mw_i, const double mu,
                  const double s) const;
    
    /**
     * @brief class III vertex correction: vector form factor 'n' for lepton
     * @param[in] l lepton
     * @param[in] Mw_i the W mass
     * @param[in] mu the renormalization scale
     * @param[in] s the momentum squared
     * @return FIII_Vnl
     */ 
    complex FIII_Anl(const StandardModel::lepton l, const double Mw_i, const double mu,
                  const double s) const;
    
    /**
     * @brief class II vertex correction: vector form factor 'n' for quark
     * @param[in] q quark
     * @param[in] Mw_i the W mass
     * @param[in] mu the renormalization scale
     * @param[in] s the momentum squared
     * @return FI_Vnq
     */ 
    complex FIII_Anq(const QCD::quark q, const double Mw_i, const double mu,
                  const double s) const;
    
    
    /**
     * @brief class III vertex correction: vector form factor 'o' for lepton
     * @param[in] l lepton
     * @param[in] Mw_i the W mass
     * @param[in] mu the renormalization scale
     * @param[in] s the momentum squared
     * @return FI_Vol
     */ 
    complex FIII_Aol(const StandardModel::lepton l, const double Mw_i, const double mu,
                  const double s) const;
    
    /**
     * @brief class II vertex correction: vector form factor 'o' for quark
     * @param[in] q quark
     * @param[in] Mw_i the W mass
     * @param[in] mu the renormalization scale
     * @param[in] s the momentum squared
     * @return FI_Voq
     */ 
    complex FIII_Aoq(const QCD::quark q, const double Mw_i, const double mu,
                  const double s) const;
    
    
    /**
     * @brief class IV vertex correction: vector form factor 'p' for lepton
     * @param[in] l lepton
     * @param[in] Mw_i the W mass
     * @param[in] mu the renormalization scale
     * @param[in] s the momentum squared
     * @return FIV_Vpl
     */ 
    complex FIV_Apl(const StandardModel::lepton l, const double Mw_i, const double mu,
                  const double s) const;
    
    /**
     * @brief class IV vertex correction: vector form factor 'p' for quark
     * @param[in] q quark
     * @param[in] Mw_i the W mass
     * @param[in] mu the renormalization scale
     * @param[in] s the momentum squared
     * @return FIV_Vpq
     */ 
    complex FIV_Apq(const QCD::quark q, const double Mw_i, const double mu,
                  const double s) const;
    
    /**
     * @brief class IV vertex correction: vector form factor 'p' for lepton
     * @param[in] l lepton
     * @param[in] Mw_i the W mass
     * @param[in] mu the renormalization scale
     * @param[in] s the momentum squared
     * @return FIV_Vpl
     */ 
    complex FIV_Aql(const StandardModel::lepton l, const double Mw_i, const double mu,
                  const double s) const;
    
    /**
     * @brief class IV vertex correction: vector form factor 'p' for quark
     * @param[in] q quark
     * @param[in] Mw_i the W mass
     * @param[in] mu the renormalization scale
     * @param[in] s the momentum squared
     * @return FIV_Vpq
     */ 
    complex FIV_Aqq(const QCD::quark q, const double Mw_i, const double mu,
                  const double s) const;
    /**
     * @brief class IV vertex correction: vector form factor 'p' for lepton
     * @param[in] l lepton
     * @param[in] Mw_i the W mass
     * @param[in] mu the renormalization scale
     * @param[in] s the momentum squared
     * @return FIV_Vpl
     */ 
    complex FIV_Arl(const StandardModel::lepton l, const double Mw_i, const double mu,
                  const double s) const;
    
    /**
     * @brief class IV vertex correction: vector form factor 'p' for quark
     * @param[in] q quark
     * @param[in] Mw_i the W mass
     * @param[in] mu the renormalization scale
     * @param[in] s the momentum squared
     * @return FIV_Vpq
     */ 
    complex FIV_Arq(const QCD::quark q, const double Mw_i, const double mu,
                  const double s) const;
    
    /**
     * @brief class V vertex correction: vector form factor 's' for lepton
     * @param[in] l lepton
     * @param[in] Mw_i the W mass
     * @param[in] mu the renormalization scale
     * @param[in] s the momentum squared
     * @return FV_Vsl
     */ 
    complex FV_Asl(const StandardModel::lepton l, const double Mw_i, const double mu,
                  const double s) const;
    
    /**
     * @brief class V vertex correction: vector form factor 's' for quark
     * @param[in] q quark
     * @param[in] Mw_i the W mass
     * @param[in] mu the renormalization scale
     * @param[in] s the momentum squared
     * @return FV_Vsq
     */ 
    complex FV_Asq(const QCD::quark q, const double Mw_i, const double mu,
                  const double s) const;
    
    /**
     * @brief class V vertex correction: vector form factor 't' for lepton
     * @param[in] l lepton
     * @param[in] Mw_i the W mass
     * @param[in] mu the renormalization scale
     * @param[in] s the momentum squared
     * @return FV_Vtl
     */ 
    complex FV_Atl(const StandardModel::lepton l, const double Mw_i, const double mu,
                  const double s) const;
    
    /**
     * @brief class V vertex correction: vector form factor 't' for quark
     * @param[in] q quark
     * @param[in] Mw_i the W mass
     * @param[in] mu the renormalization scale
     * @param[in] s the momentum squared
     * @return FV_Vtq
     */ 
    complex FV_Atq(const QCD::quark q, const double Mw_i, const double mu,
                  const double s) const;
    
    /**
     * @brief class VI vertex correction: vector form factor 'u' for lepton
     * @param[in] l lepton
     * @param[in] Mw_i the W mass
     * @param[in] mu the renormalization scale
     * @param[in] s the momentum squared
     * @return FVI_Vul
     */ 
    complex FVI_Aul(const StandardModel::lepton l, const double Mw_i, const double mu,
                  const double s) const;
    
    /**
     * @brief class VI vertex correction: vector form factor 'u' for quark
     * @param[in] q quark
     * @param[in] Mw_i the W mass
     * @param[in] mu the renormalization scale
     * @param[in] s the momentum squared
     * @return FVI_Vuq
     */ 
    complex FVI_Auq(const QCD::quark q, const double Mw_i, const double mu,
                  const double s) const;
    
    /**
     * @brief class VI vertex correction: vector form factor 'v' for lepton
     * @param[in] l lepton
     * @param[in] Mw_i the W mass
     * @param[in] mu the renormalization scale
     * @param[in] s the momentum squared
     * @return FVI_Vvl
     */ 
    complex FVI_Avl(const StandardModel::lepton l, const double Mw_i, const double mu,
                  const double s) const;
    
    /**
     * @brief class VI vertex correction: vector form factor 'v' for quark
     * @param[in] q quark
     * @param[in] Mw_i the W mass
     * @param[in] mu the renormalization scale
     * @param[in] s the momentum squared
     * @return FVI_Vvq
     */ 
    complex FVI_Avq(const QCD::quark q, const double Mw_i, const double mu,
                  const double s) const;
    
    /**
     * @brief weak vector gammaff form factor 
     * @param[in] l lepton
     * @param[in] Mw_i the W mass
     * @param[in] mu the renormalization scale
     * @param[in] s the momentum squared
     * @return FVgammal_weak
     */ 
    complex FAgammal_weak(const StandardModel::lepton l, const double Mw_i, const double mu,
                  const double s) const;
    
    /**
     * @brief weak vector gammaff form factor 
     * @param[in] q quark
     * @param[in] Mw_i the W mass
     * @param[in] mu the renormalization scale
     * @param[in] s the momentum squared
     * @return FVgammaq_weak
     */ 
    complex FAgammaq_weak(const QCD::quark q, const double Mw_i, const double mu,
                  const double s) const;
    
    /**
     * @brief weak vector gammaff form factor 
     * @param[in] l lepton
     * @param[in] Mw_i the W mass
     * @param[in] mu the renormalization scale
     * @param[in] s the momentum squared
     * @return FVZl_weak
     */ 
    complex FAZl_weak(const StandardModel::lepton l, const double Mw_i, const double mu,
                  const double s) const;
    
    /**
     * @brief weak vector gammaff form factor 
     * @param[in] l lepton
     * @param[in] Mw_i the W mass
     * @param[in] mu the renormalization scale
     * @param[in] s the momentum squared
     * @return FVZq_weak
     */ 
    complex FAZq_weak(const QCD::quark q, const double Mw_i, const double mu,
                  const double s) const;
    
    
    /**
     * @brief coefficient for weak vertex corrections ffgamma
     * @param[in] l lepton
     * @param[in] rho  the elicity
     * @param[in] s  the momentum squared
     * @param[in] Mw_i  the W mass
     * @return Chi_gammaZ
     */
    complex E1(const double mu, const double k, const double s, const double Mw_i,
                      const double W, const double X, const double Y) const;
    
    /**
     * @brief coefficient for weak vertex corrections ffZ
     * @param[in] mu renormalization scale
     * @param[in] q quark
     * @param[in] rho  the elicity
     * @param[in] s  the momentum squared
     * @param[in] Mw_i  the W mass
     * @return Chi_gammaZ
     */
    complex E2(const double mu, const double k, const double s, const double Mw_i,
                      const double W, const double X, const double Y) const;
    
    
    /**
     * @brief coefficient for weak vertex corrections ffgamma
     * @param[in] l lepton
     * @param[in] rho  the elicity
     * @param[in] s  the momentum squared
     * @param[in] Mw_i  the W mass
     * @return Chi_gammaZ
     */
    complex F1_l(const double mu,const double rho, const double s,
                              const double Mw_i,const StandardModel::lepton l) const;
    
    /**
     * @brief coefficient for weak vertex corrections ffZ
     * @param[in] mu renormalization scale
     * @param[in] q quark
     * @param[in] rho  the elicity
     * @param[in] s  the momentum squared
     * @param[in] Mw_i  the W mass
     * @return Chi_gammaZ
     */
    complex F2_l(const double mu, const double rho, const double s, const double Mw_i
                  ,const StandardModel::lepton l) const;
    
    /**
     * @brief coefficient for weak vertex corrections ffgamma
     * @param[in] l lepton
     * @param[in] rho  the elicity
     * @param[in] s  the momentum squared
     * @param[in] Mw_i  the W mass
     * @return Chi_gammaZ
     */
    complex F1_q(const double mu,const double rho, const double s,
                              const double Mw_i,const QCD::quark q) const;
    
    /**
     * @brief coefficient for weak vertex corrections ffZ
     * @param[in] mu renormalization scale
     * @param[in] q quark
     * @param[in] rho  the elicity
     * @param[in] s  the momentum squared
     * @param[in] Mw_i  the W mass
     * @return Chi_gammaZ
     */
    complex F2_q(const double mu,const double rho, const double s,
                              const double Mw_i,const QCD::quark q) const;
    
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////        WEAK BOX CONTRIBUTION        /////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////

/**
     * @brief coefficient for weak vertex corrections ffZ
     * @param[in] mu renormalization scale
     * @param[in] q quark
     * @param[in] rho  the elicity
     * @param[in] s  the momentum squared
     * @param[in] Mw_i  the W mass
     * @return Chi_gammaZ
     */
    complex A_CC(const double k,const double Mw_i) const;

/**
     * @brief coefficient for weak vertex corrections ffZ
     * @param[in] mu renormalization scale
     * @param[in] q quark
     * @param[in] rho  the elicity
     * @param[in] s  the momentum squared
     * @param[in] Mw_i  the W mass
     * @return Chi_gammaZ
     */
    complex B_CCq(const double s,const double theta,const QCD::quark q, const double Mw_i) const;
    
    
    /**
     * @brief coefficient for weak vertex corrections ffZ
     * @param[in] mu renormalization scale
     * @param[in] q quark
     * @param[in] rho  the elicity
     * @param[in] s  the momentum squared
     * @param[in] Mw_i  the W mass
     * @return Chi_gammaZ
     */
    complex B_CCl(const double s,const double theta,const StandardModel::lepton l, const double Mw_i) const;
    
    complex C_CCq(const double mu,const double s,const double theta,const QCD::quark q, const double Mw_i) const;
      
    complex C_CCl(const double mu,const double s,const double theta,const StandardModel::lepton l, const double Mw_i) const;
    
    
    complex D_CCq(const double mu,const double s,const double theta,const QCD::quark q, const double Mw_i) const;
    
    complex D_CCl(const double mu,const double s,const double theta,const StandardModel::lepton l, const double Mw_i) const;
    
    
    complex E_CCq(const double mu,const double s,const double theta,const QCD::quark q, const double Mw_i) const;
    
    complex E_CCl(const double mu,const double s,const double theta,const StandardModel::lepton l, const double Mw_i) const;
    
    complex A1_NCq(const double mu,const double s,const double theta,const QCD::quark q, const double Mw_i, const double rho, const double k) const;
    
    complex A1_NCl(const double mu,const double s,const double theta,const StandardModel::lepton l, const double Mw_i, const double rho, const double k) const;
    
    complex A2_NCq(const double mu,const double s,const double theta,const QCD::quark q, const double Mw_i, const double rho, const double k) const;
    
    complex A2_NCl(const double mu,const double s,const double theta,const StandardModel::lepton l, const double Mw_i, const double rho, const double k) const;
    
    complex A3_NCq(const double mu,const double s,const double theta,const QCD::quark q, const double Mw_i, const double rho, const double k) const;
    
    complex A3_NCl(const double mu,const double s,const double theta,const StandardModel::lepton l, const double Mw_i, const double rho, const double k) const;
    
    complex A4_NCq(const double mu,const double s,const double theta,const QCD::quark q, const double Mw_i, const double rho, const double k) const;
    
    complex A4_NCl(const double mu,const double s,const double theta,const StandardModel::lepton l, const double Mw_i, const double rho, const double k) const;
    



    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    

/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////        QED VERTEX AND BOXCONTRIBUTION        ////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////


    /**
     * @brief coefficient for eegamma-eeZ QED vertex correction
     * @param[in] s the momentum squared 
     * @return Lambda1
     */
    complex Lambda1(const double s) const;   
    
    /**
     * @brief coefficient for eegamma-eeZ QED vertex correction
     * @param[in] s  the momentum squared
     * @return E3
     */
    complex E3(const double s) const;
    
    /**
     * @brief coefficient for eegamma-eeZ QED vertex correction
     * @param[in] s the momentum squared 
     * @param[in] q quark
     * @return Lambda1
     */
    complex Lambdaq(const double s, const QCD::quark q, const double mu) const;  
    /**
     * @brief coefficient for eegamma-eeZ QED vertex correction
     * @param[in] s the momentum squared 
     * @param[in] l lepton
     * @return Lambda1
     */
    complex Lambdal(const double s,const StandardModel::lepton l, const double mu) const;  
    
    /**
     * @brief coefficient for eegamma-eeZ QED vertex correction
     * @param[in] s the momentum squared 
     * @param[in] q quark
     * @return Lambda1
     */
    complex LambdaMq(const double s,const QCD::quark q, const double mu) const;  
    /**
     * @brief coefficient for eegamma-eeZ QED vertex correction
     * @param[in] s the momentum squared 
     * @param[in] l lepton
     * @return Lambda1
     */
    complex LambdaMl(const double s,const StandardModel::lepton l, const double mu) const;  
    
    /**
     * @brief coefficient for eegamma-eeZ QED vertex correction
     * @param[in] s  the momentum squared
     * @param[in] l lepton
     * @return E3
     */
    complex E4l(const double s,const StandardModel::lepton l
               , const double mu, const double Mw_i,
                      const double W, const double X, const double Y) const;
    
    /**
     * @brief coefficient for eegamma-eeZ QED vertex correction
     * @param[in] s  the momentum squared
     * @param[in] q quark
     * @return E3
     */
    complex E4q(const double s,const QCD::quark q 
                 , const double mu, const double Mw_i,
                      const double W, const double X, const double Y) const;
    
    
    /**
     * @brief coefficient for eegamma-eeZ QED vertex correction
     * @param[in] s  the momentum squared
     * @return E3
     */
    complex E5(const double s, const double k, const double mu, const double Mw_i,
                      const double W, const double X, const double Y) const;
    
    /**
     * @brief coefficient for eegamma-eeZ QED vertex correction
     * @param[in] s  the momentum squared
     * @return E3
     */
    complex F5l(const double s, const double rho, const double mu, const double Mw_i,const StandardModel::lepton l) const;
   
    /**
     * @brief coefficient for eegamma-eeZ QED vertex correction
     * @param[in] s  the momentum squared
     * @return E3
     */
    complex F5q(const double s, const double rho, const double mu, const double Mw_i,const QCD::quark q) const;
    
    /**
     * @brief coefficient for eegamma-eeZ QED vertex correction
     * @param[in] s  the momentum squared
     * @return E3
     */
    complex G5l(const double s, const double mu, const double Mw_i,const StandardModel::lepton l) const;
    
    /**
     * @brief coefficient for eegamma-eeZ QED vertex correction
     * @param[in] s  the momentum squared
     * @return E3
     */
    complex G5q(const double s, const double mu, const double Mw_i,const QCD::quark q) const;
    
    
    /**
     * @brief function G for gammagamma QED boxes corrections
     * @param[in] s  the momentum squared
     * @return E3
     */
    complex Gfunc(const double s, const double t) const;
    
    
    /**
     * @brief function G for gammagamma QED boxes corrections
     * @param[in] s  the momentum squared
     * @return E3
     */
    complex A_gammagammaq(const double s, const QCD::quark q, const double theta) const;
    
    /**
     * @brief function G for gammagamma QED boxes corrections
     * @param[in] s  the momentum squared
     * @return E3
     */
    complex A_gammagammal(const double s, const StandardModel::lepton l, const double theta) const;
    
   
    /**
     * @brief function G for gammagamma QED boxes corrections
     * @param[in] s  the momentum squared
     * @return E3
     */
    complex V_gammagammaq(const double s, const QCD::quark q, const double theta) const;
    
    /**
     * @brief function G for gammagamma QED boxes corrections
     * @param[in] s  the momentum squared
     * @return E3
     */
    complex V_gammagammal(const double s, const StandardModel::lepton l, const double theta) const;
    
    
    /**
     * @brief function G for gammagamma QED boxes corrections
     * @param[in] s  the momentum squared
     * @return E3
     */
    complex Afunc(const double s, const double t,const double GammaZ) const;
    
    /**
     * @brief function G for gammagamma QED boxes corrections
     * @param[in] s  the momentum squared
     * @return E3
     */
    complex A_gammaZq(const double s,const QCD::quark q, const double theta,const double GammaZ) const;
    /**
     * @brief function G for gammagamma QED boxes corrections
     * @param[in] s  the momentum squared
     * @return E3
     */
    complex A_gammaZl(const double s,const StandardModel::lepton l, const double theta,const double GammaZ) const;
   
    /**
     * @brief function G for gammagamma QED boxes corrections
     * @param[in] s  the momentum squared
     * @return E3
     */
    complex V_gammaZq(const double s,const QCD::quark q, const double theta,const double GammaZ) const;
    
    /**
     * @brief function G for gammagamma QED boxes corrections
     * @param[in] s  the momentum squared
     * @return E3
     */
    complex V_gammaZl(const double s,const StandardModel::lepton l, const double theta,const double GammaZ) const;
    
    
    /* @brief coefficient for eegamma-eeZ QED vertex correction
     * @param[in] s  the momentum squared
     * @param[in] l lepton
     * @return E3
     */
    complex E6l(const double s,const StandardModel::lepton l
               , const double mu, const double Mw_i,
                      const double W, const double X, const double Y) const;
    
    /**
     * @brief coefficient for eegamma-eeZ QED vertex correction
     * @param[in] s  the momentum squared
     * @param[in] q quark
     * @return E3
     */
    complex E6q(const double s,const QCD::quark q 
                 , const double mu, const double Mw_i,
                      const double W, const double X, const double Y) const;
    
    /**
     * @brief coefficient for eegamma-eeZ QED vertex correction
     * @param[in] s  the momentum squared
     * @param[in] q quark
     * @return E3
     */
    complex F6rhoq(const double s,const double rho, const double k,const QCD::quark q, const double theta) const;
    
    /**
     * @brief coefficient for eegamma-eeZ QED vertex correction
     * @param[in] s  the momentum squared
     * @param[in] q quark
     * @return E3
     */
    complex F6rhol(const double s,const double rho, const double k,
                   const StandardModel::lepton l, const double theta) const;
    
    
    /**
     * @brief coefficient for eegamma-eeZ QED vertex correction
     * @param[in] s  the momentum squared
     * @param[in] q quark
     * @return E3
     */
    complex F7rhol(const double s,const double rho, const double k,const double Mw_i,
                   const StandardModel::lepton l, const double theta,const double GammaZ) const;
    
    
    
    /**
     * @brief coefficient for eegamma-eeZ QED vertex correction
     * @param[in] s  the momentum squared
     * @param[in] q quark
     * @return E3
     */
    complex F7rhoq(const double s,const double rho, const double k,const double Mw_i,
                   const QCD::quark q, const double theta,const double GammaZ) const;
    
    
    
/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////       SOFT PHOTON APPROXIMATION OF BREHMSTRAHLUNG CROSS SECTION     ///////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////

    
    
    
    double delta(const double s) const;
    
    
    
    double xf(const double s, const double mf) const;
    
    
    
    
    double Be(const double s) const;
    
    
    
    
    double Bf(const double s, const double mf) const;
    
    
    
    
    double Bint(const double theta, const double mf, const double s) const;
    
    
    double gammaIR_q(const double s, const QCD::quark q, const double theta) const;
    
    
    double gammaIR_l(const double s, const StandardModel::lepton l, const double theta) const;
    
    double gammadelta_q(const double s, const QCD::quark q, const double theta) const;
    
    
    double gammadelta_l(const double s, const StandardModel::lepton l, const double theta) const;
    
    complex gammadeltaINT_q(const double s, const QCD::quark q, const double theta,const double GammaZ) const;
    
    
    complex gammadeltaINT_l(const double s, const StandardModel::lepton l, const double theta,const double GammaZ) const;
    
    double gammadeltaRES_q(const double s, const QCD::quark q, const double theta,const double GammaZ) const;
    
    
    double gammadeltaRES_l(const double s, const StandardModel::lepton l, const double theta,const double GammaZ) const;
    
    double gammatail_q(const double s, const QCD::quark q,const double GammaZ) const;
    
    
    double gammatail_l(const double s, const StandardModel::lepton l,const double GammaZ) const;
    
    double gammafin_q(const double s, const QCD::quark q, const double theta) const;
    
    
    double gammafin_l(const double s, const StandardModel::lepton l,const double theta) const;
    
    
    double deltagammagamma_softq(const double s, const QCD::quark q, const double theta) const;
    
    double deltagammagamma_softl(const double s, const StandardModel::lepton l, const double theta) const;
    
    complex deltagammaZ_softq(const double s, const QCD::quark q, const double theta,const double GammaZ) const;
    
    complex deltagammaZ_softl(const double s, const StandardModel::lepton l, const double theta,const double GammaZ) const;
    
    double deltaZZ_softq(const double s, const QCD::quark q, const double theta,const double GammaZ) const;

    double deltaZZ_softl(const double s, const StandardModel::lepton l, const double theta,const double GammaZ) const;
    
    
    
/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////       CHIRALITY AMPLITUDE     //////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////
 
    
    
    
    /**
     * @brief interference of the basic chirality amplitude
     * @param[in] s  the momentum squared
     * @param[in] rho final elicity
     * @param[in] rhoprime final elicity
     * @param[in] k initial elicity
     * @param[in] l lepton
     * @return M1rhok_M1rhopk_l
     */
    complex M1rhok_M1rhopk_l(const double s,const double rho,const double rhoprime,
                             const double k,const StandardModel::lepton l,const double theta) const;
    
    
    /**
     * @brief interference of the basic chirality amplitude
     * @param[in] s  the momentum squared
     * @param[in] rho final elicity
     * @param[in] rhoprime final elicity
     * @param[in] k initial elicity
     * @param[in] q quark
     * @return M1rhok_M1rhopk_l
     */
    complex M1rhok_M1rhopk_q(const double s,const double rho,const double rhoprime,
                             const double k,const QCD::quark q,const double theta) const;
    
    /**
     * @brief interference of the basic chirality amplitude
     * @param[in] s  the momentum squared
     * @param[in] rho final elicity
     * @param[in] rhoprime final elicity
     * @param[in] k initial elicity
     * @param[in] l lepton
     * @return M2rhok_M1rhopk_l
     */
    complex M2rhok_M1rhopk_l(const double s,const double rho,const double rhoprime,
                             const double k,const StandardModel::lepton l,const double theta) const;
    
    
    /**
     * @brief interference of the basic chirality amplitude
     * @param[in] s  the momentum squared
     * @param[in] rho final elicity
     * @param[in] rhoprime final elicity
     * @param[in] k initial elicity
     * @param[in] q quark
     * @return M2rhok_M1rhopk_l
     */
    complex M2rhok_M1rhopk_q(const double s,const double rho,const double rhoprime,
                             const double k,const QCD::quark q,const double theta) const;
    
    /**
     * @brief interference of the basic chirality amplitude
     * @param[in] s  the momentum squared
     * @param[in] rho final elicity
     * @param[in] rhoprime final elicity
     * @param[in] k initial elicity
     * @param[in] l lepton
     * @return M3rhok_M1rhopk_l
     */
    complex M3rhok_M1rhopk_l(const double s,const double rho,const double rhoprime,
                             const double k,const StandardModel::lepton l,const double theta) const;
    
    
    /**
     * @brief interference of the basic chirality amplitude
     * @param[in] s  the momentum squared
     * @param[in] rho final elicity
     * @param[in] rhoprime final elicity
     * @param[in] k initial elicity
     * @param[in] q quark
     * @return M3rhok_M1rhopk_l
     */
    complex M3rhok_M1rhopk_q(const double s,const double rho,const double rhoprime,
                             const double k,const QCD::quark q,const double theta) const;
    
    
    /**
     * @brief interference of the basic chirality amplitude
     * @param[in] s  the momentum squared
     * @param[in] rho final elicity
     * @param[in] rhoprime final elicity
     * @param[in] k initial elicity
     * @param[in] l lepton
     * @return M4rhok_M1rhopk_l
     */
    complex M4rhok_M1rhopk_l(const double s,const double rho,const double rhoprime,
                             const double k,const StandardModel::lepton l,const double theta) const;
    
    
    /**
     * @brief interference of the basic chirality amplitude
     * @param[in] s  the momentum squared
     * @param[in] rho final elicity
     * @param[in] rhoprime final elicity
     * @param[in] k initial elicity
     * @param[in] q quark
     * @return M4rhok_M1rhopk_l
     */
    complex M4rhok_M1rhopk_q(const double s,const double rho,const double rhoprime,
                             const double k,const QCD::quark q,const double theta) const;
    
    
    
    
    complex ATOTq(const double mu,const QCD::quark q, const double rho,
                  const double k, const double s, const double Mw_i,const double theta,
                               const double W,const double X,const double Y,const double GammaZ) const;
    
    complex ATOTl(const double mu,const StandardModel::lepton l, const double rho,
                  const double k, const double s, const double Mw_i,const double theta,
                               const double W,const double X,const double Y,const double GammaZ) const;
    
    complex BTOTq(const double mu,const QCD::quark q, const double rho,
                               const double k, const double s, const double Mw_i,const double theta,
                               const double W,const double X,const double Y) const;
    
    complex BTOTl(const double mu,const StandardModel::lepton l, const double rho,
                  const double k, const double s, const double Mw_i,const double theta,
                               const double W,const double X,const double Y) const;
    
    
    //Total amplitude squared for q
    double MTOTq_sq(const QCD::quark q, const double k, 
                    const double s, const double Mw_i,const double theta,
                               const double W,const double X,const double Y,const double GammaZ)const;
    
    //Total amplitude squared for l
    double MTOTl_sq(const StandardModel::lepton l,const double k, const double s, 
                    const double Mw_i,const double theta,
                               const double W,const double X,const double Y,const double GammaZ)const;
    
    
    double dsigmaBrem_l(const StandardModel::lepton l, const double k,
                                 const double s, const double Mw_i, const double theta,
                      const double W, const double X, const double Y,const double GammaZ)const;
    
    double dsigmaBrem_q(const QCD::quark q, const double k,
                                 const double s, const double Mw_i, const double theta,
                      const double W, const double X, const double Y,const double GammaZ)const;
    
    
    
    double dsigma_l(const StandardModel::lepton l, const double s, const double Mw_i,
                    const double theta, const double W,const double X,const double Y,const double GammaZ) const;
                              
    
    double dsigma_q(const QCD::quark q, const double s, const double Mw_i,const double theta,
                    const double W,const double X,const double Y,const double GammaZ) const;
    
    
    
    
protected:
    
private:    
    const EWSMcache& cache; 
    const EWSMOneLoopEW& EWOL;
    const PVfunctions PV;
    const StandardModel& SM;
    double phmass;//////////////////////////////////////////////////fictious photon mass

};

#endif	/* EWSMONELOOPLEP2_H */

