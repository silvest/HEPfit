/* 
 * File:   EWPOSM.h
 * Author: aleksandr
 *
 * Created on February 23, 2011, 5:48 PM
 */

#ifndef EWPOSM_H
#define	EWPOSM_H

#include <cstring>
#include "StandardModel.h"

#define LEPTON_TMP 1.0



class EWPOSM {
public:

    /**
     * @brief EWPOSM constructor
     * @param[in] sm reference to a StandardModel object
     */
    EWPOSM(StandardModel& sm);
   
    /**
     * @brief EWPOSM destructor
     */
    virtual ~EWPOSM();


    ///////////////////////////////////////////////////////////////////////////
    
    void Check_string(const std::string ferm);
    
    double Qf(const std::string ferm);
    
        
    ///////////////////////////////////////////////////////////////////////////
    // Analytic formulas for the W boson mass and the effective weak mixing angle
    
    double mW() const;

    double sw2() const;

    double cw2() const;

    /**
     * @param ferm defines the type of the fermions: leptons -"l", charm -"c", or b quark "b"
     * @return sin theta effective for the given fermion ferm @f$\sin^2\theta_{\mathrm{eff}}^\ell$@f
     */
    double sin2thwall(const std::string ferm) const;

    double vf(const std::string ferm);
    
    
    ///////////////////////////////////////////////////////////////////////////
    // Loop functions in Bardin et al., Z.Phys.C 44,493-502 (1989)    
  
    double lambda(double Q2, double M12, double M22);
    double lambdas(double s, double Mv2);
    double fcurve(double Q2, double M12, double M22);
    gslpp::complex fcurve_c(double Q2, double M12, double M22);
    gslpp::complex fcurve2(double s, double Mv2);
    gslpp::complex fcurve3(double s, double Mv2);
    gslpp::complex LL(double Q2, double M12, double M22);
    gslpp::complex I0(double Q2, double M12, double M22);
    gslpp::complex I1(double Q2, double M12, double M22);
    gslpp::complex I3(double Q2, double M12, double M22);

    
    ///////////////////////////////////////////////////////////////////////////
    // Loop functions copied from the ZFitter codes

    gslpp::complex DI0(double Q2, double AMZ2, double M12, double M22);
    gslpp::complex DI3(double Q2, double AMZ2, double M12, double M22);
    gslpp::complex DL(double Q2, double AMZ2, double M12, double M22);

    
    ///////////////////////////////////////////////////////////////////////////
    // Self energies of the W and Z bosons 
    
    double sigmaBF0();
    gslpp::complex sigmaBFMw();
    gslpp::complex sigmaBFMz();
    gslpp::complex sigmaBFMz1();
    gslpp::complex sigmaFww(double s);
    gslpp::complex sigmaFzz();
    gslpp::complex sigmaFzz1();
    gslpp::complex sigmaFzz1a();
    
    gslpp::complex delta_rho_W_F();
    gslpp::complex delta_rho_Z_F();
    
    gslpp::complex Delta_rho();
    

    ///////////////////////////////////////////////////////////////////////////
    // Z-gamma mixing 
    
    gslpp::complex piZg();
    
    gslpp::complex XAMM1();
    
    
    ///////////////////////////////////////////////////////////////////////////
    // Vertex functions
    
    gslpp::complex v1W(double s);
    gslpp::complex v1Z(double s);
    gslpp::complex v2W(double s);
    gslpp::complex v2Wa();
    
    gslpp::complex uf(const std::string ferm);

    
    ///////////////////////////////////////////////////////////////////////////
    // delta rho
    
    gslpp::complex delta_rho();

    double delta_rho_QCD();
    
    
    ///////////////////////////////////////////////////////////////////////////
    // delta r

    double delta_r() const;    
    
    
    ///////////////////////////////////////////////////////////////////////////
    // Reminder "delta rho_rem"

    gslpp::complex delta_rho_rem(const std::string ferm);

    
    ///////////////////////////////////////////////////////////////////////////
    // Reminder "delta kappa_rem"
        
    gslpp::complex delta_k_rem(const std::string ferm);
        

    ///////////////////////////////////////////////////////////////////////////
    // Reminder "Delta r_rem" 
    
    gslpp::complex Delta_r_rem();
    
    double F1(double x);
    double F1inf(double x);
    double Delta_r_rem_QCD();
    
    gslpp::complex Delta_r_rem_tot();    
    

    ///////////////////////////////////////////////////////////////////////////
    // Electroweak form factors

    /**
     * @param ferm can be equal to "u", "d", "l"
     * @return EW form factor rho
     */
    gslpp::complex rhoZf(const std::string ferm);
    
    /**
     * @param ferm can be equal to "u", "d", "l"
     * @return EW form factor k
     */
    gslpp::complex kZf(const std::string ferm);

    double rekZf(const std::string ferm);   
    
    gslpp::complex gzf(const std::string ferm);

    gslpp::complex gzf1(const std::string ferm);

    
    ///////////////////////////////////////////////////////////////////////////
            
    gslpp::complex spence(double);
    gslpp::complex fspence(double);
    
    double ROQCD();
    double AKQCD();
    double CLQQCD();

    
    ///////////////////////////////////////////////////////////////////////////    

protected:
    double me, mmu, mtau;
    double mu, md, ms, mc, mt, mb;
    double alsMz, GF, ale, dAletotmz, mZ, mHl;


};


#endif	/* EWPOSM_H */

