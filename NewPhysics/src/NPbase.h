/*
 * Copyright (C) 2013 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef NPBASE_H
#define	NPBASE_H

#include "StandardModel.h"

/**
 * @class NPbase
 * @brief The auxiliary base model class for other model classes.
 * @ingroup NewPhysics
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 *  
 * @details This is an auxiliary Model class containing the basic structure 
 * to compute new physics (NP) contributions to the electroweak precision
 * observables. The NP contributions are described by the following quantities:
 *
 * @li @f$\Delta G@f$&nbsp;&nbsp; (with DeltaGF()),
 * @li @f$S@f$, @f$T@f$ and @f$U@f$&nbsp;&nbsp;
 * (with obliqueS(), obliqueT() and obliqueU()),
 *
 * Using these quantities, the mass and width of the @f$W@f$ boson and the NP
 * contributions to the effective vector and axial-vector couplings of the
 * @f$Z@f$-boson to leptons and quarks are computed:
 *
 * @li @f$M_W@f$, @f$c_W^2@f$ and @f$s_W^2@f$&nbsp;&nbsp; (with Mw(), cW2() and sW2()),
 * @li @f$\Gamma_W@f$&nbsp;&nbsp; (with GammaW()),
 * @li @f$\delta g_V^f@f$&nbsp;&nbsp;(with deltaGVl() and deltaGVq()),
 * @li @f$\delta g_A^f@f$&nbsp;&nbsp;(with deltaGAl() and deltaGAq()). 
 *
 * In the model classes that are inherited from the current class (see the
 * inheritance diagram above), some of these methods are reimplemented to
 * account for the details of more specific scenarios.
 *
 * 
 * @anchor NPbaseInitialization
 * <h3>Initialization</h3>
 *
 * This class is intended to be used with an inherited model class.
 *
 *
 * @anchor NPbaseParameters
 * <h3>%Model parameters</h3>
 *
 * There is no model parameter in the current class.
 *
 * 
 * @anchor NPbaseFlags
 * <h3>%Model flags</h3>
 *
 * There is no model flag in the current class.
 *
 *
 * @anchor NPbaseFunctions
 * <h3>Important member functions</h3>
 *
 * The functions are explained above. 
 *
 */
class NPbase : public StandardModel {
public:

    /**
     * @brief The default constructor.
     */
    NPbase();

    /**
     * @brief The update method for %NPbase.
     * @return a boolean that is true if the execution is successful
     */
    virtual bool Update(const std::map<std::string, double>& DPars);
    
    /**
     * @brief The postupdate method for %NPbase.
     * @return a boolean that is true if the execution is successful
     */
    virtual bool PostUpdate();
    
    ////////////////////////////////////////////////////////////////////////

    /**
     * @brief A method to return a StandardModel object from %NPbase.
     * @return an StandardModel object
     */
    virtual StandardModel getTrueSM() const
    {
        return trueSM;
    }

    ////////////////////////////////////////////////////////////////////////

    /**
     * @brief New physics contribution to the Fermi constant.
     * @details The new physics contribution @f$\Delta G@f$ is defined as
     * @f[
     * G_\mu = G_{\mu,\mathrm{SM}}(1+\Delta G)\,,
     * @f]
     * where @f$G_\mu@f$ is the experimental value measured through muon decays, 
     * and @f$G_{\mu,\mathrm{SM}}@f$ is the Fermi constant in the SM.
     * @return @f$\Delta G@f$
     */
    virtual double DeltaGF() const
    {
        return 0.;
    }

    ////////////////////////////////////////////////////////////////////////

    /**
     * @brief The oblique parameter \f$S\f$.
     * @return the value of @f$S@f$
     */
    virtual double obliqueS() const
    {
        return 0.;
    }

    /**
     * @brief The oblique parameter \f$T\f$.
     * @return the value of @f$T@f$
     */
    virtual double obliqueT() const
    {
        return 0.;
    }

    /**
     * @brief The oblique parameter \f$U\f$.
     * @return the value of @f$U@f$
     */
    virtual double obliqueU() const
    {
        return 0.;
    }
    
    /**
     * @brief The oblique parameter \f$W\f$.
     * @return the value of @f$W@f$
     */
    virtual double obliqueW() const
    {
        return 0.;
    }
    
    /**
     * @brief The oblique parameter \f$Y\f$.
     * @return the value of @f$Y@f$
     */
    virtual double obliqueY() const
    {
        return 0.;
    }

    /**
     * @brief The electromagnetic coupling at the @f$Z@f$-mass scale, 
     * @f$\alpha(M_Z^2)=\alpha/(1-\Delta\alpha(M_Z^2))@f$.
     * @return @f$\alpha(M_Z^2)@f$
     */
    virtual double alphaMz() const;

    /**
     * @brief The mass of the @f$W@f$ boson, @f$M_W@f$.
     * @details
     * The @f$W@f$-boson mass receives the new physics
     * contribution via the oblique parameters @f$S@f$, @f$T@f$ and @f$U@f$ and
     * the shift in the Fermi constant, @f$\Delta G@f$:
     * @f[
     * M_W = M_{W,\mathrm{SM}}
     * \left[
     * 1 - \frac{\alpha(M_Z^2)}{4(c_W^2-s_W^2)}
     * \left( S - 2c_W^2\,T - \frac{c_W^2-s_W^2}{2s_W^2}\,U \right)
     * - \frac{s_W^2}{2(c_W^2-s_W^2)}\,\Delta G
     * \right].
     * @f]
     *
     * See @cite Ciuchini:2013pca and references therein.
     * @return @f$M_W@f$ in GeV
     */
    virtual double Mw() const;
    
    /**
     * @brief The new physics contribution to the decay width of the @f$W@f$ boson into a given fermion pair, @f$\delta \Gamma_Z^{f}@f$.
     * @param[in] fi a lepton or quark
     * @param[in] fj a lepton or quark
     * @return @f$\delta \Gamma_W^{ff}@f$ in GeV
     */
    virtual double deltaGamma_Wff(const Particle fi, const Particle fj) const
    {
        return 0.0;
    };
    
    /**
     * @brief A partial decay width of the @f$W@f$ boson decay into a SM fermion pair.
     * @details
     * The partial @f$W@f$-boson widths receives the new physics
     * contribution via the oblique parameters @f$S@f$, @f$T@f$ and @f$U@f$ and
     * the shift in the Fermi constant, @f$\Delta G@f$:
     * @f[
     * \Gamma_W^{ij} = \Gamma_{W,\mathrm{SM}}
     * \left[ 1
     * - \frac{3\alpha(M_Z^2)}{4(c_W^2-s_W^2)}
     *  \left( S - 2c_W^2\,T - \frac{c_W^2-s_W^2}{2s_W^2}\,U \right)
     * - \frac{1+c_W^2}{2(c_W^2-s_W^2)}\, \Delta G
     * \right].
     * @f]
     * @param[in] fi a lepton or quark
     * @param[in] fj a lepton or quark
     * @return @f$\Gamma^W_{ij}@f$
     *
     * @attention Fermion masses are neglected.
     */
    virtual double GammaW(const Particle fi, const Particle fj) const;
    
    /**
     * @brief The new physics contribution to the total decay width of the @f$W@f$ boson, @f$\delta \Gamma_W@f$.
     * @return @f$\delta \Gamma_W@f$ in GeV
     */
    virtual double deltaGamma_W() const
    {
        return 0.0;
    };

    /**
     * @brief The total width of the @f$W@f$ boson, @f$\Gamma_W@f$.
     * @details
     * The @f$W@f$-boson width receives the new physics
     * contribution via the oblique parameters @f$S@f$, @f$T@f$ and @f$U@f$ and
     * the shift in the Fermi constant, @f$\Delta G@f$:
     * @f[
     * \Gamma_W = \Gamma_{W,\mathrm{SM}}
     * \left[ 1
     * - \frac{3\alpha(M_Z^2)}{4(c_W^2-s_W^2)}
     *  \left( S - 2c_W^2\,T - \frac{c_W^2-s_W^2}{2s_W^2}\,U \right)
     * - \frac{1+c_W^2}{2(c_W^2-s_W^2)}\, \Delta G
     * \right].
     * @f]
     *
     * See @cite Ciuchini:2013pca and references therein.
     * @return @f$\Gamma_W@f$ in GeV
     */
    virtual double GammaW() const;
    
    /**
     * @brief The branching ratio of the @f$W@f$ boson decaying into a SM fermion pair, @f$Br(W\to f_i f_j)@f$.
     * @return @f$Br(W\to f_i f_j)@f$ in GeV
     */
    virtual double BrW(const Particle fi, const Particle fj) const;
    
    /**
     * @brief The lepton universality ratio @f$R_{W,l_i/l_j)=\Gamma(W\to l_i \nu_i)/\Gamma(W\to l_j \nu_j)@f$.
     * @return @f$R_{W,l_i/l_j)@f$ in GeV
     */
    virtual double RWlilj(const Particle li, const Particle lj) const;
    
    /**
     * @brief The ratio @f$R_{W,c)=\Gamma(W\to c + X)/\Gamma(W\to had)@f$.
     * @return @f$R_{W,c)@f$ in GeV
     */
    virtual double RWc() const;

    /**
     * @brief 
     */    
    virtual double deltaGV_f_2(const Particle f) const {
        return 0.0;
    }  //AG:added
    
    /**
     * @brief New physics contribution to the neutral-current vector coupling @f$g_V^f@f$.
     * @details
     * The neutral-current vector coupling @f$g_V^f@f$ receives the new physics
     * contribution via the oblique parameters @f$S@f$ and @f$T@f$ and the shift
     * in the Fermi constant, @f$\Delta G@f$:
     * @f[
     * \delta g_V^f =
     * \frac{g_{V,\mathrm{SM}}^f}{2}
     * \left[ \alpha(M_Z^2)\, T - \Delta G \right]
     * +
     * \frac{\big( g_{V,\mathrm{SM}}^f - g_{A,\mathrm{SM}}^f \big)
     * \left[
     * \alpha(M_Z^2)\left( S - 4\,c_W^2s_W^2\, T \right)
     * + 4\,c_W^2s_W^2\, \Delta G
     * \right]}{4s_W^2\,(c_W^2-s_W^2)}\,.
     * @f]
     *
     * See @cite Ciuchini:2013pca and references therein.
     * @param[in] f a lepton or quark
     * @return @f$\delta g_V^f@f$
     */
    virtual double deltaGV_f(const Particle f) const;

    /**
     * @brief The total (SM+NP) contribution to the neutral-current vector coupling @f$g_V^f@f$.
     * @param[in] f a lepton or quark
     * @return @f$g_V^f@f$, including SM plus NP contributions
     */
    virtual gslpp::complex gV_f(const Particle f) const;

    /**
     * @brief 
     */    
    virtual double deltaGA_f_2(const Particle f) const {
        return 0.0;
    }  //AG:added    

    /**
     * @brief New physics contribution to the neutral-current axial-vector coupling @f$g_A^f@f$.
     * @details
     * The neutral-current axial-vector coupling @f$g_A^f@f$ receives the new
     * physics contribution via the oblique parameter @f$T@f$ and the shift in
     * the Fermi constant, @f$\Delta G@f$:
     * @f[
     * \delta g_A^f
     * = \frac{g_{A,\mathrm{SM}}^f}{2} \left[ \alpha(M_Z^2)\, T - \Delta G \right].
     * @f]
     *
     * See @cite Ciuchini:2013pca and references therein.
     * @param[in] f a lepton or quark
     * @return @f$\delta g_A^f@f$
     */
    virtual double deltaGA_f(const Particle f) const;

    /**
     * @brief The total (SM+NP) contribution to the neutral-current axial-vector coupling @f$g_A^f@f$.
     * @param[in] f a lepton or quark
     * @return @f$g_A^f@f$, including SM plus NP contributions
     */
    virtual gslpp::complex gA_f(const Particle f) const;

    /**
     * @brief The effective neutral-current coupling @f$\rho_Z^f@f$ including SM plus NP contributions.
     * @param[in] f a lepton or quark
     * @return @f$\rho_Z^f@f$, including SM plus NP contributions
     */
    virtual gslpp::complex rhoZ_f(const Particle f) const;

    /**
     * @brief The effective neutral-current coupling @f$\kappa_Z^f@f$ including SM plus NP contributions.
     * @param[in] f a lepton or quark
     * @return @f$\kappa_Z^f@f$, including SM plus NP contributions
     */
    virtual gslpp::complex kappaZ_f(const Particle f) const;
    
    /**
     * @brief The \f$\mathcal{O}(\Lambda^{-4})\f$ new physics contribution to the decay width of the \f$Z\f$ boson into a given fermion pair, @f$\Delta \Gamma_{Z,f}^{(2)}@f$.
     * @details 
     * \f$\newcommand{\nc}{\newcommand}\f$
     * \f$\nc{\gVfSM}{g_{V,f}^{SM}}\f$
     * \f$\nc{\gAfSM}{g_{A,f}^{SM}}\f$
     * \f$\nc{\gVfL}{\Delta g_{V,f}^{(1)}}\f$
     * \f$\nc{\gAfL}{\Delta g_{A,f}^{(1)}}\f$
     * \f$\nc{\gVfQ}{\Delta g_{V,f}^{(2)}}\f$
     * \f$\nc{\gAfQ}{\Delta g_{A,f}^{(2)}}\f$
     * \f[ 
     * \Delta \Gamma_{Z,f}^{(2)} = N_f \frac{\alpha(M_Z)^{SM} M_Z}{6 s_W^2 c_w^2} \left[ (\gVfL)^2 + (\gAfL)^2 + 2 \left( \gAfSM ~\gAfQ + \gVfSM ~\gVfQ \right) \right] 
     * \f]
     * @param[in] f a lepton or quark
     * @return @f$\Delta \Gamma_{Z,f}^{(2)}@f$ in GeV
     */
    virtual double deltaGamma_Zf_2(const Particle f) const;       //AG:added
    
    /**
     * @brief The new physics contribution to the decay width of the @f$Z@f$ boson into a given fermion pair, @f$\delta \Gamma_Z^{f}@f$.
     * @param[in] f a lepton or quark
     * @return @f$\delta \Gamma_Z^{f}@f$ in GeV
     */
    virtual double deltaGamma_Zf(const Particle f) const;
    
    /**
     * @brief The decay width of the @f$Z@f$ boson into a given fermion pair, @f$\Gamma_Z^{f}@f$.
     * @details
     * \f[
     * \Gamma_Z^{f} = \Gamma_{Z,f}^{SM} + \Delta \Gamma_{Z,f}^{(1)} + \Delta \Gamma_{Z,f}^{(2)}
     * \f]
     * @param[in] f a lepton or quark
     * @return @f$\Gamma_Z^{f}@f$ in GeV, including SM plus @f$\mathcal{O}(\Lambda^{-2})@f$ and @f$\mathcal{O}(\Lambda^{-4})@f$ NP contributions
     *
     * @attention This function is applicable only to the NP model classes that
     * are inherited from NPbase.
     */
    virtual double Gamma_Zf(const Particle f) const;               //AG:modified

    /**
     * @brief The @f$\mathcal{O}(\Lambda^{-4})@f$ new physics contribution to the total decay width of the @f$Z@f$ boson, @f$\Delta \Gamma_Z^{(2)}@f$.
     * @details 
     * \f[ 
     * \Delta \Gamma_{Z}^{(2)} = \Delta \Gamma_{Z,fer}^{(2)} + \Delta \Gamma_{Z,had}^{(2)}
     * \f]
     * @return @f$\Delta \Gamma_Z^{(2)}@f$ in GeV
     */
    virtual double deltaGamma_Z_2() const;                         //AG:added
    
    /**
     * @brief The new physics contribution to the total decay width of the @f$Z@f$ boson, @f$\delta \Gamma_Z@f$.
     * @return @f$\delta \Gamma_Z@f$ in GeV
     */
    virtual double deltaGamma_Z() const;

    /**
     * @brief The total decay width of the @f$Z@f$ boson, @f$\Gamma_Z@f$.
     * @details
     * \f[
     * \Gamma_Z = \Gamma_Z^{SM} + \Delta \Gamma_Z^{(1)} + \Delta \Gamma_Z^{(2)}
     * \f]
     * @return @f$\Gamma_Z@f$ in GeV, including SM plus @f$\mathcal{O}(\Lambda^{-2})@f$ and @f$\mathcal{O}(\Lambda^{-4})@f$ NP contributions
     *
     * @attention This function is applicable only to the NP model classes that
     * are inherited from NPbase.
     */
    virtual double Gamma_Z() const;                                //AG:modified
    
    /**
     * @brief The @f$\mathcal{O}(\Lambda^{-4})@f$ new physics contribution to the hadronic decay width of the @f$Z@f$ boson, @f$\Delta \Gamma_{Z,had}^{(2)}@f$.
     * @details 
     * \f[ 
     * \Delta \Gamma_{Z,had}^{(2)} = \Delta \Gamma_{Z,u}^{(2)} + \Delta \Gamma_{Z,d}^{(2)}
     * + \Delta \Gamma_{Z,c}^{(2)} + \Delta \Gamma_{Z,s}^{(2)}
     * + \Delta \Gamma_{Z,b}^{(2)}
     * \f]
     * @return @f$\Delta \Gamma_{Z,had}^{(2)}@f$ in GeV
     */
    virtual double deltaGamma_Zhad_2() const;                      //AG:added
    
    /**
     * @brief The new physics contribution to the hadronic decay width of the @f$Z@f$ boson, @f$\delta \Gamma_{Z,had}@f$.
     * @return @f$\delta \Gamma_{Z,had}@f$ in GeV
     */
    virtual double deltaGamma_Zhad() const;
    
    /**
     * @brief The hadronic decay width of the @f$Z@f$ boson, @f$\Gamma_{Z,had}@f$.
     * @details
     * \Gamma_{Z,had} = \Gamma_Z,had}^{SM} + \Delta \Gamma_{Z,had}^{(1)} + \Delta \Gamma_{Z,had}^{(2)}
     * @return @f$\Gamma_{Z,had}@f$ in GeV, including SM plus @f$\mathcal{O}(\Lambda^{-2})@f$ and @f$\mathcal{O}(\Lambda^{-4})@f$ NP contributions
     *
     * @attention This function is applicable only to the NP model classes that
     * are inherited from NPbase.
     */
    virtual double Gamma_had() const;                              //AG:modified
    
    /**
     * @brief The new physics contribution to the ratio of the @f$Z\to u\bar{u} + Z\to c\bar{c}@f$ width to the @f$Z@f$-boson hadronic width:
     * @details
     * \f[ 
     * \Delta R_{uc} = \frac{1}{2} \left(\Delta R_u^{(1)} + \Delta R_c^{(1)} \right) 
     * \f]
     * @return @f$\Delta R_{uc}^{(1)}@f$
     */
    virtual double deltaRuc() const;                               //AG:added
    
    /**
     * @brief The @f$\mathcal{O}(1/\Lambda^4)@f$ new physics contribution to the ratio of the @f$Z\to u\bar{u} + Z\to c\bar{c}@f$ width to the @f$Z@f$-boson hadronic width:
     * @details
     * \f[ 
     * \Delta R_{uc}^{(2)} = \frac{1}{2} \left(\Delta R_u^{(2)} + \Delta R_c^{(2)} \right) 
     * \f]
     * @return @f$\Delta R_{uc}^{(2)}@f$
     */
    virtual double deltaRuc_2() const;                             //AG:added
    
    /**
     * @brief The ratio of the @f$Z\to u\bar{u} + Z\to c\bar{c}@f$ width to the @f$Z@f$-boson hadronic width.
     * @return @f$\Delta R_{uc}^{0 (2)}@f$, including SM plus NP contributions
     */
    virtual double Ruc() const;                                    //AG:added
    
    /**
     * @brief The lepton universality ratio @f$R_{Z,l_i/l_j)=\Gamma(Z\to l_i^+ l_i^-)/\Gamma(Z\to l_j^+ l_j^-)@f$.
     * @return @f$R_{Z,l_i/l_j)@f$ in GeV
     */
    virtual double RZlilj(const Particle li, const Particle lj) const;

    /**
     * @brief The Branching ratio of the @f$Z@f$ boson into a given fermion pair, @f$BR_Z^{f}@f$.
     * @param[in] f a lepton or quark
     * @return @f$BR_Z^{f}@f$ including SM plus NP contributions
     *
     * @attention This function is applicable only to the NP model classes that
     * are inherited from NPbase.
     */
    virtual double BR_Zf(const Particle f) const;

    /**
     * @brief The @f$\mathcal{O}(\Lambda^{-4})@f$ new physics contribution to the cross section for the process @f$e^+ e^-\to Z\to \mathrm{hadrons}@f$
     * at the @f$Z@f$ pole, @f$\Delta \sigma_h^{0,(2)}@f$.
     * @details
     * \f$\nc{\MZt}{\widetilde{M}_Z}\f$
     * \f[
     * \Delta\sigma_{h}^{0,(2)} = \frac{12\pi}{\MZt^2} ~\frac{\Gamma_e^{SM} \Gamma_{had}^{SM}}{(\Gamma_Z^{SM})^2}
    ~ \left( \frac{\Delta\Gamma_e^{(2)}}{\Gamma_e^{SM}}
           + \frac{\Delta\Gamma_{had}^{(2)}}{\Gamma_{had}^{SM}}
           - 2 \frac{\Delta\Gamma_Z^{(2)}}{\Gamma_Z^{SM}}
           + \frac{\Delta\Gamma_e^{(1)} \Delta\Gamma_{had}^{(1)}}{\Gamma_e^{SM} \Gamma_{had}^{SM}} 
           - 2 \frac{\Delta\Gamma_e^{(1)} \Delta\Gamma_Z^{(1)}}{\Gamma_e^{SM} \Gamma_Z^{SM}} 
           - 2 \frac{\Delta\Gamma_{had}^{(1)} \Delta\Gamma_Z^{(1)}}{\Gamma_{had}^{SM} \Gamma_Z^{SM}} 
           + 3 \frac{(\Delta\Gamma_Z^{(1)})^2}{(\Gamma_Z^{SM})^2} 
     \right)
     * \f]
     * @return @f$\Delta \sigma_h^{0,(2)}@f$ in GeV@f$^{-2}@f$
     */
    virtual double deltaSigmaHadron_2() const;                     //AG:added
    
    /**
     * @brief The new physics contribution to the cross section for the process @f$e^+ e^-\to Z\to \mathrm{hadrons}@f$
     * at the @f$Z@f$ pole, @f$\delta \sigma_h^0@f$.
     * @return @f$\delta \sigma_h^0@f$ in GeV@f$^{-2}@f$
     */
    virtual double deltaSigmaHadron() const;

    /**
     * @brief The cross section for the process @f$e^+ e^-\to Z\to \mathrm{hadrons}@f$
     * at the @f$Z@f$ pole, @f$\sigma_h^0@f$.
     * @details
     * \f[
     * \sigma_h = \sigma_h^{SM} + \Delta \sigma_h^{(1)} + \Delta \sigma_h^{(2)}
     * \f]
     * @return @f$\sigma_h^0@f$ in GeV@f$^{-2}@f$, including SM plus @f$\mathcal{O}(\Lambda^{-2})@f$ and @f$\mathcal{O}(\Lambda^{-4})@f$ NP contributions
     *
     * @attention This function is applicable only to the NP model classes that
     * are inherited from NPbase.
     */
    virtual double sigma0_had() const;                             //AG:modified

    /**
     * @brief The @f$\mathcal{O}(\Lambda^{-4})@f$ new physics contribution to the effective electron weak angle @f$\Delta \sin^2\theta_{eff,e}^{(2)}@f$
     * at the @f$Z@f$ pole.
     * @details 
     * \f$\nc{\gVlSM}{g_{V,e}^{SM}}\f$
     * \f$\nc{\gAlSM}{g_{A,e}^{SM}}\f$
     * \f$\nc{\gVlL}{\Delta g_{V,e}^{(1)}}\f$
     * \f$\nc{\gAlL}{\Delta g_{A,e}^{(1)}}\f$
     * \f$\nc{\gVlQ}{\Delta g_{V,e}^{(2)}}\f$
     * \f$\nc{\gAlQ}{\Delta g_{A,e}^{(2)}}\f$
     * \f[ 
     * \Delta sin^2\theta_{eff,e}^{(2)} = \frac{1}{4}\frac{\gVlL ~\gAlL}{(\gAlSM)^2}  -\frac{1}{4}\frac{\gVlSM ~(\gAlL)^2}{(\gAlSM)^3} - \frac{1}{4} \frac{ \gAlSM ~\gVlQ - \gVlSM ~\gAlQ}{(\gAlSM)^2}
     * \f]
     * @return @f$\Delta \sin^2\theta_{eff,e}^{(2)}@f$
     */
    virtual double deltaSin2thetaEff_e_2() const;                  //AG:added
    
    /**
     * @brief The new physics contribution to the effective electron/leptonic weak angle @f$\delta \sin^2\theta_{\rm eff}^{\rm lept}@f$
     * at the @f$Z@f$ pole.
     * @return @f$\delta \sin^2\theta_{\rm eff}^{\rm lept}@f$
     */
    virtual double deltaSin2thetaEff_e() const;
    
    /**
     * @brief The @f$\mathcal{O}(\Lambda^{-4})@f$ new physics contribution to the effective muonic weak angle @f$\Delta \sin^2\theta_{eff, \mu}^{(2)}@f$
     * at the @f$Z@f$ pole.
     * @details 
     * \f$\nc{\gVlSM}{g_{V,\mu}^{SM}}\f$
     * \f$\nc{\gAlSM}{g_{A,\mu}^{SM}}\f$
     * \f$\nc{\gVlL}{\Delta g_{V,\mu}^{(1)}}\f$
     * \f$\nc{\gAlL}{\Delta g_{A,\mu}^{(1)}}\f$
     * \f$\nc{\gVlQ}{\Delta g_{V,\mu}^{(2)}}\f$
     * \f$\nc{\gAlQ}{\Delta g_{A,\mu}^{(2)}}\f$
     * \f[ 
     * \Delta sin^2\theta_{eff,\mu}^{(2)} = \frac{1}{4}\frac{\gVlL ~\gAlL}{(\gAlSM)^2}  -\frac{1}{4}\frac{\gVlSM ~(\gAlL)^2}{(\gAlSM)^3} - \frac{1}{4} \frac{ \gAlSM ~\gVlQ - \gVlSM ~\gAlQ}{(\gAlSM)^2}
     * \f]
     * @return @f$\Delta \sin^2\theta_{eff,\mu}^{(2)}@f$
     */
    virtual double deltaSin2thetaEff_mu_2() const;                 //AG:added
    
    /**
     * @brief The new physics contribution to the effective muonic weak angle @f$\delta \sin^2\theta_{\rm eff}^{\mu\mu}@f$
     * at the @f$Z@f$ pole.
     * @return @f$\delta \sin^2\theta_{\rm eff}^{\mu\mu}@f$
     */
    virtual double deltaSin2thetaEff_mu() const;

    /**
     * @brief @copybrief sin2thetaEff::computeThValue()
     * @details
     * \f[
     * \sin^2\theta_{eff} = \sin^2\theta_{eff}^{SM} + \Delta \sin^2\theta_{eff}^{(1)} + \Delta \sin^2\theta_{eff}^{(2)}
     * \f]
     * @param[in] f a lepton or quark
     * @return @f$\sin^2\theta_{\rm eff}^{\rm lept}@f$, including SM plus @f$\mathcal{O}(\Lambda^{-2})@f$ and @f$\mathcal{O}(\Lambda^{-4})@f$ NP contributions
     *
     * @attention This function is applicable only to the NP model classes that
     * are inherited from NPbase.
     */
    virtual double sin2thetaEff(const Particle f) const;           //AG:modified

    /**
     * @brief The @f$\mathcal{O}(\Lambda^{-4})@f$ new physics contribution to the left-right asymmetry in @f$e^+e^-\to Z\to f \bar{f}@f$ at the
     * @f$Z@f$-pole, @f$\Delta \mathcal{A}_f^{(2)}@f$.
     * @details
     * \f$\newcommand{\nc}{\newcommand}\f$
     * \f$\nc{\gVfSM}{g_{V,f}^{SM}}\f$
     * \f$\nc{\gAfSM}{g_{A,f}^{SM}}\f$
     * \f$\nc{\gVfL}{\Delta g_{V,f}^{(1)}}\f$
     * \f$\nc{\gAfL}{\Delta g_{A,f}^{(1)}}\f$
     * \f$\nc{\gVfQ}{\Delta g_{V,f}^{(2)}}\f$
     * \f$\nc{\gAfQ}{\Delta g_{A,f}^{(2)}}\f$
     * \f{eqnarray*}{
     * \Delta \mathcal{A}_f^{(2)} &=& 2~ \frac{ \gVfSM ~\gAfSM }{ \left[(g_{V,f}^{SM})^2+(g_{A,f}^{SM})^2\right]^3 } 
     * \left(  \left[(\gAfSM)^2 -3 (\gVfSM)^2\right] (\gAfL)^2 
     * + \left[(\gVfSM)^2 -3 (\gAfSM)^2\right] (\gVfL)^2 \right) \\
     * &-& ~ 2~ \frac{ (\gAfSM)^4 - 6(\gAfSM)^2(\gVfSM)^2 + (\gVfSM)^4 }{\left[(g_{V,f}^{SM})^2+(g_{A,f}^{SM})^2\right]^3} ~\gVfL ~\gAfL 
     * +~ 2~ \frac{ (g_{V,f}^{SM})^2-(g_{A,f}^{SM})^2 }{ \left[(g_{V,f}^{SM})^2+(g_{A,f}^{SM})^2\right]^2 } \left[ \gVfSM ~\gAfQ - \gAfSM ~\gVfQ \right]
     * \f}
     * @param[in] f a lepton or quark
     * @return @f$\Delta \mathcal{A}_f^{(2)}@f$
     */
    virtual double deltaA_f_2(const Particle f) const;             //AG:added
    
    /**
     * @brief The new physics contribution to the left-right asymmetry in @f$e^+e^-\to Z\to f \bar{f}@f$ at the
     * @f$Z@f$-pole, @f$\delta \mathcal{A}_f@f$. 
     * @param[in] f a lepton or quark
     * @return @f$\delta \mathcal{A}_f@f$
     */
    virtual double deltaA_f(const Particle f) const;

    /**
     * @brief The left-right asymmetry in @f$e^+e^-\to Z\to f \bar{f}@f$ at the
     * @f$Z@f$-pole, @f$\mathcal{A}_f@f$.
     * @details
     * \f[
     * \mathcal{A}_f = \mathcal{A}_f^{SM} + \Delta \mathcal{A}_f^{(1)} + \Delta \mathcal{A}_f^{(2)}
     * \f]
     * @param[in] f a lepton or quark
     * @return @f$\mathcal{A}_f@f$, including SM plus @f$\mathcal{O}(\Lambda^{-2})@f$ and @f$\mathcal{O}(\Lambda^{-4})@f$ NP contributions
     *
     * @attention This function is applicable only to the NP model classes that
     * are inherited from NPbase.
     */
    virtual double A_f(const Particle f) const;                    //AG:modified

    /**
     * @brief The @f$\mathcal{O}(\Lambda^{-4})@f$ new physics to the forward-backward asymmetry in @f$e^+e^-\to Z\to f \bar{f}@f$ at the
     * @f$Z@f$-pole, @f$\Delta A^f_{FB}@f$. 
     * @details
     * \f[
     * \Delta A_{FB,f}^{(2)} = \frac{3}{4} \left( \Delta A_e^{(1)} \Delta A_f^{(1)} + A_e^{SM} \Delta A_f^{(2)} + A_f^{SM} \Delta A_e^{(2)} \right)
     * \f]
     * @param[in] f a lepton or quark
     * @return @f$\Delta A^{f (2)}_{FB}@f$
     */
    virtual double deltaAFB_2(const Particle f) const;             //AG:added
    
    /**
     * @brief The new physics contribution to the forward-backward asymmetry in @f$e^+e^-\to Z\to f \bar{f}@f$ at the
     * @f$Z@f$-pole, @f$\delta A^f_{FB}@f$. 
     * @param[in] f a lepton or quark
     * @return @f$\delta A^f_{FB}@f$
     */
    virtual double deltaAFB(const Particle f) const;

    /**
     * @brief The forward-backward asymmetry in @f$e^+e^-\to Z\to f \bar{f}@f$ at the
     * @f$Z@f$-pole, @f$A^f_{FB}@f$.
     * @details
     * \f[
     * A_{FB,f} = A_{FB,f}^{SM} + \Delta A_{FB,f}^{(1)} + \Delta A_{FB,f}^{(2)}
     * \f]
     * @param[in] f a lepton or quark
     * @return @f$A^f_{FB}@f$, including SM plus @f$\mathcal{O}(\Lambda^{-2})@f$ and @f$\mathcal{O}(\Lambda^{-4})@f$ NP contributions
     */
    virtual double AFB(const Particle f) const;                    //AG:modified

    /**
     * @brief The @f$\mathcal{O}(\Lambda^{-4})@f$ new physics contribution to the ratio @f$R_\ell^0=\Gamma_{\mathrm{had}}/\Gamma_\ell@f$,
     * @f$R_q^0=\Gamma_q/\Gamma_{\mathrm{had}}@f$ and @f$R_\nu^0=\Gamma_\nu/\Gamma_{\mathrm{had}}@f$, 
     * for charged leptons, quarks and neutrinos:
     * @details
     * \f[
     * \Delta R_l^{(2)} = \frac{ \Gamma_{had}^{SM}~(\Delta \Gamma_l^{(1)})^2}{(\Gamma_l^{SM})^3}
                        - \frac{\Delta \Gamma_{had}^{(1)}~\Delta \Gamma_l^{(1)}}{(\Gamma_l^{SM})^2}
                        + \frac{ \Gamma_l^{SM} ~\Delta \Gamma_{had}^{(2)} - \Gamma_{had}^{SM} ~\Delta \Gamma_l^{(2)}}{(\Gamma_l^{SM})^2}
     * \f]
     * \f[
     * \Delta R_q^{(2)} = \frac{ \Gamma_{q}^{SM}~(\Delta \Gamma_{had}^{(1)})^2}{(\Gamma_{had}^{SM})^3}
                     - \frac{\Delta \Gamma_q^{(1)}~\Delta \Gamma_{had}^{(1)}}{(\Gamma_{had}^{SM})^2}
                     + \frac{ \Gamma_{had}^{SM} ~\Delta \Gamma_q^{(2)} - \Gamma_q^{SM} ~\Delta \Gamma_{had}^{(2)}}{(\Gamma_{had}^{SM})^2}
     * \f]
     * \f[
     * \Delta R_{\nu}^{(2)} = \frac{ \Gamma_{\nu}^{SM}~(\Delta \Gamma_{had}^{(1)})^2}{(\Gamma_{had}^{SM})^3}
                     - \frac{\Delta \Gamma_{\nu}^{(1)}~\Delta \Gamma_{had}^{(1)}}{(\Gamma_{had}^{SM})^2}
                     + \frac{ \Gamma_{had}^{SM} ~\Delta \Gamma_{\nu}^{(2)} - \Gamma_{\nu}^{SM} ~\Delta \Gamma_{had}^{(2)}}{(\Gamma_{had}^{SM})^2}
     * \f]
     * @param f a lepton or quark
     * @return @f$\Delta R_f^{0 (2)}@f$
     */
    virtual double deltaR0_f_2(const Particle f) const;            //AG:added
    
    /**
     * @brief The new physics contribution to the ratio @f$R_\ell^0=\Gamma_{\mathrm{had}}/\Gamma_\ell@f$,
     * @f$R_q^0=\Gamma_q/\Gamma_{\mathrm{had}}@f$ and @f$R_\nu^0=\Gamma_\nu/\Gamma_{\mathrm{had}}@f$, 
     * for charged leptons, quarks and neutrinos, respectively.
     * @param f a lepton or quark
     * @return @f$\delta R_f^0@f$
     */
    virtual double deltaR0_f(const Particle f) const;

    /**
     * @brief The ratio @f$R_\ell^0=\Gamma_{\mathrm{had}}/\Gamma_\ell@f$,
     * @f$R_q^0=\Gamma_q/\Gamma_{\mathrm{had}}@f$ and @f$R_\nu^0=\Gamma_\nu/\Gamma_{\mathrm{had}}@f$, 
     * for charged leptons, quarks and neutrinos, respectively. 
     * @details
     * \f[
     * R^0_f = R_f^{SM} + \Delta R_f^{(1)} + \Delta R_f^{(2)}
     * \f]
     * @param[in] f a lepton or quark
     * @return @f$R_f^0@f$, including SM plus @f$\mathcal{O}(\Lambda^{-2})@f$ and @f$\mathcal{O}(\Lambda^{-4})@f$ NP contributions
     */
    virtual double R0_f(const Particle f) const;                   //AG:modified
    
    /**
     * @brief The new physics contribution to the ratio of invisible and leptonic (electron) decay widths of the @f$Z@f$ boson, @f$\delta R_{inv}@f$.
     * @return @f$\delta R_{inv}@f$
     */
    virtual double deltaR_inv() const;

    /**
     * @brief The ratio of the invisible and leptonic (electron) decay widths of the @f$Z@f$ boson, @f$R_{inv}@f$.
     * @return @f$R_{inv}@f$, including SM plus NP contributions
     *
     * @attention This function is applicable only to the NP model classes that
     * are inherited from NPbase.
     */
    virtual double R_inv() const;
    
    /**
     * @brief The new physics contribution to the number of neutrinos dervied from the @f$Z@f$ pole measurements.
     * @return @f$\delta N_{\nu}@f$
     */
    virtual double deltaN_nu() const;

    /**
     * @brief The number of neutrinos dervied from the @f$Z@f$ pole measurements, @f$N_{\nu}@f$.
     * @return @f$N_{\nu}@f$, including SM plus NP contributions
     *
     * @attention This function is applicable only to the NP model classes that
     * are inherited from NPbase.
     */
    virtual double N_nu() const;

    ////////////////////////////////////////////////////////////////////////
    
    /**
     * @brief New physics contribution to the charged current coupling @f$W_\mu \bar{f_L}\gamma^mu f_L@f$.
     * @param[in] pbar a lepton or quark
     * @param[in] p a lepton or quark
     * @return @f$\delta g_{Wff}^{L}@f$
     */
    // no generation mixing
    virtual gslpp::complex deltaGL_Wff(const Particle pbar, const Particle p) const
    {
        return 0.0;
    };
    /**
     * @brief New physics contribution to the charged current coupling @f$W_\mu \bar{f_R}\gamma^mu f_R@f$.
     * @param[in] pbar a lepton or quark
     * @param[in] p a lepton or quark
     * @return @f$\delta g_{Wff}^{R}@f$
     */
    // no generation mixing
    virtual gslpp::complex deltaGR_Wff(const Particle pbar, const Particle p) const
    {
        return 0.0;
    };

    /**
     * @brief The new physics contribution to the coupling of the effective interaction @f$H G_{\mu\nu}^AG^{A \mu\nu}@f$.
     * @return @f$\delta g_{HGG}@f$
     */
    virtual double deltaG_hgg() const
    {
        return 0.0;
    };
    /**
     * @brief The full new physics contribution to the coupling of the effective interaction @f$H G_{\mu\nu}^AG^{A \mu\nu}@f$,
     * including new local terms and modifications on the SM-loops. Normalized to the SM value.
     * @return @f$\delta g_{HGG}/g_{HGG}^SM}@f$
     */
    virtual double deltaG_hggRatio() const
    {
        return 0.0;
    };
    /**
     * @brief The new physics contribution to the coupling of the effective interaction @f$H W_{\mu\nu}^\dagger W^{\mu\nu}@f$.
     * @return @f$\delta g_{HWW}^{(1)}@f$
     */
    virtual double deltaG1_hWW() const
    {
        return 0.0;
    };
    /**
     * @brief The new physics contribution to the coupling of the effective interaction @f$H W_{\nu}^\dagger \partial^\mu W^{\mu\nu}@f$.
     * @return @f$\delta g_{HWW}^{(2)}@f$
     */
    virtual double deltaG2_hWW() const
    {
        return 0.0;
    };
    /**
     * @brief The new physics contribution to the coupling of the effective interaction @f$H W_{\mu}^\dagger W^{\mu}@f$.
     * @return @f$\delta g_{HWW}^{(3)}@f$
     */
    virtual double deltaG3_hWW() const
    {
        return 0.0;
    };
    /**
     * @brief The new physics contribution to the coupling of the effective interaction @f$H Z_{\mu\nu} Z^{\mu\nu}@f$.
     * @return @f$\delta g_{HZZ}^{(1)}@f$
     */
    virtual double deltaG1_hZZ() const
    {
        return 0.0;
    };
    /**
     * @brief The new physics contribution to the coupling of the effective interaction @f$H Z_{\nu} \partial^\mu Z^{\mu\nu}@f$.
     * @return @f$\delta g_{HZZ}^{(2)}@f$
     */
    virtual double deltaG2_hZZ() const
    {
        return 0.0;
    };
    /**
     * @brief The new physics contribution to the coupling of the effective interaction @f$H Z_{\mu} Z^{\mu}@f$.
     * @return @f$\delta g_{HZZ}^{(3)}@f$
     */
    virtual double deltaG3_hZZ() const
    {
        return 0.0;
    };
    /**
     * @brief The new physics contribution to the coupling of the effective interaction @f$H Z_{\mu\nu} F^{\mu\nu}@f$.
     * @return @f$\delta g_{HZA}^{(1)}@f$
     */
    virtual double deltaG1_hZA() const
    {
        return 0.0;
    };
    /**
     * @brief The full new physics contribution to the coupling of the effective interaction @f$H Z_{\mu\nu} F^{A \mu\nu}@f$,
     * including new local terms and modifications on the SM-loops. Normalized to the SM value.
     * @return @f$\delta g_{HZA}^{(1)}/g_{HZA}^{(1),SM}@f$
     */
    virtual double deltaG1_hZARatio() const
    {
        return 0.0;
    };
    /**
     * @brief The new physics contribution to the coupling of the effective interaction @f$H Z_{\nu} \partial^\mu F^{\mu\nu}@f$.
     * @return @f$\delta g_{HZA}^{(2)}@f$
     */
    virtual double deltaG2_hZA() const
    {
        return 0.0;
    };
    /**
     * @brief The new physics contribution to the coupling of the effective interaction @f$H F_{\mu\nu} F^{\mu\nu}@f$.
     * @return @f$\delta g_{HAA}@f$
     */
    virtual double deltaG_hAA() const
    {
        return 0.0;
    };
    /**
     * @brief The full new physics contribution to the coupling of the effective interaction @f$H F_{\mu\nu} F^{\mu\nu}@f$,
     * including new local terms and modifications on the SM-loops. Normalized to the SM value.
     * @return @f$\delta g_{HAA}/g_{HAA}^SM}@f$
     */
    virtual double deltaG_hAARatio() const
    {
        return 0.0;
    };
    /**
     * @brief The new physics contribution to the coupling of the effective interaction @f$H f\bar{f}@f$.
     * @param[in] p a lepton or quark
     * @return @f$\delta g_{Hff}@f$
     */
    // no generation mixing
    virtual gslpp::complex deltaG_hff(const Particle p) const
    {
        return 0.0;
    };
    /**
     * @brief The new physics contribution to the Higgs self-coupling @f$ H H H@f$. Normalized to the SM value.
     * @return @f$\delta g_{HHH}/g_{HHH}^SM}@f$
     */
    virtual double deltaG_hhhRatio() const
    {
        return 0.0;
    };

    ////////////////////////////////////////////////////////////////////////

    /**
     * @brief The ratio @f$\mu_{ggH}@f$ between the gluon-gluon fusion Higgs
     * production cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ggH}@f$
     */
    virtual double muggH(const double sqrt_s) const
    {
        return 1.0;
    }
    
    /**
     * @brief The ratio @f$\mu_{ggHH}@f$ between the gluon-gluon fusion di-Higgs
     * production cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ggHH}@f$
     */
    virtual double muggHH(const double sqrt_s) const
    {
        return 1.0;
    }

    /**
     * @brief The ratio @f$\mu_{VBF}@f$ between the vector-boson fusion Higgs
     * production cross-section in the current model and in the Standard Model. 
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VBF}@f$
     */
    virtual double muVBF(const double sqrt_s) const
    {
        return 1.0;
    }
    
    /**
     * @brief The ratio @f$\mu_{VBF+\gamma}@f$ between the vector-boson fusion Higgs
     * production cross-section in association with a hard photon in the current model
     * and in the Standard Model. 
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VBF+\gamma}@f$
     */
    virtual double muVBFgamma(const double sqrt_s) const
    {
        return 1.0;
    }
    
    /**
     * @brief The ratio @f$\mu_{eeWBF}@f$ between the 
     * @f$ e^{+}e^{-}\to \nu\bar{\nu} H @f$ production
     * cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{eeWBF}@f$
     */
    virtual double mueeWBF(const double sqrt_s) const
    {
        return 1.0;
    }
    
    /**
     * @brief The ratio @f$\mu_{eeWBF}@f$ between the 
     * @f$ e^{+}e^{-}\to \nu\bar{\nu} H @f$ production
     * cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV, Pol_em and Pol_ep
     * are the polarization of electrons and positrons, respectively
     * @return @f$\mu_{eeWBF}@f$
     */
    virtual double mueeWBFPol(const double sqrt_s, const double Pol_em, const double Pol_ep) const
    {
        return 1.0;
    }
    
    /**
     * @brief The ratio @f$\mu_{e^+e^- \to H\nu\bar{\nu}}@f$ between the 
     * @f$ e^+e^- \to H\nu\bar{\nu} @f$ associated production
     * cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{e^+e^- \to H\nu\bar{\nu}}@f$
     */
    virtual double mueeHvv(const double sqrt_s) const
    {
        return 1.0;
    }
    
    /**
     * @brief The ratio @f$\mu_{e^+e^- \to H\nu\bar{\nu}}@f$ between the 
     * @f$ e^+e^- \to H\nu\bar{\nu} @f$ associated production
     * cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV, Pol_em and Pol_ep
     * are the polarization of electrons and positrons, respectively 
     * @return @f$\mu_{e^+e^- \to H\nu\bar{\nu}}@f$
     */
    virtual double mueeHvvPol(const double sqrt_s, const double Pol_em, const double Pol_ep) const
    {
        return 1.0;
    }
    
    /**
     * @brief The ratio @f$\mu_{eeZBF}@f$ between the 
     * @f$ e^{+}e^{-}\to e^{+}e^{-} H @f$ production
     * cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{eeZBF}@f$
     */
    virtual double mueeZBF(const double sqrt_s) const
    {
        return 1.0;
    }
    
    /**
     * @brief The ratio @f$\mu_{eeZBF}@f$ between the 
     * @f$ e^{+}e^{-}\to e^{+}e^{-} H @f$ production
     * cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV, Pol_em and Pol_ep
     * are the polarization of electrons and positrons, respectively 
     * @return @f$\mu_{eeZBF}@f$
     */
    virtual double mueeZBFPol(const double sqrt_s, const double Pol_em, const double Pol_ep) const
    {
        return 1.0;
    }
    
    /**
     * @brief The ratio @f$\mu_{epWBF}@f$ between the 
     * @f$ e^{-} p\to \nu j H @f$ production
     * cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{epWBF}@f$
     */
    virtual double muepWBF(const double sqrt_s) const
    {
        return 1.0;
    }
    
    /**
     * @brief The ratio @f$\mu_{epZBF}@f$ between the 
     * @f$ e^{-} p\to e^{-} j H @f$ production
     * cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{epZBF}@f$
     */
    virtual double muepZBF(const double sqrt_s) const
    {
        return 1.0;
    }

    /**
     * @brief The ratio @f$\mu_{WH}@f$ between the W-Higgs associated production
     * cross-section in the current model and in the Standard Model. 
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{WH}@f$
     */
    virtual double muWH(const double sqrt_s) const
    {
        return 1.0;
    }
    
    /**
     * @brief The ratio @f$\mu_{WH}@f$ between the W-Higgs associated production
     * cross-section in the current model and in the Standard Model, with @f$p_{T,H}>250@f$ GeV.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{WH}@f$
     */
    virtual double muWHpT250(const double sqrt_s) const
    {
        return 1.0;
    }

    /**
     * @brief The ratio @f$\mu_{ZH}@f$ between the Z-Higgs associated production
     * cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ZH}@f$
     */
    virtual double muZH(const double sqrt_s) const
    {
        return 1.0;
    }
    
    /**
     * @brief The ratio @f$\mu_{ZH}@f$ between the Z-Higgs associated production
     * cross-section in the current model and in the Standard Model, with @f$p_{T,H}>250@f$ GeV.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ZH}@f$
     */
    virtual double muZHpT250(const double sqrt_s) const
    {
        return 1.0;
    }

    /**
     * @brief The ratio @f$\mu_{eeZH}@f$ between the 
     * @f$ e^{+}e^{-}\to ZH @f$ associated production
     * cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{eeZH}@f$
     */
    virtual double mueeZH(const double sqrt_s) const
    {
        return 1.0;
    }
    
    /**
     * @brief The ratio @f$\mu_{eeZH, Z \to e^+ e^-, \mu^+ \mu^-}@f$ between the 
     * @f$ e^{+}e^{-}\to ZH, Z \to e^+ e^-, \mu^+ \mu^- @f$ associated production
     * cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{eeZH, Z \to e^+ e^-, \mu^+ \mu^-}@f$
     */
    virtual double mueeZllH(const double sqrt_s) const
    {
        return 1.0;
    }
    
    /**
     * @brief The ratio @f$\mu_{eeZH, Z \to q \bar{q}}@f$ between the 
     * @f$ e^{+}e^{-}\to ZH, Z \to q \bar{q} @f$ associated production
     * cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{eeZH, Z \to q \bar{q}}@f$
     */
    virtual double mueeZqqH(const double sqrt_s) const
    {
        return 1.0;
    }
    
    /**
     * @brief The ratio @f$\mu_{eeZH}@f$ between the 
     * @f$ e^{+}e^{-}\to ZH @f$ associated production
     * cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV, Pol_em and Pol_ep
     * are the polarization of electrons and positrons, respectively 
     * @return @f$\mu_{eeZH}@f$
     */
    virtual double mueeZHPol(const double sqrt_s, const double Pol_em, const double Pol_ep) const
    {
        return 1.0;
    }
    
    /**
     * @brief The ratio @f$\mu_{eeZH, Z \to e^+ e^-, \mu^+ \mu^-}@f$ between the 
     * @f$ e^{+}e^{-}\to ZH, Z \to e^+ e^-, \mu^+ \mu^- @f$ associated production
     * cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV, Pol_em and Pol_ep
     * are the polarization of electrons and positrons, respectively 
     * @return @f$\mu_{eeZH, Z \to e^+ e^-, \mu^+ \mu^-}@f$
     */
    virtual double mueeZllHPol(const double sqrt_s, const double Pol_em, const double Pol_ep) const
    {
        return 1.0;
    }
    
    /**
     * @brief The ratio @f$\mu_{eeZH, Z \to q \bar{q}}@f$ between the 
     * @f$ e^{+}e^{-}\to ZH, Z \to q \bar{q} @f$ associated production
     * cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV, Pol_em and Pol_ep
     * are the polarization of electrons and positrons, respectively 
     * @return @f$\mu_{eeZH, Z \to q \bar{q}}@f$
     */
    virtual double mueeZqqHPol(const double sqrt_s, const double Pol_em, const double Pol_ep) const
    {
        return 1.0;
    }
    
    
    /**
     * @brief the angular parameter @f$a@f$ from 
     * @f$\mu_{e^+e^- \to ZH}@f$ (arXiv:1708.09079 [hep-ph]).
     * @param[in] sqrt_s the center-of-mass energy in TeV, Pol_em and Pol_ep
     * are the polarization of electrons and positrons, respectively 
     * @return @f$a_{eeZH}@f$
     */
    virtual double aPskPol(const double sqrt_s, const double Pol_em, const double Pol_ep) const
    {
        return 0.0;
    }
    
    
    /**
     * @brief the angular parameter @f$b@f$ from 
     * @f$\mu_{e^+e^- \to ZH}@f$ (arXiv:1708.09079 [hep-ph]).
     * @param[in] sqrt_s the center-of-mass energy in TeV, Pol_em and Pol_ep
     * are the polarization of electrons and positrons, respectively 
     * @return @f$b_{eeZH}@f$
     */
    virtual double bPskPol(const double sqrt_s, const double Pol_em, const double Pol_ep) const
    {
        return 0.0;
    }

    
    /**
     * @brief The ratio @f$\mu_{VH}@f$ between the WH+ZH associated production
     * cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VH}@f$
     */
    virtual double muVH(const double sqrt_s) const
    {
        return 1.0;
    }
    

    /**
     * @brief The ratio @f$\mu_{VH}@f$ between the WH+ZH associated production
     * cross-section in the current model and in the Standard Model, with @f$p_{T,H}>250@f$ GeV.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VH}@f$
     */
    virtual double muVHpT250(const double sqrt_s) const
    {
        return 1.0;
    }
    

    /**
     * @brief The ratio @f$\mu_{VBF+VH}@f$ between the sum of VBF and WH+ZH associated production
     * cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VBF+VH}@f$
     */
    virtual double muVBFpVH(const double sqrt_s) const
    {
        return 1.0;
    }

    /**
     * @brief The ratio @f$\mu_{ttH}@f$ between the t-tbar-Higgs associated 
     * production cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ttH}@f$
     */
    virtual double muttH(const double sqrt_s) const
    {
        return 1.0;
    }
    
    /**
     * @brief The ratio @f$\mu_{tH}@f$ between the t-Higgs associated 
     * production cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{tH}@f$
     */
    virtual double mutH(const double sqrt_s) const                              //AG:added
    {
        return 1.0;
    }
    
    /**
     * @brief The ratio @f$\mu_{tHq}@f$ between the t-q-Higgs associated 
     * production cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{tHq}@f$
     */
    virtual double mutHq(const double sqrt_s) const
    {
        return 1.0;
    }

    /**
     * @brief The ratio @f$\mu_{ggH+ttH}@f$ between the sum of gluon-gluon fusion
     * and t-tbar-Higgs associated 
     * production cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ggH+ttH}@f$
     */
    virtual double muggHpttH(const double sqrt_s) const
    {
        return 1.0;
    }
        
    /**
     * @brief The ratio @f$\mu_{eettH}@f$ between the 
     * @f$ e^{+}e^{-}\to t\bar{t} H @f$ production
     * cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{eettH}@f$
     */
    virtual double mueettH(const double sqrt_s) const
    {
        return 1.0;
    }
    
    /**
     * @brief The ratio @f$\mu_{eettH}@f$ between the 
     * @f$ e^{+}e^{-}\to t\bar{t} H @f$ production
     * cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV, Pol_em and Pol_ep
     * are the polarization of electrons and positrons, respectively
     * @return @f$\mu_{eettH}@f$
     */
    virtual double mueettHPol(const double sqrt_s, const double Pol_em, const double Pol_ep) const
    {
        return 1.0;
    }
    
    /**
     * @brief The ratio @f$\mu_{\mu\mu H}@f$ between the @f$\sigma(\mu \mu \to H)}@f$
     * production cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{\mu\mu H}@f$
     */
    virtual double mummH(const double sqrt_s) const
    {
        return 1.0;
    }
    
    /**
     * @brief The ratio @f$\mu_{\mu\mu H}@f$ between the @f$\sigma(\mu \mu \to H)}@f$
     * production cross-section in the current model and in the Standard Model, 
     * in the narrow width approximation.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{\mu\mu H}@f$
     */
    virtual double mummHNWA(const double sqrt_s) const
    {
        return 1.0;
    }
    
    
    /**
     * @brief The ratio @f$\mu_{\mu\mu ZH}@f$ between the @f$\sigma(\mu \mu \to Z H)}@f$
     * production cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{\mu\mu ZH}@f$
     */
    virtual double mummZH(const double sqrt_s) const
    {
        return 1.0;
    }
    
    
    /**
     * @brief The ratio @f$\mu_{\mu\mu H\nu\nu}@f$ between the @f$\sigma(\mu \mu \to H \nu \nu)}@f$
     * production cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{\mu\mu H\nu\nu}@f$
     */
    virtual double mummHvv(const double sqrt_s) const
    {
        return 1.0;
    }
    
    
    /**
     * @brief The ratio @f$\mu_{\mu\mu H\mu\mu}@f$ between the @f$\sigma(\mu \mu \to H \mu \mu)}@f$
     * production cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{\mu\mu H\mu\mu}@f$
     */
    virtual double mummHmm(const double sqrt_s) const
    {
        return 1.0;
    }
    
    
    /**
     * @brief The ratio @f$\mu_{\mu\mu ttH}@f$ between the @f$\sigma(\mu \mu \to t\bar{t} H )}@f$
     * production cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{\mu\mu ttH}@f$
     */
    virtual double mummttH(const double sqrt_s) const
    {
        return 1.0;
    }
    
    
    /**
     * @brief The ratio @f$\sigma(ttH)/\sigma(ttZ)@f$ 
     * in the @f$H,Z\to b\bar{b}@f$ channel in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\sigma(ttH)/\sigma(ttZ)@f$ normalized to the SM
     */
    virtual double muttHZbbboost(const double sqrt_s) const
    {
        return 1.0;
    }
        
    virtual double muggHgaga(const double sqrt_s) const
    {
        return 1.0;
    }
    
    /**
     * @brief The ratio @f$\mu_{ggH,\gamma\gamma}@f$ between the gluon-gluon fusion Higgs
     * production cross-section with subsequent decay into 2 photons in the
     * current model and in the Standard Model. Includes interference effects
     * with the background, following arXiv:1704.08259
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ggH,\gamma\gamma}@f$
     */
    virtual double muggHgagaInt(const double sqrt_s) const
    {
        return 1.0;        
    };
    
    virtual double muVBFHgaga(const double sqrt_s) const
    {
        return 1.0;
    }
    virtual double muZHgaga(const double sqrt_s) const
    {
        return 1.0;
    }
    virtual double muWHgaga(const double sqrt_s) const
    {
        return 1.0;
    }
    virtual double muVHgaga(const double sqrt_s) const
    {
        return 1.0;
    }
    virtual double muttHgaga(const double sqrt_s) const
    {
        return 1.0;
    }   
    virtual double muggHZga(const double sqrt_s) const
    {
        return 1.0;
    }
    virtual double muVBFHZga(const double sqrt_s) const
    {
        return 1.0;
    }
    virtual double muZHZga(const double sqrt_s) const
    {
        return 1.0;
    }    
    virtual double muWHZga(const double sqrt_s) const
    {
        return 1.0;
    }
    virtual double muVHZga(const double sqrt_s) const
    {
        return 1.0;
    }
    virtual double muttHZga(const double sqrt_s) const
    {
        return 1.0;
    }  
    virtual double muggHZZ(const double sqrt_s) const
    {
        return 1.0;
    }
    virtual double muVBFHZZ(const double sqrt_s) const
    {
        return 1.0;
    }
    virtual double muZHZZ(const double sqrt_s) const
    {
        return 1.0;
    }    
    virtual double muWHZZ(const double sqrt_s) const
    {
        return 1.0;
    }
    virtual double muVHZZ(const double sqrt_s) const
    {
        return 1.0;
    }
    virtual double muttHZZ(const double sqrt_s) const
    {
        return 1.0;
    }
    
    virtual double muggHZZ4l(const double sqrt_s) const
    {
        return 1.0;
    }
    virtual double muVBFHZZ4l(const double sqrt_s) const
    {
        return 1.0;
    }
    virtual double muZHZZ4l(const double sqrt_s) const
    {
        return 1.0;
    }    
    virtual double muWHZZ4l(const double sqrt_s) const
    {
        return 1.0;
    }
    virtual double muVHZZ4l(const double sqrt_s) const
    {
        return 1.0;
    }
    virtual double muttHZZ4l(const double sqrt_s) const
    {
        return 1.0;
    }
    
    virtual double muggHWW(const double sqrt_s) const
    {
        return 1.0;
    }
    virtual double muVBFHWW(const double sqrt_s) const
    {
        return 1.0;
    }
    virtual double muZHWW(const double sqrt_s) const
    {
        return 1.0;
    }
    virtual double muWHWW(const double sqrt_s) const
    {
        return 1.0;
    }
    virtual double muVHWW(const double sqrt_s) const
    {
        return 1.0;
    }
    virtual double muttHWW(const double sqrt_s) const
    {
        return 1.0;
    }
    virtual double muggHWW2l2v(const double sqrt_s) const
    {
        return 1.0;
    }
    virtual double muVBFHWW2l2v(const double sqrt_s) const
    {
        return 1.0;
    }
    virtual double muZHWW2l2v(const double sqrt_s) const
    {
        return 1.0;
    }
    virtual double muWHWW2l2v(const double sqrt_s) const
    {
        return 1.0;
    }
    virtual double muVHWW2l2v(const double sqrt_s) const
    {
        return 1.0;
    }
    virtual double muttHWW2l2v(const double sqrt_s) const
    {
        return 1.0;
    }   
    virtual double muggHmumu(const double sqrt_s) const
    {
        return 1.0;
    }
    virtual double muVBFHmumu(const double sqrt_s) const
    {
        return 1.0;
    }
    virtual double muZHmumu(const double sqrt_s) const
    {
        return 1.0;
    }
    virtual double muWHmumu(const double sqrt_s) const
    {
        return 1.0;
    }
    virtual double muVHmumu(const double sqrt_s) const
    {
        return 1.0;
    }
    virtual double muttHmumu(const double sqrt_s) const
    {
        return 1.0;
    }
    virtual double muggHtautau(const double sqrt_s) const
    {
        return 1.0;
    }
    virtual double muVBFHtautau(const double sqrt_s) const
    {
        return 1.0;
    }
    virtual double muZHtautau(const double sqrt_s) const
    {
        return 1.0;
    }
    virtual double muWHtautau(const double sqrt_s) const
    {
        return 1.0;
    }
    virtual double muVHtautau(const double sqrt_s) const
    {
        return 1.0;
    }
    virtual double muttHtautau(const double sqrt_s) const
    {
        return 1.0;
    }
    virtual double muggHbb(const double sqrt_s) const
    {
        return 1.0;
    }
    virtual double muVBFHbb(const double sqrt_s) const
    {
        return 1.0;
    }
    virtual double muZHbb(const double sqrt_s) const
    {
        return 1.0;
    }
    virtual double muWHbb(const double sqrt_s) const
    {
        return 1.0;
    }
    virtual double muVHbb(const double sqrt_s) const
    {
        return 1.0;
    }
    virtual double muttHbb(const double sqrt_s) const
    {
        return 1.0;
    }
    virtual double muppHmumu(const double sqrt_s) const
    {
        return 1.0;
    }
    virtual double muppHZga(const double sqrt_s) const
    {
        return 1.0;
    }
        
    //AG:begin,begin
    virtual double delta_muggH_1(const double sqrt_s) const 
    { 
        return 0.0 ;
    }
    virtual double delta_muggH_2(const double sqrt_s) const 
    { 
        return 0.0 ;
    }
    
    virtual double delta_muVBF_1(const double sqrt_s) const 
    { 
        return 0.0 ;
    }
    virtual double delta_muVBF_2(const double sqrt_s) const 
    { 
        return 0.0 ;
    }
    
    virtual double delta_muWH_1(const double sqrt_s) const 
    { 
        return 0.0 ;
    }
    virtual double delta_muWH_2(const double sqrt_s) const 
    { 
        return 0.0 ;
    }
    
    virtual double delta_muZH_1(const double sqrt_s) const 
    { 
        return 0.0 ;
    }
    virtual double delta_muZH_2(const double sqrt_s) const 
    { 
        return 0.0 ;
    }
    
    virtual double delta_muVH_1(const double sqrt_s) const 
    { 
        return 0.0 ;
    }
    virtual double delta_muVH_2(const double sqrt_s) const 
    { 
        return 0.0 ;
    }
    
    virtual double delta_muttH_1(const double sqrt_s) const 
    { 
        return 0.0 ;
    }
    virtual double delta_muttH_2(const double sqrt_s) const 
    { 
        return 0.0 ;
    }
    
    virtual double delta_mutH_1(const double sqrt_s) const 
    { 
        return 0.0 ;
    }
    virtual double delta_mutH_2(const double sqrt_s) const 
    { 
        return 0.0 ;
    }
    ///////////////////////////////////////////////////////
    virtual double  deltaGammaHggRatio1() const 
    {
        return 0.0;
    }
    virtual double  deltaGammaHggRatio2() const 
    {
        return 0.0;
    }
    
    virtual double  deltaGammaHWWRatio1() const 
    {
        return 0.0;
    }
    virtual double  deltaGammaHWWRatio2() const 
    {
        return 0.0;
    }
    
    virtual double  deltaGammaHZZRatio1() const 
    {
        return 0.0;
    }
    virtual double  deltaGammaHZZRatio2() const 
    {
        return 0.0;
    }
    
    virtual double  deltaGammaHgagaRatio1() const 
    {
        return 0.0;
    }
    virtual double  deltaGammaHgagaRatio2() const 
    {
        return 0.0;
    }
    
    virtual double  deltaGammaHZgaRatio1() const 
    {
        return 0.0;
    }
    virtual double  deltaGammaHZgaRatio2() const 
    {
        return 0.0;
    }
    
    virtual double  deltaGammaHbbRatio1() const 
    {
        return 0.0;
    }
    virtual double  deltaGammaHbbRatio2() const 
    {
        return 0.0;
    }
    
    virtual double  deltaGammaHmumuRatio1() const 
    {
        return 0.0;
    }
    virtual double  deltaGammaHmumuRatio2() const 
    {
        return 0.0;
    }
    
    virtual double  deltaGammaHtautauRatio1() const 
    {
        return 0.0;
    }
    virtual double  deltaGammaHtautauRatio2() const 
    {
        return 0.0;
    }
    
    virtual double  deltaGammaHccRatio1() const 
    {
        return 0.0;
    }
    virtual double  deltaGammaHccRatio2() const 
    {
        return 0.0;
    }
    
    virtual double  deltaGammaTotalRatio1() const 
    {
        return 0.0;
    }
    virtual double  deltaGammaTotalRatio2() const 
    {
        return 0.0;
    }

    //AG:end
    
    ///////////////////////HIGGS DECAY WIDTHS/////////////////////////

    /**
     * @brief The ratio of the @f$\Gamma(H)@f$ in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H)@f$/@f$\Gamma(H)_{\mathrm{SM}}@f$
     */
    virtual double computeGammaTotalRatio() const
    {
        return 1.0;
    }
    
    ///////////////////////HIGGS BRANCHING RATIOS/////////////////////////
    
    /**
     * @brief The ratio of the Br@f$(H\to gg)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to gg)@f$/Br@f$(H\to gg)_{\mathrm{SM}}@f$
     */
    virtual double BrHggRatio() const
    {
        return 1.0;
    }

    /**
     * @brief The ratio of the Br@f$(H\to WW)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to WW)@f$/Br@f$(H\to WW)_{\mathrm{SM}}@f$
     */
    virtual double BrHWWRatio() const
    {
        return 1.0;
    }
    
    /**
     * @brief The ratio of the Br@f$(H\to ZZ)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to ZZ)@f$/Br@f$(H\to ZZ)_{\mathrm{SM}}@f$
     */
    virtual double BrHZZRatio() const
    {
        return 1.0;
    }
    
    /**
     * @brief The ratio of the Br@f$(H\to VV)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to VV)@f$/Br@f$(H\to VV)_{\mathrm{SM}}@f$
     */
    virtual double BrHVVRatio() const
    {
        return 1.0;
    }

    /**
     * @brief The ratio of the Br@f$(H\to Z\gamma)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to Z\gamma)@f$/Br@f$(H\to Z\gamma)_{\mathrm{SM}}@f$
     */
    virtual double BrHZgaRatio() const
    {
        return 1.0;
    }
    
    /**
     * @brief The ratio of the Br@f$(H\to Z\gamma\to ll\gamma)@f$ (@f$l=e,\mu @f$) in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to Z\gamma\to ll\gamma)@f$/Br@f$(H\to Z\gamma\to ll\gamma)_{\mathrm{SM}}@f$
     */
    virtual double BrHZgallRatio() const
    {
        return 1.0;
    }
    
    /**
     * @brief The ratio of the Br@f$(H\to Z\gamma\to ee\gamma)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to Z\gamma\to ee\gamma)@f$/Br@f$(H\to Z\gamma\to ee\gamma)_{\mathrm{SM}}@f$
     */
    virtual double BrHZgaeeRatio() const
    {
        return 1.0;
    }
    
    /**
     * @brief The ratio of the Br@f$(H\to Z\gamma\to \mu\mu\gamma)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to Z\gamma\to \mu\mu\gamma)@f$/Br@f$(H\to Z\gamma\to \mu\mu\gamma)_{\mathrm{SM}}@f$
     */
    virtual double BrHZgamumuRatio() const
    {
        return 1.0;
    }

    /**
     * @brief The ratio of the Br@f$(H\to \gamma\gamma)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to \gamma\gamma)@f$/Br@f$(H\to \gamma\gamma)_{\mathrm{SM}}@f$
     */
    virtual double BrHgagaRatio() const
    {
        return 1.0;
    }
    
    /////////////////////// HIGGS TO 2 FERMION DECAYS /////////////////////////

    
    /**
     * @brief The ratio of the Br@f$(H\to \mu^+\mu^-)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to \mu^+\mu^-)@f$/Br@f$(H\to \mu^+\mu^-)_{\mathrm{SM}}@f$
     */
    virtual double BrHmumuRatio() const
    {
        return 1.0;
    }

    /**
     * @brief The ratio of the Br@f$(H\to \tau^+\tau^-)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to \tau^+\tau^-)@f$/Br@f$(H\to \tau^+\tau^-)_{\mathrm{SM}}@f$
     */
    virtual double BrHtautauRatio() const
    {
        return 1.0;
    }

    /**
     * @brief The ratio of the Br@f$(H\to c\bar{c})@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to c\bar{c})@f$/Br@f$(H\to c\bar{c})_{\mathrm{SM}}@f$
     */
    virtual double BrHccRatio() const
    {
        return 1.0;
    }

    /**
     * @brief The ratio of the Br@f$(H\to b\bar{b})@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to b\bar{b})@f$/Br@f$(H\to b\bar{b})_{\mathrm{SM}}@f$
     */
    virtual double BrHbbRatio() const
    {
        return 1.0;
    }
    

    /////////////////////// HIGGS TO 4 FERMION DECAYS /////////////////////////
    
    /**
     * @brief The ratio of the Br@f$(H\to 2L2L')@f$ (@f$L,L'=e,\mu,\tau@f$) in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to 2L2L')@f$/Br@f$(H\to 2L2L')_{\mathrm{SM}}@f$
     */
    virtual double BrH2L2LRatio() const
    {
        return 1.0;
    }
    /**
     * @brief The ratio of the Br@f$(H\to 2e 2\mu)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to 2L2L)@f$/Br@f$(H\to 2e 2\mu)_{\mathrm{SM}}@f$
     */
    virtual double BrH2e2muRatio() const
    {
        return 1.0;
    }
    /**
     * @brief The ratio of the Br@f$(H\to 2v2v)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to 2v2v)@f$/Br@f$(H\to 2v2v)_{\mathrm{SM}}@f$
     */
    virtual double BrH2v2vRatio() const
    {
        return 1.0;
    }
    /**
     * @brief The ratio of the Br@f$(H\to 2L2v)@f$ (@f$L=e,\mu,\tau@f$) in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to 2L2v)@f$/Br@f$(H\to 2L2v)_{\mathrm{SM}}@f$
     */
    virtual double BrH2L2vRatio() const
    {
        return 1.0;
    }
    /**
     * @brief The ratio of the Br@f$(H\to 2L2v)@f$ (@f$L=e,\mu@f$) in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to 2L2v)@f$/Br@f$(H\to 2L2v)_{\mathrm{SM}}@f$
     */
    virtual double BrH2L2v2Ratio() const
    {
        return 1.0;
    }
    /**
     * @brief The ratio of the Br@f$(H\to 2e2v)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to 2e2v)@f$/Br@f$(H\to 2e2v)_{\mathrm{SM}}@f$
     */
    virtual double BrH2e2vRatio() const
    {
        return 1.0;
    }
    /**
     * @brief The ratio of the Br@f$(H\to 2\mu 2v)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to 2\mu 2v)@f$/Br@f$(H\to 2\mu 2v)_{\mathrm{SM}}@f$
     */
    virtual double BrH2mu2vRatio() const
    {
        return 1.0;
    }
    /**
     * @brief The ratio of the Br@f$(H\to 2u2u)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to 2u2u)@f$/Br@f$(H\to 2u2u)_{\mathrm{SM}}@f$
     */
    virtual double BrH2u2uRatio() const
    {
        return 1.0;
    }
    /**
     * @brief The ratio of the Br@f$(H\to 2d2d)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to 2d2d)@f$/Br@f$(H\to 2d2d)_{\mathrm{SM}}@f$
     */
    virtual double BrH2d2dRatio() const
    {
        return 1.0;
    }
    /**
     * @brief The ratio of the Br@f$(H\to 2u2d)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to 2u2d)@f$/Br@f$(H\to 2u2d)_{\mathrm{SM}}@f$
     */
    virtual double BrH2u2dRatio() const
    {
        return 1.0;
    }
    /**
     * @brief The ratio of the Br@f$(H\to 2L2u)@f$ (@f$L=e,\mu,\tau@f$) in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to 2L2u)@f$/Br@f$(H\to 2L2u)_{\mathrm{SM}}@f$
     */
    virtual double BrH2L2uRatio() const
    {
        return 1.0;
    }
    /**
     * @brief The ratio of the Br@f$(H\to 2L2d)@f$ (@f$L=e,\mu,\tau@f$) in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to 2L2d)@f$/Br@f$(H\to 2L2d)_{\mathrm{SM}}@f$
     */
    virtual double BrH2L2dRatio() const
    {
        return 1.0;
    }
    /**
     * @brief The ratio of the Br@f$(H\to 2v2u)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to 2v2u)@f$/Br@f$(H\to 2v2u)_{\mathrm{SM}}@f$
     */
    virtual double BrH2v2uRatio() const
    {
        return 1.0;
    }
    /**
     * @brief The ratio of the Br@f$(H\to 2v2d)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to 2v2d)@f$/Br@f$(H\to 2v2d)_{\mathrm{SM}}@f$
     */
    virtual double BrH2v2dRatio() const
    {
        return 1.0;
    }
    /**
     * @brief The ratio of the Br@f$(H\to 4L)@f$ (@f$L=e,\mu,\tau@f$) in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to 4L)@f$/Br@f$(H\to 4L)_{\mathrm{SM}}@f$
     */
    virtual double BrH4LRatio() const
    {
        return 1.0;
    }
    /**
     * @brief The ratio of the Br@f$(H\to 4L)@f$ (@f$L=e,\mu@f$) in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to 4L)@f$/Br@f$(H\to 4L)_{\mathrm{SM}}@f$
     */
    virtual double BrH4L2Ratio() const
    {
        return 1.0;
    }
    /**
     * @brief The ratio of the Br@f$(H\to 4e)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to 4e)@f$/Br@f$(H\to 4e)_{\mathrm{SM}}@f$
     */
    virtual double BrH4eRatio() const
    {
        return 1.0;
    }
    /**
     * @brief The ratio of the Br@f$(H\to 4\mu)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to 4\mu)@f$/Br@f$(H\to 4\mu)_{\mathrm{SM}}@f$
     */
    virtual double BrH4muRatio() const
    {
        return 1.0;
    }
    /**
     * @brief The ratio of the Br@f$(H\to 4v)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to 4v)@f$/Br@f$(H\to 4v)_{\mathrm{SM}}@f$
     */
    virtual double BrH4vRatio() const
    {
        return 1.0;
    }
    /**
     * @brief The ratio of the Br@f$(H\to 4u)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to 4u)@f$/Br@f$(H\to 4u)_{\mathrm{SM}}@f$
     */
    virtual double BrH4uRatio() const
    {
        return 1.0;
    }   
    /**
     * @brief The ratio of the Br@f$(H\to 4d)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to 4d)@f$/Br@f$(H\to 4d)_{\mathrm{SM}}@f$
     */
    virtual double BrH4dRatio() const
    {
        return 1.0;
    }    
    /**
     * @brief The ratio of the Br@f$(H\to LvvL)@f$ (@f$L=e,\mu,\tau@f$) in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to LvvL)@f$/Br@f$(H\to LvvL)_{\mathrm{SM}}@f$
     */
    virtual double BrHLvvLRatio() const
    {
        return 1.0;
    }
    /**
     * @brief The ratio of the Br@f$(H\to e\nu \mu\nu)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to e\nu \mu\nu)@f$/Br@f$(H\to e\nu \mu\nu)_{\mathrm{SM}}@f$
     */
    virtual double BrHevmuvRatio() const
    {
        return 1.0;
    }
    /**
     * @brief The ratio of the Br@f$(H\to uddu)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to uddu)@f$/Br@f$(H\to uddu)_{\mathrm{SM}}@f$
     */
    virtual double BrHudduRatio() const
    {
        return 1.0;
    }
    /**
     * @brief The ratio of the Br@f$(H\to Lvud)@f$ (@f$L=e,\mu,\tau@f$) in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to Lvud)@f$/Br@f$(H\to Lvud)_{\mathrm{SM}}@f$
     */
    virtual double BrHLvudRatio() const
    {
        return 1.0;
    }
    /**
     * @brief The ratio of the Br@f$(H\to 2ud)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to 2ud)@f$/Br@f$(H\to 2ud)_{\mathrm{SM}}@f$
     */
    virtual double BrH2udRatio() const
    {
        return 1.0;
    }
    /**
     * @brief The ratio of the Br@f$(H\to 2Lv)@f$ (@f$L=e,\mu,\tau@f$) in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to 2Lv)@f$/Br@f$(H\to 2Lv)_{\mathrm{SM}}@f$
     */
    virtual double BrH2LvRatio() const
    {
        return 1.0;
    }
    /**
     * @brief The ratio of the Br@f$(H\to 2Lv)@f$ (@f$L=e,\mu@f$) in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to 2Lv)@f$/Br@f$(H\to 2Lv)_{\mathrm{SM}}@f$
     */
    virtual double BrH2Lv2Ratio() const
    {
        return 1.0;
    }
    /**
     * @brief The ratio of the Br@f$(H\to 2ev)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to 2ev)@f$/Br@f$(H\to 2ev)_{\mathrm{SM}}@f$
     */
    virtual double BrH2evRatio() const
    {
        return 1.0;
    }
    /**
     * @brief The ratio of the Br@f$(H\to 2ev)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to 2\mu v)@f$/Br@f$(H\to 2\mu v)_{\mathrm{SM}}@f$
     */
    virtual double BrH2muvRatio() const
    {
        return 1.0;
    }
    /**
     * @brief The ratio of the Br@f$(H\to 4f)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to 4f)@f$/Br@f$(H\to 4f)_{\mathrm{SM}}@f$
     */
    virtual double BrH4fRatio() const
    {
        return 1.0;
    }
    
    // DECAYS INVOLVING ONLY ELECTRONS, MUONS OR NEUTRINOS IN THE FINAL STATES 
    
    /**
     * @brief The ratio of the Br@f$(H\to 4l)@f$ (@f$l=e,\mu@f$) in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to 4l)@f$/Br@f$(H\to 4l)_{\mathrm{SM}}@f$
     */
    virtual double BrH4lRatio() const
    {
        return 1.0;
    }
    /**
     * @brief The ratio of the Br@f$(H\to 2l2v)@f$ (@f$l=e,\mu@f$) in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to 2l2v)@f$/Br@f$(H\to 2l2v)_{\mathrm{SM}}@f$
     */
    virtual double BrH2l2vRatio() const
    {
        return 1.0;
    }

    
    ///////////////////////OTHER DEDICATED (SEMI-)LEPTONIC 4 FERMION DECAYS/////////////////////////
    
    /**
     * @brief The ratio of the Br@f$(H\to l l j j)@f$ (@f$l=e,\mu,~~j\not=b@f$) in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to l l j j)@f$/Br@f$(H\to l l j j)_{\mathrm{SM}}@f$
     */
    virtual double BrHlljjRatio() const
    {
        return 1.0;
    }
    /**
     * @brief The ratio of the Br@f$(H\to l \nu j j)@f$ (@f$l=e,\mu,~~j\not=b@f$) in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to l \nu j j)@f$/Br@f$(H\to l \nu j j)_{\mathrm{SM}}@f$
     */
    virtual double BrHlvjjRatio() const
    {
        return 1.0;
    }
    /**
     * @brief The ratio of the Br@f$(H\to l \nu l \nu, l \nu j j)@f$ (@f$l=e,\mu,~~j\not=b@f$) in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to l \nu l \nu, l \nu j j)@f$/Br@f$(H\to l \nu l \nu, l \nu j j)_{\mathrm{SM}}@f$
     */
    virtual double BrHlv_lvorjjRatio() const
    {
        return 1.0;
    }
    /**
     * @brief The ratio of the Br@f$(H\to l l \nu\nu, l l j j)@f$ (@f$l=e,\mu,~~j\not=b@f$) in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to l l \nu\nu, l l j j)@f$/Br@f$(H\to l l \nu\nu, l l j j)_{\mathrm{SM}}@f$
     */
    virtual double BrHll_vvorjjRatio() const
    {
        return 1.0;
    }
    
    ///////////////////////OTHER HIGGS BRANCHING RATIOS/////////////////////////    
    
    /**
     * @brief The branching ratio of the of the Higgs into exotic particles.
     * @return Br@f$(H\to exotic)@f$
     */
    virtual double Br_H_exo() const
    {
        return 0.0;
    };
    
    /**
     * @brief The branching ratio of the of the Higgs into invisible particles.
     * @return Br@f$(H\to invisible)@f$
     */
    virtual double Br_H_inv() const
    {
        return 0.0;
    };
    
    /**
     * @brief The branching ratio of the of the Higgs into invisible particles 
     * (only invisible new particles).
     * @return Br@f$(H\to invisible,NP)@f$
     */
    virtual double Br_H_inv_NP() const
    {
        return 0.0;
    };
    
    /**
     * @brief The ratio of the Br@f$(H\to visible)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to visible)@f$/Br@f$(H\to visible)_{\mathrm{SM}}@f$
     */
    virtual double BrHvisRatio() const
    {
        return 1.0;
    }
    
    /**
     * @brief The ratio of the Br@f$(H\to invisible)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to invisible)@f$/Br@f$(H\to ZZ \to invisible)_{\mathrm{SM}}@f$
     */
    virtual double BrHtoinvRatio() const
    {
        return 1.0;
    }

    /**
     * @brief 
     * @return 
     */
    virtual double UpperLimitZgammaA13(const double sqrt_s) const
    {
        return 1.0;
    }

    /**
     * @brief 
     * @return 
     */
    virtual double UpperLimitZgammaC13(const double sqrt_s) const
    {
        return 1.0;
    }

    /**
     * @brief 
     * @return 
     */
    virtual double UpperLimitZgammaA(const double sqrt_s) const
    {
        return 1.0;
    }

    /**
     * @brief 
     * @return 
     */
    virtual double UpperLimitZgammaC(const double sqrt_s) const
    {
        return 1.0;
    }

    /**
     * @brief 
     * @return 
     */
    virtual double cgplusct() const
    {
        return 1.0;
    }

    /**
     * @brief 
     * @return 
     */
    virtual double cgaplusct() const
    {
        return 1.0;
    }

    /**
     * @brief 
     * @return 
     */
    virtual double cgminuscga() const
    {
        return 0.0;
    }

    /**
     * @brief 
     * @return 
     */
    virtual double cVpluscb() const
    {
        return 2.0;
    }

    /**
     * @brief 
     * @return 
     */
    virtual double cVplusctau() const
    {
        return 2.0;
    }

    /**
     * @brief 
     * @return 
     */
    virtual double cbminuscc() const
    {
        return 0.0;
    }

    /**
     * @brief 
     * @return 
     */
    virtual double cbminusctau() const
    {
        return 0.0;
    }

    /**
     * @brief 
     * @return 
     */
    virtual double ccminusctau() const
    {
        return 0.0;
    }
    
    
    
////////////////////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------------------
//-- Special Hadron collider signal strengths with separate full TH unc U(prod x decay) ---
//-----------------------------------------------------------------------------------------
////////////////////////////////////////////////////////////////////////////////////////////
    
    virtual double muTHUggHgaga(const double sqrt_s) const
    {
        return 1.0;
    }
        
    virtual double muTHUVBFHgaga(const double sqrt_s) const
    {
        return 1.0;
    }
    virtual double muTHUZHgaga(const double sqrt_s) const
    {
        return 1.0;
    }
    virtual double muTHUWHgaga(const double sqrt_s) const
    {
        return 1.0;
    }
    virtual double muTHUVHgaga(const double sqrt_s) const
    {
        return 1.0;
    }
    virtual double muTHUttHgaga(const double sqrt_s) const
    {
        return 1.0;
    }   
    virtual double muTHUggHZga(const double sqrt_s) const
    {
        return 1.0;
    }
    virtual double muTHUVBFHZga(const double sqrt_s) const
    {
        return 1.0;
    }
    virtual double muTHUZHZga(const double sqrt_s) const
    {
        return 1.0;
    }    
    virtual double muTHUWHZga(const double sqrt_s) const
    {
        return 1.0;
    }
    virtual double muTHUVHZga(const double sqrt_s) const
    {
        return 1.0;
    }
    virtual double muTHUttHZga(const double sqrt_s) const
    {
        return 1.0;
    }  
    virtual double muTHUggHZZ(const double sqrt_s) const
    {
        return 1.0;
    }
    virtual double muTHUVBFHZZ(const double sqrt_s) const
    {
        return 1.0;
    }
    virtual double muTHUZHZZ(const double sqrt_s) const
    {
        return 1.0;
    }    
    virtual double muTHUWHZZ(const double sqrt_s) const
    {
        return 1.0;
    }
    virtual double muTHUVHZZ(const double sqrt_s) const
    {
        return 1.0;
    }
    virtual double muTHUttHZZ(const double sqrt_s) const
    {
        return 1.0;
    }
    
    virtual double muTHUggHZZ4l(const double sqrt_s) const
    {
        return 1.0;
    }
    virtual double muTHUVBFHZZ4l(const double sqrt_s) const
    {
        return 1.0;
    }
    virtual double muTHUZHZZ4l(const double sqrt_s) const
    {
        return 1.0;
    }    
    virtual double muTHUWHZZ4l(const double sqrt_s) const
    {
        return 1.0;
    }
    virtual double muTHUVHZZ4l(const double sqrt_s) const
    {
        return 1.0;
    }
    virtual double muTHUttHZZ4l(const double sqrt_s) const
    {
        return 1.0;
    }
    
    virtual double muTHUggHWW(const double sqrt_s) const
    {
        return 1.0;
    }
    virtual double muTHUVBFHWW(const double sqrt_s) const
    {
        return 1.0;
    }
    virtual double muTHUZHWW(const double sqrt_s) const
    {
        return 1.0;
    }
    virtual double muTHUWHWW(const double sqrt_s) const
    {
        return 1.0;
    }
    virtual double muTHUVHWW(const double sqrt_s) const
    {
        return 1.0;
    }
    virtual double muTHUttHWW(const double sqrt_s) const
    {
        return 1.0;
    }
    virtual double muTHUggHWW2l2v(const double sqrt_s) const
    {
        return 1.0;
    }
    virtual double muTHUVBFHWW2l2v(const double sqrt_s) const
    {
        return 1.0;
    }
    virtual double muTHUZHWW2l2v(const double sqrt_s) const
    {
        return 1.0;
    }
    virtual double muTHUWHWW2l2v(const double sqrt_s) const
    {
        return 1.0;
    }
    virtual double muTHUVHWW2l2v(const double sqrt_s) const
    {
        return 1.0;
    }
    virtual double muTHUttHWW2l2v(const double sqrt_s) const
    {
        return 1.0;
    }   
    virtual double muTHUggHmumu(const double sqrt_s) const
    {
        return 1.0;
    }
    virtual double muTHUVBFHmumu(const double sqrt_s) const
    {
        return 1.0;
    }
    virtual double muTHUZHmumu(const double sqrt_s) const
    {
        return 1.0;
    }
    virtual double muTHUWHmumu(const double sqrt_s) const
    {
        return 1.0;
    }
    virtual double muTHUVHmumu(const double sqrt_s) const
    {
        return 1.0;
    }
    virtual double muTHUttHmumu(const double sqrt_s) const
    {
        return 1.0;
    }
    virtual double muTHUggHtautau(const double sqrt_s) const
    {
        return 1.0;
    }
    virtual double muTHUVBFHtautau(const double sqrt_s) const
    {
        return 1.0;
    }
    virtual double muTHUZHtautau(const double sqrt_s) const
    {
        return 1.0;
    }
    virtual double muTHUWHtautau(const double sqrt_s) const
    {
        return 1.0;
    }
    virtual double muTHUVHtautau(const double sqrt_s) const
    {
        return 1.0;
    }
    virtual double muTHUttHtautau(const double sqrt_s) const
    {
        return 1.0;
    }
    virtual double muTHUggHbb(const double sqrt_s) const
    {
        return 1.0;
    }
    virtual double muTHUVBFHbb(const double sqrt_s) const
    {
        return 1.0;
    }
    virtual double muTHUZHbb(const double sqrt_s) const
    {
        return 1.0;
    }
    virtual double muTHUWHbb(const double sqrt_s) const
    {
        return 1.0;
    }
    virtual double muTHUVHbb(const double sqrt_s) const
    {
        return 1.0;
    }
    virtual double muTHUttHbb(const double sqrt_s) const
    {
        return 1.0;
    }
    
    virtual double muTHUVBFBRinv(const double sqrt_s) const
    {
        return 0.0;
    }
    
    virtual double muTHUVBFHinv(const double sqrt_s) const
    {
        return 1.0;
    }
    
    virtual double muTHUVHBRinv(const double sqrt_s) const
    {
        return 0.0;
    }
    
    virtual double muTHUVHinv(const double sqrt_s) const
    {
        return 1.0;
    }
       
    virtual double muTHUggHZZ4mu(const double sqrt_s) const
    {
        return 1.0;
    }
    
    virtual double muTHUggHZgamumu(const double sqrt_s) const
    {
        return 1.0;
    }
    

    ////////////////////////////////////////////////////////////////////////
      
    /**
     * @brief The new physics contribution to the anomalous triple gauge coupling @f$g_{1,Z}@f$.
     * @return @f$\delta g_{1,Z}@f$
     */
    virtual double deltag1ZNP() const
    {
        return 0.0;
    }
      
    /**
     * @brief The new physics contribution to the anomalous triple gauge coupling @f$\kappa_{\gamma}@f$.
     * @return @f$\delta \kappa_{\gamma}@f$
     */
    virtual double deltaKgammaNP() const
    {
        return 0.0;
    }
      
    /**
     * @brief The new physics contribution to the anomalous triple gauge coupling @f$\lambda_{Z}@f$.
     * @return @f$\lambda_{Z}@f$
     */
    virtual double lambdaZNP() const
    {
        return 0.0;
    }
    
    
    ////////////////////////////////////////////////////////////////////////
      
    /**
     * @brief The new physics contribution to the effective anomalous triple 
     * gauge coupling @f$g_{1,Z}^{Eff}@f$ from arXiv: 1708.09079 [hep-ph].
     * @return @f$\delta g_{1,Z}@f$
     */
    virtual double deltag1ZNPEff() const
    {
        return 0.0;
    }
      
    /**
     * @brief The new physics contribution to the effective anomalous triple 
     * gauge coupling @f$\kappa_{\gamma}^{Eff}@f$ from arXiv: 1708.09079 [hep-ph].
     * @return @f$\delta \kappa_{\gamma}@f$
     */
    virtual double deltaKgammaNPEff() const
    {
        return 0.0;
    }
    
    ////////////////////////////////////////////////////////////////////////
    
    /**
     * @brief The cross section in pb for @f$e^+ e^- \to W^+ W^- \to 4f @f$, 
     * with @f$ 4f = 0 (jjjj), 1 (e v jj), 2 (mu v jj), 3 (tau v jj),
     * 4 (e v e v), 5 (mu v mu v), 6 (tau v tau v),
     * 7 (e v mu v), 8 (e v tau v), 9 (mu v tau v), 10 (l v jj), 11 (l v l v) @f$
     * the different fermion final states for C.O.M. energies in 188-208 GeV.
     * From arXiv: 1606.06693 [hep-ph]. Defined only for the NPSMEFTd6 class.
     * @return @f$sigma@f$ [pb]
     */
    virtual double xseeWW4fLEP2(const double sqrt_s, const int fstate) const
    {
        return 0.0;
    }
    
    /**
     * @brief The new physics contribution to the cross section in pb for @f$e^+ e^- \to W^+ W^- \to 4f @f$, 
     * with @f$ 4f = 0 (jjjj), 1 (e v jj), 2 (mu v jj), 3 (tau v jj),
     * 4 (e v e v), 5 (mu v mu v), 6 (tau v tau v),
     * 7 (e v mu v), 8 (e v tau v), 9 (mu v tau v), 10 (l v jj), 11 (l v l v) @f$
     * the different fermion final states for C.O.M. energies in 188-208 GeV.
     * From arXiv: 1606.06693 [hep-ph].
     * @return @f$\delta sigma@f$ [pb]
     */
    virtual double deltaxseeWW4fLEP2(const double sqrt_s, const int fstate) const
    {
        return 0.0;
    }

    /**
     * @brief The cross section in pb for @f$e^+ e^- \to W^+ W^- \to \ell \nu \ell \nu@f$, 
     * summing over all final states for C.O.M. energies in 188-208 GeV.
     * From arXiv: 1606.06693 [hep-ph]. Defined only for the NPSMEFTd6 class.
     * @return @f$sigma@f$ [pb]
     */
    virtual double xseeWWleptLEP2(const double sqrt_s) const
    {
        return xseeWW4fLEP2(sqrt_s, 11);
    }
    
/**
     * @brief The new physics contribution to the cross section in pb for @f$e^+ e^- \to W^+ W^- \to \ell \nu \ell \nu@f$, 
     * summing over all final states for C.O.M. energies in 188-208 GeV.
     * From arXiv: 1606.06693 [hep-ph]. Defined only for the NPSMEFTd6 class.
     * @return @f$\delta sigma@f$ [pb]
     */
    virtual double deltaxseeWWleptLEP2(const double sqrt_s) const
    {
        return deltaxseeWW4fLEP2(sqrt_s, 11);
    }
    
    /**
     * @brief The cross section in pb for @f$e^+ e^- \to W^+ W^- \to \ell \nu j j@f$, 
     * summing over all final states for C.O.M. energies in 188-208 GeV.
     * From arXiv: 1606.06693 [hep-ph]. Defined only for the NPSMEFTd6 class.
     * @return @f$sigma@f$ [pb]
     */
    virtual double xseeWWsemilLEP2(const double sqrt_s) const
    {
        return xseeWW4fLEP2(sqrt_s, 10);
    }  
    
    /**
     * @brief The new physics contribution to the cross section in pb for @f$e^+ e^- \to W^+ W^- \to \ell \nu j j@f$, 
     * summing over all final states for C.O.M. energies in 188-208 GeV.
     * From arXiv: 1606.06693 [hep-ph]. Defined only for the NPSMEFTd6 class.
     * @return @f$\delta sigma@f$ [pb]
     */
    virtual double deltaxseeWWsemilLEP2(const double sqrt_s) const
    {
        return deltaxseeWW4fLEP2(sqrt_s, 10);
    } 
    
    /**
     * @brief The cross section in pb for @f$e^+ e^- \to W^+ W^- \to j j j j@f$, 
     * summing over all final states for C.O.M. energies in 188-208 GeV.
     * From arXiv: 1606.06693 [hep-ph]. Defined only for the NPSMEFTd6 class.
     * @return @f$sigma@f$ [pb]
     */
    virtual double xseeWWhadLEP2(const double sqrt_s) const
    {
        return xseeWW4fLEP2(sqrt_s, 0);
    }
    
    /**
     * @brief The new physics contribution to the cross section in pb for @f$e^+ e^- \to W^+ W^- \to j j j j@f$, 
     * summing over all final states for C.O.M. energies in 188-208 GeV.
     * From arXiv: 1606.06693 [hep-ph]. Defined only for the NPSMEFTd6 class.
     * @return @f$\delta sigma@f$ [pb]
     */
    virtual double deltaxseeWWhadLEP2(const double sqrt_s) const
    {
        return deltaxseeWW4fLEP2(sqrt_s, 0);
    }
    
    /**
     * @brief The total cross section in pb for @f$e^+ e^- \to W^+ W^-@f$, 
     * summing over all final states for C.O.M. energies in 188-208 GeV.
     * From arXiv: 1606.06693 [hep-ph]. Defined only for the NPSMEFTd6 class.
     * @return @f$sigma@f$ [pb]
     */
    virtual double xseeWWtotLEP2(const double sqrt_s) const
    {
        return 0.0;
    }
    
    /**
     * @brief The new physics contribution to the total cross section in pb for @f$e^+ e^- \to W^+ W^-@f$, 
     * summing over all final states for C.O.M. energies in 188-208 GeV.
     * From arXiv: 1606.06693 [hep-ph].
     * @return @f$\delta sigma@f$ [pb]
     */
    virtual double deltaxseeWWtotLEP2(const double sqrt_s) const
    {
        return 0.0;
    }
    
    /**
     * @brief The differential cross section in pb for @f$e^+ e^- \to W^+ W^- \to lv jj @f$, 
     * with @f$ l= e,\mu @f$ for the 4 @f$ cos{\theta}@f$ bins defined in arXiv: 1606.06693 [hep-ph].
     * for the C.O.M. energies of 182.6 and 205.9 GeV.
     * From arXiv: 1606.06693 [hep-ph]. Defined only for the NPSMEFTd6 class.
     * @return @f$d\sigma/d\cos{\theta}@f$ [pb]
     */
    virtual double dxsdcoseeWWlvjjLEP2(const double sqrt_s, const int bin) const
    {
        return 0.0;
    }
    
    /**
     * @brief The new physics contribution to the differential cross section in pb for @f$e^+ e^- \to W^+ W^- \to lv jj @f$, 
     * with @f$ l= e,\mu @f$ for the 4 @f$ cos{\theta}@f$ bins defined in arXiv: 1606.06693 [hep-ph].
     * for the C.O.M. energies of 182.6 and 205.9 GeV.
     * From arXiv: 1606.06693 [hep-ph].
     * @return @f$\delta d\sigma/d\cos{\theta}@f$ [pb]
     */
    virtual double deltadxsdcoseeWWlvjjLEP2(const double sqrt_s, const int bin) const
    {
        return 0.0;
    }
    
    /**
     * @brief The differential distribution for @f$e^+ e^- \to W^+ W^- \to jj \ell \nu@f$, 
     * with @f$\ell= e, \mu@f$, as a function of the @f$W@f$ polar angle.
     * @return @f$d\sigma/d\cos{\theta}@f$
     */
    virtual double dxseeWWdcos(const double sqrt_s, const double cos) const
    {
        return 0.0;
    }
    
    /**
     * @brief The integral of differential distribution for @f$e^+ e^- \to W^+ W^- \to jj \ell \nu@f$, 
     * with @f$\ell= e, \mu@f$ in a given bin of the @f$W@f$ polar angle.
     * @return @f$\int_{\cos{\theta_1}}^{\cos{\theta_2}} d\sigma/d\cos{\theta}@f$
     */
    virtual double dxseeWWdcosBin(const double sqrt_s, const double cos1, const double cos2) const
    {
        return 0.0;
    }
    
    /**
     * @brief Total @f$e^+ e^- \to W^+ W^- \to jj \ell \nu@f$ cross section in pb, 
     * with @f$\ell= e, \mu@f$.
     * @return @f$\sigma(e^+ e^- \to W^+ W^- \to jj \ell \nu) @f$
     */
    virtual double xseeWW(const double sqrt_s) const
    {
        return 0.0;
    }
    
    /**
     * @brief The ratio @f$\mu_{eeWW}@f$ between the 
     * @f$ e^{+}e^{-}\to W^{+}W^{-} @f$ production
     * cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{eeWW}@f$
     */
    virtual double mueeWW(const double sqrt_s) const
    {
        return 1.0;
    }
    
    /**
     * @brief The ratio @f$\mu_{eeWW}@f$ between the 
     * @f$ e^{+}e^{-}\to W^{+}W^{-} @f$ production
     * cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV, Pol_em and Pol_ep
     * are the polarization of electrons and positrons, respectively 
     * @return @f$\mu_{eeWW}@f$
     */
    virtual double mueeWWPol(const double sqrt_s, const double Pol_em, const double Pol_ep) const
    {
        return 1.0;
    }
    
    ////////////////////////////////////////////////////////////////////////
    
    //----- High Energy diboson observables at hadron colliders

    /**
     * @brief The direction constrained by @f$ p p \to Z H@f$ in the boosted regime, @f$g_p^Z@f$.
     * From arXiv:1807.01796 and the contribution to FCC CDR Vol 1. Implemented only in NPSMEFTd6 class.
     * @return @f$g_p^Z@f$
     */
    virtual double ppZHprobe(const double sqrt_s) const
    {
        return 0.0;
    }
    
    /**
     * @brief The number of events in  @f$ p p \to WZ@f$
     * in a given @f$p_{TV}@f$ bin, normalized to the SM prediction.
     * From arXiv: 1712.01310 [hep-ph] and private communication.
     * Implemented only in NPSMEFTd6 class.
     * @return @f$N_{ev}^{p_{TV}}/N_{ev,SM}^{p_{TV}}@f$
     */
    virtual double mupTVppWZ(const double sqrt_s, const double pTV1, const double pTV2) const
    {
        return 1.0;
    }
    
    
    
    ////////////////////////////////////////////////////////////////////////
    
    //----- Simplified Template Cross Sections Bins
    
    //----- Stage 0
    
    /**
     * @brief The STXS0 bin @f$pp \to H qq@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual double STXS0_qqH(const double sqrt_s) const    
    {
        return 1.0;
    }
    
    
    //----- Stage 1
    
    /**
     * @brief The STXS bin @f$gg \to H@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual double STXS_ggH_VBFtopo_j3v(const double sqrt_s) const
    {
        return 1.0;
    }
    
    
    /**
     * @brief The STXS bin @f$gg \to H@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual double STXS_ggH_VBFtopo_j3(const double sqrt_s) const
    {
        return 1.0;
    }
    
    
    /**
     * @brief The STXS bin @f$gg \to H@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual double STXS_ggH0j(const double sqrt_s) const
    {
        return 1.0;
    }
    
    /**
     * @brief The STXS bin @f$gg \to H@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual double STXS_ggH1j_pTH_0_60(const double sqrt_s) const
    {
        return 1.0;
    }
    
    /**
     * @brief The STXS bin @f$gg \to H@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual double STXS_ggH1j_pTH_60_120(const double sqrt_s) const
    {
        return 1.0;
    }
    
    /**
     * @brief The STXS bin @f$gg \to H@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual double STXS_ggH1j_pTH_120_200(const double sqrt_s) const
    {
        return 1.0;
    }
    
    /**
     * @brief The STXS bin @f$gg \to H@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual double STXS_ggH1j_pTH_200(const double sqrt_s) const
    {
        return 1.0;
    }
    
    /**
     * @brief The STXS bin @f$gg \to H@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual double STXS_ggH2j_pTH_0_200(const double sqrt_s) const
    {
        return 1.0;
    }
    
    
    /**
     * @brief The STXS bin @f$gg \to H@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual double STXS_ggH2j_pTH_0_60(const double sqrt_s) const
    {
        return 1.0;
    }
    
    
    /**
     * @brief The STXS bin @f$gg \to H@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual double STXS_ggH2j_pTH_60_120(const double sqrt_s) const
    {
        return 1.0;
    }
    
        
    /**
     * @brief The STXS bin @f$gg \to H@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual double STXS_ggH2j_pTH_120_200(const double sqrt_s) const
    {
        return 1.0;
    }
    
    
    /**
     * @brief The STXS bin @f$gg \to H@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual double STXS_ggH2j_pTH_200(const double sqrt_s) const
    {
        return 1.0;
    }
    
  
    /**
     * @brief The STXS bin @f$qq \to H qq@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual double STXS_qqHqq_VBFtopo_Rest(const double sqrt_s) const
    {
        return 1.0;
    }
    
    
    /**
     * @brief The STXS bin @f$qq \to H qq@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual double STXS_qqHqq_VBFtopo_j3v(const double sqrt_s) const
    {
        return 1.0;
    }
    
    
    /**
     * @brief The STXS bin @f$qq \to H qq@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual double STXS_qqHqq_VBFtopo_j3(const double sqrt_s) const
    {
        return 1.0;
    }
    
    /**
     * @brief The STXS bin @f$qq \to H qq@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual double STXS_qqHqq_nonVHtopo(const double sqrt_s) const
    {
        return 1.0;
    }
    
    /**
     * @brief The STXS bin @f$qq \to H qq@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual double STXS_qqHqq_VHtopo(const double sqrt_s) const
    {
        return 1.0;
    }
    
    
    /**
     * @brief The STXS bin @f$qq \to H qq@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual double STXS_qqHqq_Rest(const double sqrt_s) const
    {
        return 1.0;
    }
    
    /**
     * @brief The STXS bin @f$qq \to H qq@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual double STXS_qqHqq_pTj_200(const double sqrt_s) const
    {
        return 1.0;
    }
    
    /**
     * @brief The STXS bin @f$qq \to H \ell \nu@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual double STXS_qqHlv_pTV_0_250(const double sqrt_s) const
    {
        return 1.0;
    }
        
    
    /**
     * @brief The STXS bin @f$qq \to H \ell \nu@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual double STXS_qqHlv_pTV_0_150(const double sqrt_s) const
    {
        return 1.0;
    }
    
    
    /**
     * @brief The STXS bin @f$qq \to H \ell \nu@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual double STXS_qqHlv_pTV_150_250_0j(const double sqrt_s) const
    {
        return 1.0;
    }
    
    
    /**
     * @brief The STXS bin @f$qq \to H \ell \nu@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual double STXS_qqHlv_pTV_150_250_1j(const double sqrt_s) const
    {
        return 1.0;
    }
    
    
    /**
     * @brief The STXS bin @f$qq \to H \ell \nu@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual double STXS_qqHlv_pTV_250(const double sqrt_s) const
    {
        return 1.0;
    }
    
    /**
     * @brief The STXS bin @f$qq \to H \ell \ell@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual double STXS_qqHll_pTV_0_150(const double sqrt_s) const
    {
        return 1.0;
    }
    
    /**
     * @brief The STXS bin @f$qq \to H \ell \ell@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual double STXS_qqHll_pTV_150_250(const double sqrt_s) const
    {
        return 1.0;
    }
    
    
    /**
     * @brief The STXS bin @f$qq \to H \ell \ell@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual double STXS_qqHll_pTV_150_250_0j(const double sqrt_s) const
    {
        return 1.0;
    }
    
    
    /**
     * @brief The STXS bin @f$qq \to H \ell \ell@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual double STXS_qqHll_pTV_150_250_1j(const double sqrt_s) const
    {
        return 1.0;
    }   
    
    
    /**
     * @brief The STXS bin @f$qq \to H \ell \ell@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual double STXS_qqHll_pTV_250(const double sqrt_s) const
    {
        return 1.0;
    }
    
    
    /**
     * @brief The STXS bin @f$ ttH + tH @f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual double STXS_ttHtH(const double sqrt_s) const
    {
        return 1.0;
    }
    
    
    /**
     * @brief The STXS bin @f$ qq \to WH \to H qq @f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual double STXS_WHqqHqq_VBFtopo_j3v(const double sqrt_s) const
    {
        return 1.0;
    }
    
    /**
     * @brief The STXS bin @f$ qq \to WH \to H qq @f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual double STXS_WHqqHqq_VBFtopo_j3(const double sqrt_s) const
    {
        return 1.0;
    }
    
    /**
     * @brief The STXS bin @f$ qq \to WH \to H qq @f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual double STXS_WHqqHqq_VH2j(const double sqrt_s) const
    {
        return 1.0;
    }
    
    /**
     * @brief The STXS bin @f$ qq \to WH \to H qq @f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual double STXS_WHqqHqq_Rest(const double sqrt_s) const
    {
        return 1.0;
    }
    
    /**
     * @brief The STXS bin @f$ qq \to WH \to H qq @f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual double STXS_WHqqHqq_pTj1_200(const double sqrt_s) const
    {
        return 1.0;
    }
    
    /**
     * @brief The STXS bin @f$ qq \to ZH \to H qq @f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual double STXS_ZHqqHqq_VBFtopo_j3v(const double sqrt_s) const
    {
        return 1.0;
    }
    
    /**
     * @brief The STXS bin @f$ qq \to ZH \to H qq @f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual double STXS_ZHqqHqq_VBFtopo_j3(const double sqrt_s) const
    {
        return 1.0;
    }
    
    /**
     * @brief The STXS bin @f$ qq \to ZH \to H qq @f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual double STXS_ZHqqHqq_VH2j(const double sqrt_s) const
    {
        return 1.0;
    }
    
    /**
     * @brief The STXS bin @f$ qq \to ZH \to H qq @f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual double STXS_ZHqqHqq_Rest(const double sqrt_s) const
    {
        return 1.0;
    }
    
    /**
     * @brief The STXS bin @f$ qq \to ZH \to H qq @f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual double STXS_ZHqqHqq_pTj1_200(const double sqrt_s) const
    {
        return 1.0;
    }
    

    //----- Stage 1.2 
    // From ATLAS-CONF-2020-053
    // Expressions valid in the {G_F, M_Z, M_W} scheme
    
    /**
     * @brief The STXS BR @f$ H \to 4l @f$, @f$l=e,\mu@f$.
     */
    virtual double STXS12_BrH4lRatio() const
    {
        return 1.0;
    }
    
    /**
     * @brief The STXS BR @f$ H \to e\nu \mu\nu @f$.
     */
    virtual double STXS12_BrHevmuvRatio() const
    {
        return 1.0;
    }
    
    /**
     * @brief The STXS BR @f$ H \to \gamma \gamma @f$.
     */
    virtual double STXS12_BrHgagaRatio() const
    {
        return 1.0;
    }
    
    /**
     * @brief The STXS BR @f$ H \to bb @f$.
     */
    virtual double STXS12_BrHbbRatio() const
    {
        return 1.0;
    }
    
    /**
     * @brief The STXS bin @f$gg \to H@f$, @f$N_j\leq 1,~200<p_{TH} [GeV]<300@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual double STXS12_ggH_pTH200_300_Nj01(const double sqrt_s) const
    {
        return 1.0;
    }
    
    /**
     * @brief The STXS bin @f$gg \to H@f$, @f$N_j\leq 1,~300<p_{TH} [GeV]<450@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual double STXS12_ggH_pTH300_450_Nj01(const double sqrt_s) const
    {
        return 1.0;
    }
    
    /**
     * @brief The STXS bin @f$gg \to H@f$, @f$N_j\leq 1,~450<p_{TH} [GeV]<650@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual double STXS12_ggH_pTH450_650_Nj01(const double sqrt_s) const
    {
        return 1.0;
    }
    
    /**
     * @brief The STXS bin @f$gg \to H@f$, @f$N_j\leq 1,650<p_{TH} [GeV]@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual double STXS12_ggH_pTH650_Inf_Nj01(const double sqrt_s) const
    {
        return 1.0;
    }
    
    /**
     * @brief The STXS bin @f$gg \to H@f$, @f$N_j=0,~p_{TH} [GeV]<10@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual double STXS12_ggH_pTH0_10_Nj0(const double sqrt_s) const
    {
        return 1.0;
    }
    
    /**
     * @brief The STXS bin @f$gg \to H@f$, @f$N_j=0,~10<p_{TH} [GeV]@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual double STXS12_ggH_pTH10_Inf_Nj0(const double sqrt_s) const
    {
        return 1.0;
    }
    
    /**
     * @brief The STXS bin @f$gg \to H@f$, @f$N_j = 1,~p_{TH} [GeV]<60@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual double STXS12_ggH_pTH0_60_Nj1(const double sqrt_s) const
    {
        return 1.0;
    }
    
    /**
     * @brief The STXS bin @f$gg \to H@f$, @f$N_j = 1,~60<p_{TH} [GeV]<120@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual double STXS12_ggH_pTH60_120_Nj1(const double sqrt_s) const
    {
        return 1.0;
    }
    
    /**
     * @brief The STXS bin @f$gg \to H@f$, @f$N_j = 1,~120<p_{TH} [GeV]<200@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual double STXS12_ggH_pTH120_200_Nj1(const double sqrt_s) const
    {
        return 1.0;
    }
      
    /**
     * @brief The STXS bin @f$gg \to H@f$, @f$N_j \geq 2,~m_{jj}[GeV]<350,~p_{TH} [GeV]<60@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual double STXS12_ggH_mjj0_350_pTH0_60_Nj2(const double sqrt_s) const
    {
        return 1.0;
    }
    
    /**
     * @brief The STXS bin @f$gg \to H@f$, @f$N_j \geq 2,~m_{jj}[GeV]<350,~60<p_{TH} [GeV]<120@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual double STXS12_ggH_mjj0_350_pTH60_120_Nj2(const double sqrt_s) const
    {
        return 1.0;
    }
    
    /**
     * @brief The STXS bin @f$gg \to H@f$, @f$N_j \geq 2,~m_{jj}[GeV]<350,~120<p_{TH} [GeV]<200@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual double STXS12_ggH_mjj0_350_pTH120_200_Nj2(const double sqrt_s) const
    {
        return 1.0;
    }
    
    /**
     * @brief The STXS bin @f$gg \to H@f$, @f$N_j \geq 2,~350<m_{jj}[GeV]<700,~p_{TH} [GeV]<200,~p_{THjj}[GeV]<25@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual double STXS12_ggH_mjj350_700_pTH0_200_ptHjj0_25_Nj2(const double sqrt_s) const
    {
        return 1.0;
    }

    /**
     * @brief The STXS bin @f$gg \to H@f$, @f$N_j \geq 2,~350<m_{jj}[GeV]<700,~p_{TH} [GeV]<200,~25<p_{THjj}[GeV]@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual double STXS12_ggH_mjj350_700_pTH0_200_ptHjj25_Inf_Nj2(const double sqrt_s) const
    {
        return 1.0;
    }    

    /**
     * @brief The STXS bin @f$gg \to H@f$, @f$N_j \geq 2,~700<m_{jj}[GeV],~p_{TH} [GeV]<200,~p_{THjj}[GeV]<25@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual double STXS12_ggH_mjj700_Inf_pTH0_200_ptHjj0_25_Nj2(const double sqrt_s) const
    {
        return 1.0;
    }

    /**
     * @brief The STXS bin @f$gg \to H@f$, @f$N_j \geq 2,~700<m_{jj}[GeV],~p_{TH} [GeV]<200,~25<p_{THjj}[GeV]@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual double STXS12_ggH_mjj700_Inf_pTH0_200_ptHjj25_Inf_Nj2(const double sqrt_s) const
    {
        return 1.0;
    }
    
    /**
     * @brief The STXS bin @f$gg \to H@f$, @f$N_j \geq 2,~350<m_{jj}[GeV]<700,~p_{TH} [GeV]<200@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual double STXS12_ggH_mjj350_700_pTH0_200_Nj2(const double sqrt_s) const        //AG:added
    {
        return 1.0;
    }
    
    /**
     * @brief The STXS bin @f$gg \to H@f$, @f$N_j \geq 2,~700<m_{jj}[GeV],~p_{TH} [GeV]<200@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual double STXS12_ggH_mjj700_Inf_pTH0_200_Nj2(const double sqrt_s) const        //AG:added
    {
        return 1.0;
    }

    
    /**
     * @brief The STXS bin @f$gg \to H\ell\ell@f$, @f$p_{TV}[GeV]<75@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual double STXS12_ggHll_pTV0_75(const double sqrt_s) const
    {
        return 1.0;
    }
    
    /**
     * @brief The STXS bin @f$gg \to H\ell\ell@f$, @f$75<p_{TV}[GeV]<150@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual double STXS12_ggHll_pTV75_150(const double sqrt_s) const
    {
        return 1.0;
    }
    
    /**
     * @brief The STXS bin @f$gg \to H\ell\ell@f$, @f$N_j = 0,~150<p_{TV}[GeV]<250@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual double STXS12_ggHll_pTV150_250_Nj0(const double sqrt_s) const
    {
        return 1.0;
    }
    
    /**
     * @brief The STXS bin @f$gg \to H\ell\ell@f$, @f$N_j = 1,~150<p_{TV}[GeV]<250@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual double STXS12_ggHll_pTV150_250_Nj1(const double sqrt_s) const
    {
        return 1.0;
    }
    
    /**
     * @brief The STXS bin @f$gg \to H\ell\ell@f$, @f$250 < p_{TV}[GeV]@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual double STXS12_ggHll_pTV250_Inf(const double sqrt_s) const
    {
        return 1.0;
    }
    
    
    /**
     * @brief The STXS bin @f$qq \to Hqq@f$, @f$N_j = 0@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual double STXS12_qqHqq_Nj0(const double sqrt_s) const
    {
        return 1.0;
    }
    
    /**
     * @brief The STXS bin @f$qq \to Hqq@f$, @f$N_j = 1@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual double STXS12_qqHqq_Nj1(const double sqrt_s) const
    {
        return 1.0;
    }
    
    /**
     * @brief The STXS bin @f$qq \to Hqq@f$, @f$N_j \geq 2,~m_{jj}[GeV]<60@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual double STXS12_qqHqq_mjj0_60_Nj2(const double sqrt_s) const
    {
        return 1.0;
    }
    
    /**
     * @brief The STXS bin @f$qq \to Hqq@f$, @f$N_j \geq 2,~60<m_{jj}[GeV]<120@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual double STXS12_qqHqq_mjj60_120_Nj2(const double sqrt_s) const
    {
        return 1.0;
    }
    
    /**
     * @brief The STXS bin @f$qq \to Hqq@f$, @f$N_j \geq 2,~120<m_{jj}[GeV]<350@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual double STXS12_qqHqq_mjj120_350_Nj2(const double sqrt_s) const
    {
        return 1.0;
    }
    
    /**
     * @brief The STXS bin @f$qq \to Hqq@f$, @f$N_j \geq 2,~350<m_{jj}[GeV],~200<p_{TH}[GeV]@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual double STXS12_qqHqq_mjj350_Inf_pTH200_Inf_Nj2(const double sqrt_s) const
    {
        return 1.0;
    }
    
    /**
     * @brief The STXS bin @f$qq \to Hqq@f$, @f$N_j \geq 2,~350<m_{jj}[GeV]<700,~p_{TH}[GeV]<200,~p_{THjj}[GeV]<25@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual double STXS12_qqHqq_mjj350_700_pTH0_200_pTHjj0_25_Nj2(const double sqrt_s) const
    {
        return 1.0;
    }
    
    /**
     * @brief The STXS bin @f$qq \to Hqq@f$, @f$N_j \geq 2,~350<m_{jj}[GeV]<700,~p_{TH}[GeV]<200,~25<p_{THjj}[GeV]@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual double STXS12_qqHqq_mjj350_700_pTH0_200_pTHjj25_Inf_Nj2(const double sqrt_s) const
    {
        return 1.0;
    }    
    
    /**
     * @brief The STXS bin @f$qq \to Hqq@f$, @f$N_j \geq 2,~700<m_{jj}[GeV],~p_{TH}[GeV]<200,~p_{THjj}[GeV]<25@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual double STXS12_qqHqq_mjj700_Inf_pTH0_200_pTHjj0_25_Nj2(const double sqrt_s) const
    {
        return 1.0;
    }
    
    /**
     * @brief The STXS bin @f$qq \to Hqq@f$, @f$N_j \geq 2,~700<m_{jj}[GeV],~p_{TH}[GeV]<200,~25<p_{THjj}[GeV]@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual double STXS12_qqHqq_mjj700_Inf_pTH0_200_pTHjj25_Inf_Nj2(const double sqrt_s) const
    {
        return 1.0;
    }
    
    /**
     * @brief The STXS bin @f$qq \to Hqq@f$, @f$N_j \geq 2,~350<m_{jj}[GeV]<700,~p_{TH}[GeV]<200@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual double STXS12_qqHqq_mjj350_700_pTH0_200_Nj2(const double sqrt_s) const     //AG:added
    {
        return 1.0;
    }
    
    /**
     * @brief The STXS bin @f$qq \to Hqq@f$, @f$N_j \geq 2,~700<m_{jj}[GeV]<1000,~p_{TH}[GeV]<200@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual double STXS12_qqHqq_mjj700_1000_pTH0_200_Nj2(const double sqrt_s) const     //AG:added
    {
        return 1.0;
    }
    
    /**
     * @brief The STXS bin @f$qq \to Hqq@f$, @f$N_j \geq 2,~1000<m_{jj}[GeV]<1500,~p_{TH}[GeV]<200@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual double STXS12_qqHqq_mjj1000_1500_pTH0_200_Nj2(const double sqrt_s) const     //AG:added
    {
        return 1.0;
    }
    
    /**
     * @brief The STXS bin @f$qq \to Hqq@f$, @f$N_j \geq 2,~1500<m_{jj}[GeV],~p_{TH}[GeV]<200@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual double STXS12_qqHqq_mjj1500_Inf_pTH0_200_Nj2(const double sqrt_s) const     //AG:added
    {
        return 1.0;
    }
    
    /**
     * @brief The STXS bin @f$qq \to Hqq@f$, @f$N_j \geq 2,~350<m_{jj}[GeV]<1000,~p_{TH}[GeV]>200@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual double STXS12_qqHqq_mjj350_1000_pTH200_Inf_Nj2(const double sqrt_s) const     //AG:added
    {
        return 1.0;
    }
    
    /**
     * @brief The STXS bin @f$qq \to Hqq@f$, @f$N_j \geq 2,~1000<m_{jj}[GeV],~p_{TH}[GeV]>200@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual double STXS12_qqHqq_mjj1000_Inf_pTH200_Inf_Nj2(const double sqrt_s) const     //AG:added
    {
        return 1.0;
    }
    
    /**
     * @brief The STXS bin @f$qq \to H\ell\nu@f$, @f$p_{TV}[GeV]<75@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual double STXS12_qqHlv_pTV0_75(const double sqrt_s) const
    {
        return 1.0;
    }
    
    /**
     * @brief The STXS bin @f$qq \to H\ell\nu@f$, @f$75<p_{TV}[GeV]<150@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual double STXS12_qqHlv_pTV75_150(const double sqrt_s) const
    {
        return 1.0;
    }
    
    /**
     * @brief The STXS bin @f$qq \to H\ell\nu@f$, @f$N_j = 0,~150<p_{TV}[GeV]<250@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual double STXS12_qqHlv_pTV150_250_Nj0(const double sqrt_s) const
    {
        return 1.0;
    }
    
    /**
     * @brief The STXS bin @f$qq \to H\ell\nu@f$, @f$N_j \geq 1,~150<p_{TV}[GeV]<250@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual double STXS12_qqHlv_pTV150_250_Nj1(const double sqrt_s) const
    {
        return 1.0;
    }
    
    /**
     * @brief The STXS bin @f$qq \to H\ell\nu@f$, @f$250<p_{TV}[GeV]@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual double STXS12_qqHlv_pTV250_Inf(const double sqrt_s) const
    {
        return 1.0;
    }
    
    /**
     * @brief The STXS bin @f$qq \to H\ell\nu@f$, @f$0<p_{TV}<150[GeV]@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual double STXS12_qqHlv_pTV0_150(const double sqrt_s) const             //AG:added
    {
        return 1.0;
    }
    
    /**
     * @brief The STXS bin @f$qq \to H\ell\nu@f$, @f$250<p_{TV}<400[GeV]@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual double STXS12_qqHlv_pTV250_400(const double sqrt_s) const           //AG:added
    {
        return 1.0;
    }
    
    /**
     * @brief The STXS bin @f$qq \to H\ell\nu@f$, @f$400<p_{TV}[GeV]@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual double STXS12_qqHlv_pTV400_Inf(const double sqrt_s) const           //AG:added
    {
        return 1.0;
    }
    
    /**
     * @brief The STXS bin @f$qq \to H\ell\ell@f$, @f$p_{TV}[GeV]<75@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual double STXS12_qqHll_pTV0_75(const double sqrt_s) const
    {
        return 1.0;
    }
        
    /**
     * @brief The STXS bin @f$qq \to H\ell\ell@f$, @f$75<p_{TV}[GeV]<150@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual double STXS12_qqHll_pTV75_150(const double sqrt_s) const
    {
        return 1.0;
    }
        
    /**
     * @brief The STXS bin @f$qq \to H\ell\ell@f$, @f$N_j = 0,~150<p_{TV}[GeV]<250@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual double STXS12_qqHll_pTV150_250_Nj0(const double sqrt_s) const
    {
        return 1.0;
    }
    
    /**
     * @brief The STXS bin @f$qq \to H\ell\ell@f$, @f$N_j \geq 1,~150<p_{TV}[GeV]<250@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual double STXS12_qqHll_pTV150_250_Nj1(const double sqrt_s) const
    {
        return 1.0;
    }
    
    /**
     * @brief The STXS bin @f$qq \to H\ell\ell@f$, @f$250<p_{TV}[GeV]@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual double STXS12_qqHll_pTV250_Inf(const double sqrt_s) const
    {
        return 1.0;
    }
     
    /**
     * @brief The STXS bin @f$qq \to H\ell\ell@f$, @f$0<p_{TV}<150[GeV]@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual double STXS12_qqHll_pTV0_150(const double sqrt_s) const             //AG:added
    {
        return 1.0;
    }
    
    /**
     * @brief The STXS bin @f$qq \to H\ell\ell@f$, @f$250<p_{TV}<400[GeV]@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual double STXS12_qqHll_pTV250_400(const double sqrt_s) const           //AG:added
    {
        return 1.0;
    }
    
    /**
     * @brief The STXS bin @f$qq \to H\ell\ell@f$, @f$400<p_{TV}[GeV]@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual double STXS12_qqHll_pTV400_Inf(const double sqrt_s) const           //AG:added
    {
        return 1.0;
    }
    
    /**
     * @brief The STXS bin @f$pp \to ttH@f$, @f$p_{TH}[GeV]<60@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual double STXS12_ttH_pTH0_60(const double sqrt_s) const
    {
        return 1.0;
    }
    
    /**
     * @brief The STXS bin @f$pp \to ttH@f$, @f$60<p_{TH}[GeV]<120@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual double STXS12_ttH_pTH60_120(const double sqrt_s) const
    {
        return 1.0;
    }
    
    /**
     * @brief The STXS bin @f$pp \to ttH@f$, @f$120<p_{TH}[GeV]<200@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual double STXS12_ttH_pTH120_200(const double sqrt_s) const
    {
        return 1.0;
    }
    
    /**
     * @brief The STXS bin @f$pp \to ttH@f$, @f$200<p_{TH}[GeV]<300@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual double STXS12_ttH_pTH200_300(const double sqrt_s) const
    {
        return 1.0;
    }
    
    /**
     * @brief The STXS bin @f$pp \to ttH@f$, @f$300<p_{TH}[GeV]@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual double STXS12_ttH_pTH300_Inf(const double sqrt_s) const
    {
        return 1.0;
    }
    
    /**
     * @brief The STXS bin @f$pp \to ttH@f$, @f$300<p_{TH}[GeV]<450@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual double STXS12_ttH_pTH300_450(const double sqrt_s) const             //AG:added
    {
        return 1.0;
    }
    
     /**
     * @brief The STXS bin @f$pp \to ttH@f$, @f$450<p_{TH}[GeV]@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual double STXS12_ttH_pTH450_Inf(const double sqrt_s) const             //AG:added
    {
        return 1.0;
    }
    
    /**
     * @brief The STXS bin @f$pp \to tH@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual double STXS12_tH(const double sqrt_s) const
    {
        return 1.0;
    }

 
    ////////////////////////////////////////////////////////////////////////
    
    /**
     * @brief The effective coupling @f$\kappa_{\mu,eff}=\sqrt{\Gamma_{H\mu\mu}/\Gamma_{H\mu\mu}^{SM}}@f$.
     * @return @f$\kappa_{\mu,eff}@f$
     */
    virtual double kappamueff() const
    {
        return 1.0;
    }
    
    /**
     * @brief The effective coupling @f$\kappa_{\tau,eff}=\sqrt{\Gamma_{H\tau\tau}/\Gamma_{H\tau\tau}^{SM}}@f$.
     * @return @f$\kappa_{\tau,eff}@f$
     */
    virtual double kappataueff() const
    {
        return 1.0;
    }
    
    /**
     * @brief The effective coupling @f$\kappa_{c,eff}=\sqrt{\Gamma_{Hcc}/\Gamma_{Hcc}^{SM}}@f$.
     * @return @f$\kappa_{c,eff}@f$
     */
    virtual double kappaceff() const
    {
        return 1.0;
    }
    
    /**
     * @brief The effective coupling @f$\kappa_{b,eff}=\sqrt{\Gamma_{Hbb}/\Gamma_{Hbb}^{SM}}@f$.
     * @return @f$\kappa_{b,eff}@f$
     */
    virtual double kappabeff() const
    {
        return 1.0;
    }
    
    /**
     * @brief The effective coupling @f$\kappa_{G,eff}=\sqrt{\Gamma_{HGG}/\Gamma_{HGG}^{SM}}@f$.
     * @return @f$\kappa_{G,eff}@f$
     */
    virtual double kappaGeff() const
    {
        return 1.0;
    }
    
    /**
     * @brief The effective coupling @f$\kappa_{Z,eff}=\sqrt{\Gamma_{HZZ}/\Gamma_{HZZ}^{SM}}@f$.
     * @return @f$\kappa_{Z,eff}@f$
     */
    virtual double kappaZeff() const
    {
        return 1.0;
    }
    
    /**
     * @brief The effective coupling @f$\kappa_{W,eff}=\sqrt{\Gamma_{HWW}/\Gamma_{HWW}^{SM}}@f$.
     * @return @f$\kappa_{W,eff}@f$
     */
    virtual double kappaWeff() const
    {
        return 1.0;
    }
    
    /**
     * @brief The effective coupling @f$\kappa_{A,eff}=\sqrt{\Gamma_{HAA}/\Gamma_{HAA}^{SM}}@f$.
     * @return @f$\kappa_{A,eff}@f$
     */
    virtual double kappaAeff() const
    {
        return 1.0;
    }
    
    /**
     * @brief The effective coupling @f$\kappa_{ZA,eff}=\sqrt{\Gamma_{HZA}/\Gamma_{HZA}^{SM}}@f$.
     * @return @f$\kappa_{ZA,eff}@f$
     */
    virtual double kappaZAeff() const
    {
        return 1.0;
    }
    
    /////////////Basic interactions of the so-called Higgs basis////////////////
    
    /**
     * @brief The Higgs-basis coupling @f$\delta y_t@f$.
     * (See LHCHXSWG-INT-2015-001 document.)
     * @return @f$\delta y_t@f$
     */
    virtual double deltayt_HB() const
    {
        return 0.0;
    }
    
    /**
     * @brief The Higgs-basis coupling @f$\delta y_b@f$.
     * (See LHCHXSWG-INT-2015-001 document.)
     * @return @f$\delta y_b@f$
     */
    virtual double deltayb_HB() const
    {
        return 0.0;
    }
    
    /**
     * @brief The Higgs-basis coupling @f$\delta y_\tau@f$.
     * (See LHCHXSWG-INT-2015-001 document.)
     * @return @f$\delta y_\tau@f$
     */
    virtual double deltaytau_HB() const
    {
        return 0.0;
    }
    
    /**
     * @brief The Higgs-basis coupling @f$\delta y_c@f$.
     * (See LHCHXSWG-INT-2015-001 document.)
     * @return @f$\delta y_c@f$
     */
    virtual double deltayc_HB() const
    {
        return 0.0;
    }
    
    
    /**
     * @brief The Higgs-basis coupling @f$\delta y_\mu@f$.
     * (See LHCHXSWG-INT-2015-001 document.)
     * @return @f$\delta y_\mu@f$
     */
    virtual double deltaymu_HB() const
    {
        return 0.0;
    }
    
    /**
     * @brief The Higgs-basis coupling @f$\delta c_z@f$.
     * (See LHCHXSWG-INT-2015-001 document.)
     * @return @f$\delta c_z@f$
     */
    virtual double deltacZ_HB() const
    {
        return 0.0;
    }
    
    /**
     * @brief The Higgs-basis coupling @f$c_{z\Box}@f$.
     * (See LHCHXSWG-INT-2015-001 document.)
     * @return @f$c_{z\Box}@f$
     */
    virtual double cZBox_HB() const
    {
        return 0.0;
    }
    
    /**
     * @brief The Higgs-basis coupling @f$c_{zz}@f$.
     * (See LHCHXSWG-INT-2015-001 document.)
     * @return @f$c_{zz}@f$
     */
    virtual double cZZ_HB() const
    {
        return 0.0;
    }
    
    /**
     * @brief The Higgs-basis coupling @f$c_{z\gamma}@f$.
     * (See LHCHXSWG-INT-2015-001 document.)
     * @return @f$c_{z\gamma}@f$
     */
    virtual double cZga_HB() const
    {
        return 0.0;
    }
    
    /**
     * @brief The Higgs-basis coupling @f$c_{\gamma\gamma}@f$.
     * (See LHCHXSWG-INT-2015-001 document.)
     * @return @f$c_{\gamma\gamma}@f$
     */
    virtual double cgaga_HB() const
    {
        return 0.0;
    }
    
    /**
     * @brief The Higgs-basis coupling @f$c_{gg}@f$.
     * (See LHCHXSWG-INT-2015-001 document.)
     * @return @f$c_{gg}@f$
     */
    virtual double cgg_HB() const
    {
        return 0.0;
    }
    
    /**
     * @brief The effective Higgs-basis coupling @f$c_{gg}^{Eff}@f$. (Similar to cgg_HB but including modifications of SM loops.)
     * (See arXiv: 1505.00046 [hep-ph] document.)
     * @return @f$c_{gg}^{Eff}@f$
     */
    virtual double cggEff_HB() const
    {
        return 0.0;
    }
    
    /**
     * @brief The Higgs-basis coupling @f$\lambda_{z}@f$.
     * (See LHCHXSWG-INT-2015-001 document.)
     * @return @f$\lambda_{z}@f$
     */
    virtual double lambz_HB() const
    {
        return 0.0;
    }
    
    
    /////////////Combinations of Warsaw basis coefficients constrained by EWPO////////////////
    
    /**
     * @brief Combination of coefficients of the Warsaw basis constrained by EWPO
     * @f$(\hat{C}_{HL}^{(1)})_{11}@f$.
     * @return @f$(\hat{C}_{HL}^{(1)})_{11}@f$
     */
    virtual double CEWHL111() const
    {
        return 0.0;
    }
    
    /**
     * @brief Combination of coefficients of the Warsaw basis constrained by EWPO
     * @f$(\hat{C}_{HL}^{(1)})_{22}@f$.
     * @return @f$(\hat{C}_{HL}^{(1)})_{22}@f$
     */
    virtual double CEWHL122() const
    {
        return 0.0;
    }
    
    /**
     * @brief Combination of coefficients of the Warsaw basis constrained by EWPO
     * @f$(\hat{C}_{HL}^{(1)})_{33}@f$.
     * @return @f$(\hat{C}_{HL}^{(1)})_{33}@f$
     */
    virtual double CEWHL133() const
    {
        return 0.0;
    }
    
    /**
     * @brief Combination of coefficients of the Warsaw basis constrained by EWPO
     * @f$(\hat{C}_{HL}^{(3)})_{11}@f$.
     * @return @f$(\hat{C}_{HL}^{(3)})_{11}@f$
     */
    virtual double CEWHL311() const
    {
        return 0.0;
    }
    
    /**
     * @brief Combination of coefficients of the Warsaw basis constrained by EWPO
     * @f$(\hat{C}_{HL}^{(3)})_{22}@f$.
     * @return @f$(\hat{C}_{HL}^{(3)})_{22}@f$
     */
    virtual double CEWHL322() const
    {
        return 0.0;
    }
    
    /**
     * @brief Combination of coefficients of the Warsaw basis constrained by EWPO
     * @f$(\hat{C}_{HL}^{(3)})_{33}@f$.
     * @return @f$(\hat{C}_{HL}^{(3)})_{33}@f$
     */
    virtual double CEWHL333() const
    {
        return 0.0;
    }
    
    /**
     * @brief Combination of coefficients of the Warsaw basis constrained by EWPO
     * @f$(\hat{C}_{HQ}^{(1)})_{11}@f$.
     * @return @f$(\hat{C}_{HQ}^{(1)})_{11}@f$
     */
    virtual double CEWHQ111() const
    {
        return 0.0;
    }
    
    /**
     * @brief Combination of coefficients of the Warsaw basis constrained by EWPO
     * @f$(\hat{C}_{HQ}^{(1)})_{22}@f$.
     * @return @f$(\hat{C}_{HQ}^{(1)})_{22}@f$
     */
    virtual double CEWHQ122() const
    {
        return 0.0;
    }
    
    /**
     * @brief Combination of coefficients of the Warsaw basis constrained by EWPO
     * @f$(\hat{C}_{HQ}^{(1)})_{33}@f$.
     * @return @f$(\hat{C}_{HQ}^{(1)})_{33}@f$
     */
    virtual double CEWHQ133() const
    {
        return 0.0;
    }
    
    /**
     * @brief Combination of coefficients of the Warsaw basis constrained by EWPO
     * @f$(\hat{C}_{HQ}^{(3)})_{11}@f$.
     * @return @f$(\hat{C}_{HQ}^{(3)})_{11}@f$
     */
    virtual double CEWHQ311() const
    {
        return 0.0;
    }
    
    /**
     * @brief Combination of coefficients of the Warsaw basis constrained by EWPO
     * @f$(\hat{C}_{HQ}^{(3)})_{22}@f$.
     * @return @f$(\hat{C}_{HQ}^{(3)})_{22}@f$
     */
    virtual double CEWHQ322() const
    {
        return 0.0;
    }
    
    /**
     * @brief Combination of coefficients of the Warsaw basis constrained by EWPO
     * @f$(\hat{C}_{HQ}^{(3)})_{33}@f$.
     * @return @f$(\hat{C}_{HQ}^{(3)})_{33}@f$
     */
    virtual double CEWHQ333() const
    {
        return 0.0;
    }
    
    
    /**
     * @brief Combination of coefficients of the Warsaw basis constrained by EWPO
     * @f$(\hat{C}_{HQ}^{(d)})_{33}@f$.
     * @return @f$(\hat{C}_{HQ}^{(d)})_{33}@f$
     */
    virtual double CEWHQd33() const
    {
        return 0.0;
    }
    
        
    /**
     * @brief Combination of coefficients of the Warsaw basis constrained by EWPO
     * @f$(\hat{C}_{He})_{11}@f$.
     * @return @f$(\hat{C}_{He})_{11}@f$
     */
    virtual double CEWHe11() const
    {
        return 0.0;
    }
    
    /**
     * @brief Combination of coefficients of the Warsaw basis constrained by EWPO
     * @f$(\hat{C}_{He})_{22}@f$.
     * @return @f$(\hat{C}_{He})_{22}@f$
     */
    virtual double CEWHe22() const
    {
        return 0.0;
    }
    
    /**
     * @brief Combination of coefficients of the Warsaw basis constrained by EWPO
     * @f$(\hat{C}_{He})_{33}@f$.
     * @return @f$(\hat{C}_{He})_{33}@f$
     */
    virtual double CEWHe33() const
    {
        return 0.0;
    }
    
    
    /**
     * @brief Combination of coefficients of the Warsaw basis constrained by EWPO
     * @f$(\hat{C}_{Hu})_{11}@f$.
     * @return @f$(\hat{C}_{Hu})_{11}@f$
     */
    virtual double CEWHu11() const
    {
        return 0.0;
    }
    
    /**
     * @brief Combination of coefficients of the Warsaw basis constrained by EWPO
     * @f$(\hat{C}_{Hu})_{22}@f$.
     * @return @f$(\hat{C}_{Hu})_{22}@f$
     */
    virtual double CEWHu22() const
    {
        return 0.0;
    }
    
    /**
     * @brief Combination of coefficients of the Warsaw basis constrained by EWPO
     * @f$(\hat{C}_{Hu})_{33}@f$.
     * @return @f$(\hat{C}_{Hu})_{33}@f$
     */
    virtual double CEWHu33() const
    {
        return 0.0;
    }
    

    /**
     * @brief Combination of coefficients of the Warsaw basis constrained by EWPO
     * @f$(\hat{C}_{Hd})_{11}@f$.
     * @return @f$(\hat{C}_{Hd})_{11}@f$
     */
    virtual double CEWHd11() const
    {
        return 0.0;
    }
    
    /**
     * @brief Combination of coefficients of the Warsaw basis constrained by EWPO
     * @f$(\hat{C}_{Hd})_{22}@f$.
     * @return @f$(\hat{C}_{Hd})_{22}@f$
     */
    virtual double CEWHd22() const
    {
        return 0.0;
    }
    
    /**
     * @brief Combination of coefficients of the Warsaw basis constrained by EWPO
     * @f$(\hat{C}_{Hd})_{33}@f$.
     * @return @f$(\hat{C}_{Hd})_{33}@f$
     */
    virtual double CEWHd33() const
    {
        return 0.0;
    }  
    
    
    
    ////////////////////////////////////////////////////////////////////////
    //LEP2 Observables
    
    
//  Absolute corrections to the differential cross section        
    virtual double delta_Dsigma_f(const Particle f, const double s, const double cos) const;
    
//  Absolute corrections to the differential cross section integrated in [cos \theta_{min},cos \theta_{max}] 
//  Valid for f=/=e
    virtual double delta_sigma_f(const Particle f, const double s, const double cosmin, const double cosmax) const;
    virtual double delta_sigma_had(const double s, const double cosmin, const double cosmax) const;
    
//  Total cross sections  (full acceptance)
    virtual double delta_sigmaTot_f(const Particle f, const double s) const;
    
//  Forward-Backward asymmetry (full acceptance). Valid for f!=e !
    virtual double delta_AFB_f(const Particle f, const double s) const;
      
//   Extension of SM observable definitions
    virtual double LEP2sigmaMu(const double s) const;
    virtual double LEP2sigmaTau(const double s) const;
    virtual double LEP2sigmaHadron(const double s) const;
    virtual double LEP2sigmaCharm(const double s) const;
    virtual double LEP2sigmaBottom(const double s) const;
    
    virtual double LEP2AFBmu(const double s) const;
    virtual double LEP2AFBtau(const double s) const;
    virtual double LEP2AFBcharm(const double s) const;
    virtual double LEP2AFBbottom(const double s) const;
    virtual double LEP2Rcharm(const double s) const;
    virtual double LEP2Rbottom(const double s) const;
    
    //LEP2 Differential Observables
    virtual double LEP2dsigmadcosE(const double s, const double cos) const;
    virtual double LEP2dsigmadcosMu(const double s, const double cos) const;
    virtual double LEP2dsigmadcosTau(const double s, const double cos) const;
    
    ////////////////////////////////////////////////////////////////////////     
    // EW low-energy observables: Muon g-2
    
    /**
     * @brief The computation of the anomalous magnetic moment of the muon @f$a_\mu=(g_\mu-2)/2@f$
     * @return @f$a_\mu=(g_\mu-2)/2@f$
     */
      virtual double delta_amuon() const;
      
      
    // EW low-energy observables: Parity violation

    /**
     * @brief The computation of the electron's weak charge
     * @param[in] q2 the @f$Q^2@f$ at which the weak charge is measured
     * @param[in] y
     * @return @f$Q_{w}(e)@f$
     */
      virtual double delta_Qwemoller(const double q2, const double y) const;


    /**
     * @brief The computation of the parity violating asymmetry in Moller scattering
     * @param[in] q2 the @f$Q^2@f$ of the process
     * @param[in] y
     * @return @f$A_{LR}@f$
     */
      virtual double delta_alrmoller(const double q2, const double y) const;


    /**
     * @brief The computation of the proton weak charge: Qwp
     * @return @f$Q_{W}(p)@f$
     */
      virtual double delta_Qwp() const;

      
    /**
     * @brief The computation of the neutron weak charge: Qwn
     * @return @f$Q_{W}(n)@f$
     */
      virtual double delta_Qwn() const;
      
    // EW low-energy observables: neutrino scattering
      
    /**
     * @brief The computation of the correction to the effective neutrino nucleon LH coupling: delta_gLnuN2
     * @return @f$\Delta g_L^2(\nu N)@f$
     */
      virtual double delta_gLnuN2() const;
      
    /**
     * @brief The computation of the correction to the effective neutrino nucleon RH coupling: delta_gRnuN2
     * @return @f$\Delta g_R^2(\nu N)@f$
     */
      virtual double delta_gRnuN2() const;
      
    /**
     * @brief The computation of the correction to the effective (muon) neutrino-electron vector coupling: delta_gVnue
     * @details
     * @return @f$\Delta g_V^{\nu_\mu e}@f$
     */
      virtual double delta_gVnue() const;
      
    /**
     * @brief The computation of the correction to the effective (muon) neutrino-electron vector coupling: delta_gAnue
     * @details
     * @return @f$\Delta g_A^{\nu_\mu e}@f$
     */
      virtual double delta_gAnue() const;
      
      
//   Extension of SM observable definitions
      virtual double amuon() const;
      virtual double Qwemoller(const double q2, const double y) const;
      virtual double alrmoller(const double q2, const double y) const;
      virtual double Qwp() const;
      virtual double Qwn() const;
      virtual double gLnuN2() const;
      virtual double gRnuN2() const;
      virtual double gVnue() const;
      virtual double gAnue() const;
      
    ////////////////////////////////////////////////////////////////////////     
    // Lepton decays
      
    // Lepton Flavor universality tests in Tau decays

    /**
     * @brief The computation of the correction to the LFU ratio @f$g_\mu/ g_e @f$
     * @details
     * @return @f$\delta \Gamma(\tau \to \mu \nu \nu ) / \Gamma(\tau \to e \nu \nu ) @f$
     */
      virtual double delta_TauLFU_gmuge() const;


    /**
     * @brief The computation of the correction to the LFU ratio @f$g_\tau/ g_\mu @f$
     * @details
     * @return @f$\delta \Gamma(\tau \to e \nu \nu ) / \Gamma(\mu \to e \nu \nu ) @f$
     */
      virtual double delta_TauLFU_gtaugmu() const;


    /**
     * @brief The computation of the correction to the LFU ratio @f$g_\tau/ g_\mu @f$
     * @details
     * @return @f$\delta \Gamma(\tau \to \mu \nu \nu ) / \Gamma(\mu \to e \nu \nu ) @f$
     */
      virtual double delta_TauLFU_gtauge() const;
      
//   Extension of SM observable definitions
      virtual double TauLFU_gmuge() const;
      virtual double TauLFU_gtaugmu() const;
      virtual double TauLFU_gtauge() const;

    
    ///////////Collider observables: LHC dilepton events////////////////////////
    
    /**
     * @brief Number of di-electron events at the LHC at 13 TeV
     * @return NevLHCppee13
     */
    virtual double NevLHCppee13(const int i_bin) const
    {
        return 0.0;
    }
    
    /**
     * @brief Number of di-muon events at the LHC at 13 TeV
     * @return NevLHCppmumu13
     */
    virtual double NevLHCppmumu13(const int i_bin) const
    {
        return 0.0;
    }
    
    /**
     * @brief Number of di-tau events at the LHC at 13 TeV
     * @return NevLHCpptautau13
     */
    virtual double NevLHCpptautau13(const int i_bin) const
    {
        return 0.0;
    }
    
    ///////////Collider observables: LHC mono-lepton events////////////////////////
    
    /**
     * @brief Number of mono-electron events at the LHC at 13 TeV
     * @return NevLHCppenu13
     */
    virtual double NevLHCppenu13(const int i_bin) const
    {
        return 0.0;
    }
    
    /**
     * @brief Number of mono-muon events at the LHC at 13 TeV
     * @return NevLHCppmunu13
     */
    virtual double NevLHCppmunu13(const int i_bin) const
    {
        return 0.0;
    }
    
    /**
     * @brief Number of mono-tau events at the LHC at 13 TeV
     * @return NevLHCpptaunu13
     */
    virtual double NevLHCpptaunu13(const int i_bin) const
    {
        return 0.0;
    }
    
    /////////////Auxiliary observables////////////////
    
    /**
     * @brief Auxiliary observable AuxObs_NP1
     * @return AuxObs_NP1
     */
    virtual double AuxObs_NP1() const
    {
        return 0.0;
    }
    
    /**
     * @brief Auxiliary observable AuxObs_NP2
     * @return AuxObs_NP2
     */
    virtual double AuxObs_NP2() const
    {
        return 0.0;
    }
    
    /**
     * @brief Auxiliary observable AuxObs_NP3
     * @return AuxObs_NP3
     */
    virtual double AuxObs_NP3() const
    {
        return 0.0;
    }
    
    /**
     * @brief Auxiliary observable AuxObs_NP4
     * @return AuxObs_NP4
     */
    virtual double AuxObs_NP4() const
    {
        return 0.0;
    }
    
    /**
     * @brief Auxiliary observable AuxObs_NP5
     * @return AuxObs_NP5
     */
    virtual double AuxObs_NP5() const
    {
        return 0.0;
    }
    
    /**
     * @brief Auxiliary observable AuxObs_NP6
     * @return AuxObs_NP6
     */
    virtual double AuxObs_NP6() const
    {
        return 0.0;
    }
    
    /**
     * @brief Auxiliary observable AuxObs_NP7
     * @return AuxObs_NP7
     */
    virtual double AuxObs_NP7() const
    {
        return 0.0;
    }
    
    /**
     * @brief Auxiliary observable AuxObs_NP8
     * @return AuxObs_NP8
     */
    virtual double AuxObs_NP8() const
    {
        return 0.0;
    }
    
    /**
     * @brief Auxiliary observable AuxObs_NP9
     * @return AuxObs_NP9
     */
    virtual double AuxObs_NP9() const
    {
        return 0.0;
    }
    
    /**
     * @brief Auxiliary observable AuxObs_NP10
     * @return AuxObs_NP10
     */
    virtual double AuxObs_NP10() const
    {
        return 0.0;
    }
    
    /**
     * @brief Auxiliary observable AuxObs_NP11
     * @return AuxObs_NP11
     */
    virtual double AuxObs_NP11() const
    {
        return 0.0;
    }
    
    /**
     * @brief Auxiliary observable AuxObs_NP12
     * @return AuxObs_NP12
     */
    virtual double AuxObs_NP12() const
    {
        return 0.0;
    }
    
    /**
     * @brief Auxiliary observable AuxObs_NP13
     * @return AuxObs_NP13
     */
    virtual double AuxObs_NP13() const
    {
        return 0.0;
    }
    
    /**
     * @brief Auxiliary observable AuxObs_NP14
     * @return AuxObs_NP14
     */
    virtual double AuxObs_NP14() const
    {
        return 0.0;
    }
    
    /**
     * @brief Auxiliary observable AuxObs_NP15
     * @return AuxObs_NP15
     */
    virtual double AuxObs_NP15() const
    {
        return 0.0;
    }
    
    /**
     * @brief Auxiliary observable AuxObs_NP16
     * @return AuxObs_NP16
     */
    virtual double AuxObs_NP16() const
    {
        return 0.0;
    }
    
    /**
     * @brief Auxiliary observable AuxObs_NP17
     * @return AuxObs_NP17
     */
    virtual double AuxObs_NP17() const
    {
        return 0.0;
    }
    
    /**
     * @brief Auxiliary observable AuxObs_NP18
     * @return AuxObs_NP18
     */
    virtual double AuxObs_NP18() const
    {
        return 0.0;
    }
    
    /**
     * @brief Auxiliary observable AuxObs_NP19
     * @return AuxObs_NP19
     */
    virtual double AuxObs_NP19() const
    {
        return 0.0;
    }
    
    /**
     * @brief Auxiliary observable AuxObs_NP20
     * @return AuxObs_NP20
     */
    virtual double AuxObs_NP20() const
    {
        return 0.0;
    }
    
    /**
     * @brief Auxiliary observable AuxObs_NP21
     * @return AuxObs_NP21
     */
    virtual double AuxObs_NP21() const
    {
        return 0.0;
    }
    
    /**
     * @brief Auxiliary observable AuxObs_NP22
     * @return AuxObs_NP22
     */
    virtual double AuxObs_NP22() const
    {
        return 0.0;
    }
    
    /**
     * @brief Auxiliary observable AuxObs_NP23
     * @return AuxObs_NP23
     */
    virtual double AuxObs_NP23() const
    {
        return 0.0;
    }
    
    /**
     * @brief Auxiliary observable AuxObs_NP24
     * @return AuxObs_NP24
     */
    virtual double AuxObs_NP24() const
    {
        return 0.0;
    }
    
    /**
     * @brief Auxiliary observable AuxObs_NP25
     * @return AuxObs_NP25
     */
    virtual double AuxObs_NP25() const
    {
        return 0.0;
    }
    
    /**
     * @brief Auxiliary observable AuxObs_NP26
     * @return AuxObs_NP26
     */
    virtual double AuxObs_NP26() const
    {
        return 0.0;
    }
    
    /**
     * @brief Auxiliary observable AuxObs_NP27
     * @return AuxObs_NP27
     */
    virtual double AuxObs_NP27() const
    {
        return 0.0;
    }
    
    /**
     * @brief Auxiliary observable AuxObs_NP28
     * @return AuxObs_NP28
     */
    virtual double AuxObs_NP28() const
    {
        return 0.0;
    }
    
    /**
     * @brief Auxiliary observable AuxObs_NP29
     * @return AuxObs_NP29
     */
    virtual double AuxObs_NP29() const
    {
        return 0.0;
    }
    
    /**
     * @brief Auxiliary observable AuxObs_NP30
     * @return AuxObs_NP30
     */
    virtual double AuxObs_NP30() const
    {
        return 0.0;
    }
    
    virtual int OutputOrder() const;                                    //AG:added
      
    ////////////////////////////////////////////////////////////////////////
protected:
    StandardModel trueSM;
};

#endif	/* NPBASE_H */

