/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef MESON_H
#define	MESON_H

class QCD;

#include <stdexcept>
#include <vector>
#include <map>
#include <functional>
#include "Particle.h"
#include "boost/lexical_cast.hpp"

/**
 * @class Meson
 * @ingroup StandardModel
 * @brief A class for mesons. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to define a meson and three of its
 * characteristics: mass, lifetime and decay constant. For selected mesons,
 * also the Gegenbauer moments are defined.
 * This class inherits the public access members of the Particle class.
 * 
 *
 * @anchor MesonParameters
 * <h3>%Model parameters</h3>
 *
 * The model parameters of %Meson are summarized below:
 * <table class="model">
 * <tr>
 *   <th>Label</th>
 *   <th>LaTeX symbol</th>
 *   <th>Description</th>
 * </tr>
 * <tr>
 *   <td class="mod_name">%MP0</td>
 *   <td class="mod_symb">@f$M_{\pi^0}@f$</td>
 *   <td class="mod_desc">The mass of the \f$ \pi^0 \f$ meson in GeV.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%tP0</td>
 *   <td class="mod_symb">@f$\tau_{\pi^0}@f$</td>
 *   <td class="mod_desc">The lifetime of the \f$ \pi^0 \f$ meson in \f$\mathrm{ps}^{-1} \f$.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%FP0</td>
 *   <td class="mod_symb">@f$F_{\pi^0}@f$</td>
 *   <td class="mod_desc">The decay constant of the \f$ \pi^0 \f$ meson in GeV.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%MPp</td>
 *   <td class="mod_symb">@f$M_{\pi^\pm}@f$</td>
 *   <td class="mod_desc">The mass of the \f$ \pi^\pm \f$ meson in GeV.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%tPp</td>
 *   <td class="mod_symb">@f$\tau_{\pi^\pm}@f$</td>
 *   <td class="mod_desc">The lifetime of the \f$ \pi^\pm \f$ meson in \f$\mathrm{ps}^{-1} \f$.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%FPp</td>
 *   <td class="mod_symb">@f$F_{\pi^\pm}@f$</td>
 *   <td class="mod_desc">The decay constant of the \f$ \pi^{\pm} \f$ meson in GeV.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%MK0</td>
 *   <td class="mod_symb">@f$M_{K^0}@f$</td>
 *   <td class="mod_desc">The mass of the \f$ K^0 \f$ meson in GeV.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%tKl</td>
 *   <td class="mod_symb">@f$\tau_{K_L}@f$</td>
 *   <td class="mod_desc">The lifetime of the \f$ K_L \f$ meson in \f$\mathrm{ps}^{-1} \f$.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%MKp</td>
 *   <td class="mod_symb">@f$M_{K^\pm}@f$</td>
 *   <td class="mod_desc">The mass of the \f$ K^\pm \f$ meson in GeV.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%tKp</td>
 *   <td class="mod_symb">@f$\tau_{K^\pm}@f$</td>
 *   <td class="mod_desc">The lifetime of the \f$ K^\pm \f$ meson in \f$\mathrm{ps}^{-1} \f$.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%FK</td>
 *   <td class="mod_symb">@f$F_{K}@f$</td>
 *   <td class="mod_desc">The decay constant of the \f$ K^{0,\pm} \f$ meson in GeV.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%alpha1kp, %alpha2kp</td>
 *   <td class="mod_symb">@f$\alpha_1(K^+), \alpha_2(K^+)@f$</td>
 *   <td class="mod_desc">The Gegenbauer coefficients for the @f$K^+@f$ meson.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%MD</td>
 *   <td class="mod_symb">@f$M_{D^0}@f$</td>
 *   <td class="mod_desc">The mass of the \f$ D^0 \f$ meson in GeV.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%tD</td>
 *   <td class="mod_symb">@f$\tau_{D^0}@f$</td>
 *   <td class="mod_desc">The lifetime of the \f$ D^0 \f$ meson in \f$\mathrm{ps}^{-1} \f$.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%FD</td>
 *   <td class="mod_symb">@f$F_{D^0}@f$</td>
 *   <td class="mod_desc">The decay constant of the \f$ D^0 \f$ meson in GeV.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%MDp</td>
 *   <td class="mod_symb">@f$M_{D^\pm}@f$</td>
 *   <td class="mod_desc">The mass of the \f$ D^\pm \f$ meson in GeV.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%tDp</td>
 *   <td class="mod_symb">@f$\tau_{D^\pm}@f$</td>
 *   <td class="mod_desc">The lifetime of the \f$ D^\pm \f$ meson in \f$\mathrm{ps}^{-1} \f$.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%FDp</td>
 *   <td class="mod_symb">@f$F_{D^\pm}@f$</td>
 *   <td class="mod_desc">The decay constant of the \f$ D^\pm \f$ meson in GeV.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%MBd</td>
 *   <td class="mod_symb">@f$M_{B_d}@f$</td>
 *   <td class="mod_desc">The mass of the \f$ B_d \f$ meson in GeV.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%tBd</td>
 *   <td class="mod_symb">@f$\tau_{B_d}@f$</td>
 *   <td class="mod_desc">The lifetime of the \f$ B_d \f$ meson in \f$\mathrm{ps}^{-1} \f$.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%MBp</td>
 *   <td class="mod_symb">@f$M_{B^\pm}@f$</td>
 *   <td class="mod_desc">The mass of the \f$ B^\pm \f$ meson in GeV.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%tBp</td>
 *   <td class="mod_symb">@f$\tau_{B^\pm}@f$</td>
 *   <td class="mod_desc">The lifetime of the \f$ B^\pm \f$ meson in \f$\mathrm{ps}^{-1} \f$.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%FBsoFBd</td>
 *   <td class="mod_symb">@f$F_{B_s}/F_{B_d}@f$</td>
 *   <td class="mod_desc">The ratio of decay constants of the \f$ B_s \f$ and \f$ B_d \f$ mesons.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%MBs</td>
 *   <td class="mod_symb">@f$M_{B_s}@f$</td>
 *   <td class="mod_desc">The mass of the \f$ B_s \f$ meson in GeV.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%tBs</td>
 *   <td class="mod_symb">@f$\tau_{B_s}@f$</td>
 *   <td class="mod_desc">The lifetime of the \f$ B_s \f$ meson in \f$\mathrm{ps}^{-1} \f$.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%FBs</td>
 *   <td class="mod_symb">@f$F_{B_s}@f$</td>
 *   <td class="mod_desc">The decay constant of the \f$ B_s \f$ meson in GeV.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%lambdaB</td>
 *   <td class="mod_symb">@f$\Lambda_{B,+}@f$</td>
 *   <td class="mod_desc">The integrated leading twist light-cone distribution amplitudes of the B meson divided by the integral variable.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%DGs_Gs</td>
 *   <td class="mod_symb">@f$ \Delta \Gamma_s /\Gamma_s @f$</td>
 *   <td class="mod_desc">The oscillation parameter for the \f$B_s\f$ meson.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%MBc</td>
 *   <td class="mod_symb">@f$M_{B_c}@f$</td>
 *   <td class="mod_desc">The mass of the \f$ B_c \f$ meson in GeV.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%tBc</td>
 *   <td class="mod_symb">@f$\tau_{B_c}@f$</td>
 *   <td class="mod_desc">The lifetime of the \f$ B_c \f$ meson in \f$\mathrm{ps}^{-1} \f$.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%FBc</td>
 *   <td class="mod_symb">@f$F_{B_c}@f$</td>
 *   <td class="mod_desc">The decay constant of the \f$ B_c \f$ meson in GeV.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%Mphi</td>
 *   <td class="mod_symb">@f$M_{\phi}@f$</td>
 *   <td class="mod_desc">The mass of the \f$ \phi \f$ meson in GeV.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%tphi</td>
 *   <td class="mod_symb">@f$\tau_{\phi}@f$</td>
 *   <td class="mod_desc">The lifetime of the \f$ \phi \f$ meson in \f$\mathrm{ps}^{-1} \f$.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%Fphi</td>
 *   <td class="mod_symb">@f$F_{\phi}@f$</td>
 *   <td class="mod_desc">The decay constant of the \f$ \phi \f$ meson in GeV.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%Fphip</td>
 *   <td class="mod_symb">@f$F_{\phi}^{\perp}@f$</td>
 *   <td class="mod_desc">The decay constant of a transversely polarized \f$ \phi \f$ meson in GeV.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%alpha2phi</td>
 *   <td class="mod_symb">@f$\alpha_2(\phi)@f$</td>
 *   <td class="mod_desc">The Gegenbauer coefficient for the @f$\phi@f$ meson.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%MKstar</td>
 *   <td class="mod_symb">@f$M_{K^{*0}}@f$</td>
 *   <td class="mod_desc">The mass of the \f$ K^{*0} \f$ meson in GeV.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%MKstarP</td>
 *   <td class="mod_symb">@f$M_{K^{*\pm}}@f$</td>
 *   <td class="mod_desc">The mass of the \f$ K^{*\pm} \f$ meson in GeV.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%tKstar</td>
 *   <td class="mod_symb">@f$\tau_{K^*}@f$</td>
 *   <td class="mod_desc">The lifetime of the \f$ K^* \f$ meson in \f$\mathrm{ps}^{-1} \f$.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%FKstar</td>
 *   <td class="mod_symb">@f$F_{K^*}@f$</td>
 *   <td class="mod_desc">The decay constant of the \f$ K^* \f$ meson in GeV.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%FKstarp</td>
 *   <td class="mod_symb">@f$F_{K^*}^{\perp}@f$</td>
 *   <td class="mod_desc">The decay constant of a transversely polarized \f$ K^* \f$ meson in GeV.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%alpha1kst, %alpha2kst</td>
 *   <td class="mod_symb">@f$\alpha_1(\bar{K}^*), \alpha_2(\bar{K}^*)@f$</td>
 *   <td class="mod_desc">The Gegenbauer coefficients for the @f$\bar{K}^*@f$ meson.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%MDstarP</td>
 *   <td class="mod_symb">@f$M_{D^{*\pm}}@f$</td>
 *   <td class="mod_desc">The mass of the \f$ K^{*\pm} \f$ meson in GeV.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%tDstarP</td>
 *   <td class="mod_symb">@f$\tau_{D^{*\pm}}@f$</td>
 *   <td class="mod_desc">The lifetime of the \f$ D^{*\pm} \f$ meson in \f$\mathrm{ps}^{-1} \f$.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%FDstarP</td>
 *   <td class="mod_symb">@f$F_{D^{*\pm}}@f$</td>
 *   <td class="mod_desc">The decay constant of the \f$ D^{*\pm} \f$ meson in GeV.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%Mrho</td>
 *   <td class="mod_symb">@f$M_{\rho}@f$</td>
 *   <td class="mod_desc">The mass of the \f$ \rho \f$ meson in GeV.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%trho</td>
 *   <td class="mod_symb">@f$\tau_{\rho}@f$</td>
 *   <td class="mod_desc">The lifetime of the \f$ \rho \f$ meson in \f$\mathrm{ps}^{-1} \f$.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%Frho</td>
 *   <td class="mod_symb">@f$F_{\rho}@f$</td>
 *   <td class="mod_desc">The decay constant of the \f$ \rho \f$ meson in GeV.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%MrhoP</td>
 *   <td class="mod_symb">@f$M_{\rho^{\pm}}@f$</td>
 *   <td class="mod_desc">The mass of the \f$ \rho^{\pm} \f$ meson in GeV.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%trhoP</td>
 *   <td class="mod_symb">@f$\tau_{\rho^{\pm}}@f$</td>
 *   <td class="mod_desc">The lifetime of the \f$ \rho^{\pm} \f$ meson in \f$\mathrm{ps}^{-1} \f$.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%FrhoP</td>
 *   <td class="mod_symb">@f$F_{\rho^{\pm}}@f$</td>
 *   <td class="mod_desc">The decay constant of the \f$ \rho^{\pm} \f$ meson in GeV.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%Momega</td>
 *   <td class="mod_symb">@f$M_{\omega}@f$</td>
 *   <td class="mod_desc">The mass of the \f$ \omega \f$ meson in GeV.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%tomega</td>
 *   <td class="mod_symb">@f$\tau_{\omega}@f$</td>
 *   <td class="mod_desc">The lifetime of the \f$ \omega \f$ meson in \f$\mathrm{ps}^{-1} \f$.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%Fomega</td>
 *   <td class="mod_symb">@f$F_{\omega}@f$</td>
 *   <td class="mod_desc">The decay constant of the \f$ \omega \f$ meson in GeV.</td>
 * </tr>
 * </table>
 * 
 */
class Meson : public Particle {
public:

    /**
     * @brief The default constructor.
     */
    Meson();

    /**
     * @brief Constructor.
     * @param[in] mass the mass of the meson in GeV
     * @param[in] lifetime the lifetime of the meson in \f$ \mathrm{ps}^{-1} \f$
     * @param[in] decayconst the decay constant of the meson in GeV
     * @param[in] lambdaM the first moment of the LCDA
     * @param[in] gegenalpha1 first Gegenbauer moment of LCDA
     * @param[in] gegenalpha2 second Gegenbauer moment of LCDA
     */
    Meson(double mass, double lifetime, double decayconst, double lambdaM, double gegenalpha1, double gegenalpha2);

    /**
     * @brief The default destructor.
     */
    virtual ~Meson();
    
    void ModelParameterMapInsert(std::map< std::string, std::reference_wrapper<const double> >& ModelParamMap);
    
    std::vector<std::string> parameterList(std::string name_i);
    
    bool setParameter(std::string name_i, double value); 
    
    void initializeParameters();
    
    /**
     * @brief A get method for the lifetime of the meson.
     * @return the lifetime of the meson in \f$ \mathrm{ps}^{-1} \f$
     */
    double getLifetime() const
    {
        return lifetime;
    }

    /**
     * @brief A get method for the decay constant of the meson.
     * @return the decay constant of the meson in GeV
     */
    const double& getDecayconst() const
    {
        return decayconst;
    }
    
    /**
     * @brief A set method for the decay constant of the meson.
     * @param[in] decayconst the decay constant of the meson in GeV
     */
    void setDecayconst(double decayconst)
    {
        this->decayconst = decayconst;
    }
    
    /**
     * @brief A get method for the perpendicular decay constant of a vector meson.
     * @return the decay constant of the meson in GeV
     */
    const double& getDecayconst_p() const
    {
        return decayconst_p;
    }

    /**
     * @brief A method to compute the width of the meson from its lifetime.
     * @return the width of the meson in GeV
     */
    double computeWidth() const;
    
    /**
     * @brief A get method to get the Gegenbaur coefficient.
     * @param[in] the order of the Gegenbaur coefficient
     * @return the Gegenbaur coefficient
     */
    const double& getGegenalpha(int i) const
    {
        if (i >= 0 && i < 2)
            return gegenalpha[i];
        else 
            throw std::runtime_error("Meson::getGegenalpha(" + boost::lexical_cast<std::string>(i) + "): index out of range");
    }

    const double& getLambdaM() const
    {
        return lambdaM;
    }
    
    void setDgamma_gamma(double Dgamma_gamma){
        this->Dgamma_gamma = Dgamma_gamma;
    }
    
    const double& getDgamma_gamma() const
    {
        return Dgamma_gamma;
    }
    
    /**
     * @brief A get method to get the name of the meson
     * @return the the name of the meson
     */
    std::string getName() const
    {
        return name;
    }
    
    /**
     * @brief A set method to set the name of the meson
     * @param[in] name_i the the name of the meson
     */
    void setName(std::string name_i)
    {
        this->name = name_i;
    }
    
    double getFBsoFBd() const
    {
        return FBsoFBd;
    }

private:
    double decayconst; ///< The decay constant of the meson.
    double decayconst_p; ///< The perpendicular decay constant of a vector meson.
    double lifetime; ///< The lifetime of the meson.
    double gegenalpha[2]; ///< Gegenbauer moments 
    double lambdaM; ///< First moment of LCDA
    double Dgamma_gamma; ///< Dgamma/gamma for neutral mesons
    double FBsoFBd;
};

#endif	/* MESON_H */

