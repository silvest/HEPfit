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
 * characteristics: mass, lifetime and decay constant. 
 * This class inherits the public access members of the Particle class.
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

