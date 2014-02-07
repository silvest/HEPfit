/*
 * Copyright (C) 2013-2014 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef NPBASE_H
#define	NPBASE_H

#include <StandardModel.h>

/**
 * @class NPbase
 * @brief A base model class for general new physics contributions to electroweak 
 * precision observables. 
 * @ingroup NewPhysics
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 *  
 * @details This is a Model class containing the basic structure to compute
 * new physics contributions to the electroweak precision observables, where
 * the following quantities are handled:
 *
 * @li @f$S@f$, @f$T@f$ and @f$U@f$&nbsp;&nbsp;
 * (with obliqueS(), obliqueT() and obliqueU()),
 * @li @f$\hat{S}@f$, @f$\hat{T}@f$, @f$\hat{U}@f$,
 * @f$V@f$, @f$W@f$, @f$X@f$ and @f$Y@f$&nbsp;&nbsp;
 * (with obliqueShat(), obliqueThat(), obliqueUhat(),
 * obliqueV(), obliqueW(), obliqueX() and obliqueY()),
 * @li @f$\Delta G@f$&nbsp;&nbsp; (with DeltaGF()),
 * @li @f$\delta g_V^f@f$&nbsp;&nbsp;(with deltaGVl() and deltaGVq()),
 * @li @f$\delta g_A^f@f$&nbsp;&nbsp;(with deltaGAl() and deltaGAq()).
 *
 * In the model classes that are inherited from the current class (see the
 * inheritance diagram above), some of these methods are reimplemented to
 * account for the details of more specific scenarios.
 */
class NPbase : public StandardModel {
public:

    /**
     * @brief Th default constructor.
     */
    NPbase();

    /**
     * @brief A method to fetch the name of %NPbase.
     * @return the name of the model as a string
     */
    virtual std::string ModelName() const
    {
        return "NPbase";
    }

    /**
     * @brief A method to initialize %NPbase.
     * @return a boolean that is true if model initialization is successful
     */
    virtual bool InitializeModel();

    /**
     * @brief A method to initialize the model parameters.
     * @param[in] DPars a map of the parameters that are being updated in the Monte Carlo run
     * (including parameters that are varied and those that are held constant)
     * @return a boolean that is true if the execution is successful
     */
    virtual bool Init(const std::map<std::string, double>& DPars);
    
    /**
     * @brief The update method for %NPbase.
     * @details This method updates all the model parameters with giving DPars.
     * @param[in] DPars a map of the parameters that are being updated in the Monte Carlo run
     * (including parameters that are varied and those that are held constant)
     * @return a boolean that is true if the execution is successful
     */
    virtual bool Update(const std::map<std::string, double>& DPars);
    
    /**
     * @brief A method to check if all the mandatory parameters for %NPbase
     * have been provided in the model configuration file.
     * @param[in] DPars a map of the parameters that are being updated in the Monte Carlo run
     * (including parameters that are varied and those that are held constant)
     * @return a boolean that is true if the execution is successful
     */
    virtual bool CheckParameters(const std::map<std::string, double>& DPars);

    /**
     * @brief A method to set a flag of %NPbase.
     * @param[in] name name of a model flag
     * @param[in] value the boolean to be assigned to the flag specified by name
     * @return a boolean that is true if the execution is successful
     */
    virtual bool setFlag(const std::string name , const bool value);
    
    /**
     * @brief A method to check the sanity of the set of model flags.
     * @return a boolean that is true if the set of model flags is sane
     */
    virtual bool CheckFlags() const;

    
    ////////////////////////////////////////////////////////////////////////

    /**
     * @brief The oblique parameter \f$S\f$.
     * @return the value of the oblique parameter @f$S@f$
     */
    virtual double obliqueS() const;

    /**
     * @brief The oblique parameter \f$T\f$.
     * @return the value of the oblique parameter @f$T@f$
     */
    virtual double obliqueT() const;

    /**
     * @brief The oblique parameter \f$U\f$.
     * @return the value of the oblique parameter @f$U@f$
     */
    virtual double obliqueU() const;

    /**
     * @brief The oblique parameter \f$\hat{S}\f$.
     * @return the value of the oblique parameter 
     * \f$\displaystyle\hat{S}=\frac{\alpha}{4s^2_W}S\f$
     */
    virtual double obliqueShat() const;

    /**
     * @brief The oblique parameter \f$\hat{T}\f$.
     * @return the value of the oblique parameter \f$\hat{T}=\alpha T\f$
     */
    virtual double obliqueThat() const;

    /**
     * @brief The oblique parameter \f$\hat{U}\f$.
     * @return the value of the oblique parameter 
     * \f$\displaystyle\hat{U}=-\frac{\alpha}{4s^2_W}U\f$
     */
    virtual double obliqueUhat() const;

    /**
     * @brief The oblique parameter \f$V\f$.
     * @return the value of the oblique parameter \f$V\f$
     */
    virtual double obliqueV() const;

    /**
     * @brief The oblique parameter \f$W\f$.
     * @return the value of the oblique parameter \f$W\f$
     */
    virtual double obliqueW() const;

    /**
     * @brief The oblique parameter \f$X\f$.
     * @return the value of the oblique parameter \f$X\f$
     */
    virtual double obliqueX() const;

    /**
     * @brief The oblique parameter \f$Y\f$.
     * @return the value of the oblique parameter \f$Y\f$
     */
    virtual double obliqueY() const;
    
    /**
     * @brief New physics contribution to the Fermi constant.
     * @details The new physics contribution @f$\Delta G@f$ is defined as
     * @f[
     * G_\mu = G_{\mu,\mathrm{SM}}(1+\Delta G)\,,
     * @f]
     * where @f$G_{\mu,\mathrm{SM}}@f$ denotes the Fermi constant in the SM,
     * while @f$G_\mu@f$ is the experimentl value measured through muon decays.
     * @return @f$\Delta G@f$
     */
    virtual double DeltaGF() const;

    /**
     * @brief New physics contribution to @f$g_V^l@f$.
     * @details
     * @f[
     * \delta g_V^l
     * = g_{V,\mathrm{SM}}^l\,\frac{\alpha(M_Z^2)\, T'}{2}
     * + \big( g_{V,\mathrm{SM}}^l - g_{A,\mathrm{SM}}^l \big)\,
     * \frac{\alpha(M_Z^2)}{4s_W^2\,(c_W^2-s_W^2)}
     * \left( S' - 4\,c_W^2s_W^2\, T' \right).
     * @f]
     *
     * See, e.g., @cite Ciuchini:2013pca.
     * @param[in] l name of a lepton (see StandardModel::lepton)
     * @return @f$\delta g_V^l@f$
     */
    virtual double deltaGVl(StandardModel::lepton l) const;
    
    /**
     * @brief New physics contribution to @f$g_V^q@f$.
     * @details
     * @f[
     * \delta g_V^q
     * = g_{V,\mathrm{SM}}^q\,\frac{\alpha(M_Z^2)\, T'}{2}
     * + \big( g_{V,\mathrm{SM}}^q - g_{A,\mathrm{SM}}^q \big)\,
     * \frac{\alpha(M_Z^2)}{4s_W^2\,(c_W^2-s_W^2)}
     * \left( S' - 4\,c_W^2s_W^2\, T' \right).
     * @f]
     *
     * See, e.g., @cite Ciuchini:2013pca.
     * @param[in] q name of a quark (see QCD::quark)
     * @return @f$\delta g_V^q@f$
     */
    virtual double deltaGVq(StandardModel::quark q) const;

    /**
     * @brief New physics contribution to @f$g_A^l@f$.
     * @details
     * @f[
     * \delta g_A^l = g_{A,\mathrm{SM}}^l \frac{\alpha(M_Z^2)\, T'}{2}\,.
     * @f]
     *
     * See, e.g., @cite Ciuchini:2013pca.
     * @param[in] l name of a lepton (see StandardModel::lepton)
     * @return @f$\delta g_A^l@f$
     */   
    virtual double deltaGAl(StandardModel::lepton l) const;

    /**
     * @brief New physics contribution to @f$g_A^q@f$.
     * @details
     * @f[
     * \delta g_A^q = g_{A,\mathrm{SM}}^q \frac{\alpha(M_Z^2)\, T'}{2}\,.
     * @f]
     *
     * See, e.g., @cite Ciuchini:2013pca.
     * @param[in] q name of a quark (see QCD::quark)
     * @return @f$\delta g_A^q@f$
     */ 
    virtual double deltaGAq(StandardModel::quark q) const;

    
    ////////////////////////////////////////////////////////////////////////
protected:

    /**
     * @brief A method to set the value of a parameter of %NPbase.
     * @param[in] name name of a model parameter
     * @param[in] value the value to be assigned to the parameter specified by name
     */
    virtual void setParameter(const std::string name, const double& value);

    
};

#endif	/* NPBASE_H */

