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
 * @brief A base class for general new physics corrections to electroweak precision
 * observables
 * @ingroup NewPhysics
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details
 */
class NPbase : public StandardModel {
public:

    static const int NNPbaseflags = 1;
    static const std::string NPbaseflags[NNPbaseflags];

    /**
     * @brief Constructor. 
     */
    NPbase();

    virtual std::string ModelName() const
    {
        return "NPbase";
    }

    virtual bool InitializeModel();
    virtual void setEWSMflags(EWSM& myEWSM);

    virtual bool Init(const std::map<std::string, double>& DPars);
    virtual bool Update(const std::map<std::string, double>& DPars);
    virtual bool CheckParameters(const std::map<std::string, double>& DPars);

    virtual bool setFlag(const std::string, const bool&);
    virtual bool CheckFlags() const;

    
    ////////////////////////////////////////////////////////////////////////

    bool IsFlagFixSMcontribution() const
    {
        return FlagFixSMcontribution;
    }


    ////////////////////////////////////////////////////////////////////////

    /**
     * @return the value of the oblique parameter @f$S@f$
     */
    virtual double obliqueS() const;

    /**
     * @return the value of the oblique parameter @f$T@f$
     */
    virtual double obliqueT() const;

    /**
     * @return the value of the oblique parameter @f$U@f$
     */
    virtual double obliqueU() const;

    /**
     * @return the value of the oblique parameter @f$\hat{S}=\frac{\alpha}{4\sin^2{\theta_W}}S@f$
     */
    virtual double obliqueShat() const;

    /**
     * @return the value of the oblique parameter @f$\hat{T}=\alpha T@f$
     */
    virtual double obliqueThat() const;

    /**
     * @return the value of the oblique parameter @f$\hat{U}=-\frac{\alpha}{4\sin^2{\theta_W}}U@f$
     */
    virtual double obliqueUhat() const;

    /**
     * @return the value of the oblique parameter @f$V@f$
     */
    virtual double obliqueV() const;

    /**
     * @return the value of the oblique parameter @f$W@f$
     */
    virtual double obliqueW() const;

    /**
     * @return the value of the oblique parameter @f$X@f$
     */
    virtual double obliqueX() const;

    /**
     * @return the value of the oblique parameter @f$Y@f$
     */
    virtual double obliqueY() const;
    
    /**
     * @return the new physics contribution to the Fermi constant, @f$G_{\mu}=G_{\mu}^{\rm SM}(1+\Delta GF)@f$
     */
    virtual double DeltaGF() const;

    /**
     * @param[in] l name of a lepton
     * @return the new physics contribution to neutral-current vector coupling @f$\delta g_V^l@f$
     */
    virtual double deltaGVl(StandardModel::lepton l) const;
    
    /**
     * @param[in] q name of a quark
     * @return the new physics contribution to neutral-current vector coupling @f$\delta g_V^q@f$
     */
    virtual double deltaGVq(StandardModel::quark q) const;

    /**
     * @param[in] l name of a lepton
     * @return the new physics contribution to neutral-current axial-vector coupling @f$\delta g_A^l@f$
     */   
    virtual double deltaGAl(StandardModel::lepton l) const;

    /**
     * @param[in] q name of a quark
     * @return the new physics contribution to neutral-current axial-vector coupling @f$\delta g_A^q@f$
     */ 
    virtual double deltaGAq(StandardModel::quark q) const;

    
    ////////////////////////////////////////////////////////////////////////
protected:
    virtual void setParameter(const std::string name, const double& value);


    ////////////////////////////////////////////////////////////////////////
private:
    bool FlagFixSMcontribution;
    
};

#endif	/* NPBASE_H */

