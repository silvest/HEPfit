/*
 * Copyright (C) 2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef NPBASE_H
#define	NPBASE_H

#include <StandardModel.h>

/**
 * @class NPbase
 * @brief A base class for ...
 * @ingroup NewPhysics
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details
 */
class NPbase : public StandardModel {
public:

    static const int NNPbaseflags = 1;
    static const std::string NPbaseflags[NNPbaseflags];

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
     * @return @f$\epsilon_1@f$.
     */
    virtual double epsilon1() const;

    /**
     * @return @f$\epsilon_2@f$.
     */
    virtual double epsilon2() const;

    /**
     * @return @f$\epsilon_3@f$.
     */
    virtual double epsilon3() const;

    /**
     * @return @f$\epsilon_b@f$.
     */
    virtual double epsilonb() const;

    /**
     * @return Oblique parameter @f$S@f$.
     */
    virtual double obliqueS() const;

    /**
     * @return Oblique parameter @f$T@f$.
     */
    virtual double obliqueT() const;

    /**
     * @return Oblique parameter @f$U@f$.
     */
    virtual double obliqueU() const;

    /**
     * @return Oblique parameter @f$\hat{S}@f$.
     */
    virtual double obliqueShat() const;

    /**
     * @return Oblique parameter @f$\hat{T}@f$.
     */
    virtual double obliqueThat() const;

    /**
     * @return Oblique parameter @f$\hat{U}@f$.
     */
    virtual double obliqueUhat() const;

    /**
     * @return Oblique parameter @f$V@f$.
     */
    virtual double obliqueV() const;

    /**
     * @return Oblique parameter @f$W@f$.
     */
    virtual double obliqueW() const;

    /**
     * @return Oblique parameter @f$X@f$.
     */
    virtual double obliqueX() const;

    /**
     * @return Oblique parameter @f$Y@f$.
     */
    virtual double obliqueY() const;
    
    /**
     * @brief NP contribution to the Fermi constant @f$G_{\mu}=G_{\mu}^{\rm SM}(1+DeltaGF)@f$.
     * @return
     */
    virtual double DeltaGF() const;

    virtual double deltaGVl(StandardModel::lepton l) const;

    virtual double deltaGVq(StandardModel::quark q) const;

    virtual double deltaGAl(StandardModel::lepton l) const;

    virtual double deltaGAq(StandardModel::quark q) const;

    
    ////////////////////////////////////////////////////////////////////////
protected:
    virtual void setParameter(const std::string name, const double& value);


    ////////////////////////////////////////////////////////////////////////
private:
    bool FlagFixSMcontribution;
    
};

#endif	/* NPBASE_H */

