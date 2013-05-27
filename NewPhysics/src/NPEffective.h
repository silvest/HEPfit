/* 
 * Copyright (C) 2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef NPEFFECTIVE_H
#define	NPEFFECTIVE_H

#include <stdexcept>
#include <StandardModel.h>

/**
 * @class NPEffective
 * @brief A class for new physics with an effective Lagrangian. 
 * @ingroup NewPhysics
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details  
 */
class NPEffective : public StandardModel {
public:
    static const int NNPEffectiveVars = 11;
    static const std::string NPEffectiveVars[NNPEffectiveVars];

    /**
     * @brief NPEffective constructor. 
     */
    NPEffective();

    virtual std::string ModelName() const 
    {
        return "NPEffective";
    }

    virtual bool Update(const std::map<std::string, double>& DPars);
    virtual bool Init(const std::map<std::string, double>& DPars);    
    virtual bool CheckParameters(const std::map<std::string, double>& DPars);

    virtual bool InitializeModel();  
    virtual void SetEWSMflags(EWSM& myEWSM);    

    virtual bool SetFlag(const std::string, const bool&); 


    ////////////////////////////////////////////////////////////////////////

    virtual double v() const;

    virtual double Mw_tree() const;

    virtual double DeltaGF() const;


    ////////////////////////////////////////////////////////////////////////     

    /**
     * @return Oblique parameter S.
     */
    virtual double obliqueS() const;

    /**
     * @return Oblique parameter T.
     */
    virtual double obliqueT() const;

    /**
     * @return Oblique parameter U.
     */
    virtual double obliqueU() const;


    ////////////////////////////////////////////////////////////////////////    

    double deltaGLl(StandardModel::lepton l) const;

    double deltaGLq(StandardModel::quark q) const;
    
    double deltaGRl(StandardModel::lepton l) const;

    double deltaGRq(StandardModel::quark q) const;

    virtual double deltaGVl(StandardModel::lepton l) const;
    
    virtual double deltaGVq(StandardModel::quark q) const;
    
    virtual double deltaGAl(StandardModel::lepton l) const;
    
    virtual double deltaGAq(StandardModel::quark q) const;


    ////////////////////////////////////////////////////////////////////////

    virtual double epsilon1() const;

    virtual double epsilon2() const;

    virtual double epsilon3() const;

    virtual double epsilonb() const;


    ////////////////////////////////////////////////////////////////////////

    /**
     * @return The W boson mass.
     */
    virtual double Mw() const;

    /**
     * @return @f$M_W^2/M_Z^2@f$.
     */
    virtual double cW2() const;

    /**
     * @return @f$1-M_W^2/M_Z^2@f$.
     */
    virtual double sW2() const;

    /**
     * @return The total width of the W boson.
     */
    virtual double GammaW() const;
    
    
    ////////////////////////////////////////////////////////////////////////    

protected:
    double cWB, cH, cLL, cHLp, cHQp, cHL, cHQ, cHE, cHU, cHD, LambdaNP;
    virtual void SetParameter(const std::string name, const double& value);
    

    ////////////////////////////////////////////////////////////////////////   

private:

};

#endif	/* NPEFFECTIVE_H */

