/* 
 * Copyright (C) 2018 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef LOOPMEDIATORS_H
#define LOOPMEDIATORS_H

#include "StandardModel.h"
#include "gslpp.h"
#include "LoopMediatorsMatching.h"
#include "ThObservable.h"

/**
 * @class LoopMediators
 * @brief Model for NP contributions to flavour.
 */
class LoopMediators: public StandardModel {
public:

    static const int NLoopMediatorsvars = 10;

    static const std::string LoopMediatorsvars[NLoopMediatorsvars];
    
    /**
     * @brief FlavourWilsonCoefficient constructor
     */
    LoopMediators();
    
    /**
     * @brief FlavourWilsonCoefficient destructor
     */
    ~LoopMediators();
    
    virtual bool InitializeModel();
    
    virtual bool Init(const std::map<std::string, double>& DPars);
    
    virtual bool PreUpdate();
    
    virtual bool Update(const std::map<std::string, double>& DPars);
    
    virtual bool PostUpdate();
    
    virtual bool CheckParameters(const std::map<std::string, double>& DPars);
    
    virtual bool setFlag(const std::string name, const bool value);
    
    double F9(double x, double y);
    double F7(double x);
    double F7t(double x);
    double G7(double x);
    double G7t(double x);
    
    virtual LoopMediatorsMatching& getMatching() const
    {
        return LoopMediatorsM.getObj();
    }
    
    /**
     *
     * @return \f$ C_1 $\f
     */
    gslpp::complex getC1() const {
        return C1;
    }
    
    /**
     *
     * @return \f$ C_2 $\f
     */
    gslpp::complex getC2() const {
        return C2;
    }
    
    /**
     *
     * @return \f$ C_3 $\f
     */
    gslpp::complex getC3() const {
        return C3;
    }
    
    /**
     *
     * @return \f$ C_4 $\f
     */
    gslpp::complex getC4() const {
        return C4;
    }
    
    /**
     *
     * @return \f$ C_5 $\f
     */
    gslpp::complex getC5() const {
        return C5;
    }
    
    /**
     *
     * @return \f$ C_1' $\f
     */
    gslpp::complex getC1p() const {
        return C1p;
    }
    
    /**
     *
     * @return \f$ C_2' $\f
     */
    gslpp::complex getC2p() const {
        return C2p;
    }
    
    /**
     *
     * @return \f$ C_3' $\f
     */
    gslpp::complex getC3p() const {
        return C3p;
    }
    
    /**
     *
     * @return \f$ C_7 $\f
     */
    gslpp::complex getC7() const {
        return C7;
    }
    
    /**
     *
     * @return \f$ C_8 $\f
     */
    gslpp::complex getC8() const {
        return C8;
    }
    
    /**
     *
     * @return \f$ C_9 $\f
     */
    gslpp::complex getC9() const {
        return C9;
    }
    
    /**
     *
     * @return \f$ C_10 $\f
     */
    gslpp::complex getC10() const {
        return C10;
    }
    
    /**
     *
     * @return \f$ C_S $\f
     */
    gslpp::complex getCS() const {
        return CS;
    }
    
    /**
     *
     * @return \f$ C_P $\f
     */
    gslpp::complex getCP() const {
        return CP;
    }
    
    /**
     *
     * @return \f$ C_7' $\f
     */
    gslpp::complex getC7p() const {
        return C7p;
    }
    
    /**
     *
     * @return \f$ C_8' $\f
     */
    gslpp::complex getC8p() const {
        return C8p;
    }
    
    /**
     *
     * @return \f$ C_9' $\f
     */
    gslpp::complex getC9p() const {
        return C9p;
    }
    
    /**
     *
     * @return \f$ C_10' $\f
     */
    gslpp::complex getC10p() const {
        return C10p;
    }
    
    /**
     *
     * @return \f$ C_S' $\f
     */
    gslpp::complex getCSp() const {
        return CSp;
    }
    
    /**
     *
     * @return \f$ C_P' $\f
     */
    gslpp::complex getCPp() const {
        return CPp;
    }
    
    /**
     *
     * @return \f$ \Delta(a_\mu) $\f
     */
    double getDeltaamu() const {
        return Deltaamu;
    }
    
    /**
     *
     * @return the matching scale of the Wilson coefficients
     */
    double getWCscale() const {
        return WCscale;
    }
    
protected: 
    
    virtual void setParameter(const std::string, const double&);
    mutable Matching<LoopMediatorsMatching,LoopMediators> LoopMediatorsM;

private:
    
    gslpp::complex C1;
    gslpp::complex C2;
    gslpp::complex C3;
    gslpp::complex C4;
    gslpp::complex C5;
    
    gslpp::complex C1p;
    gslpp::complex C2p;
    gslpp::complex C3p;

    gslpp::complex C7;
    gslpp::complex C8;
    gslpp::complex C9;
    gslpp::complex C10;
    gslpp::complex CS;
    gslpp::complex CP;
    
    gslpp::complex C7p;
    gslpp::complex C8p;
    gslpp::complex C9p;
    gslpp::complex C10p;
    gslpp::complex CSp;
    gslpp::complex CPp;
    
    double Deltaamu;
    
    double GammaL;
    double GammaR;
    double GammamuL;
    double GammamuR;
    double lambdaE;
    double mphi;
    double yD;
    double yE;
    double charge;
    
    double WCscale;
    
      
};



class Deltaamu : public ThObservable {
public:
    /**
     * @brief Constructor.
     */
   Deltaamu(const StandardModel& SM_i);
     
   ~Deltaamu();

   /**
     * @brief Two positivity conditions of the Higgs potential.
     * @return
     */
    double computeThValue();
    const LoopMediators * myLM;
};

#endif /* LOOPMEDIATORS_H */

