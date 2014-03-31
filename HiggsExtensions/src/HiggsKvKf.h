/* 
 * File:   HiggsKvKf.h
 * Author: silvest
 *
 * Created on March 27, 2014, 12:01 PM
 */

#ifndef HIGGSKVKF_H
#define	HIGGSKVKF_H
#include <StandardModel.h>

class HiggsKvKf : public StandardModel {
public:
    HiggsKvKf() : StandardModel() 
    {
    };
    HiggsKvKf(const HiggsKvKf& orig);
    virtual ~HiggsKvKf();
    
    virtual double computeKW()
    {
        return Kv;
    }
    virtual double computeKZ()
    {
        return Kv;
    }

    /**
     * @brief A method to compute the ratio of the @f[HZ\gamma@f] coupling in the current model and in the SM.
     * @return the ratio of the @f[HZ\gamma@f] coupling in the current model and in the SM
     */
    virtual double computeKZga() 
    {
        double gtt = computeGammaZgatt();
        double gWW = computeGammaZgaWW();
        double gtW = computeGammaZgatW();
        return (gtt*Kf*Kf + gWW*Kv*Kv + gtW*Kf*Kv)/(gtt + gWW + gtW);
    }

     /**
     * @brief A method to compute the ratio of the @f[H\gamma\gamma@f] coupling in the current model and in the SM.
     * @return the ratio of the @f[H\gamma\gamma@f] coupling in the current model and in the SM
     */
    virtual double computeKgaga() 
    {
        double gtt = computeGammagagatt();
        double gWW = computeGammagagaWW();
        double gtW = computeGammagagatW();
        return (gtt*Kf*Kf + gWW*Kv*Kv + gtW*Kf*Kv)/(gtt + gWW + gtW);
    }

   virtual double computeKb() 
    {
        return Kf;
    }

    virtual double computeKglgl() 
    {
        return Kf;
    }

    virtual double computeKt() 
    {
        std::cout << "warning: you asked me to return the ratio of the Higgs coupling in the Standard Model over itself.\n "
              <<  "This is probably not what you intended." << std::endl;
        return 1.;
    }

    virtual double computeKtau() 
    {
        std::cout << "warning: you asked me to return the ratio of the Higgs coupling in the Standard Model over itself.\n "
              <<  "This is probably not what you intended." << std::endl;
        return 1.;
    }


};

#endif	/* HIGGSKVKF_H */

