/* 
 * File:   HiggExtensionModel.h
 * Author: enrico
 *
 * Created on 2 aprile 2014, 14.30
 */

#ifndef HIGGSEXTENSIONMODEL_H
#define	HIGGSEXTENSIONMODEL_H
#include <StandardModel.h>

/**
 * @class HiggsExtensionModel
 * @ingroup HiggsExtensions
 * @brief An abstract model class extending the %StandardModel Higgs sector.
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details This is an abstract Model class associated
 * with an extension of the %StandardModel where Higgs couplings 
 * are rescaled. 
 * This class inherits from the %StandardModel class, and provides the SM values of
 * Higgs production cross sections and individual contributions to decay rates.
 */
class HiggsExtensionModel : public StandardModel {
public:
    HiggsExtensionModel() : StandardModel() {
    };
    HiggsExtensionModel(const HiggsExtensionModel& orig);
    virtual ~HiggsExtensionModel();
    virtual double computeKW()=0;
    virtual double computeKZ()=0;
    virtual double computeKZga()=0;
    virtual double computeKgaga()=0;
    virtual double computeKglgl()=0;
    virtual double computeKb()=0;
    virtual double computeKc()=0; 
    virtual double computeKt()=0;
    virtual double computeKtau()=0;
    /**
     * @brief This method computes the top loop contribution to @f[H\to\gamma\gamma@f] in the Standard Model.
     * Currently it returns the value of tab 40 in ref. @cite Heinemeyer:2013tqa
     * @return Width of H->gamma gamma (top loop contribution squared)
     */
    double computeGammagagatt()
    {
        return 662.84;
    }

    /**
     * @brief This method computes the @f[W@f] loop contribution to @f[H\to\gamma\gamma@f] in the Standard Model.
     * Currently it returns the value of tab 40 in ref. @cite Heinemeyer:2013tqa
     * @return Width of H->gamma gamma (W loop contribution squared)
     */
    double computeGammagagaWW()
    {
        return 14731.86;
    }

    /**
     * @brief This method computes the mixed @f[t-W@f] loop contribution to @f[H\to\gamma\gamma@f] in the Standard Model.
     * Currently it returns the value of tab 40 in ref. @cite Heinemeyer:2013tqa
     * @return Width of H->gamma gamma (top W loop interference)
     */
    double computeGammagagatW()
    {
        return -6249.93;
    }
    
    /**
     * @brief This method computes the top loop contribution to @f[H\to Z\gamma@f] in the Standard Model.
     * Currently it returns the value of tab 41 in ref. @cite Heinemeyer:2013tqa
     * @return Width of H->Z gamma (top loop contribution squared)
     */
    double computeGammaZgatt()
    {
        return 21.74;
    }

    /**
     * @brief This method computes the @f[W@f] loop contribution to @f[H\to Z\gamma@f] in the Standard Model.
     * Currently it returns the value of tab 41 in ref. @cite Heinemeyer:2013tqa
     * @return Width of H->Z gamma (W loop contribution squared)
     */
    double computeGammaZgaWW()
    {
        return 7005.6;
    }

    /**
     * @brief This method computes the mixed @f[t-W@f] loop contribution to @f[H\to Z\gamma@f] in the Standard Model.
     * Currently it returns the value of tab 41 in ref. @cite Heinemeyer:2013tqa
     * @return Width of H->Z gamma (top W loop interference)
     */
    double computeGammaZgatW()
    {
        return -780.4;
    }
    /**
     * @brief This method computes the W fusion contribution @f[\sigma_{WF}f] to higgs-production
     * cross section in the Standard Model.
     * Currently it returns the value of tab 37  in ref. @cite Heinemeyer:2013tqa
     * @return W fusion contribution @f[\sigma_{WF}f] to cross section.
     */
    double computeSigmaWF()
    {
        return 1.21;
    }
    /**
     * @brief This method computes the Z fusion contribution @f[\sigma_{ZF}f] to higgs-production
     * cross section in the Standard Model.
     * Currently it returns the value of tab 37  in ref. @cite Heinemeyer:2013tqa
     * @return W fusion contribution @f[\sigma_{ZF}f] to cross section.
     */
    double computeSigmaZF()
    {
        return 0.417;
    }
    /**
     * @brief This method computes the Z W interference fusion contribution @f[\sigma_{ZWF}f] to higgs-production
     * cross section in the Standard Model.
     * Negligible (0.1%) in the Standard  model.
     * @return Z W interference fusion contribution @f[\sigma_{ZWF}f] to cross section.
     */
    double computeSigmaZWF()
    {
        return 0.;
    }
    /**
     * @brief This method computes the @f[BR(H\to WW)@f] in the Standard Model.
     * Currently it returns the value of tables in appendix A (Mh=125.5 GeV)  in ref. @cite Heinemeyer:2013tqa
     * @return The @f[BR(H\to WW)@f] in the Standard Model
     */
    double computeBRWW()
    {
        return 2.23e-1;
    }
    /**
     * @brief This method computes the @f[BR(H\to ZZ)@f] in the Standard Model.
     * Currently it returns the value of tables in appendix A (Mh=125.5 GeV)  in ref. @cite Heinemeyer:2013tqa
     * @return The @f[BR(H\to ZZ)@f] in the Standard Model
     */
    double computeBRZZ()
    {
        return 2.76e-2;
    }
    /**
     * @brief This method computes the @f[BR(H\to\gamma\gamma)@f] in the Standard Model.
     * Currently it returns the value of tables in appendix A (Mh=125.5 GeV)  in ref. @cite Heinemeyer:2013tqa
     * @return The @f[BR(H\to\gamma\gamma)@f] in the Standard Model
     */
    double computeBRgaga()
    {
        return 2.28e-3;
    }
    /**
     * @brief This method computes the @f[BR(H\to Z\gamma)@f] in the Standard Model.
     * Currently it returns the value of tables in appendix A (Mh=125.5 GeV)  in ref. @cite Heinemeyer:2013tqa
     * @return The @f[BR(H\to Z\gamma)@f] in the Standard Model
     */
    double computeBRZga()
    {
        return 1.58e-3;
    }
    /**
     * @brief This method computes the @f[BR(H\to gg)@f] in the Standard Model.
     * Currently it returns the value of tables in appendix A (Mh=125.5 GeV)  in ref. @cite Heinemeyer:2013tqa
     * @return The @f[BR(H\to gg)@f] in the Standard Model
     */
    double computeBRglgl()
    {
        return 8.52e-2;
    }
    /**
     * @brief This method computes the @f[BR(H\to bb)@f] in the Standard Model.
     * Currently it returns the value of tables in appendix A (Mh=125.5 GeV)  in ref. @cite Heinemeyer:2013tqa
     * @return The @f[BR(H\to bb)@f] in the Standard Model
     */
    double computeBRbb()
    {
        return 5.69e-1;
    }
    /**
     * @brief This method computes the @f[BR(H\to \tau\tau)@f] in the Standard Model.
     * Currently it returns the value of tables in appendix A (Mh=125.5 GeV)  in ref. @cite Heinemeyer:2013tqa
     * @return The @f[BR(H\to \tau\tau)@f] in the Standard Model
     */
    double computeBRtautau()
    {
        return 6.24e-2;
    }
    /**
     * @brief This method computes the @f[BR(H\to cc)@f] in the Standard Model.
     * Currently it returns the value of tables in appendix A (Mh=125.5 GeV)  in ref. @cite Heinemeyer:2013tqa
     * @return The @f[BR(H\to cc)@f] in the Standard Model
     */
    double computeBRcc()
    {
        return 2.87e-2;
    }
    /**
     * @brief This method computes the ratio of the total Higgs width w.r.t SM.
     * @return The he ratio of the total Higgs width w.r.t SM
     */
    virtual double computeGTotalRatio()=0;

    private:

};

#endif	/* HIGGSEXTENSIONMODEL_H */

