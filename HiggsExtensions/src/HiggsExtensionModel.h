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
class HiggsExtensionModel : virtual public StandardModel {
public:
    HiggsExtensionModel() : StandardModel() {
    };
    HiggsExtensionModel(const HiggsExtensionModel& orig);
    virtual ~HiggsExtensionModel(){};
    virtual double computeKW() const =0;
    virtual double computeKZ() const =0;
    virtual double computeKZga() const =0;
    virtual double computeKgaga() const =0;
    virtual double computeKglgl() const =0;
    virtual double computeKb() const =0;
    virtual double computeKc() const =0; 
    virtual double computeKt() const =0;
    virtual double computeKtau() const =0;

    /**
     * @brief This method computes the ratio of the total Higgs width w.r.t SM.
     * @return The he ratio of the total Higgs width w.r.t SM
     */
    virtual double computeGTotalRatio() const =0;

    private:

};

#endif	/* HIGGSEXTENSIONMODEL_H */

