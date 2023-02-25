/* 
 * Copyright (C) 2023 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef ASL_H
#define ASL_H
#include "StandardModel.h"
#include "ThObservable.h"

#include "gslpp.h"


class Asl : public ThObservable{
public:
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] lep_i final leptons of the decay
     */
    Asl(const StandardModel& SM_i, QCD::lepton lep_i);
    
    /**
     * @brief Destructor.
     */
    virtual ~Asl();
    
    /**
     * @brief The update parameter method for Asl.
     */
    void updateParameters();
    
    double computeThValue ();
    
private:
    const StandardModel& mySM;/**< Model type */
    QCD::lepton lep;/**< Final leptons type */
    
    double Md;         
    double Mc;
    double Mb;
    double Mt;
    double MW;
    double MB;
    double rhob;
    double etab;
    
    //wrong basis for now: Misiak instead of Buras
    gslpp::vector<gslpp::complex> ** allcoeff;/**<vector that contains the Wilson coeffients */
    gslpp::complex C_1;/**<Wilson coeffients @f$C_1@f$*/
    gslpp::complex C_2;/**<Wilson coeffients @f$C_2@f$*/
    gslpp::complex C_3;/**<Wilson coeffients @f$C_3@f$*/
    gslpp::complex C_4;/**<Wilson coeffients @f$C_4@f$*/
    gslpp::complex C_5;/**<Wilson coeffients @f$C_5@f$*/
    gslpp::complex C_6;/**<Wilson coeffients @f$C_6@f$*/
    
    //combinations of Wilson coeffients arxiv.org/abs/hep-ph/0202010
    gslpp::complex K_1;
    gslpp::complex K_2;
    gslpp::complex K_1prime;
    gslpp::complex K_2prime;
    gslpp::complex K_3prime;
    
    double Mb2;
    double Mt2;
    double MB2;
    double MW2;
    double z;
    gslpp::complex K12;

    double B;    
    double B_s;
    double B_sprime;
    double eta_B;
    double eta_Bb;
    double S_0;
    
    gslpp::complex kappa;
  
    //CP asymmetry in semileptonic B decays in the SM
    gslpp::complex A_sl;
};
 
#endif /* ASL_H */

