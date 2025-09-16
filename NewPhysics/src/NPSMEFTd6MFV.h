/* 
 * File:   NPSMEFTd6MFV.h
 * Author: Glioti
 *
 * Created on June 24, 2025
 */

#ifndef NPSMEFTD6MFV_H
#define NPSMEFTD6MFV_H

#include "NPSMEFTd6General.h"


class NPSMEFTd6MFV: public NPSMEFTd6General {
public:
    
    static const int NNPSMEFTd6MFVVars = 128+1;
    
    static std::string NPSMEFTd6MFVVars[NNPSMEFTd6MFVVars];
    
    NPSMEFTd6MFV();

    /**
     * @brief The post-update method for %NPSMEFTd6General.
     * @details This method runs all the procedures that are need to be executed
     * after the model is successfully updated.
     * @return a boolean that is true if the execution is successful
     */
    virtual bool PostUpdate();
	
	/**
	* @brief A method to initialize the model parameters.
	* @param[in] DPars a map of the parameters that are being updated in the Monte Carlo run
	* (including parameters that are varied and those that are held constant)
	* @return a boolean that is true if the execution is successful
	*/
	virtual bool Init(const std::map<std::string, double>& DPars);
    
    
protected:

double CG = 0.;
double CW = 0.;
double CHG = 0.;
double CHW = 0.;
double CHB = 0.;
double CHWB = 0.;
double CHD = 0.;
double CHbox = 0.;
double CH = 0.;
double CHl1 = 0.;
double CHl3 = 0.;
double CHe = 0.;
double Cll = 0.;
double Cll1 = 0.;
double Cee = 0.;
double Cle = 0.;
double CuH_0 = 0.;
double CuH_u = 0.;
double CuH_d = 0.;
double CuG_0 = 0.;
double CuG_u = 0.;
double CuG_d = 0.;
double CuW_0 = 0.;
double CuW_u = 0.;
double CuW_d = 0.;
double CuB_0 = 0.;
double CuB_u = 0.;
double CuB_d = 0.;
double CdH_0 = 0.;
double CdH_u = 0.;
double CdH_d = 0.;
double CdG_0 = 0.;
double CdG_u = 0.;
double CdG_d = 0.;
double CdW_0 = 0.;
double CdW_u = 0.;
double CdW_d = 0.;
double CdB_0 = 0.;
double CdB_u = 0.;
double CdB_d = 0.;
double CHq1_0 = 0.;
double CHq1_u = 0.;
double CHq1_d = 0.;
double CHq3_0 = 0.;
double CHq3_u = 0.;
double CHq3_d = 0.;
double CHu_0 = 0.;
double CHu_u = 0.;
double CHd_0 = 0.;
double CHd_d = 0.;
double CHud_ud = 0.;
double Clq1_0 = 0.;
double Clq1_u = 0.;
double Clq1_d = 0.;
double Clq3_0 = 0.;
double Clq3_u = 0.;
double Clq3_d = 0.;
double Cqe_0 = 0.;
double Cqe_u = 0.;
double Cqe_d = 0.;
double Clu_0 = 0.;
double Clu_u = 0.;
double Ceu_0 = 0.;
double Ceu_u = 0.;
double Cld_0 = 0.;
double Cld_d = 0.;
double Ced_0 = 0.;
double Ced_d = 0.;
double Cqq1_00 = 0.;
double Cqq1_0u = 0.;
double Cqq1_0d = 0.;
double Cqq1_u0 = 0.;
double Cqq1_uu = 0.;
double Cqq1_ud = 0.;
double Cqq1_d0 = 0.;
double Cqq1_du = 0.;
double Cqq1_dd = 0.;
double Cqq3_00 = 0.;
double Cqq3_0u = 0.;
double Cqq3_0d = 0.;
double Cqq3_u0 = 0.;
double Cqq3_uu = 0.;
double Cqq3_ud = 0.;
double Cqq3_d0 = 0.;
double Cqq3_du = 0.;
double Cqq3_dd = 0.;
double Cuu_00 = 0.;
double Cuu_0u = 0.;
double Cuu_u0 = 0.;
double Cuu_uu = 0.;
double Cdd_00 = 0.;
double Cdd_0d = 0.;
double Cdd_d0 = 0.;
double Cdd_dd = 0.;
double Cud1_00 = 0.;
double Cud1_u0 = 0.;
double Cud1_0d = 0.;
double Cud1_ud = 0.;
double Cud8_00 = 0.;
double Cud8_u0 = 0.;
double Cud8_0d = 0.;
double Cud8_ud = 0.;
double Cqu1_00 = 0.;
double Cqu1_u0 = 0.;
double Cqu1_d0 = 0.;
double Cqu1_0u = 0.;
double Cqu1_uu = 0.;
double Cqu1_du = 0.;
double Cqu8_00 = 0.;
double Cqu8_u0 = 0.;
double Cqu8_d0 = 0.;
double Cqu8_0u = 0.;
double Cqu8_uu = 0.;
double Cqu8_du = 0.;
double Cqd1_00 = 0.;
double Cqd1_u0 = 0.;
double Cqd1_d0 = 0.;
double Cqd1_0d = 0.;
double Cqd1_ud = 0.;
double Cqd1_dd = 0.;
double Cqd8_00 = 0.;
double Cqd8_u0 = 0.;
double Cqd8_d0 = 0.;
double Cqd8_0d = 0.;
double Cqd8_ud = 0.;
double Cqd8_dd = 0.;
double Cquqd1 = 0.;
double Cquqd8 = 0.;


virtual void setParameter(const std::string name, const double& value);

/**
* @brief An auxiliary method to set the WC of the general class
*/
void setNPSMEFTd6GeneralParameters();

private:

};

#endif /* NPSMEFTD6MFV_H */
