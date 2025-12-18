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
    
    static const int NNPSMEFTd6MFVVars = 178+1;
    
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

	///< Coefficients of the MFV expansion of the dimension-6 operator coefficient \f$(C_{Hl}^{(1)})_{ij}\f$.
	double CHl1_LNP = 0.;

	///< Coefficients of the MFV expansion of the dimension-6 operator coefficient \f$(C_{Hl}^{(3)})_{ij}\f$.
	double CHl3_LNP = 0.;

	///< Coefficients of the MFV expansion of the dimension-6 operator coefficient \f$(C_{He})_{ij}\f$.
	double CHe_LNP = 0.;

	///< Coefficients of the MFV expansion of the dimension-6 operator coefficient \f$(C_{ll})_{ijkm}\f$.
	double Cll_aabb_LNP = 0.;

	///< Coefficients of the MFV expansion of the dimension-6 operator coefficient \f$(C_{ll})_{ijkm}\f$.
	double Cll_abba_LNP = 0.;

	///< Coefficients of the MFV expansion of the dimension-6 operator coefficient \f$(C_{ee})_{ijkm}\f$.
	double Cee_LNP = 0.;

	///< Coefficients of the MFV expansion of the dimension-6 operator coefficient \f$(C_{le})_{ijkm}\f$.
	double Cle_LNP = 0.;

	///< Coefficients of the MFV expansion of the dimension-6 operator coefficient \f$(C_{uH})_{ij}\f$.
	double CuH_0_LNP = 0., CuH_u_LNP = 0., CuH_d_LNP = 0.;

	///< Coefficients of the MFV expansion of the dimension-6 operator coefficient \f$(C_{uG})_{ij}\f$.
	double CuG_0_LNP = 0., CuG_u_LNP = 0., CuG_d_LNP = 0.;

	///< Coefficients of the MFV expansion of the dimension-6 operator coefficient \f$(C_{uW})_{ij}\f$.
	double CuW_0_LNP = 0., CuW_u_LNP = 0., CuW_d_LNP = 0.;

	///< Coefficients of the MFV expansion of the dimension-6 operator coefficient \f$(C_{uB})_{ij}\f$.
	double CuB_0_LNP = 0., CuB_u_LNP = 0., CuB_d_LNP = 0.;

	///< Coefficients of the MFV expansion of the dimension-6 operator coefficient \f$(C_{dH})_{ij}\f$.
	double CdH_0_LNP = 0., CdH_u_LNP = 0., CdH_d_LNP = 0.;

	///< Coefficients of the MFV expansion of the dimension-6 operator coefficient \f$(C_{dG})_{ij}\f$.
	double CdG_0_LNP = 0., CdG_u_LNP = 0., CdG_d_LNP = 0.;

	///< Coefficients of the MFV expansion of the dimension-6 operator coefficient \f$(C_{dW})_{ij}\f$.
	double CdW_0_LNP = 0., CdW_u_LNP = 0., CdW_d_LNP = 0.;

	///< Coefficients of the MFV expansion of the dimension-6 operator coefficient \f$(C_{dB})_{ij}\f$.
	double CdB_0_LNP = 0., CdB_u_LNP = 0., CdB_d_LNP = 0.;

	///< Coefficients of the MFV expansion of the dimension-6 operator coefficient \f$(C_{Hq}^{(1)})_{ij}\f$.
	double CHq1_0_LNP = 0., CHq1_u_LNP = 0., CHq1_d_LNP = 0.;

	///< Coefficients of the MFV expansion of the dimension-6 operator coefficient \f$(C_{Hq}^{(3)})_{ij}\f$.
	double CHq3_0_LNP = 0., CHq3_u_LNP = 0., CHq3_d_LNP = 0.;

	///< Coefficients of the MFV expansion of the dimension-6 operator coefficient \f$(C_{Hu})_{ij}\f$.
	double CHu_0_LNP = 0., CHu_u_LNP = 0.;

	///< Coefficients of the MFV expansion of the dimension-6 operator coefficient \f$(C_{Hd})_{ij}\f$.
	double CHd_0_LNP = 0., CHd_d_LNP = 0.;

	///< Coefficients of the MFV expansion of the dimension-6 operator coefficient \f$(C_{Hud})_{ij}\f$.
	double CHud_ud_LNP = 0.;

	///< Coefficients of the MFV expansion of the dimension-6 operator coefficient \f$(C_{lq}^{(1)})_{ijkm}\f$.
	double Clq1_0_LNP = 0., Clq1_u_LNP = 0., Clq1_d_LNP = 0.;

	///< Coefficients of the MFV expansion of the dimension-6 operator coefficient \f$(C_{lq}^{(3)})_{ijkm}\f$.
	double Clq3_0_LNP = 0., Clq3_u_LNP = 0., Clq3_d_LNP = 0.;

	///< Coefficients of the MFV expansion of the dimension-6 operator coefficient \f$(C_{qe})_{ijkm}\f$.
	double Cqe_0_LNP = 0., Cqe_u_LNP = 0., Cqe_d_LNP = 0.;

	///< Coefficients of the MFV expansion of the dimension-6 operator coefficient \f$(C_{lu})_{ijkm}\f$.
	double Clu_0_LNP = 0., Clu_u_LNP = 0.;

	///< Coefficients of the MFV expansion of the dimension-6 operator coefficient \f$(C_{eu})_{ijkm}\f$.
	double Ceu_0_LNP = 0., Ceu_u_LNP = 0.;

	///< Coefficients of the MFV expansion of the dimension-6 operator coefficient \f$(C_{ld})_{ijkm}\f$.
	double Cld_0_LNP = 0., Cld_d_LNP = 0.;

	///< Coefficients of the MFV expansion of the dimension-6 operator coefficient \f$(C_{ed})_{ijkm}\f$.
	double Ced_0_LNP = 0., Ced_d_LNP = 0.;

	///< Coefficients of the MFV expansion of the dimension-6 operator coefficient \f$(C_{qq}^{(1)})_{ijkm}\f$.
	double Cqq1_00_LNP = 0., Cqq1_u0_LNP = 0., Cqq1_d0_LNP = 0., Cqq1_uu_LNP = 0., Cqq1_dd_LNP = 0., Cqq1_ud_LNP = 0., Cqq1p_00_LNP = 0., Cqq1p_u0_LNP = 0., Cqq1p_d0_LNP = 0., Cqq1p_uu_LNP = 0., Cqq1p_dd_LNP = 0., Cqq1p_ud_LNP = 0.;

	///< Coefficients of the MFV expansion of the dimension-6 operator coefficient \f$(C_{qq}^{(3)})_{ijkm}\f$.
	double Cqq3_00_LNP = 0., Cqq3_u0_LNP = 0., Cqq3_d0_LNP = 0., Cqq3_uu_LNP = 0., Cqq3_dd_LNP = 0., Cqq3_ud_LNP = 0., Cqq3p_00_LNP = 0., Cqq3p_u0_LNP = 0., Cqq3p_d0_LNP = 0., Cqq3p_uu_LNP = 0., Cqq3p_dd_LNP = 0., Cqq3p_ud_LNP = 0.;

	///< Coefficients of the MFV expansion of the dimension-6 operator coefficient \f$(C_{uu})_{ijkm}\f$.
	double Cuu_00_LNP = 0., Cuu_u0_LNP = 0., Cuu_uu_LNP = 0., Cuup_00_LNP = 0., Cuup_u0_LNP = 0., Cuup_uu_LNP = 0.;

	///< Coefficients of the MFV expansion of the dimension-6 operator coefficient \f$(C_{dd})_{ijkm}\f$.
	double Cdd_00_LNP = 0., Cdd_d0_LNP = 0., Cdd_dd_LNP = 0., Cddp_00_LNP = 0., Cddp_d0_LNP = 0., Cddp_dd_LNP = 0.;

	///< Coefficients of the MFV expansion of the dimension-6 operator coefficient \f$(C_{ud}^{(1)})_{ijkm}\f$.
	double Cud1_00_LNP = 0., Cud1_u0_LNP = 0., Cud1_0d_LNP = 0., Cud1_ud_LNP = 0., Cud1p_ud_LNP = 0.;

	///< Coefficients of the MFV expansion of the dimension-6 operator coefficient \f$(C_{ud}^{(8)})_{ijkm}\f$.
	double Cud8_00_LNP = 0., Cud8_u0_LNP = 0., Cud8_0d_LNP = 0., Cud8_ud_LNP = 0., Cud8p_ud_LNP = 0.;

	///< Coefficients of the MFV expansion of the dimension-6 operator coefficient \f$(C_{qu}^{(1)})_{ijkm}\f$.
	double Cqu1_00_LNP = 0., Cqu1_u0_LNP = 0., Cqu1_d0_LNP = 0., Cqu1_0u_LNP = 0., Cqu1_uu_LNP = 0., Cqu1_du_LNP = 0., Cqu1_y_LNP = 0., Cqu1_uy_LNP = 0., Cqu1_dy_LNP = 0., Cqu1_yu_LNP = 0., Cqu1_yd_LNP = 0.;

	///< Coefficients of the MFV expansion of the dimension-6 operator coefficient \f$(C_{qu}^{(8)})_{ijkm}\f$.
	double Cqu8_00_LNP = 0., Cqu8_u0_LNP = 0., Cqu8_d0_LNP = 0., Cqu8_0u_LNP = 0., Cqu8_uu_LNP = 0., Cqu8_du_LNP = 0., Cqu8_y_LNP = 0., Cqu8_uy_LNP = 0., Cqu8_dy_LNP = 0., Cqu8_yu_LNP = 0., Cqu8_yd_LNP = 0.;

	///< Coefficients of the MFV expansion of the dimension-6 operator coefficient \f$(C_{qd}^{(1)})_{ijkm}\f$.
	double Cqd1_00_LNP = 0., Cqd1_u0_LNP = 0., Cqd1_d0_LNP = 0., Cqd1_0d_LNP = 0., Cqd1_ud_LNP = 0., Cqd1_dd_LNP = 0., Cqd1_y_LNP = 0., Cqd1_uy_LNP = 0., Cqd1_dy_LNP = 0., Cqd1_yu_LNP = 0., Cqd1_yd_LNP = 0.;

	///< Coefficients of the MFV expansion of the dimension-6 operator coefficient \f$(C_{qd}^{(8)})_{ijkm}\f$.
	double Cqd8_00_LNP = 0., Cqd8_u0_LNP = 0., Cqd8_d0_LNP = 0., Cqd8_0d_LNP = 0., Cqd8_ud_LNP = 0., Cqd8_dd_LNP = 0., Cqd8_y_LNP = 0., Cqd8_uy_LNP = 0., Cqd8_dy_LNP = 0., Cqd8_yu_LNP = 0., Cqd8_yd_LNP = 0.;

	///< Coefficients of the MFV expansion of the dimension-6 operator coefficient \f$(C_{quqd}^{(1)})_{ijkm}\f$.
	double Cquqd1_00_LNP = 0., Cquqd1_u0_LNP = 0., Cquqd1_d0_LNP = 0., Cquqd1_0u_LNP = 0., Cquqd1_0d_LNP = 0., Cquqd1p_00_LNP = 0., Cquqd1p_u0_LNP = 0., Cquqd1p_d0_LNP = 0., Cquqd1p_0u_LNP = 0., Cquqd1p_0d_LNP = 0.;

	///< Coefficients of the MFV expansion of the dimension-6 operator coefficient \f$(C_{quqd}^{(8)})_{ijkm}\f$.
	double Cquqd8_00_LNP = 0., Cquqd8_u0_LNP = 0., Cquqd8_d0_LNP = 0., Cquqd8_0u_LNP = 0., Cquqd8_0d_LNP = 0., Cquqd8p_00_LNP = 0., Cquqd8p_u0_LNP = 0., Cquqd8p_d0_LNP = 0., Cquqd8p_0u_LNP = 0., Cquqd8p_0d_LNP = 0.;



	virtual void setParameter(const std::string name, const double& value);

	/**
	* @brief An auxiliary method to set the WC of the general class
	*/
	void setNPSMEFTd6GeneralParameters();

private:

};

#endif /* NPSMEFTD6MFV_H */
