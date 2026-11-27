/* 
 * File:   NPSMEFTd6CH_RU.h
 * Author: Glioti
 *
 * Created on May 27, 2026
 */

#ifndef NPSMEFTD6CHRU_H
#define NPSMEFTD6CHRU_H

#include "gslpp.h"
#include "ThObservable.h"
#include "NPSMEFTd6MFV.h"



class NPSMEFTd6CHRU: public NPSMEFTd6MFV {
public:
    
	// 4 composite Higgs parameters + 76 O(1) coefficients + Lambda_NP
    static const int NNPSMEFTd6CHRUVars = 4 + 76 +1; 
    
    static std::string NPSMEFTd6CHRUVars[NNPSMEFTd6CHRUVars];
    
    NPSMEFTd6CHRU();

    /**
     * @brief The post-update method for %NPSMEFTd6MFV.
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

	double get_eps_uL() const {
		double yt = getRawSMEFTCoeffEW("YuR",2,2);
		return gstar * eps_u/yt - 1;
	}

	double get_eps_dL() const {
		double yb = getRawSMEFTCoeffEW("YdR",2,2);
		return gstar * eps_d/yb - 1;
	}
    
    
protected:

	///< Composite Higgs universal parameters.
	double mstar = 0., gstar = 0.;

	///< Composite Higgs partial compositeness mixing parameters.
	double eps_u = 0., eps_d = 0.;

	///< O(1) coefficients of bosonic SILH Operators
	double cH = 0., cB = 0., cW = 0., c2B = 0., c2W = 0., c2G = 0., cHB = 0., cHW = 0., cga = 0., cg = 0., c3W = 0., c3G = 0.;


	///< O(1) coefficients of dipole Operators
	double cuG = 0., cuW = 0., cuB = 0., cdG = 0., cdW = 0., cdB = 0.;
	double cuG_u = 0., cuW_u = 0., cuB_u = 0., cdG_u = 0., cdW_u = 0., cdB_u = 0.;
	double cuG_d = 0., cuW_d = 0., cuB_d = 0., cdG_d = 0., cdW_d = 0., cdB_d = 0.;


	///< O(1) coefficients of subleading EW couplings derivative Operators
	double cqD1_u = 0., cqD1_d = 0., cqD3_u = 0., cqD3_d = 0., cuD = 0., cdD = 0.;

	///< O(1) coefficients of Higgs current Operators
	double cHu = 0., cHd = 0., cHq1_u = 0., cHq1_d = 0., cHq3_u = 0., cHq3_d = 0.;

	///< O(1) coefficients of LL LL 4-fermion Operators
	double cqq1_uu = 0., cqq1_ud = 0., cqq1_dd = 0., cqq1p_uu = 0., cqq1p_ud = 0., cqq1p_dd = 0.;
	double cqq3_uu = 0., cqq3_ud = 0., cqq3_dd = 0., cqq3p_uu = 0., cqq3p_ud = 0., cqq3p_dd = 0.;


	///< O(1) coefficients of RR RR 4-fermion Operators
	double cuu = 0., cuup = 0., cdd = 0., cddp = 0., cud1 = 0., cud8 = 0.;


	///< O(1) coefficients of LL RR 4-fermion Operators
	double cqu1_u = 0., cqu1_d = 0., cqu1_y = 0., cqu8_u = 0., cqu8_d = 0., cqu8_y = 0.;
	double cqd1_u = 0., cqd1_d = 0., cqd1_y = 0., cqd8_u = 0., cqd8_d = 0., cqd8_y = 0.;


	///< O(1) coefficients of LR LR 4-fermion Operators
	double cquqd1 = 0., cquqd1p = 0., cquqd8 = 0., cquqd8p = 0.;

	///< Flag to enable the P_LR discrete symmetry
	bool FlagPLR = true;

	///< Flag to switch to Jeffreys prior for mstar
	bool FlagJeffreysPrior = false;

	virtual bool setFlag(const std::string, const bool);

	virtual void setParameter(const std::string name, const double& value);

	/**
	* @brief An auxiliary method to set the WC of the MFV class
	*/
	void setNPSMEFTd6MFVParameters();

};





class CHRU_eps_uL : public ThObservable {
public:   

    CHRU_eps_uL(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:        
    const NPSMEFTd6CHRU& myNPSMEFTd6CHRU;
};


class CHRU_eps_dL : public ThObservable {
public:   

    CHRU_eps_dL(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFTd6CHRU& myNPSMEFTd6CHRU;
};

class oneOvermstar2 : public ThObservable {
public:

	oneOvermstar2(const StandardModel& SM_i);
	
   
	double computeThValue();

private:
	const NPSMEFTd6CHRU& myNPSMEFTd6CHRU;
};



#endif /* NPSMEFTD6CHRU_H */