/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef SUSYMASSINSERTION_H
#define	SUSYMASSINSERTION_H

#include <gslpp.h>
#include <StandardModel.h>
#include "SUSYMassInsertionMatching.h" 

class SUSYMassInsertion: public StandardModel {
/**
 * @class SUSYMassIsertion
 * @brief This class inheritates Standard Model one and contains the
 * basic code elements to study SUSY effects to Standard Model processes in
 * mass insertion approximation.
 */
public:
    static const int NSusyMIvars = 3 + 144;
    static const std::string SusyMIvars[NSusyMIvars];

    /**
     * @brief SUSYMassInsertion constructor
     */
    SUSYMassInsertion();

    /**
     * @brief SUSYMassInsertion destructor
     */
    virtual ~SUSYMassInsertion();
    
    virtual bool InitializeModel();
    
    virtual SUSYMassInsertionMatching* getMyMatching() const
    {
        return mySUSYMIA;
    }
 
     ///////////////////////////////////////////////////////////////////////////
    // Parameters 
    
    /**
     * @brief a method to check the correct assignment forthe value of
     *  all the SM and SUSY parameters in respect to the one in Model.conf
     * @param a map containing the parameters (all as double) to be updated
     * @return a boolean true or false value
     */
    virtual bool Init(const std::map<std::string, double>& DPars);

    virtual bool PreUpdate();
    /**
     * @brief a method to update SM and SUSY parameters found in the argument
     * @param a map containing the parameters (all as double) to be updated
     */
    virtual bool Update(const std::map<std::string, double>& DPars);

    virtual bool PostUpdate();
    
    virtual bool CheckParameters(const std::map<std::string, double>& DPars);
    
    
    ///////////////////////////////////////////////////////////////////////////
    // functions for the input parameters of SUSY MIA model
    /**
     * 
     * @return the gluino mass
     */
    double getM3() const {
        return m3;
    }
    
    /**
     * @brief set the gluino mass
     * @param m3 a double for  the gluino mass
     */
    void setM3(double m3) {
        this->m3 = m3;
    }

    /**
     *
     * @return the mean value squark mass
     */
    double getMsq() const {
        return Msq;
    }

    /**
     * @brief set mean value squark mass
     * @param Msq a double for the mean value squark mass
     */
    void setMsq(double Msq) {
        this->Msq = Msq;
    }

    /**
     *
     * @return the \f$ \delta^{u}_LL \f$ mass insertion
     */
    gslpp::matrix<gslpp::complex> getDu_LL() const {
        return Du_LL;
    }

    /**
     * @brief set delta^u_LL mass insertion
     * @param delta^u_LL a matrix<complex> for the up-type left-left mass insertion parameters
     */
    void setDu_LL(gslpp::matrix<gslpp::complex> Du_LL) {
        this->Du_LL = Du_LL;
    }

    /**
     *
     * @return the \f$ \delta^{u}_LR \f$ mass insertion
     */
    gslpp::matrix<gslpp::complex> getDu_LR() const {
        return Du_LR;
    }

    /**
     * @brief set \f$ \delta^{u}_LR \f$ mass insertion
     * @param delta^u_LR a matrix<complex> for the up-type left-right mass insertion parameters
     */
    void setDu_LR(gslpp::matrix<gslpp::complex> Du_LL) {
        this->Du_LR = Du_LR;
    }

    /**
     *
     * @return the \f$ \delta^{u}_RL \f$ mass insertion
     */
    gslpp::matrix<gslpp::complex> getDu_RL() const {
        return Du_RL;
    }

    /**
     * @brief set \f$ \delta^{u}_RL \f$ mass insertion
     * @param delta^u_RL a matrix<complex> for the up-type right-left mass insertion parameters
     */
    void setDu_RL(gslpp::matrix<gslpp::complex> Du_LL) {
        this->Du_RL = Du_RL;
    }

    /**
     *
     * @return the \f$ \delta^u_RR \f$ mass insertion
     */
    gslpp::matrix<gslpp::complex> getDu_RR() const {
        return Du_RR;
    }

    /**
     * @brief set \f$ \delta^{u}_RR \f$ mass insertion
     * @param delta^u_RR a matrix<complex> for the up-type right-right mass insertion parameters
     */
    void setDu_RR(gslpp::matrix<gslpp::complex> Du_LL) {
        this->Du_RR = Du_RR;
    }

    /**
     *
     * @return the \f$ \delta^{d}_LL \f$ mass insertion
     */
    gslpp::matrix<gslpp::complex> getDd_LL() const {
        return Dd_LL;
    }

    /**
     * @brief set \f$ \delta^{d}_LL \f$ mass insertion
     * @param delta^d_LL a matrix<complex> for the down-type left-left mass insertion parameters
     */
    void setDd_LL(gslpp::matrix<gslpp::complex> Dd_LL) {
        this->Dd_LL = Dd_LL;
    }

    /**
     *
     * @return the \f$ \delta^{d}_LR \f$ mass insertion
     */
    gslpp::matrix<gslpp::complex> getDd_LR() const {
        return Dd_LR;
    }

    /**
     * @brief set \f$ \delta^{d}_LR \f$ mass insertion
     * @param delta^d_LR a matrix<complex> for the down-type left-right mass insertion parameters
     */
    void setDd_LR(gslpp::matrix<gslpp::complex> Dd_LL) {
        this->Dd_LR = Dd_LR;
    }

    /**
     *
     * @return the \f$ \delta^{d}_RL \f$ mass insertion
     */
    gslpp::matrix<gslpp::complex> getDd_RL() const {
        return Dd_RL;
    }

    /**
     * @brief set \f$ \delta^{d}_RL \f$ mass insertion
     * @param delta^d_RL a matrix<complex> for the down-type right-left mass insertion parameters
     */
    void setDd_RL(gslpp::matrix<gslpp::complex> Dd_LL) {
        this->Dd_RL = Dd_RL;
    }

    /**
     *
     * @return the \f$ \delta^{d}_RR \f$ mass insertion
     */
    gslpp::matrix<gslpp::complex> getDd_RR() const {
        return Dd_RR;
    }

    /**
     * @brief set \f$ \delta^{d}_RR \f$ mass insertion
     * @param delta^d_RR a matrix<complex> for the down-type right-right mass insertion parameters
     */
    void setDd_RR(gslpp::matrix<gslpp::complex> Dd_LL) {
        this->Dd_RR = Dd_RR;
    }

    /** 
     * @return get the SUSY matching scale
     */
    double getMuM() const{
        return MuM;
    }

    /**
     * @brief set the SUSY matching scale
     * @param MuM a double for the SUSY matching scale
     */
    void setMuM(double MuM){
        this->MuM = MuM;
    }
    
protected:

    /**
     * @brief a method to set the value of all the SUSY parameters given as input in Model.conf
     */
    virtual void setParameter(const std::string, const double&);

    gslpp::matrix<gslpp::complex> Du_LL, Du_LR, Du_RL, Du_RR;
    gslpp::matrix<gslpp::complex> Dd_LL, Dd_LR, Dd_RL, Dd_RR;
    double Msq, m3, MuM;
    
    double rDULL11, rDULL12, rDULL13, rDULL21, rDULL22, rDULL23, rDULL31, rDULL32, rDULL33;
    double iDULL11, iDULL12, iDULL13, iDULL21, iDULL22, iDULL23, iDULL31, iDULL32, iDULL33;
    
    double rDURL11, rDURL12, rDURL13, rDURL21, rDURL22, rDURL23, rDURL31, rDURL32, rDURL33;
    double iDURL11, iDURL12, iDURL13, iDURL21, iDURL22, iDURL23, iDURL31, iDURL32, iDURL33;
    
    double rDULR11, rDULR12, rDULR13, rDULR21, rDULR22, rDULR23, rDULR31, rDULR32, rDULR33;
    double iDULR11, iDULR12, iDULR13, iDULR21, iDULR22, iDULR23, iDULR31, iDULR32, iDULR33;
    
    double rDURR11, rDURR12, rDURR13, rDURR21, rDURR22, rDURR23, rDURR31, rDURR32, rDURR33;
    double iDURR11, iDURR12, iDURR13, iDURR21, iDURR22, iDURR23, iDURR31, iDURR32, iDURR33;
    
    double rDDLL11, rDDLL12, rDDLL13, rDDLL21, rDDLL22, rDDLL23, rDDLL31, rDDLL32, rDDLL33;
    double iDDLL11, iDDLL12, iDDLL13, iDDLL21, iDDLL22, iDDLL23, iDDLL31, iDDLL32, iDDLL33;
    
    double rDDRL11, rDDRL12, rDDRL13, rDDRL21, rDDRL22, rDDRL23, rDDRL31, rDDRL32, rDDRL33;
    double iDDRL11, iDDRL12, iDDRL13, iDDRL21, iDDRL22, iDDRL23, iDDRL31, iDDRL32, iDDRL33;
    
    double rDDLR11, rDDLR12, rDDLR13, rDDLR21, rDDLR22, rDDLR23, rDDLR31, rDDLR32, rDDLR33;
    double iDDLR11, iDDLR12, iDDLR13, iDDLR21, iDDLR22, iDDLR23, iDDLR31, iDDLR32, iDDLR33;
    
    double rDDRR11, rDDRR12, rDDRR13, rDDRR21, rDDRR22, rDDRR23, rDDRR31, rDDRR32, rDDRR33;
    double iDDRR11, iDDRR12, iDDRR13, iDDRR21, iDDRR22, iDDRR23, iDDRR31, iDDRR32, iDDRR33;

private:    
    SUSYMassInsertionMatching* mySUSYMIA;

};        

#endif	/* SUSYMASSINSERTION_H */



