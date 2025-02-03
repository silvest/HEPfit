/* 
 * Copyright (C) 2025 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef GENERALTHDMZ2_H
#define GENERALTHDMZ2_H

#include "GeneralTHDM.h"


class GeneralTHDMZ2: public GeneralTHDM {
public:

    static const int NGeneralTHDMZ2vars = 3;

    static const std::string GeneralTHDMZ2vars[NGeneralTHDMZ2vars];

    /**
     * @brief GeneralTHDMZ2 constructor
     */
    GeneralTHDMZ2();

    /**
     * @brief GeneralTHDMZ2 destructor
     */
    ~GeneralTHDMZ2();

    /**
     * @brief The post-update method for %GeneralTHDMZ2.
     * @details This method runs all the procedures that are need to be executed
     * after the model is successfully updated.
     * @return a boolean that is true if the execution is successful
     */
    virtual bool InitializeModel();

    /**
     * @brief The pre-update method for %GeneralTHDMZ2
     * @return a boolean that is true if the execution is successful
     */
    virtual bool PreUpdate();

    /**
     * @brief The update method for %GeneralTHDMZ2.
     * @details This method updates all the model parameters with given DPars.
     * @param[in] DPars a map of the parameters that are being updated in the Monte Carlo run
     * @return a boolean that is true if the execution is successful
     */
    virtual bool Update(const std::map<std::string, double>& DPars);

    /**
     * @brief The post-update method for %GeneralTHDMZ2.
     * @details This method runs all the procedures that are need to be executed
     * after the model is successfully updated.
     * @return a boolean that is true if the execution is successful
     */
    virtual bool PostUpdate();

    /**
     * @brief A method to check if all the mandatory parameters for %GeneralTHDMZ2
     * have been provided in model initialization.
     * @param[in] DPars a map of the parameters that are being updated in the Monte Carlo run
     * (including parameters that are varied and those that are held constant)
     * @return a boolean that is true if the execution is successful
     */
    virtual bool CheckParameters(const std::map<std::string, double>& DPars);

    /**
     * @brief A method to set a string flag of %GeneralTHDMZ2.
     * @param[in] name name of a model flag
     * @param[in] value the string to be assigned to the flag specified by name
     * @return a boolean that is true if the execution is successful
     */
    virtual bool setFlagStr(const std::string name, const std::string value);


protected:

    /**
     * @brief A method to set the value of a parameter of %GeneralTHDMZ2.
     * @param[in] name name of a model parameter
     * @param[in] value the value to be assigned to the parameter specified by name
     */
    virtual void setParameter(const std::string, const double&);

    /**
     * @brief A method to check if the model type name in string form is valid.
     * @param[in] modeltype GeneralTHDMZ2 model type name
     * @return a boolean that is true if the model type name is valid
     */
    bool CheckModelType(const std::string modeltype) const;


private:

    double tanb, bma, m12_2; ///< parameters exclusively in Z2 models: tan(beta), beta-alpha, m_12^2
    double beta, cosb, cos2b, cos4b, cos6b, sinb, sin2b, sin6b, cos2bma, sin2bma;
    double vev, mh_2, mH_2, mA_2, mHp_2;
    std::string flag_model;
};

#endif /* GENERALTHDMZ2_H */
