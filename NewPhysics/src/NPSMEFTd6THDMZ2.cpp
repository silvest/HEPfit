/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/cppFiles/file.cc to edit this template
 */

#include "NPSMEFTd6THDMZ2.h"

const std::string NPSMEFTd6THDMZ2::NPSMEFTd6THDMZ2Vars[NNPSMEFTd6THDMZ2Vars] = {
    "logtb", "bma", "m12_2"
};



NPSMEFTd6THDMZ2::NPSMEFTd6THDMZ2()
: NPSMEFTd6ATHDM()
{

    setModelName("NPSMEFTd6THDMZ2");

    ModelParamMap.insert(std::make_pair("logtb", std::cref(logtb)));
    ModelParamMap.insert(std::make_pair("bma", std::cref(bma)));
    ModelParamMap.insert(std::make_pair("m12_2", std::cref(m12_2)));
    // Derived quantities, only to be used as observables
    ModelParamMap.insert(std::make_pair("cosbma", std::cref(cosbma)));
    ModelParamMap.insert(std::make_pair("Z1", std::cref(Z1)));
    ModelParamMap.insert(std::make_pair("lambda1_Z2", std::cref(lambda1_Z2)));
    ModelParamMap.insert(std::make_pair("lambda2_Z2", std::cref(lambda2_Z2)));
    ModelParamMap.insert(std::make_pair("lambda3_Z2", std::cref(lambda3_Z2)));
    ModelParamMap.insert(std::make_pair("lambda4_Z2", std::cref(lambda4_Z2)));
    ModelParamMap.insert(std::make_pair("lambda5_Z2", std::cref(lambda5_Z2)));

}



void NPSMEFTd6THDMZ2::setParameter(const std::string name, const double& value)
{

    if (name.compare("logtb") == 0) {

        logtb = value;
        tanb = pow(10., value); // Input is log10(tanb), but we work with tanb

    } else if (name.compare("bma") == 0) {

        bma = value;
        cosbma = cos(bma);

    } else if (name.compare("m12_2") == 0) {

        m12_2 = value; // Overwritten in setZ2Potential() if Decoupling = true

    } else if (name.compare("etau") == 0 || name.compare("etad") == 0 || name.compare("etae") == 0) {
        // Fixed by the Z2 symmetry in setZ2Yukawas(), any input value is ignored
    } else if (name.compare("Z2") == 0 || name.compare("Z3") == 0 || name.compare("Z4") == 0 ||
               name.compare("Z5") == 0 || name.compare("Z6") == 0 || name.compare("Z7") == 0) {
        // Computed in setZ2Potential(), any input value is ignored
    } else
        NPSMEFTd6ATHDM::setParameter(name, value);
}



void NPSMEFTd6THDMZ2::setZ2Yukawas()
{

    // Same conventions as in GeneralTHDMZ2: the eta's of NPSMEFTd6ATHDM play the role of
    // the alignment parameters sigma_u,d,l of the Higgs basis
    if (flag_model == "type1") {
        etau = 1./tanb;
        etad = 1./tanb;
        etae = 1./tanb;
    }
    else if (flag_model == "type2") {
        etau = 1./tanb;
        etad = -tanb;
        etae = -tanb;
    }
    else if (flag_model == "typeX") {
        etau = 1./tanb;
        etad = 1./tanb;
        etae = -tanb;
    }
    else if (flag_model == "typeY") {
        etau = 1./tanb;
        etad = -tanb;
        etae = 1./tanb;
    }
    else if (flag_model == "inert") {
        etau = 0.;
        etad = 0.;
        etae = 0.;
    }
    else
        throw std::runtime_error("NPSMEFTd6THDMZ2: the flag modelType has not been set (type1, type2, typeX, typeY or inert)");
}



void NPSMEFTd6THDMZ2::setZ2Potential()
{

    double v2 = v() * v();
    double mh2 = mHl * mHl;
    double M2 = sqY2 * sqY2; // mH^2 = mA^2 = mH+^2

    if (flag_model == "inert") {
        // The Z2 basis is the Higgs basis, exact alignment, and the potential does not enter
        // the tree-level matching
        bma = M_PI/2.;
        cosbma = 0.;
        m12_2 = 0.;
        Z1 = mh2 / v2;
        Z2 = 0.;
        Z3 = 0.;
        Z4 = 0.;
        Z5 = 0.;
        Z6 = 0.;
        Z7 = 0.;
        lambda1_Z2 = Z1;
        lambda2_Z2 = 0.;
        lambda3_Z2 = 0.;
        lambda4_Z2 = 0.;
        lambda5_Z2 = 0.;
        return;
    }

    double beta = atan(tanb);
    double cb = cos(beta);
    double sb = sin(beta);
    double s2b = sin(2. * beta);
    double c2b = cos(2. * beta);
    double cos2bmbma = cos(2. * (beta - bma));
    double sin2bmbma = sin(2. * (beta - bma));

    if (flag_decoupling) m12_2 = M2 * sb * cb;
    double M2aux = m12_2 / sb / cb;

    // Z2 basis, as in GeneralTHDMZ2 with mH^2 = mA^2 = mH+^2 = M2
    lambda1_Z2 = (M2 + mh2 + (M2 - mh2) * cos2bmbma - 2. * M2aux * sb * sb) / 2. / cb / cb / v2;
    lambda2_Z2 = (mh2 + M2 + (mh2 - M2) * cos2bmbma - 2. * M2aux * cb * cb) / 2. / sb / sb / v2;
    lambda3_Z2 = (2. * M2 - M2aux + (M2 - mh2) * sin2bmbma / s2b) / v2;
    lambda4_Z2 = (M2aux - M2) / v2;
    lambda5_Z2 = (M2aux - M2) / v2;

    // Rotation to the Higgs basis with lambda6 = lambda7 = 0
    double l345 = lambda3_Z2 + lambda4_Z2 + lambda5_Z2;
    double K = (lambda1_Z2 + lambda2_Z2 - 2. * l345) * s2b * s2b / 4.;
    Z1 = lambda1_Z2 * pow(cb, 4.) + lambda2_Z2 * pow(sb, 4.) + 0.5 * l345 * s2b * s2b;
    Z2 = lambda1_Z2 * pow(sb, 4.) + lambda2_Z2 * pow(cb, 4.) + 0.5 * l345 * s2b * s2b;
    Z3 = K + lambda3_Z2;
    Z4 = K + lambda4_Z2;
    Z5 = K + lambda5_Z2;
    Z6 = -0.5 * s2b * (lambda1_Z2 * cb * cb - lambda2_Z2 * sb * sb - l345 * c2b);
    Z7 = -0.5 * s2b * (lambda1_Z2 * sb * sb - lambda2_Z2 * cb * cb + l345 * c2b);
}



bool NPSMEFTd6THDMZ2::PostUpdate() {

    // The eta's and the Z's are set here and not in setParameter() since the parameters are
    // updated in alphabetical order
    setZ2Yukawas();
    setZ2Potential();

    if (!NPSMEFTd6ATHDM::PostUpdate()) return (false);

    return (true);
}



bool NPSMEFTd6THDMZ2::CheckParameters(const std::map<std::string, double>& DPars)
{

    // NPSMEFTd6General does not check its parameters, so we check here the ones that are needed by this model
    std::vector<std::string> mandatory = {"sqY2"};
    if (flag_model != "inert") {
        mandatory.push_back("logtb");
        mandatory.push_back("bma");
        if (!flag_decoupling) mandatory.push_back("m12_2");
    }
    for (const std::string& par : mandatory) {
        if (DPars.find(par) == DPars.end()) {
            std::cout << "ERROR: missing mandatory NPSMEFTd6THDMZ2 parameter " << par << std::endl;
            raiseMissingModelParameterCount();
            addMissingModelParameter(par);
        }
    }

    std::vector<std::string> ignored = {"etau", "etad", "etae", "Z2", "Z3", "Z4", "Z5", "Z6", "Z7"};
    if (flag_decoupling) ignored.push_back("m12_2");
    for (const std::string& par : ignored) {
        if (DPars.find(par) != DPars.end())
            std::cout << "WARNING: NPSMEFTd6THDMZ2 parameter " << par
                    << " is computed internally, the input value will be ignored" << std::endl;
    }

    return (NPSMEFTd6ATHDM::CheckParameters(DPars));
}



bool NPSMEFTd6THDMZ2::setFlag(const std::string name, const bool value)
{

    if (name.compare("Decoupling") == 0) {
        flag_decoupling = value;
        return true;
    }
    else
        return NPSMEFTd6ATHDM::setFlag(name, value);
}



bool NPSMEFTd6THDMZ2::setFlagStr(const std::string name, const std::string value)
{

    if (name.compare("modelType") == 0) {
        if (CheckModelType(value)) {
            flag_model = value;
            return true;
        }
        else
            throw std::runtime_error("NPSMEFTd6THDMZ2::setFlagStr(): invalid model-type flag "
                + name + "=" + value);
    }
    else
        return NPSMEFTd6ATHDM::setFlagStr(name, value);
}



bool NPSMEFTd6THDMZ2::CheckModelType(const std::string modeltype) const
{
    if (modeltype.compare("type1") == 0 ||
        modeltype.compare("type2") == 0 ||
        modeltype.compare("typeX") == 0 ||
        modeltype.compare("typeY") == 0 ||
        modeltype.compare("inert") == 0 )
        return true;
    else
        return false;
}
