/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "Meson.h"
#include "std_make_vector.h"

Meson::Meson(double mass, double lifetime = 5.e29, double decayconst = 0., 
        double lambdaM = 0., double gegenalpha1 = 0., double gegenalpha2 = 0.)
{
    this->mass = mass;
    this->lifetime = lifetime;
    this->decayconst = decayconst;
    this->lambdaM = lambdaM;
    gegenalpha[0] = gegenalpha1;
    gegenalpha[1] = gegenalpha2;
}

Meson::Meson()
{
    mass = 0.;
    lifetime = 5.e29;
    decayconst = 0.;
    lambdaM = 0.;
    gegenalpha[0] = 0.;
    gegenalpha[1] = 0.;
    name = "";
}

Meson::~Meson()
{}

void Meson::ModelParameterMapInsert(std::map< std::string, boost::reference_wrapper<const double> >& ModelParamMap)
{
    if (name.compare("P_0") == 0) {
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("MP0", boost::cref(mass)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("tP0", boost::cref(lifetime)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("FP0", boost::cref(decayconst)));
        return;
    }
    if (name.compare("P_P") == 0) {
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("MPp", boost::cref(mass)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("tPp", boost::cref(lifetime)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("FPp", boost::cref(decayconst)));
        return;
    }
    if (name.compare("K_0") == 0) {
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("MK0", boost::cref(mass)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("tKl", boost::cref(lifetime)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("FK", boost::cref(decayconst)));
        return;
    }
    if (name.compare("K_P") == 0) {
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("MKp", boost::cref(mass)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("tKp", boost::cref(lifetime)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("FK", boost::cref(decayconst)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("alpha1kp", boost::cref(gegenalpha[0])));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("alpha2kp", boost::cref(gegenalpha[1])));
        return;
    }
    if (name.compare("D_0") == 0) {
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("MD0", boost::cref(mass)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("tD0", boost::cref(lifetime)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("FD", boost::cref(decayconst)));
        return;
    }
    if (name.compare("B_D") == 0) {
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("MBd", boost::cref(mass)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("tBd", boost::cref(lifetime)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("FBd", boost::cref(decayconst)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("FBsoFBd", boost::cref(FBsoFBd)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("lambdaB", boost::cref(lambdaM)));
        return;
    }
    if (name.compare("B_P") == 0) {
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("MBp", boost::cref(mass)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("tBp", boost::cref(lifetime)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("FBp", boost::cref(decayconst)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("FBsoFBd", boost::cref(FBsoFBd)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("lambdaB", boost::cref(lambdaM)));
        return;
    }
    if (name.compare("B_S") == 0) {
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("MBs", boost::cref(mass)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("tBs", boost::cref(lifetime)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("FBs", boost::cref(decayconst)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("lambdaB", boost::cref(lambdaM)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("DGs_Gs", boost::cref(Dgamma_gamma)));
        return;
    }
    if (name.compare("PHI") == 0) {
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Mphi", boost::cref(mass)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("tphi", boost::cref(lifetime)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Fphi", boost::cref(decayconst)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Fphip", boost::cref(decayconst_p)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("alpha2phi", boost::cref(gegenalpha[1])));
        return;
    }
    if (name.compare("K_star") == 0) {
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("MKstar", boost::cref(mass)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("tKstar", boost::cref(lifetime)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("FKstar", boost::cref(decayconst)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("FKstarp", boost::cref(decayconst_p)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("alpha1kst", boost::cref(gegenalpha[0])));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("alpha2kst", boost::cref(gegenalpha[1])));
        return;
    }
    if (name.compare("K_star_P") == 0) {
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("MKstarP", boost::cref(mass)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("tKstarP", boost::cref(lifetime)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("FKstar", boost::cref(decayconst)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("FKstarp", boost::cref(decayconst_p)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("alpha1kst", boost::cref(gegenalpha[0])));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("alpha2kst", boost::cref(gegenalpha[1])));
        return;
    } 
    if (name.compare("D_star_P") == 0) {
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("MDstarP", boost::cref(mass)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("tDstarP", boost::cref(lifetime)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("FDstarP", boost::cref(decayconst)));
        return;
    } else throw std::runtime_error(name + " is not implemented in Meson class");
}

std::vector<std::string> Meson::parameterList(std::string name_i)
{
    if (name_i.compare("P_0") == 0) return make_vector<std::string>() << "MP0" << "tP0" << "FP0";
    if (name_i.compare("P_P") == 0) return make_vector<std::string>() << "MPp" << "tPp" << "FPp";
    if (name_i.compare("K_0") == 0) return make_vector<std::string>() << "MK0" << "tKl" << "FK";
    if (name_i.compare("K_P") == 0) return make_vector<std::string>() << "MKp" << "tKp" << "FK" << "alpha1kp" << "alpha2kp";
    if (name_i.compare("D_0") == 0) return make_vector<std::string>() << "MD"  << "tD"  << "FD";
    if (name_i.compare("B_D") == 0) return make_vector<std::string>() << "MBd" << "tBd" << "FBsoFBd" << "lambdaB";
    if (name_i.compare("B_P") == 0) return make_vector<std::string>() << "MBp" << "tBp" << "FBsoFBd" << "lambdaB";
    if (name_i.compare("B_S") == 0) return make_vector<std::string>() << "MBs" << "tBs" << "FBs" << "lambdaB" << "DGs_Gs";
    if (name_i.compare("PHI") == 0) return make_vector<std::string>() << "Mphi" << "tphi"  << "Fphi" << "Fphip" << "alpha2phi";
    if (name_i.compare("K_star") == 0) return make_vector<std::string>() << "MKstar"  << "tKstar"  << "FKstar" << "FKstarp" << "alpha1kst" << "alpha2kst";
    if (name_i.compare("K_star_P") == 0) return make_vector<std::string>() << "MKstarP" << "tKstar" << "FKstar" << "FKstarp" << "alpha1kst" << "alpha2kst";
    if (name_i.compare("D_star_P") == 0) return make_vector<std::string>() << "MDstarP"  << "tDstarP"  << "FDstarP";
    else throw std::runtime_error(name_i + " is not implemented in Meson class");
}

bool Meson::setParameter(std::string name_i, double value) 
{
    if (name.compare("P_0") == 0) {
        if (name_i.compare("MP0") == 0) {
            mass = value;
            return true;
        }
        if (name_i.compare("tP0") == 0) {
            lifetime = value;
            return true;
        }
        if (name_i.compare("FP0") == 0) {
            decayconst = value;
            return true;
        }
    }
    if (name.compare("P_P") == 0) {
        if (name_i.compare("MPp") == 0) {
            mass = value;
            return true;
        }
        if (name_i.compare("tPp") == 0) {
            lifetime = value;
            return true;
        }
        if (name_i.compare("FPp") == 0) {
            decayconst = value;
            return true;
        }
    }
    if (name.compare("K_0") == 0) {
        if (name_i.compare("MK0") == 0) {
            mass = value;
            return true;
        }
        if (name_i.compare("tKl") == 0) {
            lifetime = value;
            return true;
        }
        if (name_i.compare("FK") == 0) {
            decayconst = value;
            return true;
        }
    }
    if (name.compare("K_P") == 0) {
        if (name_i.compare("MKp") == 0) {
            mass = value;
            return true;
        }
        if (name_i.compare("tKp") == 0) {
            lifetime = value;
            return true;
        }
        if (name_i.compare("FK") == 0) {
            decayconst = value;
            return true;
        }
        if (name_i.compare("alpha1kp") == 0) {
            gegenalpha[0] = value;
            return true;
        }
        if (name_i.compare("alpha2kp") == 0) {
            gegenalpha[1] = value;
            return true;
        }
    }
    if (name.compare("D_0") == 0) {
        if (name_i.compare("MD") == 0) {
            mass = value;
            return true;
        }
        if (name_i.compare("tD") == 0) {
            lifetime = value;
            return true;
        }
        if (name_i.compare("FD") == 0) {
            decayconst = value;
            return true;
        }
    }
    if (name.compare("B_D") == 0) {
        if (name_i.compare("MBd") == 0) {
            mass = value;
            return true;
        }
        if (name_i.compare("tBd") == 0) {
            lifetime = value;
            return true;
        }
        if (name_i.compare("FBsoFBd") == 0) {
            FBsoFBd = value;
            return true;
        }
        if (name_i.compare("lambdaB") == 0) {
            lambdaM = value;
            return true;
        }
    }
    if (name.compare("B_P") == 0) {
        if (name_i.compare("MBp") == 0) {
            mass = value;
            return true;
        }
        if (name_i.compare("tBp") == 0) {
            lifetime = value;
            return true;
        }
        if (name_i.compare("FBsoFBd") == 0) {
            FBsoFBd = value;
            return true;
        }
        if (name_i.compare("lambdaB") == 0) {
            lambdaM = value;
            return true;
        }
    }
    if (name.compare("B_S") == 0) {
        if (name_i.compare("MBs") == 0) {
            mass = value;
            return true;
        }
        if (name_i.compare("tBs") == 0) {
            lifetime = value;
            return true;
        }
        if (name_i.compare("FBs") == 0) {
            decayconst = value;
            return true;
        }
        if (name_i.compare("lambdaB") == 0) {
            lambdaM = value;
            return true;
        }
        if (name_i.compare("DGs_Gs") == 0) {
            Dgamma_gamma = value;
            return true;
        }
    }
    if (name.compare("PHI") == 0) {
        if (name_i.compare("Mphi") == 0) {
            mass = value;
            return true;
        }
        if (name_i.compare("tphi") == 0) {
            lifetime = value;
            return true;
        }
        if (name_i.compare("Fphi") == 0) {
            decayconst = value;
            return true;
        }
        if (name_i.compare("Fphip") == 0) {
            decayconst_p = value;
            return true;
        }
        if (name_i.compare("alpha2phi") == 0) {
            gegenalpha[1] = value;
            return true;
        }
    }
    if (name.compare("K_star") == 0) {
        if (name_i.compare("MKstar") == 0) {
            mass = value;
            return true;
        }
        if (name_i.compare("tKstar") == 0) {
            lifetime = value;
            return true;
        }
        if (name_i.compare("FKstar") == 0) {
            decayconst = value;
            return true;
        }
        if (name_i.compare("FKstarp") == 0) {
            decayconst_p = value;
            return true;
        }
        if (name_i.compare("alpha1kst") == 0) {
            gegenalpha[0] = value;
            return true;
        }
        if (name_i.compare("alpha2kst") == 0) {
            gegenalpha[1] = value;
            return true;
        }
    }
    if (name.compare("K_star_P") == 0) {
        if (name_i.compare("MKstarP") == 0) {
            mass = value;
            return true;
        }
        if (name_i.compare("tKstar") == 0) {
            lifetime = value;
            return true;
        }
        if (name_i.compare("FKstar") == 0) {
            decayconst = value;
            return true;
        }
        if (name_i.compare("FKstarp") == 0) {
            decayconst_p = value;
            return true;
        }
        if (name_i.compare("alpha1kst") == 0) {
            gegenalpha[0] = value;
            return true;
        }
        if (name_i.compare("alpha2kst") == 0) {
            gegenalpha[1] = value;
            return true;
        }
    }
    if (name.compare("D_star_P") == 0) {
        if (name_i.compare("MDstarP") == 0) {
            mass = value;
            return true;
        }
        if (name_i.compare("tDstarP") == 0) {
            lifetime = value;
            return true;
        }
        if (name_i.compare("FDstarP") == 0) {
            decayconst = value;
            return true;
        }
    }
    return false;
}

double Meson::computeWidth() const
{
    return (HCUT / lifetime);
}