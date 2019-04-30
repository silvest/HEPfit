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

void Meson::ModelParameterMapInsert(std::map< std::string, std::reference_wrapper<const double> >& ModelParamMap)
{
    if (name.compare("P_0") == 0) {
        ModelParamMap.insert(std::make_pair("MP0", std::cref(mass)));
        ModelParamMap.insert(std::make_pair("tP0", std::cref(lifetime)));
        ModelParamMap.insert(std::make_pair("FP0", std::cref(decayconst)));
        return;
    }
    if (name.compare("P_P") == 0) {
        ModelParamMap.insert(std::make_pair("MPp", std::cref(mass)));
        ModelParamMap.insert(std::make_pair("tPp", std::cref(lifetime)));
        ModelParamMap.insert(std::make_pair("FPp", std::cref(decayconst)));
        return;
    }
    if (name.compare("K_0") == 0) {
        ModelParamMap.insert(std::make_pair("MK0", std::cref(mass)));
        ModelParamMap.insert(std::make_pair("tKl", std::cref(lifetime)));
        ModelParamMap.insert(std::make_pair("FK", std::cref(decayconst)));
        return;
    }
    if (name.compare("K_P") == 0) {
        ModelParamMap.insert(std::make_pair("MKp", std::cref(mass)));
        ModelParamMap.insert(std::make_pair("tKp", std::cref(lifetime)));
        ModelParamMap.insert(std::make_pair("FK", std::cref(decayconst)));
        ModelParamMap.insert(std::make_pair("alpha1kp", std::cref(gegenalpha[0])));
        ModelParamMap.insert(std::make_pair("alpha2kp", std::cref(gegenalpha[1])));
        return;
    }
    if (name.compare("D_0") == 0) {
        ModelParamMap.insert(std::make_pair("MD0", std::cref(mass)));
        ModelParamMap.insert(std::make_pair("tD0", std::cref(lifetime)));
        ModelParamMap.insert(std::make_pair("FD", std::cref(decayconst)));
        return;
    }
    if (name.compare("D_P") == 0) {
        ModelParamMap.insert(std::make_pair("MDP", std::cref(mass)));
        ModelParamMap.insert(std::make_pair("tDP", std::cref(lifetime)));
        ModelParamMap.insert(std::make_pair("FDP", std::cref(decayconst)));
        return;
    }
    if (name.compare("B_D") == 0) {
        ModelParamMap.insert(std::make_pair("MBd", std::cref(mass)));
        ModelParamMap.insert(std::make_pair("tBd", std::cref(lifetime)));
        ModelParamMap.insert(std::make_pair("FBd", std::cref(decayconst)));
        ModelParamMap.insert(std::make_pair("FBsoFBd", std::cref(FBsoFBd)));
        ModelParamMap.insert(std::make_pair("lambdaB", std::cref(lambdaM)));
        return;
    }
    if (name.compare("B_P") == 0) {
        ModelParamMap.insert(std::make_pair("MBp", std::cref(mass)));
        ModelParamMap.insert(std::make_pair("tBp", std::cref(lifetime)));
        ModelParamMap.insert(std::make_pair("FBp", std::cref(decayconst)));
        ModelParamMap.insert(std::make_pair("FBsoFBd", std::cref(FBsoFBd)));
        ModelParamMap.insert(std::make_pair("lambdaB", std::cref(lambdaM)));
        return;
    }
    if (name.compare("B_S") == 0) {
        ModelParamMap.insert(std::make_pair("MBs", std::cref(mass)));
        ModelParamMap.insert(std::make_pair("tBs", std::cref(lifetime)));
        ModelParamMap.insert(std::make_pair("FBs", std::cref(decayconst)));
        ModelParamMap.insert(std::make_pair("lambdaB", std::cref(lambdaM)));
        ModelParamMap.insert(std::make_pair("DGs_Gs", std::cref(Dgamma_gamma)));
        return;
    }
    if (name.compare("B_C") == 0) {
        ModelParamMap.insert(std::make_pair("MBc", std::cref(mass)));
        ModelParamMap.insert(std::make_pair("tBc", std::cref(lifetime)));
        ModelParamMap.insert(std::make_pair("FBc", std::cref(decayconst)));
        return;
    }
    if (name.compare("PHI") == 0) {
        ModelParamMap.insert(std::make_pair("Mphi", std::cref(mass)));
        ModelParamMap.insert(std::make_pair("tphi", std::cref(lifetime)));
        ModelParamMap.insert(std::make_pair("Fphi", std::cref(decayconst)));
        ModelParamMap.insert(std::make_pair("Fphip", std::cref(decayconst_p)));
        ModelParamMap.insert(std::make_pair("alpha2phi", std::cref(gegenalpha[1])));
        return;
    }
    if (name.compare("K_star") == 0) {
        ModelParamMap.insert(std::make_pair("MKstar", std::cref(mass)));
        ModelParamMap.insert(std::make_pair("tKstar", std::cref(lifetime)));
        ModelParamMap.insert(std::make_pair("FKstar", std::cref(decayconst)));
        ModelParamMap.insert(std::make_pair("FKstarp", std::cref(decayconst_p)));
        ModelParamMap.insert(std::make_pair("alpha1kst", std::cref(gegenalpha[0])));
        ModelParamMap.insert(std::make_pair("alpha2kst", std::cref(gegenalpha[1])));
        return;
    }
    if (name.compare("K_star_P") == 0) {
        ModelParamMap.insert(std::make_pair("MKstarP", std::cref(mass)));
        ModelParamMap.insert(std::make_pair("tKstarP", std::cref(lifetime)));
        ModelParamMap.insert(std::make_pair("FKstar", std::cref(decayconst)));
        ModelParamMap.insert(std::make_pair("FKstarp", std::cref(decayconst_p)));
        ModelParamMap.insert(std::make_pair("alpha1kst", std::cref(gegenalpha[0])));
        ModelParamMap.insert(std::make_pair("alpha2kst", std::cref(gegenalpha[1])));
        return;
    } 
    if (name.compare("D_star_P") == 0) {
        ModelParamMap.insert(std::make_pair("MDstarP", std::cref(mass)));
        ModelParamMap.insert(std::make_pair("tDstarP", std::cref(lifetime)));
        ModelParamMap.insert(std::make_pair("FDstarP", std::cref(decayconst)));
        return;
    } 
    if (name.compare("RHO") == 0) {
        ModelParamMap.insert(std::make_pair("Mrho", std::cref(mass)));
        ModelParamMap.insert(std::make_pair("trho", std::cref(lifetime)));
        ModelParamMap.insert(std::make_pair("Frho", std::cref(decayconst)));
        return;
    } 
    if (name.compare("RHO_P") == 0) {
        ModelParamMap.insert(std::make_pair("MrhoP", std::cref(mass)));
        ModelParamMap.insert(std::make_pair("trho", std::cref(lifetime)));
        ModelParamMap.insert(std::make_pair("Frho", std::cref(decayconst)));
        return;
    } 
    if (name.compare("OMEGA") == 0) {
        ModelParamMap.insert(std::make_pair("Momega", std::cref(mass)));
        ModelParamMap.insert(std::make_pair("tomega", std::cref(lifetime)));
        ModelParamMap.insert(std::make_pair("Fomega", std::cref(decayconst)));
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
    if (name_i.compare("D_P") == 0) return make_vector<std::string>() << "MDP"  << "tDP"  << "FDP";
    if (name_i.compare("B_D") == 0) return make_vector<std::string>() << "MBd" << "tBd" << "FBsoFBd" << "lambdaB";
    if (name_i.compare("B_P") == 0) return make_vector<std::string>() << "MBp" << "tBp" << "FBsoFBd" << "lambdaB";
    if (name_i.compare("B_S") == 0) return make_vector<std::string>() << "MBs" << "tBs" << "FBs" << "lambdaB" << "DGs_Gs";
    if (name_i.compare("B_C") == 0) return make_vector<std::string>() << "MBc" << "tBc" << "FBc";
    if (name_i.compare("PHI") == 0) return make_vector<std::string>() << "Mphi" << "tphi"  << "Fphi" << "Fphip" << "alpha2phi";
    if (name_i.compare("K_star") == 0) return make_vector<std::string>() << "MKstar"  << "tKstar"  << "FKstar" << "FKstarp" << "alpha1kst" << "alpha2kst";
    if (name_i.compare("K_star_P") == 0) return make_vector<std::string>() << "MKstarP" << "tKstar" << "FKstar" << "FKstarp" << "alpha1kst" << "alpha2kst";
    if (name_i.compare("D_star_P") == 0) return make_vector<std::string>() << "MDstarP"  << "tDstarP"  << "FDstarP";
    if (name_i.compare("RHO") == 0) return make_vector<std::string>() << "Mrho"  << "trho"  << "Frho";
    if (name_i.compare("RHO_P") == 0) return make_vector<std::string>() << "MrhoP"  << "trho"  << "Frho";
    if (name_i.compare("OMEGA") == 0) return make_vector<std::string>() << "Momega"  << "tomega"  << "Fomega";
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
    if (name.compare("D_P") == 0) {
        if (name_i.compare("MDP") == 0) {
            mass = value;
            return true;
        }
        if (name_i.compare("tDP") == 0) {
            lifetime = value;
            return true;
        }
        if (name_i.compare("FDP") == 0) {
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
    if (name.compare("B_C") == 0) {
        if (name_i.compare("MBc") == 0) {
            mass = value;
            return true;
        }
        if (name_i.compare("tBc") == 0) {
            lifetime = value;
            return true;
        }
        if (name_i.compare("FBc") == 0) {
            decayconst = value;
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
    if (name.compare("RHO") == 0) {
        if (name_i.compare("Mrho") == 0) {
            mass = value;
            return true;
        }
        if (name_i.compare("trho") == 0) {
            lifetime = value;
            return true;
        }
        if (name_i.compare("Frho") == 0) {
            decayconst = value;
            return true;
        }
    }
    if (name.compare("RHO_P") == 0) {
        if (name_i.compare("MrhoP") == 0) {
            mass = value;
            return true;
        }
        if (name_i.compare("trho") == 0) {
            lifetime = value;
            return true;
        }
        if (name_i.compare("Frho") == 0) {
            decayconst = value;
            return true;
        }
    }
    if (name.compare("OMEGA") == 0) {
        if (name_i.compare("Momega") == 0) {
            mass = value;
            return true;
        }
        if (name_i.compare("tomega") == 0) {
            lifetime = value;
            return true;
        }
        if (name_i.compare("Fomega") == 0) {
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