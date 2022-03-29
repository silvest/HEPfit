/*
 * Copyright (C) 2017 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include <stdexcept>
#include "BParameter.h"
#include "std_make_vector.h"

BParameter::BParameter(int n)
: bpars(n, 0.)
{};

BParameter::BParameter(int n, std::string name_i)
: bpars(n, 0.)
{
    name = name_i;
    FlagCsi = true;
};

BParameter::~BParameter()
{}

void BParameter::ModelParameterMapInsert(std::map< std::string, std::reference_wrapper<const double> >& ModelParamMap)
{
    if (name.compare("BBs") == 0) {
        ModelParamMap.insert(std::make_pair("BBs1", std::cref(bpars(0))));
        ModelParamMap.insert(std::make_pair("BBs2", std::cref(bpars(1))));
        ModelParamMap.insert(std::make_pair("BBs3", std::cref(bpars(2))));
        ModelParamMap.insert(std::make_pair("BBs4", std::cref(bpars(3))));
        ModelParamMap.insert(std::make_pair("BBs5", std::cref(bpars(4))));
        ModelParamMap.insert(std::make_pair("BBsscale", std::cref(mu)));
        if (FlagCsi) {
            ModelParamMap.insert(std::make_pair("FBsSqrtBBs1", std::cref(FBsSqrtBBs1)));
            ModelParamMap.insert(std::make_pair("FBsSqrtBBs2", std::cref(FBsSqrtBBs2)));
            ModelParamMap.insert(std::make_pair("FBsSqrtBBs3", std::cref(FBsSqrtBBs3)));
            ModelParamMap.insert(std::make_pair("FBsSqrtBBs4", std::cref(FBsSqrtBBs4)));
            ModelParamMap.insert(std::make_pair("FBsSqrtBBs5", std::cref(FBsSqrtBBs5)));
        }
        return;
    }
    if (name.compare("BBd") == 0) {
        ModelParamMap.insert(std::make_pair("BBsoBBd", std::cref(BBsoBBd)));
        ModelParamMap.insert(std::make_pair("BBd2", std::cref(bpars(1))));
        ModelParamMap.insert(std::make_pair("BBd3", std::cref(bpars(2))));
        ModelParamMap.insert(std::make_pair("BBd4", std::cref(bpars(3))));
        ModelParamMap.insert(std::make_pair("BBd5", std::cref(bpars(4))));
        ModelParamMap.insert(std::make_pair("BBdscale", std::cref(mu)));
        if (FlagCsi) {
            ModelParamMap.insert(std::make_pair("csi", std::cref(csi)));
            ModelParamMap.insert(std::make_pair("FBdSqrtBBd2", std::cref(FBdSqrtBBd2)));
            ModelParamMap.insert(std::make_pair("FBdSqrtBBd3", std::cref(FBdSqrtBBd3)));
            ModelParamMap.insert(std::make_pair("FBdSqrtBBd4", std::cref(FBdSqrtBBd4)));
            ModelParamMap.insert(std::make_pair("FBdSqrtBBd5", std::cref(FBdSqrtBBd5)));
        }
        return;
    }
    if (name.compare("BK") == 0) {
        ModelParamMap.insert(std::make_pair("BK1", std::cref(bpars(0))));
        ModelParamMap.insert(std::make_pair("BK2", std::cref(bpars(1))));
        ModelParamMap.insert(std::make_pair("BK3", std::cref(bpars(2))));
        ModelParamMap.insert(std::make_pair("BK4", std::cref(bpars(3))));
        ModelParamMap.insert(std::make_pair("BK5", std::cref(bpars(4))));
        ModelParamMap.insert(std::make_pair("BKscale", std::cref(mu)));
        return;
    }
    if (name.compare("BD") == 0) {
        ModelParamMap.insert(std::make_pair("BD1", std::cref(bpars(0))));
        ModelParamMap.insert(std::make_pair("BD2", std::cref(bpars(1))));
        ModelParamMap.insert(std::make_pair("BD3", std::cref(bpars(2))));
        ModelParamMap.insert(std::make_pair("BD4", std::cref(bpars(3))));
        ModelParamMap.insert(std::make_pair("BD5", std::cref(bpars(4))));
        ModelParamMap.insert(std::make_pair("BDscale", std::cref(mu)));
        return;
    }
    if (name.compare("BKd1") == 0) {
        ModelParamMap.insert(std::make_pair("BK(1/2)1", std::cref(bpars(0))));
        ModelParamMap.insert(std::make_pair("BK(1/2)2", std::cref(bpars(1))));
        ModelParamMap.insert(std::make_pair("BK(1/2)3", std::cref(bpars(2))));
        ModelParamMap.insert(std::make_pair("BK(1/2)4", std::cref(bpars(3))));
        ModelParamMap.insert(std::make_pair("BK(1/2)5", std::cref(bpars(4))));
        ModelParamMap.insert(std::make_pair("BK(1/2)6", std::cref(bpars(5))));
        ModelParamMap.insert(std::make_pair("BK(1/2)7", std::cref(bpars(6))));
        ModelParamMap.insert(std::make_pair("BK(1/2)8", std::cref(bpars(7))));
        ModelParamMap.insert(std::make_pair("BK(1/2)9", std::cref(bpars(8))));
        ModelParamMap.insert(std::make_pair("BK(1/2)10", std::cref(bpars(9))));
        ModelParamMap.insert(std::make_pair("BKd1_scale", std::cref(mu)));
        return;
    }
    if (name.compare("BKd3") == 0) {
        ModelParamMap.insert(std::make_pair("BK(3/2)1", std::cref(bpars(0))));
        ModelParamMap.insert(std::make_pair("BK(3/2)2", std::cref(bpars(1))));
        ModelParamMap.insert(std::make_pair("BK(3/2)3", std::cref(bpars(2))));
        ModelParamMap.insert(std::make_pair("BK(3/2)4", std::cref(bpars(3))));
        ModelParamMap.insert(std::make_pair("BK(3/2)5", std::cref(bpars(4))));
        ModelParamMap.insert(std::make_pair("BK(3/2)6", std::cref(bpars(5))));
        ModelParamMap.insert(std::make_pair("BK(3/2)7", std::cref(bpars(6))));
        ModelParamMap.insert(std::make_pair("BK(3/2)8", std::cref(bpars(7))));
        ModelParamMap.insert(std::make_pair("BK(3/2)9", std::cref(bpars(8))));
        ModelParamMap.insert(std::make_pair("BK(3/2)10", std::cref(bpars(9))));
        ModelParamMap.insert(std::make_pair("BKd3_scale", std::cref(mu)));
        return;
    } else throw std::runtime_error(name + " is not implemented in BParameter class");
}

std::vector<std::string> BParameter::parameterList(std::string name_i)
{
    if (name_i.compare("BBs") == 0 && !FlagCsi) return make_vector<std::string>() << "BBs1" << "BBs2" << "BBs3" << "BBs4" << "BBs5" << "BBsscale" << "BBsscheme";
    if (name_i.compare("BBs") == 0 && FlagCsi) return make_vector<std::string>() << "FBsSqrtBBs1" << "FBsSqrtBBs2" << "FBsSqrtBBs3" << "FBsSqrtBBs4" << "FBsSqrtBBs5" << "BBsscale" << "BBsscheme";
    if (name_i.compare("BBd") == 0 && !FlagCsi) return make_vector<std::string>() << "BBsoBBd" << "BBd2" << "BBd3" << "BBd4" << "BBd5" << "BBdscale" << "BBdscheme";
    if (name_i.compare("BBd") == 0 && FlagCsi) return make_vector<std::string>() << "csi" << "FBdSqrtBBd2" << "FBdSqrtBBd3" << "FBdSqrtBBd4" << "FBdSqrtBBd5" << "BBdscale" << "BBdscheme";
    if (name_i.compare("BK") == 0) return make_vector<std::string>() << "BK1" << "BK2" << "BK3" << "BK4" << "BK5" << "BKscale" << "BKscheme";
    if (name_i.compare("BD") == 0) return make_vector<std::string>() << "BD1" << "BD2" << "BD3" << "BD4" << "BD5" << "BDscale" << "BDscheme";
    if (name_i.compare("BKd1") == 0) return make_vector<std::string>() << "BK(1/2)1" << "BK(1/2)2" << "BK(1/2)3" << "BK(1/2)4" << "BK(1/2)5"
        << "BK(1/2)6" << "BK(1/2)7" << "BK(1/2)8" << "BK(1/2)9" << "BK(1/2)10" << "BKd1_scale" << "BKd1_scheme";
    if (name_i.compare("BKd3") == 0) return make_vector<std::string>() << "BK(3/2)1" << "BK(3/2)2" << "BK(3/2)3" << "BK(3/2)4" << "BK(3/2)5"
        << "BK(3/2)6" << "BK(3/2)7" << "BK(3/2)8" << "BK(3/2)9" << "BK(3/2)10" << "BKd3_scale" << "BKd3_scheme";
    else throw std::runtime_error(name_i + " is not implemented in BParameter class");
}

bool BParameter::setParameter(std::string name_i, double value)
{

    if (name.compare("BBs") == 0 || name.compare("BBd") == 0) {
        if (!FlagCsi) {
            if (name_i.compare("BBsoBBd") == 0) {
                BBsoBBd = value;
                return true;
            } else if (name_i.compare("BBs1") == 0) {
                setBpars(0, value);
                return true;
            } else if (name_i.compare("BBs2") == 0) {
                setBpars(1, value);
                return true;
            } else if (name_i.compare("BBs3") == 0) {
                setBpars(2, value);
                return true;
            } else if (name_i.compare("BBs4") == 0) {
                setBpars(3, value);
                return true;
            } else if (name_i.compare("BBs5") == 0) {
                setBpars(4, value);
                return true;
            } else if (name_i.compare("BBd2") == 0) {
                setBpars(1, value);
                return true;
            } else if (name_i.compare("BBd3") == 0) {
                setBpars(2, value);
                return true;
            } else if (name_i.compare("BBd4") == 0) {
                setBpars(3, value);
                return true;
            } else if (name_i.compare("BBd5") == 0) {
                setBpars(4, value);
                return true;
            }
        } else {
            if (name_i.compare("csi") == 0) {
                csi = value;
                return true;
            } else if (name_i.compare("FBsSqrtBBs1") == 0) {
                FBsSqrtBBs1 = value;
                return true;
            } else if (name_i.compare("FBsSqrtBBs2") == 0) {
                FBsSqrtBBs2 = value;
                return true;
            } else if (name_i.compare("FBsSqrtBBs3") == 0) {
                FBsSqrtBBs3 = value;
                return true;
            } else if (name_i.compare("FBsSqrtBBs4") == 0) {
                FBsSqrtBBs4 = value;
                return true;
            } else if (name_i.compare("FBsSqrtBBs5") == 0) {
                FBsSqrtBBs5 = value;
                return true;
            } else if (name_i.compare("FBdSqrtBBd2") == 0) {
                FBdSqrtBBd2 = value;
                return true;
            } else if (name_i.compare("FBdSqrtBBd3") == 0) {
                FBdSqrtBBd3 = value;
                return true;
            } else if (name_i.compare("FBdSqrtBBd4") == 0) {
                FBdSqrtBBd4 = value;
                return true;
            } else if (name_i.compare("FBdSqrtBBd5") == 0) {
                FBdSqrtBBd5 = value;
                return true;
            }
        }
        if (name_i.compare("BBdscale") == 0) {
            setMu(value);
            return true;
        } else if (name_i.compare("BBdscheme") == 0) {
            setScheme((schemes) value);
            return true;
        } else if (name_i.compare("BBsscale") == 0) {
            setMu(value);
            return true;
        } else if (name_i.compare("BBsscheme") == 0) {
            setScheme((schemes) value);
            return true;
        }
    }

    if (name.compare("BK") == 0) {
        if (name_i.compare("BK1") == 0) {
            setBpars(0, value);
            return true;
        } else if (name_i.compare("BK2") == 0) {
            setBpars(1, value);
            return true;
        } else if (name_i.compare("BK3") == 0) {
            setBpars(2, value);
            return true;
        } else if (name_i.compare("BK4") == 0) {
            setBpars(3, value);
            return true;
        } else if (name_i.compare("BK5") == 0) {
            setBpars(4, value);
            return true;
        } else if (name_i.compare("BKscale") == 0) {
            setMu(value);
            return true;
        } else if (name_i.compare("BKscheme") == 0) {
            setScheme((schemes) value);
            return true;
        }
    }
    if (name.compare("BD") == 0) {
        if (name_i.compare("BD1") == 0) {
            setBpars(0, value);
            return true;
        } else if (name_i.compare("BD2") == 0) {
            setBpars(1, value);
            return true;
        } else if (name_i.compare("BD3") == 0) {
            setBpars(2, value);
            return true;
        } else if (name_i.compare("BD4") == 0) {
            setBpars(3, value);
            return true;
        } else if (name_i.compare("BD5") == 0) {
            setBpars(4, value);
            return true;
        } else if (name_i.compare("BDscale") == 0) {
            setMu(value);
            return true;
        } else if (name_i.compare("BDscheme") == 0) {
            setScheme((schemes) value);
            return true;
        }
    }

    if (name.compare("BKd1") == 0) {
        if (name_i.compare("BK(1/2)1") == 0) {
            setBpars(0, value);
            return true;
        } else if (name_i.compare("BK(1/2)2") == 0) {
            setBpars(1, value);
            return true;
        } else if (name_i.compare("BK(1/2)3") == 0) {
            setBpars(2, value);
            return true;
        } else if (name_i.compare("BK(1/2)4") == 0) {
            setBpars(3, value);
            return true;
        } else if (name_i.compare("BK(1/2)5") == 0) {
            setBpars(4, value);
            return true;
        } else if (name_i.compare("BK(1/2)6") == 0) {
            setBpars(5, value);
            return true;
        } else if (name_i.compare("BK(1/2)7") == 0) {
            setBpars(6, value);
            return true;
        } else if (name_i.compare("BK(1/2)8") == 0) {
            setBpars(7, value);
            return true;
        } else if (name_i.compare("BK(1/2)9") == 0) {
            setBpars(8, value);
            return true;
        } else if (name_i.compare("BK(1/2)10") == 0) {
            setBpars(9, value);
            return true;
        } else if (name_i.compare("BKd1_scale") == 0) {
            setMu(value);
            return true;
        } else if (name_i.compare("BKd1_scheme") == 0) {
            setScheme((schemes) value);
            return true;
        }
    }

    if (name.compare("BKd3") == 0) {
        if (name_i.compare("BK(3/2)1") == 0) {
            setBpars(0, value);
            return true;
        } else if (name_i.compare("BK(3/2)2") == 0) {
            setBpars(1, value);
            return true;
        } else if (name_i.compare("BK(3/2)3") == 0) {
            setBpars(2, value);
            return true;
        } else if (name_i.compare("BK(3/2)4") == 0) {
            setBpars(3, value);
            return true;
        } else if (name_i.compare("BK(3/2)5") == 0) {
            setBpars(4, value);
            return true;
        } else if (name_i.compare("BK(3/2)6") == 0) {
            setBpars(5, value);
            return true;
        } else if (name_i.compare("BK(3/2)7") == 0) {
            setBpars(6, value);
            return true;
        } else if (name_i.compare("BK(3/2)8") == 0) {
            setBpars(7, value);
            return true;
        } else if (name_i.compare("BK(3/2)9") == 0) {
            setBpars(8, value);
            return true;
        } else if (name_i.compare("BK(3/2)10") == 0) {
            setBpars(9, value);
            return true;
        } else if (name_i.compare("BKd3_scale") == 0) {
            setMu(value);
            return true;
        } else if (name_i.compare("BKd3_scheme") == 0) {
            setScheme((schemes) value);
            return true;
        }
    }

    return false;

}
