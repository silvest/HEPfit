/*
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "ThObsFactory.h"
#include "ThObsFactoryConstants.h"
#include "HiggsThObservables.h"
#include "OptimizedObservablesSMEFTd6.h"
#include "EWObservables.h"
using namespace ThObsConst;

void ThObsFactory::registerHiggsLeptonObservables()
{
    //-----  Full Signal strengths per prod and decay: e+ e- colliders  ----------
    //
    // Pure WBF
    obsThFactory["mueeWBFbb230"] = [=](const StandardModel& SM) { return new mueeWBFbb(SM, sqrt_s_leptcoll_230, 0., 0.); };
    obsThFactory["mueeWBFbb240"] = [=](const StandardModel& SM) { return new mueeWBFbb(SM, sqrt_s_leptcoll_240, 0., 0.); };
    obsThFactory["mueeWBFbb250"] = [=](const StandardModel& SM) { return new mueeWBFbb(SM, sqrt_s_leptcoll_250, 0., 0.); };
    obsThFactory["mueeWBFbb350"] = [=](const StandardModel& SM) { return new mueeWBFbb(SM, sqrt_s_leptcoll_350, 0., 0.); };
    obsThFactory["mueeWBFbb365"] = [=](const StandardModel& SM) { return new mueeWBFbb(SM, sqrt_s_leptcoll_365, 0., 0.); };
    obsThFactory["mueeWBFbb380"] = [=](const StandardModel& SM) { return new mueeWBFbb(SM, sqrt_s_leptcoll_380, 0., 0.); };
    obsThFactory["mueeWBFbb500"] = [=](const StandardModel& SM) { return new mueeWBFbb(SM, sqrt_s_leptcoll_500, 0., 0.); };
    obsThFactory["mueeWBFbb550"] = [=](const StandardModel& SM) { return new mueeWBFbb(SM, sqrt_s_leptcoll_550, 0., 0.); };
    obsThFactory["mueeWBFbb1000"] = [=](const StandardModel& SM) { return new mueeWBFbb(SM, sqrt_s_leptcoll_1000, 0., 0.); };
    obsThFactory["mueeWBFbb1400"] = [=](const StandardModel& SM) { return new mueeWBFbb(SM, sqrt_s_leptcoll_1400, 0., 0.); };
    obsThFactory["mueeWBFbb1500"] = [=](const StandardModel& SM) { return new mueeWBFbb(SM, sqrt_s_leptcoll_1500, 0., 0.); };
    obsThFactory["mueeWBFbb3000"] = [=](const StandardModel& SM) { return new mueeWBFbb(SM, sqrt_s_leptcoll_3000, 0., 0.); };
    //
    obsThFactory["mueeWBFbb250_p80_m30"] = [=](const StandardModel& SM) { return new mueeWBFbb(SM, sqrt_s_leptcoll_250, pol_80, -pol_30); };
    obsThFactory["mueeWBFbb250_m80_p30"] = [=](const StandardModel& SM) { return new mueeWBFbb(SM, sqrt_s_leptcoll_250, -pol_80, pol_30); };
    obsThFactory["mueeWBFbb250_p80_0"] = [=](const StandardModel& SM) { return new mueeWBFbb(SM, sqrt_s_leptcoll_250, pol_80, pol_0); };
    obsThFactory["mueeWBFbb250_m80_0"] = [=](const StandardModel& SM) { return new mueeWBFbb(SM, sqrt_s_leptcoll_250, -pol_80, pol_0); };
    //
    obsThFactory["mueeWBFbb350_p80_m30"] = [=](const StandardModel& SM) { return new mueeWBFbb(SM, sqrt_s_leptcoll_350, pol_80, -pol_30); };
    obsThFactory["mueeWBFbb350_m80_p30"] = [=](const StandardModel& SM) { return new mueeWBFbb(SM, sqrt_s_leptcoll_350, -pol_80, pol_30); };
    obsThFactory["mueeWBFbb350_p80_0"] = [=](const StandardModel& SM) { return new mueeWBFbb(SM, sqrt_s_leptcoll_350, pol_80, pol_0); };
    obsThFactory["mueeWBFbb350_m80_0"] = [=](const StandardModel& SM) { return new mueeWBFbb(SM, sqrt_s_leptcoll_350, -pol_80, pol_0); };
    //
    obsThFactory["mueeWBFbb380_p80_m30"] = [=](const StandardModel& SM) { return new mueeWBFbb(SM, sqrt_s_leptcoll_380, pol_80, -pol_30); };
    obsThFactory["mueeWBFbb380_m80_p30"] = [=](const StandardModel& SM) { return new mueeWBFbb(SM, sqrt_s_leptcoll_380, -pol_80, pol_30); };
    obsThFactory["mueeWBFbb380_p80_0"] = [=](const StandardModel& SM) { return new mueeWBFbb(SM, sqrt_s_leptcoll_380, pol_80, pol_0); };
    obsThFactory["mueeWBFbb380_m80_0"] = [=](const StandardModel& SM) { return new mueeWBFbb(SM, sqrt_s_leptcoll_380, -pol_80, pol_0); };
    //
    obsThFactory["mueeWBFbb500_p80_m30"] = [=](const StandardModel& SM) { return new mueeWBFbb(SM, sqrt_s_leptcoll_500, pol_80, -pol_30); };
    obsThFactory["mueeWBFbb500_m80_p30"] = [=](const StandardModel& SM) { return new mueeWBFbb(SM, sqrt_s_leptcoll_500, -pol_80, pol_30); };
    obsThFactory["mueeWBFbb500_p80_0"] = [=](const StandardModel& SM) { return new mueeWBFbb(SM, sqrt_s_leptcoll_500, pol_80, pol_0); };
    obsThFactory["mueeWBFbb500_m80_0"] = [=](const StandardModel& SM) { return new mueeWBFbb(SM, sqrt_s_leptcoll_500, -pol_80, pol_0); };
    //
    obsThFactory["mueeWBFbb550_p80_m30"] = [=](const StandardModel& SM) { return new mueeWBFbb(SM, sqrt_s_leptcoll_550, pol_80, -pol_30); };
    obsThFactory["mueeWBFbb550_m80_p30"] = [=](const StandardModel& SM) { return new mueeWBFbb(SM, sqrt_s_leptcoll_550, -pol_80, pol_30); };
    obsThFactory["mueeWBFbb550_p80_0"] = [=](const StandardModel& SM) { return new mueeWBFbb(SM, sqrt_s_leptcoll_550, pol_80, pol_0); };
    obsThFactory["mueeWBFbb550_m80_0"] = [=](const StandardModel& SM) { return new mueeWBFbb(SM, sqrt_s_leptcoll_550, -pol_80, pol_0); };
    //
    obsThFactory["mueeWBFbb1000_p80_m30"] = [=](const StandardModel& SM) { return new mueeWBFbb(SM, sqrt_s_leptcoll_1000, pol_80, -pol_30); };
    obsThFactory["mueeWBFbb1000_m80_p30"] = [=](const StandardModel& SM) { return new mueeWBFbb(SM, sqrt_s_leptcoll_1000, -pol_80, pol_30); };
    obsThFactory["mueeWBFbb1000_p80_m20"] = [=](const StandardModel& SM) { return new mueeWBFbb(SM, sqrt_s_leptcoll_1000, pol_80, -pol_20); };
    obsThFactory["mueeWBFbb1000_m80_p20"] = [=](const StandardModel& SM) { return new mueeWBFbb(SM, sqrt_s_leptcoll_1000, -pol_80, pol_20); };
    obsThFactory["mueeWBFbb1000_p80_0"] = [=](const StandardModel& SM) { return new mueeWBFbb(SM, sqrt_s_leptcoll_1000, pol_80, pol_0); };
    obsThFactory["mueeWBFbb1000_m80_0"] = [=](const StandardModel& SM) { return new mueeWBFbb(SM, sqrt_s_leptcoll_1000, -pol_80, pol_0); };
    //
    obsThFactory["mueeWBFbb1400_p80_m30"] = [=](const StandardModel& SM) { return new mueeWBFbb(SM, sqrt_s_leptcoll_1400, pol_80, -pol_30); };
    obsThFactory["mueeWBFbb1400_m80_p30"] = [=](const StandardModel& SM) { return new mueeWBFbb(SM, sqrt_s_leptcoll_1400, -pol_80, pol_30); };
    obsThFactory["mueeWBFbb1400_p80_0"] = [=](const StandardModel& SM) { return new mueeWBFbb(SM, sqrt_s_leptcoll_1400, pol_80, pol_0); };
    obsThFactory["mueeWBFbb1400_m80_0"] = [=](const StandardModel& SM) { return new mueeWBFbb(SM, sqrt_s_leptcoll_1400, -pol_80, pol_0); };
    //
    obsThFactory["mueeWBFbb1500_p80_m30"] = [=](const StandardModel& SM) { return new mueeWBFbb(SM, sqrt_s_leptcoll_1500, pol_80, -pol_30); };
    obsThFactory["mueeWBFbb1500_m80_p30"] = [=](const StandardModel& SM) { return new mueeWBFbb(SM, sqrt_s_leptcoll_1500, -pol_80, pol_30); };
    obsThFactory["mueeWBFbb1500_p80_0"] = [=](const StandardModel& SM) { return new mueeWBFbb(SM, sqrt_s_leptcoll_1500, pol_80, pol_0); };
    obsThFactory["mueeWBFbb1500_m80_0"] = [=](const StandardModel& SM) { return new mueeWBFbb(SM, sqrt_s_leptcoll_1500, -pol_80, pol_0); };
    //
    obsThFactory["mueeWBFbb3000_p80_m30"] = [=](const StandardModel& SM) { return new mueeWBFbb(SM, sqrt_s_leptcoll_3000, pol_80, -pol_30); };
    obsThFactory["mueeWBFbb3000_m80_p30"] = [=](const StandardModel& SM) { return new mueeWBFbb(SM, sqrt_s_leptcoll_3000, -pol_80, pol_30); };
    obsThFactory["mueeWBFbb3000_p80_0"] = [=](const StandardModel& SM) { return new mueeWBFbb(SM, sqrt_s_leptcoll_3000, pol_80, pol_0); };
    obsThFactory["mueeWBFbb3000_m80_0"] = [=](const StandardModel& SM) { return new mueeWBFbb(SM, sqrt_s_leptcoll_3000, -pol_80, pol_0); };
    //
    obsThFactory["mueeWBFcc230"] = [=](const StandardModel& SM) { return new mueeWBFcc(SM, sqrt_s_leptcoll_230, 0., 0.); };
    obsThFactory["mueeWBFcc240"] = [=](const StandardModel& SM) { return new mueeWBFcc(SM, sqrt_s_leptcoll_240, 0., 0.); };
    obsThFactory["mueeWBFcc250"] = [=](const StandardModel& SM) { return new mueeWBFcc(SM, sqrt_s_leptcoll_250, 0., 0.); };
    obsThFactory["mueeWBFcc350"] = [=](const StandardModel& SM) { return new mueeWBFcc(SM, sqrt_s_leptcoll_350, 0., 0.); };
    obsThFactory["mueeWBFcc365"] = [=](const StandardModel& SM) { return new mueeWBFcc(SM, sqrt_s_leptcoll_365, 0., 0.); };
    obsThFactory["mueeWBFcc380"] = [=](const StandardModel& SM) { return new mueeWBFcc(SM, sqrt_s_leptcoll_380, 0., 0.); };
    obsThFactory["mueeWBFcc500"] = [=](const StandardModel& SM) { return new mueeWBFcc(SM, sqrt_s_leptcoll_500, 0., 0.); };
    obsThFactory["mueeWBFcc550"] = [=](const StandardModel& SM) { return new mueeWBFcc(SM, sqrt_s_leptcoll_550, 0., 0.); };
    obsThFactory["mueeWBFcc1000"] = [=](const StandardModel& SM) { return new mueeWBFcc(SM, sqrt_s_leptcoll_1000, 0., 0.); };
    obsThFactory["mueeWBFcc1400"] = [=](const StandardModel& SM) { return new mueeWBFcc(SM, sqrt_s_leptcoll_1400, 0., 0.); };
    obsThFactory["mueeWBFcc1500"] = [=](const StandardModel& SM) { return new mueeWBFcc(SM, sqrt_s_leptcoll_1500, 0., 0.); };
    obsThFactory["mueeWBFcc3000"] = [=](const StandardModel& SM) { return new mueeWBFcc(SM, sqrt_s_leptcoll_3000, 0., 0.); };
    //
//    obsThFactory["mueeWBFcc250_p80_m30"] = bind(boost::factory<mueeWBFcc*>(), _1, sqrt_s_leptcoll_250, pol_80, -pol_30);
//    obsThFactory["mueeWBFcc250_m80_p30"] = bind(boost::factory<mueeWBFcc*>(), _1, sqrt_s_leptcoll_250, -pol_80, pol_30);
//    obsThFactory["mueeWBFcc250_p80_0"] = bind(boost::factory<mueeWBFcc*>(), _1, sqrt_s_leptcoll_250, pol_80, pol_0);
//    obsThFactory["mueeWBFcc250_m80_0"] = bind(boost::factory<mueeWBFcc*>(), _1, sqrt_s_leptcoll_250, -pol_80, pol_0);
    //
//    obsThFactory["mueeWBFcc350_p80_m30"] = bind(boost::factory<mueeWBFcc*>(), _1, sqrt_s_leptcoll_350, pol_80, -pol_30);
//    obsThFactory["mueeWBFcc350_m80_p30"] = bind(boost::factory<mueeWBFcc*>(), _1, sqrt_s_leptcoll_350, -pol_80, pol_30);
//    obsThFactory["mueeWBFcc350_p80_0"] = bind(boost::factory<mueeWBFcc*>(), _1, sqrt_s_leptcoll_350, pol_80, pol_0);
//    obsThFactory["mueeWBFcc350_m80_0"] = bind(boost::factory<mueeWBFcc*>(), _1, sqrt_s_leptcoll_350, -pol_80, pol_0);
    //
//    obsThFactory["mueeWBFcc380_p80_m30"] = bind(boost::factory<mueeWBFcc*>(), _1, sqrt_s_leptcoll_380, pol_80, -pol_30);
//    obsThFactory["mueeWBFcc380_m80_p30"] = bind(boost::factory<mueeWBFcc*>(), _1, sqrt_s_leptcoll_380, -pol_80, pol_30);
//    obsThFactory["mueeWBFcc380_p80_0"] = bind(boost::factory<mueeWBFcc*>(), _1, sqrt_s_leptcoll_380, pol_80, pol_0);
//    obsThFactory["mueeWBFcc380_m80_0"] = bind(boost::factory<mueeWBFcc*>(), _1, sqrt_s_leptcoll_380, -pol_80, pol_0);
    //
//    obsThFactory["mueeWBFcc500_p80_m30"] = bind(boost::factory<mueeWBFcc*>(), _1, sqrt_s_leptcoll_500, pol_80, -pol_30);
//    obsThFactory["mueeWBFcc500_m80_p30"] = bind(boost::factory<mueeWBFcc*>(), _1, sqrt_s_leptcoll_500, -pol_80, pol_30);
//    obsThFactory["mueeWBFcc500_p80_0"] = bind(boost::factory<mueeWBFcc*>(), _1, sqrt_s_leptcoll_500, pol_80, pol_0);
//    obsThFactory["mueeWBFcc500_m80_0"] = bind(boost::factory<mueeWBFcc*>(), _1, sqrt_s_leptcoll_500, -pol_80, pol_0);
    //
//    obsThFactory["mueeWBFcc1000_p80_m30"] = bind(boost::factory<mueeWBFcc*>(), _1, sqrt_s_leptcoll_1000, pol_80, -pol_30);
//    obsThFactory["mueeWBFcc1000_m80_p30"] = bind(boost::factory<mueeWBFcc*>(), _1, sqrt_s_leptcoll_1000, -pol_80, pol_30);
//    obsThFactory["mueeWBFcc1000_p80_m20"] = bind(boost::factory<mueeWBFcc*>(), _1, sqrt_s_leptcoll_1000, pol_80, -pol_20);
//    obsThFactory["mueeWBFcc1000_m80_p20"] = bind(boost::factory<mueeWBFcc*>(), _1, sqrt_s_leptcoll_1000, -pol_80, pol_20);
//    obsThFactory["mueeWBFcc1000_p80_0"] = bind(boost::factory<mueeWBFcc*>(), _1, sqrt_s_leptcoll_1000, pol_80, pol_0);
//    obsThFactory["mueeWBFcc1000_m80_0"] = bind(boost::factory<mueeWBFcc*>(), _1, sqrt_s_leptcoll_1000, -pol_80, pol_0);
    //
//    obsThFactory["mueeWBFcc1400_p80_m30"] = bind(boost::factory<mueeWBFcc*>(), _1, sqrt_s_leptcoll_1400, pol_80, -pol_30);
//    obsThFactory["mueeWBFcc1400_m80_p30"] = bind(boost::factory<mueeWBFcc*>(), _1, sqrt_s_leptcoll_1400, -pol_80, pol_30);
//    obsThFactory["mueeWBFcc1400_p80_0"] = bind(boost::factory<mueeWBFcc*>(), _1, sqrt_s_leptcoll_1400, pol_80, pol_0);
//    obsThFactory["mueeWBFcc1400_m80_0"] = bind(boost::factory<mueeWBFcc*>(), _1, sqrt_s_leptcoll_1400, -pol_80, pol_0);
    //
//    obsThFactory["mueeWBFcc1500_p80_m30"] = bind(boost::factory<mueeWBFcc*>(), _1, sqrt_s_leptcoll_1500, pol_80, -pol_30);
//    obsThFactory["mueeWBFcc1500_m80_p30"] = bind(boost::factory<mueeWBFcc*>(), _1, sqrt_s_leptcoll_1500, -pol_80, pol_30);
//    obsThFactory["mueeWBFcc1500_p80_0"] = bind(boost::factory<mueeWBFcc*>(), _1, sqrt_s_leptcoll_1500, pol_80, pol_0);
//    obsThFactory["mueeWBFcc1500_m80_0"] = bind(boost::factory<mueeWBFcc*>(), _1, sqrt_s_leptcoll_1500, -pol_80, pol_0);
    //
//    obsThFactory["mueeWBFcc3000_p80_m30"] = bind(boost::factory<mueeWBFcc*>(), _1, sqrt_s_leptcoll_3000, pol_80, -pol_30);
//    obsThFactory["mueeWBFcc3000_m80_p30"] = bind(boost::factory<mueeWBFcc*>(), _1, sqrt_s_leptcoll_3000, -pol_80, pol_30);
//    obsThFactory["mueeWBFcc3000_p80_0"] = bind(boost::factory<mueeWBFcc*>(), _1, sqrt_s_leptcoll_3000, pol_80, pol_0);
//    obsThFactory["mueeWBFcc3000_m80_0"] = bind(boost::factory<mueeWBFcc*>(), _1, sqrt_s_leptcoll_3000, -pol_80, pol_0);
    //
    obsThFactory["mueeWBFgg230"] = [=](const StandardModel& SM) { return new mueeWBFgg(SM, sqrt_s_leptcoll_230, 0., 0.); };
    obsThFactory["mueeWBFgg240"] = [=](const StandardModel& SM) { return new mueeWBFgg(SM, sqrt_s_leptcoll_240, 0., 0.); };
    obsThFactory["mueeWBFgg250"] = [=](const StandardModel& SM) { return new mueeWBFgg(SM, sqrt_s_leptcoll_250, 0., 0.); };
    obsThFactory["mueeWBFgg350"] = [=](const StandardModel& SM) { return new mueeWBFgg(SM, sqrt_s_leptcoll_350, 0., 0.); };
    obsThFactory["mueeWBFgg365"] = [=](const StandardModel& SM) { return new mueeWBFgg(SM, sqrt_s_leptcoll_365, 0., 0.); };
    obsThFactory["mueeWBFgg380"] = [=](const StandardModel& SM) { return new mueeWBFgg(SM, sqrt_s_leptcoll_380, 0., 0.); };
    obsThFactory["mueeWBFgg500"] = [=](const StandardModel& SM) { return new mueeWBFgg(SM, sqrt_s_leptcoll_500, 0., 0.); };
    obsThFactory["mueeWBFgg550"] = [=](const StandardModel& SM) { return new mueeWBFgg(SM, sqrt_s_leptcoll_550, 0., 0.); };
    obsThFactory["mueeWBFgg1000"] = [=](const StandardModel& SM) { return new mueeWBFgg(SM, sqrt_s_leptcoll_1000, 0., 0.); };
    obsThFactory["mueeWBFgg1400"] = [=](const StandardModel& SM) { return new mueeWBFgg(SM, sqrt_s_leptcoll_1400, 0., 0.); };
    obsThFactory["mueeWBFgg1500"] = [=](const StandardModel& SM) { return new mueeWBFgg(SM, sqrt_s_leptcoll_1500, 0., 0.); };
    obsThFactory["mueeWBFgg3000"] = [=](const StandardModel& SM) { return new mueeWBFgg(SM, sqrt_s_leptcoll_3000, 0., 0.); };
    //
//    obsThFactory["mueeWBFgg250_p80_m30"] = bind(boost::factory<mueeWBFgg*>(), _1, sqrt_s_leptcoll_250, pol_80, -pol_30);
//    obsThFactory["mueeWBFgg250_m80_p30"] = bind(boost::factory<mueeWBFgg*>(), _1, sqrt_s_leptcoll_250, -pol_80, pol_30);
//    obsThFactory["mueeWBFgg250_p80_0"] = bind(boost::factory<mueeWBFgg*>(), _1, sqrt_s_leptcoll_250, pol_80, pol_0);
//    obsThFactory["mueeWBFgg250_m80_0"] = bind(boost::factory<mueeWBFgg*>(), _1, sqrt_s_leptcoll_250, -pol_80, pol_0);
    //
//    obsThFactory["mueeWBFgg350_p80_m30"] = bind(boost::factory<mueeWBFgg*>(), _1, sqrt_s_leptcoll_350, pol_80, -pol_30);
//    obsThFactory["mueeWBFgg350_m80_p30"] = bind(boost::factory<mueeWBFgg*>(), _1, sqrt_s_leptcoll_350, -pol_80, pol_30);
//    obsThFactory["mueeWBFgg350_p80_0"] = bind(boost::factory<mueeWBFgg*>(), _1, sqrt_s_leptcoll_350, pol_80, pol_0);
//    obsThFactory["mueeWBFgg350_m80_0"] = bind(boost::factory<mueeWBFgg*>(), _1, sqrt_s_leptcoll_350, -pol_80, pol_0);
    //
//    obsThFactory["mueeWBFgg380_p80_m30"] = bind(boost::factory<mueeWBFgg*>(), _1, sqrt_s_leptcoll_380, pol_80, -pol_30);
//    obsThFactory["mueeWBFgg380_m80_p30"] = bind(boost::factory<mueeWBFgg*>(), _1, sqrt_s_leptcoll_380, -pol_80, pol_30);
//    obsThFactory["mueeWBFgg380_p80_0"] = bind(boost::factory<mueeWBFgg*>(), _1, sqrt_s_leptcoll_380, pol_80, pol_0);
//    obsThFactory["mueeWBFgg380_m80_0"] = bind(boost::factory<mueeWBFgg*>(), _1, sqrt_s_leptcoll_380, -pol_80, pol_0);
    //
//    obsThFactory["mueeWBFgg500_p80_m30"] = bind(boost::factory<mueeWBFgg*>(), _1, sqrt_s_leptcoll_500, pol_80, -pol_30);
//    obsThFactory["mueeWBFgg500_m80_p30"] = bind(boost::factory<mueeWBFgg*>(), _1, sqrt_s_leptcoll_500, -pol_80, pol_30);
//    obsThFactory["mueeWBFgg500_p80_0"] = bind(boost::factory<mueeWBFgg*>(), _1, sqrt_s_leptcoll_500, pol_80, pol_0);
//    obsThFactory["mueeWBFgg500_m80_0"] = bind(boost::factory<mueeWBFgg*>(), _1, sqrt_s_leptcoll_500, -pol_80, pol_0);
    //
//    obsThFactory["mueeWBFgg1000_p80_m30"] = bind(boost::factory<mueeWBFgg*>(), _1, sqrt_s_leptcoll_1000, pol_80, -pol_30);
//    obsThFactory["mueeWBFgg1000_m80_p30"] = bind(boost::factory<mueeWBFgg*>(), _1, sqrt_s_leptcoll_1000, -pol_80, pol_30);
//    obsThFactory["mueeWBFgg1000_p80_m20"] = bind(boost::factory<mueeWBFgg*>(), _1, sqrt_s_leptcoll_1000, pol_80, -pol_20);
//    obsThFactory["mueeWBFgg1000_m80_p20"] = bind(boost::factory<mueeWBFgg*>(), _1, sqrt_s_leptcoll_1000, -pol_80, pol_20);
//    obsThFactory["mueeWBFgg1000_p80_0"] = bind(boost::factory<mueeWBFgg*>(), _1, sqrt_s_leptcoll_1000, pol_80, pol_0);
//    obsThFactory["mueeWBFgg1000_m80_0"] = bind(boost::factory<mueeWBFgg*>(), _1, sqrt_s_leptcoll_1000, -pol_80, pol_0);
    //
//    obsThFactory["mueeWBFgg1400_p80_m30"] = bind(boost::factory<mueeWBFgg*>(), _1, sqrt_s_leptcoll_1400, pol_80, -pol_30);
//    obsThFactory["mueeWBFgg1400_m80_p30"] = bind(boost::factory<mueeWBFgg*>(), _1, sqrt_s_leptcoll_1400, -pol_80, pol_30);
//    obsThFactory["mueeWBFgg1400_p80_0"] = bind(boost::factory<mueeWBFgg*>(), _1, sqrt_s_leptcoll_1400, pol_80, pol_0);
//    obsThFactory["mueeWBFgg1400_m80_0"] = bind(boost::factory<mueeWBFgg*>(), _1, sqrt_s_leptcoll_1400, -pol_80, pol_0);
    //
//    obsThFactory["mueeWBFgg1500_p80_m30"] = bind(boost::factory<mueeWBFgg*>(), _1, sqrt_s_leptcoll_1500, pol_80, -pol_30);
//    obsThFactory["mueeWBFgg1500_m80_p30"] = bind(boost::factory<mueeWBFgg*>(), _1, sqrt_s_leptcoll_1500, -pol_80, pol_30);
//    obsThFactory["mueeWBFgg1500_p80_0"] = bind(boost::factory<mueeWBFgg*>(), _1, sqrt_s_leptcoll_1500, pol_80, pol_0);
//    obsThFactory["mueeWBFgg1500_m80_0"] = bind(boost::factory<mueeWBFgg*>(), _1, sqrt_s_leptcoll_1500, -pol_80, pol_0);
    //
//    obsThFactory["mueeWBFgg3000_p80_m30"] = bind(boost::factory<mueeWBFgg*>(), _1, sqrt_s_leptcoll_3000, pol_80, -pol_30);
//    obsThFactory["mueeWBFgg3000_m80_p30"] = bind(boost::factory<mueeWBFgg*>(), _1, sqrt_s_leptcoll_3000, -pol_80, pol_30);
//    obsThFactory["mueeWBFgg3000_p80_0"] = bind(boost::factory<mueeWBFgg*>(), _1, sqrt_s_leptcoll_3000, pol_80, pol_0);
//    obsThFactory["mueeWBFgg3000_m80_0"] = bind(boost::factory<mueeWBFgg*>(), _1, sqrt_s_leptcoll_3000, -pol_80, pol_0);
    //
    obsThFactory["mueeWBFWW230"] = [=](const StandardModel& SM) { return new mueeWBFWW(SM, sqrt_s_leptcoll_230, 0., 0.); };
    obsThFactory["mueeWBFWW240"] = [=](const StandardModel& SM) { return new mueeWBFWW(SM, sqrt_s_leptcoll_240, 0., 0.); };
    obsThFactory["mueeWBFWW250"] = [=](const StandardModel& SM) { return new mueeWBFWW(SM, sqrt_s_leptcoll_250, 0., 0.); };
    obsThFactory["mueeWBFWW350"] = [=](const StandardModel& SM) { return new mueeWBFWW(SM, sqrt_s_leptcoll_350, 0., 0.); };
    obsThFactory["mueeWBFWW365"] = [=](const StandardModel& SM) { return new mueeWBFWW(SM, sqrt_s_leptcoll_365, 0., 0.); };
    obsThFactory["mueeWBFWW380"] = [=](const StandardModel& SM) { return new mueeWBFWW(SM, sqrt_s_leptcoll_380, 0., 0.); };
    obsThFactory["mueeWBFWW500"] = [=](const StandardModel& SM) { return new mueeWBFWW(SM, sqrt_s_leptcoll_500, 0., 0.); };
    obsThFactory["mueeWBFWW550"] = [=](const StandardModel& SM) { return new mueeWBFWW(SM, sqrt_s_leptcoll_550, 0., 0.); };
    obsThFactory["mueeWBFWW1000"] = [=](const StandardModel& SM) { return new mueeWBFWW(SM, sqrt_s_leptcoll_1000, 0., 0.); };
    obsThFactory["mueeWBFWW1400"] = [=](const StandardModel& SM) { return new mueeWBFWW(SM, sqrt_s_leptcoll_1400, 0., 0.); };
    obsThFactory["mueeWBFWW1500"] = [=](const StandardModel& SM) { return new mueeWBFWW(SM, sqrt_s_leptcoll_1500, 0., 0.); };
    obsThFactory["mueeWBFWW3000"] = [=](const StandardModel& SM) { return new mueeWBFWW(SM, sqrt_s_leptcoll_3000, 0., 0.); };
    //
//    obsThFactory["mueeWBFWW250_p80_m30"] = bind(boost::factory<mueeWBFWW*>(), _1, sqrt_s_leptcoll_250, pol_80, -pol_30);
//    obsThFactory["mueeWBFWW250_m80_p30"] = bind(boost::factory<mueeWBFWW*>(), _1, sqrt_s_leptcoll_250, -pol_80, pol_30);
//    obsThFactory["mueeWBFWW250_p80_0"] = bind(boost::factory<mueeWBFWW*>(), _1, sqrt_s_leptcoll_250, pol_80, pol_0);
//    obsThFactory["mueeWBFWW250_m80_0"] = bind(boost::factory<mueeWBFWW*>(), _1, sqrt_s_leptcoll_250, -pol_80, pol_0);
    //
//    obsThFactory["mueeWBFWW350_p80_m30"] = bind(boost::factory<mueeWBFWW*>(), _1, sqrt_s_leptcoll_350, pol_80, -pol_30);
//    obsThFactory["mueeWBFWW350_m80_p30"] = bind(boost::factory<mueeWBFWW*>(), _1, sqrt_s_leptcoll_350, -pol_80, pol_30);
//    obsThFactory["mueeWBFWW350_p80_0"] = bind(boost::factory<mueeWBFWW*>(), _1, sqrt_s_leptcoll_350, pol_80, pol_0);
//    obsThFactory["mueeWBFWW350_m80_0"] = bind(boost::factory<mueeWBFWW*>(), _1, sqrt_s_leptcoll_350, -pol_80, pol_0);
    //
//    obsThFactory["mueeWBFWW380_p80_m30"] = bind(boost::factory<mueeWBFWW*>(), _1, sqrt_s_leptcoll_380, pol_80, -pol_30);
//    obsThFactory["mueeWBFWW380_m80_p30"] = bind(boost::factory<mueeWBFWW*>(), _1, sqrt_s_leptcoll_380, -pol_80, pol_30);
//    obsThFactory["mueeWBFWW380_p80_0"] = bind(boost::factory<mueeWBFWW*>(), _1, sqrt_s_leptcoll_380, pol_80, pol_0);
//    obsThFactory["mueeWBFWW380_m80_0"] = bind(boost::factory<mueeWBFWW*>(), _1, sqrt_s_leptcoll_380, -pol_80, pol_0);
    //
//    obsThFactory["mueeWBFWW500_p80_m30"] = bind(boost::factory<mueeWBFWW*>(), _1, sqrt_s_leptcoll_500, pol_80, -pol_30);
//    obsThFactory["mueeWBFWW500_m80_p30"] = bind(boost::factory<mueeWBFWW*>(), _1, sqrt_s_leptcoll_500, -pol_80, pol_30);
//    obsThFactory["mueeWBFWW500_p80_0"] = bind(boost::factory<mueeWBFWW*>(), _1, sqrt_s_leptcoll_500, pol_80, pol_0);
//    obsThFactory["mueeWBFWW500_m80_0"] = bind(boost::factory<mueeWBFWW*>(), _1, sqrt_s_leptcoll_500, -pol_80, pol_0);
    //
//    obsThFactory["mueeWBFWW1000_p80_m30"] = bind(boost::factory<mueeWBFWW*>(), _1, sqrt_s_leptcoll_1000, pol_80, -pol_30);
//    obsThFactory["mueeWBFWW1000_m80_p30"] = bind(boost::factory<mueeWBFWW*>(), _1, sqrt_s_leptcoll_1000, -pol_80, pol_30);
//    obsThFactory["mueeWBFWW1000_p80_m20"] = bind(boost::factory<mueeWBFWW*>(), _1, sqrt_s_leptcoll_1000, pol_80, -pol_20);
//    obsThFactory["mueeWBFWW1000_m80_p20"] = bind(boost::factory<mueeWBFWW*>(), _1, sqrt_s_leptcoll_1000, -pol_80, pol_20);
//    obsThFactory["mueeWBFWW1000_p80_0"] = bind(boost::factory<mueeWBFWW*>(), _1, sqrt_s_leptcoll_1000, pol_80, pol_0);
//    obsThFactory["mueeWBFWW1000_m80_0"] = bind(boost::factory<mueeWBFWW*>(), _1, sqrt_s_leptcoll_1000, -pol_80, pol_0);
    //
//    obsThFactory["mueeWBFWW1400_p80_m30"] = bind(boost::factory<mueeWBFWW*>(), _1, sqrt_s_leptcoll_1400, pol_80, -pol_30);
//    obsThFactory["mueeWBFWW1400_m80_p30"] = bind(boost::factory<mueeWBFWW*>(), _1, sqrt_s_leptcoll_1400, -pol_80, pol_30);
//    obsThFactory["mueeWBFWW1400_p80_0"] = bind(boost::factory<mueeWBFWW*>(), _1, sqrt_s_leptcoll_1400, pol_80, pol_0);
//    obsThFactory["mueeWBFWW1400_m80_0"] = bind(boost::factory<mueeWBFWW*>(), _1, sqrt_s_leptcoll_1400, -pol_80, pol_0);
    //
//    obsThFactory["mueeWBFWW1500_p80_m30"] = bind(boost::factory<mueeWBFWW*>(), _1, sqrt_s_leptcoll_1500, pol_80, -pol_30);
//    obsThFactory["mueeWBFWW1500_m80_p30"] = bind(boost::factory<mueeWBFWW*>(), _1, sqrt_s_leptcoll_1500, -pol_80, pol_30);
//    obsThFactory["mueeWBFWW1500_p80_0"] = bind(boost::factory<mueeWBFWW*>(), _1, sqrt_s_leptcoll_1500, pol_80, pol_0);
//    obsThFactory["mueeWBFWW1500_m80_0"] = bind(boost::factory<mueeWBFWW*>(), _1, sqrt_s_leptcoll_1500, -pol_80, pol_0);
    //
//    obsThFactory["mueeWBFWW3000_p80_m30"] = bind(boost::factory<mueeWBFWW*>(), _1, sqrt_s_leptcoll_3000, pol_80, -pol_30);
//    obsThFactory["mueeWBFWW3000_m80_p30"] = bind(boost::factory<mueeWBFWW*>(), _1, sqrt_s_leptcoll_3000, -pol_80, pol_30);
//    obsThFactory["mueeWBFWW3000_p80_0"] = bind(boost::factory<mueeWBFWW*>(), _1, sqrt_s_leptcoll_3000, pol_80, pol_0);
//    obsThFactory["mueeWBFWW3000_m80_0"] = bind(boost::factory<mueeWBFWW*>(), _1, sqrt_s_leptcoll_3000, -pol_80, pol_0);
    //
    obsThFactory["mueeWBFtautau230"] = [=](const StandardModel& SM) { return new mueeWBFtautau(SM, sqrt_s_leptcoll_230, 0., 0.); };
    obsThFactory["mueeWBFtautau240"] = [=](const StandardModel& SM) { return new mueeWBFtautau(SM, sqrt_s_leptcoll_240, 0., 0.); };
    obsThFactory["mueeWBFtautau250"] = [=](const StandardModel& SM) { return new mueeWBFtautau(SM, sqrt_s_leptcoll_250, 0., 0.); };
    obsThFactory["mueeWBFtautau350"] = [=](const StandardModel& SM) { return new mueeWBFtautau(SM, sqrt_s_leptcoll_350, 0., 0.); };
    obsThFactory["mueeWBFtautau365"] = [=](const StandardModel& SM) { return new mueeWBFtautau(SM, sqrt_s_leptcoll_365, 0., 0.); };
    obsThFactory["mueeWBFtautau380"] = [=](const StandardModel& SM) { return new mueeWBFtautau(SM, sqrt_s_leptcoll_380, 0., 0.); };
    obsThFactory["mueeWBFtautau500"] = [=](const StandardModel& SM) { return new mueeWBFtautau(SM, sqrt_s_leptcoll_500, 0., 0.); };
    obsThFactory["mueeWBFtautau550"] = [=](const StandardModel& SM) { return new mueeWBFtautau(SM, sqrt_s_leptcoll_550, 0., 0.); };    
    obsThFactory["mueeWBFtautau1000"] = [=](const StandardModel& SM) { return new mueeWBFtautau(SM, sqrt_s_leptcoll_1000, 0., 0.); };
    obsThFactory["mueeWBFtautau1400"] = [=](const StandardModel& SM) { return new mueeWBFtautau(SM, sqrt_s_leptcoll_1400, 0., 0.); };
    obsThFactory["mueeWBFtautau1500"] = [=](const StandardModel& SM) { return new mueeWBFtautau(SM, sqrt_s_leptcoll_1500, 0., 0.); };
    obsThFactory["mueeWBFtautau3000"] = [=](const StandardModel& SM) { return new mueeWBFtautau(SM, sqrt_s_leptcoll_3000, 0., 0.); };
    //
//    obsThFactory["mueeWBFtautau250_p80_m30"] = bind(boost::factory<mueeWBFtautau*>(), _1, sqrt_s_leptcoll_250, pol_80, -pol_30);
//    obsThFactory["mueeWBFtautau250_m80_p30"] = bind(boost::factory<mueeWBFtautau*>(), _1, sqrt_s_leptcoll_250, -pol_80, pol_30);
//    obsThFactory["mueeWBFtautau250_p80_0"] = bind(boost::factory<mueeWBFtautau*>(), _1, sqrt_s_leptcoll_250, pol_80, pol_0);
//    obsThFactory["mueeWBFtautau250_m80_0"] = bind(boost::factory<mueeWBFtautau*>(), _1, sqrt_s_leptcoll_250, -pol_80, pol_0);
    //
//    obsThFactory["mueeWBFtautau350_p80_m30"] = bind(boost::factory<mueeWBFtautau*>(), _1, sqrt_s_leptcoll_350, pol_80, -pol_30);
//    obsThFactory["mueeWBFtautau350_m80_p30"] = bind(boost::factory<mueeWBFtautau*>(), _1, sqrt_s_leptcoll_350, -pol_80, pol_30);
//    obsThFactory["mueeWBFtautau350_p80_0"] = bind(boost::factory<mueeWBFtautau*>(), _1, sqrt_s_leptcoll_350, pol_80, pol_0);
//    obsThFactory["mueeWBFtautau350_m80_0"] = bind(boost::factory<mueeWBFtautau*>(), _1, sqrt_s_leptcoll_350, -pol_80, pol_0);
    //
//    obsThFactory["mueeWBFtautau380_p80_m30"] = bind(boost::factory<mueeWBFtautau*>(), _1, sqrt_s_leptcoll_380, pol_80, -pol_30);
//    obsThFactory["mueeWBFtautau380_m80_p30"] = bind(boost::factory<mueeWBFtautau*>(), _1, sqrt_s_leptcoll_380, -pol_80, pol_30);
//    obsThFactory["mueeWBFtautau380_p80_0"] = bind(boost::factory<mueeWBFtautau*>(), _1, sqrt_s_leptcoll_380, pol_80, pol_0);
//    obsThFactory["mueeWBFtautau380_m80_0"] = bind(boost::factory<mueeWBFtautau*>(), _1, sqrt_s_leptcoll_380, -pol_80, pol_0);
    //
//    obsThFactory["mueeWBFtautau500_p80_m30"] = bind(boost::factory<mueeWBFtautau*>(), _1, sqrt_s_leptcoll_500, pol_80, -pol_30);
//    obsThFactory["mueeWBFtautau500_m80_p30"] = bind(boost::factory<mueeWBFtautau*>(), _1, sqrt_s_leptcoll_500, -pol_80, pol_30);
//    obsThFactory["mueeWBFtautau500_p80_0"] = bind(boost::factory<mueeWBFtautau*>(), _1, sqrt_s_leptcoll_500, pol_80, pol_0);
//    obsThFactory["mueeWBFtautau500_m80_0"] = bind(boost::factory<mueeWBFtautau*>(), _1, sqrt_s_leptcoll_500, -pol_80, pol_0);
    //
//    obsThFactory["mueeWBFtautau1000_p80_m30"] = bind(boost::factory<mueeWBFtautau*>(), _1, sqrt_s_leptcoll_1000, pol_80, -pol_30);
//    obsThFactory["mueeWBFtautau1000_m80_p30"] = bind(boost::factory<mueeWBFtautau*>(), _1, sqrt_s_leptcoll_1000, -pol_80, pol_30);
//    obsThFactory["mueeWBFtautau1000_p80_m20"] = bind(boost::factory<mueeWBFtautau*>(), _1, sqrt_s_leptcoll_1000, pol_80, -pol_20);
//    obsThFactory["mueeWBFtautau1000_m80_p20"] = bind(boost::factory<mueeWBFtautau*>(), _1, sqrt_s_leptcoll_1000, -pol_80, pol_20);
//    obsThFactory["mueeWBFtautau1000_p80_0"] = bind(boost::factory<mueeWBFtautau*>(), _1, sqrt_s_leptcoll_1000, pol_80, pol_0);
//    obsThFactory["mueeWBFtautau1000_m80_0"] = bind(boost::factory<mueeWBFtautau*>(), _1, sqrt_s_leptcoll_1000, -pol_80, pol_0);
    //
//    obsThFactory["mueeWBFtautau1400_p80_m30"] = bind(boost::factory<mueeWBFtautau*>(), _1, sqrt_s_leptcoll_1400, pol_80, -pol_30);
//    obsThFactory["mueeWBFtautau1400_m80_p30"] = bind(boost::factory<mueeWBFtautau*>(), _1, sqrt_s_leptcoll_1400, -pol_80, pol_30);
//    obsThFactory["mueeWBFtautau1400_p80_0"] = bind(boost::factory<mueeWBFtautau*>(), _1, sqrt_s_leptcoll_1400, pol_80, pol_0);
//    obsThFactory["mueeWBFtautau1400_m80_0"] = bind(boost::factory<mueeWBFtautau*>(), _1, sqrt_s_leptcoll_1400, -pol_80, pol_0);
    //
//    obsThFactory["mueeWBFtautau1500_p80_m30"] = bind(boost::factory<mueeWBFtautau*>(), _1, sqrt_s_leptcoll_1500, pol_80, -pol_30);
//    obsThFactory["mueeWBFtautau1500_m80_p30"] = bind(boost::factory<mueeWBFtautau*>(), _1, sqrt_s_leptcoll_1500, -pol_80, pol_30);
//    obsThFactory["mueeWBFtautau1500_p80_0"] = bind(boost::factory<mueeWBFtautau*>(), _1, sqrt_s_leptcoll_1500, pol_80, pol_0);
//    obsThFactory["mueeWBFtautau1500_m80_0"] = bind(boost::factory<mueeWBFtautau*>(), _1, sqrt_s_leptcoll_1500, -pol_80, pol_0);
    //
//    obsThFactory["mueeWBFtautau3000_p80_m30"] = bind(boost::factory<mueeWBFtautau*>(), _1, sqrt_s_leptcoll_3000, pol_80, -pol_30);
//    obsThFactory["mueeWBFtautau3000_m80_p30"] = bind(boost::factory<mueeWBFtautau*>(), _1, sqrt_s_leptcoll_3000, -pol_80, pol_30);
//    obsThFactory["mueeWBFtautau3000_p80_0"] = bind(boost::factory<mueeWBFtautau*>(), _1, sqrt_s_leptcoll_3000, pol_80, pol_0);
//    obsThFactory["mueeWBFtautau3000_m80_0"] = bind(boost::factory<mueeWBFtautau*>(), _1, sqrt_s_leptcoll_3000, -pol_80, pol_0);
    //
    obsThFactory["mueeWBFZZ230"] = [=](const StandardModel& SM) { return new mueeWBFZZ(SM, sqrt_s_leptcoll_230, 0., 0.); };
    obsThFactory["mueeWBFZZ240"] = [=](const StandardModel& SM) { return new mueeWBFZZ(SM, sqrt_s_leptcoll_240, 0., 0.); };
    obsThFactory["mueeWBFZZ250"] = [=](const StandardModel& SM) { return new mueeWBFZZ(SM, sqrt_s_leptcoll_250, 0., 0.); };
    obsThFactory["mueeWBFZZ350"] = [=](const StandardModel& SM) { return new mueeWBFZZ(SM, sqrt_s_leptcoll_350, 0., 0.); };
    obsThFactory["mueeWBFZZ365"] = [=](const StandardModel& SM) { return new mueeWBFZZ(SM, sqrt_s_leptcoll_365, 0., 0.); };
    obsThFactory["mueeWBFZZ380"] = [=](const StandardModel& SM) { return new mueeWBFZZ(SM, sqrt_s_leptcoll_380, 0., 0.); };
    obsThFactory["mueeWBFZZ500"] = [=](const StandardModel& SM) { return new mueeWBFZZ(SM, sqrt_s_leptcoll_500, 0., 0.); };
    obsThFactory["mueeWBFZZ550"] = [=](const StandardModel& SM) { return new mueeWBFZZ(SM, sqrt_s_leptcoll_550, 0., 0.); };
    obsThFactory["mueeWBFZZ1000"] = [=](const StandardModel& SM) { return new mueeWBFZZ(SM, sqrt_s_leptcoll_1000, 0., 0.); };
    obsThFactory["mueeWBFZZ1400"] = [=](const StandardModel& SM) { return new mueeWBFZZ(SM, sqrt_s_leptcoll_1400, 0., 0.); };
    obsThFactory["mueeWBFZZ1500"] = [=](const StandardModel& SM) { return new mueeWBFZZ(SM, sqrt_s_leptcoll_1500, 0., 0.); };
    obsThFactory["mueeWBFZZ3000"] = [=](const StandardModel& SM) { return new mueeWBFZZ(SM, sqrt_s_leptcoll_3000, 0., 0.); };
    //
//    obsThFactory["mueeWBFZZ250_p80_m30"] = bind(boost::factory<mueeWBFZZ*>(), _1, sqrt_s_leptcoll_250, pol_80, -pol_30);
//    obsThFactory["mueeWBFZZ250_m80_p30"] = bind(boost::factory<mueeWBFZZ*>(), _1, sqrt_s_leptcoll_250, -pol_80, pol_30);
//    obsThFactory["mueeWBFZZ250_p80_0"] = bind(boost::factory<mueeWBFZZ*>(), _1, sqrt_s_leptcoll_250, pol_80, pol_0);
//    obsThFactory["mueeWBFZZ250_m80_0"] = bind(boost::factory<mueeWBFZZ*>(), _1, sqrt_s_leptcoll_250, -pol_80, pol_0);
    //
//    obsThFactory["mueeWBFZZ350_p80_m30"] = bind(boost::factory<mueeWBFZZ*>(), _1, sqrt_s_leptcoll_350, pol_80, -pol_30);
//    obsThFactory["mueeWBFZZ350_m80_p30"] = bind(boost::factory<mueeWBFZZ*>(), _1, sqrt_s_leptcoll_350, -pol_80, pol_30);
//    obsThFactory["mueeWBFZZ350_p80_0"] = bind(boost::factory<mueeWBFZZ*>(), _1, sqrt_s_leptcoll_350, pol_80, pol_0);
//    obsThFactory["mueeWBFZZ350_m80_0"] = bind(boost::factory<mueeWBFZZ*>(), _1, sqrt_s_leptcoll_350, -pol_80, pol_0);
    //
//    obsThFactory["mueeWBFZZ380_p80_m30"] = bind(boost::factory<mueeWBFZZ*>(), _1, sqrt_s_leptcoll_380, pol_80, -pol_30);
//    obsThFactory["mueeWBFZZ380_m80_p30"] = bind(boost::factory<mueeWBFZZ*>(), _1, sqrt_s_leptcoll_380, -pol_80, pol_30);
//    obsThFactory["mueeWBFZZ380_p80_0"] = bind(boost::factory<mueeWBFZZ*>(), _1, sqrt_s_leptcoll_380, pol_80, pol_0);
//    obsThFactory["mueeWBFZZ380_m80_0"] = bind(boost::factory<mueeWBFZZ*>(), _1, sqrt_s_leptcoll_380, -pol_80, pol_0);
    //
//    obsThFactory["mueeWBFZZ500_p80_m30"] = bind(boost::factory<mueeWBFZZ*>(), _1, sqrt_s_leptcoll_500, pol_80, -pol_30);
//    obsThFactory["mueeWBFZZ500_m80_p30"] = bind(boost::factory<mueeWBFZZ*>(), _1, sqrt_s_leptcoll_500, -pol_80, pol_30);
//    obsThFactory["mueeWBFZZ500_p80_0"] = bind(boost::factory<mueeWBFZZ*>(), _1, sqrt_s_leptcoll_500, pol_80, pol_0);
//    obsThFactory["mueeWBFZZ500_m80_0"] = bind(boost::factory<mueeWBFZZ*>(), _1, sqrt_s_leptcoll_500, -pol_80, pol_0);
    //
//    obsThFactory["mueeWBFZZ1000_p80_m30"] = bind(boost::factory<mueeWBFZZ*>(), _1, sqrt_s_leptcoll_1000, pol_80, -pol_30);
//    obsThFactory["mueeWBFZZ1000_m80_p30"] = bind(boost::factory<mueeWBFZZ*>(), _1, sqrt_s_leptcoll_1000, -pol_80, pol_30);
//    obsThFactory["mueeWBFZZ1000_p80_m20"] = bind(boost::factory<mueeWBFZZ*>(), _1, sqrt_s_leptcoll_1000, pol_80, -pol_20);
//    obsThFactory["mueeWBFZZ1000_m80_p20"] = bind(boost::factory<mueeWBFZZ*>(), _1, sqrt_s_leptcoll_1000, -pol_80, pol_20);
//    obsThFactory["mueeWBFZZ1000_p80_0"] = bind(boost::factory<mueeWBFZZ*>(), _1, sqrt_s_leptcoll_1000, pol_80, pol_0);
//    obsThFactory["mueeWBFZZ1000_m80_0"] = bind(boost::factory<mueeWBFZZ*>(), _1, sqrt_s_leptcoll_1000, -pol_80, pol_0);
    //
//    obsThFactory["mueeWBFZZ1400_p80_m30"] = bind(boost::factory<mueeWBFZZ*>(), _1, sqrt_s_leptcoll_1400, pol_80, -pol_30);
//    obsThFactory["mueeWBFZZ1400_m80_p30"] = bind(boost::factory<mueeWBFZZ*>(), _1, sqrt_s_leptcoll_1400, -pol_80, pol_30);
//    obsThFactory["mueeWBFZZ1400_p80_0"] = bind(boost::factory<mueeWBFZZ*>(), _1, sqrt_s_leptcoll_1400, pol_80, pol_0);
//    obsThFactory["mueeWBFZZ1400_m80_0"] = bind(boost::factory<mueeWBFZZ*>(), _1, sqrt_s_leptcoll_1400, -pol_80, pol_0);
    //
//    obsThFactory["mueeWBFZZ1500_p80_m30"] = bind(boost::factory<mueeWBFZZ*>(), _1, sqrt_s_leptcoll_1500, pol_80, -pol_30);
//    obsThFactory["mueeWBFZZ1500_m80_p30"] = bind(boost::factory<mueeWBFZZ*>(), _1, sqrt_s_leptcoll_1500, -pol_80, pol_30);
//    obsThFactory["mueeWBFZZ1500_p80_0"] = bind(boost::factory<mueeWBFZZ*>(), _1, sqrt_s_leptcoll_1500, pol_80, pol_0);
//    obsThFactory["mueeWBFZZ1500_m80_0"] = bind(boost::factory<mueeWBFZZ*>(), _1, sqrt_s_leptcoll_1500, -pol_80, pol_0);
    //
//    obsThFactory["mueeWBFZZ3000_p80_m30"] = bind(boost::factory<mueeWBFZZ*>(), _1, sqrt_s_leptcoll_3000, pol_80, -pol_30);
//    obsThFactory["mueeWBFZZ3000_m80_p30"] = bind(boost::factory<mueeWBFZZ*>(), _1, sqrt_s_leptcoll_3000, -pol_80, pol_30);
//    obsThFactory["mueeWBFZZ3000_p80_0"] = bind(boost::factory<mueeWBFZZ*>(), _1, sqrt_s_leptcoll_3000, pol_80, pol_0);
//    obsThFactory["mueeWBFZZ3000_m80_0"] = bind(boost::factory<mueeWBFZZ*>(), _1, sqrt_s_leptcoll_3000, -pol_80, pol_0);
    //
    obsThFactory["mueeWBFZga230"] = [=](const StandardModel& SM) { return new mueeWBFZga(SM, sqrt_s_leptcoll_230, 0., 0.); };
    obsThFactory["mueeWBFZga240"] = [=](const StandardModel& SM) { return new mueeWBFZga(SM, sqrt_s_leptcoll_240, 0., 0.); };
    obsThFactory["mueeWBFZga250"] = [=](const StandardModel& SM) { return new mueeWBFZga(SM, sqrt_s_leptcoll_250, 0., 0.); };
    obsThFactory["mueeWBFZga350"] = [=](const StandardModel& SM) { return new mueeWBFZga(SM, sqrt_s_leptcoll_350, 0., 0.); };
    obsThFactory["mueeWBFZga365"] = [=](const StandardModel& SM) { return new mueeWBFZga(SM, sqrt_s_leptcoll_365, 0., 0.); };
    obsThFactory["mueeWBFZga380"] = [=](const StandardModel& SM) { return new mueeWBFZga(SM, sqrt_s_leptcoll_380, 0., 0.); };
    obsThFactory["mueeWBFZga500"] = [=](const StandardModel& SM) { return new mueeWBFZga(SM, sqrt_s_leptcoll_500, 0., 0.); };
    obsThFactory["mueeWBFZga550"] = [=](const StandardModel& SM) { return new mueeWBFZga(SM, sqrt_s_leptcoll_550, 0., 0.); };
    obsThFactory["mueeWBFZga1000"] = [=](const StandardModel& SM) { return new mueeWBFZga(SM, sqrt_s_leptcoll_1000, 0., 0.); };
    obsThFactory["mueeWBFZga1400"] = [=](const StandardModel& SM) { return new mueeWBFZga(SM, sqrt_s_leptcoll_1400, 0., 0.); };
    obsThFactory["mueeWBFZga1500"] = [=](const StandardModel& SM) { return new mueeWBFZga(SM, sqrt_s_leptcoll_1500, 0., 0.); };
    obsThFactory["mueeWBFZga3000"] = [=](const StandardModel& SM) { return new mueeWBFZga(SM, sqrt_s_leptcoll_3000, 0., 0.); };
    //
//    obsThFactory["mueeWBFZga250_p80_m30"] = bind(boost::factory<mueeWBFZga*>(), _1, sqrt_s_leptcoll_250, pol_80, -pol_30);
//    obsThFactory["mueeWBFZga250_m80_p30"] = bind(boost::factory<mueeWBFZga*>(), _1, sqrt_s_leptcoll_250, -pol_80, pol_30);
//    obsThFactory["mueeWBFZga250_p80_0"] = bind(boost::factory<mueeWBFZga*>(), _1, sqrt_s_leptcoll_250, pol_80, pol_0);
//    obsThFactory["mueeWBFZga250_m80_0"] = bind(boost::factory<mueeWBFZga*>(), _1, sqrt_s_leptcoll_250, -pol_80, pol_0);
    //
//    obsThFactory["mueeWBFZga350_p80_m30"] = bind(boost::factory<mueeWBFZga*>(), _1, sqrt_s_leptcoll_350, pol_80, -pol_30);
//    obsThFactory["mueeWBFZga350_m80_p30"] = bind(boost::factory<mueeWBFZga*>(), _1, sqrt_s_leptcoll_350, -pol_80, pol_30);
//    obsThFactory["mueeWBFZga350_p80_0"] = bind(boost::factory<mueeWBFZga*>(), _1, sqrt_s_leptcoll_350, pol_80, pol_0);
//    obsThFactory["mueeWBFZga350_m80_0"] = bind(boost::factory<mueeWBFZga*>(), _1, sqrt_s_leptcoll_350, -pol_80, pol_0);
    //
//    obsThFactory["mueeWBFZga380_p80_m30"] = bind(boost::factory<mueeWBFZga*>(), _1, sqrt_s_leptcoll_380, pol_80, -pol_30);
//    obsThFactory["mueeWBFZga380_m80_p30"] = bind(boost::factory<mueeWBFZga*>(), _1, sqrt_s_leptcoll_380, -pol_80, pol_30);
//    obsThFactory["mueeWBFZga380_p80_0"] = bind(boost::factory<mueeWBFZga*>(), _1, sqrt_s_leptcoll_380, pol_80, pol_0);
//    obsThFactory["mueeWBFZga380_m80_0"] = bind(boost::factory<mueeWBFZga*>(), _1, sqrt_s_leptcoll_380, -pol_80, pol_0);
    //
//    obsThFactory["mueeWBFZga500_p80_m30"] = bind(boost::factory<mueeWBFZga*>(), _1, sqrt_s_leptcoll_500, pol_80, -pol_30);
//    obsThFactory["mueeWBFZga500_m80_p30"] = bind(boost::factory<mueeWBFZga*>(), _1, sqrt_s_leptcoll_500, -pol_80, pol_30);
//    obsThFactory["mueeWBFZga500_p80_0"] = bind(boost::factory<mueeWBFZga*>(), _1, sqrt_s_leptcoll_500, pol_80, pol_0);
//    obsThFactory["mueeWBFZga500_m80_0"] = bind(boost::factory<mueeWBFZga*>(), _1, sqrt_s_leptcoll_500, -pol_80, pol_0);
    //
//    obsThFactory["mueeWBFZga1000_p80_m30"] = bind(boost::factory<mueeWBFZga*>(), _1, sqrt_s_leptcoll_1000, pol_80, -pol_30);
//    obsThFactory["mueeWBFZga1000_m80_p30"] = bind(boost::factory<mueeWBFZga*>(), _1, sqrt_s_leptcoll_1000, -pol_80, pol_30);
//    obsThFactory["mueeWBFZga1000_p80_m20"] = bind(boost::factory<mueeWBFZga*>(), _1, sqrt_s_leptcoll_1000, pol_80, -pol_20);
//    obsThFactory["mueeWBFZga1000_m80_p20"] = bind(boost::factory<mueeWBFZga*>(), _1, sqrt_s_leptcoll_1000, -pol_80, pol_20);
//    obsThFactory["mueeWBFZga1000_p80_0"] = bind(boost::factory<mueeWBFZga*>(), _1, sqrt_s_leptcoll_1000, pol_80, pol_0);
//    obsThFactory["mueeWBFZga1000_m80_0"] = bind(boost::factory<mueeWBFZga*>(), _1, sqrt_s_leptcoll_1000, -pol_80, pol_0);
    //
//    obsThFactory["mueeWBFZga1400_p80_m30"] = bind(boost::factory<mueeWBFZga*>(), _1, sqrt_s_leptcoll_1400, pol_80, -pol_30);
//    obsThFactory["mueeWBFZga1400_m80_p30"] = bind(boost::factory<mueeWBFZga*>(), _1, sqrt_s_leptcoll_1400, -pol_80, pol_30);
//    obsThFactory["mueeWBFZga1400_p80_0"] = bind(boost::factory<mueeWBFZga*>(), _1, sqrt_s_leptcoll_1400, pol_80, pol_0);
//    obsThFactory["mueeWBFZga1400_m80_0"] = bind(boost::factory<mueeWBFZga*>(), _1, sqrt_s_leptcoll_1400, -pol_80, pol_0);
    //
//    obsThFactory["mueeWBFZga1500_p80_m30"] = bind(boost::factory<mueeWBFZga*>(), _1, sqrt_s_leptcoll_1500, pol_80, -pol_30);
//    obsThFactory["mueeWBFZga1500_m80_p30"] = bind(boost::factory<mueeWBFZga*>(), _1, sqrt_s_leptcoll_1500, -pol_80, pol_30);
//    obsThFactory["mueeWBFZga1500_p80_0"] = bind(boost::factory<mueeWBFZga*>(), _1, sqrt_s_leptcoll_1500, pol_80, pol_0);
//    obsThFactory["mueeWBFZga1500_m80_0"] = bind(boost::factory<mueeWBFZga*>(), _1, sqrt_s_leptcoll_1500, -pol_80, pol_0);
    //
//    obsThFactory["mueeWBFZga3000_p80_m30"] = bind(boost::factory<mueeWBFZga*>(), _1, sqrt_s_leptcoll_3000, pol_80, -pol_30);
//    obsThFactory["mueeWBFZga3000_m80_p30"] = bind(boost::factory<mueeWBFZga*>(), _1, sqrt_s_leptcoll_3000, -pol_80, pol_30);
//    obsThFactory["mueeWBFZga3000_p80_0"] = bind(boost::factory<mueeWBFZga*>(), _1, sqrt_s_leptcoll_3000, pol_80, pol_0);
//    obsThFactory["mueeWBFZga3000_m80_0"] = bind(boost::factory<mueeWBFZga*>(), _1, sqrt_s_leptcoll_3000, -pol_80, pol_0);
    //
    obsThFactory["mueeWBFgaga230"] = [=](const StandardModel& SM) { return new mueeWBFgaga(SM, sqrt_s_leptcoll_230, 0., 0.); };
    obsThFactory["mueeWBFgaga240"] = [=](const StandardModel& SM) { return new mueeWBFgaga(SM, sqrt_s_leptcoll_240, 0., 0.); };
    obsThFactory["mueeWBFgaga250"] = [=](const StandardModel& SM) { return new mueeWBFgaga(SM, sqrt_s_leptcoll_250, 0., 0.); };
    obsThFactory["mueeWBFgaga350"] = [=](const StandardModel& SM) { return new mueeWBFgaga(SM, sqrt_s_leptcoll_350, 0., 0.); };
    obsThFactory["mueeWBFgaga365"] = [=](const StandardModel& SM) { return new mueeWBFgaga(SM, sqrt_s_leptcoll_365, 0., 0.); };
    obsThFactory["mueeWBFgaga380"] = [=](const StandardModel& SM) { return new mueeWBFgaga(SM, sqrt_s_leptcoll_380, 0., 0.); };
    obsThFactory["mueeWBFgaga500"] = [=](const StandardModel& SM) { return new mueeWBFgaga(SM, sqrt_s_leptcoll_500, 0., 0.); };
    obsThFactory["mueeWBFgaga550"] = [=](const StandardModel& SM) { return new mueeWBFgaga(SM, sqrt_s_leptcoll_550, 0., 0.); };   
    obsThFactory["mueeWBFgaga1000"] = [=](const StandardModel& SM) { return new mueeWBFgaga(SM, sqrt_s_leptcoll_1000, 0., 0.); };
    obsThFactory["mueeWBFgaga1400"] = [=](const StandardModel& SM) { return new mueeWBFgaga(SM, sqrt_s_leptcoll_1400, 0., 0.); };
    obsThFactory["mueeWBFgaga1500"] = [=](const StandardModel& SM) { return new mueeWBFgaga(SM, sqrt_s_leptcoll_1500, 0., 0.); };
    obsThFactory["mueeWBFgaga3000"] = [=](const StandardModel& SM) { return new mueeWBFgaga(SM, sqrt_s_leptcoll_3000, 0., 0.); };
    //
//    obsThFactory["mueeWBFgaga250_p80_m30"] = bind(boost::factory<mueeWBFgaga*>(), _1, sqrt_s_leptcoll_250, pol_80, -pol_30);
//    obsThFactory["mueeWBFgaga250_m80_p30"] = bind(boost::factory<mueeWBFgaga*>(), _1, sqrt_s_leptcoll_250, -pol_80, pol_30);
//    obsThFactory["mueeWBFgaga250_p80_0"] = bind(boost::factory<mueeWBFgaga*>(), _1, sqrt_s_leptcoll_250, pol_80, pol_0);
//    obsThFactory["mueeWBFgaga250_m80_0"] = bind(boost::factory<mueeWBFgaga*>(), _1, sqrt_s_leptcoll_250, -pol_80, pol_0);
    //
//    obsThFactory["mueeWBFgaga350_p80_m30"] = bind(boost::factory<mueeWBFgaga*>(), _1, sqrt_s_leptcoll_350, pol_80, -pol_30);
//    obsThFactory["mueeWBFgaga350_m80_p30"] = bind(boost::factory<mueeWBFgaga*>(), _1, sqrt_s_leptcoll_350, -pol_80, pol_30);
//    obsThFactory["mueeWBFgaga350_p80_0"] = bind(boost::factory<mueeWBFgaga*>(), _1, sqrt_s_leptcoll_350, pol_80, pol_0);
//    obsThFactory["mueeWBFgaga350_m80_0"] = bind(boost::factory<mueeWBFgaga*>(), _1, sqrt_s_leptcoll_350, -pol_80, pol_0);
    //
//    obsThFactory["mueeWBFgaga380_p80_m30"] = bind(boost::factory<mueeWBFgaga*>(), _1, sqrt_s_leptcoll_380, pol_80, -pol_30);
//    obsThFactory["mueeWBFgaga380_m80_p30"] = bind(boost::factory<mueeWBFgaga*>(), _1, sqrt_s_leptcoll_380, -pol_80, pol_30);
//    obsThFactory["mueeWBFgaga380_p80_0"] = bind(boost::factory<mueeWBFgaga*>(), _1, sqrt_s_leptcoll_380, pol_80, pol_0);
//    obsThFactory["mueeWBFgaga380_m80_0"] = bind(boost::factory<mueeWBFgaga*>(), _1, sqrt_s_leptcoll_380, -pol_80, pol_0);
    //
//    obsThFactory["mueeWBFgaga500_p80_m30"] = bind(boost::factory<mueeWBFgaga*>(), _1, sqrt_s_leptcoll_500, pol_80, -pol_30);
//    obsThFactory["mueeWBFgaga500_m80_p30"] = bind(boost::factory<mueeWBFgaga*>(), _1, sqrt_s_leptcoll_500, -pol_80, pol_30);
//    obsThFactory["mueeWBFgaga500_p80_0"] = bind(boost::factory<mueeWBFgaga*>(), _1, sqrt_s_leptcoll_500, pol_80, pol_0);
//    obsThFactory["mueeWBFgaga500_m80_0"] = bind(boost::factory<mueeWBFgaga*>(), _1, sqrt_s_leptcoll_500, -pol_80, pol_0);
    //
//    obsThFactory["mueeWBFgaga1000_p80_m30"] = bind(boost::factory<mueeWBFgaga*>(), _1, sqrt_s_leptcoll_1000, pol_80, -pol_30);
//    obsThFactory["mueeWBFgaga1000_m80_p30"] = bind(boost::factory<mueeWBFgaga*>(), _1, sqrt_s_leptcoll_1000, -pol_80, pol_30);
//    obsThFactory["mueeWBFgaga1000_p80_m20"] = bind(boost::factory<mueeWBFgaga*>(), _1, sqrt_s_leptcoll_1000, pol_80, -pol_20);
//    obsThFactory["mueeWBFgaga1000_m80_p20"] = bind(boost::factory<mueeWBFgaga*>(), _1, sqrt_s_leptcoll_1000, -pol_80, pol_20);
//    obsThFactory["mueeWBFgaga1000_p80_0"] = bind(boost::factory<mueeWBFgaga*>(), _1, sqrt_s_leptcoll_1000, pol_80, pol_0);
//    obsThFactory["mueeWBFgaga1000_m80_0"] = bind(boost::factory<mueeWBFgaga*>(), _1, sqrt_s_leptcoll_1000, -pol_80, pol_0);
    //
//    obsThFactory["mueeWBFgaga1400_p80_m30"] = bind(boost::factory<mueeWBFgaga*>(), _1, sqrt_s_leptcoll_1400, pol_80, -pol_30);
//    obsThFactory["mueeWBFgaga1400_m80_p30"] = bind(boost::factory<mueeWBFgaga*>(), _1, sqrt_s_leptcoll_1400, -pol_80, pol_30);
//    obsThFactory["mueeWBFgaga1400_p80_0"] = bind(boost::factory<mueeWBFgaga*>(), _1, sqrt_s_leptcoll_1400, pol_80, pol_0);
//    obsThFactory["mueeWBFgaga1400_m80_0"] = bind(boost::factory<mueeWBFgaga*>(), _1, sqrt_s_leptcoll_1400, -pol_80, pol_0);
    //
//    obsThFactory["mueeWBFgaga1500_p80_m30"] = bind(boost::factory<mueeWBFgaga*>(), _1, sqrt_s_leptcoll_1500, pol_80, -pol_30);
//    obsThFactory["mueeWBFgaga1500_m80_p30"] = bind(boost::factory<mueeWBFgaga*>(), _1, sqrt_s_leptcoll_1500, -pol_80, pol_30);
//    obsThFactory["mueeWBFgaga1500_p80_0"] = bind(boost::factory<mueeWBFgaga*>(), _1, sqrt_s_leptcoll_1500, pol_80, pol_0);
//    obsThFactory["mueeWBFgaga1500_m80_0"] = bind(boost::factory<mueeWBFgaga*>(), _1, sqrt_s_leptcoll_1500, -pol_80, pol_0);
    //
//    obsThFactory["mueeWBFgaga3000_p80_m30"] = bind(boost::factory<mueeWBFgaga*>(), _1, sqrt_s_leptcoll_3000, pol_80, -pol_30);
//    obsThFactory["mueeWBFgaga3000_m80_p30"] = bind(boost::factory<mueeWBFgaga*>(), _1, sqrt_s_leptcoll_3000, -pol_80, pol_30);
//    obsThFactory["mueeWBFgaga3000_p80_0"] = bind(boost::factory<mueeWBFgaga*>(), _1, sqrt_s_leptcoll_3000, pol_80, pol_0);
//    obsThFactory["mueeWBFgaga3000_m80_0"] = bind(boost::factory<mueeWBFgaga*>(), _1, sqrt_s_leptcoll_3000, -pol_80, pol_0);
    //
    obsThFactory["mueeWBFmumu230"] = [=](const StandardModel& SM) { return new mueeWBFmumu(SM, sqrt_s_leptcoll_230, 0., 0.); };
    obsThFactory["mueeWBFmumu240"] = [=](const StandardModel& SM) { return new mueeWBFmumu(SM, sqrt_s_leptcoll_240, 0., 0.); };
    obsThFactory["mueeWBFmumu250"] = [=](const StandardModel& SM) { return new mueeWBFmumu(SM, sqrt_s_leptcoll_250, 0., 0.); };
    obsThFactory["mueeWBFmumu350"] = [=](const StandardModel& SM) { return new mueeWBFmumu(SM, sqrt_s_leptcoll_350, 0., 0.); };
    obsThFactory["mueeWBFmumu365"] = [=](const StandardModel& SM) { return new mueeWBFmumu(SM, sqrt_s_leptcoll_365, 0., 0.); };
    obsThFactory["mueeWBFmumu380"] = [=](const StandardModel& SM) { return new mueeWBFmumu(SM, sqrt_s_leptcoll_380, 0., 0.); };
    obsThFactory["mueeWBFmumu500"] = [=](const StandardModel& SM) { return new mueeWBFmumu(SM, sqrt_s_leptcoll_500, 0., 0.); };
    obsThFactory["mueeWBFmumu550"] = [=](const StandardModel& SM) { return new mueeWBFmumu(SM, sqrt_s_leptcoll_550, 0., 0.); };
    obsThFactory["mueeWBFmumu1000"] = [=](const StandardModel& SM) { return new mueeWBFmumu(SM, sqrt_s_leptcoll_1000, 0., 0.); };
    obsThFactory["mueeWBFmumu1400"] = [=](const StandardModel& SM) { return new mueeWBFmumu(SM, sqrt_s_leptcoll_1400, 0., 0.); };
    obsThFactory["mueeWBFmumu1500"] = [=](const StandardModel& SM) { return new mueeWBFmumu(SM, sqrt_s_leptcoll_1500, 0., 0.); };
    obsThFactory["mueeWBFmumu3000"] = [=](const StandardModel& SM) { return new mueeWBFmumu(SM, sqrt_s_leptcoll_3000, 0., 0.); };
    //
//    obsThFactory["mueeWBFmumu250_p80_m30"] = bind(boost::factory<mueeWBFmumu*>(), _1, sqrt_s_leptcoll_250, pol_80, -pol_30);
//    obsThFactory["mueeWBFmumu250_m80_p30"] = bind(boost::factory<mueeWBFmumu*>(), _1, sqrt_s_leptcoll_250, -pol_80, pol_30);
//    obsThFactory["mueeWBFmumu250_p80_0"] = bind(boost::factory<mueeWBFmumu*>(), _1, sqrt_s_leptcoll_250, pol_80, pol_0);
//    obsThFactory["mueeWBFmumu250_m80_0"] = bind(boost::factory<mueeWBFmumu*>(), _1, sqrt_s_leptcoll_250, -pol_80, pol_0);
    //
//    obsThFactory["mueeWBFmumu350_p80_m30"] = bind(boost::factory<mueeWBFmumu*>(), _1, sqrt_s_leptcoll_350, pol_80, -pol_30);
//    obsThFactory["mueeWBFmumu350_m80_p30"] = bind(boost::factory<mueeWBFmumu*>(), _1, sqrt_s_leptcoll_350, -pol_80, pol_30);
//    obsThFactory["mueeWBFmumu350_p80_0"] = bind(boost::factory<mueeWBFmumu*>(), _1, sqrt_s_leptcoll_350, pol_80, pol_0);
//    obsThFactory["mueeWBFmumu350_m80_0"] = bind(boost::factory<mueeWBFmumu*>(), _1, sqrt_s_leptcoll_350, -pol_80, pol_0);
    //
//    obsThFactory["mueeWBFmumu380_p80_m30"] = bind(boost::factory<mueeWBFmumu*>(), _1, sqrt_s_leptcoll_380, pol_80, -pol_30);
//    obsThFactory["mueeWBFmumu380_m80_p30"] = bind(boost::factory<mueeWBFmumu*>(), _1, sqrt_s_leptcoll_380, -pol_80, pol_30);
//    obsThFactory["mueeWBFmumu380_p80_0"] = bind(boost::factory<mueeWBFmumu*>(), _1, sqrt_s_leptcoll_380, pol_80, pol_0);
//    obsThFactory["mueeWBFmumu380_m80_0"] = bind(boost::factory<mueeWBFmumu*>(), _1, sqrt_s_leptcoll_380, -pol_80, pol_0);
    //
//    obsThFactory["mueeWBFmumu500_p80_m30"] = bind(boost::factory<mueeWBFmumu*>(), _1, sqrt_s_leptcoll_500, pol_80, -pol_30);
//    obsThFactory["mueeWBFmumu500_m80_p30"] = bind(boost::factory<mueeWBFmumu*>(), _1, sqrt_s_leptcoll_500, -pol_80, pol_30);
//    obsThFactory["mueeWBFmumu500_p80_0"] = bind(boost::factory<mueeWBFmumu*>(), _1, sqrt_s_leptcoll_500, pol_80, pol_0);
//    obsThFactory["mueeWBFmumu500_m80_0"] = bind(boost::factory<mueeWBFmumu*>(), _1, sqrt_s_leptcoll_500, -pol_80, pol_0);
    //
//    obsThFactory["mueeWBFmumu1000_p80_m30"] = bind(boost::factory<mueeWBFmumu*>(), _1, sqrt_s_leptcoll_1000, pol_80, -pol_30);
//    obsThFactory["mueeWBFmumu1000_m80_p30"] = bind(boost::factory<mueeWBFmumu*>(), _1, sqrt_s_leptcoll_1000, -pol_80, pol_30);
//    obsThFactory["mueeWBFmumu1000_p80_m20"] = bind(boost::factory<mueeWBFmumu*>(), _1, sqrt_s_leptcoll_1000, pol_80, -pol_20);
//    obsThFactory["mueeWBFmumu1000_m80_p20"] = bind(boost::factory<mueeWBFmumu*>(), _1, sqrt_s_leptcoll_1000, -pol_80, pol_20);
//    obsThFactory["mueeWBFmumu1000_p80_0"] = bind(boost::factory<mueeWBFmumu*>(), _1, sqrt_s_leptcoll_1000, pol_80, pol_0);
//    obsThFactory["mueeWBFmumu1000_m80_0"] = bind(boost::factory<mueeWBFmumu*>(), _1, sqrt_s_leptcoll_1000, -pol_80, pol_0);
    //
//    obsThFactory["mueeWBFmumu1400_p80_m30"] = bind(boost::factory<mueeWBFmumu*>(), _1, sqrt_s_leptcoll_1400, pol_80, -pol_30);
//    obsThFactory["mueeWBFmumu1400_m80_p30"] = bind(boost::factory<mueeWBFmumu*>(), _1, sqrt_s_leptcoll_1400, -pol_80, pol_30);
//    obsThFactory["mueeWBFmumu1400_p80_0"] = bind(boost::factory<mueeWBFmumu*>(), _1, sqrt_s_leptcoll_1400, pol_80, pol_0);
//    obsThFactory["mueeWBFmumu1400_m80_0"] = bind(boost::factory<mueeWBFmumu*>(), _1, sqrt_s_leptcoll_1400, -pol_80, pol_0);
    //
//    obsThFactory["mueeWBFmumu1500_p80_m30"] = bind(boost::factory<mueeWBFmumu*>(), _1, sqrt_s_leptcoll_1500, pol_80, -pol_30);
//    obsThFactory["mueeWBFmumu1500_m80_p30"] = bind(boost::factory<mueeWBFmumu*>(), _1, sqrt_s_leptcoll_1500, -pol_80, pol_30);
//    obsThFactory["mueeWBFmumu1500_p80_0"] = bind(boost::factory<mueeWBFmumu*>(), _1, sqrt_s_leptcoll_1500, pol_80, pol_0);
//    obsThFactory["mueeWBFmumu1500_m80_0"] = bind(boost::factory<mueeWBFmumu*>(), _1, sqrt_s_leptcoll_1500, -pol_80, pol_0);
    //
//    obsThFactory["mueeWBFmumu3000_p80_m30"] = bind(boost::factory<mueeWBFmumu*>(), _1, sqrt_s_leptcoll_3000, pol_80, -pol_30);
//    obsThFactory["mueeWBFmumu3000_m80_p30"] = bind(boost::factory<mueeWBFmumu*>(), _1, sqrt_s_leptcoll_3000, -pol_80, pol_30);
//    obsThFactory["mueeWBFmumu3000_p80_0"] = bind(boost::factory<mueeWBFmumu*>(), _1, sqrt_s_leptcoll_3000, pol_80, pol_0);
//    obsThFactory["mueeWBFmumu3000_m80_0"] = bind(boost::factory<mueeWBFmumu*>(), _1, sqrt_s_leptcoll_3000, -pol_80, pol_0);
    //
    // vvH
    // H -> bb
    obsThFactory["mueeHvvbb230"] = [=](const StandardModel& SM) { return new mueeHvvbb(SM, sqrt_s_leptcoll_230, 0., 0.); };
    obsThFactory["mueeHvvbb240"] = [=](const StandardModel& SM) { return new mueeHvvbb(SM, sqrt_s_leptcoll_240, 0., 0.); };
    obsThFactory["mueeHvvbb250"] = [=](const StandardModel& SM) { return new mueeHvvbb(SM, sqrt_s_leptcoll_250, 0., 0.); };
    obsThFactory["mueeHvvbb350"] = [=](const StandardModel& SM) { return new mueeHvvbb(SM, sqrt_s_leptcoll_350, 0., 0.); };
    obsThFactory["mueeHvvbb365"] = [=](const StandardModel& SM) { return new mueeHvvbb(SM, sqrt_s_leptcoll_365, 0., 0.); };
    obsThFactory["mueeHvvbb380"] = [=](const StandardModel& SM) { return new mueeHvvbb(SM, sqrt_s_leptcoll_380, 0., 0.); };
    obsThFactory["mueeHvvbb500"] = [=](const StandardModel& SM) { return new mueeHvvbb(SM, sqrt_s_leptcoll_500, 0., 0.); };
    obsThFactory["mueeHvvbb550"] = [=](const StandardModel& SM) { return new mueeHvvbb(SM, sqrt_s_leptcoll_550, 0., 0.); };
    obsThFactory["mueeHvvbb1000"] = [=](const StandardModel& SM) { return new mueeHvvbb(SM, sqrt_s_leptcoll_1000, 0., 0.); };
    obsThFactory["mueeHvvbb1400"] = [=](const StandardModel& SM) { return new mueeHvvbb(SM, sqrt_s_leptcoll_1400, 0., 0.); };
    obsThFactory["mueeHvvbb1500"] = [=](const StandardModel& SM) { return new mueeHvvbb(SM, sqrt_s_leptcoll_1500, 0., 0.); };
    obsThFactory["mueeHvvbb3000"] = [=](const StandardModel& SM) { return new mueeHvvbb(SM, sqrt_s_leptcoll_3000, 0., 0.); };
    //
    obsThFactory["mueeHvvbb250_p80_m30"] = [=](const StandardModel& SM) { return new mueeHvvbb(SM, sqrt_s_leptcoll_250, pol_80, -pol_30); };
    obsThFactory["mueeHvvbb250_m80_p30"] = [=](const StandardModel& SM) { return new mueeHvvbb(SM, sqrt_s_leptcoll_250, -pol_80, pol_30); };
    obsThFactory["mueeHvvbb250_p80_0"] = [=](const StandardModel& SM) { return new mueeHvvbb(SM, sqrt_s_leptcoll_250, pol_80, pol_0); };
    obsThFactory["mueeHvvbb250_m80_0"] = [=](const StandardModel& SM) { return new mueeHvvbb(SM, sqrt_s_leptcoll_250, -pol_80, pol_0); };
    //
    obsThFactory["mueeHvvbb350_p80_m30"] = [=](const StandardModel& SM) { return new mueeHvvbb(SM, sqrt_s_leptcoll_350, pol_80, -pol_30); };
    obsThFactory["mueeHvvbb350_m80_p30"] = [=](const StandardModel& SM) { return new mueeHvvbb(SM, sqrt_s_leptcoll_350, -pol_80, pol_30); };
    obsThFactory["mueeHvvbb350_p80_0"] = [=](const StandardModel& SM) { return new mueeHvvbb(SM, sqrt_s_leptcoll_350, pol_80, pol_0); };
    obsThFactory["mueeHvvbb350_m80_0"] = [=](const StandardModel& SM) { return new mueeHvvbb(SM, sqrt_s_leptcoll_350, -pol_80, pol_0); };
    //
    obsThFactory["mueeHvvbb380_p80_m30"] = [=](const StandardModel& SM) { return new mueeHvvbb(SM, sqrt_s_leptcoll_380, pol_80, -pol_30); };
    obsThFactory["mueeHvvbb380_m80_p30"] = [=](const StandardModel& SM) { return new mueeHvvbb(SM, sqrt_s_leptcoll_380, -pol_80, pol_30); };
    obsThFactory["mueeHvvbb380_p80_0"] = [=](const StandardModel& SM) { return new mueeHvvbb(SM, sqrt_s_leptcoll_380, pol_80, pol_0); };
    obsThFactory["mueeHvvbb380_m80_0"] = [=](const StandardModel& SM) { return new mueeHvvbb(SM, sqrt_s_leptcoll_380, -pol_80, pol_0); };
    //
    obsThFactory["mueeHvvbb500_p80_m30"] = [=](const StandardModel& SM) { return new mueeHvvbb(SM, sqrt_s_leptcoll_500, pol_80, -pol_30); };
    obsThFactory["mueeHvvbb500_m80_p30"] = [=](const StandardModel& SM) { return new mueeHvvbb(SM, sqrt_s_leptcoll_500, -pol_80, pol_30); };
    obsThFactory["mueeHvvbb500_p80_0"] = [=](const StandardModel& SM) { return new mueeHvvbb(SM, sqrt_s_leptcoll_500, pol_80, pol_0); };
    obsThFactory["mueeHvvbb500_m80_0"] = [=](const StandardModel& SM) { return new mueeHvvbb(SM, sqrt_s_leptcoll_500, -pol_80, pol_0); };
    //
    obsThFactory["mueeHvvbb550_p80_m30"] = [=](const StandardModel& SM) { return new mueeHvvbb(SM, sqrt_s_leptcoll_550, pol_80, -pol_30); };
    obsThFactory["mueeHvvbb550_m80_p30"] = [=](const StandardModel& SM) { return new mueeHvvbb(SM, sqrt_s_leptcoll_550, -pol_80, pol_30); };
    obsThFactory["mueeHvvbb550_p80_0"] = [=](const StandardModel& SM) { return new mueeHvvbb(SM, sqrt_s_leptcoll_550, pol_80, pol_0); };
    obsThFactory["mueeHvvbb550_m80_0"] = [=](const StandardModel& SM) { return new mueeHvvbb(SM, sqrt_s_leptcoll_550, -pol_80, pol_0); };
    //
    obsThFactory["mueeHvvbb1000_p80_m30"] = [=](const StandardModel& SM) { return new mueeHvvbb(SM, sqrt_s_leptcoll_1000, pol_80, -pol_30); };
    obsThFactory["mueeHvvbb1000_m80_p30"] = [=](const StandardModel& SM) { return new mueeHvvbb(SM, sqrt_s_leptcoll_1000, -pol_80, pol_30); };
    obsThFactory["mueeHvvbb1000_p80_m20"] = [=](const StandardModel& SM) { return new mueeHvvbb(SM, sqrt_s_leptcoll_1000, pol_80, -pol_20); };
    obsThFactory["mueeHvvbb1000_m80_p20"] = [=](const StandardModel& SM) { return new mueeHvvbb(SM, sqrt_s_leptcoll_1000, -pol_80, pol_20); };
    obsThFactory["mueeHvvbb1000_p80_0"] = [=](const StandardModel& SM) { return new mueeHvvbb(SM, sqrt_s_leptcoll_1000, pol_80, pol_0); };
    obsThFactory["mueeHvvbb1000_m80_0"] = [=](const StandardModel& SM) { return new mueeHvvbb(SM, sqrt_s_leptcoll_1000, -pol_80, pol_0); };
    //
    obsThFactory["mueeHvvbb1400_p80_m30"] = [=](const StandardModel& SM) { return new mueeHvvbb(SM, sqrt_s_leptcoll_1400, pol_80, -pol_30); };
    obsThFactory["mueeHvvbb1400_m80_p30"] = [=](const StandardModel& SM) { return new mueeHvvbb(SM, sqrt_s_leptcoll_1400, -pol_80, pol_30); };
    obsThFactory["mueeHvvbb1400_p80_0"] = [=](const StandardModel& SM) { return new mueeHvvbb(SM, sqrt_s_leptcoll_1400, pol_80, pol_0); };
    obsThFactory["mueeHvvbb1400_m80_0"] = [=](const StandardModel& SM) { return new mueeHvvbb(SM, sqrt_s_leptcoll_1400, -pol_80, pol_0); };
    //
    obsThFactory["mueeHvvbb1500_p80_m30"] = [=](const StandardModel& SM) { return new mueeHvvbb(SM, sqrt_s_leptcoll_1500, pol_80, -pol_30); };
    obsThFactory["mueeHvvbb1500_m80_p30"] = [=](const StandardModel& SM) { return new mueeHvvbb(SM, sqrt_s_leptcoll_1500, -pol_80, pol_30); };
    obsThFactory["mueeHvvbb1500_p80_0"] = [=](const StandardModel& SM) { return new mueeHvvbb(SM, sqrt_s_leptcoll_1500, pol_80, pol_0); };
    obsThFactory["mueeHvvbb1500_m80_0"] = [=](const StandardModel& SM) { return new mueeHvvbb(SM, sqrt_s_leptcoll_1500, -pol_80, pol_0); };
    //
    obsThFactory["mueeHvvbb3000_p80_m30"] = [=](const StandardModel& SM) { return new mueeHvvbb(SM, sqrt_s_leptcoll_3000, pol_80, -pol_30); };
    obsThFactory["mueeHvvbb3000_m80_p30"] = [=](const StandardModel& SM) { return new mueeHvvbb(SM, sqrt_s_leptcoll_3000, -pol_80, pol_30); };
    obsThFactory["mueeHvvbb3000_p80_0"] = [=](const StandardModel& SM) { return new mueeHvvbb(SM, sqrt_s_leptcoll_3000, pol_80, pol_0); };
    obsThFactory["mueeHvvbb3000_m80_0"] = [=](const StandardModel& SM) { return new mueeHvvbb(SM, sqrt_s_leptcoll_3000, -pol_80, pol_0); };
    //
    // H -> cc
    obsThFactory["mueeHvvcc230"] = [=](const StandardModel& SM) { return new mueeHvvcc(SM, sqrt_s_leptcoll_230, 0., 0.); };
    obsThFactory["mueeHvvcc240"] = [=](const StandardModel& SM) { return new mueeHvvcc(SM, sqrt_s_leptcoll_240, 0., 0.); };
    obsThFactory["mueeHvvcc250"] = [=](const StandardModel& SM) { return new mueeHvvcc(SM, sqrt_s_leptcoll_250, 0., 0.); };
    obsThFactory["mueeHvvcc350"] = [=](const StandardModel& SM) { return new mueeHvvcc(SM, sqrt_s_leptcoll_350, 0., 0.); };
    obsThFactory["mueeHvvcc365"] = [=](const StandardModel& SM) { return new mueeHvvcc(SM, sqrt_s_leptcoll_365, 0., 0.); };
    obsThFactory["mueeHvvcc380"] = [=](const StandardModel& SM) { return new mueeHvvcc(SM, sqrt_s_leptcoll_380, 0., 0.); };
    obsThFactory["mueeHvvcc500"] = [=](const StandardModel& SM) { return new mueeHvvcc(SM, sqrt_s_leptcoll_500, 0., 0.); };
    obsThFactory["mueeHvvcc550"] = [=](const StandardModel& SM) { return new mueeHvvcc(SM, sqrt_s_leptcoll_550, 0., 0.); };
    obsThFactory["mueeHvvcc1000"] = [=](const StandardModel& SM) { return new mueeHvvcc(SM, sqrt_s_leptcoll_1000, 0., 0.); };
    obsThFactory["mueeHvvcc1400"] = [=](const StandardModel& SM) { return new mueeHvvcc(SM, sqrt_s_leptcoll_1400, 0., 0.); };
    obsThFactory["mueeHvvcc1500"] = [=](const StandardModel& SM) { return new mueeHvvcc(SM, sqrt_s_leptcoll_1500, 0., 0.); };
    obsThFactory["mueeHvvcc3000"] = [=](const StandardModel& SM) { return new mueeHvvcc(SM, sqrt_s_leptcoll_3000, 0., 0.); };
    //
    obsThFactory["mueeHvvcc250_p80_m30"] = [=](const StandardModel& SM) { return new mueeHvvcc(SM, sqrt_s_leptcoll_250, pol_80, -pol_30); };
    obsThFactory["mueeHvvcc250_m80_p30"] = [=](const StandardModel& SM) { return new mueeHvvcc(SM, sqrt_s_leptcoll_250, -pol_80, pol_30); };
    obsThFactory["mueeHvvcc250_p80_0"] = [=](const StandardModel& SM) { return new mueeHvvcc(SM, sqrt_s_leptcoll_250, pol_80, pol_0); };
    obsThFactory["mueeHvvcc250_m80_0"] = [=](const StandardModel& SM) { return new mueeHvvcc(SM, sqrt_s_leptcoll_250, -pol_80, pol_0); };
    //
    obsThFactory["mueeHvvcc350_p80_m30"] = [=](const StandardModel& SM) { return new mueeHvvcc(SM, sqrt_s_leptcoll_350, pol_80, -pol_30); };
    obsThFactory["mueeHvvcc350_m80_p30"] = [=](const StandardModel& SM) { return new mueeHvvcc(SM, sqrt_s_leptcoll_350, -pol_80, pol_30); };
    obsThFactory["mueeHvvcc350_p80_0"] = [=](const StandardModel& SM) { return new mueeHvvcc(SM, sqrt_s_leptcoll_350, pol_80, pol_0); };
    obsThFactory["mueeHvvcc350_m80_0"] = [=](const StandardModel& SM) { return new mueeHvvcc(SM, sqrt_s_leptcoll_350, -pol_80, pol_0); };
    //
    obsThFactory["mueeHvvcc380_p80_m30"] = [=](const StandardModel& SM) { return new mueeHvvcc(SM, sqrt_s_leptcoll_380, pol_80, -pol_30); };
    obsThFactory["mueeHvvcc380_m80_p30"] = [=](const StandardModel& SM) { return new mueeHvvcc(SM, sqrt_s_leptcoll_380, -pol_80, pol_30); };
    obsThFactory["mueeHvvcc380_p80_0"] = [=](const StandardModel& SM) { return new mueeHvvcc(SM, sqrt_s_leptcoll_380, pol_80, pol_0); };
    obsThFactory["mueeHvvcc380_m80_0"] = [=](const StandardModel& SM) { return new mueeHvvcc(SM, sqrt_s_leptcoll_380, -pol_80, pol_0); };
    //
    obsThFactory["mueeHvvcc500_p80_m30"] = [=](const StandardModel& SM) { return new mueeHvvcc(SM, sqrt_s_leptcoll_500, pol_80, -pol_30); };
    obsThFactory["mueeHvvcc500_m80_p30"] = [=](const StandardModel& SM) { return new mueeHvvcc(SM, sqrt_s_leptcoll_500, -pol_80, pol_30); };
    obsThFactory["mueeHvvcc500_p80_0"] = [=](const StandardModel& SM) { return new mueeHvvcc(SM, sqrt_s_leptcoll_500, pol_80, pol_0); };
    obsThFactory["mueeHvvcc500_m80_0"] = [=](const StandardModel& SM) { return new mueeHvvcc(SM, sqrt_s_leptcoll_500, -pol_80, pol_0); };
    //
    obsThFactory["mueeHvvcc550_p80_m30"] = [=](const StandardModel& SM) { return new mueeHvvcc(SM, sqrt_s_leptcoll_550, pol_80, -pol_30); };
    obsThFactory["mueeHvvcc550_m80_p30"] = [=](const StandardModel& SM) { return new mueeHvvcc(SM, sqrt_s_leptcoll_550, -pol_80, pol_30); };
    obsThFactory["mueeHvvcc550_p80_0"] = [=](const StandardModel& SM) { return new mueeHvvcc(SM, sqrt_s_leptcoll_550, pol_80, pol_0); };
    obsThFactory["mueeHvvcc550_m80_0"] = [=](const StandardModel& SM) { return new mueeHvvcc(SM, sqrt_s_leptcoll_550, -pol_80, pol_0); };
    //
    obsThFactory["mueeHvvcc1000_p80_m30"] = [=](const StandardModel& SM) { return new mueeHvvcc(SM, sqrt_s_leptcoll_1000, pol_80, -pol_30); };
    obsThFactory["mueeHvvcc1000_m80_p30"] = [=](const StandardModel& SM) { return new mueeHvvcc(SM, sqrt_s_leptcoll_1000, -pol_80, pol_30); };
    obsThFactory["mueeHvvcc1000_p80_m20"] = [=](const StandardModel& SM) { return new mueeHvvcc(SM, sqrt_s_leptcoll_1000, pol_80, -pol_20); };
    obsThFactory["mueeHvvcc1000_m80_p20"] = [=](const StandardModel& SM) { return new mueeHvvcc(SM, sqrt_s_leptcoll_1000, -pol_80, pol_20); };
    obsThFactory["mueeHvvcc1000_p80_0"] = [=](const StandardModel& SM) { return new mueeHvvcc(SM, sqrt_s_leptcoll_1000, pol_80, pol_0); };
    obsThFactory["mueeHvvcc1000_m80_0"] = [=](const StandardModel& SM) { return new mueeHvvcc(SM, sqrt_s_leptcoll_1000, -pol_80, pol_0); };
    //
    obsThFactory["mueeHvvcc1400_p80_m30"] = [=](const StandardModel& SM) { return new mueeHvvcc(SM, sqrt_s_leptcoll_1400, pol_80, -pol_30); };
    obsThFactory["mueeHvvcc1400_m80_p30"] = [=](const StandardModel& SM) { return new mueeHvvcc(SM, sqrt_s_leptcoll_1400, -pol_80, pol_30); };
    obsThFactory["mueeHvvcc1400_p80_0"] = [=](const StandardModel& SM) { return new mueeHvvcc(SM, sqrt_s_leptcoll_1400, pol_80, pol_0); };
    obsThFactory["mueeHvvcc1400_m80_0"] = [=](const StandardModel& SM) { return new mueeHvvcc(SM, sqrt_s_leptcoll_1400, -pol_80, pol_0); };
    //
    obsThFactory["mueeHvvcc1500_p80_m30"] = [=](const StandardModel& SM) { return new mueeHvvcc(SM, sqrt_s_leptcoll_1500, pol_80, -pol_30); };
    obsThFactory["mueeHvvcc1500_m80_p30"] = [=](const StandardModel& SM) { return new mueeHvvcc(SM, sqrt_s_leptcoll_1500, -pol_80, pol_30); };
    obsThFactory["mueeHvvcc1500_p80_0"] = [=](const StandardModel& SM) { return new mueeHvvcc(SM, sqrt_s_leptcoll_1500, pol_80, pol_0); };
    obsThFactory["mueeHvvcc1500_m80_0"] = [=](const StandardModel& SM) { return new mueeHvvcc(SM, sqrt_s_leptcoll_1500, -pol_80, pol_0); };
    //
    obsThFactory["mueeHvvcc3000_p80_m30"] = [=](const StandardModel& SM) { return new mueeHvvcc(SM, sqrt_s_leptcoll_3000, pol_80, -pol_30); };
    obsThFactory["mueeHvvcc3000_m80_p30"] = [=](const StandardModel& SM) { return new mueeHvvcc(SM, sqrt_s_leptcoll_3000, -pol_80, pol_30); };
    obsThFactory["mueeHvvcc3000_p80_0"] = [=](const StandardModel& SM) { return new mueeHvvcc(SM, sqrt_s_leptcoll_3000, pol_80, pol_0); };
    obsThFactory["mueeHvvcc3000_m80_0"] = [=](const StandardModel& SM) { return new mueeHvvcc(SM, sqrt_s_leptcoll_3000, -pol_80, pol_0); };
    //
    // H -> ss
    obsThFactory["mueeHvvss230"] = [=](const StandardModel& SM) { return new mueeHvvss(SM, sqrt_s_leptcoll_230, 0., 0.); };
    obsThFactory["mueeHvvss240"] = [=](const StandardModel& SM) { return new mueeHvvss(SM, sqrt_s_leptcoll_240, 0., 0.); };
    obsThFactory["mueeHvvss250"] = [=](const StandardModel& SM) { return new mueeHvvss(SM, sqrt_s_leptcoll_250, 0., 0.); };
    obsThFactory["mueeHvvss350"] = [=](const StandardModel& SM) { return new mueeHvvss(SM, sqrt_s_leptcoll_350, 0., 0.); };
    obsThFactory["mueeHvvss365"] = [=](const StandardModel& SM) { return new mueeHvvss(SM, sqrt_s_leptcoll_365, 0., 0.); };
    obsThFactory["mueeHvvss380"] = [=](const StandardModel& SM) { return new mueeHvvss(SM, sqrt_s_leptcoll_380, 0., 0.); };
    obsThFactory["mueeHvvss500"] = [=](const StandardModel& SM) { return new mueeHvvss(SM, sqrt_s_leptcoll_500, 0., 0.); };
    obsThFactory["mueeHvvss550"] = [=](const StandardModel& SM) { return new mueeHvvss(SM, sqrt_s_leptcoll_550, 0., 0.); };
    obsThFactory["mueeHvvss1000"] = [=](const StandardModel& SM) { return new mueeHvvss(SM, sqrt_s_leptcoll_1000, 0., 0.); };
    obsThFactory["mueeHvvss1400"] = [=](const StandardModel& SM) { return new mueeHvvss(SM, sqrt_s_leptcoll_1400, 0., 0.); };
    obsThFactory["mueeHvvss1500"] = [=](const StandardModel& SM) { return new mueeHvvss(SM, sqrt_s_leptcoll_1500, 0., 0.); };
    obsThFactory["mueeHvvss3000"] = [=](const StandardModel& SM) { return new mueeHvvss(SM, sqrt_s_leptcoll_3000, 0., 0.); };
    //
    obsThFactory["mueeHvvss250_p80_m30"] = [=](const StandardModel& SM) { return new mueeHvvss(SM, sqrt_s_leptcoll_250, pol_80, -pol_30); };
    obsThFactory["mueeHvvss250_m80_p30"] = [=](const StandardModel& SM) { return new mueeHvvss(SM, sqrt_s_leptcoll_250, -pol_80, pol_30); };
    obsThFactory["mueeHvvss250_p80_0"] = [=](const StandardModel& SM) { return new mueeHvvss(SM, sqrt_s_leptcoll_250, pol_80, pol_0); };
    obsThFactory["mueeHvvss250_m80_0"] = [=](const StandardModel& SM) { return new mueeHvvss(SM, sqrt_s_leptcoll_250, -pol_80, pol_0); };
    //
    obsThFactory["mueeHvvss350_p80_m30"] = [=](const StandardModel& SM) { return new mueeHvvss(SM, sqrt_s_leptcoll_350, pol_80, -pol_30); };
    obsThFactory["mueeHvvss350_m80_p30"] = [=](const StandardModel& SM) { return new mueeHvvss(SM, sqrt_s_leptcoll_350, -pol_80, pol_30); };
    obsThFactory["mueeHvvss350_p80_0"] = [=](const StandardModel& SM) { return new mueeHvvss(SM, sqrt_s_leptcoll_350, pol_80, pol_0); };
    obsThFactory["mueeHvvss350_m80_0"] = [=](const StandardModel& SM) { return new mueeHvvss(SM, sqrt_s_leptcoll_350, -pol_80, pol_0); };
    //
    obsThFactory["mueeHvvss380_p80_m30"] = [=](const StandardModel& SM) { return new mueeHvvss(SM, sqrt_s_leptcoll_380, pol_80, -pol_30); };
    obsThFactory["mueeHvvss380_m80_p30"] = [=](const StandardModel& SM) { return new mueeHvvss(SM, sqrt_s_leptcoll_380, -pol_80, pol_30); };
    obsThFactory["mueeHvvss380_p80_0"] = [=](const StandardModel& SM) { return new mueeHvvss(SM, sqrt_s_leptcoll_380, pol_80, pol_0); };
    obsThFactory["mueeHvvss380_m80_0"] = [=](const StandardModel& SM) { return new mueeHvvss(SM, sqrt_s_leptcoll_380, -pol_80, pol_0); };
    //
    obsThFactory["mueeHvvss500_p80_m30"] = [=](const StandardModel& SM) { return new mueeHvvss(SM, sqrt_s_leptcoll_500, pol_80, -pol_30); };
    obsThFactory["mueeHvvss500_m80_p30"] = [=](const StandardModel& SM) { return new mueeHvvss(SM, sqrt_s_leptcoll_500, -pol_80, pol_30); };
    obsThFactory["mueeHvvss500_p80_0"] = [=](const StandardModel& SM) { return new mueeHvvss(SM, sqrt_s_leptcoll_500, pol_80, pol_0); };
    obsThFactory["mueeHvvss500_m80_0"] = [=](const StandardModel& SM) { return new mueeHvvss(SM, sqrt_s_leptcoll_500, -pol_80, pol_0); };
    //
    obsThFactory["mueeHvvss550_p80_m30"] = [=](const StandardModel& SM) { return new mueeHvvss(SM, sqrt_s_leptcoll_550, pol_80, -pol_30); };
    obsThFactory["mueeHvvss550_m80_p30"] = [=](const StandardModel& SM) { return new mueeHvvss(SM, sqrt_s_leptcoll_550, -pol_80, pol_30); };
    obsThFactory["mueeHvvss550_p80_0"] = [=](const StandardModel& SM) { return new mueeHvvss(SM, sqrt_s_leptcoll_550, pol_80, pol_0); };
    obsThFactory["mueeHvvss550_m80_0"] = [=](const StandardModel& SM) { return new mueeHvvss(SM, sqrt_s_leptcoll_550, -pol_80, pol_0); };
    //
    obsThFactory["mueeHvvss1000_p80_m30"] = [=](const StandardModel& SM) { return new mueeHvvss(SM, sqrt_s_leptcoll_1000, pol_80, -pol_30); };
    obsThFactory["mueeHvvss1000_m80_p30"] = [=](const StandardModel& SM) { return new mueeHvvss(SM, sqrt_s_leptcoll_1000, -pol_80, pol_30); };
    obsThFactory["mueeHvvss1000_p80_m20"] = [=](const StandardModel& SM) { return new mueeHvvss(SM, sqrt_s_leptcoll_1000, pol_80, -pol_20); };
    obsThFactory["mueeHvvss1000_m80_p20"] = [=](const StandardModel& SM) { return new mueeHvvss(SM, sqrt_s_leptcoll_1000, -pol_80, pol_20); };
    obsThFactory["mueeHvvss1000_p80_0"] = [=](const StandardModel& SM) { return new mueeHvvss(SM, sqrt_s_leptcoll_1000, pol_80, pol_0); };
    obsThFactory["mueeHvvss1000_m80_0"] = [=](const StandardModel& SM) { return new mueeHvvss(SM, sqrt_s_leptcoll_1000, -pol_80, pol_0); };
    //
    obsThFactory["mueeHvvss1400_p80_m30"] = [=](const StandardModel& SM) { return new mueeHvvss(SM, sqrt_s_leptcoll_1400, pol_80, -pol_30); };
    obsThFactory["mueeHvvss1400_m80_p30"] = [=](const StandardModel& SM) { return new mueeHvvss(SM, sqrt_s_leptcoll_1400, -pol_80, pol_30); };
    obsThFactory["mueeHvvss1400_p80_0"] = [=](const StandardModel& SM) { return new mueeHvvss(SM, sqrt_s_leptcoll_1400, pol_80, pol_0); };
    obsThFactory["mueeHvvss1400_m80_0"] = [=](const StandardModel& SM) { return new mueeHvvss(SM, sqrt_s_leptcoll_1400, -pol_80, pol_0); };
    //
    obsThFactory["mueeHvvss1500_p80_m30"] = [=](const StandardModel& SM) { return new mueeHvvss(SM, sqrt_s_leptcoll_1500, pol_80, -pol_30); };
    obsThFactory["mueeHvvss1500_m80_p30"] = [=](const StandardModel& SM) { return new mueeHvvss(SM, sqrt_s_leptcoll_1500, -pol_80, pol_30); };
    obsThFactory["mueeHvvss1500_p80_0"] = [=](const StandardModel& SM) { return new mueeHvvss(SM, sqrt_s_leptcoll_1500, pol_80, pol_0); };
    obsThFactory["mueeHvvss1500_m80_0"] = [=](const StandardModel& SM) { return new mueeHvvss(SM, sqrt_s_leptcoll_1500, -pol_80, pol_0); };
    //
    obsThFactory["mueeHvvss3000_p80_m30"] = [=](const StandardModel& SM) { return new mueeHvvss(SM, sqrt_s_leptcoll_3000, pol_80, -pol_30); };
    obsThFactory["mueeHvvss3000_m80_p30"] = [=](const StandardModel& SM) { return new mueeHvvss(SM, sqrt_s_leptcoll_3000, -pol_80, pol_30); };
    obsThFactory["mueeHvvss3000_p80_0"] = [=](const StandardModel& SM) { return new mueeHvvss(SM, sqrt_s_leptcoll_3000, pol_80, pol_0); };
    obsThFactory["mueeHvvss3000_m80_0"] = [=](const StandardModel& SM) { return new mueeHvvss(SM, sqrt_s_leptcoll_3000, -pol_80, pol_0); };    
    //
    // H -> gg
    obsThFactory["mueeHvvgg230"] = [=](const StandardModel& SM) { return new mueeHvvgg(SM, sqrt_s_leptcoll_230, 0., 0.); };
    obsThFactory["mueeHvvgg240"] = [=](const StandardModel& SM) { return new mueeHvvgg(SM, sqrt_s_leptcoll_240, 0., 0.); };
    obsThFactory["mueeHvvgg250"] = [=](const StandardModel& SM) { return new mueeHvvgg(SM, sqrt_s_leptcoll_250, 0., 0.); };
    obsThFactory["mueeHvvgg350"] = [=](const StandardModel& SM) { return new mueeHvvgg(SM, sqrt_s_leptcoll_350, 0., 0.); };
    obsThFactory["mueeHvvgg365"] = [=](const StandardModel& SM) { return new mueeHvvgg(SM, sqrt_s_leptcoll_365, 0., 0.); };
    obsThFactory["mueeHvvgg380"] = [=](const StandardModel& SM) { return new mueeHvvgg(SM, sqrt_s_leptcoll_380, 0., 0.); };
    obsThFactory["mueeHvvgg500"] = [=](const StandardModel& SM) { return new mueeHvvgg(SM, sqrt_s_leptcoll_500, 0., 0.); };
    obsThFactory["mueeHvvgg550"] = [=](const StandardModel& SM) { return new mueeHvvgg(SM, sqrt_s_leptcoll_550, 0., 0.); };
    obsThFactory["mueeHvvgg1000"] = [=](const StandardModel& SM) { return new mueeHvvgg(SM, sqrt_s_leptcoll_1000, 0., 0.); };
    obsThFactory["mueeHvvgg1400"] = [=](const StandardModel& SM) { return new mueeHvvgg(SM, sqrt_s_leptcoll_1400, 0., 0.); };
    obsThFactory["mueeHvvgg1500"] = [=](const StandardModel& SM) { return new mueeHvvgg(SM, sqrt_s_leptcoll_1500, 0., 0.); };
    obsThFactory["mueeHvvgg3000"] = [=](const StandardModel& SM) { return new mueeHvvgg(SM, sqrt_s_leptcoll_3000, 0., 0.); };
    //
    obsThFactory["mueeHvvgg250_p80_m30"] = [=](const StandardModel& SM) { return new mueeHvvgg(SM, sqrt_s_leptcoll_250, pol_80, -pol_30); };
    obsThFactory["mueeHvvgg250_m80_p30"] = [=](const StandardModel& SM) { return new mueeHvvgg(SM, sqrt_s_leptcoll_250, -pol_80, pol_30); };
    obsThFactory["mueeHvvgg250_p80_0"] = [=](const StandardModel& SM) { return new mueeHvvgg(SM, sqrt_s_leptcoll_250, pol_80, pol_0); };
    obsThFactory["mueeHvvgg250_m80_0"] = [=](const StandardModel& SM) { return new mueeHvvgg(SM, sqrt_s_leptcoll_250, -pol_80, pol_0); };
    //
    obsThFactory["mueeHvvgg350_p80_m30"] = [=](const StandardModel& SM) { return new mueeHvvgg(SM, sqrt_s_leptcoll_350, pol_80, -pol_30); };
    obsThFactory["mueeHvvgg350_m80_p30"] = [=](const StandardModel& SM) { return new mueeHvvgg(SM, sqrt_s_leptcoll_350, -pol_80, pol_30); };
    obsThFactory["mueeHvvgg350_p80_0"] = [=](const StandardModel& SM) { return new mueeHvvgg(SM, sqrt_s_leptcoll_350, pol_80, pol_0); };
    obsThFactory["mueeHvvgg350_m80_0"] = [=](const StandardModel& SM) { return new mueeHvvgg(SM, sqrt_s_leptcoll_350, -pol_80, pol_0); };
    //
    obsThFactory["mueeHvvgg380_p80_m30"] = [=](const StandardModel& SM) { return new mueeHvvgg(SM, sqrt_s_leptcoll_380, pol_80, -pol_30); };
    obsThFactory["mueeHvvgg380_m80_p30"] = [=](const StandardModel& SM) { return new mueeHvvgg(SM, sqrt_s_leptcoll_380, -pol_80, pol_30); };
    obsThFactory["mueeHvvgg380_p80_0"] = [=](const StandardModel& SM) { return new mueeHvvgg(SM, sqrt_s_leptcoll_380, pol_80, pol_0); };
    obsThFactory["mueeHvvgg380_m80_0"] = [=](const StandardModel& SM) { return new mueeHvvgg(SM, sqrt_s_leptcoll_380, -pol_80, pol_0); };
    //
    obsThFactory["mueeHvvgg500_p80_m30"] = [=](const StandardModel& SM) { return new mueeHvvgg(SM, sqrt_s_leptcoll_500, pol_80, -pol_30); };
    obsThFactory["mueeHvvgg500_m80_p30"] = [=](const StandardModel& SM) { return new mueeHvvgg(SM, sqrt_s_leptcoll_500, -pol_80, pol_30); };
    obsThFactory["mueeHvvgg500_p80_0"] = [=](const StandardModel& SM) { return new mueeHvvgg(SM, sqrt_s_leptcoll_500, pol_80, pol_0); };
    obsThFactory["mueeHvvgg500_m80_0"] = [=](const StandardModel& SM) { return new mueeHvvgg(SM, sqrt_s_leptcoll_500, -pol_80, pol_0); };
    //
    obsThFactory["mueeHvvgg550_p80_m30"] = [=](const StandardModel& SM) { return new mueeHvvgg(SM, sqrt_s_leptcoll_550, pol_80, -pol_30); };
    obsThFactory["mueeHvvgg550_m80_p30"] = [=](const StandardModel& SM) { return new mueeHvvgg(SM, sqrt_s_leptcoll_550, -pol_80, pol_30); };
    obsThFactory["mueeHvvgg550_p80_0"] = [=](const StandardModel& SM) { return new mueeHvvgg(SM, sqrt_s_leptcoll_550, pol_80, pol_0); };
    obsThFactory["mueeHvvgg550_m80_0"] = [=](const StandardModel& SM) { return new mueeHvvgg(SM, sqrt_s_leptcoll_550, -pol_80, pol_0); };
    //
    obsThFactory["mueeHvvgg1000_p80_m30"] = [=](const StandardModel& SM) { return new mueeHvvgg(SM, sqrt_s_leptcoll_1000, pol_80, -pol_30); };
    obsThFactory["mueeHvvgg1000_m80_p30"] = [=](const StandardModel& SM) { return new mueeHvvgg(SM, sqrt_s_leptcoll_1000, -pol_80, pol_30); };
    obsThFactory["mueeHvvgg1000_p80_m20"] = [=](const StandardModel& SM) { return new mueeHvvgg(SM, sqrt_s_leptcoll_1000, pol_80, -pol_20); };
    obsThFactory["mueeHvvgg1000_m80_p20"] = [=](const StandardModel& SM) { return new mueeHvvgg(SM, sqrt_s_leptcoll_1000, -pol_80, pol_20); };
    obsThFactory["mueeHvvgg1000_p80_0"] = [=](const StandardModel& SM) { return new mueeHvvgg(SM, sqrt_s_leptcoll_1000, pol_80, pol_0); };
    obsThFactory["mueeHvvgg1000_m80_0"] = [=](const StandardModel& SM) { return new mueeHvvgg(SM, sqrt_s_leptcoll_1000, -pol_80, pol_0); };
    //
    obsThFactory["mueeHvvgg1400_p80_m30"] = [=](const StandardModel& SM) { return new mueeHvvgg(SM, sqrt_s_leptcoll_1400, pol_80, -pol_30); };
    obsThFactory["mueeHvvgg1400_m80_p30"] = [=](const StandardModel& SM) { return new mueeHvvgg(SM, sqrt_s_leptcoll_1400, -pol_80, pol_30); };
    obsThFactory["mueeHvvgg1400_p80_0"] = [=](const StandardModel& SM) { return new mueeHvvgg(SM, sqrt_s_leptcoll_1400, pol_80, pol_0); };
    obsThFactory["mueeHvvgg1400_m80_0"] = [=](const StandardModel& SM) { return new mueeHvvgg(SM, sqrt_s_leptcoll_1400, -pol_80, pol_0); };
    //
    obsThFactory["mueeHvvgg1500_p80_m30"] = [=](const StandardModel& SM) { return new mueeHvvgg(SM, sqrt_s_leptcoll_1500, pol_80, -pol_30); };
    obsThFactory["mueeHvvgg1500_m80_p30"] = [=](const StandardModel& SM) { return new mueeHvvgg(SM, sqrt_s_leptcoll_1500, -pol_80, pol_30); };
    obsThFactory["mueeHvvgg1500_p80_0"] = [=](const StandardModel& SM) { return new mueeHvvgg(SM, sqrt_s_leptcoll_1500, pol_80, pol_0); };
    obsThFactory["mueeHvvgg1500_m80_0"] = [=](const StandardModel& SM) { return new mueeHvvgg(SM, sqrt_s_leptcoll_1500, -pol_80, pol_0); };
    //
    obsThFactory["mueeHvvgg3000_p80_m30"] = [=](const StandardModel& SM) { return new mueeHvvgg(SM, sqrt_s_leptcoll_3000, pol_80, -pol_30); };
    obsThFactory["mueeHvvgg3000_m80_p30"] = [=](const StandardModel& SM) { return new mueeHvvgg(SM, sqrt_s_leptcoll_3000, -pol_80, pol_30); };
    obsThFactory["mueeHvvgg3000_p80_0"] = [=](const StandardModel& SM) { return new mueeHvvgg(SM, sqrt_s_leptcoll_3000, pol_80, pol_0); };
    obsThFactory["mueeHvvgg3000_m80_0"] = [=](const StandardModel& SM) { return new mueeHvvgg(SM, sqrt_s_leptcoll_3000, -pol_80, pol_0); };
    //
    obsThFactory["mueeHvvWW230"] = [=](const StandardModel& SM) { return new mueeHvvWW(SM, sqrt_s_leptcoll_230, 0., 0.); };
    obsThFactory["mueeHvvWW240"] = [=](const StandardModel& SM) { return new mueeHvvWW(SM, sqrt_s_leptcoll_240, 0., 0.); };
    obsThFactory["mueeHvvWW250"] = [=](const StandardModel& SM) { return new mueeHvvWW(SM, sqrt_s_leptcoll_250, 0., 0.); };
    obsThFactory["mueeHvvWW350"] = [=](const StandardModel& SM) { return new mueeHvvWW(SM, sqrt_s_leptcoll_350, 0., 0.); };
    obsThFactory["mueeHvvWW365"] = [=](const StandardModel& SM) { return new mueeHvvWW(SM, sqrt_s_leptcoll_365, 0., 0.); };
    obsThFactory["mueeHvvWW380"] = [=](const StandardModel& SM) { return new mueeHvvWW(SM, sqrt_s_leptcoll_380, 0., 0.); };
    obsThFactory["mueeHvvWW500"] = [=](const StandardModel& SM) { return new mueeHvvWW(SM, sqrt_s_leptcoll_500, 0., 0.); };
    obsThFactory["mueeHvvWW550"] = [=](const StandardModel& SM) { return new mueeHvvWW(SM, sqrt_s_leptcoll_550, 0., 0.); };
    obsThFactory["mueeHvvWW1000"] = [=](const StandardModel& SM) { return new mueeHvvWW(SM, sqrt_s_leptcoll_1000, 0., 0.); };
    obsThFactory["mueeHvvWW1400"] = [=](const StandardModel& SM) { return new mueeHvvWW(SM, sqrt_s_leptcoll_1400, 0., 0.); };
    obsThFactory["mueeHvvWW1500"] = [=](const StandardModel& SM) { return new mueeHvvWW(SM, sqrt_s_leptcoll_1500, 0., 0.); };
    obsThFactory["mueeHvvWW3000"] = [=](const StandardModel& SM) { return new mueeHvvWW(SM, sqrt_s_leptcoll_3000, 0., 0.); };
    //
    obsThFactory["mueeHvvWW250_p80_m30"] = [=](const StandardModel& SM) { return new mueeHvvWW(SM, sqrt_s_leptcoll_250, pol_80, -pol_30); };
    obsThFactory["mueeHvvWW250_m80_p30"] = [=](const StandardModel& SM) { return new mueeHvvWW(SM, sqrt_s_leptcoll_250, -pol_80, pol_30); };
    obsThFactory["mueeHvvWW250_p80_0"] = [=](const StandardModel& SM) { return new mueeHvvWW(SM, sqrt_s_leptcoll_250, pol_80, pol_0); };
    obsThFactory["mueeHvvWW250_m80_0"] = [=](const StandardModel& SM) { return new mueeHvvWW(SM, sqrt_s_leptcoll_250, -pol_80, pol_0); };
    //
    obsThFactory["mueeHvvWW350_p80_m30"] = [=](const StandardModel& SM) { return new mueeHvvWW(SM, sqrt_s_leptcoll_350, pol_80, -pol_30); };
    obsThFactory["mueeHvvWW350_m80_p30"] = [=](const StandardModel& SM) { return new mueeHvvWW(SM, sqrt_s_leptcoll_350, -pol_80, pol_30); };
    obsThFactory["mueeHvvWW350_p80_0"] = [=](const StandardModel& SM) { return new mueeHvvWW(SM, sqrt_s_leptcoll_350, pol_80, pol_0); };
    obsThFactory["mueeHvvWW350_m80_0"] = [=](const StandardModel& SM) { return new mueeHvvWW(SM, sqrt_s_leptcoll_350, -pol_80, pol_0); };
    //
    obsThFactory["mueeHvvWW380_p80_m30"] = [=](const StandardModel& SM) { return new mueeHvvWW(SM, sqrt_s_leptcoll_380, pol_80, -pol_30); };
    obsThFactory["mueeHvvWW380_m80_p30"] = [=](const StandardModel& SM) { return new mueeHvvWW(SM, sqrt_s_leptcoll_380, -pol_80, pol_30); };
    obsThFactory["mueeHvvWW380_p80_0"] = [=](const StandardModel& SM) { return new mueeHvvWW(SM, sqrt_s_leptcoll_380, pol_80, pol_0); };
    obsThFactory["mueeHvvWW380_m80_0"] = [=](const StandardModel& SM) { return new mueeHvvWW(SM, sqrt_s_leptcoll_380, -pol_80, pol_0); };
    //
    obsThFactory["mueeHvvWW500_p80_m30"] = [=](const StandardModel& SM) { return new mueeHvvWW(SM, sqrt_s_leptcoll_500, pol_80, -pol_30); };
    obsThFactory["mueeHvvWW500_m80_p30"] = [=](const StandardModel& SM) { return new mueeHvvWW(SM, sqrt_s_leptcoll_500, -pol_80, pol_30); };
    obsThFactory["mueeHvvWW500_p80_0"] = [=](const StandardModel& SM) { return new mueeHvvWW(SM, sqrt_s_leptcoll_500, pol_80, pol_0); };
    obsThFactory["mueeHvvWW500_m80_0"] = [=](const StandardModel& SM) { return new mueeHvvWW(SM, sqrt_s_leptcoll_500, -pol_80, pol_0); };
    //
    obsThFactory["mueeHvvWW550_p80_m30"] = [=](const StandardModel& SM) { return new mueeHvvWW(SM, sqrt_s_leptcoll_550, pol_80, -pol_30); };
    obsThFactory["mueeHvvWW550_m80_p30"] = [=](const StandardModel& SM) { return new mueeHvvWW(SM, sqrt_s_leptcoll_550, -pol_80, pol_30); };
    obsThFactory["mueeHvvWW550_p80_0"] = [=](const StandardModel& SM) { return new mueeHvvWW(SM, sqrt_s_leptcoll_550, pol_80, pol_0); };
    obsThFactory["mueeHvvWW550_m80_0"] = [=](const StandardModel& SM) { return new mueeHvvWW(SM, sqrt_s_leptcoll_550, -pol_80, pol_0); };
    //
    obsThFactory["mueeHvvWW1000_p80_m30"] = [=](const StandardModel& SM) { return new mueeHvvWW(SM, sqrt_s_leptcoll_1000, pol_80, -pol_30); };
    obsThFactory["mueeHvvWW1000_m80_p30"] = [=](const StandardModel& SM) { return new mueeHvvWW(SM, sqrt_s_leptcoll_1000, -pol_80, pol_30); };
    obsThFactory["mueeHvvWW1000_p80_m20"] = [=](const StandardModel& SM) { return new mueeHvvWW(SM, sqrt_s_leptcoll_1000, pol_80, -pol_20); };
    obsThFactory["mueeHvvWW1000_m80_p20"] = [=](const StandardModel& SM) { return new mueeHvvWW(SM, sqrt_s_leptcoll_1000, -pol_80, pol_20); };
    obsThFactory["mueeHvvWW1000_p80_0"] = [=](const StandardModel& SM) { return new mueeHvvWW(SM, sqrt_s_leptcoll_1000, pol_80, pol_0); };
    obsThFactory["mueeHvvWW1000_m80_0"] = [=](const StandardModel& SM) { return new mueeHvvWW(SM, sqrt_s_leptcoll_1000, -pol_80, pol_0); };
    //
    obsThFactory["mueeHvvWW1400_p80_m30"] = [=](const StandardModel& SM) { return new mueeHvvWW(SM, sqrt_s_leptcoll_1400, pol_80, -pol_30); };
    obsThFactory["mueeHvvWW1400_m80_p30"] = [=](const StandardModel& SM) { return new mueeHvvWW(SM, sqrt_s_leptcoll_1400, -pol_80, pol_30); };
    obsThFactory["mueeHvvWW1400_p80_0"] = [=](const StandardModel& SM) { return new mueeHvvWW(SM, sqrt_s_leptcoll_1400, pol_80, pol_0); };
    obsThFactory["mueeHvvWW1400_m80_0"] = [=](const StandardModel& SM) { return new mueeHvvWW(SM, sqrt_s_leptcoll_1400, -pol_80, pol_0); };
    //
    obsThFactory["mueeHvvWW1500_p80_m30"] = [=](const StandardModel& SM) { return new mueeHvvWW(SM, sqrt_s_leptcoll_1500, pol_80, -pol_30); };
    obsThFactory["mueeHvvWW1500_m80_p30"] = [=](const StandardModel& SM) { return new mueeHvvWW(SM, sqrt_s_leptcoll_1500, -pol_80, pol_30); };
    obsThFactory["mueeHvvWW1500_p80_0"] = [=](const StandardModel& SM) { return new mueeHvvWW(SM, sqrt_s_leptcoll_1500, pol_80, pol_0); };
    obsThFactory["mueeHvvWW1500_m80_0"] = [=](const StandardModel& SM) { return new mueeHvvWW(SM, sqrt_s_leptcoll_1500, -pol_80, pol_0); };
    //
    obsThFactory["mueeHvvWW3000_p80_m30"] = [=](const StandardModel& SM) { return new mueeHvvWW(SM, sqrt_s_leptcoll_3000, pol_80, -pol_30); };
    obsThFactory["mueeHvvWW3000_m80_p30"] = [=](const StandardModel& SM) { return new mueeHvvWW(SM, sqrt_s_leptcoll_3000, -pol_80, pol_30); };
    obsThFactory["mueeHvvWW3000_p80_0"] = [=](const StandardModel& SM) { return new mueeHvvWW(SM, sqrt_s_leptcoll_3000, pol_80, pol_0); };
    obsThFactory["mueeHvvWW3000_m80_0"] = [=](const StandardModel& SM) { return new mueeHvvWW(SM, sqrt_s_leptcoll_3000, -pol_80, pol_0); };
    //
    obsThFactory["mueeHvvtautau230"] = [=](const StandardModel& SM) { return new mueeHvvtautau(SM, sqrt_s_leptcoll_230, 0., 0.); };
    obsThFactory["mueeHvvtautau240"] = [=](const StandardModel& SM) { return new mueeHvvtautau(SM, sqrt_s_leptcoll_240, 0., 0.); };
    obsThFactory["mueeHvvtautau250"] = [=](const StandardModel& SM) { return new mueeHvvtautau(SM, sqrt_s_leptcoll_250, 0., 0.); };
    obsThFactory["mueeHvvtautau350"] = [=](const StandardModel& SM) { return new mueeHvvtautau(SM, sqrt_s_leptcoll_350, 0., 0.); };
    obsThFactory["mueeHvvtautau365"] = [=](const StandardModel& SM) { return new mueeHvvtautau(SM, sqrt_s_leptcoll_365, 0., 0.); };
    obsThFactory["mueeHvvtautau380"] = [=](const StandardModel& SM) { return new mueeHvvtautau(SM, sqrt_s_leptcoll_380, 0., 0.); };
    obsThFactory["mueeHvvtautau500"] = [=](const StandardModel& SM) { return new mueeHvvtautau(SM, sqrt_s_leptcoll_500, 0., 0.); };
    obsThFactory["mueeHvvtautau550"] = [=](const StandardModel& SM) { return new mueeHvvtautau(SM, sqrt_s_leptcoll_550, 0., 0.); };
    obsThFactory["mueeHvvtautau1000"] = [=](const StandardModel& SM) { return new mueeHvvtautau(SM, sqrt_s_leptcoll_1000, 0., 0.); };
    obsThFactory["mueeHvvtautau1400"] = [=](const StandardModel& SM) { return new mueeHvvtautau(SM, sqrt_s_leptcoll_1400, 0., 0.); };
    obsThFactory["mueeHvvtautau1500"] = [=](const StandardModel& SM) { return new mueeHvvtautau(SM, sqrt_s_leptcoll_1500, 0., 0.); };
    obsThFactory["mueeHvvtautau3000"] = [=](const StandardModel& SM) { return new mueeHvvtautau(SM, sqrt_s_leptcoll_3000, 0., 0.); };
    //
    obsThFactory["mueeHvvtautau250_p80_m30"] = [=](const StandardModel& SM) { return new mueeHvvtautau(SM, sqrt_s_leptcoll_250, pol_80, -pol_30); };
    obsThFactory["mueeHvvtautau250_m80_p30"] = [=](const StandardModel& SM) { return new mueeHvvtautau(SM, sqrt_s_leptcoll_250, -pol_80, pol_30); };
    obsThFactory["mueeHvvtautau250_p80_0"] = [=](const StandardModel& SM) { return new mueeHvvtautau(SM, sqrt_s_leptcoll_250, pol_80, pol_0); };
    obsThFactory["mueeHvvtautau250_m80_0"] = [=](const StandardModel& SM) { return new mueeHvvtautau(SM, sqrt_s_leptcoll_250, -pol_80, pol_0); };
    //
    obsThFactory["mueeHvvtautau350_p80_m30"] = [=](const StandardModel& SM) { return new mueeHvvtautau(SM, sqrt_s_leptcoll_350, pol_80, -pol_30); };
    obsThFactory["mueeHvvtautau350_m80_p30"] = [=](const StandardModel& SM) { return new mueeHvvtautau(SM, sqrt_s_leptcoll_350, -pol_80, pol_30); };
    obsThFactory["mueeHvvtautau350_p80_0"] = [=](const StandardModel& SM) { return new mueeHvvtautau(SM, sqrt_s_leptcoll_350, pol_80, pol_0); };
    obsThFactory["mueeHvvtautau350_m80_0"] = [=](const StandardModel& SM) { return new mueeHvvtautau(SM, sqrt_s_leptcoll_350, -pol_80, pol_0); };
    //
    obsThFactory["mueeHvvtautau380_p80_m30"] = [=](const StandardModel& SM) { return new mueeHvvtautau(SM, sqrt_s_leptcoll_380, pol_80, -pol_30); };
    obsThFactory["mueeHvvtautau380_m80_p30"] = [=](const StandardModel& SM) { return new mueeHvvtautau(SM, sqrt_s_leptcoll_380, -pol_80, pol_30); };
    obsThFactory["mueeHvvtautau380_p80_0"] = [=](const StandardModel& SM) { return new mueeHvvtautau(SM, sqrt_s_leptcoll_380, pol_80, pol_0); };
    obsThFactory["mueeHvvtautau380_m80_0"] = [=](const StandardModel& SM) { return new mueeHvvtautau(SM, sqrt_s_leptcoll_380, -pol_80, pol_0); };
    //
    obsThFactory["mueeHvvtautau500_p80_m30"] = [=](const StandardModel& SM) { return new mueeHvvtautau(SM, sqrt_s_leptcoll_500, pol_80, -pol_30); };
    obsThFactory["mueeHvvtautau500_m80_p30"] = [=](const StandardModel& SM) { return new mueeHvvtautau(SM, sqrt_s_leptcoll_500, -pol_80, pol_30); };
    obsThFactory["mueeHvvtautau500_p80_0"] = [=](const StandardModel& SM) { return new mueeHvvtautau(SM, sqrt_s_leptcoll_500, pol_80, pol_0); };
    obsThFactory["mueeHvvtautau500_m80_0"] = [=](const StandardModel& SM) { return new mueeHvvtautau(SM, sqrt_s_leptcoll_500, -pol_80, pol_0); };
    //
    obsThFactory["mueeHvvtautau550_p80_m30"] = [=](const StandardModel& SM) { return new mueeHvvtautau(SM, sqrt_s_leptcoll_550, pol_80, -pol_30); };
    obsThFactory["mueeHvvtautau550_m80_p30"] = [=](const StandardModel& SM) { return new mueeHvvtautau(SM, sqrt_s_leptcoll_550, -pol_80, pol_30); };
    obsThFactory["mueeHvvtautau550_p80_0"] = [=](const StandardModel& SM) { return new mueeHvvtautau(SM, sqrt_s_leptcoll_550, pol_80, pol_0); };
    obsThFactory["mueeHvvtautau550_m80_0"] = [=](const StandardModel& SM) { return new mueeHvvtautau(SM, sqrt_s_leptcoll_550, -pol_80, pol_0); };
    //
    obsThFactory["mueeHvvtautau1000_p80_m30"] = [=](const StandardModel& SM) { return new mueeHvvtautau(SM, sqrt_s_leptcoll_1000, pol_80, -pol_30); };
    obsThFactory["mueeHvvtautau1000_m80_p30"] = [=](const StandardModel& SM) { return new mueeHvvtautau(SM, sqrt_s_leptcoll_1000, -pol_80, pol_30); };
    obsThFactory["mueeHvvtautau1000_p80_m20"] = [=](const StandardModel& SM) { return new mueeHvvtautau(SM, sqrt_s_leptcoll_1000, pol_80, -pol_20); };
    obsThFactory["mueeHvvtautau1000_m80_p20"] = [=](const StandardModel& SM) { return new mueeHvvtautau(SM, sqrt_s_leptcoll_1000, -pol_80, pol_20); };
    obsThFactory["mueeHvvtautau1000_p80_0"] = [=](const StandardModel& SM) { return new mueeHvvtautau(SM, sqrt_s_leptcoll_1000, pol_80, pol_0); };
    obsThFactory["mueeHvvtautau1000_m80_0"] = [=](const StandardModel& SM) { return new mueeHvvtautau(SM, sqrt_s_leptcoll_1000, -pol_80, pol_0); };
    //
    obsThFactory["mueeHvvtautau1400_p80_m30"] = [=](const StandardModel& SM) { return new mueeHvvtautau(SM, sqrt_s_leptcoll_1400, pol_80, -pol_30); };
    obsThFactory["mueeHvvtautau1400_m80_p30"] = [=](const StandardModel& SM) { return new mueeHvvtautau(SM, sqrt_s_leptcoll_1400, -pol_80, pol_30); };
    obsThFactory["mueeHvvtautau1400_p80_0"] = [=](const StandardModel& SM) { return new mueeHvvtautau(SM, sqrt_s_leptcoll_1400, pol_80, pol_0); };
    obsThFactory["mueeHvvtautau1400_m80_0"] = [=](const StandardModel& SM) { return new mueeHvvtautau(SM, sqrt_s_leptcoll_1400, -pol_80, pol_0); };
    //
    obsThFactory["mueeHvvtautau1500_p80_m30"] = [=](const StandardModel& SM) { return new mueeHvvtautau(SM, sqrt_s_leptcoll_1500, pol_80, -pol_30); };
    obsThFactory["mueeHvvtautau1500_m80_p30"] = [=](const StandardModel& SM) { return new mueeHvvtautau(SM, sqrt_s_leptcoll_1500, -pol_80, pol_30); };
    obsThFactory["mueeHvvtautau1500_p80_0"] = [=](const StandardModel& SM) { return new mueeHvvtautau(SM, sqrt_s_leptcoll_1500, pol_80, pol_0); };
    obsThFactory["mueeHvvtautau1500_m80_0"] = [=](const StandardModel& SM) { return new mueeHvvtautau(SM, sqrt_s_leptcoll_1500, -pol_80, pol_0); };
    //
    obsThFactory["mueeHvvtautau3000_p80_m30"] = [=](const StandardModel& SM) { return new mueeHvvtautau(SM, sqrt_s_leptcoll_3000, pol_80, -pol_30); };
    obsThFactory["mueeHvvtautau3000_m80_p30"] = [=](const StandardModel& SM) { return new mueeHvvtautau(SM, sqrt_s_leptcoll_3000, -pol_80, pol_30); };
    obsThFactory["mueeHvvtautau3000_p80_0"] = [=](const StandardModel& SM) { return new mueeHvvtautau(SM, sqrt_s_leptcoll_3000, pol_80, pol_0); };
    obsThFactory["mueeHvvtautau3000_m80_0"] = [=](const StandardModel& SM) { return new mueeHvvtautau(SM, sqrt_s_leptcoll_3000, -pol_80, pol_0); };
    //
    obsThFactory["mueeHvvZZ230"] = [=](const StandardModel& SM) { return new mueeHvvZZ(SM, sqrt_s_leptcoll_230, 0., 0.); };
    obsThFactory["mueeHvvZZ240"] = [=](const StandardModel& SM) { return new mueeHvvZZ(SM, sqrt_s_leptcoll_240, 0., 0.); };
    obsThFactory["mueeHvvZZ250"] = [=](const StandardModel& SM) { return new mueeHvvZZ(SM, sqrt_s_leptcoll_250, 0., 0.); };
    obsThFactory["mueeHvvZZ350"] = [=](const StandardModel& SM) { return new mueeHvvZZ(SM, sqrt_s_leptcoll_350, 0., 0.); };
    obsThFactory["mueeHvvZZ365"] = [=](const StandardModel& SM) { return new mueeHvvZZ(SM, sqrt_s_leptcoll_365, 0., 0.); };
    obsThFactory["mueeHvvZZ380"] = [=](const StandardModel& SM) { return new mueeHvvZZ(SM, sqrt_s_leptcoll_380, 0., 0.); };
    obsThFactory["mueeHvvZZ500"] = [=](const StandardModel& SM) { return new mueeHvvZZ(SM, sqrt_s_leptcoll_500, 0., 0.); };
    obsThFactory["mueeHvvZZ550"] = [=](const StandardModel& SM) { return new mueeHvvZZ(SM, sqrt_s_leptcoll_550, 0., 0.); };
    obsThFactory["mueeHvvZZ1000"] = [=](const StandardModel& SM) { return new mueeHvvZZ(SM, sqrt_s_leptcoll_1000, 0., 0.); };
    obsThFactory["mueeHvvZZ1400"] = [=](const StandardModel& SM) { return new mueeHvvZZ(SM, sqrt_s_leptcoll_1400, 0., 0.); };
    obsThFactory["mueeHvvZZ1500"] = [=](const StandardModel& SM) { return new mueeHvvZZ(SM, sqrt_s_leptcoll_1500, 0., 0.); };
    obsThFactory["mueeHvvZZ3000"] = [=](const StandardModel& SM) { return new mueeHvvZZ(SM, sqrt_s_leptcoll_3000, 0., 0.); };
    //
    obsThFactory["mueeHvvZZ250_p80_m30"] = [=](const StandardModel& SM) { return new mueeHvvZZ(SM, sqrt_s_leptcoll_250, pol_80, -pol_30); };
    obsThFactory["mueeHvvZZ250_m80_p30"] = [=](const StandardModel& SM) { return new mueeHvvZZ(SM, sqrt_s_leptcoll_250, -pol_80, pol_30); };
    obsThFactory["mueeHvvZZ250_p80_0"] = [=](const StandardModel& SM) { return new mueeHvvZZ(SM, sqrt_s_leptcoll_250, pol_80, pol_0); };
    obsThFactory["mueeHvvZZ250_m80_0"] = [=](const StandardModel& SM) { return new mueeHvvZZ(SM, sqrt_s_leptcoll_250, -pol_80, pol_0); };
    //
    obsThFactory["mueeHvvZZ350_p80_m30"] = [=](const StandardModel& SM) { return new mueeHvvZZ(SM, sqrt_s_leptcoll_350, pol_80, -pol_30); };
    obsThFactory["mueeHvvZZ350_m80_p30"] = [=](const StandardModel& SM) { return new mueeHvvZZ(SM, sqrt_s_leptcoll_350, -pol_80, pol_30); };
    obsThFactory["mueeHvvZZ350_p80_0"] = [=](const StandardModel& SM) { return new mueeHvvZZ(SM, sqrt_s_leptcoll_350, pol_80, pol_0); };
    obsThFactory["mueeHvvZZ350_m80_0"] = [=](const StandardModel& SM) { return new mueeHvvZZ(SM, sqrt_s_leptcoll_350, -pol_80, pol_0); };
    //
    obsThFactory["mueeHvvZZ380_p80_m30"] = [=](const StandardModel& SM) { return new mueeHvvZZ(SM, sqrt_s_leptcoll_380, pol_80, -pol_30); };
    obsThFactory["mueeHvvZZ380_m80_p30"] = [=](const StandardModel& SM) { return new mueeHvvZZ(SM, sqrt_s_leptcoll_380, -pol_80, pol_30); };
    obsThFactory["mueeHvvZZ380_p80_0"] = [=](const StandardModel& SM) { return new mueeHvvZZ(SM, sqrt_s_leptcoll_380, pol_80, pol_0); };
    obsThFactory["mueeHvvZZ380_m80_0"] = [=](const StandardModel& SM) { return new mueeHvvZZ(SM, sqrt_s_leptcoll_380, -pol_80, pol_0); };
    //
    obsThFactory["mueeHvvZZ500_p80_m30"] = [=](const StandardModel& SM) { return new mueeHvvZZ(SM, sqrt_s_leptcoll_500, pol_80, -pol_30); };
    obsThFactory["mueeHvvZZ500_m80_p30"] = [=](const StandardModel& SM) { return new mueeHvvZZ(SM, sqrt_s_leptcoll_500, -pol_80, pol_30); };
    obsThFactory["mueeHvvZZ500_p80_0"] = [=](const StandardModel& SM) { return new mueeHvvZZ(SM, sqrt_s_leptcoll_500, pol_80, pol_0); };
    obsThFactory["mueeHvvZZ500_m80_0"] = [=](const StandardModel& SM) { return new mueeHvvZZ(SM, sqrt_s_leptcoll_500, -pol_80, pol_0); };
    //
    obsThFactory["mueeHvvZZ550_p80_m30"] = [=](const StandardModel& SM) { return new mueeHvvZZ(SM, sqrt_s_leptcoll_550, pol_80, -pol_30); };
    obsThFactory["mueeHvvZZ550_m80_p30"] = [=](const StandardModel& SM) { return new mueeHvvZZ(SM, sqrt_s_leptcoll_550, -pol_80, pol_30); };
    obsThFactory["mueeHvvZZ550_p80_0"] = [=](const StandardModel& SM) { return new mueeHvvZZ(SM, sqrt_s_leptcoll_550, pol_80, pol_0); };
    obsThFactory["mueeHvvZZ550_m80_0"] = [=](const StandardModel& SM) { return new mueeHvvZZ(SM, sqrt_s_leptcoll_550, -pol_80, pol_0); };
    //
    obsThFactory["mueeHvvZZ1000_p80_m30"] = [=](const StandardModel& SM) { return new mueeHvvZZ(SM, sqrt_s_leptcoll_1000, pol_80, -pol_30); };
    obsThFactory["mueeHvvZZ1000_m80_p30"] = [=](const StandardModel& SM) { return new mueeHvvZZ(SM, sqrt_s_leptcoll_1000, -pol_80, pol_30); };
    obsThFactory["mueeHvvZZ1000_p80_m20"] = [=](const StandardModel& SM) { return new mueeHvvZZ(SM, sqrt_s_leptcoll_1000, pol_80, -pol_20); };
    obsThFactory["mueeHvvZZ1000_m80_p20"] = [=](const StandardModel& SM) { return new mueeHvvZZ(SM, sqrt_s_leptcoll_1000, -pol_80, pol_20); };
    obsThFactory["mueeHvvZZ1000_p80_0"] = [=](const StandardModel& SM) { return new mueeHvvZZ(SM, sqrt_s_leptcoll_1000, pol_80, pol_0); };
    obsThFactory["mueeHvvZZ1000_m80_0"] = [=](const StandardModel& SM) { return new mueeHvvZZ(SM, sqrt_s_leptcoll_1000, -pol_80, pol_0); };
    //
    obsThFactory["mueeHvvZZ1400_p80_m30"] = [=](const StandardModel& SM) { return new mueeHvvZZ(SM, sqrt_s_leptcoll_1400, pol_80, -pol_30); };
    obsThFactory["mueeHvvZZ1400_m80_p30"] = [=](const StandardModel& SM) { return new mueeHvvZZ(SM, sqrt_s_leptcoll_1400, -pol_80, pol_30); };
    obsThFactory["mueeHvvZZ1400_p80_0"] = [=](const StandardModel& SM) { return new mueeHvvZZ(SM, sqrt_s_leptcoll_1400, pol_80, pol_0); };
    obsThFactory["mueeHvvZZ1400_m80_0"] = [=](const StandardModel& SM) { return new mueeHvvZZ(SM, sqrt_s_leptcoll_1400, -pol_80, pol_0); };
    //
    obsThFactory["mueeHvvZZ1500_p80_m30"] = [=](const StandardModel& SM) { return new mueeHvvZZ(SM, sqrt_s_leptcoll_1500, pol_80, -pol_30); };
    obsThFactory["mueeHvvZZ1500_m80_p30"] = [=](const StandardModel& SM) { return new mueeHvvZZ(SM, sqrt_s_leptcoll_1500, -pol_80, pol_30); };
    obsThFactory["mueeHvvZZ1500_p80_0"] = [=](const StandardModel& SM) { return new mueeHvvZZ(SM, sqrt_s_leptcoll_1500, pol_80, pol_0); };
    obsThFactory["mueeHvvZZ1500_m80_0"] = [=](const StandardModel& SM) { return new mueeHvvZZ(SM, sqrt_s_leptcoll_1500, -pol_80, pol_0); };
    //
    obsThFactory["mueeHvvZZ3000_p80_m30"] = [=](const StandardModel& SM) { return new mueeHvvZZ(SM, sqrt_s_leptcoll_3000, pol_80, -pol_30); };
    obsThFactory["mueeHvvZZ3000_m80_p30"] = [=](const StandardModel& SM) { return new mueeHvvZZ(SM, sqrt_s_leptcoll_3000, -pol_80, pol_30); };
    obsThFactory["mueeHvvZZ3000_p80_0"] = [=](const StandardModel& SM) { return new mueeHvvZZ(SM, sqrt_s_leptcoll_3000, pol_80, pol_0); };
    obsThFactory["mueeHvvZZ3000_m80_0"] = [=](const StandardModel& SM) { return new mueeHvvZZ(SM, sqrt_s_leptcoll_3000, -pol_80, pol_0); };
    //
    obsThFactory["mueeHvvZga230"] = [=](const StandardModel& SM) { return new mueeHvvZga(SM, sqrt_s_leptcoll_230, 0., 0.); };
    obsThFactory["mueeHvvZga240"] = [=](const StandardModel& SM) { return new mueeHvvZga(SM, sqrt_s_leptcoll_240, 0., 0.); };
    obsThFactory["mueeHvvZga250"] = [=](const StandardModel& SM) { return new mueeHvvZga(SM, sqrt_s_leptcoll_250, 0., 0.); };
    obsThFactory["mueeHvvZga350"] = [=](const StandardModel& SM) { return new mueeHvvZga(SM, sqrt_s_leptcoll_350, 0., 0.); };
    obsThFactory["mueeHvvZga365"] = [=](const StandardModel& SM) { return new mueeHvvZga(SM, sqrt_s_leptcoll_365, 0., 0.); };
    obsThFactory["mueeHvvZga380"] = [=](const StandardModel& SM) { return new mueeHvvZga(SM, sqrt_s_leptcoll_380, 0., 0.); };
    obsThFactory["mueeHvvZga500"] = [=](const StandardModel& SM) { return new mueeHvvZga(SM, sqrt_s_leptcoll_500, 0., 0.); };
    obsThFactory["mueeHvvZga550"] = [=](const StandardModel& SM) { return new mueeHvvZga(SM, sqrt_s_leptcoll_550, 0., 0.); };
    obsThFactory["mueeHvvZga1000"] = [=](const StandardModel& SM) { return new mueeHvvZga(SM, sqrt_s_leptcoll_1000, 0., 0.); };
    obsThFactory["mueeHvvZga1400"] = [=](const StandardModel& SM) { return new mueeHvvZga(SM, sqrt_s_leptcoll_1400, 0., 0.); };
    obsThFactory["mueeHvvZga1500"] = [=](const StandardModel& SM) { return new mueeHvvZga(SM, sqrt_s_leptcoll_1500, 0., 0.); };
    obsThFactory["mueeHvvZga3000"] = [=](const StandardModel& SM) { return new mueeHvvZga(SM, sqrt_s_leptcoll_3000, 0., 0.); };
    //
    obsThFactory["mueeHvvZga250_p80_m30"] = [=](const StandardModel& SM) { return new mueeHvvZga(SM, sqrt_s_leptcoll_250, pol_80, -pol_30); };
    obsThFactory["mueeHvvZga250_m80_p30"] = [=](const StandardModel& SM) { return new mueeHvvZga(SM, sqrt_s_leptcoll_250, -pol_80, pol_30); };
    obsThFactory["mueeHvvZga250_p80_0"] = [=](const StandardModel& SM) { return new mueeHvvZga(SM, sqrt_s_leptcoll_250, pol_80, pol_0); };
    obsThFactory["mueeHvvZga250_m80_0"] = [=](const StandardModel& SM) { return new mueeHvvZga(SM, sqrt_s_leptcoll_250, -pol_80, pol_0); };
    //
    obsThFactory["mueeHvvZga350_p80_m30"] = [=](const StandardModel& SM) { return new mueeHvvZga(SM, sqrt_s_leptcoll_350, pol_80, -pol_30); };
    obsThFactory["mueeHvvZga350_m80_p30"] = [=](const StandardModel& SM) { return new mueeHvvZga(SM, sqrt_s_leptcoll_350, -pol_80, pol_30); };
    obsThFactory["mueeHvvZga350_p80_0"] = [=](const StandardModel& SM) { return new mueeHvvZga(SM, sqrt_s_leptcoll_350, pol_80, pol_0); };
    obsThFactory["mueeHvvZga350_m80_0"] = [=](const StandardModel& SM) { return new mueeHvvZga(SM, sqrt_s_leptcoll_350, -pol_80, pol_0); };
    //
    obsThFactory["mueeHvvZga380_p80_m30"] = [=](const StandardModel& SM) { return new mueeHvvZga(SM, sqrt_s_leptcoll_380, pol_80, -pol_30); };
    obsThFactory["mueeHvvZga380_m80_p30"] = [=](const StandardModel& SM) { return new mueeHvvZga(SM, sqrt_s_leptcoll_380, -pol_80, pol_30); };
    obsThFactory["mueeHvvZga380_p80_0"] = [=](const StandardModel& SM) { return new mueeHvvZga(SM, sqrt_s_leptcoll_380, pol_80, pol_0); };
    obsThFactory["mueeHvvZga380_m80_0"] = [=](const StandardModel& SM) { return new mueeHvvZga(SM, sqrt_s_leptcoll_380, -pol_80, pol_0); };
    //
    obsThFactory["mueeHvvZga500_p80_m30"] = [=](const StandardModel& SM) { return new mueeHvvZga(SM, sqrt_s_leptcoll_500, pol_80, -pol_30); };
    obsThFactory["mueeHvvZga500_m80_p30"] = [=](const StandardModel& SM) { return new mueeHvvZga(SM, sqrt_s_leptcoll_500, -pol_80, pol_30); };
    obsThFactory["mueeHvvZga500_p80_0"] = [=](const StandardModel& SM) { return new mueeHvvZga(SM, sqrt_s_leptcoll_500, pol_80, pol_0); };
    obsThFactory["mueeHvvZga500_m80_0"] = [=](const StandardModel& SM) { return new mueeHvvZga(SM, sqrt_s_leptcoll_500, -pol_80, pol_0); };
    //
    obsThFactory["mueeHvvZga550_p80_m30"] = [=](const StandardModel& SM) { return new mueeHvvZga(SM, sqrt_s_leptcoll_550, pol_80, -pol_30); };
    obsThFactory["mueeHvvZga550_m80_p30"] = [=](const StandardModel& SM) { return new mueeHvvZga(SM, sqrt_s_leptcoll_550, -pol_80, pol_30); };
    obsThFactory["mueeHvvZga550_p80_0"] = [=](const StandardModel& SM) { return new mueeHvvZga(SM, sqrt_s_leptcoll_550, pol_80, pol_0); };
    obsThFactory["mueeHvvZga550_m80_0"] = [=](const StandardModel& SM) { return new mueeHvvZga(SM, sqrt_s_leptcoll_550, -pol_80, pol_0); };
    //
    obsThFactory["mueeHvvZga1000_p80_m30"] = [=](const StandardModel& SM) { return new mueeHvvZga(SM, sqrt_s_leptcoll_1000, pol_80, -pol_30); };
    obsThFactory["mueeHvvZga1000_m80_p30"] = [=](const StandardModel& SM) { return new mueeHvvZga(SM, sqrt_s_leptcoll_1000, -pol_80, pol_30); };
    obsThFactory["mueeHvvZga1000_p80_m20"] = [=](const StandardModel& SM) { return new mueeHvvZga(SM, sqrt_s_leptcoll_1000, pol_80, -pol_20); };
    obsThFactory["mueeHvvZga1000_m80_p20"] = [=](const StandardModel& SM) { return new mueeHvvZga(SM, sqrt_s_leptcoll_1000, -pol_80, pol_20); };
    obsThFactory["mueeHvvZga1000_p80_0"] = [=](const StandardModel& SM) { return new mueeHvvZga(SM, sqrt_s_leptcoll_1000, pol_80, pol_0); };
    obsThFactory["mueeHvvZga1000_m80_0"] = [=](const StandardModel& SM) { return new mueeHvvZga(SM, sqrt_s_leptcoll_1000, -pol_80, pol_0); };
    //
    obsThFactory["mueeHvvZga1400_p80_m30"] = [=](const StandardModel& SM) { return new mueeHvvZga(SM, sqrt_s_leptcoll_1400, pol_80, -pol_30); };
    obsThFactory["mueeHvvZga1400_m80_p30"] = [=](const StandardModel& SM) { return new mueeHvvZga(SM, sqrt_s_leptcoll_1400, -pol_80, pol_30); };
    obsThFactory["mueeHvvZga1400_p80_0"] = [=](const StandardModel& SM) { return new mueeHvvZga(SM, sqrt_s_leptcoll_1400, pol_80, pol_0); };
    obsThFactory["mueeHvvZga1400_m80_0"] = [=](const StandardModel& SM) { return new mueeHvvZga(SM, sqrt_s_leptcoll_1400, -pol_80, pol_0); };
    //
    obsThFactory["mueeHvvZga1500_p80_m30"] = [=](const StandardModel& SM) { return new mueeHvvZga(SM, sqrt_s_leptcoll_1500, pol_80, -pol_30); };
    obsThFactory["mueeHvvZga1500_m80_p30"] = [=](const StandardModel& SM) { return new mueeHvvZga(SM, sqrt_s_leptcoll_1500, -pol_80, pol_30); };
    obsThFactory["mueeHvvZga1500_p80_0"] = [=](const StandardModel& SM) { return new mueeHvvZga(SM, sqrt_s_leptcoll_1500, pol_80, pol_0); };
    obsThFactory["mueeHvvZga1500_m80_0"] = [=](const StandardModel& SM) { return new mueeHvvZga(SM, sqrt_s_leptcoll_1500, -pol_80, pol_0); };
    //
    obsThFactory["mueeHvvZga3000_p80_m30"] = [=](const StandardModel& SM) { return new mueeHvvZga(SM, sqrt_s_leptcoll_3000, pol_80, -pol_30); };
    obsThFactory["mueeHvvZga3000_m80_p30"] = [=](const StandardModel& SM) { return new mueeHvvZga(SM, sqrt_s_leptcoll_3000, -pol_80, pol_30); };
    obsThFactory["mueeHvvZga3000_p80_0"] = [=](const StandardModel& SM) { return new mueeHvvZga(SM, sqrt_s_leptcoll_3000, pol_80, pol_0); };
    obsThFactory["mueeHvvZga3000_m80_0"] = [=](const StandardModel& SM) { return new mueeHvvZga(SM, sqrt_s_leptcoll_3000, -pol_80, pol_0); };
    //
    obsThFactory["mueeHvvgaga230"] = [=](const StandardModel& SM) { return new mueeHvvgaga(SM, sqrt_s_leptcoll_230, 0., 0.); };
    obsThFactory["mueeHvvgaga240"] = [=](const StandardModel& SM) { return new mueeHvvgaga(SM, sqrt_s_leptcoll_240, 0., 0.); };
    obsThFactory["mueeHvvgaga250"] = [=](const StandardModel& SM) { return new mueeHvvgaga(SM, sqrt_s_leptcoll_250, 0., 0.); };
    obsThFactory["mueeHvvgaga350"] = [=](const StandardModel& SM) { return new mueeHvvgaga(SM, sqrt_s_leptcoll_350, 0., 0.); };
    obsThFactory["mueeHvvgaga365"] = [=](const StandardModel& SM) { return new mueeHvvgaga(SM, sqrt_s_leptcoll_365, 0., 0.); };
    obsThFactory["mueeHvvgaga380"] = [=](const StandardModel& SM) { return new mueeHvvgaga(SM, sqrt_s_leptcoll_380, 0., 0.); };
    obsThFactory["mueeHvvgaga500"] = [=](const StandardModel& SM) { return new mueeHvvgaga(SM, sqrt_s_leptcoll_500, 0., 0.); };
    obsThFactory["mueeHvvgaga550"] = [=](const StandardModel& SM) { return new mueeHvvgaga(SM, sqrt_s_leptcoll_550, 0., 0.); };
    obsThFactory["mueeHvvgaga1000"] = [=](const StandardModel& SM) { return new mueeHvvgaga(SM, sqrt_s_leptcoll_1000, 0., 0.); };
    obsThFactory["mueeHvvgaga1400"] = [=](const StandardModel& SM) { return new mueeHvvgaga(SM, sqrt_s_leptcoll_1400, 0., 0.); };
    obsThFactory["mueeHvvgaga1500"] = [=](const StandardModel& SM) { return new mueeHvvgaga(SM, sqrt_s_leptcoll_1500, 0., 0.); };
    obsThFactory["mueeHvvgaga3000"] = [=](const StandardModel& SM) { return new mueeHvvgaga(SM, sqrt_s_leptcoll_3000, 0., 0.); };
    //
    obsThFactory["mueeHvvgaga250_p80_m30"] = [=](const StandardModel& SM) { return new mueeHvvgaga(SM, sqrt_s_leptcoll_250, pol_80, -pol_30); };
    obsThFactory["mueeHvvgaga250_m80_p30"] = [=](const StandardModel& SM) { return new mueeHvvgaga(SM, sqrt_s_leptcoll_250, -pol_80, pol_30); };
    obsThFactory["mueeHvvgaga250_p80_0"] = [=](const StandardModel& SM) { return new mueeHvvgaga(SM, sqrt_s_leptcoll_250, pol_80, pol_0); };
    obsThFactory["mueeHvvgaga250_m80_0"] = [=](const StandardModel& SM) { return new mueeHvvgaga(SM, sqrt_s_leptcoll_250, -pol_80, pol_0); };
    //
    obsThFactory["mueeHvvgaga350_p80_m30"] = [=](const StandardModel& SM) { return new mueeHvvgaga(SM, sqrt_s_leptcoll_350, pol_80, -pol_30); };
    obsThFactory["mueeHvvgaga350_m80_p30"] = [=](const StandardModel& SM) { return new mueeHvvgaga(SM, sqrt_s_leptcoll_350, -pol_80, pol_30); };
    obsThFactory["mueeHvvgaga350_p80_0"] = [=](const StandardModel& SM) { return new mueeHvvgaga(SM, sqrt_s_leptcoll_350, pol_80, pol_0); };
    obsThFactory["mueeHvvgaga350_m80_0"] = [=](const StandardModel& SM) { return new mueeHvvgaga(SM, sqrt_s_leptcoll_350, -pol_80, pol_0); };
    //
    obsThFactory["mueeHvvgaga380_p80_m30"] = [=](const StandardModel& SM) { return new mueeHvvgaga(SM, sqrt_s_leptcoll_380, pol_80, -pol_30); };
    obsThFactory["mueeHvvgaga380_m80_p30"] = [=](const StandardModel& SM) { return new mueeHvvgaga(SM, sqrt_s_leptcoll_380, -pol_80, pol_30); };
    obsThFactory["mueeHvvgaga380_p80_0"] = [=](const StandardModel& SM) { return new mueeHvvgaga(SM, sqrt_s_leptcoll_380, pol_80, pol_0); };
    obsThFactory["mueeHvvgaga380_m80_0"] = [=](const StandardModel& SM) { return new mueeHvvgaga(SM, sqrt_s_leptcoll_380, -pol_80, pol_0); };
    //
    obsThFactory["mueeHvvgaga500_p80_m30"] = [=](const StandardModel& SM) { return new mueeHvvgaga(SM, sqrt_s_leptcoll_500, pol_80, -pol_30); };
    obsThFactory["mueeHvvgaga500_m80_p30"] = [=](const StandardModel& SM) { return new mueeHvvgaga(SM, sqrt_s_leptcoll_500, -pol_80, pol_30); };
    obsThFactory["mueeHvvgaga500_p80_0"] = [=](const StandardModel& SM) { return new mueeHvvgaga(SM, sqrt_s_leptcoll_500, pol_80, pol_0); };
    obsThFactory["mueeHvvgaga500_m80_0"] = [=](const StandardModel& SM) { return new mueeHvvgaga(SM, sqrt_s_leptcoll_500, -pol_80, pol_0); };
    //
    obsThFactory["mueeHvvgaga550_p80_m30"] = [=](const StandardModel& SM) { return new mueeHvvgaga(SM, sqrt_s_leptcoll_550, pol_80, -pol_30); };
    obsThFactory["mueeHvvgaga550_m80_p30"] = [=](const StandardModel& SM) { return new mueeHvvgaga(SM, sqrt_s_leptcoll_550, -pol_80, pol_30); };
    obsThFactory["mueeHvvgaga550_p80_0"] = [=](const StandardModel& SM) { return new mueeHvvgaga(SM, sqrt_s_leptcoll_550, pol_80, pol_0); };
    obsThFactory["mueeHvvgaga550_m80_0"] = [=](const StandardModel& SM) { return new mueeHvvgaga(SM, sqrt_s_leptcoll_550, -pol_80, pol_0); };
    //
    obsThFactory["mueeHvvgaga1000_p80_m30"] = [=](const StandardModel& SM) { return new mueeHvvgaga(SM, sqrt_s_leptcoll_1000, pol_80, -pol_30); };
    obsThFactory["mueeHvvgaga1000_m80_p30"] = [=](const StandardModel& SM) { return new mueeHvvgaga(SM, sqrt_s_leptcoll_1000, -pol_80, pol_30); };
    obsThFactory["mueeHvvgaga1000_p80_m20"] = [=](const StandardModel& SM) { return new mueeHvvgaga(SM, sqrt_s_leptcoll_1000, pol_80, -pol_20); };
    obsThFactory["mueeHvvgaga1000_m80_p20"] = [=](const StandardModel& SM) { return new mueeHvvgaga(SM, sqrt_s_leptcoll_1000, -pol_80, pol_20); };
    obsThFactory["mueeHvvgaga1000_p80_0"] = [=](const StandardModel& SM) { return new mueeHvvgaga(SM, sqrt_s_leptcoll_1000, pol_80, pol_0); };
    obsThFactory["mueeHvvgaga1000_m80_0"] = [=](const StandardModel& SM) { return new mueeHvvgaga(SM, sqrt_s_leptcoll_1000, -pol_80, pol_0); };
    //
    obsThFactory["mueeHvvgaga1400_p80_m30"] = [=](const StandardModel& SM) { return new mueeHvvgaga(SM, sqrt_s_leptcoll_1400, pol_80, -pol_30); };
    obsThFactory["mueeHvvgaga1400_m80_p30"] = [=](const StandardModel& SM) { return new mueeHvvgaga(SM, sqrt_s_leptcoll_1400, -pol_80, pol_30); };
    obsThFactory["mueeHvvgaga1400_p80_0"] = [=](const StandardModel& SM) { return new mueeHvvgaga(SM, sqrt_s_leptcoll_1400, pol_80, pol_0); };
    obsThFactory["mueeHvvgaga1400_m80_0"] = [=](const StandardModel& SM) { return new mueeHvvgaga(SM, sqrt_s_leptcoll_1400, -pol_80, pol_0); };
    //
    obsThFactory["mueeHvvgaga1500_p80_m30"] = [=](const StandardModel& SM) { return new mueeHvvgaga(SM, sqrt_s_leptcoll_1500, pol_80, -pol_30); };
    obsThFactory["mueeHvvgaga1500_m80_p30"] = [=](const StandardModel& SM) { return new mueeHvvgaga(SM, sqrt_s_leptcoll_1500, -pol_80, pol_30); };
    obsThFactory["mueeHvvgaga1500_p80_0"] = [=](const StandardModel& SM) { return new mueeHvvgaga(SM, sqrt_s_leptcoll_1500, pol_80, pol_0); };
    obsThFactory["mueeHvvgaga1500_m80_0"] = [=](const StandardModel& SM) { return new mueeHvvgaga(SM, sqrt_s_leptcoll_1500, -pol_80, pol_0); };
    //
    obsThFactory["mueeHvvgaga3000_p80_m30"] = [=](const StandardModel& SM) { return new mueeHvvgaga(SM, sqrt_s_leptcoll_3000, pol_80, -pol_30); };
    obsThFactory["mueeHvvgaga3000_m80_p30"] = [=](const StandardModel& SM) { return new mueeHvvgaga(SM, sqrt_s_leptcoll_3000, -pol_80, pol_30); };
    obsThFactory["mueeHvvgaga3000_p80_0"] = [=](const StandardModel& SM) { return new mueeHvvgaga(SM, sqrt_s_leptcoll_3000, pol_80, pol_0); };
    obsThFactory["mueeHvvgaga3000_m80_0"] = [=](const StandardModel& SM) { return new mueeHvvgaga(SM, sqrt_s_leptcoll_3000, -pol_80, pol_0); };
    //
    obsThFactory["mueeHvvmumu230"] = [=](const StandardModel& SM) { return new mueeHvvmumu(SM, sqrt_s_leptcoll_230, 0., 0.); };
    obsThFactory["mueeHvvmumu240"] = [=](const StandardModel& SM) { return new mueeHvvmumu(SM, sqrt_s_leptcoll_240, 0., 0.); };
    obsThFactory["mueeHvvmumu250"] = [=](const StandardModel& SM) { return new mueeHvvmumu(SM, sqrt_s_leptcoll_250, 0., 0.); };
    obsThFactory["mueeHvvmumu350"] = [=](const StandardModel& SM) { return new mueeHvvmumu(SM, sqrt_s_leptcoll_350, 0., 0.); };
    obsThFactory["mueeHvvmumu365"] = [=](const StandardModel& SM) { return new mueeHvvmumu(SM, sqrt_s_leptcoll_365, 0., 0.); };
    obsThFactory["mueeHvvmumu380"] = [=](const StandardModel& SM) { return new mueeHvvmumu(SM, sqrt_s_leptcoll_380, 0., 0.); };
    obsThFactory["mueeHvvmumu500"] = [=](const StandardModel& SM) { return new mueeHvvmumu(SM, sqrt_s_leptcoll_500, 0., 0.); };
    obsThFactory["mueeHvvmumu550"] = [=](const StandardModel& SM) { return new mueeHvvmumu(SM, sqrt_s_leptcoll_550, 0., 0.); };
    obsThFactory["mueeHvvmumu1000"] = [=](const StandardModel& SM) { return new mueeHvvmumu(SM, sqrt_s_leptcoll_1000, 0., 0.); };
    obsThFactory["mueeHvvmumu1400"] = [=](const StandardModel& SM) { return new mueeHvvmumu(SM, sqrt_s_leptcoll_1400, 0., 0.); };
    obsThFactory["mueeHvvmumu1500"] = [=](const StandardModel& SM) { return new mueeHvvmumu(SM, sqrt_s_leptcoll_1500, 0., 0.); };
    obsThFactory["mueeHvvmumu3000"] = [=](const StandardModel& SM) { return new mueeHvvmumu(SM, sqrt_s_leptcoll_3000, 0., 0.); };
    //
    obsThFactory["mueeHvvmumu250_p80_m30"] = [=](const StandardModel& SM) { return new mueeHvvmumu(SM, sqrt_s_leptcoll_250, pol_80, -pol_30); };
    obsThFactory["mueeHvvmumu250_m80_p30"] = [=](const StandardModel& SM) { return new mueeHvvmumu(SM, sqrt_s_leptcoll_250, -pol_80, pol_30); };
    obsThFactory["mueeHvvmumu250_p80_0"] = [=](const StandardModel& SM) { return new mueeHvvmumu(SM, sqrt_s_leptcoll_250, pol_80, pol_0); };
    obsThFactory["mueeHvvmumu250_m80_0"] = [=](const StandardModel& SM) { return new mueeHvvmumu(SM, sqrt_s_leptcoll_250, -pol_80, pol_0); };
    //
    obsThFactory["mueeHvvmumu350_p80_m30"] = [=](const StandardModel& SM) { return new mueeHvvmumu(SM, sqrt_s_leptcoll_350, pol_80, -pol_30); };
    obsThFactory["mueeHvvmumu350_m80_p30"] = [=](const StandardModel& SM) { return new mueeHvvmumu(SM, sqrt_s_leptcoll_350, -pol_80, pol_30); };
    obsThFactory["mueeHvvmumu350_p80_0"] = [=](const StandardModel& SM) { return new mueeHvvmumu(SM, sqrt_s_leptcoll_350, pol_80, pol_0); };
    obsThFactory["mueeHvvmumu350_m80_0"] = [=](const StandardModel& SM) { return new mueeHvvmumu(SM, sqrt_s_leptcoll_350, -pol_80, pol_0); };
    //
    obsThFactory["mueeHvvmumu380_p80_m30"] = [=](const StandardModel& SM) { return new mueeHvvmumu(SM, sqrt_s_leptcoll_380, pol_80, -pol_30); };
    obsThFactory["mueeHvvmumu380_m80_p30"] = [=](const StandardModel& SM) { return new mueeHvvmumu(SM, sqrt_s_leptcoll_380, -pol_80, pol_30); };
    obsThFactory["mueeHvvmumu380_p80_0"] = [=](const StandardModel& SM) { return new mueeHvvmumu(SM, sqrt_s_leptcoll_380, pol_80, pol_0); };
    obsThFactory["mueeHvvmumu380_m80_0"] = [=](const StandardModel& SM) { return new mueeHvvmumu(SM, sqrt_s_leptcoll_380, -pol_80, pol_0); };
    //
    obsThFactory["mueeHvvmumu500_p80_m30"] = [=](const StandardModel& SM) { return new mueeHvvmumu(SM, sqrt_s_leptcoll_500, pol_80, -pol_30); };
    obsThFactory["mueeHvvmumu500_m80_p30"] = [=](const StandardModel& SM) { return new mueeHvvmumu(SM, sqrt_s_leptcoll_500, -pol_80, pol_30); };
    obsThFactory["mueeHvvmumu500_p80_0"] = [=](const StandardModel& SM) { return new mueeHvvmumu(SM, sqrt_s_leptcoll_500, pol_80, pol_0); };
    obsThFactory["mueeHvvmumu500_m80_0"] = [=](const StandardModel& SM) { return new mueeHvvmumu(SM, sqrt_s_leptcoll_500, -pol_80, pol_0); };
    //
    obsThFactory["mueeHvvmumu550_p80_m30"] = [=](const StandardModel& SM) { return new mueeHvvmumu(SM, sqrt_s_leptcoll_550, pol_80, -pol_30); };
    obsThFactory["mueeHvvmumu550_m80_p30"] = [=](const StandardModel& SM) { return new mueeHvvmumu(SM, sqrt_s_leptcoll_550, -pol_80, pol_30); };
    obsThFactory["mueeHvvmumu550_p80_0"] = [=](const StandardModel& SM) { return new mueeHvvmumu(SM, sqrt_s_leptcoll_550, pol_80, pol_0); };
    obsThFactory["mueeHvvmumu550_m80_0"] = [=](const StandardModel& SM) { return new mueeHvvmumu(SM, sqrt_s_leptcoll_550, -pol_80, pol_0); };
    //
    obsThFactory["mueeHvvmumu1000_p80_m30"] = [=](const StandardModel& SM) { return new mueeHvvmumu(SM, sqrt_s_leptcoll_1000, pol_80, -pol_30); };
    obsThFactory["mueeHvvmumu1000_m80_p30"] = [=](const StandardModel& SM) { return new mueeHvvmumu(SM, sqrt_s_leptcoll_1000, -pol_80, pol_30); };
    obsThFactory["mueeHvvmumu1000_p80_m20"] = [=](const StandardModel& SM) { return new mueeHvvmumu(SM, sqrt_s_leptcoll_1000, pol_80, -pol_20); };
    obsThFactory["mueeHvvmumu1000_m80_p20"] = [=](const StandardModel& SM) { return new mueeHvvmumu(SM, sqrt_s_leptcoll_1000, -pol_80, pol_20); };
    obsThFactory["mueeHvvmumu1000_p80_0"] = [=](const StandardModel& SM) { return new mueeHvvmumu(SM, sqrt_s_leptcoll_1000, pol_80, pol_0); };
    obsThFactory["mueeHvvmumu1000_m80_0"] = [=](const StandardModel& SM) { return new mueeHvvmumu(SM, sqrt_s_leptcoll_1000, -pol_80, pol_0); };
    //
    obsThFactory["mueeHvvmumu1400_p80_m30"] = [=](const StandardModel& SM) { return new mueeHvvmumu(SM, sqrt_s_leptcoll_1400, pol_80, -pol_30); };
    obsThFactory["mueeHvvmumu1400_m80_p30"] = [=](const StandardModel& SM) { return new mueeHvvmumu(SM, sqrt_s_leptcoll_1400, -pol_80, pol_30); };
    obsThFactory["mueeHvvmumu1400_p80_0"] = [=](const StandardModel& SM) { return new mueeHvvmumu(SM, sqrt_s_leptcoll_1400, pol_80, pol_0); };
    obsThFactory["mueeHvvmumu1400_m80_0"] = [=](const StandardModel& SM) { return new mueeHvvmumu(SM, sqrt_s_leptcoll_1400, -pol_80, pol_0); };
    //
    obsThFactory["mueeHvvmumu1500_p80_m30"] = [=](const StandardModel& SM) { return new mueeHvvmumu(SM, sqrt_s_leptcoll_1500, pol_80, -pol_30); };
    obsThFactory["mueeHvvmumu1500_m80_p30"] = [=](const StandardModel& SM) { return new mueeHvvmumu(SM, sqrt_s_leptcoll_1500, -pol_80, pol_30); };
    obsThFactory["mueeHvvmumu1500_p80_0"] = [=](const StandardModel& SM) { return new mueeHvvmumu(SM, sqrt_s_leptcoll_1500, pol_80, pol_0); };
    obsThFactory["mueeHvvmumu1500_m80_0"] = [=](const StandardModel& SM) { return new mueeHvvmumu(SM, sqrt_s_leptcoll_1500, -pol_80, pol_0); };
    //
    obsThFactory["mueeHvvmumu3000_p80_m30"] = [=](const StandardModel& SM) { return new mueeHvvmumu(SM, sqrt_s_leptcoll_3000, pol_80, -pol_30); };
    obsThFactory["mueeHvvmumu3000_m80_p30"] = [=](const StandardModel& SM) { return new mueeHvvmumu(SM, sqrt_s_leptcoll_3000, -pol_80, pol_30); };
    obsThFactory["mueeHvvmumu3000_p80_0"] = [=](const StandardModel& SM) { return new mueeHvvmumu(SM, sqrt_s_leptcoll_3000, pol_80, pol_0); };
    obsThFactory["mueeHvvmumu3000_m80_0"] = [=](const StandardModel& SM) { return new mueeHvvmumu(SM, sqrt_s_leptcoll_3000, -pol_80, pol_0); };
    //
    // ZH
    // H -> bb
    obsThFactory["mueeZHbb230"] = [=](const StandardModel& SM) { return new mueeZHbb(SM, sqrt_s_leptcoll_230, 0., 0.); };
    obsThFactory["mueeZHbb240"] = [=](const StandardModel& SM) { return new mueeZHbb(SM, sqrt_s_leptcoll_240, 0., 0.); };
    obsThFactory["mueeZHbb250"] = [=](const StandardModel& SM) { return new mueeZHbb(SM, sqrt_s_leptcoll_250, 0., 0.); };
    obsThFactory["mueeZHbb350"] = [=](const StandardModel& SM) { return new mueeZHbb(SM, sqrt_s_leptcoll_350, 0., 0.); };
    obsThFactory["mueeZHbb365"] = [=](const StandardModel& SM) { return new mueeZHbb(SM, sqrt_s_leptcoll_365, 0., 0.); };
    obsThFactory["mueeZHbb380"] = [=](const StandardModel& SM) { return new mueeZHbb(SM, sqrt_s_leptcoll_380, 0., 0.); };
    obsThFactory["mueeZHbb500"] = [=](const StandardModel& SM) { return new mueeZHbb(SM, sqrt_s_leptcoll_500, 0., 0.); };
    obsThFactory["mueeZHbb550"] = [=](const StandardModel& SM) { return new mueeZHbb(SM, sqrt_s_leptcoll_550, 0., 0.); };
    obsThFactory["mueeZHbb1000"] = [=](const StandardModel& SM) { return new mueeZHbb(SM, sqrt_s_leptcoll_1000, 0., 0.); };
    obsThFactory["mueeZHbb1400"] = [=](const StandardModel& SM) { return new mueeZHbb(SM, sqrt_s_leptcoll_1400, 0., 0.); };
    obsThFactory["mueeZHbb1500"] = [=](const StandardModel& SM) { return new mueeZHbb(SM, sqrt_s_leptcoll_1500, 0., 0.); };
    obsThFactory["mueeZHbb3000"] = [=](const StandardModel& SM) { return new mueeZHbb(SM, sqrt_s_leptcoll_3000, 0., 0.); };
    //
    obsThFactory["mueeZHbb250_p80_m30"] = [=](const StandardModel& SM) { return new mueeZHbb(SM, sqrt_s_leptcoll_250, pol_80, -pol_30); };
    obsThFactory["mueeZHbb250_m80_p30"] = [=](const StandardModel& SM) { return new mueeZHbb(SM, sqrt_s_leptcoll_250, -pol_80, pol_30); };
    obsThFactory["mueeZHbb250_p80_0"] = [=](const StandardModel& SM) { return new mueeZHbb(SM, sqrt_s_leptcoll_250, pol_80, pol_0); };
    obsThFactory["mueeZHbb250_m80_0"] = [=](const StandardModel& SM) { return new mueeZHbb(SM, sqrt_s_leptcoll_250, -pol_80, pol_0); };
    //
    obsThFactory["mueeZHbb350_p80_m30"] = [=](const StandardModel& SM) { return new mueeZHbb(SM, sqrt_s_leptcoll_350, pol_80, -pol_30); };
    obsThFactory["mueeZHbb350_m80_p30"] = [=](const StandardModel& SM) { return new mueeZHbb(SM, sqrt_s_leptcoll_350, -pol_80, pol_30); };
    obsThFactory["mueeZHbb350_p80_0"] = [=](const StandardModel& SM) { return new mueeZHbb(SM, sqrt_s_leptcoll_350, pol_80, pol_0); };
    obsThFactory["mueeZHbb350_m80_0"] = [=](const StandardModel& SM) { return new mueeZHbb(SM, sqrt_s_leptcoll_350, -pol_80, pol_0); };
    //
    obsThFactory["mueeZHbb380_p80_m30"] = [=](const StandardModel& SM) { return new mueeZHbb(SM, sqrt_s_leptcoll_380, pol_80, -pol_30); };
    obsThFactory["mueeZHbb380_m80_p30"] = [=](const StandardModel& SM) { return new mueeZHbb(SM, sqrt_s_leptcoll_380, -pol_80, pol_30); };
    obsThFactory["mueeZHbb380_p80_0"] = [=](const StandardModel& SM) { return new mueeZHbb(SM, sqrt_s_leptcoll_380, pol_80, pol_0); };
    obsThFactory["mueeZHbb380_m80_0"] = [=](const StandardModel& SM) { return new mueeZHbb(SM, sqrt_s_leptcoll_380, -pol_80, pol_0); };
    //
    obsThFactory["mueeZHbb500_p80_m30"] = [=](const StandardModel& SM) { return new mueeZHbb(SM, sqrt_s_leptcoll_500, pol_80, -pol_30); };
    obsThFactory["mueeZHbb500_m80_p30"] = [=](const StandardModel& SM) { return new mueeZHbb(SM, sqrt_s_leptcoll_500, -pol_80, pol_30); };
    obsThFactory["mueeZHbb500_p80_0"] = [=](const StandardModel& SM) { return new mueeZHbb(SM, sqrt_s_leptcoll_500, pol_80, pol_0); };
    obsThFactory["mueeZHbb500_m80_0"] = [=](const StandardModel& SM) { return new mueeZHbb(SM, sqrt_s_leptcoll_500, -pol_80, pol_0); };
    //
    obsThFactory["mueeZHbb550_p80_m30"] = [=](const StandardModel& SM) { return new mueeZHbb(SM, sqrt_s_leptcoll_550, pol_80, -pol_30); };
    obsThFactory["mueeZHbb550_m80_p30"] = [=](const StandardModel& SM) { return new mueeZHbb(SM, sqrt_s_leptcoll_550, -pol_80, pol_30); };
    obsThFactory["mueeZHbb550_p80_0"] = [=](const StandardModel& SM) { return new mueeZHbb(SM, sqrt_s_leptcoll_550, pol_80, pol_0); };
    obsThFactory["mueeZHbb550_m80_0"] = [=](const StandardModel& SM) { return new mueeZHbb(SM, sqrt_s_leptcoll_550, -pol_80, pol_0); };
    //
    obsThFactory["mueeZHbb1000_p80_m30"] = [=](const StandardModel& SM) { return new mueeZHbb(SM, sqrt_s_leptcoll_1000, pol_80, -pol_30); };
    obsThFactory["mueeZHbb1000_m80_p30"] = [=](const StandardModel& SM) { return new mueeZHbb(SM, sqrt_s_leptcoll_1000, -pol_80, pol_30); };
    obsThFactory["mueeZHbb1000_p80_m20"] = [=](const StandardModel& SM) { return new mueeZHbb(SM, sqrt_s_leptcoll_1000, pol_80, -pol_20); };
    obsThFactory["mueeZHbb1000_m80_p20"] = [=](const StandardModel& SM) { return new mueeZHbb(SM, sqrt_s_leptcoll_1000, -pol_80, pol_20); };
    obsThFactory["mueeZHbb1000_p80_0"] = [=](const StandardModel& SM) { return new mueeZHbb(SM, sqrt_s_leptcoll_1000, pol_80, pol_0); };
    obsThFactory["mueeZHbb1000_m80_0"] = [=](const StandardModel& SM) { return new mueeZHbb(SM, sqrt_s_leptcoll_1000, -pol_80, pol_0); };
    //
    obsThFactory["mueeZHbb1400_p80_m30"] = [=](const StandardModel& SM) { return new mueeZHbb(SM, sqrt_s_leptcoll_1400, pol_80, -pol_30); };
    obsThFactory["mueeZHbb1400_m80_p30"] = [=](const StandardModel& SM) { return new mueeZHbb(SM, sqrt_s_leptcoll_1400, -pol_80, pol_30); };
    obsThFactory["mueeZHbb1400_p80_0"] = [=](const StandardModel& SM) { return new mueeZHbb(SM, sqrt_s_leptcoll_1400, pol_80, pol_0); };
    obsThFactory["mueeZHbb1400_m80_0"] = [=](const StandardModel& SM) { return new mueeZHbb(SM, sqrt_s_leptcoll_1400, -pol_80, pol_0); };
    //
    obsThFactory["mueeZHbb1500_p80_m30"] = [=](const StandardModel& SM) { return new mueeZHbb(SM, sqrt_s_leptcoll_1500, pol_80, -pol_30); };
    obsThFactory["mueeZHbb1500_m80_p30"] = [=](const StandardModel& SM) { return new mueeZHbb(SM, sqrt_s_leptcoll_1500, -pol_80, pol_30); };
    obsThFactory["mueeZHbb1500_p80_0"] = [=](const StandardModel& SM) { return new mueeZHbb(SM, sqrt_s_leptcoll_1500, pol_80, pol_0); };
    obsThFactory["mueeZHbb1500_m80_0"] = [=](const StandardModel& SM) { return new mueeZHbb(SM, sqrt_s_leptcoll_1500, -pol_80, pol_0); };
    //
    obsThFactory["mueeZHbb3000_p80_m30"] = [=](const StandardModel& SM) { return new mueeZHbb(SM, sqrt_s_leptcoll_3000, pol_80, -pol_30); };
    obsThFactory["mueeZHbb3000_m80_p30"] = [=](const StandardModel& SM) { return new mueeZHbb(SM, sqrt_s_leptcoll_3000, -pol_80, pol_30); };
    obsThFactory["mueeZHbb3000_p80_0"] = [=](const StandardModel& SM) { return new mueeZHbb(SM, sqrt_s_leptcoll_3000, pol_80, pol_0); };
    obsThFactory["mueeZHbb3000_m80_0"] = [=](const StandardModel& SM) { return new mueeZHbb(SM, sqrt_s_leptcoll_3000, -pol_80, pol_0); };
    //
    // H -> cc
    obsThFactory["mueeZHcc230"] = [=](const StandardModel& SM) { return new mueeZHcc(SM, sqrt_s_leptcoll_230, 0., 0.); };
    obsThFactory["mueeZHcc240"] = [=](const StandardModel& SM) { return new mueeZHcc(SM, sqrt_s_leptcoll_240, 0., 0.); };
    obsThFactory["mueeZHcc250"] = [=](const StandardModel& SM) { return new mueeZHcc(SM, sqrt_s_leptcoll_250, 0., 0.); };
    obsThFactory["mueeZHcc350"] = [=](const StandardModel& SM) { return new mueeZHcc(SM, sqrt_s_leptcoll_350, 0., 0.); };
    obsThFactory["mueeZHcc365"] = [=](const StandardModel& SM) { return new mueeZHcc(SM, sqrt_s_leptcoll_365, 0., 0.); };
    obsThFactory["mueeZHcc380"] = [=](const StandardModel& SM) { return new mueeZHcc(SM, sqrt_s_leptcoll_380, 0., 0.); };
    obsThFactory["mueeZHcc500"] = [=](const StandardModel& SM) { return new mueeZHcc(SM, sqrt_s_leptcoll_500, 0., 0.); };
    obsThFactory["mueeZHcc550"] = [=](const StandardModel& SM) { return new mueeZHcc(SM, sqrt_s_leptcoll_550, 0., 0.); };
    obsThFactory["mueeZHcc1000"] = [=](const StandardModel& SM) { return new mueeZHcc(SM, sqrt_s_leptcoll_1000, 0., 0.); };
    obsThFactory["mueeZHcc1400"] = [=](const StandardModel& SM) { return new mueeZHcc(SM, sqrt_s_leptcoll_1400, 0., 0.); };
    obsThFactory["mueeZHcc1500"] = [=](const StandardModel& SM) { return new mueeZHcc(SM, sqrt_s_leptcoll_1500, 0., 0.); };
    obsThFactory["mueeZHcc3000"] = [=](const StandardModel& SM) { return new mueeZHcc(SM, sqrt_s_leptcoll_3000, 0., 0.); };
    //
    obsThFactory["mueeZHcc250_p80_m30"] = [=](const StandardModel& SM) { return new mueeZHcc(SM, sqrt_s_leptcoll_250, pol_80, -pol_30); };
    obsThFactory["mueeZHcc250_m80_p30"] = [=](const StandardModel& SM) { return new mueeZHcc(SM, sqrt_s_leptcoll_250, -pol_80, pol_30); };
    obsThFactory["mueeZHcc250_p80_0"] = [=](const StandardModel& SM) { return new mueeZHcc(SM, sqrt_s_leptcoll_250, pol_80, pol_0); };
    obsThFactory["mueeZHcc250_m80_0"] = [=](const StandardModel& SM) { return new mueeZHcc(SM, sqrt_s_leptcoll_250, -pol_80, pol_0); };
    //
    obsThFactory["mueeZHcc350_p80_m30"] = [=](const StandardModel& SM) { return new mueeZHcc(SM, sqrt_s_leptcoll_350, pol_80, -pol_30); };
    obsThFactory["mueeZHcc350_m80_p30"] = [=](const StandardModel& SM) { return new mueeZHcc(SM, sqrt_s_leptcoll_350, -pol_80, pol_30); };
    obsThFactory["mueeZHcc350_p80_0"] = [=](const StandardModel& SM) { return new mueeZHcc(SM, sqrt_s_leptcoll_350, pol_80, pol_0); };
    obsThFactory["mueeZHcc350_m80_0"] = [=](const StandardModel& SM) { return new mueeZHcc(SM, sqrt_s_leptcoll_350, -pol_80, pol_0); };
    //
    obsThFactory["mueeZHcc380_p80_m30"] = [=](const StandardModel& SM) { return new mueeZHcc(SM, sqrt_s_leptcoll_380, pol_80, -pol_30); };
    obsThFactory["mueeZHcc380_m80_p30"] = [=](const StandardModel& SM) { return new mueeZHcc(SM, sqrt_s_leptcoll_380, -pol_80, pol_30); };
    obsThFactory["mueeZHcc380_p80_0"] = [=](const StandardModel& SM) { return new mueeZHcc(SM, sqrt_s_leptcoll_380, pol_80, pol_0); };
    obsThFactory["mueeZHcc380_m80_0"] = [=](const StandardModel& SM) { return new mueeZHcc(SM, sqrt_s_leptcoll_380, -pol_80, pol_0); };
    //
    obsThFactory["mueeZHcc500_p80_m30"] = [=](const StandardModel& SM) { return new mueeZHcc(SM, sqrt_s_leptcoll_500, pol_80, -pol_30); };
    obsThFactory["mueeZHcc500_m80_p30"] = [=](const StandardModel& SM) { return new mueeZHcc(SM, sqrt_s_leptcoll_500, -pol_80, pol_30); };
    obsThFactory["mueeZHcc500_p80_0"] = [=](const StandardModel& SM) { return new mueeZHcc(SM, sqrt_s_leptcoll_500, pol_80, pol_0); };
    obsThFactory["mueeZHcc500_m80_0"] = [=](const StandardModel& SM) { return new mueeZHcc(SM, sqrt_s_leptcoll_500, -pol_80, pol_0); };
    //
    obsThFactory["mueeZHcc550_p80_m30"] = [=](const StandardModel& SM) { return new mueeZHcc(SM, sqrt_s_leptcoll_550, pol_80, -pol_30); };
    obsThFactory["mueeZHcc550_m80_p30"] = [=](const StandardModel& SM) { return new mueeZHcc(SM, sqrt_s_leptcoll_550, -pol_80, pol_30); };
    obsThFactory["mueeZHcc550_p80_0"] = [=](const StandardModel& SM) { return new mueeZHcc(SM, sqrt_s_leptcoll_550, pol_80, pol_0); };
    obsThFactory["mueeZHcc550_m80_0"] = [=](const StandardModel& SM) { return new mueeZHcc(SM, sqrt_s_leptcoll_550, -pol_80, pol_0); };
    //
    obsThFactory["mueeZHcc1000_p80_m30"] = [=](const StandardModel& SM) { return new mueeZHcc(SM, sqrt_s_leptcoll_1000, pol_80, -pol_30); };
    obsThFactory["mueeZHcc1000_m80_p30"] = [=](const StandardModel& SM) { return new mueeZHcc(SM, sqrt_s_leptcoll_1000, -pol_80, pol_30); };
    obsThFactory["mueeZHcc1000_p80_m20"] = [=](const StandardModel& SM) { return new mueeZHcc(SM, sqrt_s_leptcoll_1000, pol_80, -pol_20); };
    obsThFactory["mueeZHcc1000_m80_p20"] = [=](const StandardModel& SM) { return new mueeZHcc(SM, sqrt_s_leptcoll_1000, -pol_80, pol_20); };
    obsThFactory["mueeZHcc1000_p80_0"] = [=](const StandardModel& SM) { return new mueeZHcc(SM, sqrt_s_leptcoll_1000, pol_80, pol_0); };
    obsThFactory["mueeZHcc1000_m80_0"] = [=](const StandardModel& SM) { return new mueeZHcc(SM, sqrt_s_leptcoll_1000, -pol_80, pol_0); };
    //
    obsThFactory["mueeZHcc1400_p80_m30"] = [=](const StandardModel& SM) { return new mueeZHcc(SM, sqrt_s_leptcoll_1400, pol_80, -pol_30); };
    obsThFactory["mueeZHcc1400_m80_p30"] = [=](const StandardModel& SM) { return new mueeZHcc(SM, sqrt_s_leptcoll_1400, -pol_80, pol_30); };
    obsThFactory["mueeZHcc1400_p80_0"] = [=](const StandardModel& SM) { return new mueeZHcc(SM, sqrt_s_leptcoll_1400, pol_80, pol_0); };
    obsThFactory["mueeZHcc1400_m80_0"] = [=](const StandardModel& SM) { return new mueeZHcc(SM, sqrt_s_leptcoll_1400, -pol_80, pol_0); };
    //
    obsThFactory["mueeZHcc1500_p80_m30"] = [=](const StandardModel& SM) { return new mueeZHcc(SM, sqrt_s_leptcoll_1500, pol_80, -pol_30); };
    obsThFactory["mueeZHcc1500_m80_p30"] = [=](const StandardModel& SM) { return new mueeZHcc(SM, sqrt_s_leptcoll_1500, -pol_80, pol_30); };
    obsThFactory["mueeZHcc1500_p80_0"] = [=](const StandardModel& SM) { return new mueeZHcc(SM, sqrt_s_leptcoll_1500, pol_80, pol_0); };
    obsThFactory["mueeZHcc1500_m80_0"] = [=](const StandardModel& SM) { return new mueeZHcc(SM, sqrt_s_leptcoll_1500, -pol_80, pol_0); };
    //
    obsThFactory["mueeZHcc3000_p80_m30"] = [=](const StandardModel& SM) { return new mueeZHcc(SM, sqrt_s_leptcoll_3000, pol_80, -pol_30); };
    obsThFactory["mueeZHcc3000_m80_p30"] = [=](const StandardModel& SM) { return new mueeZHcc(SM, sqrt_s_leptcoll_3000, -pol_80, pol_30); };
    obsThFactory["mueeZHcc3000_p80_0"] = [=](const StandardModel& SM) { return new mueeZHcc(SM, sqrt_s_leptcoll_3000, pol_80, pol_0); };
    obsThFactory["mueeZHcc3000_m80_0"] = [=](const StandardModel& SM) { return new mueeZHcc(SM, sqrt_s_leptcoll_3000, -pol_80, pol_0); };
    //
    // H -> ss
    obsThFactory["mueeZHss230"] = [=](const StandardModel& SM) { return new mueeZHss(SM, sqrt_s_leptcoll_230, 0., 0.); };
    obsThFactory["mueeZHss240"] = [=](const StandardModel& SM) { return new mueeZHss(SM, sqrt_s_leptcoll_240, 0., 0.); };
    obsThFactory["mueeZHss250"] = [=](const StandardModel& SM) { return new mueeZHss(SM, sqrt_s_leptcoll_250, 0., 0.); };
    obsThFactory["mueeZHss350"] = [=](const StandardModel& SM) { return new mueeZHss(SM, sqrt_s_leptcoll_350, 0., 0.); };
    obsThFactory["mueeZHss365"] = [=](const StandardModel& SM) { return new mueeZHss(SM, sqrt_s_leptcoll_365, 0., 0.); };
    obsThFactory["mueeZHss380"] = [=](const StandardModel& SM) { return new mueeZHss(SM, sqrt_s_leptcoll_380, 0., 0.); };
    obsThFactory["mueeZHss500"] = [=](const StandardModel& SM) { return new mueeZHss(SM, sqrt_s_leptcoll_500, 0., 0.); };
    obsThFactory["mueeZHss550"] = [=](const StandardModel& SM) { return new mueeZHss(SM, sqrt_s_leptcoll_550, 0., 0.); };
    obsThFactory["mueeZHss1000"] = [=](const StandardModel& SM) { return new mueeZHss(SM, sqrt_s_leptcoll_1000, 0., 0.); };
    obsThFactory["mueeZHss1400"] = [=](const StandardModel& SM) { return new mueeZHss(SM, sqrt_s_leptcoll_1400, 0., 0.); };
    obsThFactory["mueeZHss1500"] = [=](const StandardModel& SM) { return new mueeZHss(SM, sqrt_s_leptcoll_1500, 0., 0.); };
    obsThFactory["mueeZHss3000"] = [=](const StandardModel& SM) { return new mueeZHss(SM, sqrt_s_leptcoll_3000, 0., 0.); };
    //
    obsThFactory["mueeZHss250_p80_m30"] = [=](const StandardModel& SM) { return new mueeZHss(SM, sqrt_s_leptcoll_250, pol_80, -pol_30); };
    obsThFactory["mueeZHss250_m80_p30"] = [=](const StandardModel& SM) { return new mueeZHss(SM, sqrt_s_leptcoll_250, -pol_80, pol_30); };
    obsThFactory["mueeZHss250_p80_0"] = [=](const StandardModel& SM) { return new mueeZHss(SM, sqrt_s_leptcoll_250, pol_80, pol_0); };
    obsThFactory["mueeZHss250_m80_0"] = [=](const StandardModel& SM) { return new mueeZHss(SM, sqrt_s_leptcoll_250, -pol_80, pol_0); };
    //
    obsThFactory["mueeZHss350_p80_m30"] = [=](const StandardModel& SM) { return new mueeZHss(SM, sqrt_s_leptcoll_350, pol_80, -pol_30); };
    obsThFactory["mueeZHss350_m80_p30"] = [=](const StandardModel& SM) { return new mueeZHss(SM, sqrt_s_leptcoll_350, -pol_80, pol_30); };
    obsThFactory["mueeZHss350_p80_0"] = [=](const StandardModel& SM) { return new mueeZHss(SM, sqrt_s_leptcoll_350, pol_80, pol_0); };
    obsThFactory["mueeZHss350_m80_0"] = [=](const StandardModel& SM) { return new mueeZHss(SM, sqrt_s_leptcoll_350, -pol_80, pol_0); };
    //
    obsThFactory["mueeZHss380_p80_m30"] = [=](const StandardModel& SM) { return new mueeZHss(SM, sqrt_s_leptcoll_380, pol_80, -pol_30); };
    obsThFactory["mueeZHss380_m80_p30"] = [=](const StandardModel& SM) { return new mueeZHss(SM, sqrt_s_leptcoll_380, -pol_80, pol_30); };
    obsThFactory["mueeZHss380_p80_0"] = [=](const StandardModel& SM) { return new mueeZHss(SM, sqrt_s_leptcoll_380, pol_80, pol_0); };
    obsThFactory["mueeZHss380_m80_0"] = [=](const StandardModel& SM) { return new mueeZHss(SM, sqrt_s_leptcoll_380, -pol_80, pol_0); };
    //
    obsThFactory["mueeZHss500_p80_m30"] = [=](const StandardModel& SM) { return new mueeZHss(SM, sqrt_s_leptcoll_500, pol_80, -pol_30); };
    obsThFactory["mueeZHss500_m80_p30"] = [=](const StandardModel& SM) { return new mueeZHss(SM, sqrt_s_leptcoll_500, -pol_80, pol_30); };
    obsThFactory["mueeZHss500_p80_0"] = [=](const StandardModel& SM) { return new mueeZHss(SM, sqrt_s_leptcoll_500, pol_80, pol_0); };
    obsThFactory["mueeZHss500_m80_0"] = [=](const StandardModel& SM) { return new mueeZHss(SM, sqrt_s_leptcoll_500, -pol_80, pol_0); };
    //
    obsThFactory["mueeZHss550_p80_m30"] = [=](const StandardModel& SM) { return new mueeZHss(SM, sqrt_s_leptcoll_550, pol_80, -pol_30); };
    obsThFactory["mueeZHss550_m80_p30"] = [=](const StandardModel& SM) { return new mueeZHss(SM, sqrt_s_leptcoll_550, -pol_80, pol_30); };
    obsThFactory["mueeZHss550_p80_0"] = [=](const StandardModel& SM) { return new mueeZHss(SM, sqrt_s_leptcoll_550, pol_80, pol_0); };
    obsThFactory["mueeZHss550_m80_0"] = [=](const StandardModel& SM) { return new mueeZHss(SM, sqrt_s_leptcoll_550, -pol_80, pol_0); };
    //
    obsThFactory["mueeZHss1000_p80_m30"] = [=](const StandardModel& SM) { return new mueeZHss(SM, sqrt_s_leptcoll_1000, pol_80, -pol_30); };
    obsThFactory["mueeZHss1000_m80_p30"] = [=](const StandardModel& SM) { return new mueeZHss(SM, sqrt_s_leptcoll_1000, -pol_80, pol_30); };
    obsThFactory["mueeZHss1000_p80_m20"] = [=](const StandardModel& SM) { return new mueeZHss(SM, sqrt_s_leptcoll_1000, pol_80, -pol_20); };
    obsThFactory["mueeZHss1000_m80_p20"] = [=](const StandardModel& SM) { return new mueeZHss(SM, sqrt_s_leptcoll_1000, -pol_80, pol_20); };
    obsThFactory["mueeZHss1000_p80_0"] = [=](const StandardModel& SM) { return new mueeZHss(SM, sqrt_s_leptcoll_1000, pol_80, pol_0); };
    obsThFactory["mueeZHss1000_m80_0"] = [=](const StandardModel& SM) { return new mueeZHss(SM, sqrt_s_leptcoll_1000, -pol_80, pol_0); };
    //
    obsThFactory["mueeZHss1400_p80_m30"] = [=](const StandardModel& SM) { return new mueeZHss(SM, sqrt_s_leptcoll_1400, pol_80, -pol_30); };
    obsThFactory["mueeZHss1400_m80_p30"] = [=](const StandardModel& SM) { return new mueeZHss(SM, sqrt_s_leptcoll_1400, -pol_80, pol_30); };
    obsThFactory["mueeZHss1400_p80_0"] = [=](const StandardModel& SM) { return new mueeZHss(SM, sqrt_s_leptcoll_1400, pol_80, pol_0); };
    obsThFactory["mueeZHss1400_m80_0"] = [=](const StandardModel& SM) { return new mueeZHss(SM, sqrt_s_leptcoll_1400, -pol_80, pol_0); };
    //
    obsThFactory["mueeZHss1500_p80_m30"] = [=](const StandardModel& SM) { return new mueeZHss(SM, sqrt_s_leptcoll_1500, pol_80, -pol_30); };
    obsThFactory["mueeZHss1500_m80_p30"] = [=](const StandardModel& SM) { return new mueeZHss(SM, sqrt_s_leptcoll_1500, -pol_80, pol_30); };
    obsThFactory["mueeZHss1500_p80_0"] = [=](const StandardModel& SM) { return new mueeZHss(SM, sqrt_s_leptcoll_1500, pol_80, pol_0); };
    obsThFactory["mueeZHss1500_m80_0"] = [=](const StandardModel& SM) { return new mueeZHss(SM, sqrt_s_leptcoll_1500, -pol_80, pol_0); };
    //
    obsThFactory["mueeZHss3000_p80_m30"] = [=](const StandardModel& SM) { return new mueeZHss(SM, sqrt_s_leptcoll_3000, pol_80, -pol_30); };
    obsThFactory["mueeZHss3000_m80_p30"] = [=](const StandardModel& SM) { return new mueeZHss(SM, sqrt_s_leptcoll_3000, -pol_80, pol_30); };
    obsThFactory["mueeZHss3000_p80_0"] = [=](const StandardModel& SM) { return new mueeZHss(SM, sqrt_s_leptcoll_3000, pol_80, pol_0); };
    obsThFactory["mueeZHss3000_m80_0"] = [=](const StandardModel& SM) { return new mueeZHss(SM, sqrt_s_leptcoll_3000, -pol_80, pol_0); };
    //
    // H -> gg
    obsThFactory["mueeZHgg230"] = [=](const StandardModel& SM) { return new mueeZHgg(SM, sqrt_s_leptcoll_230, 0., 0.); };
    obsThFactory["mueeZHgg240"] = [=](const StandardModel& SM) { return new mueeZHgg(SM, sqrt_s_leptcoll_240, 0., 0.); };
    obsThFactory["mueeZHgg250"] = [=](const StandardModel& SM) { return new mueeZHgg(SM, sqrt_s_leptcoll_250, 0., 0.); };
    obsThFactory["mueeZHgg350"] = [=](const StandardModel& SM) { return new mueeZHgg(SM, sqrt_s_leptcoll_350, 0., 0.); };
    obsThFactory["mueeZHgg365"] = [=](const StandardModel& SM) { return new mueeZHgg(SM, sqrt_s_leptcoll_365, 0., 0.); };
    obsThFactory["mueeZHgg380"] = [=](const StandardModel& SM) { return new mueeZHgg(SM, sqrt_s_leptcoll_380, 0., 0.); };
    obsThFactory["mueeZHgg500"] = [=](const StandardModel& SM) { return new mueeZHgg(SM, sqrt_s_leptcoll_500, 0., 0.); };
    obsThFactory["mueeZHgg550"] = [=](const StandardModel& SM) { return new mueeZHgg(SM, sqrt_s_leptcoll_550, 0., 0.); };
    obsThFactory["mueeZHgg1000"] = [=](const StandardModel& SM) { return new mueeZHgg(SM, sqrt_s_leptcoll_1000, 0., 0.); };
    obsThFactory["mueeZHgg1400"] = [=](const StandardModel& SM) { return new mueeZHgg(SM, sqrt_s_leptcoll_1400, 0., 0.); };
    obsThFactory["mueeZHgg1500"] = [=](const StandardModel& SM) { return new mueeZHgg(SM, sqrt_s_leptcoll_1500, 0., 0.); };
    obsThFactory["mueeZHgg3000"] = [=](const StandardModel& SM) { return new mueeZHgg(SM, sqrt_s_leptcoll_3000, 0., 0.); };
    //
    obsThFactory["mueeZHgg250_p80_m30"] = [=](const StandardModel& SM) { return new mueeZHgg(SM, sqrt_s_leptcoll_250, pol_80, -pol_30); };
    obsThFactory["mueeZHgg250_m80_p30"] = [=](const StandardModel& SM) { return new mueeZHgg(SM, sqrt_s_leptcoll_250, -pol_80, pol_30); };
    obsThFactory["mueeZHgg250_p80_0"] = [=](const StandardModel& SM) { return new mueeZHgg(SM, sqrt_s_leptcoll_250, pol_80, pol_0); };
    obsThFactory["mueeZHgg250_m80_0"] = [=](const StandardModel& SM) { return new mueeZHgg(SM, sqrt_s_leptcoll_250, -pol_80, pol_0); };
    //
    obsThFactory["mueeZHgg350_p80_m30"] = [=](const StandardModel& SM) { return new mueeZHgg(SM, sqrt_s_leptcoll_350, pol_80, -pol_30); };
    obsThFactory["mueeZHgg350_m80_p30"] = [=](const StandardModel& SM) { return new mueeZHgg(SM, sqrt_s_leptcoll_350, -pol_80, pol_30); };
    obsThFactory["mueeZHgg350_p80_0"] = [=](const StandardModel& SM) { return new mueeZHgg(SM, sqrt_s_leptcoll_350, pol_80, pol_0); };
    obsThFactory["mueeZHgg350_m80_0"] = [=](const StandardModel& SM) { return new mueeZHgg(SM, sqrt_s_leptcoll_350, -pol_80, pol_0); };
    //
    obsThFactory["mueeZHgg380_p80_m30"] = [=](const StandardModel& SM) { return new mueeZHgg(SM, sqrt_s_leptcoll_380, pol_80, -pol_30); };
    obsThFactory["mueeZHgg380_m80_p30"] = [=](const StandardModel& SM) { return new mueeZHgg(SM, sqrt_s_leptcoll_380, -pol_80, pol_30); };
    obsThFactory["mueeZHgg380_p80_0"] = [=](const StandardModel& SM) { return new mueeZHgg(SM, sqrt_s_leptcoll_380, pol_80, pol_0); };
    obsThFactory["mueeZHgg380_m80_0"] = [=](const StandardModel& SM) { return new mueeZHgg(SM, sqrt_s_leptcoll_380, -pol_80, pol_0); };
    //
    obsThFactory["mueeZHgg500_p80_m30"] = [=](const StandardModel& SM) { return new mueeZHgg(SM, sqrt_s_leptcoll_500, pol_80, -pol_30); };
    obsThFactory["mueeZHgg500_m80_p30"] = [=](const StandardModel& SM) { return new mueeZHgg(SM, sqrt_s_leptcoll_500, -pol_80, pol_30); };
    obsThFactory["mueeZHgg500_p80_0"] = [=](const StandardModel& SM) { return new mueeZHgg(SM, sqrt_s_leptcoll_500, pol_80, pol_0); };
    obsThFactory["mueeZHgg500_m80_0"] = [=](const StandardModel& SM) { return new mueeZHgg(SM, sqrt_s_leptcoll_500, -pol_80, pol_0); };
    //
    obsThFactory["mueeZHgg550_p80_m30"] = [=](const StandardModel& SM) { return new mueeZHgg(SM, sqrt_s_leptcoll_550, pol_80, -pol_30); };
    obsThFactory["mueeZHgg550_m80_p30"] = [=](const StandardModel& SM) { return new mueeZHgg(SM, sqrt_s_leptcoll_550, -pol_80, pol_30); };
    obsThFactory["mueeZHgg550_p80_0"] = [=](const StandardModel& SM) { return new mueeZHgg(SM, sqrt_s_leptcoll_550, pol_80, pol_0); };
    obsThFactory["mueeZHgg550_m80_0"] = [=](const StandardModel& SM) { return new mueeZHgg(SM, sqrt_s_leptcoll_550, -pol_80, pol_0); };
    //
    obsThFactory["mueeZHgg1000_p80_m30"] = [=](const StandardModel& SM) { return new mueeZHgg(SM, sqrt_s_leptcoll_1000, pol_80, -pol_30); };
    obsThFactory["mueeZHgg1000_m80_p30"] = [=](const StandardModel& SM) { return new mueeZHgg(SM, sqrt_s_leptcoll_1000, -pol_80, pol_30); };
    obsThFactory["mueeZHgg1000_p80_m20"] = [=](const StandardModel& SM) { return new mueeZHgg(SM, sqrt_s_leptcoll_1000, pol_80, -pol_20); };
    obsThFactory["mueeZHgg1000_m80_p20"] = [=](const StandardModel& SM) { return new mueeZHgg(SM, sqrt_s_leptcoll_1000, -pol_80, pol_20); };
    obsThFactory["mueeZHgg1000_p80_0"] = [=](const StandardModel& SM) { return new mueeZHgg(SM, sqrt_s_leptcoll_1000, pol_80, pol_0); };
    obsThFactory["mueeZHgg1000_m80_0"] = [=](const StandardModel& SM) { return new mueeZHgg(SM, sqrt_s_leptcoll_1000, -pol_80, pol_0); };
    //
    obsThFactory["mueeZHgg1400_p80_m30"] = [=](const StandardModel& SM) { return new mueeZHgg(SM, sqrt_s_leptcoll_1400, pol_80, -pol_30); };
    obsThFactory["mueeZHgg1400_m80_p30"] = [=](const StandardModel& SM) { return new mueeZHgg(SM, sqrt_s_leptcoll_1400, -pol_80, pol_30); };
    obsThFactory["mueeZHgg1400_p80_0"] = [=](const StandardModel& SM) { return new mueeZHgg(SM, sqrt_s_leptcoll_1400, pol_80, pol_0); };
    obsThFactory["mueeZHgg1400_m80_0"] = [=](const StandardModel& SM) { return new mueeZHgg(SM, sqrt_s_leptcoll_1400, -pol_80, pol_0); };
    //
    obsThFactory["mueeZHgg1500_p80_m30"] = [=](const StandardModel& SM) { return new mueeZHgg(SM, sqrt_s_leptcoll_1500, pol_80, -pol_30); };
    obsThFactory["mueeZHgg1500_m80_p30"] = [=](const StandardModel& SM) { return new mueeZHgg(SM, sqrt_s_leptcoll_1500, -pol_80, pol_30); };
    obsThFactory["mueeZHgg1500_p80_0"] = [=](const StandardModel& SM) { return new mueeZHgg(SM, sqrt_s_leptcoll_1500, pol_80, pol_0); };
    obsThFactory["mueeZHgg1500_m80_0"] = [=](const StandardModel& SM) { return new mueeZHgg(SM, sqrt_s_leptcoll_1500, -pol_80, pol_0); };
    //
    obsThFactory["mueeZHgg3000_p80_m30"] = [=](const StandardModel& SM) { return new mueeZHgg(SM, sqrt_s_leptcoll_3000, pol_80, -pol_30); };
    obsThFactory["mueeZHgg3000_m80_p30"] = [=](const StandardModel& SM) { return new mueeZHgg(SM, sqrt_s_leptcoll_3000, -pol_80, pol_30); };
    obsThFactory["mueeZHgg3000_p80_0"] = [=](const StandardModel& SM) { return new mueeZHgg(SM, sqrt_s_leptcoll_3000, pol_80, pol_0); };
    obsThFactory["mueeZHgg3000_m80_0"] = [=](const StandardModel& SM) { return new mueeZHgg(SM, sqrt_s_leptcoll_3000, -pol_80, pol_0); };
    //
    // H -> WW
    obsThFactory["mueeZHWW230"] = [=](const StandardModel& SM) { return new mueeZHWW(SM, sqrt_s_leptcoll_230, 0., 0.); };
    obsThFactory["mueeZHWW240"] = [=](const StandardModel& SM) { return new mueeZHWW(SM, sqrt_s_leptcoll_240, 0., 0.); };
    obsThFactory["mueeZHWW250"] = [=](const StandardModel& SM) { return new mueeZHWW(SM, sqrt_s_leptcoll_250, 0., 0.); };
    obsThFactory["mueeZHWW350"] = [=](const StandardModel& SM) { return new mueeZHWW(SM, sqrt_s_leptcoll_350, 0., 0.); };
    obsThFactory["mueeZHWW365"] = [=](const StandardModel& SM) { return new mueeZHWW(SM, sqrt_s_leptcoll_365, 0., 0.); };
    obsThFactory["mueeZHWW380"] = [=](const StandardModel& SM) { return new mueeZHWW(SM, sqrt_s_leptcoll_380, 0., 0.); };
    obsThFactory["mueeZHWW500"] = [=](const StandardModel& SM) { return new mueeZHWW(SM, sqrt_s_leptcoll_500, 0., 0.); };
    obsThFactory["mueeZHWW550"] = [=](const StandardModel& SM) { return new mueeZHWW(SM, sqrt_s_leptcoll_550, 0., 0.); };
    obsThFactory["mueeZHWW1000"] = [=](const StandardModel& SM) { return new mueeZHWW(SM, sqrt_s_leptcoll_1000, 0., 0.); };
    obsThFactory["mueeZHWW1400"] = [=](const StandardModel& SM) { return new mueeZHWW(SM, sqrt_s_leptcoll_1400, 0., 0.); };
    obsThFactory["mueeZHWW1500"] = [=](const StandardModel& SM) { return new mueeZHWW(SM, sqrt_s_leptcoll_1500, 0., 0.); };
    obsThFactory["mueeZHWW3000"] = [=](const StandardModel& SM) { return new mueeZHWW(SM, sqrt_s_leptcoll_3000, 0., 0.); };
    //
    obsThFactory["mueeZHWW250_p80_m30"] = [=](const StandardModel& SM) { return new mueeZHWW(SM, sqrt_s_leptcoll_250, pol_80, -pol_30); };
    obsThFactory["mueeZHWW250_m80_p30"] = [=](const StandardModel& SM) { return new mueeZHWW(SM, sqrt_s_leptcoll_250, -pol_80, pol_30); };
    obsThFactory["mueeZHWW250_p80_0"] = [=](const StandardModel& SM) { return new mueeZHWW(SM, sqrt_s_leptcoll_250, pol_80, pol_0); };
    obsThFactory["mueeZHWW250_m80_0"] = [=](const StandardModel& SM) { return new mueeZHWW(SM, sqrt_s_leptcoll_250, -pol_80, pol_0); };
    //
    obsThFactory["mueeZHWW350_p80_m30"] = [=](const StandardModel& SM) { return new mueeZHWW(SM, sqrt_s_leptcoll_350, pol_80, -pol_30); };
    obsThFactory["mueeZHWW350_m80_p30"] = [=](const StandardModel& SM) { return new mueeZHWW(SM, sqrt_s_leptcoll_350, -pol_80, pol_30); };
    obsThFactory["mueeZHWW350_p80_0"] = [=](const StandardModel& SM) { return new mueeZHWW(SM, sqrt_s_leptcoll_350, pol_80, pol_0); };
    obsThFactory["mueeZHWW350_m80_0"] = [=](const StandardModel& SM) { return new mueeZHWW(SM, sqrt_s_leptcoll_350, -pol_80, pol_0); };
    //
    obsThFactory["mueeZHWW380_p80_m30"] = [=](const StandardModel& SM) { return new mueeZHWW(SM, sqrt_s_leptcoll_380, pol_80, -pol_30); };
    obsThFactory["mueeZHWW380_m80_p30"] = [=](const StandardModel& SM) { return new mueeZHWW(SM, sqrt_s_leptcoll_380, -pol_80, pol_30); };
    obsThFactory["mueeZHWW380_p80_0"] = [=](const StandardModel& SM) { return new mueeZHWW(SM, sqrt_s_leptcoll_380, pol_80, pol_0); };
    obsThFactory["mueeZHWW380_m80_0"] = [=](const StandardModel& SM) { return new mueeZHWW(SM, sqrt_s_leptcoll_380, -pol_80, pol_0); };
    //
    obsThFactory["mueeZHWW500_p80_m30"] = [=](const StandardModel& SM) { return new mueeZHWW(SM, sqrt_s_leptcoll_500, pol_80, -pol_30); };
    obsThFactory["mueeZHWW500_m80_p30"] = [=](const StandardModel& SM) { return new mueeZHWW(SM, sqrt_s_leptcoll_500, -pol_80, pol_30); };
    obsThFactory["mueeZHWW500_p80_0"] = [=](const StandardModel& SM) { return new mueeZHWW(SM, sqrt_s_leptcoll_500, pol_80, pol_0); };
    obsThFactory["mueeZHWW500_m80_0"] = [=](const StandardModel& SM) { return new mueeZHWW(SM, sqrt_s_leptcoll_500, -pol_80, pol_0); };
    //
    obsThFactory["mueeZHWW550_p80_m30"] = [=](const StandardModel& SM) { return new mueeZHWW(SM, sqrt_s_leptcoll_550, pol_80, -pol_30); };
    obsThFactory["mueeZHWW550_m80_p30"] = [=](const StandardModel& SM) { return new mueeZHWW(SM, sqrt_s_leptcoll_550, -pol_80, pol_30); };
    obsThFactory["mueeZHWW550_p80_0"] = [=](const StandardModel& SM) { return new mueeZHWW(SM, sqrt_s_leptcoll_550, pol_80, pol_0); };
    obsThFactory["mueeZHWW550_m80_0"] = [=](const StandardModel& SM) { return new mueeZHWW(SM, sqrt_s_leptcoll_550, -pol_80, pol_0); };
    //
    obsThFactory["mueeZHWW1000_p80_m30"] = [=](const StandardModel& SM) { return new mueeZHWW(SM, sqrt_s_leptcoll_1000, pol_80, -pol_30); };
    obsThFactory["mueeZHWW1000_m80_p30"] = [=](const StandardModel& SM) { return new mueeZHWW(SM, sqrt_s_leptcoll_1000, -pol_80, pol_30); };
    obsThFactory["mueeZHWW1000_p80_m20"] = [=](const StandardModel& SM) { return new mueeZHWW(SM, sqrt_s_leptcoll_1000, pol_80, -pol_20); };
    obsThFactory["mueeZHWW1000_m80_p20"] = [=](const StandardModel& SM) { return new mueeZHWW(SM, sqrt_s_leptcoll_1000, -pol_80, pol_20); };
    obsThFactory["mueeZHWW1000_p80_0"] = [=](const StandardModel& SM) { return new mueeZHWW(SM, sqrt_s_leptcoll_1000, pol_80, pol_0); };
    obsThFactory["mueeZHWW1000_m80_0"] = [=](const StandardModel& SM) { return new mueeZHWW(SM, sqrt_s_leptcoll_1000, -pol_80, pol_0); };
    //
    obsThFactory["mueeZHWW1400_p80_m30"] = [=](const StandardModel& SM) { return new mueeZHWW(SM, sqrt_s_leptcoll_1400, pol_80, -pol_30); };
    obsThFactory["mueeZHWW1400_m80_p30"] = [=](const StandardModel& SM) { return new mueeZHWW(SM, sqrt_s_leptcoll_1400, -pol_80, pol_30); };
    obsThFactory["mueeZHWW1400_p80_0"] = [=](const StandardModel& SM) { return new mueeZHWW(SM, sqrt_s_leptcoll_1400, pol_80, pol_0); };
    obsThFactory["mueeZHWW1400_m80_0"] = [=](const StandardModel& SM) { return new mueeZHWW(SM, sqrt_s_leptcoll_1400, -pol_80, pol_0); };
    //
    obsThFactory["mueeZHWW1500_p80_m30"] = [=](const StandardModel& SM) { return new mueeZHWW(SM, sqrt_s_leptcoll_1500, pol_80, -pol_30); };
    obsThFactory["mueeZHWW1500_m80_p30"] = [=](const StandardModel& SM) { return new mueeZHWW(SM, sqrt_s_leptcoll_1500, -pol_80, pol_30); };
    obsThFactory["mueeZHWW1500_p80_0"] = [=](const StandardModel& SM) { return new mueeZHWW(SM, sqrt_s_leptcoll_1500, pol_80, pol_0); };
    obsThFactory["mueeZHWW1500_m80_0"] = [=](const StandardModel& SM) { return new mueeZHWW(SM, sqrt_s_leptcoll_1500, -pol_80, pol_0); };
    //
    obsThFactory["mueeZHWW3000_p80_m30"] = [=](const StandardModel& SM) { return new mueeZHWW(SM, sqrt_s_leptcoll_3000, pol_80, -pol_30); };
    obsThFactory["mueeZHWW3000_m80_p30"] = [=](const StandardModel& SM) { return new mueeZHWW(SM, sqrt_s_leptcoll_3000, -pol_80, pol_30); };
    obsThFactory["mueeZHWW3000_p80_0"] = [=](const StandardModel& SM) { return new mueeZHWW(SM, sqrt_s_leptcoll_3000, pol_80, pol_0); };
    obsThFactory["mueeZHWW3000_m80_0"] = [=](const StandardModel& SM) { return new mueeZHWW(SM, sqrt_s_leptcoll_3000, -pol_80, pol_0); };
    //
    // H -> tau tau
    obsThFactory["mueeZHtautau230"] = [=](const StandardModel& SM) { return new mueeZHtautau(SM, sqrt_s_leptcoll_230, 0., 0.); };
    obsThFactory["mueeZHtautau240"] = [=](const StandardModel& SM) { return new mueeZHtautau(SM, sqrt_s_leptcoll_240, 0., 0.); };
    obsThFactory["mueeZHtautau250"] = [=](const StandardModel& SM) { return new mueeZHtautau(SM, sqrt_s_leptcoll_250, 0., 0.); };
    obsThFactory["mueeZHtautau350"] = [=](const StandardModel& SM) { return new mueeZHtautau(SM, sqrt_s_leptcoll_350, 0., 0.); };
    obsThFactory["mueeZHtautau365"] = [=](const StandardModel& SM) { return new mueeZHtautau(SM, sqrt_s_leptcoll_365, 0., 0.); };
    obsThFactory["mueeZHtautau380"] = [=](const StandardModel& SM) { return new mueeZHtautau(SM, sqrt_s_leptcoll_380, 0., 0.); };
    obsThFactory["mueeZHtautau500"] = [=](const StandardModel& SM) { return new mueeZHtautau(SM, sqrt_s_leptcoll_500, 0., 0.); };
    obsThFactory["mueeZHtautau550"] = [=](const StandardModel& SM) { return new mueeZHtautau(SM, sqrt_s_leptcoll_550, 0., 0.); };
    obsThFactory["mueeZHtautau1000"] = [=](const StandardModel& SM) { return new mueeZHtautau(SM, sqrt_s_leptcoll_1000, 0., 0.); };
    obsThFactory["mueeZHtautau1400"] = [=](const StandardModel& SM) { return new mueeZHtautau(SM, sqrt_s_leptcoll_1400, 0., 0.); };
    obsThFactory["mueeZHtautau1500"] = [=](const StandardModel& SM) { return new mueeZHtautau(SM, sqrt_s_leptcoll_1500, 0., 0.); };
    obsThFactory["mueeZHtautau3000"] = [=](const StandardModel& SM) { return new mueeZHtautau(SM, sqrt_s_leptcoll_3000, 0., 0.); };
    //
    obsThFactory["mueeZHtautau250_p80_m30"] = [=](const StandardModel& SM) { return new mueeZHtautau(SM, sqrt_s_leptcoll_250, pol_80, -pol_30); };
    obsThFactory["mueeZHtautau250_m80_p30"] = [=](const StandardModel& SM) { return new mueeZHtautau(SM, sqrt_s_leptcoll_250, -pol_80, pol_30); };
    obsThFactory["mueeZHtautau250_p80_0"] = [=](const StandardModel& SM) { return new mueeZHtautau(SM, sqrt_s_leptcoll_250, pol_80, pol_0); };
    obsThFactory["mueeZHtautau250_m80_0"] = [=](const StandardModel& SM) { return new mueeZHtautau(SM, sqrt_s_leptcoll_250, -pol_80, pol_0); };
    //
    obsThFactory["mueeZHtautau350_p80_m30"] = [=](const StandardModel& SM) { return new mueeZHtautau(SM, sqrt_s_leptcoll_350, pol_80, -pol_30); };
    obsThFactory["mueeZHtautau350_m80_p30"] = [=](const StandardModel& SM) { return new mueeZHtautau(SM, sqrt_s_leptcoll_350, -pol_80, pol_30); };
    obsThFactory["mueeZHtautau350_p80_0"] = [=](const StandardModel& SM) { return new mueeZHtautau(SM, sqrt_s_leptcoll_350, pol_80, pol_0); };
    obsThFactory["mueeZHtautau350_m80_0"] = [=](const StandardModel& SM) { return new mueeZHtautau(SM, sqrt_s_leptcoll_350, -pol_80, pol_0); };
    //
    obsThFactory["mueeZHtautau380_p80_m30"] = [=](const StandardModel& SM) { return new mueeZHtautau(SM, sqrt_s_leptcoll_380, pol_80, -pol_30); };
    obsThFactory["mueeZHtautau380_m80_p30"] = [=](const StandardModel& SM) { return new mueeZHtautau(SM, sqrt_s_leptcoll_380, -pol_80, pol_30); };
    obsThFactory["mueeZHtautau380_p80_0"] = [=](const StandardModel& SM) { return new mueeZHtautau(SM, sqrt_s_leptcoll_380, pol_80, pol_0); };
    obsThFactory["mueeZHtautau380_m80_0"] = [=](const StandardModel& SM) { return new mueeZHtautau(SM, sqrt_s_leptcoll_380, -pol_80, pol_0); };
    //
    obsThFactory["mueeZHtautau500_p80_m30"] = [=](const StandardModel& SM) { return new mueeZHtautau(SM, sqrt_s_leptcoll_500, pol_80, -pol_30); };
    obsThFactory["mueeZHtautau500_m80_p30"] = [=](const StandardModel& SM) { return new mueeZHtautau(SM, sqrt_s_leptcoll_500, -pol_80, pol_30); };
    obsThFactory["mueeZHtautau500_p80_0"] = [=](const StandardModel& SM) { return new mueeZHtautau(SM, sqrt_s_leptcoll_500, pol_80, pol_0); };
    obsThFactory["mueeZHtautau500_m80_0"] = [=](const StandardModel& SM) { return new mueeZHtautau(SM, sqrt_s_leptcoll_500, -pol_80, pol_0); };
    //
    obsThFactory["mueeZHtautau550_p80_m30"] = [=](const StandardModel& SM) { return new mueeZHtautau(SM, sqrt_s_leptcoll_550, pol_80, -pol_30); };
    obsThFactory["mueeZHtautau550_m80_p30"] = [=](const StandardModel& SM) { return new mueeZHtautau(SM, sqrt_s_leptcoll_550, -pol_80, pol_30); };
    obsThFactory["mueeZHtautau550_p80_0"] = [=](const StandardModel& SM) { return new mueeZHtautau(SM, sqrt_s_leptcoll_550, pol_80, pol_0); };
    obsThFactory["mueeZHtautau550_m80_0"] = [=](const StandardModel& SM) { return new mueeZHtautau(SM, sqrt_s_leptcoll_550, -pol_80, pol_0); };
    //
    obsThFactory["mueeZHtautau1000_p80_m30"] = [=](const StandardModel& SM) { return new mueeZHtautau(SM, sqrt_s_leptcoll_1000, pol_80, -pol_30); };
    obsThFactory["mueeZHtautau1000_m80_p30"] = [=](const StandardModel& SM) { return new mueeZHtautau(SM, sqrt_s_leptcoll_1000, -pol_80, pol_30); };
    obsThFactory["mueeZHtautau1000_p80_m20"] = [=](const StandardModel& SM) { return new mueeZHtautau(SM, sqrt_s_leptcoll_1000, pol_80, -pol_20); };
    obsThFactory["mueeZHtautau1000_m80_p20"] = [=](const StandardModel& SM) { return new mueeZHtautau(SM, sqrt_s_leptcoll_1000, -pol_80, pol_20); };
    obsThFactory["mueeZHtautau1000_p80_0"] = [=](const StandardModel& SM) { return new mueeZHtautau(SM, sqrt_s_leptcoll_1000, pol_80, pol_0); };
    obsThFactory["mueeZHtautau1000_m80_0"] = [=](const StandardModel& SM) { return new mueeZHtautau(SM, sqrt_s_leptcoll_1000, -pol_80, pol_0); };
    //
    obsThFactory["mueeZHtautau1400_p80_m30"] = [=](const StandardModel& SM) { return new mueeZHtautau(SM, sqrt_s_leptcoll_1400, pol_80, -pol_30); };
    obsThFactory["mueeZHtautau1400_m80_p30"] = [=](const StandardModel& SM) { return new mueeZHtautau(SM, sqrt_s_leptcoll_1400, -pol_80, pol_30); };
    obsThFactory["mueeZHtautau1400_p80_0"] = [=](const StandardModel& SM) { return new mueeZHtautau(SM, sqrt_s_leptcoll_1400, pol_80, pol_0); };
    obsThFactory["mueeZHtautau1400_m80_0"] = [=](const StandardModel& SM) { return new mueeZHtautau(SM, sqrt_s_leptcoll_1400, -pol_80, pol_0); };
    //
    obsThFactory["mueeZHtautau1500_p80_m30"] = [=](const StandardModel& SM) { return new mueeZHtautau(SM, sqrt_s_leptcoll_1500, pol_80, -pol_30); };
    obsThFactory["mueeZHtautau1500_m80_p30"] = [=](const StandardModel& SM) { return new mueeZHtautau(SM, sqrt_s_leptcoll_1500, -pol_80, pol_30); };
    obsThFactory["mueeZHtautau1500_p80_0"] = [=](const StandardModel& SM) { return new mueeZHtautau(SM, sqrt_s_leptcoll_1500, pol_80, pol_0); };
    obsThFactory["mueeZHtautau1500_m80_0"] = [=](const StandardModel& SM) { return new mueeZHtautau(SM, sqrt_s_leptcoll_1500, -pol_80, pol_0); };
    //
    obsThFactory["mueeZHtautau3000_p80_m30"] = [=](const StandardModel& SM) { return new mueeZHtautau(SM, sqrt_s_leptcoll_3000, pol_80, -pol_30); };
    obsThFactory["mueeZHtautau3000_m80_p30"] = [=](const StandardModel& SM) { return new mueeZHtautau(SM, sqrt_s_leptcoll_3000, -pol_80, pol_30); };
    obsThFactory["mueeZHtautau3000_p80_0"] = [=](const StandardModel& SM) { return new mueeZHtautau(SM, sqrt_s_leptcoll_3000, pol_80, pol_0); };
    obsThFactory["mueeZHtautau3000_m80_0"] = [=](const StandardModel& SM) { return new mueeZHtautau(SM, sqrt_s_leptcoll_3000, -pol_80, pol_0); };
    //
    // H -> ZZ
    obsThFactory["mueeZHZZ230"] = [=](const StandardModel& SM) { return new mueeZHZZ(SM, sqrt_s_leptcoll_230, 0., 0.); };
    obsThFactory["mueeZHZZ240"] = [=](const StandardModel& SM) { return new mueeZHZZ(SM, sqrt_s_leptcoll_240, 0., 0.); };
    obsThFactory["mueeZHZZ250"] = [=](const StandardModel& SM) { return new mueeZHZZ(SM, sqrt_s_leptcoll_250, 0., 0.); };
    obsThFactory["mueeZHZZ350"] = [=](const StandardModel& SM) { return new mueeZHZZ(SM, sqrt_s_leptcoll_350, 0., 0.); };
    obsThFactory["mueeZHZZ365"] = [=](const StandardModel& SM) { return new mueeZHZZ(SM, sqrt_s_leptcoll_365, 0., 0.); };
    obsThFactory["mueeZHZZ380"] = [=](const StandardModel& SM) { return new mueeZHZZ(SM, sqrt_s_leptcoll_380, 0., 0.); };
    obsThFactory["mueeZHZZ500"] = [=](const StandardModel& SM) { return new mueeZHZZ(SM, sqrt_s_leptcoll_500, 0., 0.); };
    obsThFactory["mueeZHZZ550"] = [=](const StandardModel& SM) { return new mueeZHZZ(SM, sqrt_s_leptcoll_550, 0., 0.); };
    obsThFactory["mueeZHZZ1000"] = [=](const StandardModel& SM) { return new mueeZHZZ(SM, sqrt_s_leptcoll_1000, 0., 0.); };
    obsThFactory["mueeZHZZ1400"] = [=](const StandardModel& SM) { return new mueeZHZZ(SM, sqrt_s_leptcoll_1400, 0., 0.); };
    obsThFactory["mueeZHZZ1500"] = [=](const StandardModel& SM) { return new mueeZHZZ(SM, sqrt_s_leptcoll_1500, 0., 0.); };
    obsThFactory["mueeZHZZ3000"] = [=](const StandardModel& SM) { return new mueeZHZZ(SM, sqrt_s_leptcoll_3000, 0., 0.); };
    //
    obsThFactory["mueeZHZZ250_p80_m30"] = [=](const StandardModel& SM) { return new mueeZHZZ(SM, sqrt_s_leptcoll_250, pol_80, -pol_30); };
    obsThFactory["mueeZHZZ250_m80_p30"] = [=](const StandardModel& SM) { return new mueeZHZZ(SM, sqrt_s_leptcoll_250, -pol_80, pol_30); };
    obsThFactory["mueeZHZZ250_p80_0"] = [=](const StandardModel& SM) { return new mueeZHZZ(SM, sqrt_s_leptcoll_250, pol_80, pol_0); };
    obsThFactory["mueeZHZZ250_m80_0"] = [=](const StandardModel& SM) { return new mueeZHZZ(SM, sqrt_s_leptcoll_250, -pol_80, pol_0); };
    //
    obsThFactory["mueeZHZZ350_p80_m30"] = [=](const StandardModel& SM) { return new mueeZHZZ(SM, sqrt_s_leptcoll_350, pol_80, -pol_30); };
    obsThFactory["mueeZHZZ350_m80_p30"] = [=](const StandardModel& SM) { return new mueeZHZZ(SM, sqrt_s_leptcoll_350, -pol_80, pol_30); };
    obsThFactory["mueeZHZZ350_p80_0"] = [=](const StandardModel& SM) { return new mueeZHZZ(SM, sqrt_s_leptcoll_350, pol_80, pol_0); };
    obsThFactory["mueeZHZZ350_m80_0"] = [=](const StandardModel& SM) { return new mueeZHZZ(SM, sqrt_s_leptcoll_350, -pol_80, pol_0); };
    //
    obsThFactory["mueeZHZZ380_p80_m30"] = [=](const StandardModel& SM) { return new mueeZHZZ(SM, sqrt_s_leptcoll_380, pol_80, -pol_30); };
    obsThFactory["mueeZHZZ380_m80_p30"] = [=](const StandardModel& SM) { return new mueeZHZZ(SM, sqrt_s_leptcoll_380, -pol_80, pol_30); };
    obsThFactory["mueeZHZZ380_p80_0"] = [=](const StandardModel& SM) { return new mueeZHZZ(SM, sqrt_s_leptcoll_380, pol_80, pol_0); };
    obsThFactory["mueeZHZZ380_m80_0"] = [=](const StandardModel& SM) { return new mueeZHZZ(SM, sqrt_s_leptcoll_380, -pol_80, pol_0); };
    //
    obsThFactory["mueeZHZZ500_p80_m30"] = [=](const StandardModel& SM) { return new mueeZHZZ(SM, sqrt_s_leptcoll_500, pol_80, -pol_30); };
    obsThFactory["mueeZHZZ500_m80_p30"] = [=](const StandardModel& SM) { return new mueeZHZZ(SM, sqrt_s_leptcoll_500, -pol_80, pol_30); };
    obsThFactory["mueeZHZZ500_p80_0"] = [=](const StandardModel& SM) { return new mueeZHZZ(SM, sqrt_s_leptcoll_500, pol_80, pol_0); };
    obsThFactory["mueeZHZZ500_m80_0"] = [=](const StandardModel& SM) { return new mueeZHZZ(SM, sqrt_s_leptcoll_500, -pol_80, pol_0); };
    //
    obsThFactory["mueeZHZZ550_p80_m30"] = [=](const StandardModel& SM) { return new mueeZHZZ(SM, sqrt_s_leptcoll_550, pol_80, -pol_30); };
    obsThFactory["mueeZHZZ550_m80_p30"] = [=](const StandardModel& SM) { return new mueeZHZZ(SM, sqrt_s_leptcoll_550, -pol_80, pol_30); };
    obsThFactory["mueeZHZZ550_p80_0"] = [=](const StandardModel& SM) { return new mueeZHZZ(SM, sqrt_s_leptcoll_550, pol_80, pol_0); };
    obsThFactory["mueeZHZZ550_m80_0"] = [=](const StandardModel& SM) { return new mueeZHZZ(SM, sqrt_s_leptcoll_550, -pol_80, pol_0); };
    //
    obsThFactory["mueeZHZZ1000_p80_m30"] = [=](const StandardModel& SM) { return new mueeZHZZ(SM, sqrt_s_leptcoll_1000, pol_80, -pol_30); };
    obsThFactory["mueeZHZZ1000_m80_p30"] = [=](const StandardModel& SM) { return new mueeZHZZ(SM, sqrt_s_leptcoll_1000, -pol_80, pol_30); };
    obsThFactory["mueeZHZZ1000_p80_m20"] = [=](const StandardModel& SM) { return new mueeZHZZ(SM, sqrt_s_leptcoll_1000, pol_80, -pol_20); };
    obsThFactory["mueeZHZZ1000_m80_p20"] = [=](const StandardModel& SM) { return new mueeZHZZ(SM, sqrt_s_leptcoll_1000, -pol_80, pol_20); };
    obsThFactory["mueeZHZZ1000_p80_0"] = [=](const StandardModel& SM) { return new mueeZHZZ(SM, sqrt_s_leptcoll_1000, pol_80, pol_0); };
    obsThFactory["mueeZHZZ1000_m80_0"] = [=](const StandardModel& SM) { return new mueeZHZZ(SM, sqrt_s_leptcoll_1000, -pol_80, pol_0); };
    //
    obsThFactory["mueeZHZZ1400_p80_m30"] = [=](const StandardModel& SM) { return new mueeZHZZ(SM, sqrt_s_leptcoll_1400, pol_80, -pol_30); };
    obsThFactory["mueeZHZZ1400_m80_p30"] = [=](const StandardModel& SM) { return new mueeZHZZ(SM, sqrt_s_leptcoll_1400, -pol_80, pol_30); };
    obsThFactory["mueeZHZZ1400_p80_0"] = [=](const StandardModel& SM) { return new mueeZHZZ(SM, sqrt_s_leptcoll_1400, pol_80, pol_0); };
    obsThFactory["mueeZHZZ1400_m80_0"] = [=](const StandardModel& SM) { return new mueeZHZZ(SM, sqrt_s_leptcoll_1400, -pol_80, pol_0); };
    //
    obsThFactory["mueeZHZZ1500_p80_m30"] = [=](const StandardModel& SM) { return new mueeZHZZ(SM, sqrt_s_leptcoll_1500, pol_80, -pol_30); };
    obsThFactory["mueeZHZZ1500_m80_p30"] = [=](const StandardModel& SM) { return new mueeZHZZ(SM, sqrt_s_leptcoll_1500, -pol_80, pol_30); };
    obsThFactory["mueeZHZZ1500_p80_0"] = [=](const StandardModel& SM) { return new mueeZHZZ(SM, sqrt_s_leptcoll_1500, pol_80, pol_0); };
    obsThFactory["mueeZHZZ1500_m80_0"] = [=](const StandardModel& SM) { return new mueeZHZZ(SM, sqrt_s_leptcoll_1500, -pol_80, pol_0); };
    //
    obsThFactory["mueeZHZZ3000_p80_m30"] = [=](const StandardModel& SM) { return new mueeZHZZ(SM, sqrt_s_leptcoll_3000, pol_80, -pol_30); };
    obsThFactory["mueeZHZZ3000_m80_p30"] = [=](const StandardModel& SM) { return new mueeZHZZ(SM, sqrt_s_leptcoll_3000, -pol_80, pol_30); };
    obsThFactory["mueeZHZZ3000_p80_0"] = [=](const StandardModel& SM) { return new mueeZHZZ(SM, sqrt_s_leptcoll_3000, pol_80, pol_0); };
    obsThFactory["mueeZHZZ3000_m80_0"] = [=](const StandardModel& SM) { return new mueeZHZZ(SM, sqrt_s_leptcoll_3000, -pol_80, pol_0); };
    //
    // H -> Z ga
    obsThFactory["mueeZHZga230"] = [=](const StandardModel& SM) { return new mueeZHZga(SM, sqrt_s_leptcoll_230, 0., 0.); };
    obsThFactory["mueeZHZga240"] = [=](const StandardModel& SM) { return new mueeZHZga(SM, sqrt_s_leptcoll_240, 0., 0.); };
    obsThFactory["mueeZHZga250"] = [=](const StandardModel& SM) { return new mueeZHZga(SM, sqrt_s_leptcoll_250, 0., 0.); };
    obsThFactory["mueeZHZga350"] = [=](const StandardModel& SM) { return new mueeZHZga(SM, sqrt_s_leptcoll_350, 0., 0.); };
    obsThFactory["mueeZHZga365"] = [=](const StandardModel& SM) { return new mueeZHZga(SM, sqrt_s_leptcoll_365, 0., 0.); };
    obsThFactory["mueeZHZga380"] = [=](const StandardModel& SM) { return new mueeZHZga(SM, sqrt_s_leptcoll_380, 0., 0.); };
    obsThFactory["mueeZHZga500"] = [=](const StandardModel& SM) { return new mueeZHZga(SM, sqrt_s_leptcoll_500, 0., 0.); };
    obsThFactory["mueeZHZga550"] = [=](const StandardModel& SM) { return new mueeZHZga(SM, sqrt_s_leptcoll_550, 0., 0.); };
    obsThFactory["mueeZHZga1000"] = [=](const StandardModel& SM) { return new mueeZHZga(SM, sqrt_s_leptcoll_1000, 0., 0.); };
    obsThFactory["mueeZHZga1400"] = [=](const StandardModel& SM) { return new mueeZHZga(SM, sqrt_s_leptcoll_1400, 0., 0.); };
    obsThFactory["mueeZHZga1500"] = [=](const StandardModel& SM) { return new mueeZHZga(SM, sqrt_s_leptcoll_1500, 0., 0.); };
    obsThFactory["mueeZHZga3000"] = [=](const StandardModel& SM) { return new mueeZHZga(SM, sqrt_s_leptcoll_3000, 0., 0.); };
    //
    obsThFactory["mueeZHZga250_p80_m30"] = [=](const StandardModel& SM) { return new mueeZHZga(SM, sqrt_s_leptcoll_250, pol_80, -pol_30); };
    obsThFactory["mueeZHZga250_m80_p30"] = [=](const StandardModel& SM) { return new mueeZHZga(SM, sqrt_s_leptcoll_250, -pol_80, pol_30); };
    obsThFactory["mueeZHZga250_p80_0"] = [=](const StandardModel& SM) { return new mueeZHZga(SM, sqrt_s_leptcoll_250, pol_80, pol_0); };
    obsThFactory["mueeZHZga250_m80_0"] = [=](const StandardModel& SM) { return new mueeZHZga(SM, sqrt_s_leptcoll_250, -pol_80, pol_0); };
    //
    obsThFactory["mueeZHZga350_p80_m30"] = [=](const StandardModel& SM) { return new mueeZHZga(SM, sqrt_s_leptcoll_350, pol_80, -pol_30); };
    obsThFactory["mueeZHZga350_m80_p30"] = [=](const StandardModel& SM) { return new mueeZHZga(SM, sqrt_s_leptcoll_350, -pol_80, pol_30); };
    obsThFactory["mueeZHZga350_p80_0"] = [=](const StandardModel& SM) { return new mueeZHZga(SM, sqrt_s_leptcoll_350, pol_80, pol_0); };
    obsThFactory["mueeZHZga350_m80_0"] = [=](const StandardModel& SM) { return new mueeZHZga(SM, sqrt_s_leptcoll_350, -pol_80, pol_0); };
    //
    obsThFactory["mueeZHZga380_p80_m30"] = [=](const StandardModel& SM) { return new mueeZHZga(SM, sqrt_s_leptcoll_380, pol_80, -pol_30); };
    obsThFactory["mueeZHZga380_m80_p30"] = [=](const StandardModel& SM) { return new mueeZHZga(SM, sqrt_s_leptcoll_380, -pol_80, pol_30); };
    obsThFactory["mueeZHZga380_p80_0"] = [=](const StandardModel& SM) { return new mueeZHZga(SM, sqrt_s_leptcoll_380, pol_80, pol_0); };
    obsThFactory["mueeZHZga380_m80_0"] = [=](const StandardModel& SM) { return new mueeZHZga(SM, sqrt_s_leptcoll_380, -pol_80, pol_0); };
    //
    obsThFactory["mueeZHZga500_p80_m30"] = [=](const StandardModel& SM) { return new mueeZHZga(SM, sqrt_s_leptcoll_500, pol_80, -pol_30); };
    obsThFactory["mueeZHZga500_m80_p30"] = [=](const StandardModel& SM) { return new mueeZHZga(SM, sqrt_s_leptcoll_500, -pol_80, pol_30); };
    obsThFactory["mueeZHZga500_p80_0"] = [=](const StandardModel& SM) { return new mueeZHZga(SM, sqrt_s_leptcoll_500, pol_80, pol_0); };
    obsThFactory["mueeZHZga500_m80_0"] = [=](const StandardModel& SM) { return new mueeZHZga(SM, sqrt_s_leptcoll_500, -pol_80, pol_0); };
    //
    obsThFactory["mueeZHZga550_p80_m30"] = [=](const StandardModel& SM) { return new mueeZHZga(SM, sqrt_s_leptcoll_550, pol_80, -pol_30); };
    obsThFactory["mueeZHZga550_m80_p30"] = [=](const StandardModel& SM) { return new mueeZHZga(SM, sqrt_s_leptcoll_550, -pol_80, pol_30); };
    obsThFactory["mueeZHZga550_p80_0"] = [=](const StandardModel& SM) { return new mueeZHZga(SM, sqrt_s_leptcoll_550, pol_80, pol_0); };
    obsThFactory["mueeZHZga550_m80_0"] = [=](const StandardModel& SM) { return new mueeZHZga(SM, sqrt_s_leptcoll_550, -pol_80, pol_0); };
    //
    obsThFactory["mueeZHZga1000_p80_m30"] = [=](const StandardModel& SM) { return new mueeZHZga(SM, sqrt_s_leptcoll_1000, pol_80, -pol_30); };
    obsThFactory["mueeZHZga1000_m80_p30"] = [=](const StandardModel& SM) { return new mueeZHZga(SM, sqrt_s_leptcoll_1000, -pol_80, pol_30); };
    obsThFactory["mueeZHZga1000_p80_m20"] = [=](const StandardModel& SM) { return new mueeZHZga(SM, sqrt_s_leptcoll_1000, pol_80, -pol_20); };
    obsThFactory["mueeZHZga1000_m80_p20"] = [=](const StandardModel& SM) { return new mueeZHZga(SM, sqrt_s_leptcoll_1000, -pol_80, pol_20); };
    obsThFactory["mueeZHZga1000_p80_0"] = [=](const StandardModel& SM) { return new mueeZHZga(SM, sqrt_s_leptcoll_1000, pol_80, pol_0); };
    obsThFactory["mueeZHZga1000_m80_0"] = [=](const StandardModel& SM) { return new mueeZHZga(SM, sqrt_s_leptcoll_1000, -pol_80, pol_0); };
    //
    obsThFactory["mueeZHZga1400_p80_m30"] = [=](const StandardModel& SM) { return new mueeZHZga(SM, sqrt_s_leptcoll_1400, pol_80, -pol_30); };
    obsThFactory["mueeZHZga1400_m80_p30"] = [=](const StandardModel& SM) { return new mueeZHZga(SM, sqrt_s_leptcoll_1400, -pol_80, pol_30); };
    obsThFactory["mueeZHZga1400_p80_0"] = [=](const StandardModel& SM) { return new mueeZHZga(SM, sqrt_s_leptcoll_1400, pol_80, pol_0); };
    obsThFactory["mueeZHZga1400_m80_0"] = [=](const StandardModel& SM) { return new mueeZHZga(SM, sqrt_s_leptcoll_1400, -pol_80, pol_0); };
    //
    obsThFactory["mueeZHZga1500_p80_m30"] = [=](const StandardModel& SM) { return new mueeZHZga(SM, sqrt_s_leptcoll_1500, pol_80, -pol_30); };
    obsThFactory["mueeZHZga1500_m80_p30"] = [=](const StandardModel& SM) { return new mueeZHZga(SM, sqrt_s_leptcoll_1500, -pol_80, pol_30); };
    obsThFactory["mueeZHZga1500_p80_0"] = [=](const StandardModel& SM) { return new mueeZHZga(SM, sqrt_s_leptcoll_1500, pol_80, pol_0); };
    obsThFactory["mueeZHZga1500_m80_0"] = [=](const StandardModel& SM) { return new mueeZHZga(SM, sqrt_s_leptcoll_1500, -pol_80, pol_0); };
    //
    obsThFactory["mueeZHZga3000_p80_m30"] = [=](const StandardModel& SM) { return new mueeZHZga(SM, sqrt_s_leptcoll_3000, pol_80, -pol_30); };
    obsThFactory["mueeZHZga3000_m80_p30"] = [=](const StandardModel& SM) { return new mueeZHZga(SM, sqrt_s_leptcoll_3000, -pol_80, pol_30); };
    obsThFactory["mueeZHZga3000_p80_0"] = [=](const StandardModel& SM) { return new mueeZHZga(SM, sqrt_s_leptcoll_3000, pol_80, pol_0); };
    obsThFactory["mueeZHZga3000_m80_0"] = [=](const StandardModel& SM) { return new mueeZHZga(SM, sqrt_s_leptcoll_3000, -pol_80, pol_0); };
    //
    // H -> ga ga
    obsThFactory["mueeZHgaga230"] = [=](const StandardModel& SM) { return new mueeZHgaga(SM, sqrt_s_leptcoll_230, 0., 0.); };
    obsThFactory["mueeZHgaga240"] = [=](const StandardModel& SM) { return new mueeZHgaga(SM, sqrt_s_leptcoll_240, 0., 0.); };
    obsThFactory["mueeZHgaga250"] = [=](const StandardModel& SM) { return new mueeZHgaga(SM, sqrt_s_leptcoll_250, 0., 0.); };
    obsThFactory["mueeZHgaga350"] = [=](const StandardModel& SM) { return new mueeZHgaga(SM, sqrt_s_leptcoll_350, 0., 0.); };
    obsThFactory["mueeZHgaga365"] = [=](const StandardModel& SM) { return new mueeZHgaga(SM, sqrt_s_leptcoll_365, 0., 0.); };
    obsThFactory["mueeZHgaga380"] = [=](const StandardModel& SM) { return new mueeZHgaga(SM, sqrt_s_leptcoll_380, 0., 0.); };
    obsThFactory["mueeZHgaga500"] = [=](const StandardModel& SM) { return new mueeZHgaga(SM, sqrt_s_leptcoll_500, 0., 0.); };
    obsThFactory["mueeZHgaga550"] = [=](const StandardModel& SM) { return new mueeZHgaga(SM, sqrt_s_leptcoll_550, 0., 0.); };
    obsThFactory["mueeZHgaga1000"] = [=](const StandardModel& SM) { return new mueeZHgaga(SM, sqrt_s_leptcoll_1000, 0., 0.); };
    obsThFactory["mueeZHgaga1400"] = [=](const StandardModel& SM) { return new mueeZHgaga(SM, sqrt_s_leptcoll_1400, 0., 0.); };
    obsThFactory["mueeZHgaga1500"] = [=](const StandardModel& SM) { return new mueeZHgaga(SM, sqrt_s_leptcoll_1500, 0., 0.); };
    obsThFactory["mueeZHgaga3000"] = [=](const StandardModel& SM) { return new mueeZHgaga(SM, sqrt_s_leptcoll_3000, 0., 0.); };
    //
    obsThFactory["mueeZHgaga250_p80_m30"] = [=](const StandardModel& SM) { return new mueeZHgaga(SM, sqrt_s_leptcoll_250, pol_80, -pol_30); };
    obsThFactory["mueeZHgaga250_m80_p30"] = [=](const StandardModel& SM) { return new mueeZHgaga(SM, sqrt_s_leptcoll_250, -pol_80, pol_30); };
    obsThFactory["mueeZHgaga250_p80_0"] = [=](const StandardModel& SM) { return new mueeZHgaga(SM, sqrt_s_leptcoll_250, pol_80, pol_0); };
    obsThFactory["mueeZHgaga250_m80_0"] = [=](const StandardModel& SM) { return new mueeZHgaga(SM, sqrt_s_leptcoll_250, -pol_80, pol_0); };
    //
    obsThFactory["mueeZHgaga350_p80_m30"] = [=](const StandardModel& SM) { return new mueeZHgaga(SM, sqrt_s_leptcoll_350, pol_80, -pol_30); };
    obsThFactory["mueeZHgaga350_m80_p30"] = [=](const StandardModel& SM) { return new mueeZHgaga(SM, sqrt_s_leptcoll_350, -pol_80, pol_30); };
    obsThFactory["mueeZHgaga350_p80_0"] = [=](const StandardModel& SM) { return new mueeZHgaga(SM, sqrt_s_leptcoll_350, pol_80, pol_0); };
    obsThFactory["mueeZHgaga350_m80_0"] = [=](const StandardModel& SM) { return new mueeZHgaga(SM, sqrt_s_leptcoll_350, -pol_80, pol_0); };
    //
    obsThFactory["mueeZHgaga380_p80_m30"] = [=](const StandardModel& SM) { return new mueeZHgaga(SM, sqrt_s_leptcoll_380, pol_80, -pol_30); };
    obsThFactory["mueeZHgaga380_m80_p30"] = [=](const StandardModel& SM) { return new mueeZHgaga(SM, sqrt_s_leptcoll_380, -pol_80, pol_30); };
    obsThFactory["mueeZHgaga380_p80_0"] = [=](const StandardModel& SM) { return new mueeZHgaga(SM, sqrt_s_leptcoll_380, pol_80, pol_0); };
    obsThFactory["mueeZHgaga380_m80_0"] = [=](const StandardModel& SM) { return new mueeZHgaga(SM, sqrt_s_leptcoll_380, -pol_80, pol_0); };
    //
    obsThFactory["mueeZHgaga500_p80_m30"] = [=](const StandardModel& SM) { return new mueeZHgaga(SM, sqrt_s_leptcoll_500, pol_80, -pol_30); };
    obsThFactory["mueeZHgaga500_m80_p30"] = [=](const StandardModel& SM) { return new mueeZHgaga(SM, sqrt_s_leptcoll_500, -pol_80, pol_30); };
    obsThFactory["mueeZHgaga500_p80_0"] = [=](const StandardModel& SM) { return new mueeZHgaga(SM, sqrt_s_leptcoll_500, pol_80, pol_0); };
    obsThFactory["mueeZHgaga500_m80_0"] = [=](const StandardModel& SM) { return new mueeZHgaga(SM, sqrt_s_leptcoll_500, -pol_80, pol_0); };
    //
    obsThFactory["mueeZHgaga550_p80_m30"] = [=](const StandardModel& SM) { return new mueeZHgaga(SM, sqrt_s_leptcoll_550, pol_80, -pol_30); };
    obsThFactory["mueeZHgaga550_m80_p30"] = [=](const StandardModel& SM) { return new mueeZHgaga(SM, sqrt_s_leptcoll_550, -pol_80, pol_30); };
    obsThFactory["mueeZHgaga550_p80_0"] = [=](const StandardModel& SM) { return new mueeZHgaga(SM, sqrt_s_leptcoll_550, pol_80, pol_0); };
    obsThFactory["mueeZHgaga550_m80_0"] = [=](const StandardModel& SM) { return new mueeZHgaga(SM, sqrt_s_leptcoll_550, -pol_80, pol_0); };
    //
    obsThFactory["mueeZHgaga1000_p80_m30"] = [=](const StandardModel& SM) { return new mueeZHgaga(SM, sqrt_s_leptcoll_1000, pol_80, -pol_30); };
    obsThFactory["mueeZHgaga1000_m80_p30"] = [=](const StandardModel& SM) { return new mueeZHgaga(SM, sqrt_s_leptcoll_1000, -pol_80, pol_30); };
    obsThFactory["mueeZHgaga1000_p80_m20"] = [=](const StandardModel& SM) { return new mueeZHgaga(SM, sqrt_s_leptcoll_1000, pol_80, -pol_20); };
    obsThFactory["mueeZHgaga1000_m80_p20"] = [=](const StandardModel& SM) { return new mueeZHgaga(SM, sqrt_s_leptcoll_1000, -pol_80, pol_20); };
    obsThFactory["mueeZHgaga1000_p80_0"] = [=](const StandardModel& SM) { return new mueeZHgaga(SM, sqrt_s_leptcoll_1000, pol_80, pol_0); };
    obsThFactory["mueeZHgaga1000_m80_0"] = [=](const StandardModel& SM) { return new mueeZHgaga(SM, sqrt_s_leptcoll_1000, -pol_80, pol_0); };
    //
    obsThFactory["mueeZHgaga1400_p80_m30"] = [=](const StandardModel& SM) { return new mueeZHgaga(SM, sqrt_s_leptcoll_1400, pol_80, -pol_30); };
    obsThFactory["mueeZHgaga1400_m80_p30"] = [=](const StandardModel& SM) { return new mueeZHgaga(SM, sqrt_s_leptcoll_1400, -pol_80, pol_30); };
    obsThFactory["mueeZHgaga1400_p80_0"] = [=](const StandardModel& SM) { return new mueeZHgaga(SM, sqrt_s_leptcoll_1400, pol_80, pol_0); };
    obsThFactory["mueeZHgaga1400_m80_0"] = [=](const StandardModel& SM) { return new mueeZHgaga(SM, sqrt_s_leptcoll_1400, -pol_80, pol_0); };
    //
    obsThFactory["mueeZHgaga1500_p80_m30"] = [=](const StandardModel& SM) { return new mueeZHgaga(SM, sqrt_s_leptcoll_1500, pol_80, -pol_30); };
    obsThFactory["mueeZHgaga1500_m80_p30"] = [=](const StandardModel& SM) { return new mueeZHgaga(SM, sqrt_s_leptcoll_1500, -pol_80, pol_30); };
    obsThFactory["mueeZHgaga1500_p80_0"] = [=](const StandardModel& SM) { return new mueeZHgaga(SM, sqrt_s_leptcoll_1500, pol_80, pol_0); };
    obsThFactory["mueeZHgaga1500_m80_0"] = [=](const StandardModel& SM) { return new mueeZHgaga(SM, sqrt_s_leptcoll_1500, -pol_80, pol_0); };
    //
    obsThFactory["mueeZHgaga3000_p80_m30"] = [=](const StandardModel& SM) { return new mueeZHgaga(SM, sqrt_s_leptcoll_3000, pol_80, -pol_30); };
    obsThFactory["mueeZHgaga3000_m80_p30"] = [=](const StandardModel& SM) { return new mueeZHgaga(SM, sqrt_s_leptcoll_3000, -pol_80, pol_30); };
    obsThFactory["mueeZHgaga3000_p80_0"] = [=](const StandardModel& SM) { return new mueeZHgaga(SM, sqrt_s_leptcoll_3000, pol_80, pol_0); };
    obsThFactory["mueeZHgaga3000_m80_0"] = [=](const StandardModel& SM) { return new mueeZHgaga(SM, sqrt_s_leptcoll_3000, -pol_80, pol_0); };
    //
    // H -> mu mu
    obsThFactory["mueeZHmumu230"] = [=](const StandardModel& SM) { return new mueeZHmumu(SM, sqrt_s_leptcoll_230, 0., 0.); };
    obsThFactory["mueeZHmumu240"] = [=](const StandardModel& SM) { return new mueeZHmumu(SM, sqrt_s_leptcoll_240, 0., 0.); };
    obsThFactory["mueeZHmumu250"] = [=](const StandardModel& SM) { return new mueeZHmumu(SM, sqrt_s_leptcoll_250, 0., 0.); };
    obsThFactory["mueeZHmumu350"] = [=](const StandardModel& SM) { return new mueeZHmumu(SM, sqrt_s_leptcoll_350, 0., 0.); };
    obsThFactory["mueeZHmumu365"] = [=](const StandardModel& SM) { return new mueeZHmumu(SM, sqrt_s_leptcoll_365, 0., 0.); };
    obsThFactory["mueeZHmumu380"] = [=](const StandardModel& SM) { return new mueeZHmumu(SM, sqrt_s_leptcoll_380, 0., 0.); };
    obsThFactory["mueeZHmumu500"] = [=](const StandardModel& SM) { return new mueeZHmumu(SM, sqrt_s_leptcoll_500, 0., 0.); };
    obsThFactory["mueeZHmumu550"] = [=](const StandardModel& SM) { return new mueeZHmumu(SM, sqrt_s_leptcoll_550, 0., 0.); };
    obsThFactory["mueeZHmumu1000"] = [=](const StandardModel& SM) { return new mueeZHmumu(SM, sqrt_s_leptcoll_1000, 0., 0.); };
    obsThFactory["mueeZHmumu1400"] = [=](const StandardModel& SM) { return new mueeZHmumu(SM, sqrt_s_leptcoll_1400, 0., 0.); };
    obsThFactory["mueeZHmumu1500"] = [=](const StandardModel& SM) { return new mueeZHmumu(SM, sqrt_s_leptcoll_1500, 0., 0.); };
    obsThFactory["mueeZHmumu3000"] = [=](const StandardModel& SM) { return new mueeZHmumu(SM, sqrt_s_leptcoll_3000, 0., 0.); };
    //
    obsThFactory["mueeZHmumu250_p80_m30"] = [=](const StandardModel& SM) { return new mueeZHmumu(SM, sqrt_s_leptcoll_250, pol_80, -pol_30); };
    obsThFactory["mueeZHmumu250_m80_p30"] = [=](const StandardModel& SM) { return new mueeZHmumu(SM, sqrt_s_leptcoll_250, -pol_80, pol_30); };
    obsThFactory["mueeZHmumu250_p80_0"] = [=](const StandardModel& SM) { return new mueeZHmumu(SM, sqrt_s_leptcoll_250, pol_80, pol_0); };
    obsThFactory["mueeZHmumu250_m80_0"] = [=](const StandardModel& SM) { return new mueeZHmumu(SM, sqrt_s_leptcoll_250, -pol_80, pol_0); };
    //
    obsThFactory["mueeZHmumu350_p80_m30"] = [=](const StandardModel& SM) { return new mueeZHmumu(SM, sqrt_s_leptcoll_350, pol_80, -pol_30); };
    obsThFactory["mueeZHmumu350_m80_p30"] = [=](const StandardModel& SM) { return new mueeZHmumu(SM, sqrt_s_leptcoll_350, -pol_80, pol_30); };
    obsThFactory["mueeZHmumu350_p80_0"] = [=](const StandardModel& SM) { return new mueeZHmumu(SM, sqrt_s_leptcoll_350, pol_80, pol_0); };
    obsThFactory["mueeZHmumu350_m80_0"] = [=](const StandardModel& SM) { return new mueeZHmumu(SM, sqrt_s_leptcoll_350, -pol_80, pol_0); };
    //
    obsThFactory["mueeZHmumu380_p80_m30"] = [=](const StandardModel& SM) { return new mueeZHmumu(SM, sqrt_s_leptcoll_380, pol_80, -pol_30); };
    obsThFactory["mueeZHmumu380_m80_p30"] = [=](const StandardModel& SM) { return new mueeZHmumu(SM, sqrt_s_leptcoll_380, -pol_80, pol_30); };
    obsThFactory["mueeZHmumu380_p80_0"] = [=](const StandardModel& SM) { return new mueeZHmumu(SM, sqrt_s_leptcoll_380, pol_80, pol_0); };
    obsThFactory["mueeZHmumu380_m80_0"] = [=](const StandardModel& SM) { return new mueeZHmumu(SM, sqrt_s_leptcoll_380, -pol_80, pol_0); };
    //
    obsThFactory["mueeZHmumu500_p80_m30"] = [=](const StandardModel& SM) { return new mueeZHmumu(SM, sqrt_s_leptcoll_500, pol_80, -pol_30); };
    obsThFactory["mueeZHmumu500_m80_p30"] = [=](const StandardModel& SM) { return new mueeZHmumu(SM, sqrt_s_leptcoll_500, -pol_80, pol_30); };
    obsThFactory["mueeZHmumu500_p80_0"] = [=](const StandardModel& SM) { return new mueeZHmumu(SM, sqrt_s_leptcoll_500, pol_80, pol_0); };
    obsThFactory["mueeZHmumu500_m80_0"] = [=](const StandardModel& SM) { return new mueeZHmumu(SM, sqrt_s_leptcoll_500, -pol_80, pol_0); };
    //
    obsThFactory["mueeZHmumu550_p80_m30"] = [=](const StandardModel& SM) { return new mueeZHmumu(SM, sqrt_s_leptcoll_550, pol_80, -pol_30); };
    obsThFactory["mueeZHmumu550_m80_p30"] = [=](const StandardModel& SM) { return new mueeZHmumu(SM, sqrt_s_leptcoll_550, -pol_80, pol_30); };
    obsThFactory["mueeZHmumu550_p80_0"] = [=](const StandardModel& SM) { return new mueeZHmumu(SM, sqrt_s_leptcoll_550, pol_80, pol_0); };
    obsThFactory["mueeZHmumu550_m80_0"] = [=](const StandardModel& SM) { return new mueeZHmumu(SM, sqrt_s_leptcoll_550, -pol_80, pol_0); };
    //
    obsThFactory["mueeZHmumu1000_p80_m30"] = [=](const StandardModel& SM) { return new mueeZHmumu(SM, sqrt_s_leptcoll_1000, pol_80, -pol_30); };
    obsThFactory["mueeZHmumu1000_m80_p30"] = [=](const StandardModel& SM) { return new mueeZHmumu(SM, sqrt_s_leptcoll_1000, -pol_80, pol_30); };
    obsThFactory["mueeZHmumu1000_p80_m20"] = [=](const StandardModel& SM) { return new mueeZHmumu(SM, sqrt_s_leptcoll_1000, pol_80, -pol_20); };
    obsThFactory["mueeZHmumu1000_m80_p20"] = [=](const StandardModel& SM) { return new mueeZHmumu(SM, sqrt_s_leptcoll_1000, -pol_80, pol_20); };
    obsThFactory["mueeZHmumu1000_p80_0"] = [=](const StandardModel& SM) { return new mueeZHmumu(SM, sqrt_s_leptcoll_1000, pol_80, pol_0); };
    obsThFactory["mueeZHmumu1000_m80_0"] = [=](const StandardModel& SM) { return new mueeZHmumu(SM, sqrt_s_leptcoll_1000, -pol_80, pol_0); };
    //
    obsThFactory["mueeZHmumu1400_p80_m30"] = [=](const StandardModel& SM) { return new mueeZHmumu(SM, sqrt_s_leptcoll_1400, pol_80, -pol_30); };
    obsThFactory["mueeZHmumu1400_m80_p30"] = [=](const StandardModel& SM) { return new mueeZHmumu(SM, sqrt_s_leptcoll_1400, -pol_80, pol_30); };
    obsThFactory["mueeZHmumu1400_p80_0"] = [=](const StandardModel& SM) { return new mueeZHmumu(SM, sqrt_s_leptcoll_1400, pol_80, pol_0); };
    obsThFactory["mueeZHmumu1400_m80_0"] = [=](const StandardModel& SM) { return new mueeZHmumu(SM, sqrt_s_leptcoll_1400, -pol_80, pol_0); };
    //
    obsThFactory["mueeZHmumu1500_p80_m30"] = [=](const StandardModel& SM) { return new mueeZHmumu(SM, sqrt_s_leptcoll_1500, pol_80, -pol_30); };
    obsThFactory["mueeZHmumu1500_m80_p30"] = [=](const StandardModel& SM) { return new mueeZHmumu(SM, sqrt_s_leptcoll_1500, -pol_80, pol_30); };
    obsThFactory["mueeZHmumu1500_p80_0"] = [=](const StandardModel& SM) { return new mueeZHmumu(SM, sqrt_s_leptcoll_1500, pol_80, pol_0); };
    obsThFactory["mueeZHmumu1500_m80_0"] = [=](const StandardModel& SM) { return new mueeZHmumu(SM, sqrt_s_leptcoll_1500, -pol_80, pol_0); };
    //
    obsThFactory["mueeZHmumu3000_p80_m30"] = [=](const StandardModel& SM) { return new mueeZHmumu(SM, sqrt_s_leptcoll_3000, pol_80, -pol_30); };
    obsThFactory["mueeZHmumu3000_m80_p30"] = [=](const StandardModel& SM) { return new mueeZHmumu(SM, sqrt_s_leptcoll_3000, -pol_80, pol_30); };
    obsThFactory["mueeZHmumu3000_p80_0"] = [=](const StandardModel& SM) { return new mueeZHmumu(SM, sqrt_s_leptcoll_3000, pol_80, pol_0); };
    obsThFactory["mueeZHmumu3000_m80_0"] = [=](const StandardModel& SM) { return new mueeZHmumu(SM, sqrt_s_leptcoll_3000, -pol_80, pol_0); };
    //
    // H -> inv
    obsThFactory["mueeZHinv230"] = [=](const StandardModel& SM) { return new mueeZHinv(SM, sqrt_s_leptcoll_230, 0., 0.); };
    obsThFactory["mueeZHinv240"] = [=](const StandardModel& SM) { return new mueeZHinv(SM, sqrt_s_leptcoll_240, 0., 0.); };
    obsThFactory["mueeZHinv250"] = [=](const StandardModel& SM) { return new mueeZHinv(SM, sqrt_s_leptcoll_250, 0., 0.); };
    obsThFactory["mueeZHinv350"] = [=](const StandardModel& SM) { return new mueeZHinv(SM, sqrt_s_leptcoll_350, 0., 0.); };
    obsThFactory["mueeZHinv365"] = [=](const StandardModel& SM) { return new mueeZHinv(SM, sqrt_s_leptcoll_365, 0., 0.); };
    obsThFactory["mueeZHinv380"] = [=](const StandardModel& SM) { return new mueeZHinv(SM, sqrt_s_leptcoll_380, 0., 0.); };
    obsThFactory["mueeZHinv500"] = [=](const StandardModel& SM) { return new mueeZHinv(SM, sqrt_s_leptcoll_500, 0., 0.); };
    obsThFactory["mueeZHinv550"] = [=](const StandardModel& SM) { return new mueeZHinv(SM, sqrt_s_leptcoll_550, 0., 0.); };
    //
    obsThFactory["mueeZHinv250_p80_m30"] = [=](const StandardModel& SM) { return new mueeZHinv(SM, sqrt_s_leptcoll_250, pol_80, -pol_30); };
    obsThFactory["mueeZHinv250_m80_p30"] = [=](const StandardModel& SM) { return new mueeZHinv(SM, sqrt_s_leptcoll_250, -pol_80, pol_30); };
    obsThFactory["mueeZHinv250_p80_0"] = [=](const StandardModel& SM) { return new mueeZHinv(SM, sqrt_s_leptcoll_250, pol_80, pol_0); };
    obsThFactory["mueeZHinv250_m80_0"] = [=](const StandardModel& SM) { return new mueeZHinv(SM, sqrt_s_leptcoll_250, -pol_80, pol_0); };
    //
    obsThFactory["mueeZHinv350_p80_m30"] = [=](const StandardModel& SM) { return new mueeZHinv(SM, sqrt_s_leptcoll_350, pol_80, -pol_30); };
    obsThFactory["mueeZHinv350_m80_p30"] = [=](const StandardModel& SM) { return new mueeZHinv(SM, sqrt_s_leptcoll_350, -pol_80, pol_30); };
    obsThFactory["mueeZHinv350_p80_0"] = [=](const StandardModel& SM) { return new mueeZHinv(SM, sqrt_s_leptcoll_350, pol_80, pol_0); };
    obsThFactory["mueeZHinv350_m80_0"] = [=](const StandardModel& SM) { return new mueeZHinv(SM, sqrt_s_leptcoll_350, -pol_80, pol_0); };
    //
    obsThFactory["mueeZHinv380_p80_m30"] = [=](const StandardModel& SM) { return new mueeZHinv(SM, sqrt_s_leptcoll_380, pol_80, -pol_30); };
    obsThFactory["mueeZHinv380_m80_p30"] = [=](const StandardModel& SM) { return new mueeZHinv(SM, sqrt_s_leptcoll_380, -pol_80, pol_30); };
    obsThFactory["mueeZHinv380_p80_0"] = [=](const StandardModel& SM) { return new mueeZHinv(SM, sqrt_s_leptcoll_380, pol_80, pol_0); };
    obsThFactory["mueeZHinv380_m80_0"] = [=](const StandardModel& SM) { return new mueeZHinv(SM, sqrt_s_leptcoll_380, -pol_80, pol_0); };
    //
    obsThFactory["mueeZHinv500_p80_m30"] = [=](const StandardModel& SM) { return new mueeZHinv(SM, sqrt_s_leptcoll_500, pol_80, -pol_30); };
    obsThFactory["mueeZHinv500_m80_p30"] = [=](const StandardModel& SM) { return new mueeZHinv(SM, sqrt_s_leptcoll_500, -pol_80, pol_30); };
    obsThFactory["mueeZHinv500_p80_0"] = [=](const StandardModel& SM) { return new mueeZHinv(SM, sqrt_s_leptcoll_500, pol_80, pol_0); };
    obsThFactory["mueeZHinv500_m80_0"] = [=](const StandardModel& SM) { return new mueeZHinv(SM, sqrt_s_leptcoll_500, -pol_80, pol_0); };
    //
    obsThFactory["mueeZHinv550_p80_m30"] = [=](const StandardModel& SM) { return new mueeZHinv(SM, sqrt_s_leptcoll_550, pol_80, -pol_30); };
    obsThFactory["mueeZHinv550_m80_p30"] = [=](const StandardModel& SM) { return new mueeZHinv(SM, sqrt_s_leptcoll_550, -pol_80, pol_30); };
    obsThFactory["mueeZHinv550_p80_0"] = [=](const StandardModel& SM) { return new mueeZHinv(SM, sqrt_s_leptcoll_550, pol_80, pol_0); };
    obsThFactory["mueeZHinv550_m80_0"] = [=](const StandardModel& SM) { return new mueeZHinv(SM, sqrt_s_leptcoll_550, -pol_80, pol_0); };
    //
    obsThFactory["mueeZHinv1000_p80_m30"] = [=](const StandardModel& SM) { return new mueeZHinv(SM, sqrt_s_leptcoll_1000, pol_80, -pol_30); };
    obsThFactory["mueeZHinv1000_m80_p30"] = [=](const StandardModel& SM) { return new mueeZHinv(SM, sqrt_s_leptcoll_1000, -pol_80, pol_30); };
    obsThFactory["mueeZHinv1000_p80_m20"] = [=](const StandardModel& SM) { return new mueeZHinv(SM, sqrt_s_leptcoll_1000, pol_80, -pol_20); };
    obsThFactory["mueeZHinv1000_m80_p20"] = [=](const StandardModel& SM) { return new mueeZHinv(SM, sqrt_s_leptcoll_1000, -pol_80, pol_20); };
    obsThFactory["mueeZHinv1000_p80_0"] = [=](const StandardModel& SM) { return new mueeZHinv(SM, sqrt_s_leptcoll_1000, pol_80, pol_0); };
    obsThFactory["mueeZHinv1000_m80_0"] = [=](const StandardModel& SM) { return new mueeZHinv(SM, sqrt_s_leptcoll_1000, -pol_80, pol_0); };
    //
    obsThFactory["mueeZHinv1400_p80_m30"] = [=](const StandardModel& SM) { return new mueeZHinv(SM, sqrt_s_leptcoll_1400, pol_80, -pol_30); };
    obsThFactory["mueeZHinv1400_m80_p30"] = [=](const StandardModel& SM) { return new mueeZHinv(SM, sqrt_s_leptcoll_1400, -pol_80, pol_30); };
    obsThFactory["mueeZHinv1400_p80_0"] = [=](const StandardModel& SM) { return new mueeZHinv(SM, sqrt_s_leptcoll_1400, pol_80, pol_0); };
    obsThFactory["mueeZHinv1400_m80_0"] = [=](const StandardModel& SM) { return new mueeZHinv(SM, sqrt_s_leptcoll_1400, -pol_80, pol_0); };
    //
    obsThFactory["mueeZHinv1500_p80_m30"] = [=](const StandardModel& SM) { return new mueeZHinv(SM, sqrt_s_leptcoll_1500, pol_80, -pol_30); };
    obsThFactory["mueeZHinv1500_m80_p30"] = [=](const StandardModel& SM) { return new mueeZHinv(SM, sqrt_s_leptcoll_1500, -pol_80, pol_30); };
    obsThFactory["mueeZHinv1500_p80_0"] = [=](const StandardModel& SM) { return new mueeZHinv(SM, sqrt_s_leptcoll_1500, pol_80, pol_0); };
    obsThFactory["mueeZHinv1500_m80_0"] = [=](const StandardModel& SM) { return new mueeZHinv(SM, sqrt_s_leptcoll_1500, -pol_80, pol_0); };
    //
    obsThFactory["mueeZHinv3000_p80_m30"] = [=](const StandardModel& SM) { return new mueeZHinv(SM, sqrt_s_leptcoll_3000, pol_80, -pol_30); };
    obsThFactory["mueeZHinv3000_m80_p30"] = [=](const StandardModel& SM) { return new mueeZHinv(SM, sqrt_s_leptcoll_3000, -pol_80, pol_30); };
    obsThFactory["mueeZHinv3000_p80_0"] = [=](const StandardModel& SM) { return new mueeZHinv(SM, sqrt_s_leptcoll_3000, pol_80, pol_0); };
    obsThFactory["mueeZHinv3000_m80_0"] = [=](const StandardModel& SM) { return new mueeZHinv(SM, sqrt_s_leptcoll_3000, -pol_80, pol_0); };
    //
    obsThFactory["mueeZHBRinv230"] = [=](const StandardModel& SM) { return new mueeZHBRinv(SM, sqrt_s_leptcoll_230, 0., 0.); };
    obsThFactory["mueeZHBRinv240"] = [=](const StandardModel& SM) { return new mueeZHBRinv(SM, sqrt_s_leptcoll_240, 0., 0.); };
    obsThFactory["mueeZHBRinv250"] = [=](const StandardModel& SM) { return new mueeZHBRinv(SM, sqrt_s_leptcoll_250, 0., 0.); };
    obsThFactory["mueeZHBRinv350"] = [=](const StandardModel& SM) { return new mueeZHBRinv(SM, sqrt_s_leptcoll_350, 0., 0.); };
    obsThFactory["mueeZHBRinv365"] = [=](const StandardModel& SM) { return new mueeZHBRinv(SM, sqrt_s_leptcoll_365, 0., 0.); };
    obsThFactory["mueeZHBRinv380"] = [=](const StandardModel& SM) { return new mueeZHBRinv(SM, sqrt_s_leptcoll_380, 0., 0.); };
    obsThFactory["mueeZHBRinv500"] = [=](const StandardModel& SM) { return new mueeZHBRinv(SM, sqrt_s_leptcoll_500, 0., 0.); };
    obsThFactory["mueeZHBRinv550"] = [=](const StandardModel& SM) { return new mueeZHBRinv(SM, sqrt_s_leptcoll_550, 0., 0.); };
    //
    obsThFactory["mueeZHBRinv250_p80_m30"] = [=](const StandardModel& SM) { return new mueeZHBRinv(SM, sqrt_s_leptcoll_250, pol_80, -pol_30); };
    obsThFactory["mueeZHBRinv250_m80_p30"] = [=](const StandardModel& SM) { return new mueeZHBRinv(SM, sqrt_s_leptcoll_250, -pol_80, pol_30); };
    obsThFactory["mueeZHBRinv250_p80_0"] = [=](const StandardModel& SM) { return new mueeZHBRinv(SM, sqrt_s_leptcoll_250, pol_80, pol_0); };
    obsThFactory["mueeZHBRinv250_m80_0"] = [=](const StandardModel& SM) { return new mueeZHBRinv(SM, sqrt_s_leptcoll_250, -pol_80, pol_0); };
    //
    obsThFactory["mueeZHBRinv350_p80_m30"] = [=](const StandardModel& SM) { return new mueeZHBRinv(SM, sqrt_s_leptcoll_350, pol_80, -pol_30); };
    obsThFactory["mueeZHBRinv350_m80_p30"] = [=](const StandardModel& SM) { return new mueeZHBRinv(SM, sqrt_s_leptcoll_350, -pol_80, pol_30); };
    obsThFactory["mueeZHBRinv350_p80_0"] = [=](const StandardModel& SM) { return new mueeZHBRinv(SM, sqrt_s_leptcoll_350, pol_80, pol_0); };
    obsThFactory["mueeZHBRinv350_m80_0"] = [=](const StandardModel& SM) { return new mueeZHBRinv(SM, sqrt_s_leptcoll_350, -pol_80, pol_0); };
    //
    obsThFactory["mueeZHBRinv380_p80_m30"] = [=](const StandardModel& SM) { return new mueeZHBRinv(SM, sqrt_s_leptcoll_380, pol_80, -pol_30); };
    obsThFactory["mueeZHBRinv380_m80_p30"] = [=](const StandardModel& SM) { return new mueeZHBRinv(SM, sqrt_s_leptcoll_380, -pol_80, pol_30); };
    obsThFactory["mueeZHBRinv380_p80_0"] = [=](const StandardModel& SM) { return new mueeZHBRinv(SM, sqrt_s_leptcoll_380, pol_80, pol_0); };
    obsThFactory["mueeZHBRinv380_m80_0"] = [=](const StandardModel& SM) { return new mueeZHBRinv(SM, sqrt_s_leptcoll_380, -pol_80, pol_0); };
    //
    obsThFactory["mueeZHBRinv500_p80_m30"] = [=](const StandardModel& SM) { return new mueeZHBRinv(SM, sqrt_s_leptcoll_500, pol_80, -pol_30); };
    obsThFactory["mueeZHBRinv500_m80_p30"] = [=](const StandardModel& SM) { return new mueeZHBRinv(SM, sqrt_s_leptcoll_500, -pol_80, pol_30); };
    obsThFactory["mueeZHBRinv500_p80_0"] = [=](const StandardModel& SM) { return new mueeZHBRinv(SM, sqrt_s_leptcoll_500, pol_80, pol_0); };
    obsThFactory["mueeZHBRinv500_m80_0"] = [=](const StandardModel& SM) { return new mueeZHBRinv(SM, sqrt_s_leptcoll_500, -pol_80, pol_0); };
    //
    obsThFactory["mueeZHBRinv550_p80_m30"] = [=](const StandardModel& SM) { return new mueeZHBRinv(SM, sqrt_s_leptcoll_550, pol_80, -pol_30); };
    obsThFactory["mueeZHBRinv550_m80_p30"] = [=](const StandardModel& SM) { return new mueeZHBRinv(SM, sqrt_s_leptcoll_550, -pol_80, pol_30); };
    obsThFactory["mueeZHBRinv550_p80_0"] = [=](const StandardModel& SM) { return new mueeZHBRinv(SM, sqrt_s_leptcoll_550, pol_80, pol_0); };
    obsThFactory["mueeZHBRinv550_m80_0"] = [=](const StandardModel& SM) { return new mueeZHBRinv(SM, sqrt_s_leptcoll_550, -pol_80, pol_0); };
    //
    obsThFactory["mueeZHBRinv1000_p80_m30"] = [=](const StandardModel& SM) { return new mueeZHBRinv(SM, sqrt_s_leptcoll_1000, pol_80, -pol_30); };
    obsThFactory["mueeZHBRinv1000_m80_p30"] = [=](const StandardModel& SM) { return new mueeZHBRinv(SM, sqrt_s_leptcoll_1000, -pol_80, pol_30); };
    obsThFactory["mueeZHBRinv1000_p80_m20"] = [=](const StandardModel& SM) { return new mueeZHBRinv(SM, sqrt_s_leptcoll_1000, pol_80, -pol_20); };
    obsThFactory["mueeZHBRinv1000_m80_p20"] = [=](const StandardModel& SM) { return new mueeZHBRinv(SM, sqrt_s_leptcoll_1000, -pol_80, pol_20); };
    obsThFactory["mueeZHBRinv1000_p80_0"] = [=](const StandardModel& SM) { return new mueeZHBRinv(SM, sqrt_s_leptcoll_1000, pol_80, pol_0); };
    obsThFactory["mueeZHBRinv1000_m80_0"] = [=](const StandardModel& SM) { return new mueeZHBRinv(SM, sqrt_s_leptcoll_1000, -pol_80, pol_0); };
    //
    obsThFactory["mueeZHBRinv1400_p80_m30"] = [=](const StandardModel& SM) { return new mueeZHBRinv(SM, sqrt_s_leptcoll_1400, pol_80, -pol_30); };
    obsThFactory["mueeZHBRinv1400_m80_p30"] = [=](const StandardModel& SM) { return new mueeZHBRinv(SM, sqrt_s_leptcoll_1400, -pol_80, pol_30); };
    obsThFactory["mueeZHBRinv1400_p80_0"] = [=](const StandardModel& SM) { return new mueeZHBRinv(SM, sqrt_s_leptcoll_1400, pol_80, pol_0); };
    obsThFactory["mueeZHBRinv1400_m80_0"] = [=](const StandardModel& SM) { return new mueeZHBRinv(SM, sqrt_s_leptcoll_1400, -pol_80, pol_0); };
    //
    obsThFactory["mueeZHBRinv1500_p80_m30"] = [=](const StandardModel& SM) { return new mueeZHBRinv(SM, sqrt_s_leptcoll_1500, pol_80, -pol_30); };
    obsThFactory["mueeZHBRinv1500_m80_p30"] = [=](const StandardModel& SM) { return new mueeZHBRinv(SM, sqrt_s_leptcoll_1500, -pol_80, pol_30); };
    obsThFactory["mueeZHBRinv1500_p80_0"] = [=](const StandardModel& SM) { return new mueeZHBRinv(SM, sqrt_s_leptcoll_1500, pol_80, pol_0); };
    obsThFactory["mueeZHBRinv1500_m80_0"] = [=](const StandardModel& SM) { return new mueeZHBRinv(SM, sqrt_s_leptcoll_1500, -pol_80, pol_0); };
    //
    obsThFactory["mueeZHBRinv3000_p80_m30"] = [=](const StandardModel& SM) { return new mueeZHBRinv(SM, sqrt_s_leptcoll_3000, pol_80, -pol_30); };
    obsThFactory["mueeZHBRinv3000_m80_p30"] = [=](const StandardModel& SM) { return new mueeZHBRinv(SM, sqrt_s_leptcoll_3000, -pol_80, pol_30); };
    obsThFactory["mueeZHBRinv3000_p80_0"] = [=](const StandardModel& SM) { return new mueeZHBRinv(SM, sqrt_s_leptcoll_3000, pol_80, pol_0); };
    obsThFactory["mueeZHBRinv3000_m80_0"] = [=](const StandardModel& SM) { return new mueeZHBRinv(SM, sqrt_s_leptcoll_3000, -pol_80, pol_0); };
    //
    // ZBF
    obsThFactory["mueeZBFbb230"] = [=](const StandardModel& SM) { return new mueeZBFbb(SM, sqrt_s_leptcoll_230, 0., 0.); };
    obsThFactory["mueeZBFbb240"] = [=](const StandardModel& SM) { return new mueeZBFbb(SM, sqrt_s_leptcoll_240, 0., 0.); };
    obsThFactory["mueeZBFbb250"] = [=](const StandardModel& SM) { return new mueeZBFbb(SM, sqrt_s_leptcoll_250, 0., 0.); };
    obsThFactory["mueeZBFbb350"] = [=](const StandardModel& SM) { return new mueeZBFbb(SM, sqrt_s_leptcoll_350, 0., 0.); };
    obsThFactory["mueeZBFbb365"] = [=](const StandardModel& SM) { return new mueeZBFbb(SM, sqrt_s_leptcoll_365, 0., 0.); };
    obsThFactory["mueeZBFbb380"] = [=](const StandardModel& SM) { return new mueeZBFbb(SM, sqrt_s_leptcoll_380, 0., 0.); };
    obsThFactory["mueeZBFbb500"] = [=](const StandardModel& SM) { return new mueeZBFbb(SM, sqrt_s_leptcoll_500, 0., 0.); };
    obsThFactory["mueeZBFbb1000"] = [=](const StandardModel& SM) { return new mueeZBFbb(SM, sqrt_s_leptcoll_1000, 0., 0.); };
    obsThFactory["mueeZBFbb1400"] = [=](const StandardModel& SM) { return new mueeZBFbb(SM, sqrt_s_leptcoll_1400, 0., 0.); };
    obsThFactory["mueeZBFbb1500"] = [=](const StandardModel& SM) { return new mueeZBFbb(SM, sqrt_s_leptcoll_1500, 0., 0.); };
    obsThFactory["mueeZBFbb3000"] = [=](const StandardModel& SM) { return new mueeZBFbb(SM, sqrt_s_leptcoll_3000, 0., 0.); };
    //
    obsThFactory["mueeZBFbb250_p80_m30"] = [=](const StandardModel& SM) { return new mueeZBFbb(SM, sqrt_s_leptcoll_250, pol_80, -pol_30); };
    obsThFactory["mueeZBFbb250_m80_p30"] = [=](const StandardModel& SM) { return new mueeZBFbb(SM, sqrt_s_leptcoll_250, -pol_80, pol_30); };
    obsThFactory["mueeZBFbb250_p80_0"] = [=](const StandardModel& SM) { return new mueeZBFbb(SM, sqrt_s_leptcoll_250, pol_80, pol_0); };
    obsThFactory["mueeZBFbb250_m80_0"] = [=](const StandardModel& SM) { return new mueeZBFbb(SM, sqrt_s_leptcoll_250, -pol_80, pol_0); };
    //
    obsThFactory["mueeZBFbb350_p80_m30"] = [=](const StandardModel& SM) { return new mueeZBFbb(SM, sqrt_s_leptcoll_350, pol_80, -pol_30); };
    obsThFactory["mueeZBFbb350_m80_p30"] = [=](const StandardModel& SM) { return new mueeZBFbb(SM, sqrt_s_leptcoll_350, -pol_80, pol_30); };
    obsThFactory["mueeZBFbb350_p80_0"] = [=](const StandardModel& SM) { return new mueeZBFbb(SM, sqrt_s_leptcoll_350, pol_80, pol_0); };
    obsThFactory["mueeZBFbb350_m80_0"] = [=](const StandardModel& SM) { return new mueeZBFbb(SM, sqrt_s_leptcoll_350, -pol_80, pol_0); };
    //
    obsThFactory["mueeZBFbb380_p80_m30"] = [=](const StandardModel& SM) { return new mueeZBFbb(SM, sqrt_s_leptcoll_380, pol_80, -pol_30); };
    obsThFactory["mueeZBFbb380_m80_p30"] = [=](const StandardModel& SM) { return new mueeZBFbb(SM, sqrt_s_leptcoll_380, -pol_80, pol_30); };
    obsThFactory["mueeZBFbb380_p80_0"] = [=](const StandardModel& SM) { return new mueeZBFbb(SM, sqrt_s_leptcoll_380, pol_80, pol_0); };
    obsThFactory["mueeZBFbb380_m80_0"] = [=](const StandardModel& SM) { return new mueeZBFbb(SM, sqrt_s_leptcoll_380, -pol_80, pol_0); };
    //
    obsThFactory["mueeZBFbb500_p80_m30"] = [=](const StandardModel& SM) { return new mueeZBFbb(SM, sqrt_s_leptcoll_500, pol_80, -pol_30); };
    obsThFactory["mueeZBFbb500_m80_p30"] = [=](const StandardModel& SM) { return new mueeZBFbb(SM, sqrt_s_leptcoll_500, -pol_80, pol_30); };
    obsThFactory["mueeZBFbb500_p80_0"] = [=](const StandardModel& SM) { return new mueeZBFbb(SM, sqrt_s_leptcoll_500, pol_80, pol_0); };
    obsThFactory["mueeZBFbb500_m80_0"] = [=](const StandardModel& SM) { return new mueeZBFbb(SM, sqrt_s_leptcoll_500, -pol_80, pol_0); };
    //
    obsThFactory["mueeZBFbb1000_p80_m30"] = [=](const StandardModel& SM) { return new mueeZBFbb(SM, sqrt_s_leptcoll_1000, pol_80, -pol_30); };
    obsThFactory["mueeZBFbb1000_m80_p30"] = [=](const StandardModel& SM) { return new mueeZBFbb(SM, sqrt_s_leptcoll_1000, -pol_80, pol_30); };
    obsThFactory["mueeZBFbb1000_p80_m20"] = [=](const StandardModel& SM) { return new mueeZBFbb(SM, sqrt_s_leptcoll_1000, pol_80, -pol_20); };
    obsThFactory["mueeZBFbb1000_m80_p20"] = [=](const StandardModel& SM) { return new mueeZBFbb(SM, sqrt_s_leptcoll_1000, -pol_80, pol_20); };
    obsThFactory["mueeZBFbb1000_p80_0"] = [=](const StandardModel& SM) { return new mueeZBFbb(SM, sqrt_s_leptcoll_1000, pol_80, pol_0); };
    obsThFactory["mueeZBFbb1000_m80_0"] = [=](const StandardModel& SM) { return new mueeZBFbb(SM, sqrt_s_leptcoll_1000, -pol_80, pol_0); };
    //
    obsThFactory["mueeZBFbb1400_p80_m30"] = [=](const StandardModel& SM) { return new mueeZBFbb(SM, sqrt_s_leptcoll_1400, pol_80, -pol_30); };
    obsThFactory["mueeZBFbb1400_m80_p30"] = [=](const StandardModel& SM) { return new mueeZBFbb(SM, sqrt_s_leptcoll_1400, -pol_80, pol_30); };
    obsThFactory["mueeZBFbb1400_p80_0"] = [=](const StandardModel& SM) { return new mueeZBFbb(SM, sqrt_s_leptcoll_1400, pol_80, pol_0); };
    obsThFactory["mueeZBFbb1400_m80_0"] = [=](const StandardModel& SM) { return new mueeZBFbb(SM, sqrt_s_leptcoll_1400, -pol_80, pol_0); };
    //
    obsThFactory["mueeZBFbb1500_p80_m30"] = [=](const StandardModel& SM) { return new mueeZBFbb(SM, sqrt_s_leptcoll_1500, pol_80, -pol_30); };
    obsThFactory["mueeZBFbb1500_m80_p30"] = [=](const StandardModel& SM) { return new mueeZBFbb(SM, sqrt_s_leptcoll_1500, -pol_80, pol_30); };
    obsThFactory["mueeZBFbb1500_p80_0"] = [=](const StandardModel& SM) { return new mueeZBFbb(SM, sqrt_s_leptcoll_1500, pol_80, pol_0); };
    obsThFactory["mueeZBFbb1500_m80_0"] = [=](const StandardModel& SM) { return new mueeZBFbb(SM, sqrt_s_leptcoll_1500, -pol_80, pol_0); };
    //
    obsThFactory["mueeZBFbb3000_p80_m30"] = [=](const StandardModel& SM) { return new mueeZBFbb(SM, sqrt_s_leptcoll_3000, pol_80, -pol_30); };
    obsThFactory["mueeZBFbb3000_m80_p30"] = [=](const StandardModel& SM) { return new mueeZBFbb(SM, sqrt_s_leptcoll_3000, -pol_80, pol_30); };
    obsThFactory["mueeZBFbb3000_p80_0"] = [=](const StandardModel& SM) { return new mueeZBFbb(SM, sqrt_s_leptcoll_3000, pol_80, pol_0); };
    obsThFactory["mueeZBFbb3000_m80_0"] = [=](const StandardModel& SM) { return new mueeZBFbb(SM, sqrt_s_leptcoll_3000, -pol_80, pol_0); };
    //
    // eettH
    obsThFactory["mueettHbb500"] = [=](const StandardModel& SM) { return new mueettHbb(SM, sqrt_s_leptcoll_500, 0., 0.); };
    obsThFactory["mueettHbb550"] = [=](const StandardModel& SM) { return new mueettHbb(SM, sqrt_s_leptcoll_550, 0., 0.); };
    obsThFactory["mueettHbb1000"] = [=](const StandardModel& SM) { return new mueettHbb(SM, sqrt_s_leptcoll_1000, 0., 0.); };
    obsThFactory["mueettHbb1400"] = [=](const StandardModel& SM) { return new mueettHbb(SM, sqrt_s_leptcoll_1400, 0., 0.); };
    obsThFactory["mueettHbb1500"] = [=](const StandardModel& SM) { return new mueettHbb(SM, sqrt_s_leptcoll_1500, 0., 0.); };
    obsThFactory["mueettHbb3000"] = [=](const StandardModel& SM) { return new mueettHbb(SM, sqrt_s_leptcoll_3000, 0., 0.); };
    //
    obsThFactory["mueettHbb500_p80_m30"] = [=](const StandardModel& SM) { return new mueettHbb(SM, sqrt_s_leptcoll_500, pol_80, -pol_30); };
    obsThFactory["mueettHbb500_m80_p30"] = [=](const StandardModel& SM) { return new mueettHbb(SM, sqrt_s_leptcoll_500, -pol_80, pol_30); };
    obsThFactory["mueettHbb500_p80_0"] = [=](const StandardModel& SM) { return new mueettHbb(SM, sqrt_s_leptcoll_500, pol_80, pol_0); };
    obsThFactory["mueettHbb500_m80_0"] = [=](const StandardModel& SM) { return new mueettHbb(SM, sqrt_s_leptcoll_500, -pol_80, pol_0); };
    //
    obsThFactory["mueettHbb550_p80_m30"] = [=](const StandardModel& SM) { return new mueettHbb(SM, sqrt_s_leptcoll_550, pol_80, -pol_30); };
    obsThFactory["mueettHbb550_m80_p30"] = [=](const StandardModel& SM) { return new mueettHbb(SM, sqrt_s_leptcoll_550, -pol_80, pol_30); };
    obsThFactory["mueettHbb550_p80_0"] = [=](const StandardModel& SM) { return new mueettHbb(SM, sqrt_s_leptcoll_550, pol_80, pol_0); };
    obsThFactory["mueettHbb550_m80_0"] = [=](const StandardModel& SM) { return new mueettHbb(SM, sqrt_s_leptcoll_550, -pol_80, pol_0); };
    //
    obsThFactory["mueettHbb1000_p80_m30"] = [=](const StandardModel& SM) { return new mueettHbb(SM, sqrt_s_leptcoll_1000, pol_80, -pol_30); };
    obsThFactory["mueettHbb1000_m80_p30"] = [=](const StandardModel& SM) { return new mueettHbb(SM, sqrt_s_leptcoll_1000, -pol_80, pol_30); };
    obsThFactory["mueettHbb1000_p80_m20"] = [=](const StandardModel& SM) { return new mueettHbb(SM, sqrt_s_leptcoll_1000, pol_80, -pol_20); };
    obsThFactory["mueettHbb1000_m80_p20"] = [=](const StandardModel& SM) { return new mueettHbb(SM, sqrt_s_leptcoll_1000, -pol_80, pol_20); };
    obsThFactory["mueettHbb1000_p80_0"] = [=](const StandardModel& SM) { return new mueettHbb(SM, sqrt_s_leptcoll_1000, pol_80, pol_0); };
    obsThFactory["mueettHbb1000_m80_0"] = [=](const StandardModel& SM) { return new mueettHbb(SM, sqrt_s_leptcoll_1000, -pol_80, pol_0); };
    //
    obsThFactory["mueettHbb1400_p80_m30"] = [=](const StandardModel& SM) { return new mueettHbb(SM, sqrt_s_leptcoll_1400, pol_80, -pol_30); };
    obsThFactory["mueettHbb1400_m80_p30"] = [=](const StandardModel& SM) { return new mueettHbb(SM, sqrt_s_leptcoll_1400, -pol_80, pol_30); };
    obsThFactory["mueettHbb1400_p80_0"] = [=](const StandardModel& SM) { return new mueettHbb(SM, sqrt_s_leptcoll_1400, pol_80, pol_0); };
    obsThFactory["mueettHbb1400_m80_0"] = [=](const StandardModel& SM) { return new mueettHbb(SM, sqrt_s_leptcoll_1400, -pol_80, pol_0); };
    //
    obsThFactory["mueettHbb1500_p80_m30"] = [=](const StandardModel& SM) { return new mueettHbb(SM, sqrt_s_leptcoll_1500, pol_80, -pol_30); };
    obsThFactory["mueettHbb1500_m80_p30"] = [=](const StandardModel& SM) { return new mueettHbb(SM, sqrt_s_leptcoll_1500, -pol_80, pol_30); };
    obsThFactory["mueettHbb1500_p80_0"] = [=](const StandardModel& SM) { return new mueettHbb(SM, sqrt_s_leptcoll_1500, pol_80, pol_0); };
    obsThFactory["mueettHbb1500_m80_0"] = [=](const StandardModel& SM) { return new mueettHbb(SM, sqrt_s_leptcoll_1500, -pol_80, pol_0); };
    //
    obsThFactory["mueettHbb3000_p80_m30"] = [=](const StandardModel& SM) { return new mueettHbb(SM, sqrt_s_leptcoll_3000, pol_80, -pol_30); };
    obsThFactory["mueettHbb3000_m80_p30"] = [=](const StandardModel& SM) { return new mueettHbb(SM, sqrt_s_leptcoll_3000, -pol_80, pol_30); };
    obsThFactory["mueettHbb3000_p80_0"] = [=](const StandardModel& SM) { return new mueettHbb(SM, sqrt_s_leptcoll_3000, pol_80, pol_0); };
    obsThFactory["mueettHbb3000_m80_0"] = [=](const StandardModel& SM) { return new mueettHbb(SM, sqrt_s_leptcoll_3000, -pol_80, pol_0); };
    //
    //-----  Full Signal strengths per prod and decay: mu+ mu- colliders  ----------
    //
    obsThFactory["mumumuHbb125"] = [=](const StandardModel& SM) { return new mummHbb(SM, sqrt_s_leptcoll_125); };
    obsThFactory["mumumuHcc125"] = [=](const StandardModel& SM) { return new mummHcc(SM, sqrt_s_leptcoll_125); };
    obsThFactory["mumumuHgg125"] = [=](const StandardModel& SM) { return new mummHgg(SM, sqrt_s_leptcoll_125); };
    obsThFactory["mumumuHWW125"] = [=](const StandardModel& SM) { return new mummHWW(SM, sqrt_s_leptcoll_125); };
    obsThFactory["mumumuHtautau125"] = [=](const StandardModel& SM) { return new mummHtautau(SM, sqrt_s_leptcoll_125); };
    obsThFactory["mumumuHZZ125"] = [=](const StandardModel& SM) { return new mummHZZ(SM, sqrt_s_leptcoll_125); };
    obsThFactory["mumumuHZga125"] = [=](const StandardModel& SM) { return new mummHZga(SM, sqrt_s_leptcoll_125); };
    obsThFactory["mumumuHgaga125"] = [=](const StandardModel& SM) { return new mummHgaga(SM, sqrt_s_leptcoll_125); };
    obsThFactory["mumumuHmumu125"] = [=](const StandardModel& SM) { return new mummHmumu(SM, sqrt_s_leptcoll_125); };
    // The same in the narrow width approximation
    obsThFactory["mumumuHbbNWA125"] = [=](const StandardModel& SM) { return new mummHbbNWA(SM, sqrt_s_leptcoll_125); };
    obsThFactory["mumumuHccNWA125"] = [=](const StandardModel& SM) { return new mummHccNWA(SM, sqrt_s_leptcoll_125); };
    obsThFactory["mumumuHggNWA125"] = [=](const StandardModel& SM) { return new mummHggNWA(SM, sqrt_s_leptcoll_125); };
    obsThFactory["mumumuHWWNWA125"] = [=](const StandardModel& SM) { return new mummHWWNWA(SM, sqrt_s_leptcoll_125); };
    obsThFactory["mumumuHtautauNWA125"] = [=](const StandardModel& SM) { return new mummHtautauNWA(SM, sqrt_s_leptcoll_125); };
    obsThFactory["mumumuHZZNWA125"] = [=](const StandardModel& SM) { return new mummHZZNWA(SM, sqrt_s_leptcoll_125); };
    obsThFactory["mumumuHZgaNWA125"] = [=](const StandardModel& SM) { return new mummHZgaNWA(SM, sqrt_s_leptcoll_125); };
    obsThFactory["mumumuHgagaNWA125"] = [=](const StandardModel& SM) { return new mummHgagaNWA(SM, sqrt_s_leptcoll_125); };
    obsThFactory["mumumuHmumuNWA125"] = [=](const StandardModel& SM) { return new mummHmumuNWA(SM, sqrt_s_leptcoll_125); };
    //
    // Signal strengths above the pole
    //
    // Hvv
    obsThFactory["mumumuHvvbb3000"] = [=](const StandardModel& SM) { return new mummHvvbb(SM, sqrt_s_leptcoll_3000); };
    obsThFactory["mumumuHvvbb10000"] = [=](const StandardModel& SM) { return new mummHvvbb(SM, sqrt_s_leptcoll_10000); };
    //
    obsThFactory["mumumuHvvcc3000"] = [=](const StandardModel& SM) { return new mummHvvcc(SM, sqrt_s_leptcoll_3000); };
    obsThFactory["mumumuHvvcc10000"] = [=](const StandardModel& SM) { return new mummHvvcc(SM, sqrt_s_leptcoll_10000); };
    //
    obsThFactory["mumumuHvvgg3000"] = [=](const StandardModel& SM) { return new mummHvvgg(SM, sqrt_s_leptcoll_3000); };
    obsThFactory["mumumuHvvgg10000"] = [=](const StandardModel& SM) { return new mummHvvgg(SM, sqrt_s_leptcoll_10000); };
    //
    obsThFactory["mumumuHvvWW3000"] = [=](const StandardModel& SM) { return new mummHvvWW(SM, sqrt_s_leptcoll_3000); };
    obsThFactory["mumumuHvvWW10000"] = [=](const StandardModel& SM) { return new mummHvvWW(SM, sqrt_s_leptcoll_10000); };
    //
    obsThFactory["mumumuHvvtautau3000"] = [=](const StandardModel& SM) { return new mummHvvtautau(SM, sqrt_s_leptcoll_3000); };
    obsThFactory["mumumuHvvtautau10000"] = [=](const StandardModel& SM) { return new mummHvvtautau(SM, sqrt_s_leptcoll_10000); };
    //
    obsThFactory["mumumuHvvZZ3000"] = [=](const StandardModel& SM) { return new mummHvvZZ(SM, sqrt_s_leptcoll_3000); };
    obsThFactory["mumumuHvvZZ10000"] = [=](const StandardModel& SM) { return new mummHvvZZ(SM, sqrt_s_leptcoll_10000); };
    //
    obsThFactory["mumumuHvvZga3000"] = [=](const StandardModel& SM) { return new mummHvvZga(SM, sqrt_s_leptcoll_3000); };
    obsThFactory["mumumuHvvZga10000"] = [=](const StandardModel& SM) { return new mummHvvZga(SM, sqrt_s_leptcoll_10000); };
    //
    obsThFactory["mumumuHvvgaga3000"] = [=](const StandardModel& SM) { return new mummHvvgaga(SM, sqrt_s_leptcoll_3000); };
    obsThFactory["mumumuHvvgaga10000"] = [=](const StandardModel& SM) { return new mummHvvgaga(SM, sqrt_s_leptcoll_10000); };
    //
    obsThFactory["mumumuHvvmumu3000"] = [=](const StandardModel& SM) { return new mummHvvmumu(SM, sqrt_s_leptcoll_3000); };
    obsThFactory["mumumuHvvmumu10000"] = [=](const StandardModel& SM) { return new mummHvvmumu(SM, sqrt_s_leptcoll_10000); };
    //
    // Hmumu
    obsThFactory["mumumuHmumubb3000"] = [=](const StandardModel& SM) { return new mummHmmbb(SM, sqrt_s_leptcoll_3000); };
    obsThFactory["mumumuHmumubb10000"] = [=](const StandardModel& SM) { return new mummHmmbb(SM, sqrt_s_leptcoll_10000); };
    //
    obsThFactory["mumumuHmumucc3000"] = [=](const StandardModel& SM) { return new mummHmmcc(SM, sqrt_s_leptcoll_3000); };
    obsThFactory["mumumuHmumucc10000"] = [=](const StandardModel& SM) { return new mummHmmcc(SM, sqrt_s_leptcoll_10000); };
    //
    obsThFactory["mumumuHmumugg3000"] = [=](const StandardModel& SM) { return new mummHmmgg(SM, sqrt_s_leptcoll_3000); };
    obsThFactory["mumumuHmumugg10000"] = [=](const StandardModel& SM) { return new mummHmmgg(SM, sqrt_s_leptcoll_10000); };
    //
    obsThFactory["mumumuHmumuWW3000"] = [=](const StandardModel& SM) { return new mummHmmWW(SM, sqrt_s_leptcoll_3000); };
    obsThFactory["mumumuHmumuWW10000"] = [=](const StandardModel& SM) { return new mummHmmWW(SM, sqrt_s_leptcoll_10000); };
    //
    obsThFactory["mumumuHmumutautau3000"] = [=](const StandardModel& SM) { return new mummHmmtautau(SM, sqrt_s_leptcoll_3000); };
    obsThFactory["mumumuHmumutautau10000"] = [=](const StandardModel& SM) { return new mummHmmtautau(SM, sqrt_s_leptcoll_10000); };
    //
    obsThFactory["mumumuHmumuZZ3000"] = [=](const StandardModel& SM) { return new mummHmmZZ(SM, sqrt_s_leptcoll_3000); };
    obsThFactory["mumumuHmumuZZ10000"] = [=](const StandardModel& SM) { return new mummHmmZZ(SM, sqrt_s_leptcoll_10000); };
    //
    obsThFactory["mumumuHmumuZga3000"] = [=](const StandardModel& SM) { return new mummHmmZga(SM, sqrt_s_leptcoll_3000); };
    obsThFactory["mumumuHmumuZga10000"] = [=](const StandardModel& SM) { return new mummHmmZga(SM, sqrt_s_leptcoll_10000); };
    //
    obsThFactory["mumumuHmumugaga3000"] = [=](const StandardModel& SM) { return new mummHmmgaga(SM, sqrt_s_leptcoll_3000); };
    obsThFactory["mumumuHmumugaga10000"] = [=](const StandardModel& SM) { return new mummHmmgaga(SM, sqrt_s_leptcoll_10000); };
    //
    obsThFactory["mumumuHmumumumu3000"] = [=](const StandardModel& SM) { return new mummHmmmumu(SM, sqrt_s_leptcoll_3000); };
    obsThFactory["mumumuHmumumumu10000"] = [=](const StandardModel& SM) { return new mummHmmmumu(SM, sqrt_s_leptcoll_10000); };
    //
    // ZH
    obsThFactory["mumumuZHbb3000"] = [=](const StandardModel& SM) { return new mummZHbb(SM, sqrt_s_leptcoll_3000); };
    obsThFactory["mumumuZHbb10000"] = [=](const StandardModel& SM) { return new mummZHbb(SM, sqrt_s_leptcoll_10000); };
    //
    obsThFactory["mumumuZHcc3000"] = [=](const StandardModel& SM) { return new mummZHcc(SM, sqrt_s_leptcoll_3000); };
    obsThFactory["mumumuZHcc10000"] = [=](const StandardModel& SM) { return new mummZHcc(SM, sqrt_s_leptcoll_10000); };
    //
    obsThFactory["mumumuZHgg3000"] = [=](const StandardModel& SM) { return new mummZHgg(SM, sqrt_s_leptcoll_3000); };
    obsThFactory["mumumuZHgg10000"] = [=](const StandardModel& SM) { return new mummZHgg(SM, sqrt_s_leptcoll_10000); };
    //
    obsThFactory["mumumuZHWW3000"] = [=](const StandardModel& SM) { return new mummZHWW(SM, sqrt_s_leptcoll_3000); };
    obsThFactory["mumumuZHWW10000"] = [=](const StandardModel& SM) { return new mummZHWW(SM, sqrt_s_leptcoll_10000); };
    //
    obsThFactory["mumumuZHtautau3000"] = [=](const StandardModel& SM) { return new mummZHtautau(SM, sqrt_s_leptcoll_3000); };
    obsThFactory["mumumuZHtautau10000"] = [=](const StandardModel& SM) { return new mummZHtautau(SM, sqrt_s_leptcoll_10000); };
    //
    obsThFactory["mumumuZHZZ3000"] = [=](const StandardModel& SM) { return new mummZHZZ(SM, sqrt_s_leptcoll_3000); };
    obsThFactory["mumumuZHZZ10000"] = [=](const StandardModel& SM) { return new mummZHZZ(SM, sqrt_s_leptcoll_10000); };
    //
    obsThFactory["mumumuZHZga3000"] = [=](const StandardModel& SM) { return new mummZHZga(SM, sqrt_s_leptcoll_3000); };
    obsThFactory["mumumuZHZga10000"] = [=](const StandardModel& SM) { return new mummZHZga(SM, sqrt_s_leptcoll_10000); };
    //
    obsThFactory["mumumuZHgaga3000"] = [=](const StandardModel& SM) { return new mummZHgaga(SM, sqrt_s_leptcoll_3000); };
    obsThFactory["mumumuZHgaga10000"] = [=](const StandardModel& SM) { return new mummZHgaga(SM, sqrt_s_leptcoll_10000); };
    //
    obsThFactory["mumumuZHmumu3000"] = [=](const StandardModel& SM) { return new mummZHmumu(SM, sqrt_s_leptcoll_3000); };
    obsThFactory["mumumuZHmumu10000"] = [=](const StandardModel& SM) { return new mummZHmumu(SM, sqrt_s_leptcoll_10000); };
    //
    // mumuttH
    obsThFactory["mumumuttHbb3000"] = [=](const StandardModel& SM) { return new mummttHbb(SM, sqrt_s_leptcoll_3000); };
    obsThFactory["mumumuttHbb10000"] = [=](const StandardModel& SM) { return new mummttHbb(SM, sqrt_s_leptcoll_10000); };
    //
    obsThFactory["mumumuttHcc3000"] = [=](const StandardModel& SM) { return new mummttHcc(SM, sqrt_s_leptcoll_3000); };
    obsThFactory["mumumuttHcc10000"] = [=](const StandardModel& SM) { return new mummttHcc(SM, sqrt_s_leptcoll_10000); };
    //
    obsThFactory["mumumuttHgg3000"] = [=](const StandardModel& SM) { return new mummttHgg(SM, sqrt_s_leptcoll_3000); };
    obsThFactory["mumumuttHgg10000"] = [=](const StandardModel& SM) { return new mummttHgg(SM, sqrt_s_leptcoll_10000); };
    //
    obsThFactory["mumumuttHWW3000"] = [=](const StandardModel& SM) { return new mummttHWW(SM, sqrt_s_leptcoll_3000); };
    obsThFactory["mumumuttHWW10000"] = [=](const StandardModel& SM) { return new mummttHWW(SM, sqrt_s_leptcoll_10000); };
    //
    obsThFactory["mumumuttHtautau3000"] = [=](const StandardModel& SM) { return new mummttHtautau(SM, sqrt_s_leptcoll_3000); };
    obsThFactory["mumumuttHtautau10000"] = [=](const StandardModel& SM) { return new mummttHtautau(SM, sqrt_s_leptcoll_10000); };
    //
    obsThFactory["mumumuttHZZ3000"] = [=](const StandardModel& SM) { return new mummttHZZ(SM, sqrt_s_leptcoll_3000); };
    obsThFactory["mumumuttHZZ10000"] = [=](const StandardModel& SM) { return new mummttHZZ(SM, sqrt_s_leptcoll_10000); };
    //
    obsThFactory["mumumuttHZga3000"] = [=](const StandardModel& SM) { return new mummttHZga(SM, sqrt_s_leptcoll_3000); };
    obsThFactory["mumumuttHZga10000"] = [=](const StandardModel& SM) { return new mummttHZga(SM, sqrt_s_leptcoll_10000); };
    //
    obsThFactory["mumumuttHgaga3000"] = [=](const StandardModel& SM) { return new mummttHgaga(SM, sqrt_s_leptcoll_3000); };
    obsThFactory["mumumuttHgaga10000"] = [=](const StandardModel& SM) { return new mummttHgaga(SM, sqrt_s_leptcoll_10000); };
    //
    obsThFactory["mumumuttHmumu3000"] = [=](const StandardModel& SM) { return new mummttHmumu(SM, sqrt_s_leptcoll_3000); };
    obsThFactory["mumumuttHmumu10000"] = [=](const StandardModel& SM) { return new mummttHmumu(SM, sqrt_s_leptcoll_10000); };
    //
    //-----  Full Signal strengths per prod and decay: Lepton-Hadron colliders  ----------
    //
    obsThFactory["muepWBFbb1200"] = [=](const StandardModel& SM) { return new muepWBFbb(SM, sqrt_s_LHeC_1_2); };
    obsThFactory["muepWBFcc1200"] = [=](const StandardModel& SM) { return new muepWBFcc(SM, sqrt_s_LHeC_1_2); };
    obsThFactory["muepWBFgg1200"] = [=](const StandardModel& SM) { return new muepWBFgg(SM, sqrt_s_LHeC_1_2); };
    obsThFactory["muepWBFWW2l2v1200"] = [=](const StandardModel& SM) { return new muepWBFWW2l2v(SM, sqrt_s_LHeC_1_2); };
    obsThFactory["muepWBFZZ4l1200"] = [=](const StandardModel& SM) { return new muepWBFZZ4l(SM, sqrt_s_LHeC_1_2); };
    obsThFactory["muepWBFgaga1200"] = [=](const StandardModel& SM) { return new muepWBFgaga(SM, sqrt_s_LHeC_1_2); };
    obsThFactory["muepWBFtautau1200"] = [=](const StandardModel& SM) { return new muepWBFtautau(SM, sqrt_s_LHeC_1_2); };
    //
    obsThFactory["muepWBFbb1300"] = [=](const StandardModel& SM) { return new muepWBFbb(SM, sqrt_s_LHeC_1_3); };
    obsThFactory["muepWBFcc1300"] = [=](const StandardModel& SM) { return new muepWBFcc(SM, sqrt_s_LHeC_1_3); };
    obsThFactory["muepWBFgg1300"] = [=](const StandardModel& SM) { return new muepWBFgg(SM, sqrt_s_LHeC_1_3); };
    obsThFactory["muepWBFWW2l2v1300"] = [=](const StandardModel& SM) { return new muepWBFWW2l2v(SM, sqrt_s_LHeC_1_3); };
    obsThFactory["muepWBFZZ4l1300"] = [=](const StandardModel& SM) { return new muepWBFZZ4l(SM, sqrt_s_LHeC_1_3); };
    obsThFactory["muepWBFgaga1300"] = [=](const StandardModel& SM) { return new muepWBFgaga(SM, sqrt_s_LHeC_1_3); };
    obsThFactory["muepWBFtautau1300"] = [=](const StandardModel& SM) { return new muepWBFtautau(SM, sqrt_s_LHeC_1_3); };
    //
    obsThFactory["muepWBFbb1800"] = [=](const StandardModel& SM) { return new muepWBFbb(SM, sqrt_s_LHeC_1_8); };
    obsThFactory["muepWBFcc1800"] = [=](const StandardModel& SM) { return new muepWBFcc(SM, sqrt_s_LHeC_1_8); };
    obsThFactory["muepWBFgg1800"] = [=](const StandardModel& SM) { return new muepWBFgg(SM, sqrt_s_LHeC_1_8); };
    obsThFactory["muepWBFWW2l2v1800"] = [=](const StandardModel& SM) { return new muepWBFWW2l2v(SM, sqrt_s_LHeC_1_8); };
    obsThFactory["muepWBFZZ4l1800"] = [=](const StandardModel& SM) { return new muepWBFZZ4l(SM, sqrt_s_LHeC_1_8); };
    obsThFactory["muepWBFgaga1800"] = [=](const StandardModel& SM) { return new muepWBFgaga(SM, sqrt_s_LHeC_1_8); };
    obsThFactory["muepWBFtautau1800"] = [=](const StandardModel& SM) { return new muepWBFtautau(SM, sqrt_s_LHeC_1_8); };
    //
    obsThFactory["muepWBFbb3500"] = [=](const StandardModel& SM) { return new muepWBFbb(SM, sqrt_s_FCCep_3_5); };
    obsThFactory["muepWBFcc3500"] = [=](const StandardModel& SM) { return new muepWBFcc(SM, sqrt_s_FCCep_3_5); };
    obsThFactory["muepWBFgg3500"] = [=](const StandardModel& SM) { return new muepWBFgg(SM, sqrt_s_FCCep_3_5); };
    obsThFactory["muepWBFWW2l2v3500"] = [=](const StandardModel& SM) { return new muepWBFWW2l2v(SM, sqrt_s_FCCep_3_5); };
    obsThFactory["muepWBFZZ4l3500"] = [=](const StandardModel& SM) { return new muepWBFZZ4l(SM, sqrt_s_FCCep_3_5); };
    obsThFactory["muepWBFgaga3500"] = [=](const StandardModel& SM) { return new muepWBFgaga(SM, sqrt_s_FCCep_3_5); };
    obsThFactory["muepWBFtautau3500"] = [=](const StandardModel& SM) { return new muepWBFtautau(SM, sqrt_s_FCCep_3_5); };
    //
    obsThFactory["muepWBFbb5000"] = [=](const StandardModel& SM) { return new muepWBFbb(SM, sqrt_s_FCCep_5); };
    obsThFactory["muepWBFcc5000"] = [=](const StandardModel& SM) { return new muepWBFcc(SM, sqrt_s_FCCep_5); };
    obsThFactory["muepWBFgg5000"] = [=](const StandardModel& SM) { return new muepWBFgg(SM, sqrt_s_FCCep_5); };
    obsThFactory["muepWBFWW2l2v5000"] = [=](const StandardModel& SM) { return new muepWBFWW2l2v(SM, sqrt_s_FCCep_5); };
    obsThFactory["muepWBFZZ4l5000"] = [=](const StandardModel& SM) { return new muepWBFZZ4l(SM, sqrt_s_FCCep_5); };
    obsThFactory["muepWBFgaga5000"] = [=](const StandardModel& SM) { return new muepWBFgaga(SM, sqrt_s_FCCep_5); };
    obsThFactory["muepWBFtautau5000"] = [=](const StandardModel& SM) { return new muepWBFtautau(SM, sqrt_s_FCCep_5); };
    //
    obsThFactory["muepZBFbb1200"] = [=](const StandardModel& SM) { return new muepZBFbb(SM, sqrt_s_LHeC_1_2); };
    obsThFactory["muepZBFcc1200"] = [=](const StandardModel& SM) { return new muepZBFcc(SM, sqrt_s_LHeC_1_2); };
    obsThFactory["muepZBFgg1200"] = [=](const StandardModel& SM) { return new muepZBFgg(SM, sqrt_s_LHeC_1_2); };
    obsThFactory["muepZBFWW2l2v1200"] = [=](const StandardModel& SM) { return new muepZBFWW2l2v(SM, sqrt_s_LHeC_1_2); };
    obsThFactory["muepZBFZZ4l1200"] = [=](const StandardModel& SM) { return new muepZBFZZ4l(SM, sqrt_s_LHeC_1_2); };
    obsThFactory["muepZBFgaga1200"] = [=](const StandardModel& SM) { return new muepZBFgaga(SM, sqrt_s_LHeC_1_2); };
    obsThFactory["muepZBFtautau1200"] = [=](const StandardModel& SM) { return new muepZBFtautau(SM, sqrt_s_LHeC_1_2); };
    //
    obsThFactory["muepZBFbb1300"] = [=](const StandardModel& SM) { return new muepZBFbb(SM, sqrt_s_LHeC_1_3); };
    obsThFactory["muepZBFcc1300"] = [=](const StandardModel& SM) { return new muepZBFcc(SM, sqrt_s_LHeC_1_3); };
    obsThFactory["muepZBFgg1300"] = [=](const StandardModel& SM) { return new muepZBFgg(SM, sqrt_s_LHeC_1_3); };
    obsThFactory["muepZBFWW2l2v1300"] = [=](const StandardModel& SM) { return new muepZBFWW2l2v(SM, sqrt_s_LHeC_1_3); };
    obsThFactory["muepZBFZZ4l1300"] = [=](const StandardModel& SM) { return new muepZBFZZ4l(SM, sqrt_s_LHeC_1_3); };
    obsThFactory["muepZBFgaga1300"] = [=](const StandardModel& SM) { return new muepZBFgaga(SM, sqrt_s_LHeC_1_3); };
    obsThFactory["muepZBFtautau1300"] = [=](const StandardModel& SM) { return new muepZBFtautau(SM, sqrt_s_LHeC_1_3); };
    //
    obsThFactory["muepZBFbb1800"] = [=](const StandardModel& SM) { return new muepZBFbb(SM, sqrt_s_LHeC_1_8); };
    obsThFactory["muepZBFcc1800"] = [=](const StandardModel& SM) { return new muepZBFcc(SM, sqrt_s_LHeC_1_8); };
    obsThFactory["muepZBFgg1800"] = [=](const StandardModel& SM) { return new muepZBFgg(SM, sqrt_s_LHeC_1_8); };
    obsThFactory["muepZBFWW2l2v1800"] = [=](const StandardModel& SM) { return new muepZBFWW2l2v(SM, sqrt_s_LHeC_1_8); };
    obsThFactory["muepZBFZZ4l1800"] = [=](const StandardModel& SM) { return new muepZBFZZ4l(SM, sqrt_s_LHeC_1_8); };
    obsThFactory["muepZBFgaga1800"] = [=](const StandardModel& SM) { return new muepZBFgaga(SM, sqrt_s_LHeC_1_8); };
    obsThFactory["muepZBFtautau1800"] = [=](const StandardModel& SM) { return new muepZBFtautau(SM, sqrt_s_LHeC_1_8); };
    //
    obsThFactory["muepZBFbb3500"] = [=](const StandardModel& SM) { return new muepZBFbb(SM, sqrt_s_FCCep_3_5); };
    obsThFactory["muepZBFcc3500"] = [=](const StandardModel& SM) { return new muepZBFcc(SM, sqrt_s_FCCep_3_5); };
    obsThFactory["muepZBFgg3500"] = [=](const StandardModel& SM) { return new muepZBFgg(SM, sqrt_s_FCCep_3_5); };
    obsThFactory["muepZBFWW2l2v3500"] = [=](const StandardModel& SM) { return new muepZBFWW2l2v(SM, sqrt_s_FCCep_3_5); };
    obsThFactory["muepZBFZZ4l3500"] = [=](const StandardModel& SM) { return new muepZBFZZ4l(SM, sqrt_s_FCCep_3_5); };
    obsThFactory["muepZBFgaga3500"] = [=](const StandardModel& SM) { return new muepZBFgaga(SM, sqrt_s_FCCep_3_5); };
    obsThFactory["muepZBFtautau3500"] = [=](const StandardModel& SM) { return new muepZBFtautau(SM, sqrt_s_FCCep_3_5); };
    //
    obsThFactory["muepZBFbb5000"] = [=](const StandardModel& SM) { return new muepZBFbb(SM, sqrt_s_FCCep_5); };
    obsThFactory["muepZBFcc5000"] = [=](const StandardModel& SM) { return new muepZBFcc(SM, sqrt_s_FCCep_5); };
    obsThFactory["muepZBFgg5000"] = [=](const StandardModel& SM) { return new muepZBFgg(SM, sqrt_s_FCCep_5); };
    obsThFactory["muepZBFWW2l2v5000"] = [=](const StandardModel& SM) { return new muepZBFWW2l2v(SM, sqrt_s_FCCep_5); };
    obsThFactory["muepZBFZZ4l5000"] = [=](const StandardModel& SM) { return new muepZBFZZ4l(SM, sqrt_s_FCCep_5); };
    obsThFactory["muepZBFgaga5000"] = [=](const StandardModel& SM) { return new muepZBFgaga(SM, sqrt_s_FCCep_5); };
    obsThFactory["muepZBFtautau5000"] = [=](const StandardModel& SM) { return new muepZBFtautau(SM, sqrt_s_FCCep_5); };
    //
    //-----  Limits  ----------
    obsThFactory["UpperLimit_ppHZgammaA13"] = [=](const StandardModel& SM) { return new UpperLimit_ppHZgammaA13(SM, sqrt_s_LHC13); };
    obsThFactory["UpperLimit_ppHZgammaC13"] = [=](const StandardModel& SM) { return new UpperLimit_ppHZgammaC13(SM, sqrt_s_LHC13); };
    obsThFactory["UpperLimit_ppHZgammaA"] = [=](const StandardModel& SM) { return new UpperLimit_ppHZgammaA(SM, sqrt_s_LHC8); };
    obsThFactory["UpperLimit_ppHZgammaC"] = [=](const StandardModel& SM) { return new UpperLimit_ppHZgammaC(SM, sqrt_s_LHC8); };
    //-----  Others  ----------
    obsThFactory["cg_plus_ct"] = [](const StandardModel& SM) { return new cg_plus_ct(SM); };
    obsThFactory["cga_plus_ct"] = [](const StandardModel& SM) { return new cga_plus_ct(SM); };
    obsThFactory["cg_minus_cga"] = [](const StandardModel& SM) { return new cg_minus_cga(SM); };
    obsThFactory["cV_plus_cb"] = [](const StandardModel& SM) { return new cV_plus_cb(SM); };
    obsThFactory["cV_plus_ctau"] = [](const StandardModel& SM) { return new cV_plus_ctau(SM); };
    obsThFactory["cb_minus_cc"] = [](const StandardModel& SM) { return new cb_minus_cc(SM); };
    obsThFactory["cb_minus_ctau"] = [](const StandardModel& SM) { return new cb_minus_ctau(SM); };
    obsThFactory["cc_minus_ctau"] = [](const StandardModel& SM) { return new cc_minus_ctau(SM); };

    //-----  e+e- -> W+ W- Optimized Observables  -----
    obsThFactory["eeWWOO"] = [](const StandardModel& SM) { return new eeWW(SM); };

    //-----  Epsilon parameters  -----
    obsThFactory["epsilon1"] = [](const StandardModel& SM) { return new Epsilon1(SM); };
    obsThFactory["epsilon2"] = [](const StandardModel& SM) { return new Epsilon2(SM); };
    obsThFactory["epsilon3"] = [](const StandardModel& SM) { return new Epsilon3(SM); };
    obsThFactory["epsilonb"] = [](const StandardModel& SM) { return new Epsilonb(SM); };


    //----- NPSMEFT6dtopquark  -----


    //Now that the map is defined these lines should be useless. Remove them and check that is fine
    /*
    obsThFactory["C_phit"] = [](const StandardModel& SM) { return new C_phit(SM); };
    obsThFactory["C_phiQ3"] = [](const StandardModel& SM) { return new C_phiQ3(SM); };
    obsThFactory["C_phiQ1"] = [](const StandardModel& SM) { return new C_phiQ1(SM); };
    obsThFactory["C_phiQm"] = [](const StandardModel& SM) { return new C_phiQm(SM); };
    obsThFactory["C_tW"] = [](const StandardModel& SM) { return new C_tW(SM); };
    obsThFactory["C_tZ"] = [](const StandardModel& SM) { return new C_tZ(SM); };
    obsThFactory["C_tB"] = [](const StandardModel& SM) { return new C_tB(SM); };
    obsThFactory["C_tphi"] = [](const StandardModel& SM) { return new C_tphi(SM); };
    */
}
