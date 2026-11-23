// Run from the repository root with the HEPfit libraries linked.
#include "NPSTUVWXY.h"
#include <array>
#include <cmath>
#include <complex>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>

namespace {
using Complex = std::complex<double>;
using Parameters = std::array<double, 7>;
const char* names[] = {"obliqueShat", "obliqueThat", "obliqueUhat", "obliqueV",
                       "obliqueW", "obliqueX", "obliqueY"};

void initialize(NPSTUVWXY& model, bool massInput=false)
{
    if (!model.InitializeModel()) throw std::runtime_error("Model initialization failed");
    std::map<std::string, double> p = {
        {"lambda",0.2251}, {"A",0.826}, {"rhob",0.160}, {"etab",0.350}, {"Mz",91.1875},
        {"AlsMz",0.1185}, {"GF",1.1663787e-5}, {"ale",7.2973525698e-3}, {"dAle5Mz",0.0275},
        {"mHl",125.09}, {"mtop",163.0}, {"mbottom",4.18}, {"mcharm",1.3}, {"mstrange",0.0938},
        {"mup",0.0023}, {"mdown",0.0048}, {"muc",1.3}, {"mub",4.8}, {"mut",164.1}, {"muw",80.385},
        {"melectron",5.109989e-4}, {"mmu",0.10565837}, {"mtau",1.77682}
    };
    for (auto key : {"delMw","delSin2th_l","delSin2th_q","delSin2th_b","delGammaZ","delsigma0H",
                     "delR0l","delR0c","delR0b","delGammaWlv","delGammaWqq",
                     "mneutrino_1","mneutrino_2","mneutrino_3"}) p[key]=0.;
    if (massInput) {
        model.setFlag("MWinput", true);
        p.erase("dAle5Mz");
        p["Mw_inp"] = 80.37;
    }
    for (auto key : names) p[key] = 0.;
    if (!model.Init(p)) throw std::runtime_error("Missing model inputs");
    model.setFlag("LEP2WeakBox", false);
}

struct Integral {
    double total, fb;
};

struct Reference {
    double s, mw, mz, width, e, gz, sw, cw;
    bool massInput;
    Reference(const StandardModel& sm, double energy) : s(energy*energy), mw(sm.Mw()), mz(sm.getMz()),
        width(sm.Gamma_Z()), e(std::sqrt(4.*M_PI*sm.ale_OS(energy, FULLNLO))),
        gz(std::sqrt(4.*std::sqrt(2.)*sm.getGF())*mz),
        sw(std::sqrt(sm.sW2())), cw(std::sqrt(sm.cW2())), massInput(sm.IsFlagMWinput()) {}

    // Neutral-current chiral amplitudes. Rotate the contact interaction from
    // (W3, B) before squaring. Eq. (10) of hep-ph/0405040 fixes the pole terms.
    Integral evaluate(const Particle& f, double m, double colors, const Parameters& p) const
    {
        const double beta = std::sqrt(1.-4.*m*m/s);
        const Complex pole = s / Complex(s-mz*mz, s*width/mz);
        const double residue = p[1]-p[4]+2.*sw/cw*p[5]-sw*sw/(cw*cw)*p[6];
        const double massInputMixing = cw/sw*(p[1]-p[2]+p[3])-sw/cw*(p[0]-p[4])-p[5]/(cw*cw);
        // Transform Eq. (10) from the mass-defined angle using the oblique
        // W-mass shift at fixed alpha, G_F, M_Z (S', T', U' convention).
        const double e2 = p[2]-p[3]-p[4]+2.*sw/cw*p[5];
        const double e3 = p[0]-p[4]+p[5]/(sw*cw)-p[6];
        const double dMw2OverMw2 = (cw*cw*residue-2.*sw*sw*e3)/(cw*cw-sw*sw)-e2;
        const double mixing = gz/e*(sw*cw*massInputMixing - (massInput ? 0. : cw*cw*dMw2OverMw2));
        const double photonShift = massInput ? (cw*cw-sw*sw)/(sw*sw)*dMw2OverMw2 : 0.;
        Integral result{};
        for (int incoming : {1,-1}) {
            const double ze = gz*((incoming == 1 ? -0.5 : 0.) + sw*sw);
            std::array<Complex,2> amp;
            for (int k=0; k<2; ++k) {
                const double qf = e*f.getCharge();
                const double zf = gz*((k == 0 ? f.getIsospin() : 0.) - sw*sw*f.getCharge());
                const double w3e = -sw*e + cw*ze, be = -cw*e - sw*ze;
                const double w3f = sw*qf + cw*zf, bf = cw*qf - sw*zf;
                amp[k] = -e*qf*(1.+photonShift) + pole*(ze*zf*(1.+residue) - mixing*(ze*qf-e*zf))
                    - s/(mw*mw)*(p[4]*w3e*w3f + p[6]*be*bf + p[5]*(w3e*bf+be*w3f));
            }
            const Complex v = (amp[0]+amp[1])/2., a = (amp[0]-amp[1])/2.;
            for (int half : {-1,1}) {
                for (int node : {-1,1}) {
                    // Two-point Gauss integration is exact on each hemisphere.
                    const double x = half*(1.+node/std::sqrt(3.))/2.;
                    const double weight = std::pow(1.+incoming*x,2)*std::norm(v+beta*a)
                            + std::pow(1.-incoming*x,2)*std::norm(v-beta*a)
                            + 2.*(1.-beta*beta)*(1.-x*x)*std::norm(v);
                    const double ds = 0.5*colors*beta*weight/(128.*M_PI*s);
                    result.total += ds;
                    result.fb += half*ds;
                }
            }
        }
        return result;
    }
};

void check(double value, double expected, double tolerance, const std::string& which)
{
    if (!(std::abs(value-expected) <= tolerance)) {
        std::ostringstream out;
        out.precision(15);
        out << which << " got " << value << ", expected " << expected;
        throw std::runtime_error(out.str());
    }
}

void testDirections(const StandardModel& sm, double energy)
{
    const Reference ref(sm, energy);
    const LEP2oblique lep(sm);
    const double s = ref.s, step = 1.e-6;
    for (unsigned j=0; j<7; ++j) {
        Parameters plus{}, minus{}, unit{};
        plus[j]=step; minus[j]=-step; unit[j]=1.;
        double hadPlus=0., hadMinus=0., bPlus=0., bMinus=0.;
        auto checkRate = [&](const Particle& f, double mass, int colors, double ds, double da) {
            const Integral hi = ref.evaluate(f,mass,colors,plus);
            const Integral lo = ref.evaluate(f,mass,colors,minus);
            const double expected = (hi.total-lo.total)/(2.*step);
            check(ds, expected, 2.e-9*(hi.total+std::abs(expected)), "cross section");
            check(da,(hi.fb/hi.total-lo.fb/lo.total)/(2.*step),2.e-8,"asymmetry");
        };
        for (auto l : {QCD::MU,QCD::TAU}) {
            const Particle& f=sm.getLeptons(l);
            checkRate(f,f.getMass(),1,lep.sigma_l_LEP2_NP(l,s,f.getMass(),unit.data()),
                    lep.AFB_l_LEP2_NP(l,s,f.getMass(),unit.data()));
        }
        for (auto q : {QCD::UP,QCD::DOWN,QCD::STRANGE,QCD::CHARM,QCD::BOTTOM}) {
            const Particle& f=sm.getQuarks(q);
            const double mass=sm.getmq(q,energy);
            checkRate(f,mass,3,lep.sigma_q_LEP2_NP(q,s,mass,unit.data()),
                    lep.AFB_q_LEP2_NP(q,s,mass,unit.data()));
            const double hi=ref.evaluate(f,mass,3,plus).total, lo=ref.evaluate(f,mass,3,minus).total;
            hadPlus+=hi; hadMinus+=lo;
            if (q==QCD::BOTTOM) { bPlus=hi; bMinus=lo; }
        }
        check(lep.R_q_LEP2_NP(QCD::BOTTOM,s,sm.getmq(QCD::BOTTOM,energy),unit.data()),
              (bPlus/hadPlus-bMinus/hadMinus)/(2.*step),2.e-8,"Rb");
    }
}
}

int main()
{
    NPSTUVWXY model;
    initialize(model);
    const StandardModel& sm=model.getTrueSM();
    for (double energy : {130.2,188.6,206.7}) testDirections(sm,energy);
    NPSTUVWXY massInput;
    initialize(massInput, true);
    for (double energy : {130.2,188.6,206.7}) testDirections(massInput.getTrueSM(),energy);
    for (NPSTUVWXY* fit : {&model,&massInput}) {
        const double s=199.5*199.5;
        if (!fit->Update({{"obliqueUhat",0.001}})) throw std::runtime_error("U update failed");
        const double shift=fit->LEP2sigmaMu(s)-fit->getTrueSM().LEP2sigmaMu(s);
        if (fit->getTrueSM().IsFlagMWinput()) {
            if (!(std::abs(shift)>1.e-6)) throw std::runtime_error("Missing M_W input conversion");
        } else {
            check(shift,0.,0.,"neutral current at fixed alpha");
        }
        if (!fit->Update({{"obliqueUhat",0.}})) throw std::runtime_error("U reset failed");
    }
    for (bool exact : {false,true}) {
        model.setFlag("NoApproximateGammaZ",exact);
        const double s=199.5*199.5;
        if (!model.Update({{"obliqueX",0.}})) throw std::runtime_error("Zero-point update failed");
        const double baseline=sm.LEP2sigmaMu(s);
        check(model.LEP2sigmaMu(s),baseline,0.,"zero shift");
        const double baselineAFB=sm.LEP2AFBmu(s);
        if (!model.Update({{"obliqueX",0.001}})) throw std::runtime_error("X update failed");
        const Parameters x={0.,0.,0.,0.,0.,0.001,0.};
        LEP2oblique lep(sm);
        const double units=exact ? 1. : sm.GeVminus2_to_nb*1000.;
        const double delta=lep.sigma_l_LEP2_NP(QCD::MU,s,sm.getLeptons(QCD::MU).getMass(),x.data());
        check(model.LEP2sigmaMu(s)-baseline,units*delta,1.e-12*std::abs(baseline),"NPSTUVWXY units");
        check(model.LEP2AFBmu(s)-baselineAFB,lep.AFB_l_LEP2_NP(QCD::MU,s,sm.getLeptons(QCD::MU).getMass(),x.data()),
              1.e-12,"NPSTUVWXY asymmetry");
        std::cout << (exact ? "exact" : "tabulated") << " 199.5 GeV: delta sigma(mu) = "
                  << delta*sm.GeVminus2_to_nb*1000. << " pb\n";
    }
}
