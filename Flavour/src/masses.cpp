#include "masses.h"
#include "StandardModel.h"

up_mass::up_mass(const StandardModel& SM_i) : ThObservable(SM_i)
{
}

double up_mass::computeThValue()
{
    return(SM.getQuarks(QCD::UP).getMass());
}

down_mass::down_mass(const StandardModel& SM_i) : ThObservable(SM_i)
{
}

double down_mass::computeThValue()
{
    return(SM.getQuarks(QCD::DOWN).getMass());
}

strange_mass::strange_mass(const StandardModel& SM_i) : ThObservable(SM_i)
{
}

double strange_mass::computeThValue()
{
    return(SM.getQuarks(QCD::STRANGE).getMass());
}

charm_mass::charm_mass(const StandardModel& SM_i) : ThObservable(SM_i)
{
}

double charm_mass::computeThValue()
{
    return(SM.getQuarks(QCD::CHARM).getMass());
}

bottom_mass::bottom_mass(const StandardModel& SM_i) : ThObservable(SM_i)
{
}

double bottom_mass::computeThValue()
{
    return(SM.getQuarks(QCD::BOTTOM).getMass());
}

top_mass::top_mass(const StandardModel& SM_i) : ThObservable(SM_i)
{
}

double top_mass::computeThValue()
{
    return(SM.getQuarks(QCD::TOP).getMass());
}

mtpole::mtpole(const StandardModel& SM_i) : ThObservable(SM_i)
{
}

double mtpole::computeThValue()
{
    return(SM.getMtpole());
}

electron_mass::electron_mass(const StandardModel& SM_i) : ThObservable(SM_i)
{
}

double electron_mass::computeThValue()
{
    return(SM.getLeptons(QCD::ELECTRON).getMass());
}

muon_mass::muon_mass(const StandardModel& SM_i) : ThObservable(SM_i)
{
}

double muon_mass::computeThValue()
{
    return(SM.getLeptons(QCD::MU).getMass());
}

tau_mass::tau_mass(const StandardModel& SM_i) : ThObservable(SM_i)
{
}

double tau_mass::computeThValue()
{
    return(SM.getLeptons(QCD::TAU).getMass());
}
