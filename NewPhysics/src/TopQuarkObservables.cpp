#include "TopQuarkObservables.h"

TopQuarkObservables::TopQuarkObservables(const NPSMEFTd6General& NP_i) : NP(NP_i) {};





/////// Helicities  /////////////////////////////////////////////////////////////////////////////////////


F0_LO::F0_LO(const StandardModel& SM_i)
: ThObservable(SM_i)
{};

double F0_LO::computeThValue()
{
    TopQuarkObservables myt(static_cast<const NPSMEFTd6General&>(SM));
    
    double F0_madgraph = 0.70381;
    
    bool   flag_Quadratic= false; //Needs to be properly defined
    double F0_SM = F0_madgraph; //Needs to be properly defined
    double C_tW = myt.ewgc("CuWR",2,2);
    double C_bW = myt.ewgc("CdWR",2,2);
    
    
    if(flag_Quadratic){
        return  F0_SM + (-0.0569231*C_tW + 0.0022956*C_tW*C_tW -
                0.00215306*(C_bW/(0.999*0.6532))*(C_bW/(0.999*0.6532))
                )*(F0_SM/F0_madgraph);
        }
        else{
            return F0_SM + (-0.0569231*(C_tW) )*(F0_SM/F0_madgraph);
        }

}


FL_LO::FL_LO(const StandardModel& SM_i)
: ThObservable(SM_i)
{};

double FL_LO::computeThValue()
{
    
    TopQuarkObservables myt(static_cast<const NPSMEFTd6General&>(SM));
    
    double FL_madgraph = 0.29619;

    bool   flag_Quadratic= false;//Needs to be properly defined
    //double FL_SM = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_FL_SM();
    double FL_SM = FL_madgraph;//Needs to be properly defined
    
    double C_tW = myt.ewgc("CuWR",2,2);
    double C_bW = myt.ewgc("CdWR",2,2);
    double C_phitb = myt.ewgc("CHudR",2,2);
    

    if(flag_Quadratic){
        return  FL_SM + ( 0.0569227*C_tW - 0.0022956*C_tW*C_tW 
                - 0.0011115*((C_bW/(0.999*0.6532)))*((C_bW/(0.999*0.6532)))
                - 0.000105115*((-C_phitb/(0.998)))*((-C_phitb/(0.998)))
                )*(FL_SM/FL_madgraph);
        }
    else{
        return FL_SM + (0.0569227*(C_tW) )*(FL_SM/FL_madgraph);
        }

    
}







Test_direct::Test_direct(const StandardModel& SM_i) 
: ThObservable(SM_i), mytopobs(static_cast<const NPSMEFTd6General&> (SM))
{ };

double Test_direct::computeThValue()
{
    
    //std::cout<<"\033[1;31m   myt.GetNP().getAle()  =  \033[0m "<< myt.GetNP().getAle()  <<std::endl;

    return ewgc("CHl3R",0,0);
    //return mytopobs.GetNP().test_direct();
    

    
}



