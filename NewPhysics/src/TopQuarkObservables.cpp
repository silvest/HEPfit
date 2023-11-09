#include "TopQuarkObservables.h"
#include "std_make_vector.h"


TopQuarkObservables::TopQuarkObservables(const NPSMEFTd6General& NP_i) : NP(NP_i) {};





/////// Helicities  /////////////////////////////////////////////////////////////////////////////////////


F0_LO::F0_LO(const StandardModel& SM_i)
: ThObservable(SM_i), mytopobs(static_cast<const NPSMEFTd6General&> (SM))
{
    setParametersForObservable(make_vector<std::string>() << "F0_SM");
};

double F0_LO::computeThValue()
{
    
    bool   flag_Quadratic= false; //Needs to be properly defined
    
    double F0_madgraph = 0.70381;
    double F0_SM = SM.getOptionalParameter("F0_SM"); 

    
    
    if(flag_Quadratic){
        return  F0_SM + (-0.0569231*getSMEFTCoeffEW("CuWR",2,2) + 0.0022956*getSMEFTCoeffEW("CuWR",2,2)*getSMEFTCoeffEW("CuWR",2,2) -
                0.00215306*(getSMEFTCoeffEW("CdWR",2,2)/(0.999*0.6532))*(getSMEFTCoeffEW("CdWR",2,2)/(0.999*0.6532))
                )*(F0_SM/F0_madgraph);
        }
        else{
            return F0_SM + (-0.0569231*getSMEFTCoeffEW("CuWR",2,2) )*(F0_SM/F0_madgraph);
        }

}



FL_LO::FL_LO(const StandardModel& SM_i)
: ThObservable(SM_i), mytopobs(static_cast<const NPSMEFTd6General&> (SM))
{
    setParametersForObservable(make_vector<std::string>() << "FL_SM");
};

double FL_LO::computeThValue()
{
    
    bool   flag_Quadratic= false;//Needs to be properly defined
    
    double FL_madgraph = 0.29619;
    double FL_SM = SM.getOptionalParameter("FL_SM");
    

    if(flag_Quadratic){
        return  FL_SM + ( 0.0569227*getSMEFTCoeffEW("CuWR",2,2) - 0.0022956*getSMEFTCoeffEW("CuWR",2,2)*getSMEFTCoeffEW("CuWR",2,2) 
                - 0.0011115*((getSMEFTCoeffEW("CdWR",2,2)/(0.999*0.6532)))*((getSMEFTCoeffEW("CdWR",2,2)/(0.999*0.6532)))
                - 0.000105115*((-getSMEFTCoeffEW("CHudR",2,2)/(0.998)))*((-getSMEFTCoeffEW("CHudR",2,2)/(0.998)))
                )*(FL_SM/FL_madgraph);
        }
    else{
        return FL_SM + (0.0569227*(getSMEFTCoeffEW("CuWR",2,2)) )*(FL_SM/FL_madgraph);
        }

    
}




Test_direct::Test_direct(const StandardModel& SM_i) 
: ThObservable(SM_i), mytopobs(static_cast<const NPSMEFTd6General&> (SM))
{ 

    setParametersForObservable(make_vector<std::string>() << "Test_direct_in");
};

double Test_direct::computeThValue()
{
    double Test_direct_in = SM.getOptionalParameter("Test_direct_in");
    
    return Test_direct_in;

    
}



