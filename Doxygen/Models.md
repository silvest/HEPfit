Models   {#PageModels}
==========================================================

Here is a list of all available models with links to the summary
tables of model parameters and flags. The parameters and flags 
are inherited from a parent model class to a child model class: 

MODEL_GRAPH_INHERITE_SVG
  
where Model is a template class of models, and NPbase is an auxiliary class.
A complete list of observables are available in the function ThObsFactory::ThObsFactory(),
where not all the observables can be used in each model. For some of the models below we list
the available observables since the list is small.

## QCD:

  - %Model parameters: [@ref QCDParameters "Summary table"]
  - %Model flags: None
  - __NOTE:__  %QCD cannot be set as a model in the model configuration file

## StandardModel:

  - %Model parameters: [@ref StandardModelParameters "Summary table"]
  - %Model flags: [@ref StandardModelFlags "Summary table"]
  - %Model description: This is a Model class containing parameters and functions associated
   with the Standard %Model. This class is inherited from the QCD class, which
   defines parameters related to %QCD
  - __NOTE:__  %StandardModel is the minimum model that needs to be set in the model configuration file
  
##  NPSTU
  
  - %Model parameters: [@ref NPSTUParameters "Summary table
  - %Model flags: None 
  - %Model description: This is a %Model class containing the necessary functions to compute
  new physics contributions to the electroweak precision observables with the
  Peskin-Takeuchi oblique parameters @cite Peskin:1990zt, @cite Peskin:1991sw.
  
## NPSTUZbbbarLR
  
  - %Model parameters: [@ref NPSTUZbbbarLRParameters "Summary table"]
  - %Model flags: None
  - %Model description: This is a %Model class containing the necessary functions to compute
   new physics contributions to the electroweak precision observables with the
   Peskin-Takeuchi oblique parameters @cite Peskin:1990zt, @cite Peskin:1991sw
   and with additional shifts to the left-handed and right-handed
   @f$Zb\bar{b}@f$ vertices.

## NPEpsilons:

  - %Model parameters: [@ref NPEpsilonsParameters "Summary table"]
  - %Model flags: [@ref NPEpsilonsFlags "Summary table"]
  - %Model description: This is a %Model class containing parameters and functions to work
    with the epsilon paramterization 
    @cite Altarelli:1990zd, @cite Altarelli:1991fk,@cite Altarelli:1993sz
    for new physics contributions to the electroweak precision observables.

## NPEpsilons_pureNP:

  - %Model parameters: [@ref NPEpsilons_pureNPParameters "Summary table"]
  - %Model flags: None
  - %Model description: This is a %Model class containing parameters and functions to work 
   with the modified epsilon parameters \f$\delta\,\varepsilon_{1,2,3,b}\f$
   which parameterize only new physics contributions to the electroweak
   precision observables. 
   See @cite Altarelli:1990zd, @cite Altarelli:1991fk,@cite Altarelli:1993sz
   for the original epsilon parameterization.

## NPZbbbar (variant: NPZbbbarLR):

  - %Model parameters: [@ref NPZbbbarParameters "Summary table"]
  - %Model flags: [@ref NPZbbbarFlags "Summary table"]
  - %Model description: This class is a %Model class containing the necessary functions to
   work with new physics contributions to electroweak precision observables,
   in the form of contributions to the neutral current couplings of the bottom
   quark. Variant: \f$\delta g_{L}^b\f$
   and  \f$\delta g_{R}^b\f$ are used for the model parameters, instead of
   \f$\delta g_{V}^b\f$ and \f$\delta g_{A}^b\f$

## NPZbbbarLinearized (variant: NPZbbbarLinearizedLR):

  - %Model parameters: [@ref NPZbbbarLinearizedParameters "Summary table"]
  - %Model flags: [@ref NPZbbbarLinearizedFlags "Summary table"]
  - %Model description: This class is a linearized %Model class containing the necessary functions to
   work with new physics contributions to electroweak precision observables,
   in the form of contributions to the neutral current couplings of the bottom
   quark. Variant: \f$\delta g_{L}^b\f$
   and  \f$\delta g_{R}^b\f$ are used for the model parameters, instead of
   \f$\delta g_{V}^b\f$ and \f$\delta g_{A}^b\f$
	
## NPSMEFTd6 (variant: NPSMEFTd6_LFU_QFU)

  - %Model parameters: [@ref NPSMEFTd6Parameters "Summary table"]
  - %Model flags: [@ref NPSMEFTd6Flags "Summary table"]
  - %Model description: This is a %Model class containing parameters and functions
   associated with the general dimension-six effective Lagrangian. 
   Use the model name "NPSMEFTd6_LFU_QFU" to assume lepton and quark flavour universality. 
   The implementation is written in the basis of \cite Grzadkowski:2010es. 
   For convenience, the parameterization also includes operators appearing in
   other common bases. In particular, the complete set of parameters contains 4
   redundancies, given by the coefficients \f$C_{2B,2W,DHB,DHW,DB,DW} \f$,
   which correspond to operators not included in the basis of \cite Grzadkowski:2010es.
   For meaningful physical results one must make sure to include only
   a complete set of interactions in a given analysis.

## HiggsKvKf

  - %Model parameters: [@ref HiggsKvKfParameters "Summary table"]
  - %Model flags: None
  - %Model description: This is a %Model class containing parameters and functions associated
   with an extension of the %StandardModel where Higgs couplings to all vector bosons
   are rescaled by @f$\kappa_v@f$ and Higgs couplings to all fermions are rescaled by @f$\kappa_f@f$.
   The invisible decay width is also parametrized independently by Br@f$(H\to invisible)@f$.
   This class inherits from the %NPbase class, which defines parameters related to generic
   extensions of the %StandardModel Higgs sector.

## HiggsKvKfgen

  - %Model parameters: [@ref HiggsKvKfgenParameters "Summary table"]
  - %Model flags: None
  - %Model description: This is a %Model class containing parameters and functions associated
   with an extension of the %StandardModel where Higgs couplings to all vector bosons
   are rescaled by @f$\kappa_v@f$ and Higgs couplings to all up, down and lepton fermions
   are rescaled by @f$\kappa_u, \kappa_d, \kappa_l@f$ respectively.
   The invisible decay width is also parametrized independently by Br@f$(H\to invisible)@f$.
   This class inherits from the %NPbase class, which defines parameters related to generic
   extensions of the %StandardModel Higgs sector.

## HiggsKvgenKfgen

  - %Model parameters: [@ref HiggsKvgenKfgenParameters "Summary table"]
  - %Model flags: None
  - %Model description: This is a %Model class containing parameters and functions associated
   with an extension of the %StandardModel where Higgs couplings to W bosons
   are rescaled by @f$\kappa_W@f$, Higgs couplings to Z bosons
   are rescaled by @f$\kappa_Z@f$ and Higgs couplings to all up, down and lepton fermions
   are rescaled by @f$\kappa_u, \kappa_d, \kappa_l@f$ respectively.
   The invisible decay width is also parametrized independently by Br@f$(H\to invisible)@f$.
   This class inherits from the %NPbase class, which defines parameters related to generic
   extensions of the %StandardModel Higgs sector.
  
## HiggsKigen

  - %Model parameters: [@ref HiggsKigenParameters "Summary table"]
  - %Model flags: [@ref HiggsKigenFlags "Summary table"]
  - %Model description: This is a %Model class containing parameters and functions associated
   to rescaling the Higgs decay into vector bosons (@f$\kappa_W,\kappa_Z@f$), 
   gluons (@f$\kappa_g@f$), photons (@f$\kappa_\gamma@f$), Z and photons (@f$\kappa_{Z\gamma}@f$) 
   and fermions (@f$\kappa_{u,c,t}, \kappa_{d,s,b}, \kappa_{e\mu\tau}@f$)
   (as well as the corresponding production mechanisms) with respect to the %StandardModel.
   The possibility of extra decay width is also parametrized independently by Br@f$(H\to invisible)@f$ and Br@f$(H\to exotic)@f$.
   This class inherits from the %NPbase class, which defines parameters related to generic
   extensions of the %StandardModel Higgs sector.
  
## HiggsChiral

  - %Model parameters: [@ref HiggsChiralParameters "Summary table"]
  - %Model flags: [@ref HiggsChiralFlags "Summary table"]
  
## FlavourWilsonCoefficient

  - %Model parameters: [@ref FlavourWilsonCoefficientParameters "Summary table"]
  - %Model flags: None
  - %Model description: A %Model for NP contributions to flavour
    through shifts to Standard %Model Wilson coefficients.
  
## RealWeakEFTLFV

  - %Model parameters: [@ref RealWeakEFTLFVParameters "Summary table"]
  - %Model flags: None
  - %Model description: A %Model for WEFT LFV contributions to @f$\Delta F=1@f$
   processes like @f$ b\to s@f$ decays.
  
## RealWeakEFTCC (variant: RealWeakEFTCCPM)

  - %Model parameters: [@ref RealWeakEFTCCParameters "Summary table"]
  - %Model flags: None
  - %Model description: A %Model for NP contributions to charged current
   processes like @f$b\to c@f$ decays. Variant: the +/- basis for the
   Wilson coefficients.
  
## THDM  

  - %Model parameters: [@ref THDMParameters "Summary table"]
  - %Model flags: [@ref THDMFlags "Summary table"]
  - %Model description: The @f$Z_2@f$ symmetric Two-Higgs-Doublet models.
   The theoretical constraints are the positivity bounds (boundedness from below)
   and the requirement that the electroweak minimum be the global minimum of the
   scalar potential as well as the LO and NLO(+) unitarity conditions.
   The experimental constraints comprise electroweak precision observables (%STU),
   Higgs signal strengths, direct searches for scalar particles and a set of flavour observables
   (@f$B \to X_s \gamma@f$, @f$B_s@f$ mixing, @f$B \to \tau \nu@f$, @f$B \to D^{(*)} \tau \nu@f$
   and the anomalous magnetic moment of the muon).
   Also the renormalization group equations are implemented at next-to-leading order.
  
## GeorgiMachacek  

  - %Model parameters: [@ref GeorgiMachacekParameters "Summary table"]
  - %Model flags: [@ref GeorgiMachacekFlags "Summary table"]
  - %Model description: The Georgi-Machacek model extends the Standard Model by two scalar triplets.
   Among the implemented theoretical constraints are positivity of the scalar
   potential and the unitarity conditions.
   As experimental constraints, Higgs signal strengths and direct searches for 
   neutral, singly and doubly charged scalars
   are available.
  
## GeneralTHDM

  - %Model parameters: [@ref GeneralTHDMParameters "Summary table"]
  - %Model flags: [@ref GeneralTHDMFlags "Summary table"]
  - %Model description: The general symmetric Two-Higgs-Doublet models.
   While the theoretical bounds (positivity, unitarity and stability of the Higgs potential)
   as well as the electroweak %STU pseudo-observables
   are available in the most general case,
   the Higgs and flavour observables are currently only implemented in the flavour aligned limiting case.
  
## THDMW    

  - %Model parameters: [@ref THDMWParameters "Summary table"]
  - %Model flags: [@ref THDMWFlags "Summary table"]
  - %Model description: The extension of either the Standard %Model or the @f$Z_2@f$ symmetric Two-Higgs-Doublet model by a scalar octet.
   The Standard Model extension by a scalar octet can be accessed setting the THDMWmodel flag to "ManoharWise".
   A custodial limiting case of this model can be obtained setting the flag to "custodialMW".
   The THDM plus octet can be address by attributing "custodial1" to the THDMWmodel flag.
   The implemented observables are positivity, unitarity, Higgs signal strengths and electroweak precision measurements.
