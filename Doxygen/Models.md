Models   {#PageModels}
==========================================================

Here is a list of all available models with links to the summary
tables of model parameters and flags. The parameters and flags 
are inherited from a parent model class to a child model class: 

MODEL_GRAPH_INHERITE_SVG

where Model is a template class of models, and NPbase is an auxiliary class.
A complete list of observables are available in the function ThObsFactory::ThObsFactory(),
where not all the observables can be used in each model. Below we also list the
available observables for each model.

## QCD:

  - %Model parameters: [@ref QCDParameters "Summary table"]
  - %Model flags: None
  - %Observables: None

## StandardModel:

  - %Model parameters: [@ref StandardModelParameters "Summary table"]
  - %Model flags: [@ref StandardModelFlags "Summary table"]
  - %Observables: Mw, GammaW, GammaZ, sigmaHadron, sin2thetaEff,
    PtauPol, Alepton, Acharm, Abottom, AFBlepton, AFBcharm, AFBbottom,
    Rlepton, Rcharm, Rbottom

## NPEpsilons:

  - %Model parameters: [@ref NPEpsilonsParameters "Summary table"]
  - %Model flags: [@ref NPEpsilonsFlags "Summary table"]
  - %Observables: Mw, GammaZ, sigmaHadron, sin2thetaEff,
    PtauPol, Alepton, Acharm, Abottom, AFBlepton, AFBcharm, AFBbottom,
    Rlepton, Rcharm, Rbottom

## NPEpsilons_pureNP:

  - %Model parameters: [@ref NPEpsilons_pureNPParameters "Summary table"]
  - %Model flags: None
  - %Observables: Mw, GammaZ, sigmaHadron, sin2thetaEff,
    PtauPol, Alepton, Acharm, Abottom, AFBlepton, AFBcharm, AFBbottom,
    Rlepton, Rcharm, Rbottom

## NPSTU:

  - %Model parameters: [@ref NPSTUParameters "Summary table"]
  - %Model flags: None
  - %Observables: Mw, GammaW, GammaZ, sigmaHadron, sin2thetaEff,
    PtauPol, Alepton, Acharm, Abottom, AFBlepton, AFBcharm, AFBbottom,
    Rlepton, Rcharm, Rbottom

## NPHiggs:

  - %Model parameters: [@ref NPHiggsParameters "Summary table"]
  - %Model flags: None
  - %Observables: Mw, GammaW, GammaZ, sigmaHadron, sin2thetaEff,
    PtauPol, Alepton, Acharm, Abottom, AFBlepton, AFBcharm, AFBbottom,
    Rlepton, Rcharm, Rbottom

## NPZbbbar:

  - %Model parameters: [@ref NPZbbbarParameters "Summary table"]
  - %Model flags: [@ref NPZbbbarFlags "Summary table"]
  - %Observables: Mw, GammaW, GammaZ, sigmaHadron, sin2thetaEff,
    PtauPol, Alepton, Acharm, Abottom, AFBlepton, AFBcharm, AFBbottom,
    Rlepton, Rcharm, Rbottom

## NPZbbbarLR (a variant of NPZbbbar): 

  - %Model parameters: [@ref NPZbbbarParameters "Summary table"],
    where deltaGLb and deltaGRb must be used, instead of deltaGVb and
    deltaGAb in %NPZbbbar above.
  - %Model flags: Same as %NPZbbbar above. 
  - %Observables: Same as %NPZbbbar above. 

## NPZbbbarLinearized:

  - %Model parameters: [@ref NPZbbbarLinearizedParameters "Summary table"]
  - %Model flags: [@ref NPZbbbarLinearizedFlags "Summary table"]
  - %Observables: Mw, GammaW, GammaZ, sigmaHadron, sin2thetaEff,
    PtauPol, Alepton, Acharm, Abottom, AFBlepton, AFBcharm, AFBbottom,
    Rlepton, Rcharm, Rbottom

## NPZbbbarLinearizedLR (a variant of NPZbbbarLinearized):

  - %Model parameters: [@ref NPZbbbarLinearizedParameters "Summary table"],
    where deltaGLb and deltaGRb must be used, instead of deltaGVb and
    deltaGAb in %NPZbbbarLinearized above.
  - %Model flags: Same as %NPZbbbarLinearized above.
  - %Observables: Same as %NPZbbbarLinearized above.

## NPEffectiveBS (variants: NPEffectiveBS_LFU, NPEffectiveBS_QFU and NPEffectiveBS_LFU_QFU):

  - %Model parameters: [@ref NPEffectiveBSParameters "Summary table"]
  - %Model flags: None
  - %Observables: Mw, GammaW, GammaZ, sigmaHadron, sin2thetaEff,
    PtauPol, Alepton, Acharm, Abottom, AFBlepton, AFBcharm, AFBbottom,
    Rlepton, Rcharm, Rbottom

## NPEffectiveGIMR (variant: NPEffectiveGIMR_LFU_QFU):

  - %Model parameters: [@ref NPEffectiveGIMRParameters "Summary table"]
  - %Model flags: None
  - %Observables: Mw, GammaW, GammaZ, sigmaHadron, sin2thetaEff,
    PtauPol, Alepton, Acharm, Abottom, AFBlepton, AFBcharm, AFBbottom,
    Rlepton, Rcharm, Rbottom

## HiggsKvKf

  - %Model parameters: [@ref HiggsKvKfParameters "Summary table"]
  - %Model flags:
  - %Observables:

## HiggsKvKfgen

  - %Model parameters: [@ref HiggsKvKfgenParameters "Summary table"]
  - %Model flags:
  - %Observables:

## HiggsKvgenKf

  - %Model parameters: [@ref HiggsKvgenKfParameters "Summary table"]
  - %Model flags:
  - %Observables:

