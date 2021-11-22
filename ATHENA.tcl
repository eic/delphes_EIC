######################################################################################################################
# ATHENA detector model. Based on parametrizations from G4 simulations made by the ATHENA collaboration
# email: miguel.arratia@ucr.edu, ssekula@mail.smu.edu
#######################################################################################################################

#######################################
# Order of execution of various modules
#######################################

set ExecutionPath {
  BeamSpotFilter
  ParticlePropagator

  ChargedHadronTrackingEfficiency
  ElectronTrackingEfficiency
  MuonTrackingEfficiency

  ChargedHadronSmearing
  ElectronSmearing
  MuonSmearing

  TrackMerger

  ECal
  HCal

  Calorimeter
  EFlowMerger
  EFlowFilter

  PhotonEfficiency
  PhotonIsolation

  ElectronFilter
  ElectronEfficiency
  ElectronIsolation

  ChargedHadronFilter
  MissingET

  NeutrinoFilter
  GenJetFinder
  GenMissingET

  FastJetFinder

  JetEnergyScale

  JetFlavorAssociation
  GenJetFlavorAssociation

  UniqueObjectFinder

  ScalarHT

  TrackCountingBTagging

  pfRICH
  barrelDIRC
  dualRICH_aerogel
  dualRICH_c2f6

  TreeWriter
}


#################################
# Propagate particles in cylinder
#################################

#######################
# GenBeamSpotFilter
# Saves a particle intended to represent the beamspot
#######################

module BeamSpotFilter BeamSpotFilter {
    set InputArray Delphes/stableParticles
    set OutputArray beamSpotParticle

}



module ParticlePropagator ParticlePropagator {
    set InputArray Delphes/stableParticles
    set OutputArray stableParticles
    set BeamSpotInputArray BeamSpotFilter/beamSpotParticle
    set ChargedHadronOutputArray chargedHadrons
    set ElectronOutputArray electrons
    set MuonOutputArray muons
    # radius of the magnetic field coverage, in m
    set Radius 1.60
    # half-length of the magnetic field coverage, in m
    set HalfLength 1.92
    # magnetic field
    set Bz 3.0
}


####################################
# Common Tracking Efficiency Model
####################################
#
set CommonTrackingEfficiency {
    (eta > -3.5 && eta <= -3.0 ) * (pt*cosh(eta) > 1.25 && pt*cosh(eta)<6.0)   * (0.875) +
    (eta > -3.5 && eta <= -3.0 ) * (pt*cosh(eta) > 6.0 )                       * (0.95) +
    (eta > -3.0 && eta <= -2.5 ) * (pt*cosh(eta) > 0.55 && pt*cosh(eta)<2.0)   * (0.875) +
    (eta > -3.0 && eta <= -2.5 ) * (pt*cosh(eta) > 2.0)                        * (0.95) +
    (eta > -2.5 && eta <= -2.0)  * (pt*cosh(eta)> 0.45 && pt*cosh(eta)<0.6)    * (0.875) +
    (eta > -2.5 && eta <= -2.0)  * (pt*cosh(eta)>0.6)                          * (0.95) +
    (eta > -2.0 && eta <= -1.5)  * (pt*cosh(eta)> 0.250 && pt*cosh(eta)<0.500)* (0.875) +
    (eta > -2.0 && eta <= -1.5)  * (pt*cosh(eta)>0.500)                        * (0.95) +    
    (eta > -1.5 && eta <= -1.0)  * (pt*cosh(eta) > 0.150 && pt*cosh(eta)<0.300)* (0.86) +
    (eta > -1.5 && eta <= -1.0)  * (pt*cosh(eta) > 0.300)                      * (0.92) +
    (eta > -1.0 && eta <= -0.5)  * (pt*cosh(eta) > 0.150 && pt*cosh(eta)<0.200)* (0.89) +
    (eta > -1.0 && eta <= -0.5)  * (pt*cosh(eta) > 0.200)                      * (0.98) +
    (eta > -0.5 && eta <= 0.0)   * (pt*cosh(eta) > 0.150 && pt*cosh(eta)<0.200)*(0.89) +
    (eta > -0.5 && eta <= 0.0)   * (pt*cosh(eta) > 0.200 )                     * (0.98) +
    (eta > 0.0 && eta <= 0.5 )   * (pt*cosh(eta) > 0.150 && pt*cosh(eta)<0.200)* (0.89) +
    (eta > 0.0 && eta <= 0.5 )   * (pt*cosh(eta) > 0.200)                      * (0.98) +
    (eta > 0.5 && eta <= 1.0 )   * (pt*cosh(eta) > 0.150 && pt*cosh(eta)<0.200)* (0.89) +
    (eta > 0.5 && eta <= 1.0 )   * (pt*cosh(eta) > 0.200)                      * (0.98) +
    (eta > 1.0 && eta <= 1.5)    * (pt*cosh(eta) > 0.150 && pt*cosh(eta)<0.200)* (0.86) +
    (eta > 1.0 && eta <= 1.5)    * (pt*cosh(eta) > 0.200)                      * (0.92) +
    (eta > 1.5 && eta <= 2.0)    * (pt*cosh(eta) > 0.250 && pt*cosh(eta)<0.500)* (0.89) +
    (eta > 1.5 && eta <= 2.0)    * (pt*cosh(eta) > 0.500)                      * (0.98) +
    (eta > 2.0 && eta <= 2.5)    * (pt*cosh(eta) > 0.350 && pt*cosh(eta)<0.700)* (0.88) +
    (eta > 2.0 && eta <= 2.5)    * (pt*cosh(eta) > 0.700)                      * (0.97) +
    (eta > 2.5 && eta <= 3.0)    * (pt*cosh(eta) > 0.550 && pt*cosh(eta)<2.0)  * (0.95) +
    (eta > 2.5 && eta <= 3.0)    * (pt*cosh(eta) > 2.0)                        * (0.87) +
    (eta > 3.0 && eta <= 3.5)    * (pt*cosh(eta) > 0.850 && pt*cosh(eta)<4.0)  * (0.87) +
    (eta > 3.0 && eta <= 3.5)    * (pt*cosh(eta) > 4.0)                        * (0.95) +
    (abs(eta) > 3.5)  * (0.00)+
    0.0
}

#ATHENA Hybrid design
set CommonTrackingResolution {
    (eta<=-3.0 && eta>-3.5)  * (sqrt( (1.841e-2)^2 + (pt*cosh(eta)*7.1e-4)^2  ) )  +
    (eta<=-2.5 && eta>-3.0)  * (sqrt( (1.080e-2)^2 + (pt*cosh(eta)*1.7e-4)^2  ) )  +
    (eta<=-2.0 && eta>-2.5)  * (sqrt( (6.33e-3)^2 + (pt*cosh(eta)*0.0)^2  ) )  +
    (eta<=-1.5 && eta>-2.0)  * (sqrt( (4.76e-3)^2 + (pt*cosh(eta)*1.1e-4)^2  ) )  +
    (eta<=-1.0 && eta>-1.5)  * (sqrt( (4.33e-3)^2 + (pt*cosh(eta)*1.6e-4)^2  ) )  +
    (eta<=-0.5 && eta>-1.0)  * (sqrt( (3.98e-3)^2 + (pt*cosh(eta)*5.0e-4)^2  ) )  +
    (eta<= 0.0 && eta>-0.5)  * (sqrt( (3.53e-3)^2 + (pt*cosh(eta)*5.9e-4)^2  ) )  +

    (eta<=0.5 && eta>0)  * (sqrt( (3.50e-3)^2 + (pt*cosh(eta)*5.9e-4)^2  ) )  +
    (eta<=1.0 && eta>0.5) * (sqrt( (4.01e-3)^2 + (pt*cosh(eta)*5.0e-4)^2   ) )  +
    (eta<=1.5 && eta>1.0) * (sqrt( (4.14e-3)^2 + (pt*cosh(eta)*1.5e-4)^2   ) )  +
    (eta<=2.0 && eta>1.5) * (sqrt( (4.66e-3)^2 + (pt*cosh(eta)*1.1e-4)^2   ) )  +
    (eta<=2.5 && eta>2.0) * (sqrt( (6.38e-3)^2 + (pt*cosh(eta)*1.3e-4)^2   ) )  +
    (eta<=3.0 && eta>2.5) * (sqrt( (1.089e-2)^2 + (pt*cosh(eta)*1.1e-4)^2   ) )  +
    (eta<=3.5 && eta>3.0) * (sqrt( (1.905e-2)^2 + (pt*cosh(eta)*3.1e-4)^2  ) )  
}


####################################
# Charged hadron tracking efficiency
####################################

module Efficiency ChargedHadronTrackingEfficiency {
  set InputArray ParticlePropagator/chargedHadrons
  set OutputArray chargedHadrons
  set EfficiencyFormula $CommonTrackingEfficiency
}

##############################
# Electron tracking efficiency
##############################

module Efficiency ElectronTrackingEfficiency {
  set InputArray ParticlePropagator/electrons
  set OutputArray electrons
  set EfficiencyFormula $CommonTrackingEfficiency

}

##############################
# Muon tracking efficiency
##############################
module Efficiency MuonTrackingEfficiency {
  set InputArray ParticlePropagator/muons
  set OutputArray muons
  set EfficiencyFormula $CommonTrackingEfficiency

}

########################################
# Smearing for charged hadrons
########################################

module TrackSmearing ChargedHadronSmearing {
  set InputArray ChargedHadronTrackingEfficiency/chargedHadrons
  set BeamSpotInputArray BeamSpotFilter/beamSpotParticle
  set OutputArray chargedHadrons
#  set ApplyToPileUp true
  # magnetic field
  set Bz 3.0
  set PResolutionFormula $CommonTrackingResolution
  set CtgThetaResolutionFormula { 0.0 }
  set PhiResolutionFormula { 0.0 }

    

  set D0ResolutionFormula "
    (eta<=0.0 && eta>-0.5)    * (sqrt( (0.00541)^2 +   (0.02679/(pt*cosh(eta)))^2   ) )  +
    (eta<=-0.5 && eta>-1.0)   * (sqrt( (0.00545)^2 +   (0.02960/(pt*cosh(eta)))^2   ) )  +
    (eta<=-1.0 && eta>-1.5)   * (sqrt( (0.00577)^2 +   (0.03561/(pt*cosh(eta)))^2   ) )  +
    (eta<=-1.5 && eta>-2.0)   * (sqrt( (0.00550)^2 +   (0.04260/(pt*cosh(eta)))^2   ) )  +
    (eta<=-2.0 && eta>-2.5)   * (sqrt( (0.00596)^2 +   (0.05891/(pt*cosh(eta)))^2   ) )  +
    (eta<=-2.5 && eta>-3.0)   * (sqrt( (0.00773)^2 +   (0.07631/(pt*cosh(eta)))^2   ) )  +
    (eta<=-3.0 && eta>-3.5)   * (sqrt( (0.01392)^2 +   (0.09182/(pt*cosh(eta)))^2   ) )  +

    (eta<=0.5 && eta>0)     * (sqrt( (0.00536)^2 +   (0.02668/(pt*cosh(eta)))^2   ) )  +                                                                                                   
    (eta<=1.0 && eta>0.5)   * (sqrt( (0.00550)^2 +   (0.02951/(pt*cosh(eta)))^2   ) )  +                                                                                                              
    (eta<=1.5 && eta>1.0)   * (sqrt( (0.00570)^2 +   (0.03565/(pt*cosh(eta)))^2   ) )  +                                                                                                             
   
    (eta<=2.0 && eta>1.5)   * (sqrt( (0.00539)^2 +   (0.04250/(pt*cosh(eta)))^2   ) )  +                                                                                                                 
    (eta<=2.5 && eta>2.0)   * (sqrt( (0.00588)^2 +   (0.05919/(pt*cosh(eta)))^2   ) )  +                                                                                                              
    (eta<=3.0 && eta>2.5)   * (sqrt( (0.00650)^2 +   (0.07622/(pt*cosh(eta)))^2   ) )  +                                                                                                               
    (eta<=3.5 && eta>3.0)   * (sqrt( (0.01119)^2 +   (0.09137/(pt*cosh(eta)))^2   ) )  
  "


  set DZResolutionFormula "
    (abs(eta)<=0.5)                   * (sqrt( (0.0034)^2 +   (0.027/(pt*cosh(eta)))^2   ) )  +
    (abs(eta)<=1.0 && abs(eta)>0.5)   * (sqrt( (0.0038)^2 +   (0.036/(pt*cosh(eta)))^2   ) )  +
    (abs(eta)<=1.5 && abs(eta)>1.0)   * (sqrt( (0.0056)^2 +   (0.061/(pt*cosh(eta)))^2   ) )  +
    (abs(eta)<=2.0 && abs(eta)>1.5)   * (sqrt( (0.0072)^2 +   (0.116/(pt*cosh(eta)))^2   ) )  +   
    (abs(eta)<=2.5 && abs(eta)>2.0)   * (sqrt( (0.0095)^2 +   (0.244/(pt*cosh(eta)))^2   ) )  + 
    (abs(eta)<=3.0 && abs(eta)>2.5)   * (sqrt( (0.0330)^2 +   (2.581/(pt*cosh(eta)))^2   ) )  +
    (abs(eta)<=3.5 && abs(eta)>3.0)   * (sqrt( (0.1890)^2 +   (8.349/(pt*cosh(eta)))^2   ) )  
  "


}

###################################
# Smearing for muons
###################################

module TrackSmearing MuonSmearing {
  set InputArray MuonTrackingEfficiency/muons
  set BeamSpotInputArray BeamSpotFilter/beamSpotParticle
  set OutputArray muons
#  set ApplyToPileUp true
  # magnetic field
  set Bz 3.0
  set PResolutionFormula $CommonTrackingResolution
  set CtgThetaResolutionFormula { 0.0 }
  set PhiResolutionFormula { 0.0 }

  set D0ResolutionFormula  "                                                                                                                                                                              
    (eta<=0.0 && eta>-0.5)    * (sqrt( (0.00541)^2 +   (0.02679/(pt*cosh(eta)))^2   ) )  +                                                                                                                
    (eta<=-0.5 && eta>-1.0)   * (sqrt( (0.00545)^2 +   (0.02960/(pt*cosh(eta)))^2   ) )  +                                                                                                               
    (eta<=-1.0 && eta>-1.5)   * (sqrt( (0.00577)^2 +   (0.03561/(pt*cosh(eta)))^2   ) )  +                                                                                                                
    (eta<=-1.5 && eta>-2.0)   * (sqrt( (0.00550)^2 +   (0.04260/(pt*cosh(eta)))^2   ) )  +                                                                                                               
    (eta<=-2.0 && eta>-2.5)   * (sqrt( (0.00596)^2 +   (0.05891/(pt*cosh(eta)))^2   ) )  +                                                                                                                
    (eta<=-2.5 && eta>-3.0)   * (sqrt( (0.00773)^2 +   (0.07631/(pt*cosh(eta)))^2   ) )  +                                                                                                               
    (eta<=-3.0 && eta>-3.5)   * (sqrt( (0.01392)^2 +   (0.09182/(pt*cosh(eta)))^2   ) )  +                                                                                                              
                                                                                                                                                                                                    
    (eta<=0.5 && eta>0)     * (sqrt( (0.00536)^2 +   (0.02668/(pt*cosh(eta)))^2   ) )  +                                                                                                                  
    (eta<=1.0 && eta>0.5)   * (sqrt( (0.00550)^2 +   (0.02951/(pt*cosh(eta)))^2   ) )  +                                                                                                                
    (eta<=1.5 && eta>1.0)   * (sqrt( (0.00570)^2 +   (0.03565/(pt*cosh(eta)))^2   ) )  +                                                                                                                 
    (eta<=2.0 && eta>1.5)   * (sqrt( (0.00539)^2 +   (0.04250/(pt*cosh(eta)))^2   ) )  +                                                                                                                 
    (eta<=2.5 && eta>2.0)   * (sqrt( (0.00588)^2 +   (0.05919/(pt*cosh(eta)))^2   ) )  +                                                                                                                
    (eta<=3.0 && eta>2.5)   * (sqrt( (0.00650)^2 +   (0.07622/(pt*cosh(eta)))^2   ) )  +                                                                                                                
    (eta<=3.5 && eta>3.0)   * (sqrt( (0.01119)^2 +   (0.09137/(pt*cosh(eta)))^2   ) )                                                                                                                     
  "

  

set DZResolutionFormula "                                                                                                                                                                                
    (abs(eta)<=0.5)                   * (sqrt( (0.0034)^2 +   (0.027/(pt*cosh(eta)))^2   ) )  +                                                                                                           
    (abs(eta)<=1.0 && abs(eta)>0.5)   * (sqrt( (0.0038)^2 +   (0.036/(pt*cosh(eta)))^2   ) )  +                                                                                                          
    (abs(eta)<=1.5 && abs(eta)>1.0)   * (sqrt( (0.0056)^2 +   (0.061/(pt*cosh(eta)))^2   ) )  +                                                                                                           
    (abs(eta)<=2.0 && abs(eta)>1.5)   * (sqrt( (0.0072)^2 +   (0.116/(pt*cosh(eta)))^2   ) )  +                                                                                                           
    (abs(eta)<=2.5 && abs(eta)>2.0)   * (sqrt( (0.0095)^2 +   (0.244/(pt*cosh(eta)))^2   ) )  +                                                                                                           
    (abs(eta)<=3.0 && abs(eta)>2.5)   * (sqrt( (0.0330)^2 +   (2.581/(pt*cosh(eta)))^2   ) )  +                                                                                                           
    (abs(eta)<=3.5 && abs(eta)>3.0)   * (sqrt( (0.1890)^2 +   (8.349/(pt*cosh(eta)))^2   ) )                                                                                                              
  "


}

###################################
# Smearing for electrons
###################################


module TrackSmearing ElectronSmearing {
  set InputArray ElectronTrackingEfficiency/electrons
  set BeamSpotInputArray BeamSpotFilter/beamSpotParticle
  set OutputArray electrons
#  set ApplyToPileUp true
  # magnetic field
  set Bz 3.0
  set PResolutionFormula $CommonTrackingResolution
  set CtgThetaResolutionFormula { 0.0 }
  set PhiResolutionFormula { 0.0 }

  set D0ResolutionFormula  "                                                                                                                                                                              
    (eta<=0.0 && eta>-0.5)    * (sqrt( (0.00541)^2 +   (0.02679/(pt*cosh(eta)))^2   ) )  +                                                                                                                
    (eta<=-0.5 && eta>-1.0)   * (sqrt( (0.00545)^2 +   (0.02960/(pt*cosh(eta)))^2   ) )  +                                                                                                                
    (eta<=-1.0 && eta>-1.5)   * (sqrt( (0.00577)^2 +   (0.03561/(pt*cosh(eta)))^2   ) )  +                                                                                                                
    (eta<=-1.5 && eta>-2.0)   * (sqrt( (0.00550)^2 +   (0.04260/(pt*cosh(eta)))^2   ) )  +                                                                                                              
    (eta<=-2.0 && eta>-2.5)   * (sqrt( (0.00596)^2 +   (0.05891/(pt*cosh(eta)))^2   ) )  +                                                                                                                
    (eta<=-2.5 && eta>-3.0)   * (sqrt( (0.00773)^2 +   (0.07631/(pt*cosh(eta)))^2   ) )  +                                                                                                                
    (eta<=-3.0 && eta>-3.5)   * (sqrt( (0.01392)^2 +   (0.09182/(pt*cosh(eta)))^2   ) )  +                                                                                                               
                                                                                                                                                                                                          
    (eta<=0.5 && eta>0)     * (sqrt( (0.00536)^2 +   (0.02668/(pt*cosh(eta)))^2   ) )  +                                                                                                                  
    (eta<=1.0 && eta>0.5)   * (sqrt( (0.00550)^2 +   (0.02951/(pt*cosh(eta)))^2   ) )  +                                                                                                                  
    (eta<=1.5 && eta>1.0)   * (sqrt( (0.00570)^2 +   (0.03565/(pt*cosh(eta)))^2   ) )  +                                                                                                                  
    (eta<=2.0 && eta>1.5)   * (sqrt( (0.00539)^2 +   (0.04250/(pt*cosh(eta)))^2   ) )  +                                                                                                                  
    (eta<=2.5 && eta>2.0)   * (sqrt( (0.00588)^2 +   (0.05919/(pt*cosh(eta)))^2   ) )  +                                                                                                                  
    (eta<=3.0 && eta>2.5)   * (sqrt( (0.00650)^2 +   (0.07622/(pt*cosh(eta)))^2   ) )  +                                                                                                                  
    (eta<=3.5 && eta>3.0)   * (sqrt( (0.01119)^2 +   (0.09137/(pt*cosh(eta)))^2   ) )                                                                                                                     
  "



set DZResolutionFormula "                                                                                                                                                                                
    (abs(eta)<=0.5)                   * (sqrt( (0.0034)^2 +   (0.027/(pt*cosh(eta)))^2   ) )  +                                                                                                          
    (abs(eta)<=1.0 && abs(eta)>0.5)   * (sqrt( (0.0038)^2 +   (0.036/(pt*cosh(eta)))^2   ) )  +                                                                                                           
    (abs(eta)<=1.5 && abs(eta)>1.0)   * (sqrt( (0.0056)^2 +   (0.061/(pt*cosh(eta)))^2   ) )  +                                                                                                           
    (abs(eta)<=2.0 && abs(eta)>1.5)   * (sqrt( (0.0072)^2 +   (0.116/(pt*cosh(eta)))^2   ) )  +                                                                                                           
    (abs(eta)<=2.5 && abs(eta)>2.0)   * (sqrt( (0.0095)^2 +   (0.244/(pt*cosh(eta)))^2   ) )  +                                                                                                           
    (abs(eta)<=3.0 && abs(eta)>2.5)   * (sqrt( (0.0330)^2 +   (2.581/(pt*cosh(eta)))^2   ) )  +                                                                                                           
    (abs(eta)<=3.5 && abs(eta)>3.0)   * (sqrt( (0.1890)^2 +   (8.349/(pt*cosh(eta)))^2   ) )                                                                                                              
  "

}

##############
# Track merger
##############

module Merger TrackMerger {
# add InputArray InputArray
  add InputArray ChargedHadronSmearing/chargedHadrons
  add InputArray ElectronSmearing/electrons
  add InputArray MuonSmearing/muons
  set OutputArray tracks
}



#############
#   ECAL
#############

module SimpleCalorimeter ECal {
  set ParticleInputArray ParticlePropagator/stableParticles
  set TrackInputArray TrackMerger/tracks

  set TowerOutputArray ecalTowers
  set EFlowTrackOutputArray eflowTracks
  set EFlowTowerOutputArray eflowPhotons

  set IsEcal true
  set EnergyMin 0.50
  set EnergySignificanceMin 1.0

  set SmearTowerCenter true

  set pi [expr {acos(-1)}]

  # lists of the edges of each tower in eta and phi
  # each list starts with the lower edge of the first tower
  # the list ends with the higher edged of the last tower

    set PhiBins {}
    for {set i -30} {$i <=30} {incr i} {
	add PhiBins [expr {$i * $pi/30.0}]
    }
    for {set i -10} {$i <=10} {incr i} {
	set eta [expr {$i * 0.1}]
	add EtaPhiBins $eta $PhiBins
    }

    
   ## assume 0.1 x 0.1 (real cell size will be smaller, so this is to represent some cluster)
    set PhiBins {}
    for {set i -30} {$i <=30} {incr i} {
	add PhiBins [expr {$i * $pi/30.0}]
    }

    foreach eta {-4.0 -3.77209196 -3.57533702 -3.41102102 -3.26996837 -3.14642305 -3.03653567 -2.93760447 -2.84766006 -2.76522251 -2.68915144 -2.61854952 -2.55269788 -2.49101173 -2.43300894 -2.3782873  -2.3265078  -2.27738197 -2.23066235 -2.18613503 -2.14361383 -2.10293569 -2.063957 -2.02655061 -1.99060337 -1.95601417 -1.92269228 -1.89055593 -1.8595312  -1.82955102 -1.80055436 -1.77248548 -1.74529337 -1.71893119 -1.69335587 -1.66852765 -1.64440978 -1.62096821 -1.59817135 -1.57598979 -1.55439612 -1.53336478 -1.51287184 -1.4928949  -1.47341295 -1.45440623 -1.43585618 -1.41774529 -1.40005705 -1.38277588 -1.36588703 -1.34937654 -1.33323117 -1.31743839 -1.30198626 -1.28686345 -1.27205918 -1.25756317 -1.24336562 -1.22945719 -1.21582897 -1.20247241 -1.18937936 -1.17654201 -1.16395288 -1.15160481 -1.13949092 -1.12760462 -1.11593955 -1.10448965 -1.09324904 -1.08221211	-1.07137342 -1.06072776 -1.0502701  -1.03999558} {
	add EtaPhiBins $eta $PhiBins
    }
    foreach eta {1.0 1.0502701  1.06072776 1.07137342 1.08221211 1.09324904 1.10448965 1.11593955 1.12760462 1.13949092 1.15160481 1.16395288 1.17654201 1.18937936 1.20247241 1.21582897 1.22945719 1.24336562 1.25756317 1.27205918 1.28686345 1.30198626 1.31743839 1.33323117 1.34937654 1.36588703 1.38277588 1.40005705 1.41774529 1.43585618 1.45440623 1.47341295 1.4928949  1.51287184 1.53336478 1.55439612 1.57598979 1.59817135 1.62096821 1.64440978 1.66852765 1.69335587 1.71893119 1.74529337 1.77248548 1.80055436 1.82955102 1.8595312 1.89055593 1.92269228 1.95601417 1.99060337 2.02655061 2.063957 2.10293569 2.14361383 2.18613503 2.23066235 2.27738197 2.3265078 2.3782873  2.43300894 2.49101173 2.55269788 2.61854952 2.68915144 2.76522251 2.84766006 2.93760447 3.03653567 3.14642305 3.26996837 3.41102102 3.57533702 3.77209196 4.0} {
   	add EtaPhiBins $eta $PhiBins
    }
    
  add EnergyFraction {0} {0.0}
  # energy fractions for e, gamma and pi0
  add EnergyFraction {11} {1.0}
  add EnergyFraction {22} {1.0}
  add EnergyFraction {111} {1.0}
  # energy fractions for muon, neutrinos and neutralinos
  add EnergyFraction {12} {0.0}
  add EnergyFraction {13} {0.0}
  add EnergyFraction {14} {0.0}
  add EnergyFraction {16} {0.0}
  add EnergyFraction {1000022} {0.0}
  add EnergyFraction {1000023} {0.0}
  add EnergyFraction {1000025} {0.0}
  add EnergyFraction {1000035} {0.0}
  add EnergyFraction {1000045} {0.0}
  # energy fractions for K0short and Lambda
  add EnergyFraction {310} {0.3}
  add EnergyFraction {3122} {0.3}

  #Crystal -4.0 to -2.3
  #SciGlass -1.5 to -2.3
  #Barrel ECAL   -1.5 to 1.2
  #pECAL 1.2 to 4.0  
  set ResolutionFormula {          (eta <= -2.3 && eta>-4.0)                          * sqrt(energy^2*0.01^2 + energy*0.020^2 )+ \
                 		   (eta <= -1.5 && eta>-2.3 )                         * sqrt(energy^2*0.02^2 + energy*0.025^2 )+ \
				   (eta <= 1.2  && eta> -1.5 )                        * sqrt(energy^2*0.01^2 + energy*0.055^2 )+ \
				   (eta <= 4.0  &&  eta>1.2 )                         * sqrt(energy^2*0.02^2 + energy*0.100^2 )}

}


#############
#   HCAL
#############

module SimpleCalorimeter HCal {
  set ParticleInputArray ParticlePropagator/stableParticles
  set TrackInputArray ECal/eflowTracks

  set TowerOutputArray hcalTowers
  set EFlowTrackOutputArray eflowTracks
  set EFlowTowerOutputArray eflowNeutralHadrons

  set IsEcal false

  set EnergyMin 0.300
  set EnergySignificanceMin 1.0

  set SmearTowerCenter true

  set pi [expr {acos(-1)}]

    set PhiBins {}
    for {set i -30} {$i <=30} {incr i} {
	add PhiBins [expr {$i * $pi/30.0}]
    }
    for {set i -10} {$i <=10} {incr i} {
	set eta [expr {$i * 0.1}]
	add EtaPhiBins $eta $PhiBins
    }

    for {set i -30} {$i <=30} {incr i} {
	add PhiBins [expr {$i * $pi/30.0}]
    }
    foreach eta {-4.0 -3.34390825 -2.95880652 -2.68264484 -2.46773612 -2.29224349 -2.14432155 -2.01681569 -1.90506801 -1.80587261 -1.71692581 -1.63651428 -1.56332731 -1.49633825 -1.43472677 -1.37782606 -1.325086   -1.27604684 -1.23031998 -1.18757364 -1.14752205 -1.10991713 -1.07454199 -1.04120583 -1.00} {
	add EtaPhiBins $eta $PhiBins
    }
    
    foreach eta {1.0 1.04 1.075 1.1099 1.14752205 1.18757364 1.23031998 1.27604684 1.325086 1.37782606 1.43472677 1.49633825 1.56332731 1.63651428 1.71692581 1.80587261 1.90506801 2.01681569 2.14432155 2.29224349 2.46773612 2.68264484 2.95880652 3.34390825 4.00} {
	add EtaPhiBins $eta $PhiBins
    }


  add EnergyFraction {0} {1.0}
  # energy fractions for e, gamma and pi0
  add EnergyFraction {11} {0.0}
  add EnergyFraction {22} {0.0}
  add EnergyFraction {111} {0.0}
  # energy fractions for muon, neutrinos and neutralinos
  add EnergyFraction {12} {0.0}
  add EnergyFraction {13} {0.0}
  add EnergyFraction {14} {0.0}
  add EnergyFraction {16} {0.0}
  add EnergyFraction {1000022} {0.0}
  add EnergyFraction {1000023} {0.0}
  add EnergyFraction {1000025} {0.0}
  add EnergyFraction {1000035} {0.0}
  add EnergyFraction {1000045} {0.0}
  # energy fractions for K0short and Lambda
  add EnergyFraction {310} {0.7}
  add EnergyFraction {3122} {0.7}

  # set HCalResolutionFormula {resolution formula as a function of eta and energy}
  set ResolutionFormula {    (eta <= -1.0 && eta>-4.0)                       * sqrt(energy^2*0.05^2 + energy*0.70^2)+
                             (eta <= 1.0 && eta>-1.0 )                       * sqrt(energy^2*0.37^2 )+
                             (eta <= 4.0 && eta>1.0 )                        * sqrt(energy^2*0.06^2 + energy*0.40^2)
  }

}


#################
# Electron filter
#################

module PdgCodeFilter ElectronFilter {
  set InputArray HCal/eflowTracks
  set OutputArray electrons
  set Invert true
  add PdgCode {11}
  add PdgCode {-11}
}

######################
# ChargedHadronFilter
######################

module PdgCodeFilter ChargedHadronFilter {
  set InputArray HCal/eflowTracks
  set OutputArray chargedHadrons

  add PdgCode {11}
  add PdgCode {-11}
  add PdgCode {13}
  add PdgCode {-13}
}


###################################################
# Tower Merger (in case not using e-flow algorithm)
###################################################

module Merger Calorimeter {
# add InputArray InputArray
  add InputArray ECal/ecalTowers
  add InputArray HCal/hcalTowers
  set OutputArray towers
}



####################
# Energy flow merger
####################

module Merger EFlowMerger {
# add InputArray InputArray
  add InputArray HCal/eflowTracks
  add InputArray ECal/eflowPhotons
  add InputArray HCal/eflowNeutralHadrons
  set OutputArray eflow
}

######################
# EFlowFilter
######################

module PdgCodeFilter EFlowFilter {
  set InputArray EFlowMerger/eflow
  set OutputArray eflow

  add PdgCode {11}
  add PdgCode {-11}
  add PdgCode {13}
  add PdgCode {-13}
}

###################
# Photon efficiency
###################

module Efficiency PhotonEfficiency {
  set InputArray ECal/eflowPhotons
  set OutputArray photons

  # set EfficiencyFormula {efficiency formula as a function of eta and pt}

  # efficiency formula for photons
    set EfficiencyFormula { 1}
}

##################
# Photon isolation
##################

module Isolation PhotonIsolation {
  set CandidateInputArray PhotonEfficiency/photons
  set IsolationInputArray EFlowFilter/eflow

  set OutputArray photons

  set DeltaRMax 0.5

  set PTMin 0.5

  set PTRatioMax 0.12
}


#####################
# Electron efficiency
#####################

module Efficiency ElectronEfficiency {
  set InputArray ElectronFilter/electrons
  set OutputArray electrons

  # set EfficiencyFormula {efficiency formula as a function of eta and pt}

  # efficiency formula for electrons
    set EfficiencyFormula {1}
}

####################
# Electron isolation
####################

module Isolation ElectronIsolation {
  set CandidateInputArray ElectronEfficiency/electrons
  set IsolationInputArray EFlowFilter/eflow

  set OutputArray electrons

  set DeltaRMax 0.5

  set PTMin 0.5

  set PTRatioMax 0.12
}

###################
# Missing ET merger
###################

module Merger MissingET {
# add InputArray InputArray
  add InputArray EFlowMerger/eflow
  set MomentumOutputArray momentum
}

##################
# Scalar HT merger
##################

module Merger ScalarHT {
# add InputArray InputArray
  add InputArray UniqueObjectFinder/jets
  add InputArray UniqueObjectFinder/electrons
  add InputArray UniqueObjectFinder/photons

  set EnergyOutputArray energy
}


#####################
# Neutrino Filter
#####################

module PdgCodeFilter NeutrinoFilter {

  set InputArray Delphes/stableParticles
  set OutputArray filteredParticles

  set PTMin 0.0

  add PdgCode {12}
  add PdgCode {14}
  add PdgCode {16}
  add PdgCode {-12}
  add PdgCode {-14}
  add PdgCode {-16}

}


#####################
# MC truth jet finder
#####################

module FastJetFinder GenJetFinder {
  set InputArray NeutrinoFilter/filteredParticles

  set OutputArray jets

  # algorithm: 1 CDFJetClu, 2 MidPoint, 3 SIScone, 4 kt, 5 Cambridge/Aachen, 6 antikt
  set JetAlgorithm 6
  set ParameterR 1.0

  set JetPTMin 3.0
}

#########################
# Gen Missing ET merger
########################

module Merger GenMissingET {
# add InputArray InputArray
  add InputArray NeutrinoFilter/filteredParticles
  set MomentumOutputArray momentum
}



############
# Jet finder
############

module FastJetFinder FastJetFinder {
#  set InputArray Calorimeter/towers
  set InputArray EFlowMerger/eflow

  set OutputArray jets

  # algorithm: 1 CDFJetClu, 2 MidPoint, 3 SIScone, 4 kt, 5 Cambridge/Aachen, 6 antikt
  set JetAlgorithm 6
  set ParameterR 1.0

  set ComputeNsubjettiness 1
  set Beta 1.0
  set AxisMode 4

  set ComputeTrimming 1
  set RTrim 0.4
  set PtFracTrim 0.20
  #set PtFracTrim 0.05

  set ComputePruning 1
  set ZcutPrun 0.1
  set RcutPrun 0.5
  set RPrun 0.8

  set ComputeSoftDrop 1
  set BetaSoftDrop 0.0
  set SymmetryCutSoftDrop 0.1
  set R0SoftDrop 0.8

  set JetPTMin 3.0}


##################
# Jet Energy Scale
##################

module EnergyScale JetEnergyScale {
  set InputArray FastJetFinder/jets
  set OutputArray jets

  # scale formula for jets (do not apply it)
  set ScaleFormula {1.0}
}

########################
# Jet Flavor Association
########################

module JetFlavorAssociation JetFlavorAssociation {

  set PartonInputArray Delphes/partons
  set ParticleInputArray Delphes/allParticles
  set ParticleLHEFInputArray Delphes/allParticlesLHEF
  set JetInputArray JetEnergyScale/jets

  set DeltaR 0.5
  set PartonPTMin 4.0
  set PartonEtaMax 4.0

}

module JetFlavorAssociation GenJetFlavorAssociation {

  set PartonInputArray Delphes/partons
  set ParticleInputArray Delphes/allParticles
  set ParticleLHEFInputArray Delphes/allParticlesLHEF
  set JetInputArray GenJetFinder/jets

  set DeltaR 0.5
  set PartonPTMin 1.0
  set PartonEtaMax 4.0

}



#####################################################
# Find uniquely identified photons/electrons/tau/jets
#####################################################

module UniqueObjectFinder UniqueObjectFinder {
# earlier arrays take precedence over later ones
# add InputArray InputArray OutputArray
  add InputArray PhotonIsolation/photons photons
  add InputArray ElectronIsolation/electrons electrons
  add InputArray JetEnergyScale/jets jets
}

############################
# b-tagging (track counting)
############################

module TrackCountingBTagging TrackCountingBTagging {
    set JetInputArray JetEnergyScale/jets
    set TrackInputArray HCal/eflowTracks
    set BitNumber 0
    # maximum distance between jet and track
    set DeltaR 0.5
    # minimum pt of tracks
    set TrackPtMin 1.0
    # maximum transverse impact parameter (in mm)
    set TrackIPMax 3
    # minimum ip significance for the track to be counted
    set SigMin 2.0
    set Use3D true
    # alternate setting for 2D IP (default)
    #  set SigMin 1.3
    #  set Use3D false
    # minimum number of tracks (high efficiency n=2, high purity n=3)
    #set Ntracks 3

}

##################
# Particle ID Systems
##################

source pfRICH_0.5mrad.tcl
source pfRICH_0.0mrad.tcl
#My name is "Barrel DIRC TR=0.5 [mrad] dT=0.1 ns QE = 27 %" and I am described as follows:
#    Eta coverage =  [-1,1]
#    Assumed time precision = 0.1 ns
#    Assumed track resolution = 0.5 mrad
#    Assumed quantum efficiency of the MCP-PMT = 27%

module IdentificationMap barrelDIRC {
  set InputArray TrackMerger/tracks
  set OutputArray tracks

  # --- kaons ---

  add EfficiencyFormula {321} {321} {
    (eta<-1.00 || eta>=1.00 || pt * cosh(eta) < 0.10 || pt * cosh(eta) >= 50.00) * ( 0.00 ) +
    (-1.00 <= eta && eta < -0.90) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) +
    (-1.00 <= eta && eta < -0.90) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 5.70) * (0.999363) +
    (-1.00 <= eta && eta < -0.90) * (5.70 <= pt * cosh(eta) && pt * cosh(eta) < 6.20) * (0.977469) +
    (-1.00 <= eta && eta < -0.90) * (6.20 <= pt * cosh(eta) && pt * cosh(eta) < 6.60) * (0.953824) +
    (-1.00 <= eta && eta < -0.90) * (6.60 <= pt * cosh(eta) && pt * cosh(eta) < 6.90) * (0.930396) +
    (-1.00 <= eta && eta < -0.90) * (6.90 <= pt * cosh(eta) && pt * cosh(eta) < 7.20) * (0.908335) +
    (-1.00 <= eta && eta < -0.90) * (7.20 <= pt * cosh(eta) && pt * cosh(eta) < 7.50) * (0.885172) +
    (-1.00 <= eta && eta < -0.90) * (7.50 <= pt * cosh(eta) && pt * cosh(eta) < 7.80) * (0.861520) +
    (-1.00 <= eta && eta < -0.90) * (7.80 <= pt * cosh(eta) && pt * cosh(eta) < 8.10) * (0.838567) +
    (-1.00 <= eta && eta < -0.90) * (8.10 <= pt * cosh(eta) && pt * cosh(eta) < 8.40) * (0.817634) +
    (-1.00 <= eta && eta < -0.90) * (8.40 <= pt * cosh(eta) && pt * cosh(eta) < 8.80) * (0.793978) +
    (-1.00 <= eta && eta < -0.90) * (8.80 <= pt * cosh(eta) && pt * cosh(eta) < 9.20) * (0.769187) +
    (-1.00 <= eta && eta < -0.90) * (9.20 <= pt * cosh(eta) && pt * cosh(eta) < 9.70) * (0.745592) +
    (-1.00 <= eta && eta < -0.90) * (9.70 <= pt * cosh(eta) && pt * cosh(eta) < 10.30) * (0.721481) +
    (-1.00 <= eta && eta < -0.90) * (10.30 <= pt * cosh(eta) && pt * cosh(eta) < 11.00) * (0.698201) +
    (-1.00 <= eta && eta < -0.90) * (11.00 <= pt * cosh(eta) && pt * cosh(eta) < 11.90) * (0.676260) +
    (-1.00 <= eta && eta < -0.90) * (11.90 <= pt * cosh(eta) && pt * cosh(eta) < 13.30) * (0.654166) +
    (-1.00 <= eta && eta < -0.90) * (13.30 <= pt * cosh(eta) && pt * cosh(eta) < 16.50) * (0.629819) +
    (-1.00 <= eta && eta < -0.90) * (16.50 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.604030) +
    (-0.90 <= eta && eta < -0.80) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) +
    (-0.90 <= eta && eta < -0.80) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 5.60) * (0.999415) +
    (-0.90 <= eta && eta < -0.80) * (5.60 <= pt * cosh(eta) && pt * cosh(eta) < 6.10) * (0.978320) +
    (-0.90 <= eta && eta < -0.80) * (6.10 <= pt * cosh(eta) && pt * cosh(eta) < 6.50) * (0.954946) +
    (-0.90 <= eta && eta < -0.80) * (6.50 <= pt * cosh(eta) && pt * cosh(eta) < 6.80) * (0.931253) +
    (-0.90 <= eta && eta < -0.80) * (6.80 <= pt * cosh(eta) && pt * cosh(eta) < 7.10) * (0.908742) +
    (-0.90 <= eta && eta < -0.80) * (7.10 <= pt * cosh(eta) && pt * cosh(eta) < 7.40) * (0.884654) +
    (-0.90 <= eta && eta < -0.80) * (7.40 <= pt * cosh(eta) && pt * cosh(eta) < 7.70) * (0.860913) +
    (-0.90 <= eta && eta < -0.80) * (7.70 <= pt * cosh(eta) && pt * cosh(eta) < 8.00) * (0.837747) +
    (-0.90 <= eta && eta < -0.80) * (8.00 <= pt * cosh(eta) && pt * cosh(eta) < 8.30) * (0.815568) +
    (-0.90 <= eta && eta < -0.80) * (8.30 <= pt * cosh(eta) && pt * cosh(eta) < 8.70) * (0.791569) +
    (-0.90 <= eta && eta < -0.80) * (8.70 <= pt * cosh(eta) && pt * cosh(eta) < 9.10) * (0.767848) +
    (-0.90 <= eta && eta < -0.80) * (9.10 <= pt * cosh(eta) && pt * cosh(eta) < 9.50) * (0.747406) +
    (-0.90 <= eta && eta < -0.80) * (9.50 <= pt * cosh(eta) && pt * cosh(eta) < 10.10) * (0.723516) +
    (-0.90 <= eta && eta < -0.80) * (10.10 <= pt * cosh(eta) && pt * cosh(eta) < 10.80) * (0.699697) +
    (-0.90 <= eta && eta < -0.80) * (10.80 <= pt * cosh(eta) && pt * cosh(eta) < 11.70) * (0.677096) +
    (-0.90 <= eta && eta < -0.80) * (11.70 <= pt * cosh(eta) && pt * cosh(eta) < 13.10) * (0.654458) +
    (-0.90 <= eta && eta < -0.80) * (13.10 <= pt * cosh(eta) && pt * cosh(eta) < 16.20) * (0.630047) +
    (-0.90 <= eta && eta < -0.80) * (16.20 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.603999) +
    (-0.80 <= eta && eta < -0.70) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) +
    (-0.80 <= eta && eta < -0.70) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 5.50) * (0.999438) +
    (-0.80 <= eta && eta < -0.70) * (5.50 <= pt * cosh(eta) && pt * cosh(eta) < 6.00) * (0.978782) +
    (-0.80 <= eta && eta < -0.70) * (6.00 <= pt * cosh(eta) && pt * cosh(eta) < 6.40) * (0.953835) +
    (-0.80 <= eta && eta < -0.70) * (6.40 <= pt * cosh(eta) && pt * cosh(eta) < 6.70) * (0.929657) +
    (-0.80 <= eta && eta < -0.70) * (6.70 <= pt * cosh(eta) && pt * cosh(eta) < 7.00) * (0.907134) +
    (-0.80 <= eta && eta < -0.70) * (7.00 <= pt * cosh(eta) && pt * cosh(eta) < 7.30) * (0.882151) +
    (-0.80 <= eta && eta < -0.70) * (7.30 <= pt * cosh(eta) && pt * cosh(eta) < 7.60) * (0.858132) +
    (-0.80 <= eta && eta < -0.70) * (7.60 <= pt * cosh(eta) && pt * cosh(eta) < 7.90) * (0.835697) +
    (-0.80 <= eta && eta < -0.70) * (7.90 <= pt * cosh(eta) && pt * cosh(eta) < 8.20) * (0.813827) +
    (-0.80 <= eta && eta < -0.70) * (8.20 <= pt * cosh(eta) && pt * cosh(eta) < 8.60) * (0.790188) +
    (-0.80 <= eta && eta < -0.70) * (8.60 <= pt * cosh(eta) && pt * cosh(eta) < 9.00) * (0.765614) +
    (-0.80 <= eta && eta < -0.70) * (9.00 <= pt * cosh(eta) && pt * cosh(eta) < 9.50) * (0.742033) +
    (-0.80 <= eta && eta < -0.70) * (9.50 <= pt * cosh(eta) && pt * cosh(eta) < 10.00) * (0.720207) +
    (-0.80 <= eta && eta < -0.70) * (10.00 <= pt * cosh(eta) && pt * cosh(eta) < 10.70) * (0.698161) +
    (-0.80 <= eta && eta < -0.70) * (10.70 <= pt * cosh(eta) && pt * cosh(eta) < 11.60) * (0.675683) +
    (-0.80 <= eta && eta < -0.70) * (11.60 <= pt * cosh(eta) && pt * cosh(eta) < 13.00) * (0.653272) +
    (-0.80 <= eta && eta < -0.70) * (13.00 <= pt * cosh(eta) && pt * cosh(eta) < 16.30) * (0.628628) +
    (-0.80 <= eta && eta < -0.70) * (16.30 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.603794) +
    (-0.70 <= eta && eta < -0.60) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) +
    (-0.70 <= eta && eta < -0.60) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 5.50) * (0.999386) +
    (-0.70 <= eta && eta < -0.60) * (5.50 <= pt * cosh(eta) && pt * cosh(eta) < 6.00) * (0.977686) +
    (-0.70 <= eta && eta < -0.60) * (6.00 <= pt * cosh(eta) && pt * cosh(eta) < 6.40) * (0.952568) +
    (-0.70 <= eta && eta < -0.60) * (6.40 <= pt * cosh(eta) && pt * cosh(eta) < 6.70) * (0.927697) +
    (-0.70 <= eta && eta < -0.60) * (6.70 <= pt * cosh(eta) && pt * cosh(eta) < 7.00) * (0.904394) +
    (-0.70 <= eta && eta < -0.60) * (7.00 <= pt * cosh(eta) && pt * cosh(eta) < 7.30) * (0.880318) +
    (-0.70 <= eta && eta < -0.60) * (7.30 <= pt * cosh(eta) && pt * cosh(eta) < 7.60) * (0.856544) +
    (-0.70 <= eta && eta < -0.60) * (7.60 <= pt * cosh(eta) && pt * cosh(eta) < 7.90) * (0.833257) +
    (-0.70 <= eta && eta < -0.60) * (7.90 <= pt * cosh(eta) && pt * cosh(eta) < 8.20) * (0.810933) +
    (-0.70 <= eta && eta < -0.60) * (8.20 <= pt * cosh(eta) && pt * cosh(eta) < 8.60) * (0.786986) +
    (-0.70 <= eta && eta < -0.60) * (8.60 <= pt * cosh(eta) && pt * cosh(eta) < 9.00) * (0.763520) +
    (-0.70 <= eta && eta < -0.60) * (9.00 <= pt * cosh(eta) && pt * cosh(eta) < 9.50) * (0.740261) +
    (-0.70 <= eta && eta < -0.60) * (9.50 <= pt * cosh(eta) && pt * cosh(eta) < 10.00) * (0.718074) +
    (-0.70 <= eta && eta < -0.60) * (10.00 <= pt * cosh(eta) && pt * cosh(eta) < 10.70) * (0.696228) +
    (-0.70 <= eta && eta < -0.60) * (10.70 <= pt * cosh(eta) && pt * cosh(eta) < 11.60) * (0.674142) +
    (-0.70 <= eta && eta < -0.60) * (11.60 <= pt * cosh(eta) && pt * cosh(eta) < 13.10) * (0.651446) +
    (-0.70 <= eta && eta < -0.60) * (13.10 <= pt * cosh(eta) && pt * cosh(eta) < 17.00) * (0.625653) +
    (-0.70 <= eta && eta < -0.60) * (17.00 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.603457) +
    (-0.60 <= eta && eta < -0.50) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) +
    (-0.60 <= eta && eta < -0.50) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 5.50) * (0.999412) +
    (-0.60 <= eta && eta < -0.50) * (5.50 <= pt * cosh(eta) && pt * cosh(eta) < 6.00) * (0.978142) +
    (-0.60 <= eta && eta < -0.50) * (6.00 <= pt * cosh(eta) && pt * cosh(eta) < 6.40) * (0.953041) +
    (-0.60 <= eta && eta < -0.50) * (6.40 <= pt * cosh(eta) && pt * cosh(eta) < 6.70) * (0.928653) +
    (-0.60 <= eta && eta < -0.50) * (6.70 <= pt * cosh(eta) && pt * cosh(eta) < 7.00) * (0.906049) +
    (-0.60 <= eta && eta < -0.50) * (7.00 <= pt * cosh(eta) && pt * cosh(eta) < 7.30) * (0.882734) +
    (-0.60 <= eta && eta < -0.50) * (7.30 <= pt * cosh(eta) && pt * cosh(eta) < 7.60) * (0.858632) +
    (-0.60 <= eta && eta < -0.50) * (7.60 <= pt * cosh(eta) && pt * cosh(eta) < 7.90) * (0.835440) +
    (-0.60 <= eta && eta < -0.50) * (7.90 <= pt * cosh(eta) && pt * cosh(eta) < 8.20) * (0.813533) +
    (-0.60 <= eta && eta < -0.50) * (8.20 <= pt * cosh(eta) && pt * cosh(eta) < 8.60) * (0.789624) +
    (-0.60 <= eta && eta < -0.50) * (8.60 <= pt * cosh(eta) && pt * cosh(eta) < 9.00) * (0.765618) +
    (-0.60 <= eta && eta < -0.50) * (9.00 <= pt * cosh(eta) && pt * cosh(eta) < 9.40) * (0.744122) +
    (-0.60 <= eta && eta < -0.50) * (9.40 <= pt * cosh(eta) && pt * cosh(eta) < 9.90) * (0.722552) +
    (-0.60 <= eta && eta < -0.50) * (9.90 <= pt * cosh(eta) && pt * cosh(eta) < 10.60) * (0.700157) +
    (-0.60 <= eta && eta < -0.50) * (10.60 <= pt * cosh(eta) && pt * cosh(eta) < 11.50) * (0.677085) +
    (-0.60 <= eta && eta < -0.50) * (11.50 <= pt * cosh(eta) && pt * cosh(eta) < 12.90) * (0.654123) +
    (-0.60 <= eta && eta < -0.50) * (12.90 <= pt * cosh(eta) && pt * cosh(eta) < 16.10) * (0.629276) +
    (-0.60 <= eta && eta < -0.50) * (16.10 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.603840) +
    (-0.50 <= eta && eta < -0.40) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) +
    (-0.50 <= eta && eta < -0.40) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 5.40) * (0.999451) +
    (-0.50 <= eta && eta < -0.40) * (5.40 <= pt * cosh(eta) && pt * cosh(eta) < 5.90) * (0.979067) +
    (-0.50 <= eta && eta < -0.40) * (5.90 <= pt * cosh(eta) && pt * cosh(eta) < 6.20) * (0.958187) +
    (-0.50 <= eta && eta < -0.40) * (6.20 <= pt * cosh(eta) && pt * cosh(eta) < 6.50) * (0.937021) +
    (-0.50 <= eta && eta < -0.40) * (6.50 <= pt * cosh(eta) && pt * cosh(eta) < 6.80) * (0.913728) +
    (-0.50 <= eta && eta < -0.40) * (6.80 <= pt * cosh(eta) && pt * cosh(eta) < 7.10) * (0.889767) +
    (-0.50 <= eta && eta < -0.40) * (7.10 <= pt * cosh(eta) && pt * cosh(eta) < 7.40) * (0.865592) +
    (-0.50 <= eta && eta < -0.40) * (7.40 <= pt * cosh(eta) && pt * cosh(eta) < 7.70) * (0.842135) +
    (-0.50 <= eta && eta < -0.40) * (7.70 <= pt * cosh(eta) && pt * cosh(eta) < 8.00) * (0.819614) +
    (-0.50 <= eta && eta < -0.40) * (8.00 <= pt * cosh(eta) && pt * cosh(eta) < 8.40) * (0.795484) +
    (-0.50 <= eta && eta < -0.40) * (8.40 <= pt * cosh(eta) && pt * cosh(eta) < 8.80) * (0.770150) +
    (-0.50 <= eta && eta < -0.40) * (8.80 <= pt * cosh(eta) && pt * cosh(eta) < 9.20) * (0.747397) +
    (-0.50 <= eta && eta < -0.40) * (9.20 <= pt * cosh(eta) && pt * cosh(eta) < 9.70) * (0.725362) +
    (-0.50 <= eta && eta < -0.40) * (9.70 <= pt * cosh(eta) && pt * cosh(eta) < 10.30) * (0.703515) +
    (-0.50 <= eta && eta < -0.40) * (10.30 <= pt * cosh(eta) && pt * cosh(eta) < 11.10) * (0.681967) +
    (-0.50 <= eta && eta < -0.40) * (11.10 <= pt * cosh(eta) && pt * cosh(eta) < 12.30) * (0.659657) +
    (-0.50 <= eta && eta < -0.40) * (12.30 <= pt * cosh(eta) && pt * cosh(eta) < 14.60) * (0.636152) +
    (-0.50 <= eta && eta < -0.40) * (14.60 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.604442) +
    (-0.40 <= eta && eta < -0.30) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) +
    (-0.40 <= eta && eta < -0.30) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 5.40) * (0.999472) +
    (-0.40 <= eta && eta < -0.30) * (5.40 <= pt * cosh(eta) && pt * cosh(eta) < 5.90) * (0.979617) +
    (-0.40 <= eta && eta < -0.30) * (5.90 <= pt * cosh(eta) && pt * cosh(eta) < 6.30) * (0.955630) +
    (-0.40 <= eta && eta < -0.30) * (6.30 <= pt * cosh(eta) && pt * cosh(eta) < 6.60) * (0.932152) +
    (-0.40 <= eta && eta < -0.30) * (6.60 <= pt * cosh(eta) && pt * cosh(eta) < 6.90) * (0.909173) +
    (-0.40 <= eta && eta < -0.30) * (6.90 <= pt * cosh(eta) && pt * cosh(eta) < 7.20) * (0.885363) +
    (-0.40 <= eta && eta < -0.30) * (7.20 <= pt * cosh(eta) && pt * cosh(eta) < 7.50) * (0.860582) +
    (-0.40 <= eta && eta < -0.30) * (7.50 <= pt * cosh(eta) && pt * cosh(eta) < 7.80) * (0.836574) +
    (-0.40 <= eta && eta < -0.30) * (7.80 <= pt * cosh(eta) && pt * cosh(eta) < 8.10) * (0.813641) +
    (-0.40 <= eta && eta < -0.30) * (8.10 <= pt * cosh(eta) && pt * cosh(eta) < 8.50) * (0.789231) +
    (-0.40 <= eta && eta < -0.30) * (8.50 <= pt * cosh(eta) && pt * cosh(eta) < 8.90) * (0.764564) +
    (-0.40 <= eta && eta < -0.30) * (8.90 <= pt * cosh(eta) && pt * cosh(eta) < 9.40) * (0.740178) +
    (-0.40 <= eta && eta < -0.30) * (9.40 <= pt * cosh(eta) && pt * cosh(eta) < 10.00) * (0.715912) +
    (-0.40 <= eta && eta < -0.30) * (10.00 <= pt * cosh(eta) && pt * cosh(eta) < 10.70) * (0.693289) +
    (-0.40 <= eta && eta < -0.30) * (10.70 <= pt * cosh(eta) && pt * cosh(eta) < 11.70) * (0.670721) +
    (-0.40 <= eta && eta < -0.30) * (11.70 <= pt * cosh(eta) && pt * cosh(eta) < 13.30) * (0.647669) +
    (-0.40 <= eta && eta < -0.30) * (13.30 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.605441) +
    (-0.30 <= eta && eta < -0.20) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) +
    (-0.30 <= eta && eta < -0.20) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 5.40) * (0.999491) +
    (-0.30 <= eta && eta < -0.20) * (5.40 <= pt * cosh(eta) && pt * cosh(eta) < 5.90) * (0.979938) +
    (-0.30 <= eta && eta < -0.20) * (5.90 <= pt * cosh(eta) && pt * cosh(eta) < 6.20) * (0.959867) +
    (-0.30 <= eta && eta < -0.20) * (6.20 <= pt * cosh(eta) && pt * cosh(eta) < 6.50) * (0.938816) +
    (-0.30 <= eta && eta < -0.20) * (6.50 <= pt * cosh(eta) && pt * cosh(eta) < 6.80) * (0.916434) +
    (-0.30 <= eta && eta < -0.20) * (6.80 <= pt * cosh(eta) && pt * cosh(eta) < 7.10) * (0.892350) +
    (-0.30 <= eta && eta < -0.20) * (7.10 <= pt * cosh(eta) && pt * cosh(eta) < 7.40) * (0.868247) +
    (-0.30 <= eta && eta < -0.20) * (7.40 <= pt * cosh(eta) && pt * cosh(eta) < 7.70) * (0.845108) +
    (-0.30 <= eta && eta < -0.20) * (7.70 <= pt * cosh(eta) && pt * cosh(eta) < 8.00) * (0.821653) +
    (-0.30 <= eta && eta < -0.20) * (8.00 <= pt * cosh(eta) && pt * cosh(eta) < 8.30) * (0.799669) +
    (-0.30 <= eta && eta < -0.20) * (8.30 <= pt * cosh(eta) && pt * cosh(eta) < 8.70) * (0.776530) +
    (-0.30 <= eta && eta < -0.20) * (8.70 <= pt * cosh(eta) && pt * cosh(eta) < 9.10) * (0.754283) +
    (-0.30 <= eta && eta < -0.20) * (9.10 <= pt * cosh(eta) && pt * cosh(eta) < 9.60) * (0.732042) +
    (-0.30 <= eta && eta < -0.20) * (9.60 <= pt * cosh(eta) && pt * cosh(eta) < 10.20) * (0.708971) +
    (-0.30 <= eta && eta < -0.20) * (10.20 <= pt * cosh(eta) && pt * cosh(eta) < 11.00) * (0.685968) +
    (-0.30 <= eta && eta < -0.20) * (11.00 <= pt * cosh(eta) && pt * cosh(eta) < 12.10) * (0.663406) +
    (-0.30 <= eta && eta < -0.20) * (12.10 <= pt * cosh(eta) && pt * cosh(eta) < 14.10) * (0.640318) +
    (-0.30 <= eta && eta < -0.20) * (14.10 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.604825) +
    (-0.20 <= eta && eta < -0.10) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) +
    (-0.20 <= eta && eta < -0.10) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 5.50) * (0.999335) +
    (-0.20 <= eta && eta < -0.10) * (5.50 <= pt * cosh(eta) && pt * cosh(eta) < 6.00) * (0.976774) +
    (-0.20 <= eta && eta < -0.10) * (6.00 <= pt * cosh(eta) && pt * cosh(eta) < 6.40) * (0.951314) +
    (-0.20 <= eta && eta < -0.10) * (6.40 <= pt * cosh(eta) && pt * cosh(eta) < 6.70) * (0.926371) +
    (-0.20 <= eta && eta < -0.10) * (6.70 <= pt * cosh(eta) && pt * cosh(eta) < 7.00) * (0.902907) +
    (-0.20 <= eta && eta < -0.10) * (7.00 <= pt * cosh(eta) && pt * cosh(eta) < 7.30) * (0.878307) +
    (-0.20 <= eta && eta < -0.10) * (7.30 <= pt * cosh(eta) && pt * cosh(eta) < 7.60) * (0.855267) +
    (-0.20 <= eta && eta < -0.10) * (7.60 <= pt * cosh(eta) && pt * cosh(eta) < 7.90) * (0.831969) +
    (-0.20 <= eta && eta < -0.10) * (7.90 <= pt * cosh(eta) && pt * cosh(eta) < 8.30) * (0.806045) +
    (-0.20 <= eta && eta < -0.10) * (8.30 <= pt * cosh(eta) && pt * cosh(eta) < 8.70) * (0.780629) +
    (-0.20 <= eta && eta < -0.10) * (8.70 <= pt * cosh(eta) && pt * cosh(eta) < 9.10) * (0.757736) +
    (-0.20 <= eta && eta < -0.10) * (9.10 <= pt * cosh(eta) && pt * cosh(eta) < 9.60) * (0.734043) +
    (-0.20 <= eta && eta < -0.10) * (9.60 <= pt * cosh(eta) && pt * cosh(eta) < 10.20) * (0.711369) +
    (-0.20 <= eta && eta < -0.10) * (10.20 <= pt * cosh(eta) && pt * cosh(eta) < 11.00) * (0.688232) +
    (-0.20 <= eta && eta < -0.10) * (11.00 <= pt * cosh(eta) && pt * cosh(eta) < 12.10) * (0.665130) +
    (-0.20 <= eta && eta < -0.10) * (12.10 <= pt * cosh(eta) && pt * cosh(eta) < 14.00) * (0.641972) +
    (-0.20 <= eta && eta < -0.10) * (14.00 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.605010) +
    (-0.10 <= eta && eta < -0.00) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) +
    (-0.10 <= eta && eta < -0.00) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 5.60) * (0.999436) +
    (-0.10 <= eta && eta < -0.00) * (5.60 <= pt * cosh(eta) && pt * cosh(eta) < 6.10) * (0.978964) +
    (-0.10 <= eta && eta < -0.00) * (6.10 <= pt * cosh(eta) && pt * cosh(eta) < 6.40) * (0.958369) +
    (-0.10 <= eta && eta < -0.00) * (6.40 <= pt * cosh(eta) && pt * cosh(eta) < 6.70) * (0.938389) +
    (-0.10 <= eta && eta < -0.00) * (6.70 <= pt * cosh(eta) && pt * cosh(eta) < 7.00) * (0.916194) +
    (-0.10 <= eta && eta < -0.00) * (7.00 <= pt * cosh(eta) && pt * cosh(eta) < 7.30) * (0.893059) +
    (-0.10 <= eta && eta < -0.00) * (7.30 <= pt * cosh(eta) && pt * cosh(eta) < 7.60) * (0.869718) +
    (-0.10 <= eta && eta < -0.00) * (7.60 <= pt * cosh(eta) && pt * cosh(eta) < 7.90) * (0.846143) +
    (-0.10 <= eta && eta < -0.00) * (7.90 <= pt * cosh(eta) && pt * cosh(eta) < 8.20) * (0.825154) +
    (-0.10 <= eta && eta < -0.00) * (8.20 <= pt * cosh(eta) && pt * cosh(eta) < 8.60) * (0.800818) +
    (-0.10 <= eta && eta < -0.00) * (8.60 <= pt * cosh(eta) && pt * cosh(eta) < 9.00) * (0.778480) +
    (-0.10 <= eta && eta < -0.00) * (9.00 <= pt * cosh(eta) && pt * cosh(eta) < 9.40) * (0.755699) +
    (-0.10 <= eta && eta < -0.00) * (9.40 <= pt * cosh(eta) && pt * cosh(eta) < 9.90) * (0.733208) +
    (-0.10 <= eta && eta < -0.00) * (9.90 <= pt * cosh(eta) && pt * cosh(eta) < 10.50) * (0.710395) +
    (-0.10 <= eta && eta < -0.00) * (10.50 <= pt * cosh(eta) && pt * cosh(eta) < 11.30) * (0.687940) +
    (-0.10 <= eta && eta < -0.00) * (11.30 <= pt * cosh(eta) && pt * cosh(eta) < 12.40) * (0.665437) +
    (-0.10 <= eta && eta < -0.00) * (12.40 <= pt * cosh(eta) && pt * cosh(eta) < 14.30) * (0.642619) +
    (-0.10 <= eta && eta < -0.00) * (14.30 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.605218) +
    (-0.00 <= eta && eta < 0.10) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) +
    (-0.00 <= eta && eta < 0.10) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 5.50) * (0.999487) +
    (-0.00 <= eta && eta < 0.10) * (5.50 <= pt * cosh(eta) && pt * cosh(eta) < 6.00) * (0.980086) +
    (-0.00 <= eta && eta < 0.10) * (6.00 <= pt * cosh(eta) && pt * cosh(eta) < 6.40) * (0.955878) +
    (-0.00 <= eta && eta < 0.10) * (6.40 <= pt * cosh(eta) && pt * cosh(eta) < 6.70) * (0.931602) +
    (-0.00 <= eta && eta < 0.10) * (6.70 <= pt * cosh(eta) && pt * cosh(eta) < 7.00) * (0.908797) +
    (-0.00 <= eta && eta < 0.10) * (7.00 <= pt * cosh(eta) && pt * cosh(eta) < 7.30) * (0.885676) +
    (-0.00 <= eta && eta < 0.10) * (7.30 <= pt * cosh(eta) && pt * cosh(eta) < 7.60) * (0.862280) +
    (-0.00 <= eta && eta < 0.10) * (7.60 <= pt * cosh(eta) && pt * cosh(eta) < 7.90) * (0.839159) +
    (-0.00 <= eta && eta < 0.10) * (7.90 <= pt * cosh(eta) && pt * cosh(eta) < 8.20) * (0.815846) +
    (-0.00 <= eta && eta < 0.10) * (8.20 <= pt * cosh(eta) && pt * cosh(eta) < 8.50) * (0.795129) +
    (-0.00 <= eta && eta < 0.10) * (8.50 <= pt * cosh(eta) && pt * cosh(eta) < 8.90) * (0.773860) +
    (-0.00 <= eta && eta < 0.10) * (8.90 <= pt * cosh(eta) && pt * cosh(eta) < 9.30) * (0.751980) +
    (-0.00 <= eta && eta < 0.10) * (9.30 <= pt * cosh(eta) && pt * cosh(eta) < 9.90) * (0.728013) +
    (-0.00 <= eta && eta < 0.10) * (9.90 <= pt * cosh(eta) && pt * cosh(eta) < 10.50) * (0.705518) +
    (-0.00 <= eta && eta < 0.10) * (10.50 <= pt * cosh(eta) && pt * cosh(eta) < 11.30) * (0.683919) +
    (-0.00 <= eta && eta < 0.10) * (11.30 <= pt * cosh(eta) && pt * cosh(eta) < 12.50) * (0.661467) +
    (-0.00 <= eta && eta < 0.10) * (12.50 <= pt * cosh(eta) && pt * cosh(eta) < 14.70) * (0.638022) +
    (-0.00 <= eta && eta < 0.10) * (14.70 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.604731) +
    (0.10 <= eta && eta < 0.20) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) +
    (0.10 <= eta && eta < 0.20) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 5.40) * (0.999390) +
    (0.10 <= eta && eta < 0.20) * (5.40 <= pt * cosh(eta) && pt * cosh(eta) < 5.90) * (0.978436) +
    (0.10 <= eta && eta < 0.20) * (5.90 <= pt * cosh(eta) && pt * cosh(eta) < 6.30) * (0.954026) +
    (0.10 <= eta && eta < 0.20) * (6.30 <= pt * cosh(eta) && pt * cosh(eta) < 6.60) * (0.928980) +
    (0.10 <= eta && eta < 0.20) * (6.60 <= pt * cosh(eta) && pt * cosh(eta) < 6.90) * (0.905261) +
    (0.10 <= eta && eta < 0.20) * (6.90 <= pt * cosh(eta) && pt * cosh(eta) < 7.20) * (0.880827) +
    (0.10 <= eta && eta < 0.20) * (7.20 <= pt * cosh(eta) && pt * cosh(eta) < 7.50) * (0.855798) +
    (0.10 <= eta && eta < 0.20) * (7.50 <= pt * cosh(eta) && pt * cosh(eta) < 7.80) * (0.831816) +
    (0.10 <= eta && eta < 0.20) * (7.80 <= pt * cosh(eta) && pt * cosh(eta) < 8.10) * (0.811017) +
    (0.10 <= eta && eta < 0.20) * (8.10 <= pt * cosh(eta) && pt * cosh(eta) < 8.40) * (0.790324) +
    (0.10 <= eta && eta < 0.20) * (8.40 <= pt * cosh(eta) && pt * cosh(eta) < 8.80) * (0.766845) +
    (0.10 <= eta && eta < 0.20) * (8.80 <= pt * cosh(eta) && pt * cosh(eta) < 9.30) * (0.742732) +
    (0.10 <= eta && eta < 0.20) * (9.30 <= pt * cosh(eta) && pt * cosh(eta) < 9.80) * (0.721120) +
    (0.10 <= eta && eta < 0.20) * (9.80 <= pt * cosh(eta) && pt * cosh(eta) < 10.50) * (0.698097) +
    (0.10 <= eta && eta < 0.20) * (10.50 <= pt * cosh(eta) && pt * cosh(eta) < 11.40) * (0.675216) +
    (0.10 <= eta && eta < 0.20) * (11.40 <= pt * cosh(eta) && pt * cosh(eta) < 12.80) * (0.652610) +
    (0.10 <= eta && eta < 0.20) * (12.80 <= pt * cosh(eta) && pt * cosh(eta) < 16.20) * (0.627707) +
    (0.10 <= eta && eta < 0.20) * (16.20 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.603627) +
    (0.20 <= eta && eta < 0.30) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) +
    (0.20 <= eta && eta < 0.30) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 5.30) * (0.999486) +
    (0.20 <= eta && eta < 0.30) * (5.30 <= pt * cosh(eta) && pt * cosh(eta) < 5.80) * (0.979519) +
    (0.20 <= eta && eta < 0.30) * (5.80 <= pt * cosh(eta) && pt * cosh(eta) < 6.20) * (0.955406) +
    (0.20 <= eta && eta < 0.30) * (6.20 <= pt * cosh(eta) && pt * cosh(eta) < 6.50) * (0.930947) +
    (0.20 <= eta && eta < 0.30) * (6.50 <= pt * cosh(eta) && pt * cosh(eta) < 6.80) * (0.906986) +
    (0.20 <= eta && eta < 0.30) * (6.80 <= pt * cosh(eta) && pt * cosh(eta) < 7.10) * (0.882189) +
    (0.20 <= eta && eta < 0.30) * (7.10 <= pt * cosh(eta) && pt * cosh(eta) < 7.40) * (0.857189) +
    (0.20 <= eta && eta < 0.30) * (7.40 <= pt * cosh(eta) && pt * cosh(eta) < 7.70) * (0.833268) +
    (0.20 <= eta && eta < 0.30) * (7.70 <= pt * cosh(eta) && pt * cosh(eta) < 8.00) * (0.811348) +
    (0.20 <= eta && eta < 0.30) * (8.00 <= pt * cosh(eta) && pt * cosh(eta) < 8.40) * (0.787352) +
    (0.20 <= eta && eta < 0.30) * (8.40 <= pt * cosh(eta) && pt * cosh(eta) < 8.80) * (0.762475) +
    (0.20 <= eta && eta < 0.30) * (8.80 <= pt * cosh(eta) && pt * cosh(eta) < 9.30) * (0.737585) +
    (0.20 <= eta && eta < 0.30) * (9.30 <= pt * cosh(eta) && pt * cosh(eta) < 9.80) * (0.716144) +
    (0.20 <= eta && eta < 0.30) * (9.80 <= pt * cosh(eta) && pt * cosh(eta) < 10.50) * (0.693021) +
    (0.20 <= eta && eta < 0.30) * (10.50 <= pt * cosh(eta) && pt * cosh(eta) < 11.50) * (0.670061) +
    (0.20 <= eta && eta < 0.30) * (11.50 <= pt * cosh(eta) && pt * cosh(eta) < 13.10) * (0.646891) +
    (0.20 <= eta && eta < 0.30) * (13.10 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.605262) +
    (0.30 <= eta && eta < 0.40) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) +
    (0.30 <= eta && eta < 0.40) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 5.30) * (0.999459) +
    (0.30 <= eta && eta < 0.40) * (5.30 <= pt * cosh(eta) && pt * cosh(eta) < 5.80) * (0.979305) +
    (0.30 <= eta && eta < 0.40) * (5.80 <= pt * cosh(eta) && pt * cosh(eta) < 6.20) * (0.955111) +
    (0.30 <= eta && eta < 0.40) * (6.20 <= pt * cosh(eta) && pt * cosh(eta) < 6.50) * (0.929350) +
    (0.30 <= eta && eta < 0.40) * (6.50 <= pt * cosh(eta) && pt * cosh(eta) < 6.80) * (0.904229) +
    (0.30 <= eta && eta < 0.40) * (6.80 <= pt * cosh(eta) && pt * cosh(eta) < 7.10) * (0.880114) +
    (0.30 <= eta && eta < 0.40) * (7.10 <= pt * cosh(eta) && pt * cosh(eta) < 7.40) * (0.855471) +
    (0.30 <= eta && eta < 0.40) * (7.40 <= pt * cosh(eta) && pt * cosh(eta) < 7.70) * (0.831596) +
    (0.30 <= eta && eta < 0.40) * (7.70 <= pt * cosh(eta) && pt * cosh(eta) < 8.00) * (0.808837) +
    (0.30 <= eta && eta < 0.40) * (8.00 <= pt * cosh(eta) && pt * cosh(eta) < 8.40) * (0.785228) +
    (0.30 <= eta && eta < 0.40) * (8.40 <= pt * cosh(eta) && pt * cosh(eta) < 8.90) * (0.756982) +
    (0.30 <= eta && eta < 0.40) * (8.90 <= pt * cosh(eta) && pt * cosh(eta) < 9.40) * (0.732645) +
    (0.30 <= eta && eta < 0.40) * (9.40 <= pt * cosh(eta) && pt * cosh(eta) < 10.00) * (0.709880) +
    (0.30 <= eta && eta < 0.40) * (10.00 <= pt * cosh(eta) && pt * cosh(eta) < 10.80) * (0.687128) +
    (0.30 <= eta && eta < 0.40) * (10.80 <= pt * cosh(eta) && pt * cosh(eta) < 11.90) * (0.663925) +
    (0.30 <= eta && eta < 0.40) * (11.90 <= pt * cosh(eta) && pt * cosh(eta) < 13.90) * (0.640354) +
    (0.30 <= eta && eta < 0.40) * (13.90 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.604740) +
    (0.40 <= eta && eta < 0.50) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) +
    (0.40 <= eta && eta < 0.50) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 5.30) * (0.999416) +
    (0.40 <= eta && eta < 0.50) * (5.30 <= pt * cosh(eta) && pt * cosh(eta) < 5.80) * (0.977866) +
    (0.40 <= eta && eta < 0.50) * (5.80 <= pt * cosh(eta) && pt * cosh(eta) < 6.10) * (0.956023) +
    (0.40 <= eta && eta < 0.50) * (6.10 <= pt * cosh(eta) && pt * cosh(eta) < 6.40) * (0.934804) +
    (0.40 <= eta && eta < 0.50) * (6.40 <= pt * cosh(eta) && pt * cosh(eta) < 6.70) * (0.910945) +
    (0.40 <= eta && eta < 0.50) * (6.70 <= pt * cosh(eta) && pt * cosh(eta) < 7.00) * (0.886157) +
    (0.40 <= eta && eta < 0.50) * (7.00 <= pt * cosh(eta) && pt * cosh(eta) < 7.30) * (0.861028) +
    (0.40 <= eta && eta < 0.50) * (7.30 <= pt * cosh(eta) && pt * cosh(eta) < 7.60) * (0.836587) +
    (0.40 <= eta && eta < 0.50) * (7.60 <= pt * cosh(eta) && pt * cosh(eta) < 7.90) * (0.813125) +
    (0.40 <= eta && eta < 0.50) * (7.90 <= pt * cosh(eta) && pt * cosh(eta) < 8.30) * (0.788762) +
    (0.40 <= eta && eta < 0.50) * (8.30 <= pt * cosh(eta) && pt * cosh(eta) < 8.70) * (0.764916) +
    (0.40 <= eta && eta < 0.50) * (8.70 <= pt * cosh(eta) && pt * cosh(eta) < 9.20) * (0.739594) +
    (0.40 <= eta && eta < 0.50) * (9.20 <= pt * cosh(eta) && pt * cosh(eta) < 9.80) * (0.715349) +
    (0.40 <= eta && eta < 0.50) * (9.80 <= pt * cosh(eta) && pt * cosh(eta) < 10.50) * (0.691540) +
    (0.40 <= eta && eta < 0.50) * (10.50 <= pt * cosh(eta) && pt * cosh(eta) < 11.50) * (0.668971) +
    (0.40 <= eta && eta < 0.50) * (11.50 <= pt * cosh(eta) && pt * cosh(eta) < 13.20) * (0.645517) +
    (0.40 <= eta && eta < 0.50) * (13.20 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.605109) +
    (0.50 <= eta && eta < 0.60) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) +
    (0.50 <= eta && eta < 0.60) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 5.40) * (0.999374) +
    (0.50 <= eta && eta < 0.60) * (5.40 <= pt * cosh(eta) && pt * cosh(eta) < 5.90) * (0.976949) +
    (0.50 <= eta && eta < 0.60) * (5.90 <= pt * cosh(eta) && pt * cosh(eta) < 6.30) * (0.951356) +
    (0.50 <= eta && eta < 0.60) * (6.30 <= pt * cosh(eta) && pt * cosh(eta) < 6.60) * (0.926120) +
    (0.50 <= eta && eta < 0.60) * (6.60 <= pt * cosh(eta) && pt * cosh(eta) < 6.90) * (0.902566) +
    (0.50 <= eta && eta < 0.60) * (6.90 <= pt * cosh(eta) && pt * cosh(eta) < 7.20) * (0.878135) +
    (0.50 <= eta && eta < 0.60) * (7.20 <= pt * cosh(eta) && pt * cosh(eta) < 7.50) * (0.853567) +
    (0.50 <= eta && eta < 0.60) * (7.50 <= pt * cosh(eta) && pt * cosh(eta) < 7.80) * (0.830007) +
    (0.50 <= eta && eta < 0.60) * (7.80 <= pt * cosh(eta) && pt * cosh(eta) < 8.10) * (0.807715) +
    (0.50 <= eta && eta < 0.60) * (8.10 <= pt * cosh(eta) && pt * cosh(eta) < 8.50) * (0.784361) +
    (0.50 <= eta && eta < 0.60) * (8.50 <= pt * cosh(eta) && pt * cosh(eta) < 8.90) * (0.760491) +
    (0.50 <= eta && eta < 0.60) * (8.90 <= pt * cosh(eta) && pt * cosh(eta) < 9.40) * (0.736736) +
    (0.50 <= eta && eta < 0.60) * (9.40 <= pt * cosh(eta) && pt * cosh(eta) < 10.00) * (0.712414) +
    (0.50 <= eta && eta < 0.60) * (10.00 <= pt * cosh(eta) && pt * cosh(eta) < 10.70) * (0.690695) +
    (0.50 <= eta && eta < 0.60) * (10.70 <= pt * cosh(eta) && pt * cosh(eta) < 11.70) * (0.668691) +
    (0.50 <= eta && eta < 0.60) * (11.70 <= pt * cosh(eta) && pt * cosh(eta) < 13.40) * (0.645647) +
    (0.50 <= eta && eta < 0.60) * (13.40 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.605221) +
    (0.60 <= eta && eta < 0.70) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) +
    (0.60 <= eta && eta < 0.70) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 5.40) * (0.999428) +
    (0.60 <= eta && eta < 0.70) * (5.40 <= pt * cosh(eta) && pt * cosh(eta) < 5.90) * (0.978093) +
    (0.60 <= eta && eta < 0.70) * (5.90 <= pt * cosh(eta) && pt * cosh(eta) < 6.30) * (0.952896) +
    (0.60 <= eta && eta < 0.70) * (6.30 <= pt * cosh(eta) && pt * cosh(eta) < 6.60) * (0.928743) +
    (0.60 <= eta && eta < 0.70) * (6.60 <= pt * cosh(eta) && pt * cosh(eta) < 6.90) * (0.904835) +
    (0.60 <= eta && eta < 0.70) * (6.90 <= pt * cosh(eta) && pt * cosh(eta) < 7.20) * (0.880239) +
    (0.60 <= eta && eta < 0.70) * (7.20 <= pt * cosh(eta) && pt * cosh(eta) < 7.50) * (0.855800) +
    (0.60 <= eta && eta < 0.70) * (7.50 <= pt * cosh(eta) && pt * cosh(eta) < 7.80) * (0.832129) +
    (0.60 <= eta && eta < 0.70) * (7.80 <= pt * cosh(eta) && pt * cosh(eta) < 8.10) * (0.809971) +
    (0.60 <= eta && eta < 0.70) * (8.10 <= pt * cosh(eta) && pt * cosh(eta) < 8.50) * (0.785533) +
    (0.60 <= eta && eta < 0.70) * (8.50 <= pt * cosh(eta) && pt * cosh(eta) < 8.90) * (0.762384) +
    (0.60 <= eta && eta < 0.70) * (8.90 <= pt * cosh(eta) && pt * cosh(eta) < 9.40) * (0.739140) +
    (0.60 <= eta && eta < 0.70) * (9.40 <= pt * cosh(eta) && pt * cosh(eta) < 10.00) * (0.715270) +
    (0.60 <= eta && eta < 0.70) * (10.00 <= pt * cosh(eta) && pt * cosh(eta) < 10.70) * (0.692367) +
    (0.60 <= eta && eta < 0.70) * (10.70 <= pt * cosh(eta) && pt * cosh(eta) < 11.70) * (0.669998) +
    (0.60 <= eta && eta < 0.70) * (11.70 <= pt * cosh(eta) && pt * cosh(eta) < 13.30) * (0.647168) +
    (0.60 <= eta && eta < 0.70) * (13.30 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.605392) +
    (0.70 <= eta && eta < 0.80) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) +
    (0.70 <= eta && eta < 0.80) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 5.40) * (0.999431) +
    (0.70 <= eta && eta < 0.80) * (5.40 <= pt * cosh(eta) && pt * cosh(eta) < 5.90) * (0.978854) +
    (0.70 <= eta && eta < 0.80) * (5.90 <= pt * cosh(eta) && pt * cosh(eta) < 6.20) * (0.957614) +
    (0.70 <= eta && eta < 0.80) * (6.20 <= pt * cosh(eta) && pt * cosh(eta) < 6.50) * (0.936746) +
    (0.70 <= eta && eta < 0.80) * (6.50 <= pt * cosh(eta) && pt * cosh(eta) < 6.80) * (0.914572) +
    (0.70 <= eta && eta < 0.80) * (6.80 <= pt * cosh(eta) && pt * cosh(eta) < 7.10) * (0.889412) +
    (0.70 <= eta && eta < 0.80) * (7.10 <= pt * cosh(eta) && pt * cosh(eta) < 7.40) * (0.863632) +
    (0.70 <= eta && eta < 0.80) * (7.40 <= pt * cosh(eta) && pt * cosh(eta) < 7.70) * (0.839744) +
    (0.70 <= eta && eta < 0.80) * (7.70 <= pt * cosh(eta) && pt * cosh(eta) < 8.10) * (0.814966) +
    (0.70 <= eta && eta < 0.80) * (8.10 <= pt * cosh(eta) && pt * cosh(eta) < 8.50) * (0.787803) +
    (0.70 <= eta && eta < 0.80) * (8.50 <= pt * cosh(eta) && pt * cosh(eta) < 8.90) * (0.762097) +
    (0.70 <= eta && eta < 0.80) * (8.90 <= pt * cosh(eta) && pt * cosh(eta) < 9.40) * (0.738550) +
    (0.70 <= eta && eta < 0.80) * (9.40 <= pt * cosh(eta) && pt * cosh(eta) < 10.00) * (0.715206) +
    (0.70 <= eta && eta < 0.80) * (10.00 <= pt * cosh(eta) && pt * cosh(eta) < 10.70) * (0.692151) +
    (0.70 <= eta && eta < 0.80) * (10.70 <= pt * cosh(eta) && pt * cosh(eta) < 11.70) * (0.669829) +
    (0.70 <= eta && eta < 0.80) * (11.70 <= pt * cosh(eta) && pt * cosh(eta) < 13.30) * (0.647051) +
    (0.70 <= eta && eta < 0.80) * (13.30 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.605381) +
    (0.80 <= eta && eta < 0.90) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) +
    (0.80 <= eta && eta < 0.90) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 5.50) * (0.999463) +
    (0.80 <= eta && eta < 0.90) * (5.50 <= pt * cosh(eta) && pt * cosh(eta) < 6.00) * (0.979480) +
    (0.80 <= eta && eta < 0.90) * (6.00 <= pt * cosh(eta) && pt * cosh(eta) < 6.40) * (0.955658) +
    (0.80 <= eta && eta < 0.90) * (6.40 <= pt * cosh(eta) && pt * cosh(eta) < 6.70) * (0.931795) +
    (0.80 <= eta && eta < 0.90) * (6.70 <= pt * cosh(eta) && pt * cosh(eta) < 7.00) * (0.909021) +
    (0.80 <= eta && eta < 0.90) * (7.00 <= pt * cosh(eta) && pt * cosh(eta) < 7.30) * (0.884511) +
    (0.80 <= eta && eta < 0.90) * (7.30 <= pt * cosh(eta) && pt * cosh(eta) < 7.60) * (0.859876) +
    (0.80 <= eta && eta < 0.90) * (7.60 <= pt * cosh(eta) && pt * cosh(eta) < 7.90) * (0.837206) +
    (0.80 <= eta && eta < 0.90) * (7.90 <= pt * cosh(eta) && pt * cosh(eta) < 8.20) * (0.815556) +
    (0.80 <= eta && eta < 0.90) * (8.20 <= pt * cosh(eta) && pt * cosh(eta) < 8.60) * (0.792098) +
    (0.80 <= eta && eta < 0.90) * (8.60 <= pt * cosh(eta) && pt * cosh(eta) < 9.00) * (0.767572) +
    (0.80 <= eta && eta < 0.90) * (9.00 <= pt * cosh(eta) && pt * cosh(eta) < 9.50) * (0.742973) +
    (0.80 <= eta && eta < 0.90) * (9.50 <= pt * cosh(eta) && pt * cosh(eta) < 10.10) * (0.719838) +
    (0.80 <= eta && eta < 0.90) * (10.10 <= pt * cosh(eta) && pt * cosh(eta) < 10.80) * (0.696957) +
    (0.80 <= eta && eta < 0.90) * (10.80 <= pt * cosh(eta) && pt * cosh(eta) < 11.70) * (0.674906) +
    (0.80 <= eta && eta < 0.90) * (11.70 <= pt * cosh(eta) && pt * cosh(eta) < 13.10) * (0.652868) +
    (0.80 <= eta && eta < 0.90) * (13.10 <= pt * cosh(eta) && pt * cosh(eta) < 16.50) * (0.628233) +
    (0.80 <= eta && eta < 0.90) * (16.50 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.603774) +
    (0.90 <= eta && eta < 1.00) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) +
    (0.90 <= eta && eta < 1.00) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 5.50) * (0.999490) +
    (0.90 <= eta && eta < 1.00) * (5.50 <= pt * cosh(eta) && pt * cosh(eta) < 6.00) * (0.980505) +
    (0.90 <= eta && eta < 1.00) * (6.00 <= pt * cosh(eta) && pt * cosh(eta) < 6.40) * (0.957300) +
    (0.90 <= eta && eta < 1.00) * (6.40 <= pt * cosh(eta) && pt * cosh(eta) < 6.70) * (0.933891) +
    (0.90 <= eta && eta < 1.00) * (6.70 <= pt * cosh(eta) && pt * cosh(eta) < 7.00) * (0.910870) +
    (0.90 <= eta && eta < 1.00) * (7.00 <= pt * cosh(eta) && pt * cosh(eta) < 7.30) * (0.887137) +
    (0.90 <= eta && eta < 1.00) * (7.30 <= pt * cosh(eta) && pt * cosh(eta) < 7.60) * (0.863393) +
    (0.90 <= eta && eta < 1.00) * (7.60 <= pt * cosh(eta) && pt * cosh(eta) < 7.90) * (0.840634) +
    (0.90 <= eta && eta < 1.00) * (7.90 <= pt * cosh(eta) && pt * cosh(eta) < 8.20) * (0.818273) +
    (0.90 <= eta && eta < 1.00) * (8.20 <= pt * cosh(eta) && pt * cosh(eta) < 8.60) * (0.794422) +
    (0.90 <= eta && eta < 1.00) * (8.60 <= pt * cosh(eta) && pt * cosh(eta) < 9.00) * (0.769639) +
    (0.90 <= eta && eta < 1.00) * (9.00 <= pt * cosh(eta) && pt * cosh(eta) < 9.40) * (0.748157) +
    (0.90 <= eta && eta < 1.00) * (9.40 <= pt * cosh(eta) && pt * cosh(eta) < 10.00) * (0.724553) +
    (0.90 <= eta && eta < 1.00) * (10.00 <= pt * cosh(eta) && pt * cosh(eta) < 10.70) * (0.699984) +
    (0.90 <= eta && eta < 1.00) * (10.70 <= pt * cosh(eta) && pt * cosh(eta) < 11.60) * (0.677138) +
    (0.90 <= eta && eta < 1.00) * (11.60 <= pt * cosh(eta) && pt * cosh(eta) < 13.00) * (0.654326) +
    (0.90 <= eta && eta < 1.00) * (13.00 <= pt * cosh(eta) && pt * cosh(eta) < 16.20) * (0.629520) +
    (0.90 <= eta && eta < 1.00) * (16.20 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.603898) +
    (1.00 <= eta && eta < 1.10) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.000000) }

  add EfficiencyFormula {321} {211} {
    (eta<-1.00 || eta>=1.00 || pt * cosh(eta) < 0.10 || pt * cosh(eta) >= 50.00) * ( 0.00 ) +
    (-1.00 <= eta && eta < -0.90) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 6.10) * (0.000629) +
    (-1.00 <= eta && eta < -0.90) * (6.10 <= pt * cosh(eta) && pt * cosh(eta) < 6.80) * (0.020818) +
    (-1.00 <= eta && eta < -0.90) * (6.80 <= pt * cosh(eta) && pt * cosh(eta) < 7.30) * (0.043238) +
    (-1.00 <= eta && eta < -0.90) * (7.30 <= pt * cosh(eta) && pt * cosh(eta) < 7.70) * (0.064987) +
    (-1.00 <= eta && eta < -0.90) * (7.70 <= pt * cosh(eta) && pt * cosh(eta) < 8.10) * (0.086643) +
    (-1.00 <= eta && eta < -0.90) * (8.10 <= pt * cosh(eta) && pt * cosh(eta) < 8.50) * (0.108246) +
    (-1.00 <= eta && eta < -0.90) * (8.50 <= pt * cosh(eta) && pt * cosh(eta) < 8.90) * (0.131038) +
    (-1.00 <= eta && eta < -0.90) * (8.90 <= pt * cosh(eta) && pt * cosh(eta) < 9.30) * (0.153265) +
    (-1.00 <= eta && eta < -0.90) * (9.30 <= pt * cosh(eta) && pt * cosh(eta) < 9.70) * (0.173948) +
    (-1.00 <= eta && eta < -0.90) * (9.70 <= pt * cosh(eta) && pt * cosh(eta) < 10.20) * (0.196026) +
    (-1.00 <= eta && eta < -0.90) * (10.20 <= pt * cosh(eta) && pt * cosh(eta) < 10.70) * (0.219035) +
    (-1.00 <= eta && eta < -0.90) * (10.70 <= pt * cosh(eta) && pt * cosh(eta) < 11.30) * (0.242086) +
    (-1.00 <= eta && eta < -0.90) * (11.30 <= pt * cosh(eta) && pt * cosh(eta) < 11.90) * (0.264737) +
    (-1.00 <= eta && eta < -0.90) * (11.90 <= pt * cosh(eta) && pt * cosh(eta) < 12.60) * (0.286490) +
    (-1.00 <= eta && eta < -0.90) * (12.60 <= pt * cosh(eta) && pt * cosh(eta) < 13.40) * (0.308418) +
    (-1.00 <= eta && eta < -0.90) * (13.40 <= pt * cosh(eta) && pt * cosh(eta) < 14.30) * (0.329730) +
    (-1.00 <= eta && eta < -0.90) * (14.30 <= pt * cosh(eta) && pt * cosh(eta) < 15.40) * (0.350771) +
    (-1.00 <= eta && eta < -0.90) * (15.40 <= pt * cosh(eta) && pt * cosh(eta) < 16.80) * (0.372200) +
    (-1.00 <= eta && eta < -0.90) * (16.80 <= pt * cosh(eta) && pt * cosh(eta) < 18.60) * (0.393662) +
    (-1.00 <= eta && eta < -0.90) * (18.60 <= pt * cosh(eta) && pt * cosh(eta) < 21.20) * (0.415409) +
    (-1.00 <= eta && eta < -0.90) * (21.20 <= pt * cosh(eta) && pt * cosh(eta) < 25.30) * (0.437641) +
    (-1.00 <= eta && eta < -0.90) * (25.30 <= pt * cosh(eta) && pt * cosh(eta) < 34.20) * (0.461287) +
    (-1.00 <= eta && eta < -0.90) * (34.20 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.480405) +
    (-0.90 <= eta && eta < -0.80) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 6.00) * (0.000601) +
    (-0.90 <= eta && eta < -0.80) * (6.00 <= pt * cosh(eta) && pt * cosh(eta) < 6.70) * (0.020257) +
    (-0.90 <= eta && eta < -0.80) * (6.70 <= pt * cosh(eta) && pt * cosh(eta) < 7.20) * (0.043042) +
    (-0.90 <= eta && eta < -0.80) * (7.20 <= pt * cosh(eta) && pt * cosh(eta) < 7.60) * (0.065328) +
    (-0.90 <= eta && eta < -0.80) * (7.60 <= pt * cosh(eta) && pt * cosh(eta) < 8.00) * (0.087115) +
    (-0.90 <= eta && eta < -0.80) * (8.00 <= pt * cosh(eta) && pt * cosh(eta) < 8.40) * (0.109928) +
    (-0.90 <= eta && eta < -0.80) * (8.40 <= pt * cosh(eta) && pt * cosh(eta) < 8.80) * (0.132996) +
    (-0.90 <= eta && eta < -0.80) * (8.80 <= pt * cosh(eta) && pt * cosh(eta) < 9.30) * (0.156649) +
    (-0.90 <= eta && eta < -0.80) * (9.30 <= pt * cosh(eta) && pt * cosh(eta) < 9.70) * (0.180600) +
    (-0.90 <= eta && eta < -0.80) * (9.70 <= pt * cosh(eta) && pt * cosh(eta) < 10.20) * (0.203148) +
    (-0.90 <= eta && eta < -0.80) * (10.20 <= pt * cosh(eta) && pt * cosh(eta) < 10.70) * (0.225705) +
    (-0.90 <= eta && eta < -0.80) * (10.70 <= pt * cosh(eta) && pt * cosh(eta) < 11.30) * (0.248439) +
    (-0.90 <= eta && eta < -0.80) * (11.30 <= pt * cosh(eta) && pt * cosh(eta) < 11.90) * (0.270714) +
    (-0.90 <= eta && eta < -0.80) * (11.90 <= pt * cosh(eta) && pt * cosh(eta) < 12.60) * (0.292054) +
    (-0.90 <= eta && eta < -0.80) * (12.60 <= pt * cosh(eta) && pt * cosh(eta) < 13.40) * (0.313521) +
    (-0.90 <= eta && eta < -0.80) * (13.40 <= pt * cosh(eta) && pt * cosh(eta) < 14.40) * (0.335443) +
    (-0.90 <= eta && eta < -0.80) * (14.40 <= pt * cosh(eta) && pt * cosh(eta) < 15.60) * (0.357624) +
    (-0.90 <= eta && eta < -0.80) * (15.60 <= pt * cosh(eta) && pt * cosh(eta) < 17.10) * (0.379405) +
    (-0.90 <= eta && eta < -0.80) * (17.10 <= pt * cosh(eta) && pt * cosh(eta) < 19.10) * (0.401051) +
    (-0.90 <= eta && eta < -0.80) * (19.10 <= pt * cosh(eta) && pt * cosh(eta) < 22.00) * (0.422825) +
    (-0.90 <= eta && eta < -0.80) * (22.00 <= pt * cosh(eta) && pt * cosh(eta) < 27.00) * (0.445280) +
    (-0.90 <= eta && eta < -0.80) * (27.00 <= pt * cosh(eta) && pt * cosh(eta) < 41.00) * (0.470610) +
    (-0.90 <= eta && eta < -0.80) * (41.00 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.484124) +
    (-0.80 <= eta && eta < -0.70) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 6.00) * (0.000757) +
    (-0.80 <= eta && eta < -0.70) * (6.00 <= pt * cosh(eta) && pt * cosh(eta) < 6.60) * (0.022320) +
    (-0.80 <= eta && eta < -0.70) * (6.60 <= pt * cosh(eta) && pt * cosh(eta) < 7.10) * (0.044193) +
    (-0.80 <= eta && eta < -0.70) * (7.10 <= pt * cosh(eta) && pt * cosh(eta) < 7.50) * (0.067162) +
    (-0.80 <= eta && eta < -0.70) * (7.50 <= pt * cosh(eta) && pt * cosh(eta) < 7.90) * (0.088655) +
    (-0.80 <= eta && eta < -0.70) * (7.90 <= pt * cosh(eta) && pt * cosh(eta) < 8.30) * (0.111401) +
    (-0.80 <= eta && eta < -0.70) * (8.30 <= pt * cosh(eta) && pt * cosh(eta) < 8.70) * (0.134165) +
    (-0.80 <= eta && eta < -0.70) * (8.70 <= pt * cosh(eta) && pt * cosh(eta) < 9.10) * (0.156920) +
    (-0.80 <= eta && eta < -0.70) * (9.10 <= pt * cosh(eta) && pt * cosh(eta) < 9.50) * (0.177718) +
    (-0.80 <= eta && eta < -0.70) * (9.50 <= pt * cosh(eta) && pt * cosh(eta) < 10.00) * (0.199864) +
    (-0.80 <= eta && eta < -0.70) * (10.00 <= pt * cosh(eta) && pt * cosh(eta) < 10.50) * (0.223340) +
    (-0.80 <= eta && eta < -0.70) * (10.50 <= pt * cosh(eta) && pt * cosh(eta) < 11.10) * (0.246613) +
    (-0.80 <= eta && eta < -0.70) * (11.10 <= pt * cosh(eta) && pt * cosh(eta) < 11.70) * (0.269390) +
    (-0.80 <= eta && eta < -0.70) * (11.70 <= pt * cosh(eta) && pt * cosh(eta) < 12.40) * (0.291176) +
    (-0.80 <= eta && eta < -0.70) * (12.40 <= pt * cosh(eta) && pt * cosh(eta) < 13.20) * (0.313049) +
    (-0.80 <= eta && eta < -0.70) * (13.20 <= pt * cosh(eta) && pt * cosh(eta) < 14.20) * (0.335333) +
    (-0.80 <= eta && eta < -0.70) * (14.20 <= pt * cosh(eta) && pt * cosh(eta) < 15.40) * (0.357820) +
    (-0.80 <= eta && eta < -0.70) * (15.40 <= pt * cosh(eta) && pt * cosh(eta) < 16.90) * (0.379834) +
    (-0.80 <= eta && eta < -0.70) * (16.90 <= pt * cosh(eta) && pt * cosh(eta) < 18.90) * (0.401637) +
    (-0.80 <= eta && eta < -0.70) * (18.90 <= pt * cosh(eta) && pt * cosh(eta) < 21.80) * (0.423486) +
    (-0.80 <= eta && eta < -0.70) * (21.80 <= pt * cosh(eta) && pt * cosh(eta) < 26.80) * (0.445919) +
    (-0.80 <= eta && eta < -0.70) * (26.80 <= pt * cosh(eta) && pt * cosh(eta) < 41.30) * (0.471428) +
    (-0.80 <= eta && eta < -0.70) * (41.30 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.484680) +
    (-0.70 <= eta && eta < -0.60) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 5.90) * (0.000627) +
    (-0.70 <= eta && eta < -0.60) * (5.90 <= pt * cosh(eta) && pt * cosh(eta) < 6.50) * (0.019743) +
    (-0.70 <= eta && eta < -0.60) * (6.50 <= pt * cosh(eta) && pt * cosh(eta) < 7.00) * (0.041036) +
    (-0.70 <= eta && eta < -0.60) * (7.00 <= pt * cosh(eta) && pt * cosh(eta) < 7.40) * (0.063090) +
    (-0.70 <= eta && eta < -0.60) * (7.40 <= pt * cosh(eta) && pt * cosh(eta) < 7.80) * (0.084632) +
    (-0.70 <= eta && eta < -0.60) * (7.80 <= pt * cosh(eta) && pt * cosh(eta) < 8.20) * (0.107944) +
    (-0.70 <= eta && eta < -0.60) * (8.20 <= pt * cosh(eta) && pt * cosh(eta) < 8.60) * (0.131352) +
    (-0.70 <= eta && eta < -0.60) * (8.60 <= pt * cosh(eta) && pt * cosh(eta) < 9.00) * (0.153268) +
    (-0.70 <= eta && eta < -0.60) * (9.00 <= pt * cosh(eta) && pt * cosh(eta) < 9.40) * (0.174594) +
    (-0.70 <= eta && eta < -0.60) * (9.40 <= pt * cosh(eta) && pt * cosh(eta) < 9.90) * (0.197491) +
    (-0.70 <= eta && eta < -0.60) * (9.90 <= pt * cosh(eta) && pt * cosh(eta) < 10.40) * (0.221492) +
    (-0.70 <= eta && eta < -0.60) * (10.40 <= pt * cosh(eta) && pt * cosh(eta) < 10.90) * (0.243056) +
    (-0.70 <= eta && eta < -0.60) * (10.90 <= pt * cosh(eta) && pt * cosh(eta) < 11.50) * (0.264487) +
    (-0.70 <= eta && eta < -0.60) * (11.50 <= pt * cosh(eta) && pt * cosh(eta) < 12.20) * (0.286988) +
    (-0.70 <= eta && eta < -0.60) * (12.20 <= pt * cosh(eta) && pt * cosh(eta) < 13.00) * (0.309562) +
    (-0.70 <= eta && eta < -0.60) * (13.00 <= pt * cosh(eta) && pt * cosh(eta) < 13.90) * (0.331387) +
    (-0.70 <= eta && eta < -0.60) * (13.90 <= pt * cosh(eta) && pt * cosh(eta) < 15.00) * (0.352810) +
    (-0.70 <= eta && eta < -0.60) * (15.00 <= pt * cosh(eta) && pt * cosh(eta) < 16.40) * (0.374492) +
    (-0.70 <= eta && eta < -0.60) * (16.40 <= pt * cosh(eta) && pt * cosh(eta) < 18.20) * (0.396058) +
    (-0.70 <= eta && eta < -0.60) * (18.20 <= pt * cosh(eta) && pt * cosh(eta) < 20.80) * (0.417745) +
    (-0.70 <= eta && eta < -0.60) * (20.80 <= pt * cosh(eta) && pt * cosh(eta) < 25.00) * (0.439965) +
    (-0.70 <= eta && eta < -0.60) * (25.00 <= pt * cosh(eta) && pt * cosh(eta) < 34.60) * (0.463862) +
    (-0.70 <= eta && eta < -0.60) * (34.60 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.481927) +
    (-0.60 <= eta && eta < -0.50) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 5.90) * (0.000606) +
    (-0.60 <= eta && eta < -0.50) * (5.90 <= pt * cosh(eta) && pt * cosh(eta) < 6.50) * (0.019510) +
    (-0.60 <= eta && eta < -0.50) * (6.50 <= pt * cosh(eta) && pt * cosh(eta) < 7.00) * (0.040196) +
    (-0.60 <= eta && eta < -0.50) * (7.00 <= pt * cosh(eta) && pt * cosh(eta) < 7.40) * (0.061442) +
    (-0.60 <= eta && eta < -0.50) * (7.40 <= pt * cosh(eta) && pt * cosh(eta) < 7.80) * (0.083297) +
    (-0.60 <= eta && eta < -0.50) * (7.80 <= pt * cosh(eta) && pt * cosh(eta) < 8.20) * (0.105846) +
    (-0.60 <= eta && eta < -0.50) * (8.20 <= pt * cosh(eta) && pt * cosh(eta) < 8.60) * (0.129010) +
    (-0.60 <= eta && eta < -0.50) * (8.60 <= pt * cosh(eta) && pt * cosh(eta) < 9.00) * (0.151227) +
    (-0.60 <= eta && eta < -0.50) * (9.00 <= pt * cosh(eta) && pt * cosh(eta) < 9.40) * (0.173018) +
    (-0.60 <= eta && eta < -0.50) * (9.40 <= pt * cosh(eta) && pt * cosh(eta) < 9.90) * (0.197095) +
    (-0.60 <= eta && eta < -0.50) * (9.90 <= pt * cosh(eta) && pt * cosh(eta) < 10.40) * (0.220621) +
    (-0.60 <= eta && eta < -0.50) * (10.40 <= pt * cosh(eta) && pt * cosh(eta) < 10.90) * (0.242223) +
    (-0.60 <= eta && eta < -0.50) * (10.90 <= pt * cosh(eta) && pt * cosh(eta) < 11.50) * (0.263700) +
    (-0.60 <= eta && eta < -0.50) * (11.50 <= pt * cosh(eta) && pt * cosh(eta) < 12.20) * (0.286256) +
    (-0.60 <= eta && eta < -0.50) * (12.20 <= pt * cosh(eta) && pt * cosh(eta) < 13.00) * (0.308892) +
    (-0.60 <= eta && eta < -0.50) * (13.00 <= pt * cosh(eta) && pt * cosh(eta) < 13.90) * (0.330782) +
    (-0.60 <= eta && eta < -0.50) * (13.90 <= pt * cosh(eta) && pt * cosh(eta) < 15.00) * (0.352274) +
    (-0.60 <= eta && eta < -0.50) * (15.00 <= pt * cosh(eta) && pt * cosh(eta) < 16.40) * (0.374028) +
    (-0.60 <= eta && eta < -0.50) * (16.40 <= pt * cosh(eta) && pt * cosh(eta) < 18.20) * (0.395670) +
    (-0.60 <= eta && eta < -0.50) * (18.20 <= pt * cosh(eta) && pt * cosh(eta) < 20.80) * (0.417435) +
    (-0.60 <= eta && eta < -0.50) * (20.80 <= pt * cosh(eta) && pt * cosh(eta) < 25.00) * (0.439737) +
    (-0.60 <= eta && eta < -0.50) * (25.00 <= pt * cosh(eta) && pt * cosh(eta) < 34.60) * (0.463725) +
    (-0.60 <= eta && eta < -0.50) * (34.60 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.481858) +
    (-0.50 <= eta && eta < -0.40) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 5.90) * (0.000754) +
    (-0.50 <= eta && eta < -0.40) * (5.90 <= pt * cosh(eta) && pt * cosh(eta) < 6.50) * (0.022057) +
    (-0.50 <= eta && eta < -0.40) * (6.50 <= pt * cosh(eta) && pt * cosh(eta) < 7.00) * (0.044626) +
    (-0.50 <= eta && eta < -0.40) * (7.00 <= pt * cosh(eta) && pt * cosh(eta) < 7.40) * (0.067546) +
    (-0.50 <= eta && eta < -0.40) * (7.40 <= pt * cosh(eta) && pt * cosh(eta) < 7.80) * (0.089482) +
    (-0.50 <= eta && eta < -0.40) * (7.80 <= pt * cosh(eta) && pt * cosh(eta) < 8.20) * (0.112586) +
    (-0.50 <= eta && eta < -0.40) * (8.20 <= pt * cosh(eta) && pt * cosh(eta) < 8.60) * (0.135083) +
    (-0.50 <= eta && eta < -0.40) * (8.60 <= pt * cosh(eta) && pt * cosh(eta) < 9.00) * (0.158756) +
    (-0.50 <= eta && eta < -0.40) * (9.00 <= pt * cosh(eta) && pt * cosh(eta) < 9.40) * (0.180504) +
    (-0.50 <= eta && eta < -0.40) * (9.40 <= pt * cosh(eta) && pt * cosh(eta) < 9.90) * (0.204102) +
    (-0.50 <= eta && eta < -0.40) * (9.90 <= pt * cosh(eta) && pt * cosh(eta) < 10.40) * (0.227320) +
    (-0.50 <= eta && eta < -0.40) * (10.40 <= pt * cosh(eta) && pt * cosh(eta) < 11.00) * (0.250611) +
    (-0.50 <= eta && eta < -0.40) * (11.00 <= pt * cosh(eta) && pt * cosh(eta) < 11.60) * (0.273341) +
    (-0.50 <= eta && eta < -0.40) * (11.60 <= pt * cosh(eta) && pt * cosh(eta) < 12.30) * (0.295022) +
    (-0.50 <= eta && eta < -0.40) * (12.30 <= pt * cosh(eta) && pt * cosh(eta) < 13.10) * (0.316735) +
    (-0.50 <= eta && eta < -0.40) * (13.10 <= pt * cosh(eta) && pt * cosh(eta) < 14.10) * (0.338798) +
    (-0.50 <= eta && eta < -0.40) * (14.10 <= pt * cosh(eta) && pt * cosh(eta) < 15.30) * (0.361005) +
    (-0.50 <= eta && eta < -0.40) * (15.30 <= pt * cosh(eta) && pt * cosh(eta) < 16.80) * (0.382690) +
    (-0.50 <= eta && eta < -0.40) * (16.80 <= pt * cosh(eta) && pt * cosh(eta) < 18.80) * (0.404114) +
    (-0.50 <= eta && eta < -0.40) * (18.80 <= pt * cosh(eta) && pt * cosh(eta) < 21.80) * (0.425866) +
    (-0.50 <= eta && eta < -0.40) * (21.80 <= pt * cosh(eta) && pt * cosh(eta) < 27.10) * (0.448463) +
    (-0.50 <= eta && eta < -0.40) * (27.10 <= pt * cosh(eta) && pt * cosh(eta) < 45.70) * (0.475401) +
    (-0.50 <= eta && eta < -0.40) * (45.70 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.486661) +
    (-0.40 <= eta && eta < -0.30) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 5.90) * (0.000728) +
    (-0.40 <= eta && eta < -0.30) * (5.90 <= pt * cosh(eta) && pt * cosh(eta) < 6.50) * (0.021373) +
    (-0.40 <= eta && eta < -0.30) * (6.50 <= pt * cosh(eta) && pt * cosh(eta) < 7.00) * (0.042787) +
    (-0.40 <= eta && eta < -0.30) * (7.00 <= pt * cosh(eta) && pt * cosh(eta) < 7.40) * (0.065124) +
    (-0.40 <= eta && eta < -0.30) * (7.40 <= pt * cosh(eta) && pt * cosh(eta) < 7.80) * (0.087923) +
    (-0.40 <= eta && eta < -0.30) * (7.80 <= pt * cosh(eta) && pt * cosh(eta) < 8.20) * (0.111671) +
    (-0.40 <= eta && eta < -0.30) * (8.20 <= pt * cosh(eta) && pt * cosh(eta) < 8.60) * (0.134974) +
    (-0.40 <= eta && eta < -0.30) * (8.60 <= pt * cosh(eta) && pt * cosh(eta) < 9.00) * (0.157718) +
    (-0.40 <= eta && eta < -0.30) * (9.00 <= pt * cosh(eta) && pt * cosh(eta) < 9.40) * (0.179824) +
    (-0.40 <= eta && eta < -0.30) * (9.40 <= pt * cosh(eta) && pt * cosh(eta) < 9.90) * (0.202629) +
    (-0.40 <= eta && eta < -0.30) * (9.90 <= pt * cosh(eta) && pt * cosh(eta) < 10.40) * (0.225577) +
    (-0.40 <= eta && eta < -0.30) * (10.40 <= pt * cosh(eta) && pt * cosh(eta) < 11.00) * (0.248957) +
    (-0.40 <= eta && eta < -0.30) * (11.00 <= pt * cosh(eta) && pt * cosh(eta) < 11.60) * (0.271791) +
    (-0.40 <= eta && eta < -0.30) * (11.60 <= pt * cosh(eta) && pt * cosh(eta) < 12.30) * (0.293585) +
    (-0.40 <= eta && eta < -0.30) * (12.30 <= pt * cosh(eta) && pt * cosh(eta) < 13.10) * (0.315422) +
    (-0.40 <= eta && eta < -0.30) * (13.10 <= pt * cosh(eta) && pt * cosh(eta) < 14.10) * (0.337621) +
    (-0.40 <= eta && eta < -0.30) * (14.10 <= pt * cosh(eta) && pt * cosh(eta) < 15.30) * (0.359975) +
    (-0.40 <= eta && eta < -0.30) * (15.30 <= pt * cosh(eta) && pt * cosh(eta) < 16.80) * (0.381809) +
    (-0.40 <= eta && eta < -0.30) * (16.80 <= pt * cosh(eta) && pt * cosh(eta) < 18.80) * (0.403386) +
    (-0.40 <= eta && eta < -0.30) * (18.80 <= pt * cosh(eta) && pt * cosh(eta) < 21.80) * (0.425299) +
    (-0.40 <= eta && eta < -0.30) * (21.80 <= pt * cosh(eta) && pt * cosh(eta) < 27.10) * (0.448066) +
    (-0.40 <= eta && eta < -0.30) * (27.10 <= pt * cosh(eta) && pt * cosh(eta) < 45.00) * (0.474824) +
    (-0.40 <= eta && eta < -0.30) * (45.00 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.486349) +
    (-0.30 <= eta && eta < -0.20) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 5.90) * (0.000713) +
    (-0.30 <= eta && eta < -0.20) * (5.90 <= pt * cosh(eta) && pt * cosh(eta) < 6.50) * (0.021188) +
    (-0.30 <= eta && eta < -0.20) * (6.50 <= pt * cosh(eta) && pt * cosh(eta) < 7.00) * (0.043291) +
    (-0.30 <= eta && eta < -0.20) * (7.00 <= pt * cosh(eta) && pt * cosh(eta) < 7.40) * (0.065636) +
    (-0.30 <= eta && eta < -0.20) * (7.40 <= pt * cosh(eta) && pt * cosh(eta) < 7.80) * (0.087346) +
    (-0.30 <= eta && eta < -0.20) * (7.80 <= pt * cosh(eta) && pt * cosh(eta) < 8.20) * (0.111189) +
    (-0.30 <= eta && eta < -0.20) * (8.20 <= pt * cosh(eta) && pt * cosh(eta) < 8.60) * (0.135124) +
    (-0.30 <= eta && eta < -0.20) * (8.60 <= pt * cosh(eta) && pt * cosh(eta) < 9.00) * (0.157289) +
    (-0.30 <= eta && eta < -0.20) * (9.00 <= pt * cosh(eta) && pt * cosh(eta) < 9.40) * (0.178343) +
    (-0.30 <= eta && eta < -0.20) * (9.40 <= pt * cosh(eta) && pt * cosh(eta) < 9.90) * (0.201490) +
    (-0.30 <= eta && eta < -0.20) * (9.90 <= pt * cosh(eta) && pt * cosh(eta) < 10.40) * (0.225445) +
    (-0.30 <= eta && eta < -0.20) * (10.40 <= pt * cosh(eta) && pt * cosh(eta) < 11.00) * (0.248831) +
    (-0.30 <= eta && eta < -0.20) * (11.00 <= pt * cosh(eta) && pt * cosh(eta) < 11.60) * (0.271673) +
    (-0.30 <= eta && eta < -0.20) * (11.60 <= pt * cosh(eta) && pt * cosh(eta) < 12.30) * (0.293475) +
    (-0.30 <= eta && eta < -0.20) * (12.30 <= pt * cosh(eta) && pt * cosh(eta) < 13.10) * (0.315321) +
    (-0.30 <= eta && eta < -0.20) * (13.10 <= pt * cosh(eta) && pt * cosh(eta) < 14.10) * (0.337531) +
    (-0.30 <= eta && eta < -0.20) * (14.10 <= pt * cosh(eta) && pt * cosh(eta) < 15.30) * (0.359896) +
    (-0.30 <= eta && eta < -0.20) * (15.30 <= pt * cosh(eta) && pt * cosh(eta) < 16.80) * (0.381742) +
    (-0.30 <= eta && eta < -0.20) * (16.80 <= pt * cosh(eta) && pt * cosh(eta) < 18.80) * (0.403330) +
    (-0.30 <= eta && eta < -0.20) * (18.80 <= pt * cosh(eta) && pt * cosh(eta) < 21.80) * (0.425256) +
    (-0.30 <= eta && eta < -0.20) * (21.80 <= pt * cosh(eta) && pt * cosh(eta) < 27.10) * (0.448036) +
    (-0.30 <= eta && eta < -0.20) * (27.10 <= pt * cosh(eta) && pt * cosh(eta) < 45.00) * (0.474810) +
    (-0.30 <= eta && eta < -0.20) * (45.00 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.486341) +
    (-0.20 <= eta && eta < -0.10) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 5.90) * (0.000670) +
    (-0.20 <= eta && eta < -0.10) * (5.90 <= pt * cosh(eta) && pt * cosh(eta) < 6.50) * (0.020397) +
    (-0.20 <= eta && eta < -0.10) * (6.50 <= pt * cosh(eta) && pt * cosh(eta) < 7.00) * (0.041717) +
    (-0.20 <= eta && eta < -0.10) * (7.00 <= pt * cosh(eta) && pt * cosh(eta) < 7.40) * (0.064356) +
    (-0.20 <= eta && eta < -0.10) * (7.40 <= pt * cosh(eta) && pt * cosh(eta) < 7.80) * (0.085635) +
    (-0.20 <= eta && eta < -0.10) * (7.80 <= pt * cosh(eta) && pt * cosh(eta) < 8.20) * (0.109110) +
    (-0.20 <= eta && eta < -0.10) * (8.20 <= pt * cosh(eta) && pt * cosh(eta) < 8.60) * (0.131602) +
    (-0.20 <= eta && eta < -0.10) * (8.60 <= pt * cosh(eta) && pt * cosh(eta) < 9.00) * (0.153486) +
    (-0.20 <= eta && eta < -0.10) * (9.00 <= pt * cosh(eta) && pt * cosh(eta) < 9.40) * (0.175586) +
    (-0.20 <= eta && eta < -0.10) * (9.40 <= pt * cosh(eta) && pt * cosh(eta) < 9.90) * (0.198903) +
    (-0.20 <= eta && eta < -0.10) * (9.90 <= pt * cosh(eta) && pt * cosh(eta) < 10.40) * (0.222061) +
    (-0.20 <= eta && eta < -0.10) * (10.40 <= pt * cosh(eta) && pt * cosh(eta) < 10.90) * (0.243601) +
    (-0.20 <= eta && eta < -0.10) * (10.90 <= pt * cosh(eta) && pt * cosh(eta) < 11.50) * (0.265002) +
    (-0.20 <= eta && eta < -0.10) * (11.50 <= pt * cosh(eta) && pt * cosh(eta) < 12.20) * (0.287466) +
    (-0.20 <= eta && eta < -0.10) * (12.20 <= pt * cosh(eta) && pt * cosh(eta) < 13.00) * (0.310000) +
    (-0.20 <= eta && eta < -0.10) * (13.00 <= pt * cosh(eta) && pt * cosh(eta) < 13.90) * (0.331781) +
    (-0.20 <= eta && eta < -0.10) * (13.90 <= pt * cosh(eta) && pt * cosh(eta) < 15.00) * (0.353161) +
    (-0.20 <= eta && eta < -0.10) * (15.00 <= pt * cosh(eta) && pt * cosh(eta) < 16.40) * (0.374795) +
    (-0.20 <= eta && eta < -0.10) * (16.40 <= pt * cosh(eta) && pt * cosh(eta) < 18.20) * (0.396312) +
    (-0.20 <= eta && eta < -0.10) * (18.20 <= pt * cosh(eta) && pt * cosh(eta) < 20.80) * (0.417947) +
    (-0.20 <= eta && eta < -0.10) * (20.80 <= pt * cosh(eta) && pt * cosh(eta) < 25.10) * (0.440351) +
    (-0.20 <= eta && eta < -0.10) * (25.10 <= pt * cosh(eta) && pt * cosh(eta) < 35.00) * (0.464506) +
    (-0.20 <= eta && eta < -0.10) * (35.00 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.482178) +
    (-0.10 <= eta && eta < -0.00) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 6.10) * (0.000742) +
    (-0.10 <= eta && eta < -0.00) * (6.10 <= pt * cosh(eta) && pt * cosh(eta) < 6.70) * (0.021676) +
    (-0.10 <= eta && eta < -0.00) * (6.70 <= pt * cosh(eta) && pt * cosh(eta) < 7.20) * (0.043000) +
    (-0.10 <= eta && eta < -0.00) * (7.20 <= pt * cosh(eta) && pt * cosh(eta) < 7.60) * (0.064909) +
    (-0.10 <= eta && eta < -0.00) * (7.60 <= pt * cosh(eta) && pt * cosh(eta) < 8.00) * (0.086553) +
    (-0.10 <= eta && eta < -0.00) * (8.00 <= pt * cosh(eta) && pt * cosh(eta) < 8.40) * (0.108414) +
    (-0.10 <= eta && eta < -0.00) * (8.40 <= pt * cosh(eta) && pt * cosh(eta) < 8.90) * (0.131561) +
    (-0.10 <= eta && eta < -0.00) * (8.90 <= pt * cosh(eta) && pt * cosh(eta) < 9.30) * (0.155629) +
    (-0.10 <= eta && eta < -0.00) * (9.30 <= pt * cosh(eta) && pt * cosh(eta) < 9.70) * (0.177181) +
    (-0.10 <= eta && eta < -0.00) * (9.70 <= pt * cosh(eta) && pt * cosh(eta) < 10.20) * (0.200013) +
    (-0.10 <= eta && eta < -0.00) * (10.20 <= pt * cosh(eta) && pt * cosh(eta) < 10.70) * (0.223039) +
    (-0.10 <= eta && eta < -0.00) * (10.70 <= pt * cosh(eta) && pt * cosh(eta) < 11.30) * (0.245902) +
    (-0.10 <= eta && eta < -0.00) * (11.30 <= pt * cosh(eta) && pt * cosh(eta) < 11.90) * (0.268329) +
    (-0.10 <= eta && eta < -0.00) * (11.90 <= pt * cosh(eta) && pt * cosh(eta) < 12.60) * (0.289835) +
    (-0.10 <= eta && eta < -0.00) * (12.60 <= pt * cosh(eta) && pt * cosh(eta) < 13.40) * (0.311487) +
    (-0.10 <= eta && eta < -0.00) * (13.40 <= pt * cosh(eta) && pt * cosh(eta) < 14.40) * (0.333613) +
    (-0.10 <= eta && eta < -0.00) * (14.40 <= pt * cosh(eta) && pt * cosh(eta) < 15.60) * (0.356016) +
    (-0.10 <= eta && eta < -0.00) * (15.60 <= pt * cosh(eta) && pt * cosh(eta) < 17.10) * (0.378024) +
    (-0.10 <= eta && eta < -0.00) * (17.10 <= pt * cosh(eta) && pt * cosh(eta) < 19.10) * (0.399906) +
    (-0.10 <= eta && eta < -0.00) * (19.10 <= pt * cosh(eta) && pt * cosh(eta) < 22.00) * (0.421924) +
    (-0.10 <= eta && eta < -0.00) * (22.00 <= pt * cosh(eta) && pt * cosh(eta) < 26.90) * (0.444432) +
    (-0.10 <= eta && eta < -0.00) * (26.90 <= pt * cosh(eta) && pt * cosh(eta) < 40.00) * (0.469406) +
    (-0.10 <= eta && eta < -0.00) * (40.00 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.483535) +
    (-0.00 <= eta && eta < 0.10) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 6.00) * (0.000700) +
    (-0.00 <= eta && eta < 0.10) * (6.00 <= pt * cosh(eta) && pt * cosh(eta) < 6.60) * (0.021257) +
    (-0.00 <= eta && eta < 0.10) * (6.60 <= pt * cosh(eta) && pt * cosh(eta) < 7.10) * (0.043040) +
    (-0.00 <= eta && eta < 0.10) * (7.10 <= pt * cosh(eta) && pt * cosh(eta) < 7.50) * (0.064435) +
    (-0.00 <= eta && eta < 0.10) * (7.50 <= pt * cosh(eta) && pt * cosh(eta) < 7.90) * (0.086042) +
    (-0.00 <= eta && eta < 0.10) * (7.90 <= pt * cosh(eta) && pt * cosh(eta) < 8.30) * (0.109607) +
    (-0.00 <= eta && eta < 0.10) * (8.30 <= pt * cosh(eta) && pt * cosh(eta) < 8.70) * (0.132930) +
    (-0.00 <= eta && eta < 0.10) * (8.70 <= pt * cosh(eta) && pt * cosh(eta) < 9.10) * (0.153931) +
    (-0.00 <= eta && eta < 0.10) * (9.10 <= pt * cosh(eta) && pt * cosh(eta) < 9.50) * (0.175705) +
    (-0.00 <= eta && eta < 0.10) * (9.50 <= pt * cosh(eta) && pt * cosh(eta) < 10.00) * (0.197877) +
    (-0.00 <= eta && eta < 0.10) * (10.00 <= pt * cosh(eta) && pt * cosh(eta) < 10.50) * (0.220292) +
    (-0.00 <= eta && eta < 0.10) * (10.50 <= pt * cosh(eta) && pt * cosh(eta) < 11.00) * (0.241706) +
    (-0.00 <= eta && eta < 0.10) * (11.00 <= pt * cosh(eta) && pt * cosh(eta) < 11.60) * (0.263021) +
    (-0.00 <= eta && eta < 0.10) * (11.60 <= pt * cosh(eta) && pt * cosh(eta) < 12.30) * (0.285436) +
    (-0.00 <= eta && eta < 0.10) * (12.30 <= pt * cosh(eta) && pt * cosh(eta) < 13.10) * (0.307965) +
    (-0.00 <= eta && eta < 0.10) * (13.10 <= pt * cosh(eta) && pt * cosh(eta) < 14.00) * (0.329785) +
    (-0.00 <= eta && eta < 0.10) * (14.00 <= pt * cosh(eta) && pt * cosh(eta) < 15.10) * (0.351244) +
    (-0.00 <= eta && eta < 0.10) * (15.10 <= pt * cosh(eta) && pt * cosh(eta) < 16.50) * (0.373004) +
    (-0.00 <= eta && eta < 0.10) * (16.50 <= pt * cosh(eta) && pt * cosh(eta) < 18.30) * (0.394692) +
    (-0.00 <= eta && eta < 0.10) * (18.30 <= pt * cosh(eta) && pt * cosh(eta) < 20.90) * (0.416547) +
    (-0.00 <= eta && eta < 0.10) * (20.90 <= pt * cosh(eta) && pt * cosh(eta) < 25.10) * (0.438994) +
    (-0.00 <= eta && eta < 0.10) * (25.10 <= pt * cosh(eta) && pt * cosh(eta) < 34.40) * (0.462884) +
    (-0.00 <= eta && eta < 0.10) * (34.40 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.481363) +
    (0.10 <= eta && eta < 0.20) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 5.80) * (0.000617) +
    (0.10 <= eta && eta < 0.20) * (5.80 <= pt * cosh(eta) && pt * cosh(eta) < 6.40) * (0.019051) +
    (0.10 <= eta && eta < 0.20) * (6.40 <= pt * cosh(eta) && pt * cosh(eta) < 6.90) * (0.040330) +
    (0.10 <= eta && eta < 0.20) * (6.90 <= pt * cosh(eta) && pt * cosh(eta) < 7.30) * (0.062862) +
    (0.10 <= eta && eta < 0.20) * (7.30 <= pt * cosh(eta) && pt * cosh(eta) < 7.70) * (0.085532) +
    (0.10 <= eta && eta < 0.20) * (7.70 <= pt * cosh(eta) && pt * cosh(eta) < 8.10) * (0.107916) +
    (0.10 <= eta && eta < 0.20) * (8.10 <= pt * cosh(eta) && pt * cosh(eta) < 8.50) * (0.131693) +
    (0.10 <= eta && eta < 0.20) * (8.50 <= pt * cosh(eta) && pt * cosh(eta) < 8.90) * (0.155745) +
    (0.10 <= eta && eta < 0.20) * (8.90 <= pt * cosh(eta) && pt * cosh(eta) < 9.30) * (0.177142) +
    (0.10 <= eta && eta < 0.20) * (9.30 <= pt * cosh(eta) && pt * cosh(eta) < 9.70) * (0.196263) +
    (0.10 <= eta && eta < 0.20) * (9.70 <= pt * cosh(eta) && pt * cosh(eta) < 10.20) * (0.218643) +
    (0.10 <= eta && eta < 0.20) * (10.20 <= pt * cosh(eta) && pt * cosh(eta) < 10.70) * (0.240921) +
    (0.10 <= eta && eta < 0.20) * (10.70 <= pt * cosh(eta) && pt * cosh(eta) < 11.30) * (0.262860) +
    (0.10 <= eta && eta < 0.20) * (11.30 <= pt * cosh(eta) && pt * cosh(eta) < 12.00) * (0.285861) +
    (0.10 <= eta && eta < 0.20) * (12.00 <= pt * cosh(eta) && pt * cosh(eta) < 12.80) * (0.308892) +
    (0.10 <= eta && eta < 0.20) * (12.80 <= pt * cosh(eta) && pt * cosh(eta) < 13.70) * (0.331107) +
    (0.10 <= eta && eta < 0.20) * (13.70 <= pt * cosh(eta) && pt * cosh(eta) < 14.80) * (0.352856) +
    (0.10 <= eta && eta < 0.20) * (14.80 <= pt * cosh(eta) && pt * cosh(eta) < 16.20) * (0.374803) +
    (0.10 <= eta && eta < 0.20) * (16.20 <= pt * cosh(eta) && pt * cosh(eta) < 18.00) * (0.396560) +
    (0.10 <= eta && eta < 0.20) * (18.00 <= pt * cosh(eta) && pt * cosh(eta) < 20.60) * (0.418357) +
    (0.10 <= eta && eta < 0.20) * (20.60 <= pt * cosh(eta) && pt * cosh(eta) < 24.90) * (0.440831) +
    (0.10 <= eta && eta < 0.20) * (24.90 <= pt * cosh(eta) && pt * cosh(eta) < 34.90) * (0.465032) +
    (0.10 <= eta && eta < 0.20) * (34.90 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.482582) +
    (0.20 <= eta && eta < 0.30) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 5.80) * (0.000736) +
    (0.20 <= eta && eta < 0.30) * (5.80 <= pt * cosh(eta) && pt * cosh(eta) < 6.40) * (0.021548) +
    (0.20 <= eta && eta < 0.30) * (6.40 <= pt * cosh(eta) && pt * cosh(eta) < 6.90) * (0.043947) +
    (0.20 <= eta && eta < 0.30) * (6.90 <= pt * cosh(eta) && pt * cosh(eta) < 7.30) * (0.067742) +
    (0.20 <= eta && eta < 0.30) * (7.30 <= pt * cosh(eta) && pt * cosh(eta) < 7.70) * (0.090454) +
    (0.20 <= eta && eta < 0.30) * (7.70 <= pt * cosh(eta) && pt * cosh(eta) < 8.10) * (0.113386) +
    (0.20 <= eta && eta < 0.30) * (8.10 <= pt * cosh(eta) && pt * cosh(eta) < 8.50) * (0.136843) +
    (0.20 <= eta && eta < 0.30) * (8.50 <= pt * cosh(eta) && pt * cosh(eta) < 8.90) * (0.160311) +
    (0.20 <= eta && eta < 0.30) * (8.90 <= pt * cosh(eta) && pt * cosh(eta) < 9.40) * (0.185362) +
    (0.20 <= eta && eta < 0.30) * (9.40 <= pt * cosh(eta) && pt * cosh(eta) < 9.80) * (0.207279) +
    (0.20 <= eta && eta < 0.30) * (9.80 <= pt * cosh(eta) && pt * cosh(eta) < 10.30) * (0.230256) +
    (0.20 <= eta && eta < 0.30) * (10.30 <= pt * cosh(eta) && pt * cosh(eta) < 10.90) * (0.253794) +
    (0.20 <= eta && eta < 0.30) * (10.90 <= pt * cosh(eta) && pt * cosh(eta) < 11.50) * (0.276522) +
    (0.20 <= eta && eta < 0.30) * (11.50 <= pt * cosh(eta) && pt * cosh(eta) < 12.20) * (0.298150) +
    (0.20 <= eta && eta < 0.30) * (12.20 <= pt * cosh(eta) && pt * cosh(eta) < 13.00) * (0.319758) +
    (0.20 <= eta && eta < 0.30) * (13.00 <= pt * cosh(eta) && pt * cosh(eta) < 14.00) * (0.341663) +
    (0.20 <= eta && eta < 0.30) * (14.00 <= pt * cosh(eta) && pt * cosh(eta) < 15.20) * (0.363659) +
    (0.20 <= eta && eta < 0.30) * (15.20 <= pt * cosh(eta) && pt * cosh(eta) < 16.70) * (0.385086) +
    (0.20 <= eta && eta < 0.30) * (16.70 <= pt * cosh(eta) && pt * cosh(eta) < 18.80) * (0.406698) +
    (0.20 <= eta && eta < 0.30) * (18.80 <= pt * cosh(eta) && pt * cosh(eta) < 21.90) * (0.428629) +
    (0.20 <= eta && eta < 0.30) * (21.90 <= pt * cosh(eta) && pt * cosh(eta) < 27.60) * (0.451291) +
    (0.20 <= eta && eta < 0.30) * (27.60 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.478658) +
    (0.30 <= eta && eta < 0.40) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 5.80) * (0.000753) +
    (0.30 <= eta && eta < 0.40) * (5.80 <= pt * cosh(eta) && pt * cosh(eta) < 6.40) * (0.021899) +
    (0.30 <= eta && eta < 0.40) * (6.40 <= pt * cosh(eta) && pt * cosh(eta) < 6.90) * (0.045558) +
    (0.30 <= eta && eta < 0.40) * (6.90 <= pt * cosh(eta) && pt * cosh(eta) < 7.30) * (0.068436) +
    (0.30 <= eta && eta < 0.40) * (7.30 <= pt * cosh(eta) && pt * cosh(eta) < 7.70) * (0.091463) +
    (0.30 <= eta && eta < 0.40) * (7.70 <= pt * cosh(eta) && pt * cosh(eta) < 8.10) * (0.115435) +
    (0.30 <= eta && eta < 0.40) * (8.10 <= pt * cosh(eta) && pt * cosh(eta) < 8.50) * (0.139105) +
    (0.30 <= eta && eta < 0.40) * (8.50 <= pt * cosh(eta) && pt * cosh(eta) < 8.90) * (0.162598) +
    (0.30 <= eta && eta < 0.40) * (8.90 <= pt * cosh(eta) && pt * cosh(eta) < 9.30) * (0.182795) +
    (0.30 <= eta && eta < 0.40) * (9.30 <= pt * cosh(eta) && pt * cosh(eta) < 9.80) * (0.205458) +
    (0.30 <= eta && eta < 0.40) * (9.80 <= pt * cosh(eta) && pt * cosh(eta) < 10.30) * (0.228047) +
    (0.30 <= eta && eta < 0.40) * (10.30 <= pt * cosh(eta) && pt * cosh(eta) < 10.80) * (0.249350) +
    (0.30 <= eta && eta < 0.40) * (10.80 <= pt * cosh(eta) && pt * cosh(eta) < 11.40) * (0.270616) +
    (0.30 <= eta && eta < 0.40) * (11.40 <= pt * cosh(eta) && pt * cosh(eta) < 12.10) * (0.292865) +
    (0.30 <= eta && eta < 0.40) * (12.10 <= pt * cosh(eta) && pt * cosh(eta) < 12.90) * (0.315110) +
    (0.30 <= eta && eta < 0.40) * (12.90 <= pt * cosh(eta) && pt * cosh(eta) < 13.80) * (0.336546) +
    (0.30 <= eta && eta < 0.40) * (13.80 <= pt * cosh(eta) && pt * cosh(eta) < 14.90) * (0.357525) +
    (0.30 <= eta && eta < 0.40) * (14.90 <= pt * cosh(eta) && pt * cosh(eta) < 16.30) * (0.378695) +
    (0.30 <= eta && eta < 0.40) * (16.30 <= pt * cosh(eta) && pt * cosh(eta) < 18.20) * (0.400236) +
    (0.30 <= eta && eta < 0.40) * (18.20 <= pt * cosh(eta) && pt * cosh(eta) < 21.00) * (0.422294) +
    (0.30 <= eta && eta < 0.40) * (21.00 <= pt * cosh(eta) && pt * cosh(eta) < 25.70) * (0.444843) +
    (0.30 <= eta && eta < 0.40) * (25.70 <= pt * cosh(eta) && pt * cosh(eta) < 38.60) * (0.469965) +
    (0.30 <= eta && eta < 0.40) * (38.60 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.484558) +
    (0.40 <= eta && eta < 0.50) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 5.70) * (0.000629) +
    (0.40 <= eta && eta < 0.50) * (5.70 <= pt * cosh(eta) && pt * cosh(eta) < 6.30) * (0.019727) +
    (0.40 <= eta && eta < 0.50) * (6.30 <= pt * cosh(eta) && pt * cosh(eta) < 6.80) * (0.041772) +
    (0.40 <= eta && eta < 0.50) * (6.80 <= pt * cosh(eta) && pt * cosh(eta) < 7.20) * (0.064756) +
    (0.40 <= eta && eta < 0.50) * (7.20 <= pt * cosh(eta) && pt * cosh(eta) < 7.60) * (0.087848) +
    (0.40 <= eta && eta < 0.50) * (7.60 <= pt * cosh(eta) && pt * cosh(eta) < 8.00) * (0.112141) +
    (0.40 <= eta && eta < 0.50) * (8.00 <= pt * cosh(eta) && pt * cosh(eta) < 8.40) * (0.135408) +
    (0.40 <= eta && eta < 0.50) * (8.40 <= pt * cosh(eta) && pt * cosh(eta) < 8.80) * (0.157790) +
    (0.40 <= eta && eta < 0.50) * (8.80 <= pt * cosh(eta) && pt * cosh(eta) < 9.20) * (0.180619) +
    (0.40 <= eta && eta < 0.50) * (9.20 <= pt * cosh(eta) && pt * cosh(eta) < 9.70) * (0.203357) +
    (0.40 <= eta && eta < 0.50) * (9.70 <= pt * cosh(eta) && pt * cosh(eta) < 10.20) * (0.227838) +
    (0.40 <= eta && eta < 0.50) * (10.20 <= pt * cosh(eta) && pt * cosh(eta) < 10.70) * (0.249584) +
    (0.40 <= eta && eta < 0.50) * (10.70 <= pt * cosh(eta) && pt * cosh(eta) < 11.30) * (0.271029) +
    (0.40 <= eta && eta < 0.50) * (11.30 <= pt * cosh(eta) && pt * cosh(eta) < 12.00) * (0.293437) +
    (0.40 <= eta && eta < 0.50) * (12.00 <= pt * cosh(eta) && pt * cosh(eta) < 12.80) * (0.315808) +
    (0.40 <= eta && eta < 0.50) * (12.80 <= pt * cosh(eta) && pt * cosh(eta) < 13.70) * (0.337334) +
    (0.40 <= eta && eta < 0.50) * (13.70 <= pt * cosh(eta) && pt * cosh(eta) < 14.80) * (0.358365) +
    (0.40 <= eta && eta < 0.50) * (14.80 <= pt * cosh(eta) && pt * cosh(eta) < 16.20) * (0.379551) +
    (0.40 <= eta && eta < 0.50) * (16.20 <= pt * cosh(eta) && pt * cosh(eta) < 18.10) * (0.401066) +
    (0.40 <= eta && eta < 0.50) * (18.10 <= pt * cosh(eta) && pt * cosh(eta) < 20.90) * (0.423052) +
    (0.40 <= eta && eta < 0.50) * (20.90 <= pt * cosh(eta) && pt * cosh(eta) < 25.70) * (0.445684) +
    (0.40 <= eta && eta < 0.50) * (25.70 <= pt * cosh(eta) && pt * cosh(eta) < 39.40) * (0.471163) +
    (0.40 <= eta && eta < 0.50) * (39.40 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.485174) +
    (0.50 <= eta && eta < 0.60) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 5.80) * (0.000656) +
    (0.50 <= eta && eta < 0.60) * (5.80 <= pt * cosh(eta) && pt * cosh(eta) < 6.40) * (0.020312) +
    (0.50 <= eta && eta < 0.60) * (6.40 <= pt * cosh(eta) && pt * cosh(eta) < 6.90) * (0.041999) +
    (0.50 <= eta && eta < 0.60) * (6.90 <= pt * cosh(eta) && pt * cosh(eta) < 7.30) * (0.064530) +
    (0.50 <= eta && eta < 0.60) * (7.30 <= pt * cosh(eta) && pt * cosh(eta) < 7.70) * (0.087046) +
    (0.50 <= eta && eta < 0.60) * (7.70 <= pt * cosh(eta) && pt * cosh(eta) < 8.10) * (0.110549) +
    (0.50 <= eta && eta < 0.60) * (8.10 <= pt * cosh(eta) && pt * cosh(eta) < 8.50) * (0.133710) +
    (0.50 <= eta && eta < 0.60) * (8.50 <= pt * cosh(eta) && pt * cosh(eta) < 8.90) * (0.156249) +
    (0.50 <= eta && eta < 0.60) * (8.90 <= pt * cosh(eta) && pt * cosh(eta) < 9.30) * (0.178337) +
    (0.50 <= eta && eta < 0.60) * (9.30 <= pt * cosh(eta) && pt * cosh(eta) < 9.70) * (0.199630) +
    (0.50 <= eta && eta < 0.60) * (9.70 <= pt * cosh(eta) && pt * cosh(eta) < 10.20) * (0.220540) +
    (0.50 <= eta && eta < 0.60) * (10.20 <= pt * cosh(eta) && pt * cosh(eta) < 10.70) * (0.242222) +
    (0.50 <= eta && eta < 0.60) * (10.70 <= pt * cosh(eta) && pt * cosh(eta) < 11.30) * (0.264089) +
    (0.50 <= eta && eta < 0.60) * (11.30 <= pt * cosh(eta) && pt * cosh(eta) < 12.00) * (0.287002) +
    (0.50 <= eta && eta < 0.60) * (12.00 <= pt * cosh(eta) && pt * cosh(eta) < 12.80) * (0.309935) +
    (0.50 <= eta && eta < 0.60) * (12.80 <= pt * cosh(eta) && pt * cosh(eta) < 13.70) * (0.332047) +
    (0.50 <= eta && eta < 0.60) * (13.70 <= pt * cosh(eta) && pt * cosh(eta) < 14.80) * (0.353688) +
    (0.50 <= eta && eta < 0.60) * (14.80 <= pt * cosh(eta) && pt * cosh(eta) < 16.20) * (0.375520) +
    (0.50 <= eta && eta < 0.60) * (16.20 <= pt * cosh(eta) && pt * cosh(eta) < 18.00) * (0.397160) +
    (0.50 <= eta && eta < 0.60) * (18.00 <= pt * cosh(eta) && pt * cosh(eta) < 20.60) * (0.418834) +
    (0.50 <= eta && eta < 0.60) * (20.60 <= pt * cosh(eta) && pt * cosh(eta) < 24.90) * (0.441179) +
    (0.50 <= eta && eta < 0.60) * (24.90 <= pt * cosh(eta) && pt * cosh(eta) < 35.10) * (0.465437) +
    (0.50 <= eta && eta < 0.60) * (35.10 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.482784) +
    (0.60 <= eta && eta < 0.70) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 5.80) * (0.000610) +
    (0.60 <= eta && eta < 0.70) * (5.80 <= pt * cosh(eta) && pt * cosh(eta) < 6.40) * (0.019574) +
    (0.60 <= eta && eta < 0.70) * (6.40 <= pt * cosh(eta) && pt * cosh(eta) < 6.90) * (0.040602) +
    (0.60 <= eta && eta < 0.70) * (6.90 <= pt * cosh(eta) && pt * cosh(eta) < 7.30) * (0.063109) +
    (0.60 <= eta && eta < 0.70) * (7.30 <= pt * cosh(eta) && pt * cosh(eta) < 7.70) * (0.085463) +
    (0.60 <= eta && eta < 0.70) * (7.70 <= pt * cosh(eta) && pt * cosh(eta) < 8.10) * (0.108771) +
    (0.60 <= eta && eta < 0.70) * (8.10 <= pt * cosh(eta) && pt * cosh(eta) < 8.50) * (0.132650) +
    (0.60 <= eta && eta < 0.70) * (8.50 <= pt * cosh(eta) && pt * cosh(eta) < 8.90) * (0.154383) +
    (0.60 <= eta && eta < 0.70) * (8.90 <= pt * cosh(eta) && pt * cosh(eta) < 9.30) * (0.175861) +
    (0.60 <= eta && eta < 0.70) * (9.30 <= pt * cosh(eta) && pt * cosh(eta) < 9.80) * (0.198600) +
    (0.60 <= eta && eta < 0.70) * (9.80 <= pt * cosh(eta) && pt * cosh(eta) < 10.30) * (0.222323) +
    (0.60 <= eta && eta < 0.70) * (10.30 <= pt * cosh(eta) && pt * cosh(eta) < 10.80) * (0.244096) +
    (0.60 <= eta && eta < 0.70) * (10.80 <= pt * cosh(eta) && pt * cosh(eta) < 11.40) * (0.265662) +
    (0.60 <= eta && eta < 0.70) * (11.40 <= pt * cosh(eta) && pt * cosh(eta) < 12.10) * (0.288269) +
    (0.60 <= eta && eta < 0.70) * (12.10 <= pt * cosh(eta) && pt * cosh(eta) < 12.90) * (0.310912) +
    (0.60 <= eta && eta < 0.70) * (12.90 <= pt * cosh(eta) && pt * cosh(eta) < 13.80) * (0.332765) +
    (0.60 <= eta && eta < 0.70) * (13.80 <= pt * cosh(eta) && pt * cosh(eta) < 14.90) * (0.354177) +
    (0.60 <= eta && eta < 0.70) * (14.90 <= pt * cosh(eta) && pt * cosh(eta) < 16.30) * (0.375807) +
    (0.60 <= eta && eta < 0.70) * (16.30 <= pt * cosh(eta) && pt * cosh(eta) < 18.10) * (0.397278) +
    (0.60 <= eta && eta < 0.70) * (18.10 <= pt * cosh(eta) && pt * cosh(eta) < 20.70) * (0.418823) +
    (0.60 <= eta && eta < 0.70) * (20.70 <= pt * cosh(eta) && pt * cosh(eta) < 25.00) * (0.441082) +
    (0.60 <= eta && eta < 0.70) * (25.00 <= pt * cosh(eta) && pt * cosh(eta) < 35.20) * (0.465310) +
    (0.60 <= eta && eta < 0.70) * (35.20 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.482651) +
    (0.70 <= eta && eta < 0.80) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 5.90) * (0.000768) +
    (0.70 <= eta && eta < 0.80) * (5.90 <= pt * cosh(eta) && pt * cosh(eta) < 6.50) * (0.022252) +
    (0.70 <= eta && eta < 0.80) * (6.50 <= pt * cosh(eta) && pt * cosh(eta) < 6.90) * (0.042122) +
    (0.70 <= eta && eta < 0.80) * (6.90 <= pt * cosh(eta) && pt * cosh(eta) < 7.30) * (0.062694) +
    (0.70 <= eta && eta < 0.80) * (7.30 <= pt * cosh(eta) && pt * cosh(eta) < 7.70) * (0.085610) +
    (0.70 <= eta && eta < 0.80) * (7.70 <= pt * cosh(eta) && pt * cosh(eta) < 8.10) * (0.107607) +
    (0.70 <= eta && eta < 0.80) * (8.10 <= pt * cosh(eta) && pt * cosh(eta) < 8.50) * (0.130626) +
    (0.70 <= eta && eta < 0.80) * (8.50 <= pt * cosh(eta) && pt * cosh(eta) < 8.90) * (0.154667) +
    (0.70 <= eta && eta < 0.80) * (8.90 <= pt * cosh(eta) && pt * cosh(eta) < 9.30) * (0.176433) +
    (0.70 <= eta && eta < 0.80) * (9.30 <= pt * cosh(eta) && pt * cosh(eta) < 9.80) * (0.198617) +
    (0.70 <= eta && eta < 0.80) * (9.80 <= pt * cosh(eta) && pt * cosh(eta) < 10.30) * (0.222620) +
    (0.70 <= eta && eta < 0.80) * (10.30 <= pt * cosh(eta) && pt * cosh(eta) < 10.80) * (0.244391) +
    (0.70 <= eta && eta < 0.80) * (10.80 <= pt * cosh(eta) && pt * cosh(eta) < 11.40) * (0.265940) +
    (0.70 <= eta && eta < 0.80) * (11.40 <= pt * cosh(eta) && pt * cosh(eta) < 12.10) * (0.288527) +
    (0.70 <= eta && eta < 0.80) * (12.10 <= pt * cosh(eta) && pt * cosh(eta) < 12.90) * (0.311148) +
    (0.70 <= eta && eta < 0.80) * (12.90 <= pt * cosh(eta) && pt * cosh(eta) < 13.80) * (0.332977) +
    (0.70 <= eta && eta < 0.80) * (13.80 <= pt * cosh(eta) && pt * cosh(eta) < 14.90) * (0.354366) +
    (0.70 <= eta && eta < 0.80) * (14.90 <= pt * cosh(eta) && pt * cosh(eta) < 16.30) * (0.375969) +
    (0.70 <= eta && eta < 0.80) * (16.30 <= pt * cosh(eta) && pt * cosh(eta) < 18.20) * (0.397969) +
    (0.70 <= eta && eta < 0.80) * (18.20 <= pt * cosh(eta) && pt * cosh(eta) < 20.90) * (0.420136) +
    (0.70 <= eta && eta < 0.80) * (20.90 <= pt * cosh(eta) && pt * cosh(eta) < 25.40) * (0.442635) +
    (0.70 <= eta && eta < 0.80) * (25.40 <= pt * cosh(eta) && pt * cosh(eta) < 36.60) * (0.467206) +
    (0.70 <= eta && eta < 0.80) * (36.60 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.483337) +
    (0.80 <= eta && eta < 0.90) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 6.00) * (0.000728) +
    (0.80 <= eta && eta < 0.90) * (6.00 <= pt * cosh(eta) && pt * cosh(eta) < 6.60) * (0.021330) +
    (0.80 <= eta && eta < 0.90) * (6.60 <= pt * cosh(eta) && pt * cosh(eta) < 7.10) * (0.042978) +
    (0.80 <= eta && eta < 0.90) * (7.10 <= pt * cosh(eta) && pt * cosh(eta) < 7.50) * (0.065797) +
    (0.80 <= eta && eta < 0.90) * (7.50 <= pt * cosh(eta) && pt * cosh(eta) < 7.90) * (0.087589) +
    (0.80 <= eta && eta < 0.90) * (7.90 <= pt * cosh(eta) && pt * cosh(eta) < 8.30) * (0.109904) +
    (0.80 <= eta && eta < 0.90) * (8.30 <= pt * cosh(eta) && pt * cosh(eta) < 8.70) * (0.132432) +
    (0.80 <= eta && eta < 0.90) * (8.70 <= pt * cosh(eta) && pt * cosh(eta) < 9.10) * (0.155030) +
    (0.80 <= eta && eta < 0.90) * (9.10 <= pt * cosh(eta) && pt * cosh(eta) < 9.50) * (0.177015) +
    (0.80 <= eta && eta < 0.90) * (9.50 <= pt * cosh(eta) && pt * cosh(eta) < 10.00) * (0.198109) +
    (0.80 <= eta && eta < 0.90) * (10.00 <= pt * cosh(eta) && pt * cosh(eta) < 10.50) * (0.220572) +
    (0.80 <= eta && eta < 0.90) * (10.50 <= pt * cosh(eta) && pt * cosh(eta) < 11.00) * (0.241975) +
    (0.80 <= eta && eta < 0.90) * (11.00 <= pt * cosh(eta) && pt * cosh(eta) < 11.60) * (0.263275) +
    (0.80 <= eta && eta < 0.90) * (11.60 <= pt * cosh(eta) && pt * cosh(eta) < 12.30) * (0.285672) +
    (0.80 <= eta && eta < 0.90) * (12.30 <= pt * cosh(eta) && pt * cosh(eta) < 13.10) * (0.308182) +
    (0.80 <= eta && eta < 0.90) * (13.10 <= pt * cosh(eta) && pt * cosh(eta) < 14.00) * (0.329981) +
    (0.80 <= eta && eta < 0.90) * (14.00 <= pt * cosh(eta) && pt * cosh(eta) < 15.10) * (0.351418) +
    (0.80 <= eta && eta < 0.90) * (15.10 <= pt * cosh(eta) && pt * cosh(eta) < 16.50) * (0.373154) +
    (0.80 <= eta && eta < 0.90) * (16.50 <= pt * cosh(eta) && pt * cosh(eta) < 18.30) * (0.394818) +
    (0.80 <= eta && eta < 0.90) * (18.30 <= pt * cosh(eta) && pt * cosh(eta) < 20.90) * (0.416648) +
    (0.80 <= eta && eta < 0.90) * (20.90 <= pt * cosh(eta) && pt * cosh(eta) < 25.10) * (0.439069) +
    (0.80 <= eta && eta < 0.90) * (25.10 <= pt * cosh(eta) && pt * cosh(eta) < 34.50) * (0.463037) +
    (0.80 <= eta && eta < 0.90) * (34.50 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.481439) +
    (0.90 <= eta && eta < 1.00) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 6.00) * (0.000687) +
    (0.90 <= eta && eta < 1.00) * (6.00 <= pt * cosh(eta) && pt * cosh(eta) < 6.60) * (0.020501) +
    (0.90 <= eta && eta < 1.00) * (6.60 <= pt * cosh(eta) && pt * cosh(eta) < 7.10) * (0.041727) +
    (0.90 <= eta && eta < 1.00) * (7.10 <= pt * cosh(eta) && pt * cosh(eta) < 7.50) * (0.063800) +
    (0.90 <= eta && eta < 1.00) * (7.50 <= pt * cosh(eta) && pt * cosh(eta) < 7.90) * (0.084994) +
    (0.90 <= eta && eta < 1.00) * (7.90 <= pt * cosh(eta) && pt * cosh(eta) < 8.30) * (0.107765) +
    (0.90 <= eta && eta < 1.00) * (8.30 <= pt * cosh(eta) && pt * cosh(eta) < 8.70) * (0.130564) +
    (0.90 <= eta && eta < 1.00) * (8.70 <= pt * cosh(eta) && pt * cosh(eta) < 9.10) * (0.152720) +
    (0.90 <= eta && eta < 1.00) * (9.10 <= pt * cosh(eta) && pt * cosh(eta) < 9.50) * (0.174241) +
    (0.90 <= eta && eta < 1.00) * (9.50 <= pt * cosh(eta) && pt * cosh(eta) < 10.00) * (0.197247) +
    (0.90 <= eta && eta < 1.00) * (10.00 <= pt * cosh(eta) && pt * cosh(eta) < 10.50) * (0.220892) +
    (0.90 <= eta && eta < 1.00) * (10.50 <= pt * cosh(eta) && pt * cosh(eta) < 11.00) * (0.242281) +
    (0.90 <= eta && eta < 1.00) * (11.00 <= pt * cosh(eta) && pt * cosh(eta) < 11.60) * (0.263564) +
    (0.90 <= eta && eta < 1.00) * (11.60 <= pt * cosh(eta) && pt * cosh(eta) < 12.30) * (0.285942) +
    (0.90 <= eta && eta < 1.00) * (12.30 <= pt * cosh(eta) && pt * cosh(eta) < 13.10) * (0.308428) +
    (0.90 <= eta && eta < 1.00) * (13.10 <= pt * cosh(eta) && pt * cosh(eta) < 14.00) * (0.330204) +
    (0.90 <= eta && eta < 1.00) * (14.00 <= pt * cosh(eta) && pt * cosh(eta) < 15.10) * (0.351616) +
    (0.90 <= eta && eta < 1.00) * (15.10 <= pt * cosh(eta) && pt * cosh(eta) < 16.50) * (0.373325) +
    (0.90 <= eta && eta < 1.00) * (16.50 <= pt * cosh(eta) && pt * cosh(eta) < 18.30) * (0.394961) +
    (0.90 <= eta && eta < 1.00) * (18.30 <= pt * cosh(eta) && pt * cosh(eta) < 20.90) * (0.416763) +
    (0.90 <= eta && eta < 1.00) * (20.90 <= pt * cosh(eta) && pt * cosh(eta) < 25.10) * (0.439153) +
    (0.90 <= eta && eta < 1.00) * (25.10 <= pt * cosh(eta) && pt * cosh(eta) < 34.50) * (0.463089) +
    (0.90 <= eta && eta < 1.00) * (34.50 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.481465) +
    (1.00 <= eta && eta < 1.10) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.000000) }

  add EfficiencyFormula {321} {2212} {
    (eta<-1.00 || eta>=1.00 || pt * cosh(eta) < 0.10 || pt * cosh(eta) >= 50.00) * ( 0.00 ) +
    (-1.00 <= eta && eta < -0.90) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 10.30) * (0.000705) +
    (-1.00 <= eta && eta < -0.90) * (10.30 <= pt * cosh(eta) && pt * cosh(eta) < 11.40) * (0.021025) +
    (-1.00 <= eta && eta < -0.90) * (11.40 <= pt * cosh(eta) && pt * cosh(eta) < 12.20) * (0.042481) +
    (-1.00 <= eta && eta < -0.90) * (12.20 <= pt * cosh(eta) && pt * cosh(eta) < 12.90) * (0.063884) +
    (-1.00 <= eta && eta < -0.90) * (12.90 <= pt * cosh(eta) && pt * cosh(eta) < 13.60) * (0.085969) +
    (-1.00 <= eta && eta < -0.90) * (13.60 <= pt * cosh(eta) && pt * cosh(eta) < 14.30) * (0.108965) +
    (-1.00 <= eta && eta < -0.90) * (14.30 <= pt * cosh(eta) && pt * cosh(eta) < 15.00) * (0.132030) +
    (-1.00 <= eta && eta < -0.90) * (15.00 <= pt * cosh(eta) && pt * cosh(eta) < 15.70) * (0.154564) +
    (-1.00 <= eta && eta < -0.90) * (15.70 <= pt * cosh(eta) && pt * cosh(eta) < 16.40) * (0.176174) +
    (-1.00 <= eta && eta < -0.90) * (16.40 <= pt * cosh(eta) && pt * cosh(eta) < 17.20) * (0.198024) +
    (-1.00 <= eta && eta < -0.90) * (17.20 <= pt * cosh(eta) && pt * cosh(eta) < 18.00) * (0.219729) +
    (-1.00 <= eta && eta < -0.90) * (18.00 <= pt * cosh(eta) && pt * cosh(eta) < 18.90) * (0.240906) +
    (-1.00 <= eta && eta < -0.90) * (18.90 <= pt * cosh(eta) && pt * cosh(eta) < 19.90) * (0.262363) +
    (-1.00 <= eta && eta < -0.90) * (19.90 <= pt * cosh(eta) && pt * cosh(eta) < 21.00) * (0.283549) +
    (-1.00 <= eta && eta < -0.90) * (21.00 <= pt * cosh(eta) && pt * cosh(eta) < 22.30) * (0.304817) +
    (-1.00 <= eta && eta < -0.90) * (22.30 <= pt * cosh(eta) && pt * cosh(eta) < 23.80) * (0.326209) +
    (-1.00 <= eta && eta < -0.90) * (23.80 <= pt * cosh(eta) && pt * cosh(eta) < 25.60) * (0.347442) +
    (-1.00 <= eta && eta < -0.90) * (25.60 <= pt * cosh(eta) && pt * cosh(eta) < 27.80) * (0.368549) +
    (-1.00 <= eta && eta < -0.90) * (27.80 <= pt * cosh(eta) && pt * cosh(eta) < 30.70) * (0.389801) +
    (-1.00 <= eta && eta < -0.90) * (30.70 <= pt * cosh(eta) && pt * cosh(eta) < 34.70) * (0.411336) +
    (-1.00 <= eta && eta < -0.90) * (34.70 <= pt * cosh(eta) && pt * cosh(eta) < 40.90) * (0.433234) +
    (-1.00 <= eta && eta < -0.90) * (40.90 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.453570) +
    (-0.90 <= eta && eta < -0.80) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 10.10) * (0.000649) +
    (-0.90 <= eta && eta < -0.80) * (10.10 <= pt * cosh(eta) && pt * cosh(eta) < 11.20) * (0.020231) +
    (-0.90 <= eta && eta < -0.80) * (11.20 <= pt * cosh(eta) && pt * cosh(eta) < 12.00) * (0.041709) +
    (-0.90 <= eta && eta < -0.80) * (12.00 <= pt * cosh(eta) && pt * cosh(eta) < 12.70) * (0.063327) +
    (-0.90 <= eta && eta < -0.80) * (12.70 <= pt * cosh(eta) && pt * cosh(eta) < 13.40) * (0.085708) +
    (-0.90 <= eta && eta < -0.80) * (13.40 <= pt * cosh(eta) && pt * cosh(eta) < 14.00) * (0.107361) +
    (-0.90 <= eta && eta < -0.80) * (14.00 <= pt * cosh(eta) && pt * cosh(eta) < 14.70) * (0.129120) +
    (-0.90 <= eta && eta < -0.80) * (14.70 <= pt * cosh(eta) && pt * cosh(eta) < 15.40) * (0.152075) +
    (-0.90 <= eta && eta < -0.80) * (15.40 <= pt * cosh(eta) && pt * cosh(eta) < 16.10) * (0.174114) +
    (-0.90 <= eta && eta < -0.80) * (16.10 <= pt * cosh(eta) && pt * cosh(eta) < 16.90) * (0.196407) +
    (-0.90 <= eta && eta < -0.80) * (16.90 <= pt * cosh(eta) && pt * cosh(eta) < 17.70) * (0.218543) +
    (-0.90 <= eta && eta < -0.80) * (17.70 <= pt * cosh(eta) && pt * cosh(eta) < 18.60) * (0.240121) +
    (-0.90 <= eta && eta < -0.80) * (18.60 <= pt * cosh(eta) && pt * cosh(eta) < 19.60) * (0.261958) +
    (-0.90 <= eta && eta < -0.80) * (19.60 <= pt * cosh(eta) && pt * cosh(eta) < 20.70) * (0.283483) +
    (-0.90 <= eta && eta < -0.80) * (20.70 <= pt * cosh(eta) && pt * cosh(eta) < 22.00) * (0.305052) +
    (-0.90 <= eta && eta < -0.80) * (22.00 <= pt * cosh(eta) && pt * cosh(eta) < 23.50) * (0.326699) +
    (-0.90 <= eta && eta < -0.80) * (23.50 <= pt * cosh(eta) && pt * cosh(eta) < 25.30) * (0.348135) +
    (-0.90 <= eta && eta < -0.80) * (25.30 <= pt * cosh(eta) && pt * cosh(eta) < 27.50) * (0.369386) +
    (-0.90 <= eta && eta < -0.80) * (27.50 <= pt * cosh(eta) && pt * cosh(eta) < 30.40) * (0.390721) +
    (-0.90 <= eta && eta < -0.80) * (30.40 <= pt * cosh(eta) && pt * cosh(eta) < 34.40) * (0.412270) +
    (-0.90 <= eta && eta < -0.80) * (34.40 <= pt * cosh(eta) && pt * cosh(eta) < 40.70) * (0.434265) +
    (-0.90 <= eta && eta < -0.80) * (40.70 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.454684) +
    (-0.80 <= eta && eta < -0.70) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 10.00) * (0.000698) +
    (-0.80 <= eta && eta < -0.70) * (10.00 <= pt * cosh(eta) && pt * cosh(eta) < 11.00) * (0.020156) +
    (-0.80 <= eta && eta < -0.70) * (11.00 <= pt * cosh(eta) && pt * cosh(eta) < 11.80) * (0.040719) +
    (-0.80 <= eta && eta < -0.70) * (11.80 <= pt * cosh(eta) && pt * cosh(eta) < 12.50) * (0.062503) +
    (-0.80 <= eta && eta < -0.70) * (12.50 <= pt * cosh(eta) && pt * cosh(eta) < 13.10) * (0.083481) +
    (-0.80 <= eta && eta < -0.70) * (13.10 <= pt * cosh(eta) && pt * cosh(eta) < 13.70) * (0.103705) +
    (-0.80 <= eta && eta < -0.70) * (13.70 <= pt * cosh(eta) && pt * cosh(eta) < 14.30) * (0.124110) +
    (-0.80 <= eta && eta < -0.70) * (14.30 <= pt * cosh(eta) && pt * cosh(eta) < 15.00) * (0.145883) +
    (-0.80 <= eta && eta < -0.70) * (15.00 <= pt * cosh(eta) && pt * cosh(eta) < 15.70) * (0.168517) +
    (-0.80 <= eta && eta < -0.70) * (15.70 <= pt * cosh(eta) && pt * cosh(eta) < 16.40) * (0.189997) +
    (-0.80 <= eta && eta < -0.70) * (16.40 <= pt * cosh(eta) && pt * cosh(eta) < 17.20) * (0.211530) +
    (-0.80 <= eta && eta < -0.70) * (17.20 <= pt * cosh(eta) && pt * cosh(eta) < 18.10) * (0.234017) +
    (-0.80 <= eta && eta < -0.70) * (18.10 <= pt * cosh(eta) && pt * cosh(eta) < 19.00) * (0.255655) +
    (-0.80 <= eta && eta < -0.70) * (19.00 <= pt * cosh(eta) && pt * cosh(eta) < 20.10) * (0.277170) +
    (-0.80 <= eta && eta < -0.70) * (20.10 <= pt * cosh(eta) && pt * cosh(eta) < 21.30) * (0.298995) +
    (-0.80 <= eta && eta < -0.70) * (21.30 <= pt * cosh(eta) && pt * cosh(eta) < 22.70) * (0.320333) +
    (-0.80 <= eta && eta < -0.70) * (22.70 <= pt * cosh(eta) && pt * cosh(eta) < 24.40) * (0.341878) +
    (-0.80 <= eta && eta < -0.70) * (24.40 <= pt * cosh(eta) && pt * cosh(eta) < 26.50) * (0.363611) +
    (-0.80 <= eta && eta < -0.70) * (26.50 <= pt * cosh(eta) && pt * cosh(eta) < 29.10) * (0.384994) +
    (-0.80 <= eta && eta < -0.70) * (29.10 <= pt * cosh(eta) && pt * cosh(eta) < 32.70) * (0.406381) +
    (-0.80 <= eta && eta < -0.70) * (32.70 <= pt * cosh(eta) && pt * cosh(eta) < 38.10) * (0.428249) +
    (-0.80 <= eta && eta < -0.70) * (38.10 <= pt * cosh(eta) && pt * cosh(eta) < 47.90) * (0.450901) +
    (-0.80 <= eta && eta < -0.70) * (47.90 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.462545) +
    (-0.70 <= eta && eta < -0.60) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 9.90) * (0.000644) +
    (-0.70 <= eta && eta < -0.60) * (9.90 <= pt * cosh(eta) && pt * cosh(eta) < 10.90) * (0.019416) +
    (-0.70 <= eta && eta < -0.60) * (10.90 <= pt * cosh(eta) && pt * cosh(eta) < 11.70) * (0.039800) +
    (-0.70 <= eta && eta < -0.60) * (11.70 <= pt * cosh(eta) && pt * cosh(eta) < 12.40) * (0.061565) +
    (-0.70 <= eta && eta < -0.60) * (12.40 <= pt * cosh(eta) && pt * cosh(eta) < 13.00) * (0.082602) +
    (-0.70 <= eta && eta < -0.60) * (13.00 <= pt * cosh(eta) && pt * cosh(eta) < 13.60) * (0.102923) +
    (-0.70 <= eta && eta < -0.60) * (13.60 <= pt * cosh(eta) && pt * cosh(eta) < 14.20) * (0.123448) +
    (-0.70 <= eta && eta < -0.60) * (14.20 <= pt * cosh(eta) && pt * cosh(eta) < 14.90) * (0.145363) +
    (-0.70 <= eta && eta < -0.60) * (14.90 <= pt * cosh(eta) && pt * cosh(eta) < 15.60) * (0.168148) +
    (-0.70 <= eta && eta < -0.60) * (15.60 <= pt * cosh(eta) && pt * cosh(eta) < 16.30) * (0.189770) +
    (-0.70 <= eta && eta < -0.60) * (16.30 <= pt * cosh(eta) && pt * cosh(eta) < 17.10) * (0.211440) +
    (-0.70 <= eta && eta < -0.60) * (17.10 <= pt * cosh(eta) && pt * cosh(eta) < 17.90) * (0.232802) +
    (-0.70 <= eta && eta < -0.60) * (17.90 <= pt * cosh(eta) && pt * cosh(eta) < 18.80) * (0.253505) +
    (-0.70 <= eta && eta < -0.60) * (18.80 <= pt * cosh(eta) && pt * cosh(eta) < 19.80) * (0.274359) +
    (-0.70 <= eta && eta < -0.60) * (19.80 <= pt * cosh(eta) && pt * cosh(eta) < 21.00) * (0.295736) +
    (-0.70 <= eta && eta < -0.60) * (21.00 <= pt * cosh(eta) && pt * cosh(eta) < 22.40) * (0.317652) +
    (-0.70 <= eta && eta < -0.60) * (22.40 <= pt * cosh(eta) && pt * cosh(eta) < 24.00) * (0.339116) +
    (-0.70 <= eta && eta < -0.60) * (24.00 <= pt * cosh(eta) && pt * cosh(eta) < 26.00) * (0.360427) +
    (-0.70 <= eta && eta < -0.60) * (26.00 <= pt * cosh(eta) && pt * cosh(eta) < 28.50) * (0.381782) +
    (-0.70 <= eta && eta < -0.60) * (28.50 <= pt * cosh(eta) && pt * cosh(eta) < 31.90) * (0.403198) +
    (-0.70 <= eta && eta < -0.60) * (31.90 <= pt * cosh(eta) && pt * cosh(eta) < 36.90) * (0.424964) +
    (-0.70 <= eta && eta < -0.60) * (36.90 <= pt * cosh(eta) && pt * cosh(eta) < 45.70) * (0.447490) +
    (-0.70 <= eta && eta < -0.60) * (45.70 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.461199) +
    (-0.60 <= eta && eta < -0.50) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 9.90) * (0.000623) +
    (-0.60 <= eta && eta < -0.50) * (9.90 <= pt * cosh(eta) && pt * cosh(eta) < 10.90) * (0.019050) +
    (-0.60 <= eta && eta < -0.50) * (10.90 <= pt * cosh(eta) && pt * cosh(eta) < 11.70) * (0.039230) +
    (-0.60 <= eta && eta < -0.50) * (11.70 <= pt * cosh(eta) && pt * cosh(eta) < 12.40) * (0.060853) +
    (-0.60 <= eta && eta < -0.50) * (12.40 <= pt * cosh(eta) && pt * cosh(eta) < 13.00) * (0.081798) +
    (-0.60 <= eta && eta < -0.50) * (13.00 <= pt * cosh(eta) && pt * cosh(eta) < 13.60) * (0.102059) +
    (-0.60 <= eta && eta < -0.50) * (13.60 <= pt * cosh(eta) && pt * cosh(eta) < 14.20) * (0.122547) +
    (-0.60 <= eta && eta < -0.50) * (14.20 <= pt * cosh(eta) && pt * cosh(eta) < 14.90) * (0.144443) +
    (-0.60 <= eta && eta < -0.50) * (14.90 <= pt * cosh(eta) && pt * cosh(eta) < 15.60) * (0.167226) +
    (-0.60 <= eta && eta < -0.50) * (15.60 <= pt * cosh(eta) && pt * cosh(eta) < 16.30) * (0.188860) +
    (-0.60 <= eta && eta < -0.50) * (16.30 <= pt * cosh(eta) && pt * cosh(eta) < 17.10) * (0.210554) +
    (-0.60 <= eta && eta < -0.50) * (17.10 <= pt * cosh(eta) && pt * cosh(eta) < 17.90) * (0.231950) +
    (-0.60 <= eta && eta < -0.50) * (17.90 <= pt * cosh(eta) && pt * cosh(eta) < 18.80) * (0.252693) +
    (-0.60 <= eta && eta < -0.50) * (18.80 <= pt * cosh(eta) && pt * cosh(eta) < 19.80) * (0.273595) +
    (-0.60 <= eta && eta < -0.50) * (19.80 <= pt * cosh(eta) && pt * cosh(eta) < 21.00) * (0.295028) +
    (-0.60 <= eta && eta < -0.50) * (21.00 <= pt * cosh(eta) && pt * cosh(eta) < 22.40) * (0.317006) +
    (-0.60 <= eta && eta < -0.50) * (22.40 <= pt * cosh(eta) && pt * cosh(eta) < 24.00) * (0.338536) +
    (-0.60 <= eta && eta < -0.50) * (24.00 <= pt * cosh(eta) && pt * cosh(eta) < 26.00) * (0.359915) +
    (-0.60 <= eta && eta < -0.50) * (26.00 <= pt * cosh(eta) && pt * cosh(eta) < 28.50) * (0.381344) +
    (-0.60 <= eta && eta < -0.50) * (28.50 <= pt * cosh(eta) && pt * cosh(eta) < 31.90) * (0.402835) +
    (-0.60 <= eta && eta < -0.50) * (31.90 <= pt * cosh(eta) && pt * cosh(eta) < 36.90) * (0.424681) +
    (-0.60 <= eta && eta < -0.50) * (36.90 <= pt * cosh(eta) && pt * cosh(eta) < 45.60) * (0.447175) +
    (-0.60 <= eta && eta < -0.50) * (45.60 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.460966) +
    (-0.50 <= eta && eta < -0.40) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 9.80) * (0.000665) +
    (-0.50 <= eta && eta < -0.40) * (9.80 <= pt * cosh(eta) && pt * cosh(eta) < 10.80) * (0.020045) +
    (-0.50 <= eta && eta < -0.40) * (10.80 <= pt * cosh(eta) && pt * cosh(eta) < 11.60) * (0.040989) +
    (-0.50 <= eta && eta < -0.40) * (11.60 <= pt * cosh(eta) && pt * cosh(eta) < 12.30) * (0.063254) +
    (-0.50 <= eta && eta < -0.40) * (12.30 <= pt * cosh(eta) && pt * cosh(eta) < 12.90) * (0.084688) +
    (-0.50 <= eta && eta < -0.40) * (12.90 <= pt * cosh(eta) && pt * cosh(eta) < 13.50) * (0.105322) +
    (-0.50 <= eta && eta < -0.40) * (13.50 <= pt * cosh(eta) && pt * cosh(eta) < 14.10) * (0.126101) +
    (-0.50 <= eta && eta < -0.40) * (14.10 <= pt * cosh(eta) && pt * cosh(eta) < 14.70) * (0.146555) +
    (-0.50 <= eta && eta < -0.40) * (14.70 <= pt * cosh(eta) && pt * cosh(eta) < 15.40) * (0.167958) +
    (-0.50 <= eta && eta < -0.40) * (15.40 <= pt * cosh(eta) && pt * cosh(eta) < 16.10) * (0.189862) +
    (-0.50 <= eta && eta < -0.40) * (16.10 <= pt * cosh(eta) && pt * cosh(eta) < 16.90) * (0.211794) +
    (-0.50 <= eta && eta < -0.40) * (16.90 <= pt * cosh(eta) && pt * cosh(eta) < 17.70) * (0.233391) +
    (-0.50 <= eta && eta < -0.40) * (17.70 <= pt * cosh(eta) && pt * cosh(eta) < 18.60) * (0.254294) +
    (-0.50 <= eta && eta < -0.40) * (18.60 <= pt * cosh(eta) && pt * cosh(eta) < 19.60) * (0.275318) +
    (-0.50 <= eta && eta < -0.40) * (19.60 <= pt * cosh(eta) && pt * cosh(eta) < 20.80) * (0.296835) +
    (-0.50 <= eta && eta < -0.40) * (20.80 <= pt * cosh(eta) && pt * cosh(eta) < 22.20) * (0.318854) +
    (-0.50 <= eta && eta < -0.40) * (22.20 <= pt * cosh(eta) && pt * cosh(eta) < 23.80) * (0.340377) +
    (-0.50 <= eta && eta < -0.40) * (23.80 <= pt * cosh(eta) && pt * cosh(eta) < 25.80) * (0.361703) +
    (-0.50 <= eta && eta < -0.40) * (25.80 <= pt * cosh(eta) && pt * cosh(eta) < 28.30) * (0.383028) +
    (-0.50 <= eta && eta < -0.40) * (28.30 <= pt * cosh(eta) && pt * cosh(eta) < 31.70) * (0.404363) +
    (-0.50 <= eta && eta < -0.40) * (31.70 <= pt * cosh(eta) && pt * cosh(eta) < 36.80) * (0.426192) +
    (-0.50 <= eta && eta < -0.40) * (36.80 <= pt * cosh(eta) && pt * cosh(eta) < 45.80) * (0.448789) +
    (-0.50 <= eta && eta < -0.40) * (45.80 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.462263) +
    (-0.40 <= eta && eta < -0.30) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 9.80) * (0.000645) +
    (-0.40 <= eta && eta < -0.30) * (9.80 <= pt * cosh(eta) && pt * cosh(eta) < 10.80) * (0.019306) +
    (-0.40 <= eta && eta < -0.30) * (10.80 <= pt * cosh(eta) && pt * cosh(eta) < 11.60) * (0.039824) +
    (-0.40 <= eta && eta < -0.30) * (11.60 <= pt * cosh(eta) && pt * cosh(eta) < 12.30) * (0.061801) +
    (-0.40 <= eta && eta < -0.30) * (12.30 <= pt * cosh(eta) && pt * cosh(eta) < 12.90) * (0.083051) +
    (-0.40 <= eta && eta < -0.30) * (12.90 <= pt * cosh(eta) && pt * cosh(eta) < 13.50) * (0.103567) +
    (-0.40 <= eta && eta < -0.30) * (13.50 <= pt * cosh(eta) && pt * cosh(eta) < 14.10) * (0.124275) +
    (-0.40 <= eta && eta < -0.30) * (14.10 <= pt * cosh(eta) && pt * cosh(eta) < 14.70) * (0.144695) +
    (-0.40 <= eta && eta < -0.30) * (14.70 <= pt * cosh(eta) && pt * cosh(eta) < 15.40) * (0.166094) +
    (-0.40 <= eta && eta < -0.30) * (15.40 <= pt * cosh(eta) && pt * cosh(eta) < 16.10) * (0.188023) +
    (-0.40 <= eta && eta < -0.30) * (16.10 <= pt * cosh(eta) && pt * cosh(eta) < 16.90) * (0.210005) +
    (-0.40 <= eta && eta < -0.30) * (16.90 <= pt * cosh(eta) && pt * cosh(eta) < 17.70) * (0.231670) +
    (-0.40 <= eta && eta < -0.30) * (17.70 <= pt * cosh(eta) && pt * cosh(eta) < 18.60) * (0.252656) +
    (-0.40 <= eta && eta < -0.30) * (18.60 <= pt * cosh(eta) && pt * cosh(eta) < 19.60) * (0.273778) +
    (-0.40 <= eta && eta < -0.30) * (19.60 <= pt * cosh(eta) && pt * cosh(eta) < 20.80) * (0.295407) +
    (-0.40 <= eta && eta < -0.30) * (20.80 <= pt * cosh(eta) && pt * cosh(eta) < 22.10) * (0.316782) +
    (-0.40 <= eta && eta < -0.30) * (22.10 <= pt * cosh(eta) && pt * cosh(eta) < 23.70) * (0.337879) +
    (-0.40 <= eta && eta < -0.30) * (23.70 <= pt * cosh(eta) && pt * cosh(eta) < 25.60) * (0.359062) +
    (-0.40 <= eta && eta < -0.30) * (25.60 <= pt * cosh(eta) && pt * cosh(eta) < 28.10) * (0.380434) +
    (-0.40 <= eta && eta < -0.30) * (28.10 <= pt * cosh(eta) && pt * cosh(eta) < 31.40) * (0.402053) +
    (-0.40 <= eta && eta < -0.30) * (31.40 <= pt * cosh(eta) && pt * cosh(eta) < 36.20) * (0.423693) +
    (-0.40 <= eta && eta < -0.30) * (36.20 <= pt * cosh(eta) && pt * cosh(eta) < 44.50) * (0.446017) +
    (-0.40 <= eta && eta < -0.30) * (44.50 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.460865) +
    (-0.30 <= eta && eta < -0.20) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 9.80) * (0.000624) +
    (-0.30 <= eta && eta < -0.20) * (9.80 <= pt * cosh(eta) && pt * cosh(eta) < 10.80) * (0.019209) +
    (-0.30 <= eta && eta < -0.20) * (10.80 <= pt * cosh(eta) && pt * cosh(eta) < 11.60) * (0.039736) +
    (-0.30 <= eta && eta < -0.20) * (11.60 <= pt * cosh(eta) && pt * cosh(eta) < 12.30) * (0.061691) +
    (-0.30 <= eta && eta < -0.20) * (12.30 <= pt * cosh(eta) && pt * cosh(eta) < 12.90) * (0.082927) +
    (-0.30 <= eta && eta < -0.20) * (12.90 <= pt * cosh(eta) && pt * cosh(eta) < 13.50) * (0.103434) +
    (-0.30 <= eta && eta < -0.20) * (13.50 <= pt * cosh(eta) && pt * cosh(eta) < 14.10) * (0.124136) +
    (-0.30 <= eta && eta < -0.20) * (14.10 <= pt * cosh(eta) && pt * cosh(eta) < 14.70) * (0.144554) +
    (-0.30 <= eta && eta < -0.20) * (14.70 <= pt * cosh(eta) && pt * cosh(eta) < 15.40) * (0.165953) +
    (-0.30 <= eta && eta < -0.20) * (15.40 <= pt * cosh(eta) && pt * cosh(eta) < 16.10) * (0.187883) +
    (-0.30 <= eta && eta < -0.20) * (16.10 <= pt * cosh(eta) && pt * cosh(eta) < 16.90) * (0.209869) +
    (-0.30 <= eta && eta < -0.20) * (16.90 <= pt * cosh(eta) && pt * cosh(eta) < 17.70) * (0.231539) +
    (-0.30 <= eta && eta < -0.20) * (17.70 <= pt * cosh(eta) && pt * cosh(eta) < 18.60) * (0.252531) +
    (-0.30 <= eta && eta < -0.20) * (18.60 <= pt * cosh(eta) && pt * cosh(eta) < 19.60) * (0.273661) +
    (-0.30 <= eta && eta < -0.20) * (19.60 <= pt * cosh(eta) && pt * cosh(eta) < 20.80) * (0.295298) +
    (-0.30 <= eta && eta < -0.20) * (20.80 <= pt * cosh(eta) && pt * cosh(eta) < 22.10) * (0.316682) +
    (-0.30 <= eta && eta < -0.20) * (22.10 <= pt * cosh(eta) && pt * cosh(eta) < 23.70) * (0.337790) +
    (-0.30 <= eta && eta < -0.20) * (23.70 <= pt * cosh(eta) && pt * cosh(eta) < 25.60) * (0.358983) +
    (-0.30 <= eta && eta < -0.20) * (25.60 <= pt * cosh(eta) && pt * cosh(eta) < 28.10) * (0.380366) +
    (-0.30 <= eta && eta < -0.20) * (28.10 <= pt * cosh(eta) && pt * cosh(eta) < 31.40) * (0.401997) +
    (-0.30 <= eta && eta < -0.20) * (31.40 <= pt * cosh(eta) && pt * cosh(eta) < 36.20) * (0.423648) +
    (-0.30 <= eta && eta < -0.20) * (36.20 <= pt * cosh(eta) && pt * cosh(eta) < 44.50) * (0.445986) +
    (-0.30 <= eta && eta < -0.20) * (44.50 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.460842) +
    (-0.20 <= eta && eta < -0.10) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 9.90) * (0.000667) +
    (-0.20 <= eta && eta < -0.10) * (9.90 <= pt * cosh(eta) && pt * cosh(eta) < 10.90) * (0.019658) +
    (-0.20 <= eta && eta < -0.10) * (10.90 <= pt * cosh(eta) && pt * cosh(eta) < 11.70) * (0.040176) +
    (-0.20 <= eta && eta < -0.10) * (11.70 <= pt * cosh(eta) && pt * cosh(eta) < 12.40) * (0.062033) +
    (-0.20 <= eta && eta < -0.10) * (12.40 <= pt * cosh(eta) && pt * cosh(eta) < 13.00) * (0.083131) +
    (-0.20 <= eta && eta < -0.10) * (13.00 <= pt * cosh(eta) && pt * cosh(eta) < 13.60) * (0.103490) +
    (-0.20 <= eta && eta < -0.10) * (13.60 <= pt * cosh(eta) && pt * cosh(eta) < 14.20) * (0.124039) +
    (-0.20 <= eta && eta < -0.10) * (14.20 <= pt * cosh(eta) && pt * cosh(eta) < 14.90) * (0.145966) +
    (-0.20 <= eta && eta < -0.10) * (14.90 <= pt * cosh(eta) && pt * cosh(eta) < 15.60) * (0.168751) +
    (-0.20 <= eta && eta < -0.10) * (15.60 <= pt * cosh(eta) && pt * cosh(eta) < 16.30) * (0.190365) +
    (-0.20 <= eta && eta < -0.10) * (16.30 <= pt * cosh(eta) && pt * cosh(eta) < 17.10) * (0.212019) +
    (-0.20 <= eta && eta < -0.10) * (17.10 <= pt * cosh(eta) && pt * cosh(eta) < 17.90) * (0.233359) +
    (-0.20 <= eta && eta < -0.10) * (17.90 <= pt * cosh(eta) && pt * cosh(eta) < 18.80) * (0.254035) +
    (-0.20 <= eta && eta < -0.10) * (18.80 <= pt * cosh(eta) && pt * cosh(eta) < 19.80) * (0.274858) +
    (-0.20 <= eta && eta < -0.10) * (19.80 <= pt * cosh(eta) && pt * cosh(eta) < 21.00) * (0.296199) +
    (-0.20 <= eta && eta < -0.10) * (21.00 <= pt * cosh(eta) && pt * cosh(eta) < 22.40) * (0.318075) +
    (-0.20 <= eta && eta < -0.10) * (22.40 <= pt * cosh(eta) && pt * cosh(eta) < 24.00) * (0.339496) +
    (-0.20 <= eta && eta < -0.10) * (24.00 <= pt * cosh(eta) && pt * cosh(eta) < 26.00) * (0.360760) +
    (-0.20 <= eta && eta < -0.10) * (26.00 <= pt * cosh(eta) && pt * cosh(eta) < 28.50) * (0.382069) +
    (-0.20 <= eta && eta < -0.10) * (28.50 <= pt * cosh(eta) && pt * cosh(eta) < 31.90) * (0.403435) +
    (-0.20 <= eta && eta < -0.10) * (31.90 <= pt * cosh(eta) && pt * cosh(eta) < 36.90) * (0.425149) +
    (-0.20 <= eta && eta < -0.10) * (36.90 <= pt * cosh(eta) && pt * cosh(eta) < 45.70) * (0.447621) +
    (-0.20 <= eta && eta < -0.10) * (45.70 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.461296) +
    (-0.10 <= eta && eta < -0.00) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 10.20) * (0.000679) +
    (-0.10 <= eta && eta < -0.00) * (10.20 <= pt * cosh(eta) && pt * cosh(eta) < 11.20) * (0.019921) +
    (-0.10 <= eta && eta < -0.00) * (11.20 <= pt * cosh(eta) && pt * cosh(eta) < 12.00) * (0.039923) +
    (-0.10 <= eta && eta < -0.00) * (12.00 <= pt * cosh(eta) && pt * cosh(eta) < 12.70) * (0.061117) +
    (-0.10 <= eta && eta < -0.00) * (12.70 <= pt * cosh(eta) && pt * cosh(eta) < 13.40) * (0.083207) +
    (-0.10 <= eta && eta < -0.00) * (13.40 <= pt * cosh(eta) && pt * cosh(eta) < 14.00) * (0.104676) +
    (-0.10 <= eta && eta < -0.00) * (14.00 <= pt * cosh(eta) && pt * cosh(eta) < 14.70) * (0.126328) +
    (-0.10 <= eta && eta < -0.00) * (14.70 <= pt * cosh(eta) && pt * cosh(eta) < 15.40) * (0.149235) +
    (-0.10 <= eta && eta < -0.00) * (15.40 <= pt * cosh(eta) && pt * cosh(eta) < 16.10) * (0.171280) +
    (-0.10 <= eta && eta < -0.00) * (16.10 <= pt * cosh(eta) && pt * cosh(eta) < 16.80) * (0.192190) +
    (-0.10 <= eta && eta < -0.00) * (16.80 <= pt * cosh(eta) && pt * cosh(eta) < 17.60) * (0.213160) +
    (-0.10 <= eta && eta < -0.00) * (17.60 <= pt * cosh(eta) && pt * cosh(eta) < 18.50) * (0.235081) +
    (-0.10 <= eta && eta < -0.00) * (18.50 <= pt * cosh(eta) && pt * cosh(eta) < 19.50) * (0.257304) +
    (-0.10 <= eta && eta < -0.00) * (19.50 <= pt * cosh(eta) && pt * cosh(eta) < 20.60) * (0.279240) +
    (-0.10 <= eta && eta < -0.00) * (20.60 <= pt * cosh(eta) && pt * cosh(eta) < 21.80) * (0.300402) +
    (-0.10 <= eta && eta < -0.00) * (21.80 <= pt * cosh(eta) && pt * cosh(eta) < 23.20) * (0.321151) +
    (-0.10 <= eta && eta < -0.00) * (23.20 <= pt * cosh(eta) && pt * cosh(eta) < 24.90) * (0.342173) +
    (-0.10 <= eta && eta < -0.00) * (24.90 <= pt * cosh(eta) && pt * cosh(eta) < 27.00) * (0.363463) +
    (-0.10 <= eta && eta < -0.00) * (27.00 <= pt * cosh(eta) && pt * cosh(eta) < 29.70) * (0.384883) +
    (-0.10 <= eta && eta < -0.00) * (29.70 <= pt * cosh(eta) && pt * cosh(eta) < 33.40) * (0.406518) +
    (-0.10 <= eta && eta < -0.00) * (33.40 <= pt * cosh(eta) && pt * cosh(eta) < 38.90) * (0.428380) +
    (-0.10 <= eta && eta < -0.00) * (38.90 <= pt * cosh(eta) && pt * cosh(eta) < 49.00) * (0.451065) +
    (-0.10 <= eta && eta < -0.00) * (49.00 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.461887) +
    (-0.00 <= eta && eta < 0.10) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 10.00) * (0.000652) +
    (-0.00 <= eta && eta < 0.10) * (10.00 <= pt * cosh(eta) && pt * cosh(eta) < 11.10) * (0.019853) +
    (-0.00 <= eta && eta < 0.10) * (11.10 <= pt * cosh(eta) && pt * cosh(eta) < 11.90) * (0.041350) +
    (-0.00 <= eta && eta < 0.10) * (11.90 <= pt * cosh(eta) && pt * cosh(eta) < 12.60) * (0.063083) +
    (-0.00 <= eta && eta < 0.10) * (12.60 <= pt * cosh(eta) && pt * cosh(eta) < 13.30) * (0.085621) +
    (-0.00 <= eta && eta < 0.10) * (13.30 <= pt * cosh(eta) && pt * cosh(eta) < 13.90) * (0.107436) +
    (-0.00 <= eta && eta < 0.10) * (13.90 <= pt * cosh(eta) && pt * cosh(eta) < 14.50) * (0.127685) +
    (-0.00 <= eta && eta < 0.10) * (14.50 <= pt * cosh(eta) && pt * cosh(eta) < 15.20) * (0.149221) +
    (-0.00 <= eta && eta < 0.10) * (15.20 <= pt * cosh(eta) && pt * cosh(eta) < 15.90) * (0.171557) +
    (-0.00 <= eta && eta < 0.10) * (15.90 <= pt * cosh(eta) && pt * cosh(eta) < 16.60) * (0.192724) +
    (-0.00 <= eta && eta < 0.10) * (16.60 <= pt * cosh(eta) && pt * cosh(eta) < 17.40) * (0.213928) +
    (-0.00 <= eta && eta < 0.10) * (17.40 <= pt * cosh(eta) && pt * cosh(eta) < 18.30) * (0.236064) +
    (-0.00 <= eta && eta < 0.10) * (18.30 <= pt * cosh(eta) && pt * cosh(eta) < 19.30) * (0.258471) +
    (-0.00 <= eta && eta < 0.10) * (19.30 <= pt * cosh(eta) && pt * cosh(eta) < 20.40) * (0.280550) +
    (-0.00 <= eta && eta < 0.10) * (20.40 <= pt * cosh(eta) && pt * cosh(eta) < 21.60) * (0.301813) +
    (-0.00 <= eta && eta < 0.10) * (21.60 <= pt * cosh(eta) && pt * cosh(eta) < 23.00) * (0.322623) +
    (-0.00 <= eta && eta < 0.10) * (23.00 <= pt * cosh(eta) && pt * cosh(eta) < 24.70) * (0.343665) +
    (-0.00 <= eta && eta < 0.10) * (24.70 <= pt * cosh(eta) && pt * cosh(eta) < 26.80) * (0.364932) +
    (-0.00 <= eta && eta < 0.10) * (26.80 <= pt * cosh(eta) && pt * cosh(eta) < 29.50) * (0.386281) +
    (-0.00 <= eta && eta < 0.10) * (29.50 <= pt * cosh(eta) && pt * cosh(eta) < 33.20) * (0.407794) +
    (-0.00 <= eta && eta < 0.10) * (33.20 <= pt * cosh(eta) && pt * cosh(eta) < 38.80) * (0.429658) +
    (-0.00 <= eta && eta < 0.10) * (38.80 <= pt * cosh(eta) && pt * cosh(eta) < 49.20) * (0.452424) +
    (-0.00 <= eta && eta < 0.10) * (49.20 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.463042) +
    (0.10 <= eta && eta < 0.20) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 9.80) * (0.000688) +
    (0.10 <= eta && eta < 0.20) * (9.80 <= pt * cosh(eta) && pt * cosh(eta) < 10.80) * (0.020283) +
    (0.10 <= eta && eta < 0.20) * (10.80 <= pt * cosh(eta) && pt * cosh(eta) < 11.60) * (0.041390) +
    (0.10 <= eta && eta < 0.20) * (11.60 <= pt * cosh(eta) && pt * cosh(eta) < 12.30) * (0.063752) +
    (0.10 <= eta && eta < 0.20) * (12.30 <= pt * cosh(eta) && pt * cosh(eta) < 12.90) * (0.085249) +
    (0.10 <= eta && eta < 0.20) * (12.90 <= pt * cosh(eta) && pt * cosh(eta) < 13.50) * (0.105921) +
    (0.10 <= eta && eta < 0.20) * (13.50 <= pt * cosh(eta) && pt * cosh(eta) < 14.10) * (0.126724) +
    (0.10 <= eta && eta < 0.20) * (14.10 <= pt * cosh(eta) && pt * cosh(eta) < 14.70) * (0.147189) +
    (0.10 <= eta && eta < 0.20) * (14.70 <= pt * cosh(eta) && pt * cosh(eta) < 15.40) * (0.168592) +
    (0.10 <= eta && eta < 0.20) * (15.40 <= pt * cosh(eta) && pt * cosh(eta) < 16.10) * (0.190487) +
    (0.10 <= eta && eta < 0.20) * (16.10 <= pt * cosh(eta) && pt * cosh(eta) < 16.90) * (0.212403) +
    (0.10 <= eta && eta < 0.20) * (16.90 <= pt * cosh(eta) && pt * cosh(eta) < 17.70) * (0.233976) +
    (0.10 <= eta && eta < 0.20) * (17.70 <= pt * cosh(eta) && pt * cosh(eta) < 18.60) * (0.254850) +
    (0.10 <= eta && eta < 0.20) * (18.60 <= pt * cosh(eta) && pt * cosh(eta) < 19.60) * (0.275841) +
    (0.10 <= eta && eta < 0.20) * (19.60 <= pt * cosh(eta) && pt * cosh(eta) < 20.80) * (0.297319) +
    (0.10 <= eta && eta < 0.20) * (20.80 <= pt * cosh(eta) && pt * cosh(eta) < 22.20) * (0.319295) +
    (0.10 <= eta && eta < 0.20) * (22.20 <= pt * cosh(eta) && pt * cosh(eta) < 23.80) * (0.340773) +
    (0.10 <= eta && eta < 0.20) * (23.80 <= pt * cosh(eta) && pt * cosh(eta) < 25.80) * (0.362051) +
    (0.10 <= eta && eta < 0.20) * (25.80 <= pt * cosh(eta) && pt * cosh(eta) < 28.30) * (0.383326) +
    (0.10 <= eta && eta < 0.20) * (28.30 <= pt * cosh(eta) && pt * cosh(eta) < 31.70) * (0.404609) +
    (0.10 <= eta && eta < 0.20) * (31.70 <= pt * cosh(eta) && pt * cosh(eta) < 36.80) * (0.426384) +
    (0.10 <= eta && eta < 0.20) * (36.80 <= pt * cosh(eta) && pt * cosh(eta) < 45.90) * (0.449033) +
    (0.10 <= eta && eta < 0.20) * (45.90 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.462443) +
    (0.20 <= eta && eta < 0.30) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 9.70) * (0.000702) +
    (0.20 <= eta && eta < 0.30) * (9.70 <= pt * cosh(eta) && pt * cosh(eta) < 10.70) * (0.021415) +
    (0.20 <= eta && eta < 0.30) * (10.70 <= pt * cosh(eta) && pt * cosh(eta) < 11.40) * (0.041984) +
    (0.20 <= eta && eta < 0.30) * (11.40 <= pt * cosh(eta) && pt * cosh(eta) < 12.10) * (0.063242) +
    (0.20 <= eta && eta < 0.30) * (12.10 <= pt * cosh(eta) && pt * cosh(eta) < 12.70) * (0.085050) +
    (0.20 <= eta && eta < 0.30) * (12.70 <= pt * cosh(eta) && pt * cosh(eta) < 13.30) * (0.106045) +
    (0.20 <= eta && eta < 0.30) * (13.30 <= pt * cosh(eta) && pt * cosh(eta) < 13.90) * (0.127172) +
    (0.20 <= eta && eta < 0.30) * (13.90 <= pt * cosh(eta) && pt * cosh(eta) < 14.50) * (0.147941) +
    (0.20 <= eta && eta < 0.30) * (14.50 <= pt * cosh(eta) && pt * cosh(eta) < 15.20) * (0.169640) +
    (0.20 <= eta && eta < 0.30) * (15.20 <= pt * cosh(eta) && pt * cosh(eta) < 15.90) * (0.191806) +
    (0.20 <= eta && eta < 0.30) * (15.90 <= pt * cosh(eta) && pt * cosh(eta) < 16.70) * (0.213956) +
    (0.20 <= eta && eta < 0.30) * (16.70 <= pt * cosh(eta) && pt * cosh(eta) < 17.50) * (0.235721) +
    (0.20 <= eta && eta < 0.30) * (17.50 <= pt * cosh(eta) && pt * cosh(eta) < 18.40) * (0.256741) +
    (0.20 <= eta && eta < 0.30) * (18.40 <= pt * cosh(eta) && pt * cosh(eta) < 19.40) * (0.277837) +
    (0.20 <= eta && eta < 0.30) * (19.40 <= pt * cosh(eta) && pt * cosh(eta) < 20.60) * (0.299379) +
    (0.20 <= eta && eta < 0.30) * (20.60 <= pt * cosh(eta) && pt * cosh(eta) < 22.00) * (0.321372) +
    (0.20 <= eta && eta < 0.30) * (22.00 <= pt * cosh(eta) && pt * cosh(eta) < 23.60) * (0.342817) +
    (0.20 <= eta && eta < 0.30) * (23.60 <= pt * cosh(eta) && pt * cosh(eta) < 25.60) * (0.364014) +
    (0.20 <= eta && eta < 0.30) * (25.60 <= pt * cosh(eta) && pt * cosh(eta) < 28.20) * (0.385555) +
    (0.20 <= eta && eta < 0.30) * (28.20 <= pt * cosh(eta) && pt * cosh(eta) < 31.70) * (0.407162) +
    (0.20 <= eta && eta < 0.30) * (31.70 <= pt * cosh(eta) && pt * cosh(eta) < 37.00) * (0.429002) +
    (0.20 <= eta && eta < 0.30) * (37.00 <= pt * cosh(eta) && pt * cosh(eta) < 46.80) * (0.451798) +
    (0.20 <= eta && eta < 0.30) * (46.80 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.464291) +
    (0.30 <= eta && eta < 0.40) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 9.70) * (0.000726) +
    (0.30 <= eta && eta < 0.40) * (9.70 <= pt * cosh(eta) && pt * cosh(eta) < 10.70) * (0.020369) +
    (0.30 <= eta && eta < 0.40) * (10.70 <= pt * cosh(eta) && pt * cosh(eta) < 11.50) * (0.041656) +
    (0.30 <= eta && eta < 0.40) * (11.50 <= pt * cosh(eta) && pt * cosh(eta) < 12.20) * (0.064296) +
    (0.30 <= eta && eta < 0.40) * (12.20 <= pt * cosh(eta) && pt * cosh(eta) < 12.80) * (0.086047) +
    (0.30 <= eta && eta < 0.40) * (12.80 <= pt * cosh(eta) && pt * cosh(eta) < 13.40) * (0.106941) +
    (0.30 <= eta && eta < 0.40) * (13.40 <= pt * cosh(eta) && pt * cosh(eta) < 14.00) * (0.127942) +
    (0.30 <= eta && eta < 0.40) * (14.00 <= pt * cosh(eta) && pt * cosh(eta) < 14.60) * (0.148575) +
    (0.30 <= eta && eta < 0.40) * (14.60 <= pt * cosh(eta) && pt * cosh(eta) < 15.30) * (0.170125) +
    (0.30 <= eta && eta < 0.40) * (15.30 <= pt * cosh(eta) && pt * cosh(eta) < 16.00) * (0.192139) +
    (0.30 <= eta && eta < 0.40) * (16.00 <= pt * cosh(eta) && pt * cosh(eta) < 16.80) * (0.214143) +
    (0.30 <= eta && eta < 0.40) * (16.80 <= pt * cosh(eta) && pt * cosh(eta) < 17.60) * (0.235774) +
    (0.30 <= eta && eta < 0.40) * (17.60 <= pt * cosh(eta) && pt * cosh(eta) < 18.50) * (0.256674) +
    (0.30 <= eta && eta < 0.40) * (18.50 <= pt * cosh(eta) && pt * cosh(eta) < 19.50) * (0.277665) +
    (0.30 <= eta && eta < 0.40) * (19.50 <= pt * cosh(eta) && pt * cosh(eta) < 20.70) * (0.299113) +
    (0.30 <= eta && eta < 0.40) * (20.70 <= pt * cosh(eta) && pt * cosh(eta) < 22.10) * (0.321029) +
    (0.30 <= eta && eta < 0.40) * (22.10 <= pt * cosh(eta) && pt * cosh(eta) < 23.70) * (0.342418) +
    (0.30 <= eta && eta < 0.40) * (23.70 <= pt * cosh(eta) && pt * cosh(eta) < 25.70) * (0.363580) +
    (0.30 <= eta && eta < 0.40) * (25.70 <= pt * cosh(eta) && pt * cosh(eta) < 28.30) * (0.385108) +
    (0.30 <= eta && eta < 0.40) * (28.30 <= pt * cosh(eta) && pt * cosh(eta) < 31.80) * (0.406727) +
    (0.30 <= eta && eta < 0.40) * (31.80 <= pt * cosh(eta) && pt * cosh(eta) < 37.10) * (0.428608) +
    (0.30 <= eta && eta < 0.40) * (37.10 <= pt * cosh(eta) && pt * cosh(eta) < 46.80) * (0.451376) +
    (0.30 <= eta && eta < 0.40) * (46.80 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.463881) +
    (0.40 <= eta && eta < 0.50) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 9.60) * (0.000671) +
    (0.40 <= eta && eta < 0.50) * (9.60 <= pt * cosh(eta) && pt * cosh(eta) < 10.60) * (0.020318) +
    (0.40 <= eta && eta < 0.50) * (10.60 <= pt * cosh(eta) && pt * cosh(eta) < 11.30) * (0.040535) +
    (0.40 <= eta && eta < 0.50) * (11.30 <= pt * cosh(eta) && pt * cosh(eta) < 12.00) * (0.061663) +
    (0.40 <= eta && eta < 0.50) * (12.00 <= pt * cosh(eta) && pt * cosh(eta) < 12.60) * (0.083458) +
    (0.40 <= eta && eta < 0.50) * (12.60 <= pt * cosh(eta) && pt * cosh(eta) < 13.20) * (0.104511) +
    (0.40 <= eta && eta < 0.50) * (13.20 <= pt * cosh(eta) && pt * cosh(eta) < 13.80) * (0.125739) +
    (0.40 <= eta && eta < 0.50) * (13.80 <= pt * cosh(eta) && pt * cosh(eta) < 14.40) * (0.146635) +
    (0.40 <= eta && eta < 0.50) * (14.40 <= pt * cosh(eta) && pt * cosh(eta) < 15.10) * (0.168484) +
    (0.40 <= eta && eta < 0.50) * (15.10 <= pt * cosh(eta) && pt * cosh(eta) < 15.80) * (0.190813) +
    (0.40 <= eta && eta < 0.50) * (15.80 <= pt * cosh(eta) && pt * cosh(eta) < 16.60) * (0.213130) +
    (0.40 <= eta && eta < 0.50) * (16.60 <= pt * cosh(eta) && pt * cosh(eta) < 17.40) * (0.235057) +
    (0.40 <= eta && eta < 0.50) * (17.40 <= pt * cosh(eta) && pt * cosh(eta) < 18.30) * (0.256227) +
    (0.40 <= eta && eta < 0.50) * (18.30 <= pt * cosh(eta) && pt * cosh(eta) < 19.30) * (0.277468) +
    (0.40 <= eta && eta < 0.50) * (19.30 <= pt * cosh(eta) && pt * cosh(eta) < 20.50) * (0.299144) +
    (0.40 <= eta && eta < 0.50) * (20.50 <= pt * cosh(eta) && pt * cosh(eta) < 21.90) * (0.321260) +
    (0.40 <= eta && eta < 0.50) * (21.90 <= pt * cosh(eta) && pt * cosh(eta) < 23.50) * (0.342810) +
    (0.40 <= eta && eta < 0.50) * (23.50 <= pt * cosh(eta) && pt * cosh(eta) < 25.50) * (0.364092) +
    (0.40 <= eta && eta < 0.50) * (25.50 <= pt * cosh(eta) && pt * cosh(eta) < 28.10) * (0.385699) +
    (0.40 <= eta && eta < 0.50) * (28.10 <= pt * cosh(eta) && pt * cosh(eta) < 31.60) * (0.407349) +
    (0.40 <= eta && eta < 0.50) * (31.60 <= pt * cosh(eta) && pt * cosh(eta) < 36.90) * (0.429205) +
    (0.40 <= eta && eta < 0.50) * (36.90 <= pt * cosh(eta) && pt * cosh(eta) < 46.70) * (0.451986) +
    (0.40 <= eta && eta < 0.50) * (46.70 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.464527) +
    (0.50 <= eta && eta < 0.60) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 9.70) * (0.000629) +
    (0.50 <= eta && eta < 0.60) * (9.70 <= pt * cosh(eta) && pt * cosh(eta) < 10.70) * (0.019027) +
    (0.50 <= eta && eta < 0.60) * (10.70 <= pt * cosh(eta) && pt * cosh(eta) < 11.50) * (0.039568) +
    (0.50 <= eta && eta < 0.60) * (11.50 <= pt * cosh(eta) && pt * cosh(eta) < 12.20) * (0.061693) +
    (0.50 <= eta && eta < 0.60) * (12.20 <= pt * cosh(eta) && pt * cosh(eta) < 12.80) * (0.083113) +
    (0.50 <= eta && eta < 0.60) * (12.80 <= pt * cosh(eta) && pt * cosh(eta) < 13.40) * (0.103800) +
    (0.50 <= eta && eta < 0.60) * (13.40 <= pt * cosh(eta) && pt * cosh(eta) < 14.00) * (0.124675) +
    (0.50 <= eta && eta < 0.60) * (14.00 <= pt * cosh(eta) && pt * cosh(eta) < 14.60) * (0.145250) +
    (0.50 <= eta && eta < 0.60) * (14.60 <= pt * cosh(eta) && pt * cosh(eta) < 15.30) * (0.166798) +
    (0.50 <= eta && eta < 0.60) * (15.30 <= pt * cosh(eta) && pt * cosh(eta) < 16.00) * (0.188860) +
    (0.50 <= eta && eta < 0.60) * (16.00 <= pt * cosh(eta) && pt * cosh(eta) < 16.80) * (0.210955) +
    (0.50 <= eta && eta < 0.60) * (16.80 <= pt * cosh(eta) && pt * cosh(eta) < 17.60) * (0.232711) +
    (0.50 <= eta && eta < 0.60) * (17.60 <= pt * cosh(eta) && pt * cosh(eta) < 18.50) * (0.253762) +
    (0.50 <= eta && eta < 0.60) * (18.50 <= pt * cosh(eta) && pt * cosh(eta) < 19.50) * (0.274929) +
    (0.50 <= eta && eta < 0.60) * (19.50 <= pt * cosh(eta) && pt * cosh(eta) < 20.70) * (0.296580) +
    (0.50 <= eta && eta < 0.60) * (20.70 <= pt * cosh(eta) && pt * cosh(eta) < 22.00) * (0.317952) +
    (0.50 <= eta && eta < 0.60) * (22.00 <= pt * cosh(eta) && pt * cosh(eta) < 23.60) * (0.339024) +
    (0.50 <= eta && eta < 0.60) * (23.60 <= pt * cosh(eta) && pt * cosh(eta) < 25.50) * (0.360155) +
    (0.50 <= eta && eta < 0.60) * (25.50 <= pt * cosh(eta) && pt * cosh(eta) < 28.00) * (0.381450) +
    (0.50 <= eta && eta < 0.60) * (28.00 <= pt * cosh(eta) && pt * cosh(eta) < 31.30) * (0.402963) +
    (0.50 <= eta && eta < 0.60) * (31.30 <= pt * cosh(eta) && pt * cosh(eta) < 36.20) * (0.424674) +
    (0.50 <= eta && eta < 0.60) * (36.20 <= pt * cosh(eta) && pt * cosh(eta) < 44.70) * (0.447121) +
    (0.50 <= eta && eta < 0.60) * (44.70 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.461667) +
    (0.60 <= eta && eta < 0.70) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 9.80) * (0.000675) +
    (0.60 <= eta && eta < 0.70) * (9.80 <= pt * cosh(eta) && pt * cosh(eta) < 10.80) * (0.019837) +
    (0.60 <= eta && eta < 0.70) * (10.80 <= pt * cosh(eta) && pt * cosh(eta) < 11.60) * (0.040690) +
    (0.60 <= eta && eta < 0.70) * (11.60 <= pt * cosh(eta) && pt * cosh(eta) < 12.30) * (0.062882) +
    (0.60 <= eta && eta < 0.70) * (12.30 <= pt * cosh(eta) && pt * cosh(eta) < 12.90) * (0.084270) +
    (0.60 <= eta && eta < 0.70) * (12.90 <= pt * cosh(eta) && pt * cosh(eta) < 13.50) * (0.104874) +
    (0.60 <= eta && eta < 0.70) * (13.50 <= pt * cosh(eta) && pt * cosh(eta) < 14.10) * (0.125635) +
    (0.60 <= eta && eta < 0.70) * (14.10 <= pt * cosh(eta) && pt * cosh(eta) < 14.70) * (0.146081) +
    (0.60 <= eta && eta < 0.70) * (14.70 <= pt * cosh(eta) && pt * cosh(eta) < 15.40) * (0.167483) +
    (0.60 <= eta && eta < 0.70) * (15.40 <= pt * cosh(eta) && pt * cosh(eta) < 16.10) * (0.189393) +
    (0.60 <= eta && eta < 0.70) * (16.10 <= pt * cosh(eta) && pt * cosh(eta) < 16.90) * (0.211339) +
    (0.60 <= eta && eta < 0.70) * (16.90 <= pt * cosh(eta) && pt * cosh(eta) < 17.70) * (0.232954) +
    (0.60 <= eta && eta < 0.70) * (17.70 <= pt * cosh(eta) && pt * cosh(eta) < 18.60) * (0.253877) +
    (0.60 <= eta && eta < 0.70) * (18.60 <= pt * cosh(eta) && pt * cosh(eta) < 19.60) * (0.274927) +
    (0.60 <= eta && eta < 0.70) * (19.60 <= pt * cosh(eta) && pt * cosh(eta) < 20.80) * (0.296472) +
    (0.60 <= eta && eta < 0.70) * (20.80 <= pt * cosh(eta) && pt * cosh(eta) < 22.20) * (0.318523) +
    (0.60 <= eta && eta < 0.70) * (22.20 <= pt * cosh(eta) && pt * cosh(eta) < 23.80) * (0.340081) +
    (0.60 <= eta && eta < 0.70) * (23.80 <= pt * cosh(eta) && pt * cosh(eta) < 25.80) * (0.361442) +
    (0.60 <= eta && eta < 0.70) * (25.80 <= pt * cosh(eta) && pt * cosh(eta) < 28.30) * (0.382805) +
    (0.60 <= eta && eta < 0.70) * (28.30 <= pt * cosh(eta) && pt * cosh(eta) < 31.70) * (0.404178) +
    (0.60 <= eta && eta < 0.70) * (31.70 <= pt * cosh(eta) && pt * cosh(eta) < 36.80) * (0.426049) +
    (0.60 <= eta && eta < 0.70) * (36.80 <= pt * cosh(eta) && pt * cosh(eta) < 45.80) * (0.448688) +
    (0.60 <= eta && eta < 0.70) * (45.80 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.462189) +
    (0.70 <= eta && eta < 0.80) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 9.80) * (0.000676) +
    (0.70 <= eta && eta < 0.80) * (9.80 <= pt * cosh(eta) && pt * cosh(eta) < 10.80) * (0.019969) +
    (0.70 <= eta && eta < 0.80) * (10.80 <= pt * cosh(eta) && pt * cosh(eta) < 11.60) * (0.040896) +
    (0.70 <= eta && eta < 0.80) * (11.60 <= pt * cosh(eta) && pt * cosh(eta) < 12.30) * (0.063138) +
    (0.70 <= eta && eta < 0.80) * (12.30 <= pt * cosh(eta) && pt * cosh(eta) < 12.90) * (0.084558) +
    (0.70 <= eta && eta < 0.80) * (12.90 <= pt * cosh(eta) && pt * cosh(eta) < 13.50) * (0.105182) +
    (0.70 <= eta && eta < 0.80) * (13.50 <= pt * cosh(eta) && pt * cosh(eta) < 14.10) * (0.125956) +
    (0.70 <= eta && eta < 0.80) * (14.10 <= pt * cosh(eta) && pt * cosh(eta) < 14.70) * (0.146407) +
    (0.70 <= eta && eta < 0.80) * (14.70 <= pt * cosh(eta) && pt * cosh(eta) < 15.40) * (0.167810) +
    (0.70 <= eta && eta < 0.80) * (15.40 <= pt * cosh(eta) && pt * cosh(eta) < 16.10) * (0.189716) +
    (0.70 <= eta && eta < 0.80) * (16.10 <= pt * cosh(eta) && pt * cosh(eta) < 16.90) * (0.211653) +
    (0.70 <= eta && eta < 0.80) * (16.90 <= pt * cosh(eta) && pt * cosh(eta) < 17.70) * (0.233255) +
    (0.70 <= eta && eta < 0.80) * (17.70 <= pt * cosh(eta) && pt * cosh(eta) < 18.60) * (0.254164) +
    (0.70 <= eta && eta < 0.80) * (18.60 <= pt * cosh(eta) && pt * cosh(eta) < 19.60) * (0.275196) +
    (0.70 <= eta && eta < 0.80) * (19.60 <= pt * cosh(eta) && pt * cosh(eta) < 20.80) * (0.296722) +
    (0.70 <= eta && eta < 0.80) * (20.80 <= pt * cosh(eta) && pt * cosh(eta) < 22.20) * (0.318751) +
    (0.70 <= eta && eta < 0.80) * (22.20 <= pt * cosh(eta) && pt * cosh(eta) < 23.80) * (0.340285) +
    (0.70 <= eta && eta < 0.80) * (23.80 <= pt * cosh(eta) && pt * cosh(eta) < 25.80) * (0.361622) +
    (0.70 <= eta && eta < 0.80) * (25.80 <= pt * cosh(eta) && pt * cosh(eta) < 28.30) * (0.382958) +
    (0.70 <= eta && eta < 0.80) * (28.30 <= pt * cosh(eta) && pt * cosh(eta) < 31.70) * (0.404306) +
    (0.70 <= eta && eta < 0.80) * (31.70 <= pt * cosh(eta) && pt * cosh(eta) < 36.80) * (0.426148) +
    (0.70 <= eta && eta < 0.80) * (36.80 <= pt * cosh(eta) && pt * cosh(eta) < 45.80) * (0.448757) +
    (0.70 <= eta && eta < 0.80) * (45.80 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.462240) +
    (0.80 <= eta && eta < 0.90) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 10.00) * (0.000666) +
    (0.80 <= eta && eta < 0.90) * (10.00 <= pt * cosh(eta) && pt * cosh(eta) < 11.10) * (0.019973) +
    (0.80 <= eta && eta < 0.90) * (11.10 <= pt * cosh(eta) && pt * cosh(eta) < 11.90) * (0.041538) +
    (0.80 <= eta && eta < 0.90) * (11.90 <= pt * cosh(eta) && pt * cosh(eta) < 12.60) * (0.063316) +
    (0.80 <= eta && eta < 0.90) * (12.60 <= pt * cosh(eta) && pt * cosh(eta) < 13.30) * (0.085883) +
    (0.80 <= eta && eta < 0.90) * (13.30 <= pt * cosh(eta) && pt * cosh(eta) < 13.90) * (0.107717) +
    (0.80 <= eta && eta < 0.90) * (13.90 <= pt * cosh(eta) && pt * cosh(eta) < 14.50) * (0.127976) +
    (0.80 <= eta && eta < 0.90) * (14.50 <= pt * cosh(eta) && pt * cosh(eta) < 15.20) * (0.149517) +
    (0.80 <= eta && eta < 0.90) * (15.20 <= pt * cosh(eta) && pt * cosh(eta) < 15.90) * (0.171852) +
    (0.80 <= eta && eta < 0.90) * (15.90 <= pt * cosh(eta) && pt * cosh(eta) < 16.60) * (0.193015) +
    (0.80 <= eta && eta < 0.90) * (16.60 <= pt * cosh(eta) && pt * cosh(eta) < 17.40) * (0.214211) +
    (0.80 <= eta && eta < 0.90) * (17.40 <= pt * cosh(eta) && pt * cosh(eta) < 18.30) * (0.236336) +
    (0.80 <= eta && eta < 0.90) * (18.30 <= pt * cosh(eta) && pt * cosh(eta) < 19.30) * (0.258728) +
    (0.80 <= eta && eta < 0.90) * (19.30 <= pt * cosh(eta) && pt * cosh(eta) < 20.40) * (0.280790) +
    (0.80 <= eta && eta < 0.90) * (20.40 <= pt * cosh(eta) && pt * cosh(eta) < 21.60) * (0.302035) +
    (0.80 <= eta && eta < 0.90) * (21.60 <= pt * cosh(eta) && pt * cosh(eta) < 23.00) * (0.322825) +
    (0.80 <= eta && eta < 0.90) * (23.00 <= pt * cosh(eta) && pt * cosh(eta) < 24.70) * (0.343847) +
    (0.80 <= eta && eta < 0.90) * (24.70 <= pt * cosh(eta) && pt * cosh(eta) < 26.80) * (0.365091) +
    (0.80 <= eta && eta < 0.90) * (26.80 <= pt * cosh(eta) && pt * cosh(eta) < 29.50) * (0.386417) +
    (0.80 <= eta && eta < 0.90) * (29.50 <= pt * cosh(eta) && pt * cosh(eta) < 33.20) * (0.407905) +
    (0.80 <= eta && eta < 0.90) * (33.20 <= pt * cosh(eta) && pt * cosh(eta) < 38.80) * (0.429743) +
    (0.80 <= eta && eta < 0.90) * (38.80 <= pt * cosh(eta) && pt * cosh(eta) < 49.20) * (0.452482) +
    (0.80 <= eta && eta < 0.90) * (49.20 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.463087) +
    (0.90 <= eta && eta < 1.00) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 10.00) * (0.000631) +
    (0.90 <= eta && eta < 1.00) * (10.00 <= pt * cosh(eta) && pt * cosh(eta) < 11.10) * (0.020111) +
    (0.90 <= eta && eta < 1.00) * (11.10 <= pt * cosh(eta) && pt * cosh(eta) < 11.90) * (0.041754) +
    (0.90 <= eta && eta < 1.00) * (11.90 <= pt * cosh(eta) && pt * cosh(eta) < 12.60) * (0.063581) +
    (0.90 <= eta && eta < 1.00) * (12.60 <= pt * cosh(eta) && pt * cosh(eta) < 13.20) * (0.084516) +
    (0.90 <= eta && eta < 1.00) * (13.20 <= pt * cosh(eta) && pt * cosh(eta) < 13.80) * (0.104653) +
    (0.90 <= eta && eta < 1.00) * (13.80 <= pt * cosh(eta) && pt * cosh(eta) < 14.40) * (0.124945) +
    (0.90 <= eta && eta < 1.00) * (14.40 <= pt * cosh(eta) && pt * cosh(eta) < 15.10) * (0.146583) +
    (0.90 <= eta && eta < 1.00) * (15.10 <= pt * cosh(eta) && pt * cosh(eta) < 15.80) * (0.169067) +
    (0.90 <= eta && eta < 1.00) * (15.80 <= pt * cosh(eta) && pt * cosh(eta) < 16.50) * (0.190404) +
    (0.90 <= eta && eta < 1.00) * (16.50 <= pt * cosh(eta) && pt * cosh(eta) < 17.30) * (0.211798) +
    (0.90 <= eta && eta < 1.00) * (17.30 <= pt * cosh(eta) && pt * cosh(eta) < 18.20) * (0.234147) +
    (0.90 <= eta && eta < 1.00) * (18.20 <= pt * cosh(eta) && pt * cosh(eta) < 19.10) * (0.255662) +
    (0.90 <= eta && eta < 1.00) * (19.10 <= pt * cosh(eta) && pt * cosh(eta) < 20.20) * (0.277068) +
    (0.90 <= eta && eta < 1.00) * (20.20 <= pt * cosh(eta) && pt * cosh(eta) < 21.40) * (0.298797) +
    (0.90 <= eta && eta < 1.00) * (21.40 <= pt * cosh(eta) && pt * cosh(eta) < 22.80) * (0.320057) +
    (0.90 <= eta && eta < 1.00) * (22.80 <= pt * cosh(eta) && pt * cosh(eta) < 24.50) * (0.341541) +
    (0.90 <= eta && eta < 1.00) * (24.50 <= pt * cosh(eta) && pt * cosh(eta) < 26.60) * (0.363233) +
    (0.90 <= eta && eta < 1.00) * (26.60 <= pt * cosh(eta) && pt * cosh(eta) < 29.20) * (0.384597) +
    (0.90 <= eta && eta < 1.00) * (29.20 <= pt * cosh(eta) && pt * cosh(eta) < 32.80) * (0.405989) +
    (0.90 <= eta && eta < 1.00) * (32.80 <= pt * cosh(eta) && pt * cosh(eta) < 38.20) * (0.427888) +
    (0.90 <= eta && eta < 1.00) * (38.20 <= pt * cosh(eta) && pt * cosh(eta) < 48.00) * (0.450605) +
    (0.90 <= eta && eta < 1.00) * (48.00 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.462220) +
    (1.00 <= eta && eta < 1.10) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.000000) }

  # --- pions ---

  add EfficiencyFormula {-211} {321} {
    (eta<-1.00 || eta>=1.00 || pt * cosh(eta) < 0.10 || pt * cosh(eta) >= 50.00) * ( 0.00 ) +
    (-1.00 <= eta && eta < -0.90) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 6.10) * (0.000629) +
    (-1.00 <= eta && eta < -0.90) * (6.10 <= pt * cosh(eta) && pt * cosh(eta) < 6.80) * (0.020818) +
    (-1.00 <= eta && eta < -0.90) * (6.80 <= pt * cosh(eta) && pt * cosh(eta) < 7.30) * (0.043238) +
    (-1.00 <= eta && eta < -0.90) * (7.30 <= pt * cosh(eta) && pt * cosh(eta) < 7.70) * (0.064987) +
    (-1.00 <= eta && eta < -0.90) * (7.70 <= pt * cosh(eta) && pt * cosh(eta) < 8.10) * (0.086643) +
    (-1.00 <= eta && eta < -0.90) * (8.10 <= pt * cosh(eta) && pt * cosh(eta) < 8.50) * (0.108246) +
    (-1.00 <= eta && eta < -0.90) * (8.50 <= pt * cosh(eta) && pt * cosh(eta) < 8.90) * (0.131038) +
    (-1.00 <= eta && eta < -0.90) * (8.90 <= pt * cosh(eta) && pt * cosh(eta) < 9.30) * (0.153265) +
    (-1.00 <= eta && eta < -0.90) * (9.30 <= pt * cosh(eta) && pt * cosh(eta) < 9.70) * (0.173948) +
    (-1.00 <= eta && eta < -0.90) * (9.70 <= pt * cosh(eta) && pt * cosh(eta) < 10.20) * (0.196026) +
    (-1.00 <= eta && eta < -0.90) * (10.20 <= pt * cosh(eta) && pt * cosh(eta) < 10.70) * (0.219035) +
    (-1.00 <= eta && eta < -0.90) * (10.70 <= pt * cosh(eta) && pt * cosh(eta) < 11.30) * (0.242086) +
    (-1.00 <= eta && eta < -0.90) * (11.30 <= pt * cosh(eta) && pt * cosh(eta) < 11.90) * (0.264737) +
    (-1.00 <= eta && eta < -0.90) * (11.90 <= pt * cosh(eta) && pt * cosh(eta) < 12.60) * (0.286490) +
    (-1.00 <= eta && eta < -0.90) * (12.60 <= pt * cosh(eta) && pt * cosh(eta) < 13.40) * (0.308418) +
    (-1.00 <= eta && eta < -0.90) * (13.40 <= pt * cosh(eta) && pt * cosh(eta) < 14.30) * (0.329730) +
    (-1.00 <= eta && eta < -0.90) * (14.30 <= pt * cosh(eta) && pt * cosh(eta) < 15.40) * (0.350771) +
    (-1.00 <= eta && eta < -0.90) * (15.40 <= pt * cosh(eta) && pt * cosh(eta) < 16.80) * (0.372200) +
    (-1.00 <= eta && eta < -0.90) * (16.80 <= pt * cosh(eta) && pt * cosh(eta) < 18.60) * (0.393662) +
    (-1.00 <= eta && eta < -0.90) * (18.60 <= pt * cosh(eta) && pt * cosh(eta) < 21.20) * (0.415409) +
    (-1.00 <= eta && eta < -0.90) * (21.20 <= pt * cosh(eta) && pt * cosh(eta) < 25.30) * (0.437641) +
    (-1.00 <= eta && eta < -0.90) * (25.30 <= pt * cosh(eta) && pt * cosh(eta) < 34.20) * (0.461287) +
    (-1.00 <= eta && eta < -0.90) * (34.20 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.480405) +
    (-0.90 <= eta && eta < -0.80) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 6.00) * (0.000601) +
    (-0.90 <= eta && eta < -0.80) * (6.00 <= pt * cosh(eta) && pt * cosh(eta) < 6.70) * (0.020257) +
    (-0.90 <= eta && eta < -0.80) * (6.70 <= pt * cosh(eta) && pt * cosh(eta) < 7.20) * (0.043042) +
    (-0.90 <= eta && eta < -0.80) * (7.20 <= pt * cosh(eta) && pt * cosh(eta) < 7.60) * (0.065328) +
    (-0.90 <= eta && eta < -0.80) * (7.60 <= pt * cosh(eta) && pt * cosh(eta) < 8.00) * (0.087115) +
    (-0.90 <= eta && eta < -0.80) * (8.00 <= pt * cosh(eta) && pt * cosh(eta) < 8.40) * (0.109928) +
    (-0.90 <= eta && eta < -0.80) * (8.40 <= pt * cosh(eta) && pt * cosh(eta) < 8.80) * (0.132996) +
    (-0.90 <= eta && eta < -0.80) * (8.80 <= pt * cosh(eta) && pt * cosh(eta) < 9.30) * (0.156649) +
    (-0.90 <= eta && eta < -0.80) * (9.30 <= pt * cosh(eta) && pt * cosh(eta) < 9.70) * (0.180600) +
    (-0.90 <= eta && eta < -0.80) * (9.70 <= pt * cosh(eta) && pt * cosh(eta) < 10.20) * (0.203148) +
    (-0.90 <= eta && eta < -0.80) * (10.20 <= pt * cosh(eta) && pt * cosh(eta) < 10.70) * (0.225705) +
    (-0.90 <= eta && eta < -0.80) * (10.70 <= pt * cosh(eta) && pt * cosh(eta) < 11.30) * (0.248439) +
    (-0.90 <= eta && eta < -0.80) * (11.30 <= pt * cosh(eta) && pt * cosh(eta) < 11.90) * (0.270714) +
    (-0.90 <= eta && eta < -0.80) * (11.90 <= pt * cosh(eta) && pt * cosh(eta) < 12.60) * (0.292054) +
    (-0.90 <= eta && eta < -0.80) * (12.60 <= pt * cosh(eta) && pt * cosh(eta) < 13.40) * (0.313521) +
    (-0.90 <= eta && eta < -0.80) * (13.40 <= pt * cosh(eta) && pt * cosh(eta) < 14.40) * (0.335443) +
    (-0.90 <= eta && eta < -0.80) * (14.40 <= pt * cosh(eta) && pt * cosh(eta) < 15.60) * (0.357624) +
    (-0.90 <= eta && eta < -0.80) * (15.60 <= pt * cosh(eta) && pt * cosh(eta) < 17.10) * (0.379405) +
    (-0.90 <= eta && eta < -0.80) * (17.10 <= pt * cosh(eta) && pt * cosh(eta) < 19.10) * (0.401051) +
    (-0.90 <= eta && eta < -0.80) * (19.10 <= pt * cosh(eta) && pt * cosh(eta) < 22.00) * (0.422825) +
    (-0.90 <= eta && eta < -0.80) * (22.00 <= pt * cosh(eta) && pt * cosh(eta) < 27.00) * (0.445280) +
    (-0.90 <= eta && eta < -0.80) * (27.00 <= pt * cosh(eta) && pt * cosh(eta) < 41.00) * (0.470610) +
    (-0.90 <= eta && eta < -0.80) * (41.00 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.484124) +
    (-0.80 <= eta && eta < -0.70) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 6.00) * (0.000757) +
    (-0.80 <= eta && eta < -0.70) * (6.00 <= pt * cosh(eta) && pt * cosh(eta) < 6.60) * (0.022320) +
    (-0.80 <= eta && eta < -0.70) * (6.60 <= pt * cosh(eta) && pt * cosh(eta) < 7.10) * (0.044193) +
    (-0.80 <= eta && eta < -0.70) * (7.10 <= pt * cosh(eta) && pt * cosh(eta) < 7.50) * (0.067162) +
    (-0.80 <= eta && eta < -0.70) * (7.50 <= pt * cosh(eta) && pt * cosh(eta) < 7.90) * (0.088655) +
    (-0.80 <= eta && eta < -0.70) * (7.90 <= pt * cosh(eta) && pt * cosh(eta) < 8.30) * (0.111401) +
    (-0.80 <= eta && eta < -0.70) * (8.30 <= pt * cosh(eta) && pt * cosh(eta) < 8.70) * (0.134165) +
    (-0.80 <= eta && eta < -0.70) * (8.70 <= pt * cosh(eta) && pt * cosh(eta) < 9.10) * (0.156920) +
    (-0.80 <= eta && eta < -0.70) * (9.10 <= pt * cosh(eta) && pt * cosh(eta) < 9.50) * (0.177718) +
    (-0.80 <= eta && eta < -0.70) * (9.50 <= pt * cosh(eta) && pt * cosh(eta) < 10.00) * (0.199864) +
    (-0.80 <= eta && eta < -0.70) * (10.00 <= pt * cosh(eta) && pt * cosh(eta) < 10.50) * (0.223340) +
    (-0.80 <= eta && eta < -0.70) * (10.50 <= pt * cosh(eta) && pt * cosh(eta) < 11.10) * (0.246613) +
    (-0.80 <= eta && eta < -0.70) * (11.10 <= pt * cosh(eta) && pt * cosh(eta) < 11.70) * (0.269390) +
    (-0.80 <= eta && eta < -0.70) * (11.70 <= pt * cosh(eta) && pt * cosh(eta) < 12.40) * (0.291176) +
    (-0.80 <= eta && eta < -0.70) * (12.40 <= pt * cosh(eta) && pt * cosh(eta) < 13.20) * (0.313049) +
    (-0.80 <= eta && eta < -0.70) * (13.20 <= pt * cosh(eta) && pt * cosh(eta) < 14.20) * (0.335333) +
    (-0.80 <= eta && eta < -0.70) * (14.20 <= pt * cosh(eta) && pt * cosh(eta) < 15.40) * (0.357820) +
    (-0.80 <= eta && eta < -0.70) * (15.40 <= pt * cosh(eta) && pt * cosh(eta) < 16.90) * (0.379834) +
    (-0.80 <= eta && eta < -0.70) * (16.90 <= pt * cosh(eta) && pt * cosh(eta) < 18.90) * (0.401637) +
    (-0.80 <= eta && eta < -0.70) * (18.90 <= pt * cosh(eta) && pt * cosh(eta) < 21.80) * (0.423486) +
    (-0.80 <= eta && eta < -0.70) * (21.80 <= pt * cosh(eta) && pt * cosh(eta) < 26.80) * (0.445919) +
    (-0.80 <= eta && eta < -0.70) * (26.80 <= pt * cosh(eta) && pt * cosh(eta) < 41.30) * (0.471428) +
    (-0.80 <= eta && eta < -0.70) * (41.30 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.484680) +
    (-0.70 <= eta && eta < -0.60) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 5.90) * (0.000627) +
    (-0.70 <= eta && eta < -0.60) * (5.90 <= pt * cosh(eta) && pt * cosh(eta) < 6.50) * (0.019743) +
    (-0.70 <= eta && eta < -0.60) * (6.50 <= pt * cosh(eta) && pt * cosh(eta) < 7.00) * (0.041036) +
    (-0.70 <= eta && eta < -0.60) * (7.00 <= pt * cosh(eta) && pt * cosh(eta) < 7.40) * (0.063090) +
    (-0.70 <= eta && eta < -0.60) * (7.40 <= pt * cosh(eta) && pt * cosh(eta) < 7.80) * (0.084632) +
    (-0.70 <= eta && eta < -0.60) * (7.80 <= pt * cosh(eta) && pt * cosh(eta) < 8.20) * (0.107944) +
    (-0.70 <= eta && eta < -0.60) * (8.20 <= pt * cosh(eta) && pt * cosh(eta) < 8.60) * (0.131352) +
    (-0.70 <= eta && eta < -0.60) * (8.60 <= pt * cosh(eta) && pt * cosh(eta) < 9.00) * (0.153268) +
    (-0.70 <= eta && eta < -0.60) * (9.00 <= pt * cosh(eta) && pt * cosh(eta) < 9.40) * (0.174594) +
    (-0.70 <= eta && eta < -0.60) * (9.40 <= pt * cosh(eta) && pt * cosh(eta) < 9.90) * (0.197491) +
    (-0.70 <= eta && eta < -0.60) * (9.90 <= pt * cosh(eta) && pt * cosh(eta) < 10.40) * (0.221492) +
    (-0.70 <= eta && eta < -0.60) * (10.40 <= pt * cosh(eta) && pt * cosh(eta) < 10.90) * (0.243056) +
    (-0.70 <= eta && eta < -0.60) * (10.90 <= pt * cosh(eta) && pt * cosh(eta) < 11.50) * (0.264487) +
    (-0.70 <= eta && eta < -0.60) * (11.50 <= pt * cosh(eta) && pt * cosh(eta) < 12.20) * (0.286988) +
    (-0.70 <= eta && eta < -0.60) * (12.20 <= pt * cosh(eta) && pt * cosh(eta) < 13.00) * (0.309562) +
    (-0.70 <= eta && eta < -0.60) * (13.00 <= pt * cosh(eta) && pt * cosh(eta) < 13.90) * (0.331387) +
    (-0.70 <= eta && eta < -0.60) * (13.90 <= pt * cosh(eta) && pt * cosh(eta) < 15.00) * (0.352810) +
    (-0.70 <= eta && eta < -0.60) * (15.00 <= pt * cosh(eta) && pt * cosh(eta) < 16.40) * (0.374492) +
    (-0.70 <= eta && eta < -0.60) * (16.40 <= pt * cosh(eta) && pt * cosh(eta) < 18.20) * (0.396058) +
    (-0.70 <= eta && eta < -0.60) * (18.20 <= pt * cosh(eta) && pt * cosh(eta) < 20.80) * (0.417745) +
    (-0.70 <= eta && eta < -0.60) * (20.80 <= pt * cosh(eta) && pt * cosh(eta) < 25.00) * (0.439965) +
    (-0.70 <= eta && eta < -0.60) * (25.00 <= pt * cosh(eta) && pt * cosh(eta) < 34.60) * (0.463862) +
    (-0.70 <= eta && eta < -0.60) * (34.60 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.481927) +
    (-0.60 <= eta && eta < -0.50) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 5.90) * (0.000606) +
    (-0.60 <= eta && eta < -0.50) * (5.90 <= pt * cosh(eta) && pt * cosh(eta) < 6.50) * (0.019510) +
    (-0.60 <= eta && eta < -0.50) * (6.50 <= pt * cosh(eta) && pt * cosh(eta) < 7.00) * (0.040196) +
    (-0.60 <= eta && eta < -0.50) * (7.00 <= pt * cosh(eta) && pt * cosh(eta) < 7.40) * (0.061442) +
    (-0.60 <= eta && eta < -0.50) * (7.40 <= pt * cosh(eta) && pt * cosh(eta) < 7.80) * (0.083297) +
    (-0.60 <= eta && eta < -0.50) * (7.80 <= pt * cosh(eta) && pt * cosh(eta) < 8.20) * (0.105846) +
    (-0.60 <= eta && eta < -0.50) * (8.20 <= pt * cosh(eta) && pt * cosh(eta) < 8.60) * (0.129010) +
    (-0.60 <= eta && eta < -0.50) * (8.60 <= pt * cosh(eta) && pt * cosh(eta) < 9.00) * (0.151227) +
    (-0.60 <= eta && eta < -0.50) * (9.00 <= pt * cosh(eta) && pt * cosh(eta) < 9.40) * (0.173018) +
    (-0.60 <= eta && eta < -0.50) * (9.40 <= pt * cosh(eta) && pt * cosh(eta) < 9.90) * (0.197095) +
    (-0.60 <= eta && eta < -0.50) * (9.90 <= pt * cosh(eta) && pt * cosh(eta) < 10.40) * (0.220621) +
    (-0.60 <= eta && eta < -0.50) * (10.40 <= pt * cosh(eta) && pt * cosh(eta) < 10.90) * (0.242223) +
    (-0.60 <= eta && eta < -0.50) * (10.90 <= pt * cosh(eta) && pt * cosh(eta) < 11.50) * (0.263700) +
    (-0.60 <= eta && eta < -0.50) * (11.50 <= pt * cosh(eta) && pt * cosh(eta) < 12.20) * (0.286256) +
    (-0.60 <= eta && eta < -0.50) * (12.20 <= pt * cosh(eta) && pt * cosh(eta) < 13.00) * (0.308892) +
    (-0.60 <= eta && eta < -0.50) * (13.00 <= pt * cosh(eta) && pt * cosh(eta) < 13.90) * (0.330782) +
    (-0.60 <= eta && eta < -0.50) * (13.90 <= pt * cosh(eta) && pt * cosh(eta) < 15.00) * (0.352274) +
    (-0.60 <= eta && eta < -0.50) * (15.00 <= pt * cosh(eta) && pt * cosh(eta) < 16.40) * (0.374028) +
    (-0.60 <= eta && eta < -0.50) * (16.40 <= pt * cosh(eta) && pt * cosh(eta) < 18.20) * (0.395670) +
    (-0.60 <= eta && eta < -0.50) * (18.20 <= pt * cosh(eta) && pt * cosh(eta) < 20.80) * (0.417435) +
    (-0.60 <= eta && eta < -0.50) * (20.80 <= pt * cosh(eta) && pt * cosh(eta) < 25.00) * (0.439737) +
    (-0.60 <= eta && eta < -0.50) * (25.00 <= pt * cosh(eta) && pt * cosh(eta) < 34.60) * (0.463725) +
    (-0.60 <= eta && eta < -0.50) * (34.60 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.481858) +
    (-0.50 <= eta && eta < -0.40) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 5.90) * (0.000754) +
    (-0.50 <= eta && eta < -0.40) * (5.90 <= pt * cosh(eta) && pt * cosh(eta) < 6.50) * (0.022057) +
    (-0.50 <= eta && eta < -0.40) * (6.50 <= pt * cosh(eta) && pt * cosh(eta) < 7.00) * (0.044626) +
    (-0.50 <= eta && eta < -0.40) * (7.00 <= pt * cosh(eta) && pt * cosh(eta) < 7.40) * (0.067546) +
    (-0.50 <= eta && eta < -0.40) * (7.40 <= pt * cosh(eta) && pt * cosh(eta) < 7.80) * (0.089482) +
    (-0.50 <= eta && eta < -0.40) * (7.80 <= pt * cosh(eta) && pt * cosh(eta) < 8.20) * (0.112586) +
    (-0.50 <= eta && eta < -0.40) * (8.20 <= pt * cosh(eta) && pt * cosh(eta) < 8.60) * (0.135083) +
    (-0.50 <= eta && eta < -0.40) * (8.60 <= pt * cosh(eta) && pt * cosh(eta) < 9.00) * (0.158756) +
    (-0.50 <= eta && eta < -0.40) * (9.00 <= pt * cosh(eta) && pt * cosh(eta) < 9.40) * (0.180504) +
    (-0.50 <= eta && eta < -0.40) * (9.40 <= pt * cosh(eta) && pt * cosh(eta) < 9.90) * (0.204102) +
    (-0.50 <= eta && eta < -0.40) * (9.90 <= pt * cosh(eta) && pt * cosh(eta) < 10.40) * (0.227320) +
    (-0.50 <= eta && eta < -0.40) * (10.40 <= pt * cosh(eta) && pt * cosh(eta) < 11.00) * (0.250611) +
    (-0.50 <= eta && eta < -0.40) * (11.00 <= pt * cosh(eta) && pt * cosh(eta) < 11.60) * (0.273341) +
    (-0.50 <= eta && eta < -0.40) * (11.60 <= pt * cosh(eta) && pt * cosh(eta) < 12.30) * (0.295022) +
    (-0.50 <= eta && eta < -0.40) * (12.30 <= pt * cosh(eta) && pt * cosh(eta) < 13.10) * (0.316735) +
    (-0.50 <= eta && eta < -0.40) * (13.10 <= pt * cosh(eta) && pt * cosh(eta) < 14.10) * (0.338798) +
    (-0.50 <= eta && eta < -0.40) * (14.10 <= pt * cosh(eta) && pt * cosh(eta) < 15.30) * (0.361005) +
    (-0.50 <= eta && eta < -0.40) * (15.30 <= pt * cosh(eta) && pt * cosh(eta) < 16.80) * (0.382690) +
    (-0.50 <= eta && eta < -0.40) * (16.80 <= pt * cosh(eta) && pt * cosh(eta) < 18.80) * (0.404114) +
    (-0.50 <= eta && eta < -0.40) * (18.80 <= pt * cosh(eta) && pt * cosh(eta) < 21.80) * (0.425866) +
    (-0.50 <= eta && eta < -0.40) * (21.80 <= pt * cosh(eta) && pt * cosh(eta) < 27.10) * (0.448463) +
    (-0.50 <= eta && eta < -0.40) * (27.10 <= pt * cosh(eta) && pt * cosh(eta) < 45.70) * (0.475401) +
    (-0.50 <= eta && eta < -0.40) * (45.70 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.486661) +
    (-0.40 <= eta && eta < -0.30) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 5.90) * (0.000728) +
    (-0.40 <= eta && eta < -0.30) * (5.90 <= pt * cosh(eta) && pt * cosh(eta) < 6.50) * (0.021373) +
    (-0.40 <= eta && eta < -0.30) * (6.50 <= pt * cosh(eta) && pt * cosh(eta) < 7.00) * (0.042787) +
    (-0.40 <= eta && eta < -0.30) * (7.00 <= pt * cosh(eta) && pt * cosh(eta) < 7.40) * (0.065124) +
    (-0.40 <= eta && eta < -0.30) * (7.40 <= pt * cosh(eta) && pt * cosh(eta) < 7.80) * (0.087923) +
    (-0.40 <= eta && eta < -0.30) * (7.80 <= pt * cosh(eta) && pt * cosh(eta) < 8.20) * (0.111671) +
    (-0.40 <= eta && eta < -0.30) * (8.20 <= pt * cosh(eta) && pt * cosh(eta) < 8.60) * (0.134974) +
    (-0.40 <= eta && eta < -0.30) * (8.60 <= pt * cosh(eta) && pt * cosh(eta) < 9.00) * (0.157718) +
    (-0.40 <= eta && eta < -0.30) * (9.00 <= pt * cosh(eta) && pt * cosh(eta) < 9.40) * (0.179824) +
    (-0.40 <= eta && eta < -0.30) * (9.40 <= pt * cosh(eta) && pt * cosh(eta) < 9.90) * (0.202629) +
    (-0.40 <= eta && eta < -0.30) * (9.90 <= pt * cosh(eta) && pt * cosh(eta) < 10.40) * (0.225577) +
    (-0.40 <= eta && eta < -0.30) * (10.40 <= pt * cosh(eta) && pt * cosh(eta) < 11.00) * (0.248957) +
    (-0.40 <= eta && eta < -0.30) * (11.00 <= pt * cosh(eta) && pt * cosh(eta) < 11.60) * (0.271791) +
    (-0.40 <= eta && eta < -0.30) * (11.60 <= pt * cosh(eta) && pt * cosh(eta) < 12.30) * (0.293585) +
    (-0.40 <= eta && eta < -0.30) * (12.30 <= pt * cosh(eta) && pt * cosh(eta) < 13.10) * (0.315422) +
    (-0.40 <= eta && eta < -0.30) * (13.10 <= pt * cosh(eta) && pt * cosh(eta) < 14.10) * (0.337621) +
    (-0.40 <= eta && eta < -0.30) * (14.10 <= pt * cosh(eta) && pt * cosh(eta) < 15.30) * (0.359975) +
    (-0.40 <= eta && eta < -0.30) * (15.30 <= pt * cosh(eta) && pt * cosh(eta) < 16.80) * (0.381809) +
    (-0.40 <= eta && eta < -0.30) * (16.80 <= pt * cosh(eta) && pt * cosh(eta) < 18.80) * (0.403386) +
    (-0.40 <= eta && eta < -0.30) * (18.80 <= pt * cosh(eta) && pt * cosh(eta) < 21.80) * (0.425299) +
    (-0.40 <= eta && eta < -0.30) * (21.80 <= pt * cosh(eta) && pt * cosh(eta) < 27.10) * (0.448066) +
    (-0.40 <= eta && eta < -0.30) * (27.10 <= pt * cosh(eta) && pt * cosh(eta) < 45.00) * (0.474824) +
    (-0.40 <= eta && eta < -0.30) * (45.00 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.486349) +
    (-0.30 <= eta && eta < -0.20) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 5.90) * (0.000713) +
    (-0.30 <= eta && eta < -0.20) * (5.90 <= pt * cosh(eta) && pt * cosh(eta) < 6.50) * (0.021188) +
    (-0.30 <= eta && eta < -0.20) * (6.50 <= pt * cosh(eta) && pt * cosh(eta) < 7.00) * (0.043291) +
    (-0.30 <= eta && eta < -0.20) * (7.00 <= pt * cosh(eta) && pt * cosh(eta) < 7.40) * (0.065636) +
    (-0.30 <= eta && eta < -0.20) * (7.40 <= pt * cosh(eta) && pt * cosh(eta) < 7.80) * (0.087346) +
    (-0.30 <= eta && eta < -0.20) * (7.80 <= pt * cosh(eta) && pt * cosh(eta) < 8.20) * (0.111189) +
    (-0.30 <= eta && eta < -0.20) * (8.20 <= pt * cosh(eta) && pt * cosh(eta) < 8.60) * (0.135124) +
    (-0.30 <= eta && eta < -0.20) * (8.60 <= pt * cosh(eta) && pt * cosh(eta) < 9.00) * (0.157289) +
    (-0.30 <= eta && eta < -0.20) * (9.00 <= pt * cosh(eta) && pt * cosh(eta) < 9.40) * (0.178343) +
    (-0.30 <= eta && eta < -0.20) * (9.40 <= pt * cosh(eta) && pt * cosh(eta) < 9.90) * (0.201490) +
    (-0.30 <= eta && eta < -0.20) * (9.90 <= pt * cosh(eta) && pt * cosh(eta) < 10.40) * (0.225445) +
    (-0.30 <= eta && eta < -0.20) * (10.40 <= pt * cosh(eta) && pt * cosh(eta) < 11.00) * (0.248831) +
    (-0.30 <= eta && eta < -0.20) * (11.00 <= pt * cosh(eta) && pt * cosh(eta) < 11.60) * (0.271673) +
    (-0.30 <= eta && eta < -0.20) * (11.60 <= pt * cosh(eta) && pt * cosh(eta) < 12.30) * (0.293475) +
    (-0.30 <= eta && eta < -0.20) * (12.30 <= pt * cosh(eta) && pt * cosh(eta) < 13.10) * (0.315321) +
    (-0.30 <= eta && eta < -0.20) * (13.10 <= pt * cosh(eta) && pt * cosh(eta) < 14.10) * (0.337531) +
    (-0.30 <= eta && eta < -0.20) * (14.10 <= pt * cosh(eta) && pt * cosh(eta) < 15.30) * (0.359896) +
    (-0.30 <= eta && eta < -0.20) * (15.30 <= pt * cosh(eta) && pt * cosh(eta) < 16.80) * (0.381742) +
    (-0.30 <= eta && eta < -0.20) * (16.80 <= pt * cosh(eta) && pt * cosh(eta) < 18.80) * (0.403330) +
    (-0.30 <= eta && eta < -0.20) * (18.80 <= pt * cosh(eta) && pt * cosh(eta) < 21.80) * (0.425256) +
    (-0.30 <= eta && eta < -0.20) * (21.80 <= pt * cosh(eta) && pt * cosh(eta) < 27.10) * (0.448036) +
    (-0.30 <= eta && eta < -0.20) * (27.10 <= pt * cosh(eta) && pt * cosh(eta) < 45.00) * (0.474810) +
    (-0.30 <= eta && eta < -0.20) * (45.00 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.486341) +
    (-0.20 <= eta && eta < -0.10) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 5.90) * (0.000670) +
    (-0.20 <= eta && eta < -0.10) * (5.90 <= pt * cosh(eta) && pt * cosh(eta) < 6.50) * (0.020397) +
    (-0.20 <= eta && eta < -0.10) * (6.50 <= pt * cosh(eta) && pt * cosh(eta) < 7.00) * (0.041717) +
    (-0.20 <= eta && eta < -0.10) * (7.00 <= pt * cosh(eta) && pt * cosh(eta) < 7.40) * (0.064356) +
    (-0.20 <= eta && eta < -0.10) * (7.40 <= pt * cosh(eta) && pt * cosh(eta) < 7.80) * (0.085635) +
    (-0.20 <= eta && eta < -0.10) * (7.80 <= pt * cosh(eta) && pt * cosh(eta) < 8.20) * (0.109110) +
    (-0.20 <= eta && eta < -0.10) * (8.20 <= pt * cosh(eta) && pt * cosh(eta) < 8.60) * (0.131602) +
    (-0.20 <= eta && eta < -0.10) * (8.60 <= pt * cosh(eta) && pt * cosh(eta) < 9.00) * (0.153486) +
    (-0.20 <= eta && eta < -0.10) * (9.00 <= pt * cosh(eta) && pt * cosh(eta) < 9.40) * (0.175586) +
    (-0.20 <= eta && eta < -0.10) * (9.40 <= pt * cosh(eta) && pt * cosh(eta) < 9.90) * (0.198903) +
    (-0.20 <= eta && eta < -0.10) * (9.90 <= pt * cosh(eta) && pt * cosh(eta) < 10.40) * (0.222061) +
    (-0.20 <= eta && eta < -0.10) * (10.40 <= pt * cosh(eta) && pt * cosh(eta) < 10.90) * (0.243601) +
    (-0.20 <= eta && eta < -0.10) * (10.90 <= pt * cosh(eta) && pt * cosh(eta) < 11.50) * (0.265002) +
    (-0.20 <= eta && eta < -0.10) * (11.50 <= pt * cosh(eta) && pt * cosh(eta) < 12.20) * (0.287466) +
    (-0.20 <= eta && eta < -0.10) * (12.20 <= pt * cosh(eta) && pt * cosh(eta) < 13.00) * (0.310000) +
    (-0.20 <= eta && eta < -0.10) * (13.00 <= pt * cosh(eta) && pt * cosh(eta) < 13.90) * (0.331781) +
    (-0.20 <= eta && eta < -0.10) * (13.90 <= pt * cosh(eta) && pt * cosh(eta) < 15.00) * (0.353161) +
    (-0.20 <= eta && eta < -0.10) * (15.00 <= pt * cosh(eta) && pt * cosh(eta) < 16.40) * (0.374795) +
    (-0.20 <= eta && eta < -0.10) * (16.40 <= pt * cosh(eta) && pt * cosh(eta) < 18.20) * (0.396312) +
    (-0.20 <= eta && eta < -0.10) * (18.20 <= pt * cosh(eta) && pt * cosh(eta) < 20.80) * (0.417947) +
    (-0.20 <= eta && eta < -0.10) * (20.80 <= pt * cosh(eta) && pt * cosh(eta) < 25.10) * (0.440351) +
    (-0.20 <= eta && eta < -0.10) * (25.10 <= pt * cosh(eta) && pt * cosh(eta) < 35.00) * (0.464506) +
    (-0.20 <= eta && eta < -0.10) * (35.00 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.482178) +
    (-0.10 <= eta && eta < -0.00) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 6.10) * (0.000742) +
    (-0.10 <= eta && eta < -0.00) * (6.10 <= pt * cosh(eta) && pt * cosh(eta) < 6.70) * (0.021676) +
    (-0.10 <= eta && eta < -0.00) * (6.70 <= pt * cosh(eta) && pt * cosh(eta) < 7.20) * (0.043000) +
    (-0.10 <= eta && eta < -0.00) * (7.20 <= pt * cosh(eta) && pt * cosh(eta) < 7.60) * (0.064909) +
    (-0.10 <= eta && eta < -0.00) * (7.60 <= pt * cosh(eta) && pt * cosh(eta) < 8.00) * (0.086553) +
    (-0.10 <= eta && eta < -0.00) * (8.00 <= pt * cosh(eta) && pt * cosh(eta) < 8.40) * (0.108414) +
    (-0.10 <= eta && eta < -0.00) * (8.40 <= pt * cosh(eta) && pt * cosh(eta) < 8.90) * (0.131561) +
    (-0.10 <= eta && eta < -0.00) * (8.90 <= pt * cosh(eta) && pt * cosh(eta) < 9.30) * (0.155629) +
    (-0.10 <= eta && eta < -0.00) * (9.30 <= pt * cosh(eta) && pt * cosh(eta) < 9.70) * (0.177181) +
    (-0.10 <= eta && eta < -0.00) * (9.70 <= pt * cosh(eta) && pt * cosh(eta) < 10.20) * (0.200013) +
    (-0.10 <= eta && eta < -0.00) * (10.20 <= pt * cosh(eta) && pt * cosh(eta) < 10.70) * (0.223039) +
    (-0.10 <= eta && eta < -0.00) * (10.70 <= pt * cosh(eta) && pt * cosh(eta) < 11.30) * (0.245902) +
    (-0.10 <= eta && eta < -0.00) * (11.30 <= pt * cosh(eta) && pt * cosh(eta) < 11.90) * (0.268329) +
    (-0.10 <= eta && eta < -0.00) * (11.90 <= pt * cosh(eta) && pt * cosh(eta) < 12.60) * (0.289835) +
    (-0.10 <= eta && eta < -0.00) * (12.60 <= pt * cosh(eta) && pt * cosh(eta) < 13.40) * (0.311487) +
    (-0.10 <= eta && eta < -0.00) * (13.40 <= pt * cosh(eta) && pt * cosh(eta) < 14.40) * (0.333613) +
    (-0.10 <= eta && eta < -0.00) * (14.40 <= pt * cosh(eta) && pt * cosh(eta) < 15.60) * (0.356016) +
    (-0.10 <= eta && eta < -0.00) * (15.60 <= pt * cosh(eta) && pt * cosh(eta) < 17.10) * (0.378024) +
    (-0.10 <= eta && eta < -0.00) * (17.10 <= pt * cosh(eta) && pt * cosh(eta) < 19.10) * (0.399906) +
    (-0.10 <= eta && eta < -0.00) * (19.10 <= pt * cosh(eta) && pt * cosh(eta) < 22.00) * (0.421924) +
    (-0.10 <= eta && eta < -0.00) * (22.00 <= pt * cosh(eta) && pt * cosh(eta) < 26.90) * (0.444432) +
    (-0.10 <= eta && eta < -0.00) * (26.90 <= pt * cosh(eta) && pt * cosh(eta) < 40.00) * (0.469406) +
    (-0.10 <= eta && eta < -0.00) * (40.00 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.483535) +
    (-0.00 <= eta && eta < 0.10) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 6.00) * (0.000700) +
    (-0.00 <= eta && eta < 0.10) * (6.00 <= pt * cosh(eta) && pt * cosh(eta) < 6.60) * (0.021257) +
    (-0.00 <= eta && eta < 0.10) * (6.60 <= pt * cosh(eta) && pt * cosh(eta) < 7.10) * (0.043040) +
    (-0.00 <= eta && eta < 0.10) * (7.10 <= pt * cosh(eta) && pt * cosh(eta) < 7.50) * (0.064435) +
    (-0.00 <= eta && eta < 0.10) * (7.50 <= pt * cosh(eta) && pt * cosh(eta) < 7.90) * (0.086042) +
    (-0.00 <= eta && eta < 0.10) * (7.90 <= pt * cosh(eta) && pt * cosh(eta) < 8.30) * (0.109607) +
    (-0.00 <= eta && eta < 0.10) * (8.30 <= pt * cosh(eta) && pt * cosh(eta) < 8.70) * (0.132930) +
    (-0.00 <= eta && eta < 0.10) * (8.70 <= pt * cosh(eta) && pt * cosh(eta) < 9.10) * (0.153931) +
    (-0.00 <= eta && eta < 0.10) * (9.10 <= pt * cosh(eta) && pt * cosh(eta) < 9.50) * (0.175705) +
    (-0.00 <= eta && eta < 0.10) * (9.50 <= pt * cosh(eta) && pt * cosh(eta) < 10.00) * (0.197877) +
    (-0.00 <= eta && eta < 0.10) * (10.00 <= pt * cosh(eta) && pt * cosh(eta) < 10.50) * (0.220292) +
    (-0.00 <= eta && eta < 0.10) * (10.50 <= pt * cosh(eta) && pt * cosh(eta) < 11.00) * (0.241706) +
    (-0.00 <= eta && eta < 0.10) * (11.00 <= pt * cosh(eta) && pt * cosh(eta) < 11.60) * (0.263021) +
    (-0.00 <= eta && eta < 0.10) * (11.60 <= pt * cosh(eta) && pt * cosh(eta) < 12.30) * (0.285436) +
    (-0.00 <= eta && eta < 0.10) * (12.30 <= pt * cosh(eta) && pt * cosh(eta) < 13.10) * (0.307965) +
    (-0.00 <= eta && eta < 0.10) * (13.10 <= pt * cosh(eta) && pt * cosh(eta) < 14.00) * (0.329785) +
    (-0.00 <= eta && eta < 0.10) * (14.00 <= pt * cosh(eta) && pt * cosh(eta) < 15.10) * (0.351244) +
    (-0.00 <= eta && eta < 0.10) * (15.10 <= pt * cosh(eta) && pt * cosh(eta) < 16.50) * (0.373004) +
    (-0.00 <= eta && eta < 0.10) * (16.50 <= pt * cosh(eta) && pt * cosh(eta) < 18.30) * (0.394692) +
    (-0.00 <= eta && eta < 0.10) * (18.30 <= pt * cosh(eta) && pt * cosh(eta) < 20.90) * (0.416547) +
    (-0.00 <= eta && eta < 0.10) * (20.90 <= pt * cosh(eta) && pt * cosh(eta) < 25.10) * (0.438994) +
    (-0.00 <= eta && eta < 0.10) * (25.10 <= pt * cosh(eta) && pt * cosh(eta) < 34.40) * (0.462884) +
    (-0.00 <= eta && eta < 0.10) * (34.40 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.481363) +
    (0.10 <= eta && eta < 0.20) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 5.80) * (0.000617) +
    (0.10 <= eta && eta < 0.20) * (5.80 <= pt * cosh(eta) && pt * cosh(eta) < 6.40) * (0.019051) +
    (0.10 <= eta && eta < 0.20) * (6.40 <= pt * cosh(eta) && pt * cosh(eta) < 6.90) * (0.040330) +
    (0.10 <= eta && eta < 0.20) * (6.90 <= pt * cosh(eta) && pt * cosh(eta) < 7.30) * (0.062862) +
    (0.10 <= eta && eta < 0.20) * (7.30 <= pt * cosh(eta) && pt * cosh(eta) < 7.70) * (0.085532) +
    (0.10 <= eta && eta < 0.20) * (7.70 <= pt * cosh(eta) && pt * cosh(eta) < 8.10) * (0.107916) +
    (0.10 <= eta && eta < 0.20) * (8.10 <= pt * cosh(eta) && pt * cosh(eta) < 8.50) * (0.131693) +
    (0.10 <= eta && eta < 0.20) * (8.50 <= pt * cosh(eta) && pt * cosh(eta) < 8.90) * (0.155745) +
    (0.10 <= eta && eta < 0.20) * (8.90 <= pt * cosh(eta) && pt * cosh(eta) < 9.30) * (0.177142) +
    (0.10 <= eta && eta < 0.20) * (9.30 <= pt * cosh(eta) && pt * cosh(eta) < 9.70) * (0.196263) +
    (0.10 <= eta && eta < 0.20) * (9.70 <= pt * cosh(eta) && pt * cosh(eta) < 10.20) * (0.218643) +
    (0.10 <= eta && eta < 0.20) * (10.20 <= pt * cosh(eta) && pt * cosh(eta) < 10.70) * (0.240921) +
    (0.10 <= eta && eta < 0.20) * (10.70 <= pt * cosh(eta) && pt * cosh(eta) < 11.30) * (0.262860) +
    (0.10 <= eta && eta < 0.20) * (11.30 <= pt * cosh(eta) && pt * cosh(eta) < 12.00) * (0.285861) +
    (0.10 <= eta && eta < 0.20) * (12.00 <= pt * cosh(eta) && pt * cosh(eta) < 12.80) * (0.308892) +
    (0.10 <= eta && eta < 0.20) * (12.80 <= pt * cosh(eta) && pt * cosh(eta) < 13.70) * (0.331107) +
    (0.10 <= eta && eta < 0.20) * (13.70 <= pt * cosh(eta) && pt * cosh(eta) < 14.80) * (0.352856) +
    (0.10 <= eta && eta < 0.20) * (14.80 <= pt * cosh(eta) && pt * cosh(eta) < 16.20) * (0.374803) +
    (0.10 <= eta && eta < 0.20) * (16.20 <= pt * cosh(eta) && pt * cosh(eta) < 18.00) * (0.396560) +
    (0.10 <= eta && eta < 0.20) * (18.00 <= pt * cosh(eta) && pt * cosh(eta) < 20.60) * (0.418357) +
    (0.10 <= eta && eta < 0.20) * (20.60 <= pt * cosh(eta) && pt * cosh(eta) < 24.90) * (0.440831) +
    (0.10 <= eta && eta < 0.20) * (24.90 <= pt * cosh(eta) && pt * cosh(eta) < 34.90) * (0.465032) +
    (0.10 <= eta && eta < 0.20) * (34.90 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.482582) +
    (0.20 <= eta && eta < 0.30) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 5.80) * (0.000736) +
    (0.20 <= eta && eta < 0.30) * (5.80 <= pt * cosh(eta) && pt * cosh(eta) < 6.40) * (0.021548) +
    (0.20 <= eta && eta < 0.30) * (6.40 <= pt * cosh(eta) && pt * cosh(eta) < 6.90) * (0.043947) +
    (0.20 <= eta && eta < 0.30) * (6.90 <= pt * cosh(eta) && pt * cosh(eta) < 7.30) * (0.067742) +
    (0.20 <= eta && eta < 0.30) * (7.30 <= pt * cosh(eta) && pt * cosh(eta) < 7.70) * (0.090454) +
    (0.20 <= eta && eta < 0.30) * (7.70 <= pt * cosh(eta) && pt * cosh(eta) < 8.10) * (0.113386) +
    (0.20 <= eta && eta < 0.30) * (8.10 <= pt * cosh(eta) && pt * cosh(eta) < 8.50) * (0.136843) +
    (0.20 <= eta && eta < 0.30) * (8.50 <= pt * cosh(eta) && pt * cosh(eta) < 8.90) * (0.160311) +
    (0.20 <= eta && eta < 0.30) * (8.90 <= pt * cosh(eta) && pt * cosh(eta) < 9.40) * (0.185362) +
    (0.20 <= eta && eta < 0.30) * (9.40 <= pt * cosh(eta) && pt * cosh(eta) < 9.80) * (0.207279) +
    (0.20 <= eta && eta < 0.30) * (9.80 <= pt * cosh(eta) && pt * cosh(eta) < 10.30) * (0.230256) +
    (0.20 <= eta && eta < 0.30) * (10.30 <= pt * cosh(eta) && pt * cosh(eta) < 10.90) * (0.253794) +
    (0.20 <= eta && eta < 0.30) * (10.90 <= pt * cosh(eta) && pt * cosh(eta) < 11.50) * (0.276522) +
    (0.20 <= eta && eta < 0.30) * (11.50 <= pt * cosh(eta) && pt * cosh(eta) < 12.20) * (0.298150) +
    (0.20 <= eta && eta < 0.30) * (12.20 <= pt * cosh(eta) && pt * cosh(eta) < 13.00) * (0.319758) +
    (0.20 <= eta && eta < 0.30) * (13.00 <= pt * cosh(eta) && pt * cosh(eta) < 14.00) * (0.341663) +
    (0.20 <= eta && eta < 0.30) * (14.00 <= pt * cosh(eta) && pt * cosh(eta) < 15.20) * (0.363659) +
    (0.20 <= eta && eta < 0.30) * (15.20 <= pt * cosh(eta) && pt * cosh(eta) < 16.70) * (0.385086) +
    (0.20 <= eta && eta < 0.30) * (16.70 <= pt * cosh(eta) && pt * cosh(eta) < 18.80) * (0.406698) +
    (0.20 <= eta && eta < 0.30) * (18.80 <= pt * cosh(eta) && pt * cosh(eta) < 21.90) * (0.428629) +
    (0.20 <= eta && eta < 0.30) * (21.90 <= pt * cosh(eta) && pt * cosh(eta) < 27.60) * (0.451291) +
    (0.20 <= eta && eta < 0.30) * (27.60 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.478658) +
    (0.30 <= eta && eta < 0.40) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 5.80) * (0.000753) +
    (0.30 <= eta && eta < 0.40) * (5.80 <= pt * cosh(eta) && pt * cosh(eta) < 6.40) * (0.021899) +
    (0.30 <= eta && eta < 0.40) * (6.40 <= pt * cosh(eta) && pt * cosh(eta) < 6.90) * (0.045558) +
    (0.30 <= eta && eta < 0.40) * (6.90 <= pt * cosh(eta) && pt * cosh(eta) < 7.30) * (0.068436) +
    (0.30 <= eta && eta < 0.40) * (7.30 <= pt * cosh(eta) && pt * cosh(eta) < 7.70) * (0.091463) +
    (0.30 <= eta && eta < 0.40) * (7.70 <= pt * cosh(eta) && pt * cosh(eta) < 8.10) * (0.115435) +
    (0.30 <= eta && eta < 0.40) * (8.10 <= pt * cosh(eta) && pt * cosh(eta) < 8.50) * (0.139105) +
    (0.30 <= eta && eta < 0.40) * (8.50 <= pt * cosh(eta) && pt * cosh(eta) < 8.90) * (0.162598) +
    (0.30 <= eta && eta < 0.40) * (8.90 <= pt * cosh(eta) && pt * cosh(eta) < 9.30) * (0.182795) +
    (0.30 <= eta && eta < 0.40) * (9.30 <= pt * cosh(eta) && pt * cosh(eta) < 9.80) * (0.205458) +
    (0.30 <= eta && eta < 0.40) * (9.80 <= pt * cosh(eta) && pt * cosh(eta) < 10.30) * (0.228047) +
    (0.30 <= eta && eta < 0.40) * (10.30 <= pt * cosh(eta) && pt * cosh(eta) < 10.80) * (0.249350) +
    (0.30 <= eta && eta < 0.40) * (10.80 <= pt * cosh(eta) && pt * cosh(eta) < 11.40) * (0.270616) +
    (0.30 <= eta && eta < 0.40) * (11.40 <= pt * cosh(eta) && pt * cosh(eta) < 12.10) * (0.292865) +
    (0.30 <= eta && eta < 0.40) * (12.10 <= pt * cosh(eta) && pt * cosh(eta) < 12.90) * (0.315110) +
    (0.30 <= eta && eta < 0.40) * (12.90 <= pt * cosh(eta) && pt * cosh(eta) < 13.80) * (0.336546) +
    (0.30 <= eta && eta < 0.40) * (13.80 <= pt * cosh(eta) && pt * cosh(eta) < 14.90) * (0.357525) +
    (0.30 <= eta && eta < 0.40) * (14.90 <= pt * cosh(eta) && pt * cosh(eta) < 16.30) * (0.378695) +
    (0.30 <= eta && eta < 0.40) * (16.30 <= pt * cosh(eta) && pt * cosh(eta) < 18.20) * (0.400236) +
    (0.30 <= eta && eta < 0.40) * (18.20 <= pt * cosh(eta) && pt * cosh(eta) < 21.00) * (0.422294) +
    (0.30 <= eta && eta < 0.40) * (21.00 <= pt * cosh(eta) && pt * cosh(eta) < 25.70) * (0.444843) +
    (0.30 <= eta && eta < 0.40) * (25.70 <= pt * cosh(eta) && pt * cosh(eta) < 38.60) * (0.469965) +
    (0.30 <= eta && eta < 0.40) * (38.60 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.484558) +
    (0.40 <= eta && eta < 0.50) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 5.70) * (0.000629) +
    (0.40 <= eta && eta < 0.50) * (5.70 <= pt * cosh(eta) && pt * cosh(eta) < 6.30) * (0.019727) +
    (0.40 <= eta && eta < 0.50) * (6.30 <= pt * cosh(eta) && pt * cosh(eta) < 6.80) * (0.041772) +
    (0.40 <= eta && eta < 0.50) * (6.80 <= pt * cosh(eta) && pt * cosh(eta) < 7.20) * (0.064756) +
    (0.40 <= eta && eta < 0.50) * (7.20 <= pt * cosh(eta) && pt * cosh(eta) < 7.60) * (0.087848) +
    (0.40 <= eta && eta < 0.50) * (7.60 <= pt * cosh(eta) && pt * cosh(eta) < 8.00) * (0.112141) +
    (0.40 <= eta && eta < 0.50) * (8.00 <= pt * cosh(eta) && pt * cosh(eta) < 8.40) * (0.135408) +
    (0.40 <= eta && eta < 0.50) * (8.40 <= pt * cosh(eta) && pt * cosh(eta) < 8.80) * (0.157790) +
    (0.40 <= eta && eta < 0.50) * (8.80 <= pt * cosh(eta) && pt * cosh(eta) < 9.20) * (0.180619) +
    (0.40 <= eta && eta < 0.50) * (9.20 <= pt * cosh(eta) && pt * cosh(eta) < 9.70) * (0.203357) +
    (0.40 <= eta && eta < 0.50) * (9.70 <= pt * cosh(eta) && pt * cosh(eta) < 10.20) * (0.227838) +
    (0.40 <= eta && eta < 0.50) * (10.20 <= pt * cosh(eta) && pt * cosh(eta) < 10.70) * (0.249584) +
    (0.40 <= eta && eta < 0.50) * (10.70 <= pt * cosh(eta) && pt * cosh(eta) < 11.30) * (0.271029) +
    (0.40 <= eta && eta < 0.50) * (11.30 <= pt * cosh(eta) && pt * cosh(eta) < 12.00) * (0.293437) +
    (0.40 <= eta && eta < 0.50) * (12.00 <= pt * cosh(eta) && pt * cosh(eta) < 12.80) * (0.315808) +
    (0.40 <= eta && eta < 0.50) * (12.80 <= pt * cosh(eta) && pt * cosh(eta) < 13.70) * (0.337334) +
    (0.40 <= eta && eta < 0.50) * (13.70 <= pt * cosh(eta) && pt * cosh(eta) < 14.80) * (0.358365) +
    (0.40 <= eta && eta < 0.50) * (14.80 <= pt * cosh(eta) && pt * cosh(eta) < 16.20) * (0.379551) +
    (0.40 <= eta && eta < 0.50) * (16.20 <= pt * cosh(eta) && pt * cosh(eta) < 18.10) * (0.401066) +
    (0.40 <= eta && eta < 0.50) * (18.10 <= pt * cosh(eta) && pt * cosh(eta) < 20.90) * (0.423052) +
    (0.40 <= eta && eta < 0.50) * (20.90 <= pt * cosh(eta) && pt * cosh(eta) < 25.70) * (0.445684) +
    (0.40 <= eta && eta < 0.50) * (25.70 <= pt * cosh(eta) && pt * cosh(eta) < 39.40) * (0.471163) +
    (0.40 <= eta && eta < 0.50) * (39.40 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.485174) +
    (0.50 <= eta && eta < 0.60) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 5.80) * (0.000656) +
    (0.50 <= eta && eta < 0.60) * (5.80 <= pt * cosh(eta) && pt * cosh(eta) < 6.40) * (0.020312) +
    (0.50 <= eta && eta < 0.60) * (6.40 <= pt * cosh(eta) && pt * cosh(eta) < 6.90) * (0.041999) +
    (0.50 <= eta && eta < 0.60) * (6.90 <= pt * cosh(eta) && pt * cosh(eta) < 7.30) * (0.064530) +
    (0.50 <= eta && eta < 0.60) * (7.30 <= pt * cosh(eta) && pt * cosh(eta) < 7.70) * (0.087046) +
    (0.50 <= eta && eta < 0.60) * (7.70 <= pt * cosh(eta) && pt * cosh(eta) < 8.10) * (0.110549) +
    (0.50 <= eta && eta < 0.60) * (8.10 <= pt * cosh(eta) && pt * cosh(eta) < 8.50) * (0.133710) +
    (0.50 <= eta && eta < 0.60) * (8.50 <= pt * cosh(eta) && pt * cosh(eta) < 8.90) * (0.156249) +
    (0.50 <= eta && eta < 0.60) * (8.90 <= pt * cosh(eta) && pt * cosh(eta) < 9.30) * (0.178337) +
    (0.50 <= eta && eta < 0.60) * (9.30 <= pt * cosh(eta) && pt * cosh(eta) < 9.70) * (0.199630) +
    (0.50 <= eta && eta < 0.60) * (9.70 <= pt * cosh(eta) && pt * cosh(eta) < 10.20) * (0.220540) +
    (0.50 <= eta && eta < 0.60) * (10.20 <= pt * cosh(eta) && pt * cosh(eta) < 10.70) * (0.242222) +
    (0.50 <= eta && eta < 0.60) * (10.70 <= pt * cosh(eta) && pt * cosh(eta) < 11.30) * (0.264089) +
    (0.50 <= eta && eta < 0.60) * (11.30 <= pt * cosh(eta) && pt * cosh(eta) < 12.00) * (0.287002) +
    (0.50 <= eta && eta < 0.60) * (12.00 <= pt * cosh(eta) && pt * cosh(eta) < 12.80) * (0.309935) +
    (0.50 <= eta && eta < 0.60) * (12.80 <= pt * cosh(eta) && pt * cosh(eta) < 13.70) * (0.332047) +
    (0.50 <= eta && eta < 0.60) * (13.70 <= pt * cosh(eta) && pt * cosh(eta) < 14.80) * (0.353688) +
    (0.50 <= eta && eta < 0.60) * (14.80 <= pt * cosh(eta) && pt * cosh(eta) < 16.20) * (0.375520) +
    (0.50 <= eta && eta < 0.60) * (16.20 <= pt * cosh(eta) && pt * cosh(eta) < 18.00) * (0.397160) +
    (0.50 <= eta && eta < 0.60) * (18.00 <= pt * cosh(eta) && pt * cosh(eta) < 20.60) * (0.418834) +
    (0.50 <= eta && eta < 0.60) * (20.60 <= pt * cosh(eta) && pt * cosh(eta) < 24.90) * (0.441179) +
    (0.50 <= eta && eta < 0.60) * (24.90 <= pt * cosh(eta) && pt * cosh(eta) < 35.10) * (0.465437) +
    (0.50 <= eta && eta < 0.60) * (35.10 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.482784) +
    (0.60 <= eta && eta < 0.70) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 5.80) * (0.000610) +
    (0.60 <= eta && eta < 0.70) * (5.80 <= pt * cosh(eta) && pt * cosh(eta) < 6.40) * (0.019574) +
    (0.60 <= eta && eta < 0.70) * (6.40 <= pt * cosh(eta) && pt * cosh(eta) < 6.90) * (0.040602) +
    (0.60 <= eta && eta < 0.70) * (6.90 <= pt * cosh(eta) && pt * cosh(eta) < 7.30) * (0.063109) +
    (0.60 <= eta && eta < 0.70) * (7.30 <= pt * cosh(eta) && pt * cosh(eta) < 7.70) * (0.085463) +
    (0.60 <= eta && eta < 0.70) * (7.70 <= pt * cosh(eta) && pt * cosh(eta) < 8.10) * (0.108771) +
    (0.60 <= eta && eta < 0.70) * (8.10 <= pt * cosh(eta) && pt * cosh(eta) < 8.50) * (0.132650) +
    (0.60 <= eta && eta < 0.70) * (8.50 <= pt * cosh(eta) && pt * cosh(eta) < 8.90) * (0.154383) +
    (0.60 <= eta && eta < 0.70) * (8.90 <= pt * cosh(eta) && pt * cosh(eta) < 9.30) * (0.175861) +
    (0.60 <= eta && eta < 0.70) * (9.30 <= pt * cosh(eta) && pt * cosh(eta) < 9.80) * (0.198600) +
    (0.60 <= eta && eta < 0.70) * (9.80 <= pt * cosh(eta) && pt * cosh(eta) < 10.30) * (0.222323) +
    (0.60 <= eta && eta < 0.70) * (10.30 <= pt * cosh(eta) && pt * cosh(eta) < 10.80) * (0.244096) +
    (0.60 <= eta && eta < 0.70) * (10.80 <= pt * cosh(eta) && pt * cosh(eta) < 11.40) * (0.265662) +
    (0.60 <= eta && eta < 0.70) * (11.40 <= pt * cosh(eta) && pt * cosh(eta) < 12.10) * (0.288269) +
    (0.60 <= eta && eta < 0.70) * (12.10 <= pt * cosh(eta) && pt * cosh(eta) < 12.90) * (0.310912) +
    (0.60 <= eta && eta < 0.70) * (12.90 <= pt * cosh(eta) && pt * cosh(eta) < 13.80) * (0.332765) +
    (0.60 <= eta && eta < 0.70) * (13.80 <= pt * cosh(eta) && pt * cosh(eta) < 14.90) * (0.354177) +
    (0.60 <= eta && eta < 0.70) * (14.90 <= pt * cosh(eta) && pt * cosh(eta) < 16.30) * (0.375807) +
    (0.60 <= eta && eta < 0.70) * (16.30 <= pt * cosh(eta) && pt * cosh(eta) < 18.10) * (0.397278) +
    (0.60 <= eta && eta < 0.70) * (18.10 <= pt * cosh(eta) && pt * cosh(eta) < 20.70) * (0.418823) +
    (0.60 <= eta && eta < 0.70) * (20.70 <= pt * cosh(eta) && pt * cosh(eta) < 25.00) * (0.441082) +
    (0.60 <= eta && eta < 0.70) * (25.00 <= pt * cosh(eta) && pt * cosh(eta) < 35.20) * (0.465310) +
    (0.60 <= eta && eta < 0.70) * (35.20 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.482651) +
    (0.70 <= eta && eta < 0.80) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 5.90) * (0.000768) +
    (0.70 <= eta && eta < 0.80) * (5.90 <= pt * cosh(eta) && pt * cosh(eta) < 6.50) * (0.022252) +
    (0.70 <= eta && eta < 0.80) * (6.50 <= pt * cosh(eta) && pt * cosh(eta) < 6.90) * (0.042122) +
    (0.70 <= eta && eta < 0.80) * (6.90 <= pt * cosh(eta) && pt * cosh(eta) < 7.30) * (0.062694) +
    (0.70 <= eta && eta < 0.80) * (7.30 <= pt * cosh(eta) && pt * cosh(eta) < 7.70) * (0.085610) +
    (0.70 <= eta && eta < 0.80) * (7.70 <= pt * cosh(eta) && pt * cosh(eta) < 8.10) * (0.107607) +
    (0.70 <= eta && eta < 0.80) * (8.10 <= pt * cosh(eta) && pt * cosh(eta) < 8.50) * (0.130626) +
    (0.70 <= eta && eta < 0.80) * (8.50 <= pt * cosh(eta) && pt * cosh(eta) < 8.90) * (0.154667) +
    (0.70 <= eta && eta < 0.80) * (8.90 <= pt * cosh(eta) && pt * cosh(eta) < 9.30) * (0.176433) +
    (0.70 <= eta && eta < 0.80) * (9.30 <= pt * cosh(eta) && pt * cosh(eta) < 9.80) * (0.198617) +
    (0.70 <= eta && eta < 0.80) * (9.80 <= pt * cosh(eta) && pt * cosh(eta) < 10.30) * (0.222620) +
    (0.70 <= eta && eta < 0.80) * (10.30 <= pt * cosh(eta) && pt * cosh(eta) < 10.80) * (0.244391) +
    (0.70 <= eta && eta < 0.80) * (10.80 <= pt * cosh(eta) && pt * cosh(eta) < 11.40) * (0.265940) +
    (0.70 <= eta && eta < 0.80) * (11.40 <= pt * cosh(eta) && pt * cosh(eta) < 12.10) * (0.288527) +
    (0.70 <= eta && eta < 0.80) * (12.10 <= pt * cosh(eta) && pt * cosh(eta) < 12.90) * (0.311148) +
    (0.70 <= eta && eta < 0.80) * (12.90 <= pt * cosh(eta) && pt * cosh(eta) < 13.80) * (0.332977) +
    (0.70 <= eta && eta < 0.80) * (13.80 <= pt * cosh(eta) && pt * cosh(eta) < 14.90) * (0.354366) +
    (0.70 <= eta && eta < 0.80) * (14.90 <= pt * cosh(eta) && pt * cosh(eta) < 16.30) * (0.375969) +
    (0.70 <= eta && eta < 0.80) * (16.30 <= pt * cosh(eta) && pt * cosh(eta) < 18.20) * (0.397969) +
    (0.70 <= eta && eta < 0.80) * (18.20 <= pt * cosh(eta) && pt * cosh(eta) < 20.90) * (0.420136) +
    (0.70 <= eta && eta < 0.80) * (20.90 <= pt * cosh(eta) && pt * cosh(eta) < 25.40) * (0.442635) +
    (0.70 <= eta && eta < 0.80) * (25.40 <= pt * cosh(eta) && pt * cosh(eta) < 36.60) * (0.467206) +
    (0.70 <= eta && eta < 0.80) * (36.60 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.483337) +
    (0.80 <= eta && eta < 0.90) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 6.00) * (0.000728) +
    (0.80 <= eta && eta < 0.90) * (6.00 <= pt * cosh(eta) && pt * cosh(eta) < 6.60) * (0.021330) +
    (0.80 <= eta && eta < 0.90) * (6.60 <= pt * cosh(eta) && pt * cosh(eta) < 7.10) * (0.042978) +
    (0.80 <= eta && eta < 0.90) * (7.10 <= pt * cosh(eta) && pt * cosh(eta) < 7.50) * (0.065797) +
    (0.80 <= eta && eta < 0.90) * (7.50 <= pt * cosh(eta) && pt * cosh(eta) < 7.90) * (0.087589) +
    (0.80 <= eta && eta < 0.90) * (7.90 <= pt * cosh(eta) && pt * cosh(eta) < 8.30) * (0.109904) +
    (0.80 <= eta && eta < 0.90) * (8.30 <= pt * cosh(eta) && pt * cosh(eta) < 8.70) * (0.132432) +
    (0.80 <= eta && eta < 0.90) * (8.70 <= pt * cosh(eta) && pt * cosh(eta) < 9.10) * (0.155030) +
    (0.80 <= eta && eta < 0.90) * (9.10 <= pt * cosh(eta) && pt * cosh(eta) < 9.50) * (0.177015) +
    (0.80 <= eta && eta < 0.90) * (9.50 <= pt * cosh(eta) && pt * cosh(eta) < 10.00) * (0.198109) +
    (0.80 <= eta && eta < 0.90) * (10.00 <= pt * cosh(eta) && pt * cosh(eta) < 10.50) * (0.220572) +
    (0.80 <= eta && eta < 0.90) * (10.50 <= pt * cosh(eta) && pt * cosh(eta) < 11.00) * (0.241975) +
    (0.80 <= eta && eta < 0.90) * (11.00 <= pt * cosh(eta) && pt * cosh(eta) < 11.60) * (0.263275) +
    (0.80 <= eta && eta < 0.90) * (11.60 <= pt * cosh(eta) && pt * cosh(eta) < 12.30) * (0.285672) +
    (0.80 <= eta && eta < 0.90) * (12.30 <= pt * cosh(eta) && pt * cosh(eta) < 13.10) * (0.308182) +
    (0.80 <= eta && eta < 0.90) * (13.10 <= pt * cosh(eta) && pt * cosh(eta) < 14.00) * (0.329981) +
    (0.80 <= eta && eta < 0.90) * (14.00 <= pt * cosh(eta) && pt * cosh(eta) < 15.10) * (0.351418) +
    (0.80 <= eta && eta < 0.90) * (15.10 <= pt * cosh(eta) && pt * cosh(eta) < 16.50) * (0.373154) +
    (0.80 <= eta && eta < 0.90) * (16.50 <= pt * cosh(eta) && pt * cosh(eta) < 18.30) * (0.394818) +
    (0.80 <= eta && eta < 0.90) * (18.30 <= pt * cosh(eta) && pt * cosh(eta) < 20.90) * (0.416648) +
    (0.80 <= eta && eta < 0.90) * (20.90 <= pt * cosh(eta) && pt * cosh(eta) < 25.10) * (0.439069) +
    (0.80 <= eta && eta < 0.90) * (25.10 <= pt * cosh(eta) && pt * cosh(eta) < 34.50) * (0.463037) +
    (0.80 <= eta && eta < 0.90) * (34.50 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.481439) +
    (0.90 <= eta && eta < 1.00) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 6.00) * (0.000687) +
    (0.90 <= eta && eta < 1.00) * (6.00 <= pt * cosh(eta) && pt * cosh(eta) < 6.60) * (0.020501) +
    (0.90 <= eta && eta < 1.00) * (6.60 <= pt * cosh(eta) && pt * cosh(eta) < 7.10) * (0.041727) +
    (0.90 <= eta && eta < 1.00) * (7.10 <= pt * cosh(eta) && pt * cosh(eta) < 7.50) * (0.063800) +
    (0.90 <= eta && eta < 1.00) * (7.50 <= pt * cosh(eta) && pt * cosh(eta) < 7.90) * (0.084994) +
    (0.90 <= eta && eta < 1.00) * (7.90 <= pt * cosh(eta) && pt * cosh(eta) < 8.30) * (0.107765) +
    (0.90 <= eta && eta < 1.00) * (8.30 <= pt * cosh(eta) && pt * cosh(eta) < 8.70) * (0.130564) +
    (0.90 <= eta && eta < 1.00) * (8.70 <= pt * cosh(eta) && pt * cosh(eta) < 9.10) * (0.152720) +
    (0.90 <= eta && eta < 1.00) * (9.10 <= pt * cosh(eta) && pt * cosh(eta) < 9.50) * (0.174241) +
    (0.90 <= eta && eta < 1.00) * (9.50 <= pt * cosh(eta) && pt * cosh(eta) < 10.00) * (0.197247) +
    (0.90 <= eta && eta < 1.00) * (10.00 <= pt * cosh(eta) && pt * cosh(eta) < 10.50) * (0.220892) +
    (0.90 <= eta && eta < 1.00) * (10.50 <= pt * cosh(eta) && pt * cosh(eta) < 11.00) * (0.242281) +
    (0.90 <= eta && eta < 1.00) * (11.00 <= pt * cosh(eta) && pt * cosh(eta) < 11.60) * (0.263564) +
    (0.90 <= eta && eta < 1.00) * (11.60 <= pt * cosh(eta) && pt * cosh(eta) < 12.30) * (0.285942) +
    (0.90 <= eta && eta < 1.00) * (12.30 <= pt * cosh(eta) && pt * cosh(eta) < 13.10) * (0.308428) +
    (0.90 <= eta && eta < 1.00) * (13.10 <= pt * cosh(eta) && pt * cosh(eta) < 14.00) * (0.330204) +
    (0.90 <= eta && eta < 1.00) * (14.00 <= pt * cosh(eta) && pt * cosh(eta) < 15.10) * (0.351616) +
    (0.90 <= eta && eta < 1.00) * (15.10 <= pt * cosh(eta) && pt * cosh(eta) < 16.50) * (0.373325) +
    (0.90 <= eta && eta < 1.00) * (16.50 <= pt * cosh(eta) && pt * cosh(eta) < 18.30) * (0.394961) +
    (0.90 <= eta && eta < 1.00) * (18.30 <= pt * cosh(eta) && pt * cosh(eta) < 20.90) * (0.416763) +
    (0.90 <= eta && eta < 1.00) * (20.90 <= pt * cosh(eta) && pt * cosh(eta) < 25.10) * (0.439153) +
    (0.90 <= eta && eta < 1.00) * (25.10 <= pt * cosh(eta) && pt * cosh(eta) < 34.50) * (0.463089) +
    (0.90 <= eta && eta < 1.00) * (34.50 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.481465) +
    (1.00 <= eta && eta < 1.10) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.000000) }

  add EfficiencyFormula {-211} {211} {
    (eta<-1.00 || eta>=1.00 || pt * cosh(eta) < 0.10 || pt * cosh(eta) >= 50.00) * ( 0.00 ) +
    (-1.00 <= eta && eta < -0.90) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) +
    (-1.00 <= eta && eta < -0.90) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 5.70) * (0.999363) +
    (-1.00 <= eta && eta < -0.90) * (5.70 <= pt * cosh(eta) && pt * cosh(eta) < 6.20) * (0.977469) +
    (-1.00 <= eta && eta < -0.90) * (6.20 <= pt * cosh(eta) && pt * cosh(eta) < 6.60) * (0.953824) +
    (-1.00 <= eta && eta < -0.90) * (6.60 <= pt * cosh(eta) && pt * cosh(eta) < 6.90) * (0.930396) +
    (-1.00 <= eta && eta < -0.90) * (6.90 <= pt * cosh(eta) && pt * cosh(eta) < 7.20) * (0.908335) +
    (-1.00 <= eta && eta < -0.90) * (7.20 <= pt * cosh(eta) && pt * cosh(eta) < 7.50) * (0.885172) +
    (-1.00 <= eta && eta < -0.90) * (7.50 <= pt * cosh(eta) && pt * cosh(eta) < 7.80) * (0.861520) +
    (-1.00 <= eta && eta < -0.90) * (7.80 <= pt * cosh(eta) && pt * cosh(eta) < 8.10) * (0.838567) +
    (-1.00 <= eta && eta < -0.90) * (8.10 <= pt * cosh(eta) && pt * cosh(eta) < 8.40) * (0.817634) +
    (-1.00 <= eta && eta < -0.90) * (8.40 <= pt * cosh(eta) && pt * cosh(eta) < 8.80) * (0.793978) +
    (-1.00 <= eta && eta < -0.90) * (8.80 <= pt * cosh(eta) && pt * cosh(eta) < 9.20) * (0.769187) +
    (-1.00 <= eta && eta < -0.90) * (9.20 <= pt * cosh(eta) && pt * cosh(eta) < 9.70) * (0.745592) +
    (-1.00 <= eta && eta < -0.90) * (9.70 <= pt * cosh(eta) && pt * cosh(eta) < 10.30) * (0.721481) +
    (-1.00 <= eta && eta < -0.90) * (10.30 <= pt * cosh(eta) && pt * cosh(eta) < 11.00) * (0.698201) +
    (-1.00 <= eta && eta < -0.90) * (11.00 <= pt * cosh(eta) && pt * cosh(eta) < 11.90) * (0.676260) +
    (-1.00 <= eta && eta < -0.90) * (11.90 <= pt * cosh(eta) && pt * cosh(eta) < 13.30) * (0.654166) +
    (-1.00 <= eta && eta < -0.90) * (13.30 <= pt * cosh(eta) && pt * cosh(eta) < 16.50) * (0.629819) +
    (-1.00 <= eta && eta < -0.90) * (16.50 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.604030) +
    (-0.90 <= eta && eta < -0.80) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) +
    (-0.90 <= eta && eta < -0.80) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 5.60) * (0.999415) +
    (-0.90 <= eta && eta < -0.80) * (5.60 <= pt * cosh(eta) && pt * cosh(eta) < 6.10) * (0.978320) +
    (-0.90 <= eta && eta < -0.80) * (6.10 <= pt * cosh(eta) && pt * cosh(eta) < 6.50) * (0.954946) +
    (-0.90 <= eta && eta < -0.80) * (6.50 <= pt * cosh(eta) && pt * cosh(eta) < 6.80) * (0.931253) +
    (-0.90 <= eta && eta < -0.80) * (6.80 <= pt * cosh(eta) && pt * cosh(eta) < 7.10) * (0.908742) +
    (-0.90 <= eta && eta < -0.80) * (7.10 <= pt * cosh(eta) && pt * cosh(eta) < 7.40) * (0.884654) +
    (-0.90 <= eta && eta < -0.80) * (7.40 <= pt * cosh(eta) && pt * cosh(eta) < 7.70) * (0.860913) +
    (-0.90 <= eta && eta < -0.80) * (7.70 <= pt * cosh(eta) && pt * cosh(eta) < 8.00) * (0.837747) +
    (-0.90 <= eta && eta < -0.80) * (8.00 <= pt * cosh(eta) && pt * cosh(eta) < 8.30) * (0.815568) +
    (-0.90 <= eta && eta < -0.80) * (8.30 <= pt * cosh(eta) && pt * cosh(eta) < 8.70) * (0.791569) +
    (-0.90 <= eta && eta < -0.80) * (8.70 <= pt * cosh(eta) && pt * cosh(eta) < 9.10) * (0.767848) +
    (-0.90 <= eta && eta < -0.80) * (9.10 <= pt * cosh(eta) && pt * cosh(eta) < 9.50) * (0.747406) +
    (-0.90 <= eta && eta < -0.80) * (9.50 <= pt * cosh(eta) && pt * cosh(eta) < 10.10) * (0.723516) +
    (-0.90 <= eta && eta < -0.80) * (10.10 <= pt * cosh(eta) && pt * cosh(eta) < 10.80) * (0.699697) +
    (-0.90 <= eta && eta < -0.80) * (10.80 <= pt * cosh(eta) && pt * cosh(eta) < 11.70) * (0.677096) +
    (-0.90 <= eta && eta < -0.80) * (11.70 <= pt * cosh(eta) && pt * cosh(eta) < 13.10) * (0.654458) +
    (-0.90 <= eta && eta < -0.80) * (13.10 <= pt * cosh(eta) && pt * cosh(eta) < 16.20) * (0.630047) +
    (-0.90 <= eta && eta < -0.80) * (16.20 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.603999) +
    (-0.80 <= eta && eta < -0.70) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) +
    (-0.80 <= eta && eta < -0.70) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 5.50) * (0.999438) +
    (-0.80 <= eta && eta < -0.70) * (5.50 <= pt * cosh(eta) && pt * cosh(eta) < 6.00) * (0.978782) +
    (-0.80 <= eta && eta < -0.70) * (6.00 <= pt * cosh(eta) && pt * cosh(eta) < 6.40) * (0.953835) +
    (-0.80 <= eta && eta < -0.70) * (6.40 <= pt * cosh(eta) && pt * cosh(eta) < 6.70) * (0.929657) +
    (-0.80 <= eta && eta < -0.70) * (6.70 <= pt * cosh(eta) && pt * cosh(eta) < 7.00) * (0.907134) +
    (-0.80 <= eta && eta < -0.70) * (7.00 <= pt * cosh(eta) && pt * cosh(eta) < 7.30) * (0.882151) +
    (-0.80 <= eta && eta < -0.70) * (7.30 <= pt * cosh(eta) && pt * cosh(eta) < 7.60) * (0.858132) +
    (-0.80 <= eta && eta < -0.70) * (7.60 <= pt * cosh(eta) && pt * cosh(eta) < 7.90) * (0.835697) +
    (-0.80 <= eta && eta < -0.70) * (7.90 <= pt * cosh(eta) && pt * cosh(eta) < 8.20) * (0.813827) +
    (-0.80 <= eta && eta < -0.70) * (8.20 <= pt * cosh(eta) && pt * cosh(eta) < 8.60) * (0.790188) +
    (-0.80 <= eta && eta < -0.70) * (8.60 <= pt * cosh(eta) && pt * cosh(eta) < 9.00) * (0.765614) +
    (-0.80 <= eta && eta < -0.70) * (9.00 <= pt * cosh(eta) && pt * cosh(eta) < 9.50) * (0.742033) +
    (-0.80 <= eta && eta < -0.70) * (9.50 <= pt * cosh(eta) && pt * cosh(eta) < 10.00) * (0.720207) +
    (-0.80 <= eta && eta < -0.70) * (10.00 <= pt * cosh(eta) && pt * cosh(eta) < 10.70) * (0.698161) +
    (-0.80 <= eta && eta < -0.70) * (10.70 <= pt * cosh(eta) && pt * cosh(eta) < 11.60) * (0.675683) +
    (-0.80 <= eta && eta < -0.70) * (11.60 <= pt * cosh(eta) && pt * cosh(eta) < 13.00) * (0.653272) +
    (-0.80 <= eta && eta < -0.70) * (13.00 <= pt * cosh(eta) && pt * cosh(eta) < 16.30) * (0.628628) +
    (-0.80 <= eta && eta < -0.70) * (16.30 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.603794) +
    (-0.70 <= eta && eta < -0.60) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) +
    (-0.70 <= eta && eta < -0.60) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 5.50) * (0.999386) +
    (-0.70 <= eta && eta < -0.60) * (5.50 <= pt * cosh(eta) && pt * cosh(eta) < 6.00) * (0.977686) +
    (-0.70 <= eta && eta < -0.60) * (6.00 <= pt * cosh(eta) && pt * cosh(eta) < 6.40) * (0.952568) +
    (-0.70 <= eta && eta < -0.60) * (6.40 <= pt * cosh(eta) && pt * cosh(eta) < 6.70) * (0.927697) +
    (-0.70 <= eta && eta < -0.60) * (6.70 <= pt * cosh(eta) && pt * cosh(eta) < 7.00) * (0.904394) +
    (-0.70 <= eta && eta < -0.60) * (7.00 <= pt * cosh(eta) && pt * cosh(eta) < 7.30) * (0.880318) +
    (-0.70 <= eta && eta < -0.60) * (7.30 <= pt * cosh(eta) && pt * cosh(eta) < 7.60) * (0.856544) +
    (-0.70 <= eta && eta < -0.60) * (7.60 <= pt * cosh(eta) && pt * cosh(eta) < 7.90) * (0.833257) +
    (-0.70 <= eta && eta < -0.60) * (7.90 <= pt * cosh(eta) && pt * cosh(eta) < 8.20) * (0.810933) +
    (-0.70 <= eta && eta < -0.60) * (8.20 <= pt * cosh(eta) && pt * cosh(eta) < 8.60) * (0.786986) +
    (-0.70 <= eta && eta < -0.60) * (8.60 <= pt * cosh(eta) && pt * cosh(eta) < 9.00) * (0.763520) +
    (-0.70 <= eta && eta < -0.60) * (9.00 <= pt * cosh(eta) && pt * cosh(eta) < 9.50) * (0.740261) +
    (-0.70 <= eta && eta < -0.60) * (9.50 <= pt * cosh(eta) && pt * cosh(eta) < 10.00) * (0.718074) +
    (-0.70 <= eta && eta < -0.60) * (10.00 <= pt * cosh(eta) && pt * cosh(eta) < 10.70) * (0.696228) +
    (-0.70 <= eta && eta < -0.60) * (10.70 <= pt * cosh(eta) && pt * cosh(eta) < 11.60) * (0.674142) +
    (-0.70 <= eta && eta < -0.60) * (11.60 <= pt * cosh(eta) && pt * cosh(eta) < 13.10) * (0.651446) +
    (-0.70 <= eta && eta < -0.60) * (13.10 <= pt * cosh(eta) && pt * cosh(eta) < 17.00) * (0.625653) +
    (-0.70 <= eta && eta < -0.60) * (17.00 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.603457) +
    (-0.60 <= eta && eta < -0.50) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) +
    (-0.60 <= eta && eta < -0.50) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 5.50) * (0.999412) +
    (-0.60 <= eta && eta < -0.50) * (5.50 <= pt * cosh(eta) && pt * cosh(eta) < 6.00) * (0.978142) +
    (-0.60 <= eta && eta < -0.50) * (6.00 <= pt * cosh(eta) && pt * cosh(eta) < 6.40) * (0.953041) +
    (-0.60 <= eta && eta < -0.50) * (6.40 <= pt * cosh(eta) && pt * cosh(eta) < 6.70) * (0.928653) +
    (-0.60 <= eta && eta < -0.50) * (6.70 <= pt * cosh(eta) && pt * cosh(eta) < 7.00) * (0.906049) +
    (-0.60 <= eta && eta < -0.50) * (7.00 <= pt * cosh(eta) && pt * cosh(eta) < 7.30) * (0.882734) +
    (-0.60 <= eta && eta < -0.50) * (7.30 <= pt * cosh(eta) && pt * cosh(eta) < 7.60) * (0.858632) +
    (-0.60 <= eta && eta < -0.50) * (7.60 <= pt * cosh(eta) && pt * cosh(eta) < 7.90) * (0.835440) +
    (-0.60 <= eta && eta < -0.50) * (7.90 <= pt * cosh(eta) && pt * cosh(eta) < 8.20) * (0.813533) +
    (-0.60 <= eta && eta < -0.50) * (8.20 <= pt * cosh(eta) && pt * cosh(eta) < 8.60) * (0.789624) +
    (-0.60 <= eta && eta < -0.50) * (8.60 <= pt * cosh(eta) && pt * cosh(eta) < 9.00) * (0.765618) +
    (-0.60 <= eta && eta < -0.50) * (9.00 <= pt * cosh(eta) && pt * cosh(eta) < 9.40) * (0.744122) +
    (-0.60 <= eta && eta < -0.50) * (9.40 <= pt * cosh(eta) && pt * cosh(eta) < 9.90) * (0.722552) +
    (-0.60 <= eta && eta < -0.50) * (9.90 <= pt * cosh(eta) && pt * cosh(eta) < 10.60) * (0.700157) +
    (-0.60 <= eta && eta < -0.50) * (10.60 <= pt * cosh(eta) && pt * cosh(eta) < 11.50) * (0.677085) +
    (-0.60 <= eta && eta < -0.50) * (11.50 <= pt * cosh(eta) && pt * cosh(eta) < 12.90) * (0.654123) +
    (-0.60 <= eta && eta < -0.50) * (12.90 <= pt * cosh(eta) && pt * cosh(eta) < 16.10) * (0.629276) +
    (-0.60 <= eta && eta < -0.50) * (16.10 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.603840) +
    (-0.50 <= eta && eta < -0.40) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) +
    (-0.50 <= eta && eta < -0.40) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 5.40) * (0.999451) +
    (-0.50 <= eta && eta < -0.40) * (5.40 <= pt * cosh(eta) && pt * cosh(eta) < 5.90) * (0.979067) +
    (-0.50 <= eta && eta < -0.40) * (5.90 <= pt * cosh(eta) && pt * cosh(eta) < 6.20) * (0.958187) +
    (-0.50 <= eta && eta < -0.40) * (6.20 <= pt * cosh(eta) && pt * cosh(eta) < 6.50) * (0.937021) +
    (-0.50 <= eta && eta < -0.40) * (6.50 <= pt * cosh(eta) && pt * cosh(eta) < 6.80) * (0.913728) +
    (-0.50 <= eta && eta < -0.40) * (6.80 <= pt * cosh(eta) && pt * cosh(eta) < 7.10) * (0.889767) +
    (-0.50 <= eta && eta < -0.40) * (7.10 <= pt * cosh(eta) && pt * cosh(eta) < 7.40) * (0.865592) +
    (-0.50 <= eta && eta < -0.40) * (7.40 <= pt * cosh(eta) && pt * cosh(eta) < 7.70) * (0.842135) +
    (-0.50 <= eta && eta < -0.40) * (7.70 <= pt * cosh(eta) && pt * cosh(eta) < 8.00) * (0.819614) +
    (-0.50 <= eta && eta < -0.40) * (8.00 <= pt * cosh(eta) && pt * cosh(eta) < 8.40) * (0.795484) +
    (-0.50 <= eta && eta < -0.40) * (8.40 <= pt * cosh(eta) && pt * cosh(eta) < 8.80) * (0.770150) +
    (-0.50 <= eta && eta < -0.40) * (8.80 <= pt * cosh(eta) && pt * cosh(eta) < 9.20) * (0.747397) +
    (-0.50 <= eta && eta < -0.40) * (9.20 <= pt * cosh(eta) && pt * cosh(eta) < 9.70) * (0.725362) +
    (-0.50 <= eta && eta < -0.40) * (9.70 <= pt * cosh(eta) && pt * cosh(eta) < 10.30) * (0.703515) +
    (-0.50 <= eta && eta < -0.40) * (10.30 <= pt * cosh(eta) && pt * cosh(eta) < 11.10) * (0.681967) +
    (-0.50 <= eta && eta < -0.40) * (11.10 <= pt * cosh(eta) && pt * cosh(eta) < 12.30) * (0.659657) +
    (-0.50 <= eta && eta < -0.40) * (12.30 <= pt * cosh(eta) && pt * cosh(eta) < 14.60) * (0.636152) +
    (-0.50 <= eta && eta < -0.40) * (14.60 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.604442) +
    (-0.40 <= eta && eta < -0.30) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) +
    (-0.40 <= eta && eta < -0.30) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 5.40) * (0.999472) +
    (-0.40 <= eta && eta < -0.30) * (5.40 <= pt * cosh(eta) && pt * cosh(eta) < 5.90) * (0.979617) +
    (-0.40 <= eta && eta < -0.30) * (5.90 <= pt * cosh(eta) && pt * cosh(eta) < 6.30) * (0.955630) +
    (-0.40 <= eta && eta < -0.30) * (6.30 <= pt * cosh(eta) && pt * cosh(eta) < 6.60) * (0.932152) +
    (-0.40 <= eta && eta < -0.30) * (6.60 <= pt * cosh(eta) && pt * cosh(eta) < 6.90) * (0.909173) +
    (-0.40 <= eta && eta < -0.30) * (6.90 <= pt * cosh(eta) && pt * cosh(eta) < 7.20) * (0.885363) +
    (-0.40 <= eta && eta < -0.30) * (7.20 <= pt * cosh(eta) && pt * cosh(eta) < 7.50) * (0.860582) +
    (-0.40 <= eta && eta < -0.30) * (7.50 <= pt * cosh(eta) && pt * cosh(eta) < 7.80) * (0.836574) +
    (-0.40 <= eta && eta < -0.30) * (7.80 <= pt * cosh(eta) && pt * cosh(eta) < 8.10) * (0.813641) +
    (-0.40 <= eta && eta < -0.30) * (8.10 <= pt * cosh(eta) && pt * cosh(eta) < 8.50) * (0.789231) +
    (-0.40 <= eta && eta < -0.30) * (8.50 <= pt * cosh(eta) && pt * cosh(eta) < 8.90) * (0.764564) +
    (-0.40 <= eta && eta < -0.30) * (8.90 <= pt * cosh(eta) && pt * cosh(eta) < 9.40) * (0.740178) +
    (-0.40 <= eta && eta < -0.30) * (9.40 <= pt * cosh(eta) && pt * cosh(eta) < 10.00) * (0.715912) +
    (-0.40 <= eta && eta < -0.30) * (10.00 <= pt * cosh(eta) && pt * cosh(eta) < 10.70) * (0.693289) +
    (-0.40 <= eta && eta < -0.30) * (10.70 <= pt * cosh(eta) && pt * cosh(eta) < 11.70) * (0.670721) +
    (-0.40 <= eta && eta < -0.30) * (11.70 <= pt * cosh(eta) && pt * cosh(eta) < 13.30) * (0.647669) +
    (-0.40 <= eta && eta < -0.30) * (13.30 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.605441) +
    (-0.30 <= eta && eta < -0.20) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) +
    (-0.30 <= eta && eta < -0.20) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 5.40) * (0.999491) +
    (-0.30 <= eta && eta < -0.20) * (5.40 <= pt * cosh(eta) && pt * cosh(eta) < 5.90) * (0.979938) +
    (-0.30 <= eta && eta < -0.20) * (5.90 <= pt * cosh(eta) && pt * cosh(eta) < 6.20) * (0.959867) +
    (-0.30 <= eta && eta < -0.20) * (6.20 <= pt * cosh(eta) && pt * cosh(eta) < 6.50) * (0.938816) +
    (-0.30 <= eta && eta < -0.20) * (6.50 <= pt * cosh(eta) && pt * cosh(eta) < 6.80) * (0.916434) +
    (-0.30 <= eta && eta < -0.20) * (6.80 <= pt * cosh(eta) && pt * cosh(eta) < 7.10) * (0.892350) +
    (-0.30 <= eta && eta < -0.20) * (7.10 <= pt * cosh(eta) && pt * cosh(eta) < 7.40) * (0.868247) +
    (-0.30 <= eta && eta < -0.20) * (7.40 <= pt * cosh(eta) && pt * cosh(eta) < 7.70) * (0.845108) +
    (-0.30 <= eta && eta < -0.20) * (7.70 <= pt * cosh(eta) && pt * cosh(eta) < 8.00) * (0.821653) +
    (-0.30 <= eta && eta < -0.20) * (8.00 <= pt * cosh(eta) && pt * cosh(eta) < 8.30) * (0.799669) +
    (-0.30 <= eta && eta < -0.20) * (8.30 <= pt * cosh(eta) && pt * cosh(eta) < 8.70) * (0.776530) +
    (-0.30 <= eta && eta < -0.20) * (8.70 <= pt * cosh(eta) && pt * cosh(eta) < 9.10) * (0.754283) +
    (-0.30 <= eta && eta < -0.20) * (9.10 <= pt * cosh(eta) && pt * cosh(eta) < 9.60) * (0.732042) +
    (-0.30 <= eta && eta < -0.20) * (9.60 <= pt * cosh(eta) && pt * cosh(eta) < 10.20) * (0.708971) +
    (-0.30 <= eta && eta < -0.20) * (10.20 <= pt * cosh(eta) && pt * cosh(eta) < 11.00) * (0.685968) +
    (-0.30 <= eta && eta < -0.20) * (11.00 <= pt * cosh(eta) && pt * cosh(eta) < 12.10) * (0.663406) +
    (-0.30 <= eta && eta < -0.20) * (12.10 <= pt * cosh(eta) && pt * cosh(eta) < 14.10) * (0.640318) +
    (-0.30 <= eta && eta < -0.20) * (14.10 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.604825) +
    (-0.20 <= eta && eta < -0.10) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) +
    (-0.20 <= eta && eta < -0.10) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 5.50) * (0.999335) +
    (-0.20 <= eta && eta < -0.10) * (5.50 <= pt * cosh(eta) && pt * cosh(eta) < 6.00) * (0.976774) +
    (-0.20 <= eta && eta < -0.10) * (6.00 <= pt * cosh(eta) && pt * cosh(eta) < 6.40) * (0.951314) +
    (-0.20 <= eta && eta < -0.10) * (6.40 <= pt * cosh(eta) && pt * cosh(eta) < 6.70) * (0.926371) +
    (-0.20 <= eta && eta < -0.10) * (6.70 <= pt * cosh(eta) && pt * cosh(eta) < 7.00) * (0.902907) +
    (-0.20 <= eta && eta < -0.10) * (7.00 <= pt * cosh(eta) && pt * cosh(eta) < 7.30) * (0.878307) +
    (-0.20 <= eta && eta < -0.10) * (7.30 <= pt * cosh(eta) && pt * cosh(eta) < 7.60) * (0.855267) +
    (-0.20 <= eta && eta < -0.10) * (7.60 <= pt * cosh(eta) && pt * cosh(eta) < 7.90) * (0.831969) +
    (-0.20 <= eta && eta < -0.10) * (7.90 <= pt * cosh(eta) && pt * cosh(eta) < 8.30) * (0.806045) +
    (-0.20 <= eta && eta < -0.10) * (8.30 <= pt * cosh(eta) && pt * cosh(eta) < 8.70) * (0.780629) +
    (-0.20 <= eta && eta < -0.10) * (8.70 <= pt * cosh(eta) && pt * cosh(eta) < 9.10) * (0.757736) +
    (-0.20 <= eta && eta < -0.10) * (9.10 <= pt * cosh(eta) && pt * cosh(eta) < 9.60) * (0.734043) +
    (-0.20 <= eta && eta < -0.10) * (9.60 <= pt * cosh(eta) && pt * cosh(eta) < 10.20) * (0.711369) +
    (-0.20 <= eta && eta < -0.10) * (10.20 <= pt * cosh(eta) && pt * cosh(eta) < 11.00) * (0.688232) +
    (-0.20 <= eta && eta < -0.10) * (11.00 <= pt * cosh(eta) && pt * cosh(eta) < 12.10) * (0.665130) +
    (-0.20 <= eta && eta < -0.10) * (12.10 <= pt * cosh(eta) && pt * cosh(eta) < 14.00) * (0.641972) +
    (-0.20 <= eta && eta < -0.10) * (14.00 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.605010) +
    (-0.10 <= eta && eta < -0.00) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) +
    (-0.10 <= eta && eta < -0.00) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 5.60) * (0.999436) +
    (-0.10 <= eta && eta < -0.00) * (5.60 <= pt * cosh(eta) && pt * cosh(eta) < 6.10) * (0.978964) +
    (-0.10 <= eta && eta < -0.00) * (6.10 <= pt * cosh(eta) && pt * cosh(eta) < 6.40) * (0.958369) +
    (-0.10 <= eta && eta < -0.00) * (6.40 <= pt * cosh(eta) && pt * cosh(eta) < 6.70) * (0.938389) +
    (-0.10 <= eta && eta < -0.00) * (6.70 <= pt * cosh(eta) && pt * cosh(eta) < 7.00) * (0.916194) +
    (-0.10 <= eta && eta < -0.00) * (7.00 <= pt * cosh(eta) && pt * cosh(eta) < 7.30) * (0.893059) +
    (-0.10 <= eta && eta < -0.00) * (7.30 <= pt * cosh(eta) && pt * cosh(eta) < 7.60) * (0.869718) +
    (-0.10 <= eta && eta < -0.00) * (7.60 <= pt * cosh(eta) && pt * cosh(eta) < 7.90) * (0.846143) +
    (-0.10 <= eta && eta < -0.00) * (7.90 <= pt * cosh(eta) && pt * cosh(eta) < 8.20) * (0.825154) +
    (-0.10 <= eta && eta < -0.00) * (8.20 <= pt * cosh(eta) && pt * cosh(eta) < 8.60) * (0.800818) +
    (-0.10 <= eta && eta < -0.00) * (8.60 <= pt * cosh(eta) && pt * cosh(eta) < 9.00) * (0.778480) +
    (-0.10 <= eta && eta < -0.00) * (9.00 <= pt * cosh(eta) && pt * cosh(eta) < 9.40) * (0.755699) +
    (-0.10 <= eta && eta < -0.00) * (9.40 <= pt * cosh(eta) && pt * cosh(eta) < 9.90) * (0.733208) +
    (-0.10 <= eta && eta < -0.00) * (9.90 <= pt * cosh(eta) && pt * cosh(eta) < 10.50) * (0.710395) +
    (-0.10 <= eta && eta < -0.00) * (10.50 <= pt * cosh(eta) && pt * cosh(eta) < 11.30) * (0.687940) +
    (-0.10 <= eta && eta < -0.00) * (11.30 <= pt * cosh(eta) && pt * cosh(eta) < 12.40) * (0.665437) +
    (-0.10 <= eta && eta < -0.00) * (12.40 <= pt * cosh(eta) && pt * cosh(eta) < 14.30) * (0.642619) +
    (-0.10 <= eta && eta < -0.00) * (14.30 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.605218) +
    (-0.00 <= eta && eta < 0.10) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) +
    (-0.00 <= eta && eta < 0.10) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 5.50) * (0.999487) +
    (-0.00 <= eta && eta < 0.10) * (5.50 <= pt * cosh(eta) && pt * cosh(eta) < 6.00) * (0.980086) +
    (-0.00 <= eta && eta < 0.10) * (6.00 <= pt * cosh(eta) && pt * cosh(eta) < 6.40) * (0.955878) +
    (-0.00 <= eta && eta < 0.10) * (6.40 <= pt * cosh(eta) && pt * cosh(eta) < 6.70) * (0.931602) +
    (-0.00 <= eta && eta < 0.10) * (6.70 <= pt * cosh(eta) && pt * cosh(eta) < 7.00) * (0.908797) +
    (-0.00 <= eta && eta < 0.10) * (7.00 <= pt * cosh(eta) && pt * cosh(eta) < 7.30) * (0.885676) +
    (-0.00 <= eta && eta < 0.10) * (7.30 <= pt * cosh(eta) && pt * cosh(eta) < 7.60) * (0.862280) +
    (-0.00 <= eta && eta < 0.10) * (7.60 <= pt * cosh(eta) && pt * cosh(eta) < 7.90) * (0.839159) +
    (-0.00 <= eta && eta < 0.10) * (7.90 <= pt * cosh(eta) && pt * cosh(eta) < 8.20) * (0.815846) +
    (-0.00 <= eta && eta < 0.10) * (8.20 <= pt * cosh(eta) && pt * cosh(eta) < 8.50) * (0.795129) +
    (-0.00 <= eta && eta < 0.10) * (8.50 <= pt * cosh(eta) && pt * cosh(eta) < 8.90) * (0.773860) +
    (-0.00 <= eta && eta < 0.10) * (8.90 <= pt * cosh(eta) && pt * cosh(eta) < 9.30) * (0.751980) +
    (-0.00 <= eta && eta < 0.10) * (9.30 <= pt * cosh(eta) && pt * cosh(eta) < 9.90) * (0.728013) +
    (-0.00 <= eta && eta < 0.10) * (9.90 <= pt * cosh(eta) && pt * cosh(eta) < 10.50) * (0.705518) +
    (-0.00 <= eta && eta < 0.10) * (10.50 <= pt * cosh(eta) && pt * cosh(eta) < 11.30) * (0.683919) +
    (-0.00 <= eta && eta < 0.10) * (11.30 <= pt * cosh(eta) && pt * cosh(eta) < 12.50) * (0.661467) +
    (-0.00 <= eta && eta < 0.10) * (12.50 <= pt * cosh(eta) && pt * cosh(eta) < 14.70) * (0.638022) +
    (-0.00 <= eta && eta < 0.10) * (14.70 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.604731) +
    (0.10 <= eta && eta < 0.20) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) +
    (0.10 <= eta && eta < 0.20) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 5.40) * (0.999390) +
    (0.10 <= eta && eta < 0.20) * (5.40 <= pt * cosh(eta) && pt * cosh(eta) < 5.90) * (0.978436) +
    (0.10 <= eta && eta < 0.20) * (5.90 <= pt * cosh(eta) && pt * cosh(eta) < 6.30) * (0.954026) +
    (0.10 <= eta && eta < 0.20) * (6.30 <= pt * cosh(eta) && pt * cosh(eta) < 6.60) * (0.928980) +
    (0.10 <= eta && eta < 0.20) * (6.60 <= pt * cosh(eta) && pt * cosh(eta) < 6.90) * (0.905261) +
    (0.10 <= eta && eta < 0.20) * (6.90 <= pt * cosh(eta) && pt * cosh(eta) < 7.20) * (0.880827) +
    (0.10 <= eta && eta < 0.20) * (7.20 <= pt * cosh(eta) && pt * cosh(eta) < 7.50) * (0.855798) +
    (0.10 <= eta && eta < 0.20) * (7.50 <= pt * cosh(eta) && pt * cosh(eta) < 7.80) * (0.831816) +
    (0.10 <= eta && eta < 0.20) * (7.80 <= pt * cosh(eta) && pt * cosh(eta) < 8.10) * (0.811017) +
    (0.10 <= eta && eta < 0.20) * (8.10 <= pt * cosh(eta) && pt * cosh(eta) < 8.40) * (0.790324) +
    (0.10 <= eta && eta < 0.20) * (8.40 <= pt * cosh(eta) && pt * cosh(eta) < 8.80) * (0.766845) +
    (0.10 <= eta && eta < 0.20) * (8.80 <= pt * cosh(eta) && pt * cosh(eta) < 9.30) * (0.742732) +
    (0.10 <= eta && eta < 0.20) * (9.30 <= pt * cosh(eta) && pt * cosh(eta) < 9.80) * (0.721120) +
    (0.10 <= eta && eta < 0.20) * (9.80 <= pt * cosh(eta) && pt * cosh(eta) < 10.50) * (0.698097) +
    (0.10 <= eta && eta < 0.20) * (10.50 <= pt * cosh(eta) && pt * cosh(eta) < 11.40) * (0.675216) +
    (0.10 <= eta && eta < 0.20) * (11.40 <= pt * cosh(eta) && pt * cosh(eta) < 12.80) * (0.652610) +
    (0.10 <= eta && eta < 0.20) * (12.80 <= pt * cosh(eta) && pt * cosh(eta) < 16.20) * (0.627707) +
    (0.10 <= eta && eta < 0.20) * (16.20 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.603627) +
    (0.20 <= eta && eta < 0.30) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) +
    (0.20 <= eta && eta < 0.30) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 5.30) * (0.999486) +
    (0.20 <= eta && eta < 0.30) * (5.30 <= pt * cosh(eta) && pt * cosh(eta) < 5.80) * (0.979519) +
    (0.20 <= eta && eta < 0.30) * (5.80 <= pt * cosh(eta) && pt * cosh(eta) < 6.20) * (0.955406) +
    (0.20 <= eta && eta < 0.30) * (6.20 <= pt * cosh(eta) && pt * cosh(eta) < 6.50) * (0.930947) +
    (0.20 <= eta && eta < 0.30) * (6.50 <= pt * cosh(eta) && pt * cosh(eta) < 6.80) * (0.906986) +
    (0.20 <= eta && eta < 0.30) * (6.80 <= pt * cosh(eta) && pt * cosh(eta) < 7.10) * (0.882189) +
    (0.20 <= eta && eta < 0.30) * (7.10 <= pt * cosh(eta) && pt * cosh(eta) < 7.40) * (0.857189) +
    (0.20 <= eta && eta < 0.30) * (7.40 <= pt * cosh(eta) && pt * cosh(eta) < 7.70) * (0.833268) +
    (0.20 <= eta && eta < 0.30) * (7.70 <= pt * cosh(eta) && pt * cosh(eta) < 8.00) * (0.811348) +
    (0.20 <= eta && eta < 0.30) * (8.00 <= pt * cosh(eta) && pt * cosh(eta) < 8.40) * (0.787352) +
    (0.20 <= eta && eta < 0.30) * (8.40 <= pt * cosh(eta) && pt * cosh(eta) < 8.80) * (0.762475) +
    (0.20 <= eta && eta < 0.30) * (8.80 <= pt * cosh(eta) && pt * cosh(eta) < 9.30) * (0.737585) +
    (0.20 <= eta && eta < 0.30) * (9.30 <= pt * cosh(eta) && pt * cosh(eta) < 9.80) * (0.716144) +
    (0.20 <= eta && eta < 0.30) * (9.80 <= pt * cosh(eta) && pt * cosh(eta) < 10.50) * (0.693021) +
    (0.20 <= eta && eta < 0.30) * (10.50 <= pt * cosh(eta) && pt * cosh(eta) < 11.50) * (0.670061) +
    (0.20 <= eta && eta < 0.30) * (11.50 <= pt * cosh(eta) && pt * cosh(eta) < 13.10) * (0.646891) +
    (0.20 <= eta && eta < 0.30) * (13.10 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.605262) +
    (0.30 <= eta && eta < 0.40) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) +
    (0.30 <= eta && eta < 0.40) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 5.30) * (0.999459) +
    (0.30 <= eta && eta < 0.40) * (5.30 <= pt * cosh(eta) && pt * cosh(eta) < 5.80) * (0.979305) +
    (0.30 <= eta && eta < 0.40) * (5.80 <= pt * cosh(eta) && pt * cosh(eta) < 6.20) * (0.955111) +
    (0.30 <= eta && eta < 0.40) * (6.20 <= pt * cosh(eta) && pt * cosh(eta) < 6.50) * (0.929350) +
    (0.30 <= eta && eta < 0.40) * (6.50 <= pt * cosh(eta) && pt * cosh(eta) < 6.80) * (0.904229) +
    (0.30 <= eta && eta < 0.40) * (6.80 <= pt * cosh(eta) && pt * cosh(eta) < 7.10) * (0.880114) +
    (0.30 <= eta && eta < 0.40) * (7.10 <= pt * cosh(eta) && pt * cosh(eta) < 7.40) * (0.855471) +
    (0.30 <= eta && eta < 0.40) * (7.40 <= pt * cosh(eta) && pt * cosh(eta) < 7.70) * (0.831596) +
    (0.30 <= eta && eta < 0.40) * (7.70 <= pt * cosh(eta) && pt * cosh(eta) < 8.00) * (0.808837) +
    (0.30 <= eta && eta < 0.40) * (8.00 <= pt * cosh(eta) && pt * cosh(eta) < 8.40) * (0.785228) +
    (0.30 <= eta && eta < 0.40) * (8.40 <= pt * cosh(eta) && pt * cosh(eta) < 8.90) * (0.756982) +
    (0.30 <= eta && eta < 0.40) * (8.90 <= pt * cosh(eta) && pt * cosh(eta) < 9.40) * (0.732645) +
    (0.30 <= eta && eta < 0.40) * (9.40 <= pt * cosh(eta) && pt * cosh(eta) < 10.00) * (0.709880) +
    (0.30 <= eta && eta < 0.40) * (10.00 <= pt * cosh(eta) && pt * cosh(eta) < 10.80) * (0.687128) +
    (0.30 <= eta && eta < 0.40) * (10.80 <= pt * cosh(eta) && pt * cosh(eta) < 11.90) * (0.663925) +
    (0.30 <= eta && eta < 0.40) * (11.90 <= pt * cosh(eta) && pt * cosh(eta) < 13.90) * (0.640354) +
    (0.30 <= eta && eta < 0.40) * (13.90 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.604740) +
    (0.40 <= eta && eta < 0.50) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) +
    (0.40 <= eta && eta < 0.50) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 5.30) * (0.999416) +
    (0.40 <= eta && eta < 0.50) * (5.30 <= pt * cosh(eta) && pt * cosh(eta) < 5.80) * (0.977866) +
    (0.40 <= eta && eta < 0.50) * (5.80 <= pt * cosh(eta) && pt * cosh(eta) < 6.10) * (0.956023) +
    (0.40 <= eta && eta < 0.50) * (6.10 <= pt * cosh(eta) && pt * cosh(eta) < 6.40) * (0.934804) +
    (0.40 <= eta && eta < 0.50) * (6.40 <= pt * cosh(eta) && pt * cosh(eta) < 6.70) * (0.910945) +
    (0.40 <= eta && eta < 0.50) * (6.70 <= pt * cosh(eta) && pt * cosh(eta) < 7.00) * (0.886157) +
    (0.40 <= eta && eta < 0.50) * (7.00 <= pt * cosh(eta) && pt * cosh(eta) < 7.30) * (0.861028) +
    (0.40 <= eta && eta < 0.50) * (7.30 <= pt * cosh(eta) && pt * cosh(eta) < 7.60) * (0.836587) +
    (0.40 <= eta && eta < 0.50) * (7.60 <= pt * cosh(eta) && pt * cosh(eta) < 7.90) * (0.813125) +
    (0.40 <= eta && eta < 0.50) * (7.90 <= pt * cosh(eta) && pt * cosh(eta) < 8.30) * (0.788762) +
    (0.40 <= eta && eta < 0.50) * (8.30 <= pt * cosh(eta) && pt * cosh(eta) < 8.70) * (0.764916) +
    (0.40 <= eta && eta < 0.50) * (8.70 <= pt * cosh(eta) && pt * cosh(eta) < 9.20) * (0.739594) +
    (0.40 <= eta && eta < 0.50) * (9.20 <= pt * cosh(eta) && pt * cosh(eta) < 9.80) * (0.715349) +
    (0.40 <= eta && eta < 0.50) * (9.80 <= pt * cosh(eta) && pt * cosh(eta) < 10.50) * (0.691540) +
    (0.40 <= eta && eta < 0.50) * (10.50 <= pt * cosh(eta) && pt * cosh(eta) < 11.50) * (0.668971) +
    (0.40 <= eta && eta < 0.50) * (11.50 <= pt * cosh(eta) && pt * cosh(eta) < 13.20) * (0.645517) +
    (0.40 <= eta && eta < 0.50) * (13.20 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.605109) +
    (0.50 <= eta && eta < 0.60) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) +
    (0.50 <= eta && eta < 0.60) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 5.40) * (0.999374) +
    (0.50 <= eta && eta < 0.60) * (5.40 <= pt * cosh(eta) && pt * cosh(eta) < 5.90) * (0.976949) +
    (0.50 <= eta && eta < 0.60) * (5.90 <= pt * cosh(eta) && pt * cosh(eta) < 6.30) * (0.951356) +
    (0.50 <= eta && eta < 0.60) * (6.30 <= pt * cosh(eta) && pt * cosh(eta) < 6.60) * (0.926120) +
    (0.50 <= eta && eta < 0.60) * (6.60 <= pt * cosh(eta) && pt * cosh(eta) < 6.90) * (0.902566) +
    (0.50 <= eta && eta < 0.60) * (6.90 <= pt * cosh(eta) && pt * cosh(eta) < 7.20) * (0.878135) +
    (0.50 <= eta && eta < 0.60) * (7.20 <= pt * cosh(eta) && pt * cosh(eta) < 7.50) * (0.853567) +
    (0.50 <= eta && eta < 0.60) * (7.50 <= pt * cosh(eta) && pt * cosh(eta) < 7.80) * (0.830007) +
    (0.50 <= eta && eta < 0.60) * (7.80 <= pt * cosh(eta) && pt * cosh(eta) < 8.10) * (0.807715) +
    (0.50 <= eta && eta < 0.60) * (8.10 <= pt * cosh(eta) && pt * cosh(eta) < 8.50) * (0.784361) +
    (0.50 <= eta && eta < 0.60) * (8.50 <= pt * cosh(eta) && pt * cosh(eta) < 8.90) * (0.760491) +
    (0.50 <= eta && eta < 0.60) * (8.90 <= pt * cosh(eta) && pt * cosh(eta) < 9.40) * (0.736736) +
    (0.50 <= eta && eta < 0.60) * (9.40 <= pt * cosh(eta) && pt * cosh(eta) < 10.00) * (0.712414) +
    (0.50 <= eta && eta < 0.60) * (10.00 <= pt * cosh(eta) && pt * cosh(eta) < 10.70) * (0.690695) +
    (0.50 <= eta && eta < 0.60) * (10.70 <= pt * cosh(eta) && pt * cosh(eta) < 11.70) * (0.668691) +
    (0.50 <= eta && eta < 0.60) * (11.70 <= pt * cosh(eta) && pt * cosh(eta) < 13.40) * (0.645647) +
    (0.50 <= eta && eta < 0.60) * (13.40 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.605221) +
    (0.60 <= eta && eta < 0.70) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) +
    (0.60 <= eta && eta < 0.70) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 5.40) * (0.999428) +
    (0.60 <= eta && eta < 0.70) * (5.40 <= pt * cosh(eta) && pt * cosh(eta) < 5.90) * (0.978093) +
    (0.60 <= eta && eta < 0.70) * (5.90 <= pt * cosh(eta) && pt * cosh(eta) < 6.30) * (0.952896) +
    (0.60 <= eta && eta < 0.70) * (6.30 <= pt * cosh(eta) && pt * cosh(eta) < 6.60) * (0.928743) +
    (0.60 <= eta && eta < 0.70) * (6.60 <= pt * cosh(eta) && pt * cosh(eta) < 6.90) * (0.904835) +
    (0.60 <= eta && eta < 0.70) * (6.90 <= pt * cosh(eta) && pt * cosh(eta) < 7.20) * (0.880239) +
    (0.60 <= eta && eta < 0.70) * (7.20 <= pt * cosh(eta) && pt * cosh(eta) < 7.50) * (0.855800) +
    (0.60 <= eta && eta < 0.70) * (7.50 <= pt * cosh(eta) && pt * cosh(eta) < 7.80) * (0.832129) +
    (0.60 <= eta && eta < 0.70) * (7.80 <= pt * cosh(eta) && pt * cosh(eta) < 8.10) * (0.809971) +
    (0.60 <= eta && eta < 0.70) * (8.10 <= pt * cosh(eta) && pt * cosh(eta) < 8.50) * (0.785533) +
    (0.60 <= eta && eta < 0.70) * (8.50 <= pt * cosh(eta) && pt * cosh(eta) < 8.90) * (0.762384) +
    (0.60 <= eta && eta < 0.70) * (8.90 <= pt * cosh(eta) && pt * cosh(eta) < 9.40) * (0.739140) +
    (0.60 <= eta && eta < 0.70) * (9.40 <= pt * cosh(eta) && pt * cosh(eta) < 10.00) * (0.715270) +
    (0.60 <= eta && eta < 0.70) * (10.00 <= pt * cosh(eta) && pt * cosh(eta) < 10.70) * (0.692367) +
    (0.60 <= eta && eta < 0.70) * (10.70 <= pt * cosh(eta) && pt * cosh(eta) < 11.70) * (0.669998) +
    (0.60 <= eta && eta < 0.70) * (11.70 <= pt * cosh(eta) && pt * cosh(eta) < 13.30) * (0.647168) +
    (0.60 <= eta && eta < 0.70) * (13.30 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.605392) +
    (0.70 <= eta && eta < 0.80) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) +
    (0.70 <= eta && eta < 0.80) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 5.40) * (0.999431) +
    (0.70 <= eta && eta < 0.80) * (5.40 <= pt * cosh(eta) && pt * cosh(eta) < 5.90) * (0.978854) +
    (0.70 <= eta && eta < 0.80) * (5.90 <= pt * cosh(eta) && pt * cosh(eta) < 6.20) * (0.957614) +
    (0.70 <= eta && eta < 0.80) * (6.20 <= pt * cosh(eta) && pt * cosh(eta) < 6.50) * (0.936746) +
    (0.70 <= eta && eta < 0.80) * (6.50 <= pt * cosh(eta) && pt * cosh(eta) < 6.80) * (0.914572) +
    (0.70 <= eta && eta < 0.80) * (6.80 <= pt * cosh(eta) && pt * cosh(eta) < 7.10) * (0.889412) +
    (0.70 <= eta && eta < 0.80) * (7.10 <= pt * cosh(eta) && pt * cosh(eta) < 7.40) * (0.863632) +
    (0.70 <= eta && eta < 0.80) * (7.40 <= pt * cosh(eta) && pt * cosh(eta) < 7.70) * (0.839744) +
    (0.70 <= eta && eta < 0.80) * (7.70 <= pt * cosh(eta) && pt * cosh(eta) < 8.10) * (0.814966) +
    (0.70 <= eta && eta < 0.80) * (8.10 <= pt * cosh(eta) && pt * cosh(eta) < 8.50) * (0.787803) +
    (0.70 <= eta && eta < 0.80) * (8.50 <= pt * cosh(eta) && pt * cosh(eta) < 8.90) * (0.762097) +
    (0.70 <= eta && eta < 0.80) * (8.90 <= pt * cosh(eta) && pt * cosh(eta) < 9.40) * (0.738550) +
    (0.70 <= eta && eta < 0.80) * (9.40 <= pt * cosh(eta) && pt * cosh(eta) < 10.00) * (0.715206) +
    (0.70 <= eta && eta < 0.80) * (10.00 <= pt * cosh(eta) && pt * cosh(eta) < 10.70) * (0.692151) +
    (0.70 <= eta && eta < 0.80) * (10.70 <= pt * cosh(eta) && pt * cosh(eta) < 11.70) * (0.669829) +
    (0.70 <= eta && eta < 0.80) * (11.70 <= pt * cosh(eta) && pt * cosh(eta) < 13.30) * (0.647051) +
    (0.70 <= eta && eta < 0.80) * (13.30 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.605381) +
    (0.80 <= eta && eta < 0.90) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) +
    (0.80 <= eta && eta < 0.90) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 5.50) * (0.999463) +
    (0.80 <= eta && eta < 0.90) * (5.50 <= pt * cosh(eta) && pt * cosh(eta) < 6.00) * (0.979480) +
    (0.80 <= eta && eta < 0.90) * (6.00 <= pt * cosh(eta) && pt * cosh(eta) < 6.40) * (0.955658) +
    (0.80 <= eta && eta < 0.90) * (6.40 <= pt * cosh(eta) && pt * cosh(eta) < 6.70) * (0.931795) +
    (0.80 <= eta && eta < 0.90) * (6.70 <= pt * cosh(eta) && pt * cosh(eta) < 7.00) * (0.909021) +
    (0.80 <= eta && eta < 0.90) * (7.00 <= pt * cosh(eta) && pt * cosh(eta) < 7.30) * (0.884511) +
    (0.80 <= eta && eta < 0.90) * (7.30 <= pt * cosh(eta) && pt * cosh(eta) < 7.60) * (0.859876) +
    (0.80 <= eta && eta < 0.90) * (7.60 <= pt * cosh(eta) && pt * cosh(eta) < 7.90) * (0.837206) +
    (0.80 <= eta && eta < 0.90) * (7.90 <= pt * cosh(eta) && pt * cosh(eta) < 8.20) * (0.815556) +
    (0.80 <= eta && eta < 0.90) * (8.20 <= pt * cosh(eta) && pt * cosh(eta) < 8.60) * (0.792098) +
    (0.80 <= eta && eta < 0.90) * (8.60 <= pt * cosh(eta) && pt * cosh(eta) < 9.00) * (0.767572) +
    (0.80 <= eta && eta < 0.90) * (9.00 <= pt * cosh(eta) && pt * cosh(eta) < 9.50) * (0.742973) +
    (0.80 <= eta && eta < 0.90) * (9.50 <= pt * cosh(eta) && pt * cosh(eta) < 10.10) * (0.719838) +
    (0.80 <= eta && eta < 0.90) * (10.10 <= pt * cosh(eta) && pt * cosh(eta) < 10.80) * (0.696957) +
    (0.80 <= eta && eta < 0.90) * (10.80 <= pt * cosh(eta) && pt * cosh(eta) < 11.70) * (0.674906) +
    (0.80 <= eta && eta < 0.90) * (11.70 <= pt * cosh(eta) && pt * cosh(eta) < 13.10) * (0.652868) +
    (0.80 <= eta && eta < 0.90) * (13.10 <= pt * cosh(eta) && pt * cosh(eta) < 16.50) * (0.628233) +
    (0.80 <= eta && eta < 0.90) * (16.50 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.603774) +
    (0.90 <= eta && eta < 1.00) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) +
    (0.90 <= eta && eta < 1.00) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 5.50) * (0.999490) +
    (0.90 <= eta && eta < 1.00) * (5.50 <= pt * cosh(eta) && pt * cosh(eta) < 6.00) * (0.980505) +
    (0.90 <= eta && eta < 1.00) * (6.00 <= pt * cosh(eta) && pt * cosh(eta) < 6.40) * (0.957300) +
    (0.90 <= eta && eta < 1.00) * (6.40 <= pt * cosh(eta) && pt * cosh(eta) < 6.70) * (0.933891) +
    (0.90 <= eta && eta < 1.00) * (6.70 <= pt * cosh(eta) && pt * cosh(eta) < 7.00) * (0.910870) +
    (0.90 <= eta && eta < 1.00) * (7.00 <= pt * cosh(eta) && pt * cosh(eta) < 7.30) * (0.887137) +
    (0.90 <= eta && eta < 1.00) * (7.30 <= pt * cosh(eta) && pt * cosh(eta) < 7.60) * (0.863393) +
    (0.90 <= eta && eta < 1.00) * (7.60 <= pt * cosh(eta) && pt * cosh(eta) < 7.90) * (0.840634) +
    (0.90 <= eta && eta < 1.00) * (7.90 <= pt * cosh(eta) && pt * cosh(eta) < 8.20) * (0.818273) +
    (0.90 <= eta && eta < 1.00) * (8.20 <= pt * cosh(eta) && pt * cosh(eta) < 8.60) * (0.794422) +
    (0.90 <= eta && eta < 1.00) * (8.60 <= pt * cosh(eta) && pt * cosh(eta) < 9.00) * (0.769639) +
    (0.90 <= eta && eta < 1.00) * (9.00 <= pt * cosh(eta) && pt * cosh(eta) < 9.40) * (0.748157) +
    (0.90 <= eta && eta < 1.00) * (9.40 <= pt * cosh(eta) && pt * cosh(eta) < 10.00) * (0.724553) +
    (0.90 <= eta && eta < 1.00) * (10.00 <= pt * cosh(eta) && pt * cosh(eta) < 10.70) * (0.699984) +
    (0.90 <= eta && eta < 1.00) * (10.70 <= pt * cosh(eta) && pt * cosh(eta) < 11.60) * (0.677138) +
    (0.90 <= eta && eta < 1.00) * (11.60 <= pt * cosh(eta) && pt * cosh(eta) < 13.00) * (0.654326) +
    (0.90 <= eta && eta < 1.00) * (13.00 <= pt * cosh(eta) && pt * cosh(eta) < 16.20) * (0.629520) +
    (0.90 <= eta && eta < 1.00) * (16.20 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.603898) +
    (1.00 <= eta && eta < 1.10) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.000000) }


  # --- protons ---

  add EfficiencyFormula {2212} {321} {
    (eta<-1.00 || eta>=1.00 || pt * cosh(eta) < 0.10 || pt * cosh(eta) >= 50.00) * ( 0.00 ) +
    (-1.00 <= eta && eta < -0.90) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 10.30) * (0.000705) +
    (-1.00 <= eta && eta < -0.90) * (10.30 <= pt * cosh(eta) && pt * cosh(eta) < 11.40) * (0.021025) +
    (-1.00 <= eta && eta < -0.90) * (11.40 <= pt * cosh(eta) && pt * cosh(eta) < 12.20) * (0.042481) +
    (-1.00 <= eta && eta < -0.90) * (12.20 <= pt * cosh(eta) && pt * cosh(eta) < 12.90) * (0.063884) +
    (-1.00 <= eta && eta < -0.90) * (12.90 <= pt * cosh(eta) && pt * cosh(eta) < 13.60) * (0.085969) +
    (-1.00 <= eta && eta < -0.90) * (13.60 <= pt * cosh(eta) && pt * cosh(eta) < 14.30) * (0.108965) +
    (-1.00 <= eta && eta < -0.90) * (14.30 <= pt * cosh(eta) && pt * cosh(eta) < 15.00) * (0.132030) +
    (-1.00 <= eta && eta < -0.90) * (15.00 <= pt * cosh(eta) && pt * cosh(eta) < 15.70) * (0.154564) +
    (-1.00 <= eta && eta < -0.90) * (15.70 <= pt * cosh(eta) && pt * cosh(eta) < 16.40) * (0.176174) +
    (-1.00 <= eta && eta < -0.90) * (16.40 <= pt * cosh(eta) && pt * cosh(eta) < 17.20) * (0.198024) +
    (-1.00 <= eta && eta < -0.90) * (17.20 <= pt * cosh(eta) && pt * cosh(eta) < 18.00) * (0.219729) +
    (-1.00 <= eta && eta < -0.90) * (18.00 <= pt * cosh(eta) && pt * cosh(eta) < 18.90) * (0.240906) +
    (-1.00 <= eta && eta < -0.90) * (18.90 <= pt * cosh(eta) && pt * cosh(eta) < 19.90) * (0.262363) +
    (-1.00 <= eta && eta < -0.90) * (19.90 <= pt * cosh(eta) && pt * cosh(eta) < 21.00) * (0.283549) +
    (-1.00 <= eta && eta < -0.90) * (21.00 <= pt * cosh(eta) && pt * cosh(eta) < 22.30) * (0.304817) +
    (-1.00 <= eta && eta < -0.90) * (22.30 <= pt * cosh(eta) && pt * cosh(eta) < 23.80) * (0.326209) +
    (-1.00 <= eta && eta < -0.90) * (23.80 <= pt * cosh(eta) && pt * cosh(eta) < 25.60) * (0.347442) +
    (-1.00 <= eta && eta < -0.90) * (25.60 <= pt * cosh(eta) && pt * cosh(eta) < 27.80) * (0.368549) +
    (-1.00 <= eta && eta < -0.90) * (27.80 <= pt * cosh(eta) && pt * cosh(eta) < 30.70) * (0.389801) +
    (-1.00 <= eta && eta < -0.90) * (30.70 <= pt * cosh(eta) && pt * cosh(eta) < 34.70) * (0.411336) +
    (-1.00 <= eta && eta < -0.90) * (34.70 <= pt * cosh(eta) && pt * cosh(eta) < 40.90) * (0.433234) +
    (-1.00 <= eta && eta < -0.90) * (40.90 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.453570) +
    (-0.90 <= eta && eta < -0.80) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 10.10) * (0.000649) +
    (-0.90 <= eta && eta < -0.80) * (10.10 <= pt * cosh(eta) && pt * cosh(eta) < 11.20) * (0.020231) +
    (-0.90 <= eta && eta < -0.80) * (11.20 <= pt * cosh(eta) && pt * cosh(eta) < 12.00) * (0.041709) +
    (-0.90 <= eta && eta < -0.80) * (12.00 <= pt * cosh(eta) && pt * cosh(eta) < 12.70) * (0.063327) +
    (-0.90 <= eta && eta < -0.80) * (12.70 <= pt * cosh(eta) && pt * cosh(eta) < 13.40) * (0.085708) +
    (-0.90 <= eta && eta < -0.80) * (13.40 <= pt * cosh(eta) && pt * cosh(eta) < 14.00) * (0.107361) +
    (-0.90 <= eta && eta < -0.80) * (14.00 <= pt * cosh(eta) && pt * cosh(eta) < 14.70) * (0.129120) +
    (-0.90 <= eta && eta < -0.80) * (14.70 <= pt * cosh(eta) && pt * cosh(eta) < 15.40) * (0.152075) +
    (-0.90 <= eta && eta < -0.80) * (15.40 <= pt * cosh(eta) && pt * cosh(eta) < 16.10) * (0.174114) +
    (-0.90 <= eta && eta < -0.80) * (16.10 <= pt * cosh(eta) && pt * cosh(eta) < 16.90) * (0.196407) +
    (-0.90 <= eta && eta < -0.80) * (16.90 <= pt * cosh(eta) && pt * cosh(eta) < 17.70) * (0.218543) +
    (-0.90 <= eta && eta < -0.80) * (17.70 <= pt * cosh(eta) && pt * cosh(eta) < 18.60) * (0.240121) +
    (-0.90 <= eta && eta < -0.80) * (18.60 <= pt * cosh(eta) && pt * cosh(eta) < 19.60) * (0.261958) +
    (-0.90 <= eta && eta < -0.80) * (19.60 <= pt * cosh(eta) && pt * cosh(eta) < 20.70) * (0.283483) +
    (-0.90 <= eta && eta < -0.80) * (20.70 <= pt * cosh(eta) && pt * cosh(eta) < 22.00) * (0.305052) +
    (-0.90 <= eta && eta < -0.80) * (22.00 <= pt * cosh(eta) && pt * cosh(eta) < 23.50) * (0.326699) +
    (-0.90 <= eta && eta < -0.80) * (23.50 <= pt * cosh(eta) && pt * cosh(eta) < 25.30) * (0.348135) +
    (-0.90 <= eta && eta < -0.80) * (25.30 <= pt * cosh(eta) && pt * cosh(eta) < 27.50) * (0.369386) +
    (-0.90 <= eta && eta < -0.80) * (27.50 <= pt * cosh(eta) && pt * cosh(eta) < 30.40) * (0.390721) +
    (-0.90 <= eta && eta < -0.80) * (30.40 <= pt * cosh(eta) && pt * cosh(eta) < 34.40) * (0.412270) +
    (-0.90 <= eta && eta < -0.80) * (34.40 <= pt * cosh(eta) && pt * cosh(eta) < 40.70) * (0.434265) +
    (-0.90 <= eta && eta < -0.80) * (40.70 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.454684) +
    (-0.80 <= eta && eta < -0.70) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 10.00) * (0.000698) +
    (-0.80 <= eta && eta < -0.70) * (10.00 <= pt * cosh(eta) && pt * cosh(eta) < 11.00) * (0.020156) +
    (-0.80 <= eta && eta < -0.70) * (11.00 <= pt * cosh(eta) && pt * cosh(eta) < 11.80) * (0.040719) +
    (-0.80 <= eta && eta < -0.70) * (11.80 <= pt * cosh(eta) && pt * cosh(eta) < 12.50) * (0.062503) +
    (-0.80 <= eta && eta < -0.70) * (12.50 <= pt * cosh(eta) && pt * cosh(eta) < 13.10) * (0.083481) +
    (-0.80 <= eta && eta < -0.70) * (13.10 <= pt * cosh(eta) && pt * cosh(eta) < 13.70) * (0.103705) +
    (-0.80 <= eta && eta < -0.70) * (13.70 <= pt * cosh(eta) && pt * cosh(eta) < 14.30) * (0.124110) +
    (-0.80 <= eta && eta < -0.70) * (14.30 <= pt * cosh(eta) && pt * cosh(eta) < 15.00) * (0.145883) +
    (-0.80 <= eta && eta < -0.70) * (15.00 <= pt * cosh(eta) && pt * cosh(eta) < 15.70) * (0.168517) +
    (-0.80 <= eta && eta < -0.70) * (15.70 <= pt * cosh(eta) && pt * cosh(eta) < 16.40) * (0.189997) +
    (-0.80 <= eta && eta < -0.70) * (16.40 <= pt * cosh(eta) && pt * cosh(eta) < 17.20) * (0.211530) +
    (-0.80 <= eta && eta < -0.70) * (17.20 <= pt * cosh(eta) && pt * cosh(eta) < 18.10) * (0.234017) +
    (-0.80 <= eta && eta < -0.70) * (18.10 <= pt * cosh(eta) && pt * cosh(eta) < 19.00) * (0.255655) +
    (-0.80 <= eta && eta < -0.70) * (19.00 <= pt * cosh(eta) && pt * cosh(eta) < 20.10) * (0.277170) +
    (-0.80 <= eta && eta < -0.70) * (20.10 <= pt * cosh(eta) && pt * cosh(eta) < 21.30) * (0.298995) +
    (-0.80 <= eta && eta < -0.70) * (21.30 <= pt * cosh(eta) && pt * cosh(eta) < 22.70) * (0.320333) +
    (-0.80 <= eta && eta < -0.70) * (22.70 <= pt * cosh(eta) && pt * cosh(eta) < 24.40) * (0.341878) +
    (-0.80 <= eta && eta < -0.70) * (24.40 <= pt * cosh(eta) && pt * cosh(eta) < 26.50) * (0.363611) +
    (-0.80 <= eta && eta < -0.70) * (26.50 <= pt * cosh(eta) && pt * cosh(eta) < 29.10) * (0.384994) +
    (-0.80 <= eta && eta < -0.70) * (29.10 <= pt * cosh(eta) && pt * cosh(eta) < 32.70) * (0.406381) +
    (-0.80 <= eta && eta < -0.70) * (32.70 <= pt * cosh(eta) && pt * cosh(eta) < 38.10) * (0.428249) +
    (-0.80 <= eta && eta < -0.70) * (38.10 <= pt * cosh(eta) && pt * cosh(eta) < 47.90) * (0.450901) +
    (-0.80 <= eta && eta < -0.70) * (47.90 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.462545) +
    (-0.70 <= eta && eta < -0.60) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 9.90) * (0.000644) +
    (-0.70 <= eta && eta < -0.60) * (9.90 <= pt * cosh(eta) && pt * cosh(eta) < 10.90) * (0.019416) +
    (-0.70 <= eta && eta < -0.60) * (10.90 <= pt * cosh(eta) && pt * cosh(eta) < 11.70) * (0.039800) +
    (-0.70 <= eta && eta < -0.60) * (11.70 <= pt * cosh(eta) && pt * cosh(eta) < 12.40) * (0.061565) +
    (-0.70 <= eta && eta < -0.60) * (12.40 <= pt * cosh(eta) && pt * cosh(eta) < 13.00) * (0.082602) +
    (-0.70 <= eta && eta < -0.60) * (13.00 <= pt * cosh(eta) && pt * cosh(eta) < 13.60) * (0.102923) +
    (-0.70 <= eta && eta < -0.60) * (13.60 <= pt * cosh(eta) && pt * cosh(eta) < 14.20) * (0.123448) +
    (-0.70 <= eta && eta < -0.60) * (14.20 <= pt * cosh(eta) && pt * cosh(eta) < 14.90) * (0.145363) +
    (-0.70 <= eta && eta < -0.60) * (14.90 <= pt * cosh(eta) && pt * cosh(eta) < 15.60) * (0.168148) +
    (-0.70 <= eta && eta < -0.60) * (15.60 <= pt * cosh(eta) && pt * cosh(eta) < 16.30) * (0.189770) +
    (-0.70 <= eta && eta < -0.60) * (16.30 <= pt * cosh(eta) && pt * cosh(eta) < 17.10) * (0.211440) +
    (-0.70 <= eta && eta < -0.60) * (17.10 <= pt * cosh(eta) && pt * cosh(eta) < 17.90) * (0.232802) +
    (-0.70 <= eta && eta < -0.60) * (17.90 <= pt * cosh(eta) && pt * cosh(eta) < 18.80) * (0.253505) +
    (-0.70 <= eta && eta < -0.60) * (18.80 <= pt * cosh(eta) && pt * cosh(eta) < 19.80) * (0.274359) +
    (-0.70 <= eta && eta < -0.60) * (19.80 <= pt * cosh(eta) && pt * cosh(eta) < 21.00) * (0.295736) +
    (-0.70 <= eta && eta < -0.60) * (21.00 <= pt * cosh(eta) && pt * cosh(eta) < 22.40) * (0.317652) +
    (-0.70 <= eta && eta < -0.60) * (22.40 <= pt * cosh(eta) && pt * cosh(eta) < 24.00) * (0.339116) +
    (-0.70 <= eta && eta < -0.60) * (24.00 <= pt * cosh(eta) && pt * cosh(eta) < 26.00) * (0.360427) +
    (-0.70 <= eta && eta < -0.60) * (26.00 <= pt * cosh(eta) && pt * cosh(eta) < 28.50) * (0.381782) +
    (-0.70 <= eta && eta < -0.60) * (28.50 <= pt * cosh(eta) && pt * cosh(eta) < 31.90) * (0.403198) +
    (-0.70 <= eta && eta < -0.60) * (31.90 <= pt * cosh(eta) && pt * cosh(eta) < 36.90) * (0.424964) +
    (-0.70 <= eta && eta < -0.60) * (36.90 <= pt * cosh(eta) && pt * cosh(eta) < 45.70) * (0.447490) +
    (-0.70 <= eta && eta < -0.60) * (45.70 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.461199) +
    (-0.60 <= eta && eta < -0.50) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 9.90) * (0.000623) +
    (-0.60 <= eta && eta < -0.50) * (9.90 <= pt * cosh(eta) && pt * cosh(eta) < 10.90) * (0.019050) +
    (-0.60 <= eta && eta < -0.50) * (10.90 <= pt * cosh(eta) && pt * cosh(eta) < 11.70) * (0.039230) +
    (-0.60 <= eta && eta < -0.50) * (11.70 <= pt * cosh(eta) && pt * cosh(eta) < 12.40) * (0.060853) +
    (-0.60 <= eta && eta < -0.50) * (12.40 <= pt * cosh(eta) && pt * cosh(eta) < 13.00) * (0.081798) +
    (-0.60 <= eta && eta < -0.50) * (13.00 <= pt * cosh(eta) && pt * cosh(eta) < 13.60) * (0.102059) +
    (-0.60 <= eta && eta < -0.50) * (13.60 <= pt * cosh(eta) && pt * cosh(eta) < 14.20) * (0.122547) +
    (-0.60 <= eta && eta < -0.50) * (14.20 <= pt * cosh(eta) && pt * cosh(eta) < 14.90) * (0.144443) +
    (-0.60 <= eta && eta < -0.50) * (14.90 <= pt * cosh(eta) && pt * cosh(eta) < 15.60) * (0.167226) +
    (-0.60 <= eta && eta < -0.50) * (15.60 <= pt * cosh(eta) && pt * cosh(eta) < 16.30) * (0.188860) +
    (-0.60 <= eta && eta < -0.50) * (16.30 <= pt * cosh(eta) && pt * cosh(eta) < 17.10) * (0.210554) +
    (-0.60 <= eta && eta < -0.50) * (17.10 <= pt * cosh(eta) && pt * cosh(eta) < 17.90) * (0.231950) +
    (-0.60 <= eta && eta < -0.50) * (17.90 <= pt * cosh(eta) && pt * cosh(eta) < 18.80) * (0.252693) +
    (-0.60 <= eta && eta < -0.50) * (18.80 <= pt * cosh(eta) && pt * cosh(eta) < 19.80) * (0.273595) +
    (-0.60 <= eta && eta < -0.50) * (19.80 <= pt * cosh(eta) && pt * cosh(eta) < 21.00) * (0.295028) +
    (-0.60 <= eta && eta < -0.50) * (21.00 <= pt * cosh(eta) && pt * cosh(eta) < 22.40) * (0.317006) +
    (-0.60 <= eta && eta < -0.50) * (22.40 <= pt * cosh(eta) && pt * cosh(eta) < 24.00) * (0.338536) +
    (-0.60 <= eta && eta < -0.50) * (24.00 <= pt * cosh(eta) && pt * cosh(eta) < 26.00) * (0.359915) +
    (-0.60 <= eta && eta < -0.50) * (26.00 <= pt * cosh(eta) && pt * cosh(eta) < 28.50) * (0.381344) +
    (-0.60 <= eta && eta < -0.50) * (28.50 <= pt * cosh(eta) && pt * cosh(eta) < 31.90) * (0.402835) +
    (-0.60 <= eta && eta < -0.50) * (31.90 <= pt * cosh(eta) && pt * cosh(eta) < 36.90) * (0.424681) +
    (-0.60 <= eta && eta < -0.50) * (36.90 <= pt * cosh(eta) && pt * cosh(eta) < 45.60) * (0.447175) +
    (-0.60 <= eta && eta < -0.50) * (45.60 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.460966) +
    (-0.50 <= eta && eta < -0.40) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 9.80) * (0.000665) +
    (-0.50 <= eta && eta < -0.40) * (9.80 <= pt * cosh(eta) && pt * cosh(eta) < 10.80) * (0.020045) +
    (-0.50 <= eta && eta < -0.40) * (10.80 <= pt * cosh(eta) && pt * cosh(eta) < 11.60) * (0.040989) +
    (-0.50 <= eta && eta < -0.40) * (11.60 <= pt * cosh(eta) && pt * cosh(eta) < 12.30) * (0.063254) +
    (-0.50 <= eta && eta < -0.40) * (12.30 <= pt * cosh(eta) && pt * cosh(eta) < 12.90) * (0.084688) +
    (-0.50 <= eta && eta < -0.40) * (12.90 <= pt * cosh(eta) && pt * cosh(eta) < 13.50) * (0.105322) +
    (-0.50 <= eta && eta < -0.40) * (13.50 <= pt * cosh(eta) && pt * cosh(eta) < 14.10) * (0.126101) +
    (-0.50 <= eta && eta < -0.40) * (14.10 <= pt * cosh(eta) && pt * cosh(eta) < 14.70) * (0.146555) +
    (-0.50 <= eta && eta < -0.40) * (14.70 <= pt * cosh(eta) && pt * cosh(eta) < 15.40) * (0.167958) +
    (-0.50 <= eta && eta < -0.40) * (15.40 <= pt * cosh(eta) && pt * cosh(eta) < 16.10) * (0.189862) +
    (-0.50 <= eta && eta < -0.40) * (16.10 <= pt * cosh(eta) && pt * cosh(eta) < 16.90) * (0.211794) +
    (-0.50 <= eta && eta < -0.40) * (16.90 <= pt * cosh(eta) && pt * cosh(eta) < 17.70) * (0.233391) +
    (-0.50 <= eta && eta < -0.40) * (17.70 <= pt * cosh(eta) && pt * cosh(eta) < 18.60) * (0.254294) +
    (-0.50 <= eta && eta < -0.40) * (18.60 <= pt * cosh(eta) && pt * cosh(eta) < 19.60) * (0.275318) +
    (-0.50 <= eta && eta < -0.40) * (19.60 <= pt * cosh(eta) && pt * cosh(eta) < 20.80) * (0.296835) +
    (-0.50 <= eta && eta < -0.40) * (20.80 <= pt * cosh(eta) && pt * cosh(eta) < 22.20) * (0.318854) +
    (-0.50 <= eta && eta < -0.40) * (22.20 <= pt * cosh(eta) && pt * cosh(eta) < 23.80) * (0.340377) +
    (-0.50 <= eta && eta < -0.40) * (23.80 <= pt * cosh(eta) && pt * cosh(eta) < 25.80) * (0.361703) +
    (-0.50 <= eta && eta < -0.40) * (25.80 <= pt * cosh(eta) && pt * cosh(eta) < 28.30) * (0.383028) +
    (-0.50 <= eta && eta < -0.40) * (28.30 <= pt * cosh(eta) && pt * cosh(eta) < 31.70) * (0.404363) +
    (-0.50 <= eta && eta < -0.40) * (31.70 <= pt * cosh(eta) && pt * cosh(eta) < 36.80) * (0.426192) +
    (-0.50 <= eta && eta < -0.40) * (36.80 <= pt * cosh(eta) && pt * cosh(eta) < 45.80) * (0.448789) +
    (-0.50 <= eta && eta < -0.40) * (45.80 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.462263) +
    (-0.40 <= eta && eta < -0.30) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 9.80) * (0.000645) +
    (-0.40 <= eta && eta < -0.30) * (9.80 <= pt * cosh(eta) && pt * cosh(eta) < 10.80) * (0.019306) +
    (-0.40 <= eta && eta < -0.30) * (10.80 <= pt * cosh(eta) && pt * cosh(eta) < 11.60) * (0.039824) +
    (-0.40 <= eta && eta < -0.30) * (11.60 <= pt * cosh(eta) && pt * cosh(eta) < 12.30) * (0.061801) +
    (-0.40 <= eta && eta < -0.30) * (12.30 <= pt * cosh(eta) && pt * cosh(eta) < 12.90) * (0.083051) +
    (-0.40 <= eta && eta < -0.30) * (12.90 <= pt * cosh(eta) && pt * cosh(eta) < 13.50) * (0.103567) +
    (-0.40 <= eta && eta < -0.30) * (13.50 <= pt * cosh(eta) && pt * cosh(eta) < 14.10) * (0.124275) +
    (-0.40 <= eta && eta < -0.30) * (14.10 <= pt * cosh(eta) && pt * cosh(eta) < 14.70) * (0.144695) +
    (-0.40 <= eta && eta < -0.30) * (14.70 <= pt * cosh(eta) && pt * cosh(eta) < 15.40) * (0.166094) +
    (-0.40 <= eta && eta < -0.30) * (15.40 <= pt * cosh(eta) && pt * cosh(eta) < 16.10) * (0.188023) +
    (-0.40 <= eta && eta < -0.30) * (16.10 <= pt * cosh(eta) && pt * cosh(eta) < 16.90) * (0.210005) +
    (-0.40 <= eta && eta < -0.30) * (16.90 <= pt * cosh(eta) && pt * cosh(eta) < 17.70) * (0.231670) +
    (-0.40 <= eta && eta < -0.30) * (17.70 <= pt * cosh(eta) && pt * cosh(eta) < 18.60) * (0.252656) +
    (-0.40 <= eta && eta < -0.30) * (18.60 <= pt * cosh(eta) && pt * cosh(eta) < 19.60) * (0.273778) +
    (-0.40 <= eta && eta < -0.30) * (19.60 <= pt * cosh(eta) && pt * cosh(eta) < 20.80) * (0.295407) +
    (-0.40 <= eta && eta < -0.30) * (20.80 <= pt * cosh(eta) && pt * cosh(eta) < 22.10) * (0.316782) +
    (-0.40 <= eta && eta < -0.30) * (22.10 <= pt * cosh(eta) && pt * cosh(eta) < 23.70) * (0.337879) +
    (-0.40 <= eta && eta < -0.30) * (23.70 <= pt * cosh(eta) && pt * cosh(eta) < 25.60) * (0.359062) +
    (-0.40 <= eta && eta < -0.30) * (25.60 <= pt * cosh(eta) && pt * cosh(eta) < 28.10) * (0.380434) +
    (-0.40 <= eta && eta < -0.30) * (28.10 <= pt * cosh(eta) && pt * cosh(eta) < 31.40) * (0.402053) +
    (-0.40 <= eta && eta < -0.30) * (31.40 <= pt * cosh(eta) && pt * cosh(eta) < 36.20) * (0.423693) +
    (-0.40 <= eta && eta < -0.30) * (36.20 <= pt * cosh(eta) && pt * cosh(eta) < 44.50) * (0.446017) +
    (-0.40 <= eta && eta < -0.30) * (44.50 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.460865) +
    (-0.30 <= eta && eta < -0.20) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 9.80) * (0.000624) +
    (-0.30 <= eta && eta < -0.20) * (9.80 <= pt * cosh(eta) && pt * cosh(eta) < 10.80) * (0.019209) +
    (-0.30 <= eta && eta < -0.20) * (10.80 <= pt * cosh(eta) && pt * cosh(eta) < 11.60) * (0.039736) +
    (-0.30 <= eta && eta < -0.20) * (11.60 <= pt * cosh(eta) && pt * cosh(eta) < 12.30) * (0.061691) +
    (-0.30 <= eta && eta < -0.20) * (12.30 <= pt * cosh(eta) && pt * cosh(eta) < 12.90) * (0.082927) +
    (-0.30 <= eta && eta < -0.20) * (12.90 <= pt * cosh(eta) && pt * cosh(eta) < 13.50) * (0.103434) +
    (-0.30 <= eta && eta < -0.20) * (13.50 <= pt * cosh(eta) && pt * cosh(eta) < 14.10) * (0.124136) +
    (-0.30 <= eta && eta < -0.20) * (14.10 <= pt * cosh(eta) && pt * cosh(eta) < 14.70) * (0.144554) +
    (-0.30 <= eta && eta < -0.20) * (14.70 <= pt * cosh(eta) && pt * cosh(eta) < 15.40) * (0.165953) +
    (-0.30 <= eta && eta < -0.20) * (15.40 <= pt * cosh(eta) && pt * cosh(eta) < 16.10) * (0.187883) +
    (-0.30 <= eta && eta < -0.20) * (16.10 <= pt * cosh(eta) && pt * cosh(eta) < 16.90) * (0.209869) +
    (-0.30 <= eta && eta < -0.20) * (16.90 <= pt * cosh(eta) && pt * cosh(eta) < 17.70) * (0.231539) +
    (-0.30 <= eta && eta < -0.20) * (17.70 <= pt * cosh(eta) && pt * cosh(eta) < 18.60) * (0.252531) +
    (-0.30 <= eta && eta < -0.20) * (18.60 <= pt * cosh(eta) && pt * cosh(eta) < 19.60) * (0.273661) +
    (-0.30 <= eta && eta < -0.20) * (19.60 <= pt * cosh(eta) && pt * cosh(eta) < 20.80) * (0.295298) +
    (-0.30 <= eta && eta < -0.20) * (20.80 <= pt * cosh(eta) && pt * cosh(eta) < 22.10) * (0.316682) +
    (-0.30 <= eta && eta < -0.20) * (22.10 <= pt * cosh(eta) && pt * cosh(eta) < 23.70) * (0.337790) +
    (-0.30 <= eta && eta < -0.20) * (23.70 <= pt * cosh(eta) && pt * cosh(eta) < 25.60) * (0.358983) +
    (-0.30 <= eta && eta < -0.20) * (25.60 <= pt * cosh(eta) && pt * cosh(eta) < 28.10) * (0.380366) +
    (-0.30 <= eta && eta < -0.20) * (28.10 <= pt * cosh(eta) && pt * cosh(eta) < 31.40) * (0.401997) +
    (-0.30 <= eta && eta < -0.20) * (31.40 <= pt * cosh(eta) && pt * cosh(eta) < 36.20) * (0.423648) +
    (-0.30 <= eta && eta < -0.20) * (36.20 <= pt * cosh(eta) && pt * cosh(eta) < 44.50) * (0.445986) +
    (-0.30 <= eta && eta < -0.20) * (44.50 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.460842) +
    (-0.20 <= eta && eta < -0.10) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 9.90) * (0.000667) +
    (-0.20 <= eta && eta < -0.10) * (9.90 <= pt * cosh(eta) && pt * cosh(eta) < 10.90) * (0.019658) +
    (-0.20 <= eta && eta < -0.10) * (10.90 <= pt * cosh(eta) && pt * cosh(eta) < 11.70) * (0.040176) +
    (-0.20 <= eta && eta < -0.10) * (11.70 <= pt * cosh(eta) && pt * cosh(eta) < 12.40) * (0.062033) +
    (-0.20 <= eta && eta < -0.10) * (12.40 <= pt * cosh(eta) && pt * cosh(eta) < 13.00) * (0.083131) +
    (-0.20 <= eta && eta < -0.10) * (13.00 <= pt * cosh(eta) && pt * cosh(eta) < 13.60) * (0.103490) +
    (-0.20 <= eta && eta < -0.10) * (13.60 <= pt * cosh(eta) && pt * cosh(eta) < 14.20) * (0.124039) +
    (-0.20 <= eta && eta < -0.10) * (14.20 <= pt * cosh(eta) && pt * cosh(eta) < 14.90) * (0.145966) +
    (-0.20 <= eta && eta < -0.10) * (14.90 <= pt * cosh(eta) && pt * cosh(eta) < 15.60) * (0.168751) +
    (-0.20 <= eta && eta < -0.10) * (15.60 <= pt * cosh(eta) && pt * cosh(eta) < 16.30) * (0.190365) +
    (-0.20 <= eta && eta < -0.10) * (16.30 <= pt * cosh(eta) && pt * cosh(eta) < 17.10) * (0.212019) +
    (-0.20 <= eta && eta < -0.10) * (17.10 <= pt * cosh(eta) && pt * cosh(eta) < 17.90) * (0.233359) +
    (-0.20 <= eta && eta < -0.10) * (17.90 <= pt * cosh(eta) && pt * cosh(eta) < 18.80) * (0.254035) +
    (-0.20 <= eta && eta < -0.10) * (18.80 <= pt * cosh(eta) && pt * cosh(eta) < 19.80) * (0.274858) +
    (-0.20 <= eta && eta < -0.10) * (19.80 <= pt * cosh(eta) && pt * cosh(eta) < 21.00) * (0.296199) +
    (-0.20 <= eta && eta < -0.10) * (21.00 <= pt * cosh(eta) && pt * cosh(eta) < 22.40) * (0.318075) +
    (-0.20 <= eta && eta < -0.10) * (22.40 <= pt * cosh(eta) && pt * cosh(eta) < 24.00) * (0.339496) +
    (-0.20 <= eta && eta < -0.10) * (24.00 <= pt * cosh(eta) && pt * cosh(eta) < 26.00) * (0.360760) +
    (-0.20 <= eta && eta < -0.10) * (26.00 <= pt * cosh(eta) && pt * cosh(eta) < 28.50) * (0.382069) +
    (-0.20 <= eta && eta < -0.10) * (28.50 <= pt * cosh(eta) && pt * cosh(eta) < 31.90) * (0.403435) +
    (-0.20 <= eta && eta < -0.10) * (31.90 <= pt * cosh(eta) && pt * cosh(eta) < 36.90) * (0.425149) +
    (-0.20 <= eta && eta < -0.10) * (36.90 <= pt * cosh(eta) && pt * cosh(eta) < 45.70) * (0.447621) +
    (-0.20 <= eta && eta < -0.10) * (45.70 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.461296) +
    (-0.10 <= eta && eta < -0.00) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 10.20) * (0.000679) +
    (-0.10 <= eta && eta < -0.00) * (10.20 <= pt * cosh(eta) && pt * cosh(eta) < 11.20) * (0.019921) +
    (-0.10 <= eta && eta < -0.00) * (11.20 <= pt * cosh(eta) && pt * cosh(eta) < 12.00) * (0.039923) +
    (-0.10 <= eta && eta < -0.00) * (12.00 <= pt * cosh(eta) && pt * cosh(eta) < 12.70) * (0.061117) +
    (-0.10 <= eta && eta < -0.00) * (12.70 <= pt * cosh(eta) && pt * cosh(eta) < 13.40) * (0.083207) +
    (-0.10 <= eta && eta < -0.00) * (13.40 <= pt * cosh(eta) && pt * cosh(eta) < 14.00) * (0.104676) +
    (-0.10 <= eta && eta < -0.00) * (14.00 <= pt * cosh(eta) && pt * cosh(eta) < 14.70) * (0.126328) +
    (-0.10 <= eta && eta < -0.00) * (14.70 <= pt * cosh(eta) && pt * cosh(eta) < 15.40) * (0.149235) +
    (-0.10 <= eta && eta < -0.00) * (15.40 <= pt * cosh(eta) && pt * cosh(eta) < 16.10) * (0.171280) +
    (-0.10 <= eta && eta < -0.00) * (16.10 <= pt * cosh(eta) && pt * cosh(eta) < 16.80) * (0.192190) +
    (-0.10 <= eta && eta < -0.00) * (16.80 <= pt * cosh(eta) && pt * cosh(eta) < 17.60) * (0.213160) +
    (-0.10 <= eta && eta < -0.00) * (17.60 <= pt * cosh(eta) && pt * cosh(eta) < 18.50) * (0.235081) +
    (-0.10 <= eta && eta < -0.00) * (18.50 <= pt * cosh(eta) && pt * cosh(eta) < 19.50) * (0.257304) +
    (-0.10 <= eta && eta < -0.00) * (19.50 <= pt * cosh(eta) && pt * cosh(eta) < 20.60) * (0.279240) +
    (-0.10 <= eta && eta < -0.00) * (20.60 <= pt * cosh(eta) && pt * cosh(eta) < 21.80) * (0.300402) +
    (-0.10 <= eta && eta < -0.00) * (21.80 <= pt * cosh(eta) && pt * cosh(eta) < 23.20) * (0.321151) +
    (-0.10 <= eta && eta < -0.00) * (23.20 <= pt * cosh(eta) && pt * cosh(eta) < 24.90) * (0.342173) +
    (-0.10 <= eta && eta < -0.00) * (24.90 <= pt * cosh(eta) && pt * cosh(eta) < 27.00) * (0.363463) +
    (-0.10 <= eta && eta < -0.00) * (27.00 <= pt * cosh(eta) && pt * cosh(eta) < 29.70) * (0.384883) +
    (-0.10 <= eta && eta < -0.00) * (29.70 <= pt * cosh(eta) && pt * cosh(eta) < 33.40) * (0.406518) +
    (-0.10 <= eta && eta < -0.00) * (33.40 <= pt * cosh(eta) && pt * cosh(eta) < 38.90) * (0.428380) +
    (-0.10 <= eta && eta < -0.00) * (38.90 <= pt * cosh(eta) && pt * cosh(eta) < 49.00) * (0.451065) +
    (-0.10 <= eta && eta < -0.00) * (49.00 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.461887) +
    (-0.00 <= eta && eta < 0.10) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 10.00) * (0.000652) +
    (-0.00 <= eta && eta < 0.10) * (10.00 <= pt * cosh(eta) && pt * cosh(eta) < 11.10) * (0.019853) +
    (-0.00 <= eta && eta < 0.10) * (11.10 <= pt * cosh(eta) && pt * cosh(eta) < 11.90) * (0.041350) +
    (-0.00 <= eta && eta < 0.10) * (11.90 <= pt * cosh(eta) && pt * cosh(eta) < 12.60) * (0.063083) +
    (-0.00 <= eta && eta < 0.10) * (12.60 <= pt * cosh(eta) && pt * cosh(eta) < 13.30) * (0.085621) +
    (-0.00 <= eta && eta < 0.10) * (13.30 <= pt * cosh(eta) && pt * cosh(eta) < 13.90) * (0.107436) +
    (-0.00 <= eta && eta < 0.10) * (13.90 <= pt * cosh(eta) && pt * cosh(eta) < 14.50) * (0.127685) +
    (-0.00 <= eta && eta < 0.10) * (14.50 <= pt * cosh(eta) && pt * cosh(eta) < 15.20) * (0.149221) +
    (-0.00 <= eta && eta < 0.10) * (15.20 <= pt * cosh(eta) && pt * cosh(eta) < 15.90) * (0.171557) +
    (-0.00 <= eta && eta < 0.10) * (15.90 <= pt * cosh(eta) && pt * cosh(eta) < 16.60) * (0.192724) +
    (-0.00 <= eta && eta < 0.10) * (16.60 <= pt * cosh(eta) && pt * cosh(eta) < 17.40) * (0.213928) +
    (-0.00 <= eta && eta < 0.10) * (17.40 <= pt * cosh(eta) && pt * cosh(eta) < 18.30) * (0.236064) +
    (-0.00 <= eta && eta < 0.10) * (18.30 <= pt * cosh(eta) && pt * cosh(eta) < 19.30) * (0.258471) +
    (-0.00 <= eta && eta < 0.10) * (19.30 <= pt * cosh(eta) && pt * cosh(eta) < 20.40) * (0.280550) +
    (-0.00 <= eta && eta < 0.10) * (20.40 <= pt * cosh(eta) && pt * cosh(eta) < 21.60) * (0.301813) +
    (-0.00 <= eta && eta < 0.10) * (21.60 <= pt * cosh(eta) && pt * cosh(eta) < 23.00) * (0.322623) +
    (-0.00 <= eta && eta < 0.10) * (23.00 <= pt * cosh(eta) && pt * cosh(eta) < 24.70) * (0.343665) +
    (-0.00 <= eta && eta < 0.10) * (24.70 <= pt * cosh(eta) && pt * cosh(eta) < 26.80) * (0.364932) +
    (-0.00 <= eta && eta < 0.10) * (26.80 <= pt * cosh(eta) && pt * cosh(eta) < 29.50) * (0.386281) +
    (-0.00 <= eta && eta < 0.10) * (29.50 <= pt * cosh(eta) && pt * cosh(eta) < 33.20) * (0.407794) +
    (-0.00 <= eta && eta < 0.10) * (33.20 <= pt * cosh(eta) && pt * cosh(eta) < 38.80) * (0.429658) +
    (-0.00 <= eta && eta < 0.10) * (38.80 <= pt * cosh(eta) && pt * cosh(eta) < 49.20) * (0.452424) +
    (-0.00 <= eta && eta < 0.10) * (49.20 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.463042) +
    (0.10 <= eta && eta < 0.20) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 9.80) * (0.000688) +
    (0.10 <= eta && eta < 0.20) * (9.80 <= pt * cosh(eta) && pt * cosh(eta) < 10.80) * (0.020283) +
    (0.10 <= eta && eta < 0.20) * (10.80 <= pt * cosh(eta) && pt * cosh(eta) < 11.60) * (0.041390) +
    (0.10 <= eta && eta < 0.20) * (11.60 <= pt * cosh(eta) && pt * cosh(eta) < 12.30) * (0.063752) +
    (0.10 <= eta && eta < 0.20) * (12.30 <= pt * cosh(eta) && pt * cosh(eta) < 12.90) * (0.085249) +
    (0.10 <= eta && eta < 0.20) * (12.90 <= pt * cosh(eta) && pt * cosh(eta) < 13.50) * (0.105921) +
    (0.10 <= eta && eta < 0.20) * (13.50 <= pt * cosh(eta) && pt * cosh(eta) < 14.10) * (0.126724) +
    (0.10 <= eta && eta < 0.20) * (14.10 <= pt * cosh(eta) && pt * cosh(eta) < 14.70) * (0.147189) +
    (0.10 <= eta && eta < 0.20) * (14.70 <= pt * cosh(eta) && pt * cosh(eta) < 15.40) * (0.168592) +
    (0.10 <= eta && eta < 0.20) * (15.40 <= pt * cosh(eta) && pt * cosh(eta) < 16.10) * (0.190487) +
    (0.10 <= eta && eta < 0.20) * (16.10 <= pt * cosh(eta) && pt * cosh(eta) < 16.90) * (0.212403) +
    (0.10 <= eta && eta < 0.20) * (16.90 <= pt * cosh(eta) && pt * cosh(eta) < 17.70) * (0.233976) +
    (0.10 <= eta && eta < 0.20) * (17.70 <= pt * cosh(eta) && pt * cosh(eta) < 18.60) * (0.254850) +
    (0.10 <= eta && eta < 0.20) * (18.60 <= pt * cosh(eta) && pt * cosh(eta) < 19.60) * (0.275841) +
    (0.10 <= eta && eta < 0.20) * (19.60 <= pt * cosh(eta) && pt * cosh(eta) < 20.80) * (0.297319) +
    (0.10 <= eta && eta < 0.20) * (20.80 <= pt * cosh(eta) && pt * cosh(eta) < 22.20) * (0.319295) +
    (0.10 <= eta && eta < 0.20) * (22.20 <= pt * cosh(eta) && pt * cosh(eta) < 23.80) * (0.340773) +
    (0.10 <= eta && eta < 0.20) * (23.80 <= pt * cosh(eta) && pt * cosh(eta) < 25.80) * (0.362051) +
    (0.10 <= eta && eta < 0.20) * (25.80 <= pt * cosh(eta) && pt * cosh(eta) < 28.30) * (0.383326) +
    (0.10 <= eta && eta < 0.20) * (28.30 <= pt * cosh(eta) && pt * cosh(eta) < 31.70) * (0.404609) +
    (0.10 <= eta && eta < 0.20) * (31.70 <= pt * cosh(eta) && pt * cosh(eta) < 36.80) * (0.426384) +
    (0.10 <= eta && eta < 0.20) * (36.80 <= pt * cosh(eta) && pt * cosh(eta) < 45.90) * (0.449033) +
    (0.10 <= eta && eta < 0.20) * (45.90 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.462443) +
    (0.20 <= eta && eta < 0.30) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 9.70) * (0.000702) +
    (0.20 <= eta && eta < 0.30) * (9.70 <= pt * cosh(eta) && pt * cosh(eta) < 10.70) * (0.021415) +
    (0.20 <= eta && eta < 0.30) * (10.70 <= pt * cosh(eta) && pt * cosh(eta) < 11.40) * (0.041984) +
    (0.20 <= eta && eta < 0.30) * (11.40 <= pt * cosh(eta) && pt * cosh(eta) < 12.10) * (0.063242) +
    (0.20 <= eta && eta < 0.30) * (12.10 <= pt * cosh(eta) && pt * cosh(eta) < 12.70) * (0.085050) +
    (0.20 <= eta && eta < 0.30) * (12.70 <= pt * cosh(eta) && pt * cosh(eta) < 13.30) * (0.106045) +
    (0.20 <= eta && eta < 0.30) * (13.30 <= pt * cosh(eta) && pt * cosh(eta) < 13.90) * (0.127172) +
    (0.20 <= eta && eta < 0.30) * (13.90 <= pt * cosh(eta) && pt * cosh(eta) < 14.50) * (0.147941) +
    (0.20 <= eta && eta < 0.30) * (14.50 <= pt * cosh(eta) && pt * cosh(eta) < 15.20) * (0.169640) +
    (0.20 <= eta && eta < 0.30) * (15.20 <= pt * cosh(eta) && pt * cosh(eta) < 15.90) * (0.191806) +
    (0.20 <= eta && eta < 0.30) * (15.90 <= pt * cosh(eta) && pt * cosh(eta) < 16.70) * (0.213956) +
    (0.20 <= eta && eta < 0.30) * (16.70 <= pt * cosh(eta) && pt * cosh(eta) < 17.50) * (0.235721) +
    (0.20 <= eta && eta < 0.30) * (17.50 <= pt * cosh(eta) && pt * cosh(eta) < 18.40) * (0.256741) +
    (0.20 <= eta && eta < 0.30) * (18.40 <= pt * cosh(eta) && pt * cosh(eta) < 19.40) * (0.277837) +
    (0.20 <= eta && eta < 0.30) * (19.40 <= pt * cosh(eta) && pt * cosh(eta) < 20.60) * (0.299379) +
    (0.20 <= eta && eta < 0.30) * (20.60 <= pt * cosh(eta) && pt * cosh(eta) < 22.00) * (0.321372) +
    (0.20 <= eta && eta < 0.30) * (22.00 <= pt * cosh(eta) && pt * cosh(eta) < 23.60) * (0.342817) +
    (0.20 <= eta && eta < 0.30) * (23.60 <= pt * cosh(eta) && pt * cosh(eta) < 25.60) * (0.364014) +
    (0.20 <= eta && eta < 0.30) * (25.60 <= pt * cosh(eta) && pt * cosh(eta) < 28.20) * (0.385555) +
    (0.20 <= eta && eta < 0.30) * (28.20 <= pt * cosh(eta) && pt * cosh(eta) < 31.70) * (0.407162) +
    (0.20 <= eta && eta < 0.30) * (31.70 <= pt * cosh(eta) && pt * cosh(eta) < 37.00) * (0.429002) +
    (0.20 <= eta && eta < 0.30) * (37.00 <= pt * cosh(eta) && pt * cosh(eta) < 46.80) * (0.451798) +
    (0.20 <= eta && eta < 0.30) * (46.80 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.464291) +
    (0.30 <= eta && eta < 0.40) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 9.70) * (0.000726) +
    (0.30 <= eta && eta < 0.40) * (9.70 <= pt * cosh(eta) && pt * cosh(eta) < 10.70) * (0.020369) +
    (0.30 <= eta && eta < 0.40) * (10.70 <= pt * cosh(eta) && pt * cosh(eta) < 11.50) * (0.041656) +
    (0.30 <= eta && eta < 0.40) * (11.50 <= pt * cosh(eta) && pt * cosh(eta) < 12.20) * (0.064296) +
    (0.30 <= eta && eta < 0.40) * (12.20 <= pt * cosh(eta) && pt * cosh(eta) < 12.80) * (0.086047) +
    (0.30 <= eta && eta < 0.40) * (12.80 <= pt * cosh(eta) && pt * cosh(eta) < 13.40) * (0.106941) +
    (0.30 <= eta && eta < 0.40) * (13.40 <= pt * cosh(eta) && pt * cosh(eta) < 14.00) * (0.127942) +
    (0.30 <= eta && eta < 0.40) * (14.00 <= pt * cosh(eta) && pt * cosh(eta) < 14.60) * (0.148575) +
    (0.30 <= eta && eta < 0.40) * (14.60 <= pt * cosh(eta) && pt * cosh(eta) < 15.30) * (0.170125) +
    (0.30 <= eta && eta < 0.40) * (15.30 <= pt * cosh(eta) && pt * cosh(eta) < 16.00) * (0.192139) +
    (0.30 <= eta && eta < 0.40) * (16.00 <= pt * cosh(eta) && pt * cosh(eta) < 16.80) * (0.214143) +
    (0.30 <= eta && eta < 0.40) * (16.80 <= pt * cosh(eta) && pt * cosh(eta) < 17.60) * (0.235774) +
    (0.30 <= eta && eta < 0.40) * (17.60 <= pt * cosh(eta) && pt * cosh(eta) < 18.50) * (0.256674) +
    (0.30 <= eta && eta < 0.40) * (18.50 <= pt * cosh(eta) && pt * cosh(eta) < 19.50) * (0.277665) +
    (0.30 <= eta && eta < 0.40) * (19.50 <= pt * cosh(eta) && pt * cosh(eta) < 20.70) * (0.299113) +
    (0.30 <= eta && eta < 0.40) * (20.70 <= pt * cosh(eta) && pt * cosh(eta) < 22.10) * (0.321029) +
    (0.30 <= eta && eta < 0.40) * (22.10 <= pt * cosh(eta) && pt * cosh(eta) < 23.70) * (0.342418) +
    (0.30 <= eta && eta < 0.40) * (23.70 <= pt * cosh(eta) && pt * cosh(eta) < 25.70) * (0.363580) +
    (0.30 <= eta && eta < 0.40) * (25.70 <= pt * cosh(eta) && pt * cosh(eta) < 28.30) * (0.385108) +
    (0.30 <= eta && eta < 0.40) * (28.30 <= pt * cosh(eta) && pt * cosh(eta) < 31.80) * (0.406727) +
    (0.30 <= eta && eta < 0.40) * (31.80 <= pt * cosh(eta) && pt * cosh(eta) < 37.10) * (0.428608) +
    (0.30 <= eta && eta < 0.40) * (37.10 <= pt * cosh(eta) && pt * cosh(eta) < 46.80) * (0.451376) +
    (0.30 <= eta && eta < 0.40) * (46.80 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.463881) +
    (0.40 <= eta && eta < 0.50) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 9.60) * (0.000671) +
    (0.40 <= eta && eta < 0.50) * (9.60 <= pt * cosh(eta) && pt * cosh(eta) < 10.60) * (0.020318) +
    (0.40 <= eta && eta < 0.50) * (10.60 <= pt * cosh(eta) && pt * cosh(eta) < 11.30) * (0.040535) +
    (0.40 <= eta && eta < 0.50) * (11.30 <= pt * cosh(eta) && pt * cosh(eta) < 12.00) * (0.061663) +
    (0.40 <= eta && eta < 0.50) * (12.00 <= pt * cosh(eta) && pt * cosh(eta) < 12.60) * (0.083458) +
    (0.40 <= eta && eta < 0.50) * (12.60 <= pt * cosh(eta) && pt * cosh(eta) < 13.20) * (0.104511) +
    (0.40 <= eta && eta < 0.50) * (13.20 <= pt * cosh(eta) && pt * cosh(eta) < 13.80) * (0.125739) +
    (0.40 <= eta && eta < 0.50) * (13.80 <= pt * cosh(eta) && pt * cosh(eta) < 14.40) * (0.146635) +
    (0.40 <= eta && eta < 0.50) * (14.40 <= pt * cosh(eta) && pt * cosh(eta) < 15.10) * (0.168484) +
    (0.40 <= eta && eta < 0.50) * (15.10 <= pt * cosh(eta) && pt * cosh(eta) < 15.80) * (0.190813) +
    (0.40 <= eta && eta < 0.50) * (15.80 <= pt * cosh(eta) && pt * cosh(eta) < 16.60) * (0.213130) +
    (0.40 <= eta && eta < 0.50) * (16.60 <= pt * cosh(eta) && pt * cosh(eta) < 17.40) * (0.235057) +
    (0.40 <= eta && eta < 0.50) * (17.40 <= pt * cosh(eta) && pt * cosh(eta) < 18.30) * (0.256227) +
    (0.40 <= eta && eta < 0.50) * (18.30 <= pt * cosh(eta) && pt * cosh(eta) < 19.30) * (0.277468) +
    (0.40 <= eta && eta < 0.50) * (19.30 <= pt * cosh(eta) && pt * cosh(eta) < 20.50) * (0.299144) +
    (0.40 <= eta && eta < 0.50) * (20.50 <= pt * cosh(eta) && pt * cosh(eta) < 21.90) * (0.321260) +
    (0.40 <= eta && eta < 0.50) * (21.90 <= pt * cosh(eta) && pt * cosh(eta) < 23.50) * (0.342810) +
    (0.40 <= eta && eta < 0.50) * (23.50 <= pt * cosh(eta) && pt * cosh(eta) < 25.50) * (0.364092) +
    (0.40 <= eta && eta < 0.50) * (25.50 <= pt * cosh(eta) && pt * cosh(eta) < 28.10) * (0.385699) +
    (0.40 <= eta && eta < 0.50) * (28.10 <= pt * cosh(eta) && pt * cosh(eta) < 31.60) * (0.407349) +
    (0.40 <= eta && eta < 0.50) * (31.60 <= pt * cosh(eta) && pt * cosh(eta) < 36.90) * (0.429205) +
    (0.40 <= eta && eta < 0.50) * (36.90 <= pt * cosh(eta) && pt * cosh(eta) < 46.70) * (0.451986) +
    (0.40 <= eta && eta < 0.50) * (46.70 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.464527) +
    (0.50 <= eta && eta < 0.60) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 9.70) * (0.000629) +
    (0.50 <= eta && eta < 0.60) * (9.70 <= pt * cosh(eta) && pt * cosh(eta) < 10.70) * (0.019027) +
    (0.50 <= eta && eta < 0.60) * (10.70 <= pt * cosh(eta) && pt * cosh(eta) < 11.50) * (0.039568) +
    (0.50 <= eta && eta < 0.60) * (11.50 <= pt * cosh(eta) && pt * cosh(eta) < 12.20) * (0.061693) +
    (0.50 <= eta && eta < 0.60) * (12.20 <= pt * cosh(eta) && pt * cosh(eta) < 12.80) * (0.083113) +
    (0.50 <= eta && eta < 0.60) * (12.80 <= pt * cosh(eta) && pt * cosh(eta) < 13.40) * (0.103800) +
    (0.50 <= eta && eta < 0.60) * (13.40 <= pt * cosh(eta) && pt * cosh(eta) < 14.00) * (0.124675) +
    (0.50 <= eta && eta < 0.60) * (14.00 <= pt * cosh(eta) && pt * cosh(eta) < 14.60) * (0.145250) +
    (0.50 <= eta && eta < 0.60) * (14.60 <= pt * cosh(eta) && pt * cosh(eta) < 15.30) * (0.166798) +
    (0.50 <= eta && eta < 0.60) * (15.30 <= pt * cosh(eta) && pt * cosh(eta) < 16.00) * (0.188860) +
    (0.50 <= eta && eta < 0.60) * (16.00 <= pt * cosh(eta) && pt * cosh(eta) < 16.80) * (0.210955) +
    (0.50 <= eta && eta < 0.60) * (16.80 <= pt * cosh(eta) && pt * cosh(eta) < 17.60) * (0.232711) +
    (0.50 <= eta && eta < 0.60) * (17.60 <= pt * cosh(eta) && pt * cosh(eta) < 18.50) * (0.253762) +
    (0.50 <= eta && eta < 0.60) * (18.50 <= pt * cosh(eta) && pt * cosh(eta) < 19.50) * (0.274929) +
    (0.50 <= eta && eta < 0.60) * (19.50 <= pt * cosh(eta) && pt * cosh(eta) < 20.70) * (0.296580) +
    (0.50 <= eta && eta < 0.60) * (20.70 <= pt * cosh(eta) && pt * cosh(eta) < 22.00) * (0.317952) +
    (0.50 <= eta && eta < 0.60) * (22.00 <= pt * cosh(eta) && pt * cosh(eta) < 23.60) * (0.339024) +
    (0.50 <= eta && eta < 0.60) * (23.60 <= pt * cosh(eta) && pt * cosh(eta) < 25.50) * (0.360155) +
    (0.50 <= eta && eta < 0.60) * (25.50 <= pt * cosh(eta) && pt * cosh(eta) < 28.00) * (0.381450) +
    (0.50 <= eta && eta < 0.60) * (28.00 <= pt * cosh(eta) && pt * cosh(eta) < 31.30) * (0.402963) +
    (0.50 <= eta && eta < 0.60) * (31.30 <= pt * cosh(eta) && pt * cosh(eta) < 36.20) * (0.424674) +
    (0.50 <= eta && eta < 0.60) * (36.20 <= pt * cosh(eta) && pt * cosh(eta) < 44.70) * (0.447121) +
    (0.50 <= eta && eta < 0.60) * (44.70 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.461667) +
    (0.60 <= eta && eta < 0.70) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 9.80) * (0.000675) +
    (0.60 <= eta && eta < 0.70) * (9.80 <= pt * cosh(eta) && pt * cosh(eta) < 10.80) * (0.019837) +
    (0.60 <= eta && eta < 0.70) * (10.80 <= pt * cosh(eta) && pt * cosh(eta) < 11.60) * (0.040690) +
    (0.60 <= eta && eta < 0.70) * (11.60 <= pt * cosh(eta) && pt * cosh(eta) < 12.30) * (0.062882) +
    (0.60 <= eta && eta < 0.70) * (12.30 <= pt * cosh(eta) && pt * cosh(eta) < 12.90) * (0.084270) +
    (0.60 <= eta && eta < 0.70) * (12.90 <= pt * cosh(eta) && pt * cosh(eta) < 13.50) * (0.104874) +
    (0.60 <= eta && eta < 0.70) * (13.50 <= pt * cosh(eta) && pt * cosh(eta) < 14.10) * (0.125635) +
    (0.60 <= eta && eta < 0.70) * (14.10 <= pt * cosh(eta) && pt * cosh(eta) < 14.70) * (0.146081) +
    (0.60 <= eta && eta < 0.70) * (14.70 <= pt * cosh(eta) && pt * cosh(eta) < 15.40) * (0.167483) +
    (0.60 <= eta && eta < 0.70) * (15.40 <= pt * cosh(eta) && pt * cosh(eta) < 16.10) * (0.189393) +
    (0.60 <= eta && eta < 0.70) * (16.10 <= pt * cosh(eta) && pt * cosh(eta) < 16.90) * (0.211339) +
    (0.60 <= eta && eta < 0.70) * (16.90 <= pt * cosh(eta) && pt * cosh(eta) < 17.70) * (0.232954) +
    (0.60 <= eta && eta < 0.70) * (17.70 <= pt * cosh(eta) && pt * cosh(eta) < 18.60) * (0.253877) +
    (0.60 <= eta && eta < 0.70) * (18.60 <= pt * cosh(eta) && pt * cosh(eta) < 19.60) * (0.274927) +
    (0.60 <= eta && eta < 0.70) * (19.60 <= pt * cosh(eta) && pt * cosh(eta) < 20.80) * (0.296472) +
    (0.60 <= eta && eta < 0.70) * (20.80 <= pt * cosh(eta) && pt * cosh(eta) < 22.20) * (0.318523) +
    (0.60 <= eta && eta < 0.70) * (22.20 <= pt * cosh(eta) && pt * cosh(eta) < 23.80) * (0.340081) +
    (0.60 <= eta && eta < 0.70) * (23.80 <= pt * cosh(eta) && pt * cosh(eta) < 25.80) * (0.361442) +
    (0.60 <= eta && eta < 0.70) * (25.80 <= pt * cosh(eta) && pt * cosh(eta) < 28.30) * (0.382805) +
    (0.60 <= eta && eta < 0.70) * (28.30 <= pt * cosh(eta) && pt * cosh(eta) < 31.70) * (0.404178) +
    (0.60 <= eta && eta < 0.70) * (31.70 <= pt * cosh(eta) && pt * cosh(eta) < 36.80) * (0.426049) +
    (0.60 <= eta && eta < 0.70) * (36.80 <= pt * cosh(eta) && pt * cosh(eta) < 45.80) * (0.448688) +
    (0.60 <= eta && eta < 0.70) * (45.80 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.462189) +
    (0.70 <= eta && eta < 0.80) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 9.80) * (0.000676) +
    (0.70 <= eta && eta < 0.80) * (9.80 <= pt * cosh(eta) && pt * cosh(eta) < 10.80) * (0.019969) +
    (0.70 <= eta && eta < 0.80) * (10.80 <= pt * cosh(eta) && pt * cosh(eta) < 11.60) * (0.040896) +
    (0.70 <= eta && eta < 0.80) * (11.60 <= pt * cosh(eta) && pt * cosh(eta) < 12.30) * (0.063138) +
    (0.70 <= eta && eta < 0.80) * (12.30 <= pt * cosh(eta) && pt * cosh(eta) < 12.90) * (0.084558) +
    (0.70 <= eta && eta < 0.80) * (12.90 <= pt * cosh(eta) && pt * cosh(eta) < 13.50) * (0.105182) +
    (0.70 <= eta && eta < 0.80) * (13.50 <= pt * cosh(eta) && pt * cosh(eta) < 14.10) * (0.125956) +
    (0.70 <= eta && eta < 0.80) * (14.10 <= pt * cosh(eta) && pt * cosh(eta) < 14.70) * (0.146407) +
    (0.70 <= eta && eta < 0.80) * (14.70 <= pt * cosh(eta) && pt * cosh(eta) < 15.40) * (0.167810) +
    (0.70 <= eta && eta < 0.80) * (15.40 <= pt * cosh(eta) && pt * cosh(eta) < 16.10) * (0.189716) +
    (0.70 <= eta && eta < 0.80) * (16.10 <= pt * cosh(eta) && pt * cosh(eta) < 16.90) * (0.211653) +
    (0.70 <= eta && eta < 0.80) * (16.90 <= pt * cosh(eta) && pt * cosh(eta) < 17.70) * (0.233255) +
    (0.70 <= eta && eta < 0.80) * (17.70 <= pt * cosh(eta) && pt * cosh(eta) < 18.60) * (0.254164) +
    (0.70 <= eta && eta < 0.80) * (18.60 <= pt * cosh(eta) && pt * cosh(eta) < 19.60) * (0.275196) +
    (0.70 <= eta && eta < 0.80) * (19.60 <= pt * cosh(eta) && pt * cosh(eta) < 20.80) * (0.296722) +
    (0.70 <= eta && eta < 0.80) * (20.80 <= pt * cosh(eta) && pt * cosh(eta) < 22.20) * (0.318751) +
    (0.70 <= eta && eta < 0.80) * (22.20 <= pt * cosh(eta) && pt * cosh(eta) < 23.80) * (0.340285) +
    (0.70 <= eta && eta < 0.80) * (23.80 <= pt * cosh(eta) && pt * cosh(eta) < 25.80) * (0.361622) +
    (0.70 <= eta && eta < 0.80) * (25.80 <= pt * cosh(eta) && pt * cosh(eta) < 28.30) * (0.382958) +
    (0.70 <= eta && eta < 0.80) * (28.30 <= pt * cosh(eta) && pt * cosh(eta) < 31.70) * (0.404306) +
    (0.70 <= eta && eta < 0.80) * (31.70 <= pt * cosh(eta) && pt * cosh(eta) < 36.80) * (0.426148) +
    (0.70 <= eta && eta < 0.80) * (36.80 <= pt * cosh(eta) && pt * cosh(eta) < 45.80) * (0.448757) +
    (0.70 <= eta && eta < 0.80) * (45.80 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.462240) +
    (0.80 <= eta && eta < 0.90) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 10.00) * (0.000666) +
    (0.80 <= eta && eta < 0.90) * (10.00 <= pt * cosh(eta) && pt * cosh(eta) < 11.10) * (0.019973) +
    (0.80 <= eta && eta < 0.90) * (11.10 <= pt * cosh(eta) && pt * cosh(eta) < 11.90) * (0.041538) +
    (0.80 <= eta && eta < 0.90) * (11.90 <= pt * cosh(eta) && pt * cosh(eta) < 12.60) * (0.063316) +
    (0.80 <= eta && eta < 0.90) * (12.60 <= pt * cosh(eta) && pt * cosh(eta) < 13.30) * (0.085883) +
    (0.80 <= eta && eta < 0.90) * (13.30 <= pt * cosh(eta) && pt * cosh(eta) < 13.90) * (0.107717) +
    (0.80 <= eta && eta < 0.90) * (13.90 <= pt * cosh(eta) && pt * cosh(eta) < 14.50) * (0.127976) +
    (0.80 <= eta && eta < 0.90) * (14.50 <= pt * cosh(eta) && pt * cosh(eta) < 15.20) * (0.149517) +
    (0.80 <= eta && eta < 0.90) * (15.20 <= pt * cosh(eta) && pt * cosh(eta) < 15.90) * (0.171852) +
    (0.80 <= eta && eta < 0.90) * (15.90 <= pt * cosh(eta) && pt * cosh(eta) < 16.60) * (0.193015) +
    (0.80 <= eta && eta < 0.90) * (16.60 <= pt * cosh(eta) && pt * cosh(eta) < 17.40) * (0.214211) +
    (0.80 <= eta && eta < 0.90) * (17.40 <= pt * cosh(eta) && pt * cosh(eta) < 18.30) * (0.236336) +
    (0.80 <= eta && eta < 0.90) * (18.30 <= pt * cosh(eta) && pt * cosh(eta) < 19.30) * (0.258728) +
    (0.80 <= eta && eta < 0.90) * (19.30 <= pt * cosh(eta) && pt * cosh(eta) < 20.40) * (0.280790) +
    (0.80 <= eta && eta < 0.90) * (20.40 <= pt * cosh(eta) && pt * cosh(eta) < 21.60) * (0.302035) +
    (0.80 <= eta && eta < 0.90) * (21.60 <= pt * cosh(eta) && pt * cosh(eta) < 23.00) * (0.322825) +
    (0.80 <= eta && eta < 0.90) * (23.00 <= pt * cosh(eta) && pt * cosh(eta) < 24.70) * (0.343847) +
    (0.80 <= eta && eta < 0.90) * (24.70 <= pt * cosh(eta) && pt * cosh(eta) < 26.80) * (0.365091) +
    (0.80 <= eta && eta < 0.90) * (26.80 <= pt * cosh(eta) && pt * cosh(eta) < 29.50) * (0.386417) +
    (0.80 <= eta && eta < 0.90) * (29.50 <= pt * cosh(eta) && pt * cosh(eta) < 33.20) * (0.407905) +
    (0.80 <= eta && eta < 0.90) * (33.20 <= pt * cosh(eta) && pt * cosh(eta) < 38.80) * (0.429743) +
    (0.80 <= eta && eta < 0.90) * (38.80 <= pt * cosh(eta) && pt * cosh(eta) < 49.20) * (0.452482) +
    (0.80 <= eta && eta < 0.90) * (49.20 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.463087) +
    (0.90 <= eta && eta < 1.00) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 10.00) * (0.000631) +
    (0.90 <= eta && eta < 1.00) * (10.00 <= pt * cosh(eta) && pt * cosh(eta) < 11.10) * (0.020111) +
    (0.90 <= eta && eta < 1.00) * (11.10 <= pt * cosh(eta) && pt * cosh(eta) < 11.90) * (0.041754) +
    (0.90 <= eta && eta < 1.00) * (11.90 <= pt * cosh(eta) && pt * cosh(eta) < 12.60) * (0.063581) +
    (0.90 <= eta && eta < 1.00) * (12.60 <= pt * cosh(eta) && pt * cosh(eta) < 13.20) * (0.084516) +
    (0.90 <= eta && eta < 1.00) * (13.20 <= pt * cosh(eta) && pt * cosh(eta) < 13.80) * (0.104653) +
    (0.90 <= eta && eta < 1.00) * (13.80 <= pt * cosh(eta) && pt * cosh(eta) < 14.40) * (0.124945) +
    (0.90 <= eta && eta < 1.00) * (14.40 <= pt * cosh(eta) && pt * cosh(eta) < 15.10) * (0.146583) +
    (0.90 <= eta && eta < 1.00) * (15.10 <= pt * cosh(eta) && pt * cosh(eta) < 15.80) * (0.169067) +
    (0.90 <= eta && eta < 1.00) * (15.80 <= pt * cosh(eta) && pt * cosh(eta) < 16.50) * (0.190404) +
    (0.90 <= eta && eta < 1.00) * (16.50 <= pt * cosh(eta) && pt * cosh(eta) < 17.30) * (0.211798) +
    (0.90 <= eta && eta < 1.00) * (17.30 <= pt * cosh(eta) && pt * cosh(eta) < 18.20) * (0.234147) +
    (0.90 <= eta && eta < 1.00) * (18.20 <= pt * cosh(eta) && pt * cosh(eta) < 19.10) * (0.255662) +
    (0.90 <= eta && eta < 1.00) * (19.10 <= pt * cosh(eta) && pt * cosh(eta) < 20.20) * (0.277068) +
    (0.90 <= eta && eta < 1.00) * (20.20 <= pt * cosh(eta) && pt * cosh(eta) < 21.40) * (0.298797) +
    (0.90 <= eta && eta < 1.00) * (21.40 <= pt * cosh(eta) && pt * cosh(eta) < 22.80) * (0.320057) +
    (0.90 <= eta && eta < 1.00) * (22.80 <= pt * cosh(eta) && pt * cosh(eta) < 24.50) * (0.341541) +
    (0.90 <= eta && eta < 1.00) * (24.50 <= pt * cosh(eta) && pt * cosh(eta) < 26.60) * (0.363233) +
    (0.90 <= eta && eta < 1.00) * (26.60 <= pt * cosh(eta) && pt * cosh(eta) < 29.20) * (0.384597) +
    (0.90 <= eta && eta < 1.00) * (29.20 <= pt * cosh(eta) && pt * cosh(eta) < 32.80) * (0.405989) +
    (0.90 <= eta && eta < 1.00) * (32.80 <= pt * cosh(eta) && pt * cosh(eta) < 38.20) * (0.427888) +
    (0.90 <= eta && eta < 1.00) * (38.20 <= pt * cosh(eta) && pt * cosh(eta) < 48.00) * (0.450605) +
    (0.90 <= eta && eta < 1.00) * (48.00 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.462220) +
    (1.00 <= eta && eta < 1.10) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.000000) }

  add EfficiencyFormula {2212} {2212} {
    (eta<-1.00 || eta>=1.00 || pt * cosh(eta) < 0.10 || pt * cosh(eta) >= 50.00) * ( 0.00 ) +
    (-1.00 <= eta && eta < -0.90) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.60) * (0.000000) +
    (-1.00 <= eta && eta < -0.90) * (0.60 <= pt * cosh(eta) && pt * cosh(eta) < 9.50) * (0.999393) +
    (-1.00 <= eta && eta < -0.90) * (9.50 <= pt * cosh(eta) && pt * cosh(eta) < 10.30) * (0.979966) +
    (-1.00 <= eta && eta < -0.90) * (10.30 <= pt * cosh(eta) && pt * cosh(eta) < 10.90) * (0.959319) +
    (-1.00 <= eta && eta < -0.90) * (10.90 <= pt * cosh(eta) && pt * cosh(eta) < 11.40) * (0.938208) +
    (-1.00 <= eta && eta < -0.90) * (11.40 <= pt * cosh(eta) && pt * cosh(eta) < 11.90) * (0.916518) +
    (-1.00 <= eta && eta < -0.90) * (11.90 <= pt * cosh(eta) && pt * cosh(eta) < 12.40) * (0.893591) +
    (-1.00 <= eta && eta < -0.90) * (12.40 <= pt * cosh(eta) && pt * cosh(eta) < 12.90) * (0.870367) +
    (-1.00 <= eta && eta < -0.90) * (12.90 <= pt * cosh(eta) && pt * cosh(eta) < 13.40) * (0.847581) +
    (-1.00 <= eta && eta < -0.90) * (13.40 <= pt * cosh(eta) && pt * cosh(eta) < 13.90) * (0.825759) +
    (-1.00 <= eta && eta < -0.90) * (13.90 <= pt * cosh(eta) && pt * cosh(eta) < 14.50) * (0.803286) +
    (-1.00 <= eta && eta < -0.90) * (14.50 <= pt * cosh(eta) && pt * cosh(eta) < 15.10) * (0.780814) +
    (-1.00 <= eta && eta < -0.90) * (15.10 <= pt * cosh(eta) && pt * cosh(eta) < 15.80) * (0.759031) +
    (-1.00 <= eta && eta < -0.90) * (15.80 <= pt * cosh(eta) && pt * cosh(eta) < 16.60) * (0.737077) +
    (-1.00 <= eta && eta < -0.90) * (16.60 <= pt * cosh(eta) && pt * cosh(eta) < 17.60) * (0.714871) +
    (-1.00 <= eta && eta < -0.90) * (17.60 <= pt * cosh(eta) && pt * cosh(eta) < 18.80) * (0.692935) +
    (-1.00 <= eta && eta < -0.90) * (18.80 <= pt * cosh(eta) && pt * cosh(eta) < 20.50) * (0.671054) +
    (-1.00 <= eta && eta < -0.90) * (20.50 <= pt * cosh(eta) && pt * cosh(eta) < 23.20) * (0.648575) +
    (-1.00 <= eta && eta < -0.90) * (23.20 <= pt * cosh(eta) && pt * cosh(eta) < 32.00) * (0.621769) +
    (-1.00 <= eta && eta < -0.90) * (32.00 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.605753) +
    (-0.90 <= eta && eta < -0.80) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.60) * (0.000000) +
    (-0.90 <= eta && eta < -0.80) * (0.60 <= pt * cosh(eta) && pt * cosh(eta) < 9.40) * (0.999362) +
    (-0.90 <= eta && eta < -0.80) * (9.40 <= pt * cosh(eta) && pt * cosh(eta) < 10.20) * (0.978828) +
    (-0.90 <= eta && eta < -0.80) * (10.20 <= pt * cosh(eta) && pt * cosh(eta) < 10.80) * (0.957392) +
    (-0.90 <= eta && eta < -0.80) * (10.80 <= pt * cosh(eta) && pt * cosh(eta) < 11.30) * (0.935603) +
    (-0.90 <= eta && eta < -0.80) * (11.30 <= pt * cosh(eta) && pt * cosh(eta) < 11.80) * (0.913368) +
    (-0.90 <= eta && eta < -0.80) * (11.80 <= pt * cosh(eta) && pt * cosh(eta) < 12.30) * (0.890006) +
    (-0.90 <= eta && eta < -0.80) * (12.30 <= pt * cosh(eta) && pt * cosh(eta) < 12.80) * (0.866467) +
    (-0.90 <= eta && eta < -0.80) * (12.80 <= pt * cosh(eta) && pt * cosh(eta) < 13.30) * (0.843486) +
    (-0.90 <= eta && eta < -0.80) * (13.30 <= pt * cosh(eta) && pt * cosh(eta) < 13.80) * (0.821573) +
    (-0.90 <= eta && eta < -0.80) * (13.80 <= pt * cosh(eta) && pt * cosh(eta) < 14.40) * (0.799102) +
    (-0.90 <= eta && eta < -0.80) * (14.40 <= pt * cosh(eta) && pt * cosh(eta) < 15.00) * (0.776721) +
    (-0.90 <= eta && eta < -0.80) * (15.00 <= pt * cosh(eta) && pt * cosh(eta) < 15.70) * (0.755112) +
    (-0.90 <= eta && eta < -0.80) * (15.70 <= pt * cosh(eta) && pt * cosh(eta) < 16.50) * (0.733415) +
    (-0.90 <= eta && eta < -0.80) * (16.50 <= pt * cosh(eta) && pt * cosh(eta) < 17.50) * (0.711555) +
    (-0.90 <= eta && eta < -0.80) * (17.50 <= pt * cosh(eta) && pt * cosh(eta) < 18.80) * (0.689217) +
    (-0.90 <= eta && eta < -0.80) * (18.80 <= pt * cosh(eta) && pt * cosh(eta) < 20.60) * (0.666846) +
    (-0.90 <= eta && eta < -0.80) * (20.60 <= pt * cosh(eta) && pt * cosh(eta) < 23.70) * (0.643904) +
    (-0.90 <= eta && eta < -0.80) * (23.70 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.609976) +
    (-0.80 <= eta && eta < -0.70) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.60) * (0.000000) +
    (-0.80 <= eta && eta < -0.70) * (0.60 <= pt * cosh(eta) && pt * cosh(eta) < 9.20) * (0.999417) +
    (-0.80 <= eta && eta < -0.70) * (9.20 <= pt * cosh(eta) && pt * cosh(eta) < 10.00) * (0.980438) +
    (-0.80 <= eta && eta < -0.70) * (10.00 <= pt * cosh(eta) && pt * cosh(eta) < 10.60) * (0.959304) +
    (-0.80 <= eta && eta < -0.70) * (10.60 <= pt * cosh(eta) && pt * cosh(eta) < 11.10) * (0.937529) +
    (-0.80 <= eta && eta < -0.70) * (11.10 <= pt * cosh(eta) && pt * cosh(eta) < 11.60) * (0.915124) +
    (-0.80 <= eta && eta < -0.70) * (11.60 <= pt * cosh(eta) && pt * cosh(eta) < 12.10) * (0.891475) +
    (-0.80 <= eta && eta < -0.70) * (12.10 <= pt * cosh(eta) && pt * cosh(eta) < 12.60) * (0.867589) +
    (-0.80 <= eta && eta < -0.70) * (12.60 <= pt * cosh(eta) && pt * cosh(eta) < 13.10) * (0.844246) +
    (-0.80 <= eta && eta < -0.70) * (13.10 <= pt * cosh(eta) && pt * cosh(eta) < 13.60) * (0.821987) +
    (-0.80 <= eta && eta < -0.70) * (13.60 <= pt * cosh(eta) && pt * cosh(eta) < 14.20) * (0.799175) +
    (-0.80 <= eta && eta < -0.70) * (14.20 <= pt * cosh(eta) && pt * cosh(eta) < 14.80) * (0.776484) +
    (-0.80 <= eta && eta < -0.70) * (14.80 <= pt * cosh(eta) && pt * cosh(eta) < 15.50) * (0.754614) +
    (-0.80 <= eta && eta < -0.70) * (15.50 <= pt * cosh(eta) && pt * cosh(eta) < 16.30) * (0.732703) +
    (-0.80 <= eta && eta < -0.70) * (16.30 <= pt * cosh(eta) && pt * cosh(eta) < 17.30) * (0.710687) +
    (-0.80 <= eta && eta < -0.70) * (17.30 <= pt * cosh(eta) && pt * cosh(eta) < 18.60) * (0.688267) +
    (-0.80 <= eta && eta < -0.70) * (18.60 <= pt * cosh(eta) && pt * cosh(eta) < 20.40) * (0.665905) +
    (-0.80 <= eta && eta < -0.70) * (20.40 <= pt * cosh(eta) && pt * cosh(eta) < 23.50) * (0.643095) +
    (-0.80 <= eta && eta < -0.70) * (23.50 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.609665) +
    (-0.70 <= eta && eta < -0.60) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.60) * (0.000000) +
    (-0.70 <= eta && eta < -0.60) * (0.60 <= pt * cosh(eta) && pt * cosh(eta) < 9.20) * (0.999371) +
    (-0.70 <= eta && eta < -0.60) * (9.20 <= pt * cosh(eta) && pt * cosh(eta) < 10.00) * (0.979255) +
    (-0.70 <= eta && eta < -0.60) * (10.00 <= pt * cosh(eta) && pt * cosh(eta) < 10.60) * (0.957127) +
    (-0.70 <= eta && eta < -0.60) * (10.60 <= pt * cosh(eta) && pt * cosh(eta) < 11.10) * (0.934823) +
    (-0.70 <= eta && eta < -0.60) * (11.10 <= pt * cosh(eta) && pt * cosh(eta) < 11.60) * (0.912064) +
    (-0.70 <= eta && eta < -0.60) * (11.60 <= pt * cosh(eta) && pt * cosh(eta) < 12.10) * (0.888193) +
    (-0.70 <= eta && eta < -0.60) * (12.10 <= pt * cosh(eta) && pt * cosh(eta) < 12.60) * (0.864204) +
    (-0.70 <= eta && eta < -0.60) * (12.60 <= pt * cosh(eta) && pt * cosh(eta) < 13.10) * (0.840855) +
    (-0.70 <= eta && eta < -0.60) * (13.10 <= pt * cosh(eta) && pt * cosh(eta) < 13.60) * (0.818667) +
    (-0.70 <= eta && eta < -0.60) * (13.60 <= pt * cosh(eta) && pt * cosh(eta) < 14.20) * (0.795994) +
    (-0.70 <= eta && eta < -0.60) * (14.20 <= pt * cosh(eta) && pt * cosh(eta) < 14.80) * (0.773499) +
    (-0.70 <= eta && eta < -0.60) * (14.80 <= pt * cosh(eta) && pt * cosh(eta) < 15.50) * (0.751868) +
    (-0.70 <= eta && eta < -0.60) * (15.50 <= pt * cosh(eta) && pt * cosh(eta) < 16.30) * (0.730241) +
    (-0.70 <= eta && eta < -0.60) * (16.30 <= pt * cosh(eta) && pt * cosh(eta) < 17.30) * (0.708551) +
    (-0.70 <= eta && eta < -0.60) * (17.30 <= pt * cosh(eta) && pt * cosh(eta) < 18.60) * (0.686502) +
    (-0.70 <= eta && eta < -0.60) * (18.60 <= pt * cosh(eta) && pt * cosh(eta) < 20.40) * (0.664546) +
    (-0.70 <= eta && eta < -0.60) * (20.40 <= pt * cosh(eta) && pt * cosh(eta) < 23.60) * (0.641866) +
    (-0.70 <= eta && eta < -0.60) * (23.60 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.609387) +
    (-0.60 <= eta && eta < -0.50) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.60) * (0.000000) +
    (-0.60 <= eta && eta < -0.50) * (0.60 <= pt * cosh(eta) && pt * cosh(eta) < 9.20) * (0.999416) +
    (-0.60 <= eta && eta < -0.50) * (9.20 <= pt * cosh(eta) && pt * cosh(eta) < 10.00) * (0.979498) +
    (-0.60 <= eta && eta < -0.50) * (10.00 <= pt * cosh(eta) && pt * cosh(eta) < 10.60) * (0.957852) +
    (-0.60 <= eta && eta < -0.50) * (10.60 <= pt * cosh(eta) && pt * cosh(eta) < 11.10) * (0.935720) +
    (-0.60 <= eta && eta < -0.50) * (11.10 <= pt * cosh(eta) && pt * cosh(eta) < 11.60) * (0.913076) +
    (-0.60 <= eta && eta < -0.50) * (11.60 <= pt * cosh(eta) && pt * cosh(eta) < 12.10) * (0.889276) +
    (-0.60 <= eta && eta < -0.50) * (12.10 <= pt * cosh(eta) && pt * cosh(eta) < 12.60) * (0.865320) +
    (-0.60 <= eta && eta < -0.50) * (12.60 <= pt * cosh(eta) && pt * cosh(eta) < 13.10) * (0.841972) +
    (-0.60 <= eta && eta < -0.50) * (13.10 <= pt * cosh(eta) && pt * cosh(eta) < 13.60) * (0.819759) +
    (-0.60 <= eta && eta < -0.50) * (13.60 <= pt * cosh(eta) && pt * cosh(eta) < 14.20) * (0.797040) +
    (-0.60 <= eta && eta < -0.50) * (14.20 <= pt * cosh(eta) && pt * cosh(eta) < 14.80) * (0.774479) +
    (-0.60 <= eta && eta < -0.50) * (14.80 <= pt * cosh(eta) && pt * cosh(eta) < 15.50) * (0.752769) +
    (-0.60 <= eta && eta < -0.50) * (15.50 <= pt * cosh(eta) && pt * cosh(eta) < 16.30) * (0.731048) +
    (-0.60 <= eta && eta < -0.50) * (16.30 <= pt * cosh(eta) && pt * cosh(eta) < 17.30) * (0.709251) +
    (-0.60 <= eta && eta < -0.50) * (17.30 <= pt * cosh(eta) && pt * cosh(eta) < 18.60) * (0.687080) +
    (-0.60 <= eta && eta < -0.50) * (18.60 <= pt * cosh(eta) && pt * cosh(eta) < 20.40) * (0.664991) +
    (-0.60 <= eta && eta < -0.50) * (20.40 <= pt * cosh(eta) && pt * cosh(eta) < 23.60) * (0.642162) +
    (-0.60 <= eta && eta < -0.50) * (23.60 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.609449) +
    (-0.50 <= eta && eta < -0.40) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.60) * (0.000000) +
    (-0.50 <= eta && eta < -0.40) * (0.60 <= pt * cosh(eta) && pt * cosh(eta) < 9.10) * (0.999377) +
    (-0.50 <= eta && eta < -0.40) * (9.10 <= pt * cosh(eta) && pt * cosh(eta) < 9.90) * (0.978482) +
    (-0.50 <= eta && eta < -0.40) * (9.90 <= pt * cosh(eta) && pt * cosh(eta) < 10.50) * (0.955953) +
    (-0.50 <= eta && eta < -0.40) * (10.50 <= pt * cosh(eta) && pt * cosh(eta) < 11.00) * (0.933137) +
    (-0.50 <= eta && eta < -0.40) * (11.00 <= pt * cosh(eta) && pt * cosh(eta) < 11.50) * (0.909949) +
    (-0.50 <= eta && eta < -0.40) * (11.50 <= pt * cosh(eta) && pt * cosh(eta) < 12.00) * (0.885719) +
    (-0.50 <= eta && eta < -0.40) * (12.00 <= pt * cosh(eta) && pt * cosh(eta) < 12.50) * (0.861459) +
    (-0.50 <= eta && eta < -0.40) * (12.50 <= pt * cosh(eta) && pt * cosh(eta) < 13.00) * (0.837929) +
    (-0.50 <= eta && eta < -0.40) * (13.00 <= pt * cosh(eta) && pt * cosh(eta) < 13.50) * (0.815641) +
    (-0.50 <= eta && eta < -0.40) * (13.50 <= pt * cosh(eta) && pt * cosh(eta) < 14.10) * (0.792940) +
    (-0.50 <= eta && eta < -0.40) * (14.10 <= pt * cosh(eta) && pt * cosh(eta) < 14.70) * (0.770487) +
    (-0.50 <= eta && eta < -0.40) * (14.70 <= pt * cosh(eta) && pt * cosh(eta) < 15.40) * (0.748964) +
    (-0.50 <= eta && eta < -0.40) * (15.40 <= pt * cosh(eta) && pt * cosh(eta) < 16.20) * (0.727511) +
    (-0.50 <= eta && eta < -0.40) * (16.20 <= pt * cosh(eta) && pt * cosh(eta) < 17.20) * (0.706066) +
    (-0.50 <= eta && eta < -0.40) * (17.20 <= pt * cosh(eta) && pt * cosh(eta) < 18.50) * (0.684338) +
    (-0.50 <= eta && eta < -0.40) * (18.50 <= pt * cosh(eta) && pt * cosh(eta) < 20.40) * (0.662237) +
    (-0.50 <= eta && eta < -0.40) * (20.40 <= pt * cosh(eta) && pt * cosh(eta) < 23.90) * (0.639038) +
    (-0.50 <= eta && eta < -0.40) * (23.90 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.608741) +
    (-0.40 <= eta && eta < -0.30) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.60) * (0.000000) +
    (-0.40 <= eta && eta < -0.30) * (0.60 <= pt * cosh(eta) && pt * cosh(eta) < 9.10) * (0.999391) +
    (-0.40 <= eta && eta < -0.30) * (9.10 <= pt * cosh(eta) && pt * cosh(eta) < 9.90) * (0.979179) +
    (-0.40 <= eta && eta < -0.30) * (9.90 <= pt * cosh(eta) && pt * cosh(eta) < 10.50) * (0.957431) +
    (-0.40 <= eta && eta < -0.30) * (10.50 <= pt * cosh(eta) && pt * cosh(eta) < 11.00) * (0.934965) +
    (-0.40 <= eta && eta < -0.30) * (11.00 <= pt * cosh(eta) && pt * cosh(eta) < 11.50) * (0.912006) +
    (-0.40 <= eta && eta < -0.30) * (11.50 <= pt * cosh(eta) && pt * cosh(eta) < 12.00) * (0.887916) +
    (-0.40 <= eta && eta < -0.30) * (12.00 <= pt * cosh(eta) && pt * cosh(eta) < 12.50) * (0.863717) +
    (-0.40 <= eta && eta < -0.30) * (12.50 <= pt * cosh(eta) && pt * cosh(eta) < 13.00) * (0.840182) +
    (-0.40 <= eta && eta < -0.30) * (13.00 <= pt * cosh(eta) && pt * cosh(eta) < 13.50) * (0.817839) +
    (-0.40 <= eta && eta < -0.30) * (13.50 <= pt * cosh(eta) && pt * cosh(eta) < 14.10) * (0.795040) +
    (-0.40 <= eta && eta < -0.30) * (14.10 <= pt * cosh(eta) && pt * cosh(eta) < 14.70) * (0.772450) +
    (-0.40 <= eta && eta < -0.30) * (14.70 <= pt * cosh(eta) && pt * cosh(eta) < 15.40) * (0.750764) +
    (-0.40 <= eta && eta < -0.30) * (15.40 <= pt * cosh(eta) && pt * cosh(eta) < 16.20) * (0.729121) +
    (-0.40 <= eta && eta < -0.30) * (16.20 <= pt * cosh(eta) && pt * cosh(eta) < 17.20) * (0.707458) +
    (-0.40 <= eta && eta < -0.30) * (17.20 <= pt * cosh(eta) && pt * cosh(eta) < 18.50) * (0.685485) +
    (-0.40 <= eta && eta < -0.30) * (18.50 <= pt * cosh(eta) && pt * cosh(eta) < 20.40) * (0.663109) +
    (-0.40 <= eta && eta < -0.30) * (20.40 <= pt * cosh(eta) && pt * cosh(eta) < 23.80) * (0.639895) +
    (-0.40 <= eta && eta < -0.30) * (23.80 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.608937) +
    (-0.30 <= eta && eta < -0.20) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.60) * (0.000000) +
    (-0.30 <= eta && eta < -0.20) * (0.60 <= pt * cosh(eta) && pt * cosh(eta) < 9.10) * (0.999403) +
    (-0.30 <= eta && eta < -0.20) * (9.10 <= pt * cosh(eta) && pt * cosh(eta) < 9.90) * (0.979783) +
    (-0.30 <= eta && eta < -0.20) * (9.90 <= pt * cosh(eta) && pt * cosh(eta) < 10.50) * (0.957542) +
    (-0.30 <= eta && eta < -0.20) * (10.50 <= pt * cosh(eta) && pt * cosh(eta) < 11.00) * (0.935103) +
    (-0.30 <= eta && eta < -0.20) * (11.00 <= pt * cosh(eta) && pt * cosh(eta) < 11.50) * (0.912162) +
    (-0.30 <= eta && eta < -0.20) * (11.50 <= pt * cosh(eta) && pt * cosh(eta) < 12.00) * (0.888083) +
    (-0.30 <= eta && eta < -0.20) * (12.00 <= pt * cosh(eta) && pt * cosh(eta) < 12.50) * (0.863888) +
    (-0.30 <= eta && eta < -0.20) * (12.50 <= pt * cosh(eta) && pt * cosh(eta) < 13.00) * (0.840353) +
    (-0.30 <= eta && eta < -0.20) * (13.00 <= pt * cosh(eta) && pt * cosh(eta) < 13.50) * (0.818007) +
    (-0.30 <= eta && eta < -0.20) * (13.50 <= pt * cosh(eta) && pt * cosh(eta) < 14.10) * (0.795199) +
    (-0.30 <= eta && eta < -0.20) * (14.10 <= pt * cosh(eta) && pt * cosh(eta) < 14.70) * (0.772600) +
    (-0.30 <= eta && eta < -0.20) * (14.70 <= pt * cosh(eta) && pt * cosh(eta) < 15.40) * (0.750902) +
    (-0.30 <= eta && eta < -0.20) * (15.40 <= pt * cosh(eta) && pt * cosh(eta) < 16.20) * (0.729244) +
    (-0.30 <= eta && eta < -0.20) * (16.20 <= pt * cosh(eta) && pt * cosh(eta) < 17.20) * (0.707564) +
    (-0.30 <= eta && eta < -0.20) * (17.20 <= pt * cosh(eta) && pt * cosh(eta) < 18.50) * (0.685572) +
    (-0.30 <= eta && eta < -0.20) * (18.50 <= pt * cosh(eta) && pt * cosh(eta) < 20.40) * (0.663176) +
    (-0.30 <= eta && eta < -0.20) * (20.40 <= pt * cosh(eta) && pt * cosh(eta) < 23.80) * (0.639938) +
    (-0.30 <= eta && eta < -0.20) * (23.80 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.608946) +
    (-0.20 <= eta && eta < -0.10) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.60) * (0.000000) +
    (-0.20 <= eta && eta < -0.10) * (0.60 <= pt * cosh(eta) && pt * cosh(eta) < 9.20) * (0.999359) +
    (-0.20 <= eta && eta < -0.10) * (9.20 <= pt * cosh(eta) && pt * cosh(eta) < 10.00) * (0.978520) +
    (-0.20 <= eta && eta < -0.10) * (10.00 <= pt * cosh(eta) && pt * cosh(eta) < 10.60) * (0.956649) +
    (-0.20 <= eta && eta < -0.10) * (10.60 <= pt * cosh(eta) && pt * cosh(eta) < 11.10) * (0.934232) +
    (-0.20 <= eta && eta < -0.10) * (11.10 <= pt * cosh(eta) && pt * cosh(eta) < 11.60) * (0.911399) +
    (-0.20 <= eta && eta < -0.10) * (11.60 <= pt * cosh(eta) && pt * cosh(eta) < 12.10) * (0.887482) +
    (-0.20 <= eta && eta < -0.10) * (12.10 <= pt * cosh(eta) && pt * cosh(eta) < 12.60) * (0.863473) +
    (-0.20 <= eta && eta < -0.10) * (12.60 <= pt * cosh(eta) && pt * cosh(eta) < 13.10) * (0.840125) +
    (-0.20 <= eta && eta < -0.10) * (13.10 <= pt * cosh(eta) && pt * cosh(eta) < 13.60) * (0.817953) +
    (-0.20 <= eta && eta < -0.10) * (13.60 <= pt * cosh(eta) && pt * cosh(eta) < 14.20) * (0.795311) +
    (-0.20 <= eta && eta < -0.10) * (14.20 <= pt * cosh(eta) && pt * cosh(eta) < 14.80) * (0.772859) +
    (-0.20 <= eta && eta < -0.10) * (14.80 <= pt * cosh(eta) && pt * cosh(eta) < 15.50) * (0.751280) +
    (-0.20 <= eta && eta < -0.10) * (15.50 <= pt * cosh(eta) && pt * cosh(eta) < 16.30) * (0.729714) +
    (-0.20 <= eta && eta < -0.10) * (16.30 <= pt * cosh(eta) && pt * cosh(eta) < 17.30) * (0.708095) +
    (-0.20 <= eta && eta < -0.10) * (17.30 <= pt * cosh(eta) && pt * cosh(eta) < 18.60) * (0.686125) +
    (-0.20 <= eta && eta < -0.10) * (18.60 <= pt * cosh(eta) && pt * cosh(eta) < 20.50) * (0.663701) +
    (-0.20 <= eta && eta < -0.10) * (20.50 <= pt * cosh(eta) && pt * cosh(eta) < 23.90) * (0.640366) +
    (-0.20 <= eta && eta < -0.10) * (23.90 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.609100) +
    (-0.10 <= eta && eta < -0.00) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.60) * (0.000000) +
    (-0.10 <= eta && eta < -0.00) * (0.60 <= pt * cosh(eta) && pt * cosh(eta) < 9.40) * (0.999439) +
    (-0.10 <= eta && eta < -0.00) * (9.40 <= pt * cosh(eta) && pt * cosh(eta) < 10.20) * (0.980491) +
    (-0.10 <= eta && eta < -0.00) * (10.20 <= pt * cosh(eta) && pt * cosh(eta) < 10.80) * (0.959593) +
    (-0.10 <= eta && eta < -0.00) * (10.80 <= pt * cosh(eta) && pt * cosh(eta) < 11.30) * (0.938333) +
    (-0.10 <= eta && eta < -0.00) * (11.30 <= pt * cosh(eta) && pt * cosh(eta) < 11.80) * (0.916456) +
    (-0.10 <= eta && eta < -0.00) * (11.80 <= pt * cosh(eta) && pt * cosh(eta) < 12.30) * (0.893323) +
    (-0.10 <= eta && eta < -0.00) * (12.30 <= pt * cosh(eta) && pt * cosh(eta) < 12.80) * (0.869897) +
    (-0.10 <= eta && eta < -0.00) * (12.80 <= pt * cosh(eta) && pt * cosh(eta) < 13.30) * (0.846930) +
    (-0.10 <= eta && eta < -0.00) * (13.30 <= pt * cosh(eta) && pt * cosh(eta) < 13.80) * (0.824955) +
    (-0.10 <= eta && eta < -0.00) * (13.80 <= pt * cosh(eta) && pt * cosh(eta) < 14.40) * (0.802353) +
    (-0.10 <= eta && eta < -0.00) * (14.40 <= pt * cosh(eta) && pt * cosh(eta) < 15.00) * (0.779781) +
    (-0.10 <= eta && eta < -0.00) * (15.00 <= pt * cosh(eta) && pt * cosh(eta) < 15.70) * (0.757937) +
    (-0.10 <= eta && eta < -0.00) * (15.70 <= pt * cosh(eta) && pt * cosh(eta) < 16.50) * (0.735958) +
    (-0.10 <= eta && eta < -0.00) * (16.50 <= pt * cosh(eta) && pt * cosh(eta) < 17.50) * (0.713770) +
    (-0.10 <= eta && eta < -0.00) * (17.50 <= pt * cosh(eta) && pt * cosh(eta) < 18.70) * (0.691898) +
    (-0.10 <= eta && eta < -0.00) * (18.70 <= pt * cosh(eta) && pt * cosh(eta) < 20.40) * (0.670134) +
    (-0.10 <= eta && eta < -0.00) * (20.40 <= pt * cosh(eta) && pt * cosh(eta) < 23.20) * (0.647473) +
    (-0.10 <= eta && eta < -0.00) * (23.20 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.610659) +
    (-0.00 <= eta && eta < 0.10) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.60) * (0.000000) +
    (-0.00 <= eta && eta < 0.10) * (0.60 <= pt * cosh(eta) && pt * cosh(eta) < 9.30) * (0.999354) +
    (-0.00 <= eta && eta < 0.10) * (9.30 <= pt * cosh(eta) && pt * cosh(eta) < 10.10) * (0.979103) +
    (-0.00 <= eta && eta < 0.10) * (10.10 <= pt * cosh(eta) && pt * cosh(eta) < 10.70) * (0.958182) +
    (-0.00 <= eta && eta < 0.10) * (10.70 <= pt * cosh(eta) && pt * cosh(eta) < 11.20) * (0.936357) +
    (-0.00 <= eta && eta < 0.10) * (11.20 <= pt * cosh(eta) && pt * cosh(eta) < 11.70) * (0.914009) +
    (-0.00 <= eta && eta < 0.10) * (11.70 <= pt * cosh(eta) && pt * cosh(eta) < 12.20) * (0.890487) +
    (-0.00 <= eta && eta < 0.10) * (12.20 <= pt * cosh(eta) && pt * cosh(eta) < 12.70) * (0.866768) +
    (-0.00 <= eta && eta < 0.10) * (12.70 <= pt * cosh(eta) && pt * cosh(eta) < 13.20) * (0.843606) +
    (-0.00 <= eta && eta < 0.10) * (13.20 <= pt * cosh(eta) && pt * cosh(eta) < 13.70) * (0.821527) +
    (-0.00 <= eta && eta < 0.10) * (13.70 <= pt * cosh(eta) && pt * cosh(eta) < 14.30) * (0.798897) +
    (-0.00 <= eta && eta < 0.10) * (14.30 <= pt * cosh(eta) && pt * cosh(eta) < 14.90) * (0.776376) +
    (-0.00 <= eta && eta < 0.10) * (14.90 <= pt * cosh(eta) && pt * cosh(eta) < 15.60) * (0.754655) +
    (-0.00 <= eta && eta < 0.10) * (15.60 <= pt * cosh(eta) && pt * cosh(eta) < 16.40) * (0.732873) +
    (-0.00 <= eta && eta < 0.10) * (16.40 <= pt * cosh(eta) && pt * cosh(eta) < 17.40) * (0.710960) +
    (-0.00 <= eta && eta < 0.10) * (17.40 <= pt * cosh(eta) && pt * cosh(eta) < 18.70) * (0.688609) +
    (-0.00 <= eta && eta < 0.10) * (18.70 <= pt * cosh(eta) && pt * cosh(eta) < 20.50) * (0.666272) +
    (-0.00 <= eta && eta < 0.10) * (20.50 <= pt * cosh(eta) && pt * cosh(eta) < 23.60) * (0.643430) +
    (-0.00 <= eta && eta < 0.10) * (23.60 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.609805) +
    (0.10 <= eta && eta < 0.20) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.60) * (0.000000) +
    (0.10 <= eta && eta < 0.20) * (0.60 <= pt * cosh(eta) && pt * cosh(eta) < 9.00) * (0.999442) +
    (0.10 <= eta && eta < 0.20) * (9.00 <= pt * cosh(eta) && pt * cosh(eta) < 9.80) * (0.980896) +
    (0.10 <= eta && eta < 0.20) * (9.80 <= pt * cosh(eta) && pt * cosh(eta) < 10.40) * (0.959254) +
    (0.10 <= eta && eta < 0.20) * (10.40 <= pt * cosh(eta) && pt * cosh(eta) < 10.90) * (0.936937) +
    (0.10 <= eta && eta < 0.20) * (10.90 <= pt * cosh(eta) && pt * cosh(eta) < 11.40) * (0.914015) +
    (0.10 <= eta && eta < 0.20) * (11.40 <= pt * cosh(eta) && pt * cosh(eta) < 11.90) * (0.889852) +
    (0.10 <= eta && eta < 0.20) * (11.90 <= pt * cosh(eta) && pt * cosh(eta) < 12.40) * (0.865505) +
    (0.10 <= eta && eta < 0.20) * (12.40 <= pt * cosh(eta) && pt * cosh(eta) < 12.90) * (0.841780) +
    (0.10 <= eta && eta < 0.20) * (12.90 <= pt * cosh(eta) && pt * cosh(eta) < 13.40) * (0.819231) +
    (0.10 <= eta && eta < 0.20) * (13.40 <= pt * cosh(eta) && pt * cosh(eta) < 14.00) * (0.796203) +
    (0.10 <= eta && eta < 0.20) * (14.00 <= pt * cosh(eta) && pt * cosh(eta) < 14.60) * (0.773383) +
    (0.10 <= eta && eta < 0.20) * (14.60 <= pt * cosh(eta) && pt * cosh(eta) < 15.30) * (0.751478) +
    (0.10 <= eta && eta < 0.20) * (15.30 <= pt * cosh(eta) && pt * cosh(eta) < 16.10) * (0.729625) +
    (0.10 <= eta && eta < 0.20) * (16.10 <= pt * cosh(eta) && pt * cosh(eta) < 17.10) * (0.707770) +
    (0.10 <= eta && eta < 0.20) * (17.10 <= pt * cosh(eta) && pt * cosh(eta) < 18.40) * (0.685626) +
    (0.10 <= eta && eta < 0.20) * (18.40 <= pt * cosh(eta) && pt * cosh(eta) < 20.20) * (0.663668) +
    (0.10 <= eta && eta < 0.20) * (20.20 <= pt * cosh(eta) && pt * cosh(eta) < 23.50) * (0.640799) +
    (0.10 <= eta && eta < 0.20) * (23.50 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.609016) +
    (0.20 <= eta && eta < 0.30) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.60) * (0.000000) +
    (0.20 <= eta && eta < 0.30) * (0.60 <= pt * cosh(eta) && pt * cosh(eta) < 8.90) * (0.999456) +
    (0.20 <= eta && eta < 0.30) * (8.90 <= pt * cosh(eta) && pt * cosh(eta) < 9.70) * (0.980442) +
    (0.20 <= eta && eta < 0.30) * (9.70 <= pt * cosh(eta) && pt * cosh(eta) < 10.20) * (0.959342) +
    (0.20 <= eta && eta < 0.30) * (10.20 <= pt * cosh(eta) && pt * cosh(eta) < 10.70) * (0.938527) +
    (0.20 <= eta && eta < 0.30) * (10.70 <= pt * cosh(eta) && pt * cosh(eta) < 11.20) * (0.915372) +
    (0.20 <= eta && eta < 0.30) * (11.20 <= pt * cosh(eta) && pt * cosh(eta) < 11.70) * (0.890867) +
    (0.20 <= eta && eta < 0.30) * (11.70 <= pt * cosh(eta) && pt * cosh(eta) < 12.20) * (0.866131) +
    (0.20 <= eta && eta < 0.30) * (12.20 <= pt * cosh(eta) && pt * cosh(eta) < 12.70) * (0.842017) +
    (0.20 <= eta && eta < 0.30) * (12.70 <= pt * cosh(eta) && pt * cosh(eta) < 13.20) * (0.819110) +
    (0.20 <= eta && eta < 0.30) * (13.20 <= pt * cosh(eta) && pt * cosh(eta) < 13.70) * (0.797763) +
    (0.20 <= eta && eta < 0.30) * (13.70 <= pt * cosh(eta) && pt * cosh(eta) < 14.30) * (0.776311) +
    (0.20 <= eta && eta < 0.30) * (14.30 <= pt * cosh(eta) && pt * cosh(eta) < 15.00) * (0.753723) +
    (0.20 <= eta && eta < 0.30) * (15.00 <= pt * cosh(eta) && pt * cosh(eta) < 15.80) * (0.731218) +
    (0.20 <= eta && eta < 0.30) * (15.80 <= pt * cosh(eta) && pt * cosh(eta) < 16.80) * (0.708760) +
    (0.20 <= eta && eta < 0.30) * (16.80 <= pt * cosh(eta) && pt * cosh(eta) < 18.10) * (0.686083) +
    (0.20 <= eta && eta < 0.30) * (18.10 <= pt * cosh(eta) && pt * cosh(eta) < 19.90) * (0.663704) +
    (0.20 <= eta && eta < 0.30) * (19.90 <= pt * cosh(eta) && pt * cosh(eta) < 23.10) * (0.640869) +
    (0.20 <= eta && eta < 0.30) * (23.10 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.608858) +
    (0.30 <= eta && eta < 0.40) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.60) * (0.000000) +
    (0.30 <= eta && eta < 0.40) * (0.60 <= pt * cosh(eta) && pt * cosh(eta) < 8.90) * (0.999408) +
    (0.30 <= eta && eta < 0.40) * (8.90 <= pt * cosh(eta) && pt * cosh(eta) < 9.70) * (0.980163) +
    (0.30 <= eta && eta < 0.40) * (9.70 <= pt * cosh(eta) && pt * cosh(eta) < 10.30) * (0.959038) +
    (0.30 <= eta && eta < 0.40) * (10.30 <= pt * cosh(eta) && pt * cosh(eta) < 10.80) * (0.936756) +
    (0.30 <= eta && eta < 0.40) * (10.80 <= pt * cosh(eta) && pt * cosh(eta) < 11.30) * (0.913584) +
    (0.30 <= eta && eta < 0.40) * (11.30 <= pt * cosh(eta) && pt * cosh(eta) < 11.80) * (0.889168) +
    (0.30 <= eta && eta < 0.40) * (11.80 <= pt * cosh(eta) && pt * cosh(eta) < 12.30) * (0.864590) +
    (0.30 <= eta && eta < 0.40) * (12.30 <= pt * cosh(eta) && pt * cosh(eta) < 12.80) * (0.840671) +
    (0.30 <= eta && eta < 0.40) * (12.80 <= pt * cosh(eta) && pt * cosh(eta) < 13.30) * (0.817972) +
    (0.30 <= eta && eta < 0.40) * (13.30 <= pt * cosh(eta) && pt * cosh(eta) < 13.90) * (0.794832) +
    (0.30 <= eta && eta < 0.40) * (13.90 <= pt * cosh(eta) && pt * cosh(eta) < 14.50) * (0.771941) +
    (0.30 <= eta && eta < 0.40) * (14.50 <= pt * cosh(eta) && pt * cosh(eta) < 15.20) * (0.750012) +
    (0.30 <= eta && eta < 0.40) * (15.20 <= pt * cosh(eta) && pt * cosh(eta) < 16.00) * (0.728180) +
    (0.30 <= eta && eta < 0.40) * (16.00 <= pt * cosh(eta) && pt * cosh(eta) < 17.00) * (0.706396) +
    (0.30 <= eta && eta < 0.40) * (17.00 <= pt * cosh(eta) && pt * cosh(eta) < 18.30) * (0.684379) +
    (0.30 <= eta && eta < 0.40) * (18.30 <= pt * cosh(eta) && pt * cosh(eta) < 20.20) * (0.662061) +
    (0.30 <= eta && eta < 0.40) * (20.20 <= pt * cosh(eta) && pt * cosh(eta) < 23.70) * (0.638751) +
    (0.30 <= eta && eta < 0.40) * (23.70 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.608560) +
    (0.40 <= eta && eta < 0.50) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.60) * (0.000000) +
    (0.40 <= eta && eta < 0.50) * (0.60 <= pt * cosh(eta) && pt * cosh(eta) < 8.90) * (0.999373) +
    (0.40 <= eta && eta < 0.50) * (8.90 <= pt * cosh(eta) && pt * cosh(eta) < 9.70) * (0.978769) +
    (0.40 <= eta && eta < 0.50) * (9.70 <= pt * cosh(eta) && pt * cosh(eta) < 10.20) * (0.957448) +
    (0.40 <= eta && eta < 0.50) * (10.20 <= pt * cosh(eta) && pt * cosh(eta) < 10.70) * (0.936499) +
    (0.40 <= eta && eta < 0.50) * (10.70 <= pt * cosh(eta) && pt * cosh(eta) < 11.20) * (0.913062) +
    (0.40 <= eta && eta < 0.50) * (11.20 <= pt * cosh(eta) && pt * cosh(eta) < 11.70) * (0.888382) +
    (0.40 <= eta && eta < 0.50) * (11.70 <= pt * cosh(eta) && pt * cosh(eta) < 12.20) * (0.863567) +
    (0.40 <= eta && eta < 0.50) * (12.20 <= pt * cosh(eta) && pt * cosh(eta) < 12.70) * (0.839454) +
    (0.40 <= eta && eta < 0.50) * (12.70 <= pt * cosh(eta) && pt * cosh(eta) < 13.20) * (0.816607) +
    (0.40 <= eta && eta < 0.50) * (13.20 <= pt * cosh(eta) && pt * cosh(eta) < 13.70) * (0.795361) +
    (0.40 <= eta && eta < 0.50) * (13.70 <= pt * cosh(eta) && pt * cosh(eta) < 14.30) * (0.774050) +
    (0.40 <= eta && eta < 0.50) * (14.30 <= pt * cosh(eta) && pt * cosh(eta) < 15.00) * (0.751650) +
    (0.40 <= eta && eta < 0.50) * (15.00 <= pt * cosh(eta) && pt * cosh(eta) < 15.80) * (0.729367) +
    (0.40 <= eta && eta < 0.50) * (15.80 <= pt * cosh(eta) && pt * cosh(eta) < 16.80) * (0.707163) +
    (0.40 <= eta && eta < 0.50) * (16.80 <= pt * cosh(eta) && pt * cosh(eta) < 18.10) * (0.684773) +
    (0.40 <= eta && eta < 0.50) * (18.10 <= pt * cosh(eta) && pt * cosh(eta) < 19.90) * (0.662704) +
    (0.40 <= eta && eta < 0.50) * (19.90 <= pt * cosh(eta) && pt * cosh(eta) < 23.20) * (0.639904) +
    (0.40 <= eta && eta < 0.50) * (23.20 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.608646) +
    (0.50 <= eta && eta < 0.60) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.60) * (0.000000) +
    (0.50 <= eta && eta < 0.60) * (0.60 <= pt * cosh(eta) && pt * cosh(eta) < 9.00) * (0.999417) +
    (0.50 <= eta && eta < 0.60) * (9.00 <= pt * cosh(eta) && pt * cosh(eta) < 9.80) * (0.979539) +
    (0.50 <= eta && eta < 0.60) * (9.80 <= pt * cosh(eta) && pt * cosh(eta) < 10.40) * (0.957986) +
    (0.50 <= eta && eta < 0.60) * (10.40 <= pt * cosh(eta) && pt * cosh(eta) < 10.90) * (0.935548) +
    (0.50 <= eta && eta < 0.60) * (10.90 <= pt * cosh(eta) && pt * cosh(eta) < 11.40) * (0.912442) +
    (0.50 <= eta && eta < 0.60) * (11.40 <= pt * cosh(eta) && pt * cosh(eta) < 11.90) * (0.888165) +
    (0.50 <= eta && eta < 0.60) * (11.90 <= pt * cosh(eta) && pt * cosh(eta) < 12.40) * (0.863766) +
    (0.50 <= eta && eta < 0.60) * (12.40 <= pt * cosh(eta) && pt * cosh(eta) < 12.90) * (0.840041) +
    (0.50 <= eta && eta < 0.60) * (12.90 <= pt * cosh(eta) && pt * cosh(eta) < 13.40) * (0.817531) +
    (0.50 <= eta && eta < 0.60) * (13.40 <= pt * cosh(eta) && pt * cosh(eta) < 14.00) * (0.794579) +
    (0.50 <= eta && eta < 0.60) * (14.00 <= pt * cosh(eta) && pt * cosh(eta) < 14.60) * (0.771863) +
    (0.50 <= eta && eta < 0.60) * (14.60 <= pt * cosh(eta) && pt * cosh(eta) < 15.30) * (0.750084) +
    (0.50 <= eta && eta < 0.60) * (15.30 <= pt * cosh(eta) && pt * cosh(eta) < 16.10) * (0.728379) +
    (0.50 <= eta && eta < 0.60) * (16.10 <= pt * cosh(eta) && pt * cosh(eta) < 17.10) * (0.706692) +
    (0.50 <= eta && eta < 0.60) * (17.10 <= pt * cosh(eta) && pt * cosh(eta) < 18.40) * (0.684739) +
    (0.50 <= eta && eta < 0.60) * (18.40 <= pt * cosh(eta) && pt * cosh(eta) < 20.30) * (0.662439) +
    (0.50 <= eta && eta < 0.60) * (20.30 <= pt * cosh(eta) && pt * cosh(eta) < 23.80) * (0.639080) +
    (0.50 <= eta && eta < 0.60) * (23.80 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.608689) +
    (0.60 <= eta && eta < 0.70) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.60) * (0.000000) +
    (0.60 <= eta && eta < 0.70) * (0.60 <= pt * cosh(eta) && pt * cosh(eta) < 9.10) * (0.999335) +
    (0.60 <= eta && eta < 0.70) * (9.10 <= pt * cosh(eta) && pt * cosh(eta) < 9.90) * (0.978710) +
    (0.60 <= eta && eta < 0.70) * (9.90 <= pt * cosh(eta) && pt * cosh(eta) < 10.50) * (0.956332) +
    (0.60 <= eta && eta < 0.70) * (10.50 <= pt * cosh(eta) && pt * cosh(eta) < 11.00) * (0.933605) +
    (0.60 <= eta && eta < 0.70) * (11.00 <= pt * cosh(eta) && pt * cosh(eta) < 11.50) * (0.910475) +
    (0.60 <= eta && eta < 0.70) * (11.50 <= pt * cosh(eta) && pt * cosh(eta) < 12.00) * (0.886279) +
    (0.60 <= eta && eta < 0.70) * (12.00 <= pt * cosh(eta) && pt * cosh(eta) < 12.50) * (0.862034) +
    (0.60 <= eta && eta < 0.70) * (12.50 <= pt * cosh(eta) && pt * cosh(eta) < 13.00) * (0.838503) +
    (0.60 <= eta && eta < 0.70) * (13.00 <= pt * cosh(eta) && pt * cosh(eta) < 13.50) * (0.816201) +
    (0.60 <= eta && eta < 0.70) * (13.50 <= pt * cosh(eta) && pt * cosh(eta) < 14.10) * (0.793474) +
    (0.60 <= eta && eta < 0.70) * (14.10 <= pt * cosh(eta) && pt * cosh(eta) < 14.70) * (0.770986) +
    (0.60 <= eta && eta < 0.70) * (14.70 <= pt * cosh(eta) && pt * cosh(eta) < 15.40) * (0.749421) +
    (0.60 <= eta && eta < 0.70) * (15.40 <= pt * cosh(eta) && pt * cosh(eta) < 16.20) * (0.727920) +
    (0.60 <= eta && eta < 0.70) * (16.20 <= pt * cosh(eta) && pt * cosh(eta) < 17.20) * (0.706419) +
    (0.60 <= eta && eta < 0.70) * (17.20 <= pt * cosh(eta) && pt * cosh(eta) < 18.50) * (0.684629) +
    (0.60 <= eta && eta < 0.70) * (18.50 <= pt * cosh(eta) && pt * cosh(eta) < 20.40) * (0.662458) +
    (0.60 <= eta && eta < 0.70) * (20.40 <= pt * cosh(eta) && pt * cosh(eta) < 23.90) * (0.639180) +
    (0.60 <= eta && eta < 0.70) * (23.90 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.608771) +
    (0.70 <= eta && eta < 0.80) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.60) * (0.000000) +
    (0.70 <= eta && eta < 0.80) * (0.60 <= pt * cosh(eta) && pt * cosh(eta) < 9.00) * (0.999463) +
    (0.70 <= eta && eta < 0.80) * (9.00 <= pt * cosh(eta) && pt * cosh(eta) < 9.80) * (0.981050) +
    (0.70 <= eta && eta < 0.80) * (9.80 <= pt * cosh(eta) && pt * cosh(eta) < 10.40) * (0.959830) +
    (0.70 <= eta && eta < 0.80) * (10.40 <= pt * cosh(eta) && pt * cosh(eta) < 10.90) * (0.937687) +
    (0.70 <= eta && eta < 0.80) * (10.90 <= pt * cosh(eta) && pt * cosh(eta) < 11.40) * (0.914867) +
    (0.70 <= eta && eta < 0.80) * (11.40 <= pt * cosh(eta) && pt * cosh(eta) < 11.90) * (0.890769) +
    (0.70 <= eta && eta < 0.80) * (11.90 <= pt * cosh(eta) && pt * cosh(eta) < 12.40) * (0.866452) +
    (0.70 <= eta && eta < 0.80) * (12.40 <= pt * cosh(eta) && pt * cosh(eta) < 12.90) * (0.842728) +
    (0.70 <= eta && eta < 0.80) * (12.90 <= pt * cosh(eta) && pt * cosh(eta) < 13.40) * (0.820157) +
    (0.70 <= eta && eta < 0.80) * (13.40 <= pt * cosh(eta) && pt * cosh(eta) < 14.00) * (0.797090) +
    (0.70 <= eta && eta < 0.80) * (14.00 <= pt * cosh(eta) && pt * cosh(eta) < 14.60) * (0.774214) +
    (0.70 <= eta && eta < 0.80) * (14.60 <= pt * cosh(eta) && pt * cosh(eta) < 15.30) * (0.752241) +
    (0.70 <= eta && eta < 0.80) * (15.30 <= pt * cosh(eta) && pt * cosh(eta) < 16.10) * (0.730308) +
    (0.70 <= eta && eta < 0.80) * (16.10 <= pt * cosh(eta) && pt * cosh(eta) < 17.10) * (0.708360) +
    (0.70 <= eta && eta < 0.80) * (17.10 <= pt * cosh(eta) && pt * cosh(eta) < 18.40) * (0.686112) +
    (0.70 <= eta && eta < 0.80) * (18.40 <= pt * cosh(eta) && pt * cosh(eta) < 20.20) * (0.664041) +
    (0.70 <= eta && eta < 0.80) * (20.20 <= pt * cosh(eta) && pt * cosh(eta) < 23.40) * (0.641356) +
    (0.70 <= eta && eta < 0.80) * (23.40 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.609149) +
    (0.80 <= eta && eta < 0.90) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.60) * (0.000000) +
    (0.80 <= eta && eta < 0.90) * (0.60 <= pt * cosh(eta) && pt * cosh(eta) < 9.20) * (0.999454) +
    (0.80 <= eta && eta < 0.90) * (9.20 <= pt * cosh(eta) && pt * cosh(eta) < 10.10) * (0.979938) +
    (0.80 <= eta && eta < 0.90) * (10.10 <= pt * cosh(eta) && pt * cosh(eta) < 10.70) * (0.957951) +
    (0.80 <= eta && eta < 0.90) * (10.70 <= pt * cosh(eta) && pt * cosh(eta) < 11.20) * (0.936071) +
    (0.80 <= eta && eta < 0.90) * (11.20 <= pt * cosh(eta) && pt * cosh(eta) < 11.70) * (0.913686) +
    (0.80 <= eta && eta < 0.90) * (11.70 <= pt * cosh(eta) && pt * cosh(eta) < 12.20) * (0.890140) +
    (0.80 <= eta && eta < 0.90) * (12.20 <= pt * cosh(eta) && pt * cosh(eta) < 12.70) * (0.866410) +
    (0.80 <= eta && eta < 0.90) * (12.70 <= pt * cosh(eta) && pt * cosh(eta) < 13.20) * (0.843247) +
    (0.80 <= eta && eta < 0.90) * (13.20 <= pt * cosh(eta) && pt * cosh(eta) < 13.70) * (0.821175) +
    (0.80 <= eta && eta < 0.90) * (13.70 <= pt * cosh(eta) && pt * cosh(eta) < 14.30) * (0.798559) +
    (0.80 <= eta && eta < 0.90) * (14.30 <= pt * cosh(eta) && pt * cosh(eta) < 14.90) * (0.776059) +
    (0.80 <= eta && eta < 0.90) * (14.90 <= pt * cosh(eta) && pt * cosh(eta) < 15.60) * (0.754363) +
    (0.80 <= eta && eta < 0.90) * (15.60 <= pt * cosh(eta) && pt * cosh(eta) < 16.40) * (0.732610) +
    (0.80 <= eta && eta < 0.90) * (16.40 <= pt * cosh(eta) && pt * cosh(eta) < 17.40) * (0.710732) +
    (0.80 <= eta && eta < 0.90) * (17.40 <= pt * cosh(eta) && pt * cosh(eta) < 18.70) * (0.688420) +
    (0.80 <= eta && eta < 0.90) * (18.70 <= pt * cosh(eta) && pt * cosh(eta) < 20.50) * (0.666126) +
    (0.80 <= eta && eta < 0.90) * (20.50 <= pt * cosh(eta) && pt * cosh(eta) < 23.60) * (0.643332) +
    (0.80 <= eta && eta < 0.90) * (23.60 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.609784) +
    (0.90 <= eta && eta < 1.00) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.60) * (0.000000) +
    (0.90 <= eta && eta < 1.00) * (0.60 <= pt * cosh(eta) && pt * cosh(eta) < 9.30) * (0.999390) +
    (0.90 <= eta && eta < 1.00) * (9.30 <= pt * cosh(eta) && pt * cosh(eta) < 10.10) * (0.979355) +
    (0.90 <= eta && eta < 1.00) * (10.10 <= pt * cosh(eta) && pt * cosh(eta) < 10.70) * (0.957686) +
    (0.90 <= eta && eta < 1.00) * (10.70 <= pt * cosh(eta) && pt * cosh(eta) < 11.20) * (0.935743) +
    (0.90 <= eta && eta < 1.00) * (11.20 <= pt * cosh(eta) && pt * cosh(eta) < 11.70) * (0.913316) +
    (0.90 <= eta && eta < 1.00) * (11.70 <= pt * cosh(eta) && pt * cosh(eta) < 12.20) * (0.889743) +
    (0.90 <= eta && eta < 1.00) * (12.20 <= pt * cosh(eta) && pt * cosh(eta) < 12.70) * (0.866001) +
    (0.90 <= eta && eta < 1.00) * (12.70 <= pt * cosh(eta) && pt * cosh(eta) < 13.20) * (0.842837) +
    (0.90 <= eta && eta < 1.00) * (13.20 <= pt * cosh(eta) && pt * cosh(eta) < 13.70) * (0.820773) +
    (0.90 <= eta && eta < 1.00) * (13.70 <= pt * cosh(eta) && pt * cosh(eta) < 14.30) * (0.798174) +
    (0.90 <= eta && eta < 1.00) * (14.30 <= pt * cosh(eta) && pt * cosh(eta) < 14.90) * (0.775697) +
    (0.90 <= eta && eta < 1.00) * (14.90 <= pt * cosh(eta) && pt * cosh(eta) < 15.60) * (0.754029) +
    (0.90 <= eta && eta < 1.00) * (15.60 <= pt * cosh(eta) && pt * cosh(eta) < 16.40) * (0.732311) +
    (0.90 <= eta && eta < 1.00) * (16.40 <= pt * cosh(eta) && pt * cosh(eta) < 17.40) * (0.710471) +
    (0.90 <= eta && eta < 1.00) * (17.40 <= pt * cosh(eta) && pt * cosh(eta) < 18.70) * (0.688204) +
    (0.90 <= eta && eta < 1.00) * (18.70 <= pt * cosh(eta) && pt * cosh(eta) < 20.50) * (0.665960) +
    (0.90 <= eta && eta < 1.00) * (20.50 <= pt * cosh(eta) && pt * cosh(eta) < 23.60) * (0.643220) +
    (0.90 <= eta && eta < 1.00) * (23.60 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.609760) +
    (1.00 <= eta && eta < 1.10) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.000000) }

# Everything else is not ID'd at all.
  add EfficiencyFormula {0} {0} { 0.00 }

}


#source dualRICH_aerogel_0.0mrad.tcl
source dualRICH_aerogel_0.5mrad.tcl
source dualRICH_c2f6_0.0mrad.tcl 
#source dualRICH_c2f6_0.0mrad.tcl



##################
# ROOT tree writer
##################

# tracks, towers and eflow objects are not stored by default in the output.
# if needed (for jet constituent or other studies), uncomment the relevant
# "add Branch ..." lines.

module TreeWriter TreeWriter {
# add Branch InputArray BranchName BranchClass
  add Branch Delphes/allParticles Particle GenParticle
  add Branch BeamSpotFilter/beamSpotParticle BeamSpot GenParticle

  add Branch TrackMerger/tracks Track Track
  add Branch Calorimeter/towers Tower Tower

  add Branch HCal/eflowTracks EFlowTrack Track
  add Branch ECal/eflowPhotons EFlowPhoton Tower
  add Branch HCal/eflowNeutralHadrons EFlowNeutralHadron Tower

  add Branch pfRICH/tracks pfRICHTrack Track
  add Branch barrelDIRC/tracks barrelDIRCTrack Track
  add Branch dualRICH_aerogel/tracks dualRICHagTrack Track
  add Branch dualRICH_c2f6/tracks dualRICHcfTrack Track

  add Branch GenJetFinder/jets GenJet Jet
  add Branch GenMissingET/momentum GenMissingET MissingET

  add Branch UniqueObjectFinder/jets Jet Jet
  add Branch UniqueObjectFinder/electrons Electron Electron
  add Branch UniqueObjectFinder/photons Photon Photon

  add Branch MissingET/momentum MissingET MissingET
  add Branch ScalarHT/energy ScalarHT ScalarHT
}
