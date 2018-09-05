# CTPPSNTupler
Based on https://github.com/forthommel/DiffractiveForwardAnalysis. Indeed, this code is a copy such package, however including CTPPS pixels, strips and timing containers. 
In addition Jets, Particle Flow for Missing Mass Searches.

```sh
cmsrel CMSSW_9_4_4
cd CMSSW_9_4_4/src
cmsenv
git clone git@github.com:ForwardGroupBrazil/CTPPSNTupler.git
cd NTuplizer 
cd CMSSW_9_4_4/src
git cms-init 
git cms-merge-topic lsoffi:CMSSW_9_4_0_pre3_TnP 
git cms-merge-topic guitargeek:ElectronID_MVA2017_940pre3 
scram b -j8 
cd $CMSSW_BASE/src 
git clone https://github.com/lsoffi/RecoEgamma-PhotonIdentification.git data/RecoEgamma/PhotonIdentification/data 
cd data/RecoEgamma/PhotonIdentification/data 
git checkout CMSSW_9_4_0_pre3_TnP 
cd $CMSSW_BASE/src 
git clone https://github.com/lsoffi/RecoEgamma-ElectronIdentification.git data/RecoEgamma/ElectronIdentification/data 
cd data/RecoEgamma/ElectronIdentification/data 
git checkout CMSSW_9_4_0_pre3_TnP 
cd $CMSSW_BASE/src 
git cms-merge-topic cms-egamma:EGM_94X_v1 
cd EgammaAnalysis/ElectronTools/data 
git clone https://github.com/ECALELFS/ScalesSmearings.git 
cd ScalesSmearings/ 
git checkout Run2017_17Nov2017_v1
cd $CMSSW_BASE/src
scram b -j8 
```
