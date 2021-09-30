# BHNLAnalyser

## Installation
Note: 
* it is preferable to install and use the tool in a bash environment
* the tool is installed within the combine framework

Set up the environment
```
cmsrel CMSSW_10_2_15
cd CMSSW_10_2_15/src
cmsenv
git cms-init
```

Install combine
```
git clone https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit
cd HiggsAnalysis/CombinedLimit
source env_standalone.sh 
make -j 8; make 
```

Install the BHNLAnalyser tool
```
git clone git@github.com:BParkHNLs/BHNLAnalyser.git
```


### After first installation
```
cd CMSSW_10_2_15/src/HiggsAnalysis/CombinedLimit/BHNAnalyser
cmsenv
```

If need to use combine:
```
cd ..
source env_standalone.sh
cd BHNAnalyser
```
