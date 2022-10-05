#!/bin/bash

cmsenv
cd ..
source env_standalone.sh &> tmp
rm tmp
cd BHNLAnalyser
