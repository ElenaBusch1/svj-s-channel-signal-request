#!/bin/bash

echo "RUNNING..."

# Set up ATLAS on this machine:
source ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh

export RUN_NUMBER=$1
export N_EVENTS=$2
export PWD=`pwd`
export STREND=".root"
export LOGEND=".log.generate"

mkdir $RUN_NUMBER
#mkdir run
cp mc*.py $RUN_NUMBER
#cd run

#1 make the EVNT
(
echo "GENERATING..."
echo "Contents of DSID folder: "
ls $RUN_NUMBER
asetup AthGeneration,21.6.96
echo "Gen_tf.py --ecmEnergy=13000. --firstEvent=1  --maxEvents=$N_EVENTS --randomSeed=111 --jobConfig=$RUN_NUMBER --outputEVNTFile=$RUN_NUMBER$STREND"
Gen_tf.py --ecmEnergy=13000. --firstEvent=1  --maxEvents=$N_EVENTS --randomSeed=111 --jobConfig=$RUN_NUMBER --outputEVNTFile=$RUN_NUMBER$STREND > $RUN_NUMBER/$RUN_NUMBER$LOGEND
)

#2 make the TRUTH DAOD 
(
echo "GETTING TRUTH DAOD..."
asetup AthDerivation,21.2.132.0
Reco_tf.py --inputEVNTFile $RUN_NUMBER$STREND --outputDAODFile $RUN_NUMBER$STREND --reductionConf TRUTH3
)

##3 make the histogram file
#(
#echo "MAKING HISTOGRAMS..."
#asetup AnalysisBase,21.2.156
#TruthDerivationTester --input DAOD_TRUTH3.$RUN_NUMBER$STREND --output hists_$RUN_NUMBER$STREND --nevents -1
#)

#--------------------------------------------------------------
# Generate (already done)
#(
#asetup AthGeneration,21.6.48
#Gen_tf.py --ecmEnergy 13000 --firstEvent 1 --maxEvents $N_EVENTS --randomSeed 10041992 --jobConfig $RUN_NUMBER --outputEVNTFile EVNT.root
#)
#
## Reconstruct to Get Truth DAOD
#(
#echo "GETTING TRUTH DAOD..."
#asetup 21.2.86.0,AthDerivation
#Reco_tf.py --inputEVNTFile EVNT.root --outputDAODFile Hino.root --reductionConf TRUTH1
#)

