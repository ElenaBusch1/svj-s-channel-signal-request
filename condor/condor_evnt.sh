#!/bin/bash

echo "RUNNING..."

# Set up ATLAS on this machine:
source ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh

export RUN_NUMBER=$1
export N_EVENTS=$2
export PWD=`pwd`
export STREND=".root"

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
Gen_tf.py --ecmEnergy=13000. --firstEvent=1  --maxEvents=$N_EVENTS --randomSeed=111 --jobConfig=$RUN_NUMBER --outputEVNTFile=$RUN_NUMBER$STREND > $RUN_NUMBER/log.generate
)

#2 make the TRUTH DAOD 
#(
#echo "GETTING TRUTH DAOD..."
#asetup AthDerivation,21.2.132.0
#Reco_tf.py --inputEVNTFile $RUN_NUMBER$STREND --outputDAODFile $RUN_NUMBER$STREND --reductionConf TRUTH1
#)

##3 make the histogram file
#(
#echo "MAKING HISTOGRAMS..."
#asetup AnalysisBase,21.2.156
#TruthDerivationTester --input DAOD_TRUTH3.$RUN_NUMBER$STREND --output hists_$RUN_NUMBER$STREND --nevents -1
#)

#4 Simulate
#(
#asetup AthSimulation,21.0.136
#Sim_tf.py --inputEvgenFile $RUN_NUMBER$STREND --outputHITSFile HITS.root  --maxEvents $N_EVENTS --physicsList 'FTFP_BERT_ATL_VALIDATION' --truthStrategy 'MC15aPlusLLP' --simulator 'FullG4' --DBRelease 'all:current' --conditionsTag 'default:OFLCOND-MC16-SDR-14' --DataRunNumber '284500' --preExec 'EVNTtoHITS:simFlags.SimBarcodeOffset.set_Value_and_Lock(200000)' 'EVNTtoHITS:simFlags.TRTRangeCut=30.0;simFlags.TightMuonStepping=True' --preInclude 'EVNTtoHITS:SimulationJobOptions/preInclude.BeamPipeKill.py,SimulationJobOptions/preInclude.FrozenShowersFCalOnly.py' --geometryVersion 'default:ATLAS-R2-2016-01-00-01_VALIDATION' --postInclude 'default:RecJobTransforms/UseFrontier.py'
#
## Reconstruct to Get Truth DAOD
#(
#echo "GETTING TRUTH DAOD..."
#asetup 21.2.86.0,AthDerivation
#Reco_tf.py --inputEVNTFile EVNT.root --outputDAODFile Hino.root --reductionConf TRUTH1
#)

