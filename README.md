# SVJ s-channel Analysis - Signal Request Framework

Job options and condor scripts for the signal request

## Condor submit
1. setupATLAS (won't work if you don't do this locally first!)
2. Open condor/submitCondor.py. Make sure you are running over the desired DSIDs (see for loop). Also set number of events
3. Commands for generation and truth DAOD creation are in the condor\_evnt.sh script. Make sure the release information is correct.
4. python submitCondor.py

## DAOD evaluation
Navigate to the TruthAnalysisTests directory 
1. setupATLAS (if not already done)
2. mkdir build
3. cd source
4. asetup AnalysisBase,21.2.156,here . Check that CMakeLists file is created
5. cd ../build
7. cmake ../source
8. make
9. source x86\*/setup.sh
10. cd ../run
11. TruthDerivationTester --input <path_to_input_DAOD> --output output.root --nevents -1
12. Use myPlotter.C to analyze histogram files

To modify the code, make changes to the source/MyAnalysis/util/TruthDerivationTester.cxx. After making changes, navigate to the build directory, and run `make`. Then run `source x86*/setup.sh`. Your changes should now be fully compiled and you are ready to run the executable again.

Examples of running the executable are in the extra\_jets.sh bash script.
