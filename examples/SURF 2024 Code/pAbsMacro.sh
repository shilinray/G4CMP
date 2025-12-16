source ~/Desktop/otherG4cmp/install/share/G4CMP/g4cmp_env.sh
source ~/Desktop/myG4/install/bin/geant4.sh

for probMount in $(seq 0.0 0.02 1.0)
do 
  sed -i '' "s/double pAbsMount = [0-9]*\.*[0-9]*/double pAbsMount = $probMount/" "/Users/matthewcantor/Desktop/otherG4cmp/work/singleKIDSim/singleKIDSim/src/RISQTutorialDetectorConstruction.cc"
  for probSilicon in $(seq 0.001 0.002 0.011)
  do
    hitFile="RISQTutorial_hits_Si${probSilicon}Mount${probMount}.txt"
    primaryFile="RISQTutorial_primary_Si${probSilicon}Mount${probMount}.txt"

    sed -i '' "s/Hit_file(getenv(\"G4CMP_HIT_FILE\")?getenv(\"G4CMP_HIT_FILE\"):\"[^\"]*\")/Hit_file(getenv(\"G4CMP_HIT_FILE\")?getenv(\"G4CMP_HIT_FILE\"):\"${hitFile}\")/" "/Users/matthewcantor/Desktop/otherG4cmp/work/othersim/othersim/src/RISQTutorialConfigManager.cc"
    sed -i '' "s/Primary_file(\"[^\"]*\")/Primary_file(\"${primaryFile}\")/" "/Users/matthewcantor/Desktop/otherG4cmp/work/singleKIDSim/singleKIDSim/src/RISQTutorialConfigManager.cc"

    sed -i '' "s/double pAbsSilicon = [0-9]*\.[0-9]*/double pAbsSilicon = $probSilicon/" "/Users/matthewcantor/Desktop/otherG4cmp/work/singleKIDSim/singleKIDSim/src/RISQTutorialDetectorConstruction.cc"
    make
    make install
    ./RISQTutorial ../singleKIDSim/G4Macros/pceStudy.mac

    cd ../othersim/AnalysisTools
    root -l <<EOF
.L pceAnalysis.cc
PCEStudy("/Users/matthewcantor/Desktop/otherG4cmp/work/othersim/build/${primaryFile}", "/Users/matthewcantor/Desktop/otherG4cmp/work/singleKIDSim/build/${hitFile}")
.q
EOF
  cd ../../build
  done
done
