#!/usr/bin/env bash
# Build gdmlcheck/gdmlprobe against LCG_110a Geant4 11.4.2 (same recipe as the
# Celeritas geocheck build).
set -euo pipefail
cd "$(dirname "$0")"
IMAGE=/cvmfs/singularity.opensciencegrid.org/opensciencegrid/osgvo-el9:latest
VIEW=/cvmfs/sft.cern.ch/lcg/views/LCG_110a/x86_64-el9-gcc13-opt
XERCES=$VIEW/lib/libxerces-c-3.3.so
apptainer exec -B /cvmfs -B /scratch "$IMAGE" bash -c "
  source $VIEW/setup.sh > /dev/null 2>&1
  set -e
  LIBS=\"\$(geant4-config --libs-without-gui | sed 's/ -lG4interfaces//') $XERCES\"
  for p in gdmlcheck gdmlprobe gdmlinside; do
    g++ -O2 \$(geant4-config --cflags) \$p.cc -o \$p \$LIBS
  done
"
