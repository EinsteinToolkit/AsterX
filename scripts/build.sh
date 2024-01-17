#!/bin/bash

set -ex

export ASTERXSPACE="$PWD"
export WORKSPACE="$PWD/../workspace"
mkdir -p "$WORKSPACE"
cd "$WORKSPACE"

cd Cactus

# Set up SimFactory
cp "$ASTERXSPACE/scripts/actions-$ACCELERATOR-$REAL_PRECISION.cfg" ./simfactory/mdb/optionlists
cp "$ASTERXSPACE/scripts/actions-$ACCELERATOR-$REAL_PRECISION.ini" ./simfactory/mdb/machines
cp "$ASTERXSPACE/scripts/actions-$ACCELERATOR-$REAL_PRECISION.run" ./simfactory/mdb/runscripts
cp "$ASTERXSPACE/scripts/actions-$ACCELERATOR-$REAL_PRECISION.sub" ./simfactory/mdb/submitscripts
cp "$ASTERXSPACE/scripts/asterx.th" .
cp "$ASTERXSPACE/scripts/defs.local.ini" ./simfactory/etc

# For Formaline
git config --global user.email "asterx@einsteintoolkit.org"
git config --global user.name "Github Actions"

# Build
# The build log needs to be stored for later.
time ./simfactory/bin/sim --machine="actions-$ACCELERATOR-$REAL_PRECISION" build -j $(nproc) sim 2>&1 |
    tee build.log

# Check whether the executable exists and is executable
test -x exe/cactus_sim
