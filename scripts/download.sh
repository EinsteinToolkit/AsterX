#!/bin/bash

set -ex

export ASTERXSPACE="$PWD"
export WORKSPACE="$PWD/../workspace"
mkdir -p "$WORKSPACE"
cd "$WORKSPACE"

# Check out Cactus
wget https://raw.githubusercontent.com/gridaphobe/CRL/master/GetComponents
chmod a+x GetComponents
./GetComponents --no-parallel --shallow "$ASTERXSPACE/scripts/asterx.th"

cd Cactus

# Create a link to the AsterX repository
ln -s "$ASTERXSPACE" repos
# Create links for the AsterX thorns
mkdir -p arrangements/AsterX
pushd arrangements/AsterX
ln -s ../../repos/AsterX/* .
popd
