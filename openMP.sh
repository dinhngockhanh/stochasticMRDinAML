#!/bin/sh
source /opt/asn/etc/asn-bash-profiles-special/modules.sh
module load intel
module load acml
make all
./program_stochastic
