# set up intel compiler 12.0, inspector and amplifier
#. $DIR/x86_64/Compiler/11.1/056/bin/iccvars.sh intel64
#. $DIR/x86_64/itt/3.1/tcheck/bin/32e/tcheckvars.sh
#. $DIR/x86_64/itt/3.1/tprofile/bin/32e/tprofilevars.sh

INTEL_LOCAL_DIR=/afs/cern.ch/sw/IntelSoftware/linux
source $INTEL_LOCAL_DIR/setup.sh
source $INTEL_LOCAL_DIR/x86_64/xe2013/bin/iccvars.sh intel64
source $INTEL_LOCAL_DIR/x86_64/xe2013/bin/ifortvars.sh intel64

export PATH="/afs/cern.ch/sw/IntelSoftware/linux/x86_64/xe2013/inspector_xe_2013/bin64/:$PATH"
export PATH="/afs/cern.ch/sw/IntelSoftware/linux/x86_64/xe2013/vtune_amplifier_xe_2013/bin64/:$PATH"

icc -V
ifort -V
inspxe-cl -V
amplxe-cl -V
