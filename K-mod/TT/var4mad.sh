sed -e 's/%START/'0'/g' \
    -e 's/%ACCEL/'LHCB1'/g' \
    -e 's/%IR/'5'/g' \
    -e 's/%B12/'B1'/g' \
    -e 's/%Q12/'Q1'/g' \
    -e 's/%deltar/'0.0009'/g' \
    -e 's/%deltal/'0.00105'/g' \
    -e 's/%deltak/'1e-5'/g' \
    -e 's/%XY/'x'/g' \
    -e 's/%MODIFIERS/\/afs\/cern.ch\/eng\/sl\/lintrack\/Beta-Beat.src\/K-mod\/modifiers.madx/g' \
</afs/cern.ch/eng/sl/lintrack/Beta-Beat.src/K-mod//job.Qmatch.mask > .//t.madx 
