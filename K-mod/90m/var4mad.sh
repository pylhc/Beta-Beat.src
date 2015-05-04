sed -e 's/%IP/'5'/g' \
    -e 's/%MPATH/\/afs\/cern.ch\/eng\/sl\/lintrack\/Beta-Beat.src\/K-mod\/90m\//g' \
    -e 's/%DELTAK/4e-5/g' \
</afs/cern.ch/eng/sl/lintrack/Beta-Beat.src/K-mod/90m//job.statistics.mask > .//t.madx 
