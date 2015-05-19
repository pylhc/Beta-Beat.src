sed -e 's/%IP/'5'/g' \
    -e 's/%MPATH/\/afs\/cern.ch\/eng\/sl\/lintrack\/Beta-Beat.src\/K-mod\//g' \
    -e 's/%DELTAK/1e-5/g' \
</afs/cern.ch/eng/sl/lintrack/Beta-Beat.src/K-mod/TESTING//job.statistics.mask > .//t.madx 
