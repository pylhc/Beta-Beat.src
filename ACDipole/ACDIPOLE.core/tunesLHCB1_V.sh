sed -e 's/%Q/'0.31'/g' \
    -e 's/%D/'0.329999954376'/g' \
    -e 's/%PATH/'".\/ACDIPOLE.core\/"'/g' \
< /afs/cern.ch/eng/sl/lintrack/Beta-Beat.src//MODEL/LHCB1/AC.Dipole/twiss.LHCB1.vac.mask >  ./ACDIPOLE.core//twiss.LHCB1.vac.madx
