
if [ -f dumpXB$1.ps ]; then
  #echo "dumpXB$1.ps exists"
  mv dumpXB$1.ps dumpXB$1tmp.ps
  psmerge -odumpXB$1.ps dumpXB$1tmp.ps XB$1.ps
  mv dumpYB$1.ps dumpYB$1tmp.ps
  psmerge -odumpYB$1.ps dumpYB$1tmp.ps YB$1.ps  
  rm -f dumpXB$1tmp.ps dumpYB$1tmp.ps
else
  #echo "dumpXB$1.ps does not exist"
  mv XB$1.ps dumpXB$1.ps 
  mv YB$1.ps dumpYB$1.ps 
fi
