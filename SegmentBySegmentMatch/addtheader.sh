npar=$#
echo NPARS = $npar

ipno=1;

if [ $npar -lt 1 ] ; then
 ipno=1;
else
 ipno=$1; 
 echo read the value $ipno
fi


if [ $ipno -lt 1 ] ; then
 ipno=1;
fi

if [ $ipno -gt 8 ] ; then
 ipno=1;
fi


if [ $? -ne 0 ] ; then
 echo "Not a number? using 1"
 ipno=1;
fi


echo Using IP $ipno


sed -i 1i'@ TYPE %05s  \"USER\"' Beam1/sbs/sbsphasext_IP${ipno}.out
sed -i 1i'@ TYPE %05s  \"USER\"' Beam1/sbs/sbsphaseyt_IP${ipno}.out
sed -i 1i'@ TYPE %05s  \"USER\"' Beam2/sbs/sbsphasext_IP${ipno}.out
sed -i 1i'@ TYPE %05s  \"USER\"' Beam2/sbs/sbsphaseyt_IP${ipno}.out
