
file="AveAngle.b2.IatS.dat"
echo "#Average and RMS angle of real and imaginary families in degrees for the Injection knobs used at Squeeze: $file"

paste Squeeze/f_cminusimag_b2.IP1.dat Squeeze/f_cminusreal_b2.IP1.dat | awk 'function acos(x){return atan2(sqrt(1-x*x), x)} {t=($3*$9+$4*$10)/sqrt(($3**2+$4**2)*($9**2+$10**2)); t=acos(t);at=at+t;t2=t2+t*t;c=c+1} END{print "IP1 angle ",at/c,sqrt(t2/c-at/c*at/c)}' > $file
paste Squeeze/f_cminusimag_b2.IP2.dat Squeeze/f_cminusreal_b2.IP2.dat | awk 'function acos(x){return atan2(sqrt(1-x*x), x)} {t=($3*$9+$4*$10)/sqrt(($3**2+$4**2)*($9**2+$10**2)); t=acos(t);at=at+t;t2=t2+t*t;c=c+1} END{print "IP2 angle ",at/c,sqrt(t2/c-at/c*at/c)}' >> $file
paste Squeeze/f_cminusimag_b2.IP3.dat Squeeze/f_cminusreal_b2.IP3.dat | awk 'function acos(x){return atan2(sqrt(1-x*x), x)} {t=($3*$9+$4*$10)/sqrt(($3**2+$4**2)*($9**2+$10**2)); t=acos(t);at=at+t;t2=t2+t*t;c=c+1} END{print "IP3 angle ",at/c,sqrt(t2/c-at/c*at/c)}' >> $file
paste Squeeze/f_cminusimag_b2.IP4.dat Squeeze/f_cminusreal_b2.IP4.dat | awk 'function acos(x){return atan2(sqrt(1-x*x), x)} {t=($3*$9+$4*$10)/sqrt(($3**2+$4**2)*($9**2+$10**2)); t=acos(t);at=at+t;t2=t2+t*t;c=c+1} END{print "IP4 angle ",at/c,sqrt(t2/c-at/c*at/c)}' >> $file
paste Squeeze/f_cminusimag_b2.IP5.dat Squeeze/f_cminusreal_b2.IP5.dat | awk 'function acos(x){return atan2(sqrt(1-x*x), x)} {t=($3*$9+$4*$10)/sqrt(($3**2+$4**2)*($9**2+$10**2)); t=acos(t);at=at+t;t2=t2+t*t;c=c+1} END{print "IP5 angle ",at/c,sqrt(t2/c-at/c*at/c)}' >> $file
paste Squeeze/f_cminusimag_b2.IP6.dat Squeeze/f_cminusreal_b2.IP6.dat | awk 'function acos(x){return atan2(sqrt(1-x*x), x)} {t=($3*$9+$4*$10)/sqrt(($3**2+$4**2)*($9**2+$10**2)); t=acos(t);at=at+t;t2=t2+t*t;c=c+1} END{print "IP6 angle ",at/c,sqrt(t2/c-at/c*at/c)}' >> $file
paste Squeeze/f_cminusimag_b2.IP7.dat Squeeze/f_cminusreal_b2.IP7.dat | awk 'function acos(x){return atan2(sqrt(1-x*x), x)} {t=($3*$9+$4*$10)/sqrt(($3**2+$4**2)*($9**2+$10**2)); t=acos(t);at=at+t;t2=t2+t*t;c=c+1} END{print "IP7 angle ",at/c,sqrt(t2/c-at/c*at/c)}' >> $file
paste Squeeze/f_cminusimag_b2.IP8.dat Squeeze/f_cminusreal_b2.IP8.dat | awk 'function acos(x){return atan2(sqrt(1-x*x), x)} {t=($3*$9+$4*$10)/sqrt(($3**2+$4**2)*($9**2+$10**2)); t=acos(t);at=at+t;t2=t2+t*t;c=c+1} END{print "IP8 angle ",at/c,sqrt(t2/c-at/c*at/c)}' >> $file





file="Avef1001.b2.IatS.dat"
echo "#Average and RMS |f_{1001}| for the Injection knobs used at Squeeze: $file"

paste Squeeze/f_cminusimag_b2.IP1.dat Squeeze/f_cminusreal_b2.IP1.dat | awk '{tr=sqrt($3**2+$4**2);ti=sqrt($9**2+$10**2);atr=atr+tr;ati=ati+ti;tr2=tr2+tr*tr;ti2=ti2+ti*ti;c=c+1} END{print "IP1 ",atr/c,sqrt(tr2/c-atr/c*atr/c), ati/c,sqrt(ti2/c-ati/c*ati/c)}' > $file
paste Squeeze/f_cminusimag_b2.IP2.dat Squeeze/f_cminusreal_b2.IP2.dat | awk '{tr=sqrt($3**2+$4**2);ti=sqrt($9**2+$10**2);atr=atr+tr;ati=ati+ti;tr2=tr2+tr*tr;ti2=ti2+ti*ti;c=c+1} END{print "IP2 ",atr/c,sqrt(tr2/c-atr/c*atr/c), ati/c,sqrt(ti2/c-ati/c*ati/c)}' >> $file
paste Squeeze/f_cminusimag_b2.IP3.dat Squeeze/f_cminusreal_b2.IP3.dat | awk '{tr=sqrt($3**2+$4**2);ti=sqrt($9**2+$10**2);atr=atr+tr;ati=ati+ti;tr2=tr2+tr*tr;ti2=ti2+ti*ti;c=c+1} END{print "IP3 ",atr/c,sqrt(tr2/c-atr/c*atr/c), ati/c,sqrt(ti2/c-ati/c*ati/c)}' >> $file
paste Squeeze/f_cminusimag_b2.IP4.dat Squeeze/f_cminusreal_b2.IP4.dat | awk '{tr=sqrt($3**2+$4**2);ti=sqrt($9**2+$10**2);atr=atr+tr;ati=ati+ti;tr2=tr2+tr*tr;ti2=ti2+ti*ti;c=c+1} END{print "IP4 ",atr/c,sqrt(tr2/c-atr/c*atr/c), ati/c,sqrt(ti2/c-ati/c*ati/c)}' >> $file
paste Squeeze/f_cminusimag_b2.IP5.dat Squeeze/f_cminusreal_b2.IP5.dat | awk '{tr=sqrt($3**2+$4**2);ti=sqrt($9**2+$10**2);atr=atr+tr;ati=ati+ti;tr2=tr2+tr*tr;ti2=ti2+ti*ti;c=c+1} END{print "IP5 ",atr/c,sqrt(tr2/c-atr/c*atr/c), ati/c,sqrt(ti2/c-ati/c*ati/c)}' >> $file
paste Squeeze/f_cminusimag_b2.IP6.dat Squeeze/f_cminusreal_b2.IP6.dat | awk '{tr=sqrt($3**2+$4**2);ti=sqrt($9**2+$10**2);atr=atr+tr;ati=ati+ti;tr2=tr2+tr*tr;ti2=ti2+ti*ti;c=c+1} END{print "IP6 ",atr/c,sqrt(tr2/c-atr/c*atr/c), ati/c,sqrt(ti2/c-ati/c*ati/c)}' >> $file
paste Squeeze/f_cminusimag_b2.IP7.dat Squeeze/f_cminusreal_b2.IP7.dat | awk '{tr=sqrt($3**2+$4**2);ti=sqrt($9**2+$10**2);atr=atr+tr;ati=ati+ti;tr2=tr2+tr*tr;ti2=ti2+ti*ti;c=c+1} END{print "IP7 ",atr/c,sqrt(tr2/c-atr/c*atr/c), ati/c,sqrt(ti2/c-ati/c*ati/c)}' >> $file
paste Squeeze/f_cminusimag_b2.IP8.dat Squeeze/f_cminusreal_b2.IP8.dat | awk '{tr=sqrt($3**2+$4**2);ti=sqrt($9**2+$10**2);atr=atr+tr;ati=ati+ti;tr2=tr2+tr*tr;ti2=ti2+ti*ti;c=c+1} END{print "IP8 ",atr/c,sqrt(tr2/c-atr/c*atr/c), ati/c,sqrt(ti2/c-ati/c*ati/c)}' >> $file

file="AveAngle.b1.IatI.dat"
echo "#Beam1 Average and RMS angle of real and imaginary families in degrees for the Injection knobs used at injection: $file"

paste f_cminusimag_b1.IP1.dat f_cminusreal_b1.IP1.dat | awk 'function acos(x){return atan2(sqrt(1-x*x), x)} {t=($3*$9+$4*$10)/sqrt(($3**2+$4**2)*($9**2+$10**2)); t=acos(t);at=at+t;t2=t2+t*t;c=c+1} END{print "IP1 angle ",at/c,sqrt(t2/c-at/c*at/c)}' > $file
paste f_cminusimag_b1.IP2.dat f_cminusreal_b1.IP2.dat | awk 'function acos(x){return atan2(sqrt(1-x*x), x)}{t=($3*$9+$4*$10)/sqrt(($3**2+$4**2)*($9**2+$10**2));  t=acos(t);at=at+t;t2=t2+t*t;c=c+1} END{print "IP2 angle ",at/c,sqrt(t2/c-at/c*at/c)}' >> $file
paste f_cminusimag_b1.IP3.dat f_cminusreal_b1.IP3.dat | awk 'function acos(x){return atan2(sqrt(1-x*x), x)}{t=($3*$9+$4*$10)/sqrt(($3**2+$4**2)*($9**2+$10**2));  t=acos(t);at=at+t;t2=t2+t*t;c=c+1} END{print "IP3 angle ",at/c,sqrt(t2/c-at/c*at/c)}' >> $file
paste f_cminusimag_b1.IP4.dat f_cminusreal_b1.IP4.dat | awk 'function acos(x){return atan2(sqrt(1-x*x), x)}{t=($3*$9+$4*$10)/sqrt(($3**2+$4**2)*($9**2+$10**2));  t=acos(t);at=at+t;t2=t2+t*t;c=c+1} END{print "IP4 angle ",at/c,sqrt(t2/c-at/c*at/c)}' >> $file
paste f_cminusimag_b1.IP5.dat f_cminusreal_b1.IP5.dat | awk 'function acos(x){return atan2(sqrt(1-x*x), x)}{t=($3*$9+$4*$10)/sqrt(($3**2+$4**2)*($9**2+$10**2));  t=acos(t);at=at+t;t2=t2+t*t;c=c+1} END{print "IP5 angle ",at/c,sqrt(t2/c-at/c*at/c)}' >> $file
paste f_cminusimag_b1.IP6.dat f_cminusreal_b1.IP6.dat | awk 'function acos(x){return atan2(sqrt(1-x*x), x)}{t=($3*$9+$4*$10)/sqrt(($3**2+$4**2)*($9**2+$10**2));  t=acos(t);at=at+t;t2=t2+t*t;c=c+1} END{print "IP6 angle ",at/c,sqrt(t2/c-at/c*at/c)}' >> $file
paste f_cminusimag_b1.IP7.dat f_cminusreal_b1.IP7.dat | awk 'function acos(x){return atan2(sqrt(1-x*x), x)}{t=($3*$9+$4*$10)/sqrt(($3**2+$4**2)*($9**2+$10**2));  t=acos(t);at=at+t;t2=t2+t*t;c=c+1} END{print "IP7 angle ",at/c,sqrt(t2/c-at/c*at/c)}' >> $file
paste f_cminusimag_b1.IP8.dat f_cminusreal_b1.IP8.dat | awk 'function acos(x){return atan2(sqrt(1-x*x), x)}{t=($3*$9+$4*$10)/sqrt(($3**2+$4**2)*($9**2+$10**2));  t=acos(t);at=at+t;t2=t2+t*t;c=c+1} END{print "IP8 angle ",at/c,sqrt(t2/c-at/c*at/c)}' >> $file

file="AveAngle.b2.IatI.dat"
echo "#Average and RMS angle of real and imaginary families in degrees for the Injection knobs used at injection: $file"

paste f_cminusimag_b2.IP1.dat f_cminusreal_b2.IP1.dat | awk 'function acos(x){return atan2(sqrt(1-x*x), x)}{t=($3*$9+$4*$10)/sqrt(($3**2+$4**2)*($9**2+$10**2));  t=acos(t);at=at+t;t2=t2+t*t;c=c+1} END{print "IP1 angle ",at/c,sqrt(t2/c-at/c*at/c)}' > $file
paste f_cminusimag_b2.IP2.dat f_cminusreal_b2.IP2.dat | awk 'function acos(x){return atan2(sqrt(1-x*x), x)}{t=($3*$9+$4*$10)/sqrt(($3**2+$4**2)*($9**2+$10**2));  t=acos(t);at=at+t;t2=t2+t*t;c=c+1} END{print "IP2 angle ",at/c,sqrt(t2/c-at/c*at/c)}' >> $file
paste f_cminusimag_b2.IP3.dat f_cminusreal_b2.IP3.dat | awk 'function acos(x){return atan2(sqrt(1-x*x), x)}{t=($3*$9+$4*$10)/sqrt(($3**2+$4**2)*($9**2+$10**2));  t=acos(t);at=at+t;t2=t2+t*t;c=c+1} END{print "IP3 angle ",at/c,sqrt(t2/c-at/c*at/c)}' >> $file
paste f_cminusimag_b2.IP4.dat f_cminusreal_b2.IP4.dat | awk 'function acos(x){return atan2(sqrt(1-x*x), x)}{t=($3*$9+$4*$10)/sqrt(($3**2+$4**2)*($9**2+$10**2));  t=acos(t);at=at+t;t2=t2+t*t;c=c+1} END{print "IP4 angle ",at/c,sqrt(t2/c-at/c*at/c)}' >> $file
paste f_cminusimag_b2.IP5.dat f_cminusreal_b2.IP5.dat | awk 'function acos(x){return atan2(sqrt(1-x*x), x)}{t=($3*$9+$4*$10)/sqrt(($3**2+$4**2)*($9**2+$10**2));  t=acos(t);at=at+t;t2=t2+t*t;c=c+1} END{print "IP5 angle ",at/c,sqrt(t2/c-at/c*at/c)}' >> $file
paste f_cminusimag_b2.IP6.dat f_cminusreal_b2.IP6.dat | awk 'function acos(x){return atan2(sqrt(1-x*x), x)}{t=($3*$9+$4*$10)/sqrt(($3**2+$4**2)*($9**2+$10**2));  t=acos(t);at=at+t;t2=t2+t*t;c=c+1} END{print "IP6 angle ",at/c,sqrt(t2/c-at/c*at/c)}' >> $file
paste f_cminusimag_b2.IP7.dat f_cminusreal_b2.IP7.dat | awk 'function acos(x){return atan2(sqrt(1-x*x), x)}{t=($3*$9+$4*$10)/sqrt(($3**2+$4**2)*($9**2+$10**2));  t=acos(t);at=at+t;t2=t2+t*t;c=c+1} END{print "IP7 angle ",at/c,sqrt(t2/c-at/c*at/c)}' >> $file
paste f_cminusimag_b2.IP8.dat f_cminusreal_b2.IP8.dat | awk 'function acos(x){return atan2(sqrt(1-x*x), x)}{t=($3*$9+$4*$10)/sqrt(($3**2+$4**2)*($9**2+$10**2));  t=acos(t);at=at+t;t2=t2+t*t;c=c+1} END{print "IP8 angle ",at/c,sqrt(t2/c-at/c*at/c)}' >> $file


file="Avef1001.b2.IatI.dat"
echo "#Average and RMS |f_{1001}| for the Injection knobs used at injection: $file"

paste f_cminusimag_b2.IP1.dat f_cminusreal_b2.IP1.dat | awk '{tr=sqrt($3**2+$4**2);ti=sqrt($9**2+$10**2);atr=atr+tr;ati=ati+ti;tr2=tr2+tr*tr;ti2=ti2+ti*ti;c=c+1} END{print "IP1 ",atr/c,sqrt(tr2/c-atr/c*atr/c), ati/c,sqrt(ti2/c-ati/c*ati/c)}' > $file
paste f_cminusimag_b2.IP2.dat f_cminusreal_b2.IP2.dat | awk '{tr=sqrt($3**2+$4**2);ti=sqrt($9**2+$10**2);atr=atr+tr;ati=ati+ti;tr2=tr2+tr*tr;ti2=ti2+ti*ti;c=c+1} END{print "IP2 ",atr/c,sqrt(tr2/c-atr/c*atr/c), ati/c,sqrt(ti2/c-ati/c*ati/c)}' >> $file
paste f_cminusimag_b2.IP3.dat f_cminusreal_b2.IP3.dat | awk '{tr=sqrt($3**2+$4**2);ti=sqrt($9**2+$10**2);atr=atr+tr;ati=ati+ti;tr2=tr2+tr*tr;ti2=ti2+ti*ti;c=c+1} END{print "IP3 ",atr/c,sqrt(tr2/c-atr/c*atr/c), ati/c,sqrt(ti2/c-ati/c*ati/c)}' >> $file
paste f_cminusimag_b2.IP4.dat f_cminusreal_b2.IP4.dat | awk '{tr=sqrt($3**2+$4**2);ti=sqrt($9**2+$10**2);atr=atr+tr;ati=ati+ti;tr2=tr2+tr*tr;ti2=ti2+ti*ti;c=c+1} END{print "IP4 ",atr/c,sqrt(tr2/c-atr/c*atr/c), ati/c,sqrt(ti2/c-ati/c*ati/c)}' >> $file
paste f_cminusimag_b2.IP5.dat f_cminusreal_b2.IP5.dat | awk '{tr=sqrt($3**2+$4**2);ti=sqrt($9**2+$10**2);atr=atr+tr;ati=ati+ti;tr2=tr2+tr*tr;ti2=ti2+ti*ti;c=c+1} END{print "IP5 ",atr/c,sqrt(tr2/c-atr/c*atr/c), ati/c,sqrt(ti2/c-ati/c*ati/c)}' >> $file
paste f_cminusimag_b2.IP6.dat f_cminusreal_b2.IP6.dat | awk '{tr=sqrt($3**2+$4**2);ti=sqrt($9**2+$10**2);atr=atr+tr;ati=ati+ti;tr2=tr2+tr*tr;ti2=ti2+ti*ti;c=c+1} END{print "IP6 ",atr/c,sqrt(tr2/c-atr/c*atr/c), ati/c,sqrt(ti2/c-ati/c*ati/c)}' >> $file
paste f_cminusimag_b2.IP7.dat f_cminusreal_b2.IP7.dat | awk '{tr=sqrt($3**2+$4**2);ti=sqrt($9**2+$10**2);atr=atr+tr;ati=ati+ti;tr2=tr2+tr*tr;ti2=ti2+ti*ti;c=c+1} END{print "IP7 ",atr/c,sqrt(tr2/c-atr/c*atr/c), ati/c,sqrt(ti2/c-ati/c*ati/c)}' >> $file
paste f_cminusimag_b2.IP8.dat f_cminusreal_b2.IP8.dat | awk '{tr=sqrt($3**2+$4**2);ti=sqrt($9**2+$10**2);atr=atr+tr;ati=ati+ti;tr2=tr2+tr*tr;ti2=ti2+ti*ti;c=c+1} END{print "IP8 ",atr/c,sqrt(tr2/c-atr/c*atr/c), ati/c,sqrt(ti2/c-ati/c*ati/c)}' >> $file
