#!/usr/bin/ksh

x=0.0000   # BPM noise in m 
y=0.0000   # BPM noise in m
turn=1000  # Starting turn of the AC dipole flat plateau
#Nfail=0   # Number of BPMs failing, not active

echo Parsing MADX file
echo and Adding Random errors to BPM data: $x $y
echo and number of BPMs failing: $Nfail

#  trackone syntax more or less
#  $2 turn number
#  $9 s location
#  when #segment $6 is bpm name, $3 is the total number of turns
#
awk 'function gauss(){pi=3.1415926;
                      x1=rand();x2=rand();
                      y1 = sqrt(-2*log(x1))*cos(2*pi*x2);
                      return y1;}
     BEGIN{k=0;myend=0;xerr='$x';yerr='$y';srand(); turn='$turn'; if (turn<1) {turn=1};}
     /segment/{name=$6; tt=$3} 
     /segment/ && $6~/end/ {name=$6;myend=1;} 
     $0!~/segment/ && $0!~/^@/ && $0!~/^\*/ && $0!~/%/ {nameb[$9":"name]=name; 
              if ($2==1) {ind[k]=$9":"name; k++;}
              s[$9":"name]=$9;
              x[$9":"name":"$2]=$3;px[$9":"name":"$2]=$4; 
              y[$9":"name":"$2]=$5; py[$9":"name":"$2]=$6} 
     END{   print "# title";
            for (kk=0;kk<k;kk++){ pos=ind[kk];
             if (match(nameb[pos],"BPM") || match(nameb[pos],"bpm") ) {printf("0 %s %e " ,nameb[pos], pos); 
            for (i=turn;i<tt;i++) printf("%e ", x[pos":"i]+gauss()*xerr);printf("\n");}
             if (match(nameb[pos],"BPM") || match(nameb[pos],"bpm")) {printf("1 %s %e ",nameb[pos], pos); 
            for (i=turn;i<tt;i++) printf("%e ", y[pos":"i]+gauss()*yerr);printf("\n");}}
            }' trackone > ALLBPMs


#awk 'BEGIN{srand();
#           Nfail='$Nfail';
#           for (i=0;i<Nfail;i++) {a[i]=int(rand()*1000)+2;}}
#	  {isin=0;
#           for (i=0;i<Nfail;i++) if (a[i]==NR) isin=1;  }    
#     isin==0{print $0}     ' tt > ./data/ALLBPMs
    

# rm tt


rm trackone
