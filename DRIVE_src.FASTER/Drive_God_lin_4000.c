#include <sys/types.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#define MAXPICK 1100
#define MAXTURNS 4000   /* Critial parameter in sussixfordiveNoO.f !!!!!!!!!!!!!!*/
#define MAXTURNS4 16000 /*Always four times  MAXTURNS*/
#define MAXRUNS 100
#define M_PI 3.14159265358979310862446895044
char *get_line(FILE *);
char *get_word(FILE *);
void  get_name(),  sussix_inp(int, char*), betabeat(),savefort300(int);
int BPMstatus(int, int), norm(int,int);
/* isin(char *, char *);*/
double  getparameter();
double lector(int, char, char *);


char namefile[2000], datafile[200], labelfile[200], noisefile[200],
 datafileaux[200],inputfile[2000], ss[1000], labelfilepath[200], *bpmname[MAXPICK],
 originalpath[600];
/*     *labelfilepath="/afs/cern.ch/project/lhcnap/nap2/user/rtomas/BPM_noise/Noise_files/";*/ /* REAL ONE!!!!!!!!!!!!!!! RAPID RUN !!!!!!!!! */

double matrix[MAXPICK][MAXTURNS], bpmpos[MAXPICK], tunex, tuney, tune[2], 
      amplitude[14], phase[14], parameter[MAXRUNS], istun, noise1, noiseAve,
         maxpeak, maxfreq, maxmin, co, co2, windowa1, windowa2, windowb1, windowb2, allfreqsx[300],
  allampsx[300], allfreqsy[300], allampsy[300], doubleToSend[MAXTURNS4+4],
  allbpmamp[MAXPICK], allbpmphase[MAXPICK];

int turns, normalisation, iparameter=0, scum=0, in=0, iir=0, labelrun=0,
  label[MAXPICK],NoVerticalData=0,hv[MAXPICK], hvt[MAXPICK],nslines ;

main(int argc, char **argv)
{
	 double amin=1e17, amax, amp, backtune[2],backphase[3],
		  backamp[14],offset,offset1,offrms,tunesum[2],tune2sum[2],zsum,
		  ysum,zysum,z2sum,amp10bpm[MAXPICK],amp30bpm[MAXPICK], kper, 
		  amp01bpm[MAXPICK],amp30ave[MAXRUNS], amp01ave[MAXRUNS], 
		  amp00ave[MAXRUNS], amp1sum[2],amp1sum2[2],amp2sum[2], 
		  amp2sum2[2], phase01ave[MAXRUNS], phase02sum[2],
		  amp00bpm[MAXPICK], amp00sum[2], amp00sum2[2], raizbeta[262], 
		  amp4sum[2], amp4sum2[2],amp5sum[2],amp6sum[2], amp40bpm[MAXPICK], 
		  amp50bpm[MAXPICK],amp60bpm[MAXPICK], amp40ave[MAXRUNS],
		  amp50ave[MAXRUNS],amp5sum2,amp60ave[MAXRUNS],amp20bpm[MAXPICK], 
		  amp20sum[2], amp20sum2[2],amp20ave[MAXRUNS], meantune, fact, 
                  upperpeak, lowerpeak, peaktopeak[MAXPICK], averamp[2], averamp2[2];
      
	 int i=0, kick, tot,  ij, finii, finij, mini, counth, countv,maxcounthv,
		  pickupaxis1[MAXPICK][2],pickupaxis2[MAXPICK][2],start, flag=0, 
		  ii=0, jj=0, j, jjx, jjy,count0[2], count1[2], count2[2], count8[2], 
		  count9[2], pickstart, pickend, kcase, ibetabeat, opticlabel=1, 
		  
lastslash, kjl, count4[2], count5[2],count6[2], iim, iformat,
		  ipick, ipickm,  backlai,  backlaij, kk, Nbpms, Nturns, avercount[2],
                  counthl=0, countvl=0;

	 char s[300],bpmfile[300], linfile[300],linfilex[2000],linfiley[2000], *opticfile="monitors-hpi40.txt",*pointer,*ps, cmd[6000], 
	   path[4000], bsfile[400];

	 FILE *df, *fc, *di, *ft, *fs, *fb, *optic, *lf, *nf, *foc, *fdec, 
	   *fdod, *kf, *linfx, *linfy, *fBS;

 /* Memory allocation*/
	 for (i=0; i<MAXPICK;i++){
	   bpmname[i] = (char *) calloc(50, sizeof(char));
	 }	 
	 	
	 
/*  Path to DrivingTerms and Drive.inp */ 
printf("\n Path to Drive input: %s\n", argv[1]);
strcpy(path, argv[1]);
 strcat(path, "/sussix_v4.inp");
strcpy(namefile, argv[1]); 
strcat(namefile,"/DrivingTerms");
strcpy(inputfile, argv[1]);
strcat(inputfile, "/Drive.inp");
printf("\n Drive.inp: %s\n", inputfile);




/* check the file namefile */
         if ( access(namefile, 0) < 0 ) {
                  /* doesn't exist --> stop*/
                  printf("\nNo file %s for reading the name of the Data file\n",namefile);
                  return;
         }




/* check the input file Drive.inp */  
	 if ( access(inputfile, 0) < 0 ) {
		  /* doesn't exist --> stop*/
		  printf("\nNo input file %s\n", inputfile);
	   return;
	 }


	 count0[0]=0;
count0[1]=0;
	       tunesum[0] = 0;
	       tune2sum[0] = 0;
        tunesum[1] = 0;
	       tune2sum[1] = 0;
 /* check the OPTIC  file and take Beta(s) in the monitors if it exists */  
	 if ( access(opticfile, 0) < 0 ) {
		  /* doesn't exist --> forget*/
	   printf("\nNO OPTIC FILE: %s\n", opticfile);
	   opticlabel=0;
	 
	 } else {
	   optic=fopen(opticfile,"r");
	   for(i=0; i<262; i++){
			pointer=get_word(optic);
			for(kjl=0; *pointer != '\0'; kjl++){	 
				 s[kjl]=*pointer;
				 pointer++;
			}
			s[kjl]='\0';
			
			/*printf("%e ",atof(s));*/
		 raizbeta[i]=sqrt(atof(s));
		 get_word(optic);get_word(optic);get_word(optic);
	     
	   }
	   printf("\nOPTIC FILE: %s USED\n", opticfile);
	   fclose(optic);
	 }

	 /* Opening files for the output */	 
	 /*fs=fopen("sextupoles","w");
	 ft=fopen("detuning","w");
	 fc=fopen("coupling", "w");
	 
	 foc=fopen("octupoles", "w");
	 fdec=fopen("decapoles", "w");
	 fdod=fopen("dodecapoles", "w");*/
	 
		  
	 i=0;
	 di=fopen(inputfile, "r"); 
	 while (i<=17){  
		  *s='\0';   
		  if ((s[0]=getc(di))== '=') {
		    i++;
		    for (jj=0; ((s[jj]=getc(di)) != EOF) & 
			   (s[jj]  !='\n') & (s[jj]  !=' ')  ; jj++);
		    s[jj]='\0';
		    switch (i){
		    case 1:
		      kick=atof(s)-1;  /* C arrays start at 0 */
		      break;
		    case 2:
		      kcase=atof(s);
		      break;
		    case 3:   
		      kper=atof(s);
		      break;
		    case 4:   
		      tunex=atof(s);
		      break;
		    case 5:   
		      tuney=atof(s);
		      break;
/*		    case 6:   
		      turnspercent=atof(s);   
		      break;*/
		    case 6:   
		      pickstart=atoi(s);
		      break;
		    case 7:		
		      pickend=atoi(s);
		      break;
		    case 8:
		      normalisation=atoi(s);
		      break;
		    case 9:
		      istun=atof(s);
		      break;
		    case 10:
		      ibetabeat=atoi(s);
		      break;
		    case 11:
		      iir=atoi(s);
		      break;
		    case 12:
		      labelrun=atoi(s);
		      break;
		    case 13:
		      iformat=atoi(s);
		      break;
		    case 14:
		      strcpy(labelfilepath,s);
		      break;
		    case 15:
		      windowa1=atof(s);
		      break;
		    case 16:
		      windowa2=atof(s);
		      break;
		    case 17:
		      windowb1=atof(s);
		      break;
		    case 18:
		      windowb2=atof(s);
		      break;
		    }
		  }
	 }
      	 fclose(di);
	 if (kick>=0) printf("Known kick in %d turn\n",kick+1);
	 if (kcase==1) printf("Horizontal case\n");
	 else if (kcase==0)  printf("Vertical case\n");
	 else {
		  printf("No proper case in Drive.inp\n");
		  return;
	 }

	 if (labelrun==1) 
	   printf("\n LABELRUN: NOISE FILES WILL BE WRITTEN TO NOISEPATH\n");
	 printf("Normalisation: %d\n", normalisation);
	 printf("pickstart: %d, pickend: %d\n", pickstart, pickend);


	 while(scum==0){
	 /* From namefile assign datafile, increase iparameter,
	  asign parameter[iparameter] and turns.               */ 
	 get_name();
	 
/* Check the file datafile */
	 if ( access(datafile, 0) < 0 ) {
		  /* doesn't exist --> stop*/
		  printf("\nNo Data file %s\n", datafile);
		  break;
	 }
	 printf("Data File: %s\n", datafile);
	 /*constructing name of BPM files and labelfile*/  
	 lastslash=0;
	 for (jj=0; datafile[jj] !='\0'; jj++){
		  if (datafile[jj]=='/') {lastslash=jj+1;}
	 }
	 for (jj=0; datafile[jj+lastslash] !='\0'; jj++) bpmfile[jj]=datafile[jj+lastslash];
	 bpmfile[jj]='\0'; 
	 for (jj=0; labelfilepath[jj] !='\0'; jj++){
		  if (labelfilepath[jj]=='/') {lastslash=jj+1;}
	 }/* If labelfilepath does not end with '/' put it */
	 if (lastslash!=jj) strcat(labelfilepath , "/");
	 strcpy(labelfile, labelfilepath);
	 strcat(labelfile, bpmfile);
	 strcpy(noisefile, labelfile);
	 strcat(labelfile, "_label");

	strcpy(linfilex, argv[1]);
	strcpy(linfiley, argv[1]);
	strcat(linfilex, "/");
        strcat(linfiley, "/");
	 strcat(linfilex, bpmfile);
	 strcat(linfiley, bpmfile);
	 strcat(bpmfile, "_bpm");
	 strcat(noisefile, "_noise");
                  strcat(linfilex, "_linx");
                  strcat(linfiley, "_liny");


	 /* For Rapid run !!! ************************* 29 Jul 2002 */
	 /*strcpy(labelfile, "labelfile");*/

	 /*fb=fopen(bpmfile, "w");*/
         linfx=fopen(linfilex, "w");
         linfy=fopen(linfiley, "w");
         fprintf(linfx,"* NAME  S   BINDEX SLABEL  TUNEX   MUX  AMPX   NOISE   PK2PK  AMP01 PHASE01 CO CORMS AMP_20  PHASE_20\n ");
         fprintf(linfx,"$ %%s  %%le   %%le    %%le  %%le  %%le  %%le  %%le  %%le  %%le  %%le  %%le  %%le    %%le     %%le");
         fprintf(linfy,"* NAME  S   BINDEX SLABEL  TUNEY  MUY  AMPY   NOISE   PK2PK AMP10 PHASE10 CO CORMS\n");
         fprintf(linfy,"$ %%s  %%le  %%le    %%le  %%le  %%le  %%le  %%le  %%le %%le  %%le %%le  %%le");
	 if (labelrun==1) {lf=fopen(labelfile, "w"); nf=fopen(noisefile, "w");}

/* Constructing a matrix with all the data from the pick-ups for iformat i=0,1*/ 		 
	 if (iformat<2){
	 jj=0;
	 i=0;
	 df=fopen(datafile, "r");
	 while((s[0]=getc(df)) !='\n'){}
	 s[0]=getc(df);
	 while(s[0] != EOF ){
		  if (s[0]=='\n') {jj++; ii=0;}
		  if ((s[0]=='\n' || s[0]==' ' || s[0]=='\t') && flag==1) flag=0;
		  if ((s[0] !='\n' && s[0] !=' ' && s[0] !='\t') && flag==0){
			   while(!isspace(s[i]) && s[i] !=EOF){
					i++;
					s[i]=getc(df);
			   }
			   s[i+1]=s[i];
			   s[i]='\0';
			   iim=ii/2; 
			   			   
			   if (opticlabel==0) fact=1.0;
			   else fact=raizbeta[ii];
					
			   if (iformat==0) {			
					if (iim == (ii/2.0)) 
						 matrix[iim+MAXPICK/2][jj]=atof(s)/fact;
					else matrix[iim][jj]=atof(s)/fact;	
			   }
			   else {
					matrix[ii][jj]=atof(s)/fact;
			   }
			   ii++;
			   flag=1;
			   s[0]=s[i+1];
			   i=0;
		  }
		  if (flag==0) s[0]=getc(df);
	 }
	 
	 /*turns=turnspercent*jj/100;  !!!!!!!!!!*//* C count */
	 fclose(df);  
	 
/* Label spare pick-ups from datafile & labelfile*/
	 for(i=0; i<MAXPICK; i++) label[i]=0;
	 df=fopen(datafile, "r");
	 tot=0;
	 for(i=0; i<MAXPICK; i++){
		  while(isspace(s[0]=getc(df))){}
		  if (/*s[0]=='N' || s[0]=='S' || s[0]=='W' || s[0]=='O' ||*/ s[0]=='B'){
			 iim=i/2;
			   if (iformat==0)
					if (iim == (i/2.0)) label[iim+MAXPICK/2]=1;
					else { label[iim]=1; 
					}
			   if (iformat==1) label[i]=1;
		  }
		  while(!isspace(s[0]=getc(df))){}
		  /*Provisional to ensure 'Flats' recognisition!!!!!!!!!!*/
		  /*if (matrix[i][1]==0.0 && matrix[i][2]==0.0) label[i]=0;*/
	 }
	 /*if (labelrun==0){
		   check the file labelfile 
		  
		  if ( access(labelfile, 0) < 0 ) {
		
			   printf("\nNo file %s for labeling \n", labelfile);
		  } else {
			   printf("File %s being used \n", labelfile);
			   lf=fopen(labelfile, "r");
			   while(1){
					i=1;
					while(isspace(s[0]=getc(lf))){}
					if(s[0]==EOF) break;
					while(!isspace(s[i]=getc(lf))){if(s[i]==EOF) break;i++;}
					s[i+1]=s[i];
					s[i]='\0';
					ipick=atoi(s)-1;
					ipickm=ipick/2;
					if (ipickm == (ipick/2.0)) 
						 {label[ipickm+MAXPICK/2]=0;
						 printf(" %d ", ipickm+MAXPICK/2);
						 }
					
					else {label[ipickm]=0;printf(" %d ", ipickm);}
					  
					if(s[0]==EOF) break;
			   }
			   fclose(lf);
		  } 
	 
	 }*/
	 

	 fclose(df);

	
     }/* Closes if iformat<2 */

/* Constructing a matrix with all the data from the pick-ups for iformat i=2, RHIC*/ 
	 flag=0;
         ii=0;
	 for(i=0; i<MAXPICK; i++) label[i]=0;
	 if (iformat==2) { 
	 jj=0;
	 jjx=-1;
	 jjy=MAXPICK/2-1;
	 i=0;
	 df=fopen(datafile, "r");
	 while((originalpath[i]=getc(df)) !='\n') i++;
	 originalpath[i]='\0'; printf("original path: %s\n", originalpath);
	 s[0]=getc(df);
	 i=0;
	 /*	 while(s[0] != EOF && (jjx<=pickend || jjy<=pickend)){*/
   	 while(s[0] != EOF){
		  if (s[0]=='\n') {jj++; ii=0;}
		  if ((s[0]=='\n' || s[0]==' ' || s[0]=='\t') && flag==1) flag=0;
		  if ((s[0] !='\n' && s[0] !=' ' && s[0] !='\t') && flag==0){
		    /*s[1]='\0';printf("Jo i=%d jj=%d ii=%d s=%s\n", i, jj, ii, s);*/
			   while(!isspace(s[i]) && s[i] !=EOF){
					i++;
					s[i]=getc(df);
					if(i>100) {
					  s[i+1]='\0';
					  printf("i bigger than 100 i=%d jj=%d ii=%d s=%s\n", i, jj, ii, s); exit(0);}
			   }
			   s[i+1]=s[i];
			   s[i]='\0';
			   			   
			   if (opticlabel==0) fact=1.0;
			   else fact=raizbeta[ii];
			   if(ii>=MAXTURNS) { printf("ii bigger than MAXTURNS, %d %d", jj, ii); exit(0);}
			   if(jj>=MAXPICK){ printf("jj bigger than MAXPICK, %d %d",jj, ii); exit(0);}
			   if(ii==0){ hv[jj]=atoi(s); if(hv[jj]==0) jjx++; else jjy++;}
			   
			   if(ii==1){ if(hv[jj]==0) {hvt[jjx]=0;strcpy(bpmname[jjx],s);label[jjx]=1;}  
			   else {hvt[jjy]=1;strcpy(bpmname[jjy],s);label[jjy]=1;}}
			   
			   if(ii==2){ if(hv[jj]==0)  bpmpos[jjx]=atof(s);
			   else   bpmpos[jjy]=atof(s);}

			   if(ii>2) { 
			     if(hv[jj]==0) { 
			       matrix[jjx][ii-3]=atof(s)/fact;
			       /*if ((ii-3)>kick && matrix[jjx][ii-3]==99) label[jjx]=0;*/
			     } else {
			       matrix[jjy][ii-3]=atof(s)/fact;
			       /*if ((ii-3)>kick && matrix[jjy][ii-3]==99) label[jjy]=0;*/
			     }
			   }
			   ii++;
			   flag=1;
			   s[0]=s[i+1];
			   i=0;
		  }
		  if (flag==0) s[0]=getc(df);
	 }
	 fclose(df);  
	 Nbpms=jj;	 
	 Nturns=ii-3;
	 counth=jjx+1;
	 countv=jjy+1;
         /*for(ii=0;ii<counth;ii++) if (label[ii]==1) counthl++;
	   for(ii=MAXPICK/2;ii<countv;ii++) if (label[ii]==1) countvl++;*/
	 if(0){ /* to produce file with pick-up names and number*/
	   df=fopen("bpmnames","w");
	   for(ii=0;ii<counth;ii++)
	     fprintf(df, "%d %s %e\n", ii, bpmname[ii], bpmpos[ii]);
	   for(ii=MAXPICK/2;ii<countv;ii++)
	     fprintf(df, "%d %s %e\n", ii, bpmname[ii], bpmpos[ii]);
	   fclose(df);
	 }

	 

	 /* SOme statistics and checks*/
	 
	 printf("Total number of pick-ups: %d Last turn number: %d\n", Nbpms, Nturns);
	 

	 printf("Horizontal pick-ups: %d   Vertical pick-ups: %d\n ",counth, -MAXPICK/2+countv);

	 printf("name of pick-up 0: %s pos: %f last turn:%e\n", bpmname[0], bpmpos[0], matrix[0][999]);
	 printf("name of pick-up 165: %s pos: %f first turn:%e\n", bpmname[165], bpmpos[165], matrix[165][0]);
	 

	 }          /* CLoses if iformat==2*/


/* searching for two working adjacent pick-ups */ 
/* after the Q-kickers for finding the kick*/
	 if (kick<0){
	   i=0;
	   start=-(kcase-1)*MAXPICK/2+2;
	   while(label[start]==0 || label[start+2]==0){
	     start=start+2;
	   } 
	   
	   printf("looking for kick in pick-up:%d\n", start+1);
	   offset=0;
	   /* Find kick here and get kick */	
	   jj=1;
	   while(kick<0 && jj < turns){
	     if (abs(matrix[start][jj]- matrix[start][jj-1])>kper) {kick=jj; break;}
	     jj++;
	   }
	   if (kick<0){
	     printf("NO KICK FOUND\n");
	     return; 	
	   }
	   else
	     printf("Found kick in turn:%d\n", kick+1); /*Natural count*/
	 }
	 
	 
	 /*Calculating, substracting and removing the offset*/
	 for(i=0; i< MAXPICK; i++){
	   if (label[i] == 1) {
	     offset=0.0;
	     /*for(j=0; j< kick; j++){  FOR NORMAL RUN 
 	       offset=offset+matrix[i][j]/(kick); */
	     /*for(j=0; j< turns; j++){    */ /*FOR ANORMAL RUN */ /*
	      offset=offset+matrix[i][j]/turns; */ 
	     
	     for(j=kick; j< turns; j++)
	       matrix[i][j-kick]=matrix[i][j]-offset;
	   }
	 }
	 turns -= kick;
	 printf("Turns to be procesed: %d\n", turns);


	 /* Measuring the peak-to-peak amplitude in the first 10 turns after kick*/
         averamp[0]=0.0;
	 averamp[1]=0.0;
         averamp2[0]=0.0;
	 averamp2[1]=0.0;
         avercount[0]=0;
         avercount[1]=0;
	 for(i=0; i< MAXPICK; i++){
	   upperpeak=-10000;
	   lowerpeak=10000;
	   if (label[i] == 1) {/* This is not really needed but ok...*/
	     for(j=0; j< 10 ; j++){
               if (matrix[i][j]==99) printf("99 found in %d %d %s\n", i, j, bpmname[i]);
	       if (matrix[i][j]>upperpeak) upperpeak=matrix[i][j];
	       if (matrix[i][j]<lowerpeak) lowerpeak=matrix[i][j];
	     }
	     peaktopeak[i]=upperpeak-lowerpeak;
             /*Compute average amplitude in the arcs*/ 
	    /* if (((bpmpos[i]>154 && bpmpos[i]<482) || (bpmpos[i]>797 && bpmpos[i]<1124) || (bpmpos[i]>1432 && bpmpos[i]<1760) || (bpmpos[i]>2074 && bpmpos[i]<2400) || (bpmpos[i]>2710 && bpmpos[i]<3036) || (bpmpos[i]>3352 && bpmpos[i]<3678) )) {averamp[hvt[i]]+=peaktopeak[i];averamp2[hvt[i]]+=peaktopeak[i]*peaktopeak[i];avercount[hvt[i]]++;*/ 
	     /*printf("%d %s %e %e\n", hvt[i],bpmname[i], bpmpos[i],peaktopeak[i]);
	     }*/
	   }
	 }
	 
	 /*printf("\n counts:%d %d\n",avercount[0],avercount[1] );

	 if (avercount[0]>2) {
	   averamp[0]=averamp[0]/avercount[0]/2;
           averamp2[0]=sqrt(averamp2[0]/avercount[0]/4.-averamp[0]*averamp[0]);
	 } else {
	   printf("\nNO Horizontal pickups in  arcs!!!\n");
	   averamp[0]=0.0;
           averamp2[0]=0.0;
	 }

	 if (avercount[1]>2) {
	   averamp[1]=averamp[1]/avercount[1]/2;
           averamp2[1]=sqrt(averamp2[1]/avercount[1]/4.-averamp[1]*averamp[1]);
	 } else {
	   printf("\nNO Vertical pickups  arcs!!!\n");
	   averamp[1]=0.0;
	   averamp2[1]=0.0;
	   }*/



	 

	 /* First part of the analysis: Determine  phase of all pick-ups and noise*/
	 sussix_inp(1, path);
         if (counth>=(countv-MAXPICK/2)) maxcounthv=counth; 
	 else maxcounthv=-MAXPICK/2+countv;
         
	 if (maxcounthv>=pickend) maxcounthv=pickend+1;
	 if (maxcounthv>=MAXPICK) { printf("\nNot enough Pick-up memory\n");exit(0);}

         printf("BPMs in loop:%d\n", maxcounthv);
	 for (i=pickstart;  i<=maxcounthv ; i++){
           ii=i; 
	   ij=i+MAXPICK/2;
           
	   if (ij >= countv) ij = countv-1;
           if (ii >= counth) ii = counth-1;
	   printf("BPM indexes (H,V):%d %d\n",ii, ij);
           	   


	   for(kk=0; kk<MAXTURNS;kk++) {
	     
	     doubleToSend[kk]=matrix[ii][kk];  
	     doubleToSend[kk+2*MAXTURNS]=0.0;
	     doubleToSend[kk+MAXTURNS]=matrix[ij][kk]; /* BUG to solve TUNES TOO CLOSE  */ 
	     doubleToSend[kk+3*MAXTURNS]=0.0;
	     /*printf("%d %f %f\n", kk,doubleToSend[kk], doubleToSend[kk+MAXTURNS]);*/
	     }
	   /*printf("Hello %d %d",i, maxcounthv);*/
	   
	   sussix4drivenoise_(&doubleToSend[0], &tune[0], &amplitude[0], &phase[0], &allfreqsx[0], &allampsx[0], &allfreqsy[0], &allampsy[0], path);
	   
	   /*for(kk=0;kk<300;kk++){
	     fprintf(lf,"%d %e %e %e %e\n", kk, allfreqsx[kk], allampsx[kk], allfreqsy[kk], allampsy[kk]);
	     }*/


	   allbpmphase[ii]=phase[0];
	   allbpmphase[ij]=phase[3];
	   /* printf("%d %e %e\n", i,  phase[3],  doubleToSend[56+MAXTURNS]);*/
	   allbpmamp[ii]=amplitude[0];
	   allbpmamp[ij]=amplitude[3];
	   /*label[ii]=BPMstatus(1,ii);*/    /*   COMMENTED FOR TEST - REMEMBER TO UNCOMMENT!!!!   */
	   if (labelrun==1) fprintf(nf, "1 %d  %e %e %e %e %e %d %d %f\n", ii, noise1, noiseAve, maxpeak, maxfreq, maxmin, nslines, label[i], phase[0]/360.); 
	   
	   /*printf("%d %d %d %d %e\n", i, label[i],ij, label[ij], phase[3]);*/

	  	   /* PRINT LINEAR FILE */
	   if (amplitude[0]>0 && label[i]==1 && ii==i){
	     fprintf(linfx,"\n\"%s\" %e %d %d %e %e %e %e %e %e %e %e %e %e %e", bpmname[ii], bpmpos[ii],ii, label[ii], tune[0], phase[0]/360., amplitude[0], noise1, maxmin, amplitude[2]/amplitude[0], phase[2]/360.,co, co2, amplitude[1]/amplitude[0], phase[1]/360.);
	     
	       count0[0]++;
	       tunesum[0] += tune[0];
	       tune2sum[0] += tune[0]*tune[0];
	     

	       /* Horizontal Spectrum output*/
	       if (i<10) {
		 sprintf(bsfile,"%s/%s.x",argv[1],bpmname[i]);
	   fBS=fopen(bsfile,"w");
               fprintf(fBS,"%s %s %s\n","*","FREQ","AMP");
	       fprintf(fBS,"%s %s %s\n","$","%le","%le");

	   for(kk=0;kk<300;kk++){
	     fprintf(fBS,"%e %e\n", allfreqsx[kk], allampsx[kk]);
	   }
	   fclose(fBS);
               }



	     /*fprintf(linf, "%e %e ", amplitude[1]/amplitude[0], phase[1]-phase[0]);
	     fprintf(linf, "%e  ", amplitude[2]/amplitude[0]);
	     fprintf(linf, "%e  ", amplitude[6]/amplitude[0]);
	     fprintf(linf, "%e  ", amplitude[7]/amplitude[0]);
	     fprintf(linf, "%e  %e ", amplitude[9]/amplitude[0], phase[9]-2*phase[0]);*/
	   }
	   /* else fprintf(linfx,"\n%s %e %d %d non non non ",  bpmname[i], bpmpos[i] , i, label[i]);*/

	   label[ij]=BPMstatus(2,ij);
	   if (labelrun==1) fprintf(nf, "2 %d  %e %e %e %e %e %d %d %f\n", ij, noise1, noiseAve, maxpeak, maxfreq, maxmin, nslines, label[ij], phase[3]/360.);


	   if (amplitude[3]>0 && label[ij]==1 && ij==i+MAXPICK/2){
	     fprintf(linfy,"\n\"%s\" %e %d %d %e %e %e %e %e %e %e %e %e", bpmname[ij], bpmpos[ij], ij, label[ij], tune[1], phase[3]/360., amplitude[3], noise1, maxmin,  amplitude[5]/amplitude[3], phase[5]/360.,co, co2);

	     

	     count0[1]++;
	       tunesum[1] += tune[1];
	       tune2sum[1] += tune[1]*tune[1];
	     

	       
	       if (ij< MAXPICK/2+10){
		 sprintf(bsfile,"%s/%s.y",argv[1],bpmname[ij]);
	   fBS=fopen(bsfile,"w");
	      fprintf(fBS,"%s %s %s\n","*","FREQ","AMP");
	      fprintf(fBS,"%s %s %s\n","$","%le","%le");
	   for(kk=0;kk<300;kk++){
	     fprintf(fBS,"%e %e \n", allfreqsy[kk], allampsy[kk]);
	   }
	   fclose(fBS);
	       }

	     /*fprintf(linf, "%e  ", amplitude[5]/amplitude[3]); 
	     fprintf(linf, "%e  ", amplitude[4]/amplitude[3]);
	     fprintf(linf, "%e  ", amplitude[11]/amplitude[3]);
	     fprintf(linf, "%e  ", amplitude[13]/amplitude[3]);*/
	   }
	   /* else fprintf(linfy,"\n%s %e %d non non non ", bpmname[ij], bpmpos[ij], label[ij]);*/
	   /* if (amplitude[0]>0){
	     fprintf(linf, "%e  ", amplitude[12]/amplitude[0]);
	   }
	   else fprintf(linf, "non  ");
	   fprintf(linf,"%e %e %e %e ", bpmpos[i], bpmpos[ij], peaktopeak[i], peaktopeak[ij]);
	   if (amplitude[0]>0){
	     fprintf(linf, "%e  ", phase[2]);
	   }
	   else fprintf(linf, "non  ");
	   if (amplitude[3]>0){
	     fprintf(linf, "%e  ", phase[5]);
	   }
	   else fprintf(linf, "non  ");*/
	 }
	 fprintf(linfx,"\n");
	 fprintf(linfy,"\n");
	 fclose(linfx); fclose(linfy);
	 /* What follows is some dirty operations to move the @ Q1 ... to the top of _linx _liny*/
	 if (count0[0]>0) { 
	   sprintf(cmd, "mv %s %s/tx", linfilex,argv[1]); 
	   system(cmd); 
	 } 
	 if (count0[1]>0) {
	   sprintf(cmd, "mv %s %s/ty", linfiley,argv[1]);
	   system(cmd);
	 }
	 

	 
	 if (count0[0]>0) {
	   linfx=fopen(linfilex, "w");
	   fprintf(linfx,"@ Q1 %%le %le\n@ Q1RMS %%le %le\n", tunesum[0]/count0[0], sqrt( tune2sum[0]/count0[0] - (tunesum[0]/count0[0])*(tunesum[0]/count0[0])));
	   fclose(linfx);
	   sprintf(cmd, "cat  %s/tx >> %s ; rm %s/tx", argv[1], linfilex, argv[1]);
	   system(cmd);
	 }
	 if (count0[1]>0) {
	   linfy=fopen(linfiley, "w");
	   fprintf(linfy,"@ Q2 %%le %le\n@ Q2RMS %%le %le\n", tunesum[1]/count0[1], sqrt( tune2sum[1]/count0[1] - (tunesum[1]/count0[1])*(tunesum[1]/count0[1])));
	   fclose(linfy);
	   sprintf(cmd, "cat  %s/ty >> %s ; rm %s/ty", argv[1], linfiley, argv[1]);
	   system(cmd);
	 }
	 
	 




/*     1     2     3       4   5      6      7   8       9      10    11   12*/
/*lin:HBPM LABEL HTUNE HPHASE HAMP -20AMP -20PH 01AMP -30AMP -40AMP 20AMP 20PH*/
/*    LABEL VTUNE VAMP VPHASE 10AMPv 0-2AMPv 0-3AMPv  -1-1AMPv   02AMPh POSH  POSV PPH PPV 01PH 10PHv*/ 
/*     13   14    15    16     17     18      19       20         21    22    23   24  25  26   27*/

	 if (labelrun==1) { fclose(nf) ; fclose(lf) ; }
	 exit(0);
	 tot=0;
 
	 for(i=0; i<counth; i++) tot += label[i];
	 for(i=MAXPICK/2; i<countv; i++) tot += label[i];
	 printf("Number of working pickups: %d\n",tot);
        




         if (0==1) { /* NEW */
	 /*  Starting Analysis */
	 for(mini=0; mini<=1; mini++){
		  
		  count0[mini]=0;
		  count1[mini]=0;
		  count2[mini]=0;
		  count4[mini]=0;
		  count5[mini]=0;
		  count6[mini]=0;
		  count8[mini]=0;
		  count9[mini]=0; 
		  tune[mini]=0.0;
		  tunesum[mini]=0.0;
		  tune2sum[mini]=0.0;
		  amp1sum[mini]=0.0;
		  amp2sum[mini]=0.0;
		  amp4sum[mini]=0.0;
		  amp4sum2[mini]=0.0;
		  amp5sum[mini]=0.0;
		  amp5sum2=0.0;
		  amp20sum[mini]=0.0;
		  phase02sum[mini]=0.0;
		  amp1sum2[mini]=0.0;
		  amp2sum2[mini]=0.0;
		  amp20sum2[mini]=0.0;
		  amp00sum[mini]=0.0;
		  amp00sum2[mini]=0.0;
	 }
	 
	 finii=0;
	 finij=0;
	 i=pickstart;
	 ij=i+MAXPICK/2;
	 printf("Starting pick-ups(H, V): %d %d\n",i , ij);
	 
	 while((finii==0 || finij==0) && i<(MAXPICK/2-1) && ij<(MAXPICK-1)){
	   for(mini=0; (label[i+mini]==0 || label[i+mini+1]==0) && 
		 (i+mini) < pickend && finii==0 && (i+mini+1)<MAXPICK/2; mini++);
	   if ((i+mini+1) >= pickend || (i+mini+1) >= MAXPICK/2) {
	     if (label[i+1]==0 && finii==0) i--; 
	     finii=1;
	   }
	   else {
	     i=i+mini;
	   }
	   for(mini=0; (label[ij+mini]==0 || label[ij+mini+1]==0) && 
		 (ij+mini-MAXPICK/2) < pickend && finij==0 
		 && (ij+mini+1)<MAXPICK; mini +=1);
	   if ((ij+mini+1-MAXPICK/2) >= pickend || (ij+mini+2) >= MAXPICK) {
	     if (label[ij+1]==0 && finij==0) ij--;
	     if (ij<MAXPICK/2) { /*No vertical data*/
	       ij=MAXPICK/2+1;     /*substituted by horizontal data*/
	       label[ij]=1;
	       label[ij+1]=1;
	       NoVerticalData=1;
	       for(mini=0; (label[mini]==0); mini++); 
	       for(jj=0;jj<turns;jj++)
		 matrix[ij][jj]=matrix[mini][jj];
	       mini++;
	       for(mini=mini; (label[mini]==0); mini++);
	       for(jj=0;jj<turns;jj++)
	       matrix[ij+1][jj]=matrix[mini][jj];
	       
	     }
	     finij=1;     
	   }
	   else {
	     ij=ij+mini;
	   }
	   printf("Pick-ups and Labels: %d %d %d %d %d %d fini: %d %d NoVertical:%d\n", 
		  i, ij, label[i],label[i+1], label[ij], label[ij+1], finii, 
		  finij,NoVerticalData);
	   backtune[0]=tune[0];
	   backtune[1]=tune[1];
	   for(mini=0; mini<=13; mini++)  
	     backamp[mini]=amplitude[mini];
	   if ((finii+finij)<2 || (label[i]+label[i+1]+label[ij]+label[ij+1])==4){
	     backlai=label[i];
	     backlaij=label[ij];
	     /*printf("Normalisation in 2: %d\n", normalisation);*/
	     if (norm(i, ij)==0) {
	       printf("\n!!!!!!!!!!! Small Phase advance between pick--ups %d,%d!!!!!!!!!!!!\n",
		    i, ij);
	       if (finii== 0) i += 1;
	       if (finij== 0) ij +=1;
	       continue;
	     }
	   
	     /*savefort300(i);*/
	     printf("Tune: %e %d %e %d %e \n",  tune[0], i, amplitude[0], ij ,tune[1]);
	     /*if (labelrun==1) {
	       if(BPMstatus(1)==0) fprintf(lf, "%d  ", i);
	       if(BPMstatus(2)==0) fprintf(lf, "%d  ", ij);
	       }*/
	     
	     
	     /* Horizontal */
/*     1     2     3       4   5      6      7   8       9      10    11   12*/
/*fb: HBPM COUNT HTUNE HPHASE HAMP -20AMP -20PH 01AMP -30AMP -40AMP 20AMP 20PH*/
/*    VBPM VTUNE VAMP VPHASE 10AMPv 0-2AMPv 0-3AMPv  -1-1AMPv   02AMPh*/ 
/*     13   14    15    16     17     18      19       20         21*/

	     if ((backtune[0] !=tune[0] || backamp[0]!=amplitude[0]) 
		 && amplitude[0]>0 && label[i]==1){
	       count0[0]++;
	       tunesum[0] += tune[0];
	       tune2sum[0] += tune[0]*tune[0];
	       meantune=tunesum[0]/count0[0];
	       fprintf(fb,"\n%d %d %e %e %e ", i, count0[0], tune[0], phase[0], amplitude[0]);
	       
	       
	       if (backamp[1] != amplitude[1] &&  amplitude[1] > 0){
		 pickupaxis1[count1[0]][0] = i;
		 amp30bpm[count1[0]] = amplitude[1]/amplitude[0];
		 amp1sum[0]=amp1sum[0]+amp30bpm[count1[0]];
		 amp1sum2[0]=amp1sum2[0]+amp30bpm[count1[0]]*amp30bpm[count1[0]];
		 fprintf(fb, "%e %e ", amp30bpm[count1[0]], phase[1]-phase[0]);
		 count1[0]++;
		 
	       }
	       else fprintf(fb,"non  non ");
	       
	       if (backamp[2] != amplitude[2] &&  amplitude[2]>0 && phase[2]<=180
		   && NoVerticalData==0){
		 pickupaxis2[count2[0]][0]=i;
		 amp01bpm[count2[0]]=amplitude[2]/amplitude[0];
		 amp2sum[0]=amp2sum[0]+amp01bpm[count2[0]];
		 amp2sum2[0]=amp2sum2[0]+amp01bpm[count2[0]]*amp01bpm[count2[0]];
		 fprintf(fb, "%e  ", amp01bpm[count2[0]]);
		 count2[0]++;		   
	       }
	       else   fprintf(fb,"non  ");
	       
	       if (backamp[6] != amplitude[6] &&  amplitude[6]>0 && phase[6]<=180){
		 
		 amp40bpm[count4[0]]=amplitude[6]/amplitude[0];
		 amp4sum[0]=amp4sum[0]+amp40bpm[count4[0]];
		 amp4sum2[0]=amp4sum2[0]+amp40bpm[count4[0]]*amp40bpm[count4[0]];
		 fprintf(fb, "%e  ", amp40bpm[count4[0]]);
		 count4[0]++;
		 
	       }
	       else fprintf(fb,"non  ");
	       if (backamp[7] != amplitude[7] &&  amplitude[7]>0 && phase[7]<=180){
		 amp50bpm[count5[0]]=amplitude[7]/amplitude[0];
		 amp5sum[0]=amp5sum[0]+amp50bpm[count5[0]];
		 amp5sum2=amp5sum2+amp50bpm[count5[0]]*amp50bpm[count5[0]];
		 fprintf(fb, "%e  ", amp50bpm[count5[0]]);
		 /*amp2sum2[0]=amp2sum2[0]+amp01bpm[count2[0]]*amp01bpm[count2[0]];*/
		 count5[0]++;
	       }
	       else fprintf(fb,"non  ");
	       if (backamp[8] != amplitude[8] &&  amplitude[8]>0 && phase[8]<=180){
		 amp60bpm[count6[0]]=amplitude[8]/amplitude[0];
		 amp6sum[0]=amp6sum[0]+amp60bpm[count6[0]];
		 /*amp2sum2[0]=amp2sum2[0]+amp01bpm[count2[0]]*amp01bpm[count2[0]];*/
		 count6[0]++;
	       }
	       if (backamp[9] != amplitude[9] &&  amplitude[9]>0 && phase[9]<=180){
		 amp20bpm[count8[0]]=amplitude[9]/amplitude[0];
		 amp20sum[0]=amp20sum[0]+amp20bpm[count8[0]];
		 amp20sum2[0]=amp20sum2[0]+amp20bpm[count8[0]]*amp20bpm[count8[0]]; 
		 fprintf(fb, "%e  %e ", amp20bpm[count8[0]], phase[9]-2*phase[0]);
		 count8[0]++;
		 
	       } else   fprintf(fb,"non  non ");
	       
	       if (backamp[10] != amplitude[10] &&  amplitude[10]>0 && phase[10]<=180){ 
		 amp00bpm[count9[0]]=amplitude[10];
		 amp00sum[0]=amp00sum[0]+amp00bpm[count9[0]];
		 amp00sum2[0]=amp00sum2[0]+amp00bpm[count9[0]]*amp00bpm[count9[0]];
		 count9[0]++;
	       }      
	     } else 
	       fprintf(fb,"\n%d  non  non  non  non  non  non  non  non  non  non  non  ", i);
	     
	     
	     /* Vertical */
	     if ((backtune[1] !=tune[1] || backamp[3] !=amplitude[3]) 
		 && amplitude[3] >0 && label[ij]==1 && NoVerticalData==0){
	       count0[1]++;
	       tunesum[1]=tune[1]+tunesum[1];
	       tune2sum[1]=tune2sum[1]+tune[1]*tune[1];
	       fprintf(fb,"%d  %e  %e  %e ", ij, tune[1], amplitude[3], phase[3]);
    
	       if (backamp[5] != amplitude[5] &&  amplitude[5]>0 ){
		 pickupaxis2[count2[1]][1]=ij;
		 amp10bpm[count2[1]]=amplitude[5]/amplitude[3];
		 amp2sum[1]=amp2sum[1]+amp10bpm[count2[1]];
		 amp2sum2[1]=amp2sum2[1]+amp10bpm[count2[1]]*amp10bpm[count2[1]];
		 fprintf(fb, "%e  ", amp10bpm[count2[1]]);
		 count2[1]++; 
	       }
	       else fprintf(fb,"non  ");
  
	       if (backamp[4] != amplitude[4] &&  amplitude[4]>0 ){
		 fprintf(fb, "%e  ", amplitude[4]/amplitude[3]);
	       }
               else fprintf(fb,"non  ");
	       
	       if (backamp[11] != amplitude[11] &&  amplitude[11]>0 ){
		 fprintf(fb, "%e  ", amplitude[11]/amplitude[3]);
	       }
               else fprintf(fb,"non  ");
       
	       if (backamp[13] != amplitude[13] &&  amplitude[13]>0 ){
		 fprintf(fb, "%e  ", amplitude[13]/amplitude[3]);
	       }
               else fprintf(fb,"non  ");
	     } 
	     else fprintf(fb,"%d  non  non  non  non  non  non  non",ij);
	     /* An extra  horizontal line at the end of the output string*/ 
	     if (backamp[12] != amplitude[12] &&  amplitude[12]>0 
               && (backtune[0] !=tune[0] || backamp[0]!=amplitude[0])
               && amplitude[0]>0 && label[i]==1){
		 fprintf(fb, "%e  ", amplitude[12]/amplitude[0]);
	     }
               else fprintf(fb,"non  ");
	     /*bpm positions and peaktopeak at the end of the line*/

             fprintf(fb,"%e %e %e %e ", bpmpos[i], bpmpos[ij], peaktopeak[i], peaktopeak[ij]);
	   }
	   if (finii== 0) i += 1;
	   if (finij== 0) ij +=1;
	 }		
	 fclose(fb);
	 
	 printf("\nQx: %e Ampx: %e ", meantune,averamp[0]);

	 /*if (labelrun==1) fclose(lf);*/
	 if (count0[0] > 0) {
	   fprintf(ft, "%e %e %e %e %e %e ", parameter[iparameter],
		   parameter[iparameter]*parameter[iparameter], 
		   amplitude[0], amplitude[3], meantune, 
		   sqrt(tune2sum[0]/count0[0]-meantune*meantune));
	 } else fprintf(ft,"%e non non non non non ", parameter[iparameter]);
	 
	 if (count0[1] > 0) { 
	   printf("Qy: %e Ampy: %e \n", tunesum[1]/count0[1],averamp[1]);
	   fprintf(ft, "%e %e ", tunesum[1]/count0[1],
		    sqrt(tune2sum[1]/count0[1]-(tunesum[1]/count0[1])*(tunesum[1]/count0[1]))); 
	 } else  fprintf(ft,"non non ");
	   
	 fprintf(ft, "%e %e %e %e\n",averamp[0], averamp2[0], averamp[1], averamp2[1]);

	 if (count1[0] > 0) {
	   amp30ave[iparameter]=amp1sum[0]/count1[0];
	   fprintf(fs, "%e %e %e ", parameter[iparameter], amp30ave[iparameter],
		   sqrt(amp1sum2[0]/count1[0]-amp30ave[iparameter]*amp30ave[iparameter]));
	 } else 
	   fprintf(fs, "%e non non ", parameter[iparameter]);
	 
	 if (count8[0] > 0) {
	   amp20ave[iparameter]=amp20sum[0]/count8[0];
	   fprintf(fs, "%e %e ", amp20ave[iparameter],
		   sqrt(amp20sum2[0]/count8[0]-amp20ave[iparameter]*amp20ave[iparameter]));
	   
	 } else 
	   fprintf(fs, "non non ");


	 
	 if (count9[0] > 0) {
	   amp00ave[iparameter]=amp00sum[0]/count9[0];
	   fprintf(fs, "%e %e ", amp00ave[iparameter],
		   sqrt(amp00sum2[0]/count9[0]-amp00ave[iparameter]*amp00ave[iparameter]));
	   
	 } else 
	   fprintf(fs, "non non ");
	 
	 fprintf(fs, "%e %e %e %e\n",averamp[0], averamp2[0], averamp[1], averamp2[1]);

	 
	 if (count2[0] > 0) {
	   amp01ave[iparameter]=amp2sum[0]/count2[0];
	   fprintf(fc, "%e  %e  %e  ", parameter[iparameter], amp01ave[iparameter], 
		   sqrt(amp2sum2[0]/count2[0]-amp01ave[iparameter]*amp01ave[iparameter]));
	 } else 
	   fprintf(fc, "%e non non ", parameter[iparameter]);
	 if (count2[1] > 0) { 
	   fprintf(fc, "%e  %e  %e\n",amp2sum[1]/count2[1], 
		   sqrt(amp2sum2[1]/count2[1]-(amp2sum[1]/count2[1])*(amp2sum[1]/count2[1])),
		   sqrt(amp01ave[iparameter]*amp2sum[1]/count2[1]));
	 } else  
	   fprintf(fc, "non non non\n");
	 
	 
	 if (count4[0] > 0) {
	   amp40ave[iparameter]=amp4sum[0]/count4[0];
	   fprintf(foc, "%e %e %e\n", parameter[iparameter], amp40ave[iparameter],
		   sqrt(amp4sum2[0]/count4[0]-amp40ave[iparameter]*amp40ave[iparameter]));  
	 }
	 if (count5[0] > 0) {
	   amp50ave[iparameter]=amp5sum[0]/count5[0];
	   fprintf(fdec, "%e %e %e\n", parameter[iparameter], amp50ave[iparameter],
		   sqrt(amp5sum2/count5[0]-amp50ave[iparameter]*amp50ave[iparameter]));
	   
	 }
	 if (count6[0] > 0) {
	   amp60ave[iparameter]=amp6sum[0]/count6[0];
	   fprintf(fdod, "%e %e\n", parameter[iparameter], amp60ave[iparameter]);
	   
	 }
	 }/*NEW*/	 
	 iparameter++;
	 }
	 fclose(fc);
	 fclose(ft);
	 fclose(fs);	 
	 fclose(foc);
	 fclose(fdec);
	 fclose(fdod);
	 if (ibetabeat==1) {
		  printf("\nBETAbeating process running\n");
		  betabeat();
	 }
}


/* ***************** */
/*    savefort300    */
/* ***************** */
void savefort300(int ir)
{
	 char *path="F300", ffile[20], number[4], cmd[100];
	 number[0]=(int) ir/100 + '0';
	 number[1]=(int) ((ir-(number[0]-'0')*100)/10) + '0';
	 number[2]=(int) (ir-(number[0]-'0')*100-(number[1]-'0')*10) + '0';
	 number[3]='\0';
	 sprintf(ffile, "%s/fort.300_%s", path, number);
	 sprintf(cmd, "cp fort.300 %s", ffile);
	 system(cmd);
	 return;
	 
}


/* ***************** */
/*    GET_NAME       */
/* ***************** */
void get_name()
{
	 char s[1000];
	 int i=0;
	 FILE *nf;
	 nf=fopen(namefile, "r");
	 while (i<iparameter){
		  if ((datafile[0]=getc(nf))== '\n') i++;
	 }
	 for (i=0; isspace(datafile[0]=getc(nf)); i++);
	 for (i=1; ((datafile[i]=getc(nf)) != EOF) && 
			   (datafile[i]  !='\n') && (datafile[i]  !='%') &&
			   (datafile[i]  !=' ')  ; i++);
	 
	 if (datafile[i]==EOF) {
		  scum=1;
		  datafile[i]='\0';
		  iparameter++;
		  return;
	 }
	 datafile[i]='\0';
	 for (i=0; isspace(s[0]=getc(nf)); i++);
	 for (i=1; ((s[i]=getc(nf)) != EOF) & 
                (s[i] !='\n') &
                (s[i] !=' ')  ; i++);
	 if (s[i]==EOF) scum=1;
	 s[i]='\0';
	 parameter[iparameter]=atof(s);
	 
	 for (i=0; isspace(s[0]=getc(nf)); i++);
	 for (i=1; ((s[i]=getc(nf)) != EOF) & 
                (s[i] !='\n') &
                (s[i] !=' ')  ; i++);
	 if (s[i]==EOF) scum=1;
	 s[i]='\0';
	 turns=atof(s); 
	 /*
	 for (i=0; isspace(s[0]=getc(nf)); i++);
	 for (i=1; ((s[i]=getc(nf)) != EOF) & 
                (s[i] !='\n') &
                (s[i] !=' ')  ; i++);
	 if (s[i]==EOF) scum=1;
	 s[i]='\0';
	 tuney=atof(s);*/


	 fclose(nf);

	 

}

/* ***************** */
/*    LECTOR         */
/* ***************** */
/*double lector(int l, char separator, char *file)
{
	int i=0;
	FILE *di;
	char s[80];
	di=fopen(file, "r"); 
  	while (i<l){
		 if ((s[0]=getc(di))== separator) i++;
	}
	 for (i=0; ((s[i]=getc(di)) != EOF) & 
                (s[i]  !='\n') &
                (s[i]  !=' ')  ; i++);
	 s[i]='\0';
     fclose(di);
	 return atof(s);
}*/
/* ***************** */
/*    GETPARAMETER   */
/* ***************** */
double getparameter()
{
	int i=0; 
    char s[80];
	printf("\n Insert parameter (value for the x axis) ");
	while ((s[i]=getchar()) != '\n' ) i++;
	s[i]='\0';
	return atof(s);
}

/* ***************** */
/*      NORM                       */
/* ***************** */
int norm(int i,int j)
{
	 int jj, kk;
	 double p1, p2, a1, a2, t1, t2, p21, pp1, pp2, aa1, aa2, pp21, sinp21, 
	   sinpp21, cosp21, cospp21;
	 FILE *f9;

/*Creating sussix_v4.inp for 1 dimension */
	 printf("Normalisation in norm: %d\n", normalisation);
	 if (normalisation == 1) {
       /*sussix_inp(1);
       Creating fort.90 as input of sussix4Drive preanalysis
       f9=fopen("fort.90", "w");
       for (jj=0;  jj < turns; jj++){
	 fprintf( f9, "%e %e %e %e\n", matrix[i][jj], matrix[i+1][jj], matrix[j][jj], 
		 matrix[j+1][jj]);
       } 
       fclose(f9); */
       /*       sussix4drive_(&tune[0], &amplitude[0], &phase[0]);*/
       /*t1=tune[0]*/
       a1=allbpmamp[i];
       p1=allbpmphase[i];
       aa1=allbpmamp[j];
       pp1=allbpmphase[j];
       
       printf("Amps: %e %e\n", a1, aa1);
       /* Creating fort.90 as input of sussix4Drive preanalysis */	 
       /*f9=fopen("fort.90", "w");
       for (jj=0; jj< turns; jj++){
	  fprintf(f9, "%e %e %e %e\n", matrix[i+1][jj], matrix[i][jj], matrix[j+1][jj], 
		  matrix[j][jj]);
       } 
       fclose(f9); */
       /*       sussix4drive_(&tune[0], &amplitude[0], &phase[0]);*/
       /*t2=tune[0];*/
       a2=allbpmamp[i+1];
       p2=allbpmphase[i+1];
       aa2=allbpmamp[j+1];
       pp2=allbpmphase[j+1];
       /* Creating fort.90 as input of sussix4Drive ANALYSIS*/
       p21=M_PI*(p2-p1)/180.0;
       pp21=M_PI*(pp2-pp1)/180.0;
       sinp21=sin(p21);
       sinpp21=sin(pp21);
       /* Cut those pick-ups with phase advance different from 90d in more than 60d */
       printf("Phase increase, amps : %e %e %e %e %e %e\n", p21, pp21, a1, a2, aa1, aa2);
       printf("First points: %e %e %e %e\n", matrix[i][0], matrix[i+1][0], 
	      matrix[j][0],matrix[j+1][0]);

       if (sinp21*sinp21<0.25 || a2<=0.0 || a1<=0.0 )  
	 label[i]=0; 
       if (sinpp21*sinpp21<0.25 || aa2<=0.0 || aa1 <=0.0)
	 label[j]=0; 
	   
       if (label[i]==0 && label[j]==0) return 0;
       
       cosp21=cos(p21);
       cospp21=cos(pp21);

       if(label[i]==0){
	 printf("X Not normalisable");
	 a2=1.0;
	 cosp21=0;
	 sinp21=1.0;
       }

       if(label[j]==0){
	 printf("\n Y Not normalisable \n");
	 aa2=1.0;
	 cospp21=0;
	 sinpp21=1.0;
       }

       for(kk=0; kk<MAXTURNS;kk++) {
	     doubleToSend[kk]=matrix[i][kk];  
	     doubleToSend[kk+2*MAXTURNS]=(matrix[i+1][kk]*a1/a2-cosp21*matrix[i][kk])/sinp21;
	     doubleToSend[kk+MAXTURNS]=matrix[j][kk]; 
	     doubleToSend[kk+3*MAXTURNS]=(matrix[j+1][kk]*aa1/aa2-cospp21*matrix[j][kk])/sinpp21;
	   }
 
	      /*      f9=fopen("fort.90", "w");
       for (jj=0; jj< turns; jj++){ 
	 if (label[i]==1 && label[j]==1)
	   fprintf(f9, "%e %e %e %e\n", matrix[i][jj],  
		   (matrix[i+1][jj]*a1/a2-cosp21*matrix[i][jj])/sinp21
		   , matrix[j][jj], 
		   (matrix[j+1][jj]*aa1/aa2-cospp21*matrix[j][jj])/sinpp21);
	 if (label[i]==1 && label[j]==0) 
	   fprintf(f9, "%e %e %e %e\n", matrix[i][jj], 
		   (matrix[i+1][jj]*a1/a2-cosp21*matrix[i][jj])/sinp21
		   , matrix[j][jj], matrix[j+1][jj]);
	 if (label[i]==0 && label[j]==1)
	   fprintf(f9, "%e %e %e %e\n", matrix[i][jj], matrix[i+1][jj], matrix[j][jj],
		   (matrix[j+1][jj]*aa1/aa2-cospp21*matrix[j][jj])/sinpp21);
       }
       fclose(f9);*/
     }
     else {
         for(kk=0; kk<MAXTURNS;kk++) {
	     doubleToSend[kk]=matrix[i][kk];  
	     doubleToSend[kk+2*MAXTURNS]=matrix[i+1][kk];
	     doubleToSend[kk+MAXTURNS]=matrix[j][kk]; 
	     doubleToSend[kk+3*MAXTURNS]=matrix[j+1][kk];
	   }

       /* f9=fopen("fort.90", "w");  
       for (jj=0; jj< turns; jj++){  
	fprintf(f9, "%e %e %e %e\n", matrix[i][jj], matrix[i+1][jj], matrix[j][jj], 
		 matrix[j+1][jj]); 
       } 
       fclose(f9);*/ 
     }
     /*Creating sussix_v4.inp for complex signal*/

	 /*	 sussix_inp(iir);     
	   sussix4drivenoise_(&doubleToSend[0], &tune[0], &amplitude[0], &phase[0], &allfreqsx[0], &allampsx[0], &allfreqsy[0], &allampsy[0]);*/

     /*     sussix4drive_(&tune[0], &amplitude[0], &phase[0]);*/
     return 1;
}
	 
/* ***************** */
/*    sussix_inp     */
/* ***************** */
void sussix_inp(int ir, char *path)
{
	 FILE *si;
	 int idam;
	 /*char thefile[4000];
	 strcpy(thefile, path);
	 strcat(thefile, "sussix_v4.inp");*/
	 si=fopen(path, "w");
	 idam=2;
	 fprintf(si, "C\nC INPUT FOR SUSSIX_V4 ---17/09/1997---\n");
	 fprintf(si,"C DETAILS ARE IN THE MAIN PROGRAM SUSSIX_V4.F\nC\n\n");
	 fprintf(si,"ISIX  = 0\nNTOT  = 1\nIANA  = 1\nICONV = 0\n");
	 fprintf(si,"TURNS = 1 %d\n", turns);
	 fprintf(si,"NARM  = 300\nISTUN = 1 %e %e\n", istun, istun);
	 fprintf(si,"TUNES = %e %e .07\n", tunex, tuney);
	 fprintf(si,"NSUS  = 0\nIDAM  = %d\n", idam);
	 fprintf(si,"NTWIX = 1\nIR    = %d\nIMETH = 2\nNRC   = 6\nEPS   = 5D-3\n",
		 ir); /* EPS is the window in the secondary lines, very imp!!!*/
	 fprintf(si,"NLINE = 0\nL,M,K = \nIDAMX = 1\nNFIN  = 500\nISME  = 1\n");
	 fprintf(si,"IUSME = 200\nINV   = 0\nIINV  = 250\nICF   = 0\nIICF  = 350\n");
     fclose(si);     
}
void get_nameb();
char *namefileb="BPMfiles", datafileb[100];
int iparameterb=0, scumb=0;

/* ***************** */
/*    BetaBeat       */
/* ***************** */
void betabeat()
{
	 char s[100];
	 int i=0,j,ii=0,jj=0, jjx=0, jjy=0, flag=0, contx[MAXPICK], conty[MAXPICK],
		  xpick[MAXPICK], ypick[MAXPICK];
	 double ampx[MAXPICK], ampy[MAXPICK], sumampx[MAXPICK], sum2ampx[MAXPICK], 
		  sumampy[MAXPICK], sum2ampy[MAXPICK], ax, ay;
	 FILE *df, *bf;
	 for (j=0;j<MAXPICK; j++){
		  contx[j]=0;
		  conty[j]=0;
		  xpick[j]=0;
          ypick[j]=0;
		  sumampx[j]=0;
		  sum2ampx[j]=0;
		  sumampy[j]=0;
		  sum2ampy[j]=0;
		  ampx[j]=0;
		  ampy[j]=0;
	 }
	 
	 if ( access(namefileb, 0) < 0 )
		  printf("no file %s", namefileb);
	 
	 while(scumb==0){
		  
		    get_nameb();
		  if ( access(datafileb, 0) < 0 ) {
			   printf("No file %s", datafileb);
			   break;
		  }else
			   printf("opening %s\n", datafileb);
		  
		  df=fopen(datafileb, "r");
		  jj=0;
		  i=0;
		  ii=0;
		  jjx=0;
		  jjy=0;
		  ax=0;
		  ay=0;
		  
		  while(isspace(s[0]=getc(df)));
		  
		  while(s[0] != EOF ){
			   if (s[0]=='\n') {jj++; ii=0;}
			   if ((s[0]=='\n' || s[0]==' ' || s[0]=='\t') && flag==1) flag=0;
			   if ((s[0] !='\n' && s[0] !=' ' && s[0] !='\t') && flag==0){
					while(!isspace(s[i]) && (s[i] !=EOF)){	 
 						 i++;
						 s[i]=getc(df);
					}
					s[i+1]=s[i];
					s[i]='\0';
					
					if (ii==0) {
						 xpick[jjx]=atof(s);
 					 
					}
					if (ii==2 && s[i-1] !='n'){
						 ampx[xpick[jjx]]=atof(s);
						 ax=ax+ampx[xpick[jjx]];
						 printf("ASSIGN %e\t", ampx[jjx]);
 						 jjx++;
					}
					if (ii==5){
						 ypick[jjy]=atoi(s);
					}
					if (ii==7 && s[i-1] !='n'){
						 ampy[ypick[jjy]]=atof(s);
						 ay=ay+ampy[ypick[jjy]];
						 jjy++;
					}
					
					ii++;
					flag=1;
					s[0]=s[i+1];
					i=0;
			   }
			   if (flag==0) s[0]=getc(df);
		  }
		  fclose(df);  
		  ay=ay/jjy;
		  ax=ax/jjx;
		  
		  for (i=0; i<jjx; i++){
			   sumampx[xpick[i]]=sumampx[xpick[i]]+ampx[xpick[i]]/ax;
			   sum2ampx[xpick[i]]=sum2ampx[xpick[i]]
					+ampx[xpick[i]]*ampx[xpick[i]]/(ax*ax);
			   contx[xpick[i]]++;
					
		  }
		  for (i=0; i<jjy; i++){
			   sumampy[ypick[i]]=sumampy[ypick[i]]+ampy[ypick[i]]/ay;
			   sum2ampy[ypick[i]]=sum2ampy[ypick[i]]
					+ampy[ypick[i]]*ampy[ypick[i]]/(ay*ay);
		  conty[ypick[i]]++;
		  }
		  iparameterb++;
	 }
	bf=fopen("betabeating","w"); 
	for (i=0; i<=jjx ; i++){
		 fprintf(bf, "%d  %e  %e  %d  %e  %e\n",  
				 xpick[i], sumampx[xpick[i]]/contx[xpick[i]], 
				 sqrt(sum2ampx[xpick[i]]/contx[xpick[i]]-
					  sumampx[xpick[i]]*sumampx[xpick[i]]/(contx[xpick[i]]*
														   contx[xpick[i]])),
				 ypick[i], sumampy[ypick[i]]/conty[ypick[i]], 
				 sqrt(sum2ampy[ypick[i]]/conty[ypick[i]]-
					  sumampy[ypick[i]]*sumampy[ypick[i]]/(conty[ypick[i]]*
														   conty[ypick[i]])));
		 
	}	
	fclose(bf); 	
}

/* ***************** */
/*    GET_NAMEB      */
/* ***************** */
void get_nameb()
{
	 char s[80];
	 int i=0;
	 FILE *nf;
	 nf=fopen(namefileb, "r");
	 while (i<iparameterb){
		  if ((datafileb[0]=getc(nf))== '\n') i++;
	 }
	 for (i=0; isspace(datafileb[0]=getc(nf)); i++);
	 for (i=1; ((datafileb[i]=getc(nf)) != EOF) & 
			   (datafileb[i]  !='\n') && (datafileb[i]  !='%') &
			   (datafileb[i]  !=' ')  ; i++);
	 if (datafileb[i]==EOF) {
		  scumb=1;
		  datafileb[i]='\0';
		  return;
	 }
	 datafileb[i]='\0';
	 fclose(nf);
}

/************   GET_WORD *************/
/* get next string in between blanks */
/*************************************/

char *get_word(FILE *ff){
	 /*char ss[100];*/
  int i=1;
  while(isspace(ss[0]=getc(ff))){}
  while(!isspace(ss[i]=getc(ff))){i++;}
  ss[i]='\0';
  return &ss[0];
}



/************   BPMstatus *************/
/* Analyse fort.300 to detect noise   */
/**************************************/
#define MINSIGNAL 0.00001
#define SIGMACUT   1.8
#define MINIMUMNOISE 0.0
#define BADPICKUP  8.0
#define MAXSIGNAL 30000
int BPMstatus(int plane, int pickupN)
{
	 double norma=0.0, weight, aux, freq=0.0, freq2=0.0,margin=0.1, ave=0, amp,
		  maxe=-500000.0, mine=500000.0;
	 int counter, il, counter2=0, counter3=0, status=1;
	 FILE *fort, *f90;
	 char line[500], ssw[80], str[80];

	 maxpeak=0; /*Initialising*/
	 co=0.0;	 
         co2=0.0;
	 /* If peak-to-peak signal smaller than MINSIGNAL reject */
	 if(plane==1){
	   for(il=0; il<turns; il++){
		 co += doubleToSend[il];
		 co2 += doubleToSend[il]*doubleToSend[il];
	     /*if(doubleToSend[il]==99) return 0;  THIS IS ONLY FOR RHIC*/
	     if(doubleToSend[il]< mine) mine=doubleToSend[il];
	     if(doubleToSend[il] > maxe) maxe=doubleToSend[il];		  
	   }
	 }
	 if(plane==2){
	   for(il=MAXTURNS; il<MAXTURNS+turns; il++){
		 co += doubleToSend[il];
		 co2 += doubleToSend[il]*doubleToSend[il];
             /*if(doubleToSend[il]==99) return 0;  THIS IS ONLY FOR RHIC*/
	     if(doubleToSend[il]< mine) mine=doubleToSend[il];
	     if(doubleToSend[il] > maxe) maxe=doubleToSend[il];		  
	   }
	 }
	 co=co/turns;
	 co2=sqrt(co2/turns-co*co);
	 maxmin=maxe-mine;
	 if (maxmin<MINSIGNAL || maxmin>MAXSIGNAL) return 0;
	 /*printf("\n1\n");*/
       
	 
	 /* Compute the spread and average in the intervals [windowa1,windowa2]
            and [windowb1,windowb2] */
	 
	 /*fort=fopen("fort.300", "r");*/
	 noise1=0;
	 
	 for(counter=0; counter<300; counter++){
	   	   
	   if(plane==1) {aux=allfreqsx[counter];amp=allampsx[counter];}
	   if(plane==2) {aux=allfreqsy[counter];amp=allampsy[counter];}

	   if (amp>maxpeak && aux>0.05) {maxpeak=amp; maxfreq=aux;}

	   if ((aux<windowa2 && aux>windowa1) || (aux<windowb2 && aux>windowb1)){
	     

	     if (amp<0) { /* Something in sussix went wrong */
	       /*fclose(fort);*/
	       noise1=100;
	       noiseAve=100;
	       maxpeak=100; maxfreq=100;
	       return 0;
	     }

	     ave=amp+ave;
	     noise1=noise1+amp*amp;
	     counter3++;
	   }
	   
	 }
         /*printf("\n2\n");*/
	 /*fclose(fort);*/
	 if (counter3>0){
	   if(counter3>1)
	     noise1=sqrt((noise1/counter3-ave*ave/(counter3*counter3)));
	   else noise1=0;
	   noiseAve=ave/counter3;
	 } 
	 else {
	   noise1=MINIMUMNOISE;
	   noiseAve=MINIMUMNOISE;
	 }
	 nslines=counter3;
	 printf("%e %e", noiseAve, noise1);
	 
	 /* If tune line isn't larger than background reject*/ 
	  
	 if ((maxfreq<windowa2 && maxfreq>windowa1) || (maxfreq<windowb2 && maxfreq>windowb1))
	   printf("NoiseWindow includes largest lines, amp %e freq %e!!!!\n", maxpeak, maxfreq);
	 /*printf("\n3\n");*/
	 if (maxpeak < noiseAve+SIGMACUT*noise1) return 0;
	 /*printf("\n4\n");*/
	 if (noise1 > BADPICKUP) return 0;

	 /* Otherwise pick-up succeded to first cuts */  
	 return 1;
	 
   
}







/************   BPMstatus *************/
/* Analyse fort.300 to detect noise   */
/**************************************/
int BPMstatus_OLD(int plane)
{
	 double norma=0.0, weight, aux, freq=0.0, freq2=0.0,margin=0.1 ;
	 int counter, il, counter2=0;
	 FILE *fort, *f90;
	 char line[500], ssw[80], str[80];

	 fort=fopen("fort.300", "r");
	 f90=fopen("fort.90", "r");
	 get_word(f90);
	 strcpy(str, get_word(f90));
	 for(il=1; il<=6; il++){
		  get_word(f90);get_word(f90);get_word(f90);
		  strcpy(ssw, get_word(f90));
		  if (atof(str) != atof(ssw)) break;
		  if (il==6) return 0;
		  
	 }
	 fclose(f90);
	 
	 while (plane>0){
		  strcpy(line, get_line(fort)); 
		  if (strstr(line, "ANALYSIS OF") != NULL ) plane-- ;
	 }
	 get_line(fort);get_line(fort);get_line(fort);get_line(fort);
	 for(counter=0; counter<149; counter++){
		  get_word(fort);
		  strcpy(ssw, get_word(fort));
		  aux=atof(ssw);
		  if (aux<tunex+margin && aux>tunex-margin){
			   strcpy(ssw, get_word(fort));
			   weight=atof(ssw)*atof(ssw);
			   norma=norma+weight;
			   freq=freq+aux*weight;
			   freq2=freq2+aux*aux*weight;
			   counter2++;
			   
		  }
		  get_line(fort);
	 }
	 printf ("------------------>> %e \n", sqrt((freq2/norma-freq*freq/(norma*norma))));
	 fclose(fort);
	 if(counter2 <= 5) return 0;
	 if (sqrt((freq2/norma-freq*freq/(norma*norma))) > 0.009 )  return 0;
	 else return 1;
}

/************   GET_LINE *************/
/* get next line of the FILE pointed */
/*************************************/

char *get_line(FILE *ff){
	 /*char ss[1000];*/
  int i=0;
  while((ss[i]=getc(ff))!='\n'){i++;}
  ss[i]='\0';
  return &ss[0];
}

/************  ISIN(A,B) ***********/
/* If B contains A ISIN=1; else 0  */
/**********************************
int isin(char *ss, char *ll){
  char dd[100], word[100], line[1000];
  int countw=0, countl=0, i=0, j=0;
  strcpy(word, ss);
  while(word[countw]!='\0'){
    countw++;}
  strcpy(line, ll);
  while(line[countl]!='\0'){
    countl++;}
  while(strcmp(dd, word)!=0 && i<(countl-countw)){
    for(j=0; j<countw; j++){
      dd[j]=line[i+j];
    }
    dd[countw]='\0';
    
    i++;
  }
  if (i>=(countl-countw)) return 0;
  else return 1;
  }*/










