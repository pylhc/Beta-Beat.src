C=============================================================
C COMPUTES THE TUNE USING FFT INTERPOLATED METHOD.
C X, XP ARE THE COORDINATES OF THE ORBIT AND MAXN IS THE 
C LENGTH OF THE ORBIT.
C
C AUTHOR:     E. TODESCO - INFN AND CERN 
C             M. GIOVANNOZZI - CERN HAS INTRODUCED SOME 
C                                   MODIFICATIONS
C

      DOUBLE PRECISION FUNCTION TUNEABT(X,XP,MAXN,ZW)              
      PARAMETER(MAXITER=100000)
      IMPLICIT DOUBLE PRECISION (A-H,O-Y)
      IMPLICIT COMPLEX*16(Z)
      COMPLEX*16 ZSING(MAXITER)
      DIMENSION X(MAXITER),XP(MAXITER)
      DIMENSION Z(MAXITER)
C.............................................................
      ZU=DCMPLX(0D0,1D0)
C..................................ESTIMATION OF TUNE WITH FFT 
      PI=DATAN(1D0)*4D0
      DUEPI=2*PI
      MFT=INT(LOG(FLOAT(MAXN))/LOG(2D0)) 
      NPOINT=2**MFT
      STEP=DUEPI/NPOINT/2D+0
C.............................................................
      SUM=0D0            !..CHECKS FOR COMPLEX OR REAL DATA
      DO MF=1,NPOINT
        Z(MF)=DCMPLX(X(MF),XP(MF))
        ZSING(MF)=Z(MF)
        SUM=SUM+XP(MF)
      ENDDO 
      CALL FFT(ZSING,NPOINT,-1)
C.......................SEARCH FOR MAXIMUM OF FOURIER SPECTRUM
      NPMIN=1
      IF (SUM.EQ.0D0) THEN
        NPMAX=NPOINT/2         !..REAL FFT ONLY HALF COEFFICIENTS
      ELSE
        NPMAX=NPOINT
      ENDIF
C.............................................................
      FTMAX=0D0
      NFTMAX=0
      DO NFT=NPMIN,NPMAX
        IF (ABS(ZSING(NFT)).GT.FTMAX) THEN
          FTMAX=ABS(ZSING(NFT))
          NFTMAX=NFT
        END IF
      ENDDO
C................................................INTERPOLATION
      CF1=ABS(ZSING(NFTMAX-1))
      CF2=ABS(ZSING(NFTMAX))
      CF3=ABS(ZSING(NFTMAX+1))
      IF (CF3.GT.CF1) THEN      
        ASSK=DFLOAT(NFTMAX)+NPOINT/PI*
     .       ATAN2(CF3*DSIN(PI/NPOINT),CF2+CF3*DCOS(PI/NPOINT))
      ELSEIF (CF3.LE.CF1) THEN                   
        ASSK=DFLOAT(NFTMAX-1)+NPOINT/PI*
     .       ATAN2(CF2*DSIN(PI/NPOINT),CF1+CF2*DCOS(PI/NPOINT))
      ENDIF
C........................................FINAL INTERPOLATION
      TUNEABT=1D+0-(ASSK-1D+0)/DFLOAT(NPOINT)
C............................................................  
      ZTUNE=EXP(-ZU*DUEPI*TUNEABT)
      CALL CALCR(ZTUNE,ZW,Z,MAXN)
C............................................................  
      RETURN 
C............................................................  
      END
