C Parallelised version using Intel IA32 32-bit  compiler. On CERN lxplus needs:
C  source /afs/cern.ch/sw/IntelSoftware/linux/all-setup.csh ia32
C 
C  For the number of cores enter: setenv OMP_NUM_THREADS n
C  where n is from 1 (use for validation runs) to the number available.
C  Tests show the optimum is about 12 - probably due to i/o issues.
C
      subroutine sussix4drivenoise(xy,tunexy,amplitude,phase,ox,ax,oy,ay
     &, path)
C
C Parallelised version of sussix4drivexxNoO.f of 15/05/2011 using Openmp.
C Has matching Drive_God_lin.c      H.Renshall & E.Maclean
C
C  xy has been critically defined of dimension 40004, x1, y1 , x2, y2 max 10000
C  and 4 for the window. See matching maxturns parameter in main and datspe.
C 
C      program sussixv4
C=======================================================================
C
C SUSSIX CONTAINS:
C
C ROUTINES FOR HIGH PRECISION TUNE CALCULATION:
C
C         TUNELASR (NO WINDOW)
C         TUNENEWT (HANNING WINDOW)
C         BOTH USING ZFUNR E CALCR
C
C ROUTINES FOR FREQUENCY ANALYSIS OF SIGNALS:
C
C         SPECTRUM      !SPECTRUM COMPUTATION
C
C         ORDRES        !ORDERING OF FREQUENCIES
C         READRES       !READING THE ORDERED FREQUENCIES
C
C ROUTINES FOR POST PROCESSING OF SIGNALS:
C
C         SUSRES        !SUBTRACTION OF NEXT TO LEADING FREQUENCIES
C         READSME       !SMEAR AS QUALITY FACTOR
C         READINV       !INVARIANT CALCULATION
C         READDT3       !THIRD ORDER CONJUGATING FUNCTION COEFFICIENTS
C         READDT4       !FOURTH ORDER CONJUGATING FUNCTION COEFFICIENTS
C
C ROUTINES INTERFACE WITH SIXTRACK OUTPUT
C
C         READRIC       !READS SIXTRACK TRACKING OUTPUT
C         WRITERIC      !WRITES POSTPROCESSED DATA IN SIXTRACK OUTPUT
C                        FORMAT
C
C AUTHOR: R.BARTOLINI CERN AND BOLOGNA UNIVERSITY
C LAST MODIFIED:08/10/1997
C
C=======================================================================
C
C  MAIN PROGRAM FOR THE SPECTRAL ANALYSIS OF TRACKING OR BEAM DATA
C
C  TWO TYPES OF DATA MAY BE ANALYSED ACCORDING TO THE ISIX OPTION:
C
C  ISIX=0 ---> TRACKING DATA FROM A USER PROVIDED ASCII FILE
C  ISIX=1 ---> SIXTRACK BINARY OUTPUT,
C
C  IF ISIX=1 IT TRANSFORMS THE SIXTRACK BINARY OUTPUT INTO AN ASCII
C  FILE (SUBROUTINE READRIC),
C
C  MORE THAN ONE FILE CAN BE TREATED:
C  NTOT IS THE TOTAL NUMBER OF FILES TO BE PROCESSED
C  THE ASCII DATA FILES MUST START WITH FORT.90 IN DESCENDING ORDER
C  AND IT MUST HAVE TWO COLUMNS FOR EACH PLANE!!!!!!!!!!!!!!!!!!!!!
C
C  ONCE THE DATA ARE READ THE SPECTRUM IS CALCULATED, ORDERED AND
C  THE OUTPUT OF ALL THE CASES IS WRITTEN IN THE FILE FORT.300
C  (SUBROUTINE DATSPE ---> SUBROUTINE RSPECTRUM, ORDRES, ETC.)
C
C  DIFFERENT KIND OF POSTPROCESSING ARE AVALIABLE AN ARE SWITCHED
C  ON WITH THE CORRESPONDING FLAG, BY READING THE FORT.300 WITH
C  THE SUBROUTINE  READRES.
C  NSUS GE 1 ---> SUBTRACTS THE NEXT TO LEADING FREQUENCIES
C  (SUBROUTINE SUSRES). IN THIS CASE THE PROGRAM ALSO WRITES THE MODIFIED
C  DATA BACK INTO THE STARTING SIXTRACK BYNARY OUTPUT
C  (SUBROUTINE WRITERIC). WARNING: IT OVERWRITES THE FILES.
C  ISME=1 ---> SMEAR CALCULATION
C  (SUBROUTINE READSME)
C  INV=1 ---> INVARIANT CALCULATION
C  (SUBROUTINE READINV)
C
C  N.B.: IF THE OUTPUT IN FORT.300 IS ALREADY AVALIABLE THE OPTION
C  IANA=0 ALLOWS TO SKIP DATA ANALYSIS, STARTING DIRECTLY WITH THE
C  FORT.300 ANALYSIS.
C
C  DESCRIPTION OF ALL INPUT ITEMS:
C
C  ISIX FLAG FOR SIXTRACK (1) OR ASCII (0) DATA
C  NTOT TOTAL NUMBER OF DATA FILES TO BE ANALYZED
C  IANA FLAG FOR FULL ANALYSIS OF DATA
C  ICONV IS THE FLAG FOR LINEAR TRASFORMATION (YES=1)
C  NT1,NT2 initial & final turn number
C  NARM THE MUNBER OF HARMONIC TO BE CALCULATED
C  ISTUNE IS THE FLAG FOR FUNDAMENTALFREQUENCIES
C  ISTUNE=0 => off
C  ISTUNE=1 uses Qx,y,z values as guess values
C  ISTUNE=2 takes the Qx,y,z values as fundamental frequencies
C  etune(3) allowed distance to Qx,y,z
C  TUNEX,TUNEY,TUNEZ guess or fundamental tunes
C  NSUS THE NUMBER OF HARMONIC TO BE SUBTRACTED
C  IDAM IS THE DIMENSION OF PHASE SPACE
C  NTWIX IS A FLAG FOR THE TWIN PARTICLES (1 or 2)
C  IR IS A FLAG FOR REAL SIGNAL (YES=1)
C  IMETH IS A FLAG FOR WINDOWING (HANNING=1, NO FILTER=0)
C  NRC IS THE MAXIMUM ORDER OF LINEAR COMBINATION OF FREQUENCIES
C  EPS IS THE TOLERANCE ON THE IDENTIFICATION OF FREQUENCIES
C  NLINE IS THE NUMBER OF LINES TO BE LOOKED FOR
C  LR MR KR SPECIFY THE LINE
C  IDAMX SELECT THE PLANE TO ANALYZE
C  IFIN IS THE UNIT FOR THE FINAL OUTPUT
C  ISME FLAG FOR SMEAR CALCULATION
C  IUSME UNIT FOR SMEAR OUTPUT
C  INV FLAG FOR INVARIANTS CALCULATION
C  IINV UNIT FOR INVARIANTS OUTPUT
C  ICF FLAG FOR INVARIANTS CALCULATION
C  IICF UNIT FOR INVARIANTS OUTPUT
C
C  THE SPECTRUM IS WRITTEN IN THE UNIT 300
C
C  AUTHORS: R.BARTOLINI & F.SCHMIDT
C  MODIFIED 16/09/1996: ADDED THE ICONV OPTION
C  MODIFIED 17/09/1996: ADDED THE INPUT FROM FILE
C  MODIFIED 20/08/1997: ADDED THE INVARIANT CALCULATION
C  MODIFIED 17/09/1997: ADDED THE TREATEMENT OF ASCII BEAM DATA
C  MODIFIED 30/05/1998: ADDED DEFAULT VALUES, AND FFT
C  MODIFIED 12/10/1999: TOTAL CLEAN-UP
C
      implicit none
      integer i,iana,icf,iconv,idam,idamx,ifi,ifin,iicf,iinv,imeth,ini,
     &inv,iouk,ir,isix,isme,istune,iunit,iusme,k,kr,lr,mr,mterm,n,narm,
     &nf,nline,nlst,nrc,nsus,nt1,nt2,ntot,ntotal,nturn,ntwin,ntwix,
     &maxturns
      parameter (maxturns=10000)
      double precision eps,etune,tunex,tuney,tunez, xy(maxturns*4+4),
     &tunexy(2),amplitude(19), phase(19), ox(300), ax(300), oy(300), 
     &ay(300)
      parameter(mterm=300)
      dimension lr(100),mr(100),kr(100),etune(3)
      character*200 ch,ch1
      character*8 filename
      character*4000 path
      integer tid,itid,OMP_GET_THREAD_NUM
      save
!$OMP THREADPRIVATE(tid,itid,
!$OMP&i,iana,icf,iconv,idam,idamx,ifi,ifin,iicf,iinv,imeth,ini,
!$OMP&inv,iouk,ir,isix,isme,istune,iunit,iusme,k,kr,lr,mr,n,narm,
!$OMP&nf,nline,nlst,nrc,nsus,nt1,nt2,ntot,ntotal,nturn,ntwin,ntwix,
!$OMP&eps,etune,tunex,tuney,tunez,
!$OMP&ch,ch1,filename)
C===================
C.....INIZIALIZATION
C===================
CHRR      call omp_set_dynamic(.false.)
C=========================
C.....DEFAULT VALUES FIRST
C=========================
      isix=1
      ntot=1
      iana=1
      iconv=1
      nt1=1
      nt2=1024
      narm=1
      istune=0
      etune(1)=1d-2
      etune(2)=1d-2
      etune(3)=1d-2
      tunex=0.28
      tuney=0.31
      tunez=0.006
      nsus=0
      idam=2
      ntwix=1
      ir=0
      imeth=2
      nrc=10
      eps=1d-6
      nline=1
      do 1 i=1,100
        lr(i)=0
        mr(i)=0
        kr(i)=0
 1    continue
      lr(1)=1
      idamx=1
      ifin=500
      isme=0
      iusme=200
      inv=0
      iinv=250
      icf=0
      iicf=350
C=========================
C.....READS THE INPUT FILE
C=========================
c     write(*,*)"Here comes the input file"
      i=index(path, 'sussix_v4.inp')
c     write(*,*) path(1:i+12), i
      TID = OMP_GET_THREAD_NUM()
c     write(*,*) tid
      itid= tid + 10
!$OMP CRITICAL
      open(itid,file = path(1:i+12), form='formatted',
     &status='unknown')
      read(itid,'(4(/))')
      read(itid,'(A)') ch
      ch1=ch(9:80)//' / '
      read(ch1,*)isix
      read(itid,'(A)') ch
      ch1=ch(9:80)//' / '
      read(ch1,*)ntot
c     write(*,*)"ntot=",ntot
      read(itid,'(A)') ch
      ch1=ch(9:80)//' / '
      read(ch1,*)iana
      read(itid,'(A)') ch
      ch1=ch(9:80)//' / '
      read(ch1,*)iconv
      read(itid,'(A)') ch
      ch1=ch(9:80)//' / '
      read(ch1,*)nt1,nt2
      read(itid,'(A)') ch
      ch1=ch(9:80)//' / '
      read(ch1,*)narm
      read(itid,'(A)') ch
      ch1=ch(9:80)//' / '
      read(ch1,*)istune,etune(1),etune(2),etune(3)
      read(itid,'(A)') ch
      ch1=ch(9:80)//' / '
      read(ch1,*)tunex,tuney,tunez
      read(itid,'(A)') ch
      ch1=ch(9:80)//' / '
      read(ch1,*)nsus
      read(itid,'(A)') ch
      ch1=ch(9:80)//' / '
      read(ch1,*)idam
      read(itid,'(A)') ch
      ch1=ch(9:80)//' / '
      read(ch1,*)ntwix
      read(itid,'(A)') ch
      ch1=ch(9:80)//' / '
      read(ch1,*)ir
      read(itid,'(A)') ch
      ch1=ch(9:80)//' / '
      read(ch1,*)imeth
      read(itid,'(A)') ch
      ch1=ch(9:80)//' / '
      read(ch1,*)nrc
      read(itid,'(A)') ch
      ch1=ch(9:80)//' / '
      read(ch1,*)eps
      read(itid,'(A)') ch
      ch1=ch(9:80)//' / '
      read(ch1,*)nline
      do k=1,nline
        read(itid,'(A)') ch
        ch1=ch(9:80)//' / '
        read(ch1,*) lr(k),mr(k),kr(k)
CHRR        write (6,*)"k,lr,mr,kr",k,lr(k),mr(k),kr(k)
      enddo
      read(itid,'(A)') ch
      ch1=ch(9:80)//' / '
      read(ch1,*)idamx
      read(itid,'(A)') ch
      ch1=ch(9:80)//' / '
      read(ch1,*)ifin
      read(itid,'(A)') ch
      ch1=ch(9:80)//' / '
      read(ch1,*)isme
CHRR      write (6,*)"ISME",isme
      read(itid,'(A)') ch
      ch1=ch(9:80)//' / '
      read(ch1,*)iusme
      read(itid,'(A)') ch
      ch1=ch(9:80)//' / '
      read(ch1,*)inv
      read(itid,'(A)') ch
      ch1=ch(9:80)//' / '
      read(ch1,*)iinv
      read(itid,'(A)') ch
      ch1=ch(9:80)//' / '
      read(ch1,*)icf
      read(itid,'(A)') ch
      ch1=ch(9:80)//' / '
      read(ch1,*)iicf
      close(itid)
!$OMP END CRITICAL
C.....CHECKS
      if(narm.le.0) then
C        write(6,*) 'NARM too small'
        close(itid)
        stop
      endif
      if(narm.gt.mterm) then
C        write(6,*) 'NARM too big => reduced to maximum: ',mterm
        narm=mterm
      endif
      if(nt1.le.0) then
C        write(6,*) 'NT1 too small'
        close(itid)
        stop
      endif
      if(nt2.le.nt1) then
C        write(6,*) 'NT2 smaller than NT1'
        close(itid)
        stop
      endif
      if(idam.ne.1.and.idam.ne.2.and.idam.ne.3) then
C        write(6,*) 'The order of phase space IDAM must be 1, 2 or 3'
        close(itid)
        stop
      endif
      if(ntot.gt.32) then
C        write(6,*) 'NTOT too big => reduced to maximum: ',32
        ntot=32
      endif
      if(ntwix.ne.1.and.ntwix.ne.2) then
C        write(6,*) 'NTWIX ill defined, set to: ',1
        ntwix=1
      endif
 
C==============================
C.....END OF THE INIZIALIZATION
C==============================
C===========================
C.....STARTING DATA ANALYSIS
C===========================
 
C.....1) CHECK IF THE FORT.300 IS ALREADY PRODUCED
C.....IF NOT (IANA=1 or 2)
C........2) CHECK THE ISIX OPTION
C........3) PROCEED WITH SPECTRUM CALCULATION
C........4) FFT IF IANA=2
C.....IF YES (IANA=0)
C........2) PROCEED WITH POSTPROCESSING
 
C      open(30,file='fort.300',form='formatted',status='unknown')
      if(iana.eq.1.or.iana.eq.2) then
        if(isix.eq.1) then
C@@@@@@@@@@@@@@@@@@@@@@
C.....SIXTRACK DATA   @
C@@@@@@@@@@@@@@@@@@@@@@
!$OMP CRITICAL
          nlst=90-ntot+1
          open(91,file='fort.91',form='formatted',status='unknown')
          open(92,file='fort.92',form='formatted',status='unknown')
          do n=90,nlst,-1
            if(n.gt.58) goto 100
            write(6,*) 'Unit for sixtrack input must be above 58'
            stop
 100  continue
            filename='fort.'
            write(filename(6:7),'(i2.2)') n
            open(n,file=filename,form='unformatted',status='unknown')
            call readric(n,idam,ntwin,iconv)
C.....NOW READRIC HAS CREATED NTWIN ASCII FILES IN UNIT 91 AND
C.....EVENTUALLY 92. NTWIN IS AN OUTPUT PARAMETER OF READRIC!!
            if(ntwix.lt.ntwin) then
C              write(6,*)'WARNING: THE TWIN PARTICLE IS IGNORED'
              if(nsus.ge.2) then
C                write(6,*)'ERROR: TWIN PARTICLE NEEDED'
                close(10)
C                close(30)
                stop
              endif
            endif
            do nf=1,ntwix
              iunit=90+nf
              call datspe(xy,eps,iunit,idam,ir,nt1,nt2,nturn,
     &             imeth,narm,nrc,iana)
C.....N.B. NTURN IS AN OUTPUT PARAMETER OF DATSPE
              call ordres(eps,narm,nrc,ir,idam,iunit,nturn,
     &             -tunex,-tuney,-tunez,istune,etune, 
     &             tunexy, amplitude, phase, ox, ax, oy,ay)
              if(nsus.ge.2) then
C.....SUBTRACTS AND OVERWRITE UNIT IUNIT
                call susres(iunit,nsus,nturn,3)
              endif
            enddo
 
            if(nsus.ge.2) then
C.....N.B. WRITERIC NEEDS BOTH THE FILES TREATED WITH SUSRES
C.....     IT OVERWRITES THE INITIAL SIXTRACK OUTPUT!!!!!!!!
              call writeric(n,ntotal,ntwin,nturn,iconv)
            endif
          enddo
!$OMP END CRITICAL
        elseif(isix.eq.0) then
C@@@@@@@@@@@@@@@@@@@@
C......ASCII DATA   @
C@@@@@@@@@@@@@@@@@@@@
          do n=1,ntot
            iunit=91-n
            filename='fort.'
            write(filename(6:7),'(i2.2)') iunit
C            open(iunit,file=filename,form='formatted',status='unknown')
            call datspe(xy, eps,iunit,idam,ir,nt1,nt2,nturn,
     &           imeth,narm,nrc,iana)
            call ordres(eps,narm,nrc,ir,idam,iunit,nturn,
     &             -tunex,-tuney,-tunez,istune,etune, 
     &           tunexy, amplitude, phase, ox, ax, oy,ay)
C            close(iunit)
          enddo
        endif
      endif
      close(91)
      close(92)
C=============================================
C.....STARTING POSTPROCESSING OF ANALYZED DATA
C=============================================
 
C......1) SELECTION OF LINES
C......2) SMEAR CALCULATION
C......3) INVARIANT CALCULATION
 
      if(isix.eq.1) then
C.......NTOT*NTWIX CASES ANALYZED FROM FORT.300
        ini=90
        ifi=ini-ntot*ntwix+1
      elseif(isix.eq.0) then
C.......NTOT CASES ANALYZED FROM FORT.300
        ini=90
        ifi=ini-ntot+1
      endif
 
C      write(6,*)' '
C      write(6,*)'*************************************** '
C      write(6,*)'STARTING THE FORT.300 ANALYSIS OF CASES:',ini,ifi
C      write(6,*)'*************************************** '
 
C.....SELECTION OF NLINE LINES SPECIFIED IN THE ARRAY LR,MR,KR
C.....THE RESULTS ARE WRITTEN IN THE UNIT IOUK, IDAMX IS THE
C.....PLANE TO BE ANALYZED
      do k=1,nline
        iouk=ifin+k
        if(iouk.gt.999) then
C          write(6,*) 'Unit for resonance output  must be below 1000'
          close(10)
C          close(30)
          stop
        endif
        filename='fort.'
        write(filename(6:8),'(i3.3)') iouk
        open(iouk,file=filename,form='formatted',status='unknown')
        call readres(ini,ifi,lr(k),mr(k),kr(k),narm,
     &                  idam,idamx,iouk)
        close(iouk)
      enddo
C.....SMEAR CALCULATION: OUTPUT IN THE UNIT IUSME
      if(isme.eq.1) then
         write(6,*)'SMEAR CALCULATION: OUTPUT IN THE UNIT',iusme
        if(iusme.gt.999) then
C          write(6,*) 'Unit for smear calculation must be below 1000'
          close(10)
C          close(30)
          stop
        endif
        filename='fort.'
        write(filename(6:8),'(i3.3)') iusme
        open(iusme,file=filename,form='formatted',status='unknown')
        call readsme(ini,ifi,narm,idam,idamx,iusme)
        close(iusme)
      endif
C.....INVARIANT CALCULATION: OUTPUT IN THE UNIT IINV
      if(inv.eq.1) then
        if(iinv.gt.999) then
C          write(6,*) 'Unit for invariant calculation must be below 1000'
          close(10)
C          close(30)
          stop
        endif
        filename='fort.'
        write(filename(6:8),'(i3.3)') iinv
        open(iinv,file=filename,form='formatted',status='unknown')
C        write(6,*)'INVARIANT CALCULATION: OUTPUT IN THE UNIT',iinv
        call readinv(ini,ifi,narm,idam,idamx,iinv)
        close(iinv)
      endif
C.....3-RD ORDER CONJUGATING FUNCTION CALCULATION: OUTPUT IN THE UNIT ICF
      if(icf.eq.1) then
C        write(6,*)'3-RD ORDER CONJUG. FUNC.: OUTPUT IN THE UNIT',iicf
        if(iicf.gt.999) then
C          write(6,*) 'Unit for conjugating function must be below 1000'
          close(10)
C          close(30)
          stop
        endif
        filename='fort.'
        write(filename(6:8),'(i3.3)') iicf
        open(iicf,file=filename,form='formatted',status='unknown')
        call readdt3(ini,ifi,narm,idam,idamx,iicf)
        call readdt4(ini,ifi,narm,idam,idamx,iicf)
        close(iicf)
      endif
c      close(10)
C      close(30)
      return
      end
      subroutine spectrum(x,xp,maxn,tune,zpesi,narm,meth)
C=======================================================================
C
C SUBROUTINE SPECTRUM
C
C COMPUTE THE MAIN FREQUENCY
C X, XP ARE THE COORDINATES OF THE ORBIT AND MAXN IS THE
C LENGTH IF THE ORBIT.
C WITHOUT ORTHOGONALIZATION OF GRAM-SCHMIDT
C
C METH SELECTS THE WINDOW:
C     1 --> HANNING WINDOW
C     2 --> RECTANGULAR WINDOW
C
C AUTHOR:  R. BARTOLINI 9/1/1996
C          A. BAZZANI
C
C=======================================================================
      implicit none
      integer maxiter,maxn,meth,mterm,n,na,narm
      external tunelasr,tunenewt
      double precision duepi,freq,pi,tune,tunelasr,tunenewt,x,xp
      double complex z,zef,zgs,zpesi,zw,zx,zz
      parameter(mterm=300)
      parameter(maxiter=100000)
      dimension x(maxiter),xp(maxiter)
      dimension z(maxiter),zz(maxiter)
      dimension tune(mterm),zpesi(mterm),zgs(maxiter)
      save
!$OMP THREADPRIVATE(z,zz,freq,zef,zw,zx,zgs,n,na,duepi,pi)
C===============
C INIZIALIZATION
C===============
 
      pi=atan(1d0)*4d0
      duepi=2*pi
      if(maxn.gt.maxiter) then
C        write(6,*) 'ERROR IN SPECTRUM: MAXN TOO LARGE'
        close(10)
C        close(30)
        stop
      endif
      if(narm.lt.2) then
C        write(6,*) 'ERROR IN SPECTRUM: NA SMALLER THAN 2'
        close(10)
C        close(30)
        stop
      endif
 
      do n=1,maxn
        z(n)=dcmplx(x(n),xp(n))
        zz(n)=z(n)
      enddo
 
      do na=1,narm
        if(meth.eq.1) then
          tune(na)=tunenewt(x,xp,maxn,zw)
        elseif(meth.eq.2) then
          tune(na)=tunelasr(x,xp,maxn,zw)
        endif
C.....beginning of subtraction procedure
        freq=tune(na)
        zpesi(na)=zw/dble(maxn)
        zef=exp(dcmplx(0d+0,freq*duepi))
        zx=1
        zgs(1)=zpesi(na)*zx
 
        do n=2,maxn
          zx=zx*zef
          zgs(n)=zpesi(na)*zx
        enddo
 
        do n=1,maxn
          z(n)=z(n)-zgs(n)
        enddo
        do n=1,maxn
         x(n)=dble(z(n))
         xp(n)=dimag(z(n))
        enddo
 
      enddo
 
C.....restore the original signal in x, xp.
      do n=1,maxn
        x(n)=dble(zz(n))
        xp(n)=dimag(zz(n))
      enddo
 
      return
      end
      double precision function tunenewt(x,xp,maxn,zw)
C=======================================================================
C
C SUBROUTINE TUNENEWT
C
C COMPUTES THE TUNE USING A DISCRETE VERSION OF LASKAR METHOD.
C IT INCLUDES A NEWTON METHOD FOR THE SEARCH OF THE FREQUENCY.
C X, XP ARE THE COORDINATES OF THE ORBIT AND MAXN IS THE
C LENGTH IF THE ORBIT.
C
C AUTHOR:     A. BAZZANI - BOLOGNA UNIVERSITY
C             R. BARTOLINI - CERN HAS INTRODUCED SOME MODIFICATIONS
C
C=======================================================================
      implicit none
      integer maxiter,maxn,maxn2,mf,mft,nft,nftmax,npoint
      double precision deltat,duepi,ftmax,step,tune,tune1,tunefou,x,xp
      double complex z,zw
      complex zsing
      parameter(maxiter=100000)
      dimension x(maxiter),xp(maxiter),zsing(maxiter)
      dimension z(maxiter)
      save
!$OMP THREADPRIVATE(z,zsing,maxn,maxn2,mf,mft,nft,nftmax,npoint,
!$OMP&              deltat,ftmax,step,tune,tune1,tunefou,duepi)
C.............................................................
C    ESTIMATION OF TUNE WITH FFT
C.............................................................
      duepi=atan(1d0)*8d0
      mft=int(log(dble(maxn))/log(2d0))
      npoint=2**mft
      maxn2=maxn/2
      step=duepi/maxn
      do mf=1,maxn
        z(mf)=dcmplx(x(mf),xp(mf))*(1d0+cos(step*(mf-maxn2)))
        zsing(mf)=z(mf)
      enddo
      call cfft(zsing,-mft)
C.............................................................
C   SEARCH FOR MAXIMUM OF FOURIER SPECTRUM
C.............................................................
      ftmax=0d0
      nftmax=0
      do nft=1,npoint
        if (abs(zsing(nft)).gt.ftmax) then
          ftmax=abs(zsing(nft))
          nftmax=nft
        end if
      enddo
      tunefou=dble(nftmax-1)/dble(npoint)
      if(tunefou.ge.0.5d+0) tunefou=-(1d+0-tunefou)
      deltat=1d0/npoint
      tune1=tunefou-deltat
      call zfunr(tune,zw,z,maxn,tune1,deltat)
      tunenewt=tune
 
C............................................................
      return
C............................................................
      end
      double precision function tunelasr(x,xp,maxn,zw)
C=======================================================================
C
C SAME AS TUNENEWT BUT NO HANNING FILTER
C
C AUTHOR:     A. BAZZANI - BOLOGNA UNIVERSITY
C             R. BARTOLINI - CERN HAS INTRODUCED SOME MODIFICATIONS
C
C=======================================================================
      implicit none
      integer maxiter,maxn,maxn2,mf,mft,nft,nftmax,npoint
      double precision deltat,duepi,ftmax,step,tune,tune1,tunefou,x,xp
      double complex z,zw
      complex zsing
      parameter(maxiter=100000)
      dimension x(maxiter),xp(maxiter),zsing(maxiter)
      dimension z(maxiter)
      save
!$OMP THREADPRIVATE(z,zsing,maxn,maxn2,mf,mft,nft,nftmax,npoint,
!$OMP&              deltat,ftmax,step,tune,tune1,tunefou,duepi)
C.............................................................
C    ESTIMATION OF TUNE WITH FFT
C.............................................................
      duepi=atan(1d0)*8d0
      mft=int(log(dble(maxn))/log(2d0))
      npoint=2**mft
      maxn2=maxn/2
      step=duepi/maxn
      do mf=1,maxn
        z(mf)=dcmplx(x(mf),xp(mf))  ! no filter
        zsing(mf)=z(mf)
      enddo
      call cfft(zsing,-mft)
C.............................................................
C   SEARCH FOR MAXIMUM OF FOURIER SPECTRUM
C.............................................................
      ftmax=0d0
      nftmax=0
      do nft=1,npoint
        if (abs(zsing(nft)).gt.ftmax) then
          ftmax=abs(zsing(nft))
          nftmax=nft
        end if
      enddo
      tunefou=dble(nftmax-1)/dble(npoint)
      if(tunefou.ge.0.5d+0) tunefou=-(1d+0-tunefou)
      deltat=1d0/npoint
      tune1=tunefou-deltat
      call zfunr(tune,zw,z,maxn,tune1,deltat)
      tunelasr=tune
 
      return
      end
      subroutine zfunr(tune,zw,z,maxn,tunea1,deltat)
C=======================================================================
C AUXILIARY ROUTINE USED BY TUNENEWT.
C
C AUTHOR:     A. BAZZANI - BOLOGNA UNIVERSITY
C
C=======================================================================
      implicit none
      integer maxiter,maxn,nc,ncont,nd,ntest,num
      double precision deltat,dtune1,dtune2,dtune3,dtunea1,dtunea2,
     &duepi,err,ratio,tune,tune1,tune2,tune3,tunea1,tunea2,tunetest,
     &tuneval,tunevmax
      double complex z,zd,zf,zfd,ztune,ztune1,ztune2,ztune3,zu,zw
      parameter(maxiter=100000)
      dimension z(*),zd(maxiter),tunetest(10),tuneval(10)
      save
!$OMP THREADPRIVATE(zd,zf,zfd,ztune,ztune1,ztune2,ztune3,zu,
!$OMP&              deltat,dtune1,dtune2,dtune3,dtunea1,dtunea2,
!$OMP&              err,ratio,tune1,tune2,tune3,tunea2,tunetest,
!$OMP&              tuneval,tunevmax,nc,ncont,nd,ntest,num,duepi)
C............................................................
      duepi=atan(1d0)*8d0
      err=1d-10
      zu=dcmplx(0d0,1d0)
C............................................................
C.... WE DIVIDE DELTAT IN 5 PARTS
C............................................................
      deltat=deltat/5.d0
C............................................................
      do nd=1,maxn
        zd(nd)=zu*nd*z(nd)
      enddo
C............................................................
      ztune1=cdexp(-zu*duepi*tunea1)
! Put calcr calls inline                                                 HRR
!      call calcr(ztune1,zf,z,maxn)
!      call calcr(ztune1,zfd,zd,maxn)
         zf= z(maxn)
         zfd= zd(maxn)
         do nd= maxn-1,1,-1
            zf= zf*ztune1 + z(nd)
            zfd= zfd*ztune1 + zd(nd)
         enddo
      dtunea1=dble(zf)*dble(zfd)+dimag(zf)*dimag(zfd)
      num=1
      do ntest=1, 10
        tunea2=tunea1+deltat
        ztune2=cdexp(-zu*duepi*tunea2)
! Put calcr calls inline                                                 HRR
!       call calcr(ztune2,zf,z,maxn)
!       call calcr(ztune2,zfd,zd,maxn)
         zf= z(maxn)
         zfd= zd(maxn)
         do nd= maxn-1,1,-1
            zf= zf*ztune2 + z(nd)
            zfd= zfd*ztune2 + zd(nd)
         enddo
        dtunea2=dble(zf)*dble(zfd)+dimag(zf)*dimag(zfd)
        if ((dtunea1.le.0d0).and.(dtunea2.ge.0d0)) then
           tune1=tunea1
           tune2=tunea2
           dtune1=dtunea1
           dtune2=dtunea2
           do ncont=1,100
c ......Insertion, to avoid /0.0
             if(dtune2.ne.0) then
               ratio=-dtune1/dtune2
             else 
               ratio=0.0
             endif
              tune3=(tune1+ratio*tune2)/(1.d0+ratio)
              ztune3=cdexp(-zu*duepi*tune3)
! Put calcr calls inline                                                 HRR
!            call calcr(ztune3,zf,z,maxn)
!            call calcr(ztune3,zfd,zd,maxn)
! Execution is slightly faster with two separate loops under -O3
             zf= z(maxn)
             zfd= zd(maxn)
             do nd= maxn-1,1,-1
                zf= zf*ztune3 + z(nd)
                zfd= zfd*ztune3 + zd(nd)
             enddo
              dtune3=dble(zf)*dble(zfd)+dimag(zf)*dimag(zfd)
              if (dtune3.le.0d0) then
                 if(tune1.eq.tune3) goto 100
                 tune1=tune3
                 dtune1=dtune3
              else
                 if(tune2.eq.tune3) goto 100
                 tune2=tune3
                 dtune2=dtune3
              endif
              if (abs(tune2-tune1).le.err) goto 100
           enddo
100        tunetest(num)=tune3
           tuneval(num)=cdabs(zf)
           num=num+1
        endif
        tunea1=tunea2
        dtunea1=dtunea2
      enddo
      tune=tunetest(1)
      tunevmax=tuneval(1)
      do nc=2, num-1
         if(tunevmax.le.tuneval(nc)) then
            tunevmax=tuneval(nc)
            tune=tunetest(nc)
         endif
      enddo
      ztune=cdexp(-zu*duepi*tune)
      call calcr(ztune,zw,z,maxn)
C............................................................
      return
C............................................................
      end
      subroutine calcr(zv,zpp,zp,maxd)
C=======================================================================
C AUXILIARY ROUTINE USED BY TUNENEWT.
C
C AUTHOR:     A. BAZZANI - BOLOGNA UNIVERSITY
C
C=======================================================================
      implicit none
      integer maxd,np
      double complex zp,zpp,zv
      dimension zp(*)
      save
!$OMP THREADPRIVATE(np)
      zpp=zp(maxd)
C............................................................
      do np=maxd-1,1, -1
        zpp=zpp*zv+zp(np)
      enddo
C............................................................
      return
C............................................................
      end
      subroutine susres(nfile,nsus,max,idams)
C=======================================================================
C
C SUBROUTINE SUSRES
C
C SUBTRACTS TO THE SIGNALS THE NEXT TO LEADING FREQUENCIES
C LEAVING THE TUNE IN.
C
C NSUS = HARMONICS TO BE SUBTRACTED
C MAX  = LENGTH OF THE SIGNAL
C NFILE = THE OUTPUT FILE WITH THE SUBTRACTED SIGNAL
C
C THIS ROUTINE MUST BE CALLED ONLY AFTER THE EXECUTION OF
C THE ROUTINE DATSPE WHICH FILLS THE COMMONS
C
C AUTHOR: R.BARTOLINI
C LAST MODIFIED: 01/04/1996
C
C=======================================================================
      implicit none
      integer idams,j,k,max,maxn,mterm,nfile,nsus
      double precision duepi,ss,ssp,xx,xxp,yy,yyp
      double complex zpots,zpotx,zpoty,zts,ztx,zty,zxs,zys,zss
      parameter(mterm=300)
      parameter(maxn=100000)
      dimension zxs(maxn),zys(maxn),zss(maxn)
      double precision s,sp,x,xp,y,yp
      common/data/x(maxn),y(maxn),xp(maxn),yp(maxn),s(maxn),sp(maxn)
      double precision tsa,txa,tya
      common/tune/txa(mterm),tya(mterm),tsa(mterm)
      double complex zspes,zxpes,zypes
      common/fcoe/zxpes(mterm),zypes(mterm),zspes(mterm)
!$OMP THREADPRIVATE(/data/,/fcoe/,/tune/,
!$OMP&              j,k,max,
!$OMP&              duepi,ss,ssp,xx,xxp,yy,yyp,
!$OMP&              zpots,zpotx,zpoty,zts,ztx,zty,zxs,zys,zss)
      save
      duepi=8*atan(1d+0)
 
C      write(6,*)'STARTING THE SUBTRACTION PROCEDURE IN FILE',nfile
 
C....BUILD THE COMPLEX SIGNAL Z = X + i PX
 
      do j=1,max
        zxs(j)=dcmplx(x(j),xp(j))
        zys(j)=dcmplx(y(j),yp(j))
        zss(j)=dcmplx(s(j),sp(j))
      enddo
 
      do j=2,nsus
        zpotx=1d+0
        ztx=exp(dcmplx(0d+0,txa(j)*duepi))
        zpoty=1d+0
        zty=exp(dcmplx(0d+0,tya(j)*duepi))
        zpots=1d+0
        zts=exp(dcmplx(0d+0,tsa(j)*duepi))
        do k=1,max
          zxs(k)=zxs(k)-zxpes(j)*zpotx
          zys(k)=zys(k)-zypes(j)*zpoty
          zss(k)=zss(k)-zspes(j)*zpots
          zpotx=zpotx*ztx
          zpoty=zpoty*zty
          zpots=zpots*zts
        enddo
      enddo
 
C....WRITES THE REMAINING SIGNAL
 
      do k=1,max
        xx=dble(zxs(k))
        xxp=dimag(zxs(k))
        yy=dble(zys(k))
        yyp=dimag(zys(k))
        ss=dble(zss(k))
        ssp=dimag(zss(k))
        if(idams.eq.1) then
          write(nfile,'(2(G21.15,1X))')xx,xxp
        elseif(idams.eq.2) then
          write(nfile,'(4(G21.15,1X))')xx,xxp,yy,yyp
        elseif(idams.eq.3) then
          write(nfile,'(6(G21.15,1X))')xx,xxp,yy,yyp,ss,ssp
        endif
      enddo
 
      return
 
      end
      subroutine ordres(eps,narm,nr,ir,idam,iunit,nturn,tunex,
     &                  tuney,tunez,istune,etune, tunexy,
     &                  amplitude, phase, ox, ax, oy, ay)
C=======================================================================
C
C  ORDERS THE HARMONICS FOUND BY SPECTRUM
C
C  INPUT PARAMETERS:
C
C  NARM : NUMBER OF HARMONICS TO BE OREDERED
C  NR   : MAXIMUM ORDER OF HARMONIC TO BE LOOKED FOR IN THE LINEAR
C         COMBINATIONS. (NR<10 ALTRIMENTI L'OUTPUT PUO INCASINARSI)
C  IR   : FLAG FOR COMPLEX OR REAL SIGNALS:
C         1 = INPUT SIGNAL IS REAL.
C         N.B. EVEN IF THE SIGNAL IS REAL THE SPECTRUM IS NOT
C              TRANSPOSED IN [0,0.5]
C  EPS  : MAXIMUM ERROR ACCEPTED TO FIND THE LINEAR COMBINATIONS.
C  TUNEX,TUNEZ,TUNES ARE THE EXPECTED TUNES WITHIN 0.01
C
C  MIND :  THE INPUT DATA COME FROM THE COMMON, SO THIS ROUTINE
C          MUST BE USED ONLY AFTER THE CALL TO DATSPE
C          THE ARRAYS TX,TY,TZ ARE USED IN ORDER NOT TO CHANGE
C          THE ARRAYS TXA,TYA,TSA.
C
C  THE OUTPUT IS PLACED IN THE FILE FORT.300 WHICH IS NOT CLOSED AT THE
C  END OF THE SUBROUTINE IN ORDER TO COLLECT THE RESULTS OF DIFFERENT
C  SIGNALS IN ONE FILE ONLY. THE DIFFERENT SIGNALS ARE NUMBERED BY IOU.
C
C  A CHECK IS PERFORMED TO SEE IF THE HARMONICS ARE CLOSER THAN 1/NTURN
C
C  IOU  : INDEX WHICH IDENTIFIES THE CASE ANALIZED
C
C  AUTHOR: R.BARTOLINI 12/02/1996
C  LAST MODIFIED: 21/09/1997
C
C=======================================================================
      implicit none
      integer idam,imiss,imissx,imissy,imissz,ir,isca,iscax,iscay,iscaz,
     &isearx,iseary,isearz,istune,iunit,j,j1,k,k1,l,l1,
     &m,m1,mterm,n,narm,narm2,nr,nt,nturn,ntx,nty,ntz, flagad(19),myint
      double precision az,checkn,dt,dtunex,dtuney,dtunez,dtxz,dty,dtyz,
     &eps,epsx,epsy,epsz,etune,ex,ey,ez,fx,fxt,fy,fyt,fz,fzt,ordc,ordcx,
     &ordcy,ordcz,pi,px,pxt,pxti,pxtr,py,pyt,pyti,pytr,pz,pzt,pzti,pztr,
     &tunex,tuney,tunez,tx,txt,ty,tyt,tz,tzt, tunexy(2), 
     &amplitude(19), phase(19), ox(300), ax(300),
     &oy(300), ay(300)
      double complex zpx,zpy,zpz
      parameter(mterm=300)
      double precision tsa,txa,tya
      common/tune/txa(mterm),tya(mterm),tsa(mterm)
      double complex zspes,zxpes,zypes
      common/fcoe/zxpes(mterm),zypes(mterm),zspes(mterm)
      dimension tx(mterm),ty(mterm),tz(mterm),etune(3)
!$OMP THREADPRIVATE(/fcoe/,/tune/,
!$OMP&imiss,imissx,imissy,imissz,isca,iscax,iscay,iscaz,
!$OMP&isearx,iseary,isearz,j,j1,k,k1,l,l1,m,m1,n,
!$OMP&narm2,nt,ntx,nty,ntz,flagad,myint,
!$OMP&pi,
!$OMP&az,checkn,dt,dtunex,dtuney,dtunez,dtxz,dty,dtyz,
!$OMP&epsx,epsy,epsz,ex,ey,ez,fx,fxt,fy,fyt,fz,fzt,ordc,ordcx,
!$OMP&ordcy,ordcz,px,pxt,pxti,pxtr,py,pyt,pyti,pytr,pz,pzt,pzti,pztr,
!$OMP&tx,txt,ty,tyt,tz,tzt,
!$OMP&zpx,zpy,zpz)
      save 
      pi=4*atan(1d0)
      checkn=1d0/dble(nturn)
      flagad(1)=0
      flagad(2)=0
      flagad(3)=0
      flagad(4)=0
      flagad(5)=0
      flagad(6)=0
      flagad(7)=0
      flagad(8)=0
      flagad(9)=0
      flagad(10)=0
      flagad(11)=0
      flagad(12)=0
      flagad(13)=0
      flagad(14)=0
      flagad(15)=0
      flagad(16)=0
      flagad(17)=0
      flagad(18)=0
      flagad(19)=0
      if(nr.gt.10) then
C        write(6,*)'ERROR IN ORDRES: NR LARGER THAN 10'
        close(10)
C        close(30)
C        close(iunit)
        stop
      endif
 
      do j=1,narm
        tx(j)=txa(j)
        ty(j)=tya(j)
        tz(j)=tsa(j)
      enddo
 
C
C.. TUNES PARAMETERS AND EVENTUAL CHECK FOR THE EXPECTED TUNES
C
 
      if(istune.ge.1) then
C.....CHECK X TUNE
        dtunex=abs(abs(tx(1))-abs(tunex))
        if(dtunex.gt.etune(1).or.(tx(1)*tunex).lt.0) then
          write(6,*)'X TUNE DIFFERENT FROM EXPECTED'
C          write(6,*)-tx(1),-tunex
          ntx=1
          do nt=2,narm
            dtunex=abs(abs(tx(nt))-abs(tunex))
            if(dtunex.le.etune(1).and.(tx(nt)*tunex).gt.0) then
              ntx=nt
              goto 7
            endif
          enddo
 7        if(ntx.gt.1) then
            write(6,*)'EXPECTED TUNE X FOUND AT LINE',ntx
          elseif(ntx.eq.1) then
            write(6,*)'EXPECTED TUNE X NOT FOUND'
C            if(istune.eq.1) write(6,*)'LINE 1 ASSUMED AS TUNE!!!'
          endif
          tx(1)=tx(ntx)
          tx(ntx)=txa(1)
          zpx=zxpes(1)
          zxpes(1)=zxpes(ntx)
          zxpes(ntx)=zpx
        endif
      endif
 
      if(istune.eq.2) then
        txt=tunex
      else
        txt=tx(1)
      endif
      pxt=abs(zxpes(1))
c..... Insertion to avoid /0.0
      if (pxt.gt.0) then
        pxtr=dble(zxpes(1))/pxt
        pxti=dimag(zxpes(1))/pxt
        fxt=atan2(pxti,pxtr)
      elseif(pxt.eq.0) then
        pxtr=0
        pxti=0
        fxt=0
      endif
      fxt=fxt/pi*1.8d+2
      tyt=9999.
      tzt=8888.
      if(idam.ge.2) then
        if(istune.ge.1) then
C.....CHECK Y TUNE
          dtuney=abs(abs(ty(1))-abs(tuney))
          if(dtuney.gt.etune(2).or.(ty(1)*tuney).lt.0) then
C            write(6,*)'Y TUNE DIFFERENT FROM EXPECTED'
C            write(6,*)-ty(1),-tuney
            nty=1
            do nt=2,narm
              dtuney=abs(abs(ty(nt))-abs(tuney))
              if(dtuney.le.etune(2).and.(ty(nt)*tuney).gt.0) then
                nty=nt
                goto 8
              endif
            enddo
 8          if(nty.gt.1) then
C              write(6,*)'EXPECTED TUNE Y FOUND AT LINE',nty
            elseif(nty.eq.1) then
C              write(6,*)'EXPECTED TUNE Y NOT FOUND'
C              if(istune.eq.1) write(6,*)'LINE 1 ASSUMED AS TUNE!!!'
            endif
            ty(1)=ty(nty)
            ty(nty)=tya(1)
            zpy=zypes(1)
            zypes(1)=zypes(nty)
            zypes(nty)=zpy
          endif
        endif
 
        if(istune.eq.2) then
          tyt=tuney
        else
          tyt=ty(1)
        endif
        pyt=abs(zypes(1))
c .... Insertion to avoid /0.0
        if(pyt.gt.0) then
          pytr=dble(zypes(1))/pyt
          pyti=dimag(zypes(1))/pyt
          fyt=atan2(pyti,pytr)
        elseif(pyt.eq.0) then
          pytr=0.0          
          pyti=0.0
          fyt=0.0
        endif
        fyt=fyt/pi*1.8d+2
        dty=abs(tyt-txt)
        if(dty.le.eps) then
          write(6,*)'TUNEX AND TUNEY ARE TOO CLOSE'
           tyt=9999.
        endif
      else if(idam.lt.2.and.istune.eq.2) then
        tyt=tuney
      endif
      if(idam.eq.3) then
        if(istune.ge.1) then
C.....CHECK Z TUNE
          dtunez=abs(abs(tz(1))-abs(tunez))
          if(dtunez.gt.etune(3).or.(tz(1)*tunez).lt.0) then
            write(6,*)'Y TUNE DIFFERENT FROM EXPECTED'
            write(6,*)-tz(1),-tunez
            ntz=1
            do nt=2,narm
              dtunez=abs(abs(tz(nt))-abs(tunez))
              if(dtunez.le.etune(3).and.(tz(nt)*tunez).gt.0) then
                ntz=nt
                goto 9
              endif
            enddo
 9          if(ntz.gt.1) then
C              write(6,*)'EXPECTED TUNE S FOUND AT LINE',ntz
            elseif(ntz.eq.1) then
              write(6,*)'EXPECTED TUNE S NOT FOUND'
C              if(istune.eq.1) write(6,*)'LINE 1 ASSUMED AS TUNE!!!'
            endif
            tz(1)=tz(ntz)
            tz(ntz)=tsa(1)
            zpz=zspes(1)
            zspes(1)=zspes(ntz)
            zspes(ntz)=zpz
          endif
        endif
 
        if(istune.eq.2) then
          tzt=tunez
        else
          tzt=tz(1)
        endif
        pzt=abs(zspes(1))
        pztr=dble(zspes(1))/pzt
        pzti=dimag(zspes(1))/pzt
        fzt=atan2(pzti,pztr)
        fzt=fzt/pi*1.8d+2
        dtxz=abs(tzt-txt)
        if(dtxz.le.eps) then
          write(6,*)'TUNEX AND TUNES ARE TOO CLOSE'
          tzt=8888.
        endif
        dtyz=abs(tzt-tyt)
        if(dtyz.le.eps) then
          write(6,*)'TUNEY AND TUNES ARE TOO CLOSE'
           tzt=8888.
        endif
      else if(idam.lt.3.and.istune.eq.2) then
        tzt=tunez
      endif
 
C
C.. X SIGNAL PROCESSING
C
 
C      write(30,*)
C      write(30,*)'ANALYSIS OF X SIGNAL, CASE:',iunit
C      write(30,*)'TUNE X = ',-txt
C      write(30,*)
C      write(30,'(3a)')'        Line Frequency      Amplitude',
C     &'             Phase       ',
C     &'      Error         mx  my  ms  p'
C      write(30,*)
 
      iscax=0
      imissx=0
      do n=1,narm
            ox(n)=-tx(n)
            ax(n)=abs(zxpes(n))
C.....CHECK WITH PREVIOUS HARMONICS
        do k=1,n-1
          dt=abs(tx(k)-tx(n))
          if(dt.le.checkn) then
            iscax=1
          endif
        enddo
C.....SEARCH LINEAR COMBINATIONS
        isearx=0
        ordcx=3000 ! large number needed here
        do l=-nr,nr
          do m=-nr,nr
            do k=-nr,nr
              do j=-nr,nr
                az=l*txt+m*tyt+k*tzt+j
                ex=abs(az-tx(n))
                if(ex.lt.eps) then
                  isearx=isearx+1
C.....check for the lowest possible order of combination
                  ordc=abs(l)+abs(m)+abs(k)
                  if(ordc.lt.ordcx) then
                    l1=l
                    m1=m
                    k1=k
                    j1=j
                    epsx=ex
                    ordcx=ordc
                  endif
                endif
              enddo
            enddo
          enddo
        enddo

        if(isearx.ge.1) then
          px=abs(zxpes(n))
          if(px.gt.0) then
            pxtr=dble(zxpes(n))/px
            pxti=dimag(zxpes(n))/px
            fx=atan2(pxti,pxtr)
            fx=fx/pi*1.8d+2
c...  Insertion, to return parameters
c            write(*,*) l1, m1, k1, j1, flagad(1), tx(n)
            if(l1.eq.1.and.m1.eq.0.and.k1.eq.0.and.
     &           j1.eq.0.and.flagad(1).eq.0) then
              amplitude(1)=px
              phase(1)=-fx
              tunexy(1)=-tx(n)
              flagad(1)=1
            endif             
            if(l1.eq.-2.and.m1.eq.0.and.k1.eq.0.and.
     &           flagad(2).eq.0) then
              amplitude(2)=px
              phase(2)=-fx
              flagad(2)=1
            endif
            
            if(flagad(3).eq.0.and.l1.eq.0.and.m1.eq.1.and.
     &           k1.eq.0)then
              amplitude(3)=px
              phase(3)=-fx
              flagad(3)=1
            endif
            
            if(flagad(7).eq.0.and.l1.eq.-3.and.m1.eq.0.and.
     &           k1.eq.0)then
              amplitude(7)=px
              phase(7)=-fx
              flagad(7)=1
            endif
            
            if(flagad(8).eq.0.and.l1.eq.-4.and.m1.eq.0.and.
     &           k1.eq.0)then
              amplitude(8)=px
              phase(8)=-fx
              flagad(8)=1
            endif
            
            if(flagad(9).eq.0.and.l1.eq.-5.and.m1.eq.0.and.
     &           k1.eq.0.and.j1.eq.0)then
              amplitude(9)=px
              phase(9)=-fx
              flagad(9)=1
            endif
            
            if(flagad(10).eq.0.and.l1.eq.2.and.m1.eq.0.and.
     &           k1.eq.0)then
              amplitude(10)=px
              phase(10)=-fx
              flagad(10)=1
            endif             
            
            if(flagad(11).eq.0.and.l1.eq.0.and.m1.eq.0.and.
     &           k1.eq.0)then
              amplitude(11)=px
              phase(11)=-fx
              flagad(11)=1
            endif

            if(flagad(13).eq.0.and.l1.eq.0.and.m1.eq.2.and.
     &           k1.eq.0.and.j1.eq.0)then
              amplitude(13)=px
              phase(13)=-fx
              flagad(13)=1
            endif 
	    
	    if(l1.eq.-1.and.m1.eq.-1.and.k1.eq.0.and.
     &           j1.eq.0.and.flagad(15).eq.0) then
              amplitude(15)=px
              phase(15)=-fx
              flagad(15)=1
            endif
	    
	    if(l1.eq.2.and.m1.eq.-2.and.k1.eq.0.and.
     &           j1.eq.0.and.flagad(17).eq.0) then
              amplitude(17)=px
              phase(17)=-fx
              flagad(17)=1
            endif
	    
	    if(l1.eq.0.and.m1.eq.-2.and.k1.eq.0.and.
     &           j1.eq.0.and.flagad(19).eq.0) then
              amplitude(19)=px
              phase(19)=-fx
              flagad(19)=1
            endif
	    



            
c...  End of Insertion	
C            write(30,100)n,-tx(n),px,-fx,epsx,l1,m1,k1,j1
          elseif(px.eq.0) then
C            write(30,100)n,-tx(n),px,0,0
          endif
        elseif(isearx.eq.0) then
          imissx=imissx+1
          px=abs(zxpes(n))
          if(px.gt.0) then
            pxtr=dble(zxpes(n))/px
            pxti=dimag(zxpes(n))/px
            fx=atan2(pxti,pxtr)
            fx=fx/pi*1.8d+2
C            write(30,100)n,-tx(n),px,-fx,eps
          elseif(px.eq.0) then
C           write(30,100)n,-tx(n),px,0,0
          endif
        endif
  10     continue
      enddo

C
C.. Y SIGNAL PROCESSING
C

       
      if(idam.lt.2) then
C        write(30,*)
C        write(30,*) 'NO Y SIGNAL'
        if(istune.eq.2) then
C          write(30,*) 'TUNE Y = ',-tyt
C          write(30,*) 'TUNE S = ',-tzt
          goto 51
        else
C          write(30,*)
C          write(30,*)
          goto 50
        endif
      endif
 
C      write(30,*)
C      write(30,*) 'ANALYSIS OF Y SIGNAL, CASE:',iunit
      if(istune.eq.2) then
C        write(30,*) 'TUNE Y = ',-tyt
      else
C        write(30,*) 'TUNE Y = ',-ty(1)
      endif
C      write(30,*)
C      write(30,'(3a)')'        Line Frequency      Amplitude',
C     &'             Phase       ',
C     &'      Error         mx  my  ms  p'
C     write(30,*)
 
      iscay=0
      imissy=0
      do n=1,narm
C.....CHECK WITH PREVIOUS HARMONICS
        do k=1,n-1
          dt=abs(ty(k)-ty(n))
          if(dt.le.checkn) then
            iscay=1
          endif
        enddo
C.....SEARCH LINEAR COMBINATIONS
        iseary=0
        ordcy=3000 ! large number needed here
        do l=-nr,nr
          do m=-nr,nr
            do k=-nr,nr
              do j=-nr,nr
                az=l*txt+m*tyt+k*tzt+j
                ey=abs(az-ty(n))
                if(ey.lt.eps) then
                  iseary=iseary+1
C.....check for the lowest possible order of combination
                  ordc=abs(l)+abs(m)+abs(k)
                  if(ordc.lt.ordcy) then
                    l1=l
                    m1=m
                    k1=k
                    j1=j
                    epsy=ey
                    ordcy=ordc
                  endif
                endif
              enddo
            enddo
          enddo
        enddo
        if(iseary.ge.1) then
          py=abs(zypes(n))
          if(py.gt.0) then
            pytr=dble(zypes(n))/py
            pyti=dimag(zypes(n))/py
            fy=atan2(pyti,pytr)
            fy=fy/pi*1.8d+2
c...  Insertion, to return parameters
            if(n.eq.1) then
              amplitude(4)=py
              write(*,*)"p3 ", fy, py, ty(n), tunexy(1)
              phase(4)=-fy
              tunexy(2)=-ty(n)
            endif
            if(l1.eq.0.and.m1.eq.1.and.k1.eq.0.and.
     &           j1.eq.0.and.flagad(4).eq.0) then
              amplitude(4)=py
              phase(4)=-fy
              tunexy(2)=-ty(n)
c              write(*,*)"p32 ", fy, py, ty(n)
              flagad(4)=1
            endif             
            if(l1.eq.0.and.m1.eq.-2.and.k1.eq.0.and.
     &           flagad(5).eq.0) then
              amplitude(5)=py
              phase(5)=-fy
              flagad(5)=1
            endif
            
            if(flagad(6).eq.0.and.l1.eq.1.and.m1.eq.0.and.
     &         k1.eq.0.and.j1.eq.0)then
              amplitude(6)=py
              phase(6)=-fy
              flagad(6)=1
            endif
            
            if(flagad(12).eq.0.and.l1.eq.0.and.m1.eq.-3.and.
     &         k1.eq.0.and.j1.eq.0)then
              amplitude(12)=py
              phase(12)=-fy
              flagad(12)=1
            endif        

            if(flagad(14).eq.0.and.l1.eq.-1.and.m1.eq.-1.and.
     &           k1.eq.0.and.j1.eq.0)then
              amplitude(14)=py
              phase(14)=-fy
              flagad(14)=1
            endif 
	    
	     if(l1.eq.-2.and.m1.eq.0.and.k1.eq.0.and.
     &           flagad(16).eq.0) then
              amplitude(16)=py
              phase(16)=-fy
              flagad(16)=1
            endif
	    
	     if(l1.eq.1.and.m1.eq.-1.and.k1.eq.0.and.
     &           flagad(18).eq.0) then
              amplitude(18)=py
              phase(18)=-fy
              flagad(18)=1
            endif
	    
	    
	    

c...  End of Insertion
            oy(n)=-ty(n)
            ay(n)=py
C            write(30,100)n,-ty(n),py,-fy,epsy,l1,m1,k1,j1
          elseif(py.eq.0) then
            oy(n)=-ty(n)
            ay(n)=py
C            write(30,100)n,-ty(n),py,0,0
          endif
        elseif(iseary.eq.0) then
          imissy=imissy+1
          py=abs(zypes(n))
          if(py.gt.0) then
            pytr=dble(zypes(n))/py
            pyti=dimag(zypes(n))/py
            fy=atan2(pyti,pytr)
            fy=fy/pi*1.8d+2
            oy(n)=-ty(n)
            ay(n)=py
C            write(30,100)n,-ty(n),py,-fy,eps
          elseif(py.eq.0) then
            oy(n)=-ty(n)
            ay(n)=py
C            write(30,100)n,-ty(n),py,0,0
          endif
        endif
      enddo
 
      
      do myint=1,19
         if(flagad(myint).eq.0) then
            amplitude(myint)=0
            phase(myint)=0
            if(myint.eq.1) then
               tunexy(1)=0
            endif
            if(myint.eq.4) then
               tunexy(2)=0
               write(*,*)"VERTICAL TUNE NOT IDENTIFIED (in sussix)"
            endif
         endif
      enddo

C
C.. S SIGNAL PROCESSING
C
 
 51   if(idam.lt.3) then
C        write(30,*)
C        write(30,*) 'NO S SIGNAL'
        if(istune.eq.2) then
C          write(30,*) 'TUNE S = ',-tzt
        else
C          write(30,*)
        endif
C        write(30,*)
        goto 50
      endif
 
C      write(30,*)
C      write(30,*) 'ANALYSIS OF S SIGNAL, CASE:',iunit
      if(istune.eq.2) then
C        write(30,*) 'TUNE S = ',-tzt
      else
C        write(30,*) 'TUNE S = ',-tz(1)
      endif
C      write(30,*)
C      write(30,'(3a)')'        Line Frequency      Amplitude',
C     &'             Phase       ',
C     &'      Error         mx  my  ms  p'
C      write(30,*)
 
      iscaz=0
      imissz=0
      do n=1,narm
C.....CHECK WITH PREVIOUS HARMONICS
        do k=1,n-1
          dt=abs(tz(k)-tz(n))
          if(dt.le.checkn) then
            iscaz=1
          endif
        enddo
C.....SEARCH LINEAR COMBINATIONS
        isearz=0
        ordcz=3000 ! large number needed here
        do l=-nr,nr
          do m=-nr,nr
            do k=-nr,nr
              do j=-nr,nr
                az=l*txt+m*tyt+k*tzt+j
                ez=abs(az-tz(n))
                if(ez.lt.eps) then
                  isearz=isearz+1
C.....check for the lowest possible order of combination
                  ordc=abs(l)+abs(m)+abs(k)
                  if(ordc.lt.ordcz) then
                    l1=l
                    m1=m
                    k1=k
                    j1=j
                    epsz=ez
                    ordcz=ordc
                  endif
                endif
              enddo
            enddo
          enddo
        enddo
        if(isearz.ge.1) then
          pz=abs(zspes(n))
          pztr=dble(zspes(n))/pz
          pzti=dimag(zspes(n))/pz
          fz=atan2(pzti,pztr)
          fz=fz/pi*1.8d+2
C          write(30,100)n,-tz(n),pz,-fz,epsz,l1,m1,k1,j1
        elseif(isearz.eq.0) then
          imissz=imissz+1
          pz=abs(zspes(n))
          pztr=dble(zspes(n))/pz
          pzti=dimag(zspes(n))/pz
          fz=atan2(pzti,pztr)
          fz=fz/pi*1.8d+2
C          write(30,100)n,-tz(n),pz,-fz,eps
        endif
30    continue
      enddo
 
 50   continue
C50    write(6,*)'================================================'
C      write(6,*)'HARMONICS NON-IDENTIFIED IN X = ',imissx
C      write(6,*)'HARMONICS NON-IDENTIFIED IN Y = ',imissy
C      write(6,*)'HARMONICS NON-IDENTIFIED IN S = ',imissz
      narm2=narm/2d+0
      imiss=imissx+imissy+imissz
      if(imiss.ge.narm) then
C       write(6,*)'WARNING: CHECK EPS'
      endif
      isca=iscax+iscay+iscaz
      if(isca.ge.1) then
C        write(6,*)'WARNING: TOO CLOSE BY HARMONICS DETECTED'
C        write(6,*)'WARNING: TRY A LARGER NUMBER OF TURNS'
      endif
C      write(6,*)'WROTE IN FORT.300 THE SPECTRAL ANALYSIS OF CASE:',iunit
C      write(6,*)'================================================'
 
100   format(i3,1x,f19.16,2(1x,e18.12),1x,e17.11,:,4(1x,i3))
 
      return
      end
      subroutine datspe(xy,eps,iunit,idam,ir,nt1,nt2,nturn,
     &                   imeth,narm,nrc,iana)
C=======================================================================
C
C  SUBROUTINE DATSPE
C
C  THIS PROGRAM CALCULATES THE SPECTRUM OF A SINGLE FILE
C  OF TRACKING DATA WRITTEN IN A STANDARD ASCII FROM IN 2*IDAM COLUMNS.
C
C  WITH THE DATA SOME INFORMATIONS ARE REQUESTED:
C
C    IUNIT = FORTRAN UNIT OF INPUT DATA
C    IDAM  = DIMENSION OF PHASE SPACE
C    IR    = FLAG FOR REAL SIGNAL (1=REAL)
C    NT1   = INITIAL TURN TO BE ANALYZED
C    NT2   = FINAL TURN TO BE ANALYZED
C
C  THE TUNE AND THE LINES ARE CALCULATED WITH THE ROUTINE SPECTRUM
C
C    NARM : THE NUMBER OF HARMONIC TO BE CALCULATED
C    IMETH: THE CHIOCE ON THE WINDOWING
C           1 HANNING WINDOW (TUNENEWT)
C           2 NO WINDOW (TUNELASR)
C    NRC   = ORDER OF LINEAR COMBINATIONS TO LOOK FOR
C    EPS   = TOLERANCE ON LINEAR COMBINATIONS
C
C  N.B.: ONLY ONE FILE IS ANALIZED AND THE OUTPUT IS PLACED IN
C        THE COMMONS TUNE AND FCOE TOGETHER WITH THE TRACKING DATA.
C  N.B.: NTURN IS AN OUTPUT PARAMETER!
C
C  AUTHOR: R.BARTOLINI 21/08/1996
C  MODIFIED 17/09/1996: ADDED THE IR OPTION
C  MODIFIED 30/05/1998: ADDED FFT
C
C=======================================================================
      implicit none
      integer iana,idam,imeth,ir,iunit,j,k,
     &maxn,mterm,narm,nrc,nt1,nt2,nturn,nturn2,maxturns
      parameter (maxturns=10000)
      double precision duepi,xy(maxturns*4+4),eps
      complex zsing
      parameter(mterm=300)
      parameter(maxn=100000)
      double precision s,sp,x,xp,y,yp
      common/data/x(maxn),y(maxn),xp(maxn),yp(maxn),s(maxn),sp(maxn)
      double precision tsa,txa,tya
      common/tune/txa(mterm),tya(mterm),tsa(mterm)
      double complex zspes,zxpes,zypes
      common/fcoe/zxpes(mterm),zypes(mterm),zspes(mterm)
!$OMP THREADPRIVATE(/data/,/fcoe/,/tune/,j,nturn2,k,zsing,duepi)
      save
      duepi=8*atan(1d+0)
 
C.....READ INPUT FILE
 
C      write(6,*)'****************************'
C      write(6,*)'ANALYZING UNIT',iunit
C      write(6,*)'****************************'
      if(idam.eq.1) then
         do j=1,maxturns
C        do j=1,4000
C        do j=1,maxn
C         read(iunit,*,end=990)x(j),xp(j)
           x(j)=xy(j)
           xp(j)=xy(j)
        enddo
      elseif(idam.eq.2) then
         do j=1,maxturns
C        do j=1,4000
c        do j=1,maxn
C          read(iunit,*,end=990)x(j),xp(j),y(j),yp(j)
           x(j)=xy(j)
           xp(j)=xy(j+maxturns*2)
           y(j)=xy(j+maxturns)
           yp(j)=xy(j+maxturns*3)
        enddo
      elseif(idam.eq.3) then
        do j=1,maxn
C          read(iunit,*,end=990)x(j),xp(j),y(j),yp(j),s(j),sp(j)
        enddo
      endif
c990   nturn=j-1       ! check this it is always larger by 1
      nturn=j-1
C      write(6,*)'NUMBER OF TURNS DETECTED IN THE INPUT',nturn
C      rewind(iunit)
 
C.....CHECK FOR REQUIRED SECTIONING OF THE SIGNAL
      nturn2=nt2-nt1+1
      if(nturn.gt.nturn2) then
        if(idam.eq.1) then
          do k=1,nturn2
            x(k)=x(nt1+k-1)
            xp(k)=xp(nt1+k-1)
          enddo
        elseif(idam.eq.2) then
          do k=1,nturn2
            x(k)=x(nt1+k-1)
            xp(k)=xp(nt1+k-1)
            y(k)=y(nt1+k-1)
            yp(k)=yp(nt1+k-1)
          enddo
        elseif(idam.eq.3) then
          do k=1,nturn2
            x(k)=x(nt1+k-1)
            xp(k)=xp(nt1+k-1)
            y(k)=y(nt1+k-1)
            yp(k)=yp(nt1+k-1)
            s(k)=s(nt1+k-1)
            sp(k)=sp(nt1+k-1)
          enddo
        endif
        nturn=nturn2
C        write(6,*)'REDUCTION OF THE SIGNAL PERFORMED.'
C        write(6,*)'NEW NUMBER OF TURNS ANALYZED =',nturn
C        write(6,*)'INTERVAL',nt1,nt2
      endif
 
C.....FLAG FOR THE ANALYSIS OF THE REAL PART ONLY
 
      if(ir.eq.1) then
        do j=1,nturn
          xp(j)=0d+0
          yp(j)=0d+0
          sp(j)=0d+0
        enddo
C        write(6,*)'WARNING: ONLY THE REAL PART OF THE SIGNAL IS'
C        write(6,*)'         CONSIDERED FOR SPECTRUM CALCULATION'
      endif
 
C.....SPECTRUM CALCULATION
 
      if(idam.eq.1) then
        call spectrum(x,xp,nturn,txa,zxpes,narm,imeth)
C        write(6,*)'TURNS ANALYZED AND TUNE'
C        write(6,*)nturn,-txa(1)
        if(iana.eq.2) then
          call fftr(x,xp,nturn,zsing,imeth)
        endif
      else if(idam.eq.2) then
        call spectrum(x,xp,nturn,txa,zxpes,narm,imeth)
        call spectrum(y,yp,nturn,tya,zypes,narm,imeth)
C        write(6,*)'TURNS ANALYZED AND TUNE'
C        write(6,*)nturn,-txa(1),-tya(1)
        if(iana.eq.2) then
          call fftr(x,xp,nturn,zsing,imeth)
          call fftr(y,yp,nturn,zsing,imeth)
        endif
      elseif(idam.eq.3) then
        call spectrum(x,xp,nturn,txa,zxpes,narm,imeth)
        call spectrum(y,yp,nturn,tya,zypes,narm,imeth)
        call spectrum(s,sp,nturn,tsa,zspes,narm,imeth)
C        write(6,*)'TURNS ANALYZED AND TUNE'
C        write(6,*)nturn,-txa(1),-tya(1),-tsa(1)
        if(iana.eq.2) then
          call fftr(x,xp,nturn,zsing,imeth)
          call fftr(y,yp,nturn,zsing,imeth)
          call fftr(s,sp,nturn,zsing,imeth)
        endif
      endif
 
      end
      subroutine fftr(x,xp,maxn,zsing,meth)
C=======================================================================
C
C SUBROUTINE FFTR
C
C COMPUTES THE FFT.
C X, XP ARE THE COORDINATES OF THE ORBIT AND MAXN IS THE
C LENGTH IF THE ORBIT.
C
C METH SELECTS THE HANNING WINDOW (1) OR NOT (2)
C
C FOURIER COEFFICIENTS ARE GIVEN IN ZSING
C AND THE SPECTRUM IS WRITTEN IN FORT.301
C
C AUTHOR:     R. BARTOLINI
C
C=======================================================================
      implicit none
      integer maxiter,maxn,maxn2,meth,mf,mft,nf,npoint
      double precision amp,duepi,omnf,pha,step,x,xp
      double complex z
      complex zsing
      parameter(maxiter=100000)
      dimension x(maxiter),xp(maxiter),zsing(maxiter)
      dimension z(maxiter)
      save
!$OMP THREADPRIVATE(z,amp,omnf,pha,step,maxn2,mf,mft,nf,npoint,duepi)
C.............................................................
C    ESTIMATION OF TUNE WITH FFT
C.............................................................
      duepi=atan(1d0)*8d0
      mft=int(log(dble(maxn))/log(2d0))
      npoint=2**mft
      maxn2=maxn/2
      step=duepi/maxn
      if (meth.eq.1) then
        do mf=1,maxn
          z(mf)=dcmplx(x(mf),xp(mf))*(1d0+cos(step*(mf-maxn2)))
          zsing(mf)=z(mf)
        enddo
      elseif (meth.eq.2) then
        do mf=1,maxn
          z(mf)=dcmplx(x(mf),xp(mf))
          zsing(mf)=z(mf)
        enddo
      endif
 
      call cfft(zsing,-mft)
 
      do nf=1,npoint
        omnf=(nf-1d0)/npoint
        amp=abs(zsing(nf))
        pha=atan2(imag(zsing(nf)),real(zsing(nf)))/duepi*360
C        write(301,*)omnf,amp/npoint,pha
      enddo
 
C............................................................
      return
C............................................................
      end
      subroutine readres(imin,imax,lr,mr,kr,narm,idam,idamx,iout)
C=======================================================================
C
C  SUBROUTINE READRES
C
C  READ THE HARMONICS FROM THE FILE FORT.300
C  AND WRITE SEPARATELY THE SELECTED ONE
C
C  INPUT PARAMETERS:
C
C  IMIN,IMAX ENUMERATE THE DIFFERENT CASES STORED IN THE FORT.300
C  NARM : NUMERO DI ARMONICHE DA ORDINARE
C  IR   : FLAG FOR COMPLEX OR REAL SIGNALS:
C         1 = INPUT SIGNAL IS REAL.
C
C  N.B. :  THE STRENGTH OF NEXT TO LEADING LINES ARE NORMALIZED AT THE
C          TUNE STRENGTH.
C  MIND :  THE ARRAYS TX,TY,TZ ARE USED IN ORDER NOT TO CHANGE
C          THE ARRAYS TUNEX,TUNEY,TUNEZ.
C
C  IOU  : FORTRAN UNIT OF THE OUTPUT, IF 0 NO OUTPUT
C
C  AUTHOR: R.BARTOLINI 21/08/1996
C  LAST MODIFIED 08/10/1997
C
C=======================================================================
      implicit none
      integer i,idam,idamx,ifoun,imax,imin,iout,istep,j,k,ki,kr,l,lr,m,
     &mr,mterm,n,narm,ni
      double precision ex,fx,pi,px,tx
      parameter(mterm=300)
      dimension tx(mterm)
      double precision dtfx,dtpx,dttx
      common/dt/dtpx(100),dtfx(100),dttx(100)
      character*200 ch,ch1
!$OMP THREADPRIVATE(/dt/,
!$OMP&              i,ifoun,istep,j,k,ki,l,m,n,ni,
!$OMP&              ex,fx,pi,px,tx,ch,ch1)
      save
      pi=4*atan(1d0)
      ch1=' '
 
      if(imin.le.imax) then
        istep=1
      else
        istep=-1
      endif
 
C      rewind(30)
C      do 3 i=imin,imax,istep
CC        write(6,*)'ANALYZING UNIT 300 FOR LINE IDENTIFICATION CASE:',i
C        do ki=1,idam
C 1        read(30,'(A)',end=2) ch
c Find start for each plane
C          if(ch(1:9).eq.' ANALYSIS') then
C            read(30,'()',end=2)
C            read(30,'()',end=2)
C            read(30,'()',end=2)
C            read(30,'()',end=2)
C            ifoun=0
C            do ni=1,narm
C              read(30,'(A)',end=2) ch
C              ch(198:200)=' / '
c Detect end of data and exit to 3 to study next plane
C              if(ch.eq.ch1) goto 3
C              l=-9999
C              m=-9999
C              k=-9999
C              j=-9999
C              read(ch,*)n,tx(n),px,fx,ex,l,m,k,j
C              if(ki.eq.idamx) then
C                if((l.eq.lr.and.m.eq.mr.and.k.eq.kr).and.
C     &            ifoun.eq.0) then
C                  if(iout.ne.0) then
C                    write(iout,101)tx(1),tx(ni),px,fx,ex,l,m,k,j,
C     &              ' CASE:',i
C                  endif
c Fills the vector for eventual DT calculation
C                  dtpx(i)=px
C                  dtfx(i)=fx
C                  dttx(i)=tx(ni)
C                  ifoun=1
C                endif
C              endif
C            enddo
C            if(ki.eq.idamx.and.ifoun.eq.0.and.iout.ne.0) write(iout,*)
C          else
C            goto 1
C          endif
C        enddo
C3    continue
C2    continue
C     rewind(30)
 
C101   format(4(1x,e18.12),1x,e17.12,:,4(1x,i3),a6,1x,i3)
 
      return
      end
      subroutine readdt3(imin,imax,narm,idam,idamx,iout)
C=======================================================================
C
C  SUBROUTINE READDT3
C
C  READ THE HARMONICS FROM THE FILE FORT.300
C  AND CONVERT THEM INTO DRIVING TERMS (SEE LHC Project Report 132)
C
C  INPUT PARAMETERS:
C
C  IMIN,IMAX ENUMERATE THE DIFFERENT CASES STORED IN THE FORT.300
C  NARM : NUMERO DI ARMONICHE DA ORDINARE
C  IR   : FLAG FOR COMPLEX OR REAL SIGNALS:
C         1 = INPUT SIGNAL IS REAL.
C
C
C  IOUT  : FORTRAN UNIT OF THE OUTPUT
C
C  AUTHOR: R.BARTOLINI 21/08/1996
C
C=======================================================================
      implicit none
      integer i,idam,idamx,imax,imin,iout,istep,narm
      double precision a0012,a0030,a1002,a1011,a1020,a1110,a1200,a2001,
     &a2010,a3000,aaux,c0012,c0030,c1002,c1011,c1020,c1110,c1200,c2001,
     &c2010,c3000,f0012,f0030,f1002,f1011,f1020,f1110,f1200,f2001,f2010,
     &f3000,faux,fx,fy,pi,px,py,s0012,s0030,s1002,s1011,s1020,s1110,
     &s1200,s2001,s2010,s3000,tx,ty,y1,y2
      double complex z1011,z2100,zfli,zline
      dimension tx(50),px(50),fx(50),ty(50),py(50),
     &fy(50),aaux(50),faux(50)
      double precision dtfx,dtpx,dttx
      common/dt/dtpx(100),dtfx(100),dttx(100)
!$OMP THREADPRIVATE(/dt/,
!$OMP&              i,istep,
!$OMP&a0012,a0030,a1002,a1011,a1020,a1110,a1200,a2001,
!$OMP&a2010,a3000,aaux,c0012,c0030,c1002,c1011,c1020,c1110,c1200,c2001,
!$OMP&c2010,c3000,f0012,f0030,f1002,f1011,f1020,f1110,f1200,f2001,f2010,
!$OMP&f3000,faux,fx,fy,pi,px,py,s0012,s0030,s1002,s1011,s1020,s1110,
!$OMP&s1200,s2001,s2010,s3000,tx,ty,y1,y2,
!$OMP&z1011,z2100,zfli,zline)
      save

      pi=4*atan(1d0)
      if(imin.le.imax) then
        istep=1
      else
        istep=-1
      endif
 
      write(iout,*)
      write(iout,*)'THIRD ORDER CONJUGATING FUNCTION COEFFICIENTS'
      write(iout,*)
      write(iout,'(2a)')'jklm          Amplitude              Phase',
     &'           Cos             Sin'
      write(iout,*)
 
      do i=imin,imax,istep
C tune lines
        dtpx(i)=0d0
        dtfx(i)=0d0
        dttx(i)=0d0
        call readres(imin,i,1,0,0,narm,idam,1,0)
        tx(i)=dttx(i)
        px(i)=dtpx(i)
        fx(i)=dtfx(i)
        if(idam.ge.2) then
 
c
c attenzione se analizzo il moto verticale ci
c va la y. cosi come e' faccio un overwrite su dtpx,dtfx,sttx
c
C inoltre la linea eccitata in orizzontale
c e' 1-j+k,m-l e non l-m
c idem per il vertical avro' k-j e non j-k
c correggere.
c
 
          dtpx(i)=0d0
          dtfx(i)=0d0
          dttx(i)=0d0
          call readres(imin,i,0,1,0,narm,idam,2,0)
          ty(i)=dttx(i)
          py(i)=dtpx(i)
          fy(i)=dtfx(i)
        endif
 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C One dimensional resonances horizontal C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
 
C -2,0 line and h3000
        dtpx(i)=0d0
        dtfx(i)=0d0
        dttx(i)=0d0
        call readres(imin,i,-2,0,0,narm,idam,1,0)
        a3000=dtpx(i)/2d0/3d0/px(i)**2
        f3000=dtfx(i)-fx(i)+90d0
        c3000=a3000*cos(f3000/180*pi)
        s3000=a3000*sin(f3000/180*pi)
C 2,0 line and h1200
        dtpx(i)=0d0
        dtfx(i)=0d0
        dttx(i)=0d0
        call readres(imin,i,2,0,0,narm,idam,1,0)
        a1200=dtpx(i)/2d0/px(i)**2
        f1200=dtfx(i)-fx(i)+90d0
        c1200=a1200*cos(f1200/180*pi)
        s1200=a1200*sin(f1200/180*pi)
C stores this data for eventual subresonance calculation
        aaux(i)=a1200
        faux(i)=f1200
 
c other lines in case of 4D motion
        if(idam.gt.1) then
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C One dimensional vertical resonances (vertical motion needed) C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C 0,-2 line and h0030
          dtpx(i)=0d0
          dtfx(i)=0d0
          dttx(i)=0d0
          call readres(imin,i,0,-2,0,narm,idam,2,0)
          a0030=dtpx(i)/2d0/3d0/py(i)**2
          f0030=dtfx(i)-fy(i)+90d0
          c0030=a0030*cos(f0030/180*pi)
          s0030=a0030*sin(f0030/180*pi)
C 0,2 line and h0012
          dtpx(i)=0d0
          dtfx(i)=0d0
          dttx(i)=0d0
          call readres(imin,i,0,2,0,narm,idam,2,0)
          a0012=dtpx(i)/2d0/py(i)**2
          f0012=dtfx(i)-fy(i)+90d0
          c0012=a0012*cos(f0012/180*pi)
          s0012=a0012*sin(f0012/180*pi)
 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Normal resonances from horizontal motion C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
 
C 0,-2 line and h1020
          dtpx(i)=0d0
          dtfx(i)=0d0
          dttx(i)=0d0
          call readres(imin,i,0,-2,0,narm,idam,1,0)
          a1020=dtpx(i)/2d0/py(i)**2
          f1020=dtfx(i)-fx(i)+90d0
          c1020=a1020*cos(f1020/180*pi)
          s1020=a1020*sin(f1020/180*pi)
C 0,2 line and h1002
          dtpx(i)=0d0
          dtfx(i)=0d0
          dttx(i)=0d0
          call readres(imin,i,0,2,0,narm,idam,1,0)
          a1002=dtpx(i)/2d0/py(i)**2
          f1002=dtfx(i)-fx(i)+90d0
          c1002=a1002*cos(f1002/180*pi)
          s1002=a1002*sin(f1002/180*pi)
C 0,0 line and h1011 --subresonance: h2100 needed---
          dtpx(i)=0d0
          dtfx(i)=0d0
          dttx(i)=0d0
          call readres(imin,i,0,0,0,narm,idam,1,0)
          zline=dtpx(i)*cdexp(dcmplx(0d0,dtfx(i)/180*pi))
C.......rebulding the contribution to the (0,0) from h2100
          zfli=4*px(i)**2*dcmplx(0d0,-1d0)
          z2100=aaux(i)*cdexp(dcmplx(0d0,-faux(i)/180*pi))*zfli
C.......subtracting the contribution from h2100 to the (0,0)
          z1011=zline-z2100
C.......building h1011 from z1011
          a1011=abs(z1011)/2d0/py(i)**2
          y1=dble(z1011)
          y2=dimag(z1011)
          f1011=atan2(y2,y1)*180/pi-fx(i)+90d0
          c1011=a1011*cos(f1011/180*pi)
          s1011=a1011*sin(f1011/180*pi)
 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Skew resonances from horizontal motion C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C -1,-1 line and h2010
          dtpx(i)=0d0
          dtfx(i)=0d0
          dttx(i)=0d0
          call readres(imin,i,-1,-1,0,narm,idam,1,0)
          a2010=dtpx(i)/2d0/2d0/px(i)/py(i)
          f2010=dtfx(i)-fx(i)+90d0
          c2010=a2010*cos(f2010/180*pi)
          s2010=a2010*sin(f2010/180*pi)
C  1,-1 line and h1110
          dtpx(i)=0d0
          dtfx(i)=0d0
          dttx(i)=0d0
          call readres(imin,i,1,-1,0,narm,idam,1,0)
          a1110=dtpx(i)/2d0/px(i)/py(i)
          f1110=dtfx(i)-fx(i)+90d0
          c1110=a1110*cos(f1110/180*pi)
          s1110=a1110*sin(f1110/180*pi)
C -1,1 line and h2001
          dtpx(i)=0d0
          dtfx(i)=0d0
          dttx(i)=0d0
          call readres(imin,i,-1,1,0,narm,idam,1,0)
          a2001=dtpx(i)/2d0/2d0/px(i)/py(i)
          f2001=dtfx(i)-fx(i)+90d0
          c2001=a2001*cos(f2001/180*pi)
          s2001=a2001*sin(f2001/180*pi)
 
        endif
        write(iout,100)'3000',a3000,f3000,c3000,s3000
        write(iout,200)'2100',a1200,-f1200
        write(iout,100)'1200',a1200,f1200,c1200,s1200
        write(iout,200)'0300',a3000,-f3000
        write(iout,*)
        write(iout,100)'1020',a1020,f1020,c1020,s1020
        write(iout,100)'1011',a1011,f1011,c1011,s1011
        write(iout,100)'1002',a1002,f1002,c1002,s1002
        write(iout,200)'0120',a1002,-f1002
        write(iout,200)'0111',a1011,-f1011
        write(iout,200)'0102',a1020,-f1020
        write(iout,*)
        write(iout,100)'0030',a0030,f0030,c0030,s0030
        write(iout,200)'0021',a0012,-f0012
        write(iout,100)'0012',a0012,f0012,c0012,s0012
        write(iout,200)'0003',a0030,-f0030
        write(iout,*)
        write(iout,100)'2010',a2010,f2010,c2010,s2010
        write(iout,100)'1110',a1110,f1110,c1110,s1110
        write(iout,200)'0210',a2001,-f2001
        write(iout,100)'2001',a2001,f2001,c2001,s2001
        write(iout,200)'1101',a1110,-f1110
        write(iout,200)'0201',a2010,-f2010
 
      enddo
 
 100  format(1x,a4,4(1x,e19.13))
 200  format(1x,a4,2(1x,e19.13))
      return
      end
      subroutine readdt4(imin,imax,narm,idam,idamx,iout)
C=======================================================================
C
C  SUBROUTINE READDT4
C
C  READ THE HARMONICS FROM THE FILE FORT.300
C  AND CONVERT THEM INTO DRIVING TERMS (SEE LHC Project Report 132)
C
C  IT IS ASSUMED THAT THE DATA ARE ALREADY CLEANED FROM THIRD ORDER
C
C  INPUT PARAMETERS:
C
C  IMIN,IMAX ENUMERATE THE DIFFERENT CASES STORED IN THE FORT.300
C  NARM : NUMERO DI ARMONICHE DA ORDINARE
C  IR   : FLAG FOR COMPLEX OR REAL SIGNALS:
C         1 = INPUT SIGNAL IS REAL.
C
C
C  IOUT  : FORTRAN UNIT OF THE OUTPUT
C
C  AUTHOR: R.BARTOLINI 21/08/1996
C
C=======================================================================
      implicit none
      integer i,idam,idamx,imax,imin,iout,istep,narm
      double precision a0013,a0040,a0301,a1003,a1012,a1021,a1030,a1120,
     &a1201,a1210,a1300,a2002,a2011,a2020,a2101,a3001,a3010,a4000,aaux,
     &aauy,c0013,c0040,c1003,c1012,c1021,c1030,c1120,c1201,c1210,c1300,
     &c2002,c2011,c2020,c3001,c3010,c4000,f0013,f0040,f1003,f1012,
     &f1021,f1030,f1120,f1201,f1210,f1300,f2002,f2011,f2020,f3001,f3010,
     &f4000,faux,fauy,fx,fy,pi,px,py,s0013,s0040,s1003,s1012,s1021,
     &s1030,s1120,s1201,s1210,s1300,s2002,s2011,s2020,s3001,s3010,s4000,
     &tx,ty,y1,y2
      double complex z1012,z1021,z2011,z2101,z2110,z3100,zfli,zline
      dimension tx(50),px(50),fx(50),ty(50),py(50),
     &fy(50),aaux(50),faux(50),aauy(50),fauy(50)
      double precision dtfx,dtpx,dttx
      common/dt/dtpx(100),dtfx(100),dttx(100)
!$OMP THREADPRIVATE(/dt/,
!$OMP&              i,istep,
!$OMP&a0013,a0040,a0301,a1003,a1012,a1021,a1030,a1120,
!$OMP&a1201,a1210,a1300,a2002,a2011,a2020,a2101,a3001,a3010,a4000,aaux,
!$OMP&aauy,c0013,c0040,c1003,c1012,c1021,c1030,c1120,c1201,c1210,c1300,
!$OMP&c2002,c2011,c2020,c3001,c3010,c4000,f0013,f0040,f1003,f1012,
!$OMP&f1021,f1030,f1120,f1201,f1210,f1300,f2002,f2011,f2020,f3001,f3010,
!$OMP&f4000,faux,fauy,fx,fy,pi,px,py,s0013,s0040,s1003,s1012,s1021,
!$OMP&s1030,s1120,s1201,s1210,s1300,s2002,s2011,s2020,s3001,s3010,s4000,
!$OMP&tx,ty,y1,y2,
!$OMP&z1012,z1021,z2011,z2101,z2110,z3100,zfli,zline)
      save
 
      pi=4*atan(1d0)
      if(imin.le.imax) then
        istep=1
      else
        istep=-1
      endif
 
      write(iout,*)
      write(iout,*)'FOURTH ORDER CONJUGATING FUNCTION COEFFICIENTS'
      write(iout,*)
      write(iout,'(2a)')'jklm          Amplitude              Phase',
     &'           Cos             Sin'
 
      write(iout,*)
 
      do i=imin,imax,istep
C tune lines
        dtpx(i)=0d0
        dtfx(i)=0d0
        dttx(i)=0d0
        call readres(imin,i,1,0,0,narm,idam,1,0)
        tx(i)=dttx(i)
        px(i)=dtpx(i)
        fx(i)=dtfx(i)
        if(idam.ge.2) then
          dtpx(i)=0d0
          dtfx(i)=0d0
          dttx(i)=0d0
          call readres(imin,i,0,1,0,narm,idam,2,0)
          ty(i)=dttx(i)
          py(i)=dtpx(i)
          fy(i)=dtfx(i)
        endif
 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C One dimensional resonances horizontal C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
 
C -3,0 line and h4000
        dtpx(i)=0d0
        dtfx(i)=0d0
        dttx(i)=0d0
        call readres(imin,i,-3,0,0,narm,idam,1,0)
        a4000=dtpx(i)/2d0/4d0/px(i)**3
        f4000=dtfx(i)-fx(i)+90d0
        c4000=a4000*cos(f4000/180*pi)
        s4000=a4000*sin(f4000/180*pi)
C 3,0 line and h1300
        dtpx(i)=0d0
        dtfx(i)=0d0
        dttx(i)=0d0
        call readres(imin,i,3,0,0,narm,idam,1,0)
        a1300=dtpx(i)/2d0/px(i)**3
        f1300=dtfx(i)-fx(i)+90d0
        c1300=a1300*cos(f1300/180*pi)
        s1300=a1300*sin(f1300/180*pi)
C stores this data for eventual subresonance calculation
        aaux(i)=a1300
        faux(i)=f1300
 
c other lines in case of 4D motion
        if(idam.gt.1) then
 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C One dimensional vertical resonances (vertical motion needed) C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C 0,-3 line and h0040
          dtpx(i)=0d0
          dtfx(i)=0d0
          dttx(i)=0d0
          call readres(imin,i,0,-3,0,narm,idam,2,0)
          a0040=dtpx(i)/2d0/4d0/py(i)**3
          f0040=dtfx(i)-fy(i)+90d0
          c0040=a0040*cos(f0040/180*pi)
          s0040=a0040*sin(f0040/180*pi)
C 0,3 line and h0013
          dtpx(i)=0d0
          dtfx(i)=0d0
          dttx(i)=0d0
          call readres(imin,i,0,3,0,narm,idam,2,0)
          a0013=dtpx(i)/2d0/py(i)**3
          f0013=dtfx(i)-fy(i)+90d0
          c0013=a0013*cos(f0013/180*pi)
          s0013=a0013*sin(f0013/180*pi)
C stores this data for eventual subresonance calculation
          aauy(i)=a0013
          fauy(i)=f0013
 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Normal resonances from horizontal motion C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
 
C -1,-2 line and h2020
          dtpx(i)=0d0
          dtfx(i)=0d0
          dttx(i)=0d0
          call readres(imin,i,-1,-2,0,narm,idam,1,0)
          a2020=dtpx(i)/2d0/2d0/px(i)/py(i)**2
          f2020=dtfx(i)-fx(i)+90d0
          c2020=a2020*cos(f2020/180*pi)
          s2020=a2020*sin(f2020/180*pi)
C -1,2 line and h2002
          dtpx(i)=0d0
          dtfx(i)=0d0
          dttx(i)=0d0
          call readres(imin,i,-1,2,0,narm,idam,1,0)
          a2002=dtpx(i)/2d0/2d0/px(i)/py(i)**2
          f2002=dtfx(i)-fx(i)+90d0
          c2002=a2002*cos(f2002/180*pi)
          s2002=a2002*sin(f2002/180*pi)
C 1,-2 line and h1120
          dtpx(i)=0d0
          dtfx(i)=0d0
          dttx(i)=0d0
          call readres(imin,i,1,-2,0,narm,idam,1,0)
          a1120=dtpx(i)/2d0/px(i)/py(i)**2
          f1120=dtfx(i)-fx(i)+90d0
          c1120=a1120*cos(f1120/180*pi)
          s1120=a1120*sin(f1120/180*pi)
C -1,0 line and h2011 --subresonance: h3100 needed---
          dtpx(i)=0d0
          dtfx(i)=0d0
          dttx(i)=0d0
          call readres(imin,i,-1,0,0,narm,idam,1,0)
          zline=dtpx(i)*cdexp(dcmplx(0d0,dtfx(i)/180*pi))
C.......rebulding the contribution to the (-1,0) from h3100
          zfli=2d0*3d0*px(i)**3*dcmplx(0d0,-1d0)
          z3100=aaux(i)*cdexp(dcmplx(0d0,-faux(i)/180*pi))*zfli
C.......subtracting the contribution from h3100 to the (-1,0)
          z2011=zline-z3100
C.......building h2011 from z2011
          a2011=abs(z2011)/2d0/2d0/px(i)/py(i)**2
          y1=dble(z2011)
          y2=dimag(z2011)
          f2011=atan2(y2,y1)*180/pi-fx(i)+90d0
          c2011=a2011*cos(f2011/180*pi)
          s2011=a2011*sin(f2011/180*pi)
 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Skew resonances from horizontal motion C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
 
C -2,-1 line and h3010
          dtpx(i)=0d0
          dtfx(i)=0d0
          dttx(i)=0d0
          call readres(imin,i,-2,-1,0,narm,idam,1,0)
          a3010=dtpx(i)/2d0/3d0/px(i)**2/py(i)
          f3010=dtfx(i)-fx(i)+90d0
          c3010=a3010*cos(f3010/180*pi)
          s3010=a3010*sin(f3010/180*pi)
C  2,1 line and h1201
          dtpx(i)=0d0
          dtfx(i)=0d0
          dttx(i)=0d0
          call readres(imin,i,2,1,0,narm,idam,1,0)
          a1201=dtpx(i)/2d0/px(i)**2/py(i)
          f1201=dtfx(i)-fx(i)+90d0
          c1201=a1201*cos(f1201/180*pi)
          s1201=a1201*sin(f1201/180*pi)
C  2,-1 line and h1210
          dtpx(i)=0d0
          dtfx(i)=0d0
          dttx(i)=0d0
          call readres(imin,i,2,-1,0,narm,idam,1,0)
          a1210=dtpx(i)/2d0/px(i)**2/py(i)
          f1210=dtfx(i)-fx(i)+90d0
          c1210=a1210*cos(f1210/180*pi)
          s1210=a1210*sin(f1210/180*pi)
C -2,1 line and h3001
          dtpx(i)=0d0
          dtfx(i)=0d0
          dttx(i)=0d0
          call readres(imin,i,-2,1,0,narm,idam,1,0)
          a3001=dtpx(i)/2d0/3d0/px(i)**2/py(i)
          f3001=dtfx(i)-fx(i)+90d0
          c3001=a3001*cos(f3001/180*pi)
          s3001=a3001*sin(f3001/180*pi)
C 0,-3 line and h1030
          dtpx(i)=0d0
          dtfx(i)=0d0
          dttx(i)=0d0
          call readres(imin,i,0,-3,0,narm,idam,1,0)
          a1030=dtpx(i)/2d0/py(i)**3
          f1030=dtfx(i)-fx(i)+90d0
          c1030=a1030*cos(f1030/180*pi)
          s1030=a1030*sin(f1030/180*pi)
C 0,-1 line and h1021 --subresonance h2110 needed--
          dtpx(i)=0d0
          dtfx(i)=0d0
          dttx(i)=0d0
          call readres(imin,i,0,-1,0,narm,idam,1,0)
          zline=dtpx(i)*cdexp(dcmplx(0d0,dtfx(i)/180*pi))
C.......rebulding the contribution to the (0,-1) from h2110
          zfli=2d0*2d0*px(i)**2*py(i)*dcmplx(0d0,-1d0)
          z2110=a1201*cdexp(dcmplx(0d0,-f1201/180*pi))*zfli
C.......subtracting the contribution from h2110 to the (0,-1)
          z1021=zline-z2110
C.......building h1021 from z1021
          a1021=abs(z1021)/2d0/py(i)**3
          y1=dble(z1021)
          y2=dimag(z1021)
          f1021=atan2(y2,y1)*180/pi-fx(i)+90d0
          c1021=a1021*cos(f1021/180*pi)
          s1021=a1021*sin(f1021/180*pi)
C 0,1 line and h1012 --subresonance h2101 needed--
          dtpx(i)=0d0
          dtfx(i)=0d0
          dttx(i)=0d0
          call readres(imin,i,0,1,0,narm,idam,1,0)
          zline=dtpx(i)*cdexp(dcmplx(0d0,dtfx(i)/180*pi))
C.......rebulding the contribution to the (0,1) from h2101
          zfli=2d0*2d0*px(i)**2*py(i)*dcmplx(0d0,-1d0)
          z2101=a1210*cdexp(dcmplx(0d0,-f1210/180*pi))*zfli
C.......subtracting the contribution from h2101 to the (0,1)
          z1012=zline-z2101
C.......building h1012 from z1012
          a1012=abs(z1012)/2d0/py(i)**3
          y1=dble(z1012)
          y2=dimag(z1012)
          f1012=atan2(y2,y1)*180/pi-fx(i)+90d0
          c1012=a1012*cos(f1012/180*pi)
          s1012=a1012*sin(f1012/180*pi)
C 0,3 line and h1003
          dtpx(i)=0d0
          dtfx(i)=0d0
          dttx(i)=0d0
          call readres(imin,i,0,3,0,narm,idam,1,0)
          a1003=dtpx(i)/2d0/py(i)**3
          f1003=dtfx(i)-fx(i)+90d0
          c1003=a1003*cos(f1003/180*pi)
          s1003=a1003*sin(f1003/180*pi)
 
        endif
        write(iout,100)'4000',a4000,f4000,c4000,s4000
        write(iout,200)'3100',a1300,-f1300
        write(iout,100)'1300',a1300,f1300,c1300,s1300
        write(iout,200)'0400',a4000,-f4000
        write(iout,*)
        write(iout,100)'2020',a2020,f2020,c2020,s2020
        write(iout,100)'2011',a2011,f2011,c2011,s2011
        write(iout,100)'2002',a2002,f2002,c2002,s2002
        write(iout,100)'1120',a1120,f1120,c1120,s1120
        write(iout,200)'1102',a1120,-f1120
        write(iout,200)'0220',a2002,-f2002
        write(iout,200)'0211',a2011,-f2011
        write(iout,200)'0202',a2020,-f2020
        write(iout,*)
        write(iout,100)'0040',a0040,f0040,c0040,s0040
        write(iout,200)'0031',a0013,-f0013
        write(iout,100)'0013',a0013,f0013,c0013,s0013
        write(iout,200)'0004',a0040,-f0040
        write(iout,*)
        write(iout,100)'3010',a3010,f3010,c3010,s3010
        write(iout,200)'2110',a1201,-f1201
        write(iout,100)'1210',a1210,f1210,c1210,s1210
        write(iout,200)'0310',a3001,-f3001
        write(iout,100)'3001',a3001,f3001,c3001,s3001
        write(iout,200)'2101',a2101,-f1210
        write(iout,100)'1201',a1201,f1201,c1201,s1201
        write(iout,200)'0301',a0301,-f3010
        write(iout,*)
        write(iout,100)'1030',a1030,f1030,c1030,s1030
        write(iout,100)'1021',a1021,f1021,c1021,s1021
        write(iout,100)'1012',a1012,f1012,c1012,s1012
        write(iout,100)'1003',a1003,f1003,c1003,s1003
        write(iout,200)'0103',a1030,-f1030
        write(iout,200)'0112',a1021,-f1021
        write(iout,200)'0121',a1012,-f1012
        write(iout,200)'0130',a1003,-f1003
      enddo
 
 100  format(1x,a4,4(1x,e19.13))
 200  format(1x,a4,2(1x,e19.13))
      return
      end
      subroutine readsme(imin,imax,narm,idam,idamx,iusm)
C=======================================================================
C
C  READ THE HARMONICS FROM FORT.300
C  GIVES THE SMEAR AS A QUALITY FACTOR FOR THE DIFFERENT CASES
C
C  INPUT PARAMETERS:
C
C  NARM : NUMERO DI ARMONICHE DA ORDINARE
C  IR   : FLAG FOR COMPLEX OR REAL SIGNALS:
C         1 = INPUT SIGNAL IS REAL.
C         IF THE SIGNAL IS REAL IT IS TRASPOSED IN [0,0.5]
C
C  N.B. :  THE STRENGTH OF NEXT TO LEADING LINES ARE NORMALIZED AT THE
C          TUNE STRENGTH.
C  MIND :  THE ARRAYS TX,TY,TZ ARE USED IN ORDER NOT TO CHANGE
C          THE ARRAYS TUNEX,TUNEY,TUNEZ.
C
C  IOU  : FORTRAN UNIT OF THE OUTPUT
C
C  AUTHOR: R.BARTOLINI 21/08/1996
C
C=======================================================================
      implicit none
      integer i,idam,idamx,imax,imin,istep,iusm,j,k,ki,l,m,mterm,n,narm,
     &ni
      double precision ex,fx,pi,px,qsme,tx
      parameter(mterm=300)
      dimension tx(mterm)
      character*200 ch,ch1
!$OMP THREADPRIVATE(i,istep,j,k,ki,l,m,n,ni,
!$OMP&ex,fx,pi,px,qsme,tx,ch,ch1)

      save
 
      pi=4*atan(1d0)
      ch1=' '
 
      if(imin.le.imax) then
        istep=1
      else
        istep=-1
      endif
 
C      rewind(30)
C      do i=imin,imax,istep
CC        write(6,*)'ANALYZING UNIT 300 FOR SMEAR CALCULATION, CASE:',i
C        qsme=0d+0
C        do ki=1,idam
C 1        read(30,'(A)',end=2) ch
Cc Find start for each plane
C          if(ch(1:9).eq.' ANALYSIS') then
C            read(30,'()',end=2)
C           read(30,'()',end=2)
C           read(30,'()',end=2)
C           read(30,'()',end=2)
C           do ni=1,narm
C              read(30,'(A)',end=2) ch
C              if(ch.eq.ch1) goto 3
C              read(ch,100,end=2)n,tx(n),px,fx,ex,l,m,k,j
C              if(ki.eq.idamx) then
C                if(ni.ge.2) then
C                  qsme=qsme+px
C                endif
C              endif
C            enddo
C          else
C            goto 1
C          endif
C        enddo
C 3      continue
C        write(iusm,*)i,qsme
C      enddo
C 2    continue
C      rewind(30)
C 
C100   format(i3,1x,f19.16,2(1x,e18.12),1x,e17.11,:,4(1x,i3))
 
      return
      end
      subroutine readinv(imin,imax,narm,idam,idamx,iinv)
C=======================================================================
C
C  READ THE HARMONICS FROM FORT.300
C  GIVES THE INVARIANTS FOR THE 3 PLANES
C
C  INPUT PARAMETERS:
C
C  IMIN : MIN
C  IMAX : &MAX NUMBER OF CASES
C
C  NARM : NUMERO DI ARMONICHE DA ORDINARE
C
C  IDAM : NUMBER OF PLANES
C
C  IDAMX: DUMMY
C
C  IINV : OUTPUT UNIT
C
C  THEORY: A.BAZZANI, L.BONGINI, G.TURCHETTI, AIP 329, PAGE 120
C  AUTHOR: R.BARTOLINI/F.Schmidt 14/08/1997
C  last change 22/09/1997
C
C=======================================================================
      implicit none
      integer i,idam,idamx,ii,iinv,imax,imin,istep,j,k,ki,kii,kkini,l,
     &lkini,m,mkini,mterm,n,narm,ni
      double precision ex,fx,pi,psinv,px,pxinv,pyinv,tx,zero
      parameter(mterm=300)
      parameter(zero=0d0)
      dimension px(3,mterm)
      dimension l(3,mterm),m(3,mterm),k(3,mterm)
      character*200 ch,ch1
!$OMP THREADPRIVATE(i,ii,istep,j,k,ki,kii,kkini,l,lkini,m,mkini,n,ni,
!$OMP&ex,fx,pi,psinv,px,pxinv,pyinv,tx,ch,ch1)
      save
 
      pi=4*atan(1d0)
 
      ch1=' '
      if(imin.le.imax) then
        istep=1
      else
        istep=-1
      endif
 
      do i=1,3
        do ii=1,mterm
          l(i,ii)=0
          m(i,ii)=0
          k(i,ii)=0
        enddo
      enddo
 
      do i=1,3
        do ii=1,mterm
           px(i,ii)=zero
        enddo
      enddo
 
C      rewind(30)
Cc study imax-imin+1 sets of data
C      do i=imin,imax,istep
CC        write(6,*)'ANALYZING UNIT 300. INVARIANT CALCULATION, CASE:',i
Cc loop over planes
C        do 3 ki=1,idam
C 1        read(30,'(A)',end=2) ch
Cc Find start for each plane
C          if(ch(1:9).eq.' ANALYSIS') then
C            read(30,'()',end=2)
C            read(30,'()',end=2)
C            read(30,'()',end=2)
C            read(30,'()',end=2)
C            do 4 ni=1,narm
C              read(30,'(A)',end=2) ch
Cc Detect end of data and exit to 3 to study next plane
C              if(ch.eq.ch1) goto 3
C              lkini=0
C              mkini=0
C              kkini=0
C              read(ch,100)n,tx,px(ki,ni),fx,ex,lkini,mkini,kkini,j
Cc Reject new line if already present in analysed (ni-1) lines
C              if(ni.gt.1) then
C                do kii=1,ni-1
C                  if(l(ki,kii).eq.lkini.and.m(ki,kii).eq.mkini.and.
C     &                 k(ki,kii).eq.kkini) then
Cc true when line already exists => exit to 4 to look for new line
C                    l(ki,ni)=0
C                    m(ki,ni)=0
C                    k(ki,ni)=0
C                    goto 4
C                  endif
C                enddo
C              endif
Cc set new coefficients
C              l(ki,ni)=lkini
C              m(ki,ni)=mkini
C              k(ki,ni)=kkini
C 4          continue
C          else
C            goto 1
C          endif
C 3      continue
 
CC.......INVARIANTS CALCULATION
C        pxinv=0d0
C        pyinv=0d0
C        psinv=0d0
C        do n=1,narm
C          pxinv=pxinv+l(1,n)*px(1,n)**2+
C     &                l(2,n)*px(2,n)**2+
C     &                l(3,n)*px(3,n)**2
C          pyinv=pyinv+m(1,n)*px(1,n)**2+
C     &                m(2,n)*px(2,n)**2+
C     &                m(3,n)*px(3,n)**2
C          psinv=psinv+k(1,n)*px(1,n)**2+
C     &                k(2,n)*px(2,n)**2+
C     &                k(3,n)*px(3,n)**2
C        enddo
C        if(pxinv.lt.zero) then
CC          write(6,*)'Warning: Horizontal Invariant negative'
C          pxinv=zero
C        else
C          pxinv=sqrt(pxinv)
C        endif
C        if(pyinv.lt.zero) then
CC          write(6,*)'Warning: Vertical Invariant negative'
C          pyinv=zero
C        else
C          pyinv=sqrt(pyinv)
C        endif
C        if(psinv.lt.zero) then
CC          write(6,*)'Warning: Longitudinal Invariant negative'
C          psinv=zero
C        else
C          psinv=sqrt(psinv)
C        endif
C 
C        write(iinv,'(3(1X,E20.12))')pxinv,pyinv,psinv
C 
C      enddo
C      rewind(30)
C 2    continue
 
C100   format(i3,1x,f19.16,2(1x,e18.12),1x,e17.11,:,4(1x,i3))
 
      return
      end
      subroutine readric(nfile,idam,ntwin,iconv)
C=======================================================================
C
C SUBROUTINE READRIC
C
C GIVEN A SIXTRACK BINARY FILE IT CONVERTS IT INTO ONE OR TWO ASCII
C ACCORDING TO NTWIN.
C
C NFILE IS THE UNIT OF THE FORT.NFILE TO BE PROCESSED
C
C THE OUTPUT ARE ALWAYS IN FORT.91 AND FORT.92
C
C N.B.: NTWIN IS AN OUTPUT PARAMETER AND THE NUMBER OF TURN
C       IS AUTOMATICALLY DETERMINED FROM READING THE WHOLE FILE.
C
C=======================================================================
      implicit none
      integer i,ia,icode,iconv,idam,ifipa,ilapa,iq,itopa,its6d,j,jq,
     &nerror,nfile,nfile0,nfile1,ntwin,numl,rdummy
      double precision b,c,c1,clo,clop,d,d1,dizu0,di0x,di0z,dip0x,dip0z,
     &dmmac,dnms,dnumlr,dummy,e,e1,f,f1,g,g1,h,h1,one,p,pi,pieni,p1,qwc,
     &t,ta,tasum,two,txyz,xyzv,zero
      parameter(zero=0d0)
      parameter(pieni=1d-17,one=1d0,two=2d0)
      character*80 sixtit,comment
      character*8 cdate,ctime,prgram
      dimension qwc(3),clo(3),clop(3)
      dimension ta(6,6),t(6,6),txyz(6),xyzv(6)
!$OMP THREADPRIVATE(i,ia,icode,ifipa,ilapa,iq,itopa,its6d,j,jq,
!$OMP&              nerror,nfile0,nfile1,numl,rdummy,
!$OMP&b,c,c1,clo,clop,d,d1,dizu0,di0x,di0z,dip0x,dip0z,
!$OMP&dmmac,dnms,dnumlr,dummy,e,e1,f,f1,g,g1,h,h1,p,pi,p1,qwc,
!$OMP&t,ta,tasum,txyz,xyzv,
!$OMP&sixtit,comment,cdate,ctime,prgram)
      save

      pi=atan(1d0)*4d0
      dummy=zero
 
C---------------------------------------------------------------------
 
      nfile0=91
      nfile1=92
      rewind(nfile)
      rewind(nfile0)
      rewind(nfile1)
 
      read(nfile) sixtit,comment,cdate,ctime,
     &        prgram,ifipa,ilapa,itopa,icode,numl,qwc(1),qwc(2),qwc(3),
     &        clo(1),clop(1),clo(2),clop(2),clo(3),clop(3),
     &        di0x,dip0x,di0z,dip0z,dummy,dummy,
     &        ta(1,1),ta(1,2),ta(1,3),
     &        ta(1,4),ta(1,5),ta(1,6),
     &        ta(2,1),ta(2,2),ta(2,3),
     &        ta(2,4),ta(2,5),ta(2,6),
     &        ta(3,1),ta(3,2),ta(3,3),
     &        ta(3,4),ta(3,5),ta(3,6),
     &        ta(4,1),ta(4,2),ta(4,3),
     &        ta(4,4),ta(4,5),ta(4,6),
     &        ta(5,1),ta(5,2),ta(5,3),
     &        ta(5,4),ta(5,5),ta(5,6),
     &        ta(6,1),ta(6,2),ta(6,3),
     &        ta(6,4),ta(6,5),ta(6,6),
     &        dmmac,dnms,dizu0,dnumlr,
     &        dummy,dummy,dummy,dummy,dummy,dummy,
     &        dummy,dummy,dummy,dummy,dummy,dummy,
     &        dummy,dummy,dummy,dummy,dummy,dummy,
     &        dummy,dummy,dummy,dummy,dummy,dummy,
     &        dummy,dummy,dummy,dummy,dummy,dummy,
     &        dummy,dummy,dummy,dummy,dummy,dummy,
     &        dummy,dummy,dummy,dummy,dummy,dummy,
     &        dummy,dummy,dummy,dummy
 
C--REPAIR PARTIAL OR TOTAL INCOMPLETE TRANSFER MATRIX
 
      if(abs(ta(1,1)).le.pieni.and.abs(ta(2,2)).le.pieni) then
        ta(1,1)=one
        ta(2,2)=one
      endif
      if(abs(ta(3,3)).le.pieni.and.abs(ta(4,4)).le.pieni) then
        ta(3,3)=one
        ta(4,4)=one
      endif
      if(abs(ta(5,5)).le.pieni.and.abs(ta(6,6)).le.pieni) then
        ta(5,5)=one
        ta(6,6)=one
      endif
      its6d=0
 
C--Convert if requested
 
      if(iconv.ne.1) then
        do i=1,6
          do j=1,6
            if(i.ne.j) then
              ta(i,j)=zero
            else
              ta(i,j)=one
            endif
          enddo
        enddo
      endif
 
C--TEST IF TRANSFER MATRIX HAS LONGITUDINAL PART
 
      do 110 i=1,6
        tasum=tasum+abs(ta(i,5))+abs(ta(i,6))
110   continue
      do 120 i=1,4
        tasum=tasum+abs(ta(5,i))+abs(ta(6,i))
120   continue
      tasum=tasum-two
      if(abs(tasum).ge.pieni) its6d=1
 
C--INVERT TRANSFER MATRIX
 
      do 130 i=1,6
        do 140 j=1,6
          t(i,j)=ta(j,i)
140     continue
130   continue
      call dinv(6,t,6,rdummy,nerror)
 
C--DETERMINE ORDER OF FREEDOM OF TRACKED CASE
 
C      IF(ICODE.EQ.1.OR.ICODE.EQ.2.OR.ICODE.EQ.4) IDAM=1
C      IF(ICODE.EQ.3.OR.ICODE.EQ.5.OR.ICODE.EQ.6) IDAM=2
C      IF(ICODE.EQ.7) IDAM=3
 
C--CHECK IF 1 OR 2 SETS OF COORDINATES ARE WRITTEN
 
      ntwin=1
      if(ilapa.ne.ifipa) ntwin=2
C      write(6,*)'ANALYZING THE BINARY SIXTRACK OUTPUT'
      do 150 i=1,numl+1
        if(ntwin.eq.1)
     &    read(nfile,end=999) ia,ifipa,b,c,d,e,f,g,h,p
        if(ntwin.eq.2)
     &    read(nfile,end=999) ia,ifipa,b,c,d,e,f,g,h,p,
     &    ilapa,b,c1,d1,e1,f1,g1,h1,p1
 
C--SUBTRACT CLOSED ORBIT
 
          c=c-clo(1)
          d=d-clop(1)
          e=e-clo(2)
          f=f-clop(2)
          g=g-clo(3)
          h=h-clop(3)
 
C--!! SUBTRACT DISPERSION ONLY IF THE LONGITUDINAL
C--PART OF THE TRANSFER MATRIX HAS NOT BEEN
C--CALCULATED BUT THERE ARE SYNCHROTRON OSCILLATIONS
 
          if(icode.ge.4.and.its6d.eq.0) then
            c=c-di0x*h
            d=d-dip0x*h
            e=e-di0z*h
            f=f-dip0z*h
          endif
          xyzv(1)=c
          xyzv(2)=d
          xyzv(3)=e
          xyzv(4)=f
          xyzv(5)=g
          xyzv(6)=h
 
C--TRANSFER FIRST SET OF DATA
 
          do 160 iq=1,6
            txyz(iq)=zero
            do 170 jq=1,6
              txyz(iq)=txyz(iq)+t(jq,iq)*xyzv(jq)
170         continue
160       continue
          c=txyz(1)
          d=txyz(2)
          e=txyz(3)
          f=txyz(4)
          g=txyz(5)
 
C--!!IMPORTANT THE MOMENTUM HAS A TRIVIAL SCALING!!
 
          h=txyz(6)*1d3
          if(ntwin.eq.2) then
             c1=c1-clo(1)
             d1=d1-clop(1)
             e1=e1-clo(2)
             f1=f1-clop(2)
             g1=g1-clo(3)
             h1=h1-clop(3)
 
C--TREAT SECOND SET OF DATA IF THEY EXIST
 
            if(icode.ge.4.and.its6d.eq.0) then
               c1=c1-di0x*h1
               d1=d1-dip0x*h1
               e1=e1-di0z*h1
               f1=f1-dip0z*h1
            endif
            xyzv(1)=c1
            xyzv(2)=d1
            xyzv(3)=e1
            xyzv(4)=f1
            xyzv(5)=g1
            xyzv(6)=h1
 
C--TRANSFER SECOND SET OF DATA
 
            do 180 iq=1,6
               txyz(iq)=zero
               do 190 jq=1,6
                  txyz(iq)=txyz(iq)+t(jq,iq)*xyzv(jq)
 190           continue
 180        continue
         c1=txyz(1)
         d1=txyz(2)
         e1=txyz(3)
         f1=txyz(4)
         g1=txyz(5)
C--!!IMPORTANT THE MOMENTUM HAS A TRIVIAL SCALING!!
         h1=txyz(6)*1d3
      endif
C--OBEY DEMANDED ORDER OF PHASE SPACE
      if(idam.lt.3) then
        g=zero
        h=zero
        g1=zero
        h1=zero
      endif
      if(idam.eq.1) then
        e=zero
        f=zero
        e1=zero
        f1=zero
      endif
      if(ntwin.eq.1) then
          write(nfile0,'(6(G21.15,1X))') c,d,e,f,g,h
      endif
      if(ntwin.eq.2) then
          write(nfile0,'(6(G21.15,1X))') c,d,e,f,g,h
          write(nfile1,'(6(G21.15,1X))') c1,d1,e1,f1,g1,h1
      endif
 150  continue
 999  continue
      rewind(nfile)
      rewind(nfile0)
      rewind(nfile1)
 
      return
      end
      subroutine writeric(nfile,ntotal,ntwin,numl,iconv)
C=======================================================================
C
C SUBROUTINE WRITERIC
C
C READS THE PROCESSED ASCII FILE FROM THE UNIT 91 AND 92 IF THERE
C IS A TWIN PARTICLE.
C
C AND PLUGS THEM INTO A SIXTRACK BINARY FILE IN UNIT NFILE.
C
C
C=======================================================================
      implicit none
      integer i,i11,i22,ia,icode,iconv,ifipa,ilapa,iq,itopa,its6d,j,jq,
     &ndiff,nfile,ntotal,ntwin,numl,num1
      double precision b,c,c1,clo,clop,d,d1,di0x,di0z,dip0x,dip0z,dizu0,
     &dmmac,dnms,dnumlr,dummy,e,e1,f,f1,g,g1,h,h1,one,p,pi,p1,
     &pieni,qwc,t,ta,tasum,two,txyz,xyzv,zero
      parameter(zero=0d0)
      parameter(pieni=1d-17,one=1d0,two=2d0)
      character*80 sixtit,comment
      character*8 cdate,ctime,prgram
      dimension qwc(3),clo(3),clop(3)
      dimension ta(6,6),t(6,6),txyz(6),xyzv(6)
!$OMP THREADPRIVATE(i,i11,i22,ia,icode,ifipa,ilapa,iq,itopa,its6d,j,jq,
!$OMP&              ndiff,num1,
!$OMP&b,c,c1,clo,clop,d,d1,di0x,di0z,dip0x,dip0z,dizu0,
!$OMP&dmmac,dnms,dnumlr,dummy,e,e1,f,f1,g,g1,h,h1,p,pi,p1,
!$OMP&qwc,t,ta,tasum,txyz,xyzv,sixtit,comment,cdate,ctime,prgram)
      save
 
C.....READS AGAIN THE HEADER AND CHECK FOR THE ACTUAL NUMBER OF TURNS
 
      rewind(nfile)
 
      read(nfile) sixtit,comment,cdate,ctime,
     &        prgram,ifipa,ilapa,itopa,icode,num1,qwc(1),qwc(2),qwc(3),
     &        clo(1),clop(1),clo(2),clop(2),clo(3),clop(3),
     &        di0x,dip0x,di0z,dip0z,dummy,dummy,
     &        ta(1,1),ta(1,2),ta(1,3),
     &        ta(1,4),ta(1,5),ta(1,6),
     &        ta(2,1),ta(2,2),ta(2,3),
     &        ta(2,4),ta(2,5),ta(2,6),
     &        ta(3,1),ta(3,2),ta(3,3),
     &        ta(3,4),ta(3,5),ta(3,6),
     &        ta(4,1),ta(4,2),ta(4,3),
     &        ta(4,4),ta(4,5),ta(4,6),
     &        ta(5,1),ta(5,2),ta(5,3),
     &        ta(5,4),ta(5,5),ta(5,6),
     &        ta(6,1),ta(6,2),ta(6,3),
     &        ta(6,4),ta(6,5),ta(6,6),
     &        dmmac,dnms,dizu0,dnumlr,
     &        dummy,dummy,dummy,dummy,dummy,dummy,
     &        dummy,dummy,dummy,dummy,dummy,dummy,
     &        dummy,dummy,dummy,dummy,dummy,dummy,
     &        dummy,dummy,dummy,dummy,dummy,dummy,
     &        dummy,dummy,dummy,dummy,dummy,dummy,
     &        dummy,dummy,dummy,dummy,dummy,dummy,
     &        dummy,dummy,dummy,dummy,dummy,dummy,
     &        dummy,dummy,dummy,dummy
      read(nfile,end=999) i11,ifipa,b,c,d,e,f,g,h,p
      read(nfile,end=999) i22
      ndiff=i22-i11
      rewind(nfile)
      comment='SUSSIX TREATED DATA'
 
      rewind(91)
      rewind(92)
      rewind(nfile)
C.....WRITES THE SAME HEADER AS THE INPUT FILE
      write(nfile) sixtit,comment,cdate,ctime,
     &        prgram,ifipa,ilapa,itopa,icode,num1,qwc(1),qwc(2),qwc(3),
     &        clo(1),clop(1),clo(2),clop(2),clo(3),clop(3),
     &        di0x,dip0x,di0z,dip0z,dummy,dummy,
     &        ta(1,1),ta(1,2),ta(1,3),
     &        ta(1,4),ta(1,5),ta(1,6),
     &        ta(2,1),ta(2,2),ta(2,3),
     &        ta(2,4),ta(2,5),ta(2,6),
     &        ta(3,1),ta(3,2),ta(3,3),
     &        ta(3,4),ta(3,5),ta(3,6),
     &        ta(4,1),ta(4,2),ta(4,3),
     &        ta(4,4),ta(4,5),ta(4,6),
     &        ta(5,1),ta(5,2),ta(5,3),
     &        ta(5,4),ta(5,5),ta(5,6),
     &        ta(6,1),ta(6,2),ta(6,3),
     &        ta(6,4),ta(6,5),ta(6,6),
     &        dmmac,dnms,dizu0,dnumlr,
     &        dummy,dummy,dummy,dummy,dummy,dummy,
     &        dummy,dummy,dummy,dummy,dummy,dummy,
     &        dummy,dummy,dummy,dummy,dummy,dummy,
     &        dummy,dummy,dummy,dummy,dummy,dummy,
     &        dummy,dummy,dummy,dummy,dummy,dummy,
     &        dummy,dummy,dummy,dummy,dummy,dummy,
     &        dummy,dummy,dummy,dummy,dummy,dummy,
     &        dummy,dummy,dummy,dummy
 
C---------------------------------------------------------------------
 
      pi=atan(1d0)*4d0
 
C.....CONVERTS THE DATA FROM NORMALIZED COORDINATE TO PHISICAL COORDINATES
 
C--REPAIR PARTIAL OR TOTAL INCOMPLETE TRANSFER MATRIX
 
      if(abs(ta(1,1)).le.pieni.and.abs(ta(2,2)).le.pieni) then
        ta(1,1)=one
        ta(2,2)=one
      endif
      if(abs(ta(3,3)).le.pieni.and.abs(ta(4,4)).le.pieni) then
        ta(3,3)=one
        ta(4,4)=one
      endif
      if(abs(ta(5,5)).le.pieni.and.abs(ta(6,6)).le.pieni) then
        ta(5,5)=one
        ta(6,6)=one
      endif
      its6d=0
 
C--Convert if requested
 
      if(iconv.ne.1) then
        do i=1,6
          do j=1,6
            if(i.ne.j) then
              ta(i,j)=zero
            else
              ta(i,j)=one
            endif
          enddo
        enddo
      endif
 
C--TEST IF TRANSFER MATRIX HAS LONGITUDINAL PART
 
      do 110 i=1,6
        tasum=tasum+abs(ta(i,5))+abs(ta(i,6))
110   continue
      do 120 i=1,4
        tasum=tasum+abs(ta(5,i))+abs(ta(6,i))
120   continue
      tasum=tasum-two
      if(abs(tasum).ge.pieni) its6d=1
 
C--TRANSPOSE MATRIX
 
      do 130 i=1,6
        do 140 j=1,6
          t(i,j)=ta(j,i)
140     continue
130   continue
 
C--CHECK IF 1 OR 2 SETS OF COORDINATES ARE WRITTEN
 
      ntwin=1
      if(ilapa.ne.ifipa) ntwin=2
C      write(6,*)'WRITING SUSSIX TREATED DATA INTO THE BINARY FILE'
      do 150 i=1,numl
        if(ntwin.eq.1)
     &    read(91,*,end=999)c,d,e,f,g,h
        if(ntwin.eq.2) then
          read(91,*,end=999) c,d,e,f,g,h
          read(92,*,end=999) c1,d1,e1,f1,g1,h1
        endif
 
C--TRANSFER FIRST SET OF DATA
        xyzv(1)=c
        xyzv(2)=d
        xyzv(3)=e
        xyzv(4)=f
        xyzv(5)=g
C--!!IMPORTANT THE MOMENTUM HAS A TRIVIAL SCALING!!
        xyzv(6)=h/1d3
 
        do 180 iq=1,6
          txyz(iq)=zero
          do 190 jq=1,6
            txyz(iq)=txyz(iq)+t(jq,iq)*xyzv(jq)
190       continue
180     continue
        c=txyz(1)
        d=txyz(2)
        e=txyz(3)
        f=txyz(4)
        g=txyz(5)
        h=txyz(6)
 
C--!! ADD DISPERSION ONLY IF THE LONGITUDINAL
C--PART OF THE TRANSFER MATRIX HAS NOT BEEN
C--CALCULATED BUT THERE ARE SYNCHROTRON OSCILLATIONS
 
        if(icode.ge.4.and.its6d.eq.0) then
          c=c+di0x*h
          d=d+dip0x*h
          e=e+di0z*h
          f=f+dip0z*h
        endif
 
C--ADD CLOSED ORBIT
 
        c=c+clo(1)
        d=d+clop(1)
        e=e+clo(2)
        f=f+clop(2)
        g=g+clo(3)
        h=h+clop(3)
 
        if(ntwin.eq.2) then
C--TRANSFER SECOND SET OF DATA IF THEY EXIST
          xyzv(1)=c1
          xyzv(2)=d1
          xyzv(3)=e1
          xyzv(4)=f1
          xyzv(5)=g1
C--!!IMPORTANT THE MOMENTUM HAS A TRIVIAL SCALING!!
          xyzv(6)=h1/1d3
 
          do 280 iq=1,6
            txyz(iq)=zero
            do 290 jq=1,6
              txyz(iq)=txyz(iq)+t(jq,iq)*xyzv(jq)
290         continue
280       continue
          c1=txyz(1)
          d1=txyz(2)
          e1=txyz(3)
          f1=txyz(4)
          g1=txyz(5)
          h1=txyz(6)
 
C--!! ADD DISPERSION ONLY IF THE LONGITUDINAL
C--PART OF THE TRANSFER MATRIX HAS NOT BEEN
C--CALCULATED BUT THERE ARE SYNCHROTRON OSCILLATIONS
 
          if(icode.ge.4.and.its6d.eq.0) then
            c1=c1+di0x*h1
            d1=d1+dip0x*h1
            e1=e1+di0z*h1
            f1=f1+dip0z*h1
          endif
 
C--ADD CLOSED ORBIT
 
          c1=c1+clo(1)
          d1=d1+clop(1)
          e1=e1+clo(2)
          f1=f1+clop(2)
          g1=g1+clo(3)
          h1=h1+clop(3)
 
        endif
 
C.....WRITES THE DATA IN BINARY FORMAT ON UNIT NFILE
 
        ia=(i-1)*ndiff
        if(ntwin.eq.1) then
          write(nfile) ia,ifipa,b,c,d,e,f,g,h,p
        elseif(ntwin.eq.2) then
          write(nfile) ia,ifipa,b,c,d,e,f,g,h,p,
     &    ilapa,b,c1,d1,e1,f1,g1,h1,p1
        endif
150   continue
999   continue
      rewind(nfile)
      rewind(91)
      rewind(92)
 
      return
      end
      subroutine kerset(ercode,lgfile,limitm,limitr)
      implicit none
      integer i,kounte,l,lgfile,limitm,limitr,log,logf
      parameter(kounte = 27)
      character*6         ercode,   code(kounte)
      logical             mflag,    rflag
      integer             kntm(kounte),       kntr(kounte)
C-----------------------------------------------------------------------
      data      logf      /  0  /
      data      code(1), kntm(1), kntr(1)  / 'C204.1', 255, 255 /
      data      code(2), kntm(2), kntr(2)  / 'C204.2', 255, 255 /
      data      code(3), kntm(3), kntr(3)  / 'C204.3', 255, 255 /
      data      code(4), kntm(4), kntr(4)  / 'C205.1', 255, 255 /
      data      code(5), kntm(5), kntr(5)  / 'C205.2', 255, 255 /
      data      code(6), kntm(6), kntr(6)  / 'C305.1', 255, 255 /
      data      code(7), kntm(7), kntr(7)  / 'C308.1', 255, 255 /
      data      code(8), kntm(8), kntr(8)  / 'C312.1', 255, 255 /
      data      code(9), kntm(9), kntr(9)  / 'C313.1', 255, 255 /
      data      code(10),kntm(10),kntr(10) / 'C336.1', 255, 255 /
      data      code(11),kntm(11),kntr(11) / 'C337.1', 255, 255 /
      data      code(12),kntm(12),kntr(12) / 'C341.1', 255, 255 /
      data      code(13),kntm(13),kntr(13) / 'D103.1', 255, 255 /
      data      code(14),kntm(14),kntr(14) / 'D106.1', 255, 255 /
      data      code(15),kntm(15),kntr(15) / 'D209.1', 255, 255 /
      data      code(16),kntm(16),kntr(16) / 'D509.1', 255, 255 /
      data      code(17),kntm(17),kntr(17) / 'E100.1', 255, 255 /
      data      code(18),kntm(18),kntr(18) / 'E104.1', 255, 255 /
      data      code(19),kntm(19),kntr(19) / 'E105.1', 255, 255 /
      data      code(20),kntm(20),kntr(20) / 'E208.1', 255, 255 /
      data      code(21),kntm(21),kntr(21) / 'E208.2', 255, 255 /
      data      code(22),kntm(22),kntr(22) / 'F010.1', 255,   0 /
      data      code(23),kntm(23),kntr(23) / 'F011.1', 255,   0 /
      data      code(24),kntm(24),kntr(24) / 'F012.1', 255,   0 /
      data      code(25),kntm(25),kntr(25) / 'F406.1', 255,   0 /
      data      code(26),kntm(26),kntr(26) / 'G100.1', 255, 255 /
      data      code(27),kntm(27),kntr(27) / 'G100.2', 255, 255 /
C-----------------------------------------------------------------------
!$OMP THREADPRIVATE(i,l,logf,kntm,kntr,
!$OMP&              code)
      save

      logf  =  lgfile
         l  =  0
      if(ercode .ne. ' ')  then
         do 10  l = 1, 6
            if(ercode(1:l) .eq. ercode)  goto 12
  10        continue
  12     continue
      endif
      do 14     i  =  1, kounte
         if(l .eq. 0)  goto 13
         if(code(i)(1:l) .ne. ercode(1:l))  goto 14
  13     if(limitm.ge.0) kntm(i)  =  limitm
         if(limitr.ge.0) kntr(i)  =  limitr
  14     continue
      return
      entry kermtr(ercode,log,mflag,rflag)
      log  =  logf
      do 20     i  =  1, kounte
         if(ercode .eq. code(i))  goto 21
  20     continue
      write(*,1000)  ercode
      call abend
      return
  21  rflag  =  kntr(i) .ge. 1
      if(rflag  .and.  (kntr(i) .lt. 255))  kntr(i)  =  kntr(i) - 1
      mflag  =  kntm(i) .ge. 1
      if(mflag  .and.  (kntm(i) .lt. 255))  kntm(i)  =  kntm(i) - 1
      if(.not. rflag)  then
         if(logf .lt. 1)  then
            write(*,1001)  code(i)
         else
            write(logf,1001)  code(i)
         endif
      endif
      if(mflag .and. rflag)  then
         if(logf .lt. 1)  then
            write(*,1002)  code(i)
         else
            write(logf,1002)  code(i)
         endif
      endif
      return
1000  format(' KERNLIB LIBRARY ERROR. ' /
     &       ' ERROR CODE ',a6,' NOT RECOGNIZED BY KERMTR',
     &       ' ERROR MONITOR. RUN ABORTED.')
1001  format(/' ***** RUN TERMINATED BY CERN LIBRARY ERROR ',
     &       'CONDITION ',a6)
1002  format(/' ***** CERN LIBRARY ERROR CONDITION ',a6)
      end
      subroutine rinv(n,a,idim,ir,ifail)
C-----------------------------------------------------------------------
C
C     ******************************************************************
C
C     REPLACES A BY ITS INVERSE.
C
C     (PARAMETERS AS FOR REQINV.)
C
C     CALLS ... RFACT, RFINV, F010PR, ABEND.
C
C     ******************************************************************
C-----------------------------------------------------------------------
      implicit none
      integer idim,ifail,ir,jfail,k,kprnt,n
      real t1,t2,t3,a,det,temp,s,c11,c12,c13,c21,c22,c23,c31,c32,c33
      character*6 name
      dimension ir(n),a(idim,n)
      data name/'RINV'/,kprnt/0/
!$OMP THREADPRIVATE(jfail,k,kprnt,
!$OMP& t1,t2,t3,det,temp,s,c11,c12,c13,c21,c22,c23,c31,c32,c33)
      save

C-----------------------------------------------------------------------
C
C  TEST FOR PARAMETER ERRORS.
C
      if((n.lt.1).or.(n.gt.idim)) goto 7
C
C  TEST FOR N.LE.3.
C
      if(n.gt.3) goto 6
      ifail=0
      if(n.lt.3) goto 4
C
C  N=3 CASE.
C
C     COMPUTE COFACTORS.
      c11=a(2,2)*a(3,3)-a(2,3)*a(3,2)
      c12=a(2,3)*a(3,1)-a(2,1)*a(3,3)
      c13=a(2,1)*a(3,2)-a(2,2)*a(3,1)
      c21=a(3,2)*a(1,3)-a(3,3)*a(1,2)
      c22=a(3,3)*a(1,1)-a(3,1)*a(1,3)
      c23=a(3,1)*a(1,2)-a(3,2)*a(1,1)
      c31=a(1,2)*a(2,3)-a(1,3)*a(2,2)
      c32=a(1,3)*a(2,1)-a(1,1)*a(2,3)
      c33=a(1,1)*a(2,2)-a(1,2)*a(2,1)
      t1=abs(a(1,1))
      t2=abs(a(2,1))
      t3=abs(a(3,1))
C
C     (SET TEMP=PIVOT AND DET=PIVOT*DET.)
      if(t1.ge.t2) goto 1
         if(t3.ge.t2) goto 2
C        (PIVOT IS A21)
            temp=a(2,1)
            det=c13*c32-c12*c33
            goto 3
    1 if(t3.ge.t1) goto 2
C     (PIVOT IS A11)
         temp=a(1,1)
         det=c22*c33-c23*c32
         goto 3
C     (PIVOT IS A31)
    2    temp=a(3,1)
         det=c23*c12-c22*c13
C
C     SET ELEMENTS OF INVERSE IN A.
    3 if(det.eq.0.) goto 8
      s=temp/det
      a(1,1)=s*c11
      a(1,2)=s*c21
      a(1,3)=s*c31
      a(2,1)=s*c12
      a(2,2)=s*c22
      a(2,3)=s*c32
      a(3,1)=s*c13
      a(3,2)=s*c23
      a(3,3)=s*c33
      return
C
    4 if(n.lt.2) goto 5
C
C  N=2 CASE BY CRAMERS RULE.
C
      det=a(1,1)*a(2,2)-a(1,2)*a(2,1)
      if(det.eq.0.) goto 8
      s=1./det
      c11   =s*a(2,2)
      a(1,2)=-s*a(1,2)
      a(2,1)=-s*a(2,1)
      a(2,2)=s*a(1,1)
      a(1,1)=c11
      return
C
C  N=1 CASE.
C
    5 if(a(1,1).eq.0.) goto 8
      a(1,1)=1./a(1,1)
      return
C
C  N.GT.3 CASES.  FACTORIZE MATRIX AND INVERT.
C
    6 call rfact(n,a,idim,ir,ifail,det,jfail)
      if(ifail.ne.0) return
      call rfinv(n,a,idim,ir)
      return
C
C  ERROR EXITS.
C
    7 ifail=+1
      call f010pr(name,n,idim,k,kprnt)
      return
C
    8 ifail=-1
      return
C
      end
      subroutine dinv(n,a,idim,ir,ifail)
C-----------------------------------------------------------------------
C
C     ******************************************************************
C
C     REPLACES A BY ITS INVERSE.
C
C     (PARAMETERS AS FOR DEQINV.)
C
C     CALLS ... DFACT, DFINV, F010PR, ABEND.
C
C     ******************************************************************
C-----------------------------------------------------------------------
      implicit none
      integer idim,ifail,jfail,k,kprnt,n
      integer ir
      real t1,t2,t3
      double precision a,det,temp,s,c11,c12,c13,c21,c22,c23,c31,c32,c33
      character*6 name
      dimension ir(n),a(idim,n)
      data name/'DINV'/,kprnt/0/
!$OMP THREADPRIVATE(jfail,k,kprnt,
!$OMP& det,temp,s,c11,c12,c13,c21,c22,c23,c31,c32,c33,t1,t2,t3,name)
      save

C-----------------------------------------------------------------------
C
C  TEST FOR PARAMETER ERRORS.
C
      if((n.lt.1).or.(n.gt.idim)) goto 7
C
C  TEST FOR N.LE.3.
C
      if(n.gt.3) goto 6
      ifail=0
      if(n.lt.3) goto 4
C
C  N=3 CASE.
C
C     COMPUTE COFACTORS.
      c11=a(2,2)*a(3,3)-a(2,3)*a(3,2)
      c12=a(2,3)*a(3,1)-a(2,1)*a(3,3)
      c13=a(2,1)*a(3,2)-a(2,2)*a(3,1)
      c21=a(3,2)*a(1,3)-a(3,3)*a(1,2)
      c22=a(3,3)*a(1,1)-a(3,1)*a(1,3)
      c23=a(3,1)*a(1,2)-a(3,2)*a(1,1)
      c31=a(1,2)*a(2,3)-a(1,3)*a(2,2)
      c32=a(1,3)*a(2,1)-a(1,1)*a(2,3)
      c33=a(1,1)*a(2,2)-a(1,2)*a(2,1)
      t1=abs(sngl(a(1,1)))
      t2=abs(sngl(a(2,1)))
      t3=abs(sngl(a(3,1)))
C
C     (SET TEMP=PIVOT AND DET=PIVOT*DET.)
      if(t1.ge.t2) goto 1
         if(t3.ge.t2) goto 2
C        (PIVOT IS A21)
            temp=a(2,1)
            det=c13*c32-c12*c33
            goto 3
    1 if(t3.ge.t1) goto 2
C     (PIVOT IS A11)
         temp=a(1,1)
         det=c22*c33-c23*c32
         goto 3
C     (PIVOT IS A31)
    2    temp=a(3,1)
         det=c23*c12-c22*c13
C
C     SET ELEMENTS OF INVERSE IN A.
    3 if(det.eq.0d0) goto 8
      s=temp/det
      a(1,1)=s*c11
      a(1,2)=s*c21
      a(1,3)=s*c31
      a(2,1)=s*c12
      a(2,2)=s*c22
      a(2,3)=s*c32
      a(3,1)=s*c13
      a(3,2)=s*c23
      a(3,3)=s*c33
      return
C
    4 if(n.lt.2) goto 5
C
C  N=2 CASE BY CRAMERS RULE.
C
      det=a(1,1)*a(2,2)-a(1,2)*a(2,1)
      if(det.eq.0d0) goto 8
      s=1d0/det
      c11   =s*a(2,2)
      a(1,2)=-s*a(1,2)
      a(2,1)=-s*a(2,1)
      a(2,2)=s*a(1,1)
      a(1,1)=c11
      return
C
C  N=1 CASE.
C
    5 if(a(1,1).eq.0d0) goto 8
      a(1,1)=1d0/a(1,1)
      return
C
C  N.GT.3 CASES.  FACTORIZE MATRIX AND INVERT.
C
    6 call dfact(n,a,idim,ir,ifail,det,jfail)
      if(ifail.ne.0) return
      call dfinv(n,a,idim,ir)
      return
C
C  ERROR EXITS.
C
    7 ifail=+1
      call f010pr(name,n,idim,k,kprnt)
      return
C
    8 ifail=-1
      return
C
      end
      subroutine f010pr(name,n,idim,k,kprnt)
C     ******************************************************************
C
C     PRINT ROUTINE FOR PARAMETER ERRORS IN MATRIX SUBROUTINES $EQINV,
C     $EQN, $INV (WHERE $ IS A LETTER SPECIFYING THE ARITHMETIC TYPE).
C
C     NAME         (CHARACTER*6) NAME OF THE CALLING ROUTINE.
C
C     N,IDIM,K     PARAMETERS OF THE CALLING ROUTINE (WITH K=0 IF K IS
C                  NOT TO BE PRINTED).
C
C     KPRNT        PRINT FLAG FOR K (K IS NOT PRINTED IF KPRNT=0).
C
C     ******************************************************************
      implicit none
      integer idim,k,kprnt,lgfile,n
      character*6 name
      logical mflag,rflag
!$OMP THREADPRIVATE(lgfile,mflag,rflag)
      save

C-----------------------------------------------------------------------
      call kermtr('F010.1',lgfile,mflag,rflag)
      if(mflag) then
         if(lgfile.eq.0)  then
            if(kprnt.eq.0) write(*,2000) name,n,idim
            if(kprnt.ne.0) write(*,2001) name,n,idim,k
         else
            if(kprnt.eq.0) write(lgfile,2000) name,n,idim
            if(kprnt.ne.0) write(lgfile,2001) name,n,idim,k
         endif
      endif
      if(.not. rflag) call abend
      return
C
 2000 format( 7x, 11hsubroutine , a6, 14h ... parameter,
     &        29h error (n.lt.1 or n.gt.idim).,
     &        6x, 3hn =, i4, 6x, 6hidim =, i4, 1h. )
 2001 format( 7x, 11hsubroutine , a6, 14h ... parameter,
     &        39h error (n.lt.1 or n.gt.idim or k.lt.1).,
     &        6x, 3hn =, i4, 6x, 6hidim =, i4, 6x, 3hk =, i4, 1h. )
      end
      subroutine rfact(n,a,idim,ir,ifail,det,jfail)
      implicit none
      integer i,idim,ifail,imposs,ipairf,ir,j,jfail,jm1,jover,jp1,
     &jrange,junder,k,l,n,normal,nxch
      real a,det,g1,g2,one,p,pivotf,q,sizef,t,tf,x,y,zero
      double precision s11,s12,dotf
      character*6 hname
      dimension ir(*),a(idim,*)
      data      g1, g2              /  1.e-37,  1.e37  /
      data      hname               /  ' RFACT'  /
      data      zero, one           /  0., 1.  /
      data      normal, imposs      /  0, -1  /
      data      jrange, jover, junder  /  0, +1, -1  /
!$OMP THREADPRIVATE(i,imposs,j,jm1,jover,jp1,
!$OMP&              jrange,junder,k,l,normal,nxch,
!$OMP& g1,g2,one,p,q,t,tf,x,y,zero,s11,s12,hname)
      save

      dotf(x,y,s11)  =  dble(x)*dble(y) + s11
      ipairf(j,k)  =  j*2**12 + k
      pivotf(x)    =  abs(x)
      sizef(x)     =  abs(x)
C-----------------------------------------------------------------------
      if(idim .ge. n  .and.  n .gt. 0)  goto 110
         call tmprnt(hname,n,idim,0)
         return
 110  ifail  =  normal
      jfail  =  jrange
      nxch   =  0
      det    =  one
      do 144    j  =  1, n
 120     k  =  j
         p  =  pivotf(a(j,j))
         if(j .eq. n)  goto 122
         jp1  =  j+1
         do 121    i  =  jp1, n
            q  =  pivotf(a(i,j))
            if(q .le. p)  goto 121
               k  =  i
               p  =  q
 121        continue
         if(k .ne. j)  goto 123
 122     if(p .gt. 0.)  goto 130
            det    =  zero
            ifail  =  imposs
            jfail  =  jrange
            return
 123     do 124    l  =  1, n
            tf      =  a(j,l)
            a(j,l)  =  a(k,l)
            a(k,l)  =  tf
 124        continue
         nxch      =  nxch + 1
         ir(nxch)  =  ipairf(j,k)
 130     det     =  det * a(j,j)
         a(j,j)  =  one / a(j,j)
         t  =  sizef(det)
         if(t .lt. g1)  then
            det    =  zero
            if(jfail .eq. jrange)  jfail  =  junder
         elseif(t .gt. g2)  then
            det    =  one
            if(jfail .eq. jrange)  jfail  =  jover
         endif
         if(j .eq. n)  goto 144
         jm1  =  j-1
         jp1  =  j+1
         do 143   k  =  jp1, n
            s11  =  -a(j,k)
            s12  =  -a(k,j+1)
            if(j .eq. 1)  goto 142
            do 141  i  =  1, jm1
               s11  =  dotf(a(i,k),a(j,i),s11)
               s12  =  dotf(a(i,j+1),a(k,i),s12)
 141           continue
 142        a(j,k)    =  -s11 * a(j,j)
            a(k,j+1)  =  -dotf(a(j,j+1),a(k,j),s12)
 143        continue
 144     continue
 150  if(mod(nxch,2) .ne. 0)  det  =  -det
      if(jfail .ne. jrange)   det  =  zero
      ir(n)  =  nxch
      return
      end
      subroutine dfact(n,a,idim,ir,ifail,det,jfail)
      implicit none
      integer i,idim,ifail,imposs,ipairf,ir,j,jfail,jm1,jover,jp1,
     &jrange,junder,k,l,n,normal,nxch
      real g1,g2,p,pivotf,q,sizef,t
      double precision a,det,dotf,zero,one,s11,s12,x,y,tf
      character*6         hname
      dimension ir(*),a(idim,*)
      data      g1, g2              /  1.e-37,  1.e37  /
      data      hname               /  ' DFACT'  /
      data      zero, one           /  0.d0, 1.d0  /
      data      normal, imposs      /  0, -1  /
      data      jrange, jover, junder  /  0, +1, -1  /
!$OMP THREADPRIVATE(i,imposs,j,jm1,jover,jp1,
!$OMP&              jrange,junder,k,l,normal,nxch,
!$OMP& g1,g2,p,q,t,
!$OMP& zero,one,s11,s12,x,y,tf,hname)
      save

C-----------------------------------------------------------------------
      ipairf(j,k)  =  j*2**12 + k
      pivotf(x)    =  abs(sngl(x))
      sizef(x)     =  abs(sngl(x))
      dotf(x,y,s11)  =  x * y + s11
      if(idim .ge. n  .and.  n .gt. 0)  goto 110
      call tmprnt(hname,n,idim,0)
      return
 110  ifail  =  normal
      jfail  =  jrange
      nxch   =  0
      det    =  one
      do 144    j  =  1, n
 120     k  =  j
         p  =  pivotf(a(j,j))
         if(j .eq. n)  goto 122
         jp1  =  j+1
         do 121    i  =  jp1, n
            q  =  pivotf(a(i,j))
            if(q .le. p)  goto 121
               k  =  i
               p  =  q
 121        continue
         if(k .ne. j)  goto 123
 122     if(p .gt. 0.)  goto 130
            det    =  zero
            ifail  =  imposs
            jfail  =  jrange
            return
 123     do 124    l  =  1, n
            tf      =  a(j,l)
            a(j,l)  =  a(k,l)
            a(k,l)  =  tf
 124        continue
         nxch      =  nxch + 1
         ir(nxch)  =  ipairf(j,k)
 130     det     =  det * a(j,j)
         a(j,j)  =  one / a(j,j)
         t  =  sizef(det)
         if(t .lt. g1)  then
            det    =  zero
            if(jfail .eq. jrange)  jfail  =  junder
         elseif(t .gt. g2)  then
            det    =  one
            if(jfail .eq. jrange)  jfail  =  jover
         endif
         if(j .eq. n)  goto 144
         jm1  =  j-1
         jp1  =  j+1
         do 143   k  =  jp1, n
            s11  =  -a(j,k)
            s12  =  -a(k,j+1)
            if(j .eq. 1)  goto 142
            do 141  i  =  1, jm1
               s11  =  dotf(a(i,k),a(j,i),s11)
               s12  =  dotf(a(i,j+1),a(k,i),s12)
 141           continue
 142        a(j,k)    =  -s11 * a(j,j)
            a(k,j+1)  =  -dotf(a(j,j+1),a(k,j),s12)
 143        continue
 144     continue
 150  if(mod(nxch,2) .ne. 0)  det  =  -det
      if(jfail .ne. jrange)   det  =  zero
      ir(n)  =  nxch
      return
      end
      subroutine rfeqn(n,a,idim,ir,k,b)
      implicit none
      integer i,idim,ij,im1,ir,j,k,l,m,n,nm1,nmi,nmjp1,nxch
      real a,b,te,x,y
      double precision dotf,s21,s22
      character*6 hname
      dimension ir(*),a(idim,*),b(idim,*)
      data      hname               /  ' RFEQN'  /
!$OMP THREADPRIVATE(i,ij,im1,j,l,m,nm1,nmi,nmjp1,nxch,
!$OMP&              te,x,y,s21,s22,hname)
      save

C-----------------------------------------------------------------------
      dotf(x,y,s21)  =  dble(x)*dble(y) + s21
      if(idim .ge. n  .and.  n .gt. 0  .and.  k .gt. 0)  goto 210
      call tmprnt(hname,n,idim,k)
      return
 210  nxch  =  ir(n)
      if(nxch .eq. 0)  goto 220
      do 212    m  =  1, nxch
         ij  =  ir(m)
         i   =  ij / 4096
         j   =  mod(ij,4096)
         do 211   l  =  1, k
            te      =  b(i,l)
            b(i,l)  =  b(j,l)
            b(j,l)  =  te
 211        continue
 212     continue
 220  do 221    l  =  1, k
         b(1,l)  =  a(1,1)*b(1,l)
 221     continue
      if(n .eq. 1)  goto 299
      do 243    l  =  1, k
         do 232   i  =  2, n
            im1  =  i-1
            s21  =  - b(i,l)
            do 231   j  =  1, im1
               s21  =  dotf(a(i,j),b(j,l),s21)
 231           continue
            b(i,l)  =  - a(i,i)*s21
 232        continue
         nm1  =  n-1
         do 242   i  =  1, nm1
            nmi  =  n-i
            s22  =  - b(nmi,l)
            do 241   j  =  1, i
               nmjp1  =  n - j+1
               s22    =  dotf(a(nmi,nmjp1),b(nmjp1,l),s22)
 241           continue
            b(nmi,l)  =  - s22
 242        continue
 243     continue
 299  continue
      return
      end
      subroutine dfeqn(n,a,idim,ir,k,b)
      implicit none
      integer i,idim,ij,im1,ir,j,k,l,m,n,nm1,nmi,nmjp1,nxch
      double precision a,b,x,y,te
      double precision dotf,s21,s22
      character*6 hname
      dimension ir(*),a(idim,*),b(idim,*)
      data      hname               /  ' DFEQN'  /
!$OMP THREADPRIVATE(i,ij,im1,j,l,m,nm1,nmi,nmjp1,nxch,
!$OMP& x,y,te,s21,s22,hname)
      save

C-----------------------------------------------------------------------
      dotf(x,y,s21)  =  x*y + s21
      if(idim .ge. n  .and.  n .gt. 0  .and.  k .gt. 0)  goto 210
      call tmprnt(hname,n,idim,k)
      return
 210  nxch  =  ir(n)
      if(nxch .eq. 0)  goto 220
      do 212    m  =  1, nxch
         ij  =  ir(m)
         i   =  ij / 4096
         j   =  mod(ij,4096)
         do 211   l  =  1, k
            te      =  b(i,l)
            b(i,l)  =  b(j,l)
            b(j,l)  =  te
 211        continue
 212     continue
 220  do 221    l  =  1, k
         b(1,l)  =  a(1,1)*b(1,l)
 221     continue
      if(n .eq. 1)  goto 299
      do 243    l  =  1, k
         do 232   i  =  2, n
            im1  =  i-1
            s21  =  - b(i,l)
            do 231   j  =  1, im1
               s21  =  dotf(a(i,j),b(j,l),s21)
 231           continue
            b(i,l)  =  - a(i,i)*s21
 232        continue
         nm1  =  n-1
         do 242   i  =  1, nm1
            nmi  =  n-i
            s22  =  - b(nmi,l)
            do 241   j  =  1, i
               nmjp1  =  n - j+1
               s22    =  dotf(a(nmi,nmjp1),b(nmjp1,l),s22)
 241           continue
            b(nmi,l)  =  - s22
 242        continue
 243     continue
 299  continue
      return
      end
      subroutine rfinv(n,a,idim,ir)
      implicit none
      integer i,idim,ij,im2,ir,j,k,m,n,nm1,nmi,nxch
      real a,ti,x,y
      double precision dotf,s31,s32,s33,s34,zero
      character*6 hname
      dimension ir(*),a(idim,*)
      data      zero      /  0.d0  /
      data      hname               /  ' RFINV'  /
!$OMP THREADPRIVATE(i,ij,im2,j,k,m,nm1,nmi,nxch,
!$OMP& ti,x,y,s31,s32,s33,s34,zero,hname)
      save

C-----------------------------------------------------------------------
      dotf(x,y,s31)  =  dble(x)*dble(y) + s31
      if(idim .ge. n  .and.  n .gt. 0)  goto 310
         call tmprnt(hname,n,idim,0)
         return
 310  if(n .eq. 1)  return
      a(2,1)  =  -a(2,2) * dotf(a(1,1),a(2,1),zero)
      a(1,2)  =  -a(1,2)
      if(n .eq. 2)  goto 330
      do 314    i  =  3, n
         im2  =  i-2
         do 312 j  =  1, im2
            s31  =  zero
            s32  =  a(j,i)
            do 311  k  =  j, im2
               s31  =  dotf(a(k,j),a(i,k),s31)
               s32  =  dotf(a(j,k+1),a(k+1,i),s32)
 311           continue
            a(i,j)  =  -a(i,i) * dotf(a(i-1,j),a(i,i-1),s31)
            a(j,i)  =  -s32
 312        continue
         a(i,i-1)  =  -a(i,i) * dotf(a(i-1,i-1),a(i,i-1),zero)
         a(i-1,i)  =  -a(i-1,i)
 314     continue
 330  nm1  =  n-1
      do 335   i  =  1, nm1
         nmi  =  n-i
         do 332   j  =  1, i
            s33  =  a(i,j)
            do 331   k  =  1, nmi
               s33  =  dotf(a(i+k,j),a(i,i+k),s33)
 331           continue
            a(i,j)  =  s33
 332        continue
         do 334   j  =  1, nmi
            s34  =  zero
            do 333   k  =  j, nmi
               s34  =  dotf(a(i+k,i+j),a(i,i+k),s34)
 333           continue
            a(i,i+j)  =  s34
 334        continue
 335     continue
      nxch  =  ir(n)
      if(nxch .eq. 0)  return
        do 342 m  =  1, nxch
         k   =  nxch - m+1
         ij  =  ir(k)
         i   =  ij / 4096
         j   =  mod(ij,4096)
         do 341  k  =  1, n
            ti      =  a(k,i)
            a(k,i)  =  a(k,j)
            a(k,j)  =  ti
 341        continue
 342     continue
      return
      end
      subroutine dfinv(n,a,idim,ir)
      implicit none
      integer i,idim,ij,im2,ir,j,k,m,n,nm1,nmi,nxch
      double precision a,dotf,s31,s32,s33,s34,ti,x,y,zero
      character*6 hname
      dimension ir(*),a(idim,*)
      data      hname               /  ' DFINV'  /
      data      zero      /  0.d0  /
!$OMP THREADPRIVATE(i,ij,im2,j,k,m,nm1,nmi,nxch,
!$OMP& s31,s32,s33,s34,ti,x,y,zero,hname)
      save

C-----------------------------------------------------------------------
      dotf(x,y,s31)  =  x*y + s31
      if(idim .ge. n  .and.  n .gt. 0)  goto 310
         call tmprnt(hname,n,idim,0)
         return
 310  if(n .eq. 1)  return
      a(2,1)  =  -a(2,2) * dotf(a(1,1),a(2,1),zero)
      a(1,2)  =  -a(1,2)
      if(n .eq. 2)  goto 330
      do 314    i  =  3, n
         im2  =  i-2
         do 312 j  =  1, im2
            s31  =  zero
            s32  =  a(j,i)
            do 311  k  =  j, im2
               s31  =  dotf(a(k,j),a(i,k),s31)
               s32  =  dotf(a(j,k+1),a(k+1,i),s32)
 311           continue
            a(i,j)  =  -a(i,i) * dotf(a(i-1,j),a(i,i-1),s31)
            a(j,i)  =  -s32
 312        continue
         a(i,i-1)  =  -a(i,i) * dotf(a(i-1,i-1),a(i,i-1),zero)
         a(i-1,i)  =  -a(i-1,i)
 314     continue
 330  nm1  =  n-1
      do 335   i  =  1, nm1
         nmi  =  n-i
         do 332   j  =  1, i
            s33  =  a(i,j)
            do 331   k  =  1, nmi
               s33  =  dotf(a(i+k,j),a(i,i+k),s33)
 331           continue
            a(i,j)  =  s33
 332        continue
         do 334   j  =  1, nmi
            s34  =  zero
            do 333   k  =  j, nmi
               s34  =  dotf(a(i+k,i+j),a(i,i+k),s34)
 333           continue
            a(i,i+j)  =  s34
 334        continue
 335     continue
      nxch  =  ir(n)
      if(nxch .eq. 0)  return
        do 342 m  =  1, nxch
         k   =  nxch - m+1
         ij  =  ir(k)
         i   =  ij / 4096
         j   =  mod(ij,4096)
         do 341  k  =  1, n
            ti      =  a(k,i)
            a(k,i)  =  a(k,j)
            a(k,j)  =  ti
 341        continue
 342     continue
      return
      end
      subroutine tmprnt(name,n,idim,k)
      implicit none
      integer idim,ifmt,k,lgfile,n
      character*6 name
      logical mflag,rflag
!$OMP THREADPRIVATE(ifmt,lgfile,mflag,rflag)
      save

C-----------------------------------------------------------------------
      if(name(2:2) .eq. 'S') then
         call kermtr('F012.1',lgfile,mflag,rflag)
      else
         call kermtr('F011.1',lgfile,mflag,rflag)
      endif
      if(name(3:6) .eq. 'FEQN') assign 1002 to ifmt
      if(name(3:6) .ne. 'FEQN') assign 1001 to ifmt
      if(mflag) then
         if(lgfile .eq. 0) then
            if(name(3:6) .eq. 'FEQN') then
               write(*,ifmt) name, n, idim, k
            else
               write(*,ifmt) name, n, idim
            endif
         else
            if(name(3:6) .eq. 'FEQN') then
               write(lgfile,ifmt) name, n, idim, k
            else
               write(lgfile,ifmt) name, n, idim
            endif
         endif
      endif
      if(.not. rflag) call abend
      return
1001  format(7x, 31h parameter error in subroutine , a6,
     &         27h ... (n.lt.1 or idim.lt.n).,
     &         5x, 3hn =, i4, 5x, 6hidim =, i4, 1h. )
1002  format(7x, 31h parameter error in subroutine , a6,
     &         37h ... (n.lt.1 or idim.lt.n or k.lt.1).,
     &         5x, 3hn =, i4, 5x, 6hidim =, i4, 5x, 3hk =, i4,1h.)
      end
      subroutine lfit(x,y,l,key,a,b,e)
C-----------------------------------------------------------------------
C
C     TO FIT A STRAIGHT LINE    Y=A*X+B    TO L POINTS WITH ERROR E
C     SEE MENZEL , FORMULAS OF PHYSICS P.116
C     POINTS WITH Y=0 ARE IGNOERD IF KEY=0
C     L IS NO. OF POINTS
C
C-----------------------------------------------------------------------
      implicit none
      integer j,key,l
      real a,b,count,e,scartx,scarty,sumx,sumxx,sumxy,sumy,sumyy,x,xmed,
     &y,ymed
      dimension x(1),y(1)
!$OMP THREADPRIVATE(j,
!$OMP& count,scartx,scarty,sumx,sumxx,sumxy,sumy,sumyy,xmed,ymed)
      save

C-----------------------------------------------------------------------
C
C     CALCULATE SUMS
C
C-----------------------------------------------------------------------
      if(l-2) 25,1,1
    1 count=0.0
      sumx=0.0
      sumy=0.0
      sumxy=0.0
      sumxx=0.0
      sumyy=0.0
      do 10 j=1,l
      if(y(j).eq.0..and.key.eq.0) goto 10
      sumx=sumx+x(j)
      sumy=sumy+y(j)
      count=count+1.0
   10 continue
      if(count.le.1.) goto 25
      ymed=sumy/count
      xmed=sumx/count
      do 20 j=1,l
      if(y(j).eq.0..and.key.eq.0) goto 20
      scartx=x(j)-xmed
      scarty=y(j)-ymed
      sumxy=sumxy+scartx   *scarty
      sumxx=sumxx+scartx   *scartx
      sumyy=sumyy+scarty   *scarty
   20 continue
C
C     FIT PARAMETERS
      if(sumxx.eq.0.) goto 25
      a=sumxy/sumxx
      b=ymed-a*xmed
      if(count.lt.3.) goto 101
      e=(sumyy-sumxy*a          )/(count-2.0)
      goto 100
C
C     ISUFFICIENT POINTS
   25 a=0.0
      b=0.0
  101 e=0.0
  100 return
      end
      subroutine lfitw(x,y,w,l,key,a,b,e)
C-----------------------------------------------------------------------
C
C     TO PERFORM A WEIGHTED STRAIGHT LINE FIT
C
C     FOR FORMULAE USED SEE MENZEL, FORMULAS OF PHYSICS P.116
C
C     FIT IS OF Y=AX+B , WITH S**2 ESTIMATOR E. WEIGHTS ARE IN W.
C     IF KEY=0, POINTS WITH Y=0 ARE IGNORED
C     L IS NO. OF POINTS
C
C-----------------------------------------------------------------------
      implicit none
      integer icnt,j,key,l
      real a,b,e,w,w2,w2x,w2x2,w2xy,w2y,w2y2,ww,wwf,wwfi,x,y
      dimension x(1),y(1),w(1)
!$OMP THREADPRIVATE(icnt,j,
!$OMP& w2,w2x,w2x2,w2xy,w2y,w2y2,ww,wwf,wwfi)
      save

C-----------------------------------------------------------------------
C
C     CALCULATE SUMS
C
C-----------------------------------------------------------------------
      if(l.le.1) goto 1
      w2=0.
      w2x=0.
      w2y=0.
      w2xy=0.
      w2x2=0.
      w2y2=0.
      icnt=0
      do 2 j=1,l
      if(y(j).eq.0..and.key.eq.0) goto 2
      ww=w(j)*w(j)
      w2=ww+w2
      wwf=ww*x(j)
      w2x=wwf+w2x
      w2x2=wwf*x(j)+w2x2
      w2xy=wwf*y(j)+w2xy
      wwfi=ww*y(j)
      w2y=wwfi+w2y
      w2y2=wwfi*y(j)+w2y2
      icnt=icnt+1
    2 continue
C
C     FIT PARAMETERS
      a=(w2xy-w2x*w2y/w2)/(w2x2-w2x**2/w2)
      b=(w2y-a*w2x)/w2
      if(icnt.le.2) goto 3
      e=(w2y2-w2y**2/w2-(w2xy-w2x*w2y/w2)**2/(w2x2-w2x**2/w2))/(icnt-2)
      goto 4
C
C     ISUFFICIENT POINTS
    1 a=0.
      b=0.
    3 e=0.
    4 return
      end
      subroutine abend
      print*,'abend called ==> problem'
      return
      end
*
* $Id: cfft.F,v 1.1.1.1 1996/02/15 17:48:48 mclareni Exp $
*
* $Log: cfft.F,v $
* Revision 1.1.1.1  1996/02/15 17:48:48  mclareni
* Kernlib
*
*
*#include "pilot.h"
      SUBROUTINE CFFT(A,MSIGN)
      COMPLEX A(1),U,W,T
      save
 !$OMP THREADPRIVATE(m,n,msign,nv2,nm1,i,j,k,le,le1,l,ip,c,s,U,W,T)
      IF(MSIGN.EQ.0) RETURN
      M=IABS(MSIGN)
      N=2**M
      NV2=N/2
      NM1=N-1
      J=1
      DO 7 I=1,NM1
      IF(I.GE.J) GO TO 5
      T=A(J)
      A(J)=A(I)
      A(I)=T
 5    K=NV2
 6    IF(K.GE.J) GO TO 7
      J=J-K
      K=K/2
      GO TO 6
 7    J=J+K
      DO 8 I=1,N,2
      T=A(I+1)
      A(I+1)=A(I)-T
 8    A(I )=A(I)+T
      IF(M.EQ.1) RETURN
      C=0.
      S=ISIGN(1,MSIGN)
      LE=2
      DO 20 L=2,M
      W=CMPLX(C,S)
      U=W
      C=SQRT(C*.5+.5)
      S=AIMAG(W)/(C+C)
      LE1=LE
      LE=LE1+LE1
      DO 9 I=1,N,LE
      IP=I+LE1
      T=A(IP)
      A(IP)=A(I)-T
 9    A(I) =A(I)+T
      DO 20 J=2,LE1
      DO 10 I=J,N,LE
      IP=I+LE1
      T=A(IP)*U
      A(IP)=A(I)-T
 10   A(I) =A(I)+T
 20   U=U*W
      RETURN
      END
