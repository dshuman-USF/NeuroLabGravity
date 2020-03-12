      subroutine qcalc (IERR,ICFLAG,sig)
c     subroutine GMINTU4(IERR)<-- derived from g's version
c     V1.01  GLG  READS *.ASC OR *.REN (I2,F10.1) FILE; CALCULATES   
c     RUNNING MEAN INTERVAL OVER 20 SPIKES FOR EACH ID;  
c     WRITES UNFORMATTED 'INPUT.GRA': ID, EVTIME, AVINT  
c     INCLUDES HEADER: IFNAME(14), IHS    
c     v1.02  GLG VERSION FOR DG   
c     v1.03  GLG 20-APR-84 ADDS TO HEADER: IFNAME, IHS, ETIME   
c     v1.04 GLG 12-jun-84 subroutine;no output files;calculates forward   
c     and backward exponential charge histories for 10 neurons;
c     output as large COMMON array CHIST(10000,20); first N    
c     columns old (later) exp; next N columns backwards charge;
c     Normalization from MINT    
c     v1.041 GLG 25-feb-85 deals with edited arrays; a) eliminates dead   
c     time between edited chunks; b) translates IDs 3xx and 4x 
c     into your choice of particle codes 1-10. These are kept  
c     as IDSET(10) in IBLOCK of GHEADER.DEF. CHIST to 12000    
c     v1.042 GLG 8-aug-85 CHIST split into QE and QAFS, each N wide. 
c     Variable N (typein from main) all in GHEADER2.DEF.  
c     v2.01  GLG 30-aug-85 Universal version like 42, but allows all 
c     three possible types of normalization: 1) FWDINC=MEANINT 
c     2) FWDTAU=MEANINT/ADJF; 3) NONE. Must be called by  
c     GFAST3 or later, so that all relevant parameters are
c     selected.   
c     v2.02  GLG 19-sep-85 calculates meanint over +/- IHS only if more   
c     than 30 spikes; otherwise MEANINT = (LastT-firstT)/number.    
c     Charges defined to ENDTIME at end of run, both for forward    
c     and backward decays.  
c     v2.03  GLG 11-nov-85 keeps stimuli at constant (typin) spacing as   
c     it concatenates edited chunks. Does not eliminate silent 
c     stimulus periods, or after end of edited chunk.
c     16-dec-85   
c     In NORM=3 uses local rather than global MEANINT. The global   
c     value calculated from TIMAT rather than ENDTIME is used only  
c     if local MEANINT is zero (as at end of TIMAT). 
c     Allows last event to be MARKE.  
c     v3.01  GLG 3-feb-86 installs PST histogram for stimulus control
c     charge. Resolution 1 MS; PST tables maximum length of    
c     1000, corresponding to 1 sec between stimuli. Actual
c     duration calculated from (1),(2) event time lag.    
c     V3.02  GLG 11-nov-86 MPI VAX version, pst charge corrections.  
c     v3.03  GLG 18-dec-87 SUN4 version. Fixes concatenation problems
c     in SECOND pass. Be SURE to give stim period that
c     is larger than any expected period!!!
c     v. gmint8  BGL add lines to allow FS =-1 and no PST convol. corr.
c     v  gmint9 skips pst charge corr.. with c'd code at various locations
c     ...also limit of time chunks to 1000.ms removed as a consequence!
c     ... also c'd is poss code for
c     pst correction for back as well as fwd charge with norm = 3
c     logic not confirmed yet so not implemented 
c     
c     
c     v10. change to qe qasf 120000,10
c     v10a. makes 100 cycles for significance tests
c     v32. allows up to 100 cycles for significance tests and 32 neurons
c     
c     
c     v bg matches mods in gbatchbg - 10000 codes, 1000 steps, only 16 units
c     v2000  matches mods in gbatch2000  - see that code for changes too
c     v2kv14  12-16-03  read in i3 integers as spike codes
c     may 04 - bgl to match gbatch_04 changes - including handling of 64 codes
c     changes in def file in ghead_04 ...
c     parameter (ncl=64)!number of codes handled
c     parameter (nsz=480000)!size of charge history array - 1st dimension
c     parameter (nsppc=10000)!max number of spikes in the busiest channel
c     parameter (lstsz=640000)!Product of nsz*nsppc
c     
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      INCLUDE 'ghead.def'  
c------------------------------------------------------------------------
     
      REAL MINT,GMINT(ncl)
      INTEGER NSEQ(ncl),NSEQL(ncl),stat
      INTEGER IDTRANS(999),ICFLAG,ETEMPL
      logical sig
      CHARACTER*1 CHTEM  
      character*30 garbage
      dimension iidlist(lstsz),eevtlist(lstsz)
      character(:),allocatable :: fmt
     
c     **********ARRAY STRUCTURES**********  
c     IDLIST,EVTLIST raw edited data from I3,I8 file  
c     TIMAT(EVTIME,1->N; MINT,N+1->2N) 
c     QE(fwd or bak)(I,1->N); QAFS(fwd or bak)(I,1->N)
c     ************************************  
     
      if (.not.allocated(fmt)) allocate (fmt,source='') ! work around gcc bug 60500
      ITIMATL=nsppc             !matches timat in header
      IQELEN = nsz    
c     !matches QE and QAFS definition in GHEADER.DEF   
      IERR = 0 
      BUFTAU = 4.0*FWDTAU     
c     !buffer time between edited chunks
      MARKE = 22      
c     !end mark for edited chunk   
      MARKB = 21      
c     !begin mark for edited chunk 
     
c     *****SET UP ID TRANSFORMATION; IDTRANS CONTAINS 0's except at  
c     64 values which are desired particle number.   
c     ***initialize IDTRANS to zero***   
     
      if (ICFLAG.LE.1) then
         DO I=1,999    
            IDTRANS(I) = 0
         end do
         
c     ***KEYBOARD ASSIGNMENT OF 'ncl' ID's***    
         CHTEM = 'n'
         do while (scan(chtem,'yY').eq.0)
            DO I=1,N  
               PRINT "(1X,'ID CODE FOR PARTICLE ',I3,' IS (I3): ',$)", I 
               read (*,'(I3)') IDTEMP    
               IDTRANS(IDTEMP) = I
               IDSET(I) = IDTEMP  
            end do
            
c     verification of assignments 
            print *, 'verify your assignments: '   
            DO I=1,N 
               print "(1x,'particle ',I3,' is ID# ',I3)", I,IDSET(I)  
            end do
            
            print *, ' is this correct? (Y/N)'
            read (*,'(A)') CHTEM
         end do
         
c     *****SET UP INPUT FILE***** 
         if (sig) then
            read (*,'(A)') garbage
         else
            print "(1X,'INPUT DATA FILE NAME IS: ',$)"
            read (*,'(A)') IFNAME
         end if
      end if
      if (sig) write (IFNAME, "('sh',I0,'.rdt')") ICFLAG
     
c     open (1,file=IFNAME,status='old',form='FORMATTED',readonly)
      open (1,file=IFNAME,status='old',form='FORMATTED')
     
c     *****FIRST PASS: READ INPUT FILE***** 
      IDLIST = 0 
      EVTLIST = 0

      itempl = 0
      do while (itempl .eq. 0) ! loop until non-blank record
         READ (1,'(I5,I8)',iostat=stat) ITEMPL,ETEMPL !read .BDT format
         if (stat .ne. 0) exit
      end do
      if (stat .eq. 0 .and. ETEMPL.eq.1111111) then !if *.BDT type file then..
         fmt = '(I5,I8)'
         READ (1,'(I5,I8)',iostat=stat) ITEMPL,ETEMPL !read .BDT format
      else
         fmt = '(I2,I8)'
         rewind 1               !.ADT file;go back to start
      end if
      if (stat .eq. 0) then
         DO I=1,lstsz-1
            itemp1 = 0
            do while ((ITEMP1.eq.0).or.(ITEMP1.gt.999))
               READ (1,fmt,iostat=stat) ITEMP1,ITEMP2 ! adt and old gdt format
               if (stat .ne. 0) exit
            end do
            if (stat .ne. 0) exit
            IDLIST(I) = ITEMP1
            if (ITEMP2.EQ.0) ITEMP2 = 2 ! deals with 0 initial time
            EVTLIST(I) = 0.5*ITEMP2     ! conversion from clock tick to MS
            iidlist(I) = idlist(i)      ! save original values for call to icel_chg
            eevtlist(I) = evtlist(i)    ! save original values for call to icel_chg
         end do
      end if
      close (unit=1)
      LAST = I - 1              ! last entry index in IDLIST,EVTLIST array
     
      if (ICFLAG.LE.1) then
         print "(1X,' define TIME (ms) between stimuli (F10.1): ',$)"
         read (*,'(F10.1)') STIMPER ! stimulus period 
      end if
          
      TIMADJ = EVTLIST(1) -0.5  !! bgl mod to prevent 0.0 in 1st TIMAT
c     ...happens sometimes otherwise
     
      LISTDX = 1    
      NSTIM = 0
c     *****SECOND PASS: Concatenate time chunks, transform IDs  
c     *****FIND NEXT MARKE, CHECK FOR MARKB, ADJUST TIMES*****  

      do while (LISTDX.LE.LAST)
         IF (IDLIST(LISTDX).NE.MARKE) THEN  
            EVTLIST(LISTDX) = EVTLIST(LISTDX) - TIMADJ
         ELSE IF ((LISTDX.LT.LAST).AND.(IDLIST(LISTDX+1).NE.MARKB)) THEN
            print *,'**no MARKB after MARKE in file: FATAL ERROR**'
            print *, LISTDX,IDLIST(LISTDX),EVTLIST(LISTDX)
            STOP    
         ELSE  
            EVTLIST(LISTDX) = EVTLIST(LISTDX) - TIMADJ
            IF (LISTDX.GT.1) THEN  
c     EVTLIST(LISTDX) = EVTLIST(LISTDX-1) 
               IF (IDLIST(LISTDX-1).NE.MARKB) THEN   
                  NSTIM = NSTIM + 1    
               END IF  
            ELSE
               EVTLIST(LISTDX) = 0
            END IF  
            TIMADJ = EVTLIST(LISTDX+1) -  NSTIM*STIMPER
         END IF   
         LISTDX = LISTDX + 1
      end do
     
c     *****TRANSFORM ID's to particles 1 to N as in IDTRANS array*****    
      DO I = 1,LAST 
         ITEMP3 = IDLIST(I) 
         IF (ITEMP3.GT.999) THEN 
            PRINT *,'Illegal ID .GT. 999; FATAL ERROR' 
            STOP    
         END IF   
         IDLIST(I) = IDTRANS(ITEMP3)  
c     print *, I, IDLIST(I), EVTLIST(I)       
c     !diagnostic   
      end do
     
c     *****THIRD PASS: Load TIMAT from transformed IDLIST and   
c     concatenated EVTLIST. 
     
      DO I=1,N  
         NSEQ(I) = 1   
      end do
     
      DO J=1,2*N
         DO I=1,ITIMATL
            TIMAT(I,J) = 0.0   
         end do
      end do
     
      DO I=1,LAST 
         ID = IDLIST(I)
         IF (ID.EQ.0) cycle
         TIMAT(NSEQ(ID),ID) = EVTLIST(I)   
         IF (NSEQ(ID).GE.ITIMATL) THEN   
            PRINT "(1X,'MORE THAN ',I0,' OF ID ',I3,'-- FILES TRUNCATED')",nsppc,ID 
            print *,'Fix your spike file'
            STOP
         END IF   
         NSEQ(ID) = NSEQ(ID) + 1 
      end do
     
      ENDTIME = EVTLIST(LAST) 
      print "(1x,'TIMAT LOADED; ENDTIME: ',F10.1)",ENDTIME
     
c     *****FOURTH PASS: CALC MINT, PUT INTO RIGHT COLS TIMAT*****    
c     *****MINT CALCULATED OVER +/- IHS SPIKES OF EACH ID*****
c     IHS = 10 
      DO J=1,N  
         NSEQ(J) = 0   
      end do
c     
c     *****FIND LAST ENTRY IN EACH ID COLUMN*****
      DO J=1,N  
         DO I=1,ITIMATL
            NSEQL(J) = I-1
            IF (TIMAT(I,J).EQ.0.0) exit
         end do
         NSTEP = int(TIMAT(NSEQL(J),J)/DELTIM)
         IF (NSTEP.LE.IQELEN) cycle    
         IERR = 1 
         PRINT *,'TIME STEP TOO SMALL FOR ARRAY SIZE; TRY AGAIN'    
         RETURN
      end do
      print *,'event counts (under respective ID) are:'
      fmt = merge ('(1X,10I4)', '(1X,32I8)', sig)
      print fmt,(IDSET(ID),ID=1,N) 
      print fmt,(NSEQL(ID),ID=1,N) 
     
c     ********bgl mod: autocalc of istpfr,lfram*****************
      MAXSTP= int(ENDTIME/DELTIM)
      LFRAM=1000
      if (MAXSTP.LE.1000)LFRAM=MAXSTP
      ISTPFR=MAXSTP/LFRAM
          
c     changes in def file in ghead_04 ...see for current values
c     parameter (ncl=64)!number of codes handled
c     parameter (nsz=480000)!size of charge history array - 1st dimension
c     parameter (nsppc=10000)!max number of spikes in the busiest channel
c     parameter (lstsz=640000)!Product of nsz*nsppc
     
c     *****CALCULATE MEAN INTERVAL*****
      DO ID=1,N 
         IF (NSEQL(ID).GT.30) THEN    
            DO I=1,NSEQL(ID)
               IF (I.GE.nsppc-1) THEN      
c     !truncation message  
c     PRINT 981,ID 
                  continue      !bgl mod of g's code...see different
c     approach to truncation around label 125
               ELSE IF (I-IHS.LT.1) THEN   
                  TIMAT(I,ID+N) = (TIMAT((2*IHS)+1,ID)-TIMAT(1,ID))/(2*IHS)   
c     Array element changed to (2*IHS)+1 so number of intervals really 2*IHS
c     Added by SCN  11/98
               ELSE IF (I+IHS.GE.NSEQL(ID)) THEN
                  TIMAT(I,ID+N) = (TIMAT(NSEQL(ID),ID)-TIMAT(NSEQL(ID)-   
     X                 (2*IHS),ID))/(2*IHS)
               ELSE    
                  TIMAT(I,ID+N) = (TIMAT(I+IHS,ID)-TIMAT(I-IHS,ID))/(2*IHS)    
               END IF    
            end do
         ELSE   
            MINT = (TIMAT(NSEQL(ID),ID) - TIMAT(1,ID))/(NSEQL(ID)-1)    
            DO I=1,NSEQL(ID)
               TIMAT(I,ID+N) = MINT
            end do
         END IF   
c     ###save global mean intervals to GMINT###  
         GMINT(ID) = (TIMAT(NSEQL(ID),ID) - TIMAT(1,ID))/(NSEQL(ID)-1)   
      end do
c     
c     *****CLEAR CHARGE HISTORY ARRAYS QE and QAFS*****    
      DO J=1,N    
         DO I=1,IQELEN  
            QE(I,J) = 0.0 
            QAFS(I,J) = 0.0    
         end do
      end do
     
      MEANINT = 0
      IF (ISHIFTE.EQ.0) THEN  
c     *****EFFECTOR FORWARD CHARGE HISTORIES *****  
c     +++use GLOBAL meanint: T/N or LOCAL meanint from TIMAT+++   
         DO ID=1,N   
            MEANINT = int(TIMAT(1,ID+N))
c     !for local meanint    
            MEANINTG = int((TIMAT(NSEQL(ID),ID))/NSEQL(ID))
c     !for global meanint    
            K = 1 
            EVTIME = TIMAT(K,ID) 
c     print *,'evtime = ',evtime
c     do 740 ick=1,100 
c     print 3100, TIMAT(ick,1),TIMAT(ick,9)
c     740       CONTINUE
c     
            EVTIMEL = TIMAT(NSEQL(ID),ID) + 2*FWDTAU 
c     print *,'evtimel =',evtimel
            FWDQ = 0.0 
            DO I=1,IQELEN    
               TIME = I*DELTIM 
               IF (TIME.GE.EVTIME) THEN  
                  MEANINT = int(TIMAT(K,ID+N))
c     !for local meanint
c     print *,'I,K,EVTIME,TIME,DELTIM = ',I,K,EVTIME,TIME,DELTIM
                  FWDQ = FWDQ*FWDDEC + FWDINC*EXP((EVTIME-TIME)/FWDTAU)    
c     print *,'time = ',TIME,' fwdq = ',fwdq
                  K = K + 1   
                  EVTIME = TIMAT(K,ID)  
                  IF (EVTIME.EQ.0.0) EVTIME = ENDTIME  
               ELSE 
                  FWDQ = FWDQ*FWDDEC    
c     print *,'FWDQ = ',FWDQ

               END IF
               IF (MEANINT.EQ.0.0) meanint = meanintg
               IF (MEANINT.EQ.0.0) then
                  print *,"Particle code ",idset(id),
     +                 " may have too few spikes."
                  stop
               else
                  QE(I,ID) = FWDQ - (FWDINC*FWDTAU/MEANINT)  
c     print *,'QE(I,ID),I,FWDQ,MEANINT ',QE(I,ID),I,FWDQ,MEANINT
               end if
c     
               IF (TIME.GT.ENDTIME) exit
c     !charge defined to end of run 
            end do
         end do
     
c     do 770 ick=1,1000 
c     print 3100, QE(ick,1)
c     print *,'meanint ',MEANINT
c77   0       CONTINUE
c     
         PRINT "(1X,'EFFECTOR FORWARD CHARGE HISTORIES COMPLETED')"
c     
      ELSE
c     *****EFFECTOR BACKWARD CHARGE HISTORIES ***** 
c     +++use GLOBAL meanint: T/N or LOCAL meanint from TIMAT+++   
         DO ID=1,N   
            MEANINTG = int((TIMAT(NSEQL(ID),ID))/NSEQL(ID))
c     !for global meanint    
            K = NSEQL(ID)   
            EVTIME = TIMAT(K,ID) 
            IL = int(ENDTIME/DELTIM)
c     !charge defined to end of run  
            BAKQ = 0.0 
            DO I=IL,1,-1
               TIME = I*DELTIM 
               IF (TIME.LE.EVTIME) THEN  
                  MEANINT = int(TIMAT(K,ID+N))
c     !for local meanint
                  BAKQ = BAKQ*BAKDEC + BAKINC*EXP((TIME-EVTIME)/BAKTAU)    
                  K = K - 1   
                  IF (K.LE.0) THEN 
                     EVTIME=0.0   
                  ELSE 
                     EVTIME = TIMAT(K,ID)   
                  END IF 
               ELSE
                  BAKQ = BAKQ*BAKDEC    
               END IF
               IF (MEANINT.EQ.0.0) meanint = meanintg
               IF (MEANINT.EQ.0.0) then
                  print *,"Particle code ",idset(id),
     +                 " may have too few spikes."
                  stop
               else
                  QE(I,ID) = BAKQ - (BAKINC*BAKTAU/MEANINT)  
               end if
            end do
         end do
         PRINT "(1X,'EFFECTOR BACKWARD CHARGE HISTORIES COMPLETED')"
      END IF     
     
      IF (ISHIFTA.EQ.0) THEN  
c     *****(ACCEPTOR FORWARD CHARGE)*FS HISTORIES *****  
c     +++use GLOBAL meanint: T/N or LOCAL meanint from TIMAT+++   
         DO ID=1,N  
            MEANINT = int(TIMAT(1,ID+N))
c     !for local meanint    
            MEANINTG = int((TIMAT(NSEQL(ID),ID))/NSEQL(ID))
c     !for global meanint    
            K = 1 
            EVTIME = TIMAT(K,ID) 
            EVTIMEL = TIMAT(NSEQL(ID),ID) + 2*FWDTAU 
            FWDQ = 0.0 
            DO I=1,IQELEN   
               TIME = I*DELTIM 
               IF (TIME.GE.EVTIME) THEN  
                  MEANINT = int(TIMAT(K,ID+N))
c     !for local meanint
                  FWDQ = FWDQ*FWDDEC + FWDINC*EXP((EVTIME-TIME)/FWDTAU)    
                  K = K + 1   
                  EVTIME = TIMAT(K,ID)  
                  IF (EVTIME.EQ.0.0) EVTIME = ENDTIME  
               ELSE
                  FWDQ = FWDQ*FWDDEC    
               END IF
               IF (MEANINT.EQ.0.0) meanint = meanintg
               IF (MEANINT.EQ.0.0) then
                  print *,"Particle code ",idset(id),
     +                 " may have too few spikes."
                  stop
               else
                  QAFS(I,ID) = (FWDQ - (FWDINC*FWDTAU/MEANINT))*FS    
               end if
               IF (TIME.GT.ENDTIME) exit      
c     !charge defined to end of run
            end do
         end do
         PRINT "(1X,'(ACCEPTOR FORWARD CHARGE)*FS HISTORIES COMPLETED')"
      ELSE
c     *****(ACCEPTOR BACKWARD CHARGE)*FS HISTORIES ***** 
c     +++use GLOBAL meanint: T/N or LOCAL meanint from TIMAT+++   
         DO ID=1,N  
            MEANINTG = int((TIMAT(NSEQL(ID),ID))/NSEQL(ID))
c     !for global meanint    
            K = NSEQL(ID)   
            EVTIME = TIMAT(K,ID) 
            IL = int(ENDTIME/DELTIM)
c     !charge defined to end of run  
            BAKQ = 0.0 
            DO I=IL,1,-1    
               TIME = int(I*DELTIM)
               if (sig) TIME = int(time) ! to preserve version 0.10.3 gmint/gmintsig difference
! I don't know that it is necessary - ROC
               IF (TIME.LE.EVTIME) THEN  
                  MEANINT = int(TIMAT(K,ID+N))
c     !for local meanint
                  BAKQ = BAKQ*BAKDEC + BAKINC*EXP((TIME-EVTIME)/BAKTAU)    
                  K = K - 1   
                  IF (K.LE.0) THEN 
                     EVTIME=0.0   
                  ELSE 
                     EVTIME = TIMAT(K,ID)   
                  END IF 
               ELSE
                  BAKQ = BAKQ*BAKDEC    
               END IF
               IF (MEANINT.EQ.0.0) meanint = meanintg
               IF (MEANINT.EQ.0.0) then
                  print *,"Particle code ",idset(id),
     +                 " may have too few spikes."
                  stop
               else
                  QAFS(I,ID) = (BAKQ - (BAKINC*BAKTAU/MEANINT))*FS     
               end if
c     
            end do
         end do
         PRINT "(1X,'(ACCEPTOR BACKWARD CHARGE)*FS HISTORIES COMPLETED')"
      END IF    
c     
c     
c     Subroutine to optionally include intracellular channel
      if(ICFLAG.le.1) then
         print *,"Calling icel_win;N= ",N
         mult_ch = 0
         call icel_win(iqelen,stimper,fwdinc,mult_ch,iidlist,eevtlist)
         qeall1 = 0.
         qeall9 = 0.
         print *,"icel_win complete,N= ",N
c     print *,'QE then QAFS'
c     do 980 ick=1,2290
c     qeall1 = qeall1 + QE(ick,1)
c     qeall9 = qeall9 + QE(ick,9)
c     print 3051,ick,(QE(ick,ID),ID=1,N)
c     print *,qeall1,'   ',qeall9
c     980       CONTINUE
      end if
      print "(1X,'NORMAL EXIT; START GRAVITY ')"
      CLOSE (UNIT=1)
c     
      RETURN   
c     
      END 


