      PROGRAM gbatch_08_v1
c***************************************************************************
c       -GRAVITATIONAL CLUSTERING
c       V1.01  GLG  09-APR-82
c       V1.02  GLG  11-APR-82 TYPES POSIT(I,J)
c       V1.03  GLG  15-APR-82 NO FORCE CORE, TYPES PDISQ(I,J)
c       V2.01  GLG  23-APR-82 DISK DATA FILE, RED CHARGE
c       V2.02  GLG  24-APR-82 CHARGE HISTORY THRESH MATRIX
c       V2.03  GLG  27-APR-82 PRINTS FINAL STATE +CHG HIST MTX 
c       V2.04  GLG  28-APR-82 TYPE IN STEP TIME,CHG CONSTANTS
c       V2.05  GLG  28-APR-82 AUTOSIZE NO FORCE CORE
c       V3.01  GLG  03-JUN-82 NO SP K;CHGM=T;EFCHG=CHG-TAU;TYP NFC
c       V4.01  GLG  12-APR-83 USES DATA IN ASCII 'REN' FORMAT
c       V4.02  GLG  14-APR-83 READY FOR SUBROUTINE STRUCTURE
c       V4.03  GLG  15-APR-83 MAKES UNFORMTD OUT FILE DS AT 40ST INT
c       V6.01  GLG  12-JUN-83 KEEPS 45 DS IN UNFORMATD OUTFILE
c       V6.02  GLG  16-JUN-83 INIT DIST 100; IMPROVES NO SPIKE BEH
c       V6.03  GLG  22-AUG-83 SAVES POSITIONS INSTEAD 45 DS; UNFRMTD OUTF
c             OUTPUT EVERY 80 STEPS (RATHER THAN 40)
c       V6.04N GLG  08-DEC-83 INVERTS SIGN OF FORCE PRODUCT FOR INHIBITION
c             INPUT REN FILE IN FORMAT I2,F10.1
c       V6.05N GLG  28-DEC-83 +/- IHS SPIKES RUNNING AVERAGE FOR REDQMs
c             USES UNFMTD 'INPUT.GRA' FILE: ID, TIME, AVREDQM(ID)
c             HEADER: IFNAME(14), IHS. PRODUCED BY 'MINT1'.
c       V6.06N GLG  13-FEB-84 SAME TRANSFERRED TO DG; MODS IN lower case.
c             tests for inactive neurons near 330 are C'd out.
c       V6.06P GLG  28-FEB-84 DG version; force product made POSITIVE
c       V6.07P GLG  13-APR-84 MAKES STEPS/FRAME AND # OF FRAMES TYPEIN
c             adds header block to output file; NOTE HEADER BLK DEFNTN
c       V6.08P GLG  12-JUN-84 forward and backward charge histories in
c             array from subroutine GMINT4.
c       V7.01P GLG  26-FEB-85 uses subroutine GMINT41 for edited real data.
c             Assumes I3,I8 file of ID,TICK; concatenates chunks; 
c             reassigns IDs to be 1-10.
c       V7.01PM GLG 29-mar-85 installs FLASH array for movies.  Note that
c             IDLIST and EVTLIST are moved to GHEADER.DEF
c       V1.01M GLG 8-AUG-85 first GFAST version; improves basic loop
c             structures in calculation for speed. Uses GMINT42;
c             CHIST divided into QE and QAFS for less array calculation.
c             Queries about acceptor/effector charge time sense and
c             force sign (exit/inhib). Uses GHEADER2.
c        V1.02M GLG 22-aug-85 removes various diagnostic printouts.
c             Variable N typed in, stored in GHEADER2.DEF.
c        V1.03M GLG 6-sep-85 Universal version, allows choice of three
c             NORMalizations (1,2,3) <=> INC (original), TAU (old 'j'),
c             or none (old 's').
c             Second step normalization for mean always done.
c             Note that type in varies with choice of NORM!!!
c             Fix typos in charge time sense queries. 10-sep-85.
c       V1.04M GLG 8-nov-85 uses GMINTU4.
c       V1.04M GLG 8-nov-85 uses GMI, longer and narrower arrays.
c             Checks for N max =6.
c       V2.0  BGL 28-MAY-86 SMALLER & VIRTUAL ARRAYS USED FOR
c               11/23 VERSION. CHANGES MADE IN GHEAD3.F77 (DEFINITIONS)
c
c       V2.01 BGL 11-AUG-86 SOME PRINTS CHANGED TO TYPES...
c       NOW SOME OUTPUTS TO TERMINAL AND SOME TO PRINTER.
c       V2.2  BGL 2-SEP-86 MODIFY SOME LOGICAL*1 I/O..FORMERLY CHARACTER
c       TYPE
c       v4.0 first version for hp350: 13 cells, large arrays, code like g's
c       v4.1 22-jul-87 minor changes :pgm optionally loops to start 
c       v6.0 30-oct-87 mods so posdsx and trydsx are subroutines to 
c       gfast6 - MUST be run on main HP console (because of graphics)
c
c       v7.0 27-jan-88 mods to run with g's new gmint-our gmint7
c       v7.01 30-mar-88 bgl change to #8 to match gmint8
c       v7.02 29-jun-88 bgl change to #9 to match gmint9
c       gfast10.f converted to gbatch10.f
c       v10.1 allows mod of IHS at run time, def is 10
c       
c       3/16/94 expand to 32 neuron capability bgl gbatch32.f
c
c       big version 1 12-26-96
c       increase qe,qafs to 480000 , spikes to 10000/cell/ 1000 steps
c       but num cells decreased to 16, in arrays now at 160000.
c
c       mod from gbatchsnw to gbatch2000 by bgl 12-20-99  - 
c       1- now only have old norm 3 - no rate normalization
c       2- flash write now only dummy values; never used at USF anyway..
c       ..... simply removes some overhead associated with code
c            
c       
c       OPTION ('N') STANDARD GRAVITY- "N-DIMENSIONAL" - NO TUNING
c       OPTION ('O') "N-DIMENSIONAL" GRAVITY - TUNING
c       OPTION ('P') "ONE DIMENSIONAL" GRAVITY - NO TUNING
c       OPTION ('Q') "ONE DIMENSIONAL" GRAVITY - TUNING

c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      implicit none
      include 'config.defs'
      INCLUDE 'ghead.def'
c-----------------------------------------------------------------
c
      REAL POSIT(ncl,ncl),PFDIR(ncl,ncl,ncl)
c       next line for v12 and beyond - pair only sub option 
      DOUBLE PRECISION PDIS2(ncl,ncl)
      REAL PDIS2F(ncl,ncl)
      REAL FLASH(ncl), TEMPD(ncl)
      INTEGER ipofset(ncl,ncl)
      character*1 instr
      character*30 posfile
      character*30 garbage
      real chpr1,chpr2,dist,flashlen,pdisq,temp,temp1,time,totfk,trans
      integer i,ibuff,iend,ierr,ifram,ioptop,irange,istart,istep,ita
      integer itemp1,itemp2,j,jta,k,m,nmax,stat,slash,len,icflag,i_savN
      integer I_N
      CHARACTER(:), ALLOCATABLE :: path
      CHARACTER(:), ALLOCATABLE :: basename
      logical sig, batch

      if(.false.)print*,fildes
      posfile = ' '
      garbage = ' '
      call get_command_argument (0, length=len, status=stat)
      if (stat.ne.0) stop "Can't get command path length"
      allocate (character(len) :: path)
      call get_command_argument (0, value=path, status=stat)
      if (stat.ne.0) stop "Can't get command path"
      slash = index (path, "/", back=.TRUE.)
      allocate (character(len - slash) :: basename)
      basename = path ((slash + 1) : len)
      progid = basename
      sig = .false.
      batch = .false.
      sig = progid .eq. 'gsig' .or. progid .eq. "gsig_05"
      batch = progid .eq. 'gbatch'

c     TIMAT only used in subroutines but in common header...
c     write to it to avoid compilation warning...
      TIMAT (1,1)=0.0           !simply
      ICFLAG=1

      if (batch) then
         instr = 'y'
         do while (scan(instr,'yY').ne.0)
            call gbatch
            print "('enter ""y"" to continue-any other char to exit..')"
            read (*,'(A)') instr
         end do
      else if (sig) then
         read (*,'(I12)') nshft
         read (*,'(A)') garbage
         do icflag=1,nshft
            call gbatch
            N = i_savN
         end do
      else
         print "(A,': ',$)", basename
         stop "due to unknown program name"
      end if

      contains

c     ****************************************************************
      SUBROUTINE GBATCH
c     ****************************************************************
      integer trycnt
c       *****BASIC CONSTANTS and PARAMETERS*****     
      NMAX = ncl                !parameter max. number of cells allowed
      if (icflag .eq. 1) then
         N= ncl    
         I_N=ncl     
         ISHIFTA = 0            !these three are reset by queries.     
         ISHIFTE = 0     
         FS = 1.0               !force sign: +1=excitation; -1=inhibition     
         FLASHLEN = 3.0         !# of frames unit will flash ! for backward comp
      end if
c       three variables used in new optimizer routine defined here
c       **POSIT(PARTICLE #,COORDINATE#)**
c       **PFDIR(ME, THE OTHERS, COORDINATE COMPONENTS)**
c       **PDISQ(ME, THE OTHERS)**
c
c       *****INITIALIZES POSITION MATRIX OF N PARTICLES*****
c       +++++CLEARS DIAGONAL IN PFDIR+++++     
      PFDIR = 0.0
      POSIT = 0.0
      DO I=1,ncl
         POSIT(I,I) = 50.0*SQRT(2.0)     
      end do

c     fill flash array with dummy values for backward compatability
c     at flash write
      FLASH = FLASHLEN
      
      print "(1X,A,' version ',A)", trim(progid),version

c     *****SET UP OUTPUT FILE OF POSIT(I,J) EVERY ISTPFR STEPS*****     

      if (batch) then
         read (*,'(A)') garbage ! read and ignore first line of cmd32 file - 
         print "(1X,'OUTPUT POSIT VS TIME FILE NAME IS : ',$)"
         read (*,'(A)') ofname
      else
        write (OFNAME, "('sh',I0,'.gout')") ICFLAG
      end if

      OPEN (UNIT=2,file=OFNAME,status='NEW',FORM='UNFORMATTED')

      ierr = 1
      trycnt = 0
      do while (ierr.eq.1)
         trycnt = trycnt + 1
         if (icflag == 1 .or. trycnt .gt. 1) call get_user_input
         FWDDEC = EXP(-DELTIM/FWDTAU)     
         BAKDEC = EXP(-DELTIM/BAKTAU)
         TRANS = SLIDE*DELTIM
         TIME = DELTIM     
         ISTEP = 1     
         M = 1     
         IFRAM = 1     
c     ***EVENT TIMES AND CHARGES*****
         i_savN = N
         CALL qcalc (IERR,icflag,sig)   !make 2-charge histories from spike data     
      end do                    ! loop if IERR.eq.1

c     +++++++++++++++++++++++++++++++++++++++++++++++++
c     -> -> ->  +++++AT THIS TIME WE  WRITE OUTPUT HEADER+++++
      WRITE (2) CBLOCK,IBLOCK,RBLOCK
c     ++++++++++++++++++++++++++++++++++++ 

c     now go to major routines for different options
      if (IOPTOP.eq.0) then
         call ndim_loop (0)     ! 'N' option std gravity
      else if (IOPTOP.eq.5) then
         call pair_loop(0)      ! 'P' option - 1-D, no tuning 
      else                      ! optimizer subroutines
         call pre_tuning
         if (IOPTOP.eq.6) then
            call pair_loop(1)   ! 'Q' option - 1-D, TUNING 
         else
            call ndim_loop(1)   ! FOR IOPTOP = 3  (3= aka 'O'; 6 = aka 'Q')
         end if
      end if

      print "(1x,'program: ',A)", PROGID
      PRINT "(1X, 'OUTPUT POSITION FILE: ',A)", OFNAME
      PRINT "(1x, 'INPUT DATA FILE: ',A)", IFNAME
      PRINT "(1X,'MEANINT  +/-',I3,' SPIKES;',   
     1'  ENDTIME: ',F10.1)", IHS,ENDTIME
      if (ISHIFTA.EQ.0) then       
         PRINT "(1X,'ACCEPTOR DECAY FORWARD')"      
      else       
         PRINT "(1X,'ACCEPTOR DECAY BACKWARD')"
      end if
      if( ISHIFTE.EQ.0) then
         PRINT "(1X,'EFFECTOR DECAY FORWARD')"
      else
         PRINT "(1x,'EFFECTOR DECAY BACKWARD')"     
      end if
      if (FS .eq. 1.0) then       
         PRINT "(1X,'Force sign =',F4.1,': excitation')", FS      
      else       
         PRINT "(1X,'Force sign =',F4.1,'; inhibition')", FS     
      end if
      PRINT "(1X,'time step=',F6.1,'; SLIDE=',E8.2,'; fwdinc=',F6.1,
     1'; bakinc=',F6.1)", DELTIM,SLIDE,FWDINC,BAKINC
      PRINT "(1X,'FWDTAU=',F6.1,'; BAKTAU=',F6.1,'; CRDI=',F6.1)",
     +     FWDTAU,BAKTAU,CRDI
      PRINT "(1X,'FINAL STATE: TIME (MS)=',F8.1,'; TIME STEPS=',I8)",
     +     TIME,ISTEP     
      print "(1x,'time between stimuli= ',F10.1,' ms.')", stimper
      print "(1x,'positions after every ',I4,' step(s) saved')", ISTPFR
      print "(1x,'total of ',I4,' frames saved')", LFRAM

      if ((IOPTOP.eq.5).or.(IOPTOP.eq.6)) then
         close (unit=2)
         print "(1X,'IOPTOP=5 or 6')"
         return
      end if

      DO  I=1,N
         PRINT "(10(1X,F6.1))", (POSIT(I,J),J=1,N)
      end do

      PRINT "(2X,' ')"
      DO I=2,N     
         DO J=1,I-1     
            PDISQ = 0.0     
            DO K=1,N     
               TEMP1 = POSIT(J,K) - POSIT(I,K)     
               PFDIR(I,J,K) = TEMP1
               PDISQ = PDISQ + TEMP1*TEMP1
            end do
            TEMPD(J) = SQRT(PDISQ)
         end do
         PRINT "(10(1X,F6.1))",(TEMPD(J),J=1,I-1)
      end do
      CLOSE (UNIT=2)
      call possub(icflag, sig)  ! calc distances between all pairs
      call postodir3d(icflag, sig) !calc 3d direction
      end subroutine gbatch

c     ****************************************************************
      SUBROUTINE GET_USER_INPUT
c     ****************************************************************
c       *****INITALIZE N, TIME, DEFINE T STP, TAUs,SP STP,WELL*****

      N = 0
      do while (N.gt.NMAX.or.N.LT.2)
         print "(1x,'NUMBER OF PARTICLES (I2) 2-',I0,' : ',$)", nmax
         read (5,"(I2)") N     
      end do
      
      print "(1X,'TIME STEP (IN MILLISECONDS, F6.1) : ',$)"    
      read (5,"(F6.1)") DELTIM 
      print "(1X,'STEP FOR UNIT FORCE, (E8.2) : ',$)"     
      read (5,"(E8.2)") SLIDE
      print "(1x,'enter IHS value, default=10: ',$)" ! now get IHS here
      read (5,"(I5)") IHS
      if (IHS.eq.0) IHS = 10
      
      itemp1 = 0
      do while (itemp1.eq.0)
         print "(1x,'ACCEPTOR DECAY FWD = +1, BAK = -1 : ',$)"     
         read (5,"(I2)") itemp1     
      end do
      ISHIFTA = merge (0, N, itemp1.gt.0)

      itemp2 = 0
      do while (itemp2.eq.0)
         print "(1x,'EFFECTOR DECAY FWD = +1, BAK = -1 : ',$)"     
         read (5,"(I2)") itemp2
      end do
      ISHIFTE = merge (0, N, itemp2.gt.0)

      print "(1x,'FORCE SIGN:excitation= +1.0,inhibition= -1.0:',$)"     
      read (5,"(F6.1)") FS
      print "(1x,'rate NORMalization: ONLY NONE (3) ALLOWED : '$)"
      read (5,"(I2)") NORM      ! Read in only for backward compatibility
      NORM = 3                  ! This is the old value for no normalization
c     NORM IS IGNORED AFTER THIS- ONLY NO NORMALIZATION IS USED
      print "(1X,'FWD CHARGE DECAY TAU (IN MILLISECONDS, F6.1) : ',$)"        
      read (5,"(F6.1)") FWDTAU        
      print "(1X, 'BAK CHARGE DECAY TAU (IN MILLISECONDS, F6.1) : ',$)"        
      read (5,"(F6.1)") BAKTAU        
      print "(1x,'FWD CHARGE INCREMENT AT AP TIME (F6.1) : ',$)"        
      read (5,"(F6.1)") FWDINC        
      print "(1x,'BAK CHARGE INCREMENT AT AP TIME (F6.1) : ',$)"        
      read (5,"(F6.1)") BAKINC     
      print "(1X,'ZERO FORCE WELL DIAMETER, (F6.1) : ',$)"     
      read (5,"(F6.1)") CRDI
c     choice of normal or "TUNED - lag optimized" gravity 
      IOPTOP = 0                ! default - no tuning, N-dimension
      print "(1x,'LAG OPTIMIZED - T, sd thresh - S, NONE - N: ',$)"
      read (*,'(A)') instr
      if((instr.eq.'o').or.(instr.eq.'O')) IOPTOP=3 ! N-dimensions - TUNED
c     *********************************************************
c     OPTION FOR PAIR ONLY (1-DIMENSION) ALGORITHM- pos file out  - no gout file
c     Uses ONLY acceptor charges - AND these MUST be set to FORWARD

      if((instr.eq.'p').or.(instr.eq.'P')) IOPTOP=5 ! PAIR ALGORITHM-NO TUNING

c     NOTE: option 'Q' for PAIR ALGORITHM WITH TUNING 
      if((instr.eq.'q').or.(instr.eq.'Q')) IOPTOP=6 ! PAIR ALGORITHM-TUNING

      end subroutine get_user_input

c     ****************************************************************
      SUBROUTINE PRE_TUNING     ! and readin
c     ****************************************************************
      irange = 150              ! number of steps (in deltim units) that
                                ! shift goes in each direction - must be
                                ! even integer. Total shift search =
                                ! 2*irange+1
      ibuff=int((4*FWDTAU)+1)   ! number of extra steps to ensure that
                                ! charge array boundaries are not
                                ! exceeded NOTE WELL: ibuff assumes that
                                ! BAKTAU=FWDTAU in absolute value
      istart = irange+ibuff
      iend=int( (ENDTIME/DELTIM) - ((istart)+(2*irange+10))) ! effectively
                                                             ! changes
                                                             ! actual
                                                             ! ENDTIME
      ISTEP = istart
      ipofset = 0               ! CLEAR ARRAY THAT HOLDS "LAGS" or
                                ! "OFFSETS" in INTEGER MSEC values in
                                ! DELTIM "units" if DELTIME=1, offset of
                                ! 50 = 50 msec;if DELTIM=2 50 = 100
                                ! msec., etc.
c READ in lags for pairs in SAME order used for gravity input
c NOTE WELL  - LAG to J relative to "fixed" I !! 
      open (9, ACCESS='SEQUENTIAL',FORM='FORMATTED',
     +     STATUS='OLD', FILE='offsets.gnew')
      DO I=2,N
         DO J=1,I-1
            read (9,"(2(i5,1x),i5)") ITA, JTA, ipofset(I,J) ! Trash ITA,JTA
            if (ipofset(I,J).gt.irange) ipofset(I,J) = irange !Limits
                                                              !range to
                                                              !prevent
                                                              !error
            if (ipofset(I,J).lt.-irange) ipofset(I,J) = -irange !Limits
                                                                !range
                                                                !to
                                                                !prevent
                                                                !error
         end do
      end do
      close (9)                 ! done reading offsets.gnew files from disk
      if (.not.sig) then
c     FOR DEBUG write out file just read in as temp.gnew
         open (8, ACCESS='SEQUENTIAL',FORM='FORMATTED',STATUS='UNKNOWN',
     +        FILE='temp.gnew')
         DO I=2,N
            DO J=1,I-1
c     next line writes shift for i,j to disk file for user info
               write (8,"(2(i5,1x),i5)") I, J, ipofset(I,J) ! to temp.gnew file
            end do
         end do
         close (8)              ! done writing temp.gnew files to disk
      end if
      end subroutine pre_tuning

c     ****************************************************************
      SUBROUTINE NDIM_LOOP (TUN)
c     ****************************************************************
      integer tun
      do while (.true.)
         DO I=2,N               ! particle
            DO J=1,I-1          ! particle
               PDISQ = 0.0     
               DO K=1,N         ! coordinate
                  TEMP1 = POSIT(J,K) - POSIT(I,K)     
                  PFDIR(I,J,K) = TEMP1     
                  PDISQ = PDISQ + TEMP1*TEMP1
               end do
               TEMP = 0.0     
               DIST = SQRT(PDISQ)     
               IF (DIST.GT.CRDI) TEMP = 1.0/DIST
               CHPR1 = TEMP*QAFS(ISTEP,I)*QE((ISTEP+tun*ipofset(I,J)),J)
               CHPR2 =-TEMP*QAFS((ISTEP+tun*ipofset(I,J)),J)*QE(ISTEP,I)
               DO K=1,N     
                  PFDIR(J,I,K) = PFDIR(I,J,K)*CHPR2     
                  PFDIR(I,J,K) = PFDIR(I,J,K)*CHPR1
               end do
            end do
         end do
c       *****GET TOTAL FORCE, NEW POSITION                   *****
c       *****NEW POS = OLD POS + (TOTF*TIME STEP)*SLIDE CONST*****
         DO I=1,N
            DO K=1,N
               TOTFK = 0.0
               DO J=1,N
                  TOTFK = TOTFK + PFDIR(I,J,K)
               end do
               POSIT(I,K) = POSIT(I,K) + TOTFK*TRANS
            end do
         end do
         ISTEP = ISTEP + 1     
         TIME = TIME + DELTIM     
         if ((tun.eq.1.and.ISTEP.EQ.iend)
     +        .or.(tun.eq.0.and.TIME.GT.ENDTIME)) return
c       *****CREATES MATRIX IN THE POSITION OUTPUT TABLE*****
c       ***includes FLASH vector***     
         M = M + 1
         IF (M.gt.ISTPFR.and.IFRAM.le.LFRAM) then
            WRITE (2) ((POSIT(I,J),J=1,N),I=1,N)
            WRITE (2) (FLASH(I),I=1,N) ! for compatibility only, not used
            IFRAM = IFRAM + 1
            M = 1
         end if
      end do
      end subroutine ndim_loop
      
c     ****************************************************************
      SUBROUTINE PAIR_LOOP(TUN)
c     ****************************************************************
      integer tun
      PRINT "(2X,'OUTPUT PAIR DIST VS TIME FILENAME IS: ',$)"
      if (sig) then
        if (ICFLAG.EQ.1) READ (*,'(A)') garbage
        write (posfile, "('sh',I0,'.pos')") ICFLAG
      else
         READ (*,'(A)') posfile
      end if

      close (unit=2)            !gout write never completed with option
                                !p - so close here
      OPEN(UNIT=2,FILE=posfile,STATUS='NEW',FORM='UNFORMATTED')
      write (2) CBLOCK,IBLOCK,RBLOCK !write header to "distance" file
      write (2) TIMAT           ! spike times for "raw data plot" over
                                ! time vs. distance plot
      K=1
      PDIS2 = 100.
      do while (.true.)
         DO I=2,N
            DO J=1,I-1
               TEMP = 0.0     
               DIST = real( PDIS2(I,J))
               IF (DIST.GT.CRDI) TEMP = 1.  
               PFDIR(I,J,K)= TEMP*QAFS(ISTEP,I)
     +              *QE(ISTEP+tun*ipofset(I,J),J)
            end do
         end do
         DO I=2,N
            DO J=1,I-1 
               TOTFK = PFDIR(I,J,K)
               PDIS2(I,J)= PDIS2(I,J) - TOTFK*TRANS
            end do
         end do
         ISTEP = ISTEP + 1     
         TIME = TIME + DELTIM     
         IF (TIME.GT.ENDTIME) then
            Close (unit=2)
            return              ! out to end of run
         end if
c     creates matrix in the position output table
         M = M + 1
         if (M.gt.ISTPFR.and.IFRAM.le.LFRAM) then
            PDIS2F = real(PDIS2)
            WRITE(2) ((PDIS2F(I,J),J=1,I-1),I=2,N)
            IFRAM = IFRAM + 1
            M = 1
         end if
      end do
      end subroutine pair_loop
      
      end program gbatch_08_v1
