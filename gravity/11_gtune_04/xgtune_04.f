        PROGRAM xgtunev_04
C       

c       ****************************************************
c       version for starbase graphics in X-windows
C       bgl January 2000
c       program reads in shifts.gnew file from gbatch2k* (or gsig2k*)
c       and displays in a histogram the sum of charge products for all 101 shifts 
c       v3 read in and plot mod record structure from gbatch2kv7 (mean and sd)
c       **********************************************************
c
C       place at top before parameters
        DIMENSION zshft(101), ICN(6)
        character*60 FNAME
c
c       *****  get starbase aliases and constants  **************************
        include 'config.defs'
        include 'sbparam.defs'
c       *********************************************************************
        print 6
6       FORMAT(2X,'..creating X window: starbase graphics')

C       NEXT 3 lines create custom window
C
        fildes=gopen (1700,420,-700,-5,'xgtune')
        ICON=0 ! cross option
        IXX=0  ! 1=full plot, 2= full plot and control overlay
        print 11
11      FORMAT(//,2X,'PROGRAM DISPLAYS charge product sum histograms',/,
     2  2X,' created by gbatch2k*',//)
c       set up ICN array - legacy problem - use only two of 6 
        ICN (1)= 0 ! this will be "ref" code
        ICN (2)= 0
        ICN (3)= 0 ! this will be "tgt" code
        ICN (4)= 0
        ICN (5)= 0
        ICN (6)= 0

c
c
c       loop point
c
c
c
5       CONTINUE
        GOTO 17
15      print 16
16      FORMAT(2X,'OPEN FILE ERROR... CORRECT FILE NAME ?')
17      print 20
20      FORMAT(2X,'File name ..OR..<CR> TO EXIT: ',$)
        read (*,'(A)') FNAME
        IF(FNAME.EQ.' ')GOTO 50
c       above goto 50 is exit program..
c
c
c
        open (UNIT=2,ACCESS='DIRECT', RECL=416, FILE=FNAME,
     +  FORM='UNFORMATTED',STATUS='OLD', ERR=15)

c
c
c       loop point
c
c
c
103     itgt = 0
        print 104
104     FORMAT(2X,'RECORD #(I5),OR <CR> FOR NEW FILE or EXIT')
        print 105
105     FORMAT(2X,'OR  Ref, Tgt (i.e., gravity loop I,J): ')
        read (*,106) IREC, itgt
106     FORMAT(I5,I5)
        
        IF((IREC.gt.0).and.(itgt.eq.0)) ioptz = 1 ! use rec num to search file
        IF((IREC.gt.0).and.(itgt.gt.0)) ioptz = 2 ! use i,j to  search file
        IF((IREC.EQ.0).and.(itgt.eq.0))CLOSE (UNIT=2)
        IF((IREC.EQ.0).and.(itgt.eq.0)) goto 5
c
c

c       next ifs  calculate either two cells or record number.. depending on option at entry
c       if gravity run handles more that 10000 neurons the loop counter will have to be increased
c       ... I wait for the day.
c
c
        if (ioptz.eq.1) then
        ict=0 ! loop control counter init.
        do I=2,10000
        do J=1,I-1
        ict=ict+1
        if (ict.eq.IREC) then
        ICN(1)=I ! ref code from gravity
        ICN(3)=J ! target code from gravity
        goto 3000 ! mission accomplished - exit loop here
        end if
        end do
        end do
        end if
c
        if (ioptz.eq.2) then
        ict=0 ! loop control counter init.
        do I=2,10000
        do J=1,I-1
        ict=ict+1
        if ((IREC.eq.I).and.(itgt.eq.J)) then
        IREC= ict
        ICN(1)=I ! ref code from gravity
        ICN(3)=J ! target code from gravity
        goto 3000 ! mission accomplished - exit loop here
        end if
        end do
        end do
        end if

3000    READ(2,REC=IREC,ERR=222) zmean,twsd,DELTIM,(zshft(I),I=1,101)
        BINW=DELTIM

c       plot zshft array as histogram
        call distune (IREC,zshft,BINW,FNAME,fildes,ICN,zmean,twsd)
        goto 103
c
222     CLOSE(UNIT=2)
        GOTO 5
50      STOP    !**********************PROGRAM EXIT*******************
        END
