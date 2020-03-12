        program lbs_04_v3       
C       
c       This program is used to generate surrogate data sets in "gdt" format
c       for Monte Carlo gravity analysis significance testing.
c
c
c
C       This program shuffles local ranges of interspike intervals
c       for all codes except markb marke and, optionally, one other code, e.g., a stimulus marker
c       or an "I code" to mark the onset of a respiratrory cycle.
c
c       Allowing one code series to remain undisturbed permits, for example,
c       average histograms of time series triggered by events of that code series - so
c       that the original averages can be compares with the locally shuffled averages.
c
c       The "original" gdt file is read in once and all the shuffled output files are derived from it-
c       i.e., each *.rdt file is independently created from the same source file.
c
c       The program steps along each time series in "nshuf" steps,
c       where nshuf is the number os sucessive interspike intervals to be locally 
c       shuffled and written out.
c       
c       The objective is to preserve the approximate times of firing rate fluctuations,
c       but remove "excess" precise short time scale synchrony.
c
c       ***************************NOTES**************************************
c
c       lbs_04_v3 bgl June 2004
c
c       Thus in this version, nshuf is a user entered variable between 3 and 10, inclusive;
c       default is 5.
c
c
c       The do loop for each code advances until the last 'nshuf' +1 events in each channel
c       as a simple way to deal 
c       with end conditions in this first version; these events are not shuffled.
c
c       Channels with fewer than 2*(nshuf+1) events are not changed
c 
c       Future could include dithering of the middle
c       spike in "two interspike interval blocks" prior to the local shuffling
c       
C      COMPILE with -K option (static storage) - not stack 
c
c
c       ************************************************************
        parameter (ncd=74,nev=30000)
c       ncd = max number of different codes allowed
c       nev = max number of each code allowed
c       ************************************************************
c
        INTEGER*4 DAT(ncd,nev),IDLIS(ncd),IBIGNUM,iplow(ncd)
        INTEGER*4 iseed, iseed2,imkbtm,imketm
        REAL*4 RND
        dimension ipt2(ncd)
        CHARACTER*128 IFNAME,OFNAME
        dimension ipoint(ncd),ital(ncd)
        integer*4 times(ncd,nev),TIMES2(ncd,nev),ICURVAL
        integer*4 ttb(ncd) ! for each code. time between markb and first spike
        integer*4 tte(ncd) ! for each code. time between marke and last spike
c
c
c
c
        ICFLAG = 1
c
c       Reentry point...
c
c
1       markb=21
        marke=22
        ibdt=0  ! subrtn rdfrm returns 1 if bdt file
        IBIGNUM= 2147483647 ! largest integer*4 integer 

C       CLEAR ARRAYs
C
        do 2 i=1,ncd
        ipoint(i)=0
        ipt2(i)=0
        ital(i)=0
        IDLIS(i)=0
        iplow(i)=0
        do 3 j=1,nev
        DAT(i,j)=IBIGNUM
        times(i,j)=IBIGNUM
        times2(i,j)=IBIGNUM
3       continue
2       continue
c
c       define seed for random # gen. transparently to user
c        if(ICFLAG.eq.1)call ransub(iseed)
c        iseed2=iseed ! in case need original value
c
c
c       DEBUG LINES
c
c       PRINT 12345,iseed
c       12345   FORMAT(2X,//,I11,' random number seed')
c
c       END DEBUG
c
c       print *,'REMEMBER - *.gdt file must begin w/ markb 21 and.. '
c       print *,'.. end with marke 22 !!!!! '
c       print *,' '
c       print *,' '
c       print *,'Program is used in conjunction with gravity:'
c       print *,'Does LOCAL random shuffles within each successive block of
c       nshuf intervals
c       print *,'...for each code; EXCLUDES last nshuf intervals
c       
c       print *,'at end of each code series.''
c       print *,'  '
c       print *,'  '
        print *,'enter input file name  '
        if (ICFLAG.EQ.1) read (*,'(A)') IFNAME
c
c
        if (ICFLAG.EQ.1) OFNAME = 'sh1.rdt'
        if (ICFLAG.EQ.2) OFNAME = 'sh2.rdt'
        if (ICFLAG.EQ.3) OFNAME = 'sh3.rdt'
        if (ICFLAG.EQ.4) OFNAME = 'sh4.rdt'
        if (ICFLAG.EQ.5) OFNAME = 'sh5.rdt'
        if (ICFLAG.EQ.6) OFNAME = 'sh6.rdt'
        if (ICFLAG.EQ.7) OFNAME = 'sh7.rdt'
        if (ICFLAG.EQ.8) OFNAME = 'sh8.rdt'
        if (ICFLAG.EQ.9) OFNAME = 'sh9.rdt'
        if (ICFLAG.EQ.10) OFNAME = 'sh10.rdt'
        if (ICFLAG.EQ.11) OFNAME = 'sh11.rdt'
        if (ICFLAG.EQ.12) OFNAME = 'sh12.rdt'
        if (ICFLAG.EQ.13) OFNAME = 'sh13.rdt'
        if (ICFLAG.EQ.14) OFNAME = 'sh14.rdt'
        if (ICFLAG.EQ.15) OFNAME = 'sh15.rdt'
        if (ICFLAG.EQ.16) OFNAME = 'sh16.rdt'
        if (ICFLAG.EQ.17) OFNAME = 'sh17.rdt'
        if (ICFLAG.EQ.18) OFNAME = 'sh18.rdt'
        if (ICFLAG.EQ.19) OFNAME = 'sh19.rdt'
        if (ICFLAG.EQ.20) OFNAME = 'sh20.rdt'
        if (ICFLAG.EQ.21) OFNAME = 'sh21.rdt'
        if (ICFLAG.EQ.22) OFNAME = 'sh22.rdt'
        if (ICFLAG.EQ.23) OFNAME = 'sh23.rdt'
        if (ICFLAG.EQ.24) OFNAME = 'sh24.rdt'
        if (ICFLAG.EQ.25) OFNAME = 'sh25.rdt'
        if (ICFLAG.EQ.26) OFNAME = 'sh26.rdt'
        if (ICFLAG.EQ.27) OFNAME = 'sh27.rdt'
        if (ICFLAG.EQ.28) OFNAME = 'sh28.rdt'
        if (ICFLAG.EQ.29) OFNAME = 'sh29.rdt'
        if (ICFLAG.EQ.30) OFNAME = 'sh30.rdt'
        if (ICFLAG.EQ.31) OFNAME = 'sh31.rdt'
        if (ICFLAG.EQ.32) OFNAME = 'sh32.rdt'
        if (ICFLAG.EQ.33) OFNAME = 'sh33.rdt'
        if (ICFLAG.EQ.34) OFNAME = 'sh34.rdt'
        if (ICFLAG.EQ.35) OFNAME = 'sh35.rdt'
        if (ICFLAG.EQ.36) OFNAME = 'sh36.rdt'
        if (ICFLAG.EQ.37) OFNAME = 'sh37.rdt'
        if (ICFLAG.EQ.38) OFNAME = 'sh38.rdt'
        if (ICFLAG.EQ.39) OFNAME = 'sh39.rdt'
        if (ICFLAG.EQ.40) OFNAME = 'sh40.rdt'
        if (ICFLAG.EQ.41) OFNAME = 'sh41.rdt'
        if (ICFLAG.EQ.42) OFNAME = 'sh42.rdt'
        if (ICFLAG.EQ.43) OFNAME = 'sh43.rdt'
        if (ICFLAG.EQ.44) OFNAME = 'sh44.rdt'
        if (ICFLAG.EQ.45) OFNAME = 'sh45.rdt'
        if (ICFLAG.EQ.46) OFNAME = 'sh46.rdt'
        if (ICFLAG.EQ.47) OFNAME = 'sh47.rdt'
        if (ICFLAG.EQ.48) OFNAME = 'sh48.rdt'
        if (ICFLAG.EQ.49) OFNAME = 'sh49.rdt'
        if (ICFLAG.EQ.50) OFNAME = 'sh50.rdt'
        if (ICFLAG.EQ.51) OFNAME = 'sh51.rdt'
        if (ICFLAG.EQ.52) OFNAME = 'sh52.rdt'
        if (ICFLAG.EQ.53) OFNAME = 'sh53.rdt'
        if (ICFLAG.EQ.54) OFNAME = 'sh54.rdt'
        if (ICFLAG.EQ.55) OFNAME = 'sh55.rdt'
        if (ICFLAG.EQ.56) OFNAME = 'sh56.rdt'
        if (ICFLAG.EQ.57) OFNAME = 'sh57.rdt'
        if (ICFLAG.EQ.58) OFNAME = 'sh58.rdt'
        if (ICFLAG.EQ.59) OFNAME = 'sh59.rdt'
        if (ICFLAG.EQ.60) OFNAME = 'sh60.rdt'
        if (ICFLAG.EQ.61) OFNAME = 'sh61.rdt'
        if (ICFLAG.EQ.62) OFNAME = 'sh62.rdt'
        if (ICFLAG.EQ.63) OFNAME = 'sh63.rdt'
        if (ICFLAG.EQ.64) OFNAME = 'sh64.rdt'
        if (ICFLAG.EQ.65) OFNAME = 'sh65.rdt'
        if (ICFLAG.EQ.66) OFNAME = 'sh66.rdt'
        if (ICFLAG.EQ.67) OFNAME = 'sh67.rdt'
        if (ICFLAG.EQ.68) OFNAME = 'sh68.rdt'
        if (ICFLAG.EQ.69) OFNAME = 'sh69.rdt'
        if (ICFLAG.EQ.70) OFNAME = 'sh70.rdt'
        if (ICFLAG.EQ.71) OFNAME = 'sh71.rdt'
        if (ICFLAG.EQ.72) OFNAME = 'sh72.rdt'
        if (ICFLAG.EQ.73) OFNAME = 'sh73.rdt'
        if (ICFLAG.EQ.74) OFNAME = 'sh74.rdt'
        if (ICFLAG.EQ.75) OFNAME = 'sh75.rdt'
        if (ICFLAG.EQ.76) OFNAME = 'sh76.rdt'
        if (ICFLAG.EQ.77) OFNAME = 'sh77.rdt'
        if (ICFLAG.EQ.78) OFNAME = 'sh78.rdt'
        if (ICFLAG.EQ.79) OFNAME = 'sh79.rdt'
        if (ICFLAG.EQ.80) OFNAME = 'sh80.rdt'
        if (ICFLAG.EQ.81) OFNAME = 'sh81.rdt'
        if (ICFLAG.EQ.82) OFNAME = 'sh82.rdt'
        if (ICFLAG.EQ.83) OFNAME = 'sh83.rdt'
        if (ICFLAG.EQ.84) OFNAME = 'sh84.rdt'
        if (ICFLAG.EQ.85) OFNAME = 'sh85.rdt'
        if (ICFLAG.EQ.86) OFNAME = 'sh86.rdt'
        if (ICFLAG.EQ.87) OFNAME = 'sh87.rdt'
        if (ICFLAG.EQ.88) OFNAME = 'sh88.rdt'
        if (ICFLAG.EQ.89) OFNAME = 'sh89.rdt'
        if (ICFLAG.EQ.90) OFNAME = 'sh90.rdt'
        if (ICFLAG.EQ.91) OFNAME = 'sh91.rdt'
        if (ICFLAG.EQ.92) OFNAME = 'sh92.rdt'
        if (ICFLAG.EQ.93) OFNAME = 'sh93.rdt'
        if (ICFLAG.EQ.94) OFNAME = 'sh94.rdt'
        if (ICFLAG.EQ.95) OFNAME = 'sh95.rdt'
        if (ICFLAG.EQ.96) OFNAME = 'sh96.rdt'
        if (ICFLAG.EQ.97) OFNAME = 'sh97.rdt'
        if (ICFLAG.EQ.98) OFNAME = 'sh98.rdt'
        if (ICFLAG.EQ.99) OFNAME = 'sh99.rdt'
        if (ICFLAG.EQ.100) OFNAME = 'sh100.rdt'
c
c
        OPEN (1,FILE=IFNAME,STATUS='OLD',FORM='FORMATTED')
        OPEN (2,FILE=OFNAME,STATUS='NEW',FORM='FORMATTED')
c
c       read in data file to be shifted
        call rdfrm (DAT,ipoint,ital,ibdt)
        if (ibdt.eq.1) then
        print *,'INPUT FILE *.bdt type'
        end if
c
c       DEFINE NUM OF CODES, CODES, AND NUM EVTS/CODE:
c
        IDT=0
        JDT=1
        DO  I=1,ncd
        IF (ipoint(I).NE.0) THEN
        IDT=IDT+1
        IDLIS(JDT)=ipoint(I)
        JDT=JDT+1
        END IF
        end do
        PRINT 11,IDT
11      FORMAT(2X,//,I3,' ID CODES WERE INPUT')
c
c
        if ((IDT.gt.16).and.(IDT.le.32)) goto 5040
        if ((IDT.gt.32).and.(IDT.le.48)) goto 5050
        if ((IDT.gt.48).and.(IDT.le.64)) goto 5060
        if (IDT.gt.64) goto 5070

        PRINT 67,(IDLIS(I),I=1,IDT)
        goto 75
5040    PRINT 67,(IDLIS(I),I=1,16)
        PRINT 67,(IDLIS(I),I=17,IDT)
        goto 75
5050    PRINT 67,(IDLIS(I),I=1,16)
        PRINT 67,(IDLIS(I),I=17,32)
        PRINT 67,(IDLIS(I),I=33,IDT)
        goto 75
5060    PRINT 67,(IDLIS(I),I=1,16)
        PRINT 67,(IDLIS(I),I=17,32)
        PRINT 67,(IDLIS(I),I=33,48)
        PRINT 67,(IDLIS(I),I=49,IDT)
        goto 75
5070    PRINT 67,(IDLIS(I),I=1,16)
        PRINT 67,(IDLIS(I),I=17,32)
        PRINT 67,(IDLIS(I),I=33,48)
        PRINT 67,(IDLIS(I),I=49,64)
        PRINT 67,(IDLIS(I),I=65,IDT)
67      FORMAT (2X,16(I3,1X))
c
c
c   ****************************************************
c       THIS SECTION DONE ONLY AT TIME OF INPUT FILE READ
c   ****************************************************
c
c
75      if (ICFLAG.eq.1) then

c
c       get numshf either 20 or 100 
        print *,'enter # suffled output files(20 or 100)' 
        read (*,80) numshf
80      format (I5)
87      format (I10)
c
c
c
c       get iseed 
        print *,'enter random number seed (odd integer.LE.10 digits)' 
        read (*,87) iseed
c
        iseed2=iseed ! in case need original value
        PRINT 12345,iseed
12345   FORMAT(2X,//,I10,' random number seed')
c
c
        print *,'enter I code - <CR> to skip' 
        read (*,80) icode
        
83      print *,'enter # intervals to shuffle locally (3-10;def=5)' 
        read (*,80) nshuf
        if(nshuf.eq.0) nshuf = 5
        if ((nshuf.lt.3).or.(nshuf.gt.10)) goto 83
        end if
c
c
c
c       *****************************
c
c       enter here after first loop
c
c       *****************************
c
c       IDT = number of different codes in file - inc markb and marke
c
        do 320 jk=1,IDT
        if (ipoint(jk).eq.markb) then
        imkbtm=DAT(jk,1) ! save to define time in cc of first spike after
c                       shuffle
        end if
        if (ipoint(jk).eq.marke) then
        imketm=DAT(jk,1) 
        end if

320     continue
c
c
c
c       *************************************************
c       loop to shuffle user selected codes
c       *************************************************
c       CONVERT SPIKE TIMES in DAT ARRAY TO INTERSPIKE INTERVALS
c       in TIMES ARRAY
c
c
        DO I=1,IDT ! outer loop - for each code
        IF ((IDLIS(I).ne.markb).and.(IDLIS(I).ne. marke)) then
            ttb(I)= DAT(I,1)-imkbtm ! store time between markb and fist spike
            if (ttb(I).eq.0) ttb(I)=1 ! can not equal 0  clock counts
            tte(I)= imketm-DAT(I,ITAL(I))! store time between last spike & marke
                do J=1,ITAL(I)-1 ! inner loop for each spike event
                TRV=(DAT(I,J+1)-DAT(I,J))
                if (TRV.eq.0.0) then
                print *, 'FATAL ERROR - suspect duplicate spike'
                stop
                end if
                TIMES(I,J)=int(TRV)
                end do
                end if
          end do
c
c       **************************
c       Correct imketm value for marke pulse relative to markb 
c       ... i.e., all codes are time shifted by this program so can not use 
c       "original" time values

        imketm=imketm-imkbtm
c       **************************
c
c
c
c
c       MOVE intervals in TIMES array to TIMES2 array 
c       **************************
c
c

        DO I=1,IDT ! outer loop - for each code
        IF ((IDLIS(I).ne.markb).and.(IDLIS(I).ne. marke)) then
                DO J=1,(ITAL(I)-1)
                TIMES2(I,J)=TIMES(I,J)
                TIMES(I,J)=0

                end do
                end if
        end do
c
c
c       **************************
c       MOVE intervals in TIMES2 array back to TIMES array
c       Move 'nshuf' intervals in random order to create 'locally'
c       shuffled data
c       **************************
c
c       skip markb, marke, icode channels and channels with small tallies...
        DO 3500 I=1,IDT ! outer loop - for each code
c
c
c
c       FIRST HANDLE channels that will NOT be shuffled here
C
C
c       skip markb, marke channels 
c
        IF ((IDLIS(I).eq.markb).or.(IDLIS(I).eq. marke)) goto 3500


c       Handle icode channel and channels with small tallies...
c
        if ((IDLIS(I).eq.icode).or.(ITAL(I).lt.2*(nshuf+1))) then

        do 17520 j=1,ITAL(I)-1
        TIMES(I,j)=TIMES2(I,j)
17520   continue
        goto 3500
        end if
c
c
c
c
c       HANDLE channels to be shuffled here...
c
C       nshuf is the number of successive intervals to be
C       randomly redistributed locally within the same region of the spike train
C       by the inner most loop 
C
c       NOTE-next line also handles end of array 
        do 3400 J=1, (ITAL(I)-(nshuf+1)),nshuf ! middle loop: intervals for one code
        icount = nshuf
c
c       ******************************
c       shuffle engine
c
c
        do 2100 K=J, J+(nshuf-1) ! inner loop: randomly place nshuf intervals
        ICURVAL= TIMES2(I,K)
2000    RND=RAN(iseed)
        IR=int((1+(RND*FLOAT(icount))))
        if (IR.gt.(icount)) goto 2000
        jcount = 1
c
c       find the IRth empty bin in the TIMES array and fill
        do 2200 L=J,J+(nshuf-1) 
        IF ((TIMES (I,L).eq.0.0).and.(IR.eq.jcount))then
        TIMES (I,L) = ICURVAL
        icount = icount-1
        goto 2100
        end if
        IF ((TIMES (I,L).eq.0.0).and.(IR.ne.jcount))then
        jcount=jcount+1
        end if
c
c       other conditions must only need to increment L so loop
c
c
2200    continue
2100    continue
c
3400    continue
c
c       end shuffle engine
c       *******************************
c
c
c
c       outer loop end point:
3500    continue

c
c       loops done for all codes except markb,marke
c
C       CONVERT intervals in TIMES array back to clock times in TIMES2
c
c       REMEMBER all times and intervals in clock counts -cc- (INTEGER*4)
c
c
c
c       FIRST HANDLE channels that will NOT be shuffled here
C
C
c
        DO 4500 I=1,IDT ! outer loop - for each code
c
c
       IF ((IDLIS(I).ne.markb).and.(IDLIS(I).ne. marke)) then
c
c           put in first interval relative to original mark b
        TIMES2(I,1)= ttb(I) ! first cc for this code
c
c      Handle icode channel and channels with small tallies...
c
        if ((IDLIS(I).eq.icode).or.(ITAL(I).lt.2*(nshuf+1))) then
        do 4450 J=2, ITAL(I) ! inner loop: intervals for each code
c       add next interval to previous cc
        TIMES2(I,J)=TIMES2(I,J-1)+TIMES(I,J-1) 
4450    continue
        goto 4500
        end if
c
c        HANDLE ALL OTHER CONDITIONS HERE
c
        do 4400 J=2, ITAL(I)-(nshuf+1) ! HANDLES END CONDITION LIKE ABOVE
c       add next interval to previous cc
        TIMES2(I,J)=TIMES2(I,J-1)+TIMES(I,J-1) 
4400    continue
        goto 4500
c
c
        END IF
c
c
c
c       take care of markb and marke
        IF (IDLIS(I).eq.markb) TIMES2(I,1) = 1 ! to ensure markb written first
        IF (IDLIS(I).eq.marke) TIMES2(I,1) = imketm
c
4500    continue
c       **********************************************
c       DEBUG ROUTINE TO SEE IF TIMES2 ARRAY LOOKS OK after all this code
c       **********************************************
c
c       DO I=1,IDT ! outer loop - for each code
c       DO J=1, ITAL(I)
c       print 6666, IDLIS(I), TIMES2(I,J)
c       end do
c       end do
c       6666    format(2x, 'code: ',I5, 'cc :', I8)
c
c
c       **********************************************
c       END OF DEBUG ROUTINE 
c       **********************************************
c
c
c       *************************************************
c       write out TIMES2 array in *.bdt format.
c       ...for given time: markb code FIRST out, marke LAST out.
c       ***************************************
c       basic algorithm
c
c       1. write out in order of temporal occurrence
c
c       2.  markb first written,
c       ..marke last written.
c
c       ********************************************
c
c       MAJOR LOOP POINT HERE
c
c       initialize pointers,ETC.
        do 2575 nw=1,IDT
        ipt2(nw)=1
2575    continue
        mkeflg=0
        IDONE1=0
        print *,'      '
        print*,'..writing randomly shuffled data to output file..'
c
c       find low event in this sweep across DAT array
c
c
2150    low=IBIGNUM
        do 2275 iw=1,IDT
        low=MIN0(low,TIMES2(iw,ipt2(iw)))
2275    continue
c       
c
c       write out markb and cc if present first
c       and check for marke event to be written last
c
c
        do 2300 jw=1,IDT
        if (low.eq.TIMES2(jw,ipt2(jw))) then
c       write markb code first if present at this clock count
        if (IDLIS(jw).eq.markb) then
        if (IDONE1.eq.1)goto 2320
c       update pointer of column containing  event. 
        ipt2(jw)=ipt2(jw)+1
c       WRITE BDT HEADER HERE JUST ONCE BEFORE MARKB
           write (2,2490) 11,1111111
           write (2,2490) 11,1111111
           write (2,2490) markb,low
c          IDONE1=1
           end if
c       set marke detection flag if = marke so can be..
c       written at end of this loop as last code with this cc
c       increment buffer to be written
2320       if (IDLIS(jw).eq.marke) then
c          update pointer of column containing  event. 
           ipt2(jw)=ipt2(jw)+1
           mkeflg=1
           end if
        end if
2300    continue
c
c
c       write out any NON markb, marke codes
c       note - compare with  iplow(kw) value is to prevent end of record errors
c       due to major loops above not going att the way to to ITAL-1 for some codes
c
c
        do 2350 kw=1,IDT
        if (low.eq.TIMES2(kw,ipt2(kw))) then
           if ((IDLIS(kw).ne.markb).and.(IDLIS(kw).ne.marke)) then
           if (low.le.iplow(kw)) goto 2323
           write (2,2490) IDLIS(kw),low
           iplow(kw)=low
c          update pointer of column containing  event. 
2323           ipt2(kw)=ipt2(kw)+1
           end if
        end if
2350    continue
c
c       write out marke code & check done flag
c
        if (mkeflg.eq.1) then
        write (2,2490) marke,low
        goto 200
        end if
        goto 2150
c
2490    format (I5,I8)
200     CLOSE (UNIT=1)
        close (unit=2)
        print 390,(ICFLAG)
390     FORMAT(2X,'CYCLE# ',I4)
        ICFLAG=ICFLAG+1
        if (numshf.eq.1) then
        if (ICFLAG.LE.1) goto 1
        end if
        if (numshf.eq.20) then
        if (ICFLAG.LE.20) goto 1
        end if
        if (numshf.eq.100) then
        if (ICFLAG.LE.100) goto 1
        end if
        END

