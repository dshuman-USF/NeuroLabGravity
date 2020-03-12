        program blkshuf_05      
c       .. source file called fs.f for convenience
C       program converts *.gdt files to frameshifted *.rdt files
c       20-mar-90 increase nev from 4000 to 10000
c       7-jul-99 increase nev from 10000 to 30000 for special project
c       ..also allows numshf = 1 in addition to 21 and 101 - for special project
        parameter (ncd=74,nev=30000)
        INTEGER*4 DAT(ncd,nev),IDLIS(ncd),IBIGNUM
        INTEGER*4 iseed, iseed2
        REAL*4 RND
        dimension ipt2(ncd)
        CHARACTER*30 IFNAME,OFNAME
        dimension ipoint(ncd),ital(ncd)
        integer*4 times(ncd,nev),TIMES2(ncd,nev)
        integer*4 ttb(ncd) ! for each code. time between markb and first spike
        integer*4 tte(ncd) ! for each code. time between marke and last spike
c
c       ********************************************************************
c
c
c       v1.1 bgl march 1989
c       v2.0 bgl jan   1990
c               now enter only orig. file name and codes to be shifted
c               once. 20 different sh*rdt files generated..PROVIDED
c                number of time blocks in *.gdt file = 
c               at least 21 times the number of codes to be shifted.
c       v1 bgl convert a copy of frameshift2 to blksh1
c       This program rotates and shifts spike codes within one block defined by
c       a markb and a marke code
c       v 2.0 adds option of 100 shifts for .01 sig test.
c       v 3.0 increases ncd to 36 for extra codes other than 32 spike trains
c       7-jul-99 increase nev from 10000 to 30000 for special project
c
c       JAN 2000 - bgl:
c       v2000 adds one more loop so first shift of only 4 tau can be overwritten
c       to eliminate first sh1.rdt with only small offsets - for gsig 2000
c
c       v2kv9 - handles larger initial shifts [in first of 20 or 100 iterations only] 
c       (ishdel+ishtau*kcode) rather than just (ishtau*kcode)
c       because of tuning in gravity... readin variable decides whether it is used or not
c
c       may 25, 2000version blksh2kv11 - debug of development beta v9 and modified so that program 
c       always uses larger offsets in first shift
c       to eliminate effects of any real, BUT LONG LAG, correlation. The first shift is now
c       a multiple of ishdel - a type in - and subsequent shifts continue to be multiples of 4*tau,
c       with tau also a type in.
c
c
c
c    new MAJOR REVISION:
C                       PROGRAM: blkshuf_v3
c       Sept 3,2000  convert shift and rotate to random shuffle of each spike train.
C       INCORPORATES code fro blkshft2kv11 and tmplt8. 
c
c       This improvement will allow many neurons, many shuffles, and shorter samples.
c
c       *******************************************************************
C                       PROGRAM: blkshuf_05
c        3-jun-05 increases ncd to match lbs_04_v3.f code - the local shuffle option
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
        ibdt=0  ! subrtn rdfrm returns 1 if bdt input file type
        IBIGNUM= 2147483647 ! largest integer*4 integer 

C       CLEAR ARRAYs
C
        do 2 i=1,ncd
        ipoint(i)=0
        ipt2(i)=0
        ital(i)=0
        IDLIS(i)=0
        do 3 j=1,nev
        DAT(i,j)=IBIGNUM
        times(i,j)=IBIGNUM
        times2(i,j)=IBIGNUM
3       continue
2       continue
c
c       define seed for random # gen. transparently to user
        if(ICFLAG.eq.1)call ransub(iseed)
        iseed2=iseed ! in case need original value
c
c
c       DEBUG LINES
c
c        PRINT 12345,iseed
c       12345   FORMAT(2X,//,I11,' random number seed')
c
c       END DEBUG
c
c       print *,'REMEMBER - *.gdt file must begin w/ markb 21 and.. '
c       print *,'.. end with marke 22 !!!!! '
c       print *,' '
c       print *,' '
c       print *,'Program is used in conjunction with gravity:'
c       print *,'shuffles selected trains for use in control calculations.'
c       print *,'  '
c       print *,'  '
c       print *,'  '
        print *,'enter input file name  '
        if (ICFLAG.EQ.1) read (*,'(A)') IFNAME
        if (ICFLAG.EQ.2) IFNAME = 'sh1.rdt'
        if (ICFLAG.EQ.3) IFNAME = 'sh2.rdt'
        if (ICFLAG.EQ.4) IFNAME = 'sh3.rdt'
        if (ICFLAG.EQ.5) IFNAME = 'sh4.rdt'
        if (ICFLAG.EQ.6) IFNAME = 'sh5.rdt'
        if (ICFLAG.EQ.7) IFNAME = 'sh6.rdt'
        if (ICFLAG.EQ.8) IFNAME = 'sh7.rdt'
        if (ICFLAG.EQ.9) IFNAME = 'sh8.rdt'
        if (ICFLAG.EQ.10) IFNAME = 'sh9.rdt'
        if (ICFLAG.EQ.11) IFNAME = 'sh10.rdt'
        if (ICFLAG.EQ.12) IFNAME = 'sh11.rdt'
        if (ICFLAG.EQ.13) IFNAME = 'sh12.rdt'
        if (ICFLAG.EQ.14) IFNAME = 'sh13.rdt'
        if (ICFLAG.EQ.15) IFNAME = 'sh14.rdt'
        if (ICFLAG.EQ.16) IFNAME = 'sh15.rdt'
        if (ICFLAG.EQ.17) IFNAME = 'sh16.rdt'
        if (ICFLAG.EQ.18) IFNAME = 'sh17.rdt'
        if (ICFLAG.EQ.19) IFNAME = 'sh18.rdt'
        if (ICFLAG.EQ.20) IFNAME = 'sh19.rdt'
        if (ICFLAG.EQ.21) IFNAME = 'sh20.rdt'
        if (ICFLAG.EQ.22) IFNAME = 'sh21.rdt'
        if (ICFLAG.EQ.23) IFNAME = 'sh22.rdt'
        if (ICFLAG.EQ.24) IFNAME = 'sh23.rdt'
        if (ICFLAG.EQ.25) IFNAME = 'sh24.rdt'
        if (ICFLAG.EQ.26) IFNAME = 'sh25.rdt'
        if (ICFLAG.EQ.27) IFNAME = 'sh26.rdt'
        if (ICFLAG.EQ.28) IFNAME = 'sh27.rdt'
        if (ICFLAG.EQ.29) IFNAME = 'sh28.rdt'
        if (ICFLAG.EQ.30) IFNAME = 'sh29.rdt'
        if (ICFLAG.EQ.31) IFNAME = 'sh30.rdt'
        if (ICFLAG.EQ.32) IFNAME = 'sh31.rdt'
        if (ICFLAG.EQ.33) IFNAME = 'sh32.rdt'
        if (ICFLAG.EQ.34) IFNAME = 'sh33.rdt'
        if (ICFLAG.EQ.35) IFNAME = 'sh34.rdt'
        if (ICFLAG.EQ.36) IFNAME = 'sh35.rdt'
        if (ICFLAG.EQ.37) IFNAME = 'sh36.rdt'
        if (ICFLAG.EQ.38) IFNAME = 'sh37.rdt'
        if (ICFLAG.EQ.39) IFNAME = 'sh38.rdt'
        if (ICFLAG.EQ.40) IFNAME = 'sh39.rdt'
        if (ICFLAG.EQ.41) IFNAME = 'sh40.rdt'
        if (ICFLAG.EQ.42) IFNAME = 'sh41.rdt'
        if (ICFLAG.EQ.43) IFNAME = 'sh42.rdt'
        if (ICFLAG.EQ.44) IFNAME = 'sh43.rdt'
        if (ICFLAG.EQ.45) IFNAME = 'sh44.rdt'
        if (ICFLAG.EQ.46) IFNAME = 'sh45.rdt'
        if (ICFLAG.EQ.47) IFNAME = 'sh46.rdt'
        if (ICFLAG.EQ.48) IFNAME = 'sh47.rdt'
        if (ICFLAG.EQ.49) IFNAME = 'sh48.rdt'
        if (ICFLAG.EQ.50) IFNAME = 'sh49.rdt'
        if (ICFLAG.EQ.51) IFNAME = 'sh50.rdt'
        if (ICFLAG.EQ.52) IFNAME = 'sh51.rdt'
        if (ICFLAG.EQ.53) IFNAME = 'sh52.rdt'
        if (ICFLAG.EQ.54) IFNAME = 'sh53.rdt'
        if (ICFLAG.EQ.55) IFNAME = 'sh54.rdt'
        if (ICFLAG.EQ.56) IFNAME = 'sh55.rdt'
        if (ICFLAG.EQ.57) IFNAME = 'sh56.rdt'
        if (ICFLAG.EQ.58) IFNAME = 'sh57.rdt'
        if (ICFLAG.EQ.59) IFNAME = 'sh58.rdt'
        if (ICFLAG.EQ.60) IFNAME = 'sh59.rdt'
        if (ICFLAG.EQ.61) IFNAME = 'sh60.rdt'
        if (ICFLAG.EQ.62) IFNAME = 'sh61.rdt'
        if (ICFLAG.EQ.63) IFNAME = 'sh62.rdt'
        if (ICFLAG.EQ.64) IFNAME = 'sh63.rdt'
        if (ICFLAG.EQ.65) IFNAME = 'sh64.rdt'
        if (ICFLAG.EQ.66) IFNAME = 'sh65.rdt'
        if (ICFLAG.EQ.67) IFNAME = 'sh66.rdt'
        if (ICFLAG.EQ.68) IFNAME = 'sh67.rdt'
        if (ICFLAG.EQ.69) IFNAME = 'sh68.rdt'
        if (ICFLAG.EQ.70) IFNAME = 'sh69.rdt'
        if (ICFLAG.EQ.71) IFNAME = 'sh70.rdt'
        if (ICFLAG.EQ.72) IFNAME = 'sh71.rdt'
        if (ICFLAG.EQ.73) IFNAME = 'sh72.rdt'
        if (ICFLAG.EQ.74) IFNAME = 'sh73.rdt'
        if (ICFLAG.EQ.75) IFNAME = 'sh74.rdt'
        if (ICFLAG.EQ.76) IFNAME = 'sh75.rdt'
        if (ICFLAG.EQ.77) IFNAME = 'sh76.rdt'
        if (ICFLAG.EQ.78) IFNAME = 'sh77.rdt'
        if (ICFLAG.EQ.79) IFNAME = 'sh78.rdt'
        if (ICFLAG.EQ.80) IFNAME = 'sh79.rdt'
        if (ICFLAG.EQ.81) IFNAME = 'sh80.rdt'
        if (ICFLAG.EQ.82) IFNAME = 'sh81.rdt'
        if (ICFLAG.EQ.83) IFNAME = 'sh82.rdt'
        if (ICFLAG.EQ.84) IFNAME = 'sh83.rdt'
        if (ICFLAG.EQ.85) IFNAME = 'sh84.rdt'
        if (ICFLAG.EQ.86) IFNAME = 'sh85.rdt'
        if (ICFLAG.EQ.87) IFNAME = 'sh86.rdt'
        if (ICFLAG.EQ.88) IFNAME = 'sh87.rdt'
        if (ICFLAG.EQ.89) IFNAME = 'sh88.rdt'
        if (ICFLAG.EQ.90) IFNAME = 'sh89.rdt'
        if (ICFLAG.EQ.91) IFNAME = 'sh90.rdt'
        if (ICFLAG.EQ.92) IFNAME = 'sh91.rdt'
        if (ICFLAG.EQ.93) IFNAME = 'sh92.rdt'
        if (ICFLAG.EQ.94) IFNAME = 'sh93.rdt'
        if (ICFLAG.EQ.95) IFNAME = 'sh94.rdt'
        if (ICFLAG.EQ.96) IFNAME = 'sh95.rdt'
        if (ICFLAG.EQ.97) IFNAME = 'sh96.rdt'
        if (ICFLAG.EQ.98) IFNAME = 'sh97.rdt'
        if (ICFLAG.EQ.99) IFNAME = 'sh98.rdt'
        if (ICFLAG.EQ.100) IFNAME = 'sh99.rdt'
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
        if (IDT.gt.32) goto 5050

        PRINT 67,(IDLIS(I),I=1,IDT)
        goto 75
5040    PRINT 67,(IDLIS(I),I=1,16)
        PRINT 67,(IDLIS(I),I=17,IDT)
        goto 75
5050    PRINT 67,(IDLIS(I),I=1,16)
        PRINT 67,(IDLIS(I),I=17,32)
        PRINT 67,(IDLIS(I),I=33,IDT)
67      FORMAT (2X,16(I3,1X))
c
c
c   ****************************************************
c       THIS SECTION DONE ONLY AT TIME OF FIRST FILE READ
c   ****************************************************
c
c
75      if (ICFLAG.eq.1) then

c
c       get numshf either 20 or 100 
        print *,'enter # suffles(20 or 100)' 
        read (*,80) numshf
80      format (I5)
        end if
c
c
c
c
c       NOTE: This initial version shuffles all codes except markb, marke
c
c       define column number or array pointer and tallies of  markb, marke codes
c       enter here after first loop
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
c
c       **************************
c       MOVE intervals in TIMES array to TIMES2 array 
c
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


c       MOVE intervals in TIMES2 array back to TIMES array IN RANDOM ORDER
c
c
        DO 3500 I=1,IDT ! outer loop - for each code
                IF ((IDLIS(I).ne.markb).and.(IDLIS(I).ne. marke)) then
c
                  do 3400 J=1, (ITAL(I)-1) ! inner loop: intervals for each code
2000                RND=RAN(iseed)
                    IR=int((RND*FLOAT(ITAL(I)-1))+1)
                    if (IR.gt.ITAL(I)-1) goto 2000 
3000                IF(TIMES2(I,IR).EQ.0.0)GOTO 3010
                    TIMES(I,J)=TIMES2(I,IR)
                    TIMES2(I,IR)=0
                    goto 3400
3010                IR=IR+1
                    IF(IR.GT.ITAL(I)-1)IR=1
                    GOTO 3000
3400              continue
                end if
3500    continue
c
c
c       loops done for all codes except markb,marke
c
C       CONVERT intervals in TIMES array back to clock times in TIMES2
c
c       REMEMBER all times and intervals in clock counts -cc- (INTEGER*4)
c
c
      DO 4500 I=1,IDT ! outer loop - for each code
           IF ((IDLIS(I).ne.markb).and.(IDLIS(I).ne. marke)) then
c           put in first interval relative to original mark b

           TIMES2(I,1)= ttb(I) ! first cc for this code
                  do 4400 J=2, ITAL(I) ! inner loop: intervals for each code
                TIMES2(I,J)=TIMES2(I,J-1)+TIMES(I,J-1) ! add next interval to previous cc
4400            continue
           end if
c       take care of markb and marke
           IF (IDLIS(I).eq.markb) TIMES2(I,1) = 1 ! to ensure markb written first
           IF (IDLIS(I).eq.marke) TIMES2(I,1) = imketm
4500    continue




c       **********************************************
c       DEBUG ROUTINE TO SEE IF TIMES2 ARRAY LOOKS OK after all this code
c       **********************************************
c
c
c
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
c
c
c
c       MAJOR LOOP POINT HERE
c
c
c       initialize pointers,ETC.
        do 2100 nw=1,IDT
        ipt2(nw)=1
2100    continue
        mkeflg=0
        IDONE1=0
        print *,'      '
        print*,'..writing randomly shuffled data to output file..'
c
c       find low event in this sweep across DAT array
2150    low=IBIGNUM
        do 2200 iw=1,IDT
        low=MIN0(low,TIMES2(iw,ipt2(iw)))
2200    continue
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
c          update pointer of column containing  event. 
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
c
        do 2350 kw=1,IDT
        if (low.eq.TIMES2(kw,ipt2(kw))) then
           if ((IDLIS(kw).ne.markb).and.(IDLIS(kw).ne.marke)) then
           write (2,2490) IDLIS(kw),low
c          update pointer of column containing  event. 
           ipt2(kw)=ipt2(kw)+1
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
