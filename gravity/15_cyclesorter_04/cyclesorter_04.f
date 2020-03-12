        Program cyclesorter_04
c
c
c
c       ***************************************************************
c       ***************************************************************
c       v 1.0 bgl 18-Dec-98
c
c       PROGRAM INPUT:  *.BDT files with I, E pulses indicating 
C       respiratory phase transitions
C       and     analog channels, including "integrated" phrenic analog channel. 
c
c       This program allows each respiratory cycle or phase, or other defined "time chunk", to be defined
c       as an event. All spike codes and analog signals in each event are treated as a single entity
C       for the purpose of sorting the events as a function of a selected parameter
C       (peak phrenic amplitude, cycle duration, etc.)
c       The sorted data time chunks are written out in a sorted order in BDT format
c       with a uniform time interval between each chunk.
c
c       For files to be converted to a *.gdt gravity input file
c       be sure that the interval between chunks is 4 * the gravity tau (time constant)
c       Thus, the relative times of events within a chunk do not change, 
c       but their relation to times in other chunks may change.
c       Defined codes define the start and end of each chunk.

c       ***************************************************************
c       ***************************************************************
c       v _04 bgl 24-Nov-2004

c       New option 1 creates chunks written out IN ORIGINAL ORDER.
c       Allows concatendated chunks in original order for gravity *.gdt input data files 
c       NOTE: Chunks that exceed user defined length are excluded!!
c
c
        PARAMETER (lines=5000000) !max bdt lines to store
        PARAMETER (icyc=10000) ! max number of cycles in calc arrays 
        DIMENSION icode(lines) 
        integer*4 iclock(lines),jtrash
        DIMENSION istart(icyc), ispt (icyc), iept (icyc),
     1  iend(icyc), idur(icyc), ilo(icyc),
     2  ihi(icyc), irange(icyc), imean(icyc), isort(icyc),
     3  index(icyc),iramp(icyc)

        CHARACTER*60 ifname,ofname
c
c       reentry for another "run"
c
1111    do ic=1,icyc
        istart(ic)=0
        ispt(ic)=0
        iept(ic)=0
        iend(ic)=0
        idur(ic)=0
        ilo(ic)=0
        ihi(ic)=0
        iramp(ic)=0
        irange(ic)=0
        imean(ic)=0
        isort(ic)=0
        index(ic)=0
        end do
        do li=1,lines
        icode(li)=0
        iclock(li)=0
        end do

        

        PRINT 100
100     FORMAT(2X,'TYPE INPUT *.bdt FILENAME ..(CR:EXIT):')
        READ (*,'(A)') ifname
        IF(ifname .EQ. ' ') GOTO 300


        OPEN(UNIT=1,file=ifname,FORM='FORMATTED',STATUS='OLD')
C

1010    FORMAT (I5,I8)
C       
c       init pointer
        li=1
c       init flag to indicate reached end of input file
        iendflg=0
25      READ (1,1010,END=2000) itrash,jtrash    !read .BDT format
        if (itrash.eq.0) goto 25                ! handle blank records
        if (jtrash.eq.1111111) then     

        read(1,1010,end=2000) itrash, jtrash ! read final header line
        end if
c
c       done with bdt header handling
c
        do li=1,lines
        read(1,1010,end=2000) icode(li), iclock(li)
        lastrd=li !last value of li actually read
        end do
c       READ IN DONE
        goto 2010 ! did not reach end of file
c       calc duration of bdt read in
c       report that and if whole bdt file accepted.
2000    iendflg=1 ! reached end of input file
c       calculate total time of bdt read in minutes
2010    CLOSE(UNIT=1)
        tottime=(((iclock(lastrd)-iclock(1))*0.5)/1000.)/60.! in minutes
c
c       Report to user duration of bdt file and whether whole file was read in here
c       
        if (iendflg.eq.0) print 713
        if (iendflg.eq.1) print 714
 713    FORMAT(2X
     +       ,'Only PARTIAL readin of BDT file...increase array size')
714     FORMAT(2X,'Complete readin of BDT file..')
        print 716, tottime
716     format(2x,'Duration of data read in is: ',F6.2, '  min.')
c
c
c

        print 125
125     FORMAT(2X,'TYPE OUTPUT *.bdt FILENAME:')
        READ (*,'(A)') ofname
        OPEN(UNIT=2,file=ofname,FORM='FORMATTED',STATUS='NEW')
c
c       enter codes that define start and end of each 'cycle'
c       ...and select analog channel
c

40      PRINT 50
50      FORMAT(2X,'TYPE: (2I5) ISCODE,IECODE')
        PRINT 60
60      FORMAT(2X,'They can be same code')

        READ (*,13) ISCODE,IECODE
13      FORMAT(2I5)
        PRINT 120
120   FORMAT(2X,'ANALOG CHANNEL for PAREMETERS (I2):',$)
        READ (*,'(I2)') ICHN1
        ICHNA=ICHN1*4096
c
c       select criterion to use for later sort and output
c
9500    call pick_04(iprim)
          IF(iprim.le.6)goto 350
            GOTO 9500
c
c       define maximum cycle length allowed - gt will be filtered out at write time-TBA
c
350     Print 80
80    FORMAT(2x,'Enter max. phase duration..')
        PRINT 90
90    FORMAT(2x,'(sec, fp). Default max. = 10.0 sec: ',$)
        READ (*,6124) vmax
c       value now in clock ticks
        maxcyc = int(vmax * 2000)
        if (vmax.eq.0.0) maxcyc = 20000 ! asumes 0.5 msc clock ticks
6124    format (f10.2)
6125    format (f10.2)
c
c       enter lenght of time to stick between each chunk at write
        Print 85
85    FORMAT(2x,'Enter INTER-CHUNK interval..')
        PRINT 95
95    FORMAT(2x,'(msec, fp). Default max. = 40.0 msec: ',$)
        READ (*,6125) atau
        if(atau.eq.0.0) atau=40. 

        itau=int(atau*2.0)! integer number of 0.5 msec ticks
c
C       CONFIRM THAT SELECTIONS ARE VALID FOR CURRENT DATA FILE:
C
        iok=0
        eok=0
        do ij=1,lastrd
        IF (ISCODE.EQ.icode(ij)) iok=1
        IF (IECODE.EQ.icode(ij)) eok=1
        IF ((iok.eq.1).and.(eok.eq.1)) then
        goto 3030
        end if
        end do
        PRINT 70
        goto 40 ! try again
70      FORMAT(2x, 'CODES NOT FOUND')
3030    CONTINUE ! codes found
c
c       fill data arrays that will be sorted next
c       pass 1 - fill start and end times of each time chunk and pointers to bdt
c
c       
c       OPTION 1: iscode equal to iecode
c
        itot=0 ! set counter for total number of cycles stored in arrays
        if (iscode.eq.iecode) then
        ipt=1 ! init main array pointer
        ij = 1
3500    if (icode(ij).eq.iscode) then
        istart(ipt)=iclock(ij) ! store first clock count
        ispt(ipt)=ij               ! store first pointer to bdt array location
        goto 4000
        end if
        ij = ij+1
        if (ij.gt.lastrd)goto 6000 ! done
        goto 3500 ! in this loop until first code found
4000    ij = ij+1
        ipt = ipt+1
        if ((ij.gt.lastrd).or.(ipt.gt.icyc)) goto 6000 ! done
4050    if (icode(ij).eq.iscode) then
        istart(ipt)=iclock(ij) ! store next start clock count
        ispt(ipt)=ij               ! store next ispt pointer to bdt array location
        iend(ipt-1)=iclock(ij-1) ! store clock count of end of last chunk
        iept(ipt-1)=ij-1                   ! store pointer end of last chunk
        itot=itot+1 ! number of cycles done
        goto 4000
        end if
        ij = ij+1
        if (ij.gt.lastrd) goto 6000! done
        goto 4050
        end if

c
c       OPTION 2:  iscode not equal to iecode
c
        ipt=1                   ! init main array pointer
        ij = 1
3600    if (icode(ij).eq.iscode) then
        istart(ipt)=iclock(ij) ! store a start code's clock count
        ispt(ipt)=ij               ! store corresponding pointer to bdt array location
        ij = ij+1
        if (ij.gt.lastrd) goto 6000 ! done
        goto 4100
        end if
        ij = ij+1
        if (ij.gt.lastrd)goto 6000 ! done
        goto 3600 ! first loop
c
4100    if (icode(ij).eq.iecode) then
        iend(ipt)=iclock(ij) ! store clock count of end of last chunk
        iept(ipt)=ij               ! store pointer end of last chunk
        itot = itot+1 ! total cycles captured
        ij = ij+1
        ipt = ipt+1
        if ((ij.gt.lastrd).or.(ipt.gt.icyc)) goto 6000! done
        goto 3600
        end if
        ij = ij+1
        if (ij.gt.lastrd)goto 6000 ! done
        goto 4100
6000    continue ! done with first phase
c
c       use itot as number of complete cycles captured in arrays
c

c
c       Parameter: cycle duration-fill idur array
c
         
        do ij=1,itot
        idur(ij)= iend(ij)-istart(ij)
        end do
         
c
c       Parameter: analog channel parameters: fill high, low, range, mean        
c       ... and time to peak (iramp)
c

        
        do ik=1,itot
c
c       (re)initialize some variables
c
        max = -128000
        min = 128000
        isum=0 ! for calculating mean
        itally=0 !for calculating mean
c
c
c
        do ij=ispt(ik),iept(ik)
        IF (ICHN1.NE.(icode(ij)/4096)) GOTO 9550
        itally=itally+1
        icurrent=icode(ij)-ICHNA
c
c       next line handles 2's compliment AD conversion
c
        IF (icurrent.GT.2047) icurrent=icurrent-4096
c
c       next line makes all values positive
c
        icurrent=icurrent+2048
        isum=isum+icurrent
        max = MAX0(icurrent,max)
        if (max.eq.icurrent)imxtim=iclock(ij)
        min = MIN0(icurrent,min)
 9550   continue
        end do
c
c       fill arrays for THIS CYCLE
c
        ihi(ik)=max
        ilo(ik)=min
        irange (ik) = max-min
c
c       debug
c
        if(itally.eq.0) then
        print 6996
6996    format (2x, 'ERROR at 6996')
        stop
        end if
c
c
        imean(ik)=isum/itally
c
c       time to peak in cycle
c
        iramp(ik)=imxtim-istart(ik)
c       
c       now loop to do next cycle
c
        end do


c
c       Select array to use to sort data and fill temp sort array
c       OPTION 1 is special case - index, NOT isort, is used

        IF(iprim.EQ.1)THEN ! 
                do ij=1,itot
                index(ij)=ij
                end do
            ELSE IF(iprim.EQ.2)THEN
                do ij=1,itot
                isort(ij)= idur(ij)
                end do
          ELSE IF(iprim.EQ.3)THEN
                do ij=1,itot
                isort(ij)=ihi(ij)
                end do
          ELSE IF(iprim.EQ.4)THEN
                do ij=1,itot
                isort(ij)=ilo(ij)
                end do
          ELSE IF(iprim.EQ.5)THEN
                do ij=1,itot
                isort(ij)=irange(ij)
                end do
          ELSE IF(iprim.EQ.6)THEN
                do ij=1,itot
                isort(ij)=iramp(ij)
                end do
          ELSE ! expects 7 here
                do ij=1,itot
                isort(ij)=imean(ij)
                end do
          ENDIF
c
c       EXCEPT for ORIGINAL ORDER OPTION 1...
c       call sort routine to process isort array
c       pass both isort and index arrays. The index array will return
c       sorted pointers to the array represented by temporary array isort
c       
        if (iprim.ne.1) call sort(itot,isort,index)

c
c       using index array pointers,
c       write out 'sorted' bdt file, with tau between chunks. 
c       Chunks that exceed user defined length are excluded!!
c




c       using index array pointers,
c       write out 'sorted' bdt file, with tau between chunks. 
c       add itau "buffer" interval between each chunk to eliminate
c       boundary effects .. these will vary with analysis used.


c
c
c       outer loop cycles through each chunk
        print*,'..writing sorted data ..'
c
c       write bdt header once


           write (2,1010) 11,1111111
           write (2,1010) 11,1111111
c       
c       set up time correction
c
        lasttime=0
        do iw =1,itot
c
c       next line filters out long cycles
c
        if((iend(index(iw)))-(istart(index(iw))).ge.maxcyc) goto 4995
c
c       with each code and clock time in a particular chunk
c       
        icorfl=0 !flag to handle first time in each new chunk
        do iz=ispt(index(iw)),iept(index(iw))
        if (icorfl.eq.0)then
        lasttime = lasttime+itau
        write (2,1010) icode(iz),lasttime
        icorfl=1
        goto 4996 ! done this iteration of inner loop
        end if
        lasttime= lasttime+(iclock(iz)-iclock(iz-1))
        write (2,1010) icode(iz),lasttime
 4996   continue
        end do
 4995   continue
        end do



        close (unit=2)

c
c       loop to start of program
c
        GOTO 1111
300     STOP
        END





