C
C       CLOCK COUNTS (MUL BY 0.5 TO GET MSEC., F.P.) 
C
C       READS IN UP TO N SPIKE CODES, NCD= MAX SYSTEM INPUTS POSSIBLE
C       version RDAT:      NCD= 14 AS OF JUNE 1987
C       version RDAT2:     NCD= 32 AS OF JUNE 1991
C       version rdfrm4:    handles adt or bdt format
C       NEV=10000;MAX NUM EVENTS PER CODE
C       *.ADT INPUT DATA HAS THE FORMAT I2,I8: EVENT CODE,CLOCK COUNTSC
C       *.BDT INPUT DATA HAS THE FORMAT I5,I8 ..FIRST 2 RECORDS
C               ARE: 11,1111111
C                    11,LSTTIM WHERE LSTTIM IS LAST CLOCK COUNT IN FILE
C                    ...THUS 1111111 IS USED TO DISTINGUISH .ADT FROM .BDT
C                    FILES.
c
C       ******* BE SURE ETEMPL,NEV,IPOINT ARE 32 BIT INTEGER VARIABLES *******
c
c
c       bgl 27-feb-89 version RDFRM: this source code modified as rdfrm for 
c       linking with FORTRAN
c       main program frameshift.f...
C
c       v.5 increased nev to 30000
        SUBROUTINE RDFRM(DAT,IPOINT,ITAL,ibdt)
        PARAMETER (NCD=74,NEV=30000)
        integer*4 ETEMPL, DAT(ncd,nev)
        DIMENSION IPOINT(NCD),ITAL(NCD)
1000    FORMAT (I2,I8)
1010    FORMAT (I5,I8)
C       
C
C       FILL DAT ARRAY, leave times as integer clock counts
c
        ibdt=0 ! default is adt format:w
25      READ (1,1010,END=2000) ITEMPL,ETEMPL    !read .BDT format
        if (ITEMPL.eq.0) goto 25                ! handle blank records
        if (ETEMPL.eq.1111111) then             !if *.BDT type file then..
        ibdt=1
        read(1,1010,end=2000) ITEMPL,ETEMPL     !read end time record
        goto 3030                               !process as .BDT file
        end if
        print 26
26      format(2x,'...reading .ADT type file..')
        rewind 1                                !.ADT file;go back to start 
30      READ (1,1000,END=2000) ITEMPL,ETEMPL    !read .ADT format
        if(ITEMPL.eq.0) goto 30 ! handle blank records
        do 100 I=1,NCD
        if (ITEMPL.eq.IPOINT(I)) goto 300
100     continue
        do 200 I=1,NCD
        if (IPOINT(I).eq.0) then
        IPOINT(I) = ITEMPL
        goto 300
        end if
200     continue
        print 260
260     format (2x,'FATAL ERROR IN rdfrm sub rtn... too many'
     1  ,/,2x,' different id codes [> NCD parameter]')
        stop
300     ITAL(I)=ITAL(I)+1
        if (ITAL(I).gt.NEV) then
        print *,'truncation error: a code has > NEV events'
        print *,'increase NEV in main and subrtns & recompile..'
        print *,'or decrease max. for each code in input file'
        stop
        end if
        DAT(I,ITAL(I))= ETEMPL
        goto 30
c
c
c
c       new block of code for  handling *.BDT files
c
c
c
3030    READ (1,1010,END=2000) ITEMPL,ETEMPL
        if((ITEMPL.eq.0).or.(ITEMPL.gt.999)) goto 3030 ! ignore ID code 
c
c       the above if statement filters out analog channel and voltage 
c       data points in the *.BDT file
c
        do 3100 I=1,NCD
        if (ITEMPL.eq.IPOINT(I)) goto 3300
3100    continue
        do 3200 I=1,NCD
        if (IPOINT(I).eq.0) then
        IPOINT(I) = ITEMPL
        goto 3300
        end if
3200    continue
        print 3260
3260    format (2x,'FATAL ERROR IN rdat3 sub rtn; excess'
     1  ,/,2x,' id codes [> NCD parameter]')
        stop
3300    ITAL(I)=ITAL(I)+1
        if (ITAL(I).gt.NEV) then
        print *,'truncation error: a code has > NEV events'
        print *,'increase NEV in main and subrtns & recompile..'
        print *,'or decrease max. for each code in input file'
        stop
        end if
        DAT(I,ITAL(I))= ETEMPL
        goto 3030
c
c
c
2000    RETURN
        END
