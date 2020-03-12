c     spkpat6kbg derived from spkpat6bg - bgl - Aug 98
c     THIS VERSION WRIES OUT FILE FOR INPUT TO FIREWORKS DISPLAY PROGRAM
c
c
c
c     previous notes during development history....
c
c
c     NOTE :DIFFERS FROM spkpat6bg:
c     1. jsig=1000 rather that 100 for better sig test.
c     occurrences is > max in any of the 1000 shifted sets
c
c
c
c     PREVIOUS DOCUMENTATION ....
c     NOTE :DIFFERS FROM GRVPAT* IN THAT IT LOOKS BACKWARD IN TIME FROM EACH
c     TEMPLATE COLUMN AS WELL AS FORWARD IN TIME LIKE GRVPAT*!!!!
c
c
c
c
c
c
c     This program detects spatiotemporal patterns of neuronal synchrony that RECUR more
c     frequently than would be expected by chance..
c
c     I. BACKGROUND & OVERVIEW...
c
c     The program currently looks for  repeated single column matches (a): abbbbbbbabbbb; i.e., each
c     single column is considered a "minimum unit to match".
c
c
c     To be developed...
c     b.  repeated 'multicolumn (successive)' patterns(aa): aabbbbbbbbaabbbb
c     c.  repeated 'multicolumn (nonsuccessive)' patterns (a_a): ababbbbbbbbbbabab
c
c
c     II. SOME KEY VARIABLES
c
c     jmptal=number of non-nulled velocities in a template column
c
c     hit = number of "matching" (non-nulled) velocities in template and target columns- i.e.,
c     present for same pairs in both columns.
c
c     null=number of matching "empties" (particle repulsion not yet included) or nulled array elements
c
c
c     twc=temporary storage of absolute diff between a particular ref and target slope
c
c
c     lastr = number of pairs (rows) in data set
c
c     perhits=0.0
c     if (jmptal.gt.0) perhits= hit/float(jmptal) ! WATCH OUT - TWO INT = INT unless float one
c
c     old perxt=ixtra/float(jmptal)
c     now  perxt=ixtra/float(lastr) ! now same fixed num for each col. rather than dep on
c     diff jmptal for each ref col.
c
c     t1 or thresh = % of a template column's non-nulled elements that needed to match in target
c     column to count as match - deals with "missing", non-nulled target column elements
c
c     t2 or thresh2 = % of xtra  allowed in comparison between template and target column
c
c     t3 or thresh3 = max % magnitude difference between a template velocity and target vel. in a
c     column...for Category 7 classification of match - "hardest to achieve" -see next section
c
c     thrnul- default=0.03125 or 1/32.
c     ... IF final array element after scaling with sclfac -derived from'biggest' var-
c     .le. (thrnul*tmpmax) THEN change to 0.0. This makes comparisons match displayed data in
c     xslope32f.
c
c     t4 = % of "hits" that must be IN t3 range to count as match. - new w/ ver. 11
c
c     III. 7  SPECIFIC "MATCH CRITERIA" (Others possible)
c
c
c
c     Category 7
c
c     Exactly same pairs show sig  velocities & magnitudes all matched +/- N%
c     ...... where N set by thresh3 parameter
c
c     if ((perxt.eq.0).and.(perhits.gt.0.99).and.(perwir.gt.t4))
c     + match array element =7
c
c
c     Categories 6-1 see line num 99999,88888,and 77777 for exact logic!!
c
c     IV. SIGNIFICANCE TEST - see iopt var. for two options: ext to control and internal comparison.
c
c     Compares otal (Original tally array) with with ctal (Control tally array), where ctal has biggest
c
c
c     V. PROGRAM SECTIONS
c
c
c     NOTE 1: In the code, "jump" may refer to a specific sig. velocity for a single pair of particles.
c     NOTE 2: The first data array read in (for otal) was generated previously by the program xslope.
c     NOTE 3:  The 100 control data sets are the existing gravity *.pos files used in other sig tests too.
c     ... therefore, the initial processing of each of these (to fill ttally)
c     is different than for the first data
c     set and gives an indication of the processing used in the xslope program
c
c     VERSION 5.0 allows inc values >1. This allows slopes calculated in xslope (and in this program)
c     to represent delta distance over inc PLOTTED time steps. Idea is to look for longer intervals
c     of  associations that repeat.
c
c     version 7 scales each data set to max slope value.. improves comparisons between orig and controls.
c
c     version 9 allows control cycles to do internal rather than with ext - see opt1 setting
c     ..also changed logic of comparison sections 99999,88888,77777
c
c     v10 adds "thresholding" to cut noise, smoothing, and REMOVES internal scaling to 100.
c     NOW  comp. loops uses BIG from sub read as scale for everything.
c
c     v12 now causes ismooth > 0 to cause search to "ignore" next "ismooth" columns in target array
c
c     v13 - various bugs introduced in v12 exterminated. .. it was a case of 1 step forward and 2 back.
c
c     v14 add write out of recurring patterns
c     .. because of possibility of near false positives due to smoothing.
c
c     v15 correct some minor bugs: (fp=i/i) must be fp=i/float(i) for correct result.
c     also added tmpspklst.txt file to go with patfor3d.spk output file in addition to main listing
c     SUMMARY: the program outputs the main text list plus a file for 3djmp use plus a matching file for
c     importation to a spread sheet for quant info on *.spk
c
c     v16 adds MASK option to take advantage of corresponding changes in xslope32g. This option
c     permits the user to optionally MASK out all consideration of data from pairs that were
c     NOT significantly clustered as defined by option 2 in xslope32g. Thus, only pairs that are
c     defined as elements of one or more assemblies in the group by this method are searched for repeating
c     motifs when the MASK option is engaged.
c     ... also changed crit 7,6,5 from .eq.1.0 to .gt. 0.99 for perhits to handle poss rounding errors
c     v17 made 3 changes.
c     1-allows .05 confidence as well as .01 - for 'frameshift' rather than blkshift'
c     ..... gravity sig test data files.
c     2- data for spark plots with steps > 1 can now be written out with last program option.
c     3- changed default for thresh2 to 0.05 from 0.20.
c
c     v. 18 add tally in end of output file for numbet of matches for each "match" criterion.
c     .. also change default for threshold1 to 0.90
c
c     revisions by bgl, Dec. 97
c
c     1. in spkpat2bg changed arrays to allow 1000 -1 sparks, but only 16 diff codes
c     2. in this version, spkpat3bg:
c     -add newsig array for all columns and criteria to be set to 1 if sig repeats..
c     ..newsig(nr,n) where nr=number of columns serving as templates, i.e., 1000-1
c     .. and n = number of criteria, say 7
c     This array is filled as text table of corresponding data is being filled and
c     written to disk (do 8000 ....)
c
c     -array newsig can then be looped through and sig columns and criteria can be used
c     ... as pointers to write out multiple spark and  corresponding text files
c     ... that can be read in and displayed by 3djmp program
c     ++++++++++++++++++++++++++++++++++++++++
c     include for header file derived from gravity for pos file read
c     +++++++++++++++++++++++++++++++++++++++++
c     $INCLUDE 'head16.defs'
c     ------------------------------------------------------------------------
      parameter (ncol=999)      ! number of steps over which slopes are calculated
      parameter (ncat=8)        ! number of categories +1
      parameter (ncell=64,npts=1000,ipairs=2016)

c     **************************************************************************
c     *****************************************************************************
c
c     adjust next line to help define selectivity of search algorithm
c     SEE LABELS 99999 and 88888 for other places to modify "selectivity" and sensitivity of
c     algorithm...

      integer iarnum,ihard,ijmp,ismooth
      real thrnul,tmpmax
c     ncell = max number of cells in data file
c     npts = max num time steps in orig. grav calc.
c     ipairs=number of unique pairs for ncell
      dimension imask(ipairs)
      dimension IRQ(ncol)
      dimension newsig (ncol,ncat-1) ! added in revision spkpat3bg
      integer*4 iseed,iseed2

      real aryin (ipairs,ncol),twc
      real aryin2 (ipairs,ncol)
      integer match2(ncol)
      logical*1 match (ncol,ncol,ncat)
      integer hit
      integer otal(ncol,ncat)
      integer ctal(ncol,ncat)
      integer ttally(ncol,ncat)
      integer ist(ncat-1)       !num of hits for each match category  - for printout
      integer mtal(ncat-1)      !sum of all repeats for each match category criterion - for printout
      real ztal(ncat-1)         ! mean of all repeats for each match category criterion - for printout
      character*30 IFNAM,OFNAM
      character*2 sig(ncat-1)
      character*11 outfilespk, outfiletxt,outfilefwk
c     new BYTE arrays (x,y) where x is first column to have element that repeats in
c     to hold patterns for each match criterion - first used in v. 14
c     elements have 0 or 1-a series of 1s in a row indicates
c     a set of  matches for the corresponding column
      logical*1 cri(ncol,ncol,ncat)

      iarnum=0
      empty=0.0
c     default value of inc = 1! this is used by xslope32 to get slopes or
c     velocities between step x and x+inc!!
c     inc=1!  NOW read in from input file header!!
      ixrange=ncol

c     define seed for random # gen.
      call ransub(iseed)
      iseed2=iseed              ! save first seed for report write-up

c     define sig level
      jsig=1000                 ! other version=100
c     **************************************************
c     read subroutine does the following tasks:
c     I. set up input and output files/filenames; confiles are default names
c     II. fill inary from xslope output file and close file input device
c     **************************************************
      call readx (ifnam,ofnam,iarnum,aryin,inc,NZ,
     +     lastc,t1,t2,t3,t4,ihard,ijmp,iopt,ismooth,thrnul,tmpmax,
     +     imask,maskop,nomsk,sclfac)
      thresh=t1
      thresh2=t2
      thresh3=t3
      lastr = iarnum            ! last row with stuff in main input array
c     **************************************************
c     III. big loop to look for patterns of jumps, each jump set already defined by xslope*
c     *********************************************************

c     mod to version 12 here..
c     variable ifix = 1+value of ismooth variable
c     .. this prevents a ref column from matching an adjacent target
c     column that may match only because of the effect of
c     gaussian smoothing

      ifix=1+ismooth

      call match_loop (aryin) ! fills match and match2 arrays
      cri = match
      
c     *******************************************
c     START OF LOOP TO FILL otal array
c     *******************************************
c     Tally repeating patterns for later comparison with sig test patterns
c     scan each row of the match array (1 less element/row as row number increases)
c     ..... for 1) each category , 2) number of occurrences for each cat., 3) complexity - number of
c     jumps in template i.e. jmptal value for each particular  template column
c     column.
c     Fill each row of otal array with template complexity and #  occurr for category 1-5
c     row1: complexity,  #occur,  #occur, #occur
c     :
c     :
c     row ne: complexity,  #occur,  #occur, #occur

c     Later, compare  with ctal (ncol,ncat) where the rows will
c     have the "biggest values" for jsig shifted data sets.
c     ..ttally holds data from each of the jsig control data sets in turn...
c     if value in ttally bigger than corresponding value in ctal, then the former value replaces
c     the latter.
c     ************************************************************************************

c     fill otal array
      otal = 0
      do nc = 1, lastc
         otal(nc,1)= match2(nc)
         do nr = 1,lastc
            do icat = 1, 7
               if (match(nr,nc,icat)) otal (nc,icat+1)=otal(nc,icat+1)+1
               end do
         end do
      end do
c     *******************************************
c     END OF LOOP TO FILL otal array
c     *******************************************
c     **************************************************
c     **************************************************
c     START OF NEW SIGNIFICANCE TEST ARRAY FILLING LOOP
c     **************************************************
c     fill ctal array by:
c     1.calculating  jsig control shuffled spark sets
c     2. repeating logic in sections above on each - filling ttally array each time
c     3. keep bigest counts for,each category in ctal
c     4. define those elements in otal for which the null hypothesis of no pattern
c     must be rejected.

c     **************************************************
c     calc jsig randomized spark arrays, one at a time.
c     use ctal to hold final results

      ctal = 0

      izz=iarnum                !     define var izz

c     **************************************************
c     DO BIG LOOP DONE JSIG TIMES
c     ******************************************************
      do jaz=1, jsig

c     MOVE each pair row in aryin TO aryin2  IN RANDOM ORDER
         do jv=1,izz            ! outer loop for pair row definition
            IRQ=1
            do J=1,lastc        ! inner loop for stepping through the array columns
               R=RAN(iseed)
               IR = int ((R*lastc)+1)
               do while (IRQ(IR).EQ.0)
                  IR=IR+1
                  IF(IR.GT.lastc)IR=1
               end do
               aryin2(jv,J) = aryin(jv,IR)
               IRQ(IR)=0
            end do
         end do

         call match_loop (aryin2) ! fills match and match2 arrays

c     fill ttally array
         ttally = 0             ! initialize ttally array here
         do nc = 1, lastc
            ttally(nc,1)= match2(nc)
            do nr = 1,lastc
               do icat = 1, 7
                  if (match(nr,nc,icat))
     +                 ttally (nc,icat+1)=ttally(nc,icat+1)+1
               end do
            end do
         end do

c     LOOP TO UPDATE ctal array with any bigger
c     corresponding values from ttally array after each of the 100 loops:
         do nr = 1, lastc-1
            do icat = 1, 8
               if (ttally(nr,icat).gt.ctal(nr,icat))
     +              ctal(nr,icat)=ttally(nr,icat)
            end do
         end do

         write (*, "(2x,I4,' of',I4,' control cycles done.')") jaz,jsig
      end do
c     **************************************************
c     END OF BIG SIGNIFICANCE TEST ARRAY FILLING LOOP..
c     NOW COMPARE otal with ctal for final sig tests:
c     ***************************************************
      do it=1,ncat-1
         ist(it)=0              ! init array for match tallies
         mtal(it)=0             ! init array for repeat tallies
         ztal(it)=0.0           ! init array for repeat tallies
      end do

c     revision bgl Dec 1997:init newsig array for later
c     use in writing out multiple spark files
      newsig = 0

c     open file for report
      OPEN (unit=3,file=ofnam,status='new',FORM='FORMATTED')
      write (3,"(2x,'Generated by program: spkpat6kbg')")
      write (3,"(2X,'Thresh:',f4.2,1x,'Thresh2:',f4.2,1x,'Thresh3:',f4.2,
     +'t4:',f4.2,2x,' jsig :',I4)")
     $     thresh,thresh2,thresh3,t4,jsig
      write (3,"(2X,'gausian smoothing = ',i2,2x, '0,1,2,4... 0 = NONE')") ismooth
      write (3,"(2X,'thrnul = ',f7.5,2x, 'Slopes .le. this % ignored.')") thrnul
      write (3,"(2X, A14,/,
     +2x,'CONTROL OPTION (0=ext-default, 1=int):',I3,/,
     +2x,'sig for xslope op 4&5:',1x,I4,4x, 'sig option for jmp:',I3/,
     $2x,'template', 2x,'complexity',2x, 'Categories 1-7:')") ifnam,iopt,ihard,ijmp
      write (3,"(2X,'MASK option 1 = ON, 0 = OFF: ',i3)") maskop
      write (3,"(2X,'nomsk (if MASK=1, this is num pairs eval',i8)") nomsk
      write (3,"(2X,'total pairs in data',i8)") iarnum
      write (3,"(2X,'RANDOM NUM SEED USED: ',I16)") iseed2

      do nr = 1, lastc-1
         sig = ' '
         do icat = 1, 7
            if((otal(nr,icat+1).gt.ctal(nr,icat+1)).and.(otal(nr,1).gt.1)) then
               ist(icat)=ist(icat)+1
               mtal(icat)=mtal(icat)+otal(nr,icat+1)
               sig(icat)=' *'
               newsig(nr,icat)=1
            end if
         end do
         write (3,"(2X,I3,6X,I3,4X,I3,1X,I3,1X,A2,1x,I3,
     +1X,I3,1x,A2,1X,I3,1X,I3,1x,A2,I3,1x,I3,1x,A2,1x,I3,1x,I3,1x,A2,
     +1x,I3,1x,I3,1x,A2,1x,I3,1x,I3,1x,A2,1x,I3,1x,I3,1x,A2)")
     $        nr,otal(nr,1),otal(nr,2),
     +        ctal(nr,2),
     +        sig1,otal(nr,3),ctal(nr,3),sig2,otal(nr,4),
     +        ctal(nr,4),sig3,otal(nr,5),ctal(nr,5),sig4,
     +        otal(nr,6),ctal(nr,6),sig5,otal(nr,7),
     +        ctal(nr,7),sig6,otal(nr,8),ctal(nr,8),sig7
      end do

      write (3,"(2x,'Matches(2+vec only)- categories 1-7:',/,
     +2X,I5,4X,I5,4X,I5,4X,I5,4X,I5,4X,I5,4X,I5)")
     $     ist(1),ist(2),ist(3),ist(4), ist(5),ist(6),ist(7)

      do icat = 1, 7
         if(ist(icat).ne.0)ztal(icat)=mtal(icat)/float(ist(icat))
      end do

      write (3,"(2x,'mean reps/match(2+vec only)- categories 1-7:',/,
     +2X,f6.1,4X,f6.1,4X,f6.1,4X,f6.1,4X,f6.1,4X,f6.1,4X,f6.1)")
     $     ztal(1),ztal(2),ztal(3), ztal(4),ztal(5),ztal(6), ztal(7)
      close (unit=3)

c     **************************************************
c     option to loop for writing out specific pattern sequence file for 3djmp* input
c     **************************************************
      if (ismooth.ne.0) print "(2x,'Pattern 3-D write illegal with smoothed data')"
      if (ismooth.ne.0) stop
      print "(2x,'Enter 1 to write spk & txt files, <cr> TO QUIT.')"
      print "(/,2x,'*************** WARNING !!! ***************')"
      print "(/,2x,'REMOVE *.spk *.txt & *.fwk files from directory')"
      print "(/,2x,'before choosing option 1 - to avoid overwrite.')"
      read (*,'(I5)') ichoicez
      if (ichoicez.ne.1) stop

c     if sig do you want to write file to be read by 3djmp program?
c     routines to write out ONE recurring pattern in
c     format that 3djmp can read
c     
c     READ CRITERION TO MATCH
      do
         icrit=0                !default = ALL sig. matches - could be alot
         print "(2x,'Enter 1 criterion number;def. is 1-7!!: ',$)"
         print "(2x,'Enter 8 for ALL fireworks ONLY: ',$)"
         read (*,'(I5)') icrit
         if ((icrit.ge.0).and.(icrit.le.8)) exit
      end do

c     ********************************************************************
c     new front end loop control to write out multiple spk and txt files
c     bgl Dec 1997
      do ny=1,ncat-1
         write(outfilefwk,"('f-',I1,'.fwk')")ny ! auto gen of fire works filename (s)
         if ((icrit.ne.0).and.(icrit.ne.ny).and.(icrit.ne.8))cycle
         OPEN (UNIT=4,file=outfilefwk,status='UNKNOWN', access='SEQUENTIAL', form='UNFORMATTED')

         do isel=1,lastc-1
            if (newsig(isel,ny).eq.0) cycle !otherwise USE IT
c     dynamically define descriptive filenames here, format=
c     f<column in part.cond. profile array>-<criterion>.<type>
            write(outfilespk,"('f',I3,'-',I1,'.spk')")isel,ny
            write(outfiletxt,"('f',I3,'-',I1,'.txt')")isel,ny
            aryin2 = 0.0        ! clear aryin2 for reuse
c     A single user defined column 'isel' is used as ref to look at all OTHER cols.
c     Only recurrences of that ONE column, and that column itself,  are written out

c     remember: cri*(m,ne)=1 target column at position m (x-axis) is a match
c     for reference col 'ne' in yaxis.
            
            call row_loop (isel) ! put REFERENCE column elements in array to be written out
            
            do m=1,lastc
               if (cri(m,isel,ny)) then
                  call row_loop (m)  ! put matching TARGET columns elements in array to be written out
               end if
            end do

c     now WRITE AS FILE THIS ONE SEQUENCE OF RECURRING 'EVENTS'
            if (icrit.ne.8) then
               call idl (aryin2,outfilespk,izz,ihard,ijmp,inc,ismooth, imask,biggest)
               call wriary (aryin2,outfiletxt,izz,NZ) !ascii version file for debugging/detailed info
               write (*,"(2x,'Progress... just wrote file: ',A11)") outfilespk
            end if
            call wrfwk (aryin2,isel,ny,izz) ! out to unit 4 A BIT OF currnt fireworks file
         end do
         close (unit=4)         !done with fwk file write for one criterion
      end do

c     END of BIG LOOP TO WRITE OUT SPARK FILES

      CONTAINS

c     ******************************************************************
      SUBROUTINE ROW_LOOP (icol)
c     ******************************************************************
      do jv=1, izz
         aryin2(jv,icol*inc)=aryin(jv,icol)
         if ((inc.gt.1).and.(aryin(jv,icol).gt.0.0)) then
            do new = 1,inc-1
               aryin2(jv,(icol*inc)-new)=aryin(jv,icol)
            end do
         end if
         if((maskop.eq.1).and.(imask(jv).eq.1)) then
            aryin2(jv,icol*inc)=0.0
            if ((inc.gt.1).and.(aryin(jv,icol).gt.0.0)) then
               do new = 1,inc-1
                  aryin2(jv,(icol*inc)-new)=0.0
               end do
            end if
         end if
      end do
      end subroutine row_loop

c     ******************************************************************
      SUBROUTINE MATCH_LOOP (aryin2)
c     ******************************************************************
      real aryin2 (:,:)
      match = .false.
      do ne = 1, lastc
         jmptal=0
         do k = 1, lastr
            if((maskop.eq.1).and.(imask(k).eq.1)) cycle
            if (aryin(k,ne).ne.empty)jmptal=jmptal+1
         end do
         match2(ne)=jmptal

         do m = 0+ifix, lastc   ! inner loop for comparing two columns
            hit=0               ! initialize three "match" category variables for loop
            ixtra=0
            nexr=0

            do j=1,lastr
               if((maskop.eq.1).and.(imask(j).eq.1)) cycle
               if ((aryin(j,ne).gt.empty).and.(aryin2(j,m).gt.empty)) then
                  hit=hit+1
                  twc= (ABS(aryin(j,ne)-aryin2(j,m)))/aryin(j,ne)
                  if (twc.le.thresh3) nexr=nexr+1
               end if
               if ((aryin(j,ne).eq.empty).and.(aryin2(j,m).gt.empty)) then
                  ixtra=ixtra+1
               end if
            end do

c     done with column loop; now summarize and fill location in match array
c     lastr  = total pairs in current data set
c     jmptal, hit and null values available for use...

c     calculate % of possible hits actually achieved

            perhits=0.0
            perxt=0.0
            perwir=0.0
            if (jmptal.gt.0) perhits= hit/float(jmptal)
            perxt=ixtra/float(lastr)
            if (maskop.eq.1) perxt=ixtra/float(nomsk)
            if (hit.gt.0)
     +           perwir= nexr/float(hit) !% magnitudes within range set by thresh3

c     *********************************************
c     THRESHOLD VALUES USED HERE!!
c     *********************************************
            match(m,ne,7)=(perxt.eq.0.0).and.(perhits.gt.0.99).and.(perwir.gt.t4)
            match(m,ne,6)=(perxt.eq.0.0).and.(perhits.gt.0.99)
            match(m,ne,5)=(perxt.lt.thresh2).and.(perhits.gt.0.99)
            match(m,ne,4)=(perhits.ge.thresh).and.(perxt.lt.thresh2).and.(perwir.gt.t4)
            match(m,ne,3)=(perhits.ge.thresh).and.(perxt.lt.thresh2)
            match(m,ne,2)=(perhits.ge.thresh).and.(perwir.gt.t4)
            match(m,ne,1)=(perhits.ge.thresh)
         end do
      end do
      end subroutine match_loop

      end program


c     ******************************************************************
c     *************************** SUBROUTINES **************************
c     ******************************************************************

c     ******************************************************************
      SUBROUTINE READX  (ifnam,ofnam,iarnum,aryin,inc,NZ,
     +     lastc,T1,T2,T3,T4,ihard,ijmp,iopt,ismooth,thrnul,tmpmax,
     +     imask,maskop,nomsk,sclfac)
c     ******************************************************************
     
c     imask = array with row (pair) numbers to use; others skipped
c     maskop = flag for using mask option; 1=use mask option
c     nomsk= value to replace variable 'lastr' when mask option engaged
c     ...but only as denominator in calcs - not as a loop controller
c     ....nomsk is determined here in this routine-see below!!
     
      parameter (ncol=999)      ! number of steps over which slopes are calculated
      parameter (ncell=64,ipairs=2016) ! max pairs for 16 cells
      integer iarnum,ihard,ijmp,iopt
      real aryin (ipairs,ncol)
      real temp (ipairs,ncol)
      dimension imask(ipairs)
      character*30 ifnam,ofnam
         
      print "(2x,'PROGRAM EXPECTS npts = 1000 in *.pos input files')"

      print "(2x,'Input data file name:',$)"
      read (*,'(A)') IFNAM
      print "(2x,'Output data file name:',$)"
      read (*,'(A)') OFNAM
      print "(2x,'Enter 1 for MASK option; def is NO MASK')"
      read (*,"(I3)") maskop
      if (maskop.ne.1) maskop = 0
      print "(2x,'DEF VAL OF thrnul used to null in XSLOPE32')"
      print "(2x,'Enter thrnul (.lt. 1.0; f7.5)...DEF=.03125: ',$)"
      read (*,"(f7.5)") thrnul
      if (thrnul.eq.0.0) thrnul= .03125 ! 1/32 of max to match nulling of
c     graphics in xslope32
      print "(2x,'Enter on 1 line T1,T2,T3,T4 def= .9, .05, .25, .75',/
     +,1x,'Bigger T1 & T4.. Smaller T2 & T3 -> harder match (4f5.2):')"
      read (*,"(4f5.2)") T1,T2,T3,T4
      if (T1.eq.0.0)then
         T1=0.9
         T2=0.05
         T3=0.25
         T4=0.75
      end if
      iopt=0

c     READ INPUT FILE ROUTINE HERE

      OPEN (unit=1,file=IFNAM,status='OLD',access='sequential',
     +     form='UNFORMATTED')
      read (1) ismooth
      read (1) inc
      read (1) ihard
      read (1) ijmp
      read (1) iarnum
      read (1) biggest
      
c     part of mask option revision

      nomsk=0
      do i=1,ipairs
         read (1) imask(i)
         if (imask(i).eq.0) nomsk=nomsk+1
      end do
      write (*,"(2x,'value of nomsk variable = ',i8)") nomsk

      do i=1,iarnum
         do j=1,ncol
            read (1) temp(i,j)
         end do
      end do
c     detour to define NZ, number of neurons represented in input
c     ..all neurons are counted irrespective of mask option status!!
      NZ = (int (sqrt (8 * iarnum + 1.)) + 1) / 2;

      write (*,"(2x,'# of different spike codes: ',I7)")NZ
      close (unit=1)

c     END READ ROUTINE

c     lastc now defined here - depends on value of inc
      lastc=ncol/inc
c     now define global scaling factor derived from xslope32*
c     ... from value of variable "biggest" = largest distance between
c     ...pairs of particles in all 20 or 100 control cycles.
     
      sclfac=100./biggest
     
c     scale here with same value as later for all 100 sh*.pos files!!!
      do i=1,iarnum
         do j=1,ncol
            temp(i,j)=sclfac*temp(i,j)
         end do
      end do
     
c     scale slopes to be used to max used slope
c     LOOKS AT ALL COLUMNS, ignores inc variable
c     ...this search for tmpmax permits use of thrnul in main pgm
c     mod so that max slope looked for has to be for a neuron pair
c     that is NOT MASKED OUT!!
     
      tmpmax=0.0
      do i=1,iarnum
         if ((maskop.eq.1).and.(imask(i).eq.1)) cycle
         do j=1,ncol
            tmpmax=AMAX1(tmpmax,temp(i,j))
         end do
      end do

c     now fill array to be passed back to main
c     ..if inc>1 then get biggest value in range - only makes a diff.
c     .... if smoothing in effect
c     NOTE THAT number of columns is now defined by var iz

      do i=1,iarnum
         iz=0
         do j=1,ncol,inc
            iz=iz+1
            aryin(i,iz)=0.0
            do ka=j,j+inc-1
               tmpka=temp(i,j)
               aryin(i,iz)=AMAX1(tmpka,aryin(i,iz))
            end do
c     next line nulls small values defined by thrnul
            if (aryin(i,iz).le.(thrnul*tmpmax))
     +           aryin(i,iz)=0.0
         end do
      end do

      end subroutine readx
      
c     ******************************************************************
      SUBROUTINE IDL (aryidl,outfile,izz,ihard,ijmp,inc,ismooth, imask,biggest)
c     ******************************************************************
      parameter (npts=1000,ipairs=2016)
c     npts = max num time steps in orig. grav calc.
c     ipairs=number of unique pairs for ncell
      character*11 outfile
      dimension aryidl(ipairs,npts-1)
      dimension imask(ipairs)
c     write out 'idl compatible files' file here - over writes file of same name!!!
      OPEN (UNIT=3,file=outfile,status='UNKNOWN',access='SEQUENTIAL',
     +     form='UNFORMATTED')
      write (3) ismooth         ! number of bins for gaus smooth option
      write (3) inc             ! number of plotted steps over which a slope calculated
      write (3) ihard           ! header to indicate sig test option used to gen data
      write (3) ijmp            ! header to indicate jmp sig test option used to gen data
      write (3) izz
      write (3) biggest
      do id=1,ipairs
         write (3) imask(id)
      end do
      do i=1,izz
         do j=1,npts-1
            write (3) aryidl(i,j)
         end do
      end do
      close (UNIT=3)
      end subroutine idl
      
c     ******************************************************************
      SUBROUTINE WRIARY (aryidl,outfile,izz,NZ)
c     ******************************************************************
      parameter (npts=1000,ipairs=2016)
c     npts = max num time steps in orig. grav calc.
c     ipairs=number of unique pairs for ncell
      character*11 outfile
      dimension aryidl(ipairs,npts-1)
      dimension iunit1(ipairs),iunit2(ipairs)
     
c     fill iunit arrays with neuron pair sequences:
c     in order array column is filled - like trydis , xslope sequence
c     e.g.,
c     2,1
c     3,1
c     3,2
c     4,1
c     4,2
c     4,3
c     .
c     .
c     N,N-1
          
      icount=0
      do i=2,NZ
         do j=1,i-1
            icount=icount+1
            iunit1(icount)=i
            iunit2(icount)=j
            if (icount.eq.izz) exit
         end do
         if (icount.eq.izz) exit
      end do

c     write out 'text file' file here - over writes file of same name!!!
      OPEN (UNIT=3,file=outfile,status='UNKNOWN',access='sequential')
      do j=1,npts-1
         do i=1,izz         ! loop through columns for matching columns
            if (aryidl(i,j).ne.0.0) exit
         end do
         if (i.gt.izz) cycle
         do i=1,izz
            write (3,"(2x,i4,',',i4,',',i4,',',f8.4)") j,iunit1(i),iunit2(i),aryidl(i,j)
         end do
      end do
      close (UNIT=3)
      end subroutine wriary
      
c     ******************************************************************
      SUBROUTINE WRFWK (aryin2,isel,ny,izz)
c     ******************************************************************
c     subroutine to write fireworks file - ONE SET PER CALL TO THIS ROUTINE
      parameter (npts=1000,ipairs=2016)
c     npts = max num time steps in orig. grav calc.
c     ipairs=number of unique pairs for ncell
      dimension aryin2(ipairs,npts-1)
      integer fireworks(npts-1)
     
      mx2 = 0                     ! max vectors per spark in this set
      fireworks = 0
c     write to output fwk file here
      do j=1,npts-1
         maxcount=0
         do i=1,izz         ! loop through columns for matching columns
            if (aryin2(i,j).ne.0.0) then
               fireworks(j)=1
               maxcount=maxcount+1
            end if
         end do
         mx2=MAX0(mx2,maxcount)
      end do
     
c     ONLY WRITE IF THIS SET HAS A SPARK SERIES IN IT...
      if (mx2.gt.0)write (4) isel,ny,mx2,(fireworks(k), k=1,npts-1)
      end subroutine wrfwk

c     ******************************************************************
      SUBROUTINE RANSUB(iseed)
c     ******************************************************************
      integer*4 iseed,t(3),sec
c     define seed for random # gen. transparently to user
      print "(2x,'1<RET> to enter iseed; <RET> for auto gen.')"
      read (*,"(I12)") iopt
      if (iopt.eq.1) then
         print "(2x,'Enter iseed (odd integer, 11 digits max: ',$)"
         read (*,"(I12)") iseed
         return
      end if
      call IDATE(t)             !get 3 integers
      sec=int(SECNDS(0.0))      !seconds since midnight
c     iseed=((sec*month*day*year)*2)+1
      iseed=(((t(1)*t(2)*t(3))/10)*2)+1
      end subroutine ransub
