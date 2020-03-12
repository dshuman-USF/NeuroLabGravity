c       spkpat5bg derived from spkpat4bg - bgl - February 1997
C       THIS VERSION WRIES OUT FILE FOR INPUT TO FIREWORKS DISPLAY PROGRAM
c
c
c
c       previous notes during development history....
c
c
c       NOTE :DIFFERS FROM BKGPAT* IN THAT IT LOOKS at 100 shuffled spark
c       sets, NOT 100 shuffled spike train data sets!!!
c
c
c
c       PREVIOUS DOCUMENTATION ....
c       NOTE :DIFFERS FROM GRVPAT* IN THAT IT LOOKS BACKWARD IN TIME FROM EACH
C       TEMPLATE COLUMN AS WELL AS FORWARD IN TIME LIKE GRVPAT*!!!!
C
C
C
C
C
C
c       This program detects spatiotemporal patterns of neuronal synchrony that RECUR more 
c       frequently than would be expected by chance..
c
c       I. BACKGROUND & OVERVIEW...
c
c       The program currently looks for  repeated single column matches (a): abbbbbbbabbbb; i.e., each
c       single column is considered a "minimum unit to match".
c
c
c       To be developed...
c        b.  repeated 'multicolumn (successive)' patterns(aa): aabbbbbbbbaabbbb
c        c.  repeated 'multicolumn (nonsuccessive)' patterns (a_a): ababbbbbbbbbbabab
c
c
c       II. SOME KEY VARIABLES
c
c       jmptal=number of non-nulled velocities in a template column
c
c       hit = number of "matching" (non-nulled) velocities in template and target columns- i.e., 
c       present for same pairs in both columns.
c
c       null=number of matching "empties" (particle repulsion not yet included) or nulled array elements
c
c
c       twc=temporary storage of absolute diff between a particular ref and target slope
c
c
c       lastr = number of pairs (rows) in data set      
c 
c       perhits=0.0
c       if (jmptal.gt.0) perhits= hit/float(jmptal) ! WATCH OUT - TWO INT = INT unless float one
c  
c       old perxt=ixtra/float(jmptal)
c       now  perxt=ixtra/float(lastr) ! now same fixed num for each col. rather than dep on 
c       diff jmptal for each ref col.
c
c       t1 or thresh = % of a template column's non-nulled elements that needed to match in target 
c       column to count as match - deals with "missing", non-nulled target column elements
c
c       t2 or thresh2 = % of xtra  allowed in comparison between template and target column
c
c       t3 or thresh3 = max % magnitude difference between a template velocity and target vel. in a 
c       column...for Category 7 classification of match - "hardest to achieve" -see next section
c
c       thrnul- default=0.03125 or 1/32.
c               ... IF final array element after scaling with sclfac -derived from'biggest' var-
c       .le. (thrnul*tmpmax) THEN change to 0.0. This makes comparisons match displayed data in
c       xslope32f. 
c
c       t4 = % of "hits" that must be IN t3 range to count as match. - new w/ ver. 11
c
c       III. 7  SPECIFIC "MATCH CRITERIA" (Others possible)
c
c
c
c        Category 7
c
c                Exactly same pairs show sig  velocities & magnitudes all matched +/- N%
c               ...... where N set by thresh3 parameter
c
c       if ((perxt.eq.0).and.(perhits.gt.0.99).and.(perwir.gt.t4)) 
c     + match array element =7
c
c
c        Categories 6-1 see line num 99999,88888,and 77777 for exact logic!!
c
c       IV. SIGNIFICANCE TEST - see iopt var. for two options: ext to control and internal comparison.
c
c       Compares otal (Original tally array) with with ctal (Control tally array), where ctal has biggest
c       
c
c       V. PROGRAM SECTIONS
c
c
c       NOTE 1: In the code, "jump" may refer to a specific sig. velocity for a single pair of particles.
c       NOTE 2: The first data array read in (for otal) was generated previously by the program xslope.
c       NOTE 3:  The 100 control data sets are the existing gravity *.pos files used in other sig tests too.
c       ... therefore, the initial processing of each of these (to fill ttally) 
c       is different than for the first data
c       set and gives an indication of the processing used in the xslope program
c
c       VERSION 5.0 allows inc values >1. This allows slopes calculated in xslope (and in this program)
c       to represent delta distance over inc PLOTTED time steps. Idea is to look for longer intervals
c       of  associations that repeat.
c
c       version 7 scales each data set to max slope value.. improves comparisons between orig and controls.
c
c       version 9 allows control cycles to do internal rather than with ext - see opt1 setting
c       ..also changed logic of comparison sections 99999,88888,77777
c
c       v10 adds "thresholding" to cut noise, smoothing, and REMOVES internal scaling to 100.
c       NOW  comp. loops uses BIG from sub read as scale for everything.
c
c       v12 now causes ismooth > 0 to cause search to "ignore" next "ismooth" columns in target array
c
c       v13 - various bugs introduced in v12 exterminated. .. it was a case of 1 step forward and 2 back.
c
c       v14 add write out of recurring patterns
c       .. because of possibility of near false positives due to smoothing.
c
c       v15 correct some minor bugs: (fp=i/i) must be fp=i/float(i) for correct result.
c               also added tmpspklst.txt file to go with patfor3d.spk output file in addition to main listing
c       SUMMARY: the program outputs the main text list plus a file for 3djmp use plus a matching file for
c       importation to a spread sheet for quant info on *.spk
c
c       v16 adds MASK option to take advantage of corresponding changes in xslope32g. This option
c       permits the user to optionally MASK out all consideration of data from pairs that were
c       NOT significantly clustered as defined by option 2 in xslope32g. Thus, only pairs that are
c       defined as elements of one or more assemblies in the group by this method are searched for repeating
c       motifs when the MASK option is engaged.
c       ... also changed crit 7,6,5 from .eq.1.0 to .gt. 0.99 for perhits to handle poss rounding errors
c       v17 made 3 changes.
c       1-allows .05 confidence as well as .01 - for 'frameshift' rather than blkshift'
c       ..... gravity sig test data files.
c       2- data for spark plots with steps > 1 can now be written out with last program option.
c       3- changed default for thresh2 to 0.05 from 0.20.
c
c       v. 18 add tally in end of output file for numbet of matches for each "match" criterion.
c       .. also change default for threshold1 to 0.90
c
c       revisions by bgl, Dec. 97
c
c       1. in spkpat2bg changed arrays to allow 1000 -1 sparks, but only 16 diff codes
c       2. in this version, spkpat3bg:
c               -add newsig array for all columns and criteria to be set to 1 if sig repeats..
c                       ..newsig(nr,n) where nr=number of columns serving as templates, i.e., 1000-1
c                       .. and n = number of criteria, say 7 
c                       This array is filled as text table of corresponding data is being filled and 
c                       written to disk (do 8000 ....)
c
c               -array newsig can then be looped through and sig columns and criteria can be used 
c                       ... as pointers to write out multiple spark and  corresponding text files 
c                       ... that can be read in and displayed by 3djmp program
c       ++++++++++++++++++++++++++++++++++++++++
c       include for header file derived from gravity for pos file read
c       +++++++++++++++++++++++++++++++++++++++++
C       $INCLUDE 'head16.defs'
C       ------------------------------------------------------------------------
      use surdat_mod
      parameter (ncol=999)      ! number of steps over which slopes are calculated
      parameter (ncat=8)        ! number of categories +1 
      parameter (ncell=64,npts=1000,ipairs=2016)
c
c       **************************************************************************
c       *****************************************************************************
c
c       adjust next line to help define selectivity of search algorithm
c       SEE LABELS 99999 and 88888 for other places to modify "selectivity" and sensitivity of
c       algorithm...
c

c
c
c
      integer iarnum,ihard,ijmp,ismooth 
      real thrnul,tmpmax
c       ncell = max number of cells in data file
c       npts = max num time steps in orig. grav calc.
c       ipairs=number of unique pairs for ncell
      dimension imask(ipairs)
      dimension IRQ(ncol)
      dimension newsig (ncol,ncat-1) ! added in revision spkpat3bg
      integer*4 iseed,iseed2

      real aryin (ipairs,ncol),twc
      real aryin2 (ipairs,ncol)
      integer match (ncol,ncol,ncat),match2(ncol)
      integer hit,null
      integer otal(ncol,ncat)
      integer ctal(ncol,ncat)
      integer ttally(ncol,ncat)
      integer ist(ncat-1)       !num of hits for each match category  - for printout 
      integer mtal(ncat-1)      !sum of all repeats for each match category criterion - for printout 
      real ztal(ncat-1)         ! mean of all repeats for each match category criterion - for printout 
      character*30 IFNAM,OFNAM
      character*2 sig1,sig2,sig3,sig4,sig5,sig6,sig7,newsur
      character*11 outfilespk, outfiletxt,outfilefwk
c       new BYTE arrays (x,y) where x is first column to have element that repeats in 
c       to hold patterns for each match criterion - first used in v. 14
c       elements have 0 or 1-a series of 1s in a row indicates 
c       a set of  matches for the corresponding column
      byte cri1(ncol,ncol)
      byte cri2(ncol,ncol)
      byte cri3(ncol,ncol)
      byte cri4(ncol,ncol)
      byte cri5(ncol,ncol)
      byte cri6(ncol,ncol)
      byte cri7(ncol,ncol)
c
      call getenv ("newsur",newsur)
      iarnum=0
      empty=0.0
c       default value of inc = 1! this is used by xslope32 to get slopes or
c       velocities between step x and x+inc!!
c       inc=1!  NOW read in from input file header!!
      ixrange=ncol
c
c
c
c
c
c
c       define seed for random # gen.
c
c

c
      call ransub(iseed)
      iseed2=iseed              ! save first seed for report write-up
c
c       define sig level
      jsig=100                  ! normally=100 for .01; 20 would = 0.05
c       **************************************************
c       read subroutine does the following tasks:
c       I. set up input and output files/filenames; confiles are default names
c       II. fill inary from xslope output file and close file input device 
c       **************************************************
      call readx (ifnam,ofnam,iarnum,aryin,inc,NZ,
     +     lastc,t1,t2,t3,t4,ihard,ijmp,iopt,ismooth,thrnul,tmpmax,
     +     imask,maskop,nomsk,sclfac)
      thresh=t1
      thresh2=t2
      thresh3=t3
      lastr = iarnum            ! last row with stuff in main input array
c       **************************************************
c       III. big loop to look for patterns of jumps, each jump set already defined by xslope*
c       *********************************************************
c       initialize match array
      do 1900 ia=1,ncol
         do 1910 ib=1,ncol
            do 1920 ic=1,ncat
               match(ia,ib,ic)=0
 1920       continue
 1910    continue
 1900 continue

c       init more arrays
      do 1930 ix=1,ncol
         do 1940 iy=1,ncol
            cri1(ix,iy)=0
            cri2(ix,iy)=0
            cri3(ix,iy)=0
            cri4(ix,iy)=0
            cri5(ix,iy)=0
            cri6(ix,iy)=0
            cri7(ix,iy)=0
 1940    continue
 1930 continue
c
c
c       **********************************************************
c       **********************************************************
c       outer loop comparing contents of column ne with ne + x
c       **********************************************************
c       **********************************************************
c
c       mod to version 12 here..
c               variable ifix = 1+value of ismooth variable
c       .. this prevents a ref column from matching an adjacent target
c       column that may match only because of the effect of 
c       gaussian smoothing
c
      ifix=1+ismooth
c
      do 1000 ne = 1, lastc
         jmptal=0
         do 1300 k = 1, lastr
            if((maskop.eq.1).and.(imask(k).eq.1)) goto 1300
            if (aryin(k,ne).ne.empty)jmptal=jmptal+1
 1300    continue
c       1339    format (i5,i5,f8.4)
         match2(ne)=jmptal
c
c
         do 1100 m = 0+ifix, lastc
c       inner loop for comparing two columns
c       initialize three "match" category variables for loop
            hit=0
            ixtra=0
            null=0
            nexr=0
c
c
c       core algorithm
c

            do 1200, j=1,lastr
               if((maskop.eq.1).and.(imask(j).eq.1)) goto 1200
               if ((aryin(j,ne).eq.empty).and.(aryin(j,m).eq.empty))then
                  null=null+1
               end if
               if ((aryin(j,ne).gt.empty).and.(aryin(j,m).gt.empty))then
                  hit=hit+1
c       note cannot have divide by zero in next line, both > empty(0.0)
                  twc= ((ABS(aryin(j,ne)-aryin(j,m)))/aryin(j,ne))
                  if (twc.le.thresh3) nexr=nexr+1
               end if
               if ((aryin(j,ne).eq.empty).and.(aryin(j,m).gt.empty))then
                  ixtra=ixtra+1
               end if
 1200       continue
c
c       done with column loop; now summarize and fill location in match array
c       lastr  = total pairs in current data set
c       jmptal, hit and null values available for use... 
c
c       calculate % of possible hits actually achieved
c       calculate % of possible xtras actually achieved
c       next 2 lines prevent divide by zero poss.
c       note these "defaults" ..logic must deal with them
            perhits=0.0
            perxt=0.0
            perwir=0.0          !% magnitudes within range set by thresh3
            if (jmptal.gt.0) perhits= hit/float(jmptal)  
            perxt=ixtra/float(lastr)
            if (maskop.eq.1) perxt=ixtra/float(nomsk)
c       write (*,1996) perxt,hit,null,ixtra
c        1996   format (2x,f6.2,i5,i5,i5)
            if (hit.gt.0) 
     +           perwir= nexr/float(hit) !% magnitudes within range set by thresh3
c
c       **********************************************
c       *********************************************
c       THRESHOLD VALUES USED HERE!! and at 88888 below
c       *********************************************
            if ((perxt.eq.0.0).and.(perhits.gt.0.99)
     +           .and.(perwir.gt.t4)) then 
               match(m,ne,7)=7  ! match(x,y)
               cri7(m,ne)=1     ! Column at position m (x-axis) is a match for ref col 'ne'.
            end if
            if ((perxt.eq.0.0).and.(perhits.gt.0.990)) then
               match(m,ne,6)=6  ! match(x,y)
               cri6(m,ne)=1
            end if
            if ((perxt.lt.thresh2).and.(perhits.gt.0.99)) then 
               match(m,ne,5)=5  ! match(x,y)
               cri5(m,ne)=1
            end if
            if ((perhits.ge.thresh).and.(perxt.lt.thresh2)
     +           .and. (perwir.gt.t4)) then
               match(m,ne,4)=4  ! match(x,y)
               cri4(m,ne)=1
            end if
            if ((perhits.ge.thresh).and.(perxt.lt.thresh2)) then 
               match(m,ne,3)=3
               cri3(m,ne)=1
            end if
            if ((perhits.ge.thresh).and.(perwir.gt.t4)) then
               match(m,ne,2)=2
               cri2(m,ne)=1
            end if
            if (perhits.ge.thresh) then 
               match(m,ne,1)=1
               cri1(m,ne)=1
            end if
c       *******************************************
 1100    continue
 1000 continue
c       **************************************************
c       *******************************************
c       END OF LOOP TO FILL match array elements
c       *******************************************
c       *******************************************
c       START OF LOOP TO FILL otal array
c       *******************************************
c       Tally repeating patterns for later comparison with sig test patterns
c       scan each row of the match array (1 less element/row as row number increases)
c       ..... for 1) each category , 2) number of occurrences for each cat., 3) complexity - number of 
c       jumps in template i.e. jmptal value for each particular  template column
c       column. 
c       Fill each row of otal array with template complexity and #  occurr for category 1-5
c       row1: complexity,  #occur,  #occur, #occur
c                       :
c                       :
c       row ne: complexity,  #occur,  #occur, #occur

c       Later, compare  with ctal (ncol,ncat) where the rows will
c       have the "bigest values" for 100 shifted data sets. 
c       ..ttally holds data from each of thew 100 control data sets in turn...
c       if value in ttally bigger than corresponding value in ctal, then the former value replaces 
c       the latter.
c       initialize otal array here
c ************************************************************************************
      do 2000 ij =1,ncat
         do 2100 ik = 1, ncol
c       otal(ik,ij)=0
            otal(ik,ij)=-1      ! to remove SELF MATCH that will occur with BACKWARD SEARCH.
 2100    continue
 2000 continue
c
c       fill otal array
c       
      do 3000 nc = 1, lastc 
         otal(nc,1)= match2(nc)
         do 3100 nr = 1,lastc
            if (match(nr,nc,1).eq. 1) 
     +           otal(nc,2)=otal(nc,2)+1
            if (match(nr,nc,2).eq. 2)
     +           otal (nc,3)=otal(nc,3)+1
            if (match(nr,nc,3).eq. 3) 
     +           otal (nc,4)=otal(nc,4)+1
            if (match(nr,nc,4).eq. 4) 
     +           otal (nc,5)=otal(nc,5)+1
            if (match(nr,nc,5).eq. 5) 
     +           otal (nc,6)=otal(nc,6)+1
            if (match(nr,nc,6).eq. 6) 
     +           otal (nc,7)=otal(nc,7)+1
            if (match(nr,nc,7).eq. 7) 
     +           otal (nc,8)=otal(nc,8)+1
 3100    continue
 3000 continue
c       *******************************************
c       END OF LOOP TO FILL otal array
c       *******************************************
c       **************************************************
c       **************************************************
c       START OF NEW SIGNIFICANCE TEST ARRAY FILLING LOOP
c       **************************************************
c       fill ctal array by:
c       1.calculating  100 control shuffled spark sets
c       2. repeating logic in sections above on each - filling ttally array each time
c       3. keep bigest counts for,each category in ctal
c       4. define those elements in otal for which the null hypothesis of no pattern
c       must be rejected. 

c       **************************************************
c       calc 100 randomized spark arrays, one at a time. 
c       use ctal to hold final results

c       initialize ctal array here
      do 5335 ij =1,ncat
         do 5236 ik = 1, ncol
            ctal(ik,ij)=0
 5236    continue
 5335 continue
c
c       define var izz
      izz=iarnum
c
      if (newsur.eq.'1') call SURDAT_INIT_RC (aryin,izz,lastc)
c
c       **************************************************
c       DO BIG LOOP DONE JSIG TIMES
c       ******************************************************
      do 10000 jaz=1, jsig 
c       reinitialize match array
         do 4900 ia=1,ncol
            do 4910 ib=1,ncol
               do 4920 ic=1,ncat
                  match (ia,ib,ic)=0
 4920          continue
 4910       continue
 4900    continue

         if (newsur.eq.'1') then
            call SURDAT_RC (aryin2,izz,lastc)
         else

c       *******************************************************
c         MOVE each pair row in aryin TO aryin2  IN RANDOM ORDER
c       ********************************************************
c       outer loop for pair row definition
c
         do 3396 jv=1,izz
            do new=1,ncol
               IRQ(new)=1
            end do
c
c       inner loop for stepping through the array columns
c
            J=1
 8500       R=RAN(iseed)
            ANUM=(R*lastc)+1
            IR=int(ANUM)
 8600       IF(IRQ(IR).EQ.0)GOTO 8510
            aryin2(jv,J) = aryin(jv,IR)
            IRQ(IR)=0
            J=J+1
            IF(J.GT.lastc)GOTO 3396
            GOTO 8500
 8510       IR=IR+1
            IF(IR.GT.lastc)IR=1
            GOTO 8600
 3396    continue

         end if

c
c       ***********************************************************
c       END OF ONE  RANDOMIZATION LOOP for one of jsig trials
c       **********************************************************
c
c       write (*, 40915) jaz
c 40915 format (2x, 'End of randomization loop: ', i5)
c
c
c
c
c
c
c       now continue with Big loop to look for patterns 
c       outer loop comparing contents of column n in aryin with n + x in aryin2
c       ***************************************************
c       *******************************************************
         do 4000 ny=1,lastc
            jmptal=0
            do 4300 k = 1,lastr
               if((maskop.eq.1).and.(imask(k).eq.1)) goto 4300
               if (aryin(k,ny).ne.empty)jmptal=jmptal+1
 4300       continue
            match2(ny)=jmptal
            do 4100 m = 0+ifix, lastc
c       inner loop for comparing two columns
c       initialize three "match" category variables for loop
               hit=0
               ixtra=0
               null=0
               nexr=0
c
               do 4200, j=1,lastr
                  if((maskop.eq.1).and.(imask(j).eq.1)) goto 4200
                  if ((aryin(j,ny).eq.empty).and.(aryin2(j,m).eq.empty))
     +                 then
                     null=null+1
                  end if
                  if ((aryin(j,ny).gt.empty).and.(aryin2(j,m).gt.empty))
     +                 then
                     hit=hit+1
                     twc= (ABS(aryin(j,ny)-aryin2(j,m)))/aryin(j,ny)
                     if (twc.le.thresh3) nexr=nexr+1
                  end if 
                  if ((aryin(j,ny).eq.empty).and.(aryin2(j,m).gt.empty))
     +                 then
                     ixtra=ixtra+1
                  end if
 4200          continue
c
c
c       done with column loop; now summarize and fill location in match array
c       lastr  = total pairs in current data set
c       jmptal, hit and null values available for use... 
c
c       calculate % of possible hits achieved
c
               perhits=0.0
               perxt=0.0
               perwir=0.0
               if (jmptal.gt.0) perhits= hit/float(jmptal)  
               perxt=ixtra/float(lastr)
               if(maskop.eq.1) perxt=ixtra/float(nomsk)
               if (hit.gt.0) 
     +              perwir= nexr/float(hit) !% magnitudes within range set by thresh3
c
c       **********************************************
c       **************************************
c       SECOND PLACE WHERE KEY DECISIONS MADE
c       **************************************
c       *********************************************
c       THRESHOLD VALUES USED HERE!! and at 77777 below
c       *********************************************
               if ((perxt.eq.0.0).and.(perhits.gt.0.99)
     +              .and.(perwir.gt.t4)) 
     +              match(m,ny,7)=7 ! match(x,y)
               if ((perxt.eq.0.0).and.(perhits.gt.0.99)) 
     +              match(m,ny,6)=6 ! match(x,y)
               if ((perxt.lt.thresh2).and.(perhits.gt.0.99)) 
     +              match(m,ny,5)=5 ! match(x,y)
               if ((perhits.ge.thresh).and.(perxt.lt.thresh2)
     +              .and. (perwir.gt.t4)) 
     +              match(m,ny,4)=4 ! match(x,y)
               if ((perhits.ge.thresh).and.(perxt.lt.thresh2)) 
     +              match(m,ny,3)=3
               if ((perhits.ge.thresh).and.(perwir.gt.t4)) 
     +              match(m,ny,2)=2
               if (perhits.ge.thresh) 
     +              match(m,ny,1)=1
c       *****************************************
 4100       continue
 4000    continue
c
c
c
c
c
c
c
c       **************************************************
c       initialize ttally array here
c
c
         do 5070 ij =1,ncat
            do 5170 ik = 1, ncol
               ttally(ik,ij)=0
 5170       continue
 5070    continue
c
c
c       fill ttally array
c       
         do 6000 nc = 1, lastc 
            ttally(nc,1)= match2(nc)
            do 6100 nr = 1,lastc
               if (match(nr,nc,1).eq. 1) 
     +              ttally(nc,2)=ttally(nc,2)+1
               if (match(nr,nc,2).eq. 2)
     +              ttally (nc,3)=ttally(nc,3)+1
               if (match(nr,nc,3).eq. 3) 
     +              ttally (nc,4)=ttally(nc,4)+1
               if (match(nr,nc,4).eq. 4) 
     +              ttally (nc,5)=ttally(nc,5)+1
               if (match(nr,nc,5).eq. 5) 
     +              ttally (nc,6)=ttally(nc,6)+1
               if (match(nr,nc,6).eq. 6) 
     +              ttally (nc,7)=ttally(nc,7)+1
               if (match(nr,nc,7).eq. 7) 
     +              ttally (nc,8)=ttally(nc,8)+1
 6100       continue
 6000    continue
c       *******************************************
c       *******************************************
c       END OF LOOP TO FILL ttally array one time
c       *******************************************

c       ************************************************
c       LOOP TO UPDATE ctal array with any bigger 
c       corresponding values from ttally array after each of the 100 loops: 
c       ************************************************
         do 7000 nr = 1, lastc-1 
            if (ttally(nr,1).gt.ctal(nr,1)) ctal(nr,1)=ttally(nr,1)
            if (ttally(nr,2).gt.ctal(nr,2)) ctal(nr,2)=ttally(nr,2)
            if (ttally(nr,3).gt.ctal(nr,3)) ctal(nr,3)=ttally(nr,3)
            if (ttally(nr,4).gt.ctal(nr,4)) ctal(nr,4)=ttally(nr,4)
            if (ttally(nr,5).gt.ctal(nr,5)) ctal(nr,5)=ttally(nr,5)
            if (ttally(nr,6).gt.ctal(nr,6)) ctal(nr,6)=ttally(nr,6)
            if (ttally(nr,7).gt.ctal(nr,7)) ctal(nr,7)=ttally(nr,7)
            if (ttally(nr,8).gt.ctal(nr,8)) ctal(nr,8)=ttally(nr,8)
 7000    continue

c       *************************************************
c       put some user feedback to screen here if this takes a long time..
c       % done , etc
c       *************************************************

         write (*, 40912) jaz
10000 continue
c       **************************************************
c       END OF BIG SIGNIFICANCE TEST ARRAY FILLING LOOP
c       **************************************************
40912 format (2x,I4,' of 100 control cycles done.')
c       **************************************************
c       END OF BIG SIGNIFICANCE TEST ARRAY FILLING LOOP..
c       NOW COMPARE otal with ctal for final sig tests:
c       ***************************************************
c       *************************************************
c       write out file and close output device
c       **************************************************
      do it=1,ncat-1
         ist(it)=0              ! init array for match tallies
         mtal(it)=0             ! init array for repeat tallies
         ztal(it)=0.0           ! init array for repeat tallies
      end do
c
c
c       revision bgl Dec 1997:init newsig array for later 
c               use in writing out multiple spark files
      do ny=1,ncat-1
         do nx=1,ncol
            newsig(nx,ny)=0
         end do
      end do
c
c
c
c       open file for report
      OPEN (unit=3,file=ofnam,status='new',FORM='FORMATTED')
      write (3,8039) 
      write (3,8040) thresh,thresh2,thresh3,t4,jsig 
      write (3,8046) ismooth 
      write (3,8045) thrnul 
      write (3,8010) ifnam,iopt,ihard,ijmp
      write (3,8011) maskop
      write (3,8013) nomsk
      write (3,8014) iarnum
      write (3,8012) iseed2
      do 8000 nr = 1, lastc-1 
         sig1=' '
         sig2=' '
         sig3=' '
         sig4=' '
         sig5=' '
         sig6=' '
         sig7=' '
         if((otal(nr,2).gt.ctal(nr,2)).and.(otal(nr,1).gt.1))
     +        then
            ist(1)=ist(1)+1
            mtal(1)=mtal(1)+otal(nr,2)
            sig1=' *'
            newsig(nr,1)=1
         end if
         if((otal(nr,3).gt.ctal(nr,3)).and.(otal(nr,1).gt.1))
     +        then
            ist(2)=ist(2)+1
            mtal(2)=mtal(2)+otal(nr,3)
            sig2=' *'
            newsig(nr,2)=1
         end if
         if((otal(nr,4).gt.ctal(nr,4)).and.(otal(nr,1).gt.1))
     +        then
            ist(3)=ist(3)+1
            mtal(3)=mtal(3)+otal(nr,4)
            sig3=' *'
            newsig(nr,3)=1
         end if
         if((otal(nr,5).gt.ctal(nr,5)).and.(otal(nr,1).gt.1))
     +        then
            ist(4)=ist(4)+1
            mtal(4)=mtal(4)+otal(nr,5)
            sig4=' *'
            newsig(nr,4)=1
         end if
         if((otal(nr,6).gt.ctal(nr,6)).and.(otal(nr,1).gt.1))
     +        then
            ist(5)=ist(5)+1
            mtal(5)=mtal(5)+otal(nr,6)
            sig5=' *'
            newsig(nr,5)=1
         end if
         if((otal(nr,7).gt.ctal(nr,7)).and.(otal(nr,1).gt.1))
     +        then
            ist(6)=ist(6)+1
            mtal(6)=mtal(6)+otal(nr,7)
            sig6=' *'
            newsig(nr,6)=1
         end if
         if((otal(nr,8).gt.ctal(nr,8)).and.(otal(nr,1).gt.1))
     +        then
            ist(7)=ist(7)+1
            mtal(7)=mtal(7)+otal(nr,8)
            sig7=' *'
            newsig(nr,7)=1
         end if
c
         write (3,8020) nr,otal(nr,1),otal(nr,2),
     +        ctal(nr,2),
     +        sig1,otal(nr,3),ctal(nr,3),sig2,otal(nr,4),
     +        ctal(nr,4),sig3,otal(nr,5),ctal(nr,5),sig4,
     +        otal(nr,6),ctal(nr,6),sig5,otal(nr,7),
     +        ctal(nr,7),sig6,otal(nr,8),ctal(nr,8),sig7
 8000 continue
      write (3,8035) ist(1),ist(2),ist(3),ist(4),
     +     ist(5),ist(6),ist(7)
      if(ist(1).ne.0)ztal(1)=mtal(1)/float(ist(1))
      if(ist(2).ne.0)ztal(2)=mtal(2)/float(ist(2))
      if(ist(3).ne.0)ztal(3)=mtal(3)/float(ist(3))
      if(ist(4).ne.0)ztal(4)=mtal(4)/float(ist(4))
      if(ist(5).ne.0)ztal(5)=mtal(5)/float(ist(5))
      if(ist(6).ne.0)ztal(6)=mtal(6)/float(ist(6))
      if(ist(7).ne.0)ztal(7)=mtal(7)/float(ist(7))
      write (3,8036) ztal(1),ztal(2),ztal(3),
     +     ztal(4),ztal(5),ztal(6),
     +     ztal(7)
      close (unit=3)

c       **************************************************
c       option to loop for writing out specific pattern sequence file for 3djmp* input
c       **************************************************
      if (ismooth.ne.0) print 3160
      if (ismooth.ne.0) goto 3159
      print 3158
 3158 format (2x,'Enter 1 to write spk & txt files, <cr> TO QUIT.')
      print 3155
 3155 format (/,2x,'*************** WARNING !!! ***************')
      print 3152
 3152 format (/,2x,'REMOVE *.spk *.txt & *.fwk files from directory')
      print 3151
 3151 format (/,2x,'before choosing option 1 - to avoid overwrite.')
 3160 format(2x,'Pattern 3-D write illegal with smoothed data')
      read (*,3960) ichoicez
      if (ichoicez.ne.1) goto 3159
c
c       if sig do you want to write file to be read by 3djmp program?
c       routines to write out ONE recurring pattern in
c       format that 3djmp can read
c
c       READ CRITERION TO MATCH
 3157 icrit=0                   !default = ALL sig. matches - could be alot
      print 3156
 3156 format (2x,'Enter 1 criterion number;def. is 1-7!!: ',$)
      print 3336
 3336 format (2x,'Enter 8 for ALL fireworks ONLY: ',$)
      read (*,3960) icrit
      if ((icrit.lt.0).or.(icrit.gt.8)) goto 3157
c
c
c       ********************************************************************
c       new front end loop control to write out multiple spk and txt files
c       bgl Dec 1997
      do 5600, ny=1,ncat-1
         write(outfilefwk,5540)ny ! auto gen of fire works filename (s)
         if ((icrit.eq.0).or.(icrit.eq.ny).or.(icrit.eq.8))then
            OPEN (UNIT=4,file=outfilefwk,status='UNKNOWN',
     +           access='SEQUENTIAL',
     +           form='UNFORMATTED')
            goto 5501
         else
            goto 5600
         end if
c
 5501    do 5500, isel=1,lastc-1
c
            if (newsig(isel,ny).eq.0)goto 5500 !otherwise USE IT
c
c       dynamically define descriptive filenames here, format= 
c       f<column in part.cond. profile array>-<criterion>.<type>

            write(outfilespk,5520)isel,ny
            write(outfiletxt,5530)isel,ny

 5520       format('f',I3,'-',I1,'.spk')
 5530       format('f',I3,'-',I1,'.txt')
 5540       format('f-',I1,'.fwk')

c
c
c
c       clear aryin2 for reuse
            do jv=1,ipairs
               do kv=1,ncol
                  aryin2(jv,kv)=0.0
               end do
            end do
c       A single user defined column 'isel' is used as ref to look at all OTHER cols.
c       Only recurrences of that ONE column, and that column itself,  are written out

c       remember: cri*(m,ne)=1 target column at position m (x-axis) is a match 
c       for reference col 'ne' in yaxis.
c
c
            if (ny.eq.1) then
               do jv=1, izz
                  aryin2(jv,isel*inc)=aryin(jv,isel)
                  if ((inc.gt.1).and.(aryin(jv,isel).gt.0.0)) then
                     do new = 1,inc-1
                        aryin2(jv,(isel*inc)-new)=aryin(jv,isel)
                     end do
                  end if
                  if((maskop.eq.1).and.(imask(jv).eq.1)) then
                     aryin2(jv,isel*inc)=0.0
                     if ((inc.gt.1).and.(aryin(jv,isel).gt.0.0)) then
                        do new = 1,inc-1
                           aryin2(jv,(isel*inc)-new)=0.0
                        end do
                     end if
                  end if
               end do 
c       REFERENCE column elements now in array to be written out.
c
               do m=1,lastc
                  if (cri1(m,isel).eq.1 ) then
                     do jv=1, izz
                        aryin2(jv,m*inc)=aryin(jv,m)
                        if ((inc.gt.1).and.(aryin(jv,m).gt.0.0)) then
                           do new = 1,inc-1
                              aryin2(jv,(m*inc)-new)=aryin(jv,m)
                           end do
                        end if
                        if((maskop.eq.1).and.(imask(jv).eq.1)) then
                           aryin2(jv,m*inc)=0.0
                           if ((inc.gt.1).and.(aryin(jv,m).gt.0.0)) then
                              do new = 1,inc-1
                                 aryin2(jv,(m*inc)-new)=0.0
                              end do
                           end if
                        end if
                     end do 
                  end if
               end do
            end if
c       matching TARGET columns elements now in array to be written out.
c
c
            if (ny.eq.2) then
               do jv=1, izz
                  aryin2(jv,isel*inc)=aryin(jv,isel)
                  if ((inc.gt.1).and.(aryin(jv,isel).gt.0.0)) then
                     do new = 1,inc-1
                        aryin2(jv,(isel*inc)-new)=aryin(jv,isel)
                     end do
                  end if
                  if((maskop.eq.1).and.(imask(jv).eq.1)) then
                     aryin2(jv,isel*inc)=0.0
                     if ((inc.gt.1).and.(aryin(jv,isel).gt.0.0)) then
                        do new = 1,inc-1
                           aryin2(jv,(isel*inc)-new)=0.0
                        end do
                     end if
                  end if
               end do 
c       REFERENCE column elements now in array to be written out.
c
               do m=1,lastc
                  if (cri2(m,isel).eq.1 ) then
                     do jv=1, izz
                        aryin2(jv,m*inc)=aryin(jv,m)
                        if ((inc.gt.1).and.(aryin(jv,m).gt.0.0)) then
                           do new = 1,inc-1
                              aryin2(jv,(m*inc)-new)=aryin(jv,m)
                           end do
                        end if
                        if((maskop.eq.1).and.(imask(jv).eq.1)) then
                           aryin2(jv,m*inc)=0.0
                           if ((inc.gt.1).and.(aryin(jv,m).gt.0.0)) then
                              do new = 1,inc-1
                                 aryin2(jv,(m*inc)-new)=0.0
                              end do
                           end if
                        end if
                     end do 
                  end if
               end do
            end if
c       matching TARGET columns elements now in array to be written out.
c
c
c
c
            if (ny.eq.3) then
               do jv=1, izz
                  aryin2(jv,isel*inc)=aryin(jv,isel)
                  if ((inc.gt.1).and.(aryin(jv,isel).gt.0.0)) then
                     do new = 1,inc-1
                        aryin2(jv,(isel*inc)-new)=aryin(jv,isel)
                     end do
                  end if
                  if((maskop.eq.1).and.(imask(jv).eq.1)) then
                     aryin2(jv,isel*inc)=0.0
                     if ((inc.gt.1).and.(aryin(jv,isel).gt.0.0)) then
                        do new = 1,inc-1
                           aryin2(jv,(isel*inc)-new)=0.0
                        end do
                     end if
                  end if
               end do 
c       REFERENCE column elements now in array to be written out.
c
               do m=1,lastc
                  if (cri3(m,isel).eq.1 ) then
                     do jv=1, izz
                        aryin2(jv,m*inc)=aryin(jv,m)
                        if ((inc.gt.1).and.(aryin(jv,m).gt.0.0)) then
                           do new = 1,inc-1
                              aryin2(jv,(m*inc)-new)=aryin(jv,m)
                           end do
                        end if
                        if((maskop.eq.1).and.(imask(jv).eq.1)) then
                           aryin2(jv,m*inc)=0.0
                           if ((inc.gt.1).and.(aryin(jv,m).gt.0.0)) then
                              do new = 1,inc-1
                                 aryin2(jv,(m*inc)-new)=0.0
                              end do
                           end if
                        end if
                     end do 
                  end if
               end do
            end if
c       matching TARGET columns elements now in array to be written out.
c
c
c
            if (ny.eq.4) then
               do jv=1, izz
                  aryin2(jv,isel*inc)=aryin(jv,isel)
                  if ((inc.gt.1).and.(aryin(jv,isel).gt.0.0)) then
                     do new = 1,inc-1
                        aryin2(jv,(isel*inc)-new)=aryin(jv,isel)
                     end do
                  end if
                  if((maskop.eq.1).and.(imask(jv).eq.1)) then
                     aryin2(jv,isel*inc)=0.0
                     if ((inc.gt.1).and.(aryin(jv,isel).gt.0.0)) then
                        do new = 1,inc-1
                           aryin2(jv,(isel*inc)-new)=0.0
                        end do
                     end if
                  end if
               end do 
c       REFERENCE column elements now in array to be written out.
c
               do m=1,lastc
                  if (cri4(m,isel).eq.1 ) then
                     do jv=1, izz
                        aryin2(jv,m*inc)=aryin(jv,m)
                        if ((inc.gt.1).and.(aryin(jv,m).gt.0.0)) then
                           do new = 1,inc-1
                              aryin2(jv,(m*inc)-new)=aryin(jv,m)
                           end do
                        end if
                        if((maskop.eq.1).and.(imask(jv).eq.1)) then
                           aryin2(jv,m*inc)=0.0
                           if ((inc.gt.1).and.(aryin(jv,m).gt.0.0)) then
                              do new = 1,inc-1
                                 aryin2(jv,(m*inc)-new)=0.0
                              end do
                           end if
                        end if
                     end do 
                  end if
               end do
            end if
c       matching TARGET columns elements now in array to be written out.
c
c
c
            if (ny.eq.5) then
               do jv=1, izz
                  aryin2(jv,isel*inc)=aryin(jv,isel)
                  if ((inc.gt.1).and.(aryin(jv,isel).gt.0.0)) then
                     do new = 1,inc-1
                        aryin2(jv,(isel*inc)-new)=aryin(jv,isel)
                     end do
                  end if
                  if((maskop.eq.1).and.(imask(jv).eq.1)) then
                     aryin2(jv,isel*inc)=0.0
                     if ((inc.gt.1).and.(aryin(jv,isel).gt.0.0)) then
                        do new = 1,inc-1
                           aryin2(jv,(isel*inc)-new)=0.0
                        end do
                     end if
                  end if
               end do 
c       REFERENCE column elements now in array to be written out.
c
               do m=1,lastc
                  if (cri5(m,isel).eq.1 ) then
                     do jv=1, izz
                        aryin2(jv,m*inc)=aryin(jv,m)
                        if ((inc.gt.1).and.(aryin(jv,m).gt.0.0)) then
                           do new = 1,inc-1
                              aryin2(jv,(m*inc)-new)=aryin(jv,m)
                           end do
                        end if
                        if((maskop.eq.1).and.(imask(jv).eq.1)) then
                           aryin2(jv,m*inc)=0.0
                           if ((inc.gt.1).and.(aryin(jv,m).gt.0.0)) then
                              do new = 1,inc-1
                                 aryin2(jv,(m*inc)-new)=0.0
                              end do
                           end if
                        end if
                     end do 
                  end if
               end do
            end if
c       matching TARGET columns elements now in array to be written out.
c
c
c
c
            if (ny.eq.6) then
               do jv=1, izz
                  aryin2(jv,isel*inc)=aryin(jv,isel)
                  if ((inc.gt.1).and.(aryin(jv,isel).gt.0.0)) then
                     do new = 1,inc-1
                        aryin2(jv,(isel*inc)-new)=aryin(jv,isel)
                     end do
                  end if
                  if((maskop.eq.1).and.(imask(jv).eq.1)) then
                     aryin2(jv,isel*inc)=0.0
                     if ((inc.gt.1).and.(aryin(jv,isel).gt.0.0)) then
                        do new = 1,inc-1
                           aryin2(jv,(isel*inc)-new)=0.0
                        end do
                     end if
                  end if
               end do 
c       REFERENCE column elements now in array to be written out.
c
               do m=1,lastc
                  if (cri6(m,isel).eq.1 ) then
                     do jv=1, izz
                        aryin2(jv,m*inc)=aryin(jv,m)
                        if ((inc.gt.1).and.(aryin(jv,m).gt.0.0)) then
                           do new = 1,inc-1
                              aryin2(jv,(m*inc)-new)=aryin(jv,m)
                           end do
                        end if
                        if((maskop.eq.1).and.(imask(jv).eq.1)) then
                           aryin2(jv,m*inc)=0.0
                           if ((inc.gt.1).and.(aryin(jv,m).gt.0.0)) then
                              do new = 1,inc-1
                                 aryin2(jv,(m*inc)-new)=0.0
                              end do
                           end if
                        end if
                     end do 
                  end if
               end do
            end if
c       matching TARGET columns elements now in array to be written out.
c
c
c
            if (ny.eq.7) then
               do jv=1, izz
                  aryin2(jv,isel*inc)=aryin(jv,isel)
                  if ((inc.gt.1).and.(aryin(jv,isel).gt.0.0)) then
                     do new = 1,inc-1
                        aryin2(jv,(isel*inc)-new)=aryin(jv,isel)
                     end do
                  end if
                  if((maskop.eq.1).and.(imask(jv).eq.1)) then
                     aryin2(jv,isel*inc)=0.0
                     if ((inc.gt.1).and.(aryin(jv,isel).gt.0.0)) then
                        do new = 1,inc-1
                           aryin2(jv,(isel*inc)-new)=0.0
                        end do
                     end if
                  end if
               end do 
c       REFERENCE column elements now in array to be written out.
c
               do m=1,lastc
                  if (cri7(m,isel).eq.1 ) then
                     do jv=1, izz
                        aryin2(jv,m*inc)=aryin(jv,m)
                        if ((inc.gt.1).and.(aryin(jv,m).gt.0.0)) then
                           do new = 1,inc-1
                              aryin2(jv,(m*inc)-new)=aryin(jv,m)
                           end do
                        end if
                        if((maskop.eq.1).and.(imask(jv).eq.1)) then
                           aryin2(jv,m*inc)=0.0
                           if ((inc.gt.1).and.(aryin(jv,m).gt.0.0)) then
                              do new = 1,inc-1
                                 aryin2(jv,(m*inc)-new)=0.0
                              end do
                           end if
                        end if
                     end do 
                  end if
               end do
            end if
c       matching TARGET columns elements now in array to be written out.
c
c
c
c       TARGET column elements now in aryin2 to be written out.
c       now WRITE AS FILE THIS ONE SEQUENCE OF RECURRING 'EVENTS'
c
            if (icrit.eq.8) goto 5998
            call idl (aryin2,outfilespk,izz,ihard,ijmp,inc,ismooth,
     +           imask,biggest)
            call wriary (aryin2,outfiletxt,izz,NZ) !ascii version file for debugging/detailed info
            write (*,55501) outfilespk
 5998       call wrfwk (aryin2,isel,ny,izz) ! out to unit 4 A BIT OF currnt fireworks file 
 5500    continue
         close (unit=4)         !done with fwk file write for one criterion
 5600 continue
c
c
c       END of BIG LOOP TO WRITE OUT SPARK FILES
c
c
 3960 format (I5)
 8010 FORMAT(2X, A14,/,
     +     2x,'CONTROL OPTION (0=ext-default, 1=int):',I3,/,
     +     2x,'sig for xslope op 4&5:',1x,I4,4x,
     +     'sig option for jmp:',I3/,2x,'template',
     +     2x,'complexity',2x, 'Categories 1-7:')
 8020 FORMAT(2X,I3,6X,I3,4X,I3,1X,I3,1X,A2,1x,I3, 1X,I3,1x,A2,1X,I3,1X,
     +     I3,1x,A2,I3,1x,I3,1x,A2,1x,I3,1x,I3,1x,A2,
     +     1x,I3,1x,I3,1x,A2,1x,I3,1x,I3,1x,A2,1x,I3,1x,I3,1x,A2)
 8039 format (2x,'Generated by program: spkpat1')
 8040 FORMAT(2X,'Thresh:',f4.2,1x,'Thresh2:',f4.2,1x,'Thresh3:',f4.2,
     +     't4:',f4.2,2x,' jsig (20=.05, 100=.01):',I4)
 8045 FORMAT(2X,'thrnul = ',f7.5,2x, 'Slopes .le. this % ignored.')
 8046 FORMAT(2X,'gausian smoothing = ',i2,2x, '0,1,2,4... 0 = NONE')
 8011 FORMAT(2X,'MASK option 1 = ON, 0 = OFF: ',i3)
 8013 FORMAT(2X,'nomsk (if MASK=1, this is num pairs eval',i8)
 8014 FORMAT(2X,'total pairs in data',i8)
 8012 FORMAT(2X,'RANDOM NUM SEED USED: ',I16)
55501 format(2x,'Progress... just wrote file: ',A11)
c
c
 8035 FORMAT(2x,'Matches(2+vec only)- categories 1-7:',/,
     +     2X,I5,4X,I5,4X,I5,4X,I5,4X,I5,4X,I5,4X,I5)
 8036 FORMAT(2x,'mean reps/match(2+vec only)- categories 1-7:',/,
     +     2X,f6.1,4X,f6.1,4X,f6.1,4X,f6.1,4X,f6.1,4X,f6.1,4X,f6.1)
 3159 end
c
c
c
c
c       ********************************************************************
c       **********************************************************************
c       *************************** SUBROUTINES ***************************
c       *******************************************************************
c
      subroutine readx  (ifnam,ofnam,iarnum,aryin,inc,NZ,
     +     lastc,T1,T2,T3,T4,ihard,ijmp,iopt,ismooth,thrnul,tmpmax,
     +     imask,maskop,nomsk,sclfac)
c
c       imask = array with row (pair) numbers to use; others skipped
c       maskop = flag for using mask option; 1=use mask option
c       nomsk= value to replace variable 'lastr' when mask option engaged
c               ...but only as denominator in calcs - not as a loop controller
c               ....nomsk is determined here in this routine-see below!!
c
      parameter (ncol=999)      ! number of steps over which slopes are calculated
      parameter (ncell=16,ipairs=2016) ! max pairs for 16 cells
      integer iarnum,ihard,ijmp,iopt
      real aryin (ipairs,ncol)
      real temp (ipairs,ncol)
      dimension imask(ipairs)
      character*30 ifnam,ofnam
c
c
c
      print 51
 51   format (2x,'PROGRAM EXPECTS npts = 1000 in *.pos input files')
c       
c       
      print 100
      read (*,'(A)') IFNAM
      print 400
      read (*,'(A)') OFNAM
      print 755
 755  format (2x,'Enter 1 for MASK option; def is NO MASK')
      read (*,560) maskop
      if (maskop.ne.1) maskop = 0
      print 155
      print 156
 155  format (2x,'DEF VAL OF thrnul used to null in XSLOPE32')
 156  format (2x,'Enter thrnul (.lt. 1.0; f7.5)...DEF=.03125: ',$)
      read (*,960) thrnul
      if (thrnul.eq.0.0) thrnul= .03125 ! 1/32 of max to match nulling of
c       graphics in xslope32
      print 150
 150  format (2x,'Enter on 1 line T1,T2,T3,T4 def= .9, .05, .25, .75',/,
     +     1x,
     +     'Bigger T1 & T4.. Smaller T2 & T3 -> harder match (4f5.2):')
      read (*,950) T1,T2,T3,T4
      if (T1.eq.0.0)then
         T1=0.9
         T2=0.05
         T3=0.25
         T4=0.75
      end if
      iopt=0
c
c
c
c       READ INPUT FILE ROUTINE HERE
c
c
c
      OPEN (unit=1,file=IFNAM,status='OLD',access='sequential',
     +     form='UNFORMATTED')
      read (1) ismooth
      read (1) inc
      read (1) ihard
      read (1) ijmp
      read (1) iarnum
      read (1) biggest
c
c       part of mask option revision
c
      nomsk=0
      do 7000 i=1,ipairs
         read (1) imask(i)
         if (imask(i).eq.0) nomsk=nomsk+1
 7000 continue
      write (*,565) nomsk
c
c
c
      do 1000 i=1,iarnum
         do 1001 j=1,ncol
            read (1) temp(i,j)
 1001    continue
 1000 continue
c       detour to define NZ, number of neurons represented in input
c       ..all neurons are counted irrespective of mask option status!!
c

 560  format(i3)
 565  format(2x,'value of nomsk variable = ',i8)
      itot=1
      do 2000, il=3,ncell
         itot=itot+il-1
         if(itot.eq.iarnum)NZ=il
         if(itot.eq.iarnum)goto 2010
 2000 continue
c       ok, NZ defined....
 2010 write (*,500)NZ
c
      close (unit=1)
c
c
c       END READ ROUTINE
c
c
c
c
c       lastc now defined here - depends on value of inc
      lastc=ncol/inc
c       now define global scaling factor derived from xslope32*
c       ... from value of variable "biggest" = largest distance between
c       ...pairs of particles in all 20 or 100 control cycles.
c
      sclfac=100./biggest
c       
c       scale here with same value as later for all 100 sh*.pos files!!!
      do i=1,iarnum
         do j=1,ncol
            temp(i,j)=sclfac*temp(i,j)
         end do
      end do
c
c       scale slopes to be used to max used slope
c       LOOKS AT ALL COLUMNS, ignores inc variable
c       ...this search for tmpmax permits use of thrnul in main pgm
c       mod so that max slope looked for has to be for a neuron pair
c       that is NOT MASKED OUT!!
c
      tmpmax=0.0
      do 5000 i=1,iarnum
         if ((maskop.eq.1).and.(imask(i).eq.1))goto 5000
         do 6000 j=1,ncol
            tmpmax=AMAX1(tmpmax,temp(i,j))
 6000    continue
 5000 continue
c
c
c
c       now fill array to be passed back to main
c       ..if inc>1 then get biggest value in range - only makes a diff.
c       .... if smoothing in effect
c       NOTE THAT number of columns is now defined by var iz
c
      do 3000 i=1,iarnum
         iz=0
         do 4000 j=1,ncol,inc
            iz=iz+1
            aryin(i,iz)=0.0
            do 4500 ka=j,j+inc-1
               tmpka=temp(i,j)
               aryin(i,iz)=AMAX1(tmpka,aryin(i,iz))
 4500       continue
c       next line nulls small values defined by thrnul
            if (aryin(i,iz).le.(thrnul*tmpmax))
     +           aryin(i,iz)=0.0
 4000    continue
 3000 continue
c
c
 100  format (2x,'Input data file name:',$)
 400  format (2x,'Output data file name:',$)
 500  format (2x,'# of different spike codes: ',I7)
 950  format(4f5.2)
 960  format(f7.5)
      return
      end
c       *****************************************************
      subroutine idl (aryidl,outfile,izz,ihard,ijmp,inc,ismooth,
     +     imask,biggest)
      parameter (npts=1000,ipairs=2016) 
c       npts = max num time steps in orig. grav calc.
c       ipairs=number of unique pairs for ncell
      character*11 outfile
      dimension aryidl(ipairs,npts-1)
      dimension imask(ipairs)
c       write out 'idl compatible files' file here - over writes file of same name!!!
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
      return
      end
c       *************************************************************
c       *****************************************************
      subroutine wriary (aryidl,outfile,izz,NZ)
      parameter (npts=1000,ipairs=2016) 
c       npts = max num time steps in orig. grav calc.
c       ipairs=number of unique pairs for ncell
      character*11 outfile
      dimension aryidl(ipairs,npts-1)
      dimension iunit1(ipairs),iunit2(ipairs)
c
c       fill iunit arrays with neuron pair sequences:
c       in order array column is filled - like trydis , xslope sequence
c       e.g.,
c       2,1
c       3,1
c       3,2
c       4,1
c       4,2
c       4,3
c       .
c       .
c       N,N-1
c
c
      icount=0
      do 20 i=2,NZ
         do 10 j=1,i-1
            icount=icount+1
            iunit1(icount)=i
            iunit2(icount)=j
            if (icount.eq.izz)goto 30
 10      continue
 20   continue
c
c
c
c       write out 'text file' file here - over writes file of same name!!!
 30   OPEN (UNIT=3,file=outfile,status='UNKNOWN',access='sequential')
      do 100 j=1,npts-1
         do 200 i=1,izz         ! loop through columns for matching columns
            if (aryidl(i,j).ne.0.0) goto 800
 200     continue
         goto 100
 800     do 300 i=1,izz
            write (3,500) j,iunit1(i),iunit2(i),aryidl(i,j)
 300     continue
 100  continue
 500  format (2x,i4,',',i4,',',i4,',',f8.4)
      close (UNIT=3)
      return
      end
c       *************************************************************
c       *************************************************************
c
c       subroutine to write fireworks file - ONE SET PER CALL TO THIS ROUTINE
      subroutine wrfwk (aryin2,isel,ny,izz)
      parameter (npts=1000,ipairs=2016) 
c       npts = max num time steps in orig. grav calc.
c       ipairs=number of unique pairs for ncell
      dimension aryin2(ipairs,npts-1)
      integer fireworks(npts-1)
c
      mx2=0                     ! max vectors per spark in this set       
c
c
      do lz=1,npts-1
         fireworks(lz)=0
      end do
c
c       write to output fwk file here 
      do 100 j=1,npts-1
         maxcount=0
         do 200 i=1,izz         ! loop through columns for matching columns
            if (aryin2(i,j).ne.0.0) then
               fireworks(j)=1
               maxcount=maxcount+1
            end if
 200     continue
         mx2=MAX0(mx2,maxcount)
         
 100  continue
c
c       ONLY WRITE IF THIS SET HAS A SPARK SERIES IN IT...
      if (mx2.gt.0)write (4) isel,ny,mx2,(fireworks(k),
     +     k=1,npts-1)
      return
      end
c
c
c
c
c       *******************************************************************
      subroutine ransub(iseed)
      integer*4 iseed,t(3),sec
c       define seed for random # gen. transparently to user
c
c
      print 755
 755  format (2x,'1<RET> to enter iseed; <RET> for auto gen.')
      read (*,960) iopt
      if (iopt.eq.1) then
         print 156
 156     format (2x,'Enter iseed (odd integer, 11 digits max: ',$)
         read (*,960) iseed
 960     format(I12)

         goto 200
      end if

      call IDATE(t)             !get 3 integers
      sec=int(SECNDS(0.0))           !seconds since midnight
c       iseed=((sec*month*day*year)*2)+1
      iseed=(((t(1)*t(2)*t(3))/10)*2)+1
      call pcg32_srandom(iseed, 54);

 200  return
      end
