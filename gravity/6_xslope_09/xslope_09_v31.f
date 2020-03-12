      program xslope_04
c       latest version of xslope bgl June 2009
C
c       *********************************************************
      INCLUDE 'head_04.def'
c       *************************************************************
C       DESECRIPTION: Plots "time vs. distance" for each pair in a
c       *.pos type gravity file together with control plots in 100, 1000 or 10K
c       existing sh*.pos files. Control plots are plotted in a different color
c       and line type. A "band" defining the max and min of each control
c       time step can be optionally graphed. 
c
c       v10.0 bgl April 1, 1989 ... first version derived from trydis10.f
c       v10.2 bgl April 3, 1989 ... add mean distance plot
c       v10.3 bgl for X windows using sox11 starbase driver
c       v11.0 bgl 2-dec-91 for X windows using sox11 starbase driver
c
c       v11 
c       1. various plots and normalizations of particle velocity
c       for all pairs.
c
c       2. option of reading in single *.pos record with out
c       
c       bgl 1-8-95:
c       xslope3 - modify color tables in newplot2 call newplot2
c       with arg iopt: 1 - cool blue to red to white
c                       2 - spectrum
c                       3 - gray scale
c
c       bgl 1-10-95 xslope4
c       modify filling of "slope arrays" so that rates of aggregation over multiple time
c       steps can be calculated... e.g., slope over every three steps will result in "wider bins.
c       call newplot3
c       bgl 1-17-95 xslope5
c       modify so histograms for each spike train show sum of slope values for all pairs involving
c       that cell at that step
c
c       xslope6 bgl allows up to 100 control files
c       xslope7 bgl corrects N variable in loop control and corrects messages
c       xslope10 bgl histograms now for all pairs rather than just cells
c       xslope11 bgl 5/2/95
c       DONE: last display set in a run shows jumps - where at 
c       a given time step or range at least 3 pairs must exceed threshold of max.
c       aggregation velocity found in 100 or 1000 shifted control sets. 
c        2) writes out array for times when  particles are significantly close: whensig.out
c       DONE: 3) writes out array for times when pairs of move closer : idlmov.out
c       DONE: 3) writes out array for times when pairs aggregation vel sig. high : idlsig.out
c       DONE: 4) writes out array for times when pairs of  particles jump: idljmp.out
c       DONE: 5) writes out array for jump times and particles involved: jump.out
c       to do these things have new particle color at each step array
c       summer 1995:
c       jj and bgl sig mods to handle 32 codes; eliminate 4d array by getting values up front
c       add "menu" so don't have to do all options.. debugging as of 15-aug-95..bl
c       8-feb-96 bgl version xslope32b - fix bugs in program for jump option
c       06-mar-96 bgl version xslope32f - correct POSSIBLE divide by 0 at definition of scn variable 
c       27-aug-96 bgl version xslope32g -  add "masking" to output to grvpat16 and later... pairs
c          without sig agg are not asessed by grvpat.
c       xslopebg.f = max 16 codes, but 10000 each and 1000 steps
c
c       xslopebgv2 bgl may1,2000: add routine direct3D as subroutine called with option 1
c       bl june 2009
c       xslope 09 add color for PDFT plot line when agg vel is sig relatve to max surrogate agg vel
c
c       in progress increase to up to 10000 control cycles
c
c       parameters in head_04.def....
C       ****************-DEFINITION OF HEADER BLOCK STRUCTURE-************
c       parameter (ncl=64)!number of codes handled
c       parameter (nsz=480000)!size of charge history array - 1st dimension
c       parameter (nsppc=10000)!max number of spikes in the busiest channel
c       parameter (lstsz=640000)!Product of nsz*nsppc
c

      EXTERNAL GETLOG_
c        ncell=ncl ! same as parameter ncl
      parameter (NCFR=10000)
      parameter (ncell=64)
      parameter (npts=1000)
      parameter (ipairs=2016)   ! =ncell*(ncell-1)/2
      parameter (sigcolor=2)    ! 2 = red
      parameter (maxcolor=13)    ! 13 = blue
c
c       NCFR =  number of control files
c       ...defines sig. level: 100 = .01 , 1000 = .001 etc.
c       ncell = max number of cells in data file
c       npts = max num time steps in orig. grav calc.
c       ipairs=number of unique pairs for ncell
c       jmpval=integer =number of concurrent "jumps" to be considered JUMP
c
c       min and max intensity levels in scaled plots..
      parameter (is_size=40,imnscl=0,imxscl=255)
c
c
c
      character*1 inst
      character*16 PLOT
      character*30    LABEL
      character*30    Strng
      character*60    title
c       Next line part of routine to let each user invoke diff window
c
      character*30 posfile
      character*30 outfile
      character*30 newfile
      character*5 c_i
      character*5 c_ii
      character*5 c_jj
      character*30 surfile 
      character*30 surfil (ipairs) ! array of NCFR file names
      character*30 confil (NCFR) ! array of NCFR file names
      DIMENSION  XXX(npts+1),YYY(npts+1)
      DIMENSION PDIS(ncell,ncell,npts)
      DIMENSION CDIS(ncell,ncell,npts)
c       new arrays:
c       sl contains "raw" slope (velocity) values for each
c       pair of particles following each step.
c       amxsl contains max. slope (velocity) values for each step
c            for each pair of CONTROL particles over all control samples
c       isig contains 1 at steps where each pair exhibits "sig" aggregation
c       according to criterion in use for that "run".
      DIMENSION sl(ncell,ncell,npts-1)
      DIMENSION vel(npts-1)
      DIMENSION amxsl(ncell,ncell,npts-1)
      DIMENSION isig(ncell,ncell,npts-1)
      DIMENSION icolr(npts)     ! color segments of PDFT - sig velocity
      dimension idisary(ipairs,npts-1)
      dimension imask (ipairs)
      dimension rmax(ipairs)
      dimension itmpary(ipairs,npts-1)
      dimension inwary(ipairs,npts-1)
      dimension aryidl(ipairs,npts-1)
      dimension yca(ncell,ncell)
c         dimension diff_arr(ncell,ncell,npts)
      dimension yyy_max(ncell,ncell,npts)
      dimension yyy_min(ncell,ncell,npts)
      dimension yyy_mean(ncell,ncell,npts)
      dimension icell1(ipairs),icell2(ipairs)
      dimension ijcell1(ipairs),ijcell2(ipairs)
      include 'config.defs'
      include 'sbparam.defs'
      real XCKB,XCORIG,XC,YC,xstp
      real top,bottom,left,right,hrange,vrange !parameters for graph out
      integer*4 Edge,EndOfLine,tick,idx,cnt(ncell)
      parameter (Edge=TRUE)     !used by "interior style"
      parameter (EndOfLine=FALSE) !used by "text2d"
      parameter (icnt=ncell*ncell*npts)
      data yyy_min / icnt*1000000000.0/ !initialize mins to max
      

c********************************************************************************
c       READ IN CONTROL *.pos FILES ...throwing away CBLOCK, IBLOCK,
c       and RBLOCK data plus TIMAT as loop through the control files... STORE
c        ONLY the data points in the 4 dimensional array CDIS.
c********************************************************************************
c
c       clear screen
 1    continue
c
c       choices for this run...
c       iline defined in gaus subroutine and passed to newplot2 to define line width
c       .. changed from default when gaus smooting selected - better image on screen than no smooth
c       .. nosmooth uses default line width.
C
C
c
c       **********************************************************************
c       options for color display - facilitate various hard copy possibilities
c       **********************************************************************
58540 print 58546
58546 format (2x,'1-color  2-gray scale')
      read (*,'(A)') inst
      if (inst.eq.'1')iopt=1
      if (inst.eq.'2')iopt=2
      if ((inst.ne.'1').and. (inst.ne.'2')) goto 58540
c       ********************************************************************
c       option for "bin" value in slope plots 
c       *******************************************************************
      inc = 1 
21854 print 28546
28546 format (2x,'Calc. slope for every Nth step, (def=1):',$)
      read (*,'(I3)') inc
      if (inc.gt.npts/2) goto 21854
      if (inc.eq.0) inc = 1
c       ***********************************************
c       ***********************************************
18540 print 18546
18546 format (2x,'1- single file  2- read controls')
      read (*,'(A)') inst
      if (inst.eq.'1') then
         icflag=1
         goto 18550
      end if
      if (inst.eq.'2') then
         icflag=2
         goto 98540
      end if
      goto 18540
c
c       **********************************************
c       How many cycles?
98540 print 98546
98546 format (2x,'1- 100 *.rdt files; 2- 1000, 3 = 10K')
      read (*,'(A)') inst
      if (inst.eq.'1') then
         ncfile=100
         goto 18548
      end if
      if (inst.eq.'2') then
         ncfile=1000
         goto 18548
      end if
      if (inst.eq.'3') then
         ncfile=10000
         goto 18548
      end if

      goto 98540
18548 print *,' Be patient. Reading surrogate/control files...'
      goto 18550
c
c
c
c
C       NEXT 4 lines create custom window
C
C        EXTERNAL GETLOG_
18550 continue
c      fildes=gopen(1600,1120,0,0,'xslope')

      if (icflag.eq.1) goto 3200

c
c       generate surrogate pos file names to read in
c
      do i = 1,NCFR
         write (c_i,'(I5)') i
         c_i=adjustl(c_i)
         newfile='sh'//trim(c_i)//'.pos'
         confil(i)=newfile
c       print '(''confil('',i5,'')='',a30)',i,confil(i)

      end do
c
c     Generate NAMES for surrogate velocity distribution files names for each pair
c  
c       print '(''ncell='',I10)',ncell
      isurplt=0
      do II=2,ncell
         do JJ=1,II-1
            isurplt=isurplt+1
            write (c_ii,'(I5)') II
            write (c_jj,'(I5)') JJ
            c_ii=adjustl(c_ii)
            c_jj=adjustl(c_jj)
            surfile='sur'//trim(c_ii)//'_'//trim(c_jj)//'.sur'
            surfil(isurplt)=surfile
c       print '(''surfil('',i5,'')='',a30)',isurplt,surfil(isurplt)
         end do
      end do
c
c       Initialize some arrays
c
c
      do II=2,ncell
         do JJ=1,II-1
            yca(II,JJ)=0.0
            do K=1,NMPTS
               yyy_max(II,JJ,K)= 0.0
               yyy_min(II,JJ,K)= 1000000000. !needs to be much bigger than 100
               yyy_mean(II,JJ,K)= 0.0
               amxsl(II,JJ,K)= 0.0
            end do
         end do
      end do
c       *******************************************************************
c       get max min and mean of sur controls at each step-for sig 'band" in PDFT plots
c       ******************************************************************
      biggest=0.0
      
      amxsl = 0
      do jz=1,ncfile
         print 2001, jz
 2001    format (2x, 'Reading in control file:',1x,I5)
         posfile= confil(jz)
         OPEN (UNIT=1,file=posfile,status='OLD',form='UNFORMATTED')
         read (1) CBLOCK,IBLOCK,RBLOCK !read in header info for this plot
         read(1) TIMAT          !read in spike times
         NMPTS=LFRAM            !# of points in file for each pair of cells 
         DO K=1,NMPTS
            READ (1,end=3130,iostat=ios) 
     +           ((CDIS(II,JJ,K),JJ=1,II-1),II=2,N)
         end do
 3130    close (1)
         do II=2,N
            do JJ=1,II-1
               DO K=1,NMPTS
                  yca(II,JJ)=AMAX1(yca(II,JJ),CDIS(II,JJ,K))
                  yyy_max(II,JJ,K)=AMAX1(yyy_max(II,JJ,K),CDIS(II,JJ,K))
                  biggest=AMAX1(biggest,yyy_max(II,JJ,K))
                  yyy_min(II,JJ,K)=AMIN1(yyy_min(II,JJ,K),CDIS(II,JJ,K))
                  yyy_mean(II,JJ,K)=yyy_mean(II,JJ,K) + 
     +                 CDIS(II,JJ,K)
               end do
            end do
         end do







C       *************************************************
c       Fill amxsl array here
c     It contains the largest slope (out of ncfile control sets)
c     for each pair at each step.
C     IF inc > 1, then inc successive bins have same value.
c     This allows slope over several steps to be calculated..
C       *************************************************
c       
c
         do II=2,N
            do JJ=1,II-1
               DO K=1,NMPTS-inc,inc
                  amxsl(II,JJ,K)=AMAX1(amxsl(II,JJ,K),
     +                 CDIS(II,JJ,K)-CDIS(II,JJ,K+inc))
                  tmpsl= amxsl(II,JJ,K)
                  do IZ=K,K+inc-1
                     if (IZ.eq. npts-1)goto 3007
                     amxsl(II,JJ,IZ)=tmpsl
                  end do
 3007             continue
               end do
            end do
         end do
      end do

c
c

      print 6
 6    format (9x,'1 - min, max & mean control plot (default)',/,
     +     9x,'2 - min only')
      read (*,'(A)') inst
      if (inst.eq.' ') iflag = 1
      if (inst.eq.' ') PLOT = 'min, max, & mean'
      if (inst.eq.'1') iflag = 1
      if (inst.eq.'1') PLOT = 'min, max, & mean'
      if (inst.eq.'2') iflag = 2
      if (inst.eq.'2') PLOT = 'min only'
c
c
c
c       ********************************************************
c
c
c       yyy arrays filled - now calculate correct mean values for yyy_mean
      do II=2,N
         do JJ=1,II-1
            DO K=1,NMPTS
               yyy_mean(II,JJ,K)=yyy_mean(II,JJ,K)/ncfile 
            end do
         end do
      end do

c       ************************************************************
c       *********************************************************

 3200 continue
      print 5
 5    format (2x,'..input file? :',$)
      read (*,'(A)') posfile
      if (posfile.eq.'')posfile = '2006-08-09P_pm10.pos'
c
      XCKB=0.0                  !init. var
      OPEN (UNIT=1,file=posfile,status='OLD',form='UNFORMATTED')
      read (1) CBLOCK,IBLOCK,RBLOCK
      read(1) TIMAT             !read in spike times

C       **************** READ DATA ***************************
      NMPTS=LFRAM               !# of points in file for each pair of cells 
      DO 120 K=1,NMPTS
         READ (1,END=130) ((PDIS(II,JJ,K),JJ=1,II-1),II=2,N)
 120  CONTINUE
 130  CLOSE (1)
C
C
      XCORIG= (ENDTIME+DELTIM)/1000. !final state time in sec. of data file
      call mode(4)
      PRINT 22,XCORIG
 22   FORMAT (2X,'default timespan will be ',F5.1,' sec. or ?:',$)
      READ (5,18) XCKB
 18   format (F5.1)
      if (XCKB.eq.0.0)then
         XC=XCORIG
      else 
         XC=XCKB
      end if
      if (XC.lt.XCORIG) NMPTS= int((LFRAM*(XC/XCORIG)))
c
c
c        call GETLOG_(USER)
c        WINDOW1='xslope_04'
      fildes=gopen(1600,1120,0,0,'xslope')
      if (fildes.eq.-1)stop
      call colors(fildes)

c**************************************************************
c       get number of lines (unique pairs) in current display
c**************************************************************
      number=0
      numlin=N-1
78888 number=numlin+number
      numlin=numlin-1
      if (numlin.eq.0) goto 78889
      goto 78888
78889 continue
      iyrange = number
c       NMPTS - 1 is number of points plotted along x- axis 
      ixrange = NMPTS - 1

c**************************************************************
c       fill sl array here
c       IF inc > 1 then inc successive bins have same value..
c       allows slope over several steps to be calculated..
c       DEFAULT INC MUST BE SET TO 1
c**************************************************************
      sl = 0
      do II=2,N
         do JJ=1,II-1
            DO 207 K=1,NMPTS-inc,inc
               cursl = PDIS(II,JJ,K)-PDIS(II,JJ,K+inc)
               do 208 IZ=K,K+inc-1
                  if (IZ.eq. npts-1)goto 207
                  sl(II,JJ,IZ) = cursl
 208           continue
 207        continue
         end do
      end do


c********************************************
c       fill isig array (has 1 less element than NMPTS)
c********************************************
      do II=2,N
         do JJ=1,II-1
            do 3810 K=1,NMPTS
               if (K.eq.1) goto 3810

               isig(II,JJ,K-1) = imnscl ! coolest color
               if (PDIS(II,JJ,K).lt.yyy_min(II,JJ,K)) then
                  isig(II,JJ,K-1) = imxscl ! hottest color
               end if
 3810       continue
         end do
      end do

 508  format (2x,'Must have control files for this option '
     +     /'         (<CR> to continue): ',$)
C*******************************************************************
C       Main Menu
C*******************************************************************
 506  call mode(4)
      print *,'CHOOSE DATA DISPLAY:'
      print *,''
      print *,'         1. Pairwise aggregation'
      print *,'         2. Times of Significant aggregation'
      print *,'           . also makes MASK array for grvpat'
      print *,'         3. All movements -> reduced distance'
      print *,'              - writes: idlmov.out to disk'
      print *,'         4. LOOP TO START of  program'
      print *,'         5. EXIT program'
      read (*,'(A)') inst
      if (inst .eq. '1') goto 501
      if (inst .eq. '2') goto 507
      if (inst .eq. '3') goto 43485
      if (inst .eq. '4') goto 1
      if (inst .ne. '5') goto 506
      call gclose(fildes)
      call exit

c*************************************************************************
c       MAIN PROGRAM LOOP : TIME vs. DISTANCE PLOT FOR EACH PAIR
c*************************************************************************
 501  call colors(fildes)
      call background_color(fildes,R(2),G(2),B(2))
      call line_color(fildes,R(1),G(1),B(1))
      call text_color(fildes,R(1),G(1),B(1))
      call perimeter_color(fildes,R(1),G(1),B(1))
16548 format (2x,'Select particle to display ',
     +     '(enter 0 for all): ',$)
16549 format (I4)
      print 16548
      read (*,16549) ij1
      isurplt=0
      do II=2,N
         do JJ=1,II-1
            isurplt=isurplt+1

            if (ij1 .ne. 0) then
               if (ij1 .ne. II .and. ij1 .ne. JJ) goto 3490
            end if

            call clear_control (fildes,CLEAR_DISPLAY_SURFACE)
            call clear_view_surface (fildes)

c        auto scaling of plot AND definition of global scaling for grvpat
            
            YC=0.0
            DO 505 K=1,NMPTS
               YC= AMAX1(YC,PDIS(II,JJ,K))
 505        end do

            if (icflag .eq. 1 .or. iflag .eq.2) THEN
               YC=YC+0.05*YC
            else
               YC = AMAX1(YC,yca(II,JJ))
               YC=YC+0.05*YC
            end if

C
C       graph output parameters. These parameters define where on the 
C       screen the graphic output will appear.  All have a range of
C       0.0 to 1.0.
C                   Right - the amount of space to the right of graph
C                   Left  - the amount ofspace to the left of graph
C                   Top   - the amount of space above the graph 
C                   Botttom - the amount of space below the graph
C
            right=0.13
            left=0.15
            top=0.07
            bottom=0.35
            hrange=1-left-right
            vrange=1-top-bottom
C
            call vdc_extent(fildes,0.0,0.0,0.0,1.25,1.0,0.0)
            call view_port(fildes,left*1.25,bottom,(1-right)*1.25,
     +           (1-top))
            call view_window(fildes,0.0,0.0,XC,YC)
            call interior_style(fildes,INT_HOLLOW,Edge) !hollow rectangles
            call rectangle (fildes,0.0,0.0,XC,YC) !frame the viewport

            call make_picture_current(fildes) !dump the buffer
C
C
C       plot the data
c       call mode(4)
c       print 170
c170    format(2x,'graph title (max. 26 char.)? :',$)
c       read (*,'(A)') LABEL
            LABEL = posfile
C
c       if ((II.eq.N).and.(JJ.eq.II-1)) then
c       print *,'got to pt 1b'
c       pause
c       end if
C       PLOT TIME vs. DISTANCE LINE
C
c       lincol=6  !rich purple
            lincol=0            !black
            m = lincol + 1
            call line_color(fildes,R(m),G(m),B(m)) 
            xstp=XC/NMPTS
c
            DO 205 K=1,NMPTS    ! successive 'inc' steps same color 
               YYY(K) = PDIS(II,JJ,K)
               XXX(K)= K*xstp
               if (icflag.eq.2.and.sl(II,JJ,K).gt.amxsl(II,JJ,K)) then
                  icolr(K)=int(sigcolor)
               else 
                  icolr(K)=1
               end if
 205        end do
c
c Make this line thicker.. because of coordinates looks distorted
c verticle part looks thicker
c
            call line_width (fildes,0.002,VDC_UNITS) ! def = 0.0
            call line_endpoint (fildes,rounded) 
            call move2d(fildes,0.0,100.)
            do 206 jz=1,NMPTS
               if (jz > 1.and.icolr(jz-1).eq.sigcolor) then
                  m = icolr(jz-1) + 1
                  call line_color(fildes,R(m),G(m),B(m))
               else
                  m = lincol + 1
                  call line_color(fildes,R(m),G(m),B(m)) 
               end if
               call draw2d(fildes,XXX(jz),YYY(jz))
 206        end do

c       PLOT MAX SURROGATE VELOCITIES WHEN VEL is SIG
c       Note - successive inc steps same color (= 1 eval step) 
c$$$            m = maxcolor + 1
c$$$            call line_color(fildes,R(m),G(m),B(m))
c$$$            ypmax=0.0
c$$$            ymax=2.0
c$$$            do  jz=1,NMPTS-1
c$$$               if (sl(II,JJ,jz).gt.amxsl(II,JJ,jz)) then
c$$$                  ypmax=AMAX1(ypmax,amxsl(II,JJ,jz))
c$$$                  call move2d(fildes,XXX(jz),ymax)
c$$$                  call draw2d(fildes,XXX(jz),ymax+amxsl(II,JJ,jz))
c$$$                  print *,'max: ', XXX(jz), amxsl(II,JJ,jz)
c$$$               end if
c$$$            end do
c
c       ABOVE IT PLOT PAIR VELOCITY AT SAME STEPS - same color as PDFT plot
c       Note - successive inc steps same color (= 1 eval step) 

            if (icflag.eq.2) then
               m = int(sigcolor + 1)
               call line_color(fildes,R(m),G(m),B(m))
               yval=ypmax+4.0
               do  jz=1,NMPTS-1
                  if (sl(II,JJ,jz).gt.amxsl(II,JJ,jz)) then
                     call move2d(fildes,XXX(jz),yval)
                     if (amxsl(II,JJ,jz).eq.0) then
                        call draw2d(fildes,XXX(jz), yval+10)
                     else
                        call draw2d(fildes,XXX(jz), yval+
     +                       sl(II,JJ,jz)/amxsl(II,JJ,jz)*10)
                     end if
                     print *,'data: ',XXX(jz), sl(II,JJ,jz)
                  end if
               end do
            end if
c Line width and color back to default
            call line_width (fildes,0.0,VDC_UNITS) ! def = 0.0
            m = lincol + 1
            call line_color(fildes,R(m),G(m),B(m)) 






c       *************************************************************
c       *************************************************************
c               A SLOPE CALC.
c       *********************************************************
c
c        if (icflag.eq.1) goto 52000 !no control
            goto 52005          ! skip line thickening routine

c       MAKE LINE THICKER..
c
            DO 4205 K=1,NMPTS
               YYY(K) = PDIS(II,JJ,K)-(0.0025*YC)
 4205       end do
            call move2d(fildes,0.0,100.)
            do 4206 jz=1,NMPTS
               call draw2d(fildes,XXX(jz),YYY(jz))
 4206       end do
c
c       thicker still for monochrome...
c
            DO 4305 K=1,NMPTS
               YYY(K) = PDIS(II,JJ,K)-(0.005*YC)
 4305       end do
            call move2d(fildes,0.0,100.)
            do 4306 jz=1,NMPTS
               call draw2d(fildes,XXX(jz),YYY(jz))
 4306       end do

            DO 4405 K=1,NMPTS
               YYY(K) = PDIS(II,JJ,K)-(0.0075*YC)
 4405       end do
            call move2d(fildes,0.0,100.)
            do 4406 jz=1,NMPTS
               call draw2d(fildes,XXX(jz),YYY(jz))
 4406       end do
c
c
c       *************************************************************
c       **************************************************************
c       PLOT SIGNIFICANCE BAND
c
C       *************************************************
c       ...first plot max point in all control files
c       at each step..
C       *************************************************
c       lincol=9 ! medium green
52005       lincol=0            ! black
            m = lincol + 1
            call line_color(fildes,R(m),G(m),B(m))
            call line_type (fildes,DASH)
c
c       min only so skip this part
c
            if (iflag.eq.2) goto 6666

            call move2d(fildes,0.0,100.)

            do 3706 jz=1,NMPTS
               call draw2d(fildes,XXX(jz),yyy_max(II,JJ,jz))
 3706       end do

c********************************************************
c       plot min points in all control files at each step
c********************************************************
 6666       call move2d(fildes,0.0,100.)
            do 3806 jz=1,NMPTS
               call draw2d(fildes,XXX(jz),yyy_min(II,JJ,jz))
 3806       end do
c
c
            if (iflag.eq.2) goto 52000
c
c********************************************
c       PLOT MEAN OF CONTROL POINTS AT EACH STEP
c********************************************
c       lincol=8 ! 
            lincol=0            ! 
            m = lincol + 1
            call line_color(fildes,R(m),G(m),B(m))
            call line_type (fildes,DOT)

            call move2d(fildes,0.0,100.)
            do 3906 jz=1,NMPTS
               call draw2d(fildes,XXX(jz),yyy_mean(II,JJ,jz))
 3906       end do
c
            call make_picture_current(fildes) !dump the buffer
c
c
c       *************************************************************
c       *************************************************************
c               A SLOPE CALC.
c       *********************************************************
c       plot x-tics
52000       call line_color(fildes,R(1),G(1),B(1))
            call line_type(fildes,SOLID)
            xmrk=0.0
            do while (xmrk.le.XC)
               call move2d(fildes,xmrk,0.0)
               call draw2d(fildes,xmrk,(.02*YC))
               xmrk=xmrk+(0.25*XC)
            end do
c       plot y -tics
            ymrk=0.0
            do while (ymrk.le.YC)
               call move2d(fildes,0.0,ymrk)
               call draw2d(fildes,(.02*XC),ymrk)
               ymrk=ymrk+(0.25*YC)
            end do
C
!-----  do main title
            call text_font_index(fildes, 4) !sans serif font
            call text_alignment(fildes,TA_LEFT,TA_BOTTOM,0.0,0.0)
            call character_height(fildes,0.05) !9% of current 0-1 Y range
c       write(unit=Strng,fmt='(a)')LABEL
c       call text2d(fildes,0.50,0.94,Strng,VDC_TEXT,
c     + EndOfLine)
            call text2d(fildes,0.5,0.94,'Particles:'//char(0),
     +           VDC_TEXT,1)
            write(unit=strng,fmt='(2x,i2,x,a1)')II,char(0)
            call append_text(fildes,strng,VDC_TEXT,1)
            write(unit=strng,fmt='(2x,i2,x,a1)')JJ,char(0)
            call append_text(fildes,strng,VDC_TEXT,0)
c
c
c
!----   do x-axis title
            call text_font_index(fildes,4)
            call text_alignment(fildes,TA_CENTER,TA_BOTTOM,0.0,0.0)
            call character_height(fildes,.04)
            call text2d(fildes,1.25*(hrange/2.+left),bottom-.14*vrange,
     +           'Time (sec.)'//char(0),VDC_TEXT, EndOfLine)
!----   do y-axis title
            call text_alignment(fildes,TA_CENTER,TA_TOP,0.0,0.0)
            call text_orientation2d(fildes,-1.0,0.0,0.0,1.0)
            call text2d(fildes,1.25*(left-.17*hrange),vrange/2.+bottom,
     +           'Distance'//char(0),VDC_TEXT, EndOfLine)
!---- label x tics
            size=0.045
            call clip_indicator(fildes,CLIP_TO_VDC)
            call text_orientation2d(fildes,0.0,1.0,1.0,0.0)
            call character_height(fildes,size)
            call text_alignment(fildes,TA_CENTER,TA_TOP,0.0,0.0)
            X=0.0
            do 221 ival=0,4
               write(unit=Strng,fmt='(F5.1,x,a1)')X,char(0)
               xaddr=(.25*ival*hrange+left)*1.25
               call text2d(fildes,xaddr,bottom-size/5.,Strng,VDC_TEXT,
     +              EndOfLine)
               X=X+(0.25*XC)
 221        end do
!-----  label y tics
            call text_alignment(fildes,TA_RIGHT,TA_HALF,0.0,0.0)
            Y=0.0
            do 222 jval=0,4
               write(unit=Strng,fmt='(f5.1,x,a1)')Y,char(0)
               yaddr=.25*jval*vrange+bottom
               call text2d(fildes,left*1.25,yaddr,Strng,VDC_TEXT,
     +              EndOfLine)
               Y=Y+(.25*YC)
 222        end do
C
            call make_picture_current(fildes) !dump the buffer
c 224   call mode(4)

            call text_alignment(fildes,TA_LEFT,TA_BOTTOM,0.,0.)
            call text_path(fildes,PATH_RIGHT)
            call text_line_path(fildes,PATH_DOWN)
            call text_orientation2d(fildes,0.0,1.0,1.0,0.0)
c        call character_height(fildes,0.020)
            call character_height(fildes,0.012)
            call text_font_index(fildes,1)
c
            call vdc_extent(fildes,0.,0.,0.,1.25,1.0,0.)
            call view_window(fildes,0.0,0.0,20.0,20.0)
            call view_port(fildes,.12*1.25,0.,.88*1.25,.2)
c
c
c       spike counter for N neurons
c
            do 3 idx=1,N
               tick=1
               do while((TIMAT(tick,idx).gt.0.0).and.
     +              (tick.le.nsppc))
                  tick=tick+1
               end do
               cnt(idx)=tick-1
 3          end do
c
c
c       particle count block output
c
c
            jtmp = min(N,18)
            call text2d(fildes,.05,.250,'particle #    : '//char(0),
     +           VDC_TEXT,1)

            if (jtmp .gt. 9) then
               do i=1,9
                  write(unit=strng,fmt='(4x,i3,x,a1)')i,char(0)
                  call append_text(fildes,strng,VDC_TEXT,1)
               end do
               do i=10,jtmp
                  write(unit=strng,fmt='(4x,i3,x,a1)')i,char(0)
                  call append_text(fildes,strng,VDC_TEXT,1)
               end do
            else
               do i=1,jtmp
                  write(unit=strng,fmt='(4x,i3,x,a1)')i,char(0)
                  call append_text(fildes,strng,VDC_TEXT,1)
               end do
            end if
            call append_text(fildes,'',VDC_TEXT,0)

            call text2d(fildes,.05,.230,'event code    : '//char(0),
     +           VDC_TEXT,1)
            do j=1,jtmp
               write(unit=strng,fmt='(4x,i3,x,a1)')IDSET(j),char(0)
               call append_text(fildes,strng,VDC_TEXT,1)
            end do
            call append_text(fildes,'',VDC_TEXT,0)
c
            call text2d(fildes,.05,.210,'# of events   : '//char(0),
     +           VDC_TEXT,1)
            do k=1,jtmp
               write(unit=strng,fmt='(3x,i4,x,a1)')cnt(k),char(0)
               call append_text(fildes,strng,VDC_TEXT,1)
            end do
            call append_text(fildes,'',VDC_TEXT,0)
c
            if (N .gt. 18) then
               jtmp = min(N,36)
               call text2d(fildes,.05,.190,'particle #    : '//char(0),
     +              VDC_TEXT,1)
               do i=19,jtmp
                  write(unit=strng,fmt='(5x,i2,x,a1)')i,char(0)
                  call append_text(fildes,strng,VDC_TEXT,1)
               end do
               call append_text(fildes,'',VDC_TEXT,0)

               call text2d(fildes,.05,.170,'event code    : '//char(0),
     +              VDC_TEXT,1)
               do  j=19,jtmp
                  write(unit=strng,fmt='(4x,i3,x,a1)')IDSET(j),char(0)
                  call append_text(fildes,strng,VDC_TEXT,1)
               end do
               call append_text(fildes,'',VDC_TEXT,0)

               call text2d(fildes,.05,.150,'# of events   : '//char(0),
     +              VDC_TEXT,1)
               do j=19,jtmp
                  write(unit=strng,fmt='(3x,i4,x,a1)')cnt(j),char(0)
                  call append_text(fildes,strng,VDC_TEXT,1)
               end do
               call append_text(fildes,'',VDC_TEXT,0)
            end if

            if (N .gt. 36) then
               jtmp = min(N,54)
               call text2d(fildes,.05,.130,'particle #    : '//char(0),
     +              VDC_TEXT,1)
               do i=37,jtmp
                  write(unit=strng,fmt='(5x,i2,x,a1)')i,char(0)
                  call append_text(fildes,strng,VDC_TEXT,1)
               end do
               call append_text(fildes,'',VDC_TEXT,0)

               call text2d(fildes,.05,.110,'event code    : '//char(0),
     +              VDC_TEXT,1)
               do j=37,jtmp
                  write(unit=strng,fmt='(4x,i3,x,a1)')IDSET(j),char(0)
                  call append_text(fildes,strng,VDC_TEXT,1)
               end do
               call append_text(fildes,'',VDC_TEXT,0)

               call text2d(fildes,.05,.090,'# of events   : '//char(0),
     +              VDC_TEXT,1)
               do j=37,jtmp
                  write(unit=strng,fmt='(3x,i4,x,a1)')cnt(j),char(0)
                  call append_text(fildes,strng,VDC_TEXT,1)
               end do
               call append_text(fildes,'',VDC_TEXT,0)
            end if

            if (N .gt. 54) then
               jtmp = min(N,72)
               call text2d(fildes,.05,.070,'particle #    : '//char(0),
     +              VDC_TEXT,1)
               do i=55,jtmp
                  write(unit=strng,fmt='(5x,i2,x,a1)')i,char(0)
                  call append_text(fildes,strng,VDC_TEXT,1)
               end do
               call append_text(fildes,'',VDC_TEXT,0)

               call text2d(fildes,.05,.050,'event code    : '//char(0),
     +              VDC_TEXT,1)
               do j=55,jtmp
                  write(unit=strng,fmt='(4x,i3,x,a1)')IDSET(j),char(0)
                  call append_text(fildes,strng,VDC_TEXT,1)
               end do
               call append_text(fildes,'',VDC_TEXT,0)

               call text2d(fildes,.05,.030,'# of events   : '//char(0),
     +              VDC_TEXT,1)
               do j=55,jtmp
                  write(unit=strng,fmt='(3x,i4,x,a1)')cnt(j),char(0)
                  call append_text(fildes,strng,VDC_TEXT,1)
               end do
               call append_text(fildes,'',VDC_TEXT,0)
            end if

c
c
c
c
c
c
c
c
c
c********************************************************
C       print settings
c********************************************************
            call text2d(fildes,1.092,.975,'PRG:'//char(0),
     +           VDC_TEXT,1)
            write(unit=strng,fmt='(a)')PROGID
            call append_text(fildes,strng,VDC_TEXT,0)
c
            call text2d(fildes,1.092,.950,'M-C:'
     +           //char(0),VDC_TEXT,1)
            write(unit=strng,fmt='(a)')PLOT
            call append_text(fildes,strng,VDC_TEXT,0)
c
            call text2d(fildes,1.092,.925,'FL:'
     +           //char(0),VDC_TEXT,1)
            call character_height(fildes,0.017)
            write(unit=strng,fmt='(a)')posfile
            call append_text(fildes,strng,VDC_TEXT,0)
            call character_height(fildes,0.020)
c
            call text2d(fildes,1.092,.900,'INT:'
     +           //char(0),VDC_TEXT,1)
            write(unit=strng,fmt='(i3,x,a6,a1)')IHS,"spikes",char(0)
            call append_text(fildes,strng,VDC_TEXT,0)
c
            call text2d(fildes,1.092,.875,'END:'
     +           //char(0),VDC_TEXT,1)
            write(unit=strng,fmt='(f10.1,x,a1)')ENDTIME,char(0)
            call append_text(fildes,strng,VDC_TEXT,0)
c
            call text2d(fildes,1.092,.850,'FSGN:'
     +           //char(0),VDC_TEXT,1)
            write(unit=strng,fmt='(f4.1,x,a1)')FS,char(0)
            call append_text(fildes,strng,VDC_TEXT,1)
c
            if(FS.eq.1.0) then
               call append_text(fildes,'excit'
     +              //char(0),VDC_TEXT,0)
            else
               call append_text(fildes,'inhib'
     +              //char(0),VDC_TEXT,0)
            end if
c
            if (ISHIFTA.EQ.0) THEN
               call text2d(fildes,1.092,.825,'ACC DECAY FWD'
     +              //char(0),VDC_TEXT,0)
            else
               call text2d(fildes,1.092,.825,'ACC DECAY BKWD'
     +              //char(0),VDC_TEXT,0)
            end if
c
            call text2d(fildes,1.092,.800,'NRM:'
     +           //char(0),VDC_TEXT,1)
            write(unit=strng,fmt='(i3,x,a1)')NORM,char(0)
            call append_text(fildes,strng,VDC_TEXT,0)
c
            if (ISHIFTE.EQ.0) THEN
               call text2d(fildes,1.092,.775,'EFF DECAY FWD'
     +              //char(0),VDC_TEXT,1)
            else
               call text2d(fildes,1.092,.775,'EFF DECAY BKWD'
     +              //char(0),VDC_TEXT,1)
            end if
c
            call text2d(fildes,1.092,.750,'NRM:'
     +           //char(0),VDC_TEXT,1)
            write(unit=strng,fmt='(i3,x,a1)')NORM,char(0)
            call append_text(fildes,strng,VDC_TEXT,0)
c
            call text2d(fildes,1.092,.725,'TM ST:'
     +           //char(0),VDC_TEXT,1)
            write(unit=strng,fmt='(f6.1,x,a1)')DELTIM,char(0)
            call append_text(fildes,strng,VDC_TEXT,0)
C
            call text2d(fildes,1.092,.700,'SLD:'
     +           //char(0),VDC_TEXT,1)
            write(unit=strng,fmt='(E8.2,x,a1)')SLIDE,char(0)
            call append_text(fildes,strng,VDC_TEXT,0)
c
            call text2d(fildes,1.092,.675,'FINC:'
     +           //char(0),VDC_TEXT,1)
            write(unit=strng,fmt='(f6.1,x,a1)')FWDINC,char(0)
            call append_text(fildes,strng,VDC_TEXT,0)
c
            call text2d(fildes,1.092,.650,'BINC:'
     +           //char(0),VDC_TEXT,1)
            write(unit=strng,fmt='(f6.1,x,a1)')BAKINC,char(0)
            call append_text(fildes,strng,VDC_TEXT,0)
C
            call text2d(fildes,1.092,.625,'FTAU:'
     +           //char(0),VDC_TEXT,1)
            write(unit=strng,fmt='(f6.1,x,a1)')FWDTAU,char(0)
            call append_text(fildes,strng,VDC_TEXT,0)
c
            call text2d(fildes,1.092,.600,'BTAU:'
     +           //char(0),VDC_TEXT,1)
            write(unit=strng,fmt='(f6.1,x,a1)')BAKTAU,char(0)
            call append_text(fildes,strng,VDC_TEXT,0)
c
            call text2d(fildes,1.092,.575,'CRDI:'
     +           //char(0),VDC_TEXT,1)
            write(unit=strng,fmt='(f6.1,x,a1)')CRDI,char(0)
            call append_text(fildes,strng,VDC_TEXT,0)
c
            call text2d(fildes,1.092,.550,'MS BT:'
     +           //char(0),VDC_TEXT,1)
            write(unit=strng,fmt='(f8.1,x,a1)')stimper,char(0)
            call append_text(fildes,strng,VDC_TEXT,0)
c
            call text2d(fildes,1.092,.525,'ST SVD:'
     +           //char(0),VDC_TEXT,1)
            write(unit=strng,fmt='(i4,x,a1)')ISTPFR,char(0)
            call append_text(fildes,strng,VDC_TEXT,0)
c
            call text2d(fildes,1.092,.500,'FR SVD:'
     +           //char(0),VDC_TEXT,1)
            write(unit=strng,fmt='(i5,x,a1)')LFRAM,char(0)
            call append_text(fildes,strng,VDC_TEXT,0)
c
            call make_picture_current(fildes) !dump the buffer
c
c
c       call mode(4)

c
c       Write surrogate file option..  
c
            print 1227
 1227       format (5x,'<cr> continue or w for sur file for plot:',$)
            read (*,'(A)') inst
            if((inst.eq.'w').or.(inst.eq.'W')) goto 1030
            goto 3490           ! Done with this pair
c Write distributions for sig aggs in this pair's PDFT plot
c
c  

 1030       continue
            if (compiler.eq.'gfortran') then
               OPEN (UNIT=8,file=surfil(isurplt),status='UNKNOWN',
     +              form='UNFORMATTED',ACCESS='STREAM')
            else
               OPEN (UNIT=8,file=surfil(isurplt),status='UNKNOWN',
     +              form='UNFORMATTED',ACCESS='DIRECT',RECL=4)
            end if
            write (8) II        ! first cell of pair
            write (8) JJ        ! second cell of pair
            write (8) npts      ! max number of plotted steps
            write (8) inc       ! num steps used to calc slope (=K-(K+inc))
            write (8) ENDTIME   ! time span of total plot in ms.
            write (8) ncfile    !write number of surrogate trials
            icount = 0
            do jz=1,NMPTS-inc, inc
               vel(jz) = sl(II,JJ,jz)
               icount = icount + 1
            end do
            write (8) icount
            write (8) vel
            do kz=1,ncfile
               posfile= confil(kz)
               OPEN (UNIT=1,file=posfile,status='OLD',
     +              form='UNFORMATTED')
               read (1) CBLOCK,IBLOCK,RBLOCK !read in header
               read(1) TIMAT    !read in spike time
               NMPTS=LFRAM      !# of points in file for each pair of cells 
               DO mz=1,NMPTS
                  READ (1,end=31301,iostat=ios) 
     +                 ((CDIS(IB,JB,mz),JB=1,IB-1),IB=2,N)
               end do
31301          close (1)
               do jz=1,NMPTS-inc, inc
                  vel(jz) = CDIS(II,JJ,jz) - CDIS(II,JJ,jz+inc)
               end do
               if (compiler.eq.'gfortran') then
                  write (8) vel
               else
                  do jz=1,NMPTS-inc, inc
                     write (8) vel(jz)
                  end do
               end if
            end do
            close (8)           ! done 


 3490       continue
         end do
      end do
      go to 506
c
c
c       OLD
c        fildes=gopen ('/dev/screen/xslopebg'//char(0),
c     +  OUTDEV,'sox11'// char(0),INIT)
c
c if necseeary, now use...
c
c       fildes=gopen ('/dev/screen/'//USER//'/'
c    + //WINDOW1//char(0),OUTDEV,'sox11'// char(0),INIT)

c       ******************************************************************
c
c               SLOPE DISPLAYS
c       **************************************************************
c
c       sequence of displays showing sig aggregation across all pairs
c       and velocity of aggregation
c
c
 507  continue
      if (icflag.eq.1) then
         print 508
         read (*,'(A)') inst
         goto 506               !no control
      end if

c****************************************************
c       PLOT: TIMES OF SIG. AGGREGATION ACROSS PAIRS
c****************************************************
      title='Significant aggregation'
      loop=0                    !no loop through rows with reordering by max
      iline=1                   ! use fat lines for this plot

c       ************************************************************
c       match pairs with idisary variables
c       ************************************************************
      j=0
      do II=2,N
         do JJ=1,II-1
            j=j+1
            ijcell1(j)=II
            ijcell2(j)=JJ
         end do
      end do

      ips=0
      do 13500 II=2,N
         do 13490 JJ=1,II-1
            ips=ips+1
            do 13480 K=1,ixrange
               itmpary(ips,K)=isig(II,JJ,K)
13480       continue
13490    continue
13500 continue

      iq = 0
      do while (iq .lt. ips)
c          clear array for histogram plots      
         do 4000 iz = 1,ipairs
            do 4010 iy= 1, npts-1
               inwary(iz,iy)=0
 4010       continue
 4000    continue

c          set up idisary and icell1/2 top 40
         do iy=1,is_size
            do K=1,ixrange
               idisary(iy,K) = itmpary(iy+iq,K)
            end do
            icell1(iy) = ijcell1(iy+iq)
            icell2(iy) = ijcell2(iy+iq)
         end do
c       sl_max - max value for color key not used in this call 
         if (iyrange - iq .lt. is_size) then !eliminates garbage on LAST screen
            call plot_op2(fildes,idisary,ixrange,iyrange,
     +           title,N,loop,iopt,iline,inwary,xc,'n',
     +           icell1,icell2,iyrange-iq)
         else
            call plot_op2(fildes,idisary,ixrange,iyrange,title,
     +           N,loop,iopt,iline,inwary,xc,'n',
     +           icell1,icell2,iyrange)
         endif

         iq = iq+is_size

         print 1346
         read (*,'(A)') inst
      end do
 1346 format (2x,'<CR> to continue...')

c       *****************************************************************
c       *****************************************************************
c       insert routine here to differentiate  pairs (row) that never showed sig. aggregation
c       in PDFT plot. This arry is written out in idl*.out headers
c       for later use by grvpat's "MASK" option..
c       NOTE: if mask option 'on', then:1=MASK LATER...... 0= NO MASK LATER
c       
      do nz=1,ipairs
         imask(nz)=1            ! default = mask this pair-DO NOT USE
      end do
c
      iz=0
      do 86403 II=2,N
         do 86423 JJ=1,II-1
            iz=iz+1
            isigflg=0
            do 77983 K=1,ixrange
               if (isig(II,JJ,K).eq.imxscl) isigflg= 1
77983       continue
c       next line: use this imask pair in grvpat if MASK on
            if (isigflg.eq.1) imask(iz)=0 !CLEAR MASK & USE IN GRVPAT*
86423    continue
86403 continue
c
c       end mod
      goto 506
c       
c
c
c       ********************************************************
c       PLOT: "RAW" neg. slopes ACROSS PAIRS - global normalization
c       ********************************************************
43485 iz=0
      do 43500 II=2,N
         do 43490 JJ=1,II-1
            iz=iz+1
            rmax(iz)=0.0
            do 43480 K=1,ixrange
               rmax(iz)=AMAX1(rmax(iz),sl(II,JJ,K))
43480       continue
43490    continue
43500 continue
c
c       get global max slope over all pairs
c
      rmaxval=0.0
      do 43470 iz=1,iyrange
         rmaxval=AMAX1(rmaxval,rmax(iz))
43470 continue
c       clear array for histogram plots 
      do 6000 iz = 1,ipairs
         do 6010 iy= 1, npts-1
            inwary(iz,iy)=0
 6010    continue
 6000 continue
c
c       fill main display and histogram arrays
c
      sl_max = 0.0
c
      iz=0
      if (rmaxval.ne.0.0)scn=(imxscl + 1)/rmaxval
      if (rmaxval.eq.0.0)scn=0.0
      do 44500 II=2,N
         do 44490 JJ=1,II-1
            iz=iz+1
            do 44480 K=1,ixrange
               sltmp=sl(II,JJ,K)
               if (sltmp.lt.0)
     +              sltmp=0.0
               idisary(iz,K) = int(sltmp*scn)
               if (idisary(iz,K).gt.imxscl) idisary(iz,K) = imxscl
               aryidl(iz,K) = sltmp

C*******************************************************************
C       save max value for color coded scale in plot* subroutine 
C*******************************************************************
               sl_max = AMAX1(sl_max,sltmp)

c*****************************************************************************
c       NOTE HISTOGRAM ARRAY INWARY FILLED HERE SAME ORDER AS MAIN DISPLAY ROWS
c*****************************************************************************
               inwary(iz,K)=int(sltmp) ! note NOT SCALED & CONVERTED TO INTEGER
c       next line replaced by line in plot routine
c       if (idisary (iz,K).eq.0) idisary (iz,K)=1 ! for color = bkgrnd
44480       continue
44490    continue
44500 continue
c
      izz = iz
c       call sub to write out aryidl array
      outfile='idlmov.out'
      call idl (aryidl,outfile,izz,inc,imask,
     +     biggest)
c
c
c
c
      title='All movements -> reduced distance.'
      call mode(4)
c
      iline=0                   ! thin lines, set to 1 for fat line
      loop=1

c       ************************************************************
c       match pairs with idisary variables
c       ************************************************************
      j=0
      do II=2,N
         do JJ=1,II-1
            j=j+1
            icell1(j)=II
            icell2(j)=JJ
         end do
      end do

      call plot_op3(fildes,idisary,ixrange,iyrange,title,
     +     N,loop,iopt,iline,inwary,sl_max,xc,'y',
     +     icell1,icell2,iyrange)

c       print 4446
c4446   format (2x,'<CR> to continue...')
c       read (*,'(A)') inst
      goto 506
c       
      end
c       *********************************************************************
c       **********************************************************************
c       *************************** SUBROUTINES ***************************
c       *******************************************************************
c
      subroutine idl (aryidl,outfile,izz,inc,imask,
     +     biggest)
      
      parameter (npts=1000,ipairs=2016) ! ipairs =ncell*(ncell-1)/2
c       npts = max num time steps in orig. grav calc.
c       ipairs=number of unique pairs for ncell
      character*30 outfile
      dimension aryidl(ipairs,npts-1)
      dimension imask(ipairs)
c
c       next thre values are place holders for file compatibility
      ismooth=0                 !to deal with legacy option no longer available
      ihard=1                   !to deal with legacy option no longer available
      ijmp=0                    !to deal with legacy option no longer available

c       write out 'idl compatible files' file here - over writes file of same name!!!
      OPEN (UNIT=3,file=outfile,status='UNKNOWN',access='SEQUENTIAL',
     +     form='UNFORMATTED')
      write (3) ismooth         ! number of bins for gaus smooth option
      write (3) inc             ! number of plotted steps over which a slope calculated
      write (3) ihard           ! header to indicate sig test option used to gen data
      write (3) ijmp            ! header to indicate jmp sig test option used to gen data
      write (3) izz
      write (3) biggest         ! largest distance between pair in all controls
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
c




