      program xslope
c       latest version of xslope bgl June 2009
c
c       *********************************************************
      INCLUDE 'head_04.def'
c       *************************************************************
c       DESECRIPTION: Plots "time vs. distance" for each pair in a
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
c       steps can be calculated... e.g., slope over every three steps will result in wider bins.
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
c       add "menu" so don''t have to do all options.. debugging as of 15-aug-95..bl
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
c       ****************-DEFINITION OF HEADER BLOCK STRUCTURE-************
c       parameter (ncl=64)!number of codes handled
c       parameter (nsz=480000)!size of charge history array - 1st dimension
c       parameter (nsppc=10000)!max number of spikes in the busiest channel
c       parameter (lstsz=640000)!Product of nsz*nsppc
c

c     ncell=ncl ! same as parameter ncl
      parameter (NCFR=10000)
      parameter (ncell=64)
      parameter (npts=1000)
      parameter (ipairs=2016)   ! =ncell*(ncell-1)/2
      parameter (sigcolor=2)    ! 2 = red
      parameter (maxcolor=13)   ! 13 = blue
c     
c     NCFR =  number of control files
c     ...defines sig. level: 100 = .01 , 1000 = .001 etc.
c     ncell = max number of cells in data file
c     npts = max num time steps in orig. grav calc.
c     ipairs=number of unique pairs for ncell
c     jmpval=integer =number of concurrent "jumps" to be considered JUMP
c     
c     min and max intensity levels in scaled plots..
      parameter (is_size=40,imnscl=1,imxscl=32)
     
      character*1 inst
      character*16 PLOT
      character*30 LABEL
      character*30 Strng
      character*60 title
      character*30 posfile
      character*30 outfile
      character*30 surfil (ipairs) ! array of NCFR file names
      character*30 confil (NCFR) ! array of NCFR file names
      DIMENSION XXX(npts+1),YYY(npts+1)
      DIMENSION PDIS(ncell,ncell,npts)
      DIMENSION CDIS(ncell,ncell,npts)
c     new arrays:
c     sl contains "raw" slope (velocity) values for each
c     pair of particles following each step.
c     amxsl contains max. slope (velocity) values for each step
c     for each pair of CONTROL particles over all control samples
c     isig contains 1 at steps where each pair exhibits "sig" aggregation
c     according to criterion in use for that "run".
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
c     dimension diff_arr(ncell,ncell,npts)
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

      do
         call main_menu
      end do
      
      CONTAINS

c     ******************************************************************
      SUBROUTINE MAIN_MENU
c     ******************************************************************
c     choices for this run...
c     iline defined in gaus subroutine and passed to newplot2 to define line width
c     .. changed from default when gaus smooting selected - better image on screen than no smooth
c     .. nosmooth uses default line width.

c     options for color display - facilitate various hard copy possibilities
      do; print "(2x,'BLACK BACKGROUND: 1-blue-red-white  2-spectrum',/,
     +2x, 'WHITE BACKGROUND: 3-gray scale  4 - spectrum2',
     +2x,/,'5 - black & spec for pc')"
         read (*,*,iostat=ios) iopt
         if (ios.eq.0.and.iopt.ge.1.and.iopt.le.5) exit
      end do
      
c     option for "bin" value in slope plots 
      do; print "(2x,'Calc. slope for every Nth step, (def=1):',$)"
         read (*,'(i9)',iostat=ios) inc
         if (ios.eq.0.and.inc.le.npts/2.and.inc.ge.0) exit
      end do
      if (inc.eq.0) inc = 1

      do; print "(2x,'1- single file  2- read controls')"
         read (*,*,iostat=ios) icflag
         if (ios.eq.0.and.icflag.ge.1.and.icflag.le.2) exit
      end do

      if (icflag.eq.2) then ! How many cycles?
         do; print "(2x,'1- 100 *.rdt files; 2- 1000, 3 = 10K')"
            read (*,*,iostat=ios) iexp
            if (ios.eq.0.and.icflag.ge.1.and.icflag.le.3) exit
         end do
         ncfile = 10 * 10**iexp
      end if

      print *,' Be patient. Reading surrogate/control files...'

      if (icflag.eq.2) call initialize_controls

      print "(2x,'..input file? :',$)"
      read (*,'(A)') posfile
      if (posfile.eq.'')posfile = '2006-08-09P_pm10.pos'

      XCKB=0.0                  !init. var
      OPEN (UNIT=1,file=posfile,status='OLD',form='UNFORMATTED')
      read (1) CBLOCK,IBLOCK,RBLOCK
      read(1) TIMAT             !read in spike times

c     **************** READ DATA ***************************
      NMPTS=LFRAM               !# of points in file for each pair of cells 
      DO  K=1,NMPTS
         READ (1,iostat=ios) ((PDIS(II,JJ,K),JJ=1,II-1),II=2,N)
         if (ios.eq.iostat_end) exit
      end do
      CLOSE (1)

      XCORIG= (ENDTIME+DELTIM)/1000. !final state time in sec. of data file
      call mode(4)
      PRINT "(2X,'default timespan will be ',F5.1,' sec. or ?:',$)",XCORIG
      READ (5,"(F5.1)") XCKB
      XC = merge (XCORIG, XCKB, XCKB.eq.0.0)
      if (XC.lt.XCORIG) NMPTS= int((LFRAM*(XC/XCORIG)))

      fildes=gopen(1600,1120,0,0,'xslope')
      if (fildes.eq.-1)stop
      call colors(fildes)

c     get number of lines (unique pairs) in current display
      iyrange = int (1_8 * N * (N - 1) / 2)
      
c     NMPTS - 1 is number of points plotted along x- axis 
      ixrange = NMPTS - 1

c     **************************************************************
c     fill sl array here
c     IF inc > 1 then inc successive bins have same value..
c     allows slope over several steps to be calculated..
c     DEFAULT INC MUST BE SET TO 1
c     **************************************************************
      sl = 0
      do II=2,N
         do JJ=1,II-1
            DO K=1,NMPTS-inc,inc
               cursl = PDIS(II,JJ,K)-PDIS(II,JJ,K+inc)
               do IZ=K,K+inc-1
                  if (IZ.eq. npts-1) exit
                  sl(II,JJ,IZ) = cursl
               end do
            end do
         end do
      end do

c     fill isig array (has 1 less element than NMPTS)
      do II=2,N
         do JJ=1,II-1
            do K=2,NMPTS
               isig(II,JJ,K-1) = imnscl ! coolest color
               if (PDIS(II,JJ,K).lt.yyy_min(II,JJ,K)) then
                  isig(II,JJ,K-1) = imxscl ! hottest color
               end if
            end do
         end do
      end do

c     *******************************************************************
c     Main Menu
c     *******************************************************************
      do
         call mode(4)
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
         select case (inst)
         case ('1'); call pairwise_aggregation
         case ('2'); call significant_aggregation
         case ('3'); call all_movements
         case ('4'); return
         case ('5'); call exit
         end select
      end do

      end subroutine main_menu

c     ********************************************************************************
      SUBROUTINE INITIALIZE_CONTROLS
c     ********************************************************************************
      
      write (confil,"('sh',i0,'.pos')") (i,i=1,NCFR) ! generate surrogate pos file names to read in
c     Generate NAMES for surrogate velocity distribution files names for each pair
      write (surfil, "('sur',i0,'_',i0,'.sur')") ((ii,jj,jj=1,ii-1),ii=2,ncell)

c     Initialize some arrays
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

c     get max min and mean of sur controls at each step-for sig 'band" in PDFT plots
      biggest=0.0
      
      amxsl = 0
      do jz=1,ncfile
         print "(2x, 'Reading in control file:',1x,I5)", jz
         posfile= confil(jz)
         OPEN (UNIT=1,file=posfile,status='OLD',form='UNFORMATTED')
         read (1) CBLOCK,IBLOCK,RBLOCK !read in header info for this plot
         read(1) TIMAT          !read in spike times
         NMPTS=LFRAM            !# of points in file for each pair of cells 
         DO K=1,NMPTS
            READ (1,iostat=ios) ((CDIS(II,JJ,K),JJ=1,II-1),II=2,N)
            if (ios.eq.iostat_end) exit
         end do
         close (1)
         do II=2,N
            do JJ=1,II-1
               DO K=1,NMPTS
                  yca(II,JJ)=AMAX1(yca(II,JJ),CDIS(II,JJ,K))
                  yyy_max(II,JJ,K)=AMAX1(yyy_max(II,JJ,K),CDIS(II,JJ,K))
                  biggest=AMAX1(biggest,yyy_max(II,JJ,K))
                  yyy_min(II,JJ,K)=AMIN1(yyy_min(II,JJ,K),CDIS(II,JJ,K))
                  yyy_mean(II,JJ,K)=yyy_mean(II,JJ,K) + CDIS(II,JJ,K)
               end do
            end do
         end do
c     *************************************************
c     Fill amxsl array here
c     It contains the largest slope (out of ncfile control sets)
c     for each pair at each step.
c     IF inc > 1, then inc successive bins have same value.
c     This allows slope over several steps to be calculated..
c     *************************************************
         do II=2,N
            do JJ=1,II-1
               DO K=1,NMPTS-inc,inc
                  amxsl(II,JJ,K)=AMAX1(amxsl(II,JJ,K),
     +                 CDIS(II,JJ,K)-CDIS(II,JJ,K+inc))
                  tmpsl= amxsl(II,JJ,K)
                  do IZ=K,K+inc-1
                     if (IZ.eq. npts-1) exit
                     amxsl(II,JJ,IZ)=tmpsl
                  end do
               end do
            end do
         end do
      end do

      do; print "(9x,'1 - min, max & mean control plot (default)',/,
     $            9x,'2 - min only')"
         read (*,'(i9)',iostat=ios) iflag
         if (ios.eq.0.and.iflag.ge.1.and.iflag.le.2) exit
      end do
      if (iflag.eq.0) iflag = 1
      PLOT = merge ('min, max, & mean', 'min only        ', iflag.eq.1)

c     yyy arrays filled - now calculate correct mean values for yyy_mean
      do II=2,N
         do JJ=1,II-1
            DO K=1,NMPTS
               yyy_mean(II,JJ,K)=yyy_mean(II,JJ,K)/ncfile 
            end do
         end do
      end do
      end subroutine initialize_controls


c     ********************************************************************************
      SUBROUTINE PAIRWISE_AGGREGATION
c     ********************************************************************************
c     MAIN PROGRAM LOOP : TIME vs. DISTANCE PLOT FOR EACH PAIR
      call colors(fildes)
      call background_color(fildes,R(2),G(2),B(2))
      call line_color(fildes,R(1),G(1),B(1))
      call text_color(fildes,R(1),G(1),B(1))
      call perimeter_color(fildes,R(1),G(1),B(1))
      print "(2x,'Select particle to display ',
     +     '(enter 0 for all): ',$)"
      read (*,'(I4)') ij1
      isurplt=0
      do II=2,N
         do JJ=1,II-1
            isurplt=isurplt+1

            if (ij1 .ne. 0) then
               if (ij1 .ne. II .and. ij1 .ne. JJ) cycle
            end if

            call clear_control (fildes,CLEAR_DISPLAY_SURFACE)
            call clear_view_surface (fildes)

c        auto scaling of plot AND definition of global scaling for grvpat
            YC=0.0
            DO K=1,NMPTS
               YC= AMAX1(YC,PDIS(II,JJ,K))
            end do

            if (icflag .eq. 1 .or. iflag .eq.2) THEN
               YC=YC+0.05*YC
            else
               YC = AMAX1(YC,yca(II,JJ))
               YC=YC+0.05*YC
            end if

c       graph output parameters. These parameters define where on the 
c       screen the graphic output will appear.  All have a range of
c       0.0 to 1.0.
c                   Right - the amount of space to the right of graph
c                   Left  - the amount ofspace to the left of graph
c                   Top   - the amount of space above the graph 
c                   Botttom - the amount of space below the graph

            right=0.13
            left=0.15
            top=0.07
            bottom=0.35
            hrange=1-left-right
            vrange=1-top-bottom

            call vdc_extent(fildes,0.0,0.0,0.0,1.25,1.0,0.0)
            call view_port(fildes,left*1.25,bottom,(1-right)*1.25,
     +           (1-top))
            call view_window(fildes,0.0,0.0,XC,YC)
            call interior_style(fildes,INT_HOLLOW,Edge) !hollow rectangles
            call rectangle (fildes,0.0,0.0,XC,YC) !frame the viewport

            call make_picture_current(fildes) !dump the buffer
            LABEL = posfile
            lincol=0            !black
            m = lincol + 1
            call line_color(fildes,R(m),G(m),B(m)) 
            xstp=XC/NMPTS

            DO K=1,NMPTS    ! successive 'inc' steps same color 
               YYY(K) = PDIS(II,JJ,K)
               XXX(K)= K*xstp
               if (icflag.eq.2.and.sl(II,JJ,K).gt.amxsl(II,JJ,K)) then
                  icolr(K)=int(sigcolor)
               else 
                  icolr(K)=1
               end if
            end do

c     Make this line thicker.. because of coordinates looks distorted
c     verticle part looks thicker

            call line_width (fildes,0.002,VDC_UNITS) ! def = 0.0
            call line_endpoint (fildes,rounded) 
            call move2d(fildes,0.0,100.)
            do jz=1,NMPTS
               if (jz > 1.and.icolr(jz-1).eq.sigcolor) then
                  m = icolr(jz-1) + 1
                  call line_color(fildes,R(m),G(m),B(m))
               else
                  m = lincol + 1
                  call line_color(fildes,R(m),G(m),B(m)) 
               end if
               call draw2d(fildes,XXX(jz),YYY(jz))
            end do

c     PLOT MAX SURROGATE VELOCITIES WHEN VEL is SIG
c     Note - successive inc steps same color (= 1 eval step) 

c     ABOVE IT PLOT PAIR VELOCITY AT SAME STEPS - same color as PDFT plot
c     Note - successive inc steps same color (= 1 eval step) 

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
c     Line width and color back to default
            call line_width (fildes,0.0,VDC_UNITS) ! def = 0.0
            m = lincol + 1
            call line_color(fildes,R(m),G(m),B(m)) 

c     *************************************************************
c     **************************************************************
c     PLOT SIGNIFICANCE BAND
c     *************************************************
c     ...first plot max point in all control files
c     at each step..
c     *************************************************
c     lincol=9 ! medium green
            lincol=0            ! black
            m = lincol + 1
            call line_color(fildes,R(m),G(m),B(m))
            call line_type (fildes,DASH)

c       min only so skip this part
            if (iflag.ne.2) then
               call move2d(fildes,0.0,100.)
               do jz=1,NMPTS
                  call draw2d(fildes,XXX(jz),yyy_max(II,JJ,jz))
               end do
            end if
            
c     ********************************************************
c       plot min points in all control files at each step
c     ********************************************************
            call move2d(fildes,0.0,100.)
            do jz=1,NMPTS
               call draw2d(fildes,XXX(jz),yyy_min(II,JJ,jz))
            end do

            if (iflag.ne.2) then ! PLOT MEAN OF CONTROL POINTS AT EACH STEP
               lincol=0
               m = lincol + 1
               call line_color(fildes,R(m),G(m),B(m))
               call line_type (fildes,DOT)
               call move2d(fildes,0.0,100.)
               do jz=1,NMPTS
                  call draw2d(fildes,XXX(jz),yyy_mean(II,JJ,jz))
               end do
               call make_picture_current(fildes) !dump the buffer
            end if

c     *************************************************************
c     *************************************************************
c     A SLOPE CALC.
c     *********************************************************
c     plot x-tics
            call line_color(fildes,R(1),G(1),B(1))
            call line_type(fildes,SOLID)
            xmrk=0.0
            do while (xmrk.le.XC)
               call move2d(fildes,xmrk,0.0)
               call draw2d(fildes,xmrk,(.02*YC))
               xmrk=xmrk+(0.25*XC)
            end do
c     plot y -tics
            ymrk=0.0
            do while (ymrk.le.YC)
               call move2d(fildes,0.0,ymrk)
               call draw2d(fildes,(.02*XC),ymrk)
               ymrk=ymrk+(0.25*YC)
            end do
c
!-----do main title
            call text_font_index(fildes, 4) !sans serif font
            call text_alignment(fildes,TA_LEFT,TA_BOTTOM,0.0,0.0)
            call character_height(fildes,0.05) !9% of current 0-1 Y range
            call text2d(fildes,0.5,0.94,'Particles:'//char(0), VDC_TEXT,1)
            write(unit=strng,fmt='(2x,i2,x,a1)')II,char(0)
            call append_text(fildes,strng,VDC_TEXT,1)
            write(unit=strng,fmt='(2x,i2,x,a1)')JJ,char(0)
            call append_text(fildes,strng,VDC_TEXT,0)
!---- do x-axis title
            call text_font_index(fildes,4)
            call text_alignment(fildes,TA_CENTER,TA_BOTTOM,0.0,0.0)
            call character_height(fildes,.04)
            call text2d(fildes,1.25*(hrange/2.+left),bottom-.14*vrange,
     +           'Time (sec.)'//char(0),VDC_TEXT, EndOfLine)
!---- do y-axis title
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
            do ival=0,4
               write(unit=Strng,fmt='(F5.1,x,a1)')X,char(0)
               xaddr=(.25*ival*hrange+left)*1.25
               call text2d(fildes,xaddr,bottom-size/5.,Strng,VDC_TEXT,
     +              EndOfLine)
               X=X+(0.25*XC)
            end do
!-----label y tics
            call text_alignment(fildes,TA_RIGHT,TA_HALF,0.0,0.0)
            Y=0.0
            do jval=0,4
               write(unit=Strng,fmt='(f5.1,x,a1)')Y,char(0)
               yaddr=.25*jval*vrange+bottom
               call text2d(fildes,left*1.25,yaddr,Strng,VDC_TEXT,
     +              EndOfLine)
               Y=Y+(.25*YC)
            end do
            call make_picture_current(fildes) !dump the buffer
            call text_alignment(fildes,TA_LEFT,TA_BOTTOM,0.,0.)
            call text_path(fildes,PATH_RIGHT)
            call text_line_path(fildes,PATH_DOWN)
            call text_orientation2d(fildes,0.0,1.0,1.0,0.0)
            call character_height(fildes,0.012)
            call text_font_index(fildes,1)
            call vdc_extent(fildes,0.,0.,0.,1.25,1.0,0.)
            call view_window(fildes,0.0,0.0,20.0,20.0)
            call view_port(fildes,.12*1.25,0.,.88*1.25,.2)

c       spike counter for N neurons
            do idx=1,N
               tick=1
               do while((TIMAT(tick,idx).gt.0.0).and.
     +              (tick.le.nsppc))
                  tick=tick+1
               end do
               cnt(idx)=tick-1
            end do

c       particle count block output

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

c     ********************************************************
c       print settings
c     ********************************************************
            call text2d(fildes,1.092,.975,'PRG:'//char(0),
     +           VDC_TEXT,1)
            write(unit=strng,fmt='(a)')PROGID
            call append_text(fildes,strng,VDC_TEXT,0)

            call text2d(fildes,1.092,.950,'M-C:'
     +           //char(0),VDC_TEXT,1)
            write(unit=strng,fmt='(a)')PLOT
            call append_text(fildes,strng,VDC_TEXT,0)

            call text2d(fildes,1.092,.925,'FL:'
     +           //char(0),VDC_TEXT,1)
            call character_height(fildes,0.017)
            write(unit=strng,fmt='(a)')posfile
            call append_text(fildes,strng,VDC_TEXT,0)
            call character_height(fildes,0.020)

            call text2d(fildes,1.092,.900,'INT:'
     +           //char(0),VDC_TEXT,1)
            write(unit=strng,fmt='(i3,x,a6,a1)')IHS,"spikes",char(0)
            call append_text(fildes,strng,VDC_TEXT,0)

            call text2d(fildes,1.092,.875,'END:'
     +           //char(0),VDC_TEXT,1)
            write(unit=strng,fmt='(f10.1,x,a1)')ENDTIME,char(0)
            call append_text(fildes,strng,VDC_TEXT,0)

            call text2d(fildes,1.092,.850,'FSGN:'
     +           //char(0),VDC_TEXT,1)
            write(unit=strng,fmt='(f4.1,x,a1)')FS,char(0)
            call append_text(fildes,strng,VDC_TEXT,1)

            if(FS.eq.1.0) then
               call append_text(fildes,'excit'
     +              //char(0),VDC_TEXT,0)
            else
               call append_text(fildes,'inhib'
     +              //char(0),VDC_TEXT,0)
            end if

            if (ISHIFTA.EQ.0) THEN
               call text2d(fildes,1.092,.825,'ACC DECAY FWD'
     +              //char(0),VDC_TEXT,0)
            else
               call text2d(fildes,1.092,.825,'ACC DECAY BKWD'
     +              //char(0),VDC_TEXT,0)
            end if

            call text2d(fildes,1.092,.800,'NRM:'
     +           //char(0),VDC_TEXT,1)
            write(unit=strng,fmt='(i3,x,a1)')NORM,char(0)
            call append_text(fildes,strng,VDC_TEXT,0)

            if (ISHIFTE.EQ.0) THEN
               call text2d(fildes,1.092,.775,'EFF DECAY FWD'
     +              //char(0),VDC_TEXT,1)
            else
               call text2d(fildes,1.092,.775,'EFF DECAY BKWD'
     +              //char(0),VDC_TEXT,1)
            end if

            call text2d(fildes,1.092,.750,'NRM:'
     +           //char(0),VDC_TEXT,1)
            write(unit=strng,fmt='(i3,x,a1)')NORM,char(0)
            call append_text(fildes,strng,VDC_TEXT,0)

            call text2d(fildes,1.092,.725,'TM ST:'
     +           //char(0),VDC_TEXT,1)
            write(unit=strng,fmt='(f6.1,x,a1)')DELTIM,char(0)
            call append_text(fildes,strng,VDC_TEXT,0)

            call text2d(fildes,1.092,.700,'SLD:'
     +           //char(0),VDC_TEXT,1)
            write(unit=strng,fmt='(E8.2,x,a1)')SLIDE,char(0)
            call append_text(fildes,strng,VDC_TEXT,0)

            call text2d(fildes,1.092,.675,'FINC:'
     +           //char(0),VDC_TEXT,1)
            write(unit=strng,fmt='(f6.1,x,a1)')FWDINC,char(0)
            call append_text(fildes,strng,VDC_TEXT,0)

            call text2d(fildes,1.092,.650,'BINC:'
     +           //char(0),VDC_TEXT,1)
            write(unit=strng,fmt='(f6.1,x,a1)')BAKINC,char(0)
            call append_text(fildes,strng,VDC_TEXT,0)

            call text2d(fildes,1.092,.625,'FTAU:'
     +           //char(0),VDC_TEXT,1)
            write(unit=strng,fmt='(f6.1,x,a1)')FWDTAU,char(0)
            call append_text(fildes,strng,VDC_TEXT,0)

            call text2d(fildes,1.092,.600,'BTAU:'
     +           //char(0),VDC_TEXT,1)
            write(unit=strng,fmt='(f6.1,x,a1)')BAKTAU,char(0)
            call append_text(fildes,strng,VDC_TEXT,0)

            call text2d(fildes,1.092,.575,'CRDI:'
     +           //char(0),VDC_TEXT,1)
            write(unit=strng,fmt='(f6.1,x,a1)')CRDI,char(0)
            call append_text(fildes,strng,VDC_TEXT,0)

            call text2d(fildes,1.092,.550,'MS BT:'
     +           //char(0),VDC_TEXT,1)
            write(unit=strng,fmt='(f8.1,x,a1)')stimper,char(0)
            call append_text(fildes,strng,VDC_TEXT,0)

            call text2d(fildes,1.092,.525,'ST SVD:'
     +           //char(0),VDC_TEXT,1)
            write(unit=strng,fmt='(i4,x,a1)')ISTPFR,char(0)
            call append_text(fildes,strng,VDC_TEXT,0)

            call text2d(fildes,1.092,.500,'FR SVD:'
     +           //char(0),VDC_TEXT,1)
            write(unit=strng,fmt='(i5,x,a1)')LFRAM,char(0)
            call append_text(fildes,strng,VDC_TEXT,0)

            call make_picture_current(fildes) !dump the buffer

c     Write surrogate file option..  
            print "(5x,'<cr> continue or w for sur file for plot:',$)"
            read (*,'(A)') inst
            if((inst.ne.'w').and.(inst.ne.'W')) cycle

c     Write distributions for sig aggs in this pair's PDFT plot
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
                  READ (1,iostat=ios) 
     +                 ((CDIS(IB,JB,mz),JB=1,IB-1),IB=2,N)
                  if (ios.eq.iostat_end) exit
               end do
               close (1)
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
         end do
      end do
      end subroutine pairwise_aggregation

c     ********************************************************************************
      SUBROUTINE SIGNIFICANT_AGGREGATION
c     ********************************************************************************
c     SLOPE DISPLAYS
c     sequence of displays showing sig aggregation across all pairs
c     and velocity of aggregation

      if (icflag.eq.1) then
         print "(2x,'Must have control files for this option ',
     +           /,'         (<CR> to continue): ',$)"

         read (*,'(A)') inst
         return
      end if

c     ****************************************************
c     PLOT: TIMES OF SIG. AGGREGATION ACROSS PAIRS
c     ****************************************************
      title='Significant aggregation'
      loop=0                    !no loop through rows with reordering by max
      iline=1                   ! use fat lines for this plot

c     ************************************************************
c     match pairs with idisary variables
c     ************************************************************
      j=0
      do II=2,N
         do JJ=1,II-1
            j=j+1
            ijcell1(j)=II
            ijcell2(j)=JJ
         end do
      end do

      ips=0
      do II=2,N
         do JJ=1,II-1
            ips=ips+1
            do K=1,ixrange
               itmpary(ips,K)=isig(II,JJ,K)
            end do
         end do
      end do

      iq = 0
      do while (iq .lt. ips)
c     clear array for histogram plots      
         do iz = 1,ipairs
            do iy= 1, npts-1
               inwary(iz,iy)=0
            end do
         end do

c     set up idisary and icell1/2 top 40
         do iy=1,is_size
            do K=1,ixrange
               idisary(iy,K) = itmpary(iy+iq,K)
            end do
            icell1(iy) = ijcell1(iy+iq)
            icell2(iy) = ijcell2(iy+iq)
         end do
c     sl_max - max value for color key not used in this call 
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

         print "(2x,'<CR> to continue...')"
         read (*,'(A)') inst
      end do

c       *****************************************************************
c       *****************************************************************
c       insert routine here to differentiate  pairs (row) that never showed sig. aggregation
c       in PDFT plot. This arry is written out in idl*.out headers
c       for later use by grvpat's "MASK" option..
c       NOTE: if mask option 'on', then:1=MASK LATER...... 0= NO MASK LATER
       
      imask = 1                 ! default = mask this pair-DO NOT USE

      iz=0
      do II=2,N
         do JJ=1,II-1
            iz=iz+1
            isigflg=0
            do K=1,ixrange
               if (isig(II,JJ,K).eq.imxscl) isigflg= 1
            end do
c       next line: use this imask pair in grvpat if MASK on
            if (isigflg.eq.1) imask(iz)=0 !CLEAR MASK & USE IN GRVPAT*
         end do
      end do
      end subroutine significant_aggregation

c     ********************************************************************************
      SUBROUTINE ALL_MOVEMENTS
c     ********************************************************************************
c     PLOT: "RAW" neg. slopes ACROSS PAIRS - global normalization
      iz=0
      do II=2,N
         do JJ=1,II-1
            iz=iz+1
            rmax(iz)=0.0
            do K=1,ixrange
               rmax(iz)=AMAX1(rmax(iz),sl(II,JJ,K))
            end do
         end do
      end do

c     get global max slope over all pairs
      rmaxval=0.0
      do iz=1,iyrange
         rmaxval=AMAX1(rmaxval,rmax(iz))
      end do

      inwary = 0                ! clear array for histogram plots 

c     fill main display and histogram arrays
      sl_max = 0.0
      iz=0
      if (rmaxval.ne.0.0)scn=imxscl/rmaxval
      if (rmaxval.eq.0.0)scn=0.0
      do II=2,N
         do JJ=1,II-1
            iz=iz+1
            do K=1,ixrange
               sltmp=sl(II,JJ,K)
               if (sltmp.lt.0)
     +              sltmp=0.0
               idisary(iz,K) = int(sltmp*scn)
               aryidl(iz,K) = sltmp
               sl_max = AMAX1(sl_max,sltmp) ! save max value for color coded scale in plot* subroutine 
c     NOTE HISTOGRAM ARRAY INWARY FILLED HERE SAME ORDER AS MAIN DISPLAY ROWS
               inwary(iz,K)=int(sltmp) ! note NOT SCALED & CONVERTED TO INTEGER
            end do
         end do
      end do

      izz = iz
c     call sub to write out aryidl array
      outfile='idlmov.out'
      call idl (aryidl,outfile,izz,inc,imask, biggest)

      title='All movements -> reduced distance.'
      call mode(4)
      iline=0                   ! thin lines, set to 1 for fat line
      loop=1

c     match pairs with idisary variables
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
      end subroutine all_movements
      
      END PROGRAM XSLOPE
      
c     *******************************************************************
      SUBROUTINE IDL (aryidl,outfile,izz,inc,imask, biggest)
c     *******************************************************************
      parameter (npts=1000,ipairs=2016) ! ipairs =ncell*(ncell-1)/2
c     npts = max num time steps in orig. grav calc.
c     ipairs=number of unique pairs for ncell
      character*30 outfile
      dimension aryidl(ipairs,npts-1)
      dimension imask(ipairs)

c     next three values are place holders for file compatibility
      ismooth=0                 !to deal with legacy option no longer available
      ihard=1                   !to deal with legacy option no longer available
      ijmp=0                    !to deal with legacy option no longer available

c     write out 'idl compatible files' file here - over writes file of same name!!!
      OPEN (UNIT=3,file=outfile,status='UNKNOWN',access='SEQUENTIAL', form='UNFORMATTED')
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
      end subroutine idl
