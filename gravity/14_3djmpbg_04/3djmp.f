c     program 3djmp_04   
      character*1 inst
      inst = 'c'
      do while ((inst.eq.'c').or.(inst.eq.'C'))
         call threedjmp
         print "(2x,'C to continue, RETURN to exit: ',$)"
         read (*,'(A)') inst
      end do
      end program

c     ****************************************************************
      SUBROUTINE THREEDJMP
c     ****************************************************************
c     v. 5.0 add back wall for phase plane    bgl
c     and read in second analog ch to gen phase plane data
c     and increase # pairs to that expected for 32 neurons.
c     v.7.0 ch. color to ana side wall ch.1 with each spark..
c     ...just like phase plane changes color with each spark.
c     v. xslopebg - 10000 spikes; 1000 steps; 16 codes max.
c     NOTE: in this program nstep = one LESS that antecedent programs 
c     NA= max # analog samples in analog arrays
c     
c     v 3djmpbg2 - bgl - jan 1, 1998 (Happy new year)
c     add multiple data arrays to hold 3 (THREE) different sets of spark files
c     ...all MUST come from same gravity analysis because all are scaled globally
c     ... also read in spark files witha separate subroutine to facilitate the
c     multiple files..
c     ..also make various changes to allow different sets of sparks to be 
c     represented in different colors,... etc
c     3djmp_04 - handle 64 neurons, 2016 distinct pars (n*(n-1)/2)
      use color_mod
      implicit none
      integer NA, nstep, ipairs, iaend, iopta, ispln, iw, j,
     +     k, ianaflg, ispf, iyrange1, iyrange2, izrange,
     +     max1, max2, min1, min2, n, inc
      PARAMETER (NA=200000,nstep=999,ipairs=2016)
c     FOR PHASE PLANE:
c     next line determines amount of smooting in slope calc!!
c     must be odd number, middle array pt is location of instantaneous
c     slope.
c     ISM2 must = (0.5 *(ISM-1))+1
c     PARAMETER(scalep=0.8,xofsetp=.1,yofsetp=.1)
C     place at top before parameters
      include 'config.defs'
      include 'sbparam.defs'
      real scalep, xofsetp, yofsetp, pseg, sclfacxp,
     +     sclfacyp, timspan, zdis, zinc, zminp, zminp2,
     +     anstep, ax, rk, x, y, ysclfac1, ysclfac2, zmaxp, zmaxp2,
     +     degrad
      PARAMETER(scalep=0.9,xofsetp=.05,yofsetp=.05)
c     PARAMETER (ISM=7,ISM2=4)
c     PARAMETER (ISM=15,ISM2=8)
      integer ism, ism2
C      PARAMETER (ISM=21,ISM2=11)
c     end of phase plane stuff
      integer*4 param
      integer arrnum, sprksetcnt, i
      real            maxarr
      real            xref,yref,zref,radius,phi,factor
      real theta
      real            camera(13),grad,sign
      real            xmrk,ymrk,zmrk
      logical iset
      integer       ispark(3)
      real       data(ipairs,nstep,3)
      DIMENSION       iset(nstep, 3)
      real       anary1(NA,2),anary2(NA,2)
      integer       ianary1(NA,2),ianary2(NA,2)
      real       ana1(NA,2) ! X,Y for  element in phase plane
      integer       iana1(NA) ! X,Y,color for that element in phase plane
      character*1       char1, inst, change
      character*80      Strng
      character*30    string
      character*20    reply
      character*80    file(3)
      CHARACTER*30 fname,title
      real pi, edge, endofline
      parameter (pi=3.1415926)                 
      parameter (Edge=TRUE)     !used by "interior style"
      parameter (EndOfLine=FALSE) !used by "text2d"
      integer iaxes,iphase,ihiltp(3),
     $     iancol,itxtc,ibkgnd, ihilta(3)
      logical do_analog, do_phase, dynamic, movie
      real sprkitvl
      
      fildes=gopen (1024,768,135,203,'3djmp') ! creates window

      if (fildes.eq.-1)stop

      inst=' '
c     SELECT FROM 1 to 3 sets of sparks to display

      file=' '
      do
         print "('  How many (1-3) spark sets',$)"
         print "(' FROM SAME XSLOPE.OUT to display? ',$)"
         read (*,'(i5)')sprksetcnt 
         if ((sprksetcnt.ge.1).and.(sprksetcnt.le.3)) exit
      end do
      
      do i = 1,sprksetcnt
         print "('  Input data file',i0,' from xslopebg: ',$)",i
         read (*,'(A)') file(i)
      end do
     
c     initialize data* arrays before readin options
      data = 0
     
c     read in at least 1 spark data set and upto max allowed by this version - 
c     .. check number of arrays called data*
      do i=1,sprksetcnt
         call readspk(file(i),data(:,:,i:i),arrnum,inc)
      end do
     
      title=' '
     
c     get maxarr value to use in scaling all apark vectors to largest in all sets
c     ...read in NB: all sets must be from same XSLOPE output FILE !!!!!
c     SCALE SPARKS HERE FOR LATER DISPLAY
      maxarr = maxval (data)
      data = 0.75 * data / maxarr

c     clear array for analog data to be plotted on side wall
      ianary1=0
      ianary2=0
      ana1=0.0
      iana1=0               ! clear "add on" color elements in array
      anary1=0.0
      anary2=0.0
c     COLOR DEFINITIONS, table set by subroutine colors32..
c     .......black = 0 for background

      iaxes=10
      iphase=12
      iancol=24
      itxtc=29
      ibkgnd = 0
      ihiltp = (/18, 28, 7/)    ! phase plane highlite for 3 spark types
      ihilta = (/18, 28, 7/)    ! analog highlite for 3 spark types
      ispark = (/18, 28, 7/)

      do_analog = .true.
      do_phase = .false.
      PRINT "(2X,'BDT File for analog sig, <cr> to skip: ',$)"
      READ (*,'(A)') fname
      if (fname.eq.' ') do_analog = .false.
      if (do_analog) then
         print "(2x,'Draw the phase plane? ',$)"
         READ (*,'(A)') reply
         if(scan(reply(1:1),'yY').ne.0) do_phase = .true.
         call get_analog
         if (do_phase) call get_phase
      end if

c     GET READY TO DISPLAY ALL...
      dynamic = .true.
      print "(/2x,'Sparks will be displayed at 25 frames/sec')"
      print "(2x,'Enter animation speed in frames per spark step: ',$)"
      read (*,*) sprkitvl
      if (sprkitvl.le.1./nstep) dynamic = .false.
      movie = .false.
      if (dynamic) then
         print "(2x,'Save a movie to movie.avi? ',$)"
         READ (*,'(A)') reply
         if(scan(reply(1:1),'yY').ne.0) movie = .true.
      end if
      print "(2x,'title of plot: ',$)"
      read (*,'(A)') title 
      if (title.eq."") title = "TITLE OF PLOT"
      call colors32()
c     initialize the camera settings
      factor=1.6
      xref=0.5
      yref=0.5
      zref=0.5
      theta=degrad(30.)
      phi=degrad(35.)
      grad=degrad(5.0)

      call set_camera
      call view_camera(fildes,camera)
      call display_all
c     camera parameter adjustment
      do while(.true.)
         print "(2x, 'Do you wish to adjust the camera parameters
     +(0 to exit 3D)? ',$)"
         read(*,'(A)')reply
         if(scan(reply(1:1),'yY').ne.0) then
            call adjust
            call display_all
         else
            call gclose(fildes)
            return
         end if
      end do

      contains

c     ****************************************************************
      SUBROUTINE SET_CAMERA
c     ****************************************************************
      radius=factor
      camera(CAM_CAMX)=radius*COS(phi)*SIN(theta)+xref
      camera(CAM_CAMY)=radius*SIN(phi)+yref
      camera(CAM_CAMZ)=-(radius*COS(phi)*COS(theta))+zref
      camera(CAM_REFX)=xref
      camera(CAM_REFY)=yref
      camera(CAM_REFZ)=zref
      camera(CAM_FIELD_OV)=60.
      camera(CAM_UPX)=-SIN(theta)*SIN(phi)
      camera(CAM_UPY)=COS(phi)
      camera(CAM_UPZ)=SIN(phi)*COS(theta)
      camera(CAM_PROJECTION)=CAM_PARALLEL
      camera(CAM_FRONT)=0.
      camera(CAM_BACK)=0.
      end subroutine set_camera

c     ****************************************************************
      SUBROUTINE TICKS_AND_TEXT
c     ****************************************************************
c     plot x-tics
      xmrk=0.0
      do while (xmrk.le.1.0)
         call move3d(fildes,xmrk,0.0,0.0)
         call draw3d(fildes,xmrk,(.02),0.0)
         xmrk=xmrk+(0.05)
      end do

c     plot y -tics
      ymrk=0.0
      do while (ymrk.le.1.0)
         call move3d(fildes,0.0,ymrk,0.0)
         call draw3d(fildes,(.02),ymrk,0.0)
         ymrk=ymrk+(0.10)
      end do

c     plot z-tics
      zmrk=0.0
      do while (zmrk.le.1.0)
         call move3d(fildes,1.0,0.0,zmrk)
         call draw3d(fildes,0.985,0.0,zmrk)
         zmrk=zmrk+(0.05)
      end do

c     do main title
      call text_color_index(fildes,itxtc)
      call text_font_index(fildes, 4) !sans serif font
      call text_alignment(fildes,TA_CENTER,TA_TOP,0.0,0.0)
      call text_orientation3d(fildes,0.0,1.0,0.0,1.0,0.0,0.0)
      call character_height(fildes,0.055) !9% of current 0-1 Y range
      write(unit=Strng,fmt='(A)') title
      call text2d(fildes,0.625,1.0,trim(Strng),VDC_TEXT,EndOfLine)
c     do data file name
      call text_color_index(fildes,itxtc)
      call text_font_index(fildes, 4)         
      call text_alignment(fildes,TA_CENTER,TA_BOTTOM,0.0,0.0)
      call character_height(fildes,0.04)   
      write(unit=Strng,fmt='(A)')file(1)
      call text2d(fildes,0.625,0.0,trim(Strng),VDC_TEXT,EndOfLine)
c     do x-axis title
      call text_color_index(fildes,itxtc)
      call text_font_index(fildes,4)
      call text_alignment(fildes,TA_CENTER,TA_BOTTOM,0.0,0.0)
      call text_orientation3d(fildes,0.0,1.0,0.0,1.0,0.0,0.0)
      call character_height(fildes,.08)
      call text3d(fildes,0.4,-.15,
     +     0.0,'            ',WORLD_COORDINATE_TEXT, EndOfLine)

c     do y-axis title
      call text_orientation3d(fildes,-1.0,0.0,0.0,0.0,1.0,0.0)
c     call text3d(fildes,-.05,0.25,0.0,'            ',
      call text3d(fildes,-.05,0.25,0.0,'            ',
     +     WORLD_COORDINATE_TEXT, EndOfLine)

c     do z-axis title
      call text_orientation3d(fildes,-1.0,0.0,0.0,0.0,0.0,1.0)
      write(string,"(F7.1,a5)") timspan," SEC."
      call text_alignment(fildes,TA_RIGHT,TA_BOTTOM,0.0,0.0)
      call text3d(fildes,1.10,0.0,0.95,trim(string),
     +     WORLD_COORDINATE_TEXT, EndOfLine)
      end subroutine ticks_and_text

c     ****************************************************************
      subroutine plot_static
c     ****************************************************************
      integer istep
      call clear_view_surface (fildes)
      call draw_axes
      call ticks_and_text
      do istep = 1, nstep
         call draw_spark(istep)
      end do
      if (do_analog) then
         call draw_analog(nstep)
         if (do_phase) call draw_phase(nstep)
      end if
      end subroutine plot_static

c     ****************************************************************
      SUBROUTINE ADJUST
c     ****************************************************************
      real raddeg
      call gclose(fildes)
      fildes=gopen (1024,768,135,203,'3djmp')
      call colors32()
      call view_volume(fildes,-.5,-.5,-.5,1.5,1.5,1.5)

      call mapping_mode(fildes,DISTORT)
      call set_p1_p2(fildes,FRACTIONAL,0.,0.,0.,1.,.9,1.)
      call vdc_extent(fildes,0.0,0.0,0.0,1.25,1.0,1.0)
      call view_camera(fildes,camera)

      call draw_axes

      char1='a'
      do while(char1.ne.'q')
         print "(x,'Enter adjustments : ',$)"

c     print present camera position
         call text_alignment(fildes,TA_LEFT,TA_BOTTOM,0.,0.)
         call character_height(fildes,.03)
         call text2d(fildes,0.,0.,'  phi = ',VDC_TEXT,1)
         write(unit=Strng,fmt='(f5.1,x,a1)')raddeg(phi),char(0)
         call append_text(fildes,trim(Strng),VDC_TEXT,1)
         call append_text(fildes,'  theta = ',VDC_TEXT,1)
         write(unit=Strng,fmt='(f5.1,x,a1)')raddeg(theta),char(0)
         call append_text(fildes,trim(Strng),VDC_TEXT,1)
         call append_text(fildes,'  distance = ',VDC_TEXT,1)
         write(unit=Strng,fmt='(f5.1,x,a1)')radius,char(0)
         call append_text(fildes,trim(Strng),VDC_TEXT,0)
c     print camera adjustment menu
         call text_alignment(fildes,TA_LEFT,TA_TOP,0.,0.)
         call character_height(fildes,.03)
         call text2d(fildes,0.,1.,'r - right  ',VDC_TEXT,1)
         call append_text(fildes,'l - left   ',VDC_TEXT,1)
         call append_text(fildes,'u - up     ',VDC_TEXT,1)
         call append_text(fildes,'d - down   ',VDC_TEXT,1)
         call text2d(fildes,0.,.95,'t - toward ',VDC_TEXT,1)
         call append_text(fildes,'a - away   ',VDC_TEXT,1)
         call append_text(fildes,'q - quit   ',VDC_TEXT,0)
         call make_picture_current(fildes)
         
         read(*,'(a)')change
         if((change.eq.'u').or.(change.eq.'U')) then
            param=1
            sign=1.
         else if((change.eq.'d').or.(change.eq.'D')) then
            param=1
            sign=-1.
         else if((change.eq.'r').or.(change.eq.'R')) then
            param=2
            sign=1.
         else if((change.eq.'l').or.(change.eq.'L')) then
            param=2
            sign=-1.
         else if((change.eq.'t').or.(change.eq.'T')) then
            param=3
            sign=-1.
         else if((change.eq.'a').or.(change.eq.'A')) then
            param=3
            sign=1.
         else if((change.eq.'q').or.(change.eq.'Q')) then
            param=4
         else
c     do nothing
         end if

         if(param.eq.1) then
            phi=phi+grad*sign
         else if (param.eq.2)then
            theta=theta+grad*sign
         else if (param.eq.3)then
            factor=factor+0.1*sign
         else if (param.eq.4)then
            char1='q'
         end if
         call set_camera
         call clear_control(fildes,CLEAR_VIEWPORT)
         call clear_view_surface(fildes)
         call view_camera(fildes,camera)
         call draw_axes
         call make_picture_current(fildes)
      end do
      end subroutine adjust

c     ****************************************************************
      SUBROUTINE DRAW_AXES
c     ****************************************************************
      call line_color_index(fildes,iaxes)
      call move3d(fildes,0.0,0.0,0.0)
      call draw3d(fildes,0.0,1.0,0.0)
      call draw3d(fildes,0.0,1.0,1.0)
      call draw3d(fildes,0.0,0.0,1.0)
      call draw3d(fildes,0.0,0.0,0.0)
      call draw3d(fildes,1.0,0.0,0.0)
      call draw3d(fildes,1.0,0.0,1.0)
      call draw3d(fildes,0.0,0.0,1.0)
      call move3d(fildes,0.0,1.0,1.0)
      call draw3d(fildes,1.0,1.0,1.0)
      call draw3d(fildes,1.0,0.0,1.0)
      end subroutine draw_axes
      
c     ****************************************************************
      SUBROUTINE GET_ANALOG
c     ****************************************************************
      implicit none
      ispf=1                    !default color analog at spark times
      PRINT"(2X,'Analog plots change color with sparks(Y-def, N)?: ',$)"
      READ (*,'(A)') inst
      if((inst.eq.' ').or.(inst.eq.'y').or.(inst.eq.'Y')) ispf=1
      if((inst.eq.'n').or.(inst.eq.'N')) ispf=0
c     adjust color map for ispf selection
      if (ispf.eq.0) then
         ihilta=24
         ihiltp(1)=12
      end if
c     
      OPEN(UNIT=1,file=fname,FORM='FORMATTED',STATUS='OLD')
c     routine to read in analog values 
c     NOTE WELL: clock counts are integer!! have not been
c     multiplied by 0.5 to get msec.
      ianaflg=0
c     pointer to last filled array element
      iaend=0 
      call rana2(ianary1,ianary2,ianaflg,iaend,iopta)
      CLOSE (unit=1)
c     code to get peak value in current analog signal
      if (ianaflg.eq.1) then
         max1 = -32000
         min1 = 10000000
         max2 = -32000
         min2 = 10000000
c     find biggest and smallest analog value for scaling & offset on y axis
         do iw=1,iaend
            max1 = MAX0(ianary1(iw,1),max1)
            min1 = MIN0(ianary1(iw,1),min1)
            max2 = MAX0(ianary2(iw,1),max2)
            min2 = MIN0(ianary2(iw,1),min2)
         end do
c     now have largest int analog value(s)
      end if

c     set up to change writing color for ana ch. 1 & phase plane color at sparks      
      iset = maxval (data, 1) > 0.0

c     SCALE analog array 1 for display on side wall down Z axis. - 
c     half size if TWO signals are to be plotted.
      iyrange1 = max1-min1
      izrange= iaend 
c     will plot sig. 1 on side wall between .05 and .45; side wall goes up to 1.0 in y axis
c     add .05 to all values after scaled.
      ysclfac1=.4/iyrange1
      do iw=1,iaend
         anary1(iw,1)=((ianary1(iw,1)-min1)*ysclfac1)+0.05 ! voltage
         anary1(iw,2)=ianary1(iw,2)*.5 ! time in msec.
      end do
c     get time span of analog data on side wall
      timspan= (anary1(iaend,2) - anary1(1,2))/1000. ! in sec
      zinc=1.0/iaend            ! inc Z distance this much per analog point
      anstep= iaend/nstep       ! number of analog points per plotted gravity step.
      zdis=1.0/nstep            ! plotting distance/grav step
c     analog array 1 now ready to plot !!
c     plot second analog signal?
      if (iopta.eq.2) then
         iyrange2 = max2-min2
         izrange= iaend 
c     will plot sig. 2 on side wall between .55 and .95; side wall goes up to 1.0 in y axis
c     add .55 to all values after scaled.
         ysclfac2=.4/iyrange2
         do iw=1,iaend
            anary2(iw,1) = ((ianary2(iw,1)-min2)*ysclfac2)+0.55 ! voltage
            anary2(iw,2) = ianary2(iw,2)*.5 ! time in msec.
         end do
      end if
      end subroutine get_analog
      
c     ****************************************************************
      SUBROUTINE GET_PHASE
c     ****************************************************************
      implicit none
      integer maxlen
      real, allocatable :: SM(:)
      integer, allocatable :: IX(:)
      ispln=1                   !default solid line for phase plane always
      PRINT "(2X,'faint ph pl at non spk times(Y, N-def)?: ',$)"
c     this will change line_type from SOLID to DOT in phase plane
c     ... except for when spark is occurring
c     
      READ (*,'(A)') inst
      if((inst.eq.' ').or.(inst.eq.'n').or.(inst.eq.'N')) ispln=1
      if((inst.eq.'y').or.(inst.eq.'Y')) ispln=0
      
      maxlen = iaend - 1
      if (modulo (maxlen, 2).eq.0) maxlen = maxlen - 1
      PRINT "(2X,'smoothing kernel length (5-',I0,', odd, def=21): ',$)",maxlen
      ism = 0
      READ (*,'(I10)') ism
      if (ism .eq. 0) ism = 21;
      if (ism .lt. 5) ism = 5
      if (ism .gt. maxlen) ism = maxlen
      if (modulo(ism,2).eq.0) ism = ism + 1
      print "(2X,'smoothing kernel length set to ',I0)",ism
      ISM2  = ((ISM-1)/2)+1
      allocate (SM(ISM))
      allocate (IX(ISM))
c     FILL ARRAYS FOR PHASE PLANE
      pseg = (1.0 * iaend)/nstep 
C     we need the size of the smoothing array
C     a bigger array makes a smoother slope
C     initialize the smoothing array
      AX=0
      DO N=0,ISM-1
         RK=N*2*3.1416/(ISM-1)
         IF ((RK.GE.0.0).AND.(RK.LE.3.1416/2.0))
     +        SM(N+1)=(2.0*RK/3.1416)
         IF ((RK.GT.3.1416/2.0).AND.(RK.LE.3.0*3.1416/2.0))
     +        SM(N+1)=(2.0-2.0*RK/3.1416)
         IF ((RK.GT.3.0*3.1416/2.0).AND.(RK.LE.2.0*3.1416))
     +        SM(N+1)=(2.0*RK/3.1416-4)
         AX=AX+abs(SM(N+1))
      end do
      SM = SM / AX

c     get second variable for phase plane plot and store in ana1
C     Y is the slope of the stuff in IX
      do i=ISM2,iaend-ISM2+1
         do j=1,ISM
            IX(j)=ianary1(i-ISM2+j,1)
         end do
         Y=0.0
         X=0.0
         do k=1,ISM
            Y=Y-float(IX(k))*SM(k)
c     X=X+float(ABS(IX(k)*SM(k)))
            X=X+ABS(IX(k)*SM(k))
         end do
         ana1(i,1)=X
         ana1(i,2)=Y
c     set plot color for this point in phase plane array
         iana1(i)= iphase       ! base color of phase plane plots
         do k=1,3
            if (iset(int((i-1)/pseg)+1,k)) then
               iana1(i)= ihiltp(k) ! set point to ihiltp color
               exit
            end if
         end do
      end do
c     get zminp, zmaxp for x and y values
c     Find min & MAX PHRENIC AMPLITUDE ARRAY
      zminp=100000.
      zmaxp=-32000.
      zminp2=100000.
      zmaxp2=-32000.
      do j=ISM2,iaend-ISM2+1
         zminp = AMIN1(ana1(j,1),zminp)
         zmaxp = AMAX1(ana1(j,1),zmaxp)
         zminp2 = AMIN1(ana1(j,2),zminp2)
         zmaxp2 = AMAX1(ana1(j,2),zmaxp2)
      end do
c     scale factor
      sclfacxp = scalep/(zmaxp-zminp)
      sclfacyp = scalep/(zmaxp2-zminp2)
      end subroutine get_phase

c     ****************************************************************
      SUBROUTINE DISPLAY_ALL
c     ****************************************************************
      use ISO_C_BINDING
      implicit none
      integer istep_first, istep_last, nfrm, idur, ifrm, istep
      integer(c_int64_t) :: frmitvl = 40000 ! microseconds
      integer(c_int64_t) :: itvlsum
      integer :: slpsec = 1, is_movie
      
      call gclose(fildes)
      fildes=gopen (1024,768,135,203,'3djmp')
      call colors32()
      call mapping_mode(fildes,DISTORT)
      call set_p1_p2(fildes,FRACTIONAL,0.,0.,0.,1.,.9,1.0)
      call vdc_extent(fildes,0.0,0.0,0.0,1.25,1.0,1.0)
c     call clip_rectangle(fildes,-1.0,1.0,0.0,1.5)
c     call clip_indicator(fildes,CLIP_TO_VIEWPORT)
      call view_camera(fildes,camera)
      call background_color_index(fildes,0)
      call draw_axes
      call ticks_and_text
      if (.not.dynamic) then
         call plot_static
         call make_picture_current (fildes)
         return
      end if
c     plot the data points - dynamic
c     frame loop
      nfrm = floor ((nstep - 1) * sprkitvl) + 1
      idur = merge (1, floor (sprkitvl), sprkitvl < 1)
      if (movie) call init_movie
      call pace (0_c_int64_t)
      do ifrm = 1, nfrm
         call clear_view_surface (fildes)
         call draw_axes
         call ticks_and_text
         istep_first = ceiling ((ifrm - idur) / sprkitvl) + 1
         istep_last =  min (ceiling (ifrm / sprkitvl), nstep)

         istep_first = (istep_first - 1) / inc * inc + 1;
         istep_last = (istep_last - 1) / inc * inc + inc;
         if (istep_last.gt.nstep) istep_last = nstep;

         do istep = istep_first, istep_last
            call draw_spark(istep)
         end do
         if (do_analog) then
            call draw_analog(istep_last)
            if (do_phase) call draw_phase(istep_last)
         end if
         call make_picture_current (fildes)
         if (movie) call movie_frame (fildes)
         call pace (frmitvl)
      end do
      call sleep (slpsec)
      if (movie) then
         itvlsum = 0
         do while (itvlsum < slpsec * 1d6)
            call movie_frame (fildes)
            itvlsum = itvlsum + frmitvl
         end do
      end if
      call plot_static
      call make_picture_current (fildes)
      if (movie) then
         itvlsum = 0
         do while (itvlsum < slpsec * 1d6)
            call movie_frame (fildes)
            itvlsum = itvlsum + frmitvl
         end do
         call end_movie(is_movie)
      end if
      end subroutine display_all

c     ****************************************************************
      subroutine draw_spark(i)
c     ****************************************************************
      integer i, j, k
      real xx, yy, zz, ang
      zz = i/float(nstep)
      do k=1,3
         do j=1,arrnum
            if (data(j,i,k).ne.0.0) then
               call line_color_index(fildes,ispark(k))
               ang=pi*(j-1)/(arrnum-1)
               xx = 0.5+data(j,i,k)*cos(ang)
               yy = data(j,i,k)*sin(ang)
               call move3d(fildes,0.5,0.0,zz)
               call draw3d(fildes,xx,yy,zz)
            end if
         end do
      end do
      end subroutine draw_spark

c     ****************************************************************
      subroutine draw_analog(i)
c     ****************************************************************
      integer i, iw, k
      real zval, zval0
      zval=0.0
      do iw=1,iaend-1   
         zval0 = zval
         zval = 1.0 * iw / iaend
         call line_color_index(fildes,iancol)
         do k=1,3
            if(iset(int(zval/zdis)+1,k)) then ! if spark k at this step
               call line_color_index(fildes,ihilta(k))
            end if
         end do
         call move3d(fildes,0.0, anary1(iw,1),zval0)
         call draw3d(fildes,0.0, anary1(iw+1,1),zval)
         if (iopta.eq.2) then
            call move3d(fildes,0.0, anary2(iw,1),zval0)
            call draw3d(fildes,0.0, anary2(iw+1,1),zval)
         end if
         if (zval.gt.(i*zdis)) exit
      end do
      end subroutine draw_analog
      
c     ****************************************************************
      subroutine draw_phase(i)
c     ****************************************************************
      integer i, k
      real pxp, pyp
      pxp= ((ana1(ISM2,1)-zminp)*sclfacxp)+xofsetp
      pyp= ((ana1(ISM2,2)-zminp2)*sclfacyp)+yofsetp
      call move3d(fildes,pxp,pyp,1.0)
      do k=ISM2+1, int(i*pseg)
         if (k > iaend-ISM2+1) exit
         pxp= ((ana1(k,1)-zminp)*sclfacxp)+xofsetp
         pyp= ((ana1(k,2)-zminp2)*sclfacyp)+yofsetp
         if ((ispln.eq.0).and.(iana1(k).eq.iphase)) then
            call line_type(fildes,DOT)
         end if
         call line_color_index(fildes,iana1(k))
         call draw3d(fildes,pxp,pyp,1.0)
         call line_type(fildes,SOLID)
      end do
      end subroutine draw_phase

      end subroutine threedjmp

c     ****************************************************************
c     FUNCTION degrad()
c     a function for converting from degrees to radians
c     ****************************************************************
      real FUNCTION degrad(degree)
      real degree, pi
      parameter (pi=3.1415927)
      degrad=degree*pi/180.
      return
      end

c     ****************************************************************
c     FUNCTION raddeg()
c     a function to convert from radians to degrees
c     ****************************************************************
      real FUNCTION raddeg(radian)
      real radian,pi
      parameter (pi=3.1415927)
      raddeg=radian*180./pi
      return
      end

c     ****************************************************************
      SUBROUTINE READSPK(filename,data,arrnum,inc)
c     ****************************************************************
      PARAMETER (nstep=999,ipairs=2016)
      DIMENSION data(ipairs,nstep)
      DIMENSION  imask(ipairs)
c     DIMENSION  icnt(nstep)
      integer arrnum, inc
      character*20 filename
c     **************************************
c     next 3 lines give info about xslope32* sig tests for options 4 & 5..
c     see source code for details.. basically ihard =1 means a sig slope had
c     to be > any value in the pair row in all shifted sets; ihard=0 means a sig slope had to only be >
c     any value at the same step time in all shifted data sets.
c     
c     ismooth= number of columns over which gaussian calc..1, 2, or 4. 0=DEFAULT - no smoothing
c     inc = number of plotted steps over which a slope is calculated
c     ijmp=1 means pair had to have sig aggregation
c     in the PDFT plot to be considered for inclusion in a jump
c     jump column; ijmp=0 means ignore this sig test
c     
c     
      OPEN (UNIT=1,file=filename,status='OLD',access='SEQUENTIAL',
     +     form='UNFORMATTED')
      read (1,err=999,iostat=ios) ismooth ! header not yet used by this program
      read (1,err=999,iostat=ios) inc ! header not yet used by this program
      read (1,err=999,iostat=ios) ihard ! header not yet used by this program
      read (1,err=999,iostat=ios) ijmp ! not yet used.
      read (1,err=999,iostat=ios) arrnum
      read (1,err=999,iostat=ios) biggest
      do i=1,ipairs
c     read (1,err=999,iostat=ios) (imask(i))
         read (1,err=999,iostat=ios)  imask(i)
      end do
      do i=1,arrnum
         do j=1,nstep
            read (1,err=999,iostat=ios) data(i,j)
         end do
      end do

      CLOSE (unit=1)

      return
c     
c     simple read input file error handling...
c     
 999  print "(2x,'FATAL READ ERROR')"
      stop
      end
