c       **************************************************************
c       program hires2kv4   
c       **************************************************************
c       v. 1 - program shows single spark in 3-D; vectors go forward and back intime
c       using info from tuning options available in version of gbatch2kv9 or above
c       bgl - jan 2000.
c
c       main modules:
c
c       1. draw one spark in 3D with sparks that go forward or backward in time relative to origin
c                       in different shades of red, blue. - based on lags read in from offsets.gnew type file from
c                       gbatch2kv9 or higher.
c
c       2. versions 2,3: add 3-D spark "plane" option - in progress
c
c
c       3. version 4: add  'surface' option shows all 'nstep' histograms; 
c               ...all pairs lag sorted and compressed into 1 histogram
c
c       ipairs= total possible pairs for current max number of
c               cells allowed in gravity; e.g. for 64 neurons, ipairs=2016 (Formula= N*(N-1)/2)
c       nstep = number of plotted gravity time steps -1. - This is the number of steps for which distinct
c       aggregation velocities can be calculated
c
C       place at top before parameters
       
c
        PARAMETER (nstep=999,ipairs=2016,timr=101) 
c       timr = num bins in offset range histogram in gbatch_04* -- for zaxis plotting
c
        integer*4 param
        integer arrnum  
        real    maxarr
        real zzend(ipairs,nstep),xx(ipairs,nstep),yy(ipairs,nstep),ang
        real      xref,yref,zref,radius,theta,phi,factor
        real    camera(13),grad,sign
        real    xmrk,ymrk,zmrk
        DIMENSION itable(ipairs)
        DIMENSION       data1(ipairs,nstep)
        character*1       char1, inst, change
        character*30      Strng
        character*30    string
        character*20    reply
        character*30    file1 
        CHARACTER*30 fname,title
        include 'sbparam.defs'
        include 'config.defs'
        parameter (pi=3.1415926)                 
        parameter (Edge=TRUE)                   !used by "interior style"
        parameter (EndOfLine=FALSE)             !used by "text2d"
        
10000   fildes = gopen (1024,768,-700,-5,'hires')
        if (fildes.eq.-1)stop
        write(6, fmt ='(a4)') char(27)//'h'//char(27)//'J' !clear the screen
        inst=' '
c       *******************************
c
c       SELECT FROM 1 to 3 sets of sparks to display
c
c
c
        file1=' '

c
      print 1
1     format (2x,'Input data file1 from xslopebg: ',$)
      read (*,'(A)') file1
        
      print 2
2     format (2x,'Input offset file1 from gbatch2kv*: ',$)
      read (*,'(A)') fname

        do i=1,ipairs
        do j=1,nstep
        data1(i,j)=0.0
        end do
        end do
c
c
c       read in spark data set MUST DO THOS BEFORE call readoff to get arrnum value
c
c
        call readspk(file1,data1,arrnum)
c
c       get maxarr value to use in scaling all spark vectors to largest in all sets
c       ...read in NB: all sets must be from same XSLOPE output FILE !!!!!
c
c
c       read in tuning offsets for all pairs derived from gbatch2kv9 or higher
        call readoff(fname,itable,arrnum)
c
        title=' '
        maxarr=0.0
        do 200 i=1,arrnum
        do 150 j=1,nstep
        if (data1(i,j).gt.maxarr) maxarr=data1(i,j)
150     continue
200     continue
c
c       SCALE SPARKS HERE FOR LATER DISPLAY
c
        do 300 i=1,arrnum
        do 250 j=1,nstep
        data1(i,j)=0.75*data1(i,j)/maxarr 
250     continue
300     continue


c       *******************************************
c       COLOR DEFINITIONS, table set by subroutine colors32..
c       .......black = 0 for background

        iaxes=10
        iphase=12
        ispark=18
        itxtc=29
        ibkgnd = 0

c       *******************************************
c
c
c       GET READY TO DISPLAY ALL...
c       idflag=1 ! default is dynamic display
c
c
c       change next line when add dynamic options
        idflag=0
c
c
c
395     print 22207
22207   format (2x,' Choose display option:',/,
     +  2x,'   1 - 3-D spark',/,
     +  2x,'   2 - 3-D plane: lags & velocities',/,
     +  2x,'   3 - 3-D surface: collapsed sorted lags & velocities',/,
     +  2x,'  <ENTER> to EXIT!!')
      read (*,'(A)') inst                                               
        call mode(3) ! clear screen
        if (inst.eq.'1') goto 400
        if (inst.eq.'2') goto 398
        if (inst.eq.'3') goto 399
        if (inst.eq.' ') goto 4000
        goto 395

398     call plane(itable, data1, arrnum, maxarr,file1)
        goto 395 ! get next choice
399     call surface(itable, data1, arrnum, maxarr,file1)
        goto 395 ! get next choice


400     print 22208
22208   format (2x,' Choose spark number to display (1-999):',/,
     +  2x,' <CR> to exit')
      read (*,22209) isppl
      if (isppl.eq.0) goto 395
22209   format (i5) 
c
        write(6, fmt ='(a4)') char(27)//'h'//char(27)//'J' !clear the screen
      print 5
5     format (2x,'title of plot:',$)
      read (*,'(A)') title 

c       delay loop to permit mouse click to change color table
c      print 5987
c       5987  format (2x,'2 sec. to change color table with mouse..')
C
c      time2=0.0
c     28547 s=secnds(0.0)
c      s=secnds(s)
c      time2=time2+s
c      if(time2.lt.2.) goto 28547
C       end timed delay to slow animation

        call colors32(fildes)
c ***********************************************************
c       initialize the camera settings here
c ***********************************************************
        factor=1.6
        xref=0.5
        yref=0.5
        zref=0.5
        theta=degrad(30.)
        phi=degrad(35.)
        grad=degrad(5.0)


c ***********************************************************
c       set camera parameters here
c ***********************************************************
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
        call view_camera(fildes,camera)
        call mode(3)

c       
c
c
c
c
c
c ************************************************************
c
c BIG LOOP to do each spark from the file in turn, each in 3-D.
c
c ************************************************************
c
c
c
        i=isppl
c
        isflg=0
c
c       any sparks at this step?
        do jn1=1,arrnum
        if (data1(jn1,i).ne.0.0)isflg=1
        end do

        if (isflg.eq.0) then ! empty spark ..try again
      print 5122
5122     format (2x,'EMPTY SPARK LOCATION-try again..')
        goto 400
        end if 
c
c
c       THIS is reentry for SAME spark after changing camera view
c
c
c
174     call gclose(fildes)
        write(6, fmt ='(a4)') char(27)//'h'//char(27)//'J' !clear the screen
        fildes = gopen (1024,768,-700,-5,'hires')

        call colors32(fildes)

        call set_p1_p2(fildes,FRACTIONAL,0.,0.,0.,1.,.9,1.0)
        call vdc_extent(fildes,0.0,0.0,0.0,1.25,1.0,1.0)
        call mapping_mode(fildes,ISOTROPIC)

c       call clip_rectangle(fildes,-1.0,1.0,0.0,1.5)
c       call clip_indicator(fildes,CLIP_TO_VIEWPORT)

        call view_camera(fildes,camera)
        call background_color_idx(fildes,0)
        call line_color_idx(fildes,iaxes)
c ***********************************************************
c       draw the axes
c ***********************************************************
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

c ***********************************************************
c       plot x-tics
c ***********************************************************
        xmrk=0.0
        do while (xmrk.le.1.0)
        call move3d(fildes,xmrk,0.0,0.0)
        call draw3d(fildes,xmrk,(.02),0.0)
        xmrk=xmrk+(0.05)
        end do


c ***********************************************************
c       plot y -tics
c ***********************************************************
        ymrk=0.0
        do while (ymrk.le.1.0)
        call move3d(fildes,0.0,ymrk,0.0)
        call draw3d(fildes,(.02),ymrk,0.0)
        ymrk=ymrk+(0.10)
        end do


c ***********************************************************
c       plot z-tics
c ***********************************************************
        zmrk=0.0
        do while (zmrk.le.1.0)
        call move3d(fildes,1.0,0.0,zmrk)
        call draw3d(fildes,0.985,0.0,zmrk)
        zmrk=zmrk+(0.02)
        end do

c ***********************************************************
c       do main title
c ***********************************************************
        call text_color_idx(fildes,itxtc)
        call text_font_index(fildes, 4)         !sans serif font
        call text_alignment(fildes,TA_CENTER,TA_TOP,0.0,0.0)
        call character_height(fildes,0.055)      !9% of current 0-1 Y range
        write(unit=Strng,fmt='(A)') title
        call text2d(fildes,0.65,0.86,Strng,VDC_TEXT,EndOfLine)

c ***********************************************************
c       do data file name
c ***********************************************************
        call text_color_idx(fildes,itxtc)
        call text_font_index(fildes, 4)         
        call text_alignment(fildes,TA_CENTER,TA_TOP,0.0,0.0)
        call character_height(fildes,0.04)   
        write(unit=Strng,fmt='(A)')file1
        call text2d(fildes,1.3,0.07,Strng,VDC_TEXT,EndOfLine)
c ***********************************************************
c       do x-axis title
c ***********************************************************
        call text_color_idx(fildes,itxtc)
        call text_font_index(fildes,4)
        call text_alignment(fildes,TA_CENTER,TA_BOTTOM,0.0,0.0)
        call text_orientation3d(fildes,0.0,1.0,0.0,1.0,0.0,0.0)
        call character_height(fildes,.08)
        call text3d(fildes,0.4,-.15,
     +  0.0,'            '//char(0),WORLD_COORDINATE_TEXT, EndOfLine)

c ***********************************************************
c       do y-axis title
c ***********************************************************
        call text_orientation3d(fildes,-1.0,0.0,0.0,0.0,1.0,0.0)
        call text3d(fildes,-.05,0.25,0.0,'            '//char(0),
     +  WORLD_COORDINATE_TEXT, EndOfLine)

c ***********************************************************
c       do z-axis title
c ***********************************************************
        call text_orientation3d(fildes,-1.0,0.0,0.0,0.0,0.0,1.0)
       write(string,2001) "                  ", timspan," SEC."
       call text2d(fildes,0.90,.93,string,VDC_TEXT,0)
2001    format (a18,F7.1,a4)
        call text3d(fildes,1.10,0.0,0.45,string//char(0),
     +  WORLD_COORDINATE_TEXT, EndOfLine)
        call make_picture_current (fildes)
c ***********************************************************
c       PLOT 1 spark in 3-D
c ***********************************************************
c
c       draw hi res 3-D spark here
c
        call line_color_idx(fildes,ispark)
c       inner loop for one spark - one color
        do 350 j=1,arrnum
        j1=j-1
        z1=itable(j)-1
        angz=pi*z1/(timr-1)
        ang=pi*j1/(arrnum-1)
        xx(j,i)=0.5+data1(j,i)*cos(ang)
        yy(j,i)=data1(j,i)*sin(ang)
c       zzend(j,i)=0.5+itable(j)*cos(angz)
        zzend(j,i)=0.5+data1(j,i)*cos(angz)

c       write(*,*) i,j,xx(j,i),yy(j,i),zzend(j,i)
        call move3d(fildes,0.5,0.0,0.5)
        call draw3d(fildes,xx(j,i),yy(j,i),zzend(j,i))
350     continue

c
c       THE SPARK IS DRAWN
c
          call make_picture_current(fildes)   !dump the buffer
c
c ***********************************************************
c       camera parameter adjustment
c ***********************************************************
171     format(2x, 'Do you wish to adjust the camera parameters 
     +  (0 to exit 3D)? ',$)


        print 171
           read(*,'(A)')reply
           if (reply .eq. '0') goto 2999
        if((reply.eq.'y').or.(reply.eq.'Y')) then

           call gclose(fildes)
           write(6, fmt ='(a4)') char(27)//'h'//char(27)//'J' !clear the screen


        fildes = gopen (1024,768,-700,-5,'hires')
         call colors32(fildes)
           call view_volume(fildes,-.5,-.5,-.5,1.5,1.5,1.5)

           call set_p1_p2(fildes,FRACTIONAL,0.,0.,0.,1.,.9,1.)
           call vdc_extent(fildes,0.0,0.0,0.0,1.25,1.0,ZC/9.0)
           call mapping_mode(fildes,ISOTROPIC)
           call view_camera(fildes,camera)

           call line_color_idx(fildes,iaxes)
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

           call make_picture_current(fildes)   !dump the buffer

           char1='a'

           do while(char1.ne.'q')
           call mode(3)
           print 950
950       format(x,'Enter adjustments : ',$)

c ***********************************************************
c       print present camera position
c ***********************************************************
          call text_alignment(fildes,TA_LEFT,TA_BOTTOM,0.,0.)
          call character_height(fildes,.03)
          call text2d(fildes,0.,0.,'  phi = '//char(0),VDC_TEXT,1)
          write(unit=Strng,fmt='(f5.1,x,a1)')raddeg(phi),char(0)
          call append_text(fildes,Strng,VDC_TEXT,1)
          call append_text(fildes,'  theta = '//char(0),VDC_TEXT,1)
          write(unit=Strng,fmt='(f5.1,x,a1)')raddeg(theta),char(0)
          call append_text(fildes,Strng,VDC_TEXT,1)
          call append_text(fildes,'  distance = '//char(0),VDC_TEXT,1)
          write(unit=Strng,fmt='(f5.1,x,a1)')radius,char(0)
          call append_text(fildes,Strng,VDC_TEXT,0)

          call make_picture_current(fildes)   !dump the buffer

c ***********************************************************
c       print camera adjustment menu
c ***********************************************************
          call text_alignment(fildes,TA_LEFT,TA_TOP,0.,0.)
          call character_height(fildes,.03)
          call text2d(fildes,0.,1.,'r - right  '//char(0),VDC_TEXT,1)
          call append_text(fildes,'l - left   '//char(0),VDC_TEXT,1)
          call append_text(fildes,'u - up     '//char(0),VDC_TEXT,1)
          call append_text(fildes,'d - down   '//char(0),VDC_TEXT,1)
          call text2d(fildes,0.,.95,'t - toward '//char(0),VDC_TEXT,1)
          call append_text(fildes,'a - away   '//char(0),VDC_TEXT,1)
          call append_text(fildes,'q - quit   '//char(0),VDC_TEXT,1)
          call append_text(fildes,'<cr> - more'//char(0),VDC_TEXT,0)
          call make_picture_current(fildes)   !dump the buffer
        
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
c           do nothing
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


c ***********************************************************
c       reset camera position variables
c ***********************************************************
          radius=factor
          camera(CAM_CAMX)=radius*COS(phi)*SIN(theta)+xref
        camera(CAM_CAMY)=radius*SIN(phi)+yref
          camera(CAM_CAMZ)=-(radius*COS(phi)*COS(theta))+zref
          camera(CAM_REFX)=xref
          camera(CAM_REFY)=yref
          camera(CAM_FIELD_OV)=60.
          camera(CAM_UPX)=-SIN(theta)*SIN(phi)
          camera(CAM_UPY)=COS(phi)
          camera(CAM_UPZ)=SIN(phi)*COS(theta)
          camera(CAM_PROJECTION)=CAM_PARALLEL
          camera(CAM_FRONT)=0.
          camera(CAM_BACK)=0.
c ***********************************************************
c       reset the port
c ***********************************************************
          call clear_control(fildes,CLEAR_VIEWPORT)
          call clear_view_surface(fildes)
        call view_camera(fildes,camera)
c ***********************************************************
c       re-draw the axes
c ***********************************************************
      call line_color_idx(fildes,iaxes)
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
      call make_picture_current(fildes)   !dump the buffer

        end do

      call make_picture_current(fildes)   !dump the buffer
          call mode(3)
          goto 174 ! LOOP TO REPLOT SAME SPARK WITH NEW CAMERA VIEW
        end if ! IS THIS THE END OF THE IF reply .eq.'Y' at label 169 ??

c
c
c
c       END OF LOOP TO PLOT ONE OF THE SPARKS IN THE FILE
c
c
        goto 400        
c
c
c
4000    call mode(3)
        
        close(1)

        call gclose(fildes)
2999     print 11000
11000      format (2x,'C to continue, RETURN to exit:',$)
        read (*,'(A)') inst
        if ((inst.eq.'c').or.(inst.eq.'C')) goto 10000
        end

c
c
c       ******************************************************
c       ******************************************************
c       FUNCTIONS & SUBROUTINES
c
c       ******************************************************
c       ******************************************************

c       ****************************************************************
c       FUNCTION degrad()
c               a function for converting from degrees to radians
c       ****************************************************************

        real FUNCTION degrad(degree)

        real degree, pi
        parameter (pi=3.1415927)

        degrad=degree*pi/180.
        return
        end

c       ****************************************************************
c       FUNCTION raddeg()
c       a function to convert from radians to degrees
c       ****************************************************************

        real FUNCTION raddeg(radian)

        real radian,pi
        parameter (pi=3.1415927)

        raddeg=radian*180./pi
        return
        end

c       subroutine debug (line)
c       write (*,100) line
c 100   format (2x,'GOT HERE',i5)
c       return
c       end
c
c       *********************************************************************
        subroutine readspk(filename,data1,arrnum)
c       arrnum = number of pairs actually in the source gravity run..
        PARAMETER (nstep=999,ipairs=120)
        DIMENSION data1(ipairs,nstep)
        DIMENSION  imask(ipairs)
        integer arrnum  
        character*30 filename
c       **************************************
c       next 3 lines give info about xslope32* sig tests for options 4 & 5..
c       see source code for details.. basically ihard =1 means a sig slope had
c       to be > any value in the pair row in all shifted sets; ihard=0 means a sig slope had to only be >
c       any value at the same step time in all shifted data sets.
c
c       ismooth= number of columns over which gaussian calc..1, 2, or 4. 0=DEFAULT - no smoothing
c       inc = number of plotted steps over which a slope is calculated
c       ijmp=1 means pair had to have sig aggregation
c       in the PDFT plot to be considered for inclusion in a jump
c       jump column; ijmp=0 means ignore this sig test
c
c
        OPEN (UNIT=1,file=filename,status='OLD',access='SEQUENTIAL',
     +  form='UNFORMATTED')
        read (1,err=999,iostat=ios) ismooth ! header not yet used by this program
        read (1,err=999,iostat=ios) inc ! header not yet used by this program
        read (1,err=999,iostat=ios) ihard ! header not yet used by this program
        read (1,err=999,iostat=ios) ijmp ! not yet used.
        read (1,err=999,iostat=ios) arrnum
        read (1,err=999,iostat=ios) biggest
        do 1998 i=1,ipairs
        read (1,err=999,iostat=ios) imask(i)
1998    continue

         do 100 i=1,arrnum
         do 101 j=1,nstep

        read (1,err=999,iostat=ios) data1(i,j)
101     continue
100     continue

        CLOSE (unit=1)
        return
c
c       simple read input file error handling...
c
999     print 3999
3999    format(2x,'FATAL READ ERROR')
        stop
        end


c       *********************************************************************
        subroutine readoff(filename,itable,arrnum)
c       arrnum = number of pairs actually in the source gravity run..
        PARAMETER (ipairs=120) ! total possible pairs for current max number of
c               cells allowed in gravity, N=total allowed cells
        DIMENSION itable(ipairs)
        integer arrnum  
        character*30 filename
c       **************************************
c     READ in successive values of formatted text file offsets.gnew from gbatch2kv9 or later
c       format is I5,I5,I5 for: I,J (std grav nested loop variables), and offset (say, -50 to 50)
c
c
        OPEN (UNIT=1,file=filename,status='OLD',access='SEQUENTIAL',
     +  form='FORMATTED')


200     format(I5,I5,I5)
c
        do 50 m=1,ipairs ! N = MAX number of cells allowed in gravity
        itable(m)=0
50      continue
c
c
c       read in data
c

        do 100 i=1,arrnum
        read (1,200,END=300) ival,jval,noff ! throw away ival,jval in this program
c
c       next line saves offset value corresponding to a particular SINGLE SPARK after converting
c       from" -50 to 50" range to 1 to 101 range - for plotting purposes in this program 
c
c
        itable(i)=noff+51
100     continue
c
c
300     CLOSE (unit=1)
        return
        end
