        subroutine plane(itable, data1, arrnum, maxarr,fname)
c
c
c       This subroutine displays a SERIES of histogram PLANES in three
c       dimensions.  Each plane represents sparks at one gravity step.
c
c       v4. add line in z-plane for o-lag - tba
c          add reorder of display based on lag values in the *.gnew file
c       v6 add sort of histograms by lag value
c       v7 mod display to allow negatative as well as positive spark values.
c
c
c       **** get starbase aliases and constants ********************
        include 'config.defs'
        include 'sbparam.defs'
c       ************************************************************
c
c       place at top before parameters
        PARAMETER (nstep=999,ipairs=2016,itimr=101,itm=309) ! itm=3*itimr+6
c
        DIMENSION itable(ipairs),itable2(ipairs)
        DIMENSION iorder(ipairs)
        DIMENSION       data1(ipairs,nstep)
        integer arrnum  ! number of pairs for this run
        real hist(ipairs,itimr),maxarr
c
        integer*4 iopt
        character*30 fname
        real vertlist(itm),camera(13)
        character*1 reply,color
        if(.false.)print *,fname !suppress unused variable warning
        if(.false.)print *,maxarr !suppress unused variable warning

c
c
c       ##### prompt user for mode #####################
c
        call mode(3)
1       print 3
3       format(2x,'Routine displays 3-D "spark plane" histograms.',
     +  /,'There are six modes of operation.  They are:      ',
     +  /,'  1-single step',
     +  /,'  2-single step (sorted)',
     +  /,'  3-auto sequence - loop',
     +  /,'  4-auto sequence - loop (sorted)',
     +  /,'  5-sum plane-watch it grow',
     +  /,'  6-sum plane-watch it grow (sorted)',/)

c
        print 2
2       format(2x,/,'Enter OPTION NUMBER : ',$)
        read (5,fmt='(i1)',err=1)iopt
c
c


        if((iopt.lt.1).or.(iopt.gt.6)) goto 1
c
        icflg=0
        if(iopt.eq.1) icflg=1
        if(iopt.eq.2) icflg=2
        if(iopt.eq.3) icflg=3
        if(iopt.eq.4) icflg=4
        if(iopt.eq.5) icflg=5
        if(iopt.eq.6) icflg=6




c
c
c       init itable2 with itable values
c
        do i=1,ipairs
        itable2(i)=itable(i) ! itable2 is used in sorting histogram by lag values
        end do


c       initialize the camera settings..................
c
c
c
c
c       **************************************************************
c       begin main body of 3d graphing code
c       **************************************************************
c
c
        call gclose(fildes)
        write(6, fmt ='(a4)') char(27)//'h'//char(27)//'J' !clear the screen

C       NEXT 4 lines create custom window
C
        fildes = gopen (1024,768,-700,-5,'hires')
        if (fildes.eq.-1)stop
        color="c"
C
C
        ZC=(arrnum-1)*1.0
C
c
c       initialize the camera settings here
c
        factor=1.5 ! = "distance" from image on screen
        xref=0.5
        yref=0.5
        zref=0.5
        theta=degrad(20.)
c       theta=degrad(15.)
        phi=degrad(45.)
        grad=degrad(5.0)
c
c
c
c       set camera parameters here
c
        radius=factor
        camera(INT(CAM_CAMX))=radius*COS(phi)*SIN(theta)+xref
        camera(INT(CAM_CAMY))=radius*SIN(phi)+yref
        camera(INT(CAM_CAMZ))=-(radius*COS(phi)*COS(theta))+zref
        camera(INT(CAM_REFX))=xref 
        camera(INT(CAM_REFY))=yref 
        camera(INT(CAM_REFZ))=zref
        camera(INT(CAM_FIELD_OV))=60.
        camera(INT(CAM_UPX))=-SIN(theta)*SIN(phi)
        camera(INT(CAM_UPY))=COS(phi)
        camera(INT(CAM_UPZ))=SIN(phi)*COS(theta)
        camera(INT(CAM_PROJECTION))=CAM_PARALLEL
        camera(INT(CAM_FRONT))=0.
        camera(INT(CAM_BACK))=0.
c
        call view_camera(fildes,camera)
c
        call mode(3)
        call gclose(fildes)
        write(6, fmt ='(a4)') char(27)//'h'//char(27)//'J' !clear the screen
        fildes = gopen (1024,768,-700,-5,'hires')
c
c
        if((color.eq."p").or.(color.eq."P")) then
          call setmono(fildes)
        else
          call colors(fildes)
        end if
c
c
        call set_p1_p2(fildes,FRACTIONAL,0.0,0.0,0.0,1.0,1.0,1.0)
        call vdc_extent(fildes,0.0,0.0,0.0,1.25,.95,1.0)
        call view_camera(fildes,camera)

c       ***********************************************************
c       GET READY FOR BIG LOOP for successive planes - one for each gravity step
c       *********************************************************************
c
c
c
c       ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c       set up DEFAULT values of iorder pointer array...
c
        do j=1,arrnum
        iorder(j)=j
        end do

c         OPTIONAL reorder of histograms from largest to smallest
c         sort on lag values in itable2 (init=itable values)
c
        if((icflg.eq.2).or.(icflg.eq.4).or.(icflg.eq.6)) then
          do 55 i=1,arrnum
            largest=-50000
            do 50 j=1,arrnum
              largest=MAX(itable2(j),largest)
50          end do
c
            do 52 k=1,arrnum
              if(itable2(k).eq.largest) ibig = k
52          end do
c
            iorder(i)=ibig
            itable2(ibig)=-50000
55        end do
        end if
c
c       ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c
c
c
c
c
c       clear hist array before loop; needed for summing option

          do j=1,arrnum
          do k=1,itimr
          hist(iorder(j),k)=0.0
          end do
          end do
c
c
        speed = 0.050 ! fast - for timed delay to slow animation
c
c       init monitor of largest hist value-used in dynamic scaling of
c       'sum'plane plot as it 'grows'
c
        if ((icflg.eq.5).or.(icflg.eq.6)) then
          big=0.0
          do ic=1,nstep
          do j=1,arrnum
          do k=1,itimr
          if (itable(iorder(j)).eq.k) then
          hist(iorder(j),k)= hist(iorder(j),k)+data1(iorder(j),ic)
          if(ABS(hist(iorder(j),k)).gt.big) big = ABS(hist(iorder(j),k))
          else
          hist(iorder(j),k)=hist(iorder(j),k)
          end if
          end do
          end do
          end do
c
c       set big - if * 1.0 then get full max scale
c               - if *2.0 get 1/2 max scale, etc
c
          big= 1.5*big
c
c       have scale factor for 'growing sum plot'
c
c
          do j=1,arrnum
          do k=1,itimr
          hist(iorder(j),k)=0.0
          end do
          end do
          end if
c

c
c       HERE IS THE BIG LOOP....
c

c       **********************************************************
        do 6000 ic=1,nstep
c       ***********************************************************
c       reset the port
c       ***********************************************************
c       call clear_control(fildes,CLEAR_VIEWPORT)
        call clear_view_surface(fildes)
c
c
        if ((icflg.eq.5).or.(icflg.eq.6)) then
c         increment sum histogram here-opt 5 or 6
          do j=1,arrnum
          do k=1,itimr
          if (itable(iorder(j)).eq.k) then
          hist(iorder(j),k)= hist(iorder(j),k)+data1(iorder(j),ic)
          else
          hist(iorder(j),k)=hist(iorder(j),k)
          end if
          end do
          end do
        else
c         fill histogram array for this loop here; opt 1 through 4

          do j=1,arrnum
          do k=1,itimr
          if (itable(iorder(j)).eq.k) then
          hist(iorder(j),k)= data1(iorder(j),ic)
          else
          hist(iorder(j),k)=0.0
          end if
          end do
          end do
        end if

c       ok - histogram updated for all particular icflg options



C
c       draw the frame
c
        call move3d(fildes,0.0,0.0,1.0)
        call draw3d(fildes,1.0,0.0,1.0)
        call draw3d(fildes,1.0,0.0,0.0)
        call draw3d(fildes,0.0,0.0,0.0)
        call draw3d(fildes,0.0,0.0,1.0)
c
        call move3d(fildes,0.0,0.0,1.0)
        call draw3d(fildes,0.0,1.0,1.0)
c       draw origin line next..
        call move3d(fildes,0.51,0.0,0.0)
        call draw3d(fildes,0.51,0.0,1.0)

c       plot the data points for ONE PLANE (velocities over one plotted gravity step)
c
        call interior_style(fildes,INT_SOLID,1)
c
        zeeax=1.0
c


        do 100 i=1,arrnum
c         if((color.eq."p").or.(color.eq."P")) then
c           shade = 0.3*MOD(i,2)
c           call fill_color(fildes,shade,shade,shade)
c         else
            call fill_color_idx(fildes,i+1)
c         end if
          vertlist(1)=0.0 ! x
          vertlist(2)=0.0 ! y
          vertlist(3)=zeeax
          do 95 j=1,itimr
            vertlist((j-1)*3+4)=j*.01
        if ((icflg.eq.5).or.(icflg.eq.6))then
            vertlist((j-1)*3+5)=(hist(iorder(i),j))/big
        else
            vertlist((j-1)*3+5)=hist(iorder(i),j)
            end if
            vertlist((j-1)*3+6)=zeeax
95        end do
          vertlist(3*itimr+4)=1.0 ! x
          vertlist(3*itimr+5)=0.0 ! y
          vertlist(3*itimr+6)=zeeax
c
          call polygon3d(fildes,vertlist,103,0)
          call make_picture_current(fildes)   !dump the buffer
          zeeax=zeeax-1.0/ZC
100     end do
c
c
c
c
        if ((icflg.eq.1).or.(icflg.eq.2)) then
        print 171, ic+1
171     format(2x,'Type ENTER for next plane (#',i5,')')
        read(*,'(A)')reply
        end if

c       *****************************************
c       PAUSE IN DISPLAY TO CONTROL MOVIE SPEED
c       ****************************************
        time2=0.0
78547   s=secnds(0.0)
        s=secnds(s)
        time2=time2+s
        if(time2.lt.speed) goto 78547
c

6000  continue

c
c       ****************************************************************
c       end of main body of code for the 3d section
c       ****************************************************************
c
c
        call clear_control(fildes,CLEAR_DISPLAY_SURFACE)
        call clear_view_surface(fildes)
c
c
        return
        end
C
c
c
c       **********************************************
c       functions used in  code .. commented out because duplicated
c       in main program
c       **********************************************
c
c       a function for converting from degrees to radians
c
c       real FUNCTION degrad(degree)
c
c       real degree, pi
c       parameter (pi=3.1415927)
c
c       degrad=degree*pi/180.
c       return
c       end
c
c
c
c
c       a function to convert from radians to degrees
c
c       real FUNCTION raddeg(radian)
c
c       real radian,pi
c
c       parameter (pi=3.1415927)
c
c       raddeg=radian*180./pi
c       return
c       end
c
c
c
        subroutine setmono(fildes)
c
c       This is a routine which sets the color table for
c       monochrome dithered output.  This is used when
c       the the display is sent to the laser printer via
c       the screendump routine
c
c
        include 'config.defs'
        COMMON R(256),G(256),B(256)
c
c
c       create monochrome color table ++++++++++++++++++++++++
c
        R(1) = 0; G(1) = 0; B(1) = 0
        do i=2,64
           R(1) = 1; G(1) = 1; B(1) = 1
        end do
           
c
c       ++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
        call fill_dither(fildes,16)
c
        return
        end
