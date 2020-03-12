        subroutine surface(itable, data1, arrnum, maxarr,fname)
c
c
c       This subroutine displays a surface of histograms
c       - each histogram is a sum of the 'lag' sorted  AND compressed histograms
c        of main program option 2.
c       Each histogram represents a spark at one gravity step.
c
c       v2 fixes some bugs
c       v3 adds a simple 'contour' plane option.
c
c
c       **** get starbase aliases and constants ********************
        include 'config.defs'
        include 'sbparam.defs'
c       ************************************************************
c
        PARAMETER (nstep=999,ipairs=120,itimr=101,itm=309) ! itm=3*itimr+6

c
        DIMENSION itable(ipairs),itable2(ipairs)
        DIMENSION iorder(ipairs)
        DIMENSION data1(ipairs,nstep)
        integer arrnum  ! number of pairs for this run
        real hist(nstep,itimr),maxarr
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
     +  /,'There is one mode of operation.  It is:      ',
     +  /,'  1-sum plane, theta=15.0',/,
     +  /,'  2-sum plane, theta=20.0',/,
     +  /,'  3-contour plane, theta=15.0',/,
     +  /,'  4-contour plane, theta=20.0',/)

c
        print 2
2       format(2x,/,'Enter OPTION NUMBER : ',$)
        read (5,fmt='(i1)',err=1)iopt
c
c
        if((iopt.lt.1).or.(iopt.gt.4)) goto 1
        icflg=0
        if(iopt.eq.1) icflg=1
        if(iopt.eq.1) thetaval=15.0
        if(iopt.eq.2) icflg=1
        if(iopt.eq.2) thetaval=20.0
        if(iopt.eq.3) icflg=2
        if(iopt.eq.3) thetaval=15.0
        if(iopt.eq.4) icflg=2
        if(iopt.eq.4) thetaval=20.0
c
c
c
c
        if(icflg.eq.2) then
        print 200
200     format(2x,'Enter contour cut off - % max',/,
     +  'e.g., .25 shows all between 25 and 100%',$)
        read (5,fmt='(F10.4)',err=1)pmax
        end if
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
        fildes = gopen (1024,768,-700,-5,'hires')
        if (fildes.eq.-1)stop
        color="c"
C
C
        ZC=(nstep-1)*1.0
C
c
c       initialize the camera settings here
c
        factor=1.6 ! = "distance" from image on screen
        xref=0.5
        yref=0.5
        zref=0.5
        theta=degrad(thetaval)
        phi=degrad(45.)
        grad=degrad(5.0)
c
c
c
c       set camera parameters here
c
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
        call vdc_extent(fildes,0.0,0.0,0.0,1.25,.95,2.0)
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
          do 55 i=1,arrnum
            largest=-50000
            do 50 j=1,arrnum
              largest=MAX(itable2(j),largest)
50          end do
            do 52 k=1,arrnum
              if(itable2(k).eq.largest) ibig = k
52          end do
            iorder(i)=ibig
            itable2(ibig)=-50000
55        end do
c
c       ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c
c
c
c
c
c
c       clear hist array before loop; needed for summing option

          do j=1,nstep
          do k=1,itimr
          hist(j,k)=0.0
          end do
          end do
c
c       fill hist array and set scale factor at the same time
c
          big=0.0
          do ic=1,nstep
          do j=1,arrnum
          do k=1,itimr
          if (itable(iorder(j)).eq.k) then
          hist(ic,k)= hist(ic,k)+data1(iorder(j),ic)
          if(ABS(hist(ic,k)).gt.big) big = ABS(hist(ic,k))
          else
          hist(ic,k)=hist(ic,k)
          end if
          end do
          end do
          end do
c
c
        sclfac=2.0 
        big=sclfac*big ! this will cut max peak size - later in plotting rtn

        if(icflg.eq.2) then
        thresh=(big/sclfac)*pmax! for contour plane plot
        end if
c
c
c       ***********************************************************
c       reset the port
c       ***********************************************************

        call clear_view_surface(fildes)

c       HERE IS THE BIG LOOP....
c

c       **********************************************************

c       draw the frame
c
        call move3d(fildes,0.0,0.0,2.0)
        call draw3d(fildes,1.0,0.0,2.0)
        call draw3d(fildes,1.0,0.0,0.0)
        call draw3d(fildes,0.0,0.0,0.0)
        call draw3d(fildes,0.0,0.0,2.0)
c
        call move3d(fildes,0.0,0.0,2.0)
        call draw3d(fildes,0.0,0.5,2.0)
c       draw origin line next..
        call move3d(fildes,0.51,0.0,0.0)
        call draw3d(fildes,0.51,0.0,2.0)

c       plot the surface for all gravity steps)
c
        call interior_style(fildes,INT_SOLID,1)
        zeeax=2.0
            call fill_color_idx(fildes,0)
c
c       LOOP HERE
          do ic=1,nstep
          vertlist(1)=0.0 ! x
          vertlist(2)=0.0 ! y
          vertlist(3)=zeeax
          do 95 j=1,itimr
            vertlist((j-1)*3+4)=j*.01
            if (icflg.eq.2)then
            if(hist(ic,j).gt.thresh) 
c     + vertlist((j-1)*3+5)=(hist(ic,j)/big)*(1.0 - pmax)
     +      vertlist((j-1)*3+5)=(hist(ic,j)-thresh)/big
            if(hist(ic,j).le.thresh) vertlist((j-1)*3+5)=0.0
c
c
            else
c           vertlist((j-1)*3+5)=((hist(ic,j))/big)
            vertlist((j-1)*3+5)=hist(ic,j)/big
            end if
c
c
            vertlist((j-1)*3+6)=zeeax
95        end do
          vertlist(3*itimr+4)=1.0 ! x
          vertlist(3*itimr+5)=0.0 ! y
          vertlist(3*itimr+6)=zeeax
          call polygon3d(fildes,vertlist,103,0)
          zeeax=zeeax-2.0/ZC
          end do
c       END LOOP HERE
c
c
c
        call make_picture_current(fildes)   !dump the buffer
        print 171
171     format(2x,'Type any key to continue')
        read(*,'(A)')reply

c
c
c       ****************************************************************
c       end of main body of code for the this subroutine
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

