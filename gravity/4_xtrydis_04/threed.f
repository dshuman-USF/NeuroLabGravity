c       **************************************************************
c       subroutine threed 
c               for displaying 3D representation of pairwise distances
c               can call with any of three options (inst):
c                       1 - auto display (5 sec between screens)
c                       2 - page (user must press enter to get next scr
c                       3 - single particle (user chooses particle to display
c       **************************************************************
 
        subroutine threed(fildes,NPTS,PDIS,XC,YC,icolor,TITLE,inst)
        INCLUDE 'head_04.def'
        include 'config.defs'
c       place at top before parameters
        EXTERNAL GETLOG_
        parameter (nstep=1000,ncell=64)
        integer*4 NPTS
        DIMENSION PDIS(ncell,ncell,nstep),icolor(32)
        real YC,XC
        real xref,yref,zref,radius,theta,phi,factor
        real camera(13),sign,grad,zeeax
c       real plot(260*3) 
        real plot((nstep+20)*3) 
        integer*4 c_cnt,plotpts,param
        DIMENSION  XXX3D(nstep+1),YYY3D(nstep+1)
        character*1 char1,inst,reply,change
        character*30 strng
        character*26 TITLE
        integer*4 EndOfLine

        INCLUDE 'sbparam.defs'

        parameter (EndOfLine=FALSE)             !used by "text2d"

        call gclose(fildes)
        call mode(0)
        fildes=gopen (1700,820,-700,-5,'xtrydis')
        if (fildes.eq.-1)stop
        call colors(fildes)
 
        ZC=N*(N-1)/2
        ZC2=.40*N*(N-1)/2
        call view_volume(fildes,-.5,-.5,-.5,1.5,1.5,ZC/9.0)

c ***********************************************************
c       initialize the camera settings here
c ***********************************************************
        factor=1.8
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

c **********************************************************
c       preliminary stuff for DO LOOP
c **********************************************************
709     format (2x,'Enter particle number, or 0 to exit 3D: ',$)
710     format (I2)

        jj = 1

c **********************************************************
c       big DO Loop for displaying ncell  particle pairs
c **********************************************************
        do while (jj .le. N)

        if (inst .eq. '3') then
           print 709
           read (*,710) jj
           if (jj .eq. 0) goto 999
        end if

174     call gclose(fildes)
        call mode (0)
        fildes = gopen (1700,820,-700,-5,'xtrydis')
        call colors(fildes)
        call set_p1_p2(fildes,FRACTIONAL,0.,0.,0.,1.,.9,1.)
        call vdc_extent(fildes,0.0,0.0,0.0,1.25,1.0,ZC/9.0)
        call mapping_mode(fildes,DISTORT)
        call view_camera(fildes,camera)
 
c ***********************************************************
c       draw the axes 
c ***********************************************************
        call move3d(fildes,0.0,0.0,0.0)
        call draw3d(fildes,1.0,0.0,0.0)

        call move3d(fildes,0.0,0.0,ZC2/9.0)
        call draw3d(fildes,1.0,0.0,ZC2/9.0)

        call move3d(fildes,0.0,0.0,0.0)
        call draw3d(fildes,0.0,1.0,0.0)
        call draw3d(fildes,0.0,1.0,ZC2/9.0)
        call draw3d(fildes,0.0,0.0,ZC2/9.0)
        call draw3d(fildes,0.0,0.0,0.0)

c ***********************************************************
c       plot x-tics
c ***********************************************************
        call line_color(fildes,R(2),G(2),B(2))
        xmrk=0.0
        do while (xmrk.le.1.0)
        call move3d(fildes,xmrk,0.0,0.0)
        call draw3d(fildes,xmrk,(.02),0.0)
        xmrk=xmrk+(0.25)
        end do
 
c ***********************************************************
c       plot y -tics
c ***********************************************************
        ymrk=0.0
        do while (ymrk.le.1.0)
        call move3d(fildes,0.0,ymrk,0.0)
        call draw3d(fildes,(.02),ymrk,0.0)
        ymrk=ymrk+(0.25)
        end do
 
c ***********************************************************
c       plot the data points
c ***********************************************************
        zeeax=1.0/9.0
        I = jj
        c_cnt = N
        DO 610 J=1,N
        if (I .ne. J) then ! don't display particle w/ itself
        if (mod(j,4) .eq. 0) c_cnt = c_cnt - 1


        m = icolor(mod(J-1,32)+1) + 1
        call line_color(fildes,R(m),G(m),B(m))

        do 600 K=1,NPTS
           if (J .gt. I) then
              i_tmp = J
              j_tmp = I
           else
              i_tmp = I
              j_tmp = J
           end if
           YYY3D(K) = PDIS(i_tmp,j_tmp,K)/YC
           XXX3D(K) = (K*1.0)/(NPTS*1.0)
600     end do
        plotpts=NPTS+1
        plot(1)=0.0
        plot(2)=100.0/YC
        plot(3)=zeeax
        DO 605 K4=0,NPTS-1
          k5=K4*3
          plot(k5+4)=XXX3D(K4+1)
          plot(k5+5)=YYY3D(K4+1)
          plot(k5+6)=zeeax
605     end do
        call polyline3d(fildes,plot,plotpts,FALSE)
 
c ***********************************************************
c       this section labels each plot at its base y and max x 
c ***********************************************************
        call character_height(fildes,.07)
        call text_color(fildes,R(2),G(2),B(2))

        write(unit=strng,fmt='(a1,a1,a1,i2,a1,i2,a1,a1)') 
     +     char(46),char(46),char(40),I,",",J,char(41),char(0)
        call text3d(fildes,1.0,0.0,zeeax,strng,
     +     WORLD_COORDINATE_TEXT,0)

        call line_type(fildes,DOT)
        call line_repeat_length(fildes,.125)
        call move3d(fildes,1.0,0.0,zeeax)
        call draw3d(fildes,1.0,YYY3D(NPTS),zeeax)
        call line_type(fildes,SOLID)
        call line_repeat_length(fildes,.03125)
 
        zeeax=zeeax+1.0/9.0
        end if ! (I .ne. J)
610     CONTINUE
 
        call text_color(fildes,R(2),G(2),B(2))
 
c ***********************************************************
c       do main title
c ***********************************************************
        call text_font_index(fildes, 4)         !sans serif font
        call text_alignment(fildes,TA_CENTER,TA_TOP,0.0,0.0)
        call character_height(fildes,0.09)      !9% of current 0-1 Y range
        write(unit=Strng,fmt='(A)')TITLE
        call text2d(fildes,0.625,0.92,Strng,VDC_TEXT,EndOfLine)

c ***********************************************************
c       do x-axis title
c ***********************************************************
        call text_font_index(fildes,4)
        call text_alignment(fildes,TA_CENTER,TA_BOTTOM,0.0,0.0)
        call text_orientation3d(fildes,0.0,1.0,0.0,1.0,0.0,0.0)
        call character_height(fildes,.08)
        call text3d(fildes,0.5,-.2,
     +  0.0,'Time (sec.)'//char(0),WORLD_COORDINATE_TEXT, EndOfLine)

c ***********************************************************
c       do y-axis title
c ***********************************************************
        call text_orientation3d(fildes,-1.0,0.0,0.0,0.0,1.0,0.0)
        call text3d(fildes,-.2,0.5,0.0,'Distance'//char(0),
     +  WORLD_COORDINATE_TEXT, EndOfLine)
 
c ***********************************************************
c       label x tics
c ***********************************************************
        call clip_indicator(fildes,CLIP_TO_VDC)
        call text_orientation2d(fildes,0.0,1.0,1.0,0.0)
        call character_height(fildes,(.10))
        call text_alignment(fildes,TA_CENTER,TA_TOP,0.0,0.0)
        X=0.0
        do while (X.le.XC)
        write(unit=Strng,fmt='(F5.1,x,a1)')X,char(0)
        call text2d(fildes,X/XC,0.0,Strng,WORLD_COORDINATE_TEXT,
     +  EndOfLine)
        X=X+(0.25*XC)
        end do
 
c ***********************************************************
c       label y tics
c ***********************************************************
        call text_alignment(fildes,TA_RIGHT,TA_HALF,0.0,0.0)
        Y=0.0
        do while (Y.le.YC)
          write(unit=Strng,fmt='(f5.1,x,a1)')Y,char(0)
          call text2d(fildes,0.0,Y/YC,Strng,WORLD_COORDINATE_TEXT,
     +  EndOfLine)
        Y=Y+(.25*YC)
        end do
        call make_picture_current(fildes)   !dump the buffer
 
c ***********************************************************
c       camera parameter adjustment !jaj
c ***********************************************************
171        format(2x,
     +'Adjust camera parameters (0 to exit 3D)?',$)

        if (inst .eq. '1') then !pause for 5 sec betw screens
           call sleep (5)
        else
           print 171
           read(*,'(A)')reply
           if (reply .eq. '0') goto 999
        end if

        if((reply.eq.'y').or.(reply.eq.'Y')) then
 
          call gclose(fildes)
          call mode (0)

          fildes=gopen (1700,820,-700,-5,'xtrydis')
          
          call colors(fildes)
          call view_volume(fildes,-.5,-.5,-.5,1.5,1.5,1.5)
 
          call set_p1_p2(fildes,FRACTIONAL,0.,0.,0.,1.,.9,1.)
          call vdc_extent(fildes,0.0,0.0,0.0,1.25,1.0,ZC/9.0)
          call mapping_mode(fildes,DISTORT)
          call view_camera(fildes,camera)
 
          call move3d(fildes,0.0,0.0,0.0)
          call draw3d(fildes,1.0,0.0,0.0)
          call move3d(fildes,0.0,0.0,0.0)
          call draw3d(fildes,0.0,1.0,0.0)
          call move3d(fildes,0.0,0.0,0.0)
          call draw3d(fildes,0.0,0.0,ZC/9.0)
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
          camera(CAM_REFZ)=zref
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
          call move3d(fildes,0.0,0.0,0.0)
          call draw3d(fildes,1.0,0.0,0.0)
          call move3d(fildes,0.0,0.0,0.0)
          call draw3d(fildes,0.0,1.0,0.0)
          call move3d(fildes,0.0,0.0,0.0)
          call draw3d(fildes,0.0,0.0,ZC/9.0)
          call make_picture_current(fildes)   !dump the buffer
 
          end do
          call mode(3)
          goto 174
        end if
        call mode(3)

        if (inst .ne. '3') jj = jj + 1
        end do !big DO loop

999     call mode(3)
        return
        end
 
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
