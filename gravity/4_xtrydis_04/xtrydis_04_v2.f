        program xtrydis _04
c****************************************************************************
C       link with mode.f and colors.f / compile with -ldd300h -lsb1 -lsb2 -lm
c
C       plots *.pos data files of gravity data on hp300h, etc
c
c       v10.0 bgl 19-mar-89
c       v32 3/16/94 bgl handle 32 codes from gbatch32 version of gravity 
c
c       mod so max Y-scale value can be keyed in to match scaling of other
c       plots.
c
c       v. xtrydis_04 6-23-04 bgl
c       increase num neurons to 64 with ncell parameter; some numbers in loops remain hard coded
c*******************************************************************************
 
      INCLUDE 'head_04.def'
      INCLUDE 'config.defs'
      parameter (nstep=1000,ncell=64)
      character*1 inst
      character*90 TITLE
      character*90 Strng
      character*90 posfile
      DIMENSION  XXX(nstep+1),YYY(nstep+1)
      DIMENSION PDIS(ncell,ncell,nstep)
      dimension icolor(32)
      data (icolor(i), i=1,32) /13,3,2,5,6,7,8,9,10,11,12,14,15,16,
     +     17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34/ !plot colors

      INCLUDE 'sbparam.defs'
      integer*4 itmp
      real XCKB,XCORIG,XC,YC,xstp
      real top,bottom,left,right,hrange,vrange !parameters for graph out
      integer*4 Edge,EndOfLine,tick,idex,cnt(80),NPTS
      parameter (Edge=TRUE)     !used by "interior style"
      parameter (EndOfLine=FALSE) !used by "text2d"
      
 2    continue

C       NEXT 3 lines create custom window
C
      fildes=gopen (1279,864,0,0,'xtrydis')
      
      if (fildes.eq.-1)stop
      call colors(fildes)
      
      print 5
 5    format (2x,'..input file? : ',$)
      read (*,'(A)') posfile

 6    format (2x,'..# analog channels? : ',$)
 7    format (I2)
      print 6
      read (*,7) itmp
      
      XCKB=0.0                  !init. var
      OPEN (UNIT=1,file=posfile,status='OLD',form='UNFORMATTED')
      read (1) CBLOCK,IBLOCK,RBLOCK !read in header info for this plot
      read(1) TIMAT             !read in spike times
      
C       ******************************************************
C       READ DATA 
C       ******************************************************
      NPTS=LFRAM                !# of points in file for each pair of cells 
      DO 120 K=1,NPTS
         READ (1,END=130) ((PDIS(I,J,K),J=1,I-1),I=2,N)
 120  CONTINUE
 130  CLOSE (1)
      
      XCORIG= (ENDTIME+DELTIM)/1000. !final state time in sec. of data file
      call mode(3)
      PRINT 22,XCORIG
 22   FORMAT (2X,'default timespan will be ',F5.1,' sec. or ?:',$)
      READ (5,18) XCKB
 18   format (F5.1)
      if (XCKB.eq.0.0)then
         XC=XCORIG
      else 
         XC=XCKB
      end if
      if (XC.lt.XCORIG) NPTS= int(LFRAM*(XC/XCORIG))
      
c       ******************************************************
c       find max y value in data for a given NPTS value
c       ... used for auto scaling of plot
c       ******************************************************
      YC=0.0
      DO 520 I=2,N
         DO 510 J=1,I-1
            DO 505 K=1,NPTS
               YC= AMAX1(YC,PDIS(I,J,K))
 505        end do
 510     CONTINUE
 520  CONTINUE
      YC=YC+0.05*YC
      
      call mode(3)
      print 507
 507  format (2x, 'Enter <CR> to continue or new Y-max.(F5.1): ',$)
      read (5,18) YCCK
      if (YCCK.ne.0) YC=YCCK

c       ******************************************************
C       graph output parameters. These parameters define where on the 
C       screen the graphic output will appear.  All have a range of
C       0.0 to 1.0.
C                   Right - the amount of space to the right of graph
C                   Left  - the amount ofspace to the left of graph
C                   Top   - the amount of space above the graph 
C                   Botttom - the amount of space below the graph
c       ******************************************************
      right=0.13
      left=0.15
      top=0.07
      bottom=0.35
      hrange=1-left-right
      vrange=1-top-bottom

      call vdc_extent(fildes,0.0,0.0,0.0,1.25,1.0,1.0)
      call view_port(fildes,left*1.25,bottom,(1-right)*1.25,(1-top))
      call view_window(fildes,0.0,0.0,XC,YC)
      call interior_style(fildes,INT_HOLLOW,Edge) !hollow rectangles
      call background_color (fildes, 0., 0., 0.)
      call perimeter_color (fildes, 1., 1., 1.)
      call line_color (fildes, 1., 1., 1.)
      call clear_view_surface (fildes)
      call rectangle (fildes,0.0,0.0,XC,YC) !frame the viewport

c       ******************************************************
C       plot the data
c       ******************************************************
      call mode(3)
      TITLE = posfile
      
      DO 220 I=2,N
         DO 210 J=1,I-1
            m = icolor(mod(I-1-1,32)+1) + 1
            call line_color(fildes,R(m),G(m),B(m)) 
            xstp=XC/NPTS
            DO 205 K=1,NPTS
               YYY(K) = PDIS(I,J,K)
               XXX(K)= K*xstp
 205        end do
            call move2d(fildes,0.0,100.)
            do 206 jz=1,NPTS
               call draw2d(fildes,XXX(jz),YYY(jz))
 206        end do
 210     CONTINUE
 220  CONTINUE

c       ******************************************************
c       plot x-tics
c       ******************************************************
      call line_color(fildes,R(2),G(2),B(2))
      xmrk=0.0
      do while (xmrk.le.XC)
         call move2d(fildes,xmrk,0.0)
         call draw2d(fildes,xmrk,(.02*YC))
         xmrk=xmrk+(0.25*XC)
      end do

c       ******************************************************
c       plot y -tics
c       ******************************************************
      ymrk=0.0
      do while (ymrk.le.YC)
         call move2d(fildes,0.0,ymrk)
         call draw2d(fildes,(.02*XC),ymrk)
         ymrk=ymrk+(0.25*YC)
      end do
      
c       ******************************************************
c         do main title
c       ******************************************************

      call text_font_index(fildes, 4) !sans serif font
      call text_alignment(fildes,TA_LEFT,TA_BOTTOM,0.0,0.0)
      call character_height(fildes,0.05) !9% of current 0-1 Y range
      write(unit=Strng,fmt='(A)')TITLE
      call text2d(fildes,0.50,0.93,Strng,VDC_TEXT,
     +     EndOfLine)
      
c       ******************************************************
c       do x-axis title
c       ******************************************************
      call text_font_index(fildes,4)
      call text_alignment(fildes,TA_CENTER,TA_BOTTOM,0.0,0.0)
      call character_height(fildes,.04)
      call text2d(fildes,1.25*(hrange/2.+left),bottom-.14*vrange,
     +     'Time (sec.)'//char(0),VDC_TEXT, EndOfLine)

c       ******************************************************
c       do y-axis title
c       ******************************************************
      call text_alignment(fildes,TA_CENTER,TA_TOP,0.0,0.0)
      call text_orientation2d(fildes,-1.0,0.0,0.0,1.0)
      call text2d(fildes,1.25*(left-.17*hrange),vrange/2.+bottom,
     +     'Distance'//char(0),VDC_TEXT, EndOfLine)

c       ******************************************************
c        label x tics
c       ******************************************************
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
     +        EndOfLine)
         X=X+(0.25*XC)
 221  end do

c       ******************************************************
c       label y tics
c       ******************************************************
      call text_alignment(fildes,TA_RIGHT,TA_HALF,0.0,0.0)
      Y=0.0
      do 222 jval=0,4
         write(unit=Strng,fmt='(f5.1,x,a1)')Y,char(0)
         yaddr=.25*jval*vrange+bottom
         call text2d(fildes,left*1.25,yaddr,Strng,VDC_TEXT,
     +        EndOfLine)
         Y=Y+(.25*YC)
 222  end do
      
      call make_picture_current(fildes) !dump the buffer
      call mode(3)
      
      call text_alignment(fildes,TA_LEFT,TA_BOTTOM,0.,0.)

c      not implemented - these are just setting the default anyway
c      call text_path(fildes,PATH_RIGHT)
c      call text_line_path(fildes,PATH_DOWN)

      call text_orientation2d(fildes,0.0,1.0,1.0,0.0)
      call character_height(fildes,0.012)
      call text_font_index(fildes,1)
      call vdc_extent(fildes,0.,0.,0.,1.25,1.0,1.)
      call view_window(fildes,0.0,0.0,20.0,20.0)
      call view_port(fildes,.12*1.25,0.,.88*1.25,.2)
      
c       ******************************************************
c       spike counter for N neurons
c       ******************************************************
      do ij = 1,N
         cnt(ij) = 0
      end do

      do 3 idex=1,N
         tick=1
         do while((TIMAT(tick,idex).gt.0.0).and.
     +        (tick.le.nsppc))
            tick=tick+1
         end do
         cnt(idex)=tick-1
 3    end do
      
c       ******************************************************
c       particle count block output
c       ******************************************************
      if (N .gt. 1) then
         if (N .lt. 18) then
            imax = N
         else
            imax =18 
         end if

         call text2d(fildes,.05,.250,'particle #    : '//char(0),
     +        VDC_TEXT,1)
         do i=1,imax
            write(unit=strng,fmt='(5x,i2,x,a1)')i,char(0)
            call append_text(fildes,strng,VDC_TEXT,1)
         end do
         call append_text(fildes,'',VDC_TEXT,0)

         call text2d(fildes,.05,.230,'event code    : '//char(0),
     +        VDC_TEXT,1)
         do i=1,imax
            if (idset(i) .eq. 0) then
               write(unit=strng,fmt='(5x,a2,x,a1)')"AN",char(0)
            else
               write(unit=strng,fmt='(4x,i3,x,a1)')IDSET(i),char(0)
            end if
            call append_text(fildes,strng,VDC_TEXT,1)
         end do
         call append_text(fildes,'',VDC_TEXT,0)

         call text2d(fildes,.05,.210,'# of events   : '//char(0),
     +        VDC_TEXT,1)
         if (imax .eq. N) then
            jmax = imax-itmp
         else
            jmax = imax
         end if
         do i=1,jmax
            write(unit=strng,fmt='(3x,i4,x,a1)')cnt(i),char(0)
            call append_text(fildes,strng,VDC_TEXT,1)
         end do
         call append_text(fildes,'',VDC_TEXT,0)
      end if

      if (N .gt. 18) then
         if (N .lt. 36) then
            imax = N
         else
            imax = 36
         end if

         call text2d(fildes,.05,.190,'particle #    : '//char(0),
     +        VDC_TEXT,1)
         do  i=19,imax
            write(unit=strng,fmt='(5x,i2,x,a1)')i,char(0)
            call append_text(fildes,strng,VDC_TEXT,1)
         end do
         call append_text(fildes,'',VDC_TEXT,0)

         call text2d(fildes,.05,.170,'event code    : '//char(0),
     +        VDC_TEXT,1)
         do  i=19,imax
            if (idset(i) .eq. 0) then
               write(unit=strng,fmt='(5x,a2,x,a1)')"AN",char(0)
            else
               write(unit=strng,fmt='(4x,i3,x,a1)')IDSET(i),char(0)
            end if
            call append_text(fildes,strng,VDC_TEXT,1)
         end do
         call append_text(fildes,'',VDC_TEXT,0)

         call text2d(fildes,.05,.150,'# of events   : '//char(0),
     +        VDC_TEXT,1)
         if (imax .eq. N) then
            jmax = imax-itmp
         else
            jmax = imax
         end if
         do  i=19,jmax
            write(unit=strng,fmt='(3x,i4,x,a1)')cnt(i),char(0)
            call append_text(fildes,strng,VDC_TEXT,1)
         end do
         call append_text(fildes,'',VDC_TEXT,0)
      end if

      if (N .gt. 36) then
         if (N .lt. 54) then
            imax = N
         else
            imax = 54
         end if

         call text2d(fildes,.05,.130,'particle #    : '//char(0),
     +        VDC_TEXT,1)
         do  i=37,imax
            write(unit=strng,fmt='(5x,i2,x,a1)')i,char(0)
            call append_text(fildes,strng,VDC_TEXT,1)
         end do
         call append_text(fildes,'',VDC_TEXT,0)

         call text2d(fildes,.05,.110,'event code    : '//char(0),
     +        VDC_TEXT,1)
         do  i=37,imax
            if (idset(i) .eq. 0) then
               write(unit=strng,fmt='(5x,a2,x,a1)')"AN",char(0)
            else
               write(unit=strng,fmt='(4x,i3,x,a1)')IDSET(i),char(0)
            end if
            call append_text(fildes,strng,VDC_TEXT,1)
         end do
         call append_text(fildes,'',VDC_TEXT,0)

         call text2d(fildes,.05,.090,'# of events   : '//char(0),
     +        VDC_TEXT,1)
         if (imax .eq. N) then
            jmax = imax-itmp
         else
            jmax = imax
         end if
         do  i=37,jmax
            write(unit=strng,fmt='(3x,i4,x,a1)')cnt(i),char(0)
            call append_text(fildes,strng,VDC_TEXT,1)
         end do
         call append_text(fildes,'',VDC_TEXT,0)
      end if
      
      if (N .gt. 54) then
         if (N .lt. 72) then
            imax = N
         else
            imax =72 
         end if
         call text2d(fildes,.05,.070,'particle #    : '//char(0),
     +        VDC_TEXT,1)
         do  i=55,imax
            write(unit=strng,fmt='(5x,i2,x,a1)')i,char(0)
            call append_text(fildes,strng,VDC_TEXT,1)
         end do
         call append_text(fildes,'',VDC_TEXT,0)

         call text2d(fildes,.05,.050,'event code    : '//char(0),
     +        VDC_TEXT,1)
         do  i=55,imax
            if (idset(i) .eq. 0) then
               write(unit=strng,fmt='(5x,a2,x,a1)')"AN",char(0)
            else
               write(unit=strng,fmt='(4x,i3,x,a1)')IDSET(i),char(0)
            end if
            call append_text(fildes,strng,VDC_TEXT,1)
         end do
         call append_text(fildes,'',VDC_TEXT,0)

         call text2d(fildes,.05,.030,'# of events   : '//char(0),
     +        VDC_TEXT,1)
         if (imax .eq. N) then
            jmax = imax-itmp
         else
            jmax = imax
         end if
         do  i=55,jmax
            write(unit=strng,fmt='(3x,i4,x,a1)')cnt(i),char(0)
            call append_text(fildes,strng,VDC_TEXT,1)
         end do
         call append_text(fildes,'',VDC_TEXT,0)
      end if
      
      call text2d(fildes,1.10,.975,'PRG: '//char(0),
     +     VDC_TEXT,1)
      write(unit=strng,fmt='(a)')PROGID
      call append_text(fildes,strng,VDC_TEXT,0)
      
      call text2d(fildes,1.10,.950,'OUTPT POS FIL:'
     +     //char(0),VDC_TEXT,0)
      write(unit=strng,fmt='(a)')OFNAME
      call character_height(fildes,0.017)
      call text2d(fildes,1.12,.925,strng
     +     //char(0),VDC_TEXT,0)
      call character_height(fildes,0.020)
      
      call text2d(fildes,1.10,.900,'INPT DATA FIL:'
     +     //char(0),VDC_TEXT,0)
      write(unit=strng,fmt='(a)')IFNAME
      call character_height(fildes,0.017)
      call text2d(fildes,1.12,.875,strng
     +     //char(0),VDC_TEXT,0)
      call character_height(fildes,0.020)
      
      call text2d(fildes,1.10,.850,'MEANINT +/-:'
     +     //char(0),VDC_TEXT,0)
      write(unit=strng,fmt='(i3,x,a1)')IHS,char(0)
      call text2d(fildes,1.12,.825,strng
     +     //char(0),VDC_TEXT,1)
      
      call append_text(fildes,' spikes'
     +     //char(0),VDC_TEXT,0)
      call text2d(fildes,1.10,.800,'ENDTIME:'
     +     //char(0),VDC_TEXT,0)
      write(unit=strng,fmt='(f10.1,x,a1)')ENDTIME,char(0)
      call text2d(fildes,1.12,.775,strng
     +     //char(0),VDC_TEXT,0)
      
      call text2d(fildes,1.10,.750,'FORCE SIGN:'
     +     //char(0),VDC_TEXT,0)
      write(unit=strng,fmt='(f4.1,x,a1)')FS,char(0)
      call text2d(fildes,1.12,.725,strng
     +     //char(0),VDC_TEXT,1)
      
      if(FS.eq.1.0) then
         call append_text(fildes,' excit'
     +        //char(0),VDC_TEXT,0)
      else
         call append_text(fildes,' inhib'
     +        //char(0),VDC_TEXT,0)
      end if
      
      if (ISHIFTA.EQ.0) THEN
         call text2d(fildes,1.10,.700,'ACC DECAY FWD'
     +        //char(0),VDC_TEXT,0)
      else
         call text2d(fildes,1.10,.700,'ACC DECAY BCKWD'
     +        //char(0),VDC_TEXT,0)
      end if
      
      call text2d(fildes,1.10,.675,'NORM: '
     +     //char(0),VDC_TEXT,1)
      write(unit=strng,fmt='(i3,x,a1)')NORM,char(0)
      call append_text(fildes,strng,VDC_TEXT,0)
      
      if (ISHIFTE.EQ.0) THEN
         call text2d(fildes,1.10,.650,'EFF DECAY FWD'
     +        //char(0),VDC_TEXT,0)
      else
         call text2d(fildes,1.10,.650,'EFF DECAY BCKWD'
     +        //char(0),VDC_TEXT,0)
      end if
      
      call text2d(fildes,1.10,.625,'NORM: '
     +     //char(0),VDC_TEXT,1)
      write(unit=strng,fmt='(i3,x,a1)')NORM,char(0)
      call append_text(fildes,strng,VDC_TEXT,0)
      
      call text2d(fildes,1.10,.600,'TM STP:'
     +     //char(0),VDC_TEXT,1)
      write(unit=strng,fmt='(f6.1,x,a1)')DELTIM,char(0)
      call append_text(fildes,strng,VDC_TEXT,0)
      
      call text2d(fildes,1.10,.575,'SLID:'
     +     //char(0),VDC_TEXT,1)
      write(unit=strng,fmt='(E8.2,x,a1)')SLIDE,char(0)
      call append_text(fildes,strng,VDC_TEXT,0)
      
      call text2d(fildes,1.10,.550,'FWDINC:'
     +     //char(0),VDC_TEXT,1)
      write(unit=strng,fmt='(f6.1,x,a1)')FWDINC,char(0)
      call append_text(fildes,strng,VDC_TEXT,0)
      
      call text2d(fildes,1.10,.525,'BAKINC:'
     +     //char(0),VDC_TEXT,1)
      write(unit=strng,fmt='(f6.1,x,a1)')BAKINC,char(0)
      call append_text(fildes,strng,VDC_TEXT,0)
      
      call text2d(fildes,1.10,.500,'FWDTAU:'
     +     //char(0),VDC_TEXT,1)
      write(unit=strng,fmt='(f6.1,x,a1)')FWDTAU,char(0)
      call append_text(fildes,strng,VDC_TEXT,0)
      
      call text2d(fildes,1.10,.475,'BAKTAU:'
     +     //char(0),VDC_TEXT,1)
      write(unit=strng,fmt='(f6.1,x,a1)')BAKTAU,char(0)
      call append_text(fildes,strng,VDC_TEXT,0)
      
      call text2d(fildes,1.10,.450,'CRDI:'
     +     //char(0),VDC_TEXT,1)
      write(unit=strng,fmt='(f6.1,x,a1)')CRDI,char(0)
      call append_text(fildes,strng,VDC_TEXT,0)
      
      write(unit=strng,fmt='(f8.1,x,a1)')stimper,char(0)
      call text2d(fildes,.05,.975,strng,VDC_TEXT,1)
      call append_text(fildes,' ms. between stimuli;'
     +     //char(0),VDC_TEXT,1)
      
      call append_text(fildes,'   every '
     +     //char(0),VDC_TEXT,1)
      write(unit=strng,fmt='(i4,x,a1)')ISTPFR,char(0)
      call append_text(fildes,strng,VDC_TEXT,1)
      call append_text(fildes,'th step saved.     '
     +     //char(0),VDC_TEXT,1)
      
      write(unit=strng,fmt='(i4,x,a1)')LFRAM,char(0)
      call append_text(fildes,strng,VDC_TEXT,1)
      call append_text(fildes,' frames saved.'
     +     //char(0),VDC_TEXT,0)
      call append_text(fildes,'',VDC_TEXT,0)

      call make_picture_current(fildes) !dump the buffer
      
      call mode(3)
      
c       ******************************************************
c       this section gives the 3D plot option
c       ******************************************************
      print 1328
      print 1327
 1328 format (2x,"Enter 3D choice...(hard-copy avail in option 3)")
 1327 format (2x,
     +     '1-auto display; 2-page; 3-single particle; <cr>-continue: '
     +     ,$)
      read (*,'(A)') inst
      call mode(3)
      if (inst.eq.'3'.or.inst.eq.'2'.or.inst.eq.'1') then
         call threed(fildes,NPTS,PDIS,XC,YC,icolor,TITLE,inst)
      end if
      call mode(3)
      
c       ********************************************
c       this section gives the spike plot option 
c       ********************************************
      print 1527
 1527 format (2x,'enter s for spike plot, <cr> to continue..',$)
      read (*,'(A)') inst
      call mode(3)
      if ((inst.eq.'s').or.(inst.eq.'S')) then
         call spike(fildes,XC,TITLE,icolor,itmp)
      end if
      
      print 225
 225  format (1x,'c to replot .. other char to exit',$)
      read (*,230) inst
 230  format (1A1)
      if ((inst.eq.'c').or.(inst.eq.'C')) goto 235
      goto 250
 235  call gclose(fildes)
      goto 2
 250  call gclose(fildes)
      end
      
