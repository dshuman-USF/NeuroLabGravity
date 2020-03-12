c       FIREWORKS v 1 bgl February 1998
c       program to read fireworks files and display times (and thus durations)
c       of all spark sets for a particular criterion. Also plots 
c       ON EACH PLOTTED PAGE with the same time scale 1 or 2
c       analog channels from The EXACTLY coresponding BDT file. 
c
c
c       Spark sets are displayed in some SORTED ORDER. For example fireworks1
c       sorts by locationin time of FIRST SPARK OF SERIES. Other schemes may be found to
c       be more useful with experience - TBA.
c
c       ver2 - add reorder/sort option 2 bgl 3/9/98
c       ************************************************************************        C       place at top before parameters

        parameter (npts=999) ! npts is 1 less than num plotted steps in orig gravity
        parameter (imxscl=32)
        parameter (NA=200000)
        parameter(iperscrn = 30) ! number of rows plotted per screen

c       *************************************************************
c       min and max intensity levels in scaled plots..
c       ncell = max number of cells in data file
c       npts = max num time steps in orig. grav calc.
c       npts=num steps ploted in gravity,nlines=num diff spark sequences plotted here
c       *************************************************************
        real left,bottom,right,top
        real xc,yc,hrange,vrange
        character*60 fwkname,bdtname
        character*2 inst
        dimension anary1(NA,2),anary2(NA,2) 
        dimension ianary1(NA,2),ianary2(NA,2) 
        dimension isel(npts),ny(npts),mx2(npts)
        dimension idxst(npts),idxend(npts)
        integer fwary(npts,npts)
        integer fwaryt(npts)
        include 'sbparam.defs'
        include 'config.defs'

c       
c
c       ******************************************************
c       ******************************************************
c       Fireworks file read in and initial processing
c       ******************************************************
c       ******************************************************
c
        numspkpat=0
        do itt = 1,npts
        idxst(itt)=0
        idxend(itt)=0
        end do
4219    PRINT 4230
4230    FORMAT(2X,'Enter fireworks filename:')
        READ (*,'(A)') fwkname
        if (fwkname.eq.' ') goto 4219 
        OPEN (UNIT=3,file=fwkname,status='OLD',access='SEQUENTIAL',
     +  form='UNFORMATTED')
c
c
c
        do 4235 mf = 1, npts
        read (3,END=4250,ERR=4260,IOSTAT=ios) is,in,
     +  im, (fwaryt(k),k=1,npts)
        if (im.gt.0)then
        isel(mf)=is
        ny(mf)=in
        mx2(mf)=im
        numspkpat=numspkpat+1 ! spark sequences in this one
        ifst=0
        do 2336 jz=1,npts
        fwary(numspkpat,jz)= fwaryt(jz)
        if (ifst.eq.0) then
        if (fwaryt(jz).eq.1) then
        idxst(numspkpat)=jz
        ifst=1
        end if
        end if
2336    continue
        end if
4235    continue
c
c
c
4250    close (UNIT=3)
        numscrn= numspkpat/iperscrn!get num screens (pages to plot)
        irem= numspkpat-(numscrn*iperscrn)
        if (irem.gt.0) then
        numscrn=numscrn+1
        lastscrn=irem
        else
        lastscrn=0
        end if
c       Define ixrange and iyrange for plotting spark part of each displayed page
        ixrange = npts ! parameter for max number of spark positions along x axis
        iyrange = iperscrn ! number of rows plotted per screen

c
c       DEFINE DISPLAY OPTION HERE
c       1= ordered by template column number
c       2= ordered by time of first spark in each spark sequence 
        iopt=2 ! fixed option for display while in development
c        - add user choice later
c
c
c       REORDER SEQUENCE WITH OPTIONAL SORTING ROUTINES HERE

c*****************************************************
c       for this column, reorder rows from hi to low
c*****************************************************
        call idex (numspkpat,idxst,idxend)

c       END OF SORTING ROUTINES
c       ***************************************************
        goto 4300 ! now do analog stuff
c       **********************************************
c       READ ERROR HANDLING
4260    stop'FATAL READ ERROR'
c       *****************************************************
c       *****************************************************
c       ANALOG DATA PROCESSING - FROM BDT FILE..
c       DATA MUST BE EXACTLY SAME AS SOURCE OF DATA USED IN GRAVITY!!
c       TO ENSURE THAT SPARK AND ANALOG TIMES ALIGN EXACTLY.
c
c       clear array for analog data to be plotted on side wall
4300    do 9150 jz=1,2
        do 9155 jy=1,NA
        ianary1(jy,jz)=0
        ianary2(jy,jz)=0
        anary1(jy,jz)=0.0
        anary2(jy,jz)=0.0
9155    continue
9150    continue

c       routine to open BDT file for later read in of analog data
        ianz=0 ! controls analog display option
        PRINT 1030
1030    FORMAT(2X,'BDT File for analog sig, <cr> to skip:')
        READ (*,'(A)') bdtname
        if (bdtname.eq.' ') ianz=1 
c       if (ianz.eq.0) pause
        if (ianz.eq.1) goto 1776
c       
c       NOTE WELL rana2 subroutine expects UNIT = 1!!!!!        
        OPEN(UNIT=1,file=bdtname,FORM='FORMATTED',STATUS='OLD')
C
c       **************************************
c
c       routine to read in analog values 
c       NOTE WELL: clock counts are interger!! have not been
c       multiplied by 0.5 to get msec.
        ianaflg=0
c       pointer to last filled array element
        iaend=0 
        call rana2(ianary1,ianary2,ianaflg,iaend,iopta)
        CLOSE (UNIT=1)
c
c       ***************************************************
c       code to get peak value in current analog signal
c       *****************************************************
c
        if (ianaflg.eq.1) then
        max1 = -32000
        min1 = 10000000
        max2 = -32000
        min2 = 10000000
c       find biggest and smallest analog value for scaling & offset on y axis
        do 62003 iw=1,iaend
        max1 = MAX0(ianary1(iw,1),max1)
        min1 = MIN0(ianary1(iw,1),min1)
        max2 = MAX0(ianary2(iw,1),max2)
        min2 = MIN0(ianary2(iw,1),min2)
62003   continue
c       now have largest int analog value(s)
        end if
c       *****************************************************

c
c
c       SCALE analog array 1 for display on bottom panel. - 
c       half size if TWO signals are to be plotted.
c       ***************************************************
c
        iyr1 = max1-min1
c       will plot sig. 1 on side wall between .01 and .09; side wall goes up to 1.0 in y axis
c       add .05 to all values after scaled.
        yscf1=.08/float(iyr1)
        do 6900 iw=1,iaend
        anary1(iw,1)=((ianary1(iw,1)-min1)*yscf1)+0.01 ! voltage
        anary1(iw,2)=ianary1(iw,2)*.5 ! time in msec.
6900    continue
c       get time span of analog data on bottom panel
        timspan= (anary1(iaend,2) - anary1(1,2))/1000. ! in sec
c
c       debug
        write(*,2001) timspan," sec",iaend
2001    format (F7.3,a4,',','iaend=',I8)
c
c
c       *****************************************
c       analog array 1 now ready to plot !!
c       ********************************************
c       plot second analog signal?
c       *********************************************
        if (iopta.eq.2) then
        iyr2 = max2-min2
c       will plot sig. 2 on side wall between .11 and .19; side wall goes up to 1.0 in y axis
c       add .55 to all values after scaled.
        yscf2=.08/float(iyr2)
        do 6905 iw=1,iaend
        anary2(iw,1) = ((ianary2(iw,1)-min2)*yscf2)+0.11 ! voltage
c       anary2(iw,2) = ianary2(iw,2)*.5 ! time in msec. ! currently do not need time for 2nd ch.
6905    continue
c
        end if
c       ******************************************************
c       ***********************************************************
C       END OF ANALOG PROCESSING OPTION
c       ************************************************************
c       ************************************************************
c       ******************************************************
c       ***********************************************************
c       *************************************************************
c       DISPLAY PLOT LOOP 
c       *************************************************************
c       ******************************************************
c       ***********************************************************

c
c       ************************************
c       NO ANALOG REENTRY POINT & setup colors
c       ************************************
1776    ibkgnd=0 ! background color
        iline=1  ! color of all drawn lines
        itxtcl=1 ! color of all text
c
c        fildes = gopen (1700,420,-700,-5,'xfire')
        fildes = gopen (1700,820,-700,-5,'xfire')
        if (fildes.eq.-1)then
        PRINT 9088
9088    FORMAT(2X,'STARBASE gopen failure..')
        stop
        end if
c
c
c
c       ************************************
c       setup colors
c       ************************************
        call colors (fildes) ! define color table
        call background_color_idx(fildes,ibkgnd)
        call line_color_idx(fildes,iline)
        call text_color_idx(fildes,itxtcl)
c       line thickness 0-thin, 1=fat
        itline=0
c
c
c
        istpyc=0 !initialize this variable here
c
c       *******************************************
c       START OF ACTUAL LOOP HERE:
c       *******************************************
        do 5000 ilp=1, numscrn
c       *************************************************************
C       graph output parameters. These parameters define where on the
C       screen the graphic output will appear.  All have a range of
C       0.0 to 1.0.
C                   Right - the amount of space to the right of graph
C                   Left  - the amount ofspace to the left of graph
C                   Top   - the amount of space above the graph
C                   Botttom - the amount of space below the graph
c       *************************************************************
        xc=ixrange
        yc=iyrange 
        istyc=istpyc+1
        istpyc= ilp*iyrange
        right=0.13
        left=0.15
        top=0.07
        bottom=0.20
        hrange=1-left-right
        vrange=1-top-bottom
        xstep=hrange/xc
        ystep=vrange/yc
        call vdc_extent(fildes,0.0,0.0,0.0,1.25,1.0,0.0)
        call view_port(fildes,left*1.25,bottom,(1-right)*1.25,(1-top))
        call view_window(fildes,0.0,0.0,hrange,vrange)
        call clear_control(fildes,CLEAR_DISPLAY_SURFACE)
        call clear_view_surface(fildes)
        left=0.0
        right=left+xstep
        top=vrange
        bottom=top-ystep
c       ***********************************************************
c       mod for line width
c       ***********************************************************
        if (itline.eq.1) call line_width (fildes,(1.25/500),VDC_UNITS)
        if ((lastscrn.gt.0).and.(ilp.eq.numscrn)) then!last screen? partial?
        istpyc=istyc+lastscrn-1 ! yes, so change end of do loop NOW
        end if
        do 10000 j=istyc,istpyc
           do 11000 k=1,ixrange
c
        if ((iopt.eq.1).and.(fwary(j,k).eq.0)) then 
                left=right
                right=right+xstep
                goto 11000
        end if
        if ((iopt.eq.1).and.(fwary(j,k).eq.1)) then
                call move2d(fildes,left,bottom)
                call draw2d(fildes,left,top)
                left=right
                right=right+xstep
                        goto 11000
                end if

        if ((iopt.eq.2).and.(fwary(idxend(j),k).eq.0)) then 
                left=right
                right=right+xstep
                goto 11000
        end if
        if ((iopt.eq.2).and.(fwary(idxend(j),k).eq.1)) then 
                call move2d(fildes,left,bottom)
                call draw2d(fildes,left,top)
                left=right
                right=right+xstep
                        goto 11000
                end if


11000      continue
           left=0.0
           right=left+xstep
           top=bottom
           bottom=bottom-ystep
10000   continue
          
        call make_picture_current(fildes)        
c       *****************************************************
c         create the frame  
c       ******************************************************
        call interior_style(fildes,INT_HOLLOW,1)
        call rectangle(fildes,0.0,0.0,hrange,vrange)
        call make_picture_current(fildes)        

c       *****************************************************
c         create the pair labels  
c       ******************************************************



        call labels(fildes,isel,ny,mx2,xc,
     +  yc,istyc,istpyc,iopt,idxend)
c       *******************************************
c       plot the analog signal(s) at the bottom
c       get xy coordimnates and scaling to match page..
c       *******************************************
c
c       any analog plot?
c
c       
        if (ianz.eq.1) goto 1812 ! NO ANALOG PLOT


c       *************************************************************
C       graph output parameters. These parameters define where on the
C       screen the graphic output will appear.  All have a range of
C       0.0 to 1.0.
C                   Right - the amount of space to the right of graph
C                   Left  - the amount ofspace to the left of graph
C                   Top   - the amount of space above the graph
C                   Botttom - the amount of space below the graph
c       *************************************************************
c
c                       ANALOG SECTION DISPLAY SETUP
c
        xc=iaend
c       yc=0.2 ! percent of 1.0 y range devoted to analog display
        right=0.13
        left=0.15
        top=0.8
        bottom=0.00
        hrange=1-left-right
        vrange=1.0-top-bottom ! total verticle range raw is 1.0 - analog already scaled
        xstep=hrange/xc
c       ystep=vrange/yc ! vertical stuff already scaled
c 
c
        call vdc_extent(fildes,0.0,0.0,0.0,1.25,1.0,0.0)
        call view_port(fildes,left*1.25,bottom,(1-right)*1.25,(1-top))
        call view_window(fildes,0.0,0.0,hrange,vrange)

c
c
c
        do 468 iw=1,iaend-1   
        zval=iw*xstep
        call move2d(fildes,zval, anary1(iw,1))
        call draw2d(fildes,zval, anary1(iw+1,1))
        if (iopta.eq.2) then ! plot 2nd channel too 
        call move2d(fildes,zval, anary2(iw,1))
        call draw2d(fildes,zval, anary2(iw+1,1))
        end if
468     continue
c         create the frame  
c       ******************************************************
        call interior_style(fildes,INT_HOLLOW,1)
        call rectangle(fildes,0.0,0.0,hrange,vrange)
        call make_picture_current(fildes)        
c
c       ****************************************************
c       come here if no analog plot:
c       ****************************************************
c
1812    call stitle (fildes,fwkname,bdtname,timspan,ianz)
        call make_picture_current(fildes) ! dump buffer to screen
        continue
c
c
c
c       **********************************************************
        PRINT 4285
4285    FORMAT(2X,'<ENTER> for next page') ! for a pause to capture or study screen
        READ (*,'(A)') inst
c       *********************************************************
c       ONE PAGE DONE .... loop to next or done....
5000    continue                
c       ************************************************************
        call gclose(fildes)
        end
c
c
c
c       *******************************************************************
c       SUBROUTINES HERE 
c       *******************************************************************
c       *******************************************************
c       SUBROUTINE stitle
c       *******************************************************
        subroutine stitle(fildes,fwkname,bdtname,timspan,ianz)
 
        include 'sbparam.defs'
        include 'config.defs'
 
        character*60 fwkname,bdtname,string
        real right,left,top,bottom
        right=0.13
        left=0.15
        top=0.07
        bottom=0.07
        hrange=1-left-right
        vrange=1-top-bottom
        call vdc_extent(fildes,0.0,0.0,0.0,1.25,1.0,0.0)
        call view_port(fildes,left*1.25,bottom,(1-right)*1.25,(1-top))
        call view_window(fildes,0.0,0.0,hrange,vrange)

c       ********************************************************
c       output the title
c       *********************************************************
        call character_expansion_factor(fildes,.75)
        call text_orientation2d(fildes,0.0,1.0,1.0,0.0)
        call text_alignment(fildes,TA_LEFT,TA_BOTTOM,0.0,0.0)
        call character_height(fildes,.04)
        write(UNIT=string,fmt='(a)') fwkname
        call text2d(fildes,0.0,0.95,string,VDC_TEXT,0)
c
        if(ianz.ne.1)then
        write(UNIT=string,fmt='(a)') bdtname
        call text2d(fildes,0.45,0.95,string,VDC_TEXT,0)
        end if
c
        write(string,2001) timspan," sec"
        call text2d(fildes,0.90,.93,string,VDC_TEXT,0)
2001    format (F7.3,a4)

        call make_picture_current(fildes)
        return
        end
 
c       *******************************************************
c       SUBROUTINE labels
c          puts pair labels on each line & adds "row and column
c          pointers when rows are rearranged.
c       *******************************************************


        subroutine labels(fildes,isel,ny,mx2,xc,
     +  yc,istyc,istpyc,iopt,idxend)
 
        include 'sbparam.defs'
        include 'config.defs'
        parameter (iperscrn=30,ncell=16,npts=999,ipairs=120)
        parameter (imxscl=32)
 

        character*11 string
        real right,left,top,bottom
        dimension isel(npts),ny(npts),mx2(npts)
        dimension idxend(npts)
        right=0.13
        left=0.15
        top=0.07
        bottom=0.20
        hrange=1-left-right
        vrange=1-top-bottom
        xstep=hrange/xc
        ystep=vrange/yc
 
        call vdc_extent(fildes,0.0,0.0,0.0,1.25,1.0,0.0)
        call view_port(fildes,left*1.25,bottom,(1-right)*1.25,(1-top))
        call view_window(fildes,0.0,0.0,hrange,vrange)
        left=0.00
        right=left+xstep
        topd=1-top
        bottom=topd-ystep



c       ********************************************************
c       output labels for each set of sparks... one set per line
c       *********************************************************
        call character_expansion_factor(fildes,.75)
        call text_orientation2d(fildes,0.0,1.0,1.0,0.0)
        call text_alignment(fildes,TA_LEFT,TA_BOTTOM,0.0,0.0)
        call character_height(fildes,ystep*.7)

2000    format(I2,'-',I3,'  ',I3)
        do 5000, j=istyc,istpyc
          if (iopt.eq.1) then
          write(string,2000) ny(j), isel(j), mx2(j)
        end if
        if (iopt.eq.2) then
        write(string,2000) ny(idxend(j)), isel(idxend(j)),
     +   mx2(idxend(j))
        end if
          call text2d(fildes,0.0,bottom,string,VDC_TEXT,0)
          bottom=bottom-ystep
5000    continue
 
          call make_picture_current(fildes)      
        return
        end


