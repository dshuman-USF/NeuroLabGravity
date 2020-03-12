        SUBROUTINE distune (IREC,zshft,BINW,FNAME,fildes,ICN,zmean,twsd)
c
c       v4 - handles charge product values; scales for bigger of max histogram or zmean+twsd value
c                               - omits cumsum plot option
c
        DIMENSION zshft(101),ICN(6),
     +  ZHIST(101),ZIH(101)
        CHARACTER*30 FNAME
        real binw,zmean,twsd
c
c
c
c       *****  get starbase aliases and constants  **************************
        include 'config.defs'
        include 'sbparam.defs'
c       *********************************************************************
c
c       *****  set screen parameters  ***************************************
c       call colors(fildes)           ! reset hardware color table
        call vdc_extent(fildes,0.0,0.0,0.0,1.25,1.0,1.0)
        call view_port(fildes,.15*1.25,.15,.85*1.25,.75)
        call view_window(fildes,85.,85.,590.,390.)
c       *********************************************************************
c
c
c       ZHIST ARRAY WILL ALWAYS HAVE original zshft values..
c       zshft values may change - depending upon options.
c
c

        DO 2001 I=1,101
        ZHIST(I)=zshft(I)
2001    continue
c
c       scale before plot
c
        CALL DISSCL(ZHIST,ZIH,0,zero,zmean,twsd)
c
c
c       plot
c
        CALL VECHIS(ZIH,fildes,zero,zmean,twsd)
        call labels (fildes,binw,ICN,IREC,
     +  FNAME) ! label axes, etc
        call title (fildes)
c
c
c       HSTOGRAM HAS BEEN PLOTTED
c
c
c       *******************************************************
c       *******************************************************
c                       OPTIONS
c       *******************************************************
c       ******************************************************
c
        PRINT 2010
2010    FORMAT(36X,'OPTIONS:',//,2X,
     +  '<CR>--NEXT HISTOGRAM')
        READ(*, 2020) INST1
2020    FORMAT(1A2)
        goto 9000       ! leave dishisx
c       *********************************************
c       LEAVE dishisx subroutine
c       **********************************************
9000    CALL CLEAR(fildes)
        RETURN
        END
c       ******************************************************
c       ******************************************************
c
c
c       SOME SUBROUTINES 
c
c
c       *******************************************************
c       *******************************************************
        SUBROUTINE labels(fildes,binw,ICN,IREC,
     +  FNAME)
c
c
c       ****** a subroutine to create x-axis labels  ***********
c
c       *****  get starbase aliases and constants  **************************
        include 'config.defs'
        include 'sbparam.defs'
c       *********************************************************************
c
        character*10 rtic,mtic,ltic,tscale
        character*38 params
        character*30 FNAME
        real scale,binw
        integer*4 ICN(6)
c
c
c
c       ++ draw 3 tics  +++++++++++++++++++++++++++++++++
        call line_color(fildes,0.,1.,1.)    !cyan tics
        call move2d(fildes,85.,85.)
        call draw2d(fildes,85.,90.)
        call move2d(fildes,338.,85.)
        call draw2d(fildes,338.,90.)
        call move2d(fildes,590.,85.)
        call draw2d(fildes,590.,90.)
c       +++++++++++++++++++++++++++++++++++++++++++++++++
c
c       ++ set time scale  ++++++++++++++++++++++++++++++
        scale=1.0
        if(binw.gt.1000.) scale=.001
c       +++++++++++++++++++++++++++++++++++++++++++++++++
c
c       ++ label the tics  ++++++++++++++++++++++++++++++
          write(unit=ltic,fmt='(f8.1,a1)')-50.*scale*binw,char(0)
          mtic='0'//char(0)
          write(unit=rtic,fmt='(f8.1,a1)') 50.*scale*binw,char(0)
        call text_alignment(fildes,TA_CENTER,TA_TOP,0.,0.)
        call character_height(fildes,.04)
        call character_width(fildes,.04*9./16.)
        call text_orientation2d(fildes,0.,1.,1.0,0.0)
        call text_color(fildes,0.,1.,1.)
        call text2d(fildes,.15*1.25,.14,ltic,VDC_TEXT,0)
        call text2d(fildes,.5*1.25,.14,mtic,VDC_TEXT,0)
        call text2d(fildes,.83*1.25,.14,rtic,VDC_TEXT,1)
c
c       previous line changed to 1 in last var from 0.
c       ++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c       ++ label x-axis  +++++++++++++++++++++++++++++++++++
        call character_height(fildes,.045)
        call character_width(fildes,.045*9./16.)
        if(binw.gt.1000.) then
         tscale=' sec'//char(0)
        else
         tscale=' ms'//char(0)
        end if
        call append_text(fildes,tscale,VDC_TEXT,1)
c       ++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c
c       ++ ref, targ, stim labels ++++++++++++++++++++++++++
        call character_height(fildes,.03)
        call character_width(fildes,.03*9./16.)
        call text_alignment(fildes,TA_LEFT,TA_BOTTOM,0.,0.)
        call text2d(fildes,0.0,.04,'ref: '//char(0),VDC_TEXT,1)
        write(unit=params,fmt='(i3,x,a1)')ICN(1),char(0)
        call append_text(fildes,params,VDC_TEXT,1)
        call append_text(fildes,' targ: '//char(0),VDC_TEXT,1)
        write(unit=params,fmt='(i3,x,a1)')ICN(3),char(0)
        call append_text(fildes,params,VDC_TEXT,1)
        call text2d(fildes,0.0,.0,'rec#: '//char(0),VDC_TEXT,1)
        write(unit=params,fmt='(i5,x,a1)')IREC,char(0)
        call append_text(fildes,params,VDC_TEXT,1)
        call make_picture_current(fildes)
c
c       add ons... check
c
c
        call append_text(fildes,'  File: '//char(0),VDC_TEXT,1)
        write(unit=params,fmt='(a30,x,a1)')FNAME,char(0)
        call append_text(fildes,params,VDC_TEXT,1)


c       ++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c
        call make_picture_current(fildes)
        return
        end
c
c
c
c
        SUBROUTINE blockout(fildes)
c
c       ------------------------------------------------------
c       |     THIS IS A SUBROUTINE WHICH DELETES A BLOCK OF  |
C       |     TEXT WRITTEN ON THE TOP OF THE SCREEN.         |
C       ------------------------------------------------------
c
c
c       *****  get starbase aliases and constants  **************************
        include 'config.defs'
        include 'sbparam.defs'
c       *********************************************************************
c
c
c
c       set the window to the text block settings
c       and clear the window.
c
        call vdc_extent(fildes,0.0,0.0,0.0,1.25,1.0,1.0)
        call view_port(fildes,0.0,.76,1.25,1.0)
        call view_window(fildes,0.0,10.0,100.0,0.0)
        call clear_control(fildes,CLEAR_VIEWPORT)
        call clear_view_surface(fildes)
        call make_picture_current(fildes)
c
c       reset the window to the histogram plot settings
c
        call view_port(fildes,.15*1.25,.15,.85*1.25,.75)
        call view_window(fildes,85.,85.,590.,390.)
        return
        end
c
c
        SUBROUTINE blocktext(fildes,string,x,y)
c
c
c       >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c       |  This subroutine writes the string in input parameter |
c       |  "string" to the text block at the position "x","y"   |
c       |  (also input parameters).  The block will hold 10     |
c       |  lines, at about 110 characters per line.  The input  |
c       |  parameters can range al follows:                     |
c       |                                                       |
c       |    "string" - an array of characters up to 110        |
c       |               characters in length.                   |
c       |                                                       |
c       |    "x"      - a real value such that.. 0.0 < x < 100. |
c       |                                                       |
c       |    "y"      - a real value such that.. 0.0 < y < 10.  |
c       |                                                       |
c       <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c
c
c       *****  get starbase aliases and constants  **************************
        include 'config.defs'
        include 'sbparam.defs'
c       *********************************************************************
c
        real x,y
        character*130 string
c
c
c       set the window to the text block settings
c       and write the string to the window
c
        call vdc_extent(fildes,0.0,0.0,0.0,1.25,1.0,1.0)
        call view_port(fildes,0.0,.76,1.25,1.0)
        call view_window(fildes,0.0,10.0,100.0,0.0)
        call text_alignment(fildes,TA_LEFT,TA_BOTTOM,0.,0.)
        call text_orientation2d(fildes,0.,-1.,1.0,0.0)
        call character_height(fildes,1.20)
        call character_width(fildes,.98)
        call text2d(fildes,x,y,string,1,0)
        call make_picture_current(fildes)
c
c       reset the window to the histogram plot settings
        call view_port(fildes,.15*1.25,.15,.85*1.25,.75)
        call view_window(fildes,85.,85.,590.,390.)
        return
        end
c
c
c
c
c
        SUBROUTINE title(fildes)
c
c       This is a routine to allow user to input a title
c       for the plot to be printed on the laser printer.
c       The title will appear at the top of the hardcopy
c       within the text block area and can be erased with
c       the blockout routine.
c
c
c       *****  get starbase aliases and constants  **************************
        include 'config.defs'
        include 'sbparam.defs'
c       *********************************************************************
c
        character*41 string,strtitle
c
c       get plot title from user
c
c        print 100
c100     format(2x,'Graph title (Max. 25 characters)?  : ',$)
c        read(*,'(A)') strtitle
          strtitle='SUM Ref A * Tgt E Charges' 
        write(unit=string,fmt='(A)')strtitle
c
c
c       set the window to the text block settings
c       and write the string to the window
c
        call view_port(fildes,0.0,.76,1.25,1.0)
        call view_window(fildes,0.0,10.0,100.0,0.0)
        call text_alignment(fildes,TA_CENTER,TA_BOTTOM,0.,0.)
        call text_orientation2d(fildes,0.,-1.,1.0,0.0)
        call character_height(fildes,1.28*2.)
        call character_width(fildes,.98*2.)
        call text2d(fildes,50.,8.,string//char(0),1,0)
        call make_picture_current(fildes)
c
c       reset the window to the histogram plot settings
c
        call view_port(fildes,.15*1.25,.15,.85*1.25,.75)
        call view_window(fildes,85.,85.,590.,390.)
        return
        end
c
c       *******************************************************
c
c
        SUBROUTINE VECHIS(ZIH,fildes,zero,zmean,twsd)
c
c       get starbase aliases and constants 
        include 'config.defs'
        include 'sbparam.defs'
        DIMENSION ZIH(101)
c
c       *****  clear the screen  ***
        call clear(fildes)
C       *****  draw axes from 85,390 to 85,90 to 590,90 ***
        call line_type(fildes,SOLID)
        call move2d(fildes,85.,390.)
        call draw2d(fildes,85.,90.)
        call draw2d(fildes,590.,90.)
c         draw zero charge product line - ONLY IF zmean = zero, so no SD option
        if (zero.eq.zmean) call move2d(fildes,85.,(zero+90.))
        if (zero.eq.zmean)call draw2d(fildes,590.,(zero+90.))
c       draw mean of non peak and trough bins
        call move2d(fildes,85.,(zmean+90.))
        call draw2d(fildes,590.,(zmean+90.))
c       draw +/- 2sd of non peak,trough bins
        call line_type(fildes,DOT)
        call move2d(fildes,85.,(zmean+twsd+90.))
        call draw2d(fildes,590.,(zmean+twsd+90.))
        call move2d(fildes,85.,(zmean-twsd+90.))
        call draw2d(fildes,590.,(zmean-twsd+90.))

c       *****  draw a broken line from 338,90 to 338,390  ****
        call line_type(fildes,DOT)
        call move2d(fildes,338.,90.)
        call draw2d(fildes,338.,390.)
        call line_type(fildes,SOLID)
C
        linetype=0
        call dohis(ZIH,fildes,linetype)
        call make_picture_current(fildes)
        RETURN
        END
c
c
c
c       ************************************************** 
c
C       THIS SUBROUTINE DRAWS BARS (POLYGONS) IN HISTOGRAM
        SUBROUTINE dohis(ZIH,fildes,linetype)
c
c       get starbase aliases and constants  
        include 'config.defs'
        include 'sbparam.defs'
        DIMENSION ZIH(101)
        call line_type(fildes,linetype)
        call move2d(fildes,85,90)
        IZ=101
        IX=85
        DO 100 I=1,IZ
        ZIH(I)=ZIH(I)+90.
        J=I+1
        IF(J.EQ.102)GOTO 100
        ZTEMP=ZIH(J)
        ZTEMP=ZTEMP+90.
        IX2=IX+5
        call move2d(fildes,float(IX),ZIH(I))
        call draw2d(fildes,float(IX2),ZIH(I))
        call draw2d(fildes,float(IX2),ZTEMP)
        IX=IX2
100     CONTINUE
        call make_picture_current(fildes)
        RETURN
        END
c
c       ************************************************** 
c
C       ROUTINE DRAWS ONLY OVERLAY IN DIFFERENT linetype
        SUBROUTINE VECOVL(ZIH,fildes)
        include 'config.defs'
        include 'sbparam.defs'
        DIMENSION ZIH(101)
        linetype=2
        call dohis(ZIH,fildes,linetype)
        RETURN
        END
c
c       ************************************************** 
c
C       ROUTINE TO CARRY OUT FINAL SCALING BEFORE PLOT
        SUBROUTINE DISSCL(ZHIST,ZIH,ISTAT,zero,zmean,twsd)
        DIMENSION ZHIST(101),ZIH(101)
        IZ=101
        TOP=-9000000000000.
        CNBOT=+9000000000000.
        BOT=+9000000000000.
        DO 10 I=1,IZ
c       CNTOP=AMAX1(CNTOP,ZHIST(I))
        TOP=AMAX1(TOP,ZHIST(I),zmean+twsd)
        CNBOT=AMIN1(CNBOT,ZHIST(I))
        BOT=AMIN1(BOT,ZHIST(I),zmean-twsd)
10      CONTINUE
c
c
        zero=0.0
        if (CNBOT.lt.0) zero= abs(CNBOT) ! takes "top" of most negative bin as zero
c
c
        zofset=abs(CNBOT)+ 0.025*(TOP-CNBOT)
        if (BOT.lt.CNBOT) zofset= abs(BOT)+ 0.025*(TOP-BOT)
c
c
c
        DO 50 I=1,IZ
        ZHIST(I)=ZHIST(I)+ zofset
50      CONTINUE

        CNTOP=TOP+zofset
        CNBOT=0.0

        IF(CNTOP.EQ.0.)SCLFC=1.0
        IF(CNTOP.EQ.0.)GOTO 30
        IF(ISTAT.EQ.1)GOTO 30
        SCLFC=300.0/CNTOP
30      DO 20 J=1,IZ
        ZVAL=ZHIST(J)*SCLFC
        ZIH(J)=ZVAL
        IF(ZIH(J).GT.300.)GOTO 100
        IF(ZIH(J).LT.0.)GOTO 100
        GOTO 20
100     PRINT 110
110     FORMAT(2X,'ZIH VALUE ERROR')
20      CONTINUE
c       add y axis value = to "0" after shift up.
        zero=zero*SCLFC
        zmean=(zmean+zofset)*SCLFC
        twsd=twsd*SCLFC
        RETURN
        END
c
        SUBROUTINE clear (fildes)
c
c       get starbase aliases and constants 
        include 'config.defs'
        include 'sbparam.defs'
        call clear_control(fildes,CLEAR_DISPLAY_SURFACE)
        call clear_view_surface(fildes)
        call make_picture_current(fildes)
        RETURN
        END
c
