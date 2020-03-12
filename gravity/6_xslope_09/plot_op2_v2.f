c*************************************************************
c       This subroutine plots pair info for sig and slope calcs
c       ie, *.pos files from gravity with control calcs.
c*************************************************************
c       plot_op2_v1.f - derived from newplot in "'04" series
c       bgl 7/10/2004
        SUBROUTINE plot_op2(fildes,idisary,ixrange,iyrange,title,
     +  N_,loop,iopt,iline,inwary,xc_sec,prsch,
     +  icell1,icell2,il_cnt)
 
        include 'head_04.def'
        include 'sbparam.defs'
        include 'config.defs'

        parameter (is_size=40,ncell=64,npts=1000,ipairs=2016)
        parameter (imxscl=32)

c*************************************************************
c       min and max intensity levels in scaled plots..
c       ncell = max number of cells in data file
c       npts = max num time steps in orig. grav calc.
c       ipairs= max poss. number of unique pairs for ncell
c*************************************************************

        real left,bottom,right,top
        integer*4 xc,yc,istyc,istpyc
        dimension idisary(ipairs,npts-1)
        dimension irow(ipairs)
        dimension icell1(ipairs),icell2(ipairs)
        dimension inwary(ipairs,npts-1)
        character*60 title
        character*1 inst6,prsch

        if(.false.)print *,inwary
        if(.false.)print *,loop
        if(.false.)print *,N_
        if(.false.)print *,prsch

c*************************************************************
c       COLOR SLIDE OPTIONS 
c*************************************************************
 
        ibkgnd=1
        itxtcl=256
        
        if (iopt.eq.1 )call hot()
        
        if (iopt.eq.2 )call graybb()
        
        m = ibkgnd
        call background_color(fildes,R(m),G(m),B(m))
        call clear_control(fildes,CLEAR_DISPLAY_SURFACE)
        call clear_view_surface(fildes)
 
c*************************************************************
C       graph output parameters. These parameters define where on the
C       screen the graphic output will appear.  All have a range of
C       0.0 to 1.0.
C                   Right - the amount of space to the right of graph
C                   Left  - the amount ofspace to the left of graph
C                   Top   - the amount of space above the graph
C                   Botttom - the amount of space below the graph
c*************************************************************
        icont=1 
        xc=ixrange
        yc=MIN0(is_size,iyrange)
        istyc=1
        istpyc=MIN0(yc,il_cnt) !eliminates garbage on last screen of sig aggr
        right=0.13
        left=0.15
        top=0.07
        bottom=0.15
        hrange=1-left-right
        vrange=1-top-bottom
        xstep=hrange/xc
        ystep=vrange/yc
 
        call vdc_extent(fildes,0.0,0.0,0.0,1.25,1.0,0.0)
        call view_port(fildes,left*1.25,bottom,(1-right)*1.25,(1-top))
        call view_window(fildes,0.0,0.0,hrange,vrange)

        m = ibkgnd
        call background_color(fildes,R(m),G(m),B(m))
        call clear_control(fildes,CLEAR_DISPLAY_SURFACE)
        call clear_view_surface(fildes)

c***********************************************************
C       MAIN PLOT
c          ... 1 line/pair starting across the top
C          scaling already done; color given by idisary variables
c          row content defined by value of irow(ipairs) !!!!!
c          ...this pointer defines icell1,icells pair labels!!
c***********************************************************
        left=0.0
        right=left+xstep
        top=vrange
        bottom=top-ystep
 
c***********************************************************
c       mod for line width
c       if iline flag set to one then use gaus 1 smooth..
c***********************************************************
        if (iline.eq.1) call line_width (fildes,(1.25/500),VDC_UNITS)

        do 10000 j=istyc,istpyc
           irow(j)=j            !default values for first run
           do 11000 k=1,ixrange
           if (idisary(irow(j),k).eq.0) then
              call line_color(fildes,
     +        R(ibkgnd),G(ibkgnd),B(ibkgnd)) ! BACKGROUND COLOR
        else
           m = idisary (irow(j),k) + 1
              call line_color(fildes,
     +        R(m),G(m),B(m))
        end if
                call move2d(fildes,left,bottom)
                call draw2d(fildes,left,top)
                left=right
                right=right+xstep
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
        m = itxtcl
        call perimeter_color(fildes,R(m),G(m),B(m))
        call rectangle(fildes,0.0,0.0,hrange,vrange)
        call make_picture_current(fildes)        

c*****************************************************
c       print color scheme + time off to side
c******************************************************
c        if (prsch .eq. 'y') then
c           atop=vrange
c           abottom=top - 0.028
c           do i=32,1,-1
c              m = i + 1
c              call line_color(fildes,R(m),G(m),B(m))
c              call move2d(fildes,0.74,abottom)
c              call draw2d(fildes,0.74,atop)
c              atop=atop - 0.028
c              abottom=abottom - 0.028
c           end do
c           call label_cs2(fildes,vrange,sl_max,itxtcl)
c        end if

        ixp=1
        iyp=1           

        call labels2(fildes,icell1,icell2,irow,ixrange,
     +  yc,ixp,iyp,istyc,istpyc,'a',itxtcl)
        call make_picture_current(fildes)        
        call stitle2(fildes,title,xc_sec,itxtcl)
c       Query user for more screens 
c       ******************************************************
10027   format (2x,'CAPTURE or <CR> to continue..')
           print 10027
           read (*,'(A)') inst6
c**************************************************** 
        call clear_control(fildes,CLEAR_DISPLAY_SURFACE)
        call clear_view_surface(fildes)

        return
        end
 
c*******************************************************
c       SUBROUTINE stitle2
c*******************************************************
        subroutine stitle2(fildes,title,xc_sec,itxtcl)
 
        include 'head_04.def'
        include 'sbparam.defs'
        include 'config.defs'
 
        character*60 title,string
        integer itxtcl
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

c********************************************************
c       output the title
c*********************************************************
        call character_expansion_factor(fildes,.75)
        call text_orientation2d(fildes,0.0,1.0,1.0,0.0)
        call text_alignment(fildes,TA_LEFT,TA_BOTTOM,0.0,0.0)
        call character_height(fildes,.04)
        m = itxtcl
        call text_color(fildes,R(m),G(m),B(m))
        write(unit=string,fmt='(a)') title
        call text2d(fildes,0.0,0.95,string,VDC_TEXT,0)
 
        write(string,2001) xc_sec," sec"
        call text2d(fildes,0.90,.93,string,VDC_TEXT,0)
2001    format (F7.3,a4)

        call make_picture_current(fildes)
        return
        end
 
c*******************************************************
c       SUBROUTINE labels2
c          puts pair labels on each line & adds "row and column
c          pointers when rows are rearranged.
c*******************************************************
        subroutine labels2(fildes,icell1,icell2,irow,ixrange,
     +  irg,ixp,iyp,istyc,istpyc,h_flag,itxtcl)
 
        include 'head_04.def'
        include 'sbparam.defs'
        include 'config.defs'
        parameter (is_size=40,ncell=16,npts=1000,ipairs=120)
        parameter (imxscl=32)
 
        character*9 string
        character*1 h_flag
        integer itxtcl
        integer*4 istyc,istpyc
        real right,left,top,bottom
        dimension icell1(ipairs),icell2(ipairs)
        dimension irow(ipairs)

        if(.false.)print *,iyp
 
        xc=ixrange
        yc=MIN0(is_size,irg)
        right=0.13
        left=0.15
        if (h_flag .eq. 'h') then ! currently not used option
           top=0.10
           bottom=0.03
        else
           top=0.07
           bottom=0.15
        end if
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

c       *********************************************************
c       output selected column tic mark
c       *********************************************************
        if (h_flag .ne. 'h') then
           zleft=0.0-1*xstep
           m = itxtcl
           call line_color(fildes,R(m),G(m),B(m))
           call move2d(fildes,zleft+ixp*xstep,vrange)
c          call draw2d(fildes,zleft+ixp*xstep,vrange+ystep)
           call draw2d(fildes,zleft+ixp*xstep,vrange+1.-topd)
        end if

c       ********************************************************
c       output pair labels for lines
c       *********************************************************
        call character_expansion_factor(fildes,.75)
        call text_orientation2d(fildes,0.0,1.0,1.0,0.0)
        call text_alignment(fildes,TA_LEFT,TA_BOTTOM,0.0,0.0)
        call character_height(fildes,ystep*.7)
        m = itxtcl
        call text_color(fildes,R(m),G(m),B(m))

2000    format(I3,',',I3)
        do 5000, j=istyc,istpyc
          write(string,2000) icell1(irow(j)),icell2(irow(j))
          call text2d(fildes,0.0,bottom,string,VDC_TEXT,0)
          bottom=bottom-ystep
5000    continue
 
          call make_picture_current(fildes)      
        return
        end

c*******************************************************
c       SUBROUTINE label_cs2
c          labels color scheme at right with min and max values
c*******************************************************
        subroutine label_cs2(fildes,vrange,sl_max,itxtcl)
 
        include 'head_04.def'
        include 'sbparam.defs'
        include 'config.defs'

        parameter (imxscl=32)
        character*60 string
        integer itxtcl

        call character_height(fildes,0.028*.7)
        m = itxtcl
        call text_color(fildes,R(m),G(m),B(m))

        tmpmax = sl_max
3000    format(F16.7)
        do j=imxscl,1,-1
          write(string,3000) tmpmax
          call text2d(fildes,1.09,vrange+0.135,string,VDC_TEXT,0)
          vrange=vrange - .028
          tmpmax = tmpmax - sl_max/imxscl
        end do

        call make_picture_current(fildes)        
        return
        end

