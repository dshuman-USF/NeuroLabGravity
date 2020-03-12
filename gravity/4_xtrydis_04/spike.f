c       **************************************************************
c       subroutine spike()
c               for displaying spike plot representation with xtrydis
c       **************************************************************
c
        subroutine spike(fildes,XC,TITLE,icolor,num_an)
        INCLUDE 'head_04.def'
        include 'config.defs'

C       place at top before parameters
        parameter (nstep=1000, ncell=64)
        integer*4 EndOfLine,num_an
        real XC,right,left,hrange
        character*30 Strng
        character*26 TITLE
        real spikeplot(nstep+1,ncell)
        DIMENSION icolor(32)

        include 'sbparam.defs'

        parameter (EndOfLine=FALSE)             !used by "text2d"

        right=0.13
        left=0.15
        hrange=1-left-right

        call gclose(fildes)

C       NEXT 3 lines create custom window
C

        fildes=gopen (1700,820,-700,-5,'xtrydis')
        call colors(fildes)
        call vdc_extent(fildes,0.0,0.0,0.0,1.25,1.0,1.0)
        call view_port(fildes,left*1.25,0.1,(1-right)*1.25,.95)
        call view_window(fildes,0.0,-.1,XC,1.0)
 
c   *************************************************
c       do main title
c   *************************************************
        call text_font_index(fildes, 4)         !sans serif font
        call text_alignment(fildes,TA_LEFT,TA_TOP,0.0,0.0)
        call character_height(fildes,0.05)      !4% of current 0-1 Y range
        write(unit=Strng,fmt='(A)')TITLE
        call text2d(fildes,0.5,.99,Strng,VDC_TEXT,
     +  EndOfLine)

c   *************************************************
c       do y-axis title
c   *************************************************
        call character_height(fildes,.06)
        call text_alignment(fildes,TA_CENTER,TA_TOP,0.0,0.0)
        call text_orientation2d(fildes,-1.0,0.0,0.0,1.0)
        call text2d(fildes,1.25*(left-.10),.55,
     + 'Neuron #'//char(0),VDC_TEXT, EndOfLine)
 
c   *************************************************
c       this section draws the time scale on the bottom of the
c       plot.
c   *************************************************
        call move2d(fildes,0.0,-.05)
        call draw2d(fildes,XC,-.05)
 
        call line_color(fildes,R(2),G(2),B(2))
        xmrk=0.0
        do while (xmrk.le.XC)
        call move2d(fildes,xmrk,-.05)
        call draw2d(fildes,xmrk,-.075)
        xmrk=xmrk+(0.25*XC)
        end do
 
c   *************************************************
c       this section labels the time scale tics
c   *************************************************
        size=0.05
        call clip_indicator(fildes,CLIP_TO_VDC)
        call text_orientation2d(fildes,0.0,1.0,1.0,0.0)
        call character_height(fildes,size)
        call text_alignment(fildes,TA_CENTER,TA_TOP,0.0,0.0)
        X=0.0
        do 921 ival=0,4
          write(unit=Strng,fmt='(F5.1,x,a1)')X,char(0)
          xaddr=(.25*ival*hrange+left)*1.25
          call text2d(fildes,xaddr,.1,Strng,VDC_TEXT,
     +    EndOfLine)
          X=X+(0.25*XC)
921     end do
        
c   *************************************************
c       this section titles the time scale
c   *************************************************
        call text_font_index(fildes,4)
        call text_alignment(fildes,TA_CENTER,TA_BOTTOM,0.0,0.0)
        call character_height(fildes,.06)
        call text2d(fildes,1.25*(hrange/2.+left),.003,
     +  'Time (sec.)'//char(0),VDC_TEXT, EndOfLine)
 
c   **************************************************************
c       this section of code plots the spike data on a new screen
c       it gives the option of a dump to laser printer.
c   **************************************************************
        binwidth=XC/float(nstep)
        binwidthms=binwidth*1000.
 
        do 3100 j5=1,N-num_an
          k5=1
          base=0.0
          do 3000 i5=1,LFRAM
            spikes=0.0
         do while((timat(k5,j5).le.(base+binwidthms)).and.(k5.le.nsppc))
              if(timat(k5,j5).gt.0.0) then
                spikes=spikes+1.0
              end if
              k5=k5+1
            end do
            base=base+binwidthms
            spikeplot(i5,j5)=spikes
3000      end do
3100    end do
 
c   *************************************************
c       the next block allows for proportional scaling
c   *************************************************
        yscale=0.0
        do 4100 i6=1,N-num_an
          do 4000 j6=1,LFRAM
            yscale=AMAX1(yscale,spikeplot(j6,i6))
4000      end do
4100    end do
 
        scale=yscale/binwidth
 
        xcoord=((1.0-right)+.02)*1.25
        do 6000 i3=1,N-num_an
          x=0.0 
          base=1.0-1.0*i3/N
          ymax2 = base + (1.0/(N*1.4))
 
          call text_alignment(fildes,TA_RIGHT,TA_BASE,0.,0.)
          call character_height(fildes,0.02) !jaj
          write(unit=strng,fmt='(i2,x,a1)')i3,char(0)
          call wc_to_vdc(fildes,0.0,base,0.0,xlabel,ylabel,zlabel)
          call text2d(fildes,xlabel,ylabel,strng,VDC_TEXT,0)
 
          m = icolor(mod(i3-1,32)+1) + 1
          call perimeter_color(fildes,R(m),G(m),B(m)) 
          call fill_color(fildes,R(m),G(m),B(m)) !jaj1
          do 5500 j3=1,LFRAM
            binval=base+(spikeplot(j3,i3)/(N*1.4*yscale))
            call rectangle(fildes,x,base,x+binwidth,binval)
            x=x+binwidth
5500      end do

c   *************************************************
c         this section creates the scale line for each plot
c   *************************************************
          call line_color(fildes,R(2),G(2),B(2)) !white only
          call move2d(fildes,XC,base)
          call draw2d(fildes,XC,ymax2)
          call draw2d(fildes,XC-XC/100.,ymax2)
 
6000    end do

c   *************************************************
c         this block labels a scale line (above)
c   *************************************************
          call text_alignment(fildes,TA_LEFT,TA_HALF,0.,0.)
          call character_height(fildes,0.05)
          call wc_to_vdc(fildes,XC*1.01,ymax2,0.0,xlabel,ylabel,zlabel)
          write(unit=strng,fmt='(F5.1,x,a1)')scale,char(0)
          call text2d(fildes,xlabel,ylabel,strng,VDC_TEXT,0)
          call text2d(fildes,xlabel,ylabel-.04,'spikes/'
     +     //char(0),VDC_TEXT,0)
          call text2d(fildes,xlabel,ylabel-.08,'second'
     +     //char(0),VDC_TEXT,0)
 
        call make_picture_current(fildes)   !dump the buffer
        return
        end
