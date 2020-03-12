        INCLUDE 'head_04.def'
        include 'config.defs'

C       place at top before parameters
        EXTERNAL GETLOG_
        parameter(nstep=1000,ncell=64)
        DIMENSION POSIT(ncell,ncell),PRJPOS(ncell,2)
        DIMENSION ORIG(ncell),P1(ncell),P2(ncell)
        DIMENSION V1(ncell),V2(ncell),X(ncell)
        DIMENSION Y(ncell),TEMP1(ncell),TEMP2(ncell)
        dimension icolor(34)
        dimension jmpary(ncell,nstep-1)!ncell, nstep-1
        REAL*4 NDIMOR(ncell)
        REAL*4 MAG,DOT1,DOT2,TEMP3,FLASH(ncell)
        real path(nstep*2,ncell),lastpos(ncell,2)
c       real display(260,32,2),fintens,finsave,speed
        real display(nstep+20,ncell,2),speed ! try this for bg
        real*4 xmax,ymax,xmin,ymin,kbtim
        integer*4 count,texdex,back
        integer*4 npts,upper,lower,step
        character*30 JNAME,ANAME,strng,title
        character*1 inst,resp,dummy,bigger,circ
        data(icolor(i), i=1,34) 
     +  /3,4,5,6,7,8,9,10,11,12,13,14,15,16,
     +  17,18,19,20,21,22,23,24,
     +   25,26,27,28,29,30,31,32,34,35,1,2/
        include 'sbparam.defs'
2       call mode (0)

C       NEXT 3 lines create custom window
C
        fildes=gopen(1024,768,-300,-5,'xprojtm')

        if (fildes.eq.-1) stop
        call colors(fildes)
c       *********************************************
        call mode(3)                                                   
        print 10
10      FORMAT (2x,'Program can NOT display more than 16 neurons')
        print 1122
1122      FORMAT (2x,' INPUT POSITION FILE NAME IS : ',$)
        read (*,'(A)') ANAME
c       clear new jump array - read to change particle color
c       at indicated time steps - times of "jumps" - see xslopebg
c
        do i=1,ncell
        do j= 1,nstep-1
        jmpary(i,j)=0
        end do
        end do
        
        print 11107
        print 11108
11107      FORMAT (2x,' OPTION: JUMP.OUT FILE NAME IS')
11108      FORMAT (2x,' (DEFAULT = jump.out) : ',$)
        read (*,'(A)') JNAME
        if (JNAME.eq.' ') JNAME='jump.out'
        OPEN (UNIT=1,file=ANAME,status='OLD',FORM='UNFORMATTED')
        OPEN (UNIT=2,file=JNAME,status='OLD',
     +  err= 22222,FORM='FORMATTED')                       
11118   read (2,11109,end = 11110)icnmb,istpp
11109   format(I3,I4)
c       fill jump arrary here
        jmpary(icnmb,istpp)= 1
        goto 11118
11110   continue ! DONE reading jump.out file
        close (UNIT= 2)
c
22222   continue ! come here if there is no jump.out file

        call mode(3) ! clear screen
c       determine speed of animation sequence
        print 22207
22207   format (2x,' Choose animation speed:',/,
     +  2x,'   1 - slow',/,
     +  2x,'   2 - standard',/,
     +  2x,'   3 - fast',/,
     +  2x,' ')
        read (*,'(A)') inst
        call mode(3) ! clear screen
        if (inst.eq.' ') speed = 0.025 ! default speed = medium
        if (inst.eq.'1') speed = 0.1 ! slow
        if (inst.eq.'2') speed = 0.025 ! default speed = medium
        if (inst.eq.'3') speed = 0.0025 ! fast
c
        read(1) CBLOCK,IBLOCK,RBLOCK  !read in header info for this plot
c
c
c       check if N > 16 ---> NO GO for this 2-D display of projection
        if (N.gt.16) then
        print 44449
44449   format (2x, 'NO PLOT: Number of neurons > 16... EXIT')
        STOP
        end if
c       *****  set time span for data display  *****************
        xckb=0.0 
        npts=LFRAM
        XCORIG=(ENDTIME+DELTIM)/1000.
        totim=XCORIG
        call mode(3)
        print 11,XCORIG
11      format(2x,'default timespan will be ',f5.1,' sec. or ?: ',$)
        read (*,'(f5.1)')xckb
        if((xckb.gt.0.0).and.(xckb.lt.XCORIG)) then
          npts=int(LFRAM*xckb/XCORIG)
          totim=xckb
        end if
c       *****  set step size for data display  *****************
        step=npts 
        call mode(3)
        print 115,totim
115     format(2x,'animation will stop after every ',f5.1,' sec. or',
     +  ' after every ??? sec. ? ',$)
        read (*,'(f5.1)')kbtim
        if((kbtim.gt.0.).and.(kbtim.lt.totim))step=int(npts*kbtim/totim)
c       ********************************************************
        back=0
c       *****  get plot title  **************************************
        call mode(3)                                                   
        title = ANAME
        call mode(3)                                                   
c       *************************************************************
C       *****fills NDIMOR with zeroes****                                             
        DO 5 J = 1,N
          NDIMOR(J) = 0.
5       CONTINUE
c       *****  initialize parameter maximums and minimums  ***********
        xmin=100000.
        ymin=100000.
        xmax=0.
        ymax=0.
c       **************************************************************
        count=0
C       *****READS INPUT BLOCK*****
200     READ (1,END=900)((POSIT(I,J),J=1,N),I=1,N)
        READ (1,END=900)(FLASH(I),I=1,N)
c
c       
C                                                                               
        IF((N .ge.3).and.(N.le.4)) then 
                isplit = N / 3
C               ***** DEFINE NEW ORIGIN AS CG (N-1)(1)((N)/2) *****                          
                DO 130 J=1,N
                ORIG(J) = 0.0
                do 129 k =1,isplit + 1
                  ORIG(J) = ORIG(J) + POSIT(MOD(k,N)+1,J)/(1. *isplit)
129             continue
130             CONTINUE
C               *****POSIT RE NEW ORIGIN*****
                DO 145 I=1,N
                DO 140 J=1,N
                POSIT(I,J) = POSIT(I,J)-ORIG(J)
140             CONTINUE
145             CONTINUE
C               ***** NOTE ALL POSITIONS ARE NOW RELATIVE TO NEW ORIG*****
C               ***** DEFINE P2 AS CG (N-1)/2 (N) (N-2)/3  *****
                DO 155 J=1,N
                P2(J) = 0.0
                do 149 k =isplit + 1, 2*isplit+ 1
                  P2(J) = P2(J) + POSIT(MOD(k,N)+1,J)/isplit
149             continue
155             CONTINUE
C               ***** NOTE THIS POSITION NOW RELATIVE TO NEW 2D ORIG*****
C               ***** DEFINE P1 AS CG (N)/3. (1)/2. (N-1)/4. *****
                DO 158 J=1,N
                  P1(J) = 0.0
                    do 157 k =2*isplit + 1, 3*isplit + 1
                          P1(J) = P1(J) + POSIT(MOD(k,N)+1,J)/isplit
157                 continue
158             CONTINUE
C               *******************************************************
        ELSE
C
C               ***** DEFINE NEW ORIGIN AS CG (15)(13)
C               (11)(9)(7)(5)(2)(3)((1)/2) *****
                div = 9.
                if (N.lt.15) div = 8.
                if (N.lt.13) div = 7.
                if (N.lt.11) div = 6.
                if (N.lt.9) div = 5.
                if (N.lt.7) div = 4.
                if (N.lt.5) div = 3.
                DO 30 J=1,N
                ORIG(J)=(POSIT(15,J)+
     +          POSIT(13,J)+POSIT(11,J)+POSIT(9,J)+POSIT(5,J)+
     +          POSIT(7,J)+POSIT(2,J)+POSIT(3,J)+(0.5*POSIT(1,J)))/div
30              CONTINUE
C               *******************************************************
C
C               *****POSIT RE NEW ORIGIN*****
                DO 45 I=1,N
                DO 40 J=1,N
                POSIT(I,J) = POSIT(I,J)-ORIG(J)
40              CONTINUE
45              CONTINUE
C               ***** NOTE ALL POSITIONS ARE NOW RELATIVE TO NEW ORIG*****
C
C               *******************************************************
C               ***** DEFINE P2 AS CG (14)
C               (10)(5)(6)(7)((8)/2) *****
                div = 6.
                if (N.lt.14) div = 5.
                if (N.lt.10) div = 4.
                if (N.lt.8) div = 3.
                if (N.lt.7) div = 2.
                DO 55 J=1,N
                P2(J)=
     +          (POSIT(14,J)+POSIT(10,J)+
     +          POSIT(5,J)+POSIT(6,J)+POSIT(7,J)+(0.5*POSIT(8,J)))/div
55              CONTINUE
C               ***** NOTE THIS POSITION NOW RELATIVE TO NEW 2D ORIG*****
C               *******************************************************
C
C               *******************************************************
C               ***** DEFINE P1 AS CG (16)((15)/2)(14)(12)(11)
C               (10)(9)(3)(4)(5) *****
                div = 10.
                if (N.lt.16) div = 9.
                if (N.lt.15) div = 8.
                if (N.lt.14) div = 7.
                if (N.lt.12) div = 6.
                if (N.lt.11) div = 5.
                if (N.lt.10) div = 4.
                if (N.lt.9) div = 3.
                if (N.lt.5) div = 2.
                DO 58 J=1,N
                P1(J)=((0.5*POSIT(15,J))
     +          +POSIT(16,J)+
     +          POSIT(14,J)+POSIT(12,J)+POSIT(11,J)+POSIT(10,J)+
     +          POSIT(9,J)+
     +          POSIT(3,J)+POSIT(4,J)+POSIT(5,J))/div
58              CONTINUE
C               *******************************************************
C
                 END IF
c
C       *****CALCULATE V2 = OR->P2*****
        DO 56 J=1,N
           V2(J) = P2(J)
56      CONTINUE
C       *******************************************************
C       *****CALCULATE V1, NORM X=V1/|V1| *****
        DO 60 J=1,N
        V1(J) = P1(J)
60      CONTINUE
        MAG = 0.0
        DO 62 J=1,N
        MAG = MAG + (V1(J))**2
62      CONTINUE
        DOT2 = MAG              !SAVE IT FOR LATER
        MAG = SQRT(MAG)
        DO 65 J=1,N
        X(J) = V1(J)/MAG
65      CONTINUE
C       *****USE Ax(BxC) = B(A.C) - C(A.B) ;
C       *****CALCULATE Y AND NORMALIZED Y
        DOT1 = 0.0
        DO 70 J=1,N
        DOT1 = DOT1 + V1(J)*V2(J)
70      CONTINUE
        DO 74 J=1,N
        TEMP1(J) = V1(J)*DOT1
74      CONTINUE
        DO 76 J=1,N
        TEMP2(J) = V2(J)*DOT2
76      CONTINUE
        DO 78 J=1,N
        Y(J) = TEMP1(J) - TEMP2(J)
78      CONTINUE
C       *** MAGNITUDE FOR NORMALIZATION ***
        MAG = 0.0
        DO 80 J = 1,N
        MAG = MAG + Y(J)**2
80      CONTINUE
        MAG = SQRT(MAG)
C       *** NORMALIZE Y ***
        DO 82 J=1,N
        Y(J) = Y(J)/MAG
82      CONTINUE
C       *****CALCULATE PROJECTIONS; RESULT IN PRJPOS(I,J)*****
        DO 94 I=1,N
        TEMP3 = 0.0
        DO 90 K=1,N
        TEMP3 = TEMP3 + POSIT(I,K)*X(K)
90      CONTINUE
        PRJPOS(I,1) = TEMP3
        TEMP3 = 0.0
        DO 92 K=1,N
        TEMP3 = TEMP3 + POSIT(I,K)*Y(K)
92      CONTINUE
        PRJPOS(I,2) = TEMP3
94      CONTINUE
c       ******  fill the display array  ******************
        count=count+1
        do 103 k2=1,N
          display(count,k2,1)=PRJPOS(k2,1)
          display(count,k2,2)=PRJPOS(k2,2)
103     end do
c       **************************************************
        goto 200
c       *****  get minimum and maximum x and y values  ***

900     do 1045 i2=1,npts
          do 104 j2=1,N
            xmin=AMIN1(xmin,display(i2,j2,1))
            xmax=AMAX1(xmax,display(i2,j2,1))
            ymin=AMIN1(ymin,display(i2,j2,2))
            ymax=AMAX1(ymax,display(i2,j2,2))
104       end do
1045    end do

c       **************************************************
c       *****  initialize the previous positions  ***
        do 105 i2=1,N
          lastpos(i2,1)=display(1,i2,1)
          lastpos(i2,2)=display(1,i2,2)
105     end do
c       ***************************************************
c       *****  define starbase parameters  ****************
        call vdc_extent(fildes,0.0,0.0,0.0,1.25,1.0,1.0)
        call view_port(fildes,.1*1.25,.1,.9*1.25,.9)

        call view_window(fildes,xmin-10.,ymin-20.,xmax+30.,ymax+10.)
        m = back + 1
        call fill_color(fildes,R(m),G(m),B(m))
        call interior_style(fildes,INT_SOLID,1)
        call rectangle(fildes,xmin-10.,ymin-10.,xmax+10.,ymax+10.)

c       **** show 0.0 time ************************
        call rectangle(fildes,xmax-15.,ymin-10.,xmax+10.,ymin)
        call text_color(fildes,R(2),G(2),B(2))
        write(unit=strng,fmt='(f5.1,a4,x,a1)') 0.0,
     +          " sec",char(0)
        call wc_to_vdc(fildes,xmax-12.0,ymin-7.0,
     +          0.0,wcx,wcy,wcz)
        call character_height(fildes,.03)
        call text2d(fildes,wcx,wcy,strng,VDC_TEXT,0)

        call make_picture_current(fildes)

C       size=.02*(ymax-ymin)
C       sizex=.02*(xmax-xmin)
C       need bigger "particle" circles to hold 2 digits...
        size=.030*(ymax-ymin)
        sizex=.030*(xmax-xmin)

c       *****  plot the data  ******************************
        call interior_style(fildes,INT_HOLLOW,1)

c       +++++  plot initial positions and prompt user to continue  ++++
        call text_alignment(fildes,TA_CENTER,TA_HALF,0.,0.)
        do 106 j=1,N
           m = icolor(mod(j-1,34)+1) + 1
          call fill_color(fildes,R(m),G(m),B(m))
          call perimeter_color(fildes,R(m),G(m),B(m))

          call ellipse(fildes,sizex,size,
     +    display(1,j,1),display(1,j,2),0.)
          call character_height(fildes,.025*(ymax-ymin))
          texdex=icolor(mod(j-1,34)+1)
          m = texdex + 1
          call text_color(fildes,R(m),G(m),B(m))
          write(unit=strng,fmt='(i2,x,a1)')j,char(0)
          call text2d(fildes,display(1,j,1),
     +        display(1,j,2),strng,1,0)
106     end do

c        *****  print option  *******************************
        call make_picture_current(fildes)
        call mode(3)
c        here also the code was removed
          print 997
          read(*,'(a)')dummy
          call mode(3)
c       ****************************************************
c
c       +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

c       ****** erase initial time *********************
          m = back + 1
        call text_color(fildes,R(m),G(m),B(m))
        call text_alignment(fildes,TA_LEFT,TA_NORMAL_VERTICAL,0.,0.)
        write(unit=strng,fmt='(f5.1,a4,x,a1)') 0.0,
     +          " sec",char(0)
        call wc_to_vdc(fildes,xmax-12.0,ymin-7.0,
     +          0.0,wcx,wcy,wcz)
        call character_height(fildes,.03)
        call text2d(fildes,wcx,wcy,strng,VDC_TEXT,0)

        time=xmin-10.
        eltim=0.0
        del=(xmax-xmin+20.)/(1.0*npts)
        count=1
        lower=1
        upper=step

        do while(count.le.npts)
          do 108 i3=lower,upper
C
C       delay loop to "slow" animation on fast computers
C
        time2=0.0
18547   s=secnds(0.0)
        s=secnds(s)
        time2=time2+s
        if(time2.lt.speed) goto 18547
C       end timed delay to slow animation
c
            call character_height(fildes,.028*(ymax-ymin))
            call text_alignment(fildes,TA_CENTER,TA_HALF,0.,0.)

            m = back + 1
            call fill_color(fildes,R(m),G(m),B(m))
            call interior_style(fildes,INT_SOLID,1)
            call rectangle(fildes,xmin-10.,ymin-10.,xmax+10.,ymax+10.)
            call interior_style(fildes,INT_HOLLOW,1)

            do 107 j3=1,N
C             +++++++  erase last position  ++++++++++++++++++
c               m = back + 1
c              call perimeter_color(fildes,R(m),G(m),B(m))
c              call text_color(fildes,R(m),G(m),B(m))
c              write(unit=strng,fmt='(i2,x,a1)')j3,char(0)
c              call text2d(fildes,lastpos(j3,1),lastpos(j3,2)
c     +        ,strng,1,0)
c              call ellipse(fildes,sizex,size,lastpos(j3,1)
c     +        ,lastpos(j3,2),0.)
              texdex=icolor(mod(j3-1,34)+1)
c             ++++++++++++++++++++++++++++++++++++++++++++++++
c             +++++ write current position  +++++++++++++++++
        if (jmpary(j3,count).eq.1) then
           m = icolor(33) + 1
              call fill_color(fildes,R(m),G(m),B(m))
              call perimeter_color(fildes,R(m),G(m),B(m))
              call ellipse(fildes,sizex,size,
     +        display(i3,j3,1),display(i3,j3,2),0.)
              call text_color(fildes,R(m),G(m),B(m))
              write(unit=strng,fmt='(i2,x,a1)')j3,char(0)
              call text2d(fildes,display(i3,j3,1),display(i3,j3,2)
     +        ,strng,1,0)
                else
                   m = icolor(mod(j3-1,34)+1) + 1
              call fill_color(fildes,R(m),G(m),B(m))
              call perimeter_color(fildes,R(m),G(m),B(m))
              call ellipse(fildes,sizex,size,
     +        display(i3,j3,1),display(i3,j3,2),0.)
              m = texdex + 1
              call text_color(fildes,R(m),G(m),B(m))
              write(unit=strng,fmt='(i2,x,a1)')j3,char(0)
              call text2d(fildes,display(i3,j3,1),display(i3,j3,2)
     +        ,strng,1,0)
        end if
c             ++++++++++++++++++++++++++++++++++++++++++++++++
c             +++++  record  neuron paths ++++++++++++++++++++
              path(2*i3-1,j3)=display(i3,j3,1)
              path(2*i3,j3)=display(i3,j3,2)
c             ++++++++++++++++++++++++++++++++++++++++++++++++
c             +++++  update last position array  +++++++++++++
              lastpos(j3,1)=display(i3,j3,1)
              lastpos(j3,2)=display(i3,j3,2)
c             ++++++++++++++++++++++++++++++++++++++++++++++++
107         end do
            count=count+1
            call interior_style(fildes,INT_SOLID,1)
            call fill_color(fildes,R(3),G(3),B(3))
            call perimeter_color(fildes,R(3),G(3),B(3))
            call rectangle(fildes,time,ymin-18.,time+del,ymin-16.)
            time=time+del
            call interior_style(fildes,INT_HOLLOW,1)
            call character_height(fildes,.03)
            call text_alignment(fildes,TA_RIGHT,TA_BOTTOM,0.,0.)

            call text_color(fildes,R(2),G(2),B(2))
            call fill_color(fildes,R(3),G(3),B(3))
            call perimeter_color(fildes,R(2),G(2),B(2))
            call rectangle(fildes,xmax-15.,ymin-10.,xmax+10.,ymin)

            call text_alignment(fildes,TA_LEFT,TA_BOTTOM,0.,0.)
            call character_height(fildes,.03)

            ! ********* erase previous time printout
            m = back + 1
            call text_color(fildes,R(m),G(m),B(m))
            write(unit=strng,fmt='(f5.1,a4,x,a1)') eltim,
     +           " sec",char(0)
            call wc_to_vdc(fildes,xmax-12.0,ymin-7.0,
     +          0.0,wcx,wcy,wcz)
            call text2d(fildes,wcx,wcy,strng,VDC_TEXT,0)

            !********* print new time
            eltim=eltim+totim/npts
            call text_color(fildes,R(2),G(2),B(2))
            write(unit=strng,fmt='(f5.1,a4,x,a1)') eltim,
     +          " sec",char(0)
            call wc_to_vdc(fildes,xmax-12.0,ymin-7.0,
     +          0.0,wcx,wcy,wcz)
            call text2d(fildes,wcx,wcy,strng,VDC_TEXT,0)
            m = back + 1
            call fill_color(fildes,R(m),G(m),B(m))

            call text_alignment(fildes,TA_LEFT,TA_BOTTOM,0.,0.)
          call make_picture_current(fildes)
108       end do

          lower=lower+step
          upper=upper+step
          if(upper.gt.npts) upper=npts

c         *****  print option  *******************************
          call make_picture_current(fildes)
          call mode(3)
            print 997
            read(*,'(a)')dummy
            call mode(3)
c         ****************************************************
        end do
        call mode(3)
        print 95
95      FORMAT (2x,'do you want to see the particle paths ? ',$)
        read (*,'(A)') inst
        call mode(3)
        if((inst.eq.'y').or.(inst.eq.'Y')) then
          do 97 k=1,N
             m = icolor(mod(k-1,34)+1) + 1
            call line_color(fildes,R(m),G(m),B(m))
            call move2d(fildes,path(1,k),path(2,k))
            do 96 j=2,npts
              call draw2d(fildes,path(j*2-1,k),path(j*2,k))
96          end do
97        end do
c         +++++  rewrite last positions for clarity  +++++++++
          call interior_style(fildes,INT_SOLID,1)
          call character_height(fildes,.028*(ymax-ymin))
          call text_alignment(fildes,TA_CENTER,TA_HALF,0.,0.)
          do 98 l=1,N
             m = icolor(mod(l-1,34)+1) + 1
            call perimeter_color(fildes,R(m),G(m),B(m))
            m = back + 1
            call fill_color(fildes,R(m),G(m),B(m))
             m = icolor(mod(l-1,34)+1) + 1
            call text_color(fildes,R(m),G(m),B(m))
            xpansion=(l/10)*1.
            call ellipse(fildes,sizex,size,
     +      display(npts,l,1),display(npts,l,2),0.)
            write(unit=strng,fmt='(i2,x,a1)')l,char(0)
            call text2d(fildes,display(npts,l,1),display(npts,l,2)
     +      ,strng,1,0)
98        end do
c         ++++++++++++++++++++++++++++++++++++++++++++++++++++
c         *****  print option  *******************************
          call make_picture_current(fildes)
          call mode(3)
            print 997
            read(*,'(a)')dummy
            call mode(3)
c         ****************************************************
        end if
        call mode(3)
        print 301
301     FORMAT (2x,'enter c for circle plot, <cr> to continue ? ',$)
        read (*,'(A)') circ
        call mode(3)
        if((circ.eq.'C').or.(circ.eq.'c')) then
c       *****  here is the circle draw version  **************
        call mode(3)
        print 300
300     FORMAT (2x,'do you want to see the particle paths ? ',$)
        read (*,'(A)') resp
        call mode(3)
        bigger = 'n'
        print 303
303     FORMAT (2x,'do you want larger circles (video) ?',$)
        read (*,'(A)') bigger
        call mode(3)
        call gclose(fildes)

C       NEXT 3 lines create custom window
C
        fildes=gopen(1024,768,-300,-5,'xprojtm')
        if (fildes.eq.-1) stop
        call colors(fildes)
c       *****  initialize the previous positions  ***
        do 305 i2=1,N
          lastpos(i2,1)=display(1,i2,1)
          lastpos(i2,2)=display(1,i2,2)
305     end do
c       ***************************************************
c       *****  define starbase parameters  ****************
        call vdc_extent(fildes,0.0,0.0,0.0,1.25,1.0,1.0)
        call view_port(fildes,.1*1.25,.1,.9*1.25,.9)
        call view_window(fildes,xmin-10.,ymin-20.,xmax+30.,ymax+10.)
        call interior_style(fildes,INT_HOLLOW,1)
        call rectangle(fildes,xmin-10.,ymin-10.,xmax+10.,ymax+10.)
        size=.03*(ymax-ymin)
        sizex=.03*(xmax-xmin)
c       ****************************************************
c       *****  write title on screen  **********************
        call character_height(fildes,.06)
        call text_alignment(fildes,TA_CENTER,TA_BOTTOM,0.,0.)
        write(unit=strng,fmt='(A)')title
        call text2d(fildes,.625,.91,strng,VDC_TEXT,0)
c       ****************************************************
c       ***** draw the elapsed time box and sec. label *****
        call rectangle(fildes,xmax-15.,ymin-10.,xmax+10.,ymin)
        call text_alignment(fildes,TA_LEFT,TA_BOTTOM,0.,0.)
        call character_height(fildes,.03)

        write(unit=strng,fmt='(f5.1,a4,x,a1)') 0.0," sec",char(0)
        call wc_to_vdc(fildes,xmax-12.0,ymin-7.0,
     +          0.0,wcx,wcy,wcz)
        call text2d(fildes,wcx,wcy,strng,VDC_TEXT,0)
        call make_picture_current(fildes)

c       ****************************************************
        if((bigger.eq.'y').or.(bigger.eq.'Y')) then
          size = size * 2.0
          sizex = sizex * 2.0
        end if
c       +++++  plot initial positions and prompt user to continue  ++++
        m = back + 1
        call fill_color(fildes,R(m),G(m),B(m))
        call character_height(fildes,1.*size)
        call text_alignment(fildes,TA_CENTER,TA_HALF,0.,0.)
        do 306 j=1,N
           m = icolor(mod(j-1,34)+1) + 1
          call perimeter_color(fildes,R(m),G(m),B(m))
          call ellipse(fildes,sizex,size,path(1,j)
     +    ,path(2,j),0.)
          call text_color(fildes,R(m),G(m),B(m))
          write(unit=strng,fmt='(i2,x,a1)')j,char(0)
          call text2d(fildes,display(1,j,1),display(1,j,2)
     +    ,strng,1,0)
306     end do
c        *****  print option  *******************************
        call make_picture_current(fildes)
        call mode(3)
          print 997
          read(*,'(a)')dummy
          call mode(3)
c       ****************************************************
C       +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        time=xmin-10.
        eltim=0.0
        count=1
        lower=1
        upper=step
        do while(count.le.npts)
          do 310 j=lower,upper
C
C       delay loop to "slow" animation on fast computers
C
        time2=0.0
18548   s=secnds(0.0)
        s=secnds(s)
        time2=time2+s
        if(time2.lt.speed) goto 18548
C       end timed delay to slow animation
            call character_height(fildes,1.*size)
            call text_alignment(fildes,TA_CENTER,TA_HALF,0.,0.)
            do 308 k=1,N
C             +++++++  erase last position  ++++++++++++++++++
              if((resp.ne.'y').and.(resp.ne.'Y')) then
                 m = back + 1
                call perimeter_color(fildes,R(m),G(m),B(m))
              else
                 m = icolor(k) + 1
                call perimeter_color(fildes,R(m),G(m),B(m))
              end if
              m = back + 1
              call text_color(fildes,R(m),G(m),B(m))
              write(unit=strng,fmt='(i2,x,a1)')k,char(0)
              call text2d(fildes,lastpos(k,1),lastpos(k,2)
     +        ,strng,1,0)
              call ellipse(fildes,sizex,size,lastpos(k,1)
     +        ,lastpos(k,2),0.)
c             ++++++++++++++++++++++++++++++++++++++++++++++++
c             +++++  write present position  ++++++++++++++++
              m = icolor(k) + 1
              call perimeter_color(fildes,R(m),G(m),B(m))
              call ellipse(fildes,sizex,size,path(j*2-1,k)
     +        ,path(j*2,k),0.)
              call text_color(fildes,R(m),G(m),B(m))
              write(unit=strng,fmt='(i2,x,a1)')k,char(0)
              call text2d(fildes,display(j,k,1),display(j,k,2)
     +        ,strng,1,0)
c             +++++++++++++++++++++++++++++++++++++++++++++++
c             +++++  update last position array  +++++++++++++
              lastpos(k,1)=path(j*2-1,k)
              lastpos(k,2)=path(j*2,k)
308         end do
            count=count+1
            call interior_style(fildes,INT_SOLID,1)
C           if((bigger.ne.'y').and.(bigger.ne.'Y')) then
              call fill_color(fildes,R(3),G(3),B(3))
              call perimeter_color(fildes,R(3),G(3),B(3))

              call rectangle(fildes,time,ymin-18.,time+del,ymin-16.)

              call text_alignment(fildes,TA_LEFT,TA_BOTTOM,0.,0.)
              call character_height(fildes,.03)

              ! ********* erase previous time printout
              m = back + 1
              call text_color(fildes,R(m),G(m),B(m))
              write(unit=strng,fmt='(f5.1,a4,x,a1)') eltim,
     +          " sec",char(0)
              call wc_to_vdc(fildes,xmax-12.0,ymin-7.0,
     +          0.0,wcx,wcy,wcz)
              call text2d(fildes,wcx,wcy,strng,VDC_TEXT,0)

              !********* print new time
              time=time+del
              eltim=eltim+totim/npts
              call text_color(fildes,R(2),G(2),B(2))
              write(unit=strng,fmt='(f5.1,a4,x,a1)') eltim,
     +          " sec",char(0)
              call wc_to_vdc(fildes,xmax-12.0,ymin-7.0,
     +          0.0,wcx,wcy,wcz)
              call text2d(fildes,wcx,wcy,strng,VDC_TEXT,0)
              m = back + 1
              call fill_color(fildes,R(m),G(m),B(m))
C           end if

            call character_height(fildes,.03)
            call text_alignment(fildes,TA_LEFT,TA_BOTTOM,0.,0.)
            call make_picture_current(fildes)
310       end do
          lower=lower+step
          upper=upper+step
          if(upper.gt.npts) upper=npts
c         *****  print option  *******************************
          call make_picture_current(fildes)
          call mode(3)
            print 997
            read(*,'(a)')dummy
            call mode(3)
c         ****************************************************
        end do
c       ******************************************************
        end if
c       *****  loop for more data option  ********************
        call mode(3)
        print 998
998     format(2x,'c to replot, other char to exit : ',$)
        read(*,'(a)')inst
        call mode(3)
        if((inst.eq.'c').or.(inst.eq.'C')) goto 902
        goto 904
902     CLOSE (UNIT=1)
        call gclose(fildes)
        goto 2
c       ******************************************************
c       *****  print option format  **************************
997     format(2x,'press <cr> to continue ...',$)
904     CLOSE (UNIT=1)
        call gclose(fildes)
        STOP
        END
