c       another subroutine for gausian smoothing
c
c       derived from gaus4.f
c       ..for smoothing the data used by grvpat10.f, 
c
c
        subroutine smooth(ISMPD,aryidl,ixrange,izz)
c       variable inc not used ...yet
        parameter (npts=1000,ipairs=2016,imxscl=32)
c       npts=max poss. number of steps in gravity plot
c       ipairs = max. poss num ber of pairs in calc
c       imxscl= max no of colors in plot (see newplot rtn.)
        dimension aryidl(ipairs,npts-1)
        dimension temp(ipairs,npts-1)
        dimension GTAB1(19),GTAB2(19),GTAB4(19)
C
        DATA GTAB1 / 0., 0., 0., 0., 0., 0., 0., 5.39910E-02,
     1  2.41971E-01, 3.98942E-01, 2.41971E-01, 5.39910E-02,
     2  0., 0., 0., 0., 0., 0., 0./
C
        DATA GTAB2 / 0., 0., 0., 0., 0., 2.69955E-02, 6.47588E-02,
     1  1.20985E-01, 1.76033E-01, 1.99471E-01, 1.76033E-01,
     2  1.20985E-01, 6.47588E-02, 2.69955E-02, 0., 0., 0., 0., 0./
C
        DATA GTAB4 / 7.9349130E-03, 1.3497742E-02, 2.1569330E-02,
     1  3.2379400E-02, 4.5662273E-02, 6.0492679E-02, 7.5284354E-02,
     2  8.8016331E-02, 9.6667029E-02, 9.9735573E-02, 9.6667029E-02,
     3  8.8016331E-02, 7.5284354E-02, 6.0492679E-02, 4.5662273E-02,
     4  3.2379400E-02, 2.1569330E-02, 1.3497742E-02, 7.9349130E-03/
C
c
        do 130 J=1,izz
        DO 120 IA=1,ixrange
        ICL = IA - 9
        ICR = IA + 9
        IF (ICL.LT.1) ICL = 1
        IF (ICR.GT.ixrange) ICR = ixrange
        temp(J,IA) = 0.0
        DO 110 I=ICL,ICR
        if (ISMPD.eq.1) WFAC = GTAB1(I-IA+10)
        if (ISMPD.eq.2) WFAC = GTAB2(I-IA+10)
        if (ISMPD.eq.4) WFAC = GTAB4(I-IA+10)
        temp(J,IA) = temp(J,IA) + WFAC*aryidl(J,I)
110     CONTINUE
120     CONTINUE
130     continue
c       put smoothed values back into aryidl
        do 170 j=1,izz
        do 160 k=1,ixrange
        aryidl(j,k)= temp(j,k)
160     continue
170     continue
c
c
        RETURN
        END
