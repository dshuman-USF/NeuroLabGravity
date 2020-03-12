        SUBROUTINE MODE(IVT)
C
C       HP-UX version
C
C       * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C       * THIS SUBROUTINE WILL ERASE THE SCREEN AND SEND THE CURSOR TO THE  *
C       * HOME POSITION (UPPER LEFT CORNER OF THE SCREEN) ON A VT52 OR A    *
C       * VT100 TERMINAL
C       * ... call with IVT=3 homes cursor and blanks that line on hp300h   * 
C       * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
C
        CHARACTER*1 ESC,BR,TWO,FIVE,EEE,JJJ,HHH
        CHARACTER*1 KKK,MMM,CCC
        PARAMETER (ESC=CHAR(27))        !escape
        PARAMETER (BR=CHAR(91))         ![
        PARAMETER (TWO=CHAR(50))        !2
        PARAMETER (FIVE=CHAR(53))       !5
        PARAMETER (EEE=CHAR(69))        !E
        PARAMETER (HHH=CHAR(72))        !H
        PARAMETER (JJJ=CHAR(74))        !J
        parameter (KKK=CHAR(75))        !K
        parameter (MMM=char(77))        !M
        parameter (CCC=char(67))        !C
        return
        IF(IVT.EQ.1)print 102,ESC,EEE,ESC,HHH           !IVT=1 --> VT52
        IF(IVT.EQ.2)print 103,ESC,BR,TWO,JJJ,ESC,BR,HHH !IVT=2 --> VT100
C       IF(IVT.eq.3)print 104,ESC,KKK,ESC,HHH,ESC,MMM   !IVT=3 --> hp300h
        IF(IVT.eq.3)THEN
          print 104,ESC,KKK,ESC,HHH,ESC,MMM
          print 105
          print 104,ESC,KKK,ESC,HHH,ESC,MMM
        ENDIF
        IF(IVT.EQ.4)then
        print 106,ESC,HHH       !IVT=4 hp300h 
        do i=1,15
         print *,"                                    ",
     +       "               "
        end do
        print 106,ESC,HHH       !IVT=4 hp300h 
        print 105
        print 106,ESC,HHH       !IVT=4 hp300h 
        end if
102     FORMAT(/,2X,4A1)
106     FORMAT(/,2X,2A1)
103     FORMAT(/,2X,7A1)
104     format(/,2x,6a1)
105     FORMAT(120X,$)
        RETURN
        END
