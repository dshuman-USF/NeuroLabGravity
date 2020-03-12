c
c       Derived from NUMERICAL RECIPIES
c
        SUBROUTINE INDEXX(N,IARRIN,INDX,icol)
        parameter (npts=1000,ipairs=120)
         DIMENSION IARRIN(ipairs,npts-1),INDX(ipairs)
        DO 11 J=1,N
        INDX(J)=J
11      CONTINUE
        L=N/2+1
        IR =N
10      CONTINUE
        IF(L.GT.1)THEN
          L=L-1
          INDXT=INDX(L)
          IQ=IARRIN(INDXT,icol)
        ELSE
          INDXT=INDX(IR)
          IQ=IARRIN(INDXT,icol)
          INDX(IR)=INDX(1)
          IR=IR-1
          IF(IR.EQ.1)THEN
            INDX(1)=INDXT
            RETURN
          ENDIF
        ENDIF
        I=L
        J=L+L
20      IF(J.LE.IR)THEN
          IF(J.LT.IR)THEN
            IF(IARRIN(INDX(J),icol).LT.IARRIN(INDX(J+1),icol))J=J+1
          ENDIF
          IF(IQ.LT.IARRIN(INDX(J),icol))THEN
            INDX(I)=INDX(J)
            I=J
            J=J+J
          ELSE
            J=IR+1
          ENDIF
        GO TO 20
        ENDIF
        INDX(I)=INDXT
      GO TO 10
      END
