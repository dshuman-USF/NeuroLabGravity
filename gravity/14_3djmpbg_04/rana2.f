      subroutine rana2 (ianary,ianary2,ianaflg,iaend,iopta)
      PARAMETER (NA= 200000)
      DIMENSION ianary (NA,2)
      DIMENSION ianary2 (NA,2)
 1    FORMAT (I5,I8)
 105  print 110
 110  FORMAT(2X,'1 - analog or 2 - analog chs. #?;  <cr> = skip: ',$)
      READ (*,'(I2)') iopta
      if (iopta.eq.0) goto 5000
      if ((iopta.eq.1).or.(iopta.eq.2)) then
         ianaflg=1
      else
         goto 105
      end if
c     ************************************************
c     handle 1 OR 2 channel
c     ************************************************
      if ((iopta.eq.1).or.(iopta.eq.2)) then
         rewind (1) 
         PRINT 120
 120     FORMAT(2X,'ENTER ONE channel # (I2) now: ',$)
         READ (*,'(I2)') ICHN1
         ICHNA=ICHN1*4096
         N=1
c     
c     read loop
c     
 130     READ(1,1,end=5000) I,J
         IF (ICHN1.NE.I/4096) GOTO 130
         ianary(N,1)=I-ICHNA
c     
c     next line handles 2's compliment AD conversion
c     
         IF (ianary(N,1).GT.2047) ianary(N,1)=ianary(N,1)-4096
c     next line makes all values  positive
         ianary(N,1)= ianary(N,1) + 2048
         ianary(N,2)=J
         N=N+1
         IF (N.LE.NA) GOTO 130
c     
c     end of loop to handle 1 channel
c     
      end if
c     *****************************************
c     handle second channel only - get here from line 5000
c     *****************************************
 6000 rewind (1) 

      PRINT 122
 122  FORMAT(2X,'ENTER second channel #: ',$)
      READ (*,'(I2)') ICHN2
      ICHNB=ICHN2*4096

      N=1
c     
c     read loop
c     
 133  READ(1,1,end=7000) I,J
      IF (ICHN2.NE.I/4096) GOTO 133
      ianary2(N,1)=I-ICHNB


c     
c     next line handles 2's compliment AD conversion
c     
      IF(ianary2(N,1).GT.2047)
     +     ianary2(N,1)=ianary2(N,1)-4096
c     next line makes all values  positive
      ianary2(N,1)= ianary2(N,1) + 2048
      ianary2(N,2)=J
      N=N+1
      IF (N.LE.NA) GOTO 133
c     
c     end of loop to handle 2nd channel
c     

 5000 if (iopta.eq.2) goto 6000 ! READ 1 DONE, IS THERE a SECOND?
 7000 iaend = N - 1
      return
      END

