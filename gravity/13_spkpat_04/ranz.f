c *******************************************************
c        subroutine to insert in spkpat6bg to randomize original
c data set in anary before any analysis
c ********************************************************
c outer loop for pair row definition
c
      subroutine ranz(iseed,lastc,izz,aryin,aryin2)
      parameter (ipairs=120, ncol=999)
      real aryin(ipairs,ncol)
      real aryin2(ipairs,ncol)
      integer*4 iseed
      dimension IRQ(ncol)
c      izz= iarnum in calling program

      do 3396 jv=1,izz
            do new=1,ncol
            IRQ(new)=1
            end do
c
c      inner loop for stepping through the array columns
c
      J=1
8500    R=RAN(iseed)
      ANUM=(R*lastc)+1
      IR=ANUM
8600    IF(IRQ(IR).EQ.0)GOTO 8510
          aryin(jv,J) = aryin2(jv,IR)
      IRQ(IR)=0
      J=J+1
      IF(J.GT.lastc)GOTO 3396
      GOTO 8500
8510    IR=IR+1
      IF(IR.GT.lastc)IR=1
      GOTO 8600
3396      continue

      return
      end

