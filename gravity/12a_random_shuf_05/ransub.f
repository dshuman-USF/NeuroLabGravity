        subroutine ransub(iseed)
        integer*4 iseed,t(3),sec
c       define seed for random # gen. transparently to user
c
c
        print 755
755      format (2x,'1<RET> to enter iseed; <RET> for auto gen.')
        read (*,960) iopt
        if (iopt.eq.1) then
        print 156
156     format (2x,'Enter iseed (odd integer, 11 digits max: ',$)
        read (*,960) iseed
960     format(I12)

        goto 200
        end if

        call IDATE(t)!get 3 integers
        sec=int(SECNDS(0.0))!seconds since midnight
c       iseed=((sec*month*day*year)*2)+1
        iseed=(((t(1)*t(2)*t(3))/10)*2)+1

200     return
        end

