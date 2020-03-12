

        subroutine pick_04(iprim)
c
c       options
c
        integer iprim
150     print 310
200     format(I5)
310     format (5X,'Sort Choices:',
     3/,2x,'1 -- in ORIGINAL Order',
     4/,2x,'2 -- by Interval Duration',
     5/,2x,'3 -- by Analog Channel Hi Peak',
     6/,2x,'4 -- by Analog Channel Low',
     7/,2x,'5 -- by Analog Channel Voltage Range',
     8/,2x,'6 -- by Analog Channel Time to Peak',
     9/,2x,'7 -- by Analog Channel Mean Value',
     9//,4x,'>> ',$)
        READ (*,200) iprim
        if ((iprim.lt.1).or.(iprim.gt.7)) goto 150
        return
        end



