      module color_mod
      real :: R(32), G(32), B(32)
      contains

      SUBROUTINE colors32()
c
c       ********************************************************
c       *  This subroutine resets the color table entries      *
c       *  for the first 32 entries in the hardware color      *
c       *  table.    The device corresponding to 'fildes'      *
c       *  must be open when this subroutine is called.        *
c       ********************************************************
c
c
c
c
c
c
C black
      R( 1)= 0.000; G( 1)= 0.000; B( 1)= 0.000
c purple to blue
      R( 2)= 0.500; G( 2)= 0.000; B( 2)= 0.500
      R( 3)= 0.400; G( 3)= 0.000; B( 3)= 0.600
      R( 4)= 0.300; G( 4)= 0.100; B( 4)= 0.700
      R( 5)= 0.200; G( 5)= 0.200; B( 5)= 0.800
      R( 6)= 0.100; G( 6)= 0.300; B( 6)= 0.900
c blue to green
      R( 7)= 0.000; G( 7)= 0.400; B( 7)= 1.000
      R( 8)= 0.000; G( 8)= 0.400; B( 8)= 0.800
      R( 9)= 0.000; G( 9)= 0.400; B( 9)= 0.600
      R(10)= 0.000; G(10)= 0.750; B(10)= 0.700
      R(11)= 0.000; G(11)= 0.450; B(11)= 0.200
      R(12)= 0.000; G(12)= 0.475; B(12)= 0.100
c green to yellow
      R(13)= 0.000; G(13)= 0.500; B(13)= 0.000
      R(14)= 0.333; G(14)= 0.500; B(14)= 0.000
      R(15)= 0.444; G(15)= 0.600; B(15)= 0.000
      R(16)= 0.555; G(16)= 0.700; B(16)= 0.000
      R(17)= 0.666; G(17)= 0.800; B(17)= 0.000
C yellow to orange
      R(18)= 1.000; G(18)= 1.000; B(18)= 0.000
      R(19)= 1.000; G(19)= 0.900; B(19)= 0.000
      R(20)= 1.000; G(20)= 0.800; B(20)= 0.000
      R(21)= 1.000; G(21)= 0.700; B(21)= 0.000
      R(22)= 1.000; G(22)= 0.600; B(22)= 0.000
C orange to red
      R(23)= 1.000; G(23)= 0.500; B(23)= 0.000
      R(24)= 1.000; G(24)= 0.400; B(24)= 0.000
      R(25)= 1.000; G(25)= 0.300; B(25)= 0.000
      R(26)= 1.000; G(26)= 0.200; B(26)= 0.000
      R(27)= 1.000; G(27)= 0.100; B(27)= 0.000
C red to  white
      R(28)= 1.000; G(28)= 0.000; B(28)= 0.000
      R(29)= 1.000; G(29)= 0.250; B(29)= 0.250
      R(30)= 1.000; G(30)= 0.500; B(30)= 0.500
      R(31)= 1.000; G(31)= 0.750; B(31)= 0.750
      R(32)= 1.000; G(32)= 1.000; B(32)= 1.000
c
c
      RETURN
      end subroutine colors32

      SUBROUTINE line_color_index(fildes,i)
      include 'config.defs'
      call line_color (fildes,R(i+1),G(i+1),B(i+1))
      end SUBROUTINE line_color_index

      SUBROUTINE background_color_index(fildes,i)
      include 'config.defs'
      call background_color (fildes,R(i+1),G(i+1),B(i+1))
      end SUBROUTINE background_color_index

      SUBROUTINE text_color_index(fildes,i)
      include 'config.defs'
      call text_color (fildes,R(i+1),G(i+1),B(i+1))
      end SUBROUTINE text_color_index

      end module color_mod
