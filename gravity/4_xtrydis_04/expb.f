
      INCLUDE 'config.defs'
      INCLUDE 'sbparam.defs'

      integer*4 Edge
      parameter (Edge=TRUE)     !used by "interior style"
      real clist(12)
c     fildes=gopen3d (1279,864,0,0,'xtrydis')
      fildes=gopen (1700,820,-700,-5,'xtrydis')
c      fildes=gopen3d (320,240,0,0,'xtrydis')
      call make_picture_current(fildes);

      if (.true.) then
         call fill_color (fildes, 1., 1., 1.)
         clist = (/0.,0.,0.,.5,0.,0.,.5,.5,0.,0.,.5,0./)
         call rectangle (fildes,0.0,0.0,.5,.5) !frame the viewport
c         call polygon3d (fildes, clist, 4, 0);
         call make_picture_current(fildes);
      end if

      print '(/,T10,''Press <cr> to continue  >> '',$)'
      read (*,'(A)')
      
      end
