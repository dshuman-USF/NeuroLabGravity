
      INCLUDE 'config.defs'
      INCLUDE 'sbparam.defs'

      integer*4 Edge
      parameter (Edge=TRUE)     !used by "interior style"
      real clist(24)
c     fildes=gopen3d (1279,864,0,0,'xtrydis')
      fildes=gopen3d (600,600,600,0,'test_ogl_sb')
      call background_color (fildes, 1., 1., 1.)
      call clear_view_surface ( fildes )
      if (.true.) then
         call fill_color (fildes, 0., 0., 0.)
         call line_type (fildes, DOT);
         call line_color (fildes, 0., 0., 0.)
         clist = (/
     +        .0,.0,.0,
     +        .2,.1,.0,
     +        .5,.0,.0,
     +        .4,.2,.0,
     +        .5,.5,.0,
     +        .3,.4,.0,
     +        .0,.5,.0,
     +        .1,.3,.0/)

c     +        .0,.0,.0,
c     +        .5,.0,.0,
c     +        .2,.2,.0,
c     +        .5,.5,.0,
c     +        .0,.5,.0/)
         call polyline3d (fildes, clist, 8, 0);
         call make_picture_current(fildes);
      else
         call text_color (fildes, 1., 1., 1.)
         call character_height (fildes, .1)
         call text_font_index (fildes, 6)
c         call view_window (fildes, 0., 0., 1., 1.)
         call text_orientation3d (fildes, 0., 1., 0., 1., 0., 0.)
         call text_alignment (fildes, TA_LEFT, TA_BASE, 0, 0)
         call text3d (fildes, 0.01, .04, .0, "Hello, World",
     +        WORLD_COORDINATE_TEXT, 0)
      end if

      call make_picture_current(fildes)
      print '(/,T10,''Press <cr> to continue  >> '',$)'
      read (*,'(A)')

      
      end
