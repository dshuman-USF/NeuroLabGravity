      parameter (ncol=999)      ! number of steps over which slopes are calculated
      parameter (ncell=64,npts=1000,ipairs=2016)
      integer iarnum,ihard,ijmp,ismooth,maskop,inc,nomsk,iseed, iseed2, lastc
      dimension imask(ipairs)
      character*256 IFNAM,OFNAM,env
      double precision, allocatable :: aryin (:,:)
      double precision, allocatable :: aryin2 (:,:)
      double precision, allocatable :: scm (:,:)
      real biggest
      double precision, allocatable :: meanobs(:)
      double precision fdr_thresh
c      character*11 outfilespk, outfiletxt,outfilefwk
      integer isigtot, isigcnt(ncol)
      double precision dmin(ncol), pwdist(ncol,ncol)

      double precision, allocatable :: U(:,:), VT(:,:), S(:,:), WORK(:)
      interface
         integer(c_int) function pcg32_boundedrand (bound)
         use iso_c_binding, only: c_int
         integer(c_int), VALUE :: bound
         end function pcg32_boundedrand
      end interface

      call getenv ("svddemo",env)

      fdr_thresh = 0.05
      
      call ransub(iseed)
      iseed2=iseed              ! save first seed for report write-up
      call readx
      
      allocate (aryin2 (iarnum, lastc))
      allocate (meanobs (iarnum))


      n = iarnum
      meanobs = sum (aryin, 2) / lastc
      do i = 1, lastc
         aryin(:,i) = aryin(:,i) - meanobs
      end do
      
      if (env.eq.'shuffled') then
         call rotate_rows
         aryin = aryin2
         call print_data
         stop
      end if

      if (env.eq.'correlated') then
         call print_data
         stop
      end if

      allocate (scm (n, n))
      scm = matmul (aryin, transpose (aryin)) / (n - 1)

      allocate (U (n, n))
      allocate (VT (n, n))
      allocate (S (n, n))
      allocate (WORK (10 * n))
      call dgesvd_ ('A','A',n,n,scm,n,S,U,n,VT,n,WORK,10*n,info)
      aryin = matmul (VT, aryin)
      
      if (env.eq.'rotshuf') then
         call rotate_rows
         aryin = aryin2
         call print_data
         stop
      end if
      
      if (env.eq.'rotback') then
         call rotate_rows
         aryin = aryin2
         aryin = matmul (U, aryin)
         call print_data
         stop
      end if
      
      if (env.eq.'rotated') then
         call print_data
         stop
      end if
      
      call pairwise_distances
      call significant_distances
      print *, 'isigtot: ', isigtot
      print *, 'dmin(1): ', dmin(1)
      do i = 2, lastc
         print *, 'isigcnt	', isigcnt(i)
      end do
      do i = 2, lastc
         print *, 'pwdist	', pwdist(1,i)
      end do

      
      CONTAINS

c     ******************************************************************
      SUBROUTINE PRINT_DATA
c     ******************************************************************
      do ir = 1, iarnum
         do ic = 1, lastc
            print "(f21.17,$)", aryin(ir,ic)
         end do
         print *
      end do
      end subroutine print_data
      

c     ******************************************************************
      SUBROUTINE ROTATE_ROWS
c     ******************************************************************
      integer irow, icol, ic0
      do irow = 1, iarnum
         if((maskop.eq.1).and.(imask(k).eq.1)) cycle
         ic0 = pcg32_boundedrand (lastc)
         do icol = 1, lastc
            aryin2(irow, modulo (ic0 + icol, lastc) + 1) = aryin(irow, icol)
         end do
      end do
      end subroutine rotate_rows

c     ******************************************************************
      SUBROUTINE PAIRWISE_DISTANCES
c     ******************************************************************
      integer iref, itar
      double precision d
      do iref = 1, lastc
         pwdist(iref, iref) = 0
         do itar = iref + 1, lastc
            d = distance (aryin, iref, itar)
            pwdist (iref, itar) = d
            pwdist (itar, iref) = d
         end do
      end do
      end subroutine pairwise_distances

c     ******************************************************************
      SUBROUTINE UPDATE_SIGCNT (iref)
c     ******************************************************************
      isigcnt(iref) = 0
      do itar = 1, lastc
         if (itar.eq.iref) cycle
         if (pwdist (iref, itar).lt.dmin(iref)) then
            isigcnt(iref) = isigcnt(iref) + 1
         end if
      end do
      end subroutine update_sigcnt

c     ******************************************************************
      FUNCTION DISTANCE(ary, iref, itar)
c     ******************************************************************
      double precision distance, sum
      double precision ary (:,:)
      sum = 0
      do irow = 1, iarnum
         if((maskop.eq.1).and.(imask(k).eq.1)) cycle
         sum = sum + (aryin(irow,iref) - ary(irow,itar))**2
      end do
      distance = sum
      end function distance

c     ******************************************************************
      SUBROUTINE SIGNIFICANT_DISTANCES
c     ******************************************************************
      double precision d, fdr
      integer idata, itar, iref, ireps
      idata = lastc * (lastc - 1)
      isigtot = idata
      fdr = 1
      ireps = 0
      isigcnt = lastc - 1
      dmin = huge(dmin(1))
      do while (isigtot.gt.0.and.fdr.gt.fdr_thresh)
         call rotate_rows
         do itar = 1, lastc
            isigtot = 0
            do iref = 1, lastc
               if (isigcnt(iref).eq.0) cycle
               d = distance (aryin2, iref, itar)
               if (iref.eq.1) print *,'surdist	', d
               if (d.lt.dmin(iref)) then
                  dmin(iref) = d
                  call update_sigcnt (iref)
               end if
               isigtot = isigtot + isigcnt(iref)
            end do
            ireps = ireps + 1
            fdr = idata / (ireps + 1.) / isigtot
            if (fdr <= fdr_thresh) exit
            print *, "fdr: ", fdr, "isigtot: ", isigtot
         end do
      end do
      end subroutine significant_distances

c     ******************************************************************
      SUBROUTINE RANSUB(iseed)
c     ******************************************************************
      integer*4 iseed,t(3),sec, iopt
c     define seed for random # gen. transparently to user
      print "(2x,'1<RET> to enter iseed; <RET> for auto gen.')"
      read (*,"(I12)") iopt
      if (iopt.eq.1) then
         print "(2x,'Enter iseed (integer, 11 digits max: ',$)"
         read (*,"(I12)") iseed
      else
         call IDATE(t)          !get 3 integers
         sec=int(SECNDS(0.0))   !seconds since midnight
         iseed=(((t(1)*t(2)*t(3))/10)*2)+1
      end if
      call pcg32_srandom(iseed, 54);
      end subroutine ransub

c     ******************************************************************
      SUBROUTINE READX
c     ******************************************************************
      double precision sclfac,tmpmax,tmpka,thrnul
      real temp (ipairs,ncol)
      integer iz, ka
c     imask = array with row (pair) numbers to use; others skipped
c     maskop = flag for using mask option; 1=use mask option
c     nomsk= value to replace variable 'lastr' when mask option engaged
c     ...but only as denominator in calcs - not as a loop controller
c     ....nomsk is determined here in this routine-see below!!
     
      thrnul= .03125
      print "(2x,'PROGRAM EXPECTS npts = 1000 in *.pos input files')"

      print "(2x,'Input data file name:',$)"
      read (*,'(A)') IFNAM
      print "(2x,'Output data file name:',$)"
      read (*,'(A)') OFNAM
      print "(2x,'Enter 1 for MASK option; def is NO MASK')"
      read (*,"(I3)") maskop
      if (maskop.ne.1) maskop = 0

      OPEN (unit=1,file=IFNAM,status='OLD',access='sequential', form='UNFORMATTED')
      read (1) ismooth
      read (1) inc
      read (1) ihard
      read (1) ijmp
      read (1) iarnum
      read (1) biggest
      lastc=ncol/inc

      nomsk=0
      do i=1,ipairs
         read (1) imask(i)
         if (imask(i).eq.0) nomsk=nomsk+1
      end do
      write (*,"(2x,'value of nomsk variable = ',i8)") nomsk

      do i=1,iarnum
         do j=1,ncol
            read (1) temp(i,j)
         end do
      end do
c     NZ is the number of neurons represented in input
c     ..all neurons are counted irrespective of mask option status!!
      NZ = (int (sqrt (8 * iarnum + 1.)) - 1) / 2;

      write (*,"(2x,'# of different spike codes: ',I7)")NZ
      close (unit=1)

c     now define global scaling factor
c     ... from value of variable "biggest" = largest distance between
c     ...pairs of particles in all 20 or 100 control cycles.
     
      sclfac=100./biggest
     
c     scale here with same value as later for all 100 sh*.pos files!!!
      do i=1,iarnum
         do j=1,ncol
            temp(i,j)=real(sclfac*temp(i,j))
         end do
      end do
     
c     scale slopes to be used to max used slope
c     LOOKS AT ALL COLUMNS, ignores inc variable
c     max slope looked for has to be for a neuron pair
c     that is NOT MASKED OUT!!
     
      tmpmax=0.0
      do i=1,iarnum
         if ((maskop.eq.1).and.(imask(i).eq.1)) cycle
         do j=1,ncol
            tmpmax=AMAX1(tmpmax,temp(i,j))
         end do
      end do

c     now fill array to be passed back to main
c     ..if inc>1 then get biggest value in range - only makes a diff.
c     .... if smoothing in effect
c     NOTE THAT number of columns is now defined by var iz

      print *, "tmpmax: ", tmpmax

      allocate (aryin  (2, lastc))
      
      k = 0
      do i=1,iarnum
         if(i.ne.116.and.i.ne.119) cycle
         k=k+1
         iz=0
         do j=1,ncol,inc
            iz=iz+1
            if (iz.gt.lastc) exit
            aryin(k,iz)=0.0
            do ka=j,j+inc-1
               tmpka=temp(i,j)
               aryin(k,iz)=AMAX1(tmpka,aryin(k,iz))
            end do

         end do
      end do
      iarnum = 2
      
      end subroutine readx
      end program
      
