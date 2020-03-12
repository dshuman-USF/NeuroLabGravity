      module surdat_mod
      
      private
      save
      double precision, allocatable :: X0 (:,:), X(:,:)
      double precision, allocatable :: scm(:,:), E(:,:), Er(:,:)
      double precision, allocatable :: V(:), W(:), meanobs(:)
      integer, allocatable :: shift(:)
      integer :: nrows = 0, ncols, n
      public surdat, surdat_init
      
      interface
         integer(c_int) function pcg32_boundedrand (bound)
         use iso_c_binding, only: c_int
         integer(c_int), VALUE :: bound
         end function pcg32_boundedrand
      end interface
      
      contains
      
c     ******************************************************************
      SUBROUTINE SURDAT_INIT (aryin)
c     ******************************************************************
      double precision aryin (:,:)
      if (size(aryin,1).ne.nrows.or.size(aryin,2).ne.ncols) then
         if (allocated (X)) then
            deallocate (X0,X,scm,E,Er,V,W,meanobs,shift)
         end if
         nrows = size(aryin,1)
         ncols = size(aryin,2)
         n = nrows
         allocate (X(nrows, ncols), meanobs(nrows), scm(n,n), E(n,n))
         allocate (X0(nrows, ncols), V(n), shift(nrows), W(n), Er (n,n))
      end if

c     center the data
      meanobs = sum (aryin, 2) / ncols
      do i = 1, ncols
         X0(:,i) = aryin(:,i) - meanobs
      end do

c     1. calculate the SCM (sample covariance matrix)
      scm = matmul (X0, transpose (X0)) / (ncols - 1)

c     2. calculate the eigendecomposition of the SCM
      call evd (scm, E, V)

      end subroutine surdat_init
      
c     ******************************************************************
      SUBROUTINE SURDAT (Y)
c     ******************************************************************
      double precision Y (:,:)
c     3. randomly rotate the rows of the centered data
      do irow = 1, nrows
         shift(irow) = pcg32_boundedrand (ncols)
      end do
      X = cshift (X0, shift, 2)

c     4. calculate the SCM of the rotated data
      scm = matmul (X, transpose (X)) / (ncols - 1)

c     5. calculate the eigendecomposition of the SCM of the rotated data
      call evd (scm, Er, W)
      
c     6. Premultiply Xᵣ by Eᵣ' to get uncorrelated data Xᵤ with an SCM of Λᵣ (= W)
      X = matmul (transpose(Er), X)
      
c     7. scale each row of Xᵤ to make the variances V
      nzero = nrows - ncols + 1
      if (nzero < 0) nzero = 0
      do i = nzero + 1, nrows
         X(i,:) = X(i,:) * sqrt (V(i) / W(i))
      end do

c     8. Premultiply the rotated scaled data by E to get data with the original SCM
      Y = matmul (E, X)
      
c     put the mean back
      do i = 1, ncols
         Y(:,i) = Y(:,i) + meanobs
      end do
      end subroutine surdat

c     ******************************************************************
      SUBROUTINE EVD(scm,E,V)
c     ******************************************************************
      double precision :: scm(:,:), E(:,:), V(:)
      integer :: n = 0, LWORK, LIWORK, info
      double precision, allocatable, save :: WORK(:)
      integer, allocatable, save :: IWORK(:)

      if (n.ne.size(scm,1)) then
         if (allocated(WORK)) deallocate (WORK,IWORK)
         n = size (scm, 1)
         LWORK = 1 + 6*n + 2*n**2
         allocate (WORK (LWORK))
         LIWORK = 3 + 5*n
         allocate (IWORK (LIWORK))
      end if

      call dsyevd_ ('V', 'U', n, scm, n, V, WORK, LWORK, IWORK, LIWORK, info)
      if (info.ne.0) then
         print *,"info = ",info
         stop
      end if
      E = scm                   ! dsyevd overwrites scm with the eigenvectors.  Return them as E.
      end subroutine evd

      end module surdat_mod
      
      program spkpatdist
      use surdat_mod
      parameter (ncol=999)      ! number of steps over which slopes are calculated
      parameter (ncell=64,npts=1000,ipairs=2016)
      integer iarnum,ihard,ijmp,ismooth,maskop,inc,nomsk,iseed, iseed2, lastc
      dimension imask(ipairs)
      character*256 IFNAM,OFNAM,env
      double precision, allocatable :: aryin (:,:), Y(:,:)
      double precision, allocatable :: aryin2 (:,:)
      real biggest
      double precision fdr_thresh
c      character*11 outfilespk, outfiletxt,outfilefwk
      integer isigtot, isigcnt(ncol)
      double precision dmin(ncol), pwdist(ncol,ncol)

      
      interface
         integer(c_int) function pcg32_boundedrand (bound)
         use iso_c_binding, only: c_int
         integer(c_int), VALUE :: bound
         end function pcg32_boundedrand
      end interface

      call getenv ("spkpatdist",env)

      fdr_thresh = 0.05
      
      call ransub(iseed)
      iseed2=iseed              ! save first seed for report write-up
      call readx

      if (env.eq.'aryin') then
         call print_mat (aryin)
         stop
      end if
      if (env(1:6).eq.'surdat') then
         allocate (Y(size(aryin,1), size(aryin,2)))
         call surdat_init (aryin)
         call surdat (Y)
         call print_mat (Y)
         stop
      end if

      allocate (aryin2 (iarnum, lastc))
      
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
      SUBROUTINE PRINT_MAT(MAT)
c     ******************************************************************
      double precision mat(:,:)
      integer nrows, ncols
      nrows = size(mat,1)
      ncols = size(mat,2)
      do ir = 1, nrows
         do ic = 1, ncols
            print "(e25.17,$)", mat(ir,ic)
         end do
         print *
      end do
      end subroutine print_mat

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
      call surdat_init (aryin)
      nzero = 0
      do while (isigtot.gt.0.and.fdr.gt.fdr_thresh)
         call surdat (aryin2)


      nzerotot = 0
      if (inc.eq.1) then
         do j=1,ncol,inc
            do i=1,iarnum
               if (aryin2(i, j).eq.0) nzerotot = nzerotot + 1
            end do
         end do
      end if
      print *, "nzerotot: ", nzerotot
         
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
            if (sum(aryin2(:,itar)).eq.0) nzero = nzero + 1
            print *, "fdr: ", fdr, "isigtot: ", isigtot, "ireps: ", ireps, "idata: ", idata, "nzero: ", nzero
            if (fdr <= fdr_thresh) exit
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
      integer iz, ka, iskip
      character*10 cskip
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
      print *, "inc: ", inc
      print *, "iarnum: ", iarnum

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
      NZ = (int (sqrt (8 * iarnum + 1.)) + 1) / 2;


      itot=1
      do 2000, il=3,ncell
         itot=itot+il-1
         if(itot.eq.iarnum)NZ=il
         if(itot.eq.iarnum)goto 2010
 2000 continue
 2010 continue


      write (*,"(2x,'# of different spike codes: ',I0)")NZ
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
      
      call getenv ("spkpatdist_skip", cskip)
      read (cskip,'(i10)') iskip
      
      allocate (aryin  (iarnum - iskip, lastc))

      k = 0
      do i=1,iarnum
         if (i.le.iskip) cycle
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
c            if (aryin(i,iz).le.(thrnul*tmpmax)) aryin(i,iz)=0.0
         end do
      end do
      iarnum = iarnum - iskip

      nzero = 0
      if (inc.eq.1) then
         do j=1,ncol,inc
            do i=1,iarnum
               if (aryin(i, j).ne.0) exit
            end do
            if (i.gt.iarnum) nzero = nzero + 1
         end do
      end if
      print *, "nzero: ", nzero

      nzero = 0
      if (inc.eq.1) then
         do j=1,ncol,inc
            if (sum(aryin(:,j)).eq.0) nzero = nzero + 1
         end do
      end if
      print *, "nzero: ", nzero

      nzerotot = 0
      if (inc.eq.1) then
         do j=1,ncol,inc
            do i=1,iarnum
               if (aryin(i, j).eq.0) nzerotot = nzerotot + 1
            end do
         end do
      end if
      print *, "nzerotot: ", nzerotot
      

      end subroutine readx
      end program spkpatdist
      
