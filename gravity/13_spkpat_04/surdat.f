      module surdat_mod
      
      private
      save
      double precision, allocatable :: X0 (:,:), X(:,:)
      double precision, allocatable :: scm(:,:), E(:,:), Er(:,:)
      double precision, allocatable :: V(:), W(:), meanobs(:)
      integer, allocatable :: shift(:)
      integer :: nrows = 0, ncols, n
      public surdat, surdat_init, surdat_init_rc, surdat_rc
      
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
      SUBROUTINE SURDAT_INIT_RC (aryin,nrow,ncol)
c     ******************************************************************
      real aryin (:,:)
      double precision, allocatable :: temp(:,:)
      if (allocated (temp)) deallocate (temp)
      allocate (temp(nrow, ncol))
      temp = aryin(1:nrow,1:ncol)
      call surdat_init (temp)
      end subroutine SURDAT_INIT_RC
      
c     ******************************************************************
      SUBROUTINE SURDAT_RC (aryin2,nrow,ncol)
c     ******************************************************************
      real aryin2 (:,:)
      double precision, allocatable :: temp(:,:)
      if (allocated (temp)) deallocate (temp)
      allocate (temp(nrow, ncol))
      call surdat (temp)
      aryin2(1:nrow, 1:ncol) = real(temp)
      end subroutine SURDAT_RC

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

