      subroutine postodir3d (ICFLAG,sig)
      
c     FOR GRAVITATIONAL CLUSTERING PROGRAMS
c     CONVERTS POSIT VS TIME TO DIST VS TIME
c     INPUT FROM GRAV6C OR D (GPOS6D OR GNEG6D) AND UP (E,F,...)
c     
c     V1.01 12-DEC-83 GLG USE FOR PLOTTER PDIS6C INPUT
c     V2.01 14-FEB-84 GLG VERSION FOR DG
c     V2.02 1-MAY-85  GLG EATS AND IGNORES FLASH VECTOR;
c     FOR GEDTEST2 AND LATER
c     V2.03 6-DEC-85 GLG MAKES N A TYPEIN. NEEDED FOR GFASTX AND GMINTX
c     V2.1  11-AUG-86 BGL TYPE IN VERSION OF PGM FOR INIT DEV ON 11/23
c     V3.0  25-AUG-86 BGL CHANGE READ AND WRITE TO DIRECT ACCESS,BUFERED I/O
c     V4.0  16-JUL-87 FIRST VERSION FOR HP350..REINSTATE G'S CODE
c     v5.0  30 oct 87 convert to subroutine
c     3/16/94 mod to possub32 to handle 32 codes bgl
c     
c     mods for first big version bgl 12-26-96
c     
c     changes in def file in 2004 mods...
c     parameter (ncl=64)!number of codes handled
c     parameter (nsz=480000)!size of charge history array - 1st dimension
c     parameter (nsppc=10000)!max number of spikes in the busiest channel
c     parameter (lstsz=640000)!Product of nsz*nsppc
c     
c     VERSION TO GO WITH G*2000.F FILES BGL 12-20-99
c     *********************************************************************
      INCLUDE 'ghead.def'                                                        
c     **********************************************************************

      DIMENSION POSIT(ncl,ncl),PFDIR(ncl,ncl,ncl)
      DIMENSION PDIS(ncl,ncl),FLASH(ncl)
      dimension dis(ncl),sumdsq(ncl)
      character*30 dirfile, garbage
      real origin
      logical sig
      integer stat
      
      OPEN(UNIT=1,FILE=OFNAME,STATUS='OLD',FORM='UNFORMATTED')
      read (1) CBLOCK,IBLOCK,RBLOCK !read in header info
     
c     ****SET UP OUTPUT FILE
c     PRINT 35
c     35     FORMAT(2X,'OUTPUT *.dir FILENAME IS: ',$)
c     READ (*,'(A)') dirfile
     
      if (sig) then
         if (ICFLAG.EQ.1) READ (*,'(A)') garbage
         write (dirfile, "('sh',I0,'.dir')") ICFLAG
      else
         PRINT "(2X,'OUTPUT PAIR DIST VS TIME FILENAME IS: ',$)"
         READ (*,'(A)') dirfile
      end if
      OPEN(UNIT=2,FILE=dirfile,STATUS='NEW',FORM='UNFORMATTED')
      write (2) CBLOCK,IBLOCK,RBLOCK !write header to "dir file" file
     
c     calculate the value of the non-zero origin coordinate
c     in N dimension space
      origin = (SQRT (2.0)/2.0)*100.0
     
      do while (.true.)
         READ(1,iostat=stat) ((POSIT(I,J),J=1,N),I=1,N)
         if (stat.ne.0) exit
         READ(1,iostat=stat) (FLASH(I),I=1,N)
         if (stat.ne.0) exit
         
c     GET DISTANCE OF EACH PARTICLE FROM ITS ORIGIN
c     AT THIS TIME STEP
         do  I=1,N
            sumdsq(I)=0.0
            do  J=1,N
               if (I.eq.J) then
                  sumdsq(I)=sumdsq(I)+(POSIT(I,J)-origin)**2
               else
                  sumdsq(I)=sumdsq(I)+POSIT(I,J)**2
               end if
            end do
            dis(I)=SQRT(sumdsq(I)) ! = distance of particle I from its origin
         end do
     
c     DIRECTION BETWEEN POINT OF PAIR
         DO I=2,N
            DO J=1,I-1
               DO K=1,N
                  PFDIR(I,J,K)=POSIT(J,K)-POSIT(I,K)
               end do
            end do
         end do
     
c     CALC DISTANCE BETWEEN POINTS OF PAIR
         DO I=2,N
            DO J=1,I-1
               PDIS(I,J)=0.0
               DO K=1,N
                  PDIS(I,J)=PDIS(I,J)+PFDIR(I,J,K)**2 ! SQUARED DISTANCE
               end do
               PDIS(I,J)=SQRT(PDIS(I,J)) ! ROOT IT
            end do
         end do
     
c     WRITE DIST MATRIX TO OUTPUT FILE
         WRITE(2) ((PDIS(I,J),J=1,I-1),I=2,N)
         WRITE(2) (dis(I),I=1,N)
      end do
     
      PRINT "(2X,'NORMAL EXIT..postodir3d..')"
      CLOSE (UNIT=1)
      CLOSE (UNIT=2)
     
      return
      END

