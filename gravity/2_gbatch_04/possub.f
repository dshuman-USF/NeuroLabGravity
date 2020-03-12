      subroutine possub (ICFLAG,sig)
      
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
c     v10a  bgl outputs 100 files  
c     v32  bgl outputs 100 files for upto 32 particles
c     changes in def file in 2004 mods...
c     parameter (ncl=64)!number of codes handled
c     parameter (nsz=480000)!size of charge history array - 1st dimension
c     parameter (nsppc=10000)!max number of spikes in the busiest channel
c     parameter (lstsz=640000)!Product of nsz*nsppc
c     
c     changes in *_04 version to handle 64 units,particles
c     *********************************************************************
      INCLUDE 'ghead.def'
c     **********************************************************************
      DIMENSION POSIT(ncl,ncl),PFDIR(ncl,ncl,ncl),PDIS(ncl,2*ncl), ! 2*ncl from possubsig v. 0.10.3
     +     FLASH(ncl)
      character*30 posfile,garbage
      integer ICFLAG,stat
      logical sig
c     
c     
      OPEN(UNIT=1,FILE=OFNAME,STATUS='OLD',FORM='UNFORMATTED')
      read (1) CBLOCK,IBLOCK,RBLOCK !read in header info
c     
c     ****SET UP OUTPUT FILE
c     PRINT 35
c     35    FORMAT(2X,'OUTPUT PAIR DIST VS TIME FILENAME IS: ',$)
      if (sig) then
         if (ICFLAG.EQ.1) READ (*,'(A)') garbage
         write (posfile, "('sh',I0,'.pos')") ICFLAG
      else
         PRINT "(2X,'OUTPUT PAIR DIST VS TIME FILENAME IS: ',$)"
         READ (*,'(A)') posfile
      end if

c     
      OPEN(UNIT=2,FILE=posfile,STATUS='NEW',FORM='UNFORMATTED')
      write (2) CBLOCK,IBLOCK,RBLOCK !write header to "distance" file
      write (2) TIMAT           ! spike times for "raw data plot" over time vs. 
c     distance plot

      do while (.true.)
         READ(1,iostat=stat) ((POSIT(I,J),J=1,N),I=1,N)
         if (stat.ne.0) exit
         READ(1,iostat=stat) (FLASH(I),I=1,N)
         if (stat.ne.0) exit
          
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
         
c     **********WRITE DIST MATRIX TO OUTPUT FILE
         WRITE(2) ((PDIS(I,J),J=1,I-1),I=2,N)
      end do
     
c     ****NORMAL EXIT
      PRINT "(2X,'NORMAL EXIT..possub..')"
      CLOSE (UNIT=1)
      CLOSE (UNIT=2)
c     
      return
      END

