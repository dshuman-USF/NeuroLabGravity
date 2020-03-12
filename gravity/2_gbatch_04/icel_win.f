      subroutine icel_win(iiqelen,sper,ifwdinc,mult_ch,iidlist,
     + eevtlist)
C JA Jacques 01/27/94 Version 2.0: removes action potentials b4 computing
C Try to increase size of vlt_arr  - in doing so, got rid of separate
C subroutines and made voltage array allocatable    9/9/98 SC Nuding
C
C 1/99 - Decided to make mean charge calculation over entire block, 
C allow different options to calculate local voltage average, calculate 
C voltage around "jump" value - not just arbitrary array element
C Also - change way charge histories are calculated - see proportion
C of voltage change over entire max-min traverse
C  SC Nuding
C *****************************************************************
C This subroutine loads the N+1 th particle portion of the charge 
C table with values derived from an intracellular recording stored
C in an .HDT file.
C
C This version uses .8 of max voltage to identify an action pot.
C
C       changes in def file in 2004 mods...
c       parameter (ncl=64)!number of codes handled
c       parameter (nsz=480000)!size of charge history array - 1st dimension
c       parameter (nsppc=10000)!max number of spikes in the busiest channel
c       parameter (lstsz=640000)!Product of nsz*nsppc

C       *********************************************************************
        INCLUDE 'ghead.def'                                                 

C Data type declarations

        integer*4 IV_MAX,i_cnt,I,id_idx,time_ptr,e_blk_rec,sav_idx

        parameter (IV_MAX=300000) !maximum size of vlt_arr

        integer*2 vlt_arr(:) !Allocate size later
        Allocatable vlt_arr

        character*15 f_filnam

        integer blk_cnt,arr_max
        integer*2 ap_width,ep_ip,ch_idx,i_chtord,i_style,i_apdel,holder
        
        ! assumes 10 sec max on blk size
        integer*2 s_rate,num_chs,tmp1,tmp2,P,frst_block,sec_block 
        integer*2 offset,loc_win,min_v,max_v,ic_ch,a_minv,a_maxv
        integer*2 diff   !,bit

        real b_blk_tim,e_blk_tim,loc_avg,ifwdinc,ap_recs,loc_max
        real volt_avg,loc_min,qafsall,jump_avg

        dimension iidlist(lstsz),eevtlist(lstsz)


C Format statements
100     format(I2)
115     format(I3)
110     format(I4)

C Read user variables from bcm1 file and open .hdt file
        if (mult_ch .eq. 0) then
           print *,"Enter .hdt file name: "
           read(*,'(A)') f_filnam ! .hdt file
           print *,"Enter # channels to read: "
           read(*,100) i_chtord ! number of channels in .hdt file to read
           print *,"Enter channel to use: "
           read(*,100) ch_idx ! channel in .hdt file to use
           print *,"Enter 1 for EPSPs: "
           read(*,100) ep_ip    ! whether to use EPSP's or IPSP's
           print *,"Enter AP width in msec: "
           read(*,100) ap_width ! width of action potential in msec
           print *,"Enter signal processing option: "
           read(*,100) i_style  ! signal processing option
           print *,"Enter 1 to call ap_del: "
           read(*,100) i_apdel  ! whether to call ap_del
           print *,"Enter # records to offset in block 1 "
           read(*,115) offset   ! # records to skip
           print *,"Enter # for array entries around AP peak "
           read(*,110) loc_win  ! # recs to get better estimate of voltage
           print *,"Enter channel number for intracellular spike train:"
           read(*,100) ic_ch    ! channel carrying spikes for intracellular data
           if (f_filnam .eq. "N") return
        end if
        open (8,file=f_filnam,ACTION='READ',status='old',ACCESS='DIRECT'
     +       ,RECL=2)
c       open (3,file='ap_ok.dat',status='new',ACCESS='SEQUENTIAL')
        if (i_apdel .eq. 1) then 
           open (4,file='ap_del.hdt',status='new',ACCESS='DIRECT',
     +        RECL=2)
           open (5,file='ap_del.dat',status='new',ACCESS='SEQUENTIAL')
        end if 

C Read header info from file (see hdtfmt.dat file for .hdt layout)
        read(8,rec=1) s_rate
        rdg_msec = s_rate/1000 ! gives sampling rate in msec

        IV_HERE = NINT(sper*rdg_msec)
        
              if (IV_HERE .gt. IV_MAX) then !jaj
                print *,"Calculated Array bounds (IV_MAX) exceeded!"
     +             //" Program terminated."
                print *,IV_HERE
                print *,IV_MAX
                goto 999
              end if

        Allocate (vlt_arr (2*IV_HERE))

        read(8,rec=2) num_chs
        ! read number of blocks
        read(8,rec=2*num_chs+3) frst_block
        read(8,rec=2*num_chs+4) sec_block

        num_blks = sec_block + frst_block * 32768
        blk_cnt = 0 

C Now set values for the record number of the beginning of 
C the current block and for the end of the current block. 
        time_ptr = 2*num_chs+5 ! first recno for block times 

        read(8,rec=time_ptr)   tmp1 ! get beginning of first block time
        read(8,rec=time_ptr+1) tmp2
        b_blk_tim = tmp2/2 + (tmp1 * 32768)/2 ! convert to msec

        read(8,rec=time_ptr+2) tmp1 ! get end of block time
        read(8,rec=time_ptr+3) tmp2
        e_blk_tim = tmp2/2 + (tmp1 * 32768)/2
c       print *,'begin, end blks ',b_blk_tim,' ',e_blk_tim

        ! set the number of records to jump for each snapshot to qafs
        rec_jmp = deltim*rdg_msec

        ! number of records per block
        recs_blk = (e_blk_tim - b_blk_tim)*rdg_msec*num_chs

        ! number of 0 values to place in qafs between blocks
        gap_cnt = (sper - (e_blk_tim - b_blk_tim))/deltim

        ! number of records comprising an action potential
        ap_recs = ap_width*rdg_msec

C ************** main program loop **********************************
        id_idx = 0 ! index of iidlist array
c       do i_chcnt = 1,i_chtord ! do for each .hdt channel to read
        
C       add a new particle
        N = N + 1

        ! skip past begin and end time blocks and set record number
        ! to beginning of voltage area for chosen channel
        st_volt = 4 + num_chs*2 + num_blks*4 + ch_idx
        b_blk_rec = st_volt
       print *, "This is using 1000 records for mean, max, min"

        e_blk_rec = 0
        qafsall = 0.0
        arr_max = 0
        L = 1
        do while (L .le. iiqelen .and. blk_cnt .le. num_blks)
           !*************************************
           ! load chg array with 0's for gap
           ! load next block and find min and max voltage
           !*************************************

           if (blk_cnt .gt. 0) then 
           ! don't load at beginning 
              do i = 1,INT(gap_cnt)
                 qafs(L+i,N) = 0
                 qe(L+i,N) = 0
              end do
              L = int(L + gap_cnt)
              b_blk_rec = b_blk_rec + recs_blk
           end if

           e_blk_rec = int(b_blk_rec + recs_blk - num_chs)

           min_v = 32000
           max_v = -32000
           volt_avg = 0
           J = 1

           !*******************************************
           ! to take care of invalid values at
           ! beginning of block 1 (sets offset) 
           !*******************************************
           i_cnt = int(b_blk_rec)
           if (blk_cnt .eq. 0) then
                i_cnt = int(b_blk_rec + offset*num_chs)
           end if

           !*****************************************
           ! read the block into the array vlt_arr()
           !*****************************************
c               print *,'i_cnt = ', i_cnt
c               print *,'e_blk_rec', e_blk_rec
c               print *,'blk_cnt = ',blk_cnt
           do I = i_cnt, e_blk_rec, num_chs 
              read(8,rec=I,err=300) vlt_arr(J)
              if (min_v .gt. vlt_arr(J)) min_v = vlt_arr(J)
              if (max_v .lt. vlt_arr(J)) max_v = vlt_arr(J)
              volt_avg = volt_avg + vlt_arr(J)
              J = J + 1
                      if (J .gt. IV_HERE) then !jaj
                print *,"Array bounds (IV_HERE) exceeded!"
     +             //" Program terminated."
                print *,J
                print *,IV_HERE
                print *,i_cnt
                print *,e_blk_rec
                print *,num_chs
                print *,I
                goto 990
                      end if
            end do
300        continue
c          print *,e_blk_rec
           arr_max = J - 1
           volt_avg = volt_avg/arr_max
c          print *,'arr_max = ', arr_max
           blk_cnt = blk_cnt + 1

c             do bit = 1,arr_max
c               write(3,110) vlt_arr(bit)
c             end do
c          print 301,min_v, max_v, volt_avg 
c           print *,'blk_cnt = ',blk_cnt
c301        Format ('min_v, max_v, volt_avg ',I5,' , ',I5,' ,',F10.3)
c*******************************************************************
           ! start process to eliminate action potentials
C*******************************************************************

        If (i_apdel.eq.1) then 
           write(4,rec=1) s_rate
           write(4,rec=2) num_chs

        ! read min and max voltages for each channel
           do i = 3,2*num_chs+2
              read(8,rec=i) tmp1
              write(4,rec=i) tmp1
           end do
        ! write number of blocks
           write(4,rec=2*num_chs+3) frst_block
           write(4,rec=2*num_chs+4) sec_block
        
        !**************************************************
        ! read each begin and end block time for each block
        ! each value is stored in 2 2-byte integers
        ! this is for creating the output .hdt file only
        ! doesn't affect the gravity algorithm
        !**************************************************
           i_cnt  = 2*num_chs+5
           i_cnt2 = 2*num_chs + 4*num_blks + 4
           do i = i_cnt,i_cnt2
              read(8,rec=i) tmp1
              write(4,rec=i) tmp1
           end do
           
        !***************************************************
        ! write all the voltage values for all channels 
        ! ap_del will overwrite those records which are 
        ! changed
        ! Modified 1/99 to handle hdt files of max size SCN
        !***************************************************
           i_cnt = i_cnt2 + 1
           holder = 0
           do while (holder.eq.0)
              read(8,rec=i_cnt,err=200) tmp1
              write(4,rec=i_cnt) tmp1
              i_cnt = i_cnt + 1
           end do
           
200     continue
        
C ************** old subroutine ap_del() ******************************
C Deletes action potentials from the array and replaces them with a
C straight line connecting the ends (using local window average voltage)
C *********************************************************************
           
           id_idx = id_idx + 1
           st_time = eevtlist(id_idx) !value at ID 21
        do while (iidlist(id_idx) .ne. 22) ! end of block
           if (iidlist(id_idx) .eq. ic_ch) then ! site of AP in .gdt file
              !*****************************************************************
              ! find the peak in vlt_arr corresponding to the spike in eevtlist.
              ! Must use rdg_msec to calculate offset.
              !*****************************************************************
              I = int(rdg_msec * (eevtlist(id_idx) - st_time))

              !*****************************************************************
              ! The above gets us close to the action potential
              ! Now find peak within ap_recs window
              !*****************************************************************
              i_st = int(I-ap_recs)
              if (i_st .lt. 1) i_st = 1

              i_end = int(I+ap_recs)
              if (i_end .gt. arr_max) i_end = arr_max

              sav_idx = I
              do J= i_st,i_end
                 if (vlt_arr(J) .gt. vlt_arr(sav_idx))sav_idx = J
              end do

              ! ********************************************************
              ! The peak is at index sav_idx. Delete the left and right
              ! halves of the action potential. Then delete the after
              ! hyperpolarization section until the mean is reached.
              ! Connect the two loose ends with a straight line of value
              ! "loc_avg" . Do this only if sav_idx is a peak. If 2 APs
              ! are very close, the first go-round eliminates the peak, so
              ! don't do it twice.
              ! ********************************************************
              if (vlt_arr(sav_idx) .gt. vlt_arr(i_st) .and.
     +            vlt_arr(sav_idx) .gt. vlt_arr(i_end)) then

                 i_st = int(sav_idx - ap_recs/2)
                 i_end = int(i_st + 2*ap_recs)
                 if (i_st .lt. 1) i_st = 1 ! don't get out of array bounds
                 if (i_end .gt. arr_max) i_end = arr_max

                 ! set new values to local mean
C*********************** old subroutine local_mn but still in ap-del************
C Finds local mean in vlt_arr()
C *******************************************************************
                 i_st2 = sav_idx - loc_win/2
                 i_end2 = sav_idx + loc_win/2
                 if (i_st2 .lt. 1) then
                   i_st2 = 1
                   i_end2 = loc_win
                 else 
                   if (i_end2 .gt. arr_max) then
                     i_st2 = arr_max - loc_win
                     i_end2 = arr_max
                   end if
                 end if

                 loc_max = -32000
                 loc_min = 32000
                 loc_avg = 0.0
                 do K = i_st2,i_end2
                    loc_avg = loc_avg + vlt_arr(K)
                    if (loc_max .lt. vlt_arr(K)) loc_max = vlt_arr(K)
                    if (loc_min .gt. vlt_arr(K)) loc_min = vlt_arr(K)
                 end do
                 loc_avg = loc_avg/loc_win
c       print *,"loc_avg,loc_win,loc_max,loc_min ",loc_avg," ",loc_win
c     +         ," ",loc_max," ",loc_min
C*********************** old subroutine local_mn ********************
C end of old subroutine local_mn .... within ap_del
C *******************************************************************
        
                 do NP = i_st,i_end
                    vlt_arr(NP) = int(loc_avg,kind(vlt_arr(NP)))
                 end do
                 NP = i_end + 1
                 do while (vlt_arr(NP) .lt. loc_avg 
     +              .and. NP .le. arr_max)
                    vlt_arr(NP) = int(loc_avg,kind(vlt_arr(NP)))
                    NP = NP + 1
                 end do
               end if ! If within limits
               end if ! If ic_ch has spike
            id_idx = id_idx + 1
        end do ! While not end of block
        end if ! If action potential delete
C*********************** old subroutine ap_del **********************
C end of old subroutine ap_del ....
C *******************************************************************

        !********************************!
        ! find new max, min, mean values !
        !********************************!
           a_minv = 32000
           a_maxv = -32000
           a_tot = 0.0

           do J = 1,arr_max
c            print *,vlt_arr(J)
             a_tot = a_tot + vlt_arr(J)
             if (a_minv .gt. vlt_arr(J)) a_minv = vlt_arr(J)
             if (a_maxv .lt. vlt_arr(J)) a_maxv = vlt_arr(J)
           end do
           a_mean = a_tot/arr_max ! a_mean is the mean voltage over block
c          print 302, a_minv, a_maxv, a_mean 
c302        Format ('a_minv, a_maxv, a_mean ',I4,' , ',I4,' ,',F10.3)
c           print *,'a_tot,arr_max = ',a_tot,' ',arr_max

           if (i_style .ne. 1) then
c          print *,'into rectify'   
C ************** old subroutine chc_rect() **************************
C Half wave or full wave rectifies
C       i_style = 2 : top half wave rectify
C       i_style = 3 : bottom half wave rectify
C       i_style = 4 : full wave rectify
C *******************************************************************

        !*******************!
        ! half wave rectify !
        !*******************!
           if (i_style .eq. 2) then ! top rectify
              do J = 1,arr_max
                if (vlt_arr(J) .le. a_mean) vlt_arr(J) = 0
              end do
           else if (i_style .eq. 3) then ! bottom rectify
              do J = 1,arr_max
                 if (vlt_arr(J).ge.a_mean) then
                          vlt_arr(J)=0 
                 else 
                      vlt_arr(J)
     +                   =int(abs(a_mean-vlt_arr(J)),kind(vlt_arr(J)))
                 end if
              end do
           else ! full wave rectify
             do J = 1,arr_max
             if (vlt_arr(J).lt.a_mean) vlt_arr(J)
     +               =int(abs(a_mean-vlt_arr(J)),kind(vlt_arr(J)))
             end do
            end if  ! different choices
         end if ! i-style choices
c         print *,'past rectify'
C ************** end of old subroutine chc_rect() *******************

           ! write modified data to ap_del.hdt
           if (i_apdel .eq. 1) then
              j = 1
              do I = INT(b_blk_rec),e_blk_rec,num_chs
                write(4,rec=I) vlt_arr(j)
                write(5,110) vlt_arr(j)
                j = j + 1
              end do
           end if

           ! **************************************************
           ! calc local avg and scale to gmint32.f chg values
           ! **************************************************
        do P=INT((rec_jmp)/2,kind(P)),INT(arr_max-(rec_jmp/2),kind(P))
     +          ,INT(rec_jmp,kind(P))

        i_st = int((P - rec_jmp/2)+1)
        i_end = int(P + rec_jmp/2)
          if (i_end .gt. arr_max) then
           i_st = int(arr_max - rec_jmp+1)
           i_end = arr_max
          end if

        jump_avg = 0.0
        do K = i_st,i_end
           jump_avg = jump_avg + vlt_arr(K)
c          print *,vlt_arr(K)
        end do
         jump_avg = jump_avg/rec_jmp
c          print 305,P,vlt_arr(P),jump_avg
c305        format (i5,1x,i5,1x,(f6.1))


        i_st = P - loc_win/2
        i_end = P + loc_win/2
        if (i_st .lt. 1) then
             i_st = 1
             i_end = loc_win
        else 
            if (i_end .gt. arr_max) then
             i_st = arr_max - loc_win +1
             i_end = arr_max
            end if
        end if
c        print *,'i_st,i_end ',i_st,' ',i_end

        loc_max = -32000
        loc_min = 32000
        loc_avg = 0
        do K = i_st,i_end
           loc_avg = loc_avg + vlt_arr(K)
           if (loc_max .lt. vlt_arr(K)) loc_max = vlt_arr(K)
           if (loc_min .gt. vlt_arr(K)) loc_min = vlt_arr(K)
        end do
         loc_avg = loc_avg/loc_win
c       print *,"loc_avg,max,min ",loc_avg," ",loc_max," ",loc_min

C      Decide if want local average for gravity mean charge zero, or use a
C      fraction of it, or use entire block??
C      Run with next three lines in or commented out and compare
       a_mean = loc_avg
       a_minv = int(loc_min,kind(a_minv))
       a_maxv = int(loc_max,kind(a_maxv))

              !**************************************************
              ! determine whether IPSP's or EPSP's and scale vals
              !Change way charges are calculated - look at diff.
              !between voltage and max-min travers ... SCN
              !**************************************************
              diff = a_maxv - a_minv
              if (diff.ne.0) then
                qafs(L,N) = ifwdinc*
     +          (jump_avg-a_mean)/(a_maxv-a_minv)
              else
                qafs(L,N) = qafs(L-1,N)
              end if
c             if (jump_avg .le. a_mean) qafs(L,N) = -qafs(L,N)

              !***************************************
              ! if IPSP's desired, reverse the charges
              !***************************************
              if (ep_ip .eq. 0) then 
                 qafs(L,N) = -qafs(L,N)
              end if
              qafsall = qafsall + qafs(L,N)
c             print *,'diff, jump_avg ', diff,' ',jump_avg
c             print *,'qafs,qafsall = ',L,' ',qafs(L,N),' ',qafsall
              qe(L,N) = qafs(L,N)
              L = L + 1
              if (L .gt. iiqelen) return
           end do ! for P
c       print *,'L = ',L
        end do ! for L
c       blk_cnt = 1
c       print *,"to 990"
c       print *,blk_cnt-1
c       end do ! for each .hdt channel to read
C ************** end main program loop ******************************

990     if (i_apdel .eq. 1) then
            close(4)
            close(5)
        end if
c       close(3)
        Deallocate (vlt_arr)
999     close(8)
        return
        end



