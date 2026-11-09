C*********************************************************
      SUBROUTINE printderivl (CST,NPP, maxvt, maxk)                                               
      implicit none
      save

C     
C*********************************************************

      include 'common_dim.f'
      include 'common_eigen.f'
      include 'common_eigen1.f'
      include 'common_tor.f'
      include 'common_dat.f'
      include 'common_rotor.f'
      include 'common_big2.f'
      include 'common_expdat.f'
      include 'common_symb.f'
      include 'common_quant.f'
      include 'common_chara.f'
      include 'common_cst.f'
      include 'common_int.f'

      integer maxvt, maxk
      integer ipr, ivv, ink,  npp, k,III
      integer ik, vt
      double precision ee
      dimension ee(2)

C     loop through energy levels

      do isig = 0,1
         do vt = 0, maxvt
            do k = 0, maxk
               do n = k, jmax
C                 A
                  if  ( isig .eq. 0 ) then

                     if ( k .ne. 0 ) then
                        ipr = 0
                        ik = k
                        ivv = vt+1
                        ink = n**2 + n + ik + 1 - ipr*n
                        ee(2) =  ener(ivv, ink, 1)
                     else
                        ee(2) = 0.0d0
                     endif
                     
                     ipr = 1
                     ik = k
                     ivv = vt+1
                     ink = n**2 + n + ik + 1 - ipr*n
                     ee(1) = ener(ivv, ink, 1)

                     do iii = 1,npp
                     
                        if ( k .ne. 0 ) then
                           write(6, 561) 'A', vt, n, k, ee(1), ee(2)                  
                        else
                           write(6, 562) 'A', vt, n, k, ee(1)
                        endif
                     end do
C     E
                  else if ( isig.eq. 1 ) then
                     ik = k
                     ivv = vt+6
                     ink = n**2 + n + ik + 1
                     ee(1) = ener(ivv, ink, 1)
                     
                     ik = -k
                     ivv = vt+6
                     ink = n**2 + n + ik + 1
                     ee(2) = ener(ivv, ink, 1)

                     write(6, 561) 'E', vt, n, k, ee(1), ee(2)
                  endif

               end do
            end do
         end do
      end do
      
 561  FORMAT  (A1, ' VT = ', I3, ' N = ', I3, ' K = ',I3,
     &     ' E = ', F16.10, ' +   E = ', F16.10, ' - ') 

 562  FORMAT  (A1, ' VT = ', I3, ' N = ', I3, ' K = ',I3,
     &     ' E = ', F16.10, ' +  ')
      
      return
      end

