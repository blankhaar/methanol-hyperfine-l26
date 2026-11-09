C     ******************************************************
      subroutine expec(vt, sigma, K)
      implicit none
      save
C     ******************************************************

      include 'common_dim.f'
      include 'common_quant.f'
      include 'common_expdat.f'
      include 'common_rotor.f'
      include 'common_dat.f'
      include 'common_tor59.f'
      
      integer vt, k, nn, kind, sigma

      kind = k+Jmax+1

      expc3(vt, kind, sigma+1) = 0.0d0
      expc6(vt, kind, sigma+1) = 0.0d0

      do nn = 1, ndimto-1
         expc3(vt, kind, sigma+1) = expc3(vt, kind, sigma+1) + 
     &        A59(vt, nn, kind, Jmax+1, sigma+1)*
     &        A59(vt, nn+1, kind, Jmax+1, sigma+1)
      end do

      do nn = 1, ndimto-2
         expc6(vt, kind, sigma+1) = expc6(vt, kind, sigma+1) + 
     &        A59(vt, nn, kind, Jmax+1, sigma+1)*
     &        A59(vt, nn+2, kind, Jmax+1, sigma+1)
      end do

      expc3(vt, kind, sigma+1) = 1.0d0 - expc3(vt, kind, sigma+1)
      expc6(vt, kind, sigma+1) = 1.0d0 - expc6(vt, kind, sigma+1)

      return
      end


C     ******************************************************      
      subroutine printexp(m)
      implicit none
      save
C     ******************************************************

      include 'common_dim.f'
      include 'common_quant.f'
      include 'common_expdat.f'
      include 'common_rotor.f'
      include 'common_dat.f'

      integer kind, vt, sigma, k, m

      do sigma = 0,1
         do K = -jmax, jmax

            kind = k+jmax+1

            if ( m .eq. 3 ) then
               write(6, 103) sigma, K, 
     &              expc3(1, kind, sigma+1),
     &              expc3(2, kind, sigma+1),
     &              expc3(3, kind, sigma+1),
     &              expc3(4, kind, sigma+1)
            else if ( m .eq. 6 ) then
               write(6, 103) sigma, K,
     &              expc6(1, kind, sigma+1),
     &              expc6(2, kind, sigma+1),
     &              expc6(3, kind, sigma+1),
     &              expc6(4, kind, sigma+1)
            end if

 103        format(i3, ' ', i3, ' ',
     &           f24.16, ' ', f24.16, ' ', f24.16, ' ', f24.16)
            
         end do
      end do
      
      end
