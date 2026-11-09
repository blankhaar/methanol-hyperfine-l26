      subroutine rule(cst)

C     optional routine to be called at the start of each iteration

      implicit none
      save

      include 'common_dim.f'
      include 'common_dat.f'
      include 'common_chara.f'
      include 'common_rotor.f'
      include 'common_symb.f'
      include 'common_quant.f'
      include 'common_cst.f'

      double precision AK7K
      integer lo

      do lo = 1,2
         AK7K=CST(LO,60)
         
         if ( AK7K .gt. 1.0d-5 ) then
            AK7K = AK7K - 1.0d-5
            CST(LO, 60) = AK7K
            
            print *, 'AK7K = ', AK7K
         end if 
      end do
      
      end
      
