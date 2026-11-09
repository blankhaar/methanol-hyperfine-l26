C***********************************************************                    
      FUNCTION ENCAL2(E,LO,EGVL,EGVC,A, SK)                                   
      implicit none
      save
C
C     J M Fisher
C     New labelling scheme (not mine!), see Illyushin
C
C***********************************************************                    
      
      include 'common_dim.f'
      include 'common_eigen.f'
      include 'common_rotor.f'
      include 'common_quant.f'
      include 'common_chara.f'
      include 'common_expdat.f'
      include 'common_cntl.f'
      
      double precision egvl, egvc, egvc4, a, e
      integer ivtorp, iu, ii, lo, j, i, ivt, ik, sk
      integer idx1, idx2
      integer encal2
      double precision sf, sfmax, deltaE, predE, factE
      double precision factor, sfEmax

      double precision N1, N2, N3, NN, E1, E2, E3, EN
      
      DIMENSION EGVL(NDMX),EGVC(NDMX,NDMX),A(9,NTORMX,KDMX)                     
      DIMENSION EGVC4(KDMX,NTORMX,NDMX)

      sfmax = 0.0d0
      sfEmax = 0.0d0
      j = 0

C     EGVLJ has the last three energies (in J)
C     use a 2nd order Lagrange interpolation to predict this energy
C     calculate a scaling factor exp(-|(Epred-E)/E|) to multiply the SF by

      N1 = N-3.d0
      N2 = N-2.d0
      N3 = N-1.d0
      NN = N

      N1 = N1*(N1+1.d0)
      N2 = N2*(N2+1.d0)
      N3 = N3*(N3+1.d0)
      NN = NN*(NN+1.d0)

      E1 = EGVLJ(IV+1, ISIG+1, SK+JMAX+1, 1)
      E2 = EGVLJ(IV+1, ISIG+1, SK+JMAX+1, 2)
      E3 = EGVLJ(IV+1, ISIG+1, SK+JMAX+1, 3)

      EN = E1*(NN-N2)*(NN-N3) / ( (N1-N2)*(N1-N3) ) 
      EN = EN + E2*(NN-N1)*(NN-N3) / ( (N2-N1)*(N2-N3) )
      EN = EN + E3*(NN-N1)*(NN-N2) / ( (N3-N1)*(N3-N2) )

C     loop through all eigenvectors
      do i = 1, nrotor
         E = EGVL(i)

         factE = exp( -abs( (EN-E)/E ) )

         sf = 0.0d0
C     get the similarity factor
         do ivt = 1, 9
            do ii = -N+1, N-1

               idx1 = (ivt-1)*(2*(N-1)+1) + ii+(N-1)+1
               idx2 = (ivt-1)*(2*N+1) + ii+N+1
               
               sf = sf + EGVC3J(IV+1, ISIG+1,
     &              SK+JMAX+1, idx1) * EGVC(idx2, i)
            end do
         end do

         sf = abs(sf)

C     larger SF?
         if ( label .eq. 2 .and. (abs(sk) <= (N-3)) ) then
            if ( sf*factE .gt. sfEmax ) then
               sfmax = sf
               sfEmax = sf*factE
               j = i
               factor = factE
            end if
         else
            if ( sf .gt. sfmax ) then
               sfmax = sf
               j = i
               factor = factE
               sfEmax = sf*factE
            end if
         end if
      end do


 162  format('VT = ', I2, ' ISIG = ', I2, ' N = ', I2, ' K = ', I3, 
     &' SF = ', F24.16, ' EN = ', F24.16)

C      write(6, 162) iv, isig, n, kk, sfmax, factor

      encal2 = j

      RETURN                                                                    
      END                                                                       


