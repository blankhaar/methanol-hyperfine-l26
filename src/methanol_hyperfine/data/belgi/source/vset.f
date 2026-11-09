C     *************************************************************************
      subroutine vset1(lo,v,cst)                                               
C     *************************************************************************
      
C     *************************************************************************  
      implicit none
      save
 
      include 'common_cst.f'
      include 'common_dim.f'
      include 'common_eigen.f'
      include 'common_eigen1.f'
      include 'common_tor.f'
      include 'common_dat.f'
      include 'common_rotor.f'
      include 'common_par.f'
      include 'common_big2.f'
      include 'common_chara.f'
      include 'common_expdat.f'
      include 'common_symb.f'
      include 'common_quant.f'
      include 'common_vset.f'

      double precision v

      integer lo, iii, i, k
      DIMENSION V(maxparam)          

      include 'common_vrtd.f'

                                                   
      DO I = 1, 2 * N + 1                                                    
         K = I - N - 1                                                          
         DO III = 1, NDIMTO                                                  
            RHOEFF = RHO                                                        
            R2A(III,I) = 3.d0*(III-KTRON1) + ISIG + RHOEFF*K                    
            R2A9(III,I) = R2A(III,I)**3

            R2AC(III,I) = 3.d0*(III-KTRON1) + ISIG + RHOEFF*(K+1)                    

            R2X        = 3.d0*(III-KTRON1) + ISIG + RHOEFF*(K+2)                
            R2B(III,I) = R2A(III,I)*R2A(III,I) + R2X*R2X                        
            R2BK(III,I) = R2A(III,I)**2*K**2 + R2X**2*(K+2)**2                       
            R2BKK(III,I) = R2A(III,I)**2*K**4 + R2X**2*(K+2)**4                        
            R2C(III,I) = R2A(III,I)**4 + R2X**4                                 
            R2CK(III,I) = R2A(III,I)**4*K**2 + R2X**4*(K+2)**2
            R2CKK(III,I) = R2A(III,I)**4*K**4 + R2X**4*(K+2)**4
            R2D(III,I) = R2X*(K+2) + R2A(III,I)*K                             
*                                                                               
            R4X = 3.d0*(III-KTRON1) + ISIG + RHOEFF*(K+1)                       
            R4A(III,I) = R4X*(K+1)*(K+1)+R2A(III,I)*K*K                         
            R8A(III,I) = R4X*(K+1)**4+R2A(III,I)*K**4                         
            R4AB(III,I)=R4X**3*(K+1)*(K+1)+R2A(III,I)**3*K*K
            R4B(III,I) = R2A(III,I)*R2A(III,I)*K + R4X*R4X*(K+1)                
            R7B(III,I) = R2A(III,I)*R2A(III,I)*K**3 + R4X*R4X*(K+1)**3                
            R8B(III,I) = R2A(III,I)*R2A(III,I)*K**5 + R4X*R4X*(K+1)**5                
            R4BB(III,I)=(R2A(III,I))**4*K+(R4X)**4*(K+1)
            R4BC(III,I)=(R2A(III,I))**4*K**3+(R4X)**4*(K+1)**3
            R2E(III,I)=R2A(III,I)**3*(K)+R2X**3*(K+2)
            R2EK(III,I)=R2A(III,I)**3*(K**3)+R2X**3*(K+2)**3
            R2EKK(III,I)=R2A(III,I)**3*(K**5)+R2X**3*(K+2)**5
C            R2F(III,I)=R2A(III,I)*K*(K+2)**2+R2X*(K+2)
C            *(K)**2

            R2F(III,I)=R2A(III,I)*K**3+R2X*(K+2)**3
            R2FK(III,I)=R2A(III,I)*K**5+R2X*(K+2)**5
                                                                               
            R5A(III,I) = R2A(III,I)*R2A(III,I)                                  
            R5A8(III,I) = R2A(III,I)**4

            R5B(III,I) = R2A(III,I)*R2A(III,I)*R2A(III,I)                       
            R5B2(III,I) = R2A(III,I)*R2A(III,I)*R5B(III,I)                      
            R5C(III,I) = R2A(III,I)*R2A(III,I)*R2A(III,I)*R2A(III,I)            
            R5C2(III,I) = R2A(III,I)*R2A(III,I)*R5C(III,I)                      
            R5D(III,I) = 3.d0*(III-KTRON1+1) + ISIG + RHOEFF*K                  
            R5D9(III,I) = R5D(III, I)**3

            R5DD(III,I) = 3.d0*(III-KTRON1+2) + ISIG + RHOEFF*K                 

            R5E(III,I) = R5D(III,I)*R5D(III,I)                                  
            R5E8(III,I) = R5D(III,I)**4                                  

            R5F(III,I) = 3.d0*(III-KTRON1-1) + ISIG + RHOEFF*K                  
            R5F9(III,I) = R5F(III, I)**3

            R5FF(III,I) = 3.d0*(III-KTRON1-2) + ISIG + RHOEFF*K                 

            R5G(III,I) = R5F(III,I)*R5F(III,I)                                 
            R5G8(III,I) = R5F(III,I)**4                   
*                                                                               
            RVA(III,I) = 3.d0*(III-KTRON1) + ISIG + RHORHO*K                    
            RVB(III,I) = RVA(III,I)*RVA(III,I)                                  

            R2B1(III, I) = 3.d0*(III-KTRON1) + ISIG + RHOEFF*K
            R2B6(III, I) = 3.d0*(III-KTRON1) + ISIG + RHOEFF*(K+1)
            R2B2(III, I) = 3.d0*(III-KTRON1) + ISIG + RHOEFF*(K+2)

            R2B5(III, I) = 3.d0*(III-KTRON1+1) + ISIG + RHOEFF*K
            R2B4(III, I) = 3.d0*(III-KTRON1+1) + ISIG + RHOEFF*(K+1)
            R2B3(III, I) = 3.d0*(III-KTRON1+1) + ISIG + RHOEFF*(K+2)

            R2B1(III, I) = R2B1(III, I)*R2B1(III, I)
            R2B2(III, I) = R2B2(III, I)*R2B2(III, I)
            R2B3(III, I) = R2B3(III, I)*R2B3(III, I)
            R2B4(III, I) = R2B4(III, I)*R2B4(III, I)
            R2B5(III, I) = R2B5(III, I)*R2B5(III, I)
            R2B6(III, I) = R2B6(III, I)*R2B6(III, I)
           
            R2G1(III,I) = R2A(III,I)**6 + R2X**6                        
            R2G2(III,I) = K*R2A(III,I)**5 + R2X**5 * (K+2)                        
         end do
      end do
*                      

      return 
      end
C     *************************************************************************  


C     *************************************************************************  
      subroutine vset2
C     *************************************************************************  

C     A little bit of profiling revealed that this function takes about 80%
C     of the total runtime. It looks like it could probably parallelize, which
C     would give a boost in performance at the cost of portability.
C
C     JM Fisher - By splitting the largest loop into two loops I was able to
C     get it down to about 60%. The smaller loop bodies seem to have a large
C     impact on performance. 
C     *************************************************************************  
      implicit none
      save

      include 'common_cstdim.f'
      include 'common_dim.f'
      include 'common_eigen.f'
      include 'common_eigen1.f'
      include 'common_tor.f'
      include 'common_dat.f'
      include 'common_rotor.f'
      include 'common_par.f'
      include 'common_big2.f'
      include 'common_chara.f'
      include 'common_expdat.f'
      include 'common_symb.f'
      include 'common_quant.f'
      include 'common_vset.f'

c      double precision v

c      integer lo
      integer ivtorp, ivtor1, i
      integer i130, i140, i150, i160, i170, i180, i190
c      DIMENSION V(maxparam)           

      include 'common_vrtd.f'

                                                         
      do IVTORP = 1,9                                                      
         DO  IVTOR1 = 1,9                                                   
            DO  I = 1, 2 * N + 1                        

*     Pour VRTDI1                                                                   
               do i130 = 1, ndimto
                  AAA = A(IVTORP,I130,M-1+I)*A(IVTOR1,I130,M-1+I)                 
                  A1T1(I,IVTOR1,IVTORP) = A1T1(I,IVTOR1,IVTORP) + AAA           
                  AK7T1(I,IVTOR1,IVTORP) = AK7T1(I,IVTOR1,IVTORP) 
     &                 + AAA * (R2A(I130, I)**2 + R2AC(I130, I)**2)
               end do

*     Pour VRTDI2                                                                   
               do i130 = 1, ndimto
                  A1T31(I,IVTOR1,IVTORP) = A1T31(I,IVTOR1,IVTORP) +
     &                 A(IVTORP,I130,M-1+I)*A(IVTOR1,I130,M+I+3)

                  A1T32(I,IVTOR1,IVTORP) = A1T32(I,IVTOR1,IVTORP) +
     &                 A(IVTORP,I130,M+I+3)*A(IVTOR1,I130,M-1+I)

                  AA2 = A(IVTORP,I130,M-1+I)*A(IVTOR1,I130,M+I+1)                 

                  A1T2(I,IVTOR1,IVTORP) = A1T2(I,IVTOR1,IVTORP) + AA2           
                  A2T2(I,IVTOR1,IVTORP) = A2T2(I,IVTOR1,IVTORP)                 
     &                 + AA2 * R2B(I130,I)                                       
                  A2T2K(I,IVTOR1,IVTORP) = A2T2K(I,IVTOR1,IVTORP)                 
     &                 + AA2 * R2BK(I130,I)                                       
                  A2T2KK(I,IVTOR1,IVTORP) = A2T2KK(I,IVTOR1,IVTORP)                 
     &                 + AA2 * R2BKK(I130,I)                                       
                  A2TG1(I,IVTOR1,IVTORP) = A2TG1(I,IVTOR1,IVTORP)
     &                 + AA2 * R2G1(I130,I)
                  A2TG2(I,IVTOR1,IVTORP) = A2TG2(I,IVTOR1,IVTORP)
     &                 + AA2 * R2G2(I130,I)
                  A3T2(I,IVTOR1,IVTORP) = A3T2(I,IVTOR1,IVTORP)                 
     &                 + AA2 * R2C(I130,I)                                       
                  A3T2K(I,IVTOR1,IVTORP) = A3T2K(I,IVTOR1,IVTORP)                 
     &                 + AA2 * R2CK(I130,I)                                       
                  A3T2KK(I,IVTOR1,IVTORP) = A3T2KK(I,IVTOR1,IVTORP)                 
     &                 + AA2 * R2CKK(I130,I)                                       
                  A7T2(I,IVTOR1,IVTORP) = A7T2(I,IVTOR1,IVTORP)                 
     &                 + AA2 * R2D(I130,I)                                       
                  A10T2(I,IVTOR1,IVTORP) =A10T2(I,IVTOR1,IVTORP)                
     &                 + AA2 * R2E(I130,I)                                       
                  A10T2K(I,IVTOR1,IVTORP) =A10T2K(I,IVTOR1,IVTORP)                
     &                 + AA2 * R2EK(I130,I)                                       
                  A10T2KK(I,IVTOR1,IVTORP) =A10T2KK(I,IVTOR1,IVTORP)                
     &                 + AA2 * R2EKK(I130,I)                                       
                  A11T2(I,IVTOR1,IVTORP) =A11T2(I,IVTOR1,IVTORP)                
     &                 + AA2 * R2F(I130,I)                                       
                  A11T2K(I,IVTOR1,IVTORP) =A11T2K(I,IVTOR1,IVTORP)                
     &                 + AA2 * R2FK(I130,I)                                       
               end do

*     Pour VRTDI4                                                                   
               do i130 = 1, ndimto
                  AA4 = A(IVTORP,I130,M-1+I)*A(IVTOR1,I130,M+I)                   
                  A1T4(I,IVTOR1,IVTORP) = A1T4(I,IVTOR1,IVTORP) + AA4           
                  A1T4B(I,IVTOR1,IVTORP) = A1T4B(I,IVTOR1,IVTORP)
     &                 + AA4 * ( R2A(I130, I) + R2AC(I130, I) )
                  A5T4(I,IVTOR1,IVTORP) = A5T4(I,IVTOR1,IVTORP)                 
     &                 + AA4 * R4A(I130,I)                                       
                  A8T4(I,IVTOR1,IVTORP) = A8T4(I,IVTOR1,IVTORP)                 
     &                 + AA4 * R8A(I130,I)                                       
                  A5T4B(I,IVTOR1,IVTORP) = A5T4B(I,IVTOR1,IVTORP)               
     &                 + AA4 * R4AB(I130,I)                                      
                  A6T4(I,IVTOR1,IVTORP) = A6T4(I,IVTOR1,IVTORP)                 
     &                 + AA4 * R4B(I130,I)                                       
                  A7T4(I,IVTOR1,IVTORP) = A7T4(I,IVTOR1,IVTORP)                 
     &                 + AA4 * R7B(I130,I)                                       
                  A9T4(I,IVTOR1,IVTORP) = A9T4(I,IVTOR1,IVTORP)                 
     &                 + AA4 * R8B(I130,I)                                       
                  A6T4B(I,IVTOR1,IVTORP) = A6T4B(I,IVTOR1,IVTORP)               
     &                 + AA4 * R4BB(I130,I)                                      
                  A6T4C(I,IVTOR1,IVTORP) = A6T4C(I,IVTOR1,IVTORP)               
     &                 + AA4 * R4BC(I130,I)                                      
               end do

*     Pour VRTDI5                                                                                  
               do i130 = 1, ndimto
                  AAA = A(IVTORP,I130,M-1+I)*A(IVTOR1,I130,M-1+I)                 
                  A2T5(I,IVTOR1,IVTORP) = A2T5(I,IVTOR1,IVTORP)                 
     &                 + AAA * R5A(I130,I)                                       
                  A2T58(I,IVTOR1,IVTORP) = A2T58(I,IVTOR1,IVTORP)                 
     &                 + AAA * R5A8(I130,I)                                       
                  A3T5(I,IVTOR1,IVTORP) = A3T5(I,IVTOR1,IVTORP)                 
     &                 + AAA * R2A(I130,I)                                       
                  A3T59(I,IVTOR1,IVTORP) = A3T59(I,IVTOR1,IVTORP)                 
     &                 + AAA * R2A9(I130,I)                                       
                  A4T5(I,IVTOR1,IVTORP) = A4T5(I,IVTOR1,IVTORP)                 
     &                 + AAA * R5B(I130,I)                                       
                  A4T5B(I,IVTOR1,IVTORP) = A4T5B(I,IVTOR1,IVTORP)               
     &                 + AAA * R5B2(I130,I)                                      
                  A5T5(I,IVTOR1,IVTORP) = A5T5(I,IVTOR1,IVTORP)                 
     &                 + AAA * R5C(I130,I)                                       
                  A5T5B(I,IVTOR1,IVTORP) = A5T5B(I,IVTOR1,IVTORP)               
     &                 + AAA * R5C2(I130,I)                                      
                  A4T5C(I,IVTOR1,IVTORP) = A4T5C(I,IVTOR1,IVTORP)               
     &                 + AAA * R5C2(I130,I)*R2A(I130,I)**2                       
                  A3T5C(I,IVTOR1,IVTORP) = A3T5C(I,IVTOR1,IVTORP)               
     &                 + AAA * R5B2(I130,I)*R2A(I130,I)**2
                  A4T5D(I,IVTOR1,IVTORP) = A4T5D(I,IVTOR1,IVTORP)               
     &                 + AAA * R5C2(I130,I)*R2A(I130,I)**4                       
                  A3T5D(I,IVTOR1,IVTORP) = A3T5D(I,IVTOR1,IVTORP)               
     &                 + AAA * R5B2(I130,I)*R2A(I130,I)**4
               end do

*     Pour VTORS                                                                    
                  do i130 = 1, ndimto
                  AAA = A(IVTORP,I130,M-1+I)*A(IVTOR1,I130,M-1+I)                 
                  A1TS(I,IVTOR1,IVTORP) = A1TS(I,IVTOR1,IVTORP)                 
     &                 + AAA * RVA(I130,I)                                       
                  A2TS(I,IVTOR1,IVTORP) = A2TS(I,IVTOR1,IVTORP)                 
     &                 + AAA * RVB(I130,I)  
               enddo

*     Pour VRTDI2                                                                                  
               do i140 = 1, ndimto-1
                  A4T2(I,IVTOR1,IVTORP) = A4T2(I,IVTOR1,IVTORP)                 
     &                 + A(IVTORP,I140+1,M-1+I)*A(IVTOR1,I140,M+I+1)              
                  A5T2(I,IVTOR1,IVTORP) = A5T2(I,IVTOR1,IVTORP)                 
     &                 + A(IVTORP,I140,M-1+I)*A(IVTOR1,I140+1,M+I+1)              

                  A4T2B(I,IVTOR1,IVTORP) = A4T2B(I,IVTOR1,IVTORP)                 
     &                 + A(IVTORP,I140+1,M-1+I)*A(IVTOR1,I140,M+I+1)              
     &                 * ( R2B3(I140, I) + R2B5(I140, I) +
     &                     R2B1(I140, I) + R2B2(I140, I) )
                  A5T2B(I,IVTOR1,IVTORP) = A5T2B(I,IVTOR1,IVTORP)                 
     &                 + A(IVTORP,I140,M-1+I)*A(IVTOR1,I140+1,M+I+1)              
     &                 * ( R2B3(I140, I) + R2B5(I140, I) + 
C     &                 * ( R2B4(I140, I) + R2B6(I140, I) + 
     &                     R2B1(I140, I) + R2B2(I140, I) )
               end do

*     Pour VRTDI4                   
               do i140 = 1, ndimto-1                                                
                  A2T4(I,IVTOR1,IVTORP) = A2T4(I,IVTOR1,IVTORP)                 
     &                 + A(IVTORP,I140+1,M-1+I)*A(IVTOR1,I140,M+I)                
                  A3T4(I,IVTOR1,IVTORP) = A3T4(I,IVTOR1,IVTORP)                 
     &                 + A(IVTORP,I140,M-1+I)*A(IVTOR1,I140+1,M+I)                

                  AK7T2(I,IVTOR1,IVTORP) = AK7T2(I,IVTOR1,IVTORP)  
     &                 + A(IVTORP, I140, M-1+I)*A(IVTOR1,I140+1,M+I)
     &                 * ( R2B1(i140,I) + R2B5(i140, I)
     &                 + R2B6(i140,I) + R2B4(i140, I))

                  AK7T3(I,IVTOR1,IVTORP) = AK7T2(I,IVTOR1,IVTORP)  
     &                 + A(IVTORP, I140+1, M-1+I)*A(IVTOR1,I140,M+I)
     &                 * ( R2B1(i140,I) + R2B5(i140, I)
     &                 + R2B6(i140,I) + R2B4(i140, I))

               end do

*     Pour VRTDI5                                                                   
               do i140 = 1, ndimto-1
                  AA5 = A(IVTORP,I140,M-1+I)*A(IVTOR1,I140+1,M-1+I)               
                  A1T5(I,IVTOR1,IVTORP) = A1T5(I,IVTOR1,IVTORP) + AA5           
                  A9T5(I,IVTOR1,IVTORP) = A9T5(I,IVTOR1,IVTORP)                 
     &                 + AA5 * R2A(I140,I)                                       
                  A9T59(I,IVTOR1,IVTORP) = A9T59(I,IVTOR1,IVTORP)                 
     &                 + AA5 * R2A9(I140,I)                                       
                  A11T5(I,IVTOR1,IVTORP) = A11T5(I,IVTOR1,IVTORP)               
     &                 + AA5 * R5A(I140,I)                                       
                  A11T58(I,IVTOR1,IVTORP) = A11T58(I,IVTOR1,IVTORP)               
     &                 + AA5 * R5A8(I140,I)                                       
                  A14T5(I,IVTOR1,IVTORP) = A14T5(I,IVTOR1,IVTORP)               
     &                 + AA5 * R5D(I140,I)                                       
                  A14T59(I,IVTOR1,IVTORP) = A14T59(I,IVTOR1,IVTORP)               
     &                 + AA5 * R5D9(I140,I)                                       
                  A16T5(I,IVTOR1,IVTORP) = A16T5(I,IVTOR1,IVTORP)               
     &                 + AA5 * R5E(I140,I)                             
                  A16T58(I,IVTOR1,IVTORP) = A16T58(I,IVTOR1,IVTORP)               
     &                 + AA5 * R5E8(I140,I)                             
               end do
               do i190 = 1, ndimto-2
                  A8T2(I,IVTOR1,IVTORP) = A8T2(I,IVTOR1,IVTORP)                 
     &                 + A(IVTORP,I190+2,M-1+I)*A(IVTOR1,I190,M+I+1)              
                  A9T2(I,IVTOR1,IVTORP) = A9T2(I,IVTOR1,IVTORP)                 
     &                 + A(IVTORP,I190,M-1+I)*A(IVTOR1,I190+2,M+I+1)              
                  A20T4(I,IVTOR1,IVTORP) = A20T4(I,IVTOR1,IVTORP)               
     &                 + A(IVTORP,I190+2,M-1+I)*A(IVTOR1,I190,M+I)                
                  A21T4(I,IVTOR1,IVTORP) = A21T4(I,IVTOR1,IVTORP)               
     &                 + A(IVTORP,I190,M-1+I)*A(IVTOR1,I190+2,M+I)                
               end do

*     Pour VRTDI4                                                                                  
               do i150 = 2, ndimto
                  A4T4(I,IVTOR1,IVTORP) = A4T4(I,IVTOR1,IVTORP)                 
     &                 + A(IVTORP,I150,M-1+I)*A(IVTOR1,I150-1,M+I)                
               end do

               
*     Pour VRTDI5                                                                   
               do i150 = 2, ndimto
                  AA5 = A(IVTORP,I150,M-1+I)*A(IVTOR1,I150-1,M-1+I)               
                  A7T5(I,IVTOR1,IVTORP) = A7T5(I,IVTOR1,IVTORP) + AA5           
                  A8T5(I,IVTOR1,IVTORP) = A8T5(I,IVTOR1,IVTORP)                 
     &                 + AA5 * R2A(I150,I)                                       
                  A8T59(I,IVTOR1,IVTORP) = A8T59(I,IVTOR1,IVTORP)                 
     &                 + AA5 * R2A9(I150,I)                                       
                  A10T5(I,IVTOR1,IVTORP) = A10T5(I,IVTOR1,IVTORP)               
     &                 + AA5 * R5A(I150,I)                                       
                  A10T58(I,IVTOR1,IVTORP) = A10T58(I,IVTOR1,IVTORP)               
     &                 + AA5 * R5A8(I150,I)                                       
                  A15T5(I,IVTOR1,IVTORP) = A15T5(I,IVTOR1,IVTORP)               
     &                 + AA5 * R5F(I150,I)                                       
                  A15T59(I,IVTOR1,IVTORP) = A15T59(I,IVTOR1,IVTORP)               
     &                 + AA5 * R5F9(I150,I)                                       
                  A17T5(I,IVTOR1,IVTORP) = A17T5(I,IVTOR1,IVTORP)               
     &                 + AA5 * R5G(I150,I)                                       
                  A17T58(I,IVTOR1,IVTORP) = A17T58(I,IVTOR1,IVTORP)               
     &                 + AA5 * R5G8(I150,I)                                       
               enddo
                  
*     Pour VRTDI5                                                                   
               do i160 = 1, ndimto-2
                  AA5B = A(IVTORP,I160,M-1+I)*A(IVTOR1,I160+2,M-1+I)               
                  A12T5(I,IVTOR1,IVTORP) = A12T5(I,IVTOR1,IVTORP)               
     &                 + A(IVTORP,I160,M-1+I)*A(IVTOR1,I160+2,M-1+I)              
                  A8T5B(I,IVTOR1,IVTORP) = A8T5B(I,IVTOR1,IVTORP)               
     &                 + AA5B* R2A(I160,I)                                       
                  A15T5B(I,IVTOR1,IVTORP) = A15T5B(I,IVTOR1,IVTORP)             
     &                 + AA5B* R5DD(I160,I)                                      

                  A8T5C(I,IVTOR1,IVTORP) = A8T5C(I,IVTOR1,IVTORP)               
     &                 + AA5B* R2A(I160,I)**2                                       
                  A15T5C(I,IVTOR1,IVTORP) = A15T5C(I,IVTOR1,IVTORP)             
     &                 + AA5B* R5DD(I160,I)**2                                      

               enddo

*     Pour VRTDI5                   
                  do i170 = 3, ndimto
                  AA5B = A(IVTORP,I170,M-1+I)*A(IVTOR1,I170-2,M-1+I)               
                     A13T5(I,IVTOR1,IVTORP) = A13T5(I,IVTOR1,IVTORP)               
     &                    + A(IVTORP,I170,M-1+I)*A(IVTOR1,I170-2,M-1+I)              
                  A9T5B(I,IVTOR1,IVTORP) = A9T5B(I,IVTOR1,IVTORP)               
     &                 + AA5B* R2A(I170,I)                                       
                  A14T5B(I,IVTOR1,IVTORP) = A14T5B(I,IVTOR1,IVTORP)             
     &                 + AA5B* R5FF(I170,I)                                      
                  A9T5C(I,IVTOR1,IVTORP) = A9T5C(I,IVTOR1,IVTORP)               
     &                 + AA5B* R2A(I170,I)**2                                       
                  A14T5C(I,IVTOR1,IVTORP) = A14T5C(I,IVTOR1,IVTORP)             
     &                 + AA5B* R5FF(I170,I)**2                                      
                  enddo
                  
*     Pour VTORS       
                  do i180 = 1, ndimto-3
                     A5TS(I,IVTOR1,IVTORP) = A5TS(I,IVTOR1,IVTORP)                 
     &                    + A(IVTORP,I180,M-1+I)*A(IVTOR1,I180+3,M-1+I)              
     &                    + A(IVTORP,I180+3,M-1+I)*A(IVTOR1,I180,M-1+I)
                  A22T4(I,IVTOR1,IVTORP) = A22T4(I,IVTOR1,IVTORP)               
     &                 + A(IVTORP,I180+3,M-1+I)*A(IVTOR1,I180,M+I)                
                  A23T4(I,IVTOR1,IVTORP) = A23T4(I,IVTOR1,IVTORP)               
     &                 + A(IVTORP,I180,M-1+I)*A(IVTOR1,I180+3,M+I)              
                  A12T2(I,IVTOR1,IVTORP) = A12T2(I,IVTOR1,IVTORP)                 
     &                 + A(IVTORP,I180+3,M-1+I)*A(IVTOR1,I180,M+I+1)              
                  A13T2(I,IVTOR1,IVTORP) = A13T2(I,IVTOR1,IVTORP)                 
     &                 + A(IVTORP,I180,M-1+I)*A(IVTOR1,I180+3,M+I+1)              
                  enddo

                  do i180 = 1, ndimto-4
                     A12TS(I,IVTOR1,IVTORP) = A12TS(I,IVTOR1,IVTORP)                 
     &                    + A(IVTORP,I180,M-1+I)*A(IVTOR1,I180+4,M-1+I)              
     &                    + A(IVTORP,I180+4,M-1+I)*A(IVTOR1,I180,M-1+I)
                  end do

               end do
            end do
         end do
         
      return 
      end
C     *************************************************************************  




C     *************************************************************************  
      SUBROUTINE VSET3(LO,V,CST)                                                 
C     *************************************************************************  
C     *************************************************************************  

      implicit none
      save

      include 'common_cst.f'
      include 'common_dim.f'
      include 'common_eigen.f'
      include 'common_eigen1.f'
      include 'common_tor.f'
      include 'common_dat.f'
      include 'common_rotor.f'
      include 'common_par.f'
      include 'common_big2.f'
      include 'common_chara.f'
      include 'common_expdat.f'
      include 'common_symb.f'
      include 'common_quant.f'
      include 'common_vset.f'

      double precision v, vrtdi1, vrtdi2, vrtdi4, vrtdi5, vtors, vrtdi6

      integer lo, l
      DIMENSION V(maxparam)           

      include 'common_vrtd.f'

      DO L=NA,NB                                                           
C     PRINT*,'PARA DANS VSET=',PARA(LO,L)
         V(L)=0.d0                                                              
         IF(PARA(LO,L).EQ.'OA    '.OR.PARA(LO,L).EQ.'DJ    '.OR.                
     &        PARA(LO,L).EQ.'DJK   '.OR.PARA(LO,L).EQ.'DK    '.OR.              
     &        PARA(LO,L).EQ.'HJ    '.OR.PARA(LO,L).EQ.'HJK   '.OR.              
     &        PARA(LO,L).EQ.'HKJ   '.OR.PARA(LO,L).EQ.'HK    '.OR.              
     &        PARA(LO,L).eq.'LJ    '.or.PARA(LO,L).eq.'BC4   '.or.
     &        PARA(LO,L).eq.'LK    '.or.PARA(LO,L).eq.'LJKKK '.or.
     &        PARA(LO,L).eq.'LJJKK '.or.PARA(LO,L).eq.'LJJJK '.or.
     &        PARA(LO,L).eq.'OBC4  ') THEN                                      
            V(L)=V(L)+VRTDI1(L,LO)                                              
C            IF(N.EQ.1) PRINT*,'V(L)=',V(L)
         ELSE IF(PARA(LO,L).EQ.'B     '.OR.PARA(LO,L).EQ.'C     ')THEN          
            V(L)=V(L)+VRTDI1(L,LO)+VRTDI2(L,LO)                                 
         ELSE IF(PARA(LO,L).EQ.'DAB   '.OR.PARA(LO,L).EQ.'ODAB  '.OR.           
     &           PARA(LO,L).EQ.'DABJ  '.OR.PARA(LO,L).EQ.'DABK  '.OR.           
     &           PARA(LO,L).EQ.'DAC   '.OR.PARA(LO,L).EQ.'ODELTA'.OR.           
     &           PARA(LO,L).EQ.'DELTAB'.OR.PARA(LO,L).EQ.'ODELTB'.OR.           
     &           PARA(LO,L).EQ.'DELTA '.OR.PARA(LO,L).EQ.'ODAB6 '.or.
     &           PARA(LO,L).eq.'HABK  '.or.PARA(LO,L).eq.'HABJ  '.or.
     &           PARA(LO,L).eq.'HABJK '.or.PARA(LO,L).eq.'ODABK '.or.
     &           para(lo,l).eq.'ODABJ '.or.para(lo,l).eq.'DELTAJ'.or.
     &           para(lo,l).eq.'DELTAK'.or.para(lo,l).eq.'DAGJ  '.or.
     &           para(lo,l).eq.'DAGK  '.or.para(lo,l).eq.'CABK  '.or.
     &           para(lo,l).eq.'CABJ  '.or.para(lo,l).eq.'DACK  '.or.
     &           para(lo,l).eq.'DACJ  '.or.para(lo,l).eq.'ODAC  '.or.
     &           para(lo,l).eq.'DAC6  '.or.para(lo,l).eq.'GB    '.or.
     &           para(lo,l).eq.'DELTJK'.or.para(lo,l).eq.'DELTKK'.or.
     &           para(lo,l).eq.'AK7AB '.or.para(lo,l).eq.'DACJJ '.or.
     &           para(lo,l).eq.'DACJK '.or.para(lo,l).eq.'DACKK '.or.
     &           para(lo,l).eq.'ODAB9 '.or.para(lo,l).eq.'DG4J  '.or.
     &           para(lo,l).eq.'DG4K  '.or.para(lo,l).eq.'ODABKK')THEN       
            V(L)=V(L)+VRTDI4(L,LO)                                              
         ELSE IF(PARA(LO,L).EQ.'ODELN '.OR.PARA(LO,L).EQ.'ODELK '.OR.           
     &           PARA(LO,L).EQ.'OHJ   '.OR.PARA(LO,L).EQ.'OHJK  '.OR.           
     &           PARA(LO,L).EQ.'OHK   '.OR.PARA(LO,L).EQ.'C1    '.OR.           
     &           PARA(LO,L).EQ.'C2    '.OR.PARA(LO,L).EQ.'C3    '.OR.           
     &           PARA(LO,L).EQ.'C11   '.OR.para(lo,l).eq.'C9BC  '.or.                                     
     &           PARA(LO,L).EQ.'C2J   '.OR.PARA(LO,L).EQ.'C1J   '.OR.           
     &           PARA(LO,L).EQ.'C1K   '.OR.PARA(LO,L).EQ.'C2K   '.OR.           
     &           PARA(LO,L).EQ.'C4K   '.OR.PARA(LO,L).EQ.'C4J   '.OR.           
     &           PARA(LO,L).EQ.'C12   '.OR.para(lo,l).eq.'CBCJ  '.or.
     &           para(lo,l).eq.'CBCK  '.or.para(lo,l).eq.'AK4BC '.or.
     &           para(lo,l).eq.'AK3BC '.or.para(lo,l).eq.'DBC6  '.or.
     &           PARA(LO,L).EQ.'DBC   '.OR.PARA(LO,L).EQ.'C4    '.or.
     &           PARA(LO,L).eq.'OLKKK '.or.PARA(LO,L).eq.'OLJKK '.or.
     &           PARA(LO,L).eq.'OLJJK '.or.PARA(LO,L).eq.'OLJJJ '.or.
     &           PARA(LO,L).eq.'AK7BC '.or.para(lo,l).eq.'DBCJ  '.or.
     &           PARA(LO,L).eq.'C6BCJ '.or.para(lo,l).eq.'AG3BCJ'.or.
     &           PARA(LO,L).eq.'G4BCJ '.or.parA(lo,l).eq.'C6BCK '.or.
     &           PARA(LO,L).eq.'AG3BCK'.or.PARA(LO,L).eq.'G4BCK '.or.
     &           PARA(LO,L).eq.'DBCK  '.or.para(lo,l).eq.'DBCJK '.or.
     &           PARA(LO,L).eq.'C2JK  '.or.para(lo,l).eq.'C2KK  '.or.
     &           PARA(LO,L).eq.'C1JK  '.or.para(lo,l).eq.'C1KK  '.or.
     &           PARA(LO,L).eq.'C4JK  '.or.para(lo,l).eq.'C4KK  '.or.
     &           PARA(LO,L).eq.'C6BCJK'.or.para(lo,l).eq.'C6BCKK'.or.
     &           PARA(LO,L).eq.'G3BCJK'.or.para(lo,l).eq.'G3BCKK'.or.
     &           PARA(LO,L).eq.'G4BCJK'.or.para(lo,l).eq.'G4BCKK') THEN         
            V(L)=V(L)+VRTDI2(L,LO)                                              
         ELSE IF(PARA(LO,L).EQ.'FV    '.OR.PARA(LO,L).EQ.'OFV   '.OR.           
     &           PARA(LO,L).EQ.'GV    '.OR.PARA(LO,L).EQ.'AK1   '.OR.           
     &           PARA(LO,L).EQ.'AK2   '.OR.PARA(LO,L).EQ.'AK3   '.OR.           
     &           PARA(LO,L).EQ.'AK4   '.OR.PARA(LO,L).EQ.'AK5   '.OR.           
     &           PARA(LO,L).EQ.'AK6   '.OR.PARA(LO,L).EQ.'AK7   '.OR.           
     &           PARA(LO,L).EQ.'AK6K  '.OR.PARA(LO,L).EQ.'AK7K  '.OR.           
     &           PARA(LO,L).EQ.'AK6KK '.OR.PARA(LO,L).EQ.'AK6JK '.OR.           
     &           PARA(LO,L).EQ.'AK6JJ '.OR.PARA(LO,L).EQ.'AK7KK '.OR.
     &           PARA(LO,L).EQ.'AK7JK '.OR.PARA(LO,L).EQ.'AK7JJ '.or.
     &           PARA(LO,L).EQ.'AK66  '.OR.PARA(LO,L).EQ.'AK6J  '.OR.           
     &           PARA(LO,L).EQ.'ALV   '.OR.PARA(LO,L).EQ.'OLV   '.OR.           
     &           PARA(LO,L).EQ.'AK5J  '.OR.PARA(LO,L).EQ.'GVJ   '.OR.           
     &           PARA(LO,L).EQ.'AK2J  '.OR.PARA(LO,L).EQ.'ANV   '.OR.           
     &           PARA(LO,L).EQ.'C6J4  '.OR.PARA(LO,L).EQ.'C6J2K2'.OR.           
     &           PARA(LO,L).EQ.'C6K4  '.OR.PARA(LO,L).EQ.'G4J4  '.OR.           
     &           PARA(LO,L).EQ.'G4J2K2'.OR.PARA(LO,L).EQ.'G4K4  '.OR.           
     &           PARA(LO,L).EQ.'V9J   '.OR.PARA(LO,L).EQ.'V9K   '.OR.           
     &           para(lo,l).eq.'V12J  '.or.para(lo,l).eq.'AK8   '.or.
     &           para(lo,l).eq.'AK9   '.or.para(lo,l).eq.'AK8J  '.or.
     &           para(lo,l).eq.'AK8K  '.or.para(lo,l).eq.'AK9J  '.or.
     &           para(lo,l).eq.'AK9K  '.or.para(lo,l).eq.'AK76  '.or.
     &           PARA(LO,L).EQ.'BK2   '.OR.PARA(LO,L).EQ.'AMV   '.OR.           
     &           PARA(LO,L).EQ.'BK1   '.OR.PARA(LO,L).EQ.'AK7J  '.OR.           
     &           PARA(LO,L).EQ.'AK1J  '.OR.PARA(LO,L).EQ.'AK5K  '.OR.           
     &           PARA(LO,L).EQ.'AK2K  '.OR.PARA(LO,L).EQ.'AK1K  '.OR.           
     &           PARA(LO,L).EQ.'AK3B  '.or.PARA(LO,L).eq.'AK3BK '.OR.
     &           PARA(LO,L).EQ.'AK4B  '.OR.PARA(LO,L).eq.'GAJ4  '.OR.
     &           PARA(LO,L).EQ.'GAJ2K2'.OR.PARA(LO,L).EQ.'GAK4  '.OR.
     &           PARA(LO,L).EQ.'AK3J  '.OR.PARA(LO,L).EQ.'AK3K  '.OR.
     &           PARA(LO,L).EQ.'GJ6   '.OR.PARA(LO,L).EQ.'GJ4K2 '.OR.
     &           PARA(LO,L).EQ.'GJ2K4 '.OR.PARA(LO,L).EQ.'GK6   '.or.
     &           PARA(LO,L).eq.'FJ6   '.or.PARA(LO,L).eq.'FJ4K2 '.or.
     &           PARA(LO,L).eq.'FJ2K4 '.or.PARA(LO,L).eq.'FK6   '.or.
     &           PARA(LO,L).EQ.'AGJ2K4'.OR.PARA(LO,L).EQ.'AGJ6  '.or.
     &           PARA(LO,L).EQ.'AGJ4K2'.OR.PARA(LO,L).EQ.'AGK6  '.or.
     &           para(lo,l).eq.'AK3BJ '.or.para(lo,l).eq.'AK4BJ '.or.
     &           para(lo,l).eq.'AK4BK '.or.para(lo,l).eq.'AK4C  '.or.
     &           para(lo,l).eq.'AK3C  '.or.para(lo,l).eq.'AK3CK4'.or.
     &           para(lo,l).eq.'AK3CK6'.or.para(lo,l).eq.'AK3BK4'.or.
     &           para(lo,l).eq.'AK3BK6'.or.para(lo,l).eq.'AK4BK4'.or.
     &           para(lo,l).eq.'AK4D  '.or.para(lo,l).eq.'AK3D  '.or.
     &           para(lo,l).eq.'AK4CK '.or.para(lo,l).eq.'AK3CK '.or.
     &           para(lo,l).eq.'AK4BK6'.or.para(lo,l).eq.'AK4CK4')
     &           then

C        PRINT*,'PARA BEFORE VRTDI5=',PARA(LO,L)
            V(L)=V(L)+VRTDI5(L,LO)                                              
C        PRINT*,'PARA AFTER VRTDI5=',PARA(LO,L),'LO=',LO,'L=',L
C        IF(N.EQ.1.AND.KK.EQ.1) THEN
C        PRINT*,'V=',V(L),'VRTDI5=',VRTDI5(L,LO),'KK=',KK,'IV=',IV
C    1   ,'LO=',LO,'L=',L,'ISIG=',ISIG
C        ENDIF
         ELSE IF(PARA(LO,L).EQ.'FPARA '.OR.PARA(LO,L).EQ.'V3    '.OR.           
     &           PARA(LO,L).EQ.'V6    '.OR.PARA(LO,L).EQ.'V9    '.OR.           
     &           PARA(LO,L).EQ.'RHORHO'.or.para(lo,l).eq.'V12   ') THEN                                   
            V(L)=V(L)+VTORS(L,CST,LO)                                           
C           PRINT*,' A= ',A(1,1,M),' M= ',M                                     
C        ELSE IF(PARA(LO,L).EQ.'RHO   '                                         
C     &     .OR.PARA(LO,L).EQ.'RHOJ  '.OR. PARA(LO,L).EQ.'RHOK  ')THEN          
         ELSE IF(PARA(LO,L).EQ.'RHO   ') THEN                                   
           V(L)=V(L)+VTORS(L,CST,LO)+VRTDI6(L,CST,LO,A(1,1,M))                  
         ENDIF                                                                  
      end do

      return 
      end
C     *************************************************************************  


C     *************************************************************************  
      subroutine vset(lo,v,cst)
C     ************************************************************************* 
C     Sets up the V matrix for each parameter to be fitted by taking the 
C     derivative of the energy. This subroutine is now only calculating the
C     derivatives and does not call ROTTOR (diagonal elements of the matrices
C     and energy calculation). Note that once the variables are initialized, 
C     this function calls vset1, vset2, and vset3. Originally these were
C     contained in the body of this function, but this function was split up
C     into more manageable chunks. Note that the bulk of the calculation is 
C     done in vset2.
C     *************************************************************************  
      implicit none
      save
C                  
      include 'common_cst.f'
      include 'common_dim.f'
      include 'common_eigen.f'
      include 'common_eigen1.f'
      include 'common_tor.f'
      include 'common_dat.f'
      include 'common_rotor.f'
      include 'common_par.f'
      include 'common_big2.f'
      include 'common_chara.f'
      include 'common_expdat.f'
      include 'common_symb.f'
      include 'common_quant.f'
      include 'common_vset.f'

      double precision v

      integer lo, ivtorp, ivtor1, i
      DIMENSION V(maxparam)           
*
      include 'common_vrtd.f'
*                                                                               
C     DEFINITION OF N,IV,KK,PR,ISIG CORRESPONDING TO THE CURRENT                
C     DATA K AND TO 1 STATE                                                     
C                                                                               
C     N=NOBS(LO,K)                                                              
C     IV=IVOBS(LO,K)                                                            
C     KK=KOBS(LO,K)                                                             
C     PR=PROBS(LO,K)                                                            
C     IF(K.GT.NADAT) THEN                                                       
C         ISIG=1                                                                
C     ELSE                                                                      
C        ISIG=0                                                                 
C     ENDIF                                                                     
*                                                                               
      IF(LO.EQ.1) THEN                                                          
         NI=1                                                                   
         M=1                                                                    
      ELSE                                                                      
         NI=NDMX+1                                                              
         M=KDMX+1                                                               
      END IF                                                                    
*                                                                               
      LOO=LO                                                                    
      IF(IRMW.EQ.1.AND.LO.EQ.2) LO=1                                            
      IF(LO.EQ.1) THEN                                                          
         NA=1                                                                   
         NB=NP(1)                                                               
      ELSE                                                                      
         NA=NP(1)+1                                                             
         NB=NP(1)+NP(2)                                                         
C        PRINT*,' NA= ',NA,' NB= ',NB                                           
      END IF                                                                    
*                                                                               
C     IBLK:DNUMBER OF THE N BLOCK                                               
C     FROM HERE,THE MATRIX EGVC,EGVL,A,ETOR WILL BE SPLITTED                    
C     INTO 2 PARTS,ONE FOR THE LOWER STATE(LO=1),ONE FOR THE                    
C     UPPER STATE(LO=2)                                                         
C     IF(IBLK(LOO,K).NE.IBLSTO(LOO)) THEN                                       
C        IBLSTO(LOO)=IBLK(LOO,K)                                                
C        WRITE(6,*) ' K = ',K,';  N = ',N                                       
C        CALL ROTTOR(CST,LO,EGVL(NI),EGVC(1,NI),ETOR(M,1)                       
C    1        ,A(1,1,M))                                                        
C     ENDIF                                                                     
C     NROTOR=(2*N+1)*7                                                          
C     CALL ENCAL(E,LOO,EGVL(NI),EGVC(1,NI))                                     
C     SUBROUTINE ENCAL TAKES THE GOOD EIGENVALUE AND EIGENVECT.                 
C     CORRESPONDING TO THE K CURRENT DATA                                       
C     DO LOOP ON THE PARAMETER L TO BE FITTED                                   
C                                                                               
      RHO=CST(LO,73)                                                        
      FPARA=CST(LO,72)                                                          
      RHORHO=CST(LO,77)                                                     

      DO IVTORP = 1,9                                                       
         DO IVTOR1 = 1,9                                                    
            DO I = 1, 2 * N + 1                                              
               AK7T1(I,IVTOR1,IVTORP) = 0.d0
               AK7T2(I,IVTOR1,IVTORP) = 0.d0
               AK7T3(I,IVTOR1,IVTORP) = 0.d0
               A1T1(I,IVTOR1,IVTORP) = 0.d0                                     
               A1T31(I,IVTOR1,IVTORP) = 0.d0
               A1T32(I,IVTOR1,IVTORP) = 0.d0
               A1T2(I,IVTOR1,IVTORP) = 0.d0                                     
               A2T2(I,IVTOR1,IVTORP) = 0.d0                                     
               A2T2K(I,IVTOR1,IVTORP) = 0.d0  
               A2T2KK(I,IVTOR1,IVTORP) = 0.d0                                     
               A3T2(I,IVTOR1,IVTORP) = 0.d0                                     
               A3T2K(I,IVTOR1,IVTORP) = 0.d0                                     
               A3T2KK(I,IVTOR1,IVTORP) = 0.d0                                     
               A4T2(I,IVTOR1,IVTORP) = 0.d0                                     
               A5T2(I,IVTOR1,IVTORP) = 0.d0                                     
               A4T2B(I,IVTOR1,IVTORP) = 0.d0
               A5T2B(I,IVTOR1,IVTORP) = 0.d0
               A7T2(I,IVTOR1,IVTORP) = 0.d0                                     
               A8T2(I,IVTOR1,IVTORP) = 0.d0                                     
               A9T2(I,IVTOR1,IVTORP) = 0.d0                                     
               A12T2(I,IVTOR1,IVTORP) = 0.d0                                     
               A13T2(I,IVTOR1,IVTORP) = 0.d0                                     
               A1T4(I,IVTOR1,IVTORP) = 0.d0                                     
               A1T4B(I,IVTOR1,IVTORP) = 0.d0
               A2T4(I,IVTOR1,IVTORP) = 0.d0                                     
               A3T4(I,IVTOR1,IVTORP) = 0.d0                                     
               A4T4(I,IVTOR1,IVTORP) = 0.d0                                     
               A5T4(I,IVTOR1,IVTORP) = 0.d0                                     
               A6T4(I,IVTOR1,IVTORP) = 0.d0                                     
               A7T4(I,IVTOR1,IVTORP) = 0.d0                                     
               A8T4(I,IVTOR1,IVTORP) = 0.d0
               A9T4(I,IVTOR1,IVTORP) = 0.d0
               A20T4(I,IVTOR1,IVTORP) = 0.d0                                    
               A21T4(I,IVTOR1,IVTORP) = 0.d0                                    
               A22T4(I,IVTOR1,IVTORP) = 0.d0                                    
               A23T4(I,IVTOR1,IVTORP) = 0.d0                                    
               A1T5(I,IVTOR1,IVTORP) = 0.d0                                     
               A2T5(I,IVTOR1,IVTORP) = 0.d0                                     
               A2T58(I,IVTOR1,IVTORP) = 0.d0                                     
               A3T5(I,IVTOR1,IVTORP) = 0.d0                                     
               A3T5(I,IVTOR1,IVTORP) = 0.d0                                     
               A4T5(I,IVTOR1,IVTORP) = 0.d0                                     
               A4T5B(I,IVTOR1,IVTORP) = 0.d0                                    
               A5T5(I,IVTOR1,IVTORP) = 0.d0                                     
               A5T5B(I,IVTOR1,IVTORP) = 0.d0                                    
               A7T5(I,IVTOR1,IVTORP) = 0.d0                                     
               A8T5(I,IVTOR1,IVTORP) = 0.d0                                     
               A8T59(I,IVTOR1,IVTORP) = 0.d0                                     
               A8T5B(I,IVTOR1,IVTORP) = 0.d0                                    
               A8T5C(I,IVTOR1,IVTORP) = 0.d0                                    
               A9T5(I,IVTOR1,IVTORP) = 0.d0                                     
               A9T59(I,IVTOR1,IVTORP) = 0.d0                                     
               A9T5B(I,IVTOR1,IVTORP) = 0.d0                                    
               A9T5C(I,IVTOR1,IVTORP) = 0.d0                                    
               A10T5(I,IVTOR1,IVTORP) = 0.d0                                   
               A11T5(I,IVTOR1,IVTORP) = 0.d0                                    
               A10T58(I,IVTOR1,IVTORP) = 0.d0                                   
               A11T58(I,IVTOR1,IVTORP) = 0.d0                                    
               A12T5(I,IVTOR1,IVTORP) = 0.d0                                    
               A13T5(I,IVTOR1,IVTORP) = 0.d0                                    
               A14T5(I,IVTOR1,IVTORP) = 0.d0                                    
               A14T59(I,IVTOR1,IVTORP) = 0.d0                                    
               A14T5B(I,IVTOR1,IVTORP) = 0.d0                                   
               A14T5C(I,IVTOR1,IVTORP) = 0.d0                                    
               A15T5(I,IVTOR1,IVTORP) = 0.d0                                    
               A15T59(I,IVTOR1,IVTORP) = 0.d0                                    
               A15T5B(I,IVTOR1,IVTORP) = 0.d0                                   
               A15T5C(I,IVTOR1,IVTORP) = 0.d0                                    
               A16T5(I,IVTOR1,IVTORP) = 0.d0                                    
               A17T5(I,IVTOR1,IVTORP) = 0.d0                                    
               A16T58(I,IVTOR1,IVTORP) = 0.d0                                    
               A17T58(I,IVTOR1,IVTORP) = 0.d0                                    
               A1TS(I,IVTOR1,IVTORP) = 0.d0                                     
               A2TS(I,IVTOR1,IVTORP) = 0.d0                                     
               A5TS(I,IVTOR1,IVTORP) = 0.d0                                     
               A12TS(I,IVTOR1,IVTORP) = 0.d0
               A5T4B(I,IVTOR1,IVTORP) = 0.d0                                    
               A6T4B(I,IVTOR1,IVTORP) = 0.d0                                    
               A6T4C(I,IVTOR1,IVTORP) = 0.d0                                    
               A10T2(I,IVTOR1,IVTORP) = 0.d0                                    
               A10T2K(I,IVTOR1,IVTORP) = 0.d0                                    
               A10T2KK(I,IVTOR1,IVTORP) = 0.d0                                    
               A11T2(I,IVTOR1,IVTORP) = 0.d0                                    
               A11T2K(I,IVTOR1,IVTORP) = 0.d0                                    
               A4T5C(I,IVTOR1,IVTORP) = 0.d0
               A3T5C(I,IVTOR1,IVTORP) = 0.d0
               A4T5D(I,IVTOR1,IVTORP) = 0.d0
               A3T5D(I,IVTOR1,IVTORP) = 0.d0
            end do
         end do
      end do

      call vset1(lo, v, cst)
      call vset2
      call vset3(lo, v, cst)
                                                         
      end                                                                       
C     *************************************************************************  


