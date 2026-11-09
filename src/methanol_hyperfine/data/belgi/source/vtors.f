C*************************************************************                  
      FUNCTION VTORS (L, CST, LO)                                               
      implicit none
      save
C*************************************************************                  
C     DERIVATIVES OF THE MATRIX ELEMENTS FOR THE TOSION HAMILTONIAN             
C     F,V3,V6,V9                                                                
C                  

      include 'common_dim.f'

      include 'common_eigen.f'
      include 'common_chara.f'
      include 'common_quant.f'
      include 'common_vrtd.f'
      include 'common_cst.f'
      
      double precision vtors, vtorsa, vtorsb, vtorsc, vtorsd, 
     &     vtorse, fpara
      integer lo, ivtorp, ivtor1, i, l
      double precision EGT(KDMX,9,9)                                                      

      vtors = 0.0d0
*                                                                               
      FPARA=CST(LO,72)                                                          
*                                                                               
      DO IVTORP = 1,9                                                     
         DO IVTOR1 = 1,9                                                  
            DO I = 1, 2 * N + 1                                            
               EGT(I,IVTOR1,IVTORP) = EGVC3(I,IVTORP) * EGVC3(I,IVTOR1)         
            end do
         end do
      end do
*                                                                               
      IF (PARA (LO, L) .EQ. 'RHORHO') THEN                                      
       VTORSA = 0.d0                                                            
         DO IVTORP = 1,9                                                  
            DO IVTOR1 = 1,9                                               
               DO I = 1, 2 * N + 1                                         
                VTORSA = VTORSA + EGT(I,IVTOR1,IVTORP)                          
     &                 * (I - N - 1) * A1TS(I,IVTOR1,IVTORP)                    
             end do
          end do
       end do
         VTORS = VTORSA * 2.d0 * FPARA                                              
      ELSE IF (PARA (LO, L) .EQ. 'FPARA ') THEN                                 
         VTORSB = 0.d0                                                          
         DO IVTORP = 1,9                                                  
            DO IVTOR1 = 1,9                                               
               DO I = 1, 2 * N + 1                                         
                  VTORSB = VTORSB + EGT(I,IVTOR1,IVTORP)                        
     &                 * A2TS(I,IVTOR1,IVTORP)                                  
               end do
            end do
         end do
         VTORS = VTORSB                                                         
      ELSE IF (PARA (LO, L) .EQ. 'V3    ') THEN                                 
         VTORSC = 0.d0                                                          
         DO IVTORP = 1,9                                                  
            DO IVTOR1 = 1,9                                               
               DO I = 1, 2 * N + 1                                         
C Jonathan
C                  VTORSC = VTORSC + EGT(I,IVTOR1,IVTORP)                        
C*    &                 * (A1T1(I,IVTOR1,IVTORP) - A4TS(I,IVTOR1,IVTORP))        
C     &                 * (A1T1(I,IVTOR1,IVTORP) - A1T5(I,IVTOR1,IVTORP))        

                  VTORSC = VTORSC + EGT(I,IVTOR1,IVTORP)*
     &             (A1T1(I,IVTOR1,IVTORP)
     &              -0.5d0*A1T5(I,IVTOR1,IVTORP)
     &              -0.5d0*A7T5(I,IVTOR1,IVTORP))

               end do
            end do
         end do
         VTORS = VTORSC / 2.d0                                                  
      ELSE IF (PARA (LO, L) .EQ. 'V6    ') THEN                                 
         VTORSD = 0.d0                                                          
         DO IVTORP = 1,9                                                  
            DO IVTOR1 = 1,9                                               
               DO I = 1, 2 * N + 1                            
C JONATHAN
C                  VTORSD = VTORSD + EGT(I,IVTOR1,IVTORP)                        
C     &                * (A1T1(I,IVTOR1,IVTORP)  - A12T5(I,IVTOR1,IVTORP)
                  VTORSD = VTORSD + EGT(I,IVTOR1,IVTORP)                        
     &                * (A1T1(I,IVTOR1,IVTORP) 
     &                 - 0.5d0*A12T5(I,IVTOR1,IVTORP)
     &                 - 0.5d0*A13T5(I,IVTOR1,IVTORP))        
               end do
            end do
         end do
         VTORS = VTORSD / 2.d0                                                  
      ELSE IF (PARA (LO, L) .EQ. 'V9    ') THEN                                 
         VTORSE = 0.d0                                                          
         DO IVTORP = 1,9                                                  
            DO IVTOR1 = 1,9                                               
               DO I = 1, 2 * N + 1                                         
C Jonathan
C                  VTORSE = VTORSE + EGT(I,IVTOR1,IVTORP)                        
C     &                 * (A1T1(I,IVTOR1,IVTORP) - A5TS(I,IVTOR1,IVTORP))        
                  VTORSE = VTORSE + EGT(I,IVTOR1,IVTORP)                        
     &                 *(A1T1(I,IVTOR1,IVTORP)
     &                 - 0.5d0*A5TS(I,IVTOR1,IVTORP))
               end do
            end do
         end do
         VTORS = VTORSE / 2.d0                                                  
      ELSE IF (PARA (LO, L) .EQ. 'V12   ') THEN                                 
         VTORSE = 0.d0                                                          
         DO IVTORP = 1,9                                                  
            DO IVTOR1 = 1,9                                               
               DO I = 1, 2 * N + 1                                         
                  VTORSE = VTORSE + EGT(I,IVTOR1,IVTORP)                        
     &                 * (A1T1(I,IVTOR1,IVTORP) 
C     &                 - A12TS(I,IVTOR1,IVTORP))        
     &                 - 0.5d0*A12TS(I,IVTOR1,IVTORP))        
               end do
            end do
         end do
         VTORS = VTORSE / 2.d0                                                  
      END IF                                                                    
      END                                                                       
