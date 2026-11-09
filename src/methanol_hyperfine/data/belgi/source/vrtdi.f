C***********************************************************                    
      FUNCTION VRTDI1 (L, LO)                                                   
      implicit none
      save
C***********************************************************                    
C     DERIVATIVES OF THE ENERGY INVOLVING THE DIAGONAL ELEMENTS OF              
C     HAMILTONIAN(PARAMETERS:DA,B,C,DJ,DJK,DK,HJ,HJK,HKJ,HK)                    
C         

      include 'common_dim.f'
      include 'common_eigen.f'
      include 'common_chara.f'
      include 'common_quant.f'
      include 'common_vrtd.f'

      double precision vrtdi1, vrtdia, vrtdib, vrtdic, vrtdid
      double precision vrtdie, pm4fact, k1, k2, k3, k4
      integer ivtorp, ivtor1, lo, k, i, l
*                                                                               
      double precision EGTA(KDMX,9,9)                                                     
      double precision EGT(KDMX,9,9)
*                          
      xn = dble(n)
      xnn = dble(n*(n+1))

      vrtdi1 = 0.0d0

      DO IVTORP = 1,9                                                     
         DO IVTOR1 = 1,9                                                  
            DO I = 1, 2 * N + 1                                            
               EGTA(I,IVTOR1,IVTORP) = EGVC3(I,IVTORP) * EGVC3(I,IVTOR1)        
     &                                           * A1T1(I,IVTOR1,IVTORP)        

               EGT(I,IVTOR1,IVTORP) = EGVC3(I,IVTORP)*EGVC3(I+4,IVTOR1)        
            end do
         end do
      end do

C      if (PARA(LO, L) .eq. 'LJ    ') then
C         vrtdia = 0.0d0
C      end if
*                                                                               
      IF (PARA(LO, L).EQ.'B     ' .OR. PARA(LO, L).EQ.'C     ' .OR.             
     &    PARA(LO, L).EQ.'DJ    ' .OR. PARA(LO, L).EQ.'HJ    ' .or. 
     &    PARA(LO, L).eq.'LJ    ') THEN            
         VRTDIA = 0.d0                                                          
         DO IVTORP = 1,9                                                  
            DO IVTOR1 = 1,9                                               
               DO I = 1, 2 * N + 1                                         
                  VRTDIA = VRTDIA + EGTA(I,IVTOR1,IVTORP)                       
               end do
            end do
         end do
      END IF                                                                    
*                                                                               
      IF (PARA(LO, L).EQ.'OA    ' .OR. PARA(LO, L).EQ.'B     ' .OR.             
     &    PARA(LO, L).EQ.'C     ' .OR. PARA(LO, L).EQ.'DJK   ' .OR.             
     &    PARA(LO, L).EQ.'HJK   ' .or. PARA(LO, L).eq.'LJJJK ') THEN                                         
         VRTDIB = 0.d0                                                          
         DO IVTORP = 1,9                                                  
            DO IVTOR1 = 1,9                                               
               DO I = 1, 2 * N + 1                                         
                  K = I - N - 1                                                 
                  XK = dble(K)**2
                  VRTDIB = VRTDIB + EGTA(I,IVTOR1,IVTORP) * XK
               end do
            end do
         end do
      END IF                                                                    
*                           

      IF (PARA(LO,L).EQ.'BC4   ') then
         VRTDIB = 0.d0                                                          
         DO IVTORP = 1,9                                                  
            DO IVTOR1 = 1,9                                               
               DO I = 1, 2 * N - 1                                         
                  K = I - N - 1                                                 
                  VRTDIB = VRTDIB + EGTA(I,IVTOR1,IVTORP)                        
     &                 * dble(N*(N+1)-K*(K-1))*dble(N*(N+1)-(K-1)*(K-2))
     &                 * A1T1(I, IVTOR1, IVTORP)

                  VRTDIB = VRTDIB + EGTA(I,IVTOR1,IVTORP)                        
     &                 * dble(N*(N+1)-K*(K+1))*dble(N*(N+1)-(K+1)*(K+2))
     &                 * A1T1(I, IVTOR1, IVTORP)

                  k1 = dble(N*(N+1)-K*(K+1))
                  k2 = dble(N*(N+1)-(K+1)*(K+2))
                  k3 = dble(N*(N+1)-(K+2)*(K+3))
                  k4 = dble(N*(N+1)-(K+3)*(K+4))

                  pm4fact = k1*k2*k3*k4

C                  if ( pm4fact .gt. 0.0 ) pm4fact = SQRT(pm4fact)
C                  if ( pm4fact .lt. 0.0 ) pm4fact = 0.0

                  if ( pm4fact .gt. 0.0 ) pm4fact = 2.0d0*SQRT(pm4fact)
                  if ( pm4fact .lt. 0.0 ) pm4fact = 0.0

                  VRTDIB = VRTDIB + pm4fact * EGT(I, IVTOR1, IVTORP)
     &                 * A1T31(I, IVTOR1, IVTORP)

C                  k1 = dble(N*(N+1)-K*(K-1))
C                  k2 = dble(N*(N+1)-(K-1)*(K-2))
C                  k3 = dble(N*(N+1)-(K-2)*(K-3))
C                  k4 = dble(N*(N+1)-(K-3)*(K-4))

C                  pm4fact = k1*k2*k3*k4

C                  if ( pm4fact .gt. 0.0 ) pm4fact = SQRT(pm4fact)
C                  if ( pm4fact .lt. 0.0 ) pm4fact = 0.0

C                  VRTDIB = VRTDIB + pm4fact * EGT(I, IVTOR1, IVTORP)
C     &                 * A1T32(I, IVTOR1, IVTORP)

               end do
            end do
         end do
      END IF                                                                    

      IF (PARA(LO,L).EQ.'OBC4  ') then
         VRTDIB = 0.d0                                                          
         DO IVTORP = 1,9                                                  
            DO IVTOR1 = 1,9                                               
               DO I = 1, 2 * N - 1                                         
                  K = I - N - 1                                                 

                  k1 = dble(N*(N+1)-K*(K+1))
                  k2 = dble(N*(N+1)-(K+1)*(K+2))
                  k3 = dble(N*(N+1)-(K+2)*(K+3))
                  k4 = dble(N*(N+1)-(K+3)*(K+4))

                  pm4fact = k1*k2*k3*k4

                  if ( pm4fact .gt. 0.0 ) pm4fact = 2.0d0*SQRT(pm4fact)
                  if ( pm4fact .lt. 0.0 ) pm4fact = 0.0

                  VRTDIB = VRTDIB + pm4fact * EGT(I, IVTOR1, IVTORP)
     &                 * A1T31(I, IVTOR1, IVTORP)

               end do
            end do
         end do
      END IF                                                                    
                                                    
      IF (PARA(LO, L).EQ.'DK    ' .OR. PARA(LO, L).EQ.'HKJ   ' .or.
     &    PARA(LO, L).eq.'LJJKK ') THEN            
         VRTDIC = 0.d0                                                          
         DO IVTORP = 1,9                                                  
            DO IVTOR1 = 1,9                                               
               DO I = 1, 2 * N + 1                                         
                  K = I - N - 1                               
                  XK = dble(K)**4
                  VRTDIC = VRTDIC + EGTA(I,IVTOR1,IVTORP) * XK
               end do
            end do
         end do
      END IF                                                                    
*                                                                               
      IF (PARA(LO, L).EQ.'HK    ' .OR.
     &    PARA(LO, L).eq.'LJKKK ') THEN            
         VRTDID = 0.d0                                                          
         DO IVTORP = 1,9                                                  
            DO IVTOR1 = 1,9                                               
               DO I = 1, 2 * N + 1                                         
                  K = I - N - 1       
                  xk = dble(k)**6
                  VRTDID = VRTDID + EGTA(I,IVTOR1,IVTORP) * xk
               end do
            end do
         end do
      END IF                                                                    

      IF (PARA(LO, L).EQ.'LK    ') then
         VRTDIe = 0.d0                                                          
         DO  IVTORP = 1,9                                                  
            DO  IVTOR1 = 1,9                                               
               DO  I = 1, 2 * N + 1                                         
                  K = I - N - 1           
                  xk = dble(k)**8
                  VRTDIe = VRTDIe + EGTA(I,IVTOR1,IVTORP) *
     &                       xk
                  end do
               end do
            end do
      END IF                                                                    

*                                                                               
      IF (PARA(LO, L).EQ. 'OA    ') THEN                                        
         VRTDI1 = VRTDIB                                                        
C         IF(N.EQ.1) PRINT*,'VRTDI1=',VRTDI1
      ELSE IF (PARA(LO, L).EQ.'B     '.OR.PARA(LO, L).EQ.'C     ') THEN         
         VRTDI1 = 0.5d0 * (VRTDIA * N*(N+1) - VRTDIB)                             
      ELSE IF (PARA(LO, L).EQ.'DJ    ') THEN                                    
         VRTDI1 = - VRTDIA * xnn**2                                    
      ELSE IF (PARA(LO, L).EQ.'DJK   ') THEN                                    
         VRTDI1 = - VRTDIB * xnn                                            
      ELSE IF (PARA(LO, L).EQ.'DK    ') THEN                                    
         VRTDI1 = - VRTDIC                                                      
      ELSE IF (PARA(LO, L).EQ.'HJ    ') THEN                                    
         VRTDI1 = VRTDIA * (xnn) ** 3                                       
      ELSE IF (PARA (LO, L) .EQ. 'HJK   ') THEN                                 
         VRTDI1 = VRTDIB * xnn**2                                      
      ELSE IF (PARA (LO, L) .EQ. 'HKJ   ') THEN                                 
         VRTDI1 = VRTDIC * xnn                                              
      ELSE IF (PARA (LO, L) .EQ. 'HK    ') THEN                                 
         VRTDI1 = VRTDID                                                        
      else if (PARA(LO, L) .eq. 'LJ    ') THEN
         vrtdi1 = vrtdia * (xnn) ** 4
      else if (PARA(LO,L) .eq. 'LK    ') then
         vrtdi1 = vrtdie
      else if (PARA(LO,L) .eq. 'LJKKK ') then
         vrtdi1 = vrtdid* xnn   
      else if (PARA(LO,L) .eq. 'LJJKK ') then
         vrtdi1 = vrtdic * (xnn)**2
      else if (PARA(LO,L) .eq. 'LJJJK ') then
         vrtdi1 = vrtdib * (xnn)**3
      ELSE IF (PARA(LO,L).EQ.'BC4   ') THEN                                     
         VRTDI1 = VRTDIB * 0.25d0                                              
      ELSE IF (PARA(LO,L).EQ.'OBC4  ') THEN                                     
         VRTDI1 = VRTDIB * 0.25d0                                              
      END IF                         
                                           

      END                                                                       


C**********************************************************                     
      FUNCTION VRTDI2 (L, LO)                                                   
      implicit none
      save
C************************************************************                   
C     DERIVATIVES INVOLVING THE K/K+ OU - 2 MATRIX ELEMENTS OF                  
C     THE HAMILTONIAN(B,C,ODELN,ODELK,OHJ,OHJK,OHK,C1,C2,C3)                    
C         

      include 'common_dim.f'
      include 'common_eigen.f'
      include 'common_quant.f'
      include 'common_chara.f'
      include 'common_vrtd.f'

      double precision vrtdi2, vrtdia, vrtdib, vrtdic, vrtdid
      double precision vrtdii, vrtdie, vrtdij, vrtdik, vrtdif
      double precision vrtdig, vrtdih, vrtdil, vrtdim, vrtdin
      integer ivtorp, ivtor1, lo, i, l, k
*                                                                               
      double precision EGT(KDMX,9,9),SQT(KDMX)                                            
*                                                                               
      xn = dble(n)
      xnn = dble(n*(n+1))

      vrtdi2 = 0.0d0
      DO I = 1, 2 * N - 1                                                   
         SQT(I) = 2. * SQRT (dble((N*(N+1) - (I-N-1)*(I-N))                   
     &                         * (N*(N+1) - (I-N)*(I-N+1))))                    
      end do
      DO IVTORP = 1,9                                                     
         DO IVTOR1 = 1,9                                                  
            DO I = 1, 2 * N - 1                                            
               EGT(I,IVTOR1,IVTORP) = SQT(I) * EGVC3(I,IVTORP)                  
     &                 * EGVC3(I+2,IVTOR1)                                      
            end do
         end do
      end do
*                                                                               
      IF (PARA(LO,L).EQ.'B     ' .OR. PARA(LO,L).EQ.'C     ' .OR.               
     &    PARA(LO,L).EQ.'ODELN ' .OR. PARA(LO,L).EQ.'OHJ   '.or.
     &    PARA(LO,L).eq.'OLJJJ ') THEN              
         VRTDIA = 0.d0                                                          
         DO IVTORP = 1,9                                                  
            DO IVTOR1 = 1,9                                               
               DO I = 1, 2 * N - 1                                         
                  VRTDIA = VRTDIA + EGT(I,IVTOR1,IVTORP)                        
     &                 * A1T2(I,IVTOR1,IVTORP)                                  
               end do
            end do
         end do
      END IF                                                                    
*                                                                               
      IF (PARA(LO,L).EQ.'ODELK '.OR. PARA(LO,L).EQ.'OHJK  ' .or.
     &    PARA(LO,L).eq.'OLJJK ') THEN               
         VRTDIB = 0.d0                                                          
         DO IVTORP = 1,9                                                  
            DO IVTOR1 = 1,9                                               
               DO I = 1, 2 * N - 1                                         
                  K = I - N - 1                                                 
                  VRTDIB = VRTDIB + EGT(I,IVTOR1,IVTORP)                        
     &                 * (dble(K*K) + dble((K+2)*(K+2))) * 
     &                 A1T2(I,IVTOR1,IVTORP)            

               end do
            end do
         end do
      END IF                                                                    

      IF (PARA(LO,L).EQ.'OLKKK ') THEN               
         VRTDIn = 0.d0                                                          
         DO  IVTORP = 1,9                                                  
            DO  IVTOR1 = 1,9                                               
               DO  I = 1, 2 * N - 1                                         
                 K = I - N - 1      
                 xk = dble(k)
                  VRTDIn = VRTDIn + EGT(I,IVTOR1,IVTORP)                        
     &                 * (xK**6 + (xk+2)**6) * A1T2(I,IVTOR1,IVTORP)            
               end do
            end do
         end do
      END IF                  
*                                                                               
      IF (PARA(LO,L).EQ.'OHK   '.or. PARA(LO,L).eq.'OLJKK ') THEN                                          
         VRTDIC = 0.d0                                                          
         DO IVTORP = 1,9                                                  
            DO IVTOR1 = 1,9                                               
               DO I = 1, 2 * N - 1                                         
                  K = I - N - 1      
                  xk = dble(k)
                  VRTDIC = VRTDIC + EGT(I,IVTOR1,IVTORP)                        
     &                 * (xK**4 + (xK+2)**4) * A1T2(I,IVTOR1,IVTORP)              
               end do
            end do
         end do
      END IF                                                                    
*                                                                               
      IF (PARA(LO,L).EQ.'C1    ' .OR. PARA(LO,L).EQ.'C1J   ') THEN              
         VRTDID = 0.d0                                                          
         DO IVTORP = 1,9                                                  
            DO IVTOR1 = 1,9                                               
               DO I = 1, 2 * N - 1                                         
                  VRTDID = VRTDID + EGT(I,IVTOR1,IVTORP)                        
     &                 * A2T2(I,IVTOR1,IVTORP)                                  
               end do
            end do
         end do
      END IF                                                                    
*                                                                               
      IF (PARA(LO,L).EQ.'C1K   ')THEN                                           
         VRTDII = 0.d0                                                          
         DO IVTORP = 1,9                                                  
            DO IVTOR1 = 1,9                                               
               DO I = 1, 2 * N - 1                                         
               K=I-N-1
                  VRTDII = VRTDII + EGT(I,IVTOR1,IVTORP)                        
     &                 * A2T2K(I,IVTOR1,IVTORP)                                  
               end do
            end do
         end do
      END IF                                                                    

      IF (PARA(LO,L).EQ.'C1JK  ')THEN                                           
         VRTDII = 0.d0                                                          
         DO IVTORP = 1,9                                                  
            DO IVTOR1 = 1,9                                               
               DO I = 1, 2 * N - 1                                         
               K=I-N-1
                  VRTDII = VRTDII + EGT(I,IVTOR1,IVTORP)                        
     &                 * A2T2K(I,IVTOR1,IVTORP)*xnn
               end do
            end do
         end do
      END IF                                                                    

      IF (PARA(LO,L).EQ.'C1KK  ')THEN                                           
         VRTDII = 0.d0                                                          
         DO IVTORP = 1,9                                                  
            DO IVTOR1 = 1,9                                               
               DO I = 1, 2 * N - 1                                         
               K=I-N-1
                  VRTDII = VRTDII + EGT(I,IVTOR1,IVTORP)                        
     &                 * A2T2KK(I,IVTOR1,IVTORP)                                  
               end do
            end do
         end do
      END IF                                                                    
*
      IF (PARA(LO,L).EQ.'C2    ' .OR. PARA(LO,L).EQ.'C2J   ') THEN              
         VRTDIE = 0.d0                                                          
         DO IVTORP = 1,9                                                  
            DO IVTOR1 = 1,9                                               
               DO I = 1, 2 * N - 1                                         
                  VRTDIE = VRTDIE + EGT(I,IVTOR1,IVTORP)                        
     &                 * (A1T2(I,IVTOR1,IVTORP)                                 
     &                     - 0.5d0 * (A4T2(I,IVTOR1,IVTORP)                     
     &                              + A5T2(I,IVTOR1,IVTORP)))                   
               end do
            end do
         end do
      END IF                                                                    

      IF (PARA(LO,L).EQ.'AK4BC ') THEN              
         VRTDIE = 0.d0                                                          
         DO IVTORP = 1,9                                                  
            DO IVTOR1 = 1,9                                               
               DO I = 1, 2 * N - 1                                         
                  VRTDIE = VRTDIE + EGT(I,IVTOR1,IVTORP)                        
     &                 * A2TG1(I, IVTOR1, IVTORP)

               end do
            end do
         end do
      END IF                                                                    

      IF (PARA(LO,L).EQ.'AK3BC ') THEN              
         VRTDIE = 0.d0                                                          
         DO IVTORP = 1,9                                                  
            DO IVTOR1 = 1,9                                               
               DO I = 1, 2 * N - 1                                         
                  VRTDIE = VRTDIE + EGT(I,IVTOR1,IVTORP)                        
     &                 * A2TG2(I,IVTOR1,IVTORP)

               end do
            end do
         end do
      END IF                                                                    

      IF (PARA(LO,L).EQ.'AK7BC ') THEN              
         VRTDIE = 0.d0                                                          
         DO  IVTORP = 1,9                                                  
            DO  IVTOR1 = 1,9                                               
               DO I = 1, 2 * N - 1                                         
                  VRTDIE = VRTDIE + EGT(I,IVTOR1,IVTORP)                        
     &                 * (2.d0 * A2T2(I,IVTOR1,IVTORP)                                 
     &                     - 0.5d0 * (A4T2B(I,IVTOR1,IVTORP)                     
     &                              + A5T2B(I,IVTOR1,IVTORP)))                   
               end do
            end do
         end do
      END IF                                                                    
*                                                                               
      IF (PARA(LO,L).EQ.'C2K   ') THEN                                          
         VRTDIJ = 0.d0                                                          
         DO IVTORP = 1,9                                                  
            DO IVTOR1 = 1,9                                               
               DO I = 1, 2 * N - 1                                         
               K=I-N-1
               xk = dble(k)
                  VRTDIJ = VRTDIJ + EGT(I,IVTOR1,IVTORP)                        
     &                 * (A1T2(I,IVTOR1,IVTORP)                                 
     &                     - 0.5d0 * (A4T2(I,IVTOR1,IVTORP)                     
     &                              + A5T2(I,IVTOR1,IVTORP)))                   
     &                 *(xK**2+(xK+2)**2)
               end do
            end do
         end do
      END IF                                                                    

      IF (PARA(LO,L).EQ.'C2JK  ') THEN                                          
         VRTDIJ = 0.d0                                                          
         DO IVTORP = 1,9                                                  
            DO IVTOR1 = 1,9                                               
               DO I = 1, 2 * N - 1                                         
               K=I-N-1
               xk = dble(k)
                  VRTDIJ = VRTDIJ + EGT(I,IVTOR1,IVTORP)                        
     &                 * (A1T2(I,IVTOR1,IVTORP)                                 
     &                     - 0.5d0 * (A4T2(I,IVTOR1,IVTORP)                     
     &                              + A5T2(I,IVTOR1,IVTORP)))                   
     &                 *(xK**2+(xK+2)**2)*xnn
               end do
            end do
         end do
      END IF                                                                    

      IF (PARA(LO,L).EQ.'C2KK  ') THEN                                          
         VRTDIJ = 0.d0                                                          
         DO IVTORP = 1,9                                                  
            DO IVTOR1 = 1,9                                               
               DO I = 1, 2 * N - 1                                         
               K=I-N-1
               xk = dble(k)
                  VRTDIJ = VRTDIJ + EGT(I,IVTOR1,IVTORP)                        
     &                 * (A1T2(I,IVTOR1,IVTORP)                                 
     &                     - 0.5d0 * (A4T2(I,IVTOR1,IVTORP)                     
     &                              + A5T2(I,IVTOR1,IVTORP)))                   
     &                 *(xK**4+(xK+2)**4)
               end do
            end do
         end do
      END IF                                                                    
*                                                                               
      IF (PARA(LO,L).EQ.'C11   '.or.PARA(LO,L).eq.'C6BCJ ') THEN              
         VRTDIK = 0.d0                                                          
         DO IVTORP = 1,9                                                  
            DO IVTOR1 = 1,9                                               
               DO I = 1, 2 * N - 1                                         
                  VRTDIK = VRTDIK + EGT(I,IVTOR1,IVTORP)                        
     &                 * (A1T2(I,IVTOR1,IVTORP)                                 
     &                     - 0.5d0 * (A8T2(I,IVTOR1,IVTORP)                     
     &                              + A9T2(I,IVTOR1,IVTORP)))                   
               end do
            end do
         end do
      END IF                                                                    

      IF (PARA(LO,L).EQ.'C6BCK ') THEN              
         VRTDIK = 0.d0                                                          
         DO IVTORP = 1,9                                                  
            DO IVTOR1 = 1,9                                               
               DO I = 1, 2 * N - 1                                         
                  K = I - N - 1  
                  xk = dble(k)
                  VRTDIK = VRTDIK + EGT(I,IVTOR1,IVTORP)                        
     &                 * (A1T2(I,IVTOR1,IVTORP)                                 
     &                     - 0.5d0 * (A8T2(I,IVTOR1,IVTORP)                     
     &                              + A9T2(I,IVTOR1,IVTORP)))                   
     &                 * (xK**2 + (xK+2)**2)
               end do
            end do
         end do
      END IF                                                                    

      IF (PARA(LO,L).EQ.'C6BCJK') THEN              
         VRTDIK = 0.d0                                                          
         DO IVTORP = 1,9                                                  
            DO IVTOR1 = 1,9                                               
               DO I = 1, 2 * N - 1                                         
                  K = I - N - 1     
                  xk = dble(k)
                  VRTDIK = VRTDIK + EGT(I,IVTOR1,IVTORP)                        
     &                 * (A1T2(I,IVTOR1,IVTORP)                                 
     &                     - 0.5d0 * (A8T2(I,IVTOR1,IVTORP)                     
     &                              + A9T2(I,IVTOR1,IVTORP)))                   
     &                 * (xK**2 + (xK+2)**2)*xnn
               end do
            end do
         end do
      END IF                                                                    

      IF (PARA(LO,L).EQ.'C6BCKK') THEN              
         VRTDIK = 0.d0                                                          
         DO IVTORP = 1,9                                                  
            DO IVTOR1 = 1,9                                               
               DO I = 1, 2 * N - 1                                         
                  K = I - N - 1                                                 
                  xk = dble(k)
                  VRTDIK = VRTDIK + EGT(I,IVTOR1,IVTORP)                        
     &                 * (A1T2(I,IVTOR1,IVTORP)                                 
     &                     - 0.5d0 * (A8T2(I,IVTOR1,IVTORP)                     
     &                              + A9T2(I,IVTOR1,IVTORP)))                   
     &                 * (xK**4 + (xK+2)**4)
               end do
            end do
         end do
      END IF                                                                    

      IF (PARA(LO,L).EQ.'CBCJ  ') THEN              
         VRTDIK = 0.d0                                                          
         DO IVTORP = 1,9                                                  
            DO IVTOR1 = 1,9                                               
               DO I = 1, 2 * N - 1                                         
                  VRTDIK = VRTDIK + EGT(I,IVTOR1,IVTORP)                        
     &                 * (A1T2(I,IVTOR1,IVTORP)                                 
     &                     - 0.5d0 * (A8T2(I,IVTOR1,IVTORP)                     
     &                              + A9T2(I,IVTOR1,IVTORP)))                   
               end do
            end do
         end do
      END IF                                                                    

      IF (PARA(LO,L).EQ.'CBCK  ') THEN              
         VRTDIK = 0.d0                                                          
         DO IVTORP = 1,9                                                  
            DO IVTOR1 = 1,9                                               
               DO I = 1, 2 * N - 1                                         

                  K=I-N-1
                  xk = dble(k)

                  VRTDIK = VRTDIK + EGT(I,IVTOR1,IVTORP)                        
     &                 * (A1T2(I,IVTOR1,IVTORP)                                 
     &                     - 0.5d0 * (A8T2(I,IVTOR1,IVTORP)                     
     &                              + A9T2(I,IVTOR1,IVTORP)))                   
     &                 * (xK**2 + (xK+2)**2)
               end do
            end do
         end do
      END IF                                                                    

      IF (PARA(LO,L).EQ.'C9BC  ') THEN              
         VRTDIK = 0.d0                                                          
         DO IVTORP = 1,9                                                  
            DO IVTOR1 = 1,9                                               
               DO I = 1, 2 * N - 1                                         
                  VRTDIK = VRTDIK + EGT(I,IVTOR1,IVTORP)                        
     &                 * (A1T2(I,IVTOR1,IVTORP)                                 
     &                     - 0.5d0 * (A12T2(I,IVTOR1,IVTORP)                     
     &                              + A13T2(I,IVTOR1,IVTORP)))                   
               end do
            end do
         end do
      END IF                                                                    
*                                                                               
      IF (PARA(LO,L).EQ.'C3    '.or.PARA(LO,l).eq.'G4BCJ ') THEN                                          
         VRTDIF = 0.d0                                                          
         DO IVTORP = 1,9                                                  
            DO IVTOR1 = 1,9                                               
               DO I = 1, 2 * N - 1                                         
                  VRTDIF = VRTDIF + EGT(I,IVTOR1,IVTORP)                        
     &                 * A3T2(I,IVTOR1,IVTORP)                                  
               end do
            end do
         end do
      END IF   

      IF (PARA(LO,L).EQ.'G4BCK ') THEN                                          
         VRTDIF = 0.d0                                                          
         DO IVTORP = 1,9                                                  
            DO IVTOR1 = 1,9                                               
               DO I = 1, 2 * N - 1                                         
                  K = I - N - 1                                                 
                  VRTDIF = VRTDIF + EGT(I,IVTOR1,IVTORP)                        
     &                 * A3T2K(I,IVTOR1,IVTORP)                                  
               end do
            end do
         end do
      END IF                                                                    

      IF (PARA(LO,L).EQ.'G4BCJK') THEN                                          
         VRTDIF = 0.d0                                                          
         DO IVTORP = 1,9                                                  
            DO IVTOR1 = 1,9                                               
               DO I = 1, 2 * N - 1                                         
                  K = I - N - 1                                                 
                  VRTDIF = VRTDIF + EGT(I,IVTOR1,IVTORP)                        
     &                 * A3T2K(I,IVTOR1,IVTORP)*xnn
               end do
            end do
         end do
      END IF                                                                    

      IF (PARA(LO,L).EQ.'G4BCKK') THEN                                          
         VRTDIF = 0.d0                                                          
         DO IVTORP = 1,9                                                  
            DO IVTOR1 = 1,9                                               
               DO I = 1, 2 * N - 1                                         
                  K = I - N - 1                                                 
                  VRTDIF = VRTDIF + EGT(I,IVTOR1,IVTORP)                        
     &                 * A3T2KK(I,IVTOR1,IVTORP)                                  
               end do
            end do
         end do
      END IF                                                                    

      IF (PARA(LO,L).EQ.'DBCK  '.or.para(lo,l).eq.'DBCJK ') THEN                                          
         VRTDIG = 0.d0                                                          
         DO IVTORP = 1,9                                                  
            DO IVTOR1 = 1,9                                               
               DO I = 1, 2 * N - 1                                         
                  K = I - N - 1                                                 
                  xk = dble(k)
                  VRTDIG = VRTDIG + EGT(I,IVTOR1,IVTORP)                        
     &                 * (A4T2(I,IVTOR1,IVTORP) - A5T2(I,IVTOR1,IVTORP))        
     &                 * (xK**2 + (xK+2)**2)
               end do
            end do
         end do
      END IF                                                                    

*                                                                               
      IF (PARA(LO,L).EQ.'DBC   '.or.para(lo,l).eq.'DBCJ  ') THEN                                          
         VRTDIG = 0.d0                                                          
         DO IVTORP = 1,9                                                  
            DO IVTOR1 = 1,9                                               
               DO I = 1, 2 * N - 1                                         
                  VRTDIG = VRTDIG + EGT(I,IVTOR1,IVTORP)                        
     &                 * (A4T2(I,IVTOR1,IVTORP) - A5T2(I,IVTOR1,IVTORP))        
               end do
            end do
         end do
      END IF                                                                    

      IF (PARA(LO,L).EQ.'DBC6  ') THEN                                          
         VRTDIG = 0.d0                                                          
         DO IVTORP = 1,9                                                  
            DO IVTOR1 = 1,9                                               
               DO I = 1, 2 * N - 1                                         
                  VRTDIG = VRTDIG + EGT(I,IVTOR1,IVTORP)                        
     &                 * (A9T2(I,IVTOR1,IVTORP) - A8T2(I,IVTOR1,IVTORP))        
               end do
            end do
         end do
      END IF                                                                    
*                                                                               
      IF (PARA(LO,L).EQ.'C4    '.OR.PARA(LO,L).EQ.'C4J   ')THEN                 
         VRTDIH = 0.d0                                                          
         DO IVTORP = 1,9                                                  
            DO IVTOR1 = 1,9                                               
               DO I = 1, 2 * N - 1                                         
                  VRTDIH = VRTDIH + EGT(I,IVTOR1,IVTORP)                        
     &                 * A7T2(I,IVTOR1,IVTORP)                                  
               end do
            end do
         end do
      END IF                                                                    
*                                                                               
      IF (PARA(LO,L).EQ.'C12   '.or.PARA(LO,L).eq.'AG3BCJ')THEN                                           
        VRTDIL=0.d0
         DO IVTORP = 1,9                                                  
            DO IVTOR1 = 1,9                                               
               DO I = 1, 2 * N - 1                                         
                  VRTDIL = VRTDIL + EGT(I,IVTOR1,IVTORP)                        
     &            *A10T2(I,IVTOR1,IVTORP)
               end do
            end do
         end do
      END IF                                                                    

      IF (PARA(LO,L).EQ.'AG3BCK')THEN                                           
        VRTDIL=0.d0
         DO IVTORP = 1,9                                                  
            DO IVTOR1 = 1,9                                               
               DO I = 1, 2 * N - 1                                         
                  K = I - N - 1                                                 
                  VRTDIL = VRTDIL + EGT(I,IVTOR1,IVTORP)                        
     &            *A10T2K(I,IVTOR1,IVTORP)
               end do
            end do
         end do
      END IF                                                                    

      IF (PARA(LO,L).EQ.'G3BCKK')THEN                                           
        VRTDIL=0.d0
         DO IVTORP = 1,9                                                  
            DO IVTOR1 = 1,9                                               
               DO I = 1, 2 * N - 1                                         
                  K = I - N - 1                                                 
                  VRTDIL = VRTDIL + EGT(I,IVTOR1,IVTORP)                        
     &            *A10T2KK(I,IVTOR1,IVTORP)
               end do
            end do
         end do
      END IF                                                                    

      IF (PARA(LO,L).EQ.'G3BCJK')THEN                                           
        VRTDIL=0.d0
         DO IVTORP = 1,9                                                  
            DO IVTOR1 = 1,9                                               
               DO I = 1, 2 * N - 1                                         
                  K = I - N - 1                                                 
                  VRTDIL = VRTDIL + EGT(I,IVTOR1,IVTORP)                        
     &            *A10T2K(I,IVTOR1,IVTORP)*xnn
               end do
            end do
         end do
      END IF                                                                    
*                                                                               
      IF (PARA(LO,L).EQ.'C4K   ')THEN                                           
         VRTDIM = 0.d0                                                          
         DO IVTORP = 1,9                                                  
            DO IVTOR1 = 1,9                                               
               DO I = 1, 2 * N - 1                                         
                  VRTDIM = VRTDIM + EGT(I,IVTOR1,IVTORP)                        
     &                 * A11T2(I,IVTOR1,IVTORP)                                 
               end do
            end do
         end do
      END IF                                                                    

      IF (PARA(LO,L).EQ.'C4KJ  ')THEN                                           
         VRTDIM = 0.d0                                                          
         DO IVTORP = 1,9                                                  
            DO IVTOR1 = 1,9                                               
               DO I = 1, 2 * N - 1                                         
                  VRTDIM = VRTDIM + EGT(I,IVTOR1,IVTORP)                        
     &                 * A11T2(I,IVTOR1,IVTORP)*xnn
               end do
            end do
         end do
      END IF                                                                    

      IF (PARA(LO,L).EQ.'C4KK  ')THEN                                           
         VRTDIM = 0.d0                                                          
         DO IVTORP = 1,9                                                  
            DO IVTOR1 = 1,9                                               
               DO I = 1, 2 * N - 1                                         
                  VRTDIM = VRTDIM + EGT(I,IVTOR1,IVTORP)                        
     &                 * A11T2K(I,IVTOR1,IVTORP)                                 
               end do
            end do
         end do
      END IF                                                                    
*                                                                               
      IF (PARA(LO,L).EQ.'B     ') THEN                                          
         VRTDI2 = 0.25d0 * VRTDIA                                               
      ELSE IF (PARA(LO,L).EQ.'C     ') THEN                                     
          VRTDI2 = -0.25d0 * VRTDIA                                             
      ELSE IF (PARA(LO,L).EQ.'ODELN ') THEN                                     
         VRTDI2 = - VRTDIA * xnn                                            
      ELSE IF (PARA(LO,L).EQ.'ODELK ') THEN                                     
         VRTDI2 = -0.5d0 * VRTDIB                                               
      ELSE IF (PARA(LO,L).EQ.'OHJ   ') THEN                                     
         VRTDI2 = VRTDIA * xnn**2
      ELSE IF (PARA(LO,L).EQ.'OHJK  ') THEN                                     
         VRTDI2 = 0.5d0 * VRTDIB * xnn
      ELSE IF (PARA(LO,L).EQ.'OHK   ') THEN                                     
         VRTDI2 = 0.5d0 * VRTDIC                                                
      ELSE IF (PARA(LO,L).EQ.'C1    ') THEN                                     
         VRTDI2 = 0.5d0 * VRTDID                                                
      ELSE IF (PARA(LO,L).EQ.'C1J   ') THEN                                     
         VRTDI2 = 0.5d0 * VRTDID * xnn
      ELSE IF (PARA(LO,L).EQ.'C2    ') THEN                                     
         VRTDI2 = 0.5d0 * VRTDIE                                                
      ELSE IF (PARA(LO,L).EQ.'AK4BC ') THEN                                     
         VRTDI2 = 0.5d0 * VRTDIE                                                
      ELSE IF (PARA(LO,L).EQ.'AK3BC ') THEN                                     
         VRTDI2 = 0.5d0 * VRTDIE                                                
      ELSE IF (PARA(LO,L).EQ.'AK7BC ') THEN                                     
         VRTDI2 = 0.5d0 * VRTDIE                                                
      ELSE IF (PARA(LO,L).EQ.'C2J   ') THEN                                     
         VRTDI2 = 0.5d0 * VRTDIE * xnn
      ELSE IF (PARA(LO,L).EQ.'C11   ') THEN
         VRTDI2 = 0.5d0 * VRTDIK
      ELSE IF (PARA(LO,L).EQ.'C6BCK ') THEN
         VRTDI2 = 0.5d0 * VRTDIK
      ELSE IF (PARA(LO,L).EQ.'C6BCJK') THEN
         VRTDI2 = 0.5d0 * VRTDIK
      ELSE IF (PARA(LO,L).EQ.'C6BCKK') THEN
         VRTDI2 = 0.5d0 * VRTDIK
      ELSE IF (PARA(LO,L).EQ.'C6BCJ ') THEN
         VRTDI2 = 0.5d0 * VRTDIK * xnn
      ELSE IF (PARA(LO,L).EQ.'CBCJ  ') THEN
         VRTDI2 = 0.5d0 * VRTDIK * xnn
      ELSE IF (PARA(LO,L).EQ.'CBCK  ') THEN
         VRTDI2 = 0.5d0 * VRTDIK
      ELSE IF (PARA(LO,L).EQ.'C9BC  ') THEN
         VRTDI2 = 0.5d0 * VRTDIK
      ELSE IF (PARA(LO,L).EQ.'C3    ') THEN                                     
         VRTDI2 = 0.5d0 * VRTDIF                                                
      ELSE IF (PARA(LO,L).EQ.'G4BCK ') THEN                                     
         VRTDI2 = 0.5d0 * VRTDIF                                                
      ELSE IF (PARA(LO,L).EQ.'G4BCJK') THEN                                     
         VRTDI2 = 0.5d0 * VRTDIF                                                
      ELSE IF (PARA(LO,L).EQ.'G4BCKK') THEN                                     
         VRTDI2 = 0.5d0 * VRTDIF                                                
      ELSE IF (PARA(LO,L).EQ.'G4BCJ ') THEN                                     
         VRTDI2 = 0.5d0 * VRTDIF * xnn
      ELSE IF (PARA(LO,L).EQ.'DBC   ') THEN                                     
         VRTDI2 = -0.25d0 * VRTDIG                                              
      ELSE IF (PARA(LO,L).EQ.'DBC6  ') THEN                                     
         VRTDI2 = -0.25d0 * VRTDIG                                              
      ELSE IF (PARA(LO,L).EQ.'DBCK  ') THEN                                     
         VRTDI2 = -0.25d0 * VRTDIG                                              
      ELSE IF (PARA(LO,L).EQ.'DBCJK ') THEN                                     
         VRTDI2 = -0.25d0 * VRTDIG *xnn
      ELSE IF (PARA(LO,L).EQ.'DBCJ  ') THEN                                     
         VRTDI2 = -0.25d0 * VRTDIG * xnn
      ELSE IF (PARA(LO,L).EQ.'C4    ') THEN                                     
         VRTDI2 = 0.5d0 * VRTDIH                                                
      ELSE IF (PARA(LO,L).EQ.'C4J   ') THEN                                     
         VRTDI2 = 0.5d0 * VRTDIH *xnn
      ELSE IF (PARA(LO,L).EQ.'C1K   ') THEN                                     
         VRTDI2 = 0.5d0 * VRTDII                                                
      ELSE IF (PARA(LO,L).EQ.'C1JK  ') THEN                                     
         VRTDI2 = 0.5d0 * VRTDII                                                
      ELSE IF (PARA(LO,L).EQ.'C1KK  ') THEN                                     
         VRTDI2 = 0.5d0 * VRTDII                                                
      ELSE IF (PARA(LO,L).EQ.'C2K   ') THEN                                     
         VRTDI2 = 0.5d0 * VRTDIJ                                                
      ELSE IF (PARA(LO,L).EQ.'C2KK  ') THEN                                     
         VRTDI2 = 0.5d0 * VRTDIJ                                                
      ELSE IF (PARA(LO,L).EQ.'C2JK  ') THEN                                     
         VRTDI2 = 0.5d0 * VRTDIJ                                                
      ELSE IF (PARA(LO,L).EQ.'C12   ') THEN                                     
         VRTDI2 = 0.5d0 * VRTDIL                                                
      ELSE IF (PARA(LO,L).EQ.'AG3BCK') THEN                                     
         VRTDI2 = 0.5d0 * VRTDIL                                                
      ELSE IF (PARA(LO,L).EQ.'G3BCJK') THEN                                     
         VRTDI2 = 0.5d0 * VRTDIL                                                
      ELSE IF (PARA(LO,L).EQ.'G3BCKK') THEN                                     
         VRTDI2 = 0.5d0 * VRTDIL                                                
      ELSE IF (PARA(LO,L).EQ.'AG3BCJ') THEN                                     
         VRTDI2 = 0.5d0 * VRTDIL * xnn
      ELSE IF (PARA(LO,L).EQ.'C4K   ') THEN                                     
         VRTDI2 = 0.5d0 * VRTDIM                                                
      ELSE IF (PARA(LO,L).EQ.'C4JK  ') THEN                                     
         VRTDI2 = 0.5d0 * VRTDIM                            
      ELSE IF (PARA(LO,L).EQ.'C4KK  ') THEN                                     
         VRTDI2 = 0.5d0 * VRTDIM                                                              
      else if (PARA(LO,L).eq.'OLKKK ') then
         vrtdi2 = vrtdin
      else if (PARA(LO,L).eq.'OLJKK ') then
         vrtdi2 = vrtdic * xnn
      else if (PARA(LO,L).eq.'OLJJK ') then
         vrtdi2 = vrtdib * (xnn)**2
      else if (PARA(LO,L).eq.'OLJJJ ') then
         vrtdi2 = vrtdia * (xnn)**3
      END IF                                                                    
      END                                                                       
C***********************************************************                    
      FUNCTION VRTDI4 (L, LO)                                                   
      implicit none
      save
C**********************************************************                     
C     DERIVATIVES INVOLVING THE K/K+OR-1 MATRIX ELEMENTS OF THE                 
C     MATRIX ELEMENT(DAB,ODAB,DABJ,DABK)                                        
C         

      include 'common_dim.f'
      include 'common_eigen.f'
      include 'common_chara.f'
      include 'common_quant.f'
      include 'common_vrtd.f'

      double precision vrtdi4
      integer lo, ivtorp, ivtor1, l, i, k
*                                                                               
      double precision EGT(KDMX,9,9),SQT(KDMX)                                            
*             
      xn = dble(n)
      xnn = dble(n*(n+1))
                                                                  
      DO  I = 1, 2 * N                                                       
         SQT(I) = SQRT (dble(N*(N+1) - (I-N-1)*(I-N)))                         
      end do
      DO IVTORP = 1,9                                                     
         DO IVTOR1 = 1,9                                                  
            DO I = 1, 2 * N                                                
               EGT(I,IVTOR1,IVTORP) = SQT(I) * EGVC3(I,IVTORP)                  
     &                 * EGVC3(I+1,IVTOR1)                                      
            end do
         end do
      end do
*     
      VRTDI4 = 0.d0                                                             
*                                                                               
      IF (PARA (LO, L) .EQ. 'DAB   ') THEN                                      
         DO IVTORP = 1,9                                                  
            DO IVTOR1 = 1,9                                               
               DO I = 1, 2 * N                                             
                  VRTDI4 = VRTDI4 + 2.d0 * EGT(I,IVTOR1,IVTORP)                 
     &                 * (I - N - 0.5d0) * A1T4(I,IVTOR1,IVTORP)                
               end do
            end do
         end do
      else  IF (PARA (LO, L) .EQ. 'GB    ') THEN                                      
         DO IVTORP = 1,9                                                  
            DO IVTOR1 = 1,9                                               
               DO I = 1, 2 * N                                             
                  VRTDI4 = VRTDI4 + 2.d0 * EGT(I,IVTOR1,IVTORP)                 
     &                 * A1T4B(I,IVTOR1,IVTORP) * 0.5d0
               end do
            end do
         end do

      ELSE IF (PARA (LO, L) .EQ. 'ODAB  ') THEN                                 
         DO IVTORP = 1,9                                                  
            DO IVTOR1 = 1,9                                               
               DO I = 1, 2 * N                                             
                  VRTDI4 = VRTDI4 + 2. * EGT(I,IVTOR1,IVTORP)                   
     &                 * (I - N - 0.5d0)                                        
     &                 * (A1T4(I,IVTOR1,IVTORP)                                 
     &                 - 0.5d0 * (A2T4(I,IVTOR1,IVTORP)                      
     &                 + A3T4(I,IVTOR1,IVTORP)))                  
               end do
            end do
         end do

      ELSE IF (PARA (LO, L) .EQ. 'AK7AB ') THEN                                 
         DO IVTORP = 1,9                                                  
            DO IVTOR1 = 1,9                                               
               DO I = 1, 2 * N                                             
                  VRTDI4 = VRTDI4 + 2. * EGT(I,IVTOR1,IVTORP)                   
     &                 * (I - N - 0.5d0)                                        
     &                 * (2.0d0*AK7T1(I,IVTOR1,IVTORP)                                 
     &                 - 0.5d0 * (AK7T2(I,IVTOR1,IVTORP) +
     &                 AK7T3(I, IVTOR1, IVTORP)))
               end do
            end do
         end do

      ELSE IF (PARA (LO, L) .EQ. 'ODABJ ') THEN                                 
         DO  IVTORP = 1,9                                                  
            DO  IVTOR1 = 1,9                                               
               DO  I = 1, 2 * N                                             
                  VRTDI4 = VRTDI4 + 2. * EGT(I,IVTOR1,IVTORP)                   
     &                 * 0.5*((I - N - 1) + (I - N))                                                              
     &                 * (A1T4(I,IVTOR1,IVTORP)                                 
     &                 - 0.5d0 * (A2T4(I,IVTOR1,IVTORP)                      
     &                 + A3T4(I,IVTOR1,IVTORP)))
     &                 * xnn
               end do
            end do
         end do

      ELSE IF (PARA (LO, L) .EQ. 'ODABK ') THEN                                 
         DO  IVTORP = 1,9                                                  
            DO  IVTOR1 = 1,9                                               
               DO  I = 1, 2 * N                                             
                  VRTDI4 = VRTDI4 + 2. * EGT(I,IVTOR1,IVTORP)                   
     &                 * 0.5*((I - N - 1)**3 + (I - N)**3)    
     &                 * (A1T4(I,IVTOR1,IVTORP)                                 
     &                 - 0.5d0 * (A2T4(I,IVTOR1,IVTORP)                      
     &                 + A3T4(I,IVTOR1,IVTORP)))                  
               end do
            end do
         end do

      ELSE IF (PARA (LO, L) .EQ. 'ODABKK') THEN                                 
         DO  IVTORP = 1,9                                                  
            DO  IVTOR1 = 1,9                                               
               DO  I = 1, 2 * N                                             
                  VRTDI4 = VRTDI4 + 2. * EGT(I,IVTOR1,IVTORP)                   
     &                 * 0.5*((I - N - 1)**5 + (I - N)**5)    
     &                 * (A1T4(I,IVTOR1,IVTORP)                                 
     &                 - 0.5d0 * (A2T4(I,IVTOR1,IVTORP)                      
     &                 + A3T4(I,IVTOR1,IVTORP)))                  
               end do
            end do
         end do

      ELSE IF (PARA (LO, L) .EQ. 'ODAB6 ') THEN                                 
         DO IVTORP = 1,9                                                  
            DO IVTOR1 = 1,9                                               
               DO I = 1, 2 * N                                             
                  VRTDI4 = VRTDI4 + 2.d0 * EGT(I,IVTOR1,IVTORP)                 
     &                 * (I - N - 0.5d0)                                        
     &                 * (A1T4(I,IVTOR1,IVTORP)                                 
     &                    - 0.5d0 * (A20T4(I,IVTOR1,IVTORP)                     
     &                               + A21T4(I,IVTOR1,IVTORP)))                 
               end do
            end do
         end do
      ELSE IF (PARA (LO, L) .EQ. 'ODAB9 ') THEN                                 
         DO IVTORP = 1,9                                                  
            DO IVTOR1 = 1,9                                               
               DO I = 1, 2 * N                                             
                  VRTDI4 = VRTDI4 + 2.d0 * EGT(I,IVTOR1,IVTORP)                 
     &                 * (I - N - 0.5d0)                                        
     &                 * (A1T4(I,IVTOR1,IVTORP)                                 
     &                    - 0.5d0 * (A22T4(I,IVTOR1,IVTORP)                     
     &                               + A23T4(I,IVTOR1,IVTORP)))                 
               end do
            end do
         end do
      ELSE IF (PARA (LO, L) .EQ. 'CABJ  ') THEN                                 
         DO IVTORP = 1,9                                                  
            DO IVTOR1 = 1,9                                               
               DO I = 1, 2 * N                                             
                  VRTDI4 = VRTDI4 + 2.d0 * EGT(I,IVTOR1,IVTORP)                 
     &                 * (I - N - 0.5d0)                                        
     &                 * (A1T4(I,IVTOR1,IVTORP)                                 
     &                    - 0.5d0 * (A20T4(I,IVTOR1,IVTORP)                     
     &                               + A21T4(I,IVTOR1,IVTORP)))                 
     &                 * xnn
               end do
            end do
         end do
      ELSE IF (PARA (LO, L) .EQ. 'CABK  ') THEN                                 
         DO IVTORP = 1,9                                                  
            DO IVTOR1 = 1,9                                               
               DO I = 1, 2 * N           

                  K = I - N - 1
                            
                  xk = dble(k)
                  VRTDI4 = VRTDI4 + 2.d0 * EGT(I,IVTOR1,IVTORP)                 
     &                 * 0.5d0 * ( xK**3 + (xK+1)**3 ) 
     &                 * (A1T4(I,IVTOR1,IVTORP)                                 
     &                    - 0.5d0 * (A20T4(I,IVTOR1,IVTORP)                     
     &                               + A21T4(I,IVTOR1,IVTORP)))                 
               end do
            end do
         end do
      ELSE IF (PARA (LO, L) .EQ. 'DABJ  ') THEN                                 
         DO IVTORP = 1,9                                                  
            DO IVTOR1 = 1,9                                               
               DO I = 1, 2 * N                                             
                  VRTDI4 = VRTDI4 + 2.d0 * EGT(I,IVTOR1,IVTORP)                 
     &                 * (I - N - 0.5d0)                                        
     &                 * A1T4(I,IVTOR1,IVTORP) * xnn
               end do
            end do
         end do
      ELSE IF (PARA (LO, L) .EQ. 'DABK  ') THEN                                 
         DO IVTORP = 1,9                                                  
            DO IVTOR1 = 1,9                                               
               DO I = 1, 2 * N                         
                  xk = dble(I - N - 1)
C                 VRTDI4 = VRTDI4 + 2.d0 * EGT(I,IVTOR1,IVTORP)                 
C    &                 * (I - N - 0.5d0)                                        
C    &                 * A1T4(I,IVTOR1,IVTORP)                                  
C    &                 * xk * xk
C*******DABK DERIVATIVE MADE HERMITIAN BY JTH 18 APR 92                         
                  VRTDI4 = VRTDI4 + 2.d0 * EGT(I,IVTOR1,IVTORP)                 
     &                 * A1T4(I,IVTOR1,IVTORP)                                  
C    &                 * (I - N - 0.5d0)                                        
C    &                 * xk * xk
     &                 * 0.5*(xk**3 + dble(I - N)**3)                      
               end do
            end do
         end do

      ELSE IF (PARA (LO, L) .EQ. 'HABK  ') THEN                                 
         DO IVTORP = 1,9                                                  
            DO IVTOR1 = 1,9                                               
               DO I = 1, 2 * N                                             

                  VRTDI4 = VRTDI4 + 2.d0 * EGT(I,IVTOR1,IVTORP)                 
     &                 * A1T4(I,IVTOR1,IVTORP)                                  
     &                 * (dble(I - N - 1)**5 + dble(I - N)**5)*0.5d0
               end do
            end do
         end do

      ELSE IF (PARA (LO, L) .EQ. 'HABJ  ') THEN                                 
         DO IVTORP = 1,9                                                  
            DO IVTOR1 = 1,9                                               
               DO I = 1, 2 * N                                             
                  VRTDI4 = VRTDI4 + 2.d0 * EGT(I,IVTOR1,IVTORP)                 
     &                 * A1T4(I,IVTOR1,IVTORP)                                  
     &                 * ((I - N - 1)**1 + (I - N)**1)*(xnn)**2                      
     &                 * 0.5d0
               end do
            end do
         end do

      ELSE IF (PARA (LO, L) .EQ. 'HABJK ') THEN                                 
         DO IVTORP = 1,9                                                  
            DO IVTOR1 = 1,9                                               
               DO I = 1, 2 * N                                             
                  VRTDI4 = VRTDI4 + 2.d0 * EGT(I,IVTOR1,IVTORP)                 
     &                 * A1T4(I,IVTOR1,IVTORP)                                  
     &                 * (dble(I - N - 1)**3 + dble(I - N)**3)*(xnn)                      
     &                 * 0.5d0
               end do
            end do
         end do

      ELSE IF (PARA (LO, L) .EQ. 'DAC   ') THEN                                 
         DO IVTORP = 1,9                                                  
            DO IVTOR1 = 1,9                                               
               DO I = 1, 2 * N                                             
C ********I BELIEVE THERE IS AN EXTRA FACTOR OF 2 IN THE NEXT LINE, AND         
C ********I HAVE REMOVED IT                                                     
C ********VRTDIS4=VRTDIS4+2.*EG*SQ*(K+0.5)*(-1.)*(A4-A3)                        
                  VRTDI4 = VRTDI4 + EGT(I,IVTOR1,IVTORP)                        
     &                 * (I - N - 0.5d0)                                        
     &                 * (A3T4(I,IVTOR1,IVTORP) - A4T4(I,IVTOR1,IVTORP))        
               end do
            end do
         end do

      ELSE IF (PARA (LO, L) .EQ. 'DACJJ ') THEN                                 
         DO IVTORP = 1,9                                                  
            DO IVTOR1 = 1,9                                               
               DO I = 1, 2 * N                                             
                  VRTDI4 = VRTDI4 + EGT(I,IVTOR1,IVTORP)* xnn**2
     &                 * (I - N - 0.5d0)                                        
     &                 * (A3T4(I,IVTOR1,IVTORP) - A4T4(I,IVTOR1,IVTORP))        
               end do
            end do
         end do

      ELSE IF (PARA (LO, L) .EQ. 'DAC6  ') THEN                                 
         DO IVTORP = 1,9                                                  
            DO IVTOR1 = 1,9                                               
               DO I = 1, 2 * N                                             
                  VRTDI4 = VRTDI4 + EGT(I,IVTOR1,IVTORP)                        
     &                 * (I - N - 0.5d0)                                        
     &                 * (A20T4(I,IVTOR1,IVTORP) 
     &                 - A21T4(I,IVTOR1,IVTORP))        
               end do
            end do
         end do

      ELSE IF (PARA (LO, L) .EQ. 'ODAC  ') THEN                                 
         DO IVTORP = 1,9                                                  
            DO IVTOR1 = 1,9                                               
               DO I = 1, 2 * N                                             
                  K = I - N - 1
                  VRTDI4 = VRTDI4 + EGT(I,IVTOR1,IVTORP)                        
     &                 * (K + (K+1))
     &                 * (A3T4(I,IVTOR1,IVTORP)/4.d0
     &                 - A4T4(I,IVTOR1,IVTORP)/4.d0
     &                 + A12T5(I,IVTOR1,IVTORP)/8.d0
     &                 - A13T5(I,IVTOR1,IVTORP)/8.d0)        
               end do
            end do
         end do

      ELSE IF (PARA (LO, L) .EQ. 'DACK  ') THEN                                 
         DO IVTORP = 1,9                                                  
            DO IVTOR1 = 1,9                                               
               DO I = 1, 2 * N                                             
                  VRTDI4 = VRTDI4 + EGT(I,IVTOR1,IVTORP)                        
     &                 * 0.5d0 * ( dble(I - N)**3 + dble(I - N - 1)**3)                                       
     &                 * (A3T4(I,IVTOR1,IVTORP) - A4T4(I,IVTOR1,IVTORP))
               end do
            end do
         end do

      ELSE IF (PARA (LO, L) .EQ. 'DACJK ') THEN                                 
         DO IVTORP = 1,9                                                  
            DO IVTOR1 = 1,9                                               
               DO I = 1, 2 * N                                             
                  VRTDI4 = VRTDI4 + EGT(I,IVTOR1,IVTORP) *xnn
     &                 * 0.5d0 * ( dble(I - N)**3 + dble(I - N - 1)**3)                                       
     &                 * (A3T4(I,IVTOR1,IVTORP) - A4T4(I,IVTOR1,IVTORP))
               end do
            end do
         end do

      ELSE IF (PARA (LO, L) .EQ. 'DACKK ') THEN                                 
         DO IVTORP = 1,9                                                  
            DO IVTOR1 = 1,9                                               
               DO I = 1, 2 * N                                             
                  VRTDI4 = VRTDI4 + EGT(I,IVTOR1,IVTORP)                        
     &                 * 0.5d0 * ( dble(I - N)**5 + dble(I - N - 1)**5)                                       
     &                 * (A3T4(I,IVTOR1,IVTORP) - A4T4(I,IVTOR1,IVTORP))
               end do
            end do
         end do

      ELSE IF (PARA (LO, L) .EQ. 'DACJ  ') THEN                                 
         DO IVTORP = 1,9                                                  
            DO IVTOR1 = 1,9                                               
               DO I = 1, 2 * N                                             
                  VRTDI4 = VRTDI4 + EGT(I,IVTOR1,IVTORP)                        
     &                 * (I - N - 0.5d0)                                        
     &                 * (A3T4(I,IVTOR1,IVTORP) - A4T4(I,IVTOR1,IVTORP))
     &                 * xnn
               end do
            end do
         end do

      ELSE IF (PARA (LO, L) .EQ. 'ODELTA') THEN                                 
         DO IVTORP = 1,9                                                  
            DO IVTOR1 = 1,9                                               
               DO I = 1, 2 * N                                             
                  VRTDI4 = VRTDI4 + EGT(I,IVTOR1,IVTORP)                        
     &                 * A5T4(I,IVTOR1,IVTORP)                                  
               end do
            end do
         end do

      ELSE IF (PARA (LO, L) .EQ. 'DAGJ  ') THEN                                 
         DO IVTORP = 1,9                                                  
            DO IVTOR1 = 1,9                                               
               DO I = 1, 2 * N                                             
                  VRTDI4 = VRTDI4 + EGT(I,IVTOR1,IVTORP)                        
     &                 * A5T4(I,IVTOR1,IVTORP) * xnn
               end do
            end do
         end do

      ELSE IF (PARA (LO, L) .EQ. 'DAGK  ') THEN                                 
         DO IVTORP = 1,9                                                  
            DO IVTOR1 = 1,9                                               
               DO I = 1, 2 * N                                             
                  VRTDI4 = VRTDI4 + EGT(I,IVTOR1,IVTORP)                        
     &                 * A8T4(I,IVTOR1,IVTORP)                                  
               end do
            end do
         end do

      ELSE IF (PARA (LO, L) .EQ. 'DELTA ') THEN                                 
         DO IVTORP = 1,9                                                  
            DO IVTOR1 = 1,9                                               
               DO I = 1, 2 * N                                            
                  VRTDI4 = VRTDI4 + EGT(I,IVTOR1,IVTORP)                        
     &                 * A6T4(I,IVTOR1,IVTORP)                                  
               end do
            end do
         end do

      ELSE IF (PARA (LO, L) .EQ. 'DELTAJ') THEN                                 
         DO IVTORP = 1,9                                                  
            DO IVTOR1 = 1,9                                               
               DO I = 1, 2 * N                                            
                  VRTDI4 = VRTDI4 + EGT(I,IVTOR1,IVTORP)                        
     &                 * A6T4(I,IVTOR1,IVTORP) * xnn
               end do
            end do
         end do

      ELSE IF (PARA (LO, L) .EQ. 'DELTAK') THEN                                 
         DO IVTORP = 1,9                                                  
            DO IVTOR1 = 1,9                                               
               DO I = 1, 2 * N                                            
                  VRTDI4 = VRTDI4 + EGT(I,IVTOR1,IVTORP)                        
     &                 * A7T4(I,IVTOR1,IVTORP)
               end do
            end do
         end do

      ELSE IF (PARA (LO, L) .EQ. 'DELTKK') THEN                                 
         DO IVTORP = 1,9                                                  
            DO IVTOR1 = 1,9                                               
               DO I = 1, 2 * N                                            
                  VRTDI4 = VRTDI4 + EGT(I,IVTOR1,IVTORP)                        
     &                 * A9T4(I,IVTOR1,IVTORP)
               end do
            end do
         end do

      ELSE IF (PARA (LO, L) .EQ. 'DELTJK') THEN                                 
         DO IVTORP = 1,9                                                  
            DO IVTOR1 = 1,9                                               
               DO I = 1, 2 * N                                            
                  VRTDI4 = VRTDI4 + EGT(I,IVTOR1,IVTORP)                        
     &                 * A7T4(I,IVTOR1,IVTORP) * xnn
               end do
            end do
         end do
*                                                                               
      ELSE IF (PARA (LO, L) .EQ. 'ODELTB') THEN                                 
         DO IVTORP = 1,9                                                  
            DO IVTOR1 = 1,9                                               
               DO I = 1, 2 * N                                            
                  VRTDI4 = VRTDI4 + EGT(I,IVTOR1,IVTORP)                        
     &                 * A5T4B(I,IVTOR1,IVTORP)                                 
               end do
            end do
         end do
*                                                                               
      ELSE IF (PARA (LO, L) .EQ. 'DELTAB') THEN                                 
         DO IVTORP = 1,9                                                  
            DO IVTOR1 = 1,9                                               
               DO I = 1, 2 * N                                            
                  VRTDI4 = VRTDI4 + EGT(I,IVTOR1,IVTORP)                        
     &                 * A6T4B(I,IVTOR1,IVTORP)                                  
               end do
            end do
         end do

      ELSE IF (PARA (LO, L) .EQ. 'DG4J  ') THEN                                 
         DO IVTORP = 1,9                                                  
            DO IVTOR1 = 1,9                                               
               DO I = 1, 2 * N                                            
                  VRTDI4 = VRTDI4 + EGT(I,IVTOR1,IVTORP)                        
     &                 * A6T4B(I,IVTOR1,IVTORP) *xnn
               end do
            end do
         end do

      ELSE IF (PARA (LO, L) .EQ. 'DG4K  ') THEN                                 
         DO IVTORP = 1,9                                                  
            DO IVTOR1 = 1,9                                               
               DO I = 1, 2 * N                                            
                  VRTDI4 = VRTDI4 + EGT(I,IVTOR1,IVTORP)                        
     &                 * A6T4C(I,IVTOR1,IVTORP)                                  
               end do
            end do
         end do

      END IF                                                                    
*                                                                               
      END                                                                       
C********************************************************************           
      FUNCTION VRTDI5 (L, LO)                                                   
      implicit none
      save
C********************************************************************           
C     DERIVATIVES INVOLVING THE K/K MATRIX ELEMENTS OF THE                      
C     HAMILTONIAN BUT WITH THE SUM OVER THE A COEFICIENTS                       
C     (FV,OFV,GV,AK1-AK7,ALV,OLV)                                               
C         

      include 'common_dim.f'
      include 'common_eigen.f'
      include 'common_chara.f'
      include 'common_quant.f'
      include 'common_vrtd.f'

      double precision vrtdi5, vrtdia, vrtdib, vrtdic, vrtdid
      double precision vrtdie, vrtdif, vrtdiz, vrtdiq, vrtdig
      double precision vrtdiy, vrtdix, vrtdih, vrtdij, vrtdiv
      double precision vrtdik, vrtdil, vrtdim, vrtdiu, vrtdiw
      double precision vrtdin, vrtdio, vrtdip
      integer lo, ivtorp, ivtor1, l, i, k
*                                                                               
      double precision EGT(KDMX,9,9)                   

      xn = dble(n)
      xnn = dble(n*(n+1))

C     JMF
      vrtdi5 = 0.0d0
*                                                                               
      DO IVTORP = 1,9                                                     
         DO IVTOR1 = 1,9                                                  
            DO I = 1, 2 * N + 1                                            
               EGT(I,IVTOR1,IVTORP) = EGVC3(I,IVTORP) * EGVC3(I,IVTOR1)         
            end do
         end do
      end do
*                                                                               
      IF (PARA(LO,L).EQ.'FV    ' .OR. PARA(LO,L).EQ.'OFV   ') THEN              
         VRTDIA = 0.d0                                                          
         DO IVTORP = 1,9                                                  
            DO IVTOR1 = 1,9                                               
               DO I = 1, 2 * N + 1                                         
                  VRTDIA = VRTDIA + EGT(I,IVTOR1,IVTORP)                        
     &                 * (A1T1(I,IVTOR1,IVTORP)                                 
     &                    - 0.5d0*A7T5(I,IVTOR1,IVTORP)                         
     &                    - 0.5d0*A1T5(I,IVTOR1,IVTORP))                        
               end do
            end do
         end do
      END IF                                                                    

*                                                                               
      IF (PARA(LO,L).EQ.'GV    ' .OR. PARA(LO,L).EQ.'GVJ   '
     &.or.PARA(LO,L).EQ.'GJ6   ') THEN              
         VRTDIB = 0.d0                                                          
         DO IVTORP = 1,9                                                  
            DO IVTOR1 = 1,9                                               
               DO I = 1, 2 * N + 1                                         
                  VRTDIB = VRTDIB + EGT(I,IVTOR1,IVTORP)                        
     &                 * A2T5(I,IVTOR1,IVTORP)                                  
               end do
            end do
         end do
      END IF                                                                    
*                                                                               
      IF (PARA(LO,L).EQ.'ALV   '.OR.PARA(LO,L).EQ.'OLV   '
     &.or.PARA(LO,L).eq.'AGJ6  '.or.PARA(LO,L).eq.'AGJ4K2'
     &.or.PARA(LO,L).eq.'AGJ2K4'.or.PARA(LO,L).eq.'AGK6  ')THEN                 
         VRTDIC = 0.d0                                    
         VRTDID = 0.d0
         VRTDIE = 0.d0
         VRTDIF = 0.d0
         DO IVTORP = 1,9                                                  
            DO IVTOR1 = 1,9                                               
               DO I = 1, 2 * N + 1                                         
                  K = I - N - 1                                                 
                  xk = dble(k)
                  VRTDIC = VRTDIC + EGT(I,IVTOR1,IVTORP) * xk
     &                 * A3T5(I,IVTOR1,IVTORP)                                  

                  VRTDID = VRTDID + EGT(I,IVTOR1,IVTORP) * xk**3
     &                 * A3T5(I,IVTOR1,IVTORP)                                  

                  VRTDIE = VRTDIE + EGT(I,IVTOR1,IVTORP) * xk**5
     &                 * A3T5(I,IVTOR1,IVTORP)                                  

                  VRTDIF = VRTDIF + EGT(I,IVTOR1,IVTORP) * xk**7
     &                 * A3T5(I,IVTOR1,IVTORP)                                  
               end do
            end do
         end do
      END IF                                                                    
*                                                                               
      IF (PARA(LO,L).EQ.'AK1   ' .OR. PARA(LO,L).EQ.'AK1J  ' )THEN              
C    &     PARA(LO,L).EQ.'OLV   ') THEN                                         
         VRTDID = 0.d0                                                          
         DO IVTORP = 1,9                                                  
            DO IVTOR1 = 1,9                                               
               DO I = 1, 2 * N + 1                                         
                  K = I - N - 1                                                 
                  xk = dble(k)**3
                  VRTDID = VRTDID + EGT(I,IVTOR1,IVTORP) * xk
     &                 * A3T5(I,IVTOR1,IVTORP)                                  
               end do
            end do
         end do
      END IF                                                                    
*                                                                               
      IF (PARA(LO,L).EQ.'AK2   ' .OR. PARA(LO,L).EQ.'AK2J  ' .OR.
     &    PARA(LO,L).EQ.'GJ4K2 ') THEN              
         VRTDIE = 0.d0                                                          
         DO IVTORP = 1,9                                                  
            DO IVTOR1 = 1,9                                               
               DO I = 1, 2 * N + 1                                         
                  K = I - N - 1                            
                  xk = dble(k)**2
                  VRTDIE = VRTDIE + EGT(I,IVTOR1,IVTORP) * xk
     &                 * A2T5(I,IVTOR1,IVTORP)                                  
               end do
            end do
         end do
      END IF                                                                    
*                                                                               
      IF (PARA(LO,L).EQ.'AK3   ' .OR. PARA(LO,L).EQ.'AK3J  ' .or.
     &     PARA(LO, L).eq.'GAJ4  ' .or. PARA(LO, L).eq. 'GAJ2K2' .or.
     &     PARA(LO, L).eq.'GAK4  ')THEN              
         VRTDIF = 0.d0                                                          
         VRTDIG = 0.d0                                                          
         VRTDIH = 0.d0                                                          
         DO  IVTORP = 1,9                                                  
            DO IVTOR1 = 1,9                                               
               DO  I = 1, 2 * N + 1                                         
                  K = I - N - 1                              
                  xk = dble(k)
                  VRTDIF = VRTDIF + EGT(I,IVTOR1,IVTORP) * xk                    
     &                 * A4T5(I,IVTOR1,IVTORP)                                  
                  VRTDIG = VRTDIG + EGT(I,IVTOR1,IVTORP) * xk**3
     &                 * A4T5(I,IVTOR1,IVTORP)                                  
                  VRTDIH = VRTDIH + EGT(I,IVTOR1,IVTORP) * xk**5                 
     &                 * A4T5(I,IVTOR1,IVTORP)                                  
               end do
            end do
         end do
      END IF                                                                    
*                                                                               
      IF (PARA(LO,L).EQ.'AK3K  ' ) THEN                                         
         VRTDIZ = 0.d0                                                          
         DO IVTORP = 1,9                                                  
            DO IVTOR1 = 1,9                                               
               DO I = 1, 2 * N + 1                                         
                  K = I - N - 1                                                 
                  xk = dble(k)**3
                  VRTDIZ = VRTDIZ + EGT(I,IVTOR1,IVTORP) * xk
     &                 * A4T5(I,IVTOR1,IVTORP)                                  
               end do
            end do
         end do
      END IF                                                                    
*                                                                               
      IF (PARA(LO,L).EQ.'AK3B  '.or.PARA(LO,L).eq.'AK3BK '.or.
     &     para(lo, l).eq.'AK3BJ ') THEN                                          
         VRTDIQ = 0.d0                                                          
         VRTDIW = 0.d0                                                          
         DO IVTORP = 1,9                                                  
            DO IVTOR1 = 1,9                                               
               DO I = 1, 2 * N + 1                                         
                  K = I - N - 1                                                 
                  xk = dble(k)
                  VRTDIQ = VRTDIQ + EGT(I,IVTOR1,IVTORP) * xk
     &                 * A4T5B(I,IVTOR1,IVTORP)                                 

                  VRTDIW = VRTDIW + EGT(I,IVTOR1,IVTORP) * xk**3
     &                 * A4T5B(I,IVTOR1,IVTORP)                                 
               end do
            end do
         end do
      END IF                                                                    

      IF (PARA(LO,L).EQ.'AK3BK6') then
         VRTDIQ = 0.d0                                                          
         DO IVTORP = 1,9                                                  
            DO IVTOR1 = 1,9                                               
               DO I = 1, 2 * N + 1                                         
                  K = I - N - 1                                 
                  xk = dble(k)
                  VRTDIQ = VRTDIQ + EGT(I,IVTOR1,IVTORP) * xk**7                    
     &                 * A4T5B(I,IVTOR1,IVTORP)
               end do
            end do
         end do
      END IF                                                                    

      IF (PARA(LO,L).EQ.'AK3BK4') then
         VRTDIQ = 0.d0                                                          
         DO IVTORP = 1,9                                                  
            DO IVTOR1 = 1,9                                               
               DO I = 1, 2 * N + 1                                         
                  K = I - N - 1                
                  xk = dble(k)
                  VRTDIQ = VRTDIQ + EGT(I,IVTOR1,IVTORP) * xk**5                    
     &                 * A4T5B(I,IVTOR1,IVTORP)
               end do
            end do
         end do
      END IF                                                                    

      IF (PARA(LO,L).EQ.'AK3C  ') THEN                                          
         VRTDIQ = 0.d0                                                          
         DO IVTORP = 1,9                                                  
            DO IVTOR1 = 1,9                                               
               DO I = 1, 2 * N + 1                                         
                  K = I - N - 1                                                 
                  VRTDIQ = VRTDIQ + EGT(I,IVTOR1,IVTORP) * K                    
     &                 * A3T5C(I,IVTOR1,IVTORP)                                 
               end do
            end do
         end do
      END IF                                                                    

      IF (PARA(LO,L).EQ.'AK3D  ') THEN                                          
         VRTDIQ = 0.d0                                                          
         DO IVTORP = 1,9                                                  
            DO IVTOR1 = 1,9                                               
               DO I = 1, 2 * N + 1                                         
                  K = I - N - 1                                                 
                  VRTDIQ = VRTDIQ + EGT(I,IVTOR1,IVTORP) * K                    
     &                 * A3T5D(I,IVTOR1,IVTORP)                                 
               end do
            end do
         end do
      END IF                                                                    

      IF (PARA(LO,L).EQ.'AK3CK  ') THEN                                          
         VRTDIQ = 0.d0                                                          
         DO IVTORP = 1,9                                                  
            DO IVTOR1 = 1,9                                               
               DO I = 1, 2 * N + 1                                         
                  K = I - N - 1                
                  xk = dble(k)
                  VRTDIQ = VRTDIQ + EGT(I,IVTOR1,IVTORP) * K                    
     &                 * A3T5C(I,IVTOR1,IVTORP) * K**2                                 
               end do
            end do
         end do
      END IF                                                                    

      IF (PARA(LO,L).EQ.'AK3CK4') THEN                                          
         VRTDIQ = 0.d0                                                          
         DO IVTORP = 1,9                                                  
            DO IVTOR1 = 1,9                                               
               DO I = 1, 2 * N + 1                                         
                  K = I - N - 1                            
                  xk = dble(k)
                  VRTDIQ = VRTDIQ + EGT(I,IVTOR1,IVTORP) * xk**5                    
     &                 * A3T5C(I,IVTOR1,IVTORP) 
               end do
            end do
         end do
      END IF                                                                    
*                                                                               
      IF (PARA(LO,L).EQ.'AK4   ' .OR. PARA(LO,L).EQ.'AMV   '.or.
     &     PARA(LO, L).eq.'G4J4  '.or.PARA(LO,L).eq.'G4J2K2'.or.
     &     PARA(LO, L).eq.'G4K4  ') THEN              
         VRTDIG = 0.d0                                                          
         VRTDIH = 0.d0                                                          
         VRTDIL = 0.d0                                                          
         DO IVTORP = 1,9                                                  
            DO IVTOR1 = 1,9                                               
               DO I = 1, 2 * N + 1                                         

                 K = I - N - 1                                                 
                 xk = dble(k)

                  VRTDIG = VRTDIG + EGT(I,IVTOR1,IVTORP)                        
     &                 * A5T5(I,IVTOR1,IVTORP)                                  
                  VRTDIH = VRTDIH + EGT(I,IVTOR1,IVTORP)                        
     &                 * A5T5(I,IVTOR1,IVTORP) * xK**2
                  VRTDIL = VRTDIL + EGT(I,IVTOR1,IVTORP)                        
     &                 * A5T5(I,IVTOR1,IVTORP) * xK**4                           
               end do
            end do
         end do
      END IF                                                                    
*                                                                               
      IF (PARA(LO,L).EQ.'AK4BK4') then
          VRTDIY = 0.d0            
         DO IVTORP = 1,9                                                  
            DO IVTOR1 = 1,9                                               
               DO I = 1, 2 * N + 1                                         
                 K = I - N - 1                                                 
                 xk = dble(k)
                  VRTDIY = VRTDIY + EGT(I,IVTOR1,IVTORP)                        
     &                 * A5T5B(I,IVTOR1,IVTORP) * xK**4
               end do
            end do
         end do
      END IF                                                                    

      IF (PARA(LO,L).EQ.'AK4BK6') then
          VRTDIY = 0.d0            
         DO IVTORP = 1,9                                                  
            DO IVTOR1 = 1,9                                               
               DO I = 1, 2 * N + 1                                         
                 K = I - N - 1                     
                 xk = dble(k)
                  VRTDIY = VRTDIY + EGT(I,IVTOR1,IVTORP)                        
     &                 * A5T5B(I,IVTOR1,IVTORP) * xK**6
               end do
            end do
         end do
      END IF                                                                    

      IF (PARA(LO,L).EQ.'AK4B  '.or.para(lo,l).eq.'AK4BJ '.or.
     &     para(lo, l).eq.'AK4BK ') THEN                                          
         VRTDIY = 0.d0            
         VRTDIQ = 0.d0
         DO IVTORP = 1,9                                                  
            DO IVTOR1 = 1,9                                               
               DO I = 1, 2 * N + 1                                         

                 K = I - N - 1                     
                 xk = dble(k)

                  VRTDIY = VRTDIY + EGT(I,IVTOR1,IVTORP)                        
     &                 * A5T5B(I,IVTOR1,IVTORP)                                 
                  VRTDIQ = VRTDIQ + EGT(I,IVTOR1,IVTORP)                        
     &                 * A5T5B(I,IVTOR1,IVTORP) * xK**2                                 
               end do
            end do
         end do
      END IF                                                                    

      IF (PARA(LO,L).EQ.'AK4C  ') THEN                                          
         VRTDIY = 0.d0            
         DO IVTORP = 1,9                                                  
            DO IVTOR1 = 1,9                                               
               DO I = 1, 2 * N + 1                                         
                 K = I - N - 1                                                 
                  VRTDIY = VRTDIY + EGT(I,IVTOR1,IVTORP)                        
     &                 * A4T5C(I,IVTOR1,IVTORP)                                 
               end do
            end do
         end do
      END IF                                                                    

      IF (PARA(LO,L).EQ.'AK4D  ') THEN                                          
         VRTDIY = 0.d0            
         DO IVTORP = 1,9                                                  
            DO IVTOR1 = 1,9                                               
               DO I = 1, 2 * N + 1                                         
                 K = I - N - 1                                                 
                  VRTDIY = VRTDIY + EGT(I,IVTOR1,IVTORP)                        
     &                 * A4T5D(I,IVTOR1,IVTORP)                                 
               end do
            end do
         end do
      END IF                                                                    

      IF (PARA(LO,L).EQ.'AK4CK ') THEN                                          
         VRTDIY = 0.d0            
         DO IVTORP = 1,9                                                  
            DO IVTOR1 = 1,9                                               
               DO I = 1, 2 * N + 1                                         
                 K = I - N - 1                     
                 xk = dble(k)
                  VRTDIY = VRTDIY + EGT(I,IVTOR1,IVTORP)                        
     &                 * A4T5C(I,IVTOR1,IVTORP) *xK**2
               end do
            end do
         end do
      END IF                                                                    

      IF (PARA(LO,L).EQ.'AK4CK4') THEN                                          
         VRTDIY = 0.d0            
         DO IVTORP = 1,9                                                  
            DO IVTOR1 = 1,9                                               
               DO I = 1, 2 * N + 1                                         
                 K = I - N - 1                                                 
                 xk = dble(k)
                  VRTDIY = VRTDIY + EGT(I,IVTOR1,IVTORP)                        
     &                 * A4T5C(I,IVTOR1,IVTORP) *xK**4
               end do
            end do
         end do
      END IF                                                                    
*                                                                               
C     PRINT*,'PARA=',PARA(LO,L)
      IF (PARA(LO,L).EQ.'BK1   ')THEN                                           
         VRTDIX = 0.d0                                                          
         DO IVTORP = 1,9                                                  
            DO IVTOR1 = 1,9                                               
               DO I = 1, 2 * N + 1                                         
                  K=I-N-1
                  xk = dble(k)
                  VRTDIX = VRTDIX + EGT(I,IVTOR1,IVTORP)                        
     &                 * A5T5(I,IVTOR1,IVTORP)*xK**2                             
               end do
            end do
         end do
      END IF                                                                    
*                                                                               
      IF (PARA(LO,L).EQ.'AK5   ' .OR. PARA(LO,L).EQ.'AK5J  ') THEN              
         VRTDIH = 0.d0                                                          
         DO IVTORP = 1,9                                                  
            DO IVTOR1 = 1,9                                               
               DO I = 1, 2 * N + 1                                         
                  K = I - N - 1                 
                  xk = dble(k)
                  VRTDIH = VRTDIH + EGT(I,IVTOR1,IVTORP) * xk * xk               
     &                 * (A1T1(I,IVTOR1,IVTORP)                                 
     &                    - 0.5d0*A7T5(I,IVTOR1,IVTORP)                         
     &                    - 0.5d0*A1T5(I,IVTOR1,IVTORP))                        
               end do
            end do
         end do
      END IF                                                                    
*                                                                               
      IF (PARA(LO,L).EQ.'AK7   '.OR.PARA(LO,L).EQ.'AK7J  ')THEN                 
         VRTDIJ = 0.d0                                                          
         DO IVTORP = 1,9                                                  
            DO IVTOR1 = 1,9                                               
               DO I = 1, 2 * N + 1                                         
                  K = I - N - 1                                                 
                  VRTDIJ = VRTDIJ + EGT(I,IVTOR1,IVTORP)                        
     &                 * (A2T5(I,IVTOR1,IVTORP)                                 
     &                    - 0.5d0*A10T5(I,IVTOR1,IVTORP)                        
     &                    - 0.5d0*A11T5(I,IVTOR1,IVTORP)                        
     &                    + A2T5(I,IVTOR1,IVTORP)                               
     &                    - 0.5d0*A16T5(I,IVTOR1,IVTORP)                        
     &                    - 0.5d0*A17T5(I,IVTOR1,IVTORP))                       
               end do
            end do
         end do
      END IF                                                                    

      IF (PARA(LO,L).EQ.'AK8   '.or.para(lo,l).eq.'AK8J  '.or.
     &     para(lo,l).eq.'AK8K  ')THEN                 
         VRTDIJ = 0.d0           
         VRTDIK = 0.d0
         DO IVTORP = 1,9                                                  
            DO IVTOR1 = 1,9                                               
               DO I = 1, 2 * N + 1                                         
                  K = I - N - 1                                                 
                  VRTDIJ = VRTDIJ + EGT(I,IVTOR1,IVTORP)                        
     &                 * (A2T58(I,IVTOR1,IVTORP)                                 
     &                    - 0.5d0*A10T58(I,IVTOR1,IVTORP)                        
     &                    - 0.5d0*A11T58(I,IVTOR1,IVTORP)                        
     &                    + A2T58(I,IVTOR1,IVTORP)                               
     &                    - 0.5d0*A16T58(I,IVTOR1,IVTORP)                        
     &                    - 0.5d0*A17T58(I,IVTOR1,IVTORP))                       

                  VRTDIK = VRTDIK + EGT(I,IVTOR1,IVTORP)                        
     &                 * (A2T58(I,IVTOR1,IVTORP)                                 
     &                    - 0.5d0*A10T58(I,IVTOR1,IVTORP)                        
     &                    - 0.5d0*A11T58(I,IVTOR1,IVTORP)                        
     &                    + A2T58(I,IVTOR1,IVTORP)                               
     &                    - 0.5d0*A16T58(I,IVTOR1,IVTORP)                        
     &                    - 0.5d0*A17T58(I,IVTOR1,IVTORP))   
     &                 * K**2
               end do
            end do
         end do
      END IF                                                                    
*                                                                               
      IF (PARA(LO,L).EQ.'AK7K  '.or.para(lo,l).eq.'AK7KK '
     &     .or.para(lo,l).eq.'AK7JK '.or.para(lo,l).eq.'AK7JJ ') THEN                                          
         VRTDIV = 0.d0                                         
         VRTDIU = 0.d0
         VRTDIW = 0.d0
         DO IVTORP = 1,9                                                  
            DO IVTOR1 = 1,9                                               
               DO I = 1, 2 * N + 1                                         
                  K = I - N - 1                                   
                  xk = dble(k)
                  VRTDIV = VRTDIV + EGT(I,IVTOR1,IVTORP)*xk**2
     &                 * (A2T5(I,IVTOR1,IVTORP)                          
     &                    - 0.5d0*A10T5(I,IVTOR1,IVTORP)                        
     &                    - 0.5d0*A11T5(I,IVTOR1,IVTORP)                        
     &                    + A2T5(I,IVTOR1,IVTORP)                               
     &                    - 0.5d0*A16T5(I,IVTOR1,IVTORP)                        
     &                    - 0.5d0*A17T5(I,IVTOR1,IVTORP))                       

                  VRTDIU = VRTDIU + EGT(I,IVTOR1,IVTORP)                  
     &                 * (A2T5(I,IVTOR1,IVTORP)                          
     &                    - 0.5d0*A10T5(I,IVTOR1,IVTORP)                        
     &                    - 0.5d0*A11T5(I,IVTOR1,IVTORP)                        
     &                    + A2T5(I,IVTOR1,IVTORP)                               
     &                    - 0.5d0*A16T5(I,IVTOR1,IVTORP)                        
     &                    - 0.5d0*A17T5(I,IVTOR1,IVTORP))                       

                  VRTDIW = VRTDIW + EGT(I,IVTOR1,IVTORP)*xK**4                  
     &                 * (A2T5(I,IVTOR1,IVTORP)                          
     &                    - 0.5d0*A10T5(I,IVTOR1,IVTORP)                        
     &                    - 0.5d0*A11T5(I,IVTOR1,IVTORP)                        
     &                    + A2T5(I,IVTOR1,IVTORP)                               
     &                    - 0.5d0*A16T5(I,IVTOR1,IVTORP)                        
     &                    - 0.5d0*A17T5(I,IVTOR1,IVTORP))                       

               end do
            end do
         end do
         ENDIF
*
         IF (PARA(LO,L).EQ.'ANV   '.or. PARA(LO, L).eq. 'C6J4  '.or.
     &        PARA(LO,L).eq.'C6J2K2'.or.PARA(LO,L).eq.'C6K4  ') THEN                                          
            VRTDIK = 0.d0                                        
            VRTDIL = 0.d0                                        
            VRTDIM = 0.d0                                        
         
         DO IVTORP = 1,9                                                  
            DO IVTOR1 = 1,9                                               
               DO I = 1, 2 * N + 1                    

                 K = I - N - 1                                                 
                 xk = dble(k)    
                 
                  VRTDIK = VRTDIK + EGT(I,IVTOR1,IVTORP)                        
     &                 * (A1T1(I,IVTOR1,IVTORP)                                 
     &                    - 0.5d0*A12T5(I,IVTOR1,IVTORP)                        
     &                    - 0.5d0*A13T5(I,IVTOR1,IVTORP))                       

                  VRTDIL = VRTDIL + EGT(I,IVTOR1,IVTORP)                        
     &                 * (A1T1(I,IVTOR1,IVTORP)                                 
     &                    - 0.5d0*A12T5(I,IVTOR1,IVTORP)                        
     &                    - 0.5d0*A13T5(I,IVTOR1,IVTORP)) * xK**2

                  VRTDIM = VRTDIM + EGT(I,IVTOR1,IVTORP)                        
     &                 * (A1T1(I,IVTOR1,IVTORP)                                 
     &                    - 0.5d0*A12T5(I,IVTOR1,IVTORP)                        
     &                    - 0.5d0*A13T5(I,IVTOR1,IVTORP)) * xK**4
               end do
            end do
         end do
      END IF                                                                    

      IF (PARA(lo, l).eq.'V9J   '.or.PARA(lo,l).eq.'V9K   ') then
         VRTDIK = 0.d0                                                          
         VRTDIL = 0.d0                                                          
         DO  IVTORP = 1,9                                                  
            DO  IVTOR1 = 1,9                                               
               DO  I = 1, 2 * N + 1                                         

                 K = I - N - 1                               
                 xk = dble(k)

                  VRTDIK = VRTDIK + EGT(I,IVTOR1,IVTORP)                        
     &                 * (A1T1(I,IVTOR1,IVTORP) 
     &                - 0.5d0*A5TS(I,IVTOR1,IVTORP))

                  VRTDIL = VRTDIL + EGT(I,IVTOR1,IVTORP)                        
     &                 * (A1T1(I,IVTOR1,IVTORP) 
     &                 - 0.5d0 * A5TS(I,IVTOR1,IVTORP)) * xK**2
               end do
            end do
         end do
      end if

      IF (PARA(lo, l).eq.'V12J  ') then
         VRTDIK = 0.d0                                                          
         DO  IVTORP = 1,9                                                  
            DO  IVTOR1 = 1,9                                               
               DO  I = 1, 2 * N + 1                                         

                 K = I - N - 1                                                 

                  VRTDIK = VRTDIK + EGT(I,IVTOR1,IVTORP)                        
     &                * (A1T1(I,IVTOR1,IVTORP) - 
     &                0.5d0 *A12TS(I,IVTOR1,IVTORP))

               end do
            end do
         end do
      end if

*                                                                               
      IF (PARA(LO,L).EQ.'BK2   ') THEN                                          
         VRTDIL = 0.d0                                                          
         DO IVTORP = 1,9                                                  
            DO IVTOR1 = 1,9                                               
               DO I = 1, 2 * N + 1                                         
                  K = I - N - 1                            
                  xk = dble(k)
                  VRTDIL = VRTDIL + EGT(I,IVTOR1,IVTORP) * xK * xK                
     &                 * (A1T1(I,IVTOR1,IVTORP)                                 
     &                    - 0.5d0*A12T5(I,IVTOR1,IVTORP)                        
     &                    - 0.5d0*A13T5(I,IVTOR1,IVTORP))                       
               end do
            end do
         end do
      END IF                                                                    
*                                                                               
      IF (PARA(LO,L).EQ.'AK6   '.OR.PARA(LO,L).EQ.'AK6K  '.OR.  
     1     PARA(LO,L).EQ.'AK6J  '.or.PARA(LO,L).eq.'AK6KK '.or.
     2     PARA(LO,L).eq.'AK6JK '.or.PARA(LO,L).eq.'AK6JJ ') THEN
         VRTDIM = 0.d0                                                          
         VRTDIU = 0.d0                                                          
         vrtdiw = 0.d0
         DO IVTORP = 1,9                                                  
            DO IVTOR1 = 1,9                                               
               DO I = 1, 2 * N + 1                                         
                  K = I - N - 1                                                 
                  xk = dble(k)
                  VRTDIM = VRTDIM + EGT(I,IVTOR1,IVTORP) * xK                    
     &                 * (A3T5(I,IVTOR1,IVTORP)                                 
     &                    - 0.5d0*A8T5(I,IVTOR1,IVTORP)                         
     &                    - 0.5d0*A9T5(I,IVTOR1,IVTORP)                         
     &                    + A3T5(I,IVTOR1,IVTORP)                               
     &                    - 0.5d0*A14T5(I,IVTOR1,IVTORP)                        
     &                    - 0.5d0*A15T5(I,IVTOR1,IVTORP))                       
                  VRTDIU = VRTDIU + EGT(I,IVTOR1,IVTORP)*xk**3
     &                 * (A3T5(I,IVTOR1,IVTORP)                                 
     &                    - 0.5d0*A8T5(I,IVTOR1,IVTORP)                         
     &                    - 0.5d0*A9T5(I,IVTOR1,IVTORP)                         
     &                    + A3T5(I,IVTOR1,IVTORP)                               
     &                    - 0.5d0*A14T5(I,IVTOR1,IVTORP)                        
     &                    - 0.5d0*A15T5(I,IVTOR1,IVTORP))                       
                  VRTDIw = VRTDIw + EGT(I,IVTOR1,IVTORP)*xk**5
     &                 * (A3T5(I,IVTOR1,IVTORP)                                 
     &                    - 0.5d0*A8T5(I,IVTOR1,IVTORP)                         
     &                    - 0.5d0*A9T5(I,IVTOR1,IVTORP)                         
     &                    + A3T5(I,IVTOR1,IVTORP)                               
     &                    - 0.5d0*A14T5(I,IVTOR1,IVTORP)                        
     &                    - 0.5d0*A15T5(I,IVTOR1,IVTORP))                       
               end do
            end do
         end do
      END IF                                                                    

      IF (PARA(LO,L).EQ.'AK9   '.OR.PARA(LO,L).EQ.'AK9K  '.OR.  
     1     PARA(LO,L).EQ.'AK9J  ') THEN
         VRTDIM = 0.d0                                                          
         VRTDIU = 0.d0                                                          
         vrtdiw = 0.d0
         DO IVTORP = 1,9                                                  
            DO IVTOR1 = 1,9                                               
               DO I = 1, 2 * N + 1                                         
                  K = I - N - 1                               
                  xk = dble(k)
                  VRTDIM = VRTDIM + EGT(I,IVTOR1,IVTORP) * xK                    
     &                 * (A3T59(I,IVTOR1,IVTORP)                                 
     &                    - 0.5d0*A8T59(I,IVTOR1,IVTORP)                         
     &                    - 0.5d0*A9T59(I,IVTOR1,IVTORP)                         
     &                    + A3T59(I,IVTOR1,IVTORP)                               
     &                    - 0.5d0*A14T59(I,IVTOR1,IVTORP)                        
     &                    - 0.5d0*A15T59(I,IVTOR1,IVTORP))                       
                  VRTDIU = VRTDIU + EGT(I,IVTOR1,IVTORP)*xk**3
     &                 * (A3T59(I,IVTOR1,IVTORP)                                 
     &                    - 0.5d0*A8T59(I,IVTOR1,IVTORP)                         
     &                    - 0.5d0*A9T59(I,IVTOR1,IVTORP)                         
     &                    + A3T59(I,IVTOR1,IVTORP)                               
     &                    - 0.5d0*A14T59(I,IVTOR1,IVTORP)                        
     &                    - 0.5d0*A15T59(I,IVTOR1,IVTORP))                       

               end do
            end do
         end do
      END IF                                                                    

*                                                                               
      IF (PARA(LO,L).EQ.'AK66  ') THEN                                          
         VRTDIW = 0.d0                                                          
         DO IVTORP = 1,9                                                  
            DO IVTOR1 = 1,9                                               
               DO I = 1, 2 * N + 1                                         
                  K = I - N - 1                                                 
                  VRTDIW = VRTDIW + EGT(I,IVTOR1,IVTORP) * K                    
     &                 * (A3T5(I,IVTOR1,IVTORP)                                 
     &                    - 0.5d0*A8T5B(I,IVTOR1,IVTORP)                         
     &                    - 0.5d0*A9T5B(I,IVTOR1,IVTORP)                         
     &                    + A3T5(I,IVTOR1,IVTORP)                               
     &                    - 0.5d0*A14T5B(I,IVTOR1,IVTORP)                        
     &                    - 0.5d0*A15T5B(I,IVTOR1,IVTORP))                       
               end do
            end do
         end do
      END IF                                                                    

      IF (PARA(LO,L).EQ.'AK76  ') THEN                                          
         VRTDIW = 0.d0                                                          
         DO IVTORP = 1,9                                                  
            DO IVTOR1 = 1,9                                               
               DO I = 1, 2 * N + 1                                         
                  K = I - N - 1                                                 
                  VRTDIW = VRTDIW + EGT(I,IVTOR1,IVTORP)
     &                 * (A3T5(I,IVTOR1,IVTORP)                                 
     &                    - 0.5d0*A8T5C(I,IVTOR1,IVTORP)                         
     &                    - 0.5d0*A9T5C(I,IVTOR1,IVTORP)                         
     &                    + A3T5(I,IVTOR1,IVTORP)                               
     &                    - 0.5d0*A14T5C(I,IVTOR1,IVTORP)                        
     &                    - 0.5d0*A15T5C(I,IVTOR1,IVTORP))                       
               end do
            end do
         end do
      END IF                                                                    
*                                                                               
      IF (PARA(LO,L).EQ.'AK5K  ') THEN                                          
         VRTDIN = 0.d0                                                          
         DO IVTORP = 1,9                                                  
            DO IVTOR1 = 1,9                                               
               DO I = 1, 2 * N + 1                                         
                  K = I - N - 1                              
                  xk = dble(k)
                  VRTDIN = VRTDIN + EGT(I,IVTOR1,IVTORP) * xk**4
     &                 * (A1T1(I,IVTOR1,IVTORP)                                 
     &                    - 0.5d0*A7T5(I,IVTOR1,IVTORP)                         
     &                    - 0.5d0*A1T5(I,IVTOR1,IVTORP))                        
               end do
            end do
         end do
      END IF                                                                    

      IF (PARA(LO,L).EQ.'FJ6   '.or.PARA(LO,L).eq.'FJ4K2 '.or.
     &     PARA(LO,L).eq.'FJ2K4 '.or.PARA(LO,L).eq.'FK6   ') THEN                                          
         VRTDIK = 0.d0                                                          
         VRTDIL = 0.d0                                                          
         VRTDIM = 0.d0                                                          
         VRTDIN = 0.d0                                                          
         DO IVTORP = 1,9                                                  
            DO IVTOR1 = 1,9                                               
               DO I = 1, 2 * N + 1                                         
                  K = I - N - 1                                 
                  xk = dble(k)
                  VRTDIK = VRTDIK + EGT(I,IVTOR1,IVTORP)
     &                 * (A1T1(I,IVTOR1,IVTORP)                                 
     &                 - 0.5d0*A7T5(I,IVTOR1,IVTORP)                         
     &                 - 0.5d0*A1T5(I,IVTOR1,IVTORP))                        

                  VRTDIL = VRTDIL + EGT(I,IVTOR1,IVTORP) * xK**2
     &                 * (A1T1(I,IVTOR1,IVTORP)                                 
     &                 - 0.5d0*A7T5(I,IVTOR1,IVTORP)                         
     &                 - 0.5d0*A1T5(I,IVTOR1,IVTORP))                        

                  VRTDIM = VRTDIM + EGT(I,IVTOR1,IVTORP) * xK**4        
     &                 * (A1T1(I,IVTOR1,IVTORP)                                 
     &                 - 0.5d0*A7T5(I,IVTOR1,IVTORP)                         
     &                 - 0.5d0*A1T5(I,IVTOR1,IVTORP))                        

                  VRTDIN = VRTDIN + EGT(I,IVTOR1,IVTORP) * xK**6        
     &                 * (A1T1(I,IVTOR1,IVTORP)                                 
     &                 - 0.5d0*A7T5(I,IVTOR1,IVTORP)                         
     &                 - 0.5d0*A1T5(I,IVTOR1,IVTORP))                        
               end do
            end do
         end do
      END IF                                                                    
                                                                               
      IF (PARA(LO,L).EQ.'AK2K  ' .or. PARA(LO,L).eq.'GJ2K4 ') THEN                                          
         VRTDIO = 0.d0                                                          
         DO IVTORP = 1,9                                                  
            DO IVTOR1 = 1,9                                               
               DO I = 1, 2 * N + 1                                         
                  K = I - N - 1                             
                  xk = dble(k)
                  VRTDIO = VRTDIO + EGT(I,IVTOR1,IVTORP) * xK**4 
     &                 * A2T5(I,IVTOR1,IVTORP)                                  
               end do
            end do
         end do
      END IF                                                                    
*

      IF (PARA(LO,L).EQ.'GK6   ') THEN                                          
         VRTDIO = 0.d0                                                          
         DO IVTORP = 1,9                                                  
            DO IVTOR1 = 1,9                                               
               DO I = 1, 2 * N + 1                                         
                  K = I - N - 1                                                 
                  xk = dble(k)
                  VRTDIO = VRTDIO + EGT(I,IVTOR1,IVTORP) * xK**6
     &                 * A2T5(I,IVTOR1,IVTORP)                                  
               end do
            end do
         end do
      END IF                                                                    
                                                                               
      IF (PARA(LO,L).EQ.'AK1K  ') THEN                                          
         VRTDIP = 0.d0                                                          
         DO IVTORP = 1,9                                                  
            DO IVTOR1 = 1,9                                               
               DO I = 1, 2 * N + 1                                         
                  K = I - N - 1                             
                  xk = dble(k)
                  VRTDIP = VRTDIP + EGT(I,IVTOR1,IVTORP) * (xK ** 5)             
     &                 * A3T5(I,IVTOR1,IVTORP)                                  
               end do
            end do
         end do
      END IF                                                                    
*                                                                               
      IF (PARA(LO,L).EQ.'FV    ') THEN                                          
         VRTDI5 = VRTDIA * xnn
      ELSE IF (PARA(LO,L).EQ.'ANV   ') THEN                                     
         VRTDI5 = VRTDIK * xnn
      ELSE IF (PARA(LO,L).EQ.'C6J4  ') THEN                                     
         VRTDI5 = VRTDIK * xnn**2
      ELSE IF (PARA(LO,L).EQ.'C6J2K2') THEN                                     
         VRTDI5 = VRTDIL * xnn
      ELSE IF (PARA(LO,L).EQ.'C6K4  ') THEN                                     
         VRTDI5 = VRTDIM                                            
      ELSE IF (PARA(LO,L).EQ.'V9J   ') THEN                                     
         VRTDI5 = VRTDIK * xnn
      ELSE IF (PARA(LO,L).EQ.'V12J  ') THEN                                     
         VRTDI5 = VRTDIK * xnn
      ELSE IF (PARA(LO,L).EQ.'V9K   ') THEN                                     
         VRTDI5 = VRTDIL                                              
      ELSE IF (PARA(LO,L).EQ.'BK2   ') THEN                                     
         VRTDI5 = VRTDIL                                                        
      ELSE IF (PARA(LO,L).EQ.'OFV   ') THEN                                     
         VRTDI5 = VRTDIA * xnn**2
      ELSE IF (PARA(LO,L).EQ.'GV    ') THEN                                     
         VRTDI5 = VRTDIB * xnn                                              
      ELSE IF (PARA(LO,L).EQ.'GVJ   ') THEN                                     
         VRTDI5 = VRTDIB * xnn**2
      ELSE IF (PARA(LO,L).EQ.'GJ6   ') THEN                                     
         VRTDI5 = VRTDIB * xnn**3
      ELSE IF (PARA(LO,L).EQ.'AK1   ') THEN                                     
         VRTDI5 = VRTDID                                                        
      ELSE IF (PARA(LO,L).EQ.'AK1J  ') THEN                                     
         VRTDI5 = VRTDID * xnn
      ELSE IF (PARA(LO,L).EQ.'AK1K  ') THEN                                     
         VRTDI5 = VRTDIP                                                        
      ELSE IF (PARA(LO,L).EQ.'AK2   ') THEN                                     
         VRTDI5 = VRTDIE                                                        
      ELSE IF (PARA(LO,L).EQ.'GJ4K2 ') THEN                                     
         VRTDI5 = VRTDIE * xnn**2
      ELSE IF (PARA(LO,L).EQ.'AK2K  ') THEN                                     
         VRTDI5 = VRTDIO                                                        
      ELSE IF (PARA(LO,L).EQ.'GJ2K4 ') THEN                                     
         VRTDI5 = VRTDIO * xnn
      ELSE IF (PARA(LO,L).EQ.'GK6   ') THEN                                     
         VRTDI5 = VRTDIO                                                                                            
      ELSE IF (PARA(LO,L).EQ.'AK2J  ') THEN                                     
         VRTDI5 = VRTDIE * xnn
      ELSE IF (PARA(LO,L).EQ.'AK3   ') THEN                                     
         VRTDI5 = VRTDIF                                                        
      ELSE IF (PARA(LO,L).EQ.'AK3B  ') THEN                                     
         VRTDI5 = VRTDIQ                                                        
      ELSE IF (PARA(LO,L).EQ.'AK3BK6') THEN                                     
         VRTDI5 = VRTDIQ                                                        
      ELSE IF (PARA(LO,L).EQ.'AK3BK4') THEN                                     
         VRTDI5 = VRTDIQ                                                        
      ELSE IF (PARA(LO,L).EQ.'AK3C  ') THEN                                     
         VRTDI5 = VRTDIQ                                                        
      ELSE IF (PARA(LO,L).EQ.'AK3D  ') THEN                                     
         VRTDI5 = VRTDIQ                                                        
      ELSE IF (PARA(LO,L).EQ.'AK3CK ') THEN                                     
         VRTDI5 = VRTDIQ                                                        
      ELSE IF (PARA(LO,L).EQ.'AK3CK4') THEN                                     
         VRTDI5 = VRTDIQ                                                        
      ELSE IF (PARA(LO,L).EQ.'AK3BJ ') THEN                                     
         VRTDI5 = VRTDIQ * xnn
      ELSE IF (PARA(LO,L).EQ.'AK3BK ') THEN                                     
         VRTDI5 = VRTDIW                                                        
      ELSE IF (PARA(LO,L).EQ.'AK3J  ') THEN                                     
         VRTDI5 = VRTDIF * xnn
      ELSE IF (PARA(LO,L).EQ.'GAJ4  ') THEN                                     
         VRTDI5 = VRTDIF * xnn**2
      ELSE IF (PARA(LO,L).EQ.'GAJ2K2') THEN                                     
         VRTDI5 = VRTDIG * xnn
      ELSE IF (PARA(LO,L).EQ.'GAK4  ') THEN                                     
         VRTDI5 = VRTDIH                                              
      ELSE IF (PARA(LO,L).EQ.'AK3K  ') THEN                                     
* Probleme : K * K ,problem solved (30 janvier 92)                              
         VRTDI5 = VRTDIZ                                                        
      ELSE IF (PARA(LO,L).EQ.'AK4   ') THEN                                     
         VRTDI5 = VRTDIG                                                        
      ELSE IF (PARA(LO,L).EQ.'AK4B  ') THEN                                     
         VRTDI5 = VRTDIY                                                        
      ELSE IF (PARA(LO,L).EQ.'AK4BK4') THEN                                     
         VRTDI5 = VRTDIY                                                        
      ELSE IF (PARA(LO,L).EQ.'AK4BK6') THEN                                     
         VRTDI5 = VRTDIY                                                        
      ELSE IF (PARA(LO,L).EQ.'AK4C  ') THEN                                     
         VRTDI5 = VRTDIY                                                        
      ELSE IF (PARA(LO,L).EQ.'AK4D  ') THEN                                     
         VRTDI5 = VRTDIY                                                        
      ELSE IF (PARA(LO,L).EQ.'AK4CK ') THEN                                     
         VRTDI5 = VRTDIY                                                        
      ELSE IF (PARA(LO,L).EQ.'AK4CK4') THEN                                     
         VRTDI5 = VRTDIY                                                        
      ELSE IF (PARA(LO,L).EQ.'AK4BJ ') THEN                                     
         VRTDI5 = VRTDIY * xnn
      ELSE IF (PARA(LO,L).EQ.'AK4BK ') THEN                                     
         VRTDI5 = VRTDIQ                                                        
      ELSE IF (PARA(LO,L).EQ.'AMV   ') THEN                                     
         VRTDI5 = VRTDIG * xnn
      ELSE IF (PARA(LO,L).EQ.'G4J4  ') THEN                                     
         VRTDI5 = VRTDIG * xnn**2
      ELSE IF (PARA(LO,L).EQ.'G4J2K2') THEN                                     
         VRTDI5 = VRTDIH * xnn
      ELSE IF (PARA(LO,L).EQ.'G4K4  ') THEN                                     
         VRTDI5 = VRTDIL                                              
      ELSE IF (PARA(LO,L).EQ.'AK5   ') THEN                                     
         VRTDI5 = VRTDIH                                                        
      ELSE IF (PARA(LO,L).EQ.'AK5J  ') THEN                                     
         VRTDI5 = VRTDIH * xnn
      ELSE IF (PARA(LO,L).EQ.'AK5K  ') THEN                                     
         VRTDI5 = VRTDIN
      else if (PARA(LO,L).eq.'FJ6   ') THEN
         VRTDI5 = VRTDIK * xnn**3
      else if (PARA(LO,L).eq.'FJ4K2 ') THEN
         VRTDI5 = VRTDIL * xnn**2
      else if (PARA(LO,L).eq.'FJ2K4 ') THEN
         VRTDI5 = VRTDIM * xnn
      else if (PARA(LO,L).eq.'FK6   ') THEN
         VRTDI5 = VRTDIN
      ELSE IF (PARA(LO,L).EQ.'AK6   ') THEN                                     
         VRTDI5 = VRTDIM                                                        
      ELSE IF (PARA(LO,L).EQ.'AK6J  ') THEN                                     
         VRTDI5 = VRTDIM*xnn
      ELSE IF (PARA(LO,L).EQ.'AK9   ') THEN                                     
         VRTDI5 = VRTDIM                                                        
      ELSE IF (PARA(LO,L).EQ.'AK9J  ') THEN                                     
         VRTDI5 = VRTDIM*xnn
      ELSE IF (PARA(LO,L).EQ.'AK66  ') THEN                                     
         VRTDI5 = VRTDIW
      ELSE IF (PARA(LO,L).EQ.'AK76  ') THEN                                     
         VRTDI5 = VRTDIW
      ELSE IF (PARA(LO,L).EQ.'AK6K  ') THEN
         VRTDI5 =VRTDIU
      ELSE IF (PARA(LO,L).EQ.'AK9K  ') THEN
         VRTDI5 =VRTDIU
      ELSE IF (PARA(LO,L).EQ.'AK6KK ') THEN
         VRTDI5 =VRTDIW
      ELSE IF (PARA(LO,L).EQ.'AK6JK ') THEN
         VRTDI5 =VRTDIU*xnn
      ELSE IF (PARA(LO,L).EQ.'AK6JJ ') THEN
         VRTDI5 =VRTDIm * xnn**2
      ELSE IF (PARA(LO,L).EQ.'AK7   ') THEN                                     
         VRTDI5 = VRTDIJ                                                        
      ELSE IF (PARA(LO,L).EQ.'AK8   ') THEN                                     
         VRTDI5 = VRTDIJ                                                        
      ELSE IF (PARA(LO,L).EQ.'AK8J  ') THEN                                     
         VRTDI5 = VRTDIJ * xnn
      ELSE IF (PARA(LO,L).EQ.'AK8K  ') THEN                                     
         VRTDI5 = VRTDIK                                                        
      ELSE IF (PARA(LO,L).EQ.'AK7K  ') THEN
         VRTDI5 = VRTDIV
      ELSE IF (PARA(LO,L).EQ.'AK7KK ') THEN
         VRTDI5 = VRTDIW
      ELSE IF (PARA(LO,L).EQ.'AK7JK ') THEN
         VRTDI5 = VRTDIV * xnn
      ELSE IF (PARA(LO,L).EQ.'AK7JJ ') THEN
         VRTDI5 = VRTDIU * xnn**2
      ELSE IF (PARA(LO,L).EQ.'AK7J  ') THEN
         VRTDI5 = VRTDIJ*xnn
      ELSE IF (PARA(LO,L).EQ.'BK1   ') THEN
         VRTDI5 = VRTDIX
C     PRINT*,'VRTDIX=',VRTDIX
      ELSE IF (PARA(LO,L).EQ.'ALV   ') THEN                                     
         VRTDI5 = VRTDIC * xnn                                              
      ELSE IF (PARA(LO,L).EQ.'OLV   ') THEN                                     
         VRTDI5 = VRTDIC * xnn**2

      ELSE IF (PARA(LO,L).EQ.'AGJ6  ') THEN                                     
         VRTDI5 = VRTDIC * xnn**3

      ELSE IF (PARA(LO,L).EQ.'AGJ4K2') THEN                                     
         VRTDI5 = VRTDID * xnn**2

      ELSE IF (PARA(LO,L).EQ.'AGJ2K4') THEN                                     
         VRTDI5 = VRTDIE * xnn

      ELSE IF (PARA(LO,L).EQ.'AGK6  ') THEN                                     
         VRTDI5 = VRTDIF
      END IF                                                                    
      END                                                                       


C**************************************************************                 
      FUNCTION VRTDI6 (L,CST,LO,A)                                              
      implicit none
      save
C**************************************************************                 
C     DERIVATIVES OF THE ELEMENTS OF THE MATRIX HAMILTONIAN                     
C     FOR THE RHO PARAMETER                                                     
C     THE DERIVATIVE FOR C4 IS MISSING HERE.                                    
C     THIS SUBROUTINE IS NOT USEFUL ANYMORE SINCE THE ONLY RHO                  
C     TERM IS IN THE TORSIONAL HAMILTONIEN                                      
C         

      include 'common_dim.f'
      include 'common_eigen.f'
      include 'common_dat.f'
      include 'common_rotor.f'
      include 'common_par.f'
      include 'common_symb.f'
      include 'common_chara.f'
      include 'common_quant.f'
      include 'common_cst.f'

      double precision vrtdi6, a

      DIMENSION A(9,NTORMX,KDMX)                                      

      integer lo, ivtorp, iii, iplus2, ivtor1, i, k, l
      double precision rho, gv, alv, ak1, ak2, ak3, ak4, ak5, ak6
      double precision ak7, c1, c3, c2j, gvj, ak2j, c1j, amv, vrtdia
      double precision vrtdib, a1, a2, a3, a4, rhoeff, rhoksg, rhoks1
      double precision eg, ann, sq

      double precision gj6, gj4k2, gj2k4

      xn = dble(n)
      xnn = dble(n*(n+1))

      RHO=CST(LO,73)                                                            
C     RHOJ=CST(LO,60)                                                           
C     RHOK=CST(LO,61)                                                           
      GV=CST(LO,18)                                                             
      ALV=CST(LO,19)                                                            
      AK1=CST(LO,20)                                                            
      AK2=CST(LO,21)                                                            
      AK3=CST(LO,22)                                                            
      AK4=CST(LO,23)                                                            
      AK5=CST(LO,24)                                                            
      AK6=CST(LO,25)                                                            
      AK7=CST(LO,26)                                                            
      C1=CST(LO,32)                                                             
      C3=CST(LO,34)                                                             
      C2J=CST(LO,37)                                                            
      GVJ=CST(LO,38)                                                            
      AK2J=CST(LO,39)                                                           
      C1J=CST(LO,40)                                                            
      AMV=CST(LO,43)                                                            

      GJ6 = CST(LO, 90)
      GJ4K2 = CST(LO, 91)
      GJ2K4 = CST(LO, 92)
C      GK6 = CST(LO, 93)

C     RHOJ=CST(LO,60)                                                           
C     RHOK=CST(LO,61)                                                           

      VRTDIA =0.d0                                                              
      VRTDIB =0.d0                                                              
      DO IVTORP=1,9                                                           
         DO I=1,2*N+1                                                         
            K=I-N-1     
            xk = dble(k)
            A1=0.d0                                                             
            A2=0.d0                                                             
            A3=0.d0                                                             
            A4=0.d0                                                           
            DO III=1,NDIMTO                                                   
C     RHOEFF=RHO+RHOJ*N*(N+1)+RHOK*xK**2                                         
               RHOEFF=RHO                                                       
               RHOKSG =3.*(III-KTRON1)+ISIG+RHOEFF*K                            
               RHOKS1 =RHOKSG +3.                                               
               IF (III.NE.NDIMTO ) THEN                                         
                  A4=A4+A(IVTORP,III,I)*A(IVTORP,III+1,I)                       
     &                 *(RHOKSG+RHOKS1)                                         
               ENDIF                                                            
               A1=A1+A(IVTORP,III,I)*A(IVTORP,III,I)*RHOKSG                     
               A2=A2+A(IVTORP,III,I)*A(IVTORP,III,I)*RHOKSG*RHOKSG              
               A3=A3+A(IVTORP,III,I)*A(IVTORP,III,I)*RHOKSG*RHOKSG              
     &              *RHOKSG                                                     
            end do
            EG=EGVC3(I,IVTORP)**2                                               
            VRTDIA =VRTDIA + EG*(2.d0*(GV+GVJ*xnn)*K*A1+                    
     *           ALV*xK**2)*xnn                                              
            VRTDIA =VRTDIA +(AK1*xK**4+2.d0*(AK2+AK2J*ann)                   
     *           *xK**3*A1)*EG                                                   
            VRTDIA =VRTDIA +(3.d0*AK3*xK**2*A2 + AK4*K*A3*4.d0+                  
     *           AMV*xnn*K*A3*4.d0)*EG                                      
            VRTDIA =VRTDIA +(AK6*xK**2 +AK7*2.d0*(K*A1-0.5d0*K*A4))*EG           
            
C           VRTDIA = VRTDIA + EG*A1*2.d0 * GJ6   * K    * xnn**3
C           VRTDIA = VRTDIA + EG*A1*2.d0 * GJ4K2 * xK**3 * xnn**2
C           VRTDIA = VRTDIA + EG*A1*2.d0 * GJ2K4 * xK**5 * xnn

         end do
      end do
      
      DO IVTORP=1,9                                                           
         DO I=1,2*N-1                                                        
            IPLUS2=I+2                                                          
            K=I-N-1                                                             
            xk = dble(k)
            DO IVTOR1 =1,9                                                   
               A1=0.d0                                                          
               A2=0.d0                                                          
               DO III=1,NDIMTO                                               
                  RHOKSG =3.d0*(III-KTRON1)+ISIG+RHOEFF*K                       
                  RHOKS1 =RHOKSG +2.d0*RHOEFF                                   
                  A1=A1+A(IVTORP,III,I)*A(IVTOR1,III,IPLUS2)                    
     &                 *(K*RHOKSG +(K+2)*RHOKS1 )                               
                  A2=A2+A(IVTORP,III,I)*A(IVTOR1,III,IPLUS2)                    
     &                 *(2.d0*K*(RHOKSG **3)+2.d0*(K+2)*(RHOKS1 **3))           
               end do
               ANN=dble((N*(N+1)-K*(K+1))*(N*(N+1)-(K+1)*(K+2)))               
               SQ=SQRT(ANN)                                                     
               VRTDIB =VRTDIB +EGVC3(I,IVTORP)*EGVC3(IPLUS2,IVTOR1 )            
     *              *SQ*(C1*A1+C3*A2+C1J*A1*N*(N+1))                            
            end do
         end do
      end do
      VRTDI6 = VRTDIA + VRTDIB                                                  
      RETURN                                                                    
      END                                                                       





