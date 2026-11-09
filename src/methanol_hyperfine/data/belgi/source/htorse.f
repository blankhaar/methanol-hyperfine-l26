C********************************************************                       
      SUBROUTINE HTORSE (K,CST,LO)                                              
      implicit none
      save
C                                                                               
C     SETS UP THE TORSION HAMILTANIAN                                           
C                                                                               
C     THIS SUBROUTINE SETS UP THE TORSION HAMILTONIAN                           
C         

      include 'common_dim.f'
      include 'common_eigen.f'
      include 'common_dat.f'
      include 'common_symb.f'
      include 'common_quant.f'
      include 'common_rotor.f'
      include 'common_chara.f'
      include 'common_cst.f'


      integer lo, ktor, ktorp, i, j, k
      double precision rho, rhorho, v3, v6, v9, rhoeff, fpara, v12
      double precision CORREC
C                                                                               
C     LOOP ON THE BASIS FUNCTIONS                                               
      FPARA=CST(LO,72)                                                              
      RHO=CST(LO,73)                                                            
C     RHOJ=CST(LO,60)                                                           
C     RHOK=CST(LO,61)                                                           
      RHORHO=CST(LO,77)                                                         
      V3=CST(LO,74)                                                             
      V6=CST(LO,75)                                                             
      V9=CST(LO,76)                                                             
      V12 = CST(LO, 146)

C     debug
C      write(6, *) 'V12 ', V12

C     INITIALISATION                                                            
      DO 3 I=1,NDIMTO                                                           
         DO 4 J=1,I                                                             
            H(I,J)=0.d0                                                         
 4       CONTINUE                                                               
 3    CONTINUE                                                                  
C                                                                               
C     KTOR:DLEFT SIDE OF THE MATRIX ELEMENT                                     
C     KTORP:DRIGHT SIDE OF THE MATRIX ELEMENT                                   
C                                                                               
      DO 10 I=1,NDIMTO                                                          
         DO 11 J=1,I                                                            
            KTOR=I - KTRON1                                                     
            KTORP=J - KTRON1                                                    
C     RHOEFF=RHORHO+RHOJ*N*(N+1)+RHOK*K**2                                      
            RHOEFF=RHORHO                                                       
C                                                                               
C     KTOR=KTORP                                                                
C                                                                               
C     PRINT*,' KTOR= ',KTOR,' KTORP= ',KTORP                                    
            IF (KTOR.EQ.KTORP)  THEN                                            
               H(I,J) = FPARA*(3.d0*KTOR+ISIG+RHOEFF*K)**2+V3/2.d0                  
     &              +V6/2.d0+V9/2.d0
               H(I,J) = H(I, J) + V12/2.d0                                           

C     PRINT*,' H= ',H(I,J)                                                      
            ELSE IF(IABS(KTOR-KTORP).EQ.1)  THEN                                
C                                                                               
C     KTOR=KTORP+1 OR KTOR=KTORP-1                                              
C                                                                               
               H(I,J) = -V3/4.d0                                                
            ELSE IF (IABS(KTOR-KTORP).EQ.2) THEN                                
C                                                                               
C     KTOR=KTORP+2 OR KTOR=KTORP-2                                              
C                                                                              
               H(I,J) = -V6/4.d0                                                
            ELSE IF(IABS(KTOR-KTORP).EQ.3) THEN                                 
               H(I,J)=-V9/4.d0                                                  
            else if (IABS(KTOR-KTORP).eq.4) then
               H(I,J) = H(I, J) -V12/4.d0
            ENDIF                                                               
 11      CONTINUE                                                               
 10   CONTINUE                                                                  
C     SYMETRISES THE MATRIX 
      DO 5 I=2,NDIMTO                                                           
         DO 6 J=1,I-1                                                           
            H(J,I) = H(I,J) 
 6       CONTINUE                                                               
 5    CONTINUE                                                                  


C     when doing the H^2 vs. eigenvalue^2 check, there is a small
C     problem with the upper-left and lower-right corners, due to
C     the m truncation. CORREC is a term that fixes this problem
      CORREC = CST(LO, 173)

      H(1, 1) = H(1, 1) + CORREC
      H(NDIMTO, NDIMTO) = H(NDIMTO, NDIMTO) + CORREC

      RETURN                                                                    
      END                                                                       

