C************************************************************                   
      SUBROUTINE HRTSET(CST,LO,ETOR,A)                                          
      implicit none
      save
C************************************************************                   
C     SETS UP THE ROTATION-TORSION HAMILTONIAN                                  


      include 'common_dim.f'
      include 'common_eigen.f'
      include 'common_dat.f'
      include 'common_rotor.f'
      include 'common_symb.f'
      include 'common_quant.f'
      include 'common_chara.f'
      include 'common_cst.f'

      double precision etor, h0, h1, h2, h3, h4, a, h5
      integer lo, i, k, j

      DIMENSION A(9,NTORMX,KDMX),ETOR(KDMX,9)                         

      integer ivtor, ivtorp, kp


      DO 3 I=1,NROTOR                                                           
         DO 4 J=1,I                                                             
C     INITIALISATION                                                            
            H(I,J)=0.d0                                                         
 4       CONTINUE                                                               
 3    CONTINUE                                                                  
      DO 7 I=1,NROTOR                                                           
         DO 8 J=1,I                                                             
C     LEFT-SIDE MATRIX ELEMENT                                                  
            IVTOR=(I-1)/(2*N+1)+1                                               
            K=I+N*(1 -2*IVTOR)-IVTOR                                            
C     RIGHT-SIDE MATRIX ELEMENT                                                 
            IVTORP=(J-1)/(2*N+1)+1                                              
            KP=J+N*(1 -2*IVTORP)-IVTORP                                         
C                                                                               
C     K=KP                                                                      
C     PRINT*,' K= ',K,' KP= ',KP,' IVTOR= ',IVTOR,' IVTORP= ',                  
C     1IVTORP                                                                   
            IF (K.EQ.KP) THEN                                                   
C     IF(IVTOR.EQ.IVTORP) THEN                                                  
               H(I,J)=       H0(K,IVTOR,KP,IVTORP,CST,LO,ETOR,A)                
     1              +       H4(K,IVTOR,KP,IVTORP,CST,LO,ETOR,A)                 
C     PRINT*,' (I,J)= ',H(I,J),' I= ',I,' J= ',J                                
C     PRINT*,'H0=',H0(K,IVTOR,KP,IVTORP,CST,LO,ETOR,A),                         
C     1'K=',K,'IVTOR=',IVTOR,'KP=',KP,'IVTORP=',IVTORP,                         
C     2'LO=',LO                                                                 
C     PRINT*,'H4=',H4(K,IVTOR,KP,IVTORP,CST,LO,ETOR,A)                          
C     ENDIF                                                                     
            ELSE IF (IABS(K-KP).EQ.1) THEN                                      
C              IF     (IABS(IVTOR-IVTORP).LT.9) THEN                            
               H(I,J)=       H2(K,IVTOR,KP,IVTORP,CST,LO,ETOR,A)             
     1              +       H1(K,IVTOR,KP,IVTORP,CST,LO,ETOR,A)              
C     PRINT*,' H(I,J)= ',H(I,J),' I= ',I,' J= ',J                               
C     PRINT*,'H2=',H2(K,IVTOR,KP,IVTORP,CST,LO,ETOR,A)                          
C     PRINT*,'H1=',H1(K,IVTOR,KP,IVTORP,CST,LO,ETOR,A)                          
C     ENDIF                                                            
            ELSE IF(IABS(K-KP).EQ.2) THEN                                       
C     IF (IABS(IVTOR-IVTORP).LT.9) THEN                                
               H(I,J)=       H3(K,IVTOR,KP,IVTORP,CST,LO,ETOR,A)             
C     PRINT*,' H(I,J)= ',H(I,J),' I= ',I,' J= ',J                               
C     PRINT*,'H3=',H3(K,IVTOR,KP,IVTORP,CST,LO,ETOR,A)                          
C     ENDIF                                                            
               
            ELSE IF (IABS(K-KP).eq.4) then
               H(I,J) = H5(K,IVTOR,KP,IVTORP,CST,LO,ETOR,A)
               
            ENDIF                                                               
 8       CONTINUE                                                               

C         call expec(IVTOR, K, A)

 7    CONTINUE                                                                  
      DO 20 I=2,NROTOR                                                          
         DO 21 J=1,I-1                                                          
            H(J,I)=H(I,J)                                                       
 21      CONTINUE                                                               
 20   CONTINUE                                                                  
      RETURN                                                                    
      END                                                                       
