                                                                
      SUBROUTINE ASET(K,EGVL,EGVC,ETOR,A)                                       
      implicit none
      save
C************************************************************                   
C     STORAGE OF THE COEFICIENTS OF THE TORSION EIGENVECTORS                    
C         

      include 'common_dim.f'
      double precision egvl, egvc, etor, a
      DIMENSION EGVL(NDMX),EGVC(NDMX,NDMX),ETOR(KDMX,9)                         

      include 'common_dat.f'
      include 'common_rotor.f'
      include 'common_quant.f'
      include 'common_expdat.f'
      include 'common_chara.f'
      include 'common_tor59.f'

      DIMENSION A(9,NTORMX,KDMX)                                                
      integer k, i, j, ik, np1


C      call A59SET(K,EGVL,EGVC,ETOR,A59)

      IK=K+N   +1           
      NP1 = N + 1
      IISIG = ISIG+1
                                                    
      DO 1 I=1,(09)                                                             
         DO 2 J=1,NDIMTO                                                        
            A(I,J,IK)=EGVC(J,I)                                                 
            A59(I, J, IK, NP1, IISIG) = EGVC(J, I)
 2       CONTINUE                                                               
         ETOR(IK,I)=EGVL(I)                                                     
C     PRINT*,'ETOR=',ETOR(IK,I),'IK=',IK,'I=',I                                 
 1    CONTINUE                                                                  
      IF(K.EQ.-N)    return                                                   
      DO 3 I=1,(09)                                                             
         IF((A(I,10,IK)/A(I,10,IK-1)).LT.0.D0) THEN                             
            DO 4 J=1,NDIMTO                                                     
               A(I,J,IK)=-A(I,J,IK)                                             
               A59(I,J,IK,NP1,IISIG) = -A59(I,J,IK,NP1,IISIG)
 4          CONTINUE                                                            
         ENDIF                                                                  
 3    CONTINUE                                                                  
 16   RETURN                                                                    
      END                                                                       

