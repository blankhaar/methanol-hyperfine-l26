C***********************************************************                    
      SUBROUTINE RESUM1(LO,EGVL)                                                
      implicit none
      save
C**********************************************************                     
C     WRITES THE ENERGIES OF TORSION LEVELS                                     
C                                                                               
C                                                                               

      include 'common_dim.f'
      include 'common_quant.f'
      include 'common_chara.f'
      include 'common_rotor.f'

      integer lo, i
      double precision egvl
      DIMENSION EGVL(NDMX)                                            

      DO 1 I=1,NDIMTO                                                           
         WRITE(6,53) I-1,EGVL(I)                                                
 53      FORMAT(/,' VT= ',I3,' E= ',F16.10,' CM-1')                             
 1    CONTINUE                                                                  
      RETURN                                                                    
      END                                                                       
C***********************************************************                    
      SUBROUTINE RESUM2(LO,EGVL,EGVC)                                           
      implicit none
      save
C***********************************************************                    
C     WRITES THE ENERGIES OF THE ROTATION-TORSION LEVELS(E STATES)              
C         

      include 'common_dim.f'
      include 'common_dat.f'
      include 'common_quant.f'
      include 'common_chara.f'
      
      integer lo, ipst, nn, nnt, ivt, in, kkk, i, k, j
      double precision egvl, egvc, cst, cc, c
      character pur

      DIMENSION IPST(3*kdmx),EGVL(NDMX),EGVC(NDMX,NDMX)                   
      dimension pur(3*kdmx)


      WRITE(6,*)' N= ',N,' ISIG= ',ISIG                                         
      WRITE(6,50)                                                               
 50   FORMAT(/,' TYPE E ',/)                                                    
      NN=2*N+1                                                                  
      NNT=NN*4                                                                  
      DO J=1,NNT                                                              
         CST=0.d0                                                               
         DO I=1,NNT                                                           
            C=EGVC(I,J)**2                                                      
            IF(C.GT.CST) THEN                                                   
               CST=C                                                            
               IPST(J)=I                                                        
            ENDIF                                                               
         end do
C     TEST THE PURITY                                                           
         PUR(J)=' '                                                             
         DO K=1,NNT                                                           
            IF(K.EQ.IPST(J)) cycle                                           
            CC=EGVC(K,J)**2                                                     
            IF(CC.GE.(CST*0.7D0)) THEN                                          
               PUR(J)='*'                                                       
               exit                                                        
            ENDIF                                                               
         end do
         IVT=0                                                                  
         DO KKK=1,(04)                                                        
            IN=NN*KKK                                                           
            IF (IPST(J).GT.IN) IVT=IVT+1                                        
         end do
         K=IPST(J)-(N+1)-IVT*NN                                                 
         EGVL(J)=EGVL(J)+EVIB(LO)                                               
         WRITE(6,55) IVT,K,EGVL(J),N,PUR(J)                                     
         WRITE(LO+32,55) IVT,K,EGVL(J),N,PUR(J)                                 
 55      FORMAT  (' VT= ',I3,' K= ',I3,' E= ',F16.10,' CM-1 ',' N= ',I3,        
     &        '  PURITY= ',A1)                                                  
      end do
      RETURN                                                                    
      END                                                                       
C*************************************************************                  
      SUBROUTINE RESUM3(LO,EGVL,EGVC)                                           
      implicit none
      save
C**************************************************************                 
C     WRITES THE ENERGIES OF THE ROT-TORSION LEVELS (A STATES)                  
C         

      include 'common_dim.f'
      include 'common_dat.f'
      include 'common_quant.f'
      include 'common_chara.f'
      
      double precision egvl, egvc
      integer nn, nnt, ivt, kkk, in, lo, k, j
      character par

      DIMENSION EGVL(NDMX),EGVC(NDMX,NDMX)                                      
      WRITE(6,*) ' N= ',N,' ISIG= ',ISIG                                        
      WRITE(6,50)                                                               
 50   FORMAT(/,' TYPE A ',/)                                                    
      NN=2*N+1                                                                  
      NNT=NN*4                                                                  
      DO 1 J=1,NNT                                                              
         IVT=0                                                                  
         DO 3 KKK=1,(04)                                                        
            IN=NN*KKK                                                           
            IF (J.GT.IN)  IVT=IVT+1                                             
 3       CONTINUE                                                               
         K=(N*J+N-IVT*N*NN)/NN                                                  
         IF (K.EQ.0) THEN                                                       
            PAR= '+'                                                            
         ELSE IF (MOD(K,2).NE.0) THEN                                           
            IF(MOD(IVT,2).EQ.0) THEN                                            
               IF(MOD(J,2).EQ.0) THEN                                           
                  PAR='+'                                                       
               ELSE                                                             
                  PAR='-'                                                       
               ENDIF                                                            
            ELSE                                                                
               IF(MOD(J,2).EQ.0) THEN                                           
                  PAR='-'                                                       
               ELSE                                                             
                  PAR='+'                                                       
               ENDIF                                                            
            ENDIF                                                               
         ELSE IF(MOD(K,2).EQ.0) THEN                                            
            IF(MOD(IVT,2).EQ.0) THEN                                            
               IF(MOD(J,2).EQ.0) THEN                                           
                  PAR='-'                                                       
               ELSE                                                             
                  PAR='+'                                                       
               ENDIF                                                            
            ELSE                                                                
               IF(MOD(J,2).EQ.0) THEN                                           
                  PAR='+'                                                       
               ELSE                                                             
                  PAR='-'                                                       
               ENDIF                                                            
            ENDIF                                                               
         ENDIF                                                                  
         EGVL(J)=EGVL(J)+EVIB(LO)                                               
         WRITE(6,56) IVT,K,EGVL(J),N,PAR                                        
         WRITE(LO+30,56) IVT,K,EGVL(J),N,PAR                                    
 56      FORMAT  (' VT= ',I3,' K= ',I3,' E= ',F16.10,' CM-1 ',' N= ',I3,        
     &        '  PAR= ',A1)                                                     
 1    CONTINUE                                                                  
      RETURN                                                                    
      END                                                                       
