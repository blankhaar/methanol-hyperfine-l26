C *********************************************************                 
      SUBROUTINE REJEC(REJ)                                                     
      implicit none
      save
C                                                                               
C     THIS SUBROUTINE REJECTS AN EXPERIMENTAL VALUE IF THE                      
C     CORRESPONDING TME IS GREATER THAN 3*SIGMA (S).                            
C *********************************************************                 

      include 'common_dim.f'
      include 'common_dat.f'
      include 'common_par.f'
      include 'common_expdat.f'
      include 'common_cntl.f'
      include 'common_quant.f'
      include 'common_chara.f'

      INTEGER REJ                                                               
      double precision rejlim
      integer irejl, itme, i
C                                                                               
      REJLIM = 3*S                                                              
      WRITE(6,803) REJLIM                                                       
      IREJL = 10000 * REJLIM                                                    
      DO 88 I = 1, NDATA                                                        
         IF (IAST(I) .EQ. 0) THEN                                               
            ITME = 10000 * ABS(TME(I))                                          
            IF (ITME .GT. IREJL) THEN                                           
               IAST(I) = 2                                                      
               WRITE(6,800) I,ETRANS(I)                                         
               REJ = 1                                                          
            END IF                                                              
         END IF                                                                 
 88   CONTINUE                                                                  
 800  FORMAT(' LINE NO ',I5,' AT',F14.4,' HAS BEEN REJECTED')                   
 803  FORMAT(' THE REJECTION LIMIT IS ',F10.4)                                  
      RETURN                                                                    
      END                                                                       
