C************************************************************                   
C     SUBROUTINE CNTRL(NPP,DELX,DELTAXL,CST)                                    
C************************************************************                   
C     MONITORS THE FITTING(CONCONVERGENCE,REJECTION OF TRANSITIONS,             
C     DIVERGENCE)                                                               
C                                                                               
C************************************************************                   
      SUBROUTINE CNTRL(NPP,DELX,DELTAX ,CST)                                    
      implicit none
      save
C************************************************************                   
C     MONITORS THE FITTING(CONCONVERGENCE,REJECTION OF TRANSITIONS,             
C     DIVERGENCE)                                                               
C         

      include 'common_dim.f'
                                            

      INTEGER REJ                                                               

      include 'common_dat.f'
      include 'common_par.f'
      include 'common_expdat.f'
      include 'common_cntl.f'
      include 'common_chara.f'
      include 'common_cst.f'

      double precision delx, deltax

      DIMENSION DELX(maxparam),DELTAX (maxparam)

      double precision rs, ri
      integer npp, ii, i, j


      IF(NIT.EQ.NITT) THEN                                                      
         IEND=0                                                                 
C     WRITE(6,'("FINISHED",I4," PRESET ITERATIONS")') NIT                       
         RETURN                                                                 
      ENDIF                                                                     

C     THE CONVERGENCE CRITERION HAS BEEN CHANGED (NOV.1990).                    
CC    I = 1                                                                     

      IF(NIT.EQ.1) THEN                                                         
         IEND=1                                                                 
         WRITE(6,307) NIT                                                       
         NIT=NIT+1                                                              
         RETURN                                                                 
      ENDIF                                                                     


      RS=DELTAS/S                                                               
      DO 6 I=1,NPP                                                              
         RI=ABS(DELX(I))/DELTAX (I)                                             
         IF(RI.GT.1.E-3) THEN                                                   
            IF(RS.GT.1.E-1) THEN                                                
CC6   IF (ABS(DELX(I)) .GE. DELTAXL(I)) THEN                                    
C     IF (DELTAS .GT. 0 .AND. NIT .GE. 4) THEN                                  
               WRITE(6,*) ' PROGRAM IS DIVERGING ...'                           
C     commented by JMF
C               IF(NITT.LE.8) GO TO 54                                           
C               STOP                                                             
 54            CONTINUE                                                         
            END IF                                                              
CC    IF(NIT.GT.NITT) STOP                                                      
            IEND = 1                                                            
            WRITE(6,307) NIT                                                    
            NIT = NIT + 1                                                       
            RETURN                                                              
         END IF                                                                 
CC    I = I + 1                                                                 
CC    IF (I .LE. NPP) GO TO 6                                                   
C                                                                               
 6    CONTINUE                                                                  
      REJ = 0                                                                   
C     CALL REJEC(TME,S,REJ)                                                     
      IF (REJ .EQ. 0) THEN                                                      
         WRITE(6,*) ' PROGRAM HAS    CONVERGED AFTER ',NIT,                     
     1        ' ITERATIONS'                                                     
         IEND = 0                                                               
         RETURN                                                                 
      ELSE                                                                      
         DO 94 I = 1,NP12                                                       
            IF (I .GT. NP(1)) THEN                                              
               II = 2                                                           
               J = I - NP(1)                                                    
            ELSE                                                                
               II = 1                                                           
               J = I                                                            
            END IF                                                              
            CST(II,IND(II,J)) = X(I)                                            
 94      CONTINUE                                                               
         TERM = TERMST                                                          
         IEND = -1                                                              
         WRITE(6,807)                                                           
         RETURN                                                                 
      END IF                                                                    
 307  FORMAT(/,' *** END OF CYCLE NO. ',I3,' *** ')                             
 807  FORMAT(//,' *** NEW LEAST SQUARE FITTING PROCEDURE ***',                  
     1     /,' ==========================================')                     
      END                                                                       
