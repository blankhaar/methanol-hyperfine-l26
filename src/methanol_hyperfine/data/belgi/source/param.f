C************************************************************                   
      SUBROUTINE PARAM(I,CST)                                                   
      implicit none
      save
C************************************************************                   
C     READING OF THE PARAMETERS TO BE FITTED AND STORING THEM IN                
C     X VECTOR(ONLY ONE VECTOR FOR THE 2 STATES) PARA(LO,L)                     
C     RETAINS WHICH PARAMETER                                                   
C         

      include 'common_dim.f'
      include 'common_par.f'
      include 'common_chara.f'
      include 'common_symb.f'
      include 'common_cst.f'

      integer jj, ik, irep, i, k, j

      READ(5,20) NP(I)                                                          
 20   FORMAT(I3)                                                                
      WRITE(6,*) ' NBRE OF PARAMETERS TO FIT ',NP(I)                            
      IF(NP(I).NE.0) THEN                                                       
         K=0                                                                    
         IF(I.EQ.2) K=NP(1)                                                     
         WRITE(6,*) ' PARAMETER TO BE FIT '                                     
         DO 2 J=1,NP(I)                                                         
            IF(I.EQ.2.AND.NP(1).NE.0) THEN                                      
               JJ=J+NP(1)                                                       
            ELSE                                                                
               JJ=J                                                             
            ENDIF                                                               
            READ(5,'(A6)') PARA(I,JJ)                                           
            WRITE(6,*) ' PARA = ',PARA(I,JJ)                                    
            IK=IREP(PARA(I,JJ))                                                 
            IND(I,J)=IK                                                         
            X(K+J)=CST(I,IK)                                                    
 2       CONTINUE                                                               
      ENDIF                                                                     
C     PRINT*,' NP(1)',NP(1),' NP(2)= ',NP(2)                                    
      RETURN                                                                    
      END                                                                       
