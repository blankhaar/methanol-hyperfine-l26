C***********************************************************
      SUBROUTINE ENCAL(E,LO,EGVL,EGVC,A)                                   
      implicit none
      save
C***********************************************************

      include 'common_dim.f'
      include 'common_eigen.f'
      include 'common_rotor.f'
      include 'common_quant.f'
      include 'common_chara.f'
      include 'common_cntl.f'
      include 'common_expdat.f'

      double precision egvl, egvc, egvc4, a, e
      integer ivtorp, iu, ii, lo, ipose, iposa, j, i, sk
      integer encal2
      integer pickup

      DIMENSION EGVL(NDMX),EGVC(NDMX,NDMX),A(9,NTORMX,KDMX)                     
      DIMENSION EGVC4(KDMX,NTORMX,NDMX)

      pickup = 0

C     sk is "signed K"
      if ( isig .eq. 0 ) then 
         if ( pr .eq. '+' ) then
            sk = kk
         else
            sk = -kk
         end if
      else
         sk = kk
      end if

C     the values of "label" tells us which labelling scheme to use
      if ( (label .eq. 1 .or. label .eq. 2) .and. (N .gt. 5) .and.
     *     (kk .ne. N) .and. (kk .ne. -N) ) then
         nbr(lo) = encal2(E, LO, EGVL, EGVC, A, SK)
      else

C     manual fix of the pickup problem
         if ( label .eq. 3 ) then
            if ( iv .eq. 2 .and. kk .eq. -9 .and. n .ge. 19) then
               kk = -kk 
               pickup = 1
            end if

C            if ( iv .eq. 1 .and. kk .eq. -8 .and. n .ge. 18) then
C               kk = -kk 
C               pickup = 1
C            end if
         end if
         
C     TAKE THE GOOD EIGENVALUE(E) AND ASSOCIATE A NUMBER(NBR)                   
C     TO THE EIGENVECTOR CORRESPONDING TO THE K DATA                            
C     
         IF(ISIG.EQ.0) THEN                                                        
            DO 5 J=1,NROTOR
               IVTORP=1                                                      
               DO 4 I=1,NROTOR                                           
                  IU=(2*N+1)*IVTORP+1                                                    
                  IF(I.EQ.IU)                IVTORP=IVTORP+1                             
                  II=I-(2*N+1)*(IVTORP-1)                                                
                  EGVC4(II,IVTORP,J)=EGVC(I,J)                                           
 4             CONTINUE
 5          CONTINUE
            NBR(LO)=IPOSA(EGVC,EGVL,A,EGVC4)                                       
         ELSE                                                                      
            NBR(LO)=IPOSE(EGVC,EGVL)                                               
         ENDIF                                                                     
      endif
C     PRINT*,' NBR= ',NBR(LO),' LO= ',LO                                        
      E=EGVL(NBR(LO))                                                           
      IVTORP=1                                                                  
      DO 3 I=1,NROTOR                                                           
         IU=(2*N+1)*IVTORP+1                                                    
         IF(I.EQ.IU) IVTORP=IVTORP+1                             
         II=I-(2*N+1)*(IVTORP-1)                                                
         EGVC3(II,IVTORP)=EGVC(I,NBR(LO))                                       
         IF(IVTORP.GE.9) EGVC3(II,IVTORP)=0.D0                                  

         EGVC3J(IV+1, ISIG+1, SK+JMAX+1, I) = EGVC3(II, IVTORP)

 3    CONTINUE                                                                  

      EGVLJ(IV+1, ISIG+1, SK+JMAX+1, 1) = 
     &     EGVLJ(IV+1, ISIG+1, SK+JMAX+1, 2)

      EGVLJ(IV+1, ISIG+1, SK+JMAX+1, 2) = 
     &     EGVLJ(IV+1, ISIG+1, SK+JMAX+1, 3)

      EGVLJ(IV+1, ISIG+1, SK+JMAX+1, 3) = E

C     manual fix of the pickup problem
      if ( label .eq. 3 .and. pickup .eq. 1) then
            if ( iv .eq. 2 .and. kk .eq. 9 .and. n .ge. 19) then
               kk = -kk
            end if

C            if ( iv .eq. 1 .and. kk .eq. 8 .and. n .ge. 18) then
C               kk = -kk 
C            end if
      end if

      RETURN                                                                    

      END                                                                       
      

