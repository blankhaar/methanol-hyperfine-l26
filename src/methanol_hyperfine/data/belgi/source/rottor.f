C***********************************************************                    
      SUBROUTINE ROTTOR(CST,LO,EGVL,EGVC,ETOR,A)                                
      implicit none
      save
C***********************************************************                    
C     CALCULATES THE ENERGIES OF LEVELS AND TRANSITIONS FOR EACH                
C     N,LO BLOCK.SAME SUBROUTINE AS THE ROTOR PROGRAM                           
C                  


      include 'common_dim.f'

      include 'common_eigen.f'
      include 'common_dat.f'
      include 'common_rotor.f'
      include 'common_symb.f'
      include 'common_quant.f'
      include 'common_chara.f'
      include 'common_cst.f'

      double precision egvl, egvc, etor, a
      integer lo, ierr, k, kkk

      DIMENSION EGVL(NDMX)                  
      DIMENSION EGVC(NDMX,NDMX),ETOR(KDMX,9),A(9,NTORMX,KDMX)                   


C                                                                               
      KTRON1 =KTRONC+1                                                          
      NDIMTO =2*KTRONC+1                                                        
C     NDIMTOR:DEFFECTIVE DIMENSION OF THE TORSION MATRIX                        
C     PRINT*,' KTRONC1= ',KTRONC1,' NDIMTOR= ',NDIMTOR                          
      IF (NDIMTO .GT.NTORMX) THEN                                               
         WRITE(6,*) ' DIMENSION NTORMX TOO LOW '                                
         STOP                                                                   
      ENDIF                                                                     
      DO K=-N,N                                                               
         IF(IBUG3.NE.0) WRITE(6,60) ISIG,K                                      
         CALL HTORSE (K,CST,LO)                                                 
 60      FORMAT(/,' SIGMA= ',I3,' K= ',I3,/)                                    
         IF (IBUG1.NE.0) THEN                                                   
            WRITE(6,63)                                                         
 63         FORMAT(/,' TORSION MATRIX TO BE DIAGONALISED',/)                    
            CALL MATOUT(H,NDMX,NDMX,NDIMTO ,NDIMTO )                            
         ENDIF                                                                  
C                                                                               
C     PERFORMS DIAGONALISATION                                                  
C                                                                               

C         CALL TRED2(NDMX, NDIMTO,H,EGVL,WORK,EGVC)                              
C         CALL TQL2(NDMX, NDIMTO,EGVL,WORK,EGVC,IERR)                 

C         print *, 'debug 1: ndimto'
         call dcopy(ndmx*ndmx, H, 1, EGVC, 1)

 123     format('DSYEVD ', I1, 1X, F20.0, 1X, I20, 1X, I20, 1X, I20)

C         call dsyevd('V', 'U', ndimto, EGVC, ndmx, egvl, 
C     *        work, lwork, iwork, liwork, ierr)

C         if ( work(1) .gt. lwork .or. iwork(1) .gt. liwork ) then
C            write(6, 123) 1, work(1), iwork(1), lwork, liwork
C         endif

         call dsyev('V', 'U', ndimto, EGVC, ndmx, egvl, 
     *        work, lwork, ierr)

         IF (IERR.NE.0) THEN                                                    
            WRITE(6,62)                                                         
 62         FORMAT(/,' DIAGONALISATION FAILED FOR TORSION ',/)                  
            write(6,*) 'IERR ', IERR
            STOP    

         
C     debug
C         else
C            write (6, *) 'DIAGONALISATION SUCCEEDED FOR TORSION '
         ENDIF                                                                  



         
         IF(IBUG3.NE.0) WRITE(6,81)                                             
 81      FORMAT(//,' TORSION EIGENVALUES AND EIGENVECTORS ',//)                 
         IF (IBUG1.NE.0) THEN                                                   
            CALL MATOU2 (EGVC,EGVL,NDMX,NDMX,NDIMTO ,NDIMTO )                   
         ENDIF                                                                  
!         IF(IBUG3.NE.0) CALL RESUM1(LO,EGVL)                                    
C                                                                               
C     STOCKS THE EIGENVECTORS IN A 3-DIMENSION MATRIX A                         
C                                                                               
         CALL ASET(K,EGVL,EGVC,ETOR,A)                                          
C                                                                               
C     END OF THE LOOP ON K                                                      
C     
      end do
      IF(IBUG2.NE.0) WRITE(6,70) ISIG,N                                         
 70   FORMAT(/,' SIGMA= ',I3,' N= ',I3,/)                                       
      NROTOR=(2*N+1)*9                                                          
C     PRINT*,'NROTOR=',NROTOR                                                   
      CALL HRTSET(CST,LO,ETOR,A)                                                
C     IF (IBUG2.NE.0) THEN                                                      
C        WRITE(6,73)                                                            
C73      FORMAT(/,' ROTATION-TORSION MATRIX TO BE                               
C    1        DIAGONALISED ',/)                                                 
C        CALL MATOUT(H,NDMX,NDMX,NROTOR,NROTOR)                                 
C     ENDIF                                                                     
C                                                                               
C     PERFORMS DIAGONALISATION                                                  
C                                                                               

C      CALL TRED2(NDMX,NROTOR,H,EGVL,WORK,EGVC)                                  
C      CALL TQL2(NDMX,NROTOR,EGVL,WORK,EGVC,IERR)                                


C     print *, 'debug 2: nrotor'
      call dcopy(ndmx*ndmx, H, 1, EGVC, 1)

C     dsyev can calculate the optimal work size
C     for debugging/speed purposes, bail out
C     if lwork is too small, since this can be fixed
C     by resizing and then compiling
C      call dsyev('V', 'U', nrotor, EGVC, ndmx, egvl, 
C     *     work, -1, ierr)

C     optimal LWORK is stored in work(1)
C      if ( lwork .lt. work(1) ) then
C         write (*,*) 'LWORK is too small'
C         stop
C      end if

C     continue as scheduled
      call dsyev('V', 'U', nrotor, EGVC, ndmx, egvl, 
     *     work, lwork, ierr)

C     LAPACK has an alternate divide and conquer routine
C     this should be faster, but if this causes any error
C     (floating point) then we will hace to go back to
C     dsyev
C      call dsyevd('V', 'U', nrotor, EGVC, ndmx, egvl, 
C     *     work, lwork, iwork, liwork, ierr)

C      if ( work(1) .gt. lwork .or. iwork(1) .gt. liwork ) then
C         write(6, 123) 2, work(1), iwork(1), lwork, liwork
C      endif

      IF (IERR.NE.0) THEN                                                       
         WRITE(6,72)                                                            
 72      FORMAT(/,' DIAGONALISATION FAILED FOR ROTATION-                        
     1        TORSION ',/)                                                      
         STOP                                                                   
      ENDIF                                                                     


C     energy without assignments
      if ( ibug5 .eq. 2 ) then
         do kkk = 1, nrotor
            write(6, 163) N, isig, kkk, egvl(kkk)
 163        format('N =  ', I4, 'ISIG = ', I4, 
     &           '  kkk = ', I4, '     ', F24.16)
         end do
      endif
      
C     IF(IBUG3.NE.0) WRITE(6,71)                                                
C**** ibug3 -> ibug2  jth march 1994
      IF(IBUG2.NE.0) WRITE(6,71)                                                
C     IF (IBUG2.NE.0) THEN                                                      
      IF (IBUG2.NE.0.AND.N.EQ.1) THEN                                          
      PRINT*,' N= ',N,'ISIG=',ISIG                                              
         CALL MATOU2 (EGVC,EGVL,NDMX,NDMX,NROTOR,NROTOR)                        
      ENDIF                                                                     
 71   FORMAT(//,' ROTATION-TORSION EIGENVALUES AND                              
     1     EIGENVECTORS ',//)                                                   
      IF(IBUG3.NE.0.OR.IBUG4.NE.0) THEN                                         
         IF(ISIG.EQ.1) THEN                                                     
            CALL RESUM2(LO,EGVL,EGVC)                                           
         ELSE                                                                   
            CALL RESUM3(LO,EGVL,EGVC)                                           
         ENDIF                                                                  
      ENDIF                                                                     
      REWIND 31                                                                 
      REWIND 32                                                                 
      REWIND 33                                                                 
      REWIND 34                                                                 



      END                                                                                 
