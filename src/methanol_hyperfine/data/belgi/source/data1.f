      subroutine data1(cst)                                                     
C     *************************************************************************
C     This subroutine reads the experimental values and calls two
C     subroutines: data11, which reads the values of the constants, and param,
C     which reads which parameters must be fitted.
C     *************************************************************************
      implicit none
      save


C     READ THE CONSTANT VALUES AND DECIDE WHICH PARAMETERS TO FIT               
C     I=1  LOWER STATE;I=2  UPPER STATE                                         
C     

      include 'common_dim.f'
      include 'common_dat.f'
      include 'common_par.f'
      include 'common_expdat.f'
      include 'common_symb.f'
      include 'common_chara.f'
      include 'common_cst.f'
      include 'common_int.f'
      include 'common_cntl.f'

      integer ii, i, ni, iquit
      character*64 infile

c      call getarg(1, infile)
c      open(unit=5, file=infile)
       open(unit=5,file='input.txt',access='sequential',status='old') 

      WRITE(6,*) ' This output is from the program `globalfit`'
      READ(5,*) EPSI,IBUG1,IBUG2,IBUG3,IBGTME,IBGVTD ,IBUG4, ibug5
      WRITE(6,*) ' EPSI= ',EPSI,' IBUG1= ',IBUG1,' IBUG2= ',                    
     *     IBUG2,' IBUG3= ',IBUG3,' IBGTME= ',IBGTME,' IBGVTDV= ',              
     *     IBGVTD ,'IBUG4= ',IBUG4,'IBUG5=',IBUG5                                              
c      READ(5,20) KTRONC                                                         
c      READ(5,20) IRMW  
       read (5,*) ktronc
       read (5,*) irmw                                                      
 20   FORMAT(I3)                                                                
c      WRITE(6,21)KTRONC,IRMW  
      write(6,*) KTRONC,IRMW                                                  
 21   FORMAT(' KTRONC = ',I3,' IRMW= ',I3)                                      
c      IF(KTRONC.GT.(ktroncmax)) THEN  
       IF(KTRONC.GT.(10)) THEN                                                 
         WRITE(6,*) ' DIMENSION KTRONC TOO LOW '                                
         stop                                                                   
      endif                                                                     
c     the next card was added for nbs computer.                                 
      evib(1)=0.d0                                                              
      evib(2)=0.d0                                                              
      if(irmw.eq.1) then                                                        
         do i=1,1                                                             
            call data11(i,cst)                                                  
            read (5,*) evib(i)                                                  
            do  ii=1,cstm
               cst(2,ii)=cst(1,ii)
            end do
            evib(2)=evib(1)
            call param(i,cst)                                                   
c     print*,' np apres param= ',np(i)                                          
         end do
      elseif(irmw.eq.2) then                                                    
         do i=1,2                                                             
            CALL DATA11(I,CST)                                                  
C     PRINT*, 'I= ',I,' IRMW= ',IRMW                                            
            READ (5,*) EVIB(I)                                                  
            WRITE(6,*) ' Evib= ',EVIB(I),'I= ',I                                
            PRINT*,' Evib(1)= ',EVIB(1),' Evib(2)= ',EVIB(2)                    
            CALL PARAM(I,CST)                                                   
C     PRINT*,' NP APRES PARAM= ',NP(I)                                          
         end do
      ENDIF                                                                     
      int = 0
      unc = 0
      dxmul = 1.0
C      READ(5,*) ITO,NITT,IKAKC,JMAX,vtmax
      READ(5,*) ITO
      READ(5,*) NITT
      READ(5,*) IKAKC
      READ(5,*) JMAX
      READ(5,*) vtmax
      READ(5,*) INT
      READ(5,*) TPR
      READ(5,*) unc
      READ(5,*) dxmul
      READ(5,*) en0
      READ(5,*) label
C      label = 0

      WRITE(6,*) ' Vibrational energy to be fitted? ITO = ', ITO                          
      WRITE(6,*)' Number of iterations: ',NITT                                  
      WRITE(6,*)' Maximum value of J: ',JMAX                                    
      write(6,*)' Maxmum value of vt: ', vtmax
      WRITE(6,*)' Intensity? ',INT
      write(6,*)' Uncertainty? ', UNC
      WRITE(6,*)' Temperature = ',TPR
      write(6,*)' en0 = ', en0

      if ( int .eq. 1 .and. unc .eq. 1 ) then
         write(6,*) ' ERROR: int = 1 and unc = 1 '
         stop
      endif

C     
      I=1                                                                       
C     
C     IBLK:DNUMBER OF THE N BLOCK;NBR:DPLACE IN THIS BLOCK                      
C     

      iquit = 0

      do while ( iquit .eq. 0 )
         
         IF(IRMW.EQ.1) THEN                                                        
            READ(5,239,iostat=iquit) ETRANS(I),
     &           IVOBS(1,I),NOBS(1,I),KOBS(1,I),PROBS(1,I),IBLK(1,I),
     &           IVOBS(2,I),NOBS(2,I),KOBS(2,I),PROBS(2,I),IBLK(2,I),
     &           IAST(I),IW(I),(REF(I, NI), NI=1,reflen)

            
            if ( iquit .eq. 0 ) then

C******THE NEXT WRITE IS USEFUL FOR DEBUGGING TYPOS IN THE DATA FILE.  JTH      
C     WRITE(6,24) ETRANS(I),IVOBS(1,I),NOBS(1,I),KOBS(1,I)                      
C     A    ,PROBS(1,I),IBLK(1,I),IVOBS(2,I),NOBS(2,I),KOBS(2,I)                  
C     B    ,PROBS(2,I),IBLK(2,I),IAST(I),I,IW(I),REF(I)                          
C     THE WEIGHT (W) IS EQUAL TO THE INVERSE OF THE SQUARE OF                   
C     EXPERIMENTAL UNCERTAINTY                                                  
               IF(IW(I).EQ.1005) THEN                                                 
                  ETRANS(I)=ETRANS(I)/29979.2458D0                                
                  W(I)=3.59502D+13                                                    
               ELSEIF(IW(I).EQ.1010) THEN                                             
                  ETRANS(I)=ETRANS(I)/29979.2458D0                                    
                  W(I)=8.98755D+12                                                    
               ELSEIF(IW(I).EQ.1040) THEN                                             
                  ETRANS(I)=ETRANS(I)/29979.2458D0                                    
                  W(I)=5.61722D+11                                                    
               ELSEIF(IW(I).EQ.1004) THEN                                             
                  ETRANS(I)=ETRANS(I)/29979.2458D0                                    
                  W(I)=5.61722D+13                                                    
               ELSEIF(IW(I).EQ.1020) THEN                                             
                  ETRANS(I)=ETRANS(I)/29979.2458D0                                    
                  W (I)=2.24688D+12                                                   
               ELSEIF(IW(I).EQ.1050) THEN                                             
                  ETRANS(I)=ETRANS(I)/29979.2458D0                                    
                  W(I)=3.59502D+11                                                    
               ELSEIF(IW(I).EQ.1030) THEN                                             
                  ETRANS(I)=ETRANS(I)/29979.2458D0                                    
                  W(I)=9.98617D+11                                                    
               ELSEIF(IW(I).EQ.1008) THEN                                             
                  ETRANS(I)=ETRANS(I)/29979.2458D0                                    
                  W(I)=1.4043D+13                                                     
               ELSEIF(IW(I).EQ.1100) THEN                                             
                  ETRANS(I)=ETRANS(I)/29979.2458D0                                    
                  W(I)=8.98755D+10                                                    
               ELSEIF(IW(I).EQ.1990) THEN                                             
                  ETRANS(I)=ETRANS(I)/29979.2458D0                                    
                  W(I)=8.98755D+8                                                     
               ELSEIF(IW(I).EQ.1180) THEN                                             
                  ETRANS(I)=ETRANS(I)/29979.2458D0                                    
                  W(I)=2.77393D+10                                                    
               ELSEIF(IW(I).EQ.1130) THEN                                             
                  ETRANS(I)=ETRANS(I)/29979.2458D0                                    
                  W(I)=5.31807D+10                                                    
               ELSEIF(IW(I).EQ.1150) THEN                                             
                  ETRANS(I)=ETRANS(I)/29979.2458D0                                    
                  W(I)=3.51076D+10                                                    
               ELSEIF(IW(I).EQ.1160) THEN                                             
                  ETRANS(I)=ETRANS(I)/29979.2458D0                                    
                  W(I)=3.99446D+10                                                    
               ELSEIF(IW(I).EQ.1200) THEN                                             
                  ETRANS(I)=ETRANS(I)/29979.2458D0                                    
                  W(I)=2.24688D+10                                                    
               ELSEIF(IW(I).EQ.1300) THEN                                             
                  ETRANS(I)=ETRANS(I)/29979.2458D0                                    
                  W(I)=9.98617D+9                                                    
               ELSEIF(IW(I).EQ.1080) THEN                                             
                  ETRANS(I)=ETRANS(I)/29979.2458D0                                    
                  W(I)=1.4043D+11                                                     
               ELSEIF(IW(I).EQ.1070) THEN                                             
                  ETRANS(I)=ETRANS(I)/29979.2458D0                                    
                  W(I)=1.83419D+11                                                    
               ELSEIF(IW(I).EQ.2000) THEN                                             
                  ETRANS(I)=ETRANS(I)/29979.2458D0                                    
                  W(I)=2.24688D+8                                                     
               ELSEIF(IW(I).EQ.1000.AND.IVOBS(1,I).EQ.1) THEN                         
                  W(I)=5.61722D+11                                                    
                  ETRANS(I)=ETRANS(I)/29979.2458D0                                    
               ELSE IF(IW(I).EQ.10) THEN                                              
                  W(I)=4000000.D0                                                     
               else if (iw(i).eq.1) then
                  w(i)=25000000.d0   
               ELSE IF(IW(I).EQ.2) THEN                                               
C                  W(I)=8000000.D0
                  w(i)=8163265.d0
               else if (iw(i).eq.3) then
                  w(i)=4000000.d0
               else if (iw(i).eq.4) then
                  w(i)=2777777.d0
               else
C     W(I)=8.0D+6                                                         
C***  25,000,000=(1/0.0002)**2 FOR METHANOL IR DATA.
                  W(I)=25.0D+6                   
               ENDIF                                                                  
               IF(ETRANS(I).EQ.-900.0000)  THEN                                       
                  NADAT=I-1                                                           
                  cycle
               ENDIF                                                                  
               I=I+1                                                                  
            end if
         ELSE                                                                      
            W(I)=1.D0                                                              
C***  OVERTONE TRANSITIONS WEIGHTED AS 0.002 ACCURACY
            W(I)=25.0D+4                   
C     NEXT READ STATEMENT IS NOW CONSISTANT WITH IRMW = 1, LHX                 
            READ(5,23,iostat=iquit) ETRANS(I),IVOBS(1,I),NOBS(1,I),                      
     A           KOBS(1,I),PROBS(1,I),IBLK(1,I),IVOBS(2,I),NOBS(2,I),                    
     B           KOBS(2,I),PROBS(2,I),IBLK(2,I),IAST(I),IW(I)                            
C***  NEXT WRITE STATEMENT COMMENTED BY LHX
C     WRITE(6,226) ETRANS(I),IVOBS(1,I),NOBS(1,I),KOBS(1,I)                     
C     A,PROBS(1,I),IBLK(1,I),IVOBS(2,I),NOBS(2,I),KOBS(2,I)                      
C     B,PROBS(2,I),IBLK(2,I),IAST(I),I,IW(I)                                     
 225        FORMAT  (F12.5,1X,I3,1X,I3,1X,I3,1X,A1,1X,I3,1X,I3,1X,I3,1X,              
     1           I3,1X,A1,1X,I3,1X,I3,1X,I4)                                               
 226        FORMAT  (1X,F12.5,1X,I3,1X,I3,1X,I3,1X,A1,1X,I3,1X,I3,1X,I3,
     1           1X,I3,1X,A1,1X,I3,1X,I3,1X,I4,1X,I4)                                         
C     THE FOLLOWING CARDS TO BE TURNED IN|                                      
C     READ(5,23,END=11) ETRANS(I),IVOBS(1,I),NOBS(1,I),                      
C     A        KOBS(1,I),PROBS(1,I),IBLK(1,I),IVOBS(2,I),NOBS(2,I),              
C     B        KOBS(2,I),PROBS(2,I),IBLK(2,I),IAST(I),IW(I),REF(I)               
C     WRITE(6,26) ETRANS(I),IVOBS(1,I),NOBS(1,I),KOBS(1,I)                      
C     A ,PROBS(1,I),IBLK(1,I),IVOBS(2,I),NOBS(2,I),KOBS(2,I)                     
C     B ,PROBS(2,I),IBLK(2,I),IAST(I),I,IW(I),REF(I)                             
C     23   FORMAT  (F12.5,I2,2X,I2,2X,I3,2X,A1,1X,I2,2X,I1,2X,I2,2X,                 
C     &     I3,2X,A1,1X,I2,2X,I1,3X,I4,3X,A1)                                    
C     25   FORMAT  (F12.5,1X,I3,1X,I3,1X,I3,1X,A1,1X,I3,1X,I3,1X,I3,1X,              
C     1     I3,1X,A1,1X,I3,1X,I3,1X,I4,3X,A1)                                    

            if ( iquit .eq. 0 ) then
               IF(ETRANS(I).EQ.-900.)THEN                                             
                  NADAT=I-1                                                           
                  cycle
               ENDIF                                                                  
               I=I+1                                                                  
               cycle
            ENDIF                                                                     
         end if
      end do

      NDATA=I-1                                                              
      PRINT*,NDATA,NADAT                                                     


C     
C     NADAT:DNBREOF DATA A TYPE(SIGMA=0);NDATA:DNBRE OF TOTAL DATA              
C     
      IF(NDATA.GT.(NDTA))THEN                                                  
         WRITE(6,*) ' TOO MANY DATA...'                                         
         STOP                                                                   
      ENDIF                                                                     
C     3    FORMAT  (F12.5,1X,I3,1X,I3,1X,I3,1X,A1,1X,I3,1X,I3,1X,I3,1X,              
C     1I3,1X,A1,1X,I3,1X,I3,1X,I4)                                              

C     reflen
 239  FORMAT  (F12.5,I2,2X,I2,2X,I3,2X,A1,1X,I2,2X,I1,2X,I2,2X,                 
     &     I3,2X,A1,1X,I2,2X,I1,3X,I4,3X,40A)                                    

 23   FORMAT  (F12.5,I2,2X,I2,2X,I3,2X,A1,1X,I2,2X,I1,2X,I2,2X,                 
     &     I3,2X,A1,1X,I2,2X,I1,3X,I4,3X,A1)                                    
 25   FORMAT  (F12.5,1X,I3,1X,I3,1X,I3,1X,A1,1X,I3,1X,I3,1X,I3,1X,              
     1     I3,1X,A1,1X,I3,1X,I3,1X,I4,3X,A1)                                    
 26   FORMAT  (1X,F12.5,1X,I3,1X,I3,1X,I3,1X,A1,1X,I3,1X,I3,1X,I3,1X,           
     1     I3,1X,A1,1X,I3,1X,I3,1X,I4,1X,I4,3X,A1)                              
 24   FORMAT  (1X,F12.5,1X,I3,1X,I3,1X,I3,1X,A1,1X,I3,1X,I3,1X,I3,1X,           
     1     I3,1X,A1,1X,I3,1X,I3,1X,I4,1X,I4,3X,A1)                              
      IF(IRMW.EQ.1) NP(2)=0                                                     
C     PRINT*, 'NP(1) A LA FIN DE DATA ',NP(1),'NP(2)= ',NP(2)                   
      RETURN                                                                    
      END                                                                       
C     *************************************************************************



