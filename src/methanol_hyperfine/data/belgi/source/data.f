
      subroutine data1(cst)                                                     
C     *************************************************************************
C     This subroutine reads the experimentals values and calls two
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

      integer ii, i


      WRITE(6,*) ' This output is from the program `globalfit`'
      READ(5,*) EPSI,IBUG1,IBUG2,IBUG3,IBGTME,IBGVTD ,IBUG4                     
      WRITE(6,*) ' EPSI= ',EPSI,' IBUG1= ',IBUG1,' IBUG2= ',                    
     *     IBUG2,' IBUG3= ',IBUG3,' IBGTME= ',IBGTME,' IBGVTDV= ',              
     *     IBGVTD ,'IBUG4= ',IBUG4                                              
      READ(5,20) KTRONC                                                         
      READ(5,20) IRMW                                                           
 20   FORMAT(I3)                                                                
      WRITE(6,21)KTRONC,IRMW                                                    
 21   FORMAT(' KTRONC = ',I3,' IRMW= ',I3)                                      
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
            do  ii=1,(77)
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
      READ(5,*) ITO,NITT,IKAKC,JMAX                                             
      WRITE(6,*) ' Vibrational energy to be fitted? ITO = ', ITO                          
      WRITE(6,*)' Number of iterations: ',NITT                                  
      WRITE(6,*)' Maximum value of J: ',JMAX                                    
C                                                                               
      I=1                                                                       
C                                                                               
C     IBLK:DNUMBER OF THE N BLOCK;NBR:DPLACE IN THIS BLOCK                      
C                                                                               
 1    IF(IRMW.EQ.1) THEN                                                        
         READ(5,23,END=10) ETRANS(I),IVOBS(1,I),NOBS(1,I),                      
     A        KOBS(1,I),PROBS(1,I),IBLK(1,I),IVOBS(2,I),NOBS(2,I),              
     B        KOBS(2,I),PROBS(2,I),IBLK(2,I),IAST(I),IW(I),REF(I)               
C******THE NEXT WRITE IS USEFUL FOR DEBUGGING TYPOS IN THE DATA FILE.  JTH      
C     WRITE(6,24) ETRANS(I),IVOBS(1,I),NOBS(1,I),KOBS(1,I)                      
C    A    ,PROBS(1,I),IBLK(1,I),IVOBS(2,I),NOBS(2,I),KOBS(2,I)                  
C    B    ,PROBS(2,I),IBLK(2,I),IAST(I),I,IW(I),REF(I)                          
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
         ELSE IF(IW(I).EQ.2) THEN                                               
            W(I)=8000000.D0                                                     
         ELSE                                                                   
C           W(I)=8.0D+6                                                         
C***     25,000,000=(1/0.0002)**2 FOR METHANOL IR DATA.
            W(I)=25.0D+6                   
         ENDIF                                                                  
         IF(ETRANS(I).EQ.-900.0000)  THEN                                       
            NADAT=I-1                                                           
            GO TO 1                                                             
         ENDIF                                                                  
         I=I+1                                                                  
         GO TO 1                                                                
 10      NDATA=I-1                                                              
         PRINT*,NDATA,NADAT                                                     
      ELSE                                                                      
         W(I)=1.D0                                                              
C***     OVERTONE TRANSITIONS WEIGHTED AS 0.002 ACCURACY
            W(I)=25.0D+4                   
C      NEXT READ STATEMENT IS NOW CONSISTANT WITH IRMW = 1, LHX                 
        READ(5,23,END=11) ETRANS(I),IVOBS(1,I),NOBS(1,I),                      
     A  KOBS(1,I),PROBS(1,I),IBLK(1,I),IVOBS(2,I),NOBS(2,I),                    
     B  KOBS(2,I),PROBS(2,I),IBLK(2,I),IAST(I),IW(I)                            
C*** NEXT WRITE STATEMENT COMMENTED BY LHX
C     WRITE(6,226) ETRANS(I),IVOBS(1,I),NOBS(1,I),KOBS(1,I)                     
C    A,PROBS(1,I),IBLK(1,I),IVOBS(2,I),NOBS(2,I),KOBS(2,I)                      
C    B,PROBS(2,I),IBLK(2,I),IAST(I),I,IW(I)                                     
225   FORMAT  (F12.5,1X,I3,1X,I3,1X,I3,1X,A1,1X,I3,1X,I3,1X,I3,1X,              
     1I3,1X,A1,1X,I3,1X,I3,1X,I4)                                               
226   FORMAT  (1X,F12.5,1X,I3,1X,I3,1X,I3,1X,A1,1X,I3,1X,I3,1X,I3,1X,           
     1I3,1X,A1,1X,I3,1X,I3,1X,I4,1X,I4)                                         
C     THE FOLLOWING CARDS TO BE TURNED IN|                                      
C        READ(5,23,END=11) ETRANS(I),IVOBS(1,I),NOBS(1,I),                      
C    A        KOBS(1,I),PROBS(1,I),IBLK(1,I),IVOBS(2,I),NOBS(2,I),              
C    B        KOBS(2,I),PROBS(2,I),IBLK(2,I),IAST(I),IW(I),REF(I)               
C     WRITE(6,26) ETRANS(I),IVOBS(1,I),NOBS(1,I),KOBS(1,I)                      
C    A ,PROBS(1,I),IBLK(1,I),IVOBS(2,I),NOBS(2,I),KOBS(2,I)                     
C    B ,PROBS(2,I),IBLK(2,I),IAST(I),I,IW(I),REF(I)                             
C23   FORMAT  (F12.5,I2,2X,I2,2X,I3,2X,A1,1X,I2,2X,I1,2X,I2,2X,                 
C    &     I3,2X,A1,1X,I2,2X,I1,3X,I4,3X,A1)                                    
C25   FORMAT  (F12.5,1X,I3,1X,I3,1X,I3,1X,A1,1X,I3,1X,I3,1X,I3,1X,              
C    1     I3,1X,A1,1X,I3,1X,I3,1X,I4,3X,A1)                                    
         IF(ETRANS(I).EQ.-900.)THEN                                             
            NADAT=I-1                                                           
            GO TO 1                                                             
         ENDIF                                                                  
         I=I+1                                                                  
         GO TO 1                                                                
 11      NDATA=I-1                                                              
         PRINT*,NDATA,NADAT                                                     
      ENDIF                                                                     
C                                                                               
C     NADAT:DNBREOF DATA A TYPE(SIGMA=0);NDATA:DNBRE OF TOTAL DATA              
C                                                                               
      IF(NDATA.GT.(NDTA))THEN                                                  
         WRITE(6,*) ' TOO MANY DATA...'                                         
         STOP                                                                   
      ENDIF                                                                     
C3    FORMAT  (F12.5,1X,I3,1X,I3,1X,I3,1X,A1,1X,I3,1X,I3,1X,I3,1X,              
C     1I3,1X,A1,1X,I3,1X,I3,1X,I4)                                              
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


C     *************************************************************************
      subroutine data11(lo,cst)                                                 
C     *************************************************************************
C     Read the constants for each state and put them in the CST vector.
C     *************************************************************************
      implicit none
      save

      include 'common_dim.f'
                                                                      
      CHARACTER*4 ANSWER     
      integer ii, iii, lo

      include 'common_symb.f'
      include 'common_chara.f'
      include 'common_dat.f'
      include 'common_cst.f'
      include 'common_cste.f'

*     MODIFICATION DENIS GIROU (CIRCE)                                          
*     ABC = 'OA    B     C     DJ    DJK   DK    ODELN ODELK HJ    HJK          
*     1 HKJ   HK    OHJ   OHJK  OHK   FV    OFV   GV    ALV   AK1   AK2         
*     2 AK3   AK4   AK5   AK6   AK7   DAB   ODAB  DABJ  DABK  OLV   C1          
*     3 C2    C3    HLKJ  AK5J  C2J   GVJ   AK2J  C1J   ANV   BK2   AMV         
*     4 DAC   DBC   C4    ODELTADELTA AK1J  AK5K  AK2K  AK1K  AK3J  AK3K        
*     5 F     RHO   V3    V6    V9    RHORHO'                                   
      ABC(1:210)=                                                               
     1      'OA    B     C     DJ    DJK   DK    ODELN ODELK HJ    HJK          
     2 HKJ   HK    OHJ   OHJK  OHK   FV    OFV   GV    ALV   AK1   AK2          
     3 AK3   AK4   AK5   AK6   AK7   DAB   ODAB  DABJ  DABK  OLV   C1           
     4 C2    C3    HLKJ  '                                                      
      ABC(211:462)=     'AK5J  C2J   GVJ   AK2J  C1J   ANV   BK2   AMV          
     4 DAC   DBC   C4    ODELTADELTA AK1J  AK5K  AK2K  AK1K  AK3J  AK3K         
     5 BK1   AK7J  C11   ODAB6 AK6K  AK7K  AK66  AK6J  AK3B  AK4B  C2K
     6 C1K   C4J   C4K   DELTABODELTBC12   F     RHO   V3    V6    V9
     7 RHORHO'                              
C     INITIALISATION                                                            
C                                                                               
      OA=0.D0                                                                   
      B=0.D0                                                                    
      C=0.D0                                                                    
      DJ=0.D0                                                                   
      DJK=0.D0                                                                  
      DK=0.D0                                                                   
      ODELN=0.D0                                                                
      ODELK=0.D0                                                                
      HJ=0.D0                                                                   
      HJK=0.D0                                                                  
      HKJ=0.D0                                                                  
      HK=0.D0                                                                   
      OHJ=0.D0                                                                  
      OHJK=0.D0                                                                 
      OHK=0.D0                                                                  
      FV=0.D0                                                                   
      OFV=0.D0                                                                  
      GV=0.D0                                                                   
      ALV=0.D0                                                                  
      AK1=0.D0                                                                  
      AK2=0.D0                                                                  
      AK3=0.D0                                                                  
      AK4=0.D0                                                                  
      AK5=0.D0                                                                  
      AK6=0.D0                                                                  
      AK7=0.D0                                                                  
      DAB=0.D0                                                                  
      DBC=0.D0                                                                  
      C4=0.D0                                                                   
      ODELTA=0.D0                                                               
      AK1J=0.D0                                                                 
      AK5K=0.D0                                                                 
      AK2K=0.D0                                                                 
      AK1K=0.D0                                                                 
      AK3J=0.D0                                                                 
      AK3K=0.D0                                                                 
      DELTA=0.D0                                                                
      ODAB=0.D0                                                                 
      DABJ=0.D0                                                                 
      DABK=0.D0                                                                 
      OLV=0.D0                                                                  
      C1=0.D0                                                                   
      C2=0.D0                                                                   
      C3=0.D0                                                                   
      HLKJ=0.D0                                                                 
      AK5J=0.D0                                                                 
      C2J=0.D0                                                                  
      GVJ=0.D0                                                                  
      AK2J=0.D0                                                                 
      C1J=0.D0                                                                  
      ANV=0.D0                                                                  
      BK2=0.D0                                                                  
      AMV=0.D0                                                                  
      DAC=0.D0                                                                  
      BK1=0.D0
      AK7J=0.D0
      C11=0.D0
      ODAB6=0.D0
      AK6K=0.D0
      AK7K=0.D0
      AK66=0.D0
      AK6J=0.D0
      AK3B=0.D0
      AK4B=0.D0
      C2K=0.D0
      C1K=0.D0
      C4J=0.D0
      C4K=0.D0
      DELTAB=0.D0
      ODELTB=0.D0
      C12=0.D0
      F=0.D0                                                                    
      RHO=0.D0                                                                  
      V3=0.D0                                                                   
      V6=0.D0                                                                   
      V9=0.D0                                                                   
C     RHOJ=0.D0                                                                 
C     RHOK=0.D0                                                                 
      RHORHO=0.D0                                                               
C                                                                               
C     INITIALISATION                                                            
C                                                                               
      DO 3 II=1,(77)                                                            
         CST(LO,II)=0.D0                                                        
 3    CONTINUE                                                                  
C     ASK THE USER FOR THE CONSTANT VALUES                                      
      WRITE(6,132)                                                              
 132  FORMAT('   CONSTANTS IN CM-1 OR IN MHZ?  ',':DA')                         
      READ(5,'(A4)') ANSWER                                                     
      WRITE(6,'(A4)') ANSWER                                                    
      CONV=1.D0                                                                 
      IF(ANSWER.EQ.'MHZ') CONV=29979.2458D0                                    
C     READ FROM NAMELIST                                                        
      READ(5,CSTE)                                                              
      WRITE(6,CSTE)                                                             
C     PUT THE CONSTANT VALUES IN CST VECTOR                                     
      WRITE(35) OA,B,C,DJ,DJK,DK,ODELN,ODELK,HJ,HJK,HKJ,HK,                     
     1     OHJ,OHJK,OHK,FV,OFV,GV,ALV,AK1,AK2,AK3,AK4,AK5,AK6,AK7,              
     2     DAB,ODAB,DABJ,DABK,OLV,C1,C2,C3,HLKJ,AK5J,C2J,GVJ,AK2J,C1J           
     3     ,ANV,BK2,AMV,DAC,DBC,C4,ODELTA,DELTA,AK1J,AK5K,AK2K,AK1K,AK3J        
     4     ,AK3K,BK1,AK7J,C11,ODAB6,AK6K,AK7K,AK66,AK6J,AK3B,AK4B               
     5     ,C2K,C1K,C4J,C4K,DELTAB,ODELTB,C12                              
     6     ,F,RHO,V3,V6,V9,RHORHO                                          
      REWIND 35                                                                 
      READ(35) (CST(LO,II),II=1,(77))                                           
      REWIND 35                                                                 
      DO 18 III=1,71                                                            
         CST(LO,III)=CST(LO,III)/CONV                                           
 18   CONTINUE                                                                  
C THE NEXT FOUR LINES ARE FOR CHECKING PRECISION
      DO 19 III=1,(77)
      WRITE(6,1000) CST(LO,III)
1000  FORMAT(D30.20)
 19   CONTINUE
C     PRINT*,CST(LO,38)                                                         
C     PRINT*,RHO,CST(LO,1)                                                      
      RETURN                                                                    
      END                                                                       
C     *************************************************************************


