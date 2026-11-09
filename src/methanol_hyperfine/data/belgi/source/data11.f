

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

*     constant table (names)
      ABC(   1:  36)= 'OA    B     C     DJ    DJK   DK    '
      ABC(  37:  72)= 'ODELN ODELK HJ    HJK   HKJ   HK    '
      ABC(  73: 108)= 'OHJ   OHJK  OHK   FV    OFV   GV    '
      ABC( 109: 144)= 'ALV   AK1   AK2   AK3   AK4   AK5   '
      ABC( 145: 180)= 'AK6   AK7   DAB   ODAB  DABJ  DABK  '
      ABC( 181: 216)= 'OLV   C1    C2    C3    DUMMY AK5J  '
      ABC( 217: 252)= 'C2J   GVJ   AK2J  C1J   ANV   BK2   '
      ABC( 253: 288)= 'AMV   DAC   DBC   C4    ODELTADELTA '
      ABC( 289: 324)= 'AK1J  AK5K  AK2K  AK1K  AK3J  AK3K  '
      ABC( 325: 360)= 'BK1   AK7J  C11   ODAB6 AK6K  AK7K  '
      ABC( 361: 396)= 'AK66  AK6J  AK3B  AK4B  C2K   C1K   '
      ABC( 397: 432)= 'C4J   C4K   DELTABODELTBC12   FPARA '
      ABC( 433: 468)= 'RHO   V3    V6    V9    RHORHOLJ    '
      ABC( 469: 504)= 'LK    LJKKK LJJKK LJJJK OLKKK OLJKK '
      ABC( 505: 540)= 'OLJJK OLJJJ HABK  HABJ  HABJK GJ6   '
      ABC( 541: 576)= 'GJ4K2 GJ2K4 GK6   FK6   AK6KK AK6JK '
      ABC( 577: 612)= 'AK6JJ AK7KK AK7JK AK7JJ AK3BK AMUA0 '
      ABC( 613: 648)= 'AMUA3 AMUB0 AMUB3 AMUC3 AGJ6  AGJ4K2'
      ABC( 649: 684)= 'AGJ2K4AGK6  FJ6   FJ4K2 FJ2K4 AK7BC '
      ABC( 685: 720)= 'V9J   V9K   ODABK C6J4  C6J2K2C6K4  '
      ABC( 721: 756)= 'G4J4  G4J2K2G4K4  GAJ4  GAJ2K2GAK4  '
      ABC( 757: 792)= 'AK3BJ AK4BJ AK4BK ODABJ ODAC  ODBC  '
      ABC( 793: 828)= 'DELTAJDELTAKDAGJ  DAGK  CABK  CABJ  '
      ABC( 829: 864)= 'DUMMY DUMMY DG4J  DG4K  DACK  DBCJ  '
      ABC( 865: 900)= 'DACJ  V12   V12J  AK4C  AK3C  DAC6  '
      ABC( 901: 936)= 'DBC6  GB    AK8   AK9   AK8J  AK8K  '
      ABC( 937: 972)= 'AK9J  AK9K  DELTJKDELTKKODAB9 C9BC  '
      ABC( 973:1008)= 'C6BCJ C6BCK AG3BCJAG3BCKG4BCJ G4BCK '
      ABC(1009:1044)= 'AK4BC AK3BC AK7AB BC4   CORRECAK76  '
      ABC(1045:1080)= 'OBC4  AK3BK6AK4BK4AK3CK4AK4D  AK3D  '
      ABC(1081:1116)= 'AK4CK AK3CK AK3BK4AK4BK6AK4CK4DACJJ '
      ABC(1117:1152)= 'DACJK DACKK DBCK  DBCJK ODABKKC2JK  '
      ABC(1153:1188)= 'C1JK  C4JK  C6BCJKG3BCJKG4BCJKC2KK  '
      ABC(1189:1224)= 'C1KK  C4KK  C6BCKKG3BCKKG4BCKKDG3J  '
      ABC(1225:1260)= 'DG3K                                '
      ABC(1261:1298)= '                                    '
      ABC(1297:1332)= '                                    '


C     INITIALISATION                                                         
C     
      C2JK = 0.d0
      C2KK = 0.d0

      C1JK = 0.d0
      C1KK = 0.d0

      C2JK = 0.d0
      C2KK = 0.d0

      C6BCJK = 0.d0
      C6BCKK = 0.d0
      
      G3BCJK = 0.d0
      G3BCKK = 0.d0
      
      G4BCJK = 0.d0
      G4BCKK = 0.d0

      DG3J = 0.d0
      DG3K = 0.d0

      DBCK  = 0.d0
      DBCJK  = 0.d0

      DACJJ = 0.d0
      DACJK = 0.d0
      DACKK = 0.d0


      AK4BK6 = 0.d0
      AK4CK4 = 0.d0

      AK4D = 0.d0
      AK3D = 0.d0
      AK4CK = 0.d0
      AK3CK = 0.d0

      AK3BK6 = 0.d0
      AK4BK4 = 0.d0
      AK3BK4 = 0.d0
      AK3CK4 = 0.d0

      AK7AB = 0.d0
      AK76 = 0.d0

      AK4BC = 0.d0
      AK3BC = 0.d0

      dummy1 = 0.d0

      C6BCJ = 0.d0
      C6BCK = 0.d0

      AG3BCJ = 0.d0
      AG3BCK = 0.d0

      G4BCJ = 0.d0
      G4BCK = 0.d0

      ODAB9 = 0.d0
      C9BC = 0.d0

      DELTJK = 0.d0
      DELTKK = 0.d0

      AK8 = 0.d0
      AK8J = 0.d0
      AK8K = 0.d0

      AK9 = 0.d0
      AK9J = 0.d0
      AK9K = 0.d0

      GB = 0.d0

      DAC6 = 0.d0
      DBC6 = 0.d0

      AK4C = 0.d0
      AK3C = 0.d0

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
      FPARA=0.D0                                                                    
      RHO=0.D0                                                                  
      V3=0.D0                                                                   
      V6=0.D0                                                                   
      V9=0.D0                                                                   
      
      V12 = 0.d0
      V12J = 0.d0

C     RHOJ=0.D0                                                                 
C     RHOK=0.D0                                                                 
      RHORHO=0.D0                                                               
C     
      LJ=0.0D0
      LK=0.0D0
      LJKKK=0.0D0
      LJJKK = 0.0d0
      LJJJK = 0.0d0
      OLKKK = 0.0d0
      OLJKK = 0.0d0
      OLJJK = 0.0d0
      OLJJJ = 0.0d0
      HABK = 0.0d0
      HABJ = 0.0d0
      HABJK = 0.0d0

      GJ6 = 0.0d0
      GJ4K2 = 0.0d0
      GJ2K4 = 0.0d0
      GK6 = 0.0d0

      FJ6 = 0.0d0
      FJ4K2 = 0.0d0
      FJ2K4 = 0.0d0
      FK6 = 0.0d0

      AK6KK = 0.0d0
      AK6JK = 0.0d0
      AK6JJ = 0.0d0

      BC4 = 0.d0
      OBC4 = 0.d0
 
      AK7BC = 0.d0
      AK7KK = 0.d0
      AK7JK = 0.d0
      AK7JJ = 0.d0

      CORREC = 0.d0

      AK3BK = 0.d0
      AK3BJ = 0.d0
      AK4BJ = 0.d0
      AK4BK = 0.d0

      AMUA0 = 0.d0
      AMUA3 = 0.d0
      AMUB0 = 0.d0
      AMUA3 = 0.d0
      AMUC3 = 0.d0

      AGJ6 = 0.d0
      AGJ4K2 = 0.d0
      AGJ2K4 = 0.d0
      AGK6 = 0.d0

      V9J = 0.d0
      V9K = 0.d0

      ODABJ = 0.d0
      ODABK = 0.d0  
      ODABKK = 0.d0
      
      ODAC = 0.d0
      ODBC = 0.d0

      DACK = 0.d0
      DACJ = 0.d0
      DBCJ = 0.d0

      C6J4 = 0.d0
      C6J2K2 = 0.d0
      C6K4 = 0.d0

      G4J4 = 0.d0
      G4J2K2 = 0.d0
      G4K4 = 0.d0

      GAJ4 = 0.d0
      GAJ2K2 = 0.d0
      GAK4 = 0.d0

      DELTAJ = 0.d0
      DELTAK = 0.d0
      
      DAGJ = 0.d0
      DAGK = 0.d0

      DG4J = 0.d0
      DG4K = 0.d0

      DABK = 0.d0

      CABK = 0.d0
      CABJ = 0.d0
      CBCK = 0.d0
      CBCJ = 0.d0

C     INITIALISATION                                                            
C     
      DO 3 II=1,cstm                                                            
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
     2     DAB,ODAB,DABJ,DABK,OLV,C1,C2,C3,dummy1,AK5J,C2J,GVJ,AK2J,C1J           
     3     ,ANV,BK2,AMV,DAC,DBC,C4,ODELTA,DELTA,AK1J,AK5K,AK2K,AK1K,AK3J        
     4     ,AK3K,BK1,AK7J,C11,ODAB6,AK6K,AK7K,AK66,AK6J,AK3B,AK4B               
     5     ,C2K,C1K,C4J,C4K,DELTAB,ODELTB,C12                              
     6     ,FPARA,RHO,V3,V6,V9,RHORHO,LJ,LK,LJKKK,LJJKK,LJJJK,
     7     OLKKK,OLJKK, OLJJK, OLJJJ,
     8     HABK, HABJ, HABJK, GJ6, GJ4K2, GJ2K4, GK6, FK6,
     9     AK6KK, AK6JK, AK6JJ,
     *     AK7KK, AK7JK, AK7JJ, AK3BK,
     1     AMUA0, AMUA3, AMUB0, AMUB3, AMUC3,
     2     AGJ6, AGJ4K2, AGJ2K4, AGK6,
     3     FJ6, FJ4K2, FJ2K4,
     4     AK7BC,
     5     V9J, V9K,
     6     ODABK,
     7     C6J4, C6J2K2, C6K4,
     8     G4J4, G4J2K2, G4K4,
     9     GAJ4, GAJ2K2, GAK4,
     *     AK3BJ, AK4BJ, AK4BK,
     1     ODABJ, ODAC, ODBC, DELTAJ, DELTAK, DAGJ, DAGK,
     2     CABK, CABJ, CBCK, CBCJ,
     3     DG4J, DG4K, DACK, DBCJ, DACJ, V12, V12J,
     4     AK4C, AK3C, DAC6, DBC6, GB, AK8, AK9,
     5     AK8J, AK8K, AK9J, AK9K, DELTJK, DELTKK, ODAB9,
     6     C9BC, C6BCJ, C6BCK, AG3BCJ, AG3BCK, G4BCJ, G4BCK,
     7     AK4BC, AK3BC, AK7AB, BC4, CORREC, AK76, OBC4,
     8     AK3BK6, AK4BK4, AK3CK4,
     9     AK4D, AK3D, AK4CK, AK3CK, AK3BK4,
     *     AK4BK6, AK4CK4, DACJJ, DACJK, DACKK, DBCK,
     1     DBCJK, ODABKK, C2JK, C1JK, C4JK, C6BCJK,
     2     G3BCJK, G4BCJK, C2KK, C1KK,C4KK,C6BCKK,
     3     G3BCKK, G4BCKK, DG3J, DG3K

      REWIND 35                                                                 
      READ(35) (CST(LO,II),II=1,cstm)                                           
      REWIND 35                     

      print *, 'CONV is ', CONV
                                            
      DO 18 III=1,71                                                            
         CST(LO,III)=CST(LO,III)/CONV                                           
 18   CONTINUE                                                                  


C     jmf be careful about conversion here
      do III = 78, cstm
         CST(LO,III) = CST(LO,III)/CONV
      end do

C     THE NEXT FOUR LINES ARE FOR CHECKING PRECISION
      DO 19 III=1,cstm
         WRITE(6,1000) CST(LO,III)
 1000    FORMAT(D30.20)
 19   CONTINUE
C     PRINT*,CST(LO,38)                                                         
C     PRINT*,RHO,CST(LO,1)                                                      
      RETURN                                                                    
      END                                                                       
C     *************************************************************************


