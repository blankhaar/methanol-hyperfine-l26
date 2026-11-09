                                                                 
      double precision dummy1
      double precision OA
      double precision B     
      double precision C  
      double precision DJ    
      double precision DJK                                                                  
      double precision DK                                                                   
      double precision ODELN                                                                
      double precision ODELK                                                                
      double precision HJ                                                                   
      double precision HJK                                                                  
      double precision HKJ                                                                  
      double precision HK                                                                   
      double precision OHJ                                                                  
      double precision OHJK                                                                 
      double precision OHK                                                                  
      double precision FV                                                                   
      double precision OFV                                                                  
      double precision GV                                                                   
      double precision ALV                                                                  
      double precision AK1                                                                  
      double precision AK2                                                                  
      double precision AK3                                                                  
      double precision AK4                                                                  
      double precision AK5                                                                  
      double precision AK6                                                                  
      double precision AK7                                                                  
      double precision DAB  
      double precision DBC  
      double precision C4   
      double precision ODELTA
      double precision AK1J 
      double precision AK5K 
      double precision AK2K 
      double precision AK1K 
      double precision AK3J 
      double precision AK3K 
      double precision DELTA
      double precision ODAB 
      double precision DABJ 
      double precision DABK 
      double precision OLV  
      double precision C1   
      double precision C2   
      double precision C3   
      double precision AK5J 
      double precision C2J  
      double precision GVJ  
      double precision AK2J 
      double precision C1J  
      double precision ANV  
      double precision BK2  
      double precision AMV  
      double precision DAC  
      double precision BK1
      double precision AK7J
      double precision C11
      double precision ODAB6
      double precision AK6K
      double precision AK7K
      double precision AK66
      double precision AK6J
      double precision AK3B
      double precision AK4B
      double precision C2K
      double precision C1K
      double precision C4J
      double precision C4K
      double precision DELTAB
      double precision ODELTB
      double precision C12
      double precision FPARA 
      double precision RHO 
      double precision V3   
      double precision V6   
      double precision V9   
C     RHOJ 
C     RHOK 
      double precision RHORHO       
C                                                                               
      double precision LJ
      double precision LK
      double precision LJKKK
      double precision LJJKK
      double precision LJJJK
      double precision OLKKK, OLJKK, OLJJK, OLJJJ
      double precision HABK, HABJ, HABJK
      double precision GJ6, GJ4K2, GJ2K4, GK6
      double precision FK6
      double precision AK6KK, AK6JK, AK6JJ
      double precision AK7KK, AK7JK, AK7JJ
      double precision AK3BK, AK3BJ
      double precision AK4BK, AK4BJ
      double precision AMUA0, AMUA3, AMUB0, AMUB3, AMUC3
      double precision AGJ6, AGJ4K2, AGJ2K4, AGK6
      double precision FJ6, FJ4K2, FJ2K4
      double precision AK7BC
      double precision V9J, V9K
      double precision ODABJ, ODABK, ODBC, ODAC, DELTAJ, DELTAK
      double precision C6J4, C6J2K2, C6K4
      double precision G4J4, G4J2K2, G4K4
      double precision GAJ4, GAJ2K2, GAK4
      double precision DAGJ, DAGK
      double precision CABK, CABJ, CBCK, CBCJ
      double precision DG4J, DG4K, DACK, DBCJ, DACJ
      double precision V12, V12J, AK4C, AK3C
      double precision DAC6, DBC6, GB
      double precision AK8, AK9, AK8J, AK8K, AK9J, AK9K
      double precision DELTJK, DELTKK
      double precision ODAB9, C9BC
      double precision C6BCJ, C6BCK, AG3BCJ, AG3BCK, G4BCJ, G4BCK
      double precision AK4BC, AK3BC, AK7AB, BC4, CORREC
      double precision AK76, OBC4
      double precision AK3BK6, AK4BK4, AK3CK4
      double precision AK4D, AK3D, AK4CK, AK3CK, AK3BK4
      double precision AK4BK6, AK4CK4
      double precision DACJJ, DACJK, DACKK, DBCK, DBCJK, ODABKK
      double precision C2JK, C1JK, C4JK, C6BCJK, G3BCJK, G4BCJK
      double precision C2KK, C1KK, C4KK, C6BCKK, G3BCKK, G4BCKK
      double precision DG3J, DG3K
      

      NAMELIST/CSTE/OA,B,C,DJ,DJK,DK,ODELN,ODELK,HJ,HJK                         
     1     ,HKJ,HK,OHJ,OHJK,OHK,FV,OFV,GV,ALV,AK1,AK2,AK3,                      
     2     AK4,AK5,AK6,AK7,DAB,ODAB,DABJ,DABK,OLV,C1,C2,C3,dummy1,AK5J,           
     &     C2J,                                                                 
     3     GVJ,AK2J,C1J,ANV,BK2,AMV,DAC,DBC,C4,ODELTA,DELTA,AK1J,AK5K,          
     4     AK2K,AK1K,AK3J,AK3K,BK1,AK7J,C11,ODAB6,AK6K,AK7K,AK66,AK6J,          
     5     AK3B,AK4B,C2K,C1K,C4J,C4K,DELTAB,ODELTB,C12,               
     6     FPARA,RHO,V3,V6,V9,RHORHO,
     7     LJ,LK,LJKKK,LJJKK,LJJJK,OLKKK,
     8     OLJKK, OLJJK, OLJJJ, HABK, HABJ, HABJK, 
     9     GJ6, GJ4K2, GJ2K4, GK6,
     *     FK6,
     1     AK6KK, AK6JK, AK6JJ,
     2     AK7KK, AK7JK, AK7JJ,
     3     AK3BK,
     4     AMUA0, AMUA3, AMUB0, AMUB3, AMUC3,
     5     AGJ6, AGJ4K2, AGJ2K4, AGK6,
     6     FJ6, FJ4K2, FJ2K4,
     7     AK7BC,
     8     V9J, V9K,
     9     ODABK,
     *     C6J4, C6J2K2, C6K4,
     1     G4J4, G4J2K2, G4K4,
     2     GAJ4, GAJ2K2, GAK4,
     3     AK3BJ, AK4BJ, AK4BK,
     4     ODABJ, ODAC, ODBC, DELTAJ, DELTAK, DAGJ, DAGK, 
     5     CABK, CABJ, CBCK, CBCJ,
     6     DG4J, DG4K, DACK, DBCJ, DACJ, V12, V12J, AK4C, AK3C,
     7     DAC6, DBC6, GB, AK8, AK9, AK8J, AK8K, AK9J, AK9K,
     8     DELTJK, DELTKK, ODAB9, C9BC,
     9     C6BCJ, C6BCK, AG3BCJ, AG3BCK, G4BCJ, G4BCK, AK4BC, AK3BC,
     *     AK7AB, BC4, CORREC, AK76, OBC4, AK3BK6, AK4BK4, AK3CK4,
     1     AK4D, AK3D, AK4CK, AK3CK, AK3BK4, AK4BK6, AK4CK4,
     2     DACJJ, DACJK, DACKK, DBCK, DBCJK,ODABKK, C2JK, 
     3     C1JK, C4JK, C6BCJK, G3BCJK, G4BCJK, C2KK, C1KK,
     4     C4KK, C6BCKK, G3BCKK, G4BCKK, DG3J, DG3K
      
