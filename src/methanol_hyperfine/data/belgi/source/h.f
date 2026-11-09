C**********************************************************                     
      FUNCTION H0 (K,IVTOR,KP,IVTORP,CST,LO,ETOR,A)                             
      implicit none
      save
C***********************************************************                    
C     SETS UP THE DIAGONAL ELEMENTS OF THE HAMILTONIAN                          
      

      include 'common_dim.f'
      include 'common_dat.f'
      include 'common_chara.f'
      include 'common_rotor.f'
      include 'common_symb.f'
      include 'common_quant.f'
      include 'common_cst.f'

      integer ivtor, kp, ivtorp, lo, ik, kl, ii, k
      double precision h0, etor
      double precision oa, dj, djk, dk, odeln, odelk, hj, hjk, hkj, hk
      double precision ohj, ohjk, ohk, fv, ofv, gv, alv, ak1, ak2, ak3
      double precision ak4, ak5, ak6, ak7, dab, odab, olv, c1, c2, c3
      double precision ak5j, c2j, gvj, ak2j, c1j, anv, bk2, amv
      double precision ak1j, ak5k, ak2k, ak1k, ak3j, bk1, ak7j, ak6k
      double precision ak7k, ak66, ak6j, ak3b, ak4b, rho, rhorho, ak3k
      double precision a1, a2, a3, a4, a4b, a5, a5b, a6, a7, a8, a8b
      double precision a9b, a14b, a15b, a9, a10, a11, a12, a13, a14
      double precision a15, a16, a17, coefa, cofaaa, coefaa, coefa4
      double precision coefklp3, coefklm3, av91, av92
      double precision coefa5, coefa6, rhoeff, rhoksg, rhokk
      double precision rhokk2, rhokm1, rhokm2, b, c, a
      double precision coefklp4, coefklm4, av121, av122, bc4

      double precision lj,lk,ljkkk,ljjkk,ljjjk
      double precision gj6, gj4k2, gj2k4, gk6
      double precision fj6, fj4k2, fj2k4, fk6
      double precision ak6kk, ak6jk, ak6jj
      double precision ak7kk, ak7jk, ak7jj
      double precision ak3bj, ak3bk, ak4bj, ak4bk

      double precision agj6, agj4k2, agj2k4, agk6
      double precision v9j, v9k
      double precision odabk, odabj
      
      double precision c6j4, c6j2k2, c6k4
      double precision g4j4, g4j2k2, g4k4
      double precision gaj4, gaj2k2, gak4

      double precision ak3bk6, ak4bk4, ak3ck4, ak3bk4
      double precision ak4bk6, ak4ck4
      
      double precision v12j, a4c, a3c, ak4c, ak3c
      double precision ak8, ak9, aak8, aak9, ak8j, ak8k
      double precision ak9k, ak9j, ak76, a8c, a9c, a14c, a15c
      
      double precision AK4D, AK3D, AK4CK, AK3CK, A4D, A3D

      DIMENSION A(9,NTORMX,KDMX),ETOR(KDMX,9)                         

      XN = N * 1.0d0
      XNN = (N*(N+1)) * 1.0d0
      XK = K * 1.0d0
      

      AK4BK6 = CST(LO, 184)
      AK4CK4 = CST(LO, 185)

      AK4D = CST(LO, 179)
      AK3D = CST(LO, 180)
      AK4CK = CST(LO, 181)
      AK3CK = CST(LO, 182)

      AK3BK6 = CST(LO, 176)
      AK4BK4 = CST(LO, 177)
      AK3CK4 = CST(LO, 178)
      AK3BK4 = CST(LO, 183)

      AK8 = CST(LO, 153)
      AK8J = CST(LO, 155)
      AK8K = CST(LO, 156)

      AK9 = CST(LO, 154)
      AK9J = CST(LO, 157)
      AK9K = CST(LO, 158)

      BC4 = CST(LO, 172)

      AK76 = CST(LO, 174)

      OA=CST(LO,1)                                                              
      B=CST(LO,2)                                                               
      C=CST(LO,3)                                                               
      DJ=CST(LO,4)                                                              
      DJK=CST(LO,5)                                                             
      DK=CST(LO,6)                                                              
      ODELN=CST(LO,7)                                                           
      ODELK=CST(LO,8)                                                           
      HJ=CST(LO,9)                                                              
      HJK=CST(LO,10)                                                            
      HKJ=CST(LO,11)                                                            
      HK=CST(LO,12)                                                             
      OHJ=CST(LO,13)                                                            
      OHJK=CST(LO,14)                                                           
      OHK=CST(LO,15)                                                            
      FV=CST(LO,16)                                                             
      OFV=CST(LO,17)                                                            
      GV=CST(LO,18)                                                             
      ALV=CST(LO,19)                                                            
      AK1=CST(LO,20)                                                            
      AK2=CST(LO,21)                                                            
      AK3=CST(LO,22)                                                            
      AK4=CST(LO,23)                                                            
      AK5=CST(LO,24)                                                            
      AK6=CST(LO,25)                                                            
      AK7=CST(LO,26)                                                            
      DAB=CST(LO,27)                                                            
      ODAB=CST(LO,28)                                                           
      OLV=CST(LO,31)                                                            
      C1=CST(LO,32)                                                             
      C2=CST(LO,33)                                                             
      C3=CST(LO,34)                                                             

      AK5J=CST(LO,36)                                                           
      C2J=CST(LO,37)                                                            
      GVJ=CST(LO,38)                                                            
      AK2J=CST(LO,39)                                                           
      C1J=CST(LO,40)                                                            
      ANV=CST(LO,41)                                                            
      BK2=CST(LO,42)                                                            
      AMV=CST(LO,43)                                                            
      AK1J=CST(LO,49)                                                           
      AK5K=CST(LO,50)                                                           
      AK2K=CST(LO,51)                                                           
      AK1K=CST(LO,52)                                                           
      AK3J=CST(LO,53)                                                           
      AK3K=CST(LO,54)                                                           
      BK1=CST(LO,55)
      AK7J=CST(LO,56)
      AK6K=CST(LO,59)
      AK7K=CST(LO,60)
      AK66=CST(LO,61)
      AK6J=CST(LO,62)
      AK3B=CST(LO,63)
      AK4B=CST(LO,64)
      RHO=CST(LO,73)                                                            
      RHORHO=CST(LO,77)

      AK4C = CST(LO, 148)
      AK3C = CST(LO, 149)

      LJ = CST(LO, 78)
      LK = CST(LO, 79)
      LJKKK = CST(LO,80)
      LJJKK = CST(LO,81)
      LJJJK = CST(LO,82)

      GJ6 = CST(LO, 90)
      GJ4K2 = CST(LO, 91)
      GJ2K4 = CST(LO, 92)
      GK6 = CST(LO, 93)

      FK6 = CST(LO, 94)
      FJ6 = CST(LO, 111)
      FJ4K2 = CST(LO, 112)
      FJ2K4 = CST(LO, 113)
      
      AK6KK = CST(LO, 95)
      AK6JK = CST(LO, 96)
      AK6JJ = CST(LO, 97)

      AK7KK = CST(LO, 98)
      AK7JK = CST(LO, 99)
      AK7JJ = CST(LO, 100)

      AK3BK = CST(LO, 101)
      AK3BJ = CST(LO, 127)
      AK4BJ = CST(LO, 128)
      AK4BK = CST(LO, 129)

      AGJ6 = CST(LO, 107)
      AGJ4K2 = CST(LO, 108)
      AGJ2K4 = CST(LO, 109)
      AGK6 = CST(LO, 110)
      
      V9J = CST(LO, 115)
      V9K = CST(LO, 116)

      ODABK = CST(Lo, 117)
      ODABJ = CST(lo, 130)

      C6J4 = CST(LO, 118)
      C6J2K2 = CST(LO, 119)
      C6K4 = CST(LO, 120)

      G4J4 = CST(LO, 121)
      G4J2K2 = CST(LO, 122)
      G4K4 = CST(LO, 123)

      GAJ4 = CST(LO, 124)
      GAJ2K2 = CST(LO, 125)
      GAK4 = CST(LO, 126)

      V12J = CST(LO, 147)

C     RHOJ=CST(LO,60)                                                           
C     RHOK=CST(LO,61)                                                           
      IK=K+N   +1                                                               
C     INITIALISATION                                                            
      A1=0.d0                                                                   
      A2=0.d0                                                                   
      A3=0.d0                                                                   
      A4=0.d0                                                                   
      A4B=0.d0                                                                  
      A5=0.d0                                                                   
      A5B=0.d0                                                                  
      A6=0.d0                                                                   
      A7=0.d0                                                                   
      A7=0.d0                                                                   
      A8=0.d0                                                                   
      A8B=0.d0
      A9B=0.d0
      A8C=0.d0
      A9C=0.d0
      A14B=0.d0
      A15B=0.d0
      A14C=0.d0
      A15C=0.d0
      A9=0.d0                                                                   
      A10=0.d0                                                                  
      A11=0.d0                                                                  
      A12=0.d0                                                                  
      A13=0.d0                                                                  
      A14=0.d0                                                                  
      A15=0.d0                                                                  
      A16=0.d0                                                                  
      A17=0.d0                                                                  

      AV91 = 0.d0
      AV92 = 0.d0

      AV121 = 0.d0
      AV122 = 0.d0

      A4C = 0.d0
      A3C = 0.d0
      A4D = 0.d0
      A3D = 0.d0

      AAK8 = 0.d0
      AAK9 = 0.d0
      
      DO II=1,NDIMTO                                                          
         KL=II-KTRON1                                                           
         COEFA= A(IVTOR,II,IK)                                                  
         COFAAA =A(IVTORP,II,IK)                                                
         IF (II.EQ.NDIMTO ) THEN                                                
            COEFAA=0.d0                                                         
         ELSE                                                                   
            COEFAA=A(IVTORP,II+1,IK)                                            
         ENDIF                                                                  
         IF(II.EQ.1) THEN                                                       
            COEFA4=0.d0                                                         
         ELSE                                                                   
            COEFA4=A(IVTORP,II-1,IK)                                            
         ENDIF                                                                  
         IF(II.GT.(NDIMTO -2)) THEN                                             
            COEFA5=0.d0                                                         
         ELSE                                                                   
            COEFA5=A(IVTORP,II+2,IK)                                            
         ENDIF                                                                  
         IF(II.LE.2) THEN                                                       
            COEFA6=0.d0                                                         
         ELSE                                                                   
            COEFA6=A(IVTORP,II-2,IK)                                            
         ENDIF                                                                  

         if (ii.le.3) then
            coefklm3 = 0.d0
         else
            coefklm3 = A(IVTORP, II-3, IK)
         endif
           
         if (ii.gt.(ndimto-3)) then
            coefklp3 = 0.d0
         else
            coefklp3 = A(IVTORP, II+3, IK)
         endif

         if (ii.gt.4) then
            coefklm4 = A(IVTORP, II-4, IK)
         else
            coefklm4 = 0.d0
         endif
           
         if (ii.lt.(ndimto-4)) then
            coefklp4 = 0.d0
         else
            coefklp4 = A(IVTORP, II+4, IK)
         endif

C     RHOEFF=RHO+RHOJ*XNN+RHOK*K**2                                         
         RHOEFF=RHO                                                             
         RHOKSG =3.d0*KL+ISIG+RHOEFF*K                                          
         RHOKK=3.d0*(KL+1)+ISIG+RHOEFF*K                                        
         RHOKK2=3.d0*(KL+2)+ISIG+RHOEFF*K                                       
         RHOKM1=3.d0*(KL-1)+ISIG+RHOEFF*K                                      
         RHOKM2=3.d0*(KL-2)+ISIG+RHOEFF*K                                       
         A1=A1+COEFA*COEFAA                                                     
         A2=A2+(COEFA*COFAAA )*(RHOKSG **2)                                     
         A3=A3+COEFA*COFAAA *RHOKSG                                             
         A4=A4+COEFA*COFAAA *RHOKSG **3                                         
         A4B=A4B+COEFA*COFAAA *RHOKSG **5                                       
         A5=A5+COEFA*COFAAA *RHOKSG **4                                         
         A5B=A5B+COEFA*COFAAA *RHOKSG **6                                       
         A6=A6+COEFA*COFAAA                                                     
         A7=A7+COEFA*COEFA4                                                     
         A8=A8+COEFA*COEFA4*RHOKSG                                              
         A9=A9+COEFA*COEFAA*RHOKSG                                              

         A8B=A8B+COEFA*COEFA6*RHOKSG                                            
         A9B=A9B+COEFA*COEFA5*RHOKSG                                            
         A14B=A14B+COEFA*COEFA5*RHOKK2                                          
         A15B=A15B+COEFA*COEFA6*RHOKM2                                          

         A8C=A8C+COEFA*COEFA6*RHOKSG**2                                            
         A9C=A9C+COEFA*COEFA5*RHOKSG**2                                            
         A14C=A14C+COEFA*COEFA5*RHOKK2**2                                          
         A15C=A15C+COEFA*COEFA6*RHOKM2**2                                          

         A10=A10+COEFA*COEFA4*RHOKSG **2                                        
         A11=A11+COEFA*COEFAA*RHOKSG **2                                        
         A12=A12+COEFA*COEFA5                                                   
         A13=A13+COEFA*COEFA6                                                   
         A14=A14+COEFA*COEFAA*RHOKK                                             

         A15=A15+COEFA*COEFA4*RHOKM1                                            

         A16=A16+COEFA*COEFAA*RHOKK**2                                          
         A17=A17+COEFA*COEFA4*RHOKM1**2                                         
         
         AV91 = AV91 + COEFA * COEFKLP3
         AV92 = AV92 + COEFA * COEFKLM3

         AV121 = AV121 + COEFA * COEFKLP4
         AV122 = AV122 + COEFA * COEFKLM4

         A4C=A4C+COEFA*COFAAA *RHOKSG **8                                       
         A3C=A3C+COEFA*COFAAA *RHOKSG **7                                       

         A4D=A4D+COEFA*COFAAA *RHOKSG **10                                       
         A3D=A3D+COEFA*COFAAA *RHOKSG **9                                      

         AAK8 = AAK8 + 2.d0 * COEFA*COFAAA * RHOKSG**4
         AAK8 = AAK8 - 0.5d0 * COEFA*COEFAA * RHOKSG**4
         AAK8 = AAK8 - 0.5d0 * COEFA*COEFAA * RHOKK**4
         AAK8 = AAK8 - 0.5d0 * COEFA*COEFA4 * RHOKSG**4
         AAK8 = AAK8 - 0.5d0 * COEFA*COEFA4 * RHOKM1**4

         AAK9 = AAK9 + 2.d0 * COEFA*COFAAA * RHOKSG**3
         AAK9 = AAK9 - 0.5d0 * COEFA*COEFA4 * RHOKSG**3
         AAK9 = AAK9 - 0.5d0 * COEFA*COEFAA * RHOKSG**3
         AAK9 = AAK9 - 0.5d0 * COEFA*COEFA4 * RHOKM1**3
         AAK9 = AAK9 - 0.5d0 * COEFA*COEFAA * RHOKK**3
      end do
C     PRINT*,' A1= ',A1,' A2= ',A2,' A3= ',A3,' A4= ',A4,' A5= ',               
C     1A5                                                                       
      H0=(0.5d0*(B+C)*(XNN-XK**2)+OA*XK**2)*A6                                
      H0=H0+(FV+OFV*XNN)*(A6-0.5d0*A7-0.5d0*A1)*XNN                     
      H0=H0+(ANV*XNN)*(A6-0.5d0*A13-0.5d0*A12)                              
      H0=H0+(A6-0.5d0*A12-0.5d0*A13)*(XK**2)*BK2                                 
      H0=H0+GV*A2*XNN                                                      
      H0=H0+GVJ*A2*(XNN**2)
      H0=H0+ALV*K*A3*XNN
      H0=H0-XK**2*DJK*XNN*A6                                                 
      H0=H0-XNN*DJ*XNN*A6                                               
      H0=H0+AK1*XK**3*A3                                                         
      H0=H0+AK1K*XK**5*A3                                                        
      H0=H0+AK1J*XK**3*A3*XNN
      H0=H0+AK2*XK**2*A2                                                         
      H0=H0+AK2K*XK**4*A2                                                        
      H0=H0+AK2J*XK**2*A2*XNN
C     PRINT*,' H0= ',H0                                                         
      H0=H0+AK3*XK*A4                                                            
      H0=H0+AK3B*XK*A4B                                                          
      H0=H0+AK3J*XK*A4*XNN
      H0=H0+AK3K*XK*A4*(K**2)                                                    
      H0=H0+AK4*A5                                                              
      H0=H0+AK4B*A5B                                                            
      H0=H0+AMV*A5*XNN
      H0=H0+BK1*A5*XK**2
      H0=H0+AK5*XK**2*(A6-0.5d0*A7-0.5d0*A1)                                     
      H0=H0+AK5K*XK**4*(A6-0.5d0*A7-0.5d0*A1)                                    
      H0=H0+AK5J*XK**2*(A6-0.5d0*A7-0.5d0*A1)*XNN
      H0=H0+AK6*XK*(A3-0.5d0*A8-0.5d0*A9+A3-0.5d0*A14-0.5d0*A15)               
      H0=H0+AK66*XK*(A3-0.5d0*A8B-0.5d0*A9B+A3-0.5d0*A14B-0.5d0*A15B)           
      H0=H0+AK6K*XK*(A3-0.5d0*A8-0.5d0*A9+A3-0.5d0*A14-0.5d0*A15)*XK**2           
      H0=H0+AK6J*XK*(A3-0.5d0*A8-0.5d0*A9+A3-0.5d0*A14-0.5d0*A15)*XNN
      H0=H0+AK7*(A2-0.5d0*A10-0.5d0*A11+A2-0.5d0*A16-0.5d0*A17)                 
      H0=H0+AK7K*(A2-0.5d0*A10-0.5d0*A11+A2-0.5d0*A16-0.5d0*A17)*XK**2           
      H0=H0+AK7J*(A2-0.5d0*A10-0.5d0*A11+A2-0.5d0*A16-0.5d0*A17)*XNN        
C     PRINT*,' AK7=',AK7,' A1=',A1,'RHOKSIG=',RHOKSIG,'RHOKK=',                 
C     1RHOKK,'N=',N,'K=',K,'IVTOR=',IVTOR                                       
      H0=H0-DK*XK**4*A6                                                          
      H0=H0+HJ*(XNN)**3*A6                                                  
      H0=H0+HJK*(XNN*XK)**2*A6                                               
      H0=H0+HKJ*XNN*XK**4*A6                                                 

      H0=H0+HK*XK**6*A6                                                         

      H0=H0+AK76*(A3-0.5d0*A8C-0.5d0*A9C+A3-0.5d0*A14C-0.5d0*A15C)            

      H0 = H0 + LJ*A6*(XNN)**4
      H0 = H0 + LK*A6*XK**8
      H0 = H0 + LJKKK*A6*(XNN)*XK**6
      H0 = H0 + LJJKK*A6*(XNN)**2*XK**4
      H0 = H0 + LJJJK*A6*(XNN)**3*XK**2

      H0 = H0 + GJ6 * A2 * XNN**3
      H0 = H0 + GJ4K2 * A2 * XNN**2 * XK**2                                                
      H0 = H0 + GJ2K4 * A2 * XNN * XK**4                                         
      H0 = H0 + GK6 * A2 * XK**6

      H0 = H0 + FJ6 *(A6-0.5d0*A7-0.5d0*A1) * XNN**3
      H0 = H0 + FJ4K2 * XK**2 *(A6-0.5d0*A7-0.5d0*A1) * XNN**2
      H0 = H0 + FJ2K4 * XK**4 *(A6-0.5d0*A7-0.5d0*A1) * XNN
      H0 = H0 + FK6 * XK**6 *(A6-0.5d0*A7-0.5d0*A1)                                    

      H0= H0 + AK6KK * XK * 
     & (A3-0.5d0*A8-0.5d0*A9+A3-0.5d0*A14-0.5d0*A15)
     &     * XK**4           

      H0= H0 + AK6JK * K * (A3-0.5d0*A8-0.5d0*A9+A3-0.5d0*A14-0.5d0*A15)
     &     * XK**2 * XNN

      H0= H0 + AK6JJ * K * (A3-0.5d0*A8-0.5d0*A9+A3-0.5d0*A14-0.5d0*A15)
     &     * XNN**2

      H0 = H0 + AK7KK * (A2-0.5d0*A10-0.5d0*A11+A2-0.5d0*A16-0.5d0*A17)
     &     * XK**4           

      H0 = H0 + AK7JK * (A2-0.5d0*A10-0.5d0*A11+A2-0.5d0*A16-0.5d0*A17)
     &     * XK**2 * XNN

      H0 = H0 + AK7JJ * (A2-0.5d0*A10-0.5d0*A11+A2-0.5d0*A16-0.5d0*A17)
     &     * XNN**2

      H0 = H0 + AK3BK * XK * A4B * K**2
      H0 = H0 + AK3BJ * XK * A4B * XNN
      H0 = H0 + AK3BK6 * XK * A4B * K**6
      
      H0 = H0 + AK4BJ * A5B * XNN
      H0 = H0 + AK4BK * A5B * XK**2
      H0 = H0 + AK4BK4 * A5B * XK**4
      H0 = H0 + AK4BK6 * A5B * XK**6
  
      H0 = H0 + AK3BK4 * XK * A4B * XK**4

      H0=H0 + K*A3 * AGJ6 * XNN**3
      H0=H0 + K*A3 * AGJ4K2 * XK**2 * XNN**2
      H0=H0 + K*A3 * AGJ2K4 * XK**4 * XNN
      H0=H0 + K*A3 * AGK6 * XK**6

      H0 = H0 + (A6 - 0.5d0*AV91 - 0.5d0*AV92) * V9J * XNN
      H0 = H0 + (A6 - 0.5d0*AV91 - 0.5d0*AV92) * V9K * XK**2

      H0=H0 + (A6-0.5d0*A13-0.5d0*A12) * C6J4 * XNN**2
      H0=H0 + (A6-0.5d0*A13-0.5d0*A12) * C6J2K2 * XNN * XK**2
      H0=H0 + (A6-0.5d0*A13-0.5d0*A12) * C6K4 * XK**4

      H0=H0 + A5 * G4J4 * XNN**2
      H0=H0 + A5 * G4J2K2 * XNN * XK**2
      H0=H0 + A5 * G4K4 * XK**4

      H0=H0 + XK * A4 * GAJ4 * XNN**2
      H0=H0 + XK * A4 * GAJ2K2 * XNN * XK**2
      H0=H0 + XK * A4 * GAK4 * XK**4

      H0 = H0 + XNN * V12J * (A6 - 0.5d0*(AV121 + AV122) )

      H0 = H0 + AK4C * A4C
      H0 = H0 + AK4CK * A4C * XK**2
      H0 = H0 + AK4CK4 * A4C * XK**4

      H0 = H0 + AK3C * A3C * XK
      H0 = H0 + AK3CK * A3C * XK**3
      H0 = H0 + AK3CK4 * A3C * XK * XK**4

      H0  = H0 + AK4D * A4D
      H0  = H0 + AK3D * A3D * XK


      H0 = H0 + AK8 * AAK8
      H0 = H0 + AK8J * AAK8 * XNN
      H0 = H0 + AK8K * AAK8 * XK**2

      H0 = H0 + AK9 * AAK9 * XK
      H0 = H0 + AK9J * AAK9 * XK * XNN
      H0 = H0 + AK9K * AAK9 * XK * XK**2

      H0 = H0 + 0.25d0*A6*BC4*(XNN-XK*(XK+1))*(XNN-(XK+1)*(XK+2))
      H0 = H0 + 0.25d0*A6*BC4*(XNN-XK*(XK-1))*(XNN-(XK-1)*(XK-2))

C     PRINT*,' CONV   = ',CONV                                                  
C     H0 = H0/CONV                                                              
C     PRINT*,' H0= ',H0                                                         
      IF(IVTOR.EQ.IVTORP) H0=H0+ETOR(IK,IVTOR)                                  
C     PRINT*,ETOR(IK,IVTOR)                                                     
C     PRINT*,' H0= ',H0                                                         
      RETURN                                                                    
      END                                                                       


C**********************************************************                     
      FUNCTION H1(K,IVTOR,KP,IVTORP,CST,LO,ETOR,A)                              
      implicit none
      save
C*********************************************************                      
C     SETS UP THE ELEMENT K/K+OR-1,VT/VT' OF THE HAMILTONIAN                    
C     (ODAB)                                                                    


      include 'common_dim.f'
      include 'common_dat.f'
      include 'common_rotor.f'
      include 'common_chara.f'
      include 'common_symb.f'
      include 'common_quant.f'
      include 'common_cst.f'

      double precision etor, a

      DIMENSION ETOR(KDMX,9),A(9,NTORMX,KDMX)                         

      integer ii, ik, isign, lo, kp, ivtorp, ivtor, nn, k
      double precision h1
      double precision rho, odab, dac, odab6, a1, a2, a3, a4, a5, a6
      double precision coefa, coefa1, coefaa, cofaa1, coefa3, coefa4
      double precision ann, coefa2, odab9, coef91, coef92, a7, a8

      double precision odabk, odabj, odac, odbc, cabk, cabj
      double precision dack, dacj, dac6, dacjj, dacjk, dackk
      double precision odabkk

      XN = N * 1.0d0
      XNN = (N*(N+1)) * 1.0d0
      XK = K * 1.0d0

      ODABKK = CST(LO, 191)

      ODAB9 = CST(LO, 161)

      DAC6 = CST(LO, 150)

      DACK = CST(LO, 143)
      DACJ = CST(LO, 145)

      DACJJ = CST(LO, 186)
      DACJK = CST(LO, 187)
      DACKK = CST(LO, 188)

      RHO=CST(LO,73)                                                            
C     RHOJ=CST(LO,60)                                                           
C     RHOK=CST(LO,61)                                                           
      ODAB=CST(LO,28)                                                           
      DAC=CST(LO,44)                                                            
      ODAB6=CST(LO,58)
      
      ODABK = CST(LO, 117)
      ODABJ = CST(LO, 130)

      ODAC = CST(LO, 131)
      ODBC = CST(LO, 132)

      CABK = CST(LO, 137)
      CABJ = CST(LO, 138)

      
      IK=K+N   +1                                                               
      ISIGN=1                                                                   
      IF (KP.EQ.(K-1)) ISIGN=-1                                                 
      A1=0.d0                                                                  
      A2=0.d0                                                                  
      A3=0.d0                                                                  
      A4=0.d0                                                                  
      A5=0.d0
      A6=0.d0
      A7=0.d0
      A8=0.d0
      DO II=1,NDIMTO                                                         
         COEFA=A(IVTOR,II,IK)                                                   
         COEFA1=A(IVTORP,II,IK+ISIGN)                                           
         IF (II.EQ.NDIMTO ) THEN                                                
            COEFAA=0.d0                                                         
            COFAA1 =0.d0                                                        
         ELSE                                                                   
            COEFAA=A(IVTOR,II+1,IK)                                             
            COFAA1 =A(IVTORP,II+1,IK+ISIGN)                                     
         ENDIF                                                                  
         IF(II.GT.(NDIMTO-2)) THEN
            COEFA3=0.d0
            COEFA4=0.d0
         ELSE
            COEFA3=A(IVTOR,II+2,IK)
            COEFA4=A(IVTORP,II+2,IK+ISIGN)
         ENDIF

         if (ii.gt.(ndimto-3)) then
            coef91 = 0.d0
            coef92 = 0.d0
         else
            coef91=a(ivtor, ii+3, ik)
            coef92=a(ivtorp, ii+3, ik+isign)
         endif

         IF(II.EQ.1) THEN                                                       
            COEFA2 =0.d0                                                        
         ELSE                                                                   
            COEFA2 =A(IVTORP,II-1,IK+ISIGN)                                     
         ENDIF                                                                  
         A1=A1+COEFA*COEFA1                                                     
         A2=A2+COEFAA*COEFA1                                                    
         A3=A3+COEFA*COFAA1                                                     
         A4=A4+COEFA*COEFA2                                                     
         A5=A5+COEFA1*COEFA3
         A6=A6+COEFA*COEFA4
         A7=A7+COEFA1*COEF91
         A8=A8+COEFA*COEF92

      end do
C     PRINT*,' A1= ',A1,' A2= ',A2,' A3= ',A3                                   
      NN=N*(N+1)-K*(K+ISIGN)                                                    
      ANN=dble(NN)                                                             
      H1=ODAB*SQRT(ANN)*(K+0.5d0*ISIGN)                                         
      H1=H1*(A1-0.5d0*(A2+A3))                                                 
      H1=H1+ODAB6*SQRT(ANN)*(XK+0.5d0*ISIGN)*(A1-0.5d0*(A5+A6))
      H1=H1+(-ISIGN)*(DAC/2.d0)*SQRT(ANN)*(XK+0.5d0*ISIGN)*(A4-A3)               

      H1=H1 + ODABK * SQRT(ANN) * (A1-0.5d0*(A2+A3)) 
     &     * 0.5d0 * (XK**3 + (XK+ISIGN)**3)

      H1=H1 + ODABKK * SQRT(ANN) * (A1-0.5d0*(A2+A3)) 
     &     * 0.5d0 * (XK**5 + (XK+ISIGN)**5)

      H1=H1 + ODABJ * SQRT(ANN) * (A1-0.5d0*(A2+A3)) 
     &     * 0.5d0 * (XK + (XK+ISIGN)) * XNN

      H1=H1 + CABK * SQRT(ANN) * 0.5d0 * (XK**3 + (XK+ISIGN)**3)
     &     * (A1-0.5d0*(A5+A6))

      H1=H1 + CABJ * SQRT(ANN) * 0.5d0 * (XK + (XK+ISIGN))
     &     * (A1-0.5d0*(A5+A6)) * XNN

      H1=H1+(-ISIGN)*(DACK/2.d0)*SQRT(ANN)*
     &     0.5d0 * (XK**3 + (XK+ISIGN)**3)*(A4-A3)               

      H1=H1+(-ISIGN)*(DACKK/2.d0)*SQRT(ANN)*
     &     0.5d0 * (XK**5 + (XK+ISIGN)**5)*(A4-A3)               

      H1=H1+(-ISIGN)*(DACJK/2.d0)*SQRT(ANN)*
     &     0.5d0 * (XK**3 + (XK+ISIGN)**3)*(A4-A3) * XNN

      H1=H1+(-ISIGN)*(DACJ/2.d0)*SQRT(ANN)*
     &     (XK+0.5d0*ISIGN)*(A4-A3) * XNN

      H1=H1+(-ISIGN)*(DACJJ/2.d0)*SQRT(ANN)*
     &     (XK+0.5d0*ISIGN)*(A4-A3) * XNN**2

      H1=H1+(-ISIGN)*(ODAC)*SQRT(ANN)*
     &     (XK + (XK+ISIGN))*(A4/4.d0 - A3/4.d0 + A6/4.d0 - A5/4.d0)

      H1=H1+(-ISIGN)*(DAC6/2.d0)*SQRT(ANN)*(XK+0.5d0*ISIGN)*(A6-A5)               

      H1=H1+ODAB9*SQRT(ANN)*(XK+0.5d0*ISIGN)*(A1-0.5d0*(A7+A8))

C     H1 = H1/CONV                                                              
      RETURN                                                                    
      END                                                                       


C***********************************************************                    
      FUNCTION H2(K,IVTOR,KP,IVTORP,CST,LO,ETOR,A)                              
      implicit none
      save
C                                                                               
C     SETS UP THE K/K+OR-1 ELEMENTS(DAB)                                        


      include 'common_dim.f'
      include 'common_dat.f'
      include 'common_chara.f'
      include 'common_rotor.f'
      include 'common_symb.f'
      include 'common_quant.f'
      include 'common_cst.f'

      double precision h2, etor, a
      integer ivtor, kp, ivtorp, lo, ik, isign, ii, kl, nn, k

      DIMENSION ETOR(KDMX,9),A(9,NTORMX,KDMX)                         

      double precision rho, dab, dabj, dabk, odelta, delta, deltab
      double precision odeltb, a1, a2, a3, a2b, a3b, coefa, coefa1
      double precision rhoeff, rhoksg, rhoksh, rhoks2, rhoks3, rhoksi
      double precision rhoks4, ann
      double precision HABK, HABJ, HABJK, deltaj, deltak, a4, a5
      double precision dagj, dagk, dg4j, dg4k, gb, a6
      double precision deltjk, deltkk, a7, ak7ab, a8, a9, a10
      double precision rhoks5, rhoks6, rhoks7, rhoks8
      double precision coefa2, coefa3, a3c

      XN = N * 1.0d0
      XNN = (N*(N+1)) * 1.0d0
      XK = K * 1.0d0

      GB = CST(LO, 152)

      AK7AB = CST(LO, 171)

      RHO=CST(LO,73)                                                            
      DAB=CST(LO,27)                                                            
      DABJ=CST(LO,29)                                                           
      DABK=CST(LO,30)                                                           
      ODELTA=CST(LO,47)                                                         
      DELTA=CST(LO,48)                                                          
      DELTAB=CST(LO,69)
      ODELTB=CST(LO,70)

      HABK = CST(LO, 87)
      HABJ = CST(LO, 88)
      HABJK = CST(LO, 89)

      DELTAJ = CST(LO, 133)
      DELTAK = CST(LO, 134)

      DELTJK = CST(LO, 159)
      DELTKK = CST(LO, 160)

      DAGJ = CST(LO, 135)
      DAGK = CST(lo, 136)

      DG4J = CST(LO, 141)
      DG4K = CST(LO, 142)

C     RHOJ=CST(LO,60)                                                           
C     RHOK=CST(LO,61)                                                           
      IK=K+N   +1                                                               
      ISIGN=1                                                                   
      IF (KP.EQ.(K-1)) ISIGN=-1                                                 
      A1=0.d0                                                                   
      A2=0.d0                                                                   
      A3=0.d0                                                                   
      A2B=0.d0
      A3B=0.d0
      A3C=0.d0
      A4 = 0.d0
      A5 = 0.d0
      A6 = 0.d0
      A7 = 0.d0
      A8 = 0.d0
      A9 = 0.d0
      A10 = 0.d0
      DO II=1,NDIMTO                                                         
         KL=II-KTRON1                                                           
         COEFA=A(IVTOR,II,IK)                                                   
         COEFA1=A(IVTORP,II,IK+ISIGN)                                           

         if (ii.le.1) then
            coefa2 = 0.d0
         else
            coefa2 = A(IVTORP,II-1,IK+ISIGN)
         endif
         
         if ( ii.ge.ndimto ) then
            coefa3 = 0.d0
         else
            coefa3 = A(IVTORP, II+1, IK+ISIGN)
         endif

C     RHOEFF=RHO+RHOJ*XNN+RHOK*K**2                                         
         RHOEFF=RHO                                                             
         RHOKSG =(3.d0*KL+ISIG+RHOEFF*(XK+ISIGN))*((XK+ISIGN)**2)                 
         RHOKSH =(3.d0*KL+ISIG+RHOEFF*(XK+ISIGN))**3*((XK+ISIGN)**2)              
         RHOKS2 =(3.d0*KL+ISIG+RHOEFF*K)*(XK**2)                                 
         RHOKSI =(3.d0*KL+ISIG+RHOEFF*K)**3*(XK**2)                              
         RHOKS3 =3.d0*KL+ISIG+RHOEFF*K                                         
         RHOKS4 =3.d0*KL+ISIG+RHOEFF*(K+ISIGN)                                  
         RHOKS5 =3.d0*(KL+1)+ISIG+RHOEFF*(K)
         RHOKS6 =3.d0*(KL+1)+ISIG+RHOEFF*(K+ISIGN)                                  
         RHOKS7 =3.d0*(KL-1)+ISIG+RHOEFF*(K)
         RHOKS8 =3.d0*(KL-1)+ISIG+RHOEFF*(K+ISIGN)                                  

         A1=A1+COEFA*COEFA1                                                     
         A2=A2+COEFA*COEFA1*(RHOKSG +RHOKS2 )                                   
         A2B=A2B+COEFA*COEFA1*(RHOKSH +RHOKSI )                                 
         A3=A3+COEFA*COEFA1*(RHOKS3**2 * XK+RHOKS4 **2*(XK+ISIGN))                 
         A3B=A3B+COEFA*COEFA1*(RHOKS3**4 * XK+RHOKS4 **4*(XK+ISIGN))               
         A3C=A3C+COEFA*COEFA1*
     &     (RHOKS3**4 * XK**3+RHOKS4 **4*(XK+ISIGN)**3)               
         A4=A4+COEFA*COEFA1*(RHOKS3**2 * XK**3+RHOKS4 **2*(XK+ISIGN)**3)           
         A5=A5+COEFA*COEFA1*(RHOKSG*(XK+ISIGN)**2 + RHOKS2*XK**2 )                 
         A6 = A6 + COEFA*COEFA1*(RHOKS3 + RHOKS4)
         A7=A7+COEFA*COEFA1*
     &     (RHOKS3**2 * XK**5 + RHOKS4**2 * (XK+ISIGN)**5)           

         A8 = A8 + COEFA*COEFA1*(RHOKS3**2*K + RHOKS4**2*(K+ISIGN))
         A9 = A9 + COEFA*coefa3*(rhoks3**2*K + rhoks4**2*(K+ISIGN))
         A9 = A9 + COEFA*coefa3*(rhoks5**2*K + rhoks6**2*(K+ISIGN))
         A10 = A10 + COEFA*coefa2*(rhoks3**2*K + rhoks4**2*(K+ISIGN))
         A10 = A10 + COEFA*coefa2*(rhoks7**2*K + rhoks8**2*(K+ISIGN))
      end do
C     PRINT*,' A1= ',A1                                                         
C     PRINT*,' A1= ',A1                                                         
      NN=N*(N+1)-K*(K+ISIGN)                                                    
      ANN=dble(NN)                                                             
C     H2=DAB*SQRT(ANN)*(K+0.5*ISIGN)*A1+DABJ*SQRT(ANN)*(K+0.5*ISIGN)            
C    *     *A1*XNN+DABK*SQRT(ANN)*(K+0.5*ISIGN)*A1*XK**2                     
C******DABK TERM MADE HERMITIAN BY JTH 18 APR 92.                               
C     H2=DAB*SQRT(ANN)*(K+0.5*ISIGN)*A1+DABJ*SQRT(ANN)*(K+0.5*ISIGN)            
C    *     *A1*XNN+DABK*SQRT(ANN)*(K+0.5*ISIGN)*A1*K**2                     
      H2=DAB*SQRT(ANN)*(K+0.5d0*ISIGN)*A1                                       
     * +DABJ*SQRT(ANN)*(K+0.5d0*ISIGN)*A1*XNN                               
     * +DABK*SQRT(ANN)*0.5d0*(XK**3+(K+1.d0*ISIGN)**3)*A1                        
      H2=H2+ODELTA*A2*SQRT(ANN)*0.5d0                                           
      H2=H2+ODELTB*A2B*SQRT(ANN)*0.5d0                                          
      H2=H2+DELTA*A3*SQRT(ANN)*0.5d0                                            
      H2=H2+DELTAB*A3B*SQRT(ANN)*0.5d0                                          

      H2=H2+DELTAK*A4*SQRT(ANN)*0.5d0                                            
      H2=H2+DELTAJ*A3*SQRT(ANN)*0.5d0*XNN        

      H2=H2+DELTKK*A7*SQRT(ANN)*0.5d0                                            
      H2=H2+DELTJK*A4*SQRT(ANN)*0.5d0*XNN

      H2=H2+DAGK*A5*SQRT(ANN)*0.5d0                                           
      H2=H2+DAGJ*A2*SQRT(ANN)*0.5d0*XNN

      H2 = H2 + A1*0.5d0*HABK*SQRT(ANN)*(XK**5 + (XK+ISIGN)**5)
      H2 = H2 + A1*0.5d0*HABJ*SQRT(ANN)*(XNN)**2 * (XK + (XK+ISIGN))
      H2 = H2 + A1*0.5d0*HABJK*SQRT(ANN)*(XNN)*(XK**3 + (XK+ISIGN)**3)

      H2 = H2 + GB*SQRT(ANN)*A6 * 0.5d0

      H2 = H2 + AK7AB*(2.d0*A8-0.5d0*A9-0.5d0*A10)*SQRT(ANN)*0.5d0

      H2=H2+DG4J*A3B*SQRT(ANN)*0.5d0*XNN
      H2=H2+DG4K*A3C*SQRT(ANN)*0.5d0

C     H2 = H2/CONV                                                              
C     PRINT*,' H2= ',H2                                                         
      RETURN                                                                    
      END                                                                       


C************************************************************                   
      FUNCTION H3(K,IVTOR,KP,IVTORP,CST,LO,ETOR,A)                              
      implicit none
      save
C                                                                               
C     SETS UP THE K/K+-2 ELEMENTS                                               


      include 'common_dim.f'
      include 'common_dat.f'
      include 'common_chara.f'
      include 'common_rotor.f'
      include 'common_symb.f'
      include 'common_quant.f'
      include 'common_cst.f'

      double precision h3, a
      integer ivtor, kp, ivtorp, lo, k
      double precision etor

      DIMENSION A(9,NTORMX,KDMX),ETOR(KDMX,9)                         

      integer ik, isign, ii, kl
      double precision oa, dj, djk, dk, hk, ohj, ohjk, ohk, fv, ofv, gv
      double precision ak1, ak2, odeln, odelk, hj, hjk, hkj, alv, olv
      double precision c1, c2, c3, c11, rho, c2j, gvj, ak2j, c1j
      double precision dbc, c4, c2k, c1k, c4j, c4k, c12, ak3, ak4, ak5
      double precision ak6, ak7, dab, odab, a1, a2, a3, a4, a5, a6, a7
      double precision a8, a9, a10, a11, coefa, coefa2, coefa1, cofaa2
      double precision coefa3, coefa4, coefaa4, rhoeff, rhoksg, rhoks2
      double precision rhoks3, rhoks4, rhoks5, rhoks6, rhoks7, rhok2
      double precision ann, all, b, c
      double precision OLKKK, OLJKK, OLJJK, OLJJJ
      double precision AK7BC
      double precision rhokl1, rhokl2, rhokl3, rhokl4
      double precision dbcj, c9bc, a12, a13, coefa5, coefaa5
      double precision c6bcj, ag3bcj, g4bcj
      double precision c6bck, ag3bck, g4bck
      double precision c6bcjk, g3bcjk, g4bcjk
      double precision c6bckk, g3bckk, g4bckk
      double precision ak4bc, ak3bc
      double precision ab1, ab2, ab3, ab4, ab5
      double precision ab6, ab7, ab8, ab9, ab10, a10bk, a5bk
      double precision ab11, ab12, a10b, a5b, dbck, dbcjk, dbc6
      double precision c1jk, c2jk, c4jk, c1kk, c2kk, c4kk, a11k
      double precision a2k, a2kk

      XN = N * 1.0d0
      XNN = (N*(N+1)) * 1.0d0
      XK = K * 1.0d0

      C6BCJK = CST(LO, 195)
      G3BCJK = CST(LO, 196)
      G4BCJK = CST(LO, 197)
      C6BCKK = CST(LO, 201)
      G3BCKK = CST(LO, 202)
      G4BCKK = CST(LO, 203)

      C2JK = CST(LO, 192)
      C1JK = CST(LO, 193)    
      C4JK = CST(LO, 194)
      C2KK = CST(LO, 198)
      C1KK = CST(LO, 199)
      C4KK = CST(LO, 200)

      DBC6 = CST(LO, 151)

      DBCK = CST(LO, 189)
      DBCJK = CST(LO, 190)

      C6BCJ = CST(LO, 163)
      AG3BCJ = CST(LO, 165)
      G4BCJ = CST(LO, 167)

      C6BCK = CST(LO, 164)
      AG3BCK = CST(LO, 166)
      G4BCK = CST(LO, 168)
      
      C9BC = CST(LO, 162)
      DBCJ = CST(LO,  144)

      OA=CST(LO,1)                                                              
      B=CST(LO,2)                                                               
      C=CST(LO,3)                                                               
      DJ=CST(LO,4)                                                              
      DJK=CST(LO,5)                                                             
      DK=CST(LO,6)                                                              
      ODELN=CST(LO,7)                                                           
      ODELK=CST(LO,8)                                                           
      HJ=CST(LO,9)                                                              
      HJK=CST(LO,10)                                                            
      HKJ=CST(LO,11)                                                            
      HK=CST(LO,12)                                                             
      OHJ=CST(LO,13)                                                            
      OHJK=CST(LO,14)                                                           
      OHK=CST(LO,15)                                                            
      FV=CST(LO,16)                                                             
      OFV=CST(LO,17)                                                            
      GV=CST(LO,18)                                                             
      ALV=CST(LO,19)                                                            
      AK1=CST(LO,20)                                                            
      AK2=CST(LO,21)                                                            
      AK3=CST(LO,22)                                                            
      AK4=CST(LO,23)                                                            
      AK5=CST(LO,24)                                                            
      AK6=CST(LO,25)                                                            
      AK7=CST(LO,26)                                                            
      DAB=CST(LO,27)                                                            
      ODAB=CST(LO,28)                                                           
      OLV=CST(LO,31)                                                            
      C1=CST(LO,32)                                                             
      C2=CST(LO,33)                                                             
      C3=CST(LO,34)                                                             
      C11=CST(LO,57)

      RHO=CST(LO,73)                                                            
C     RHOJ=CST(LO,60)                                                           
C     RHOK=CST(LO,61)                                                           
      C2J=CST(LO,37)                                                            
      GVJ=CST(LO,38)                                                            
      AK2J=CST(LO,39)                                                           
      C1J=CST(LO,40)                                                            
      DBC=CST(LO,45)                                                            
      C4=CST(LO,46)                                                             
      C2K=CST(LO,65)
      C1K=CST(LO,66)
      C4J=CST(LO,67)
      C4K=CST(LO,68)
      C12=CST(LO,71)

      OLKKK = CST(LO, 83)
      OLJKK = CST(LO, 84)
      OLJJK = CST(LO, 85)
      OLJJJ = CST(LO, 86)

      AK7BC = CST(LO, 114)

      AK4BC = CST(LO, 169)
      AK3BC = CST(LO, 170)

      IK=K+N   +1                                                               
      ISIGN=2                                                                   
C     PRINT*,' C1= ',C1,' C2= ',C2,' C3= ',C3,' OHJK= ',                        
C    1 OHJK,' ODELN= ',ODELN,' ODELK= ',ODELK,' OHJ= ',OHJ,' OHK= '             
C    2 ,OHK                                                                     
      IF (KP.EQ.(K-2)) ISIGN=-2                                                 
      A1=0.d0                                                                   
      A2=0.d0                                                                   
      A2K=0.d0                                                                   
      A2KK=0.d0                                                                   
      A3=0.d0                                                                   
      A4=0.d0                                                                   
      A5=0.d0                                                                   
      A6=0.d0                                                                   
      A7=0.d0                                                                   
      A8=0.d0                                                                   
      A9=0.d0                                                                   
      A10=0.d0
      A11=0.d0
      A11K=0.d0
      A12 = 0.d0
      A13 = 0.d0

      ab1 = 0.d0
      ab2 = 0.d0
      ab3 = 0.d0
      ab4 = 0.d0
      ab5 = 0.d0
      ab6 = 0.d0
      ab7 = 0.d0
      ab8 = 0.d0
      ab9 = 0.d0
      ab10 = 0.d0
      ab11 = 0.d0
      ab12 = 0.d0
      a10b = 0.d0
      a5b = 0.d0
      a10bk = 0.d0
      a5bk = 0.d0

      DO II=1,NDIMTO                                                         
         KL=II-KTRON1                                                           
         COEFA=A(IVTOR,II,IK)                                                   
         COEFA2=A(IVTORP,II,IK+ISIGN)                                           
         IF (II.EQ.NDIMTO ) THEN                                                
            COEFA1=0.d0                                                         
            COFAA2 =0.d0                                                        
         ELSE                                                                   
            COEFA1=A(IVTOR,II+1,IK)                                             
            COFAA2 =A(IVTORP,II+1,IK+ISIGN)                                     
         ENDIF                                                                  
         IF(II.EQ.1) THEN                                                       
            COEFA3 =0.d0                                                        
         ELSE                                                                   
            COEFA3 =A(IVTORP,II-1,IK+ISIGN)                                     
         ENDIF                                                                  
         IF(II.GT.(NDIMTO-2)) THEN
            COEFA4=0.d0
            COEFAA4=0.d0
         ELSE
            COEFA4=A(IVTOR,II+2,IK)
            COEFAA4=A(IVTORP,II+2,IK+ISIGN)
         ENDIF

         IF(II.GT.(NDIMTO-3)) THEN
            COEFA5=0.d0
            COEFAA5=0.d0
         ELSE
            COEFA5=A(IVTOR,II+3,IK)
            COEFAA5=A(IVTORP,II+3,IK+ISIGN)
         ENDIF

C     RHOEFF=RHO+RHOJ*XNN+RHOK*K**2                                         
         RHOEFF=RHO                                                             
         RHOKSG =3.d0*KL+ISIG+RHOEFF*K                                          
         RHOKS2 =(3.d0*KL+ISIG+RHOEFF*(K+ISIGN))*(K+ISIGN)                      
         RHOKS3 =(3.d0*KL+ISIG+RHOEFF*K)*(K)                                    
         RHOKS4 =(3.d0*KL+ISIG+RHOEFF*K)**3                                     
         RHOKS4 =RHOKS4*K                                                       
         RHOKS5 =(3.d0*KL+ISIG+RHOEFF*(K+ISIGN))**3                             
         RHOKS5=RHOKS5*(K+ISIGN)
C         RHOKS6 =(3.d0*KL+ISIG+RHOEFF*K)*(K)                                    
C         RHOKS6=RHOKS6*(XK+ISIGN)**2
C         RHOKS7 =(3.d0*KL+ISIG+RHOEFF*(XK+ISIGN))*(XK+ISIGN)                      
C         RHOKS7=RHOKS7*(XK)**2

         RHOKS6 =(3.d0*KL+ISIG+RHOEFF*K) * XK**3
         RHOKS7 =(3.d0*KL+ISIG+RHOEFF*(XK+ISIGN)) * (XK+ISIGN)**3

         RHOK2=3.d0*KL+ISIG+RHOEFF*(XK+ISIGN)                                    

         rhokl1 = 3.d0*(KL+1) + isig + rhoeff*(K)
         rhokl2 = 3.d0*(KL-1) + isig + rhoeff*(K)
         rhokl3 = 3.d0*(KL+1) + isig + rhoeff*(K+ISIGN)
         rhokl4 = 3.d0*(KL-1) + isig + rhoeff*(K+ISIGN)

         ab1 = ab1 + rhoksg**2 * coefa * coefa2
         ab2 = ab2 + rhok2**2 * coefa * coefa2
         ab3 = ab3 + rhoksg**2 * coefa * cofaa2
         ab4 = ab4 + rhoksg**2 * coefa * coefa3
         ab5 = ab5 + rhok2**2 * coefa * cofaa2
         ab6 = ab6 + rhok2**2 * coefa * coefa3
         ab7 = ab7 + rhokl1**2 * coefa * cofaa2
         ab8 = ab8 + rhokl2**2 * coefa * coefa3
         ab9 = ab9 + rhokl3**2 * coefa * cofaa2
         ab10 = ab10 + rhokl4**2 * coefa * coefa3

         ab11 = ab11 + (rhoksg**6 + rhok2**6) * coefa * coefa2
         ab12 = ab12 + 
     &        (rhoksg**5*k + rhok2**5*(k+isign)) * coefa * coefa2

         A1=A1+COEFA*COEFA2                                                     
C     PRINT*,' COEFA= ',COEFA,' COEFA1= ',COEFA1,' COEFA2= ',                   
C     1COEFA2,' COEFAA2= ',COEFAA2                                              
         A2=A2+COEFA*COEFA2*(RHOKSG**2+RHOK2**2)                               
         A2K=A2K+COEFA*COEFA2*
     &        (RHOKSG**2*XK**2+RHOK2**2*(XK+ISIGN)**2) 
         A2KK=A2KK+COEFA*COEFA2*
     &        (RHOKSG**2*XK**4+RHOK2**2*(XK+ISIGN)**4) 

         A3=A3+COEFA1*COEFA2                                                    
         A4=A4+COEFA*COFAA2                                                     
         A5=A5+COEFA*COEFA2*(RHOKSG **4+RHOK2**4)                               
         A6=A6+COEFA*COEFA3                                                     
         A7=A7+COEFA*COEFA2*(RHOKS2 +RHOKS3 )                                   
         A8=A8+COEFA2*COEFA4
         A9=A9+COEFA*COEFAA4
         A10=A10+COEFA*COEFA2*(RHOKS4+RHOKS5)
         A11=A11+COEFA*COEFA2*(RHOKS6+RHOKS7)
         A11K=A11K+COEFA*COEFA2
     &        *(RHOKS6*K**2+RHOKS7*(K+ISIGN)**2)


         A12=A12+COEFA2*COEFA5
         A13=A13+COEFA*COEFAA5         

         A5B  = A5B  + COEFA * COEFA2 
     &        *(RHOKSG**4*K**2 + RHOK2**4 * (XK+ISIGN)**2)                               
         A10B = A10B + COEFA * COEFA2 
     &        *(RHOKS4*K**2 + RHOKS5 * (XK+ISIGN)**2)

         A5BK  = A5BK + COEFA * COEFA2 
     &        *(RHOKSG**4 * XK**4 + RHOK2**4 * (XK+ISIGN)**4)                               
         A10BK = A10BK+ COEFA * COEFA2 
     &        *(RHOKS4 * XK**4 + RHOKS5 * (XK+ISIGN)**4)

      end do
C     PRINT*,' A1= ',A1,' A2= ',A2,' A3= ',A3,' A4= ',A4                        
C     PRINT*,' K= ',K,' KP= ',KP,' IVTOR= ',IVTOR,' IVTORP= ',                  
C     1IVTORP,' ISIGN= ',ISIGN                                                  
      H3=(B-C)/4.d0-ODELN*XNN-0.5d0*ODELK*(XK**2+(XK+ISIGN)**2)+              
     *     OHJ*(XNN)**2+0.5d0*OHJK*XNN*(XK**2+(XK+ISIGN)**2)              
     *     +0.5d0*OHK*(XK**4+(XK+ISIGN)**4)                                       
C     PRINT*,' H3= ',H3                                                         
      H3=A1*H3                                                                  
      H3=H3+0.5d0*C1*A2                                                         
      H3=H3+0.5d0*C1J*A2*XNN                                                
C     PRINT*,' H3= ',H3                                                         
      H3=H3+0.5d0*C2*(A1-0.5d0*(A3+A4))                                         
      H3=H3+0.5d0*C2J*(A1-0.5d0*(A3+A4))*(XNN)                              
      H3=H3+0.5d0*C3*A5                                                         
      H3=H3+0.5d0*C4*A7                                                         
      H3=H3+0.5d0*C12*A10
      H3=H3+0.5d0*C4J*A7*XNN
      H3=H3+0.5d0*C1K*A2K
      H3=H3+0.5d0*C2K*(A1-0.5d0*(A3+A4))*(XK**2+(XK+ISIGN)**2)
      H3=H3+0.5d0*C4K*A11
      H3=H3+0.5d0*C11*(A1-0.5d0*(A8+A9))


      H3=H3+0.5d0*C1JK*A2K*XNN
      H3=H3+0.5d0*C2JK*(A1-0.5d0*(A3+A4))*(XK**2+(XK+ISIGN)**2)
     &     *XNN
      H3=H3+0.5d0*C4JK*A11*XNN

      H3=H3+0.5d0*C1KK*A2KK
      H3=H3+0.5d0*C2KK*(A1-0.5d0*(A3+A4))*(XK**4+(XK+ISIGN)**4)
      H3=H3+0.5d0*C4KK*A11K

      H3 = H3 + 0.5d0*A1*OLKKK*(XK**6 + (XK+ISIGN)**6)
      H3 = H3 + 0.5d0*A1*OLJKK*(XK**4 + (XK+ISIGN)**4)*(XNN)
      H3 = H3 + 0.5d0*A1*OLJJK*(XK**2 + (XK+ISIGN)**2)*(XNN)**2
      H3 = H3 + 0.5d0*A1*OLJJJ*(XNN)**3

      H3 = H3 + 0.5d0 * AK7BC * ( 2.0d0*(ab1 + ab2) -
     &     0.5d0*(ab3 + ab4 + ab5 + ab6 + ab7 + ab8 + ab9 + ab10))

      ANN=XNN-K*(K+0.5d0*ISIGN)                                             
      ALL=XNN-(K+0.5d0*ISIGN)*(K+ISIGN)                                     

      H3=H3+0.5d0 * C6BCJ * (A1-0.5d0*(A8+A9)) * XNN
      H3=H3+0.5d0 * AG3BCJ * A10 * XNN
      H3=H3+0.5d0 * G4BCJ * A5 * XNN                                              

      H3=H3+0.5d0 * C6BCK * (A1-0.5d0*(A8+A9)) * (XK**2 + (XK+ISIGN)**2)
      H3=H3+0.5d0 * AG3BCK * A10B
      H3=H3+0.5d0 * G4BCK * A5B

      H3=H3+0.5d0 * C6BCJK * (A1-0.5d0*(A8+A9)) * 
     &   (XK**2 + (XK+ISIGN)**2)
     &     *XNN
      H3=H3+0.5d0 * G3BCJK * A10B *XNN
      H3=H3+0.5d0 * G4BCJK * A5B *XNN

      H3=H3+0.5d0 * C6BCKK * (A1-0.5d0*(A8+A9)) * 
     &  (XK**4 + (XK+ISIGN)**4)
      H3=H3+0.5d0 * G3BCKK * A10BK
      H3=H3+0.5d0 * G4BCKK * A5BK

      H3 = H3 + 0.5d0 * AK4BC * AB11
      H3 = H3 + 0.5d0 * AK3BC * AB12

      H3=H3*SQRT(ANN*ALL)                                                       

C     PRINT*,' H3= ',H3                                                         
      H3=H3+((A6-A4)*(-ISIGN)*0.5d0*DBC*(1.d0/4.d0))*(SQRT(ANN*ALL)) 
      H3=H3+((A9-A8)*(-ISIGN)*0.5d0*DBC6*(1.d0/4.d0))*(SQRT(ANN*ALL)) 

      H3 = H3 + ((A6-A4) * (-ISIGN) * 0.5d0 * 
     &     DBCJ * (1.d0/4.d0))*(SQRT(ANN*ALL)) * XNN

      H3=H3+0.5d0*C9BC*(A1-0.5d0*(A12+A13))*SQRT(ANN*ALL)

      H3=H3+((A6-A4)*(-ISIGN)*0.5d0*DBCK*(1.d0/4.d0))*(SQRT(ANN*ALL)) 
     &     * (XK**2 + (XK+ISIGN)**2)

      H3=H3+((A6-A4)*(-ISIGN)*0.5d0*DBCJK*(1.d0/4.d0))*(SQRT(ANN*ALL)) 
     &     * (K**2 + (K+ISIGN)**2)*XNN



      RETURN                                                                    
      END                                                                       



C**************************************************************                 
      FUNCTION H4(K,IVTOR,KP,IVTORP,CST,LO,ETOR,A)                              
      implicit none
      save
C**************************************************************                 
C     SETS UP THE ELEMENT K/K,VT,VT'(OLV)                                       
     

      include 'common_dim.f'
      include 'common_dat.f'
      include 'common_chara.f'
      include 'common_rotor.f'
      include 'common_symb.f'
      include 'common_quant.f'
      include 'common_cst.f'

      double precision h4, etor, a
      integer ivtor, kp, ivtorp, lo, ik, ii, kl, k

      double precision rho, olv, a1, coefa, coefaa, rhoeff, rhoksg
      DIMENSION ETOR(KDMX,9),A(9,NTORMX,KDMX)                         

      XN = N * 1.0d0
      XNN = (N*(N+1)) * 1.0d0
      XK = K * 1.0d0

      RHO=CST(LO,73)                                                            
C     RHOJ=CST(LO,60)                                                           
C     RHOK=CST(LO,61)                                                           
      OLV=CST(LO,31)                                                            
      IK=K+N   +1                                                               
      A1=0.d0                                                                   
      DO II=1,NDIMTO                                                         
         KL=II-KTRON1                                                           
         COEFA=A(IVTOR,II,IK)                                                   
         COEFAA=A(IVTORP,II,IK)                                                 
C     RHOEFF=RHO+RHOJ*N*(N+1)+RHOK*XK**2                                         
         RHOEFF=RHO                                                             
         RHOKSG =3.d0*KL+ISIG+RHOEFF*XK                                         
         A1=A1+COEFA*COEFAA*RHOKSG                                              
      end do
C     H4=OLV*XNN*K**3*A1                                                    
C     OLV SHOULD BE A J(J+1) DEPENDANCE OF ALV
      H4=OLV*XNN*K*A1*XNN                                            
C     H4=H4/CONV                                                                
      RETURN                                                                    
      END                                                                       


C***********************************************************
      FUNCTION H5(K,IVTOR,KP,IVTORP,CST,LO,ETOR,A)                             
      implicit none
      save
C***********************************************************                    
C     +-4 off diagonal elements (!)                    
      

      include 'common_dim.f'
      include 'common_dat.f'
      include 'common_chara.f'
      include 'common_rotor.f'
      include 'common_symb.f'
      include 'common_quant.f'
      include 'common_cst.f'

      integer ivtor, kp, ivtorp, lo, ik, kl, ii, k
      integer isign, isign1, isign2, isign3, isign4
      double precision h5, etor, a, coefa, cofaaa
      double precision rhoeff, rho, rhoksg, rhorho
      double precision bc4, pm4fact, obc4

      double precision a1

      DIMENSION A(9,NTORMX,KDMX),ETOR(KDMX,9)                         

      XN = N * 1.0d0
      XNN = (N*(N+1)) * 1.0d0
      XK = K * 1.0d0

      RHO=CST(LO,73)                                                            
      RHORHO=CST(LO,77)
      BC4 = CST(LO,172)
      OBC4 = CST(LO, 175)
      
      IK = K+N+1
      ISIGN=1                                                                   
      IF (KP.EQ.(K-4)) ISIGN=-1                    

      isign1 = isign
      isign2 = 2*isign
      isign3 = 3*isign
      isign4 = 4*isign

      a1 = 0.d0

      DO II=1,NDIMTO                                                          
         KL=II-KTRON1                                                           
         COEFA= A(IVTOR,II,IK)                                                  
         COFAAA =A(IVTORP,II,IK+isign4)                                                

         RHOEFF=RHO                                                             
         RHOKSG =3.d0*KL+ISIG+RHOEFF*K                                          

         a1 = a1 + coefa*cofaaa
      end do

      
      H5 = 0.0d0

      H5 = H5 + BC4*A1
      H5 = H5 + OBC4*A1

C     common factors
      pm4fact = XNN-dble(K*(K+ISIGN1))
      pm4fact = pm4fact *  (XNN - dble((K+ISIGN1)*(K+ISIGN2)))
      pm4fact = pm4fact *  (XNN - dble((K+ISIGN2)*(K+ISIGN3)))
      pm4fact = pm4fact *  (XNN - dble((K+ISIGN3)*(K+ISIGN4)))

      if ( pm4fact .gt. 0.0 ) pm4fact = SQRT(pm4fact)
      if ( pm4fact .lt. 0.0 ) pm4fact = 0.0

      H5 = H5 * 0.25d0 * pm4fact

      RETURN                                                                    
      END                                                                       








