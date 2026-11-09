
      double precision rho, rhorho, rhoeff, r2x, r4x
      double precision aaa, aa2, aa4, aa5, aa5b, fpara

      integer m, ni, loo, na, nb
                           
      double precision R2A,R2B,R2C,                
     &       R2D,                                                  
     &       R2E,                                                  
     &       R2F,                                                  
     &       R4A,R4B,                                 
     &       R4AB,R4BB,                               
     &       R5A,R5B,R5C,                
     &       R5B2,R5C2,                
     &       R5D,R5E,R5F,                
     &       R5DD,R5FF,                
     &       R5G,r2ck, r2ckk,
     &       RVA,RVB,R2EK,R2EKK,
     &     R2B1, R2B2, R2B3, R2B4, R2B5, R2B6, r7b, r8a,
     &     R2AC, R5A8, R5E8, R5G8,
     &     R2A9, R5D9, R5F9, r8b, r2g1, r2g2, r4bc, r2fk,
     &     R2BK,R2BKK
 
                           
      common /vsetdata/ 
     &     R2A(NTORMX,KDMX),R2B(NTORMX,KDMX),R2C(NTORMX,KDMX),                
     &     R2D(NTORMX,KDMX),                                                  
     &     R2E(NTORMX,KDMX),                                                  
     &     R2F(NTORMX,KDMX),                                                  
     &     R4A(NTORMX,KDMX),R4B(NTORMX,KDMX), R7B(NTORMX,KDMX),  
     &     R4AB(NTORMX,KDMX),R4BB(NTORMX,KDMX),                               
     &     R5A(NTORMX,KDMX),R5B(NTORMX,KDMX),R5C(NTORMX,KDMX),                
     &     R5B2(NTORMX,KDMX),R5C2(NTORMX,KDMX),                
     &     R5D(NTORMX,KDMX),R5E(NTORMX,KDMX),R5F(NTORMX,KDMX),                
     &     R5DD(NTORMX,KDMX),R5FF(NTORMX,KDMX),                
     &     R5G(NTORMX,KDMX),                                                  
     &     RVA(NTORMX,KDMX),RVB(NTORMX,KDMX),                                  
     &     R2B1(NTORMX,KDMX),R2B2(NTORMX,KDMX),                                  
     &     R2B3(NTORMX,KDMX),R2B4(NTORMX,KDMX),                                  
     &     R2B5(NTORMX,KDMX),R2B6(NTORMX,KDMX),                                  
     &     R8A(NTORMX,KDMX),R2AC(NTORMX, KDMX),
     &     R5A8(NTORMX, KDMX),R5E8(NTORMX, KDMX),
     &     R5G8(NTORMX, KDMX),R2A9(NTORMX, KDMX),
     &     R5D9(NTORMX, KDMX),R5F9(NTORMX, KDMX),
     &     r8b(NTORMX, KDMX), r2g1(NTORMX, KDMX), 
     &     r2g2(NTORMX, KDMX),R4BC(NTORMX,KDMX),
     &     R2FK(NTORMX,KDMX),R2EK(NTORMX,KDMX),
     &     R2EKK(NTORMX,KDMX),R2CK(NTORMX,KDMX),
     &     R2CKK(NTORMX,KDMX),R2BK(NTORMX,KDMX),
     &     R2BKK(NTORMX,KDMX),
     &     rho, rhorho, rhoeff, r2x, r4x,
     &     aaa, aa2, aa4, aa5, aa5b, fpara, 
     &     m, ni, loo, na, nb
