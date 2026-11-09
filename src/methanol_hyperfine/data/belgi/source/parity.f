C*************************************************************
      SUBROUTINE PARITY(A,EGVC4,II)
C     THIS SUBROUTINE CALCULATES THE TRUE PARITY OF THE CURRENT STATE
C     (+or-).(-1)J+VT = "true" parity
C****************************************************************
      
      implicit none
      save


      include 'common_dim.f'
      include 'common_rotor.f'
      include 'common_dat.f'
      include 'common_c2c.f'
      include 'common_eigen.f'
      include 'common_quant.f'
      include 'common_chara.f'

      double precision egvc4, a
      integer ii, nlow, nlow1, knp2, ivvlo, klow, ivtorp, ivtor1, ilow
      integer kbaslo, iopplo, npklo, iii, iiiopp
      double precision an00a, an00b, an00c, an01a, an01b, an01c, an0m1a
      double precision an0m1b, an0m1c, an10a, an10b, an10c, an11a, an11b
      double precision an11c, an1m1a, an1m1b, an1m1c, anm10a, anm10b
      double precision anm10c, anm11a, anm11b, anm11c, anm1ma, anm1mb
      double precision anm1mc, expplo, a2lo, a3lo, coefal, coefaal
      double precision coef4l

      DIMENSION A(9,NTORMX,KDMX)              
      DIMENSION EGVC4(KDMX,NTORMX,NDMX)
C     
      NLOW=N        
      NLOW1 = NLOW+1
C     P OR PRIME STANDS FOR LEFT OR LOWER
      KNP2=2*NLOW
      KLOW=KK  
      IVVLO=IV 
      PERC2L=0.d0
      PER23L=0.d0
C     INITIALIZATION OF THE SUM ON TORSIONAL EIGENFUNCTIONS
      AN00A=0.d0
      AN00B=0.d0
      AN00C=0.d0
      AN01A=0.d0
      AN01B=0.d0
      AN01C=0.d0
      AN0M1A=0.d0
      AN0M1B=0.d0
      AN0M1C=0.d0
      AN10A=0.d0
      AN10B=0.d0
      AN10C=0.d0
      AN11A=0.d0
      AN11B=0.d0
      AN11C=0.d0
      AN1M1A=0.d0
      AN1M1B=0.d0
      AN1M1C=0.d0
      ANM10A=0.d0
      ANM10B=0.d0
      ANM10C=0.d0
      ANM11A=0.d0
      ANM11B=0.d0
      ANM11C=0.d0
      ANM1MA=0.d0
      ANM1MB=0.d0
      ANM1MC=0.d0
C     
C     EXPECTATION VALUE OF C2(C) AND OF (23)* FOR CURRENT STATE
C     
      DO IVTORP=1,9
         DO IVTOR1=1,9
            DO ILOW=1,2*NLOW+1
               KBASLO=ILOW-NLOW-1
               IOPPLO=-KBASLO+NLOW+1
               NPKLO=NLOW+KBASLO
               EXPPLO=1.d0
               IF(MOD(NPKLO,2).NE.0) EXPPLO=-1.d0
               A2LO=0.d0
               A3LO=0.d0
               DO III=1,NDIMTO
                  IIIOPP=2*KTRONC+1-(III-1)
                  COEFAL=A(IVTORP,III,ILOW)
                  COEFAAL=A(IVTOR1,III,IOPPLO)
                  IF(ISIG.EQ.0) THEN
                     COEF4L=A(IVTOR1,IIIOPP,IOPPLO)
                  ELSE
                     COEF4L=0.d0
                  ENDIF
                  A2LO=A2LO+COEFAL*COEFAAL
                  A3LO=A3LO+COEFAL*COEF4L
               end do
               PERC2L=PERC2L+A2LO*EXPPLO*
     1              EGVC4(ILOW,IVTORP,II)*EGVC4(IOPPLO,IVTOR1,II)
               PER23L=PER23L+A3LO*EXPPLO*
     1              EGVC4(ILOW,IVTORP,II)*EGVC4(IOPPLO,IVTOR1,II)
C     IF(N.EQ.1) PRINT*,'EGVC4(K)=',EGVC4(ILOW,IVTORP,II), 
C     1   'EGVC4(-K)=',EGVC4(IOPPLO,IVTOR1,II),'ILOW=',ILOW,
C     2   'IVTORP=',IVTORP,'II=',II
C     PRINT*,'PER23L=',PER23L,'A3LO=',A3LO,'NLOW=',NLOW,'KK=',KK,'IV=',IV
C     1,'EGVC3=',EGVC3(ILOW,IVTORP),'EGVC3(-K)=',EGVC3(IOPPLO,IVTOR1)
            end do
         end do
      end do

C     PRINT*,'PER23L=',PER23L,'N=',N,'IV=',IV,'KK=',KK
      RETURN
      END
