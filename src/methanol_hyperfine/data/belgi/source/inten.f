
C****************************************************************
      SUBROUTINE INTEN(RINT,CST,K,IVV,INK,MUABC,ICTT,CF)
C     THIS SUBROUTINE CALCULATES THE LINE STRENGTH OF THE CURRENT
C     LINE 
C****************************************************************
C      implicit double precision (a-h, o-z)
      implicit none
      save

 1337 format(4(f6.2, 1x))
      
      character*6 muabc

      include 'common_dim.f'
      include 'common_cst.f'
      include 'common_eigen.f'
      include 'common_tor.f'
      include 'common_rotor.f'
      include 'common_dat.f'
      include 'common_expdat.f'
      include 'common_c2c.f'
      include 'common_quant.f'
      include 'common_chara.f'
      include 'common_int.f'
      include 'common_symb.f'
      include 'common_tor59.f'

      DIMENSION IVV(2),INK(2),ICTT(6),CTT(6)
      double precision ctt, cf
      integer ictt, ivv, ink, nlow, k, nup, nlow1, nup1
      integer knp2, kn2, klow, kup, ivvlo, ivvup, inklo, inkup, m
      integer iopplo, npklo, ivtorp, ivtor1, ilow, kbaslo
      integer iii, iiiopp, iup, kbasup, ioppup, npkup
      integer kk2, kkp2, isign
      double precision amua0, amua3, amub0, amub3, amuc3
      double precision perc2s, per23s, eps, an00a, an00b, an00c
      double precision an01a, an01b, an01c, an0m1a, an0m1b, an0m1c
      double precision an10a, an10b, an10c, an11a, an11b, an11c
      double precision an1m1a, an1m1b, an1m1c, anm10a, anm10b, anm10c
      double precision anm11a, anm11b, anm11c, anm1ma, anm1mb, anm1mc
      double precision q, s, s2, s4, expplo, a2lo, a3lo
      double precision coefal, coefaal, coef4l, exppup, a2up, a3up
      double precision coefas, coefaas, coef4s, coefa, coefa0, coefa0b
      double precision coefa0c, cprod0, tmua, tmua3, ct, ad3j, cprod1
      double precision tmub, tmub3, tmuc3, cprodm1, s3, rint
      

      iisig = isig+1

C
C    Set-up of MUABS = 'ABCabc'
C    A for mu_a    / B for mu_b    / C for mu_c   (zero in the present case)
C    a for mu_a(3) / b for mu_b(3) / c for mu_c(3)
C
      CF=2.0
      MUABC='      '
      AMUA0=CST(1,102)
      AMUA3=CST(1,103)
      AMUB0=CST(1,104)
      AMUB3=CST(1,105)
      AMUC3=CST(1,106)
      NLOW=NOBS(1,K)
      NUP=NOBS(2,K)
      NLOW1 = NLOW+1
      NUP1 = NUP+1
C     P OR PRIME STANDS FOR LEFT OR LOWER
      KNP2=2*NLOW
      KN2=2*NUP
      Q=(KNP2+1.)*(KN2+1.)
      KLOW=KOBS(1,K)
      KUP=KOBS(2,K)
      IVVLO=IVV(1)
      IVVUP=IVV(2)
      INKLO=INK(1)
      INKUP=INK(2)
      S=0.0
      S2=0.0
      S4=0.0
      DO 1 M = 1,(6)
        CTT(M) = 0.0
        ICTT(M) = 0
    1 CONTINUE
      PERC2L=0.0
      PERC2S=0.0
      PER23L=0.
      PER23S=0.
      EPS=1.E-10
C     INITIALIZATION OF THE SUM ON TORSIONAL EIGENFUNCTIONS
         AN00A=0.0
         AN00B=0.0
         AN00C=0.0
         AN01A=0.0
         AN01B=0.0
         AN01C=0.0
         AN0M1A=0.0
         AN0M1B=0.0
         AN0M1C=0.0
         AN10A=0.0
         AN10B=0.0
         AN10C=0.0
         AN11A=0.0
         AN11B=0.0
         AN11C=0.0
         AN1M1A=0.0
         AN1M1B=0.0
         AN1M1C=0.0
         ANM10A=0.0
         ANM10B=0.0
         ANM10C=0.0
         ANM11A=0.0
         ANM11B=0.0
         ANM11C=0.0
         ANM1MA=0.0
         ANM1MB=0.0
         ANM1MC=0.0
C
C     EXPECTATION VALUE OF C2(C) AND OF (23)* FOR LOWER STATE
C
      DO 304 IVTORP=1,9
       DO 302 IVTOR1=1,9
        DO 300 ILOW=1,2*NLOW+1
         KBASLO=ILOW-NLOW-1
         IOPPLO=-KBASLO+NLOW+1
         NPKLO=NLOW+KBASLO
         EXPPLO=1.
         IF(MOD(NPKLO,2).NE.0) EXPPLO=-1.
         A2LO=0.0
         A3LO=0.0
          DO 298 III=1,NDIMTO
          IIIOPP=2*KTRONC+1-(III-1)
          COEFAL=A59(IVTORP,III,ILOW,NLOW1,IISIG)
          COEFAAL=A59(IVTOR1,III,IOPPLO,NLOW1,IISIG)

          IF(ISIG.EQ.0) THEN
             COEF4L=A59(IVTOR1,IIIOPP,IOPPLO,NLOW1,IISIG)
          ELSE
             COEF4L=0.0
          ENDIF
          
          A2LO=A2LO+COEFAL*COEFAAL
          A3LO=A3LO+COEFAL*COEF4L
 298   CONTINUE
       PERC2L=PERC2L+A2LO*EXPPLO*
     1      VECTPR(IVVLO,INKLO,ILOW,IVTORP)*
     2      VECTPR(IVVLO,INKLO,IOPPLO,IVTOR1)
       PER23L=PER23L+A3LO*EXPPLO*
     1      VECTPR(IVVLO,INKLO,ILOW,IVTORP)*
     2      VECTPR(IVVLO,INKLO,IOPPLO,IVTOR1)
C     PRINT*,'PER23L=',PER23L,'A3LO=',A3LO,'NLOW=',NLOW

300      CONTINUE
302     CONTINUE
304   CONTINUE
C         print 1337, perc2l, per23l, a2lo, a3lo


C      print *, perc2l, ' ', per23l, ' ', a2lo, ' ' , a3lo


C
C     EXPECTATION VALUE OF C2(C) AND OF (23)* FOR UPPER STATE
C
      DO 404 IVTORP=1,9
       DO 402 IVTOR1=1,9
        DO 400 IUP=1,2*NUP+1
         KBASUP=IUP-NUP-1
         IOPPUP=-KBASUP+NUP+1
         NPKUP=NUP+KBASUP
         EXPPUP=1.
         IF(MOD(NPKUP,2).NE.0) EXPPUP=-1.
         A2UP=0.0
         A3UP=0.0
          DO 398 III=1,NDIMTO
          IIIOPP=2*KTRONC+1-(III-1)
          COEFAS=A59(IVTORP,III,IUP,NUP1,IISIG)
          COEFAAS=A59(IVTOR1,III,IOPPUP,NUP1,IISIG)
          IF(ISIG.EQ.0) THEN
            COEF4S=A59(IVTOR1,IIIOPP,IOPPUP,NUP1,IISIG)
          ELSE
            COEF4S=0.0
          ENDIF
          A2UP=A2UP+COEFAS*COEFAAS
          A3UP=A3UP+COEFAS*COEF4S
398       CONTINUE
         PERC2S=PERC2S+A2UP*EXPPUP*
     1   VECTPR(IVVUP,INKUP,IUP,IVTORP)*
     2   VECTPR(IVVUP,INKUP,IOPPUP,IVTOR1)
         PER23S=PER23S+A3UP*EXPPUP*
     1   VECTPR(IVVUP,INKUP,IUP,IVTORP)*
     2   VECTPR(IVVUP,INKUP,IOPPUP,IVTOR1)
C     PRINT*,'PER23S=',PER23S,'A3UP=',A3UP,'NUP=',NUP
400      CONTINUE
402     CONTINUE
404   CONTINUE

C      print 1337, perc2s, per23s, a2up, a3up
C
C  DELTA N=0   N'=N   Nlow=Nup
C
      IF(NUP.EQ.NLOW) THEN
C     SUM ON VT,VT',K,K'
C     EQUIVALENT OF THE DO 1 I=1,NDM(1) AND DO 2 J=1,NDM(2) OF
C     SUBROUTINE LNSTR (GODEFROID)
      DO 1004 IVTORP=1,9
       DO 1002 IVTOR1=1,9
        DO 1000 ILOW=1,2*NLOW+1
         DO 998 IUP=1,2*NUP+1
         KBASLO=ILOW-NLOW-1
         KBASUP=IUP-NUP-1
C        KBASLO=VALUE OF K FOR THE BASIS FUNCTION ON THE LEFT
C        SIDE OF THE MATRIX ELEMENT; TAKEN AS THE LOWER (INITIAL)
C        STATE OF THE TRANSITION
C        KBASUP= SAME FOR THE RIGHT SIDE OF THE MATRIX ELEMENT
C        KK2,RIGHT K FROM THE BASIS WAVEFUNCTION
C        KKP2= LEFT
         KK2=2*KBASUP
         KKP2=2*KBASLO
         ISIGN=+1
         IF(MOD(KNP2-KKP2,4).NE.0)ISIGN=-1
C
C        DELTA K=0    K'=K    Klow=Kup
C
         IF(IUP.EQ.ILOW) THEN
C         SUM ON THE TORSIONAL EIGENFUNCTIONS
          AN00A=0.
          AN00B=0.
          AN00C=0.
          DO 996 III=1,NDIMTO
          COEFA=A59(IVTORP,III,ILOW,NLOW1,IISIG)
          COEFA0=A59(IVTOR1,III,ILOW,NUP1,IISIG)
            IF(III.EQ.1) THEN
             COEFA0B=0.0
            ELSE
             COEFA0B=A59(IVTOR1,III-1,ILOW,NUP1,IISIG)
            ENDIF
            IF(III.EQ.NDIMTO) THEN
             COEFA0C=0.0
            ELSE
             COEFA0C=A59(IVTOR1,III+1,ILOW,NUP1,IISIG)
            ENDIF
          AN00A=AN00A + COEFA*COEFA0
          AN00B=AN00B + COEFA*COEFA0B
          AN00C=AN00C + COEFA*COEFA0C
996       CONTINUE
          CPROD0=VECTPR(IVVLO,INKLO,ILOW,IVTORP)*
     1    VECTPR(IVVUP,INKUP,ILOW,IVTOR1)
C
          IF(ABS(CPROD0).LT.EPS) GO TO 998
          TMUA=AN00A*AMUA0
          TMUA3=(AN00A-0.5*AN00B-0.5*AN00C)*AMUA3
          CT=CPROD0*AD3J(KNP2,2,KN2,-KKP2,0)*ISIGN
          CTT(1) = CTT(1) + CT*TMUA
          CTT(4) = CTT(4) + CT*TMUA3
          S=S+CT*(TMUA+TMUA3)
          S4=S4+ABS(CT*(TMUA+TMUA3))
C
C        DELTA K=+1   K=K'+1   Kup=Klow+1
C
         ELSEIF(IUP.EQ.(ILOW+1)) THEN
C
          AN01A=0.
          AN01B=0.
          AN01C=0.
          DO 1006 III=1,NDIMTO
            COEFA=A59(IVTORP,III,ILOW,NLOW1,IISIG)
            COEFA0=A59(IVTOR1,III,ILOW+1,NUP1,IISIG)
            IF(III.EQ.1) THEN
             COEFA0B=0.0
            ELSE
             COEFA0B=A59(IVTOR1,III-1,ILOW+1,NUP1,IISIG)
            ENDIF
            IF(III.EQ.NDIMTO) THEN
             COEFA0C=0.0
            ELSE
             COEFA0C=A59(IVTOR1,III+1,ILOW+1,NUP1,IISIG)
            ENDIF
            AN01A=AN01A + COEFA*COEFA0
            AN01B=AN01B + COEFA*COEFA0B
            AN01C=AN01C + COEFA*COEFA0C
1006      CONTINUE
          CPROD1=VECTPR(IVVLO,INKLO,ILOW,IVTORP)*
     1    VECTPR(IVVUP,INKUP,ILOW+1,IVTOR1)
          IF(ABS(CPROD1).LT.EPS) GO TO 998
          TMUB=AN01A*AMUB0
          TMUB3=(AN01A-0.5*(AN01B+AN01C))*AMUB3
          TMUC3=0.5*(-AN01B+AN01C)*AMUC3
          CT=CPROD1*AD3J(KNP2,2,KN2,-KKP2,-2)*ISIGN/SQRT(2.)
          CTT(2) = CTT(2) + CT*TMUB
          CTT(5) = CTT(5) + CT*TMUB3
          CTT(6) = CTT(6) + CT*TMUC3
          S=S+CT*(TMUB+TMUB3+TMUC3)
          S4=S4+ABS(CT*(TMUB+TMUB3+TMUC3))
C
C        DELTA K=-1   K=K'-1   Kup=Klow-1
C
         ELSEIF(IUP.EQ.(ILOW-1)) THEN
         AN0M1A=0.
         AN0M1B=0.
         AN0M1C=0.
          DO 1016 III=1,NDIMTO
            COEFA=A59(IVTORP,III,ILOW,NLOW1,IISIG)
            COEFA0=A59(IVTOR1,III,ILOW-1,NUP1,IISIG)
            IF(III.EQ.1) THEN
             COEFA0B=0.0
            ELSE
             COEFA0B=A59(IVTOR1,III-1,ILOW-1,NUP1,IISIG)
            ENDIF
            IF(III.EQ.NDIMTO) THEN
             COEFA0C=0.0
            ELSE
             COEFA0C=A59(IVTOR1,III+1,ILOW-1,NUP1,IISIG)
            ENDIF
            AN0M1A=AN0M1A + COEFA*COEFA0
            AN0M1B=AN0M1B + COEFA*COEFA0B
            AN0M1C=AN0M1C + COEFA*COEFA0C
1016      CONTINUE
          CPRODM1=VECTPR(IVVLO,INKLO,ILOW,IVTORP)*
     1    VECTPR(IVVUP,INKUP,ILOW-1,IVTOR1)
          IF(ABS(CPRODM1).LT.EPS) GO TO 998
          TMUB=-AN0M1A*AMUB0
          TMUB3=(-AN0M1A+0.5*(AN0M1B+AN0M1C))*AMUB3
          TMUC3=0.5*(-AN0M1B+AN0M1C)*AMUC3
          CT=CPRODM1*AD3J(KNP2,2,KN2,-KKP2,2)*ISIGN/SQRT(2.)
          CTT(2) = CTT(2) + CT*TMUB
          CTT(5) = CTT(5) + CT*TMUB3
          CTT(6) = CTT(6) + CT*TMUC3
          S=S+CT*(TMUB+TMUB3+TMUC3)
          S4=S4+ABS(CT*(TMUB+TMUB3+TMUC3))
         ENDIF
998      CONTINUE
1000    CONTINUE
1002   CONTINUE
1004  CONTINUE
C
C  DELTA N=+1     N=N'+1     Nup=Nlow+1
C
      ELSEIF(NUP.EQ.(NLOW+1)) THEN
      DO 1034 IVTORP=1,9
       DO 1032 IVTOR1=1,9
        DO 1030 ILOW=1,2*NLOW+1
         DO 1028 IUP=1,2*NUP+1
C     print*,' ivtorp = ',ivtorp,' ivtor1 = ',ivtor1
C     print*,' ilow = ',ilow,' iup = ',iup
         KBASLO=ILOW-NLOW-1
         KBASUP=IUP-NUP-1
         KKP2=2*KBASLO
         KK2=2*KBASUP
C     print*,' kbaslo = ',kbaslo,' kbasup = ',kbasup
         ISIGN=+1
         IF(MOD(KNP2-KKP2,4).NE.0)ISIGN=-1
C
C        DELTA K=0    K'=K    Klow=Kup
C
         IF(IUP.EQ.(ILOW+1)) THEN
         AN10A=0.
         AN10B=0.
         AN10C=0.
C        print*,'AN10A avant calc=',AN10A
          DO 1026 III=1,NDIMTO
          COEFA=A59(IVTORP,III,ILOW,NLOW1,IISIG)
          COEFA0=A59(IVTOR1,III,ILOW+1,NUP1,IISIG)
            IF(III.EQ.1) THEN
             COEFA0B=0.0
            ELSE
             COEFA0B=A59(IVTOR1,III-1,ILOW+1,NUP1,IISIG)
            ENDIF
            IF(III.EQ.NDIMTO) THEN
             COEFA0C=0.0
            ELSE
             COEFA0C=A59(IVTOR1,III+1,ILOW+1,NUP1,IISIG)
            ENDIF
          AN10A=AN10A + COEFA*COEFA0
          AN10B=AN10B + COEFA*COEFA0B
          AN10C=AN10C + COEFA*COEFA0C
1026      CONTINUE
          CPROD0=VECTPR(IVVLO,INKLO,ILOW,IVTORP)*
     1    VECTPR(IVVUP,INKUP,ILOW+1,IVTOR1)
C         print*,' cprod0 = ',cprod0
          IF(ABS(CPROD0).LT.EPS) GO TO 1028
          TMUA=AN10A*AMUA0
C         print*,'AN10A=',AN10A,'TMUA=',TMUA,'AMUA0=',AMUA0
          TMUA3=(AN10A-0.5*AN10B-0.5*AN10C)*AMUA3
          CT=CPROD0*AD3J(KNP2,2,KN2,-KKP2,0)*ISIGN
          CTT(1) = CTT(1) + CT*TMUA
          CTT(4) = CTT(4) + CT*TMUA3
          S=S+CT*(TMUA+TMUA3)
          S4=S4+ABS(CT*(TMUA+TMUA3))
C         PRINT*,'AD3J=',AD3J(KNP2,2,KN2,-KKP2,0),'KNP2=',KNP2,
C    1    'KN2=',KN2,'KKP2=',KKP2
C         print*,' in delta k = 0 : ct = ',ct,' s = ',s
C
C        DELTA K=+1   K=K'+1   Kup=Klow+1
C
         ELSEIF(IUP.EQ.(ILOW+2)) THEN
         AN11A=0.
         AN11B=0.
         AN11C=0.
          DO 1036 III=1,NDIMTO
          COEFA=A59(IVTORP,III,ILOW,NLOW1,IISIG)
          COEFA0=A59(IVTOR1,III,ILOW+2,NUP1,IISIG)
            IF(III.EQ.1) THEN
             COEFA0B=0.0
            ELSE
             COEFA0B=A59(IVTOR1,III-1,ILOW+2,NUP1,IISIG)
            ENDIF
            IF(III.EQ.NDIMTO) THEN
             COEFA0C=0.0
            ELSE
             COEFA0C=A59(IVTOR1,III+1,ILOW+2,NUP1,IISIG)
            ENDIF
          AN11A=AN11A + COEFA*COEFA0
          AN11B=AN11B + COEFA*COEFA0B
          AN11C=AN11C + COEFA*COEFA0C
1036      CONTINUE
          CPROD1=VECTPR(IVVLO,INKLO,ILOW,IVTORP)*
     1    VECTPR(IVVUP,INKUP,ILOW+2,IVTOR1)
C          print*,' cprod1 = ',cprod1
          IF(ABS(CPROD1).LT.EPS) GO TO 1028
          TMUB=AN11A*AMUB0
          TMUB3=(AN11A-0.5*(AN11B+AN11C))*AMUB3
          TMUC3=0.5*(-AN11B+AN11C)*AMUC3
          CT=CPROD1*AD3J(KNP2,2,KN2,-KKP2,-2)*ISIGN/SQRT(2.)
          CTT(2) = CTT(2) + CT*TMUB
          CTT(5) = CTT(5) + CT*TMUB3
          CTT(6) = CTT(6) + CT*TMUC3
          S=S+CT*(TMUB+TMUB3+TMUC3)
          S4=S4+ABS(CT*(TMUB+TMUB3+TMUC3))
C         PRINT*,'AD3J=',AD3J(KNP2,2,KN2,-KKP2,-2),'KNP2=',KNP2,
C    1    'KN2=',KN2,'KKP2=',KKP2
C         print*,' in delta k = +1 : ct = ',ct,' s = ',s
C
C        DELTA K=-1   K=K'-1   Kup=Klow-1
C
         ELSEIF(IUP.EQ.ILOW) THEN
C        PRINT*,'IUP=',IUP,'ILOW=',ILOW
         AN1M1A=0.
         AN1M1B=0.
         AN1M1C=0.
          DO 1046 III=1,NDIMTO
          COEFA=A59(IVTORP,III,ILOW,NLOW1,IISIG)
          COEFA0=A59(IVTOR1,III,ILOW,NUP1,IISIG)
C         PRINT*,'COEFA=',COEFA,'COEFA0=',COEFA0
            IF(III.EQ.1) THEN
             COEFA0B=0.0
            ELSE
             COEFA0B=A59(IVTOR1,III-1,ILOW,NUP1,IISIG)
            ENDIF
            IF(III.EQ.NDIMTO) THEN
             COEFA0C=0.0
            ELSE
             COEFA0C=A59(IVTOR1,III+1,ILOW,NUP1,IISIG)
            ENDIF
          AN1M1A=AN1M1A + COEFA*COEFA0
          AN1M1B=AN1M1B + COEFA*COEFA0B
          AN1M1C=AN1M1C + COEFA*COEFA0C
1046      CONTINUE
          CPRODM1=VECTPR(IVVLO,INKLO,ILOW,IVTORP)* 
     1    VECTPR(IVVUP,INKUP,ILOW,IVTOR1)
C         print*,' cprodm1 = ',cprodm1
          IF(ABS(CPRODM1).LT.EPS) GO TO 1028
          TMUB=-AN1M1A*AMUB0
          TMUB3=(-AN1M1A+0.5*(AN1M1B+AN1M1C))*AMUB3
          TMUC3=0.5*(-AN1M1B+AN1M1C)*AMUC3
          CT=CPRODM1*AD3J(KNP2,2,KN2,-KKP2,2)*ISIGN/SQRT(2.)
C         PRINT*,'AD3J=',AD3J(KNP2,2,KN2,-KKP2,2),'KNP2=',KNP2,
C    1    'KN2=',KN2,'KKP2=',KKP2
          CTT(2) = CTT(2) + CT*TMUB
          CTT(5) = CTT(5) + CT*TMUB3
          CTT(6) = CTT(6) + CT*TMUC3
          S=S+CT*(TMUB+TMUB3+TMUC3)
          S4=S4+ABS(CT*(TMUB+TMUB3+TMUC3))
C         PRINT*,'AD3J=',AD3J(KNP2,2,KN2,-KKP2,2),'KNP2=',KNP2,
C    1    'KN2=',KN2,'KKP2=',KKP2
C         print*,' in delta k = -1 : ct = ',ct,' s = ',s
         ENDIF
1028     CONTINUE
1030    CONTINUE
1032   CONTINUE
1034  CONTINUE
C
C  DELTA N=-1
C
      ELSEIF(NUP.EQ.(NLOW-1)) THEN
      DO 1064 IVTORP=1,9
       DO 1062 IVTOR1=1,9
        DO 1060 ILOW=1,2*NLOW+1
         DO 1058 IUP=1,2*NUP+1
         KBASLO=ILOW-NLOW-1
         KBASUP=IUP-NUP-1
         KKP2=2*KBASLO
         KK2=2*KBASUP
         ISIGN=+1
         IF(MOD(KNP2-KKP2,4).NE.0)ISIGN=-1
C
C        DELTA K=0    K'=K    Klow=Kup
C
         IF(IUP.EQ.(ILOW-1)) THEN
          ANM10A=0.                         
          ANM10B=0.                    
          ANM10C=0.                     
          DO 1056 III=1,NDIMTO
          COEFA=A59(IVTORP,III,ILOW,NLOW1,IISIG)
          COEFA0=A59(IVTOR1,III,ILOW-1,NUP1,IISIG)
            IF(III.EQ.1) THEN
             COEFA0B=0.0
            ELSE
             COEFA0B=A59(IVTOR1,III-1,ILOW-1,NUP1,IISIG)
            ENDIF
            IF(III.EQ.NDIMTO)THEN
             COEFA0C=0.0
            ELSE
             COEFA0C=A59(IVTOR1,III+1,ILOW-1,NUP1,IISIG)
            ENDIF
          ANM10A=ANM10A + COEFA*COEFA0
          ANM10B=ANM10B + COEFA*COEFA0B
          ANM10C=ANM10C + COEFA*COEFA0C
1056      CONTINUE
          CPROD0=VECTPR(IVVLO,INKLO,ILOW,IVTORP)*
     1    VECTPR(IVVUP,INKUP,ILOW-1,IVTOR1)
          IF(ABS(CPROD0).LT.EPS) GO TO 1058
          TMUA=ANM10A*AMUA0
          TMUA3=(ANM10A-0.5*ANM10B-0.5*ANM10C)*AMUA3
          CT=CPROD0*AD3J(KNP2,2,KN2,-KKP2,0)*ISIGN
          CTT(1) = CTT(1) + CT*TMUA
          CTT(4) = CTT(4) + CT*TMUA3
          S=S+CT*(TMUA+TMUA3)
          S4=S4+ABS(CT*(TMUA+TMUA3))
C
C        DELTA K=+1   K=K'+1   Kup=Klow+1
C
         ELSEIF(IUP.EQ.ILOW) THEN
          ANM11A=0.                            
          ANM11B=0.                   
          ANM11C=0.                    
          DO 1066 III=1,NDIMTO
          COEFA=A59(IVTORP,III,ILOW,NLOW1,IISIG)
          COEFA0=A59(IVTOR1,III,ILOW,NUP1,IISIG)
          IF(III.EQ.1)THEN
             COEFA0B=0.0
            ELSE
             COEFA0B=A59(IVTOR1,III-1,ILOW,NUP1,IISIG)
            ENDIF
          IF(III.EQ.NDIMTO)THEN
             COEFA0C=0.0
            ELSE
             COEFA0C=A59(IVTOR1,III+1,ILOW,NUP1,IISIG)
            ENDIF
          ANM11A=ANM11A + COEFA*COEFA0
          ANM11B=ANM11B+ COEFA*COEFA0B
          ANM11C=ANM11C + COEFA*COEFA0C
1066      CONTINUE
          CPROD1=VECTPR(IVVLO,INKLO,ILOW,IVTORP)*
     1    VECTPR(IVVUP,INKUP,ILOW,IVTOR1)
          IF(ABS(CPROD1).LT.EPS) GO TO 1058
          TMUB=ANM11A*AMUB0
          TMUB3=(ANM11A-0.5*(ANM11B+ANM11C))*AMUB3
          TMUC3=0.5*(-ANM11B+ANM11C)*AMUC3
          CT=CPROD1*AD3J(KNP2,2,KN2,-KKP2,-2)*ISIGN/SQRT(2.)
          CTT(2) = CTT(2) + CT*TMUB
          CTT(5) = CTT(5) + CT*TMUB3
          CTT(6) = CTT(6) + CT*TMUC3
          S=S+CT*(TMUB+TMUB3+TMUC3)
          S4=S4+ABS(CT*(TMUB+TMUB3+TMUC3))
C
C        DELTA K=-1   K=K'-1   Kup=Klow-1
C
         ELSEIF(IUP.EQ.(ILOW-2)) THEN
          ANM1MA=0.                   
          ANM1MB=0.                    
          ANM1MC=0.                    
          DO 1076 III=1,NDIMTO
          COEFA=A59(IVTORP,III,ILOW,NLOW1,IISIG)
          IF(ILOW.EQ.1.OR.ILOW.EQ.2) THEN
          COEFA0=0.0
          ELSE
          COEFA0=A59(IVTOR1,III,ILOW-2,NUP1,IISIG)
          ENDIF
            IF(III.EQ.1)THEN
             COEFA0B=0.0
            ELSE
             COEFA0B=A59(IVTOR1,III-1,ILOW-2,NUP1,IISIG)
            ENDIF
            IF(III.EQ.NDIMTO)THEN
             COEFA0C=0.0
            ELSE
             COEFA0C=A59(IVTOR1,III+1,ILOW-2,NUP1,IISIG)
            ENDIF
          ANM1MA=ANM1MA + COEFA*COEFA0
          ANM1MB=ANM1MB + COEFA*COEFA0B
          ANM1MC=ANM1MC + COEFA*COEFA0C
1076      CONTINUE
          CPRODM1=VECTPR(IVVLO,INKLO,ILOW,IVTORP)*
     1    VECTPR(IVVUP,INKUP,ILOW-2,IVTOR1)
          IF(ABS(CPRODM1).LT.EPS) GO TO 1058
          TMUB=-ANM1MA*AMUB0
          TMUB3=(-ANM1MA+0.5*(ANM1MB+ANM1MC))*AMUB3
          TMUC3=0.5*(-ANM1MB+ANM1MC)*AMUC3
          CT=CPRODM1*AD3J(KNP2,2,KN2,-KKP2,2)*ISIGN/SQRT(2.)
          CTT(2) = CTT(2) + CT*TMUB
          CTT(5) = CTT(5) + CT*TMUB3
          CTT(6) = CTT(6) + CT*TMUC3
          S=S+CT*(TMUB+TMUB3+TMUC3)
          S4=S4+ABS(CT*(TMUB+TMUB3+TMUC3))
         ENDIF
1058     CONTINUE
1060    CONTINUE
1062   CONTINUE
1064  CONTINUE
      ENDIF
C
      IF(ABS(CTT(1)).GT.EPS) MUABC(1:1)='A'
      IF(ABS(CTT(2)).GT.EPS) MUABC(2:2)='B'
      IF(ABS(CTT(4)).GT.EPS) MUABC(4:4)='a'
      IF(ABS(CTT(5)).GT.EPS) MUABC(5:5)='b'
      IF(ABS(CTT(6)).GT.EPS) MUABC(6:6)='c'
C
      S2 = 0.0
      S3 = 0.0
      DO 2 KK = 1,(6)
         S2 = S2 + CTT(KK)
         S3 = S3 + ABS(CTT(KK))
    2 CONTINUE
      IF(ABS(S-S2).GT.EPS) stop
C
      S=S**2
C
C If the line strength is larger than eps, then...
C
      
      IF(S.LE.EPS) GO TO 3001
C 
C ...calculate the percentages of contributions A,B,C,a,b,c to S^{1/2}
C    and keep the signs
C
      DO 3 KK = 1,(6)
         ICTT(KK) = NINT(CTT(KK)*100.0/S3)
    3 CONTINUE
C
C ...evaluate the cancellation factor (Cowan, formula 14.107)
C
      CF=(S2/S4)**2
C
C add the degeneracy and temperature dependence factors
C
 3001 CONTINUE
      SSQ = S*Q
C      BLTZ2=1.
      RINT=BLTZ2*S*Q
      RETURN
      END
