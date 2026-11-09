
C  *********************************************************************
C
C
C       CALCULATE THE LOGS OF FACTORIALS REQUIRED BY THE RACAH
C       COEFFICIENT ROUTINE RAC ( N.S. SCOTT)
C
C  ****************************************************************



      SUBROUTINE FACTT
      implicit double precision (a-h, o-z)
      save
      
      include 'common_facts.f'
      include 'common_consts.f'

      DATA FIFTY/50.0E0/
      GAM(1) = ONE
      GAM(2) = ONE
      X = TWO
      DO 10 I = 3,50
      GAM(I) = GAM(I-1)*X
      X = X + ONE
10    CONTINUE
      DO 20 I = 1,50
      GAM(I) = LOG(GAM(I))
20    CONTINUE
      X = FIFTY
      DO 30 I = 51,500
      GAM(I) = GAM(I-1) + LOG(X)
      X = X + ONE
30    CONTINUE
      RETURN
      END
