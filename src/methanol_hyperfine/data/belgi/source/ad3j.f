
      FUNCTION AD3J(J1,J2,J3,M1,M2)
      implicit double precision (a-h, o-z)
      save

      include 'common_facts.f'
      DIMENSION MC(12)
      AD3J = 0.
      IF(IABS(J1-J2).GT.J3.OR.(J1+J2).LT.J3.OR.IABS(J2-J3).GT.J1.OR.(J2+
     1J3).LT.J1.OR.IABS(J3-J1).GT.J2.OR.(J3+J1).LT.J2) GO TO 99
      IF(IABS(M1).GT.J1.OR.IABS(M2).GT.J2.OR.IABS(M1+M2).GT.J3) GO TO 99
      IF  ((IABS(M1)+IABS(M2)).EQ.0.AND.MOD(J1+J2+J3,4).NE.0) GO TO 99
      MC(1) = J1 + J2 - J3
      MC(2) = J1 - J2 + J3
      MC(3) = -J1 + J2 + J3
      MC(4) = J1 + M1
      MC(5) = J1 - M1
      MC(6) = J2 + M2
      MC(7) = J2 - M2
      MC(8) = J3 - M1 - M2
      MC(9) = J3 + M1 + M2
      MC(10) = J1 + J2 + J3 + 2
      MC(11) = -J3 + J2 - M1
      MC(12) = -J3 + J1 + M2
      DO 6 I =1,12
      IF (MOD (MC(I),2).NE.0) GO TO 98
      MC(I) = MC(I)/2
6     CONTINUE
      IZD = MAX0(MC(11),MC(12),0)
      IZM = MIN0(MC(1),MC(5),MC(6))
      ISIG = 1
      IF (MOD(IZD,2).NE.0) ISIG = -ISIG
      FCT = GAM(MC(1)+1) + GAM(MC(2)+1) + GAM(MC(3)+1)
     1    + GAM(MC(5)+1) + GAM(MC(4)+1) + GAM(MC(7)+1)
     2    + GAM(MC(6)+1) + GAM(MC(8)+1) + GAM(MC(9)+1)
     3    - GAM(MC(10)+1)
      FCT = FCT/2.
      AD3J = 0.0E0
      DO 60 IZ = IZD,IZM
      ADEN = GAM(IZ+1) + GAM(MC(1)-IZ+1) + GAM(MC(5)-IZ+1)
     1     + GAM(MC(6)-IZ+1) + GAM(-MC(11)+IZ+1)
     2     + GAM(-MC(12)+IZ+1)
      ADEN = dble(ISIG)/EXP(ADEN - FCT)
      AD3J = AD3J + ADEN
      ISIG = -ISIG
60    CONTINUE
      IF (MOD(J1-J2+M1+M2,4) .NE. 0) AD3J = -AD3J
      GO TO 99
98    WRITE(6,101) J1,J2,J3,M1,M2
101   FORMAT(60H ERROR IN ATTEMPTING TO CALCULATE 3-J SYMBOL WITH J - VA
     1LUES,3I3,9H M-VALUES,2I3)
99    RETURN
      END
