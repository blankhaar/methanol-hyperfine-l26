      integer lwork, liwork, iwork
      parameter (lwork=2*ndmx*ndmx+6*ndmx+1)
      parameter (liwork=3+5*ndmx+1)
C      parameter (lwork=ndmx)
      double precision work, egvc3, h, vectpr
      double precision egvc3j, egvlj
      integer nbr
      COMMON/EIGEN/H(NDMX,NDMX),WORK(lwork),EGVC3(KDMX,9),NBR(2),
     *     VECTPR(10, NPNP, KDMX, 9),
     *     EGVC3J(5, 2, kdmx, NDMX),
     *     EGVLJ(5, 2, kdmx, 3),
     *     IWORK(liwork)
      
