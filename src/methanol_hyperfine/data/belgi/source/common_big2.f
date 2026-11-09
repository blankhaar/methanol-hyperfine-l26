      double precision derivl, ener, derivs
      integer jsq
      double precision ener2, deriv
      parameter (jsq = (46+1)**2)

      common/big2/derivl(10,jsq,maxparam),ener(10,jsq,2),
     *     derivs(10,jsq,maxparam),ener2(10,jsq),
     *     deriv(10, jsq, maxparam)
