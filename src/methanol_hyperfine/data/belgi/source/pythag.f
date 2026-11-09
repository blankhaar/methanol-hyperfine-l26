********************************************                                  

      double precision function pythag(a,b)
      implicit none 
      save
      double precision a,b
      
c
c     finds dsqrt(a**2+b**2) without overflow or destructive underflow
c
      double precision p,r,s,t,u
      p = dmax1(dabs(a),dabs(b))

      if (p .eq. 0.0d0) then
         pythag = p
         return
      endif

      r = (dmin1(dabs(a),dabs(b))/p)**2
      t = 0.0d0

      do while (t .ne. 4.0d0)
         t = 4.0d0 + r
         s = r/t
         u = 1.0d0 + 2.0d0*s
         p = u*p
         r = (s/u)**2 * r
         t = 4.0d0 + r
      end do
      pythag = p
      return
      end


