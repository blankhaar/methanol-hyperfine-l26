      program belgi_wavefunctions
C     Project-local exporter, derived from methanol_lamda's frozen driver.
C     Uses the unchanged accepted BELGI Hamiltonian and its state labels.
      implicit none
      include 'common_dim.f'
      include 'common_cst.f'
      include 'common_dat.f'
      include 'common_par.f'
      include 'common_expdat.f'
      include 'common_big2.f'
      include 'common_cntl.f'
      include 'common_eigen.f'
      include 'common_tor59.f'
      integer output_jmax, ii, ios, j, k, vt, ipr, ink
      integer isig, ividx, nk, ik, it, nt, ib, nreq, ir
      integer maxreq
      parameter (maxreq=6*(46+1)**2)
      integer rj(maxreq), rk(maxreq), rv(maxreq), rp(maxreq),
     &        rs(maxreq)
      character*512 arg, request_path
      character*64 vector_header
      double precision phi, energy
      if(command_argument_count().ne.3) stop 2
      call get_command_argument(1,arg)
      read(arg,*) output_jmax
      call get_command_argument(2,arg)
      read(arg,*) ktronc
      if(output_jmax.gt.46.or.ktronc.gt.20) stop 2
      call get_command_argument(3,request_path)
      open(17,file=trim(request_path),status='old')
      nreq=0
 10   continue
      read(17,*,iostat=ios) isig,vt,j,k,ipr
      if(ios.ne.0) goto 20
      nreq=nreq+1
      if(nreq.gt.maxreq) stop 2
      rs(nreq)=isig
      rv(nreq)=vt
      rj(nreq)=j
      rk(nreq)=k
      rp(nreq)=ipr
      goto 10
 20   close(17)
      read(5,'(A)',iostat=ios) vector_header
      if(vector_header.ne.'BELGI_FROZEN_CST_205') stop 3
      do ii=1,cstm
         read(5,*,iostat=ios) cst(1,ii)
         if(ios.ne.0) stop 3
         cst(2,ii)=cst(1,ii)
      enddo
      np(1)=0
      np(2)=0
      ito=0
      ktron1=ktronc+1
      conv=1.d0
      irmw=1
      jmax=output_jmax
      vtmax=2
      int=1
      unc=0
      ibug1=0
      ibug2=0
      ibug3=0
      ibug4=0
      ibug5=0
      ibgtme=0
      ibgvtd=0
      evib(1)=0.d0
      evib(2)=0.d0
      term=0.d0
      en0=0.d0
      label=1
      call setup(cst)
      write(6,'(A)') 'BEGIN_WAVEFUNCTIONS'
      do ir=1,nreq
         isig=rs(ir)
         vt=rv(ir)
         j=rj(ir)
         k=rk(ir)
         ipr=rp(ir)
         ink=j*j+j+k+1
         ividx=vt+6
         if(isig.eq.0) then
            ink=ink-ipr*j
            ividx=vt+1
         endif
         energy=ener(ividx,ink,1)
         write(6,100) ir,isig,vt,j,k,ipr,energy
         do ik=1,2*j+1
            do it=1,2*ktronc+1
               phi=0.d0
               do ib=1,9
                  phi=phi+vectpr(ividx,ink,ik,ib)*
     &                 a59(ib,it,ik,j+1,isig+1)
               enddo
               write(6,101) ik-j-1,it-ktronc-1,phi
            enddo
         enddo
      enddo
      write(6,'(A)') 'END_WAVEFUNCTIONS'
 100  format('STATE ',6(I0,1X),ES25.17E3)
 101  format(I0,1X,I0,1X,ES25.17E3)
      end
