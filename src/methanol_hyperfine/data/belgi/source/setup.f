c ************************************************************
      subroutine setup(cst)                                                     
      implicit none
      save

c ************************************************************
c     this subroutine sets up the deriv matrix with all the                     
c     energies and derivatives.                                                 

      include 'common_dim.f'
      include 'common_cst.f'
      include 'common_eigen.f'
      include 'common_eigen1.f'
      include 'common_dat.f'
      include 'common_tor.f'
      include 'common_rotor.f'
      include 'common_expdat.f'
      include 'common_par.f'
      include 'common_big2.f'
      include 'common_symb.f'
      include 'common_quant.f'
      include 'common_chara.f'

      double precision v1, v2, e
      dimension v1(maxparam),v2(maxparam)

      integer jjmax, np12, jj, lo, loo, ink, ni, ipr
      integer i, j, m, in, inn, max
      
      integer ivtorp, iu, ii

      max = 2*nn2+1
c                                                                               
c     initialization                                                            
      jjmax=(jmax+1)**2                                                         
      np12=np(1)+np(2)
      do i=1,10                                                              
         do j=1,jjmax                                                       
            do jj=1,np12                                                    
               ener(i,j,1)=0.d0                                                 
               ener(i,j,2)=0.d0                                                 
               derivl(i,j,jj)=0.d0                                              
               derivs(i,j,jj)=0.0                                              
            end do

            do in = 1, max
               do inn=1,9
                  vectpr(i, j, in, inn) = 0.0
               enddo
            enddo

         end do
      end do

C     test
C     call vset2
C     test

      if(irmw.eq.2) then
         do lo=1,2
            if(lo.eq.1) then
               ni=1
               m=1
            else
               ni=ndmx+1
               m=kdmx+1
            endif
            loo=lo
            do isig=0,1                                                            
               iisig = isig+1

               if(isig.eq.1) then                                                     
                  do  n=0,jmax                                                     
                     call rottor(cst,lo,egvl(ni),egvc(1,ni),etor(m,1),                
     1                    a(1,1,m))                                                   
                     nrotor=(2*n+1)*9                                                 
                     do iv=0,4                                                     
                        do kk=-n,n                                                 
                           call encal(e,loo,egvl(ni),egvc(1,ni),
     2                          a(1,1,m))             
                           if(lo.eq.1) call vset(1,v1,cst)
                           if(lo.eq.2) call vset(2,v2,cst)
                           ink=n**2+n+kk+1                                            
                           ener(iv+6,ink,loo)=e                                       
                           
                           IF(INT.EQ.1 .AND. LO.EQ.1) THEN
                              IVTORP=1
                              DO I=1,NROTOR
                                 IU=(2*N+1)*IVTORP+1
                                 IF(I.EQ.IU) IVTORP=IVTORP+1
                                 II=I-(2*N+1)*(IVTORP-1)
                                 VECTPR(IV+6,INK,II,IVTORP)=
     &                                EGVC3(II,IVTORP)
                              enddo
                           ENDIF

                           do i=1,np12                                             
                              if(loo.eq.1) then
                                 derivl(iv+6,ink,i)=v1(i) 
                 if(iv.eq.0) then
c            write(6,*),'N=',N,'KK=',KK,'isig=',isig,'derivl=',
c     1    derivl(iv+6,ink,i),'i=',i
                    endif                           
                              else
                                 derivs(iv+6,ink,i)=v2(i)
                  if(iv.eq.0) then
c            write(6,*)'N=',N,'KK=',KK,'isig=',isig,'derivs=',
c     1    derivs(iv+6,ink,i),'i=',i
                    endif                                     
                              endif
                           end do
                        end do
                     end do
                  end do
               else                                                                   
                  do n=0,jmax                                                      
                     call rottor(cst,lo,egvl(ni),egvc(1,ni),etor(m,1)                 
     1                    ,a(1,1,m))                                                  
                     nrotor=(2*n+1)*9                                                 
                     do iv=0,4                                                     
                        do kk=0,n                                                  
                           if(kk.eq.0) then                                           
                              pr='+'                                                  
                              ipr=1                                                   
                              call encal(e,loo,egvl(ni),egvc(1,ni), 
     3                             a(1,1,m))          
                              if(lo.eq.1) call vset(1,v1,cst)
                              if(lo.eq.2) call vset(2,v2,cst)
                              ink=n**2+n+kk+1-ipr*n                                   
                              ener(iv+1,ink,loo)=e                                    
                           
                              IF(INT.EQ.1 .AND. LO.EQ.1) THEN
                                 IVTORP=1
                                 DO I=1,NROTOR
                                    IU=(2*N+1)*IVTORP+1
                                    IF(I.EQ.IU) IVTORP=IVTORP+1
                                    II=I-(2*N+1)*(IVTORP-1)
                                    VECTPR(IV+1,INK,II,IVTORP)=
     &                                   EGVC3(II,IVTORP)
                                 enddo
                              ENDIF

                              do i=1,np12                                          
                                 if(lo.eq.1) then
                                    derivl(iv+1,ink,i)=v1(i) 
                     if(iv.eq.0) then
c            write(6,*)'N=',N,'KK=',KK,'isig=',isig,'derivl=',
c     1    derivl(iv+1,ink,i),'i=',i
                    endif                                  
                                 else
                                    derivs(iv+1,ink,i)=v2(i)
              if(iv.eq.0) then
c           write(6,*)'N=',N,'KK=',KK,'isig=',isig,'derivs=',
c     1    derivs(iv+1,ink,i),'i=',i
                    endif                                   
                                 endif
                              end do
                           else                                                       
                              do ipr=0,1                                           
                                 if(ipr.eq.0) pr='-'                                  
                                 if(ipr.eq.1) pr='+'                                  
                                 ink=n**2+n+kk+1-ipr*n                                
                                 call encal(e,loo,egvl(ni),egvc(1,ni),                
     1                                a(1,1,m))
                                 if(lo.eq.1) call vset(1,v1,cst)
                                 if(lo.eq.2) call vset(2,v2,cst)
                                 ener(iv+1,ink,loo)=e                                 

                                 IF(INT.EQ.1 .AND. LO.EQ.1) THEN
                                    IVTORP=1
                                    DO I=1,NROTOR
                                       IU=(2*N+1)*IVTORP+1
                                       IF(I.EQ.IU) IVTORP=IVTORP+1
                                       II=I-(2*N+1)*(IVTORP-1)
                                       VECTPR(IV+1,INK,II,IVTORP)=
     &                                      EGVC3(II,IVTORP)
                                    enddo
                                 ENDIF

                                 do i=1,np12                                       
                                    if(lo.eq.1) then
                                       derivl(iv+1,ink,i)=v1(i)
                     if(iv.eq.0) then
c             write(6,*)'N=',N,'KK=',KK,'isig=',isig,'derivl=',
c     1    derivl(iv+1,ink,i),'i=',i
                    endif                                
                                    else 
                                       derivs(iv+1,ink,i)=v2(i)  
               if(iv.eq.0) then
c            write(6,*)'N=',N,'KK=',KK,'isig=',isig,'derivs=',
c     1    derivs(iv+1,ink,i),'i=',i
                    endif                              
                                    endif
                                 end do
                              end do
                           endif                                                      
                        end do
                     end do
                  end do
               endif                                                                  
            end do
         end do
      elseif(irmw.eq.1) then
         ni=1
         m=1
         lo=1
         loo=lo
         do  isig=0,1                                                           
            iisig = isig+1
            if(isig.eq.1) then                                                     
               do n=0,jmax                                                     
                  call rottor(cst,lo,egvl(ni),egvc(1,ni),etor(m,1),                
     1                 a(1,1,m))                                                   
                  nrotor=(2*n+1)*9                                                 
                  do iv=0,4                                                    
                     do kk=-n,n                                                
                        call encal(e,loo,egvl(ni),egvc(1,ni),a(1,1,m))             
                        call vset(1,v1,cst)
                        ink=n**2+n+kk+1                                            
                        ener(iv+6,ink,loo)=e                                       

                        IF(INT.EQ.1 .AND. LO.EQ.1) THEN
                           IVTORP=1
                           DO I=1,NROTOR
                              IU=(2*N+1)*IVTORP+1
                              IF(I.EQ.IU) IVTORP=IVTORP+1
                              II=I-(2*N+1)*(IVTORP-1)
                              VECTPR(IV+6,INK,II,IVTORP)=
     &                             EGVC3(II,IVTORP)
                           enddo
                        ENDIF

                        do i=1,np12                                            
                           derivl(iv+6,ink,i)=v1(i)                                
                        end do
                     end do
                  end do
               end do
            else                                                                   
               do n=0,jmax                                                     
                  call rottor(cst,lo,egvl(ni),egvc(1,ni),etor(m,1)                 
     1                 ,a(1,1,m))                                                  
                  nrotor=(2*n+1)*9                                                 
                  do iv=0,4                                                    
                     do kk=0,n                                                 
                        if(kk.eq.0) then                                           
                           pr='+'                                                  
                           ipr=1                                                   
                           call encal(e,loo,egvl(ni),egvc(1,ni),
     3                          a(1,1,m))          
                           call vset(1,v1,cst)
                           ink=n**2+n+kk+1-ipr*n                                  
                           ener(iv+1,ink,loo)=e                                    

                           IF(INT.EQ.1 .AND. LO.EQ.1) THEN
                              IVTORP=1
                              DO I=1,NROTOR
                                 IU=(2*N+1)*IVTORP+1
                                 IF(I.EQ.IU) IVTORP=IVTORP+1
                                 II=I-(2*N+1)*(IVTORP-1)
                                 VECTPR(IV+1,INK,II,IVTORP)=
     &                                EGVC3(II,IVTORP)
                              enddo
                           ENDIF

                           do i=1,np12                                         
                              derivl(iv+1,ink,i)=v1(i)                             
                           end do
                        else                                                       
                           do ipr=0,1                                          
                              if(ipr.eq.0) pr='-'                                  
                              if(ipr.eq.1) pr='+'                                  
                              ink=n**2+n+kk+1-ipr*n                                
                              call encal(e,loo,egvl(ni),egvc(1,ni),               
     1                             a(1,1,m))
                              call vset(1,v1,cst)
                              ener(iv+1,ink,loo)=e                                 
                              
                              IF(INT.EQ.1 .AND. LO.EQ.1) THEN
                                 IVTORP=1
                                 DO I=1,NROTOR
                                    IU=(2*N+1)*IVTORP+1
                                    IF(I.EQ.IU) IVTORP=IVTORP+1
                                    II=I-(2*N+1)*(IVTORP-1)
                                    VECTPR(IV+1,INK,II,IVTORP)=
     &                                   EGVC3(II,IVTORP)
                                 enddo
                              ENDIF

                              do  i=1,np12                                      
                                 derivl(iv+1,ink,i)=v1(i)                          
                              end do
                           end do                                                
                        endif                                                      
                     end do
                  end do
               end do
            endif                                                                  
         end do
      endif

C      if ( int .eq. 1 ) then
C         call setup59(cst)
C         return
C      endif

      return                                                                    
      end                                                                       
c     here begins a test for the case of a important mixing of the              
c     e states between +k and -k basis functions                                
c     i2=i-2*kk                                                                 
c     ratio=egvc(i,ipose)/egvc(i2,ipose)                                        
c     isgn=' '                                                                  
c     if(abs(ratio).gt.eps) then                                                
c     isgn='+'                                                                  
c     if(ratio.lt.0.d0) isgn='-'                                                
c     endif                                                                     
c     return                                                                    
c     end                                                                       
c     *************************************************************************  

