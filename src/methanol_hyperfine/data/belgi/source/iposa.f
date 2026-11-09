c**********************************************************
c     iposa for methanol by li-hong xu, march 25, 1994                         
c**********************************************************
      function iposa(egvc,egvl,a,egvc4)                             
      implicit none
      save
c***********************************************************
c     
c     determines the position of eigenvectors and values associated             
c     with the k data and for isig=0                                            
c     

      include 'common_dim.f'
      include 'common_quant.f'
      include 'common_chara.f'
      include 'common_rotor.f'
      include 'common_c2c.f'
      include 'common_eigen.f'
      
      integer iposa, i, j, m
      double precision egvc, sum, egvl, egvc4, a
      integer nvt, kkk, ipr, ii, ikt, kt, nin, jj, mm, niv, iprity
      double precision xmax

      dimension egvc(ndmx,ndmx),sum(9,ndmx),nvt(ndmx),kkk(ndmx),                
     &     egvl(ndmx)

      dimension a(9,ntormx,kdmx)
      dimension egvc4(kdmx,ntormx,ndmx)
      
      iposa = 0

c     parity assignment                                                            
      if (pr.eq.'+') ipr=1                                                      
      if (pr.eq.'-') ipr=2                                                      
c     torsional labeling                                                           
      do i=1,ndmx
         do j=1,9
            sum(j,i)=0.d0
         end do
      end do

      do i=1,ndmx
         nvt(i)=0
         kkk(i)=0
      end do

c     loop through eigenvectors
      do ii=1,9*(2*n+1)                                                       
         ikt=0                 
c     loop through vt blocks
         do kt=1,9*(2*n+1),2*n+1                                               
            if (kt.gt.9*(2*n+1)) cycle
            ikt=ikt+1                                                             
            do i=kt,2*n+kt                                                      
c     print *, 'egvc',egvc(i,ii)                                          
c     magnitude of the projection of the vector onto this
c     vt block
               sum(ikt,ii)=egvc(i,ii)*egvc(i,ii)+sum(ikt,ii)                       
            end do
         end do
      end do

c     loop through vt blocks
      do i=1,9*(2*n+1)
         xmax=0.d0
c     loop through vt
         do j=1,9
c     find vector with largest vt component
            if (sum(j,i).gt.xmax) then
               xmax=sum(j,i)
               nin=j
            end if
         end do
c     label it
         nvt(i)=nin-1
      end do

c     rotational labeling                                                          
c     loop through eigenvectors
      do jj=1,9*(2*n+1)                                                       
         xmax=0.d0
         kkk(jj)=n
c     loop through the vt block containing jj
         do ii=(2*n+1)*nvt(jj)+1,(2*n+1)*nvt(jj)+(n+1)
            m=ii-(n+1)-nvt(jj)*(2*n+1)
            mm=ii-2*m
            if (ii.ne.mm) then 
               if ((egvc(ii,jj)**2+egvc(mm,jj)**2).gt.xmax) then
                  xmax=egvc(ii,jj)**2+egvc(mm,jj)**2
                  kkk(jj)=iabs(m)
               end if
            else
               if (egvc(ii,jj)**2.gt.xmax) then
                  xmax=egvc(ii,jj)**2
                  kkk(jj)=0
               end if
            end if
         end do
      end do
      
c     loop through eigenvectors
      do ii=1,9*(2*n+1)                                                       
         if (kk.eq.kkk(ii).and.iv.eq.nvt(ii).and.isig.eq.0)then                  
            call parity(a,egvc4,ii)
c     print*,'ii apres parity',ii
            niv=n+nvt(ii)
            if(mod(niv,2).eq.0) then
               if(per23l.ge.0.)then
                  iprity=1
               else
                  iprity=2
               endif
            else  
               if(per23l.ge.0.)then
                  iprity=2
               else
                  iprity=1
               endif
            endif

            if(ipr.eq.iprity) then                                                 
               iposa=ii                                                              

            end if                                                                  
         end if                                                                  
      end do

      return                                                                    
      end                                                                       
      
