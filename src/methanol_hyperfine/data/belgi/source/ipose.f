C     *************************************************************************
      integer function ipose(egvc,egvl)                                         
C     *************************************************************************
C     Determines the position of the eigenvalue associated the k data and for
C     the state lo in the case isig = 1.
C     *************************************************************************
      implicit none
      save

      include 'common_dim.f'
      include 'common_quant.f'
      include 'common_chara.f'
      include 'common_rotor.f'

      double precision egvc, egvl, sum, xmax
      integer nvt, kkk, ii, ikt, kt, jj, nin, i, j
      
      dimension egvc(ndmx,ndmx),nvt(iposedim),kkk(iposedim),egvl(ndmx)
      dimension sum(iposedim,9)                                                     
      
      ipose = 0
      
c     torsional labeling                                                           
      do i=1,iposedim                                                              
         do j=1,9                                                              
            sum(i,j)=0.d0                                                           
         end do
      end do
      
      do i=1,iposedim                                                              
         nvt(i)=0                                                                
         kkk(i)=0                                                                
      end do

C     loop through the 9*(2*n+1) eigenvectors
      do ii=1,9*(2*n+1)                                                       
         ikt=0                                                                     
         kt = 1

C     loop through the 9 vt blocks
         do kt=1,9*(2*n+1),2*n+1                                               
            ikt=ikt+1             

C     block begins at index kt and ends at 2*n+kt
            do i=kt,2*n+kt                                                      
C     sum contains the portion of the modulus that is just
C     in this torsional subblock
               sum(ii,ikt)=egvc(i,ii)*egvc(i,ii)+sum(ii,ikt)                       
            end do
         end do
      end do
c     print *, nvt(ii)                                                        
      
C     loop through vt
      do j=1,9                                                               
C     loop through K
         do jj=1,2*n+1                                                        
            xmax=0.d0                                                             
C     find the vector with the largest component for this vt
C     store it in nin
            do i=1,9*(2*n+1)                                                   
               if (sum(i,j).gt.xmax) then                                          
                  xmax=sum(i,j)                                                     
                  nin=i                                                             
               end if                                                              
            end do
C     label nin with the current vt
            nvt(nin)=j-1           
C     remove nin from the sum, since it has already been labeled
            do ii=1,9                                                          
               sum(nin,ii)=0.d0                                                    
            end do
         end do
      end do

c     rotational labeling. the same scheme as isabelle's in resum2.               

c     loop through all of the eigenvectors
      do jj=1,9*(2*n+1)                                                       
         xmax=0.d0                                                               
         kkk(jj)=-n                       
c     loop through the vt block containing jj
         do ii=(2*n+1)*nvt(jj)+1,(2*n+1)*nvt(jj)+(2*n+1)                       
C     find the vector with the largest K component
            if (egvc(ii,jj)*egvc(ii,jj).gt.xmax) then                             
               xmax=egvc(ii,jj)*egvc(ii,jj)                                        
               kkk(jj)=ii-(n+1)-nvt(jj)*(2*n+1)                                    
            end if                                                                
         end do
      end do
      
      do i=1,9*(2*n+1)                                                        
         if (nvt(i).eq.2) then                                                   
            if (n.ge.19.and.abs(kkk(i)).eq.9) then                                
               kkk(i)=-kkk(i)                                                      
            end if                                                                
         end if                                                                  
         if (kk.eq.kkk(i).and.iv.eq.nvt(i)) then                                 
            ipose=i                                                               
            
         end if                                                                  
      end do
      return                                                                    
      end    
C     *************************************************************************


