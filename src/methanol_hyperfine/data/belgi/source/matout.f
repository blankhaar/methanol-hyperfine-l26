C     *************************************************************************
      subroutine matou2(x,z,m,n,mm,nn)                                          
C     *************************************************************************
C     *************************************************************************
      implicit none
      save

      double precision x, z
      integer m, n, mm, nn
      integer iflg, ilower, iupper, i, j

      dimension x(m,n),z(m)                                                     
 1004 format(/)                                                                 
 1003 format(4x,9f14.7)                                                         
 1000 format(9(11x,i3))                                                         
 1002 format(i4,9f14.6)                                                         
c                                                                               
      iflg=1                                                                    
      ilower=1                                                                  
      do while (iflg .ne. 0)
         iupper=ilower+8                                                           
         if( (iupper-nn) .ge. 0) then                                                        
            iupper=nn                                                                 
            iflg=0                      
         endif
         write(6,1000)(j,j=ilower,iupper)                                          
         write(6,1004)                                                             
         write(6,1003) (z(j),j=ilower,iupper)                                      
         write(6,1004)                                                             
         do i=1,mm                                                               
            write(6,1002)i,(x(i,j),j=ilower,iupper)                                
         end do
         write(6,1004)                                                             
         ilower=ilower+9                                                           
      end do
      return                                                                    
      end                                                                      
C     *************************************************************************


C     ************************************************************************* 
C*************************************************************                  
      subroutine matout(x,m,n,mm,nn)                                           
C     *************************************************************************
C     ************************************************************************* 
      implicit none
      save

      double precision x
      integer m, n, mm, nn
      integer iflg, ilower, iupper, i, j

      dimension x(m,n)                                                          
 1004 format(/)                                                                 
 1000 format(9(11x,i3))                                                         
 1002 format(i4,9f14.6)                                                         
                                                                               
      iflg=1                                                                    
      ilower=1                                                                  
      iupper=ilower+8                                                           
      do while (iflg .ne. 0)
         if ( (iupper-nn) .ge. 0 ) then                                                        
            iupper=nn                                                                 
            iflg=0                        
         end if
         write(6,1000)(j,j=ilower,iupper)                                          
         write(6,1004)                                                             
         do i=1,mm                                                               
            write(6,1002)i,(x(i,j),j=ilower,iupper)                                
         end do                                                                  
         write(6,1004)                                                             
        ilower=ilower+9                                                           
      end do                                                                   
      return                                                                    
      end                                                                       
C     *************************************************************************

