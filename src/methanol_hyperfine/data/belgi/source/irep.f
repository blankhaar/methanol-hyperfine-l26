C***********************************************************                    
      FUNCTION IREP(SYMBOL)                                                     
      implicit none
      save
C     THIS FUNCTION RECOGNIZE THE SYMBOLIC CONSTANT ENCOUNTERED                 
C     ON TAPE 5 AND RETURNS ITS POSITION IN THE ABC STRING                      


      include 'common_symb.f'
      include 'common_cstdim.f'

      integer irep, k
      CHARACTER SYMBOL*6
      DO K=1,cstm                                                    
         IF(SYMBOL.EQ.ABC(6*(K-1)+1:6*K))  THEN                                 
            IREP=K                                                              
            RETURN                                                              
         END IF                                                                 
      end do
      WRITE(6,100) SYMBOL                                                       
 100  FORMAT(//,' NO MATCHING CONSTANT SYMBOL FOR ',A6,/)                       
      STOP                                                                      
      END                                                                       
