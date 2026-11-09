C*********************************************************
      SUBROUTINE LINECA (CST,NPP)                                               
      implicit none
      save
C     THIS SUBROUTINE CALCULATES AND WRITES                                     
C     1) THE OBS-CALC FOR ALL THE EXPERIMENTAL DATA                             
C     2) THE ENERGY OF THE LINES GIVEN                                          
C     
C*********************************************************

      include 'common_dim.f'
      include 'common_eigen.f'
      include 'common_eigen1.f'
      include 'common_tor.f'
      include 'common_dat.f'
      include 'common_rotor.f'
      include 'common_big2.f'
      include 'common_expdat.f'
      include 'common_symb.f'
      include 'common_quant.f'
      include 'common_chara.f'
      include 'common_cst.f'
      include 'common_int.f'

      character*6 muabc

      integer ni, i

      integer iblsto, kc, lo, ipr, ivv, ink, kt1k, kt2k, npp, k
      double precision ee, smw, smw0, smw1, smw2, smw3, sir, sir10
      double precision sir21, sig004, sig080, sig010, sig100, sig020
      double precision sig200, sig045, sig1m0, wmw, wmw0, wmw1, wmw2
      double precision wmw3, wir, wir10, wir21, wig004, wig080, wig010
      double precision wig100, wig020, wig200, wig045, wig1m0, wig070
      double precision sgst0, ecalc, eobs, tme2, wgst0, tme3, ccc, s
      character ast

      double precision sft00, sft01, sft02, sft11, sft12, sft22
      double precision wft00, wft01, wft02, wft11, wft12, wft22
      double precision wft00a, wft01a, wft02a, wft11a, wft12a, wft22a
      double precision wft00e, wft01e, wft02e, wft11e, wft12e, wft22e
      double precision sft, wft
      double precision sft1, sft2, sft3, sft4
      double precision sft1a, sft2a, sft3a, sft4a
      double precision sft1e, sft2e, sft3e, sft4e

      double precision smw00, smw01, smw02, smw11, smw12, smw22
      double precision wmw00, wmw01, wmw02, wmw11, wmw12, wmw22
      double precision wmw00a, wmw01a, wmw02a, wmw11a, wmw12a, wmw22a
      double precision wmw00e, wmw01e, wmw02e, wmw11e, wmw12e, wmw22e
      double precision smwj, wmwj

      double precision sig004a, sig080a, sig010a, sig100a, sig020a
      double precision sig200a, sig045a, sig1m0a, sw070a

      double precision sig004e, sig080e, sig010e, sig100e, sig020e
      double precision sig200e, sig045e, sig1m0e, sw070e

      integer nft, nmw, inkf

      double precision G, BETA, BLNU, RINT, cf
      integer ictt(6)

      DIMENSION IBLSTO(2),EE(2),KC(2,NDTA), ivv(2), ink(2), inkf(2)

      S=0.d0                                                                    
      SMW=0.d0                                                                  
      SMW0=0.d0                                                                 
      SMW1=0.d0                                                                 
      SMW2=0.d0                                                                 
      SMW3=0.d0                                                                 
      SIR=0.d0                                                                  
      SIR10=0.d0                                                                
      SIR21=0.d0                                                                
      SIG004=0.d0                                                               
      SIG080=0.d0                                                              
      SIG010=0.d0                                                              
      SIG100=0.d0                                                              
      SIG020=0.d0                                                              
      SIG200=0.d0                                                              
      SIG045=0.d0                                                              
      SIG1M0=0.d0                                                              
      SW070=0.d0                                                               
      WMW=0.d0                                                                  
      WMW0=0.d0                                                                 
      WMW1=0.d0                                                                 
      WMW2=0.d0                                                                 
      WMW3=0.d0                                                                 
      WIR=0.d0                                                                  
      WIR10=0.d0                                                                
      WIR21=0.d0                                                                
      WIG004=0.d0                                                               
      WIG080=0.d0                                                              
      WIG010=0.d0                                                              
      WIG100=0.d0                                                               
      WIG020=0.d0                                                               
      WIG200=0.d0                                                               
      WIG045=0.d0                                                               
      WIG1M0=0.d0                                                               
      WIG070=0.d0                                                               
      SGST0=0.d0                                                                


      SFT00 = 0.d0
      SFT01 = 0.d0
      SFT02 = 0.d0
      SFT11 = 0.d0
      SFT12 = 0.d0
      SFT22 = 0.d0

      WFT00 = 0.d0
      WFT01 = 0.d0
      WFT02 = 0.d0
      WFT11 = 0.d0
      WFT12 = 0.d0
      WFT22 = 0.d0

      WFT00A = 0.d0
      WFT01A = 0.d0
      WFT02A = 0.d0
      WFT11A = 0.d0
      WFT12A = 0.d0
      WFT22A = 0.d0

      WFT00E = 0.d0
      WFT01E = 0.d0
      WFT02E = 0.d0
      WFT11E = 0.d0
      WFT12E = 0.d0
      WFT22E = 0.d0

      SFT = 0.d0
      WFT = 0.d0
      
      SFT1 = 0.d0
      SFT2 = 0.d0
      SFT3 = 0.d0
      SFT4 = 0.d0

      SFT1A = 0.d0
      SFT2A = 0.d0
      SFT3A = 0.d0
      SFT4A = 0.d0

      SFT1E = 0.d0
      SFT2E = 0.d0
      SFT3E = 0.d0
      SFT4E = 0.d0

      SMW00 = 0.d0
      SMW01 = 0.d0
      SMW02 = 0.d0
      SMW11 = 0.d0
      SMW12 = 0.d0
      SMW22 = 0.d0

      WMW00 = 0.d0
      WMW01 = 0.d0
      WMW02 = 0.d0
      WMW11 = 0.d0
      WMW12 = 0.d0
      WMW22 = 0.d0

      WMW00A = 0.d0
      WMW01A = 0.d0
      WMW02A = 0.d0
      WMW11A = 0.d0
      WMW12A = 0.d0
      WMW22A = 0.d0

      WMW00E = 0.d0
      WMW01E = 0.d0
      WMW02E = 0.d0
      WMW11E = 0.d0
      WMW12E = 0.d0
      WMW22E = 0.d0

      SMWJ = 0.d0
      WMWJ = 0.d0

      SIG004A=0.d0                                                               
      SIG080A=0.d0                                                              
      SIG010A=0.d0                                                              
      SIG100A=0.d0                                                              
      SIG020A=0.d0                                                              
      SIG200A=0.d0                                                              
      SIG045A=0.d0                                                              
      SIG1M0A=0.d0                                                              

      SIG004E=0.d0                                                               
      SIG080E=0.d0                                                              
      SIG010E=0.d0                                                              
      SIG100E=0.d0                                                              
      SIG020E=0.d0                                                              
      SIG200E=0.d0                                                              
      SIG045E=0.d0                                                              
      SIG1M0E=0.d0                                                              

      IBLSTO(1)=0                                                               
      IBLSTO(2)=0                                                               
      WRITE(6,602)                                                              

      DO 22 K=1,NDATA                                                           
         IF(NOBS(1,K).GT.JMAX.OR.NOBS(2,K).GT.JMAX) GO TO 22                       
         if(ivobs(1,k).gt.vtmax.or.ivobs(2,k).gt.vtmax) go to 22
C     IF(K.GT.NADAT) THEN                                                       
C     ISIG=1                                                                    
C     ELSE                                                                      
C     ISIG=0                                                                    
C     ENDIF                                                                     
C     DO 24 II=1,2                                                              
C     IF(II.EQ.1) THEN                                                          
C     NI=1                                                                      
C     M=1                                                                       
C     ELSE                                                                      
C     NI=NDMX+1                                                                 
C     M=KDMX+1                                                                  
C     ENDIF                                                                     
C     N=NOBS(II,K)                                                              
C     IV=IVOBS(II,K)                                                            
C     KK=KOBS(II,K)                                                             
C     PR=PROBS(II,K)                                                            
C     LOO=II                                                                    
C     III=II                                                                    
C     IF(IRMW.EQ.1.AND.II.EQ.2) III=1                                           
C     IF(IBLK(LOO,K).NE.IBLSTO(LOO)) THEN                                       
C     IBLSTO(LOO)=IBLK(LOO,K)                                                   
C     CALL ROTTOR(CST,III,EGVL(NI),EGVC(1,NI),ETOR(M,1)                         
C     1        ,A(1,1,M))                                                       
C     ENDIF                                                                     
C     NROTOR=(2*N+1)*9                                                          
C     CALL ENCAL(E,LOO,EGVL(NI),EGVC(1,NI))                                     
C     EE(LOO)=E                                                                 
C     24   CONTINUE                                                                  
         DO 31 LO=1,2                                                           
            N=NOBS(LO,K)                                                        
            KK=KOBS(LO,K)                                                       
            IV=IVOBS(LO,K)                                                      
            PR=PROBS(LO,K)                                                      
            IF(PR.EQ.'+') IPR=1                                                 
            IF(PR.EQ.'-') IPR=0                                                 
            IF(K.GT.NADAT) THEN                                                 
               ISIG=1                                                           

               IVV(lo)=IV+6                                                         
               INK(lo)=N**2+N+KK+1                                                  
               inkf(1) = ink(1)
            ELSE                                                                
               ISIG=0                                                           
               IVV(lo)=IV+1                                                         
               INK(lo)=N**2+N+KK+1-IPR*N                                            
            ENDIF                                                               
            IF(IRMW.EQ.2) THEN
               EE(LO)=ENER(IVV(lo),INK(lo),LO) - en0   
C               EE(LO)=ENER(IVV(lo),INK(lo),LO)                                             
            ELSEIF(IRMW.EQ.1) THEN
               EE(LO)=ENER(IVV(lo),INK(lo),1) - en0                                       
C               EE(LO)=ENER(IVV(lo),INK(lo),1)                                              
            ENDIF
 31      CONTINUE                                                               
         ECALC=TERM+EE(2)-EE(1)                                                 
         iisig = isig+1

C
C     Calculate the Boltzman and frequency factor of the
C     intensity
C     
         G=1.0
         IF (INT.EQ.1)THEN
            BETA = (CONV2/TPR) * (-EE(1) - EVIB(1)) 
            BLTZ = EXP(BETA) 
            BLNU = 1 - EXP(-ECALC*(CONV2/TPR))
            BLTZ2 = BLTZ*ECALC*BLNU*G
C     print*,'BLTZ=',BLTZ
C     PRINT*,'ECALC=',ECALC,'CONV2=',CONV2,'BLTZ2=',BLTZ2
         END IF
C     HERE, THE USER HAS THE POSSIBILITY OF HAVING A PRINT OUT                  
C     INCLUDING THE JKAKC NOTATION AND THE INTERNAL ROTOR NOTATION              
C     BUT, FOR THE E STATES,SINCE THERE IS NO SYMMETRY RULES ALLOW              
C     THE CONVERSION, IT IS BASED ON ENERGY LEVELS ORDER,M WHICH                
C     IS MOLECULE DEPENDANT:DHERE IT IS FOR CH3CHO ONLY!                        
         IF(IKAKC.NE.1) GO TO 60                                                
         IF(ISIG.EQ.0) THEN                                                     
            IF(MOD(KOBS(1,K),2).EQ.0) THEN                                      
               IF(PROBS(1,K).EQ.'+') THEN                                       
                  KC(1,K)=NOBS(1,K)-KOBS(1,K)                                   
               ELSE                                                             
                  KC(1,K)=NOBS(1,K)-KOBS(1,K)+1                                 
               ENDIF                                                            
            ELSEIF(MOD(KOBS(1,K),2).NE.0) THEN                                  
               IF(PROBS(1,K).EQ.'+') THEN                                       
                  KC(1,K)=NOBS(1,K)-KOBS(1,K)+1                                 
               ELSE                                                             
                  KC(1,K)=NOBS(1,K)-KOBS(1,K)                                   
               ENDIF                                                            
            ENDIF                                                               
            IF(MOD(KOBS(2,K),2).EQ.0) THEN                                      
               IF(PROBS(2,K).EQ.'+') THEN                                       
                  KC(2,K)=NOBS(2,K)-KOBS(2,K)                                   
               ELSE                                                             
                  KC(2,K)=NOBS(2,K)-KOBS(2,K)+1                                 
               ENDIF                                                            
            ELSEIF(MOD(KOBS(2,K),2).NE.0) THEN                                  
               IF(PROBS(2,K).EQ.'+') THEN                                       
                  KC(2,K)=NOBS(2,K)-KOBS(2,K)+1                                 
               ELSE                                                             
                  KC(2,K)=NOBS(2,K)-KOBS(2,K)                                   
               ENDIF                                                            
            ENDIF                                                               
         ELSE                                                                   
            KT1K=KOBS(1,K)*(-1)**IVOBS(1,K)                                     
            IF(KT1K.EQ.0.OR.KT1K.EQ.-1) THEN                                    
               KC(1,K)=NOBS(1,K)                                                
            ELSEIF(KT1K.EQ.1.OR.KT1K.EQ.-2) THEN                                
               KC(1,K)=NOBS(1,K)-1                                              
            ELSEIF(KT1K.EQ.2.OR.KT1K.EQ.-3) THEN                                
               KC(1,K)=NOBS(1,K)-2                                              
            ELSEIF(KT1K.EQ.3.OR.KT1K.EQ.-4) THEN                                
               KC(1,K)=NOBS(1,K)-3                                              
            ELSEIF(KT1K.EQ.4.OR.KT1K.EQ.5) THEN                                 
               KC(1,K)=NOBS(1,K)-4                                              
            ELSEIF(KT1K.EQ.-5.OR.KT1K.EQ.6) THEN                                
               KC(1,K)=NOBS(1,K)-5                                              
            ELSEIF(KT1K.EQ.-6.OR.KT1K.EQ.7) THEN                                
               KC(1,K)=NOBS(1,K)-6                                              
            ELSEIF(KT1K.EQ.-7.OR.KT1K.EQ.8) THEN                                
               KC(1,K)=NOBS(1,K)-7                                              
            ELSEIF(KT1K.EQ.-8.OR.KT1K.EQ.9) THEN                                
               KC(1,K)=NOBS(1,K)-8                                              
            ELSEIF(KT1K.EQ.-9.OR.KT1K.EQ.10) THEN                               
               KC(1,K)=NOBS(1,K)-9                                              
            ELSEIF(KT1K.EQ.-10.OR.KT1K.EQ.11) THEN                              
               KC(1,K)=NOBS(1,K)-10                                             
            ELSEIF(KT1K.EQ.-11.OR.KT1K.EQ.12) THEN                              
               KC(1,K)=NOBS(1,K)-11                                             
            ELSEIF(KT1K.EQ.-12.OR.KT1K.EQ.13) THEN                              
               KC(1,K)=NOBS(1,K)-12                                             
            ELSEIF(KT1K.EQ.-13.OR.KT1K.EQ.14) THEN                              
               KC(1,K)=NOBS(1,K)-13                                             
            ELSEIF(KT1K.EQ.-14.OR.KT1K.EQ.15) THEN                              
               KC(1,K)=NOBS(1,K)-14                                             
            ELSEIF(KT1K.EQ.-15.OR.KT1K.EQ.16) THEN                              
               KC(1,K)=NOBS(1,K)-15                                             
            ENDIF                                                               
            KT2K=KOBS(2,K)*(-1)**IVOBS(2,K)                                     
            IF(KT2K.EQ.0.OR.KT2K.EQ.-1) THEN                                    
               KC(2,K)=NOBS(2,K)                                                
            ELSEIF(KT2K.EQ.1.OR.KT2K.EQ.-2) THEN                                
               KC(2,K)=NOBS(2,K)-1                                              
            ELSEIF(KT2K.EQ.2.OR.KT2K.EQ.-3) THEN                                
               KC(2,K)=NOBS(2,K)-2                                              
            ELSEIF(KT2K.EQ.3.OR.KT2K.EQ.-4) THEN                                
               KC(2,K)=NOBS(2,K)-3                                              
            ELSEIF(KT2K.EQ.4.OR.KT2K.EQ.5) THEN                                 
               KC(2,K)=NOBS(2,K)-4                                              
            ELSEIF(KT2K.EQ.-5.OR.KT2K.EQ.6) THEN                                
               KC(2,K)=NOBS(2,K)-5                                              
            ELSEIF(KT2K.EQ.-6.OR.KT2K.EQ.7) THEN                                
               KC(2,K)=NOBS(2,K)-6                                              
            ELSEIF(KT2K.EQ.-7.OR.KT2K.EQ.8) THEN                                
               KC(2,K)=NOBS(2,K)-7                                              
            ELSEIF(KT2K.EQ.-8.OR.KT2K.EQ.9) THEN                                
               KC(2,K)=NOBS(2,K)-8                                              
            ELSEIF(KT2K.EQ.-9.OR.KT2K.EQ.10) THEN                               
               KC(2,K)=NOBS(2,K)-9                                              
            ELSEIF(KT2K.EQ.-10.OR.KT2K.EQ.11) THEN                              
               KC(2,K)=NOBS(2,K)-10                                             
            ELSEIF(KT2K.EQ.-11.OR.KT2K.EQ.12) THEN                              
               KC(2,K)=NOBS(2,K)-11                                             
            ELSEIF(KT2K.EQ.-12.OR.KT2K.EQ.13) THEN                              
               KC(2,K)=NOBS(2,K)-12                                             
            ELSEIF(KT2K.EQ.-13.OR.KT2K.EQ.14) THEN                              
               KC(2,K)=NOBS(2,K)-13                                             
            ELSEIF(KT2K.EQ.-14.OR.KT2K.EQ.15) THEN                              
               KC(2,K)=NOBS(2,K)-14                                             
            ELSEIF(KT2K.EQ.-15.OR.KT2K.EQ.16) THEN                              
               KC(2,K)=NOBS(2,K)-15                                             
            ENDIF                                                               
         ENDIF                                                                  
 60      IF(IAST(K).EQ.1) THEN                                                  
            EOBS=ETRANS(K)                                                      
            TME2=EOBS-ECALC                                                     
            S=S+TME2**2*W(K)                                                    

            if ( iw(k) .lt. 1000 ) then

               SFT = SFT + TME2**2
               WFT = WFT + TME2**2 * W(K)

               if ( iw(k) .eq. 1 ) SFT1 = SFT1 + TME2**2
               if ( iw(k) .eq. 2 ) SFT2 = SFT2 + TME2**2
               if ( iw(k) .eq. 3 ) SFT3 = SFT3 + TME2**2
               if ( iw(k) .eq. 4 ) SFT4 = SFT4 + TME2**2

               if ( probs(1, k) .eq. ' ' ) then
                  if ( iw(k) .eq. 1 ) SFT1E = SFT1E + TME2**2
                  if ( iw(k) .eq. 2 ) SFT2E = SFT2E + TME2**2
                  if ( iw(k) .eq. 3 ) SFT3E = SFT3E + TME2**2
                  if ( iw(k) .eq. 4 ) SFT4E = SFT4E + TME2**2
               else
                  if ( iw(k) .eq. 1 ) SFT1A = SFT1A + TME2**2
                  if ( iw(k) .eq. 2 ) SFT2A = SFT2A + TME2**2
                  if ( iw(k) .eq. 3 ) SFT3A = SFT3A + TME2**2
                  if ( iw(k) .eq. 4 ) SFT4A = SFT4A + TME2**2
               endif

               if ( (ivobs(1, k).eq.0) .and. (ivobs(2,k).eq.0)) then
                  WFT00=WFT00+TME2**2*W(K)                                             
                  SFT00=SFT00+TME2**2                                                  
               else if ( ( (ivobs(1, k).eq.0) 
     &                 .and. (ivobs(2,k).eq.1) )
     &                 .or. ( (ivobs(1,k).eq.1) 
     &                 .and. (ivobs(2,k).eq.0)) ) then
                  WFT01=WFT01+TME2**2*W(K)                                             
                  SFT01=SFT01+TME2**2                                                  
               else if ( ( (ivobs(1, k).eq.0) 
     &                 .and. (ivobs(2,k).eq.2) )
     &                 .or. ( (ivobs(1,k).eq.2) 
     &                 .and. (ivobs(2,k).eq.0)) ) then
                  WFT02=WFT02+TME2**2*W(K)                                             
                  SFT02=SFT02+TME2**2                                                  
               else if ((ivobs(1, k).eq.1).and.(ivobs(2,k).eq.1)) then
                  WFT11=WFT11+TME2**2*W(K)                                             
                  SFT11=SFT11+TME2**2                                                  
               else if ( ( (ivobs(1, k).eq.1) 
     &                 .and. (ivobs(2,k).eq.2) )
     &                 .or. ( (ivobs(1,k).eq.2) 
     &                 .and. (ivobs(2,k).eq.1)) ) then
                  WFT12=WFT12+TME2**2*W(K)                                             
                  SFT12=SFT12+TME2**2                                                  
               else if ((ivobs(1, k).eq.2).and.(ivobs(2,k).eq.2)) then
                  WFT22=WFT22+TME2**2*W(K)                                             
                  SFT22=SFT22+TME2**2                                                  
               endif
            endif

            if ( probs(1, k) .eq. ' ' ) then
               if ( iw(k) .lt. 1000 ) then
                  if ( (ivobs(1, k).eq.0) .and. (ivobs(2,k).eq.0)) then
                     WFT00E=WFT00E+TME2**2*W(K)                                             
                  else if ( ( (ivobs(1, k).eq.0) 
     &                    .and. (ivobs(2,k).eq.1) )
     &                    .or. ( (ivobs(1,k).eq.1) 
     &                    .and. (ivobs(2,k).eq.0)) ) then
                     WFT01E=WFT01E+TME2**2*W(K)                                             
                  else if ( ( (ivobs(1, k).eq.0) 
     &                    .and. (ivobs(2,k).eq.2) )
     &                    .or. ( (ivobs(1,k).eq.2) 
     &                    .and. (ivobs(2,k).eq.0)) ) then
                     WFT02E=WFT02E+TME2**2*W(K)                                             
                  else if ((ivobs(1, k).eq.1).and.(ivobs(2,k).eq.1))
     &                    then
                     WFT11E=WFT11E+TME2**2*W(K)                                             
                  else if ( ( (ivobs(1, k).eq.1) 
     &                    .and. (ivobs(2,k).eq.2) )
     &                    .or. ( (ivobs(1,k).eq.2) 
     &                    .and. (ivobs(2,k).eq.1)) ) then
                     WFT12E=WFT12E+TME2**2*W(K)                                             
                  else if ((ivobs(1, k).eq.2).and.(ivobs(2,k).eq.2)) 
     &                    then
                     WFT22E=WFT22E+TME2**2*W(K)                                             
                  endif
               endif
            else
               if ( iw(k) .lt. 1000 ) then
                  if ( (ivobs(1, k).eq.0) .and. (ivobs(2,k).eq.0)) then
                     WFT00A=WFT00A+TME2**2*W(K)                                             
                  else if ( ( (ivobs(1, k).eq.0) 
     &                    .and. (ivobs(2,k).eq.1) )
     &                    .or. ( (ivobs(1,k).eq.1) 
     &                    .and. (ivobs(2,k).eq.0)) ) then
                     WFT01A=WFT01A+TME2**2*W(K)                                             
                  else if ( ( (ivobs(1, k).eq.0) 
     &                    .and. (ivobs(2,k).eq.2) )
     &                    .or. ( (ivobs(1,k).eq.2) 
     &                    .and. (ivobs(2,k).eq.0)) ) then
                     WFT02A=WFT02A+TME2**2*W(K)                                             
                  else if ((ivobs(1, k).eq.1).and.(ivobs(2,k).eq.1))
     &                    then
                     WFT11A=WFT11A+TME2**2*W(K)                                             
                  else if ( ( (ivobs(1, k).eq.1) 
     &                    .and. (ivobs(2,k).eq.2) )
     &                    .or. ( (ivobs(1,k).eq.2) 
     &                    .and. (ivobs(2,k).eq.1)) ) then
                     WFT12A=WFT12A+TME2**2*W(K)                                             
                  else if ((ivobs(1, k).eq.2).and.(ivobs(2,k).eq.2))
     &                    then
                     WFT22A=WFT22A+TME2**2*W(K)                                             
                  endif
               endif
            endif


            if ( iw(k) .ge. 1000 ) then

               SMWJ = SMWJ + TME2**2
               WMWJ = WMWJ + TME2**2 * W(K)

               if ( (ivobs(1, k).eq.0) .and. (ivobs(2,k).eq.0)) then
                  WMW00=WMW00+TME2**2*W(K)                                             
                  SMW00=SMW00+TME2**2                                                  
               else if ( ( (ivobs(1, k).eq.0) 
     &                 .and. (ivobs(2,k).eq.1) )
     &                 .or. ( (ivobs(1,k).eq.1) 
     &                 .and. (ivobs(2,k).eq.0)) ) then
                  WMW01=WMW01+TME2**2*W(K)                                             
                  SMW01=SMW01+TME2**2                                                  
               else if ( ( (ivobs(1, k).eq.0) 
     &                 .and. (ivobs(2,k).eq.2) )
     &                 .or. ( (ivobs(1,k).eq.2) 
     &                 .and. (ivobs(2,k).eq.0)) ) then
                  WMW02=WMW02+TME2**2*W(K)                                             
                  SMW02=SMW02+TME2**2                                                  
               else if ((ivobs(1, k).eq.1).and.(ivobs(2,k).eq.1)) then
                  WMW11=WMW11+TME2**2*W(K)                                             
                  SMW11=SMW11+TME2**2                                                  
               else if ( ( (ivobs(1, k).eq.1) 
     &                 .and. (ivobs(2,k).eq.2) )
     &                 .or. ( (ivobs(1,k).eq.2) 
     &                 .and. (ivobs(2,k).eq.1)) ) then
                  WMW12=WMW12+TME2**2*W(K)                                             
                  SMW12=SMW12+TME2**2                                                  
               else if ((ivobs(1, k).eq.2).and.(ivobs(2,k).eq.2)) then
                  WMW22=WMW22+TME2**2*W(K)                                             
                  SMW22=SMW22+TME2**2                                                  
               endif
            endif

            if ( probs(1, k) .eq. ' ' ) then
               if ( iw(k) .ge. 1000 ) then
                  if ( (ivobs(1, k).eq.0) .and. (ivobs(2,k).eq.0)) then
                     WMW00E=WMW00E+TME2**2*W(K)                                             
                  else if ( ( (ivobs(1, k).eq.0) 
     &                    .and. (ivobs(2,k).eq.1) )
     &                    .or. ( (ivobs(1,k).eq.1) 
     &                    .and. (ivobs(2,k).eq.0)) ) then
                     WMW01E=WMW01E+TME2**2*W(K)                                             
                  else if ( ( (ivobs(1, k).eq.0) 
     &                    .and. (ivobs(2,k).eq.2) )
     &                    .or. ( (ivobs(1,k).eq.2) 
     &                    .and. (ivobs(2,k).eq.0)) ) then
                     WMW02E=WMW02E+TME2**2*W(K)                                             
                  else if ((ivobs(1, k).eq.1).and.(ivobs(2,k).eq.1))
     &                    then
                     WMW11E=WMW11E+TME2**2*W(K)                                             
                  else if ( ( (ivobs(1, k).eq.1) 
     &                    .and. (ivobs(2,k).eq.2) )
     &                    .or. ( (ivobs(1,k).eq.2) 
     &                    .and. (ivobs(2,k).eq.1)) ) then
                     WMW12E=WMW12E+TME2**2*W(K)                                             
                  else if ((ivobs(1, k).eq.2).and.(ivobs(2,k).eq.2)) 
     &                    then
                     WMW22E=WMW22E+TME2**2*W(K)                                             
                  endif
               endif
            else
               if ( iw(k) .ge. 1000 ) then
                  if ( (ivobs(1, k).eq.0) .and. (ivobs(2,k).eq.0)) then
                     WMW00A=WMW00A+TME2**2*W(K)                                             
                  else if ( ( (ivobs(1, k).eq.0) 
     &                    .and. (ivobs(2,k).eq.1) )
     &                    .or. ( (ivobs(1,k).eq.1) 
     &                    .and. (ivobs(2,k).eq.0)) ) then
                     WMW01A=WMW01A+TME2**2*W(K)                                             
                  else if ( ( (ivobs(1, k).eq.0) 
     &                    .and. (ivobs(2,k).eq.2) )
     &                    .or. ( (ivobs(1,k).eq.2) 
     &                    .and. (ivobs(2,k).eq.0)) ) then
                     WMW02A=WMW02A+TME2**2*W(K)                                             
                  else if ((ivobs(1, k).eq.1).and.(ivobs(2,k).eq.1))
     &                    then
                     WMW11A=WMW11A+TME2**2*W(K)                                             
                  else if ( ( (ivobs(1, k).eq.1) 
     &                    .and. (ivobs(2,k).eq.2) )
     &                    .or. ( (ivobs(1,k).eq.2) 
     &                    .and. (ivobs(2,k).eq.1)) ) then
                     WMW12A=WMW12A+TME2**2*W(K)                                             
                  else if ((ivobs(1, k).eq.2).and.(ivobs(2,k).eq.2))
     &                    then
                     WMW22A=WMW22A+TME2**2*W(K)                                             
                  endif
               endif
            endif

            IF(IW(K).GE.1000) THEN                                              
               WMW=WMW+TME2**2*W(K)                                             
               SMW=SMW+TME2**2                                                  

               if ( probs(1, k) .eq. ' ' ) then
                  IF(IW(K).EQ.1004.OR.IW(K).EQ.1005) THEN                          
                     SIG004E=SIG004E+TME2**2                                         
                  ELSEIF(IW(K).EQ.1080) THEN                                       
                     SIG080E=SIG080E+TME2**2                                         
                  ELSEIF(IW(K).EQ.1070) THEN                                       
                     SW070E=SW070E+TME2**2                                           
                  ELSEIF(IW(K).EQ.1010.OR.IW(K).EQ.1008) THEN                      
                     SIG010E=SIG010E+TME2**2                                         
                  ELSEIF(IW(K).EQ.1100) THEN                                       
                     SIG100E=SIG100E+TME2**2                                         
                  ELSEIF(IW(K).EQ.1020) THEN                                       
                     SIG020E=SIG020E+TME2**2                                         
                  ELSEIF(IW(K).LE.1300.AND.IW(K).GE.1130) THEN                 
                     SIG200E=SIG200E+TME2**2                                         
                  ELSEIF(IW(K).EQ.1030.OR.IW(K).EQ.1040.OR.
     &                    IW(K).EQ.1050) THEN                                                     
                     SIG045E=SIG045E+TME2**2                                         
                  ELSEIF(IW(K).EQ.1990) THEN                                       
                     SIG1M0E=SIG1M0E+TME2**2                                         
                  ENDIF                           
               else
                  IF(IW(K).EQ.1004.OR.IW(K).EQ.1005) THEN                          
                     SIG004A=SIG004A+TME2**2                                         
                  ELSEIF(IW(K).EQ.1080) THEN                                       
                     SIG080A=SIG080A+TME2**2                                         
                  ELSEIF(IW(K).EQ.1070) THEN                                       
                     SW070A=SW070A+TME2**2                                           
                  ELSEIF(IW(K).EQ.1010.OR.IW(K).EQ.1008) THEN                      
                     SIG010A=SIG010A+TME2**2                                         
                  ELSEIF(IW(K).EQ.1100) THEN                                       
                     SIG100A=SIG100A+TME2**2                                         
                  ELSEIF(IW(K).EQ.1020) THEN                                       
                     SIG020A=SIG020A+TME2**2                                         
                  ELSEIF(IW(K).LE.1300.AND.IW(K).GE.1130) THEN                 
                     SIG200A=SIG200A+TME2**2                                         
                  ELSEIF(IW(K).EQ.1030.OR.IW(K).EQ.1040.OR.
     &                    IW(K).EQ.1050) THEN                                                     
                     SIG045A=SIG045A+TME2**2                                         
                  ELSEIF(IW(K).EQ.1990) THEN                                       
                     SIG1M0A=SIG1M0A+TME2**2                                         
                  ENDIF                    
               endif                                 
               
               IF(IVOBS(1,K).EQ.0.AND.IW(K).GE.1000) THEN                       
                  WMW0=WMW0+TME2**2*W(K)                                        
                  SMW0=SMW0+TME2**2                                             
               ELSEIF(IVOBS(1,K).EQ.1.AND.IW(K).GE.1000) THEN                   
                  WMW1=WMW1+TME2**2*W(K)                                        
                  SMW1=SMW1+TME2**2                                             
               ELSEIF(IVOBS(1,K).EQ.2.AND.IW(K).GE.1000) THEN                   
                  WMW2=WMW2+TME2**2*W(K)                                        
                  SMW2=SMW2+TME2**2                                             
               ELSEIF(IVOBS(1,K).EQ.3.AND.IW(K).GE.1000) THEN                   
                  WMW3=WMW3+TME2**2*W(K)                                        
                  SMW3=SMW3+TME2**2                                             
C     ELSEIF(IW(K).EQ.2000) THEN                                                
C     SWGOOD=SWGOOD+TME2**2*W(K)                                                
               ENDIF                                                            
               IF(IW(K).EQ.1004.OR.IW(K).EQ.1005) THEN                          
                  WIG004=WIG004+TME2**2*W(K)                                    
                  SIG004=SIG004+TME2**2                                         
               ELSEIF(IW(K).EQ.1080) THEN                                       
                  WIG080=WIG080+TME2**2*W(K)                                    
                  SIG080=SIG080+TME2**2                                         
               ELSEIF(IW(K).EQ.1070) THEN                                       
                  WIG070=WIG070+TME2**2*W(K)                                    
                  SW070=SW070+TME2**2                                           
               ELSEIF(IW(K).EQ.1010.OR.IW(K).EQ.1008) THEN                      
                  WIG010=WIG010+TME2**2*W(K)                                    
                  SIG010=SIG010+TME2**2                                         
               ELSEIF(IW(K).EQ.1100) THEN                                       
                  WIG100=WIG100+TME2**2*W(K)                                    
                  SIG100=SIG100+TME2**2                                         
               ELSEIF(IW(K).EQ.1020) THEN                                       
                  WIG020=WIG020+TME2**2*W(K)                                    
                  SIG020=SIG020+TME2**2                                         
               ELSEIF(IW(K).LE.1300.AND.IW(K).GE.1130) THEN                 
                  WIG200=WIG200+TME2**2*W(K)                                    
                  SIG200=SIG200+TME2**2                                         
               ELSEIF(IW(K).EQ.1030.OR.IW(K).EQ.1040.OR.IW(K).EQ.1050)          
     &                 THEN                                                     
                  WIG045=WIG045+TME2**2*W(K)                                    
                  SIG045=SIG045+TME2**2                                         
               ELSEIF(IW(K).EQ.1990) THEN                                       
                  WIG1M0=WIG1M0+TME2**2*W(K)                                    
                  SIG1M0=SIG1M0+TME2**2                                         
               ENDIF                                                            
               EOBS=EOBS*29979.2458d0                                           
               ECALC=ECALC*29979.2458d0                                         
               TME2=TME2*29979.2458d0                                           
            ELSEIF(IW(K).EQ.1) THEN                                             
               WIR10=WIR10+TME2**2*W(K)                                         
               SIR10=SIR10+TME2**2                                              
               WIR=WIR+TME2**2*W(K)                                             
               SIR=SIR+TME2**2                                                  
            ELSEIF(IW(K).EQ.2) THEN                                             
               WIR21=WIR21+TME2**2*W(K)                                         
               SIR21=SIR21+TME2**2                                              
               WIR=WIR+TME2**2*W(K)                                             
               SIR=SIR+TME2**2                                                  
            ELSEIF(IW(K).EQ.10) THEN                                            
               WGST0=WGST0+TME2**2*W(K)                                         
               SGST0=SGST0+TME2**2                                              
            ENDIF                                                               


C
C            EE(1) = EE(1) - en0
C            EE(2) = EE(2) - en0
C

C     IF(IW(K).GE.1000) IW(K)=IW(K)-1000                                        
C     IF(IW(K).EQ.10) IW(K)=10000                                               
C     WRITE(6,604) IVOBS(2,K),NOBS(2,K),KOBS(2,K),KC(2,K)                       
C     1  ,PROBS(2,K),IVOBS(1,K),NOBS(1,K),KOBS(1,K),KC(1,K),                    
C     2  PROBS(1,K),EOBS,ECALC,TME2,IW(K),IAST(K),REF(K,1)                        

            IF(IW(K).GE.1000) THEN                                              
               IF(IW(K).EQ.2000) IW(K)=2                                        
               IF(IW(K).LT.2000) IW(K)=IW(K)-1000                               
C     
C     HERE THE USER HAS THE CHOICE BETWEEN A COMPACT FORMAT (THE ONE            
C     WHICH IS NOW PRECEDED BY "C" OR A LESS DENSE FORMAT                       
C     
               if ( int.ne.1 .and. unc .eq. 0) then
                  WRITE(6,704) IVOBS(2,K),NOBS(2,K),KOBS(2,K),KC(2,K)              
     1                 ,PROBS(2,K),IVOBS(1,K),NOBS(1,K),KOBS(1,K),
     2                 KC(1,K),PROBS(1,K),EOBS,IW(K),ECALC,TME2,
     3                 (REF(K, NI), NI = 1,reflen) 

               elseif ( int .eq. 1 ) then
                  CALL INTEN(RINT,CST,K,IVV,INK,MUABC,ICTT,CF)
                  
                  WRITE(6,834) IVOBS(2,K),NOBS(2,K),KOBS(2,K),KC(2,K)
     1                 ,PROBS(2,K),IVOBS(1,K),NOBS(1,K),KOBS(1,K),
     2                 KC(1,K), PROBS(1,K),EOBS,IW(K),ECALC,TME2,RINT,
     3                 MUABC,REF(K,1),EE(1),EE(2),BLTZ,SSQ,
     4                 (ICTT(I),I=1,6),CF
               elseif ( unc.eq.1) then
                  WRITE(6,504) IVOBS(2,K),NOBS(2,K),KOBS(2,K),KC(2,K)              
     1                 ,PROBS(2,K),IVOBS(1,K),NOBS(1,K),KOBS(1,K),
     2                 KC(1,K),PROBS(1,K),EOBS,IW(K),ECALC,TME2,
     3                 aincer(k), (REF(K, NI), NI = 1,2)                                     
               endif


C     WRITE(6,604) IVOBS(2,K),NOBS(2,K),KOBS(2,K),KC(2,K)                       
C     1  ,PROBS(2,K),IVOBS(1,K),NOBS(1,K),KOBS(1,K),KC(1,K),                    
C     2  PROBS(1,K),EOBS,IW(K),ECALC,TME2,IAST(K),REF(K,1)                        

            ELSEIF(IW(K).EQ.10) THEN                                            
               IW(K)=5                                                          
               if (int.ne.1 .and. unc .eq. 0) then
               WRITE(6,734) IVOBS(2,K),NOBS(2,K),KOBS(2,K),KC(2,K)              
     1              ,PROBS(2,K),IVOBS(1,K),NOBS(1,K),KOBS(1,K),KC(1,K),         
     2              PROBS(1,K),EOBS,IW(K),ECALC,TME2,IAST(K),
     3              (ref(K, ni), ni=1,reflen)

               elseif ( int .eq. 1 ) then
                  CALL INTEN(RINT,CST,K,IVV,INK,MUABC,ICTT,CF)
               WRITE(6,834) IVOBS(2,K),NOBS(2,K),KOBS(2,K),KC(2,K)
     1              ,PROBS(2,K),IVOBS(1,K),NOBS(1,K),KOBS(1,K),KC(1,K),
     2              PROBS(1,K),EOBS,IW(K),ECALC,TME2,RINT,MUABC,REF(K,1)
     3              ,EE(1),EE(2),BLTZ,SSQ,(ICTT(I),I=1,6),CF
               else if ( unc .eq. 1 ) then
                  WRITE(6,534) IVOBS(2,K),NOBS(2,K),KOBS(2,K),KC(2,K)              
     1              ,PROBS(2,K),IVOBS(1,K),NOBS(1,K),KOBS(1,K),KC(1,K),         
     2              PROBS(1,K),EOBS,IW(K),ECALC,TME2,IAST(K),
     3              aincer(k), (ref(K, ni), ni=1,2)
               endif
            ELSE                                              
C     IW(K)=4                                                          
C**** uncertainty changed to 2 for good methanol ir data!
C               IW(K)=2                                                
C     
C     HERE AGAIN THE USER CAN CHOICE A FORMAT                                   
C     
C     WRITE(6,635) IVOBS(2,K),NOBS(2,K),KOBS(2,K),KC(2,K)                       
C     1  ,PROBS(2,K),IVOBS(1,K),NOBS(1,K),KOBS(1,K),KC(1,K),                    
C     2  PROBS(1,K),EOBS,IW(K),ECALC,TME2,IAST(K),REF(K,1)                        

               if ( int.ne.1 .and. unc .eq. 0) then
               WRITE(6,735) IVOBS(2,K),NOBS(2,K),KOBS(2,K),KC(2,K)           
     1              ,PROBS(2,K),IVOBS(1,K),NOBS(1,K),KOBS(1,K),KC(1,K),         
     2              PROBS(1,K),EOBS,IW(K),ECALC,TME2,
     3              (ref(k,ni),ni=1,reflen)

               elseif ( int .eq. 1 ) then
                  CALL INTEN(RINT,CST,K,IVV,INK,MUABC,ICTT,CF)

                  WRITE(6,835) IVOBS(2,K),NOBS(2,K),KOBS(2,K),KC(2,K)
     1              ,PROBS(2,K),IVOBS(1,K),NOBS(1,K),KOBS(1,K),KC(1,K),
     2              PROBS(1,K),EOBS,IW(K),ECALC,TME2,RINT,MUABC,REF(K,1)
     3              ,EE(1),EE(2),BLTZ,SSQ,(ICTT(I),I=1,6),CF

               else if ( unc .eq. 1 ) then
                  WRITE(6,535) IVOBS(2,K),NOBS(2,K),KOBS(2,K),KC(2,K)           
     1               ,PROBS(2,K),IVOBS(1,K),NOBS(1,K),KOBS(1,K),KC(1,K),         
     2               PROBS(1,K),EOBS,IW(K),ECALC,TME2,
     3               aincer(k), (ref(k,ni),ni=1,2)
               endif
C           ELSEIF(IW(K).EQ.2)  THEN                                            
C               IW(K)=5                                                          
C     
C     HERE AGAIN THE USER CAN CHOICE A FORMAT                                   
C     
C     WRITE(6,635) IVOBS(2,K),NOBS(2,K),KOBS(2,K),KC(2,K)                       
C     1  ,PROBS(2,K),IVOBS(1,K),NOBS(1,K),KOBS(1,K),KC(1,K),                    
C     2  PROBS(1,K),EOBS,IW(K),ECALC,TME2,IAST(K),REF(K,1)                        
C               WRITE(6,735) IVOBS(2,K),NOBS(2,K),KOBS(2,K),KC(2,K)           
C     1              ,PROBS(2,K),IVOBS(1,K),NOBS(1,K),KOBS(1,K),KC(1,K),         
C     2              PROBS(1,K),EOBS,IW(K),ECALC,TME2,
C     3              (REF(K, ni), ni=1,reflen)
            ENDIF                                                            
         ELSEIF(IAST(K).NE.1) THEN                                           
            EOBS=ETRANS(K)                                                   
            TME3=EOBS-ECALC                                                  
            IF(IW(K).GE.1000) THEN                                           
               EOBS=EOBS*29979.2458d0
               ECALC=ECALC*29979.2458d0                                      
               TME3=TME3*29979.2458d0                                        
               AST='*'    

               if ( int.ne.1 .and. unc .eq. 0) then
                                                   
               WRITE(6,705) 
     1              IVOBS(2,K),NOBS(2,K),KOBS(2,K),KC(2,K),
     2              PROBS(2,K),
     3              IVOBS(1,K),NOBS(1,K),KOBS(1,K),KC(1,K),
     4              PROBS(1,K),
     5              EOBS,ECALC,TME3,AST,
     6              (REF(K,ni),ni=1,reflen)
               else if ( int .eq. 1 ) then
                  CALL INTEN(RINT,CST,K,IVV,INK,MUABC,ICTT,CF)
                   WRITE(6,805) IVOBS(2,K),NOBS(2,K),KOBS(2,K),KC(2,K),
     1              PROBS(2,K),IVOBS(1,K),NOBS(1,K),KOBS(1,K),KC(1,K),
     2              PROBS(1,K),EOBS,ECALC,TME3,RINT,MUABC,AST,REF(K,1)
     3              ,EE(1),EE(2),BLTZ,SSQ,(ICTT(I),I=1,6),CF
                else if ( unc .eq. 1 ) then
                   WRITE(6,505) 
     1              IVOBS(2,K),NOBS(2,K),KOBS(2,K),KC(2,K),
     2                  PROBS(2,K),
     3                  IVOBS(1,K),NOBS(1,K),KOBS(1,K),KC(1,K),
     4                  PROBS(1,K),
     5                  EOBS,ECALC,TME3,AST,
     6                  aincer(k), (REF(K,ni),ni=1,2)
               end if
            ELSE                                                             
               TME2=EOBS-ECALC                                               
               AST='*'    

               if ( int.ne.1 .and. unc .eq. 0) then
                                                   
               WRITE(6,708) IVOBS(2,K),NOBS(2,K),KOBS(2,K),KC(2,K),          
     1              PROBS(2,K),IVOBS(1,K),NOBS(1,K),KOBS(1,K),
     2              KC(1,K),PROBS(1,K),EOBS,ECALC,TME2,AST,
     3              (REF(K,ni), ni=1,reflen)
               else if ( int .eq. 1 ) then
                  CALL INTEN(RINT,CST,K,IVV,INK,MUABC,ICTT,CF)
                  WRITE(6,808) IVOBS(2,K),NOBS(2,K),KOBS(2,K),KC(2,K),
     1              PROBS(2,K),IVOBS(1,K),NOBS(1,K),KOBS(1,K),KC(1,K),
     2              PROBS(1,K),EOBS,ECALC,TME2,RINT,MUABC,AST,REF(K,1)
     3              ,EE(1),EE(2),BLTZ,SSQ,(ICTT(I),I=1,6),CF
               else if ( unc .eq. 1 ) then
               WRITE(6,508) IVOBS(2,K),NOBS(2,K),KOBS(2,K),KC(2,K),          
     1              PROBS(2,K),IVOBS(1,K),NOBS(1,K),KOBS(1,K),
     2              KC(1,K),PROBS(1,K),EOBS,ECALC,TME2,AST,
     3              aincer(k), (REF(K,ni), ni=1,2)
               endif
            ENDIF                                                            


         ENDIF                                                               
 22   CONTINUE                                                                  
      IF((NOFIT-NPP).LE.0) GO TO 23                                             
      S=SQRT(S/(NOFIT-NPP))                                                     
 23   IF(NFITMW.EQ.0) THEN                                                      
         WMW=0.d0                                                               
         SMW=0.d0                                                               
      ELSE                                                                      
         WMW=SQRT(WMW/NFITMW)                                                   
         SMW=SQRT(SMW/NFITMW)                                                   
      ENDIF                                                                     
      IF(NFITIR.EQ.0) THEN                                                      
         WIR=0.d0                                                               
         SIR=0.d0                                                               
      ELSE                                                                      
         WIR=SQRT(WIR/NFITIR)                                                   
         SIR=SQRT(SIR/NFITIR)                                                   
      ENDIF                                                                     
      IF(NFIR10.EQ.0) THEN                                                      
         WIR10=0.d0                                                             
         SIR10=0.d0                                                             
      ELSE                                                                      
         WIR10=SQRT(WIR10/NFIR10)                                               
         SIR10=SQRT(SIR10/NFIR10)                                               
      ENDIF                                                                     
      IF(NFIR21.EQ.0) THEN                                                      
         WIR21=0.d0                                                             
         SIR21=0.d0                                                             
      ELSE                                                                      
         WIR21=SQRT(WIR21/NFIR21)                                               
         SIR21=SQRT(SIR21/NFIR21)                                               
      ENDIF                                                                     
      IF(NFIMW0 .EQ.0) THEN                                                     
         WMW0=0.d0                                                              
         SMW0=0.d0                                                              
      ELSE                                                                      
         WMW0=SQRT(WMW0/NFIMW0 )                                                
         SMW0=SQRT(SMW0/NFIMW0 )                                                
      ENDIF                                                                     
      IF(NFI070 .EQ.0) THEN                                                     
         SW070=0.d0                                                             
         WIG070=0.d0                                                            
      ELSE                                                                      
         SW070=SQRT(SW070/NFI070 )                                              
         WIG070=SQRT(WIG070/NFI070 )                                            
      ENDIF                                                                     
      IF(NFIMW1 .EQ.0) THEN                                                     
         WMW1=0.d0                                                              
         SMW1=0.d0                                                              
      ELSE                                                                      
         WMW1=SQRT(WMW1/NFIMW1 )                                                
         SMW1=SQRT(SMW1/NFIMW1 )                                                
      ENDIF                                                                     
      IF(NFIMW2.EQ.0) THEN                                                      
         WMW2=0.d0                                                              
         SMW2=0.d0                                                              
      ELSE                                                                      
         WMW2=SQRT(WMW2/NFIMW2)                                                 
         SMW2=SQRT(SMW2/NFIMW2)                                                 
      ENDIF                                                                     
      IF(NFIMW3.EQ.0) THEN                                                      
         WMW3=0.d0                                                              
         SMW3=0.d0                                                              
      ELSE
         WMW3=SQRT(WMW3/NFIMW3)                                                 
         SMW3=SQRT(SMW3/NFIMW3)                                                 
      ENDIF
      IF(NFI004 .EQ.0) THEN                                                     
         WIG004=0.d0                                                            
         SIG004=0.d0                                                            
      ELSE                                                                      
         WIG004=SQRT(WIG004/NFI004 )                                            
         SIG004=SQRT(SIG004/NFI004 )                                            
      ENDIF                                                                     
      IF(NFI080 .EQ.0) THEN                                                     
         WIG080=0.d0                                                            
         SIG080=0.d0                                                            
      ELSE                                                                      
         WIG080=SQRT(WIG080/NFI080 )                                            
         SIG080=SQRT(SIG080/NFI080 )                                            
      ENDIF                                                                     
      IF(NFI010 .EQ.0) THEN                                                     
         WIG010=0.d0                                                            
         SIG010=0.d0                                                            
      ELSE                                                                      
         WIG010=SQRT(WIG010/NFI010 )                                            
         SIG010=SQRT(SIG010/NFI010 )                                            
      ENDIF                                                                     
      IF(NFI100 .EQ.0) THEN                                                     
         WIG100=0.d0                                                            
         SIG100=0.d0                                                            
      ELSE                                                                      
         WIG100=SQRT(WIG100/NFI100 )                                            
         SIG100=SQRT(SIG100/NFI100 )                                            
      ENDIF                                                                     
      IF(NFI020 .EQ.0) THEN                                                     
         WIG020=0.d0                                                            
         SIG020=0.d0                                                            
      ELSE                                                                      
         WIG020=SQRT(WIG020/NFI020 )                                            
         SIG020=SQRT(SIG020/NFI020 )                                            
      ENDIF                                                                     
      IF(NFI200 .EQ.0) THEN                                                     
         WIG200=0.d0                                                            
         SIG200=0.d0                                                            
      ELSE                                                                      
         WIG200=SQRT(WIG200/NFI200 )                                            
         SIG200=SQRT(SIG200/NFI200 )                                            
      ENDIF                                                                     
      IF(NFI045 .EQ.0) THEN                                                     
         WIG045=0.d0                                                            
         SIG045=0.d0                                                            
      ELSE                                                                      
         WIG045=SQRT(WIG045/NFI045 )                                            
         SIG045=SQRT(SIG045/NFI045 )                                            
      ENDIF                                                                     
      IF(NFI1M0 .EQ.0) THEN                                                     
         WIG1M0=0.d0                                                            
         SIG1M0=0.d0                                                            
      ELSE                                                                      
         WIG1M0=SQRT(WIG1M0/NFI1M0 )                                            
         SIG1M0=SQRT(SIG1M0/NFI1M0 )                                            
      ENDIF                                                                     
      IF(NGST0.EQ.0) THEN                                                       
         WGST0=0.d0                                                             
         SGST0=0.d0                                                             
      ELSE                                                                      
         WGST0=SQRT(WGST0/NGST0)                                                
         SGST0=SQRT(SGST0/NGST0)                                                
      ENDIF                                                                     

      IF (NFT00.eq.0) THEN
         WFT00 = 0.d0
         SFT00 = 0.d0
      ELSE
         WFT00 = SQRT(WFT00 / NFT00)
         SFT00 = SQRT(SFT00 / NFT00)
      ENDIF

      IF (NFT01.eq.0) THEN
         WFT01 = 0.d0
         SFT01 = 0.d0
      ELSE
         WFT01 = SQRT(WFT01 / NFT01)
         SFT01 = SQRT(SFT01 / NFT01)
      ENDIF

      IF (NFT02.eq.0) THEN
         WFT02 = 0.d0
         SFT02 = 0.d0
      ELSE
         WFT02 = SQRT(WFT02 / NFT02)
         SFT02 = SQRT(SFT02 / NFT02)
      ENDIF

      IF (NFT11.eq.0) THEN
         WFT11 = 0.d0
         SFT11 = 0.d0
      ELSE
         WFT11 = SQRT(WFT11 / NFT11)
         SFT11 = SQRT(SFT11 / NFT11)
      ENDIF

      IF (NFT12.eq.0) THEN
         WFT12 = 0.d0
         SFT12 = 0.d0
      ELSE
         WFT12 = SQRT(WFT12 / NFT12)
         SFT12 = SQRT(SFT12 / NFT12)
      ENDIF

      IF (NFT22.eq.0) THEN
         WFT22 = 0.d0
         SFT22 = 0.d0
      ELSE
         WFT22 = SQRT(WFT22 / NFT22)
         SFT22 = SQRT(SFT22 / NFT22)
      ENDIF

      IF (NFT00E.eq.0) THEN
         WFT00E = 0.d0
      ELSE
         WFT00E = SQRT(WFT00E / NFT00E)
      ENDIF

      IF (NFT01E.eq.0) THEN
         WFT01E = 0.d0
      ELSE
         WFT01E = SQRT(WFT01E / NFT01E)
      ENDIF

      IF (NFT02E.eq.0) THEN
         WFT02E = 0.d0
      ELSE
         WFT02E = SQRT(WFT02E / NFT02E)
      ENDIF

      IF (NFT11E.eq.0) THEN
         WFT11E = 0.d0
      ELSE
         WFT11E = SQRT(WFT11E / NFT11E)
      ENDIF

      IF (NFT12E.eq.0) THEN
         WFT12E = 0.d0
      ELSE
         WFT12E = SQRT(WFT12E / NFT12E)
      ENDIF

      IF (NFT22E.eq.0) THEN
         WFT22E = 0.d0
      ELSE
         WFT22E = SQRT(WFT22E / NFT22E)
      ENDIF

      IF (NFT00A.eq.0) THEN
         WFT00A = 0.d0
      ELSE
         WFT00A = SQRT(WFT00A / NFT00A)
      ENDIF

      IF (NFT01A.eq.0) THEN
         WFT01A = 0.d0
      ELSE
         WFT01A = SQRT(WFT01A / NFT01A)
      ENDIF

      IF (NFT02A.eq.0) THEN
         WFT02A = 0.d0
      ELSE
         WFT02A = SQRT(WFT02A / NFT02A)
      ENDIF

      IF (NFT11A.eq.0) THEN
         WFT11A = 0.d0
      ELSE
         WFT11A = SQRT(WFT11A / NFT11A)
      ENDIF

      IF (NFT12A.eq.0) THEN
         WFT12A = 0.d0
      ELSE
         WFT12A = SQRT(WFT12A / NFT12A)
      ENDIF

      IF (NFT22A.eq.0) THEN
         WFT22A = 0.d0
      ELSE
         WFT22A = SQRT(WFT22A / NFT22A)
      ENDIF

      NFT = NFT00 + NFT01 + NFT02 + NFT11 + NFT12 + NFT22

      IF ( NFT .eq. 0 ) then
         SFT = 0.d0
         WFT = 0.d0
      ELSE
         SFT = SQRT(SFT / NFT)
         WFT = SQRT(WFT / NFT)
      ENDIF

      IF ( NFTW1 .eq. 0 ) then
         SFT1 = 0.d0
      ELSE
         SFT1 = SQRT( SFT1 / NFTW1 )
      END IF

      IF ( NFTW2 .eq. 0 ) then
         SFT2 = 0.d0
      ELSE
         SFT2 = SQRT( SFT2 / NFTW2 )
      END IF

      IF ( NFTW3 .eq. 0 ) then
         SFT3 = 0.d0
      ELSE
         SFT3 = SQRT( SFT3 / NFTW3 )
      END IF

      IF ( NFTW4 .eq. 0 ) then
         SFT4 = 0.d0
      ELSE
         SFT4 = SQRT( SFT4 / NFTW4 )
      END IF

      IF ( NFTW1A .eq. 0 ) then
         SFT1A = 0.d0
      ELSE
         SFT1A = SQRT( SFT1A / NFTW1A )
      END IF

      IF ( NFTW2A .eq. 0 ) then
         SFT2A = 0.d0
      ELSE
         SFT2A = SQRT( SFT2A / NFTW2A )
      END IF

      IF ( NFTW3A .eq. 0 ) then
         SFT3A = 0.d0
      ELSE
         SFT3A = SQRT( SFT3A / NFTW3A )
      END IF

      IF ( NFTW4A .eq. 0 ) then
         SFT4A = 0.d0
      ELSE
         SFT4A = SQRT( SFT4A / NFTW4A )
      END IF


      IF ( NFTW1E .eq. 0 ) then
         SFT1E = 0.d0
      ELSE
         SFT1E = SQRT( SFT1E / NFTW1E )
      END IF

      IF ( NFTW2E .eq. 0 ) then
         SFT2E = 0.d0
      ELSE
         SFT2E = SQRT( SFT2E / NFTW2E )
      END IF

      IF ( NFTW3E .eq. 0 ) then
         SFT3E = 0.d0
      ELSE
         SFT3E = SQRT( SFT3E / NFTW3E )
      END IF

      IF ( NFTW4E .eq. 0 ) then
         SFT4E = 0.d0
      ELSE
         SFT4E = SQRT( SFT4E / NFTW4E )
      END IF







      IF (NMW00.eq.0) THEN
         WMW00 = 0.d0
         SMW00 = 0.d0
      ELSE
         WMW00 = SQRT(WMW00 / NMW00)
         SMW00 = SQRT(SMW00 / NMW00)
      ENDIF

      IF (NMW01.eq.0) THEN
         WMW01 = 0.d0
         SMW01 = 0.d0
      ELSE
         WMW01 = SQRT(WMW01 / NMW01)
         SMW01 = SQRT(SMW01 / NMW01)
      ENDIF

      IF (NMW02.eq.0) THEN
         WMW02 = 0.d0
         SMW02 = 0.d0
      ELSE
         WMW02 = SQRT(WMW02 / NMW02)
         SMW02 = SQRT(SMW02 / NMW02)
      ENDIF

      IF (NMW11.eq.0) THEN
         WMW11 = 0.d0
         SMW11 = 0.d0
      ELSE
         WMW11 = SQRT(WMW11 / NMW11)
         SMW11 = SQRT(SMW11 / NMW11)
      ENDIF

      IF (NMW12.eq.0) THEN
         WMW12 = 0.d0
         SMW12 = 0.d0
      ELSE
         WMW12 = SQRT(WMW12 / NMW12)
         SMW12 = SQRT(SMW12 / NMW12)
      ENDIF

      IF (NMW22.eq.0) THEN
         WMW22 = 0.d0
         SMW22 = 0.d0
      ELSE
         WMW22 = SQRT(WMW22 / NMW22)
         SMW22 = SQRT(SMW22 / NMW22)
      ENDIF

      IF (NMW00E.eq.0) THEN
         WMW00E = 0.d0
      ELSE
         WMW00E = SQRT(WMW00E / NMW00E)
      ENDIF

      IF (NMW01E.eq.0) THEN
         WMW01E = 0.d0
      ELSE
         WMW01E = SQRT(WMW01E / NMW01E)
      ENDIF

      IF (NMW02E.eq.0) THEN
         WMW02E = 0.d0
      ELSE
         WMW02E = SQRT(WMW02E / NMW02E)
      ENDIF

      IF (NMW11E.eq.0) THEN
         WMW11E = 0.d0
      ELSE
         WMW11E = SQRT(WMW11E / NMW11E)
      ENDIF

      IF (NMW12E.eq.0) THEN
         WMW12E = 0.d0
      ELSE
         WMW12E = SQRT(WMW12E / NMW12E)
      ENDIF

      IF (NMW22E.eq.0) THEN
         WMW22E = 0.d0
      ELSE
         WMW22E = SQRT(WMW22E / NMW22E)
      ENDIF

      IF (NMW00A.eq.0) THEN
         WMW00A = 0.d0
      ELSE
         WMW00A = SQRT(WMW00A / NMW00A)
      ENDIF

      IF (NMW01A.eq.0) THEN
         WMW01A = 0.d0
      ELSE
         WMW01A = SQRT(WMW01A / NMW01A)
      ENDIF

      IF (NMW02A.eq.0) THEN
         WMW02A = 0.d0
      ELSE
         WMW02A = SQRT(WMW02A / NMW02A)
      ENDIF

      IF (NMW11A.eq.0) THEN
         WMW11A = 0.d0
      ELSE
         WMW11A = SQRT(WMW11A / NMW11A)
      ENDIF

      IF (NMW12A.eq.0) THEN
         WMW12A = 0.d0
      ELSE
         WMW12A = SQRT(WMW12A / NMW12A)
      ENDIF

      IF (NMW22A.eq.0) THEN
         WMW22A = 0.d0
      ELSE
         WMW22A = SQRT(WMW22A / NMW22A)
      ENDIF

      NMW = NMW00 + NMW01 + NMW02 + NMW11 + NMW12 + NMW22

      IF ( NMW .eq. 0 ) then
         SMWJ = 0.d0
         WMWJ = 0.d0
      ELSE
         SMWJ = SQRT(SMWJ / NMW)
         WMWJ = SQRT(WMWJ / NMW)
      ENDIF



      IF(NFI004A .EQ.0) THEN                                                     
         SIG004A=0.d0                                                            
      ELSE                                                                      
         SIG004A=SQRT(SIG004A/NFI004A )                                            
      ENDIF                                                                     
      IF(NFI080A .EQ.0) THEN                                                     
         SIG080A=0.d0                                                            
      ELSE                                                                      
         SIG080A=SQRT(SIG080A/NFI080A )                                            
      ENDIF                                                                     
      IF(NFI010A .EQ.0) THEN                                                     
         SIG010A=0.d0                                                            
      ELSE                                                                      
         SIG010A=SQRT(SIG010A/NFI010A )                                            
      ENDIF                                                                     
      IF(NFI100A .EQ.0) THEN                                                     
         SIG100A=0.d0                                                            
      ELSE                                                                      
         SIG100A=SQRT(SIG100A/NFI100A )                                            
      ENDIF                                                                     
      IF(NFI020A .EQ.0) THEN                                                     
         SIG020A=0.d0                                                            
      ELSE                                                                      
         SIG020A=SQRT(SIG020A/NFI020A )                                            
      ENDIF                                                                     
      IF(NFI200A .EQ.0) THEN                                                     
         SIG200A=0.d0                                                            
      ELSE                                                                      
         SIG200A=SQRT(SIG200A/NFI200A )                                            
      ENDIF                                                                     
      IF(NFI045A .EQ.0) THEN                                                     
         SIG045A=0.d0                                                            
      ELSE                                                                      
         SIG045A=SQRT(SIG045A/NFI045A )                                            
      ENDIF                                                                     
      IF(NFI1M0A .EQ.0) THEN                                                     
         SIG1M0A=0.d0                                                            
      ELSE                                                                      
         SIG1M0A=SQRT(SIG1M0A/NFI1M0A )                                            
      ENDIF                                                                     

      IF(NFI004E .EQ.0) THEN                                                     
         SIG004E=0.d0                                                            
      ELSE                                                                      
         SIG004E=SQRT(SIG004E/NFI004E )                                            
      ENDIF                                                                     
      IF(NFI080E .EQ.0) THEN                                                     
         SIG080E=0.d0                                                            
      ELSE                                                                      
         SIG080E=SQRT(SIG080E/NFI080E )                                            
      ENDIF                                                                     
      IF(NFI010E .EQ.0) THEN                                                     
         SIG010E=0.d0                                                            
      ELSE                                                                      
         SIG010E=SQRT(SIG010E/NFI010E )                                            
      ENDIF                                                                     
      IF(NFI100E .EQ.0) THEN                                                     
         SIG100E=0.d0                                                            
      ELSE                                                                      
         SIG100E=SQRT(SIG100E/NFI100E )                                            
      ENDIF                                                                     
      IF(NFI020E .EQ.0) THEN                                                     
         SIG020E=0.d0                                                            
      ELSE                                                                      
         SIG020E=SQRT(SIG020E/NFI020E )                                            
      ENDIF                                                                     
      IF(NFI200E .EQ.0) THEN                                                     
         SIG200E=0.d0                                                            
      ELSE                                                                      
         SIG200E=SQRT(SIG200E/NFI200E )                                            
      ENDIF                                                                     
      IF(NFI045E .EQ.0) THEN                                                     
         SIG045E=0.d0                                                            
      ELSE                                                                      
         SIG045E=SQRT(SIG045E/NFI045E )                                            
      ENDIF                                                                     
      IF(NFI1M0E .EQ.0) THEN                                                     
         SIG1M0E=0.d0                                                            
      ELSE                                                                      
         SIG1M0E=SQRT(SIG1M0E/NFI1M0E )                                            
      ENDIF                                                                     

      IF(NFI070A .EQ.0) THEN                                                     
         SW070A=0.d0                                                             
      ELSE                                                                      
         SW070A=SQRT(SW070A/NFI070A )                                              
      ENDIF                                                                     

      IF(NFI070E .EQ.0) THEN                                                     
         SW070E=0.d0                                                             
      ELSE                                                                      
         SW070E=SQRT(SW070E/NFI070E )                                              
      ENDIF                                                                     


      WRITE(6,900) S,NOFIT                                                      
      WRITE(6,901) WMW,NFITMW                                                   
      SMW0=SMW0*29979.2458d0                                                    
      SMW1=SMW1*29979.2458d0                                                    
      SMW2=SMW2*29979.2458d0                                                    
      SMW3=SMW3*29979.2458d0                                                    
      WRITE(6,903) WMW0,NFIMW0                                                  
      WRITE(6,904) WMW1,NFIMW1                                                  
      WRITE(6,700) WMW2,NFIMW2                                                  
      WRITE(6,600) WMW3,NFIMW3                                                  
      WRITE(6,902) WIR,NFITIR                                                   
      WRITE(6,912) WIR10,NFIR10                                                 
      WRITE(6,913) WIR21,NFIR21                                                 
      WRITE(6,899) WGST0,NGST0                                                  
      WRITE(6,920) WIG004,NFI004                                                
      WRITE(6,931) WIG010,NFI010                                                
      WRITE(6,933) WIG020,NFI020                                                
      WRITE(6,935) WIG045,NFI045                                                
      WRITE(6,936) WIG070,NFI070                                                
      WRITE(6,930) WIG080,NFI080                                                
      WRITE(6,932) WIG100,NFI100                                                
      WRITE(6,934) WIG200,NFI200                                                
      WRITE(6,937) WIG1M0,NFI1M0                                                


C     SIR=SIR*0.000353553d0                                                     
C     SMW0=SMW0*0.078d0                                                         
C     SMW1=SMW1*0.04d0                                                          
      CCC=29979.2458d0                                                          
      SIG004=SIG004*CCC                                                         
      SIG080=SIG080*CCC                                                         
      SIG010=SIG010*CCC                                                         
      SIG100=SIG100*CCC                                                         
      SIG020=SIG020*CCC                                                         
      SIG200=SIG200*CCC                                                         
      SIG045=SIG045*CCC                                                         
      SIG1M0=SIG1M0*CCC                                                         
      SW070=SW070*CCC                                                           
C     SGST0=SGST0*0.0005d0                                                      

      SMW00 = SMW00 * CCC
      SMW01 = SMW01 * CCC
      SMW02 = SMW02 * CCC
      SMW11 = SMW11 * CCC
      SMW12 = SMW12 * CCC
      SMW22 = SMW22 * CCC

      SIG004A=SIG004A*CCC                                                         
      SIG080A=SIG080A*CCC                                                         
      SIG010A=SIG010A*CCC                                                         
      SIG100A=SIG100A*CCC                                                         
      SIG020A=SIG020A*CCC                                                         
      SIG200A=SIG200A*CCC                                                         
      SIG045A=SIG045A*CCC                                                         
      SIG1M0A=SIG1M0A*CCC                                                         
      SW070A=SW070A*CCC                                                           

      SIG004E=SIG004E*CCC                                                         
      SIG080E=SIG080E*CCC                                                         
      SIG010E=SIG010E*CCC                                                         
      SIG100E=SIG100E*CCC                                                         
      SIG020E=SIG020E*CCC                                                         
      SIG200E=SIG200E*CCC                                                         
      SIG045E=SIG045E*CCC                                                         
      SIG1M0E=SIG1M0E*CCC                                                         
      SW070E=SW070E*CCC                                                           


      write(6,*) '*** New Summary ***'

      WRITE(6,1000) S,NOFIT 
      write(6, *)
                                                     
C      WRITE(6,938) 0.0, NFT1
C      WRITE(6,939) 0.0, NFT2
C      WRITE(6,940) 0.0, NFT1A
C      WRITE(6,941) 0.0, NFT2A
C      WRITE(6,942) 0.0, NFT1E
C      WRITE(6,943) 0.0, NFT2E

      write(6, 1001) WFT00, NFT00
      write(6, 1002) WFT01, NFT01
      write(6, 1003) WFT02, NFT02
      write(6, 1004) WFT11, NFT11
      write(6, 1005) WFT12, NFT12
      write(6, 1006) WFT22, NFT22
      write(6, *)

      write(6, 1011) SFT00, NFT00
      write(6, 1012) SFT01, NFT01
      write(6, 1013) SFT02, NFT02
      write(6, 1014) SFT11, NFT11
      write(6, 1015) SFT12, NFT12
      write(6, 1016) SFT22, NFT22
      write(6, *)

      write(6, 1021) WFT00A, NFT00A
      write(6, 1022) WFT01A, NFT01A
      write(6, 1023) WFT02A, NFT02A
      write(6, 1024) WFT11A, NFT11A
      write(6, 1025) WFT12A, NFT12A
      write(6, 1026) WFT22A, NFT22A
      write(6, *)

      write(6, 1031) WFT00E, NFT00E
      write(6, 1032) WFT01E, NFT01E
      write(6, 1033) WFT02E, NFT02E
      write(6, 1034) WFT11E, NFT11E
      write(6, 1035) WFT12E, NFT12E
      write(6, 1036) WFT22E, NFT22E
      write(6, *)

      write(6, 1041) WFT, NFT
      write(6, 1042) SFT, NFT
      write(6, *)

      write(6, 1051) SFT1, NFTW1
      write(6, 1052) SFT2, NFTW2
      write(6, 1053) SFT3, NFTW3
      write(6, 1054) SFT4, NFTW4
      write(6, *)

      write(6, 1061) SFT1A, NFTW1A
      write(6, 1062) SFT2A, NFTW2A
      write(6, 1063) SFT3A, NFTW3A
      write(6, 1064) SFT4A, NFTW4A
      write(6, *)

      write(6, 1071) SFT1E, NFTW1E
      write(6, 1072) SFT2E, NFTW2E
      write(6, 1073) SFT3E, NFTW3E
      write(6, 1074) SFT4E, NFTW4E
      write(6, *)

      write(6, 1101) WMW00, NMW00
      write(6, 1102) WMW01, NMW01
      write(6, 1103) WMW02, NMW02
      write(6, 1104) WMW11, NMW11
      write(6, 1105) WMW12, NMW12
      write(6, 1106) WMW22, NMW22
      write(6, *)

      write(6, 1111) SMW00, NMW00
      write(6, 1112) SMW01, NMW01
      write(6, 1113) SMW02, NMW02
      write(6, 1114) SMW11, NMW11
      write(6, 1115) SMW12, NMW12
      write(6, 1116) SMW22, NMW22
      write(6, *)

      write(6, 1121) WMW00A, NMW00A
      write(6, 1122) WMW01A, NMW01A
      write(6, 1123) WMW02A, NMW02A
      write(6, 1124) WMW11A, NMW11A
      write(6, 1125) WMW12A, NMW12A
      write(6, 1126) WMW22A, NMW22A
      write(6, *)

      write(6, 1131) WMW00E, NMW00E
      write(6, 1132) WMW01E, NMW01E
      write(6, 1133) WMW02E, NMW02E
      write(6, 1134) WMW11E, NMW11E
      write(6, 1135) WMW12E, NMW12E
      write(6, 1136) WMW22E, NMW22E
      write(6, *)

      write(6, 1151) SIG004, NFI004
      write(6, 1152) SIG010, NFI010
      write(6, 1153) SIG020, NFI020
      write(6, 1154) SIG045, NFI045
      write(6, 1155) SW070,  NFI070
      write(6, 1156) SIG080, NFI080
      write(6, 1157) SIG100, NFI100
      write(6, 1158) SIG200, NFI200
      write(6, 1159) SIG1M0, NFI1M0
      write(6, *)

      write(6, 1161) SIG004A, NFI004A
      write(6, 1162) SIG010A, NFI010A
      write(6, 1163) SIG020A, NFI020A
      write(6, 1164) SIG045A, NFI045A
      write(6, 1165) SW070A,  NFI070A
      write(6, 1166) SIG080A, NFI080A
      write(6, 1167) SIG100A, NFI100A
      write(6, 1168) SIG200A, NFI200A
      write(6, 1169) SIG1M0A, NFI1M0A
      write(6, *)

      write(6, 1171) SIG004E, NFI004E
      write(6, 1172) SIG010E, NFI010E
      write(6, 1173) SIG020E, NFI020E
      write(6, 1174) SIG045E, NFI045E
      write(6, 1175) SW070E,  NFI070E
      write(6, 1176) SIG080E, NFI080E
      write(6, 1177) SIG100E, NFI100E
      write(6, 1178) SIG200E, NFI200E
      write(6, 1179) SIG1M0E, NFI1M0E
      write(6, *)


      write(6, 1141) WMWJ, NMW
      write(6, 1142) SMWJ*CCC, NMW

      write(6,*) '*** End New Summary ***'


      WRITE(6,906) SMW0                                                         
      WRITE(6,907) SMW1                                                         
      WRITE(6,701) SMW2                                                         
      WRITE(6,601) SMW3                                                         
      WRITE(6,905) SIR                                                          
      WRITE(6,915) SIR10                                                        
      WRITE(6,916) SIR21                                                        
      WRITE(6,888) SGST0                                                        
      WRITE(6,890) SIG004                                                       
      WRITE(6,892) SIG010                                                       
      WRITE(6,894) SIG020                                                       
      WRITE(6,896) SIG045                                                       
      WRITE(6,897) SW070                                                        
      WRITE(6,891) SIG080                                                       
      WRITE(6,893) SIG100                                                       
      WRITE(6,895) SIG200                                                       
      WRITE(6,898) SIG1M0                                                       

 900  FORMAT(/,' STD.DEV.(UNITLESS)  = ',F20.8,' NB.DATA= ',I6)                 
 901  FORMAT(/,' RMS.DEV.MW(UNITLESS)     = ',F20.4,2X,I4)                      
 902  FORMAT(/,' RMS.DEV.IR(UNITLESS)       = ',F20.6,2X,I6)                    
 912  FORMAT(/,' RMS.DEV.IR(UNITLESS) V=1-0 = ',F20.6,2X,I6)                    
 913  FORMAT(/,' RMS.DEV.IR(UNITLESS) V=2-1 = ',F20.6,2X,I6)                    
 903  FORMAT(/,' RMS.DEV.MW(UNITLESS)VT=0 = ',F20.4,2X,I4)                      
 904  FORMAT(/,' RMS.DEV.MW(UNITLESS)VT=1 = ',F20.4,2X,I4)                      
 700  FORMAT(/,' RMS.DEV.MW(UNITLESS)VT=2 = ',F20.4,2X,I4)                      
 600  FORMAT(/,' RMS.DEV.MW(UNITLESS)VT=3 = ',F20.4,2X,I4)                      
 899  FORMAT(/,' RMS.DEV.GROUND STATE DIFF. VT=0) = ',F20.4,2X,I4)              
 920  FORMAT(/,' STD.MW. W=0.004MHZ',F20.6,2X,I4)                               
 930  FORMAT(/,' STD.MW. W=0.080MHZ',F20.6,2X,I4)                               
 931  FORMAT(/,' STD.MW. W=0.010MHZ',F20.6,2X,I4)                               
 932  FORMAT(/,' RMS.MW. W=0.100MHZ',F20.6,2X,I4)                               
 933  FORMAT(/,' RMS.MW. W=0.020MHZ',F20.6,2X,I4)                               
 934  FORMAT(/,' RMS.MW. W=0.200MHZ',F20.6,2X,I4)                               
 935  FORMAT(/,' RMS.MW. W=0.045MHZ',F20.6,2X,I4)                               
 936  FORMAT(/,' RMS.MW. W=0.070MHZ',F20.6,2X,I4)                               
 937  FORMAT(/,' RMS.MW. W=1.000MHZ',F20.6,2X,I4)                               


 905  FORMAT(/,' RMS.DEV.IR(CM-1)       = ',F20.6)                              
 915  FORMAT(/,' RMS.DEV.IR(CM-1) V=1-0 = ',F20.6)                              
 916  FORMAT(/,' RMS.DEV.IR(CM-1) V=2-1 = ',F20.6)                              
 906  FORMAT(/,' RMS.DEV.MICRONDE VT=0(MHZ)= ',F20.4)                           
 907  FORMAT(/,' RMS.DEV.MICRONDE VT=1(MHZ)= ',F20.4)                           
 701  FORMAT(/,' RMS.DEV.MICRONDE VT=2(MHZ)= ',F20.4)                           
 601  FORMAT(/,' RMS.DEV.MICRONDE VT=3(MHZ)= ',F20.4)                           
 888  FORMAT(/,' RMS.DEV.GROUND STATE DIFF.VT=0 (CM-1)) = ',F20.6)              
 890  FORMAT(/,' RMS.DEV.MICROWAVE,WEIGHT=0.004MHZ',F20.4,' MHZ')               
 891  FORMAT(/,' RMS.DEV.MICROWAVE,WEIGHT=0.080MHZ',F20.4,' MHZ')               
 892  FORMAT(/,' RMS.DEV.MICROWAVE,WEIGHT=0.010MHZ',F20.4,' MHZ')               
 893  FORMAT(/,' RMS.DEV.MICROWAVE,WEIGHT=0.100MHZ',F20.4,' MHZ')               
 894  FORMAT(/,' RMS.DEV.MICROWAVE,WEIGHT=0.020MHZ',F20.4,' MHZ')               
 895  FORMAT(/,' RMS.DEV.MICROWAVE,WEIGHT=0.200MHZ',F20.4,' MHZ')               
 896  FORMAT(/,' RMS.DEV.MICROWAVE,WEIGHT=0.045MHZ',F20.4,' MHZ')               
 897  FORMAT(/,' RMS.DEV.MICROWAVE,WEIGHT=0.070MHZ',F20.4,' MHZ')               
 898  FORMAT(/,' RMS.DEV.MICROWAVE,WEIGHT=1.000MHZ',F20.4,' MHZ')               


 938  FORMAT(' RMS.IR. W=0.0002 CM-1             ', F20.6, 2X, I8)                               
 939  FORMAT(' RMS.IR. W=0.0004 CM-1             ', F20.6, 2X, I8)                               

 940  FORMAT(' RMS.IR. A W=0.0002 CM-1           ', F20.6, 2X, I8)                               
 941  FORMAT(' RMS.IR. A W=0.0004 CM-1           ', F20.6, 2X, I8)                               

 942  FORMAT(' RMS.IR. E W=0.0002 CM-1           ', F20.6, 2X, I8)                               
 943  FORMAT(' RMS.IR. E W=0.0004 CM-1           ', F20.6, 2X, I8)                               

 1000 FORMAT(' STD.DEV.(UNITLESS)                ', F20.6, 2X, I8)                 

 1001 FORMAT(' RMS.DEV.IR(UNITLESS) VT=0-0       ', F20.6, 2X, I8)
 1002 FORMAT(' RMS.DEV.IR(UNITLESS) VT=1-0       ', F20.6, 2X, I8)
 1003 FORMAT(' RMS.DEV.IR(UNITLESS) VT=2-0       ', F20.6, 2X, I8)
 1004 FORMAT(' RMS.DEV.IR(UNITLESS) VT=1-1       ', F20.6, 2X, I8)
 1005 FORMAT(' RMS.DEV.IR(UNITLESS) VT=2-1       ', F20.6, 2X, I8)
 1006 FORMAT(' RMS.DEV.IR(UNITLESS) VT=2-2       ', F20.6, 2X, I8)

 1011 FORMAT(' RMS.DEV.IR(CM-1) VT=0-0           ', F20.6, 2X, I8)
 1012 FORMAT(' RMS.DEV.IR(CM-1) VT=1-0           ', F20.6, 2X, I8)
 1013 FORMAT(' RMS.DEV.IR(CM-1) VT=2-0           ', F20.6, 2X, I8)
 1014 FORMAT(' RMS.DEV.IR(CM-1) VT=1-1           ', F20.6, 2X, I8)
 1015 FORMAT(' RMS.DEV.IR(CM-1) VT=2-1           ', F20.6, 2X, I8)
 1016 FORMAT(' RMS.DEV.IR(CM-1) VT=2-2           ', F20.6, 2X, I8)

 1021 FORMAT('     RMS.DEV.IR(UNITLESS) VT=0-0 A ', F20.6, 2X, I8)
 1022 FORMAT('     RMS.DEV.IR(UNITLESS) VT=1-0 A ', F20.6, 2X, I8)
 1023 FORMAT('     RMS.DEV.IR(UNITLESS) VT=2-0 A ', F20.6, 2X, I8)
 1024 FORMAT('     RMS.DEV.IR(UNITLESS) VT=1-1 A ', F20.6, 2X, I8)
 1025 FORMAT('     RMS.DEV.IR(UNITLESS) VT=2-1 A ', F20.6, 2X, I8)
 1026 FORMAT('     RMS.DEV.IR(UNITLESS) VT=2-2 A ', F20.6, 2X, I8)

 1031 FORMAT('     RMS.DEV.IR(UNITLESS) VT=0-0 E ', F20.6, 2X, I8)
 1032 FORMAT('     RMS.DEV.IR(UNITLESS) VT=1-0 E ', F20.6, 2X, I8)
 1033 FORMAT('     RMS.DEV.IR(UNITLESS) VT=2-0 E ', F20.6, 2X, I8)
 1034 FORMAT('     RMS.DEV.IR(UNITLESS) VT=1-1 E ', F20.6, 2X, I8)
 1035 FORMAT('     RMS.DEV.IR(UNITLESS) VT=2-1 E ', F20.6, 2X, I8)
 1036 FORMAT('     RMS.DEV.IR(UNITLESS) VT=2-2 E ', F20.6, 2X, I8)

 1041 FORMAT(' RMS.DEV.IR(UNITLESS)              ', F20.6, 2X, I8)
 1042 FORMAT(' RMS.DEV.IR(CM-1)                  ', F20.6, 2X, I8)

 1051 FORMAT('   RMS.DEV.WEIGHT=0.00020 CM-1     ', F20.6, 2X, I8)
 1052 FORMAT('   RMS.DEV.WEIGHT=0.00035 CM-1     ', F20.6, 2X, I8)
 1053 FORMAT('   RMS.DEV.WEIGHT=0.00050 CM-1     ', F20.6, 2X, I8)
 1054 FORMAT('   RMS.DEV.WEIGHT=0.00060 CM-1     ', F20.6, 2X, I8)

 1061 FORMAT('     RMS.DEV.WEIGHT=0.00020 CM-1 A ', F20.6, 2X, I8)
 1062 FORMAT('     RMS.DEV.WEIGHT=0.00035 CM-1 A ', F20.6, 2X, I8)
 1063 FORMAT('     RMS.DEV.WEIGHT=0.00050 CM-1 A ', F20.6, 2X, I8)
 1064 FORMAT('     RMS.DEV.WEIGHT=0.00060 CM-1 A ', F20.6, 2X, I8)

 1071 FORMAT('     RMS.DEV.WEIGHT=0.00020 CM-1 E ', F20.6, 2X, I8)
 1072 FORMAT('     RMS.DEV.WEIGHT=0.00035 CM-1 E ', F20.6, 2X, I8)
 1073 FORMAT('     RMS.DEV.WEIGHT=0.00050 CM-1 E ', F20.6, 2X, I8)
 1074 FORMAT('     RMS.DEV.WEIGHT=0.00060 CM-1 E ', F20.6, 2X, I8)


 1101 FORMAT(' RMS.DEV.MW(UNITLESS) VT=0-0       ', F20.6, 2X, I8)
 1102 FORMAT(' RMS.DEV.MW(UNITLESS) VT=1-0       ', F20.6, 2X, I8)
 1103 FORMAT(' RMS.DEV.MW(UNITLESS) VT=2-0       ', F20.6, 2X, I8)
 1104 FORMAT(' RMS.DEV.MW(UNITLESS) VT=1-1       ', F20.6, 2X, I8)
 1105 FORMAT(' RMS.DEV.MW(UNITLESS) VT=2-1       ', F20.6, 2X, I8)
 1106 FORMAT(' RMS.DEV.MW(UNITLESS) VT=2-2       ', F20.6, 2X, I8)

 1111 FORMAT(' RMS.DEV.MW(MHZ ) VT=0-0           ', F20.6, 2X, I8)
 1112 FORMAT(' RMS.DEV.MW(MHZ ) VT=1-0           ', F20.6, 2X, I8)
 1113 FORMAT(' RMS.DEV.MW(MHZ ) VT=2-0           ', F20.6, 2X, I8)
 1114 FORMAT(' RMS.DEV.MW(MHZ ) VT=1-1           ', F20.6, 2X, I8)
 1115 FORMAT(' RMS.DEV.MW(MHZ ) VT=2-1           ', F20.6, 2X, I8)
 1116 FORMAT(' RMS.DEV.MW(MHZ ) VT=2-2           ', F20.6, 2X, I8)

 1121 FORMAT('     RMS.DEV.MW(UNITLESS) VT=0-0 A ', F20.6, 2X, I8)
 1122 FORMAT('     RMS.DEV.MW(UNITLESS) VT=1-0 A ', F20.6, 2X, I8)
 1123 FORMAT('     RMS.DEV.MW(UNITLESS) VT=2-0 A ', F20.6, 2X, I8)
 1124 FORMAT('     RMS.DEV.MW(UNITLESS) VT=1-1 A ', F20.6, 2X, I8)
 1125 FORMAT('     RMS.DEV.MW(UNITLESS) VT=2-1 A ', F20.6, 2X, I8)
 1126 FORMAT('     RMS.DEV.MW(UNITLESS) VT=2-2 A ', F20.6, 2X, I8)

 1131 FORMAT('     RMS.DEV.MW(UNITLESS) VT=0-0 E ', F20.6, 2X, I8)
 1132 FORMAT('     RMS.DEV.MW(UNITLESS) VT=1-0 E ', F20.6, 2X, I8)
 1133 FORMAT('     RMS.DEV.MW(UNITLESS) VT=2-0 E ', F20.6, 2X, I8)
 1134 FORMAT('     RMS.DEV.MW(UNITLESS) VT=1-1 E ', F20.6, 2X, I8)
 1135 FORMAT('     RMS.DEV.MW(UNITLESS) VT=2-1 E ', F20.6, 2X, I8)
 1136 FORMAT('     RMS.DEV.MW(UNITLESS) VT=2-2 E ', F20.6, 2X, I8)

 1141 FORMAT(' RMS.DEV.MW(UNITLESS)              ', F20.6, 2X, I8)
 1142 FORMAT(' RMS.DEV.MW(MHZ )                  ', F20.6, 2X, I8)

 1151 FORMAT('   RMS.DEV.WEIGHT=0.004   MHZ      ', F20.6, 2X, I8)
 1152 FORMAT('   RMS.DEV.WEIGHT=0.010   MHZ      ', F20.6, 2X, I8)
 1153 FORMAT('   RMS.DEV.WEIGHT=0.020   MHZ      ', F20.6, 2X, I8)
 1154 FORMAT('   RMS.DEV.WEIGHT=0.045   MHZ      ', F20.6, 2X, I8)
 1155 FORMAT('   RMS.DEV.WEIGHT=0.070   MHZ      ', F20.6, 2X, I8)
 1156 FORMAT('   RMS.DEV.WEIGHT=0.080   MHZ      ', F20.6, 2X, I8)
 1157 FORMAT('   RMS.DEV.WEIGHT=0.100   MHZ      ', F20.6, 2X, I8)
 1158 FORMAT('   RMS.DEV.WEIGHT=0.200   MHZ      ', F20.6, 2X, I8)
 1159 FORMAT('   RMS.DEV.WEIGHT=1.000   MHZ      ', F20.6, 2X, I8)

 1161 FORMAT('     RMS.DEV.WEIGHT=0.004   MHZ A  ', F20.6, 2X, I8)
 1162 FORMAT('     RMS.DEV.WEIGHT=0.010   MHZ A  ', F20.6, 2X, I8)
 1163 FORMAT('     RMS.DEV.WEIGHT=0.020   MHZ A  ', F20.6, 2X, I8)
 1164 FORMAT('     RMS.DEV.WEIGHT=0.045   MHZ A  ', F20.6, 2X, I8)
 1165 FORMAT('     RMS.DEV.WEIGHT=0.070   MHZ A  ', F20.6, 2X, I8)
 1166 FORMAT('     RMS.DEV.WEIGHT=0.080   MHZ A  ', F20.6, 2X, I8)
 1167 FORMAT('     RMS.DEV.WEIGHT=0.100   MHZ A  ', F20.6, 2X, I8)
 1168 FORMAT('     RMS.DEV.WEIGHT=0.200   MHZ A  ', F20.6, 2X, I8)
 1169 FORMAT('     RMS.DEV.WEIGHT=1.000   MHZ A  ', F20.6, 2X, I8)

 1171 FORMAT('     RMS.DEV.WEIGHT=0.004   MHZ E  ', F20.6, 2X, I8)
 1172 FORMAT('     RMS.DEV.WEIGHT=0.010   MHZ E  ', F20.6, 2X, I8)
 1173 FORMAT('     RMS.DEV.WEIGHT=0.020   MHZ E  ', F20.6, 2X, I8)
 1174 FORMAT('     RMS.DEV.WEIGHT=0.045   MHZ E  ', F20.6, 2X, I8)
 1175 FORMAT('     RMS.DEV.WEIGHT=0.070   MHZ E  ', F20.6, 2X, I8)
 1176 FORMAT('     RMS.DEV.WEIGHT=0.080   MHZ E  ', F20.6, 2X, I8)
 1177 FORMAT('     RMS.DEV.WEIGHT=0.100   MHZ E  ', F20.6, 2X, I8)
 1178 FORMAT('     RMS.DEV.WEIGHT=0.200   MHZ E  ', F20.6, 2X, I8)
 1179 FORMAT('     RMS.DEV.WEIGHT=1.000   MHZ E  ', F20.6, 2X, I8)

 800           FORMAT(' <(23)*>l = ',F7.4,' <(23)*>s = ',F7.4,
     1                ' <C2(c)>l = ',F7.4,' <C2(c)>s = ',F7.4)

 805          FORMAT(1X,I3,1X,I3,1X,I3,1X,I3,1X,A1,1X,I3,1X,I3,1X,I3,1X,
     1    I3,1X,A1,1X,F12.3,6X,F12.3,2X,F7.3,1X,E12.6,1X,A6,1X,A1,1X,A1,
     2             ' El = ',F12.6,' Eu = ',F12.6,' BLTZ = ',E12.6,
C    2             ' El = ',E12.6,' Eu = ',E12.6,' BLTZ = ',E12.6,
     3             ' SSQ = ',E12.6,' % = ', 6(I4,1H/),1H/,' CF = ',f4.2)

 808          FORMAT(1X,I3,1X,I3,1X,I3,1X,I3,1X,A1,1X,I3,1X,I3,1X,I3,1X,
     1 I3,1X,A1,11X,F12.4,1X,F16.4,1X,F16.4,1X,E12.6,1X,A6,1X,A1,1X,A1,
     2             ' El = ',F12.6,' Eu = ',F12.6,' BLTZ = ',E12.6,
C    2             ' El = ',E12.6,' Eu = ',E12.6,' BLTZ = ',E12.6,
     3             ' SSQ = ',E12.6,' % = ', 6(I4,1H/),1H/,' CF = ',f4.2)

 835             FORMAT(1X,I3,1X,I3,1X,I3,1X,I3,1X,A1,1X,I3,1X,I3,1X,I3,
     1              1X,I3,1X,A1,1X,F12.4,'(',I1,')',1X,F12.4,1X,F7.4,
     2              1X,E12.6,1X,A6,1X,A1,         
     3              ' El = ',F12.6,' Eu = ',F12.6,' BLTZ = ',E12.6,
C    3              ' El = ',E12.6,' Eu = ',E12.6,' BLTZ = ',E12.6,
     4           ' SSQ = ',E12.6,' % = ', 6(I4,1H/),1H/,' CF = ',f4.2)

 834           FORMAT(1X,I3,1X,I3,1X,I3,1X,I3,1X,A1,1X,I3,1X,I3,1X,I3,
     1              1X,I3,1X,A1,1X,F12.3,'(',I3,')',1X,F12.3,2X,F7.3,
     2              1X,E12.6,1X,A6,1X,A1,          
     3              ' El = ',F12.6,' Eu = ',F12.6,' BLTZ = ',E12.6,
C     3              ' El = ',E12.6,' Eu = ',E12.6,' BLTZ = ',E12.6,
     4              ' SSQ = ',E12.6,' % = ', 6(I4,1H/),1H/,' CF = ',
     5              f4.2)

 602  FORMAT(/,4X,' UPPER STATE',3X,' LOWER STATE',5X,                          
     1     ' OBSERVATION ',8X,' CALCULATION',6X,' TME ',4X,' WEIGHT',           
     2     4X,' FITTED OR NO'                                                   
     3     /,2(2X,' VT J KA KC PARITY'))                                        

 604              FORMAT(
     1                 1X,I1,1X,I3,1X,I3,1X,I3,1X,A1,1X,I1,1X,I3,1X,I3,          
     2                 1X,I3,1X,A1,1X,F12.3,'(',I3,')',1X,F12.3,2X,F7.3,           
     3                 1X,A1)

 605        FORMAT(1X,I1,1X,I3,1X,I3,1X,I3,1X,A1,1X,I1,1X,I3,1X,I3,1X,        
     1           I3,1X,A1,1X,F12.3,6X,F12.3,2X,F7.3,1X,A1,1X,A1)             

 608        FORMAT(1X,I1,1X,I3,1X,I3,1X,I3,1X,A1,1X,I1,1X,I3,1X,I3,1X,        
     1           I3,1X,A1,1X,F12.5,1X,F15.5,1X,F12.5,1X,A1,1X,A1)           

 634           FORMAT(1X,I1,1X,I3,1X,I3,1X,I3,1X,A1,1X,I1,1X,I3,1X,I3,          

     1              1X,I3,1X,A1,1X,F12.5,'(',I1,')',1X,F16.5,1X,F17.5,
     2               1X, I5, 1X,A1)             

 635           FORMAT(1X,I1,1X,I3,1X,I3,1X,I3,1X,A1,1X,I1,1X,I3,1X,I3,        
     1              1X,I3,1X,A1,1X,F12.5,'(',I1,')',1X,F12.5,1X,F8.5,           
     2              1X,A1)

 704              FORMAT(
     1                 1X,I1,1X,I3,1X,I3,1X,I3,1X,A1,1X,I1,1X,I3,1X,I3,          
     2                 1X,I3,1X,A1,1X,F12.3,'(',I3,')',1X,F12.3,2X,F7.3,           
     3                 1X,40A)        

 705        FORMAT(1X,I1,1X,I3,1X,I3,1X,I3,1X,A1,1X,I1,1X,I3,1X,I3,1X,        
     1           I3,1X,A1,1X,F12.3,6X,F12.3,2X,F7.3,1X,A1,1X,40A)             

 708        FORMAT(1X,I1,1X,I3,1X,I3,1X,I3,1X,A1,1X,I1,1X,I3,1X,I3,1X,        
     1           I3,1X,A1,1X,F12.5,1X,F15.5,1X,F12.5,1X,A1,1X,40A)           

 734           FORMAT(1X,I1,1X,I3,1X,I3,1X,I3,1X,A1,1X,I1,1X,I3,1X,I3,          
     1              1X,I3,1X,A1,1X,F12.5,'(',I1,')',1X,F16.5,1X,F17.5,
     2              1X, I5,1X,40A)                                                            

 735           FORMAT(1X,I1,1X,I3,1X,I3,1X,I3,1X,A1,1X,I1,1X,I3,1X,I3,        
     1              1X,I3,1X,A1,1X,F12.5,'(',I1,')',1X,F12.5,1X,F8.5,           
     2              1X,40A)                                                      

 504              FORMAT(
     1                 1X,I1,1X,I3,1X,I3,1X,I3,1X,A1,1X,I1,1X,I3,1X,I3,          
     2                 1X,I3,1X,A1,1X,F12.3,'(',I3,')',1X,F12.3,2X,F7.3,           
     3                 1X,F10.7,2X,2A)        

 505        FORMAT(1X,I1,1X,I3,1X,I3,1X,I3,1X,A1,1X,I1,1X,I3,1X,I3,1X,        
     1           I3,1X,A1,1X,F12.3,6X,F12.3,2X,F7.3,1X,A1,1X,F10.7,
     2                 2X,2A)             

 508        FORMAT(1X,I1,1X,I3,1X,I3,1X,I3,1X,A1,1X,I1,1X,I3,1X,I3,1X,        
     1           I3,1X,A1,1X,F12.5,1X,F15.5,1X,F12.5,1X,A1,1X,F10.7,
     2           2X,2A)           

 534           FORMAT(1X,I1,1X,I3,1X,I3,1X,I3,1X,A1,1X,I1,1X,I3,1X,I3,          
     1              1X,I3,1X,A1,1X,F12.5,'(',I1,')',1X,F16.5,1X,F17.5,
     2              1X, I5,1X,F10.7,2X,2A)                                                            

 535           FORMAT(1X,I1,1X,I3,1X,I3,1X,I3,1X,A1,1X,I1,1X,I3,1X,I3,        
     1              1X,I3,1X,A1,1X,F12.5,'(',I1,')',1X,F12.5,1X,F8.5,           
     2              1X,F10.7,2X,2A)                                                      

      RETURN                                                                    
      END                                                                       

