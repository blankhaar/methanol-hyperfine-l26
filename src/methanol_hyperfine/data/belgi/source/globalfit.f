C     
C     Updated by Jonathan Fisher (q41aq@unb.ca)
C     
C     The program is now dependent on BLAS and LAPACK. The
C     appropriate libraries must be linked when compiling.
C     
C     The program has now been split into multiple source files;
C     makefiles for GNU Fortran and Portland Group F77 are 
C     included. Common block definitions have been placed into
C     their own sources files; procedures can simply include the
C     appropriate files (i.e. common_xxxxx.f).
C     
C     I have consistently added implicit none to every procedure;
C     although this is not standard Fortran 77, it is a common
C     extension and is a very good practice, as it forces variables
C     to be explicitly typed.
C     
C     I am currently in the process of going through every procedure
C     and cleaning up the code, adding comments, etc. In its current
C     state it is very difficult to follow.
C     

C*******************************************************************
C     I.Kleiner and M.Godefroid (Free University of Brussels,                   
C     50, Av. F-D. Roosevelt, 1050 Brussels, Belgium, 1987-91)                  
*********************************************************************


C     20 aout 1992                                                              
C     
C     This program is a least-squares procedure able to fit                     
C     data from a molecule with torsion moving. It is based                     
C     on the least-squares program written by M. Godefroid                       
C     "rtfit" and on the rotator-torsion program "rottor".                       
C     The V matrix elements are calculated as the derivatives                   
C     of the energy to the parameter to be fitted.                               
C     This version can trait microawve or infrared data.                          
C     This version is the copy of torfit but was modified in                    
C     order to run in the cdc nbs computer and is a weighted 
C     least-squares.
C     This program includes v9 term,dabj*j*(j+1),dabk*k**2,rhoj and             
C     rhok.                                                                      
C     fitv9b can trait high j values thanks a special procedure of              
C     picking up eigenvalues that involves a sum of the square of               
C     coeffiecients c(k) et c(-k) (see function iposa and ipose)                
C     in this version only the weights of the mw has changed                    
C     (each mw line can receive a different weight).                            
C     in this version the output has changed (see linecal)                      
C     new terms have been added:dc1j,c2j,ak2j,gvj,amv,anv,bk2                   
C     
C     dac,dbc,c4,odelta,ak3j,ak3k                                               
C     rho has been replaced by rhorho and is only fitted for htors.             
C     fitwgt9 is a version where the diagonal matrix elements in k              
C     and off-diag. in vt have been added,contrarly to herbst.                  
C     new parameter in fitwg10:d delta                                          
C     ftjon10:d jon has modified the program...                                 
C     and i have add ak1j,ak1k,ak2k,ak5k                                        
C     fastfit:version that has been modified in order to gain time              
C     (see subroutines setup and vset).                                         
C     the energies and eigenvectors are calculated at the very beginni          
C     of each cycle for j=0 to 14 for every k and vt states and then            
C     those energies and eigenvectors are stored in 2 big matrices              
C     (deriv and ener)                                                          
C     this is the last (improved) version of girou (circe):iksrc18              
C     for vt=2+ipose,iposa modified (iksrc20)                                   
C     belgi28: changes by ortigoso, aout 1992                                   
C     belgi30 : changes june 1993; the iposa has now a parity calc
C     form (23)* operator and the basis set has been extended from 7 to 9
C     can treat up to vt=4
C     
C     march, 1994       li-hong xu
C     (i) the iposa and ipose have been rewritten for methanol.
C     (ii) parameters programed in have been increased to 80.
c     
c************************************************************************

      program globalfit
      implicit none
      save

      include 'common_dim.f'
      include 'common_eigen.f'
      include 'common_eigen1.f'
      include 'common_tor.f'
      include 'common_dat.f'
      include 'common_rotor.f'
      include 'common_par.f'
      include 'common_big2.f'
      include 'common_expdat.f'

      include 'common_cntl.f'
      include 'common_symb.f'
      include 'common_quant.f'
      include 'common_chara.f'
      include 'common_cst.f'
      include 'common_int.f'

      double precision vtldv, vtldtm, deltax, correl
      double precision v, v1, v2, unit
      double precision delx, delxst
      integer iblsto, ipiv, npp, vt, sigma
      double precision sst, sir10, sir21, sig010, r
      integer lo, ipr, ivv, ink, ii, i, j, k, jj, jjj
      double precision elower, eupper, acc
      double precision xx, rjth, viben
      

      character*64 paramfmt
      character*((maxparam+1)*6) fitparam
      dimension vtldv(maxparam,maxparam),vtldtm (maxparam)
      dimension deltax(maxparam), unit(maxparam,maxparam),v(maxparam)
      dimension v1(maxparam),v2(maxparam), delx(maxparam)
      dimension delxst(maxparam),iblsto(2),ipiv(maxparam)           
      integer ifail, ia, kdat,ipath,irank
      double precision worki(lwork), varcov, prod, prod1
      dimension correl(maxparam, maxparam)
      double precision v1u(maxparam, ndta), v2u(maxparam, ndta)
      double precision vu(maxparam, ndta), vnou(maxparam, ndta)
      double precision v1nou(maxparam, ndta), v2nou(maxparam, ndta)
      double precision svd,vr,ul,tol

	  

C     JMF - these are used to calculate kappa
C     double precision kappa, wr(maxparam), wi(maxparam)
C     double precision lambda(2), lmb
C     double precision det
C     JMF

C     
C     this program may trait (ndta) observed data and fit (maxparam)                
C     parameters (upper and/or lower).                                          
C     
C     call the subroutine which read the data'                                  
C     
      write(6,*) ' BEGIN PROGRAM   '                                            
      call data1(cst)                                                           
C     
      term= evib(2)-evib(1)                                                     
      termst=term                                                               
C     
C     beginning of the fitting procedure                                        
C     
      np12=np(1)+np(2)                                                          
      npp=np12                                                                  
C     print*, 'npp= ',npp,' np12= ',np12                                        
      if(ito.ne.0) npp=np12+1                                                   
C     
C     initialisations                                                           
C     
      ipr = 0
      iend = -1
      sst = 0.d0

      do while ( iend .ne. 0 )

         if  ( ibug5 .eq. 5 ) call rule(cst)
         
         if ( iend .lt. 0 ) then
 23         nit=1                                                                     
            sst=0.d0                                                                  
            do i=1,npp                                                             
               delxst(i)=0.d0                                                         
            end do
         end if
         
 25      do i=1,npp                                                              
            v(i)=0.d0
            
            if ( unc .eq. 1 ) then
               do jjj = 1, ndata
                  vu(i, jjj) = 0.d0
                  vnou(i, jjj) = 0.d0
               end do
            end if
            
            vtldtm (i)=0.d0                                                        
            do j=1,npp                                                           
               vtldv(i,j)=0.d0                                                     
            end do
         end do

         if ( unc .eq. 1 ) then 
            do ia = 1, ndata
               aincer(ia) = 0.0
            end do
         end if 
         
         call factt

         iblsto(1)=0.d0                                                            
         iblsto(2)=0.d0                                                            
         s=0.d0                                                                    
         smw=0.d0                                                                 
         sir=0.d0                                                                  
         sir10=0.d0                                                                
         sir21=0.d0                                                                
         smw0=0.d0                                                                 
         smw1=0.d0                                                                 
         smw2=0.d0                                                                 
         smw3=0.d0
         sig004=0.d0                                                               
         sig080=0.d0                                                               
         sig010=0.d0                                                               
         sig100=0.d0                                                               
         sig020=0.d0                                                               
         sig200=0.d0                                                               
         sig045=0.d0                                                               
         sig1m0=0.d0                                                               
         sgst0=0.d0                                                                
         sw070=0.d0                                                                
         nofit=0                                                                   
         nfitmw=0                                                                  
         nfitir=0                                                                  
         nfir10=0                                                                  
         nfir21=0                                                                  
         nfimw0 =0                                                                 
         nfimw1 =0                                                                 
         nfimw2=0                                                                  
         nfi004 =0                                                                 
         nfi080 =0                                                                 
         nfi010 =0                                                                 
         nfi100 =0                                                                 
         nfi020 =0                                                                 
         nfi200 =0                                                                 
         nfi045 =0                                                                 
         nfi1m0 =0                                                                 
         ngst0=0                                                                   
         nfi070 =0                                                                 
         nfimw3=0

         NFTW1 = 0
         NFTW2 = 0
         NFTW3 = 0
         NFTW4 = 0

         NFTW1A = 0
         NFTW2A = 0
         NFTW3A = 0
         NFTW4A = 0

         NFTW1E = 0
         NFTW2E = 0
         NFTW3E = 0
         NFTW4E = 0

         nft1 = 0
         nft2 = 0
         nft1a = 0
         nft1e = 0
         nft2a = 0
         nft2e = 0

         NFT00 = 0
         NFT01 = 0
         NFT02 = 0
         NFT11 = 0
         NFT12 = 0
         NFT22 = 0

         NFT00A = 0
         NFT01A = 0
         NFT02A = 0
         NFT11A = 0
         NFT12A = 0
         NFT22A = 0

         NFT00E = 0
         NFT01E = 0
         NFT02E = 0
         NFT11E = 0
         NFT12E = 0
         NFT22E = 0


         NMW00 = 0
         NMW01 = 0
         NMW02 = 0
         NMW11 = 0
         NMW12 = 0
         NMW22 = 0

         NMW00A = 0
         NMW01A = 0
         NMW02A = 0
         NMW11A = 0
         NMW12A = 0
         NMW22A = 0

         NMW00E = 0
         NMW01E = 0
         NMW02E = 0
         NMW11E = 0
         NMW12E = 0
         NMW22E = 0

         nfi004a =0                                                                 
         nfi080a =0                                                                 
         nfi010a =0                                                                 
         nfi100a =0                                                                 
         nfi020a =0                                                                 
         nfi200a =0                                                                 
         nfi045a =0                                                                 
         nfi1m0a =0  
         
         nfi004e =0                                                                 
         nfi080e =0                                                                 
         nfi010e =0                                                                 
         nfi100e =0                                                                 
         nfi020e =0                                                                 
         nfi200e =0                                                                 
         nfi045e =0                                                                 
         nfi1m0e =0                                                                 

         if(np12.gt.(maxparam)) then                                                    
            write(6,*) 'TOO MANY PARAMETERS TO FIT'                                
            stop                                                                   
         endif                                                                     
C     
C     beginning of the do loop on data,if one doesn't want to                   
C     fit a data,iast can't be equal to 1.ndata is the total number             
C     of data including sigma=0 and 1.                                           

C     we have modified the procedure here:din the matrix deriv, onea            
C     is calculating the energies and derivatives for all the j,k levels        
C     at once and then one will pick up the "good" values.                      
         call setup(cst)                                                           
         do  k=1,ndata              


            
            if (ivobs(1,k).gt.vtmax.or.ivobs(2,k).gt.vtmax) cycle
            if (nobs(1,k).gt.jmax.or.nobs(2,k).gt.jmax) cycle                  

            if (iast(k).eq.1) then

               nofit=nofit+1              

               if ( iw(k) .lt. 1000) then
                  if ( iw(k) .eq. 1 ) NFTW1 = NFTW1 + 1
                  if ( iw(k) .eq. 2 ) NFTW2 = NFTW2 + 1
                  if ( iw(k) .eq. 3 ) NFTW3 = NFTW3 + 1
                  if ( iw(k) .eq. 4 ) NFTW4 = NFTW4 + 1

                  if ( probs(1, k) .eq. ' ') then
                     if ( iw(k) .eq. 1 ) NFTW1E = NFTW1E + 1
                     if ( iw(k) .eq. 2 ) NFTW2E = NFTW2E + 1
                     if ( iw(k) .eq. 3 ) NFTW3E = NFTW3E + 1
                     if ( iw(k) .eq. 4 ) NFTW4E = NFTW4E + 1
                  else
                     if ( iw(k) .eq. 1 ) NFTW1A = NFTW1A + 1
                     if ( iw(k) .eq. 2 ) NFTW2A = NFTW2A + 1
                     if ( iw(k) .eq. 3 ) NFTW3A = NFTW3A + 1
                     if ( iw(k) .eq. 4 ) NFTW4A = NFTW4A + 1
                  endif
               endif

               if ( iw(k) .lt. 1000 ) then
                  if ( (ivobs(1, k).eq.0) .and. (ivobs(2,k).eq.0)) then
                     NFT00 = NFT00 + 1
                  else if ( ( (ivobs(1, k).eq.0) 
     &                    .and. (ivobs(2,k).eq.1) )
     &                    .or. ( (ivobs(1,k).eq.1) 
     &                    .and. (ivobs(2,k).eq.0)) ) then
                     NFT01 = NFT01 + 1
                  else if ( ( (ivobs(1, k).eq.0) 
     &                    .and. (ivobs(2,k).eq.2) )
     &                    .or. ( (ivobs(1,k).eq.2) 
     &                    .and. (ivobs(2,k).eq.0)) ) then
                     NFT02 = NFT02 + 1
                  else if ((ivobs(1, k).eq.1).and.(ivobs(2,k).eq.1)) 
     &                    then
                     NFT11 = NFT11 + 1
                  else if ( ( (ivobs(1, k).eq.1) 
     &                    .and. (ivobs(2,k).eq.2) )
     &                    .or. ( (ivobs(1,k).eq.2) 
     &                    .and. (ivobs(2,k).eq.1)) ) then
                     NFT12 = NFT12 + 1
                  else if ((ivobs(1, k).eq.2).and.(ivobs(2,k).eq.2)) 
     &                    then
                     NFT22 = NFT22 + 1
                  endif
               endif

               if ( probs(1, k) .eq. ' ' ) then
                  if ( iw(k) .lt. 1000 ) then
                     if ( (ivobs(1, k).eq.0) .and. (ivobs(2,k).eq.0))
     &                    then
                        NFT00E = NFT00E + 1
                     else if ( ( (ivobs(1, k).eq.0) 
     &                       .and. (ivobs(2,k).eq.1) )
     &                       .or. ( (ivobs(1,k).eq.1) 
     &                       .and. (ivobs(2,k).eq.0)) ) then
                        NFT01E = NFT01E + 1
                     else if ( ( (ivobs(1, k).eq.0) 
     &                       .and. (ivobs(2,k).eq.2) )
     &                       .or. ( (ivobs(1,k).eq.2) 
     &                       .and. (ivobs(2,k).eq.0)) ) then
                        NFT02E = NFT02E + 1
                     else if ((ivobs(1, k).eq.1).and.(ivobs(2,k).eq.1))
     &                       then
                        NFT11E = NFT11E + 1
                     else if ( ( (ivobs(1, k).eq.1) 
     &                       .and. (ivobs(2,k).eq.2) )
     &                       .or. ( (ivobs(1,k).eq.2) 
     &                       .and. (ivobs(2,k).eq.1)) ) then
                        NFT12E = NFT12E + 1
                     else if ((ivobs(1, k).eq.2).and.(ivobs(2,k).eq.2))
     &                       then
                        NFT22E = NFT22E + 1
                     endif
                  endif

               else

                  if ( iw(k) .lt. 1000 ) then
                     if ( (ivobs(1, k).eq.0) .and. (ivobs(2,k).eq.0)) 
     &                    then
                        NFT00A = NFT00A + 1
                     else if ( ( (ivobs(1, k).eq.0) 
     &                       .and. (ivobs(2,k).eq.1) )
     &                       .or. ( (ivobs(1,k).eq.1) 
     &                       .and. (ivobs(2,k).eq.0)) ) then
                        NFT01A = NFT01A + 1
                     else if ( ( (ivobs(1, k).eq.0) 
     &                       .and. (ivobs(2,k).eq.2) )
     &                       .or. ( (ivobs(1,k).eq.2) 
     &                       .and. (ivobs(2,k).eq.0)) ) then
                        NFT02A = NFT02A + 1
                     else if ((ivobs(1, k).eq.1).and.(ivobs(2,k).eq.1))
     &                       then
                        NFT11A = NFT11A + 1
                     else if ( ( (ivobs(1, k).eq.1) 
     &                       .and. (ivobs(2,k).eq.2) )
     &                       .or. ( (ivobs(1,k).eq.2) 
     &                       .and. (ivobs(2,k).eq.1)) ) then
                        NFT12A = NFT12A + 1
                     else if ((ivobs(1, k).eq.2).and.(ivobs(2,k).eq.2))
     &                       then
                        NFT22A = NFT22A + 1
                     endif
                  endif

               endif

               if ( iw(k) .ge. 1000 ) then
                  if ( (ivobs(1, k).eq.0) .and. (ivobs(2,k).eq.0)) then
                     NMW00 = NMW00 + 1
                  else if ( ( (ivobs(1, k).eq.0) 
     &                    .and. (ivobs(2,k).eq.1) )
     &                    .or. ( (ivobs(1,k).eq.1) 
     &                    .and. (ivobs(2,k).eq.0)) ) then
                     NMW01 = NMW01 + 1
                  else if ( ( (ivobs(1, k).eq.0) 
     &                    .and. (ivobs(2,k).eq.2) )
     &                    .or. ( (ivobs(1,k).eq.2) 
     &                    .and. (ivobs(2,k).eq.0)) ) then
                     NMW02 = NMW02 + 1
                  else if ((ivobs(1, k).eq.1).and.(ivobs(2,k).eq.1)) 
     &                    then
                     NMW11 = NMW11 + 1
                  else if ( ( (ivobs(1, k).eq.1) 
     &                    .and. (ivobs(2,k).eq.2) )
     &                    .or. ( (ivobs(1,k).eq.2) 
     &                    .and. (ivobs(2,k).eq.1)) ) then
                     NMW12 = NMW12 + 1
                  else if ((ivobs(1, k).eq.2).and.(ivobs(2,k).eq.2)) 
     &                    then
                     NMW22 = NMW22 + 1
                  endif
               endif

               if ( probs(1, k) .eq. ' ' ) then
                  if ( iw(k) .ge. 1000 ) then
                     if ( (ivobs(1, k).eq.0) .and. (ivobs(2,k).eq.0))
     &                    then
                        NMW00E = NMW00E + 1
                     else if ( ( (ivobs(1, k).eq.0) 
     &                       .and. (ivobs(2,k).eq.1) )
     &                       .or. ( (ivobs(1,k).eq.1) 
     &                       .and. (ivobs(2,k).eq.0)) ) then
                        NMW01E = NMW01E + 1
                     else if ( ( (ivobs(1, k).eq.0) 
     &                       .and. (ivobs(2,k).eq.2) )
     &                       .or. ( (ivobs(1,k).eq.2) 
     &                       .and. (ivobs(2,k).eq.0)) ) then
                        NMW02E = NMW02E + 1
                     else if ((ivobs(1, k).eq.1).and.(ivobs(2,k).eq.1))
     &                       then
                        NMW11E = NMW11E + 1
                     else if ( ( (ivobs(1, k).eq.1) 
     &                       .and. (ivobs(2,k).eq.2) )
     &                       .or. ( (ivobs(1,k).eq.2) 
     &                       .and. (ivobs(2,k).eq.1)) ) then
                        NMW12E = NMW12E + 1
                     else if ((ivobs(1, k).eq.2).and.(ivobs(2,k).eq.2))
     &                       then
                        NMW22E = NMW22E + 1
                     endif
                  endif

               else

                  if ( iw(k) .ge. 1000 ) then
                     if ( (ivobs(1, k).eq.0) .and. (ivobs(2,k).eq.0))
     & then
                        NMW00A = NMW00A + 1
                     else if ( ( (ivobs(1, k).eq.0) 
     &                       .and. (ivobs(2,k).eq.1) )
     &                       .or. ( (ivobs(1,k).eq.1) 
     &                       .and. (ivobs(2,k).eq.0)) ) then
                        NMW01A = NMW01A + 1
                     else if ( ( (ivobs(1, k).eq.0) 
     &                       .and. (ivobs(2,k).eq.2) )
     &                       .or. ( (ivobs(1,k).eq.2) 
     &                       .and. (ivobs(2,k).eq.0)) ) then
                        NMW02A = NMW02A + 1
                     else if ((ivobs(1, k).eq.1).and.(ivobs(2,k).eq.1))
     &                       then
                        NMW11A = NMW11A + 1
                     else if ( ( (ivobs(1, k).eq.1) 
     &                       .and. (ivobs(2,k).eq.2) )
     &                       .or. ( (ivobs(1,k).eq.2) 
     &                       .and. (ivobs(2,k).eq.1)) ) then
                        NMW12A = NMW12A + 1
                     else if ((ivobs(1, k).eq.2).and.(ivobs(2,k).eq.2))
     &                       then
                        NMW22A = NMW22A + 1
                     endif
                  endif

               endif

               if(iw(k).ge.1000) nfitmw=nfitmw+1                                      
               if(iw(k).lt.1000) nfitir=nfitir+1                           

               if(iw(k).eq.1) then
                  nfir10=nfir10+1                                         
                  nft1 = nft1 + 1

                  if ( probs(1, k) .eq. ' ' ) then
                     nft1e = nft1e + 1
                  else
                     nft1a = nft1a + 1
                  endif
                  
               endif

               if(iw(k).eq.2) then 
                  nfir21=nfir21+1                                         
                  nft2 = nft2 + 1

                  if ( probs(1, k) .eq. ' ' ) then
                     nft2e = nft2e + 1
                  else
                     nft2a = nft2a + 1
                  endif

               end if

               if(iw(k).ge.1000.and.ivobs(1,k).eq.0) nfimw0 =nfimw0 +1                
               if(iw(k).ge.1000.and.ivobs(1,k).eq.1) nfimw1 =nfimw1 +1                
               if(iw(k).ge.1000.and.ivobs(1,k).eq.2) nfimw2 =nfimw2 +1                
               if(iw(k).ge.1000.and.ivobs(1,k).eq.3) nfimw3 =nfimw3 +1                
C     if(iw(k).eq.2000) nfitgd=nfitgd+1                                         
               if(iw(k).eq.10) ngst0=ngst0+1                                          
               if(iw(k).eq.1004.or.iw(k).eq.1005) nfi004 =nfi004 +1                   
               if(iw(k).eq.1080) nfi080 =nfi080 +1                                    
               if(iw(k).eq.1010.or.iw(k).eq.1008)                                     
     *              nfi010 =nfi010 +1                                                 
               if(iw(k).eq.1070) nfi070 =nfi070 +1                                    
               if(iw(k).eq.1100) nfi100 =nfi100 +1                                    
               if(iw(k).eq.1020) nfi020 =nfi020 +1                                    
               if(iw(k).le.1300.and.iw(k).ge.1130) nfi200 =nfi200 +1                 
               if(iw(k).eq.1030.or.iw(k).eq.1040.or.iw(k).eq.1050)                    
     *              nfi045 =nfi045 +1                                                 
               if(iw(k).eq.1990) nfi1m0 =nfi1m0 +1                                    

               if ( probs(1, k) .eq. ' ' ) then
                  if(iw(k).eq.1004.or.iw(k).eq.1005) 
     *                 nfi004e =nfi004e +1                   
                  if(iw(k).eq.1080) nfi080e =nfi080e +1                                    
                  if(iw(k).eq.1010.or.iw(k).eq.1008)                                     
     *                 nfi010e =nfi010e +1                                                 
                  if(iw(k).eq.1070) nfi070e =nfi070e +1                                    
                  if(iw(k).eq.1100) nfi100e =nfi100e +1                                    
                  if(iw(k).eq.1020) nfi020e =nfi020e +1                                    
                  if(iw(k).le.1300.and.iw(k).ge.1130)
     *                 nfi200e =nfi200e +1                 
                  if(iw(k).eq.1030.or.iw(k).eq.1040.or.
     *                 iw(k).eq.1050) nfi045e =nfi045e +1                                                 
                  if(iw(k).eq.1990) nfi1m0e =nfi1m0e +1                                    
               else
                  if(iw(k).eq.1004.or.iw(k).eq.1005) 
     *                 nfi004a =nfi004a +1                   
                  if(iw(k).eq.1080) nfi080a =nfi080a +1                                    
                  if(iw(k).eq.1010.or.iw(k).eq.1008)                                     
     *                 nfi010a =nfi010a +1                                                 
                  if(iw(k).eq.1070) nfi070a =nfi070a +1                                    
                  if(iw(k).eq.1100) nfi100a =nfi100a +1                                    
                  if(iw(k).eq.1020) nfi020a =nfi020a +1                                    
                  if(iw(k).le.1300.and.iw(k).ge.1130)
     *                 nfi200a =nfi200a +1                 
                  if(iw(k).eq.1030.or.iw(k).eq.1040.or.
     *                 iw(k).eq.1050) nfi045a =nfi045a +1                                                 
                  if(iw(k).eq.1990) nfi1m0a =nfi1m0a +1                                    
               endif
            endif

C      if(nitt.eq.0 .and. unc .eq. 0) cycle
C     setup the matrix v,for current data                                       
            do i=1,np12                                                       
               v1(i)=0.d0                                                          
               v2(i)=0.d0                                                          

               if ( unc .eq. 1 ) then
                  v1u(i, k) = 0.d0
                  v1nou(i, k) = 0.d0
                  v2u(i, k) = 0.d0
                  v2nou(i, k) = 0.d0
               end if
            end do
C     
C     here we are calculating the correct indices in order to fetcha            
C     the energies and derivatives corresponding to the observed transit        
C     

            do lo=1,2                                                  
               n=nobs(lo,k)                                                        
               kk=kobs(lo,k)                                                       
               iv=ivobs(lo,k)                                                      
               pr=probs(lo,k)                                                      
               if(pr.eq.'+') ipr=1                                                 
               if(pr.eq.'-') ipr=0                                                 
               if(k.gt.nadat) then                                                 
                  isig=1                                                           
                  ivv=iv+6                                                         
                  ink=n**2+n+kk+1                                                  
               else                                                                
                  isig=0                                                           
                  ivv=iv+1                                                         
                  ink=n**2+n+kk+1-ipr*n                                            
               endif                                                               
               if(irmw.eq.2) then
                  if(lo.eq.1) then                                                    
                     elower=ener(ivv,ink,1)                                     
                  else                                                                
                     eupper=ener(ivv,ink,2)                                  
                  endif                                                               
               elseif(irmw.eq.1) then
                  if(lo.eq.1) then                                                    
                     elower=ener(ivv,ink,1)                                  
                  else                                                                
                     eupper=ener(ivv,ink,1)                                      
                  endif                                                               
               endif                                                               

               do i=1,np12
                  if(irmw.eq.2) then
                     if(lo.eq.1) then
                        v1(i)= derivl(ivv,ink,i)
c       if(iv.EQ.0) print*,'N=',N,'KK=',KK,'v1=',derivl(ivv,ink,i),'i=',i
c     & ,'Ener=',ener(ivv,ink,1)
C     if(n.eq.1) print*,'v1 dans main=',v1(i)
                     else
                        v2(i)=derivs(ivv,ink,i)
c      if(iv.EQ.0) print*,'N=',N,'KK=',KK,'v2=',derivs(ivv,ink,i),'i=',i
c     & ,'Ener=',ener(ivv,ink,1)
C     if(n.eq.1) print*,'v2 dans main=',v2(i)
                     endif
                  elseif(irmw.eq.1) then
                     if(lo.eq.1) then
                        v1(i)= derivl(ivv,ink,i)
c       if(iv.EQ.0) print*,'N=',N,'KK=',KK,'v1=',derivl(ivv,ink,i),'i=',i
c     & ,'Ener=',ener(ivv,ink,1)
                     else
                        v2(i)= derivl(ivv,ink,i)
c       if(iv.EQ.0) print*,'N=',N,'KK=',KK,'v2=',derivl(ivv,ink,i),'i=',i
c     & ,'Ener=',ener(ivv,ink,1)
                     endif
                  endif

                  if ( iast(k) .ne. 0 ) then
                  v(i)=v2(i)-v1(i)
       print*,'N=',N,'KK=',KK,'v=',v(i),'v2=',v2(i),'v1=',v1(i),'i=',i
     & ,'Ener Lower=',ener(ivv,ink,1),'Ener Upper=',ener(ivv,ink,2),
     &'data number=',K,'etrans=',etrans(K)*29979.2458           
                  endif
               end do

               if ( unc .eq. 1 ) then
                  do i=1,np12
                     if(irmw.eq.2) then
                        if(lo.eq.1) then
                           if ( iast(k) .ne. 0 ) then
                              v1u(i, k) = v1(i)
                           else
                              v1nou(i, k) = v1(i)
                           endif
C     if(n.eq.1) print*,'v1 dans main=',v1(i)
                        else
                           if ( iast(k) .ne. 0 ) then
                              v2u(i, k) = v2(i)
                           else
                              v2nou(i, k) = v2(i)
                           endif                           
C     if(n.eq.1) print*,'v2 dans main=',v2(i)
                        endif
                     elseif(irmw.eq.1) then
                        if(lo.eq.1) then
                           if ( iast(k) .ne. 0 ) then
                              v1u(i, k) = v1(i)
                           else
                              v1nou(i, k) = v1(i)
                           endif
                        else
                           if ( iast(k) .ne. 0 ) then
                              v2u(i, k) = v2(i)
                           else
                              v2nou(i, k) = v2(i)
                           endif
                        endif
                     endif
                     if ( iast(k) .ne. 0 ) then
                        vu(i, k) = v2u(i, k) - v1u(i, k)
                     else
                        vnou(i, k) = v2nou(i, k) - v1nou(i, k)
                     endif
                  end do
               end if
            end do
C     
C     if the vib.energy has to be fit,the diag element of the                   
C     v matrix will contain 1.                                                  
C     here,it seems to be an error on the previous way of calculation!          
C     
            if(ito.ne.0) then                                                      
               v(npp)=1.d0                                                         

               if ( iast(k) .ne. 0 ) then
                  vu(npp, k) = 1.d0
               else
                  vnou(npp, k) = 1.d0
               endif
C     v(npp)=0.d0                                                               
C     do 8 i=1,nrotor                                                           
C     v(npp)=v(npp)+egvc(i,nbr(2))**2                                           
C     print*,' v(npp)= ',v(npp)                                                 
            endif                                                                  
C     
C     set up the vtilde v matrix(np*np)                                         
C     
            if ( iast(k) .eq. 1 ) then
               do i=1,npp                                                          
                  do j=1,npp                                                      
                     vtldv(i,j)=vtldv(i,j)+v(i)*v(j)*w(k)                             
                     correl(i, j) = sqrt( vtldv(i, i) * vtldv(j, j) );
                     
                     if ( correl(i, j) .ne. 0.0 ) then
                        correl(i, j) = vtldv(i, j) / correl(i, j)
                     end if
                     
                  end do
               end do


C     
C     observed minus calculated                                                 
C     
               tme(k)=etrans(k)-term-eupper   +elower                                 
C     standard deviation for each line,weighted                                 
               s=s+tme(k)**2*w(k)                                                     
               if(iw(k).ge.1000) then                                                 
                  smw=smw+tme(k)**2*w(k)                                              
                  if(ivobs(1,k).eq.0.and.iw(k).ge.1000) then                          
                     smw0=smw0+tme(k)**2*w(k)                                         
                  elseif(ivobs(1,k).eq.1.and.iw(k).ge.1000) then                      
                     smw1=smw1+tme(k)**2*w(k)                                         
                  elseif(ivobs(1,k).eq.2.and.iw(k).ge.1000) then                            
                     smw2=smw2+tme(k)**2*w(k)                                                 
                  elseif(ivobs(1,k).eq.3.and.iw(k).ge.1000) then                            
                     smw3=smw3+tme(k)**2*w(k)                                                 
                  endif                                                               
C     elseif (iw(k).eq.2000) then                                               
C     swgood=swgood+tme(k)**2*w(k)                                              
                  if(iw(k).eq.1004.or.iw(k).eq.1005) then                             
                     sig004=sig004+tme(k)**2*w(k)                                     
                  elseif(iw(k).eq.1080) then                                          
                     sig080=sig080+tme(k)**2*w(k)                                     
                  elseif(iw(k).eq.1010.or.iw(k).eq.1008) then                         
                     sig010=sig010+tme(k)**2*w(k)                                     
                  elseif(iw(k).eq.1100) then                                          
                     sig100=sig100+tme(k)**2*w(k)                                     
                  elseif(iw(k).eq.1020) then                                          
                     sig020=sig020+tme(k)**2*w(k)                                     
                  elseif(iw(k).eq.1070) then                                          
                     sw070=sw070+tme(k)**2*w(k)                                       
                  elseif(iw(k).le.1300.and.iw(k).ge.1130) then             
                     sig200=sig200+tme(k)**2*w(k)                                     
                  elseif(iw(k).eq.1030.or.iw(k).eq.1040.or.iw(k).eq.                  
     1                    1050) then                                               
                     sig045=sig045+tme(k)**2*w(k)                                     
                  elseif(iw(k).eq.1990) then                                          
                     sig1m0=sig1m0+tme(k)**2*w(k)                                     
                  endif                                                               
               elseif(iw(k).eq.1) then                                                
                  sir10=sir10+tme(k)**2*w(k)                                          
C     sir=sir+tme(k)**2*w(k)                                              
               elseif(iw(k).eq.2) then                                                
                  sir21=sir21+tme(k)**2*w(k)                                          
C     sir=sir+tme(k)**2*w(k)                                              
               elseif(iw(k).eq.10) then                                               
                  sgst0=sgst0+tme(k)**2*w(k)                                          
               endif                                                                  

               if ( iw(k) .lt. 1000) sir = sir + tme(k)**2 * w(k)

C     if(iw(k).ge.1000) tme(k)=tme(k)*29979.2458d0                             
               if(ibgtme.eq.1) write(6,300) k,tme(k),eupper,elower,s                  
 300           format(' K= ',i5,2x,' TME= ',f17.9,2x,' EUPPER= ',f16.8,            
     1              2x,' ELOWER= ',f16.8,' S= ',f20.8)                             
               do i=1,npp                                                          
                  vtldtm (i)=vtldtm (i)+v(i)*tme(k)*w(k)                              
               end do

            endif
C     
C     end of the loop on data                                                   
C     
         end do
C     sec=tcpu()                                                                
C     write(6,*) ' end   do 4   ',seC                                           
C     write(6,*) ' END   do 4   '                                               
C         if(nitt.eq.0) exit

         if ( nitt .ne. 0 ) then
            print*,'S(BEFORE FITTING)= ',s                                            
         endif
         s=sqrt(s/(nofit-npp))                                                     
         if(nfitmw.eq.0) then                                                      
            smw=0.d0                                                               
         else                                                                      
            smw=sqrt(smw/nfitmw)                                                   
         endif                                                                     
         if    (nfitir.eq.0) then                                                  
            sir=0.d0                                                               
         else                                                                      
            sir=sqrt(sir/nfitir)                                                   
         endif                                                                     
         if    (nfir10.eq.0) then                                                  
            sir10=0.d0                                                             
         else                                                                      
            sir10=sqrt(sir10/nfir10)                                               
         endif                                                                     
         if    (nfir21.eq.0) then                                                  
            sir21=0.d0                                                             
         else                                                                      
            sir21=sqrt(sir21/nfir21)                                               
         endif                                                                     
         if(nfimw0 .eq.0) then                                                     
            smw0=0.d0                                                              
         else                                                                      
            smw0=sqrt(smw0/nfimw0 )                                                
         endif                                                                     
         if(nfimw1 .eq.0) then                                                     
            smw1=0.d0                                                              
         else                                                                      
            smw1=sqrt(smw1/nfimw1 )                                                
         endif                                                                     
         if(nfimw2.eq.0) then                                                      
            smw2=0.d0                                                              
         else                                                                      
            smw2=sqrt(smw2/nfimw2)                                                 
         endif                                                                     
         if(nfimw3.eq.0) then                                                      
            smw3=0.d0                                                              
         else                                                                      
            smw3=sqrt(smw3/nfimw3)                                                 
         endif
         if(nfi004 .eq.0) then                                                     
            sig004=0.d0                                                            
         else                                                                      
            sig004=sqrt(sig004/nfi004 )                                            
         endif                                                                     
         if(nfi070 .eq.0) then                                                     
            sw070=0.d0                                                             
         else                                                                      
            sw070=sqrt(sw070/nfi070 )                                              
         endif                                                                     
         if(nfi080 .eq.0) then                                                     
            sig080=0.d0                                                            
         else                                                                      
            sig080=sqrt(sig080/nfi080 )                                            
         endif                                                                     
         if(nfi010 .eq.0) then                                                     
            sig010=0.d0                                                            
         else                                                                      
            sig010=sqrt(sig010/nfi010 )                                            
         endif                                                                     
         if(nfi100 .eq.0) then                                                     
            sig100=0.d0                                                            
         else                                                                      
            sig100=sqrt(sig100/nfi100 )                                            
         endif                                                                     
         if(nfi020 .eq.0) then                                                     
            sig020=0.d0                                                            
         else                                                                      
            sig020=sqrt(sig020/nfi020 )                                            
         endif                                                                     
         if(nfi200 .eq.0) then                                                     
            sig200=0.d0                                                            
         else                                                                      
            sig200=sqrt(sig200/nfi200 )                                            
         endif                                                                     
         if(nfi045 .eq.0) then                                                     
            sig045=0.d0                                                            
         else                                                                      
            sig045=sqrt(sig045/nfi045 )                                            
         endif                                                                     
         if(nfi1m0 .eq.0) then                                                     
            sig1m0=0.d0                                                            
         else                                                                      
            sig1m0=sqrt(sig1m0/nfi1m0 )                                            
         endif                                                                     
         if(ngst0.eq.0) then                                                       
            sgst0=0.d0                                                             
         else                                                                      
            sgst0=sqrt(sgst0/ngst0)                                                
         endif                                                                     


            
 123        format(' S = ', F15.10, 
     &           ' NOFIT = ', I6, ' NPP =  ', I4, 
     &           ' SMW =  ', F15.10, ' SIR =  ', F15.10)
 124        format(A10, I6)
 125        format(A10, '  ', F15.10)
            write(6, 123) S, NOFIT, NPP, SMW, SIR
            write(6, 124) ' NFITMW= ',NFITMW
            write(6, 124) ' NFITIR= ',NFITIR
            write(6, 124) ' NFIMW1= ',NFIMW1 
            write(6, 125) ' SMW0= ',SMW0
            write(6, 125) ' SMW1= ',SMW1
            write(6, 124) ' NFITMW0= ',NFIMW0
            write(6, 124) ' NFITMW1= ',NFIMW1
            write(6, 125) ' SIG004= ',SIG004
            write(6, 124) ' NFIT004= ', NFI004
            write(6, 125) ' SIG080= ',SIG080
            write(6, 124) ' NFIT080= ',NFI080 
            write(6, 125) ' SIG010= ',SIG010
            write(6, 124) ' NFIT010= ',NFI010 
            write(6, 125) ' SIG100= ', SIG100
            write(6, 124) ' NFIT100= ',NFI100
            write(6, 125) ' SIG020= ',SIG020
            write(6, 124) ' NFIT020= ', NFI020
            write(6, 125) ' SIG200= ',SIG200
            write(6, 124) ' NFIT200= ',NFI200
            write(6, 125) ' SIG045= ',SIG045
            write(6, 124) ' NFIT045= ',NFI045
            write(6, 125) ' SGST0= ',SGST0
            write(6, 124) ' NGST0= ',NGST0
            write(6, 125) ' SW070= ',SW070
            write(6, 124) ' NFIT070= ',NFI070
            write(6, 125) ' SMW2= ',SMW2
            write(6, 124) ' NFIMW2= ',NFIMW2
            write(6, 125) ' SIG1M0= ',SIG1M0              
            write(6, 124) ' NFIT1M0= ',NFI1M0
            write(6, 125) ' SIG10= ',SIR10
            write(6, 124) ' NFIT10= ',NFIR10
            write(6, 125) ' SIG21= ',SIR21
            write(6, 124) ' NFIT21= ',NFIR21
            write(6, 125) ' SMW3= ',SMW3            
            write(6, 124) ' NFIMW3= ',NFIMW3

            WRITE(6,901) NOFIT                                                        
 901        FORMAT(/,' NUMBER OF DATA INCLUDED IN THE FIT= ',I6)                      
            IF(IBGVTD .NE.0) THEN                                                     
               WRITE(6,51)                                                            
               CALL MATOUT(VTLDV,maxparam,maxparam,NPP,NPP)                                     
            ENDIF                                                                     
 51         FORMAT(/,' VTLDV MATRIX BEFORE INVERSION',/)                              


C     
C     SET UP THE VTLDV MATRIX WITH ONE ON THE DIAGONAL.                         
C     DELX   IS THE SQUARE ROOT OF THE DIAGONAL                                 
C     
         do I=1,NPP                                                             
            DELX(I)=  SQRT(VTLDV(I,I))                                             
         end do
         do I=1,NPP                                                             
            do J=1,NPP                                                          
               VTLDV(I,J)=VTLDV(I,J)/(DELX(I)*DELX(J))                             
            end do
         end do
         IF(IBGVTD .NE.0) THEN                                                     
            WRITE(6,51)                                                            
            CALL MATOUT(VTLDV,maxparam,maxparam,NPP,NPP)                                     
         ENDIF   


C     JMF debug
C     
C     Curious to see the condition of VTLDV
C     kappa = ||A^-1|| * ||A||
C     for the usual operator norm ||A|| = sigma_1,
C     and for a square matrix, this is lambda_max / lambda_min

C     compute eigenvalues
C     call DGEEV(JOBVL, JOBVR, N, A, LDA, WR, WI, VL, LDVL, VR,
C     LDVR, WORK, LWORK, INFO)

C     do I=1,NPP
C     do J=1,NPP
C     UNIT(I,J)=VTLDV(I,J)
C     end do
C     end do


C     call DGEEV('N', 'N', npp, unit, maxparam, wr, wi, 0, 1, 0, 1,
C     1        WORKi, LWORK, ifail)

C     print *, '%%%%%%%%%%%%%%%%%%%%%%%%%%%'
C     print *, 'VTLDV'
         
C     if (ifail .ne. 0) print *, 'error'
         
C     lambda(1) = sqrt(wr(1)*wr(1) + wi(2)*wi(2))
C     lambda(2) = lambda(1)

C     det = 1.0
         
C     do i = 1, npp
C     lmb = sqrt(wr(i)*wr(i) + wi(i)*wi(i))
C     if (lambda(1) .lt. lmb) lambda(1) = lmb
C     if (lambda(2) .gt. lmb) lambda(2) = lmb

C     det = det * lmb
C     end do
         
C     kappa = lambda(1) / lambda(2)
c     print *, 'kappa = ', kappa
C     print *, 'det   = ', det
         
C     print *, '%%%%%%%%%%%%%%%%%%%%%%%%%%%'
C     JMF debug
         

C     TODO JMF
C     do we actually need to invert anything? Check to make sure
C     if we only need to calculate the dx vector, then inversion is NOT
C     required, so we can solve the system MUCH faster using other routines

C     
C     INVERSION OF THE VTLDV MATRIX.FROM THIS POINT,VTLDV                       
C     BECOMES THE INVERSE MATRIX                                                
C     
C     NAG16 CALL F01AAE(VTLDV,130,NPP,UNIT,130,WKSPCE,IFAIL)                         
C     NAG  BELOW ARE TWO AMUR ROUTINES                         
C     CALL F07ADF(NPP,NPP,VTLDV,130,IPIV,INFO)                          
C     CALL F07AJF(NPP,VTLDV,130,IPIV,WKSPCE,130,IFAIL) 
      
                         
C     UNIT=VTLDV
         do i=1,npp
            do j=1,npp
               unit(i,j)=vtldv(i,j)
            end do
         end do

C     print out the derivative matrix
         if ( ibug5 .eq. 4 ) then

            write(6, *) 'VTLDV Matrix'
            do i = 1, npp
               
               write(6, *) 
     &              '{', 
     &              vtldv(i, 1),
     &              (', ', vtldv(i, j), j=2, npp),
     &              '},'
               
            enddo
            write(6, *)
            
         end if

c     lworki=800         
c     two lines below are for calling two routines in lapack                     
c         call dgetrf(npp,npp,unit,maxparam,ipiv,ifail)
c         if ( ifail.ne. 0 ) print *, 'DGETRF ', IFAIL

c         call dgetri(npp,unit,maxparam,ipiv,worki,lwork,ifail)
c         if ( ifail.ne. 0 ) print *, 'DGETRI ', IFAIL
cc ik modif 12dec 2010: to use on PC on has to turn on the dsytrf and dsytri subroutines
C     VTLDV is symmetric
         call dsytrf('U', npp, unit, maxparam, 
     &        ipiv, worki, lwork, ifail)
C
         if ( ifail .ne. 0 ) print *, 'DSYTRF ', ifail
C
         call dsytri('U', npp, unit, maxparam,
     &        ipiv, worki, ifail)
c  ik modif
c         call dpotrf('U', npp, unit, maxparam, ifail)
         if ( ifail.ne. 0 ) print *, 'DPOTRF ', IFAIL
c         call dpotri('U', npp, unit, maxparam, ifail)
         if ( ifail.ne. 0 ) print *, 'DPOTRF ', IFAIL

C     copy the upper diagonal to the lower diagonal
C         do j=1,npp
C            do i=1,j
C               unit(j,i)=unit(i,j)
C            end do
C         end do
         
C         if ( ifail .ne. 0 ) print *, 'DSYTRI ', ifail

         print*,' IFAIL= ',ifail                                                   
         if(ibgvtd .ne.0) then                                                     
            print*,' UNIT AFTER INVERSION...'                                      
            call matout(unit,maxparam,maxparam,npp,npp)                                      
         endif                                                                     
C         do i=1,npp                                                             
C            do j=1,npp                                                          
C               vtldv(i,j)=unit(i,j)/(delx(i)*delx(j))                              
C            end do
C         end do

        do j=1,npp                                                             
            do i=1,j                                                          
               vtldv(i,j)=unit(i,j)/(delx(i)*delx(j))                              
               vtldv(j, i) = vtldv(i, j)
            end do
         end do

         if(ibgvtd .ne.0) then                                                     
            print*,' VTLDV AFTER INVERSION'                                        
            call matout(vtldv,maxparam,maxparam,npp,npp)                                     
         endif                                                                     
c     
c     from this point,delx   becomes the variation of                           
c     the parameters                                                            
c     
         do i=1,npp                                                             
            delx(i)=0.d0                                                           

            if ( nitt .ne. 0 ) then
               do k=1,npp                                                         
                  delx(i)=  vtldv(i,k)*vtldtm (k)+delx(i)                             
               end do
            endif
         end do
C         acc=1.d0                                                                  
         acc = dxmul

C     testing, variable acc
C     try playing around with this until it works
         if ( dxmul .lt. 0.0 ) then
            acc = 1 / ( (s-1)**2 + 1 )
         endif

         if ( unc .eq. 1 ) then
            do kdat=1,ndata 
               do ii=1,npp
                  do jj=1,npp
                     varcov=vtldv(ii,jj)*(s**2)
                     if(iast(kdat).ne.0) then
                        prod=vu(ii,kdat)*varcov*vu(jj,kdat)
                     else
                        prod=vnou(ii,kdat)*varcov*vnou(jj,kdat)
                     endif
                     if(iw(kdat).ge.1000) then
                        prod1=29979.2458**2*prod   
                     else
                        prod1=prod       
                     endif
                     aincer(kdat)=aincer(kdat)+prod1                     
                  end do
               end do
               aincer(kdat)=sqrt(aincer(kdat))
            end do
         end if

C         if ( nitt .eq. 0 ) cycle
         if ( nitt .eq. 0 ) exit


C     another variable ACC
         if (  dxmul .lt. -10.0 ) then
            acc = 1.0

            do i=1,npp                                                             
               deltax(i)=sqrt(vtldv(i,i))*s                                          
               r=  abs(delx(i)/deltax(i))
               if ( 1/r .lt. acc ) acc = 1/r
            end do
         endif

         write(6,*) 'ACC= ',acc                                                    
         write(6,302)                                                              
c     ***** format changes by jth
 302     FORMAT(/,'STATE',2X,'CONSTANT',4X,'OLD VALUE',13X,                        
     1        'VARIATION',13X,'NEW VALUE',13X,'STD.DEV.',12X,
     2        'VARIA/STV',        
     2        2X,'CONST/STV')




         do i=1,npp                                                             
            deltax (i)=sqrt(vtldv(i,i))*s                                          
            if(i.le.np12) then                                                     
               if(i.gt.np(1)) then                                                 
                  ii=2                                                             
                  j=i-np(1)                                                        
               else                                                                
                  ii=1                                                             
                  j=i                                                              
               endif                                                               
               k=ind(ii,j)                                                         
               xx=cst(ii,k)+delx(i)*acc                                            
               r=  delx(i)/deltax (i)                                              
               rjth=xx/deltax(i)                                              
               write(6,301) ii,abc(6*(k-1)+1:6*k),cst(ii,k),delx(i),               
     &              xx,deltax (i),r,rjth                                   
               cst(ii,k)=xx                                                        
               if(irmw.eq.1) cst(2,k)=cst(1,k)
 301           format(i4,2x,a6,4f22.15,f12.7,d11.2)      
            else
               ii=2                                                                
               viben=term+delx(i)*acc                                              
               r=  delx(i)/deltax (i)                                              
               write(6,301) ii,' TERM ',term,delx(i),viben,deltax (i),r            
               term=viben                                                          
            endif                                                                  
            delxst(i)=delx(i)                                                      
         end do

         call flush(6)

         
C******NEW MATERIAL BELOW                                                       
         do  i=1,npp                                                            
            if(i.le.np12) then                                                     
               if(i.gt.np(1)) then                                                 
                  ii=2                                                             
                  j=i-np(1)                                                        
               else                                                                
                  ii=1                                                             
                  j=i                                                              
               endif                                                               
               k=ind(ii,j)                                                         
               write(6,323) abc(6*(k-1)+1:6*k),cst(ii,k)                                 
c     323        format(2x,a6,' = ',f24.17,',')                                      
 323           format(2x, a6, ' = ', e28.21, ',')
            endif                                                                  
         end do
c******newmaterial above                                                        


         write(6,900) s                                                            
 900     format(/,' STANDARD DEVIATION BEFORE= ',f20.6)                            
         deltas=s-sst                                                              
         sst=s                                                                     
         call cntrl(npp,delx,deltax ,cst)                                          
c     
c     call chekptx(ifiles,1)                                                    
c     

      end do

C     print the correlation coefficient matrix
      if ( nitt .ne. 0 ) then
         do i=1,npp                                                             
            if(i.le.np12) then                                                     
               if(i.gt.np(1)) then                                                 
                  ii=2                                                             
                  j=i-np(1)                                                        
               else                                                                
                  ii=1                                                             
                  j=i                                                              
               endif                                                               
               k=ind(ii,j)                                                         
               
               do kk = 1,6
                  fitparam(6*i+1:6*(i+1)) = abc(6*(k-1)+1:6*k)
               end do
            else
               ii=2                                                                
            endif                                                                  
         end do
         
         fitparam(1:6) = 'para  '

         write(6, *) ''
         write(6, *) 'CORRELATION COEFFICIENTS'

         rewind(44)
         write(44, *) '(', npp+1, '(A6, 10X))'
         rewind(44)
         read(44, '(A64)') paramfmt

         write(6, paramfmt) 
     &        ( fitparam((ii)*6+1:(ii+1)*6), ii=0, npp)

         rewind(44)
         write(44, *) '(A6, 10X, ', npp, '(D15.5, 1X))'
         rewind(44)
         read(44, '(A64)') paramfmt
         
         do ii = 1, npp
            write (6, paramfmt) 
     &           fitparam(ii*6+1:(ii+1)*6), 
     &           (correl(ii, jj), jj = 1, npp)
         end do
      end if

      write(6, *) ''


      call flush(6)
C     

      
C     if(iend) 23,24,25                                                         
c     24   sec=tcpu()                                                                
c     write(6,*) ' begin setup linecal   ',sec                                  
      write(6, *) ' BEGIN SETUP '
      call flush(6)
 24   call setup(cst)                                                           
      write(6, *) ' END SETUP '
      call flush(6)

C     print energy
      if ( ibug5 .eq. 1 ) then
         write (6, *)
         write (6, *) ' ENERGY '
         call printener(cst, npp, 3, jmax)
         write (6, *)
         call flush(6)
      endif

C     cos3g, cos6g expectation values
      if ( ibug5 .eq. 3 ) then
         write(6, *)
         write(6, *) 'EXPECTATION VALUES'

         do vt = 1, 9
            do K = -Jmax, Jmax
               do sigma = 0,1
                  call expec(vt, sigma, K)
               end do
            end do
         end do

         write(6, *) 'isig K <1-cos3g> (vt = 0, 1, 2, 3)'
         call printexp(3)
         write(6, *) 'isig K <1-cos6g> (vt = 0, 1, 2, 3)'
         call printexp(6)
         write(6, *)
      endif

c******thenext statement used to be number 24.                                  
      write(6, *) ' BEGIN LINECAL   '
      call flush(6)
      call lineca (cst,npp)                                                     
      write(6,*) ' END LINECAL   '    
      call flush(6)
      end                                                                       



