
      double precision etrans, epsi, sw070, w, dxmul, en0
      integer ivobs, nobs, kobs, iblk, iw
      integer ndata, nadat, nofit, nfitmw, nfitir, nfimw0
      integer nfimw1, nitt, iast, nfi004, nfi080, nfi100, nfi020
      integer nfi200, nfi045, ngst0, nfi010
      integer nfi070, nfimw2, ikakc, jmax, nfi1m0, nfir10
      integer nfir21, nfimw3
      integer vtmax
      integer nft1, nft2, nft1a, nft1e, nft2a, nft2e
      integer nft00, nft01, nft02, nft11, nft12, nft22
      integer nft00a, nft01a, nft02a, nft11a, nft12a, nft22a
      integer nft00e, nft01e, nft02e, nft11e, nft12e, nft22e
      integer nftw1, nftw2, nftw3, nftw4
      integer nftw1a, nftw2a, nftw3a, nftw4a
      integer nftw1e, nftw2e, nftw3e, nftw4e
      integer nmw00, nmw01, nmw02, nmw11, nmw12, nmw22
      integer nmw00a, nmw01a, nmw02a, nmw11a, nmw12a, nmw22a
      integer nmw00e, nmw01e, nmw02e, nmw11e, nmw12e, nmw22e
      integer nfi004a, nfi080a, nfi100a, nfi020a, nfi1m0a
      integer nfi200a, nfi045a, nfi010a, nfi070a
      integer nfi004e, nfi080e, nfi100e, nfi020e, nfi1m0e
      integer nfi200e, nfi045e, nfi010e, nfi070e
      integer int, unc

      double precision aincer

      common/expdat /  dxmul,      etrans(ndta),ivobs(2,ndta), en0,
     *     nobs(2,ndta),kobs(2,ndta),iblk(2,ndta),iw(ndta),w(ndta),        
     *     ndata,nadat,epsi,nofit,nfitmw,nfitir, nfimw0,
     *     nfimw1, nitt, iast(ndta), nfi004, nfi080, nfi100, nfi020,
     *     nfi200, nfi045, ngst0, nfi010,sw070,nfi070 ,nfimw2,ikakc,
     *     jmax,nfi1m0,               
     &     nfir10,nfir21,nfimw3,
     &     vtmax, nft1, nft2, nft1a, nft1e, nft2a, nft2e,
     &     nft00, nft01, nft02, nft11, nft12, nft22,
     &     nft00a, nft01a, nft02a, nft11a, nft12a, nft22a,
     &     nft00e, nft01e, nft02e, nft11e, nft12e, nft22e,
     &     nftw1, nftw2, nftw3, nftw4,
     &     nftw1a, nftw2a, nftw3a, nftw4a,
     &     nftw1e, nftw2e, nftw3e, nftw4e,
     &     nmw00, nmw01, nmw02, nmw11, nmw12, nmw22,
     &     nmw00a, nmw01a, nmw02a, nmw11a, nmw12a, nmw22a,
     &     nmw00e, nmw01e, nmw02e, nmw11e, nmw12e, nmw22e,
     &     nfi004a, nfi080a, nfi100a, nfi020a, nfi1m0a,
     &     nfi200a, nfi045a, nfi010a, nfi070a,
     &     nfi004e, nfi080e, nfi100e, nfi020e, nfi1m0e,
     &     nfi200e, nfi045e, nfi010e, nfi070e, int, aincer(ndta), unc

