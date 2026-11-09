      character*6 para
      character*1 pr, probs, ref
      integer reflen
      parameter (reflen=40)
      common/chara/para(2,maxparam),pr,probs(2,ndta),
     *     ref(ndta, reflen)
