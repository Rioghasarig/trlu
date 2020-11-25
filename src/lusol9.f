      subroutine lu9mod( m, n, beta, v, w,
     $                   lena,lsize, luparm, parmlu,
     $                   a, indc, indr, 
     $                   lv, li, lj, ip, iq,
     $                   lenc, lenr, locc, locr,
     $                   inform )

      implicit           double precision (a-h,o-z)
      integer            luparm(30)
      double precision   parmlu(30), a(lena), v(m), w(n),lv(lsize)
      integer            indc(lena), indr(lena), ip(m), iq(n)
      integer            lenc(n), lenr(m), li(lsize), lj(lsize)
      integer            locc(n), locr(m)
      logical            singlr

      nout   = luparm(1)
      lprint = luparm(2)
      nrank  = luparm(16)
      nrank1 = luparm(17)
      lenl   = 0
      lenu   = luparm(24)
      lrow   = luparm(25)
      small  = parmlu(3)
      utol1  = parmlu(4)
      nrank0 = nrank
   
!     ------------------------------------------------------------------
!     Find the first nonzero in  w  (in pivotal column order).
!     ------------------------------------------------------------------
      do k = 1, n
         kfirst = k
         j      = iq(k)
         if (abs( w(j) ) .gt. small) go to 120
      end do
      go to 900

!     ------------------------------------------------------------------
!     Eliminate any nonzeros in  v  below the trapezoid.
!     ------------------------------------------------------------------
  120 if (nrank .lt. m) then
         jelm   = 0

         call lu9elm( m, n, jelm, v,
     $                lena,lsize, luparm, parmlu,
     $                lenl, lenu, lrow, nrank1,
     $                a, indc, indr, lv, li, lj,
     $                ip, iq, lenr, locc, locr,
     $                inform, diag )

         if (inform .eq. 1) nrank1 = nrank1+1;
         if (inform .eq. 7) go to 970

      end if

!     ------------------------------------------------------------------
!     Find the last nonzero in  v  (in pivotal row order).
!     ------------------------------------------------------------------
      do k = nrank, 1, -1
         klast  = k
         i      = ip(k)
         if (abs( v(i) ) .gt. small) go to 220
      end do
      go to 900

!     ------------------------------------------------------------------
!     Perform a backward sweep of eliminations to reduce part of  v
!     to a multiple of the unit vector  e(iw),  where  iw = ip(klast).
!     Elements  ip(kfirst+1),  ip(kfirst+2),  ...,  ip(klast)  of  v
!     are involved.
!     L, U  and  ip  are updated accordingly.
!     U  will then be trapezoidal except for row  iw = ip(klast).
!     ------------------------------------------------------------------
  220 if (kfirst+1 .lt.  klast) then
         call lu9bak( m, n, kfirst, klast, v,
     $                lena, lsize,luparm, parmlu,
     $                lenl, lenu, lrow,
     $                a, indc, indr, lv,li,lj,
     $                ip, iq, lenr, locc, locr,
     $                inform )
         if (inform .ne. 0) go to 970
      end if
!     ------------------------------------------------------------------
!     Pack the nonzeros of  w  in pivotal order in front of  L.
!     (We will treat the packed  w  much like a normal row of  U.)
!     Set markers on  w  and initialize the corresponding
!     elements of  indc(*),  which are later used by  lu7asv.
!     ------------------------------------------------------------------

      lfree  = lena - lenl
      minfre = n + 1 - kfirst
      nfree  = lfree - lrow
      if (nfree .ge. minfre) go to 310
      call lu1rec( m, .true., luparm, lrow, lena, a, indr, lenr, locr )
      nfree  = lfree - lrow
      if (nfree .lt. minfre) go to 970

  310 lw     = lfree + 1

      do 320   k  = n, kfirst, -1
         j        = iq(k)
         if (abs( w(j) ) .le. small) go to 320
         lw       = lw - 1
         a(lw)    = w(j)
         indr(lw) = j
         indc(lw) = 0
         locc(j)  = lw
  320 continue

      lw1    = lw
      lw2    = lfree
      lenw   = lw2 + 1 - lw1
      lfree  = lfree - lenw

!     ------------------------------------------------------------------
!     Add multiples of  w  to the first  kfirst  rows of  U.
!     (This does not alter the trapezoidal form of  U.)
!     ------------------------------------------------------------------

      do 450 k  = 1, kfirst
         iv     = ip(k)
         cv     = v(iv)
         if (abs( cv ) .le. small) go to 450

!        ===============================================================
!        Compress storage if necessary, so there will be room if
!        row  iv  has to be moved to the end.
!        ===============================================================
         minfre = n
         nfree  = lfree - lrow
         if (nfree .ge. minfre) go to 420
         call lu1rec( m, .true., luparm, lrow, lena, a,indr,lenr,locr )
         nfree  = lfree - lrow
         if (nfree .lt. minfre) go to 970

!        ===============================================================
!        Set  v  =  v  +  wmult * w.
!        ===============================================================
  420    wmult  = beta * cv
         call lu7asv( m, n, iv, lenw, lw1, lw2, k, wmult,
     $                lena, luparm, parmlu,
     $                lenu, lrow,
     $                a, indc, indr, lenr, locc, locr )
  450 continue

!     ------------------------------------------------------------------
!     Add a multiple of  w  to row  iw  of  U.
!     ------------------------------------------------------------------

      if (kfirst .lt. klast) then
         minfre = n
         nfree  = lfree - lrow
         if (nfree .ge. minfre) go to 500
         call lu1rec( m, .true., luparm, lrow, lena, a,indr,lenr,locr )
         nfree  = lfree - lrow
         if (nfree .lt. minfre) go to 970

  500    iw     = ip(klast)
         marker = m + 1
         wmult  = beta * v(iw)
         call lu7asv( m, n, iw, lenw, lw1, lw2, marker, wmult,
     $                lena, luparm, parmlu,
     $                lenu, lrow,
     $                a, indc, indr, lenr, locc, locr )
      end if


!    -----------------------------------------------------------------
!    Add a multiple of w to rows nrank+1,...,nrank1 of U
!    -----------------------------------------------------------------
      do 460 k  = nrank+1,nrank1
         iv     = ip(k)
         cv     = v(iv)
         if (abs( cv ) .le. small) go to 460

!        ===============================================================
!        Compress storage if necessary, so there will be room if
!        row  iv  has to be moved to the end.
!        ===============================================================
         minfre = n
         nfree  = lfree - lrow
         if (nfree .ge. minfre) go to 430
         call lu1rec( m, .true., luparm, lrow, lena, a,indr,lenr,locr )
         nfree  = lfree - lrow
         if (nfree .lt. minfre) go to 970

!        ===============================================================
!        Set  v  =  v  +  wmult * w.
!        ===============================================================
  430    wmult  = beta * cv
         call lu7asv( m, n, iv, lenw, lw1, lw2, k, wmult,
     $                lena, luparm, parmlu,
     $                lenu, lrow,
     $                a, indc, indr, lenr, locc, locr )
  460 continue

!     ------------------------------------------------------------------
!     Cancel the markers on  w.
!     ------------------------------------------------------------------
      do lw  = lw1, lw2
         jw       = indr(lw)
         locc(jw) = 0
      end do
!     ------------------------------------------------------------------
!     Apply a forward sweep to eliminate the nonzeros in row  iw.
!     ------------------------------------------------------------------

      if (kfirst .gt. klast) then
         if (klast .lt. nrank) go to 900
      else
         call lu9for( m, n, kfirst, klast,
     $                lena, lsize,luparm, parmlu,
     $                lenl, lenu, lrow,
     $                a, indc, indr, lv,li,lj,
     $                ip, iq, lenr, locc, locr,
     $                inform, diag )
         if (inform .eq. 7) go to 970
      end if


!     ------------------------------------------------------------------
!     Set inform for exit.
!     ------------------------------------------------------------------

  900 if (nrank .eq. nrank0) then
         inform =  0
      else if (nrank .lt. nrank0) then
         inform = -1
      else
         inform =  1
      end if
      go to 990

      ! Not enough storage.

  970 inform = 7
      if (nout. gt. 0  .and.  lprint .ge. 0)
     &     write(nout, 1700) lena

      ! Exit.

  990 luparm(10) = inform
      luparm(15) = luparm(15) + 1
      luparm(16) = nrank
      luparm(17) = nrank1
      luparm(23) = lenl      !lenL <- lenl
      luparm(24) = lenu
      luparm(25) = lrow
      return

 1700 format(/ ' lu8mod  error...  Insufficient storage.',
     $         '    lena =', i8)

      end ! subroutine lu9mod
      

      subroutine lu9elm( m, n, jelm, v,
     $                   lena, lsize, luparm, parmlu,
     $                   lenL, lenU, lrow, nrank1,
     $                   a, indc, indr, lv, li, lj,
     $                   ip, iq, lenr, locc, locr,
     $                   inform, diag )

      implicit           double precision (a-h,o-z)
      integer            luparm(30)
      double precision   parmlu(30), a(lena), v(m), lv(lsize)
      integer            indc(lena), indr(lena), ip(m), iq(n), lenr(m)
      integer            locc(n), locr(m), li(lsize), lj(lsize)

      !-----------------------------------------------------------------
      ! lu9elm  eliminates the subdiagonal elements of a vector  v(*),
      ! where  L*v = y  for some vector y.
      ! If  jelm > 0,  y  has just become column  jelm  of the matrix  A.
      ! lu9elm  should not be called unless  m  is greater than  nrank.
      !
      ! inform = 0 if y contained no subdiagonal nonzeros to eliminate.
      ! inform = 1 if y contained at least one nontrivial subdiagonal.
      ! inform = 7 if there is insufficient storage.
      !
      ! 09 May 1988: First f77 version.
      !              No longer calls lu9for at end.  lu8rpc, lu8mod do so.
      ! 20 Dec 2015: ilast is now output by lu1rec.
      !-----------------------------------------------------------------

      integer          jw, lw, lw1, lw2, lenw
      parameter        ( zero = 0.0d+0 )

      small  = parmlu(3)
      nrank2 = nrank1 + 1
      diag   = zero


      ! Pack the subdiagonals of  v  into  L,  and find the largest.

  100 vmax   = zero
      kmax   = 0
      l      = lsize - lenL + 1

      do 200 k = nrank2, m
         i       = ip(k)
         vi      = abs( v(i) )
         if (vi .le. small) go to 200
         l       = l - 1
         lv(l)    = v(i)
         li(l) = i
         if (vmax .ge. vi ) go to 200
         vmax    = vi
         kmax    = k
         lmax    = l
  200 continue

      if (kmax .eq. 0) go to 900

      !-----------------------------------------------------------------
      ! Remove  vmax  by overwriting it with the last packed  v(i).
      ! Then set the multipliers in  L  for the other elements.
      !-----------------------------------------------------------------
      imax       = ip(kmax)
      vmax       = lv(lmax)
      lv(lmax)   = lv(l)
      li(lmax)   = li(l)
      l1         = l + 1
      l2         = lsize - lenL
      lenL       = lenL + (l2 - l)
       
      
      do 300 l = l1, l2
         lv(l)    = - lv(l) / vmax
         lj(l) =   imax
  300 continue
      
      ! Move the row containing vmax to pivotal position nrank + 1.

      ip(kmax  ) = ip(nrank2)
      ip(nrank2) = imax
      diag       = vmax

      !-----------------------------------------------------------------
      ! If jelm is positive, insert  vmax  into a new row of  U.
      ! This is now the only subdiagonal element.
      !-----------------------------------------------------------------

      if (jelm .gt. 0) then
         lrow       = lrow + 1
         locr(imax) = lrow
         lenr(imax) = 1
         a(lrow)    = vmax
         indr(lrow) = jelm
      end if

      inform = 1
      go to 990

      ! No elements to eliminate.

  900 inform = 0
      go to 990

      ! Not enough storage.

  970 inform = 7

  990 return

      end ! subroutine lu9elm

      subroutine lu9bak( m, n, kfirst, klast, c,
     $                   lena, lsize, luparm, parmlu,
     $                   lenl, lenu, lrow,
     $                   a, indc, indr, lv, li, lj, 
     $                   ip, iq, lenr, locc, locr,
     $                   inform )

      implicit           double precision (a-h,o-z)
      integer            luparm(30)
      double precision   parmlu(30), a(lena), c(m), lv(lsize)
      integer            indc(lena), indr(lena), ip(m), iq(n), lenr(m)
      integer            locc(n), locr(m), li(lsize), lj(lsize)

!     ------------------------------------------------------------------
!     lu9bak  (Backward Sweep) updates the LU factorization  A = L*U
!     when part of the column vector  c  is reduced to a multiple of
!     a unit vector by a backward sweep of stabilized row operations
!     of the form
!
!        LSWEEP * c  =  c(iw) * e(iw).
!
!     Elements  ip(kfirst),  ip(kfirst+2),  ...,  ip(klast)  of  c
!     are involved.  It is assumed that  (kfirst + 1) < klast.
!     The matrices  L  and  U  are updated to become
!
!        L(new)  =  L * LSWEEP(inverse)    and    U(new)  =  LSWEEP * U
!
!     where  ip * U(new) * iq  is upper triangular except for the
!     row  iw = ip(klast).  The row permutation  ip  is updated
!     accordingly.  The column permutation  iq  is not altered.
!
!     Row  iw  is held in compact form like the other rows of  U.
!
!     On entry, all elements of  locc  are assumed to be zero.
!     On a successful exit (inform = 0) they will again be zero.
!
!     -- Feb 1985: First version.
!     10 May 1988: First F77 version.
!     ------------------------------------------------------------------

      double precision   lmax
      logical            cmprss, first, swap
      
      nrank  = luparm(15)
      lmax   = parmlu(2)
      small  = parmlu(3)
      lfree  = lsize - lenl
      iw     = ip(klast)
      cw     = c(iw)
      lenw   = lenr(iw)

!     ==================================================================
!     Eliminate the nonzeros of  c  backwards.
!     ==================================================================
      first  = .true.

      do 500 k  = klast-1, kfirst+1, -1
         iv     = ip(k)
         cv     = c(iv)
         if (abs( cv ) .le. small) go to 500

         ! cv  has to be eliminated.
         
         cmprss = .false.
         swap   = .false.

!        ---------------------------------------------------------------
!        Compress storage if necessary, so there will be room if
!        v  has to be moved to the end.  (If we waited until  v  is
!        actually moved, the markers on  w  would have to be reset.)
!        ---------------------------------------------------------------
         minfre = n + 1
         nfree  = lena - lrow
         if (nfree .ge. minfre) go to 200
         call lu1rec( m, .true., luparm, lrow, lena, a,indr,lenr,locr )
         cmprss = .true.
         nfree  = lena - lrow
         if (nfree .lt. minfre) go to 970

!        ---------------------------------------------------------------
!        We are about to set  v  =  v  +  ( - cv/cw ) * w.
!        See if rows  v  and  w  should be interchanged first.
!        ---------------------------------------------------------------
  200    if (klast .gt. nrank)  go to 210 
         if (lmax * abs( cv )  .lt.  abs( cw )) go to 210
         if (lmax * abs( cw )  .lt.  abs( cv )) go to 220
         if (lenr(iv)          .lt.  lenw     ) go to 220

  210    if (first             .or. cmprss ) go to 250
         go to 300

!        ---------------------------------------------------------------
!        Interchange rows  iv  and  iw.
!        ---------------------------------------------------------------
  220    swap      = .true.
         i         = iv
         iv        = iw
         iw        = i
         w         = cv
         cv        = cw
         cw        = w
         ip(k)     = iv
         ip(klast) = iw
         if (first  .or.  lenw .eq. 0) go to 250

         ! Cancel the markers on the old row  iw.

         do 230 l = lw1, lw2
            jw       = indr(l)
            locc(jw) = 0
  230    continue

!        ---------------------------------------------------------------
!        Set markers on the new row  iw.
!        This happens if  first,  cmprss,  or  swap.
!        ---------------------------------------------------------------
  250    lenw   = lenr(iw)
         if (lenw .eq. 0) go to 300
         lw1    = locr(iw)
         lw2    = lw1 + lenw - 1

         do 260 lw = lw1, lw2
            jw       = indr(lw)
            locc(jw) = lw
            indc(lw) = 0
  260    continue

!        ---------------------------------------------------------------
!        Form the multiplier and store it in the L file.
!        ---------------------------------------------------------------
  300    first       = .false.
         wmult       = - cv / cw
         lv(lfree)   = wmult
         lj(lfree)   = iw
         li(lfree)   = iv
         lenl        = lenl  + 1
         lfree       = lfree - 1

!        ---------------------------------------------------------------
!        Set  v  =  v  +  wmult * w.
!        ---------------------------------------------------------------
         if (lenw .gt. 0) then
            call lu7asv( m, n, iv, lenw, lw1, lw2, k, wmult,
     $                   lena, luparm, parmlu,
     $                   lenu, lrow,
     $                   a, indc, indr, lenr, locc, locr )
         end if

!        ---------------------------------------------------------------
!        If there was a swap, make sure the diagonal of  v  is first.
!        ---------------------------------------------------------------
         if (.not. swap ) go to 500
         lenv   = lenr(iv)
         if (lenv .eq. 0) go to 500
         if ( k   .gt. n) go to 500
         jdiag  = iq(k)
         lv1    = locr(iv)
         lv2    = lv1 + lenv - 1

         do l = lv1, lv2
            if (indr(l) .eq. jdiag) go to 470
         end do
         go to 500

  470    indr(l)   = indr(lv1)
         indr(lv1) = jdiag
         diag      = a(l)
         a(l)      = a(lv1)
         a(lv1)    = diag
  500 continue

!     ==================================================================
!     End of main elimination loop.
!     ==================================================================

      ! Cancel the markers on  w,  unless they were never set.

      if (first  .or.  lenw .eq. 0) go to 900
      do lw = lw1, lw2
         jw       = indr(lw)
         locc(jw) = 0
      end do

      ! Successful update.

  900 inform = 0
      go to 990

      ! Not enough storage.

  970 inform = 7

      ! Exit.

  990 return

      end ! subroutine lu9bak

      subroutine lu9for( m, n, kfirst, klast,
     $                   lena, lsize, luparm, parmlu,
     $                   lenL, lenU, lrow,
     $                   a, indc, indr, lv,li,lj, 
     $                   ip, iq, lenr, locc, locr,
     $                   inform, diag )

      implicit           double precision (a-h,o-z)
      integer            luparm(30)
      double precision   parmlu(30), a(lena), lv(lsize)
      integer            indc(lena), indr(lena), ip(m), iq(n), lenr(m)
      integer            locc(n), locr(m), li(lsize), lj(lsize)

      !-----------------------------------------------------------------
      ! lu9for  (forward sweep) updates the LU factorization  A = L*U
      ! when row  iw = ip(klast)  of  U  is eliminated by a forward
      ! sweep of stabilized row operations, leaving  ip * U * iq  upper
      ! triangular.
      !
      ! The row permutation  ip  is updated to preserve stability and/or
      ! sparsity.  The column permutation  iq  is not altered.
      !
      ! kfirst  is such that row  ip(kfirst)  is the first row involved
      ! in eliminating row  iw.  (Hence,  kfirst  marks the first nonzero
      ! in row  iw  in pivotal order.)  If  kfirst  is unknown it may be
      ! input as  1.
      !
      ! klast   is such that row  ip(klast)  is the row being eliminated.
      ! klast   is not altered.
      !
      ! lu9for  should be called only if  kfirst .le. klast.
      ! If  kfirst = klast,  there are no nonzeros to eliminate, but the
      ! diagonal element of row  ip(klast)  may need to be moved to the
      ! front of the row.
      !
      ! On entry,  locc(*)  must be zero.
      !
      ! On exit:
      ! inform = 0  if row iw has a nonzero diagonal (could be small).
      ! inform = 1  if row iw has no diagonal.
      ! inform = 7  if there is not enough storage to finish the update.
      !
      ! On a successful exit (inform le 1),  locc(*)  will again be zero.
      !
      !    Jan 1985: Final f66 version.
      ! 09 May 1988: First f77 version.
      ! 20 Dec 2015: ilast is now output by lu1rec.
      ! ------------------------------------------------------------------

      parameter        ( zero = 0.0d+0 )

      double precision   Ltol
      logical            swappd
      
      nrank  = luparm(16)
      Ltol   = parmlu(2)
      small  = parmlu(3)
      uspace = parmlu(6)
      kbegin = kfirst
      swappd = .false.

      ! We come back here from below if a row interchange is performed.

  100 iw     = ip(klast)
      lenw   = lenr(iw)
      if (lenw   .eq.   0  ) go to 910
      lw1    = locr(iw)
      lw2    = lw1 + lenw - 1
      jfirst = iq(kbegin)
      if (kbegin .ge. klast) go to 700

      ! Make sure there is room at the end of the row file
      ! in case row  iw  is moved there and fills in completely.

      minfre = n + 1
      nfree  = lena - lrow
      if (nfree .lt. minfre) then
         call lu1rec( m, .true., luparm, lrow, ilast, lena,
     $                a, indr, lenr, locr )
         lw1    = locr(iw)
         lw2    = lw1 + lenw - 1
         nfree  = lena - lrow
         if (nfree .lt. minfre) go to 970
      end if

      ! Set markers on row  iw.

      do 120 l = lw1, lw2
         j       = indr(l)
         locc(j) = l
  120 continue


      !=================================================================
      ! Main elimination loop.
      !=================================================================
      kstart = kbegin
      kstop  = min( klast, n )

      do 500 k  = kstart, kstop
         jfirst = iq(k)
         lfirst = locc(jfirst)
         if (lfirst .eq. 0) go to 490

         ! Row  iw  has its first element in column  jfirst.

         wj     = a(lfirst)
         if (k .eq. klast) go to 490

         !--------------------------------------------------------------
         ! We are about to use the first element of row iv
         ! to eliminate the first element of row  iw.
         ! However, we may wish to interchange the rows instead,
         ! to preserve stability and/or sparsity.
         !--------------------------------------------------------------
         iv     = ip(k)
         lenv   = lenr(iv)
         lv1    = locr(iv)
         vj     = zero

         if (lenv      .eq.   0   ) go to 150
         if (indr(lv1) .ne. jfirst) go to 150
         vj     = a(lv1)
         if (            swappd               ) go to 200
         if (Ltol * abs( wj )  .lt.  abs( vj )) go to 200
         if (Ltol * abs( vj )  .lt.  abs( wj )) go to 150
         if (            lenv  .le.  lenw     ) go to 200

         !--------------------------------------------------------------
         ! Interchange rows  iv  and  iw.
         !--------------------------------------------------------------
  150    ip(klast) = iv
         ip(k)     = iw
         kbegin    = k
         swappd    = .true.
         go to 600

         !--------------------------------------------------------------
         ! Delete the eliminated element from row  iw
         ! by overwriting it with the last element.
         !--------------------------------------------------------------
  200    a(lfirst)    = a(lw2)
         jlast        = indr(lw2)
         indr(lfirst) = jlast
         indr(lw2)    = 0
         locc(jlast)  = lfirst
         locc(jfirst) = 0
         lenw         = lenw - 1
         lenU         = lenU - 1
         if (lrow .eq. lw2) lrow = lrow - 1
         lw2          = lw2  - 1

         !--------------------------------------------------------------
         ! Form the multiplier and store it in the  L  file.
         !--------------------------------------------------------------
         if (abs( wj ) .le. small) go to 490
         amult   = - wj / vj
         l       = lsize - lenL
         lv(l)   = amult
         lj(l)   = iv
         li(l)   = iw
         lenL    = lenL + 1

         !--------------------------------------------------------------
         ! Add the appropriate multiple of row  iv  to row  iw.
         ! We use two different inner loops.  The first one is for the
         ! case where row  iw  is not at the end of storage.
         !--------------------------------------------------------------
         if (lenv .eq. 1) go to 490
         lv2    = lv1 + 1
         lv3    = lv1 + lenv - 1
         if (lw2 .eq. lrow) go to 400

         !..............................................................
         ! This inner loop will be interrupted only if
         ! fill-in occurs enough to bump into the next row.
         !..............................................................
         do 350 lv0 = lv2, lv3
            jv     = indr(lv0)
            lw     = locc(jv)
            if (lw .gt. 0) then   ! No fill-in.
               a(lw)  = a(lw)  +  amult * a(lv0)
               if (abs( a(lw) ) .le. small) then ! Delete small element
                  a(lw)     = a(lw2)
                  j         = indr(lw2)
                  indr(lw)  = j
                  indr(lw2) = 0
                  locc(j)   = lw
                  locc(jv)  = 0
                  lenU      = lenU - 1
                  lenw      = lenw - 1
                  lw2       = lw2  - 1
               end if
            else  ! Row iw doesn't have an element in column jv yet
                  ! so there is a fill-in.
               if (indr(lw2+1) .ne. 0) go to 360
               lenU      = lenU + 1
               lenw      = lenw + 1
               lw2       = lw2  + 1
               a(lw2)    = amult * a(lv0)
               indr(lw2) = jv
               locc(jv)  = lw2
            end if
  350    continue

         go to 490

         ! Fill-in interrupted the previous loop.
         ! Move row iw to the end of the row file.

  360    lv2      = lv0
         locr(iw) = lrow + 1

         do 370 l = lw1, lw2
            lrow       = lrow + 1
            a(lrow)    = a(l)
            j          = indr(l)
            indr(l)    = 0
            indr(lrow) = j
            locc(j)    = lrow
  370    continue

         lw1    = locr(iw)
         lw2    = lrow

         !..............................................................
         ! Inner loop with row iw at the end of storage.
         !..............................................................
  400    do 450 lv0 = lv2, lv3
            jv     = indr(lv0)
            lw     = locc(jv)
            if (lw .gt. 0) then   ! No fill-in.
               a(lw)  = a(lw)  +  amult * a(lv0)
               if (abs( a(lw) ) .le. small) then ! Delete small element.
                  a(lw)     = a(lw2)
                  j         = indr(lw2)
                  indr(lw)  = j
                  indr(lw2) = 0
                  locc(j)   = lw
                  locc(jv)  = 0
                  lenU      = lenU - 1
                  lenw      = lenw - 1
                  lw2       = lw2  - 1
               end if
            else   ! Row iw doesn't have an element in column jv yet
                   ! so there is a fill-in.
               lenU      = lenU + 1
               lenw      = lenw + 1
               lw2       = lw2  + 1
               a(lw2)    = amult * a(lv0)
               indr(lw2) = jv
               locc(jv)  = lw2
            end if
  450    continue

         lrow   = lw2

         ! The k-th element of row iw has been processed.
         ! Reset swappd before looking at the next element.

  490    swappd = .false.
  500 continue

      !=================================================================
      ! End of main elimination loop.
      !=================================================================

      ! Cancel markers on row iw.

  600 lenr(iw) = lenw
      if (lenw .eq. 0) go to 910
      do 620 l = lw1, lw2
         j       = indr(l)
         locc(j) = 0
  620 continue

      ! Move the diagonal element to the front of row iw.
      ! At this stage, lenw gt 0 and klast le n.

  700 do 720 l = lw1, lw2
         ldiag = l
         if (indr(l) .eq. jfirst) go to 730
  720 continue
      go to 910

  730 diag        = a(ldiag)
      a(ldiag)    = a(lw1)
      a(lw1)      = diag
      indr(ldiag) = indr(lw1)
      indr(lw1)   = jfirst

      ! If an interchange is needed, repeat from the beginning with the
      ! new row iw, knowing that the opposite interchange cannot occur.

      if (swappd) go to 100
      inform = 0
      go to 950

      ! Singular.

  910 diag   = zero
      inform = 1

      ! Force a compression if the file for U is much longer than the
      ! no. of nonzeros in U (i.e. if lrow is much bigger than lenU).
      ! This should prevent memory fragmentation when there is far more
      ! memory than necessary (i.e. when lena is huge).

  950 limit  = int(uspace*real(lenU)) + m + n + 1000
      if (lrow .gt. limit) then
         call lu1rec( m, .true., luparm, lrow, ilast, lena,
     $                a, indr, lenr, locr )
      end if
      go to 990

      ! Not enough storage.

  970 inform = 7

      ! Exit.

  990 return

      end ! subroutine lu9for

      subroutine lu9clr( m, n,
     $                   lena, lsize,luparm, parmlu, 
     $                   a, indc, indr, lv, li, lj,
     $                   ip, iq, lenr, locc, locr,
     $                   c, inform)

      implicit           double precision (a-h,o-z)
      integer            luparm(30)
      double precision   parmlu(30), a(lena), c(m), lv(lsize)
      integer            indc(lena), indr(lena), ip(m), iq(n), lenr(m)
      integer            locc(n), locr(m), li(lsize), lj(lsize)


      parameter        ( zero = 0.0d+0 )

      
      nrank  = luparm(16)
      nrank1 = luparm(17)
      if (nrank1 .eq. nrank) go to 990
      lenL   = 0
      lenU   = luparm(24)
      lrow   = luparm(25)
      small  = parmlu(3)
      uspace = parmlu(6)
      kbegin = 1;
 

      iw     = ip(nrank1)
      lenw   = lenr(iw)
      if (lenw   .eq.   0  ) go to 900
      lw1    = locr(iw)
      lw2    = lw1 + lenw - 1
      jfirst = iq(kbegin)

      ! Make sure there is room at the end of the row file
      ! in case row  iw  is moved there and fills in completely.

      minfre = n + 1
      nfree  = lena - lrow
      if (nfree .lt. minfre) then
         call lu1rec( m, .true., luparm, lrow, ilast, lena,
     $                a, indr, lenr, locr )
         lw1    = locr(iw)
         lw2    = lw1 + lenw - 1
         nfree  = lena - lrow
         if (nfree .lt. minfre) go to 970
      end if

      ! Set markers on row  iw.

      do 120 l = lw1, lw2
         j       = indr(l)
         locc(j) = l
  120 continue


      !=================================================================
      ! Main elimination loop.
      !=================================================================
     

      do 500 k  = 1, nrank
         jfirst = iq(k)
         lfirst = locc(jfirst)
         if (lfirst .eq. 0) go to 500

         ! Row  iw  has its first element in column  jfirst.

         wj     = a(lfirst)
        
         iv     = ip(k)
         lenv   = lenr(iv)
         lv1    = locr(iv)
         vj     = zero

         if (lenv      .eq.  0    ) go to 980
         if (indr(lv1) .ne. jfirst) go to 980
         vj = a(lv1)
         !--------------------------------------------------------------
         ! Delete the eliminated element from row  iw
         ! by overwriting it with the last element.
         !--------------------------------------------------------------
         a(lfirst)    = a(lw2)
         jlast        = indr(lw2)
         indr(lfirst) = jlast
         indr(lw2)    = 0
         locc(jlast)  = lfirst
         locc(jfirst) = 0
         lenw         = lenw - 1
         lenU         = lenU - 1
         if (lrow .eq. lw2) lrow = lrow - 1
         lw2          = lw2  - 1

         !--------------------------------------------------------------
         ! Form the multiplier and store it in the  L  file.
         !--------------------------------------------------------------
         if (abs( wj ) .le. small) go to 500
         amult   = - wj / vj
         l       = lsize - lenL
         lv(l)   = amult
         lj(l)   = iv
         li(l)   = iw
         lenL    = lenL + 1

         !--------------------------------------------------------------
         ! Add the appropriate multiple of row  iv  to row  iw.
         ! We use two different inner loops.  The first one is for the
         ! case where row  iw  is not at the end of storage.
         !--------------------------------------------------------------
         if (lenv .eq. 1) go to 500
         lv2    = lv1 + 1
         lv3    = lv1 + lenv - 1
         if (lw2 .eq. lrow) go to 400

         !..............................................................
         ! This inner loop will be interrupted only if
         ! fill-in occurs enough to bump into the next row.
         !..............................................................
         do 350 lv0 = lv2, lv3
            jv     = indr(lv0)
            lw     = locc(jv)
            if (lw .gt. 0) then   ! No fill-in.
               a(lw)  = a(lw)  +  amult * a(lv0)
               if (abs( a(lw) ) .le. small) then ! Delete small element
                  a(lw)     = a(lw2)
                  j         = indr(lw2)
                  indr(lw)  = j
                  indr(lw2) = 0
                  locc(j)   = lw
                  locc(jv)  = 0
                  lenU      = lenU - 1
                  lenw      = lenw - 1
                  lw2       = lw2  - 1
               end if
            else  ! Row iw doesn't have an element in column jv yet
                  ! so there is a fill-in.
               if (indr(lw2+1) .ne. 0) go to 360
               lenU      = lenU + 1
               lenw      = lenw + 1
               lw2       = lw2  + 1
               a(lw2)    = amult * a(lv0)
               indr(lw2) = jv
               locc(jv)  = lw2
            end if
  350    continue

         go to 500

         ! Fill-in interrupted the previous loop.
         ! Move row iw to the end of the row file.

  360    lv2      = lv0
         locr(iw) = lrow + 1

         do 370 l = lw1, lw2
            lrow       = lrow + 1
            a(lrow)    = a(l)
            j          = indr(l)
            indr(l)    = 0
            indr(lrow) = j
            locc(j)    = lrow
  370    continue

         lw1    = locr(iw)
         lw2    = lrow

         !..............................................................
         ! Inner loop with row iw at the end of storage.
         !..............................................................
  400    do 450 lv0 = lv2, lv3
            jv     = indr(lv0)
            lw     = locc(jv)
            if (lw .gt. 0) then   ! No fill-in.
               a(lw)  = a(lw)  +  amult * a(lv0)
               if (abs( a(lw) ) .le. small) then ! Delete small element.
                  a(lw)     = a(lw2)
                  j         = indr(lw2)
                  indr(lw)  = j
                  indr(lw2) = 0
                  locc(j)   = lw
                  locc(jv)  = 0
                  lenU      = lenU - 1
                  lenw      = lenw - 1
                  lw2       = lw2  - 1
               end if
            else   ! Row iw doesn't have an element in column jv yet
                   ! so there is a fill-in.
               lenU      = lenU + 1
               lenw      = lenw + 1
               lw2       = lw2  + 1
               a(lw2)    = amult * a(lv0)
               indr(lw2) = jv
               locc(jv)  = lw2
            end if
  450    continue

         lrow   = lw2

         ! The k-th element of row iw has been processed.
         ! Reset swappd before looking at the next element.

  
  500 continue

      !=================================================================
      ! End of main elimination loop.
      !=================================================================

      ! Cancel markers on row iw.

      lenr(iw) = 0
      do 620 l = lw1, lw2
         j       = indr(l)
         locc(j) = 0
  620 continue

      !================================================================== 
      ! Clear the extra column of L 
      !==================================================================
     
      ! Compress row file if necessary.

      minfre = m - nrank1
      nfree  = lena  - lrow
      if (nfree .ge. minfre) go to 630
      call lu1rec( m, .true., luparm, lrow, ilast,
     &             lena, a, indr, lenr, locr )
      nfree  = lena - lrow
      if (nfree .lt. minfre) go to 970

  630 imax = ip(nrank1)
      vmax = c(imax) 
      l = lsize - lenL
      do 700 k = nrank1+1, m
        i = ip(k)
        ci = abs(c(i))
        if (ci .le. small) go to 700 
        lv(l) = c(i)/vmax
        lj(l) = imax
        li(l) = i
        l = l - 1;
        lenL = lenL + 1;
  700 continue
  900 inform = 0
      luparm(17) = nrank1 - 1; 
      luparm(23) = lenL; 
      luparm(24) = lenU; 
      go to 950

 

      ! Force a compression if the file for U is much longer than the
      ! no. of nonzeros in U (i.e. if lrow is much bigger than lenU).
      ! This should prevent memory fragmentation when there is far more
      ! memory than necessary (i.e. when lena is huge).

  950 limit  = int(uspace*real(lenU)) + m + n + 1000
      if (lrow .gt. limit) then
         call lu1rec( m, .true., luparm, lrow, ilast, lena,
     $                a, indr, lenr, locr )
      end if
      go to 990

      ! Not enough storage.

  970 inform = 7

      ! Singular U11
  980 inform = 8
      ! Exit.
 
  990 return

      end !e

