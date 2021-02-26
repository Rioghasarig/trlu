! --- TRQRCP --- ==============================================================
! Truncated low-rank randomized QRCP.
subroutine trqrcp(m,n,k,blkMax,A,ldA,jpvt,Y,ldY,tau,R,ldR,iseed,info)
   implicit none

   ! .. Scalar Arguments ..
   integer, intent(in) :: m, n, k, blkMax, ldA, ldY, ldR
   ! m: rows in A
   ! n: cols in A
   ! k: approximation rank
   ! blkMax: max columns to permute per sample block
   integer, intent(out) :: info
   ! info: error info

   ! .. Array Arguments ..
   real*8, intent(inout) :: A(ldA,n)
   ! A[in]: matrix to factor
   ! A[out]=A[in](:,jpvt) \approx Q R
   integer, intent(out) :: jpvt(n)
   ! jpvt: column pivots;
   real*8, intent(out) :: tau(k), Y(ldY,k), R(ldR,n)
   ! tau: reflector coefficients
   integer, intent(inout) :: iseed(4)
   ! iseed: randomization seed

   ! .. Local Scalars ..
   integer :: blk, l, jc, i, j
   integer :: nThd, maxThd, dline, iline
   real*8 :: dtmp

   ! .. Local work arrays, leading dimensions ..
   real*8, allocatable, dimension(:,:) :: Omega, B, Wt, WS, BW, B0
   integer, allocatable, dimension(:,:) :: iWS
   ! Omega: compression matrix
   ! B: Omega*A
   ! Wt: T'*Y'*A
   ! WS: workspace with column-based processor affinity.
   ! BW: BW'=Tb'*Yb'*B. This should make qrcp(B) slightly faster.
   ! iWS: integer workspace. Also column processor affinity.
   integer :: ldB, ldOmega, ldWt, ldWS, ldiWS, ldBW
   real*8, allocatable, dimension(:) :: Btau, BN
   integer, allocatable, dimension(:) :: Bjpvt
   ! Btau, BN, Bjpvt : B QRCP info.

   ! .. Parameters ..
   integer, parameter :: samplePad=8, numaMax=12, cacheLine=64
   real*8, parameter :: negOne=-1.0, zero=0.0, one=1.0
   logical, parameter :: checkError=.false., timeSubroutines=.false.

   ! .. Functions ..
   integer, external :: omp_get_max_threads
   real*8, external :: sampleError, qrpError, QR_Error

   ! ===========================
   ! .. Executable Statements ..

   ! Check arguments: trqrcp(m,n,k,blkMax,A,ldA,jpvt,Y,ldY,tau,R,ldR,iseed,info)
   info=0
   if (m<0) then
      info=-1
   elseif (n<0) then
      info=-2
   elseif (k<0.or.k>min(m,n)) then
      info=-3
   elseif (blkMax<=0.or.blkMax>k) then
      info=-4
   elseif (ldA<min(m,1)) then
      info=-6
   elseif (ldY<min(m,1)) then
      info=-9
   elseif (ldR<min(k,1)) then
      info=-12
   endif
   if (info/=0) then
      call xerbla('trqrcp',-info)
      return
   endif

   ! Quick return if possible.
   if (m==0.or.n==0.or.k==0) return

   ! cores, threads
   nThd=omp_get_max_threads()
   maxThd=nThd

   ! cachelines
   iline=cacheLine/kind(jpvt(1))
   dline=cacheLine/kind(A(1,1))

   ! array dimensions
   l=blkMax+samplePad  ! rows in compression and sample
   ldOmega=(1+(l-1)/dline)*dline ! round to cacheline
   ldB=ldOmega;
   ldWt=(1+(k-1)/dline)*dline ! Wt=T'Y'A
   ldWS=(1+(blkMax*k-1)/dline)*dline ! each column needs to store Y(:,jc+1:jc+blkMax)' * Y(:,1:jc) = blkMax*(k-blkMax)
   ldiWS=iline
   ldBW=(1+(n-1)/dline)*dline !BW'=Tb'*Yb'*B   ->   BW = B'*Yb*Tb

   !D write(*,'(a,i,a,i,a,i,a,i)') 'm=',m, 'n=',n, 'k=',k, 'l=',l
   !D write(*,'(a,i,a,i,a,i,a,i)') 'dline=',dline,'iline=',iline,'ldOmega=',ldOmega,'ldWt=',ldWt
   !D write(*,'(a,i,a,i,a,i)') 'ldWS=',ldWS,'ldiWS=',ldiWS,'ldBW=',ldBW
   ! Allocate local arrays.
   allocate(Omega(ldOmega,m),B(ldB,n),Bjpvt(n),Btau(blkMax),BN(n),Wt(ldWt,n),WS(ldWS,maxThd),iWS(ldiWS,maxThd),BW(ldBW,blkMax))
   if (checkError) allocate(B0(ldB,n))
   !! If everything is quite successful then we can probably put the omega allocation into Y as long as l<k

   ! Construct original permutation order, zero R.
   do j=1,n
      jpvt(j)=j
      do i=1,k
         R(i,j)=zero !! We could just zero below diagonal.
      enddo
   enddo

   jc=0 ! completed columns

   ! 0: Construct initial compression matrices.
   if (checkError) write(*,'(a)')'getB_2'
   if (timeSubroutines) call tic('getB_2')
   call getB_2(l,m,n,jc,nThd,Omega,ldOmega,A,ldA,B,ldB,iseed)
   if (timeSubroutines) call toc('getB_2')

   do
!D      if (jc==32) then
!D         call printMatrix('A0',m,n,A,ldA)
!D         call printMatrix('Omega0',l,m,Omega,ldOmega)
!D         call printMatrix('B0',l,n,B,ldB)
!D      endif

      blk=min(blkMax,k-jc) ! number of columns to move in this block.
      
      ! A: Check sample error
      if (checkError) then
         dtmp = sampleError(l,m,n,jc,B,ldB,Omega,ldOmega,A,ldA)
         write(*,'(a,i5,a,i3,a,ES)')'jc=',jc,' blk=',blk,' Sample error: ', dtmp
         !D if (jc>32) stop
      endif

      ! B1: Precopy for B factorization error
      if (checkError) call dlacpy('A',l,n-jc,B(1,jc+1),ldB,B0(1,jc+1),ldB)

      ! 1: qrcp(B)
      if (checkError) write(*,'(a)')'qrcpB_3'
      if (timeSubroutines) call tic('qrcpB_3')
      call qrcpB_3(l,n,jc,blk,nThd,B,ldB,Bjpvt,Btau,BN,BW,ldBW,iWS,ldiWS)
      if (timeSubroutines) call toc('qrcpB_3')

      ! B2: Check B factor error. Needs B0.
      if (checkError) then
         do j=1,n-jc
            Bjpvt(jc+j)=Bjpvt(jc+j)-jc
         enddo
         dtmp=qrpError(l,n-jc,blk,B0(1,jc+1),ldB,B(1,jc+1),ldB,Btau,Bjpvt(jc+1))
         write (*,'(a,es9.2)') 'QR(B) error:', dtmp
         do j=1,n-jc
            Bjpvt(jc+j)=Bjpvt(jc+j)+jc
         enddo
      endif

      ! 2: Permute A,R,Wt,jpvt
      if (checkError) write (*,'(a)') 'permuteBlock'
      if (timeSubroutines) call tic('permuteBlock')
      call permuteBlock(m,n,jc,A,ldA,R,ldR,Wt,ldWt,jpvt,Bjpvt)
      if (timeSubroutines) call toc('permuteBlock')

      ! 3,4: Reconstruct A_J in Y_J and make reflectors. T in upper right Y_J, R_JJ in R.
      if (checkError) write (*,'(a)') 'makeReflectors'
      if (timeSubroutines) call tic('makeReflectors')
      call restrictNuma(maxThd,numaMax,nThd)
      call makeReflectors(m,jc,blk,nThd,A,ldA,Y,ldY,Wt,ldWt,R,ldR,WS,ldWS)
      call unrestrictNuma(maxThd,nThd)
      if (timeSubroutines) call toc('makeReflectors')

      ! 5,6: Construct Wt_J: and R_J:
      if (checkError) write (*,'(a)') 'updateR'
      if (timeSubroutines) call tic('updateR')
      call updateR(m,n,jc,blk,nThd,A,ldA,Y,ldY,Wt,ldWt,R,ldR,WS,ldWS)
      if (timeSubroutines) call toc('updateR')

      if (checkError) then
         dtmp = QR_Error(m,n,jc+blk,jc+blk,A,ldA,Y,ldY,R,ldR)
         write (*, '(a,es9.2)') 'R error: ', dtmp
      endif

      if (jc+blk==k) exit !We are done. No need to update the sample.

      ! 7: Update Omega <- Q_b' Omega (I-Q_a*Q_a') if error checking.
      if (checkError) then
         write (*,'(a)') 'updateOmega'
         call updateOmega(l,m,jc,blk,Omega,ldOmega,B,ldB,Btau,Y,ldY)
      endif

!D      if (jc==32) then
!D         write(*,'(a)',advance='no')'Bjpvt=['
!D         do j=1,jc
!D            write(*,'(i4,a)',advance='no')Bjpvt(j),' '
!D         enddo
!D         do j=jc+1,n
!D            write(*,'(i4,a)',advance='no')Bjpvt(j)+jc,' '
!D         enddo
!D         write(*,'(a)')'];'
!D         call printMatrix('A1',m,blk,A(1,jc+1),ldA)
!D         call printMatrix('Y1',m,jc+blk,Y,ldY)
!D         call printMatrix('Ra1',jc+blk,n,R,ldR)
!D         call printMatrix('Rb1',l,n,B,ldB)
!D         call printMatrix('Omega1',l,m,Omega,ldOmega)
!D      endif

      ! 8: Sample update using projection. B(1:blk,jc+blk+1:n) -= ( B(1:blk,jc+1:jc+blk)/R(jc+1:jc+blk,jc+1:jc+blk) ) * R(jc+1:jc+blk,jc+blk+1:n)
      if (checkError) write (*,'(a)') 'updateB_2'
      if (timeSubroutines) call tic('updateB_2')
      call updateB_2(blk,n,jc,nThd,B,ldB,R,ldR)
      if (timeSubroutines) call toc('updateB_2')

!D      if (jc==32) call printMatrix('B2',l,n,B,ldB)

      jc=jc+blk

   enddo

   ! 9: Set tau from diagonal of Y.
   !!! Note that it would be more beneficial to the user to skip tau and just set the upper triangle of Y to T.
   do j=1,k
      tau(j)=Y(j,j)
   enddo

   deallocate(Omega,B,Bjpvt,Btau,BN,Wt,WS,iWS,BW)

   if (checkError) deallocate(B0)

end subroutine trqrcp
! =============================================================================






! --- GETB_2 --- ==============================================================
! Get sample matrix B = Omega A.
subroutine getB_2(l,m,n,jc,nThd,Omega,ldOmega,A,ldA,B,ldB,iseed)

   ! .. Scalar Arguments ..
   integer, intent(in) :: l,m,n,jc,nThd,ldOmega,ldA,ldB

   ! .. Array Arguments ..
   real*8, intent(in) :: A(ldA,n)
   real*8, intent(out) :: Omega(ldOmega,m), B(ldB,n)
   ! Omega: random compression matrix
   ! B: compressed matrix, B=Omega*A

   integer, intent(inout) :: iseed(4)
   ! iseed: integer seed for randomized Omega

   ! .. Local Scalars ..
   integer :: iThd, j1, nP

   ! .. Parameters ..
   real*8, parameter :: zero=0.0_8, one=1.0_8

   ! .. External Functions ..
   integer, external :: omp_get_thread_num

   ! ===========================
   ! .. Executable Statements ..

   if (jc>0) then
      write(*,'(a)') 'Error: Resampling not implemented.'
      ! We have to cleverly remove projections. I'm not sure how to do this efficiently.
      ! Probably project out Q before multiplying A.
      stop
   endif

   ! Construct random compression matrix Omega.
   ! Note: This is non-trivial to parallelize due to causal seed dependency!
   !call tic('genOmega')
   call dlarnv(3, iseed, ldOmega*m, Omega)
   !call toc('genOmega')

   ! Construct compressed matrix. B=Omega*A.
   ! dgemm would autoparallelize, but I want to control memory affinity.
!$OMP PARALLEL PRIVATE(iThd,j1,nP) SHARED(l,m,n,jc,nThd,Omega,A,B,ldOmega,ldA,ldB)
   iThd=omp_get_thread_num()
   call getFullPanel(n,jc,iThd,nThd,j1,nP)
   if (nP>0) then ! PROCESS JP=j1:j1+nP-1
      ! B(1:l,JP)=Omega(1:l,1:m)*A(1:m,JP)
      call dgemm('N','N',l,nP,m,one,Omega,ldOmega,A(1,j1),ldA,zero,B(1,j1),ldB)
      ! B=Omega*(I-Q*Q')*A
      ! B -= (Omega*Q)*R
   endif
!$OMP END PARALLEL

end subroutine getB_2
! =============================================================================






! --- GETFULLPANEL --- ========================================================
subroutine getFullPanel(n,jc,iThd,nThd,j1,nP)

   ! .. Arguments ..
   integer, intent(in) :: n, jc, iThd, nThd
   integer, intent(out) :: j1, nP
   integer, parameter :: blockMin=32

   ! .. Local scalars ..
   integer :: nUsed, jLast

   ! .. Executable statements ..
   nUsed = max(1,min(nThd, (n-jc)/blockMin))

   j1=jc+1+(iThd*(n-jc))/nUsed
   jLast=jc+((iThd+1)*(n-jc))/nUsed
   if (j1>n) then !This thread does not participate in paneling
      j1=n+1
      !jLast=n
      nP=0
   elseif (jLast>n) then !This shouldn't happen
      write(0,'(a)') 'Error: jFirst<=n, jLast>n'
   else
      nP=jLast-j1+1
   endif

end subroutine getFullPanel
! =============================================================================





! --- SAMPLEERROR --- =========================================================
! FUNCTION: Return relative Frobenius error, ||B - Omega(I-Q Q')A|| / ||B||.
!D function sampleError(l,m,n,jc,B,ldB,Omega,ldOmega,Y,ldY,tau,A,ldA) result(relErr)
function sampleError(l,m,n,jc,B,ldB,Omega,ldOmega,A,ldA) result(relErr)
   implicit none

   ! .. Scalar arguments ..
!D   integer, intent(in) :: l, m, n, jc, ldB, ldOmega, ldY, ldA
   integer, intent(in) :: l, m, n, jc, ldB, ldOmega, ldA
   ! l: rows in Omega, B
   ! m: rows in A
   ! n: cols in A
   ! jc: number of reflectors.

   ! .. Array arguments ..
!D   real*8, intent(in) :: B(ldB,n), Omega(ldOmega,m), Y(ldY,jc), A(ldA,n), tau(jc)
   real*8, intent(in) :: B(ldB,n), Omega(ldOmega,m), A(ldA,n)

   ! .. Scalar result ..
   real*8 :: relErr

   ! .. Local scalars ..
!D   integer :: i,j,info
   real*8 :: dtmp

   ! .. Local arrays ..
!D   real*8, allocatable, dimension(:,:) :: Q1, POmega, BErr, work
   real*8, allocatable, dimension(:,:) :: BErr

   ! .. Parameters ..
   real*8, parameter :: negOne=-1.0_8, zero=0.0_8, one=1.0_8

   ! .. Functions ..
   real*8 :: dlange

   ! ===========================
   ! .. Executable Statements ..
!D   allocate(POmega(l,m),BErr(l,n-jc),Q1(m,jc),work(l,jc))
   allocate(BErr(l,n-jc))

   !D I don't think we need the projection code here. Omega should be updated to remove
   !D I believe tau is the problem. It isn't initialized. Remove this projection all together.
   ! Set Q1 to E1
   !do j=1,jc
   !   do i=1,j-1
   !      Q1(i,j)=zero
   !   enddo
   !   Q1(j,j)=one
   !   do i=j+1,m
   !      Q1(i,j)=zero
   !   enddo
   !enddo
   !if (jc>0) then! Apply Q to obtain Q1
   !   call dormqr('L','N',m,jc,jc,Y,ldY,tau,Q1,m,work,l*jc,info)
   !endif

   !D Copy Omega to POmega
   !call dlacpy('A',l,m,Omega,ldOmega,POmega,l)

   !if (jc>0) then
   !   ! Multiply POmega*Q1 -> work
   !   call dgemm('N','N',l,jc,m,one,POmega,l,Q1,m,zero,work,l)
   !   ! POmega-= work*Q1'
   !   call dgemm('N','T',l,m,jc,negOne,work,l,Q1,m,one,POmega,l)
   !endif

   ! Copy B to BErr
   call dlacpy('A',l,n-jc,B(1,jc+1),ldB,BErr,l)
   !D BErr-=POmega*A
   !call dgemm('N','N',l,n-jc,m,negOne,POmega,l,A(1,jc+1),ldA,one,BErr,l)
   ! BErr-=Omega*A
   call dgemm('N','N',l,n-jc,m,negOne,Omega,ldOmega,A(1,jc+1),ldA,one,BErr,l)

   ! Relative error ||BErr||_F / ||B||_F
   relErr=dlange('F',l,n-jc,BErr,l,dtmp)/dlange('F',l,n-jc,B(1,jc+1),ldB,dtmp)
!D   deallocate(POmega,BErr,Q1,work)
   deallocate(BErr)

end function sampleError
! =============================================================================





! --- QRCPB_3 --- =============================================================
subroutine qrcpB_3(l,n,jc,blk,nThd,B,ldB,Bjpvt,Btau,BN,BW,ldBW,iWS,ldiWS)
   implicit none

   ! .. Scalar Arguments ..
   integer, intent(in) :: l,n,jc,blk,nThd,ldB,ldBW,ldiWS
   ! l: rows of B
   ! n: cols of B
   ! jc: cols completed in A. Start factoring at jc+1.
   ! blk: columns to permute to front.

   ! .. Array Arguments ..
   real*8, intent(inout) :: B(ldB,n)
   ! B: sample matrix to factorize
   integer, intent(out) :: Bjpvt(n), iWS(ldiWS,nThd)
   ! Bjpvt: column indices for B pivots.
   ! iWS: small integer workspace for column decisions
   real*8, intent(out) :: Btau(blk), BW(ldBW,blk), BN(n)
   ! Btau : B reflector coefficients
   ! BW'=T'*Y'*B.

   ! .. Local Scalars ..
   integer :: iPrc, jPrc, iThd, j1, nP, i, maxInd
   real*8 :: maxNrm, bjj, dtmp

   ! .. Parameters ..
   real*8, parameter :: negOne=-1.0_8, zero=0.0_8, one=1.0_8
   integer, parameter :: cacheLine=64

   ! .. External Functions ..
   real*8 :: dnrm2
   integer, external :: omp_get_thread_num, omp_get_max_threads

   ! ===========================
   ! .. Executable Statements ..

   if (blk>n-jc) then
      write(*,'(a,i,a,i,a)') 'QRCPB Error: ',blk,' columns requested but only ',n-jc,' columns remain.'
      stop
   endif

   iPrc=0 !Last row processed
   jPrc=jc !Last column processed

!$OMP PARALLEL PRIVATE(iThd,i,j1,nP,maxNrm,dtmp) SHARED(l,n,jc,blk,iPrc,jPrc,maxInd,nThd,iWS,Bjpvt,BN,B,ldB,BW,ldBW,bjj,Btau)
   iThd=omp_get_thread_num()

   ! 0: Get initial norms, column ordering, candidates.
   iWS(1,iThd+1)=0 ! No valid column candidate.
   maxNrm=zero
   call getFullPanel(n,jc,iThd,nThd,j1,nP)
   !D write(*,'(a,i2,a,i4,a,i4)')'T',iThd,' j1=',j1,' nP=',nP
   if (nP>0) then
      ! Process panel j1:j1+nP-1
      do i=j1,j1+nP-1
         Bjpvt(i)=i ! Original permutation ordering
         BN(i)=dnrm2(l,B(1,i),1) ! Compute 2-norm
         if (BN(i)>maxNrm) then
            maxNrm=BN(i)
            iWS(1,iThd+1)=i
         endif
      enddo
   endif

   i=iPrc
   !INVARIANT:
   !  A: i is last row processed. Temp private before iPrc is altered to avoid race condition.
   !  B: BN(:) are correct, iWS(1,:) have max column indices
   !  C: j1:j1+nP have correct affinity with last operation on B.
!$OMP BARRIER

   ! Make branch decisions using i since iPrc could change.
   do while (i<blk)

      ! 1: Construct next reflector.
      !D write (*,'(a,i2,a,i4,a,i4,a,i4,a)')'T',iThd,' ',jc+i,'?[',j1,':',j1+nP-1,']'
      if (jc+i+1>=j1.and.jc+i+1<j1+nP) then ! This thread owns the next column to process.
         !D write(*,'(a,i2)')'In T',iThd
         ! Update loop variables:
         iPrc=iPrc+1 ! Row to process.
         jPrc=jPrc+1 ! Col to process.
         ! Note: jPrc=jc+iPrc, but both are used so much this is better.

         ! Get best candidate:
         maxNrm=0
         maxInd=0
         do i=1,nThd
            if (iWS(1,i)>0.and.BN(iWS(1,i))>maxNrm) then
               maxInd=iWS(1,i)
               maxNrm=BN(maxInd)
            endif
         enddo

         !write(*,'(a,es9.2)')'Pivot norm = ',maxNrm

         if (maxInd==0) then
            write(0,'(a)')'Error: zero pivot.' !!! Probably switch maxInd to jPrc, maxNrm to BN(jPrc).
            !!! Do this after it is working. It prevents an error when trailing matrix is identically zero.
            write(0,'(a,i,a,i,a,i)')'jc=',jc,' iPrc=',iPrc,' jPrc=',jPrc
            stop
         endif

         ! Pivot on B, BW, BN, Bjpvt.
         if (maxInd/=jPrc) then
            ! Pivot B columns.
            call dswap(l,B(1,jPrc),1,B(1,maxInd),1)
            ! Pivot BW rows. BW'=T'*Y'*B.
            call dswap(iPrc-1,BW(jPrc,1),ldBW,BW(maxInd,1),ldBW)
            ! Pivot column norm. New jPrc is never referenced so don't bother.
            BN(maxInd)=BN(jPrc)
            ! Pivot column pivot index
            i=Bjpvt(jPrc)
            Bjpvt(jPrc)=Bjpvt(maxInd)
            Bjpvt(maxInd)=i
         endif

         ! Construct jPrc column which will then become reflector.
         ! B(iPrc:l,jPrc)-=Y(iPrc:l,1:iPrc-1)*BW'(1:iPrc-1,jPrc)
         ! Y(iPrc:l,1:iPrc-1)=B(iPrc:l,jc+1:jc+iPrc-1)
         ! BW'(1:iPrc-1,jPrc)=BW(jPrc,1:iPrc-1)'
         if (iPrc>1) then
            call dgemv('N',l-iPrc+1,iPrc-1,negOne,B(iPrc,jc+1),ldB,BW(jPrc,1),ldBW,one,B(iPrc,jPrc),1)
         endif
         !D write(*,'(a,i2,a,i,a,i,a,i)')'T',iThd,' makeRefl: l=',l,' iPrc=',iPrc,' jPrc=',jPrc
         !D write(*,'(a,i2,a,i,a,i,a,i)')'T',iThd,' jc=',jc, ' ldB=',ldB, ' ldBW=', ldBW
         ! Construct reflector.
         call dlarfg(l-iPrc+1,B(iPrc,jPrc),B(iPrc+1,jPrc),1,Btau(iPrc))
         bjj=B(iPrc,jPrc)
         B(iPrc,jPrc)=one
         ! Construct temporary vector wrk=-Btau*Y(iPrc:l,1:iPrc-1)'*Y(iPrc:l,iPrc)
         ! This is part ofthe BW update formula to add row iPrc to BW'
         ! wrk=BW(1,blk) This should have at least blk-1 rows free.
         if (iPrc>1) then
            call dgemv('T',l-iPrc+1,iPrc-1,-Btau(iPrc),B(iPrc,jc+1),ldB,B(iPrc,jPrc),1,zero,BW(1,blk),1)
         endif
      endif

!$OMP BARRIER
      ! iPrc, jPrc, wrk, reflector in B(iPrc:l,jPrc) must be current before use next.

      ! 2: Apply new reflector to construct new row. Update norms, candidates.
      iWS(1,iThd+1)=0 ! No valid column candidates for next iteration.
      maxNrm=zero
      call getFullPanel(n,jPrc,iThd,nThd,j1,nP)
      if (nP>0) then
         ! Process panel J = j1:j1+nP-1

         ! Construct Btau(iPrc) * Y(iPrc:l,iPrc)' * B(iPrc:l,J) -> BW'(iPrc,J)
         ! Y(iPrc:l,iPrc) = B(iPrc:l,jc+iPrc=jPrc), inc 1
         ! BW'(iPrc,J) = BW(J,iPrc), inc 1
         call dgemv('T',l-iPrc+1,nP,Btau(iPrc),B(iPrc,j1),ldB,B(iPrc,jPrc),1,zero,BW(j1,iPrc),1)

         ! BW'(iPrc,JP)+=wrk*BW'(1:iPrc-1,JP)
         ! BW'(1:iPrc-1,JP)=BW(JP,1:iPrc-1)
         ! The lower left subtriangle of BW' is incorrect, but unreferenced.
         ! wrk=BW(1,blk)
         if (iPrc>1) then
            call dgemv('N',nP,iPrc-1,one,BW(j1,1),ldBW,BW(1,blk),1,one,BW(j1,iPrc),1)
         endif

         ! Update B(iPrc,JP) -= Y(iPrc,1:iPrc)*BW'(1:iPrc,JP)
         ! B(iPrc,JP), inc ldB
         ! Y(iPrc,1:iPrc)=B(iPrc,jc+1:jc+iPrc), inc ldB
         ! BW'(1:iPrc,JP)=BW(JP,1:iPrc)'
         call dgemv('N',nP,iPrc,negOne,BW(j1,1),ldBW,B(iPrc,jc+1),ldB,one,B(iPrc,j1),ldB)

         ! Update column norms if needed.
         if (iPrc<blk) then
            do i=j1,j1+nP-1
               dtmp=abs(B(iPrc,i))/BN(i)
               BN(i)=BN(i)*sqrt(max(zero, (one+dtmp)*(one-dtmp)))
               if (BN(i)>maxNrm) then
                  maxNrm=BN(i)
                  iWS(1,iThd+1)=i
               endif
            enddo
         endif

      endif

      i=iPrc
      !INVARIANT:
      !  A: i is last row processed. Temp private before iPrc is altered to avoid race condition.
      !  B: BN(:) are correct, iWS(1,:) have max column indices
      !  C: j1:j1+nP have correct affinity with last operation on B.
!$OMP BARRIER
      ! Don't update B(iPrc,jPrc)=bjj until reflectors are finished being used.
      ! Don't decide with iPrc or jPrc in this section since they may change.

      ! One thread needs to update B(iPrc,jPrc).
      if (iThd==0) then
         B(i,jc+i)=bjj ! =B(iPrc,jPrc). Note iPrc and jPrc could update before this is set!
      endif

   enddo

   ! There doesn't seem to be a need to wait here. BW may not be finished being
   ! computed, however, the only thread to use it will be the one that computes
   ! it first.

   if (iPrc<l) then
      ! 3: Finish updating trailing matrix.
      call getFullPanel(n,jPrc,iThd,nThd,j1,nP)
      if (nP>0) then
         ! Process panel JP=j1:j1+nP-1
         ! Update B(iPrc+1:l,JP) -= Y(iPrc+1:l,1:iPrc)*BW'(1:iPrc,JP)
         ! Y(iPrc+1:l,1:iPrc) = B(iPrc+1:l,jc+1:jc+iPrc)
         ! BW'(1:iPrc,JP) = BW(JP, 1:iPrc)'
         call dgemm('N','T',l-iPrc,nP,iPrc,negOne,B(iPrc+1,jc+1),ldB,BW(j1,1),ldBW,one,B(iPrc+1,j1),ldB)
      endif
   endif

!$OMP END PARALLEL

end subroutine qrcpB_3
! =============================================================================





! --- PERMUTEBLOCK --- ========================================================
! Apply Bjpvt to A, R, Wt, jpvt
subroutine permuteBlock(m,n,jc,A,ldA,R,ldR,Wt,ldWt,jpvt,Bjpvt)
   implicit none

   ! .. Scalar Arguments ..
   integer, intent(in) :: m, n, jc, ldA, ldR, ldWt

   ! .. Array Arguments ..
   real*8, intent(inout) :: A(ldA,n), R(ldR,n), Wt(ldWt,n)
   integer, intent(inout) :: jpvt(n), Bjpvt(n)

   ! .. Local Scalars ..
   integer :: j, t0, jt, js

   ! ===========================
   ! .. Executable Statements ..

   !! All this nonsense is just to update the column indices of jpvt.
   !! At least it should be fast since it doesn't move much data.
   ! Adjust permutation array.
   do j=jc+1,n
      if (Bjpvt(j)/=j) then
         if (jpvt(j)>0) then

            ! This column has not been processed. Do it now.
            t0 = jpvt(j) ! Store contents of position j. Create hole.
            jt = j ! Target column; This is where the hole is located.
            js = Bjpvt(jt) ! Source column.
            jpvt(jt)=jpvt(js) ! Move source to target.

            do

               jt=js ! Update target/hole position.
               js=Bjpvt(jt) ! Update source.
               if (js==j) exit !When we complete the orbit, finish.
               jpvt(jt)=-jpvt(js) ! Move source to target.

            enddo

            jpvt(jt)=-t0

         else

            ! This column was processed in a prior orbit. Fix sign.
            jpvt(j)=-jpvt(j)

         endif
      endif
   enddo

   ! Modification for having Bjpvt start at jc+1
   do j=jc+1,n
      Bjpvt(j)=Bjpvt(j)-jc
   enddo

   ! Permute columns of A.
   call dlapmt(.true.,m,n-jc,A(1,jc+1),ldA,Bjpvt(jc+1))

   ! Permute columns of R, restricted to nontrivial rows.
   call dlapmt(.true.,jc,n-jc,R(1,jc+1),ldR,Bjpvt(jc+1))

   ! Permute columns of Wt, restricted to nontrivial rows.
   call dlapmt(.true.,jc,n-jc,Wt(1,jc+1),ldWt,Bjpvt(jc+1))

end subroutine permuteBlock
! =============================================================================






! --- MAKEREFLECTORS --- ======================================================
! Connection matrix T is stored in upper triangle: Y(jc+1:jc+blk,jc+1:jc+blk)
! new reflectors with implied unit-heads: Y(jc+1:m,jc+1:jc+blk)
! new upper triangle of R: R(jc+1:jc+blk,jc+1:jc+blk)
! T is needed to efficiently construct Wt for next step.
! tau is unnecessary since it is on the diagonal of T (Y).
subroutine makeReflectors(m,jc,blk,nThd,A,ldA,Y,ldY,Wt,ldWt,R,ldR,WS,ldWS)
   implicit none

   ! .. Scalar Arguments ..
   integer, intent(in) :: m,jc,blk,nThd,ldA,ldY,ldWt,ldR,ldWS
   ! m: rows of A,Y
   ! jc: rows / cols completed
   ! blk: number of new columns to process
   ! nThd: number of threads

   ! .. Array Arguments ..
   real*8, intent(in) :: A(ldA,*), Wt(ldWt,*)
   real*8, intent(inout) :: Y(ldY,*), R(ldR,*)
   real*8, intent(out) :: WS(ldWS,nThd)
   ! Wt(1:jc,remaining) = T'*Y'*A(:,remaining)
   ! Y[out]: lower triangle with implied unit diagonal, reflectors
   !       : upper triangle is connection matrix T; Q=(I - Y*T*Y').
   ! WS: workspace with column thread affinity

   ! .. Local Scalars ..
   integer :: iPrc, jPrc, iThd, i, j1, nP

   ! .. Local Arrays ..

   ! .. Parameters ..
   real*8, parameter :: negOne=-1.0_8, zero=0.0_8, one=1.0_8, two=2.0_8

   ! .. External Functions ..
   integer, external :: omp_get_thread_num, omp_get_max_threads

   ! ===========================
   ! .. Executable Statements ..


   iPrc=0 !Last row/col processed within this block.
   jPrc=jc !Last column in matrix processed. Use for column affinity.

   ! 1: Construct reflected columns in Y which become basis for next reflection block.
!$OMP PARALLEL PRIVATE(iThd,j1,nP,i) SHARED(m,jc,blk,iPrc,jPrc,nThd,A,ldA,Y,ldY,Wt,ldWt,WS,ldWS,R,ldR)
   iThd=omp_get_thread_num()

   ! Y(jc+1:m,jc+1:jc+blk) = A(jc+1:m,jc+1:jc+blk) - Y(jc+1:m,1:jc)*Wt(1:jc,jc+1:jc+blk)
   ! Y(JP,jc+1:jc+blk) = A(JP,jc+1:jc+blk) - Y(JP,1:jc)*Wt
   call getFullPanel(m,jc,iThd,nThd,j1,nP) ! row blocking.
   if (nP>0) then
      call dlacpy('A',nP,blk,A(j1,jc+1),ldA,Y(j1,jc+1),ldY)
      call dgemm('N','N',nP,blk,jc,negOne,Y(j1,1),ldY,Wt(1,jc+1),ldWt,one,Y(j1,jc+1),ldY)
   endif

!$OMP BARRIER
   ! 2: QRF(Ynew)
   ! INVARIENTS:
   !  iPrc gives last row / column processed within the active block.
   !  Columns in Y(jc+iPrc+1:m,jc+iPrc+1:jc+blk) must be current with all prior reflections.
   do while (iPrc<blk)
      ! 2a: Begin computing partial inner products with processing column.
      ! We require y_2'* [Y_1 y_2 Y_3].
      ! y_2[in]'*y_2[in]:
      !     y_2[in] = [alpha; x]
      !     [alpha; x] -> [beta; 0] where beta = -sign(alpha)*sqrt(alpha**2 + x'*x)
      !     y_2[out] = [1; x*gamma] where gamma = 1/(alpha - beta)
      !     t_22 = 2/(y_2[out]'*y_2[out]) = 2/( 1 + x'*x*gamma**2)
      !          = -1/(beta*gamma)
      ! y_2[in]'*Y_1:
      !     T_12 = -T_11*(Y_1'*y_2[out])*t_22
      ! y_2[in]'*Y_3:
      !    Y_3 -= y_2[out]*t_22'*y_2[out]'*Y_3
      call getFullPanel(m,jPrc+1,iThd,nThd,j1,nP) ! row blocking.
      ! This blocking excludes row jPrc+1, which is the current row being processed.
      ! This allows the reflector head to be computed and scaled after.
      ! We adjust inner products accordingly but retain efficient blocking.
      if (nP>0) then
         call dgemv('T',nP,blk,one,Y(j1,jc+1),ldY,Y(j1,jPrc+1),1,zero,WS(1,iThd+1),1)
      endif

      ! 2b: Gather partial results.
      nP=nThd !Each processor privately knows how many copies exist to accumulate
      do while (nP>1)
!$OMP BARRIER
         ! Wait until all active copies in WS are current.
         if (iThd<nP/2) then
            j1=iThd+1+((nP+1)/2) ! Source workspace column
            do i=1,blk
               WS(i,iThd+1)=WS(i,iThd+1)+WS(i,j1)
            enddo
         endif
         nP=(nP+1)/2 !New number of active copies
      enddo

!$OMP BARRIER
      ! WS(1:blk,1) = [Y_1 y_2 Y_3](jc+i+1:m,:)' * y_2(jc+i+1:m)

      ! 2c: Generate reflectors, new column of Y, row of R, row of Wt.
      if (iThd==0) then

         iPrc=iPrc+1 ! Now iPrc is current row / column being processed in block.
         jPrc=jPrc+1 ! jPrc=jc+iPrc. Row / column in Y,R,A,Wt.

         ! beta=-s1*sqrt(alpha**2 + x'*x).
         ! beta**2-alpha**2 = x'*x
         R(jPrc,jPrc)=-sign(sqrt(Y(jPrc,jPrc)**2+WS(iPrc,1)),Y(jPrc,jPrc)) ! beta

         ! gamma = 1/(alpha - beta). y[out]=[1; x*gamma]
         WS(iPrc,1)=one/(Y(jPrc,jPrc)-R(jPrc,jPrc)) ! gamma

         ! tau = 2/(y_2[out]'*y_2[out])
         ! = 2/(1+x'x*gamma**2)
         ! = 2/(1+(beta**2 - alpha**2)/(alpha-beta)**2)
         ! = 2/((alpha-beta-alpha-beta)/(alpha-beta))
         ! = -(alpha-beta)/beta
         ! = -1/(gamma*beta)
         Y(jPrc,jPrc)=-one/(WS(iPrc,1)*R(jPrc,jPrc)) ! tau=t_22

         ! T_12 = -T_11 * (Y_1'*y_2[out]) * t_22
         ! -(Y_1' * [1; x*gamma])*t_22 = -(Y(jPrc,jc+1:jc+iPrc-1)' + WS(1:iPrc-1,1)*gamma)*t_22
         do i=1,iPrc-1
            Y(jc+i,jPrc) = -(Y(jPrc,jc+i) + WS(i,1)*WS(iPrc,1))*Y(jPrc,jPrc)
         enddo

         ! Finish T12 = T11 * (-Y1'*y2*t22)
         call dtrmv('U','N','N',iPrc-1,Y(jc+1,jc+1),ldY,Y(jc+1,jPrc),1)

         ! Complete inner products in WS(iPrc+1:blk,1) for rank 1 update.
         ! WS(iPrc+1:blk,1) = Y_3'*y_2[out]*t_22
         ! = Y_3'*[1; x*gamma]*t_22
         ! = ( Y(jPrc,jc+iPrc+1:jc+blk)' + WS(iPrc+1:blk,1)*gamma)*Y(jPrc,jPrc)
         ! Also finish R(jPrc,jc+iPrc+1:jc+blk) = Y(jPrc,jc+iPrc+1:jc+blk) - WS(iPrc+1:blk,1)'
         do i=iPrc+1,blk
            WS(i,1) = (Y(jPrc,jc+i)+WS(i,1)*WS(iPrc,1))*Y(jPrc,jPrc)
            R(jPrc,jc+i) = Y(jPrc,jc+i) - WS(i,1)
         enddo

      endif

!$OMP BARRIER
      ! Now that reflector inner product is complete, apply rank 1 update to Y_32

      ! 2d: Apply rank-1 update to trailing rows / cols of Y
      call getFullPanel(m,jPrc,iThd,nThd,j1,nP) ! row blocking.
      ! This blocking excludes processed row jPrc, which is already done.
      if (nP>0) then

         ! Finish scaling reflector by gamma:
         do i=j1,j1+nP-1
            Y(i,jPrc)=Y(i,jPrc)*WS(iPrc,1)
         enddo

         ! Y(JP,jc+iPrc+1:jc+blk) -= Y(JP,jPrc) * Wt(jPrc,jc+iPrc+1:jc+blk)
         call dger(nP,blk-iPrc,negOne,Y(j1,jPrc),1,WS(iPrc+1,1),1,Y(j1,jPrc+1),ldY)

      endif

!$OMP BARRIER
   ! INVARIENTS:
   !  iPrc gives last row / column finished being processed within the active block.
   !  Columns in Y(jc+iPrc+1:m,jc+iPrc+1:jc+blk) must be current with all prior reflections.
   !  Don't write to iPrc since it is used for parallel branch decisions.

   enddo

!$OMP END PARALLEL

end subroutine makeReflectors
! =============================================================================





! --- UPDATER --- ============================================================
subroutine updateR(m,n,jc,blk,nThd,A,ldA,Y,ldY,Wt,ldWt,R,ldR,WS,ldWS)
   implicit none

   ! .. Scalar Arguments ..
   integer, intent(in) :: m,n,jc,blk,nThd,ldA,ldY,ldWt,ldR,ldWS

   ! .. Array Arguments ..
   real*8, intent(in) :: A(ldA,n), Y(ldY,*)
   real*8, intent(inout) :: Wt(ldWt,n), R(ldR,n)
   real*8, intent(out) :: WS(ldWS,nThd)
   ! Y(jc+1:m,jc+1:jc+blk) has newest reflectors on ltri and T on upper triangle.
   ! R(jc+1:jc+blk,jc+blk+1:n) is constructed
   ! Wt(1:jc,remaining) = T'*Y'*A(:,remaining)
   ! WS: workspace with column thread affinity

   ! .. Local Scalars ..
   integer :: iThd,j1,nP,i,j

   ! .. Local Arrays ..

   ! .. Parameters ..
   real*8, parameter :: negOne=-1.0_8, zero=0.0_8, one=1.0_8, two=2.0_8

   ! .. External Functions ..
   integer, external :: omp_get_thread_num

   ! ===========================
   ! .. Executable Statements ..

!$OMP PARALLEL PRIVATE(iThd,j1,nP,i,j) SHARED(m,n,jc,blk,nThd,Y,ldY,WS,ldWS,Wt,ldWt,A,ldA,R,ldR)
   iThd=omp_get_thread_num()

   ! 1: Wt(jc+1:jc+blk,1:jc) = Y(jc+1:m,jc+1:jc+blk)' * Y(jc+1:m,1:jc)
   !     = Y(jc+1:jc+blk,jc+1:jc+blk)' * Y(jc+1:jc+blk,1:jc)
   !     + Y(jc+blk+1:m,jc+1:jc+blk)' * Y(jc+blk+1:m,1:jc)
   ! Row panelling on jc+blk+1:m
   call getFullPanel(m,jc+blk,iThd,nThd,j1,nP)
   if (iThd==0) then
      call dlacpy('A',blk,jc,Y(jc+1,1),ldY,WS(1,1),blk)
      call dtrmm('L','L','T','U',blk,jc,one,Y(jc+1,jc+1),ldY,WS(1,1),blk)
      call dgemm('T','N',blk,jc,nP,one,Y(j1,jc+1),ldY,Y(j1,1),ldY,one,WS(1,1),blk)
   elseif (nP>0) then
      call dgemm('T','N',blk,jc,nP,one,Y(j1,jc+1),ldY,Y(j1,1),ldY,zero,WS(1,iThd+1),blk)
   endif

   nP=nThd !Each processor privately knows how many copies exist to accumulate
   do while (nP>1)
!$OMP BARRIER
      ! Wait until all active copies in WS are current.
      if (iThd<nP/2) then
         j1=iThd+1+((nP+1)/2) ! Source workspace column
         do i=1,blk*jc
            WS(i,iThd+1)=WS(i,iThd+1)+WS(i,j1)
         enddo
      endif
      nP=(nP+1)/2 !New number of active copies
   enddo
!$OMP BARRIER

   ! Col panelling on jc+blk+1:n
   call getFullPanel(n,jc+blk,iThd,nThd,j1,nP)
   if (nP>0) then
      ! 2: Wt(jc+1:jc+blk,jc+blk+1:n) = Y(jc+1:m,jc+1:jc+blk)' * A(jc+1:m,jc+blk+1:n)
      call dlacpy('A',blk,nP,A(jc+1,j1),ldA,Wt(jc+1,j1),ldWt)
      call dtrmm('L','L','T','U',blk,nP,one,Y(jc+1,jc+1),ldY,Wt(jc+1,j1),ldWt)
      call dgemm('T','N',blk,nP,m-jc-blk,one,Y(jc+blk+1,jc+1),ldY,A(jc+blk+1,j1),ldA,one,Wt(jc+1,j1),ldWt)
      ! 3: Wt(jc+1:jc+blk,jc+blk+1:n) -= Wt(jc+1:jc+blk,1:jc) * Wt(1:jc,jc+blk+1:n)
      call dgemm('N','N',blk,nP,jc,negOne,WS(1,1),blk,Wt(1,j1),ldWt,one,Wt(jc+1,j1),ldWt)
      ! 4: utri(Y(jc+1:jc+blk,jc+1:jc+blk))' *-> Wt(jc+1:jc+blk,jc+blk+1:n)
      call dtrmm('L','U','T','N',blk,nP,one,Y(jc+1,jc+1),ldY,Wt(jc+1,j1),ldWt)
      ! 5: R(jc+1:jc+blk,jc+blk+1:n) = A(jc+1:jc+blk,jc+blk+1:n) - Y(jc+1:jc+blk,1:jc+blk)*Wt(1:jc+blk,jc+blk+1:n)
      !     = -Y(jc+1:jc+blk,jc+1:jc+blk)*Wt(jc+1:jc+blk,JP)
      !     - Y(jc+1:jc+blk,1:jc)*Wt(1:jc,JP)
      !     + A(jc+1:jc+blk,JP)
      call dlacpy('A',blk,nP,Wt(jc+1,j1),ldWt,R(jc+1,j1),ldR)
      call dtrmm('L','L','N','U',blk,nP,negOne,Y(jc+1,jc+1),ldY,R(jc+1,j1),ldR)
      call dgemm('N','N',blk,nP,jc,negOne,Y(jc+1,1),ldY,Wt(1,j1),ldWt,one,R(jc+1,j1),ldR)
      do j=j1,j1+nP-1
         do i=jc+1,jc+blk
            R(i,j)=R(i,j)+A(i,j)
         enddo
      enddo
   endif

!$OMP END PARALLEL

end subroutine updateR
! =============================================================================






! --- QR_ERROR --- =============================================================
! FUNCTION: Return relative Frobenius error, || Q(:,1:l)'*A - R(1:l,:) ||_F / || R(1:l,:) ||_F
function QR_Error(m,n,k,l,A,ldA,Y,ldY,R,ldR) result(relErr)
   implicit none

   ! .. Scalar arguments ..
   integer, intent(in) :: m, n, k, l, ldA, ldY, ldR
   ! m: rows of A, Q
   ! n: cols of A, R
   ! k: number of reflectors in Q, Y
   !    reflector coefficients tau are on the diagonal of Y.
   ! l: rows of R to reconstruct

   ! .. Array arguments ..
   real*8, intent(in) :: A(ldA,n), Y(ldY,k), R(ldR,n)

   ! .. Scalar result ..
   real*8 :: relErr

   ! .. Local scalars ..
   integer :: ldError, ldQ, i, j, lwork, info
   real*8 :: dtmp

   ! .. Local arrays ..
   real*8, allocatable, dimension(:) :: tau, work
   real*8, allocatable, dimension(:,:) :: Error, Q

   ! .. Parameters ..
   real*8, parameter :: negOne=-1.0_8, zero=0.0_8, one=1.0_8

   ! .. Functions ..
   real*8 :: dlange

   ! ===========================
   ! .. Executable Statements ..

   ! Allocate Error, Q(:,1:l).
   ldError = l
   ldQ = m
   ! Get workspace requirement in dtmp.
      ! Apply reflectors Y, tau => Q
   call dormqr('L','N',m,l,k,Y,ldY,dtmp,Q,ldQ,dtmp,-1,info)
   lwork=dtmp
   allocate(Q(ldQ,l),Error(ldError,n),tau(k),work(lwork))

   ! Initialize Q=I(1:m,1:l) to identity block.
   do j=1,l
      do i=1,j-1
         Q(i,j)=zero
      enddo
      Q(j,j)=one
      do i=j+1,m
         Q(i,j)=zero
      enddo
   enddo

   ! Initialize tau.
   do j=1,k
      tau(j)=Y(j,j)
   enddo

   ! Apply reflectors Y, tau => Q
   call dormqr('L','N',m,l,k,Y,ldY,tau,Q,ldQ,work,lwork,info)

   ! dlacpy R -> Error
   call dlacpy('A',l,n,R,ldR,Error,ldError)

   ! dgemm: Error -= Q' * A
   call dgemm('T','N',l,n,m,negOne,Q,ldQ,A,ldA,one,Error,ldError)

   ! return ||Error||_F / ||R||_F
   relErr=dlange('F',l,n,Error,ldError,dtmp)/dlange('F',l,n,R,ldR,dtmp)
   deallocate(Q,Error,tau,work)

end function QR_Error
! =============================================================================





! --- UPDATEB_2 --- =============================================================
! B(1:blk,jc+blk+1:n) -= ( B(1:blk,jc+1:jc+blk)/R(jc+1:jc+blk,jc+1:jc+blk) ) * R(jc+1:jc+blk,jc+blk+1:n)
subroutine updateB_2(blk,n,jc,nThd,B,ldB,R,ldR)
! B(1:blk,jc+1:jc+blk)/R(...) will be utri. How do we efficienty apply without too much memory cost.
! If ldR is rounded up to integer multiple of blkMax, there will be room in adjacent rows.
! Otherwise produce an error.
   implicit none

   ! .. Scalar Arguments ..
   integer, intent(in) :: blk,n,jc,nThd,ldB,ldR

   ! .. Array Arguments ..
   real*8, intent(inout) :: B(ldB,n), R(ldR,n)
   ! R(jc+blk+1:jc+2*blk,jc+blk+1:n) is used for workspace.

   ! .. Local Scalars ..
   integer :: iThd,j1,nP,i,j

   ! .. Parameters ..
   real*8, parameter :: negOne=-1.0_8, zero=0.0_8, one=1.0_8, two=2.0_8

   ! .. External Functions ..
   integer, external :: omp_get_thread_num

   ! ===========================
   ! .. Executable Statements ..

   if (ldR<jc+2*blk) then
      write(*,'(a)') 'Error: ldR is too small for necessary work.'
      stop
   endif

!$OMP PARALLEL PRIVATE(iThd,j1,nP,i,j) SHARED(n,jc,blk,nThd,R,ldR,B,ldB)
   iThd=omp_get_thread_num()

   ! 1: B(1:blk,jc+1:jc+blk)/R(jc+1:jc+blk,jc+1:jc+blk) -> B(1:blk,jc+1:jc+blk)
   if (iThd==0) then
      ! Way too small to parallelize.
      ! Zero B below diagonal; it contains reflector information.
      do j=1,blk-1
         do i=j+1,blk
            B(i,jc+j)=zero
         enddo
      enddo
      call dtrsm ('R', 'U', 'N', 'N', blk, blk, one, R(jc+1,jc+1), ldR, B(1,jc+1), ldB)
   endif

!$OMP BARRIER

   ! 2: B(1:blk,jc+blk+1:n) -= (R_b1/R_a1)(1:blk,jc+1:jc+blk) * R(jc+1:jc+blk:jc+blk+1:n)
   ! Col panelling on jc+blk+1:n
   call getFullPanel(n,jc+blk,iThd,nThd,j1,nP)
   if (nP>0) then
      ! 2a: Copy R row block into next row block for triangle multiply.
      call dlacpy('A',blk,nP,R(jc+1,j1),ldR,R(jc+blk+1,j1),ldR)
      ! 2b: Apply triangle from above.
      call dtrmm('L','U','N','N',blk,nP,one,B(1,jc+1),ldB,R(jc+blk+1,j1),ldR)
      ! 2c: Subtract result from B.
      do j=j1,j1+nP-1
         do i=1,blk
            B(i,j)=B(i,j)-R(jc+blk+i,j)
         enddo
      enddo
   endif

!$OMP BARRIER

   if (iThd==0) then
      do j=jc+blk+1,jc+2*blk-1
         do i=j+1,jc+2*blk
            R(i,j)=zero
         enddo
      enddo
   endif

!$OMP END PARALLEL

end subroutine updateB_2
! =============================================================================





! --- UPDATEOMEGA --- ========================================================
subroutine updateOmega(l,m,jc,blk,Omega,ldOmega,B,ldB,Btau,Y,ldY)
   implicit none

   ! .. Scalar arguments ..
   integer, intent(in) :: l, m, jc, blk, ldOmega, ldB, ldY
   ! l: rows of Omega
   ! m: cols of Omega, rows of Q
   ! jc: completed columns in Y prior to this block.
   ! blk: number of columns in Q block to project out of Omega.

   ! .. Array arguments ..
   real*8, intent(in) :: B(ldB,jc+blk), Btau(blk), Y(ldY,jc+blk)
   real*8, intent(inout) :: Omega(ldOmega,m)

   ! .. Local scalars ..
   integer :: ldQ, i, j, lwork, info
   real*8 :: dtmp

   ! .. Local arrays ..
   real*8, allocatable, dimension(:) :: tau, work
   real*8, allocatable, dimension(:,:) :: Q

   ! .. Parameters ..
   real*8, parameter :: negOne=-1.0_8, zero=0.0_8, one=1.0_8

   ! ===========================
   ! .. Executable Statements ..

   ! Allocate Q(:,jc+1:jc+blk).
   ldQ = m
   ! Get workspace requirement in dtmp.
   call dormqr('L','N',m,blk,jc+blk,Y,ldY,dtmp,Q,ldQ,dtmp,-1,info)
   lwork=dtmp !Convert double to integer
   lwork=max(lwork,l*blk)
   allocate(Q(ldQ,l),tau(jc+blk),work(lwork))

   ! Rotate Omega = Qb'*Omega
   call dormqr('L','T',l,m,blk,B(1,jc+1),ldB,Btau,Omega,ldOmega,work,lwork,info)

   !D call printMatrix('QbtOmega0',l,m,Omega,ldOmega)

   ! Initialize Q=I(1:m,jc+1:jc+blk) to identity block.
   do j=1,blk
      do i=1,jc+j-1
         Q(i,j)=zero
      enddo
      Q(jc+j,j)=one
      do i=jc+j+1,m
         Q(i,j)=zero
      enddo
   enddo

   ! Initialize tau.
   do j=1,jc+blk
      tau(j)=Y(j,j)
   enddo

   ! Apply reflectors Y, tau => Q
   call dormqr('L','N',m,blk,jc+blk,Y,ldY,tau,Q,ldQ,work,lwork,info)

   !D call printMatrix('QaBlk',m,blk,Q,ldQ)

   ! dgemm: work = Omega*Q
   call dgemm('N','N',l,blk,m,one,Omega,ldOmega,Q,ldQ,zero,work,l)

   ! dgemm: Omega -= work*Q'
   call dgemm('N','T',l,m,blk,negOne,work,l,Q,ldQ,one,Omega,ldOmega)
   deallocate(Q,tau,work)

end subroutine updateOmega
! =============================================================================
