! --- qrpERROR ---
! FUNCTION: Return relative Frobenius error, || Q R - A P ||_F / || A ||_F.
function qrpError(m,n,k,A,ldA,QR,ldQR,tau,jpvt) result(relErr)
   implicit none

   ! .. Scalar arguments ..
   integer, intent(in) :: m, n, k, ldA, ldQR
   ! m: rows
   ! n: cols
   ! k: factorization rank.

   ! .. Array arguments ..
   real*8, intent(in) :: A(ldA,n), QR(ldQR,n), tau(min(m,n))
   integer, intent(in) :: jpvt(n)
   
   ! .. Scalar result ..
   real*8 :: relErr

   ! .. Local scalars ..
   integer :: lwork, ldErr, j, info !i
   real*8 :: dtmp

   ! .. Local arrays ..
   real*8, allocatable, dimension(:) :: work
   real*8, allocatable, dimension(:,:) :: Err

   ! .. Parameters ..
   real*8, parameter :: negOne=-1.0_8

   ! .. Functions ..
   real*8 :: dlange

   ! ===========================
   ! .. Executable Statements ..

   ldErr=m
   allocate(Err(ldErr,n))
   Err=0.0
   ! Copy R to Err
   call dlacpy('U',m,k,QR,ldQR,Err,ldErr) ! leading k columns are upper tri
   call dlacpy('A',m,n-k,QR(1,k+1),ldQR,Err(1,k+1),ldErr) ! rest is full
   ! Get lwork and allocate
   call dormqr('L','N',m,n,k,QR,ldQR,tau,Err,ldErr,dtmp,-1,info)
   lwork=dtmp
   allocate(work(lwork))

   ! Apply Q to Err
   call dormqr('L','N',m,n,k,QR,ldQR,tau,Err,ldErr,work,lwork,info)
 
   ! Subtract A(:,jpvt)
   if (jpvt(1)==-1) then
!$OMP PARALLEL PRIVATE(j)
!$OMP DO SCHEDULE(DYNAMIC)
      do j=1,n
         call daxpy(m,negOne,A(1,j),1,Err(1,j),1)
      enddo
!$OMP END DO
!$OMP END PARALLEL
   else
!$OMP PARALLEL PRIVATE(j)
!$OMP DO SCHEDULE(DYNAMIC)
      do j=1,n
         call daxpy(m,negOne,A(1,jpvt(j)),1,Err(1,j),1)
      enddo
!$OMP END DO
!$OMP END PARALLEL
   endif

   relErr=dlange('F',m,n,Err,ldErr,dtmp)/dlange('F',m,n,A,ldA,dtmp)
   deallocate(Err,work)

end function qrpError

