! --- SVD ---
! Truncated SVD
subroutine svd(A0,m,n,rank)

   ! .. Scalar arguments ..
   integer, intent(in) :: m, n, rank

   ! .. Array arguments ..
   real*8, intent(inout) :: A0(m,n)

   ! .. Local scalars ..
   integer :: i, info, LWORK

   ! .. Local arrays ..
   real*8, allocatable, dimension(:,:) :: U, VT
   real*8, allocatable, dimension(:) :: S, WORK

   ! .. Parameters ..
   real*8, parameter :: negOne=-1.0_8, zero=0.0_8, one=1.0_8

   ! ===========================
   ! .. Executable statements ..

   ! Allocate U, S, VT
   allocate(U(m,n),S(n),VT(n,n),WORK(1))

   ! Run dgesvd
   call dgesvd('S','S',m,n,A0,m,S,U,m,VT,n,WORK,-1,info)
   LWORK=WORK(1)
   deallocate(WORK)
   allocate(WORK(LWORK))
   call dgesvd('S','S',m,n,A0,m,S,U,m,VT,n,WORK,LWORK,info)
   !do i=1,10
   !   write(*,'(es9.2)') S(i)
   !enddo

   ! Zero out singular values
   S(rank+1:)=0

   ! Reconstruct A0
   do i=1,n
      call dlascl('G',1,1,one,S(i),m,1,U(1,i),m,info)
   enddo
   call dgemm('n','n',m,n,rank,one,U,m,VT,n,zero,A0,m)
   !call dgesvd('S','S',m,n,A0,m,S,U,m,VT,n,WORK,LWORK,info)

   ! Deallocate
   deallocate(U,S,VT,WORK)

end subroutine svd
