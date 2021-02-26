! --- RESTRICTNUMA ---
subroutine restrictNuma(currentMax,numaMax,nThd)

   ! .. Scalar Arguments ..
   integer, intent(out) :: currentMax, nThd
   integer, intent(in) :: numaMax

   ! .. Functions ..
   integer, external :: omp_get_max_threads

   ! .. Subroutines ..
   external :: omp_set_num_threads
   ! ===========================
   ! .. Executable Statements ..

   currentMax=omp_get_max_threads()
   if (currentMax>numaMax) then
      call omp_set_num_threads(numaMax)
      nThd=numaMax
   else
      nThd=currentMax
   endif

   ! ===========================
end subroutine restrictNuma




! --- UNRESTRICTNUMA ---
subroutine unrestrictNuma(priorMax,nThd)

   ! .. Scalar Arguments ..
   integer, intent(in) :: priorMax
   integer, intent(out) :: nThd

   ! .. Subroutines ..
   external :: omp_set_num_threads

   ! ===========================
   ! .. Executable Statements ..

   call omp_set_num_threads(priorMax)
   nThd=priorMax

   ! ===========================
end subroutine unrestrictNuma
