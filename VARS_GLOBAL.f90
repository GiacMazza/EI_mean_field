MODULE VARS_GLOBAL
  implicit none
  !
  integer :: Lk,Norb,Nspin
  real(8),dimension(:),allocatable :: wtk
  real(8),dimension(:,:),allocatable :: Hk 
  real(8),dimension(:,:),allocatable :: Umat
  real(8) :: Ndens


END MODULE VARS_GLOBAL
