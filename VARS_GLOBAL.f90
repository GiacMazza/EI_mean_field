MODULE VARS_GLOBAL
  USE DMFT_VECTORS
!  USE DMFT_TIGHT_BINDING
  USE SF_IOTOOLS
  USE SF_PARSE_INPUT
  implicit none
  !
  integer :: Lk,Norb,Nspin,Nso
  real(8),dimension(:),allocatable :: wtk
  complex(8),dimension(:,:,:),allocatable :: Hk 
  real(8),dimension(:,:),allocatable :: Umat_loc,Umat
  real(8) :: Ndens
  real(8) :: beta
  real(8) :: temp

  type(vect2D),dimension(:),allocatable :: k_bz
  
  integer :: comm,rank,ierr
  logical :: master
  character(len=100) :: init_HF,init_Hsb


contains

  subroutine get_global_vars
    !
    call parse_input_variable(Norb,"NORB","input.conf",default=2)
    call parse_input_variable(Nspin,"Nspin","input.conf",default=2)
    call parse_input_variable(Ndens,"Ndens","input.conf",default=1.d0)
    call parse_input_variable(temp,"TEMP","input.conf",default=1.d-2)
    call parse_input_variable(init_HF,"init_HF_file","input.conf",default='delta_HF.conf')
    call parse_input_variable(init_Hsb,"init_Hsb_file","input.conf",default='Hsb.conf')
    if(master) call save_input_file("input.conf")
    beta=1.d0/temp
    !
  end subroutine get_global_vars


END MODULE VARS_GLOBAL
