MODULE VARS_GLOBAL
  !USE DMFT_VECTORS
!  USE DMFT_TIGHT_BINDING
  USE SF_IOTOOLS
  USE SF_PARSE_INPUT
  USE SF_CONSTANTS
  implicit none
  !
  integer :: Lk,Norb,Nspin,Nso,Lkr,N_cut_off
  real(8),dimension(:),allocatable :: wtk,wtk_rl
  complex(8),dimension(:,:,:),allocatable :: Hk 
  real(8),dimension(:,:),allocatable :: Umat_loc,Umat
  real(8) :: Ndens
  real(8) :: beta
  real(8) :: temp
  !

  real(8),dimension(:,:),allocatable :: kpt_latt

  real(8),dimension(:,:),allocatable :: k_bz,krl
  real(8) :: dk_mesh
  integer(8),dimension(:,:,:),allocatable :: ik2ii,igr2ik,igr2ik_rlm
  integer(8),dimension(:,:),allocatable  :: ikrl2ii
  !integer(8),dimension(:,:),allocatable :: ik_stride

  real(8),dimension(:),allocatable :: kxgrid,kygrid,kzgrid,kxx,kyy,kzz
  integer :: Nk_x,Nk_y,Nk_z



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
    call parse_input_variable(N_cut_off,"Nk_cutoff","input.conf",default=4)
    if(master) call save_input_file("input.conf")

    !temp=kelvin_electron_volt_relationship*temp
    temp=kelvin_electron_volt_relationship*temp
    beta=1.d0/temp
    !
  end subroutine get_global_vars


END MODULE VARS_GLOBAL
