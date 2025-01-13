MODULE VARS_GLOBAL
  !USE DMFT_VECTORS
  !  USE DMFT_TIGHT_BINDING
  USE SF_IOTOOLS
  USE SF_PARSE_INPUT
  USE SF_CONSTANTS
  implicit none
  !
  integer :: Lk,Norb,Nspin,Nso,Lkr,N_cut_off,Lk_aux
  real(8),dimension(:),allocatable :: wtk,wtk_rl
  complex(8),dimension(:,:,:),allocatable :: Hk 
  real(8),dimension(:,:),allocatable :: Umat_loc,Umat
  real(8) :: Ndens
  real(8) :: beta
  real(8) :: temp
  !

  real(8),dimension(3,3) :: Bkinv
  real(8),dimension(3) :: R1,R2,R3
  real(8),dimension(3) :: Bk1,Bk2,Bk3


  real(8),dimension(:,:),allocatable :: kpt_latt
  real(8),dimension(:,:),allocatable :: kpt_latt_aux

  integer,dimension(:),allocatable :: ndegen

  real(8),dimension(:,:),allocatable :: rpt_latt
  integer :: nrpts
  integer,dimension(:,:),allocatable :: irvec

  real(8),dimension(:,:),allocatable :: k_bz,krl
  real(8) :: dk_mesh
  integer(8),dimension(:,:,:),allocatable :: ik2ii,igr2ik,igr2ik_aux
  integer(8),dimension(:,:),allocatable  :: ikrl2ii
  integer(8),dimension(:,:,:),allocatable :: ir_stride

  real(8),dimension(:),allocatable :: kxgrid,kygrid,kzgrid,kxx,kyy,kzz
  real(8),dimension(:),allocatable :: kxgrid_aux,kygrid_aux,kzgrid_aux
  integer(8),dimension(:),allocatable :: ix_aux,iy_aux,iz_aux

  integer :: Nk_x,Nk_y,Nk_z
  integer,allocatable :: ik_stride(:,:),ik_stride_aux(:,:)  
  
  integer :: comm,rank,ierr,mpiERR
  logical :: master
  character(len=100) :: init_HF,init_Hsb

  logical :: whartree,only_LOC,wfock


  integer :: ir0,irL,irR,irU,irD
  integer,dimension(:),allocatable :: ir_mirror

  real(8) :: hbarc_ev_nm,hbar_ev_ps

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
    call parse_input_variable(whartree,"whartree","input.conf",default=.false.)
    call parse_input_variable(wfock,"wfock","input.conf",default=.true.)
    call parse_input_variable(only_loc,"only_loc","input.conf",default=.false.)


    if(master) call save_input_file("input.conf")

    !temp=kelvin_electron_volt_relationship*temp
    temp=kelvin_electron_volt_relationship*temp
    beta=1.d0/temp
    Ndens=Ndens*dble(Nspin)
    hbarc_ev_nm = Planck_constant_in_eV_s/2.d0/pi*speed_of_light_in_vacuum*1.d9
    hbar_ev_ps=Planck_constant_in_eV_s/2.d0/pi*1.d12
    write(*,*) 'beta-temp eV',temp,beta
    ! stop
    !
  end subroutine get_global_vars
  
  subroutine stops(msg)
    character(len=*),optional :: msg
    character(len=200) :: msg_
    
    msg_='stoppati!'
    if(present(msg)) then
       if(len(msg).gt.200) then
          msg_='stoppati!'
       else          
          msg_=msg
       end if
    end if
    call MPI_BARRIER(comm,mpiERR)
    if(master) write(*,*) msg_
    stop
  end subroutine stops
  
  

END MODULE VARS_GLOBAL
