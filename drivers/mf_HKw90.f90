program officina
  !
  USE SCIFOR
  USE DMFT_TOOLS
  !
  USE VARS_GLOBAL
  USE DMFT_VECTORS
  USE HF
  !
  !USE DMFT_TIGHT_BINDING
  !USE DMFT_MISC
  !USE SF_IOTOOLS
  !
  USE MPI
  !
  implicit none
  integer :: Nk_x,Nk_y,Nk_z

  real(8) :: KKx,KKy
  complex(8),dimension(:,:,:),allocatable :: delta_hf,H_hf,delta_hf_
  real(8) :: mu_fix
  real(8) :: Uintra,Uinter,Delta_CF
  real(8) :: Uq,thop,err_hf
  real(8),dimension(3) :: R1,R2,R3
  real(8),dimension(3) :: Bk1,Bk2,Bk3

  complex(8),dimension(:,:),allocatable :: Hsb

  complex(8),dimension(:),allocatable :: Uft

  complex(8),dimension(:,:),allocatable :: Hloc
  complex(8),dimension(:,:,:),allocatable :: Hk_w90
  complex(8),dimension(:,:,:,:,:),allocatable :: Hk_w90_reshape
  real(8),dimension(:,:),allocatable :: kpts

  integer :: iso,jso,ispin,jspin,iorb,jorb,ik
  integer :: i,j,k,idim
  integer :: Nhf,ihf,unit_err,unit_obs,unit_in,uio
  integer :: flen,iread

  integer,allocatable :: ik_stride(:,:)

  complex(8),dimension(:),allocatable :: obs_loc
  real(8) :: wmix

  character(len=100) :: fileRlattice

  integer,allocatable :: Nkvect(:)

  real(8),dimension(:,:),allocatable :: kpt_latt
  real(8),dimension(:),allocatable :: kx,ky,kz,kx_,ky_,kz_

  !+- START MPI -+!
  call init_MPI()
  comm = MPI_COMM_WORLD
  call StartMsg_MPI(comm)
  rank = get_Rank_MPI(comm)
  master = get_Master_MPI(comm)
  !

  !+- PARSE INPUT DRIVER -+!
  call parse_input_variable(Nk_x,"Nk_x","input.conf",default=10)
  call parse_input_variable(Nk_y,"Nk_y","input.conf",default=10)
  call parse_input_variable(Nk_z,"Nk_z","input.conf",default=10)
  call parse_input_variable(Uintra,"Uintra","input.conf",default=1.d0)
  call parse_input_variable(Uinter,"Uinter","input.conf",default=1.d0)
  call parse_input_variable(thop,"thop","input.conf",default=0.25d0)
  call parse_input_variable(Delta_CF,"delta","input.conf",default=1.d0)
  call parse_input_variable(Nhf,"Nhf","input.conf",default=100)
  call parse_input_variable(wmix,"wmix","input.conf",default=0.5d0)
  call parse_input_variable(fileRlattice,"R_unit_cell","input.conf",default='R_unit_cell.conf')
  !
  call get_global_vars
  !
  Nso=Norb*Nspin
  !
  !
  !+- read primitive cell lattice vectors -+!
  flen=file_length(trim(fileRlattice))
  if(flen.ne.3) stop "error in reading unit cell vectors"
  unit_in=free_unit(); 
  open(unit=unit_in,file=trim(fileRlattice),status="old",action="read")  
  read(unit_in,*) R1(:)
  read(unit_in,*) R2(:)
  read(unit_in,*) R3(:)
  close(unit_in)
  write(*,*) R1
  write(*,*) R2
  write(*,*) R3

  !+- build a monkhorst-pack grid -+!
  call TB_set_ei(R1,R2,R3)
  call TB_get_bk(Bk1,Bk2,Bk3)
  call build_mp_grid(Nk_x,Nk_y,Nk_z)
  Lk=Nk_x*Nk_y*Nk_z
  allocate(kpt_latt(Lk,3),ik_stride(Lk,3))
  ik=0
  do i=1,Nk_x
     do j=1,Nk_y
        do k=1,Nk_z           
           ik=ik+1
           ik_stride(ik,1) = i
           ik_stride(ik,2) = j
           ik_stride(ik,3) = k
           !
           do idim=1,3
              kpt_latt(ik,idim) = kpt_latt(ik,idim) + kxgrid(i)*Bk1(idim)
              kpt_latt(ik,idim) = kpt_latt(ik,idim) + kygrid(j)*Bk2(idim)
              kpt_latt(ik,idim) = kpt_latt(ik,idim) + kzgrid(k)*Bk3(idim)
           end do
           write(300,*) kpt_latt(ik,:)
        end do
     end do
  end do

  uio=free_unit()
  open(unit=uio,file='mp_kgrid.f90')
  do ik=1,Lk
     write(uio,'(10F18.10)') kpt_latt(ik,:)
  end do
  close(uio)
  !+---------------------------------+!  
  allocate(Hk_w90(Nso,Nso,Lk),Hloc(Nso,Nso))
  call TB_hr_to_hk(R1,R2,R3,Hk_w90,Hloc,'TNS_hr.dat',1,6,1,kpt_latt,Hkfile='Hk_grid.out')

  allocate(kx(Nk_x),ky(Nk_x),kz(Nk_x))
  allocate(Hk_w90_reshape(Nso,Nso,Nk_x,Nk_y,Nk_z))
  do ik=1,Lk
     i=ik_stride(ik,1)
     j=ik_stride(ik,2)
     k=ik_stride(ik,3)
     Hk_w90_reshape(:,:,i,j,k) = Hk_w90(:,:,ik)
     kx(i) = kpt_latt(ik,1)
     ky(j) = kpt_latt(ik,2)
     kz(k) = kpt_latt(ik,3)
  end do

  do i=1,Nk_x
     do j=1,Nk_y
        do k=1,Nk_z
           write(200,'(10F18.10)') kx(i),ky(j),kz(k)
        end do
     end do
  end do
  
  ! !+- old(petocchi?) version -+!
  ! allocate(Nkvect(3))
  ! Nkvect(1)=Nk_x
  ! Nkvect(2)=Nk_y
  ! Nkvect(3)=Nk_z
  ! call TB_hr_to_hk(R1,R2,R3,Hk_w90,Hloc,'TNS_hr.dat',1,6,1,Nkvect,kpt_latt=kpt_latt,Hkfile='Hk.out')



  !+- now to plot the bands the interpolation of Hk should be invoked;
  !+- see the other driver

  stop


  ! subroutine hk_from_w90_hr(R1,R2,R3,ham_k,ham_loc,w90_file,Nspin,Norb,Nlat,Nkvec,kpt_latt,Hkfile,Kpointfile)
  !   implicit none
  !   real(8)               ,intent(in)            ::   R1(:),R2(:),R3(:)
  !  complex(8),allocatable,intent(inout)         ::   ham_k(:,:,:)
  !   complex(8),allocatable,intent(inout)         ::   ham_loc(:,:)
  !  character(len=*)      ,intent(in)            ::   w90_file
  !  integer               ,intent(in)            ::   Nspin,Norb,Nlat
  ! integer(4),allocatable,intent(in)            ::   Nkvec(:)             ![Nkx,Nky,Nkz]
  ! real(8)   ,allocatable,intent(out),optional  ::   kpt_latt(:,:)        ![ik,3]
  ! character(len=*)      ,intent(in) ,optional  ::   Hkfile
  ! character(len=*)      ,intent(in) ,optional  ::   Kpointfile
  !



  ! call build_2d_bz(Nx)
  ! allocate(Umat(Nso,Nso),Umat_loc(Nso,Nso))

  ! Umat=0.d0
  ! Umat_loc=0.d0
  ! unit_in=free_unit()
  ! open(unit=unit_in,file='Umat_loc.in')
  ! do ispin=1,Nspin
  !    do iorb=1,Norb
  !       iso=(ispin-1)*Nspin+iorb
  !       do jspin=1,Nspin
  !          do jorb=1,Norb
  !             jso=(jspin-1)*Nspin+jorb
  !             !
  !             if(iorb.eq.jorb) then
  !                Umat(iso,jso) = Uintra
  !                if(ispin.ne.jspin) Umat_loc(iso,jso)=Uintra
  !             else
  !                Umat(iso,jso) = Uinter
  !                Umat_loc(iso,jso) = Uinter
  !             end if
  !             !
  !          end do
  !       end do
  !       write(unit_in,*) Umat_loc(iso,:),'',Umat(iso,:)
  !    end do
  ! end do
  ! close(unit_in)

  ! open(unit=unit_in,file='Hk.in')

  ! allocate(Hk(Nso,Nso,Lk))
  ! Hk=0.d0
  ! do ik=1,Lk
  !    !
  !    iorb=1
  !    do ispin=1,Nspin
  !       iso=(ispin-1)*Nspin+iorb
  !       Hk(iso,iso,ik) = -Delta_CF*0.5d0+2.d0*thop*(cos(k_bz(ik)%x)+cos(k_bz(ik)%y))
  !    end do
  !    iorb=2
  !    do ispin=1,Nspin
  !       iso=(ispin-1)*Nspin+iorb
  !       Hk(iso,iso,ik) = Delta_CF*0.5d0-2.d0*thop*(cos(k_bz(ik)%x)+cos(k_bz(ik)%y))
  !    end do
  !    write(unit_in,'(10F18.10)') k_bz(ik)%x,k_bz(ik)%y,dreal(Hk(1,1,ik)),dreal(Hk(2,2,ik))
  !    !     
  ! end do


contains



  function get_HK(kpoint,N) result(Hk_out)
    implicit none
    real(8),dimension(:) :: kpoint
    integer              :: N
    complex(8),dimension(N,N) :: Hk_out
    integer :: io,jo
    real(8),dimension(Nk_x,Nk_y,Nk_z) :: Hkr
    real(8) :: Hr,dH

    if(size(kpoint).ne.3) stop "kpoint.ne.3"    
    if(N.ne.Nso)  stop "N/=Nso"
    Hk_out=0.d0
    do io=1,Nso
       do jo=io,Nso
          
          
          ! Hkr=dreal(Hk_w90_reshape(io,jo,:,:,:))
          ! Hr=0.d0
          ! !call polin2(kxx,kyy,Hkr,kpoint(1),kpoint(2),Hr,dH)
          ! call polin3(kxx,kyy,Hkr,kpoint(1),kpoint(2),Hr,dH)
          ! Hk_Hf(io,jo) = Hk_Hf(io,jo) + Hr
          ! !
          ! Hkr=dimag(Hhf_reshape(io,jo,:,:))
          ! Hr=0.d0
          ! call polin2(kxx,kyy,Hkr,kpoint(1),kpoint(2),Hr,dH)
          ! Hk_Hf(io,jo) = Hk_Hf(io,jo) + xi*Hr
          ! Hk_Hf(jo,io) = conjg(Hk_Hf(io,jo))
          !
       end do
    end do

    !write(*,*) Hk_Hf(1,1),kpoint,dH

  end function Get_HK



  subroutine build_mp_grid(Nx,Ny,Nz)
    implicit none
    integer :: Nx,Ny,Nz
    integer :: i,ix,iy,ik,unitk
    integer :: Nx_
    real(8) :: kmax
    real(8),dimension(:),allocatable :: kxr,kyr
    !
    !    
    allocate(kxgrid(Nx),kygrid(Ny),kzgrid(Nz))
    kmax=dble(Nx-1)*0.5d0/dble(Nx)
    kxgrid = linspace(-kmax,kmax,Nx)!,mesh=dk_mesh)
    !
    kmax=dble(Ny-1)*0.5d0/dble(Ny)
    kygrid = linspace(-kmax,kmax,Ny)!,mesh=dk_mesh)
    !
    kmax=dble(Nz-1)*0.5d0/dble(Nz)
    kzgrid = linspace(-kmax,kmax,Nz)
    !
    !
    if(mod(Nx,2).eq.0) then
       kxgrid=kxgrid-0.5d0/dble(Nx)
    end if   
    if(mod(Ny,2).eq.0) then
       kygrid=kygrid-0.5d0/dble(Ny)
    end if   
    if(mod(Nz,2).eq.0) then
       kzgrid=kzgrid-0.5d0/dble(Nz)
    end if   


    !
    ! Lk=Nx*Nx
    ! allocate(k_bz(Lk,2),wtk(Lk),ik2ii(Lk,2))
    ! allocate(igr2ik(Nx,Nx))
    ! wtk=1.d0/dble(Lk)
    ! !
    ! KK(1)=2.d0*pi
    ! KK(2)=2.d0*pi
    ! dk_mesh=dk_mesh*2.d0*pi
    ! !
    ! unitk=free_unit()
    ! open(unit=unitk,file='k_points.out')
    ! ik=0
    ! do ix=1,Nx
    !    do iy=1,Nx
    !       ik = ik + 1          
    !       ik2ii(ik,1)=ix
    !       ik2ii(ik,2)=iy
    !       igr2ik(ix,iy) = ik
    !       k_bz(ik,1)=kxgrid(ix)*KK(1)
    !       k_bz(ik,2)=kygrid(iy)*KK(2)
    !       write(unitk,'(10F18.10)') k_bz(ik,1),k_bz(ik,2)
    !    end do
    ! end do
    ! close(unitk)
    ! kxgrid=kxgrid*KK(1)
    ! kygrid=kxgrid*KK(2)
    ! !
    ! !+-Next brillouin zones! -+!
    ! Nx_=Nx*N_cut_off
    ! !
    ! allocate(kxr(Nx_),kyr(Nx_))
    ! kmax=dble(Nx_-1)*0.5d0/dble(Nx)
    ! kxr = linspace(-kmax,kmax,Nx_)
    ! kyr = kxr    
    ! if(mod(Nx,2).eq.0) then
    !    kxr=kxr-0.5d0/dble(Nx)
    !    kyr=kyr-0.5d0/dble(Nx)
    ! end if
    ! !
    ! Lkr=Nx_*Nx_
    ! allocate(krl(Lkr,2),wtk_rl(Lkr),ikrl2ii(Lkr,2))
    ! wtk_rl=1.d0/dble(Lkr)
    ! !
    ! KK(1)=2.d0*pi
    ! KK(2)=2.d0*pi
    ! !
    ! unitk=free_unit()
    ! open(unit=unitk,file='kRL_points.out')
    ! ik=0
    ! do ix=1,Nx_
    !    do iy=1,Nx_
    !       ik = ik + 1          
    !       krl(ik,1)=kxr(ix)*KK(1)
    !       krl(ik,2)=kyr(iy)*KK(2)
    !       ikrl2ii(ik,1)=ix
    !       ikrl2ii(ik,2)=iy
    !       write(unitk,'(10F18.10)') krl(ik,1),krl(ik,2)
    !    end do
    ! end do
    ! close(unitk)
    ! kxr=kxr*KK(1)
    ! kyr=kyr*KK(2)
    ! !

    ! ik=13
    ! write(*,*) krl(ik,:)
    ! call shift_BZ(ik,ix)
    ! write(*,*) k_bz(ix,:)


    ! ik=18*16
    ! write(*,*) krl(ik,:)
    ! call shift_BZ(ik,ix)
    ! write(*,*) k_bz(ix,:)
    
    !stop
  end subroutine build_mp_grid



  ! subroutine get_U_ft
    

  ! contains
  !   !
  !   function U_R
  !     !
  !   end function U_R
  !   !
  ! end subroutine get_U_ft



  
  function Uq_jellium(q)
    USE DMFT_VECTORS  
    USE VARS_GLOBAL
    implicit none
    type(vect2D) :: q
    real(8),dimension(Nso,Nso) :: Uq_jellium
    !
    Uq_jellium = 0.d0
    if(modulo(q).ne.0.d0) then
       !Uq_jellium = Umat_loc + Umat/modulo(q)
       Uq_jellium = Umat/modulo(q)
    end if
    !
  end function Uq_jellium
   
  function Uloc(q)
    USE DMFT_VECTORS  
    USE VARS_GLOBAL
    implicit none
    type(vect2D) :: q
    real(8),dimension(Nso,Nso) :: Uloc
    !
    Uloc = Umat_loc
    !
  end function Uloc
  
  ! subroutine Uq_dens_dens_loc(k,U)
  !   type(vect2d) :: k
  !   real(8),dimension(Nso,Nso) :: U
  !   U = Umat_loc
  !   if(modulo(k).ne.0.d0) then
  !      U = Umat/modulo(k)
  !   end if
  ! end subroutine Uq_dens_dens_loc


  ! subroutine build_2d_bz(Nx)
  !   implicit none
  !   real(8),dimension(:),allocatable :: kxgrid,kygrid
  !   integer :: i,Nx,ix,iy,ik,unitk
  !   real(8) :: kmax
  !   !
  !   !    
  !   allocate(kxgrid(Nx),kygrid(Nx))
  !   kmax=dble(Nx-1)*0.5d0/dble(Nx)
  !   kxgrid = linspace(-kmax,kmax,Nx)
  !   kygrid = kxgrid    
  !   if(mod(Nx,2).eq.0) then
  !      kxgrid=kxgrid-0.5d0/dble(Nx)
  !      kygrid=kygrid-0.5d0/dble(Nx)
  !   end if
  !   !
  !   Lk=Nx*Nx
  !   allocate(k_bz(Lk),wtk(Lk))
  !   wtk=1.d0/dble(Lk)
  !   !
  !   KK%x=2.d0*pi
  !   KK%y=2.d0*pi
  !   !
  !   unitk=free_unit()
  !   open(unit=unitk,file='k_points.out')
  !   ik=0
  !   do ix=1,Nx
  !      do iy=1,Nx
  !         ik = ik + 1
  !         k_bz(ik)%x=kxgrid(ix)*KK%x
  !         k_bz(ik)%y=kygrid(iy)*KK%y
  !         write(unitk,'(10F18.10)') k_bz(ik)%x,k_bz(ik)%y
  !      end do
  !   end do
  !   close(unitk)
  !   !
  ! end subroutine build_2d_Bz



  ! subroutine build_flat_dos  
  !   implicit none
  !   real(8) :: kx
  !   integer :: i,ir
  !   real(8) :: check,R
  !   complex(8) :: tmp
  !   !
  !   !
  !   allocate(epsik(Lk),wtk(Lk))
  !   allocate(vk_orb(2,Lk))

  !   epsik=linspace(-2.d0*thop,2.d0*thop,Lk)
  !   wtk = 1.d0/dble(Lk)       

  !   check=0.d0
  !   do i=1,Lk
  !      ! kx = pi/dble(Lk+1)*dble(i)
  !      ! epsik(i) = !-2.d0*thop*(dble(kx)-dble(Lk)/2.0)/dble(Lk)
  !      vk_orb(1,i) = 2.d0*thop*sin(kx)*alpha_hop
  !      vk_orb(2,i) = 2.d0*thop*sin(kx)       
  !      ! wtk(i) = 1.d0/dble(Lk)       
  !      ! !write(373,'(10F18.10)') kx,1.d0+epsik(i)*(1.d0-alpha_hop)
  !      ! check=check+1.d0/(1.d0+epsik(i)*(1.d0-alpha_hop))*wtk(i)
  !   end do
  !   ! write(373,'(10F18.10)') check
  !   ! call get_free_dos(epsik,wtk,file='DOS_free.kgrid')

  !   do ir=1,Lk
  !      tmp=0.d0
  !      R=dble(ir)-dble((Lk+1)/2)
  !      do i=1,Lk
  !         tmp=tmp+exp(xi*pi/dble(Lk+1)*dble(i)*R)*sqrt(1.0-2.0*thop*(1.d0-alpha_hop)*cos(pi/dble(Lk+1)*dble(i)))
  !      end do
  !      !write(390,'(10F18.10)') R,tmp
  !   end do
  !   !
  ! end subroutine build_flat_dos













end program Officina



!AMOEBA TEST


