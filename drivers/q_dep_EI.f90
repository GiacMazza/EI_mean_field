program officina
  USE SCIFOR
  USE DMFT_TOOLS
  !
  USE VARS_GLOBAL
  USE HF
  !
  USE MPI
  !
  implicit none

  !type(vect2D) :: KK


  real(8) :: KKx,KKy
  complex(8),dimension(:,:,:),allocatable :: delta_hf,H_hf,delta_hf_,Hhf_grid
  complex(8),dimension(:,:,:,:),allocatable :: Hhf_reshape
  real(8) :: mu_fix
  real(8) :: Uintra,Uinter,Delta_CF
  real(8) :: Uq,thop,err_hf
  complex(8),dimension(:,:),allocatable :: Hsb,Hkint
  complex(8),dimension(:,:,:),allocatable :: Hsbk

  complex(8),dimension(:),allocatable :: Uft

  integer :: iso,jso,ispin,jspin,iorb,jorb,ik,ix,iy,ii,jj,i,j
  integer :: Nhf,ihf,unit_err,unit_obs,unit_in

  complex(8),dimension(:),allocatable :: obs_loc
  real(8) :: wmix

  real(8),dimension(2) :: ktest_int
  real(8) :: kint,dk

  !+- START MPI -+!
  call init_MPI()
  comm = MPI_COMM_WORLD
  call StartMsg_MPI(comm)
  rank = get_Rank_MPI(comm)
  master = get_Master_MPI(comm)
  !

  !+- PARSE INPUT DRIVER -+!
  call parse_input_variable(Nx,"Nx","input.conf",default=10)
  call parse_input_variable(Uintra,"Uintra","input.conf",default=1.d0)
  call parse_input_variable(Uinter,"Uinter","input.conf",default=1.d0)
  call parse_input_variable(thop,"thop","input.conf",default=0.25d0)
  call parse_input_variable(Delta_CF,"delta","input.conf",default=1.d0)
  call parse_input_variable(Nhf,"Nhf","input.conf",default=100)
  call parse_input_variable(wmix,"wmix","input.conf",default=0.5d0)

  call get_global_vars
  !
  Nso=Norb*Nspin
  !
  !
  call build_2d_bz(Nx)!;stop
  allocate(Umat(Nso,Nso),Umat_loc(Nso,Nso))
  
  Umat=0.d0
  Umat_loc=0.d0
  unit_in=free_unit()
  open(unit=unit_in,file='Umat_loc.in')
  do ispin=1,Nspin
     do iorb=1,Norb
        iso=(ispin-1)*Nspin+iorb
        do jspin=1,Nspin
           do jorb=1,Norb
              jso=(jspin-1)*Nspin+jorb
              !
              if(iorb.eq.jorb) then
                 Umat(iso,jso) = Uintra
                 if(ispin.ne.jspin) Umat_loc(iso,jso)=Uintra
              else
                 Umat(iso,jso) = Uinter
                 Umat_loc(iso,jso) = Uinter
              end if
              !
           end do
        end do
        write(unit_in,*) Umat_loc(iso,:),'',Umat(iso,:)
     end do
  end do
  close(unit_in)

  open(unit=unit_in,file='Hk.in')

  allocate(Hk(Nso,Nso,Lk))
  Hk=0.d0
  do ik=1,Lk
     !
     iorb=1
     do ispin=1,Nspin
        iso=(ispin-1)*Nspin+iorb
        Hk(iso,iso,ik) = -Delta_CF*0.5d0+2.d0*thop*(cos(k_bz(ik,1))+cos(k_bz(ik,2)))
     end do
     iorb=2
     do ispin=1,Nspin
        iso=(ispin-1)*Nspin+iorb
        Hk(iso,iso,ik) = Delta_CF*0.5d0-2.d0*thop*(cos(k_bz(ik,1))+cos(k_bz(ik,2)))
     end do
     write(unit_in,'(10F18.10)') k_bz(ik,1),k_bz(ik,2),dreal(Hk(1,1,ik)),dreal(Hk(2,2,ik))
     !     
  end do
  close(unit_in)

  !+- loop hartree-fock -+!
  allocate(delta_hf(Nso,Nso,Lk),H_hf(Nso,Nso,Lk),delta_hf_(Nso,Nso,Lk))
  !
  call init_var_params(delta_hf)
  delta_hf_=delta_hf

  unit_err=free_unit()
  open(unit=unit_err,file='err_q.err')
  unit_obs=free_unit()
  open(unit=unit_obs,file='obs_hf_q.loop')
  err_hf=1.d0  
  do ihf=1,Nhf
     write(unit_err,'(10F18.10)') dble(ihf),err_hf,mu_fix          
     call local_single_particle_observables(delta_hf,obs_loc)
     write(unit_obs,'(20F18.10)') dble(ihf),dreal(obs_loc)
     call build_HF_hamiltonian(H_hf,delta_hf,Uq_jellium)     
     
     call find_chem_pot(H_hf,delta_hf,mu_fix)
     delta_hf=wmix*delta_hf+(1.d0-wmix)*delta_hf_
     !
     err_hf=check_conv(delta_hf,delta_hf_)     
 
     !
     delta_hf_=delta_hf
     
  end do
  close(unit_err)
  close(unit_obs)
  !
  call store_HF_hamiltonian_BZgrid(Hhf_grid,delta_hf,mu_fix,Uq_jellium)
  !
  !
  !
  call init_var_params(delta_hf)
  delta_hf_=delta_hf
  !
  unit_err=free_unit()
  open(unit=unit_err,file='err_loc.err')
  unit_obs=free_unit()
  open(unit=unit_obs,file='obs_hf_loc.loop')
  err_hf=1.d0  
  do ihf=1,Nhf
     write(unit_err,'(10F18.10)') dble(ihf),err_hf,mu_fix          
     call local_single_particle_observables(delta_hf,obs_loc)
     write(unit_obs,'(20F18.10)') dble(ihf),dreal(obs_loc)
     call build_HF_hamiltonian(H_hf,delta_hf,Uloc)     
     !
     call find_chem_pot(H_hf,delta_hf,mu_fix)
     delta_hf=wmix*delta_hf+(1.d0-wmix)*delta_hf_
     !
     err_hf=check_conv(delta_hf,delta_hf_)      
     !
     delta_hf_=delta_hf
     !
  end do
  close(unit_err)
  close(unit_obs)

  call save_array('delta.out',delta_hf)  
  call store_HF_hamiltonian_BZgrid(Hhf_grid,delta_hf,mu_fix,Uloc)
  ! allocate(Hhf_reshape(Nso,Nso,3*Nx,3*Nx))
  ! allocate(kxx(3*Nx),kyy(3*Nx))
  ! do ii=1,3
  !    do ix=1,Nx
  !       i=(ii-1)*Nx+ix
  !       kxx(i) = kxgrid(ix)+(ii-2)*2*pi
  !    end do
  ! end do
  ! kyy=kxx
  ! !
  ! do ii=1,3
  !    do jj=1,3
  !       do ik=1,Lk
  !          ix=ik2ii(ik,1)
  !          iy=ik2ii(ik,2)
  !          i=(ii-1)*Nx+ix
  !          j=(jj-1)*Nx+iy
  !          Hhf_reshape(:,:,i,j) = Hhf_grid(:,:,ik)     
  !       end do
  !    end do
  ! end do

  ! do ii=1,3*Nx
  !    do jj=1,3*Nx
  !       write(200,'(10F18.10)') kxx(ii),kyy(jj),Hhf_reshape(1,1,ii,jj)
  !    end do
  ! end do

  ! !+- TEST inperpolate
  ! do ik=1,Lk
  !    write(204,'(10F18.10)') k_bz(ik,:),dreal(Hhf_grid(1,1,ik)),dreal(Hhf_grid(2,2,ik)),dreal(Hhf_grid(1,2,ik))
  ! end do
 
  ! !+- interpolation -+!
  ! allocate(Hkint(Nso,Nso))
  ! dk=pi/dble(50)
  ! kint=-dk
  ! do ix=1,50
  !    kint=kint+dk
  !    ktest_int(1) = kint
  !    ktest_int(2) = kint
  !    Hkint=Hk_HF(ktest_int,Nso)
  !    write(205,'(10F18.10)') ktest_int(:),dreal(Hkint(1,1)),dreal(Hkint(2,2)),dreal(Hkint(1,2))
  ! end do

  ! do ix=1,Nx
  !    ik=igr2ik(ix,ix)
  !    write(206,'(10F18.10)') k_bz(ik,:),dreal(Hhf_grid(1,1,ik)) 
  ! end do


  
  !+- reshape -+!
  !allocate(Hhf_reshape(Nso,Nso,3*Nx,3*Nx))
  
!  allocate(kxx(3*Nx),kyy(3*Nx))
  ! do ii=1,3
  !    do ix=1,Nx
  !       i=(ii-1)*Nx+ix
  !       kxx(i) = kxgrid(ix)+(ii-2)*2*pi
  !       !        write(*,*) kxx(i)
  !    end do
  ! end do
  ! kyy=kxx
  ! !
  ! do ii=1,3
  !    do jj=1,3
  !       do ik=1,Lk
  !          ix=ik2ii(ik,1)
  !          iy=ik2ii(ik,2)
  !          i=(ii-1)*Nx+ix
  !          j=(jj-1)*Nx+iy
  !          Hhf_reshape(:,:,i,j) = Hhf_grid(:,:,ik)     
  !       end do
  !    end do
  ! end do

  ! do ik=1,Lk
  !    ix=ik2ii(ik,1)
  !    iy=ik2ii(ik,2)
  !    do ii=1,3
  !       do jj=1,3
  !          i=(ii-1)*3+ix
  !          j=(jj-1)*3+iy
  !          Hhf_reshape(:,:,i,j) = Hhf_grid(:,:,ik)     
  !       end do
  !    end do
  ! end do
  ! stop
  ! do ii=1,3*Nx
  !    do jj=1,3*Nx
  !       write(300,'(10F18.10)') kxx(ii),kyy(jj),Hhf_reshape(1,1,ii,jj)
  !    end do
  ! end do

  ! !+- TEST inperpolate
  ! do ik=1,Lk
  !    write(304,'(10F18.10)') k_bz(ik,:),dreal(Hhf_grid(1,1,ik)),dreal(Hhf_grid(2,2,ik)),dreal(Hhf_grid(1,2,ik))
  ! end do
 
  ! !+- interpolation -+!
  ! !allocate(Hkint(Nso,Nso))
  ! dk=pi/dble(50)
  ! kint=-dk
  ! do ix=1,50
  !    kint=kint+dk
  !    ktest_int(1) = kint
  !    ktest_int(2) = kint
  !    Hkint=Hk_HF(ktest_int,Nso)
  !    write(305,'(10F18.10)') ktest_int(:),dreal(Hkint(1,1)),dreal(Hkint(2,2)),dreal(Hkint(1,2))
  ! end do

  ! do ix=1,Nx
  !    ik=igr2ik(ix,ix)
  !    write(306,'(10F18.10)') k_bz(ik,:),dreal(Hhf_grid(1,1,ik)) 
  ! end do

  

  !+- here I should consider a generic "solve along a given BZ-path a given k-dep Hamiltonian"
  !call print_HF_bands(delta_hf)

  
  call save_array('delta.out',delta_hf)

  delta_hf=0.d0

  call read_array('delta.out',delta_hf)
  call save_array('delta_.out',delta_hf)
  stop


  H_hf=Hk;mu_fix=0.1d0
  !forall(ik=1:Lk) H_hf(:,:,ik) =  H_hf(:,:,ik) + Hsb
  call find_chem_pot(H_hf,delta_hf,mu_fix)
  delta_hf_=delta_hf
  unit_err=free_unit()
  open(unit=unit_err,file='err_loc_normal.err')
  unit_obs=free_unit()
  open(unit=unit_obs,file='obs_hf_loc_normal.loop')
  err_hf=1.d0  
  do ihf=1,Nhf
     write(unit_err,'(10F18.10)') dble(ihf),err_hf,mu_fix          
     call local_single_particle_observables(delta_hf,obs_loc)
     write(unit_obs,'(20F18.10)') dble(ihf),dreal(obs_loc)
     call build_HF_hamiltonian(H_hf,delta_hf,Uloc)          
     call find_chem_pot(H_hf,delta_hf,mu_fix)
     delta_hf=wmix*delta_hf+(1.d0-wmix)*delta_hf_
     !
     err_hf=check_conv(delta_hf,delta_hf_)     
     !
     delta_hf_=delta_hf     
  end do
  close(unit_err)
  close(unit_obs)



  H_hf=Hk;mu_fix=0.1d0
  !forall(ik=1:Lk) H_hf(:,:,ik) =  H_hf(:,:,ik) + Hsb
  call find_chem_pot(H_hf,delta_hf,mu_fix)
  delta_hf_=delta_hf
  unit_err=free_unit()
  open(unit=unit_err,file='err_q_normal.err')
  unit_obs=free_unit()
  open(unit=unit_obs,file='obs_hf_q_normal.loop')
  err_hf=1.d0  
  do ihf=1,Nhf
     write(unit_err,'(10F18.10)') dble(ihf),err_hf,mu_fix          
     call local_single_particle_observables(delta_hf,obs_loc)
     write(unit_obs,'(20F18.10)') dble(ihf),dreal(obs_loc)
     call build_HF_hamiltonian(H_hf,delta_hf,Uq_jellium)     
     
     call find_chem_pot(H_hf,delta_hf,mu_fix)
     delta_hf=wmix*delta_hf+(1.d0-wmix)*delta_hf_
     !
     err_hf=check_conv(delta_hf,delta_hf_)     
 
     !
     delta_hf_=delta_hf
     
  end do
  close(unit_err)
  close(unit_obs)


contains
  


  function Hk_HF(kpoint,N)
    implicit none
    real(8),dimension(:) :: kpoint
    integer              :: N
    complex(8),dimension(N,N) :: Hk_HF    
    integer :: io,jo
    real(8),dimension(3*Nx,3*Nx) :: Hkr
    real(8) :: Hr,dH

    if(size(kpoint).ne.2) stop "kpoint.ne.2"    
    Hk_HF=0.d0
    do io=1,Nso
       do jo=io,Nso
          Hkr=dreal(Hhf_reshape(io,jo,:,:))
          Hr=0.d0
          call polin2(kxx,kyy,Hkr,kpoint(1),kpoint(2),Hr,dH)
          Hk_Hf(io,jo) = Hk_Hf(io,jo) + Hr
          !
          Hkr=dimag(Hhf_reshape(io,jo,:,:))
          Hr=0.d0
          call polin2(kxx,kyy,Hkr,kpoint(1),kpoint(2),Hr,dH)
          Hk_Hf(io,jo) = Hk_Hf(io,jo) + xi*Hr
          Hk_Hf(jo,io) = conjg(Hk_Hf(io,jo))
          !
       end do
    end do
    
    write(*,*) Hk_Hf(1,1),kpoint,dH

  end function Hk_HF



  function Uq_jellium(q)
    USE DMFT_VECTORS  
    USE VARS_GLOBAL
    implicit none
    real(8),dimension(:) :: q
    real(8),dimension(Nso,Nso) :: Uq_jellium
    real(8) :: modq
    integer :: i
    !
    if(size(q).ne.2) stop "size q/=2; jellium" 
    modq=0.d0
    do i=1,size(q)
       modq=modq+q(i)**2.d0
    end do
    Uq_jellium = 0.d0
    if(modq.gt.1.d-10) then
       Uq_jellium = Umat/sqrt(modq)
    end if
    !
  end function Uq_jellium

  function Uloc(q)
    USE DMFT_VECTORS  
    USE VARS_GLOBAL
    implicit none
    real(8),dimension(:) :: q
    !type(vect2D) :: q
    real(8),dimension(Nso,Nso) :: Uloc
    !
    if(size(q).ne.2) stop "size q/=2; Uloc" 
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


  subroutine build_2d_bz(Nx)
    implicit none
    integer :: i,Nx,ix,iy,ik,unitk
    integer :: Nx_
    real(8) :: kmax
    real(8),dimension(:),allocatable :: kxr,kyr
    !
    !    
    allocate(kxgrid(Nx),kygrid(Nx))
    kmax=dble(Nx-1)*0.5d0/dble(Nx)
    kxgrid = linspace(-kmax,kmax,Nx,mesh=dk_mesh)
    kygrid = kxgrid    
    if(mod(Nx,2).eq.0) then
       kxgrid=kxgrid-0.5d0/dble(Nx)
       kygrid=kygrid-0.5d0/dble(Nx)
    end if
    !
    Lk=Nx*Nx
    allocate(k_bz(Lk,2),wtk(Lk),ik2ii(Lk,2))
    allocate(igr2ik(Nx,Nx))
    wtk=1.d0/dble(Lk)
    !
    KK(1)=2.d0*pi
    KK(2)=2.d0*pi
    dk_mesh=dk_mesh*2.d0*pi
    !
    unitk=free_unit()
    open(unit=unitk,file='k_points.out')
    ik=0
    do ix=1,Nx
       do iy=1,Nx
          ik = ik + 1          
          ik2ii(ik,1)=ix
          ik2ii(ik,2)=iy
          igr2ik(ix,iy) = ik
          k_bz(ik,1)=kxgrid(ix)*KK(1)
          k_bz(ik,2)=kygrid(iy)*KK(2)
          write(unitk,'(10F18.10)') k_bz(ik,1),k_bz(ik,2)
       end do
    end do
    close(unitk)
    kxgrid=kxgrid*KK(1)
    kygrid=kxgrid*KK(2)
    !
    !+- continue from here -+!
    Nx_=Nx*N_cut_off
    !
    allocate(kxr(Nx_),kyr(Nx_))
    kmax=dble(Nx_-1)*0.5d0/dble(Nx)
    kxr = linspace(-kmax,kmax,Nx_)
    kyr = kxr    
    if(mod(Nx,2).eq.0) then
       kxr=kxr-0.5d0/dble(Nx)
       kyr=kyr-0.5d0/dble(Nx)
    end if
    !
    Lkr=Nx_*Nx_
    allocate(krl(Lkr,2),wtk_rl(Lkr),ikrl2ii(Lkr,2))
    wtk_rl=1.d0/dble(Lkr)
    !
    KK(1)=2.d0*pi
    KK(2)=2.d0*pi
    !
    unitk=free_unit()
    open(unit=unitk,file='kRL_points.out')
    ik=0
    do ix=1,Nx_
       do iy=1,Nx_
          ik = ik + 1          
          krl(ik,1)=kxr(ix)*KK(1)
          krl(ik,2)=kyr(iy)*KK(2)
          ikrl2ii(ik,1)=ix
          ikrl2ii(ik,2)=iy
          write(unitk,'(10F18.10)') krl(ik,1),krl(ik,2)
       end do
    end do
    close(unitk)
    kxr=kxr*KK(1)
    kyr=kyr*KK(2)
    !
    
    ik=13
    write(*,*) krl(ik,:)
    call shift_BZ(ik,ix)
    write(*,*) k_bz(ix,:)


    ik=18*16
    write(*,*) krl(ik,:)
    call shift_BZ(ik,ix)
    write(*,*) k_bz(ix,:)

    !stop
  end subroutine build_2d_Bz






end program



!AMOEBA TEST


