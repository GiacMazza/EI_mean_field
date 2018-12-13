program officina
  USE SCIFOR
  USE DMFT_TOOLS
  !
  USE VARS_GLOBAL
  USE DMFT_VECTORS
  USE HF
  !
  USE MPI
  !
  implicit none
  integer :: Nx
  type(vect2D) :: KK

  real(8) :: KKx,KKy
  complex(8),dimension(:,:,:),allocatable :: delta_hf,H_hf,delta_hf_
  real(8) :: mu_fix
  real(8) :: Uintra,Uinter,Delta_CF
  real(8) :: Uq,thop,err_hf
  complex(8),dimension(:,:),allocatable :: Hsb

  complex(8),dimension(:),allocatable :: Uft

  integer :: iso,jso,ispin,jspin,iorb,jorb,ik
  integer :: Nhf,ihf,unit_err,unit_obs,unit_in

  complex(8),dimension(:),allocatable :: obs_loc
  real(8) :: wmix

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
  call build_2d_bz(Nx)
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
        Hk(iso,iso,ik) = -Delta_CF*0.5d0+2.d0*thop*(cos(k_bz(ik)%x)+cos(k_bz(ik)%y))
     end do
     iorb=2
     do ispin=1,Nspin
        iso=(ispin-1)*Nspin+iorb
        Hk(iso,iso,ik) = Delta_CF*0.5d0-2.d0*thop*(cos(k_bz(ik)%x)+cos(k_bz(ik)%y))
     end do
     write(unit_in,'(10F18.10)') k_bz(ik)%x,k_bz(ik)%y,dreal(Hk(1,1,ik)),dreal(Hk(2,2,ik))
     !     
  end do
  close(unit_in)

  !+- loop hartree-fock -+!
  allocate(delta_hf(Nso,Nso,Lk),H_hf(Nso,Nso,Lk),delta_hf_(Nso,Nso,Lk))
  !
  allocate(Hsb(Nso,Nso)); Hsb(1,2)=1.d-2; Hsb(2,1)=1.d-2
  H_hf=Hk;mu_fix=0.1d0
  forall(ik=1:Lk) H_hf(:,:,ik) =  H_hf(:,:,ik) + Hsb
  call find_chem_pot(H_hf,delta_hf,mu_fix)
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
  !
  !
  !
  H_hf=Hk;mu_fix=0.1d0
  forall(ik=1:Lk) H_hf(:,:,ik) =  H_hf(:,:,ik) + Hsb
  call find_chem_pot(H_hf,delta_hf,mu_fix)
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


  subroutine build_2d_bz(Nx)
    implicit none
    real(8),dimension(:),allocatable :: kxgrid,kygrid
    integer :: i,Nx,ix,iy,ik,unitk
    real(8) :: kmax
    !
    !    
    allocate(kxgrid(Nx),kygrid(Nx))
    kmax=dble(Nx-1)*0.5d0/dble(Nx)
    kxgrid = linspace(-kmax,kmax,Nx)
    kygrid = kxgrid    
    if(mod(Nx,2).eq.0) then
       kxgrid=kxgrid-0.5d0/dble(Nx)
       kygrid=kygrid-0.5d0/dble(Nx)
    end if
    !
    Lk=Nx*Nx
    allocate(k_bz(Lk),wtk(Lk))
    wtk=1.d0/dble(Lk)
    !
    KK%x=2.d0*pi
    KK%y=2.d0*pi
    !
    unitk=free_unit()
    open(unit=unitk,file='k_points.out')
    ik=0
    do ix=1,Nx
       do iy=1,Nx
          ik = ik + 1
          k_bz(ik)%x=kxgrid(ix)*KK%x
          k_bz(ik)%y=kygrid(iy)*KK%y
          write(unitk,'(10F18.10)') k_bz(ik)%x,k_bz(ik)%y
       end do
    end do
    close(unitk)
    !
  end subroutine build_2d_Bz



  ! function exciton_ins_photon(x) result(out_x)
  !   implicit none
  !   real(8),dimension(:)   :: x
  !   real(8),dimension(size(x)) :: out_x
  !   real(8) :: Iexc,dpol,phi0,EK,wp,wm,V_int_,df
  !   real(8) :: dpol_,Iexc_,xmu,phi0_
  !   real(8) :: e1,e2,ep,em
  !   integer :: ik
  !   real(8) :: g_,gc,w_gap_,w_ph_,Delta,w0,gk
  !   !
  !   !+-> bogoliubov trans for photons <-+!
  !   w_gap_=w_gap
  !   w_ph_=w_ph
  !   g_=g
  !   !
  !   Delta=Asq*g_*g_/w_gap_
  !   w0=w_ph_+2.d0*Delta
  !   w_ph_ = sqrt(w_ph_*w_ph_+4.d0*Delta*w_ph_)
  !   gc=g_*(sqrt(w0+w_ph_)-sqrt(w0-w_ph_))/sqrt(2.d0*w_ph_)
  !   !
  !   w_ph_= w_ph_-4.d0*dw0 ! coming from the bath
  !   !
  !   dpol=x(1)
  !   Iexc=x(2)
  !   phi0=x(3)
  !   !
  !   !
  !   w_gap_ = w_gap
  !   V_int_ = V_int
  !   xmu = -V_int_*0.5d0
  !   !
  !   dpol_ = 0.d0
  !   Iexc_ = 0.d0
  !   phi0_ = 0.d0
  !   do ik=1,Lk
  !      !+- bloch electrons: Dipolar couplings -+!
  !      e1 = alpha_hop*epsik(ik)-w_gap_*0.5d0
  !      e2 = epsik(ik)+w_gap_*0.5d0       
  !      gk=gc
  !      if(gck) gk=gc*sqrt(abs(e2-e1)/w_gap_)
  !      !+--------------------------------------+!

  !      !+- get back to the 2x2 Hamiltonian -+!
  !      e1 = e1 + xmu + V_int_*0.5d0*(1.d0-dpol)
  !      e2 = e2 + xmu + V_int_*0.5d0*(1.d0+dpol)       
  !      vk = -V_int_*Iexc + 2.d0*phi0*gk 
  !      !
  !      wp = 0.5d0*(e1+e2)
  !      wm = 0.5d0*(e1-e2)
  !      Ek = sqrt(wm**2.d0+vk**2.d0)
  !      !
  !      df=(fermi(wp+Ek,beta)-fermi(wp-Ek,beta))
  !      !
  !      dpol_ = dpol_ + df*wm/Ek*wtk(ik)
  !      Iexc_ = Iexc_ + 0.5d0*df*vk/Ek*wtk(ik)
  !      phi0_ = phi0_ - df*gk*vk/w_ph_/Ek*wtk(ik)
  !      !
  !   end do
  !   !
  !   out_x(1) = dpol_
  !   out_x(2) = Iexc_
  !   out_x(3) = phi0_
  !   !
  ! end function exciton_ins_photon




  ! function delta_exciton_ins_photon(x) result(out_x)
  !   implicit none
  !   real(8),dimension(:)   :: x
  !   real(8),dimension(size(x)) :: out_x
  !   real(8) :: Iexc,dpol,phi0,EK,wp,wm,V_int_,df
  !   real(8) :: dpol_,Iexc_,xmu,phi0_
  !   real(8) :: e1,e2
  !   integer :: ik
  !   real(8) :: g_,gc,w_gap_,w_ph_,Delta,w0,gk
  !   !
  !   !+-> bogoliubov trans for photons <-+!
  !   w_gap_=w_gap
  !   w_ph_=w_ph
  !   g_=g
  !   !
  !   Delta=Asq*g_*g_
  !   w0=w_ph_+2.d0*Delta
  !   w_ph_ = sqrt(w_ph_*w_ph_+4.d0*Delta*w_ph_)
  !   gc=g_*(sqrt(w0+w_ph_)-sqrt(w0-w_ph_))/sqrt(2.d0*w_ph_)
  !   !
  !   w_ph_= w_ph_-4.d0*dw0 ! coming from the bath
  !   !
  !   dpol=x(1)
  !   Iexc=x(2)
  !   phi0=x(3)
  !   !
  !   !
  !   w_gap_ = w_gap
  !   V_int_ = V_int
  !   xmu = -V_int_*0.5d0
  !   !
  !   dpol_ = 0.d0
  !   Iexc_ = 0.d0
  !   phi0_ = 0.d0
  !   do ik=1,Lk
  !      !+- bloch electrons: Dipolar couplings -+!
  !      e1 = alpha_hop*epsik(ik)-w_gap_*0.5d0
  !      e2 = epsik(ik)+w_gap_*0.5d0
  !      gk=gc
  !      if(gck) gk=gc*sqrt(abs(e2-e1))
  !      !+--------------------------------------+!

  !      !+- 2x2 Hamiltonian -+!
  !      e1 = e1 + xmu + V_int_*0.5d0*(1.d0-dpol) 
  !      e2 = e2 + xmu + V_int_*0.5d0*(1.d0+dpol) 
  !      vk = -V_int_*Iexc + 2.d0*phi0*gk
  !      !
  !      wp = 0.5d0*(e1+e2)
  !      wm = 0.5d0*(e1-e2)
  !      Ek = sqrt(wm**2.d0+vk**2.d0)
  !      !
  !      df=(fermi(wp+Ek,beta)-fermi(wp-Ek,beta))
  !      !
  !      dpol_ = dpol_ + df*wm/Ek*wtk(ik)
  !      Iexc_ = Iexc_ + 0.5d0*df*vk/Ek*wtk(ik)
  !      phi0_ = phi0_ - df*gk*vk/w_ph_/Ek*wtk(ik)
  !   end do
  !   !
  !   out_x(1) = dpol_-dpol
  !   out_x(2) = Iexc_-Iexc
  !   out_x(3) = phi0_-phi0
  !   !
  ! end function delta_exciton_ins_photon
  ! !
  ! !

  ! !
  ! function critical_coupling_met_ins(x) result(out_x)
  !   implicit none
  !   real(8),intent(in)   :: x
  !   real(8) :: out_x,beta_
  !   real(8) :: Iexc,dpol,phi0,EK,wp,wm,V_int_,df
  !   real(8) :: dpol_,Iexc_,xmu,phi0_
  !   real(8) :: e1,e2
  !   integer :: ik,iloop
  !   real(8) :: g_,gc,w_gap_,w_ph_,Delta,w0,gk
  !   real(8) :: M_gap(2,2)
  !   !
  !   beta_=beta
  !   w_gap_=w_gap
  !   w_ph_=w_ph
  !   g_=0.d0
  !   !
  !   Delta=Asq*g_*g_
  !   w0=w_ph_+2.d0*Delta
  !   w_ph_ = sqrt(w_ph_*w_ph_+4.d0*Delta*w_ph_)
  !   gc=g_*(sqrt(w0+w_ph_)-sqrt(w0-w_ph_))/sqrt(2.d0*w_ph_)
  !   !
  !   w_ph_= w_ph_-4.d0*dw0 ! coming from the bath
  !   !
  !   w_gap_ = w_gap
  !   V_int_ = x
  !   xmu = -V_int_*0.5d0
  !   !+- compute the polarization -+!
  !   dpol=0.d0
  !   do iloop=1,100       
  !      !
  !      dpol_=0.d0
  !      do ik=1,Lk
  !         e1 = alpha_hop*epsik(ik)-w_gap_*0.5d0
  !         e2 = epsik(ik)+w_gap_*0.5d0
  !         !
  !         e1 = e1 + xmu + V_int_*0.5d0*(1.d0-dpol)
  !         e2 = e2 + xmu + V_int_*0.5d0*(1.d0+dpol)       
  !         dpol_=dpol_+(fermi(e1,beta_)-fermi(e2,beta_))*wtk(ik)
  !         !
  !      end do
  !      if(abs(dpol-dpol_).lt.1.d-8) exit
  !      dpol=dpol_
  !      !
  !   end do
  !   !
  !   out_x=V_int_*dpol-(4.d0*thop-w_gap_)
  !   !
  ! end function critical_coupling_met_ins



  ! function critical_coupling_exciton_ins_photon(x) result(out_x)
  !   implicit none
  !   real(8),intent(in)   :: x
  !   real(8) :: out_x,beta_
  !   real(8) :: Iexc,dpol,phi0,EK,wp,wm,V_int_,df
  !   real(8) :: dpol_,Iexc_,xmu,phi0_
  !   real(8) :: e1,e2
  !   integer :: ik,iloop
  !   real(8) :: g_,gc,w_gap_,w_ph_,Delta,w0,gk
  !   real(8) :: M_gap(2,2)
  !   !
  !   beta_=beta
  !   w_gap_=w_gap
  !   w_ph_=w_ph
  !   g_=x
  !   !
  !   Delta=Asq*g_*g_
  !   w0=w_ph_+2.d0*Delta
  !   w_ph_ = sqrt(w_ph_*w_ph_+4.d0*Delta*w_ph_)
  !   gc=g_*(sqrt(w0+w_ph_)-sqrt(w0-w_ph_))/sqrt(2.d0*w_ph_)
  !   !
  !   w_ph_= w_ph_-4.d0*dw0 ! coming from the bath
  !   !
  !   w_gap_ = w_gap
  !   V_int_ = V_int
  !   xmu = -V_int_*0.5d0
  !   !+- compute the polarization -+!
  !   dpol=0.d0
  !   do iloop=1,100       
  !      !
  !      dpol_=0.d0
  !      do ik=1,Lk
  !         e1 = alpha_hop*epsik(ik)-w_gap_*0.5d0
  !         e2 = epsik(ik)+w_gap_*0.5d0
  !         !
  !         e1 = e1 + xmu + V_int_*0.5d0*(1.d0-dpol)
  !         e2 = e2 + xmu + V_int_*0.5d0*(1.d0+dpol)       
  !         dpol_=dpol_+(fermi(e1,beta_)-fermi(e2,beta_))*wtk(ik)
  !         !
  !      end do
  !      if(abs(dpol-dpol_).lt.1.d-12) exit
  !      dpol=dpol_
  !      !
  !   end do
  !   !
  !   Iexc=0.d0
  !   phi0=0.d0
  !   !
  !   M_gap=0.d0
  !   M_gap(1,1)=-1.d0
  !   M_gap(2,2)=-1.d0
  !   do ik=1,Lk       
  !      e1 = alpha_hop*epsik(ik)-w_gap_*0.5d0
  !      e2 = epsik(ik)+w_gap_*0.5d0
  !      gk=gc
  !      if(gck) gk=gc*sqrt(abs(e2-e1))
  !      !+- 2x2 Hamiltonian -+!
  !      e1 = e1 + xmu + V_int_*0.5d0*(1.d0-dpol)
  !      e2 = e2 + xmu + V_int_*0.5d0*(1.d0+dpol)       
  !      !+------------------+!
  !      vk = -V_int_*Iexc + 2.d0*phi0*gk
  !      !
  !      wp = 0.5d0*(e1+e2)
  !      wm = 0.5d0*(e1-e2)
  !      Ek = sqrt(wm**2.d0+vk**2.d0)
  !      !
  !      df=(fermi(wp+Ek,beta_)-fermi(wp-Ek,beta_))
  !      !
  !      M_gap(1,1) = M_gap(1,1) - 2.d0*gk*gk/w_ph_/Ek*df*wtk(ik)
  !      M_gap(1,2) = M_gap(1,2) + V_int_*gk/w_ph_/Ek*df*wtk(ik)       
  !      M_gap(2,1) = M_gap(2,1) + 2.d0*gk*0.5d0/Ek*df*wtk(ik)
  !      M_gap(2,2) = M_gap(2,2) - V_int_*0.5d0/Ek*df*wtk(ik)
  !      !
  !      !
  !   end do
  !   !
  !   out_x=M_gap(1,1)*M_gap(2,2)-M_gap(1,2)*M_gap(2,1)
  !   !if(mpiID==0) write(150,*) x,M_gap(1,1),M_gap(2,2),M_gap(1,2),M_gap(2,1)
  !   !
  ! end function critical_coupling_exciton_ins_photon
  ! !
  ! !
  ! function critical_temperature_exciton_ins_photon(x) result(out_x)
  !   implicit none
  !   real(8),intent(in)   :: x
  !   real(8) :: out_x,beta_
  !   real(8) :: Iexc,dpol,phi0,EK,wp,wm,V_int_,df
  !   real(8) :: dpol_,Iexc_,xmu,phi0_
  !   real(8) :: e1,e2
  !   integer :: ik,iloop
  !   real(8) :: g_,gc,w_gap_,w_ph_,Delta,w0,gk
  !   real(8) :: M_gap(2,2)
  !   !
  !   beta_=1.d0/x
  !   w_gap_=w_gap
  !   w_ph_=w_ph
  !   g_=g
  !   !
  !   Delta=Asq*g_*g_
  !   w0=w_ph_+2.d0*Delta
  !   w_ph_ = sqrt(w_ph_*w_ph_+4.d0*Delta*w_ph_)
  !   gc=g_*(sqrt(w0+w_ph_)-sqrt(w0-w_ph_))/sqrt(2.d0*w_ph_)
  !   !
  !   w_ph_= w_ph_-4.d0*dw0 ! coming from the bath
  !   !
  !   w_gap_ = w_gap
  !   V_int_ = V_int
  !   xmu = -V_int_*0.5d0
  !   !+- compute the polarization -+!
  !   dpol=0.d0
  !   do iloop=1,100       
  !      !
  !      dpol_=0.d0
  !      do ik=1,Lk
  !         e1 = alpha_hop*epsik(ik)-w_gap_*0.5d0
  !         e2 = epsik(ik)+w_gap_*0.5d0
  !         !
  !         e1 = e1 + xmu + V_int_*0.5d0*(1.d0-dpol)
  !         e2 = e2 + xmu + V_int_*0.5d0*(1.d0+dpol)       
  !         dpol_=dpol_+(fermi(e1,beta_)-fermi(e2,beta_))*wtk(ik)
  !         !
  !      end do
  !      if(abs(dpol-dpol_).lt.1.d-8) exit
  !      dpol=dpol_
  !      !
  !   end do
  !   !
  !   Iexc=0.d0
  !   phi0=0.d0
  !   !
  !   M_gap=0.d0
  !   M_gap(1,1)=-1.d0
  !   M_gap(2,2)=-1.d0
  !   do ik=1,Lk       
  !      e1 = alpha_hop*epsik(ik)-w_gap_*0.5d0
  !      e2 = epsik(ik)+w_gap_*0.5d0
  !      gk=gc
  !      if(gck) gk=gc*sqrt(abs(e2-e1))
  !      !+- 2x2 Hamiltonian -+!
  !      e1 = e1 + xmu + V_int_*0.5d0*(1.d0-dpol)
  !      e2 = e2 + xmu + V_int_*0.5d0*(1.d0+dpol)       
  !      !+------------------+!
  !      vk = -V_int_*Iexc + 2.d0*phi0*gk
  !      !
  !      wp = 0.5d0*(e1+e2)
  !      wm = 0.5d0*(e1-e2)
  !      Ek = sqrt(wm**2.d0+vk**2.d0)
  !      !
  !      df=(fermi(wp+Ek,beta_)-fermi(wp-Ek,beta_))
  !      !
  !      M_gap(1,1) = M_gap(1,1) - 2.d0*gk*gk/w_ph_/Ek*df*wtk(ik)
  !      M_gap(1,2) = M_gap(1,2) + V_int_*gk/w_ph_/Ek*df*wtk(ik)       
  !      M_gap(2,1) = M_gap(2,1) + 2.d0*gk*0.5d0/Ek*df*wtk(ik)
  !      M_gap(2,2) = M_gap(2,2) - V_int_*0.5d0/Ek*df*wtk(ik)
  !      !
  !      !
  !   end do
  !   !
  !   out_x=M_gap(1,1)*M_gap(2,2)-M_gap(1,2)*M_gap(2,1)
  !   !
  ! end function critical_temperature_exciton_ins_photon
  ! !
  ! !
  ! function bound_states_ei(x) result(out_x)
  !   implicit none
  !   real(8),intent(in)   :: x
  !   real(8) :: out_x,beta_
  !   real(8) :: Iexc,dpol,phi0,EK,wp,wm,V_int_,df
  !   real(8) :: dpol_,Iexc_,xmu,phi0_
  !   real(8) :: e1,e2
  !   integer :: ik,iloop
  !   real(8) :: g_,gc,w_gap_,w_ph_,Delta,w0,gk
  !   real(8) :: M_gap(2,2)
  !   real(8) :: omega_chi
  !   complex(8) ::chi0_tmp
  !   !
  !   beta_=beta
  !   w_gap_=w_gap
  !   w_ph_=w_ph
  !   g_=0.d0
  !   !
  !   Delta=Asq*g_*g_
  !   w0=w_ph_+2.d0*Delta
  !   w_ph_ = sqrt(w_ph_*w_ph_+4.d0*Delta*w_ph_)
  !   gc=g_*(sqrt(w0+w_ph_)-sqrt(w0-w_ph_))/sqrt(2.d0*w_ph_)
  !   !
  !   w_ph_= w_ph_-4.d0*dw0 ! coming from the bath
  !   !
  !   w_gap_ = w_gap
  !   V_int_ = V_int
  !   xmu = -V_int_*0.5d0
  !   !+- compute the polarization -+!
  !   dpol=0.d0
  !   do iloop=1,100       
  !      !
  !      dpol_=0.d0
  !      do ik=1,Lk
  !         e1 = alpha_hop*epsik(ik)-w_gap_*0.5d0
  !         e2 = epsik(ik)+w_gap_*0.5d0
  !         !
  !         e1 = e1 + xmu + V_int_*0.5d0*(1.d0-dpol)
  !         e2 = e2 + xmu + V_int_*0.5d0*(1.d0+dpol)       
  !         dpol_=dpol_+(fermi(e1,beta_)-fermi(e2,beta_))*wtk(ik)
  !         !
  !      end do
  !      if(abs(dpol-dpol_).lt.1.d-8) exit
  !      dpol=dpol_
  !      !
  !   end do
  !   !
  !   chi0_tmp=0.d0
  !   omega_chi=x
  !   do ik=1,Lk
  !      e1 = alpha_hop*epsik(ik)-w_gap_*0.5d0
  !      e2 = epsik(ik)+w_gap_*0.5d0
  !      !
  !      e1 = e1 + xmu + V_int_*0.5d0*(1.d0-dpol)
  !      e2 = e2 + xmu + V_int_*0.5d0*(1.d0+dpol)       
  !      !
  !      df= fermi(e2,beta) - fermi(e1,beta)
  !      chi0_tmp = chi0_tmp + df/(omega_chi + (e2-e1)+xi*eps)*wtk(ik)     
  !   end do
  !   !
  !   write(*,*) x,dreal(chi0_tmp),1.d0/dreal(chi0_tmp)+V_int_
  !   out_x= 1.d0/dreal(chi0_tmp)+V_int_
  !   !
  ! end function bound_states_ei



  ! !
  ! function tot_energy(x) result(x_out)
  !   implicit none
  !   real(8),dimension(:)   :: x
  !   real(8),dimension(size(x)) :: x_out
  !   real(8) :: Iexc,dpol,phi0,wp,wm,V_int_,df
  !   real(8) :: dpol_,Iexc_,xmu,phi0_
  !   real(8) :: e1,e2
  !   real(8),dimension(2,2) :: Hk
  !   real(8),dimension(2) :: Ek,tk
  !   integer :: ik,iorb,jorb,korb
  !   real(8) :: g_,gc,w_gap_,w_ph_,Delta,w0,gk
  !   real(8) :: E_tot,E_kin,E_pot
  !   !
  !   !+-> bogoliubov trans for photons <-+!
  !   w_gap_=w_gap
  !   w_ph_=w_ph
  !   g_=g
  !   !
  !   Delta=Asq*g_*g_
  !   w0=w_ph_+2.d0*Delta
  !   w_ph_ = sqrt(w_ph_*w_ph_+4.d0*Delta*w_ph_)
  !   gc=g_*(sqrt(w0+w_ph_)-sqrt(w0-w_ph_))/sqrt(2.d0*w_ph_)
  !   !
  !   w_ph_= w_ph_-4.d0*dw0 ! coming from the bath
  !   !
  !   dpol=x(1)
  !   Iexc=x(2)
  !   phi0=x(3)
  !   !
  !   !
  !   w_gap_ = w_gap
  !   V_int_ = V_int
  !   xmu = -V_int_*0.5d0
  !   !
  !   E_kin=0.d0
  !   E_tot=0.d0
  !   do ik=1,Lk
  !      !+- bloch electrons: Dipolar couplings -+!
  !      e1 = alpha_hop*epsik(ik)-w_gap_*0.5d0
  !      e2 = epsik(ik)+w_gap_*0.5d0
  !      !
  !      tk(1) = alpha_hop*epsik(ik)
  !      tk(2) = epsik(ik)
  !      !
  !      gk=gc
  !      if(gck) gk=gc*sqrt(abs(e2-e1))
  !      !+--------------------------------------+!

  !      !+- 2x2 Hamiltonian -+!
  !      e1 = e1 + xmu + V_int_*0.5d0*(1.d0-dpol)
  !      e2 = e2 + xmu + V_int_*0.5d0*(1.d0+dpol)       
  !      vk = -V_int_*Iexc + 2.d0*phi0*gk
  !      !
  !      Hk(1,1) = e1
  !      Hk(2,2) = e2
  !      Hk(1,2) = vk
  !      Hk(2,1) = vk
  !      !
  !      call matrix_diagonalize(Hk,Ek)
  !      !
  !      E_tot=E_tot+Ek(1)*fermi(Ek(1),beta)*wtk(ik)
  !      E_tot=E_tot+Ek(2)*fermi(Ek(2),beta)*wtk(ik)       
  !      !
  !      do iorb=1,2
  !         do jorb=1,2
  !            E_kin = E_kin + tk(iorb)*Hk(iorb,jorb)*Hk(iorb,jorb)*fermi(Ek(jorb),beta)*wtk(ik)
  !            !write(455,*) Ekin,tk(iorb)*Hk(iorb,jorb)*Hk(iorb,jorb)*fermi(Ek(jorb),beta)*wtk(ik)
  !         end do
  !      end do
  !      !
  !   end do
  !   ! write(455,*)
  !   ! write(455,*) Ekin
  !   !
  !   !+- double counting terms -+!
  !   E_tot=E_tot-V_int_*(1.d0-dpol*dpol)*0.25d0
  !   E_tot=E_tot+V_int*Iexc*Iexc
  !   E_pot = E_tot-E_kin
  !   !
  !   x_out(1) = E_tot
  !   x_out(2) = E_kin
  !   x_out(3) = E_pot
  !   !write(*,*) x_out
  !   !
  ! end function tot_energy
  ! !
  ! !
  ! !+- old self consistency -+!
  ! !
  ! !
  ! function self_consistency(hyb) result(f)
  !   real(8),intent(in) :: hyb
  !   real(8)  :: f
  !   integer :: ik

  !   real(8) :: e1,e2,ep,em,ek,hyb_


  !   f=0.d0
  !   do ik=1,Lk
  !      !
  !      e1 = alpha_hop*epsik(ik)-w_gap*0.5d0 -xmu
  !      e2 = epsik(ik)+w_gap*0.5d0 - xmu
  !      ep=0.5d0*(e1+e2)
  !      em=0.5d0*(e1-e2)
  !      hyb_=hyb*g*g/w_ph
  !      !
  !      ek=sqrt(em**2 + hyb_**2)       
  !      !
  !      f= f + 0.5d0*g*g/w_ph * (fermi(ep+ek,beta)-fermi(ep-ek,beta))/ek*wtk(ik)
  !      !
  !   end do
  !   f=f+1.d0
  !   !write(405,*) hyb,f 
  ! end function self_consistency

  ! function self_consistency_CR(hyb) result(f)
  !   real(8),intent(in) :: hyb
  !   real(8)  :: f
  !   integer :: ik

  !   real(8) :: e1,e2,ep,em,ek,hyb_
  !   real(8) :: g_,gc,w_gap_,w_ph_,Delta,w0

  !   !+-> bogoliubov trans
  !   w_gap_=w_gap
  !   w_ph_=w_ph
  !   g_=g
  !   !
  !   Delta=Asq*g_*g_/w_gap_
  !   w0=w_ph_+2.d0*Delta
  !   w_ph_ = sqrt(w_ph_*w_ph_+4.d0*Delta*w_ph_)
  !   gc=g_*(sqrt(w0+w_ph_)-sqrt(w0-w_ph_))/sqrt(2.d0*w_ph_)
  !   !
  !   w_ph_= w_ph_-4.d0*dw0


  !   f=0.d0
  !   do ik=1,Lk
  !      !
  !      e1 = alpha_hop*epsik(ik)-w_gap_*0.5d0 -xmu
  !      e2 = epsik(ik)+w_gap_*0.5d0 - xmu
  !      ep=0.5d0*(e1+e2)
  !      em=0.5d0*(e1-e2)
  !      hyb_=2.d0*hyb*gc*gc/w_ph_
  !      !
  !      ek=sqrt(em**2 + hyb_**2)       
  !      !
  !      f= f + 2.d0*gc*gc/w_ph_ * (fermi(ep+ek,beta)-fermi(ep-ek,beta))/ek*wtk(ik)
  !      !
  !   end do
  !   f=f+1.d0
  !   !write(405,*) hyb,f 
  ! end function self_consistency_CR





  ! function phase_boundary(g_) result(f)
  !   real(8),intent(in) :: g_
  !   real(8)  :: f
  !   integer :: ik

  !   real(8) :: e1,e2,ep,em,ek,hyb_

  !   f=0.d0
  !   do ik=1,Lk
  !      !
  !      e1 = alpha_hop*epsik(ik)-w_gap*0.5d0 -xmu
  !      e2 = epsik(ik)+w_gap*0.5d0 - xmu
  !      ep=0.5d0*(e1+e2)
  !      em=0.5d0*(e1-e2)
  !      hyb_=0.d0
  !      !
  !      ek=sqrt(em**2 + hyb_**2)       
  !      !
  !      f= f + 0.5d0*g_*g_/w_ph * (fermi(ep+ek,beta)-fermi(ep-ek,beta))/ek*wtk(ik)
  !      !
  !   end do
  !   f=f+1.d0
  !   !write(405,*) hyb,f 
  ! end function phase_boundary


  ! function phase_boundary_CR(g_) result(f)
  !   real(8),intent(in) :: g_
  !   real(8)  :: f
  !   integer :: ik
  !   real(8) :: e1,e2,ep,em,ek,hyb_,kcoup,gtmp
  !   real(8) :: gc,w_gap_,w_ph_,Delta,w0

  !   !+-> bogoliubov trans
  !   w_gap_=w_gap
  !   w_ph_=w_ph
  !   !
  !   Delta=Asq*g_*g_/w_gap_
  !   w0=w_ph_+2.d0*Delta
  !   w_ph_ = sqrt(w_ph_*w_ph_+4.d0*Delta*w_ph_)
  !   gc=g_*(sqrt(w0+w_ph_)-sqrt(w0-w_ph_))/sqrt(2.d0*w_ph_)
  !   !
  !   w_ph_= w_ph_-4.d0*dw0
  !   !
  !   kcoup=0.d0
  !   if(gck) kcoup=1.d0
  !   !
  !   f=0.d0
  !   do ik=1,Lk
  !      !
  !      e1 = alpha_hop*epsik(ik)-w_gap_*0.5d0 -xmu
  !      e2 = epsik(ik)+w_gap_*0.5d0 - xmu
  !      ep=0.5d0*(e1+e2)
  !      em=0.5d0*(e1-e2)
  !      hyb_=0.d0
  !      !
  !      ek=sqrt(em**2 + hyb_**2)       
  !      !
  !      gtmp=gc*sqrt(w_gap_+kcoup*epsik(ik)*(1.d0-alpha_hop))
  !      f= f + 2.d0*gtmp*gtmp/w_ph_*(fermi(ep+ek,beta)-fermi(ep-ek,beta))/ek*wtk(ik)
  !      !
  !      !write(405,*) ik,gtmp,2.d0*gtmp*gtmp/w_ph_*(fermi(ep+ek,beta)-fermi(ep-ek,beta))/ek
  !   end do
  !   !stop
  !   f=f+1.d0
  !   !
  ! end function phase_boundary_CR


  ! subroutine build_lattice_model  
  !   implicit none
  !   real(8) :: kx
  !   integer :: i,ir
  !   real(8) :: check,R
  !   complex(8) :: tmp
  !   !
  !   !
  !   allocate(epsik(Lk),wtk(Lk))
  !   allocate(vk_orb(2,Lk))

  !   check=0.d0
  !   do i=1,Lk
  !      kx = pi/dble(Lk+1)*dble(i)
  !      epsik(i) = -2.d0*thop*cos(kx)
  !      vk_orb(1,i) = 2.d0*thop*sin(kx)*alpha_hop
  !      vk_orb(2,i) = 2.d0*thop*sin(kx)       
  !      wtk(i) = 1.d0/dble(Lk)       
  !      !write(373,'(10F18.10)') kx,1.d0+epsik(i)*(1.d0-alpha_hop)
  !      check=check+1.d0/(1.d0+epsik(i)*(1.d0-alpha_hop))*wtk(i)
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
  ! end subroutine build_lattice_model



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










  !   function delta_hyb(alpha) result(delta)
  !     implicit none
  !     real(8),dimension(:) :: alpha
  !     real(8),dimension(size(alpha)) :: delta
  !     real(8) :: alpha_,nph,cond
  !     complex(8),dimension(:,:),allocatable :: Hk,Hph
  !     complex(8),dimension(:,:),allocatable :: fermi_ave
  !     complex(8),dimension(:,:,:,:),allocatable :: Gk
  !     complex(8),dimension(:,:,:),allocatable :: Gloc
  !     complex(8),dimension(2,2) :: Gktmp
  !     real(8),dimension(:),allocatable :: ek,Eph
  !     real(8) :: ph_gs,ph_gs_old
  !     integer :: ph_max
  !     allocate(Hk(2,2),ek(2),fermi_ave(2,2))    
  !     !+- BOSON part -+!
  !     ph_gs_old=1.d0
  !     do ph_max=Nmax,Nmax
  !        if(allocated(Hph)) deallocate(Hph)
  !        if(allocated(Eph)) deallocate(Eph)
  !        allocate(Hph(ph_max,ph_max),Eph(ph_max))
  !        Hph=zero
  !        do in=1,ph_max
  !           !+- build Hph
  !           Hph(in,in) = w_ph*dble(in-1)
  !           if(in.lt.ph_max) then
  !              Hph(in,in+1) = g*sqrt(dble(in))*alpha(1)!fermi_ave(2,1)
  !              Hph(in+1,in) = g*sqrt(dble(in))*alpha(1)!fermi_ave(1,2)
  !           end if
  !        end do
  !        !
  !        call matrix_diagonalize(Hph,Eph)
  !        ph_gs=Eph(1)       
  !        if(abs(ph_gs-ph_gs_old).lt.10.d-10) exit
  !        ph_gs_old=ph_gs
  !        !
  !     end do
  !     ph_max=ph_max-1

  !     ! !+- compute new_condensate
  !     cond=0.d0
  !     nph=0.d0
  !     do in=1,ph_max
  !        nph = nph + conjg(Hph(1,in))*Hph(1,in)*dble(in-1)
  !        if(in.lt.ph_max) then
  !           cond = cond + conjg(Hph(in,1))*Hph(in+1,1)*sqrt(dble(in+1))
  !        end if
  !     end do
  !     !
  !     fermi_ave = zero     
  !     Ekin = 0.d0
  !     do ik=1,Lk
  !        !
  !        Hk(1,1) = alpha_hop*epsik(ik)-w_gap*0.5d0
  !        Hk(2,2) = epsik(ik)+w_gap*0.5d0
  !        Hk(1,2) = g*cond
  !        Hk(2,1) = g*cond
  !        !
  !        call matrix_diagonalize(Hk,ek)
  !        !
  !        do is=1,2
  !           do js=1,2              
  !              do ks=1,2
  !                 fermi_ave(is,js) = fermi_ave(is,js) + &
  !                      conjg(Hk(is,ks))*Hk(js,ks)*fermi(ek(ks),1000.d0)*wtk(ik)                 
  !              end do
  !           end do
  !        end do
  !     end do

  !     !
  !     delta(1) = fermi_ave(1,2)!alpha(1)-alpha_
  ! !    iself=iself+1
  !     write(unit_err,'(10F18.10)') dble(iself),abs(delta(1)),dble(ph_max)
  !     write(unit_slater,'(10F18.10)') dble(iself),dreal(fermi_ave(1,1)),dreal(fermi_ave(2,2)),dreal(fermi_ave(1,2)),dreal(fermi_ave(2,1))
  !     write(unit_ph,'(10F18.10)') dble(iself),cond,nph
  !     write(666,*) alpha,cond
  !     !
  !   end function delta_hyb








end program Officina



!AMOEBA TEST


