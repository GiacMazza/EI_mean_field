program officina  !
  USE SCIFOR
  USE DMFT_TOOLS
  !
  USE VARS_GLOBAL
  USE HF_real
  USE MPI
  !
  implicit none
  !integer :: Nint
  real(8) :: KKx,KKy
  complex(8),dimension(:,:,:),allocatable :: delta_hf,H_hf,delta_hf_,Hhf_grid,delta_hfr,delta_hfr_,H_hf_dyn
  complex(8),dimension(:,:,:),allocatable :: Hr_w90,Hr_w90_tmp,Hr_toy  
  real(8) :: mu_fix
  real(8) :: Uintra,Uinter,Delta_CF
  real(8) :: Uq,thop,err_hf,hf_conv,ktmp
  logical :: plot_w90_bands
  integer :: iconv
  !
  real(8),dimension(3) :: ktest
  real(8),dimension(:),allocatable :: ek_out,ek_out_up,ek_out_dw
  !
  complex(8),dimension(:),allocatable :: Uft
  complex(8),dimension(:,:),allocatable :: Hloc,Hktmp
  complex(8),dimension(:,:,:),allocatable :: Hk_w90,Hk_w90_tmp,Hk_toy
  complex(8),dimension(:,:,:,:,:),allocatable :: Hk_w90_reshape,Hk_hf_reshape
  real(8),dimension(:,:),allocatable :: kpts,kpt_path
  !
  integer,dimension(:,:),allocatable :: irvec2d,itmp
  integer,dimension(:),allocatable :: stride2D,stride2D_
  integer,dimension(:,:),allocatable :: irvec1d
  integer,dimension(:),allocatable :: stride1D,stride1D_
  !
  integer :: nr2d,nr1d
  integer :: iso,jso,ispin,jspin,iorb,jorb,ik,isys,jk,iik
  integer :: i,j,k,idim
  integer :: Nhf,Nhf_,ihf,unit_err,unit_obs,unit_in,uio,unit_io,Nobs,jhf
  integer :: Nhf_opt,it_print
  integer,dimension(:),allocatable :: iorb_to_ihf
  integer,dimension(:),allocatable :: ihf_to_iorb
  integer,dimension(:,:),allocatable :: ijorb_to_ihf
  complex(8),dimension(:,:),allocatable :: ni_orb
  
  integer,dimension(:),allocatable :: units_loc_obs
  integer :: flen,iread

  integer :: ir,jr,kr,iir
  real(8),dimension(3) :: Rlat 

  complex(8),dimension(:),allocatable :: obs,obs_loc
  real(8) :: wmix
  complex(8),dimension(:,:,:),allocatable :: Ur_TNS,Vr_TNS
  complex(8),dimension(:,:),allocatable :: Uss_vs_R,UsAs_vs_R
  complex(8),dimension(:,:),allocatable :: Uss_vs_q,UsAs_vs_q
  
  real(8) :: Uloc,Vloc
  complex(8),dimension(:,:,:),allocatable :: Uq_TNS,Vq_TNS,Uq_log,Vq_log
  complex(8),dimension(:,:),allocatable :: check_Uloc

  real(8),dimension(2) :: Uq_read,Vq_read

  character(len=100) :: fileRlattice,file_w90_hr,file_UV,read_tns
  character(len=200) :: file_name

  integer,allocatable :: Nkvect(:)
  integer :: Nkint
  
  real(8),dimension(:),allocatable :: kx,ky,kz
  !
  real(8),dimension(:,:),allocatable :: kpath
  real(8) ::      delta_kpath(3)
  type(rgb_color),dimension(2) :: color_bands
  character(len=1),dimension(4) :: kpoint_name

  real(8) :: Eout,EoutHF,EoutLgr,E_dc,Ephn,Ekin_bare,sout,Sphn,Nphn
  real(8) :: Fout !+- free energy




  real(8) :: modk,tmpi,tmpj,checkR
  real(8) :: kmod(3)
  real(8) :: Ncell
  complex(8),dimension(:,:),allocatable :: ccij
  logical :: hartree
  real(8),allocatable,dimension(:,:) :: Ta_fat,Ni_fat
  integer :: ihyb_,jhyb_

  integer :: ix,iy,iz
  real(8),dimension(:),allocatable :: read_tmp
  real(8) :: cf_ext
  real(8) :: Ucut_off
  integer :: Ncut

  integer :: Lreal
  real(8),dimension(:),allocatable :: wreal
  complex(8),dimension(:,:),allocatable :: G_loc
  real(8),dimension(:),allocatable :: Aw

  logical :: Uradius
  logical :: H1d
  real(8) :: alphaU,deltar,hybloc
  real(8),dimension(:),allocatable :: tk
  complex(8),dimension(:,:,:),allocatable :: x_iter,x_iter_,x_iter_ir,x_iter_bare
  complex(8),dimension(:,:,:),allocatable :: hf_self_fock
  complex(8),dimension(:,:,:),allocatable :: hf_self_hartree

  complex(8),dimension(:),allocatable :: xtmp,xtmp_
  real(8) :: xphi(2)
  complex(8),dimension(2) :: xpi,xpi_
  real(8) :: CDSBseed,TRSBseed,CDSB_breaking_field,PHN_seed
  real(8) :: phn_energy
  real(8) :: ephn_g
  real(8) :: phn_ell0,Xphn_iter(4),gphn(4)
  real(8) :: phn_dyn(16),Xphn_dyn(4),Pphn_dyn(4),X2phn_dyn(4),P2phn_dyn(4)
  real(8),dimension(14) :: xr_iter
  real(8),dimension(28) :: xr_tmp

  real(8) :: Ucell,Vcell,Wcell,xi_int,UTa,UNi
  real(8) :: Evalence,Econduction,tconduction,tvalence,tt_hyb,nn_hyb,tn_hyb
  real(8) :: w0gap
  logical :: use_fsolve  
  real(8) :: fs_tol
  logical :: hf_symm  !+-> HF calculation w/ symmetric order parameter <-+!
  logical :: HF_solve !+-> HF calculation with a root-finder routine  <-+
  !logical :: hf_in
  integer :: lgr_fix_verbose
  real(8) :: lgr_fix_phase(2)

  real(8) :: phi_phase,dphi,ntot,test,err
  real(8) :: phi_abs_start,phi_abs_stop,theta_phi
  real(8),dimension(:),allocatable :: mod_phi,varphi
  integer :: ixr(3),Nop
  !
  real(8),dimension(4) :: phi_sgn
  complex(8) :: lgr_iter(2)
  real(8),dimension(2) :: lgr_iter_tmp
  real(8),dimension(1) :: lgr_iter_tmp_

  real(8),dimension(:),allocatable :: phi_list
  real(8) :: gfactor,Bfield,hop_phase,TNTN_flux,spin_seed

  !

  !+- the dynamics
  complex(8),dimension(:),allocatable :: psit
  integer :: Nt_aux,Ndyn,it,Nt
  real(8) :: tstart,tstop,time,tstep
  real(8),dimension(:),allocatable :: t_grid,t_grid_aux


  integer,dimension(:,:),allocatable :: ivec2idelta
  integer,dimension(:,:,:),allocatable :: idelta2ivec
  logical ::tns_toy,op_symm,spin_deg,enforce_cds
  
  
  !+- START MPI -+!
  call init_MPI()
  comm = MPI_COMM_WORLD
  call StartMsg_MPI(comm)
  rank = get_Rank_MPI(comm)
  master = get_Master_MPI(comm)
  
  !+- PARSE INPUT DRIVER -+!
  call parse_input_variable(Nk_x,"Nk_x","input.conf",default=51)
  call parse_input_variable(Nk_y,"Nk_y","input.conf",default=1)
  call parse_input_variable(Nhf,"Nhf","input.conf",default=100)
  call parse_input_variable(it_print,"it_print","input.conf",default=20)
  ! call parse_input_variable(Nhf_,"Nhf_","input.conf",default=10)
  call parse_input_variable(Lreal,"Lreal","input.conf",default=2000)
  call parse_input_variable(wmix,"wmix","input.conf",default=0.5d0)
  call parse_input_variable(fileRlattice,"R_unit_cell","input.conf",default='R_unit_cell.conf')
  call parse_input_variable(file_w90_hr,"w90_hr_input","input.conf",default='wannier90_hr.dat')
  call parse_input_variable(file_UV,"UVq_input","input.conf",default='VU.dat')  
  call parse_input_variable(read_tns,"READ_TNS","input.conf", &
       default='/home/mazza/project_data/TNS_EI/MALTE_DATA/K_12x12x3_B_120/')
  call parse_input_variable(plot_w90_bands,"plot_w90","input.conf",default=.true.)
  !
  call parse_input_variable(Vcell,"V","input.conf",default=1.d0)
  call parse_input_variable(UNi,"UNi","input.conf",default=1.d0)
  call parse_input_variable(UTa,"UTa","input.conf",default=1.d0)

  call parse_input_variable(Wcell,"W","input.conf",default=0.d0)
  !
  call parse_input_variable(Econduction,"Econduction","input.conf",default=1.7d0)
  call parse_input_variable(Evalence,"Evalence","input.conf",default=-0.9d0)
  call parse_input_variable(tconduction,"tconduction","input.conf",default=-0.8d0)
  call parse_input_variable(tvalence,"tvalence","input.conf",default=0.4d0)

  call parse_input_variable(tt_hyb,"tt_hyb","input.conf",default=0.0d0)
  call parse_input_variable(nn_hyb,"nn_hyb","input.conf",default=0.0d0)
  call parse_input_variable(tn_hyb,"tn_hyb","input.conf",default=0.0d0)

  
  call parse_input_variable(w0gap,"w0gap","input.conf",default=0.0d0)
  call parse_input_variable(hybloc,"hybloc","input.conf",default=0.d0)
  call parse_input_variable(use_fsolve,"use_fsolve","input.conf",default=.true.)
  call parse_input_variable(fs_tol,"fs_tol","input.conf",1.d-8)
  !
  call parse_input_variable(theta_phi,"theta_phi","input.conf",default=0d0)
  call parse_input_variable(Nop,"Nop","input.conf",default=1)
  call parse_input_variable(phi_abs_start,"phi_abs_start","input.conf",default=1.d-3)
  call parse_input_variable(phi_abs_stop,"phi_abs_end","input.conf",default=0.2d0)
  !
  call parse_input_variable(phi_sgn,"phi_sign","input.conf",default=[1.d0,-1.d0,-1.d0,1.d0])
  !
  call parse_input_variable(hf_symm,"hf_symm","input.conf",default=.false.)
  call parse_input_variable(hf_solve,"hf_solve","input.conf",default=.false.)

  call parse_input_variable(hf_conv,"hf_conv","input.conf",default=1.d-10)
  call parse_input_variable(lgr_fix_verbose,"lgr_fix_verbose","input.conf",default=0)
  call parse_input_variable(lgr_fix_phase,"lgr_fix_phase","input.conf",default=[0d0,0d0])
  !
  !call parse_input_variable(hop_phase,"hop_phase","input.conf",default=1.d-3)
  call parse_input_variable(gfactor,"gfactor","input.conf",default=0.d0)
  call parse_input_variable(Bfield,"Bfield","input.conf",default=0.d0,comment='B-field in Tesla')
  !
  call parse_input_variable(TRSBseed,"TRSBseed","input.conf",default=0.d0)
  call parse_input_variable(CDSBseed,"CDSBseed","input.conf",default=0.d0)
  call parse_input_variable(CDSB_breaking_field,"CDSB_breaking_field","input.conf",default=0.d0)
  call parse_input_variable(PHN_seed,"PHN_seed","input.conf",default=0.d0)

  
  
  !+- monoclinic distortion
  call parse_input_variable(phn_energy,"phn_energy","input.conf",default=0.01d0,comment='phonon energy in eV')
  call parse_input_variable(ephn_g,"ephn_g","input.conf",default=0.0d0,comment='electron-phonon energy in eV')
  call parse_input_variable(phn_ell0,"phn_ell0","input.conf",default=0.05d0,comment='Ta-distortion in armstrong')
  
  
  call parse_input_variable(ucut_off,"ucut_off","input.conf",default=1d0)
  call parse_input_variable(xi_int,"xi_int","input.conf",default=10000d0)

  call parse_input_variable(tns_toy,"tns_toy","input.conf",default=.false.)
  call parse_input_variable(op_symm,"op_symm","input.conf",default=.true.)
  call parse_input_variable(spin_deg,"spin_deg","input.conf",default=.true.)
  call parse_input_variable(enforce_cds,"enforce_cds","input.conf",default=.false.)
  !
  call parse_input_variable(tstart,"TSTART","input.conf",default=0.d0,comment='initial time [ps]')
  call parse_input_variable(tstop,"TSTOP","input.conf",default=10d0,comment='final time [ps]')
  call parse_input_variable(tstep,"TSTEP","input.conf",default=1.d-3)  
  !
  call get_global_vars
  !
  Econduction = Econduction + w0gap*0.5d0
  Evalence = Evalence - w0gap*0.5d0
  !
  Nso=Norb!*Nspin
  !
  !+- read primitive cell lattice vectors -+!
  file_name=reg(read_tns)//reg(fileRlattice)  
  flen=file_length(trim(file_name))
  if(flen.ne.3) stop "error in reading unit cell vectors"
  unit_in=free_unit(); 
  open(unit=unit_in,file=trim(file_name),status="old",action="read")  
  read(unit_in,*) R1(:)
  read(unit_in,*) R2(:)
  read(unit_in,*) R3(:)
  close(unit_in)
  call TB_set_ei(R1,R2,R3)
  !
  R2=0.d0
  Bk1=0.d0;Bk2=0.d0;Bk3=0.d0
  Bk1(1)=2.d0*pi/R1(1)
  Bk3(3)=2.d0*pi/R3(3)
  !+- build a monkhorst-pack grid -+!
  call build_mp_grid_2d(Nk_x,Nk_y)


  TNTN_flux=Bfield*sqrt(dot_product(R3,R3))*0.25*sqrt(dot_product(R1,R1))
  write(*,*) "TNTN_flux",TNTN_flux
  write(*,*) "TNTN_flux_over_flux_quantum",TNTN_flux/flux_quantum_TAng2
  write(*,*) "zeeman field",Bohr_magneton_in_eVoT*gfactor*0.5d0*Bfield
  hop_phase=TNTN_flux/4./(flux_quantum_TAng2/2./pi)
  write(*,*) "hop-phase",hop_phase


  !
  Lk=Nk_x*Nk_y
  Nk_z=1
  allocate(Nkvect(3));
  Nkvect(1)=Nk_x
  Nkvect(2)=Nk_y
  Nkvect(3)=Nk_z
  !
  allocate(kpt_latt(Lk,3),ik_stride(Lk,3),wtk(Lk),igr2ik(Nk_x,Nk_y,Nk_z))
  wtk=1.d0/dble(Lk)
  ik=0
  kpt_latt=0.d0
  do i=1,Nk_x
     do k=1,Nk_z           
        do j=1,Nk_y        
           !
           ik=ik+1
           !
           ik_stride(ik,1) = i
           ik_stride(ik,2) = j
           ik_stride(ik,3) = k
           igr2ik(i,j,k) = ik
           !
           do idim=1,3
              kpt_latt(ik,idim) = kpt_latt(ik,idim) + kxgrid(i)*Bk1(idim)+kygrid(j)*Bk3(idim)
           end do
           modk=sqrt(kpt_latt(ik,1)**2.d0+kpt_latt(ik,2)**2.d0+kpt_latt(ik,3)**2.d0)
           if(modk.lt.1.d-12) ik0=ik

           write(300,*) kpt_latt(ik,:)
        end do
     end do
  end do
  !
  nrpts=Nk_x*Nk_y
  write(*,*) 'nrpts',nrpts
  !
  if(allocated(irvec)) deallocate(irvec); allocate(irvec(nrpts,3))
  if(allocated(ndegen)) deallocate(ndegen); allocate(ndegen(nrpts));ndegen=1
  !
  allocate(rpt_latt(nrpts,3))
  irvec=0;ir=0
  do ix=1,Nk_x
     do iy=1,Nk_y
        ir=ir+1
        !
        irvec(ir,1) = - (Nk_x-1)/2 + (ix-1)
        irvec(ir,3) = - (Nk_y-1)/2 + (iy-1)
        rpt_latt(ir,:) = irvec(ir,1)*R1 + irvec(ir,3)*R3
        modk=sqrt(rpt_latt(ir,1)**2.d0+rpt_latt(ir,2)**2.d0+rpt_latt(ir,3)**2.d0)
        if(modk.lt.1.d-12) ir0=ir
        modk=sqrt((rpt_latt(ir,1)-R1(1))**2.d0+rpt_latt(ir,2)**2.d0+rpt_latt(ir,3)**2.d0)
        if(modk.lt.1.d-12) irR=ir
        modk=sqrt((rpt_latt(ir,1)+R1(1))**2.d0+rpt_latt(ir,2)**2.d0+rpt_latt(ir,3)**2.d0)
        if(modk.lt.1.d-12) irL=ir
        !
     end do
  end do  
  ixr(1)=ir0;ixr(2)=irR;ixr(3)=irL  
  !
  if(allocated(itmp)) deallocate(itmp);allocate(itmp(nrpts,1))  
  allocate(stride1D_(nrpts))
  i=0
  do ir=1,nrpts
     if(irvec(ir,2).eq.0.and.irvec(ir,3).eq.0) then
        i=i+1
        itmp(i,1) = irvec(ir,1)
        stride1D_(i) = ir
     end if
  end do
  nr1d=i
  write(*,*) "nr-points on x chain",nr1d
  allocate(irvec1d(nr1d,1));   allocate(stride1D(nr1d)); 
  irvec1d(:,1)=itmp(1:nr1d,1)
  stride1D=stride1D_(1:nr1d)
  !
  if(allocated(kpath)) deallocate(kpath)
  allocate(kpath(4,3))
  !+- M-point
  kpath(1,:)=0.5d0*Bk1+0.5d0*Bk3
  !+- Z-point
  kpath(2,:)=0.5d0*Bk3
  !+- G-point
  kpath(3,:)=0.d0
  !+- X-point
  kpath(4,:)=0.5d0*Bk1     
  !
  Norb=6
  Nspin=2
  Nso=Norb*Nspin
  allocate(Hk_toy(Nso,Nso,Lk),delta_hf(Nso,Nso,Lk),delta_hf_(Nso,Nso,Lk))
  allocate(tk(Norb))
  tk(1:2) = tconduction; tk(4:5) = tconduction
  tk(3) = tvalence; tk(6) = tvalence
  Hk_toy=0.d0
  do ik=1,Lk
     do ispin=1,Nspin
        !
        !+- upper chain -+!
        !
        iorb=1
        iso=(ispin-1)*Norb+iorb
        Hk_toy(iso,iso,ik) = Econduction  + 2.d0*tk(iorb)*dcos(kpt_latt(ik,1)*R1(1)+kpt_latt(ik,2)*R1(2)+kpt_latt(ik,3)*R1(3)+4d0*hop_phase)
        Hk_toy(iso,iso,ik) = Hk_toy(iso,iso,ik) + (-1d0)**dble(ispin)*Bohr_magneton_in_eVoT*gfactor*0.5d0*Bfield  !+- Zeeman field
        !
        iorb=2
        iso=(ispin-1)*Norb+iorb
        Hk_toy(iso,iso,ik) = Econduction  + 2.d0*tk(iorb)*dcos(kpt_latt(ik,1)*R1(1)+kpt_latt(ik,2)*R1(2)+kpt_latt(ik,3)*R1(3)-4d0*hop_phase)
        Hk_toy(iso,iso,ik) = Hk_toy(iso,iso,ik) + (-1d0)**dble(ispin)*Bohr_magneton_in_eVoT*gfactor*0.5d0*Bfield  !+- Zeeman field
        !
        iorb=3
        iso=(ispin-1)*Norb+iorb
        Hk_toy(iso,iso,ik) = Evalence  + 2.d0*tk(iorb)*dcos(kpt_latt(ik,1)*R1(1)+kpt_latt(ik,2)*R1(2)+kpt_latt(ik,3)*R1(3))
        Hk_toy(iso,iso,ik) = Hk_toy(iso,iso,ik) + (-1d0)**dble(ispin)*Bohr_magneton_in_eVoT*gfactor*0.5d0*Bfield  !+- Zeeman field
        !
        !C*_{Ta} C_{Ni}
        iorb=1
        jorb=3
        iso=(ispin-1)*Norb+iorb
        jso=(ispin-1)*Norb+jorb
        Rlat=R1
        Hk_toy(iso,jso,ik) = hybloc*(1.d0*exp(-xi*hop_phase)-exp(xi*dot_product(Rlat,kpt_latt(ik,:)))*exp(xi*hop_phase))
        Hk_toy(iso,jso,ik) = Hk_toy(iso,jso,ik) + CDSB_breaking_field*(1.d0+exp(xi*dot_product(Rlat,kpt_latt(ik,:))))
        Hk_toy(jso,iso,ik) = conjg(Hk_toy(iso,jso,ik))!hybloc*(1.d0-exp(-xi*dot_product(Rlat,kpt_latt(ik,:))))
        !
        iorb=2
        jorb=3
        iso=(ispin-1)*Norb+iorb
        jso=(ispin-1)*Norb+jorb
        Rlat=R1
        Hk_toy(iso,jso,ik) = hybloc*(1.d0*exp(xi*hop_phase)-exp(xi*dot_product(Rlat,kpt_latt(ik,:)))*exp(-xi*hop_phase))
        Hk_toy(iso,jso,ik) = Hk_toy(iso,jso,ik) - CDSB_breaking_field*(1.d0+exp(xi*dot_product(Rlat,kpt_latt(ik,:))))
        Hk_toy(jso,iso,ik) = conjg(Hk_toy(iso,jso,ik))!hybloc*(1.d0-exp(-xi*dot_product(Rlat,kpt_latt(ik,:))))
        !
        !+- lower chain -+!
        !
        iorb=4
        iso=(ispin-1)*Norb+iorb
        Hk_toy(iso,iso,ik) = Econduction  + 2.d0*tk(iorb)*dcos(kpt_latt(ik,1)*R1(1)+kpt_latt(ik,2)*R1(2)+kpt_latt(ik,3)*R1(3)+4d0*hop_phase)
        Hk_toy(iso,iso,ik) = Hk_toy(iso,iso,ik) + (-1d0)**dble(ispin)*Bohr_magneton_in_eVoT*gfactor*0.5d0*Bfield  !+- Zeeman field
        !
        iorb=5
        iso=(ispin-1)*Norb+iorb
        Hk_toy(iso,iso,ik) = Econduction  + 2.d0*tk(iorb)*dcos(kpt_latt(ik,1)*R1(1)+kpt_latt(ik,2)*R1(2)+kpt_latt(ik,3)*R1(3)-4d0*hop_phase)
        Hk_toy(iso,iso,ik) = Hk_toy(iso,iso,ik) + (-1d0)**dble(ispin)*Bohr_magneton_in_eVoT*gfactor*0.5d0*Bfield  !+- Zeeman field
        !
        iorb=6
        iso=(ispin-1)*Norb+iorb
        Hk_toy(iso,iso,ik) = Evalence  + 2.d0*tk(iorb)*dcos(kpt_latt(ik,1)*R1(1)+kpt_latt(ik,2)*R1(2)+kpt_latt(ik,3)*R1(3))
        Hk_toy(iso,iso,ik) = Hk_toy(iso,iso,ik) + (-1d0)**dble(ispin)*Bohr_magneton_in_eVoT*gfactor*0.5d0*Bfield  !+- Zeeman field
        !
        iorb=4
        jorb=6
        iso=(ispin-1)*Norb+iorb
        jso=(ispin-1)*Norb+jorb
        Rlat=-R1
        Hk_toy(iso,jso,ik) = hybloc*(1.d0*exp(xi*hop_phase)-exp(-xi*dot_product(Rlat,kpt_latt(ik,:)))*exp(-xi*hop_phase))
        Hk_toy(iso,jso,ik) = Hk_toy(iso,jso,ik) - CDSB_breaking_field*(1.d0+exp(-xi*dot_product(Rlat,kpt_latt(ik,:))))
        Hk_toy(jso,iso,ik) = conjg(Hk_toy(iso,jso,ik))
        !
        iorb=5
        jorb=6
        iso=(ispin-1)*Norb+iorb
        jso=(ispin-1)*Norb+jorb
        Rlat=-R1
        Hk_toy(iso,jso,ik) = hybloc*(1.d0*exp(-xi*hop_phase)-exp(-xi*dot_product(Rlat,kpt_latt(ik,:)))*exp(xi*hop_phase))
        Hk_toy(iso,jso,ik) = Hk_toy(iso,jso,ik) + CDSB_breaking_field*(1.d0+exp(-xi*dot_product(Rlat,kpt_latt(ik,:))))        
        Hk_toy(jso,iso,ik) = conjg(Hk_toy(iso,jso,ik))
        !
        !+- Ta-Ta chain hybridization -+!
        !
        iorb = 2
        jorb = 4
        iso  = (ispin-1)*Norb+iorb
        jso  = (ispin-1)*Norb+jorb
        Rlat = R1
        Hk_toy(iso,jso,ik) = Hk_toy(iso,jso,ik) + tt_hyb*(1.d0+exp(xi*dot_product(Rlat,kpt_latt(ik,:))))
        Hk_toy(jso,iso,ik) = conjg(Hk_toy(iso,jso,ik))
        ! !
        iorb = 1
        jorb = 5
        iso  = (ispin-1)*Norb+iorb
        jso  = (ispin-1)*Norb+jorb
        Rlat = -R3
        Hk_toy(iso,jso,ik) = Hk_toy(iso,jso,ik) + tt_hyb*(exp(xi*dot_product(Rlat,kpt_latt(ik,:))))
        !
        Rlat = -R3+R1
        Hk_toy(iso,jso,ik) = Hk_toy(iso,jso,ik) + tt_hyb*exp(xi*dot_product(Rlat,kpt_latt(ik,:)))
        Hk_toy(jso,iso,ik) = conjg(Hk_toy(iso,jso,ik))
        !
        
        !+- Interchain Ni-Ni hopping -+!
        iorb = 3
        jorb = 6
        iso  = (ispin-1)*Norb+iorb
        jso  = (ispin-1)*Norb+jorb
        Rlat=-R1
        Hk_toy(iso,jso,ik) = Hk_toy(iso,jso,ik) + nn_hyb*(1.d0+exp(xi*dot_product(Rlat,kpt_latt(ik,:))))
        Rlat=-R3
        Hk_toy(iso,jso,ik) = Hk_toy(iso,jso,ik) + nn_hyb*(exp(xi*dot_product(Rlat,kpt_latt(ik,:))))
        Rlat=-R3-R1
        Hk_toy(iso,jso,ik) = Hk_toy(iso,jso,ik) + nn_hyb*(exp(xi*dot_product(Rlat,kpt_latt(ik,:))))
        Hk_toy(jso,iso,ik) = conjg(Hk_toy(iso,jso,ik))
        !

        !+- Interchain Ta-Ni hopping -+!
        iorb = 2
        jorb = 6
        iso  = (ispin-1)*Norb+iorb
        jso  = (ispin-1)*Norb+jorb
        Rlat=R1
        Hk_toy(iso,jso,ik) = Hk_toy(iso,jso,ik) + tn_hyb*exp(xi*dot_product(Rlat,kpt_latt(ik,:)))
        Rlat=-R1
        Hk_toy(iso,jso,ik) = Hk_toy(iso,jso,ik) - tn_hyb*exp(xi*dot_product(Rlat,kpt_latt(ik,:)))
        Hk_toy(jso,iso,ik) = conjg(Hk_toy(iso,jso,ik))
        !
        iorb = 1
        jorb = 6
        iso  = (ispin-1)*Norb+iorb
        jso  = (ispin-1)*Norb+jorb
        Rlat=-R3+R1
        Hk_toy(iso,jso,ik) = Hk_toy(iso,jso,ik) + tn_hyb*exp(xi*dot_product(Rlat,kpt_latt(ik,:)))
        Rlat=-R3-R1
        Hk_toy(iso,jso,ik) = Hk_toy(iso,jso,ik) - tn_hyb*exp(xi*dot_product(Rlat,kpt_latt(ik,:)))
        Hk_toy(jso,iso,ik) = conjg(Hk_toy(iso,jso,ik))
                
        !+- Interchain Ta-Ni hopping -+!
        iorb = 4
        jorb = 3
        iso  = (ispin-1)*Norb+iorb
        jso  = (ispin-1)*Norb+jorb
        Rlat=  R1
        Hk_toy(iso,jso,ik) = Hk_toy(iso,jso,ik) - tn_hyb*exp(xi*dot_product(Rlat,kpt_latt(ik,:)))
        Rlat= -R1
        Hk_toy(iso,jso,ik) = Hk_toy(iso,jso,ik) + tn_hyb*exp(xi*dot_product(Rlat,kpt_latt(ik,:)))
        Hk_toy(jso,iso,ik) = conjg(Hk_toy(iso,jso,ik))
        
        iorb = 5
        jorb = 3
        iso  = (ispin-1)*Norb+iorb
        jso  = (ispin-1)*Norb+jorb
        Rlat = R1+R3
        Hk_toy(iso,jso,ik) = Hk_toy(iso,jso,ik) - tn_hyb*exp(xi*dot_product(Rlat,kpt_latt(ik,:)))
        Rlat= -R1+R3
        Hk_toy(iso,jso,ik) = Hk_toy(iso,jso,ik) + tn_hyb*exp(xi*dot_product(Rlat,kpt_latt(ik,:)))
        Hk_toy(jso,iso,ik) = conjg(Hk_toy(iso,jso,ik))
        !
     end do
     !
     !
  end do
  !
  !
  allocate(Hr_toy(Nso,Nso,nrpts))
  do ir=1,nrpts
     call FT_q2r(rpt_latt(ir,:),Hr_toy(:,:,ir),Hk_toy)
  end do
  !
  mu_fix=0.d0
  !+- change this too
  call fix_mu(Hk_toy,delta_hf,mu_fix,eout,sout)
  write(*,*) 'mu_fix',mu_fix


  allocate(kpt_path(400,3))
  allocate(ek_out_up(Norb))
  allocate(ek_out_dw(Norb))  
  allocate(Hktmp(Nso,Nso))

  uio=free_unit()
  open(unit=uio,file='tns_bare_bands.out')
  unit_in=free_unit()
  open(unit_in,file='tns_fat_bare.tns')
  allocate(Ta_fat(Norb,Nspin),Ni_fat(Norb,Nspin))
  !
  modk=0.d0
  do i=1,3
     delta_kpath=kpath(i+1,:)-kpath(i,:)
     do ik=1,100
        j=(i-1)*100 + ik
        kpt_path(j,:) = kpath(i,:) + dble(ik-1)/100.d0*delta_kpath
        modk=modk+sqrt(dot_product(1.d0/100.d0*delta_kpath,1.d0/100.d0*delta_kpath))
        !
        Hktmp=0.d0
        call FT_r2q(kpt_path(j,:),Hktmp,Hr_toy)
        !
        ek_out_up=0.0d0
        call eigh(Hktmp(1:Norb,1:Norb),ek_out_up)
        !
        ek_out_dw=0.0d0
        call eigh(Hktmp(Norb+1:2*Norb,Norb+1:2*Norb),ek_out_dw)
        !
        Ta_fat=0.d0
        do ispin=1,2
           do iorb=1,2
              iso=(ispin-1)*Norb+iorb
              Ta_fat(:,ispin)=Ta_fat(:,ispin)+abs(Hktmp(iso,(ispin-1)*Norb+1:(ispin-1)*Norb+Norb))**2.d0
           end do
           do iorb=4,5
              iso=(ispin-1)*Norb+iorb
              Ta_fat(:,ispin)=Ta_fat(:,ispin)+abs(Hktmp(iso,(ispin-1)*Norb+1:(ispin-1)*Norb+Norb))**2.d0
           end do
        end do
        Ni_fat=0.d0
        do ispin=1,2
           iorb=3
           iso=(ispin-1)*Norb+iorb
           Ni_fat(:,ispin)=Ni_fat(:,ispin)+abs(Hktmp(iso,(ispin-1)*Norb+1:(ispin-1)*Norb+Norb))**2.d0
           iorb=6
           iso=(ispin-1)*Norb+iorb
           Ni_fat(:,ispin)=Ni_fat(:,ispin)+abs(Hktmp(iso,(ispin-1)*Norb+1:(ispin-1)*Norb+Norb))**2.d0
        end do
        !
        write(uio,'(30F18.10)') modk,ek_out_up-mu_fix,ek_out_dw-mu_fix
        write(unit_in,'(30F18.10)') modk,Ta_fat(:,1),Ta_fat(:,2),Ni_fat(:,1),Ni_fat(:,2)        
        !
     end do
     !
  end do
  close(uio)
  close(unit_in)


  
  !
  allocate(delta_hfr(Nso,Nso,nrpts))
  do ir=1,nrpts
     !+- this choice is consistent w/ the interaction Hamiltonian written as (n_{R (1,2)} + n_{R+d (1,2)}) n_{R 3}
     delta_hfr(:,:,ir)=0.d0
     do ik=1,Lk
        delta_hfr(:,:,ir)=delta_hfr(:,:,ir) + &
             delta_hf(:,:,ik)*exp(xi*dot_product(rpt_latt(ir,:),kpt_latt(ik,:)))*wtk(ik)
     end do
  end do
  !
  uio=free_unit()
  open(unit=uio,file='hyb_info.out')
  iso=0
  do i=1,Norb
     do j=1,Norb
        iso=iso+1
        write(uio,*) i,j,iso+1
     end do
  end do
  close(uio)
  !+-  this changeds big times, cause the x_iter to be optimized self-consistently
  !    now becomes (~huge) of the order of Nkx (k-space)
  !    in particular, the number of variational parameters are 10*Nk*Nspin
  !    2chains. For each chain and spin there are 5 \Delta(k)
  !    \Delta1(k) = <C*_{T1}C_{T1}>
  !    \Delta2(k) = <C*_{T2}C_{T2}>
  !    \Delta3(k) = <C*_{N}C_{N}>
  !    \Delta4(k) = <C*_{N}C_{T1}>
  !    \Delta5(k) = <C*_{N}C_{T2}>
  !  
  !+- allocate and print the bare order parameters
  Nhf_opt=10
  allocate(x_iter(Lk,Nhf_opt,Nspin)); x_iter=0d0
  allocate(x_iter_(Lk,Nhf_opt,Nspin)); x_iter_=0d0

  allocate(x_iter_bare(Lk,Nhf_opt,Nspin)); x_iter_bare=0d0
  call deltak_to_xiter(delta_hf,x_iter)
  call print_xiter(x_iter,filename='bare_xiter')

  call enforce_inv_hf(x_iter,op_symm=op_symm,spin_symm=spin_deg)
  call print_xiter(x_iter,filename='bare_after_inv_x_iter_BLS')

  
  x_iter_bare=x_iter
  call init_ihfopt_strides
  call print_hyb(x_iter,filename='bare_TNShyb')

  
  !+- print the bare energy
  ntot = 0.0
  do iso=1,Nso
     ntot = ntot + delta_hfr(iso,iso,ir0)
  end do
  !
  uio=free_unit()
  open(unit=uio,file='bare_energy.out')
  write(uio, '(20F18.10)') Eout+mu_fix*ntot,Eout,sout
  close(uio)
  
  !allocate HF-hamiltonian
  allocate(H_Hf(Nso,Nso,Lk))  

  !+- define the long-range interaction
  unit_io=free_unit()
  allocate(Uss_Vs_R(Nhf_opt,nrpts),UsAs_Vs_R(Nhf_opt,nrpts))
  Uss_VS_R =0d0  !+- same-spin interaction
  UsAs_VS_R=0d0  !+- opposite-spin interaction
  xi_int=xi_int*R1(1)
  do ir=1,nrpts
     !
     Uss_VS_R(1,ir) = UTa*ucut_off*R1(1)/(sqrt(dot_product(rpt_latt(ir,:),rpt_latt(ir,:)))+ucut_off*R1(1))*exp(-abs(rpt_latt(ir,1)/xi_int))
     Uss_VS_R(2,ir) = UTa*ucut_off*R1(1)/(sqrt(dot_product(rpt_latt(ir,:),rpt_latt(ir,:)))+ucut_off*R1(1))*exp(-abs(rpt_latt(ir,1)/xi_int))
     Uss_VS_R(3,ir) = UNi*ucut_off*R1(1)/(sqrt(dot_product(rpt_latt(ir,:),rpt_latt(ir,:)))+ucut_off*R1(1))*exp(-abs(rpt_latt(ir,1)/xi_int))
     Uss_VS_R(4,ir) = Vcell*R1(1)/(abs(rpt_latt(ir,1))+abs(rpt_latt(ir,1)+R1(1)))*exp(-(abs(rpt_latt(ir,1))+abs(rpt_latt(ir,1)+R1(1))-R1(1))/xi_int)
     Uss_VS_R(5,ir) = Vcell*R1(1)/(abs(rpt_latt(ir,1))+abs(rpt_latt(ir,1)+R1(1)))*exp(-(abs(rpt_latt(ir,1))+abs(rpt_latt(ir,1)+R1(1))-R1(1))/xi_int)
     !
     Uss_VS_R(6,ir) = UTa*ucut_off*R1(1)/(sqrt(dot_product(rpt_latt(ir,:),rpt_latt(ir,:)))+ucut_off*R1(1))*exp(-abs(rpt_latt(ir,1)/xi_int))
     Uss_VS_R(7,ir) = UTa*ucut_off*R1(1)/(sqrt(dot_product(rpt_latt(ir,:),rpt_latt(ir,:)))+ucut_off*R1(1))*exp(-abs(rpt_latt(ir,1)/xi_int))
     Uss_VS_R(8,ir) = UNi*ucut_off*R1(1)/(sqrt(dot_product(rpt_latt(ir,:),rpt_latt(ir,:)))+ucut_off*R1(1))*exp(-abs(rpt_latt(ir,1)/xi_int))
     Uss_VS_R(9,ir) = Vcell*R1(1)/(abs(rpt_latt(ir,1))+abs(rpt_latt(ir,1)-R1(1)))*exp(-(abs(rpt_latt(ir,1))+abs(rpt_latt(ir,1)-R1(1))-R1(1))/xi_int)
     Uss_VS_R(10,ir) = Vcell*R1(1)/(abs(rpt_latt(ir,1))+abs(rpt_latt(ir,1)-R1(1)))*exp(-(abs(rpt_latt(ir,1))+abs(rpt_latt(ir,1)-R1(1))-R1(1))/xi_int)
     !
  end do
  !
  if(tns_toy) then
     Uss_VS_R=0d0
     !+- 
     Uss_VS_R(1,ir0) = UTa
     Uss_VS_R(2,ir0) = UTa
     Uss_VS_R(3,ir0) = UNi
     !+- 
     Uss_VS_R(4,ir0) = Vcell
     Uss_VS_R(5,ir0) = Vcell
     Uss_VS_R(4,irL) = Vcell
     Uss_VS_R(5,irL) = Vcell
     !+-
     Uss_VS_R(6,ir0) = UTa
     Uss_VS_R(7,ir0) = UTa
     Uss_VS_R(8,ir0) = UNi
     !+- 
     Uss_VS_R(9,ir0) = Vcell
     Uss_VS_R(10,ir0) = Vcell
     Uss_VS_R(9,irR) = Vcell
     Uss_VS_R(10,irR) = Vcell
     !
  end if
  UsAs_VS_R = Uss_VS_R
  !
  Uss_VS_R(1,ir0) = 0d0
  Uss_VS_R(2,ir0) = 0d0
  Uss_VS_R(3,ir0) = 0d0  
  !
  Uss_VS_R(6,ir0) = 0d0
  Uss_VS_R(7,ir0) = 0d0
  Uss_VS_R(8,ir0) = 0d0
  !
  open(unit=unit_io,file='U_VS_rpt.out') !+-
  do ir=1,nrpts
     write(unit_io,'(20F18.10)') rpt_latt(ir,1),dreal(Uss_VS_R(1:5,ir)),dreal(UsAs_VS_R(1:5,ir))
  end do
  close(unit_io)
  !
  allocate(Uss_VS_q(Nhf_opt,Lk)); Uss_VS_q=0d0
  allocate(UsAs_VS_q(Nhf_opt,Lk)); UsAs_VS_q=0d0
  open(unit=unit_io,file='U_VS_kpt.out') !+- 
  do ik=1,Lk
     do ihf=1,Nhf_opt
        call FTr2q(kpt_latt(ik,:),Uss_Vs_q(ihf,ik),Uss_VS_R(ihf,:))
        call FTr2q(kpt_latt(ik,:),UsAs_Vs_q(ihf,ik),UsAs_VS_R(ihf,:))
     end do
     write(unit_io,'(30F18.10)') kpt_latt(ik,1),dreal(Uss_VS_q(1:5,ik)),dreal(UsAs_VS_q(1:5,ik)), &
          dimag(Uss_VS_q(1:5,ik)),dimag(UsAs_VS_q(1:5,ik))
  end do
  close(unit_io)
  !+- THE UNIT CELL IS -+!
  !+  |        Ta|
  !   |-N+       |
  !   |        Ta|
  !   |   (*)    |  <---(*) the inversion point
  !   |Ta        |
  !   !      +Ni-|
  !   |Ta        |
  !
  !+- initialise the X-iter
  call init_xiter_loop(x_iter_bare,x_iter,printseed=.true.)

  !+- phonon couplings -+!
  gphn(1)=ephn_g
  gphn(2)=ephn_g
  gphn(3)=-1d0*ephn_g
  gphn(4)=-1d0*ephn_g
  
  !+- initialize phonons accordingly to init-electrons
  call xiter_ik2ir(x_iter,x_iter_ir)
  XPHN_iter=0d0
  do ispin=1,Nspin
     XPHN_iter(1) = XPHN_iter(1) - 4d0*gphn(1)/phn_energy*dreal(x_iter_ir(ir0,4,ispin)+x_iter_ir(irL,4,ispin))
     ! XPHN_iter(2) = XPHN_iter(2) - 4d0*gphn(2)/phn_energy*dreal(x_iter_ir(ir0,5,ispin)+x_iter_ir(irL,5,ispin))
     ! XPHN_iter(3) = XPHN_iter(3) - 4d0*gphn(3)/phn_energy*dreal(x_iter_ir(ir0,9,ispin)+x_iter_ir(irR,9,ispin))
     ! XPHN_iter(4) = XPHN_iter(4) - 4d0*gphn(4)/phn_energy*dreal(x_iter_ir(ir0,10,ispin)+x_iter_ir(irR,10,ispin))
  end do
  XPHN_iter(2) = -1d0*XPHN_iter(1)
  XPHN_iter(3) =  1d0*XPHN_iter(1)
  XPHN_iter(4) = -1d0*XPHN_iter(1)
  
  !+- time-dep loop
  ! uio=free_unit()
  ! open(unit=uio,file='hf_dyn_BLS_hybs_chain1_ns1.out')
  ! close(uio)
  ! open(unit=uio,file='hf_dyn_BLS_hybs_chain1_ns2.out')
  ! close(uio)
  ! open(unit=uio,file='hf_dyn_BLS_hybs_chain2_ns1.out')
  ! close(uio)
  ! open(unit=uio,file='hf_dyn_BLS_hybs_chain2_ns2.out')
  ! close(uio)

  !+-
  open(unit=uio,file='hf_dyn_order_parameter_chain1_plus.out')
  close(uio)
  open(unit=uio,file='hf_dyn_order_parameter_chain1_minus.out')
  close(uio)
  open(unit=uio,file='hf_dyn_order_parameter_chain2_plus.out')
  close(uio)
  open(unit=uio,file='hf_dyn_order_parameter_chain2_minus.out')
  close(uio)
  open(unit=uio,file='hf_dyn_energy.out')
  close(uio)
  open(unit=uio,file='hf_dyn_phonons.out')
  close(uio)
  
  !+- initialize dynamics -+!
  !+- remember - I measure time in ps.
  Nt=nint((tstop-tstart)/tstep)
  if(master) write(*,*) 'number of time steps Nt',Nt
  Nt_aux=2*Nt+1
  allocate(t_grid(Nt),t_grid_aux(Nt_aux))
  !
  t_grid = linspace(tstart,tstep*real(Nt-1,8),Nt)
  t_grid_aux = linspace(tstart,0.5d0*tstep*real(Nt_aux-1,8),Nt_aux)
  !
  Ndyn=size(delta_hf) + 16 !+- here the + 16 = 4(X,P,X^2,P^2)x4  for the phonons
  !
  !+- init tdHF strides -+!
  call init_tdHF_strides

  !+- init the dynamical vector -+!
  allocate(psit(Ndyn)); psit=0d0

  !+- solve the equilibrium hamiltonian after reading the solution -+!
  call get_double_counting_energy(x_iter,E_dc);
  call get_ni_loc(x_iter,ni_orb,ntot)
  write(500,*) E_dc,Xphn_iter
  
  H_Hf = HF_hamiltonian(x_iter,xphn_=XPHN_iter)
  H_Hf = H_Hf + Hk_toy
  call fix_mu(H_Hf,delta_hf,mu_fix,eout,sout)
  call deltak_to_xiter(delta_hf,x_iter)
  
  
  call xiter_ik2ir(x_iter,x_iter_ir)
  call get_ni_loc(x_iter,ni_orb,ntot)
  !
  call print_xiter(x_iter,filename='before_inv_x_iter_BLS')
  call print_hyb(x_iter,filename='before_hyb')

  call xiter_ik2ir(x_iter,x_iter_ir)
  write(700,'(30F18.10)') x_iter_ir(ir0,4,1),x_iter_ir(irL,4,1)
       !dreal(x_iter_ir(ir0,4,1)+x_iter_ir(irL,4,1)),dimag(x_iter_ir(ir0,4,1)),dimag(x_iter_ir(irL,4,1)), &
       !dreal(x_iter_ir(ir0,4,2)+x_iter_ir(irL,4,2)),dimag(x_iter_ir(ir0,4,2)),dimag(x_iter_ir(irL,4,2))
  write(700,'(30F18.10)') x_iter_ir(ir0,5,1),x_iter_ir(irL,5,1)
       ! dreal(x_iter_ir(ir0,5,1)+x_iter_ir(irL,5,1)),dimag(x_iter_ir(ir0,5,1)),dimag(x_iter_ir(irL,5,1)), &
       ! dreal(x_iter_ir(ir0,5,2)+x_iter_ir(irL,5,2)),dimag(x_iter_ir(ir0,5,2)),dimag(x_iter_ir(irL,5,2))

  write(700,*)
  call enforce_inv_hf(x_iter,op_symm=op_symm,spin_symm=spin_deg)
  call print_xiter(x_iter,filename='after_inv_x_iter_BLS')
  call xiter_ik2ir(x_iter,x_iter_ir)
  write(700,'(30F18.10)') x_iter_ir(ir0,4,1),x_iter_ir(irL,4,1)
             ! dreal(x_iter_ir(ir0,4,1)+x_iter_ir(irL,4,1)),dimag(x_iter_ir(ir0,4,1)),dimag(x_iter_ir(irL,4,1)), &
             ! dreal(x_iter_ir(ir0,4,2)+x_iter_ir(irL,4,2)),dimag(x_iter_ir(ir0,4,2)),dimag(x_iter_ir(irL,4,2))
  write(700,'(30F18.10)') x_iter_ir(ir0,5,1),x_iter_ir(irL,5,1)
             ! dreal(x_iter_ir(ir0,5,1)+x_iter_ir(irL,5,1)),dimag(x_iter_ir(ir0,5,1)),dimag(x_iter_ir(irL,5,1)), &
             ! dreal(x_iter_ir(ir0,5,2)+x_iter_ir(irL,5,2)),dimag(x_iter_ir(ir0,5,2)),dimag(x_iter_ir(irL,5,2))

  call get_double_counting_energy(x_iter,E_dc); eout = eout - E_dc
  write(500,*) eout+mu_fix*ntot,E_dc,eout
  stop
  phn_dyn=0d0 ![1:4  X, 5:8 P , 9:12 X^2, 13:16 P^2]
  ! NBB
  ! Xphn \equiv     (b+b*)/\sqrt(2)
  ! Pphn \equiv -xi*(b-b*)/\sqrt(2)
  !init <Xphn>
  do ispin=1,Nspin
     phn_dyn(1)  =  phn_dyn(1) - 4d0*gphn(1)/phn_energy*dreal(x_iter_ir(ir0,4,ispin)+x_iter_ir(irL,4,ispin))/sqrt(2d0)
     phn_dyn(2)  =  phn_dyn(2) - 4d0*gphn(2)/phn_energy*dreal(x_iter_ir(ir0,5,ispin)+x_iter_ir(irL,5,ispin))/sqrt(2d0)
     phn_dyn(3)  =  phn_dyn(3) - 4d0*gphn(3)/phn_energy*dreal(x_iter_ir(ir0,9,ispin)+x_iter_ir(irR,9,ispin))/sqrt(2d0)
     phn_dyn(4)  =  phn_dyn(4) - 4d0*gphn(4)/phn_energy*dreal(x_iter_ir(ir0,10,ispin)+x_iter_ir(irR,10,ispin))/sqrt(2d0)
  end do
  !init <Pphn>
  phn_dyn(5:8) = 0d0
  !init <X2phn> and <P2phn>
  do i=1,4
     phn_dyn(i+8) =   0.5d0*phn_dyn(i)**2d0 + 0.5d0*phn_dyn(i)**2d0 + 1./(exp(beta*phn_energy)-1.0) + 0.5d0  !+- <X2phn> -+!
     phn_dyn(i+12) = -0.5d0*phn_dyn(i)**2d0 + 0.5d0*phn_dyn(i)**2d0 + 1./(exp(beta*phn_energy)-1.0) + 0.5d0  !+- <P2phn> -+!
  end do
  call delta2psi(delta_hf,phn_dyn,psit)
  
  ! do ik=1,Lk  ->>>NB this is now done inside delta2psi
  !    do iso=1,Nso
  !       do jso=1,Nso
  !          isys=idelta2ivec(iso,jso,ik)
  !          psit(isys) = delta_hf(iso,jso,ik)
  !       end do
  !    end do
  ! end do
  !
  !
  !+- allocate the dynamical HF hamiltonian
  allocate(H_hf_dyn(Nso,Nso,Lk)); H_hf_dyn=0d0  
  !+-
  do it=1,Nt-1
     !
     !+- compute the dynamical HF hamiltonian -+!
     !call psi2delta(psit,delta_hf,xphn_iter)
     call psi2delta(psit,delta_hf,phn_dyn)
     xphn_iter=phn_dyn(1:4)*sqrt(2d0)
     H_Hf_dyn=HF_hamiltonian(x_iter,xphn_=xphn_iter)
     H_Hf_dyn=H_Hf_dyn+Hk_toy     
     ! !
     ! call fix_mu(H_Hf,delta_hf,mu_fix,eout,sout)
     H_HF_dyn=H_HF
     eout=0d0     !
     do ik=1,Lk
        do iso=1,Nso
           eout = eout + delta_hf(iso,iso,ik)*H_HF_dyn(iso,iso,ik)*wtk(ik)
           do jso=iso+1,Nso
              eout = eout + 2d0*dreal(delta_hf(iso,jso,ik)*H_HF_dyn(iso,jso,ik))*wtk(ik)
           end do
        end do
     end do
     !     
     ! call deltak_to_xiter(delta_hf,x_iter)
     !+- here @equ ---> enforce_cds + enforce_inv
     ! call enforce_inv_hf(x_iter,op_symm=op_symm,spin_symm=spin_deg)
     call deltak_to_xiter(delta_hf,x_iter)     
     call get_double_counting_energy(x_iter,E_dc); eout = eout - E_dc
     !+- phonon energy
     Ephn=0d0
     do i=1,4
        !Nphn=0.5*(<X2> + <P2>)
        Nphn = 0.5d0*phn_dyn(i+8)+0.5d0*phn_dyn(i+12) 
        !Nphn = 0.25*XPHN_iter(i)**2d0 + 1./(exp(beta*phn_energy)-1.0)
        Ephn = Ephn + phn_energy*Nphn
     end do
     Eout = Eout + Ephn     
     call get_ni_loc(x_iter,ni_orb,ntot)
     
     !PRINTING
     if(mod(it-1,it_print).eq.0) then
        !+- get real space for printing 
        call xiter_ik2ir(x_iter,x_iter_ir)
        uio=free_unit()
        open(unit=uio,file='hf_dyn_phonons.out',status='old',position='append')
        write(uio,'(30F18.10)') XPHN_iter(1:4),XPHN_iter(1)*phn_ell0
        close(uio)
        !
        !
        ! open(unit=uio,file='hf_dyn_BLS_hybs_chain1_ns1.out',status='old',position='append')
        ! write(uio,'(30F18.10)') x_iter_ir(ir0,1:5,1),x_iter_ir(irL,1:5,1)
        ! close(uio)
        ! open(unit=uio,file='hf_dyn_BLS_hybs_chain1_ns2.out',status='old',position='append')
        ! write(uio,'(30F18.10)') x_iter_ir(ir0,1:5,2),x_iter_ir(irL,1:5,2)
        ! close(uio)
        ! open(unit=uio,file='hf_dyn_BLS_hybs_chain2_ns1.out',status='old',position='append')
        ! write(uio,'(30F18.10)') x_iter_ir(ir0,6:10,2),x_iter_ir(irL,6:10,2)
        ! close(uio)
        ! open(unit=uio,file='hf_dyn_BLS_hybs_chain2_ns2.out',status='old',position='append')
        ! write(uio,'(30F18.10)') x_iter_ir(ir0,6:10,2),x_iter_ir(irL,6:10,2)
        ! close(uio)
        !
        !
        open(unit=uio,file='hf_dyn_order_parameter_chain1_plus.out',status='old',position='append')
        write(uio,'(30F18.10)') t_grid(it), &
             dreal(x_iter_ir(ir0,4,1)+x_iter_ir(irL,4,1)),dimag(x_iter_ir(ir0,4,1)),dimag(x_iter_ir(irL,4,1)), &
             dreal(x_iter_ir(ir0,4,2)+x_iter_ir(irL,4,2)),dimag(x_iter_ir(ir0,4,2)),dimag(x_iter_ir(irL,4,2))
        close(uio)
        
        open(unit=uio,file='hf_dyn_order_parameter_chain1_minus.out',status='old',position='append')
        write(uio,'(30F18.10)') t_grid(it), &
             dreal(x_iter_ir(ir0,5,1)+x_iter_ir(irL,5,1)),dimag(x_iter_ir(ir0,5,1)),dimag(x_iter_ir(irL,5,1)), &
             dreal(x_iter_ir(ir0,5,2)+x_iter_ir(irL,5,2)),dimag(x_iter_ir(ir0,5,2)),dimag(x_iter_ir(irL,5,2))     
        close(uio)
        
        open(unit=uio,file='hf_dyn_order_parameter_chain2_plus.out',status='old',position='append')
        write(uio,'(30F18.10)') t_grid(it), &
             dreal(x_iter_ir(ir0,9,1)+x_iter_ir(irR,9,1)),dimag(x_iter_ir(ir0,9,1)),dimag(x_iter_ir(irR,9,1)), &
             dreal(x_iter_ir(ir0,9,2)+x_iter_ir(irL,9,2)),dimag(x_iter_ir(ir0,9,2)),dimag(x_iter_ir(irL,9,2))
        close(uio)
        
        open(unit=uio,file='hf_dyn_order_parameter_chain2_minus.out',status='old',position='append')
        write(uio,'(30F18.10)') t_grid(it), &
             dreal(x_iter_ir(ir0,10,1)+x_iter_ir(irR,10,1)),dimag(x_iter_ir(ir0,10,1)),dimag(x_iter_ir(irR,10,1)), &
             dreal(x_iter_ir(ir0,10,2)+x_iter_ir(irR,10,2)),dimag(x_iter_ir(ir0,10,2)),dimag(x_iter_ir(irR,10,2))
        close(uio)
        !
        open(unit=uio,file='hf_dyn_energy.out',status='old',position='append')
        write(uio,'(30F18.10)') t_grid(it),Eout,E_dc,Ephn
        close(uio)
        !
     end if
     !
     !psit = RK_step(Ndyn,4,tstep,t_grid(it),psit,HF_eqs_of_motion)
     !     
  end do
  !+- next the post-processing for computing the response function
  

  !+- equilibrium HF loop -+!
  !   !+-
  ! unit_err=free_unit()
  ! open(unit=unit_err,file='err_BLS.out')
  ! close(unit_err)
  
  ! !+-
  ! err_hf=1.d0;iconv=0
  ! do jhf=1,Nhf
  !    !
  !    unit_err=free_unit()
  !    open(unit=unit_err,file='err_BLS.out',status='old',position='append')
  !    write(unit_err,*) err_hf           
  !    close(unit_err)
  !    !
  !    x_iter_=x_iter
  !    !
  !    !+- here I should compute the HF hamiltonian -+! 
  !    H_Hf=HF_hamiltonian(x_iter,xphn_=XPHN_iter)
  !    H_Hf=H_Hf+Hk_toy
  !    !
  !    call fix_mu(H_Hf,delta_hf,mu_fix,eout,sout)
  !    eoutHF=eout     !
  !    call deltak_to_xiter(delta_hf,x_iter)
  !    if(enforce_cds) then        
  !       call xiter_ik2ir(x_iter,x_iter_ir)
  !       do ispin=1,Nspin
  !          x_iter_ir(ir0,4,ispin) = -dreal(x_iter_ir(irL,4,ispin))   + xi*dimag(x_iter_ir(ir0,4,ispin))
  !          x_iter_ir(ir0,5,ispin) = -dreal(x_iter_ir(irL,5,ispin))   + xi*dimag(x_iter_ir(ir0,5,ispin)) 
  !          x_iter_ir(ir0,9,ispin) = -dreal(x_iter_ir(irR,9,ispin))   + xi*dimag(x_iter_ir(ir0,9,ispin))
  !          x_iter_ir(ir0,10,ispin) = -dreal(x_iter_ir(irR,10,ispin)) + xi*dimag(x_iter_ir(ir0,10,ispin))
  !       end do
  !       call xiter_ir2ik(x_iter_ir,x_iter)
  !    end if

  !    call enforce_inv_hf(x_iter,op_symm=op_symm,spin_symm=spin_deg)
  !    !
  !    call get_double_counting_energy(x_iter,E_dc); Eout = Eout - E_dc
  !    !+- phonon energy
  !    Ephn=0d0
  !    Sphn = 0d0
  !    do i=1,4        
  !       Nphn = 0.25*XPHN_iter(i)**2d0 + 1./(exp(beta*phn_energy)-1.0)
  !       Ephn = Ephn + phn_energy*Nphn
  !       Sphn = Sphn + (1+Nphn)*log(1+Nphn)-Nphn*log(Nphn) 
  !    end do
  !    Eout = Eout + Ephn
  !    Fout = Eout - temp*(Sphn+sout)
     
  !    call get_ni_loc(x_iter,ni_orb,ntot)
  !    !
  !    !+- 
     
  !    !PRINTING
  !    if(mod(jhf-1,ihf_print).eq.0) then
  !       !+- get real space for printing 
  !       call xiter_ik2ir(x_iter,x_iter_ir)
  !       uio=free_unit()
  !       open(unit=uio,file='hf_dyn_phonons.out',status='old',position='append')
  !       write(uio,'(30F18.10)') XPHN_iter(1:4),XPHN_iter(1)*phn_ell0
  !       close(uio)
  !       !
  !       open(unit=uio,file='hf_dyn_BLS_hybs_chain1_ns1.out',status='old',position='append')
  !       write(uio,'(30F18.10)') x_iter_ir(ir0,1:5,1),x_iter_ir(irL,1:5,1)
  !       close(uio)
  !       open(unit=uio,file='hf_dyn_BLS_hybs_chain1_ns2.out',status='old',position='append')
  !       write(uio,'(30F18.10)') x_iter_ir(ir0,1:5,2),x_iter_ir(irL,1:5,2)
  !       close(uio)
  !       open(unit=uio,file='hf_dyn_BLS_hybs_chain2_ns1.out',status='old',position='append')
  !       write(uio,'(30F18.10)') x_iter_ir(ir0,6:10,2),x_iter_ir(irL,6:10,2)
  !       close(uio)
  !       open(unit=uio,file='hf_dyn_BLS_hybs_chain2_ns2.out',status='old',position='append')
  !       write(uio,'(30F18.10)') x_iter_ir(ir0,6:10,2),x_iter_ir(irL,6:10,2)
  !       close(uio)
  !       !
  !       open(unit=uio,file='hf_dyn_order_parameter_chain1_plus.out',status='old',position='append')
  !       write(uio,'(30F18.10)') dreal(x_iter_ir(ir0,4,1)+x_iter_ir(irL,4,1)),dimag(x_iter_ir(ir0,4,1)),dimag(x_iter_ir(irL,4,1)), &
  !            dreal(x_iter_ir(ir0,4,2)+x_iter_ir(irL,4,2)),dimag(x_iter_ir(ir0,4,2)),dimag(x_iter_ir(irL,4,2))
  !       close(uio)
  !       open(unit=uio,file='hf_dyn_order_parameter_chain1_minus.out',status='old',position='append')
  !       write(uio,'(30F18.10)') dreal(x_iter_ir(ir0,5,1)+x_iter_ir(irL,5,1)),dimag(x_iter_ir(ir0,5,1)),dimag(x_iter_ir(irL,5,1)), &
  !            dreal(x_iter_ir(ir0,5,2)+x_iter_ir(irL,5,2)),dimag(x_iter_ir(ir0,5,2)),dimag(x_iter_ir(irL,5,2))     
  !       close(uio)
  !       open(unit=uio,file='hf_dyn_order_parameter_chain2_plus.out',status='old',position='append')
  !       write(uio,'(30F18.10)') dreal(x_iter_ir(ir0,9,1)+x_iter_ir(irR,9,1)),dimag(x_iter_ir(ir0,9,1)),dimag(x_iter_ir(irR,9,1)), &
  !            dreal(x_iter_ir(ir0,9,2)+x_iter_ir(irL,9,2)),dimag(x_iter_ir(ir0,9,2)),dimag(x_iter_ir(irL,9,2))
  !       close(uio)
  !       open(unit=uio,file='hf_dyn_order_parameter_chain2_minus.out',status='old',position='append')
  !       write(uio,'(30F18.10)') dreal(x_iter_ir(ir0,10,1)+x_iter_ir(irR,10,1)),dimag(x_iter_ir(ir0,10,1)),dimag(x_iter_ir(irR,10,1)), &
  !            dreal(x_iter_ir(ir0,10,2)+x_iter_ir(irR,10,2)),dimag(x_iter_ir(ir0,10,2)),dimag(x_iter_ir(irR,10,2))
  !       close(uio)
        
  !       Ekin_bare=0d0
  !       do ik=1,Lk
  !          do iso=1,Nso
  !             Ekin_bare=Ekin_bare+Hk_toy(iso,iso,ik)*delta_hf(iso,iso,ik)*wtk(ik)
  !             do jso=iso+1,Nso
  !                Ekin_bare=Ekin_bare+2d0*dreal(Hk_toy(iso,jso,ik)*delta_hf(iso,jso,ik))*wtk(ik)
  !             end do
  !          end do
  !       end do        
  !    end if
  !    open(unit=uio,file='hf_dyn_energy.out',status='old',position='append')
  !    write(uio,'(30F18.10)') Eout+mu_fix*ntot,Fout+mu_fix*ntot,sout,Sphn,Eout,Fout,E_dc,Ephn,Ekin_bare
  !    close(uio)
  !    !
  !    !update Xiter
  !    x_iter=x_iter*wmix+(1.d0-wmix)*x_iter_
  !    call xiter_ik2ir(x_iter,x_iter_ir)
  !    !update Xph_iter -> the phonons
  !    XPHN_iter=0d0
  !    do ispin=1,Nspin
  !       XPHN_iter(1) = XPHN_iter(1) - gphn(1)/phn_energy*dreal(x_iter_ir(ir0,4,ispin)+x_iter_ir(irL,4,ispin))
  !       XPHN_iter(2) = XPHN_iter(2) - gphn(2)/phn_energy*dreal(x_iter_ir(ir0,5,ispin)+x_iter_ir(irL,5,ispin))
  !       XPHN_iter(3) = XPHN_iter(3) - gphn(3)/phn_energy*dreal(x_iter_ir(ir0,9,ispin)+x_iter_ir(irR,9,ispin))
  !       XPHN_iter(4) = XPHN_iter(4) - gphn(4)/phn_energy*dreal(x_iter_ir(ir0,10,ispin)+x_iter_ir(irR,10,ispin))
  !    end do
  !    err_hf=get_hf_err(x_iter,x_iter_)
  !    if(err_hf.lt.hf_conv) iconv=iconv+1
  !    if(iconv.eq.2) exit
  ! end do

  ! call print_xiter(x_iter,filename='final_x_iter_BLS')
  ! call print_hyb(x_iter,filename='final_hyb')

  ! H_Hf=HF_hamiltonian(x_iter,xphn_=XPHN_iter)
  ! H_Hf=H_Hf+Hk_toy
  ! call fix_mu(H_Hf,delta_hf,mu_fix,eout)
  
  ! !+- print bands
  ! unit_in=free_unit()
  ! open(unit=unit_in,file='TNS_bands_BLS.out')
  ! uio=free_unit()
  ! open(uio,file='TNS_fat_BLS.out')  
  ! do ir=1,nrpts
  !    call FT_q2r(rpt_latt(ir,:),Hr_toy(:,:,ir),H_hf)
  ! end do
  ! modk=0.d0
  ! do i=1,3
  !    delta_kpath=kpath(i+1,:)-kpath(i,:)
  !    do ik=1,100
  !       j=(i-1)*100 + ik
  !       kpt_path(j,:) = kpath(i,:) + dble(ik-1)/100.d0*delta_kpath
  !       modk=modk+sqrt(dot_product(1.d0/100.d0*delta_kpath,1.d0/100.d0*delta_kpath))
  !       !
  !       call FT_r2q(kpt_path(j,:),Hktmp,Hr_toy)
  !       !
  !       ek_out_up=0.0d0
  !       call eigh(Hktmp(1:Norb,1:Norb),ek_out_up)
  !       !
  !       ek_out_dw=0.0d0
  !       call eigh(Hktmp(Norb+1:2*Norb,Norb+1:2*Norb),ek_out_dw)
  !       !
  !       Ta_fat=0.d0
  !       do ispin=1,2
  !          do iorb=1,2
  !             iso=(ispin-1)*Norb+iorb
  !             Ta_fat(:,ispin)=Ta_fat(:,ispin)+abs(Hktmp(iso,(ispin-1)*Norb+1:(ispin-1)*Norb+Norb))**2.d0
  !             !Ta_fat=Ta_fat+abs(Hktmp(iso,:))**2.d0
  !          end do
  !          do iorb=4,5
  !             iso=(ispin-1)*Norb+iorb
  !             Ta_fat(:,ispin)=Ta_fat(:,ispin)+abs(Hktmp(iso,(ispin-1)*Norb+1:(ispin-1)*Norb+Norb))**2.d0
  !             !Ta_fat=Ta_fat+abs(Hktmp(iso,:))**2.d0
  !          end do
  !       end do
  !       Ni_fat=0.d0
  !       do ispin=1,2
  !          iorb=3
  !          iso=(ispin-1)*Norb+iorb
  !          Ni_fat(:,ispin)=Ni_fat(:,ispin)+abs(Hktmp(iso,(ispin-1)*Norb+1:(ispin-1)*Norb+Norb))**2.d0
  !          ! Ni_fat=Ni_fat+abs(Hktmp(iso,:))**2.d0
  !          iorb=6
  !          iso=(ispin-1)*Norb+iorb
  !          Ni_fat(:,ispin)=Ni_fat(:,ispin)+abs(Hktmp(iso,(ispin-1)*Norb+1:(ispin-1)*Norb+Norb))**2.d0
  !          ! Ni_fat=Ni_fat+abs(Hktmp(iso,:))**2.d0           
  !       end do
  !       !
  !       write(unit_in,'(30F18.10)') modk,ek_out_up-mu_fix,ek_out_dw-mu_fix
  !       write(uio,'(30F18.10)') modk,Ta_fat(:,1),Ta_fat(:,2),Ni_fat(:,1),Ni_fat(:,2)        
  !       !
  !    end do
  !    !
  ! end do
  ! close(unit_in)
  ! close(uio)

  
  stop

  

  ! if(hf_in) then
  !    !allocate(x_iter(Nspin,Lk,Nhf_opt)); x_iter=0d0
  !    !
  !    flen=file_length('hf_BLS_free_final.out')
  !    if(flen.eq.14) then
  !       call read_array('hf_BLS_free_final.out',x_iter)
  !    end if
  !    !
  ! else
  !    x_iter(1) = 0.1d0
  !    x_iter(2) = 0.1d0
  !    x_iter(3) = 1.8d0
  !    !
  !    x_iter(4) = 0.1d0+xi*0.1
  !    x_iter(5) = 0.1d0-xi*0.1
  !    !  
  !    x_iter(6) = -0.1d0+xi*0.1
  !    x_iter(7) = -0.1d0-xi*0.1
  !    !
  !    x_iter(8) = 0.1d0
  !    x_iter(9) = 0.1d0
  !    x_iter(10) = 1.8d0
  !    !
  !    x_iter(11) = 0.1d0+xi*0.1
  !    x_iter(12) = 0.1d0-xi*0.1
  !    !  
  !    x_iter(13) = -0.1d0+xi*0.1
  !    x_iter(14) = -0.1d0-xi*0.1     
  ! end if
  ! !
  ! !
  ! call save_array('hf_BLS_free.init',x_iter)

  
  stop
  !
  !+- HF loop with unconstrained order parameter -+!
  
  

  
  
  !+-   
  !stop

  ! !
  ! do ihf=1,size(phi_list,1)
  !    !
  !    !+- HERE: set the target OP -+!
  !    xphi(1)=phi_list(ihf)*cos(theta_phi*pi)
  !    xphi(2)=phi_list(ihf)*sin(theta_phi*pi)
  !    write(*,*) 'fixed phi loop',ihf,xphi
  !    !
  !    lgr_iter=0d0
  !    !+-
  !    err_hf=1.d0
  !    do jhf=1,Nhf
  !       !
  !       unit_err=free_unit()
  !       open(unit=unit_err,file='err_symm.out',status='old',position='append')
  !       write(unit_err,*) err_hf           
  !       close(unit_err)
  !       !
  !       x_iter_=x_iter
  !       !
  !       lgr_iter_tmp=lgr_fix_phase
  !       if(use_fsolve) then
  !          call fsolve(fix_lgr_params,lgr_iter_tmp,tol=fs_tol,check=.false.)
  !       else
  !          call broyden1(fix_lgr_params,lgr_iter_tmp,tol=1.d-8)
  !       end if
  !       write(*,*) 'fsolve - ihf,jhf',ihf,jhf
  !       lgr_iter=lgr_iter_tmp(1)+xi*lgr_iter_tmp(2)
  !       !
  !       H_Hf=HF_hamiltonian(x_iter,lgr_iter)
  !       H_Hf=H_Hf+Hk_toy
  !       !
  !       call fix_mu(H_Hf,delta_hf,mu_fix,eout)
  !       eoutHF=eout
  !       !
  !       do i=1,3
  !          ir=ixr(i)
  !          delta_hfr(:,:,ir)=0.d0
  !          do ik=1,Lk
  !             delta_hfr(:,:,ir)=delta_hfr(:,:,ir) + &
  !                  delta_hf(:,:,ik)*exp(xi*dot_product(rpt_latt(ir,:),kpt_latt(ik,:)))*wtk(ik)
  !          end do
  !       end do
  !       !
  !       !+ -enforce the symmetry
  !       x_iter(1) = delta_hfr(1,1,ir0)+delta_hfr(1+Norb,1+Norb,ir0)
  !       x_iter(2) = delta_hfr(2,2,ir0)+delta_hfr(2+Norb,2+Norb,ir0)
  !       x_iter(3) = delta_hfr(3,3,ir0)+delta_hfr(3+Norb,3+Norb,ir0)
  !       !
  !       x_iter(4) = delta_hfr(1,3,ir0)                   !+- computed
  !       x_iter(6) = delta_hfr(1,3,irR)                   !+- computed
  !       !+- enforce symmetries -+!
  !       x_iter(5) = -dreal(x_iter(6))-xi*dimag(x_iter(4))! delta_hfr(2,3,ir0)
  !       x_iter(7) = -dreal(x_iter(4))-xi*dimag(x_iter(6))! delta_hfr(2,3,irR)
  !       !
  !       x_iter(8) = x_iter(1)   !delta_hfr(4,4,ir0)+delta_hfr(4+Norb,4+Norb,ir0)
  !       x_iter(9) = x_iter(2)   !delta_hfr(5,5,ir0)+delta_hfr(5+Norb,5+Norb,ir0)
  !       x_iter(10) = x_iter(3)  !delta_hfr(6,6,ir0)+delta_hfr(6+Norb,6+Norb,ir0)
  !       !
  !       x_iter(11) = x_iter(5) !delta_hfr(4,6,ir0)
  !       x_iter(12) = x_iter(4) !delta_hfr(5,6,ir0)
  !       !
  !       x_iter(13) = x_iter(7) !delta_hfr(4,6,irL)
  !       x_iter(14) = x_iter(6)  !delta_hfr(5,6,irL)
  !       !           
  !       Eout=Eout-Ucell*0.25d0*(dreal(x_iter(1))**2.d0+dreal(x_iter(2))**2.d0+dreal(x_iter(3))**2.d0)
  !       Eout=Eout-2d0*Vcell*(dreal(x_iter(1))*dreal(x_iter(3))+dreal(x_iter(2))*dreal(x_iter(3)))
  !       Eout=Eout + 2.d0*Vcell*(abs(x_iter(4))**2.d0+abs(x_iter(5))**2.d0)
  !       Eout=Eout + 2.d0*Vcell*(abs(x_iter(6))**2.d0+abs(x_iter(7))**2.d0)
  !       !
  !       Eout=Eout-Ucell*0.25d0*(dreal(x_iter(8))**2.d0+dreal(x_iter(9))**2.d0+dreal(x_iter(10))**2.d0)
  !       Eout=Eout-2d0*Vcell*(dreal(x_iter(8))*dreal(x_iter(10))+dreal(x_iter(9))*dreal(x_iter(10)))
  !       Eout=Eout + 2.d0*Vcell*(abs(x_iter(11))**2.d0+abs(x_iter(12))**2.d0)
  !       Eout=Eout + 2.d0*Vcell*(abs(x_iter(13))**2.d0+abs(x_iter(14))**2.d0)
  !       !
  !       !+- lagrange parameter term -+!
  !       EoutLgr = 0d0
  !       EoutLgr = EoutLgr-4d0*dreal(lgr_iter(1))*phi_sgn(1)*dreal(x_iter(4)+x_iter(6))
  !       EoutLgr = EoutLgr-4d0*dreal(lgr_iter(1))*phi_sgn(2)*dreal(x_iter(5)+x_iter(7)) 
  !       EoutLgr = EoutLgr-4d0*dreal(lgr_iter(1))*phi_sgn(3)*dreal(x_iter(11)+x_iter(13))
  !       EoutLgr = EoutLgr-4d0*dreal(lgr_iter(1))*phi_sgn(4)*dreal(x_iter(12)+x_iter(14))
  !       !
  !       EoutLgr = EoutLgr+4d0*dimag(lgr_iter(1))*phi_sgn(1)*dimag(x_iter(4)+x_iter(6))
  !       EoutLgr = EoutLgr+4d0*dimag(lgr_iter(1))*phi_sgn(2)*dimag(x_iter(5)+x_iter(7)) 
  !       EoutLgr = EoutLgr+4d0*dimag(lgr_iter(1))*phi_sgn(3)*dimag(x_iter(11)+x_iter(13))
  !       EoutLgr = EoutLgr+4d0*dimag(lgr_iter(1))*phi_sgn(4)*dimag(x_iter(12)+x_iter(14)) 
  !       !
  !       Eout=Eout+EoutLgr
  !       !
  !       ntot=dreal(x_iter(1))+dreal(x_iter(2))+dreal(x_iter(3))
  !       ntot=ntot+dreal(x_iter(8))+dreal(x_iter(9))+dreal(x_iter(10))
  !       !
  !       uio=free_unit()
  !       open(unit=uio,file='loop_fixed_order_parameter.out',status='old',position='append')
  !       write(uio, '(50F18.10)') dreal(x_iter(1:14)),dimag(x_iter(1:14)),Eout+mu_fix*ntot,EoutHF+mu_fix*ntot,EoutLgr,dreal(lgr_iter(1)),dimag(lgr_iter(2))
  !       close(uio)        
  !       !
  !       x_iter=x_iter*wmix+(1.d0-wmix)*x_iter_
  !       err_hf=0.d0
  !       do i=1,14
  !          err_hf=err_hf+abs(x_iter(i)-x_iter_(i))**2.d0
  !       end do
  !       if(err_hf.lt.hf_conv) exit
  !    end do
  !    open(unit=uio,file='loop_fixed_order_parameter.out',status='old',position='append')
  !    write(uio, '(30F18.10)')
  !    write(uio, '(30F18.10)')
  !    close(uio)
  !    !
  !    unit_io=free_unit()
  !    open(unit=unit_io,file='ENE_VS_OP.out',status='old',position='append')
  !    write(unit_io, '(10F18.10)') xphi,Eout+mu_fix*ntot
  !    close(unit_io)

     
  !    unit_io=free_unit()
  !    open(unit=unit_io,file='XITER_VS_OP.out',status='old',position='append')
  !    write(unit_io, '(50F18.10)') phi_list(ihf),dreal(x_iter(1:14)),dimag(x_iter(1:14)),Eout+mu_fix*ntot,EoutHF+mu_fix*ntot,EoutLgr,dreal(lgr_iter(1)),dimag(lgr_iter(2))
  !    close(unit_io)
  !    !        
  !    !
  ! end do
  ! close(uio)
  ! close(unit_in)     
  ! !
  ! stop
  
  ! !+- the dynamics -+!

  ! !+- here put the option to go to the HF equilibrium   
  ! Ndyn=size(delta_hf)
  ! allocate(psit(Ndyn)); psit=0d0
  
  ! call init_tdHF_strides
  ! do ik=1,Lk
  !    do iso=1,Nso
  !       do jso=1,Nso
  !          isys=idelta2ivec(iso,jso,ik)
  !          psit(isys) = delta_hf(iso,jso,ik)
  !       end do
  !    end do
  ! end do

  
  ! unit_io=free_unit()
  ! if(master) open(unit=unit_io,file='TNS_dynamics.out')

  ! do it=1,Nt-1
  !    !
  !    call set_HF_get_obs(psit,x_iter,eout)
  !    !!+- set the HF hamiltonian     

  !    !
  !    ! psit = RK_step(Ndyn,4,tstep,t_grid(it),psit,HF_eqs_of_motion)
  !    !     
  ! end do
  !
contains
  
  !+- initialize X-iter with some symmetry breaking seed
  subroutine init_xiter_loop(x_iter_in,x_iter_out,printseed)
    complex(8),dimension(Lk,Nhf_opt,Nspin),intent(in) :: x_iter_in
    complex(8),dimension(:,:,:),allocatable,intent(out) :: x_iter_out
    complex(8),dimension(:,:,:),allocatable :: x_iter_out_ir
    character(len=200) :: filein
    logical,optional :: printseed
    logical :: printseed_
    integer :: ihf,ispin,flenin
    logical :: hf_in
    real(8),dimension(3) :: tmpread
    !
    printseed_=.false.
    if(present(printseed)) printseed_=printseed
    if(allocated(x_iter_out)) deallocate(x_iter_out)
    allocate(x_iter_out(Lk,Nhf_opt,Nspin)); x_iter_out=0d0
    !
    x_iter_out=x_iter_in
    call xiter_ik2ir(x_iter_out,x_iter_out_ir)
    x_iter_out_ir(ir0,4,:) =  x_iter_out_ir(ir0,4,:) + xi*TRSBseed
    x_iter_out_ir(irL,4,:) = -dreal(x_iter_out_ir(ir0,4,:)) + CDSBseed + xi*TRSBseed
    call xiter_ir2ik(x_iter_out_ir,x_iter_out)
    call enforce_inv_hf(x_iter_out,op_symm=op_symm,spin_symm=spin_deg)
    !
    if(printseed_) call print_hyb(x_iter_out,filename='seed_TNShyb')
    !+- if present read some previous solution
    do ihf=1,Nhf_opt
       do ispin=1,Nspin
          filein="final_x_iter_BLS_spin"//reg(txtfy(ispin))//"_ihf"//reg(txtfy(ihf))//".out"
          inquire(file=filein,exist=hf_in)
          flenin=file_length(filein)
          if(hf_in.and.flenin.eq.Lk) then
             write(*,*) "reading x_iter"
             !read it
             uio=free_unit()
             open (unit=uio,file=filein)
             do ik=1,Lk
                read (uio,*) tmpread(1),tmpread(2),tmpread(3)
                x_iter_out(ik,ihf,ispin) = tmpread(2)+xi*tmpread(3)
             end do
             close(uio)
          end if
       end do
    end do
    if(printseed_) call print_xiter(x_iter_out,filename='init_x_iter_BLS')
    !
  end subroutine init_xiter_loop
  

  
  subroutine get_double_counting_energy(x_iter_in,E_dc)
    implicit none
    complex(8),dimension(Lk,Nhf_opt,Nspin),intent(in) :: x_iter_in
    real(8) :: E_dc
    complex(8) :: Ecmplx
    complex(8),dimension(:,:,:),allocatable :: hf_self_fock
    integer :: ihf,ik,ispin
    !
    E_dc=0d0
    Ecmplx=0d0
    !
    !call get_hf_self_fock(x_iter,hf_self_fock,iprint=.false.)
    call get_hf_self_fock(x_iter_in,hf_self_fock,iprint=.false.)
    do ik=1,Lk
       do ihf=1,Nhf_opt
          do ispin=1,Nspin
             Ecmplx = Ecmplx + hf_self_fock(ik,ihf,ispin)*conjg(x_iter_in(ik,ihf,ispin))*wtk(ik)
             !Ecmplx = Ecmplx + conjg(hf_self_fock(ik,ihf,ispin))*(x_iter_in(ik,ihf,ispin))*wtk(ik)
          end do
       end do
    end do
    if(abs(dimag(Ecmplx)).gt.1E-6) write(*,*) "WARNING abs(dimag(Ecmplx)).gt.1E-6",Ecmplx
    E_dc=dreal(Ecmplx)
    !
  end subroutine get_double_counting_energy
  
  function get_hf_err(x_iter_in,x_iter_inn) result(err)
    complex(8),dimension(Lk,Nhf_opt,Nspin),intent(in) ::  x_iter_in
    complex(8),dimension(Lk,Nhf_opt,Nspin),intent(in) ::  x_iter_inn
    real(8) :: err
    integer :: ik,ihf,ispin,i
    err=0d0
    do ik=1,Lk
       do ihf=1,Nhf_opt
          do ispin=1,Nspin
             err = err + abs(x_iter_in(ik,ihf,ispin)-x_iter_inn(ik,ihf,ispin))**2d0
          end do
       end do
    end do
    err=err/dble(Lk*Nspin*Nhf_opt)
  end function get_hf_err

  
  subroutine init_ihfopt_strides
    !
    allocate(iorb_to_ihf(Norb)); iorb_to_ihf=0
    iorb_to_ihf(1) = 1
    iorb_to_ihf(2) = 2
    iorb_to_ihf(3) = 3
    iorb_to_ihf(4) = 6
    iorb_to_ihf(5) = 7
    iorb_to_ihf(6) = 8
    !
    allocate(ihf_to_iorb(Nhf_opt)); ihf_to_iorb=0
    ihf_to_iorb(1) = 1
    ihf_to_iorb(2) = 2
    ihf_to_iorb(3) = 3
    ihf_to_iorb(6) = 4
    ihf_to_iorb(7) = 5
    ihf_to_iorb(8) = 6
    !
    allocate(ijorb_to_ihf(Norb,Norb)); ijorb_to_ihf=0
    ijorb_to_ihf(1,1) = 1
    ijorb_to_ihf(2,2) = 2
    ijorb_to_ihf(3,3) = 3
    ijorb_to_ihf(1,3) = 4
    ijorb_to_ihf(3,1) = 4
    ijorb_to_ihf(2,3) = 5
    ijorb_to_ihf(3,2) = 5
    !
    ijorb_to_ihf(4,4) = 6
    ijorb_to_ihf(5,5) = 7
    ijorb_to_ihf(6,6) = 8  
    ijorb_to_ihf(4,6) = 9
    ijorb_to_ihf(6,4) = 9
    ijorb_to_ihf(5,6) = 10
    ijorb_to_ihf(6,5) = 10
    !
  end subroutine init_ihfopt_strides



  subroutine enforce_inv_hf(x_iter_in,op_symm,spin_symm)
    complex(8),dimension(Lk,Nhf_opt,Nspin),intent(inout) :: x_iter_in
    logical,optional :: op_symm,spin_symm
    logical :: op_symm_,spin_symm_
    real(8) :: tmpkphase
    complex(8),dimension(:,:,:),allocatable :: x_iter_out
    integer :: ispin,ihf,ik
    !
    allocate(x_iter_out(Lk,Nhf_opt,Nspin)); 
    x_iter_out=x_iter_in
    !
    op_symm_=.false.
    if(present(op_symm)) op_symm_=op_symm
    spin_symm_=.false.
    if(present(spin_symm)) spin_symm_=spin_symm
    !
    if(op_symm_) then
       !+- enforce the order parameter symmetry of the upper chains
       do ik=1,Lk
          do ispin=1,Nspin
             ihf=2
             x_iter_out(ik,ihf,ispin) = x_iter_in(Lk+1-ik,1,ispin)
             tmpkphase=R1(1)*kpt_latt(ik,1)
             x_iter_out(ik,5,ispin) = -1.d0*exp(-xi*tmpkphase)*x_iter_in(Lk+1-ik,4,ispin)
          end do
       end do
    end if
    !+- enforce the inversion symmetry of the double chains
    do ik=1,Lk
       do ispin=1,Nspin
          !+- diagonal terms
          !+- ihf=6 !<Ta2+Ta2+>          
          !+- ihf=7 !<Ta2-Ta2->          
          !+- ihf=8 !<Ni2Ni2>
          do ihf=6,8
             if(op_symm_) then
                x_iter_out(ik,ihf,ispin) = x_iter_out(Lk+1-ik,ihf-5,ispin)
             else
                x_iter_out(ik,ihf,ispin) = x_iter_in(Lk+1-ik,ihf-5,ispin)
             end if
          end do
          !+- off-diagonal terms
          !+- ihf=9 !<Ta2+Ni2>
          ihf=9
          if(op_symm_) then
             x_iter_out(ik,ihf,ispin) = x_iter_out(Lk+1-ik,5,ispin)
          else
             x_iter_out(ik,ihf,ispin) = x_iter_in(Lk+1-ik,5,ispin)             
          end if
          !+- ihf=10 !<Ta2-Ni2>
          ihf=10
          x_iter_out(ik,ihf,ispin) = x_iter_in(Lk+1-ik,4,ispin)          
       end do
    end do
    !
    x_iter_in = x_iter_out
    if(spin_symm_) x_iter_in(:,:,2)=x_iter_in(:,:,1)
    
    ! call print_xiter(x_iter_in,'INinv_test_xiter')
    ! call print_xiter(x_iter_out,'OUTinv_test_xiter')
    
  end subroutine enforce_inv_hf

  
  subroutine xiter_ik2ir(x_iter_ik,x_iter_ir)
    implicit none
    complex(8),dimension(Lk,Nhf_opt,Nspin),intent(in) :: x_iter_ik
    complex(8),dimension(:,:,:),allocatable,intent(out) :: x_iter_ir
    integer :: uio,ihf,ispin,ir    
    if(allocated(x_iter_ir)) deallocate(x_iter_ir) 
    allocate(x_iter_ir(nrpts,Nhf_opt,Nspin)); x_iter_ir=0d0    
    do ihf=1,Nhf_opt
       do ispin=1,Nspin
          do ir=1,nrpts
             call FTq2r(rpt_latt(ir,:),x_iter_ir(ir,ihf,ispin),x_iter_ik(:,ihf,ispin))             
          end do
       end do
    end do
  end subroutine xiter_ik2ir
  subroutine xiter_ir2ik(x_iter_ir,x_iter_ik)
    implicit none
    complex(8),dimension(nrpts,Nhf_opt,Nspin),intent(in) :: x_iter_ir
    complex(8),dimension(:,:,:),allocatable,intent(out) :: x_iter_ik
    integer :: uio,ihf,ispin,ir,ik    
    if(allocated(x_iter_ik)) deallocate(x_iter_ik) 
    allocate(x_iter_ik(Lk,Nhf_opt,Nspin)); x_iter_ik=0d0    
    do ihf=1,Nhf_opt
       do ispin=1,Nspin
          do ik=1,Lk
             call FTr2q(kpt_latt(ik,:),x_iter_ik(ik,ihf,ispin),x_iter_ir(:,ihf,ispin))             
          end do
       end do
    end do
  end subroutine xiter_ir2ik


  !+- print the X-iter
  subroutine print_xiter(x_iter_in,filename)
    complex(8),dimension(Lk,Nhf_opt,Nspin),intent(in) :: x_iter_in
    integer :: ispin,ihf,ik,uio
    character(len=*),optional :: filename
    character(len=100) :: filename_    
    filename_="x_iter"
    if(present(filename)) filename_=filename    
    do ispin=1,Nspin
       do ihf=1,Nhf_opt
          uio=free_unit()
          open(unit=uio,file=reg(filename_)//"_spin"//reg(txtfy(ispin))//"_ihf"//reg(txtfy(ihf))//".out")
          do ik=1,Lk
             write(uio, '(20F18.10)') kpt_latt(ik,1),dreal(x_iter(ik,ihf,ispin)),dimag(x_iter(ik,ihf,ispin))
          end do
          close(uio)
       end do
    end do
  end subroutine print_xiter

  
  !+- write down the bare hybridizations
  subroutine print_hyb(x_iter_in,filename)
    implicit none
    complex(8),dimension(Lk,Nhf_opt,Nspin),intent(in) :: x_iter_in
    complex(8),dimension(:,:,:),allocatable :: x_iter_ir
    character(len=*),optional :: filename
    character(len=100) :: filename_
    character(len=100) :: filehyb
    integer :: uio,ihf,ispin,ir

    
    filename_="TNShyb"
    if(present(filename)) filename_=filename

    allocate(x_iter_ir(nrpts,Nhf_opt,Nspin)); x_iter_ir=0d0
    ! write(*,*) size(x_iter_in,1),Lk,nrpts
    ! write(*,*) size(x_iter_in,2),Nhf_opt
    ! write(*,*) size(x_iter_in,3),Nspin    
    ! write(*,*) size(x_iter_ir,1),Lk,nrpts
    ! write(*,*) size(x_iter_ir,2),Nhf_opt
    ! write(*,*) size(x_iter_ir,3),Nspin
    ! write(*,*) size(rpt_latt)

    do ihf=1,Nhf_opt
       do ispin=1,Nspin
          uio=free_unit()
          filehyb=reg(filename_)//"_vs_ir_ihf"//reg(txtfy(ihf))//"_ispin"//reg(txtfy(ispin))//".out"          
          open(unit=uio,file=filehyb)
          do ir=1,nrpts
             !call FT_q2r(rpt_latt(ir,:),Hr_toy(:,:,ir),Hk_toy)
             !x_iter(ik,ihf,ispin))
             call FTq2r(rpt_latt(ir,:),x_iter_ir(ir,ihf,ispin),x_iter_in(:,ihf,ispin))             
             Rlat=irvec1d(ir,1)*R1
             write(uio,'(10F18.10)') Rlat(1),x_iter_ir(ir,ihf,ispin)
          end do
          close(uio)
       end do
    end do

  end subroutine print_hyb

  

  subroutine get_hf_self_hartree(x_iter_in,hf_self_hartree_out,iprint)
    implicit none
    complex(8),dimension(Lk,Nhf_opt,Nspin),intent(in) :: x_iter_in
    complex(8),dimension(:,:,:),allocatable,intent(out) :: hf_self_hartree_out
    logical,optional :: iprint
    logical :: iprint_        
    complex(8),dimension(:,:),allocatable :: ni_orb_tmp
    integer :: ispin,ihf,ik,iik,jk
    real(8) :: ktmp
    !
    iprint_=.false.
    if(present(iprint)) iprint_=iprint
    if(allocated(hf_self_hartree_out)) deallocate(hf_self_hartree_out)
    allocate(hf_self_hartree_out(Lk,Norb,Nspin)); hf_self_hartree_out=0d0
    call get_ni_loc(x_iter,ni_orb_tmp)
    !
    do ispin=1,Nspin
       do iorb=1,Norb
          uio=free_unit()
          if(iprint_) open(unit=uio,file="hf_self_hartree_spin"//reg(txtfy(ispin))//"_iorb"//reg(txtfy(iorb))//".out")        
          do ik=1,Lk
             hf_self_hartree_out(ik,iorb,ispin)=0d0
             do jspin=1,Nspin
                do jorb=1,Norb
                   ihf=ijorb_to_ihf(iorb,jorb)
                   if(ihf.gt.0) then
                      if(ispin.eq.jspin) then
                         hf_self_hartree_out(ik,iorb,ispin) = hf_self_hartree_out(ik,iorb,ispin) + &
                              Uss_VS_q(ihf,ik0)*ni_orb_tmp(jspin,jorb)
                      else
                         hf_self_hartree_out(ik,iorb,ispin) = hf_self_hartree_out(ik,iorb,ispin) + &
                              UsAs_VS_q(ihf,ik0)*ni_orb_tmp(jspin,jorb)
                      end if
                   end if
                end do
             end do
             if(iprint_) write(uio, '(40F18.10)') kpt_latt(ik,1),dreal(hf_self_hartree_out(ik,iorb,ispin)),dimag(hf_self_hartree_out(ik,iorb,ispin))              
          end do
          if(iprint_) close(uio)
       end do
    end do
  end subroutine get_hf_self_hartree
  !
  subroutine get_hf_self_fock(x_iter_in,hf_self_fock_out,iprint)
    implicit none
    complex(8),dimension(Lk,Nhf_opt,Nspin),intent(in) :: x_iter_in
    complex(8),dimension(:,:,:),allocatable,intent(out) :: hf_self_fock_out
    logical,optional :: iprint
    logical :: iprint_
    integer :: ispint,ihf,ik,iik,jk
    real(8) :: ktmp
    !
    iprint_=.false.
    if(present(iprint)) iprint_=iprint
    if(allocated(hf_self_fock_out)) deallocate(hf_self_fock_out)
    allocate(hf_self_fock_out(Lk,Nhf_opt,Nspin)); hf_self_fock_out=0d0
    
    do ispin=1,Nspin
       do ihf=1,Nhf_opt          
          uio=free_unit() 
          if(iprint_) open(unit=uio,file="hf_self_fock_spin"//reg(txtfy(ispin))//"_ihf"//reg(txtfy(ihf))//".out")
          do ik=1,Lk        
             !
             hf_self_fock_out(ik,ihf,ispin)=0d0
             !           
             do jk=1,Lk              
                if(ik.ne.jk.or.tns_toy) then
                   ktmp=kxgrid(ik)-kxgrid(jk)
                   do while(ktmp.lt.-0.5d0)
                   !do while(ktmp.lt.kxgrid(1))
                      ktmp=ktmp+1.d0
                   end do
                   do while(ktmp.gt.0.5d0)
                   !do while(ktmp.gt.kxgrid(Lk))
                      ktmp=ktmp-1.d0
                   end do
                   iik = 1+nint((ktmp-kxgrid(1))/(kxgrid(2)-kxgrid(1)))
                   if(iik.lt.1.or.iik.gt.Lk) then
                      write(*,*) iik,ktmp,kxgrid(1),Lk
                      stop "(iik.lt.1.or.iik.gt.Lk)"
                   end if
                   !
                   hf_self_fock_out(ik,ihf,ispin) = hf_self_fock_out(ik,ihf,ispin) - &
                        Uss_VS_q(ihf,iik)*x_iter(jk,ihf,ispin)*wtk(jk)
                end if
             end do
             if(iprint_) write(uio, '(40F18.10)') kpt_latt(ik,1),dreal(hf_self_fock_out(ik,ihf,ispin)),dimag(hf_self_fock_out(ik,ihf,ispin))
          end do
          if(iprint_) close(uio)
       end do
    end do    
    !
  end subroutine get_hf_self_fock
  !
  subroutine get_ni_loc(x_iter_in,ni_out,ntot)
    complex(8),dimension(Lk,Nhf_opt,Nspin),intent(in) :: x_iter_in
    complex(8),dimension(:,:),allocatable,intent(out) :: ni_out
    integer :: iorb,ispin,ihf
    real(8),optional :: ntot
    real(8) :: ntot_
    if(allocated(ni_out)) deallocate(ni_out)
    allocate(ni_out(Nspin,Norb)); ni_out=0d0
    !
    ntot_=0d0
    do iorb=1,Norb
       do ispin=1,Nspin
          do ik=1,Lk
             ihf=iorb_to_ihf(iorb)
             ni_out(ispin,iorb) = ni_out(ispin,iorb) + x_iter_in(ik,ihf,ispin)*wtk(ik) !
          end do
          ntot_=ntot_+ni_out(ispin,iorb) 
       end do
    end do
    if(present(ntot)) ntot=ntot_
    !
  end subroutine get_ni_loc

  
  subroutine deltak_to_xiter(deltak_in,x_iter_out)
    complex(8),dimension(Nso,Nso,Lk),intent(in) :: deltak_in
    complex(8),dimension(Lk,Nhf_opt,Nspin),intent(out) :: x_iter_out
    integer :: i,j,ih,ik
    !
    !x_iter(jk,ihf,ispin)
    x_iter_out = 0d0
    !
    x_iter_out(:,1,1) = deltak_in(1,1,:)
    x_iter_out(:,1,2) = deltak_in(1+Norb,1+Norb,:)
    !
    x_iter_out(:,2,1) = deltak_in(2,2,:)
    x_iter_out(:,2,2) = deltak_in(2+Norb,2+Norb,:)
    !
    x_iter_out(:,3,1) = deltak_in(3,3,:)
    x_iter_out(:,3,2) = deltak_in(3+Norb,3+Norb,:)
    !
    x_iter_out(:,4,1) = deltak_in(1,3,:)
    x_iter_out(:,4,2) = deltak_in(1+Norb,3+Norb,:)
    !
    x_iter_out(:,5,1) = deltak_in(2,3,:)
    x_iter_out(:,5,2) = deltak_in(2+Norb,3+Norb,:)
    !
    x_iter_out(:,6,1) = deltak_in(4,4,:)
    x_iter_out(:,6,2) = deltak_in(4+Norb,4+Norb,:)
    !
    x_iter_out(:,7,1) = deltak_in(5,5,:)
    x_iter_out(:,7,2) = deltak_in(5+Norb,5+Norb,:)
    !
    x_iter_out(:,8,1) = deltak_in(6,6,:)
    x_iter_out(:,8,2) = deltak_in(6+Norb,6+Norb,:)
    !
    x_iter_out(:,9,1) = deltak_in(4,6,:)
    x_iter_out(:,9,2) = deltak_in(4+Norb,6+Norb,:)
    !
    x_iter_out(:,10,1) = deltak_in(5,6,:)
    x_iter_out(:,10,2) = deltak_in(5+Norb,6+Norb,:)
  end subroutine deltak_to_xiter


  

  ! subroutine set_HF_get_obs(y,x_iter,eout)
  !   implicit none
  !   complex(8),intent(inout)  :: x_iter(14)
  !   real(8),intent(inout)     :: eout
  !   complex(8),dimension(:),allocatable,intent(in) :: y
  !   complex(8),dimension(:,:,:),allocatable :: delta_hf_dyn
  !   integer :: iso,jso,ik,isys,i,ir
  !   !
  !   if(.not.allocated(y)) stop "psi2delta not allocated"
  !   if(size(y).ne.Ndyn)  stop "psi2delta size(y).ne.Ndyn"
  !   !
  !   call psi2delta(y,delta_hf_dyn)
  !   !
  !   do i=1,3
  !      ir=ixr(i)
  !      delta_hfr(:,:,ir)=0.d0
  !      do ik=1,Lk
  !         delta_hfr(:,:,ir)=delta_hfr(:,:,ir) + &
  !              delta_hf_dyn(:,:,ik)*exp(xi*dot_product(rpt_latt(ir,:),kpt_latt(ik,:)))*wtk(ik)
  !      end do
  !   end do
  !   !
  !   x_iter(1) = delta_hfr(1,1,ir0)+delta_hfr(1+Norb,1+Norb,ir0)
  !   x_iter(2) = delta_hfr(2,2,ir0)+delta_hfr(2+Norb,2+Norb,ir0)
  !   x_iter(3) = delta_hfr(3,3,ir0)+delta_hfr(3+Norb,3+Norb,ir0)
  !   !
  !   x_iter(4) = delta_hfr(1,3,ir0)
  !   x_iter(5) = delta_hfr(2,3,ir0)
  !   !
  !   x_iter(6) = delta_hfr(1,3,ir0)+delta_hfr(1,3,irR)
  !   x_iter(7) = delta_hfr(2,3,ir0)+delta_hfr(2,3,irR)
  !   !
  !   x_iter(8) = delta_hfr(4,4,ir0)+delta_hfr(4+Norb,4+Norb,ir0)
  !   x_iter(9) = delta_hfr(5,5,ir0)+delta_hfr(5+Norb,5+Norb,ir0)
  !   x_iter(10) = delta_hfr(6,6,ir0)+delta_hfr(6+Norb,6+Norb,ir0)
  !   !
  !   x_iter(11) = delta_hfr(4,6,ir0)
  !   x_iter(12) = delta_hfr(5,6,ir0)
  !   !
  !   x_iter(13) = delta_hfr(4,6,ir0)+delta_hfr(4,6,irL)
  !   x_iter(14) = delta_hfr(5,6,ir0)+delta_hfr(5,6,irL)
  !   !
  !   ! set the dynamical HF hamiltonian
  !   H_hf_dyn=HF_hamiltonian(x_iter)
  !   H_hf_dyn=H_hf_dyn+Hk_toy
  !   !
  !   eout=0d0
  !   do ik=1,Lk
  !      do iso=1,Nso
  !         eout = eout + H_hf_dyn(iso,iso,ik)*delta_hf_dyn(iso,iso,ik)*wtk(ik) 
  !         do jso=iso+1,Nso
  !            eout = eout + 2d0*dreal(H_hf_dyn(iso,jso,ik)*delta_hf_dyn(iso,jso,ik))*wtk(ik) 
  !         end do
  !      end do
  !   end do
  !   Eout=Eout-Ucell*0.25d0*(dreal(x_iter(1))**2.d0+dreal(x_iter(2))**2.d0+dreal(x_iter(3))**2.d0)
  !   Eout=Eout-2d0*Vcell*(dreal(x_iter(1))*dreal(x_iter(3))+dreal(x_iter(2))*dreal(x_iter(3)))
  !   Eout=Eout + 2.d0*Vcell*(abs(x_iter(4))**2.d0+abs(x_iter(5))**2.d0)
  !   Eout=Eout + 2.d0*Vcell*(abs(x_iter(6)-x_iter(4))**2.d0+abs(x_iter(7)-x_iter(5))**2.d0)
  !   !
  !   Eout=Eout-Ucell*0.25d0*(dreal(x_iter(8))**2.d0+dreal(x_iter(9))**2.d0+dreal(x_iter(10))**2.d0)
  !   Eout=Eout-2d0*Vcell*(dreal(x_iter(8))*dreal(x_iter(10))+dreal(x_iter(9))*dreal(x_iter(10)))
  !   !
  !   Eout=Eout + 2.d0*Vcell*(abs(x_iter(11))**2.d0+abs(x_iter(12))**2.d0)
  !   Eout=Eout + 2.d0*Vcell*(abs(x_iter(13)-x_iter(11))**2.d0+abs(x_iter(14)-x_iter(12))**2.d0)
  !   !
  ! end subroutine set_HF_get_obs
  !
  !  
  !
  ! function fix_lgr_params(x) result(out_x)
  !   implicit none
  !   real(8),dimension(:),intent(in)   :: x
  !   real(8),dimension(size(x)) :: out_x    
  !   complex(8),dimension(:),allocatable :: x_hf
  !   complex(8),dimension(:,:,:),allocatable :: Hhf_lgr,delta_hfr_lgr
  !   complex(8),dimension(:,:),allocatable :: dHK
  !   real(8) :: xmu_lgr
  !   complex(8) :: lgrx(2)
  !   real(8) :: phi_out(2),phi_CDSB_tmp(4),phi_TRSB_tmp(4,2)
  !   integer :: ik,iorb,jorb
  !   real(8) :: eout
    
  !   if(size(x).ne.2) stop "size(x).ne.2  fix lgr params"
  !   allocate(x_hf(size(x_iter))); x_hf=x_iter
  !   if(.not.allocated(Hk_toy)) stop "fix_lgr: Hk_toy not allocated"        
  !   allocate(Hhf_lgr(Nso,Nso,Lk)); Hhf_lgr=0.d0           
  !   !
  !   lgrx=x(1)+xi*x(2)
  !   Hhf_lgr = HF_hamiltonian(x_iter,lgrx)
  !   Hhf_lgr = Hhf_lgr + Hk_toy
  !   !
  !   call fix_mu(Hhf_lgr,delta_hf,mu_fix)
  !   !
  !   allocate(delta_hfr_lgr(Nso,Nso,nrpts)); delta_hfr_lgr=0d0
  !   do i=1,3
  !      ir=ixr(i)
  !      delta_hfr_lgr(:,:,ir)=0.d0
  !      do ik=1,Lk
  !         delta_hfr_lgr(:,:,ir)=delta_hfr_lgr(:,:,ir) + &
  !              delta_hf(:,:,ik)*exp(xi*dot_product(rpt_latt(ir,:),kpt_latt(ik,:)))*wtk(ik)
  !      end do
  !   end do
  !   !
  !   phi_CDSB_tmp(1) = dreal(delta_hfr_lgr(1,3,ir0)+delta_hfr_lgr(1,3,irR))
  !   phi_TRSB_tmp(1,1) = dimag(delta_hfr_lgr(1,3,ir0))
  !   phi_TRSB_tmp(1,2) = dimag(delta_hfr_lgr(1,3,irR))
  !   !
  !   phi_CDSB_tmp(2) = dreal(delta_hfr_lgr(2,3,ir0)+delta_hfr_lgr(2,3,irR))
  !   phi_TRSB_tmp(2,1) = dimag(delta_hfr_lgr(2,3,ir0))
  !   phi_TRSB_tmp(2,2) = dimag(delta_hfr_lgr(2,3,irR))
  !   !
  !   phi_CDSB_tmp(3) = dreal(delta_hfr_lgr(4,6,ir0)+delta_hfr_lgr(4,6,irL))
  !   phi_TRSB_tmp(3,1) = dimag(delta_hfr_lgr(4,6,ir0))
  !   phi_TRSB_tmp(3,2) = dimag(delta_hfr_lgr(4,6,irL))
  !   !
  !   phi_CDSB_tmp(4) = dreal(delta_hfr_lgr(5,6,ir0)+delta_hfr_lgr(5,6,irL))
  !   phi_TRSB_tmp(4,1) = dimag(delta_hfr_lgr(5,6,ir0))
  !   phi_TRSB_tmp(4,2) = dimag(delta_hfr_lgr(5,6,irL))
  !   !
  !   phi_out(1) = phi_CDSB_tmp(1)
  !   phi_out(2) = phi_TRSB_tmp(1,1)
  !   !
  !   if(master.and.lgr_fix_verbose>0) then
  !      write(*,'(A6,10F18.10)') 'lgr',lgrx
  !      write(*,'(A6,10F18.10)') 'CDSB',phi_CDSB_tmp
  !      write(*,'(A6,10F18.10)') 'TRSB',phi_TRSB_tmp(1,:),phi_TRSB_tmp(2,:),phi_TRSB_tmp(3,:),phi_TRSB_tmp(4,:)
  !      write(*,'(A6,10F18.10)') 'out',phi_out,xphi,sqrt((phi_out(1)-xphi(1))**2.d0+(phi_out(2)-xphi(2))**2.d0)
  !      write(*,*)
  !      write(*,*)
  !   end if
  !   !
  !   phi_out = phi_out - xphi
  !   out_x=phi_out
  !   ! out_x(1) = dreal(phi_out)
  !   ! out_x(2) = dimag(phi_out)
  !   !    
  ! end function fix_lgr_params


  subroutine delta2psi(delta,xphn,y)
    integer :: ik,iso,jso,isys
    complex(8),dimension(:,:,:),allocatable,intent(in) :: delta
    real(8),dimension(16),intent(in) :: xphn
    complex(8),dimension(:),allocatable,intent(inout) :: y

    if(.not.allocated(delta)) stop "delta2psi not allocated"
    if(size(delta,1).ne.Nso)  stop "delta2psi size(delta,1).ne.Nso"
    if(size(delta,2).ne.Nso)  stop "delta2psi size(delta,2).ne.Nso"
    if(size(delta,3).ne.Lk)   stop "delta2psi size(delta,3).ne.Lk"
    !
    if(allocated(y)) deallocate(y)
    allocate(y(Ndyn));y=0d0    
    do ik=1,Lk
       do iso=1,Nso
          do jso=1,Nso
             isys=idelta2ivec(iso,jso,ik)
             y(isys) = delta(iso,jso,ik)
          end do
       end do
    end do
    if(isys.ne.Ndyn-16) stop "if(isys.ne.Ndyn-16)"
    y(isys+1:Ndyn) = xphn
  end subroutine delta2psi
  subroutine psi2delta(y,delta,xphn)
    integer :: ik,iso,jso,isys
    complex(8),dimension(:,:,:),allocatable,intent(inout) :: delta
    real(8),dimension(16),intent(inout) :: xphn
    complex(8),dimension(:),allocatable,intent(in) :: y

    if(.not.allocated(y)) stop "psi2delta not allocated"
    if(size(y).ne.Ndyn)  stop "psi2delta size(y).ne.Ndyn"
    !
    if(allocated(delta)) deallocate(delta)
    allocate(delta(Nso,Nso,Lk));delta=0d0
    do isys=1,Ndyn-16
       iso=ivec2idelta(isys,1)
       jso=ivec2idelta(isys,2)
       ik=ivec2idelta(isys,3)
       delta(iso,jso,ik) = y(isys)       
    end do
    !write(*,*) isys,isys+1-Ndyn,Lk*Nso*Nso
    xphn(1:16) = dreal(y(isys:Ndyn))
  end subroutine psi2delta
  
  subroutine init_tdHF_strides
    integer :: ik,iso,jso,isys
    !
    allocate(ivec2idelta(Lk*Nso*Nso,3)) ; ivec2idelta=0
    allocate(idelta2ivec(Nso,Nso,Lk)) ; idelta2ivec=0
    !
    isys=0
    do ik=1,Lk
       do iso=1,Nso
          do jso=1,Nso
             isys=isys+1
             !
             ivec2idelta(isys,1) = iso 
             ivec2idelta(isys,2) = jso 
             ivec2idelta(isys,3) = ik
             !
             idelta2ivec(iso,jso,ik) = isys
             !
          end do
       end do
    end do
    !
  end subroutine init_tdHF_strides



  function HF_eqs_of_motion(time,yin,Nsys) result(yout)
    integer :: Nsys
    real(8) :: time
    complex(8),dimension(Nsys) :: yin,yout
    integer :: ik,jk,iso,jso,iso_,jso_,isys,jsys,jsys_
    integer ::  checkNsys,unit,it_aux
    do isys=1,Ndyn
       !
       iso = ivec2idelta(isys,1)
       jso = ivec2idelta(isys,2)
       ik  = ivec2idelta(isys,3)       
       !
       yout(iso) = 0d0
       do jso_=1,Nso
          !
          jsys=idelta2ivec(iso,jso_,ik)
          yout(iso) = yout(iso) + H_hf_dyn(jso,jso_,ik)*yin(jsys)
          !
          jsys=idelta2ivec(jso_,jso,ik)          
          yout(iso) = yout(iso) - H_hf_dyn(jso_,iso,ik)*yin(jsys)
          !
       end do
    end do
    yout=-xi*yout/(hbar_ev_ps) !+- here the time is measured in ps 
    !
  end function HF_eqs_of_motion
  
  !
  !call get_hf_self_fock(x_iter,hf_self_fock,iprint=.true.)
  !subroutine get_hf_self_fock(x_iter_in,hf_self_fock_out,iprint)

  function HF_hamiltonian(x_iter,xphn_,lgr_) result(Hhf)
    implicit none
    complex(8),dimension(Lk,Nhf_opt,Nspin),intent(in) :: x_iter
    complex(8),dimension(:,:,:),allocatable :: hf_self_fock
    !
    complex(8),dimension(Nso,Nso,Lk) :: Hhf
    real(8),optional :: xphn_(4)
    real(8)          :: xphn(4)
    
    complex(8),optional :: lgr_(2)
    complex(8)          :: lgr(2)
    integer :: iso,jso,iorb,jorb
    !
    lgr=0d0;
    if(present(lgr_)) lgr=lgr_
    xphn=0d0;
    if(present(xphn_)) xphn=xphn_
    call get_hf_self_fock(x_iter,hf_self_fock,iprint=.false.)
    Hhf=0d0
    do ik=1,Lk
       do ispin=1,Nspin
          do iorb=1,Norb
             do jorb=iorb,Norb
                iso=(ispin-1)*Norb+iorb
                jso=(ispin-1)*Norb+jorb
                !
                ihf = ijorb_to_ihf(iorb,jorb)
                if(ihf.ne.0) then
                   Hhf(jso,iso,ik) = hf_self_fock(ik,ihf,ispin)
                   Hhf(iso,jso,ik) = conjg(Hhf(jso,iso,ik))
                end if
                !
             end do
          end do
          !+- e-ph coupling -+!
          iorb=1
          jorb=3          
          iso=(ispin-1)*Norb+iorb
          jso=(ispin-1)*Norb+jorb
          !+-
          Rlat=R1
          Hhf(iso,jso,ik) = Hhf(iso,jso,ik) + gphn(1)*xphn(1)*(1.d0+exp(xi*dot_product(Rlat,kpt_latt(ik,:))))
          Hhf(jso,iso,ik) = conjg(Hhf(iso,jso,ik))
          !
          iorb=2
          jorb=3          
          iso=(ispin-1)*Norb+iorb
          jso=(ispin-1)*Norb+jorb
          Rlat=R1
          Hhf(iso,jso,ik) = Hhf(iso,jso,ik) + gphn(2)*xphn(2)*(1.d0+exp(xi*dot_product(Rlat,kpt_latt(ik,:))))
          Hhf(jso,iso,ik) = conjg(Hhf(iso,jso,ik))
          !
          iorb=4
          jorb=6          
          iso=(ispin-1)*Norb+iorb
          jso=(ispin-1)*Norb+jorb
          Rlat=R1
          Hhf(iso,jso,ik) = Hhf(iso,jso,ik) + gphn(3)*xphn(3)*(1.d0+exp(-xi*dot_product(Rlat,kpt_latt(ik,:))))
          Hhf(jso,iso,ik) = conjg(Hhf(iso,jso,ik))
          !
          iorb=5
          jorb=6          
          iso=(ispin-1)*Norb+iorb
          jso=(ispin-1)*Norb+jorb
          Rlat=R1
          Hhf(iso,jso,ik) = Hhf(iso,jso,ik) + gphn(4)*xphn(4)*(1.d0+exp(-xi*dot_product(Rlat,kpt_latt(ik,:))))
          Hhf(jso,iso,ik) = conjg(Hhf(iso,jso,ik))
          !          
       end do
    end do
    !
  end function HF_hamiltonian
  !
  !
  subroutine build_mp_grid_2d(Nx,Ny)
    implicit none
    integer :: Nx,Ny
    integer :: Nx_,Ny_
    integer :: i,ix,iy,ik,iz,unitk
    real(8) :: kmax,kxmax,kymax,kzmax,dx,dy,dz
    real(8),dimension(:),allocatable :: kxr,kyr,kzr
    !
    !    
    if(allocated(kxgrid)) deallocate(kxgrid)
    if(allocated(kygrid)) deallocate(kygrid)
    allocate(kxgrid(Nx),kygrid(Ny))
    kmax=dble(Nx-1)*0.5d0/dble(Nx)
    kxgrid = linspace(-kmax,kmax,Nx,mesh=dx)
    kxmax=kmax+dx*0.5d0
    !
    kmax=dble(Ny-1)*0.5d0/dble(Ny)
    kygrid = linspace(-kmax,kmax,Ny,iend=.false.,mesh=dy)
    kymax=kmax+dy*0.5d0
    !
    ! kmax=dble(Nz-1)*0.5d0/dble(Nz)
    ! kzgrid = linspace(-kmax,kmax,Nz,mesh=dz)
    ! kzmax=kmax+dz*0.5d0
    !
    !
    if(mod(Nx,2).eq.0) then
       kxgrid=kxgrid-0.5d0/dble(Nx)
    end if
    if(mod(Ny,2).eq.0) then
       kygrid=kygrid-0.5d0/dble(Ny)
    end if
    !
    unitk=free_unit()
    open(unit=unitk,file='mp_grid.out')
    do ix=1,Nx
       do iy=1,Ny
          write(unitk,'(10F18.10)') kxgrid(ix),kygrid(iy)
       end do
    end do
    close(unitk)
  end subroutine build_mp_grid_2d



end program Officina

