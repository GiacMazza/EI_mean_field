program officina
  !
  USE SCIFOR
  USE DMFT_TOOLS
  !
  USE VARS_GLOBAL
  USE HF_real
  !
  !
  USE MPI
  !
  !
  implicit none
  integer :: Nint

  real(8) :: KKx,KKy
  complex(8),dimension(:,:,:),allocatable :: delta_hf,H_hf,delta_hf_,Hhf_grid,delta_hfr,delta_hfr_
  complex(8),dimension(:,:,:),allocatable :: Hr_w90,Hr_w90_tmp,Hr_toy  
  real(8) :: mu_fix
  real(8) :: Uintra,Uinter,Delta_CF
  real(8) :: Uq,thop,err_hf,hf_conv
  logical :: plot_w90_bands
  !
  real(8),dimension(3) :: ktest
  real(8),dimension(:),allocatable :: ek_out,ek_out_up,ek_out_dw

  complex(8),dimension(:,:),allocatable :: Hsb

  complex(8),dimension(:),allocatable :: Uft

  complex(8),dimension(:,:),allocatable :: Hloc,Hktmp
  complex(8),dimension(:,:,:),allocatable :: Hk_w90,Hk_w90_tmp,Hk_toy
  complex(8),dimension(:,:,:,:,:),allocatable :: Hk_w90_reshape,Hk_hf_reshape
  real(8),dimension(:,:),allocatable :: kpts,kpt_path
  integer,dimension(:,:),allocatable :: irvec2d,itmp
  integer,dimension(:),allocatable :: stride2D,stride2D_

  integer,dimension(:,:),allocatable :: irvec1d
  integer,dimension(:),allocatable :: stride1D,stride1D_
  

  integer :: nr2d,nr1d

  integer :: iso,jso,ispin,jspin,iorb,jorb,ik
  integer :: i,j,k,idim
  integer :: Nhf,Nhf_,ihf,unit_err,unit_obs,unit_in,uio,Nobs,jhf
  integer,dimension(:),allocatable :: units_loc_obs
  integer :: flen,iread

  integer :: ir,jr,kr,iir
  real(8),dimension(3) :: Rlat 

  complex(8),dimension(:),allocatable :: obs,obs_loc
  real(8) :: wmix
  real(8),dimension(:,:),allocatable :: Uloc_TNS
  complex(8),dimension(:,:,:),allocatable :: Uq_TNS,Vq_TNS,Ur_TNS,Ur_tmp
  complex(8),dimension(:,:),allocatable :: check_Uloc

  real(8),dimension(2) :: Uq_read,Vq_read

  character(len=100) :: fileRlattice,file_w90_hr,file_UV,read_tns
  character(len=200) :: file_name

  integer,allocatable :: Nkvect(:)


  real(8),dimension(:),allocatable :: kx,ky,kz
  !



  real(8),dimension(:,:),allocatable :: kpath
  real(8) ::      delta_kpath(3)
  type(rgb_color),dimension(2) :: color_bands
  character(len=1),dimension(4) :: kpoint_name

  real(8) :: Eout,E_dc
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
  complex(8),dimension(:),allocatable :: x_iter,x_iter_,x_iter_bare
  complex(8),dimension(9) :: xtmp,xtmp_
  complex(8) :: xphi,xphi_
  complex(8),dimension(2) :: xpi,xpi_
  real(8),dimension(14) :: xr_iter
  real(8),dimension(18) :: xr_tmp

  real(8) :: Ucell,Vcell,Wcell
  real(8) :: Evalence,Econduction,tconduction,tvalence,tt_hyb,nn_hyb,tn_hyb
  real(8) :: w0gap
  logical :: fix_phi  !+-> HF calculation with fixed order parameter <-+!
  logical :: hf_symm  !+-> HF calculation w/ symmetric order parameter <-+!
  logical :: HF_solve !+-> HF calculation with a root-finder routine  <-+
  logical :: hf_in_up,hf_in_dw

  real(8) :: phi_start,phi_end,dphi,ntot,test,err
  integer :: Nphi
  integer :: ixr(3)
  real(8) :: gfactor,Bfield,hop_phase,TNTN_flux,spin_seed
  logical :: enforce_spin,hf_hartree
  !

  !+- START MPI -+!
  call init_MPI()
  comm = MPI_COMM_WORLD
  call StartMsg_MPI(comm)
  rank = get_Rank_MPI(comm)
  master = get_Master_MPI(comm)
  !+- PARSE INPUT DRIVER -+!
  call parse_input_variable(Nk_x,"Nk_x","input.conf",default=51)
  call parse_input_variable(Nk_y,"Nk_y","input.conf",default=51)
  call parse_input_variable(Nhf,"Nhf","input.conf",default=100)
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
  call parse_input_variable(Ucell,"U","input.conf",default=1.d0)
  call parse_input_variable(Wcell,"W","input.conf",default=0.d0)
  !
  call parse_input_variable(Econduction,"Econduction","input.conf",default=1.7d0)
  call parse_input_variable(Evalence,"Evalence","input.conf",default=-0.9d0)
  call parse_input_variable(tconduction,"tconduction","input.conf",default=-0.8d0)
  call parse_input_variable(tvalence,"tvalence","input.conf",default=0.4d0)
  !
  call parse_input_variable(tt_hyb,"tt_hyb","input.conf",default=0.0d0)
  call parse_input_variable(nn_hyb,"nn_hyb","input.conf",default=0.0d0)
  call parse_input_variable(tn_hyb,"tn_hyb","input.conf",default=0.0d0)
  !
  call parse_input_variable(w0gap,"w0gap","input.conf",default=0.0d0)
  call parse_input_variable(hybloc,"hybloc","input.conf",default=0.d0)
  call parse_input_variable(fix_phi,"fix_phi","input.conf",default=.false.)
  call parse_input_variable(phi_start,"phi_start","input.conf",default=0.5d0)
  call parse_input_variable(phi_end,"phi_end","input.conf",default=0.d0)
  call parse_input_variable(Nphi,"Nphi","input.conf",default=50)
  !
  call parse_input_variable(hf_symm,"hf_symm","input.conf",default=.false.)
  call parse_input_variable(hf_solve,"hf_solve","input.conf",default=.false.)
  call parse_input_variable(hf_conv,"hf_conv","input.conf",default=1.d-10)
  !
  call parse_input_variable(gfactor,"gfactor","input.conf",default=2.d0)
  call parse_input_variable(Bfield,"Bfield","input.conf",default=1.d0,comment='B-field in Tesla')
  call parse_input_variable(spin_seed,"spin_seed","input.conf",default=0.d-2)
  call parse_input_variable(enforce_spin,"enforce_spin","input.conf",default=.true.)
  call parse_input_variable(hf_hartree,"hf_hartree","input.conf",default=.true.)


  !
  !
  call get_global_vars

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
  call TB_get_bk(Bk1,Bk2,Bk3)
  Bkinv(:,1) = Bk1
  Bkinv(:,2) = Bk2
  Bkinv(:,3) = Bk3
  call inv(Bkinv)
  !+ k-space path
  allocate(kpath(3,3))
  !+- X-point
  kpath(1,:)=0.5d0*Bk1
  !+- G-point
  kpath(2,:)=0.d0
  !+- X-point
  kpath(3,:)=0.5d0*Bk1
  !+- build a monkhorst-pack grid -+!
  call build_mp_grid_2d(Nk_x,Nk_y)
  !
  Lk=Nk_x*Nk_y
  Nk_z=1
  allocate(Nkvect(3));
  Nkvect(1)=Nk_x
  Nkvect(2)=Nk_y
  Nkvect(3)=Nk_z

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
           write(300,*) kpt_latt(ik,:)
        end do
     end do
  end do
  !
  file_name=reg(read_tns)//reg(file_w90_hr)  
  !+-  read the w90 output -> this is just to have the w90 hamiltonian in hand <-  -+!
  allocate(Hloc(Nso,Nso))
  call TB_hr_to_hk(R1,R2,R3,hk_w90,Hloc,trim(file_name),1,6,1,Nkvect)  
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
  
  allocate(Hr_w90(Nso,Nso,nrpts))
  do ir=1,nrpts
     call FT_q2r(rpt_latt(ir,:),Hr_w90(:,:,ir),Hk_w90)
  end do

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
  R2=0.d0
  Bk1=0.d0;Bk2=0.d0;Bk3=0.d0
  Bk1(1)=2.d0*pi/R1(1)
  Bk3(3)=2.d0*pi/R3(3)
  deallocate(kpath); allocate(kpath(4,3))
  !+- M-point
  kpath(1,:)=0.5d0*Bk1+0.5d0*Bk3
  !+- Z-point
  kpath(2,:)=0.5d0*Bk3
  !+- G-point
  kpath(3,:)=0.d0
  !+- X-point
  kpath(4,:)=0.5d0*Bk1
  !
  !estimate flux through the Ta-Ni-Ta-Ni loops
  TNTN_flux=Bfield*sqrt(dot_product(R3,R3))*0.25*sqrt(dot_product(R1,R1))
  write(*,*) "TNTN_flux",TNTN_flux
  write(*,*) "TNTN_flux_over_flux_quantum",TNTN_flux/flux_quantum_TAng2
  write(*,*) "zeeman field",Bohr_magneton_in_eVoT*gfactor*0.5d0*Bfield
  hop_phase=TNTN_flux/4./(flux_quantum_TAng2/2./pi)
  write(*,*) "hop-phase",hop_phase
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
        Hk_toy(iso,iso,ik) = Econduction  + 2.d0*tk(iorb)*dcos(kpt_latt(ik,1)*R1(1)+kpt_latt(ik,2)*R1(2)+kpt_latt(ik,3)*R1(3)+4.d0*hop_phase)
        Hk_toy(iso,iso,ik) = Hk_toy(iso,iso,ik) + (-1d0)**dble(ispin)*Bohr_magneton_in_eVoT*gfactor*0.5d0*Bfield  !+- Zeeman field
        !
        iorb=2
        iso=(ispin-1)*Norb+iorb
        Hk_toy(iso,iso,ik) = Econduction  + 2.d0*tk(iorb)*dcos(kpt_latt(ik,1)*R1(1)+kpt_latt(ik,2)*R1(2)+kpt_latt(ik,3)*R1(3)-4.d0*hop_phase)
        Hk_toy(iso,iso,ik) = Hk_toy(iso,iso,ik) + (-1d0)**dble(ispin)*Bohr_magneton_in_eVoT*gfactor*0.5d0*Bfield  !+- Zeeman field
        !
        iorb=3
        iso=(ispin-1)*Norb+iorb
        Hk_toy(iso,iso,ik) = Evalence  + 2.d0*tk(iorb)*dcos(kpt_latt(ik,1)*R1(1)+kpt_latt(ik,2)*R1(2)+kpt_latt(ik,3)*R1(3))
        Hk_toy(iso,iso,ik) = Hk_toy(iso,iso,ik) + (-1d0)**dble(ispin)*Bohr_magneton_in_eVoT*gfactor*0.5d0*Bfield  !+- Zeeman field
        !
        iorb=1
        jorb=3
        iso=(ispin-1)*Norb+iorb
        jso=(ispin-1)*Norb+jorb
        Rlat=R1
        Hk_toy(iso,jso,ik) = hybloc*(1.d0*exp(-xi*hop_phase)-exp(xi*dot_product(Rlat,kpt_latt(ik,:)))*exp(xi*hop_phase)) 
        Hk_toy(jso,iso,ik) = conjg(Hk_toy(iso,jso,ik))!hybloc*(1.d0-exp(-xi*dot_product(Rlat,kpt_latt(ik,:))))
        !
        iorb=2
        jorb=3
        iso=(ispin-1)*Norb+iorb
        jso=(ispin-1)*Norb+jorb
        Rlat=R1
        Hk_toy(iso,jso,ik) = hybloc*(1.d0*exp(xi*hop_phase)-exp(xi*dot_product(Rlat,kpt_latt(ik,:)))*exp(-xi*hop_phase)) 
        Hk_toy(jso,iso,ik) = conjg(Hk_toy(iso,jso,ik))!hybloc*(1.d0-exp(-xi*dot_product(Rlat,kpt_latt(ik,:))))
        !
        !+- lower chain -+!
        !
        iorb=4
        iso=(ispin-1)*Norb+iorb
        Hk_toy(iso,iso,ik) = Econduction  + 2.d0*tk(iorb)*dcos(kpt_latt(ik,1)*R1(1)+kpt_latt(ik,2)*R1(2)+kpt_latt(ik,3)*R1(3)+4.d0*hop_phase)
        Hk_toy(iso,iso,ik) = Hk_toy(iso,iso,ik) + (-1d0)**dble(ispin)*Bohr_magneton_in_eVoT*gfactor*0.5d0*Bfield  !+- Zeeman field
        iorb=5
        iso=(ispin-1)*Norb+iorb
        Hk_toy(iso,iso,ik) = Econduction  + 2.d0*tk(iorb)*dcos(kpt_latt(ik,1)*R1(1)+kpt_latt(ik,2)*R1(2)+kpt_latt(ik,3)*R1(3)-4.d0*hop_phase)
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
        Hk_toy(iso,jso,ik) = hybloc*(1.d0*exp(xi*hop_phase)-exp(xi*dot_product(Rlat,kpt_latt(ik,:)))*exp(-xi*hop_phase))
        Hk_toy(jso,iso,ik) = conjg(Hk_toy(iso,jso,ik))
        !
        iorb=5
        jorb=6
        iso=(ispin-1)*Norb+iorb
        jso=(ispin-1)*Norb+jorb
        Rlat=-R1
        Hk_toy(iso,jso,ik) = hybloc*(1.d0*exp(-xi*hop_phase)-exp(xi*dot_product(Rlat,kpt_latt(ik,:)))*exp(xi*hop_phase))
        Hk_toy(jso,iso,ik) = conjg(Hk_toy(iso,jso,ik))
        !
        !
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
  allocate(Hr_toy(Nso,Nso,nrpts))
  do ir=1,nrpts
     call FT_q2r(rpt_latt(ir,:),Hr_toy(:,:,ir),Hk_toy)
  end do
  !
  mu_fix=0.d0
  call fix_mu(Hk_toy,delta_hf,mu_fix,eout)
  write(*,*) 'mu_fix',mu_fix
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
  uio=free_unit()
  open(unit=uio,file='hyb_bareR.out')
  allocate(obs(Norb*Norb))

  modk=0.d0
  do ir=1,nr1d
     !
     Rlat=irvec1d(ir,1)*R1
     !
     iso=0
     do i=1,Norb
        do j=1,Norb
           iso=iso+1
           obs(iso) = delta_hfr(i,j,stride1D(ir))
        end do
     end do
     !
     write(unit_in,'(100F10.5)') Rlat(1),dreal(obs(:)),dimag(obs(:))
     !
  end do
  ! 
  !+- plot bare bands -+!
  allocate(kpt_path(400,3))
  allocate(ek_out_up(Norb))
  allocate(ek_out_dw(Norb))  
  allocate(Hktmp(Nso,Nso))
  !
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
  call get_xiter(delta_hfr,x_iter)
  open(unit=uio,file='bare_xiter.out')
  do iso=1,14
     write(uio,*) x_iter(iso),x_iter(iso+(Nspin-1)*14)     
  end do
  allocate(x_iter_bare(size(x_iter))); x_iter_bare=x_iter
  close(uio)
  call get_ntot(x_iter,ntot)
  uio=free_unit()
  open(unit=uio,file='bare_energy.out')
  write(uio, '(20F18.10)') Eout+mu_fix*ntot,Eout,ntot,mu_fix
  close(uio)
  !
  !
  allocate(H_Hf(Nso,Nso,Lk))
  !
  write(*,*) "starting the symmetric HF-loop"
  inquire(file='hf_symm_final_up.out',exist=hf_in_up)
  inquire(file='hf_symm_final_dw.out',exist=hf_in_dw)  
  if(hf_in_up.and.hf_in_dw) then
     !
     flen=file_length('hf_symm_final_up.out')
     if(flen.eq.14) then
        call read_array('hf_symm_final_up.out',x_iter(1:14))
     end if
     flen=file_length('hf_symm_final_dw.out')
     if(flen.eq.14) then
        call read_array('hf_symm_final_dw.out',x_iter(15:28))
     end if
     !
  else
     !
     !+- seeds
     x_iter(1) = 0.05d0+spin_seed
     x_iter(1+14) = 0.05d0-spin_seed
     !
     x_iter(2) = 0.05d0+spin_seed
     x_iter(2+14) = 0.05d0-spin_seed
     !
     x_iter(3) = 0.9d0+spin_seed
     x_iter(3+14) = 0.9d0-spin_seed
     
     x_iter(4) = 0.1d0+spin_seed   !<c*_{Ta1,i} c_{Ni,1}>
     x_iter(4+14) = 0.1d0-spin_seed   !<c*_{Ta1,i} c_{Ni,1}>
     x_iter(6) = -0.1d0-spin_seed  !
     x_iter(6+14) = -0.1d0+spin_seed  !
     call enforce_inv_hf(x_iter)   
     !
  end if
  allocate(x_iter_(size(x_iter))); x_iter_=x_iter
  !
  call save_array('hf_symm_up.init',x_iter(1:14))
  call save_array('hf_symm_dw.init',x_iter(1+14:14+14))
  !
  uio=free_unit()
  open(unit=uio,file='loop_phi_symm_ns1.out')
  close(uio)
  open(unit=uio,file='loop_phi_symm_ns2.out')
  close(uio)
  open(unit=uio,file='loop_energy_symm.out')
  close(uio)
  !
  unit_err=free_unit()
  open(unit=unit_err,file='err_symm.out')
  close(unit_err)
  !
  err_hf=1.d0
  do ihf=1,Nhf
     write(*,*) "symmetric HF: loop nr",ihf,'/',Nhf
     unit_err=free_unit()
     open(unit=unit_err,file='err_symm.out',status='old',position='append')
     write(unit_err,*) err_hf
     close(unit_err)
     !
     x_iter_=x_iter
     !
     H_Hf=HF_hamiltonian(x_iter)
     H_Hf=H_Hf+Hk_toy
     !
     call fix_mu(H_Hf,delta_hf,mu_fix)
     !
     do i=1,3
        ir=ixr(i)
        delta_hfr(:,:,ir)=0.d0
        do ik=1,Lk
           delta_hfr(:,:,ir)=delta_hfr(:,:,ir) + &
                delta_hf(:,:,ik)*exp(xi*dot_product(rpt_latt(ir,:),kpt_latt(ik,:)))*wtk(ik)
        end do
     end do     
     call get_xiter(delta_hfr,x_iter)
     ! call get_ntot(x_iter,ntot)
     ! write(*,*) 'ntot before symmetry',ntot
     ! do i=1,14
     !    write(*,*) x_iter(i),x_iter(i+14)
     ! end do
     ! write(*,*)
     ! write(*,*)
     ! write(*,*)
     call enforce_inv_hf(x_iter)
     ! do i=1,14
     !    write(*,*) x_iter(i),x_iter(i+14)
     ! end do
     ! !
     ! call get_ntot(x_iter,ntot)
     ! write(*,*) 'ntot before mixing',ntot
     ! stop
     x_iter=x_iter*wmix+(1.d0-wmix)*x_iter_             
     !     
     H_Hf=HF_hamiltonian(x_iter) 
     H_Hf=H_Hf+Hk_toy
     call fix_mu(H_Hf,delta_hf,mu_fix,eout,iprint=.false.) !+- this is here because I mix the x_iter before writing - better writing before mixing. fix later      
     !
     call get_double_counting_energy(x_iter,E_dc); Eout = Eout + E_dc
     call get_ntot(x_iter,ntot)
     ! write(*,*) 'ntot',ntot
     ! stop
     !
     uio=free_unit()
     open(unit=uio,file='loop_phi_symm_ns1.out',status='old',position='append')     
     write(uio, '(50F18.10)') dreal(x_iter(1:14)),dimag(x_iter(1:14))
     close(uio)
     open(unit=uio,file='loop_phi_symm_ns2.out',status='old',position='append')     
     write(uio, '(50F18.10)') dreal(x_iter(1+14:14+14)),dimag(x_iter(1+14:14+14))
     close(uio)
     open(unit=uio,file='loop_energy_symm.out',status='old',position='append')     
     write(uio, '(50F18.10)') Eout+mu_fix*ntot,Eout,E_dc
     close(uio)
     err_hf=0.d0
     do i=1,28
        err_hf=err_hf+abs(x_iter(i)-x_iter_(i))**2.d0
     end do
     if(err_hf.lt.hf_conv) exit
     !stop
     ! 
  end do
  call save_array('hf_symm_final_up.out',x_iter(1:14))
  call save_array('hf_symm_final_dw.out',x_iter(15:28))
  !stop
  !
  !
  !
  !
  H_Hf=HF_hamiltonian(x_iter)
  H_Hf=H_Hf+Hk_toy
  call fix_mu(H_Hf,delta_hf,mu_fix,eout)  
  !+- plot real space hybridizations -+!
  uio=free_unit()
  open(unit=uio,file='hyb_TNS_VS_r_symm.out')
  obs=0.d0
  do ir=1,nrpts
     delta_hfr(:,:,ir)=0.d0
     do ik=1,Lk
        delta_hfr(:,:,ir)=delta_hfr(:,:,ir) + &
             delta_hf(:,:,ik)*exp(xi*dot_product(rpt_latt(ir,:),kpt_latt(ik,:)))*wtk(ik)
     end do
  end do
  do ir=1,nr1d
     Rlat=irvec1d(ir,1)*R1
     iso=0
     do i=1,Norb
        do j=1,Norb
           iso=iso+1
           obs(iso) = delta_hfr(i,j,stride1D(ir))
        end do
     end do
     !
     write(uio,'(100F10.5)') Rlat(1),dreal(obs(:)),dimag(obs(:))
     !
  end do
  close(uio)
  !
  !
  unit_in=free_unit()
  open(unit=unit_in,file='TNS_bands_symm.out')
  uio=free_unit()
  open(uio,file='TNS_fat_symm.out')  
  do ir=1,nrpts
     call FT_q2r(rpt_latt(ir,:),Hr_toy(:,:,ir),H_hf)
  end do
  modk=0.d0
  do i=1,3
     delta_kpath=kpath(i+1,:)-kpath(i,:)
     do ik=1,100
        j=(i-1)*100 + ik
        kpt_path(j,:) = kpath(i,:) + dble(ik-1)/100.d0*delta_kpath
        modk=modk+sqrt(dot_product(1.d0/100.d0*delta_kpath,1.d0/100.d0*delta_kpath))
        !
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
              !Ta_fat=Ta_fat+abs(Hktmp(iso,:))**2.d0
           end do
           do iorb=4,5
              iso=(ispin-1)*Norb+iorb
              Ta_fat(:,ispin)=Ta_fat(:,ispin)+abs(Hktmp(iso,(ispin-1)*Norb+1:(ispin-1)*Norb+Norb))**2.d0
              !Ta_fat=Ta_fat+abs(Hktmp(iso,:))**2.d0
           end do
        end do
        Ni_fat=0.d0
        do ispin=1,2
           iorb=3
           iso=(ispin-1)*Norb+iorb
           Ni_fat(:,ispin)=Ni_fat(:,ispin)+abs(Hktmp(iso,(ispin-1)*Norb+1:(ispin-1)*Norb+Norb))**2.d0
           ! Ni_fat=Ni_fat+abs(Hktmp(iso,:))**2.d0
           iorb=6
           iso=(ispin-1)*Norb+iorb
           Ni_fat(:,ispin)=Ni_fat(:,ispin)+abs(Hktmp(iso,(ispin-1)*Norb+1:(ispin-1)*Norb+Norb))**2.d0
           ! Ni_fat=Ni_fat+abs(Hktmp(iso,:))**2.d0           
        end do
        !
        write(unit_in,'(30F18.10)') modk,ek_out_up-mu_fix,ek_out_dw-mu_fix
        write(uio,'(30F18.10)') modk,Ta_fat(:,1),Ta_fat(:,2),Ni_fat(:,1),Ni_fat(:,2)        
        !
        !
     end do
     !
  end do
  close(unit_in)
  close(uio)
  !
  !
  !
  x_iter=x_iter_bare
  inquire(file='hf_BLS_final_up.out',exist=hf_in_up)
  inquire(file='hf_BLS_final_dw.out',exist=hf_in_dw)  
  if(hf_in_up.and.hf_in_dw) then
     !
     flen=file_length('hf_BLS_final_up.out')
     if(flen.eq.14) then
        call read_array('hf_BLS_final_up.out',x_iter(1:14))
     end if
     flen=file_length('hf_BLS_final_dw.out')
     if(flen.eq.14) then
        call read_array('hf_BLS_final_dw.out',x_iter(15:28))
     end if
     !
  else
     !+- symmetry breaking seeds -+!
     x_iter(1)    = 0.05d0+spin_seed
     x_iter(1+14) = 0.05d0-spin_seed
     x_iter(2)    = 0.05d0+spin_seed
     x_iter(2+14) = 0.05d0-spin_seed
     !
     x_iter(3)    = 0.9d0+spin_seed
     x_iter(3+14) = 0.9d0-spin_seed
     !
     x_iter(4)    =  0.1d0+xi*(0.1+spin_seed)
     x_iter(4+14) =  0.1d0+xi*(0.1-spin_seed)
     x_iter(6)    = -0.1d0+xi*(0.1+spin_seed)
     x_iter(6+14) = -0.1d0+xi*(0.1-spin_seed)
     call enforce_inv_hf(x_iter)   
     !
  end if
  !
  call save_array('hf_BLS_up.init',x_iter(1:14))
  call save_array('hf_BLS_dw.init',x_iter(1+14:14+14))  
  !
  !
  uio=free_unit()
  open(unit=uio,file='loop_phi_BLS_ns1.out')
  close(uio)
  open(unit=uio,file='loop_phi_BLS_ns2.out')
  close(uio)
  open(unit=uio,file='loop_energy_BLS.out')
  close(uio)
  !
  unit_err=free_unit()
  open(unit=unit_err,file='err_BLS.out')
  close(unit_err)
  !
  err_hf=1.0
  do ihf=1,Nhf
     !
     write(*,*) "BLS HF: loop nr",ihf,'/',Nhf
     unit_err=free_unit()
     open(unit=unit_err,file='err_BLS.out',status='old',position='append')
     write(unit_err,*) err_hf
     close(unit_err)
     x_iter_=x_iter
     !
     H_Hf=HF_hamiltonian(x_iter)
     H_Hf=H_Hf+Hk_toy
     !
     call fix_mu(H_Hf,delta_hf,mu_fix)
     !
     do i=1,3!
        ir=ixr(i)
        delta_hfr(:,:,ir)=0.d0
        do ik=1,Lk
           delta_hfr(:,:,ir)=delta_hfr(:,:,ir) + &
                delta_hf(:,:,ik)*exp(xi*dot_product(rpt_latt(ir,:),kpt_latt(ik,:)))*wtk(ik)
        end do
     end do
     call get_xiter(delta_hfr,x_iter)
     call enforce_inv_hf(x_iter)   
     !
     x_iter=x_iter*wmix+(1.d0-wmix)*x_iter_             
     !     
     H_Hf=HF_hamiltonian(x_iter)
     H_Hf=H_Hf+Hk_toy
     call fix_mu(H_Hf,delta_hf,mu_fix,eout)   !+- to be fixed later
     
     call get_double_counting_energy(x_iter,E_dc); Eout = Eout + E_dc
     call get_ntot(x_iter,ntot)
     uio=free_unit()
     open(unit=uio,file='loop_phi_BLS_ns1.out',status='old',position='append')     
     write(uio, '(50F18.10)') dreal(x_iter(1:14)),dimag(x_iter(1:14))
     close(uio)
     open(unit=uio,file='loop_phi_BLS_ns2.out',status='old',position='append')     
     write(uio, '(50F18.10)') dreal(x_iter(1+14:14+14)),dimag(x_iter(1+14:14+14))
     close(uio)
     open(unit=uio,file='loop_energy_BLS.out',status='old',position='append')     
     write(uio, '(50F18.10)') Eout+mu_fix*ntot,Eout,E_dc
     close(uio)
     err_hf=0.d0
     do i=1,14
        err_hf=err_hf+abs(x_iter(i)-x_iter_(i))**2.d0
     end do
     if(err_hf.lt.hf_conv) exit
     !
  end do
  !
  close(unit_err)
  close(uio)
  !
  call save_array('hf_BLS_final_up.out',x_iter(1:14))
  call save_array('hf_BLS_final_dw.out',x_iter(15:28))
  !stop
  !
  H_Hf=HF_hamiltonian(x_iter)
  H_Hf=H_Hf+Hk_toy
  call fix_mu(H_Hf,delta_hf,mu_fix,eout)

  !+- plot real space hybridizations -+!
  uio=free_unit()
  open(unit=uio,file='hyb_TNS_VS_r_BLS_free.out')
  !allocate(obs(Nso*Nso))
  obs=0.d0
  do ir=1,nrpts
     delta_hfr(:,:,ir)=0.d0
     do ik=1,Lk
        delta_hfr(:,:,ir)=delta_hfr(:,:,ir) + &
             delta_hf(:,:,ik)*exp(xi*dot_product(rpt_latt(ir,:),kpt_latt(ik,:)))*wtk(ik)
     end do
  end do
  do ir=1,nr1d
     Rlat=irvec1d(ir,1)*R1
     iso=0
     do i=1,Norb
        do j=1,Norb
           iso=iso+1
           obs(iso) = delta_hfr(i,j,stride1D(ir))
        end do
     end do
     !
     write(uio,'(100F10.5)') Rlat(1),dreal(obs(:)),dimag(obs(:))
     !
  end do
  close(uio)

  unit_in=free_unit()
  open(unit=unit_in,file='TNS_bands_BLS.out')
  uio=free_unit()
  open(uio,file='TNS_fat_BLS.out')  
  do ir=1,nrpts
     call FT_q2r(rpt_latt(ir,:),Hr_toy(:,:,ir),H_hf)
  end do
  modk=0.d0
  do i=1,3
     delta_kpath=kpath(i+1,:)-kpath(i,:)
     do ik=1,100
        j=(i-1)*100 + ik
        kpt_path(j,:) = kpath(i,:) + dble(ik-1)/100.d0*delta_kpath
        modk=modk+sqrt(dot_product(1.d0/100.d0*delta_kpath,1.d0/100.d0*delta_kpath))
        !
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
              !Ta_fat=Ta_fat+abs(Hktmp(iso,:))**2.d0
           end do
           do iorb=4,5
              iso=(ispin-1)*Norb+iorb
              Ta_fat(:,ispin)=Ta_fat(:,ispin)+abs(Hktmp(iso,(ispin-1)*Norb+1:(ispin-1)*Norb+Norb))**2.d0
              !Ta_fat=Ta_fat+abs(Hktmp(iso,:))**2.d0
           end do
        end do
        Ni_fat=0.d0
        do ispin=1,2
           iorb=3
           iso=(ispin-1)*Norb+iorb
           Ni_fat(:,ispin)=Ni_fat(:,ispin)+abs(Hktmp(iso,(ispin-1)*Norb+1:(ispin-1)*Norb+Norb))**2.d0
           ! Ni_fat=Ni_fat+abs(Hktmp(iso,:))**2.d0
           iorb=6
           iso=(ispin-1)*Norb+iorb
           Ni_fat(:,ispin)=Ni_fat(:,ispin)+abs(Hktmp(iso,(ispin-1)*Norb+1:(ispin-1)*Norb+Norb))**2.d0
           ! Ni_fat=Ni_fat+abs(Hktmp(iso,:))**2.d0           
        end do
        !
        write(unit_in,'(30F18.10)') modk,ek_out_up-mu_fix,ek_out_dw-mu_fix
        write(uio,'(30F18.10)') modk,Ta_fat(:,1),Ta_fat(:,2),Ni_fat(:,1),Ni_fat(:,2)        
        !
     end do
     !
  end do
  close(unit_in)
  close(uio)
  

  
  stop
  
  !
contains
  !
  subroutine get_double_counting_energy(xin,E_dc)
    complex(8),dimension(:),allocatable,intent(in) :: xin
    real(8),intent(out) ::  E_dc
    integer :: ispin,iorb
    if(.not.allocated(xin)) stop "get_ntot"
    if(size(xin).ne.14*Nspin)stop "get_ntot"
    E_dc=0d0
    if(hf_hartree) then
       do i=1,3
          E_dc=E_dc-Ucell*dreal(xin(i))*dreal(xin(i+14))
       end do
       E_dc = E_dc - 2*Vcell*( dreal(xin(1))+dreal(xin(15)) )*( dreal(xin(3))+dreal(xin(17)) )
       E_dc = E_dc - 2*Vcell*( dreal(xin(2))+dreal(xin(16)) )*( dreal(xin(3))+dreal(xin(17)) )
    end if
    E_dc = E_dc + Vcell*(abs(xin(4))**2.d0+abs(xin(5))**2.d0)+Vcell*(abs(xin(18))**2.d0+abs(xin(19))**2.d0)
    E_dc = E_dc + Vcell*(abs(xin(6))**2.d0+abs(xin(7))**2.d0)+Vcell*(abs(xin(20))**2.d0+abs(xin(21))**2.d0)
    !
    if(hf_hartree) then
       do i=8,10
          E_dc=E_dc-Ucell*dreal(xin(i))*dreal(xin(i+14))
       end do
       E_dc = E_dc - 2*Vcell*( dreal(xin(8))+dreal(xin(22)) )*( dreal(xin(10))+dreal(xin(24)) )
       E_dc = E_dc - 2*Vcell*( dreal(xin(9))+dreal(xin(23)) )*( dreal(xin(10))+dreal(xin(24)) )
    end if
    !
    E_dc=E_dc + Vcell*(abs(xin(11))**2.d0+abs(xin(12))**2.d0)+Vcell*(abs(xin(25))**2.d0+abs(xin(26))**2.d0)
    E_dc=E_dc + Vcell*(abs(xin(13))**2.d0+abs(xin(14))**2.d0)+Vcell*(abs(xin(27))**2.d0+abs(xin(28))**2.d0)


    !+- double counting term -+!
    !Eout=Eout-Ucell*0.25d0*(dreal(x_iter(1))**2.d0+dreal(x_iter(2))**2.d0+dreal(x_iter(3))**2.d0)
    !Eout=Eout-2*Vcell*(dreal(x_iter(1))*dreal(x_iter(3))+dreal(x_iter(2))*dreal(x_iter(3)))
    !Eout=Eout + 2.d0*Vcell*(abs(x_iter(4))**2.d0+abs(x_iter(5))**2.d0)
    !Eout=Eout + 2.d0*Vcell*(abs(x_iter(6))**2.d0+abs(x_iter(7))**2.d0)
    !Eout=Eout-Ucell*0.25d0*(dreal(x_iter(8))**2.d0+dreal(x_iter(9))**2.d0+dreal(x_iter(10))**2.d0)
    !Eout=Eout-2*Vcell*(dreal(x_iter(8))*dreal(x_iter(10))+dreal(x_iter(9))*dreal(x_iter(10)))
    !Eout=Eout + 2.d0*Vcell*(abs(x_iter(11))**2.d0+abs(x_iter(12))**2.d0)
    !Eout=Eout + 2.d0*Vcell*(abs(x_iter(13))**2.d0+abs(x_iter(14))**2.d0)
    ! do i=1,3
    !    Eout=Eout-Ucell*dreal(x_iter(i))*dreal(x_iter(i+14))
    ! end do     
    ! Eout = Eout - 2*Vcell*(dreal(x_iter(1)+dreal(x_iter(15))*(dreal(x_iter(3)+dreal(x_iter(17))
    ! Eout = Eout - 2*Vcell*(dreal(x_iter(2)+dreal(x_iter(16))*(dreal(x_iter(3)+dreal(x_iter(17))
    ! Eout = Eout + Vcell*(abs(x_iter(4))**2.d0+abs(x_iter(5))**2.d0)+Vcell*(abs(x_iter(18))**2.d0+abs(x_iter(19))**2.d0))
    ! Eout = Eout + Vcell*(abs(x_iter(6))**2.d0+abs(x_iter(7))**2.d0)+Vcell*(abs(x_iter(20))**2.d0+abs(x_iter(21))**2.d0)
    ! !
    ! do i=8,10
    !    Eout=Eout-Ucell*dreal(x_iter(i))*dreal(x_iter(i+14))
    ! end do
    ! Eout = Eout - 2*Vcell*(dreal(x_iter(8)+dreal(x_iter(22))*(dreal(x_iter(10)+dreal(x_iter(24))
    ! Eout = Eout - 2*Vcell*(dreal(x_iter(9)+dreal(x_iter(23))*(dreal(x_iter(10)+dreal(x_iter(24))
    ! !
    ! Eout=Eout + Vcell*(abs(x_iter(11))**2.d0+abs(x_iter(12))**2.d0)+Vcell*(abs(x_iter(25))**2.d0+abs(x_iter(26))**2.d0)
    ! Eout=Eout + Vcell*(abs(x_iter(13))**2.d0+abs(x_iter(14))**2.d0)+Vcell*(abs(x_iter(27))**2.d0+abs(x_iter(28))**2.d0)
  end subroutine get_double_counting_energy
  
  subroutine enforce_inv_hf(xin)
    complex(8),dimension(:),allocatable,intent(inout) :: xin
    integer :: ispin,iorb
    if(.not.allocated(xin))   stop "enforce_inv_hf"
    if(size(xin).ne.14*Nspin) stop "enforce_inv_hf"
    !+- diagonal terms
    do ispin=1,Nspin
       do iorb=8,10
          xin(iorb+(ispin-1)*14) = xin(iorb-7+(ispin-1)*14)
       end do
    end do
    !+- hybridization
    do ispin=1,Nspin
       xin(5+(ispin-1)*14) = -dreal(xin(6+(ispin-1)*14))-xi*dimag(xin(4+(ispin-1)*14))
       xin(7+(ispin-1)*14) = -dreal(xin(4+(ispin-1)*14))-xi*dimag(xin(6+(ispin-1)*14))       
       xin(11+(ispin-1)*14) = xin(5+(ispin-1)*14)
       xin(12+(ispin-1)*14) = xin(4+(ispin-1)*14)
       xin(13+(ispin-1)*14) = xin(7+(ispin-1)*14)
       xin(14+(ispin-1)*14) = xin(6+(ispin-1)*14)
    end do
    if(enforce_spin) xin(15:28)=xin(1:14)
  end subroutine enforce_inv_hf
  
  subroutine get_ntot(xin,nout)
    complex(8),dimension(:),allocatable,intent(in) :: xin
    real(8),intent(out) ::  nout
    integer :: ispin,iorb
    if(.not.allocated(xin)) stop "get_ntot"
    if(size(xin).ne.14*Nspin)stop "get_ntot"
    nout=0d0
    do ispin=1,Nspin
       do iorb=1,3
          nout=nout+dreal(xin(iorb+(ispin-1)*14))
       end do
       do iorb=8,10
          nout=nout+dreal(xin(iorb+(ispin-1)*14))
       end do       
    end do
  end subroutine get_ntot
  
  subroutine get_xiter(dhfr_in,xout)
    complex(8),dimension(:),allocatable,intent(out) :: xout
    complex(8),dimension(Nso,Nso,nrpts),intent(in) :: dhfr_in
    integer ::  ispin
    if(allocated(xout)) deallocate(xout)
    allocate(xout(14*Nspin));xout=0d0
    do ispin=1,Nspin
       xout(1+(ispin-1)*14) = dhfr_in(1+(ispin-1)*Norb,1+(ispin-1)*Norb,ir0)
       xout(2+(ispin-1)*14) = dhfr_in(2+(ispin-1)*Norb,2+(ispin-1)*Norb,ir0)
       xout(3+(ispin-1)*14) = dhfr_in(3+(ispin-1)*Norb,3+(ispin-1)*Norb,ir0)
       !
       xout(4+(ispin-1)*14) = dhfr_in(1+(ispin-1)*Norb,3+(ispin-1)*Norb,ir0)
       xout(5+(ispin-1)*14) = dhfr_in(2+(ispin-1)*Norb,3+(ispin-1)*Norb,ir0)
       !
       xout(6+(ispin-1)*14) = dhfr_in(1+(ispin-1)*Norb,3+(ispin-1)*Norb,irR)
       xout(7+(ispin-1)*14) = dhfr_in(2+(ispin-1)*Norb,3+(ispin-1)*Norb,irR)
       !
       xout(8+(ispin-1)*14)  = dhfr_in(4+(ispin-1)*Norb,4+(ispin-1)*Norb,ir0)
       xout(9+(ispin-1)*14)  = dhfr_in(5+(ispin-1)*Norb,5+(ispin-1)*Norb,ir0)
       xout(10+(ispin-1)*14) = dhfr_in(6+(ispin-1)*Norb,6+(ispin-1)*Norb,ir0)     
       !
       xout(11+(ispin-1)*14) = dhfr_in(4+(ispin-1)*Norb,6+(ispin-1)*Norb,ir0)
       xout(12+(ispin-1)*14) = dhfr_in(5+(ispin-1)*Norb,6+(ispin-1)*Norb,ir0)
       !
       xout(13+(ispin-1)*14) = dhfr_in(4+(ispin-1)*Norb,6+(ispin-1)*Norb,irL)
       xout(14+(ispin-1)*14) = dhfr_in(5+(ispin-1)*Norb,6+(ispin-1)*Norb,irL)
       !
    end do
  end subroutine get_xiter  
  !
  function HF_hamiltonian(x_iter) result(Hhf) !+da cambiare con il caso generale
    implicit none
    complex(8) :: x_iter(28)
    complex(8),dimension(Nso,Nso,Lk) :: Hhf
    integer :: iso,jso,iorb,jorb
    !
    Hhf=0.d0
    do ik=1,Lk
       do ispin=1,Nspin
          !
          iorb=1
          iso=(ispin-1)*Norb+iorb
          if(hf_hartree) then
             Hhf(iso,iso,ik) = Ucell*x_iter(1+(3-ispin-1)*14)
             Hhf(iso,iso,ik) = Hhf(iso,iso,ik) + 2.d0*Vcell*(x_iter(3)+x_iter(17)) !+2.d0*Vcell*x_iter(3)
          end if
          !
          iorb=2
          iso=(ispin-1)*Norb+iorb
          if(hf_hartree) then
             Hhf(iso,iso,ik) = Ucell*x_iter(2+(3-ispin-1)*14)
             Hhf(iso,iso,ik) = Hhf(iso,iso,ik) + 2.d0*Vcell*(x_iter(3)+x_iter(17)) !+2.d0*Vcell*x_iter(3)
          end if
          !
          iorb=3
          iso=(ispin-1)*Norb+iorb
          if(hf_hartree) then
             Hhf(iso,iso,ik) = Ucell*x_iter(3+(3-ispin-1)*14)
             Hhf(iso,iso,ik) = Hhf(iso,iso,ik) + 2.d0*Vcell*(x_iter(1)+x_iter(2)+x_iter(15)+x_iter(16))!+2.d0*Vcell*(x_iter(2)+x_iter(1))
          end if
          !
          iorb=1; iso=(ispin-1)*Norb+iorb          
          jorb=3; jso=(ispin-1)*Norb+jorb
          !Hhf(jso,iso,ik) = -Vcell*(x_iter(4)+x_iter(6)*exp(-xi*dot_product(R1,kpt_latt(ik,:))))
          Hhf(jso,iso,ik) = -Vcell*(x_iter(4+(ispin-1)*14)+x_iter(6+(ispin-1)*14)*exp(-xi*dot_product(R1,kpt_latt(ik,:))))
          Hhf(iso,jso,ik) = conjg(Hhf(jso,iso,ik))
          !
          iorb=2; iso=(ispin-1)*Norb+iorb          
          jorb=3; jso=(ispin-1)*Norb+jorb
          !Hhf(jso,iso,ik) = -Vcell*(x_iter(5)+x_iter(7)*exp(-xi*dot_product(R1,kpt_latt(ik,:))))
          Hhf(jso,iso,ik) = -Vcell*(x_iter(5+(ispin-1)*14)+x_iter(7+(ispin-1)*14)*exp(-xi*dot_product(R1,kpt_latt(ik,:))))
          Hhf(iso,jso,ik) = conjg(Hhf(jso,iso,ik))
          !
          !
          iorb=4
          iso=(ispin-1)*Norb+iorb
          if(hf_hartree) then
             Hhf(iso,iso,ik) = Ucell*x_iter(8+(3-ispin-1)*14)
             Hhf(iso,iso,ik) = Hhf(iso,iso,ik) + 2.d0*Vcell*(x_iter(10)+x_iter(24)) !+2.d0*Vcell*x_iter(10)
          end if
          !
          iorb=5
          iso=(ispin-1)*Norb+iorb
          if(hf_hartree) then
             Hhf(iso,iso,ik) = Ucell*x_iter(9+(3-ispin-1)*14)
             Hhf(iso,iso,ik) = Hhf(iso,iso,ik) + 2.d0*Vcell*(x_iter(10)+x_iter(24)) !+ 2.d0*Vcell*x_iter(10)
          end if
          !
          iorb=6
          iso=(ispin-1)*Norb+iorb
          if(hf_hartree) then
             Hhf(iso,iso,ik) = Ucell*x_iter(10+(3-ispin-1)*14) 
             Hhf(iso,iso,ik) = Hhf(iso,iso,ik) + 2.d0*Vcell*(x_iter(8)+x_iter(9)+x_iter(22)+x_iter(23)) !+ 2.d0*Vcell*(x_iter(9)+x_iter(8))
          end if
          !
          iorb=4; iso=(ispin-1)*Norb+iorb          
          jorb=6; jso=(ispin-1)*Norb+jorb
          !Hhf(jso,iso,ik) = -Vcell*(x_iter(11)+x_iter(13)*exp(xi*dot_product(R1,kpt_latt(ik,:))))
          Hhf(jso,iso,ik) = -Vcell*(x_iter(11+(ispin-1)*14)+x_iter(13+(ispin-1)*14)*exp(xi*dot_product(R1,kpt_latt(ik,:))))
          Hhf(iso,jso,ik) = conjg(Hhf(jso,iso,ik))
          !
          iorb=5; iso=(ispin-1)*Norb+iorb          
          jorb=6; jso=(ispin-1)*Norb+jorb
          !Hhf(jso,iso,ik) = -Vcell*(x_iter(12)+x_iter(14)*exp(xi*dot_product(R1,kpt_latt(ik,:))))
          Hhf(jso,iso,ik) = -Vcell*(x_iter(12+(ispin-1)*14)+x_iter(14+(ispin-1)*14)*exp(xi*dot_product(R1,kpt_latt(ik,:))))
          Hhf(iso,jso,ik) = conjg(Hhf(jso,iso,ik))
          !
       end do
    end do
    !
  end function HF_hamiltonian
  


  function HF_hamiltonian_ei(x_iter) result(Hhf)
    implicit none
    complex(8) :: x_iter(11)
    complex(8),dimension(Nso,Nso,Lk) :: Hhf
    integer :: iso,jso,iorb,jorb

    !
    Hhf=0.d0
    do ik=1,Lk
       do ispin=1,Nspin
          !
          iorb=1
          iso=(ispin-1)*Norb+iorb
          Hhf(iso,iso,ik) = 0.5d0*Ucell*x_iter(1)+2.d0*Vcell*x_iter(3)
          !
          iorb=2
          iso=(ispin-1)*Norb+iorb
          Hhf(iso,iso,ik) = 0.5d0*Ucell*x_iter(2)+2.d0*Vcell*x_iter(3)          
          !
          iorb=3
          iso=(ispin-1)*Norb+iorb
          Hhf(iso,iso,ik) = 0.5d0*Ucell*x_iter(3)+2.d0*Vcell*(x_iter(2)+x_iter(1))
          !
          
          !
          iorb=1; iso=(ispin-1)*Norb+iorb          
          jorb=3; jso=(ispin-1)*Norb+jorb
          Hhf(jso,iso,ik) = -Vcell*x_iter(4)*(1.d0+exp(-xi*dot_product(R1,kpt_latt(ik,:))))
          ! Hhf(jso,iso,ik) = Hhf(jso,iso,ik) - Vcell*x_iter(10)*exp(-xi*dot_product(R1,kpt_latt(ik,:)))
          ! Hhf(iso,jso,ik) = conjg(Hhf(jso,iso,ik))
          !
          iorb=2; iso=(ispin-1)*Norb+iorb          
          jorb=3; jso=(ispin-1)*Norb+jorb
          Hhf(jso,iso,ik) = -Vcell*x_iter(5)*(1.d0+exp(-xi*dot_product(R1,kpt_latt(ik,:))))
          ! Hhf(jso,iso,ik) = Hhf(jso,iso,ik) - Vcell*x_iter(11)*exp(-xi*dot_product(R1,kpt_latt(ik,:)))
          ! Hhf(iso,jso,ik) = conjg(Hhf(jso,iso,ik))
          !
       end do
    end do
    !
  end function HF_hamiltonian_ei






  subroutine get_Hk_w90(kpoint,Hk) 
    implicit none
    real(8),dimension(3) :: kpoint
    complex(8),dimension(Nso,Nso) :: Hk
    integer :: ir,nrpts
    real(8),dimension(3) :: Rlat
    real(8) :: dotRk
    !
    nrpts=size(irvec,1)    
    Hk=0.d0
    do ir=1,nrpts
       Rlat=irvec(ir,1)*R1+irvec(ir,2)*R2+irvec(ir,3)*R3
       dotRk=dot_product(Rlat,kpoint)
       Hk=Hk+Hr_w90(:,:,ir)*exp(xi*dotRK)/dble(ndegen(ir))
    end do
    !
  end subroutine get_Hk_w90


  ! function get_HK(kpoint,N,Hk_in) result(Hk_out)
  !   implicit none
  !   real(8),dimension(:) :: kpoint
  !   integer              :: N
  !   complex(8),dimension(N,N) :: Hk_out
  !   complex(8),dimension(N,N,Nk_x,Nk_y,Nk_z) :: Hk_in
  !   integer :: io,jo,Nint_,ik_target
  !   complex(8),dimension(:,:,:,:,:),allocatable :: Hk_inter
  !   real(8),dimension(:,:,:),allocatable :: Hkr
  !   real(8),dimension(:),allocatable :: kix,kiy,kiz
  !   real(8) :: Hr,dH,delta_k,diffk
  !   integer :: i0,j0,k0
  !   integer :: i,j,k,ii,jj,kk
  !   integer :: ix,ik
  !   real(8),dimension(:,:),allocatable :: kpt_int
  !   real(8),dimension(3) :: ktarget,ktmp
  !   integer :: Lkint,ikint
  !   real(8),dimension(:),allocatable :: yout,yout_,hk_tmp
  !   real(8),dimension(:),allocatable :: grid_aux

  !   Nint_=2*Nint+1
  !   allocate(Hk_inter(Nso,Nso,Nint_,Nint_,Nint_))
  !   allocate(kix(Nint_),kiy(Nint_),kiz(Nint_))


  !   if(size(kpoint).ne.3) stop "kpoint.ne.3"    
  !   if(N.ne.Nso)  stop "N/=Nso"

  !   !+- find the closest point to the target one -+!
  !   diffk=100.d0
  !   ik_target=1
  !   do ik=1,Lk
  !      delta_k=0.d0
  !      delta_k=delta_k+(kpoint(1)-kpt_latt(ik,1))**2.d0
  !      delta_k=delta_k+(kpoint(2)-kpt_latt(ik,2))**2.d0
  !      delta_k=delta_k+(kpoint(3)-kpt_latt(ik,3))**2.d0
  !      if(delta_k.le.diffk) then
  !         ik_target=ik
  !         diffk=delta_k
  !      end if
  !   end do
  !   !
  !   ktarget=matmul(Bkinv,kpoint)
  !   !
  !   !+- get the *centered* index of the aux lattice -+!
  !   i0=ik_stride(ik_target,1)+Nk_x
  !   j0=ik_stride(ik_target,2)+Nk_y
  !   k0=ik_stride(ik_target,3)+Nk_z
  !   !
  !   !
  !   Lkint=Nint_*Nint_*Nint_
  !   allocate(kpt_int(Lkint,3))

  !   ikint=0
  !   do i=1,Nint_       
  !      ii=i0-(Nint+1)+i
  !      kix(i)= kxgrid_aux(ii)
  !      do j=1,Nint_
  !         jj=j0-(Nint+1)+j
  !         kiy(j)= kygrid_aux(jj)
  !         do k=1,Nint_
  !            kk=k0-(Nint+1)+k
  !            kiz(k)=kzgrid_aux(kk)
  !            ikint=ikint+1
  !            !
  !            kpt_int(ikint,:)=kix(i)*Bk1+ kiy(j)*Bk2 + kiz(k)*Bk3
  !            !
  !            Hk_inter(:,:,i,j,k)=Hk_in(:,:,ix_aux(ii),iy_aux(jj),iz_aux(kk)) 
  !            !
  !         end do
  !      end do
  !   end do
  !   !

  !   allocate(Hkr(Nint_,Nint_,Nint_))    

  !   Hk_out=0.d0
  !   do io=1,Nso
  !      do jo=1,Nso

  !         !+-get Hkr to be interpolated 
  !         Hkr=dreal(Hk_inter(io,jo,:,:,:))
  !         Hr=0.d0
  !         call polin3(kix,kiy,kiz,Hkr,ktarget(1),ktarget(2),ktarget(3),Hr,dH)          
  !         Hk_out(io,jo)=Hk_out(io,jo) + Hr
  !         !
  !         Hkr=dimag(Hk_inter(io,jo,:,:,:))
  !         Hr=0.d0
  !         call polin3(kix,kiy,kiz,Hkr,ktarget(1),ktarget(2),ktarget(3),Hr,dH)          
  !         Hk_out(io,jo)=Hk_out(io,jo) + xi*Hr          
  !         !
  !      end do
  !   end do
  !   !
  !   write(*,'(20F10.5)') kpoint,kpt_latt(ik_target,:),ktarget,Hk_w90(1,1,ik_target),Hk_out(1,1)
  !   !write(*,*) 'test_inter polin3',Hk_out(1,1)
  !   !write(*,*) Hk_out(1,1),Hk_out(2,2)
  ! end function Get_HK
  
  ! subroutine build_mp_grid_1d(Nx,Ny,Nz)
  !   implicit none
  !   integer :: Nx,Ny,Nz
  !   integer :: Nx_,Ny_,Nz_
  !   integer :: i,ix,iy,ik,iz,unitk
  !   real(8) :: kmax,kxmax,kymax,kzmax,dx,dy,dz
  !   real(8),dimension(:),allocatable :: kxr,kyr,kzr
  !   !
  !   !    
  !   if(allocated(kxgrid)) deallocate(kxgrid)
  !   allocate(kxgrid(Nx))
  !   kmax=dble(Nx-1)*0.5d0/dble(Nx)
  !   kxgrid = linspace(-kmax,kmax,Nx,mesh=dx)
  !   kxmax=kmax+dx*0.5d0
  !   !
  !   if(mod(Nx,2).eq.0) then
  !      kxgrid=kxgrid-0.5d0/dble(Nx)
  !   end if
  !   !
  !   unitk=free_unit()
  !   open(unit=unitk,file='mp_grid.out')
  !   do ix=1,Nx
  !      do iy=1,Ny
  !         do iz=1,Nz
  !            write(unitk,'(10F18.10)') kxgrid(ix)
  !         end do
  !      end do
  !   end do
  !   close(unitk)
  ! end subroutine build_mp_grid_1d


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
    kygrid = linspace(-kmax,kmax,Ny,mesh=dy)
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

