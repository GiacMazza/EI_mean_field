program officina
  !
  USE SCIFOR
  USE DMFT_TOOLS
  !
  USE VARS_GLOBAL
  USE HF_real
  !
  USE RK_IDE
  USE MPI
  !
  !
  implicit none
  integer :: Nint

  real(8) :: KKx,KKy
  complex(8),dimension(:,:,:),allocatable :: delta_hf,H_hf,delta_hf_,Hhf_grid,delta_hfr,delta_hfr_,H_hf_dyn
  complex(8),dimension(:,:,:),allocatable :: Hr_w90,Hr_w90_tmp,Hr_toy  
  real(8) :: mu_fix
  real(8) :: Uintra,Uinter,Delta_CF
  real(8) :: Uq,thop,err_hf,hf_conv
  logical :: plot_w90_bands
  !
  real(8),dimension(3) :: ktest
  real(8),dimension(:),allocatable :: ek_out

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

  integer :: iso,jso,ispin,jspin,iorb,jorb,ik,isys
  integer :: i,j,k,idim
  integer :: Nhf,Nhf_,ihf,unit_err,unit_obs,unit_in,uio,unit_io,Nobs,jhf
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

  real(8) :: Eout,EoutHF,EoutLgr




  real(8) :: modk,tmpi,tmpj,checkR
  real(8) :: kmod(3)
  real(8) :: Ncell
  complex(8),dimension(:,:),allocatable :: ccij
  logical :: hartree
  real(8),allocatable,dimension(:) :: Ta_fat,Ni_fat
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
  complex(8),dimension(14) :: x_iter,x_iter_
  complex(8),dimension(14) :: xtmp,xtmp_
  complex(8) :: xphi,xphi_,xphi_vec(2)
  complex(8),dimension(2) :: xpi,xpi_
  complex(8),dimension(4) :: op_TNS
  real(8),dimension(14) :: xr_iter
  real(8),dimension(28) :: xr_tmp

  real(8) :: Ucell,Vcell,Wcell
  real(8) :: Evalence,Econduction,tconduction,tvalence,tt_hyb,nn_hyb,tn_hyb
  real(8) :: w0gap
  logical :: fix_phi  !+-> HF calculation with fixed order parameter <-+!
  logical :: hf_symm  !+-> HF calculation w/ symmetric order parameter <-+!
  logical :: HF_solve !+-> HF calculation with a root-finder routine  <-+
  logical :: hf_in
  integer :: lgr_fix_verbose
  logical :: lgr_fix_phase

  real(8) :: phi_start,phi_end,phi_phase,dphi,ntot,test,err
  real(8) :: phase_start,phase_end,dphase
  real(8),dimension(:),allocatable :: mod_phi,varphi
  integer :: Nphi,Nphi_phase
  integer :: ixr(3)
  !
  real(8),dimension(4) :: phi_sgn
  complex(8) :: lgr_iter(2)
  real(8),dimension(2) :: lgr_iter_tmp
  real(8),dimension(1) :: lgr_iter_tmp_

  complex(8),dimension(:),allocatable :: phi_list
  !

  !+- the dynamics
  complex(8),dimension(:),allocatable :: psit
  integer :: Nt_aux,Ndyn,it,Nt
  real(8) :: tstart,tstop,tstep,time

  integer,dimension(:,:),allocatable :: ivec2idelta
  integer,dimension(:,:,:),allocatable :: idelta2ivec
  real(8),dimension(:),allocatable :: t_grid,t_grid_aux
  
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
  ! call parse_input_variable(Nint,"Nint","input.conf",default=1)
  ! call parse_input_variable(hartree,"Hartree","input.conf",default=.true.)  

  ! call parse_input_variable(Ucut_off,"U_CUT","input.conf",default=20.d0)
  ! call parse_input_variable(H1d,"H1D","input.conf",default=.false.)
  ! call parse_input_variable(alphaU,"alphaU","input.conf",default=1.d0)

  call parse_input_variable(Vcell,"V","input.conf",default=1.d0)
  call parse_input_variable(Ucell,"U","input.conf",default=1.d0)
  call parse_input_variable(Wcell,"W","input.conf",default=0.d0)


  call parse_input_variable(Econduction,"Econduction","input.conf",default=1.7d0)
  call parse_input_variable(Evalence,"Evalence","input.conf",default=-0.9d0)
  call parse_input_variable(tconduction,"tconduction","input.conf",default=-0.8d0)
  call parse_input_variable(tvalence,"tvalence","input.conf",default=0.4d0)

  call parse_input_variable(tt_hyb,"tt_hyb","input.conf",default=0.0d0)
  call parse_input_variable(nn_hyb,"nn_hyb","input.conf",default=0.0d0)
  call parse_input_variable(tn_hyb,"tn_hyb","input.conf",default=0.0d0)

  
  call parse_input_variable(w0gap,"w0gap","input.conf",default=0.0d0)
  call parse_input_variable(hybloc,"hybloc","input.conf",default=0.d0)
  call parse_input_variable(fix_phi,"fix_phi","input.conf",default=.false.)
  call parse_input_variable(phi_start,"phi_start","input.conf",default=0.01d0)
  call parse_input_variable(phi_end,"phi_end","input.conf",default=0.2d0)

  call parse_input_variable(phase_start,"phase_start","input.conf",default=0.01d0)
  call parse_input_variable(phase_end,"phase_end","input.conf",default=0.2d0)

  
  call parse_input_variable(phi_phase,"phi_phase","input.conf",default=0d0)

  call parse_input_variable(Nphi,"Nphi","input.conf",default=1)
  call parse_input_variable(Nphi_phase,"Nphi_phase","input.conf",default=1)

  call parse_input_variable(phi_sgn,"phi_sign","input.conf",default=[1.d0,-1.d0,-1.d0,1.d0])
  
  
  call parse_input_variable(hf_symm,"hf_symm","input.conf",default=.false.)
  call parse_input_variable(hf_solve,"hf_solve","input.conf",default=.false.)

  call parse_input_variable(hf_conv,"hf_conv","input.conf",default=1.d-10)
  call parse_input_variable(lgr_fix_verbose,"lgr_fix_verbose","input.conf",default=0)
  call parse_input_variable(lgr_fix_phase,"lgr_fix_phase","input.conf",default=.false.)

  call parse_input_variable(tstart,"TSTART","input.conf",default=0.d0)
  call parse_input_variable(tstop,"TSTOP","input.conf",default=0.1d0)
  call parse_input_variable(tstep,"TSTEP","input.conf",default=1.d-3)
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
        Hk_toy(iso,iso,ik) = Econduction  + 2.d0*tk(iorb)*dcos(kpt_latt(ik,1)*R1(1)+kpt_latt(ik,2)*R1(2)+kpt_latt(ik,3)*R1(3))
        iorb=2
        iso=(ispin-1)*Norb+iorb
        Hk_toy(iso,iso,ik) = Econduction  + 2.d0*tk(iorb)*dcos(kpt_latt(ik,1)*R1(1)+kpt_latt(ik,2)*R1(2)+kpt_latt(ik,3)*R1(3))
        !
        iorb=3
        iso=(ispin-1)*Norb+iorb
        Hk_toy(iso,iso,ik) = Evalence  + 2.d0*tk(iorb)*dcos(kpt_latt(ik,1)*R1(1)+kpt_latt(ik,2)*R1(2)+kpt_latt(ik,3)*R1(3))
        !
        iorb=1
        jorb=3
        iso=(ispin-1)*Norb+iorb
        jso=(ispin-1)*Norb+jorb
        Rlat=R1
        Hk_toy(iso,jso,ik) = hybloc*(1.d0-exp(xi*dot_product(Rlat,kpt_latt(ik,:)))) 
        Hk_toy(jso,iso,ik) = hybloc*(1.d0-exp(-xi*dot_product(Rlat,kpt_latt(ik,:))))
        !
        iorb=2
        jorb=3
        iso=(ispin-1)*Norb+iorb
        jso=(ispin-1)*Norb+jorb
        Rlat=R1
        Hk_toy(iso,jso,ik) = hybloc*(1.d0-exp(xi*dot_product(Rlat,kpt_latt(ik,:)))) 
        Hk_toy(jso,iso,ik) = hybloc*(1.d0-exp(-xi*dot_product(Rlat,kpt_latt(ik,:))))
        !
        !+- lower chain -+!
        !
        iorb=4
        iso=(ispin-1)*Norb+iorb
        Hk_toy(iso,iso,ik) = Econduction  + 2.d0*tk(iorb)*dcos(kpt_latt(ik,1)*R1(1)+kpt_latt(ik,2)*R1(2)+kpt_latt(ik,3)*R1(3))
        iorb=5
        iso=(ispin-1)*Norb+iorb
        Hk_toy(iso,iso,ik) = Econduction  + 2.d0*tk(iorb)*dcos(kpt_latt(ik,1)*R1(1)+kpt_latt(ik,2)*R1(2)+kpt_latt(ik,3)*R1(3))
        !
        iorb=6
        iso=(ispin-1)*Norb+iorb
        Hk_toy(iso,iso,ik) = Evalence  + 2.d0*tk(iorb)*dcos(kpt_latt(ik,1)*R1(1)+kpt_latt(ik,2)*R1(2)+kpt_latt(ik,3)*R1(3))
        !
        iorb=4
        jorb=6
        iso=(ispin-1)*Norb+iorb
        jso=(ispin-1)*Norb+jorb
        Rlat=-R1
        Hk_toy(iso,jso,ik) = hybloc*(1.d0-exp(xi*dot_product(Rlat,kpt_latt(ik,:))))
        Hk_toy(jso,iso,ik) = conjg(Hk_toy(iso,jso,ik))
        !
        iorb=5
        jorb=6
        iso=(ispin-1)*Norb+iorb
        jso=(ispin-1)*Norb+jorb
        Rlat=-R1
        Hk_toy(iso,jso,ik) = hybloc*(1.d0-exp(xi*dot_product(Rlat,kpt_latt(ik,:))))
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
  !
  close(uio)
  uio=free_unit()
  open(unit=uio,file='hyb_bareR.out')
  allocate(obs(Norb*Norb))
  !
  modk=0.d0
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
     write(unit_in,'(100F10.5)') Rlat(1),dreal(obs(:)),dreal(obs(:))
     !
  end do
  ! 
  !+- plot bare bands -+!
  allocate(kpt_path(400,3))
  allocate(ek_out(Nso))
  allocate(Hktmp(Nso,Nso))
  !
  uio=free_unit()
  open(unit=uio,file='tns_bare_bands.out')
  unit_in=free_unit()
  open(unit_in,file='tns_fat_bare.tns')
  allocate(Ta_fat(Nso),Ni_fat(Nso))
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
        ek_out=0.0d0
        call eigh(Hktmp,ek_out)
        !
        Ta_fat=0.d0
        do ispin=1,2
           do iorb=1,2
              iso=(ispin-1)*Norb+iorb
              Ta_fat=Ta_fat+abs(Hktmp(iso,:))**2.d0
           end do
           do iorb=4,5
              iso=(ispin-1)*Norb+iorb
              Ta_fat=Ta_fat+abs(Hktmp(iso,:))**2.d0
           end do
        end do
        Ni_fat=0.d0
        do ispin=1,2
           iorb=3
           iso=(ispin-1)*Norb+iorb
           Ni_fat=Ni_fat+abs(Hktmp(iso,:))**2.d0
           iorb=6
           iso=(ispin-1)*Norb+iorb
           Ni_fat=Ni_fat+abs(Hktmp(iso,:))**2.d0           
        end do
        !
        write(uio,'(30F18.10)') modk,ek_out-mu_fix
        write(unit_in,'(30F18.10)') modk,Ta_fat,Ni_fat        
        !
     end do
     !
  end do
  close(uio)
  close(unit_in)
  !  
  x_iter(1) = delta_hfr(1,1,ir0)+delta_hfr(1+Norb,1+Norb,ir0)
  x_iter(2) = delta_hfr(2,2,ir0)+delta_hfr(2+Norb,2+Norb,ir0)
  x_iter(3) = delta_hfr(3,3,ir0)+delta_hfr(3+Norb,3+Norb,ir0)
  !
  x_iter(4) = delta_hfr(1,3,ir0)
  x_iter(5) = delta_hfr(2,3,ir0)
  !
  x_iter(6) = dreal(delta_hfr(1,3,ir0)+delta_hfr(1,3,irR))
  x_iter(7) = dreal(delta_hfr(2,3,ir0)+delta_hfr(2,3,irR))
  !
  x_iter(8) = delta_hfr(4,4,ir0)+delta_hfr(4+Norb,4+Norb,ir0)
  x_iter(9) = delta_hfr(5,5,ir0)+delta_hfr(5+Norb,5+Norb,ir0)
  x_iter(10) = delta_hfr(6,6,ir0)+delta_hfr(6+Norb,6+Norb,ir0)
  !
  x_iter(11) = delta_hfr(4,6,ir0)
  x_iter(12) = delta_hfr(5,6,ir0)
  !
  x_iter(13) = dreal(delta_hfr(4,6,ir0)+delta_hfr(4,6,irL))
  x_iter(14) = dreal(delta_hfr(5,6,ir0)+delta_hfr(5,6,irL))
  !
  write(*,*) "HF optimization: bare parameters"
  do iso=1,14
     write(*,*) x_iter(iso)
  end do


  ntot=dreal(x_iter(1))+dreal(x_iter(2))+dreal(x_iter(3))
  ntot=ntot+dreal(x_iter(8))+dreal(x_iter(9))+dreal(x_iter(10))

  uio=free_unit()
  open(unit=uio,file='bare_energy.out')
  write(uio, '(20F18.10)') dreal(x_iter(1:14)),Eout+mu_fix*ntot,Eout
  close(uio)
  !
  ! x_iter(6)=0.1d0
  ! x_iter(7)=0.1d0
  ! x_iter(13)=0.1d0
  ! x_iter(14)=0.1d0
  !
  allocate(H_Hf(Nso,Nso,Lk))

  ! x_iter_=x_iter
  ! allocate(varphi(50)); varphi=linspace(0d0,0.5d0,50)
  ! do j=1,50
  !    x_iter_(6)=0.1d0*exp(xi*varphi(j))*phi_sgn(1)
  !    x_iter_(7)=0.1d0*exp(xi*varphi(j))*phi_sgn(2)
  !    x_iter_(13)=0.1d0*exp(xi*varphi(j))*phi_sgn(3)
  !    x_iter_(14)=0.1d0*exp(xi*varphi(j))*phi_sgn(4)
  !    !
  !    H_Hf=HF_hamiltonian(x_iter_)
  !    H_Hf=H_Hf+Hk_toy
  !    !
  !    call fix_mu(H_Hf,delta_hf,mu_fix,eout)


     
  !    do i=1,3
  !       ir=ixr(i)
  !       delta_hfr(:,:,ir)=0.d0
  !       do ik=1,Lk
  !          delta_hfr(:,:,ir)=delta_hfr(:,:,ir) + &
  !               delta_hf(:,:,ik)*exp(xi*dot_product(rpt_latt(ir,:),kpt_latt(ik,:)))*wtk(ik)
  !       end do
  !    end do
  !    !
  !    x_iter(1) = delta_hfr(1,1,ir0)+delta_hfr(1+Norb,1+Norb,ir0)
  !    x_iter(2) = delta_hfr(2,2,ir0)+delta_hfr(2+Norb,2+Norb,ir0)
  !    x_iter(3) = delta_hfr(3,3,ir0)+delta_hfr(3+Norb,3+Norb,ir0)
  !    !
  !    x_iter(4) = delta_hfr(1,3,ir0)
  !    x_iter(5) = delta_hfr(2,3,ir0)
  !    !
  !    x_iter(6) = delta_hfr(1,3,ir0)+delta_hfr(1,3,irR)
  !    x_iter(7) = delta_hfr(2,3,ir0)+delta_hfr(2,3,irR)
  !    !
  !    x_iter(8) = delta_hfr(4,4,ir0)+delta_hfr(4+Norb,4+Norb,ir0)
  !    x_iter(9) = delta_hfr(5,5,ir0)+delta_hfr(5+Norb,5+Norb,ir0)
  !    x_iter(10) = delta_hfr(6,6,ir0)+delta_hfr(6+Norb,6+Norb,ir0)
  !    !
  !    x_iter(11) = delta_hfr(4,6,ir0)
  !    x_iter(12) = delta_hfr(5,6,ir0)
  !    !
  !    x_iter(13) = delta_hfr(4,6,ir0)+delta_hfr(4,6,irL)
  !    x_iter(14) = delta_hfr(5,6,ir0)+delta_hfr(5,6,irL)
  !    !
  !    write(700,'(20F18.10)') varphi(j),eout+mu_fix*4d0,x_iter(6),x_iter(7),x_iter(13),x_iter(14)
  ! end do
  ! stop
  !
  !
  write(*,*) "HF optimization @ fixed phi: initial paramteres"
  do ihf=1,size(x_iter)
     write(*,*) x_iter(ihf)
  end do
  
  !+- create the list of OP -+!
  allocate(phi_list(Nphi*Nphi_phase));phi_list=0d0
  if(allocated(varphi)) deallocate(varphi); allocate(varphi(Nphi_phase));
  varphi=phase_start
  if(Nphi_phase>1) varphi=linspace(phase_start,phase_end,Nphi_phase)
  
  allocate(mod_phi(Nphi));mod_phi=0d0
  uio=free_unit()
  open(unit=uio,file='list_phi.out')     
  ihf=0
  do i=1,Nphi_phase        
     if((-1d0)**i.lt.0d0) then
        mod_phi=phi_start
        if(Nphi>1) mod_phi=linspace(phi_start,phi_end,Nphi)
     else
        mod_phi=phi_start
        if(Nphi>1) mod_phi=linspace(phi_end,phi_start,Nphi)
     end if
     do j=1,Nphi
        ihf=ihf+1
        phi_list(ihf) = mod_phi(j)*exp(xi*varphi(i)*pi)
        write(uio,'(5F18.10)') phi_list(ihf)
     end do
  end do
  close(uio)
  !
  !
  uio=free_unit()
  open(unit=uio,file='loop_fixed_order_parameter.out')
  close(uio)
  !
  !
  unit_err=free_unit()
  open(unit=unit_err,file='err_symm.out')
  close(unit_err)
  !
  unit_io=free_unit()
  open(unit=unit_io,file='ENE_VS_OP.out')
  close(unit_io)
  !
  dphi = (phi_end-phi_start)/dble(Nphi)
  xphi=(phi_start-dphi)*exp(xi*phi_phase*pi)
  do ihf=1,size(phi_list)
     !
     xphi=phi_list(ihf)
     write(*,*) 'fixed phi loop',ihf,xphi
     !
     lgr_iter=0d0
     !+-
     err_hf=1.d0
     do jhf=1,Nhf
        
        unit_err=free_unit()
        open(unit=unit_err,file='err_symm.out',status='old',position='append')
        write(unit_err,*) err_hf           
        close(unit_err)
        !
        x_iter_=x_iter
        !
        !+- here I should fix the lgr params
        if(lgr_fix_phase) then        
           lgr_iter_tmp_(1)=abs(lgr_iter(1))
           call fsolve(fix_lgr_params_,lgr_iter_tmp_,tol=1.d-8)                   
           lgr_iter=lgr_iter_tmp_(1)*dreal(xphi)/abs(xphi)+xi*lgr_iter_tmp_(1)*dimag(xphi)/abs(xphi)
        else
           lgr_iter_tmp(1)=dreal(lgr_iter(1))
           lgr_iter_tmp(2)=dimag(lgr_iter(1))
           call fsolve(fix_lgr_params,lgr_iter_tmp,tol=1.d-8)
           lgr_iter=lgr_iter_tmp(1)+xi*lgr_iter_tmp(2)
        end if
        !
        H_Hf=HF_hamiltonian(x_iter,lgr_iter)
        H_Hf=H_Hf+Hk_toy
        !
        call fix_mu(H_Hf,delta_hf,mu_fix,eout)
        eoutHF=eout        
        !
        do i=1,3
           ir=ixr(i)
           delta_hfr(:,:,ir)=0.d0
           do ik=1,Lk
              delta_hfr(:,:,ir)=delta_hfr(:,:,ir) + &
                   delta_hf(:,:,ik)*exp(xi*dot_product(rpt_latt(ir,:),kpt_latt(ik,:)))*wtk(ik)
           end do
        end do
        !
        x_iter(1) = delta_hfr(1,1,ir0)+delta_hfr(1+Norb,1+Norb,ir0)
        x_iter(2) = delta_hfr(2,2,ir0)+delta_hfr(2+Norb,2+Norb,ir0)
        x_iter(3) = delta_hfr(3,3,ir0)+delta_hfr(3+Norb,3+Norb,ir0)
        !
        x_iter(4) = delta_hfr(1,3,ir0)
        x_iter(5) = delta_hfr(2,3,ir0)
        !
        x_iter(6) = delta_hfr(1,3,ir0)+delta_hfr(1,3,irR)
        x_iter(7) = delta_hfr(2,3,ir0)+delta_hfr(2,3,irR)
        !
        x_iter(8) = delta_hfr(4,4,ir0)+delta_hfr(4+Norb,4+Norb,ir0)
        x_iter(9) = delta_hfr(5,5,ir0)+delta_hfr(5+Norb,5+Norb,ir0)
        x_iter(10) = delta_hfr(6,6,ir0)+delta_hfr(6+Norb,6+Norb,ir0)
        !
        x_iter(11) = delta_hfr(4,6,ir0)
        x_iter(12) = delta_hfr(5,6,ir0)
        !
        x_iter(13) = delta_hfr(4,6,ir0)+delta_hfr(4,6,irL)
        x_iter(14) = delta_hfr(5,6,ir0)+delta_hfr(5,6,irL)
        !
        !
        !
        ntot=dreal(x_iter(1))+dreal(x_iter(2))+dreal(x_iter(3))
        ntot=ntot+dreal(x_iter(8))+dreal(x_iter(9))+dreal(x_iter(10))
        !
        !
        !
        Eout=Eout-Ucell*0.25d0*(dreal(x_iter(1))**2.d0+dreal(x_iter(2))**2.d0+dreal(x_iter(3))**2.d0)
        Eout=Eout-2d0*Vcell*(dreal(x_iter(1))*dreal(x_iter(3))+dreal(x_iter(2))*dreal(x_iter(3)))
        Eout=Eout + 2.d0*Vcell*(abs(x_iter(4))**2.d0+abs(x_iter(5))**2.d0)
        Eout=Eout + 2.d0*Vcell*(abs(x_iter(6)-x_iter(4))**2.d0+abs(x_iter(7)-x_iter(5))**2.d0)
        !
        Eout=Eout-Ucell*0.25d0*(dreal(x_iter(8))**2.d0+dreal(x_iter(9))**2.d0+dreal(x_iter(10))**2.d0)
        Eout=Eout-2d0*Vcell*(dreal(x_iter(8))*dreal(x_iter(10))+dreal(x_iter(9))*dreal(x_iter(10)))
        !
        Eout=Eout + 2.d0*Vcell*(abs(x_iter(11))**2.d0+abs(x_iter(12))**2.d0)
        Eout=Eout + 2.d0*Vcell*(abs(x_iter(13)-x_iter(11))**2.d0+abs(x_iter(14)-x_iter(12))**2.d0)
        !
        !+- lagrange parameter term -+!
        EoutLgr = 0d0
        EoutLgr = EoutLgr-4d0*dreal(phi_sgn(1)*x_iter(6)*lgr_iter(1)) 
        EoutLgr = EoutLgr-4d0*dreal(phi_sgn(2)*x_iter(7)*lgr_iter(2)) 
        EoutLgr = EoutLgr-4d0*dreal(phi_sgn(3)*x_iter(13)*lgr_iter(2)) 
        EoutLgr = EoutLgr-4d0*dreal(phi_sgn(4)*x_iter(14)*lgr_iter(1))        
        Eout=Eout+EoutLgr
        ! write(546,'(10F18.10)') phi_sgn(1)*x_iter(6)*lgr_iter,phi_sgn(2)*x_iter(7)*lgr_iter
        ! write(547,'(10F18.10)') phi_sgn(3)*x_iter(13)*lgr_iter,phi_sgn(4)*x_iter(14)*lgr_iter
        !

        !
        uio=free_unit()
        open(unit=uio,file='loop_fixed_order_parameter.out',status='old',position='append')
        write(uio, '(40F18.10)') x_iter(1:14),Eout+mu_fix*ntot,EoutHF+mu_fix*ntot,EoutLgr,dreal(lgr_iter),dimag(lgr_iter)
        close(uio)        
        !
        err_hf=0.d0
        do i=1,14
           err_hf=err_hf+abs(x_iter(i)-x_iter_(i))**2.d0
        end do
        if(err_hf.lt.hf_conv) exit
        x_iter=x_iter*wmix+(1.d0-wmix)*x_iter_
     end do
     open(unit=uio,file='loop_fixed_order_parameter.out',status='old',position='append')
     write(uio, '(30F18.10)')
     write(uio, '(30F18.10)')
     close(uio)
     !
     unit_io=free_unit()
     open(unit=unit_io,file='ENE_VS_OP.out',status='old',position='append')
     write(unit_io, '(10F18.10)') xphi,Eout+mu_fix*ntot
     close(unit_io)
     !        
     !
  end do
  close(uio)
  close(unit_in)     
  !
  if(fix_phi) stop
  
  !+- the dynamics -+!  
  Nt=nint((tstop-tstart)/tstep)
  if(master) write(*,*) 'number of time steps Nt',Nt
  Nt_aux=2*Nt+1
  allocate(t_grid(Nt),t_grid_aux(Nt_aux))
  !
  t_grid = linspace(tstart,tstep*real(Nt-1,8),Nt)
  t_grid_aux = linspace(tstart,0.5d0*tstep*real(Nt_aux-1,8),Nt_aux)


  
  !+- here put the option to go to the HF equilibrium   
  Ndyn=size(delta_hf)
  allocate(psit(Ndyn)); psit=0d0
  
  call init_tdHF_strides
  do ik=1,Lk
     do iso=1,Nso
        do jso=1,Nso
           isys=idelta2ivec(iso,jso,ik)
           psit(isys) = delta_hf(iso,jso,ik)
        end do
     end do
  end do
  ! eoutHF=0d0
  ! do ik=1,Lk
  !    do iso=1,Nso
  !       eoutHF = eoutHF + H_hf(iso,iso,ik)*delta_hf(iso,iso,ik)*wtk(ik) 
  !       do jso=iso+1,Nso
  !          eoutHF = eoutHF + 2d0*dreal(H_hf(iso,jso,ik)*delta_hf(iso,jso,ik))*wtk(ik) 
  !       end do
  !    end do
  ! end do
  ! write(436,*) eoutHF

  allocate(H_hf_dyn(Nso,Nso,Lk)); H_hf_dyn=0d0
  
  unit_io=free_unit()
  if(master) then
     open(unit=unit_io,file='TNS_dynamics.out')
     close(unit_io)
  end if

  do it=1,Nt-1
     !
     call set_HF_get_obs(psit,x_iter,eout)
     !
     if(master.and.(mod(it,20).eq.0)) then
        unit_io=free_unit()
        open(unit=unit_io,file='TNS_dynamics.out',status='old',position='append')
        write(unit_io,'(20F18.10)') t_grid(it),x_iter(1),x_iter(6),x_iter(7),eout
        close(unit_io)        
     end if

     !
     psit = RK_step(Ndyn,4,tstep,t_grid(it),psit,HF_eqs_of_motion)
     !     
  end do
  
  
  !
contains

  subroutine set_HF_get_obs(y,x_iter_dyn,eout)
    implicit none
    complex(8),intent(out)  :: x_iter_dyn(14)
    real(8),intent(out)     :: eout
    complex(8),dimension(:),allocatable,intent(in) :: y
    complex(8),dimension(:,:,:),allocatable :: delta_hf_dyn
    integer :: iso,jso,ik,isys,i,ir
    !
    if(.not.allocated(y)) stop "psi2delta not allocated"
    if(size(y).ne.Ndyn)  stop "psi2delta size(y).ne.Ndyn"
    !
    call psi2delta(y,delta_hf_dyn)    
    !
    do i=1,3
       ir=ixr(i)
       delta_hfr(:,:,ir)=0.d0
       do ik=1,Lk
          delta_hfr(:,:,ir)=delta_hfr(:,:,ir) + &
               delta_hf_dyn(:,:,ik)*exp(xi*dot_product(rpt_latt(ir,:),kpt_latt(ik,:)))*wtk(ik)
       end do
    end do
    !
    x_iter_dyn(1) = delta_hfr(1,1,ir0)+delta_hfr(1+Norb,1+Norb,ir0)
    x_iter_dyn(2) = delta_hfr(2,2,ir0)+delta_hfr(2+Norb,2+Norb,ir0)
    x_iter_dyn(3) = delta_hfr(3,3,ir0)+delta_hfr(3+Norb,3+Norb,ir0)
    !
    x_iter_dyn(4) = delta_hfr(1,3,ir0)
    x_iter_dyn(5) = delta_hfr(2,3,ir0)
    !
    x_iter_dyn(6) = delta_hfr(1,3,ir0)+delta_hfr(1,3,irR)
    x_iter_dyn(7) = delta_hfr(2,3,ir0)+delta_hfr(2,3,irR)
    !
    x_iter_dyn(8) = delta_hfr(4,4,ir0)+delta_hfr(4+Norb,4+Norb,ir0)
    x_iter_dyn(9) = delta_hfr(5,5,ir0)+delta_hfr(5+Norb,5+Norb,ir0)
    x_iter_dyn(10) = delta_hfr(6,6,ir0)+delta_hfr(6+Norb,6+Norb,ir0)
    !
    x_iter_dyn(11) = delta_hfr(4,6,ir0)
    x_iter_dyn(12) = delta_hfr(5,6,ir0)
    !
    x_iter_dyn(13) = delta_hfr(4,6,ir0)+delta_hfr(4,6,irL)
    x_iter_dyn(14) = delta_hfr(5,6,ir0)+delta_hfr(5,6,irL)
    !
    ! set the dynamical HF hamiltonian
    H_hf_dyn=HF_hamiltonian(x_iter_dyn)
    H_hf_dyn=H_hf_dyn+Hk_toy
    !
    eout=0d0
    do ik=1,Lk
       do iso=1,Nso
          eout = eout + H_hf_dyn(iso,iso,ik)*delta_hf_dyn(iso,iso,ik)*wtk(ik) 
          do jso=iso+1,Nso
             eout = eout + 2d0*dreal(H_hf_dyn(iso,jso,ik)*delta_hf_dyn(iso,jso,ik))*wtk(ik) 
          end do
       end do
    end do
    Eout=Eout-Ucell*0.25d0*(dreal(x_iter_dyn(1))**2.d0+dreal(x_iter_dyn(2))**2.d0+dreal(x_iter_dyn(3))**2.d0)
    Eout=Eout-2d0*Vcell*(dreal(x_iter_dyn(1))*dreal(x_iter_dyn(3))+dreal(x_iter_dyn(2))*dreal(x_iter_dyn(3)))
    Eout=Eout + 2.d0*Vcell*(abs(x_iter_dyn(4))**2.d0+abs(x_iter_dyn(5))**2.d0)
    Eout=Eout + 2.d0*Vcell*(abs(x_iter_dyn(6)-x_iter_dyn(4))**2.d0+abs(x_iter_dyn(7)-x_iter_dyn(5))**2.d0)
    !
    Eout=Eout-Ucell*0.25d0*(dreal(x_iter_dyn(8))**2.d0+dreal(x_iter_dyn(9))**2.d0+dreal(x_iter_dyn(10))**2.d0)
    Eout=Eout-2d0*Vcell*(dreal(x_iter_dyn(8))*dreal(x_iter_dyn(10))+dreal(x_iter_dyn(9))*dreal(x_iter_dyn(10)))
    !
    Eout=Eout + 2.d0*Vcell*(abs(x_iter_dyn(11))**2.d0+abs(x_iter_dyn(12))**2.d0)
    Eout=Eout + 2.d0*Vcell*(abs(x_iter_dyn(13)-x_iter_dyn(11))**2.d0+abs(x_iter_dyn(14)-x_iter_dyn(12))**2.d0)
    
  end subroutine set_HF_get_obs




  function fix_lgr_params2(x) result(out_x)
    implicit none
    real(8),dimension(:),intent(in)   :: x
    real(8),dimension(size(x)) :: out_x    
    complex(8),dimension(:),allocatable :: x_hf
    complex(8),dimension(:,:,:),allocatable :: Hhf_lgr,delta_hfr_lgr
    complex(8),dimension(:,:),allocatable :: dHK
    real(8) :: xmu_lgr
    complex(8) :: phi_tmp(4),lgrx(2),phi_out(2)
    integer :: ik,iorb,jorb
    real(8) :: eout
    
    if(size(x).ne.4) stop "size(x).ne.4  fix_lgr_params2"
    allocate(x_hf(size(x_iter))); x_hf=x_iter
    if(.not.allocated(Hk_toy)) stop "fix_lgr: Hk_toy not allocated"        
    allocate(Hhf_lgr(Nso,Nso,Lk)); Hhf_lgr=0.d0           
    !
    lgrx(1)=x(1)+xi*x(2)
    lgrx(2)=x(3)+xi*x(4)
    !
    Hhf_lgr = HF_hamiltonian(x_iter,lgrx)
    Hhf_lgr = Hhf_lgr + Hk_toy
    !
    call fix_mu(Hhf_lgr,delta_hf,mu_fix)
    !
    allocate(delta_hfr_lgr(Nso,Nso,nrpts)); delta_hfr_lgr=0d0
    do i=1,3
       ir=ixr(i)
       delta_hfr_lgr(:,:,ir)=0.d0
       do ik=1,Lk
          delta_hfr_lgr(:,:,ir)=delta_hfr_lgr(:,:,ir) + &
               delta_hf(:,:,ik)*exp(xi*dot_product(rpt_latt(ir,:),kpt_latt(ik,:)))*wtk(ik)
       end do
    end do
    !
    phi_tmp(1) = delta_hfr_lgr(1,3,ir0)+delta_hfr_lgr(1,3,irR)
    phi_tmp(2) = delta_hfr_lgr(2,3,ir0)+delta_hfr_lgr(2,3,irR)
    !
    phi_tmp(3) = delta_hfr_lgr(4,6,ir0)+delta_hfr_lgr(4,6,irL)
    phi_tmp(4) = delta_hfr_lgr(5,6,ir0)+delta_hfr_lgr(5,6,irL)
    !
    phi_out(1) = phi_tmp(1)
    phi_out(2) = phi_tmp(2)
    !
    phi_out = phi_out - xphi_vec
    if(master.and.lgr_fix_verbose>0) then
       write(*,'(10F18.10)') phi_tmp
       write(*,'(10F18.10)') phi_out,lgrx
       write(*,'(10F18.10)')
    end if
    !
    out_x(1) = dreal(phi_out(1))
    out_x(2) = dimag(phi_out(1))
    out_x(3) = dreal(phi_out(2))
    out_x(4) = dimag(phi_out(2))
    !    
  end function fix_lgr_params2
  
  !
  function fix_lgr_params(x) result(out_x)
    implicit none
    real(8),dimension(:),intent(in)   :: x
    real(8),dimension(size(x)) :: out_x    
    complex(8),dimension(:),allocatable :: x_hf
    complex(8),dimension(:,:,:),allocatable :: Hhf_lgr,delta_hfr_lgr
    complex(8),dimension(:,:),allocatable :: dHK
    real(8) :: xmu_lgr
    complex(8) :: phi_tmp(4),lgrx(2),phi_out
    integer :: ik,iorb,jorb
    real(8) :: eout
    
    if(size(x).ne.2) stop "size(x).ne.2  fix lgr params"
    allocate(x_hf(size(x_iter))); x_hf=x_iter
    if(.not.allocated(Hk_toy)) stop "fix_lgr: Hk_toy not allocated"        
    allocate(Hhf_lgr(Nso,Nso,Lk)); Hhf_lgr=0.d0           
    !
    lgrx=x(1)+xi*x(2)
    Hhf_lgr = HF_hamiltonian(x_iter,lgrx)
    Hhf_lgr = Hhf_lgr + Hk_toy
    !
    call fix_mu(Hhf_lgr,delta_hf,mu_fix)
    !
    allocate(delta_hfr_lgr(Nso,Nso,nrpts)); delta_hfr_lgr=0d0
    do i=1,3
       ir=ixr(i)
       delta_hfr_lgr(:,:,ir)=0.d0
       do ik=1,Lk
          delta_hfr_lgr(:,:,ir)=delta_hfr_lgr(:,:,ir) + &
               delta_hf(:,:,ik)*exp(xi*dot_product(rpt_latt(ir,:),kpt_latt(ik,:)))*wtk(ik)
       end do
    end do
    !
    phi_tmp(1) = delta_hfr_lgr(1,3,ir0)+delta_hfr_lgr(1,3,irR)
    phi_tmp(2) = delta_hfr_lgr(2,3,ir0)+delta_hfr_lgr(2,3,irR)
    !
    phi_tmp(3) = delta_hfr_lgr(4,6,ir0)+delta_hfr_lgr(4,6,irL)
    phi_tmp(4) = delta_hfr_lgr(5,6,ir0)+delta_hfr_lgr(5,6,irL)
    !
    phi_out = phi_tmp(1)!dot_product(phi_sgn,phi_tmp)/4d0
    phi_out = phi_out - xphi
    if(master.and.lgr_fix_verbose>0) then
       write(*,'(10F18.10)') phi_tmp
       write(*,'(10F18.10)') phi_out,lgrx
       write(*,'(10F18.10)')
    end if
    !
    out_x(1) = dreal(phi_out)
    out_x(2) = dimag(phi_out)
    !    
  end function fix_lgr_params


  function fix_lgr_params_(x) result(out_x)
    implicit none
    real(8),dimension(:),intent(in)   :: x
    real(8),dimension(size(x)) :: out_x    
    complex(8),dimension(:),allocatable :: x_hf
    complex(8),dimension(:,:,:),allocatable :: Hhf_lgr,delta_hfr_lgr
    complex(8),dimension(:,:),allocatable :: dHK
    real(8) :: xmu_lgr
    complex(8) :: phi_tmp(4),lgrx(2),phi_out
    integer :: ik,iorb,jorb
    real(8) :: eout
    
    if(size(x).ne.1) stop "size(x).ne.1  fix lgr params_"
    allocate(x_hf(size(x_iter))); x_hf=x_iter
    if(.not.allocated(Hk_toy)) stop "fix_lgr: Hk_toy not allocated"        
    allocate(Hhf_lgr(Nso,Nso,Lk)); Hhf_lgr=0.d0           
    !
    lgrx = x(1)*dreal(xphi)/abs(xphi)+xi*x(1)*dimag(xphi)/abs(xphi)
    !
    Hhf_lgr = HF_hamiltonian(x_iter,lgrx)
    Hhf_lgr = Hhf_lgr + Hk_toy
    !
    call fix_mu(Hhf_lgr,delta_hf,mu_fix)
    !
    allocate(delta_hfr_lgr(Nso,Nso,nrpts)); delta_hfr_lgr=0d0
    do i=1,3
       ir=ixr(i)
       delta_hfr_lgr(:,:,ir)=0.d0
       do ik=1,Lk
          delta_hfr_lgr(:,:,ir)=delta_hfr_lgr(:,:,ir) + &
               delta_hf(:,:,ik)*exp(xi*dot_product(rpt_latt(ir,:),kpt_latt(ik,:)))*wtk(ik)
       end do
    end do
    !
    phi_tmp(1) = delta_hfr_lgr(1,3,ir0)+delta_hfr_lgr(1,3,irR)
    phi_tmp(2) = delta_hfr_lgr(2,3,ir0)+delta_hfr_lgr(2,3,irR)
    !
    phi_tmp(3) = delta_hfr_lgr(4,6,ir0)+delta_hfr_lgr(4,6,irL)
    phi_tmp(4) = delta_hfr_lgr(5,6,ir0)+delta_hfr_lgr(5,6,irL)
    !
    phi_out = phi_tmp(1)!dot_product(phi_sgn,phi_tmp)/4d0
    phi_out = phi_out - xphi
    if(master.and.lgr_fix_verbose>0) then
       write(*,'(10F18.10)') phi_tmp
       write(*,'(10F18.10)') phi_out,lgrx
       write(*,'(10F18.10)')
    end if
    !
    out_x(1) = dreal(phi_out)
    !out_x(2) = dimag(phi_out)
    !    
  end function fix_lgr_params_

  subroutine delta2psi(delta,y)
    integer :: ik,iso,jso,isys
    complex(8),dimension(:,:,:),allocatable,intent(in) :: delta
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
  end subroutine delta2psi
  subroutine psi2delta(y,delta)
    integer :: ik,iso,jso,isys
    complex(8),dimension(:,:,:),allocatable,intent(inout) :: delta
    complex(8),dimension(:),allocatable,intent(in) :: y

    if(.not.allocated(y)) stop "psi2delta not allocated"
    if(size(y).ne.Ndyn)  stop "psi2delta size(y).ne.Ndyn"
    !
    if(allocated(delta)) deallocate(delta)
    allocate(delta(Nso,Nso,Lk));delta=0d0
    do isys=1,Ndyn
       iso=ivec2idelta(isys,1)
       jso=ivec2idelta(isys,2)
       ik=ivec2idelta(isys,3)
       delta(iso,jso,ik) = y(isys)
    end do
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
       yout(isys) = 0d0
       do jso_=1,Nso
          !
          jsys=idelta2ivec(iso,jso_,ik)
          yout(isys) = yout(isys) + H_hf_dyn(jso,jso_,ik)*yin(jsys)
          !
          jsys=idelta2ivec(jso_,jso,ik)          
          yout(isys) = yout(isys) - H_hf_dyn(jso_,iso,ik)*yin(jsys)
          !
       end do
    end do
    yout=-xi*yout/(hbar_ev_ps*1d3) !+- here the time is measured in fs
    !
  end function HF_eqs_of_motion

  
  !
  function HF_hamiltonian(x_iter,lgr_) result(Hhf)
    implicit none
    complex(8) :: x_iter(14)
    complex(8),dimension(Nso,Nso,Lk) :: Hhf
    complex(8),optional :: lgr_(2)
    complex(8)          :: lgr(2)
    integer :: iso,jso,iorb,jorb
    !
    lgr=0d0
    if(present(lgr_)) lgr=lgr_
    
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
          iorb=1; iso=(ispin-1)*Norb+iorb          
          jorb=3; jso=(ispin-1)*Norb+jorb
          Hhf(jso,iso,ik) = -Vcell*x_iter(4)*(1.d0-exp(-xi*dot_product(R1,kpt_latt(ik,:))))
          Hhf(jso,iso,ik) = Hhf(jso,iso,ik) - Vcell*x_iter(6)*exp(-xi*dot_product(R1,kpt_latt(ik,:)))
          Hhf(jso,iso,ik) = Hhf(jso,iso,ik) + conjg(lgr(1)*phi_sgn(1))*(1.d0+exp(-xi*dot_product(R1,kpt_latt(ik,:))))
          Hhf(iso,jso,ik) = conjg(Hhf(jso,iso,ik))
          !
          iorb=2; iso=(ispin-1)*Norb+iorb          
          jorb=3; jso=(ispin-1)*Norb+jorb
          Hhf(jso,iso,ik) = -Vcell*x_iter(5)*(1.d0-exp(-xi*dot_product(R1,kpt_latt(ik,:))))
          Hhf(jso,iso,ik) = Hhf(jso,iso,ik) - Vcell*x_iter(7)*exp(-xi*dot_product(R1,kpt_latt(ik,:)))
          Hhf(jso,iso,ik) = Hhf(jso,iso,ik) + conjg(lgr(2)*phi_sgn(2))*(1.d0+exp(-xi*dot_product(R1,kpt_latt(ik,:))))
          Hhf(iso,jso,ik) = conjg(Hhf(jso,iso,ik))
          !
          !
          iorb=4
          iso=(ispin-1)*Norb+iorb
          Hhf(iso,iso,ik) = 0.5d0*Ucell*x_iter(8)+2.d0*Vcell*x_iter(10)
          !
          iorb=5
          iso=(ispin-1)*Norb+iorb
          Hhf(iso,iso,ik) = 0.5d0*Ucell*x_iter(9)+2.d0*Vcell*x_iter(10)          
          !
          iorb=6
          iso=(ispin-1)*Norb+iorb
          Hhf(iso,iso,ik) = 0.5d0*Ucell*x_iter(10)+2.d0*Vcell*(x_iter(9)+x_iter(8))
          !
          iorb=4; iso=(ispin-1)*Norb+iorb          
          jorb=6; jso=(ispin-1)*Norb+jorb
          Hhf(jso,iso,ik) = -Vcell*x_iter(11)*(1.d0-exp(xi*dot_product(R1,kpt_latt(ik,:))))
          Hhf(jso,iso,ik) = Hhf(jso,iso,ik) - Vcell*x_iter(13)*exp(xi*dot_product(R1,kpt_latt(ik,:)))
          Hhf(jso,iso,ik) = Hhf(jso,iso,ik) + conjg(lgr(2)*phi_sgn(3))*(1.d0+exp(xi*dot_product(R1,kpt_latt(ik,:))))
          Hhf(iso,jso,ik) = conjg(Hhf(jso,iso,ik))
          !
          iorb=5; iso=(ispin-1)*Norb+iorb          
          jorb=6; jso=(ispin-1)*Norb+jorb
          Hhf(jso,iso,ik) = -Vcell*x_iter(12)*(1.d0-exp(xi*dot_product(R1,kpt_latt(ik,:))))
          Hhf(jso,iso,ik) = Hhf(jso,iso,ik) - Vcell*x_iter(14)*exp(xi*dot_product(R1,kpt_latt(ik,:)))
          Hhf(jso,iso,ik) = Hhf(jso,iso,ik) + conjg(lgr(1)*phi_sgn(4))*(1.d0+exp(xi*dot_product(R1,kpt_latt(ik,:))))
          !
          Hhf(iso,jso,ik) = conjg(Hhf(jso,iso,ik))
          !
       end do
    end do
    !
  end function HF_hamiltonian
  !
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

