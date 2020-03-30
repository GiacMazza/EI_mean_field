program officina
  !
  USE SCIFOR
  USE DMFT_TOOLS
  !
  USE VARS_GLOBAL
  USE DMFT_VECTORS
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
  !real(8),dimension(3,3) :: Bkinv

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

  real(8) :: Eout




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
  complex(8),dimension(14) :: xhf,xhf_
  complex(8),dimension(9) :: xtmp,xtmp_
  complex(8) :: xphi,xphi_
  complex(8),dimension(2) :: xpi,xpi_
  real(8),dimension(14) :: xr_iter
  real(8),dimension(18) :: xr_tmp
  complex(8),dimension(4) :: phi_target


  real(8) :: Ucell,Vcell,Wcell
  real(8) :: Evalence,Econduction,tconduction,tvalence,tt_hyb,nn_hyb,tn_hyb
  real(8) :: w0gap
  logical :: fix_phi  !+-> HF calculation with fixed order parameter <-+!
  logical :: hf_symm  !+-> HF calculation w/ symmetric order parameter <-+!
  logical :: HF_solve !+-> HF calculation with a root-finder routine  <-+
  logical :: hf_in

  real(8) :: phi_start,phi_end,dphi,ntot,test,err
  integer :: Nphi
  integer :: ixr(3)
  !
  real(8) :: op_phase
  real(8),dimension(2) :: lgr_phase_min
  real(8),dimension(8) :: lgr_phase,lgr_phase_
  complex(8),dimension(4) :: lgr_op
  real(8),dimension(8) :: lgr_tmp
  real(8) :: check_lgr
  integer :: info_lgr
  
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
  call parse_input_variable(Nhf_,"Nhf_","input.conf",default=10)
  call parse_input_variable(Lreal,"Lreal","input.conf",default=2000)
  call parse_input_variable(wmix,"wmix","input.conf",default=0.5d0)
  call parse_input_variable(fileRlattice,"R_unit_cell","input.conf",default='R_unit_cell.conf')
  call parse_input_variable(file_w90_hr,"w90_hr_input","input.conf",default='wannier90_hr.dat')
  call parse_input_variable(file_UV,"UVq_input","input.conf",default='VU.dat')  
  call parse_input_variable(read_tns,"READ_TNS","input.conf", &
       default='/home/mazza/project_data/TNS_EI/MALTE_DATA/K_12x12x3_B_120/')
  call parse_input_variable(Nint,"Nint","input.conf",default=1)
  call parse_input_variable(hartree,"Hartree","input.conf",default=.true.)  
  call parse_input_variable(cf_ext,"cf_ext","input.conf",default=0.d0)
  call parse_input_variable(Ucut_off,"U_CUT","input.conf",default=20.d0)
  call parse_input_variable(H1d,"H1D","input.conf",default=.false.)
  call parse_input_variable(alphaU,"alphaU","input.conf",default=1.d0)

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
  call parse_input_variable(phi_start,"phi_start","input.conf",default=0.5d0)
  call parse_input_variable(phi_end,"phi_end","input.conf",default=0.d0)
  call parse_input_variable(Nphi,"Nphi","input.conf",default=50)


  call parse_input_variable(hf_symm,"hf_symm","input.conf",default=.false.)
  call parse_input_variable(hf_solve,"hf_solve","input.conf",default=.false.)

  call parse_input_variable(hf_conv,"hf_conv","input.conf",default=1.d-06)
  call parse_input_variable(op_phase,"op_phase","input.conf",default=0.d0)

  
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
  call read_w90_hr(R1,R2,R3,Hr_w90,Hloc,irvec,ndegen,trim(file_name),1,6,1)
  !
  nrpts=Nk_x*Nk_y
  write(*,*) 'nrpts',nrpts
  deallocate(irvec); allocate(irvec(nrpts,3))
  deallocate(ndegen); allocate(ndegen(nrpts));ndegen=1
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
  allocate(Hr_toy(Nso,Nso,nrpts))
  do ir=1,nrpts
     call FT_q2r(rpt_latt(ir,:),Hr_toy(:,:,ir),Hk_toy)
  end do
  !
  Hk_toy=0.d0
  do ik=1,Lk
     call FT_r2q(kpt_latt(ik,:),Hk_toy(:,:,ik),Hr_toy)
  end do  
  mu_fix=0.d0
  call fix_mu(Hk_toy,delta_hf,mu_fix)  
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
     write(unit_in,'(100F10.5)') Rlat(1),dreal(obs(:)),dreal(obs(:))

     !
  end do
  ! 
  !
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

  !stop  
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
  !stop
  !
  allocate(H_Hf(Nso,Nso,Lk))
  !

  !+- first optimize order-parameter in the ground state -+!
  write(*,*) "starting the optimization of the ground state order parameter"
  !
  inquire(file='hf_BLS_free_final.out',exist=hf_in)  
  if(hf_in) then
     !
     flen=file_length('hf_BLS_free_final.out')
     if(flen.eq.14) then
        call read_array('hf_BLS_free_final.out',x_iter)
     end if
     !
  else

     x_iter(1) = 0.1d0
     x_iter(2) = 0.1d0
     x_iter(3) = 1.8d0
     !
     x_iter(4) =  0.1d0*exp(0.5d0*pi*xi)
     x_iter(5) =  0.1d0*exp(0.5d0*pi*xi)
     !  
     x_iter(6) =  0.4d0*exp(0.5d0*pi*xi)
     x_iter(7) = -0.4d0*exp(0.5d0*pi*xi)
     !
     x_iter(8) =  0.1d0
     x_iter(9) =  0.1d0
     x_iter(10) = 1.8d0
     !
     x_iter(11) = 0.1d0*exp(0.5d0*pi*xi)
     x_iter(12) = 0.1d0*exp(0.5d0*pi*xi)
     !  
     x_iter(13) = -0.4d0*exp(0.5d0*pi*xi)
     x_iter(14) =  0.4d0*exp(0.5d0*pi*xi)     
  end if
  !
  !
  call save_array('hf_BLS_free.init',x_iter)
  !
  !
  uio=free_unit()
  open(unit=uio,file='loop_phi_BLS_free.out')
  unit_err=free_unit()
  open(unit=unit_err,file='err_BLS.out')
  !
  err_hf=1.0
  do ihf=1,Nhf
     !
     write(*,*) "BLS HF: loop nr",ihf,'/',Nhf
     write(unit_err,*) err_hf
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
     x_iter=x_iter*wmix+(1.d0-wmix)*x_iter_             
     !     
     H_Hf=HF_hamiltonian(x_iter)
     H_Hf=H_Hf+Hk_toy
     call fix_mu(H_Hf,delta_hf,mu_fix,eout)
     
     !+- double counting term -+!
     Eout=Eout-Ucell*0.25d0*(dreal(x_iter(1))**2.d0+dreal(x_iter(2))**2.d0+dreal(x_iter(3))**2.d0)
     Eout=Eout-2*Vcell*(dreal(x_iter(1))*dreal(x_iter(3))+dreal(x_iter(2))*dreal(x_iter(3)))
     !
     Eout=Eout + 2.d0*Vcell*(abs(x_iter(4))**2.d0+abs(x_iter(5))**2.d0)
     Eout=Eout + 2.d0*Vcell*(abs(x_iter(6)-x_iter(4))**2.d0+abs(x_iter(7)-x_iter(5))**2.d0)
     !
     Eout=Eout-Ucell*0.25d0*(dreal(x_iter(8))**2.d0+dreal(x_iter(9))**2.d0+dreal(x_iter(10))**2.d0)
     Eout=Eout-2*Vcell*(dreal(x_iter(8))*dreal(x_iter(10))+dreal(x_iter(9))*dreal(x_iter(10)))
     !
     Eout=Eout + 2.d0*Vcell*(abs(x_iter(11))**2.d0+abs(x_iter(12))**2.d0)
     Eout=Eout + 2.d0*Vcell*(abs(x_iter(13)-x_iter(11))**2.d0+abs(x_iter(14)-x_iter(12))**2.d0)
     !
     ntot=dreal(x_iter(1))+dreal(x_iter(2))+dreal(x_iter(3))
     ntot=ntot+dreal(x_iter(8))+dreal(x_iter(9))+dreal(x_iter(10))
     write(uio, '(40F18.10)') dreal(x_iter(1:14)),dimag(x_iter(1:14)),Eout+mu_fix*ntot,Eout
     err_hf=0.d0
     do i=1,14
        err_hf=err_hf+abs(x_iter(i)-x_iter_(i))**2.d0
     end do
     if(err_hf.lt.hf_conv) exit
     !
     !
  end do
  !
  close(unit_err)
  close(uio)
  !stop
  !
  call save_array('hf_BLS_free_final.out',x_iter)
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
     write(uio,'(100F10.5)') Rlat(1),dreal(obs(:)),dreal(obs(:))
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
        write(unit_in,'(30F18.10)') modk,ek_out-mu_fix
        write(uio,'(30F18.10)') modk,Ta_fat,Ni_fat        
        !
     end do
     !
  end do
  close(unit_in)
  close(uio)

  !+- then do the fixed phi optimzation by adding a lagrang multiplier to keep a phase on the order parameter -+!
  !
  phi_target(1) = x_iter(6)*exp(2.d0*pi*xi*op_phase)
  phi_target(2) = x_iter(7)*exp(2.d0*pi*xi*op_phase)
  phi_target(3) = x_iter(13)*exp(2.d0*pi*xi*op_phase)
  phi_target(4) = x_iter(14)*exp(2.d0*pi*xi*op_phase)


  ! x_iter(6:7) = x_iter(6:7)*exp(2.d0*pi*xi*op_phase)
  ! x_iter(13:14) = x_iter(13:14)*exp(2.d0*pi*xi*op_phase)
  
  write(*,*) 'phi_target',phi_target

  
  
  lgr_phase=0.d0
  lgr_phase_min=0.d0
  call fsolve(hf_fix_lgr_params_,lgr_phase_min,tol=1.d-10,info=info_lgr)
  !
  lgr_phase(1)= lgr_phase_min(1)
  lgr_phase(2)=-lgr_phase_min(1)
  lgr_phase(3)=-lgr_phase_min(1)
  lgr_phase(4)= lgr_phase_min(1)
  !
  lgr_phase(5)= lgr_phase_min(2)
  lgr_phase(6)=-lgr_phase_min(2)
  lgr_phase(7)=-lgr_phase_min(2)
  lgr_phase(8)= lgr_phase_min(2)
  !
  do i=1,4
     lgr_op(i) = lgr_phase(i)+xi*lgr_phase(i)
  end do     
  !     
  H_Hf=HF_hamiltonian(x_iter)
  H_Hf=H_Hf+Hk_toy
  !+- add the lgr_part of the hamiltonian -+!
  do ik=1,Lk
     do ispin=1,Nspin           
        !+- x_iter(6)
        iorb=1; iso=(ispin-1)*Norb+iorb          
        jorb=3; jso=(ispin-1)*Norb+jorb          
        H_hf(iso,jso,ik) = H_hf(iso,jso,ik) + lgr_op(1)*(1.d0+exp(-xi*dot_product(R1,kpt_latt(ik,:))))
        H_hf(jso,iso,ik) = H_hf(jso,iso,ik) + conjg(lgr_op(1)*(1.d0+exp(-xi*dot_product(R1,kpt_latt(ik,:)))))           
        ! !+- x_iter(7)
        iorb=2; iso=(ispin-1)*Norb+iorb          
        jorb=3; jso=(ispin-1)*Norb+jorb          
        H_hf(iso,jso,ik) = H_hf(iso,jso,ik) + lgr_op(2)*(1.d0+exp(-xi*dot_product(R1,kpt_latt(ik,:))))
        H_hf(jso,iso,ik) = H_hf(jso,iso,ik) + conjg(lgr_op(2)*(1.d0+exp(-xi*dot_product(R1,kpt_latt(ik,:)))))          
        !+- x_iter(13)
        iorb=4; iso=(ispin-1)*Norb+iorb          
        jorb=6; jso=(ispin-1)*Norb+jorb          
        H_hf(iso,jso,ik) = H_hf(iso,jso,ik) + lgr_op(3)*(1.d0+exp(+xi*dot_product(R1,kpt_latt(ik,:))))
        H_hf(jso,iso,ik) = H_hf(jso,iso,ik) + conjg(lgr_op(3)*(1.d0+exp(+xi*dot_product(R1,kpt_latt(ik,:)))))           
        !+- x_iter(14)
        iorb=5; iso=(ispin-1)*Norb+iorb          
        jorb=6; jso=(ispin-1)*Norb+jorb          
        H_hf(iso,jso,ik) = H_hf(iso,jso,ik) + lgr_op(4)*(1.d0+exp(+xi*dot_product(R1,kpt_latt(ik,:))))
        H_hf(jso,iso,ik) = H_hf(jso,iso,ik) + conjg(lgr_op(4)*(1.d0+exp(+xi*dot_product(R1,kpt_latt(ik,:)))))          
     end do
  end do
  call fix_mu(H_Hf,delta_hf,mu_fix,eout)
  
  !+- double counting term -+!
  Eout=Eout-Ucell*0.25d0*(dreal(x_iter(1))**2.d0+dreal(x_iter(2))**2.d0+dreal(x_iter(3))**2.d0)
  Eout=Eout-2*Vcell*(dreal(x_iter(1))*dreal(x_iter(3))+dreal(x_iter(2))*dreal(x_iter(3)))
  !
  Eout=Eout + 2.d0*Vcell*(abs(x_iter(4))**2.d0+abs(x_iter(5))**2.d0)
  Eout=Eout + 2.d0*Vcell*(abs(x_iter(6)-x_iter(4))**2.d0+abs(x_iter(7)-x_iter(5))**2.d0)
  !
  Eout=Eout-Ucell*0.25d0*(dreal(x_iter(8))**2.d0+dreal(x_iter(9))**2.d0+dreal(x_iter(10))**2.d0)
  Eout=Eout-2*Vcell*(dreal(x_iter(8))*dreal(x_iter(10))+dreal(x_iter(9))*dreal(x_iter(10)))
  !
  Eout=Eout + 2.d0*Vcell*(abs(x_iter(11))**2.d0+abs(x_iter(12))**2.d0)
  Eout=Eout + 2.d0*Vcell*(abs(x_iter(13)-x_iter(11))**2.d0+abs(x_iter(14)-x_iter(12))**2.d0)
  !
  !+- here I should remove the energy of the lagrange multiplier
  Eout=Eout-2.d0*dreal(lgr_op(1)*x_iter(6))
  Eout=Eout-2.d0*dreal(lgr_op(2)*x_iter(7))
  Eout=Eout-2.d0*dreal(lgr_op(3)*x_iter(13))
  Eout=Eout-2.d0*dreal(lgr_op(4)*x_iter(14))

  uio=free_unit()
  open(unit=uio,file='energy_vs_op_phase.out')
  write(uio,'(10F18.10)') x_iter(6),lgr_op(1),Eout
  close(uio)

  stop
  
  !
  uio=free_unit()
  open(unit=uio,file='loop_fixed_phi_phase.out')
  unit_err=free_unit()
  open(unit=unit_err,file='err_fixed_phi_phase.out')
  !
  err_hf=1.0
  do ihf=1,Nhf
     !
     write(*,*) "fixed_phi_phase HF: loop nr",ihf,'/',Nhf
     write(unit_err,*) err_hf
     x_iter_=x_iter
     lgr_phase_=lgr_phase
     !

     !+- here call fsolve
     call fsolve(fix_lgr_params,lgr_phase,tol=1.d-10,info=info_lgr)
     lgr_tmp=fix_lgr_params(lgr_phase)
     check_lgr=0.d0
     do i=1,8
        check_lgr=check_lgr+abs(lgr_tmp(i))**2.d0
     end do
     check_lgr=sqrt(check_lgr)
     write(222,*) check_lgr
     x_iter=xhf
     !
     lgr_phase=lgr_phase*wmix+(1.d0-wmix)*lgr_phase_     
     x_iter=x_iter*wmix+(1.d0-wmix)*x_iter_
     do i=1,4
        lgr_op(i) = lgr_phase(i)+xi*lgr_phase(i)
     end do
     
     !     
     H_Hf=HF_hamiltonian(x_iter)
     H_Hf=H_Hf+Hk_toy
     !+- add the lgr_part of the hamiltonian -+!
     do ik=1,Lk
        do ispin=1,Nspin           
           !+- x_iter(6)
           iorb=1; iso=(ispin-1)*Norb+iorb          
           jorb=3; jso=(ispin-1)*Norb+jorb          
           H_hf(iso,jso,ik) = H_hf(iso,jso,ik) + lgr_op(1)*(1.d0+exp(-xi*dot_product(R1,kpt_latt(ik,:))))
           H_hf(jso,iso,ik) = H_hf(jso,iso,ik) + conjg(lgr_op(1)*(1.d0+exp(-xi*dot_product(R1,kpt_latt(ik,:)))))           
           ! !+- x_iter(7)
           iorb=2; iso=(ispin-1)*Norb+iorb          
           jorb=3; jso=(ispin-1)*Norb+jorb          
           H_hf(iso,jso,ik) = H_hf(iso,jso,ik) + lgr_op(2)*(1.d0+exp(-xi*dot_product(R1,kpt_latt(ik,:))))
           H_hf(jso,iso,ik) = H_hf(jso,iso,ik) + conjg(lgr_op(2)*(1.d0+exp(-xi*dot_product(R1,kpt_latt(ik,:)))))          
           !+- x_iter(13)
           iorb=4; iso=(ispin-1)*Norb+iorb          
           jorb=6; jso=(ispin-1)*Norb+jorb          
           H_hf(iso,jso,ik) = H_hf(iso,jso,ik) + lgr_op(3)*(1.d0+exp(+xi*dot_product(R1,kpt_latt(ik,:))))
           H_hf(jso,iso,ik) = H_hf(jso,iso,ik) + conjg(lgr_op(3)*(1.d0+exp(+xi*dot_product(R1,kpt_latt(ik,:)))))           
           !+- x_iter(14)
           iorb=5; iso=(ispin-1)*Norb+iorb          
           jorb=6; jso=(ispin-1)*Norb+jorb          
           H_hf(iso,jso,ik) = H_hf(iso,jso,ik) + lgr_op(4)*(1.d0+exp(+xi*dot_product(R1,kpt_latt(ik,:))))
           H_hf(jso,iso,ik) = H_hf(jso,iso,ik) + conjg(lgr_op(4)*(1.d0+exp(+xi*dot_product(R1,kpt_latt(ik,:)))))          
        end do
     end do
     call fix_mu(H_Hf,delta_hf,mu_fix,eout)
     
     !+- double counting term -+!
     Eout=Eout-Ucell*0.25d0*(dreal(x_iter(1))**2.d0+dreal(x_iter(2))**2.d0+dreal(x_iter(3))**2.d0)
     Eout=Eout-2*Vcell*(dreal(x_iter(1))*dreal(x_iter(3))+dreal(x_iter(2))*dreal(x_iter(3)))
     !
     Eout=Eout + 2.d0*Vcell*(abs(x_iter(4))**2.d0+abs(x_iter(5))**2.d0)
     Eout=Eout + 2.d0*Vcell*(abs(x_iter(6)-x_iter(4))**2.d0+abs(x_iter(7)-x_iter(5))**2.d0)
     !
     Eout=Eout-Ucell*0.25d0*(dreal(x_iter(8))**2.d0+dreal(x_iter(9))**2.d0+dreal(x_iter(10))**2.d0)
     Eout=Eout-2*Vcell*(dreal(x_iter(8))*dreal(x_iter(10))+dreal(x_iter(9))*dreal(x_iter(10)))
     !
     Eout=Eout + 2.d0*Vcell*(abs(x_iter(11))**2.d0+abs(x_iter(12))**2.d0)
     Eout=Eout + 2.d0*Vcell*(abs(x_iter(13)-x_iter(11))**2.d0+abs(x_iter(14)-x_iter(12))**2.d0)
     !
     !+- here I should remove the energy of the lagrange multiplier
     Eout=Eout-2.d0*dreal(lgr_op(1)*x_iter(6))
     Eout=Eout-2.d0*dreal(lgr_op(2)*x_iter(7))
     Eout=Eout-2.d0*dreal(lgr_op(3)*x_iter(13))
     Eout=Eout-2.d0*dreal(lgr_op(4)*x_iter(14))
     !
     
     ntot=dreal(x_iter(1))+dreal(x_iter(2))+dreal(x_iter(3))
     ntot=ntot+dreal(x_iter(8))+dreal(x_iter(9))+dreal(x_iter(10))
     write(uio, '(40F18.10)') dreal(x_iter(1:14)),dimag(x_iter(1:14)),Eout+mu_fix*ntot,Eout
     err_hf=0.d0
     do i=1,14
        err_hf=err_hf+abs(x_iter(i)-x_iter_(i))**2.d0
     end do
     do i=1,8
        err_hf=err_hf+abs(lgr_phase(i)-lgr_phase_(i))**2.d0
     end do

     !write(*,*) 'err',err_hf,hf_conv 
     if(err_hf.lt.hf_conv) exit
     !
     !
  end do
  !
  close(unit_err)
  close(uio)
  
  


  stop
  
  !
contains
  !
  function hf_fix_lgr_params_(xlgr) result(xdelta_op)
    real(8),dimension(:),intent(in) :: xlgr
    !real(8),dimension(:) :: xlgr
    real(8),dimension(size(xlgr)) :: xdelta_op
    complex(8),dimension(:,:,:),allocatable :: Hhf
    complex(8),dimension(:,:,:),allocatable :: deltahfr,deltahf
    real(8)   :: mu
    complex(8),dimension(4) :: lgr,delta_op
    integer :: i,ir,ik,iorb,jorb,ispin,iso,jso
    !

    !+-> for a given value of the lagrange parameters solve HF @ self-consistency -+!    
    if(size(xlgr).ne.2) stop "wrong sized in xlgr"
    ! do i=1,4
    !    lgr(i) = xlgr(i)+xi*xlgr(i+4)
    ! end do

    lgr(1) = xlgr(1)+xi*xlgr(2)
    lgr(2) = -lgr(1)
    lgr(3) = lgr(2)
    lgr(4) = lgr(1)
    !
    allocate(Hhf(Nso,Nso,Lk),deltahf(Nso,Nso,Lk),deltahfr(Nso,Nso,nrpts))



    !+- initialize x_iter -+!
    !+- initialize to the HF solution at zero lgr;
    xhf=x_iter
    err_hf=1.0
    do ihf=1,Nhf
       !
       write(444,'(20F18.10)') xhf(6),xhf(7),xhf(13),xhf(14),lgr
       write(555,'(20F18.10)') err_hf,xhf(6),xhf(7),xhf(13),xhf(14)
       xhf_=xhf
       !    !
       HHf = HF_hamiltonian(xhf)
       HHf = HHf+Hk_toy
       !
       !+-  here add the part of the lgr multipliers to the operator of the order parameter -+!
       do ik=1,Lk
          do ispin=1,Nspin          
             !+- x_iter(6)
             iorb=1; iso=(ispin-1)*Norb+iorb          
             jorb=3; jso=(ispin-1)*Norb+jorb          
             Hhf(iso,jso,ik) = Hhf(iso,jso,ik) + lgr(1)*(1.d0+exp(-xi*dot_product(R1,kpt_latt(ik,:))))
             Hhf(jso,iso,ik) = Hhf(jso,iso,ik) + conjg(lgr(1)*(1.d0+exp(-xi*dot_product(R1,kpt_latt(ik,:)))))
             
             ! !+- x_iter(7)
             iorb=2; iso=(ispin-1)*Norb+iorb          
             jorb=3; jso=(ispin-1)*Norb+jorb          
             Hhf(iso,jso,ik) = Hhf(iso,jso,ik) + lgr(2)*(1.d0+exp(-xi*dot_product(R1,kpt_latt(ik,:))))
             Hhf(jso,iso,ik) = Hhf(jso,iso,ik) + conjg(lgr(2)*(1.d0+exp(-xi*dot_product(R1,kpt_latt(ik,:)))))
             
             !+- x_iter(13)
             iorb=4; iso=(ispin-1)*Norb+iorb          
             jorb=6; jso=(ispin-1)*Norb+jorb          
             Hhf(iso,jso,ik) = Hhf(iso,jso,ik) + lgr(3)*(1.d0+exp(+xi*dot_product(R1,kpt_latt(ik,:))))
             Hhf(jso,iso,ik) = Hhf(jso,iso,ik) + conjg(lgr(3)*(1.d0+exp(+xi*dot_product(R1,kpt_latt(ik,:)))))
             
             !+- x_iter(14)
             iorb=5; iso=(ispin-1)*Norb+iorb          
             jorb=6; jso=(ispin-1)*Norb+jorb          
             Hhf(iso,jso,ik) = Hhf(iso,jso,ik) + lgr(4)*(1.d0+exp(+xi*dot_product(R1,kpt_latt(ik,:))))
             Hhf(jso,iso,ik) = Hhf(jso,iso,ik) + conjg(lgr(4)*(1.d0+exp(+xi*dot_product(R1,kpt_latt(ik,:)))))
             
          end do
       end do
       !
       mu=0.d0
       call fix_mu(HHf,deltahf,mu) ; mu_fix=mu
       !
       deltahfr=0.d0
       do i=1,3!
          ir=ixr(i)
          do ik=1,Lk
             deltahfr(:,:,ir)=deltahfr(:,:,ir) + &
                  deltahf(:,:,ik)*exp(xi*dot_product(rpt_latt(ir,:),kpt_latt(ik,:)))*wtk(ik)
          end do
       end do
       !
       xhf(1) = deltahfr(1,1,ir0)+deltahfr(1+Norb,1+Norb,ir0)
       xhf(2) = deltahfr(2,2,ir0)+deltahfr(2+Norb,2+Norb,ir0)
       xhf(3) = deltahfr(3,3,ir0)+deltahfr(3+Norb,3+Norb,ir0)
       !
       xhf(4) = deltahfr(1,3,ir0)
       xhf(5) = deltahfr(2,3,ir0)
       !
       xhf(6) = deltahfr(1,3,ir0)+deltahfr(1,3,irR)
       xhf(7) = deltahfr(2,3,ir0)+deltahfr(2,3,irR)
       !
       xhf(8) = deltahfr(4,4,ir0)+deltahfr(4+Norb,4+Norb,ir0)
       xhf(9) = deltahfr(5,5,ir0)+deltahfr(5+Norb,5+Norb,ir0)
       xhf(10) = deltahfr(6,6,ir0)+deltahfr(6+Norb,6+Norb,ir0)
       !
       xhf(11) = deltahfr(4,6,ir0)
       xhf(12) = deltahfr(5,6,ir0)
       !
       xhf(13) = deltahfr(4,6,ir0)+deltahfr(4,6,irL)
       xhf(14) = deltahfr(5,6,ir0)+deltahfr(5,6,irL)
       !
       xhf=xhf*wmix+(1.d0-wmix)*xhf_
       !
       !
       err_hf=0.d0
       do i=1,14
          err_hf=err_hf+abs(xhf(i)-xhf_(i))**2.d0
       end do

       if(err_hf.lt.hf_conv) exit
       !
       !
    end do
    write(555,*)
    write(555,*)
    !
    delta_op(1) = xhf(6)  - phi_target(1)
    delta_op(2) = xhf(7)  - phi_target(2)
    delta_op(3) = xhf(13) - phi_target(3)
    delta_op(4) = xhf(14) - phi_target(4)    
    !
    do i=1,1
       xdelta_op(i) = dreal(delta_op(i))
       xdelta_op(i+1) = dimag(delta_op(i))
    end do
  end function hf_fix_lgr_params_




    function hf_fix_lgr_params(xlgr) result(xdelta_op)
    real(8),dimension(:),intent(in) :: xlgr
    !real(8),dimension(:) :: xlgr
    real(8),dimension(size(xlgr)) :: xdelta_op
    complex(8),dimension(:,:,:),allocatable :: Hhf
    complex(8),dimension(:,:,:),allocatable :: deltahfr,deltahf
    real(8)   :: mu
    complex(8),dimension(4) :: lgr,delta_op
    integer :: i,ir,ik,iorb,jorb,ispin,iso,jso
    !

    !+-> for a given value of the lagrange parameters solve HF @ self-consistency -+!    
    if(size(xlgr).ne.8) stop "wrong sized in xlgr"
    do i=1,4
       lgr(i) = xlgr(i)+xi*xlgr(i+4)
    end do
    !
    allocate(Hhf(Nso,Nso,Lk),deltahf(Nso,Nso,Lk),deltahfr(Nso,Nso,nrpts))



    !+- initialize x_iter -+!
    !+- initialize to the HF solution at zero lgr;
    xhf=x_iter
    err_hf=1.0
    do ihf=1,Nhf
       !
       write(444,'(20F18.10)') xhf(6),xhf(7),xhf(13),xhf(14),lgr
       write(555,'(20F18.10)') err_hf,xhf(6),xhf(7),xhf(13),xhf(14)
       xhf_=xhf
       !    !
       HHf = HF_hamiltonian(xhf)
       HHf = HHf+Hk_toy
       !
       !+-  here add the part of the lgr multipliers to the operator of the order parameter -+!
       do ik=1,Lk
          do ispin=1,Nspin          
             !+- x_iter(6)
             iorb=1; iso=(ispin-1)*Norb+iorb          
             jorb=3; jso=(ispin-1)*Norb+jorb          
             Hhf(iso,jso,ik) = Hhf(iso,jso,ik) + lgr(1)*(1.d0+exp(-xi*dot_product(R1,kpt_latt(ik,:))))
             Hhf(jso,iso,ik) = Hhf(jso,iso,ik) + conjg(lgr(1)*(1.d0+exp(-xi*dot_product(R1,kpt_latt(ik,:)))))
             
             ! !+- x_iter(7)
             iorb=2; iso=(ispin-1)*Norb+iorb          
             jorb=3; jso=(ispin-1)*Norb+jorb          
             Hhf(iso,jso,ik) = Hhf(iso,jso,ik) + lgr(2)*(1.d0+exp(-xi*dot_product(R1,kpt_latt(ik,:))))
             Hhf(jso,iso,ik) = Hhf(jso,iso,ik) + conjg(lgr(2)*(1.d0+exp(-xi*dot_product(R1,kpt_latt(ik,:)))))
             
             !+- x_iter(13)
             iorb=4; iso=(ispin-1)*Norb+iorb          
             jorb=6; jso=(ispin-1)*Norb+jorb          
             Hhf(iso,jso,ik) = Hhf(iso,jso,ik) + lgr(3)*(1.d0+exp(+xi*dot_product(R1,kpt_latt(ik,:))))
             Hhf(jso,iso,ik) = Hhf(jso,iso,ik) + conjg(lgr(3)*(1.d0+exp(+xi*dot_product(R1,kpt_latt(ik,:)))))
             
             !+- x_iter(14)
             iorb=5; iso=(ispin-1)*Norb+iorb          
             jorb=6; jso=(ispin-1)*Norb+jorb          
             Hhf(iso,jso,ik) = Hhf(iso,jso,ik) + lgr(4)*(1.d0+exp(+xi*dot_product(R1,kpt_latt(ik,:))))
             Hhf(jso,iso,ik) = Hhf(jso,iso,ik) + conjg(lgr(4)*(1.d0+exp(+xi*dot_product(R1,kpt_latt(ik,:)))))
             
          end do
       end do
       !
       mu=0.d0
       call fix_mu(HHf,deltahf,mu) ; mu_fix=mu
       !
       deltahfr=0.d0
       do i=1,3!
          ir=ixr(i)
          do ik=1,Lk
             deltahfr(:,:,ir)=deltahfr(:,:,ir) + &
                  deltahf(:,:,ik)*exp(xi*dot_product(rpt_latt(ir,:),kpt_latt(ik,:)))*wtk(ik)
          end do
       end do
       !
       xhf(1) = deltahfr(1,1,ir0)+deltahfr(1+Norb,1+Norb,ir0)
       xhf(2) = deltahfr(2,2,ir0)+deltahfr(2+Norb,2+Norb,ir0)
       xhf(3) = deltahfr(3,3,ir0)+deltahfr(3+Norb,3+Norb,ir0)
       !
       xhf(4) = deltahfr(1,3,ir0)
       xhf(5) = deltahfr(2,3,ir0)
       !
       xhf(6) = deltahfr(1,3,ir0)+deltahfr(1,3,irR)
       xhf(7) = deltahfr(2,3,ir0)+deltahfr(2,3,irR)
       !
       xhf(8) = deltahfr(4,4,ir0)+deltahfr(4+Norb,4+Norb,ir0)
       xhf(9) = deltahfr(5,5,ir0)+deltahfr(5+Norb,5+Norb,ir0)
       xhf(10) = deltahfr(6,6,ir0)+deltahfr(6+Norb,6+Norb,ir0)
       !
       xhf(11) = deltahfr(4,6,ir0)
       xhf(12) = deltahfr(5,6,ir0)
       !
       xhf(13) = deltahfr(4,6,ir0)+deltahfr(4,6,irL)
       xhf(14) = deltahfr(5,6,ir0)+deltahfr(5,6,irL)
       !
       xhf=xhf*wmix+(1.d0-wmix)*xhf_
       !
       !
       err_hf=0.d0
       do i=1,14
          err_hf=err_hf+abs(xhf(i)-xhf_(i))**2.d0
       end do

       if(err_hf.lt.hf_conv) exit
       !
       !
    end do
    write(555,*)
    write(555,*)
    !
    delta_op(1) = xhf(6)  - phi_target(1)
    delta_op(2) = xhf(7)  - phi_target(2)
    delta_op(3) = xhf(13) - phi_target(3)
    delta_op(4) = xhf(14) - phi_target(4)    
    !
    do i=1,4
       xdelta_op(i) = dreal(delta_op(i))
       xdelta_op(i+4) = dimag(delta_op(i))
    end do
  end function hf_fix_lgr_params





  
  function fix_lgr_params(xlgr) result(xdelta_op)
    !real(8),dimension(:),intent(in) :: xlgr
    real(8),dimension(:) :: xlgr
    real(8),dimension(size(xlgr)) :: xdelta_op
    complex(8),dimension(:,:,:),allocatable :: Hhf
    complex(8),dimension(:,:,:),allocatable :: deltahfr,deltahf
    real(8)   :: mu
    complex(8),dimension(4) :: lgr,delta_op
    integer :: i,ir,ik,iorb,jorb,ispin,iso,jso
    !
    if(size(xlgr).ne.8) stop "wrong sized in xlgr"
    do i=1,4
       lgr(i) = xlgr(i)+xi*xlgr(i+4)
    end do
    !
    allocate(Hhf(Nso,Nso,Lk),deltahf(Nso,Nso,Lk),deltahfr(Nso,Nso,nrpts))
    xhf=x_iter
    !
    HHf = HF_hamiltonian(xhf)
    HHf = HHf+Hk_toy
    !
    !+-  here add the part of the lgr multipliers to the operator of the order parameter -+!
    do ik=1,Lk
       do ispin=1,Nspin
          
          !+- x_iter(6)
          iorb=1; iso=(ispin-1)*Norb+iorb          
          jorb=3; jso=(ispin-1)*Norb+jorb          
          Hhf(iso,jso,ik) = Hhf(iso,jso,ik) + lgr(1)*(1.d0+exp(-xi*dot_product(R1,kpt_latt(ik,:))))
          Hhf(jso,iso,ik) = Hhf(jso,iso,ik) + conjg(lgr(1)*(1.d0+exp(-xi*dot_product(R1,kpt_latt(ik,:)))))
          
          ! !+- x_iter(7)
          iorb=2; iso=(ispin-1)*Norb+iorb          
          jorb=3; jso=(ispin-1)*Norb+jorb          
          Hhf(iso,jso,ik) = Hhf(iso,jso,ik) + lgr(2)*(1.d0+exp(-xi*dot_product(R1,kpt_latt(ik,:))))
          Hhf(jso,iso,ik) = Hhf(jso,iso,ik) + conjg(lgr(2)*(1.d0+exp(-xi*dot_product(R1,kpt_latt(ik,:)))))
          
          !+- x_iter(13)
          iorb=4; iso=(ispin-1)*Norb+iorb          
          jorb=6; jso=(ispin-1)*Norb+jorb          
          Hhf(iso,jso,ik) = Hhf(iso,jso,ik) + lgr(3)*(1.d0+exp(+xi*dot_product(R1,kpt_latt(ik,:))))
          Hhf(jso,iso,ik) = Hhf(jso,iso,ik) + conjg(lgr(3)*(1.d0+exp(+xi*dot_product(R1,kpt_latt(ik,:)))))

          !+- x_iter(14)
          iorb=5; iso=(ispin-1)*Norb+iorb          
          jorb=6; jso=(ispin-1)*Norb+jorb          
          Hhf(iso,jso,ik) = Hhf(iso,jso,ik) + lgr(4)*(1.d0+exp(+xi*dot_product(R1,kpt_latt(ik,:))))
          Hhf(jso,iso,ik) = Hhf(jso,iso,ik) + conjg(lgr(4)*(1.d0+exp(+xi*dot_product(R1,kpt_latt(ik,:)))))
          
       end do
    end do    
    !
    mu=0.d0
    call fix_mu(HHf,deltahf,mu) ; mu_fix=mu
    !
    deltahfr=0.d0
    do i=1,3!
       ir=ixr(i)
       do ik=1,Lk
          deltahfr(:,:,ir)=deltahfr(:,:,ir) + &
               deltahf(:,:,ik)*exp(xi*dot_product(rpt_latt(ir,:),kpt_latt(ik,:)))*wtk(ik)
       end do
    end do
    !
    xhf(1) = deltahfr(1,1,ir0)+deltahfr(1+Norb,1+Norb,ir0)
    xhf(2) = deltahfr(2,2,ir0)+deltahfr(2+Norb,2+Norb,ir0)
    xhf(3) = deltahfr(3,3,ir0)+deltahfr(3+Norb,3+Norb,ir0)
    !
    xhf(4) = deltahfr(1,3,ir0)
    xhf(5) = deltahfr(2,3,ir0)
    !
    xhf(6) = deltahfr(1,3,ir0)+deltahfr(1,3,irR)
    xhf(7) = deltahfr(2,3,ir0)+deltahfr(2,3,irR)
    !
    xhf(8) = deltahfr(4,4,ir0)+deltahfr(4+Norb,4+Norb,ir0)
    xhf(9) = deltahfr(5,5,ir0)+deltahfr(5+Norb,5+Norb,ir0)
    xhf(10) = deltahfr(6,6,ir0)+deltahfr(6+Norb,6+Norb,ir0)
    !
    xhf(11) = deltahfr(4,6,ir0)
    xhf(12) = deltahfr(5,6,ir0)
    !
    xhf(13) = deltahfr(4,6,ir0)+deltahfr(4,6,irL)
    xhf(14) = deltahfr(5,6,ir0)+deltahfr(5,6,irL)
    !
    delta_op(1) = xhf(6)  - phi_target(1)
    delta_op(2) = xhf(7)  - phi_target(2)
    delta_op(3) = xhf(13) - phi_target(3)
    delta_op(4) = xhf(14) - phi_target(4)
    !
    do i=1,4
       xdelta_op(i) = dreal(delta_op(i))
       xdelta_op(i+4) = dimag(delta_op(i))
    end do

    !x_iter=xhf
    ! xdelta_op(1) = dreal(delta_op(i))
    ! xdelta_op(2) = dimag(delta_op(i))
    
    !write(*,'(20F18.10)') delta_op!,xlgr
    ! write(111,'(10F18.10)') xhf(6),phi_target(1)
    ! write(112,'(10F18.10)') xhf(7),phi_target(2)
    ! write(113,'(10F18.10)') xhf(13),phi_target(3)
    ! write(114,'(10F18.10)') xhf(14),phi_target(4)
    ! write(115,'(10F18.10)') xlgr


    !write(*,'(20F18.10)') xhf(6),xhf(7),xhf(13),xhf(14)
    
  end function fix_lgr_params



  function min_fix_lgr_params(xlgr) result(xdelta_op)
    real(8),dimension(:) :: xlgr
    real(8) :: xdelta_op
    complex(8),dimension(14) :: xhf
    complex(8),dimension(:,:,:),allocatable :: Hhf
    complex(8),dimension(:,:,:),allocatable :: deltahfr,deltahf
    real(8)   :: mu
    complex(8),dimension(4) :: lgr,delta_op
    integer :: i,ir,ik,iorb,jorb,ispin,iso,jso
    !
    if(size(xlgr).ne.2) stop "wrong sized in xlgr"
    ! do i=1,4
    !    lgr(i) = xlgr(i)+xi*xlgr(i+4)
    ! end do
    lgr(1) = xlgr(1)+xi*xlgr(2)
    lgr(2) = -lgr(2)
    lgr(3) = -lgr(3)
    lgr(4) = lgr(1)
    
    allocate(Hhf(Nso,Nso,Lk),deltahf(Nso,Nso,Lk),deltahfr(Nso,Nso,nrpts))
    xhf=x_iter
    !
    HHf = HF_hamiltonian(xhf)
    HHf = HHf+Hk_toy
    !
    !+-  here add the part of the lgr multipliers to the operator of the order parameter -+!
    do ik=1,Lk
       do ispin=1,Nspin
          
          !+- x_iter(6)
          iorb=1; iso=(ispin-1)*Norb+iorb          
          jorb=3; jso=(ispin-1)*Norb+jorb          
          Hhf(iso,jso,ik) = Hhf(iso,jso,ik) + lgr(1)*(1.d0+exp(-xi*dot_product(R1,kpt_latt(ik,:))))
          Hhf(jso,iso,ik) = Hhf(jso,iso,ik) + conjg(lgr(1)*(1.d0+exp(-xi*dot_product(R1,kpt_latt(ik,:)))))
          
          !+- x_iter(7)
          iorb=2; iso=(ispin-1)*Norb+iorb          
          jorb=3; jso=(ispin-1)*Norb+jorb          
          Hhf(iso,jso,ik) = Hhf(iso,jso,ik) + lgr(2)*(1.d0+exp(-xi*dot_product(R1,kpt_latt(ik,:))))
          Hhf(jso,iso,ik) = Hhf(jso,iso,ik) + conjg(lgr(2)*(1.d0+exp(-xi*dot_product(R1,kpt_latt(ik,:)))))
          
          !+- x_iter(13)
          iorb=4; iso=(ispin-1)*Norb+iorb          
          jorb=6; jso=(ispin-1)*Norb+jorb          
          Hhf(iso,jso,ik) = Hhf(iso,jso,ik) + lgr(3)*(1.d0+exp(xi*dot_product(R1,kpt_latt(ik,:))))
          Hhf(jso,iso,ik) = Hhf(jso,iso,ik) + conjg(lgr(3)*(1.d0+exp(xi*dot_product(R1,kpt_latt(ik,:)))))

          !+- x_iter(14)
          iorb=5; iso=(ispin-1)*Norb+iorb          
          jorb=6; jso=(ispin-1)*Norb+jorb          
          Hhf(iso,jso,ik) = Hhf(iso,jso,ik) + lgr(4)*(1.d0+exp(xi*dot_product(R1,kpt_latt(ik,:))))
          Hhf(jso,iso,ik) = Hhf(jso,iso,ik) + conjg(lgr(4)*(1.d0+exp(xi*dot_product(R1,kpt_latt(ik,:)))))
          
       end do
    end do    
    !
    call fix_mu(HHf,deltahf,mu) ; mu_fix=mu
    !
    do i=1,3!
       ir=ixr(i)
       deltahfr(:,:,ir)=0.d0
       do ik=1,Lk
          deltahfr(:,:,ir)=deltahfr(:,:,ir) + &
               deltahf(:,:,ik)*exp(xi*dot_product(rpt_latt(ir,:),kpt_latt(ik,:)))*wtk(ik)
       end do
    end do
    !

    xhf(1) = deltahfr(1,1,ir0)+deltahfr(1+Norb,1+Norb,ir0)
    xhf(2) = deltahfr(2,2,ir0)+deltahfr(2+Norb,2+Norb,ir0)
    xhf(3) = deltahfr(3,3,ir0)+deltahfr(3+Norb,3+Norb,ir0)
    !
    xhf(4) = deltahfr(1,3,ir0)
    xhf(5) = deltahfr(2,3,ir0)
    !
    xhf(6) = deltahfr(1,3,ir0)+deltahfr(1,3,irR)
    xhf(7) = deltahfr(2,3,ir0)+deltahfr(2,3,irR)
    !
    xhf(8) = deltahfr(4,4,ir0)+deltahfr(4+Norb,4+Norb,ir0)
    xhf(9) = deltahfr(5,5,ir0)+deltahfr(5+Norb,5+Norb,ir0)
    xhf(10) = deltahfr(6,6,ir0)+deltahfr(6+Norb,6+Norb,ir0)
    !
    xhf(11) = deltahfr(4,6,ir0)
    xhf(12) = deltahfr(5,6,ir0)
    !
    xhf(13) = deltahfr(4,6,ir0)+deltahfr(4,6,irL)
    xhf(14) = deltahfr(5,6,ir0)+deltahfr(5,6,irL)
    !
    delta_op(1) = xhf(6)  - phi_target(1)
    delta_op(2) = xhf(7)  - phi_target(2)
    delta_op(3) = xhf(13) - phi_target(3)
    delta_op(4) = xhf(14) - phi_target(4)
    !
    xdelta_op=0.d0
    do i=1,4
       xdelta_op=xdelta_op+abs(delta_op(i))**2.d0
    end do

    write(*,'(20F18.10)') xdelta_op,xhf(6),xhf(7),xhf(13),xhf(14)
    
  end function min_fix_lgr_params


  


  function root_find_HF(x) result(out_x)
    implicit none
    real(8),dimension(:) :: x
    real(8),dimension(size(x)) :: out_x
    complex(8),dimension(14) :: xhf,xhf_
    complex(8),dimension(Nso,Nso,Lk) :: H_Hf
    complex(8),dimension(Nso,Nso,Lk) :: delta_hf
    complex(8),dimension(Nso,Nso,nrpts) :: delta_hfr
    integer :: i,j
    real(8) :: mu_fix

    if(size(x).ne.14) stop 'delta HF s(x)/= 14'
    !
    do i=1,14
       xhf(i) = x(i) + xi*0.d0
    end do
    !
    H_Hf=HF_hamiltonian(xhf)
    H_Hf=H_Hf+Hk_toy    
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
    !
    xhf_(1) = delta_hfr(1,1,ir0)+delta_hfr(1+Norb,1+Norb,ir0)
    xhf_(2) = delta_hfr(2,2,ir0)+delta_hfr(2+Norb,2+Norb,ir0)
    xhf_(3) = delta_hfr(3,3,ir0)+delta_hfr(3+Norb,3+Norb,ir0)
    !
    xhf_(4) = delta_hfr(1,3,ir0)
    xhf_(5) = delta_hfr(2,3,ir0)
    !
    xhf_(6) = delta_hfr(1,3,ir0) + delta_hfr(1,3,irR)
    xhf_(7) = delta_hfr(2,3,ir0) + delta_hfr(2,3,irR)
    !

    xhf_(8) = delta_hfr(4,4,ir0)+delta_hfr(4+Norb,4+Norb,ir0)
    xhf_(9) = delta_hfr(5,5,ir0)+delta_hfr(5+Norb,5+Norb,ir0)
    xhf_(10) = delta_hfr(6,6,ir0)+delta_hfr(6+Norb,6+Norb,ir0)
    !
    xhf_(11) = delta_hfr(4,6,ir0)
    xhf_(12) = delta_hfr(5,6,ir0)
    !
    xhf_(13) = delta_hfr(4,6,ir0) + delta_hfr(1,3,irL)
    xhf_(14) = delta_hfr(5,6,ir0) + delta_hfr(2,3,irL)
    xhf_=xhf_-xhf
    do i=1,14
       out_x(i)=dreal(xhf_(i))
    end do
    !
  end function root_find_HF

  !
  function HF_hamiltonian(x_iter) result(Hhf)
    implicit none
    complex(8) :: x_iter(14)
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
          iorb=1; iso=(ispin-1)*Norb+iorb          
          jorb=3; jso=(ispin-1)*Norb+jorb
          Hhf(jso,iso,ik) = -Vcell*x_iter(4)*(1.d0-exp(-xi*dot_product(R1,kpt_latt(ik,:))))
          Hhf(jso,iso,ik) = Hhf(jso,iso,ik) - Vcell*x_iter(6)*exp(-xi*dot_product(R1,kpt_latt(ik,:)))
          Hhf(iso,jso,ik) = conjg(Hhf(jso,iso,ik))
          !
          iorb=2; iso=(ispin-1)*Norb+iorb          
          jorb=3; jso=(ispin-1)*Norb+jorb
          Hhf(jso,iso,ik) = -Vcell*x_iter(5)*(1.d0-exp(-xi*dot_product(R1,kpt_latt(ik,:))))
          Hhf(jso,iso,ik) = Hhf(jso,iso,ik) - Vcell*x_iter(7)*exp(-xi*dot_product(R1,kpt_latt(ik,:)))
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
          Hhf(iso,jso,ik) = conjg(Hhf(jso,iso,ik))
          !
          iorb=5; iso=(ispin-1)*Norb+iorb          
          jorb=6; jso=(ispin-1)*Norb+jorb
          Hhf(jso,iso,ik) = -Vcell*x_iter(12)*(1.d0-exp(xi*dot_product(R1,kpt_latt(ik,:))))
          Hhf(jso,iso,ik) = Hhf(jso,iso,ik) - Vcell*x_iter(14)*exp(xi*dot_product(R1,kpt_latt(ik,:)))
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

    nrpts=size(irvec,1)    
    Hk=0.d0
    do ir=1,nrpts
       Rlat=irvec(ir,1)*R1+irvec(ir,2)*R2+irvec(ir,3)*R3
       dotRk=dot_product(Rlat,kpoint)
       Hk=Hk+Hr_w90(:,:,ir)*exp(xi*dotRK)/dble(ndegen(ir))
    end do

  end subroutine get_Hk_w90


  function get_HK(kpoint,N,Hk_in) result(Hk_out)
    implicit none
    real(8),dimension(:) :: kpoint
    integer              :: N
    complex(8),dimension(N,N) :: Hk_out
    complex(8),dimension(N,N,Nk_x,Nk_y,Nk_z) :: Hk_in
    integer :: io,jo,Nint_,ik_target
    complex(8),dimension(:,:,:,:,:),allocatable :: Hk_inter
    real(8),dimension(:,:,:),allocatable :: Hkr
    real(8),dimension(:),allocatable :: kix,kiy,kiz
    real(8) :: Hr,dH,delta_k,diffk
    integer :: i0,j0,k0
    integer :: i,j,k,ii,jj,kk
    integer :: ix,ik
    real(8),dimension(:,:),allocatable :: kpt_int
    real(8),dimension(3) :: ktarget,ktmp
    integer :: Lkint,ikint
    real(8),dimension(:),allocatable :: yout,yout_,hk_tmp
    real(8),dimension(:),allocatable :: grid_aux

    Nint_=2*Nint+1
    allocate(Hk_inter(Nso,Nso,Nint_,Nint_,Nint_))
    allocate(kix(Nint_),kiy(Nint_),kiz(Nint_))


    if(size(kpoint).ne.3) stop "kpoint.ne.3"    
    if(N.ne.Nso)  stop "N/=Nso"

    !+- find the closest point to the target one -+!
    diffk=100.d0
    ik_target=1
    do ik=1,Lk
       delta_k=0.d0
       delta_k=delta_k+(kpoint(1)-kpt_latt(ik,1))**2.d0
       delta_k=delta_k+(kpoint(2)-kpt_latt(ik,2))**2.d0
       delta_k=delta_k+(kpoint(3)-kpt_latt(ik,3))**2.d0
       if(delta_k.le.diffk) then
          ik_target=ik
          diffk=delta_k
       end if
    end do
    !
    ktarget=matmul(Bkinv,kpoint)
    !
    !+- get the *centered* index of the aux lattice -+!
    i0=ik_stride(ik_target,1)+Nk_x
    j0=ik_stride(ik_target,2)+Nk_y
    k0=ik_stride(ik_target,3)+Nk_z
    !
    !
    Lkint=Nint_*Nint_*Nint_
    allocate(kpt_int(Lkint,3))

    ikint=0
    do i=1,Nint_       
       ii=i0-(Nint+1)+i
       kix(i)= kxgrid_aux(ii)
       do j=1,Nint_
          jj=j0-(Nint+1)+j
          kiy(j)= kygrid_aux(jj)
          do k=1,Nint_
             kk=k0-(Nint+1)+k
             kiz(k)=kzgrid_aux(kk)
             ikint=ikint+1
             !
             kpt_int(ikint,:)=kix(i)*Bk1+ kiy(j)*Bk2 + kiz(k)*Bk3
             !
             Hk_inter(:,:,i,j,k)=Hk_in(:,:,ix_aux(ii),iy_aux(jj),iz_aux(kk)) 
             !
          end do
       end do
    end do
    !

    allocate(Hkr(Nint_,Nint_,Nint_))    

    Hk_out=0.d0
    do io=1,Nso
       do jo=1,Nso

          !+-get Hkr to be interpolated 
          Hkr=dreal(Hk_inter(io,jo,:,:,:))
          Hr=0.d0
          call polin3(kix,kiy,kiz,Hkr,ktarget(1),ktarget(2),ktarget(3),Hr,dH)          
          Hk_out(io,jo)=Hk_out(io,jo) + Hr
          !
          Hkr=dimag(Hk_inter(io,jo,:,:,:))
          Hr=0.d0
          call polin3(kix,kiy,kiz,Hkr,ktarget(1),ktarget(2),ktarget(3),Hr,dH)          
          Hk_out(io,jo)=Hk_out(io,jo) + xi*Hr          
          !
       end do
    end do
    !
    write(*,'(20F10.5)') kpoint,kpt_latt(ik_target,:),ktarget,Hk_w90(1,1,ik_target),Hk_out(1,1)
    !write(*,*) 'test_inter polin3',Hk_out(1,1)
    !write(*,*) Hk_out(1,1),Hk_out(2,2)
  end function Get_HK
  
  subroutine build_mp_grid_1d(Nx,Ny,Nz)
    implicit none
    integer :: Nx,Ny,Nz
    integer :: Nx_,Ny_,Nz_
    integer :: i,ix,iy,ik,iz,unitk
    real(8) :: kmax,kxmax,kymax,kzmax,dx,dy,dz
    real(8),dimension(:),allocatable :: kxr,kyr,kzr
    !
    !    
    if(allocated(kxgrid)) deallocate(kxgrid)
    allocate(kxgrid(Nx))
    kmax=dble(Nx-1)*0.5d0/dble(Nx)
    kxgrid = linspace(-kmax,kmax,Nx,mesh=dx)
    kxmax=kmax+dx*0.5d0
    !
    if(mod(Nx,2).eq.0) then
       kxgrid=kxgrid-0.5d0/dble(Nx)
    end if
    !
    unitk=free_unit()
    open(unit=unitk,file='mp_grid.out')
    do ix=1,Nx
       do iy=1,Ny
          do iz=1,Nz
             write(unitk,'(10F18.10)') kxgrid(ix)
          end do
       end do
    end do
    close(unitk)
  end subroutine build_mp_grid_1d


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

