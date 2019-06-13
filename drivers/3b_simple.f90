program officina
  !
  USE SCIFOR
  USE DMFT_TOOLS
  !
  USE VARS_GLOBAL
  USE DMFT_VECTORS
  USE HF
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
  complex(8),dimension(:,:,:),allocatable :: Hr_w90,Hr_w90_tmp
  real(8) :: mu_fix
  real(8) :: Uintra,Uinter,Delta_CF
  real(8) :: Uq,thop,err_hf
  !real(8),dimension(3,3) :: Bkinv

  real(8),dimension(3) :: ktest
  real(8),dimension(:),allocatable :: ek_out

  complex(8),dimension(:,:),allocatable :: Hsb

  complex(8),dimension(:),allocatable :: Uft

  complex(8),dimension(:,:),allocatable :: Hloc,Hktmp
  complex(8),dimension(:,:,:),allocatable :: Hk_w90,Hk_w90_tmp
  complex(8),dimension(:,:,:,:,:),allocatable :: Hk_w90_reshape,Hk_hf_reshape
  real(8),dimension(:,:),allocatable :: kpts,kpt_path
  integer,dimension(:,:),allocatable :: irvec2d,itmp
  integer,dimension(:),allocatable :: stride2D,stride2D_

  integer :: nr2d

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
  complex(8),dimension(7) :: x_iter,x_iter_
  complex(8),dimension(5) :: xtmp,xtmp_
  complex(8) :: xphi,xphi_
  complex(8),dimension(2) :: xpi,xpi_
  real(8),dimension(14) :: xR
  real(8),dimension(10) :: xr_tmp

  real(8) :: Ucell,Vcell
  real(8) :: Evalence,Econduction,tconduction,tvalence
  real(8) :: w0gap
  logical :: fix_phi
  !

  !+- START MPI -+!
  call init_MPI()
  comm = MPI_COMM_WORLD
  call StartMsg_MPI(comm)
  rank = get_Rank_MPI(comm)
  master = get_Master_MPI(comm)
  !

  !+- PARSE INPUT DRIVER -+!
  call parse_input_variable(Nk_x,"Nk_x","input.conf",default=12)
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


  call parse_input_variable(Econduction,"Econduction","input.conf",default=1.7d0)
  call parse_input_variable(Evalence,"Evalence","input.conf",default=-0.9d0)

  call parse_input_variable(tconduction,"tconduction","input.conf",default=-0.8d0)
  call parse_input_variable(tvalence,"tvalence","input.conf",default=0.4d0)


  call parse_input_variable(w0gap,"w0gap","input.conf",default=0.0d0)
  call parse_input_variable(hybloc,"hybloc","input.conf",default=0.d0)
  call parse_input_variable(fix_phi,"fix_phi","input.conf",default=.false.)
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
  call build_mp_grid_1d(Nk_x,Nk_y,Nk_z)
  !
  Lk=Nk_x
  Nk_y=1
  Nk_z=1
  ! R2=0.d0
  ! R3=0.d0
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
              kpt_latt(ik,idim) = kpt_latt(ik,idim) + kxgrid(i)*Bk1(idim)
           end do
           write(300,*) kpt_latt(ik,:)
        end do
     end do
  end do
  !


  file_name=reg(read_tns)//reg(file_w90_hr)  
  !+- read the w90 output -+!
  allocate(Hloc(Nso,Nso))
  call read_w90_hr(R1,R2,R3,Hr_w90,Hloc,irvec,ndegen,trim(file_name),1,6,1)
  nrpts=size(irvec,1)
  nrpts=Nk_x
  write(*,*) 'nrpts',nrpts
  deallocate(irvec); allocate(irvec(nrpts,3))
  deallocate(ndegen); allocate(ndegen(nrpts));ndegen=1

  allocate(rpt_latt(nrpts,3))
  irvec=0
  do ir=1,nrpts
     irvec(ir,1) = - (Nk_x-1)/2 + (ir-1)
     rpt_latt(ir,:)=irvec(ir,1)*R1
     modk=sqrt(rpt_latt(ir,1)**2.d0+rpt_latt(ir,2)**2.d0+rpt_latt(ir,3)**2.d0)
     if(modk.lt.1.d-12) ir0=ir
  end do
  deallocate(Hr_w90)
  R2=0.d0
  R3=0.d0
  Bk1=0.d0;Bk2=0.d0;Bk3=0.d0
  Bk1(1)=2.d0*pi/R1(1)
  !+- X-point
  kpath(1,:)=0.5d0*Bk1
  !+- G-point
  kpath(2,:)=0.d0
  !+- X-point
  kpath(3,:)=0.5d0*Bk1


  Norb=3
  Nspin=2
  Nso=Norb*Nspin
  allocate(Hk_w90(Nso,Nso,Lk),delta_hf(Nso,Nso,Lk),delta_hf_(Nso,Nso,Lk))
  allocate(tk(Norb))
  tk(1:2) = tconduction
  tk(3) = tvalence
  Hk_w90=0.d0
  do ik=1,Lk
     do ispin=1,Nspin
        !
        iorb=1
        iso=(ispin-1)*Norb+iorb
        Hk_w90(iso,iso,ik) = 2.d0*tk(iorb)*dcos(kpt_latt(ik,1)*R1(1)+kpt_latt(ik,2)*R1(2)+kpt_latt(ik,3)*R1(3))+Econduction
        iorb=2
        iso=(ispin-1)*Norb+iorb
        Hk_w90(iso,iso,ik) = 2.d0*tk(iorb)*dcos(kpt_latt(ik,1)*R1(1)+kpt_latt(ik,2)*R1(2)+kpt_latt(ik,3)*R1(3))+Econduction
        !
        iorb=3
        iso=(ispin-1)*Norb+iorb
        Hk_w90(iso,iso,ik) = 2.d0*tk(iorb)*dcos(kpt_latt(ik,1)*R1(1)+kpt_latt(ik,2)*R1(2)+kpt_latt(ik,3)*R1(3))+Evalence
        !
        iorb=1
        jorb=3
        iso=(ispin-1)*Norb+iorb
        jso=(ispin-1)*Norb+jorb
        Rlat=R1
        Hk_w90(iso,jso,ik) = hybloc*(1.d0-exp(xi*dot_product(Rlat,kpt_latt(ik,:)))) 
        Hk_w90(jso,iso,ik) = hybloc*(1.d0-exp(-xi*dot_product(Rlat,kpt_latt(ik,:))))
        !
        iorb=2
        jorb=3
        iso=(ispin-1)*Norb+iorb
        jso=(ispin-1)*Norb+jorb
        Rlat=R1
        Hk_w90(iso,jso,ik) = hybloc*(1.d0-exp(xi*dot_product(Rlat,kpt_latt(ik,:)))) 
        Hk_w90(jso,iso,ik) = hybloc*(1.d0-exp(-xi*dot_product(Rlat,kpt_latt(ik,:))))
        !
     end do
  end do
  !
  allocate(Hr_w90(Nso,Nso,nrpts))
  do ir=1,nrpts
     call FT_q2r(rpt_latt(ir,:),Hr_w90(:,:,ir),Hk_w90)
  end do
  !
  mu_fix=0.d0
  call fix_mu(Hk_w90,delta_hf,mu_fix)

  !
  allocate(delta_hfr(Nso,Nso,nrpts))
  do ir=1,nrpts
     !call FT_q2r(rpt_latt(ir,:),delta_hfr(:,:,ir),delta_hf) !+- be careful with FT; make a choice for the FT of operators!!
     !+- this choice is consistent w/ the interaction Hamiltonian written as (n_{R (1,2)} + n_{R+d (1,2)}) n_{R 3}
     delta_hfr(:,:,ir)=0.d0
     do ik=1,Lk
        delta_hfr(:,:,ir)=delta_hfr(:,:,ir) + &
             delta_hf(:,:,ik)*exp(xi*dot_product(rpt_latt(ir,:),kpt_latt(ik,:)))*wtk(ik)
     end do
  end do

  uio=free_unit()
  open(unit=uio,file='hyb_info.out')
  iso=0
  do i=1,Nso                
     do j=1,Nso
        iso=iso+1
        write(uio,*) i,j,iso+1
     end do
  end do

  close(uio)
  uio=free_unit()
  open(unit=uio,file='hyb_bareR.out')
  allocate(obs(Nso*Nso))
  do ir=1,Nrpts
     iso=0
     do i=1,Nso                
        do j=1,Nso
           iso=iso+1
           obs(iso) = delta_hfr(i,j,ir)
        end do
     end do
     write(uio,'(100F18.10)') rpt_latt(ir,1),dreal(obs(:)),dreal(obs(:))     
  end do
  close(uio)
  !
  uio=free_unit()
  open(unit=uio,file='k_hyb_bare.out')
  do ik=1,Lk
     write(uio,'(40F18.10)') kpt_latt(ik,1),dreal(delta_hf(:,:,ik)),dimag(delta_hf(:,:,ik))
  end do
  close(uio)
  !
  allocate(kpt_path(300,3))
  allocate(ek_out(Nso))
  allocate(Hktmp(Nso,Nso))
  !+- tmp
  ! do ik=1,Lk
  !    Hktmp=Hk_w90(:,:,ik)
  !    call eigh(Hktmp,ek_out)
  !    write(567,'(30F18.10)') ek_out-mu_fix     
  ! end do
  !+- tmp

  uio=free_unit()
  open(unit=uio,file='tns_bare_bands.out')
  !
  modk=0.d0
  do i=1,2
     delta_kpath=kpath(i+1,:)-kpath(i,:)
     do ik=1,100
        j=(i-1)*100 + ik
        kpt_path(j,:) = kpath(i,:) + dble(ik-1)/100.d0*delta_kpath
        modk=modk+sqrt(dot_product(1.d0/100.d0*delta_kpath,1.d0/100.d0*delta_kpath))
        !
        call FT_r2q(kpt_path(j,:),Hktmp,Hr_w90)
        !
        call eigh(Hktmp,ek_out)
        write(uio,'(30F18.10)') modk,ek_out-mu_fix
        !
     end do
     !
  end do
  close(uio)

  x_iter(1) = delta_hfr(1,1,ir0)+delta_hfr(1+Norb,1+Norb,ir0)
  x_iter(2) = delta_hfr(2,2,ir0)+delta_hfr(2+Norb,2+Norb,ir0)
  x_iter(3) = delta_hfr(3,3,ir0)+delta_hfr(3+Norb,3+Norb,ir0)
  !
  x_iter(4) = delta_hfr(1,3,ir0)
  x_iter(5) = delta_hfr(2,3,ir0)
  !
  x_iter(6) = delta_hfr(1,3,ir0)+delta_hfr(1,3,ir0+1)
  x_iter(7) = delta_hfr(2,3,ir0)+delta_hfr(2,3,ir0+1)
  !
  write(*,*) "HF optimization: bare parameters"
  write(*,*) x_iter(1)
  write(*,*) x_iter(2)
  write(*,*) x_iter(3)
  write(*,*) x_iter(4)
  write(*,*) x_iter(5)
  write(*,*) x_iter(6)
  write(*,*) x_iter(7)
  write(*,*)
  write(*,*)
  !
  x_iter(1) = 0.d0
  x_iter(2) = 0.d0
  x_iter(3) = 2.d0
  !
  x_iter(4) = 0.d0
  x_iter(5) = 0.d0
  !
  x_iter(6) = 0.5
  x_iter(7) = 0.5
  write(*,*) "HF optimization: initial paramteres"
  write(*,*) x_iter(1)
  write(*,*) x_iter(2)
  write(*,*) x_iter(3)
  write(*,*) x_iter(4)
  write(*,*) x_iter(5)
  write(*,*) x_iter(6)
  write(*,*) x_iter(7)
  write(*,*)
  write(*,*)

  allocate(H_Hf(Nso,Nso,Lk))

  xphi=0.50 
  x_iter(6:7) = xphi
  xtmp=x_iter(1:5)  
  do i=1,5
     xr_tmp(i) = dreal(xtmp(i))
     xr_tmp(i+5) = dimag(xtmp(i))
  end do

  !+- TO DO THINGS -+!

  if(fix_phi) then
     uio=free_unit()
     open(unit=uio,file='loop_fixed_order_parameter.out')
     unit_in=free_unit()
     open(unit=unit_in,file='bands_VS_order_parameter.out')
     xphi=-0.01d0
     do ihf=1,50
        !
        xphi=xphi+0.01d0
        write(*,*) 'fixed phi loop',ihf,xphi
        x_iter(6:7)=xphi
        do jhf=1,Nhf_
           x_iter_=x_iter
           !
           H_Hf=HF_hamiltonian(x_iter)
           H_Hf=H_Hf+Hk_w90
           !
           call fix_mu(H_Hf,delta_hf,mu_fix)
           !
           do ir=ir0,ir0+1
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
           x_iter=x_iter*wmix+(1.d0-wmix)*x_iter_             
        end do
        xtmp=x_iter(1:5)  
        do i=1,5
           xr_tmp(i) = dreal(xtmp(i))
           xr_tmp(i+5) = dimag(xtmp(i))
        end do
        call fsolve(root_find_inner_loop,xr_tmp,tol=1.d-10)
        do i=1,5
           xtmp(i) = xr_tmp(i)+xi*xr_tmp(i+5)
        end do
        x_iter(1:5)=xtmp;x_iter(6:7)=xphi
        H_Hf=HF_hamiltonian(x_iter)
        H_Hf=H_Hf+Hk_w90
        call fix_mu(H_Hf,delta_hf,mu_fix,eout)
        !+- plot bands for the fixed value of the order parameter -+!
        do ir=1,nrpts
           call FT_q2r(rpt_latt(ir,:),Hr_w90(:,:,ir),H_hf)
        end do
        modk=0.d0
        do i=1,2
           delta_kpath=kpath(i+1,:)-kpath(i,:)
           do ik=1,100
              j=(i-1)*100 + ik
              kpt_path(j,:) = kpath(i,:) + dble(ik-1)/100.d0*delta_kpath
              modk=modk+sqrt(dot_product(1.d0/100.d0*delta_kpath,1.d0/100.d0*delta_kpath))
              !
              call FT_r2q(kpt_path(j,:),Hktmp,Hr_w90)
              !
              call eigh(Hktmp,ek_out)
              write(unit_in,'(30F18.10)') modk,ek_out-mu_fix
              !
           end do
           !
        end do
        write(unit_in,*)
        write(unit_in,*)          
        !
        !
        !+- double counting term -+!
        Eout=Eout-Ucell*0.25d0*(dreal(x_iter(1))**2.d0+dreal(x_iter(2))**2.d0+dreal(x_iter(3))**2.d0)
        Eout=Eout-2*Vcell*(dreal(x_iter(1))*dreal(x_iter(3))+dreal(x_iter(2))*dreal(x_iter(3)))
        !
        Eout=Eout+Vcell*abs(x_iter(4))**2.d0+abs(x_iter(5))**2.d0
        Eout=Eout+Vcell*abs(x_iter(6)-x_iter(4))**2.d0+abs(x_iter(7)-x_iter(5))**2.d0
        !
        write(uio, '(10F18.10)') dreal(xphi),dreal(delta_hfr(1,3,ir0)+delta_hfr(1,3,ir0+1)),dreal(x_iter(1:5)),Eout+mu_fix*(dreal(x_iter(1))+dreal(x_iter(2))+dreal(x_iter(3)))
        !
     end do
     close(uio)
     close(unit_in)     
  end if
  !
  xphi=0.5d0
  !
  x_iter(1) = 0.d0
  x_iter(2) = 0.d0
  x_iter(3) = 2.d0
  !
  x_iter(4) = 0.1d0
  x_iter(5) = 0.1d0
  !
  !
  uio=free_unit()
  open(unit=uio,file='loop_phi_12symm.out')
  !
  do ihf=1,Nhf
     !
     xphi_=xphi
     x_iter(6:7)= xphi
     !x_iter(7)=-xphi
     do jhf=1,Nhf_
        !
        write(*,*) 'self-consistency phi loop',ihf,xphi
        x_iter_=x_iter
        !
        H_Hf=HF_hamiltonian_symm(x_iter)
        H_Hf=H_Hf+Hk_w90
        !
        call fix_mu(H_Hf,delta_hf,mu_fix)
        !
        do ir=ir0,ir0+1
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
        ! write(235,'(10F18.10)') x_iter(4),x_iter(5)
        !
        x_iter=x_iter*wmix+(1.d0-wmix)*x_iter_             
        !
     end do
     xtmp=x_iter(1:5)  
     do i=1,5
        xr_tmp(i) = dreal(xtmp(i))
        xr_tmp(i+5) = dimag(xtmp(i))
     end do
     call fsolve(root_find_inner_loop_symm,xr_tmp,tol=1.d-10)
     do i=1,5
        xtmp(i)=xr_tmp(i)+xi*xr_tmp(i+5)
     end do
     x_iter(1:5)=xtmp;
     x_iter(6:7) =  xphi
     !x_iter(7) = -xphi
     !
     H_Hf=HF_hamiltonian_symm(x_iter)
     H_Hf=H_Hf+Hk_w90
     !
     call fix_mu(H_Hf,delta_hf,mu_fix,eout)
     do ir=ir0,ir0+1
        delta_hfr(:,:,ir)=0.d0
        do ik=1,Lk
           delta_hfr(:,:,ir)=delta_hfr(:,:,ir) + &
                delta_hf(:,:,ik)*exp(xi*dot_product(rpt_latt(ir,:),kpt_latt(ik,:)))*wtk(ik)
        end do
     end do
     xphi=dreal(delta_hfr(1,3,ir0)+delta_hfr(1,3,ir0+1))
     xphi=wmix*xphi+(1-wmix)*xphi_
     xphi_=xphi
     !
     !+- double counting term -+!
     Eout=Eout-Ucell*0.25d0*(dreal(x_iter(1))**2.d0+dreal(x_iter(2))**2.d0+dreal(x_iter(3))**2.d0)
     Eout=Eout-2*Vcell*(dreal(x_iter(1))*dreal(x_iter(3))+dreal(x_iter(2))*dreal(x_iter(3)))
     !
     Eout=Eout+Vcell*abs(x_iter(4))**2.d0+abs(x_iter(5))**2.d0
     Eout=Eout+Vcell*abs(x_iter(6)-x_iter(4))**2.d0+abs(x_iter(7)-x_iter(5))**2.d0
     !
     write(uio, '(20F18.10)') dreal(xphi),dreal(x_iter(1:5)),Eout+mu_fix*(dreal(x_iter(1))+dreal(x_iter(2))+dreal(x_iter(3))), &
          dreal(delta_hfr(1,3,ir0)+delta_hfr(1,3,ir0+1)),dreal(delta_hfr(2,3,ir0)+delta_hfr(2,3,ir0+1))
     !
     !
  end do
  close(uio)
  
  H_Hf=HF_hamiltonian(x_iter)
  H_Hf=H_Hf+Hk_w90
  call fix_mu(H_Hf,delta_hf,mu_fix,eout)

  !+- plot real space hybridizations -+!
  uio=free_unit()
  open(unit=uio,file='hyb_TNS_VS_r_12symm.out')
  !allocate(obs(Nso*Nso))
  obs=0.d0
  do ir=1,nrpts
     delta_hfr(:,:,ir)=0.d0
     do ik=1,Lk
        delta_hfr(:,:,ir)=delta_hfr(:,:,ir) + &
             delta_hf(:,:,ik)*exp(xi*dot_product(rpt_latt(ir,:),kpt_latt(ik,:)))*wtk(ik)
     end do
  end do
  do ir=1,Nrpts
     iso=0
     do i=1,Nso                
        do j=1,Nso
           iso=iso+1
           obs(iso) = delta_hfr(i,j,ir)
        end do
     end do
     write(uio,'(100F18.10)') rpt_latt(ir,1),dreal(obs(:)),dreal(obs(:))     
  end do
  close(uio)

  
  xphi=0.5d0
  !
  x_iter(1) = 0.d0
  x_iter(2) = 0.d0
  x_iter(3) = 2.d0
  !
  x_iter(4) = 0.d0
  x_iter(5) = x_iter(4)+xphi
  !
  !
  uio=free_unit()
  open(unit=uio,file='loop_phi_12minus.out')
  !
  do ihf=1,Nhf
     !
     xphi_=xphi
     x_iter(6)= xphi
     x_iter(7)=-xphi
     do jhf=1,Nhf_
        !
        write(*,*) 'self-consistency phi loop',ihf,xphi
        x_iter_=x_iter
        !
        H_Hf=HF_hamiltonian(x_iter)
        H_Hf=H_Hf+Hk_w90
        !
        call fix_mu(H_Hf,delta_hf,mu_fix)
        !
        do ir=ir0,ir0+1
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
        ! write(235,'(10F18.10)') x_iter(4),x_iter(5)
        !
        x_iter=x_iter*wmix+(1.d0-wmix)*x_iter_             
        !
     end do
     xtmp=x_iter(1:5)  
     do i=1,5
        xr_tmp(i) = dreal(xtmp(i))
        xr_tmp(i+5) = dimag(xtmp(i))
     end do
     call fsolve(root_find_inner_loop,xr_tmp,tol=1.d-10)
     do i=1,5
        xtmp(i)=xr_tmp(i)+xi*xr_tmp(i+5)
     end do
     x_iter(1:5)=xtmp;
     x_iter(6) =  xphi
     x_iter(7) = -xphi
     !
     H_Hf=HF_hamiltonian(x_iter)
     H_Hf=H_Hf+Hk_w90
     !
     call fix_mu(H_Hf,delta_hf,mu_fix,eout)
     do ir=ir0,ir0+1
        delta_hfr(:,:,ir)=0.d0
        do ik=1,Lk
           delta_hfr(:,:,ir)=delta_hfr(:,:,ir) + &
                delta_hf(:,:,ik)*exp(xi*dot_product(rpt_latt(ir,:),kpt_latt(ik,:)))*wtk(ik)
        end do
     end do
     xphi=dreal(delta_hfr(1,3,ir0)+delta_hfr(1,3,ir0+1))
     xphi=wmix*xphi+(1-wmix)*xphi_
     xphi_=xphi
     !
     !+- double counting term -+!
     Eout=Eout-Ucell*0.25d0*(dreal(x_iter(1))**2.d0+dreal(x_iter(2))**2.d0+dreal(x_iter(3))**2.d0)
     Eout=Eout-2*Vcell*(dreal(x_iter(1))*dreal(x_iter(3))+dreal(x_iter(2))*dreal(x_iter(3)))
     !
     Eout=Eout+Vcell*abs(x_iter(4))**2.d0+abs(x_iter(5))**2.d0
     Eout=Eout+Vcell*abs(x_iter(6)-x_iter(4))**2.d0+abs(x_iter(7)-x_iter(5))**2.d0
     !
     write(uio, '(10F18.10)') dreal(xphi),dreal(x_iter(1:5)),Eout+mu_fix*(dreal(x_iter(1))+dreal(x_iter(2))+dreal(x_iter(3))), &
          dreal(delta_hfr(1,3,ir0)+delta_hfr(1,3,ir0+1)),dreal(delta_hfr(2,3,ir0)+delta_hfr(2,3,ir0+1))
     !
     !
  end do
  close(uio)




  !
  !
  H_Hf=HF_hamiltonian(x_iter)
  H_Hf=H_Hf+Hk_w90
  call fix_mu(H_Hf,delta_hf,mu_fix,eout)


  !+- plot real space hybridizations -+!
  uio=free_unit()
  open(unit=uio,file='hyb_TNS_VS_r_12minus.out')
  !allocate(obs(Nso*Nso))
  obs=0.d0
  do ir=1,nrpts
     delta_hfr(:,:,ir)=0.d0
     do ik=1,Lk
        delta_hfr(:,:,ir)=delta_hfr(:,:,ir) + &
             delta_hf(:,:,ik)*exp(xi*dot_product(rpt_latt(ir,:),kpt_latt(ik,:)))*wtk(ik)
     end do
  end do
  do ir=1,Nrpts
     iso=0
     do i=1,Nso                
        do j=1,Nso
           iso=iso+1
           obs(iso) = delta_hfr(i,j,ir)
        end do
     end do
     write(uio,'(100F18.10)') rpt_latt(ir,1),dreal(obs(:)),dreal(obs(:))     
  end do
  close(uio)










  xpi(1)=0.5d0
  xpi(2)=-0.2d0
  !
  x_iter(1) = 0.05d0 !+- Delta_{11} (0)
  x_iter(2) = 0.15d0 !+- Delta_{22} (0)
  x_iter(3) = 1.8d0  !+- Delta_{33} (0)
  !
  x_iter(4) = 0.20   !+- Delta_{13} (0)
  x_iter(5) = 0.1    !+- Delta_{23} (0)
  !
  uio=free_unit()
  open(unit=uio,file='loop_phi_12free.out')
  !
  do ihf=1,Nhf
     !
     xpi_=xpi
     x_iter(6)= xpi(1)
     x_iter(7)= xpi(2)
     !
     do jhf=1,Nhf_
        !
        write(*,*) 'self-consistency phi loop',ihf,xpi
        x_iter_=x_iter
        !
        H_Hf=HF_hamiltonian(x_iter)
        H_Hf=H_Hf+Hk_w90
        !
        call fix_mu(H_Hf,delta_hf,mu_fix)
        !
        do ir=ir0,ir0+1
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
        x_iter=x_iter*wmix+(1.d0-wmix)*x_iter_             
        !
     end do
     xtmp=x_iter(1:5)  
     do i=1,5
        xr_tmp(i) = dreal(xtmp(i))
        xr_tmp(i+5) = dimag(xtmp(i))
     end do
     call fsolve(root_find_inner_loop_free,xr_tmp,tol=1.d-10)
     do i=1,5
        xtmp(i)=xr_tmp(i)+xi*xr_tmp(i+5)
     end do
     x_iter(1:5)=xtmp;
     x_iter(6) =  xpi(1)
     x_iter(7) =  xpi(2)
     !
     H_Hf=HF_hamiltonian(x_iter)
     H_Hf=H_Hf+Hk_w90
     !
     call fix_mu(H_Hf,delta_hf,mu_fix,eout)
     do ir=ir0,ir0+1
        delta_hfr(:,:,ir)=0.d0
        do ik=1,Lk
           delta_hfr(:,:,ir)=delta_hfr(:,:,ir) + &
                delta_hf(:,:,ik)*exp(xi*dot_product(rpt_latt(ir,:),kpt_latt(ik,:)))*wtk(ik)
        end do
     end do
     !
     xpi(1)=dreal(delta_hfr(1,3,ir0)+delta_hfr(1,3,ir0+1))
     xpi(2)=dreal(delta_hfr(2,3,ir0)+delta_hfr(2,3,ir0+1))
     !
     xpi  = wmix*xpi+(1-wmix)*xpi_
     xpi_ = xpi
     !
     !+- double counting term -+!
     Eout=Eout-Ucell*0.25d0*(dreal(x_iter(1))**2.d0+dreal(x_iter(2))**2.d0+dreal(x_iter(3))**2.d0)
     Eout=Eout-2*Vcell*(dreal(x_iter(1))*dreal(x_iter(3))+dreal(x_iter(2))*dreal(x_iter(3)))
     !
     Eout=Eout+Vcell*abs(x_iter(4))**2.d0+abs(x_iter(5))**2.d0
     Eout=Eout+Vcell*abs(x_iter(6)-x_iter(4))**2.d0+abs(x_iter(7)-x_iter(5))**2.d0
     !
     write(uio, '(20F18.10)') dreal(xpi),dreal(x_iter(1:5)),Eout+mu_fix*(dreal(x_iter(1))+dreal(x_iter(2))+dreal(x_iter(3))), &
          dreal(delta_hfr(1,3,ir0)+delta_hfr(1,3,ir0+1)),dreal(delta_hfr(2,3,ir0)+delta_hfr(2,3,ir0+1))
     !
     !
  end do
  close(uio)




  !
  !
  H_Hf=HF_hamiltonian(x_iter)
  H_Hf=H_Hf+Hk_w90
  call fix_mu(H_Hf,delta_hf,mu_fix,eout)


  !+- plot real space hybridizations -+!
  uio=free_unit()
  open(unit=uio,file='hyb_TNS_VS_r_12free.out')
  !allocate(obs(Nso*Nso))
  obs=0.d0
  do ir=1,nrpts
     delta_hfr(:,:,ir)=0.d0
     do ik=1,Lk
        delta_hfr(:,:,ir)=delta_hfr(:,:,ir) + &
             delta_hf(:,:,ik)*exp(xi*dot_product(rpt_latt(ir,:),kpt_latt(ik,:)))*wtk(ik)
     end do
  end do
  do ir=1,Nrpts
     iso=0
     do i=1,Nso                
        do j=1,Nso
           iso=iso+1
           obs(iso) = delta_hfr(i,j,ir)
        end do
     end do
     write(uio,'(100F18.10)') rpt_latt(ir,1),dreal(obs(:)),dreal(obs(:))     
  end do
  close(uio)








  !+- plot bands for the fixed value of the order parameter -+!
  unit_in=free_unit()
  open(unit=unit_in,file='TNS_bands.out')
  do ir=1,nrpts
     call FT_q2r(rpt_latt(ir,:),Hr_w90(:,:,ir),H_hf)
  end do
  modk=0.d0
  do i=1,2
     delta_kpath=kpath(i+1,:)-kpath(i,:)
     do ik=1,100
        j=(i-1)*100 + ik
        kpt_path(j,:) = kpath(i,:) + dble(ik-1)/100.d0*delta_kpath
        modk=modk+sqrt(dot_product(1.d0/100.d0*delta_kpath,1.d0/100.d0*delta_kpath))
        !
        call FT_r2q(kpt_path(j,:),Hktmp,Hr_w90)
        !
        call eigh(Hktmp,ek_out)
        write(unit_in,'(30F18.10)') modk,ek_out-mu_fix
        !
     end do
     !
  end do
  close(unit_in)

  !+- plot final bands -+!

  ! xphi=-0.3d0
  ! do ihf=1,Nhf
  !    !
  !    xphi_=xphi
  !    x_iter(6:7)=xphi
  !    do jhf=1,Nhf
  !       x_iter_=x_iter
  !       !
  !       H_Hf=HF_hamiltonian(x_iter)
  !       H_Hf=H_Hf+Hk_w90
  !       !
  !       call fix_mu(H_Hf,delta_hf,mu_fix)
  !       !
  !       do ir=1,nrpts
  !          delta_hfr(:,:,ir)=0.d0
  !          do ik=1,Lk
  !             delta_hfr(:,:,ir)=delta_hfr(:,:,ir) + &
  !                  delta_hf(:,:,ik)*exp(xi*dot_product(rpt_latt(ir,:),kpt_latt(ik,:)))*wtk(ik)
  !          end do
  !       end do
  !       !
  !       x_iter(1) = delta_hfr(1,1,ir0)+delta_hfr(1+Norb,1+Norb,ir0)
  !       x_iter(2) = delta_hfr(2,2,ir0)+delta_hfr(2+Norb,2+Norb,ir0)
  !       x_iter(3) = delta_hfr(3,3,ir0)+delta_hfr(3+Norb,3+Norb,ir0)
  !       !
  !       x_iter(4) = delta_hfr(1,3,ir0)
  !       x_iter(5) = delta_hfr(2,3,ir0)
  !       !
  !       ! x_iter(6) = delta_hfr(1,3,ir0)+delta_hfr(1,3,ir0+1)
  !       ! x_iter(7) = delta_hfr(2,3,ir0)+delta_hfr(2,3,ir0+1)
  !       !
  !       !
  !       ! x_iter=x_iter*wmix+(1.d0-wmix)*x_iter_             
  !       ! H_Hf=HF_hamiltonian(x_iter)
  !       ! H_Hf=H_Hf+Hk_w90     
  !       ! call fix_mu(H_Hf,delta_hf,mu_fix,Eout)
  !       write(*,*) ihf,jhf
  !       !write(478,'(10F18.10)') dreal(x_iter(1:7))
  !    end do
  !    !
  !    xtmp=x_iter(1:5)  
  !    do i=1,5
  !       xr_tmp(i) = dreal(xtmp(i))
  !       xr_tmp(i+5) = dimag(xtmp(i))
  !    end do
  !    call fsolve(root_find_inner_loop,xr_tmp,tol=1.d-10)
  !    H_Hf=HF_hamiltonian(x_iter)
  !    H_Hf=H_Hf+Hk_w90
  !    !    
  !    call fix_mu(H_Hf,delta_hf,mu_fix,eout)

  !    do ir=1,nrpts
  !       delta_hfr(:,:,ir)=0.d0
  !       do ik=1,Lk
  !          delta_hfr(:,:,ir)=delta_hfr(:,:,ir) + &
  !               delta_hf(:,:,ik)*exp(xi*dot_product(rpt_latt(ir,:),kpt_latt(ik,:)))*wtk(ik)
  !       end do
  !    end do
  !    xphi=delta_hfr(1,3,ir0)+delta_hfr(1,3,ir0+1)
  !    xphi=xphi*wmix+(1.d0-wmix)*xphi_

  !    !
  !    ! Eout=0.d0
  !    ! do ik=1,Lk
  !    !    do iso=1,Nso
  !    !       do jso=1,Nso
  !    !          Eout=Eout+H_Hf(iso,jso,ik)*delta_hf(iso,jso,ik)*wtk(ik)
  !    !       end do
  !    !    end do
  !    ! end do

  !    !
  !    !+- double counting term -+!
  !    !write(480, '(10F18.10)') dreal(xphi),Eout,mu_fix
  !    Eout=Eout-Ucell*0.25d0*(dreal(x_iter(1))**2.d0+dreal(x_iter(2))**2.d0+dreal(x_iter(3))**2.d0)
  !    Eout=Eout-2*Vcell*(dreal(x_iter(1))*dreal(x_iter(3))+dreal(x_iter(2))*dreal(x_iter(3)))
  !    !
  !    Eout=Eout+Vcell*abs(x_iter(4))**2.d0+abs(x_iter(5))**2.d0
  !    Eout=Eout+Vcell*abs(x_iter(6)-x_iter(4))**2.d0+abs(x_iter(7)-x_iter(5))**2.d0
  !    !
  !    write(uio,'(10F18.10)') dreal(xphi),Eout,Eout+mu_fix*(dreal(x_iter(1))+dreal(x_iter(2))+dreal(x_iter(3)))
  !    !write(479, '(10F18.10)') dreal(xphi),dreal(x_iter)
  !    !
  ! end do
  ! close(uio)


  stop

  !
  uio=free_unit()
  open(unit=uio,file='hf_loop.out')
  !
  do ihf=1,Nhf
     write(uio,'(20F18.10)') dreal(x_iter),dimag(x_iter),mu_fix
     write(346,'(20F18.10)') xphi
     !
     xphi_ = xphi
     !
     do i=1,5
        xr_tmp(i) = dreal(xtmp(i))
        xr_tmp(i+5) = dimag(xtmp(i))
     end do
     call fsolve(root_find_inner_loop,xr_tmp,tol=1.d-10)
     do i=1,5
        xtmp(i)=xr_tmp(i)+xi*xr_tmp(i+5)
     end do
     x_iter(1:5)=xtmp;x_iter(6:7)=xphi
     H_Hf=HF_hamiltonian(x_iter)
     H_Hf=H_Hf+Hk_w90
     call fix_mu(H_Hf,delta_hf,mu_fix,eout)
     do ir=1,nrpts
        delta_hfr(:,:,ir)=0.d0
        do ik=1,Lk
           delta_hfr(:,:,ir)=delta_hfr(:,:,ir) + &
                delta_hf(:,:,ik)*exp(xi*dot_product(rpt_latt(ir,:),kpt_latt(ik,:)))*wtk(ik)
        end do
     end do
     ! !
     x_iter(1) = delta_hfr(1,1,ir0)+delta_hfr(1+Norb,1+Norb,ir0)
     x_iter(2) = delta_hfr(2,2,ir0)+delta_hfr(2+Norb,2+Norb,ir0)
     x_iter(3) = delta_hfr(3,3,ir0)+delta_hfr(3+Norb,3+Norb,ir0)
     !
     x_iter(4) = delta_hfr(1,3,ir0)
     x_iter(5) = delta_hfr(2,3,ir0)
     !
     xphi = delta_hfr(1,3,ir0)+delta_hfr(1,3,ir0+1);x_iter(6:7)=xphi


     Eout=Eout-Ucell*0.25d0*(dreal(x_iter(1))**2.d0+dreal(x_iter(2))**2.d0+dreal(x_iter(3))**2.d0)
     Eout=Eout-2*Vcell*(dreal(x_iter(1))*dreal(x_iter(3))+dreal(x_iter(2))*dreal(x_iter(3)))
     !
     Eout=Eout+Vcell*abs(x_iter(4))**2.d0+abs(x_iter(5))**2.d0
     Eout=Eout+Vcell*abs(x_iter(6)-x_iter(4))**2.d0+abs(x_iter(7)-x_iter(5))**2.d0


     write(345,'(20F18.10)') dreal(x_iter(6)),eout,eout+mu_fix*(dreal(x_iter(1))+dreal(x_iter(2))+dreal(x_iter(3)))
     !
     xphi=xphi*wmix+(1.d0-wmix)*xphi_     
     !
  end do

  close(uio)
  !stop

  x_iter(1) = delta_hfr(1,1,ir0)+delta_hfr(1+Norb,1+Norb,ir0)
  x_iter(2) = delta_hfr(2,2,ir0)+delta_hfr(2+Norb,2+Norb,ir0)
  x_iter(3) = delta_hfr(3,3,ir0)+delta_hfr(3+Norb,3+Norb,ir0)
  !
  x_iter(4) = delta_hfr(1,3,ir0)!-0.3!+0.1*xi
  x_iter(5) = delta_hfr(2,3,ir0)!-0.2!-0.1*xi
  !
  x_iter(6) = delta_hfr(1,3,ir0)+delta_hfr(1,3,ir0+1)+0.5d0
  x_iter(7) = delta_hfr(2,3,ir0)+delta_hfr(2,3,ir0+1)+0.3d0
  !
  write(*,*) x_iter(1)
  write(*,*) x_iter(2)
  write(*,*) x_iter(3)
  write(*,*) x_iter(4)
  write(*,*) x_iter(5)
  write(*,*) x_iter(6)
  write(*,*) x_iter(7)


  H_Hf=HF_hamiltonian(x_iter)
  H_Hf=H_Hf+Hk_w90

  !
  call fix_mu(H_Hf,delta_hf,mu_fix)
  iorb=1
  jorb=3
  ! do ik=1,Lk
  !    write(569,'(10F18.10)') kpt_latt(ik,1),delta_hf(iorb,jorb,ik)
  ! end do


  do ir=1,nrpts
     delta_hfr(:,:,ir)=0.d0
     do ik=1,Lk
        delta_hfr(:,:,ir)=delta_hfr(:,:,ir) + &
             delta_hf(:,:,ik)*exp(xi*dot_product(rpt_latt(ir,:),kpt_latt(ik,:)))*wtk(ik)
     end do
  end do

  ! x_iter(1) = delta_hfr(1,1,ir0)+delta_hfr(1+Norb,1+Norb,ir0)
  ! x_iter(2) = delta_hfr(2,2,ir0)+delta_hfr(2+Norb,2+Norb,ir0)
  ! x_iter(3) = delta_hfr(3,3,ir0)+delta_hfr(3+Norb,3+Norb,ir0)
  ! !
  ! x_iter(4) = delta_hfr(1,3,ir0)
  ! x_iter(5) = delta_hfr(2,3,ir0)
  ! !
  ! x_iter(6) = delta_hfr(1,3,ir0)+delta_hfr(1,3,ir0+1)
  ! x_iter(7) = delta_hfr(2,3,ir0)+delta_hfr(2,3,ir0+1)

  !+- full loop -+!
  uio=free_unit()
  open(unit=uio,file='full_hf_loop.out')
  write(uio,'(20F18.10)') dreal(x_iter),dimag(x_iter)
  !

  do ihf=1,Nhf

     write(uio,'(20F18.10)') dreal(x_iter),dimag(x_iter),mu_fix

     x_iter_=x_iter
     !
     H_Hf=HF_hamiltonian(x_iter)
     H_Hf=H_Hf+Hk_w90
     !
     call fix_mu(H_Hf,delta_hf,mu_fix)
     !
     do ir=1,nrpts
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
     x_iter(6) = delta_hfr(1,3,ir0)+delta_hfr(1,3,ir0+1)
     x_iter(7) = delta_hfr(2,3,ir0)+delta_hfr(2,3,ir0+1)
     !
     x_iter=x_iter*wmix+(1.d0-wmix)*x_iter_     
     !
  end do
  do i=1,7
     xr(i)=dreal(x_iter(i))
     xr(i+7)=dimag(x_iter(i))
  end do
  call fsolve(root_find_HF,xr,tol=1.d-10)
  write(uio,'(20F18.10)') xr
  close(uio)


  stop




  !




  stop

contains


  function root_find_inner_loop(x) result(out_x)
    implicit none
    real(8),dimension(:) :: x
    real(8),dimension(size(x)) :: out_x
    complex(8),dimension(7) :: xhf,xhf_
    complex(8),dimension(Nso,Nso,Lk) :: H_Hf
    complex(8),dimension(Nso,Nso,Lk) :: delta_hf
    complex(8),dimension(Nso,Nso,nrpts) :: delta_hfr
    integer :: i,j
    real(8) :: mu_fix

    if(size(x).ne.10) stop 'delta HF s(x)/= 10'
    !
    do i=1,5
       xhf(i) = x(i) + xi*x(i+5)
    end do
    xhf(6) =  xphi
    xhf(7) = -xphi
    !
    H_Hf=HF_hamiltonian(xhf)
    H_Hf=H_Hf+Hk_w90    
    !
    call fix_mu(H_Hf,delta_hf,mu_fix)
    !
    do ir=ir0,ir0+1
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
    !
    xhf_=xhf_-xhf
    do i=1,5
       out_x(i)=dreal(xhf_(i))
       out_x(i+5)=dimag(xhf_(i))
    end do
    !
  end function root_find_inner_loop


  function root_find_inner_loop_free(x) result(out_x)
    implicit none
    real(8),dimension(:) :: x
    real(8),dimension(size(x)) :: out_x
    complex(8),dimension(7) :: xhf,xhf_
    complex(8),dimension(Nso,Nso,Lk) :: H_Hf
    complex(8),dimension(Nso,Nso,Lk) :: delta_hf
    complex(8),dimension(Nso,Nso,nrpts) :: delta_hfr
    integer :: i,j
    real(8) :: mu_fix

    if(size(x).ne.10) stop 'delta HF s(x)/= 10'
    !
    do i=1,5
       xhf(i) = x(i) + xi*x(i+5)
    end do
    xhf(6) =  xpi(1)
    xhf(7) =  xpi(2)
    !
    H_Hf=HF_hamiltonian(xhf)
    H_Hf=H_Hf+Hk_w90    
    !
    call fix_mu(H_Hf,delta_hf,mu_fix)
    !
    do ir=ir0,ir0+1
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
    !
    xhf_=xhf_-xhf
    do i=1,5
       out_x(i)=dreal(xhf_(i))
       out_x(i+5)=dimag(xhf_(i))
    end do
    !
  end function root_find_inner_loop_free




  function root_find_inner_loop_symm(x) result(out_x)
    implicit none
    real(8),dimension(:) :: x
    real(8),dimension(size(x)) :: out_x
    complex(8),dimension(7) :: xhf,xhf_
    complex(8),dimension(Nso,Nso,Lk) :: H_Hf
    complex(8),dimension(Nso,Nso,Lk) :: delta_hf
    complex(8),dimension(Nso,Nso,nrpts) :: delta_hfr
    integer :: i,j
    real(8) :: mu_fix

    if(size(x).ne.10) stop 'delta HF s(x)/= 10'
    !
    do i=1,5
       xhf(i) = x(i) + xi*x(i+5)
    end do
    xhf(6) =  xphi
    xhf(7) =  xphi
    !
    H_Hf=HF_hamiltonian(xhf)
    H_Hf=H_Hf+Hk_w90    
    !
    call fix_mu(H_Hf,delta_hf,mu_fix)
    !
    do ir=ir0,ir0+1
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
    !
    xhf_=xhf_-xhf
    do i=1,5
       out_x(i)=dreal(xhf_(i))
       out_x(i+5)=dimag(xhf_(i))
    end do
    !
  end function root_find_inner_loop_symm




  function root_find_HF(x) result(out_x)
    implicit none
    real(8),dimension(:) :: x
    real(8),dimension(size(x)) :: out_x
    complex(8),dimension(7) :: xhf,xhf_
    complex(8),dimension(Nso,Nso,Lk) :: H_Hf
    complex(8),dimension(Nso,Nso,Lk) :: delta_hf
    complex(8),dimension(Nso,Nso,nrpts) :: delta_hfr
    integer :: i,j
    real(8) :: mu_fix

    if(size(x).ne.14) stop 'delta HF s(x)/= 14'

    do i=1,7
       xhf(i) = x(i) + xi*x(i+7)
    end do

    !
    H_Hf=HF_hamiltonian(xhf)
    H_Hf=H_Hf+Hk_w90
    !
    call fix_mu(H_Hf,delta_hf,mu_fix)
    !
    do ir=1,nrpts
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
    xhf_(6) = delta_hfr(1,3,ir0)+delta_hfr(1,3,ir0+1)
    xhf_(7) = delta_hfr(2,3,ir0)+delta_hfr(2,3,ir0+1)
    !
    ! xhf_(1) = delta_hfr(1,1,ir0)
    ! xhf_(2) = delta_hfr(2,2,ir0)
    ! xhf_(3) = delta_hfr(1,2,ir0)
    ! xhf_(4) = delta_hfr(1,2,ir0)+delta_hfr(1,2,ir0+1)
    !
    xhf_=xhf_-xhf
    do i=1,7
       out_x(i)=dreal(xhf_(i))
       out_x(i+7)=dimag(xhf_(i))
    end do

    !
  end function root_find_HF






  function HF_hamiltonian_symm(x_iter) result(Hhf)
    implicit none
    complex(8) :: x_iter(7)
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
       end do
    end do
    !
  end function HF_hamiltonian_symm
  !
  function HF_hamiltonian(x_iter) result(Hhf)
    implicit none
    complex(8) :: x_iter(7)
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
       end do
    end do
    !
  end function HF_hamiltonian
    





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


  subroutine build_mp_grid(Nx,Ny,Nz)
    implicit none
    integer :: Nx,Ny,Nz
    integer :: Nx_,Ny_,Nz_
    integer :: i,ix,iy,ik,iz,unitk
    real(8) :: kmax,kxmax,kymax,kzmax,dx,dy,dz
    real(8),dimension(:),allocatable :: kxr,kyr,kzr
    !
    !    
    if(allocated(kxgrid)) deallocate(kxgrid)
    if(allocated(kygrid)) deallocate(kygrid)
    if(allocated(kzgrid)) deallocate(kzgrid)
    allocate(kxgrid(Nx),kygrid(Ny),kzgrid(Nz))
    kmax=dble(Nx-1)*0.5d0/dble(Nx)
    kxgrid = linspace(-kmax,kmax,Nx,mesh=dx)
    kxmax=kmax+dx*0.5d0
    !
    kmax=dble(Ny-1)*0.5d0/dble(Ny)
    kygrid = linspace(-kmax,kmax,Ny,mesh=dy)
    kymax=kmax+dy*0.5d0
    !
    kmax=dble(Nz-1)*0.5d0/dble(Nz)
    kzgrid = linspace(-kmax,kmax,Nz,mesh=dz)
    kzmax=kmax+dz*0.5d0
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
    ! allocate(kxgrid_aux(3*Nx),kygrid_aux(3*Ny),kzgrid_aux(3*Nz))
    ! allocate(ix_aux(3*Nx),iy_aux(3*Ny),iz_aux(3*Nz))
    ! !
    ! do ix=1,Nx
    !    !
    !    kxgrid_aux(ix) = kxgrid(ix) - 2*kxmax
    !    kxgrid_aux(ix+Nx) = kxgrid(ix) 
    !    kxgrid_aux(ix+2*Nx) = kxgrid(ix) + 2*kxmax
    !    !
    !    ix_aux(ix) = ix
    !    ix_aux(ix+Nx) = ix
    !    ix_aux(ix+2*Nx) = ix
    !    !
    ! end do
    ! do iy=1,Ny
    !    !
    !    kygrid_aux(iy) = kygrid(iy) - 2*kymax
    !    kygrid_aux(iy+Ny) = kygrid(iy) 
    !    kygrid_aux(iy+2*Ny) = kygrid(iy) + 2*kymax
    !    !
    !    iy_aux(iy) = iy
    !    iy_aux(iy+Ny) = iy
    !    iy_aux(iy+2*Ny) = iy
    !    !
    ! end do
    ! do iz=1,Nz
    !    !
    !    kzgrid_aux(iz) = kzgrid(iz) - 2*kzmax
    !    kzgrid_aux(iz+Nz) = kzgrid(iz) 
    !    kzgrid_aux(iz+2*Nz) = kzgrid(iz) + 2*kzmax
    !    !
    !    iz_aux(iz) = iz
    !    iz_aux(iz+Nz) = iz
    !    iz_aux(iz+2*Nz) = iz
    !    !
    ! end do

    unitk=free_unit()
    open(unit=unitk,file='mp_grid.out')
    do ix=1,Nx
       do iy=1,Ny
          do iz=1,Nz
             write(unitk,'(10F18.10)') kxgrid(ix),kygrid(iy),kzgrid(iz)
          end do
       end do
    end do
    close(unitk)
    ! unitk=free_unit()
    ! open(unit=unitk,file='mp_grid_aux.out')    
    ! do ix=1,3*Nx
    !    do iy=1,3*Ny
    !       do iz=1,3*Nz
    !          write(unitk,'(10F18.10)') kxgrid_aux(ix),kygrid_aux(iy),kzgrid_aux(iz)
    !       end do
    !    end do
    ! end do
    ! close(unitk)
    !


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
    ! Ny_=Ny*N_cut_off
    ! Nz_=Nz*N_cut_off

    ! allocate(kxr(Nx),kyr(Ny),kzr(Nz))
    ! kmax=dble(Nx_-1)*0.5d0/dble(Nx_)
    ! kxr = linspace(-kmax,kmax,Nx_)
    ! !
    ! kmax=dble(Ny_-1)*0.5d0/dble(Ny_)
    ! kyr = linspace(-kmax,kmax,Ny_)
    ! !
    ! kmax=dble(Nz_-1)*0.5d0/dble(Nz_)
    ! kzr = linspace(-kmax,kmax,Nz_)
    ! !
    ! if(mod(Nx_,2).eq.0) then
    !    kxr=kxr-0.5d0/dble(Nx_)
    ! end if
    ! if(mod(Nx_,2).eq.0) then
    !    kyr=kyr-0.5d0/dble(Ny_)
    ! end if
    ! if(mod(Nz_,2).eq.0) then
    !    kzr=kzr-0.5d0/dble(Nz_)
    ! end if
    ! !
    ! Lkr = Nx_*Ny_*Nz_
    ! allocate(krl(Lkr,3),wtk_rl(Lkr),ikrl2ii(Lkr,3))
    ! wtk_rl=1.d0/dble(Lkr)
    ! !
    ! !
    ! unitk=free_unit()
    ! open(unit=unitk,file='kRL_points.out')
    ! ik=0
    ! krl=0
    ! do ix=1,Nx_
    !    do iy=1,Nx_
    !       do iz=1,Nz_
    !          ik = ik + 1          
    !          !
    !          krl(ik,:)=krl(ik,:)+kxr(ix)*Bk1(:)
    !          krl(ik,:)=krl(ik,:)+kyr(iy)*Bk2(:)
    !          krl(ik,:)=krl(ik,:)+kzr(iz)*Bk3(:)
    !          !
    !          ikrl2ii(ik,1)=ix
    !          ikrl2ii(ik,2)=iy
    !          ikrl2ii(ik,3)=iz
    !          write(unitk,'(10F18.10)') krl(ik,:)
    !       end do
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
  function Uq_ij(iq)
    USE DMFT_VECTORS  
    USE VARS_GLOBAL
    implicit none
    !
    integer :: iq,jq
    complex(8),dimension(Nso,Nso) :: Uq_ij
    real(8) :: mod
    integer :: idelta(3),iqvec(3),jqvec(3),ix,iy,iz
    real(8) ::kx_tmp,ky_tmp,kz_tmp
    !
    Uq_ij=Uq_TNS(:,:,iq)
    if(only_loc) Uq_ij=check_Uloc
    !
  end function Uq_ij




  function Uq_jellium(q)
    USE DMFT_VECTORS  
    USE VARS_GLOBAL
    implicit none
    !
    real(8),dimension(:)       :: q
    real(8),dimension(Nso,Nso) :: Uq_jellium
    real(8) :: mod
    integer :: i
    !
    mod=0.d0
    do i=1,size(q)
       mod=mod+q(i)**2.0
    end do
    mod=sqrt(mod)
    Uq_jellium=Uloc_TNS
    if(.not.hartree) Uq_jellium = 0.d0
    if(mod.ne.0.d0) then
       Uq_jellium=Uloc_TNS
       !Uq_jellium = Umat/mod
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


  !+- get the truncated interaction -+!
  ! open(unit_in,file='Ur_truncated.tns')
  ! Ur_TNS=0.d0
  ! do i=1,2*Ncut+1
  !    do j=1,2*Ncut+1
  !       do k=1,2*Ncut+1
  !          Rlat=(i-(Ncut+1))*R1+(j-(Ncut+1))*R2+(k-(Ncut+1))*R3           
  !          modk=sqrt(Rlat(1)**2.d0+Rlat(2)**2.d0+Rlat(3)**2.d0)
  !          do ir=1,nrpts
  !             rpt_latt(ir,:)=irvec(ir,1)*R1+irvec(ir,2)*R2+irvec(ir,3)*R3              
  !             CheckR=sqrt(rpt_latt(ir,1)**2.d0+rpt_latt(ir,2)**2.d0+rpt_latt(ir,3)**2.d0)
  !             if(abs(modk-CheckR).lt.1.d-10) then
  !                write(*,*) i-(Ncut+1),j-(Ncut+1),k-(Ncut+1)
  !                call FT_q2r(rpt_latt(ir,:),Ur_TNS(:,:,ir),Uq_TNS)
  !                exit
  !             end if
  !          end do
  !       end do
  !    end do
  ! end do
  ! do ir=1,nrpts
  !    rpt_latt(ir,:)=irvec(ir,1)*R1+irvec(ir,2)*R2+irvec(ir,3)*R3     
  !    modk=rpt_latt(ir,1)**2.d0+rpt_latt(ir,2)**2.d0+rpt_latt(ir,3)**2.d0
  !    !     
  !    modk=sqrt(modk)
  !    write(unit_in,'(50F10.5)') modk,rpt_latt(ir,1:3),Ur_TNS(1,1:6,ir),Ur_TNS(2,2:6,ir),Ur_TNS(3,3:6,ir),Ur_TNS(4,4:6,ir),Ur_TNS(5,5:6,ir),Ur_TNS(6,6,ir)
  ! end do
  ! close(unit_in)
  ! stop
