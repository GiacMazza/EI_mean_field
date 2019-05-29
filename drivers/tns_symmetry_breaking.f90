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
  integer :: Nhf,ihf,unit_err,unit_obs,unit_in,uio,Nobs
  integer,dimension(:),allocatable :: units_loc_obs
  integer :: flen,iread

  integer :: ir,jr,kr,iir
  real(8),dimension(3) :: Rlat 

  complex(8),dimension(:),allocatable :: obs_loc
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
  !
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
  real(8) :: alphaU,deltar
  real(8),dimension(6) :: tk


  integer :: Nkpath
  real(8) :: vhyb
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
  call parse_input_variable(Nk_y,"Nk_y","input.conf",default=12)
  call parse_input_variable(Nk_z,"Nk_z","input.conf",default=3)
  call parse_input_variable(Nhf,"Nhf","input.conf",default=100)
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

  call parse_input_variable(Nkpath,"Nkpath","input.conf",default=300)
  call parse_input_variable(vhyb,"vhyb","input.conf",default=0.05d0)


  !
  call get_global_vars
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
  allocate(kpath(4,3))
  ! !+- M-point
  ! kpath(1,:)=0.5d0*Bk1+0.5d0*Bk3
  ! !+- Z-point
  ! kpath(2,:)=0.5d0*Bk3
  ! !+- G-point
  ! kpath(3,:)=0.d0
  ! !+- X-point
  ! kpath(4,:)=0.5d0*Bk1

  !+- X-point
  kpath(1,:)=0.5d0*Bk1;write(*,*) kpath(1,:)
  !+- G-point
  kpath(2,:)=0.d0;write(*,*) kpath(2,:)
  !+- X-point
  kpath(3,:)=0.5d0*Bk1;write(*,*) kpath(3,:)
  !+- M-point
  kpath(4,:)=0.5d0*Bk1+0.5d0*Bk3;write(*,*) kpath(4,:)

  

  !+- build a monkhorst-pack grid -+!  
  call build_mp_grid(Nk_x,Nk_y,Nk_z)
  !
  Lk=Nk_x*Nk_y*Nk_z
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
              kpt_latt(ik,idim) = kpt_latt(ik,idim) + kygrid(j)*Bk2(idim)
              kpt_latt(ik,idim) = kpt_latt(ik,idim) + kzgrid(k)*Bk3(idim)
           end do
           write(300,*) kpt_latt(ik,:)
        end do
     end do
  end do
  !
  ! Lk_aux=3*3*3*Lk
  ! allocate(kpt_latt_aux(Lk_aux,3),ik_stride_aux(Lk_aux,3),igr2ik_aux(3*Nk_x,3*Nk_y,3*Nk_z))
  ! ik=0
  ! kpt_latt_aux=0.d0
  ! do i=1,3*Nk_x
  !    do k=1,3*Nk_z
  !       do j=1,3*Nk_y
  !          !
  !          ik=ik+1
  !          !
  !          ik_stride_aux(ik,1)=i
  !          ik_stride_aux(ik,2)=j
  !          ik_stride_aux(ik,3)=k
  !          igr2ik_aux(i,j,k) = ik
  !          !
  !          do idim=1,3
  !             kpt_latt_aux(ik,idim) = kpt_latt_aux(ik,idim) + kxgrid_aux(i)*Bk1(idim)
  !             kpt_latt_aux(ik,idim) = kpt_latt_aux(ik,idim) + kygrid_aux(j)*Bk2(idim)
  !             kpt_latt_aux(ik,idim) = kpt_latt_aux(ik,idim) + kzgrid_aux(k)*Bk3(idim)
  !          end do
  !          !
  !       end do
  !    end do
  ! end do
  !  
  uio=free_unit()
  open(unit=uio,file='mp_kgrid.f90')
  do ik=1,Lk
     write(uio,'(10F18.10)') kpt_latt(ik,:)
  end do
  close(uio)

  ! uio=free_unit()
  ! open(unit=uio,file='mp_kgrid_aux.f90')
  ! do ik=1,Lk_aux
  !    write(uio,'(10F18.10)') kpt_latt_aux(ik,:)
  ! end do
  ! close(uio)


  !+---------------------------------+!  
  file_name=reg(read_tns)//reg(file_w90_hr)  
  !+- read the w90 output -+!
  allocate(Hk_w90(Nso,Nso,Lk),Hloc(Nso,Nso))
  call read_w90_hr(R1,R2,R3,Hr_w90,Hloc,irvec,ndegen,trim(file_name),1,6,1)
  nrpts=size(irvec,1)
  allocate(rpt_latt(nrpts,3))
  do ir=1,nrpts
     rpt_latt(ir,:)=irvec(ir,1)*R1+irvec(ir,2)*R2+irvec(ir,3)*R3     
  end do
  !
  allocate(itmp(nrpts,2))  
  allocate(stride2D_(nrpts))
  i=0
  do ir=1,nrpts
     if(irvec(ir,2).eq.0) then
        i=i+1
        itmp(i,1) = irvec(ir,1)
        itmp(i,2) = irvec(ir,3)
        stride2D_(i) = ir
     end if
  end do
  nr2d=i
  write(*,*) "nr-points on xz plane",nr2d
  allocate(irvec2d(nr2d,2));   allocate(stride2D(nr2d)); 
  irvec2d(:,1)=itmp(1:nr2d,1)
  irvec2d(:,2)=itmp(1:nr2d,2)
  stride2D=stride2D_(1:nr2d)
  !
  !+- read local interaction term -+!
  allocate(Uloc_TNS(Nso,Nso))
  flen=file_length('TNS_Uloc.dat')
  if(flen.ne.Nso*Nso) stop "error in reading Uloc"
  unit_in=free_unit(); 
  open(unit=unit_in,file='TNS_Uloc.dat',status="old",action="read")  
  do i=1,Nso
     do j=1,Nso
        read(unit_in,*) tmpi,tmpj,Uloc_TNS(j,i)
     end do
  end do
  close(unit_in)
  
  !+- read U(q) -+!
  file_name=reg(read_tns)//reg(file_UV)
  open(unit=unit_in,file=trim(file_name),status="old",action="read")  
  read(unit_in,*) 
  read(unit_in,*) 
  allocate(read_tmp(9))
  !
  allocate(Uq_TNS(Nso,Nso,Lk),Vq_TNS(Nso,Nso,Lk),check_Uloc(Nso,Nso))
  check_Uloc=0.d0
  ik=0
  do ix=1,Nk_x
     do iz=1,Nk_z
        do iy=1,Nk_y
           ik=ik+1
           !
           do i=1,Nso
              do j=1,Nso
                 read(unit_in,*) read_tmp(1:9),Vq_read,Uq_read
                 Uq_TNS(i,j,ik) = Uq_read(1)+xi*Uq_read(2)
                 Vq_TNS(i,j,ik) = Vq_read(1)+xi*Vq_read(2)                 
                 check_Uloc(i,j) = check_Uloc(i,j) + Uq_TNS(i,j,ik)*wtk(ik)
              end do
           end do
           !
        end do
     end do
  end do
  close(unit_in)

  open(unit_in,file='read_Uq.tns')
  do ik=1,Lk
     !+- here print out Uq vs |q| -+!
     modk=kpt_latt(ik,1)**2.d0+kpt_latt(ik,2)**2.d0+kpt_latt(ik,3)**2.d0
     modk=sqrt(modk)     
     write(unit_in,'(50F10.5)') modk,Uq_TNS(1,1:6,ik),Uq_TNS(2,2:6,ik),Uq_TNS(3,3:6,ik),Uq_TNS(4,4:6,ik),Uq_TNS(5,5:6,ik),Uq_TNS(6,6,ik)
  end do
  close(unit_in)
  allocate(Ur_TNS(Nso,Nso,nrpts))
  open(unit_in,file='Ur_full_space.tns')
  !+- get the full real space interaction -+!
  Ur_TNS=0.d0
  do ir=1,nrpts
     !
     call FT_q2r(rpt_latt(ir,:),Ur_TNS(:,:,ir),Uq_TNS)
     modk=rpt_latt(ir,1)**2.d0+rpt_latt(ir,2)**2.d0+rpt_latt(ir,3)**2.d0
     modk=sqrt(modk)
     write(unit_in,'(50F10.5)') modk,rpt_latt(ir,1:3),Ur_TNS(1,1:6,ir),Ur_TNS(2,2:6,ir),Ur_TNS(3,3:6,ir),Ur_TNS(4,4:6,ir),Ur_TNS(5,5:6,ir),Ur_TNS(6,6,ir)
     !
  end do
  close(unit_in)
  Rlat=-1000.0
  do ir=1,nrpts
     do i=1,3
        if(rpt_latt(ir,i).gt.Rlat(i)) Rlat(i)=rpt_latt(ir,i)
     end do
  end do
  write(*,*) 'largest components of R vectors in the supercell',Rlat
  checkR=1000.d0
  do i=1,3
     if(Rlat(i).lt.checkR) checkR=Rlat(i)
  end do
  write(*,*) 'largest radius sphere in the supercell',checkR,Rlat
  if(checkR.lt.Ucut_off) then
     write(*,*) 'cut-off larger than the largest radius sphere in the supercell' 
     Ucut_off=checkR
     write(*,*) 'New cut-off set to',Ucut_off 
  end if
  open(unit_in,file='Ur_truncated.tns')
  Ur_TNS=0.d0
  do ir=1,nrpts
     modk=rpt_latt(ir,1)**2.d0+rpt_latt(ir,2)**2.d0+rpt_latt(ir,3)**2.d0
     !     
     modk=sqrt(modk)
     if(modk.lt.Ucut_off) then
        call FT_q2r(rpt_latt(ir,:),Ur_TNS(:,:,ir),Uq_TNS)
     end if
     write(unit_in,'(50F10.5)') modk,rpt_latt(ir,1:3),Ur_TNS(1,1:6,ir),Ur_TNS(2,2:6,ir),Ur_TNS(3,3:6,ir),Ur_TNS(4,4:6,ir),Ur_TNS(5,5:6,ir),Ur_TNS(6,6,ir)
  end do
  close(unit_in)  
  !
  !+- transform back to momentum space -+!
  !
  open(unit_in,file='Uq_truncated.tns')
  Uq_TNS=0.d0
  do ik=1,Lk
     call FT_r2q(kpt_latt(ik,:),Uq_TNS(:,:,ik),Ur_TNS)
     modk=kpt_latt(ik,1)**2.d0+kpt_latt(ik,2)**2.d0+kpt_latt(ik,3)**2.d0
     modk=sqrt(modk)
     write(unit_in,'(50F9.5)') modk,Uq_TNS(1,1:6,ik),Uq_TNS(2,2:6,ik),Uq_TNS(3,3:6,ik),Uq_TNS(4,4:6,ik),Uq_TNS(5,5:6,ik),Uq_TNS(6,6,ik)
  end do
  close(unit_in)
  ! 

  !+- new MP grid -+!
  Nk_x=1000
  Nk_x=2
  Nk_x=2
  call build_mp_grid(Nk_x,Nk_y,Nk_z)
  !
  deallocate(kpt_latt,ik_stride,wtk,igr2ik)
  Lk=Nk_x*Nk_y*Nk_z
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
              kpt_latt(ik,idim) = kpt_latt(ik,idim) + kygrid(j)*Bk2(idim)
              kpt_latt(ik,idim) = kpt_latt(ik,idim) + kzgrid(k)*Bk3(idim)
           end do
           write(301,*) kpt_latt(ik,:)
        end do
     end do
  end do
  
  
  !
  allocate(Hr_w90_tmp(Nso,Nso,nrpts)); Hr_w90_tmp=Hr_w90
  deallocate(Hr_w90)
  allocate(Ur_tmp(Nso,Nso,nrpts)) ; Ur_tmp=Ur_TNS
  deallocate(Ur_TNS)            !
  !
  Nso=Norb*Nspin
  allocate(Hr_w90(Nso,Nso,nrpts))
  !
  do ispin=1,Nspin
     do iorb=1,Norb
        iso=(ispin-1)*Norb+iorb
        do jorb=1,Norb
           jso=(ispin-1)*Norb+jorb
           Hr_w90(iso,jso,:) = Hr_w90_tmp(iorb,jorb,:)
        end do
     end do
  end do
  deallocate(Hr_w90_tmp)
  !

  Hr_w90=dreal(Hr_w90)

  !+- bare bands -+!
  allocate(delta_hfr(Nso,Nso,nrpts),delta_hfr_(Nso,Nso,nrpts),H_hf(Nso,Nso,nrpts))
  call init_var_params_latt(delta_hfr,Hr_w90)
  
  ! do ispin=1,Nspin
  !    do iorb=1,4
  !       iso=(ispin-1)*Norb+iorb
  !       Hr_w90(iso,iso,ir0) = Hr_w90(iso,iso,ir0) + cf_ext
  !    end do
  ! end do
  ! do ispin=1,Nspin
  !    do iorb=5,6
  !       iso=(ispin-1)*Norb+iorb
  !       Hr_w90(iso,iso,ir0) = Hr_w90(iso,iso,ir0) - cf_ext
  !    end do
  ! end do


  if(H1d) then
     !
     tk(1:4) = -0.8d0
     tk(5:6) = 0.4d0
     !
     Hr_w90=0.d0
     do ir=1,nrpts
        Rlat=irvec(ir,1)*R1+irvec(ir,2)*R2+irvec(ir,3)*R3
        Rlat=Rlat+R1
        deltar=sqrt(Rlat(1)**2.d0+Rlat(2)**2.d0+Rlat(3)**2.d0)
        if(deltar.lt.1.d-10) then
           do ispin=1,Nspin
              do iorb=1,Norb
                 iso=(ispin-1)*Norb+iorb        
                 Hr_w90(iso,iso,ir) = tk(iorb)
              end do
              !
              jso=(ispin-1)*Norb+1
              do iorb=5,5
                 iso=(ispin-1)*Norb+iorb         
                 Hr_w90(iso,jso,ir) = vhyb
              end do
              jso=(ispin-1)*Norb+2
              do iorb=5,5
                 iso=(ispin-1)*Norb+iorb         
                 Hr_w90(iso,jso,ir) = vhyb
              end do
              !
              jso=(ispin-1)*Norb+3
              do iorb=6,6
                 iso=(ispin-1)*Norb+iorb         
                 Hr_w90(iso,jso,ir) = vhyb
              end do
              jso=(ispin-1)*Norb+4
              do iorb=6,6
                 iso=(ispin-1)*Norb+iorb         
                 Hr_w90(iso,jso,ir) = vhyb
              end do
              !
           end do
        end if
        Rlat=irvec(ir,1)*R1+irvec(ir,2)*R2+irvec(ir,3)*R3
        Rlat=Rlat-R1
        deltar=sqrt(Rlat(1)**2.d0+Rlat(2)**2.d0+Rlat(3)**2.d0)
        if(deltar.lt.1.d-10) then
           do ispin=1,Nspin
              do iorb=1,Norb
                 iso=(ispin-1)*Norb+iorb        
                 Hr_w90(iso,iso,ir) = tk(iorb)
              end do
              !
              jso=(ispin-1)*Norb+5
              do iorb=1,2
                 iso=(ispin-1)*Norb+iorb         
                 Hr_w90(iso,jso,ir) = vhyb
              end do
              !
              jso=(ispin-1)*Norb+6
              do iorb=3,4
                 iso=(ispin-1)*Norb+iorb         
                 Hr_w90(iso,jso,ir) = vhyb
              end do

           end do
        end if
     end do
     !
     do ispin=1,Nspin
        do iorb=1,4
           iso=(ispin-1)*Norb+iorb
           Hr_w90(iso,iso,ir0) = Hr_w90(iso,iso,ir0) + 1.7
        end do
        !
        jso=(ispin-1)*Norb+1
        do iorb=5,5
           iso=(ispin-1)*Norb+iorb         
           Hr_w90(iso,jso,ir0) = -vhyb
           Hr_w90(jso,iso,ir0) = -vhyb
        end do
        jso=(ispin-1)*Norb+2
        do iorb=5,5
           iso=(ispin-1)*Norb+iorb         
           Hr_w90(iso,jso,ir0) = -vhyb
           Hr_w90(jso,iso,ir0) = -vhyb
        end do
        !
        jso=(ispin-1)*Norb+3
        do iorb=6,6
           iso=(ispin-1)*Norb+iorb         
           Hr_w90(iso,jso,ir0) = -vhyb
           Hr_w90(jso,iso,ir0) = -vhyb
        end do
        jso=(ispin-1)*Norb+4
        do iorb=6,6
           iso=(ispin-1)*Norb+iorb         
           Hr_w90(iso,jso,ir0) = -vhyb
           Hr_w90(jso,iso,ir0) = -vhyb
        end do
        !
     end do
     do ispin=1,Nspin
        do iorb=5,6
           iso=(ispin-1)*Norb+iorb
           Hr_w90(iso,iso,ir0) = Hr_w90(iso,iso,ir0) - 0.9
        end do
     end do
  end if
  !
  mu_fix=4.d0
  call find_chem_pot_latt(Hr_w90,delta_hfr,mu_fix)
  
  allocate(kpt_path(3*Nkpath,3))
  allocate(Hktmp(Nso,Nso))
  allocate(ek_out(Nso))

  uio=free_unit()
  open(unit_in,file='H_Gamma_point.tns')
  Hktmp=0.d0
  do ir=1,nrpts
     Hktmp = Hktmp + Hr_w90(:,:,ir)
  end do
  do iso=1,Nso
     do jso=1,Nso
        write(unit_in,*) iso,jso,Hktmp(iso,jso)
     end do
  end do
  close(unit_in)

  uio=free_unit()
  open(unit=uio,file='tns_bare_bands.out')
  unit_in=free_unit()
  open(unit_in,file='tns_fat_bare.tns')
  !
  allocate(Ta_fat(Nso),Ni_fat(Nso))
  modk=0.d0
  delta_kpath=kpath(2,:)-kpath(1,:)
  modk=modk-sqrt(dot_product(delta_kpath/dble(Nkpath),delta_kpath/dble(Nkpath)))
  do i=1,3
     delta_kpath=kpath(i+1,:)-kpath(i,:)
     do ik=1,Nkpath
        j=(i-1)*Nkpath + ik
        kpt_path(j,:) = kpath(i,:) + dble(ik-1)/dble(Nkpath)*delta_kpath
        modk=modk+sqrt(dot_product(delta_kpath/dble(Nkpath),delta_kpath/dble(Nkpath)))
        ! modk=modk+sqrt(dot_product(1.d0/100.d0*delta_kpath,1.d0/100.d0*delta_kpath))
        !
        Hktmp=0.d0
        call FT_r2q(kpt_path(j,:),Hktmp,Hr_w90)
        !
        !
        write(264,'(30F18.10)') modk,Hktmp(1,1:6)
        write(265,'(30F18.10)') modk,Hktmp(7,7:12)
        !
        call eigh(Hktmp,ek_out)
        ! Ta_fat=abs(Hktmp(1,:))**2.d0+abs(Hktmp(2,:))**2.d0+abs(Hktmp(3,:))**2.d0+abs(Hktmp(4,:))**2.d0
        ! Ni_fat=abs(Hktmp(5,:))**2.d0+abs(Hktmp(6,:))**2.d0

        Ta_fat=abs(Hktmp(1,:))**2.d0+abs(Hktmp(2,:))**2.d0+abs(Hktmp(3,:))**2.d0+abs(Hktmp(4,:))**2.d0
        Ni_fat=abs(Hktmp(5,:))**2.d0+abs(Hktmp(6,:))**2.d0
        !
        do iso=1,Nso
           if(Ta_fat(iso).lt.1.d-5) Ta_fat(iso)=0.d0;Ni_fat(iso)=1.d0-Ta_fat(iso)
           if(Ni_fat(iso).lt.1.d-5) Ni_fat(iso)=0.d0;Ta_fat(iso)=1.d0-Ni_fat(iso)
        end do
        
        write(uio,'(30F18.10)') modk,ek_out-mu_fix,Ta_fat,Ni_fat
        write(unit_in,'(30F18.10)') modk,Ta_fat,Ni_fat
        !
     end do
     write(*,*) modk
  end do
  close(uio)
  close(unit_in)  

  

  !
  open(unit_in,file='Deltar_bare_real.tns')
  modk=0.d0
  do ir=1,nrpts
     modk=sqrt(rpt_latt(ir,1)**2.d0+rpt_latt(ir,2)**2.d0+rpt_latt(ir,3)**2.d0)
     !     
     write(unit_in,'(50F10.5)') modk,rpt_latt(ir,1:3), &
          dreal(delta_hfr(1,1,ir)), &
          dreal(delta_hfr(5,5,ir))

          ! dreal(delta_hfr(2,2:6,ir)), &
          ! dreal(delta_hfr(3,3:6,ir)), &
          ! dreal(delta_hfr(4,4:6,ir)), &
          ! dreal(delta_hfr(5,5:6,ir)), &
          ! dreal(delta_hfr(6,6,ir))
  end do
  close(unit_in)
  open(unit_in,file='Deltar_bare_imag.tns')
  modk=0.d0
  do ir=1,nrpts
     modk=sqrt(rpt_latt(ir,1)**2.d0+rpt_latt(ir,2)**2.d0+rpt_latt(ir,3)**2.d0)
     !     
     write(unit_in,'(50F10.5)') modk,rpt_latt(ir,1:3), &
          dimag(delta_hfr(1,1:6,ir)), &
          dimag(delta_hfr(2,2:6,ir)), &
          dimag(delta_hfr(3,3:6,ir)), &
          dimag(delta_hfr(4,4:6,ir)), &
          dimag(delta_hfr(5,5:6,ir)), &
          dimag(delta_hfr(6,6,ir))
  end do
  close(unit_in)
  
  open(unit_in,file='2d_deltar_bare_real.tns')
  modk=0.d0
  do ir=1,nr2d
     Rlat=irvec2d(ir,1)*R1+irvec2d(ir,2)*R3
     modk=sqrt(Rlat(1)**2.d0+Rlat(2)**2.d0+Rlat(3)**2.d0)
     !     
     write(unit_in,'(50F10.5)') modk,Rlat(1:3), &
          dreal(delta_hfr(1,1:6,stride2D(ir))), &
          dreal(delta_hfr(2,2:6,stride2D(ir))), &
          dreal(delta_hfr(3,3:6,stride2D(ir))), &
          dreal(delta_hfr(4,4:6,stride2D(ir))), &
          dreal(delta_hfr(5,5:6,stride2D(ir))), &
          dreal(delta_hfr(6,6,stride2D(ir)))
  end do
  close(unit_in)
  open(unit_in,file='2d_deltar_bare_imag.tns')
  modk=0.d0
  do ir=1,nr2d
     Rlat=irvec2d(ir,1)*R1+irvec2d(ir,2)*R3
     modk=sqrt(Rlat(1)**2.d0+Rlat(2)**2.d0+Rlat(3)**2.d0)
     !     
     write(unit_in,'(50F10.5)') modk,Rlat(1:3), &
          dimag(delta_hfr(1,1:6,stride2D(ir))), &
          dimag(delta_hfr(2,2:6,stride2D(ir))), &
          dimag(delta_hfr(3,3:6,stride2D(ir))), &
          dimag(delta_hfr(4,4:6,stride2D(ir))), &
          dimag(delta_hfr(5,5:6,stride2D(ir))), &
          dimag(delta_hfr(6,6,stride2D(ir)))
  end do
  close(unit_in)
  
  stop

  allocate(Ur_TNS(Nso,Nso,nrpts)); 
  do ispin=1,Nspin
     do iorb=1,Norb
        iso=(ispin-1)*Norb+iorb
        do jspin=1,Nspin
           do jorb=1,Norb
              jso=(jspin-1)*Norb+jorb              
              Ur_TNS(iso,jso,:) = Ur_tmp(iorb,jorb,:)
           end do
        end do
     end do
  end do
  if(only_loc) then
     Ur_TNS=0.d0
     write(*,*) 'LOCAL COULOMB INTERACTION'
     do ispin=1,Nspin
        do iorb=1,Norb
           iso=(ispin-1)*Norb+iorb
           do jspin=1,Nspin
              do jorb=1,Norb
                 jso=(jspin-1)*Norb+jorb              
                 Ur_TNS(iso,jso,ir0) = Uloc_TNS(iorb,jorb)
                 write(*,*) iso,jso,Ur_TNS(iso,jso,ir0)
              end do
           end do
        end do
     end do
  end if
  do ispin=1,Nspin
     do iorb=1,Norb
        iso=(ispin-1)*Norb+iorb
        Ur_TNS(iso,iso,ir0) = 0.d0
     end do
  end do
  !
  Ur_TNS=alphaU*Ur_TNS
  !  
  mu_fix=4.d0
  call find_chem_pot_latt(Hr_w90,delta_hfr,mu_fix)
  ncell=0.d0
  do iso=1,Nso
     write(*,*) delta_hfr(iso,iso,ir0)
     ncell=ncell+delta_hfr(iso,iso,ir0)
  end do
  write(*,*) mu_fix,ncell,ndens
  !
  Nobs=Nso*(Nso+1)/2
  allocate(units_loc_obs(Nobs))
  units_loc_obs=free_units(Nobs)
  iso=0
  do i=1,Nso
     do j=i,Nso
        iso=iso+1
        open(unit=units_loc_obs(iso),file="bare_obs_hf_"//reg(txtfy(i)//reg(txtfy(j))))        
        write(units_loc_obs(iso),'(10F18.10)') 0.d0,delta_hfr(i,j,ir0)
        close(units_loc_obs(iso))
     end do
  end do  
  iso=0
  do i=1,Nso
     do j=i,Nso
        iso=iso+1
        open(unit=units_loc_obs(iso),file="obs_hf_"//reg(txtfy(i)//reg(txtfy(j))))        
     end do
  end do
  
  do ispin=1,Nspin
     !do iorb=1,Norb
        !do jorb=iorb+1,Norb
     iorb=1
     jorb=5
     iso=(ispin-1)*Norb+iorb
     jso=(ispin-1)*Norb+jorb
     delta_hfr(iso,jso,ir0) = 0.1d0!delta_hfr(iso,jso,ir0) + 0.1d0
     delta_hfr(jso,iso,ir0) = 0.1d0!delta_hfr(jso,iso,ir0) + 0.1d0           

     iorb=2
     jorb=5
     iso=(ispin-1)*Norb+iorb
     jso=(ispin-1)*Norb+jorb
     delta_hfr(iso,jso,ir0) = 0.1d0!delta_hfr(iso,jso,ir0) + 0.1d0
     delta_hfr(jso,iso,ir0) = 0.1d0!delta_hfr(jso,iso,ir0) + 0.1d0           


     iorb=3
     jorb=6
     iso=(ispin-1)*Norb+iorb
     jso=(ispin-1)*Norb+jorb
     delta_hfr(iso,jso,ir0) = 0.1d0!delta_hfr(iso,jso,ir0) + 0.1d0
     delta_hfr(jso,iso,ir0) = 0.1d0!delta_hfr(jso,iso,ir0) + 0.1d0           


     iorb=4
     jorb=6
     iso=(ispin-1)*Norb+iorb
     jso=(ispin-1)*Norb+jorb
     delta_hfr(iso,jso,ir0) = 0.1d0!delta_hfr(iso,jso,ir0) + 0.1d0
     delta_hfr(jso,iso,ir0) = 0.1d0!delta_hfr(jso,iso,ir0) + 0.1d0           


        !end do
     !end do
  end do
  ! do iso=1,Nso
  !    do jso=iso+1,Nso
  !       delta_hfr(iso,jso,ir0) = delta_hfr(iso,jso,ir0) + 0.1d0
  !       delta_hfr(jso,iso,ir0) = delta_hfr(jso,iso,ir0) + 0.1d0
  !    end do
  ! end do
  !
  !delta_hfr=0.d0
  delta_hfr_=delta_hfr

  unit_obs=free_unit()
  open(unit=unit_obs,file='ndens_hf.loop')
 
  unit_err=free_unit()
  open(unit=unit_err,file='err_q.err')  
  err_hf=1.d0  
  do ihf=1,Nhf

     
     write(unit_err,'(10F18.10)') dble(ihf),err_hf,mu_fix               
     iso=0
     Ncell=0.d0
     do i=1,Nso           
        Ncell=Ncell+delta_hfr(i,i,ir0)
        do j=i,Nso
           iso=iso+1
           write(units_loc_obs(iso),'(10F18.10)') dble(ihf),delta_hfr(i,j,ir0)
        end do
     end do
     write(unit_obs,'(10F18.10)') dble(ihf),ncell

     H_hf=Hr_w90
     call build_HF_hamiltonian_latt(H_hf,delta_hfr,Ur_TNS)          
     call find_chem_pot_latt(H_hf,delta_hfr,mu_fix)

     ! if(Nspin.eq.2) then
     !    do iorb=1,Norb
     !       do jorb=1,Norb
     !          iso=iorb+Norb
     !          jso=jorb
     !          delta_hfr(iso,jso,:) = 0.d0
     !          iso=iorb
     !          jso=jorb+Norb
     !          delta_hfr(iso,jso,:) = 0.d0        
     !       end do
     !    end do
     ! end if
     
     delta_hfr=wmix*delta_hfr+(1.d0-wmix)*delta_hfr_


     !if(ihf.eq.2) stop
     !
     !
     err_hf=check_conv(delta_hfr,delta_hfr_)     
     !
     delta_hfr_=delta_hfr
     !
  end do
  close(unit_err)
  close(unit_obs)

  

  open(unit_in,file='Deltar_HF_real.tns')
  modk=0.d0
  do ir=1,nrpts
     Rlat=irvec(ir,1)*R1+irvec(ir,2)*R2+irvec(ir,3)*R3     
     modk=sqrt(Rlat(1)**2.d0+Rlat(2)**2.d0+Rlat(3)**2.d0)
     !     
     write(unit_in,'(50F10.5)') modk,Rlat(1:3), &
          dreal(delta_hfr(1,1,ir))
          ! dreal(delta_hfr(1,1:6,ir)), &
          ! dreal(delta_hfr(2,2:6,ir)), &
          ! dreal(delta_hfr(3,3:6,ir)), &
          ! dreal(delta_hfr(4,4:6,ir)), &
          ! dreal(delta_hfr(5,5:6,ir)), &
          ! dreal(delta_hfr(6,6,ir))
  end do
  close(unit_in)
  open(unit_in,file='Deltar_HF_imag.tns')
  modk=0.d0
  do ir=1,nrpts
     Rlat=irvec(ir,1)*R1+irvec(ir,2)*R2+irvec(ir,3)*R3     
     modk=sqrt(Rlat(1)**2.d0+Rlat(2)**2.d0+Rlat(3)**2.d0)
     !     
     write(unit_in,'(50F10.5)') modk,Rlat(1:3), &
          dimag(delta_hfr(1,1:6,ir)), &
          dimag(delta_hfr(2,2:6,ir)), &
          dimag(delta_hfr(3,3:6,ir)), &
          dimag(delta_hfr(4,4:6,ir)), &
          dimag(delta_hfr(5,5:6,ir)), &
          dimag(delta_hfr(6,6,ir))
  end do
  close(unit_in)

  open(unit_in,file='2d_deltar_HF_real.tns')
  modk=0.d0
  do ir=1,nr2d
     Rlat=irvec2d(ir,1)*R1+irvec2d(ir,2)*R3     
     modk=sqrt(Rlat(1)**2.d0+Rlat(2)**2.d0+Rlat(3)**2.d0)
     !     
     write(unit_in,'(50F10.5)') modk,Rlat(1:3), &
          dreal(delta_hfr(1,1:6,stride2D(ir))), &
          dreal(delta_hfr(2,2:6,stride2D(ir))), &
          dreal(delta_hfr(3,3:6,stride2D(ir))), &
          dreal(delta_hfr(4,4:6,stride2D(ir))), &
          dreal(delta_hfr(5,5:6,stride2D(ir))), &
          dreal(delta_hfr(6,6,stride2D(ir)))
  end do
  close(unit_in)
  open(unit_in,file='2d_deltar_HF_imag.tns')
  modk=0.d0
  do ir=1,nr2d
     Rlat=irvec2d(ir,1)*R1+irvec2d(ir,2)*R3     
     modk=sqrt(Rlat(1)**2.d0+Rlat(2)**2.d0+Rlat(3)**2.d0)
     !     
     write(unit_in,'(50F10.5)') modk,Rlat(1:3), &
          dimag(delta_hfr(1,1:6,stride2D(ir))), &
          dimag(delta_hfr(2,2:6,stride2D(ir))), &
          dimag(delta_hfr(3,3:6,stride2D(ir))), &
          dimag(delta_hfr(4,4:6,stride2D(ir))), &
          dimag(delta_hfr(5,5:6,stride2D(ir))), &
          dimag(delta_hfr(6,6,stride2D(ir)))
  end do
  close(unit_in)


  !+- spin-symmetriz delta_hfr -+!
  if(Nspin.eq.2) then
     do iorb=1,Norb
        do jorb=1,Norb
           iso=iorb+Norb
           jso=jorb
           ! delta_hfr(iso,jso,:) = 0.d0
           ! delta_hfr(jso,iso,:) = 0.d0        
        end do
     end do     
     do iorb=1,Norb
        do jorb=1,Norb
           do ir=1,nrpts
              iso=iorb+Norb
              jso=jorb+Norb
              !delta_hfr(iso,jso,ir) = delta_hfr(iorb,jorb,ir)
           end do
        end do
     end do
  end if

  !H_hf=Hr_w90
  !call build_HF_hamiltonian_latt(H_hf,delta_hfr,Ur_TNS)
  call find_chem_pot_latt(H_hf,delta_hfr,mu_fix)
  write(*,*) 'mu fix',mu_fix
  !do ir=1,nrpts
  do iso=1,Nso
     H_hf(iso,iso,ir0) =  H_hf(iso,iso,ir0) - mu_fix
  end do
  !end do
  !+- here plot bands -+! 
  !+- problems: bands are not spin symmetric -+!
  uio=free_unit()
  open(unit=uio,file='tns_bands.out')
  !
  Hktmp=0.d0
  modk=0.d0
  do i=1,3
     delta_kpath=kpath(i+1,:)-kpath(i,:)
     do ik=1,100
        j=(i-1)*100 + ik
        kpt_path(j,:) = kpath(i,:) + dble(ik-1)/100.d0*delta_kpath
        modk=modk+sqrt(dot_product(1.d0/100.d0*delta_kpath,1.d0/100.d0*delta_kpath))
        !
        call FT_r2q(kpt_path(j,:),Hktmp,H_hf)
        !
        ! Hktmp=0.d0;ek_out=0.d0
        ! Hktmp=get_HK(kpt_path(j,:),Nso,Hk_w90_reshape)
        call eigh(Hktmp,ek_out)
        write(uio,'(30F18.10)') modk,ek_out
        !
     end do
     !
  end do
  close(uio)

  
  stop
  











  !+- old stuff in momentum space -+!
  !
  
  !+- get Hk_w90 for each k-point -+!
  Hk_w90=0.d0
  do ik=1,Lk
     call get_Hk_w90(kpt_latt(ik,:),Hk_w90(:,:,ik))
     !
     do i=1,4
        Hk_w90(i,i,ik) = Hk_w90(i,i,ik) - cf_ext
     end do
     do i=5,6
        Hk_w90(i,i,ik) = Hk_w90(i,i,ik) + cf_ext
     end do
     !
  end do

  !+- compute number of particles -+!
  allocate(ccij(Nso,Nso))
  allocate(Hktmp(Nso,Nso))
  allocate(ek_out(Nso))
  Ncell=0.d0
  ccij=0.d0
  allocate(delta_hf(Nso,Nso,Lk),H_hf(Nso,Nso,Lk),delta_hf_(Nso,Nso,Lk))


  unit_in=free_unit(); 
  open(unit=unit_in,file='Ta_Ni_bare_vk.dat')    
  do ik=1,Lk     
     Hktmp=Hk_w90(:,:,ik)
     call eigh(Hktmp,ek_out)
     do i=1,Nso
        do j=1,Nso
           Ncell = Ncell + conjg(Hktmp(i,j))*Hktmp(i,j)*fermi(ek_out(j)-4.d0,beta)*wtk(ik)

           delta_hf(i,j,ik) = 0.d0
           do k=1,Nso
              delta_hf(i,j,ik) = delta_hf(i,j,ik) + conjg(Hktmp(i,k))*Hktmp(j,k)*fermi(ek_out(k)-4.d0,beta)*wtk(ik)           
              ccij(i,j) = ccij(i,j) + conjg(Hktmp(i,k))*Hktmp(j,k)*fermi(ek_out(k)-4.d0,beta)*wtk(ik)           
           end do

        end do
     end do
     write(unit_in,*) dble(ik),dreal(delta_hf(1,1,ik)),dreal(delta_hf(2,2,ik)),dreal(delta_hf(3,3,ik)),dreal(delta_hf(4,4,ik)),dreal(delta_hf(6,6,ik)),dreal(delta_hf(6,6,ik)),dreal(Hk_w90(1,6,ik)),dreal(Hk_w90(2,6,ik)),dreal(Hk_w90(3,6,ik)),dreal(Hk_w90(3,5,ik))
  end do
  close(unit_in)
  !
  unit_in=free_unit(); 
  open(unit=unit_in,file='Ta_Ni_bare_occupations.dat')    
  do i=1,Nso
     do j=1,Nso
        write(unit_in,*) i,j,dreal(ccij(i,j)),dimag(ccij(i,j))
     end do
     write(unit_in,*)
  end do
  close(unit_in)
  write(*,*) 'Ncell',Ncell





  do i=1,Nso
     do j=1,Nso
        write(768,'(10F18.10)') check_Uloc(i,j)
     end do
     check_Uloc(i,i) = 0.d0
  end do

  !+- here obtain the Fourier transform of the density density Coulomb interaction U(R-R') -+!  
  !+- TEST of THE FT -+!
  !
  Ur_TNS=0.d0
  iir=0
  do ir=1,6
     ! do jr=1,3
     !    do kr=1,3
     iir=iir+1
     !
     Rlat=dble(ir-2)*R1!+dble(jr-2)*R2+dble(kr-2)*R3
     call FT_q2r(Rlat,Ur_TNS(:,:,iir),Uq_TNS)
     modk=Rlat(1)**2.d0+Rlat(2)**2.d0+Rlat(3)**2.d0
     modk=sqrt(modk)           
     write(701,'(20F18.10)') modk,Ur_TNS(1,1,iir),Ur_TNS(1,2,iir),Ur_TNS(1,6,iir),Ur_TNS(1,5,iir)
     write(702,'(20F18.10)') Rlat,Ur_TNS(1,5,iir),Ur_TNS(1,6,iir),Ur_TNS(4,5,iir),Ur_TNS(5,4,iir)

     !
     !    end do
     ! end do
  end do

  Ur_TNS=0.d0
  iir=0
  do ir=1,Nk_z
     iir=iir+1
     !
     Rlat=dble(ir-2)*R3-1.d0*R1!+dble(jr-2)*R2+dble(kr-2)*R3
     call FT_q2r(Rlat,Ur_TNS(:,:,iir),Uq_TNS)
     modk=Rlat(1)**2.d0+Rlat(2)**2.d0+Rlat(3)**2.d0
     modk=sqrt(modk)           
     write(711,'(20F18.10)') modk,Ur_TNS(1,6,iir),Ur_TNS(6,1,iir),Ur_TNS(4,5,iir),Ur_TNS(5,4,iir)
     write(712,'(20F18.10)') Rlat,Ur_TNS(1,6,iir),Ur_TNS(6,1,iir),Ur_TNS(3,1,iir),Ur_TNS(1,3,iir)

     Rlat=dble(ir-2)*R3+0.d0*R1!+dble(jr-2)*R2+dble(kr-2)*R3
     call FT_q2r(Rlat,Ur_TNS(:,:,iir),Uq_TNS)
     modk=Rlat(1)**2.d0+Rlat(2)**2.d0+Rlat(3)**2.d0
     modk=sqrt(modk)           
     write(721,'(20F18.10)') modk,Ur_TNS(1,6,iir),Ur_TNS(6,1,iir),Ur_TNS(4,5,iir),Ur_TNS(5,4,iir)
     write(722,'(20F18.10)') Rlat,Ur_TNS(1,6,iir),Ur_TNS(6,1,iir),Ur_TNS(3,1,iir),Ur_TNS(1,3,iir)
     write(723,'(20F18.10)') Rlat,Ur_TNS(2,4,iir),Ur_TNS(4,2,iir)
     write(724,'(20F18.10)') Rlat,Ur_TNS(2,5,iir),Ur_TNS(5,2,iir)

  end do

  ! Ur_TNS=0.d0
  ! iir=0
  ! do ir=1,3
  !    do jr=1,3
  !       do kr=1,3
  !          iir=iir+1
  !          !
  !          Rlat=dble(ir-2)*R1+dble(jr-2)*R2+dble(kr-2)*R3
  !          call FT_q2r(Rlat,Ur_TNS(:,:,iir),Uq_TNS)
  !          modk=Rlat(1)**2.d0+Rlat(2)**2.d0+Rlat(3)**2.d0
  !          modk=sqrt(modk)           
  !          write(801,'(20F18.10)') modk,Ur_TNS(1,1,iir),Ur_TNS(1,2,iir),Ur_TNS(1,6,iir),Ur_TNS(1,5,iir)
  !    !
  !       end do
  !    end do
  ! end do
  
  
  




  stop


  delta_hf_=delta_hf

  unit_err=free_unit()
  open(unit=unit_err,file='err_q.err')


  !Uloc_TNS
  err_hf=1.d0  
  do ihf=1,Nhf
     write(unit_err,'(10F18.10)') dble(ihf),err_hf,mu_fix          
     call local_single_particle_observables(delta_hf,obs_loc)
     iso=0
     Ncell=0.d0
     do i=1,Nso           
        Ncell=Ncell+obs_loc(i)
        do j=i,Nso
           iso=iso+1
           write(units_loc_obs(iso),'(10F18.10)') dble(ihf),obs_loc(iso)
        end do
     end do
     write(unit_obs,'(10F18.10)') dble(ihf),ndens
     !
     call build_HF_hamiltonian_ij(H_hf,delta_hf,Uq_ij)          
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

  !
  call store_HF_hamiltonian_BZgrid_ij(Hhf_grid,delta_hf,mu_fix,Uq_ij)


  !+- here now I should plot the HF bands; see wheather the small hybridizations actually open a gap


  !+ M-point
  allocate(kpath(4,3))
  kpath(1,:)=0.5d0*Bk1+0.5d0*Bk3
  !+- Z-point
  kpath(2,:)=0.5d0*Bk3
  !+- G-point
  kpath(3,:)=0.d0
  !+- X-point
  kpath(4,:)=0.5d0*Bk1
  !
  color_bands=black
  kpoint_name(1)='M'
  kpoint_name(2)='Z'
  kpoint_name(3)='G'
  kpoint_name(4)='X'
  !
  allocate(kpt_path(300,3))
  allocate(kx(Lk),ky(Lk),kz(Lk))

  allocate(Hk_w90_reshape(Nso,Nso,Nk_x,Nk_y,Nk_z))
  Hk_w90_reshape=0.d0
  do ik=1,Lk
     i=ik_stride(ik,1)
     j=ik_stride(ik,2)
     k=ik_stride(ik,3)
     Hk_w90_reshape(:,:,i,j,k) = Hk_w90(:,:,ik)
     kx(ik) = kpt_latt(ik,1)
     ky(ik) = kpt_latt(ik,2)
     kz(ik) = kpt_latt(ik,3)
  end do
  modk=0.d0
  !



  allocate(Hk_hf_reshape(Nso,Nso,Nk_x,Nk_y,Nk_z))
  Hk_hf_reshape=0.d0
  do ik=1,Lk
     i=ik_stride(ik,1)
     j=ik_stride(ik,2)
     k=ik_stride(ik,3)
     Hk_hf_reshape(:,:,i,j,k) = Hhf_grid(:,:,ik)
     kx(ik) = kpt_latt(ik,1)
     ky(ik) = kpt_latt(ik,2)
     kz(ik) = kpt_latt(ik,3)
  end do
  modk=0.d0


  !+- compute spectral density -+!
  allocate(wreal(Lreal),G_loc(Nso,Lreal),Aw(Lreal))
  wreal=linspace(-5.0d0,5.0d0,Lreal)
  G_loc = 0.d0
  do ik=1,Lk
     Hktmp=0.d0;ek_out=0.d0     
     Hktmp=Hhf_grid(:,:,ik)
     call eigh(Hktmp,ek_out)
     do j=1,Nso     
        do i=1,Lreal
           G_loc(j,i) = G_loc(j,i) + 1.d0/(wreal(i)+xi*0.01d0-ek_out(j))*wtk(ik)
        end do
     end do
  end do
  Aw=0.d0
  do i=1,Nso
     Aw=Aw-1.d0/pi*dimag(G_loc(i,:))
  end do

  uio=free_unit()
  open(unit=uio,file='Aw.out')
  do i=1,Lreal
     write(uio,'(10F18.10)') wreal(i),Aw(i),-1.d0/pi*dimag(G_loc(:,i))
  end do
  close(uio)


  !
  !+- PLOTTING W90 BANDS -+!
  uio=free_unit()
  open(unit=uio,file='tns_HF_bands.out')
  !
  kmod=0.d0
  do i=1,3
     delta_kpath=kpath(i+1,:)-kpath(i,:)
     kmod = kmod + sqrt(dot_product(delta_kpath,delta_kpath))
     do ik=1,100
        j=(i-1)*100 + ik
        kpt_path(j,:) = kpath(i,:) + dble(ik-1)/100.d0*delta_kpath
        modk=modk+sqrt(dot_product(1.d0/100.d0*delta_kpath,1.d0/100.d0*delta_kpath))
        !
        Hktmp=0.d0;ek_out=0.d0
        Hktmp=get_HK(kpt_path(j,:),Nso,Hk_hf_reshape)
        call eigh(Hktmp,ek_out)
        write(uio,*) modk,ek_out
        !
     end do
     write(*,*) i,kmod(1)
     !
  end do
  close(uio)


  allocate(Hk_w90_tmp(Nso,Nso,Lk))
  Hk_w90_tmp=Hk_w90


  do i=1,2
     do j=3,6
        Hk_w90_tmp(i,j,:) = Hk_w90_tmp(i,j,:) + 0.10
        Hk_w90_tmp(j,i,:) = Hk_w90_tmp(j,i,:) + 0.10
     end do
  end do

  ! Hk_w90_tmp(5,1,:) = Hk_w90_tmp(5,1,:) + 0.10
  ! Hk_w90_tmp(1,6,:) = Hk_w90_tmp(1,6,:) + 0.10
  ! Hk_w90_tmp(6,1,:) = Hk_w90_tmp(6,1,:) + 0.10
  ! Hk_w90_tmp(2,5,:) = Hk_w90_tmp(2,5,:) + 0.10
  ! Hk_w90_tmp(5,2,:) = Hk_w90_tmp(5,2,:) + 0.10
  ! Hk_w90_tmp(3,6,:) = Hk_w90_tmp(3,6,:) + 0.10
  ! Hk_w90_tmp(6,3,:) = Hk_w90_tmp(6,3,:) + 0.10
  ! Hk_w90_tmp(4,6,:) = Hk_w90_tmp(4,6,:) + 0.10
  ! Hk_w90_tmp(6,4,:) = Hk_w90_tmp(6,4,:) + 0.10


  Hk_w90_reshape=0.d0
  do ik=1,Lk
     i=ik_stride(ik,1)
     j=ik_stride(ik,2)
     k=ik_stride(ik,3)
     Hk_w90_reshape(:,:,i,j,k) = Hk_w90_tmp(:,:,ik)
     kx(ik) = kpt_latt(ik,1)
     ky(ik) = kpt_latt(ik,2)
     kz(ik) = kpt_latt(ik,3)
  end do
  modk=0.d0
  !
  !+- PLOTTING W90 BANDS -+!
  uio=free_unit()
  open(unit=uio,file='tns_15_25_36_46_bands.out')
  !
  do i=1,3
     delta_kpath=kpath(i+1,:)-kpath(i,:)
     do ik=1,100
        j=(i-1)*100 + ik
        kpt_path(j,:) = kpath(i,:) + dble(ik-1)/100.d0*delta_kpath
        modk=modk+sqrt(dot_product(1.d0/100.d0*delta_kpath,1.d0/100.d0*delta_kpath))
        !
        Hktmp=0.d0;ek_out=0.d0
        Hktmp=get_HK(kpt_path(j,:),Nso,Hk_w90_reshape)
        call eigh(Hktmp,ek_out)
        write(uio,'(30F18.10)') modk,ek_out

        !
     end do
     !
  end do
  close(uio)



  Hk_w90_tmp=Hk_w90

  Hk_w90_tmp(1,6,:) = Hk_w90_tmp(1,6,:) + 0.10
  Hk_w90_tmp(6,1,:) = Hk_w90_tmp(6,1,:) + 0.10

  Hk_w90_tmp(3,5,:) = Hk_w90_tmp(3,5,:) + 0.10
  Hk_w90_tmp(5,3,:) = Hk_w90_tmp(5,3,:) + 0.10

  Hk_w90_tmp(4,5,:) = Hk_w90_tmp(4,5,:) + 0.10
  Hk_w90_tmp(5,4,:) = Hk_w90_tmp(5,4,:) + 0.10

  Hk_w90_tmp(6,2,:) = Hk_w90_tmp(6,2,:) + 0.10
  Hk_w90_tmp(2,6,:) = Hk_w90_tmp(2,6,:) + 0.10


  Hk_w90_reshape=0.d0
  do ik=1,Lk
     i=ik_stride(ik,1)
     j=ik_stride(ik,2)
     k=ik_stride(ik,3)
     Hk_w90_reshape(:,:,i,j,k) = Hk_w90_tmp(:,:,ik)
     kx(ik) = kpt_latt(ik,1)
     ky(ik) = kpt_latt(ik,2)
     kz(ik) = kpt_latt(ik,3)
  end do
  modk=0.d0
  !
  !+- PLOTTING W90 BANDS -+!
  uio=free_unit()
  open(unit=uio,file='tns_16_26_35_45_bands.out')
  !
  do i=1,3
     delta_kpath=kpath(i+1,:)-kpath(i,:)
     do ik=1,100
        j=(i-1)*100 + ik
        kpt_path(j,:) = kpath(i,:) + dble(ik-1)/100.d0*delta_kpath
        modk=modk+sqrt(dot_product(1.d0/100.d0*delta_kpath,1.d0/100.d0*delta_kpath))
        !
        Hktmp=0.d0;ek_out=0.d0
        Hktmp=get_HK(kpt_path(j,:),Nso,Hk_w90_reshape)
        call eigh(Hktmp,ek_out)
        write(uio,'(30F18.10)') modk,ek_out

        !
     end do
     !
  end do
  close(uio)






  !+- STILL TO DO: 

  !+- Do the fourier transform of the Interaction:
  !   for this I may need a larger real space lattice
  !   


  !+ Just do the mean field; need to include more than-one brillouin zone. 
  !  In 3d this could be annoying...maybe parallel integration could be considered...
  !  


  stop




contains

  !



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
