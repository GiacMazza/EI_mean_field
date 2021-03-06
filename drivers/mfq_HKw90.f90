program officina
  !
  USE SCIFOR
  USE DMFT_TOOLS
  !
  USE VARS_GLOBAL
  USE DMFT_VECTORS
  USE HF
  !
  !
  USE MPI
  !
  !
  implicit none


  integer :: Nint

  real(8) :: KKx,KKy
  complex(8),dimension(:,:,:),allocatable :: delta_hf,H_hf,delta_hf_,Hhf_grid
  complex(8),dimension(:,:,:),allocatable :: hr_w90
  real(8) :: mu_fix
  real(8) :: Uintra,Uinter,Delta_CF
  real(8) :: Uq,thop,err_hf
  real(8),dimension(3) :: R1,R2,R3
  real(8),dimension(3) :: Bk1,Bk2,Bk3
  !real(8),dimension(3,3) :: Bkinv

  real(8),dimension(3) :: ktest
  real(8),dimension(:),allocatable :: ek_out

  complex(8),dimension(:,:),allocatable :: Hsb

  complex(8),dimension(:),allocatable :: Uft

  complex(8),dimension(:,:),allocatable :: Hloc,Hktmp
  complex(8),dimension(:,:,:),allocatable :: Hk_w90,Hk_w90_tmp
  complex(8),dimension(:,:,:,:,:),allocatable :: Hk_w90_reshape,Hk_hf_reshape
  real(8),dimension(:,:),allocatable :: kpts,kpt_path


  integer :: iso,jso,ispin,jspin,iorb,jorb,ik
  integer :: i,j,k,idim
  integer :: Nhf,ihf,unit_err,unit_obs,unit_in,uio,Nobs
  integer,dimension(:),allocatable :: units_loc_obs
  integer :: flen,iread

  integer :: ir,jr,kr,iir,nrpts  
  real(8),dimension(3) :: Rlat 

  complex(8),dimension(:),allocatable :: obs_loc
  real(8) :: wmix
  real(8),dimension(:,:),allocatable :: Uloc_TNS
  complex(8),dimension(:,:,:),allocatable :: Uq_TNS,Vq_TNS,Ur_TNS
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


  integer,dimension(:,:),allocatable :: irvec

  integer,dimension(:),allocatable :: ndegen

  real(8) :: modk,tmpi,tmpj
  real(8) :: kmod(3)
  real(8) :: Ncell
  complex(8),dimension(:,:),allocatable :: ccij
  logical :: hartree
  real(8),allocatable,dimension(:) :: Ta_fat,Ni_fat
  integer :: ihyb_,jhyb_

  integer :: ix,iy,iz
  real(8),dimension(:),allocatable :: read_tmp
  real(8) :: cf_ext

  integer :: Lreal
  real(8),dimension(:),allocatable :: wreal
  complex(8),dimension(:,:),allocatable :: G_loc
  real(8),dimension(:),allocatable :: Aw

  logical :: Uradius

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
  !
  call get_global_vars
  !
  Nso=Norb*Nspin
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

  Lk_aux=3*3*3*Lk
  allocate(kpt_latt_aux(Lk_aux,3),ik_stride_aux(Lk_aux,3),igr2ik_aux(3*Nk_x,3*Nk_y,3*Nk_z))
  ik=0
  kpt_latt_aux=0.d0
  do i=1,3*Nk_x
     do k=1,3*Nk_z
        do j=1,3*Nk_y
           !
           ik=ik+1
           !
           ik_stride_aux(ik,1)=i
           ik_stride_aux(ik,2)=j
           ik_stride_aux(ik,3)=k
           igr2ik_aux(i,j,k) = ik
           !
           do idim=1,3
              kpt_latt_aux(ik,idim) = kpt_latt_aux(ik,idim) + kxgrid_aux(i)*Bk1(idim)
              kpt_latt_aux(ik,idim) = kpt_latt_aux(ik,idim) + kygrid_aux(j)*Bk2(idim)
              kpt_latt_aux(ik,idim) = kpt_latt_aux(ik,idim) + kzgrid_aux(k)*Bk3(idim)
           end do
           !
        end do
     end do
  end do
  !  
  uio=free_unit()
  open(unit=uio,file='mp_kgrid.f90')
  do ik=1,Lk
     write(uio,'(10F18.10)') kpt_latt(ik,:)
  end do
  close(uio)

  uio=free_unit()
  open(unit=uio,file='mp_kgrid_aux.f90')
  do ik=1,Lk_aux
     write(uio,'(10F18.10)') kpt_latt_aux(ik,:)
  end do
  close(uio)

  !+---------------------------------+!  

  file_name=reg(read_tns)//reg(file_w90_hr)  
  !+- read the w90 output -+!
  allocate(Hk_w90(Nso,Nso,Lk),Hloc(Nso,Nso))
  call read_w90_hr(R1,R2,R3,Hr_w90,Hloc,irvec,ndegen,trim(file_name),1,6,1)
  nrpts=size(irvec,1)
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

  mu_fix=4.d0
  !call find_chem_pot(Hk_w90,ccij,mu_fix)
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
  if(Nspin.eq.1) then
     do i=1,Nso
        Uloc_TNS(i,i) = 0.d0
     end do
  end if

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


  do ik=1,Lk
     !+- here print out Uq vs |q| -+!
     modk=kpt_latt(ik,1)**2.d0+kpt_latt(ik,2)**2.d0+kpt_latt(ik,3)**2.d0
     modk=sqrt(modk)
     write(601,'(20F18.10)') modk,Uq_TNS(1,1,ik),Uq_TNS(1,2,ik),Uq_TNS(1,6,ik),Uq_TNS(3,5,ik)
     ! write(501,'(20F18.10)') modk,check_Uloc(1,:)
     ! write(602,'(20F18.10)') modk,Uq_TNS(2,:,ik)
     ! write(603,'(20F18.10)') modk,Uq_TNS(3,:,ik)
     ! write(604,'(20F18.10)') modk,Uq_TNS(4,:,ik)
     ! write(605,'(20F18.10)') modk,Uq_TNS(5,:,ik)
     ! write(606,'(20F18.10)') modk,Uq_TNS(6,:,ik)
     do i=1,Nso
        ! Uq_TNS(i,i,ik)=Uq_TNS(i,i,ik)-check_Uloc(i,i)
        ! Vq_TNS(i,i,ik)=Vq_TNS(i,i,ik)-check_Uloc(i,i)
     end do
  end do
  !stop

  do i=1,Nso
     do j=1,Nso
        write(768,'(10F18.10)') check_Uloc(i,j)
     end do
     check_Uloc(i,i) = 0.d0
  end do

  !+- here obtain the Fourier transform of the density density Coulomb interaction U(R-R') -+!  
  !+- TEST of THE FT -+!
  allocate(rpt_latt(nrpts,3))
  allocate(Ur_TNS(Nso,Nso,nrpts))
  Ur_TNS=0.d0
  do ir=1,nrpts
     !
     rpt_latt(ir,:)=irvec(ir,1)*R1+irvec(ir,2)*R2+irvec(ir,3)*R3     
     call FT_q2r(rpt_latt(ir,:),Ur_TNS(:,:,ir),Uq_TNS)
     modk=rpt_latt(ir,1)**2.d0+rpt_latt(ir,2)**2.d0+rpt_latt(ir,3)**2.d0
     modk=sqrt(modk)
     write(501,'(20F18.10)') modk,Ur_TNS(1,1,ir),Ur_TNS(1,2,ir),Ur_TNS(1,6,ir),Ur_TNS(1,5,ir)
     write(511,'(20F18.10)') rpt_latt(ir,1:3),Ur_TNS(4,5,iir),Ur_TNS(5,4,iir),Ur_TNS(1,6,ir),Ur_TNS(6,1,ir)
     !
  end do
  ! Uq_TNS=0.d0
  ! do ik=1,Lk
  !    call FT_r2q(kpt_latt(ik,:),Uq_TNS(:,:,ik),Ur_TNS)
  !    modk=kpt_latt(ik,1)**2.d0+kpt_latt(ik,2)**2.d0+kpt_latt(ik,3)**2.d0
  !    modk=sqrt(modk)
  !    write(502,'(20F18.10)') modk,Uq_TNS(1,1,ik),Uq_TNS(1,2,ik),Uq_TNS(1,6,ik),Uq_TNS(3,5,ik)     
  ! end do
  
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
  
  
  Ur_TNS=0.d0
  do ir=1,nrpts
     rpt_latt(ir,:)=irvec(ir,1)*R1+irvec(ir,2)*R2+irvec(ir,3)*R3     
     modk=rpt_latt(ir,1)**2.d0+rpt_latt(ir,2)**2.d0+rpt_latt(ir,3)**2.d0
     !
     Rlat=dble(Nk_x*0.5d0)*R1+dble(Nk_y*0.5d0)*R2+dble(Nk_z*0.5d0)*R3
     Uradius=.true.
     Uradius=Uradius.and.(abs(rpt_latt(ir,1)).lt.abs(Rlat(1)))
     Uradius=Uradius.and.(abs(rpt_latt(ir,2)).lt.abs(Rlat(2)))
     Uradius=Uradius.and.(abs(rpt_latt(ir,3)).lt.abs(Rlat(3)))
     !     
     modk=sqrt(modk)
     if(modk.lt.20.0) then        
        call FT_q2r(rpt_latt(ir,:),Ur_TNS(:,:,ir),Uq_TNS)
        write(802,*) rpt_latt(ir,:),Rlat
     end if
     write(801,'(20F18.10)') modk,Ur_TNS(1,1,ir),Ur_TNS(1,2,ir),Ur_TNS(1,6,ir),Ur_TNS(1,5,ir)
  end do
  


  Uq_TNS=0.d0
  do ik=1,Lk
     call FT_r2q(kpt_latt(ik,:),Uq_TNS(:,:,ik),Ur_TNS)
     modk=kpt_latt(ik,1)**2.d0+kpt_latt(ik,2)**2.d0+kpt_latt(ik,3)**2.d0
     modk=sqrt(modk)
     write(502,'(20F18.10)') modk,Uq_TNS(1,1,ik),Uq_TNS(1,2,ik),Uq_TNS(1,6,ik),Uq_TNS(3,5,ik)     
  end do

  
  stop






  allocate(Hk(Nso,Nso,Lk)); Hk=Hk_w90  
  call init_var_params(delta_hf)
  do i=1,4
     do j=5,6
        delta_hf(i,j,:) = 0.10
        delta_hf(j,i,:) = 0.10
     end do
  end do



  delta_hf_=delta_hf

  unit_err=free_unit()
  open(unit=unit_err,file='err_q.err')

  Nobs=Nso*(Nso+1)/2
  allocate(units_loc_obs(Nobs))
  units_loc_obs=free_units(Nobs)
  iso=0
  do i=1,Nso
     do j=i,Nso
        iso=iso+1
        open(unit=units_loc_obs(iso),file="obs_hf_"//reg(txtfy(i)//reg(txtfy(j))))        
     end do
  end do
  unit_obs=free_unit()
  open(unit=unit_obs,file='ndens_hf.loop')

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

  !+- PLOTTING W90 BANDS -+!
  uio=free_unit()
  open(unit=uio,file='tns_bare_bands.out')
  !
  allocate(Ta_fat(Nso),Ni_fat(Nso))
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
        Ta_fat=abs(Hktmp(1,:))**2.d0+abs(Hktmp(2,:))**2.d0+abs(Hktmp(3,:))**2.d0+abs(Hktmp(4,:))**2.d0
        Ni_fat=abs(Hktmp(5,:))**2.d0+abs(Hktmp(6,:))**2.d0
        write(uio,'(30F18.10)') modk,ek_out,Ta_fat,Ni_fat
        write(234,'(30F18.10)') modk,fermi(ek_out-4.0,beta)

        !
     end do
     !
  end do
  close(uio)


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


  subroutine FT_q2r(Rlat,Fr,Fk) 
    implicit none
    real(8),dimension(3) :: kpoint
    complex(8),dimension(Nso,Nso) :: Fr
    complex(8),dimension(Nso,Nso,Lk) :: Fk
    integer :: ik
    real(8),dimension(3) :: Rlat
    real(8) :: dotRk
    !
    Fr=0.d0
    do ik=1,Lk
       kpoint=kpt_latt(ik,:)
       dotRk=dot_product(Rlat,kpoint)
       Fr=Fr+exp(-xi*dotRk)*Fk(:,:,ik)*wtk(ik)
    end do
  end subroutine FT_q2r

  subroutine FT_r2q(kpoint,Fk,Fr) 
    implicit none
    real(8),dimension(3) :: kpoint
    complex(8),dimension(Nso,Nso) :: Fk
    complex(8),dimension(Nso,Nso,nrpts) :: Fr
    integer :: ir
    real(8),dimension(3) :: Rlat
    real(8) :: dotRk
    !
    Fk=0.d0
    do ir=1,nrpts
       Rlat=irvec(ir,1)*R1+irvec(ir,2)*R2+irvec(ir,3)*R3
       dotRk=dot_product(Rlat,kpoint)
       Fk = Fk + Fr(:,:,ir)*exp(xi*dotRK)/dble(ndegen(ir))
    end do
    !
  end subroutine FT_r2q

  !
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
    allocate(kxgrid_aux(3*Nx),kygrid_aux(3*Ny),kzgrid_aux(3*Nz))
    allocate(ix_aux(3*Nx),iy_aux(3*Ny),iz_aux(3*Nz))
    !
    do ix=1,Nx
       !
       kxgrid_aux(ix) = kxgrid(ix) - 2*kxmax
       kxgrid_aux(ix+Nx) = kxgrid(ix) 
       kxgrid_aux(ix+2*Nx) = kxgrid(ix) + 2*kxmax
       !
       ix_aux(ix) = ix
       ix_aux(ix+Nx) = ix
       ix_aux(ix+2*Nx) = ix
       !
    end do
    do iy=1,Ny
       !
       kygrid_aux(iy) = kygrid(iy) - 2*kymax
       kygrid_aux(iy+Ny) = kygrid(iy) 
       kygrid_aux(iy+2*Ny) = kygrid(iy) + 2*kymax
       !
       iy_aux(iy) = iy
       iy_aux(iy+Ny) = iy
       iy_aux(iy+2*Ny) = iy
       !
    end do
    do iz=1,Nz
       !
       kzgrid_aux(iz) = kzgrid(iz) - 2*kzmax
       kzgrid_aux(iz+Nz) = kzgrid(iz) 
       kzgrid_aux(iz+2*Nz) = kzgrid(iz) + 2*kzmax
       !
       iz_aux(iz) = iz
       iz_aux(iz+Nz) = iz
       iz_aux(iz+2*Nz) = iz
       !
    end do

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
    unitk=free_unit()
    open(unit=unitk,file='mp_grid_aux.out')    
    do ix=1,3*Nx
       do iy=1,3*Ny
          do iz=1,3*Nz
             write(unitk,'(10F18.10)') kxgrid_aux(ix),kygrid_aux(iy),kzgrid_aux(iz)
          end do
       end do
    end do
    close(unitk)
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
    Nx_=Nx*N_cut_off
    Ny_=Ny*N_cut_off
    Nz_=Nz*N_cut_off

    allocate(kxr(Nx),kyr(Ny),kzr(Nz))
    kmax=dble(Nx_-1)*0.5d0/dble(Nx_)
    kxr = linspace(-kmax,kmax,Nx_)
    !
    kmax=dble(Ny_-1)*0.5d0/dble(Ny_)
    kyr = linspace(-kmax,kmax,Ny_)
    !
    kmax=dble(Nz_-1)*0.5d0/dble(Nz_)
    kzr = linspace(-kmax,kmax,Nz_)
    !
    if(mod(Nx_,2).eq.0) then
       kxr=kxr-0.5d0/dble(Nx_)
    end if
    if(mod(Nx_,2).eq.0) then
       kyr=kyr-0.5d0/dble(Ny_)
    end if
    if(mod(Nz_,2).eq.0) then
       kzr=kzr-0.5d0/dble(Nz_)
    end if
    !
    Lkr = Nx_*Ny_*Nz_
    allocate(krl(Lkr,3),wtk_rl(Lkr),ikrl2ii(Lkr,3))
    wtk_rl=1.d0/dble(Lkr)
    !
    !
    unitk=free_unit()
    open(unit=unitk,file='kRL_points.out')
    ik=0
    krl=0
    do ix=1,Nx_
       do iy=1,Nx_
          do iz=1,Nz_
             ik = ik + 1          
             !
             krl(ik,:)=krl(ik,:)+kxr(ix)*Bk1(:)
             krl(ik,:)=krl(ik,:)+kyr(iy)*Bk2(:)
             krl(ik,:)=krl(ik,:)+kzr(iz)*Bk3(:)
             !
             ikrl2ii(ik,1)=ix
             ikrl2ii(ik,2)=iy
             ikrl2ii(ik,3)=iz
             write(unitk,'(10F18.10)') krl(ik,:)
          end do
       end do
    end do
    close(unitk)
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


