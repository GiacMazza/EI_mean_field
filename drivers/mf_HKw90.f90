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
  complex(8),dimension(:,:,:),allocatable :: delta_hf,H_hf,delta_hf_
  complex(8),dimension(:,:,:),allocatable :: hr_w90
  real(8) :: mu_fix
  real(8) :: Uintra,Uinter,Delta_CF
  real(8) :: Uq,thop,err_hf
  real(8),dimension(3) :: R1,R2,R3
  real(8),dimension(3) :: Bk1,Bk2,Bk3
  real(8),dimension(3,3) :: Bkinv
  
  real(8),dimension(3) :: ktest
  real(8),dimension(:),allocatable :: ek_out

  complex(8),dimension(:,:),allocatable :: Hsb

  complex(8),dimension(:),allocatable :: Uft

  complex(8),dimension(:,:),allocatable :: Hloc,Hktmp
  complex(8),dimension(:,:,:),allocatable :: Hk_w90
  complex(8),dimension(:,:,:,:,:),allocatable :: Hk_w90_reshape
  real(8),dimension(:,:),allocatable :: kpts,kpt_path
  

  integer :: iso,jso,ispin,jspin,iorb,jorb,ik
  integer :: i,j,k,idim
  integer :: Nhf,ihf,unit_err,unit_obs,unit_in,uio,Nobs
  integer,dimension(:),allocatable :: units_loc_obs
  integer :: flen,iread

  integer,allocatable :: ik_stride(:,:)

  complex(8),dimension(:),allocatable :: obs_loc
  real(8) :: wmix
  real(8),dimension(:,:),allocatable :: Uloc_TNS

  character(len=100) :: fileRlattice

  integer,allocatable :: Nkvect(:)


  real(8),dimension(:),allocatable :: kx,ky,kz
  !
  real(8),dimension(:),allocatable :: kxgrid_aux,kygrid_aux,kzgrid_aux
  integer(8),dimension(:),allocatable :: ix_aux,iy_aux,iz_aux

  real(8),dimension(:,:),allocatable :: kpath
  real(8) ::      delta_kpath(3)
  type(rgb_color),dimension(2) :: color_bands
  character(len=1),dimension(4) :: kpoint_name


  integer,dimension(:,:),allocatable :: irvec

  integer,dimension(:),allocatable :: ndegen

  real(8) :: modk,tmpi,tmpj
  real(8) :: Ncell
  complex(8),dimension(:,:),allocatable :: ccij
  !

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
  call parse_input_variable(Nint,"Nint","input.conf",default=5)
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
  call TB_set_ei(R1,R2,R3)
  call TB_get_bk(Bk1,Bk2,Bk3)
  Bkinv(:,1) = Bk1
  Bkinv(:,2) = Bk2
  Bkinv(:,3) = Bk3
  call inv(Bkinv)
  !+- build a monkhorst-pack grid -+!  
  call build_mp_grid(Nk_x,Nk_y,Nk_z)
  Lk=Nk_x*Nk_y*Nk_z
  allocate(kpt_latt(Lk,3),ik_stride(Lk,3),wtk(Lk),igr2ik(Nk_x,Nk_y,Nk_z))
  wtk=1.d0/dble(Lk)
  ik=0
  do i=1,Nk_x
     do j=1,Nk_y
        do k=1,Nk_z           
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
  
  uio=free_unit()
  open(unit=uio,file='mp_kgrid.f90')
  do ik=1,Lk
     write(uio,'(10F18.10)') kpt_latt(ik,:)
  end do
  close(uio)
  !+---------------------------------+!  

  !+- read the w90 output -+!
  allocate(Hk_w90(Nso,Nso,Lk),Hloc(Nso,Nso))
  call read_w90_hr(R1,R2,R3,Hr_w90,Hloc,irvec,ndegen,'TNS_hr.dat',1,6,1)
  !+- get Hk_w90 for each k-point -+!
  Hk_w90=0.d0
  do ik=1,Lk
     call get_Hk_w90(kpt_latt(ik,:),Hk_w90(:,:,ik))
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
     write(unit_in,*) dble(ik),dreal(Hk_w90(1,6,ik)),dreal(Hk_w90(2,6,ik)),dreal(Hk_w90(3,6,ik)),dreal(Hk_w90(3,5,ik))
  end do
  close(unit_in)

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





  allocate(Hk(Nso,Nso,Lk)); Hk=Hk_w90  
  call init_var_params(delta_hf)
  ! delta_hf(1,:,:)=delta_hf(1,:,:)+0.1d0
  ! delta_hf(:,1,:)=delta_hf(:,1,:)+0.1d0
  delta_hf=delta_hf+0.1
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
     call build_HF_hamiltonian(H_hf,delta_hf,Uq_jellium)          
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
  !call store_HF_hamiltonian_BZgrid(Hhf_grid,delta_hf,mu_fix,Uq_jellium)


  !+- here now I should plot the HF bands; see wheather the small hybridizations actually open a gap


  stop



  !
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
  modk=0.d0
  !
  !+- PLOTTING W90 BANDS -+!
  uio=free_unit()
  open(unit=uio,file='tns_bands.out')
  !
  do i=1,3
     delta_kpath=kpath(i+1,:)-kpath(i,:)
     do ik=1,100
        j=(i-1)*100 + ik
        kpt_path(j,:) = kpath(i,:) + dble(ik-1)/100.d0*delta_kpath
        modk=modk+sqrt(dot_product(1.d0/100.d0*delta_kpath,1.d0/100.d0*delta_kpath))
        !
        Hktmp=0.d0;ek_out=0.d0
        Hktmp=get_HK(kpt_path(j,:),Nso)
        call eigh(Hktmp,ek_out)
        write(uio,*) modk,ek_out
        !
     end do
     !
  end do
  

  !+- STILL TO DO: 

  !+- Do the fourier transform of the Interaction:
  !   for this I may need a larger real space lattice
  !   


  !+ Just do the mean field; need to include more than-one brillouin zone. 
  !  In 3d this could be annoying...maybe parallel integration could be considered...
  !  


  stop




contains


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


  function get_HK(kpoint,N) result(Hk_out)
    implicit none
    real(8),dimension(:) :: kpoint
    integer              :: N
    complex(8),dimension(N,N) :: Hk_out
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
    ! ktarget(1)=ktarget(1)!-dble(Nk_x-1)*0.5d0/dble(Nk_x)
    ! ktarget(2)=ktarget(2)!-dble(Nk_y-1)*0.5d0/dble(Nk_y)
    ! ktarget(3)=ktarget(3)!-dble(Nk_z-1)*0.5d0/dble(Nk_z)
    !write(*,'(20F10.5)') kpoint,kpt_latt(ik_target,:),ktarget,Hk_w90(1,1,ik_target)
    ! write(*,*) 'closest k',kpt_latt(ik_target,:),Hk_w90(1,1,ik_target)
    ! write(*,*) 'k target',ktarget
    !stop
    !
    !+- get the *centered* index of the aux lattice -+!
    i0=ik_stride(ik_target,1)+Nk_x
    j0=ik_stride(ik_target,2)+Nk_y
    k0=ik_stride(ik_target,3)+Nk_z
    !
    !+- some temporary check
    ! write(251,*) kxgrid(i0),kxgrid_aux(i0+Nk_x)
    ! write(252,*) kygrid(j0),kygrid_aux(j0+Nk_y)
    ! write(253,*) kzgrid(k0),kzgrid_aux(k0+Nk_z)
    ! !
    ! do ix=1,3*Nk_x
    !    write(254,*) kxgrid_aux(ix),kxgrid(ix_aux(ix))
    ! end do
    ! do ix=1,3*Nk_y
    !    write(255,*) kygrid_aux(ix),kygrid(iy_aux(ix))
    ! end do
    ! do ix=1,3*Nk_z
    !    write(256,*) kzgrid_aux(ix),kzgrid(iz_aux(ix))
    ! end do

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

             kpt_int(ikint,:)=kix(i)*Bk1+ kiy(j)*Bk2 + kiz(k)*Bk3
             !
             Hk_inter(:,:,i,j,k)=Hk_w90_reshape(:,:,ix_aux(ii),iy_aux(jj),iz_aux(kk)) 
             !
             ! write(350,*) kix(i),kiy(j),kiz(k),ktarget!,kpoint
             ! write(353,*) kpt_int(ikint,:)
             !
             ! if(k==Nint+1) then
             !    write(450,*) kpt_int(ikint,1:2),dreal(Hk_inter(1,1,i,j,k))
             ! end if
             !
          end do
       end do
    end do
    !
    
    ! +- a lot of confusion, I needed to do some test to be confortable w/ interpolation
    !
    
    ! allocate(yout(Nint_),yout_(Nint_),hk_tmp(Nint_))
    
    ! allocate(grid_aux(Nint_))

    ! io=1;jo=1   
    ! do k=1,Nint_
    !    do j=1,Nint_          
          
    !       !+- interpolate along the BK1 direction -+!

    !       !+- take values of Hamiltonian to interpolate -+!
    !       hk_tmp(:) = dreal(Hk_inter(io,jo,:,j,k))          
          
    !       !+- get grid for the interpolation -+!
    !       do i=1,Nint_
    !          !ktmp=kiy(j)*Bk2+kiz(k)*Bk3
    !          !ktmp=ktmp+kix(i)*Bk1 
    !          grid_aux(i)=dot_product(ktmp,Bk1)
    !       end do          
    !       ! !+- get the target point -+!
    !       ! ktmp=ktarget(1)*Bk1+kiy(j)*Bk2+kiz(k)*Bk3
    !       ! ktmp=dot_product(ktmp,Bk1)
          
    !       call polint(kix,hk_tmp,ktmp,yout(j),dH)          
    !    end do       
    !    !+- interpolate along the BK2 direction -+!
       
    !    hk_tmp=yout !+- take previous results as new interpolation points

    !    !+- get grid for the interpolation -+!
    !    ! do i=1,Nint_
    !    !    ktmp=kiy(j)*Bk2+kiz(k)*Bk3
    !    !    ktmp=ktmp+kix(i)*Bk1 
    !    !    grid_aux(i)=dot_product(ktmp,Bk1)
    !    ! end do
       



    !    call polint(kiy,hk_tmp,ktarget(2),yout_(k),dH)       
    ! end do
    ! !+- interpolate along the z axis -+!
    ! hk_tmp=yout_
    ! call polint(kiz,hk_tmp,ktarget(3),Hr,dH)
    ! ! write(*,*) kiz
    ! ! write(*,*) hk_tmp
    ! ! !write(*,*) ktarget(3)
    ! ! write(*,*) 'test_inter',Hr
    ! stop

    ! io=1
    ! jo=1
    ! allocate(Hkr(Nint_,Nint_,Nint_))    
    ! Hkr=dreal(Hk_inter(io,jo,:,:,:))    
    ! call polin3(kix,kiy,kiz,Hkr,ktarget(1),ktarget(2),ktarget(3),Hr,dH)    
    ! write(*,*) 'test_inter polin3',Hr
    
    !stop
    
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
    ! Uq_jellium = 0.d0
    ! if(mod.ne.0.d0) then
    !    !Uq_jellium = Umat_loc + Umat/modulo(q)
    !    Uq_jellium = Umat/mod
    ! end if
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


