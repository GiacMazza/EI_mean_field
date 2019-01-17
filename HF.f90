MODULE HF
  USE VARS_GLOBAL
  USE SF_LINALG
  USE SF_SPECIAL
  USE SF_OPTIMIZE
  USE SF_IOTOOLS
  USE SF_CONSTANTS
  !USE DMFT_VECTORS
  implicit none
  private

  public :: solve_HF_hamiltonian
  public :: build_HF_hamiltonian
  public :: find_chem_pot
  public :: check_conv
  public :: local_single_particle_observables
  public :: init_var_params
  public :: store_HF_hamiltonian_BZgrid

  public :: shift_BZ

contains

  subroutine shift_BZ(ik_in,ik_out)
    integer :: ik_in,ik_out
    integer :: ix_out,iy_out,ix,iy
    !
    !
    ix=ikrl2ii(ik_in,1)
    iy=ikrl2ii(ik_in,2)
    !
    !
    ix_out=ix; 
!    write(*,*) 'ix',ix_out,ix_out-Nx,Nx
    do while(ix_out.gt.Nx) 
       ix_out=ix_out-Nx
    end do

    iy_out=iy
!    write(*,*) 'iy',iy_out,iy_out-Nx,Nx
    do while(iy_out.gt.Nx) 
       iy_out=iy_out-Nx
    end do
!    write(*,*) ix_out,iy_out
    ik_out = igr2ik(ix_out,iy_out)

  end subroutine shift_BZ


  subroutine init_var_params(delta_hf)
    complex(8),dimension(Nso,Nso,Lk),intent(inout) :: delta_hf
    complex(8),dimension(Nso,Nso,Lk) :: Hhf        
    logical   :: IOfile
    integer   :: flen
    real(8)   :: mu_fix
    
    delta_hf=0.d0

    Hhf = Hk
    inquire(file=trim(init_Hsb),exist=IOfile)    
    if(IOfile) then
       flen = file_length(trim(init_Hsb))
       !
       if(flen.eq.Nso*Nso*Lk) then    
          write(*,*) 'reading Hsb from file'
          call read_array(trim(init_Hsb),Hhf)
          Hhf = Hhf + Hk
       end if
       !       
       !
    end if
    !
    mu_fix=0.1d0
    call find_chem_pot(Hhf,delta_hf,mu_fix)       
    !
    inquire(file=trim(init_HF),exist=IOfile)
    if(IOfile) then
       flen = file_length(trim(init_HF))
       if(flen.eq.Nso*Nso*Lk) then
          write(*,*) "reading HF variational parameters from file"
          call read_array(trim(init_HF),delta_hf)
       end if
    end if
    !
  end subroutine init_var_params



  function check_conv(deltas_new,deltas_old) result(err)
    complex(8),dimension(Nso,Nso,Lk),intent(in) :: deltas_new,deltas_old
    real(8) :: err    
    integer :: ik,iso,jso
    err=0.d0
    do ik=1,Lk
       do iso=1,Nso
          do jso=1,iso
             err=err+abs(deltas_new(iso,jso,ik)-deltas_old(iso,jso,ik))*wtk(ik)
          end do
       end do
    end do
  end function check_conv

  subroutine find_chem_pot(Hhf,delta_hf,mu)
    complex(8),dimension(Nso,Nso,Lk),intent(in) :: Hhf
    complex(8),dimension(Nso,Nso,Lk),intent(inout) :: delta_hf
    complex(8),dimension(Nso,Nso,Lk) :: deltas
    real(8),intent(inout) :: mu    
    real(8),dimension(1) :: mu_,Nout
    integer :: iter    
    !
    deltas=delta_hf
    call fsolve(delta_N,mu_,tol=1.d-10,info=iter)
    mu=mu_(1)
    Nout=delta_N(mu_)
    delta_hf=deltas
    write(530,'(10F18.10)') Nout
    !
  contains
    !
    function delta_N(xmu) result(deltaN)
      real(8),dimension(:) :: xmu
      real(8),dimension(size(xmu)) :: deltaN
      complex(8),dimension(Nso,Nso,Lk) :: Htmp
      integer :: ik,iso
      do ik=1,Lk
         Htmp(:,:,ik)=Hhf(:,:,ik)-xmu(1)*zeye(Nso)        
      end do
      call solve_HF_hamiltonian(Htmp,deltas) 
      deltaN=0.d0
      do ik=1,Lk
         do iso=1,Nso
            deltaN = deltaN + deltas(iso,iso,ik)*wtk(ik)
         end do
      end do
      deltaN(1)=deltaN(1)-Ndens
    end function delta_N
  end subroutine find_chem_pot

  subroutine solve_HF_hamiltonian(Hhf,Deltas)
    complex(8),dimension(:,:,:),intent(in) :: Hhf
    complex(8),dimension(:,:,:),intent(inout) :: Deltas
    real(8),dimension(size(Hhf,1)) :: Ehf
    complex(8),dimension(size(Hhf,1),size(Hhf,1)) :: Htmp    
    integer :: ik,i,j,ii,jj
    !
    if(size(Hhf,1).ne.Nso) stop "error in Hhf1"
    if(size(Hhf,2).ne.Nso) stop "error in Hhf2"
    if(size(Hhf,3).ne.Lk) stop "error in Hhf3"
    !
    if(size(Deltas,1).ne.Nso) stop "error in Deltas1"
    if(size(Deltas,2).ne.Nso) stop "error in Deltas2"
    if(size(Deltas,3).ne.Lk) stop "error in Deltas3"
    !
    do ik=1,Lk
       !
       Htmp=Hhf(:,:,ik)
       call eigh(Htmp,Ehf)
       !
       do i=1,Nso
          do j=1,i
             Deltas(i,j,ik)=0.d0             
             do jj=1,Nso
                Deltas(i,j,ik) = Deltas(i,j,ik) + &
                     conjg(Htmp(i,jj))*Htmp(j,jj)*fermi(Ehf(jj),beta)
             end do
             Deltas(j,i,ik)=conjg(Deltas(i,j,ik))
          end do
       end do
       !
    end do
  end subroutine solve_HF_hamiltonian
  
  subroutine build_HF_hamiltonian(Hhf,Deltas,Uq)
    complex(8),dimension(:,:,:),intent(inout) :: Hhf
    complex(8),dimension(:,:,:),intent(in) :: Deltas    
    integer :: ik,i,j,ii,jj,jk,iso,jso,jjk
    real(8),dimension(:),allocatable :: deltaK,deltak_
    real(8),dimension(Nso,Nso) :: Uh,Uf
    complex(8),dimension(Nso,Nso) :: Htest
    interface
       function Uq(q)
         USE DMFT_VECTORS  
         USE VARS_GLOBAL
         implicit none
         real(8),dimension(:) :: q
         real(8),dimension(Nso,Nso) :: Uq
       end function Uq
    end interface
    !
    if(size(Hhf,1).ne.Nso) stop "error in Hhf1"
    if(size(Hhf,2).ne.Nso) stop "error in Hhf2"
    if(size(Hhf,3).ne.Lk) stop "error in Hhf3"
    !
    if(size(Deltas,1).ne.Nso) stop "error in Deltas1"
    if(size(Deltas,2).ne.Nso) stop "error in Deltas2"
    if(size(Deltas,3).ne.Lk) stop "error in Deltas3"
    !
    allocate(deltak(size(k_bz,2)),deltak_(size(k_bz,2)))
    !
    Hhf=0.d0
    do ik=1,Lk
       Hhf(:,:,ik) = Hk(:,:,ik)       
       do jk=1,Lkr
          !
          deltaK = k_bz(ik,:) - krl(jk,:)
          !
          Uf=Uq(deltaK)
          deltaK=0.d0
          Uh=Uq(deltak)
          !
          call shift_BZ(jk,jjk)
          !
          do iso=1,Nso
             do jso=1,Nso
                Hhf(iso,jso,ik) = Hhf(iso,jso,ik) - &
                     Uf(iso,jso)*Deltas(iso,jso,jjk)*wtk_rl(jk)                     
                Hhf(iso,iso,ik) = Hhf(iso,iso,ik) + &
                     Uh(iso,jso)*Deltas(jso,jso,jjk)*wtk_rl(jk)
             end do
          end do
          !
       end do
    end do
    !
  end subroutine build_HF_hamiltonian


  ! subroutine build_HF_hamiltonian(Hhf,Deltas,Uq)
  !   complex(8),dimension(:,:,:),intent(inout) :: Hhf
  !   complex(8),dimension(:,:,:),intent(in) :: Deltas    
  !   integer :: ik,i,j,ii,jj,jk,iso,jso
  !   real(8),dimension(:),allocatable :: deltaK
  !   real(8),dimension(Nso,Nso) :: Uh,Uf
  !   complex(8),dimension(Nso,Nso) :: Htest
  !   interface
  !      function Uq(q)
  !        USE DMFT_VECTORS  
  !        USE VARS_GLOBAL
  !        implicit none
  !        !type(vect2D) :: q
  !        real(8),dimension(:) :: q
  !        real(8),dimension(Nso,Nso) :: Uq
  !      end function Uq
  !   end interface
  !   !
  !   if(size(Hhf,1).ne.Nso) stop "error in Hhf1"
  !   if(size(Hhf,2).ne.Nso) stop "error in Hhf2"
  !   if(size(Hhf,3).ne.Lk) stop "error in Hhf3"
  !   !
  !   if(size(Deltas,1).ne.Nso) stop "error in Deltas1"
  !   if(size(Deltas,2).ne.Nso) stop "error in Deltas2"
  !   if(size(Deltas,3).ne.Lk) stop "error in Deltas3"
  !   !
  !   !Htest=0.d0
    
  !   allocate(deltak(size(k_bz,2)))

  !   Hhf=0.d0
  !   do ik=1,Lk
  !      Hhf(:,:,ik) = Hk(:,:,ik)       
  !      do jk=1,Lk
  !         !
  !         deltaK = k_bz(ik,:) - k_bz(jk,:)
  !         Uf=Uq(deltaK)
  !         deltaK=0.d0
  !         Uh=Uq(deltak)
  !         !
  !         do iso=1,Nso
  !            do jso=1,Nso
  !               Hhf(iso,jso,ik) = Hhf(iso,jso,ik) - &
  !                    Uf(iso,jso)*Deltas(iso,jso,jk)*wtk(jk)                     
  !               Hhf(iso,iso,ik) = Hhf(iso,iso,ik) + &
  !                    Uh(iso,jso)*Deltas(jso,jso,jk)*wtk(jk)
  !            end do
  !         end do
  !         !
  !      end do
  !   end do
  !   !
  ! end subroutine build_HF_hamiltonian



  ! subroutine print_HF_bands(Deltas,xmu,Uq)
  !   complex(8),dimension(:,:,:),intent(inout) :: Deltas
  !   complex(8),dimension(Nso,Nso,Lk) :: Hhf
  !   real(8),dimension(size(Deltas,1)) :: Ehf
  !   complex(8),dimension(size(Hhf,1),size(Hhf,1)) :: Htmp    
  !   integer :: ik,i,j,ii,jj
  !   interface
  !      function Uq(q)
  !        USE DMFT_VECTORS  
  !        USE VARS_GLOBAL
  !        implicit none
  !        type(vect2D) :: q
  !        real(8),dimension(Nso,Nso) :: Uq
  !      end function Uq
  !   end interface

  !   !
  !   if(size(Hhf,1).ne.Nso) stop "error in Hhf1"
  !   if(size(Hhf,2).ne.Nso) stop "error in Hhf2"
  !   if(size(Hhf,3).ne.Lk) stop "error in Hhf3"
  !   !
  !   if(size(Deltas,1).ne.Nso) stop "error in Deltas1"
  !   if(size(Deltas,2).ne.Nso) stop "error in Deltas2"
  !   if(size(Deltas,3).ne.Lk) stop "error in Deltas3"
  !   !

  !   call build_HF_hamiltonian

  !   do ik=1,Lk
  !      !
  !      Htmp=Hhf(:,:,ik)
  !      call eigh(Htmp,Ehf)
  !      !
  !      do i=1,Nso
  !         do j=1,i
  !            Deltas(i,j,ik)=0.d0             
  !            do jj=1,Nso
  !               Deltas(i,j,ik) = Deltas(i,j,ik) + &
  !                    conjg(Htmp(i,jj))*Htmp(j,jj)*fermi(Ehf(jj),beta)
  !            end do
  !            Deltas(j,i,ik)=conjg(Deltas(i,j,ik))
  !         end do
  !      end do
  !      !
  !   end do
  ! end subroutine print_HF_bands

  subroutine store_HF_hamiltonian_BZgrid(Hhf_grid,deltas,xmu,Uq)
    complex(8),dimension(:,:,:),allocatable :: Hhf_grid
    complex(8),dimension(Nso,Nso,Lk),intent(inout) :: Deltas    
    real(8) :: xmu
    integer :: ik
    interface
       function Uq(q)
         USE DMFT_VECTORS  
         USE VARS_GLOBAL
         implicit none
         real(8),dimension(:) :: q
         real(8),dimension(Nso,Nso) :: Uq
       end function Uq
    end interface
    !
    if(allocated(Hhf_grid)) deallocate(Hhf_grid)
    allocate(Hhf_grid(Nso,Nso,Lk)); Hhf_grid=0.d0
    call build_HF_hamiltonian(Hhf_grid,deltas,Uq)    
    do ik=1,Lk
       Hhf_grid(:,:,ik) = Hhf_grid(:,:,ik)-xmu*zeye(Nso)
    end do
    !
  end subroutine store_HF_hamiltonian_BZgrid



  !
  subroutine local_single_particle_observables(deltas,vec_obs)
    complex(8),dimension(Nso,Nso,Lk),intent(inout) :: Deltas    
    complex(8),dimension(:),allocatable :: vec_obs
    integer :: ik,iso,jso,ii
    !
    ii=Nso*(Nso+1)/2
    if(allocated(vec_obs)) deallocate(vec_obs)
    !
    allocate(vec_obs(ii))
    vec_obs=0.d0
    do ik=1,Lk
       ii=0
       do iso=1,Nso
          ii=ii+1
          vec_obs(ii) = vec_obs(ii)+deltas(iso,iso,ik)*wtk(ik)
       end do
       do iso=1,Nso
          do jso=iso+1,Nso
             ii=ii+1
             vec_obs(ii) = vec_obs(ii)+deltas(iso,jso,ik)*wtk(ik)
          end do
       end do
    end do
    !
  end subroutine local_single_particle_observables
  !
END MODULE HF
