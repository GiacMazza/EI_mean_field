MODULE HF
  USE VARS_GLOBAL
  USE SF_LINALG
  USE SF_SPECIAL
  USE SF_OPTIMIZE
  USE SF_IOTOOLS
  USE DMFT_VECTORS
  implicit none
  private

  public :: solve_HF_hamiltonian
  public :: build_HF_hamiltonian
  public :: find_chem_pot
  public :: check_conv
  public :: local_single_particle_observables
  public :: init_var_params

contains

  subroutine init_var_params(delta_hf)
    complex(8),dimension(Nso,Nso,Lk),intent(inout) :: delta_hf
    complex(8),dimension(Nso,Nso,Lk) :: Hhf        
    logical   :: IOfile
    integer   :: flen
    real(8)   :: mu_fix

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
  
  subroutine find_chem_pot(Hhf,deltas,mu)
    complex(8),dimension(Nso,Nso,Lk),intent(in) :: Hhf
    complex(8),dimension(Nso,Nso,Lk),intent(inout) :: deltas
    real(8),intent(inout) :: mu    
    real(8),dimension(1) :: mu_,Nout
    integer :: iter    
    call fsolve(delta_N,mu_,tol=1.d-10,info=iter)
    mu=mu_(1)
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
    integer :: ik,i,j,ii,jj,jk,iso,jso
    type(vect2D) :: deltaK
    real(8),dimension(Nso,Nso) :: Uh,Uf

    complex(8),dimension(Nso,Nso) :: Htest
    interface
       function Uq(q)
         USE DMFT_VECTORS  
         USE VARS_GLOBAL
         implicit none
         type(vect2D) :: q
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
    !Htest=0.d0
    Hhf=0.d0
    do ik=1,Lk
       Hhf(:,:,ik) = Hk(:,:,ik)       
       do jk=1,Lk
          !
          deltaK = k_bz(ik) - k_bz(jk)
          Uf=Uq(deltaK)
          Uh=Uq(Vzero)
          !
          do iso=1,Nso
             do jso=1,Nso
                Hhf(iso,jso,ik) = Hhf(iso,jso,ik) - &
                     Uf(iso,jso)*Deltas(iso,jso,jk)*wtk(jk)                     
                Hhf(iso,iso,ik) = Hhf(iso,iso,ik) + &
                     Uh(iso,jso)*Deltas(jso,jso,jk)*wtk(jk)
             end do
          end do
          !
       end do
       !Htest=Htest+Hhf(:,:,ik)*wtk(ik)
    end do

    ! write(*,*) Htest;
    ! stop

    !
  end subroutine build_HF_hamiltonian
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
