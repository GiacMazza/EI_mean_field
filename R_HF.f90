MODULE HF_real
  USE VARS_GLOBAL
  USE SF_LINALG
  USE SF_SPECIAL
  USE SF_OPTIMIZE
  USE SF_IOTOOLS
  USE SF_CONSTANTS
  ! !USE DMFT_VECTORS

  implicit none
  private

  public :: solve_HF_hamiltonian
  public :: build_HF_hamiltonian_latt
  public :: find_chem_pot_latt
  public :: check_conv_latt
  public :: init_var_params_latt
  !
  public :: FT_r2q,FTr2q
  public :: FT_q2r,FTq2r
  !
  public :: fix_mu


  complex(8),dimension(:,:,:),allocatable :: Hhf_tmp,delta_hf_tmp
  logical :: tmp_print

contains

  function check_conv_latt(deltas_new,deltas_old) result(err)
    complex(8),dimension(Nso,Nso,nrpts),intent(in) :: deltas_new,deltas_old
    real(8) :: err    
    integer(8):: ir,iso,jso
    err=0.d0
    do ir=1,nrpts
       do iso=1,Nso
          do jso=1,iso
             err=err+abs(deltas_new(iso,jso,ir)-deltas_old(iso,jso,ir))/dble(nrpts)
          end do
       end do
    end do
  end function check_conv_latt


  subroutine init_var_params_latt(delta_hf,Hr_in)
    complex(8),dimension(Nso,Nso,nrpts),intent(inout) :: delta_hf
    complex(8),dimension(Nso,Nso,nrpts),intent(in) :: Hr_in    
    complex(8),dimension(Nso,Nso,nrpts) :: Hhf,Htmp,delta_hf_tmp        
    logical   :: IOfile
    integer(8)  :: flen
    real(8)   :: mu_fix
    real(8),dimension(3) :: ri,rj
    real(8) :: deltar
    integer :: unit,ir,jr
    !
    Hhf = Hr_in
    inquire(file=trim(init_Hsb),exist=IOfile)    
    if(IOfile) then
       flen = file_length(trim(init_Hsb))
       !
       if(flen.eq.Nso*Nso*Lk) then    
          write(*,*) 'reading Hsb from file'
          call read_array(trim(init_Hsb),Htmp)
          Hhf = Hhf + Htmp
       end if
       !       
       !
    end if
    !
    mu_fix=0.1d0
    call find_chem_pot_latt(Hhf,delta_hf_tmp,mu_fix)       
    !
    inquire(file=trim(init_HF),exist=IOfile)
    if(IOfile) then
       flen = file_length(trim(init_HF))
       if(flen.eq.Nso*Nso*Lk) then
          write(*,*) "reading HF variational parameters from file"
          call read_array(trim(init_HF),delta_hf)
       end if
    else
       delta_hf=delta_hf_tmp
    end if
    !
    unit=free_unit()
    open(unit=unit,file='ir_mirror.out')
    allocate(ir_mirror(nrpts))
    do ir=1,nrpts
       do jr=1,nrpts
          ri=irvec(ir,1)*R1+irvec(ir,2)*R2+irvec(ir,3)*R3
          ri=ri+(irvec(jr,1)*R1+irvec(jr,2)*R2+irvec(jr,3)*R3)
          deltar=sqrt(ri(1)**2.d0+ri(2)**2.d0+ri(3)**2.d0)
          if(abs(deltar).lt.1.d-12) then
             ir_mirror(ir) = jr
             if(ir.eq.jr) ir0=jr
             exit
          end if
       end do
       write(unit,*) ir,ir_mirror(ir),ir0,rpt_latt(ir,:),rpt_latt(ir_mirror(ir),:)
    end do
    close(unit)
  end subroutine init_var_params_latt
  
  
  subroutine find_chem_pot_latt(Hhf,delta_hf,mu)
    complex(8),dimension(Nso,Nso,nrpts),intent(in) :: Hhf
    complex(8),dimension(Nso,Nso,nrpts),intent(inout) :: delta_hf
    complex(8),dimension(Nso,Nso,Lk) :: Hhf_k,delta_hf_k
    real(8),intent(inout) :: mu    
    real(8),dimension(1) :: mu_,Nout
    integer :: iter,ik,ir
    ! real(8),dimension(3) :: Rlat
    ! real(8) :: dotRk
    !
    !+- transform the real-space Hamiltonian into a k-space one -+!
    !    
    do ik=1,Lk
       call FT_r2q(kpt_latt(ik,:),Hhf_k(:,:,ik),Hhf)
    end do
    ! Hhf_k=0.d0
    ! do ir=1,nrpts
    !    Rlat=irvec(ir,1)*R1+irvec(ir,2)*R2+irvec(ir,3)*R3
    !    dotRk=dot_product(Rlat,kpoint)
    !    Hk=Hk+Hr_w90(:,:,ir)*exp(xi*dotRK)/dble(ndegen(ir))
    ! end do
    !
    !call fsolve(delta_N,mu_,tol=1.d-10,info=iter)
    !mu=mu_(1)
    !Nout=delta_N(mu_)
    mu=brentq(delta_Nb,-100.d0,5000.d0)
    Nout=delta_Nb(mu)
    !
    !
    do ir=1,nrpts
       call FT_q2r(rpt_latt(ir,:),delta_hf(:,:,ir),delta_hf_k)
    end do
    !write(530,*) Nout,Ndens,mu_,iter
    !
  contains
    !
    function delta_N(xmu) result(deltaN)
      real(8),dimension(:) :: xmu
      real(8),dimension(size(xmu)) :: deltaN
      complex(8),dimension(Nso,Nso,Lk) :: Htmp
      integer(8):: ir,iso
      do ik=1,Lk
         Htmp(:,:,ik)=Hhf_k(:,:,ik)-xmu(1)*zeye(Nso)        
      end do
      call solve_HF_hamiltonian(Htmp,delta_hf_k) 
      deltaN=0.d0
      do ik=1,Lk
         do iso=1,Nso
            deltaN = deltaN + dreal(delta_hf_k(iso,iso,ik)*wtk(ik))
         end do
      end do
      !write(531,'(10F18.10)') deltaN,Ndens
      deltaN(1)=deltaN(1)-Ndens
    end function delta_N


    function delta_Nb(xmu) result(deltaN)
      real(8),intent(in) :: xmu
      real(8) :: deltaN
      complex(8),dimension(Nso,Nso,Lk) :: Htmp
      integer(8):: ir,iso
      do ik=1,Lk
         Htmp(:,:,ik)=Hhf_k(:,:,ik)-xmu*zeye(Nso)        
      end do
      call solve_HF_hamiltonian(Htmp,delta_hf_k) 
      deltaN=0.d0
      do ik=1,Lk
         do iso=1,Nso
            deltaN = deltaN + dreal(delta_hf_k(iso,iso,ik)*wtk(ik))
         end do
      end do
      !write(531,'(10F18.10)') deltaN,Ndens,xmu
      deltaN=deltaN-Ndens
    end function delta_Nb

  end subroutine find_chem_pot_latt

  subroutine solve_HF_hamiltonian(Hhf,Deltas,eout)
    complex(8),dimension(:,:,:),intent(in) :: Hhf
    complex(8),dimension(:,:,:),intent(inout) :: Deltas
    real(8),dimension(size(Hhf,1)) :: Ehf
    complex(8),dimension(size(Hhf,1),size(Hhf,1)) :: Htmp    
    integer(8):: ik,i,j,ii,jj
    integer :: iorb,jorb
    real(8) :: eout_
    real(8),optional :: eout
    !
    if(size(Hhf,1).ne.Nso) stop "error in Hhf1"
    if(size(Hhf,2).ne.Nso) stop "error in Hhf2"
    if(size(Hhf,3).ne.Lk) stop "error in Hhf3"
    !
    if(size(Deltas,1).ne.Nso) stop "error in Deltas1"
    if(size(Deltas,2).ne.Nso) stop "error in Deltas2"
    if(size(Deltas,3).ne.Lk) stop "error in Deltas3"
    !
    eout_=0.d0
    Deltas=0.d0             
    do ik=1,Lk
       !
       Htmp=Hhf(:,:,ik)
       call eigh(Htmp,Ehf)
       !
       do i=1,Nso
          eout_=eout_+Ehf(i)*fermi(Ehf(i),beta)*wtk(ik)
          do j=i,Nso
             ! do jj=1,Nso
             !    Deltas(j,i,ik) = Deltas(j,i,ik) + &
             !         conjg(Htmp(i,jj))*Htmp(j,jj)*fermi(Ehf(jj),beta)
             ! end do
             ! Deltas(i,j,ik)=conjg(Deltas(j,i,ik))
             do jj=1,Nso
                Deltas(i,j,ik) = Deltas(i,j,ik) + &
                     conjg(Htmp(i,jj))*Htmp(j,jj)*fermi(Ehf(jj),beta)
             end do
             Deltas(j,i,ik)=conjg(Deltas(i,j,ik))
          end do
       end do
       !
    end do
    if(present(eout)) eout=eout_
  end subroutine solve_HF_hamiltonian


  subroutine build_HF_hamiltonian_latt(Hhf,Delta_Hf,Ur)
    complex(8),dimension(Nso,Nso,nrpts),intent(inout) :: Hhf
    complex(8),dimension(Nso,Nso,nrpts),intent(in) :: Delta_HF
    complex(8),dimension(Nso,Nso,nrpts),intent(in) :: Ur
    !complex(8),dimension(Nso,Nso,nrpts) :: Hhf    
    integer(8):: ir,jr,iso,jso
    complex(8),dimension(Nso) :: Uhartree
    ! real(8),dimension(:),allocatable :: deltaK,deltak_
    ! real(8),dimension(Nso,Nso) :: Uh,Uf
    ! complex(8),dimension(Nso,Nso) :: Htest
    !Hhf=Hr_w90
    !+- 
    Uhartree = 0.d0
    do iso=1,Nso
       do ir=1,nrpts
          do jso=1,Nso
             Uhartree(iso) = Uhartree(iso) + Ur(iso,jso,ir)*Delta_HF(jso,jso,ir0)
          end do
       end do
    end do
    !+- 
    do iso=1,Nso
       !+- Hartree term -+!
       if(whartree) then
          Hhf(iso,iso,ir0) = Hhf(iso,iso,ir0) + Uhartree(iso)
       end if
       !+----------------+!
       do jso=1,Nso
          do ir=1,nrpts       
             !+- Fock Term -+!
             if(wfock) then
                !Hhf(iso,jso,ir) = Hhf(iso,jso,ir) - Ur(iso,jso,ir)*Delta_HF(jso,iso,ir_mirror(ir))
                Hhf(iso,jso,ir) = Hhf(iso,jso,ir) - Ur(iso,jso,ir)*Delta_HF(iso,jso,ir)
             end if
             !+-------------+!
          end do
       end do
    end do
    
    ! if(Nspin.eq.2) then
    !    do ir=1,nrpts
    !       write(345,*) Hhf(1,7,ir)
    !    end do
    ! end if
  end subroutine build_HF_hamiltonian_latt

  
  ! !
  ! subroutine observables(deltas,vec_obs)
  !   complex(8),dimension(Nso,Nso,Lk),intent(inout) :: Deltas    
  !   complex(8),dimension(:,:,:),allocatable :: vec_obs
  !   integer(8):: ik,iso,jso,ii
  !   !
  !   ! ii=Nso*(Nso+1)/2
  !   if(allocated(vec_obs)) deallocate(vec_obs)
  !   allocate(vec_obs(Nso,Nso,9))
  !   vec_obs=0.d0
  !   do ik=1,Lk
  !      ii=0
  !      do iso=1,Nso
  !         do jso=iso,Nso
  !            ii=ii+1
  !            vec_obs(ii) = vec_obs(ii)+deltas(iso,jso,ik)*wtk(ik)             
  !         end do
  !      end do
  !   end do
  !   !
  ! end subroutine observables











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
    complex(8),dimension(:,:) :: Fk
    complex(8),dimension(:,:,:) :: Fr
    integer :: ir
    real(8),dimension(3) :: Rlat
    real(8) :: dotRk
    !
    if(size(Fr,3).ne.nrpts) stop "(size(Fr,3).ne.nrpts)" 
    Fk=0.d0
    do ir=1,nrpts
       Rlat=irvec(ir,1)*R1+irvec(ir,2)*R2+irvec(ir,3)*R3
       dotRk=dot_product(Rlat,kpoint)
       Fk = Fk + Fr(:,:,ir)*exp(xi*dotRK)/dble(ndegen(ir))
    end do
    !
  end subroutine FT_r2q


  subroutine FTq2r(Rlat,Fr,Fk) 
    implicit none
    real(8),dimension(3) :: kpoint
    complex(8) :: Fr
    complex(8),dimension(Lk) :: Fk
    integer :: ik
    real(8),dimension(3) :: Rlat
    real(8) :: dotRk
    !
    if(size(Fk).ne.Lk) stop "(size(Fr,3).ne.nrpts)" 
    Fr=0.d0
    do ik=1,Lk
       kpoint=kpt_latt(ik,:)
       dotRk=dot_product(Rlat,kpoint)
       Fr=Fr+exp(-xi*dotRk)*Fk(ik)*wtk(ik)
    end do
  end subroutine FTq2r

  

  subroutine FTr2q(kpoint,Fk,Fr) 
    implicit none
    real(8),dimension(3) :: kpoint
    complex(8) :: Fk
    complex(8),dimension(:) :: Fr
    integer :: ir
    real(8),dimension(3) :: Rlat
    real(8) :: dotRk
    !
    if(size(Fr).ne.nrpts) stop "(size(Fr,3).ne.nrpts)" 
    Fk=0.d0
    do ir=1,nrpts
       Rlat=irvec(ir,1)*R1+irvec(ir,2)*R2+irvec(ir,3)*R3
       dotRk=dot_product(Rlat,kpoint)
       Fk = Fk + Fr(ir)*exp(xi*dotRK)/dble(ndegen(ir))
    end do
    !
  end subroutine FTr2q

  !



  ! subroutine build_HF_hamiltonian_ij(Hhf,Deltas,Uq)
  !   complex(8),dimension(:,:,:),intent(inout) :: Hhf
  !   complex(8),dimension(:,:,:),intent(in) :: Deltas    
  !   integer:: ik,i,j,ii,jj,kk,jk,iso,jso,jjk
  !   real(8),dimension(:),allocatable :: deltaK,deltak_
  !   complex(8),dimension(Nso,Nso) :: Uh,Uf
  !   complex(8),dimension(Nso,Nso) :: Hhartree

  !   real(8),dimension(3) :: deltak_target

  !   interface
  !      function Uq(iq)
  !        USE DMFT_VECTORS  
  !        USE VARS_GLOBAL
  !        implicit none
  !        !real(8),dimension(:) :: q
  !        integer :: iq
  !        complex(8),dimension(Nso,Nso) :: Uq
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
  !   allocate(deltak(size(kpt_latt,2)),deltak_(size(kpt_latt,2)))
  !   !
  !   Hhf=0.d0
  !   do ik=1,Lk
  !      Hhf(:,:,ik) = Hk(:,:,ik)       

  !      Hhartree=0.d0

  !      do jk=1,Lk !+- U(jk)*Delta(ik-jk)
  !         !
  !         deltak = kpt_latt(ik,:) + kpt_latt(jk,:)
  !         !+- find deltaK in the aux lattice -+!
  !         deltak_target=matmul(Bkinv,deltak)
  !         !
  !         do ii=1,3*Nk_x
  !            if(abs(deltak_target(1)-kxgrid_aux(ii)).lt.1.d-10) exit              
  !         end do
  !         do jj=1,3*Nk_y
  !            if(abs(deltak_target(2)-kygrid_aux(jj)).lt.1.d-10) exit              
  !         end do
  !         do kk=1,3*Nk_z
  !            if(abs(deltak_target(3)-kzgrid_aux(kk)).lt.1.d-10) exit              
  !         end do
  !         !
  !         Uf=Uq(jk)
  !         if(abs(kpt_latt(jk,1))**2.d0+abs(kpt_latt(jk,2))**2.d0+abs(kpt_latt(jk,3))**2.d0.lt.1.d-8) then
  !            Uf=0.d0
  !            Uh=0.d0
  !            if(whartree) then
  !               Uf=Uq(jk)
  !               Uh=Uq(jk)
  !            end if
  !         end if
  !         Uh=0.d0
  !         !
  !         !
  !         ii=ix_aux(ii)
  !         jj=iy_aux(jj)
  !         kk=iz_aux(kk)
  !         !
  !         jjk=igr2ik(ii,jj,kk)
  !         !
  !         do iso=1,Nso
  !            do jso=1,Nso                
  !               Hhf(iso,jso,ik) = Hhf(iso,jso,ik)- &
  !                    0.5d0*Uf(iso,jso)*Deltas(jso,iso,jjk)*wtk(jk)- &
  !                    0.5d0*Uf(iso,jso)*Deltas(iso,jso,jjk)*wtk(jk)                
  !            end do
  !            Hhartree(iso,iso) = Hhartree(iso,iso) + Deltas(iso,iso,jk)*wtk(jk)
  !         end do
  !         !          
  !      end do

  !      do iso=1,Nso       
  !         do jso=1,Nso
  !            Hhf(iso,iso,ik) = Hhf(iso,iso,ik) + Hhartree(iso,iso)*Uh(iso,jso)
  !         end do
  !      end do

  !   end do

  !   !
  ! end subroutine build_HF_hamiltonian_ij


  ! subroutine store_HF_hamiltonian_BZgrid(Hhf_grid,deltas,xmu,Uq)
  !   complex(8),dimension(:,:,:),allocatable :: Hhf_grid
  !   complex(8),dimension(Nso,Nso,Lk),intent(inout) :: Deltas    
  !   real(8) :: xmu
  !   integer(8):: ik
  !   interface
  !      function Uq(q)
  !        USE DMFT_VECTORS  
  !        USE VARS_GLOBAL
  !        implicit none
  !        real(8),dimension(:) :: q
  !        real(8),dimension(Nso,Nso) :: Uq
  !      end function Uq
  !   end interface
  !   !
  !   if(allocated(Hhf_grid)) deallocate(Hhf_grid)
  !   allocate(Hhf_grid(Nso,Nso,Lk)); Hhf_grid=0.d0
  !   call build_HF_hamiltonian(Hhf_grid,deltas,Uq)    
  !   do ik=1,Lk
  !      Hhf_grid(:,:,ik) = Hhf_grid(:,:,ik)-xmu*zeye(Nso)
  !   end do
  !   !
  ! end subroutine store_HF_hamiltonian_BZgrid


  ! subroutine store_HF_hamiltonian_BZgrid_ij(Hhf_grid,deltas,xmu,Uq)
  !   complex(8),dimension(:,:,:),allocatable :: Hhf_grid
  !   complex(8),dimension(Nso,Nso,Lk),intent(inout) :: Deltas    
  !   real(8) :: xmu
  !   integer(8):: ik
  !   interface
  !      function Uq(iq)
  !        USE DMFT_VECTORS  
  !        USE VARS_GLOBAL
  !        implicit none
  !        integer :: iq
  !        complex(8),dimension(Nso,Nso) :: Uq
  !      end function Uq
  !   end interface
  !   !
  !   if(allocated(Hhf_grid)) deallocate(Hhf_grid)
  !   allocate(Hhf_grid(Nso,Nso,Lk)); Hhf_grid=0.d0
  !   call build_HF_hamiltonian_ij(Hhf_grid,deltas,Uq)    
  !   do ik=1,Lk
  !      Hhf_grid(:,:,ik) = Hhf_grid(:,:,ik)-xmu*zeye(Nso)
  !   end do
  !   !
  ! end subroutine store_HF_hamiltonian_BZgrid_ij



  !+- K-SPACE ROUTINES -+!
  subroutine fix_mu(Hhf,delta_hf,mu,eout,iprint) 
    complex(8),dimension(Nso,Nso,Lk),intent(in) :: Hhf
    complex(8),dimension(Nso,Nso,Lk),intent(inout) :: delta_hf
    complex(8),dimension(Nso,Nso,Lk) :: Hhf_k,delta_hf_k
    complex(8),dimension(Nso,Nso,Lk) :: Htmp

    real(8),intent(inout) :: mu    
    real(8),dimension(1) :: mu_,Nout
    integer :: iter,ik,ir,im
    real(8) :: x1,x2
    real(8) :: mu1,mu2
    real(8) :: eout_
    real(8),optional :: eout
    logical,optional :: iprint
    !
    allocate(Hhf_tmp(Nso,Nso,Lk)); Hhf_tmp=Hhf
    allocate(delta_hf_tmp(Nso,Nso,Lk)); delta_hf_tmp=delta_hf

    tmp_print=.false.
    if(present(iprint)) tmp_print=iprint

    mu1=-10.d0
    mu2= 10.d0
    ! do im=1,100
    !    mu1=mu1-0.1
    !    mu2=mu2+0.1
    !    x1=deltaN_fix(mu1)
    !    x2=deltaN_fix(mu2)
    !    if(x1*x2.lt.0.d0) exit
    ! end do
    
    !write(*,*) x1,x2
    mu=brentq(deltaN_fix,mu1,mu2)
    Nout=deltaN_fix(mu)
    write(530,*) Nout,Ndens,mu

    !+- compute energy -+!
    if(present(eout)) then
       eout=0.d0
       do ik=1,Lk
          Htmp(:,:,ik)=Hhf_tmp(:,:,ik)-mu*zeye(Nso)
          !          write(987,'(10F18.10)') mu
       end do
       call solve_HF_hamiltonian(Htmp,delta_hf_tmp,eout)
    end if
    !+------------------+!
    delta_hf=delta_hf_tmp
    deallocate(delta_hf_tmp)
    deallocate(Hhf_tmp)
    !
  end subroutine fix_mu
  ! 
  function deltaN_fix(xmu) result(deltaN)
    real(8),intent(in) :: xmu
    real(8) :: deltaN
    complex(8),dimension(Nso,Nso,Lk) :: Htmp
    integer(8):: ik,iso
    do ik=1,Lk
       Htmp(:,:,ik)=Hhf_tmp(:,:,ik)-xmu*zeye(Nso)
    end do
    call solve_HF_hamiltonian(Htmp,delta_hf_tmp) 
    deltaN=0.d0
    do ik=1,Lk
       do iso=1,Nso
          deltaN = deltaN + dreal(delta_hf_tmp(iso,iso,ik)*wtk(ik))
       end do
    end do
    if(tmp_print) write(531,'(10F18.10)') deltaN,Ndens,xmu
    deltaN=deltaN-Ndens
  end function deltaN_fix
    





  ! subroutine HF_solver(Hhf,delta_hf,mu,eout,iprint) 
  !   complex(8),dimension(Nso,Nso,Lk),intent(in) :: Hhf
  !   complex(8),dimension(Nso,Nso,Lk),intent(inout) :: delta_hf
  !   complex(8),dimension(Nso,Nso,Lk) :: Hhf_k,delta_hf_k
  !   complex(8),dimension(Nso,Nso,Lk) :: Htmp

  !   real(8),intent(inout) :: mu    
  !   real(8),dimension(1) :: mu_,Nout
  !   integer :: iter,ik,ir,im
  !   real(8) :: x1,x2
  !   real(8) :: mu1,mu2
  !   real(8) :: eout_
  !   real(8),optional :: eout
  !   logical,optional :: iprint
  !   !
  !   allocate(Hhf_tmp(Nso,Nso,Lk)); Hhf_tmp=Hhf
  !   allocate(delta_hf_tmp(Nso,Nso,Lk)); delta_hf_tmp=delta_hf

  !   tmp_print=.false.
  !   if(present(iprint)) tmp_print=iprint

  !   mu1=-10.d0
  !   mu2= 10.d0
  !   ! do im=1,100
  !   !    mu1=mu1-0.1
  !   !    mu2=mu2+0.1
  !   !    x1=deltaN_fix(mu1)
  !   !    x2=deltaN_fix(mu2)
  !   !    if(x1*x2.lt.0.d0) exit
  !   ! end do
    
  !   !write(*,*) x1,x2
  !   mu=brentq(deltaN_fix,mu1,mu2)
  !   Nout=deltaN_fix(mu)
  !   write(530,*) Nout,Ndens,mu

  !   !+- compute energy -+!
  !   if(present(eout)) then
  !      eout=0.d0
  !      do ik=1,Lk
  !         Htmp(:,:,ik)=Hhf_tmp(:,:,ik)-mu*zeye(Nso)
  !         !          write(987,'(10F18.10)') mu
  !      end do
  !      call solve_HF_hamiltonian(Htmp,delta_hf_tmp,eout)
  !   end if
  !   !+------------------+!
  !   delta_hf=delta_hf_tmp
  !   deallocate(delta_hf_tmp)
  !   deallocate(Hhf_tmp)
  !   !
  ! end subroutine HF_solver
  ! ! 
  ! function deltaN_fix(xmu) result(deltaN)
  !   real(8),intent(in) :: xmu
  !   real(8) :: deltaN
  !   complex(8),dimension(Nso,Nso,Lk) :: Htmp
  !   integer(8):: ik,iso
  !   do ik=1,Lk
  !      Htmp(:,:,ik)=Hhf_tmp(:,:,ik)-xmu*zeye(Nso)
  !   end do
  !   call solve_HF_hamiltonian(Htmp,delta_hf_tmp) 
  !   deltaN=0.d0
  !   do ik=1,Lk
  !      do iso=1,Nso
  !         deltaN = deltaN + dreal(delta_hf_tmp(iso,iso,ik)*wtk(ik))
  !      end do
  !   end do
  !   if(tmp_print) write(531,'(10F18.10)') deltaN,Ndens,xmu
  !   deltaN=deltaN-Ndens
  ! end function deltaN_fix





  !
END MODULE HF_REAL
