subroutine RHF(nBas,nO,S,T,V,Hc,ERI,X,ENuc,EHF,e,c)

! Perform a restricted Hartree-Fock calculation

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nBas

  integer,intent(in)            :: nO
  double precision,intent(in)   :: S(nBas,nBas)
  double precision,intent(in)   :: T(nBas,nBas)
  double precision,intent(in)   :: V(nBas,nBas)
  double precision,intent(in)   :: Hc(nBas,nBas) 
  double precision,intent(in)   :: X(nBas,nBas) 
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: ENuc

! Local variables

  integer,parameter             :: maxSCF = 64
  double precision,parameter    :: thresh = 1d-5
  integer                       :: nSCF, mu, nu, lambda, sigma, i
  double precision              :: Conv
  double precision              :: Gap
  double precision              :: ET,EV,EJ
  double precision              :: EK
  double precision,allocatable  :: cp(:,:)
  double precision,allocatable  :: P(:,:)
  double precision,allocatable  :: J(:,:)
  double precision,allocatable  :: K(:,:)
  double precision,allocatable  :: F(:,:),Fp(:,:), G(:,:)
  double precision,allocatable  :: error(:,:)
  double precision,external     :: trace_matrix

! Output variables

  double precision,intent(out)  :: EHF
  double precision,intent(out)  :: e(nBas)
  double precision,intent(out)  :: c(nBas,nBas)

! Hello world

  write(*,*)
  write(*,*)'************************************************'
  write(*,*)'|      Restricted Hartree-Fock calculation     |'
  write(*,*)'************************************************'
  write(*,*)

! Memory allocation

  allocate(cp(nBas,nBas),P(nBas,nBas),      &
           J(nBas,nBas),K(nBas,nBas),F(nBas,nBas),Fp(nBas,nBas), &
           error(nBas,nBas), g(nbas,nbas))

! Guess coefficients and eigenvalues

  F(:,:) = Hc(:,:)
  g(:,:) = 0.d0

! Initialization

  nSCF = 0
  Conv = 1d0
  P(:, :) = 0.d0 ! we start guessing that the density is 0 

!------------------------------------------------------------------------
! Main SCF loop
!------------------------------------------------------------------------

  write(*,*)
  write(*,*)'----------------------------------------------------'
  write(*,*)'| RHF calculation                                  |'
  write(*,*)'----------------------------------------------------'
  write(*,'(1X,A1,1X,A3,1X,A1,1X,A16,1X,A1,1X,A10,1X,A1,1X,A10,1X,A1,1X)') &
            '|','#','|','HF energy','|','Conv','|','HL Gap','|'
  write(*,*)'----------------------------------------------------'

  do while(Conv > thresh .and. nSCF < maxSCF)
  

      ! definition of g and F 
      F = Hc + G 
      
      ! Transform F to Fp 
      Fp = matmul(transpose(x), matmul(F, x))

      ! diagonalize Fp to get cp and epsilon
      call diagonalize_matrix(nbas, Fp, e)

      cp(:,:) = fp(:,:)

      ! Calculate C = X @ cp
      c = matmul(x, cp) 

      ! Calculate new P 
      p(:,:) = 2.0d0 * matmul(c(:,1:no), transpose(c(:,1:no)))

      ! Calculate gap 
      gap = abs(e(no) - e(no+1))


      !   Increment 
      nSCF = nSCF + 1

      ! recalculate F 
      g(:,:) = 0

      j(:,:) = 0.d0
      k(:,:) = 0.d0

      do mu = 1, nBas
          do nu = 1, nBas
              do lambda = 1, nBas
                  do sigma = 1, nBas
                      j(mu, nu) = j(mu, nu) + P(lambda, sigma) * eri(mu,lambda,nu,sigma)
                      k(mu, nu) = k(mu, nu) - 0.5d0 * P(lambda, sigma) * eri(mu,lambda,sigma,nu)
                      G(mu, nu) = G(mu, nu) +  P(lambda, sigma) * (eri(mu,lambda,nu,sigma) - 0.5d0 * eri(mu,lambda,sigma,nu))
                  enddo
              enddo
          enddo
      end do

      F = Hc + G

      ! commutator criterion
      conv = maxval(abs(matmul(f, matmul(p, s)) - matmul(s, matmul(p, f))))

      ! calculate HF energy
      EHF =  0.5d0 * trace_matrix(nbas, matmul(P, Hc+F))

      !   Dump results
      write(*,'(1X,A1,1X,I3,1X,A1,1X,F16.10,1X,A1,1X,F10.6,1X,A1,1X,F10.6,1X,A1,1X)') &
      '|',nSCF,'|',EHF+ENuc,'|',Conv,'|',Gap,'|'
 
      enddo
    
      write(*,*)'----------------------------------------------------'
!------------------------------------------------------------------------
! End of SCF loop
!------------------------------------------------------------------------

! Did it actually converge?

  if(nSCF == maxSCF) then

    write(*,*)
    write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    write(*,*)'                 Convergence failed                 '
    write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    write(*,*)

    stop

  endif

! Compute final HF energy

  ET = trace_matrix(nBas, matmul(p, T))
  EV = trace_matrix(nBas, matmul(p, V))
  EJ = 0.5d0 * trace_matrix(nBas, matmul(p, J))
  EK = 0.5d0 * trace_matrix(nBas, matmul(p, K))



  call print_RHF(nBas,nO,e,C,ENuc,ET,EV,EJ,EK,EHF)

end subroutine RHF
