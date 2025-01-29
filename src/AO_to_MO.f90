subroutine AO_to_MO(nBas,c,ERI_AO,ERI_MO)

! Expression of bi-electronic integrals in the MO basis set

  implicit none
!  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nBas
  double precision,intent(in)   :: c(nBas,nBas)
  double precision,intent(in)   :: ERI_AO(nBas,nBas,nBas,nBas)

! Local variables

  integer                       :: mu,nu,la,si
  integer                       :: p,q,r,s
  double precision,allocatable  :: scr(:,:,:,:)

! Output variables

  double precision,intent(out)  :: ERI_MO(nBas,nBas,nBas,nBas)

! Memory allocation
  
!--------------------------------------
! AO to MO transformation starts here !
!--------------------------------------

ERI_MO(:,:,:,:) = 0.d0

! Brief reminder that the integral (pq|rs) = <pr|qs>. In the previous suborutine 
! we have used the notation <pr|qs>. Therefore we have to consider the transformation
! of the equation 4 to (pq|rs) = \sum c c c c <mu lambda | nu sigma>

! do mu = 1, nBas
!  do la = 1, nBas
!   do nu = 1, nBas
!    do si = 1, nBas
!     do p = 1, nbas   
!      do q = 1, nbas   
!       do r = 1, nbas   
!        do s = 1, nbas  
!            ERI_MO(p,r,q,s) = ERI_MO(p,r,q,s) + c(mu,p) * c(nu,q) * c(la,r) * c(si,s) * ERI_AO(mu,la,nu,si)! saasdfasdf
!        end do 
!       end do 
!      end do 
!     end do 
!    end do 
!   end do 
!  end do 
! end do 

!write(*,*) 'First approach'
!write(*,*) ERI_mo(1,1,1,:)

! This approach is quite costly as the number of basis functions increases as it has a computational cost of O(nbas**8)


! A better approach is to do the transformations stepwise. 

! We do then the following transformations:  < mu la | nu si> -> <mu la | nu s> -> <mu la | q s > -> < mu r | q s > -> < p r | q s >

! With this transformation instead of O(nbas**8) we get 4 * O(nbas**5)

allocate(scr(nbas,nbas,nbas,nbas))

scr(:,:,:,:) = 0.d0
eri_mo(:,:,:,:) = 0.d0

! For this we want to transform < mu la | nu si> -> <mu la | nu s> first 

do mu = 1, nBas
 do la = 1, nbas
  do nu = 1, nbas
   do s = 1, nbas
    do si = 1, nbas
     scr(mu, la, nu, s) = scr(mu,la,nu,s) + c(si, s) * eri_ao(mu, la, nu, si)
    end do 
   end do
  end do
 end do 
end do

! For this we want to transform < mu la | nu s > -> <mu la | q s >  now 
eri_mo(:,:,:,:) = 0.d0 
do mu = 1, nBas
 do q = 1, nbas
  do la = 1, nbas
   do s = 1, nbas
    do nu = 1, nbas
     eri_mo(mu, la, q, s) = eri_mo(mu, la, q, s) + c(nu, q) * scr(mu, la, nu, s)
    end do 
   end do
  end do
 end do 
end do

scr = eri_mo

! The third is to transform < mu la | q s > -> < mu r | q s >  now 
eri_mo = 0.d0

do q = 1, nBas
 do r = 1, nbas
  do mu = 1, nbas
   do s = 1, nbas
    do la = 1, nbas
     eri_mo(mu, r, q, s) = eri_mo(mu, r, q, s) + c(la, r) * scr(mu, la, q, s)
    end do 
   end do
  end do
 end do 
end do

scr = eri_mo

! Lastly we transform < mu r | q s > -> < p r | q s >  now 
eri_mo = 0.d0
do p = 1, nbas
 do r = 1, nBas
  do q = 1, nbas
   do s = 1, nbas
    do mu = 1, nbas
     eri_mo(p, r, q, s) = eri_mo(p, r, q, s) + c(mu, p) * scr(mu, r, q, s)
    end do
   end do
  end do
 end do 
end do

!write(*,*) 'second approach'
!write(*,*) ERI_mo(1,1,1,:)

end subroutine AO_to_MO
