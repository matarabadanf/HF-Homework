subroutine TDHF(nbas, no, nv, e, eri_mo, TDHF_energies)
implicit none 

! intent in 
integer, intent(in) :: nbas 
integer, intent(in) :: nv
integer, intent(in) :: no
double precision, intent(in) :: eri_mo(nbas, nbas, nbas, nbas)
double precision, intent(in) :: e(nbas)

! local variables 
integer                       :: ov
double precision, allocatable :: a_mat(:,:)
double precision, allocatable :: b_mat(:,:)
double precision, allocatable :: c_mat(:,:)
double precision, allocatable :: a_p_b(:,:)
double precision, allocatable :: a_m_b(:,:)
double precision, allocatable :: a_m_b_eigen(:)
double precision, allocatable :: a_m_b_r(:,:)
double precision, allocatable :: c_array(:)
double precision, allocatable :: cc_array(:)
integer                       :: i, j, mu, nu, la, si, p, q, r, s, ex, ia, jb, a, b 
integer                       :: ex_vec(no * nv, 2)
integer                       :: identity(no + nv, no + nv)
double precision, allocatable :: diag_a_m_b(:,:)

! intent out
double precision, intent(out)  ::  TDHF_energies(no*nv)
!double precision, intent(out) :: 



ov = no * nv 
allocate(a_mat(ov, ov), b_mat(ov, ov), c_mat(ov,ov), a_p_b(ov, ov),&
        a_m_b(ov, ov), a_m_b_r(ov, ov), a_m_b_eigen(ov), diag_a_m_b(ov,ov),&
        c_array(ov), cc_array(ov))

write(*, *)

ex = 0 
do i = 1, no 
        do j = no + 1, no + nv  
                ex = ex + 1 
                ex_vec(ex, 1) = i 
                ex_vec(ex, 2) = j 
        end do 
end do 

write(*,'(A16)') '----------------'
write(*,'(A16)') 'TDHF EXCITATIONS'
write(*,'(A16)') '----------------'
write(*, *) 
do i = 1, ov 
    write(*, '(A11, i4, a2, i4, a4, i4)') 'Excitation ', i, ': ', ex_vec(i,1), ' -> ', ex_vec(i, 2)
end do 


identity(:,:) = 0 

do i = 1, no + nv
        identity(i, i) = 1
end do 

write(*,*)
write(*,*) 'Constructing A matrix ...'

! Filling the A matrix
do ia = 1, ex
        do jb = 1, ex
                i = ex_vec(ia, 1)
                a = ex_vec(ia, 2)
                j = ex_vec(jb, 1)
                b = ex_vec(jb, 2)
                a_mat(ia, jb) = (e(a)-e(i)) * identity(i,j) * identity(a,b) + &
                       2.d0*(eri_mo(i,b,a,j)) - eri_mo(i,b,j,a) 
        end do 
end do 

write(*,*)
write(*,*) 'Constructing B matrix ...'
write(*,*)

! Filling the B matrix
do ia = 1, ex
        do jb = 1, ex
                i = ex_vec(ia, 1)
                a = ex_vec(ia, 2)
                j = ex_vec(jb, 1)
                b = ex_vec(jb, 2)
                b_mat(ia, jb) = 2.d0 * eri_mo(i,j,a,b) - eri_mo(i,j,b,a) 
        end do 
end do 

a_p_b = a_mat + b_mat 
a_m_b = a_mat - b_mat 

! evaluate square root of A-B to obtain the C matrix

! to do this we will use the symmetric diagonalization 
! S = U @ s @ U.T
! Then 
! S**0.5 = U @ s**0.5 @ U.T
! Therefore first we need the eigenvalues and eigenvectors of a_m_b

! write(*,*)
! write(*,'(A24)') 'A-B matrix:'
! write(*,*)
! 
! do i = 1, ov
!         write(*, '(*(f6.2))') a_m_b(i,:)
! end do 
! write(*,*)

write(*,*) 'Calculating (A-B)**0.5 ...'
write(*,*)

call diagonalize_matrix(ov, a_m_b, a_m_b_eigen)


! write(*,*)
! write(*,'(A24)') 'a_m_b eigenvectors:'
! write(*,*)
! 
! do i = 1, ov
!         write(*, '(*(f6.2))') a_m_b(i,:)
! end do 
! 
! 
! write(*,*)
! write(*,'(A24)') 'a_m_b eigenvalues:'
! write(*,*)
! write(*, '(*(f6.2))') a_m_b_eigen(:)
! write(*,*)
! write(*,'(A24)') 'a_m_b eigenvalues sqrt:'
! write(*,*)
! write(*, '(*(f6.2))') sqrt(a_m_b_eigen(:))
! write(*,*)

! with these, we obtain S**0.5 = U @ s**0.5 @ U.T

diag_a_m_b = 0.d0

do i = 1, ov
        diag_a_m_b(i,i) = a_m_b_eigen(i)
end do 

! a_m_b_r = matmul(a_m_b, matmul(diag_a_m_b, transpose(a_m_b)))
! write(*,*) 'Testing diagonalization, the next should be 0:'
! do i = 1, ov
!         write(*, '(*(f6.2))') a_mat(i,:) - b_mat(i,:) - a_m_b_r(i,:)
! end do 

a_m_b_r = matmul(a_m_b, matmul(sqrt(diag_a_m_b), transpose(a_m_b)))

! write(*,*)
! write(*,*) '(A-B)**0.5:'
! write(*,*)
! do i = 1, ov
!         write(*, '(*(f6.2))') a_m_b_r(i,:)
! end do 
! write(*,*)

write(*,*) 'Building C matrix ...'
write(*,*)

c_mat = matmul(a_m_b_r, matmul(a_p_b, a_m_b_r))

! write(*,'(A15)') '---------------'
! write(*,'(A15)') 'C MATRIX'
! write(*,'(A15)') '---------------'
! write(*, *) 
! 
! do i = 1, ov
!         write(*, '(*(f6.2))') c_mat(i,:)
! end do 

write(*,*) 'Diagonalizing C matrix ...'
write(*,*)

call diagonalize_matrix(ov, c_mat, TDHF_energies)

tdhf_energies = sqrt(tdhf_energies)

write(*,'(A24)') '------------------------'
write(*,'(A24)') 'TDHF EXCITATION ENERGIES'
write(*,'(A24)') '------------------------'
write(*, *) 

! do i = 1, ov
!         write(*, '((A11, i4, a2, i4, a4, i4), a2, (f10.5), a13, f10.5, a3)')&
!                 'Excitation ', i, ': ', ex_vec(i,1), ' -> ', ex_vec(i, 2), &
!                 ', ', TDHF_energies(i), ' hartree     ',TDHF_energies(i) * 27.2114, ' eV'
! end do 
! 
! write(*,*)
! ----------------------------------
do i = 1, ov
        write(*, '(A6, i4, a2, f10.5, a13, f10.5, a3)')&
                'State ', i, ': ', tdhf_energies(i), ' hartree     ',tdhf_energies(i) * 27.2114, ' eV'
        c_array = c_mat(i,:)
        cc_array = c_mat(i,:)
        c_array = c_array ** 2
        C_array = c_array / norm2(c_array)
        do j = 1, ov
                if ( c_array(j) > 0.005d0 ) then 
                        write(*, '(i2, a4,i2,a2, f10.5, a6, f10.5, a2)') ex_vec(j,1),&
                                ' -> ', ex_vec(j, 2), ': ', c_array(j), ' ( c=', cc_array(j), ')'
                end if 
        end do 
        write(*,*)
end do 
end subroutine TDHF

