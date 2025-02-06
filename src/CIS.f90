subroutine CIS(nbas, no, nv, e, eri_mo,cis_energies)
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
double precision, allocatable :: c_array(:)
double precision, allocatable :: cc_array(:)
integer                       :: i, j, mu, nu, la, si, p, q, r, s, ex, ia, jb, a, b 
integer                       :: ex_vec(no * nv, 2)
integer                       :: identity(no + nv, no + nv)

! intent out
double precision, intent(out)  ::  cis_energies(no*nv)
!double precision, intent(out) :: 



ov = no * nv 
allocate(a_mat(ov, ov), c_array(ov), cc_array(ov))

write(*, *)
write(*, '(A65, i5, A65, i5)') 'The number of occupied orbitals is: ',&
        no, '. The number of virtual orbitals is: ', nv
write(*, *)

ex = 0 
do i = 1, no 
        do j = no + 1, no + nv  
                ex = ex + 1 
                ex_vec(ex, 1) = i 
                ex_vec(ex, 2) = j 
        end do 
end do 

write(*,'(A15)') '---------------'
write(*,'(A15)') 'CIS EXCITATIONS'
write(*,'(A15)') '---------------'
write(*, *) 
do i = 1, ov 
    write(*, '(A11, i4, a2, i4, a4, i4)') 'Excitation ', i, ': ', ex_vec(i,1), ' -> ', ex_vec(i, 2)
end do 


identity(:,:) = 0 

do i = 1, no + nv
        identity(i, i) = 1
end do 

write(*,*)
write(*,*) 'Building CIS matrix...'
write(*,*)

do ia = 1, ex
        do jb = 1, ex
                i = ex_vec(ia, 1)
                a = ex_vec(ia, 2)
                j = ex_vec(jb, 1)
                b = ex_vec(jb, 2)
                a_mat(ia, jb) = (e(a)-e(i)) * identity(i,j) * identity(a,b) + &
                       2*(eri_mo(i,b,a,j)) - eri_mo(i,b,j,a) 
        end do 
end do 


! write(*,'(A15)') '---------------'
! write(*,'(A15)') 'CIS MATRIX'
! write(*,'(A15)') '---------------'
! write(*, *) 
! 
! do i = 1, ov
!         write(*, '(*(f6.2))') a_mat(i,:)
! end do 
! 
! write(*, *) 


write(*,*) 'Diagonalizing CIS matrix...'
write(*,*)

call diagonalize_matrix(ov, a_mat, cis_energies)

write(*,'(A23)') '-----------------------'
write(*,'(A23)') 'CIS EXCITATION ENERGIES'
write(*,'(A23)') '-----------------------'
write(*, *) 

do i = 1, ov
        write(*, '(A6, i4, a2, f10.5, a13, f10.5, a3)')&
                'State ', i, ': ', cis_energies(i), ' hartree     ',cis_energies(i) * 27.2114, ' eV'
        c_array = a_mat(i,:)
        cc_array = a_mat(i,:)
        c_array = c_array ** 2
        C_array = c_array / norm2(c_array)
        do j = 1, ov
                if ( c_array(j) > 0.01d0 ) then 
                        write(*, '(i2, a4,i2,a2, f10.5, a6, f10.5, a2)') ex_vec(j,1),&
                                ' -> ', ex_vec(j, 2), ': ', c_array(j), ' ( c=', cc_array(j), ')'
                end if 
        end do 
        write(*,*)
end do 

write(*,*)

end subroutine CIS 

