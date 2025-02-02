subroutine CIS(nbas, no, nv, e, eri_mo,cis_energies)
implicit none 

! intent in 
integer, intent(in) :: nbas 
integer, intent(in) :: nv
integer, intent(in) :: no
double precision, intent(in) :: eri_mo(nbas, nbas, nbas, nbas)
double precision, intent(in) :: e(nbas)

! local variables 
integer, allocatable          :: identity_matrix(:, :)
integer                       :: ov
double precision, allocatable :: a_mat(:,:)
integer                       :: i, j, mu, nu, la, si, p, q, r, s, ex
integer                       :: ex_vec(no * nv, 2)

! intent out
double precision, intent(out)  ::  cis_energies(ov)
!double precision, intent(out) :: 



ov = no * nv 
allocate(identity_matrix(ov, ov), a_mat(ov, ov))

identity_matrix = 0

do i = 1, ov
    identity_matrix(i,i) = 1
end do 

do i = 1, ov
        write(*, '(*(i2))') identity_matrix(i,:)
end do 

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

end subroutine CIS 

