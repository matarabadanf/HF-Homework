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

ex = 0 
do i = 1, no 
        do j = no + 1, no + nv  
                ex = ex + 1 
                ex_vec(ex, 1) = i 
                ex_vec(ex, 2) = j 
        end do 
end do 
        
write(*, '(*(i2))') ex_vec(:,:)


end subroutine CIS 

