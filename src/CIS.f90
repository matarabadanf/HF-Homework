subroutine CIS 
implicit none 

! intent in 
integer, intent(in) :: nbas 
integer, intent(in) :: nbas 
integer, intent(in) :: ho 

! local variables 
integer, allocatable          :: identity_matrix(:, :)
! intent out
double precision, intent(out) ::
double precision, intent(out) :: 



ov = no * nv 
allocate(identity_matrix(ov, ov), a_matrix(ov, ov), b_matrix(ov, ov), a_plus_b(ov, ov), a_minus_n(ov, ov))

identity_matrix = 0

do i = 1, ov
    identity_matrix(i,i) = 0
end do 

end subroutine CIS 

