subroutine orthogonalization_matrix(number_of_basis, overlap_matrix, x_matrix)
        implicit none 
        
        integer, intent(in) :: number_of_basis
        double precision, intent(in) :: overlap_matrix(number_of_basis,number_of_basis)
        double precision, intent(out) :: x_matrix(number_of_basis, number_of_basis)

        integer :: i 
        double precision :: s_eigenvalues(number_of_basis)                         ! vector de autovalores de s 
        double precision :: s_diagonal(number_of_basis, number_of_basis)           ! matriz u 
        double precision :: identity(number_of_basis, number_of_basis)             ! matriz u 
        double precision :: s_minuscula(number_of_basis, number_of_basis)          ! matriz s minuscula
        double precision :: s_minuscula_un_medio(number_of_basis, number_of_basis) ! matriz s minuscula elevada a un medio
        ! diagonalize S 

        s_diagonal(:,:) = overlap_matrix(:,:)

        s_minuscula(:,:) = 0.d0
        s_minuscula_un_medio(:,:) = 0.d0

        call diagonalize_matrix(number_of_basis, s_diagonal, s_eigenvalues)        
       
        ! generate the diagonal matrix

        do i = 1, number_of_basis
                s_minuscula(i,i) = s_eigenvalues(i)
        enddo
        
        ! power the diagonal matrix

        do i = 1, number_of_basis
                s_minuscula_un_medio(i,i) = 1.d0 / sqrt(s_eigenvalues(i))
        enddo
        

        ! multiply the matrices X @ s**0.5 @ X.T 

        x_matrix = matmul(s_diagonal , matmul(s_minuscula_un_medio, transpose(s_diagonal)))

        identity = matmul(X_matrix, matmul(overlap_matrix, transpose(x_matrix)))


        open(42, access='append', file='matrices.dat')

        write(42, *)
        write(42, *) 'S matrix after diagonalization'
        write(42, *)
   
        do i = 1, number_of_basis
                write(42, '(*(f16.6))') s_minuscula(i,:)
        enddo

        write(42, *)
        write(42, *) 'Identity matrix'
        write(42, *)

        do i = 1, number_of_basis
                write(42, '(*(f16.6))') identity(i,:)
        enddo
        
        
        write(42, *)
        write(42, *) 'S 1/2 matrix '
        write(42, *)

        do i = 1, number_of_basis
                write(42, '(*(f16.6))') s_minuscula_un_medio(i,:)
        enddo

        close(42)

end subroutine
