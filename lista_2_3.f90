include 'methods.f90'

program lista_2_3
implicit none
    integer :: n, m, i, j, methodId
    real, allocatable, dimension(:) :: x, b
    real, allocatable, dimension(:,:) :: a
    Character(len = 100) :: fileName
    
    print*,"Digite o nome do arquivo com o sistema a ser resolvido:"
    read(*,*) fileName
    print*,"Escolha o método de solução:"
    print*,"[1] Jacobi"
    print*,"[2] Gauss-Seide"
    
    read(*,*) methodId
    

    open (1, file=fileName, status='old', action='read')
    read(1, *) n, m
    allocate(a(n,m))
    do i=1, m
        read(1, *) ( a(i,j) , j= 1 , n)
    end do
    allocate(b(n))
    read(1, *) ( b(i) , i= 1 , n)
    close(1)

    if (methodId == 1) then
        call interativoJacobi(a, b, x, n)
    else if (methodId == 2) then
        call GaussSeidel(a, b, x, n)
    else
        print *, "Método Inválido"
        call sleep(10000)
        stop
    end if

    open(2, file='RESUL.txt', status='replace')
    do i=1, n
        write(2,*) (a(i, j), j=1, n)
    end do
    write(2, *) (x(i), i=1, n)
    close(2)

    call sleep(10000)
    
end program lista_2_3