include 'methods.f90'

program gauss
implicit none
    integer :: n, m, i, j
    real, allocatable, dimension(:) :: b, x
    real, allocatable, dimension(:,:) :: a
    Character(len = 100) :: fileName


    print*,"Digite o nome do arquivo com o sistema a ser resolvido:"
    read(*,*) fileName
    open (1, file=fileName, status='old', action='read')

     read(1,*) n, m
     allocate(a(n,m))

     do i=1, m
         read(1,*) ( a(j,i) , j=1 , n)
     end do
     allocate(b(n))
     read(1,*) ( b(i) , i= 1 , n)
    close(1)


     allocate(x(n))

     call eliminacaoGauss(a, b, x, n)

     print*, "x = ", x
     
     open(2, file='RESUL_' // fileName //'.txt', status='replace')
     write(2, *) "A"
     do i=1, n
         write(2,*) (a(j,i), j=1, n)
     end do
     write(2, *) "A^-1"
     write(2, *) "x"
     write(2, *) (x(i), i=1, n)
     close(2)

    call sleep(10000)

end program gauss