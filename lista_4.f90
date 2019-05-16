include 'methods.f90'
include 'functions.f90'


program lista_4
    integer n, i, j
    real :: a, b, tol, y
    real lambda

    tol = 0.0001
    a = 0
    b = 10

    call bissecao(tol, a, b)
    !call f(b, y)

    write (*,*) a

end program lista_4