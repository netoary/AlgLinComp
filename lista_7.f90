include 'methods.f90'
include 'functions.f90'


program lista_7
    use functions
    IMPLICIT NONE
    real :: x1, x2, x3, deltaX, deltaX1, deltaX2, p
    real, allocatable, dimension(:) :: fLinha1, fLinha2, fLinha3

    allocate(fLinha1(4))
    allocate(fLinha2(4))
    allocate(fLinha3(4))
    deltaX = 0.001
    deltaX1 = 0.5
    deltaX2 = 0.25
    x1 = 3
    x2 = 2
    x3 = 6
    p = 1
    fLinha1 = 0.0
    fLinha2 = 0.0
    fLinha3 = 0.0

    call diferencaCentral(f1, x1, deltaX, fLinha1(1))
    call passoFrente(f1, x1, deltaX, fLinha1(2))
    call passoTras(f1, x1, deltaX, fLinha1(3))
    call interpolacaoRichard(f1, x1, deltaX1, deltaX2, fLinha1(4), p)

    call diferencaCentral(f2, x2, deltaX, fLinha2(1))
    call passoFrente(f2, x2, deltaX, fLinha2(2))
    call passoTras(f2, x2, deltaX, fLinha2(3))
    call interpolacaoRichard(f2, x2, deltaX1, deltaX2, fLinha2(4), p)

    call diferencaCentral(f3, x3, deltaX, fLinha3(1))
    call passoFrente(f3, x3, deltaX, fLinha3(2))
    call passoTras(f3, x3, deltaX, fLinha3(3))
    call interpolacaoRichard(f3, x3, deltaX1, deltaX2, fLinha3(4), p)

    write (*,*) fLinha1
    write (*,*) fLinha2
    write (*,*) fLinha3

contains
    real function f1(x) result(y)
        real, intent(in) :: x
        y = x**3 + 1/exp(x)
    end function

    real function f2(x) result(y)
        real, intent(in) :: x
        y = x**(1/3) + log(x)
    end function

    real function f3(x) result(y)
        real, intent(in) :: x
        y = 1 - exp(-(x/5)**2)
    end function
end program lista_7