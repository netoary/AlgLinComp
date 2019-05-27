include 'methods.f90'
include 'functions.f90'


program lista_6
    use functions
    IMPLICIT NONE
    real :: xZero, tZero, deltaT, deltaT2, deltaT3, xLZero
    integer :: k, i
    real, allocatable, dimension(:) :: xK, xK2, xK3

    allocate(xK(3))
    allocate(xK2(2))
    allocate(xK3(2))

    xK = 0.0
    xZero = 1.0
    tZero = 0.0
    deltaT = 0.1
    deltaT2 = 0.1
    deltaT3 = 0.1
    k = 2
    xLZero = 0.0

    call euler(f1, xK(1), xZero, tZero, k, deltaT)
    tZero = 0.0
    xZero = 1.0
    call rungeKutta2(f1, xK(2), xZero, tZero, k, deltaT2)
    tZero = 0.0
    xZero = 1.0
    call rungeKutta4(f1, xK(3), xZero, tZero, k, deltaT3)

    do i=1, 4
        xK2 = 0.0
        xZero = 0.0
        tZero = 0.0
        deltaT = 1.0/(10**i)
        write(*,*) "delta"
        write(*,*) deltaT
        k = 100
        xLZero = 0.0

        call taylor2(f2, xK2(1), xZero, xLZero, tZero, k, deltaT)
        write(*,*) "taylor"
        write(*,*) xK2(1)
        xZero = 0.0
        tZero = 0.0
        call rungeKuttaNystrom(f2, xK2(2), xZero, xLZero, tZero, k, deltaT)
        write(*,*) "runge Kutta Nystrom"
        write(*,*) xK2(2)
    end do

    xK3 = 0.0
    xZero = 0.0
    tZero = 0.0
    deltaT = 0.1
    k = 20
    xLZero = 0.0

    call taylor2(f3, xK3(1), xZero, xLZero, tZero, k, deltaT)
    xZero = 0.0
    tZero = 0.0
    xLZero = 0.0
    call rungeKuttaNystrom(f3, xK3(2), xZero, xLZero, tZero, k, deltaT)

    write (*,*) xK
    write (*,*) xK2
    write (*,*) xK3

contains
    real function f1(x, t) result(y)
        real, intent(in) :: x, t
        y = -2*t*(x**2)
    end function

    real function f2(x, t, xL) result(y)
        real, intent(in) :: x, t, xL
        y = 2*sin(0.5*t) + sin(2*0.5*t) + cos(3*0.5*t) - 0.2*xL - x
        end function

    real function f3(x, t, xL) result(y)
        real, intent(in) :: x, t, xL
        y = - 9.81 - (1 * xL * abs(xL))
    end function


end program lista_6