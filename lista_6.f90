include 'methods.f90'
include 'functions.f90'


program lista_6
    use functions
    IMPLICIT NONE
    real :: xK, xZero, tZero, deltaT, xLZero
    integer :: k

    xK = 0.0
    xZero = 1.0
    tZero = 0.0
    deltaT = 0.1
    k = 2
    xLZero = 0.0

    call euler(f1, xK, xZero, tZero, k, deltaT)
    !call rungeKutta2(xK, xZero, tZero, k, deltaT)
    !call rungeKutta4(xK, xZero, tZero, k, deltaT)

    !call taylor2(xK, xZero, xLZero, tZero, k, deltaT)
    !call rungeKuttaNystrom(xK, xZero, xLZero, tZero, k, deltaT)

    write (*,*) xK

contains
    real function f1(x, t) result(y)
        real, intent(in) :: x, t
        y = -2*t*(x**2)
    end function

end program lista_6