include 'methods.f90'
include 'functions.f90'


program lista_6
    real :: xK, xZero, tZero, deltaT, xLZero
    integer :: k

    xK = 0.0
    xZero = 0.0
    tZero = 0.0
    deltaT = 0.1
    k = 10
    xLZero = 0.0

    !call euler(xK, xZero, tZero, k, deltaT)
    !call rungeKutta2(xK, xZero, tZero, k, deltaT)
    !call rungeKutta4(xK, xZero, tZero, k, deltaT)

    !TESTAR ESSAS DUAS
    !call taylor2(xK, xZero, xLZero, tZero, k, deltaT)
    !call rungeKuttaNystrom(xK, xZero, xLZero, tZero, k, deltaT)

    write (*,*) xK

end program lista_6