include 'methods.f90'
include 'functions.f90'


program lista_7
    real :: x, deltaX, deltaX1, deltaX2, fLinha, p

    deltaX = 0.001
    deltaX1 = 0.5
    deltaX2 = 0.25
    x = 2
    p = 1
    fLinha = 0

    !call diferencaCentral(x, deltaX, fLinha)
    !call passoFrente(x, deltaX, fLinha)
    !call passoTras(x, deltaX, fLinha)
    call interpolacaoRichard(x, deltaX1, deltaX2, fLinha, p)
    !call f(b, y)

    write (*,*) fLinha

end program lista_7