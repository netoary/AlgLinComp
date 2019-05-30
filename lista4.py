import numpy as np

def f1(X):
    return [16*(X[0]**4) + 16*(X[1]**4) + (X[2]**4) - 16,
             X[0]**2 + X[1]**2 + X[2]**2 - 3,
             X[0]**3 - X[1] + X[2] - 1]

def jac1(X):
    return [[64*(X[0]**3),  64*(X[1]**3),   4*(X[2]**3)],
           [2*X[0],         2*X[1],         2*X[2]],
           [3*X[0]**2,      -1,             1]]

def newton(funcs, jacobian, X, tol, nIter):
    k = 0
    tolK = 100000000
    success = False
    while (k < nIter):
        k = k+1
        F = funcs(X)
        J = np.matrix(jacobian(X))
        J_1 = np.array(J.I)
        deltaX = - J_1 @ F#-np.matmul(J, F)
        X = X + deltaX

        tolK = np.linalg.norm(deltaX) / np.linalg.norm(X**2)
        
        if (tolK <= tol):
            success = True
            break
    return (success, X, k, tolK)


def broyden(funcs, jacobian, X, tol, nIter):
    k = 0
    tolK = 100000000
    J = np.matrix(jacobian(X))
    F_1 = funcs(X)
    while (k < nIter):
        k = k+1
        J_1 = np.array(J.I)
        
        deltaX = - J_1 @ F_1
       
        X = X + deltaX
        F = funcs(X)
        Y = np.array(F) - np.array(F_1)
        tolK = np.linalg.norm(deltaX) / np.linalg.norm(X**2)
        
        if (tolK <= tol):
            success = True
            break
        else:
            J = J + ((Y - J @ deltaX) @ deltaX.T) / (deltaX.T @ deltaX)
            F_1 = F
            
    return (success, X, k, tolK)
    
X = [1, 1, 1]
success, X, k, tolK = newton(f1, jac1, X, 0.000001, 100)
print ("Newton")
print ("Convergencia atingida: " + str(success))
print ("X = " + str(X))
print ("Total de iterações: " + str(k))
print ("Tolerância final: " + str(tolK))

print ("")

X = [1, 1, 1]
success, X, k, tolK = broyden(f1, jac1, X, 0.000001, 100)
print ("Broyden")
print ("Convergencia atingida: " + str(success))
print ("X = " + str(X))
print ("Total de iterações: " + str(k))
print ("Tolerância final: " + str(tolK))
