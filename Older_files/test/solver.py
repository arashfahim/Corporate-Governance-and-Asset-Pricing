import numpy as np
from scipy import integrate
import math


class ODE(object):
    """The ODE is solved here."""
    def __init__(self):
        self.cond = None
        self.param = None
        self.t = None
        self.x = None



    def nonlinear(X, param):
        raise NotImplementedError

    def solve(self):
        raise NotImplementedError

class my_ode(ODE):
    def __init__(self):
        super(my_ode, self).__init__()
        self.cond = [1, 0]
        self.param = [0, 0.5, 1, 0, [1.1], [0.1]] #matrix [a_11, a_12, a_21, a_22]
        self.t = np.linspace(0,1,100)
        self.x = np.zeros((2,100))

    def nonlinear(self,X):
        """Nonlinearity H: ddF = H(dF,F,w)"""
        sgm = self.param[4]
        rho = self.param[5]
        if len(sgm) != len(rho):
            print("sigma and rho must me the same length.")
        H_inv = -2*(X-rho[0])/(sgm[0]*sgm[0])
        for i in range(1,len(sgm)):
            H_inv = np.min(H_inv, -2*(X-rho[i])/(sgm[i]*sgm[i]), axis=0)
        return H_inv

    def solve(self):
        def RHS(X, t, param):
            F = [X[1], self.nonlinear(param[2] * X[0] + param[3] * X[1])]
            return F
        param = self.param
        print(param)
        cond = self.cond
        t = self.t
        return integrate.odeint(RHS, cond, t, args=(param,))
