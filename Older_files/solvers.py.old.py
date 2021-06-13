import numpy as np
from scipy import integrate
from scipy.integrate import solve_bvp
import math



class ODESolver(object):
    """The ODE is solved here."""
    def __init__(self, ode_config):
        self.config = ode_config
        self.Upper_dy_lim = ode_config.Upper_dy_lim
        self.Lower_dy_lim = ode_config.Lower_dy_lim
        self.it = ode_config.num_iterations
        self.eps = ode_config.error
        self.N = ode_config.N
        self.dx = None
        self.df0 = None
        self.func = None
        self.message = None
        self.__ddw__ = None
        self.__dw__ = None
        self.__w__ = None
        self.y = np.empty(self.N)
        self.dy = np.empty(self.N)
        self.ddy = np.empty(self.N)
        self.x_p = None
        self.x_s = None
        self.x = None


    def test(self):
        """Nonlinearity H: ddy = H(dy,y,x)"""
        raise NotImplementedError     # def solve(self):
    #     """Nonlinearity H: ddy = H(dy,y,x)"""
    #     raise NotImplementedError

    def solvr(Y, t, self):
        """function to solve the ode via integrate"""
        raise NotImplementedError

    #can we use argmin and min inside this funciton?
    def findMin(X):
        x=np.amin(X)
        for index, i in enumerate(X):
            if i == x:
                return index, i

class solve_F_w(ODESolver):
    def __init__(self, ode_config):
        super(solve_F_w, self).__init__(ode_config)
        print(r'F^{\prime}[0]: ', self.config.dy_init)

    def test(self):
        print(self.N)

    def solvr(Y, t, param):
        dotY = [Y[1], eqn.nonlinear(param[0] - param[2] * Y[0] + param[1]*t*Y[1], param[3], param[4], param[5])]
        return dotY


    def solve(self,eqn):
        self.dx = eqn.range/self.N
        self.x = np.linspace(0,eqn.range,self.N)
        print(self.dx)
        df_0 = self.df0
        #end point in bisection method
        Udf_tmp=self.Upper_dy_lim
        Ldf_tmp=self.Lower_dy_lim
        param = eqn.param
        eqn.Gamma_order()
        a_t = self.x
        b_t = (param[0] - param[1]*a_t)/param[2]
        #solving the ode(second order is change to a system of two first orders inside solvr)
        cond = [0, df_0]
        solvr = self.solvr
        asol = integrate.odeint(solvr, cond, a_t, args=(param,))
        dF=asol[:,1]
        F=asol[:,0]
        """[mu, gamma, r ,lmb ,sgm, rho]:"""
        ddF = self.eqn.nonlin(param[0] - param[2] * F + param[1]*np.multiply(a_t,dF), param[4], param[5])
        #finding point w_p using three methods
        __w__=findMin(dF)
        #error of the three evaluations of w_p
        var=np.absolute(F[__w__[0]]-b_t[__w__[0]])+np.absolute(dF[__w__[0]]+1) + np.absolute(ddF[__w__[0]])
        #iteration to minimize the error (brute force)
        k=0
        while k<self.it and var>self.eps:
            asol = integrate.odeint(self.solvr, [0, df_0], a_t, args=(param,))
            F=asol[:,0]
            dF=asol[:,1]
            """[mu, gamma, r ,lmb ,sgm, rho]:"""
            ddF = Nonlin(param[0] - param[2] * F + param[1]*np.multiply(a_t,dF), param[4], param[5])
            __w__=self.findMin(dF)
            if __w__[1]<-1:
                Ldf_tmp=df_0
                df_0 += (Udf_tmp-df_0)/2
            else:
                Udf_tmp=df_0
                df_0 -= df_0/2
                k += 1
                var=np.absolute(dF[__w__[0]]+1)+np.absolute(ddF[__w__[0]])+np.absolute(F[__w__[0]]-b_t[__w__[0]])
        print('Inside the method dF/dw(0) = ', dF[0])
        setattr(self,__ddw__,self.findMin(np.absolute(ddF)))
        self.__ddw__ = self.findMin(np.absolute(ddF))
        self.__dw__ = self.findMin(np.absolute(F-b_t))
        self.__w__ = self.findMin(dF)
        self.message = 'If the three numbers below are equal, the scheme is working.'\
        +'Reset dF(0) and rerun the code.'\
        +'Calculation of payment boundary by three methods:'\
        +'Minimum of abs(F-(mu-gamma*x)/r:'+str(__w__[0])+'\n'\
        +'Minimum of abs(dF):'+str(__dw__[0])+'\n'\
        +'Minimum of abs(ddF):'+str(__ddw__[0])
        self.y = F
        self.dy = dF
        self.ddy = ddF
        return F

    # What does function do? Find the value and index of the minimum of a float vector.
