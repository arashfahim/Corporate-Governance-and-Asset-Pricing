import numpy as np
import math
from scipy import integrate
from scipy.integrate import solve_bvp
import datetime as dt


def findMin(X):
    x=np.amin(X)
    for index, i in enumerate(X):
        if i == x:
            return index, i


class Equation(object):
    """Base class for a boundary value problem"""

    def __init__(self,config):
        self.config = config
        self.range = config.eqn_config.range
        """param=[mu, gamma, r ,lmb ,sgm, rho]"""
        self.param = config.eqn_config.param
        # self.param_alt = config.eqn_config.param_alt
        self.y_init = 0
        self.Gamma = self.Gamma_order()
        """Parameters for the solver"""
        self.N = self.config.ode_config.N
        self.dx = self.range/self.N
        self.x = np.linspace(0,self.range,self.N)
        self.y = np.empty(self.N)
        self.line = [(self.param[0] - self.param[1] * x)/self.param[2] for x in self.x]
        self.dy = np.empty(self.N)
        self.ddy = np.empty(self.N)
        self.Upper_dy_lim = config.ode_config.Upper_dy_lim
        self.Lower_dy_lim = config.ode_config.Lower_dy_lim
        self.it = config.ode_config.num_iterations
        self.eps = config.ode_config.eps
        self.N = config.ode_config.N
        self.df0 = config.ode_config.df0
        self.y = np.empty(self.N)
        self.dy = np.empty(self.N)
        self.ddy = np.empty(self.N)
        self.num_paths =  config.bvp_config.num_paths
        self.max_length = config.bvp_config.max_length
        self.num_samples =  self.num_paths*self.max_length




    def nonlinear(self,X,param):
        """Nonlinearity H: ddy = H(dy,y,x)"""
        raise NotImplementedError

    def optimal_sigma(self,X,param):
        """Nonlinearity H: ddy = H(dy,y,x)"""
        raise NotImplementedError

    def Gamma_order(self):
        """Filtering sigma and rho for """
        raise NotImplementedError

    def solver(self):
        """Solver for the equation"""
        raise NotImplementedError

    def findD(self,X):
        x=np.amin(X)
        for index, i in enumerate(X):
            if i == x:
                return index, i


class eqn_Ff(Equation):
    """The equation for principal as a function of the agents value"""
    def __init__(self,config):
        super().__init__(config)


    def Gamma_order(self):
        sgm_ = self.param[4]
        rho_ = self.param[5]
        self.param[0] = self.param[0] - rho_[0]
        rho_ = [r - rho_[0] for r in rho_]
        Gamma = []
        self.swtch = []
        _sgm_ = [sgm_[0]]
        _rho_ = [rho_[0]]
        for i in range(0,len(sgm_)):
            for j in range(1,len(sgm_)):
                if i < j:
                    gamma = 2*(rho_[i]-rho_[j])/(sgm_[i]*sgm_[i]-sgm_[j]*sgm_[j])
                    val = -gamma*(sgm_[j]*sgm_[j])/2 + rho_[j]
                    val_tmp = val
                    for k in range(0,len(sgm_)):
                        if k != i and k != j:
                            val_tmp =  np.minimum(val_tmp,-gamma*(sgm_[k]*sgm_[k])/2 + rho_[k])
                    if val_tmp == val:
                        gamma_ = [gamma , i , j]
                        Gamma.append(gamma_)
        Gamma = sorted(Gamma, reverse =True)
        # print('sgm = ',sgm_[0], '      AND rho = ', rho_[0])
        print('Switching values and active regimes:\n ')
        for index, i in enumerate(Gamma):
                _sgm_.append(sgm_[Gamma[index][2]])
                _rho_.append(rho_[Gamma[index][2]])
                self.swtch.append(Gamma[index][0])#/(self.param[3]*self.param[3]))
                print(Gamma[index][1], ' to ', Gamma[index][2], ' at Gamma_',index+1,' = ', round(Gamma[index][0],3),' \n')
        # _sgm_.reverse()
        self.param[4] = _sgm_
        self.param[5] = _rho_


    def nonlinear(self,X):
        """Nonlinearity H: ddF = H(dF,F,w)"""
        lmb = self.param[3]
        sgm = self.param[4]
        rho = self.param[5]
        if len(sgm) != len(rho):
            raise print("sigma and rho must me the same length.")
        H_inv = -2*(X-rho[0])/(sgm[0]*sgm[0]*lmb*lmb)
        for i in range(1,len(sgm)):
            H_inv = np.min((H_inv, -2*(X-rho[i])/(sgm[i]*sgm[i]*lmb*lmb)), axis=0)
        return H_inv


    def optimal_sigma(self,t):
        sgm_ = self.param[4]
        y = np.ones(len(t))
        if len(sgm_) == 1:
            y = y * sgm_[0]
        else:
            temp_swtch = [0]
            temp_swtch.extend(self.swtch)
            temp_swtch.append(self.swtch[-1]+np.min(self.ddf)-1)
            for indexm, i in enumerate(t):
                #solver of BVP: scipy.integrate.solve_bvp requires
                #this to work for t of any size. This is the reason
                #behind the loop.
                i__ = self.findD(np.abs(t[indexm]-self.m))
                y[indexm] = self.ddf[i__[0]]
            for index, g in enumerate(temp_swtch[0:-1]):
                y = np.where((y <= g) & (y > temp_swtch[index+1]), sgm_[index], y)
        return y



    def bvp_solve(self):
        def RHS(X, t, param):
            H_ = [X[1], self.nonlinear(param[0] - param[2] * X[0] + param[1]*np.multiply(t,X[1]))]
            return H_
        df0 = self.df0
        cond = [0, df0]
        Udf_tmp=self.Upper_dy_lim
        Ldf_tmp=self.Lower_dy_lim
        param = self.param
        # self.Gamma_order()
        a_t = self.x
        b_t = self.line
        asol = integrate.odeint(RHS, cond, a_t, args=(param,))
        F=asol[:,0]
        dF=asol[:,1]
        """[mu, gamma, r ,lmb ,sgm, rho]:"""
        ddF = self.nonlinear(param[0] - param[2] * F + param[1]*np.multiply(a_t,dF))
        __w__=findMin(dF)
        if __w__[1]<-1:
            Ldf_tmp=df0
            df0 += (Udf_tmp-df0)/2
        else:
            Udf_tmp=df0
            df0 -= df0/2
        self.df0 = df0
        self.Upper_dy_lim = Udf_tmp
        self.Lower_dy_lim = Ldf_tmp
        self.y = F
        self.dy = dF
        self.ddy = ddF
        self.__ddw__ = findMin(np.absolute(self.ddy))
        self.__w__ = findMin(np.absolute(self.y-self.line))
        self.__dw__ = findMin(self.dy)



    def solver(self):
        var = 1
        k=0
        while k<self.it and var>self.eps:
            self.bvp_solve()
            k += 1
            temp_k = k/50
            var=np.absolute(self.dy[self.__w__[0]]+1)*np.absolute(self.dy[self.__w__[0]]+1)\
+np.absolute(self.ddy[self.__w__[0]])*np.absolute(self.ddy[self.__w__[0]])\
+np.absolute(self.y[self.__w__[0]]-self.line[self.__w__[0]])*\
np.absolute(self.y[self.__w__[0]]-self.line[self.__w__[0]])
            if (temp_k.is_integer()):
                print('var = ',round(var,5),'       Iteration = ',k)
        print('var = ',round(var,5),'       Iteration = ',k)




        __x_s__ = findMin(np.absolute(self.dy))
        self.x_s = self.x[__x_s__[0]]
        print(self.__w__[0]+1)
        self.x_p = self.x[self.__w__[0]+1]
        self.m = self.x/self.param[3]

        self.y = np.where(self.ddy <= 0, self.y, self.y[self.__w__[0]]-(self.x - self.x_p))
        self.dy = np.where(self.ddy <= 0, self.dy, -1)
        self.ddy = np.where(self.ddy <= 0, self.ddy, 0)



        self.m_s = self.x_s/self.param[3]
        self.m_p = self.x_p/self.param[3]
        self.f = self.y
        self.df = self.dy*self.param[3]
        self.ddf = self.ddy*self.param[3]*self.param[3]
        self.fline = [(self.param[0] - self.param[1] * self.param[3] * x)/self.param[2] for x in self.m]
        self.sol= {'domain':self.x, 'solution':self.y,'1st_deriv':self.dy,'2nd_deriv':self.dy,'3rd_deriv':self.ddy,'pay_point':self.x_p, 'negot_point':self.x_s}
        opt_sig = self.optimal_sigma(self.m)
        self.swtch_pos = np.where(np.abs(np.diff(opt_sig)) > 0)[0]
        swtch = [str(round(self.m[I],2)) for I in self.swtch_pos]
        swtch_ = 'The switching point(s) are given by '+ ', '.join(swtch)
        self.message = f'''\nThis solves the value function as a function of w, the agent's utility.\n\n'''\
'''If the three numbers below, approximations of payment boundary by three methods,\
are almost equal, the scheme is working. \
\nOtherwise, reset df0 in the json file and run the code again.\n\
Calculation of payment boundary by three methods: \n\
Minimum of abs(F-(mu-gamma*x)/r,  Minimum of abs(dF+1), and Minimum of abs(ddF) are, respectively, at\n'''\
+str(round(self.x[self.__w__[0]],5))+', '\
+str(round(self.x[self.__dw__[0]],5))+', and '\
+str(round(self.x[self.__ddw__[0]],5))+'. '\
+'''\nThe value of dF(0) is '''+ str(round(self.dy[0],5))\
+'''. \nIf the principal has the bargaining power, the minimum utility of the agent is '''\
+str(round(self.x_s,5))+'.'\
+'''\nThe initial capital to start implemeting the contract is '''\
+str(round(self.x_s/self.param[3],5))+'. '+''' \nThe divident boundary is given by '''\
+str(round(self.m_p,5))+'. ' + '\n'\
+swtch_+'.'

        ####################



class eqn_FfS(eqn_Ff):
    """The equation for principal as a function of the agents value"""
    def __init__(self,config):
        super().__init__(config)



    def bvp_S(self):
        self.solver()
        def eqn_form(t, X):
            dis_cnt_coeff = 2/(self.optimal_sigma(t)*self.optimal_sigma(t))
            H_ = np.vstack((X[1], dis_cnt_coeff*(self.param[2]*X[0] -self.param[1]*t * X[1])))
            return H_
        def bc(ya, yb):
            return np.array([ya[0], yb[1]-1])
        self.mS = self.x[0:(self.__w__[0]+1)]
        init_y  = np.ones((2,len(self.mS)))
        init_y[0] = np.power(self.mS,4/5)
        a = dt.datetime.now()
        solution = solve_bvp(eqn_form, bc, self.mS, init_y)
        b = dt.datetime.now() - a
        print('BVP for S is solved in',b.total_seconds(),'seconds.\n')
        self.S = solution.sol(self.mS)[0]
        self.dS = solution.sol(self.mS)[1]
        dis_cnt_coeff = 2/(self.optimal_sigma(self.mS)*self.optimal_sigma(self.mS))
        self.ddS = (self.param[2] * self.S - self.param[1] * self.mS * self.dS)\
*dis_cnt_coeff
        for i in self.swtch_pos:
            j = self.findD(np.abs(self.mS - self.m[i]))
            self.ddS[j[0]] = np.nan
            self.ddS[j[0]+1] = np.nan


    # def sample(self,x):
    #     self.Dlt_B = np.random.normal(0, self.dx, size=(self.num_samples))
    #     dt = self.dx*self.dx
    #     i = 0
    #     self.paths = []
    #     for j in range(0,self.num_paths):
    #         path_ = [np.random.uniform(0,x_p,1)[0]]
    #         while ((i < self.max_length) & (path_[-1] > 0)):
    #             tmp = np.amin([x_p, path_[-1] + Dlt_B[i]])
    #             path_.append(tmp)
    #             i = i + 1
    #             if (path_[-1] < 0):
    #                 path_[-1] = 0
    #                 self.paths.append(path_)
    #     self.input = np.array([p[0] for p in self.paths.reshape(-1, 1)

    # def base(x,y,x_input,K):
    #     I = []
    #     y_output = []
    #     a = np.empty(K)
    #     b = np.empty(K)
    #     X_ind = np.empty(K)
    #     X = np.empty(K)
    #     X_sq = np.empty(K)
    #     Y = np.empty(K)
    #     YX = np.empty(K)
    #     x_y = np.stack((x,y.reshape(-1,1)),axis=1)
    #     for x0 in x_input:
    #         for i in range(1,K):
    #             if ((x0 > (i-1)*x_p/K) & (x0 <= i*x_p/K)):
    #                 I = i
    #         X_ind = np.sum([1 if ((z[0] > (I-1)*x_p/K) & (z[0] <= I*x_p/K)) else 0 for z in x_y])
    #         X = np.sum([z[0] if ((z[0] > (I-1)*x_p/K) & (z[0] <= I*x_p/K)) else 0 for z in x_y])
    #         X_sq = np.sum([z[0]*z[0] if ((z[0] > (I-1)*x_p/K) & (z[0] <= I*x_p/K)) else 0 for z in x_y])
    #         Y = np.sum([z[1] if ((z[0] > (I-1)*x_p/K) & (z[0] <= I*x_p/K)) else 0 for z in x_y])
    #         YX = np.sum([z[1]*z[0] if ((z[0] > (I-1)*x_p/K) & (z[0] <= I*x_p/K)) else 0 for z in x_y])
    #         det = X_sq*X_ind-X*X
    #         a = (X_ind*YX-X*Y)/det
    #         b = (-X*YX+X_sq*Y)/det
    #         y_output.append(np.multiply(a,x0) + b)
    #     return y_output
    #
    #










#END OF THE CODE
