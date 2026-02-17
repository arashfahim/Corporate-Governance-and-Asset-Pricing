import numpy as np
import math
from scipy import integrate
from scipy.integrate import solve_bvp
import datetime as dt
import os

def findMin(X):
    x=np.amin(X)
    for index, i in enumerate(X):
        if i == x:
            return index, i


class Equation(object):
    """Base class for a boundary value problem"""

    def __init__(self,config):
        self.config = config
        self.range = config.Numerical_settings.range
        """param=[mu, gamma, r ,lmb ,sgm, rho]"""
        self.param = list(config.Parameters.values())
        self.y_init = 0
        self.Gamma_order()
        """Parameters for the solver"""
        self.N = config.Numerical_settings.Number_of_points
        self.dx = self.range/self.N
        self.x = np.linspace(0,self.range,self.N)
        self.y = np.empty(self.N)
        self.line = [(self.param[0] - self.param[1] * x)/self.param[2] for x in self.x]
        self.dy = np.empty(self.N)
        self.ddy = np.empty(self.N)
        self.Upper_dy_lim = config.Numerical_settings.Upper_dy_lim
        self.Lower_dy_lim = config.Numerical_settings.Lower_dy_lim
        self.it = config.Numerical_settings.num_iterations
        self.eps = config.Numerical_settings.Stop_criterion_for_F_ODE
        self.df0 = config.Numerical_settings.df0
        self.path = os.path.dirname(__file__) #path to save the output file
        # self.y = np.empty(self.N)
        # self.dy = np.empty(self.N)
        # self.ddy = np.empty(self.N)
        # self.num_paths =  config.Simulation_settings.num_paths
        # self.max_length = config.Simulation_settings.max_length
        # self.num_samples =  self.num_paths*self.max_length
        self.solver_tag = True




    def nonlinear(self,X,param):
        """Nonlinearity H: ddy = H(dy,y,x)"""
        raise NotImplementedError

    def optimal(self,t):
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
        lmbd_ = self.param[3]
        self.param[0] = self.param[0] - rho_[0]
        rho_ = [r - rho_[0] for r in rho_]
        self.Gamma = []
        self.swtch = []
        _sgm_ = [sgm_[0]]
        _rho_ = [rho_[0]]
        for i in range(0,len(sgm_)):
            for j in range(1,len(sgm_)):
                if i < j:
                    gamma = 2*(rho_[i]-rho_[j])/(lmbd_*lmbd_*(sgm_[i]*sgm_[i]-sgm_[j]*sgm_[j]))
                    val = -gamma*(lmbd_*lmbd_*sgm_[j]*sgm_[j])/2 + rho_[j]
                    val_tmp = val
                    for k in range(0,len(sgm_)):
                        if k != i and k != j:
                            val_tmp =  np.minimum(val_tmp,-gamma*(lmbd_*lmbd_*sgm_[k]*sgm_[k])/2 + rho_[k])
                    if val_tmp == val:
                        gamma_ = [gamma , i , j]
                        self.Gamma.append(gamma_)
        self.Gamma = sorted(self.Gamma, reverse =True)
        # print('sgm = ',sgm_[0], '      AND rho = ', rho_[0])
        print('Switching values and active regimes:\n ')
        for index, i in enumerate(self.Gamma):
                _sgm_.append(sgm_[self.Gamma[index][2]])
                _rho_.append(rho_[self.Gamma[index][2]])
                self.swtch.append(self.Gamma[index][0])#holds the values of ddF where switching happends
                print(self.Gamma[index][1], ' to ', self.Gamma[index][2], ' at Gamma_',index+1,' = ', round(self.Gamma[index][0],3),' \n')
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


    def optimal(self,t):#def optimal_sigma_rho(self,t):
        sgm_ = self.param[4]
        rho_ = self.param[5]
        y = np.ones(len(t))
        z = np.ones(len(t))
        if len(sgm_) == 1:
            y = y * sgm_[0]
            z = z * rho_[0]
        else:
            temp_swtch = [0]
            temp_swtch.extend(self.swtch)
            temp_swtch.append(self.swtch[-1]+np.min(self.ddf)-1) # the switch cannot happen at values of Gamma than the minimum of ddf
            for indexm, i in enumerate(t):
                #solver of BVP: scipy.integrate.solve_bvp requires
                #this to work for t of any size. This is the reason
                #behind the loop.
                i__ = self.findD(np.abs(i-self.m))#closest point on the m-grid to t
                y[indexm] = self.ddf[i__[0]]
                z[indexm] = self.ddf[i__[0]]
            for index, g in enumerate(temp_swtch[0:-1]):
                y = np.where((y <= g) & (y > temp_swtch[index+1]), sgm_[index], y)
                z = np.where((z <= g) & (z > temp_swtch[index+1]), rho_[index], z)
        return y,z
    

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
        self.__dw__ = findMin(np.absolute(self.dy+1))



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
        self.x_s = self.x[__x_s__[0]]#switching points
        ind = self.N #np.min([int((self.__w__[0]+1)*1.2),self.N])
        self.x = self.x[0:ind]
        # print(self.__w__[0]+1)
        self.x_p = self.x[self.__w__[0]+1]#payment point
        self.m = self.x/self.param[3]# cash reserve variable
        self.mS = self.m[0:(self.__w__[0]+1)]# cash reserve variable

        self.y = np.where(self.ddy[0:ind] <= 0, self.y[0:ind], self.y[self.__w__[0]]-(self.x - self.x_p))
        self.dy = np.where(self.ddy[0:ind] <= 0, self.dy[0:ind], -1)
        self.ddy = np.where(self.ddy[0:ind] <= 0, self.ddy[0:ind], 0)
        self.line = self.line[0:ind]


        self.m_s = self.x_s/self.param[3]#switching points in cash reserve variable m
        self.m_p = self.x_p/self.param[3]#payment point in cash reserve variable m
        self.f = self.y
        self.df = self.dy*self.param[3]
        self.ddf = self.ddy*self.param[3]*self.param[3]
        self.fline = [(self.param[0] - self.param[1] * self.param[3] * x)/self.param[2] for x in self.m]
        data = np.column_stack((self.x, self.y, self.dy, self.ddy, self.line, self.m, self.f, self.df, self.ddf, self.fline))
        header_text = 'x, F, dF, ddF, Fline, m, f, df, ddf, fline'
        np.savetxt(self.path+'output.dat',  # Filename
                data,           # Array data
                fmt='%.4f',     # Format of data values (4 decimal places)
                delimiter=',',  # String or character separating columns
                header=header_text, # Text written at the beginning of the file
                comments='# '   # String prefixed to header to mark as a comment
        )        


        self.sol= {'domain':self.x, 'solution':self.y,'1st_deriv':self.dy,'2nd_deriv':self.dy,'3rd_deriv':self.ddy,'pay_point':self.x_p, 'negot_point':self.x_s}
        opt_sig, opt_rho = self.optimal(self.m)
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

        self.solver_tag = False
        ####################



class eqn_FfS(eqn_Ff):
    """The equation for principal as a function of the agents value"""
    def __init__(self,config):
        super().__init__(config)



    def bvp_S(self):
        def eqn_form(t, X):# defining equation for S
            opt_sig, opt_rho = self.optimal(t)
            dis_cnt_coeff = 2/(opt_sig*opt_sig)
            H_ = np.vstack((X[1], dis_cnt_coeff*(self.param[2]*X[0] -self.param[1]*t * X[1])))
            return H_
        def bc(ya, yb):# defining boundary values for S
            return np.array([ya[0], yb[1]-1])
        
        init_y  = np.ones((2,len(self.mS)))
        init_y[0] = np.power(self.mS,4/5)
        a = dt.datetime.now()
        self.solution_S = solve_bvp(eqn_form, bc, self.mS, init_y)
        b = dt.datetime.now() - a
        print('BVP for S is solved in',b.total_seconds(),'seconds.\n')
        self.S = self.solution_S.sol(self.mS)[0]
        self.dS = self.solution_S.sol(self.mS)[1]
        opt_sig, _ = self.optimal(self.mS)
        dis_cnt_coeff = 2/(opt_sig*opt_sig)
        self.ddS = (self.param[2] * self.S - self.param[1] * self.mS * self.dS)\
*dis_cnt_coeff
        for i in self.swtch_pos:
            j = self.findD(np.abs(self.mS - self.mS[i]))
            self.ddS[j[0]] = np.nan
            self.ddS[j[0]+1] = np.nan



    def bvp_T(self):
        def eqn_form(t, X):# defining equation for T
            opt_sig, _ = self.optimal(t)
            dis_cnt_coeff = 2/(opt_sig*opt_sig)
            H_ = np.vstack((X[1], dis_cnt_coeff*(-1 + self.param[2]*X[0] -self.param[1]*t * X[1])))
            return H_
        def bc(ya, yb):# defining boundary values for T
            return np.array([ya[0], yb[1]])

        init_y  = np.ones((2,len(self.mS)))
        init_y[0] = np.power(self.mS,4/5)
        a = dt.datetime.now()
        self.solution_T = solve_bvp(eqn_form, bc, self.mS, init_y)
        b = dt.datetime.now() - a
        print('BVP for T is solved in',b.total_seconds(),'seconds.\n')
        self.T = self.solution_T.sol(self.mS)[0]
        self.dT = self.solution_T.sol(self.mS)[1]
        opt_sig , _ = self.optimal(self.mS)
        dis_cnt_coeff = 2/(opt_sig*opt_sig)
        self.ddT = (-1 + self.param[2] * self.T - self.param[1] * self.mS * self.dT)\
*dis_cnt_coeff
        for i in self.swtch_pos:#setting value at place of jumps avoids ploting jumps as a continuous function
            j = self.findD(np.abs(self.mS - self.mS[i]))
            self.ddT[j[0]] = np.nan
            self.ddT[j[0]+1] = np.nan

    def bvp_C(self):
        def eqn_form(t, X):# defining equation for T
            opt_sig, opt_rho = self.optimal(t)
            dis_cnt_coeff = 2/(opt_sig*opt_sig)
            H_ = np.vstack((X[1], dis_cnt_coeff*(-opt_rho + self.param[2]*X[0] -self.param[1]*t * X[1])))
            return H_
        def bc(ya, yb):# defining boundary values for T
            return np.array([ya[0], yb[1]])

        init_y  = np.ones((2,len(self.mS)))
        init_y[0] = np.power(self.mS,4/5)
        a = dt.datetime.now()
        self.solution_C = solve_bvp(eqn_form, bc, self.mS, init_y)
        b = dt.datetime.now() - a
        print('BVP for C is solved in',b.total_seconds(),'seconds.\n')
        self.C = self.solution_C.sol(self.mS)[0]
        self.dC = self.solution_C.sol(self.mS)[1]
        opt_sig, opt_rho = self.optimal(self.mS) 
        dis_cnt_coeff = 2/(opt_sig*opt_sig)
        self.ddC = (- opt_rho + self.param[2] * self.C - self.param[1] * self.mS * self.dC)\
*dis_cnt_coeff
        for i in self.swtch_pos:#setting value at place of jumps avoids plotting jumps as a continuous function
            j = self.findD(np.abs(self.mS - self.mS[i]))
            self.ddC[j[0]] = np.nan
            self.ddC[j[0]+1] = np.nan

    def bvp(self):# solves both equations for T and S
        if self.solver_tag:#if bvp solved tag is false
            self.solver()
        self.bvp_S()
        self.bvp_T()
        self.bvp_C()
        # mf=self.f[0:(self.__w__[0]+1)]#f(m) in cash reserve range: until the payment point only
        self.C1=self.param[0]*self.T-self.param[3]*self.S-self.f[0:(self.__w__[0]+1)]#monitoring cost obtained from the equality f(m)=\mu/r T(m) - \lambda*S(m)-monit_cost(m)
        self.dC1=self.param[0]*self.dT-self.param[3]*self.dS-self.df[0:(self.__w__[0]+1)]
        self.ddC1=self.param[0]*self.ddT-self.param[3]*self.ddS-self.ddf[0:(self.__w__[0]+1)]
        data = np.column_stack((self.S, self.T, self.C))
        header_text = 'S, T, C'
        np.savetxt(self.path+'/output_STC.dat',  # Filename
                data,           # Array data
                fmt='%.4f',     # Format of data values (4 decimal places)
                delimiter=',',  # String or character separating columns
                header=header_text, # Text written at the beginning of the file
                comments='# '   # String prefixed to header to mark as a comment
        )        




#END OF THE CODE
