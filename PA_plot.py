import matplotlib.pyplot as plt
plt.rcParams.update({
    "text.usetex": True,
    # "font.family": "sans-serif",
    # "font.sans-serif": ["Helvetica"],
    'text.latex.preamble' : r'\usepackage{amssymb} \usepackage{amsmath}' #for \text command
    })
import matplotlib.mlab as mlab
import matplotlib.gridspec as gridspec
import tikzplotlib
import os

class myfigures(object):
    def __init__(self,PA,string):
        '''obj is a list of solved equations with different parameter sets'''
        self.path = os.path.dirname(__file__)
        self.string = string
        self.lab = []
        self.obj = PA
        self.ind = len(PA)
        self.fig = plt.figure(figsize=(10, 5),constrained_layout=True)
        self.spec = gridspec.GridSpec(ncols=3, nrows=1, figure=self.fig)
        for indx, obj in enumerate(self.obj):
            # for r in obj.param[5]:
            #     rho += [str(r)]
            # for s in obj.param[4]:
            #     sgm += [str(s)]
            # rho = ','.join(rho)
            # sgm = ','.join(sgm)
            rho = str(obj.param[5]).rstrip("]").lstrip("[")
            sgm = str(obj.param[4]).rstrip("]").lstrip("[")
            mu = str(round(obj.param[0],4))
            gamma = str(round(obj.param[1],4))
            R = str(round(obj.param[2],4))
            lmb = str(round(obj.param[3],4))
            lab = 'Parameters: ' + r'$\mu = $' + mu + '\n' +\
r'$\gamma = $' + gamma + ' $r = $' + R + '\n' + r'$ \lambda = $' + lmb + '\n' +\
r'$\sigma = ($' + sgm + '$)$' + '\n' + r' ${\rho} = ($' + rho + '$)$' + '\n' +\
'$w_p = $' +  str(round(obj.x_p,4)) + '\n' +\
'$m_p = $' +  str(round(obj.m_p,4)) + '\n ' +  r' $F^{\prime}(0) = $' + \
str(round(obj.dy[0],4))
            self.lab.append(lab)        
        
        self.ax = self.fig.add_subplot(self.spec[0, 0])
        self.ax.set_title(self.string)
        self.ax_d = self.fig.add_subplot(self.spec[0, 1])
        self.ax_dd = self.fig.add_subplot(self.spec[0, 2])

    def plot_iT(self):
            if self.string == '$F(w)$':
                string1 = '$F^{\prime}(w)$'
                self.ax_d.set_title(string1)
                string2 = '$F^{\prime\prime}(w)$'
                self.ax_dd.set_title(string2)
                for indx, obj in enumerate(self.obj):
                    plt.gca().set_prop_cycle(None)
                    color = next(self.ax._get_lines.prop_cycler)['color']
                    self.ax.plot(self.obj[indx].x,self.obj[indx].y, label = self.lab[indx], color = color)
                    self.ax.plot(self.obj[indx].x,self.obj[indx].line, 'k:')
                    self.ax_d.plot(self.obj[indx].x,self.obj[indx].dy, label = self.lab[indx], color = color)
                    self.ax_dd.plot(self.obj[indx].x,self.obj[indx].ddy, label = self.lab[indx], color = color)
                self.ax.legend(loc='upper center', \
    shadow = True, ncol = 1, bbox_to_anchor=(1.5, 0.7))
                tikzplotlib.save(self.path+ "/valueFw.tex")

            if self.string == '$f(m)$':
                string1 = '$f^{\prime}(m)$'
                self.ax_d.set_title(string1)
                string2 = '$f^{\prime\prime}(m)$'
                self.ax_dd.set_title(string2)
                for indx, obj in enumerate(self.obj):
                    plt.gca().set_prop_cycle(None)
                    color = next(self.ax._get_lines.prop_cycler)['color']
                    self.ax.plot(self.obj[indx].m,self.obj[indx].f, label = self.lab[indx], color = color)
                    self.ax.plot(self.obj[indx].m,self.obj[indx].fline, 'k:')
                    self.ax_d.plot(self.obj[indx].m,self.obj[indx].df, label = self.lab[indx], color = color)
                    self.ax_dd.plot(self.obj[indx].m,self.obj[indx].ddf, label = self.lab[indx], color = color)
                self.ax.legend(loc='upper center', \
    shadow = True, ncol = 1, bbox_to_anchor=(1.5, 0.7))
                tikzplotlib.save(self.path+"/valuefm.tex")

            if self.string == '$S(m)$':
                for indx, obj in enumerate(self.obj):
                    plt.gca().set_prop_cycle(None)
                    color = next(self.ax._get_lines.prop_cycler)['color']
                    self.ax.plot(self.obj[indx].mS,self.obj[indx].S, label = self.lab[indx], color = color)
                    self.ax_d.plot(self.obj[indx].mS,self.obj[indx].dS, label = self.lab[indx], color = color)
                    self.ax_dd.plot(self.obj[indx].mS,self.obj[indx].ddS, label = self.lab[indx], color = color)
                self.ax.legend(loc='upper center', \
    shadow = True, ncol = 1, bbox_to_anchor=(1.5, 0.7))
                tikzplotlib.save(self.path+"/asset.tex")

            if self.string == '$T(m)$':
                for indx, obj in enumerate(self.obj):
                    plt.gca().set_prop_cycle(None)
                    color = next(self.ax._get_lines.prop_cycler)['color']
                    self.ax.plot(self.obj[indx].mS,self.obj[indx].T, label = self.lab[indx], color = color)
                    self.ax_d.plot(self.obj[indx].mS,self.obj[indx].dT, label = self.lab[indx], color = color)
                    self.ax_dd.plot(self.obj[indx].mS,self.obj[indx].ddT, label = self.lab[indx], color = color)
                self.ax.legend(loc='upper center', \
    shadow = True, ncol = 1, bbox_to_anchor=(1.5, 0.7))
                tikzplotlib.save(self.path+"/survival.tex")
            
            if self.string == '$C(m)$':                
                for indx, obj in enumerate(self.obj):
                    plt.gca().set_prop_cycle(None)
                    color = next(self.ax._get_lines.prop_cycler)['color']
                    self.ax.plot(self.obj[indx].mS,self.obj[indx].C, label = self.lab[indx], color = color)
                    self.ax_d.plot(self.obj[indx].mS,self.obj[indx].dC, label = self.lab[indx], color = color)
                    self.ax_dd.plot(self.obj[indx].mS,self.obj[indx].ddC, label = self.lab[indx], color = color)
                self.ax.legend(loc='upper center', \
    shadow = True, ncol = 1, bbox_to_anchor=(1.5, 0.7))
                tikzplotlib.save(self.path+"/cost.tex")

            if self.string == '$C(m)=\mu*T(m)-\lambda S(m)-f(m)$': 
                for indx, obj in enumerate(self.obj):
                    plt.gca().set_prop_cycle(None)
                    color = next(self.ax._get_lines.prop_cycler)['color']
                    self.ax.plot(self.obj[indx].mS,self.obj[indx].C1, label = self.lab[indx], color = color)
                    self.ax_d.plot(self.obj[indx].mS,self.obj[indx].dC1, label = self.lab[indx], color = color)
                    self.ax_dd.plot(self.obj[indx].mS,self.obj[indx].ddC1, label = self.lab[indx], color = color)
                self.ax.legend(loc='upper center', \
    shadow = True, ncol = 1, bbox_to_anchor=(1.5, 0.7))
                tikzplotlib.save(self.path+"/cost1.tex")            

            if self.string == '$C(m)+\lambda S(m)$': 
                for indx, obj in enumerate(self.obj):
                    plt.gca().set_prop_cycle(None)
                    color = next(self.ax._get_lines.prop_cycler)['color']
                    self.ax.plot(self.obj[indx].mS,self.obj[indx].C+self.obj[indx].param[3]*self.obj[indx].S, label = self.lab[indx], color = color)
                    self.ax_d.plot(self.obj[indx].mS,self.obj[indx].dC+self.obj[indx].param[3]*self.obj[indx].dS, label = self.lab[indx], color = color)
                    self.ax_dd.plot(self.obj[indx].mS,self.obj[indx].ddC+self.obj[indx].param[3]*self.obj[indx].ddS, label = self.lab[indx], color = color)
                self.ax.legend(loc='upper center', \
    shadow = True, ncol = 1, bbox_to_anchor=(1.5, 0.7))
                tikzplotlib.save(self.path+"/cost_lambda_asset.tex")     

            if self.string == '$\lambda S(m)$': 
                for indx, obj in enumerate(self.obj):
                    plt.gca().set_prop_cycle(None)
                    color = next(self.ax._get_lines.prop_cycler)['color']
                    self.ax.plot(self.obj[indx].mS,self.obj[indx].param[3]*self.obj[indx].S, label = self.lab[indx], color = color)
                    self.ax_d.plot(self.obj[indx].mS,self.obj[indx].param[3]*self.obj[indx].dS, label = self.lab[indx], color = color)
                    self.ax_dd.plot(self.obj[indx].mS,self.obj[indx].param[3]*self.obj[indx].ddS, label = self.lab[indx], color = color)
                self.ax.legend(loc='upper center', \
    shadow = True, ncol = 1, bbox_to_anchor=(1.5, 0.7))
                tikzplotlib.save(self.path+"/lambda_asset.tex")   

            plt.show(block=False)
            plt.pause(1)
    
    def close_all(self):
            input("PRESS ENTER TO CLOSE.")
            plt.close('all')    

#######################
# class plotF(object):
#     def __init__(self,PA):
#         self.path = os.path.dirname(__file__)
#         print(self.path)
#         self.obj = PA
#         self.ind = len(PA)
#         self.figF = plt.figure(figsize=(10, 5),constrained_layout=True)
#         self.specF = gridspec.GridSpec(ncols=3, nrows=1, figure=self.figF)
#         self.figf = plt.figure(figsize=(10, 5),constrained_layout=True)
#         self.specf = gridspec.GridSpec(ncols=3, nrows=1, figure=self.figf)        
#         self.figS = plt.figure(figsize=(10, 5),constrained_layout=True)
#         self.specS = gridspec.GridSpec(ncols=3, nrows=1, figure=self.figS)
#         self.figT = plt.figure(figsize=(10, 5),constrained_layout=True)
#         self.specT = gridspec.GridSpec(ncols=3, nrows=1, figure=self.figT)  
#         self.figC = plt.figure(figsize=(10, 5),constrained_layout=True)
#         self.specC = gridspec.GridSpec(ncols=3, nrows=1, figure=self.figC)  
#         self.figC1 = plt.figure(figsize=(10, 5),constrained_layout=True)
#         self.specC1 = gridspec.GridSpec(ncols=3, nrows=1, figure=self.figC1)  
#         self.lab = []
#         rho = []
#         sgm = []
#         for indx, obj in enumerate(self.obj):
#             # for r in obj.param[5]:
#             #     rho += [str(r)]
#             # for s in obj.param[4]:
#             #     sgm += [str(s)]
#             # rho = ','.join(rho)
#             # sgm = ','.join(sgm)
#             rho = str(obj.param[5]).rstrip("]").lstrip("[")
#             sgm = str(obj.param[4]).rstrip("]").lstrip("[")
#             mu = str(round(obj.param[0],4))
#             gamma = str(round(obj.param[1],4))
#             R = str(round(obj.param[2],4))
#             lmb = str(round(obj.param[3],4))
#             lab = 'Parameters: ' + r'$\mu = $' + mu + '\n' +\
# r'$\gamma = $' + gamma + ' $r = $' + R + '\n' + r'$ \lambda = $' + lmb + '\n' +\
# r'$\sigma = ($' + sgm + '$)$' + '\n' + r' ${\rho} = ($' + rho + '$)$' + '\n' +\
# '$w_p = $' +  str(round(obj.x_p,4)) + '\n' +\
# '$m_p = $' +  str(round(obj.m_p,4)) + '\n ' +  r' $F^{\prime}(0) = $' + \
# str(round(obj.dy[0],4))
#             self.lab.append(lab)

#         self.ax_F = self.figF.add_subplot(self.specF[0, 0])
#         self.ax_F.set_title('$F(w)$')
#         self.ax_dF = self.figF.add_subplot(self.specF[0, 1])
#         self.ax_dF.set_title('$F^{\prime}(w)$')
#         self.ax_ddF = self.figF.add_subplot(self.specF[0, 2])
#         self.ax_ddF.set_title('$F^{\prime\prime}(w)$')

#         self.ax_f = self.figf.add_subplot(self.specf[0, 0])
#         self.ax_f.set_title('$f(m)=F(\lambda m)$')
#         self.ax_df = self.figf.add_subplot(self.specf[0, 1])
#         self.ax_df.set_title('$f^{\prime}(m)$')
#         self.ax_ddf = self.figf.add_subplot(self.specf[0, 2])
#         self.ax_ddf.set_title('$f^{\prime\prime}(m)$')

#         self.ax_S = self.figS.add_subplot(self.specS[0, 0])
#         self.ax_S.set_title('$S(m)$')
#         self.ax_dS = self.figS.add_subplot(self.specS[0, 1])
#         self.ax_dS.set_title('$S^{\prime}(m)$')
#         self.ax_ddS = self.figS.add_subplot(self.specS[0, 2])
#         self.ax_ddS.set_title('$S^{\prime\prime}(m)$')

#         self.ax_T = self.figT.add_subplot(self.specT[0, 0])
#         string = '$T(m)=\\mathbb{E}[\\int_0^{\\tau} e^{rt} dt]$'
#         self.ax_T.set_title(string)
#         self.ax_dT = self.figT.add_subplot(self.specT[0, 1])
#         self.ax_dT.set_title('$T^{\prime}(m)$')
#         self.ax_ddT = self.figT.add_subplot(self.specT[0, 2])
#         self.ax_ddT.set_title('$T^{\prime\prime}(m)$')        

#         self.ax_C = self.figC.add_subplot(self.specC[0, 0])
#         self.ax_C.set_title('$C(m)$')
#         self.ax_dC = self.figC.add_subplot(self.specC[0, 1])
#         self.ax_dC.set_title('$C^{\prime}(m)$')
#         self.ax_ddC = self.figC.add_subplot(self.specC[0, 2])
#         self.ax_ddC.set_title('$C^{\prime\prime}(m)$')     

        
#         self.ax_C1 = self.figC1.add_subplot(self.specC1[0, 0])
#         string = '$C(m)={\mu}*T(m)-{\lambda}S(m)-f(m)$'
#         self.ax_C1.set_title(string)
#         self.ax_dC1 = self.figC1.add_subplot(self.specC1[0, 1])
#         self.ax_dC1.set_title('$C^{\prime}(m)$')
#         self.ax_ddC1 = self.figC1.add_subplot(self.specC1[0, 2])
#         self.ax_ddC1.set_title('$C^{\prime\prime}(m)$')            

# #string tells the code what to plot.
# #For example  string contains 'F' or 'f', it only plots all F(w) and f(m).
# #If it contains 'S' otrr 's', it only plots S(m). Otherwise, it plots them all.
#     def plot_iT(self,string):
#         if ('F' in string) | ('f' in string):
#             for indx, obj in enumerate(self.obj):
#                 plt.gca().set_prop_cycle(None)
#                 color = next(self.ax_F._get_lines.prop_cycler)['color']
#                 self.ax_F.plot(self.obj[indx].x,self.obj[indx].y, label = self.lab[indx], color = color)
#                 self.ax_F.plot(self.obj[indx].x,self.obj[indx].line, 'k:')
#                 self.ax_dF.plot(self.obj[indx].x,self.obj[indx].dy, label = self.lab[indx], color = color)
#                 self.ax_ddF.plot(self.obj[indx].x,self.obj[indx].ddy, label = self.lab[indx], color = color)
#             self.ax_F.legend(loc='upper center', \
# shadow = True, ncol = 1, bbox_to_anchor=(1.5, 0.7))
#             tikzplotlib.save(self.path+ "/valueF.tex")

#             for indx, obj in enumerate(self.obj):
#                 plt.gca().set_prop_cycle(None)
#                 self.ax_f.plot(self.obj[indx].m,self.obj[indx].f, label = self.lab[indx], color = color)
#                 self.ax_f.plot(self.obj[indx].m,self.obj[indx].fline, 'k:')
#                 self.ax_df.plot(self.obj[indx].m,self.obj[indx].df, label = self.lab[indx], color = color)
#                 self.ax_ddf.plot(self.obj[indx].m,self.obj[indx].ddf, label = self.lab[indx], color = color)
#             self.ax_f.legend(loc='upper center', \
# shadow = True, ncol = 1, bbox_to_anchor=(1.5, 0.7))
#             tikzplotlib.save(self.path+"/valuef.tex")

#         if ('S' in string) | ('s' in string):
#             for indx, obj in enumerate(self.obj):
#                 plt.gca().set_prop_cycle(None)
#                 color = next(self.ax_S._get_lines.prop_cycler)['color']
#                 self.ax_S.plot(self.obj[indx].mS,self.obj[indx].S, label = self.lab[indx], color = color)
#                 self.ax_dS.plot(self.obj[indx].mS,self.obj[indx].dS, label = self.lab[indx], color = color)
#                 self.ax_ddS.plot(self.obj[indx].mS,self.obj[indx].ddS, label = self.lab[indx], color = color)
#             self.ax_S.legend(loc='upper center', \
# shadow = True, ncol = 1, bbox_to_anchor=(1.5, 0.7))
#             tikzplotlib.save(self.path+"/asset.tex")

#             for indx, obj in enumerate(self.obj):
#                 plt.gca().set_prop_cycle(None)
#                 color = next(self.ax_T._get_lines.prop_cycler)['color']
#                 self.ax_T.plot(self.obj[indx].mS,self.obj[indx].T, label = self.lab[indx], color = color)
#                 self.ax_dT.plot(self.obj[indx].mS,self.obj[indx].dT, label = self.lab[indx], color = color)
#                 self.ax_ddT.plot(self.obj[indx].mS,self.obj[indx].ddT, label = self.lab[indx], color = color)
#             self.ax_T.legend(loc='upper center', \
# shadow = True, ncol = 1, bbox_to_anchor=(1.5, 0.7))
#             tikzplotlib.save(self.path+"/survival.tex")
            
#             for indx, obj in enumerate(self.obj):
#                 plt.gca().set_prop_cycle(None)
#                 color = next(self.ax_C._get_lines.prop_cycler)['color']
#                 self.ax_C.plot(self.obj[indx].mS,self.obj[indx].C, label = self.lab[indx], color = color)
#                 self.ax_dC.plot(self.obj[indx].mS,self.obj[indx].dC, label = self.lab[indx], color = color)
#                 self.ax_ddC.plot(self.obj[indx].mS,self.obj[indx].ddC, label = self.lab[indx], color = color)
#             self.ax_C.legend(loc='upper center', \
# shadow = True, ncol = 1, bbox_to_anchor=(1.5, 0.7))
#             tikzplotlib.save(self.path+"/cost.tex")

#             for indx, obj in enumerate(self.obj):
#                 plt.gca().set_prop_cycle(None)
#                 color = next(self.ax_C1._get_lines.prop_cycler)['color']
#                 self.ax_C1.plot(self.obj[indx].mS,self.obj[indx].C1, label = self.lab[indx], color = color)
#                 self.ax_dC1.plot(self.obj[indx].mS,self.obj[indx].dC1, label = self.lab[indx], color = color)
#                 self.ax_ddC1.plot(self.obj[indx].mS,self.obj[indx].ddC1, label = self.lab[indx], color = color)
#             self.ax_C1.legend(loc='upper center', \
# shadow = True, ncol = 1, bbox_to_anchor=(1.5, 0.7))
#             tikzplotlib.save(self.path+"/cost12.tex")            

# #             f = self.obj.y; lab = '$m_p = $' +  \
# # str(round(self.obj.m_p,4)) + '\n ' +  '$f^{\prime}(0) = $' + \
# # str(round(self.obj.df[0],4))
# #             F_ = self.obj_.y; lab_ = '$m_p = $' +  \
# # str(round(self.obj_.m_p,4)) + '\n ' +  '$f^{\prime}(0) = $' + \
# # str(round(self.obj_.df[0],4))

# #             ax_f.legend(loc='upper center', \
# # shadow = True, ncol = 1, bbox_to_anchor=(1.1, 0.7))

#         if ('S' not in string) & ('s' not in string) & \
# ('F' not in string) & ('f' not in string):
#             self.plot_iT('fs')

#         # plt.draw()
#         plt.show(block=False)
#         plt.pause(1)
#         input("PRESS ENTER TO CLOSE.")
#         plt.close('all')
#         # plt.close(self.figS)
#         # plt.show()
