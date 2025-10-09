import matplotlib.pyplot as plt #version 3.2 consistancy fails on 3.5
plt.rcParams.update({
    "text.usetex": True,
    # "font.family": "sans-serif",
    # "font.sans-serif": ["Helvetica"],
    'text.latex.preamble' : r'\usepackage{amssymb} \usepackage{amsmath}' #for \text command
    })
import matplotlib.mlab as mlab
import matplotlib.gridspec as gridspec
# import tikzplotlib
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
                plt.savefig(self.path+ "/valueFw.png")

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
                plt.savefig(self.path+"/valuefm.png")

            if self.string == '$S(m)$':
                for indx, obj in enumerate(self.obj):
                    plt.gca().set_prop_cycle(None)
                    color = next(self.ax._get_lines.prop_cycler)['color']
                    self.ax.plot(self.obj[indx].mS,self.obj[indx].S, label = self.lab[indx], color = color)
                    self.ax_d.plot(self.obj[indx].mS,self.obj[indx].dS, label = self.lab[indx], color = color)
                    self.ax_dd.plot(self.obj[indx].mS,self.obj[indx].ddS, label = self.lab[indx], color = color)
                self.ax.legend(loc='upper center', \
    shadow = True, ncol = 1, bbox_to_anchor=(1.5, 0.7))
                plt.savefig(self.path+"/asset.png")

            if self.string == '$T(m)$':
                for indx, obj in enumerate(self.obj):
                    plt.gca().set_prop_cycle(None)
                    color = next(self.ax._get_lines.prop_cycler)['color']
                    self.ax.plot(self.obj[indx].mS,self.obj[indx].T, label = self.lab[indx], color = color)
                    self.ax_d.plot(self.obj[indx].mS,self.obj[indx].dT, label = self.lab[indx], color = color)
                    self.ax_dd.plot(self.obj[indx].mS,self.obj[indx].ddT, label = self.lab[indx], color = color)
                self.ax.legend(loc='upper center', \
    shadow = True, ncol = 1, bbox_to_anchor=(1.5, 0.7))
                plt.savefig(self.path+"/survival.png")
            
            if self.string == '$C(m)$':                
                for indx, obj in enumerate(self.obj):
                    plt.gca().set_prop_cycle(None)
                    color = next(self.ax._get_lines.prop_cycler)['color']
                    self.ax.plot(self.obj[indx].mS,self.obj[indx].C, label = self.lab[indx], color = color)
                    self.ax_d.plot(self.obj[indx].mS,self.obj[indx].dC, label = self.lab[indx], color = color)
                    self.ax_dd.plot(self.obj[indx].mS,self.obj[indx].ddC, label = self.lab[indx], color = color)
                self.ax.legend(loc='upper center', \
    shadow = True, ncol = 1, bbox_to_anchor=(1.5, 0.7))
                plt.savefig(self.path+"/cost.png")

            if self.string == '$C(m)=\mu*T(m)-\lambda S(m)-f(m)$': 
                for indx, obj in enumerate(self.obj):
                    plt.gca().set_prop_cycle(None)
                    color = next(self.ax._get_lines.prop_cycler)['color']
                    self.ax.plot(self.obj[indx].mS,self.obj[indx].C1, label = self.lab[indx], color = color)
                    self.ax_d.plot(self.obj[indx].mS,self.obj[indx].dC1, label = self.lab[indx], color = color)
                    self.ax_dd.plot(self.obj[indx].mS,self.obj[indx].ddC1, label = self.lab[indx], color = color)
                self.ax.legend(loc='upper center', \
    shadow = True, ncol = 1, bbox_to_anchor=(1.5, 0.7))
                plt.savefig(self.path+"/cost1.png")            

            if self.string == '$C(m)+\lambda S(m)$': 
                for indx, obj in enumerate(self.obj):
                    plt.gca().set_prop_cycle(None)
                    color = next(self.ax._get_lines.prop_cycler)['color']
                    self.ax.plot(self.obj[indx].mS,self.obj[indx].C+self.obj[indx].param[3]*self.obj[indx].S, label = self.lab[indx], color = color)
                    self.ax_d.plot(self.obj[indx].mS,self.obj[indx].dC+self.obj[indx].param[3]*self.obj[indx].dS, label = self.lab[indx], color = color)
                    self.ax_dd.plot(self.obj[indx].mS,self.obj[indx].ddC+self.obj[indx].param[3]*self.obj[indx].ddS, label = self.lab[indx], color = color)
                self.ax.legend(loc='upper center', \
    shadow = True, ncol = 1, bbox_to_anchor=(1.5, 0.7))
                plt.savefig(self.path+"/cost_lambda_asset.png")     

            if self.string == '$\lambda S(m)$': 
                for indx, obj in enumerate(self.obj):
                    plt.gca().set_prop_cycle(None)
                    color = next(self.ax._get_lines.prop_cycler)['color']
                    self.ax.plot(self.obj[indx].mS,self.obj[indx].param[3]*self.obj[indx].S, label = self.lab[indx], color = color)
                    self.ax_d.plot(self.obj[indx].mS,self.obj[indx].param[3]*self.obj[indx].dS, label = self.lab[indx], color = color)
                    self.ax_dd.plot(self.obj[indx].mS,self.obj[indx].param[3]*self.obj[indx].ddS, label = self.lab[indx], color = color)
                self.ax.legend(loc='upper center', \
    shadow = True, ncol = 1, bbox_to_anchor=(1.5, 0.7))
                plt.savefig(self.path+"/lambda_asset.png")   

            plt.show(block=False)
            plt.pause(1)
    
    def close_all(self):
            input("PRESS ENTER TO CLOSE.")
            plt.close('all')    

#######################
