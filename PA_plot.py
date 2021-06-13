import matplotlib.pyplot as plt
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica"]})
import matplotlib.mlab as mlab
import matplotlib.gridspec as gridspec


class plotF(object):
    def __init__(self,PA):
        self.obj = PA
        self.ind = len(PA)
        self.figF = plt.figure(figsize=(15, 7.5),constrained_layout=True)
        self.specF = gridspec.GridSpec(ncols=3, nrows=2, figure=self.figF)
        self.figS = plt.figure(figsize=(10, 5),constrained_layout=True)
        self.specS = gridspec.GridSpec(ncols=3, nrows=1, figure=self.figS)
        self.lab = []
        rho = []
        sgm = []
        for indx, obj in enumerate(self.obj):
            # for r in obj.param[5]:
            #     rho += [str(r)]
            # for s in obj.param[4]:
            #     sgm += [str(s)]
            # rho = ','.join(rho)
            # sgm = ','.join(sgm)
            rho = str(obj.param[5]).rstrip("]").lstrip("[")
            sgm = str(obj.param[4]).rstrip("]").lstrip("[")
            mu = str(round(obj.param[0],2))
            gamma = str(round(obj.param[1],2))
            R = str(round(obj.param[2],2))
            lmb = str(round(obj.param[3],2))
            lab = 'Parameters: ' + r'$\mu = $' + mu + '\n' +\
r'$\gamma = $' + gamma + ' $r = $' + R + '\n' + r'$ \lambda = $' + lmb + '\n' +\
r'$\sigma = ($' + sgm + '$)$' + '\n' + r' ${\rho} = ($' + rho + '$)$' + '\n' +\
'$w_p = $' +  str(round(obj.x_p,2)) + '\n' +\
'$m_p = $' +  str(round(obj.m_p,2)) + '\n ' +  r' $F^{\prime}(0) = $' + \
str(round(obj.dy[0],2))
            self.lab.append(lab)

        self.ax_F = self.figF.add_subplot(self.specF[0, 0])
        self.ax_F.set_title('$F(w)$')
        self.ax_dF = self.figF.add_subplot(self.specF[0, 1])
        self.ax_dF.set_title('$F^{\prime}(w)$')
        self.ax_ddF = self.figF.add_subplot(self.specF[0, 2])
        self.ax_ddF.set_title('$F^{\prime\prime}(w)$')
        self.ax_f = self.figF.add_subplot(self.specF[1, 0])
        self.ax_f.set_title('$f(m)$')
        self.ax_df = self.figF.add_subplot(self.specF[1, 1])
        self.ax_df.set_title('$f^{\prime}(m)$')
        self.ax_ddf = self.figF.add_subplot(self.specF[1, 2])
        self.ax_ddf.set_title('$f^{\prime\prime}(m)$')

        self.ax_S = self.figS.add_subplot(self.specS[0, 0])
        self.ax_S.set_title('$S(m)$')
        self.ax_dS = self.figS.add_subplot(self.specS[0, 1])
        self.ax_dS.set_title('$S^{\prime}(m)$')
        self.ax_ddS = self.figS.add_subplot(self.specS[0, 2])
        self.ax_ddS.set_title('$S^{\prime\prime}(m)$')

#string tells the code what to plot.
#For example  string contains 'F' or 'f', it only plots all F(w) and f(m).
#If it contains 'S' otrr 's', it only plots S(m). Otherwise, it plots them all.
    def plot_iT(self,string):
        if ('F' in string) | ('f' in string):
            for indx, obj in enumerate(self.obj):
                plt.gca().set_prop_cycle(None)
                color = next(self.ax_F._get_lines.prop_cycler)['color']
                self.ax_F.plot(self.obj[indx].x,self.obj[indx].y, label = self.lab[indx], color = color)
                self.ax_F.plot(self.obj[indx].x,self.obj[indx].line, 'k:')
                self.ax_dF.plot(self.obj[indx].x,self.obj[indx].dy, label = self.lab[indx], color = color)
                self.ax_ddF.plot(self.obj[indx].x,self.obj[indx].ddy, label = self.lab[indx], color = color)
                self.ax_f.plot(self.obj[indx].m,self.obj[indx].f, label = self.lab[indx], color = color)
                self.ax_f.plot(self.obj[indx].m,self.obj[indx].fline, 'k:')
                self.ax_df.plot(self.obj[indx].m,self.obj[indx].df, label = self.lab[indx], color = color)
                self.ax_ddf.plot(self.obj[indx].m,self.obj[indx].ddf, label = self.lab[indx], color = color)
            self.ax_f.legend(loc='upper center', \
shadow = True, ncol = 1, bbox_to_anchor=(1.5, 1.4))

        if ('S' in string) | ('s' in string):
            for indx, obj in enumerate(self.obj):
                plt.gca().set_prop_cycle(None)
                color = next(self.ax_S._get_lines.prop_cycler)['color']
                self.ax_S.plot(self.obj[indx].mS,self.obj[indx].S, label = self.lab[indx], color = color)
                self.ax_dS.plot(self.obj[indx].mS,self.obj[indx].dS, label = self.lab[indx], color = color)
                self.ax_ddS.plot(self.obj[indx].mS,self.obj[indx].ddS, label = self.lab[indx], color = color)
            self.ax_S.legend(loc='upper center', \
shadow = True, ncol = 1, bbox_to_anchor=(1.5, 0.7))

#             f = self.obj.y; lab = '$m_p = $' +  \
# str(round(self.obj.m_p,2)) + '\n ' +  '$f^{\prime}(0) = $' + \
# str(round(self.obj.df[0],2))
#             F_ = self.obj_.y; lab_ = '$m_p = $' +  \
# str(round(self.obj_.m_p,2)) + '\n ' +  '$f^{\prime}(0) = $' + \
# str(round(self.obj_.df[0],2))

#             ax_f.legend(loc='upper center', \
# shadow = True, ncol = 1, bbox_to_anchor=(1.1, 0.7))

        if ('S' not in string) & ('s' not in string) & \
('F' not in string) & ('f' not in string):
            self.plot_iT('fs')

        # plt.draw()
        plt.show(block=False)
        plt.pause(1)
        input("PRESS ENTER TO CLOSE.")
        plt.close('all')
        # plt.close(self.figS)
        # plt.show()
