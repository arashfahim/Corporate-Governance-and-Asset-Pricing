# A finite regime singular control problem in infinite horizon
 
This program aims to solve the singular control problem below
$\sup\mathbb{E}\Big[\int_0^\tau e^{-rt}\big(\mu-\rho_{I_t}-dP_t\big)+e^{-r\tau}X_\tau\Big]$
where $dX_t=\gamma X_t -dP_t + \theta_{I_t}dB_t$, $P_t$ is a singular control process, $I_t$ is a control process taking values in 
a finite set $\{1,...,N\}. There are $N$ possible values for the pair $(\theta_i,\rho_i)$ which are ordered by $0<\theta_1<\cdots<\theta_N$
$\mu>\rho_1>\cdots>\rho_N=0$.
The code solved the variational HJB
$0=\min\Big(\inf_i{\{}-\frac{\theta^2_i}{2}V''-\rho_i{\}}-\gamma x V' + r V - \mu, V' + 1\Big)$

To run the file, you need to run main.py. A UI window appears and you choose up to 6 sets of parameters. Then, you just close the UI window and the program runs to create some verbose and plots the solution and some other related plots.
