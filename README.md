# A finite regime singular control problem in infinite horizon
 
This program aims to solve the singular control problem below
$\sup\mathbb{E}\Big[\int_0^\tau e^{-rt}\big(\mu-\rho_{I_t}-dP_t\big)+e^{-r\tau}X_\tau\Big]$
where $dX_t=\gamma X_t -dP_t + \theta_{I_t}dB_t$, and the supremum is over all nondecreasing adapted processes $P_t$, possibly singular, and all adapted processes $I_t$ taking values in a finite set $[N]:=\lbrace1,...,N\rbrace$. There are $N$ possible values for the pair $(\theta_i,\rho_i)$ which are ordered by $0<\theta_1<\cdots<\theta_N$
$\mu>\rho_1>\cdots>\rho_N=0$.

This problem arises in the optimal case of the agency problem with monitoring in theoretical economics. The principal obseves te process $X$ and have freedom to choose $\theta_i$ for $i\in [N]$. The small $\theta_i$ is, the higher the level of monitoring is. While monitoring level $\theta_i$ is applied, the principal pays a sunning cost at rate $\rho_i$.  Process $P_t$ is the cumulative payment made by principal to the agent. 

The code first solves the variational HJB
$0=\min\Big(\inf_i\lbrace-\frac{\theta^2_i}{2}V^{\prime\prime}-\rho_i\rbrace-\gamma x V^{\prime} + r V - \mu, V^{\prime} + 1\Big)$ by turing it into a systm of nonlinear ODEs.


$$\begin{cases}
U= V^{\prime}\\
U^{\prime}= H^{-1}(\gamma x V^{\prime} - r V + \mu)
\end{cases}$$
where $H(\Gamma)=\inf_i\lbrace-\frac{\theta^2_i}{2}\Gamma-\rho_i\rbrace$ 

To run the file, you need to run main.py. A UI window appears and you choose up to 6 sets of parameters. Then, you just close the UI window and the program runs to create some verbose and plots the solution and some other related plots.
