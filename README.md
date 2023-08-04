# A finite regime singular control problem in infinite horizon
 
This program aims to solve the singular control problem below
$\sup\mathbb{E}\Big[\int_0^\tau e^{-rt}\big(\mu-\rho_{I_t}-dP_t\big)+e^{-r\tau}X_\tau\Big]$
where $dX_t=\gamma X_t -dP_t + \theta_{I_t}dB_t$, and the supremum is over all nondecreasing adapted processes $P_t$, possibly singular, and all adapted processes $I_t$ taking values in a finite set $[N]:=\lbrace1,...,N\rbrace$. There are $N$ possible values for the pair $(\theta_i,\rho_i)$ which are ordered by $0<\theta_1<\cdots<\theta_N$
$\mu>\rho_1>\cdots>\rho_N=0$.

This problem arises in the optimal case of the agency problem with monitoring in theoretical economics. The principal obseves te process $X$ and have freedom to choose $\theta_i$ for $i\in [N]$. The small $\theta_i$ is, the higher the level of monitoring is. While monitoring level $\theta_i$ is applied, the principal pays a sunning cost at rate $\rho_i$.  Process $P_t$ is the cumulative payment made by principal to the agent. 

The code first solves the variational HJB

$$0=\min\Big(\inf_i\lbrace-\frac{\theta^2_i}{2}V^{\prime\prime}-\rho_i\rbrace-\gamma x V^{\prime} + r V - \mu, V^{\prime} + 1\Big)$$ 

by turing it into a systm of nonlinear ODEs.

$$\begin{cases}
V^{\prime}=U\\
U^{\prime}= H^{-1}(\gamma x U - r V + \mu)\\
V(0)=0\\
U(0)= u_0
\end{cases}$$

where $H(\Gamma)=\inf_i\lbrace-\frac{\theta^2_i}{2}\Gamma-\rho_i\rbrace$. Then, we determine triplet $(x_0,x_1,x_2)$ such that $V(x_0)=\mu/r-\gamma/r x_0$, $U(x_1)+1=0$, and $U^{\prime}(x_2)=0$. By iterating over $u_0$, we converge towrad achieving $x_0=x_1=x_2$. Upon achieving this goal, $V$ is the solution to the variational  HJB and $U$ and $U^\prime$ are its first and second derivatives. 

# How to run the program:

To run the file, you need to run main.py. A UI window appears and you choose up to 6 sets of parameters. Then, you just close the UI window and the program runs to create some verbose and plots the solution and some other related plots.
