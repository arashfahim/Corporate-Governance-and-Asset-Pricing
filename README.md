# A finite regime singular control problem in infinite horizon
 
This program aims to solve the singular control problem below

$$\sup\mathbb{E}\Big[\int_0^\tau e^{-rt}\big(\mu-\rho_{I_t}-dP_t\big)+e^{-r\tau}X_\tau\Big]$$

where $dX_t=\gamma X_t -dP_t + \theta_{I_t}dB_t$, and the supremum is over all nondecreasing adapted processes $P_t$, possibly singular, and all adapted processes $I_t$ taking values in a finite set $[N]:=\lbrace1,...,N\rbrace$. There are $N$ possible values for the pair $(\theta_i,\rho_i)$ which are ordered by $0<\theta_1<\cdots<\theta_N$
$\mu>\rho_1>\cdots>\rho_N=0$. The loser the level of noise, $\theta_i$, is, the higher $\rho_i$. $\rho_i$ are considered a cost associated to level $i$.

This problem arises in the optimal case of the agency problem with monitoring in theoretical economics. The principal obseves te process $X$ and have freedom to choose $\theta_i$ for $i\in [N]$. The small $\theta_i$ is, the higher the level of monitoring is. While monitoring level $\theta_i$ is applied, the principal pays a running cost at rate $\rho_i$. In the problem $\theta_i=\lambda\sigma_i$, where $1-\lambda$ is a proportional  cost of stealing from the company by the agent. The higher $\lambda$ is, the lower the agent benefits from stealing.

Process $P_t$ is the cumulative payment made by principal to the agent. 

The code first solves the variational HJB

$$0=\min\Big(\inf_i\lbrace-\frac{\theta^2_i}{2}V^{\prime\prime}-\rho_i\rbrace-\gamma x V^{\prime} + r V - \mu, V^{\prime} + 1\Big)$$ 

by turing it into a systm of nonlinear ODEs.

$$\begin{cases}
V^{\prime}=U\\
U^{\prime}= H^{-1}(\gamma x U - r V + \mu)\\
V(0)=0\\
U(0)= u_0
\end{cases}$$

where $H(\Gamma)=\inf_i\lbrace-\frac{\theta^2_i}{2}\Gamma-\rho_i\rbrace$. Then, we determine triplet $(x_0,x_1,x_2)$ such that $V(x_0)=\mu/r-\gamma/r x_0$, $U(x_1)+1=0$, and $U^{\prime}(x_2)=0$. By iterating over $u_0$, we converge towrad achieving $x_0=x_1=x_2$. Upon achieving this goal, $V$ is the solution to the variational  HJB and $U$ and $U^\prime$ are its first and second derivatives of $V$. 

After finding $V$ and $V^{\prime\prime}$, the optimal contract is described as below: 

1. Optimal monitoring strategy: $I(x)=argmin_i \lbrace -\frac{\theta^2_i}{2} V^{\prime\prime}(x)-\rho_i\rbrace $, and, therefore, $\theta(x)=\lambda\sigma_{I(x)}$ and  $\rho(x)=\rho_{I(x)}$
2. Optimal payment: $P_t$ is the local time of $dX_t=\gamma X_t + \theta(X_t)dB_t$ at point $x_0$, described by $V(x_0)=\mu/r-\gamma/r x_0$, $V^{\prime}(x_0)+1=0$, and $V^{\prime\prime}(x_0)=0$.

# Implementation of the contract
To implement the contract, principal create a cash reserve given by $M_t=\frac{X_t}{\lambda}$ and transfers a $\lambda$ portion of ownership of the company to the agent by issuing a non-transferable security that pays dividend when $M_t$ hits  $m_0=\frac{x_0}{\lambda}$ at rate equal to the local time $\lambda^{-1}P_t$.

The value of the principal as a function of $M_t$ is given by

$$f(m):=\mathbb{E}\Big[\int_0^\tau e^{-rt}\big(\mu-\rho(M_t)-dP_t\big)\Big]$$

where $\rho(m)=\rho(\lambda m)$, $\sigma(m)=\sigma(\lambda m)$, and $dM_t=\gamma M_t dt -\lambda^{-1}dP_t + \sigma(M_t)dB_t$

## Price of the asset of the company

The value of the security, issued by the principal is given by 

$$S(m)=\lambda^{-1}\mathbb{E}\Big[\int_0^\tau e^{-rt} dP_t\Big]$$

and satisfies the ODE

$$\begin{cases}
-\frac{\sigma^2(m)}{2}S^{\prime\prime}-\gamma m S^{\prime} + r S=0\\
S(0)=0\\
S^\prime(m_0)= 1
\end{cases}$$

## Credit worthiness of the company

The credit worthiness of the company can be measured by 

$$\mathbb{E}[e^{-r\tau} ]$$

The function $T(m)=\mathbb{E}\Big[\int_0^\tau e^{-rt} dt\Big]=r^{-1}(1-\mathbb{E}[e^{-r\tau}])$, measurng survival of the company, satisfiess the ODE 

$$\begin{cases}
-\frac{\sigma^2(m)}{2}T^{\prime\prime}-\gamma m T^{\prime} + r T -1=0\\
T(0)=0\\
T^\prime(m_0)= 0
\end{cases}$$


## Monitoring cost

The cost of monitoring is measured by 

$$T(m)=\mathbb{E}\Big[\int_0^\tau e^{-rt} \rho(M_t) dt\Big]$$



# How to run the program:

To run the file, you need to run main.py. A UI window appears and you choose up to 6 sets of parameters. Then, you just close the UI window and the program runs to create some verbose and plots the solution and some other related plots.
