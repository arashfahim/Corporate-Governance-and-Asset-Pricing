# A finite regime singular control problem in infinite horizon
 
This program aims to solve the singular control problem below

$$\sup\mathbb{E}\Big[\int_0^\tau e^{-rt}\big(\mu-\rho_{I_t}-dP_t\big)+e^{-r\tau}X_\tau\Big]$$

where $dX_t=\gamma X_t -dP_t + \theta_{I_t}dB_t$, and the supremum is over all nondecreasing adapted processes $P_t$, possibly singular, and all adapted processes $I_t$ taking values in a finite set $[N]:=\lbrace1,...,N\rbrace$. There are $N$ possible values for the pair $(\theta_i,\rho_i)$ which are ordered by $0<\theta_1<\cdots<\theta_N$
$\mu>\rho_1>\cdots>\rho_N=0$. The loser the level of noise, $\theta_i$, is, the higher $\rho_i$ is. $\rho_i$ are considered a cost associated to level $i$.

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

In each iteration of the code, the nonlinear system of ODEs is solved by ```from scipy import integrate``` .

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

$$C(m)=\mathbb{E}\Big[\int_0^\tau e^{-rt} \rho(M_t) dt\Big]$$

$$\begin{cases}
-\frac{\sigma^2(m)}{2}C^{\prime\prime}-\gamma m C^{\prime} + r C -\rho(m)=0\\
C(0)=0\\
C^\prime(m_0)= 0
\end{cases}$$


## An identity

While the program solves each of the above functions separately, it also passes the verification step that the identity below holds true.

$$f(m)=\mu T(m) - C(m) -\lambda S(m)$$

## Solving the boundary value problems for $S$, $T$, and $C$

We use the package ```from scipy.integrate import solve_bvp``` to solve the ODEs above. Note that the ODEs are are linear but with   discontinuous coefficients.




# How to run the program:

These codes are checked with `Python 3.12.4`, `numpy 1.26.4`, `scipy 1.13.1`, `matplotlib 3.10`, `json 2.0.9`, and `munch 2.5.0`. 

To run the file, you need to run ```python main.py``` in a location which has a copy of```equation.py```, ```gui.py```, and ```PA_plot.py```. A UI window appears and you choose up to 6 sets of parameters. Then, you must click on the 'Save Parameters' to save your inserted parameters as a `json` file. To run the code, close the UI window. Then, the program runs to solve the problems and show verbose at each step of the solution. Finally, it plots the solutions using ```PA_plot.py```. Note that the  ```PA_plot.myfigures``` allows for ploting many combined functions such as $\lambda S(m)$, or $C(m)+\lambda S(m)$.

If you do not need to plot the solutions, you exclude ```pap.myfigures``` in ```main.py```. The solutions to the boundary value problems are also written in the file ```output.dat```.

Package ```tikzplotlib``` that was used in the earlier versions are deprecated and the code yields png files. To plot in higher resolution in $\LaTeX$, you need to use ```output.dat```.


