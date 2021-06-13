import solver as slv
from absl import app
from absl import flags
from absl import logging as absl_logging
from scipy import integrate


def main(argv):
    ode = slv.my_ode() #getattr(slv,'my_ode')()
    print(r"The ODE is:\n dX(1)/dt = [",ode.param[0], "X(1) + ", ode.param[1], "X(2)")
    print(r"              dX(2)/dt = [",ode.param[2], "X(1) + ", ode.param[3], "X(2)")
    param = ode.param
    cond = ode.cond
    t = ode.t
    ode.x = ode.solve()
    print(ode.x)

if __name__ == '__main__':
    app.run(main)
