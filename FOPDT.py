import time

import math
import numpy as np
import matplotlib.pyplot as plt

class ClosedLoopFOPDT:
    def __init__(
        self,

        DC_gain: float = None, 
        time_constant: float = None, 
        dead_time: float = None,
        
        Kp: float = None, 
        Ki: float = None, 
        Kd: float = None,
        
        PID_min = -np.inf, PID_max = +np.inf,
        integral_min = -np.inf, integral_max = +np.inf):

        # FOPDT parameters
        self._dc_gain = DC_gain
        self._time_constant = time_constant
        self._dead_time = dead_time

        # PID gains
        self._Kp = Kp
        self._Ki = Ki
        self._Kd = Kd

        # Lower and upper bounds for the
        # output of the PID controller
        self._pid_min = PID_min
        self._pid_max = PID_max

        # Lower and upper bounds for the
        # integral of the PID controller
        self._integral_min = integral_min
        self._integral_max = integral_max

        # A dictionary which associates
        # every attribute with a more
        # convenient name
        self._attr_names = {
            'DC_gain': '_dc_gain',
            'time_constant': '_time_constant',
            'dead_time': '_dead_time',
            
            'Kp': '_Kp',
            'Ki': '_Ki',
            'Kd': '_Kd',

            'PID_min': '_pid_min',
            'PID_max': '_pid_max',

            'integral_min': '_integral_min',
            'integral_max': '_integral_max',
        }

    def set_param(self, **kwargs) -> None:
        for key,val in kwargs.items():
            if not hasattr(self, self._attr_names[key]):
                raise AttributeError('Unknown Attribute')
            
            setattr(self, self._attr_names[key], val)

    def get_param(self, key: str) -> float:
        if not hasattr(self, self._attr_names[key]):
            raise AttributeError('Unknown Attribute')

        return getattr(self, self._attr_names[key])

    def autotune_pid_gains(self, method, *args):
        K = self._dc_gain
        τ = self._time_constant
        θ = self._dead_time

        # Internal Model Method
        if method == 'IMC':
            # The desired closed loop time constant
            tau_c = args[0]

            # Proportional gain
            self._Kp = (1/K)*(τ+θ/2)/(tau_c+θ/2) 
        
            # PID gains in standard form
            tau_I = τ + θ/2
            tau_D = (τ*θ)/(2*τ+θ)

            # PID gains in parallel form
            self._Ki = self._Kp / tau_I
            self._Kd = self._Kp * tau_D

        # Cohen Coon Method
        elif method == 'Cohen_Coon':
            # Proportional gain
            self._Kp = (τ/(K*θ))*(4/3+θ/(4*τ))

            # PID gains in standard form
            tau_I = θ*(32+6*θ/τ)/(13+8*θ/τ)
            tau_D = 4*θ/(11+2*θ/τ)

            # PID gains in parallel form
            self._Ki = self._Kp / tau_I
            self._Kd = self._Kp * tau_D

        # Ziegler-Nichols Method
        elif method == 'Ziegler_Nichols':
            # Proportional gain
            self._Kp = 1.1*τ/(K*θ)

            # PID gains in standard form
            tau_I = 2.0*θ
            tau_D = 0.5*θ

            # PID gains in parallel form
            self._Ki = self._Kp / tau_I
            self._Kd = self._Kp * tau_D

        # Chien-Hrones-Reswick Method
        elif method == 'CHR':
            # Proportional gain
            self._Kp = 0.95*τ/(K*θ)

            # PID gains in standard form
            tau_I = 2.40*θ
            tau_D = 0.42*θ

            # PID gains in parallel form
            self._Ki = self._Kp / tau_I
            self._Kd = self._Kp / tau_D

        else:
            raise ValueError('Unknown PID tuner')
            
    def simulate_step_response(
        self,
        t0: float, dt: float, tf: float,
        y0: float, x: float) -> tuple:

        # Arguments check
        if t0 >= tf:
            raise ValueError('Invalid time interval')
        elif dt <= 0:
            raise ValueError('Invalid sampling period')

        # Unpack the FOPDT process parameters
        G = self._dc_gain
        τ = self._time_constant
        θ = self._dead_time

        # Unpack the gains of the PID controller
        Kp, Ki, Kd = self._Kp, self._Ki, self._Kd

        # Unpack the lower and upper bounds of
        # the PID controller output
        pid_min = self._pid_min
        pid_max = self._pid_max

        # Unpack the lower and upper bounds of
        # the PID integral term
        int_min = self._integral_min
        int_max = self._integral_max

        # Parameters Check:
        if pid_min > pid_max:
            raise ValueError('Invalid PID output thresholds')

        elif int_min > int_max:
            raise ValueError('Invalid PID integral thresholds')

        elif τ <= 0:
            raise ValueError('Invalid time constant')

        elif θ <= 0:
            raise ValueError('Invalid dead time')
        
        # If necessary, adjust the sampling period
        # to achieve adequate resolution in the
        # time axis:
        #   -> lag:  index offset for dead time
        #   -> dt:   adjusted sampling period
        #   -> N:    number of time steps
        lag  = max(2,math.ceil(θ/dt))
        dt   = θ/lag
        N    = math.floor((tf-t0)/dt)

        # t: time axis
        # e: feedback-error
        # s: PID integral term
        # u: PID output
        # y: FOPDT process output
        t = [None]*N
        e = [None]*N
        s = [None]*N
        u = [None]*N
        y = [None]*N

        # Initial conditions
        t[0] = t0
        y[0] = y0
        s[0] = 0.
        e[0] = x - y0
        
        u[0] = Kp * e[0] + Ki * s[0]
        u[0] = max(pid_min, min(pid_max, u[0]))

        # A function handle for calculating the
        # gradients of the state variables
        odefun = lambda y,s,u: (-y/τ+(G/τ)*u, x-y)

        # Fourth order Runge-Kutta solver 
        def rk4(i:int)->tuple:
            
            # Retrieve the previous outputs of
            # the PID controller
            u1,u2,u3,u4 = 0.0, 0.0, 0.0, 0.0

            if i>=lag:
                u1 = u[i-lag]           # u(t[i]-θ)
                u4 = u[i+1-lag]         # u(t[i]-θ+dt)
                u2 = 0.5*(u1+u4)        # u(t[i]-θ+dt/2)
                u3 = u2                 # u(t[i]-θ+dt/2)
            
            # interpolation at the beginning
            # of the interval (t[i])
            t1 = t[i]
            y1 = y[i]
            s1 = s[i]
            dy1, ds1 = odefun(y1, s1, u1)

            # 1st interpolation at the midpoint
            # of the interval (t[i]+dt/2)
            t2 = t[i]+0.5*dt
            y2 = y[i]+0.5*dy1*dt
            s2 = s[i]+0.5*ds1*dt
            dy2, ds2 = odefun(y2, s2, u2)

            # 2nd interpolation at the midpoint
            # of the interval (t[i]+dt/2)
            t3 = t2
            y3 = y[i]+0.5*dy2*dt
            s3 = y[i]+0.5*ds2*dt
            dy3, ds3 = odefun(y3, s3, u3)

            # interpolation at the end 
            # of the interval (t[i]+dt)
            t4 = t[i]+1.0*dt
            y4 = y[i]+1.0*dy3*dt
            s4 = y[i]+1.0*ds3*dt
            dy4, ds4 = odefun(y4, s4, u4)

            # Weighted average of the interpolated slopes
            dy = (1*dy1 + 2*dy2 + 2*dy3 + 1*dy4) / 6
            ds = (1*ds1 + 2*ds2 + 2*ds3 + 1*ds4) / 6
            
            return (dy*dt, ds*dt)

        for i in range(N-1):
            t[i+1] = t[i] + dt

            # Numerical integration with the
            # 4th order Runge-Kutta method
            dY, dS = rk4(i)

            # Update the state variables
            # (PID integral and FOPDT output)
            y[i+1] = y[i] + dY
            s[i+1] = s[i] + dS

            # Feedback-error
            e[i+1] = x - y[i+1]

            # Update the PID output
            de     = y[i+1]/τ - (G/τ)*(u[i+1-lag] if i>=lag else 0.0) 
            u[i+1] = Kp*e[i+1] + Ki*s[i+1] + Kd*de 
            u[i+1] = max(pid_min, min(pid_max, u[i+1]))

            # Disable the PID integral term,
            # if the PID output is saturated (anti-reset windup)
            if u[i+1] <= pid_min or u[i+1] >= pid_max:
                s[i+1] = s[i]

        t = np.array(t)
        e = np.array(e)
        s = np.array(s)
        u = np.array(u)
        y = np.array(y)

        return t, (e,s,u,y)
        
    def sweep_param(self, key: str, values: list):
        pass

if __name__ == '__main__':
    # Create a Closed Loop Object
    sys = ClosedLoopFOPDT()

    # Set the parameters of the FOPDT process
    K,tau,theta=1000,450,5

    sys.set_param(DC_gain=K)
    sys.set_param(time_constant=tau)
    sys.set_param(dead_time=theta)

    # Set the gains of the PID controller
    sys.autotune_pid_gains('IMC',5)

    kp = sys.get_param('Kp')
    ki = sys.get_param('Ki')
    kd = sys.get_param('Kd')

    # Some manual fine-tuning is usually required
    ki = 0.28/255
    sys.set_param(Ki=ki)

    print("PID gains:")
    print("Proportional gain: {:.2f}/255".format(kp*255))
    print("Integral gain: {:.2f}/255".format(ki*255))
    print("Derivative gain: {:.2f}/255".format(kd*255))
    print("\n")
    
    # Set the upper and lower bounds for
    # the output of the PID controller
    sys.set_param(PID_min=0.0, PID_max=1.0)

    # Simulation parameters
    # (initial conditions, time step etc...)
    y0 = 0.0
    t0 = 0.0
    dt = 0.25
    tf = 200.0
    x  = 175.0

    # Run the simulation
    t1 = time.time()
    t, (e,s,u,y) = sys.simulate_step_response(t0,dt,tf,y0,x)
    t2 = time.time()

    Mp  = max(0,100*(max(y)-x)/x)
    yss = np.mean(y[t>0.8*tf])

    print("Overshoot: {:.2f}%".format(Mp))
    print("Steady state error: {:.2f}".format(yss-x))
    print("\n")

    # Measure the elapsed time (performance benchmark)
    elapsed_time_ms = 1000*(t2-t1)
    print("Elapsed time: {:.1f} milliseconds".format(elapsed_time_ms))

    # Waveform plots
    fig, ax = plt.subplots(2,2)

    # Plot of feedback error
    ax[0][0].plot(t,e)
    ax[0][0].set(ylabel='Feedback Error')
    ax[0][0].grid()

    # Plot of PID error-integral
    ax[1][1].plot(t,s)
    ax[1][1].set(ylabel='PID integral',xlabel='time')
    ax[1][1].grid()

    # Plot of PID output
    ax[0][1].plot(t,u)
    ax[0][1].set(ylabel='PID output')
    ax[0][1].grid()

    # Plot of FOPDT output
    ax[1][0].plot(t,y)
    ax[1][0].set(ylabel='FOPDT output',xlabel='time')
    ax[1][0].grid()

    plt.show()
