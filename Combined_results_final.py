import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

# Each stage of the process
def stage_1(v0,S0,J,a_6):
    # initial condition

    t1 = np.sqrt((2 * S0) / J)
    # time points
    t = np.linspace(0, t1)
    def velocity(v,t):
        dvdt = a_6 - (a_6*t)/t1
        # plot results
        return dvdt

    # solve ODE
    v = odeint(velocity, v0/3.1, t)
    v1 = v[-1]
    # plot results
    #plt.plot(t, v)
    #plt.xlabel('time')
    #plt.ylabel('v(t)')
    Q_1 = w*A*(v0*t1 + a_6*t1**2/3)
    #plt.show()
    return v1, t1,Q_1
def stage_2(v1, H_instant, H_static, g, p, k, D,E,e,w,Av,A,Ev):
    a = np.sqrt((1)/(p/k)+ (D/E*e))
    delta_H = H_instant - H_static
    v2 = v1 - (delta_H *g)/a
    Z = (a*w*Av**2)/(A*g*Ev)
    t2 = Z*np.log(v1/v2)
    return v2, t2, delta_H, Z,a
def stage_3(t2,v1,L,m,a,h,H,Z):
    tb = 2*L/a - t2*1.1315
    delta_V = (4*g*(h-H) + m*v1)/(4*a+m) - 1.1315
    vp = (v1 - delta_V)
    tp = Z * np.log((2*v1)/(vp))
    t3 = tb + tp
    N = 6
    #Water pumped in this section
    Q_p = w*Ac*(vp**2/20)
    #Q_p = w*A*(N*delta_V*tb +(N-1)*2*Z*delta_V +Z*vp - N**2*delta_V*tb -(N-1)**2*delta_V*t2 -(N-1)*v1*t2-(2*v1-vp)*tp)
    return vp, t3, delta_V,Q_p
def stage_4x(delta_V,a,Z,v3):
    v4 = np.sqrt(((a*Z)/L)*(delta_V**2))
    t4 = 2*Z/(v4-v3) + 2*L/a
    return v4,t4
def stage_4(delta_V,a,Z,v3):
    v4 = np.sqrt(((a*Z)/L)*(v3**2))
    t4 = -2*Z/(v4-v3) + 2*L/a
    return v4,t4
def stage_5(X,j,v4,v0):
    t5 = (L/(j*X))*np.log(((X+v0)/(X-v0))*((X-v4)/(X+v4)))
    v5 = v0
    #Water pumped in final section
    Q_5 = (w*A*L/j) * np.log((X**2 - v4**2)/(X**2 - v0**2))
    return v5 ,t5,Q_5
def stage_3x(t2,v1,L,m,a,h,H,Z):
    tb = 2*L/a - t2
    delta_V = (4*g*(h-H) + m*v1)/(4*a+m)
    vp = (v1 - delta_V)
    tp = Z* np.log((2*v1)/(2*v1 - vp))
    t3 = tb + tp
    N = 6
    #Water pumped in this section
    Q_p = w*Ac*(vp**2/20)
    #Q_p = w*A*(N*delta_V*tb +(N-1)*2*Z*delta_V +Z*vp - N**2*delta_V*tb -(N-1)**2*delta_V*t2 -(N-1)*v1*t2-(2*v1-vp)*tp)
    return vp , t3, delta_V,Q_p
def Overall_results(t1,t2,t3,t4,t5,Q_1,Q_5,Q_p):
    W_wasted = Q_1 + Q_5
    W_pumped = Q_p
    Efficiency = W_pumped/(W_wasted+W_pumped) * 100
    #print(W_wasted)
    return W_wasted,t1+t2+t3+t4+t5
def full_process(v0,a0):
    #Stage 1 processes
    v1,t1,Q_1 = stage_1(v0,S0,J,a0)

    #Stage 2 process
    v2,t2,delta_H,Z,a = stage_2(v1, H_instant, H_static, g, p, k, D,E,e,w,Av,A,Ev)

    #Stage 3 process
    v3, t3, delta_V,Q_p = stage_3(t2,v1,L,m,a,h,H_static,Z)
    #Stage 4 process
    v4,t4 = stage_4(delta_V,a,Z,v3)

    #Stage 5 process

    v5,t5,Q_5 = stage_5(X,j,v4,v0)
    return v1,t1,v2,t2,v3,t3,v4,t4,v5,t5, Q_1,Q_5,Q_p
def full_processx(v0,a0):
    #Stage 1 processes
    v1,t1,Q_1 = stage_1(v0,S0,J,a0)

    #Stage 2 process
    v2,t2,delta_H,Z,a = stage_2(v1, H_instant, H_static, g, p, k, D,E,e,w,Av,A,Ev)

    #Stage 3 process
    v3, t3, delta_V,Q_p = stage_3x(t2,v1,L,m,a,h,H_static,Z)
    #Stage 4 process
    v4,t4 = stage_4x(delta_V,a,Z,v3)

    #Stage 5 process

    v5,t5,Q_5 = stage_5(X,j,v4,v0)
    return v1,t1,v2,t2,v3,t3,v4,t4,v5,t5, Q_1,Q_5,Q_p


#Stage 1 parameters
S0 = 0.0161
J = 4.0

#Stage 2 Parameters
H_static = 9.2
H_instant = 65
g = 32.2
p = 62.4
k = 2200000000
D = 0.166667
E = 103
e = 0.153937
w = 62.43
Av = 0.1043
Ac = 0.018241
A = 0.0233
Ev = 3870000

#Stage 3 Parameters
L = 54.8
m = 817
#Stage 5 parameters
j = 15.5
X = np.sqrt(2*g*H_static/j)

def rankin_eff():

    #Rankine Efficiency
    v0 = 1.64
    a0 = (X ** 2 -v0 ** 2) / (2 * L / j)

    # Testing involving varying h to get the quantity pumped.

    H_array1 = []
    rankine_array = []
    for h in range(0, 300, 10):
        v1, t1, v2, t2, v3, t3, v4, t4, v5, t5, Q_1, Q_5, Q_p = full_process(v0, a0)
        w_waste1, tx = Overall_results(t1, t2, t3, t4, t5, Q_1, Q_5, Q_p)
        rankine_eff = (Q_p*(h-H_static))/(w_waste1*H_static)*10
        H_array1.append(h/3.281)
        rankine_array.append(rankine_eff)


    #Rankine Efficiency
    vx = 3.28
    ax = (X ** 2 - vx ** 2) / (2 * L / j)

    # Testing involving varying h to get the quantity pumped.

    H_array2 = []
    rankine_array2 = []
    for h in range(0, 600, 10):
        v1, t1, v2, t2, v3, t3, v4, t4, v5, t5, Q_1, Q_5, Q_p = full_process(vx, ax)
        w_waste2, tx = Overall_results(t1, t2, t3, t4, t5, Q_1, Q_5, Q_p)
        rankine_eff = (Q_p*(h-H_static))/(w_waste2*H_static)*10
        H_array2.append(h/3.281)
        rankine_array2.append(rankine_eff)



    #Rankine Efficiency
    vy = 4.92
    ay = (X ** 2 - vy ** 2) / (2 * L / j)

    # Testing involving varying h to get the quantity pumped.

    H_array3 = []
    rankine_array3 = []
    for h in range(0, 800, 10):
        v1, t1, v2, t2, v3, t3, v4, t4, v5, t5, Q_1, Q_5, Q_p = full_process(vy, ay)
        w_waste3, tx = Overall_results(t1, t2, t3, t4, t5, Q_1, Q_5, Q_p)
        rankine_eff = (Q_p*(h-H_static))/(w_waste3*H_static)*10
        H_array3.append(h/3.281)

        rankine_array3.append(rankine_eff)

    plt.plot(H_array1,rankine_array,label = 'v0 = 0.5 m/sec')
    plt.plot(H_array2,rankine_array2,label = 'v0 = 1.0 m/sec')
    plt.plot(H_array3,rankine_array3,label = 'v0 = 1.5 m/sec')
    plt.grid()
    plt.legend()
    plt.xlabel('The static delivery head, h, (metres)')
    plt.ylabel('The Rankine Efficiency of the hydraulic ram pump (%)')
    plt.show()
h = 65
#Tesing the initial velocity against the maximum water pumped
v_array1 = []
v_vary_waste_array1 = []
v_vary_pumped_array = []
def v_vary_waste():
    for i in range(20,80):
        v0 = 0.1*i
        a0 = (X ** 2 - v0 ** 2) / (2 * L / j)
        v1, t1, v2, t2, v3, t3, v4, t4, v5, t5, Q_1, Q_5, Q_p = full_process(v0,a0)
        v1x, t1x, v2x, t2x, v3x, t3x, v4x, t4x, v5x, t5x, Q_1x, Q_5x, Q_px = full_processx(v0,a0)

        # Overall Process
        w_waste1,tx = Overall_results(t1, t2, t3, t4, t5, Q_1, Q_5, Q_p)
        v_array1.append(v0)
        v_vary_pumped_array.append(Q_p/2.204623)
        v_vary_waste_array1.append(w_waste1/2.204623)
        print('the value of v0 is', v0, 'and water wasted is ',w_waste1 )
    plt.plot(v_array1,v_vary_waste_array1,label = 'Water wasted per cycle,Q_w')
#    plt.plot(v_array1,v_vary_pumped_array, label = 'Water pumped per cycle,Q_p')
    plt.xlabel('Velocity threshold of the waste valve before movement occurs, v0')
    plt.ylabel('Water wasted/pumped per cycle (litres)')
    plt.legend()
    plt.grid()
    plt.show()
def pumped_wasted_etc():
    v0 = 1.64
    a0 = (X ** 2 - v0 ** 2) / (2 * L / j)
    time5 = []
    # Testing involving varying h to get the quantity pumped.

    H_array1 = []
    Pumped_array1 = []
    wasted_array1 = []
    time_array1 = []
    for h in range(20,350,10):
        v1,t1,v2,t2,v3,t3,v4,t4,v5,t5,Q_1,Q_5,Q_p = full_process(v0,a0)
        w_waste1, tx = Overall_results(t1, t2, t3, t4, t5, Q_1, Q_5, Q_p)
        H_array1.append(h/3.281)
        Pumped_array1.append(Q_p/2.204623)
        wasted_array1.append(w_waste1/2.204623)
        time_array1.append(tx)

    #second for graphs
    vx = 3.28
    ax = (X ** 2 - v0 ** 2) / (2 * L / j)

    # Testing involving varying h to get the quantity pumped.
    #Second variable
    H_array2 = []
    Pumped_array2 = []
    wasted_array2 = []
    time_array2 = []
    for h in range(20,700,10):
        v1,t1,v2,t2,v3,t3,v4,t4,v5,t5,Q_1,Q_5,Q_p = full_process(vx,ax)
        w_waste2, ty = Overall_results(t1, t2, t3, t4, t5, Q_1, Q_5, Q_p)
        H_array2.append(h/3.281)
        Pumped_array2.append(Q_p/2.204623)
        wasted_array2.append(w_waste2/2.204623)
        time_array2.append(ty)
    vy = 4.92
    ay = (X ** 2 - v0 ** 2) / (2 * L / j)

        # Testing involving varying h to get the quantity pumped.
    #Third variable
    H_array3 = []
    Pumped_array3 = []
    wasted_array3 = []
    time_array3 = []
    for h in range(20,1000,10):
        v1,t1,v2,t2,v3,t3,v4,t4,v5,t5,Q_1,Q_5,Q_p = full_process(vy,ay)
        w_waste3, tz = Overall_results(t1, t2, t3, t4, t5, Q_1, Q_5, Q_p)

        H_array3.append(h/3.281)
        Pumped_array3.append(Q_p/2.204623)
        wasted_array3.append(w_waste3/2.204623)
        time_array3.append(tz)
        print(t5)
#plt.plot(H_array1,Pumped_array1,label = 'v0 = 0.5 m/sec')
#plt.plot(H_array2,Pumped_array2,label = 'v0 = 1.0 m/sec')
#plt.plot(H_array3,Pumped_array3,label = 'v0 = 1.5 m/sec')
#plt.xlabel('The static delivery head, h, (metres)')
#plt.ylabel('The quantity of water pumped per cycle, Q_p, (litres)')
#plt.legend()
#plt.grid()


#plt.plot(H_array1,wasted_array1
#,label = 'v0 = 0.5 m/sec')
#plt.plot(H_array2,wasted_array2,label = 'v0 = 1.0 m/sec')
#plt.plot(H_array3,wasted_array3,label = 'v0 = 1.5 m/sec')
#plt.xlabel('The static delivery head, h, (metres)')
#plt.ylabel('The quantity of water wasted per cycle, Q_w, (litres)')
#plt.legend()
#plt.grid()

#plt.plot(H_array1,time_array1,label = 'v0 = 0.5 m/sec')
#plt.plot(H_array2,time_array2,label = 'v0 = 1.0 m/sec')
#plt.plot(H_array3,time_array3,label = 'v0 = 1.5 m/sec')
#plt.xlabel('The static delivery head, h (m)')
#plt.ylabel('The total time of each ram pump cycle (seconds)')
#plt.grid()
#plt.show()


def water_wasted_testing():
    h_array = []
    wasted_array = []
    for h in range(50,500,1):
        v1,t1,v2,t2,v3,t3,v4,t4,v5,t5,Q_1,Q_5,Q_p = full_process()


        # Overall Process
        w_waste = Overall_results(t1,t2,t3,t4,t5,Q_1,Q_5,Q_p)
        h_array.append(h)
        wasted_array.append(w_waste)
    plt.plot(h_array,wasted_array)
    plt.show()
    return
def time_of_process():
    v0 = 3.1
    a0 = (X ** 2 - v0 ** 2) / (2 * L / j)

    h_array = []
    time_array = []
    for h in range(50,500,1):
        v1,t1,v2,t2,v3,t3,v4,t4,v5,t5,Q_1,Q_5,Q_p = full_process()


        # Overall Process
        w_waste,t = Overall_results(t1,t2,t3,t4,t5,Q_1,Q_5,Q_p)
        h_array.append(h)
        time_array.append(t)
    plt.plot(h_array,time_array)
    plt.show()

#second for graphs
vx = 3.281
ax = (X ** 2 - vx ** 2) / (2 * L / j)

# Testing involving varying h to get the quantity pumped.
#Second variable
H_array2 = []
Pumped_array2 = []
wasted_array2 = []
time_array2 = []
H_array1 = []
Pumped_array1 = []
wasted_array1 = []
time_array1 = []
for h in range(20,450,5):
    v1,t1,v2,t2,v3,t3,v4,t4,v5,t5,Q_1,Q_5,Q_p = full_process(vx,ax)
    w_waste2, ty = Overall_results(t1, t2, t3, t4, t5, Q_1, Q_5, Q_p)
    v1x, t1x, v2x, t2x, v3x, t3x, v4x, t4x, v5x, t5x, Q_1x, Q_5x, Q_px = full_processx(vx, ax)
    w_waste2x, tyx = Overall_results(t1x, t2x, t3x, t4x, t5x, Q_1x, Q_5x, Q_px)
    H_array2.append(h/3.281)
    Pumped_array2.append(Q_p/2.204623)
    wasted_array2.append(w_waste2/2.204623)
    time_array2.append(ty)
    H_array1.append(h/3.281)
    Pumped_array1.append(Q_px/2.204623)
    wasted_array1.append(w_waste2x/2.204623)
    time_array1.append(tyx)
#plt.plot(H_array1,Pumped_array1,label = 'v0 = 0.5 m/sec')
plt.plot(H_array2,time_array2,'r',label = 'Combined transient flow water hammer model')
#plt.plot(H_array1,Pumped_array1,'b',label = 'Initial model')
#plt.plot(H_array3,Pumped_array3,label = 'v0 = 1.5 m/sec')
plt.xlabel('The static delivery head, h, (metres)')
plt.ylabel('The total cycle time for the hydraulic ram pump (seconds)')
plt.legend()
plt.grid()


import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

# Each stage of the process
def stage_1(v0,S0,J,a_6):
    # initial condition

    t1 = np.sqrt((2 * S0) / J)
    # time points
    t = np.linspace(0, t1)
    def velocity(v,t):
        dvdt = a_6 - (a_6*t)/t1
        # plot results
        return dvdt

    # solve ODE
    v = odeint(velocity, v0, t)
    v1 = v[-1]
    # plot results
#    plt.plot(t, v,'r',label = 'Simulation Results from using parameter values in table 1')
#    plt.xlabel('time')
#    plt.ylabel('v(t)')
    Q_1 = w*A*(v0*t1 + a_6*t1**2/3)
    #plt.show()
    return v1, t1,Q_1
def stage_2(v1, H_instant, H_static, g, p, k, D,E,e,w,Av,A,Ev):
    a = np.sqrt((1)/(p/k)+ (D/E*e))
    delta_H = H_instant - H_static
    v2 = v1 - (delta_H *g)/a
    Z = (a*w*Av**2)/(A*g*Ev)
    t2 = Z*np.log(v1/v2)
    return v2, t2, delta_H, Z,a
def stage_3(t2,v1,L,m,a,h,H,Z):
    tb = 2*L/a - t2
    delta_V = (4*g*(h-H) + m*v1)/(4*a+m)
    vp = (v1 - delta_V)
    tp = Z* np.log((2*v1)/(2*v1 - vp))
    t3 = tb + tp
    N = 6
    #Water pumped in this section
    Q_p = w*Ac*(vp**2/20)
    #Q_p = w*A*(N*delta_V*tb +(N-1)*2*Z*delta_V +Z*vp - N**2*delta_V*tb -(N-1)**2*delta_V*t2 -(N-1)*v1*t2-(2*v1-vp)*tp)
    return vp , t3, delta_V,Q_p
def stage_4(delta_V,a,Z,v3):
    v4 = np.sqrt(((a*Z)/L)*(delta_V**2))
    t4 = -2*Z/(v4-v3) + 2*L/a
    return v4,t4
def stage_5(X,j,v4,v0):
    t5 = (L/(j*X))*np.log(((X+v0)/(X-v0))*((X-v4)/(X+v4)))
    v5 = v0
    #Water pumped in final section
    Q_5 = (w*A*L/j) * np.log((X**2 - v4**2)/(X**2 - v0**2))
    return v5 ,t5,Q_5

def Overall_results(t1,t2,t3,t4,t5,Q_1,Q_5,Q_p):
    t = t1, t1+ t2, t1+ t2 + t3,t1+t2+t3+t4,t1+t2+t3 + t4 + t5
    v = v1,v2,v3,v4,v5
#    plt.plot(t, v,'r')
#    plt.xlabel('Time (seconds)')
#    plt.ylabel('Velocity of flowing water at the waste valve, v(t)')
#    plt.legend(loc='upper right')
 #   plt.grid()

    W_wasted = Q_1 + Q_5
    W_pumped = Q_p
    Efficiency = W_pumped/(W_wasted+W_pumped) * 100
    #print(W_wasted)
    return W_wasted,t1+t2+t3+t4+t5
def full_process(v0,a0):
    #Stage 1 processes
    v1,t1,Q_1 = stage_1(v0,S0,J,a0)

    #Stage 2 process
    v2,t2,delta_H,Z,a = stage_2(v1, H_instant, H_static, g, p, k, D,E,e,w,Av,A,Ev)

    #Stage 3 process
    v3, t3, delta_V,Q_p = stage_3(t2,v1,L,m,a,h,H_static,Z)
    #Stage 4 process
    v4,t4 = stage_4(delta_V,a,Z,v3)

    #Stage 5 process

    v5,t5,Q_5 = stage_5(X,j,v4,v0)
    return v1,t1,v2,t2,v3,t3,v4,t4,v5,t5, Q_1,Q_5,Q_p


#Stage 1 parameters
S0 = 0.0161
J = 4.0

#Stage 2 Parameters
H_static = 9.2
H_instant = 65
g = 32.2
p = 62.4
k = 2200000000
D = 0.166667
E = 103
e = 0.153937
w = 62.43
Av = 0.1043
Ac = 0.018241
A = 0.0233
Ev = 3870000

#Stage 3 Parameters
L = 54.8
m = 817
#Stage 5 parameters
j = 15.5
X = np.sqrt(2*g*H_static/j)

def rankin_eff():

    #Rankine Efficiency
    v0 = 1.64
    a0 = (X ** 2 -v0 ** 2) / (2 * L / j)

    # Testing involving varying h to get the quantity pumped.

    H_array1 = []
    rankine_array = []
    for h in range(0, 300, 10):
        v1, t1, v2, t2, v3, t3, v4, t4, v5, t5, Q_1, Q_5, Q_p = full_process(v0, a0)
        w_waste1, tx = Overall_results(t1, t2, t3, t4, t5, Q_1, Q_5, Q_p)
        rankine_eff = (Q_p*(h-H_static))/(w_waste1*H_static)*10
        H_array1.append(h/3.281)
        rankine_array.append(rankine_eff)


    #Rankine Efficiency
    vx = 3.28
    ax = (X ** 2 - vx ** 2) / (2 * L / j)

    # Testing involving varying h to get the quantity pumped.

    H_array2 = []
    rankine_array2 = []
    for h in range(0, 600, 10):
        v1, t1, v2, t2, v3, t3, v4, t4, v5, t5, Q_1, Q_5, Q_p = full_process(vx, ax)
        w_waste2, tx = Overall_results(t1, t2, t3, t4, t5, Q_1, Q_5, Q_p)
        rankine_eff = (Q_p*(h-H_static))/(w_waste2*H_static)*10
        H_array2.append(h/3.281)
        rankine_array2.append(rankine_eff)



    #Rankine Efficiency
    vy = 4.92
    ay = (X ** 2 - vy ** 2) / (2 * L / j)

    # Testing involving varying h to get the quantity pumped.

    H_array3 = []
    rankine_array3 = []
    for h in range(0, 800, 10):
        v1, t1, v2, t2, v3, t3, v4, t4, v5, t5, Q_1, Q_5, Q_p = full_process(vy, ay)
        w_waste3, tx = Overall_results(t1, t2, t3, t4, t5, Q_1, Q_5, Q_p)
        rankine_eff = (Q_p*(h-H_static))/(w_waste3*H_static)*10
        H_array3.append(h/3.281)

        rankine_array3.append(rankine_eff)

    plt.plot(H_array1,rankine_array,label = 'v0 = 0.5 m/sec')
    plt.plot(H_array2,rankine_array2,label = 'v0 = 1.0 m/sec')
    plt.plot(H_array3,rankine_array3,label = 'v0 = 1.5 m/sec')
    plt.grid()
    plt.legend()
    plt.xlabel('The static delivery head, h, (metres)')
    plt.ylabel('The Rankine Efficiency of the hydraulic ram pump (%)')
    plt.show()
h = 65
#Tesing the initial velocity against the maximum water pumped
v_array1 = []
v_vary_waste_array1 = []
v_vary_pumped_array = []
def v_vary_waste():
    for i in range(20,80):
        v0 = 0.1*i
        a0 = (X ** 2 - v0 ** 2) / (2 * L / j)
        v1, t1, v2, t2, v3, t3, v4, t4, v5, t5, Q_1, Q_5, Q_p = full_process(v0,a0)

        # Overall Process
        w_waste1,tx = Overall_results(t1, t2, t3, t4, t5, Q_1, Q_5, Q_p)
        v_array1.append(v0)
        v_vary_pumped_array.append(Q_p/2.204623)
        v_vary_waste_array1.append(w_waste1/2.204623)
        print('the value of v0 is', v0, 'and water wasted is ',w_waste1 )
    plt.plot(v_array1,v_vary_waste_array1,label = 'Water wasted per cycle,Q_w')
    plt.plot(v_array1,v_vary_pumped_array, label = 'Water pumped per cycle,Q_p')
    plt.xlabel('Velocity threshold of the waste valve before movement occurs, v0')
    plt.ylabel('Water wasted/pumped per cycle (litres)')
    plt.legend()
    plt.grid()
    plt.show()
def pumped_wasted_etc():
    v0 = 1.64
    a0 = (X ** 2 - v0 ** 2) / (2 * L / j)
    time5 = []
    # Testing involving varying h to get the quantity pumped.

    H_array1 = []
    Pumped_array1 = []
    wasted_array1 = []
    time_array1 = []
    for h in range(20,350,10):
        v1,t1,v2,t2,v3,t3,v4,t4,v5,t5,Q_1,Q_5,Q_p = full_process(v0,a0)
        w_waste1, tx = Overall_results(t1, t2, t3, t4, t5, Q_1, Q_5, Q_p)
        H_array1.append(h/3.281)
        Pumped_array1.append(Q_p/2.204623)
        wasted_array1.append(w_waste1/2.204623)
        time_array1.append(tx)

    #second for graphs
    vx = 3.28
    ax = (X ** 2 - v0 ** 2) / (2 * L / j)

    # Testing involving varying h to get the quantity pumped.
    #Second variable
    H_array2 = []
    Pumped_array2 = []
    wasted_array2 = []
    time_array2 = []
    for h in range(20,700,10):
        v1,t1,v2,t2,v3,t3,v4,t4,v5,t5,Q_1,Q_5,Q_p = full_process(vx,ax)
        w_waste2, ty = Overall_results(t1, t2, t3, t4, t5, Q_1, Q_5, Q_p)
        H_array2.append(h/3.281)
        Pumped_array2.append(Q_p/2.204623)
        wasted_array2.append(w_waste2/2.204623)
        time_array2.append(ty)
    vy = 4.92
    ay = (X ** 2 - v0 ** 2) / (2 * L / j)

        # Testing involving varying h to get the quantity pumped.
    #Third variable
    H_array3 = []
    Pumped_array3 = []
    wasted_array3 = []
    time_array3 = []
    for h in range(20,1000,10):
        v1,t1,v2,t2,v3,t3,v4,t4,v5,t5,Q_1,Q_5,Q_p = full_process(vy,ay)
        w_waste3, tz = Overall_results(t1, t2, t3, t4, t5, Q_1, Q_5, Q_p)

        H_array3.append(h/3.281)
        Pumped_array3.append(Q_p/2.204623)
        wasted_array3.append(w_waste3/2.204623)
        time_array3.append(tz)
        print(t5)

#second for graphs
vx = 3.28
ax = (X ** 2 - vx ** 2) / (2 * L / j)

# Testing involving varying h to get the quantity pumped.
#Second variable
H_array2 = []
Pumped_array2 = []
wasted_array2 = []
time_array2 = []
for h in range(20,650,10):
    v1,t1,v2,t2,v3,t3,v4,t4,v5,t5,Q_1,Q_5,Q_p = full_process(vx,ax)
    w_waste2, ty = Overall_results(t1, t2, t3, t4, t5, Q_1, Q_5, Q_p)
    H_array2.append(h/3.281)
    Pumped_array2.append(Q_p/2.204623)
    wasted_array2.append(w_waste2/2.204623)
    time_array2.append(ty)

#plt.plot(H_array1,Pumped_array1,label = 'v0 = 0.5 m/sec')
plt.plot(H_array2,time_array2,label = 'Initial lumped mass model')
#plt.plot(H_array3,Pumped_array3,label = 'v0 = 1.5 m/sec')
plt.xlabel('The static delivery head, h, (metres)')
#plt.ylabel('The quantity of water pumped per cycle, Q_p, (litres)')
plt.legend()
plt.grid()
plt.show()