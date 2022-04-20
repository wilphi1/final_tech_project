import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate
#Steady flow rate for the interior of each section of a pipe;
def inner_node_steady(H0, V0, dt, g,L,D,a):
    HP = np.zeros(mx)
    VP = np.zeros(mx)
    # property of current pipe
    f = 4.72   # unitless      # m   # m/s
    #A = np.pi * D**2. / 4.  # m^2
    theta = 0.23
    ga = g/1000
    for i in range(1,L-1):
        V1 = V0[i-1]; H1 = H0[i-1]
        V2 = V0[i+1]; H2 = H0[i+1]
        C = np.zeros((2,2), dtype=np.float64)

        J1 = f*dt/2./D*V1*abs(V1)
        C[0,0] = V1 + ga*H1 - J1 + ga* dt *V1*theta
        C[0,1] = ga

        J2 = f*dt/2./D*V2*abs(V2)
        C[1,0] = -V2+ ga*H2 + J2 - ga* dt *V2*theta
        C[1,1] = ga

        HP[i] = (C[0,0] + C[1,0]) / (C[0,1] + C[1,1])
        VP[i] = np.float64(-C[1,0]+ C[1,1]*HP[i])

    return HP[1:-1], VP[1:-1]
#Each pipe section of the process:
def pipe_flow_supply_resevoir(H, V,HR,VR, H0, V0, dt,g, LL,LR,a):
    D1 = D2 = 0.00968  #
    H[1:-1], V[1:-1] = inner_node_steady(H0, V0, dt, g,LL,D1,a)
    f1 =f2 = 4.72  # unitless m
    a1 = a2 = a  # m/s
    g = 9.81
    theta1 = theta2 = 0.23
    #Resevoir side
    # Pipe start (outer boundayr conditions)
    V2 = V0[1]
    H2 = H0[1]

    #Resevoir side boundary conditions:
    H[0], V[0] = rev_end(H2,V2,H0[0],f1,D1,g,a1,dt)

    s1 = 1
    s2 = -1
    # Pipe end  (inner boundary conditions)
    V1 = V0[-2]
    H1 = H0[-2]     # upstream node
    KL_inv = 0.000000215
    H[-1],V[-1] = valve_node(KL_inv,H,V,HR,VR,dt,g,s1,s2,LL,LR,f1,f2,D1,D2,a1,a2,theta1,theta2)
    #V2 = V20;     H2 = H20         # downstream nodes
    return H,V
def pipe_flow_delivery_resevoir(H, V,HR,VR, H0, V0, dt,g, LL,LR,a):
    D1 = D2 = 0.0026  # m
    H[1:-1], V[1:-1] = inner_node_steady(H0, V0, dt, g,LR,D1,a)
    f1 =f2 = 4.72 # unitless
    a1 = a2 = a  # m/s
    g = 9.81
    theta1 = theta2 = 0.2
    #Resevoir side
    # Pipe start (outer boundayr conditions)
    s1 = 1
    s2 = -1
    # Pipe end  (inner boundary conditions)
    V1 = VR[-2]
    H1 = HR[-2]     # upstream node
    a3 = 1300
    a4 = 1300
    #Resevoir side boundary conditions:
    H[0],V[0],Qs,vs = air_chamber(H, V, HR, VR, dt, 0, s1, s2,LL,LR,f1,f2,D1,D2,a1,a2,theta1,theta2)
    H[-1],V[-1] = rev_end(H1,V1,H0[0],f1,D1,g,a1,dt)
    #V2 = V20;     H2 = H20         # downstream nodes
    return H,V
def supply_to_waste_valve(H, V,HL,VL, H0, V0, dt,g, LL,LR,a):
    f1 = f2 = 4.72  # unitless
    D1 = D2 = 0.00968 # m
    a1= a2 = a # m/s#
    theta1=theta2=0.2
    g = 9.81
    s1 = -1
    s2 = 1
    KL_inv = 251
    V1 = V0[-2]
    H1 = H0[-2]
    H[-1],V[-1] = waste_valve_end(H1,V1,H[-1],V[-1],f1,D1,g,a1,dt)
    # Pipe end  (inner boundary conditions)
    # upstream node

    H[0],V[0] =valve_node(KL_inv,HL,VL,H,V,dt,g,s1,s2,LL,LR,f1,f2,D1,D2,a1,a2,theta1,theta2)
    # V2 = V20;     H2 = H20         # downstream nodes

    return H, V
def ram_entrance_to_check_valve(H, V,HL,VL, H0, V0, dt,g, LL,LR,a):
    f1 = f2 = 4.72  # unitless
    D1 = D2 = 0.00968 # m
    a1= a2 = a # m/s#
    theta1=theta2=0.2
    g = 9.81
    s1 = -1
    s2 = 1
    KL_inv = 251
    #Check valve inner boundary junction:
    H[0], V[0] = valve_node(KL_inv, HL, VL, H, V, dt, g, s1, s2, LL, LR, f1, f2, D1, D2, a1, a2, theta1, theta2)

    #Air chamber boundary:
    V1 = V0[-2]
    H1 = H0[-2]
    H[-1],V[-1],Qs,vs = air_chamber(HL, VL, H, V, dt, LR, s1, s2,LL,LR,f1,f2,D1,D2,a1,a2,theta1,theta2)
    # Pipe end  (inner boundary conditions)
    # upstream node

    #Now other end of the pipe, boundary condition of air chamber node.

    return H,V


#Boundary Conditions:
def valve_node(KL_inv, H1, V1,H2,V2, dt, g, s1, s2,L_left,L_right,f1,f2,D1,D2,a1,a2,theta1,theta2):
    nn = mx
    A,A1, C1, C2 = cal_Cs(H1,V1,H2,V2,s1,s2,dt,L_left,L_right,f1,f2,D1,D2,a1,a2,theta1,theta2)
    # parameters of the quadratic polynomial
    aq = 1
    bq = 2*g*KL_inv* (1/C1[0,1] + 1/C2[0,1])
    cq = 2*g*KL_inv* (C2[0,0]/C2[0,1] - C1[0,0]/C1[0,1])

    # solve the quadratic equation
    delta = bq**2 - 4*aq*cq

    if delta >= 0:
        VP = (-bq + np.sqrt(delta))/(2*aq)
    elif delta > -1.0e-7 and delta <0 :
        VP = (-bq)/(2*aq)
    else:
        VP = (-bq)/(2*aq)

    if VP >=0 : # positive flow
        if nn == 0:  # pipe start
            VP = VP
            HP = (C2[0,0] + VP) / C2[0,1]
        else:        # pipe end
            VP = VP
            HP = (C1[0,0] - VP) / C1[0,1]

    else : # reverse flow
        # reconstruct the quadratic equation
        # parameters of the quadratic polynomial
        aq = 1
        bq = 2*g*KL_inv* (-A/C2[0,1]-1/C1[0,1])
        cq = 2*g*KL_inv* (-C2[0,0]/C2[0,1]+C1[0,0]/C1[0,1])

        # solve the quadratic equation
        delta = bq**2 - 4*aq*cq

        if delta >= 0:
            VP = (-bq - np.sqrt(delta))/(2*aq)
        elif delta > -1.0e-7 and delta <0 :
            VP = (-bq)/(2*aq)
        else:
            VP = (-bq)/(2*aq)


        if nn == 0:  # pipe start
            VP = VP
            HP = (C2[0,0] + VP ) / C2[0,1]
        else:        # pipe end
            VP = VP
            HP = (C1[0,0] - VP) / C1[0,1]
    return HP, VP
def waste_valve_end(H1,V1,H,V,f,D,g,a,dt):
    # Resevoir side
    J = cal_friction(f,D,V,dt)
    HP = H1 + a / g * (V - V1) + a / g * J
    VP = V1 - g/a*(H-H1) - J
    return HP, VP
def rev_end(H2,V2,H,f,D,g,a,dt):
    J = cal_friction(f, D, V2, dt)
    VP = V2 - g/a*(H - H2) - J
    HP = H
    return HP, VP
def air_chamber(H1, V1, H2, V2, dt, nn, s1, s2,L_left,L_right,f1,f2,D1,D2,a1,a2,theta1,theta2):
    A1, A2, C1, C2 = cal_Cs(H1, V1,H2,V2, s1, s2, dt,L_left,L_right,f1,f2,D1,D2,a1,a2,theta1,theta2)
    # parameters
    Hb = 10.5 # barometric pressure head
    m = 1.2
    tank = [0.283,0.19,0.287058,0.095,0.5]
    As, ht, C, z, Qs = tank  # tank properties and results at last time step
    at = 2.* As/dt
    Va = (ht-z)*As  # air volume at last time step
    Cor = 0
    a = np.dot(C1[-1,0], A1) + np.dot(C2[-1,0],A2)
    b = np.dot(C1[-1,1], A1) + np.dot(C2[-1,1],A2)

    def tank_flow(QPs, Qs, a, b, As, C, z, at, Va, Cor, m, Hb):
        return (((a-QPs)/b + Hb - z - (Qs+QPs)/at - Cor*QPs*np.abs(QPs))
                 * (Va- (Qs+QPs)*As/at)**m - C)

    def tank_flow_prime(QPs, Qs, a, b, As, C, z, at, Va, Cor, m, Hb):
        p1 = (-m*As/at * (Va- (Qs+QPs)*As/at)**(m-1)*
            ((a-QPs)/b + Hb - z - (Qs+QPs)/at - Cor*QPs*np.abs(QPs)))
        p2 = (-1/b -1/at - Cor*2.*QPs*np.sign(QPs))* (Va- (Qs+QPs)*As/at)**m
        return p1+p2

    # solve nonlinear equation for tank flow at this time step
    from scipy import optimize
    QPs = optimize.newton(tank_flow,Qs,fprime=tank_flow_prime, args=(Qs, a, b, As, C, z, at, Va, Cor, m, Hb))
        #optimize.newton(tank_flow, Qs, fprime=tank_flow_prime,

        #    args=(Qs, a, b, As, ht, C, z, at, Va, Cor, m, Hb),
       #     tol=1e-10)

    zp = z + (Qs+QPs)/at

    HP = (a - QPs)/b
    VP2 = - C2[:,0] + C2[:,1]*HP
    VP1 =  C1[:,0] - C1[:,1]*HP
    HP = np.float64(HP)
    if nn == 0:  # pipe start
        VP = np.float64(VP2)
    else:        # pipe end
        VP = np.float64(VP1)
    VP = VP[0]
    return HP, VP, QPs, zp

#Calculations
def cal_friction(f, D, V, dt):
    Ju = 0
    Js = f * dt/ 2. / D * V * abs(V)  # steady frictio
    return Ju + Js
def cal_Cs(H1, V1,H2,V2, s1, s2, dt,L_left,L_right,f1,f2,D1,D2,a1,a2,theta1,theta2):
    g = 9.81
    # property of left adjacent pipe
    A1 = np.pi * D1**2. / 4  # m^2
    A2 = np.pi * D2**2. / 4
    C1 = np.zeros((2,2))
    C2 = np.zeros((2, 2))
    # J = f1[i]*dt/2./D1[i]*V1[i]*abs(V1[i])
    for i in range(0,L_left):
        J = cal_friction(f1, D1, V1[i], dt)

        C1[:,0] = s1*V1[i] + g/a1*H1[i] - s1*J + g/a1* dt *V1[i]*theta1
        C1[:,1] = g/a1
    for i in range(0,L_right):
        J = cal_friction(f2, D2, V1[i], dt)
        C2[:,0] = s2*V2[i] + g/a2*H2[i] - s2 *J + g/a2* dt *V2[i]*theta2
        C2[:,1] = g/a1
    return A1,A2, C1, C2
# Set numerical parameters
mx = 17
my = 2
mz = 2 # number of gridpoints in space
mw = 17
mt = 70  # number of gridpoints in time
dt = 1/200
L1 = mx
L2 = 2
L3 = 2
L4 = mw
# print("deltax=", deltax)
# print("deltat=", deltat)
# print("lambda=", lmbda)

# results from last time step
H_supply = [0] * mx
H_waste = [0] * my
H_check = [0]* mz
H_deliver = [0]*mw
V_supply = [0] * mx
V_waste = [0] * my
V_check = [0] * mz
V_deliver = [0]*mw
# results at current time step
HN_supply = [0] * mx
HN_waste = [0] * my
HN_check = [0] *mz
HN_deliver = [0]*mw
VN_supply = [0] * mx
VN_waste = [0] * my
VN_check = [0] * mz
VN_deliver = [0] *mw

Hb = 10.3  # barometric head
g = 9.81
#Initial Conditions
for i in range(0,mx):
    H_supply[i] = 2.8014 - 2.80148*(i/mx)
    V_supply[i] = -0.8488
for i in range(0,my):
    #Waste Valve to supply pipe
    H_waste[i] = 0
    V_waste[i] = 0
    H_check[i] = 0
    V_check[i] = 0
for i in range(0,mw):
    H_deliver[i] = 65*(i/mx)
    V_deliver[i] = 0
V_at_resevoir = []
V_middle_of_supply = []
V_end_of_supply = []
V_waste_valve = [0]
V_check_valve = [0]
V_delivery_pipe = [0]

H_at_resevoir = []
H_middle_of_supply = []
H_end_of_supply = []
H_waste_valve = [0]
H_check_valve = [0]
H_delivery_pipe = []
#Start Calculation
for ts in range(1,mt):
    for pn in range(0,mx):
        HN_supply[pn] = np.zeros_like(H_supply[pn])
        VN_supply[pn] = np.zeros_like(V_supply[pn])
    for py in range(0,my):
        HN_waste[py] = np.zeros_like(H_waste[py])
        VN_waste[py] = np.zeros_like(V_waste[py])
        HN_check[py] = np.zeros_like(H_check[py])
        VN_check[py] = np.zeros_like(V_check[py])
    for pw in range(0,mw):
        HN_deliver[pw] = np.zeros_like(H_deliver[pw])
        VN_deliver[pw] = np.zeros_like(V_deliver[pw])

    a = 1383
    for pn in range(0,mx):
        HN_supply,VN_supply = pipe_flow_supply_resevoir(HN_supply, VN_supply,H_waste,V_waste, H_supply, V_supply,dt,g,L1,L2,a)
    for pw in range(0,mw):
        HN_deliver,VN_deliver = pipe_flow_delivery_resevoir(HN_deliver, VN_deliver,H_check,V_check, H_deliver, V_deliver,dt,g,L4,L3,a)
    for py in range(0,my):
        HN_waste,VN_waste = supply_to_waste_valve(HN_waste, VN_waste,H_supply,V_supply, H_waste, V_waste,dt,g,L1,L2,a)
        HN_check,VN_check = ram_entrance_to_check_valve(HN_check, VN_check,H_waste,V_waste, H_waste, V_waste,dt,g,L2,L3,a)
        # record results
    # march in time
    for pn in range(0,mx):
        H_supply[pn] = HN_supply[pn]
        V_supply[pn] = VN_supply[pn]
    for py in range(0,my):
        H_waste[py] = HN_waste[py]
        V_waste[py] = VN_waste[py]
        H_check[py] = HN_check[py]
        V_check[py] = VN_check[py]
    for pw in range(0,mw):
        H_deliver[pw] = HN_deliver[pw]
        V_deliver[pw] = VN_deliver[pw]

#    print('The head of the supply pipe is', H_supply)
#    print('The velocity of the supply pipe is ', V_supply)
#    print('The head of the waste pipe is', H_waste)
#    print('The velocity of the waste pipe is ', V_waste)
#    print('The head of the check pipe is', H_check)
#    print('The velocity of the check pipe is ', V_check)
#    print('The head of the delivery pipe is ',H_deliver)
#    print('The velocity of the delivery pipe is  ', V_deliver)
    #Some kind of plotting function to record results
#    plt.plot(mx,H_waste)
    V_at_resevoir.append(V_supply[1])
    V_middle_of_supply.append(V_supply[10])
    V_end_of_supply.append(V_supply[-1])
    V_waste_valve.append(V_waste[-1])
    V_check_valve.append(V_check[-1])
    V_delivery_pipe.append(V_deliver[0]*0.1*0.026)
    H_at_resevoir.append(H_supply[1])
    H_middle_of_supply.append(H_supply[10])

    H_end_of_supply.append(H_supply[-1])
    H_waste_valve.append(H_waste[-1])
    H_check_valve.append(H_check[0]/2.1)
    H_delivery_pipe.append(H_deliver[0])
    #Head_at_waste_valve.append()
    print(np.max(H_at_resevoir))
    #    plt.plot(mx,H_check)

#print(Head_at_resevoir)
#tt = np.linspace(0,mt,mt-1)
tt = np.linspace(0,mt*dt, mt-1)
#Head close to the resevoir:
plt.plot(tt,H_at_resevoir,'g',label ='Pressure head close to the supply resevoir after water hammer')
plt.ylabel('Pressure head (m)')
plt.xlabel('Time (seconds)')
plt.legend()
plt.show()
plt.plot(tt,V_at_resevoir,'g',label ='Flow velocity close to the supply resevoir after water hammer')
plt.ylabel('Flow Velocity (m/s)')
plt.xlabel('Time (seconds)')
plt.legend()
plt.show()
#Head middle of the supply pipe
#plt.plot(tt,H_middle_of_supply,label='Head at supply resevoir')
#plt.legend()
#plt.show()
#Head end of the supply pipe
#plt.plot(tt,H_end_of_supply,label='Head at supply resevoir')
#plt.legend()
#plt.show()
tx = np.linspace(0,mt*dt, mt)
#Head at the waste valve
plt.plot(tx,V_waste_valve,'r',label ='Flow velocity at the waste valve after water hammer')
plt.ylabel('Flow Velocity (m/s)')
plt.xlabel('Time (seconds)')
plt.legend()
plt.show()
plt.plot(tx,V_check_valve,'r',label ='Flow velocity at the delivery valve after water hammer')
plt.ylabel('Flow Velocity (m/s)')
plt.xlabel('Time (seconds)')
plt.legend()
plt.show()

plt.plot(tx,H_waste_valve,'r',label ='Pressure head at the waste valve after water hammer')
plt.ylabel('Pressure head (m)')
plt.xlabel('Time (seconds)')
plt.legend()
#Head before the delivery valve
plt.plot(tx,H_check_valve,'b',label ='Pressure head at the delivery valve after water hammer')
plt.ylabel('Pressure head (m)')
plt.xlabel('Time (seconds)')
plt.legend()
plt.show()
#Velocity at the Delivery pipe
plt.plot(tx,V_delivery_pipe,'y',label ='Flow Velocity at the beginning of the delivery pipe after water hammer')
plt.ylabel('Flow Velocity (m/s)')
plt.xlabel('Time (seconds)')
plt.legend()
plt.show()

#tt = np.linspace(0,mt,mt-1)



