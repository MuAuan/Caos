from scipy.integrate import odeint, simps
import numpy as np
import matplotlib.pyplot as plt

"""
def pend(y, t, b, c):
    theta, omega = y
    dydt = [omega, -b*omega - c*np.sin(theta)]
    return dydt

b = 0.025
c = 5.0
y0 = [np.pi - 0.1, 0.0]
t = np.linspace(0, 100, 1001)
sol = odeint(pend, y0, t, args=(b, c))

x, p = sol.T[0], sol.T[1]
plt.plot(x, p, ".", markersize=4)
plt.pause(3)
plt.savefig('plot_xp.png', dpi=60)
plt.close()

plt.plot(t, sol[:, 0], 'b', label='theta(t)')
plt.xlabel('t')
plt.grid()
plt.legend(loc='best')
plt.pause(3)
plt.plot(t, sol[:, 1], 'g', label='omega(t)')
plt.legend(loc='best')
plt.pause(3)
plt.savefig('plot_thomegavst.png', dpi=60)
plt.close()
"""

def dualPendulum(y, t, L, omegaJ,M,M1):
    theta1, omega1,theta2, omega2 = y
    dydt = [omega1, 
            (-omegaJ*(np.sin(theta1)+M*np.sin(theta1-2*theta2))-2*np.sin(theta1-theta2)*M*(omega2*omega2*L+omega1*omega1*np.cos(theta1-theta2)))/(1*(1-M*np.cos(2*(theta1-theta2)))), 
            omega2, 
            (2*np.sin(theta1-theta2)*(omega1*omega1*M1+omegaJ*M1*np.cos(theta1)+omega2*omega2*L*M1*np.cos(theta1-theta2)))/(L*(1 - M*np.cos(2*(theta1-theta2))))]
    return dydt

def doublePendulum(y, t, L, omegaJ,M):
    theta1, omega1,theta2, omega2 = y
    theta12=theta1 -theta2
    dydt = [omega1,
            (omegaJ*(-np.sin(theta1)+M*np.cos(theta12)*np.sin(theta2))-M*(omega1*omega1*np.cos(theta12)+L*omega2*omega2)*np.sin(theta12))/(1-M*np.cos(theta12)*np.cos(theta12)),
            omega2,
            (omegaJ*(-np.sin(theta2)+np.cos(theta12)*np.sin(theta1))+(M*L*omega2*omega2*np.cos(theta12)+omega1*omega1)*np.sin(theta12))/(L-M*L*np.cos(theta12)*np.cos(theta12))]
    return dydt
  
def doublePendulum2(y, t, L, omegaJ,M,M1):
    theta1, omega1,theta2, omega2 = y
    theta12=theta1 -theta2
    dydt = [omega1,
            (omegaJ*(-M1*np.sin(theta1)+M*np.cos(theta12)*np.sin(theta2))-M*(2*omega1*omega1*np.cos(theta12)+L*omega2*omega2)*np.sin(theta12))/(1-2*M*np.cos(theta12)*np.cos(theta12)),
            omega2,
            (omegaJ*(-np.sin(theta2)+2*M1*np.cos(theta12)*np.sin(theta1))+(2*M*L*omega2*omega2*np.cos(theta12)+2*omega1*omega1)*np.sin(theta12))/(L-2*M*L*np.cos(theta12)*np.cos(theta12))]
    return dydt    


L1=7
L2=1
g=9.8
m1=100
m2=20
y01=0.5
for i in range(0,10,1):
    y01=(2.1115*0.0001)*np.exp(-i)  #3.14*np.exp(i)/np.exp(10)
    M=m2/(2*m1+m2)
    M1=(m1+m2)/(2*m1+m2)
    omegaJ=g/L1
    L=L2/L1
    y0 = [y01, 0.0,0, 0.0]
    t = np.linspace(0, 1000, 4001)
    sol = odeint(dualPendulum, y0, t, args=(L, omegaJ,M,M1))

    x1, p1,x2,p2 = sol.T[0],sol.T[1], sol.T[2], sol.T[3]
    plt.plot(x1, p1, ".", markersize=4)
    plt.pause(3)
    plt.savefig('plot_xp1dual_L1{}L2{}m1{}m2{}y01{}.png'.format(L1, L2,m1,m2,y01), dpi=180)
    plt.close()
    plt.plot(x2, p2, ".", markersize=4)
    plt.pause(3)
    plt.savefig('plot_xp2dual_L1{}L2{}m1{}m2{}y01{}.png'.format(L1, L2,m1,m2,y01),dpi=180)
    plt.close()
    """
    plt.plot(t, sol[:, 0], 'b', label='theta1(t)')
    plt.xlabel('t')
    plt.xlim(0,100)
    #plt.ylim(-0.5,0.5)
    plt.grid()
    plt.legend(loc='best')
    plt.pause(3)
    plt.plot(t, sol[:, 2], 'g', label='theta2(t)')
    plt.legend(loc='best')
    plt.pause(3)
    plt.savefig('plot_theta2vstdual_L1{}L2{}m1{}m2{}y01{}.png'.format(L1, L2,m1,m2,y01), dpi=180)
    plt.close()

    plt.plot(t, sol[:, 1], 'b', label='omega1(t)')
    plt.xlabel('t')
    plt.xlim(0,100)
    #plt.ylim(-1,1)
    plt.grid()
    plt.legend(loc='best')
    plt.pause(3)
    plt.plot(t, sol[:, 3], 'g', label='omega2(t)')
    plt.legend(loc='best')
    plt.pause(3)
    plt.savefig('plot_omega2vstdual_L1{}L2{}m1{}m2{}y01{}.png'.format(L1, L2,m1,m2,y01), dpi=180)
    plt.close()
    """
    
    M=m2/(m1+m2)
    L=L2/L1
    omegaJ=g/L1  
    #y01=0.0020
    y0 = [y01, 0.0,0, 0.0]
    t = np.linspace(0, 1000, 4001)
    #sol = odeint(dualPendulum, y0, t, args=(L1, L2,g,m1,m2))
    #sol = odeint(doublePendulum, y0, t, args=(L1, L2,g,m1,m2))
    sol = odeint(doublePendulum, y0, t, args=(L, omegaJ,M))

    x1, p1,x2,p2 = sol.T[0],sol.T[1], sol.T[2], sol.T[3]
    plt.plot(x1, p1, ".", markersize=4)
    plt.pause(3)
    #plt.savefig('plot_xp1pen_L{}M{}omegaJ{}y01{}n.png'.format(L, omegaJ,M,y01), dpi=180)
    plt.savefig('plot_xp1pen_L1{}L2{}m1{}m2{}y01{}n.png'.format(L1, L2,m1,m2,y01), dpi=180)
    plt.close()
    plt.plot(x2, p2, ".", markersize=4)
    plt.pause(3)
    plt.savefig('plot_xp2pen_L{}M{}omegaJ{}y01{}n.png'.format(L, omegaJ,M,y01),dpi=180)
    plt.close()
    """
    plt.plot(t, sol[:, 0], 'b', label='theta1(t)')
    plt.xlabel('t')
    plt.xlim(0,100)
    #plt.ylim(-0.5,0.5)
    plt.grid()
    plt.legend(loc='best')
    plt.pause(3)
    plt.plot(t, sol[:, 2], 'g', label='theta2(t)')
    plt.legend(loc='best')
    plt.pause(3)
    plt.savefig('plot_theta2vstpen_L1{}L2{}m1{}m2{}y01{}n.png'.format(L1, L2,m1,m2,y01), dpi=180)
    plt.close()

    plt.plot(t, sol[:, 1], 'b', label='omega1(t)')
    plt.xlabel('t')
    plt.xlim(0,100)
    #plt.ylim(-1,1)
    plt.grid()
    plt.legend(loc='best')
    plt.pause(3)
    plt.plot(t, sol[:, 3], 'g', label='omega2(t)')
    plt.legend(loc='best')
    plt.pause(3)
    plt.savefig('plot_omega2vstpen_L1{}L2{}m1{}m2{}y01{}n.png'.format(L1, L2,m1,m2,y01), dpi=180)
    plt.close()
    """

    L1=L1
    L2=L2
    g=9.8
    m1=m1
    m2=m2
    M=2*m2/(m1+4*m2)
    M1=(m1+2*m2)/(m1+4*m2)
    L=L2/L1
    omegaJ=g/L1  
    y01=y01
    y0 = [y01, 0.0,0, 0.0]
    t = np.linspace(0, 1000, 4001)
    #sol = odeint(dualPendulum, y0, t, args=(L1, L2,g,m1,m2))
    #sol = odeint(doublePendulum, y0, t, args=(L1, L2,g,m1,m2))
    sol = odeint(doublePendulum2, y0, t, args=(L, omegaJ,M,M1))

    x1, p1,x2,p2 = sol.T[0],sol.T[1], sol.T[2], sol.T[3]
    plt.plot(x1, p1, ".", markersize=4)
    plt.pause(3)
    #plt.savefig('plot_xp1Pen2_L{}M{}omegaJ{}y01{}n.png'.format(L, omegaJ,M,M1,y01), dpi=180)
    plt.savefig('plot_xp1Pen2_L1{}L2{}m1{}m2{}y01{}n.png'.format(L1, L2,m1,m2,y01), dpi=180)
    plt.close()
    plt.plot(x2, p2, ".", markersize=4)
    plt.pause(3)
    plt.savefig('plot_xp2Pen2_L{}M{}omegaJ{}y01{}n.png'.format(L, omegaJ,M,y01),dpi=180)
    plt.close()
    """
    plt.plot(t, sol[:, 0], 'b', label='theta1(t)')
    plt.xlabel('t')
    plt.xlim(0,100)
    #plt.ylim(-0.5,0.5)
    plt.grid()
    plt.legend(loc='best')
    plt.pause(3)
    plt.plot(t, sol[:, 2], 'g', label='theta2(t)')
    plt.legend(loc='best')
    plt.pause(3)
    plt.savefig('plot_theta2vstPen2_L{}M{}omegaJ{}y01{}n.png'.format(L, omegaJ,M,y01), dpi=180)
    plt.close()

    plt.plot(t, sol[:, 1], 'b', label='omega1(t)')
    plt.xlabel('t')
    plt.xlim(0,100)
    #plt.ylim(-1,1)
    plt.grid()
    plt.legend(loc='best')
    plt.pause(3)
    plt.plot(t, sol[:, 3], 'g', label='omega2(t)')
    plt.legend(loc='best')
    plt.pause(3)
    plt.savefig('plot_omega2vstPen2_L{}M{}omegaJ{}y01{}n.png'.format(L, omegaJ,M,y01), dpi=180)
    plt.close()
    """