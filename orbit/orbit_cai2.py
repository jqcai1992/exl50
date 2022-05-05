# -*- coding: utf-8 -*-
"""
Created on Wed Dec  9 10:25:11 2020

@author: 83022
"""

from math import *
import numpy as np
import json
from scipy import interpolate
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation
from tqdm.std import trange


shot_no=6935
time=4.45
file_name = 'EXL%06.d_%05.fms.json'%(shot_no,time*1000)
with open(file_name) as file:
    exl=json.load(file)
rgrid=exl['r_grid']
zgrid=exl['z_grid']
Br=exl['Br']
Bt=exl['Bt']
Bz=exl['Bz']
R0=0.6
K0=5e3  #6keV的电子
pitch_angle=85 #65
direction= -1
    
    
Br_func = interpolate.interp2d(rgrid, zgrid,Br, kind='linear')
Bt_func = interpolate.interp2d(rgrid, zgrid,Bt, kind='linear')
Bz_func = interpolate.interp2d(rgrid, zgrid,Bz, kind='linear')
rmin=min(rgrid)
rmax=max(rgrid)
zmin=min(zgrid)
zmax=max(zgrid)
q = 1.602176565e-19; # Elementary charge (Coulomb)
m = 1.672621777e-27; # Proton mass (kg)
m_el = 9.10938291e-31; # Electron mass (kg)
c = 299792458; # speed of light (m/s)


def getb(x,y,z):
    fR=np.sqrt(x**2+y**2)
    fz=z
    if rmin<=fR<=rmax and zmin<=fz<=zmax:
        bR=Br_func(fR,fz)
        bt=Bt_func(fR,fz)
        bz=Bz_func(fR,fz)
            
        by=bR*y/fR+bt*x/fR
        bx=x/fR*bR-bt*y/fR
    else:
        bx=[0]
        by=[0]
        bz=[0]
        
    return bx[0],by[0],bz[0]  
        


def getb1(x,y,z):
    r=sqrt((sqrt(x**2+y**2)-R0)**2+z**2)+1e-10
    R=sqrt(x**2+y**2)
    q=1+r**2
    bt=B0*R0/R
    bp=bt*r/q/R0
    bx=-y*bt/R-bp*z*x/r/R
    by=bt*x/R-bp*z/r*y/R
    bz=bp*(R-R0)/r
    return bx,by,bz

def Lonz(y):
    dy=np.zeros(len(y))
    dy[0]=y[3]
    dy[1]=y[4]
    dy[2]=y[5]
    bx,by,bz=getb(y[0],y[1],y[2])
    dy[3]=q/m*(y[4]*bz-y[5]*by)
    dy[4]=q/m*(y[5]*bx-y[3]*bz)
    dy[5]=q/m*(y[3]*by-y[4]*bx)
    return dy

def rk4(func,t_seq=[],y0=[]):
    y=np.zeros([len(t_seq)+1,len(y0)])
    y[0]=y0
    for i in trange(len(t_seq)):              
        k1 = func(y[i])
        k2 = func(y[i]+k1*tau/2)
        k3 = func(y[i]+k2*tau/2)
        k4 = func(y[i]+k3*tau) 
        y[i+1]=y[i]+1/6*tau*(k1+2*k2+2*k3+k4)
    return y

def animate(i):
    line.set_xdata(xx[:100*i])
    line.set_ydata(yy[:100*i])
    line.set_3d_properties(zz[:100*i])
    point.set_xdata(xx[100*i])
    point.set_ydata(yy[100*i])
    point.set_3d_properties(zz[100*i])
    


v=sqrt(K0*q/m)
vpara=direction*v*cos(pitch_angle*pi/180)
vpen=v*sin(pitch_angle*pi/180)
bx0,by0,bz0=getb(R0,0,0)
B=sqrt(bx0**2+by0**2+bz0**2)

vx0=(vpara*bx0+vpen*bx0*bz0/sqrt(bx0**2+by0**2))/B
vy0=(vpara*by0+vpen*by0*bz0/sqrt(bx0**2+by0**2))/B
vz0=(vpara*bz0-vpen*sqrt(bx0**2+by0**2))/B
y0=[R0,0,0,vx0,vy0,vz0]
T=2*pi*m/q/B
tau=T/64
t_win=[i*tau for i in range(2000*16)]
y=rk4(Lonz,t_win,y0)

plt.figure(1)
plt.plot(np.sqrt(y[:,0]**2+y[:,1]**2),y[:,2])


fig = plt.figure(figsize=[10,10])
ax = plt.axes(projection='3d')

[u,v]=np.meshgrid(np.arange(0,2*pi,2*pi/50),np.arange(0,1.5*pi,1.5*pi/50));
X=(0.5+0.5*np.cos(u))*np.cos(v); Y=(0.5+0.5*np.cos(u))*np.sin(v);
Z=0.5*np.sin(u)

xx=y[:,0]
yy=y[:,1]
zz=y[:,2]

theta=np.arange(0,2*pi,pi/40)

x1=3.7*np.sin(theta)
y1=3.7*np.cos(theta)
x2=2.3*np.sin(theta)
y2=2.3*np.cos(theta)
z1=np.zeros(len(theta))

#ax.plot3D(xx,yy,zz,color='r')
#ax.plot3D(x1,y1,z1)
#ax.plot3D(x2,y2,z1)
ax.plot_wireframe(X,Y,Z,rstride=2,cstride=2,color=[0.4,0.7,0.7])
#ax.plot_surface(X,Y,Z,rstride=1, cstride=1,cmap='viridis', edgecolor='none')
#t=np.linspace(0,800*T,1000)
point,=ax.plot([xx[0]],[yy[0]],[zz[0]],'yo',markersize=10)
line, = ax.plot([xx[0]], [yy[0]], [zz[0]], lw=3)
ani = animation.FuncAnimation(fig=fig,
                              func=animate,
                              frames=len(xx),
                              interval=200,
                              repeat=True,
                              blit=False)
plt.show()
