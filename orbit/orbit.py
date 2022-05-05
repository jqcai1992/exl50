# -*- coding: utf-8 -*-
"""
Created on Wed Dec  9 10:25:11 2020

@author: 83022
"""

from math import *
import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation
import json
from scipy import interpolate
import tqdm
from tqdm.std import trange


class Orbit(object):
    
    def __init__(self,shot_no,time,R0,K0,pitch_angle,direction):
        self.file_name = 'EXL%06.d_%05.fms.json'%(shot_no,time*1000)    
        with open(self.file_name) as file:
            self.exl=json.load(file)
        self.rgrid=self.exl['r_grid']
        self.zgrid=self.exl['z_grid']
        self.Br=self.exl['Br']
        self.Bt=self.exl['Bt']
        self.Bz=self.exl['Bz']
        self.R0=R0  #入射径向位置
        self.K0= K0  #6keV的电子
        self.pitch_angle = pitch_angle #65
        self.direction= direction
    
    
        self.Br_func = interpolate.interp2d(self.rgrid, self.zgrid,self.Br, kind='linear')
        self.Bt_func = interpolate.interp2d(self.rgrid, self.zgrid,self.Bt, kind='linear')
        self.Bz_func = interpolate.interp2d(self.rgrid, self.zgrid,self.Bz, kind='linear')
        self.rmin=min(self.rgrid)
        self.rmax=max(self.rgrid)
        self.zmin=min(self.zgrid)
        self.zmax=max(self.zgrid)
        
    def getb(self,x,y,z):
        
        fR=np.sqrt(x**2+y**2)
        fz=z
        if self.rmin<=fR<=self.rmax and self.zmin<=fz<=self.zmax:
            bR=self.Br_func(fR,fz)
            bt=self.Bt_func(fR,fz)
            bz=self.Bz_func(fR,fz)
            
            by=bR*y/fR+bt*x/fR
            bx=x/fR*bR-bt*y/fR
        else:
            bx=[0]
            by=[0]
            bz=[0]
        
        return bx[0],by[0],bz[0]   


    def getb1(self,x,y,z):
        R0=2
        B0=1
        r=sqrt((sqrt(x**2+y**2)-R0)**2+z**2)+1e-10
        R=sqrt(x**2+y**2)
        q=1+r**2
        bt=B0*R0/R
        bp=bt*r/q/R0
        bx=-y*bt/R-bp*z*x/r/R
        by=bt*x/R-bp*z/r*y/R
        bz=bp*(R-R0)/r
        return bx,by,bz    

    def Lonz(self,y):
        q = 1.602176565e-19; # Elementary charge (Coulomb)
        m = 1.672621777e-27; # Proton mass (kg)
        m_el = 9.10938291e-31
        
        dy=np.zeros(len(y))
        dy[0]=y[3]
        dy[1]=y[4]
        dy[2]=y[5]
        bx,by,bz=self.getb(y[0],y[1],y[2])
        dy[3]=q/m*(y[4]*bz-y[5]*by)
        dy[4]=q/m*(y[5]*bx-y[3]*bz)
        dy[5]=q/m*(y[3]*by-y[4]*bx)
        return dy

    def rk4(self,func,tau,t_seq=[],y0=[]):
        y=np.zeros([len(t_seq)+1,len(y0)])
        y[0]=y0
        for i in trange(len(t_seq)):              
            k1 = func(y[i])
            k2 = func(y[i]+k1*tau/2)
            k3 = func(y[i]+k2*tau/2)
            k4 = func(y[i]+k3*tau) 
            y[i+1]=y[i]+1/6*tau*(k1+2*k2+2*k3+k4)
        return y


    def trace(self):

        q = 1.602176565e-19; # Elementary charge (Coulomb)
        m = 1.672621777e-27; # Proton mass (kg)
        m_el = 9.10938291e-31; # Electron mass (kg)
        c = 299792458; # speed of light (m/s)
            
        v=sqrt(self.K0*q/m)
        vpara=self.direction*v*cos(self.pitch_angle*pi/180)  #调节正负
        vpen=v*sin(self.pitch_angle*pi/180)
        bx0,by0,bz0=self.getb(self.R0,0,0)
        B0=sqrt(bx0**2+by0**2+bz0**2)
        
        vx0=(vpara*bx0+vpen*bx0*bz0/sqrt(bx0**2+by0**2))/B0
        vy0=(vpara*by0+vpen*by0*bz0/sqrt(bx0**2+by0**2))/B0
        vz0=(vpara*bz0-vpen*sqrt(bx0**2+by0**2))/B0
        y0=[self.R0,0,0,vx0,vy0,vz0]
        self.T=2*pi*m/q/B0
        self.tau=self.T/16
        t_win=[i*self.tau for i in range(6000*16)]
        self.y=self.rk4(self.Lonz,self.tau,t_win,y0)
        
               
def animate(i):
    line.set_xdata(xx[:1000*i])
    line.set_ydata(yy[:1000*i])
    line.set_3d_properties(zz[:1000*i])
    point.set_xdata(xx[1000*i])
    point.set_ydata(yy[1000*i])
    point.set_3d_properties(zz[1000*i])   

if __name__== "__main__": 
    R0=0.4
    K0=2e6
    pitch_angle=15
    direction=1
    a=Orbit(6935,4.45,R0,K0,pitch_angle,direction)
    a.trace()
    psi=np.array(a.exl['psi'])
    Psi_lcfs_Wb_rad=a.exl['Psi_lcfs_Wb_rad']
    r_grid=a.exl['r_grid']
    z_grid=a.exl['z_grid']
    psi_norm=(psi-np.min(psi))/(Psi_lcfs_Wb_rad-np.min(psi))
    
    plt.figure(1)
    plt.contour(r_grid,z_grid,psi_norm,[0.01,0.04,0.15,0.3,0.5,0.7,0.9,1.5,2.5,4,6,8],colors='b')
    plt.contour(r_grid,z_grid,psi_norm,[1],colors='b',linewidths=2.5)
    plt.xlabel('$R(m)$',fontfamily='Times New Roman',fontsize=14)
    plt.ylabel('$Z(m)$',fontfamily='Times New Roman',fontsize=14)
    plt.axis('scaled')
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'
    plt.tick_params(labelsize=14)
    plt.xlim((0.165,1.3))
    plt.ylim((-1.2,1.2))  
        
    xx=a.y[:,0]
    yy=a.y[:,1]
    zz=a.y[:,2]
    plt.plot(np.sqrt(xx**2+yy**2),zz,'r')
    fig = plt.figure(2)
    ax = plt.axes(projection='3d')
    point, = ax.plot([xx[0]],[yy[0]],[zz[0]],'yo')
    line, = ax.plot([xx[0]], [yy[0]], [zz[0]], lw=1)


    [u,v]=np.meshgrid(np.arange(0,2*pi,2*pi/100),np.arange(0,1.5*pi,1.5*pi/50));
    X=(1+0.6*np.cos(u))*np.cos(v); Y=(1+0.6*np.cos(u))*np.sin(v);

    Z=0.6*np.sin(u)
    ax.plot_wireframe(X,Y,Z,rstride=2,cstride=2,color=[0.4,0.7,0.7])

    ani = animation.FuncAnimation(fig=fig,
                              func=animate,
                              frames=len(yy),
                              interval=200,
                              repeat=False,
                              blit=False)


    plt.show()