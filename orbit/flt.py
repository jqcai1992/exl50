# -*- coding: utf-8 -*-
"""
Created on Tue May 25 14:57:59 2021

@author: 83022
"""
import numpy as np
import json
import matplotlib.pyplot as plt
from scipy import interpolate
from mpl_toolkits.mplot3d import Axes3D
import tqdm
from tqdm.std import trange
from orbit import Orbit


class Flt(Orbit):
    
    def __init__(self,shot_no,time):
        self.file_name = 'EXL%06.d_%05.fms.json'%(shot_no,time*1000)
        Orbit.__init__(self,shot_no,time,1,1,1,1)
       
        self.steps=100000   #追踪步数
        self.dphi=2*np.pi/360   #每步循环dphi
        self.r0=[self.exl['R_lcfs_out_m']]   #中平面上起始追踪位置的大半径

        
    def dr_dz(self,y):
        """
        计算dR/dphi以及dZ/dphi
        """
        dy=np.zeros(2)
        Br1=self.Br_func(y[0],y[1]) #fr=y[0] fz=y[1]
        Bt1=self.Bt_func(y[0],y[1])
        Bz1=self.Bz_func(y[0],y[1])
        
        dy[0]=Br1*y[0]/Bt1   #dr_dphi
        dy[1]=Bz1*y[0]/Bt1    #dz_dphi
        
        return dy
    
    
    def run_flt(self):
        
        phi_window=[i*self.dphi for i in range(self.steps)]
        self.x=np.zeros((len(self.r0),len(phi_window)))
        self.y=np.zeros((len(self.r0),len(phi_window)))
        self.z=np.zeros((len(self.r0),len(phi_window)))
        for i in range(len(a.r0)):
            y0=[self.r0[i],0]
            self.yy=self.rk4(self.dr_dz,self.dphi,phi_window,y0)
            self.x[i]=self.yy[1:,0]*np.cos(phi_window)
            self.y[i]=self.yy[1:,0]*np.sin(phi_window)
            self.z[i]=self.yy[1:,1]
        
        
    def plot_flt(self,number=None):
        
        ax1 = plt.axes(projection='3d')  
        if number is None:
            for i in range(len(self.r0)):
                ax1.plot3D(a.x[i],a.y[i],a.z[i])
        elif isinstance(number,int):
            ax1.plot(a.x[:,number],a.y[:,number],a.z[:,number])
        else:
            for i in number:
                ax1.plot(a.x[:,i],a.y[:,i],a.z[:,i])
        ax1.set_xlabel('X(m)')
        ax1.set_ylabel('Y(m)')
        ax1.set_zlabel('Z(m)') 

        
        
if __name__== "__main__":         
    a=Flt(6935,4.45)   
    a.steps=100000  # 追踪步数
    a.r0 = [0.8]  #追踪点数
    #a.rk44()   
    a.run_flt()
    a.plot_flt()


            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
