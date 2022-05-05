# -*- coding: utf-8 -*-
"""
Created on Thu Feb  4 09:59:22 2021
输入shot，time
exl_data.txt2json() txt数据转为js格式
exl_data.readjson() 读取json格式
exl_data.txt2gfile() 将txt数据转为gfile文本格式

@author: 蔡剑青
"""
import numpy as np
import json
import matplotlib.pyplot as plt
from scipy import interpolate

shot_no=6935
time_point=4.45

class FunctionJson:
    
    def __init__(self,shot,time):
        self.shot=shot
        self.time=time
        self.file_path='MultiFluid_output/shot%06d/%06d.%05.f/'%(shot,shot,1000*time)

        self.file_name = 'EXL%06.d_%05.fms.json'%(shot,time*1000)
        self.name=['Br','Bt','Bz','density_eh','density_el',
                   'density_p','Jeht','Jelt','Jpt','Jt_total',
                   'psi','Teh','Tel','Tp','Uehr','Ueht','Uehz',
                   'Uelr','Uelt','Uelz','Upr','Upt','Upz','VE_el']
    
    def txt2json(self):

        with open(self.file_path+self.file_name,'w') as fload:
            temp_dict={}
    
            with open (self.file_path+'equilibrium_data.txt') as f1:
                temp2=[]
                g1=f1.readlines()
                for line in g1:
                    a=[i for i in line.split('\t')]
                    try:
                        if (a[0]=='shot_No') or (a[0]=='nw') or (a[0]=='nh'):
                            temp_dict[a[0]]=int(a[1])
                        else:
                            temp_dict[a[0]]=float(a[1])
                            
                    except:
                        0
            for i in range(len(self.name)):
                with open (self.file_path+self.name[i]+'.txt') as f2:
                    temp=[]
                    g2=f2.readlines()
                    for line in g2[3:]:
                        line1=[float(i) for i in line.split()]
                        temp.append(line1)
                        temp_dict[self.name[i]]=temp
            
            with open (self.file_path+'psi.txt') as f3:
                temp2=[]
                g3=f3.readlines()
                z_grid=[float(i) for i in g3[1].split()]
                r_grid=[float(i) for i in g3[2].split()]
                temp_dict['z_grid']=z_grid
                temp_dict['r_grid']=r_grid
                
            json.dump(temp_dict,fload)

    def readjson(self):   
        with open(self.file_path+self.file_name) as file:
            exl=json.load(file)
            for k,v in exl.items():
                if isinstance(v,list):
                    globals()[k]=np.array(v)
                else: 
                    globals()[k]=v
            print(exl.keys())
            
    def txt2gfile(self):
        
        #self.txt2json()
        self.readjson()
        
        f=open(self.file_path+'EXL%06.d_%05.fms'%(self.shot,self.time*1000),'w')
        

        #写入第一行
        f.write("%7s%14s%5s"%("EFITD", "06/21/2015", '#')) 
        f.write('%6g%8s%12g%4g%4g\n'%(shot_No,'%dms'%(1000*self.time),3,nw,nh))
        
        #写入2-4行
        rzero=r_grid[0]/2+r_grid[-1]/2
        B=np.sqrt(Br**2+Bt**2+Bz**2)
        cpasma=abs(Iplasma)*1000
        b_interpolate=interpolate.interp2d(r_grid,z_grid,B.transpose())
        bcenter=b_interpolate(rzero,zmid)
        f.write('%16.9E%16.9E%16.9E%16.9E%16.9E\n'%(xdim,zdim,rzero,rmin,zmid))
        f.write('%16.9E%16.9E%16.9E%16.9E%16.9E\n'%(R_mag_m,0,Psi_mag_Wb_rad,Psi_lcfs_Wb_rad,bcenter))
        f.write('%16.9E%16.9E%16.9E%16.9E%16.9E\n'%(cpasma,Psi_mag_Wb_rad,0,R_mag_m,0))
        f.write('%16.9E%16.9E%16.9E%16.9E%16.9E\n'%(0,0,Psi_lcfs_Wb_rad,0,0))
        
        #写入f(psi)
        nprof=int(nw)
        dpsi=(Psi_lcfs_Wb_rad-Psi_mag_Wb_rad)/nprof
        #归一化psi
        psi_nprof=np.linspace(Psi_mag_Wb_rad,Psi_lcfs_Wb_rad,nw)
        f_psi=[]

        for i in range(int(nw)):
            psi_contour=plt.contour(r_grid,z_grid,psi,[psi_nprof[i]])
            data0= psi_contour.allsegs[0][0]
            r_prof=data0[:,0];z_prof=data0[:,1]
            newfunc = interpolate.interp2d(r_grid, z_grid, Bt, kind='cubic')
            temp = [float(newfunc(r_prof[j],z_prof)[j]) for j in range(len(r_prof))]*r_prof
            f_psi.append(temp.mean())
            
        ##求ffprime
        ffprime=list(np.diff(f_psi)/dpsi*f_psi[:-1])
        ffprime.append(0)
      
        for i in range(nprof):
            if (i+1)%5==0:
                f.write('%16.9E\n'%f_psi[i])
            else:
                f.write('%16.9E'%f_psi[i])
        f.write('\n')     
        
        #写入P(psi)
        ev=1.6e-19
        
        P=Tel*density_el*ev*1000+Tp*density_p*ev*1000
        P_psi=[]
        for i in range(50):
            psi_contour=plt.contour(r_grid,z_grid,psi,[psi_nprof[i]])
            data0= psi_contour.allsegs[0][0]
            r_prof=data0[:,0];z_prof=data0[:,1]
            newfunc = interpolate.interp2d(r_grid, z_grid, P, kind='cubic')
            temp = [float(newfunc(r_prof[j],z_prof)[j]) for j in range(len(r_prof))]*r_prof
            P_psi.append(temp.mean())
            
        
        for i in range(nprof):
            if (i+1)%5==0:
                f.write('%16.9E\n'%f_psi[i])
            else:
                f.write('%16.9E'%f_psi[i])      
        f.write('\n')
        
         #写入ffprim
        for i in range(nprof):
            if (i+1)%5==0:
                f.write('%16.9E\n'%ffprime[i])
            else:
                f.write('%16.9E'%ffprime[i])       
        f.write('\n')
        
        #写入pprime
        for i in range(nprof):
            if (i+1)%5==0:
                f.write('%16.9E\n'%f_psi[i])
            else:
                f.write('%16.9E'%f_psi[i])          
        f.write('\n')
        
        
        #写入磁通psi(R,Z)    
        psi2=[ i for item in psi for i in item]
        
        
        for i in range(int(nh*nw)):
            if (i+1)%5==0:
                f.write('%16.9E\n'%psi2[i])
            else:
                f.write('%16.9E'%psi2[i])          
        f.write('\n')
        
        #写入q(psi)
        for i in range(nprof):
            if (i+1)%5==0:
                f.write('%16.9E\n'%f_psi[i])
            else:
                f.write('%16.9E'%f_psi[i])          
        f.write('\n')        
                
        #写入最外闭合磁面位置和限制器位置
        #plt.contour(r_grid,z_grid,psi,20)
        lcfs_contour=plt.contour(r_grid,z_grid,psi,[Psi_lcfs_Wb_rad])
        data0= lcfs_contour.allsegs[0][0]
        rbbbs=data0[:,0]
        zbbbs=data0[:,1]
        plt.plot(rbbbs,zbbbs,'r')
        nbbbs=np.zeros(2*len(rbbbs))
        for i in range(len(rbbbs)):
            nbbbs[2*i]=rbbbs[i]
            nbbbs[2*i+1]=zbbbs[i]

        nlim=7
        f.write('%3i'%len(rbbbs))
        f.write('%3i\n'%nlim)
        for i in range(len(nbbbs)):
            if (i+1)%5==0:
                f.write('%16.9E\n'%nbbbs[i])
            else:
                f.write('%16.9E'%nbbbs[i])          
        f.write('\n')
        data_lim=[]
        data_lim.extend([0.17,0,0.17,1.245,1.525,1.245,1.525])
        data_lim.extend([-1.245,0.17,-1.245,0.17,0,0.17,0])
        
        for i in range(2*nlim):
            if (i+1)%5==0:
                f.write('%16.9E\n'%data_lim[i])
            else:
                f.write('%16.9E'%data_lim[i])          
        f.write('\n')        
        return temp

if __name__== "__main__":
    exl_data=FunctionJson(shot_no,time_point)
    #exl_data.txt2json()
    #exl_data.readjson()
    aa=exl_data.txt2gfile()
    plt.contour(r_grid,z_grid,psi,40)      
