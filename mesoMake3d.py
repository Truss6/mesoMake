import numpy as np
import matplotlib
matplotlib.use('agg')
from matplotlib import pyplot as plt
from math import *
import pandas as pd
import os,csv
# from numba import jit
from sphere import plot_S

##### DOMAIN #####
class domain:
    x=None
    y=None
    z=None
    area=None
    system=None
    def __init__(self,x,y,z):
        self.x=x
        self.y=y
        self.z=z
        self.area=(x[1]-x[0])*(y[1]-y[0])*(z[1]-z[0])
    def plot(obj,N):
        f = plt.figure(dpi=150)
        ax = f.add_subplot(projection='3d')
        ax.set_xlim(obj.x)
        ax.set_ylim(obj.y)
        ax.set_zlim(obj.z)
        obj.system.plot(ax)
        plt.savefig('meso'+str(N).zfill(3)+'.png')
        plt.close()
    def rand(obj,N):
        x=np.random.rand(N)*(obj.x[1]-obj.x[0])+obj.x[0]
        y=np.random.rand(N)*(obj.y[1]-obj.y[0])+obj.y[0]
        z=np.random.rand(N)*(obj.z[1]-obj.z[0])+obj.z[0]
        return x,y,z
    ##### UPDATE DOMAIN #####
    def update(obj,h=0.125):
        exit_flag=True
        err=0
        N=obj.system.N_points
        F=[[0,0,0] for _ in range(N)]
        for i in range(N): #point of interest
            x0=obj.system.points[i].x
            y0=obj.system.points[i].y
            z0=obj.system.points[i].z
            r0=obj.system.points[i].r
            for j in range(N): # point of comparison
                if i!=j:
                    dx=x0-obj.system.points[j].x # seperation in the x
                    dy=y0-obj.system.points[j].y
                    dz=z0-obj.system.points[j].z# seperation in the y
                    d=(dx**2+dy**2+dz**2)**(1/2)
                    r=(r0+obj.system.points[j].r)
                    if d<r:
                        if (r-d)/r > 1e-1:
                            exit_flag=False
                            err+=1
                        # print(d)
                        n=[dx,dy,dz]/d
                        F[i]+=0.5*h*(r-d)*n
            if (x0-r0)<obj.x[0]:
                exit_flag=False;err+=1
                F[i][0]+=(obj.x[0]-(x0-r0))
            elif (x0+r0)>obj.x[1]:
                exit_flag=False;err+=1
                F[i][0]+=-((x0+r0)-obj.x[1])
            if (y0-r0)<obj.y[0]:
                exit_flag=False;err+=1
                F[i][1]+=(obj.y[0]-(y0-r0))
            elif (y0+r0)>obj.y[1]:
                exit_flag=False;err+=1
                F[i][1]+=-((y0+r0)-obj.y[1])
            if (z0-r0)<obj.z[0]:
                exit_flag=False;err+=1
                F[i][2]+=(obj.z[0]-(z0-r0))
            elif (z0+r0)>obj.z[1]:
                exit_flag=False;err+=1
                F[i][2]+=-((z0+r0)-obj.z[1])                
        # for i in range(N):
            # print(F[i])
            obj.system.points[i].x+=F[i][0]
            obj.system.points[i].y+=F[i][1]
            obj.system.points[i].z+=F[i][2]
            # obj.plot(i+1)
        return exit_flag,err

##### SYSTEM #####
class system:
    N_points=None
    points=[]
    def __init__(self,points):
        self.points=points
        self.N_points=len(points)
    def plot(obj,ax):
        for point in obj.points:
            point.plot(ax)
    def write(obj,fname):
        with open(fname,'w') as out:
            out.write('x,y,z,r,id\n')
            for point in obj.points:
                line=str(point.x)+','+str(point.y)+','+str(point.z)+','+str(point.r)+','+str(point.i)+'\n'
                out.write(line)
##### POINTS #####
class point:
    x=None
    y=None
    z=None
    r=None
    i=None
    c=None
    def __init__(self,x,y,z,r,i,c=None):
        self.x=x
        self.y=y
        self.z=z
        self.r=r
        self.i=i
        if c==None:
            self.c=np.random.rand(3)
        else:
            self.c=c
    def plot(obj,ax):
        circle=plot_S(obj.x,obj.y,obj.z,obj.r,obj.c,ax)
        # ax.add_patch(circle)
        
def makePoints(x,y,z,r,i=None,c=None): # this initializes a set of points
    points=[]
    if len(x)!=len(y):
        print('ERROR1')
        return None
    else:
        N_points=len(x)
    if len(r)==1:
        r=[r for _ in range(N_points)]
    if c==None:
        c=[None for _ in range(N_points)]
    for X,Y,Z,R,I,C in zip(x,y,z,r,i,c):
        points.append(point(X,Y,Z,R,I,C)) 
    return points

Acirc=lambda r: (4/3)*pi*(r**3) # area of a circle

##### MAIN #####
if __name__=='__main__':
    
    ##### POINTS PARAMS #####
    ## Species 1 ##
    vfrac1=0.5
    r1=0.02
    # A1=Acirc(r1)
    
    ## Species 2 ##
    vfrac2=0.2
    r2=0.05
    # A2=Acirc(r2)

    D=domain([-0.5,0.5],[0,1],[-0.5,0.5]) ## init domain ##
    # N1=int(D.area*vfrac1/A1)
    # N2=int(D.area*vfrac2/A2)
    
    A1=D.area*vfrac1
    A2=D.area*vfrac2
    ##### INIT POINTS #####
    r=[];ID=[];c=[];
    while A1>0:
        # r0=np.random.chisquare(10)*r1/10
        r0=r1
        A1-=Acirc(r0)
        r.append(r0)
        ID.append(0)
        c.append([1,0,0])
    # while A2>0:
    #     # r0=np.random.chisquare(5)*r2/5
    #     r0=r2
    #     A2-=Acirc(r0)
    #     r.append(r0)
    #     ID.append(2)
    #     c.append([0,0,1])
        
    x,y,z=D.rand(len(ID)) ## randomly creats points inside the domain ##
    
    ##### INIT SYS #####
    points=makePoints(x,y,z,r,ID,c) ## init set of points ##
    D.system=system(points) ## assign points to domain ##
    I=0;h=1;A=False
    D.plot(I)
    exit_flag=False
    
    f=plt.figure(dpi=100);ax=f.gca() ## init plotting ##
    
    for I in range(10):
        
        ## update ##
        exit_flag,err=D.update(h)
        D.plot(I+1)
        
        ## plot convergence ##
        print(err)
        ax.semilogy(I,err,'ko')
        
        ## check convergence or divergence ##
        if I>1:
            if err/ERR<0.03:
                print('Converged')
                break
            elif err>err0:
                h=h*0.9
                # print('Unable to Converge')
                # brek
            elif err==err0:
                if h<1:
                    h=h/0.9
                else:
                    if A:
                        print('Unable to Converge. Remaining Error: '+str(int(100*err/ERR))+'%')
                        break
                    else:
                        A=True
                        h=h*0.9
        else:
            ERR=err
        err0=err
    ## end loop ##
    f.savefig('convergence.png') ## save convergence plot ##
    D.system.write('meso.csv')