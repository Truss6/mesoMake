'''
MesoMake5 is the fifth iteration of the mesoMake program
this version is 3D and includes zonging

Zoning subdivides the domain into NX by NY by NZ subdomains

points are assigned to subdomains

points are solved against other points in the same domain and neighboring domains

only points in boarder domains are checked for boundary conditions

points undergo an eularian remap after the lagrangian update

Zoning is effective for thousands of points
    Zone sizing is limited by the largest diameter

'''

import numpy as np
import matplotlib
matplotlib.use('agg')
from matplotlib import pyplot as plt
from math import *
import pandas as pd
import os,csv,time
import matplotlib.patches as Patches
from sphere import plot_S


##### DOMAIN #####
class domain:
    x=None
    y=None
    z=None
    # X=None
    # Y=None
    # Z=None
    NX=None
    NY=None
    NZ=None
    area=None
    system=None
    subd=None
    def __init__(self,x,y,z):
        self.x=x
        self.y=y
        self.z=z
        self.area=(x[1]-x[0])*(y[1]-y[0])*(z[1]-z[0])

    def subDiv(self,NX,NY,NZ):
        self.NX=NX
        self.NY=NY
        self.NZ=NZ
        X=np.linspace(self.x[0],self.x[1],NX+1)
        Y=np.linspace(self.y[0],self.y[1],NY+1)
        Z=np.linspace(self.z[0],self.z[1],NY+1)
        self.subd=[[[None for _ in range(NX)] for _ in range(NY)] for _ in range(NZ)]
        for i in range(NX):
            for j in range(NY):
                for k in range(NZ):
                    self.subd[i][j][k]=subdomain([X[i],X[i+1]],[Y[j],Y[j+1]],[Z[k],[Z[k+1]]])
        return self
    def plot(obj,N):
        f = plt.figure(dpi=150)
        ax = f.add_subplot(projection='3d')
        # for i in range(obj.NX):
        #     for j in range(obj.NY):
        #         obj.subd[i][j].plot(ax)
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
        
        for I in range(obj.NX):
            if I==0:
                RII=[I,I+1]
            elif I==obj.NX-1:
                RII=[I-1,I]
            else:
                RII=[I-1,I,I+1]
            for J in range(obj.NY):
                # print("I = "+str(I)+" J = "+str(J))
                if J==0:
                    RJJ=[J,J+1]
                elif J==obj.NY-1:
                    RJJ=[J-1,J]
                else:
                    RJJ=[J-1,J,J+1]
                for K in range(obj.NZ):
                    if K==0:
                        RKK=[K,K+1]
                    elif K==obj.NZ-1:
                        RKK=[K-1,K]
                    else:
                        RKK=[K-1,K,K+1]
                    for II in RII:
                        for JJ in RJJ:
                            for KK in RKK:
                            
                                for i in obj.subd[I][J][K].points:
                                    F=[0,0,0]
                                    x0=obj.system.points[i].x
                                    y0=obj.system.points[i].y
                                    z0=obj.system.points[i].z
                                    r0=obj.system.points[i].r
                                    for j in obj.subd[II][JJ][KK].points: # point of comparison
                                        if i!=j:
                                            dx=x0-obj.system.points[j].x # seperation in the x
                                            dy=y0-obj.system.points[j].y # seperation in the y
                                            dz=z0-obj.system.points[j].z
                                            d=(dx**2+dy**2+dz**2)**(1/2)
                                            r=(r0+obj.system.points[j].r)
                                            if d<r:
                                                if (r-d)/r > 1e-1:
                                                    exit_flag=False
                                                    err+=1
                                                  
                                                n=[dx,dy,dz]/d
                                                F+=0.5*h*(r-d)*n
                                    if I==0:
                                        if (x0-r0)<obj.x[0]:
                                            exit_flag=False;err+=1
                                            F[0]+=(obj.x[0]-(x0-r0))
                                    elif I==obj.NX-1:
                                        if (x0+r0)>obj.x[1]:
                                            exit_flag=False;err+=1
                                            F[0]+=-((x0+r0)-obj.x[1])
                                    if J==0:
                                        if (y0-r0)<obj.y[0]:
                                            exit_flag=False;err+=1
                                            F[1]+=(obj.y[0]-(y0-r0))
                                    elif J==obj.NY-1:
                                        if (y0+r0)>obj.y[1]:
                                            exit_flag=False;err+=1
                                            F[1]+=-((y0+r0)-obj.y[1])
                                    if K==0:
                                        if (z0-r0)<obj.y[0]:
                                            exit_flag=False;err+=1
                                            F[2]+=(obj.z[0]-(z0-r0))
                                    elif K==obj.NZ-1:
                                        if (z0+r0)>obj.z[1]:
                                            exit_flag=False;err+=1
                                            F[2]+=-((z0+r0)-obj.z[1])
                                    obj.system.points[i].x+=F[0]
                                    obj.system.points[i].y+=F[1]
                                    obj.system.points[i].z+=F[2]
            # obj.plot(i+1)
        return exit_flag,err
    def assign(obj):
        points=obj.system.points
        for i in range(obj.NX):
            for j in range(obj.NY):
                for k in range(obj.NZ):
                # print("i = "+str(i)+" j = "+str(j))
                    obj.subd[i][j][k].assign(points,i==0,i==obj.NX-1,j==0,j==obj.NY-1,k==0,k==obj.NZ-1)
##### SUBDOMAIN #####
class subdomain:
    x=None
    y=None
    z=None
    points=[]
    def __init__(self,x,y,z):
        self.x=x
        self.y=y
        self.z=z
        # return self
    def IN(obj,X,Y,Z,xb,xt,yb,yt,zb,zt):
        x=obj.x;y=obj.y;z=obj.z
        return ((X>x[0]) or xb)*((X<x[1]) or xt)*((Y>y[0]) or yb)*((Y<y[1]) or yt)*((Z>z[0]) or zb)*((Z<z[1]) or zt)
    def assign(obj,points,xb,xt,yb,yt,zb,zt):
        obj.points=[]
        for i in range(len(points)):
            if obj.IN(points[i].x,points[i].y,points[i].z,xb,xt,yb,yt,zb,zt):
                obj.points.append(i)
    # def plot(obj,ax):
    #     R=Patches.Rectangle((obj.x[0],obj.y[0]),obj.x[1]-obj.x[0],obj.y[1]-obj.y[0],facecolor='w',edgecolor='k')
    #     ax.add_patch(R)
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
        # circle=plt.Circle((obj.x,obj.y),obj.r,color=obj.c,alpha=0.5)
        # ax.add_patch(circle)
        circle=plot_S(obj.x,obj.y,obj.z,obj.r,obj.c,ax)

        
def makePoints(x,y,z,r,i=None,c=None): # this initializes a set of points
    points=[]
    # if len(x)!=len(y):
    #     print('ERROR1')
    #     return None
    # else:
    N_points=len(x)
    if len(r)==1:
        r=[r for _ in range(N_points)]
    if c==None:
        c=[None for _ in range(N_points)]
    for X,Y,Z,R,I,C in zip(x,y,z,r,i,c):
        points.append(point(X,Y,Z,R,I,C)) 
    return points

Acirc=lambda r: (4/3)*pi*r**3 # area of a circle

##### MAIN #####
if __name__=='__main__':
    tic=time.time()
    np.random.seed(0)

    ##### POINTS PARAMS #####
    ## Species 1 ##
    vfrac1=0.5
    r1=0.02
    # A1=Acirc(r1)   
    
    ## Specie 2 ##
    vfrac2=0.5
    r2=0.05
    # A2=Acirc(r2)

    D=domain([0,1],[0,1],[0,1]) ## init domain ##
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
        ID.append(3)
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
    
    
    dMax=max(r)*2
    NX=floor(1*(D.x[1]-D.x[0])/dMax)
    NY=floor(1*(D.y[1]-D.y[0])/dMax)
    NZ=floor(1*(D.z[1]-D.z[0])/dMax)
    
    D=D.subDiv(NX,NY,NZ)
    D.assign()
    
    I=0;h=1.5;A=False
    D.plot(I)
    exit_flag=False
    
    f=plt.figure(dpi=100);ax=f.gca() ## init plotting ##
    
    for I in range(10):
        
        ## update ##
        
        ## lagrangian update ##
        exit_flag,err=D.update(h)
        
        ## eulerian update ##
        D.assign()
        
        ## plot ##
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
    print(time.time()-tic)
    f.savefig('convergence.png') ## save convergence plot ##
    D.system.write('meso.csv')
    plt.close()
    