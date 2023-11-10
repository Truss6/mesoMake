'''
MesoMake4 is the fourth iteration of the mesoMake program
this version is 2D and includes zongin

Zoning subdivides the domain into NX by NY subdomains

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

##### DOMAIN #####
class domain:
    x=None
    y=None
    X=None
    Y=None
    NX=None
    NY=None
    area=None
    system=None
    subd=None
    def __init__(self,x,y):
        self.x=x
        self.y=y
        self.area=(x[1]-x[0])*(y[1]-y[0])

    def subDiv(self,NX,NY):
        self.NX=NX
        self.NY=NY
        X=np.linspace(self.x[0],self.x[1],NX+1)
        Y=np.linspace(self.y[0],self.y[1],NY+1)
        self.subd=[[None for _ in range(NY)] for _ in range(NX)]
        for i in range(NX):
            for j in range(NY):
                self.subd[i][j]=subdomain([X[i],X[i+1]],[Y[j],Y[j+1]])
        return self
    def plot(obj,N):
        f = plt.figure(dpi=150)
        ax = plt.gca()
        plt.xlim(obj.x)
        plt.ylim(obj.y)
        for i in range(obj.NX):
            for j in range(obj.NY):
                obj.subd[i][j].plot(ax)     
        obj.system.plot(ax)
        plt.savefig('meso'+str(N).zfill(3)+'.png')
        plt.close()
    def rand(obj,N):
        x=np.random.rand(N)*(obj.x[1]-obj.x[0])+obj.x[0]
        y=np.random.rand(N)*(obj.y[1]-obj.y[0])+obj.y[0]
        return x,y
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
                for II in RII:
                    for JJ in RJJ:
                        for i in obj.subd[I][J].points:
                            F=[0,0]
                            x0=obj.system.points[i].x
                            y0=obj.system.points[i].y
                            r0=obj.system.points[i].r
                            for j in obj.subd[II][JJ].points: # point of comparison
                                if i!=j:
                                    dx=x0-obj.system.points[j].x # seperation in the x
                                    dy=y0-obj.system.points[j].y # seperation in the y
                                    d=(dx**2+dy**2)**(1/2)
                                    r=(r0+obj.system.points[j].r)
                                    if d<r:
                                        if (r-d)/r > 1e-1:
                                            exit_flag=False
                                            err+=1
                                          
                                        n=[dx,dy]/d
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
                            obj.system.points[i].x+=F[0]
                            obj.system.points[i].y+=F[1]
            # obj.plot(i+1)
        return exit_flag,err
    def assign(obj):
        points=obj.system.points
        for i in range(obj.NX):
            for j in range(obj.NY):
                # print("i = "+str(i)+" j = "+str(j))
                obj.subd[i][j].assign(points,i==0,i==obj.NX-1,j==0,j==obj.NY-1)
##### SUBDOMAIN #####
class subdomain:
    x=None
    y=None
    points=[]
    def __init__(self,x,y):
        self.x=x
        self.y=y

        # return self
    def IN(obj,X,Y,xb,xt,yb,yt):
        x=obj.x;y=obj.y
        return ((X>x[0]) or xb)*((X<x[1]) or xt)*((Y>y[0]) or yb)*((Y<y[1]) or yt)
    def assign(obj,points,xb,xt,yb,yt):
        obj.points=[]
        for i in range(len(points)):
            if obj.IN(points[i].x,points[i].y,xb,xt,yb,yt):
                obj.points.append(i)
    def plot(obj,ax):
        R=Patches.Rectangle((obj.x[0],obj.y[0]),obj.x[1]-obj.x[0],obj.y[1]-obj.y[0],facecolor='w',edgecolor='k')
        ax.add_patch(R)
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
            out.write('x,y,r,id\n')
            for point in obj.points:
                line=str(point.x)+','+str(point.y)+','+str(point.r)+','+str(point.i)+'\n'
                out.write(line)
##### POINTS #####
class point:
    x=None
    y=None
    r=None
    i=None
    c=None
    def __init__(self,x,y,r,i,c=None):
        self.x=x
        self.y=y
        self.r=r
        self.i=i
        if c==None:
            self.c=np.random.rand(3)
        else:
            self.c=c
    def plot(obj,ax):
        circle=plt.Circle((obj.x,obj.y),obj.r,color=obj.c,alpha=0.5)
        ax.add_patch(circle)
        
def makePoints(x,y,r,i=None,c=None): # this initializes a set of points
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
    for X,Y,R,I,C in zip(x,y,r,i,c):
        points.append(point(X,Y,R,I,C)) 
    return points

Acirc=lambda r: pi*r**2 # area of a circle

##### MAIN #####
if __name__=='__main__':
    tic=time.time()
    np.random.seed(0)

    ##### POINTS PARAMS #####
    ## Species 1 ##
    vfrac1=0.5
    r1=0.02
    # A1=Acirc(r1)   
    
    ## Species 2 ##
    vfrac2=0.2
    r2=0.01
    # A2=Acirc(r2)

    D=domain([-2,2],[0,1]) ## init domain ##
    # N1=int(D.area*vfrac1/A1)
    # N2=int(D.area*vfrac2/A2)
    
    A1=D.area*vfrac1
    A2=D.area*vfrac2
    ##### INIT POINTS #####
    r=[];ID=[];c=[];id=2;
    while A1>0:
        # r0=np.random.chisquare(10)*r1/10
        r0=r1
        A1-=Acirc(r0)
        r.append(r0)
        ID.append(id)
        c.append([1,0,0])
        id+=1
    while A2>0:
        # r0=np.random.chisquare(5)*r2/5
        r0=r2
        A2-=Acirc(r0)
        r.append(r0)
        ID.append(2)
        c.append([0,0,1])
        
    x,y=D.rand(len(ID)) ## randomly creats points inside the domain ##
    
    ##### INIT SYS #####
    points=makePoints(x,y,r,ID,c) ## init set of points ##
    D.system=system(points) ## assign points to domain ##
    
    
    dMax=max(r)*4
    NX=floor(1*(D.x[1]-D.x[0])/dMax)
    NY=floor(1*(D.y[1]-D.y[0])/dMax)
        
    D=D.subDiv(NX,NY) # this is throwing an error
    D.assign()
    
    I=0;h=1;A=False
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
    