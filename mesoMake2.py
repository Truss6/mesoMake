import numpy as np
import matplotlib
matplotlib.use('agg')
from matplotlib import pyplot as plt
from math import *
import pandas as pd
import os,csv,time

##### DOMAIN #####
class domain:
    x=None
    y=None
    area=None
    system=None
    def __init__(self,x,y):
        self.x=x
        self.y=y
        self.area=(x[1]-x[0])*(y[1]-y[0])
    def plot(obj,N):
        f = plt.figure(dpi=150)
        ax = plt.gca()
        plt.xlim(obj.x)
        plt.ylim(obj.y)
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
        F=[[0,0] for _ in range(N)]
        for i in range(N): #point of interest
            x0=obj.system.points[i].x
            y0=obj.system.points[i].y
            r0=obj.system.points[i].r
            for j in range(N): # point of comparison
                if i!=j:
                    dx=x0-obj.system.points[j].x # seperation in the x
                    dy=y0-obj.system.points[j].y # seperation in the y
                    d=(dx**2+dy**2)**(1/2)
                    r=(r0+obj.system.points[j].r)
                    if d<r:
                        if (r-d)/r > 1e-1:
                            exit_flag=False
                            err+=1
                            # print(d)
                        n=[dx,dy]/d
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
        # for i in range(N):
            # print(F[i])
            obj.system.points[i].x+=F[i][0]
            obj.system.points[i].y+=F[i][1]
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
    vfrac1=0.3
    r1=0.125
    # A1=Acirc(r1)   
    
    ## Species 2 ##
    vfrac2=0.7
    r2=0.1
    # A2=Acirc(r2)
    # A2=Acirc(r2)

    D=domain([0,1],[0,1]) ## init domain ##
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
    I=0;h=1;A=False
    D.plot(I)
    exit_flag=False
    
    f=plt.figure(dpi=100);ax=f.gca() ## init plotting ##
    
    for I in range(100):
        
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
    print(time.time()-tic)
    f.savefig('convergence.png') ## save convergence plot ##
    D.system.write('meso.csv')
    plt.close()