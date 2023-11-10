'''
    This is an example mesoMake input deck
    The required information is: 
        vfrac, R,id -- these must have the same lengts, one per spicies
        X,Y -- these must have length 2, and be in ascending order
    
    Units are arbitrary
'''
vfrac=[0.125,1,0.2] #   volume fraction of each species
R=[0.2,0.05,0.125]   #   radius of each species 
id=[2,-1,1]       #   identity of each species

# vfrac=[0.25]
# R=[0.01]
# id=[0]

X=[-2,2]        #   X limits
Y=[0,1]         #   Y limits

'''
    Optional information
'''
color=['r',[0,0,1],'#00ff00'] #   user defined colors, must be in proper python color formats
#   if you dont define your own colors, it will default to random ugly colors

n_iter_max=25 # maximum number of solver iterations, this is a recommended input, but default is 10
earlyExit=False
# seed=0 # you may define a random seed for reproducable results

sobol=True     # if you do not want sobol sampling, set to false. default is true
# sobol sampling prevents clumps at the initiliazation, It is most helpful when constituents are similarily sized
# resave=True   # if you want a meso.csv file saved on every loop, set to true. default is false

ratio=1

imgFMT='png'