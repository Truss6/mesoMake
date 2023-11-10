import numpy as np
'''
    These are the default values
'''
n_iter_max = 10 # max number of iterations
seed = 0    # initial random seed
color = [list(np.random.rand(3)) for _ in range(100)]   # random colors
sobol=True      # sobol sample?
plotting=True   # plotting? Adds a significant ammount of time
resave=False    # save every iteration?
earlyExit=True  # exit on converg?
fileName='meso.csv'     # save to file name
ratio=1     # zoning ratio, there must be greater than 2 zones per direction
imgFMT='png'

dpi=150     # plotting resolution

idCount=False # increment id's

def reassign(i):
    if i<5: return True  # remap the first 5 loops
    if (i%10==0): return True   # and every 10 thereafter
    print("No Remap")
    return False
    # define which cycles to execute eularian remap
    # after first couple of cycles, there is little need to frequently remap
# reassign=lambda i: True # always reassign