# ---------------------- importing modules

import warnings
warnings.filterwarnings('ignore')

import sys, os
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from matplotlib import cm
from itertools import permutations, combinations

from mpl_toolkits.axes_grid1.inset_locator import (inset_axes, mark_inset)
from shapely.validation import make_valid
from shapely.ops import unary_union
from shapely.geometry import Point, Polygon

from psc.g_space import g, hsurf_g, F, hsurf_F, hsurf_F2
from psc.x2Dlinearize import double_segment_EPA, double_segment_nEPA, single_segment_EPA, single_segment_nEPA 
from psc.x2Dpolygon import multistrip, getploygons_EPA_SS, getploygons_EPA_DS, getploygons_nEPA, polyintersect, get_error
from psc.x2Drepetition import linrep_DS, linrep_SS, writedata, getsigncom, mesh
from psc.x2Dplot import plot_segment, plotisosurf_EPA, plotisosurf_nEPA
from psc.x2Dwritesolution import pseudosolution, realsolution, isInside, analyzesolution




# ---------------------- nEPA MODEL



def solve2Dstrucutre_nEPA(totalRO:int, scatteringfactors:list, structure:list,
                                      SingleorDoublesegement:str='single',
                                      isosurfaceplot: bool=False,
                                      plotresult: bool=False) -> None:
    
    """
    This routine solve structure within nEPA framework.
    
    Args:
        totalRO (int): number of reflection orders to consider in calculation. should be >=2 
        scatteringfactors (list): list of atomic structure factors of each atom. defult 
        strucutre (list): structure to solve. it should of of the form [x1, x2].
        SingleorDoublesegement (str, optional): Single- or Double-segment approach. Defaults to 'single'. for Double-segment
                                                use 'double'
        isosurfaceplot (bool, optional): decides to plot generated isosurface of each reflection order. Defaults to False.
        plotresult (bool, optional): decides to plot final result. Defaults to False.
    """    
    
       
    f = scatteringfactors if scatteringfactors else [10, 4]
    xcoor = structure if structure else np.random.uniform(0.0, 0.5, size=(1, 2))[0]
    
    if len(xcoor) >2:
        print(f"---> The structure should be a list of length 2. Given has length {len(xcoor)}. I am stooping now")
        sys.exit()
    elif len(f) >2:
        print(f"---> The atomic scattering factor should be a list of length 2. Given has length {len(f)}. I am stooping now")
        sys.exit()
    else:
        print(f"\x1b[1;31m--> The given or assumed coordinate \x1b[1;32m{xcoor}\x1b[1;31m is treated with atomic scattering factors \x1b[1;32m{f}\x1b[1;31m with reflection orders upto \x1b[1;32m{totalRO}\x1b[0m")
    
    
    # ---> assume we solve for second atomic coordinate 
    j     = 1

    # ---> defining number of reflection orders to consider
    h   = totalRO
    hss = np.arange(1,totalRO+1)

    # ---> generating possible atomic corrdinates in 2D PS. hence grid has dimension of 2x500. 
    isos = np.linspace(0., 0.5, 500)
    grid = [isos]*len(f)

    if isosurfaceplot:
        # ---> define plot continer
        plt.rc('xtick', labelsize=16) ; plt.rc('ytick', labelsize=16)
        fig, ax = plt.subplots(nrows=1, ncols=1, figsize = (12,4), subplot_kw = dict(aspect = 1.0))
        ax.set_xlim(isos.min(), isos.max()) ; ax.set_ylim(isos.min(), isos.max())
        ax.set_xlabel(r'$z_1$' ,fontsize=21,labelpad=12) ; ax.set_ylabel(r'$z_2$' ,fontsize=21,labelpad=12)

    ### open a dat file to store data
    fn=open("temppnts_%g.dat"%(h), "wt+")  

    ### staring amplitude calculation for each RO h and plot them
    for ic, l in enumerate(hss):
        xcoor  = np.sort(xcoor)[::-1]
        
        # calculate amplitude for given strucutre and RO
        I = F(l, xcoor, f)**2
        
        if isosurfaceplot:
            # calculate the isosurface over entrie PS using above gi for s=+1 and s=-1
            y1 = hsurf_F2(I, l, grid, f, j, s =  1, s2 =  1)  # c= k-
            y2 = hsurf_F2(I, l, grid, f, j, s =  1, s2 = -1)  # c= r-
            y3 = hsurf_F2(I, l, grid, f, j, s = -1, s2 =  1)  # c= b--
            y4 = hsurf_F2(I, l, grid, f, j, s = -1, s2 = -1)  # c= g--

            
            # plot calculated isosurfcae. define cc='k' if same is wanted. or isosurface colour will change automatically
            r = np.random.uniform(0.0, 0.8, 3) ; cc = (r[0],r[1],r[2],1)
            plotisosurf_nEPA(l, h, I, ax, isos, y1, y2, y3, y4, cc, lw=2, imax=0.5)
            
        #---> inearization process with error of err 
        errr = 0
        meshlist = mesh(l, f, imax=0.5)
        
        if SingleorDoublesegement == 'single':
            #---> single segment method - EPA
            pnts = single_segment_nEPA(np.sqrt(I), l, f, j, error=0)
            plist=linrep_SS(l, f, pnts, meshlist, imin=0, imax=0.5)
        
        elif SingleorDoublesegement == 'double':
            #---> double segment method - EPA
            pnts  = double_segment_nEPA(np.sqrt(I), l, f, j, error=0)
            plist = linrep_SS(l, f, pnts, meshlist, imin=0, imax=0.5)
        
        #---> plot segments
        if isosurfaceplot:
            plot_segment(ax, plist, cc)
        
        #---> storing segment data
        writedata(fn, plist)
        
    fn.close()
    
    if isosurfaceplot:
        fig.tight_layout()
        plt.legend(loc=2, prop={'size': 9})#, bbox_to_anchor=(1.0, 1.020))
        plt.show()
    
    fpath = os.path.join(os.getcwd(),'temppnts_%g.dat'%(h))
    with open(fpath) as finput:
        pointlst = [[float(num) for num in line.split()] for line in finput]

    poly = getploygons_nEPA(h, pointlst, imax=0.5)

    with open ('info.dat', "at+") as finfo:
        solution = polyintersect(h, poly, finfo)
    
    analyzesolution(solution, xcoor, plotting=plotresult)

    return



# ---------------------- EPA MODEL

def solve2Dstrucutre_EPA(totalRO:int, structure:list, IorG='intensity',
                                      SingleorDoublesegement:str='single',
                                      isosurfaceplot: bool=False,
                                      plotresult: bool=False, imax=0.5) -> None:
    """
    This routine solve structure within EPA framework. This module can apply 'intensity' or 'amplitude' methods.
    
    Args:
        totalRO (int): number of reflection orders to consider in calculation. should be >=2
        strucutre (list): structure to solve. it should of of the form [x1, x2].
        IorG (str, optional): decide which method ('intensity' or 'amplitude') to use.Defaults to 'intensity'.
        SingleorDoublesegement (str, optional): Single- or Double-segment approach. Defaults to 'single'. for Double-segment use 'double'
        isosurfaceplot (bool, optional): decides to plot generated isosurface of each reflection order. Defaults to False.
        plotresult (bool, optional): decides to plot final result. Defaults to False.
    """
    
    # EPA model
    f = [1, 1]
    
    # ---> The structure:
    if len(structure) >2:
        print(f"---> The structure should be a list of length 2. Given has length {len(structure)}. I am stooping now")
        sys.exit()
    else:
        print(f"\x1b[1;31m--> The given or assumed coordinate \x1b[1;32m{structure}\x1b[1;31m is treated with atomic scattering factors \x1b[1;32m{f}\x1b[1;31m with reflection orders upto \x1b[1;32m{totalRO}\x1b[0m")
    
    # ---> implementing permutation combination and origin of choice
    
    Asym = Polygon([[0,0],[0.5,0 ],[0.25,0.25]])
    
    xcoor = structure if structure else np.random.uniform(0.0, 0.5, size=(1, 2))[0]
    xcoor  = np.sort(xcoor)[::-1] 
    xcoor = xcoor if Asym.contains(Point(xcoor)) else np.sort(0.5-xcoor)[::-1]
    
    # ---> assume we solve for second atomic corrdinate 
    j = 1

    # ---> defining number of reflection orders to consider
    h   = totalRO
    hss = np.arange(1,totalRO+1)
    
    if isosurfaceplot:
        ### generating possible atomic corrdinates in 2D PS. hence grid has dimension of 2x500. 
        isos = np.linspace(0., imax/1, 500)  # !!! to change the precision in isosurface change 500 to higher number
        grid = [isos]*len(f)

        ### define plot container
        fig, ax = plt.subplots(nrows=1, ncols=1, figsize = (12,4), subplot_kw = dict(aspect = 1.0))
        fig.tight_layout()
        ax.set_xlim(isos.min(), isos.max()) ; ax.set_ylim(isos.min(), isos.max())
        ax.set_xlabel(r'$z_1$' ,fontsize=21,labelpad=12) ; ax.set_ylabel(r'$z_2$' ,fontsize=21,labelpad=12)
        plt.rc('xtick', labelsize=16) ; plt.rc('ytick', labelsize=16)

    ### open a dat file to store data
    fn=open("xpnts_%g.dat"%(h), "wt+")  

    ### staring amplitude calculation for each RO h and plot them
    for ic, l in enumerate(hss):
        
        # calculate amplitude for given structure and RO
        gi    = np.abs(g(l, xcoor, f))
        
        
        # calculate the isosurface over entire PS and plot
        if isosurfaceplot:
            
            giso1 = hsurf_g(l, grid, f, gi, j, s=1)
            giso2 = hsurf_g(l, grid, f, gi, j, s=-1)
            
            r = np.random.uniform(0.0, 0.8, 3) ; cc = (r[0],r[1],r[2],1)
            plotisosurf_EPA(l, h, gi, ax, isos, giso1, giso2, cc, lw=2, imax=0.5)
        
        #---> inearization process with error of err 
        errr = 0
        meshlist = mesh(l, f, imax)

        #---> Select which approach you want
        if SingleorDoublesegement == 'double':
            #---> double-segment method - EPA
            pnts = double_segment_EPA(gi, l, f, error=0)
            plist=linrep_DS(l, f, pnts, meshlist, IorG, signofIorG=np.sign(g(l, xcoor, f)), imin=0, imax=0.5)
        
        elif SingleorDoublesegement == 'single':
            #---> single-segment method - EPA
            pnts = single_segment_EPA(gi, l, f, error=0)
            plist=linrep_SS(l, f, pnts, meshlist, IorG, signofIorG=np.sign(g(l, xcoor, f)),imin=0, imax=0.5)

        #---> plot segments
        if isosurfaceplot:
            plot_segment(ax, plist, cc)

        #---> storing segment data
        writedata(fn, plist)

    fn.close()
    if isosurfaceplot:
        plt.legend(loc=2, prop={'size': 9})#, bbox_to_anchor=(1.0, 1.020))
        plt.show()
    
    #fpath= Path("C:/Users/pearl/Desktop/2021_Freiberg/")
    pointlst = np.loadtxt(os.path.join(os.getcwd(),'xpnts_%g.dat'%(h)),delimiter='\t')
    
    if SingleorDoublesegement == 'double':
        poly = getploygons_EPA_DS(h, pointlst, IorG, imax=0.5)
    elif SingleorDoublesegement == 'single':
        poly = getploygons_EPA_SS(h, pointlst, IorG, imax=0.5)
    
    with open ('info.dat', "at+") as finfo:
        solution = polyintersect(h, poly, finfo)

    analyzesolution(solution, xcoor, plotting=plotresult)
    
    return
