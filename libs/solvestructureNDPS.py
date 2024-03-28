'''
Module created for structure determination in parameter space of dimension 3 or more.
'''

# ============== Load Modules for EPA ==============

import os, sys
import numpy as np
import polytope as pc

from datetime import datetime

from psc.g_space import g, F, hsurf_F2, hsurf_g
from psc.x3Dlinearizazion import linearizenD_EPA, linearizenD_nEPA
from psc.x3Drepetition import getpolytope_EPA , getpolytope_nEPA
from psc.x3Dintersection import find_intersection
from psc.x3Dreadwrite import wrtcoor, wrtdata, wrtallsolution, wrttime_mc



# ============== Modules for EPA ==============

def solvenDstrucutre_EPA(dimension: int, structure: list=[], totalRO: int=9, 
                         IorG: str='amplitude', imax=0.5) -> None:
    
    TS0=datetime.now()
    sn = []
    polyFirstRO = []
    
    
    #---> define asymmetric part of PS
    temp = np.tril(np.ones(shape=(dimension, dimension)) , 0 )
    temp = imax*np.vstack([[0]*dimension, temp])
    asym = pc.qhull(np.array(temp))
    
    dimension = dimension if dimension >=3 else 3
    
    # ---> Generating required no of random positions in list
    #rp = structure if structure else np.random.uniform(0.0, imax, size=(1, dimension))[0]
    if (isinstance(structure, list)) and structure:
        if len(structure)<3:
            print(f"This module is for structure of dimension >= 3. Your structure has dimension {len(structure)}. I am exciting now ")
            sys.exit()
        else:
            rp = structure
    elif not structure:
        rp = np.random.uniform(0.0, imax, size=(1, dimension))[0]
    else:
        print("===> I want structure as a list of dimension >= 3. I am exciting now")
        sys.exit()
    
    
    # ---> define result folder to save results
    fpath  = os.path.join(os.path.join(os.environ['USERPROFILE']), 'Desktop', datetime.now().strftime('NDstructure-'+'%Y-%m-%d-%H%M%S')) 
    
    if os.path.isdir(fpath):
        print("\x1b[0;34m===> Output files will be saved to \x1b[1;33m", fpath, "\x1b[0;34m location\n" )
    else:
        os.mkdir(fpath)
        print(f"\x1b[0;34m===> Dir \x1b[1;33mresults \x1b[0;34mis created. Output files will be saved to {fpath} location" )
    
    
    
    if (isinstance(totalRO, int) and totalRO <= 2):
        hss = np.arange(2, 3)
        fname = os.path.join(fpath,'pnew_%g.h5'%(totalRO))
    elif (isinstance(totalRO, int) and totalRO >2):
        hss = np.arange(2, totalRO+1)
        fname = os.path.join(fpath,'pnew_%g.h5'%(totalRO))
    elif (isinstance(totalRO, list) and len(totalRO) <= 2):
        totalRO.append(totalRO[-1]+1)
        hss = totalRO
        fname = os.path.join(fpath,'pnew_%g.h5'%(len(totalRO)))
    elif (isinstance(totalRO, list) and len(totalRO) > 2):
        hss = totalRO
        fname = os.path.join(fpath,'pnew_%g.h5'%(len(totalRO)))
    else:
        print(f"totalRO should be either int of value >=2 or list of length >=2")
    
    wrtcoor(fname, [rp])
    
    for rc, h in enumerate(hss): #range(2, totalRO+1):
        
        # ---> Step 0: Initialization
        
        TS=datetime.now()
               
        print(f"\n\x1b[1;31m===> for h = {h}\x1b[0m")    
        
        tinfo = []
        
        f     = [1.0]*dimension
        j     = len(f)-1
        
        xexp  = np.sort(rp)[::-1]
        ll    = 1
        if np.sign(g(ll, xexp, f))>0:
            xexp = xexp
        else:
            xexp = np.sort(0.5-xexp)[::-1]
        
        print(f"\x1b[1;34m===> The given/generated structure {rp} is taken as {xexp} \x1b[0m")
        
        
        
        # ---> Step 1: Get isomatten 
        
        tinfo.append(h) #; tinfo.append(rc+1)
        t1=datetime.now()
        tr=datetime.now().strftime("%H:%M:%S")
        print(f"\x1b[1;32mstarting: {tr}\x1b[0m",end=" ")
        
        #nar = isosurfs_EPA(h, xexp, f)
        
        nar  = [] 
        if h <= 2:
            for l in range(1,h+1):
                gi = g(l, xexp, f)
                normal, dist = linearizenD_EPA(l, f, np.abs(gi))
                nar.append([l, normal, dist, np.sign(gi)])
        else:
            gi = g(h, xexp, f)
            normal, dist = linearizenD_EPA(h, f, np.abs(gi))
            nar.append([h, normal, dist, np.sign(gi)])
        
        tlinearize=datetime.now()
        tinfo.append(tlinearize.timestamp() - t1.timestamp())
        
        tr=datetime.now().strftime("%H:%M:%S")
        print(f"\x1b[1;32miso:\x1b[0m{tr}",end=" ")
        
        
        
        # ---> Step 2: Get repetition
        polylist = []
        
        if h <= 2:
            for i in nar: # here i[0]->l,i[1]->normal,i[2]->distance,i[3]->amplitudesign
                o = getpolytope_EPA(i[0], i[1], i[2], i[3], IorG, imax=0.5)
                
                if i[0]==1:
                    polyFirstRO.append(pc.Region([o.intersect(asym)]))
                    polylist.append(o)
                else:
                    pplist=[]
                    for ipoly in o:
                        if (asym.intersect(ipoly)):
                            if (polyFirstRO[0].intersect(ipoly)):
                                pplist.append(ipoly)
                    polylist.append(pc.Region(pplist))
        else:
            pplist=[]
            i = nar[0]
            o = getpolytope_EPA(i[0], i[1], i[2], i[3], IorG, imax=0.5)
            for ipoly in o:
                if (asym.intersect(ipoly)):
                    if polyFirstRO[0].intersect(ipoly):
                        pplist.append(ipoly)
            polylist=pc.Region(pplist)
        
        tgetpoly=datetime.now()
        tinfo.append(tgetpoly.timestamp() - tlinearize.timestamp())
        
        tr=datetime.now().strftime("%H:%M:%S")
        print(f"\x1b[1;32mpoly: \x1b[0m{tr}",end=" ")
        
        
        
        # ---> Step 3: Get intersection , write all solution
        if h == 2:
            solun = find_intersection(polyFirstRO[0], polylist[1])
            wrtallsolution(fname, rc, solun)
            sn.append(solun)
        else:
            solun   = find_intersection(sn[-1], polylist)
            
            wrtallsolution(fname, rc, solun)
            sn.append(solun)
        
        tintersect=datetime.now()
        tinfo.append(tintersect.timestamp() - tgetpoly.timestamp())
        
        tr=datetime.now().strftime("%H:%M:%S")
        print(f"\x1b[1;32mintersection: \x1b[0m{tr}",end=" ")
        
        # 3.2 Sorting solution
        tr=datetime.now().strftime("%H:%M:%S")
        print(f"\x1b[1;32msort sol:\x1b[0m{tr}",end=" ")
        
        # 3.3 Final step
        sol_all, vol_all  = [], []
        volAsym , sel1    = 0, -1
        
        for jc, poly in enumerate(solun):
            if xexp in poly:
                
                sol_all.append(poly)
                if not vol_all:
                    sel1   = jc
                    vol_all.append(poly.volume)
                    volume = poly.volume
                    
                else:
                    sel1   = sel1   if vol_all[-1] <= poly.volume else jc
                    volume = volume if vol_all[-1] <= poly.volume else poly.volume 
                    vol_all.append(poly.volume)
            
            volAsym += pc.volume(poly)
        
        print("L(sol)= ", len(solun)," s =",sel1) 
        
        rr=((3/4)*(1/np.pi)*(volAsym))**(1/3)
        
        sortorder = np.argsort(rp)[::-1]
        localmat  = pc.extreme(solun[sel1])
        
        for zi, zj in enumerate(localmat):
            localmat[zi]=zj[sortorder.argsort()]
                
        dmax = np.max(localmat, axis=0)
        dmin = np.min(localmat, axis=0)
        err  = np.abs(dmax-dmin)/2
        wrtdata(fname, rc, pc.volume(solun[sel1]), err, solun[sel1], localmat, [volAsym, rr], len(solun))
        
        twrite=datetime.now()
        tinfo.append(twrite.timestamp()-tintersect.timestamp())
        
        t2  = datetime.now()
        tr=datetime.now().strftime("%H:%M:%S")
        print(f"\x1b[1;32mwrt:\x1b[0m{tr}\x1b[1;32m Total: \x1b[0m{t2-t1}")
        
        ttotal=datetime.now()
        tinfo.append(ttotal.timestamp()-t1.timestamp())
        
        wrttime_mc(rc, fname, tinfo)
        
        m = np.mean(localmat,0)
        print(f"---> The structure is predicted to be : {m} or {np.sort(0.5-m)[::-1]}")    
        
    TE=datetime.now()    
    print("Total: \x1b[0m{}".format(TE-TS0))
    
    return


# ============== Modules for non EPA ==============

def solvenDstrucutre_nEPA(dimension: int, structure: list=[], totalRO: int=9,
                          structurefactor: list=[], IorG: str='intensity',
                          imax=0.5) -> None:
    
    TS0=datetime.now()
    sn = []
    
    #---> define asymmetric part of PS
    temp = np.tril(np.ones(shape=(dimension, dimension)) , 0 )
    temp = imax*np.vstack([[0]*dimension, temp])
    asym = pc.qhull(np.array(temp))
    
    dimension = dimension if dimension >=3 else 3
    
    # ---> Generating required no of random positions in list
    if (isinstance(structure, list)) and structure:
        if len(structure)<3:
            print(f"This module is for structure of dimension >= 3. Your structure has dimension {len(structure)}. I am exciting now ")
            sys.exit()
        else:
            rp = structure
    elif not structure:
        rp = np.random.uniform(0.0, imax, size=(1, dimension))[0]
    else:
        print("===> I want structure as a list of dimension >= 3. I am exciting now")
        sys.exit()
    
    
    # ---> define result folder to save results
    fpath  = os.path.join(os.path.join(os.environ['USERPROFILE']), 'Desktop', datetime.now().strftime('NDstructure-'+'%Y%m%d-%H%M%S')) 
    
    if os.path.isdir(fpath):
        print(f"\x1b[0;34m===> Output files will be saved to \x1b[1;33m{fpath} \x1b[0;34m location\n" )
    else:
        os.mkdir(fpath)
        print(f"\x1b[0;34m===> Output files will be saved in the directory {fpath}" )
    
    
    
    if (isinstance(totalRO, int) and totalRO <= 2):
        hss = np.arange(2, 3)
        fname = os.path.join(fpath,'pnew_%g.h5'%(totalRO))
    elif (isinstance(totalRO, int) and totalRO >2):
        hss = np.arange(2, totalRO+1)
        fname = os.path.join(fpath,'pnew_%g.h5'%(totalRO))
    elif (isinstance(totalRO, list) and len(totalRO) <= 2):
        totalRO.append(totalRO[-1]+1)
        hss = totalRO
        fname = os.path.join(fpath,'pnew_%g.h5'%(len(totalRO)))
    elif (isinstance(totalRO, list) and len(totalRO) > 2):
        hss = totalRO
        fname = os.path.join(fpath,'pnew_%g.h5'%(len(totalRO)))
    else:
        print(f"totalRO should be either int of value >=2 or list of length >=2")
    
    wrtcoor(fname, [rp])
    
    for rc, h in enumerate(hss): #range(2, totalRO+1):
        
        # ---> Step 0: Initialization
        
        TS=datetime.now()
               
        print(f"\n\x1b[1;31m===> Processing Reflection order: {h}\x1b[0m") 
        
        tinfo = []
        
        if structurefactor:
            if len(structurefactor) == dimension:
                f = structurefactor
            else:
                print(f"\x1b[1;31m===> ERROR FOUND - \x1b[0;32mThe length of structurefactor {len(structurefactor)} is different from dimension({dimension}). Check your data. I am exciting now ")
                sys.exit()            
        else:
            print(f"\x1b[1;31m===> ERROR FOUND - \x1b[0;32mThe structurefactor is not acceptable. The structurefactor should be a list of dimension {dimension}. I am exciting now")
            sys.exit()
        
        j     = len(f)-1
        
        xexp  = np.sort(rp)[::-1]
        ll    = 1
        if np.sign(g(ll, xexp, f))>0:
            xexp = xexp
        else:
            xexp = 0.5-xexp #np.sort(0.5-xexp)[::-1]
        
        print(f"\x1b[1;34m===> The given/generated structure {rp} is taken as {xexp} \x1b[0m")
        
        
        
        # ---> Step 1: Get isomatten 
        
        tinfo.append(h) #; tinfo.append(rc+1)
        t1=datetime.now()
        tr=datetime.now().strftime("%H:%M:%S")
        print(f"\x1b[1;32mstarting: {tr}\x1b[0m",end=" ")
        
        #nar = isosurfs_EPA(h, xexp, f)
        
        nar  = [] 
        if h <= 2:
            for l in range(1,h+1):
                gi = g(l, xexp, f)
                I  = F(l, xexp, f)**2
                normal, dist = linearizenD_nEPA(l, f, np.sqrt(I))
                nar.append([l, normal, dist, np.sign(gi)])
        else:
            gi = g(h, xexp, f)
            I  = F(h, xexp, f)**2
            normal, dist = linearizenD_nEPA(h, f, np.sqrt(I))
            nar.append([h, normal, dist, np.sign(gi)])
                    
        tlinearize=datetime.now()
        tinfo.append(tlinearize.timestamp() - t1.timestamp())
        
        tr=datetime.now().strftime("%H:%M:%S")
        print(f"\x1b[1;32miso:\x1b[0m{tr}",end=" ")
        
        
        
        # ---> Step 2: Get repetition
        polylist = []
        
        if h <= 2:
            for i in nar: # here i[0]->l,i[1]->normal,i[2]->distance,i[3]->amplitudesign
                o = getpolytope_nEPA(i[0], i[1], i[2], i[3], IorG, imax=0.5)
                
                if i[0]==1:
                    polyFirstRO = o
                    polylist.append(o)
                    #print(f"\n---> len(o) : {len(o)} {type(polyFirstRO)} {pc.extreme(o[0])}")
                else:
                    pplist=[]
                    
                    for ipoly in o:
                        if (polyFirstRO.intersect(ipoly)):
                            pplist.append(ipoly)
                    #print(f"\n---> len(o) : {len(o)} len(pplist) : {len(pplist)} ")
                    polylist.append(pc.Region(pplist))
        else:
            pplist=[]
            o = getpolytope_nEPA(nar[0][0], nar[0][1], nar[0][2], nar[0][3], IorG, imax=0.5)
            
            for ipoly in o:
                if polyFirstRO.intersect(ipoly):
                    pplist.append(ipoly)
            #print(f"\n---> len(o) : {len(o)} len(pplist) : {len(pplist)}")
            polylist=pc.Region(pplist)
        
        tgetpoly=datetime.now()
        tinfo.append(tgetpoly.timestamp() - tlinearize.timestamp())
        
        tr=datetime.now().strftime("%H:%M:%S")
        print(f"\x1b[1;32mpoly: \x1b[0m{tr}",end=" ")
        
        
        
        # ---> Step 3: Get intersection , write all solution
        if h == 2:
            solun = find_intersection(polyFirstRO, polylist[1])
            wrtallsolution(fname, rc, solun)
            sn.append(solun)
        else:
            solun = find_intersection(sn[-1], polylist)
            
            wrtallsolution(fname, rc, solun)
            sn.append(solun)
        
        tintersect=datetime.now()
        tinfo.append(tintersect.timestamp() - tgetpoly.timestamp())
        
        tr=datetime.now().strftime("%H:%M:%S")
        print(f"\x1b[1;32mintersection: \x1b[0m{tr}",end=" ")
        
        # 3.2 Sorting solution
        tr=datetime.now().strftime("%H:%M:%S")
        print(f"\x1b[1;32msort sol:\x1b[0m{tr}",end=" ")
        
        # 3.3 Final step
        sol_all, vol_all  = [], []
        volAsym , sel1    = 0, -1
        
        for jc, poly in enumerate(solun):
            if xexp in poly:
                
                sol_all.append(poly)
                if not vol_all:
                    sel1   = jc
                    vol_all.append(poly.volume)
                    volume = poly.volume
                    
                else:
                    sel1   = sel1   if vol_all[-1] <= poly.volume else jc
                    volume = volume if vol_all[-1] <= poly.volume else poly.volume 
                    vol_all.append(poly.volume)
            
            volAsym += pc.volume(poly)
        
        print(f"L(sol)= {len(solun)} s = {sel1}") 
        
        rr=((3/4)*(1/np.pi)*(volAsym))**(1/3)
        
        sortorder = np.argsort(rp)[::-1]
        localmat  = pc.extreme(solun[sel1])
        
        for zi, zj in enumerate(localmat):
            localmat[zi]=zj[sortorder.argsort()]
                
        dmax = np.max(localmat, axis=0)
        dmin = np.min(localmat, axis=0)
        err  = np.abs(dmax-dmin)/2
        wrtdata(fname, rc, pc.volume(solun[sel1]), err, solun[sel1], localmat, [volAsym, rr], len(solun))
        
        twrite=datetime.now()
        tinfo.append(twrite.timestamp()-tintersect.timestamp())
        
        t2  = datetime.now()
        tr=datetime.now().strftime("%H:%M:%S")
        print(f"\x1b[1;32mwrt:\x1b[0m{tr}\x1b[1;32m Total: \x1b[0m{t2-t1}")
        
        ttotal=datetime.now()
        tinfo.append(ttotal.timestamp()-t1.timestamp())
        
        wrttime_mc(rc, fname, tinfo)
        
        m = np.mean(localmat,0)
        print(f"---> The structure is predicted to be : {m} or {np.sort(0.5-m)[::-1]}")    
        
    TE=datetime.now()    
    print("Total: \x1b[0m{}".format(TE-TS0))
    
    return





# def solvenDstrucutre_nEPA(dimension: int, structure: list=[], totalRO: int=9,
#                           structurefactor: list=[], IorG: str='intensity',
#                           imax=0.5) -> None:
    
#     TS0=datetime.now()
#     sn = []
#     polyFirstRO = []
    
    
#     #---> define asymmetric part of PS
#     temp = np.tril(np.ones(shape=(dimension, dimension)) , 0 )
#     temp = imax*np.vstack([[0]*dimension, temp])
#     asym = pc.qhull(np.array(temp))
    
#     dimension = dimension if dimension >=3 else 3
    
#     # ---> Generating required no of random positions in list
#     #rp = structure if structure else np.random.uniform(0.0, imax, size=(1, dimension))[0]
#     if (isinstance(structure, list)) and structure:
#         if len(structure)<3:
#             print(f"This module is for structure of dimension >= 3. Your structure has dimension {len(structure)}. I am exciting now ")
#             sys.exit()
#         else:
#             rp = structure
#     elif not structure:
#         rp = np.random.uniform(0.0, imax, size=(1, dimension))[0]
#     else:
#         print("===> I want structure as a list of dimension >= 3. I am exciting now")
#         sys.exit()
    
    
#     # ---> define result folder to save results
#     fpath  = os.path.join(os.path.join(os.environ['USERPROFILE']), 'Desktop', datetime.now().strftime('NDstructure-'+'%Y-%m-%d-%H%M%S')) 
    
#     if os.path.isdir(fpath):
#         print("\x1b[0;34m===> Output files will be saved to \x1b[1;33m", fpath, "\x1b[0;34m location\n" )
#     else:
#         os.mkdir(fpath)
#         print(f"\x1b[0;34m===> Dir \x1b[1;33mresults \x1b[0;34mis created. Output files will be saved to {fpath} location" )
    
    
    
#     if (isinstance(totalRO, int) and totalRO <= 2):
#         hss = np.arange(2, 3)
#         fname = os.path.join(fpath,'pnew_%g.h5'%(totalRO))
#     elif (isinstance(totalRO, int) and totalRO >2):
#         hss = np.arange(2, totalRO+1)
#         fname = os.path.join(fpath,'pnew_%g.h5'%(totalRO))
#     elif (isinstance(totalRO, list) and len(totalRO) <= 2):
#         totalRO.append(totalRO[-1]+1)
#         hss = totalRO
#         fname = os.path.join(fpath,'pnew_%g.h5'%(len(totalRO)))
#     elif (isinstance(totalRO, list) and len(totalRO) > 2):
#         hss = totalRO
#         fname = os.path.join(fpath,'pnew_%g.h5'%(len(totalRO)))
#     else:
#         print(f"totalRO should be either int of value >=2 or list of length >=2")
    
#     wrtcoor(fname, [rp])
    
#     for rc, h in enumerate(hss): #range(2, totalRO+1):
        
#         # ---> Step 0: Initialization
        
#         TS=datetime.now()
               
#         print(f"\n\x1b[1;31m===> for h = {h}\x1b[0m") 
        
#         tinfo = []
        
#         if structurefactor:
#             f = structurefactor
#         else:
#             print(f"\x1b[1;31m===> ERROR FOUND - The structurefactor is not acceptable. The structurefactor should be a list of dimension {dimension}. I am exciting now")
#             sys.exit()
        
#         j     = len(f)-1
        
#         xexp  = np.sort(rp)[::-1]
#         ll    = 1
#         if np.sign(g(ll, xexp, f))>0:
#             xexp = xexp
#         else:
#             xexp = np.sort(0.5-xexp)[::-1]
        
#         print(f"\x1b[1;34m===> The given/generated structure {rp} is taken as {xexp} \x1b[0m")
        
        
        
#         # ---> Step 1: Get isomatten 
        
#         tinfo.append(h) #; tinfo.append(rc+1)
#         t1=datetime.now()
#         tr=datetime.now().strftime("%H:%M:%S")
#         print(f"\x1b[1;32mstarting: {tr}\x1b[0m",end=" ")
        
#         #nar = isosurfs_EPA(h, xexp, f)
        
#         nar  = [] 
#         if h <= 2:
#             for l in range(1,h+1):
#                 gi = g(l, xexp, f)
#                 normal, dist = linearizenD_EPA(l, f, np.abs(gi))
#                 nar.append([l, normal, dist, np.sign(gi)])
#         else:
#             gi = g(h, xexp, f)
#             normal, dist = linearizenD_EPA(h, f, np.abs(gi))
#             nar.append([h, normal, dist, np.sign(gi)])
        
#         tlinearize=datetime.now()
#         tinfo.append(tlinearize.timestamp() - t1.timestamp())
        
#         tr=datetime.now().strftime("%H:%M:%S")
#         print(f"\x1b[1;32miso:\x1b[0m{tr}",end=" ")
        
        
        
#         # ---> Step 2: Get repetition
#         polylist = []
        
#         if h <= 2:
#             for i in nar: # here i[0]->l,i[1]->normal,i[2]->distance,i[3]->amplitudesign
#                 o = getpolytope_EPA(i[0], i[1], i[2], i[3], IorG, imax=0.5)
                
#                 if i[0]==1:
#                     polyFirstRO.append(pc.Region([o.intersect(asym)]))
#                     polylist.append(o)
#                 else:
#                     pplist=[]
#                     for ipoly in o:
#                         if (asym.intersect(ipoly)):
#                             if (polyFirstRO[0].intersect(ipoly)):
#                                 pplist.append(ipoly)
#                     polylist.append(pc.Region(pplist))
#         else:
#             pplist=[]
#             i = nar[0]
#             o = getpolytope_EPA(i[0], i[1], i[2], i[3], IorG, imax=0.5)
#             for ipoly in o:
#                 if (asym.intersect(ipoly)):
#                     if polyFirstRO[0].intersect(ipoly):
#                         pplist.append(ipoly)
#             polylist=pc.Region(pplist)
        
#         tgetpoly=datetime.now()
#         tinfo.append(tgetpoly.timestamp() - tlinearize.timestamp())
        
#         tr=datetime.now().strftime("%H:%M:%S")
#         print(f"\x1b[1;32mpoly: \x1b[0m{tr}",end=" ")
        
        
        
#         # ---> Step 3: Get intersection , write all solution
#         if h == 2:
#             solun = find_intersection(polyFirstRO[0], polylist[1])
#             wrtallsolution(fname, rc, solun)
#             sn.append(solun)
#         else:
#             solun   = find_intersection(sn[-1], polylist)
            
#             wrtallsolution(fname, rc, solun)
#             sn.append(solun)
        
#         tintersect=datetime.now()
#         tinfo.append(tintersect.timestamp() - tgetpoly.timestamp())
        
#         tr=datetime.now().strftime("%H:%M:%S")
#         print(f"\x1b[1;32mintersection: \x1b[0m{tr}",end=" ")
        
#         # 3.2 Sorting solution
#         tr=datetime.now().strftime("%H:%M:%S")
#         print(f"\x1b[1;32msort sol:\x1b[0m{tr}",end=" ")
        
#         # 3.3 Final step
#         sol_all, vol_all  = [], []
#         volAsym , sel1    = 0, -1
        
#         for jc, poly in enumerate(solun):
#             if xexp in poly:
                
#                 sol_all.append(poly)
#                 if not vol_all:
#                     sel1   = jc
#                     vol_all.append(poly.volume)
#                     volume = poly.volume
                    
#                 else:
#                     sel1   = sel1   if vol_all[-1] <= poly.volume else jc
#                     volume = volume if vol_all[-1] <= poly.volume else poly.volume 
#                     vol_all.append(poly.volume)
            
#             volAsym += pc.volume(poly)
        
#         print("L(sol)= ", len(solun)," s =",sel1) 
        
#         rr=((3/4)*(1/np.pi)*(volAsym))**(1/3)
        
#         sortorder = np.argsort(rp)[::-1]
#         localmat  = pc.extreme(solun[sel1])
        
#         for zi, zj in enumerate(localmat):
#             localmat[zi]=zj[sortorder.argsort()]
                
#         dmax = np.max(localmat, axis=0)
#         dmin = np.min(localmat, axis=0)
#         err  = np.abs(dmax-dmin)/2
#         wrtdata(fname, rc, pc.volume(solun[sel1]), err, solun[sel1], localmat, [volAsym, rr], len(solun))
        
#         twrite=datetime.now()
#         tinfo.append(twrite.timestamp()-tintersect.timestamp())
        
#         t2  = datetime.now()
#         tr=datetime.now().strftime("%H:%M:%S")
#         print(f"\x1b[1;32mwrt:\x1b[0m{tr}\x1b[1;32m Total: \x1b[0m{t2-t1}")
        
#         ttotal=datetime.now()
#         tinfo.append(ttotal.timestamp()-t1.timestamp())
        
#         wrttime_mc(rc, fname, tinfo)
        
#         m = np.mean(localmat,0)
#         print(f"---> The structure is predicted to be : {m} or {np.sort(0.5-m)[::-1]}")    
        
#     TE=datetime.now()    
#     print("Total: \x1b[0m{}".format(TE-TS0))
    
#     return
