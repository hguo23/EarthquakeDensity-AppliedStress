"""
Created on Wed Apr 20 11:49:23 2022

@author: hguo23
"""

from joblib import Parallel, delayed
import multiprocessing

from okada_wrapper import dc3d0wrapper, dc3dwrapper
import numpy as np
import scipy.io
from math import sin, cos, sqrt, atan2, radians
from random import random, seed

from joblib import Parallel, delayed
from tqdm import tqdm
from scipy.io import savemat

def get_params():  
    mu = 0.5 #coefficient of friction
    poisson_ratio = 0.25
    G = 3e10
    lmda = 2 * G * poisson_ratio / (1 - 2 * poisson_ratio)
    alpha = (lmda + G) / (lmda + 2 * G)
    return mu, poisson_ratio, G, lmda, alpha

def get_distance(lat1, lon1, lat2, lon2):    
    # approximate radius of earth in m
    R = 6378137
    lat1 = radians(lat1)
    lat2 = radians(lat2)
    lon1 = radians(lon1)
    lon2 = radians(lon2)
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = sin(dlat / 2)**2 + cos(lat1) * cos(lat2) * sin(dlon / 2)**2
    c = 2 * atan2(sqrt(a), sqrt(1 - a))
    return R * c


def Stresses(S, dip):
    cd = cos(radians(dip))
    sd = sin(radians(dip))
    normal = [sd, 0, cd]
    along = [cd, 0, -sd]
    
    S = np.array([S[0], S[1], S[2], S[1], S[3], S[4], S[2], S[4], S[5]]) 
    nml = np.array([normal[0], normal[1], normal[2], normal[0], normal[1],
                    normal[2], normal[0], normal[1], normal[2]])
    acons = S*nml
    acons = [acons[0]+acons[1]+acons[2], acons[3]+acons[4]+acons[5], acons[6]+acons[7]+acons[8]]         
    shr_stress = (np.array(along) * np.array(acons)).sum()
    nml_stress = (np.array(normal) * np.array(acons)).sum()
    return shr_stress, nml_stress

def CoordVec(strike, dip):
    ss = sin(radians(strike))
    cs = cos(radians(strike))
    sd = sin(radians(dip))
    cd = cos(radians(dip))
    nn = np.array([-ss*sd, cs*sd, -cd])
    ns = np.array([cs, ss, 0])
    nd = np.array([-ss*cd, cs*cd, sd])
    return nn, ns, nd

def StressesM(S, strike, dip):
    nn, ns, nd = CoordVec(strike, dip)
    t = np.array([(S[0]*nn).sum(),(S[1]*nn).sum(),(S[2]*nn).sum()])
    tn = (t*nn).sum()
    td = (t*nd).sum()
    ts = (t*ns).sum()
    tau = np.sqrt(td**2 + ts**2)
    return tau, tn

def dc3d(loopnum, strikevec1, dipvec1, rakevec1, strikevec2, dipvec2, rakevec2, MAG,lat0, lon0, dep0, lat, lon, dep):
    outfile = '../output/qtm_2D_CSC{}.mat'.format(str(int(loopnum)))
    N = len(MAG)
    
    mu, poisson_ratio, G, lmda, alpha = get_params()
    errxout = []
    erryout = []
    errzout = []
    stc1 = []
    stc2 = []
    R1 = []
    R2 = []
   
    for n in range(0, N):
        if (dep[n][0]).size:
            Nd = len(dep[n][0][0])
            s1 = float(strikevec1[n])
            d1 = float(dipvec1[n])
            r1 = float(rakevec1[n])
            s2 = float(strikevec2[n])
            d2 = float(dipvec2[n])
            r2 = float(rakevec2[n])
            
            M0 = 10**(float(MAG[n])*1.5+9.1)
            l  = 10*10**(0.44*float(MAG[n]))
            w = l
            slip = M0/(w*l*G)
                        
            ss1 = slip*cos(radians(r1))
            ds1 = slip*sin(radians(r1))
            ss2 = slip*cos(radians(r2))
            ds2 = slip*sin(radians(r2))
            
            al1 = l*np.random.rand(1)
            al2 = l - al1
            aw1 = w*np.random.rand(1)
            aw2 = w - aw1
            
            slat = float(lat0[n])
            slon = float(lon0[n])
            sdep = float(dep0[n])
            
            errh0 = 0.1343*3*1e3*np.random.randn(1)
            errz0 = 0.3464*3*np.random.randn(1)
            theta0 = 2*np.pi*np.random.randn(1)
            errx0 = errh0*cos(theta0)
            erry0 = sqrt(errh0**2-errx0**2)
                
            errxout.append(-errh0*sin(theta0-radians(s1))/1e3)
            erryout.append(errh0*cos(theta0-radians(s1))/1e3)
            errzout.append(errz0)

            for nd in range(Nd):
                alat = lat[n][0][0][nd]
                alon = lon[n][0][0][nd]
                adep = dep[n][0][0][nd]
                x = np.sign(alon-slon)*get_distance(slat, slon, slat, alon)
                y = np.sign(alat-slat)*get_distance(slat, slon, alat, slon)
                dis = float(np.sqrt(x**2+y**2))
                theta1 = atan2(y,x)-(np.pi-radians(s1))
                theta2 = atan2(y,x)-(np.pi-radians(s2))
                x1 = dis*cos(theta1)
                y1 = dis*sin(theta1)
                x2 = dis*cos(theta2)
                y2 = dis*sin(theta2)

                errh = 0.1343*3*1e3*np.random.randn(1)
                theta1 = 2*np.pi*np.random.randn(1)
                errz1 = 0.3464*3*np.random.randn(1)
                errx1 = errh*cos(theta1)
                erry1 = sqrt(errh**2-errx1**2)
                
                errz2 = 0.3464*3*np.random.randn(1)
                theta2 = 2*np.pi*np.random.randn(1)
                errx2 = errh*cos(theta2)
                erry2 = sqrt(errh**2-errx2**2)
                
                x1 = x1 + errx1 - errx0
                y1 = y1 + erry1 - erry0
                x2 = x2 + errx2 - errx0
                y2 = y2 + erry2 - erry0

                dis1 = sqrt(x1**2 + y1**2)
                dis3D1 = sqrt((dis1/1e3)**2+(adep-sdep+errz1-errz0)**2)
                dis2 = sqrt(x2**2 + y2**2)
                dis3D2 = sqrt((dis2/1e3)**2+(adep-sdep+errz2-errz0)**2)
                R1.append(np.round(dis3D1, 6))
                R2.append(np.round(dis3D2, 6))             
                                 
                success1, disp1, grad_u1 = dc3dwrapper(alpha, [x1, y1, -(adep+errz1)*1e3],
                                                    (sdep+errz0)*1e3, d1, [-al1, al2], [-aw1, aw2],
                                                   [ss1, ds1, 0.0])
                success2, disp2, grad_u2 = dc3dwrapper(alpha, [x2, y2, -(adep+errz2)*1e3],
                                                    (sdep+errz0)*1e3, d2, [-al1, al2], [-aw1, aw2],
                                                   [ss2, ds2, 0.0])

                if success1 == 0 and success2 == 0:
                    ud1 = np.array(grad_u1)
                    ud1 = np.reshape(ud1,[3, 3])
                    udtmp = np.transpose(ud1)
                    ud1 = (ud1 + udtmp)/2 #symmetrize
                    stress1 = 2*G*ud1 + lmda*(np.diagonal(ud1)*np.eye(3)).sum()
                    stmp1 = 90 + s1 - s1
                    stmp2 = 90 + s2 - s1
                    tau1, tn1 = StressesM(stress1, stmp1, d1)
                    tau2, tn2 = StressesM(stress1, stmp2, d2)
                    
                    ud2 = np.array(grad_u2)
                    ud2 = np.reshape(ud2,[3, 3])
                    udtmp = np.transpose(ud2)
                    ud2 = (ud2 + udtmp)/2 #symmetrize
                    stress2 = 2*G*ud2 + lmda*(np.diagonal(ud2)*np.eye(3)).sum()
                    stmp1 = 90 + s1 - s2
                    stmp2 = 90 + s2 - s2
                    tau3, tn3 = StressesM(stress2, stmp1, d1)
                    tau4, tn4 = StressesM(stress2, stmp2, d2)
                    stc1.append(np.round(tau1+mu*tn1, 4))
                    stc1.append(np.round(tau2+mu*tn2, 4))
                    stc2.append(np.round(tau3+mu*tn3, 4))
                    stc2.append(np.round(tau4+mu*tn4, 4))
                    
                else:
                    stc1.append(0)
                    stc1.append(0)
                    stc2.append(0)
                    stc2.append(0)    
    stc1 = np.array(stc1)
    stc2 = np.array(stc2)
    R1 = np.array(R1)
    R2 = np.array(R2)
    errxout = np.array(errxout)
    erryout = np.array(erryout)
    errzout = np.array(errzout)
    stas = {"Errx": errxout, "Erry": erryout, "Errz": errzout,
            "R1": R1, "R2": R2, "stc1": stc1, "stc2": stc2}
    savemat(outfile, stas)    
    
def importfile():
    mat = scipy.io.loadmat('../output/qtm_2D.mat')
    N = mat['Rec']['MainDep'].shape[1]
    dep0 = np.array_split(mat['Rec']['MainDep'][0],N)
    lat0 = np.array_split(mat['Rec']['MainLat'][0],N)
    lon0 = np.array_split(mat['Rec']['MainLon'][0],N)
    MAG = np.array_split(mat['Rec']['MainMag'][0],N)
    lat = np.array_split(mat['Rec']['Lat'][0],N)
    lon = np.array_split(mat['Rec']['Lon'][0],N)
    dep = np.array_split(mat['Rec']['Dep'][0],N)
    strikevec1 = np.array_split(mat['Rec']['s1'][0],N)
    dipvec1 = np.array_split(mat['Rec']['d1'][0],N)
    rakevec1 = np.array_split(mat['Rec']['r1'][0],N)
    strikevec2 = np.array_split(mat['Rec']['s2'][0],N)
    dipvec2 = np.array_split(mat['Rec']['d2'][0],N)
    rakevec2 = np.array_split(mat['Rec']['r2'][0],N)
    return MAG, lat0, lon0, dep0, lat, lon, dep, strikevec1, dipvec1, rakevec1, strikevec2, dipvec2, rakevec2

if __name__ == '__main__':
    num_cores = 30
    Nboot = 5e2
    MAG, lat0, lon0, dep0,lat, lon, dep, strikevec1, dipvec1, rakevec1, strikevec2, dipvec2, rakevec2 = importfile()

#    dc3d(1, strikevec1, dipvec1, rakevec1,strikevec2, dipvec2, rakevec2, MAG,lat0, lon0, dep0, lat, lon, dep)

    Parallel(n_jobs = num_cores)(delayed(dc3d)(loopnum, strikevec1, dipvec1, rakevec1, 
                                               strikevec2, dipvec2, rakevec2, MAG,
                                                lat0, lon0, dep0, lat, lon, dep) 
                                 for loopnum in tqdm(np.arange(Nboot)))


"""
Created on Mon Apr 18 11:13:56 2022

@author: hguo23
"""
