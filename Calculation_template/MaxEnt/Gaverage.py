import numpy as np
import collections as cllct
import os
import os.path
import sys
import errno


# --------------------------------------------------------------------------- #

def path_exists(path):
    try:
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise

# --------------------------------------------------------------------------- #


El="Mo"
orbs=[1]
spins=[1]

#
dirs = [d for d in os.listdir('./') if os.path.isdir(os.path.join('./', d))]
dirs.remove("0")
if "avg" in dirs: dirs.remove("avg")
Nit = len(dirs)

#
taulat = np.genfromtxt('%s/Convergence/Glat/Glat_t_o%s_s%s.DAT'%(dirs[0],orbs[0],spins[0]), dtype='double', usecols=(0), unpack=True, comments='#')
tauimp = np.genfromtxt('%s/Convergence/Gqmc_%s/Gqmc_%s_t_o%s_s%s.DAT'%(dirs[0],El,El,orbs[0],spins[0]), dtype='double', usecols=(0), unpack=True, comments='#')
wm = np.genfromtxt('%s/Convergence/Wlat/Wlat_w_(%s,%s)(%s,%s).DAT'%(dirs[0],orbs[0],orbs[0],orbs[0],orbs[0]), dtype='double', usecols=(0), unpack=True, comments='#')

#
Glat = cllct.defaultdict(dict)
Wlat = cllct.defaultdict(dict)
Gqmc = cllct.defaultdict(dict)

idir = 1
dir = dirs[0]
print(dir,idir,Nit)
for iorb in orbs:
    for ispin in spins:
        Glat[iorb][ispin] = np.genfromtxt('%s/Convergence/Glat/Glat_t_o%s_s%s.DAT'%(dir,iorb,ispin), dtype='double', usecols=(1), unpack=True, comments='#')/Nit
        Gqmc[iorb][ispin] = np.genfromtxt('%s/Convergence/Gqmc_%s/Gqmc_%s_t_o%s_s%s.DAT'%(dir,El,El,iorb,ispin), dtype='double', usecols=(1), unpack=True, comments='#')/Nit
    Wlat[iorb] = np.genfromtxt('%s/Convergence/Wlat/Wlat_w_(%s,%s)(%s,%s).DAT'%(dir,iorb,iorb,iorb,iorb), dtype='double', usecols=(1), unpack=True, comments='#')/Nit
idir+=1

for dir in dirs[1:]:
    print(dir,idir,Nit)
    for iorb in orbs:
        for ispin in spins:
            Glat[iorb][ispin] += np.genfromtxt('%s/Convergence/Glat/Glat_t_o%s_s%s.DAT'%(dir,iorb,ispin), dtype='double', usecols=(1), unpack=True, comments='#')/Nit
            Gqmc[iorb][ispin] += np.genfromtxt('%s/Convergence/Gqmc_%s/Gqmc_%s_t_o%s_s%s.DAT'%(dir,El,El,iorb,ispin), dtype='double', usecols=(1), unpack=True, comments='#')/Nit
        Wlat[iorb] += np.genfromtxt('%s/Convergence/Wlat/Wlat_w_(%s,%s)(%s,%s).DAT'%(dir,iorb,iorb,iorb,iorb), dtype='double', usecols=(1), unpack=True, comments='#')/Nit
    idir+=1


path_exists("./avg/Convergence/Glat")
path_exists("./avg/Convergence/Wlat")
path_exists("./avg/Convergence/Gqmc_%s"%El)

for iorb in orbs:
    for ispin in spins:
        np.savetxt('./avg/Convergence/Glat/Glat_t_o%s_s%s.DAT'%(iorb,ispin), np.c_[ taulat, Glat[iorb][ispin] ], delimiter='\t')
        np.savetxt('./avg/Convergence/Gqmc_%s/Gqmc_%s_t_o%s_s%s.DAT'%(El,El,iorb,ispin), np.c_[ tauimp, Gqmc[iorb][ispin] ], delimiter='\t')
    np.savetxt('./avg/Convergence/Wlat/Wlat_w_(%s,%s)(%s,%s).DAT'%(iorb,iorb,iorb,iorb), np.c_[ wm, Wlat[iorb] ], delimiter='\t')





exit()
# --------------------------------------------------------------------------- #

flattenOrbndx=[o for subset in SiteList.values() for o in subset]

f_in=open('GWinput/Hk.dat','r')
dimes=f_in.readline().split()
N_Orb=int(dimes[2])
print('Norb=',N_Orb)
f_in.close()

tau   = np.genfromtxt('%s/GGWt_up/0_0.dat'%options.start, dtype='double', usecols=(0), unpack=True, comments='#')
wmatsG = np.genfromtxt('%s/GGWw_up/0_0.dat'%options.start, dtype='double', usecols=(0), unpack=True, comments='#')
wmatsW = np.genfromtxt('%s/WGWloc/(0, 0)_(0, 0).dat'%options.start, dtype='double', usecols=(0), unpack=True, comments='#')

GS = cllct.defaultdict(dict)
GGW = cllct.defaultdict(dict)
for ispin in ["up","dn"]:
    for axis in ["t","w"]:
        ndx="%s_%s"%(axis,ispin)
        l = len(tau) if axis=="t" else len(wmatsG)
        tp = np.float64 if axis=="t" else np.complex128
        GS[ndx]=np.zeros((N_Orb,N_Orb,l),dtype=tp)
        GGW[ndx]=np.zeros((N_Orb,N_Orb,l),dtype=tp)

impSw = cllct.defaultdict(dict)
impXC = cllct.defaultdict(dict)
for ispin in ["up","dn"]:
    l = len(wmatsG)
    impSw[ispin]=np.zeros((N_Orb,N_Orb,l),dtype=np.complex128)
    impXC[ispin]=np.zeros((N_Orb,N_Orb,l),dtype=np.complex128)

US=np.zeros((N_Orb,len(wmatsW)))
WS=np.zeros((N_Orb,len(wmatsW)))
WGW=np.zeros((N_Orb,len(wmatsW)))

Mu_vec={}
for a,v in SiteList.items():
    Mu_vec[a]=np.zeros(2*len(v))

HartreeU=np.zeros([N_Orb,N_Orb],dtype=numpy.float64)

rho_up=np.zeros([N_Orb,N_Orb],dtype=numpy.float64)
rho_dn=np.zeros([N_Orb,N_Orb],dtype=numpy.float64)

# --------------------------------------------------------------------------- #

for iteration in range(int(options.start),int(options.end)+1):
    #
    mess="it: %s"%iteration
    if str(iteration) in options.exceptl.split("_"):
        print(mess+" --> skipped")
        continue
    print(mess)
    # Gfs
    for ispin in ["up","dn"]:
        for axis in ["t","w"]:
            ndx="%s_%s"%(axis,ispin)
            for iorb in flattenOrbndx:
                for jorb in flattenOrbndx:
                    #
                    if axis=="w":
                        Re, Im = np.genfromtxt('%s/GS%s_%s/%s_%s.dat'%(iteration,axis,ispin,iorb,jorb) , dtype='double', usecols=(1,2), unpack=True, comments='#')
                        GS[ndx][iorb,jorb,:] += (Re+1j*Im)/Nfunct
                        Re, Im = np.genfromtxt('%s/GGW%s_%s/%s_%s.dat'%(iteration,axis,ispin,iorb,jorb) , dtype='double', usecols=(1,2), unpack=True, comments='#')
                        GGW[ndx][iorb,jorb,:] += (Re+1j*Im)/Nfunct
                    else:
                        GS[ndx][iorb,jorb,:] += np.genfromtxt('%s/GS%s_%s/%s_%s.dat'%(iteration,axis,ispin,iorb,jorb) , dtype='double', usecols=(1), unpack=True, comments='#')/Nfunct
                        GGW[ndx][iorb,jorb,:] += np.genfromtxt('%s/GGW%s_%s/%s_%s.dat'%(iteration,axis,ispin,iorb,jorb) , dtype='double', usecols=(1), unpack=True, comments='#')/Nfunct
    # self-energy
    for ispin in ["up","dn"]:
        for iorb in flattenOrbndx:
            for jorb in flattenOrbndx:
                #
                Re, Im = np.genfromtxt('%s/SigmaSw_%s/%s_%s.dat'%(iteration,ispin,iorb,jorb) , dtype='double', usecols=(1,2), unpack=True, comments='#')
                impSw[ispin][iorb,jorb,:] +=  (Re+1j*Im)/Nfunct
                Re, Im = np.genfromtxt('%s/SigmaXC_%s/%s_%s.dat'%(iteration,ispin,iorb,jorb) , dtype='double', usecols=(1,2), unpack=True, comments='#')
                impXC[ispin][iorb,jorb,:] +=  (Re+1j*Im)/Nfunct
    # interactions
    for iorb in flattenOrbndx:
        US[iorb,:] += np.genfromtxt('%s/USw/(%s, %s)_(%s, %s).dat'%(iteration,iorb,iorb,iorb,iorb) , dtype='double', usecols=(1), unpack=True, comments='#')/Nfunct
        WS[iorb,:] += np.genfromtxt('%s/WSw/(%s, %s)_(%s, %s).dat'%(iteration,iorb,iorb,iorb,iorb) , dtype='double', usecols=(1), unpack=True, comments='#')/Nfunct
        WGW[iorb,:] += np.genfromtxt('%s/WGWloc/(%s, %s)_(%s, %s).dat'%(iteration,iorb,iorb,iorb,iorb) , dtype='double', usecols=(1), unpack=True, comments='#')/Nfunct
    # local level
    for a,v in SiteList.items():
        f=open('%s/SolverOutputs_%s/MU_VECTOR'%(iteration,a),'r')
        muline=f.readline().split()
        for iorb in range(2*len(v)):
            Mu_vec[a][iorb] += float(muline[iorb])/Nfunct
        f.close()
    # Hartree
    HartreeU += np.loadtxt("%s/HartreeU"%iteration)/Nfunct
    # Gw density matrix
    rho_up += np.loadtxt("%s/GW_density_matrix_1.DAT"%iteration)/Nfunct
    rho_dn += np.loadtxt("%s/GW_density_matrix_2.DAT"%iteration)/Nfunct

# --------------------------------------------------------------------------- #

for axis in ["t","w"]:
    l = len(tau) if axis=="t" else len(wmatsG)
    for ispin in ["up","dn"]:
        ndx="%s_%s"%(axis,ispin)
        for a,v in G_equivalent.items():
            GSeq=np.zeros(l)
            GGWeq=np.zeros(l)
            Sigmaeq=np.zeros(l)
            for iorb in v:
                GSeq+=GS[ndx][iorb,iorb,:]/len(v)
                GGWeq+=GGW[ndx][iorb,iorb,:]/len(v)
                if axis=="w":Sigmaeq+=impSw[ispin][iorb,iorb,:]/len(v)
        for a,v in G_equivalent.items():
            for iorb in v:
                GS[ndx][iorb,iorb,:]=GSeq.copy()
                GGW[ndx][iorb,iorb,:]=GGWeq.copy()
                if axis=="w":impSw[ispin][iorb,iorb,:]=Sigmaeq.copy()


for a,v in W_equivalent.items():
    USeq=np.zeros(len(wmatsW))
    WSeq=np.zeros(len(wmatsW))
    WGWeq=np.zeros(len(wmatsW))
    for iorb in v:
        USeq+=US[iorb,:]/len(v)
        WSeq+=WS[iorb,:]/len(v)
        WGWeq+=WGW[iorb,:]/len(v)
for a,v in W_equivalent.items():
    for iorb in v:
        US[iorb,:]=USeq.copy()
        WS[iorb,:]=WSeq.copy()
        WGW[iorb,:]=WGWeq.copy()

# --------------------------------------------------------------------------- #

if(paramagnet):
    for iorb in flattenOrbndx:
        for jorb in flattenOrbndx:
            for axis in ["t","w"]:
                GS["%s_%s"%(axis,"up")][iorb,jorb,:]=0.5*(GS["%s_%s"%(axis,"up")][iorb,jorb,:]+GS["%s_%s"%(axis,"dn")][iorb,jorb,:])
                GS["%s_%s"%(axis,"dn")][iorb,jorb,:]=GS["%s_%s"%(axis,"up")][iorb,jorb,:]
                GGW["%s_%s"%(axis,"up")][iorb,jorb,:]=0.5*(GGW["%s_%s"%(axis,"up")][iorb,jorb,:]+GGW["%s_%s"%(axis,"dn")][iorb,jorb,:])
                GGW["%s_%s"%(axis,"dn")][iorb,jorb,:]=GGW["%s_%s"%(axis,"up")][iorb,jorb,:]
            impSw["up"][iorb,jorb,:]=0.5*(impSw["up"][iorb,jorb,:]+impSw["dn"][iorb,jorb,:])
            impSw["dn"][iorb,jorb,:]=impSw["up"][iorb,jorb,:]
            impXC["up"][iorb,jorb,:]=0.5*(impXC["up"][iorb,jorb,:]+impXC["dn"][iorb,jorb,:])
            impXC["dn"][iorb,jorb,:]=impXC["up"][iorb,jorb,:]

# --------------------------------------------------------------------------- #

print("Last iteration: ",iteration)
outpath = str(startit)+'_'+str(endit)
path_exists(outpath)
# Gfs
for ispin in ["up","dn"]:
    for axis in ["t","w"]:
        absis = tau if axis=="t" else wmatsG
        path_exists('%s/GS%s_%s'%(outpath,axis,ispin))
        path_exists('%s/GGW%s_%s'%(outpath,axis,ispin))
        for iorb in flattenOrbndx:
            for jorb in flattenOrbndx:
                if axis=="w":
                    np.savetxt('%s/GS%s_%s/%s_%s.dat'%(outpath,axis,ispin,iorb,jorb),  np.c_[absis,GS["%s_%s"%(axis,ispin)][iorb,jorb,:].real,GS["%s_%s"%(axis,ispin)][iorb,jorb,:].imag], delimiter='\t')
                    np.savetxt('%s/GGW%s_%s/%s_%s.dat'%(outpath,axis,ispin,iorb,jorb), np.c_[absis,GGW["%s_%s"%(axis,ispin)][iorb,jorb,:].real,GGW["%s_%s"%(axis,ispin)][iorb,jorb,:].imag], delimiter='\t')
                else:
                    np.savetxt('%s/GS%s_%s/%s_%s.dat'%(outpath,axis,ispin,iorb,jorb),  np.c_[absis,GS["%s_%s"%(axis,ispin)][iorb,jorb,:].real], delimiter='\t')
                    np.savetxt('%s/GGW%s_%s/%s_%s.dat'%(outpath,axis,ispin,iorb,jorb), np.c_[absis,GGW["%s_%s"%(axis,ispin)][iorb,jorb,:].real], delimiter='\t')
# self-energy
for ispin in ["up","dn"]:
    path_exists('%s/SigmaSw_%s'%(outpath,ispin))
    path_exists('%s/SigmaXC_%s'%(outpath,ispin))
    for iorb in flattenOrbndx:
        for jorb in flattenOrbndx:
            np.savetxt('%s/SigmaSw_%s/%s_%s.dat'%(outpath,ispin,iorb,jorb), np.c_[wmatsG,impSw[ispin][iorb,jorb,:].real,impSw[ispin][iorb,jorb,:].imag], delimiter='\t')
            np.savetxt('%s/SigmaXC_%s/%s_%s.dat'%(outpath,ispin,iorb,jorb), np.c_[wmatsG,impXC[ispin][iorb,jorb,:].real,impXC[ispin][iorb,jorb,:].imag], delimiter='\t')
# interactions
path_exists('%s/USw'%outpath)
path_exists('%s/WSw'%outpath)
path_exists('%s/WGWloc'%outpath)
for a,v in SiteList.items():
    for iorb in v:
        np.savetxt('%s/USw/(%s, %s)_(%s, %s).dat'%(outpath,iorb,iorb,iorb,iorb), np.c_[wmatsW,US[iorb,:]], delimiter='\t')
        np.savetxt('%s/WSw/(%s, %s)_(%s, %s).dat'%(outpath,iorb,iorb,iorb,iorb), np.c_[wmatsW,WS[iorb,:]], delimiter='\t')
        np.savetxt('%s/WGWloc/(%s, %s)_(%s, %s).dat'%(outpath,iorb,iorb,iorb,iorb), np.c_[wmatsW,WGW[iorb,:]], delimiter='\t')
# local level
for a,v in SiteList.items():
    path_exists('%s/SolverOutputs_%s'%(outpath,a))
    f=open('%s/SolverOutputs_%s/MU_VECTOR'%(outpath,a),'w')
    for iorb in range(2*len(v)):
        f.write("%s "%Mu_vec[a][iorb])
    f.close()
# Hartree
np.savetxt("%s/HartreeU"%outpath,HartreeU)
# Gw density matrix
np.savetxt("%s/GW_density_matrix_1.DAT"%outpath,rho_up)
np.savetxt("%s/GW_density_matrix_2.DAT"%outpath,rho_dn)


exit()
