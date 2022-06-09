import numpy as np
import collections as cllct
import os
import os.path
import sys
import errno
import glob
from optparse import OptionParser


# --------------------------------------------------------------------------- #

def path_exists(path):
    try:
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise

def average(G,tlimit,window,log):
    Gin = np.log(np.abs(G-(1e-15))) if log else G.copy()
    Ntau = np.shape(G)[0]; Norb = np.shape(G)[1]
    tbound = tlimit + int(window/2)
    for iorb in range(1,Norb):
        Gout = np.convolve(Gin[tlimit:Ntau-1-tlimit,iorb], np.ones(window), 'valid') / window
        G[tbound:tbound+len(Gout),iorb] = -np.exp(Gout) if log else Gout.copy()

# --------------------------------------------------------------------------- #

# Default arguments
spins = [1]; orbs = [1]; dirs=[]
# Optional parsing
parser = OptionParser()
parser.add_option("--El"    , dest="El"     , default="Mo" )
parser.add_option("--orbs"  , dest="orbs"   , default="1"  )
parser.add_option("--itlist", dest="itlist" , default="0"  )
parser.add_option("--mag"   , dest="mag", default=False, action="store_true")
parser.add_option("--Akw"   , dest="Akw", default=False, action="store_true")
parser.add_option("--avg"   , dest="avg", default="0" )
parser.add_option("--log"   , dest="log", default=False, action="store_true")
parser.add_option("--dlt"   , dest="dlt", default=False, action="store_true")
(options, args) = parser.parse_args()
#
El = options.El
spins = [1,2] if options.mag else [1]
orbs = [ int(o) for o in options.orbs.split(",") ] if options.orbs != "1" else [1]
#
alldirs = [d for d in os.listdir('./') if os.path.isdir(os.path.join('./', d))]
if "0" in alldirs: alldirs.remove("0")
if "avg" in alldirs: alldirs.remove("avg")
if options.itlist != "0":
    dirs = [ i for i in alldirs if i in options.itlist.split(",") ]
else:
    dirs = alldirs
#
if options.Akw:
    removelist=[]
    for dir in dirs:
        if not os.path.isdir(dir+"/K_resolved/MaxEnt_Gk_path_t_s1"): removelist+=[dir]
    dirs = [ i for i in dirs if i not in removelist ]
#
Nit = len(dirs)
#
print("El=%s, spins=%s, orbs=%s"%(El,spins,orbs))
print("itlist=%s, Nit=%s"%(dirs,Nit))

# --------------------------------------------------------------------------- #
# Here I choose to either perform the avergave on the local observables
# OR on the K-resolved ones. They're too fdifferent in kind
# --------------------------------------------------------------------------- #
if options.Akw:
    #
    # Check for Nkpt consistency between iterations
    Gk = glob.glob(dirs[0]+"/K_resolved/MaxEnt_Gk_path_t_s1/Gk_t_k*.DAT")
    Nkpt = len(Gk)
    print(Nkpt)
    for dir in dirs:
        Gk = glob.glob(dir+"/K_resolved/MaxEnt_Gk_path_t_s1/Gk_t_k*.DAT")
        print("Check for Nkpt in folder: %s ,Nkpt=%s "%(dir,len(Gk)))
        if len(Gk) != Nkpt: sys.exit("Wrong number of K-points")
    #
    path_exists("./avg/K_resolved/MaxEnt_Gk_path_t_s1")
    #
    for ik in range(1,Nkpt+1):
        dir = dirs[0]
        Gpath = dir+"/K_resolved/MaxEnt_Gk_path_t_s1/Gk_t_k%s.DAT"%ik
        Gtauk = np.loadtxt(Gpath)
        if options.avg !="0": average(Gtauk,int(options.avg.split(",")[0]),int(options.avg.split(",")[1]),options.log)
        Gtau = Gtauk.copy()/Nit
        for dir in dirs[1:]:
            Gpath = dir+"/K_resolved/MaxEnt_Gk_path_t_s1/Gk_t_k%s.DAT"%ik
            Gtauk = np.loadtxt(Gpath)
            if options.avg !="0": average(Gtauk,int(options.avg.split(",")[0]),int(options.avg.split(",")[1]),options.log)
            Gtau += Gtauk/Nit
        if options.dlt: Gtau = np.delete(Gtau, slice(None, None, 2), 0)
        np.savetxt('./avg/K_resolved/MaxEnt_Gk_path_t_s1/Gk_t_k%s.DAT'%ik, Gtau )
        print("Akw average ik=%s, Nkpt=%s"%(ik,Nkpt))
else:
    #
    # Default meshes
    taulat = np.genfromtxt('%s/Convergence/Glat/Glat_t_o%s_s%s.DAT'%(dirs[0],orbs[0],spins[0]), dtype='double', usecols=(0), unpack=True, comments='#')
    tauimp = np.genfromtxt('%s/Convergence/Gqmc_%s/Gqmc_%s_t_o%s_s%s.DAT'%(dirs[0],El,El,orbs[0],spins[0]), dtype='double', usecols=(0), unpack=True, comments='#')
    wm = np.genfromtxt('%s/Convergence/Wlat/Wlat_w_(%s,%s)(%s,%s).DAT'%(dirs[0],orbs[0],orbs[0],orbs[0],orbs[0]), dtype='double', usecols=(0), unpack=True, comments='#')
    #
    # Dictionaries
    Glat = cllct.defaultdict(dict)
    Wlat = cllct.defaultdict(dict)
    Gqmc = cllct.defaultdict(dict)
    #
    # Data initialization
    idir = 1
    dir = dirs[0]
    print(dir,idir,Nit)
    for iorb in orbs:
        for ispin in spins:
            Glat[iorb][ispin] = np.genfromtxt('%s/Convergence/Glat/Glat_t_o%s_s%s.DAT'%(dir,iorb,ispin), dtype='double', usecols=(1), unpack=True, comments='#')/Nit
            Gqmc[iorb][ispin] = np.genfromtxt('%s/Convergence/Gqmc_%s/Gqmc_%s_t_o%s_s%s.DAT'%(dir,El,El,iorb,ispin), dtype='double', usecols=(1), unpack=True, comments='#')/Nit
        Wlat[iorb] = np.genfromtxt('%s/Convergence/Wlat/Wlat_w_(%s,%s)(%s,%s).DAT'%(dir,iorb,iorb,iorb,iorb), dtype='double', usecols=(1), unpack=True, comments='#')/Nit
    idir+=1
    #
    # Adding up to average
    for dir in dirs[1:]:
        print("local average ",dir,idir,Nit)
        for iorb in orbs:
            for ispin in spins:
                Glat[iorb][ispin] += np.genfromtxt('%s/Convergence/Glat/Glat_t_o%s_s%s.DAT'%(dir,iorb,ispin), dtype='double', usecols=(1), unpack=True, comments='#')/Nit
                Gqmc[iorb][ispin] += np.genfromtxt('%s/Convergence/Gqmc_%s/Gqmc_%s_t_o%s_s%s.DAT'%(dir,El,El,iorb,ispin), dtype='double', usecols=(1), unpack=True, comments='#')/Nit
            Wlat[iorb] += np.genfromtxt('%s/Convergence/Wlat/Wlat_w_(%s,%s)(%s,%s).DAT'%(dir,iorb,iorb,iorb,iorb), dtype='double', usecols=(1), unpack=True, comments='#')/Nit
        idir+=1
    #
    # Print output
    path_exists("./avg/Convergence/Glat")
    path_exists("./avg/Convergence/Wlat")
    path_exists("./avg/Convergence/Gqmc_%s"%El)
    for iorb in orbs:
        for ispin in spins:
            np.savetxt('./avg/Convergence/Glat/Glat_t_o%s_s%s.DAT'%(iorb,ispin), np.c_[ taulat, Glat[iorb][ispin] ], delimiter='\t')
            np.savetxt('./avg/Convergence/Gqmc_%s/Gqmc_%s_t_o%s_s%s.DAT'%(El,El,iorb,ispin), np.c_[ tauimp, Gqmc[iorb][ispin] ], delimiter='\t')
        np.savetxt('./avg/Convergence/Wlat/Wlat_w_(%s,%s)(%s,%s).DAT'%(iorb,iorb,iorb,iorb), np.c_[ wm, Wlat[iorb] ], delimiter='\t')
