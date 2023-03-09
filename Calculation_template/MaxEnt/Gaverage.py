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

def read_Convergence(path,):
    data = np.genfromtxt(path, dtype='double', usecols=(1), unpack=True, comments='#')
    return data

def average(G,tlimit,window,log,orblist=None):
    Gin = np.log(np.abs(G-(1e-15))) if log else G.copy()
    Ntau = np.shape(G)[0]; Norb = np.shape(G)[1]
    tbound = tlimit + int(window/2)
    orblist_ = range(Norb) if orblist==None else [ o-1 for o in orblist]
    print("tlimit=%s, window=%s, log=%s, orblist=%s"%(tlimit,window,log,orblist_))
    for iorb in orblist_:
        Gout = np.convolve(Gin[tlimit:Ntau-1-tlimit,iorb], np.ones(window), 'valid') / window
        G[tbound:tbound+len(Gout),iorb] = -np.exp(Gout) if log else Gout.copy()

# ================================= EXAMPLES ================================= #
#
# to print observable average:
# python Gaverage.py --field Obs --itlist it_start..it_end
#
# to print Glat average:
# python Gaverage.py --field Gw --pad lat --orbs 1,2,3,.. --itlist it_start..it_end
#
# to print Gqmc_El average:
# python Gaverage.py --field Gw --pad qmc_El --orbs 1,2,3,.. --itlist it_start..it_end
#
# to print Wlat average:
# python Gaverage.py --field Ww --pad lat --orbs 1,2,3,.. --itlist it_start..it_end
#
# to print curlyUimp average:
# python Gaverage.py --field curlyUw --pad imp --orbs 1,2,3,.. --itlist it_start..it_end
#
# to print Gk average:
# python Gaverage.py --field Gkw --itlist it_start..it_endv
#
# ============================================================================ #

# Default arguments
spins = [1]; orbs = [1]; dirs=[]
# Optional parsing
parser = OptionParser()
parser.add_option("--field" , dest="field"  , default="Gw"    )                 # [G,W]kw, [G,W]w, Obs
parser.add_option("--pad"   , dest="pad"    , default="None"  )                 # local: lat, imp, qmc_* kresolved: _Tr, _Hetero
parser.add_option("--orbs"  , dest="orbs"   , default="1"     )                 # list of orbitals separated by comma
parser.add_option("--itlist", dest="itlist" , default="0"     )                 # list of iterations separated by comma or initial and final iteration separated by ".."
parser.add_option("--skip"  , dest="skip"   , default="0"     )                 # list of iteration to be ignored separated by comma
parser.add_option("--spin"  , dest="spin"   , default="1"     )                 # spin index
parser.add_option("--Krol"  , dest="Krol"   , default="None"  )                 # string of 3 parameters separated by comma in dicating tau_max and tau_window for rolling average, if last one is ==1 then the rolling average is done in log scale
(options, args) = parser.parse_args()
#
#
if options.itlist != "None":
    #
    checklist = options.itlist
    checkrang = options.itlist
    #
    if len(checklist.split(","))>1:
        dirs = [ str(i) for i in options.itlist.split(",") if os.path.isdir(str(i)) ]
    elif len(checkrang.split(".."))==2:
        dirs = [ str(i) for i in range( int(options.itlist.split("..")[0]) , int(options.itlist.split("..")[1])+1 ) if os.path.isdir(str(i)) ]
    else:
        sys.exit("--itlist error. Exiting.")
else:
    sys.exit("--itlist error. Exiting.")
#
#
if options.skip != "0":
    skip_dirs = list(options.skip.split(","))
    print("skip_dirs: %s"%skip_dirs)
    for dir in skip_dirs:
        dirs.remove(dir)
#
#
if options.orbs != "1":
    checklist = options.orbs
    checkrang = options.orbs
    #check if the user provided the list
    if len(checklist.split(","))>1:
        orbs = [ int(o) for o in options.orbs.split(",") ]
    elif len(checkrang.split(".."))==2:
        orbs = range( int(options.orbs.split("..")[0]) , int(options.orbs.split("..")[1])+1 )
    else:
        sys.exit("--orbs error. Exiting.")
else:
    orbs = [1]
#
#
spin = int(options.spin)
#
#
kresolved = "k" in options.field
local =  not kresolved
#
#
Field = options.field.rstrip("kw")
#
#
if(Field=="G" or Field=="S"): stat = "Fermion"
if(Field=="W" or Field=="curlyU" or Field=="C" or Field=="M"): stat = "Boson"
if(Field=="Obs"): stat = "Observables"
#
Kfolder = "K_resolved"
Kfolder_out = "K_resolved"
#
#
# --------------------------------------------------------------------------- #
#                                  RECAP                                      #
# --------------------------------------------------------------------------- #
Nit = len(dirs)
if Nit==0: sys.exit("No folder found. Exiting.")
print("itlist=%s, Nit=%s"%(dirs,Nit))
print("Field=%s (%s), spin=%s, orbs=%s, kresolved=%s"%(Field,stat,spin,orbs,kresolved))
#
#
#
if (stat == "Observables"):
    #
    path_exists("./avg/")
    #
    # OBSERVABLES: Zqpsc, Zdmft full orbital basis
    dir = dirs[0]
    Zmat = cllct.defaultdict(dict)
    for st in ['dmft','qpsc']:
        File = '%s/Z_%s_s%s.DAT'%(dir,st,spin)
        if os.path.isfile(File): Zmat[st] = np.loadtxt(File)/Nit
    for dir in dirs[1:]:
        for st in ['dmft','qpsc']:
            File = '%s/Z_%s_s%s.DAT'%(dir,st,spin)
            if os.path.isfile(File): Zmat[st] += np.loadtxt(File)/Nit
    for st in ['dmft','qpsc']:np.savetxt('./avg/Z_%s_s%s.DAT'%(st,spin), Zmat[st])
    #
    #
    # OBSERVABLES: Nlat, Nimp full orbital basis
    dir = dirs[0]
    rho = cllct.defaultdict(dict)
    for st in ['imp','lat']:
        rho[st] = np.loadtxt('%s/N%s_s%s.DAT'%(dir,st,spin))/Nit
    for dir in dirs[1:]:
        for st in ['imp','lat']:
            rho[st] += np.loadtxt('%s/N%s_s%s.DAT'%(dir,st,spin))/Nit
    for st in ['imp','lat']:np.savetxt('./avg/N%s_s%s.DAT'%(st,spin), rho[st])
    #
    #
    # OBSERVABLES: matrices in Solver_*
    dir = dirs[0]
    SolverDirs = [ S.lstrip("%s/"%dir) for S in glob.glob("%s/Solver_*"%dir)]
    print("SolverDirs: %s"%SolverDirs)
    Umat = cllct.defaultdict(dict)
    Eloc = cllct.defaultdict(dict)
    for st in SolverDirs:
        Umat[st] = np.loadtxt('%s/%s/Umat.DAT'%(dir,st))/Nit
        Eloc[st] = np.loadtxt('%s/%s/Eloc.DAT'%(dir,st))/Nit
    for dir in dirs[1:]:
        for st in SolverDirs:
            Umat[st] += np.loadtxt('%s/%s/Umat.DAT'%(dir,st))/Nit
            Eloc[st] += np.loadtxt('%s/%s/Eloc.DAT'%(dir,st))/Nit
    for st in SolverDirs:
        np.savetxt('./avg/Umat_%s.DAT'%st.lstrip("Solver_"), Umat[st])
        np.savetxt('./avg/Eloc_%s.DAT'%st.lstrip("Solver_"), Eloc[st])
#
#
#
if (stat != "Observables") and local:
    #
    if options.pad=="None": sys.exit("--pad unspecified. Exiting.")
    #
    Folder = Field + options.pad
    path_exists("./avg/Convergence/%s/"%Folder)
    #
    for o in orbs:
        #
        File = Folder + "_t_o%s"%o + "_s%s.DAT"%spin if (stat == "Fermion") else Folder + "_w_(%s,%s)(%s,%s).DAT"%(o,o,o,o)
        print("Averaging file: %s"%File)
        #
        dir = dirs[0]
        #
        path = '%s/Convergence/%s/%s'%(dir,Folder,File)
        axis = np.genfromtxt(path, dtype='double', usecols=(0), unpack=True, comments='#')
        data = read_Convergence(path)/Nit
        for dir in dirs[1:]:
            path = '%s/Convergence/%s/%s'%(dir,Folder,File)
            data += read_Convergence(path)/Nit
        #
        np.savetxt('./avg/Convergence/%s/%s'%(Folder,File), np.c_[ axis, data ] , delimiter='\t')
        data*=0.0
#
#
#
if (stat != "Observables") and kresolved:
    #
    pad = "" if options.pad == "None" else options.pad
    Folder = "MaxEnt_%sk_path_s%s"%(Field,spin) if (stat == "Fermion") else "MaxEnt_%sk_path%s"%Field
    path_exists("./avg/%s/%s/"%(Kfolder_out,Folder))
    #
    # Check for Nkpt consistency between iterations
    File = Field + "k_%s_k"%("t" if (stat == "Fermion") else "w" )
    path = "./%s/%s/%s"%(dirs[0],Kfolder,Folder)
    Nkpt = 0
    for ik in range(1,1000):
        if os.path.exists("%s/%s%s%s.DAT"%(path,File,ik,pad)): Nkpt += 1
    print("File: %s, pad: %s, Nkpt: %s"%(File,pad,Nkpt))
    for dir in dirs:
        path = "./%s/%s/%s"%(dir,Kfolder,Folder)
        Nkpt_ = 0
        for ik in range(1,1000):
            if os.path.exists("%s/%s%s%s.DAT"%(path,File,ik,pad)): Nkpt_ += 1
        if(Nkpt_!=Nkpt): sys.exit("Wrong number of K-points")
    #
    path_exists("./avg/%s/%s/"%(Kfolder_out,Folder))
    #
    # Averaging one K-point at a time
    for ik in range(1,Nkpt+1):
        #
        Filek = "%s%s%s.DAT"%(File,ik,pad)
        #
        #initialization
        dir = dirs[0]
        #
        path = "./%s/%s/%s/%s"%(dir,Kfolder,Folder,Filek)
        axis = np.loadtxt(path)[:,0]
        Datak_read = np.loadtxt(path)[:,1:]
        if options.Krol != "None":
            print("Rolling average dir=%s, ik=%s"%(dir,ik))
            average(Datak_read, int(options.Krol.split(",")[0]), int(options.Krol.split(",")[1]), int(options.Krol.split(",")[1])==1, orbs if orbs!=[1] else None )
        Datak = Datak_read.copy()/Nit
        #
        for dir in dirs[1:]:
            #
            path = "./%s/%s/%s/%s"%(dir,Kfolder,Folder,Filek)
            Datak_read = np.loadtxt(path)[:,1:]
            if options.Krol != "None":
                print("Rolling average dir=%s, ik=%s"%(dir,ik))
                average(Datak_read, int(options.Krol.split(",")[0]), int(options.Krol.split(",")[1]), int(options.Krol.split(",")[1])==1, orbs if orbs!=[1] else None )
            Datak += Datak_read/Nit
            #
        #
        # remove positive G(tau)
        if (stat == "Fermion"):
            Datak[ Datak>0.0 ] = -1e-9
        #
        #store the averaged K-point
        np.savetxt("./avg/%s/%s/%s"%(Kfolder_out,Folder,Filek),np.c_[axis,Datak])
        print("Done average over ik=%s"%ik)
        #


















exit()
