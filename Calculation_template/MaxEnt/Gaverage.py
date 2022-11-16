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
parser.add_option("--mode"  , dest="mode"   , default="None" )                      # [G,W,C]kw, loc, models(create a model for maxent starting from an existing Akw)
parser.add_option("--El"    , dest="El"     , default="Mo" )
parser.add_option("--orbs"  , dest="orbs"   , default="1"  )
parser.add_option("--itlist", dest="itlist" , default="0"  )                        # list of iteration separated by comma
parser.add_option("--itrng" , dest="itrng"  , default="0"  )                        # initial and final iteration separated by ..
parser.add_option("--skip"  , dest="skip"   , default="0"  )                        # list of iteration to be ignored separated by comma
parser.add_option("--mag"   , dest="mag"    , default=False, action="store_true")   # consider also the spin
parser.add_option("--rol"   , dest="rol"    , default="0" )                         # string of 2 parameters separated by comma in dicating tau_max and tau_window for rolling average
parser.add_option("--log"   , dest="log"    , default=False, action="store_true")   # this is to take the rolling average in log scale
parser.add_option("--dlt"   , dest="dlt"    , default=False, action="store_true")   # this is to delete half of the points
parser.add_option("--Kpad"  , dest="Kpad"   , default=""  )                         # initial and final iteration separated by ..
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
elif options.itrng != "0":
    itrng = range( int(options.itrng.split("..")[0]) , int(options.itrng.split("..")[1])+1 )
    dirs = [ i for i in alldirs if i in map(str,itrng) ]
else:
    dirs = alldirs
#
if options.skip != "0":
    skip_dirs = list(options.skip.split(","))
    print("skip_dirs: %s"%skip_dirs)
    for dir in skip_dirs:
        dirs.remove(dir)
#
Kfolder = "K_resolved"+options.Kpad
print("Kfolder: %s"%Kfolder)
#
if options.mode=="Akw":
    removelist=[]
    for dir in dirs:
        if not os.path.isdir("%s/%s/MaxEnt_Gk_path_s1"%(dir,Kfolder)): removelist+=[dir]
    dirs = [ i for i in dirs if i not in removelist ]
#
Nit = len(dirs)
#
print("El=%s, spins=%s, orbs=%s"%(El,spins,orbs))
print("itlist=%s, Nit=%s"%(dirs,Nit))

if Nit==0:
    print("No folder found. Exiting")
    exit()

# --------------------------------------------------------------------------- #
# Here I choose to either perform the avergave on the local observables
# OR on the K-resolved ones. They're too different in kind
# --------------------------------------------------------------------------- #
if options.mode.rstrip("kw") in [ "G", "W", "C" ]:
    #
    field = options.mode.rstrip("kw")
    if field == "G":
        funct = "Gk_t"
    else:
        funct = "%sk_w"%field
        spins=[1]
    #
    for spin in spins:
        #
        pad = "_s%s"%spin if field=="G" else ""
        folder = "%s/MaxEnt_%sk_path_s%s"%(Kfolder,field,spin) if field=="G" else "%s/MaxEnt_%sk_path"%(Kfolder,field)
        #
        # Check for Nkpt consistency between iterations
        if field == "G":
            Gk = glob.glob(dirs[0]+"/%s/%s_k*.DAT"%(folder,funct))
            Nkpt = len(Gk)
            print(Nkpt)
            for dir in dirs:
                Gk = glob.glob("%s/%s/%s_k*.DAT"%(dir,folder,funct))
                print("Check for Nkpt in folder: %s ,Nkpt=%s "%(dir,len(Gk)))
                if len(Gk) != Nkpt: sys.exit("Wrong number of K-points")
        else:
            # this is temporarty as long as Im keeping the usual maxent
            for orb in orbs:
                Gk = glob.glob(dirs[0]+"/%s/%s_k*_o%s.DAT"%(folder,funct,orb))
                Nkpt = len(Gk)
                print(Nkpt)
                for dir in dirs:
                    Gk = glob.glob("%s/%s/%s_k*_o%s.DAT"%(dir,folder,funct,orb))
                    print("Check for Nkpt in folder: %s, orb=%s, Nkpt=%s "%(dir,orb,len(Gk)))
                    if len(Gk) != Nkpt: sys.exit("Wrong number of K-points")
        #
        path_exists("./avg/%s"%folder)
        #
        for ik in range(1,Nkpt+1):
            #
            if field == "G":
                #
                dir = dirs[0]
                Gpath = "%s/%s/%s_k%s.DAT"%(dir,folder,funct,ik)
                Gtauk = np.loadtxt(Gpath)
                #
                if options.rol !="0":
                    print("Rolling average dir=%s, ik=%s"%(dir,ik))
                    average(Gtauk,int(options.rol.split(",")[0]),int(options.rol.split(",")[1]),options.log)
                #
                Gtau = Gtauk.copy()/Nit
                for dir in dirs[1:]:
                    Gpath = "%s/%s/%s_k%s.DAT"%(dir,folder,funct,ik)
                    Gtauk = np.loadtxt(Gpath)
                    #
                    if options.rol !="0":
                        print("Rolling average dir=%s, ik=%s"%(dir,ik))
                        average(Gtauk,int(options.rol.split(",")[0]),int(options.rol.split(",")[1]),options.log)
                    #
                    Gtau += Gtauk/Nit
                #
                if options.dlt:
                    print("Deleting tau points, ik=%s"%ik)
                    Gtau = np.delete(Gtau, slice(None, None, 2), 0)
                #
                np.savetxt('./avg/%s/%s_k%s.DAT'%(folder,funct,ik), Gtau )
                print("Akw average ik=%s, Nkpt=%s"%(ik,Nkpt))
                #
            else:
                # this is temporarty as long as Im keeping the usual maxent
                for orb in orbs:
                    #
                    dir = dirs[0]
                    Gpath = "%s/%s/%s_k%s_o%.DAT"%(dir,folder,funct,ik,orb)
                    Gtauk = np.loadtxt(Gpath)
                    #
                    if options.rol !="0":
                        print("Rolling average dir=%s, orb=%s, ik=%s"%(dir,orb,ik))
                        average(Gtauk,int(options.rol.split(",")[0]),int(options.rol.split(",")[1]),options.log)
                    #
                    Gtau = Gtauk.copy()/Nit
                    for dir in dirs[1:]:
                        Gpath = "%s/%s/%s_k%s_o%.DAT"%(dir,folder,funct,ik,orb)
                        Gtauk = np.loadtxt(Gpath)
                        #
                        if options.rol !="0":
                            print("Rolling average dir=%s, orb=%s, ik=%s"%(dir,orb,ik))
                            average(Gtauk,int(options.rol.split(",")[0]),int(options.rol.split(",")[1]),options.log)
                        #
                        Gtau += Gtauk/Nit
                    #
                    if options.dlt:
                        print("Deleting tau points, orb=%s, ik=%s"%(orb,ik))
                        Gtau = np.delete(Gtau, slice(None, None, 2), 0)
                    #
                    np.savetxt('./avg/%s/%s_k%s_o%.DAT'%(folder,funct,ik,orb), Gtau )
                    print("Akw average ik=%s, orb=%s, Nkpt=%s"%(ik,orb,Nkpt))
                    #
    #
elif options.mode=="loc":
    #
    # Default meshes
    taulat = np.genfromtxt('%s/Convergence/Glat/Glat_t_o%s_s%s.DAT'%(dirs[0],orbs[0],spins[0]), dtype='double', usecols=(0), unpack=True, comments='#')
    tauimp = np.genfromtxt('%s/Convergence/Gqmc_%s/Gqmc_%s_t_o%s_s%s.DAT'%(dirs[0],El,El,orbs[0],spins[0]), dtype='double', usecols=(0), unpack=True, comments='#')
    wmB = np.genfromtxt('%s/Convergence/Wlat/Wlat_w_(%s,%s)(%s,%s).DAT'%(dirs[0],orbs[0],orbs[0],orbs[0],orbs[0]), dtype='double', usecols=(0), unpack=True, comments='#')
    wmF = np.genfromtxt('%s/Convergence/Glat/Glat_w_o%s_s%s.DAT'%(dirs[0],orbs[0],spins[0]), dtype='double', usecols=(0), unpack=True, comments='#')
    #
    # Dictionaries
    Glat = cllct.defaultdict(dict)
    Wlat = cllct.defaultdict(dict)
    Gqmc = cllct.defaultdict(dict)
    Sful_Gamma = cllct.defaultdict(dict)
    curlyUimp = cllct.defaultdict(dict)
    Zmat = cllct.defaultdict(dict)
    #
    # Data initialization
    idir = 1
    dir = dirs[0]
    print(dir,idir,Nit)
    for iorb in orbs:
        for ispin in spins:
            Glat[iorb][ispin] = np.genfromtxt('%s/Convergence/Glat/Glat_t_o%s_s%s.DAT'%(dir,iorb,ispin), dtype='double', usecols=(1), unpack=True, comments='#')/Nit
            Gqmc[iorb][ispin] = np.genfromtxt('%s/Convergence/Gqmc_%s/Gqmc_%s_t_o%s_s%s.DAT'%(dir,El,El,iorb,ispin), dtype='double', usecols=(1), unpack=True, comments='#')/Nit
            ReS, ImS = np.genfromtxt('%s/Convergence/Sful_Gamma/Sful_Gamma_w_o%s_s%s.DAT'%(dir,iorb,ispin), dtype='double', usecols=(1,2), unpack=True, comments='#')
            Sful_Gamma[iorb][ispin] = (ReS+1j*ImS)/Nit
        Wlat[iorb] = np.genfromtxt('%s/Convergence/Wlat/Wlat_w_(%s,%s)(%s,%s).DAT'%(dir,iorb,iorb,iorb,iorb), dtype='double', usecols=(1), unpack=True, comments='#')/Nit
        curlyUimp[iorb] = np.genfromtxt('%s/Convergence/curlyUimp/curlyUimp_w_(%s,%s)(%s,%s).DAT'%(dir,iorb,iorb,iorb,iorb), dtype='double', usecols=(1), unpack=True, comments='#')/Nit
    Umat = np.loadtxt('%s/Solver_%s/Umat.DAT'%(dir,El))/Nit
    Eloc = np.loadtxt('%s/Solver_%s/Eloc.DAT'%(dir,El))/Nit
    for st in ['dmft','qpsc']: Zmat[st] = np.loadtxt('%s/Z_%s_s1.DAT'%(dir,st))/Nit
    idir+=1
    #
    # Adding up to average
    for dir in dirs[1:]:
        print("local average ",dir,idir,Nit)
        for iorb in orbs:
            for ispin in spins:
                Glat[iorb][ispin] += np.genfromtxt('%s/Convergence/Glat/Glat_t_o%s_s%s.DAT'%(dir,iorb,ispin), dtype='double', usecols=(1), unpack=True, comments='#')/Nit
                Gqmc[iorb][ispin] += np.genfromtxt('%s/Convergence/Gqmc_%s/Gqmc_%s_t_o%s_s%s.DAT'%(dir,El,El,iorb,ispin), dtype='double', usecols=(1), unpack=True, comments='#')/Nit
                ReS, ImS = np.genfromtxt('%s/Convergence/Sful_Gamma/Sful_Gamma_w_o%s_s%s.DAT'%(dir,iorb,ispin), dtype='double', usecols=(1,2), unpack=True, comments='#')
                Sful_Gamma[iorb][ispin] += (ReS+1j*ImS)/Nit
            Wlat[iorb] += np.genfromtxt('%s/Convergence/Wlat/Wlat_w_(%s,%s)(%s,%s).DAT'%(dir,iorb,iorb,iorb,iorb), dtype='double', usecols=(1), unpack=True, comments='#')/Nit
            curlyUimp[iorb] += np.genfromtxt('%s/Convergence/curlyUimp/curlyUimp_w_(%s,%s)(%s,%s).DAT'%(dir,iorb,iorb,iorb,iorb), dtype='double', usecols=(1), unpack=True, comments='#')/Nit
        Umat += np.loadtxt('%s/Solver_%s/Umat.DAT'%(dir,El))/Nit
        Eloc += np.loadtxt('%s/Solver_%s/Eloc.DAT'%(dir,El))/Nit
        for st in ['dmft','qpsc']: Zmat[st] += np.loadtxt('%s/Z_%s_s1.DAT'%(dir,st))/Nit
        idir+=1
    #
    # Print output
    path_exists("./avg/Convergence/Glat")
    path_exists("./avg/Convergence/Sful_Gamma")
    path_exists("./avg/Convergence/Gqmc_%s"%El)
    path_exists("./avg/Convergence/Wlat")
    path_exists("./avg/Convergence/curlyUimp")
    path_exists("./avg/Solver_%s"%El)
    for iorb in orbs:
        for ispin in spins:
            np.savetxt('./avg/Convergence/Glat/Glat_t_o%s_s%s.DAT'%(iorb,ispin), np.c_[ taulat, Glat[iorb][ispin] ], delimiter='\t')
            np.savetxt('./avg/Convergence/Gqmc_%s/Gqmc_%s_t_o%s_s%s.DAT'%(El,El,iorb,ispin), np.c_[ tauimp, Gqmc[iorb][ispin] ], delimiter='\t')
            np.savetxt('./avg/Convergence/Sful_Gamma/Sful_Gamma_w_o%s_s%s.DAT'%(iorb,ispin), np.c_[ wmF, Sful_Gamma[iorb][ispin].real, Sful_Gamma[iorb][ispin].imag ], delimiter='\t')
        np.savetxt('./avg/Convergence/Wlat/Wlat_w_(%s,%s)(%s,%s).DAT'%(iorb,iorb,iorb,iorb), np.c_[ wmB, Wlat[iorb] ], delimiter='\t')
        np.savetxt('./avg/Convergence/curlyUimp/curlyUimp_w_(%s,%s)(%s,%s).DAT'%(iorb,iorb,iorb,iorb), np.c_[ wmB, curlyUimp[iorb] ], delimiter='\t')
    np.savetxt('./avg/Solver_%s/Umat.DAT'%El, Umat)
    np.savetxt('./avg/Solver_%s/Eloc.DAT'%El, Eloc)
    for st in ['dmft','qpsc']:np.savetxt('./avg/Z_%s_s1.DAT'%st, Zmat[st])
    #
elif options.mode=="models":
    #
    # Check for Nkpt consistency between iterations
    Gk = glob.glob(dirs[0]+"/%s/MaxEnt_Gk_path_s1/Gk_t_k*.DAT_dos.dat"%Kfolder)
    Nkpt = len(Gk)
    for dir in dirs:
        Gk = glob.glob("%s/%s/MaxEnt_Gk_path_s1/Gk_t_k*.DATT_dos.dat"%(dir,Kfolder))
        print("Check for Nkpt in folder: %s ,Nkpt=%s "%(dir,len(Gk)))
        if len(Gk) != Nkpt: sys.exit("Wrong number of K-points")
    #
    for ik in range(1,Nkpt+1):
        for dir in dirs:
            #
            Gpath = "%s/%s/MaxEnt_Gk_path_s1/"%(dir,Kfolder)
            print("dir: %s, ik: %s"%(dir,ik))
            #
            Gtauk = np.loadtxt( Gpath+"Gk_t_k%s.DAT"%ik )
            if options.rol !="0": average(Gtauk,int(options.rol.split(",")[0]),int(options.rol.split(",")[1]),options.log)
            Grealk = np.loadtxt( Gpath+"Gk_t_k%s.DAT_dos.dat"%ik )
            if(np.shape(Gtauk)[1]==np.shape(Grealk)[1]):
                Norb=np.shape(Gtauk)[1]
            else:
                print("Norb_Gt=",np.shape(Gtauk)[1],"Norb_Gr=",np.shape(Grealk)[1])
            #
            if(ik==1): path_exists("./%s/%s_models/MaxEnt_Gk_path_s1/Models"%(dir,Kfolder))
            #
            # I'm using the trace as model
            #np.savetxt("%s/%s_models/MaxEnt_Gk_path_s1/Gk_t_k%s_o%s.DAT"%(dir,Kfolder,ik,io), np.c_[Gtauk[:,0],Gtauk[:,io]] )
            #np.savetxt("%s/%s_models/MaxEnt_Gk_path_s1/Models/Gk_t_k%s_o%s.DAT_dos.dat"%(dir,Kfolder,ik,io), np.c_[Grealk[:,0],np.sum(Grealk,axis=1)/Norb )
            #
            #for io in range(1,Norb):
            #    np.savetxt("%s/K_resolved_models/MaxEnt_Gk_path_s1/Gk_t_k%s_o%s.DAT"%(dir,ik,io), np.c_[Gtauk[:,0],Gtauk[:,io]] )
            #    np.savetxt("%s/K_resolved_models/MaxEnt_Gk_path_s1/Models/Gk_t_k%s_o%s.DAT_dos.dat"%(dir,ik,io), np.c_[Grealk[:,0],Grealk[:,io]] )










exit()
