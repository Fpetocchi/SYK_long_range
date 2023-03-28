#!/usr/bin/env python
import glob
import numpy as np
import collections as cllct
import math
from optparse import OptionParser

# ----------------------------------------------------------------------------#

parser = OptionParser()
#Mandatory input
parser.add_option("--Gt_folder", dest="Gt_folder", default="./MaxEnt_Gk_path_t_s1" )
parser.add_option("--Ak_folder", dest="Ak_folder", default="./MaxEnt_Gk_path_t_s1" )
#optional
(options, args) = parser.parse_args()

#
#Check number of k files within Gt_folder
NkptG=len( glob.glob("%s/Gk_t_k*.DAT"%options.Gt_folder) )
print("Nkp(G) = %s"%NkptG)

#
#Check number of k files within Ak_folder
NkptA=len( glob.glob("%s/*.DAT"%options.Ak_folder) )
print("Nkp(A) = %s"%NkptA)

#
assert(NkptG==NkptA)
Nkpt=NkptG

Aprefix=glob.glob("%s/*.DAT"%options.Ak_folder)[0].split("/")[1].split("_k")[0]
Aprefix+="_k"
print("Akw prefix: %s"%Aprefix)

# Read in reference G(tau)
Gtau = cllct.defaultdict(dict)
err=np.zeros(Nkpt)
for ik in range(1,Nkpt+1):
    #
    Gtau=np.loadtxt("%s/Gk_t_k%s.DAT"%(options.Gt_folder,ik))
    Akw_read=np.loadtxt("%s/%s%s.DAT"%(options.Ak_folder,Aprefix,ik))
    #
    Norb=len(Gtau[0,:])-1
    tau=Gtau[:,0]  ;Ntau=len(tau); Beta=Gtau[Ntau-1,0]
    #
    thresh=1230.0/Beta
    #
    wr=Akw_read[abs(Akw_read[:,0])<thresh,0]    ;Nwr=len(wr)  ; dw=abs(wr[10]-wr[9])
    Akw=Akw_read[abs(Akw_read[:,0])<thresh,:]
    #
    if ik==1:
        K=dw*np.array([[np.exp((Beta/2-t)*w)/(np.exp(Beta/2*w)+np.exp(-Beta/2*w)) for w in wr] for t in tau], dtype=np.float64)
        print("Dimension of kernel:",K.shape)
        print("Max real frequency:",thresh)
        print(wr)
        print(tau)
    #
    print("data at ik=%s are read. Norb=%s Ntau=%s Nwr=%s dw=%s Beta=%s"%(ik,Norb,Ntau,Nwr,dw,Beta))
    #
    Grebuilt=np.zeros_like(Gtau)
    Grebuilt[:,0]=tau
    #
    for iorb in range(Norb):
        for itau in range(Ntau):
            Grebuilt[itau,1+iorb] = -np.dot(K[itau,:],Akw[:,1+iorb])
            err[ik-1]+=abs(Grebuilt[itau,1+iorb]-Gtau[itau,1+iorb])/(Norb*Ntau)
    print("error= %s"%err[ik-1])
    #
    np.savetxt("%s/%s%s.DAT_rebuilt"%(options.Ak_folder,Aprefix,ik),Grebuilt)
np.savetxt("%s/Errors.DAT"%options.Ak_folder,np.c_[range(1,Nkpt+1),err])
