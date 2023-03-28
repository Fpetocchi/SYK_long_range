#!/usr/bin/env python
import numpy as np
import collections as cllct
import math
from optparse import OptionParser

# ----------------------------------------------------------------------------#

def Kramers_Kronig_imag2real(ImPart,wr,wlimit=1000,screenval=None,bareval=None):
    #
    # I should take tha bare but it depends on HighestMatsuFreq
    #
    dw=abs(wr[10]-wr[9])
    Lreal=len(wr)
    RePart=np.zeros(Lreal)
    #
    for iw in np.argwhere(wr<wlimit):
        print(wr[iw])
        freqdep=0.0
        for jw in range(Lreal):
            if iw!=jw: freqdep += 2.0/(math.pi) *dw * ( wr[jw] * ImPart[jw] ) / ( wr[jw]**2 - wr[iw]**2  )
        RePart[iw] = freqdep if bareval is not None else (bareval + freqdep)
    if screenval is not None : RePart -= (RePart[0] - screenval)
    #
    return RePart

# ----------------------------------------------------------------------------#


parser = OptionParser()
#Mandatory input
parser.add_option("--Bands", dest="Bands", default="Bands.DAT"   )
parser.add_option("--Norb" , dest="Norb" , default="3"      )
parser.add_option("--Emax" , dest="Emax" , default="5"      )
parser.add_option("--mode" , dest="mode" , default="G"      )
#optional
parser.add_option("--s"    , dest="spin" , default="1"      )
parser.add_option("--eta"  , dest="eta"  , default="0.05"   )
(options, args) = parser.parse_args()

Norb  = int(options.Norb)
ispin = int(options.spin)
Ethrs = float(options.Emax)
eta   = float(options.eta)
MaxEntFolder = "./%skt_s%s"%(options.mode,options.spin)


#readin the Bands
BandsRead = np.loadtxt( options.Bands )
NKtot = len( BandsRead[:,0] )-1 #usually Im not doing MaxEnt on the last Kpoint
kaxis = BandsRead[:,1]
Bands = np.array(BandsRead[0,2:Norb+2]).copy()


if options.mode=="G":
    #
    #Read the spectral functions directly
    Akw = cllct.defaultdict(dict)
    for ik in range(NKtot):
        path = MaxEntFolder+"/Gkt_k%s.DAT_dos.dat"%(ik+1)
        print(path)
        try:
            Akw[str(ik)] = np.loadtxt( path )
            if(len(Akw[str(ik)][0,:])-1 != Norb):
                print("Wrong number of orbitals in ik: %s"%ik)
                print("len(Akw[str(ik)][0,:])-1: %s"%(len(Akw[str(ik)][0,:])-1))
                print("len( Bands[0,:] )-2: %s"%Norb)
                exit()
            Akw["mask"][str(ik)] = True
        except Exception as e:
            Akw["mask"][str(ik)] = False
            print(e)
            pass
    #
    #print some info
    print("NKtot: %s"%NKtot     )
    print("Kmin : %s"%kaxis[0]  )
    print("Kmax : %s"%kaxis[-1] )
    print("Nktot: %s"%NKtot     )
    print("Norb : %s"%Norb      )
    #
    #join the spectral functions
    Akw_f = open("Akw_s%s.DAT"%ispin,"w")
    for ik in range(NKtot):
        #
        if Akw["mask"][str(ik)]:
            #
            wr=Akw[str(ik)][:,0]
            #
            for iw in range(len(wr)):
                #
                if(abs(wr[iw])<Ethrs):
                    #
                    Akw_f.write("%4d %10.5f %12.8f "%(ik+1,kaxis[ik],wr[iw]) )
                    for iorb in range(Norb): Akw_f.write("%30.20e "%(Akw[str(ik)][iw,iorb+1]/np.sum(Akw[str(ik)][:,iorb+1])))
                    Akw_f.write("\n" )
                    #
                #
            Akw_f.write("\n" )

    Akw_f.close()
    #
elif options.mode=="S":
    #
    #Read the self-energy on the real axis
    Skw = cllct.defaultdict(dict)
    for ik in range(NKtot):
        path = MaxEntFolder+"/Skt_k%s.DAT_dos.dat"%(ik+1)
        print(path)
        try:
            Skw["Im"][str(ik)] = np.loadtxt( path )
            Skw["Re"][str(ik)] = np.zeros_like(Skw["Im"][str(ik)])
            if(len(Skw["Im"][str(ik)][0,:])-1 != Norb):
                print("Wrong number of orbitals in ik: %s"%ik)
                print("len(Skw[str(ik)][0,:])-1: %s"%(len(Skw["Im"][str(ik)][0,:])-1))
                print("len( Bands[0,:] )-2: %s"%Norb)
                exit()
            Skw["mask"][str(ik)] = True
        except Exception as e:
            Skw["mask"][str(ik)] = False
            print(e)
            pass
    #
    #print some info
    print("NKtot: %s"%NKtot     )
    print("Kmin : %s"%kaxis[0]  )
    print("Kmax : %s"%kaxis[-1] )
    print("Nktot: %s"%NKtot     )
    print("Norb : %s"%Norb      )
    #
    #Print it
    Skw_f = open("Skw_s%s.DAT"%ispin,"w")
    for ik in range(NKtot):
        #
        if Skw["mask"][str(ik)]:
            #
            wr=Skw["Im"][str(ik)][:,0]
            #
            for iw in range(len(wr)):
                #
                if(abs(wr[iw])<Ethrs):
                    #
                    Skw_f.write("%4d %10.5f %12.8f "%(ik+1,kaxis[ik],wr[iw]) )
                    for iorb in range(Norb): Skw_f.write("%30.20e "%(Skw["Im"][str(ik)][iw,iorb+1] ) )
                    Skw_f.write("\n" )
                    #
                #
            Skw_f.write("\n" )
    Skw_f.close()
    #
    #
    #--------------------------------------------------------------------------#
    # THE FOLLOWING IS TOO HEAVY TO BE PERFORMED WITH SERIAL PYTHON
    # A PROPER EXECUTABLE IS CONTAINED IN THE REPOSITORY
    #--------------------------------------------------------------------------#
    exit()
    #
    #Read the correction coefficient and the bare real limit
    file_Z=open("Spath_vars/Spath_Rot_k%s.DAT"%ik,'r')
    Params = cllct.defaultdict(dict)
    Params["Bare"] = np.zeros((NKtot+1,Norb))
    Params["Coef"] = np.zeros((NKtot+1,Norb))
    shift=2*Norb*(ispin-1)
    for iorb in range(Norb):
        Params["Bare"][:,iorb] = np.genfromtxt("Spath_vars/Spath_Params.DAT", dtype='double', comments='#', usecols=(iorb+shift), unpack=True)
        Params["Coef"][:,iorb] = np.genfromtxt("Spath_vars/Spath_Params.DAT", dtype='double', comments='#', usecols=(Norb+iorb+shift), unpack=True)
    #
    #Manipulations
    for ik in range(NKtot):
        #
        if Skw["mask"][str(ik)]:
            #
            for iorb in range(Norb):
                #
                #Correct the magnitude
                Skw["Im"][str(ik)][:,iorb+1] *= Params["Coef"][ik,iorb]
                #
                #Flip the poles
                Skw["Im"][str(ik)][:,iorb+1] = Skw["Im"][str(ik)][::-1,iorb+1]
                #
                #Build the real part
            #    Skw["Re"][str(ik)][:,iorb+1] = Kramers_Kronig_imag2real(Skw["Im"][str(ik)][:,iorb+1],Skw["Im"][str(ik)][:,0],wlimit=Ethrs,screenval=None,bareval=Params["Bare"][ik,iorb])
            #    print("Done KK on ik=%s and iorb=%s"%(ik,iorb))
