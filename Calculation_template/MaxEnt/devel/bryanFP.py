#!/usr/bin/env python
# Written by Francesco Petocchi (francesco.petocchi@gmail.com), 2022-2023
import numpy as np
import sys
import os
from optparse import OptionParser
from scipy.interpolate import interp1d
#?import warnings
#?import itertools
#?import time

parser = OptionParser()
parser.add_option("--data"      , dest="data"      , default="None" , help="Path to datafile" )
parser.add_option("--kernel"    , dest="kernel"    , default="tau"  , help="Imaginary time (tau-deafult) or imaginary frequency (mats) kernel." )
parser.add_option("--S"         , dest="statistics", default="F"    , help="Fermionic (F-deafult) Bosonic (B) continuation." )
parser.add_option("--W"         , dest="Wrange"    , default="None" , help="Frequncy range. Symmetric for -S F, positive for -S X." )
parser.add_option("--Nw"        , dest="Nreal"     , default="None" , help="Number of real frequency points." )
parser.add_option("--sigma"     , dest="sigma"     , default="None" , help="Standard deviation of the input data. Orbital independent." )
parser.add_option("--sigmafile" , dest="sigmafile" , default="None" , help="Path to standard deviation file. Orbital dependent." )
parser.add_option("--moments"   , dest="moments"   , default="None" , help="M1,M2 of the default model. Comma separated. Orbital independent." )
parser.add_option("--modelfile" , dest="modelfile" , default="None" , help="Path to model file. Orbital dependent." )
parser.add_option("--verb"      , dest="verb"      , default=False  , help="Verbosity", action="store_true")
(options, args) = parser.parse_args()

# Few preliminary checks
if options.data == "None":
    print >> sys.stderr, "--data parameter not provided."
    sys.exit(1)
if options.Wrange == "None":
    print >> sys.stderr, "--W parameter not provided."
    sys.exit(1)
if options.Nreal == "None":
    print >> sys.stderr, "--Nw parameter not provided."
    sys.exit(1)
if options.sigma == "None" and options.sigmafile == "None":
    print >> sys.stderr, "neither --sigma nor --sigmafile parameter provided."
    sys.exit(1)
if options.statistics != "F" and options.statistics != "B":
    print >> sys.stderr, "--S wrong entry. Only F and B allowed."
    sys.exit(1)
if options.kernel != "tau" and options.kernel != "mats":
    print >> sys.stderr, "--kernel wrong entry. Only tau and mats allowed."
    sys.exit(1)
if options.sigma != "None" and options.sigmafile != "None":
    print >> "--sigma and --sigmafile provided. --sigma overrides."
    options.sigmafile = "None"
if options.moments != "None" and options.modelfile != "None":
    print >> "--moments and --modelfile provided. --moments overrides."
    options.modelfile == "None":

# Link to variables
Wrange = float(options.Wrange)
Nreal = int(options.Nreal)

# global variable
precD = np.float128
precC = np.complex128
M0 = 1.0 #mettimi nella definizione dei kernel
verb = options.verb
import bryan_lib
#
#
################################################################################
#                             DATA PRE-PROCESSING                              #
################################################################################
#
#
# Setting up real frequency mesh
wreal = np.linspace(-Wrange/2.0,Wrange/2.0,num=Nreal)
#
#
# Setting up input field
data = np.loadtxt( options.data, dtype=precD )
if options.kernel == "tau":
    inpMesh = data[:,0] #(tau)
    Beta = inpMesh[0]+inpMesh[-1]
    # Both statistics are real in tau
    Field = data[:,1:]
elif options.kernel == "mats":
    inpMesh = data[:,0] #(wmats)
    Beta = np.pi/inpMesh[0]
    # Only fermion are assumed to be complex
    if options.statistics == "F":
        Field = data[:,1::2] + 1j*data[:,2::2]
    elif options.statistics == "B":
        Field = data[:,1:]
#
#
# Input dimensions
Npts = np.shape(Field)[0] # input Field points
Comp = np.shape(Field)[1] # input Field components
#
#
# Setting up statistical error
if options.sigma != "None":
    # Component independent error
    sigmapm2 = np.ones_like( Field, dtype=precD ) / float(options.sigma)**2
elif options.sigmafile != "None":
    # Component dependent error only from file
    sigmapm2 = 1.0 / np.power( np.loadtxt( options.sigmafile, dtype=precD ), 2 )
    # check shape
    if np.shape(sigmapm2) != np.shape(Field):
        print >> sys.stderr, "wrong shape in "+options.sigmafile
        sys.exit(1)
#
#
# Setting up deafult model
model = np.ones( (Nreal,Comp), dtype=precD )
# Wait to code the interpolation before coding the model for bosons
if options.statistics == "F":
    # User-provided default model. Component dependent.
    if options.modelfile != "None":
        # Read the component-dependent model from file
        wreal_model_read = np.loadtxt( options.modelfile )[:,0]
        model_read = np.abs(np.loadtxt( options.modelfile )[:,1:])
        # check shape
        if np.shape(model_read)[1] != Comp:
            print >> sys.stderr, "wrong shape in "+options.modelfile
            sys.exit(1)
        # interpolate to user-provided real frequency mesh
        model *= 0.0
        for c in range(Comp):
            f = interp1d( wreal_model_read, model_read[:,c], kind='linear', fill_value=1e-15, bounds_error=False )
            model[:,c] = f(wreal)
        #
    # Gaussian default model. Component independent
    if options.moments != "None":
        # Get the moments
        Mvec = np.array( options.moments.split(","), dtype=precD )
        M1 = Mvec[0]
        M2 = Mvec[1]
        # Create a gaussian model for all the components
        model *= 0.0
        for c in range(Comp):
            model[:,c] = M0 * np.exp( -(wreal-M1)**2 / (2.0*M2**2) ) / np.sqrt(2*np.pi*M2**2)+1e-8
        #
    # Non-normalized moments (only for Experts)
    if M0 != 1.0:
        for c in range(Comp):
            model[:,c] *= M0 / np.trapz( model[:,c], wreal )
    #
#
#
#
################################################################################
#                               MAXIMUM ENTROPY                                #
################################################################################
#
#
# Setting up kernels
B = Bryan( inpMesh, wreal, Beta, options.statistics, options.kernel )
#
# MaxEnt on each component of the input Field
for c in range(Comp):
    # Initialize M matrix
    B.initField( Field[:,c], sigmapm2[:,c], model[:,c] )
    #
    ( aMax, AMax, S, L, lam ) = B.approxMax()
    #
    X = ( 0.5*np.sum(np.log(aMax/(aMax+lam))) + aMax*S - L ) - 2
    F = MeasureFunctor( inpMesh, wreal, X, aMax, AMax, S, L, lam, Field[:,c] )
    #
    B.screen( np.log(aMax), +0.1, AMax, F.add )
    B.screen( np.log(aMax), -0.1, AMax, F.add )
    #
    ABryan = F.finalize()



class Bryan:

    def approxMax(self,aMax=1.0):
        AMax=self.model.copy()
        while True:
            (AMax,S,L,lam,_)=self.maxent(aMax,AMax)
            tmp=np.sum(lam/(aMax+lam))/(-2*aMax*S)
            if abs(tmp-1)<1e-4: break
            print(tmp)
            aMax*=(tmp-1)*1.0+1
        return(aMax,AMax,S,L,lam)


    def maxent(self,alpha,initA=None):


        print("Calculating maxent with alpha=%s"%alpha)
        A=np.copy(initA)
        A+=1e-8
        u=np.dot(self.TU.transpose(),np.log(A*self.oneovermodel))
        #
        mu=1e-2
        #muthres=1.0
        #muthres=0.001
        #
        iter=0
        while True:
            #
            G_=np.dot(self.TVt.transpose(),np.dot(np.diag(self.TSigma),np.dot(self.TU.transpose(),A))) # G_=K.dot(A)
            dG=G_+self.G
            g=np.dot(np.diag(self.TSigma),np.dot(self.TVt,dG*self.sigmapm2)) # (10) and page 166 left column upper half
            K=np.dot(np.dot(self.TU.transpose(),np.diag(A)),self.TU) # before (B11)
            MK=np.dot(self.M,K)
            du=np.linalg.solve((alpha)*np.eye(MK.shape[0])+MK, -alpha*u-g) # (B11)
            mdu=np.dot(np.dot(du,K),du) # after (B12)
            #
            if verbose: print("mdu=%s"%(mdu))
            #
            if mdu<1e-7: break # there is an alternate criterion in Bryan's paper
            if mdu>muthres:
                du=np.linalg.solve((alpha+mu)*np.eye(MK.shape[0])+MK, -alpha*u-g) # (B12)
                mdu=np.dot(np.dot(du,K),du) #after (B12)
                while (abs(mdu-muthres)>1e-3*muthres):
                    mu*=1+(mdu-muthres)
                    du=np.linalg.solve((alpha+mu)*np.eye(MK.shape[0])+MK, -alpha*u-g) # (B12)
                    mdu=np.dot(np.dot(du,K),du) # after (B12)
            #
            u+=du
            A=self.model*np.exp(np.dot(self.TU,u)) # (9)
            iter=iter+1
            print("mdu=%s it=%s"%(mdu,iter))
            #



            
        lam=np.linalg.eigvalsh(MK)
        tmp=np.log(np.where(A>1e-6,A*(self.oneovermodel+1e-8),1e-6*self.oneovermodel))
        S=np.sum(A-model-(A*tmp)[tmp>-np.inf])
        L=0.5*np.sum((G_+self.G)**2*self.deltaTau*self.sigmapm2) # page 166 left column upper half
        print("Q=%s"%(alpha*S-L))
        return (A,S,L,lam,-G_)



    def screen (self, log_a, step, initA, functor):
        A=np.copy(initA)
        while True:
            log_a+=step
            a=np.exp(log_a)
            (A,S,L,lam,G)=self.maxent(a,A)
            if not (functor(a,A,self.model,S,L,lam,G)):
                break

class MeasureFunctor:
    def __init__(self, tauMesh, omegaMesh, thresh, a, A, S, L, lam,G):
        if gnuplot:
            self.plot2("set logscale y")
        self.G=G
        self.omegaMesh=omegaMesh
        self.tauMesh=tauMesh
        self.thresh=thresh
        logP=0.5*np.sum(np.log(a/(a+lam)))+a*S-L
        self.logAlpha=[np.log(a)]
        self.logPrAlpha=[logP]
        self.sumPrAlpha=np.exp(logP)
        self.sumPrAlphaA=np.exp(logP)*np.copy(A)
    def add (self, a, A, model, S, L, lam, G):
        logP=0.5*np.sum(np.log(a/(a+lam)))+a*S-L
        self.logAlpha.append(np.log(a))
        self.logPrAlpha.append(logP)
        self.sumPrAlpha+=np.exp(logP)
        self.sumPrAlphaA+=np.exp(logP)*A
        print("logP=%s, thresh=%s"%(logP, self.thresh))
        return (logP>self.thresh)
    def finalize (self):
        return self.sumPrAlphaA/self.sumPrAlpha

if __name__=="__main__":

    # Read doc string from README
    path = os.path.dirname(__file__) # find module path
    filename = os.path.join(path, 'README') # path to module README
    __doc__ = ''.join(open(filename, 'r').readlines()[:16]) # get __doc__ string
    from docopt import docopt
    arguments=docopt(__doc__)



    ########### READ FILE ###########
    if arguments['--statistics']=='X':
        from scipy.optimize import leastsq
        Ns=1
        Ww=np.loadtxt(arguments['<datafile>'])[:,:2]
        Ww[:,1]*=-1
        Beta=2*np.pi/Ww[1,0]
        p,c,d=leastsq(lambda p,m,d: [ np.sum( [(pi/((1j*mi)**(2*i))).real for i,pi in enumerate(p)] )-di for mi,di in zip(m,d) ], [1,-2,-4] , args=(Ww[-5:,0],Ww[-15:,1]) )[0]
        Ww[:,1]-=p
        tauMesh=np.linspace(0,Beta,int(Ww.shape[0]/2-1))
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            oneoverWwmesh=np.where(Ww[:,0]>0.0,1/Ww[:,0],0.0)
        G=(np.array([np.sum(np.cos(Ww[:,0]*t)*(Ww[:,1]+c*oneoverWwmesh**2-d*oneoverWwmesh**4)) for t in tauMesh])*2.0-Ww[0,1])/Beta+c*(-1.0/12*Beta+1.0/2*tauMesh-tauMesh**2/2.0/Beta)+d*(Beta**3-30*Beta*tauMesh**2+60*tauMesh**3-30./Beta*tauMesh**4)/720.
        G=G.reshape(G.shape[0],1)


        omegaMesh=np.linspace(0,float(arguments['--frequencyrange']),int(arguments['--frequencies']))
        if arguments['--modelfile']:
            origomegaMesh=np.loadtxt(arguments['--modelfile'])[:,0]
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
            origmodel=-np.loadtxt(arguments['--modelfile'])[0:,-1]
            f=interpolate.interp1d(origomegaMesh,origmodel,kind='linear',fill_value=1e-15,bounds_error=False)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                oneoveromegaMesh=np.where(omegaMesh==0.0,0.0,1/omegaMesh)
            model=-f(omegaMesh)*oneoveromegaMesh*2.0/np.pi/Ww[0,1]+1e-8
        #else:
        model=omegaMesh*0.0+1.0

        model*=float(arguments["--modelnorm"])/np.trapz(model,omegaMesh)
        B=Bryan(tauMesh,omegaMesh,Beta,model,kernel='realbosonic',norm=-Ww[0,1])



    elif arguments['--statistics']=='F':

        ######### MAXENTING #############

        B=Bryan(tauMesh, omegaMesh, Beta, model)
    else:
        print >> sys.stderr, "only real bosonic (X) and fermionic (F) continuation for now"
        sys.exit(1)








    outputA=[]
    print("shape(G)=%s, shape(sigmapm2=%s)"%(G.shape,sigmapm2.shape))
    for i in range(G.shape[1]):

        B.initG(G[:,i], sigmapm2[:,i])

        (aMax,AMax,S,L,lam)=B.approxMax()

        F=MeasureFunctor(tauMesh,omegaMesh,(0.5*np.sum(np.log(aMax/(aMax+lam)))+aMax*S-L)-2, aMax, AMax, S, L, lam,G[:,i])


        B.screen(np.log(aMax), 0.1, AMax, F.add)
        B.screen(np.log(aMax), -0.1, AMax, F.add)
        ABryan=F.finalize()

        if arguments['--statistics']=='F':
            outputA.append(np.copy(ABryan))
        elif arguments['--statistics']=='X':
            outputA.append(ABryan*B.omegaMesh*np.pi/2.0*Ww[0,1])
        elif arguments['--statistics']=='Xt':
            outputA.append(ABryan*B.omegaMesh*np.pi/2.0*np.trapz(G[:,0],tauMesh)/Beta)

    ####### WRTINTING FILE #########

    print("writing output")
    outfile=arguments['--outfile'] if arguments['--outfile'] else arguments['<datafile>']+'_dos.dat'
    np.savetxt(outfile,np.hstack([B.omegaMesh.reshape(B.omegaMesh.shape[0],1),np.array(outputA).transpose()]))
    f=open(outfile,'a')
    f.write("#"+" ".join(sys.argv)+"\n")
    f.close()

    GtauRepFile=open(outfile+"_reproduced.dat", 'w')
    for i in range(len(tauMesh)):
        s="%s  "%tauMesh[i]
        ## for j in range (Ns):
        s+="%s  "%(np.dot(B.K[i,:],ABryan))
        s+="\n"
        GtauRepFile.write(s)
    GtauRepFile.close()
    print("done")
