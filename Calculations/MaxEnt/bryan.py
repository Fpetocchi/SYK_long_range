#!/usr/bin/env python
# Written by Lewin Boehnke (lboehnke@physnet.uni-hamburg.de), 2009-2015, consulted by Hugo Strand, 2015

try:
    import Gnuplot
    from Gnuplot import Data
    gnuplot=True
except ImportError:
    gnuplot=False

import numpy as np
import warnings
import itertools
import os
import sys
import time
#from scipy import interpolate

class Bryan:
    def __init__(self, tauMesh, omegaMesh, Beta, model, kernel='fermionic',norm=1.0):
        self.tauMesh=tauMesh
        self.omegaMesh=omegaMesh
        self.Beta=Beta*1.0
        self.deltaTau=tauMesh[1]-tauMesh[0]
        self.deltaOmega=omegaMesh[1]-omegaMesh[0]
        if kernel=='fermionic':
            self.K=self.deltaOmega*np.array([[np.exp((Beta/2-tau)*omega)/(np.exp(Beta/2*omega)+np.exp(-Beta/2*omega)) for omega in omegaMesh] for tau in tauMesh], np.float64)
        elif kernel=='realbosonic':
            self.K=self.deltaOmega*np.array([[(np.cosh((Beta/2-tau)*omega)*(omega/np.sinh(Beta/2*omega) if abs(omega)>1e-8 else 2.0/Beta)*norm/2.0) if Beta*omega<700.0 else (np.exp(-tau*omega)+np.exp(-(Beta-tau)*omega))*omega*norm/2.0 for omega in omegaMesh] for tau in tauMesh], np.float64)
        print "Calculating SVD"
        (V,Sigma,Ut)=np.linalg.svd(self.K)
        print "Done"
        self.model=np.array(model,np.float64)
        self.oneovermodel=np.where(model>0.0,1.0/model,0.0)
        good=[i for i in range (Sigma.shape[0]) if Sigma[i]>Sigma[-1]*100]
        self.s=len(good)
        self.TSigma=Sigma[good]
        self.TU=Ut.transpose()[:,good]
        self.TVt=V.transpose()[good,:]
        print "Dimension of kernel:%s; Dimension of singular space:%s"%(self.K.shape, self.s)

    def initG(self,G,sigmapm2):
        self.G=np.array(G,np.float64)
        self.sigmapm2=np.array(sigmapm2,np.float64)
        self.M=np.dot(np.dot(np.diag(self.TSigma),np.dot(np.dot(self.TVt,np.diag(self.sigmapm2)),self.TVt.transpose())),np.diag(self.TSigma)) #before (11)

    def maxent(self,alpha,initA=None):
        if initA==None:
            initA=np.array([1.0/(np.sqrt(2*np.pi)*0.5)*np.exp(-((x)/0.5)**2*0.5)+1e-3 for x in self.omegaMesh],np.float64)
        print "Calculating maxent with alpha=%s"%alpha
        A=np.copy(initA)
        A+=1e-8
        u=np.dot(self.TU.transpose(),np.log(A*self.oneovermodel))
        mu=1e-2
        while True:
            G_=np.dot(self.TVt.transpose(),np.dot(np.diag(self.TSigma),np.dot(self.TU.transpose(),A))) # G_=K.dot(A)
            dG=G_+self.G
            g=np.dot(np.diag(self.TSigma),np.dot(self.TVt,dG*self.sigmapm2)) # (10) and page 166 left column upper half
            K=np.dot(np.dot(self.TU.transpose(),np.diag(A)),self.TU) # before (B11)
            MK=np.dot(self.M,K)
            du=np.linalg.solve((alpha)*np.eye(MK.shape[0])+MK, -alpha*u-g) # (B11)
            mdu=np.dot(np.dot(du,K),du) # after (B12)
            if mdu<1e-7: break # there is an alternate criterion in Bryan's paper
            if mdu>1:
                du=np.linalg.solve((alpha+mu)*np.eye(MK.shape[0])+MK, -alpha*u-g) # (B12)
                mdu=np.dot(np.dot(du,K),du) #after (B12)
                while (abs(mdu-1)>1e-4):
                    mu*=1+(mdu-1)
                    du=np.linalg.solve((alpha+mu)*np.eye(MK.shape[0])+MK, -alpha*u-g) # (B12)
                    mdu=np.dot(np.dot(du,K),du) # after (B12)
            u+=du
            A=self.model*np.exp(np.dot(self.TU,u)) # (9)
        lam=np.linalg.eigvalsh(MK)
        tmp=np.log(np.where(A>1e-6,A*(self.oneovermodel+1e-8),1e-6*self.oneovermodel))
        S=np.sum(A-model-(A*tmp)[tmp>-np.inf])
        L=0.5*np.sum((G_+self.G)**2*self.deltaTau*self.sigmapm2) # page 166 left column upper half
        print "Q=%s"%(alpha*S-L)
        return (A,S,L,lam,-G_)

    def approxMax(self,aMax=1.0):
        AMax=self.model.copy()
        while True:
            (AMax,S,L,lam,_)=self.maxent(aMax,AMax)
            tmp=np.sum(lam/(aMax+lam))/(-2*aMax*S)
            if abs(tmp-1)<1e-4: break
            print tmp
            aMax*=(tmp-1)*1.0+1
        return(aMax,AMax,S,L,lam)

    def screen (self, log_a, step, initA, functor):
        A=np.copy(initA)
        while True:
            log_a+=step
            a=np.exp(log_a)
            (A,S,L,lam,G)=self.maxent(a,A)
            if not (functor(a,A,self.model,S,L,lam,G)):
                break

class MeasureFunctor:
    if gnuplot:
        plot=Gnuplot.Gnuplot()
        plot2=Gnuplot.Gnuplot()
        plot3=Gnuplot.Gnuplot()
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
        if gnuplot:
            self.plot.plot(Data(self.logAlpha,self.logPrAlpha, with_="l"), Data(self.logAlpha,np.exp(self.logPrAlpha), with_="l", axes="x1y2"))
            self.plot2.plot(Data(self.omegaMesh, A, with_="lp"), Data(self.omegaMesh, model, with_='l ls 3',inline=True))
            self.plot3.plot(Data(self.tauMesh,self.G,with_="p"),Data(self.tauMesh,G,with_="l ls 3",inline=True))
        print "logP=%s, thresh=%s"%(logP, self.thresh)
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

    if not(arguments['--gnuplot']):
        gnuplot=False

    ########### READ FILE ###########
    if arguments['--statistics']=='X':
        from scipy.optimize import leastsq
        Ns=1
        Ww=np.loadtxt(arguments['<datafile>'])[:,:2]
        Ww[:,1]*=-1
        Beta=2*np.pi/Ww[1,0]
        p,c,d=leastsq(lambda p,m,d: [ np.sum( [(pi/((1j*mi)**(2*i))).real for i,pi in enumerate(p)] )-di for mi,di in itertools.izip(m,d) ], [1,-2,-4] , args=(Ww[-5:,0],Ww[-15:,1]) )[0]
        Ww[:,1]-=p
        tauMesh=np.linspace(0,Beta,Ww.shape[0]/2-1)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            oneoverWwmesh=np.where(Ww[:,0]>0.0,1/Ww[:,0],0.0)
        G=(np.array([np.sum(np.cos(Ww[:,0]*t)*(Ww[:,1]+c*oneoverWwmesh**2-d*oneoverWwmesh**4)) for t in tauMesh])*2.0-Ww[0,1])/Beta+c*(-1.0/12*Beta+1.0/2*tauMesh-tauMesh**2/2.0/Beta)+d*(Beta**3-30*Beta*tauMesh**2+60*tauMesh**3-30./Beta*tauMesh**4)/720.
        G=G.reshape(G.shape[0],1)
        if arguments["--sigma"]:
            sigmapm2=G*0.0+1.0/float(arguments['--sigma'])**2
        else:
            print >> sys.stderr, "--sigma has to be specified for -S X"
            sys.exit(1)

        if arguments['--moments']:
            print >> sys.stderr, "--moments not implemented for -S X"
            sys.exit(1)
        omegaMesh=np.linspace(0,float(arguments['--frequencyrange']),int(arguments['--frequencies']))
        #if arguments['--modelfile']:
        #    origomegaMesh=np.loadtxt(arguments['--modelfile'])[:,0]
        #    with warnings.catch_warnings():
        #        warnings.simplefilter("ignore")
        #    origmodel=-np.loadtxt(arguments['--modelfile'])[0:,-1]
        #    f=interpolate.interp1d(origomegaMesh,origmodel,kind='linear',fill_value=1e-15,bounds_error=False)
        #    with warnings.catch_warnings():
        #        warnings.simplefilter("ignore")
        #        oneoveromegaMesh=np.where(omegaMesh==0.0,0.0,1/omegaMesh)
        #    model=-f(omegaMesh)*oneoveromegaMesh*2.0/np.pi/Ww[0,1]+1e-8
        #else:
        model=omegaMesh*0.0+1.0

        model*=float(arguments["--modelnorm"])/np.trapz(model,omegaMesh)
        B=Bryan(tauMesh,omegaMesh,Beta,model,kernel='realbosonic',norm=-Ww[0,1])

    elif arguments['--statistics']=='Xt':
        from scipy.optimize import leastsq
        Ns=1
        tauMesh=np.loadtxt(arguments['<datafile>'])[:,0]
        G=-np.loadtxt(arguments['<datafile>'])[:,1].reshape(tauMesh.shape[0],1)
        Beta=tauMesh[0]+tauMesh[-1]

        if arguments["--sigma"]:
            sigmapm2=G*0.0+1.0/float(arguments['--sigma'])**2
        else:
            sigmapm2=(1/np.loadtxt(arguments['<datafile>'])[:,2]**2).reshape(tauMesh.shape[0],1)

        if arguments['--moments']:
            print >> sys.stderr, "--moments not implemented for -S Xt"
            sys.exit(1)

        omegaMesh=np.linspace(0,float(arguments['--frequencyrange']),int(arguments['--frequencies']))
        #if arguments['--modelfile']:
        #    origomegaMesh=np.loadtxt(arguments['--modelfile'])[:,0]
        #    with warnings.catch_warnings():
        #        warnings.simplefilter("ignore")
        #    origmodel=-np.loadtxt(arguments['--modelfile'])[0:,-1]
        #    f=interpolate.interp1d(origomegaMesh,origmodel,kind='linear',fill_value=1e-15,bounds_error=False)
        #    with warnings.catch_warnings():
        #        warnings.simplefilter("ignore")
        #        oneoveromegaMesh=np.where(omegaMesh==0.0,0.0,1/omegaMesh)
        #    model=-f(omegaMesh)*oneoveromegaMesh*2.0+1e-8
        #else:
        model=omegaMesh*0.0+1.0

        model*=float(arguments["--modelnorm"])/np.trapz(model,omegaMesh)
        B=Bryan(tauMesh,omegaMesh,Beta,model,kernel='realbosonic',norm=-np.trapz(G[:,0],tauMesh)/Beta)

    elif arguments['--statistics']=='F':
        omegaMesh=np.linspace(-float(arguments["--frequencyrange"])/2.0,float(arguments["--frequencyrange"])/2.0,int(arguments["--frequencies"]))
        Gtau=np.loadtxt(arguments["<datafile>"],np.float128)
        tauMesh=Gtau[:,0]
        Beta=tauMesh[0]+tauMesh[-1]

        if arguments["--sigma"]:
            Ns=Gtau.shape[1]-1
            sigmaFromFile=False
        else:
            Ns=(Gtau.shape[1]-1)/2
            sigmapm2=[[] for i in range(Ns)]
            sigmaFromFile=True
        G=Gtau[:,1::(2 if sigmaFromFile else 1)]
        if sigmaFromFile:
            sigmapm2=1./Gtau[:,2::2]**2
        else:
            sigmapm2=G*0.0+1.0/float(arguments['--sigma'])**2

        #if arguments['--modelfile']:
        #    origomegaMesh=np.loadtxt(arguments['--modelfile'])[:,0]
        #    with warnings.catch_warnings():
        #        warnings.simplefilter("ignore")
        #    origmodel=-np.loadtxt(arguments['--modelfile'])[0:,-1]
        #    f=interpolate.interp1d(origomegaMesh,origmodel,kind='linear',fill_value=1e-15,bounds_error=False)
        #    with warnings.catch_warnings():
        #        warnings.simplefilter("ignore")
        #        oneoveromegaMesh=np.where(omegaMesh==0.0,0.0,1/omegaMesh)
        #    model=f(omegaMesh)+1e-6
        #el
        if arguments['--moments']:
            model=0.0*omegaMesh
            for momentset in arguments['--moments']:
                m=np.array(momentset.split(','),np.float)
                if len(m)<3:
                    m=np.hstack([[1.0],m])
                print "using moments %s"%m
                model+=m[0]*np.exp(-(omegaMesh-m[1])**2/2.0/m[2]**2)/np.sqrt(2*np.pi*m[2]**2)+1e-8
        else:
            model=np.ones(omegaMesh.shape, np.float64)
        model*=float(arguments["--modelnorm"])/np.trapz(model,omegaMesh)
        print "Done"

        ######### MAXENTING #############

        B=Bryan(tauMesh, omegaMesh, Beta, model)
    else:
        print >> sys.stderr, "only real bosonic (X) and fermionic (F) continuation for now"
        sys.exit(1)

    outputA=[]
    print "shape(G)=%s, shape(sigmapm2=%s)"%(G.shape,sigmapm2.shape)
    for i in range(G.shape[1]):
        B.initG(G[:,i], sigmapm2[:,i])
        (aMax,AMax,S,L,lam)=B.approxMax()
        F=MeasureFunctor(tauMesh,omegaMesh,(0.5*np.sum(np.log(aMax/(aMax+lam)))+aMax*S-L)-2, aMax, AMax, S, L, lam,G[:,i])
        if gnuplot:
            F.plot("set arrow 1 from %s, graph 0 to %s, graph 1 nohead"%(np.log(aMax), np.log(aMax)))


        B.screen(np.log(aMax), 0.1, AMax, F.add)
        B.screen(np.log(aMax), -0.1, AMax, F.add)
        ABryan=F.finalize()
        if gnuplot:
            F.plot2.plot(Data(omegaMesh,ABryan, with_="l", title="Bryan's MaxEnt"), Data([0],[1e400/1e-99]), Data(omegaMesh, AMax, with_="l", title="Classic MaxEnt", inline=True))
        if arguments['--statistics']=='F':
            outputA.append(np.copy(ABryan))
        elif arguments['--statistics']=='X':
            outputA.append(ABryan*B.omegaMesh*np.pi/2.0*Ww[0,1])
        elif arguments['--statistics']=='Xt':
            outputA.append(ABryan*B.omegaMesh*np.pi/2.0*np.trapz(G[:,0],tauMesh)/Beta)
        if gnuplot:
            print "Press return to go on"
            sys.stdin.readline()

    ####### WRTINTING FILE #########

    print "writing output"
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
    print "done"
