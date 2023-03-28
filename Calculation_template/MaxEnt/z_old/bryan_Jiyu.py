#!/usr/bin/env python
# Written by Lewin Boehnke (lboehnke@physnet.uni-hamburg.de), 2009-2015, consulted by Hugo Strand, 2015

try:
    import Gnuplot
    from Gnuplot import Data
    gnuplot=True
except ImportError:
    gnuplot=False
    print('-sdnifos-')

import numpy as np
import warnings
import itertools
import os
import sys
import time
import copy


print('start')

def fermion_function(X,beta):
    n_mesh = []
    for x in X:
        if beta*x>20:
            n=0
        elif beta*x<-20:
            n=1
        else:
            n = 1/(np.exp(beta*x)+1)
        n_mesh.append(n)
        #print(x,beta*x,n)
    return (np.array(n_mesh))

def kernel_function(omega,tau,beta,kernel_mode):
    if kernel_mode == 'fermionic':
        return np.exp(-tau*omega)/(1.+np.exp(-Beta*omega)) if omega>0.0 else np.exp((Beta-tau)*omega)/(np.exp(Beta*omega)+1.)  
    if kernel_mode == 'realbosonic':
        return 1/beta if abs(omega)<1e-7 else omega/2*(np.exp(-omega*tau)+np.exp(-omega*(beta-tau)))/(1-np.exp(-omega*beta)) if omega>0.0 else omega/2*(np.exp(omega*tau)+np.exp(omega*(beta-tau)))/(np.exp(omega*beta)-1) 
    else:
        print('error: unknown mode')
        
        
        
class Bryan:
    def __init__(self, tauMesh, omegaMesh, Beta, model,K,TSigma,TU,TVt,s,mesh_scale, kernel='fermionic',norm=1.0):
        self.tauMesh=tauMesh
        self.omegaMesh=omegaMesh
        self.Beta=Beta*1.0
        self.deltaTau=tauMesh[1]-tauMesh[0]
        self.deltaOmega=omegaMesh[1]-omegaMesh[0]

        self.K=K
        self.TSigma=TSigma
        self.TU=TU
        self.TVt=TVt  
        self.s=s
        self.model=np.array(model,np.float64) # 1/range
        self.oneovermodel=np.where(model>0.0,1.0/model,0.0)

    def initG(self,G,sigmapm2):
        self.G=np.array(G,np.float64)
        self.sigmapm2=np.array(sigmapm2,np.float64)
        self.M=np.dot(np.dot(np.diag(self.TSigma),np.dot(np.dot(self.TVt,np.diag(self.sigmapm2)),self.TVt.transpose())),np.diag(self.TSigma)) #before (11)



    def maxent(self,alpha,initA=None): ## Q = alpha *S-L
        if any(initA==None):
            initA=np.array([1.0/(np.sqrt(2*np.pi)*0.5)*np.exp(-((x)/0.5)**2/2)+1e-3 for x in self.omegaMesh],np.float64) ###
        print("Calculating maxent with alpha=%s"%alpha)
        A=np.copy(initA)
        A+=1e-8
        u=np.dot(self.TU.transpose(),np.log(A*self.oneovermodel))   ## u = TU log(A/model)
        mu=1e-2
        time = 0
        while True:
            print('in loop')
            G_=np.dot(self.TVt.transpose(),np.dot(np.diag(self.TSigma),np.dot(self.TU.transpose(),A))) # G_=K.dot(A)
            dG=G_+self.G
            g=np.dot(np.diag(self.TSigma),np.dot(self.TVt,dG*self.sigmapm2)) # (10) and page 166 left column upper half
            K=np.dot(np.dot(self.TU.transpose(),np.diag(A)),self.TU) # before (B11)
            MK=np.dot(self.M,K)
            du=np.linalg.solve((alpha)*np.eye(MK.shape[0])+MK, -alpha*u-g) # (B11)
            mdu=np.dot(np.dot(du,K),du) # after (B12)
            if time%100==0:
                print('time=',time,mdu<1e-4,'mdu=',mdu)
            if mdu<1e-7: 
                break # there is an alternate criterion in Bryan's paper
            if mdu>1:
                du=np.linalg.solve((alpha+mu)*np.eye(MK.shape[0])+MK, -alpha*u-g) # (B12)
                mdu=np.dot(np.dot(du,K),du) #after (B12)
                time2=0
                while (abs(mdu-1)>1e-4):
                    mu*=1+(mdu-1)
                    du=np.linalg.solve((alpha+mu)*np.eye(MK.shape[0])+MK, -alpha*u-g) # (B12)
                    mdu=np.dot(np.dot(du,K),du) # after (B12)
                    time2+=1
                    if time2%100==0:
                        print('time=',time,'time2=',time2,' mdu=',mdu, 'abs',abs(mdu-1))
            u+=du
            A=self.model*np.exp(np.dot(self.TU,u)) # (9)
            time+=1

        lam=np.linalg.eigvalsh(MK)
        tmp=np.log(np.where(A>1e-6,A*(self.oneovermodel+1e-8),1e-6*self.oneovermodel))
        S=np.sum(A-model-(A*tmp)[tmp>-np.inf])
        L=0.5*np.sum((G_+self.G)**2*self.deltaTau*self.sigmapm2) # page 166 left column upper half
        print("Q=%s"%(alpha*S-L))
        return (A,S,L,lam,-G_)


    def approxMax(self,aMax=1.0):
        AMax=self.model.copy()
        while True:
            (AMax,S,L,lam,_)=self.maxent(aMax,AMax)
            tmp=np.sum(lam/(aMax+lam))/(-2*aMax*S)
            if abs(tmp-1)<1e-4: 
                print('break!',tmp)
                break
            aMax*=(tmp-1)*1.0+1
            print('----',aMax,np.sum(lam/(aMax+lam))/(-2*S))
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

    def __init__(self, tauMesh, omegaMesh, thresh, a, A, S, L, lam,G):

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
    print(arguments)

    if arguments['--statistics']=='F':
        kernel='fermionic'
        #### Read Green's function ####
        Gtau=np.loadtxt(arguments["<datafile>"],np.float64)
        G=copy.deepcopy(Gtau[:,1::1])  
        
        
        np.savetxt('test.out', abs(Gtau))
        print('path name: ',path)
        print('file name:',filename)
        print('G dimension', G.shape[1])
        
        #### define tau mesh ####
        tauMesh=Gtau[:,0]
        Beta=tauMesh[0]+tauMesh[-1]
        Ns=Gtau.shape[1]-1 # 1

        #### define omega mesh ####
        try:
            omega_scale=arguments["--omega_scale"]
        except:
            omega_scale="homogeneous"
            
        if omega_scale=='sinh':
            print('sinh omega mesh')
            omegaMesh=np.sinh(4*np.linspace(-1,1,int(arguments['--frequencies'])))/np.sinh(4)*float(arguments['--frequencyrange'])/2.0  
        elif omega_scale=='homogeneous':
            print('homogeneous omega mesh')
            omegaMesh=np.linspace(-float(arguments['--frequencyrange'])/2.0,float(arguments['--frequencyrange'])/2.0,int(arguments['--frequencies']))
        else:
            print("unknown omega mesh mode")


        #### decide sigma ####
        if arguments["--sigma"]:
            Ns=Gtau.shape[1]-1
            sigmaFromFile=False
        else:
            Ns=(Gtau.shape[1]-1)/2
            sigmapm2=[[] for i in range(Ns)]
            sigmaFromFile=True


        if sigmaFromFile:
            sigmapm2=1./Gtau[:,2::2]**2
        else:
            sigmapm2=G*0.0+1.0/float(arguments['--sigma'])**2

        #### assign K ####
        print("Calculating K...")
        time_start=time.time()
        K=np.array([[kernel_function(omega,tau,Beta,kernel) for omega in omegaMesh] for tau in tauMesh], np.float64)
        
        if omega_scale=='homogeneous':
            K*=omegaMesh[1]-omegaMesh[0]

        if omega_scale=='sinh':
            Delta=[omegaMesh[i+1]-omegaMesh[i] for i in range(len(omegaMesh)-1)]
            Delta.append(Delta[0])
    
            K*=Delta
        
        
        time_end=time.time()
        print('time spent assigning K:',time_end-time_start)

        time_spent=time_end-time_start
        f=open('timeitK','a')
        f.write('time spend '+str(time_spent)+"(s)\n")
        f.close()

        ############### SVD ###################
        print("Calculating SVD...")
        time_start=time.time()
        (V,Sigma,Ut)=np.linalg.svd(K)
        time_end=time.time()
        print('time spent SVD:',time_end-time_start)
    
        time_spent=time_end-time_start
        f=open('timeitSVD','a')
        f.write('time spend '+str(time_spent)+"(s)\n")
        f.close()

        print("Done")
        
        good=[i for i in range (Sigma.shape[0]) if Sigma[i]>Sigma[-1]*100]
        s=len(good)
        TSigma=Sigma[good]
        TU=Ut.transpose()[:,good]
        TVt=V.transpose()[good,:] 
    
    
        print("Done")
        print("Dimension of kernel:%s; Dimension of singular space:%s"%(K.shape, s))
        print('dimension TVt', TVt.shape)
        print('dimension TSigma', TSigma.shape)
        print('dimension TU', TU.shape)  


        #### smoothing Green's function ####
        error_cut=5e-4
        for dim in range(G.shape[1]):
            G_smooth=[(G[i-2,dim]+G[i-1,dim]+G[i,dim]+G[i+1,dim]+G[i+2,dim])/5 if 20<i<G.shape[0]-20 else G[i,dim] for i in range(G.shape[0])]
            G[:,dim]=G_smooth        
            if max(G[:,dim]) > -error_cut:
                index_start= [x>-error_cut for x in G[:,dim]].index(True)
                index_end  = len(G[:,dim])-[x>-error_cut for x in G[:,dim]][::-1].index(True)
                print('For dim ',dim,' index of extrapolation: ',index_start,index_end)
                k1=(np.log(abs(G[index_start,dim]))-np.log(abs(G[0,dim])))/(index_start-0)
                k2=(np.log(abs(G[-1,dim]))-np.log(abs(G[index_end,dim])))/(len(G[:,dim])-1-index_end)
                logy1 = np.log(abs(G[0,dim]))
                logy2 = np.log(abs(G[-1,dim]))
                x1=0
                x2=len(G[:,dim])-1
                x=(k1*x1-k2*x2+logy1-logy2)/(k1-k2)
                print('For dim ',dim,' index of intersection: ',x)
                for i in range(index_start,round(x)):
                    G[i,dim]=-np.exp(k1*(i-x1)+logy1)
                for i in range(round(x),index_end):
                    G[i,dim]=-np.exp(k2*(i-x2)+logy2)   
            else:
                print('For dim ',dim,' no extrapolation')
        
            

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
                print("using moments %s"%m)
                model+=m[0]*np.exp(-(omegaMesh-m[1])**2/2.0/m[2]**2)/np.sqrt(2*np.pi*m[2]**2)+1e-8
        else:
            model=np.ones(omegaMesh.shape, np.float64)
        model*=float(arguments["--modelnorm"])/np.trapz(model,omegaMesh)
        print("Done")

        ######### MAXENTING #############

        B=Bryan(tauMesh, omegaMesh, Beta, model,K,TSigma,TU,TVt,s,omega_scale)  
        
    else:
        print >> sys.stderr, "only fermionic (F) continuation for the new version"
        sys.exit(1)

    outputA=[]
    outputAMax=[]
    print("shape(G)=%s, shape(sigmapm2=%s)"%(G.shape,sigmapm2.shape))
    for i in range(G.shape[1]):
        print('#############addressing dimension',i)
        B.initG(G[:,i], sigmapm2[:,i])
        (aMax,AMax,S,L,lam)=B.approxMax()
        outputAMax.append(AMax)      
        
        print('Found max aMax,AMax,S,L,lam','aMax=',aMax)
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
            print("Press return to go on")
            sys.stdin.readline()

    ####### WRTINTING FILE #########

    print("writing output")
    outfile=arguments['--outfile'] if arguments['--outfile'] else arguments['<datafile>']+'_dos.dat'
    np.savetxt(outfile,np.hstack([B.omegaMesh.reshape(B.omegaMesh.shape[0],1),np.array(outputA).transpose()]))
    f=open(outfile,'a')
    f.write("#"+" ".join(sys.argv)+"\n")
    f.close()

    outfile=arguments['--outfile'] if arguments['--outfile'] else arguments['<datafile>']+'_dos.dat_max'
    np.savetxt(outfile,np.hstack([B.omegaMesh.reshape(B.omegaMesh.shape[0],1),np.array(outputAMax).transpose()]))
    f=open(outfile,'a')
    f.write("#"+" ".join(sys.argv)+"\n")
    f.close()


    outfile=arguments['--outfile'] if arguments['--outfile'] else arguments['<datafile>']+'_dos.dat_error'+arguments["--sigma"]
    np.savetxt(outfile,np.hstack([B.omegaMesh.reshape(B.omegaMesh.shape[0],1),np.array(outputA).transpose()]))
    f=open(outfile,'a')
    f.write("#"+" ".join(sys.argv)+"\n")
    f.close()


    #if arguments['--statistics']=='F':
    outfile2=arguments['--outfile'] if arguments['--outfile'] else arguments['<datafile>']+'_dos.dat_fd'
    n_mesh = fermion_function(B.omegaMesh,Beta)
    temp1= (np.array(outputA)*n_mesh).transpose()
    np.savetxt(outfile2,np.hstack([B.omegaMesh.reshape(B.omegaMesh.shape[0],1),temp1]))
    f=open(outfile2,'a')
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
    print(arguments['--statistics'])
    print(arguments)
