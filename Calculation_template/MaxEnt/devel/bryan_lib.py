def Kernel_F_tau( tau, wr, Beta ):
    #
    Nt = len(tau)
    Nw = len(wr)
    K = np.zeros( (Ntau,Nw), dtype=precD )
    #
    for it in range(Nt):
        for iw in range(Nw):
            #
            t = tau[it]
            w = wr[iw]
            #
            K[it,iw] = np.exp((Beta/2-t)*w) / ( np.exp(Beta/2*w) + np.exp(-Beta/2*w) )
            #
    return K

def Kernel_B_tau( tau, wr, Beta ):
    #
    Nt = len(tau)
    Nw = len(wr)
    K = np.zeros( (Ntau,Nw), dtype=precD )
    #
    for it in range(Nt):
        for iw in range(Nw):
            #
            t = tau[it]
            w = wr[iw]
            #
            if abs(w) < 1e-8:
                K[it,iw] = ( np.cosh((Beta/2-t)*w) / Beta )
            elif Beta*w >700.0:
                K[it,iw] = ( np.exp(-t*w) + np.exp(-w*(Beta-t)) ) * (w/2.0)
            else:
                K[it,iw] = ( np.cosh((Beta/2-t)*w) / np.sinh(Beta/2*w) ) * (w/2.0)
            #
    return K

def Kernel_F_mats( wm, wr, Beta ):
    #
    Nm = len(wm)
    Nw = len(wr)
    K = np.zeros( (2*Nm-1,Nw), dtype=precC )
    #
    for im in range(Nm):
        for iw in range(Nw):
            #
            m = wm[im]
            w = wr[iw]
            #
            K[Nm-1+im,iw] = 1.0 / ( 1j*m - w  )
            #
        if im>0: K[Nm-1-im,:] = np.conjugate(K[Nm-1+im,:])
    return K

def Kernel_B_mats( wm, wr, Beta ):
    #
    Nm = len(wm)
    Nw = len(wr)
    K = np.zeros( (2*Nm-1,Nw), dtype=precC )
    #
    for im in range(Nm):
        for iw in range(Nw):
            #
            m = wm[im]
            w = wr[iw]
            #
            K[Nm-1+im,iw] = w / ( 1j*m +s w  )
            #
        if im>0: K[Nm-1-im,:] = np.conjugate(K[Nm-1+im,:])
    return K

class Bryan:
    #
    # Initialization of the Kernel depending on Field and statistics
    def __init__(self, axis, wr, Beta, statistics, kernel, convergence ):
        #
        self.convergence = convergence
        self.axis = axis            # input axis
        self.wr = wr                # output real frequency mesh
        self.dw = abs(wr[1]-wr[0])
        self.Beta = Beta*1.0
        self.dTau = abs(axis[1]-axis[0]) if (kernel=="tau") else 0.0
        #
        # Setting up Kernel
        if statistics == "F":
            if kernel == "tau":
                self.K = self.dw * Kernel_F_tau( axis, wr, Beta )
            elif kernel == "mats":
                self.K = self.dw * Kernel_F_mats( axis, wr, Beta )
            else:
                print >> sys.stderr, "Bryan.__init__ kernel not allowed "
                sys.exit(1)
        elif statistics == "B":
            if kernel == "tau":
                self.K = self.dw * Kernel_B_tau( axis, wr, Beta )
            elif kernel == "mats":
                self.K = self.dw * Kernel_B_mats( axis, wr, Beta )
            else:
                print >> sys.stderr, "Bryan.__init__ kernel not allowed "
                sys.exit(1)
        else:
            print >> sys.stderr, "Bryan.__init__ statistics not allowed "
            sys.exit(1)
        #
        # Computing SVD. Last line of pg.167 left column
        ( V, Sigma, Ut ) = np.linalg.svd( self.K )
        #
        # Reducing the number of singular values. Top lines of pg.167 right column
        good = [i for i in range (Sigma.shape[0]) if Sigma[i]>Sigma[-1]*100 ]
        self.s = len(good)
        self.TSigma = Sigma[good]
        self.TU = Ut.transpose()[:,good]
        self.TVt = V.transpose()[good,:]
        print("Dimension of kernel:%s; Dimension of singular space:%s"%(self.K.shape, self.s))
    #
    # Initialization of the model and M matrix
    def initField( self, Field, sigma, model ):
        #
        self.model = np.array( model, dtype=precD )
        self.oneovermodel = np.where( model>0.0, 1.0/model, 0.0 )
        self.Field = np.array( Field, dtype=precD )
        self.sigmapm2 = np.array( sigmapm2, dtype=precD )
        # Definition of M before Eq.11 of pg.168 left column
        self.M = np.dot(np.dot(np.diag(self.TSigma),np.dot(np.dot(self.TVt,np.diag(self.sigmapm2)),self.TVt.transpose())),np.diag(self.TSigma))
    #
    # MaxEnt algorithm
    def maxent( self, alpha, initA ):
        #
        print("Calculating maxent with alpha=%s"%alpha)
        #
        # Initialization of the spectral function
        A = np.copy(initA) + 1e-8
        # Definition of new variable "u" before Eq.9 of pg.167 right column
        u = np.dot( self.TU.transpose(), np.log(A*self.oneovermodel) )
        #
        # Loop over the solution of Eq.11 of pg.168 left column
        err = 1.0
        iter = 0
        while err > self.convergence:
            #
            # Field derived with the Kernel with fewer singular values (K=V*Sigma*Ut) and the trial spectral function
            Field_ = np.dot(self.TVt.transpose(),np.dot(np.diag(self.TSigma),np.dot(self.TU.transpose(),A)))
            # Variation (because A is postive?) between the real and trial Field. dL/dF of the paper.
            dF = (Field_ + self.Field) * self.sigmapm2
            # Definition of new variable "g" in Eq.10 of pg.168 left column
            g = np.dot(np.diag(self.TSigma),np.dot(self.TVt,dF))
            # Definition of new variables "K" and "MK" before Eq.11 of pg.168 left column
            K = np.dot(np.dot(self.TU.transpose(),np.diag(A)),self.TU)
            MK = np.dot(self.M,K)
            # "du" is the solution of Eq.11 of pg.168 left column
            du = np.linalg.solve( (alpha)*np.eye(MK.shape[0]) + MK, -alpha*u-g )
            # control variable which sets the convergence criteria. Right after Eq.12 of pg.168 right column
            err = np.dot(np.dot(du,K),du)


            #
            if(verb): print("iteration # %s, error %s"%(iter,err))
            #
            # increment in the "u" veriable
            u += du
            # update of the trial spectral function of Eq.9 of pg.167 right column
            A = self.model*np.exp(np.dot(self.TU,u))
            iter = iter+1
            print("mdu=%s it=%s"%(mdu,iter))
            #





































#
