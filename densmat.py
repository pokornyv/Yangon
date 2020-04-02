################################################################
# Yangon                                                       #
# transport properties from ab-initio DFT results              #
# Copyright (C) 2019-2020 Vladislav Pokorny; pokornyv@fzu.cz   #
# homepage: github.com/pokornyv/Yangon                         #
# method described in J. Chem. Phys. 126, 174101 (2007).       #
################################################################

## list of functions
# DensmatEq
# FpqMatrices
# Densmat
# ChargeDensmat
# SpinDensmat
# MagMoment
# Nelectron
# LowdinCharges
# PlotSigmaDep
# PlotFermiDep
# FindRSigma
# FindEFermi
# FindRSigmaB
# FindEFermiB
# FindRSigmaN
# FindEFermiN
# DensmatDiff
# PulayMixer


from config import *
from selib import TransformGammaSOCdm,SaveDensmat_tm2ait,UpdateSelfEnergy
from scipy.linalg.lapack import zgeev

def DensmatEq(occn_A,mos_A,Spin=''):
	''' constructs the equilibrium density matrix from TM output '''
	Dmat_A = multi_dot([mos_A,sp.diag(occn_A),sp.conj(mos_A.T)])
	if not sp.allclose(sp.conj(Dmat_A.T), Dmat_A):
		print(' - Error: Dmat0 is not Hermitean!')
	print(' - Tr Dmat0'+Spin+' = {0: .8f}'.format(sp.real(sp.trace(Dmat_A))),flush='True')
	return Dmat_A


def FpqMatrices(mu,bias,Eext_A):
	''' Fpq matrices as defined in Ref. 
	mu is the avg. chemical potential
	sp.log correctly works out the imaginary part arctan2(y,x) '''
	N = len(Eext_A)
	muL = mu-0.5*bias
	muR = mu+0.5*bias
	# I_ij = 1/(z_i-z*_j), Ln_ij = Log1[i]-Log2[j]
	I_A = 1.0/sp.add(sp.outer(Eext_A,sp.ones(N)),-sp.outer(sp.ones(N),sp.conj(Eext_A)))
	Log1L_A = sp.log(-Eext_A+muL)
	Log2L_A = sp.log(sp.conj(Eext_A)-muL)
	LnL_A = sp.add(sp.outer(Log1L_A,sp.ones(N)),-sp.outer(sp.ones(N),Log2L_A))
	FpqL_A = I_A*(LnL_A-1.0j*sp.pi)/(2.0*sp.pi)
	if not sp.allclose(sp.conj(FpqL_A.T), FpqL_A): 
		print(' - Error: Fpq(L) is not Hermitean')
	if bias != 0.0:
		Log1R_A = sp.log(-Eext_A+muR)
		Log2R_A = sp.log(sp.conj(Eext_A)-muR)
		LnR_A = sp.add(sp.outer(Log1R_A,sp.ones(N)),-sp.outer(sp.ones(N),Log2R_A))
		FpqR_A = I_A*(LnR_A-1.0j*sp.pi)/(2.0*sp.pi)
		if not sp.allclose(sp.conj(FpqR_A.T), FpqR_A): 
			print(' - Error: Fpq(R) is not Hermitean')
	else: FpqR_A = sp.copy(FpqL_A)
	return FpqL_A,FpqR_A


def Densmat(avgmu,bias,Eext_A,Bmat_A,GammaL_A,GammaR_A,FullMat=True,Spin=''):
	''' construct the non-equilibrium density matrix 
	if FullMat=False, saves some time by calculating only the number of electrons '''
	t = time()
	#FpqL_A = FpqMatrix2(avgmu-0.5*bias,Eext_A)
	#FpqR_A = FpqMatrix2(avgmu+0.5*bias,Eext_A)
	FpqL_A,FpqR_A = FpqMatrices(avgmu,bias,Eext_A)
	#print('Fpq calculation: {0: .8f}s'.format(time()-t))
	#print(' - Tr FpqL',sp.trace(FpqL_A))
	#print(' - Tr FpqR',sp.trace(FpqR_A))
	if cmos:
		# arrays of two matrices 2Nx2N, up and dn, we need their sum for Mpq
		# this is a bottleneck of the calculation of the density metrix
		# although the real bottleneck of any calculation is the Hext diagonalization
		GammaL2_A, GammaR2_A = TransformGammaSOCdm(GammaL_A,GammaR_A,Bmat_A)
		#print('Gammas transformed in ',time()-t)
		## this is done element-by-element:
		MpqL_A = GammaL2_A*FpqL_A
		MpqR_A = GammaR2_A*FpqR_A
	else:
		Binv_A  = inv(Bmat_A)
		BinvH_A = sp.conj(Binv_A.T)
		## Gamma2*Fpq, done element-by-element:
		MpqL_A  = multi_dot([Binv_A*GammaL_A,BinvH_A])*FpqL_A
		MpqR_A  = multi_dot([Binv_A*GammaR_A,BinvH_A])*FpqR_A
	Mpq_A = MpqL_A + MpqR_A
	#print(' - Tr Mpq', sp.trace(Mpq_A))
	if not sp.allclose(sp.conj(Mpq_A.T), Mpq_A):
		print(' - Warning:  Mpq is not Hermitean!')
		DR = sp.amax(sp.fabs(sp.real(sp.conj(Mpq_A.T)-Mpq_A)))
		DI = sp.amax(sp.fabs(sp.imag(sp.conj(Mpq_A.T)-Mpq_A)))
		print(' Maximum deviation: real: {0: 5e}, imag: {1: 5e}'.format(DR,DI))
		exit(1)
	if cmos:
		NBF = int(len(Bmat_A)/2)
		Bup_A = Bmat_A[:NBF,:] 	## spin-up block of B (top), Nx2N
		Bdn_A = Bmat_A[NBF:,:]	## spin-dn block of B (bottom), Nx2N
		BupH_A = sp.conj(Bup_A.T)
		BdnH_A = sp.conj(Bdn_A.T)
		if FullMat: ## we calculate the whole matrix
			Dmatuu_A =  multi_dot([Bup_A,Mpq_A,BupH_A])
			Dmatdd_A =  multi_dot([Bdn_A,Mpq_A,BdnH_A])
			Dmatud_A =  multi_dot([Bup_A,Mpq_A,BdnH_A])
			Dmatdu_A =  multi_dot([Bdn_A,Mpq_A,BupH_A])
			Densmat_A = sp.block([[Dmatuu_A,Dmatud_A],[Dmatdu_A,Dmatdd_A]])
			TrD = sp.real_if_close(sp.trace(Dmatuu_A)+sp.trace(Dmatdd_A))
		else: ## calculate only the trace, Tr(AB)=sum(A*B.T)
			Densmat_A = 0
			Aup_A = sp.dot(Bup_A,Mpq_A)
			Adn_A = sp.dot(Bdn_A,Mpq_A)
			TrDup = sp.sum(Aup_A*sp.conj(Bup_A))
			TrDdn = sp.sum(Adn_A*sp.conj(Bdn_A))
			TrD = sp.real_if_close(TrDup+TrDdn)
		#print(' - Tr Duu',sp.trace(Dmatuu_A))
		#print(' - Tr Dud',sp.trace(Dmatud_A))
		#print(' - Tr Ddu',sp.trace(Dmatdu_A))
		#print(' - Tr Ddd',sp.trace(Dmatdd_A))
	else: ## not cmos
		Densmat_A = multi_dot([Bmat_A,Mpq_A,sp.conj(Bmat_A.T)])
		if NSpin == 1: Densmat_A *= 2.0 ## two spin channels
		TrD = sp.real_if_close(sp.trace(Densmat_A))
		if NSpin == 2: 
			NBF = int(len(Bmat_A)/2)
			print(' - Tr DmatUp: {0: .8f}, Tr DmatDn: {1: .8f}'\
			.format(sp.real(sp.trace(Densmat_A[:NBF,:NBF])),sp.real(sp.trace(Densmat_A[NBF:,NBF:]))))
	print(' - Tr Dmat'+Spin+' = {0: .8f}'.format(sp.real(TrD)),flush='True')
	#if not sp.allclose(sp.conj(Densmat_A.T), Densmat_A):
	#	print(" - Error: Density matrix is not Hermitean!")
	#	exit(1)
	print(' - Densmat: Dmat calculation: {0: .8f}s'.format(time()-t))
	if FullMat: return Densmat_A, TrD
	else: return TrD


def ChargeDensmat(Dmat_A):
	''' calculates the charge density matrix Duu+Ddd
	spin-orbit case only '''
	N = int(len(Dmat_A)/2)
	return Dmat_A[:N,:N]+Dmat_A[N:,N:]


def SpinDensmat(Dmat_A):
	''' calculates the spin density matrices Dud+Ddu, i(Dud-Ddu), Duu-Ddd
	spin-orbit case only  '''
	N = int(len(Dmat_A)/2)
	Mx_A = Dmat_A[N:,:N]+Dmat_A[:N,N:]
	My_A = 1.0j*(Dmat_A[N:,:N]-Dmat_A[:N,N:])
	Mz_A = Dmat_A[:N,:N]-Dmat_A[N:,N:]
	print(' - Tr(DmatUp): {0: .8f}'.format(sp.real(sp.trace(Dmat_A[:N,:N]))))
	print(' - Tr(DmatDn): {0: .8f}'.format(sp.real(sp.trace(Dmat_A[N:,N:]))))
	return Mx_A, My_A, Mz_A


def MagMoment(Dmat_A,tag=''):
	''' calculate the magnetic moment from density matrix 
	the one on the end of the output line prohibits Turbomole's 
	cgnce script to go crazy as it looks for lines with 6 entries to find Etot '''
	if cmos:
		## caculating the spin expectation value:
		Mx_A, My_A, Mz_A = SpinDensmat(Dmat_A)
		Sx = sp.real(sp.trace(Mx_A))/2.0
		Sy = sp.real(sp.trace(My_A))/2.0
		Sz = sp.real(sp.trace(Mz_A))/2.0
		S2 = sp.sqrt(Sx**2+Sy**2+Sz**2)
		print('\nSpin expectation value:')
		print('  <Sx> = {0: .5f}  <Sy> = {1: .5f}  <Sz> = {2: .5f}  S2 = {3: .5f}'.format(Sx,Sy,Sz,S2))
		if tm2trans: print('{0: 5d}\t{1: .6f}\t{2: .6f}\t{3: .6f}\t{4: .6f}\t:OUT_SPIN'\
		.format(tm2iter+1,Sx,Sy,Sz,S2)+tag+' 1')
		return Sx, Sy, Sz
	elif NSpin == 2:
		N = int(len(Dmat_A)/2)
		Sz = sp.real(sp.trace(Dmat_A[0:N,0:N])-sp.trace(Dmat_A[N:,N:]))/2.0
		print('\nSpin expectation value:')
		print('  <Sz> = {0: .6f}'.format(Sz))
		if tm2trans: print('{0: 5d}\t{1: .6f}\t{2: .6f}\t{3: .6f}\t{4: .6f}\t:OUT_SPIN'\
		.format(tm2iter+1,0.0,0.0,Sz,Sz)+tag+' 1')
		return 0.0, 0.0, Sz
	else:
		return 0.0, 0.0, 0.0


def Nelectron(ReSigma,Ham_A,SEl_A,SEr_A,EFermi,bias):
	''' function that calculates number of electrons from the non-equilibrium
	density matrix as a function of ReSigma 
	used as a lambda function for root finders 
	diagonalization of Hext is the bottleneck of the calculation '''
	SEl_A, SEr_A = UpdateSelfEnergy(ReSigma,SEl_A,SEr_A)
	GammaL_A = -2.0*sp.imag(SEl_A)	## these do not change with ReSigma, should be global
	GammaR_A = -2.0*sp.imag(SEr_A)
	ExtHam_A = Ham_A + sp.diag(SEl_A+SEr_A)
	#Eext_A, BmatL_A, BmatR_A = eig(ExtHam_A,left=True,right=True)
	#Eext_A, BmatL_A, BmatR_A, zout = zgeev(ExtHam_A,compute_vl=0, compute_vr=1, overwrite_a=0)
	t = time()
	Eext_A, BmatR_A = eig(ExtHam_A,left=False,check_finite=True)
	print(' - Nelectron: Eigval calculation: {0: .8f}s'.format(time()-t))
	TrD = Densmat(EFermi,bias,Eext_A,BmatR_A,GammaL_A,GammaR_A,FullMat=False)
	return sp.real(TrD) ## sp.real to change the data type so root solver does not complain


def LowdinCharges(Dmat_A,Atoms_L,AtomicSeg_A,NAtom):
	''' performs a Lowdin population analysis '''
	NBF = sp.sum(AtomicSeg_A) ## number of basis functions
	charges_A = sp.zeros(NAtom)
	spins_A   = sp.zeros(NAtom)
	ddiag_A   = sp.real(sp.diag(Dmat_A))
	NL = 0
	for i in range(NAtom):
		charges_A[i] = sp.sum(ddiag_A[NL:NL+AtomicSeg_A[i]])
		if cmos or NSpin == 2:
			cdn = sp.sum(ddiag_A[NL+NBF:NL+NBF+AtomicSeg_A[i]])
			spins_A[i] = charges_A[i] - cdn
			charges_A[i] += cdn
		NL += AtomicSeg_A[i]
	return charges_A, spins_A


def PlotSigmaDep(Ham_A,SEleft_A,SEright_A,EFermi,bias,ReSigma):
	''' plot the dependence of total electron number on Re Sigma '''
	Smax = 10.0
	dS = 0.1
	ReS = -Smax
	while ReS < Smax:
		NE = Nelectron(ReS,Ham_A,SEleft_A,SEright_A,EFermi,bias)-NElec
		print('{0: .8f}\t{1: .8f}\t{2: .8f}'.format(ReS,sp.real(NE),sp.imag(NE)))
		ReS += dS
	

def PlotFermiDep(Eext_A,BmatR_A,GammaL_A,GammaR_A,EF,bias,ReSigma):
	''' plot the dependence of total electron number on Fermi energy '''
	Emax = 10.0
	dE = 0.1
	EF = -Emax
	while EF < Emax:
		NE = Densmat(EF,bias,Eext_A,BmatR_A,GammaL_A,GammaR_A,FullMat=False)-NElec
		print('{0: .8f}\t{1: .8f}\t{2: .8f}'.format(EF,sp.real(NE),sp.imag(NE)))
		EF += dE


##### MINPACK root solvers #######
def FindRSigma(Ham_A,SEleft_A,SEright_A,EFermi,bias,ReSigma):
	''' search for ReSigma to tune TrD to number of electrons '''
	eqn = lambda x: Nelectron(x,Ham_A,SEleft_A,SEright_A,EFermi,bias)-NElec
	if sp.fabs(eqn(ReSigma))>1e-5:
		sol = root(eqn,[ReSigma],method='hybr',options={'xtol':1e-5})
		if sol.success:
			ReS = sol.x[0]
			print('   ReSigma1 = {0: .8f} H, dn: {1: .5e}'\
			.format(ReS*ImSigma1,ReS-NElec),flush=True)
		else: 
			if sp.fabs(eqn(sol.x[0]))<1e-5:
				print('Warning: failed to calculate ReSigma, MINPACK error message:')
				print(sol.message)
				print('Residue is small so we continue...')
			else:
				print('Error: failed to calculate ReSigma, MINPACK error message:')
				print(sol.message)
				if tm2trans: SaveDensmat_tm2ait(sp.zeros([2,2]),'ait2tm')
				exit(1)
	else: 
		ReS = ReSigma
	return ReS


def FindEFermi(Eext_A,BmatR_A,GammaL_A,GammaR_A,EF,bias,ReSigma):
	''' search for EFermi to tune TrD to number of electrons, uses root() '''
	eqn = lambda x: Densmat(EF,bias,Eext_A,BmatR_A,GammaL_A,GammaR_A,FullMat=False)-NElec
	if sp.fabs(eqn(EF))>1e-6:
		sol = root(eqn,[EF],method='hybr',options={'xtol':1e-5})
		if sol.success:
			EF = sol.x[0]
			print('   EFermi = {0: .8f} H, dn: {1: .5e}'.format(EF,sol.fun))
		else: 
			if sp.fabs(eqn(sol.x[0]))<1e-5:
				print('Warning: failed to calculate EFermi, MINPACK error message:')
				print(sol.message)
				print('Residue is small so we continue...')
			else:
				print('Error: failed to calculate EFermi, MINPACK error message:')
				print(sol.message)
				if tm2trans: SaveDensmat_tm2ait(sp.zeros([2,2]),'ait2tm')
				exit(1)
	return EF


###### bisections ################
def FindRSigmaB(Ham_A,SEleft_A,SEright_A,EFermi,bias,ReSigma):
	''' search for ReSigma to tune TrD to number of electrons using bisection
	we assume there is only one solution '''
	eqn = lambda x: Nelectron(x,Ham_A,SEleft_A,SEright_A,EFermi,bias)-NElec
	'''
	Emax = 2.0
	a = b = 0.0
	while True:
		if eqn(ReSigma)*eqn(ReSigma+Emax)<0.0: 
			a,b = ReSigma,ReSigma+Emax
			break
		elif eqn(ReSigma-Emax)*eqn(ReSigma)<0.0: 
			a,b = ReSigma-Emax,ReSigma
			break
		else:
			Emax += 1.0
		print('Emax = {0: .5f}, [{1: .5f} :{2: .5f}]'.format(Emax,a,b))
	'''
	a = -50.0
	b =  10.0
	print('   Bracketing interval: [{0: .3f} :{1: .3f}]'.format(a,b))
	if sp.fabs(eqn(ReSigma))>1e-4:
		ReS = brentq(eqn,a,b,xtol=1e-4,rtol=1e-4)
	else: 
		ReS = ReSigma
	return ReS


def FindEFermiB(Eext_A,BmatR_A,GammaL_A,GammaR_A,EF,bias,ReSigma):
	''' search for EFermi to tune TrD to number of electrons using bisection
	we assume there is only one solution '''
	eqn = lambda x: Densmat(EF,bias,Eext_A,BmatR_A,GammaL_A,GammaR_A,FullMat=False)-NElec
	'''
	Emax = 1.0
	a = b = 0.0
	while True:
		if eqn(EF)*eqn(EF+Emax)<0.0: 
			a,b = EF,EF+Emax
			break
		elif eqn(EF-Emax)*eqn(EF)<0.0: 
			a,b = EF-Emax,EF
			break
		else:
			Emax += 1.0
		print('Emax = {0: .5f}, [{1: .5f} :{2: .5f}]'.format(Emax,a,b))
	'''
	a = -50.0
	b =  50.0
	print(eqn(a),eqn(b))
	print('   Bracketing interval: [{0: .3f} :{1: .3f}]'.format(a,b))
	if sp.fabs(eqn(EF))>1e-5:
		EF = brentq(eqn,a,b,xtol=1e-5,rtol=1e-5)
	return EF


###### Newton solvers ############
def FindRSigmaN(Ham_A,SEleft_A,SEright_A,EFermi,bias,ReS):
	''' search for ReSigma to tune TrD to number of electrons '''
	eqn = lambda x: Nelectron(x,Ham_A,SEleft_A,SEright_A,EFermi,bias)-NElec
	if sp.fabs(eqn(ReS))>1e-4:
		ReS = newton(eqn,ReS,fprime=None,tol=1e-4,maxiter=50)
	return ReS


def FindEFermiN(Eext_A,BmatR_A,GammaL_A,GammaR_A,EF,bias,ReSigma):
	''' search for EFermi to tune TrD to number of electrons, uses brentq() '''
	eqn = lambda x: Densmat(EF,bias,Eext_A,BmatR_A,GammaL_A,GammaR_A,FullMat=False)-NElec
	if sp.fabs(eqn(EF))>1e-5:
		EF = newton(eqn,EF,fprime=None,tol=1e-5,maxiter=50)
	return EF

##################################

def DensmatDiff(Dmat_A,DmatOld_A):
	''' calculate deviation of two matrices '''
	Dev_A = (Dmat_A-DmatOld_A)**2
	return sp.sqrt(sp.sum(Dev_A))/len(Dmat_A)


def PulayMixer(DMats_L,NMat):
	''' use Pulay extrapolation to calculate a new guess for the density matrix 
	the density matrices are in the Turbomole's non-orthogonal basis '''
	N = len(DMats_L)
	if N != MaxMix+1:
		print(' - Error: insufficient data for Pulay mixer, returning the last matrix')
		return DMats_L[0]
	ErrorVectors_A = sp.zeros([NMat**2,MaxMix],dtype='complex')
	for i in range(MaxMix):
		## construct error vectors N(i+1)-N(i) (begin from zero)
		ErrorVectors_A[:,i] = sp.array(sp.ravel(DMats_L[i+1]-DMats_L[i]))
	## last row and last column belong to the Lagrange multiplier
	## that keeps the sum of coeffs at one
	Overlap_A = sp.ones([MaxMix+1,MaxMix+1],dtype='complex')
	for i,j in product(range(MaxMix),repeat=2):
		Overlap_A[i,j] = sp.dot(ErrorVectors_A[:,i],ErrorVectors_A[:,j])
	Overlap_A[-1,-1] = 0.0 
	#print('Overlap of error vectors:')
	#for i in range(MaxMix+1):
	#	for j in range(MaxMix+1):
	#		print('{0: .6f}  '.format(Overlap_A[i,j]),end='')
	#	print('')
	RHS_A = sp.zeros(MaxMix+1)
	RHS_A[-1] = 1.0
	try:
		Coeffs_A = solve(Overlap_A,RHS_A)
	except LinAlgError:
		print(' - Cannot find Pulay weights, matrix is probably singular.')
		Coeffs_A = sp.zeros(MaxMix+1)
		Coeffs_A[0] = 1.0
	print('  Pulay weights:  ',end='')
	for i in range(MaxMix): print('{0: .6f}  '.format(sp.real(Coeffs_A[i])),end='')
	print('       Sum: {0: .3f}'.format(sp.real(sp.sum(Coeffs_A[:-1]))))
	print('  Lagrange miltiplier: Re L = {0: .6f}, Im L = {1: .6f}'.format(sp.real(Coeffs_A[-1]),sp.imag(Coeffs_A[-1])))
	DmatPulay_A = sp.zeros_like(DMats_L[0]) ##<class 'numpy.ndarray'>
	for i in range(MaxMix):
		DmatPulay_A+=Coeffs_A[i]*DMats_L[i]
	return DmatPulay_A, Coeffs_A


"""
def FpqMatrix(mu,Eext_A):
	''' the Fpq matrix as defined in Ref. [x] 
	mu is the chemical potential of the lead 
	sp.log correctly works out the imaginary part arctan2(y,x) '''
	N = len(Eext_A)
	Log1_A = sp.log(      -(Eext_A -mu))
	Log2_A = sp.log(sp.conj(Eext_A)-mu )
	# I_ij = 1/(z_i-z*_j), Ln_ij = Log1[i]-Log2[j]
	I_A = 1.0/sp.add(sp.outer(Eext_A,sp.ones(N)),-sp.outer(sp.ones(N),sp.conj(Eext_A)))
	Ln_A = sp.add(sp.outer(Log1_A,sp.ones(N)),-sp.outer(sp.ones(N),Log2_A))
	Fpq_A = I_A*(Ln_A-1.0j*sp.pi)/(2.0*sp.pi)
	if not sp.allclose(sp.conj(Fpq_A.T), Fpq_A): 
		print(" - Error: Fpq is not Hermitean!")
	return Fpq_A

def FpqMatrix2(mu,Eext_A):
	''' the Fpq matrix as defined in Ref. [x] 
	mu is the chemical potential of the lead 
	backup that calculates separately real and imaginary parts of the complex logarithm '''
	N = len(Eext_A)
	eps_A = mu-sp.real(Eext_A)
	eta_A =   -sp.imag(Eext_A)
	LnArg_A = eps_A**2 + eta_A**2
	Ln_A = sp.log(sp.outer(LnArg_A,1.0/LnArg_A))
	AtanP_A = sp.arctan2(eta_A, eps_A)
	AtanM_A = sp.arctan2(eta_A,-eps_A)
	# I_ij = 1/(z_i-z*_j), A_ij = AtanP_A[i]+AtanM_A[j]
	A_A = sp.add(sp.outer(AtanP_A,sp.ones(N)),-sp.outer(sp.ones(N),AtanM_A))
	I_A = 1.0/sp.add(sp.outer(Eext_A,sp.ones(N)),-sp.outer(sp.ones(N),sp.conj(Eext_A)))
	Fpq_A = I_A*(0.5*Ln_A-1.0j*(sp.pi - A_A))/(2.0*sp.pi)
	if not sp.allclose(sp.conj(Fpq_A.T), Fpq_A): 
		print(" - Error: Fpq is not Hermitean!")
	return Fpq_A

"""

## densmat.py END

