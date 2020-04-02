################################################################
# Yangon                                                       #
# transport properties from ab-initio DFT results              #
# Copyright (C) 2019-2020 Vladislav Pokorny; pokornyv@fzu.cz   #
# homepage: github.com/pokornyv/Yangon                         #
# method described in J. Chem. Phys. 126, 174101 (2007).       #
################################################################

## list of functions
# FindAtomsInPlane
# SelfEnergyVec
# ReadSEfile
# UpdateSelfEnergy
# MagneticVec
# TransformGamma
# TransformGammaSOCdm


from config import *
from iolib import SaveDensmat_tm2ait
#from scipy.linalg.blas import cgemm,zgemm

def FindAtomsInPlane(a1,a2,a3,coord_A):
	''' constructs a plane from two vectors and finds atoms laying in that plane '''
	atoms1_L = []
	atoms2_L = []
	atoms3_L = []
	## atoms are numbered from one, coord_A from zero
	Atom1_A = coord_A[a1-1]
	Atom2_A = coord_A[a2-1]
	Atom3_A = coord_A[a3-1]
	Vec1_A = Atom1_A - Atom2_A
	Vec2_A = Atom1_A - Atom3_A
	normal_A = sp.cross(Vec1_A,Vec2_A) ## normal of a plane (n1,n2,n3)
	normal_A /= sp.sqrt(sp.sum(normal_A**2)) ## normalize
	if sp.sum(sp.dot(normal_A,normal_A)) < 1e-2:
		print(' - Error: atoms {0: 2d}, {1: 2d}, {2: 2d} lie almost in a line.'.format(a1,a2,a3))
		print(' - Please select different atoms to define a plane.')
		if tm2trans: SaveDensmat_tm2ait(sp.zeros([2,2]),'ait2tm')
		exit(1)
	d = -sp.sum(normal_A*Atom1_A) ## fix of the plane. The plane is n1.x + n2.y + n3.z = d
	dist_A = sp.zeros(len(coord_A))
	for atom in range(len(coord_A)):
		## find distance of all atoms from the plane
		dist_A[atom] = sp.fabs(sp.sum(normal_A*coord_A[atom])+d) 	## normal_A is normalized
		#print(atom,dist_A[atom]) 
		## 1st layer, the defined plane
		if dist_A[atom] < small_dist: 
			atoms1_L.append(atom)
			dist_A[atom] = 1e8
	## TODO: add variable number of layers
	#for layer in range(NLayers-1): ...
	## find atoms in 2nd layer
	if NLayers > 1:
		distance = sp.amin(dist_A)
		for atom in range(len(coord_A)):
			if dist_A[atom] < distance+small_dist: 
				atoms2_L.append(atom)
				dist_A[atom] = 1e8
	## find atoms in 3rd layer
	if NLayers > 2:
		distance = sp.amin(dist_A)
		for atom in range(len(coord_A)):
			if dist_A[atom] < distance+small_dist: 
				atoms3_L.append(atom)
				dist_A[atom] = 1e8
	return atoms1_L, atoms2_L, atoms3_L


def SelfEnergyVec(AtomicSeg_A,AtomsSE1_L,AtomsSE2_L,AtomsSE3_L,ImSE1,ImSE2,ImSE3,NBF):
	''' constructs the self-energy vector (left or right) '''
	NAtom = len(AtomicSeg_A)
	SE_A  = sp.zeros(NBF,dtype='complex')
	k = 0
	for atom in range(NAtom):
		NL = AtomicSeg_A[atom]
		if   atom in AtomsSE1_L: SE_A[k:k+NL] = -1.0j*sp.ones(NL)*ImSE1
		elif atom in AtomsSE2_L: SE_A[k:k+NL] = -1.0j*sp.ones(NL)*ImSE2
		elif atom in AtomsSE3_L: SE_A[k:k+NL] = -1.0j*sp.ones(NL)*ImSE3
		k += NL
	## TODO: polarized electrodes
	## for polarized electrodes, this must be updated
	if cmos or NSpin == 2: SE_A = sp.concatenate([SE_A,SE_A])
	return SE_A


def ReadSEfile(fname,AtomicSeg_A,NBF):
	''' read self-energy from file 
	second columns is the lead, 0 = left, 1 = right 
	real part is not written to nor read from this file,
	read RSFactor tag from trans.in to construct the real part '''
	if fname not in listdir('.'):
		print(' - Error: self-energy file '+str(fname)+' missing.')
		exit(1)
	## structure of fname is [atom,lead,ImSE]
	SEdata_A = sp.loadtxt(fname)
	N = len(SEdata_A)
	## caution, atoms are numbered from one, arrays from zero
	AtomsSE_A = SEdata_A[:,0]-1
	Lead_A = SEdata_A[:,1]
	ImSE_A = SEdata_A[:,2]
	NAtom = len(AtomicSeg_A)
	SEleft_A  = sp.zeros(NBF,dtype='complex')
	SEright_A = sp.zeros(NBF,dtype='complex')
	k = 0
	for atom in range(NAtom):
		NL = AtomicSeg_A[atom]
		if atom in AtomsSE_A: 
			index = sp.where(AtomsSE_A==atom)[0][0]
			if Lead_A[index] == 0: ## left lead
				SEleft_A[k:k+NL] = -sp.ones(NL)*1.0j*ImSE_A[index]
			elif Lead_A[index] == 1: ## right lead
				SEright_A[k:k+NL] = -sp.ones(NL)*1.0j*ImSE_A[index]
		k += NL
	## TODO: polarized electrodes
	## for polarized electrodes, this must be updated
	if cmos or NSpin == 2: 
		SEleft_A  = sp.concatenate([SEleft_A ,SEleft_A ])
		SEright_A = sp.concatenate([SEright_A,SEright_A])
	return SEleft_A, SEright_A


def UpdateSelfEnergy(ReSE,SEL_A,SER_A):
	''' updates the real part of the self-energy vectors '''
	SEL_A = ReSE*sp.imag(SEL_A)+1.0j*sp.imag(SEL_A)
	SER_A = ReSE*sp.imag(SER_A)+1.0j*sp.imag(SER_A)
	return SEL_A,SER_A


def MagneticVec(AtomicSeg_A,MagAtoms_A,dh,NBF):
	''' constructs a vector of local magnetic fields '''
	NAtom = len(AtomicSeg_A)
	hloc_A  = sp.zeros(NBF)
	k = 0
	for atom in range(NAtom):
		NL = AtomicSeg_A[atom]
		if atom in MagAtoms_A: hloc_A[k:k+NL] = dh
		k += NL
	return hloc_A
	

def TransformGamma(GammaL_A,GammaR_A,Bmat_A):
	''' transforms coupling vectors Gamma to basis of eigenvectors of H(ext) 
	to calculate transmission function using diagonal Green function
	resulting Gammas are matrices, not vectors 
	Gamma matrices for complex density matrix are different! '''
	Binv_A = inv(Bmat_A)
	## this is the bottleneck of complex density matrix calculation
	if cmos: ## complex case, Gammas are ordered as {|alpha>,|beta>}
		NBF = int(len(GammaL_A)/2)
		#t = time()
		Bup_A =  Bmat_A[:NBF,:] ## spin-up block of B (top)
		Bdn_A =  Bmat_A[NBF:,:] ## spin-dn block of B (bottom)
		BIup_A = Binv_A[:,:NBF] ## spin-up block of inv(B) (left)
		BIdn_A = Binv_A[:,NBF:] ## spin-dn block of inv(B) (right)
		#GammaL2up_A = zgemm(1.0,sp.conj(Bup_A.T)*GammaL_A[:NBF],Bup_A)
		#GammaL2dn_A = zgemm(1.0,sp.conj(Bdn_A.T)*GammaL_A[NBF:],Bdn_A)
		GammaL2up_A = multi_dot([sp.conj(Bup_A.T)*GammaL_A[:NBF],Bup_A])
		GammaL2dn_A = multi_dot([sp.conj(Bdn_A.T)*GammaL_A[NBF:],Bdn_A])
		#print('time multUp',time()-t)
		#t = time()
		BICUp_A = sp.conj(Binv_A.T[:NBF,:])
		BICDn_A = sp.conj(Binv_A.T[NBF:,:])
		GammaR2up_A = multi_dot([BIup_A*GammaR_A[:NBF],BICUp_A])
		GammaR2dn_A = multi_dot([BIdn_A*GammaR_A[NBF:],BICDn_A])
		#print('time multDn',time()-t)
		## for complex case, we have four Gamma matrices, 
		## we export them as pairs L(up/dn), R(up/dn)
		GammaL2_A = sp.array([GammaL2up_A,GammaL2dn_A])
		GammaR2_A = sp.array([GammaR2up_A,GammaR2dn_A])
	elif NSpin == 1: ## closed-shell
		GammaL2_A = multi_dot([sp.conj(Bmat_A.T)*sp.diag(GammaL_A),Bmat_A])
		GammaR2_A = multi_dot([Binv_A*sp.diag(GammaR_A),sp.conj(Binv_A.T)])
	else: ## NSpin = 2
		NBF = int(len(GammaL_A)/2)
		GammaL2_A = multi_dot([sp.conj(Bmat_A.T)*sp.diag(GammaL_A),Bmat_A])
		GammaR2_A = multi_dot([Binv_A*sp.diag(GammaR_A),sp.conj(Binv_A.T)])	
		## for open shell case, we have four Gamma matrices, 
		## we export them as pairs L(up/dn), R(up/dn)
		GammaL2_A = sp.array([GammaL2_A[:NBF,:NBF],GammaL2_A[NBF:,NBF:]])
		GammaR2_A = sp.array([GammaR2_A[:NBF,:NBF],GammaR2_A[NBF:,NBF:]])
	return GammaL2_A, GammaR2_A


def TransformGammaSOCdm(GammaL_A,GammaR_A,Bmat_A):
	''' transforms coupling vectors Gamma to basis of eigenvectors of H(ext) 
	their sum is used for calculation of the complex density matrix
	Gamma matrices for complex transmission calcualtion are different! '''
	t = time()
	NBF = int(len(GammaL_A)/2)
	Binv_A = inv(Bmat_A)
	BIup_A = Binv_A[:,:NBF] ## spin-up block of inv(B) (left)
	BIdn_A = Binv_A[:,NBF:] ## spin-dn block of inv(B) (right)
	BICUp_A = sp.conj(Binv_A.T[:NBF,:])
	BICDn_A = sp.conj(Binv_A.T[NBF:,:])
	GammaL2up_A = multi_dot([BIup_A*GammaL_A[:NBF],BICUp_A])
	GammaL2dn_A = multi_dot([BIdn_A*GammaL_A[NBF:],BICDn_A])
	GammaR2up_A = multi_dot([BIup_A*GammaR_A[:NBF],BICUp_A])
	GammaR2dn_A = multi_dot([BIdn_A*GammaR_A[NBF:],BICDn_A])
	GammaL2_A = GammaL2up_A+GammaL2dn_A
	GammaR2_A = GammaR2up_A+GammaR2dn_A
	print('time multGamma',time()-t)
	return GammaL2_A, GammaR2_A

## selib.py END

