################################################################
# Yangon                                                       #
# transport properties from ab-initio DFT results              #
# Copyright (C) 2019-2020 Vladislav Pokorny; pokornyv@fzu.cz   #
# homepage: github.com/pokornyv/Yangon                         #
# method described in J. Chem. Phys. 126, 174101 (2007).       #
################################################################

from config import *
from iolib import *

def AtomicSegments(Atoms_L, AtomTypes_L, BasisSize_A, NBF):
	''' find segments of the basis that belong to individual atoms '''
	NAtom = len(Atoms_L)
	Segments_A = sp.zeros(NAtom,dtype='int16')
	for atom in range(NAtom):
		Segments_A[atom] = BasisSize_A[AtomTypes_L.index(Atoms_L[atom])]
	if sp.sum(Segments_A) != NBF:
		print(' - Error: AtomicSegments: sum of segment length differs from basis size.')
		if tm2trans: SaveDensmat_tm2ait(sp.zeros([2,2]),'ait2tm')
		exit(1)
	return Segments_A


def LowdinSmatrix(smat_A):
	''' construct the Lowdin S(1/2) matrix from the overlap matrix smat
	to orthogonalize the molecular orbitals '''
	nsaos = len(smat_A)
	Small_A = sp.ones([nsaos,nsaos])*small
	eigvals_A, V_A = eigh(smat_A) ## eigensystem of the overlap matrix
	Sdiag12_A    = sp.sqrt(eigvals_A)
	Sdiag12inv_A = 1.0/sp.sqrt(eigvals_A)
	Smat12_A     = sp.dot(V_A*Sdiag12_A,V_A.T)
	Smat12inv_A  = sp.dot(V_A*Sdiag12inv_A,V_A.T)
	## test if S(1/2).S(-1/2) = I
	Zero_A = sp.dot(Smat12_A,Smat12inv_A)-sp.eye(nsaos)
	#Zero_A = sp.dot(Smat12_A,Smat12_A)-smat_A
	#Zero_A = sp.dot(Smat12inv_A,Smat12inv_A)-inv(smat_A)
	if sp.amax(sp.fabs(Zero_A))>small:
		print(" - Error, test failed: S(1/2).S(-1/2) is NOT a unit matrix!")
		if tm2trans: SaveDensmat_tm2ait(sp.zeros([2,2]),'ait2tm')
		exit(1)
	else:
		if chat: print(" - test passed: S(1/2).S(-1/2) is a unit matrix")
	return Smat12_A, Smat12inv_A


def OrthoMos(mos_A,S12_A):
	''' orthogonalize molecular orbitals using the S(1/2) matrix '''
	return sp.dot(S12_A,mos_A)


def CheckOrth(mos1_A,mos2_A):
	''' chcek if columns of mos1_A are orthogonal to columns of mos2_A
	for mos1=mos2 it checks orthogonality of a set of vectors '''
	nsaos1 = len(mos1_A)
	nsaos2 = len(mos2_A)
	if nsaos1!=nsaos2:
		print(' - dim(A1) and dim(A2) differ, cannot check orthogonality')
	else:
		orth_A = sp.zeros_like(mos1_A)
		for i,j in product(range(nsaos1),repeat=2):
			orth_A[i,j] = sp.dot(sp.conj(mos1_A[:,i]),mos2_A[:,j])
		#print('(debug) First 10x10 block of the overlap matrix:')
		#PrintMatrix(orth_A[0:10,0:10])
		Res = sp.amax(sp.fabs(sp.real(orth_A-sp.eye(nsaos1))))
		Ims = sp.amax(sp.fabs(sp.imag(orth_A-sp.eye(nsaos1))))
		if chat: print(' - largest deviation: real: {0: .4e}, imag: {1: .4e}i'.format(Res,Ims))
		if sp.allclose(orth_A,sp.eye(nsaos1),atol=1e-4):
			if chat: print(' - Test passed: vectors are orthonormal')
		elif sp.allclose(orth_A,sp.diag(sp.diag(orth_A)),atol=1e-4):
			SaveVector(sp.diag(orth_A),'norm_mos.dat')
			print(65*'#')
			print(' - Warning, test failed: vectors are orthogonal but not orthonormal!')
			print(65*'#')
		else:
			SaveMatrix(orth_A,'ortho_mos.dat')
			print(60*'#')
			print(' - Warning, test failed: vectors are not orthogonal!')
			print(60*'#')


def SortEigvals(en_A,mos_A,spinblock = False, vec = 'lines'):
	''' sort the eigenvalues and shuffle eigenvectors accordingly 
	if spinblock = True, reshuffle them to create spin-up and spin-dn blocks 
	vec marks wheter eigenvectors are lines or columns of mos_A '''
	perm_A = sp.argsort(sp.real(en_A))	## permutation that sorts energies by real part
	enSort_A = en_A[perm_A]
	if vec == 'lines': mosSort_A  = mos_A[perm_A]   ## sorts lines
	else:              mosSort_A  = mos_A[:,perm_A] ## sorts columns
	if spinblock:
		Nhalf = int(len(mosSort_A)/2)
		enSort2_A = sp.zeros_like(enSort_A)
		mosSort2_A = sp.zeros_like(mosSort_A)
		for i in range(Nhalf):
			enSort2_A[i] = enSort_A[2*i]
			enSort2_A[i+Nhalf] = enSort_A[2*i+1]
			mosSort2_A[:,i] = mosSort_A[:,2*i]
			mosSort2_A[:,i+Nhalf] = mosSort_A[:,2*i+1]
		mosSort_A = sp.copy(mosSort2_A)
		enSort_A = sp.copy(enSort2_A)
	return enSort_A, mosSort_A


def Hamiltonian(En_A,MOOrth_A,test_eig = True):
	''' construct the Hamiltonian from eigenvectors and eigenvalues 
	can be used to construct density matrix, then set test_eig to False '''
	Ham_A = sp.dot(MOOrth_A*En_A,sp.conj(MOOrth_A.T))
	if sp.allclose(sp.conj(Ham_A.T), Ham_A):
		if chat: print(" - Hamiltonian is Hermitean")
	else:
		print(60*'#')
		print(" - Warning: Hamiltonian is NOT Hermitean!")
		print(60*'#')
	if test_eig:	## test if new eigenvalues match the input ones
		En2_A = eigh(Ham_A)[0]		## eigenvalues of the new Hamiltonian
		perm_A = sp.argsort(En_A)	## permutation that sorts input eigenenergies by size
		En_A = En_A[perm_A]
		if any(sp.fabs(En_A-En2_A)>small):
			print(" - Warning: test failed:")
			print(" - Eigenvalues of the Hamiltonian differ from the input ones.")
			SaveVector(En_A-En2_A,'en_diff.dat')
			print(" - largest deviation: {0: .6e}".format(sp.amax(sp.fabs(sp.real(En_A-En2_A)))))
		else:
			if chat: print(" - Test passed, largest deviation of eigenvalues from input: {0: .6e}"\
			.format(sp.amax(sp.fabs(En_A-En2_A))))
	return Ham_A


def GreenFunctionZero(w,Ham_A,izero):
	''' calculates the retarded Green function (resolvent) G(w)=1/(w-H) '''
	nsaos = len(Ham_A)
	Ginv_A = sp.eye(nsaos)*w+1.0j*izero-Ham_A
	return inv(Ginv_A)

"""
def NormalizeMos(mos1_A,mos2_A):
	''' normalize a set of orthogonal but not orthonormal eigenvectors '''
	nsaos1 = len(mos1_A)
	nsaos2 = len(mos2_A)
	if nsaos1!=nsaos2:
		print(' - dim(A1) and dim(A2) differ, cannot check orthogonality')
	else:
		norm_A = sp.zeros(nsaos1,dtype=complex)
		for i in range(nsaos1):
			norm_A[i] = sp.dot(sp.conj(mos1_A[:,i]),mos2_A[:,i])
			mos1_A[:,i] /= 1.0/sp.sqrt(norm_A[i])
			mos2_A[:,i] /= 1.0/sp.sqrt(norm_A[i])
		SaveVector(norm_A,'norm1.dat')
		## check 
		for i in range(nsaos1):
			norm_A[i] = sp.dot(sp.conj(mos1_A[:,i]),mos2_A[:,i])
		SaveVector(norm_A,'norm2.dat')
	return mos1_A,mos2_A	
"""	

## tmlib.py END

