################################################################
# Yangon                                                       #
# transport properties from ab-initio DFT results              #
# Copyright (C) 2019-2020 Vladislav Pokorny; pokornyv@fzu.cz   #
# homepage: github.com/pokornyv/Yangon                         #
# method described in J. Chem. Phys. 126, 174101 (2007).       #
################################################################

from config import *
from tmlib import GreenFunctionZero


def runTrans(Emin,Emax,Estep,EnExt,GammaL,GammaR,EF,ReSigma1):
	''' runs the cycle over energy points and calculates transmission function 
	Gammas are matrices and EnExt is a vector of eigvalues '''
	NPoints = int((Emax-Emin)/Estep)
	En_A = sp.linspace(Emin,Emax,NPoints)
	NP = len(En_A)
	NPP = int(NP/10)
	if cmos:  Trans_A = sp.zeros([NP,4],dtype=complex)
	else:     Trans_A = sp.zeros([NP,2],dtype=complex)
	k = 0
	for w in En_A:
		if cmos: Trans_A[k] = TransmissionSOC([w,EnExt,GammaL,GammaR])
		else:    Trans_A[k] = Transmission([w,EnExt,GammaL,GammaR])
		k += 1
		if int(k%NPP) == 0: print('{0: 2d}%'.format(int(100.0*k/NP)),end='',flush=True)
	print()
	return Trans_A, En_A


def runTransPara(Emin,Emax,Estep,EnExt,GammaL,GammaR,EF,ReSigma1):
	''' runs the cycle over energy points and calculates transmission function 
	Gammas are matrices and EnExt is a vector of eigvalues '''
	try:             NCores = int(environ['OMP_NUM_THREADS'])
	except KeyError: NCores = 1
	NPoints = int((Emax-Emin)/Estep)
	En_A = sp.linspace(Emin,Emax,NPoints)
	NP  = len(En_A)
	NPP = int(NP/10)
	if cmos:  Trans_A = sp.zeros([NP,4],dtype=complex)
	else:     Trans_A = sp.zeros([NP,2],dtype=complex)
	k = 0
	pool = Pool(NCores)
	print(NCores)
	params_L = [(w,EnExt,GammaL,GammaR) for w in En_A]
	if cmos: 
		Trans_A = pool.map(TransmissionSOC,params_L)
	else:
		Trans_A = pool.map(Transmission,params_L)
	return Trans_A, En_A


def Transmission(params_L):
	''' calculates transmission using diagonalized Green function 
	GL/R are matrices, G is diagonal vector 1/(w-h) '''
	w,EnExt,GammaL2,GammaR2 = params_L
	if NSpin == 1:
		N = len(EnExt)
		GF      = 1.0/(sp.ones(N)*(w+1.0j*izero)-EnExt)
		GFcross = sp.conj(GF)  ## GF_A is diagonal
		#T = sp.trace(sp.dot(GammaL2*GF,GammaR2*GFcross))
		AL_A = GammaL2*GF
		AR_A = GammaR2*GFcross
		Tup = Tdn = sp.sum(AL_A*AR_A.T)
	else: ## NSpin = 2:
		N2 = int(len(EnExt)/2)
		## spin-up:
		GF      = 1.0/(sp.ones(N2)*(w+1.0j*izero)-EnExt[:N2])
		GFcross = sp.conj(GF)  ## GF_A is diagonal
		AL_A = GammaL2[0]*GF
		AR_A = GammaR2[0]*GFcross
		Tup = sp.sum(AL_A*AR_A.T)
		## spin-dn:
		GF      = 1.0/(sp.ones(N2)*(w+1.0j*izero)-EnExt[N2:])
		GFcross = sp.conj(GF)  ## GF_A is diagonal
		AL_A = GammaL2[1]*GF
		AR_A = GammaR2[1]*GFcross
		Tdn = sp.sum(AL_A*AR_A.T)
	return Tup, Tdn


def TransmissionSOC(params_L):
	''' calculates transmission as Tr[GL G GR G+] 
	GL/R are arrays of two matrices, up/dn, G is a vector '''
	#t = time()
	w,EnExt,GammaL2,GammaR2 = params_L
	T = sp.zeros(4,dtype=complex)
	N = len(EnExt)
	GF      = 1.0/(sp.ones(N)*(w+1.0j*izero)-EnExt)
	GFcross = sp.conj(GF)  ## GF is a vector
	k = 0
	for s1 in [0,1]:
		for s2 in [0,1]:
			#T[k] = sp.trace(sp.dot(GammaL2[s1]*GF,GammaR2[s2]*GFcross))
			AL_A = GammaL2[s1]*GF
			AR_A = GammaR2[s2]*GFcross
			T[k] = sp.sum(AL_A*AR_A.T)
			k +=1
	#print(T)
	#print(time()-t)
	return T

## old routines used as a benchmark
"""
def TransmissionNodiag(w,Hext,GammaL,GammaR):
	''' calculates transmission as Tr[GL G GR G+] 
	GL/R are vectors, G is a matrix '''
	GF = GreenFunctionZero(w,Hext,izero)
	GFcross = sp.conj(GF.T)
	T = sp.trace(sp.dot(sp.diag(GammaL),sp.dot(GF,sp.dot(sp.diag(GammaR),GFcross))))
	return T

def TransmissionSOCnodiag(w,Hext,GammaL,GammaR):
	''' calculates transmission as Tr[GL G GR G+] 
	GL/R are vectors, G is a matrix '''
	T = sp.zeros(4,dtype=complex)
	N = int(len(GammaL)/2)
	k = 0
	GF = GreenFunctionZero(w,Hext,izero)
	for s1 in [0,1]:
		GL = GammaL[s1*N:(s1+1)*N]
		for s2 in [0,1]:
			GR = GammaR[s2*N:(s2+1)*N]
			GFblock      = GF[s1*N:(s1+1)*N,s2*N:(s2+1)*N]
			GFcrossblock = sp.conj(GF.T)[s2*N:(s2+1)*N,s1*N:(s1+1)*N]
			T[k] = sp.trace(sp.dot(sp.diag(GL),sp.dot(GFblock,sp.dot(sp.diag(GR),GFcrossblock))))
			k +=1
	return T
"""

## landauer.py END

