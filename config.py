################################################################
# Yangon                                                       #
# transport properties from ab-initio DFT results              #
# Copyright (C) 2019-2020 Vladislav Pokorny; pokornyv@fzu.cz   #
# homepage: github.com/pokornyv/Yangon                         #
# method described in J. Chem. Phys. 126, 174101 (2007).       #
################################################################

from sys import argv,exit,path,stdout
from os import path as ospath
from joblib import Parallel, delayed, parallel_backend
from multiprocessing import Pool, cpu_count
#path.insert(0,ospath.abspath('/home/pokorny/bin/numpy_mkl/lib64/python3.7/site-packages'))

import scipy as sp
import numpy as np
from scipy.linalg import eigh,eig,inv,solve,LinAlgError
from scipy.io import FortranFile
from scipy.optimize import brentq,fixed_point,root_scalar,newton,root
from numpy.linalg import multi_dot

from time import time,ctime
from os import listdir,environ,uname,rename
from itertools import product
from configparser import SafeConfigParser
from subprocess import call,Popen,PIPE

izero = 1e-8        ## small imaginary part if the complex energy
small = 1e-8        ## energy difference that we consider zero (in Hartree)
small_n = 1e-6		## charge difference that we consider zero
small_dist = 1e-1   ## distance that we consider zero (in a.u.)

hashes = 80*'#'

###########################################################
## reading config file ####################################
cfile = 'trans.in'

if cfile not in listdir('.'): 
	print('- Parameter file '+cfile+' missing. Exit.')
	exit(1)

config = SafeConfigParser()
config.read(cfile)

## default values
ef_in_input    = False
cmos           = False
NSpin          = 1
## NElectron has no default
TuneRSigma     = True
TuneFermi      = True
PrintLowdin    = True
Pulay          = True
FixMix         = False
MaxMix         = 5              ## number of matrices to store for Pulay
dmix           = 0.0
bias           = 0.0
biasUnit       = 'H'
esol           = 'newton'
TMoutputFile   = 'ridft.out'
TMcoordFile    = 'coord'
TMeigerFile    = 'eiger.out'
TMmosFile      = 'mos'
TMalpha        = 'alpha'
TMbeta         = 'beta'
TMspinorR      = 'spinor.r'
TMspinorI      = 'spinor.i'
SEfile         = 'selfenergy.in'
##LeftAtoms RightAtoms have no defaults
UseSE          =  True
ReadSEfromFile =  False
NLayers        =  1
ReSigma        =  0.0
ImSigma1       =  0.100
ImSigma2       =  0.050
ImSigma3       =  0.025
CalcTrans      =  True
Emin           = -0.5
Emax           =  0.5
Estep          =  0.01
tm2trans       =  False
tm2iter        =  0
UseMag         =  False
dh_A           =  sp.empty(0)
mag_atoms_A    =  sp.empty(0)

## Params block
if config.has_option('Params','complex'):
	cmos = bool(int(config.get('Params','complex')))
if config.has_option('Params','spin_channels'):
	NSpin = int(config.get('Params','spin_channels'))
if config.has_option('Params','NElectron'):
	NElec = int(config.get('Params','NElectron')) ## no default for this one
else:
	print(' - Missing essetial input tag "NElectron". Exit.')
	exit(1)
if config.has_option('Params','TuneRSigma'):
	TuneRSigma = bool(int(config.get('Params','TuneRSigma')))
if config.has_option('Params','TuneFermi'):
	TuneFermi = bool(int(config.get('Params','TuneFermi')))
if config.has_option('Params','PrintLowdin'):
	PrintLowdin = bool(int(config.get('Params','PrintLowdin')))
if config.has_option('Params','Pulay'):
	Pulay = bool(int(config.get('Params','Pulay')))
if config.has_option('Params','FixMix'):
	FixMix = bool(int(config.get('Params','FixMix')))
if config.has_option('Params','MaxMix'):
	MaxMix = int(config.get('Params','MaxMix'))
if config.has_option('Params','dmix'):
	dmix = float(config.get('Params','dmix'))
if config.has_option('Params','bias'):
	bias = float(config.get('Params','bias'))
if config.has_option('Params','biasUnit'):
	biasUnit = config.get('Params','biasUnit').lower()
if config.has_option('Params','EFermi'):
	EFermi0 = float(config.get('Params','EFermi'))
	ef_in_input = True
if config.has_option('Params','solver'):
	esol = config.get('Params','solver').lower()

## TurboInput block
if config.has_option('TurboInput','TMoutputFile'):
	TMoutputFile = config.get('TurboInput','TMoutputFile')
if config.has_option('TurboInput','TMcoordFile'):
	TMcoordFile = config.get('TurboInput','TMcoordFile')
if config.has_option('TurboInput','TMeigerFile'):
	TMeigerFile = config.get('TurboInput','TMeigerFile')
if config.has_option('TurboInput','TMmosFile'):
	TMmosFile = config.get('TurboInput','TMmosFile')
if config.has_option('TurboInput','TMalpha'):
	TMalpha = config.get('TurboInput','TMalpha')
if config.has_option('TurboInput','TMbeta'):
	TMbeta = config.get('TurboInput','TMbeta')
if config.has_option('TurboInput','TMspinorR'):
	TMspinorR = config.get('TurboInput','TMspinorR')
if config.has_option('TurboInput','TMspinorI'):
	TMspinorI = config.get('TurboInput','TMspinorI')
if config.has_option('TurboInput','SEfile'):
	SEfile = config.get('TurboInput','SEfile')

## SelfEnergy block
if config.has_option('SelfEnergy','UseSE'):
	UseSE      = bool(int(config.get('SelfEnergy','UseSE')))
if UseSE:
	if config.has_option('SelfEnergy','RSFactor'):
		ReSigma = float(config.get('SelfEnergy','RSFactor'))
	if config.has_option('SelfEnergy','ReadSEfromFile'):
		ReadSEfromFile = bool(int(config.get('SelfEnergy','ReadSEfromFile')))
	if not ReadSEfromFile: ## building self-energy layer by layer
		if config.has_option('SelfEnergy','LeftAtoms'):
			[LeftA1, LeftA2, LeftA3 ] = \
			[int(i) for i in config.get('SelfEnergy','LeftAtoms' ).split()]
		else:
			print(' - Missing essetial input tag "LeftAtoms". Exit.')
			exit(1)	
		if config.has_option('SelfEnergy','RightAtoms'):
			[RightA1,RightA2,RightA3] = \
			[int(i) for i in config.get('SelfEnergy','RightAtoms').split()]
		else:
			print(' - Missing essetial input tag "RightAtoms". Exit.')
			exit(1)	
		if config.has_option('SelfEnergy','NLayers'):
			NLayers = int(config.get('SelfEnergy','NLayers'))
		if config.has_option('SelfEnergy','ImSigma1'):
			ImSigma1 = float(config.get('SelfEnergy','ImSigma1'))
		if config.has_option('SelfEnergy','ImSigma1'):
			ImSigma2 = float(config.get('SelfEnergy','ImSigma2'))
		if config.has_option('SelfEnergy','ImSigma1'):
			ImSigma3 = float(config.get('SelfEnergy','ImSigma3'))

## Trans block
if config.has_option('Trans','CalcTrans'):
	CalcTrans = bool(int(config.get('Trans','CalcTrans')))
if config.has_option('Trans','Emin'):
	Emin = float(config.get('Trans','Emin'))
if config.has_option('Trans','Emax'):
	Emax = float(config.get('Trans','Emax'))
if config.has_option('Trans','Estep'):
	Estep = float(config.get('Trans','Estep'))

## Tmait block
if config.has_option('tm2ait','tm2trans'):
	tm2trans = bool(int(config.get('tm2ait','tm2trans')))
if config.has_option('tm2ait','iter'):
	tm2iter = int(config.get('tm2ait','iter'))

## Mag block
if config.has_option('Mag','UseMag'):
	UseMag = bool(int(config.get('Mag','UseMag')))
if config.has_option('Mag','dh'):
	dh_L = config.get('Mag','dh').split()
	dh_A = sp.array([float(i) for i in dh_L])
if config.has_option('Mag','mag_atoms'):
	MagAtoms_L = config.get('Mag','mag_atoms').split()
	MagAtoms_A = sp.array([int(i) for i in MagAtoms_L])

###########################################################
## global variables #######################################
if cmos: blocks_T = ['uu','ud','du','dd']

if biasUnit == 'ev': bias /= 27.211399

## print info to standard output?
chat = True
if tm2trans: chat = False 
if tm2iter < 2 or CalcTrans: chat = True 

## mixing of density matrices, fixed numbers for MaxMix
if FixMix: 
	Pulay = True
	MaxMix = 4

## list of elements
elements_L = [
'H','He','Li','Be','B','C','N','O','F','Ne',\
'Na','Mg','Al','Si','P','S','Cl','Ar',\
'K','Ca','Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr',\
'Rb','Sr','Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te','I','Xe',\
'Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu',\
'Hf','Ta','W','Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn','Fr','Ra',\
'Ac','Th','Pa','U','Np','Pu','Am','Cm','Bk','Cf','Es','Fm','Md','No','Lr',\
'Rf','Db','Sg','Bh','Hs','Mt','Ds','Rg','Cn','Nh','Fl','Mc','Lv','Ts','Og']

unit_bohr    = 1.889725989 ## Bohr radius in Angstroms
unit_hartree = 27.211399   ## Hartree in eV

## config.py END

