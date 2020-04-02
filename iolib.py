################################################################
# Yangon                                                       #
# transport properties from ab-initio DFT results              #
# Copyright (C) 2019-2020 Vladislav Pokorny; pokornyv@fzu.cz   #
# homepage: github.com/pokornyv/Yangon                         #
# method described in J. Chem. Phys. 126, 174101 (2007).       #
################################################################

from config import *

def MakeMatrix(sin_A):
	''' make a symmetric matrix from an array representing lower triangle '''
	ls = len(sin_A)
	ln = int(0.5*(-1+sp.sqrt(1+8*ls)))	## ls  = ln*(ln+1)/2
	sout_A = sp.zeros([ln,ln])
	il = sp.tril_indices(ln) ## indices of a lower triangular matrix
	for i in range(ls):
		sout_A[il[0][i],il[1][i]] = sin_A[i]
	sout_A = sout_A+sout_A.T-sp.eye(ln)*sp.diag(sout_A)
	return sout_A

## reading Turbomole output files

def ReadAtomTypes(fname):
	''' find atom types and size of the basis for each atom
	from Turbomole output file '''
	textf='type   atoms  prim   cont   basis'
	Natom = 0
	Atoms_L = []	     ## atomic symbols, e.g. C, H, Au....
	AtomsNumber_L = []	## number of atoms of a given element
	BasisSize_L = []    ## basis size (number of CGTO) of a given element
	if chat: print(' - reading '+fname+' file...',end='',flush=True)
	if fname not in listdir('.'):
		print(' - Error: ridft output file '+str(fname)+' missing.')
		exit(1)
	f = open(fname,'r')
	while True:
		line = f.readline()
		if len(line) == 0: break	## EOF
		if line.strip() == textf:
			if chat: print(' atomic species and basis data found. ')
			f.readline() ## reading the line with ----------------------....	
			while True:
				line = f.readline()
				if line.strip()[0] == '-': 
					end = True
					break
				else: 
					line2 = line.split()
					Natom += 1
					Atoms_L.append(line2[0])
					AtomsNumber_L.append(int(line2[1]))
					BasisSize_L.append(int(line2[3]))
			if end: break	## TM output contains the basis data twice
	f.close()
	return Atoms_L, sp.array(AtomsNumber_L), sp.array(BasisSize_L)


def ReadCoordTM(fname,NAtom):
	''' read coord file with Turbomole geometry data '''
	f = open(fname,'r')
	if fname not in listdir('.'):
		print(' - Error: geometry file '+str(fname)+' missing.')
		exit(1)
	if chat: print(' - reading '+fname+' file...',end='',flush=True)
	Coord_A = sp.zeros([NAtom,3])
	Atoms_L = []
	line = f.readline()
	if line.strip() != '$coord':
		print(' - Error: '+fname+' is not a Turbomole coord file.')
		exit(1)
	atom = 0
	while True:
		line = f.readline()
		if len(line) == 0: break	## EOF, just in case
		if line.strip()[0] == '$': break
		line_L = line.split()
		Coord_A[atom,:] = [float(line_L[0]),float(line_L[1]),float(line_L[2])]
		Atoms_L.append(line_L[3])
		atom += 1
	if atom != NAtom: 
		print(' - Error, we missed some atoms from coord.')
		exit(1)
	if chat: print(' done.')
	f.close()
	return Coord_A, Atoms_L


def ReadSmatTM(fname,NBF):
	''' reads the overlap matrix from the Turbomole output file
	requires $intsdebug sao tag in control file '''
	textf = 'OVERLAP(SAO)'
	nelem = int(NBF*(NBF+1)/2) ## number of elements in
	nlines = int(nelem/3+nelem%3)	 ## number of lines with coefficients
	smat_array_L = []
	if fname not in listdir('.'):
		print(' - Error: ridft output file '+str(fname)+' missing.')
		exit(1)
	with open(fname,'r') as f:
		if chat: print(' - reading '+fname+' file...',end='',flush=True)
		t = time()
		s = 0
		while True:
			line = f.readline()
			if len(line) == 0: break	## EOF
			if line.strip() == textf:
				if chat: print(' overlap matrix found... ',end='',flush=True)
				f.readline() ## reading the line with ----------------------
				s = 1
				for i in range(nlines):
					line_L = f.readline().split()
					smat_array_L.extend(line_L)
				break
	if not s:
		print(' - Error: Overlap data not found. Put "$intsdebug sao" in control file and run ridft again.')
		exit(1)
	if chat: print(' done in {0: .2f}s.'.format(time()-t))
	smat_array_A = sp.array(smat_array_L,dtype=float)
	smat_A = MakeMatrix(smat_array_A)
	if cmos: ## we must double the overlap matrix
		smat_A = sp.block([[smat_A,sp.zeros([NBF,NBF])],[sp.zeros([NBF,NBF]),smat_A]])
	return smat_A


def ReadMosTM(fname,nsaos, sort = False):
	''' reads the mos/alpha/beta/spinor.r/spinor.i file in Turbomole format 
	!!! columns of mos_A are the eigenvectors !!!
	spinblock sorts the vectors to create spin-up and spin-down blocks
	this is a new implementation, the old one using fortranformat was super slow '''
	if fname not in listdir('.'):
		print(' - Error: molecular orbital file '+str(fname)+' missing.')
		exit(1)
	en_A  = sp.zeros(nsaos)
	mos_A = sp.zeros([nsaos,nsaos])
	#if nsaos%4 != 0: addline = 1
	nlines = int(nsaos/4)+int(bool(nsaos%4))	## number of lines with coefficients
	if chat: print(' reading '+fname+' file... ',end='',flush=True)
	t = time()
	with open(fname,'r') as f:
		content = f.readlines()
	if chat: print(' file loaded in {0: .2f}s... '.format(time()-t),end='',flush=True)
	n = 0
	t = time()
	for line in content:
		## Turbomole sometimes writes strage stuff in MOS file during calculation:
		if line[:9] == '$restartd': break
		## modified ridft code writes null characters to spinor files, \00\00\00\00
		## one can also use sed -i 's/\x00//g' spinor.r and sed -i 's/\x00//g' spinor.i
		line = line.replace('\x00','')
		## ignore header and comments
		if line[0] == '#': continue
		if line[0] == '$': continue
		if 'eigenvalue' in line:
			## if the error ValueError: could not convert string to float:
			## apperas, it means TM wrote NaN to mos files.
			en_A[n] = float(line.split()[2][11:].replace('D','E'))
			k = 0
			coeff_A = sp.zeros(nsaos)
		else:
			ncoeff = int(len(line[:-1])/20)	## last symbol is '\n'
			for i in range(ncoeff):
				coeff_A[k] = float(line[i*20:(i+1)*20].replace('D','E'))
				k += 1
			## reading coefficients, python list is faster than numpy array
			if k >= nsaos:
				mos_A[n] = coeff_A[coeff_A != sp.array(None)]
				if n == nsaos: break
				n += 1
	## sort the eigenvalues
	if sort:
		perm_A = sp.argsort(en_A)	## permutation matrix to sort energies
		en_A   = en_A[perm_A]
		mos_A  = mos_A[perm_A]
	'''
	if spinblock: ## create a spin-block structure [|alpha>,|beta>], 
		## useless, kept for historical reasons
		## eigenvectors are now lines of mos_A
		print(' - Molecular orbitals will be ordered into alpha and beta')
		mos2_A = sp.zeros_like(mos_A)
		en2_A = sp.zeros_like(en_A)
		Nhalf = int(nsaos/2)
		for i in range(Nhalf):
			mos2_A[i] = mos_A[2*i]
			en2_A[i] = en_A[2*i]
			mos2_A[i+Nhalf] = mos_A[2*i+1]
			en2_A[i+Nhalf] = en_A[2*i+1]
		mos_A = sp.copy(mos2_A)
		en_A = sp.copy(en2_A)
	print(en_A)
	'''
	if chat: print(' file processed in {0: .2f}s.'.format(time()-t),flush=True)
	return [en_A,mos_A.T]


def ReadFileTM(fname,nsaos):
	''' read lower triangle of a square matrix from file '''
	if fname not in listdir('.'):
		print(' - Error: input file '+str(fname)+' missing.')
		exit(1)
	nelem = int(nsaos*(nsaos+1)/2) ## number of elements in lower triangle
	mat_array_L = []
	if chat: print(' reading '+fname+' file...')
	t = time()
	for line in open(fname,'r'):
		if line[0] == '#': continue
		mat_array_L.extend(line.split())
	if len(mat_array_L)>nelem:
		print(" - Warning, input array too long.")
	## arrays from TM output are sometimes too long. Why???
	mat_array_A = sp.array(mat_array_L[:nelem],dtype=float)
	mat_A = MakeMatrix(mat_array_A)
	if sp.allclose(sp.conj(mat_A.T), mat_A):
		print(" - matrix is Hermitean")
	return mat_A


def ReadDMcomplexTM(NBF):
	''' read the reaa, reab, rebb, imaa, imab, imbb files '''
	Files = 1
	for fname in ['reaa', 'reab', 'rebb', 'imaa', 'imab', 'imbb']:
		if fname not in listdir('.'):
			print(' - Warning: input file '+str(fname)+' missing.')
			print(' - Turbomole density matrix will not be read')
			Files = 0
			break
	if Files:
		with FortranFile('reaa','r') as f:
			ReAA = sp.array(f.read_reals(dtype='float64')).reshape([NBF,NBF])
		with FortranFile('imaa','r') as f:
			ImAA = sp.array(f.read_reals(dtype='float64')).reshape([NBF,NBF])
		with FortranFile('reab','r') as f:
			ReAB = sp.array(f.read_reals(dtype='float64')).reshape([NBF,NBF])
		with FortranFile('imab','r') as f:
			ImAB = sp.array(f.read_reals(dtype='float64')).reshape([NBF,NBF])
		with FortranFile('rebb','r') as f:
			ReBB = sp.array(f.read_reals(dtype='float64')).reshape([NBF,NBF])
		with FortranFile('imbb','r') as f:
			ImBB = sp.array(f.read_reals(dtype='float64')).reshape([NBF,NBF])
		Dmat_A = sp.block([[ReAA+1.0j*ImAA,ReAB+1.0j*ImAB],[ReAB-1.0j*ImAB,ReBB+1.0j*ImBB]])
	else: Dmat_A = sp.zeros([2*NBF,2*NBF])
	return Dmat_A


def ReadEigerTM(fname,nsaos):
	''' read the eiger output to extract occupation numbers '''
	if fname not in listdir('.'):
		print(' - Error: eiger output file '+str(fname)+' missing.')
		exit(1)
	if NSpin == 1:
		occn_A = sp.zeros(nsaos)
	else:
		occn_A = sp.zeros([nsaos,2])
	if chat: print(' - reading '+fname+' file')
	reading = False
	k1 = k2 = 0
	for line in open(fname,'r'):
		if len(line.split()) < 1 and not reading: continue
		if len(line.split()) < 1 and reading: break
		if line.split()[0] == 'Nr.': 
			reading = True
			continue
		if reading:
			if NSpin == 1:
				## orbitals are numbered from one, arrays from zero
				orb = int(float(line.split()[0]))-1 
				if len(line.split()) == 8: occn_A[orb] = 0.0
				else: occn_A[orb] = float(line.split()[3])
			else:
			## line is longer by one column in open shell case
				orb = int(line.split()[2])-1 ## orbitals are numbered from one, arrays from zero
				if len(line.split()) == 9: 
					if line.split()[1] == 'a':	occn_A[orb][0] = 0.0
					else: 					occn_A[orb][1] = 0.0
				elif line.split()[1] == 'a': occn_A[orb][0] = float(line.split()[4])
				elif line.split()[1] == 'b': occn_A[orb][1] = float(line.split()[4])
				else: print('Error reading '+fname+' file.')
#	if k1 < nsaos: ## eiger file does not contain all virtual orbitals
#		occn_A = sp.concatenate([sp.zeros(nsaos-k1),occn_A[:k1]])
	return occn_A,int(sp.sum(occn_A))

## output functions

def WriteEiger(fname):
	''' use Turbomole's eiger to write eiger file '''
	with open(fname,'w') as f: p = Popen(['eiger','-a'],stdout=f)
	## we need this line otherwise f remains open for writing
	out = p.communicate()[0]


def PrintMatrix(M_A):
	''' prints matrix in a reasonable manner '''
	M,N = M_A.shape
	print()
	for i in range(M):
		for j in range(N):
			print('  {0: .5f}'.format(sp.real(M_A[i,j])),end='')
		print()


def SaveVector(V_A,fname):
	''' prints vector to a file '''
	f = open(fname,'w')
	f.write('# File written on '+ctime()+'\n')
	M = len(V_A)
	with open(fname,'w') as f:
		for i in range(M):
			f.write('{0: 3d}\t{1: .12f}\t{2: .12f}\n'\
			.format(i+1,sp.real(V_A[i]),sp.imag(V_A[i])))
	print(' file '+str(fname)+' saved.')


def SaveMatrix(M_A,fname):
	''' prints matrix to a file '''
	f = open(fname,'w')
	f.write('# File written on '+ctime()+'\n')
	M,N = M_A.shape
	with open(fname,'w') as f:
		for i in range(M):
			for j in range(N):
				f.write('{0: 3d}\t{1: 3d}\t{2: .12f}\t{3: .12f}\n'\
				.format(i+1,j+1,sp.real(M_A[i,j]),sp.imag(M_A[i,j])))
			f.write('\n')
	print(' file '+str(fname)+' saved.')


def PrintLowdinCharges(fname,Atoms_L,Lcharges_A,Lspins_A,Coord_A,tag = ''):
	print('\nLowdin population analysis:')
	NAtom = len(Coord_A)
	with open(fname,'w') as f:
		for atom in range(NAtom):
			line = tag+'  {0: 3d} {1:<2}    N= {2:7.4f}   NS= {3:7.4f} z= {4: .6f}'\
			.format(atom+1,Atoms_L[atom].capitalize(),Lcharges_A[atom],Lspins_A[atom],Coord_A[atom][2])
			print(line)
			f.write(line+'\n')
	print('\ntotal:     N = {0: .4f}\tNS = {1: .4f}'\
	.format(sp.sum(Lcharges_A),sp.sum(Lspins_A)))


def SaveDensmat_tm2ait(Dmat_A,fname):
	''' prints the density matrix to a file ait2tm
	for charge self-consistent calculation with Turbomole (tm2ait) '''
	t = time()
	N  = len(Dmat_A)
	N2 = int(N/2)
	if cmos: ## dmat has a block structure!
		DmatUU_A = Dmat_A[:N2,:N2]
		DmatUD_A = Dmat_A[N2:,:N2]
		DmatDD_A = Dmat_A[N2:,N2:]
		with open(fname,'w') as f:
			for i,j in product(range(int(N/2)),repeat=2):
				f.write('{0: 20.14e}\n'.format(float(sp.real(DmatUU_A[i][j]))).replace('e','D'))
			for i,j in product(range(int(N/2)),repeat=2):
				f.write('{0: 20.14e}\n'.format(float(sp.imag(DmatUU_A[i][j]))).replace('e','D'))
			for i,j in product(range(int(N/2)),repeat=2):
				f.write('{0: 20.14e}\n'.format(float(sp.real(DmatUD_A[i][j]))).replace('e','D'))
			for i,j in product(range(int(N/2)),repeat=2):
				f.write('{0: 20.14e}\n'.format(float(sp.imag(DmatUD_A[i][j]))).replace('e','D'))
			for i,j in product(range(int(N/2)),repeat=2):
				f.write('{0: 20.14e}\n'.format(float(sp.real(DmatDD_A[i][j]))).replace('e','D'))
			for i,j in product(range(int(N/2)),repeat=2):
				f.write('{0: 20.14e}\n'.format(float(sp.imag(DmatDD_A[i][j]))).replace('e','D'))
	elif NSpin == 1: ## closed shell, only lower triangle of Dmat is written
		with open(fname,'w') as f:
			for i in range(N):
				for j in range(i+1):
					f.write('{0: 20.14e}\n'.format(float(sp.real(Dmat_A[i][j]))).replace('e','D'))
	else: ## NSpin = 2, open shell, lower triangles of Dmat_uu and Dmat_dd are written 
		with open(fname,'w') as f:
			for i in range(N2): ## spin-up
				for j in range(i+1):
					f.write('{0: 20.14e}\n'.format(float(sp.real(Dmat_A[i][j]))).replace('e','D'))
			for i in range(N2): ## spin-dn
				for j in range(i+1):
					f.write('{0: 20.14e}\n'.format(float(sp.real(Dmat_A[i+N2][j+N2]))).replace('e','D'))
	print(' - file '+fname+' written in {0: .3f}s'.format(time()-t))


def LoadDensmat_tm2ait(NBF,fname):
	''' reads the density matrix from a file ait2tm to mix with the new one
	for charge self-consistent calculation with Turbomole (tm2ait) '''
	t = time()
	if cmos: Dmat_A = sp.zeros([2*NBF,2*NBF],dtype='complex')
	else:    Dmat_A = sp.zeros([NBF,NBF])
	print(' - reading '+fname+' file...',end='',flush=True)
	if fname not in listdir('.'):
		print(' - '+fname+' file is missing, we return a zero matrix')
	else:
		try:
			if cmos: ## dmat has a block structure!
				DmatUU_A = sp.zeros([NBF,NBF],dtype='complex')
				DmatUD_A = sp.zeros([NBF,NBF],dtype='complex')
				DmatDD_A = sp.zeros([NBF,NBF],dtype='complex')
				with open(fname,'r') as f:
					for i,j in product(range(NBF),repeat=2):
						DmatUU_A[i][j] = float(f.readline().replace('D','e'))
					for i,j in product(range(NBF),repeat=2):
						DmatUU_A[i][j] += 1.0j*float(f.readline().replace('D','e'))
					for i,j in product(range(NBF),repeat=2):
						DmatUD_A[i][j] = float(f.readline().replace('D','e'))
					for i,j in product(range(NBF),repeat=2):
						DmatUD_A[i][j] += 1.0j*float(f.readline().replace('D','e'))
					for i,j in product(range(NBF),repeat=2):
						DmatDD_A[i][j] = float(f.readline().replace('D','e'))
					for i,j in product(range(NBF),repeat=2):
						DmatDD_A[i][j] += 1.0j*float(f.readline().replace('D','e'))
				Dmat_A = sp.block([[DmatUU_A,DmatUD_A],[sp.conj(DmatUD_A),DmatDD_A]])
			elif NSpin == 1: ## closed shell, only lower triangle is written
				Dmat_A = sp.zeros([NBF,NBF])
				with open(fname,'r') as f:
					for i in range(NBF):
						for j in range(i+1):
							elem = float(f.readline().replace('D','e'))
							Dmat_A[i][j] = elem
							Dmat_A[j][i] = elem
			else: ## NSpin = 2, open shell, lower triangles of Dmat_uu and Dmat_dd are written 
				Dmat_A = sp.zeros([2*NBF,2*NBF])
				with open(fname,'r') as f:
					for i in range(NBF):
						for j in range(i+1):
							elem = float(f.readline().replace('D','e'))
							Dmat_A[i][j] = elem
							Dmat_A[j][i] = elem
					for i in range(NBF):
						for j in range(i+1):
							elem = float(f.readline().replace('D','e'))
							Dmat_A[i+NBF][j+NBF] = elem
							Dmat_A[j+NBF][i+NBF] = elem
		except ValueError:
			print(' ... file '+fname+' seems corrupted or empty, we return zero matrix.')
	print(' done in {0: .3f}s'.format(time()-t))
	return Dmat_A


def WriteTransmission(En_A,T_A,EF,bias,ReSigma1,ImSigma1):
	''' writes files with the transmission function 
	we write two files, one in Hartrees and other in electronvolts, energy is w.r.t. Fermi energy '''
	NP = len(En_A)
	#NPP = int(NP/10)
	fnameH  = 'transmission_H.dat'
	fnameeV = 'transmission_eV.dat'
	with open(fnameH,'w') as f1:
		with open(fnameeV,'w') as f2:
			f1.write('# File written on '+ctime()+'\n')
			f2.write('# File written on '+ctime()+'\n')
			f1.write('# EFermi = {0: .6f} H, ReSigma1 = {1: .6f} H, ImSigma1 = {2: .6f} H, bias = {3: .6f} H\n'\
			.format(EF,ReSigma1,ImSigma1,bias))
			f2.write('# EFermi = {0: .6f} eV, ReSigma1 = {1: .6f} eV, ImSigma1 = {2: .6f} eV, bias = {3: .6f} eV\n'\
			.format(EF*unit_hartree,ReSigma1*unit_hartree,ImSigma1*unit_hartree,bias*unit_hartree))
			if cmos:	
				f1.write('# E-EF [H]\t\tTuu(E)\t\tTud(E)\t\tTdu(E)\t\tTdd(E)\n')
				f2.write('# E-EF [eV]\t\tTuu(E)\t\tTud(E)\t\tTdu(E)\t\tTdd(E)\n')
			else:	
				f1.write('# E-EF [H]\t\tTup(E)\t\tTdn(E)\n')
				f2.write('# E-EF [eV]\t\tTup(E)\t\tTdn(E)\n')
			for i in range(NP):
				if cmos: 
					f1.write('{0: .6f}\t{1: .12f}\t{2: .12f}\t{3: .12f}\t{4: .12f}\n'\
					.format(En_A[i]-EF,sp.real(T_A[i][0]),sp.real(T_A[i][1]),sp.real(T_A[i][2]),sp.real(T_A[i][3])))
					f2.write('{0: .6f}\t{1: .12f}\t{2: .12f}\t{3: .12f}\t{4: .12f}\n'\
					.format((En_A[i]-EF)*unit_hartree,sp.real(T_A[i][0]),sp.real(T_A[i][1]),sp.real(T_A[i][2]),sp.real(T_A[i][3])))
				elif NSpin == 1:
					f1.write('{0: .12f}\t{1: .12f}\t{2: .12f}\n'\
					.format(En_A[i]-EF,sp.real(T_A[i][0]),sp.real(T_A[i][0])))
					f2.write('{0: .12f}\t{1: .12f}\t{2: .12f}\n'\
					.format((En_A[i]-EF)*unit_hartree,sp.real(T_A[i][0]),sp.real(T_A[i][0])))
				else:
					f1.write('{0: .6f}\t{1: .12f}\t{2: .12f}\n'\
					.format(En_A[i]-EF,sp.real(T_A[i][0]),sp.real(T_A[i][1])))
					f2.write('{0: .6f}\t{1: .12f}\t{2: .12f}\n'\
					.format((En_A[i]-EF)*unit_hartree,sp.real(T_A[i][0]),sp.real(T_A[i][1])))
				#if int(i%NPP) == 0: print('{0: 2d}%'.format(int(100.0*i/NP)),end='',flush=True)
	#print()
	print('- File '+str(fnameH)+' saved.')
	print('- File '+str(fnameeV)+' saved.')


def UpdateInfile(tag,group,value,form='float'):
	''' updates the value of a tag in main input file cfile (trans.in) '''
	cgroup = '['+group+']'
	with open(cfile,'r') as f:  lines_L = f.readlines()
	if form == 'float':	
		wl = tag+(16-len(tag))*' '+': {0: .10f}'.format(float(sp.real(value)))+'\n'	
	else: ## integer
		wl = tag+(16-len(tag))*' '+': {0: 5d}'.format(value)+'\n'	
	## if the tag is already in the input file, is updated
	line2_L = [wl if len(line) > 1 and line.split()[0] == tag else line for line in lines_L]
	## if the tag was originally missing from input file, it is added
	if not wl in line2_L: 
		for line in lines_L:
			if len(line) > 1  and line.split()[0] == cgroup:
				wl = cgroup+'\n'+wl
		line2_L = [wl if len(line) > 1 and line.split()[0] == cgroup else line for line in lines_L]
	if not wl in line2_L: ## just in case...
		print(' - Warning, tag '+tag+' was not written to '+cfile)
	print(' - tag '+tag+' in '+cfile +' is updated')
	with open(cfile,'w') as f:  f.writelines(line2_L)

## iolib.py END

