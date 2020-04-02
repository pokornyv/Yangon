Main input file description
=====
Here we provide short description for the tags in the main input file *trans.in*  
The file is separated into sections according to different tasks.  
Yes/no values are coded as 0/1, value in brackets is the default.  

### Params 
- *EFermi* (no default) - Fermi energy in Hartrees. Written by yangon during calculation.
- *complex*  (0) - Complex calculation (read *spinor.r/spinor.i* instead of *mos* or *alpha/beta*).
- *spin_channels* (1) - Number of spin channels. 1 or 2 means closed shell or open shell calcualtion. Ignored in complex case.
- *NElectron* (no default) - Number of electrons. Copy it from eiger ouput.
- *TuneRSigma* (1) - Tune real part of the self-energy?
- *TuneFermi* (1) - Tune chemical potential?
- *PrintLowdin* (1) - Print lowdin charges to *lowdin.dat*?
- *bias* (0.0) - Voltage bias in *biasUnit* units
- *biasUnit* (H) - Unit for bias, H or eV
- *Pulay* (1) - Use Pulay mixing of density matrices? Used in self-consistent calculations (if *tm2ait* on)
- *FixMix* (0) - Use fixed weights in mixing. Debug only, do not use.
- *MaxMix* (5) - Number of old iterations to use in Pulay mixer.
- *dmix* (0.1) - Linear mixing parameter. Used if *Pulay* os off.

### TurboInput
- *TMoutputFile*     (proper.out) - file with basic output of Turbomole and the overlap matrix
- *TMcoordFile*      (coord) - geometry file
- *TMmosFile*        (mos) - molecular orbital file read if complex of off and *spin_channels* = 1
- *TMalpha*          (alpha) - molecular orbital file read if complex of off and *spin_channels* = 2
- *TMbeta*           (beta) - molecular orbital file read if complex of off and *spin_channels* = 2
- *TMspinorR*        (spinor.r) - molecular orbital file read if *complex* is on
- *TMspinorI*        (spinor.i) - molecular orbital file read if *complex* is on
- *SEfile*          (selfenergy.in) - file with imaginary part of self-energy for given atoms, used if *ReadSEfromFile* is on.

### SelfEnergy
- *RSFactor*        (no default) - Prefactor of the real part of the self-energy. Written by yangon during calculation. ReSigma=RSFactor*ImSigma.
- *UseSE*           (1) - Apply self-energy to boundary regions?
- *ReadSEfromFile*  (0) - Read self-energy from file *SEfile*?
- *LeftAtoms*       (no default) - Three atoms to define left boundary plane. Ignored if *ReadSEfromFile* is on.
- *RightAtoms*      (no default) - Three atoms to define right boundary plane. Ignored if *ReadSEfromFile* is on.
- *NLayers*         (2) - Number of players with self-energy. Ignored if *ReadSEfromFile* is on.
- *ImSigma1*        (0.1) - ImSigma in first outer layer. Ignored if *ReadSEfromFile* is on.
- *ImSigma2*        (0.05) - ImSigma in second outer layer. Ignored if *ReadSEfromFile* is on.
- *ImSigma3*        (0.025) - ImSigma in third outer layer. Ignored if *ReadSEfromFile* is on.

### Trans
- *CalcTrans*       (1) - Calculate transmission function?
- *Emin*            (-0.5) - Minimum of energy interval for transmission calculation, in Hartrees.
- *Emax*            ( 0.5) - Maximum of energy interval for transmission calculation, in Hartrees.
- *Estep*           (0.001) - Step in energy interval, in Hartrees.

### tm2ait
- *tm2trans*        (0) - Run self-consistent loop with Turbomole?
- *iter*            (0) - Actual iteration in self-consistent loop.

### Mag 
- *UseMag*          (0) - Use magnetic field?
- *dh*              (0) - Magnetic field, in Hartree units (times g.mu_B). One number for field in z direction, three values for arbitrary direction.
- *mag_atoms*       (no default) - List of atoms which feel the magnetic field. Use 0 to apply field to all atoms.

