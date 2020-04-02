# Yangon
=====
#### Description:

Yangon is a code for calculating transport properties of molecular junctions. It is interfaced
with TURBOMOLE DFT package. The implemented method is described in [1].

Yangon is in very early stage of development and not yet ready for obtaining reliable results. Use on your own risk.

Yangon is a free software distributed under the GPL license.

#### Project homepage:
[github.com/pokornyv/Yangon](https://github.com/pokornyv/Yangon)

#### List of files:
- *yangon.py* - main code  
- *config.py* - reads the configutaion file *trans.in*, sets up the parameters  
- *iolib.py* - I/O functions  
- *tmlib.py* - functions for processing TURBOMOLE input  
- *selib.py* - functions for the self-energy caculations  
- *densmat.py* - functions for non-equilibrium density matrix calculations  
- *landauer.py* - functions for calculating transmission functions  
- *trans.in* -  configutaion file template  
- *LICENSE* - a copy of the GNU General Public License  
- *README.md* - this document  
- *trans.md* - configuration file description  

#### References:
1. [A. Arnold, F. Weigend, F. Evers, *J. Chem. Phys.* **126**, 174101 (2007).](https://doi.org/10.1063/1.2716664)

