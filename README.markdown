Non-Linear Mapping
================================

Content
-------------------------

This is a simple example of the multicenter nonlinear mapping described in "Nonlinear grid mapping applied to an FDTD-based, multi-center 3D Schrödinger equation solver".


Download
-------------------------

Clone the git repository:

``git clone git://github.com/nbigaouette/nonlinearmapping.git``

You can also download an archive from a tag.


Attribution
-------------------------

You are kindly asked to cite the following publication if you use this work:

* N. Bigaouette, E. Ackad, L. Ramunno. Nonlinear grid mapping applied to an FDTD-based, multi-center 3D Schrödinger equation solver. Submitted to
Computer Physics Communications, 2011.


Requirements
-------------------------

The code is written in Python on [ArchLinux](http://www.archlinux.org/). It should work on Pyton 2.6 and 2.7.

The exact requirements are:

* Python 2.6 or 2.7.
* [Numpy](http://numpy.scipy.org/)
* [Matplotlib](http://matplotlib.sourceforge.net/)

Due to Numpy and Matplotlib not (yet) ported to Python 3, the requirement over Python is 2.7 or less.

Windows users can try [Python(x,y)](http://www.pythonxy.com/) as it includes all requirements. [PortablePython](http://www.portablepython.com/) version [2.5.4](http://www.portablepython.com/wiki/PortablePython1.1Py2.5.4) could also be used.
Note that the scripts here have not been tested with either Python(x,y) nor PortablePython.


Usage
-------------------------
./nonlinear_mapping.py


License
-------------------------

This code is distributed under the terms of the [GNU General Public License v3 (GPLv3)](http://www.gnu.org/licenses/gpl.html) and is Copyright 2010 Nicolas Bigaouette.
