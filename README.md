# HEDGES

A package for encoding and decoding arbitrary byte data to and from strands of DNA using a robust an error-correcting code (ECC).

### HEDGES Error-Correcting Code for DNA Storage Corrects Indels and Allows Sequence Constraints

**William H. Press, John A. Hawkins, Stephen Knox Jones Jr, Jeffrey M. Schaub, and Ilya J. Finkelstein**

*Proc Natl Acad Sci*. accepted for publication (June, 2020)

**Updated and improved by Karam Abu Hanna, Melad Shinnawi**

### Installation

The following instructions should work across platforms, except that installing virtualenv with apt-get is Ubuntu specific. For other platforms, install virtualenv appropriately if desired.

First, clone the repository to a local directory:

```
git clone https://github.com/MeladSh/HedgesProject.git
```

Optionally, you can install into a virtual environment (recommended):

```
sudo apt-get install -y virtualenv
cd hedges
virtualenv envhedges
. envhedges/bin/activate
```

Now install required packages:

```
sudo apt-get install python2.7-dev
sudo apt-get install python2.7-numpy
```

### What is supplied
Supplied is not a single program, but a kit for variable user applications.  The kit consists of

1. C++ source code that compiles (in Linux) to the Python-includable module `NRpyDNAcode`.  Precompiled binaries are supplied for Python 2.7 in Linux, but recompilation may be necessary if these don't work.  This module implements the HEDGES "inner code" as described in the paper.

2.  C++ source code that compiles (in Linux) to the Python-includable module `NRpyRS`.  Precompiled binaries are supplied for Python 2.7 in Linux, but recompilation may be necessary if these don't work.  This module implements the Schifra Reed-Solomon Error Correcting Code Library.  See http://www.schifra.com  for details and license restrictions.  This module is not needed for the HEDGES inner code, but is needed only to implement the "outer code" as described in the paper.  Some users will instead want to utilize their own outer codes.
 
3.  Python program `print_module_test_files.py`, which verifies that the above modules can be loaded and prints their usage.  Most users will not need to use any of the routines in these files directly, but should instead use the Python functions in the following file:
 
4. Python program `test_program.py` .  This defines various user-level functions for implementing the HEDGES inner and Reed-Solomon outer codes as described in the paper.  The example inputs arbitrary bytes from the file `WizardOfOzInEsperanto.txt`, encodes a specified number of packets (each with 255 DNA strands), corrupts the strands with a specified level of random substitutions, insertions, and deletions, decodes the strands, and verifies the error correction.  To better validate the installation, the code rate and corruption level set by default are chosen to be stressful to HEDGES and is greater than that in an intended use case. 

### Testing and familiarization

Run the program `test_program.py` using the command `python2.7 test_program.py`.
When running the program you can take the default substitution, deletion, insertion rates or you can set custom rates.
Also, default data input file is WizardOfOzInEsperanto.txt, it also can be changed to a custom input file.
Same for strand length, custom DNA strand length can be provided.
Same for defining output path, custom output path can be provided.

### Recompiling the C++ modules

The modules are built using the Numerical Recipes C++ class library `nr3python.h` . This is included here and also freely available for unlimited distribution at http://numerical.recipes/nr3python.h .  Generally, you will not need to understand this library, but, if you are curious, a tutorial on its use is at http://numerical.recipes/nr3_python_tutorial.html .  You should also consult this tutorial if you have difficulty recompiling the modules.  Note that while other Numerical Recipes routines are copyright and require a license, no restricted routines are used in the two modules here supplied.

In Linux, go to the directory `LinuxC++Compile` containing the source code and run the script `compile_all.sh` .  Then copy the two files produced, `NRpyDNAcode.so` and `NRpyRS.so`, to the directory containing `test_program.py`.  The most common source of errors is the compiler's inability to find required Python and Numpy include and library files that are part of your Python installation.  Unfortunately, we can't help you with that.

### LICENSE (MIT License)

Copyright 2020 by William H. Press, John A. Hawkins, Stephen K. Jones jr, and Ilya J. Finkelstein

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
