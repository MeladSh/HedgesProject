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
sudo chmod u+x setup.sh
./setup.sh
```

### What is supplied
Supplied is not a single program, but a kit for variable user applications.  The kit consists of

1. C++ source code that compiles (in Linux) to the Python-includable module `NRpyDNAcode`.  This module implements the HEDGES "inner code" as described in the paper.

2. C++ source code that compiles (in Linux) to the Python-includable module `NRpyRS`.  This module implements the Schifra Reed-Solomon Error Correcting Code Library.  See http://www.schifra.com  for details and license restrictions.  This module is not needed for the HEDGES inner code, but is needed only to implement the "outer code" as described in the paper.  Some users will instead want to utilize their own outer codes.
 
3. Python program `project.py` .  This defines various user-level functions for implementing the HEDGES inner and Reed-Solomon outer codes as described in the paper.  The default example inputs arbitrary bytes from the file `WizardOfOzInEsperanto.txt`, encodes a specified number of packets (each with 255 DNA strands), corrupts the strands with a specified level of random substitutions, insertions, and deletions, decodes the strands, and verifies the error correction.  To better validate the installation, the code rate and corruption level set by default are chosen to be stressful to HEDGES and is greater than that in an intended use case. 

### Testing and familiarization

Compile and run the program `project.py` using the running script `./run.sh`.
When running the program you can take choose the following: 
1.	Choose custom strand length, default value is 300.
2.	Choose custom output path, default output path is stdout.
3.	Choose custom code rate of the following: [0.75, 0.6, 0.5, 0.33, 0.25, 0.166], default code rate is 0.5.
4.	Choose custom substitution/deletion/insertion rates, default rates are [s:0.0357, d:0.0123, i:0.00585].
5.	Choose custom input file, default input file is WizardOfOzInEsperanto.txt.

### LICENSE (MIT License)

Copyright 2020 by William H. Press, John A. Hawkins, Stephen K. Jones jr, and Ilya J. Finkelstein

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
