The cpp file contains the routines for pricing and calibration by Yiran Cui (Cui et al., 2017). Only minors adjustments were made to cope with the interest rate and dividend yield term structures, and a method for cumputing delta was added. 

'heston.cpp' requires the LAPACK and LAVMAR libraries to work. Below we provide a short step-by-step tutorial for using these libraries in Visual Studio. Alternativelly, we provide a pre-built dynamic-link library (/MD) 'hestonCalibrator.dll' so the C++ code can easily be called from R. See 'masterThesis/R/source/SV/heston_wrappers.R' for the respective R wrapping functions. 

Using the C++ code in MS Visual Studio 
---
#### 1. LAPACK Library (w/o Fortana compiler) ####
Requires: Visual Studio, CMAKE. 
* Download CLAPAK from: http://icl.cs.utk.edu/lapack-for-windows/clapack/index.html and follow the instructions under 'Easy Windows Build'
* Add the libraries libf2c.lib, lapack.lib, and blas.lib to your 'external libraries directory' 
* Ad the headers blaswrap.h, clapack.h, and f2c.h to your 'external include directory' 

#### 2. LEVMAR Llibrary ####
Requires: Visual Studio, CMAKE. 
* Download last version from http://users.ics.forth.gr/~lourakis/levmar/index.html
* open CMAKE 
    * select directories 
    * manually change the name and location of the libraries according to step 1
    * configure and generate. 
* open and build LEVMAR.sln
* add levmar.lib and levmar.h to the relevant directories

#### 3. Using the c++ code ####
* Add the directory of the .h files in: Properties > C/C++ > General > Additional Include Directories
* Add the directory of the .lib files in: Properties > Linker > General > Additional Library Directories
* add the 4 libraries names in: Properties > Linker > Input > Additional Dependencies 

NB: Be consistent with 32bit vs 64bit, and Debug vs Release versions!

Using rhe dll in R 
--- 
#### Debugging #### 
* Build a Debug dll. Make sure it has debug info: General > C/C++ > Debug Information format
* Load the dll in R.
* in VS: Debug > attach to process, and select 'rsession.exe'. 
* Tools > Options > Debugger > Uncheck "Enable Just my Code", for debugging inside external libraries 
