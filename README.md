# Spectral correlation v2
this repository is dedicated to provide c++ implementation of spectral correlation and spectral coherence for cyclostationary signal analysis.
Using visual studio 2017.
The project is based on c++ console application. 
the application calculates spectral correlation/coherence using ACP, fast ACP, FAM, Antoni et al Fast SC, 
Borghesani and Antoni fast dirichlet based SC, and Aggregated FAM.

Please cite the code as:
-Aggregated FFT Accumulation Method to Compute Spectral Correlation and Coherence and its Application in Bearing Fault Diagnosis 


# Project Dependency
- FFTW library is required to compile and run the project. Depending on your buil configuration CPU type, 
The header file, dll and intermediate library are already contained in the project.
Copy the file libfftw3-3.dll from projectdir/SigProcess/fftw32 for 32-bit building or projectdir/SigProcess/fftw64
for 64-bit building and put it in your application .exe folder

Current build configuration uses x64.
Once, you extract the downloaded zip file, open the project and select x86 build setting.
The latest c++ toolset is required to build the project.


# License
Please read [https://github.com/jaafarkh/spectcorrV2/blob/add-license-1/LICENSE](https://github.com/jaafarkh/spectcorrV2/blob/master/Academic%20License.txt)
