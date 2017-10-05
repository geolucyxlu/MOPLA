======================================================================================
readme.txt for the MATLAB package

This MATLAB package is provided as a supplementary material for

An optimal scheme for numerical evaluation of Eshelby tensors and its 
implementation in a MATLAB package for simulating the motion of viscous 
ellipsoids in slow flows
Mengmeng Qu, Dazhi Jiang, Xi Lu

Please let us know about any bug or error:
       mqu5@uwo.ca
       djiang3@uwo.ca
       xlu245@uwo.ca 

May 2016
=======================================================================================
LICENSE
Permission is hereby granted, free of charge, to any person obtaining a
copy of this software and associated documentation files (the "Software"), 
to deal in the Software without restriction, including without limitation 
the rights to use, copy, modify, merge, publish, distribute, sublicense, 
and/or sell copies of the Software, and to permit persons to whom the 
Software is furnished to do so, subject to the following conditions:

    *The above copyright notice, this permission notice and the following 
     disclaimer shall be included in all copies or substantial portions of 
     the Software.
 
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE 
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, 
WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN 
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
=========================================================================================
NOTES
1) The MATLAB functions are wrriten in MATLAB version R2015a;

2) Before executing these functions, one needs to add all folders and subfolders to 
   current path in MATLAB;

3) "partTGL.mexw64" and "partTLeb.mexw64" in "Routines" are two MEX(MATLAB Executable) 
   files, written in C language and then compield using the MATLAB's built-in mex 
   function and the Microsoft Visual C++ 2015 compiler on Win64 platform using MATLAB 
   version R2016a. 
   So if your default C/C++ compiler in MATLAB and/or the platform and/or the MATLAB 
   version are different from above, the compiled files (i.e.,"partTGL.mexw64" and 
   "partTLeb.mexw64") may not run properly. You need to re-compile the C codes:
   
   Step 1. Open MATLAB, and type the following command in MATLAB command window,
   
           >> mex -setup
		   
           mex will locate the installed compilers and choose a default compiler. 
           If you do not have an avialable compiler installed on your machine, you can
           visit the Mathworks website and check the supported compilers for your MATLAB
           (google "supported compilers for MATLAB version####"). Install a compiler, and
           then re-do Step 1.
		   
   Step 2. Set the sudfolder "C codes" in "Routines" as the current folder in the MATLAB
           toolbar. You will see the two C codes (partTGL.c, partTLeb.c) in the MATLAB 
	   current folder. 
			
   Step 3. Compile the two C codes by running the following command in the command window,

           >> mex partTGL.c
           >> mex partTLeb.c

           The c-mex invokes the selected compiler to compile, link, and generate the 
           binaries, "partTGL.mex###" and "partTLeb.mex###". 	
		   
   Step 4. Copy the new mex files, "partTGL.mex###" and "partTLeb.mex###", to the folder"Routines",
           and delete the previous "partTGL.mexw64" and "partTLeb.mexw64".

   Step 5. Add all folders and subfolders to current path in MATLAB, and have fun! :)

4) You can run the two functions ("validation_deformable.m", "validation_rigid.m") in the folder 
   "Validations", and obtain the Figs.7&8 in the paper.
			
   

