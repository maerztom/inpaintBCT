This C++ code implements "Fast Image Inpainting Based on Coherence Transport"
with a mex-interface to Matlab

In order to compile do:
1. start Matlab
2. setup mex by typing "mex -setup" on Matlab's command line and choose C/C++ compiler
   (you can gcc under Linux and Mac OS X, and Visual C++ under Windows)
3. compile by running "compile_mex", output file is called "inpaintBCT"
4. test it by running "inpaint_test" (in folder "inpaint_test", make sure that it finds "inpaintBCT")  
