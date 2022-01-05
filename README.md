# Standard-DSP-
MATLAB code of typical functions and algorithms used in digital signal processing

I wrote these matlab files as simple self contained references on typical tools used in DSP

the first file (makingFFTs.m)is all about the fast fourier transform written recursively. it's compared with the DFT and matlabs built in fft functions. This file also includes a convolution function built using the FFT. For these to work properly they need specific inputs. I didnt want to implement any error checking logic as I thought it would distract from the heart of these algorithms. the files should work on their own or as long as you change the inputs so that N is powers of 2 (4,8,16,32,64,128,256,512,1024 etc) and corresponds to the actual length of the inputs. with the convolution function you also should make sure both vector inputs are the same length

the next file (makingFilters.m) is all about digital filtering. I implement some simple pole and zero filters using the Z transform as a mathematical basis. This is followed by a more complex design process for a Biquad low pass filter with variable Q, cutoff frequency and gain. There are some links and derivations in the comments

the third file (ParametricSignalModeling.m) is a set of functions written to create a signal model using inderect least squares modeling (sometimes known as inverse modeling). I have class notes that I referenced in making these and havent been able to find a good online analog to explain the mathematical derivations yet but some insight is provided in the comments.

Hopefully these are helpful for any of you looking to study or implement some barebones digital signal processing!
