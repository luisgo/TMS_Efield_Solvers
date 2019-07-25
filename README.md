# TMS_Efield_Solvers
Introduction to the code

This depository contains implementations of all FEM, FDM, and BEM formulations used in “Conditions for numerically accurate TMS electric field simulation,” published in Brain Stimulation Journal. The 'accTMS-codes-manual.pdf' document introduces the reader into the code requirements, as well as, provides some mathematical background on the different simulation methods.

System requirements

Codes require the use of a 64-bit system and Windows, MacOS, or Linux operating systems. All codes require a Matlab installation, and FDM and FEM require the Matlab PDE toolbox. Finally, these codes should work with most versions of Matlab. Note these codes use OMP, which currently is not supported by matlab, as such if you have problems I might be able to help with the linking process. If you experience any problems please let me know the operating system and Matlab version, so that I can recommend how to proceed. 