# TMS_Efield_Solvers
Introduction to the code

This depository contains implementations of all 1st-3rd order FEM, 27 point stencil FDM, and 0th order BEM formulations used in “Conditions for numerically accurate TMS electric field simulation,” published in Brain Stimulation Journal. The "acc-TMS-codes-manual.pdf" document introduces the reader into the code requirements, as well as, provides some mathematical background on the different simulation methods.

System requirements

All codes require a 64-bit system and either Windows , MacOS, or Linux operating systems. All codes require a Matlab installation, and FDM and FEM require the Matlab PDE toolbox. Finally, these codes should work with most versions of Matlab. Note that these codes use OpenMP, which is not supported by matlab, as a result, the code might not work. If you experience problems, please E-mail me the operating system, error message, and matlab version at luis.gomez@duke.edu, and I will try to help you link to the correct libraries.

This code uses subroutines from the following repositories:

Matlab interface generated using:
http://www.cs.cornell.edu/~bindel/sw/mwrap/
https://github.com/ahbarnett/mwrapdemo

Fast 1/R^n evaluations using FMM from:
https://github.com/zgimbutas/fmmlib3d
