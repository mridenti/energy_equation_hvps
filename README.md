# energy_equation_hvps

1) Description of Mathematica notebook aer_equil.nb

In the Mathematica notebook called aer_equil.nb we present the basic theory of the energy balance equation of particle grains immersed in a plasma medium. It was run and tested in Mathematica 10.0 . 

The equations are solved using physical parameters that are consistent with the conditions found on a High Velocity Plasma Spray (HVPS) plasma torch (Journal of Materials Processing Technology, v. 237, p. 351-360, 2016). 

For more information on the physical motivation of this problem, wait for the paper "Temperature Measurement by Optical Emission Spectroscopy of the Plasma Jet Produced by a High Velociy Plasma Spray (HVPS)", to be published soon. 

2) Description of Matlab code N2_CN_rot_vib_fit_2T_VibDist.m

This code fits any given molecular spectra in the wavelength range 376 nm to 393 nm containing the following molecular system: CN violet, N2+ first negative system and N2 second positive system. The molecular state distribution model is specific for the case reported in the paper  "Temperature Measurement by Optical Emission Spectroscopy of the Plasma Jet Produced by a High Velociy Plasma Spray (HVPS)", where the vibrational states of CN were overpoplated. 

The input files absoluto.txt and resposta.txt are calibrating data related to the spectrometer response function. The file containing the experimental spectrum is named medida15mm.txt. It contains the data corresponging to one of the conditions analysed in the above mentioned paper. 

For more information about the spectra fitting procedure and the physical parameters used in the code, please refer to the following paper: 
M. Ridenti et al, Plasma Chemistry and Plasma Processing, 38(2), 311-329   
