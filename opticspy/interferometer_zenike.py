import numpy as __np__
import matplotlib.pyplot as __plt__
import zernike as __zernike__
import tools as __tools__
from matplotlib import cm as __cm__



def twyman_green(coefficients, lambda_1 = 632, PR = 1):
	"""
	Genertate Twyman_Green Interferogram based on zernike polynomials
	=============================================
	
	input
	----------------------------------------------
	
	Class zernike polynomials coefficients in wavenumber
	 
	see Class:opticspy.zernike.Coefficients
	
	lambda_1: wavelength in nanometer, default = 632nm
	PR: pupil radius, default = 1mm

	output
	----------------------------------------------
	Interferogram of aberration
	"""
	lambda_1 = lambda_1*(10**-9)
	coefficients = coefficients.__coefficients__
	r = __np__.linspace(-PR, PR, 400)
	x, y = __np__.meshgrid(r,r) 
	rr = __np__.sqrt(x**2 + y**2)
	OPD = 	__zernike__.__zernikecartesian__(coefficients,x,y)*2/PR

	ph = 2 * __np__.pi * OPD
	I1 = 1
	I2 = 1
	Ixy = I1 + I2 + 2 * __np__.sqrt(I1*I2) * __np__.cos(ph)
	__tools__.makecircle(Ixy, r, PR)
#======================================================
	fig = __plt__.figure(figsize=(11, 9), dpi=60)
	__plt__.imshow(-Ixy, extent=[-PR,PR,-PR,PR])
	ax = fig.gca()
	ax.set_aspect('equal', 'datalim')
	__plt__.set_cmap('Greys')
	__plt__.title('Twyman Green Interferogram',fontsize=18)
	__plt__.colorbar()
	fig.set_tight_layout(True)
	return fig

################################################################