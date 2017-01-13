import numpy as __np__
import matplotlib.pyplot as __plt__
#from opticspy.zernike import __zernikecartesian__
from opticspy import tools as __tools__
from matplotlib import cm as __cm__
from numpy import cos as __cos__
from numpy import sin as __sin__
from numpy import sqrt as __sqrt__


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
	OPD = 	__zernikecartesian__(coefficients,x,y)*2/PR

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


def __zernikecartesian__(coefficient,x,y):
	"""
	------------------------------------------------
	__zernikecartesian__(coefficient,x,y):

	Return combined aberration

	Zernike Polynomials Caculation in Cartesian coordinates

	coefficient: Zernike Polynomials Coefficient from input
	x: x in Cartesian coordinates
	y: y in Cartesian coordinates
	------------------------------------------------
	"""
	Z = [0]+coefficient
	r = __sqrt__(x**2 + y**2)
	Z1  =  Z[1]  * 1
	Z2  =  Z[2]  * 2*x
	Z3  =  Z[3]  * 2*y
	Z4  =  Z[4]  * __sqrt__(3)*(2*r**2-1)
	Z5  =  Z[5]  * 2*__sqrt__(6)*x*y
	Z6  =  Z[6]  * __sqrt__(6)*(x**2-y**2)
	Z7  =  Z[7]  * __sqrt__(8)*y*(3*r**2-2)
	Z8  =  Z[8]  * __sqrt__(8)*x*(3*r**2-2)
	Z9  =  Z[9]  * __sqrt__(8)*y*(3*x**2-y**2)
	Z10 =  Z[10] * __sqrt__(8)*x*(x**2-3*y**2)
	Z11 =  Z[11] * __sqrt__(5)*(6*r**4-6*r**2+1)
	Z12 =  Z[12] * __sqrt__(10)*(x**2-y**2)*(4*r**2-3)
	Z13 =  Z[13] * 2*__sqrt__(10)*x*y*(4*r**2-3)
	Z14 =  Z[14] * __sqrt__(10)*(r**4-8*x**2*y**2)
	Z15 =  Z[15] * 4*__sqrt__(10)*x*y*(x**2-y**2)
	Z16 =  Z[16] * __sqrt__(12)*x*(10*r**4-12*r**2+3)
	Z17 =  Z[17] * __sqrt__(12)*y*(10*r**4-12*r**2+3)
	Z18 =  Z[18] * __sqrt__(12)*x*(x**2-3*y**2)*(5*r**2-4)
	Z19 =  Z[19] * __sqrt__(12)*y*(3*x**2-y**2)*(5*r**2-4)
	Z20 =  Z[20] * __sqrt__(12)*x*(16*x**4-20*x**2*r**2+5*r**4)
	Z21 =  Z[21] * __sqrt__(12)*y*(16*y**4-20*y**2*r**2+5*r**4)
	Z22 =  Z[22] * __sqrt__(7)*(20*r**6-30*r**4+12*r**2-1)
	Z23 =  Z[23] * 2*__sqrt__(14)*x*y*(15*r**4-20*r**2+6)
	Z24 =  Z[24] * __sqrt__(14)*(x**2-y**2)*(15*r**4-20*r**2+6)
	Z25 =  Z[25] * 4*__sqrt__(14)*x*y*(x**2-y**2)*(6*r**2-5)
	Z26 =  Z[26] * __sqrt__(14)*(8*x**4-8*x**2*r**2+r**4)*(6*r**2-5)
	Z27 =  Z[27] * __sqrt__(14)*x*y*(32*x**4-32*x**2*r**2+6*r**4)
	Z28 =  Z[28] * __sqrt__(14)*(32*x**6-48*x**4*r**2+18*x**2*r**4-r**6)
	Z29 =  Z[29] * 4*y*(35*r**6-60*r**4+30*r**2+10)
	Z30 =  Z[30] * 4*x*(35*r**6-60*r**4+30*r**2+10)
	Z31 =  Z[31] * 4*y*(3*x**2-y**2)*(21*r**4-30*r**2+10)
	Z32 =  Z[32] * 4*x*(x**2-3*y**2)*(21*r**4-30*r**2+10)
	Z33 =  Z[33] * 4*(7*r**2-6)*(4*x**2*y*(x**2-y**2)+y*(r**4-8*x**2*y**2))
	Z34 =  Z[34] * (4*(7*r**2-6)*(x*(r**4-8*x**2*y**2)-4*x*y**2*(x**2-y**2)))
	Z35 =  Z[35] * (8*x**2*y*(3*r**4-16*x**2*y**2)+4*y*(x**2-y**2)*(r**4-16*x**2*y**2))
	Z36 =  Z[36] * (4*x*(x**2-y**2)*(r**4-16*x**2*y**2)-8*x*y**2*(3*r**4-16*x**2*y**2))
	Z37 =  Z[37] * 3*(70*r**8-140*r**6+90*r**4-20*r**2+1)
	ZW = 	Z1 + Z2 +  Z3+  Z4+  Z5+  Z6+  Z7+  Z8+  Z9+ \
			Z10+ Z11+ Z12+ Z13+ Z14+ Z15+ Z16+ Z17+ Z18+ Z19+ \
			Z20+ Z21+ Z22+ Z23+ Z24+ Z25+ Z26+ Z27+ Z28+ Z29+ \
			Z30+ Z31+ Z32+ Z33+ Z34+ Z35+ Z36+ Z37
	return ZW
