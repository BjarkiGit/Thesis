import numpy as np



def gaussian17(x, z, wid, amp0, cen0, amp1, cen1, amp2, cen2, amp3, cen3, amp4, cen4,
            amp5, cen5, amp6, cen6, amp7, cen7, amp8, cen8, amp9, cen9, amp10,
             cen10, amp11, cen11, amp12, cen12, amp13, cen13, amp14, cen14, amp15, cen15,
             amp16, cen16, amp17, cen17):
             
    return amp0 * np.exp(-(x-cen0*(1+z))**2 / (2*wid**2)) +\
           amp1 * np.exp(-(x-cen1*(1+z))**2 / (2*wid**2)) +\
           amp2 * np.exp(-(x-cen2*(1+z))**2 / (2*wid**2)) +\
           amp3 * np.exp(-(x-cen3*(1+z))**2 / (2*wid**2)) +\
           amp4 * np.exp(-(x-cen4*(1+z))**2 / (2*wid**2)) +\
           amp5 * np.exp(-(x-cen5*(1+z))**2 / (2*wid**2)) +\
           amp6 * np.exp(-(x-cen6*(1+z))**2 / (2*wid**2)) +\
           amp7 * np.exp(-(x-cen7*(1+z))**2 / (2*wid**2)) +\
           amp8 * np.exp(-(x-cen8*(1+z))**2 / (2*wid**2)) +\
           amp9 * np.exp(-(x-cen9*(1+z))**2 / (2*wid**2)) +\
           amp10 * np.exp(-(x-cen10*(1+z))**2 / (2*wid**2)) +\
           amp11 * np.exp(-(x-cen11*(1+z))**2 / (2*wid**2)) +\
           amp12 * np.exp(-(x-cen12*(1+z))**2 / (2*wid**2)) +\
           amp13 * np.exp(-(x-cen13*(1+z))**2 / (2*wid**2)) +\
           amp14 * np.exp(-(x-cen14*(1+z))**2 / (2*wid**2)) +\
           amp15 * np.exp(-(x-cen15*(1+z))**2 / (2*wid**2)) +\
           amp16 * np.exp(-(x-cen16*(1+z))**2 / (2*wid**2)) +\
           amp17 * np.exp(-(x-cen17*(1+z))**2 / (2*wid**2)) 
           
           
           
def gaussian(x, z, wid, strength, cen):
	i = 0
	st = 100 #Strengths are given as a fraction of Halpha strength,
		 #so this is initial guess for Halpha
	amp = np.zeros_like(strength)
	for i in range(len(strength)):
		a = st*strength[i]
		vars()['amp' + str(i)] = a
		amp[i] = vars()
		#cen = wl_vac[i]
		#vars()['cen' + str(i)] = cen
		
	gauss = np.zeros_like(amp)
	for i in range(len(amp)):
		a = amp[i]
		c = cen[i]
		g =  a * np.exp(-(x-c*(1+z))**2 / (2*wid**2))
		gauss[i] = g
		i += 1
			
	return gauss
	
def gaussTest(x, z, wid, amp, cen):
	i = 0
	for i in range(len(amp)):
		return amp[i] * np.exp(-(x-cen[i]*(1+z))**2 / (2*wid**2))
		i += 1
		
		
def gaussTest2(pars, x, data=None):
	vals = pars.values(dict)
	z = vals["z"]
	wid = vals["wid"]
	amp = vals["amp"]
	cen = vals["cen"]
	model =  amp * np.exp(-(x-cen*(1+z))**2 / (2*wid**2))
	return model - data
	


def gaussTest3(x, z, wid, amp, cen):
	return amp * np.exp(-(x-cen*(1+z))**2 / (2*wid**2))
	
	
def gauss0(x, z, wid, amp, cen, flux):
	
	i = 0
	G = np.array([])
	for a in amp:
		c = cen[i]
		g = gaussTest3(x, z, wid, a, c)
		if len(G) == 0:
			G = np.append(G, g)
		else:
			G += g
		i += 1
	return G-flux
	

def gauss(pars, x, f=None, lines=None):
	cen = lines["wl"].values
	names = lines["line"]
	vals = pars.valuesdict()
	z = vals["z"]
	wid = vals["wid"]
	i = 0
	G = np.array([])	
	if len(cen) == 0:
		print("No lines to model")
		return 0

	for key, value in pars.items():
		
		if key != "z" and key != "wid":
			c = cen[i]
			a = value.value			
			g = gaussTest3(x, z, wid, a, c)
			i += 1
			if len(G) == 0:
				G = np.append(G, g)
			else:
				# G = np.maximum(g, G)
				G += g
				# or G += g, not sure which would be better

		else:
			pass
		
	if f is None:
		print("No flux given, returning model")
		return G
			
	return f-G