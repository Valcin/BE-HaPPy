import numpy as np
import scipy.optimize as op
import matplotlib.pyplot as plt
import emcee
import sys
sys.path.append('/home/david/codes/FAST-PT')
import myFASTPT as FPT
import math
import time
from scipy.optimize import curve_fit
from scipy.special import erf



########################################################################
########### POWER LAW even
########################################################################

def coeffit_pl2 (kstop,lb1, errlb1, pop, k ,b ,errb):
	#~ lim = np.where(k < kstop)[0]
	lim = np.where((k < kstop)&(k > 1e-2))[0]
	def lnlike(theta, x, y, yerr):
		b1, b2, b4 = theta
		model = b1 + b2*x[lim]**2 + b4*x[lim]**4 
		inv_sigma2 = 1.0/(yerr[lim]**2)
		return -0.5*(np.sum((y[lim]-model)**2*inv_sigma2 - np.log(inv_sigma2)))
	
	def lnprior(theta):
		b1, b2, b4 = theta
		if lb1 - 3*errlb1 < b1 < lb1 + 3*errlb1  and b1 > 0:
			return 0.0
		return -np.inf
	
	def lnprob(theta, x, y, yerr):
		lp = lnprior(theta)
		if not np.isfinite(lp):
			return -np.inf
		return lp + lnlike(theta, x, y, yerr)
		
		
		
	

	nll = lambda *args: -lnlike(*args)
	result = op.minimize(nll, [pop],  method='Nelder-Mead', args=(k, b ,errb ),  options={'maxfev': 2000} )
	b1_ml, b2_ml, b4_ml = result["x"]
	#~ print pop
	#~ print result
	#~ max_l = lnlike(result["x"], k, b, errb )
	#~ AIC = 2*4. - 2 * max_l
	#~ print 'maximum likelihood is '+str(max_l)
	#~ print 'AIC = '+str(AIC)
	
	#~ ndim, nwalkers = len(pop), 300
	#~ pos = [result["x"] + 1e-4*np.random.randn(ndim) for i in range(nwalkers)]
	
	#~ sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(k, b, errb))
	#~ sampler.run_mcmc(pos, 1000)
	
	#~ samples = sampler.chain[:, 50:, :].reshape((-1, ndim))

	
	#~ import corner
	#~ fig = corner.corner(samples, labels=["$b1$", "$b2$", "$b3$", "$b4$" ], truths=[b1_ml, b2_ml, b3_ml, b4_ml])
	#~ fig.savefig("/home/david/triangle.png")
	

	#~ b1_mcmc, b2_mcmc, b4_mcmc = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]), zip(*np.percentile(samples, [16, 50, 84], axis=0)))
	
	#~ print b1_mcmc, b2_mcmc, b4_mcmc
	

	#~ plt.figure()
	#~ ax1 = plt.subplot(221)
	#~ ax1.set_title('b1')
	#~ for i in xrange(0,nwalkers):
		#~ ax1.plot(np.arange(1000), sampler.chain[i,:,0])
	#~ ax2 = plt.subplot(222)
	#~ ax2.set_title('b2')
	#~ for i in xrange(0,nwalkers):
		#~ ax2.plot(np.arange(1000), sampler.chain[i,:,1])
	#~ ax3 = plt.subplot(223)
	#~ ax3.set_title('b3')
	#~ for i in xrange(0,nwalkers):
		#~ ax3.plot(np.arange(1000), sampler.chain[i,:,2])
	#~ plt.show()
	
	#~ return b1_mcmc, b2_mcmc, b4_mcmc
	return b1_ml, b2_ml, b4_ml
########################################################################
########### POWER LAW odd
########################################################################

def coeffit_pl (kstop,lb1, errlb1, pop, k ,b ,errb):
	#~ lim = np.where(k < kstop)[0]
	lim = np.where((k < kstop)&(k > 1e-2))[0]
	def lnlike(theta, x, y, yerr):
		b1, b2, b3, b4 = theta
		model = b1 + b2*x[lim]**2 + b3*x[lim]**3 + b4*x[lim]**4 
		inv_sigma2 = 1.0/(yerr[lim]**2)
		#~ print theta
		#~ print -0.5*(np.sum((y[lim]-model)**2*inv_sigma2 - np.log(inv_sigma2)))
		return -0.5*(np.sum((y[lim]-model)**2*inv_sigma2 - np.log(inv_sigma2)))
	
	def lnprior(theta):
		b1, b2, b3, b4 = theta
		if lb1 - 3*errlb1 < b1 < lb1 + 3*errlb1  and b1 > 0:
			return 0.0
		return -np.inf
	
	def lnprob(theta, x, y, yerr):
		lp = lnprior(theta)
		if not np.isfinite(lp):
			return -np.inf
		#~ if lnlike(theta, x, y, yerr) < 0:
			#~ print ('-- ') + str(lp + lnlike(theta, x, y, yerr))
			#~ print ('-- ') + str(lp + lnlike(theta, x, y, yerr))
			#~ from termcolor import colored
			#~ print colored(theta, 'red')
		#~ else:
			#~ print '!! '+ str(lp + lnlike(theta, x, y, yerr))
			#~ print '!! '+ str(lp + lnlike(theta, x1,x2,x3,x4,x5,x6,x7,x8, y, yerr))
			#~ from termcolor import colored
			#~ print colored(theta, 'blue')
		return lp + lnlike(theta, x, y, yerr)
		
		
		
	

	nll = lambda *args: -lnlike(*args)
	result = op.minimize(nll, [pop],  method='Nelder-Mead', args=(k, b ,errb ),  options={'maxfev': 2000} )
	b1_ml, b2_ml, b3_ml, b4_ml = result["x"]
	#~ print pop
	#~ print result
	#~ max_l = lnlike(result["x"], k, b, errb )
	#~ AIC = 2*4. - 2 * max_l
	#~ print 'maximum likelihood is '+str(max_l)
	#~ print 'AIC = '+str(AIC)
	
	#~ ndim, nwalkers = len(pop), 300
	#~ pos = [result["x"] + 1e-3*np.random.randn(ndim) for i in range(nwalkers)]
	
	#~ posnum = 1000
	#~ sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(k, b, errb))
	#~ sampler.run_mcmc(pos, posnum)
	
	#~ samples = sampler.chain[:, 50:, :].reshape((-1, ndim))

	#~ print("Mean acceptance fraction: {0:.3f}"
                #~ .format(np.mean(sampler.acceptance_fraction)))
	#~ import corner
	#~ cmap = ['w','b','g','c']
	#~ fig = corner.corner(samples, labels=[r"$b_1$", r"$b_2$", r"$b_3$", r"$b_4$" ], truths=[b1_ml, b2_ml, b3_ml, b4_ml], fontsize = 14)
	#~ fig.savefig("/home/david/triangle.png")
	

	#~ b1_mcmc, b2_mcmc, b3_mcmc, b4_mcmc = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]), zip(*np.percentile(samples, [16, 50, 84], axis=0)))
	
	#~ print b1_mcmc, b2_mcmc, b3_mcmc, b4_mcmc
	
	#~ import seaborn as sns
	#~ import pandas as pd
	#~ sns.set(style="white")
	#~ iris = sns.load_dataset("iris")
	#~ print iris
	#~ wnv_approaches = pd.DataFrame({"b1": samples[:,0], "b2": samples[:,1], "b3": samples[:,2],"b4": samples[:,3]})
	#~ print wnv_approaches
	#~ g = sns.PairGrid(wnv_approaches)
	#~ g = g.map_upper(plt.scatter)
	#~ g = g.map_lower(sns.kdeplot, cmap="Blues_d")
	#~ g = g.map_diag(sns.kdeplot, lw=3, legend=False)
	
	
	#~ kill
	#~ plt.figure()
	#~ ax1 = plt.subplot(221)
	#~ ax1.set_title('b1')
	#~ for i in xrange(0,nwalkers):
		#~ ax1.plot(np.arange(1000), sampler.chain[i,:,0])
	#~ ax1.plot(np.arange(posnum), np.mean(sampler.chain[:,:,0], axis =0))
	#~ ax2 = plt.subplot(222)
	#~ ax2.set_title('b2')
	#~ for i in xrange(0,nwalkers):
		#~ ax2.plot(np.arange(1000), sampler.chain[i,:,1])
	#~ ax2.plot(np.arange(posnum), np.mean(sampler.chain[:,:,1], axis =0))
	#~ ax3 = plt.subplot(223)
	#~ ax3.set_title('b3')
	#~ for i in xrange(0,nwalkers):
		#~ ax3.plot(np.arange(1000), sampler.chain[i,:,2])
	#~ ax3.plot(np.arange(posnum), np.mean(sampler.chain[:,:,2], axis =0))
	#~ ax4 = plt.subplot(224)
	#~ ax4.set_title('b4')
	#~ for i in xrange(0,nwalkers):
		#~ ax4.plot(np.arange(1000), sampler.chain[i,:,3])
	#~ ax4.plot(np.arange(posnum), np.mean(sampler.chain[:,:,3], axis =0))
	#~ plt.show()
	
	#~ return b1_mcmc, b2_mcmc, b3_mcmc, b4_mcmc
	#~ return b1_ml, b2_ml, b3_ml, b4_ml, b1_mcmc, b2_mcmc, b3_mcmc, b4_mcmc
	return b1_ml, b2_ml, b3_ml, b4_ml
########################################################################
######### bias expansion 2nd order
########################################################################
def coeffit_exp1(kstop, Pmm, A, B, C, D, E, lb1, errlb1, pop, k ,b ,errb, noise):
	#~ lim = np.where(k < kstop)[0]
	lim = np.where((k < kstop)&(k > 1e-2))[0]

	#~ from scipy.signal import savgol_filter
	#~ b = savgol_filter(b, 51, 3) 
	#~ errb = savgol_filter(errb, 51, 6) 
	
	def lnlike(theta, x, y, yerr):
		b1, b2, bs, N = theta
		model = np.sqrt((b1**2 * Pmm[lim]+ b1*b2*A[lim] + 1/4.*b2**2*B[lim] + b1*bs*C[lim] + 1/2.*b2*bs*D[lim] + \
		1/4.*bs**2*E[lim] + N)/Pmm[lim])
		inv_sigma2 = 1.0/(yerr[lim]**2)
		return -0.5*(np.sum((y[lim]-model)**2*inv_sigma2 - np.log(inv_sigma2))) 
	
	def lnprior(theta):
		b1, b2, bs, N = theta
		#~ if lb1 - 3*errlb1 < b1 < lb1 + 3*errlb1  and b1 > 0:
		if  0 < b1  and -3*noise < N < 3*noise:
			return 0.0
		#~ return 0.0
		return -np.inf
	
	def lnprob(theta, x, y, yerr):
		lp = lnprior(theta)
		if not np.isfinite(lp):
			return -np.inf
		#~ if lnlike(theta, x, y, yerr) < 0:
			#~ print ('-- ') + str(lp + lnlike(theta, x, y, yerr))
			#~ print ('-- ') + str(lp + lnlike(theta, x1,x2,x3,x4,x5,x6,x7,x8, y, yerr))
			#~ from termcolor import colored
			#~ print colored(theta, 'red')
		#~ else:
			#~ print '!! '+ str(lp + lnlike(theta, x, y, yerr))
			#~ print '!! '+ str(lp + lnlike(theta, x1,x2,x3,x4,x5,x6,x7,x8, y, yerr))
			#~ from termcolor import colored
			#~ print colored(theta, 'blue')
		return lp + lnlike(theta, x, y, yerr)


	nll = lambda *args: -lnlike(*args)
	#~ result = op.minimize(nll, [pop],  method='Nelder-Mead', args=(k, b ,errb ),  options={'maxfev': 2000} )
	result = op.minimize(nll, [pop],method='SLSQP', args=(k, b ,errb ),  bounds=((0, None), (None,None), (None,None), (-3*noise,3*noise)),  options={'maxfev': 5000} )
	b1_ml, b2_ml, bs_ml, N_ml = result["x"]
	#~ print pop
	#~ print(result)
	
	#~ max_l = lnlike(result["x"], k, b, errb )
	#~ AIC = 2*4. - 2 * max_l
	#~ print 'maximum likelihood is '+str(max_l)
	#~ print 'AIC = '+str(AIC)
	
	#~ ndim, nwalkers = len(pop), 100
	#~ pos = [result["x"] + 1e-4*np.random.randn(ndim) for i in range(nwalkers)]

	#~ posnum = 1000
	#~ sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(k, b, errb), a = 2.)
	#~ sampler.run_mcmc(pos, posnum)
	
	#~ samples = sampler.chain[:, 200:, :].reshape((-1, ndim))

	#~ print("Mean acceptance fraction: {0:.3f}"
                #~ .format(np.mean(sampler.acceptance_fraction)))
	#~ import corner
	#~ fig = corner.corner(samples, labels=["$b1$", "$b2$", "$bs$" , "$N$" ], truths=[b1_ml, b2_ml, bs_ml, N_ml])
	#~ fig.savefig("/home/david/triangle.png")
	

	#~ b1_mcmc, b2_mcmc, bs_mcmc, N_mcmc = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]), zip(*np.percentile(samples, [16, 50, 84], axis=0)))
	#~ plt.figure()
	#~ ax1 = plt.subplot(221)
	#~ ax1.set_title('b1')
	#~ for i in xrange(0,nwalkers):
		#~ ax1.plot(np.arange(posnum), sampler.chain[i,:,0])
	#~ ax1.plot(np.arange(posnum), np.mean(sampler.chain[:,:,0], axis =0),color='k')
	#~ ax1.axhline(b1_mcmc[0], color='k', linestyle='--')
	#~ ax1.axhline(b1_ml, color='k')
	#~ ax2 = plt.subplot(222)
	#~ ax2.set_title('b2')
	#~ for i in xrange(0,nwalkers):
		#~ ax2.plot(np.arange(posnum), sampler.chain[i,:,1])
	#~ ax2.plot(np.arange(posnum), np.mean(sampler.chain[:,:,1], axis =0),color='k')
	#~ ax2.axhline(b2_mcmc[0], color='k', linestyle='--')
	#~ ax2.axhline(b2_ml, color='k')
	#~ ax3 = plt.subplot(223)
	#~ ax3.set_title('bs')
	#~ for i in xrange(0,nwalkers):
		#~ ax3.plot(np.arange(posnum), sampler.chain[i,:,2])
	#~ ax3.plot(np.arange(posnum), np.mean(sampler.chain[:,:,2], axis =0),color='k')
	#~ ax3.axhline(bs_mcmc[0], color='k', linestyle='--')
	#~ ax3.axhline(bs_ml, color='k')

	#~ plt.show()

	#~ return b1_mcmc, b2_mcmc, bs_mcmc, N_mcmc, b1_ml, b2_ml, bs_ml, N_ml
	#~ return b1_ml, b2_ml, bs_ml, N_ml, b1_mcmc, b2_mcmc, bs_mcmc, N_mcmc
	return b1_ml, b2_ml, bs_ml, N_ml
########################################################################
######### bias expansion 3rd order free
########################################################################
def coeffit_exp2(kstop, Pmm, A, B, C, D, E, F, lb1, errlb1, pop, k ,b ,errb, noise):
	
	#~ lim = np.where(k < kstop)[0]
	lim = np.where((k < kstop)&(k > 1e-2))[0]
	#~ from scipy.signal import savgol_filter
	#~ b = savgol_filter(b, 51, 3) 
	#~ errb = savgol_filter(errb, 51, 6) 
	
	def lnlike(theta, x, y, yerr):
		b1, b2, bs, b3nl, N = theta
		model = np.sqrt((b1**2 * Pmm[lim]+ b1*b2*A[lim] + 1/4.*b2**2*B[lim] + b1*bs*C[lim] + 1/2.*b2*bs*D[lim] + \
		1/4.*bs**2*E[lim] + 2*b1*b3nl*F[lim] + N)/Pmm[lim])
		inv_sigma2 = 1.0/(yerr[lim]**2)
		return -0.5*(np.sum((y[lim]-model)**2*inv_sigma2 - np.log(inv_sigma2)))
	
	def lnprior(theta, x, y, yerr):
		b1, b2, bs, b3nl, N = theta
		if  0 < b1 and  -3*noise < N < 3*noise :
			return 0.0
		return -np.inf
		
	def lnprob(theta, x, y, yerr):
		lp = lnprior(theta, x, y, yerr)
		if not np.isfinite(lp):
			return -np.inf
		#~ if lnlike(theta, x, y, yerr) < 0:
			#~ print ('-- ') + str(lp + lnlike(theta, x, y, yerr))
			#~ print ('-- ') + str(lp + lnlike(theta, x1,x2,x3,x4,x5,x6,x7,x8, y, yerr))
			#~ from termcolor import colored
			#~ print colored(theta, 'red')
		#~ else:
			#~ print '!! '+ str(lp + lnlike(theta, x, y, yerr))
			#~ print '!! '+ str(lp + lnlike(theta, x1,x2,x3,x4,x5,x6,x7,x8, y, yerr))
			#~ from termcolor import colored
			#~ print colored(theta, 'blue')
		return lp + lnlike(theta, x, y, yerr)

	
	nll = lambda *args: -lnlike(*args)
	#~ result = op.minimize(nll, [pop],  method='Nelder-Mead', args=(k, b ,errb ),  options={'maxfev': 5000} )
	result = op.minimize(nll, [pop],method='SLSQP', args=(k, b ,errb ),  bounds=((0, None), (None,None), (None,None), (None,None), (-3*noise,3*noise)),  options={'maxfev': 5000} )
	b1_ml, b2_ml, bs_ml, b3_ml, N_ml = result["x"]
	#~ print pop
	#~ print(result)
	
	#~ max_l = lnlike(result["x"], k, b, errb )
	#~ AIC = 2*5. - 2 * max_l
	#~ print 'maximum likelihood is '+str(max_l)
	#~ print 'AIC = '+str(AIC)
	
	#~ ndim, nwalkers = len(pop), 300
	#~ pos = [result["x"] + 1e-4*np.random.randn(ndim) for i in range(nwalkers)]

	#~ posnum = 1000
	#~ sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(k, b, errb), a = 2.)
	#~ sampler.run_mcmc(pos, posnum)
	
	#~ samples = sampler.chain[:, 200:, :].reshape((-1, ndim))

	#~ print("Mean acceptance fraction: {0:.3f}"
                #~ .format(np.mean(sampler.acceptance_fraction)))
	#~ import corner
	#~ samplesbis = np.delete(samples, 4, 1)
	#~ fig = corner.corner(samplesbis, labels=[r"$b_1$", r"$b_2$", r"$b_{s2}$" ,r"$b_{3nl}$"],\
	#~ truths=[b1_ml, b2_ml, bs_ml, b3_ml], fontsize = 14)
	
	#~ fig.savefig("/home/david/triangle.png")
	

	#~ b1_mcmc, b2_mcmc, bs_mcmc, b3_mcmc, N_mcmc = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]), zip(*np.percentile(samples, [16, 50, 84], axis=0)))
	#~ print(result)
	#~ print b1_mcmc[0], b2_mcmc[0], bs_mcmc[0], b3_mcmc[0], N_mcmc[0]
	#~ print pop
	#~ print b1_ml, b2_ml, bs_ml, b3_ml, N_ml
	

	#~ plt.figure()
	#~ ax1 = plt.subplot(221)
	#~ ax1.set_title('b1')
	#~ for i in xrange(0,nwalkers):
		#~ ax1.plot(np.arange(posnum), sampler.chain[i,:,0])
	#~ ax1.plot(np.arange(posnum), np.mean(sampler.chain[:,:,0], axis =0),color='k')
	#~ ax1.axhline(b1_mcmc[0], color='k', linestyle='--')
	#~ ax1.axhline(b1_ml, color='k')
	#~ ax2 = plt.subplot(222)
	#~ ax2.set_title('b2')
	#~ for i in xrange(0,nwalkers):
		#~ ax2.plot(np.arange(posnum), sampler.chain[i,:,1])
	#~ ax2.plot(np.arange(posnum), np.mean(sampler.chain[:,:,1], axis =0),color='k')
	#~ ax2.axhline(b2_mcmc[0], color='k', linestyle='--')
	#~ ax2.axhline(b2_ml, color='k')
	#~ ax3 = plt.subplot(223)
	#~ ax3.set_title('bs')
	#~ for i in xrange(0,nwalkers):
		#~ ax3.plot(np.arange(posnum), sampler.chain[i,:,2])
	#~ ax3.plot(np.arange(posnum), np.mean(sampler.chain[:,:,2], axis =0),color='k')
	#~ ax3.axhline(bs_mcmc[0], color='k', linestyle='--')
	#~ ax3.axhline(bs_ml, color='k')
	#~ ax4 = plt.subplot(224)
	#~ ax4.set_title('b3nl')
	#~ for i in xrange(0,nwalkers):
		#~ ax4.plot(np.arange(posnum), sampler.chain[i,:,3])
	#~ ax4.plot(np.arange(posnum), np.mean(sampler.chain[:,:,3], axis =0),color='k')
	#~ ax4.axhline(b3_mcmc[0], color='k', linestyle='--')
	#~ ax4.axhline(b3_ml, color='k')

	#~ plt.show()
	#~ kill
	#~ return b1_mcmc, b2_mcmc, bs_mcmc, b3_mcmc, N_mcmc
	#~ return b1_ml, b2_ml, bs_ml, b3_ml, N_ml, b1_mcmc, b2_mcmc, bs_mcmc, b3_mcmc, N_mcmc
	return b1_ml, b2_ml, bs_ml, b3_ml, N_ml
########################################################################
######### bias expansion 3rd order fixed
########################################################################
def coeffit_exp3(kstop, Pmm, A, B, C, D, E, F, lb1, errlb1, pop, k ,b ,errb, noise):
	
	#~ lim = np.where(k < kstop)[0]
	lim = np.where((k < kstop)&(k > 1e-2))[0]
	def lnlike(theta, x, y, yerr):
		b1, b2, N = theta
		bs = -4/7.*(b1-1)
		b3nl = 32/315.*(b1-1)
		model = np.sqrt((b1**2 * Pmm[lim]+ b1*b2*A[lim] + 1/4.*b2**2*B[lim] + b1*bs*C[lim] + 1/2.*b2*bs*D[lim] + 1/4.*bs**2*E[lim]\
		+ 2*b1*b3nl*F[lim] + N)/Pmm[lim])
		inv_sigma2 = 1.0/(yerr[lim]**2)
		return -0.5*(np.sum((y[lim]-model)**2*inv_sigma2 - np.log(inv_sigma2)))
	
	def lnprior(theta):
		b1, b2, N = theta
		if  0 < b1 and  -3*noise < N < 3*noise :
			return 0.0
		#~ return 0.0
		return -np.inf
	
	def lnprob(theta, x, y, yerr):
		lp = lnprior(theta)
		if not np.isfinite(lp):
			return -np.inf
		if lnlike(theta, x, y, yerr) < 0:
			print ('-- ') + str(lp + lnlike(theta, x, y, yerr))
			#~ print ('-- ') + str(lp + lnlike(theta, x1,x2,x3,x4,x5,x6,x7,x8, y, yerr))
			from termcolor import colored
			print colored(theta, 'red')
		else:
			print '!! '+ str(lp + lnlike(theta, x, y, yerr))
			#~ print '!! '+ str(lp + lnlike(theta, x1,x2,x3,x4,x5,x6,x7,x8, y, yerr))
			from termcolor import colored
			print colored(theta, 'blue')
		return lp + lnlike(theta, x, y, yerr)
		
		
		
	

	nll = lambda *args: -lnlike(*args)
	#~ result = op.minimize(nll, [pop],  method='Nelder-Mead', args=(k, b ,errb ),  options={'maxfev': 2000} )
	#~ result = op.minimize(nll, [pop],  method='Nelder-Mead', args=(k, b ,errb ),  options={'maxfev': 2000} )
	result = op.minimize(nll, [pop],method='SLSQP', args=(k, b ,errb ),  bounds=((0, None), (None,None), (-3*noise,3*noise)),  options={'maxfev': 5000} )
	b1_ml, b2_ml, N_ml = result["x"]
	
	
	#~ max_l = lnlike(result["x"], k, b, errb )
	#~ AIC = 2*3. - 2 * max_l
	#~ print 'maximum likelihood is '+str(max_l)
	#~ print 'AIC = '+str(AIC)
	
	#~ ndim, nwalkers = len(pop), 300
	#~ pos = [result["x"] + 1e-4*np.random.randn(ndim) for i in range(nwalkers)]

	#~ posnum = 1000
	#~ sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(k, b, errb))
	#~ sampler.run_mcmc(pos, posnum)
	
	#~ samples = sampler.chain[:, 200:, :].reshape((-1, ndim))

	
	#~ import corner
	#~ fig = corner.corner(samples, labels=["$b1$", "$b2$", "$N$" ], truths=[b1_ml, b2_ml, N_ml])
	#~ fig.savefig("/home/david/triangle.png")
	

	#~ b1_mcmc, b2_mcmc, N_mcmc = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]), zip(*np.percentile(samples, [16, 50, 84], axis=0)))
	
	#~ print("Mean acceptance fraction: {0:.3f}"
                #~ .format(np.mean(sampler.acceptance_fraction)))
	#~ print b1_mcmc, b2_mcmc, N_mcmc
	#~ print pop
	#~ print result

	#~ plt.figure()
	#~ ax1 = plt.subplot(221)
	#~ ax1.set_title('b1')
	#~ for i in xrange(0,nwalkers):
		#~ ax1.plot(np.arange(1000), sampler.chain[i,:,0])
	#~ ax2 = plt.subplot(222)
	#~ ax2.set_title('b2')
	#~ for i in xrange(0,nwalkers):
		#~ ax2.plot(np.arange(1000), sampler.chain[i,:,1])
	#~ ax3 = plt.subplot(223)
	#~ ax3.set_title('N')
	#~ for i in xrange(0,nwalkers):
		#~ ax3.plot(np.arange(1000), sampler.chain[i,:,2])
	

	#~ plt.show()
	#~ return b1_mcmc, b2_mcmc, N_mcmc
	return b1_ml, b2_ml, N_ml

########################################################################
#### Linear Kaiser
########################################################################
def coeffit_Kaiser(j, fcc, kstop, Pmm, lb1, k ,b ,errb):
	

	start = time.time()
	#~ lim = np.where(k < kstop)[0]
	lim = np.where((k < kstop)&(k > 1e-2))[0]

	def lnlike(theta, x, y, yerr):
		sigma = theta
		#~ kappa = x[lim]*sigma
		#~ coeffA = math.sqrt(math.pi)/2. * erf(kappa)/kappa
		#~ coeffB = 3./2./kappa**2*(coeffA - np.exp(-kappa**2))
		#~ coeffC = 5./2./kappa**2*(coeffB - np.exp(-kappa**2))
		kappa = x[lim]*sigma*fcc[lim]*Dz[ind]
		coeffA = np.arctan(kappa/math.sqrt(2))/(math.sqrt(2)*kappa) + 1/(2+kappa**2)
		coeffB = 6/kappa**2*(coeffA - 2/(2+kappa**2))
		coeffC = -10/kappa**2*(coeffB - 2/(2+kappa**2))
		model = Pmm[lim]*(lb1[lim]**2*coeffA +  2/3.*lb1[lim]*fcc[lim]*coeffB + 1/5.*fcc[lim]**2*coeffC)
		#~ model = Pmm[lim]*(lb1**2*coeffA +  2/3.*lb1*fcc[lim]*coeffB + 1/5.*fcc[lim]**2)
		inv_sigma2 = 1.0/(yerr[lim]**2)
		return -0.5*(np.sum((y[lim]-model)**2*inv_sigma2 - np.log(inv_sigma2)))
	
	def lnprior(theta):
		sigma = theta
		#~ if lb1 - 3*errlb1 < b1 < lb1 + 3*errlb1 :
		if 0 < sigma < 100:
			return 0.0
		return -np.inf
	
	def lnprob(theta, x, y, yerr):
		lp = lnprior(theta)
		if not np.isfinite(lp):
			return -np.inf
		return lp + lnlike(theta, x, y, yerr)
		
	z = [0.0,0.5,1.0,2.0]
	red = ['0.0','0.5','1.0','2.0','3.0']
	ind = red.index(str(z[j]))
	f = [0.518,0.754,0.872,0.956,0.98]
	Dz = [ 1.,0.77,0.61,0.42]

	nll = lambda *args: -lnlike(*args)
	#~ result = op.minimize(nll, [pop], bounds= [(1,1000)], args=(k, b ,errb ))
	result = op.minimize_scalar(nll, bounds=(1,1000), method='bounded',args=(k, b ,errb ))
	b1_ml = result["x"]
	print(result)
	
	#~ max_l = lnlike(result["x"], k, b, errb )
	#~ AIC = 2*4. - 2 * max_l
	#~ print 'maximum likelihood is '+str(max_l)
	#~ print 'AIC = '+str(AIC)
	
	#~ ndim, nwalkers = 1, 200
	#~ pos = [result["x"] + 1e-4*np.random.randn(ndim) for i in range(nwalkers)]



	#~ sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(k, b, errb))
	#~ sampler.run_mcmc(pos, 1000)
	#~ samples = sampler.chain[:, 200:, :].reshape((-1, ndim))

	
	#~ import corner
	#~ fig = corner.corner(samples, labels=["$b1$", "$b2$", "$bs$", "$b3nl$" ], truths=[b1_ml, b2_ml, bs_ml, b3nl])
	#~ fig.savefig("/home/david/triangle.png")
	

	#~ b1_mcmc = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]), zip(*np.percentile(samples, [16, 50, 84], axis=0)))
	
	#~ print b1_mcmc
	

	#~ plt.figure()
	#~ plt.title('sigma kaiser')
	#~ for i in xrange(0,nwalkers):
		#~ plt.plot(np.arange(1000), sampler.chain[i,:,0])
	#~ plt.show()
	
	end = time.time()
	print 'time is '+str((end - start))
	#~ return b1_mcmc
	return b1_ml
########################################################################
#### Linear Kaiser 2
########################################################################
def coeffit_Kaiser2(j, fcc, kstop, Pmm,  k ,b ,errb):
	

	start = time.time()
	#~ lim = np.where(k < kstop)[0]
	lim = np.where((k < kstop)&(k > 1e-2))[0]

	def lnlike(theta, x, y, yerr):
		b1, b2, b3, b4, sigma = theta
		#~ kappa = x[lim]*sigma
		#~ coeffA = math.sqrt(math.pi)/2. * erf(kappa)/kappa
		#~ coeffB = 3./2./kappa**2*(coeffA - np.exp(-kappa**2))
		#~ coeffC = 5./2./kappa**2*(coeffB - np.exp(-kappa**2))
		bcc = b1 + b2*x[lim]**2 + b3*x[lim]**3 + b4*x[lim]**4 
		kappa = x[lim]*sigma*fcc[lim]*Dz[ind]
		coeffA = np.arctan(kappa/math.sqrt(2))/(math.sqrt(2)*kappa) + 1/(2+kappa**2)
		coeffB = 6/kappa**2*(coeffA - 2/(2+kappa**2))
		coeffC = -10/kappa**2*(coeffB - 2/(2+kappa**2))
		model = Pmm[lim]*(bcc**2*coeffA +  2/3.*bcc*fcc[lim]*coeffB + 1/5.*fcc[lim]**2*coeffC)
		#~ model = Pmm[lim]*(lb1**2*coeffA +  2/3.*lb1*fcc[lim]*coeffB + 1/5.*fcc[lim]**2)
		inv_sigma2 = 1.0/(yerr[lim]**2)
		return -0.5*(np.sum((y[lim]-model)**2*inv_sigma2 - np.log(inv_sigma2)))
	
	def lnprior(theta):
		b1, b2, b3, b4, sigma = theta
		#~ if lb1 - 3*errlb1 < b1 < lb1 + 3*errlb1 :
		if 0 < sigma < 100:
			return 0.0
		return -np.inf
	
	def lnprob(theta, x, y, yerr):
		lp = lnprior(theta)
		if not np.isfinite(lp):
			return -np.inf
		return lp + lnlike(theta, x, y, yerr)
		
	z = [0.0,0.5,1.0,2.0]
	red = ['0.0','0.5','1.0','2.0','3.0']
	ind = red.index(str(z[j]))
	f = [0.518,0.754,0.872,0.956,0.98]
	Dz = [ 1.,0.77,0.61,0.42]

	pop = [1,1,1,1,6]
	#~ nll = lambda *args: -lnlike(*args)
	#~ result = op.minimize_scalar(nll, bounds=(1,1000), method='bounded',args=(k, b ,errb ))
	#~ b1_ml = result["x"]
	nll = lambda *args: -lnlike(*args)
	result = op.minimize(nll, [pop],  method='Nelder-Mead', args=(k, b ,errb ),  options={'maxfev': 2000} )
	b1_ml, b2_ml, b3_ml, b4_ml, sig_ml = result["x"]
	#~ print(result)
	
	#~ max_l = lnlike(result["x"], k, b, errb )
	#~ AIC = 2*4. - 2 * max_l
	#~ print 'maximum likelihood is '+str(max_l)
	#~ print 'AIC = '+str(AIC)
	
	#~ ndim, nwalkers = 1, 200
	#~ pos = [result["x"] + 1e-4*np.random.randn(ndim) for i in range(nwalkers)]



	#~ sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(k, b, errb))
	#~ sampler.run_mcmc(pos, 1000)
	#~ samples = sampler.chain[:, 200:, :].reshape((-1, ndim))

	
	#~ import corner
	#~ fig = corner.corner(samples, labels=["$b1$", "$b2$", "$bs$", "$b3nl$" ], truths=[b1_ml, b2_ml, bs_ml, b3nl])
	#~ fig.savefig("/home/david/triangle.png")
	

	#~ b1_mcmc = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]), zip(*np.percentile(samples, [16, 50, 84], axis=0)))
	
	#~ print b1_mcmc
	

	#~ plt.figure()
	#~ plt.title('sigma kaiser')
	#~ for i in xrange(0,nwalkers):
		#~ plt.plot(np.arange(1000), sampler.chain[i,:,0])
	#~ plt.show()
	
	end = time.time()
	print 'time is '+str((end - start))
	#~ return b1_mcmc
	return b1_ml, b2_ml, b3_ml, b4_ml, sig_ml
	
########################################################################
#### Scoccimaro
########################################################################
def coeffit_Scocci(j, fcc, kstop,Pmod_dd, Pmod_dt, Pmod_tt, lb1, k ,b ,errb):
	start = time.time()
	#~ lim = np.where(k < kstop)[0]
	lim = np.where((k < kstop)&(k > 1e-2))[0]
	
	def lnlike(theta, x, y, yerr):
		sigma = theta
		#~ kappa = x[lim]*sigma
		#~ coeffA = math.sqrt(math.pi)/2. * erf(kappa)/kappa
		#~ coeffB = 3./2./kappa**2*(coeffA - np.exp(-kappa**2))
		#~ coeffC = 5./2./kappa**2*(coeffB - np.exp(-kappa**2))
		kappa = x[lim]*sigma*fcc[lim][ind]*Dz[ind]
		coeffA = np.arctan(kappa/math.sqrt(2))/(math.sqrt(2)*kappa) + 1/(2+kappa**2)
		coeffB = 6/kappa**2*(coeffA - 2/(2+kappa**2))
		coeffC = -10/kappa**2*(coeffB - 2/(2+kappa**2))
		model = lb1[lim]**2*Pmod_dd[lim]*coeffA + 2/3.*lb1[lim]*fcc[lim]*Pmod_dt[lim]*coeffB + 1/5.*fcc[lim]**2*Pmod_tt[lim]*coeffC
		inv_sigma2 = 1.0/(yerr[lim]**2)
		return -0.5*(np.sum((y[lim]-model)**2*inv_sigma2 - np.log(inv_sigma2)))
	
	def lnprior(theta):
		sigma = theta
		#~ if lb1 - 3*errlb1 < b1 < lb1 + 3*errlb1 :
		if 0 < sigma < 100:
			return 0.0
		return -np.inf
	
	def lnprob(theta, x, y, yerr):
		lp = lnprior(theta)
		if not np.isfinite(lp):
			return -np.inf
		return lp + lnlike(theta, x, y, yerr)
		
	z = [0.0,0.5,1.0,2.0]
	red = ['0.0','0.5','1.0','2.0','3.0']
	ind = red.index(str(z[j]))
	f = [0.518,0.754,0.872,0.956,0.98]
	Dz = [ 1.,0.77,0.61,0.42]

	nll = lambda *args: -lnlike(*args)
	#~ result = op.minimize(nll, [pop], bounds= [(1,1000)], args=(k, b ,errb ))
	result = op.minimize_scalar(nll, bounds=(1,1000), method='bounded',args=(k, b ,errb ))
	b1_ml = result["x"]
	print(result)
	
	#~ max_l = lnlike(result["x"], k, b, errb )
	#~ AIC = 2*4. - 2 * max_l
	#~ print 'maximum likelihood is '+str(max_l)
	#~ print 'AIC = '+str(AIC)
	
	#~ ndim, nwalkers = 1, 200
	#~ pos = [result["x"] + 1e-4*np.random.randn(ndim) for i in range(nwalkers)]



	#~ sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(k, b, errb))
	#~ sampler.run_mcmc(pos, 1000)
	
     

	#~ samples = sampler.chain[:, 200:, :].reshape((-1, ndim))

	
	#~ import corner
	#~ fig = corner.corner(samples, labels=["$b1$", "$b2$", "$bs$", "$b3nl$" ], truths=[b1_ml, b2_ml, bs_ml, b3nl])
	#~ fig.savefig("/home/david/triangle.png")
	

	#~ b1_mcmc = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]), zip(*np.percentile(samples, [16, 50, 84], axis=0)))
	
	#~ print b1_mcmc
	

	#~ plt.figure()
	#~ plt.title('sigma scocci')
	#~ for i in xrange(0,nwalkers):
		#~ plt.plot(np.arange(1000), sampler.chain[i,:,0])
	#~ plt.show()
	
	end = time.time()
	print 'time is '+str((end - start))
	#~ return b1_mcmc
	return b1_ml
	
########################################################################
#### Scoccimaro 2
########################################################################
def coeffit_Scocci2(j, fcc, kstop,Pmod_dd, Pmod_dt, Pmod_tt, k ,b ,errb):
	start = time.time()
	#~ lim = np.where(k < kstop)[0]
	lim = np.where((k < kstop)&(k > 1e-2))[0]
	
	def lnlike(theta, x, y, yerr):
		b1, b2, b3, b4, sigma = theta
		#~ kappa = x[lim]*sigma
		#~ coeffA = math.sqrt(math.pi)/2. * erf(kappa)/kappa
		#~ coeffB = 3./2./kappa**2*(coeffA - np.exp(-kappa**2))
		#~ coeffC = 5./2./kappa**2*(coeffB - np.exp(-kappa**2))
		bcc = b1 + b2*x[lim]**2 + b3*x[lim]**3 + b4*x[lim]**4 
		kappa = x[lim]*sigma*fcc[lim][ind]*Dz[ind]
		coeffA = np.arctan(kappa/math.sqrt(2))/(math.sqrt(2)*kappa) + 1/(2+kappa**2)
		coeffB = 6/kappa**2*(coeffA - 2/(2+kappa**2))
		coeffC = -10/kappa**2*(coeffB - 2/(2+kappa**2))
		model = bcc**2*Pmod_dd[lim]*coeffA + 2/3.*bcc*fcc[lim]*Pmod_dt[lim]*coeffB + 1/5.*fcc[lim]**2*Pmod_tt[lim]*coeffC
		inv_sigma2 = 1.0/(yerr[lim]**2)
		return -0.5*(np.sum((y[lim]-model)**2*inv_sigma2 - np.log(inv_sigma2)))
	
	def lnprior(theta):
		b1, b2, b3, b4, sigma = theta
		#~ if lb1 - 3*errlb1 < b1 < lb1 + 3*errlb1 :
		if 0 < sigma < 100:
			return 0.0
		return -np.inf
	
	def lnprob(theta, x, y, yerr):
		lp = lnprior(theta)
		if not np.isfinite(lp):
			return -np.inf
		return lp + lnlike(theta, x, y, yerr)
		
	z = [0.0,0.5,1.0,2.0]
	red = ['0.0','0.5','1.0','2.0','3.0']
	ind = red.index(str(z[j]))
	f = [0.518,0.754,0.872,0.956,0.98]
	Dz = [ 1.,0.77,0.61,0.42]

	pop = [1,1,1,1,6]
	#~ nll = lambda *args: -lnlike(*args)
	#~ result = op.minimize_scalar(nll, bounds=(1,1000), method='bounded',args=(k, b ,errb ))
	#~ b1_ml = result["x"]
	nll = lambda *args: -lnlike(*args)
	result = op.minimize(nll, [pop],  method='Nelder-Mead', args=(k, b ,errb ),  options={'maxfev': 2000} )
	b1_ml, b2_ml, b3_ml, b4_ml, sig_ml = result["x"]
	#~ print(result)
	
	#~ max_l = lnlike(result["x"], k, b, errb )
	#~ AIC = 2*4. - 2 * max_l
	#~ print 'maximum likelihood is '+str(max_l)
	#~ print 'AIC = '+str(AIC)
	
	#~ ndim, nwalkers = 1, 200
	#~ pos = [result["x"] + 1e-4*np.random.randn(ndim) for i in range(nwalkers)]



	#~ sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(k, b, errb))
	#~ sampler.run_mcmc(pos, 1000)
	
     

	#~ samples = sampler.chain[:, 200:, :].reshape((-1, ndim))

	
	#~ import corner
	#~ fig = corner.corner(samples, labels=["$b1$", "$b2$", "$bs$", "$b3nl$" ], truths=[b1_ml, b2_ml, bs_ml, b3nl])
	#~ fig.savefig("/home/david/triangle.png")
	

	#~ b1_mcmc = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]), zip(*np.percentile(samples, [16, 50, 84], axis=0)))
	
	#~ print b1_mcmc
	

	#~ plt.figure()
	#~ plt.title('sigma scocci')
	#~ for i in xrange(0,nwalkers):
		#~ plt.plot(np.arange(1000), sampler.chain[i,:,0])
	#~ plt.show()
	
	end = time.time()
	print 'time is '+str((end - start))
	#~ return b1_mcmc
	return b1_ml, b2_ml, b3_ml, b4_ml, sig_ml

########################################################################
#### TNS
########################################################################
def coeffit_TNS(j, fcc, kstop,Pmod_dd, Pmod_dt, Pmod_tt, lb1, k ,b ,errb,AB2,AB4,AB6,AB8 ):
		
	start = time.time()
	#~ lim = np.where(k < kstop)[0]
	lim = np.where((k < kstop)&(k > 1e-2))[0]
	
	def lnlike(theta, x, y, yerr):
		sigma = theta
		#~ kappa = x[lim]*sigma
		#~ coeffA = math.sqrt(math.pi)/2. * erf(kappa)/kappa
		#~ coeffB = 3./2./kappa**2*(coeffA - np.exp(-kappa**2))
		#~ coeffC = 5./2./kappa**2*(coeffB - np.exp(-kappa**2))
		#~ coeffD = 7./2./kappa**2*(coeffC - np.exp(-kappa**2))
		#~ coeffE = 9./2./kappa**2*(coeffD - np.exp(-kappa**2))
		kappa = x[lim]*sigma*fcc[ind]*Dz[ind]
		coeffA = np.arctan(kappa/math.sqrt(2))/(math.sqrt(2)*kappa) + 1/(2+kappa**2)
		coeffB = 6/kappa**2*(coeffA - 2/(2+kappa**2))
		coeffC = -10/kappa**2*(coeffB - 2/(2+kappa**2))
		coeffD = -2/3./kappa**2*(coeffC - 2/(2+kappa**2))
		coeffE = -4/10./kappa**2*(7.*coeffD - 2/(2+kappa**2))
		model = lb1[lim]**2*Pmod_dd[lim]*coeffA + 2/3.*lb1[lim]*fcc[lim]*Pmod_dt[lim]*coeffB + 1/5.*fcc[lim]**2*Pmod_tt[lim]*coeffC \
		+ (1/3.*AB2[lim]*coeffB+ 1/5.*AB4[lim]*coeffC+ 1/7.*AB6[lim]*coeffD+ 1/9.*AB8[lim]*coeffE)
		inv_sigma2 = 1.0/(yerr[lim]**2)
		return -0.5*(np.sum((y[lim]-model)**2*inv_sigma2 - np.log(inv_sigma2)))
	def lnprior(theta):
		sigma = theta
		#~ if lb1 - 3*errlb1 < b1 < lb1 + 3*errlb1 :
		if 0 < sigma < 100:
			return 0.0
		return -np.inf
	
	def lnprob(theta, x, y, yerr):
		lp = lnprior(theta)
		if not np.isfinite(lp):
			return -np.inf
		return lp + lnlike(theta, x, y, yerr)
		
	z = [0.0,0.5,1.0,2.0]
	red = ['0.0','0.5','1.0','2.0','3.0']
	ind = red.index(str(z[j]))
	f = [0.518,0.754,0.872,0.956,0.98]
	Dz = [ 1.,0.77,0.61,0.42]

	nll = lambda *args: -lnlike(*args)
	#~ result = op.minimize(nll, [pop], bounds= [(1,1000)], args=(k, b ,errb ))
	result = op.minimize_scalar(nll, bounds=(1,1000), method='bounded',args=(k, b ,errb ))
	b1_ml = result["x"]
	print(result)
	
	#~ max_l = lnlike(result["x"], k, b, errb )
	#~ AIC = 2*4. - 2 * max_l
	#~ print 'maximum likelihood is '+str(max_l)
	#~ print 'AIC = '+str(AIC)
	
	#~ ndim, nwalkers = 1, 200
	#~ pos = [result["x"] + 1e-4*np.random.randn(ndim) for i in range(nwalkers)]



	#~ sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(k, b, errb))
	#~ sampler.run_mcmc(pos, 1000)
	
     

	#~ samples = sampler.chain[:, 200:, :].reshape((-1, ndim))

	
	#~ import corner
	#~ fig = corner.corner(samples, labels=["$b1$", "$b2$", "$bs$", "$b3nl$" ], truths=[b1_ml, b2_ml, bs_ml, b3nl])
	#~ fig.savefig("/home/david/triangle.png")
	

	#~ b1_mcmc = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]), zip(*np.percentile(samples, [16, 50, 84], axis=0)))
	
	#~ print b1_mcmc
	

	#~ plt.figure()
	#~ plt.title('sigma tns')
	#~ for i in xrange(0,nwalkers):
		#~ plt.plot(np.arange(1000), sampler.chain[i,:,0])
	#~ plt.show()
	
	end = time.time()
	print 'time is '+str((end - start))
	#~ return b1_mcmc
	return b1_ml
	
########################################################################
#### TNS 2
########################################################################
def coeffit_TNS2(j, fcc, kstop,Pmod_dd, Pmod_dt, Pmod_tt, k ,b ,errb,AB2,AB4,AB6,AB8 ):
		
	start = time.time()
	#~ lim = np.where(k < kstop)[0]
	lim = np.where((k < kstop)&(k > 1e-2))[0]
	
	def lnlike(theta, x, y, yerr):
		b1, b2, b3, b4, sigma = theta
		#~ kappa = x[lim]*sigma
		#~ coeffA = math.sqrt(math.pi)/2. * erf(kappa)/kappa
		#~ coeffB = 3./2./kappa**2*(coeffA - np.exp(-kappa**2))
		#~ coeffC = 5./2./kappa**2*(coeffB - np.exp(-kappa**2))
		bcc = b1 + b2*x[lim]**2 + b3*x[lim]**3 + b4*x[lim]**4 
		kappa = x[lim]*sigma*fcc[ind]*Dz[ind]
		coeffA = np.arctan(kappa/math.sqrt(2))/(math.sqrt(2)*kappa) + 1/(2+kappa**2)
		coeffB = 6/kappa**2*(coeffA - 2/(2+kappa**2))
		coeffC = -10/kappa**2*(coeffB - 2/(2+kappa**2))
		coeffD = -2/3./kappa**2*(coeffC - 2/(2+kappa**2))
		coeffE = -4/10./kappa**2*(7.*coeffD - 2/(2+kappa**2))
		model = bcc**2*Pmod_dd[lim]*coeffA + 2/3.*bcc*fcc[lim]*Pmod_dt[lim]*coeffB + 1/5.*fcc[lim]**2*Pmod_tt[lim]*coeffC \
		+ (1/3.*AB2[lim]*coeffB+ 1/5.*AB4[lim]*coeffC+ 1/7.*AB6[lim]*coeffD+ 1/9.*AB8[lim]*coeffE)
		inv_sigma2 = 1.0/(yerr[lim]**2)
		return -0.5*(np.sum((y[lim]-model)**2*inv_sigma2 - np.log(inv_sigma2)))
	def lnprior(theta):
		b1, b2, b3, b4, sigma = theta
		#~ if lb1 - 3*errlb1 < b1 < lb1 + 3*errlb1 :
		if 0 < sigma < 100:
			return 0.0
		return -np.inf
	
	def lnprob(theta, x, y, yerr):
		lp = lnprior(theta)
		if not np.isfinite(lp):
			return -np.inf
		return lp + lnlike(theta, x, y, yerr)
		
	z = [0.0,0.5,1.0,2.0]
	red = ['0.0','0.5','1.0','2.0','3.0']
	ind = red.index(str(z[j]))
	f = [0.518,0.754,0.872,0.956,0.98]
	Dz = [ 1.,0.77,0.61,0.42]

	pop = [1,1,1,1,6]
	#~ nll = lambda *args: -lnlike(*args)
	#~ result = op.minimize_scalar(nll, bounds=(1,1000), method='bounded',args=(k, b ,errb ))
	#~ b1_ml = result["x"]
	nll = lambda *args: -lnlike(*args)
	result = op.minimize(nll, [pop],  method='Nelder-Mead', args=(k, b ,errb ),  options={'maxfev': 2000} )
	b1_ml, b2_ml, b3_ml, b4_ml, sig_ml = result["x"]
	#~ print(result)
	
	#~ max_l = lnlike(result["x"], k, b, errb )
	#~ AIC = 2*4. - 2 * max_l
	#~ print 'maximum likelihood is '+str(max_l)
	#~ print 'AIC = '+str(AIC)
	
	#~ ndim, nwalkers = 1, 200
	#~ pos = [result["x"] + 1e-4*np.random.randn(ndim) for i in range(nwalkers)]



	#~ sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(k, b, errb))
	#~ sampler.run_mcmc(pos, 1000)
	
     

	#~ samples = sampler.chain[:, 200:, :].reshape((-1, ndim))

	
	#~ import corner
	#~ fig = corner.corner(samples, labels=["$b1$", "$b2$", "$bs$", "$b3nl$" ], truths=[b1_ml, b2_ml, bs_ml, b3nl])
	#~ fig.savefig("/home/david/triangle.png")
	

	#~ b1_mcmc = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]), zip(*np.percentile(samples, [16, 50, 84], axis=0)))
	
	#~ print b1_mcmc
	

	#~ plt.figure()
	#~ plt.title('sigma tns')
	#~ for i in xrange(0,nwalkers):
		#~ plt.plot(np.arange(1000), sampler.chain[i,:,0])
	#~ plt.show()
	
	end = time.time()
	print 'time is '+str((end - start))
	#~ return b1_mcmc
	return b1_ml, b2_ml, b3_ml, b4_ml, sig_ml
	
########################################################################
#### TNS + saito or eTns
########################################################################
def coeffit_eTNS(j, fcc, kstop, b1, b2, bs, b3nl, Pmod_dd, Pmod_dt, Pmod_tt, A, B, C, D, E, F, G, H, k ,b ,\
	errb,AB2,AB4,AB6,AB8,N, sca=None ):
	
	start = time.time()
	#~ lim = np.where(k < kstop)[0]
	lim = np.where((k < kstop)&(k > 1e-2))[0]
	
	def lnlike(theta, x, y, yerr):
		#~ b1, b2, bs, b3nl, sigma = theta
		sigma = theta
		PsptD1z = b1**2*Pmod_dd[lim] + b1*b2*A[lim] + 1/4.*b2**2*B[lim] + b1*bs*C[lim] + 1/2.*b2*bs*D[lim] + 1/4.*bs**2*E[lim] \
		+ 2*b1*b3nl*F[lim] +N
		PsptT = b1* Pmod_dt[lim] + b2*G[lim] + bs*H[lim] + b3nl*F[lim]
		kappa = x[lim]*sigma
		#~ coeffA = math.sqrt(math.pi)/2. * erf(kappa)/kappa
		#~ coeffB = 3./2./kappa**2*(coeffA - np.exp(-kappa**2))
		#~ coeffC = 5./2./kappa**2*(coeffB - np.exp(-kappa**2))
		#~ coeffD = 7./2./kappa**2*(coeffC - np.exp(-kappa**2))
		#~ coeffE = 9./2./kappa**2*(coeffD - np.exp(-kappa**2))
		kappa = x[lim]*sigma*fcc[lim][ind]*Dz[ind]
		coeffA = np.arctan(kappa/math.sqrt(2))/(math.sqrt(2)*kappa) + 1/(2+kappa**2)
		coeffB = 6/kappa**2*(coeffA - 2/(2+kappa**2))
		coeffC = -10/kappa**2*(coeffB - 2/(2+kappa**2))
		coeffD = -2/3./kappa**2*(coeffC - 2/(2+kappa**2))
		coeffE = -4/10./kappa**2*(7.*coeffD - 2/(2+kappa**2))
		if sca:
			model = PsptD1z*coeffA*sca**2 + 2/3.*fcc[lim]*PsptT*coeffB*sca + 1/5.*fcc[lim]**2*Pmod_tt[lim]*coeffC \
			+ (1/3.*AB2[lim]*coeffB+ 1/5.*AB4[lim]*coeffC+ 1/7.*AB6[lim]*coeffD+ 1/9.*AB8[lim]*coeffE)
			inv_sigma2 = 1.0/(yerr[lim]**2)
		else:
			model = PsptD1z*coeffA + 2/3.*fcc[lim]*PsptT*coeffB + 1/5.*fcc[lim]**2*Pmod_tt[lim]*coeffC \
			+ (1/3.*AB2[lim]*coeffB+ 1/5.*AB4[lim]*coeffC+ 1/7.*AB6[lim]*coeffD+ 1/9.*AB8[lim]*coeffE)
			inv_sigma2 = 1.0/(yerr[lim]**2)
		return -0.5*(np.sum((y[lim]-model)**2*inv_sigma2 - np.log(inv_sigma2)))
		
	def lnprior(theta):
		sigma = theta
		#~ if lb1 - 3*errlb1 < b1 < lb1 + 3*errlb1 :
		if 0 < sigma < 100:
			return 0.0
		return -np.inf
	
	def lnprob(theta, x, y, yerr):
		lp = lnprior(theta)
		if not np.isfinite(lp):
			return -np.inf
		return lp + lnlike(theta, x, y, yerr)
		
	z = [0.0,0.5,1.0,2.0]
	red = ['0.0','0.5','1.0','2.0','3.0']
	ind = red.index(str(z[j]))
	f = [0.518,0.754,0.872,0.956,0.98]
	Dz = [ 1.,0.77,0.61,0.42]

	nll = lambda *args: -lnlike(*args)
	#~ result = op.minimize(nll, [pop], bounds= [(1,1000)], args=(k, b ,errb ))
	result = op.minimize_scalar(nll, bounds=(1,1000), method='bounded',args=(k, b ,errb ))
	b1_ml = result["x"]
	#~ print(result)
	
	#~ max_l = lnlike(result["x"], k, b, errb )
	#~ AIC = 2*4. - 2 * max_l
	#~ print 'maximum likelihood is '+str(max_l)
	#~ print 'AIC = '+str(AIC)
	
	#~ ndim, nwalkers = 1, 200
	#~ pos = [result["x"] + 1e-4*np.random.randn(ndim) for i in range(nwalkers)]



	#~ sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(k, b, errb))
	#~ sampler.run_mcmc(pos, 1000)
	
     

	#~ samples = sampler.chain[:, 200:, :].reshape((-1, ndim))

	
	#~ import corner
	#~ fig = corner.corner(samples, labels=["$b1$", "$b2$", "$bs$", "$b3nl$" ], truths=[b1_ml, b2_ml, bs_ml, b3nl])
	#~ fig.savefig("/home/david/triangle.png")
	

	#~ b1_mcmc = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]), zip(*np.percentile(samples, [16, 50, 84], axis=0)))
	
	#~ print b1_mcmc
	

	#~ plt.figure()
	#~ plt.title('sigma tns')
	#~ for i in xrange(0,nwalkers):
		#~ plt.plot(np.arange(1000), sampler.chain[i,:,0])
	#~ plt.show()
	end = time.time()
	print 'time is '+str((end - start))
	#~ return b1_mcmc
	return b1_ml
	
########################################################################
#### TNS + saito or eTns 2
########################################################################
def coeffit_eTNS2(j, fcc, kstop, Pmod_dd, Pmod_dt, Pmod_tt, A, B, C, D, E, F, G, H, k ,b ,\
	errb,AB2,AB4,AB6,AB8,N, sca=None ):
	
	start = time.time()
	#~ lim = np.where(k < kstop)[0]
	lim = np.where((k < kstop)&(k > 1e-2))[0]
	
	def lnlike(theta, x, y, yerr):
		b1, b2, N, sigma = theta
		#~ sigma = theta
		bs = -4/7.*(b1-1)
		b3nl = 32/315.*(b1-1)
		PsptD1z = b1**2*Pmod_dd[lim] + b1*b2*A[lim] + 1/4.*b2**2*B[lim] + b1*bs*C[lim] + 1/2.*b2*bs*D[lim] + 1/4.*bs**2*E[lim] \
		+ 2*b1*b3nl*F[lim] +N
		PsptT = b1* Pmod_dt[lim] + b2*G[lim] + bs*H[lim] + b3nl*F[lim]
		kappa = x[lim]*sigma
		#~ coeffA = math.sqrt(math.pi)/2. * erf(kappa)/kappa
		#~ coeffB = 3./2./kappa**2*(coeffA - np.exp(-kappa**2))
		#~ coeffC = 5./2./kappa**2*(coeffB - np.exp(-kappa**2))
		#~ coeffD = 7./2./kappa**2*(coeffC - np.exp(-kappa**2))
		#~ coeffE = 9./2./kappa**2*(coeffD - np.exp(-kappa**2))
		kappa = x[lim]*sigma*fcc[lim][ind]*Dz[ind]
		coeffA = np.arctan(kappa/math.sqrt(2))/(math.sqrt(2)*kappa) + 1/(2+kappa**2)
		coeffB = 6/kappa**2*(coeffA - 2/(2+kappa**2))
		coeffC = -10/kappa**2*(coeffB - 2/(2+kappa**2))
		coeffD = -2/3./kappa**2*(coeffC - 2/(2+kappa**2))
		coeffE = -4/10./kappa**2*(7.*coeffD - 2/(2+kappa**2))
		if sca:
			model = PsptD1z*coeffA*sca**2 + 2/3.*fcc[lim]*PsptT*coeffB*sca + 1/5.*fcc[lim]**2*Pmod_tt[lim]*coeffC \
			+ (1/3.*AB2[lim]*coeffB+ 1/5.*AB4[lim]*coeffC+ 1/7.*AB6[lim]*coeffD+ 1/9.*AB8[lim]*coeffE)
			inv_sigma2 = 1.0/(yerr[lim]**2)
		else:
			model = PsptD1z*coeffA + 2/3.*fcc[lim]*PsptT*coeffB + 1/5.*fcc[lim]**2*Pmod_tt[lim]*coeffC \
			+ (1/3.*AB2[lim]*coeffB+ 1/5.*AB4[lim]*coeffC+ 1/7.*AB6[lim]*coeffD+ 1/9.*AB8[lim]*coeffE)
			inv_sigma2 = 1.0/(yerr[lim]**2)
		return -0.5*(np.sum((y[lim]-model)**2*inv_sigma2 - np.log(inv_sigma2)))
		
	def lnprior(theta):
		b1, b2,N, sigma = theta
		#~ if lb1 - 3*errlb1 < b1 < lb1 + 3*errlb1 :
		if 0 < sigma < 100:
			return 0.0
		return -np.inf
	
	def lnprob(theta, x, y, yerr):
		lp = lnprior(theta)
		if not np.isfinite(lp):
			return -np.inf
		return lp + lnlike(theta, x, y, yerr)
		
	z = [0.0,0.5,1.0,2.0]
	red = ['0.0','0.5','1.0','2.0','3.0']
	ind = red.index(str(z[j]))
	f = [0.518,0.754,0.872,0.956,0.98]
	Dz = [ 1.,0.77,0.61,0.42]

	pop = [1,1,1,6]
	nll = lambda *args: -lnlike(*args)
	result = op.minimize(nll, [pop],  method='Nelder-Mead', args=(k, b ,errb ),  options={'maxfev': 5000} )
	b1_ml, b2_ml,N_ml, sigma_ml  = result["x"]
	print(result)
	
	#~ max_l = lnlike(result["x"], k, b, errb )
	#~ AIC = 2*4. - 2 * max_l
	#~ print 'maximum likelihood is '+str(max_l)
	#~ print 'AIC = '+str(AIC)
	
	#~ ndim, nwalkers = 1, 200
	#~ pos = [result["x"] + 1e-4*np.random.randn(ndim) for i in range(nwalkers)]



	#~ sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(k, b, errb))
	#~ sampler.run_mcmc(pos, 1000)
	
     

	#~ samples = sampler.chain[:, 200:, :].reshape((-1, ndim))

	
	#~ import corner
	#~ fig = corner.corner(samples, labels=["$b1$", "$b2$", "$bs$", "$b3nl$" ], truths=[b1_ml, b2_ml, bs_ml, b3nl])
	#~ fig.savefig("/home/david/triangle.png")
	

	#~ b1_mcmc = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]), zip(*np.percentile(samples, [16, 50, 84], axis=0)))
	
	#~ print b1_mcmc
	

	#~ plt.figure()
	#~ plt.title('sigma tns')
	#~ for i in xrange(0,nwalkers):
		#~ plt.plot(np.arange(1000), sampler.chain[i,:,0])
	#~ plt.show()
	end = time.time()
	print 'time is '+str((end - start))
	#~ return b1_mcmc
	return b1_ml, b2_ml,N_ml, sigma_ml 
