import matplotlib.pyplot as plt
import numpy as np
import pylab as pl
from crpropa import *
from auger_data import *
import os

#Functions to compute quantities
def f_cut(e, zr):
	if (e < zr or zr == 0):
		return 1.
	else:
		return np.exp(1.-e/zr)

def delta_n(e, frac, gamma, zeta, Rcut):
	return frac*e**(-gamma)*f_cut(e*(10**18.),zeta*Rcut)


def not_empty_file(filename):
	with open(filename, 'r') as readobj:
		one_char = readobj.read(1)
		if not one_char:
			return False
	return True

def output_dir():
	if not(os.path.exists("/output")):
		os.mkdir("/output")

# Simulation for Power Law and distance distribution
def simulate_dd_1D_pl(params,distances,num = 1000000,title = "new_simulation",model = "Gilmore12",sim_seed = 0):

	output_dir()

	#Energy distributions for each particle: H, He, N, Si, Fe
	bins = 10000
	energies = np.linspace(1.0,1000.0,bins)
	g = params[0]
	f_h = params[1]
	f_he = params[2]
	f_n = params[3]
	f_si = params[4]
	f_fe = params[5]
	n_d = len(distances[0])	#for-iteration number for different distances
	
	dist = distances[0]  # distance in Mpc
	probd = distances[1]
	
	
	sim = ModuleList()
	if sim_seed:
		Random_seedThreads(sim_seed)

	#simulation steps
	sim.add(SimplePropagation(1*kpc, 10*Mpc))


	#simulation processes
	sim.add(PhotoPionProduction(CMB()))
	sim.add(ElectronPairProduction(CMB()))
	sim.add(PhotoDisintegration(CMB()))
	
	if (model == "Gilmore12"):
		sim.add(PhotoPionProduction(IRB_Gilmore12()))
		sim.add(ElectronPairProduction(IRB_Gilmore12()))
		sim.add(PhotoDisintegration(IRB_Gilmore12()))
	
	if (model == "Dominguez11"): 
		sim.add(PhotoPionProduction(IRB_Dominguez11()))
		sim.add(ElectronPairProduction(IRB_Dominguez11()))
		sim.add(PhotoDisintegration(IRB_Dominguez11()))
	
	sim.add(NuclearDecay())

	#stop if particle reaches this energy level or below
	sim.add(MinimumEnergy(1*EeV))

	#define observer
	obs = Observer()
	#observer at x=0
	obs.add(ObserverPoint()) 

	# name of the file
	filename = 'output/'+title+'.dat'

	#write output on detection of event
	output = TextOutput(filename, Output.Event1D)
	#output.disableAll()
	#output.enable(Output.CurrentEnergyColumn)
	#output.enable(Output.SourceEnergyColumn)
	#output.enable(Output.CreatedPositionColumn)
	#output.enable(Output.SourcePositionColumn)
	#output.enable(Output.CurrentPositionColumn)
	output.enable(Output.SerialNumberColumn)
	#output.enable(Output.SourceIdColumn)
	#output.enable(Output.CurrentIdColumn)
	obs.onDetection(output)

	sim.add(obs)

	sourcelist = SourceList()

	#composition
	composition= SourceMultipleParticleTypes()
	
	composition.add(nucleusId(1,1),f_h)
	composition.add(nucleusId(4,2),f_he)
	composition.add(nucleusId(14,7),f_n)
	composition.add(nucleusId(28,14),f_si)
	composition.add(nucleusId(56,26),f_fe)

	#Load different sources (for different distances and energies)
	for k in range(n_d):
		source = Source()
		source.add(composition)
		source.add(SourceRedshift1D())
		source.add(SourcePosition(Vector3d(dist[k], 0, 0) * Mpc))
		source.add(SourcePowerLawSpectrum(1*EeV,1000.*EeV,-g))
		sourcelist.add(source, probd[k])
	
		
	# run simulation
	sim.setShowProgress(True)

	sim.run(sourcelist,num)

	# load events
	output.close()

def simulate_dd_1D_separate(gamma,distances,num = 1000000,title = "new_simulation",model = "Gilmore12",sim_seed = 0):
	
	output_dir()
	#Energy distributions for each particle: H, He, N, Si, Fe
	bins = 10000
	g = gamma
	n_d = len(distances[0])	#for-iteration number for different distances
	
	dist = distances[0]  # distance in Mpc
	probd = distances[1]
	
	elems = ['H','He','N','Si','Fe']
	
	for i in range(len(elems)):
				
		sim = ModuleList()
		if sim_seed:
			Random_seedThreads(sim_seed)

		#simulation steps
		sim.add(SimplePropagation(1*kpc, 10*Mpc))


		#simulation processes
		sim.add(PhotoPionProduction(CMB()))
		sim.add(ElectronPairProduction(CMB()))
		sim.add(PhotoDisintegration(CMB()))
		
		if (model == "Gilmore12"):
			sim.add(PhotoPionProduction(IRB_Gilmore12()))
			sim.add(ElectronPairProduction(IRB_Gilmore12()))
			sim.add(PhotoDisintegration(IRB_Gilmore12()))
		
		if (model == "Dominguez11"): 
			sim.add(PhotoPionProduction(IRB_Dominguez11()))
			sim.add(ElectronPairProduction(IRB_Dominguez11()))
			sim.add(PhotoDisintegration(IRB_Dominguez11()))
		
		sim.add(NuclearDecay())

		#stop if particle reaches this energy level or below
		sim.add(MinimumEnergy(1*EeV))

		#define observer
		obs = Observer()
		#observer at x=0
		obs.add(ObserverPoint()) 

		# name of the file
		filename = 'output/'+elems[i]+'_'+title+'.dat'

		#write output on detection of event
		output = TextOutput(filename, Output.Event1D)
		#output.disableAll()
		#output.enable(Output.CurrentEnergyColumn)
		#output.enable(Output.SourceEnergyColumn)
		#output.enable(Output.CreatedPositionColumn)
		#output.enable(Output.SourcePositionColumn)
		#output.enable(Output.CurrentPositionColumn)
		output.enable(Output.SerialNumberColumn)
		#output.enable(Output.SourceIdColumn)
		#output.enable(Output.CurrentIdColumn)
		obs.onDetection(output)

		sim.add(obs)

		sourcelist = SourceList()

		#composition
		composition= SourceMultipleParticleTypes()
		
		if elems[i]=='H':
			composition.add(nucleusId(1,1),1.0)
		elif elems[i]=='He':
			composition.add(nucleusId(4,2),1.0)
		elif elems[i]=='N':
			composition.add(nucleusId(14,7),1.0)
		elif elems[i]=='Si':
			composition.add(nucleusId(28,14),1.0)
		elif elems[i]=='Fe':
			composition.add(nucleusId(56,26),1.0)

		#Load different sources (for different distances and energies)
		for k in range(n_d):
			source = Source()
			source.add(composition)
			source.add(SourceRedshift1D())
			source.add(SourcePosition(Vector3d(dist[k], 0, 0) * Mpc))
			source.add(SourcePowerLawSpectrum(1*EeV,1000.*EeV,-g))
			sourcelist.add(source, probd[k])
		
			
		# run simulation
		sim.setShowProgress(True)

		sim.run(sourcelist,num)

		# load events
		output.close()
	

def simulate_dd_1D_parts(params,distances,num = 1000000,title = "new_simulation",model = "Gilmore12",parts = 1,n_e = 10000,seed = 0,sim_seed = 0):

	output_dir()

	#Energy distributions for each particle: H, He, N, Si, Fe
	bins = 10000
	energies = np.linspace(1.0,1000.0,bins)
	g = params[0]
	rcut = 10**params[1]
	f_h = params[2]
	f_he = params[3]
	f_n = params[4]
	f_si = params[5]
	f_fe = params[6]
	num_tot = num

	if seed:
		np.random.seed(seed)
	
	
	if(f_h != 0.):
		probs_h = np.array([delta_n(energies[j], f_h, g, 1., rcut) for j in range(bins)])
		n_h = probs_h.sum()
		probs_h /= n_h
		energy_h = np.random.choice(energies,n_e,p=probs_h)

	if(f_he != 0.):
		probs_he = np.array([delta_n(energies[j], f_he, g, 2., rcut) for j in range(bins)])
		n_he = probs_he.sum()
		probs_he /= n_he
		energy_he = np.random.choice(energies,n_e,p=probs_he)

	if(f_n != 0.):
		probs_n = np.array([delta_n(energies[j], f_n, g, 7., rcut) for j in range(bins)])
		n_n = probs_n.sum()
		probs_n /= n_n
		energy_n = np.random.choice(energies,n_e,p=probs_n)

	if(f_si != 0.):
		probs_si = np.array([delta_n(energies[j], f_si, g, 14., rcut) for j in range(bins)])
		n_si = probs_si.sum()
		probs_si /= n_si
		energy_si = np.random.choice(energies,n_e,p=probs_si)

	if(f_fe != 0.):
		probs_fe = np.array([delta_n(energies[j], f_fe, g, 26., rcut) for j in range(bins)])
		n_fe = probs_fe.sum()
		probs_fe /= n_fe
		energy_fe = np.random.choice(energies,n_e,p=probs_fe)
	
	dist_tot = distances[0]
	probd_tot = distances[1]
	num_parts = len(dist_tot) / parts

	for i_parts in range(parts):

		if (i_parts != (parts - 1)):
			dist = dist_tot[num_parts * i_parts:num_parts * (i_parts+1)]  # distance in Mpc
			probd = probd_tot[num_parts * i_parts:num_parts * (i_parts+1)]
			n_d = len(dist)	#for-iteration number for different distances
			num = int(num_tot * probd.sum())
		else:
			dist = dist_tot[num_parts * i_parts:]  # distance in Mpc
			probd = probd_tot[num_parts * i_parts:]
			n_d = len(dist)	#for-iteration number for different distances
			num = int(num_tot * probd.sum())
		
		sim = ModuleList()
		if sim_seed:
			Random_seedThreads(sim_seed)

		#simulation steps
		sim.add(SimplePropagation(1*kpc, 10*Mpc))


		#simulation processes
		sim.add(PhotoPionProduction(CMB()))
		sim.add(ElectronPairProduction(CMB()))
		sim.add(PhotoDisintegration(CMB()))
		
		if (model == "Gilmore12"):
			sim.add(PhotoPionProduction(IRB_Gilmore12()))
			sim.add(ElectronPairProduction(IRB_Gilmore12()))
			sim.add(PhotoDisintegration(IRB_Gilmore12()))
		
		if (model == "Dominguez11"): 
			sim.add(PhotoPionProduction(IRB_Dominguez11()))
			sim.add(ElectronPairProduction(IRB_Dominguez11()))
			sim.add(PhotoDisintegration(IRB_Dominguez11()))
		
		sim.add(NuclearDecay())

		#stop if particle reaches this energy level or below
		sim.add(MinimumEnergy(1*EeV))

		#define observer
		obs = Observer()
		#observer at x=0
		obs.add(ObserverPoint()) 

		# name of the file
		filename = 'output/'+title+str(i_parts)+'.dat'

		#write output on detection of event
		output = TextOutput(filename, Output.Event1D)
		#output.disableAll()
		#output.enable(Output.CurrentEnergyColumn)
		#output.enable(Output.SourceEnergyColumn)
		#output.enable(Output.CreatedPositionColumn)
		#output.enable(Output.SourcePositionColumn)
		#output.enable(Output.CurrentPositionColumn)
		output.enable(Output.SerialNumberColumn)
		#output.enable(Output.SourceIdColumn)
		#output.enable(Output.CurrentIdColumn)
		obs.onDetection(output)

		sim.add(obs)

		sourcelist = SourceList()

		#composition list
		composition_h= SourceMultipleParticleTypes()
		composition_h.add(nucleusId(1,1),1)

		composition_he= SourceMultipleParticleTypes()
		composition_he.add(nucleusId(4,2),1)

		composition_n= SourceMultipleParticleTypes()
		composition_n.add(nucleusId(14,7),1)

		composition_si= SourceMultipleParticleTypes()
		composition_si.add(nucleusId(28,14),1)

		composition_fe= SourceMultipleParticleTypes()
		composition_fe.add(nucleusId(56,26),1)

		#Load different sources (for different distances and energies)
		for k in range(n_d):
			for i in range(n_e):
				if(f_h != 0.):
					source_h = Source()
					source_h.add(SourcePosition(Vector3d(dist[k], 0, 0) * Mpc))
					source_h.add(SourceIsotropicEmission())
					source_h.add(composition_h)
					source_h.add( SourceEnergy(energy_h[i] * EeV) )
					source_h.add(SourceRedshift1D())
					sourcelist.add(source_h, f_h*probd[k])
				
				if(f_he != 0.):
					source_he = Source()
					source_he.add(SourcePosition(Vector3d(dist[k], 0, 0) * Mpc))
					source_he.add(SourceIsotropicEmission())
					source_he.add(composition_he)
					source_he.add( SourceEnergy(energy_he[i] * EeV) )
					source_he.add(SourceRedshift1D())
					sourcelist.add(source_he, f_he*probd[k])
				
				if(f_n != 0.):	
					source_n = Source()
					source_n.add(SourcePosition(Vector3d(dist[k], 0, 0) * Mpc))
					source_n.add(SourceIsotropicEmission())
					source_n.add(composition_n)
					source_n.add( SourceEnergy(energy_n[i] * EeV) )
					source_n.add(SourceRedshift1D())
					sourcelist.add(source_n, f_n*probd[k])
				
				if(f_si != 0.):	
					source_si = Source()
					source_si.add(SourcePosition(Vector3d(dist[k], 0, 0) * Mpc))
					source_si.add(SourceIsotropicEmission())
					source_si.add(composition_si)
					source_si.add( SourceEnergy(energy_si[i] * EeV) )
					source_si.add(SourceRedshift1D())
					sourcelist.add(source_si, f_si*probd[k])
					
				if(f_fe != 0.):	
					source_fe = Source()
					source_fe.add(SourcePosition(Vector3d(dist[k], 0, 0) * Mpc))
					source_fe.add(SourceIsotropicEmission())
					source_fe.add(composition_fe)
					source_fe.add( SourceEnergy(energy_fe[i] * EeV) )
					source_fe.add(SourceRedshift1D())
					sourcelist.add(source_fe, f_fe*probd[k])

		# run simulation
		sim.setShowProgress(True)

		sim.run(sourcelist,num)

		# load events
		output.close()

def plot_power(title = "new_simulation",plotfile = "new_simulation_plot",plottitle = "My simulation",power = 2., logemin = 18, logemax = 20.4):

	output_dir()
	
	mask = (ecens >= logemin) & (ecens <= logemax)
	mask_bins = (ebins >= logemin) & (ebins <= logemax)
	ebins_ = ebins[mask_bins]
	ecens_ = ecens[mask]
	auger_ = auger[mask]
	sauger_ = sauger[mask]
	nauger_ = nauger[mask]
	sigma_auger_ = sigma_auger[mask]

	# load events
	filename = 'output/'+title+'.dat'
	fileout = 'output/'+plotfile+'.png'
	d = pl.genfromtxt(filename, names=True)

	# observed quantities
	Z = pl.array([chargeNumber(id) for id in d['ID'].astype(int)])  # element
	A = pl.array([massNumber(id) for id in d['ID'].astype(int)])  # atomic mass number
	lE = pl.log10(d['E']) + 18  # energy in log10(E/eV))

	lEbins = ebins_  # logarithmic bins
	lEcens = (lEbins[1:] + lEbins[:-1]) / 2  # logarithmic bin centers
	dE = 10**lEbins[1:] - 10**lEbins[:-1]  # bin widths

	# identify mass groups
	idx1 = A == 1
	idx2 = (A > 1) * (A <= 4)
	idx3 = (A > 4) * (A <= 22)
	idx4 = (A > 22) * (A <= 38)
	idx5 = (A > 38)

	# calculate spectrum: J(E) = dN/dE
	J  = pl.histogram(lE, bins=lEbins)[0] / dE
	J1 = pl.histogram(lE[idx1], bins=lEbins)[0] / dE
	J2 = pl.histogram(lE[idx2], bins=lEbins)[0] / dE
	J3 = pl.histogram(lE[idx3], bins=lEbins)[0] / dE
	J4 = pl.histogram(lE[idx4], bins=lEbins)[0] / dE
	J5 = pl.histogram(lE[idx5], bins=lEbins)[0] / dE
	
	# Compute J * E^{power}
	J = J * (10**lEcens)**power
	J1 = J1 * (10**lEcens)**power
	J2 = J2 * (10**lEcens)**power
	J3 = J3 * (10**lEcens)**power
	J4 = J4 * (10**lEcens)**power
	J5 = J5 * (10**lEcens)**power
	
	
	pl.figure(figsize=(10,7))
	pl.plot(lEcens, J,  color='SaddleBrown', label='Total')
	pl.plot(lEcens, J1, color='blue', label='A = 1')
	pl.plot(lEcens, J2, color='grey', label='A = 2-4')
	pl.plot(lEcens, J3, color='green', label='A = 5-22')
	pl.plot(lEcens, J4, color='purple', label='A = 23-38')
	pl.plot(lEcens, J5, color='red', label='A $>$ 38')

	pl.legend(fontsize=15, frameon=True)
	pl.title(plottitle)
	pl.grid()
	pl.ylabel('$J(E) E^$'+str(power))
	pl.xlabel('$\log_{10}$(E/eV)')
	pl.savefig(fileout)


def plot_parts(title = "new_simulation",plotfile = "new_simulation",plottitle = "My simulation",parts = 1, logemin = 18, logemax = 20.4):

	output_dir()
	
	mask = (ecens >= logemin) & (ecens <= logemax)
	mask_bins = (ebins >= logemin) & (ebins <= logemax)
	ebins_ = ebins[mask_bins]
	ecens_ = ecens[mask]
	auger_ = auger[mask]
	sauger_ = sauger[mask]
	nauger_ = nauger[mask]
	sigma_auger_ = sigma_auger[mask]

	# load events
	filename = ['output/'+title+str(i)+'.dat' for i in range(parts)]
	fileout = 'output/'+plotfile+'.png'
	
	
	lEbins = ebins_  # logarithmic bins
	lEcens = (lEbins[1:] + lEbins[:-1]) / 2  # logarithmic bin centers
	dE = 10**lEbins[1:] - 10**lEbins[:-1]  # bin widths
	
	count = 0

	for i_parts in range(parts):
	
		if not_empty_file(filename[i_parts]):
			d = pl.genfromtxt(filename[i_parts], names=True)
			
			if (d.ndim == 2):
				# observed quantities
				Z = pl.array([chargeNumber(id) for id in d['ID'].astype(int)])  # element
				A = pl.array([massNumber(id) for id in d['ID'].astype(int)])  # atomic mass number
				
				lE = pl.log10(d['E']) + 18  # energy in log10(E/eV))
				
				# identify mass groups
				idx1 = A == 1
				idx2 = (A > 1) * (A <= 4)
				idx3 = (A > 4) * (A <= 22)
				idx4 = (A > 22) * (A <= 38)
				idx5 = (A > 38)

				
				
				if (count == 0):
					# calculate spectrum: J(E) = dN/dE
					J  = pl.histogram(lE, bins=lEbins)[0] / dE
					J1 = pl.histogram(lE[idx1], bins=lEbins)[0] / dE
					J2 = pl.histogram(lE[idx2], bins=lEbins)[0] / dE
					J3 = pl.histogram(lE[idx3], bins=lEbins)[0] / dE
					J4 = pl.histogram(lE[idx4], bins=lEbins)[0] / dE
					J5 = pl.histogram(lE[idx5], bins=lEbins)[0] / dE
				else:
					J  = J + pl.histogram(lE, bins=lEbins)[0] / dE
					J1 = J1 + pl.histogram(lE[idx1], bins=lEbins)[0] / dE
					J2 = J2 + pl.histogram(lE[idx2], bins=lEbins)[0] / dE
					J3 = J3 + pl.histogram(lE[idx3], bins=lEbins)[0] / dE
					J4 = J4 + pl.histogram(lE[idx4], bins=lEbins)[0] / dE
					J5 = J5 + pl.histogram(lE[idx5], bins=lEbins)[0] / dE
					
				count += 1

	
	# normalize
	J1 /= J[0]
	J2 /= J[0]
	J3 /= J[0]
	J4 /= J[0]
	J5 /= J[0]
	J /= J[0]

	# calculate spectrum: J(E) = dN/dE for AUGER
	Ja1 = plt.hist(ecens_,bins=ebins_,weights=auger1_)[0] / dE
	Ja2 = plt.hist(ecens_,bins=ebins_,weights=auger2_)[0] / dE
	Ja3 = plt.hist(ecens_,bins=ebins_,weights=auger3_)[0] / dE
	Ja4 = plt.hist(ecens_,bins=ebins_,weights=auger4_)[0] / dE
	Ja5 = plt.hist(ecens_,bins=ebins_,weights=auger5_)[0] / dE
	Ja = plt.hist(ecens_,bins=ebins_,weights=auger_)[0] / dE

	
	# normalize AUGER
	Ja1 /= Ja[0]
	Ja2 /= Ja[0]
	Ja3 /= Ja[0]
	Ja4 /= Ja[0]
	Ja5 /= Ja[0]
	Ja /= Ja[0]

	pl.figure(figsize=(10,7))
	pl.plot(lEcens, J,  color='SaddleBrown', label='Total')
	pl.plot(lEcens, J1, color='blue', label='A = 1')
	pl.plot(lEcens, J2, color='grey', label='A = 2-4')
	pl.plot(lEcens, J3, color='green', label='A = 5-22')
	pl.plot(lEcens, J4, color='purple', label='A = 23-38')
	pl.plot(lEcens, J5, color='red', label='A $>$ 38')

	pl.plot(lEcens, Ja,  color='SaddleBrown',linestyle='dashed')
	pl.plot(lEcens, Ja1, color='blue',linestyle='dashed')
	pl.plot(lEcens, Ja2, color='grey',linestyle='dashed')
	pl.plot(lEcens, Ja3, color='green',linestyle='dashed')
	pl.plot(lEcens, Ja4, color='purple',linestyle='dashed')
	pl.plot(lEcens, Ja5, color='red',linestyle='dashed')

	pl.legend(fontsize=15, frameon=True)
	pl.semilogy()
	pl.ylim(1e-5)
	pl.title(plottitle)
	pl.grid()
	pl.ylabel('$J(E)$ [a.u.]')
	pl.xlabel('$\log_{10}$(E/eV)')
	pl.savefig(fileout)



def plot_power_parts(title = "new_simulation",plotfile = "new_simulation",plottitle = "My simulation",parts = 1, power = 2., logemin = 18, logemax = 20.4):

	output_dir()
	
	mask = (ecens >= logemin) & (ecens <= logemax)
	mask_bins = (ebins >= logemin) & (ebins <= logemax)
	ebins_ = ebins[mask_bins]
	ecens_ = ecens[mask]
	auger_ = auger[mask]
	sauger_ = sauger[mask]
	nauger_ = nauger[mask]
	sigma_auger_ = sigma_auger[mask]

	# load events
	filename = ['output/'+title+str(i)+'.dat' for i in range(parts)]
	fileout = 'output/'+plotfile+'.png'
	
	
	lEbins = ebins_  # logarithmic bins
	lEcens = (lEbins[1:] + lEbins[:-1]) / 2  # logarithmic bin centers
	dE = 10**lEbins[1:] - 10**lEbins[:-1]  # bin widths
	
	count = 0

	for i_parts in range(parts):
	
		if not_empty_file(filename[i_parts]):
			d = pl.genfromtxt(filename[i_parts], names=True)
			
			if (d.ndim == 2):
				# observed quantities
				Z = pl.array([chargeNumber(id) for id in d['ID'].astype(int)])  # element
				A = pl.array([massNumber(id) for id in d['ID'].astype(int)])  # atomic mass number
				
				lE = pl.log10(d['E']) + 18  # energy in log10(E/eV))
				
				# identify mass groups
				idx1 = A == 1
				idx2 = (A > 1) * (A <= 4)
				idx3 = (A > 4) * (A <= 22)
				idx4 = (A > 22) * (A <= 38)
				idx5 = (A > 38)

				
				
				if (count == 0):
					# calculate spectrum: J(E) = dN/dE
					J  = pl.histogram(lE, bins=lEbins)[0] / dE
					J1 = pl.histogram(lE[idx1], bins=lEbins)[0] / dE
					J2 = pl.histogram(lE[idx2], bins=lEbins)[0] / dE
					J3 = pl.histogram(lE[idx3], bins=lEbins)[0] / dE
					J4 = pl.histogram(lE[idx4], bins=lEbins)[0] / dE
					J5 = pl.histogram(lE[idx5], bins=lEbins)[0] / dE
				else:
					J  = J + pl.histogram(lE, bins=lEbins)[0] / dE
					J1 = J1 + pl.histogram(lE[idx1], bins=lEbins)[0] / dE
					J2 = J2 + pl.histogram(lE[idx2], bins=lEbins)[0] / dE
					J3 = J3 + pl.histogram(lE[idx3], bins=lEbins)[0] / dE
					J4 = J4 + pl.histogram(lE[idx4], bins=lEbins)[0] / dE
					J5 = J5 + pl.histogram(lE[idx5], bins=lEbins)[0] / dE
				
				count += 1

	
	
	# Compute J * E^{power}
	J = J * (10**lEcens)**power
	J1 = J1 * (10**lEcens)**power
	J2 = J2 * (10**lEcens)**power
	J3 = J3 * (10**lEcens)**power
	J4 = J4 * (10**lEcens)**power
	J5 = J5 * (10**lEcens)**power
	
	pl.figure(figsize=(10,7))
	pl.plot(lEcens, J,  color='SaddleBrown', label='Total')
	pl.plot(lEcens, J1, color='blue', label='A = 1')
	pl.plot(lEcens, J2, color='grey', label='A = 2-4')
	pl.plot(lEcens, J3, color='green', label='A = 5-22')
	pl.plot(lEcens, J4, color='purple', label='A = 23-38')
	pl.plot(lEcens, J5, color='red', label='A $>$ 38')

	
	pl.legend(fontsize=15, frameon=True)
	pl.title(plottitle)
	pl.grid()
	pl.ylabel('$J(E) E^$'+str(power))
	pl.xlabel('$\log_{10}$(E/eV)')
	pl.savefig(fileout)


def plot_power_sources(title = "new_simulation",plotfile = "new_simulation",plottitle = "My simulation",parts = 1, power = 2., logemin = 18, logemax = 20.4):

	output_dir()
	
	mask = (ecens >= logemin) & (ecens <= logemax)
	mask_bins = (ebins >= logemin) & (ebins <= logemax)
	ebins_ = ebins[mask_bins]
	ecens_ = ecens[mask]
	auger_ = auger[mask]
	sauger_ = sauger[mask]
	nauger_ = nauger[mask]
	sigma_auger_ = sigma_auger[mask]

	# load events
	filename = ['output/'+title+str(i)+'.dat' for i in range(parts)]
	fileout = 'output/'+plotfile+'.png'
	
	
	lEbins = ebins_  # logarithmic bins
	lEcens = (lEbins[1:] + lEbins[:-1]) / 2  # logarithmic bin centers
	dE = 10**lEbins[1:] - 10**lEbins[:-1]  # bin widths
	
	count = 0

	for i_parts in range(parts):
	
		if not_empty_file(filename[i_parts]):
			d = pl.genfromtxt(filename[i_parts], names=True)
			db = pl.genfromtxt(filename[i_parts])
			if (db.ndim == 2):
				# observed quantities
				Z = pl.array([chargeNumber(id) for id in d['ID0'].astype(int)])  # element
				A = pl.array([massNumber(id) for id in d['ID0'].astype(int)])  # atomic mass number
				
				lE = pl.log10(d['E0']) + 18  # energy in log10(E/eV))
				
				# identify mass groups
				idx1 = A == 1
				idx2 = (A > 1) * (A <= 4)
				idx3 = (A > 4) * (A <= 22)
				idx4 = (A > 22) * (A <= 38)
				idx5 = (A > 38)

				
				
				if (count == 0):
					# calculate spectrum: J(E) = dN/dE
					J  = pl.histogram(lE, bins=lEbins)[0] / dE
					J1 = pl.histogram(lE[idx1], bins=lEbins)[0] / dE
					J2 = pl.histogram(lE[idx2], bins=lEbins)[0] / dE
					J3 = pl.histogram(lE[idx3], bins=lEbins)[0] / dE
					J4 = pl.histogram(lE[idx4], bins=lEbins)[0] / dE
					J5 = pl.histogram(lE[idx5], bins=lEbins)[0] / dE
				else:
					J  = J + pl.histogram(lE, bins=lEbins)[0] / dE
					J1 = J1 + pl.histogram(lE[idx1], bins=lEbins)[0] / dE
					J2 = J2 + pl.histogram(lE[idx2], bins=lEbins)[0] / dE
					J3 = J3 + pl.histogram(lE[idx3], bins=lEbins)[0] / dE
					J4 = J4 + pl.histogram(lE[idx4], bins=lEbins)[0] / dE
					J5 = J5 + pl.histogram(lE[idx5], bins=lEbins)[0] / dE
				
				count += 1

	
	
	# Compute J * E^{power}
	J = J * (10**lEcens)**power
	J1 = J1 * (10**lEcens)**power
	J2 = J2 * (10**lEcens)**power
	J3 = J3 * (10**lEcens)**power
	J4 = J4 * (10**lEcens)**power
	J5 = J5 * (10**lEcens)**power
	
	pl.figure(figsize=(10,7))
	pl.plot(lEcens, J,  color='SaddleBrown', label='Total')
	pl.plot(lEcens, J1, color='blue', label='A = 1')
	pl.plot(lEcens, J2, color='grey', label='A = 2-4')
	pl.plot(lEcens, J3, color='green', label='A = 5-22')
	pl.plot(lEcens, J4, color='purple', label='A = 23-38')
	pl.plot(lEcens, J5, color='red', label='A $>$ 38')

	
	pl.legend(fontsize=15, frameon=True)
	pl.title(plottitle)
	pl.grid()
	pl.ylabel('$J(E) E**$'+str(power))
	pl.xlabel('$\log_{10}$(E/eV)')
	pl.savefig(fileout)


def plot_errors_rcut(title = "new_simulation",plotfile = "new_simulation",plottitle = "My simulation",rcut = 21., logemin = 18, logemax = 20.4):

	output_dir()
	
	mask = (ecens >= logemin) & (ecens <= logemax)
	mask_bins = (ebins >= logemin) & (ebins <= logemax)
	ebins_ = ebins[mask_bins]
	ecens_ = ecens[mask]
	auger_ = auger[mask]
	sauger_ = sauger[mask]
	nauger_ = nauger[mask]
	sigma_auger_ = sigma_auger[mask]

	# load events
	filename = 'output/'+title+'.dat'
	fileout = 'output/'+plotfile+'.png'
	d = pl.genfromtxt(filename, names=True)

	# observed quantities
	Z = pl.array([chargeNumber(id) for id in d['ID'].astype(int)])  # element
	Z0 = pl.array([chargeNumber(id) for id in d['ID0'].astype(int)])  # element at source
	A = pl.array([massNumber(id) for id in d['ID'].astype(int)])  # atomic mass number
	lE = pl.log10(d['E']) + 18  # energy in log10(E/eV))
	num = len(Z0)
	energy0 = d['E0']

	lEbins = ebins_  # logarithmic bins
	lEcens = (lEbins[1:] + lEbins[:-1]) / 2  # logarithmic bin centers
	dE = 10**lEbins[1:] - 10**lEbins[:-1]  # bin widths

	# identify mass groups
	idx1 = A == 1
	idx2 = (A > 1) * (A <= 4)
	idx3 = (A > 4) * (A <= 22)
	idx4 = (A > 22) * (A <= 38)
	idx5 = (A > 38)
	
	# Modifying distribution according to rcut
	weight_rcut = np.array([f_cut(energy0[i]*(10**18.),Z0[i]*(10**rcut)) for i in range(num)])

	# calculate spectrum: J(E) = dN/dE
	J  = pl.histogram(lE, bins=lEbins,weights=weight_rcut)[0] / dE
	J1 = pl.histogram(lE[idx1], bins=lEbins, weights=weight_rcut[idx1])[0] / dE
	J2 = pl.histogram(lE[idx2], bins=lEbins, weights=weight_rcut[idx2])[0] / dE
	J3 = pl.histogram(lE[idx3], bins=lEbins, weights=weight_rcut[idx3])[0] / dE
	J4 = pl.histogram(lE[idx4], bins=lEbins, weights=weight_rcut[idx4])[0] / dE
	J5 = pl.histogram(lE[idx5], bins=lEbins, weights=weight_rcut[idx5])[0] / dE

	# calculate spectrum: J(E) = dN/dE for AUGER
	Ja = plt.hist(ecens_,bins=ebins_,weights=auger_)[0] / dE
	Ja = auger_ / dE
	Jerrors = sigma_auger_ / dE

	# normalize
	J1 /= J[0]
	J2 /= J[0]
	J3 /= J[0]
	J4 /= J[0]
	J5 /= J[0]
	J /= J[0]

	# normalize AUGER
	Jerrors /= Ja[0]
	Ja /= Ja[0]
	

	pl.figure(figsize=(10,7))
	pl.plot(lEcens, J,  color='SaddleBrown', label='Total')
	pl.plot(lEcens, J1, color='blue', label='A = 1')
	pl.plot(lEcens, J2, color='grey', label='A = 2-4')
	pl.plot(lEcens, J3, color='green', label='A = 5-22')
	pl.plot(lEcens, J4, color='purple', label='A = 23-38')
	pl.plot(lEcens, J5, color='red', label='A $>$ 38')

	pl.plot(lEcens, Ja, "ok")
	pl.errorbar(lEcens, Ja, yerr=Jerrors, fmt = "ok")

	pl.legend(fontsize=15, frameon=True)
	pl.semilogy()
	pl.ylim(1e-5)
	pl.title(plottitle)
	pl.grid()
	pl.ylabel('$J(E)$ [a.u.]')
	pl.xlabel('$\log_{10}$(E/eV)')
	pl.savefig(fileout)

def plot_errors_rcut_separate(fractions,title = "new_simulation",plotfile = "new_simulation",plottitle = "My simulation",rcut = 21., logemin = 18, logemax = 20.4):

	mask = (ecens >= logemin) & (ecens <= logemax)
	mask_bins = (ebins >= logemin) & (ebins <= logemax)
	ebins_ = ebins[mask_bins]
	ecens_ = ecens[mask]
	auger_ = auger[mask]
	sauger_ = sauger[mask]
	nauger_ = nauger[mask]
	sigma_auger_ = sigma_auger[mask]
	
	fileout = 'output/'+plotfile+'.png'
	first = True
	elems = ['H','He','N','Si','Fe']
	
	for i in range(len(elems)):
		
		if round(fractions[i],7)==0:
			continue
		
		# load events
		filename = 'output/'+elems[i]+'_'+title+'.dat'
		d = pl.genfromtxt(filename, names=True)

		# observed quantities
		Z = pl.array([chargeNumber(id) for id in d['ID'].astype(int)])  # element
		Z0 = pl.array([chargeNumber(id) for id in d['ID0'].astype(int)])  # element at source
		A = pl.array([massNumber(id) for id in d['ID'].astype(int)])  # atomic mass number
		lE = pl.log10(d['E']) + 18  # energy in log10(E/eV))
		num = len(Z0)
		energy0 = d['E0']

		lEbins = ebins_  # logarithmic bins
		lEcens = (lEbins[1:] + lEbins[:-1]) / 2  # logarithmic bin centers
		dE = 10**lEbins[1:] - 10**lEbins[:-1]  # bin widths

		# identify mass groups
		idx1 = A == 1
		idx2 = (A > 1) * (A <= 4)
		idx3 = (A > 4) * (A <= 22)
		idx4 = (A > 22) * (A <= 38)
		idx5 = (A > 38)
		
		# Modifying distribution according to rcut
		weight_rcut = np.array([f_cut(energy0[w]*(10**18.),Z0[w]*(10**rcut)) for w in range(num)])
		weight_rcut *= fractions[i]

		# calculate spectrum: J(E) = dN/dE
		if first:
			J  = pl.histogram(lE, bins=lEbins,weights=weight_rcut)[0] / dE
			J1 = pl.histogram(lE[idx1], bins=lEbins, weights=weight_rcut[idx1])[0] / dE
			J2 = pl.histogram(lE[idx2], bins=lEbins, weights=weight_rcut[idx2])[0] / dE
			J3 = pl.histogram(lE[idx3], bins=lEbins, weights=weight_rcut[idx3])[0] / dE
			J4 = pl.histogram(lE[idx4], bins=lEbins, weights=weight_rcut[idx4])[0] / dE
			J5 = pl.histogram(lE[idx5], bins=lEbins, weights=weight_rcut[idx5])[0] / dE
			first = False
		else:
			J  += pl.histogram(lE, bins=lEbins,weights=weight_rcut)[0] / dE
			J1 += pl.histogram(lE[idx1], bins=lEbins, weights=weight_rcut[idx1])[0] / dE
			J2 += pl.histogram(lE[idx2], bins=lEbins, weights=weight_rcut[idx2])[0] / dE
			J3 += pl.histogram(lE[idx3], bins=lEbins, weights=weight_rcut[idx3])[0] / dE
			J4 += pl.histogram(lE[idx4], bins=lEbins, weights=weight_rcut[idx4])[0] / dE
			J5 += pl.histogram(lE[idx5], bins=lEbins, weights=weight_rcut[idx5])[0] / dE

	# calculate spectrum: J(E) = dN/dE for AUGER
	Ja = plt.hist(ecens_,bins=ebins_,weights=auger_)[0] / dE
	Ja = auger_ / dE
	Jerrors = sigma_auger_ / dE

	# normalize
	J1 /= J[0]
	J2 /= J[0]
	J3 /= J[0]
	J4 /= J[0]
	J5 /= J[0]
	J /= J[0]

	# normalize AUGER
	Jerrors /= Ja[0]
	Ja /= Ja[0]
	

	pl.figure(figsize=(10,7))
	pl.plot(lEcens, J,  color='SaddleBrown', label='Total')
	pl.plot(lEcens, J1, color='blue', label='A = 1')
	pl.plot(lEcens, J2, color='grey', label='A = 2-4')
	pl.plot(lEcens, J3, color='green', label='A = 5-22')
	pl.plot(lEcens, J4, color='purple', label='A = 23-38')
	pl.plot(lEcens, J5, color='red', label='A $>$ 38')

	pl.plot(lEcens, Ja, "ok")
	pl.errorbar(lEcens, Ja, yerr=Jerrors, fmt = "ok")

	pl.legend(fontsize=15, frameon=True)
	pl.semilogy()
	pl.ylim(1e-5)
	pl.title(plottitle)
	pl.grid()
	pl.ylabel('$J(E)$ [a.u.]')
	pl.xlabel('$\log_{10}$(E/eV)')
	pl.savefig(fileout)

def plot_errors_parts(title = "new_simulation",plotfile = "new_simulation",plottitle = "My simulation",parts = 1, logemin = 18, logemax = 20.4):

	output_dir()
	
	mask = (ecens >= logemin) & (ecens <= logemax)
	mask_bins = (ebins >= logemin) & (ebins <= logemax)
	ebins_ = ebins[mask_bins]
	ecens_ = ecens[mask]
	auger_ = auger[mask]
	sauger_ = sauger[mask]
	nauger_ = nauger[mask]
	sigma_auger_ = sigma_auger[mask]

	# load events
	filename = ['output/'+title+str(i)+'.dat' for i in range(parts)]
	fileout = 'output/'+plotfile+'.png'
	
	lEbins = ebins_  # logarithmic bins
	lEcens = (lEbins[1:] + lEbins[:-1]) / 2  # logarithmic bin centers
	dE = 10**lEbins[1:] - 10**lEbins[:-1]  # bin widths
	
	count = 0

	for i_parts in range(parts):
	
		if not_empty_file(filename[i_parts]):
			d = pl.genfromtxt(filename[i_parts], names=True)
			d2 = pl.genfromtxt(filename[i_parts])
			
			if (d2.ndim == 2):
				# observed quantities
				Z = pl.array([chargeNumber(id) for id in d['ID'].astype(int)])  # element
				A = pl.array([massNumber(id) for id in d['ID'].astype(int)])  # atomic mass number
				lE = pl.log10(d['E']) + 18  # energy in log10(E/eV))
				
				# identify mass groups
				idx1 = A == 1
				idx2 = (A > 1) * (A <= 4)
				idx3 = (A > 4) * (A <= 22)
				idx4 = (A > 22) * (A <= 38)
				idx5 = (A > 38)

				if (count == 0):
					# calculate spectrum: J(E) = dN/dE
					J  = pl.histogram(lE, bins=lEbins)[0] / dE
					J1 = pl.histogram(lE[idx1], bins=lEbins)[0] / dE
					J2 = pl.histogram(lE[idx2], bins=lEbins)[0] / dE
					J3 = pl.histogram(lE[idx3], bins=lEbins)[0] / dE
					J4 = pl.histogram(lE[idx4], bins=lEbins)[0] / dE
					J5 = pl.histogram(lE[idx5], bins=lEbins)[0] / dE
				else:
					J  = J + pl.histogram(lE, bins=lEbins)[0] / dE
					J1 = J1 + pl.histogram(lE[idx1], bins=lEbins)[0] / dE
					J2 = J2 + pl.histogram(lE[idx2], bins=lEbins)[0] / dE
					J3 = J3 + pl.histogram(lE[idx3], bins=lEbins)[0] / dE
					J4 = J4 + pl.histogram(lE[idx4], bins=lEbins)[0] / dE
					J5 = J5 + pl.histogram(lE[idx5], bins=lEbins)[0] / dE
					
				count += 1

	
	# calculate spectrum: J(E) = dN/dE for AUGER
	Ja = plt.hist(ecens_,bins=ebins_,weights=auger_)[0] / dE
	Ja = auger_ / dE
	Jerrors = sigma_auger_ / dE
	
	# normalize
	J1 /= J[0]
	J2 /= J[0]
	J3 /= J[0]
	J4 /= J[0]
	J5 /= J[0]
	J /= J[0]
	
	# calculate spectrum: J(E) = dN/dE for AUGER
	Ja = plt.hist(ecens_,bins=ebins_,weights=auger_)[0] / dE
	Ja = auger_ / dE
	Jerrors = sigma_auger_ / dE
	
	# normalize AUGER
	Jerrors /= Ja[0]
	Ja /= Ja[0]
	

	pl.figure(figsize=(10,7))
	pl.plot(lEcens, J,  color='SaddleBrown', label='Total')
	pl.plot(lEcens, J1, color='blue', label='A = 1')
	pl.plot(lEcens, J2, color='grey', label='A = 2-4')
	pl.plot(lEcens, J3, color='green', label='A = 5-22')
	pl.plot(lEcens, J4, color='purple', label='A = 23-38')
	pl.plot(lEcens, J5, color='red', label='A $>$ 38')

	pl.plot(lEcens, Ja, "ok")
	pl.errorbar(lEcens, Ja, yerr=Jerrors, fmt = "ok")

	pl.legend(fontsize=15, frameon=True)
	pl.semilogy()
	pl.ylim(1e-6)
	pl.title(plottitle)
	pl.grid()
	pl.ylabel('$J(E)$ [a.u.]')
	pl.xlabel('$\log_{10}$(E/eV)')
	pl.savefig(fileout)



def chi2_auger(title = "new_simulation0", logemin = 18, logemax = 20.4):

	output_dir()
	
	mask = (ecens >= logemin) & (ecens <= logemax)
	mask_bins = (ebins >= logemin) & (ebins <= logemax)
	ebins_ = ebins[mask_bins]
	ecens_ = ecens[mask]
	auger_ = auger[mask]
	sauger_ = sauger[mask]
	nauger_ = nauger[mask]
	
	sauger1_ = sauger1[mask]
	nauger1_ = nauger1[mask]
	sauger2_ = sauger2[mask]
	nauger2_ = nauger2[mask]
	sauger3_ = sauger3[mask]
	nauger3_ = nauger3[mask]
	sauger4_ = sauger4[mask]
	nauger4_ = nauger4[mask]
	sauger5_ = sauger5[mask]
	nauger5_ = nauger5[mask]


	# load events
	filename = 'output/'+title+'.dat'
	data = np.genfromtxt(filename, names=True)
	
	# observed quantities
	Z = pl.array([chargeNumber(id) for id in data['ID'].astype(int)])  # element
	A = pl.array([massNumber(id) for id in data['ID'].astype(int)])  # atomic mass number
	
	# identify mass groups
	idx1 = A == 1
	idx2 = (A > 1) * (A <= 4)
	idx3 = (A > 4) * (A <= 22)
	idx4 = (A > 22) * (A <= 38)
	idx5 = (A > 38)
	
	#Computing histograms for each kind of particle
	logE  = np.log10(data['E']) + 18
	hfull_test1 = plt.hist(logE[idx1],  bins=ebins_, histtype='stepfilled', alpha=0.5, label='Observed')
	test1 = hfull_test1[0]
	hfull_test2 = plt.hist(logE[idx2],  bins=ebins_, histtype='stepfilled', alpha=0.5, label='Observed')
	test2 = hfull_test2[0]
	hfull_test3 = plt.hist(logE[idx3],  bins=ebins_, histtype='stepfilled', alpha=0.5, label='Observed')
	test3 = hfull_test3[0]
	hfull_test4 = plt.hist(logE[idx4],  bins=ebins_, histtype='stepfilled', alpha=0.5, label='Observed')
	test4 = hfull_test4[0]
	hfull_test5 = plt.hist(logE[idx5],  bins=ebins_, histtype='stepfilled', alpha=0.5, label='Observed')
	test5 = hfull_test5[0]
	
	#Normalization of test
	n_test = test1.sum()+test2.sum()+test3.sum()+test4.sum()+test5.sum()
	n_test = float(n_test)
	test1 = test1/n_test
	test2 = test2/n_test
	test3 = test3/n_test
	test4 = test4/n_test
	test5 = test5/n_test
	stest1 = test1/n_test
	stest2 = test2/n_test
	stest3 = test3/n_test
	stest4 = test4/n_test
	stest5 = test5/n_test
	
	#Compute arrays with chi2 value for each bin and each kind of particle
	h1 = np.array(((nauger1_-test1)**2)/(sauger1_+stest1))
	h2 = np.array(((nauger2_-test2)**2)/(sauger2_+stest2))
	h3 = np.array(((nauger3_-test3)**2)/(sauger3_+stest3))
	h4 = np.array(((nauger4_-test4)**2)/(sauger4_+stest4))
	h5 = np.array(((nauger5_-test5)**2)/(sauger5_+stest5))
	
	h1 = h1[np.isfinite(h1)]
	h2 = h2[np.isfinite(h2)]
	h3 = h3[np.isfinite(h3)]
	h4 = h4[np.isfinite(h4)]
	h5 = h5[np.isfinite(h5)]
	
	chi2sum = h1.sum() + h2.sum() + h3.sum() + h4.sum() + h5.sum()
	
	return chi2sum


def chi2_particles_auger(title = "new_simulation0", logemin = 18, logemax = 20.4):

	output_dir()
	
	mask = (ecens >= logemin) & (ecens <= logemax)
	mask_bins = (ebins >= logemin) & (ebins <= logemax)
	ebins_ = ebins[mask_bins]
	ecens_ = ecens[mask]
	auger_ = auger[mask]
	sauger_ = sauger[mask]
	nauger_ = nauger[mask]
	
	sauger1_ = sauger1[mask]
	nauger1_ = nauger1[mask]
	sauger2_ = sauger2[mask]
	nauger2_ = nauger2[mask]
	sauger3_ = sauger3[mask]
	nauger3_ = nauger3[mask]
	sauger4_ = sauger4[mask]
	nauger4_ = nauger4[mask]
	sauger5_ = sauger5[mask]
	nauger5_ = nauger5[mask]

	# load events
	filename = 'output/'+title+'.dat'
	data = np.genfromtxt(filename, names=True)
	
	# observed quantities
	Z = pl.array([chargeNumber(id) for id in data['ID'].astype(int)])  # element
	A = pl.array([massNumber(id) for id in data['ID'].astype(int)])  # atomic mass number
	
	# identify mass groups
	idx1 = A == 1
	idx2 = (A > 1) * (A <= 4)
	idx3 = (A > 4) * (A <= 22)
	idx4 = (A > 22) * (A <= 38)
	idx5 = (A > 38)
	
	#Computing histograms for each kind of particle
	logE  = np.log10(data['E']) + 18
	hfull_test1 = plt.hist(logE[idx1],  bins=ebins_, histtype='stepfilled', alpha=0.5, label='Observed')
	test1 = hfull_test1[0]
	hfull_test2 = plt.hist(logE[idx2],  bins=ebins_, histtype='stepfilled', alpha=0.5, label='Observed')
	test2 = hfull_test2[0]
	hfull_test3 = plt.hist(logE[idx3],  bins=ebins_, histtype='stepfilled', alpha=0.5, label='Observed')
	test3 = hfull_test3[0]
	hfull_test4 = plt.hist(logE[idx4],  bins=ebins_, histtype='stepfilled', alpha=0.5, label='Observed')
	test4 = hfull_test4[0]
	hfull_test5 = plt.hist(logE[idx5],  bins=ebins_, histtype='stepfilled', alpha=0.5, label='Observed')
	test5 = hfull_test5[0]
	
	#Normalization of test
	n_test = test1.sum()+test2.sum()+test3.sum()+test4.sum()+test5.sum()
	n_test = float(n_test)
	test1 = test1/n_test
	test2 = test2/n_test
	test3 = test3/n_test
	test4 = test4/n_test
	test5 = test5/n_test
	stest1 = test1/n_test
	stest2 = test2/n_test
	stest3 = test3/n_test
	stest4 = test4/n_test
	stest5 = test5/n_test
	
	#Compute arrays with chi2 value for each bin and each kind of particle
	h1 = np.array(((nauger1_-test1)**2)/(sauger1_+stest1))
	h2 = np.array(((nauger2_-test2)**2)/(sauger2_+stest2))
	h3 = np.array(((nauger3_-test3)**2)/(sauger3_+stest3))
	h4 = np.array(((nauger4_-test4)**2)/(sauger4_+stest4))
	h5 = np.array(((nauger5_-test5)**2)/(sauger5_+stest5))
	
	h1 = h1[np.isfinite(h1)]
	h2 = h2[np.isfinite(h2)]
	h3 = h3[np.isfinite(h3)]
	h4 = h4[np.isfinite(h4)]
	h5 = h5[np.isfinite(h5)]
	
	chi2sum = h1.sum() + h2.sum() + h3.sum() + h4.sum() + h5.sum()
	
	return [h1.sum(),h2.sum(),h3.sum(),h4.sum(),h5.sum(),chi2sum]


def chi2_global_auger(title = "new_simulation", logemin = 18, logemax = 20.4):

	output_dir()
	
	mask = (ecens >= logemin) & (ecens <= logemax)
	mask_bins = (ebins >= logemin) & (ebins <= logemax)
	ebins_ = ebins[mask_bins]
	ecens_ = ecens[mask]
	auger_ = auger[mask]
	sauger_ = sauger[mask]
	nauger_ = nauger[mask]

	# load events
	filename = 'output/'+title+'.dat'
	data = np.genfromtxt(filename, names=True)
	
	# observed quantities
	Z = pl.array([chargeNumber(id) for id in data['ID'].astype(int)])  # element
	A = pl.array([massNumber(id) for id in data['ID'].astype(int)])  # atomic mass number
	
	# identify mass groups
	idx1 = A == 1
	idx2 = (A > 1) * (A <= 4)
	idx3 = (A > 4) * (A <= 22)
	idx4 = (A > 22) * (A <= 38)
	idx5 = (A > 38)
	
	#Computing histograms for each kind of particle
	logE  = np.log10(data['E']) + 18
	hfull_test = plt.hist(logE,  bins=ebins_, histtype='stepfilled', alpha=0.5, label='Observed')
	test = hfull_test[0]
	
	#Normalization of test
	n_test = test.sum()
	n_test = float(n_test)
	test = test/n_test
	stest = test/n_test
	
	#Compute arrays with chi2 value for each bin and each kind of particle
	h = np.array(((nauger_-test)**2)/(sauger_+stest))
	
	h = h[np.isfinite(h)]

	chi2sum = h.sum()
	
	return chi2sum

def chi2_global_auger_parts(title = "new_simulation", parts = 1, logemin = 18, logemax = 20.4):

	output_dir()
	
	mask = (ecens >= logemin) & (ecens <= logemax)
	mask_bins = (ebins >= logemin) & (ebins <= logemax)
	ebins_ = ebins[mask_bins]
	ecens_ = ecens[mask]
	auger_ = auger[mask]
	sauger_ = sauger[mask]
	nauger_ = nauger[mask]
	
	# load events
	filename = ['output/'+title+str(i)+'.dat' for i in range(parts)]
	count = 0
	test = np.array([0. for i in range(len(auger_))])
	
	for i_parts in range(parts):
	
		if not_empty_file(filename[i_parts]):
			data = np.genfromtxt(filename[i_parts], names=True)
			datab = np.genfromtxt(filename[i_parts])
			
			if (datab.ndim == 2):
				# observed quantities
				Z = pl.array([chargeNumber(id) for id in data['ID'].astype(int)])  # element
				A = pl.array([massNumber(id) for id in data['ID'].astype(int)])  # atomic mass number
				
				# identify mass groups
				idx1 = A == 1
				idx2 = (A > 1) * (A <= 4)
				idx3 = (A > 4) * (A <= 22)
				idx4 = (A > 22) * (A <= 38)
				idx5 = (A > 38)
						
				#Computing histograms for each kind of particle
				logE  = np.log10(data['E']) + 18
				hfull_test = plt.hist(logE,  bins=ebins_, histtype='stepfilled', alpha=0.5, label='Observed')
				
				if (count == 0):
					test = hfull_test[0]
				else:
					test = test + hfull_test[0]
					
				count += 1
	
	#Normalization of test
	n_test = test.sum()
	n_test = float(n_test)
	test = test/n_test
	stest = test/n_test
	
	sigmas = sauger_ + stest
	test = test[sigmas != 0]
	nauger_ = nauger_[sigmas != 0]
	sigmas = sigmas[sigmas != 0]
	
	#Compute arrays with chi2 value for each bin and each kind of particle
	h = np.array(((nauger_-test)**2)/(sigmas))
	
	#h = h[np.isfinite(h)]

	chi2sum = h.sum()
	
	return chi2sum


def chi2_global_auger_rcut(title = "new_simulation",rcut = 21., logemin = 18, logemax = 20.4):

	output_dir()
	
	mask = (ecens >= logemin) & (ecens <= logemax)
	mask_bins = (ebins >= logemin) & (ebins <= logemax)
	ebins_ = ebins[mask_bins]
	ecens_ = ecens[mask]
	auger_ = auger[mask]
	sauger_ = sauger[mask]
	nauger_ = nauger[mask]

	# load events
	filename = 'output/'+title+'.dat'
	data = np.genfromtxt(filename, names=True)
	
	# observed quantities
	Z = pl.array([chargeNumber(id) for id in data['ID'].astype(int)])  # element
	Z0 = pl.array([chargeNumber(id) for id in data['ID0'].astype(int)])  # element
	A = pl.array([massNumber(id) for id in data['ID'].astype(int)])  # atomic mass number
	num = len(Z0)
	energy0 = data['E0']
	
	# identify mass groups
	idx1 = A == 1
	idx2 = (A > 1) * (A <= 4)
	idx3 = (A > 4) * (A <= 22)
	idx4 = (A > 22) * (A <= 38)
	idx5 = (A > 38)
	
	# Modifying distribution according to rcut
	weight_rcut = np.array([f_cut(energy0[i]*(10**18.),Z0[i]*(10**rcut)) for i in range(num)])
	
	#Computing histograms for each kind of particle
	logE  = np.log10(data['E']) + 18
	hfull_test = plt.hist(logE,  bins=ebins_, weights=weight_rcut, histtype='stepfilled', alpha=0.5, label='Observed')
	test = hfull_test[0]
	
	#Normalization of test
	n_test = test.sum()
	n_test = float(n_test)
	test = test/n_test
	stest = test/n_test
	
	sigmas = sauger_ + stest
	test = test[sigmas != 0]
	nauger_ = nauger_[sigmas != 0]
	sigmas = sigmas[sigmas != 0]
	
	#Compute arrays with chi2 value for each bin and each kind of particle
	h = np.array(((nauger_-test)**2)/(sigmas))
	
	chi2sum = h.sum()
	
	return chi2sum

def chi2_global_auger_rcut_separate(fractions,title = "new_simulation",rcut = 21., logemin = 18, logemax = 20.4):

	mask = (ecens >= logemin) & (ecens <= logemax)
	mask_bins = (ebins >= logemin) & (ebins <= logemax)
	ebins_ = ebins[mask_bins]
	ecens_ = ecens[mask]
	auger_ = auger[mask]
	sauger_ = sauger[mask]
	nauger_ = nauger[mask]
	
	first = True
	
	elems = ['H','He','N','Si','Fe']
	
	for i in range(len(elems)):

		
		if round(fractions[i],7)==0:
			continue
				
		# load events
		filename = 'output/'+elems[i]+'_'+title+'.dat'
		data = np.genfromtxt(filename, names=True)
		
		# observed quantities
		Z = pl.array([chargeNumber(id) for id in data['ID'].astype(int)])  # element
		Z0 = pl.array([chargeNumber(id) for id in data['ID0'].astype(int)])  # element
		A = pl.array([massNumber(id) for id in data['ID'].astype(int)])  # atomic mass number
		num = len(Z0)
		energy0 = data['E0']
		
		# identify mass groups
		idx1 = A == 1
		idx2 = (A > 1) * (A <= 4)
		idx3 = (A > 4) * (A <= 22)
		idx4 = (A > 22) * (A <= 38)
		idx5 = (A > 38)
		
		# Modifying distribution according to rcut
		weight_rcut = np.array([f_cut(energy0[w]*(10**18.),Z0[w]*(10**rcut)) for w in range(num)])
		weight_rcut *= fractions[i]
		
		#Computing histograms for each kind of particle
		logE  = np.log10(data['E']) + 18
		hfull_test = plt.hist(logE,  bins=ebins_, weights=weight_rcut, histtype='stepfilled', alpha=0.5, label='Observed')
		
		if first:
			test = hfull_test[0]
			first = False
		else:
			test += hfull_test[0]
	
	#Normalization of test
	n_test = test.sum()
	n_test = float(n_test)
	test = test/n_test
	stest = test/n_test
	
	sigmas = sauger_ + stest
	test = test[sigmas != 0]
	nauger_ = nauger_[sigmas != 0]
	sigmas = sigmas[sigmas != 0]
	
	#Compute arrays with chi2 value for each bin and each kind of particle
	h = np.array(((nauger_-test)**2)/(sigmas))
	
	chi2sum = h.sum()
	
	return chi2sum
