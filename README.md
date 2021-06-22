# Welcome to sim_functions, a software to run crpropa simulations with particular energy and distance distributions

Here we describe the available functions via the header "sim_functions"

Requirements:
- python 2.7
- numpy, pandas, matplotlib
- CRPropa3

# DISTANCES

There are 2 lists of sources:

1. distances_v1.dat
	* Taken from: Véron-Cetty, M. -P. ; Véron, P. A catalogue of quasars and active nuclei: 13th edition (https://ui.adsabs.harvard.edu/abs/2010A&A...518A..10V)
		
2. distances_v2.dat
	* Provided by the Radio Fundamental Catalog (RFC). (http://astrogeo.org/sol/rfc/rfc_2021a/) Just sources that coincide with distances_v1.dat were taken.
		

# INSTALL CRPROPA3
Run:
````bash
	sudo bash install_crpropa.sh
```

Make sure that the make test works perfectly.

# FUNCTIONS

````python
	sim_functions.simulate_dd_1D_pl(params,distances,num = 1000000,title = "new_simulation",model = "Gilmore12",sim_seed = 0)
```
	
	It simulates the propagation of Cosmic Rays using a powerlaw distribution of energy from 1 to 1000 EeV without an R_cut.

	Parameters:	params: 1-dimensional ndarray or list
				The list must contain in order: the negative of the exponent gamma, the fractions of
				hydrogen, helium, nitrogen, silicum and iron, respectively (maximum value: 1.).
				
			distances: 2-dimensional array or list
				It must contain a list with 2 rows: the distance in Mpc and the weight of that source.
			
			num: int
				Number of particles at sources.
			
			title: string
				Name of the document in which the simulated data will be printed. No extension needed.
			
			model: string
				Model to do the simulation. Options: "Gilmore12" or "Dominguez11".
			
			sim_seed: int
				An integer to fix the seed for the simulation. If 0, no seed is fixed.


````python 
	sim_functions.simulate_dd_1D_parts(params,distances,num = 1000000,title = "new_simulation",model = "Gilmore12",parts = 1,n_e = 10000,seed = 0,sim_seed = 0)
```
	
	It simulates the propagation of Cosmic Rays using a powerlaw distribution of energy from 1 to 1000 EeV with an R_cut.

	Parameters:	params: 1-dimensional ndarray or list
				The list must contain in order: the negative of the exponent gamma, the log10(R_cut/1 eV), the fractions of
				hydrogen, helium, nitrogen, silicum and iron, respectively (maximum value: 1.).
				
			distances: 2-dimensional array or list
				It must contain a list with 2 rows: the distance in Mpc and the weight of that source.
			
			num: int
				Number of particles at sources.
			
			title: string
				Name of the document in which the simulated data will be printed. No extension needed.
			
			model: string
				Model to do the simulation. Options: "Gilmore12" or "Dominguez11".
			
			parts: int
				The number of parts in which the simulation will be done in order to the process not to be killed.
			
			n_e: int
				Number of energy bins considered for the simulation.
			
			seed: int
				An integer to fix the seed for the energies to be used for the distribution. If 0, no seed is fixed.
			
			sim_seed: int
				An integer to fix the seed for the simulation. If 0, no seed is fixed.


````python
	sim_functions.plot_power(title = "new_simulation",plotfile = "new_simulation_plot",plottitle = "My simulation",power = 2.)
```
	
	It plots the events of the simulated propagation of Cosmic Rays using CRPropa. X_axis: energy in EeV. Y_axis: dN/dE * E**power.
	The events are normalized according to the first bin of the total events.

	Parameters:	title: string
				Name of the document in which the simulated data is. No extension needed.
			
			plotfile: string
				Name of the .png where the plot will be saved. No extension needed.
			
			plottile: string
				Title of the plot. It can use LaTeX.
			
			power: float
				Power considered for the Y_axis of the plot. See main description of the function.


````python
	sim_functions.plot_parts(title = "new_simulation",plotfile = "new_simulation",plottitle = "My simulation",parts = 1)
```
	
	It plots the events of the simulated propagation of Cosmic Rays produced by the function sim_functions.simulate_dd_1D_parts. X_axis: energy in EeV. Y_axis: dN/dE.
	The events are normalized according to the first bin of the total events.

	Parameters:	title: string
				Name of the document in which the simulated data is. No extension needed.
			
			plotfile: string
				Name of the .png where the plot will be saved. No extension needed.
			
			plottile: string
				Title of the plot. It can use LaTeX.
			
			parts: int
				Number of part in which the simulation was divided.
				It must coincide with the "parts" parameter of the original function that produced the data.



````python
	sim_functions.plot_power_parts(title = "new_simulation",plotfile = "new_simulation",plottitle = "My simulation",parts = 1, power = 2.)
```
	
	It plots the events of the simulated propagation of Cosmic Rays produced by the function sim_functions.simulate_dd_1D_parts. X_axis: energy in EeV. Y_axis: dN/dE * E**power.
	The events are normalized according to the first bin of the total events.

	Parameters:	title: string
				Name of the document in which the simulated data is. No extension needed.
			
			plotfile: string
				Name of the .png where the plot will be saved. No extension needed.
			
			plottile: string
				Title of the plot. It can use LaTeX.
			
			parts: int
				Number of part in which the simulation was divided.
				It must coincide with the "parts" parameter of the original function that produced the data.

			power: float
				Power considered for the Y_axis of the plot. See main description of the function.



````python
	sim_functions.plot_errors_rcut(title = "new_simulation",plotfile = "new_simulation",plottitle = "My simulation",rcut = 21.)
```
	
	It plots the events of the simulated propagation of Cosmic Rays produced by the function sim_functions.simulate_dd_1D_pl. It automatically filters the output in order to add the R_cut effect. X_axis: energy in EeV. Y_axis: dN/dE.
	The events are normalized according to the first bin of the total events.

	Parameters:	title: string
				Name of the document in which the simulated data is. No extension needed.
			
			plotfile: string
				Name of the .png where the plot will be saved. No extension needed.
			
			plottile: string
				Title of the plot. It can use LaTeX.
			
			rcut: float
				log10(R_cut/1 eV). R_cut of the source energy distribution.



````python
	sim_functions.plot_errors_parts(title = "new_simulation",plotfile = "new_simulation",plottitle = "My simulation",parts = 1)
```
	
	It plots the events of the simulated propagation of Cosmic Rays produced by the function sim_functions.simulate_dd_1D_parts, including the errors with respect to the Auger data. X_axis: energy in EeV. Y_axis: dN/dE.
	The events are normalized according to the first bin of the total events.

	Parameters:	title: string
				Name of the document in which the simulated data is. No extension needed.
			
			plotfile: string
				Name of the .png where the plot will be saved. No extension needed.
			
			plottile: string
				Title of the plot. It can use LaTeX.
			
			parts: int
				Number of part in which the simulation was divided.
				It must coincide with the "parts" parameter of the original function that produced the data.


````python
	sim_functions.chi2_auger(title = "new_simulation0")
```

	It computes the chi2 of a simulation done in 1 part and distinguishes between different types of particles.
	Output: float, chi2.
	
	Parameters:	title: string
				Name of the document in which the simulated data is. No extension needed.

````python
	sim_functions.chi2_particles_auger(title = "new_simulation0")
```

	It computes the chi2 of a simulation done in 1 part and distinguishes between different types of particles.
	Output: list of floats, chi2 of each element and total chi2.
	
	Parameters:	title: string
				Name of the document in which the simulated data is. No extension needed.

````python
	sim_functions.chi2_global_auger(title = "new_simulation")
```

	It computes the chi2 of a simulation done in 1 part and doesn't distinguish between different types of particles.
	Output: float, chi2.
	
	Parameters:	title: string
				Name of the document in which the simulated data is. No extension needed.

````python
	sim_functions.chi2_global_auger_parts(title = "new_simulation", parts = 1)
```

	It computes the chi2 of a simulation done in several parts and doesn't distinguish between different types of particles.
	Output: float, chi2.
	
	Parameters:	title: string
				Name of the document in which the simulated data is. No extension needed.
			
			parts: int
				Number of part in which the simulation was divided.
				It must coincide with the "parts" parameter of the original function that produced the data.

````python
	sim_functions.chi2_global_auger_rcut(title = "new_simulation",rcut = 21.)
```

	It computes the chi2 considering the events of the simulated propagation of Cosmic Rays produced by the function sim_functions.simulate_dd_1D_pl. It automatically filters the output in order to add the R_cut effect.
	Output: float, chi2.
	
	Parameters:	title: string
				Name of the document in which the simulated data is. No extension needed.
			
			rcut: float
				log10(R_cut/1 eV). R_cut of the source energy distribution.

