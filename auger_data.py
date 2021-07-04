import numpy as np

#Energy bins for the histograms Auger
ebins = np.array([18.4019, 18.5004, 18.601, 18.702, 18.801, 18.8996, 18.9986, 19.0996, 19.2006, 19.3002, 19.3997, 19.5002, 19.6011, 19.7001, 19.7983, 19.8976, 19.9979, 20.0993, 20.1998, 20.2986, 20.3973])

ecens = np.array([18.45115, 18.5507 , 18.6515 , 18.7515 , 18.8503 , 18.9491 , 19.0491 , 19.1501 , 19.2504 , 19.34995, 19.44995, 19.55065, 19.6506 , 19.7492 , 19.84795, 19.94775, 20.0486 , 20.14955, 20.2492 , 20.34795])

#Data collected from Auger
auger1 = np.array([61975, 29143, 12885, 5338, 2039, 666, 181, 38, 6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])			#A=1
auger2 = np.array([20332, 17495, 13859, 10550, 8036, 5566, 3442, 1759, 630, 135, 14, 0, 0, 0, 0, 0, 0, 0, 0, 0])		#2<A<4
auger3 = np.array([1351, 1454, 1655, 1919, 2307, 2463, 2373, 2132, 1773, 1337, 796, 455, 232, 70, 21, 3, 0, 0, 0, 0])	#5<A<22
auger4 = np.array([0, 40, 45, 53, 69, 89, 114, 146, 174, 188, 157, 134, 114, 67, 44, 12, 3, 1, 0, 0])				#23<A<38
auger5 = np.array([0, 0, 0, 0, 0, 0, 0, 8, 11, 16, 19, 22, 24, 18, 16, 7, 4, 4, 4, 2])						#39<A<56
#auger = auger1+auger2+auger3+auger4+auger5
auger = np.array([76176,44904,26843,16970,12109,8515,5939,4048,2567,1664,979,619,373,152,80,23,9,6,0,0])

#Normalization of auger data
n_auger = auger.sum()
n_auger = float(n_auger)
nauger = auger/n_auger
nauger1 = auger1/n_auger
nauger2 = auger2/n_auger
nauger3 = auger3/n_auger
nauger4 = auger4/n_auger
nauger5 = auger5/n_auger

sauger = nauger/n_auger
sigma_auger = np.sqrt(auger)
sauger1 = nauger1/n_auger
sauger2 = nauger2/n_auger
sauger3 = nauger3/n_auger
sauger4 = nauger4/n_auger
sauger5 = nauger5/n_auger
