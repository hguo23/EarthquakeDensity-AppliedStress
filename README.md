# EarthquakeDensity-AppliedStress
# Density-AppliedStress    
This project is for studying the relation between earthquake density and the applied stress perturbations    


##
WORKFLOW: ./Mainshock-Aftershock Clusters. 

FOLDERS  

	data:  
		qtm_b0.8_mc1_d1.6.txt (declustering result)
			NUMBER YEAR MONTH DAY HOUR MINUTE SECOND LATITUDE LONTITUDE DEPTH MAGNITUDE TYPE ClusterID EventID. 
	
		sc_FM.txt (focal mechanism information). 
			Time (in MATLAB datenum), Lat, Lon, Dep, Mag, Strike1, Strike2, Dip1, Dip2, Rake1, Rake2  
			
	

	code:  
		lldistkm.m (compute epicentral distance)
		SelectFM.m
		RadialThreeDNND.m (compute 3D radial density)
		calP_window.m (bin stress, fit density with Gaussian distritbuion)
		ConfidenceInter.m (determine confidence interval)
		FitSlope_error.m (fit slopes). 


STEPs  
1: FMD analysis for magnitude distribution (codes description in the following)

2: NND analysis for earthquake clustering (codes description in the following)  

	outputfile example: data/qtm_b0.8_mc1_d1.6.txt

3: Data arrangement:
	
	PreData.m  
	outputfile: output/qtm_2D.mat  

4: Calculate static stresses caused by mainshocks (codes description in the following and in the manuscript)
	
	qtm2D.py     
	outputfile example: output/qtm_2D_CFS0.mat  

5: Dynamic stress, density   
	
	ComputeDensity.m
	outputfile example: output/qtm_dens_0.mat

6: Bin the applied stresses, gaussian fit for data in each bin, display the result  
	
	run_fit_qtm.m



##
WORKFLOW: ./indSeism

FOLDERS  

	data:
		info_Paralana.mat
			Time Lon Lat Dep Mag	

	code:
		Pp_iso.m (pore pressure changes in the isotropic model)
		poro_strain_iso.m (poroelastic strain in the isotropic model)
		strain_to_stress_iso.m
		RadialThreeDNND_t.m (compute 3D radial density)
		calP_window.m (bin stress, fit density with Gaussian distritbuion)
		ConfidenceInter.m (determine confidence interval)
		FitSlope_error.m (fit slopes)

Steps:   
1: FMD analysis for magnitude distribution   
2: Calculate pore pressure and poroelastic stress based on injection history and locations
	see code/process.m

3: k neareat-neighbors distances for determining earthquake density
	see code/RadialThreeDNND_t.m

4: Bin the applied stresses, gaussian fit for data in each bin, display the result
	run_Paralana.m
   

## METHODS DESCRIPTION   
###########################################################################   
FMD analysis for magnitude distribution that follow Gutenberg_Richter relationship: log(N) = a - bM   

References:   
1. Aki, K., 1965, Maximum likelihood estimate of b in the formula log N = a - bM and its confidence limits: Bull. Earthquake Res. Inst., Tokyo Univ., v. 43, p. 237–239.   
2. Clauset, A., Shalizi, C.R., and Newmann, M.E.J., 2009, Power-law distributions in empirical data: SIAM review, v. 51, no. 4, p. 661–703.   
3. Goebel, T.H.W., Kwiatek, G., Becker, T.W., Brodsky, E.E., and Dresen, G., 2017, What allows seismic events to grow big?: Insights from b-value and fault roughness analysis in laboratory stick-slip experiments: Geology, v. 45, no. 9, p. 815–818, doi: 10.1130/G39147.1.   

Codes available at:    
https://github.com/tgoebel/magnitude-distribution    


###########################################################################    
Seismicity Clustering Analysis Based on nearest neighbor distances in a space-time-magnitude domain    

References:    
1. Zaliapin, I., and Ben-Zion, Y., 2013, Earthquake clusters in southern California I: Identification and stability: Journal of Geophysical Research: Solid Earth, v. 118, no. 6, p. 2847–2864, doi: 10.1002/jgrb.50179.    
2. Goebel, T.H.W., Rosson, Z., Brodsky, E.E., and Walter, J.I., 2019, Aftershock deficiency of induced earthquake sequences during rapid mitigation efforts in Oklahoma: Earth and Planetary Science Letters, v. 522, p. 135–143, doi: 10.1016/j.epsl.2019.06.036.    

Codes availabel at:      
https://github.com/tgoebel/clustering-analysis   


###########################################################################    
Static stress calculation    

Reference:     
Okada, Y., 1992, Internal deformation due to shear and tensile faults in a half-space, Bull. Seism. Soc. Am., 82, 1018-1040.     

Codes available at:     
1. disloc3d available on Paul Segall's website (https://pangea.stanford.edu/research/CDFM/software/index.html).
2. https://github.com/tbenthompson/okada_wrapper.git
