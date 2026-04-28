# R CODE and DATABASE for the paper: "Metabolic responses of mammals to heat and cold waves vary across climates"

Abstract

Aim: Extreme temperature fluctuations are becoming more frequent under climate change, which may have an impact on the energy and water balance of organisms. Yet, the magnitude of these impacts and whether species from different climatic regions differ in their metabolic responses to cold and heat waves remain poorly understood. Here, we examine how short-term temperature fluctuations influence metabolic requirements of mammals across broad climatic gradients.

Location: Global.

Time period: 1966 - 2007 

Major taxa studied: Mammals

Methods: We compiled field metabolic rate (FMR) data from 390 individuals representing 41 mammal species worldwide. For each measurement, we extracted ambient temperature at the time and location of FMR collection using historical weather data and quantified temperature anomalies relative to long-term mean conditions. We modelled FMR as a function of temperature anomaly, long-term mean temperature, and their interaction.

Results: Metabolic responses to temperature anomalies varied systematically across climatic regimes. Mammals from colder environments increased FMR during both negative (cold) and positive (heat) anomalies. In contrast, species from warmer regions showed stronger increases in FMR during negative anomalies but exhibited little change or even decreases in FMR during positive anomalies.

Main conclusions: Mammalian metabolic responses to short-term temperature fluctuations are strongly climate dependent, likely reflecting divergent adaptations to local thermal regimes. These findings highlight the importance of considering climatic context when predicting the energetic consequences of increasing thermal extremes under climate change.

 # Data description (columns)
"Class": Taxonomic classification (string)       
"Order": Taxonomic classification (string)        
"Family": Taxonomic classification (string)  
"Genus": Taxonomic classification (string)  
"Species": Taxonomic classification (string)  
"Mass": Body mass (g)
"Lat_deg": Latitude (deg) 
"Lon_deg": Longitude (deg)      
"Elevation": Elevation (m)    
"Date_original": Date of measurement
"Date_start": Starting date of DLW measurement
"Date_end": End date of DLW measurement    
"Jan": Month included in DWL average (0 - 1)         
"Feb"          
"Mar"     
"Apr"      
"May"         
"Jun"      
"Jul"     
"Aug"          
"Sep"   
"Oct"    
"Nov"      
"Dec"      
"Year1": Year included in DWL average
"Year2"        
"Year3" 
"Year4"    
"FMR_Watt": FMR value (W)
"FMR_kJday": FMR value (kJ/day)
"TAREFmean": Air temperature at the time and location of FMR measurement (ºC) - reference height (1.2m)
"TALOCmean": Air temperature at the time and location of FMR measurement (ºC) - individual height 
"SOLRAD": Solar radiation at the time and location of FMR measurement (W/m2) 
"RHREFmean": Relative humidity at the time and location of FMR measurement (%) - reference height (1.2m) 
"RHLOCmean": Relative humidity at the time and location of FMR measurement (%) - individual height  
"FMR_M": Mass specific FMR (residuals) 
"Tanomalies": Temperature anomaly (ºC) 
"TMEAN": Mean temperature of the month of FMR measurement (ºC)        
"Reference": Complete reference for FMR data



 

   


 


