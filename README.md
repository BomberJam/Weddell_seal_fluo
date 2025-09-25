Weddell seal fluo project

Mandatory data :

bathy folder

amsr2/ 
        2019 
        2020
        2023
        
fas-ice/ 
        mertz_sara_akiko_19.nc
        mertz_sara_akiko_20.nc
        mertz_sara_akiko_23.nc
        
oceanographic_data/
        wd11
        wd12
        wd20
        
Scripts :
1_data_correction_select.Rmd : Correcting lat/lon coordinates, light and Fluo chl-a. Create 1 rds for each individual.
To select an individual, modify the lines where "SealOfApproval" is written (Ctrl+F : SealOfApproval).
2_deploy_dataset.Rmd : Gather all the previously created Rds to generate Rds by deployment or a global Rds, with all the data.
3_general_analysis.Rmd : Analysis on the number of individuals, number of profiles, etc...
4_2019_Analysis.Rmd : 2019 data analysis (Time series, MLD, Maps, T/S diagram).
5_2020_Analysis.Rmd : 2020 data analysis (Time series, MLD, Maps, T/S diagram).
6_2023_Analysis.Rmd : 2023 data analysis (Time series, MLD, Maps, T/S diagram).
7_ACP.R : Statistical analysis of the data.
