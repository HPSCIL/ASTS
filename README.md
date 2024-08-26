This code aims to provide the source code for a spatiotemporal adaptive sampling algorithm model.  
The model is applicable to research aimed at minimizing long-term cumulative prediction errors.  
The minimum requirements for model construction are: projection coordinates utmX, utmY (used for building the inla-spde mesh) and the target variable.If there are covariates involved in the model building, they can also be applied (just modify the expression as needed).  
This software package is compiled and executed in a Windows environment. 
The Mac environment requires independent configuration of packages.
>[!NOTE]
>The INLA package and related papers referenced for use are as follows.
>_https://www.r-inla.org_

The dataset is a processed dataframe type (UTM_X, UTM_Y, PM2017...2021).  
>[!NOTE]
>The raw data collection can be found at _[doi.org/10.5281/zenodo.3539349](url)_

The detailed usage steps for the code are written in the comments, which can be found in cube1 to cube5. 
The entry point for the program is at
>[!Important]
>**line 156 :asts_start(para.combinations, iterations = 5)**.
    
>[!Caution]
>It is highly recommended to read the content in the comments before using the code.














 
