### Summary

Repository of R and Stan code to simulate and analyse foraging trip data (patch choice and harvest success). 


### Accompanying manuscript

Socio-economic predictors of Inuit hunting strategies and their implications for climate change adaptation

Phil. Trans. R. Soc. B 378: 20220395. https://doi.org/10.1098/rstb.2022.0395

F. Hillemann, B. A. Beheim, E. Ready

This manuscript is part of the special issue ‘Climate change adaptation needs a science of culture,’ published in Philosophical Transactions of the Royal Society B in 2023.

**Manuscript Abstract:**  
In the Arctic, seasonal variation in the accessibility of the land, sea ice, and open waters influences which resources can be harvested safely and efficiently. Climate stressors are also increasingly affecting access to subsistence resources. Within Inuit communities, people differ in their involve- ment with subsistence activities, but little is know about how engagement in the cash economy (time and money available) and other socio-economic factors shape the food production choices of Inuit harvesters, and their ability to adapt to rapid ecological change. We analyse 281 forag- ing trips involving 23 Inuit harvesters from Kangiqsujuaq, Nunavik, using a Bayesian approach modelling both patch choice and within-patch success. Gender and income predict Inuit har- vest strategies: while men, especially men from low-income households, often visit patches with a relatively low success probability, women and high-income hunters generally have a higher propensity to choose low-risk patches. Inland hunting, marine hunting, and fishing differ in the required equipment and effort, and hunters may have to shift their subsistence activities if certain patches become less profitable or less safe due to high costs of transportation or climate change (e.g., navigate larger areas inland instead of targeting seals on the sea ice). Our finding that household income predicts patch choice suggests that the capacity to maintain access to country foods depends on engagement with the cash economy.

**Manuscript Keywords:**  
Arctic Canada, food security, hunting, Inuit, risk-sensitive foraging, socio-ecological systems


### Description of the file structure

The repository includes R code files to simulate and analyse harvest trip data, and Stan files to reproduce the patch choice model and the harvest success model presented in the accompanying manuscript.

The .R file also contains information on code prerequisites, and code sections to generate model summaries (table) and visualisations used in the manuscript.

**Content of the simulated data**  
The simulated dataset contains the following variable:

| Variable Name    | Meaning                                                |
|:-----------------|:-------------------------------------------------------|
| trip_id          | unique harvest episode                                 |
| season/season_id | season (1: snow/ice season, 2: ice-free season)        |
| Nhunt/Nhunt_s    | number of hunters during episode (standardised)        |
| j_id             | focal harvester's ID                                   |
| age_cat          | j_id's age category                                    |
| gender/gender_id | j_id's gender (1: male, 2: female)                     |
| income_s         | j_id's income, standardised                            |
| indegree_s       | j_id's indegree in food-sharing network, standardised  |
| outdegree_s      | j_id's outdegree in food-sharing network, standardised |
| patch_id         | harvest episode patch choice ID                        |
| harvest          | harvest episode success (1: yes, 0: no)                |
 


### Access information

Published paper: https://doi.org/10.1098/rstb.2022.0395

Preprint of the manuscript: https://doi.org/10.31235/osf.io/sm683

Code repository on Data Dryad: https://doi.org/10.5061/dryad.k3j9kd5dv 

Code repository on Github: https://github.com/fhillemann/MSrepo_harvest_patch_choice.git


### Citation

We value transparency and reproducibility in our research! Providing our methods and code allows others to validate and build upon our findings.

If you use our code or findings, please consider citing our manuscript and this repository (see CITATION.cff file).

Thank you for your interest in our work!
