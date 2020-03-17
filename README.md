# mesocosm
mesocosm experiment and associated individual based disease model

Starting in March 2020, I will simulate community disassembly by creating a series of nested artificial plant communities to empirically test how changes in species composition impact the spread of a fungal plant pathogen. I will be using greenhouse mesocosms consisting of crop seedlings contained in mesocosms inoculated with the soilborne pathogen, Rhizoctonia solani. My primary goal is to examine how changes in host competency and density across a richness gradient affects disease. I also will predict infection probability as a function of a community’s average host competency, density, and richness to estimate the relative effects of species composition and richness per se. Scripts and outputs related to the empirical data generation can be found in folders starting with 'gh_'  

In efforts to futher explore how community composition affects disease risk, I wrote an individual-based disease model to accompany my mesocosm experiment. My data from the mesocosm experiment will allow me to validate my model so I can explore additional starting parameters. The model is a spatially explicit discrete time model that tracks state changes for all individuals over the duration of the epidemic. For details on how the model works, see https://github.com/lisamr/mesocosm/blob/master/IBM/scripts/math_for_IBM.Rmd
