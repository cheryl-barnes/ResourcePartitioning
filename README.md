# Overview
These files detail analyses used to estimate the degree of resource partitioning (as a proxy for the potential for competition) between Pacific Halibut and Arrowtooth Flounder (30 to 69 cm fork length) in the Gulf of Alaska. 

<b>Code Citation</b>: Barnes, C.L., A.H. Beaudreau, M.E. Hunsicker, and L. Ciannelli. In Press. Assessing the potential for competition between Pacific Halibut (<i>Hippoglossus stenolepis</i>) and Arrowtooth Flounder (<i>Atheresthes stomias</i>) in the Gulf of Alaska. PLOS ONE. Forth-coming.

# Data Sources
We used standardized survey data procured from the Alaska Fisheries Science Center (National Marine Fisheries Service, National Oceanic and Atmospheric Association [NOAA]). Bottom trawl survey data (1990 to 2017) were collected by the Resource Assessment and Conservation Engineering (RACE) Division and are publicly accessible here: https://www.afsc.noaa.gov/RACE/groundfish/survey_data/data.htm.
Food habits data (1990 to 2013) were provided by the Resource Ecology and Ecosystem Modeling (REEM) Program and are publicly accessible here: https://access.afsc.noaa.gov/REEM/WebDietData/DietDataIntro.php. All the data necessary to complete the following analyses can be found in the 'Data' folder. See von Szalay and Raring (2016) and Livingston et al. (2017) for data collection methods.

# Spatial Overlap
The 'SpatialAnalyses' script file includes the code necessary to construct a delta model for estimating spatial overlap between Pacific Halibut and Arrowtooth Flounder in the Gulf of Alaska. Methods were modified from Hunsicker et al. (2013) and Shelton et al. (2017).

We used generalized additive models to quantify and predict the probability of occurrence and relative abundance of Pacific Halibut and Arrowtooth Flounder across a uniform grid spanning the study area. We multiplied probabilities of occurrence and relative abundances to estimate overall abundance in each survey year-grid cell combination. We then multiplied standardized abundances of each species to estimate spatial overlap in each combination of survey year and grid cell. 

# Dietary Overlap
The 'DietaryAnalyses' script file includes the code necessary to calculate species-specific proportions of prey by weight in each survey year and grid cell. We then used proportions of prey by weight to calculate Schoener's index of dietary overlap.

# Resource Partitioning
The 'ResourcePartitioningAnalyses' script file combines spatial overlap and dietary overlap to quantify the correlation between the two measures and thus the degree of resource partitioning between Pacific Halibut and Arrowtooth Flounder in the Gulf of Alaska. Descriptions of resource partitioning can be found in Schoener (1974) and Ross (1986).

# References
Chipps SR, Garvey JE. Assessment of diets and feeding patterns. In: Analysis and interpretation of freshwater fisheries data. Guy CS, Brown ML, editors. Bethesda: American Fisheries Society; 2007. pp. 473–514. <br>

Hunsicker ME, Ciannelli L, Bailey KM, Zador S, Stige L. Climate and demography dictate the strength of predator-prey overlap in a subarctic marine ecosystem. PLOS ONE. 2013;8(6):e66025. doi:10.1371/journal.pone.006602 <br>

Livingston PA, Aydin K, Buckley TW, Lang GM, Yang MS, Miller BS. Quantifying food web interactions in the North Pacific – a data-based approach. Environ Biol Fishes. 2017;100(4):443–470. <br>

Ross ST. Resource partitioning in fish assemblage: review of field studies. Copeia. 1986;1986(2):352–388.<br>

Schoener TW. Resource partitioning in ecological communities. Science. 1974;185(4145):27–39.<br>

Shelton AO, Hunsicker ME, Ward EJ, Feist BE, Blake R, Ward CL, et al. Spatio-temporal models reveal subtle changes to demersal communities following the Exxon Valdez oil spill. ICES J Mar Sci. 2017. doi: 10.1093/icesjms/fsx079 <br>

von Szalay PG, Raring NW. Data report: 2015 Gulf of Alaska bottom trawl survey. Seattle (WA): National Oceanic and Atmospheric Administration; 2016. Technical Memorandum: NMFS-AFSC-325. Sponsored by the US Department of Commerce.

# Acknowledgments
Wayne Palsson and Kirstin Holsman queried and provided guidance regarding the use of AFSC databases. Analytical recommendations from Franz Mueter and Jordan Watson were incorporated throughout.
