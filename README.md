# PNR_Beetles
Data and analysis code associated with ground beetle projects at Powdermill Nature Reserve. 
Part of this folder relates to Aaron Tayal's M.S. thesis, titled "The effects of a tornado 
and salvage-logging on ground beetles".

We sampled ground beetles in a forest that had been affected by a natural disturbance (a 
tornado) in 2012, and an anthropogenic disturbance (salvage-logging) in 2013. We sampled
ground beetles in both 2015 and 2022, using barrier pitfall traps. Ground beetles were
identified to species, and functional traits of these species were measured.

Here are some folders in this repository:

PNR_raw_data: folder contains raw data as it was entered initially
-->"PNR2022_Aaron_pinned_InvertebrateCommunity.xslx" contains raw counts of carabids and
also silphids in 2022. The scarab and other beetle counts are not complete - see the Perry
Lab Microsoft Teams for the complete data.
-->"PNR_SpeciesTraits_2022.xlsx" contains the data used for the functional trait analysis
in Aaron's MS thesis. The first page is metadata, and the second page has each specimen
listed as a row, and each column as a trait that was measured.

Aaron_PNR_formatted_data: folder contains spreadsheets that are cleaned or summarized 
versions of the raw data.
-->"PNR2015_2022_carabid_counts.csv" has the ground beetle counts from 2015 and 2022 
organized into the same spreadsheet. This file is the input data for the R script I used
for my MS thesis.
-->"PNR2015_2022_carabid_counts_by_plot_standardized.csv" has ground beetle counts that
have been standardized for sampling effort. See Aaron's MS thesis for the method
used to standardize counts for sampling effort. Importantly, this table also has the 
alpha-diversity indices and community-weighted-mean (CWM) trait values listed as columns.
It also has values for the NMDS, which was used to investigate beta-diversity. The
"PC1" column is the CWM for principal component axis 1. Note that in this spreadsheet,
counts have been summed across all sampling intervals in a season.
-->"Aaron_formatted_PNR_PitfallTrapLocations_2015.csv" has information about the locations
of each plot within the study site. This includes info on transects, treatment groups,
GPS coordinates, and elevations. Note that in 2022, Plot 63 was unnaccessible, so the
pitfall trap was placed at Plot 65, which is very close to Plot 63 from 2015.

R_scripts:
-->"2015_2022_carabids_analysis.R" is the main R script. It standardizes the counts and sums
across sampling intervals. It runs species accumulation curves. It calculates taxonomic and 
functional alpha- and beta-diversity. It calculates the activity-abundance of open-habitat
and forest specialists. It calculates community-weighted-mean trait values. It exports 
spreadsheets.
-->"UPDATED_carabid_traits.R" is the script which cleans and summarizes the functional trait 
data. Traits are standardized to body length, correlation between traits is calculated,
principal component analyses are run, and a Gower distance matrix between species is
calculated. The Gower distance matrix is written into a spreadsheet, 
"carabid_dist_in_trait_space.csv" which is used by the main R script ("2015_2022_carabids_
analysis.R") to calculate metrics of functional diversity.
-->"carabids_stats.R" is the script which runs stats. Specifically, it runs stats to test
for differences based on forest management treatment and/or year of sampling. The response
variables include activity-abundance, species richness, Shannon diversity, functional
alpha-diversity, community-weighted mean traits, and environmental variables such as 
canopy openness and vegetation percentage cover.
-->"which_carabids_were_caught_in_2015_41_65.R" In 2015, more plots were sampled than just 
the 24 included in Aaron's MS thesis. This script simply subsets the 2015 data so that it
only includes data for the 24 plots that were also sampled in 2022.
-->"PNR_environmental_vars_exploratory.R" In 2015 and 2022, environmental variables were
recorded in multiple datasheets. This script combines and summarizes the environmental data,
typically by averaging the observations across the season. It also makes some graphs
to look at changes to env. variables across the season, such as soil moisture. It also
examines the correlations between different variables.
-->"weather_data.R" This script uses downloaded weather data from the NOAA 
<https://www.ncei.noaa.gov/access/past-weather/40.31510524031136,-79.72542724948556,40.10676305409062,-79.22526387171922>
in order to look at temperature and precipitation trends at Powdermill Nature Reserve over
the study period.









