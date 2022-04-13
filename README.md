# BiodiversityChangeInQuebec
Code and data for the paper:

Crockett ETH, M Vellend, and EM Bennett. Tree biodiversity in northern forests shows temporal stability over 35 years at different scales, levels, and dimensions. 

## Code

The R code in the file "BiodiversityChange_MainAnalyses.R" runs the analyses described in the above mentioned paper. This R script calls functions created in the file "BiodiversityChange_Functions.R" and reads in data from the DataFiles folder.

## Data

Data files called into the code are based in the "Data folder". These files include cleaned versions of the data sources referenced below, and list values of the variables for the forest plots used in this analysis.


'CommMatrixFiles' folders provide the community data matrices for the Plot, 50km, 100km, 200km spatial scales, with 100 replicates generated through the resampling procedure.

'Explanatory_variables' provides data values for climate variables, disturbance, and land designations for each of the spatial scales.

'ForestChange_byBioDomain' indicates the percent of harvesting and fire disturbances across the different biodomains.

'Placette_UTM19_Included' provides coordinates for the forest plots.

'QC_50km', 'QC_100km', and 'QC_200km' files provide the grid cells generated.

'QCphyloCodes' provides phylogenetic information for the tree species in this study.

'infoFirst_tiges' details information about each forest plot during the first sampling period (1970s).

'traits17' gives information on the 17 functional traits for the tree species found in this study.


## References

Données Quebec. (2017). Placettes-échantillons temporaires du quatrième inventaire. Retrieved from: https://www.donneesquebec.ca/recherche/fr/dataset/placettes-echantillons-temporaires-4e-inventaire-2004-a-aujourd-hui

Données Quebec. (2019). Couche des territoires récréatifs du Québec à l'échelle de 1/100 000. Retrieved from: https://www.donneesquebec.ca/recherche/dataset/couche-des-territoires-recreatifs-du-quebec-a-l-echelle-de-1-100-000

McKenney, D. W., Hutchinson, M. F., Papadopol, P., Lawrence, K., Pedlar, J., Campbell, K., Milewska, E., Hopkinson, R. F., Price, D., & Owen, T. (2011). Customized spatial climate models for North America. Bulletin of the American Meteorological Society, 92(12), 1611–1622.

McKenney, D. W., Pedlar, J. H., Papadopol, P., & Hutchinson, M. F. (2006). The development of 1901–2000 historical monthly climate models for Canada and the United States. Agricultural and Forest Meteorology, 138(1–4), 69–81. https://doi.org/10.1016/j.agrformet.2006.03.012

Paquette, A., & Messier, C. (2011). The effect of biodiversity on tree productivity: From temperate to boreal forests: The effect of biodiversity on the productivity. Global Ecology and Biogeography, 20(1), 170–180. https://doi.org/10.1111/j.1466-8238.2010.00592.x

Pedlar, J. H., McKenney, D. W., Lawrence, K., Papadopol, P., Hutchinson, M. F., & Price, D. (2015). A Comparison of Two Approaches for Generating Spatial Models of Growing-Season Variables for Canada. Journal of Applied Meteorology and Climatology, 54(2), 506–518. https://doi.org/10.1175/JAMC-D-14-0045.1

Qian, H., & Jin, Y. (2016). An updated megaphylogeny of plants, a tool for generating plant phylogenies and an analysis of phylogenetic community structure. Journal of Plant Ecology, 9(2), 233–239. https://doi.org/10.1093/jpe/rtv047

White, J. C., Wulder, M. A., Hermosilla, T., Coops, N. C., & Hobart, G. W. (2017). A nationwide annual characterization of 25 years of forest disturbance and recovery for Canada using Landsat time series. Remote Sensing of Environment, 194, 303–321. https://doi.org/10.1016/j.rse.2017.03.035
