# Pan-cancer_Lodietal
This repository contains code related to data processing and downstream analysis associated with the study "A pan-cancer dissection of the tumor immune microenvironment transcriptomic landscape", under submission at Cell Reports Medicine.

# Data Availability 
## Processed & filtered scRNAseq data
All processed scRNAseq data are available for in-browser exploration and download through the Data Access official portal https://lambrechtslab.sites.vib.be/en/dataaccess 

## Shiny Apps
We created interactive and publicly available Shiny Apps designed to provide the research community with an intuitive way of visualizing, exploring, and downloading our data, while using the annotation from our pan-cancer dataset (temporarely at https://marine-lab.shinyapps.io/PanCancer-Lodi/). 
Besides being able to examine expression of individual genes or customized gene signatures in each cell type across cancer types, these apps allow us to rank individual tumors based on a gene signature in any given cell type and correlate subcluster abundance with this rank.

# Contacts
All other relevant data and analysis are available from the authors upon request. For further enquires, please either raise an issue via GitHub or email Bram Boeckx (bram.boeckx@kuleuven.be), Diether Lambrechts (diether.lambrechts@kuleuven.be) or Francesca Lodi (francesca.lodi@kuleuven.be). 

# List of scripts 
1) Heatmap showing the (scaled) expression of curated gene signatures (columns) across subclusters (rows)
2) Box plots displaying the fractions of major cell types detected across cancer types
