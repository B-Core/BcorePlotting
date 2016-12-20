library(testthat)
require(limma)
require(data.table)
baseDir_v <- "~/stable_repos_11_17/"
source(paste0(baseDir_v, "BcorePlotting/tests/test_data.R"))
source(paste0(baseDir_v, "BcorePlotting/ClusteringPlots.R"))


# Get Attributes
myattribs <- list("Treatment" = metadata_dt$Treatment)
print("Treatments used for plotting: "); myattribs

makeHeatmap(sampleData_mat,
            list("Treatment" = sampleAttribs_ls$Treatment),
            plottitle = "This is a Heatmap",
            subtitle = "Bottom of the Heatmap",
            setColv = "Treatment",
            colOrder_v = unique(sampleAttribs_ls$Treatment))

makeHeatmap(sampleData_mat,
            list("Treatment" = sampleAttribs_ls$Treatment),
            plottitle = "This is a Heatmap",
            subtitle = "Bottom of the Heatmap")