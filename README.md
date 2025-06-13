# Introduction

This is the figure code repository for the MFG project.


## Open Commend

The merge request is open for code review and maintenance.


## Repository Structure


The code files are organized as follows in this repository for convenient access:


```
.
├── fig1
│   ├── fig1c_250120.R
│   ├── fig1d_250120.R
│   ├── fig1e_250120.R
│   ├── fig1f_250120.R
│   ├── fig1g_250120.R
│   ├── fig1h_250120.R
│   ├── fig1i_250120.R
│   └── plotfgsea.R
├── fig2
│   ├── fig2b2c_250120.R
│   ├── fig2d_250121.R
│   ├── fig2e_250121.R
│   └── fig2f_250121.R
├── fig3
│   ├── fig3D_241120.R
│   ├── fig3e_250121.R
│   ├── fig3f_3g_250121.R
│   ├── fig3h_250121.R
│   ├── fig3i_250121.R
│   └── plotfgsea.R
├── fig4
│   ├── fig4a_250121.R
│   ├── fig4c_250121.R
│   ├── fig4d_250121.R
│   ├── fig4e_4f_250121.R
│   ├── fig4e_4f_250312.R
│   ├── fig4g_250121.R
│   ├── fig4h_250121.R
│   ├── fig4i_250121.R
│   └── plotfgsea.R
├── fig5
│   ├── fig5a3_241125.R
│   ├── fig5a5_241125.R
│   ├── fig5a6_241125.R
│   ├── fig5b_241125.R
│   └── plotfgsea_241125.R
├── LICENSE
├── README.md
├── sfig1
│   ├── matchSimWithReal220106_heatmap.R
│   ├── sfig1a.R
│   ├── sfig1b.R
│   ├── sfig1c.R
│   ├── sfig1f.R
│   ├── sfig1h.R
│   ├── sfig1i.R
│   ├── sfig1j.R
│   └── sfig1k.R
├── sfig2
│   ├── sfig2b.R
│   └── sfig2c.R
├── sfig3
│   ├── sfig3a.R
│   └── sfig3b.R
└── sfig4
    ├── plotfgsea.R
    ├── sfig4a.R
    ├── sfig4b.R
    ├── sfig4c.R
    ├── sfig4d.R
    └── sfig4e.R
```


## Session information

To run the code, please ensure the following prerequisites are met before you begin:

```
> sessionInfo()
R version 4.4.1 (2024-06-14)
Platform: x86_64-apple-darwin20
Running under: macOS 15.5

Matrix products: default
BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
LAPACK: /Library/Frameworks/R.framework/Versions/4.4-x86_64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.0

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

time zone: Asia/Hong_Kong
tzcode source: internal

attached base packages:
[1] stats4    parallel  stats     graphics  grDevices utils     datasets 
[8] methods   base     

other attached packages:
 [1] qusage_2.38.0        limma_3.60.6         GSEABase_1.66.0     
 [4] graph_1.82.0         annotate_1.82.0      XML_3.99-0.18       
 [7] AnnotationDbi_1.66.0 IRanges_2.38.1       S4Vectors_0.42.1    
[10] Biobase_2.64.0       BiocGenerics_0.50.0  Seurat_5.3.0        
[13] SeuratObject_5.1.0   sp_2.2-0             AUCell_1.26.0       
[16] ggpubr_0.6.0         pheatmap_1.0.13      igraph_2.1.4        
[19] RColorBrewer_1.1-3   viridis_0.6.5        viridisLite_0.4.2   
[22] SPATA2_3.0.1         ggplot2_3.5.2        data.table_1.17.4   

loaded via a namespace (and not attached):
  [1] RcppAnnoy_0.0.22            splines_4.4.1              
  [3] later_1.4.2                 bitops_1.0-9               
  [5] R.oo_1.27.1                 tibble_3.2.1               
  [7] polyclip_1.10-7             fastDummies_1.7.5          
  [9] lifecycle_1.0.4             rstatix_0.7.2              
 [11] globals_0.18.0              lattice_0.22-7             
 [13] MASS_7.3-65                 backports_1.5.0            
 [15] magrittr_2.0.3              plotly_4.10.4              
 [17] rmarkdown_2.29              yaml_2.3.10                
 [19] httpuv_1.6.16               glmGamPoi_1.16.0           
 [21] sctransform_0.4.2           spam_2.11-1                
 [23] spatstat.sparse_3.1-0       reticulate_1.42.0          
 [25] DBI_1.2.3                   cowplot_1.1.3              
 [27] pbapply_1.7-2               lubridate_1.9.4            
 [29] multcomp_1.4-28             abind_1.4-8                
 [31] zlibbioc_1.50.0             Rtsne_0.17                 
 [33] GenomicRanges_1.56.2        R.utils_2.13.0             
 [35] purrr_1.0.4                 RCurl_1.98-1.17            
 [37] TH.data_1.1-3               sandwich_3.1-1             
 [39] GenomeInfoDbData_1.2.12     ggrepel_0.9.6              
 [41] irlba_2.3.5.1               listenv_0.9.1              
 [43] spatstat.utils_3.1-4        units_0.8-7                
 [45] goftest_1.2-3               RSpectra_0.16-2            
 [47] spatstat.random_3.4-1       fitdistrplus_1.2-2         
 [49] parallelly_1.45.0           DelayedMatrixStats_1.26.0  
 [51] codetools_0.2-20            DelayedArray_0.30.1        
 [53] tidyselect_1.2.1            UCSC.utils_1.0.0           
 [55] farver_2.1.2                SPATAData_1.0.0            
 [57] matrixStats_1.5.0           spatstat.explore_3.4-3     
 [59] jsonlite_2.0.0              Formula_1.2-5              
 [61] progressr_0.15.1            emmeans_1.11.1             
 [63] ggridges_0.5.6              survival_3.8-3             
 [65] tools_4.4.1                 ica_1.0-3                  
 [67] Rcpp_1.0.14                 glue_1.8.0                 
 [69] gridExtra_2.3               SparseArray_1.4.8          
 [71] xfun_0.52                   MatrixGenerics_1.16.0      
 [73] fftw_1.0-9                  GenomeInfoDb_1.40.1        
 [75] EBImage_4.46.0              dplyr_1.1.4                
 [77] withr_3.0.2                 fastmap_1.2.0              
 [79] digest_0.6.37               estimability_1.5.1         
 [81] timechange_0.3.0            R6_2.6.1                   
 [83] mime_0.13                   colorspace_2.1-1           
 [85] scattermore_1.2             tensor_1.5                 
 [87] RSQLite_2.4.0               jpeg_0.1-11                
 [89] spatstat.data_3.1-6         R.methodsS3_1.8.2          
 [91] tidyr_1.3.1                 generics_0.1.4             
 [93] httr_1.4.7                  htmlwidgets_1.6.4          
 [95] S4Arrays_1.4.1              uwot_0.2.3                 
 [97] pkgconfig_2.0.3             gtable_0.3.6               
 [99] rsconnect_1.4.1             blob_1.2.4                 
[101] lmtest_0.9-40               SingleCellExperiment_1.26.0
[103] XVector_0.44.0              htmltools_0.5.8.1          
[105] carData_3.0-5               dotCall64_1.2              
[107] fftwtools_0.9-11            scales_1.4.0               
[109] png_0.1-8                   spatstat.univar_3.1-3      
[111] knitr_1.50                  rstudioapi_0.17.1          
[113] tzdb_0.5.0                  reshape2_1.4.4             
[115] coda_0.19-4.1               nlme_3.1-168               
[117] cachem_1.1.0                zoo_1.8-14                 
[119] stringr_1.5.1               KernSmooth_2.23-26         
[121] miniUI_0.1.2                pillar_1.10.2              
[123] grid_4.4.1                  vctrs_0.6.5                
[125] RANN_2.6.2                  promises_1.3.3             
[127] car_3.1-3                   xtable_1.8-4               
[129] cluster_2.1.8.1             evaluate_1.0.3             
[131] readr_2.1.5                 mvtnorm_1.3-3              
[133] cli_3.6.5                   locfit_1.5-9.12            
[135] compiler_4.4.1              rlang_1.1.6                
[137] crayon_1.5.3                future.apply_1.11.3        
[139] ggsignif_0.6.4              confuns_1.0.3              
[141] labeling_0.4.3              plyr_1.8.9                 
[143] stringi_1.8.7               deldir_2.0-4               
[145] Biostrings_2.72.1           lazyeval_0.2.2             
[147] tiff_0.1-12                 spatstat.geom_3.4-1        
[149] Matrix_1.7-3                RcppHNSW_0.6.0             
[151] hms_1.1.3                   patchwork_1.3.0            
[153] bit64_4.6.0-1               sparseMatrixStats_1.16.0   
[155] future_1.58.0               statmod_1.5.0              
[157] KEGGREST_1.44.1             shiny_1.10.0               
[159] SummarizedExperiment_1.34.0 ROCR_1.0-11                
[161] memoise_2.0.1               broom_1.0.8                
[163] bit_4.6.0   
```

Notes: 

- To ensure the code runs smoothly, please confirm that all required packages are available and loaded at the beginning.

- To install the necessary packages, use `install.packages(c("package1", "package2", ...))` or `BiocManager::install(c("package1", "package2", ...))` as appropriate. This will help you build a similar development environment.


