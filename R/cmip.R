#' Clean up CMIP3 and CMIP5 meta data
#'
#' @aliases cmip5.model_id
#'
#' @export
cmip3.model_id <- function(txt) {
  # filtering the names
  # remove the first 8 characters
  txt <- substr(txt,8,nchar(txt))
  txt[txt=="bccr_bcm2_0_144.nc"] <- 			"bcc_bcm2_1"                           
  txt[txt=="sresa1b_00_cgcm3.1_t47_144.nc"] <- 		"cgcm3.1(t47)"                  
  txt[txt=="sresa1b_01_cgcm3.1_t47_144.nc"] <- 		"cgcm3.1(t47)"                
  txt[txt=="sresa1b_02_cgcm3.1_t47_144.nc"] <- 		"cgcm3.1(t47)"                  
  txt[txt=="sresa1b_03_cgcm3.1_t47_144.nc"] <- 		"cgcm3.1(t47)"                  
  txt[txt=="sresa1b_04_cgcm3.1_t47_144.nc"] <- 		"cgcm3.1(t47)"                  
  txt[txt=="sresa1b_1_cgcm3.1_t63_1850_2300_144.nc"] <- "cgcm3.1(t63)"          
  txt[txt=="cnrm_cm3_144.nc"] <- 			"cnrm_cm3"                                
  txt[txt=="csiro_mk3_0_144.nc"] <- 			"csiro_mk3"                              
  txt[txt=="csiro_mk3_5_144.nc"] <-  			"csiro_mk3_5"                             
  txt[txt=="gfdl_cm2_0_144.nc"] <-  			"gfdl_cm2_0"                               
  txt[txt=="gfdl_cm2_1_144.nc"] <- 			"gfdl_cm2"                              
  txt[txt=="giss_aom_00_144.nc"] <- 			"giss_aom"                              
  txt[txt=="giss_aom_01_144.nc"] <-  			"giss_aom"                             
  txt[txt=="GISS3.SRESA1B.00_144.nc"] <- 		"GISS3"                         
  txt[txt=="GISS3.SRESA1B.01_144.nc"] <- 		"GISS3"                         
  txt[txt=="GISS3.SRESA1B.02_144.nc"] <-  		"GISS3"                        
  txt[txt=="GISS1.SRESA1B.00_144.nc"] <-  		"GISS1"                        
  txt[txt=="GISS1.SRESA1B.01_144.nc"] <-  		"GISS1"                        
  txt[txt=="GISS1.SRESA1B.02_144.nc" ] <-  		"GISS1"                        
  txt[txt=="GISS1.SRESA1B.03_144.nc"] <-  		"GISS1-3"                        
  txt[txt=="GISS1.SRESA1B.04_144.nc"] <-   		"GISS1"                       
  txt[txt=="ingv_echam4_144.nc"] <- 		      	"ingv_echam4"                              
  txt[txt=="inmcm3_0_144.nc"] <-  			"inmcm3-1"                                
  txt[txt=="ipsl_cm4_144.nc"] <- 			"ipsl_cm4"                                 
  txt[txt=="miroc3_2_medres_00_144.nc"] <- 		"miroc3_2-medres"                       
  txt[txt=="miroc3_2_medres_01_144.nc"] <-		"miroc3_2-medres"                        
  txt[txt=="miroc3_2_medres_02_144.nc"] <- 		"miroc3_2_medres"                        
  txt[txt=="miroc3_2_hires_144.nc"] <- 			"miroc3_2_hires"                           
  txt[txt=="miub_echo_g_00_144.nc"] <-  			"miub_echo_g"                          
  txt[txt=="miub_echo_g_01_144.nc"] <- 			"miub_echo_g"                           
  txt[txt=="miub_echo_g_02_144.nc"] <- 			"miub_echo_g"                           
  txt[txt=="mpi_echam5_00_144.nc"] <-			"mpi_echam5"                          
  txt[txt=="mpi_echam5_01_144.nc"] <-  			"mpi_echam5"                         
  txt[txt=="mpi_echam5_02_144.nc"] <-  			"mpi_echam5"                         
  txt[txt=="mpi_echam5_03_144.nc"] <-  			"mpi_echam5"                       
  txt[txt=="mri_cgcm2_3_2a_00_144.nc"] <-  		"mri_cgcm2_3_2a"                       
  txt[txt=="mri_cgcm2_3_2a_01_144.nc"] <-  		"mri_cgcm2_3_2a"                        
  txt[txt=="mri_cgcm2_3_2a_02_144.nc"] <-  		"mri_cgcm2_3_2a"                       
  txt[txt=="mri_cgcm2_3_2a_03_144.nc"] <-   		"mri_cgcm2_3_2a"                       
  txt[txt=="mri_cgcm2_3_2a_04_144.nc"] <-  		"mri_cgcm2_3_2a"                     
  txt[txt=="SRESA1B_00.CCSM.atmm.2000-01_cat_2099-12_144.nc"] <- "CCSM" 
  txt[txt=="SRESA1B_01.CCSM.atmm.2000-01_cat_2099-12_144.nc"] <- "CCSM" 
  txt[txt=="SRESA1B_02.CCSM.atmm.2000-01_cat_2099-12_144.nc"] <- "CCSM" 
  txt[txt=="SRESA1B_03.CCSM.atmm.2000-01_cat_2099-12_144.nc"] <- "CCSM" 
  txt[txt=="SRESA1B_04.CCSM.atmm.2000-01_cat_2099-12_144.nc"] <- "CCSM" 
  txt[txt=="SRESA1B_05.CCSM.atmm.2000-01_cat_2099-12_144.nc"] <- "CCSM" 
  txt[txt=="SRESA1B_06.CCSM.atmm.2000-01_cat_2099-12_144.nc"] <- "CCSM" 
  txt[txt=="SRESA1B_00.PCM1_144.nc"] <- 		"PCM1"                          
  txt[txt=="SRESA1B_01.PCM1_144.nc"] <- 		"PCM1"                       
  txt[txt=="SRESA1B_02.PCM1_144.nc"] <- 		"PCM1"                          
  txt[txt=="SRESA1B_03.PCM1_144.nc"] <- 		"PCM1"                          
  txt[txt=="ukmo_hadgem1_144.nc"    ] <- 		"ukmo_hadgem1"                          
  txt[txt=="ukmo_hadcm3_144.nc"     ] <- 		"ukmo_hadcm3"
  invisible(txt)
} 

#' @export
cmip5.model_id <- function(txt) {
  txt <- substr(txt,10,nchar(txt))
  txt2 <- unlist(strsplit(txt,split="_"))
  txt <- txt2[1]
  invisible(txt)
}         
