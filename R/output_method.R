# Description: Create next generation from the 50 individuals example
# robustOCS

library("AlphaSimR")
library("tidyverse")
library("pedigreemm")
library("MASS")
library("optiSel")
library("bayestestR")
library("data.table")
library("gridExtra")
library("VennDiagram")


load("/Users/s2745345/Documents/GitHub/jortiz_robust_ocs/jortiz_robust_ocs_fork/rObjects/breedingValues.RData")

output_method <- read.table("/Users/s2745345/Documents/GitHub/jortiz_phd_review/data/output_method.txt", 
                            header = TRUE, 
                            sep = "\t")

summary(output_method[output_method$sex == "M",])

nextGenerationTable = NULL
generation = 11
year <- year + 1
nProgeny = 1


# Correlation
cor(phenotypes$EBV,phenotypes$tbv)




###### -------- Truncation ------- ######

# Males selection
  selectCriterion = output_method %>%
  dplyr::filter(sex == "M") %>%
  dplyr::filter(w_truncation > 0.0) %>% 
  dplyr::select(idd) 
  
  selSires = pop_males@id %in% selectCriterion$idd

# Females selection
  out_females = output_method %>%
  dplyr::filter(sex == "F") %>%
  dplyr::select(idd)
  
  selFemales = pop_females@id %in% out_females$idd

# Mating plan
  pop_fem = pop_females[selFemales, ]
  pop_mal = pop_males[selSires, ]
    
# Next generation

  next_gen_truncation = randCross2(females = pop_fem,
                        males = pop_mal,
                        nCrosses = nrow(out_females))
                        
  
  (realised_genetic_gain_progeny = mean(next_gen_truncation@gv)) 


  next_gen_truncation@sex = rep("MF", next_gen_truncation@nInd)
  (nextGenerationTable <- PullSumm(nextGenerationTable, c(next_gen_truncation), type = "TR_PEBV"))
  
## Construyendo la matrix A
  
  ind = next_gen_truncation@iid
  father = as.numeric(next_gen_truncation@father)
  mother = as.numeric(next_gen_truncation@mother)
  
  df_rm = data.frame(
    ind, 
    father, 
    mother
  )
  
  truncation_pedigree = rbind(pedigree, df_rm)
  
# get A
  
  ped2_truncation = pedigreemm::pedigree(sire = truncation_pedigree$father, dam = truncation_pedigree$mother, label = truncation_pedigree$ind)
  
  Amatrix = getA(ped2_truncation)
  dim(Amatrix)
  
  nSelect = 25
  
  A = as.matrix(Amatrix[(nrow(Amatrix)-nSelect+1):nrow(Amatrix), (ncol(Amatrix)-nSelect+1):ncol(Amatrix)]) 
  dim(A)
  
  (group_coancestry_progeny_truncation = mean(A))
  

# Group coancestry  
  
  group_coancestry_sel_truncation = 0.3257680763244627
  (deltaC = (group_coancestry_progeny_truncation - group_coancestry_sel_truncation) / (1.0 - group_coancestry_sel_truncation))

# Effective population size  
  
  (Ne = 1 / (2 * deltaC))
  

###### -------- Conservation ------- ######
  
# Male selection
  sireID = output_method %>%
    dplyr::filter(sex == "M") %>%
    dplyr::filter(w_conservation > 0.0) %>% 
    dplyr::select(idd) 
  
  sireOC = output_method %>%
    dplyr::filter(sex == "M") %>%
    dplyr::select(w_conservation)
    
  
  (sires = rep(x = sireID[[1]], times = floor(2*25*(sireOC[[1]]))))
  
# Female selection
  
  damID = output_method %>%
    dplyr::filter(sex == "F") %>%
    dplyr::select(idd)
  
  damOC = output_method %>%
    dplyr::filter(sex == "F") %>%
    dplyr::select(w_conservation)
  
  females = pop_females@id %in% out_females$idd
  
  (dams = rep(x = damID[[1]], times = floor(2*25*(damOC[[1]]))))
  
# Mating plan

  (matingPlan = cbind(as.character(dams), 
                     as.character(sample(sireID[[1]], size = length(dams), replace = TRUE, prob = sireOC[[1]]))))
  
# Next generation
  
  next_gen_conservation = makeCross2(females = pop_females,
                                     males = pop_males,
                                     crossPlan = matingPlan)
                                   
  
  (realised_genetic_gain_progeny = mean(next_gen_conservation@gv)) 
  
  next_gen_conservation@sex = rep("MF", next_gen_conservation@nInd)
  (nextGenerationTable <- PullSumm(nextGenerationTable, c(next_gen_conservation), type = "TR_PEBV"))
  
  ## Construyendo la matrix A
  
  ind = next_gen_conservation@iid
  father = as.numeric(next_gen_conservation@father)
  mother = as.numeric(next_gen_conservation@mother)
  
  df_rm = data.frame(
    ind, 
    father, 
    mother
  )
  
  conservation_pedigree = rbind(pedigree, df_rm)
  
  # get A
  
  ped2_conservation = pedigreemm::pedigree(sire = conservation_pedigree$father, dam = conservation_pedigree$mother, label = truncation_pedigree$ind)
  
  Amatrix = getA(ped2_conservation)
  dim(Amatrix)
  
  nSelect = 25
  
  A = as.matrix(Amatrix[(nrow(Amatrix)-nSelect+1):nrow(Amatrix), (ncol(Amatrix)-nSelect+1):ncol(Amatrix)]) 
  dim(A)
  
  (group_coancestry_progeny_conservation = mean(A))
  
  # Group coancestry  
  
  group_coancestry_sel_conservation = 0.0885547903075112
  (deltaC = (group_coancestry_progeny_conservation - group_coancestry_sel_conservation) / (1.0 - group_coancestry_sel_conservation))
  
  # Effective population size  
  
  (Ne = 1 / (2 * deltaC))
  
  ###### -------- Uncertainty ------- ######
  
  # Male selection
  sireID = output_method %>%
    dplyr::filter(sex == "M") %>%
    dplyr::filter(w_uncertainty > 0.0) %>% 
    dplyr::select(idd) 
  
  sireOC = output_method %>%
    dplyr::filter(sex == "M") %>%
    dplyr::select(w_uncertainty)
  
  
  (sires = rep(x = sireID[[1]], times = floor(2*25*(sireOC[[1]]))))
  
  # Female selection
  
  damID = output_method %>%
    dplyr::filter(sex == "F") %>%
    dplyr::select(idd)
  
  damOC = output_method %>%
    dplyr::filter(sex == "F") %>%
    dplyr::select(w_uncertainty)
  
  females = pop_females@id %in% out_females$idd
  
  (dams = rep(x = damID[[1]], times = floor(2*25*(damOC[[1]]))))
  
  # Mating plan
  
  (matingPlan = cbind(as.character(dams), 
                      as.character(sample(sireID[[1]], size = length(dams), replace = TRUE, prob = sireOC[[1]]))))
  
  # Next generation
  
  next_gen_uncertainty = makeCross2(females = pop_females,
                                     males = pop_males,
                                     crossPlan = matingPlan)
  
  
  (realised_genetic_gain_progeny = mean(next_gen_uncertainty@gv)) 
  
  next_gen_uncertainty@sex = rep("MF", next_gen_conservation@nInd)
  (nextGenerationTable <- PullSumm(nextGenerationTable, c(next_gen_uncertainty), type = "TR_PEBV"))
  
  ## Construyendo la matrix A
  
  ind = next_gen_uncertainty@iid
  father = as.numeric(next_gen_uncertainty@father)
  mother = as.numeric(next_gen_uncertainty@mother)
  
  df_rm = data.frame(
    ind, 
    father, 
    mother
  )
  
  uncertainty_pedigree = rbind(pedigree, df_rm)
  
  # get A
  
  ped2_uncertainty = pedigreemm::pedigree(sire = uncertainty_pedigree$father, dam = uncertainty_pedigree$mother, label = uncertainty_pedigree$ind)
  
  Amatrix = getA(ped2_uncertainty)
  dim(Amatrix)
  
  nSelect = 25
  
  A = as.matrix(Amatrix[(nrow(Amatrix)-nSelect+1):nrow(Amatrix), (ncol(Amatrix)-nSelect+1):ncol(Amatrix)]) 
  dim(A)
  
  (group_coancestry_progeny_uncertainty = mean(A))
  
  # Group coancestry  
  
  group_coancestry_sel_uncertainty = 0.0893782291869913
  (deltaC = (group_coancestry_progeny_uncertainty - group_coancestry_sel_uncertainty) / (1.0 - group_coancestry_sel_uncertainty))
  
  # Effective population size  
  
  (Ne = 1 / (2 * deltaC))
  
  
  
  
  
  