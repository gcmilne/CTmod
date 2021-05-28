############################################
# Diagnostic sensitivities & specificities #
############################################

## Unless otherwise specified, data are from the following ref:
## https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8078101/

## Agglutination assay (AA) ##

#sensitivity
se_aa <- c(0.988,  #Pastorex
           0.968,  #Toxocell
           0.937   #Toxolatex
) 

#specificity
sp_aa <- c(0.988,  #Pastorex
           0.976,  #Toxocell
           0.971   #Toxolatex
) 

#sample size weights
n_aa <- c(589,  #Pastorex
          589,  #Toxocell
          589   #Toxolatex
)

#weighted mean sensitivity
se_aa <- weighted.mean(se_aa, weights = n_aa)

#weighted mean specificity
sp_aa <- weighted.mean(sp_aa, weights = n_aa)

#sum sample size
n_aa <- sum(n_aa)


## Chemiluminescence assay (CLIA) ##

#sensitivity
se_clia <- 0.934  #Vitros ECiQ Toxoplasma IgG

#specificity
sp_clia <- 1      #Vitros ECiQ Toxoplasma IgG

#sample size
n_clia <- 719


## Electrochemiluminescence assay (ECLIA) ##

#sensitivity
se_eclia <- 0.986  #Elecsys Toxo IgG

#specificity
sp_eclia <- 0.996  #Elecsys Toxo IgG

#sample size
n_eclia <- 2214


## Enzyme-linked fluorescence assay (ELFA) ##

#sensitivity
se_elfa <- 0.955  #Vidas Toxo IgG

#specificity
sp_elfa <- 0.998  #Vidas Toxo IgG

#sample size
n_elfa  <- 3368


## Enzyme-linked immunofluorescence assay (ELISA) ##

#sensitivity
se_elisa <- 0.964  #Platelia Toxo IgG (microplate ELISA)

#specificity
sp_elisa <- 0.994  #Platelia Toxo IgG (microplate ELISA)

#sample size
n_elisa  <- 1436


## Haemagglutination assay (HA) ##

#sensitivity
se_ha <- 1      #Toxo HAI

#specificity
sp_ha <- 0.992  #Toxo HAI

#sample size
n_ha  <- 589


## Indirect immunofluoresence assay (IFA) ##

#sensitivity
se_ifa <- c(0.9304,  #ref = https://www.sciencedirect.com/science/article/pii/S0732889305001033
            0.909,
            1
)

#specificity
sp_ifa <- c(0.9945,  #ref = https://www.sciencedirect.com/science/article/pii/S0732889305001033
            1, 
            0.98
)

#sample size
n_ifa  <- c(656,    #ref = https://www.sciencedirect.com/science/article/pii/S0732889305001033
            NA,
            NA
)

#mean sensitivity
se_ifa <- mean(se_ifa)

#mean specificity
sp_ifa <- mean(sp_ifa)


## Microparticle enzyme immunoassay (MEIA) ##

#sensitivity
se_meia <- c(0.897,  #Access Toxo IgG II
             0.907,  #Architect Toxo IgG
             0.961   #Axsym Toxo IgG
)

#specificity
sp_meia <- c(1,      #Access Toxo IgG II
             0.998,  #Architect Toxo IgG
             0.997   #Axsym Toxo IgG
)

#sample size weights
n_meia  <- c(760,    #Access Toxo IgG II
             2992,   #Architect Toxo IgG
             1555    #Axsym Toxo IgG
)

#weighted mean sensitivity
se_meia <- weighted.mean(se_meia, weights = n_meia)

#weighted mean specificity
sp_meia <- weighted.mean(sp_meia, weights = n_meia)

#sum sample size
n_meia <- sum(n_meia)


## Western blot (WB) ##

#sensitivity
se_wb <- 1  #WB Toxo IgGII

#specificity
sp_wb <- 1  #WB Toxo IgGII

#sample size
n_wb <- 569


#########################################################
## Are values significantly different between methods? ##
#########################################################

##sensitivity
se_values <- data.frame("se" = c(se_aa, se_clia, se_eclia, se_elfa, se_elisa, se_ha, se_meia, se_wb), 
                        method = c("AA", "CLIA", "ECLIA", "ELFA", "ELISA", "HA", "MEIA", "WB"), 
                        n = c(n_aa, n_clia, n_eclia, n_elfa, n_elisa, n_ha, n_meia, n_wb))

#plot
# par(mfrow=c(2,1))
# boxplot(se_values$se ~ se_values$method, xlab = "Method", ylab = "Sensitivity")

#logistic regression
# se_mod <- glm(se ~ method, family=binomial(logit), data=se_values, weights = n)
# summary(se_mod)

##specificity
sp_values <- data.frame("sp" = c(sp_aa, sp_clia, sp_eclia, sp_elfa, sp_elisa, sp_ha, sp_meia, sp_wb), 
                        method = se_values$method, 
                        n = c(n_aa, n_clia, n_eclia, n_elfa, n_elisa, n_ha, n_meia, n_wb))

#plot
# boxplot(sp_values$sp ~ sp_values$method, xlab = "Method", ylab = "Specificity")

#logistic regression
# sp_mod <- glm(sp ~ method, family=binomial(logit), data=sp_values, weights = n)
# summary(sp_mod)


# rm(c(se_mod, sp_mod))
