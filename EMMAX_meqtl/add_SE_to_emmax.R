
args <- commandArgs(trailingOnly = TRUE)
j <- as.numeric(args[1])
pop <- as.character(args[2])


if(pop == "Baka") {
best_covs <- read.delim(paste0("/share/hennlab/users/glmeeks/age_methylation/EMMAX_meqtl/fake_perms/", j, "beta_fake_meqtl_regress_Baka.ps"), header=FALSE)
print(paste0("/share/hennlab/users/glmeeks/age_methylation/EMMAX_meqtl/fake_perms/fake_perms/", j, "beta_fake_meqtl_regress_Baka.ps"))
#best_covs <- read.delim("Baka_sex_methPC23_Neutro_genoPC12345_meqtl_regressed.ps", header=FALSE)
#best_covs <- read.delim("Baka_sex_methPC23_Neutro_genoPC12345_fusion_regressed.ps", header=FALSE)

n = 38
}
    
if(pop == "KHS") {
best_covs <- read.delim(paste0("/share/hennlab/users/glmeeks/age_methylation/EMMAX_meqtl/fake_perms/", j, "beta_fake_meqtl_regress_KHS.ps"), header=FALSE)
print(paste0("/share/hennlab/users/glmeeks/age_methylation/EMMAX_meqtl/fake_perms/fake_perms/", j, "beta_fake_meqtl_regress_KHS.ps"))
#best_covs <- read.delim("KHS_sex_batch_PC13_meth_Neutro_genoPCs_meqtl_regressed.ps", header=FALSE)
#best_covs <- read.delim("KHS_sex_batch_PC13_meth_Neutro_genoPCs_fusion_regressed.ps", header=FALSE)    
n = 52
}
    
if(pop == "Himba"){
best_covs <- read.delim(paste0("/share/hennlab/users/glmeeks/age_methylation/EMMAX_meqtl/fake_perms/", j, "beta_fake_meqtl_regress_Himba.ps"), header=FALSE)
print(paste0("/share/hennlab/users/glmeeks/age_methylation/EMMAX_meqtl/fake_perms/fake_perms/", j, "beta_fake_meqtl_regress_Himba.ps"))
#best_covs <- read.delim("Himba_sex_cpPCs_Neutro_genoPCs_meqtl_regressed.ps", header=FALSE)
#best_covs <- read.delim("Himba_sex_cpPCs_Neutro_genoPCs_fusion_regressed.ps", header=FALSE)    
n = 51
}


# if(pop == "all"){
# best_covs <- read.delim(paste0("fake_perms/", j, "all_covs_meqtl_regressed.ps "), header=FALSE)}

colnames(best_covs)<-c("CPG.Labels","V2","P.value")


#example: se.from.p(effect.size = 0.71, p = 0.013, N = 75,effect.size.type= "difference", calculate.g = TRUE)

##Heterogeneity of effects
se.from.p = function(effect.size, p, N, effect.size.type = "difference", calculate.g = FALSE) {
# Define helper funcs
    sssbc = function(totaln){return(1 - (3/(4 * totaln - 9)))}
    hedges_g = function(d, totaln){mapply(function(.x, .y) .x * sssbc(.y), d, totaln)}
# Set params
    ES = effect.size
    p = p
    N = N
    ES.type = effect.size.type
    calculate.g = calculate.g
if (is.numeric(ES) == FALSE) {
        stop("'effect.size' is not of type numeric().")
    }
if (is.numeric(p) == FALSE) {
        stop("'p' is not of type numeric().")
    }
if (sum(p < 0 | p > 1) > 0){
        stop("values of 'p' must range between 0 and 1.")
    }
if (is.numeric(N) == FALSE) {
        stop("'N' is not of type numeric().")
    }
if (ES.type %in% c("difference", "ratio") == FALSE) {
        stop("'effect.size.type' must be either 'difference' or 'ratio'.")
    }
if (ES.type == "ratio" & sum(ES < 0) > 0) {
        stop("when 'effect.size.type' is 'ratio', values of 'effect.size' must be equal to or greater than 0.")
    }
# Difference vs. Ratio
# Difference
    if (ES.type == "difference") {
        if (calculate.g == TRUE) {
            ES = hedges_g(d = ES, totaln = N)
            z = -0.862 + sqrt(0.743 - 2.404 * log(p))
            SE = ES/z
            SD = SE * sqrt(N)
            LLCI = ES - qnorm(0.975) * SE
            ULCI = ES + qnorm(0.975) * SE
            data = data.frame(ES, abs(SE), abs(SD), LLCI, ULCI)
            colnames(data) = c("Hedges.g", "StandardError", "StandardDeviation", "LLCI", "ULCI")

        } else {
            z = -0.862 + sqrt(0.743 - 2.404 * log(p))
            SE = ES/z
            SD = SE * sqrt(N)
            LLCI = ES - qnorm(0.975) * SE
            ULCI = ES + qnorm(0.975) * SE
            data = data.frame(ES, abs(SE), abs(SD), LLCI, ULCI)
            colnames(data) = c("EffectSize", "StandardError", "StandardDeviation", "LLCI", "ULCI")
        }
    }

    if (ES.type == "ratio") {
        if (calculate.g == TRUE) {
            stop("Hedges' g cannot be calculated for ratios using this function; set 'calculate.g=FALSE'.")
        } else {
            z = -0.862 + sqrt(0.743 - 2.404 * log(p))
            ES = log(ES)
            SE = abs(ES/z)
            SD = SE * sqrt(N)
            LLCI = ES - qnorm(0.975) * SE
            ULCI = ES + qnorm(0.975) * SE
            # Exponentiate to get original scale
            expES = exp(ES)
            expLLCI = exp(LLCI)
            expULCI = exp(ULCI)
            data = data.frame(ES, abs(SE), abs(SD), LLCI, ULCI, expES, expLLCI, expULCI)
            colnames(data) = c("logEffectSize", "logStandardError", "logStandardDeviation", "logLLCI", "logULCI", "EffectSize",
                "LLCI", "ULCI")
        }
    }
    #return(data)
    return(abs(SE))                                      
}                                       
best_covs$SE <- apply(best_covs, 1, function(x) se.from.p(effect.size=as.numeric(x[2]), p=as.numeric(x[3]), N=n, effect.size.type="difference", calculate.g=FALSE))   
#save(list=c("best_covs"),file=paste0(pop, "best_covs_with_SE_meqtl_regressed.RData"))      
#save(list=c("best_covs"),file=paste0(pop, "best_covs_with_SE_fusion_regressed.RData"))      
save(list=c("best_covs"),file=paste0("fake_perms/", pop, j, "best_covs_with_SE_fake_meqtl_regressed.RData"))                         

