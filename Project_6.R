# Identifying genetic markers associated with Chronic Fatigue Syndrome (CFS)

# Aliaksei Shaukunou
# Mateusz Zacharecki
# Julia Girtler


##################################################################

### 1. Descriptive analysis

# Installing packages (!don't run if already installed!)

install.packages("dplyr")
install.packages("pROC")
install.packages("OptimalCutpoints")
install.packages('HardyWeinberg')

# Importing libraries

library(dplyr)
library(pROC)
library(OptimalCutpoints)
library('HardyWeinberg')

# Importing data

data <- read.csv("https://raw.githubusercontent.com/immune-stats/Biostatistics_PW_2025_2026/main/project_6/data_project_6.csv")

head(data)
str(data)

# Splitting dataset based on group of examined patients

hc_data <- data %>% filter(Group == "HC")
cfs_data <- data %>% filter(Group == "CFS")

# Dimensions and NAs

cat("Dimensions:", dim(data), "\n")
cat("Total NAs:", sum(is.na(data)), "\n\n")
cat("NAs by Column:\n")
print(colSums(is.na(data)))

# Table frequencies

cat("Group distribution:\n")
table(data$Group, useNA = 'ifany')
cat("\nGender distribution:\n")
table(data$Gender, useNA = 'ifany')
cat("\nrs2476601 distribution:\n")
table(data$rs2476601, useNA = 'ifany')
cat("\nrs3087243 distribution:\n")
table(data$rs3087243, useNA = 'ifany')
cat("\nrs3807306 distribution:\n")
table(data$rs3807306, useNA = 'ifany')
cat("\nrs1800629 distribution:\n")
table(data$rs1800629, useNA = 'ifany')
cat("\nrs1799724 distribution:\n")
table(data$rs1799724, useNA = 'ifany')

# Plots of frequencies

par(mfrow = c(3, 3), mar = c(4, 4, 2, 1))

for (col in names(data)) {
  counts <- table(data[[col]], useNA = "ifany")
  names(counts)[is.na(names(counts))] <- "NA"
  
  bp <- barplot(counts,
                main = paste(col, "distribution"),
                ylab = "Count",
                ylim = c(0, max(counts) * 1.2),
                col = "#b8dee8")
  
  text(x = bp,
       y = counts,
       labels = counts,
       pos = 3)
}


##################################################################

### 2. Information for each genetic marker from the 1000HGP

par(mfrow = c(1, 1))

### **rs2476601**
  
# Chromosome: **1**
#   
#   Location: **113834946 (GRCh38.p14)**
#   
#   Alleles:
#   
#   
#   *   Referential Allele: **A**
#   *   Alternative Allele: **G**
#   
#   
#   
#   
#   
#   
#   
#   | Population | Ref Allele | Alt Allele
# | :--- | :--- | :--- |
#   | **Global** | A=0.0274 | G=0.9726 |
#   | African | A=0.0030 | G=0.9970 |
#   | East Asian | A=0.0000 | G=1.0000 |
#   | Europe | A=0.0944 | G=0.9056 |
#   | South Asian | A=0.013 | G=0.987 |
#   | American | A=0.036 | G=0.964 |

populations <- c("Global", "African", "East Asian", "Europe", "South Asian", "American")
A <- c(0.0274, 0.0030, 0.0000, 0.0944, 0.013, 0.036)
G <- c(0.9726, 0.9970, 1.0000, 0.9056, 0.987, 0.964)

data_matrix <- rbind(A = A, G = G)
colnames(data_matrix) <- populations

barplot(data_matrix,
        main = "Allele Distribution of rs2476601",
        xlab = "Population",
        ylab = "Frequency",
        col = c("#8ed087", "#f4d059"),
        cex.names = 0.8,
        legend.text = TRUE,
        args.legend = list(x = "topright", bty = "n", inset = c(-0.05, -0.15)))

### **rs3087243**

# Chromosome: **2**
#   
#   Location: **203874196 (GRCh38.p14)**
#   
#   Alleles:
#   
#   
#   *   Referential Allele: **G**
#   *   Alternative Allele: **A**
#   
#   
#   
#   
#   
#   
#   
#   | Population | Ref Allele | Alt Allele
# | :--- | :--- | :--- |
#   | **Global** | G=0.6310 | A=0.3690 |
#   | African | G=0.8238 | A=0.1762 |
#   | East Asian | G=0.7371 | A=0.2629 |
#   | Europe | G=0.5298 | A=0.4702 |
#   | South Asian | G=0.373 | A=0.627 |
#   | American | G=0.620 | A=0.380 |

populations <- c("Global", "African", "East Asian", "Europe", "South Asian", "American")
G <- c(0.6310, 0.8238, 0.7371, 0.5298, 0.373, 0.620)
A <- c(0.3690, 0.1762, 0.2629, 0.4702, 0.627, 0.380)

data_matrix <- rbind(G = G, A = A)
colnames(data_matrix) <- populations

barplot(data_matrix,
        main = "Allele Distribution of rs3087243",
        xlab = "Population",
        ylab = "Frequency",
        col = c("#8ed087", "#f4d059"),
        cex.names = 0.8,
        legend.text = TRUE,
        args.legend = list(x = "topright", bty = "n", inset = c(-0.05, -0.15)))

### **rs3807306**

# Chromosome: **7**
#   
#   Location: **128940626 (GRCh38.p14)**
#   
#   Alleles:
#   
#   
#   *   Referential Allele: **G**
#   *   Alternative Allele: **T**
#   
#   
#   
#   
#   
#   
#   
#   | Population | Ref Allele | Alt Allele
# | :--- | :--- | :--- |
#   | **Global** | G=0.6540 | T=0.3460 |
#   | African | G=0.7035 | T=0.2965 |
#   | East Asian | G=0.8472 | T=0.1528 |
#   | Europe | G=0.4871 | T=0.5129 |
#   | South Asian | G=0.553 | T=0.447 |
#   | American | G=0.663 | T=0.337 |

### **rs3807306**

populations <- c("Global", "African", "East Asian", "Europe", "South Asian", "American")
G <- c(0.6540, 0.7035, 0.8472, 0.4871, 0.553, 0.663)
T <- c(0.3460, 0.2965, 0.1528, 0.5129, 0.447, 0.337)

data_matrix <- rbind(G = G, T = T)
colnames(data_matrix) <- populations

barplot(data_matrix,
        main = "Allele Distribution of rs3807306",
        xlab = "Population",
        ylab = "Frequency",
        col = c("#8ed087", "#f4d059"),
        cex.names = 0.8,
        legend.text = TRUE,
        args.legend = list(x = "topright", bty = "n", inset = c(-0.05, -0.15)))

### **rs1800629**

# Chromosome: **6**
#   
#   Location: **31575254 (GRCh38.p14)**
#   
#   Alleles:
#   
#   
#   *   Referential Allele: **G**
#   *   Alternative Allele: **A**
#   
#   
#   
#   
#   
#   
#   
#   | Population | Ref Allele | Alt Allele
# | :--- | :--- | :--- |
#   | **Global** | G=0.9097 | A=0.0903 |
#   | African | G=0.8805 | A=0.1195 |
#   | East Asian | G=0.9415 | A=0.0585 |
#   | Europe | G=0.8658 | A=0.1342 |
#   | South Asian | G=0.947 | A=0.053 |
#   | American | G=0.931 | A=0.069 |

populations <- c("Global", "African", "East Asian", "Europe", "South Asian", "American")
G <- c(0.9097, 0.8805, 0.9415, 0.8658, 0.947, 0.931)
A <- c(0.0903, 0.1195, 0.0585, 0.1342, 0.053, 0.069)

data_matrix <- rbind(G = G, A = A)
colnames(data_matrix) <- populations

barplot(data_matrix,
        main = "Allele Distribution of rs1800629",
        xlab = "Population",
        ylab = "Frequency",
        col = c("#8ed087", "#f4d059"),
        cex.names = 0.8,
        legend.text = TRUE,
        args.legend = list(x = "topright", bty = "n", inset = c(-0.05, -0.15)))

### **rs1799724**

# Chromosome: **6**
#   
#   Location: **31574705 (GRCh38.p14)**
#   
#   Alleles:
#   
#   
#   *   Referential Allele: **C**
#   *   Alternative Allele: **T**
#   
#   
#   
#   
#   
#   
#   
#   | Population | Ref Allele | Alt Allele
# | :--- | :--- | :--- |
#   | **Global** | C=0.9010 | T=0.0990 |
#   | African | C=0.9758 | T=0.0242 |
#   | East Asian | C=0.8750 | T=0.1250 |
#   | Europe | C=0.9056 | T=0.0944 |
#   | South Asian | C=0.881 | T=0.119 |
#   | American | C=0.817 | T=0.183 |

populations <- c("Global", "African", "East Asian", "Europe", "South Asian", "American")
C <- c(0.9010, 0.9758, 0.8750, 0.9056, 0.881, 0.817)
T <- c(0.0990, 0.0242, 0.1250, 0.0944, 0.119, 0.183)

data_matrix <- rbind(C = C, T = T)
colnames(data_matrix) <- populations

barplot(data_matrix,
        main = "Allele Distribution of rs1799724",
        xlab = "Population",
        ylab = "Frequency",
        col = c("#8ed087", "#f4d059"),
        cex.names = 0.8,
        legend.text = TRUE,
        args.legend = list(x = "topright", bty = "n", inset = c(-0.05, -0.15)))


##################################################################

### 3. Allele distributions of each biomarker for HC and comparison against distributions from the 1000HGP

# rs2476601

A <- 2 * 0 + table(hc_data$rs2476601, useNA = 'ifany')[["AG"]]
G <- 2 * table(hc_data$rs2476601, useNA = 'ifany')[["GG"]] + table(hc_data$rs2476601, useNA = 'ifany')[["AG"]]

cat("Allele distribution:\n")
as.table(c(A = A, G = G))

cat("\nComparison of distribution:\n")

cat("\nAfrican:\n")
binom.test(A, A + G, p = 0.0030)

cat("\nEast Asian:\n")
binom.test(A, A + G, p = 0.0000)

cat("\nEurope:\n")
binom.test(A, A + G, p = 0.0944)

cat("\nSouth Asian:\n")
binom.test(A, A + G, p = 0.013)

cat("\nAmerican:\n")
binom.test(A, A + G, p = 0.036)

# rs3087243

G <- 2 * table(hc_data$rs3087243, useNA = 'ifany')[["GG"]] + table(hc_data$rs3087243, useNA = 'ifany')[["AG"]]
A <- 2 * table(hc_data$rs3087243, useNA = 'ifany')[["AA"]] + table(hc_data$rs3087243, useNA = 'ifany')[["AG"]]

cat("Allele distribution:\n")
as.table(c(G = G, A = A))

cat("\nComparison of distribution:\n")

cat("\nAfrican:\n")
binom.test(G, G + A, p = 0.8238)

cat("\nEast Asian:\n")
binom.test(G, G + A, p = 0.7371)

cat("\nEurope:\n")
binom.test(G, G + A, p = 0.5298)

cat("\nSouth Asian:\n")
binom.test(G, G + A, p = 0.373)

cat("\nAmerican:\n")
binom.test(G, G + A, p = 0.620)

# rs3807306

G <- 2 * table(hc_data$rs3807306, useNA = 'ifany')[["GG"]] + table(hc_data$rs3807306, useNA = 'ifany')[["GT"]]
T <- 2 * table(hc_data$rs3807306, useNA = 'ifany')[["TT"]] + table(hc_data$rs3807306, useNA = 'ifany')[["GT"]]

cat("Allele distribution:\n")
as.table(c(G = G, T = T))

cat("\nComparison of distribution:\n")

cat("\nAfrican:\n")
binom.test(G, G + T, p = 0.7035)

cat("\nEast Asian:\n")
binom.test(G, G + T, p = 0.8472)

cat("\nEurope:\n")
binom.test(G, G + T, p = 0.4871)

cat("\nSouth Asian:\n")
binom.test(G, G + T, p = 0.553)

cat("\nAmerican:\n")
binom.test(G, G + T, p = 0.663)

# rs1800629

G <- 2 * table(hc_data$rs1800629, useNA = 'ifany')[["GG"]] + table(hc_data$rs1800629, useNA = 'ifany')[["AG"]]
A <- 2 * table(hc_data$rs1800629, useNA = 'ifany')[["AA"]] + table(hc_data$rs1800629, useNA = 'ifany')[["AG"]]

cat("Allele distribution:\n")
as.table(c(G = G, A = A))

cat("\nComparison of distribution:\n")

cat("\nAfrican:\n")
binom.test(G, G + A, p = 0.8805)

cat("\nEast Asian:\n")
binom.test(G, G + A, p = 0.9415)

cat("\nEurope:\n")
binom.test(G, G + A, p = 0.8658)

cat("\nSouth Asian:\n")
binom.test(G, G + A, p = 0.947)

cat("\nAmerican:\n")
binom.test(G, G + A, p = 0.931)

# rs1799724

C <- 2 * table(hc_data$rs1799724, useNA = 'ifany')[["CC"]] + table(hc_data$rs1799724, useNA = 'ifany')[["CT"]]
T <- 2 * table(hc_data$rs1799724, useNA = 'ifany')[["TT"]] + table(hc_data$rs1799724, useNA = 'ifany')[["CT"]]

cat("Allele distribution:\n")
as.table(c(C = C, T = T))

cat("\nComparison of distribution:\n")

cat("\nAfrican:\n")
binom.test(C, C + T, p = 0.9758)

cat("\nEast Asian:\n")
binom.test(C, C + T, p = 0.8750)

cat("\nEurope:\n")
binom.test(C, C + T, p = 0.9056)

cat("\nSouth Asian:\n")
binom.test(C, C + T, p = 0.881)

cat("\nAmerican:\n")
binom.test(C, C + T, p = 0.817)

cat('European population is the closets one to the HC, the closest allele is rs3087243')
# rs3087243: HC allele % vs 1000HGP allele % (A and G)

G <- 2 * table(hc_data$rs3087243, useNA = "no")[["GG"]] + table(hc_data$rs3087243, useNA = "no")[["AG"]]
A <- 2 * table(hc_data$rs3087243, useNA = "no")[["AA"]] + table(hc_data$rs3087243, useNA = "no")[["AG"]]

total <- A + G

hc_pct   <- c(A = 100 * A / total, G = 100 * G / total)
hgp_pct  <- c(A = 100 * 0.4702,     G = 100 * 0.5298)

mat <- rbind(HC = hc_pct,
             `1000HGP` = hgp_pct)
bp <- barplot(mat,
              beside = TRUE,
              col = c("red", "#5170ff"),
              ylab = "Allele frequency (%)",
              ylim = c(0, 70),
              main = "rs3087243 allele distribution: control group vs European population")

legend("topright",
       legend = c("HC", "1000HGP"),
       fill = c("red", "#5170ff"),
       bty = "n")

text(x = bp, y = mat, labels = paste0(sprintf("%.1f", mat), ' %'), pos = 3)


##################################################################

### 4. Genotype distributions of each biomarker for both HC and CFS and Hardy-Weinberg Equilibrium

# Table frequencies for Healthy Control

cat("\nrs2476601 distribution for Healthy Control:\n")
table(hc_data$rs2476601, useNA = 'ifany')
cat("\nrs3087243 distribution for Healthy Control:\n")
table(hc_data$rs3087243, useNA = 'ifany')
cat("\nrs3807306 distribution for Healthy Control:\n")
table(hc_data$rs3807306, useNA = 'ifany')
cat("\nrs1800629 distribution for Healthy Control:\n")
table(hc_data$rs1800629, useNA = 'ifany')
cat("\nrs1799724 distribution for Healthy Control:\n")
table(hc_data$rs1799724, useNA = 'ifany')

# Testing Hardy-Weinberg Equilibrium for Healthy Control using exact test

identificators <- c('rs2476601', 'rs3087243', 'rs3807306', 'rs1800629', 'rs1799724')
levels <- list('rs2476601' = c('AA', 'AG', 'GG'),
               'rs3087243' = c('AA', 'AG', 'GG'),
               'rs3807306' = c('GG', 'GT', 'TT'),
               'rs1800629' = c('AA', 'AG', 'GG'),
               'rs1799724' = c('CC', 'CT', 'TT'))
for(rs in identificators) {
  t <- factor(hc_data[[rs]], levels = levels[[rs]])
  t <- table(t)
  test <- HWExact(t, verbose = F)
  pval <- test$pval
  cat("\nExact test on HC's ", rs, " biomarker:\n")
  cat("p-value: ", pval, '\n')
  cat("\n--------------------------------------------\n")
}

# Testing Hardy-Weinberg Equilibrium for Healthy Control

cat("\nchi-square test on HC's rs2476601 biomarker:\n")

AA <- 0
AG <- table(hc_data$rs2476601, useNA = 'ifany')[['AG']]
GG <- table(hc_data$rs2476601, useNA = 'ifany')[['GG']]

pi <- (2 * AA + AG) / (2 * (AA + AG + GG))

cat("p-value: ", chisq.test(table(factor(hc_data$rs2476601, levels = c('AA', 'AG', 'GG'))),
                            p = c(pi**2, 2 * pi * (1 - pi), (1 - pi)**2))$p.value, "\n")

cat("\n--------------------------------------------\n")
cat("\nchi-square test on HC's rs3087243 biomarker:\n")

AA <- table(hc_data$rs3087243, useNA = 'ifany')[['AA']]
AG <- table(hc_data$rs3087243, useNA = 'ifany')[['AG']]
GG <- table(hc_data$rs3087243, useNA = 'ifany')[['GG']]

pi <- (2 * AA + AG) / (2 * (AA + AG + GG))

cat("p-value: ", chisq.test(table(factor(hc_data$rs3087243, levels = c('AA', 'AG', 'GG'))),
                            p = c(pi**2, 2 * pi * (1 - pi), (1 - pi)**2))$p.value, "\n")

cat("\n--------------------------------------------\n")
cat("\nchi-square test on HC's rs3807306 biomarker:\n")

GG <- table(hc_data$rs3807306, useNA = 'ifany')[['GG']]
GT <- table(hc_data$rs3807306, useNA = 'ifany')[['GT']]
TT <- table(hc_data$rs3807306, useNA = 'ifany')[['TT']]

pi <- (2 * GG + GT) / (2 * (GG + GT + TT))

cat("p-value: ", chisq.test(table(factor(hc_data$rs3807306, levels = c('GG', 'GT', 'TT'))),
                            p = c(pi**2, 2 * pi * (1 - pi), (1 - pi)**2))$p.value, "\n")

cat("\n--------------------------------------------\n")
cat("\nchi-square test on HC's rs1800629 biomarker:\n")

AA <- table(hc_data$rs1800629, useNA = 'ifany')[['AA']]
AG <- table(hc_data$rs1800629, useNA = 'ifany')[['AG']]
GG <- table(hc_data$rs1800629, useNA = 'ifany')[['GG']]

pi <- (2 * AA + AG) / (2 * (AA + AG + GG))

cat("p-value: ", chisq.test(table(factor(hc_data$rs1800629, levels = c('AA', 'AG', 'GG'))),
                            p = c(pi**2, 2 * pi * (1 - pi), (1 - pi)**2))$p.value, "\n")

cat("\n--------------------------------------------\n")
cat("\nchi-square test on HC's rs1799724 biomarker:\n")

CC <- table(hc_data$rs1799724, useNA = 'ifany')[['CC']]
CT <- table(hc_data$rs1799724, useNA = 'ifany')[['CT']]
TT <- table(hc_data$rs1799724, useNA = 'ifany')[['TT']]

pi <- (2 * CC + CT) / (2 * (CC + CT + TT))

cat("p-value: ", chisq.test(table(factor(hc_data$rs1799724, levels = c('CC', 'CT', 'TT'))),
                            p = c(pi**2, 2 * pi * (1 - pi), (1 - pi)**2))$p.value, "\n")

# Table frequencies for patients with Chronic Fatigue Syndrome with an infection trigger at their disease onset

cat("\nrs2476601 distribution for patients with Chronic Fatigue Syndrome:\n")
table(cfs_data$rs2476601, useNA = 'ifany')
cat("\nrs3087243 distribution for patients with Chronic Fatigue Syndrome:\n")
table(cfs_data$rs3087243, useNA = 'ifany')
cat("\nrs3807306 distribution for patients with Chronic Fatigue Syndrome:\n")
table(cfs_data$rs3807306, useNA = 'ifany')
cat("\nrs1800629 distribution for patients with Chronic Fatigue Syndrome:\n")
table(cfs_data$rs1800629, useNA = 'ifany')
cat("\nrs1799724 distribution for patients with Chronic Fatigue Syndrome:\n")
table(cfs_data$rs1799724, useNA = 'ifany')

# Testing Hardy-Weinberg Equilibrium for patients with Chronic Fatigue Syndrome with an infection trigger at their disease onset using exact tst

identificators <- c('rs2476601', 'rs3087243', 'rs3807306', 'rs1800629', 'rs1799724')
levels <- list('rs2476601' = c('AA', 'AG', 'GG'),
               'rs3087243' = c('AA', 'AG', 'GG'),
               'rs3807306' = c('GG', 'GT', 'TT'),
               'rs1800629' = c('AA', 'AG', 'GG'),
               'rs1799724' = c('CC', 'CT', 'TT'))
for(rs in identificators) {
  t <- factor(cfs_data[[rs]], levels = levels[[rs]])
  t <- table(t)
  test <- HWExact(t, verbose = F)
  pval <- test$pval
  cat("\nExact test on HC's ", rs, " biomarker:\n")
  cat("p-value: ", pval, '\n')
  cat("\n--------------------------------------------\n")
}

# Testing Hardy-Weinberg Equilibrium for patients with Chronic Fatigue Syndrome with an infection trigger at their disease onset

cat("\nchi-square test on patients with CFS's rs2476601 biomarker:\n")

AA <- table(cfs_data$rs2476601, useNA = 'ifany')[['AA']]
AG <- table(cfs_data$rs2476601, useNA = 'ifany')[['AG']]
GG <- table(cfs_data$rs2476601, useNA = 'ifany')[['GG']]

pi <- (2 * AA + AG) / (2 * (AA + AG + GG))

cat("p-value: ", chisq.test(table(factor(cfs_data$rs2476601, levels = c('AA', 'AG', 'GG'))),
                            p = c(pi**2, 2 * pi * (1 - pi), (1 - pi)**2))$p.value, "\n")

cat("\n--------------------------------------------\n")
cat("\nchi-square test on patients with CFS's rs3087243 biomarker:\n")

AA <- table(cfs_data$rs3087243, useNA = 'ifany')[['AA']]
AG <- table(cfs_data$rs3087243, useNA = 'ifany')[['AG']]
GG <- table(cfs_data$rs3087243, useNA = 'ifany')[['GG']]

pi <- (2 * AA + AG) / (2 * (AA + AG + GG))

cat("p-value: ", chisq.test(table(factor(cfs_data$rs3087243, levels = c('AA', 'AG', 'GG'))),
                            p = c(pi**2, 2 * pi * (1 - pi), (1 - pi)**2))$p.value, "\n")

cat("\n--------------------------------------------\n")
cat("\nchi-square test on patients with CFS's rs2476601 biomarker:\n")

GG <- table(cfs_data$rs3807306, useNA = 'ifany')[['GG']]
GT <- table(cfs_data$rs3807306, useNA = 'ifany')[['GT']]
TT <- table(cfs_data$rs3807306, useNA = 'ifany')[['TT']]

pi <- (2 * GG + GT) / (2 * (GG + GT + TT))

cat("p-value: ", chisq.test(table(factor(cfs_data$rs3807306, levels = c('GG', 'GT', 'TT'))),
                            p = c(pi**2, 2 * pi * (1 - pi), (1 - pi)**2))$p.value, "\n")

cat("\n--------------------------------------------\n")
cat("\nchi-square test on patients with CFS's rs2476601 biomarker:\n")

AA <- table(cfs_data$rs1800629, useNA = 'ifany')[['AA']]
AG <- table(cfs_data$rs1800629, useNA = 'ifany')[['AG']]
GG <- table(cfs_data$rs1800629, useNA = 'ifany')[['GG']]

pi <- (2 * AA + AG) / (2 * (AA + AG + GG))

cat("p-value: ", chisq.test(table(factor(cfs_data$rs1800629, levels = c('AA', 'AG', 'GG'))),
                            p = c(pi**2, 2 * pi * (1 - pi), (1 - pi)**2))$p.value, "\n")

cat("\n--------------------------------------------\n")
cat("\nchi-square test on patients with CFS's rs1799724 biomarker:\n")

CC <- table(cfs_data$rs1799724, useNA = 'ifany')[['CC']]
CT <- table(cfs_data$rs1799724, useNA = 'ifany')[['CT']]
TT <- table(cfs_data$rs1799724, useNA = 'ifany')[['TT']]

pi <- (2 * CC + CT) / (2 * (CC + CT + TT))

cat("p-value: ", chisq.test(table(factor(cfs_data$rs1799724, levels = c('CC', 'CT', 'TT'))),
                            p = c(pi**2, 2 * pi * (1 - pi), (1 - pi)**2))$p.value, "\n")


##################################################################

### 5.  Comparison of  the genotype distribution of CFS and HC for each genetic marker

#rs2476601
AA_hc <- 0
AG_hc <- table(hc_data$rs2476601, useNA = "ifany")[["AG"]]
GG_hc <- table(hc_data$rs2476601, useNA = "ifany")[["GG"]]
AA_cfs <- table(cfs_data$rs2476601, useNA = "ifany")[["AA"]]
AG_cfs <- table(cfs_data$rs2476601, useNA = "ifany")[["AG"]]
GG_cfs <- table(cfs_data$rs2476601, useNA = "ifany")[["GG"]]
genotype_rs2476601 <- rbind(
  HC  = c(AA = AA_hc,  AG = AG_hc,  GG = GG_hc),
  CFS = c(AA = AA_cfs, AG = AG_cfs, GG = GG_cfs)
)

genotype_rs2476601

chi_rs2476601 <- chisq.test(genotype_rs2476601, correct = FALSE)
chi_rs2476601

fisher_rs2476601 <- fisher.test(genotype_rs2476601)
fisher_rs2476601

#rs3087243
AA_hc <- table(hc_data$rs3087243, useNA = "ifany")[["AA"]]
AG_hc <- table(hc_data$rs3087243, useNA = "ifany")[["AG"]]
GG_hc <- table(hc_data$rs3087243, useNA = "ifany")[["GG"]]

AA_cfs <- table(cfs_data$rs3087243, useNA = "ifany")[["AA"]]
AG_cfs <- table(cfs_data$rs3087243, useNA = "ifany")[["AG"]]
GG_cfs <- table(cfs_data$rs3087243, useNA = "ifany")[["GG"]]

genotype_rs3087243 <- rbind(
  HC  = c(AA = AA_hc,  AG = AG_hc,  GG = GG_hc),
  CFS = c(AA = AA_cfs, AG = AG_cfs, GG = GG_cfs)
)

genotype_rs3087243

chi_rs3087243<-chisq.test(genotype_rs3087243, correct = FALSE)
fisher_rs3087243 <- fisher.test(genotype_rs3087243)

#rs3807306
GG_hc <- table(hc_data$rs3807306, useNA = "ifany")[["GG"]]
GT_hc <- table(hc_data$rs3807306, useNA = "ifany")[["GT"]]
TT_hc <- table(hc_data$rs3807306, useNA = "ifany")[["TT"]]

GG_cfs <- table(cfs_data$rs3807306, useNA = "ifany")[["GG"]]
GT_cfs <- table(cfs_data$rs3807306, useNA = "ifany")[["GT"]]
TT_cfs <- table(cfs_data$rs3807306, useNA = "ifany")[["TT"]]

genotype_rs3807306 <- rbind(
  HC  = c(GG = GG_hc,  GT = GT_hc,  TT = TT_hc),
  CFS = c(GG = GG_cfs, GT = GT_cfs, TT = TT_cfs)
)

genotype_rs3807306

chi_rs3807306 <- chisq.test(genotype_rs3807306, correct = FALSE)
fisher_rs3807306 <- fisher.test(genotype_rs3807306)

#rs1800629
AA_hc <- table(hc_data$rs1800629, useNA = "ifany")[["AA"]]
AG_hc <- table(hc_data$rs1800629, useNA = "ifany")[["AG"]]
GG_hc <- table(hc_data$rs1800629, useNA = "ifany")[["GG"]]

AA_cfs <- table(cfs_data$rs1800629, useNA = "ifany")[["AA"]]
AG_cfs <- table(cfs_data$rs1800629, useNA = "ifany")[["AG"]]
GG_cfs <- table(cfs_data$rs1800629, useNA = "ifany")[["GG"]]

genotype_rs1800629 <- rbind(
  HC  = c(AA = AA_hc,  AG = AG_hc,  GG = GG_hc),
  CFS = c(AA = AA_cfs, AG = AG_cfs, GG = GG_cfs)
)

genotype_rs1800629

chi_rs1800629<- chisq.test(genotype_rs1800629, correct = FALSE)
fisher_rs1800629 <- fisher.test(genotype_rs1800629)

#rs1799724
CC_hc <- table(hc_data$rs1799724, useNA = "ifany")[["CC"]]
CT_hc <- table(hc_data$rs1799724, useNA = "ifany")[["CT"]]
TT_hc <- table(hc_data$rs1799724, useNA = "ifany")[["TT"]]

CC_cfs <- table(cfs_data$rs1799724, useNA = "ifany")[["CC"]]
CT_cfs <- table(cfs_data$rs1799724, useNA = "ifany")[["CT"]]
TT_cfs <- table(cfs_data$rs1799724, useNA = "ifany")[["TT"]]

genotype_rs1799724 <- rbind(
  HC  = c(CC = CC_hc,  CT = CT_hc,  TT = TT_hc),
  CFS = c(CC = CC_cfs, CT = CT_cfs, TT = TT_cfs)
)

genotype_rs1799724

chi_rs1799724 <- chisq.test(genotype_rs1799724, correct = FALSE)
fisher_rs1799724 <- fisher.test(genotype_rs1799724)

# Benjamini-Hochberg correction
p_values_chi <- c(
  rs2476601 = chi_rs2476601$p.value,
  rs3087243 = chi_rs3087243$p.value,
  rs3807306 = chi_rs3807306$p.value,
  rs1800629 = chi_rs1800629$p.value,
  rs1799724 = chi_rs1799724$p.value
)

p_adjusted_chi <- p.adjust(p_values_chi, method = "BH")

print(round(p_adjusted_chi, 5))
significant_chi <- names(p_adjusted_chi[p_adjusted_chi < 0.05])
print(significant_chi)

# Benjamini-Hochberg correction

p_values_list <- c(
  rs2476601 = fisher_rs2476601$p.value,
  rs3087243 = fisher_rs3087243$p.value,
  rs3807306 = fisher_rs3807306$p.value,
  rs1800629 = fisher_rs1800629$p.value,
  rs1799724 = fisher_rs1799724$p.value
)

p_adjusted <- p.adjust(p_values_list, method = "BH")
print(round(p_adjusted, 5))

significant_fisher <- names(p_adjusted[p_adjusted < 0.05])
print(significant_fisher)


##################################################################

### 6. Additive models

#mapping variables to make glm model

# HC = 0, CFS = 1
data6 <- data
data6$Group <- ifelse(data$Group == "CFS", 1, 0)


# from point 2 which allele is rarer
map_rs2476601 <- c(
  'GG' = 0,
  'GA' = 1,
  'AG' = 1,
  'AA' = 2
)
map_rs3087243 <-c(
  'GG' = 0,
  'GA' = 1,
  'AG' = 1,
  'AA' = 2
)

map_rs3807306 <- c(
  'TT' = 2,
  'TG' = 1,
  'GT' = 1,
  'GG' = 0
)

map_rs1800629 <-c(
  'GG' = 0,
  'GA' = 1,
  'AG' = 1,
  'AA' = 2
)

map_rs1799724 <-c(
  'CC' = 0,
  'CT' = 1,
  'TC' = 1,
  'TT' = 2
)

data6$rs2476601 <- map_rs2476601[data$rs2476601]
data6$rs3087243 <- map_rs3087243[data$rs3087243]
data6$rs3807306 <- map_rs3807306[data$rs3807306]
data6$rs1800629 <- map_rs1800629[data$rs1800629]
data6$rs1799724 <- map_rs1799724[data$rs1799724]


model_rs2476601 <- glm(Group ~ rs2476601 + Gender, data = data6, family = binomial("logit"))
summary(model_rs2476601)

model_rs3087243 <- glm(Group ~ rs3087243 + Gender, data = data6, family = binomial("logit"))
summary(model_rs3087243)

model_rs3807306 <- glm(Group ~ rs3807306 + Gender, data = data6, family = binomial("logit"))
summary(model_rs3807306)

model_rs1800629 <- glm(Group ~ rs1800629 + Gender, data = data6, family = binomial("logit"))
summary(model_rs1800629)

model_rs1799724 <- glm(Group ~ rs1799724 + Gender, data = data6, family = binomial("logit"))
summary(model_rs1799724)

# Benjamini-Hochberg correction
p_values <- c(
  summary(model_rs2476601)$coefficients[2, 4],
  summary(model_rs3087243)$coefficients[2, 4],
  summary(model_rs3807306)$coefficients[2, 4],
  summary(model_rs1800629)$coefficients[2, 4],
  summary(model_rs1799724)$coefficients[2, 4]
)

names(p_values) <- c("rs2476601", "rs3087243", "rs3807306", "rs1800629", "rs1799724")

p_adjusted_BH <- p.adjust(p_values, method = "BH")

print("Adjusted P-values (Benjamini-Hochberg):")
print(round(p_adjusted_BH, 5))

significant_markers <- names(p_adjusted_BH[p_adjusted_BH < 0.05])
print(paste("Significant markers after correction:", paste(significant_markers, collapse = ", ")))

final_model <- glm(Group ~ rs3087243 + rs2476601 + Gender,
                   data = data6,
                   family = binomial("logit"))
summary(final_model)
#Both genetic markers remained statistically significant in the combined model
#Gender was not a significant predictor in this model

#checking the full model
full_model <- glm(Group ~ rs3087243 + rs2476601 + rs3807306 +rs1800629+ rs1799724+ Gender,
                  data = data6,
                  family = "binomial")
summary(full_model)

#confusion matrix and accuracy of final model
actual_values <- final_model$y
probabilities <- fitted(final_model)
predicted_classes <- ifelse(probabilities > 0.5, 1, 0)


conf_matrix <- table(Actual = actual_values, Predicted = predicted_classes)
accuracy <- sum(diag(conf_matrix)) / sum(conf_matrix)

print("Confusion Matrix")
print(conf_matrix)
print(paste("Final Model Accuracy:", round(accuracy * 100, 2), "%"))

# genetic variants are statistically significant risk factors, they have limited predictive power as a standalone diagnostic too



##################################################################

### 7. Can these genetic data be used as a diagnostic tool for CFS?

dim(na.omit(data6[, c('rs2476601', 'rs3087243', 'Gender')]))

# logistic regression with logit linkage function
model_logit <- glm(Group ~ rs2476601 + rs3087243 + Gender,
                   data = data6, family = binomial("logit"), na.action = na.exclude)

#backward_logit <- step(model_logit, direction='backward')

data6$Prob_logit <- predict(model_logit, type = "response")
roc_logit <- roc(Group ~ Prob_logit, data = data6)

# AUC
print(roc_logit$auc)

# ROC Curve
plot(roc_logit, col = "#8ed087", main = "ROC Curve for CFS diagnosis (logit)")

# Optimal cutpoints for logistic regression with logit linkage function

ROC01_logit <- optimal.cutpoints(Prob_logit ~ Group, tag.healthy = 0, data = data6, methods = "ROC01")
SpEqualSe_logit <- optimal.cutpoints(Prob_logit ~ Group, tag.healthy = 0, data = data6, methods = "SpEqualSe")
MaxEfficiency_logit <- optimal.cutpoints(Prob_logit ~ Group, tag.healthy = 0, data = data6, methods = "MaxEfficiency")
summary(ROC01_logit)
summary(SpEqualSe_logit)
summary(MaxEfficiency_logit)

# logistic regression with probit linkage function
model_probit <- glm(Group ~ rs2476601 + rs3087243 + Gender,
                    data = data6, family = binomial("probit"), na.action = na.exclude)

data6$Prob_probit <- predict(model_probit, type = "response")
roc_probit <- roc(Group ~ Prob_probit, data = data6)

# AUC
print(roc_probit$auc)

# ROC Curve
plot(roc_probit, col = "#f4d059", main = "ROC Curve for CFS diagnosis (probit)")

# Optimal cutpoints for logistic regression with probit linkage function

ROC01_probit <- optimal.cutpoints(Prob_probit ~ Group, tag.healthy = 0, data = data6, methods = "ROC01")
SpEqualSe_probit <- optimal.cutpoints(Prob_probit ~ Group, tag.healthy = 0, data = data6, methods = "SpEqualSe")
MaxEfficiency_probit <- optimal.cutpoints(Prob_probit ~ Group, tag.healthy = 0, data = data6, methods = "MaxEfficiency")
summary(ROC01_probit)
summary(SpEqualSe_probit)
summary(MaxEfficiency_probit)