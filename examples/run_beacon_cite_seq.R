#Load required packages
library(Seurat)
library(ggplot2)
library(fitdistrplus)
library(uwot)

#Set seed
set.seed(1234)

#Read in Kotliarov data
seurat <- readRDS(".../RDS/H1_day0_demultilexed_singlets.RDS")

########################
#### Fit Background ####
########################

# Keep the Cite-Seq transcriptomes
cells_keep <- rownames(seurat@meta.data)
rna <- seurat@raw.data
rna <- rna[,cells_keep]
adt <- seurat@assay$CITE@raw.data
data <- as.matrix(adt)

#### Learn a "universal" background, with truncated negative binomial ####

# First visualize the entire count distribution
x <- as.vector(log10(data + 1))
h <- hist(x, breaks = 1000, plot = FALSE)
df <- data.frame(
  x = h$mids,
  prob = h$counts/sum(h$counts)
)
df_plot <- subset(df, prob > 0)
ggplot(df_plot, aes(x, prob)) +
  geom_point(size = 1) +
  scale_y_log10() +
  labs(
    x = expression(log[10](ADT~counts+1)),
    y = "Probability"
  ) +
  theme_classic()

#Learn a universal background
Tmax <- 10
data_back <- data[data<=Tmax]

#Define the PDF
dtrunc_nbinom <- function(x, size, mu, log = FALSE) {
  dens <- dnbinom(x, size = size, mu = mu) /
    pnbinom(Tmax, size = size, mu = mu)
  dens[x > Tmax] <- 0
  if (log)
    return(log(dens))
  dens
}

#Define the CDF
ptrunc_nbinom <- function(q, size, mu, lower.tail = TRUE, log.p = FALSE) {
  p <- pnbinom(pmin(q, Tmax), size = size, mu = mu) /
    pnbinom(Tmax, size = size, mu = mu)
  p[q < 0] <- 0
  p[q >= Tmax] <- 1
  if (!lower.tail)
    p <- 1 - p
  if (log.p)
    p <- log(p)
  p
}

#Fit truncated NB by maximum likelihood
fit_back_nb <- fitdist(
  data_back,
  distr = "trunc_nbinom",
  start = list(
    size = 1,
    mu = mean(data_back)
  ),
  lower = c(1e-6, 1e-6)
)

#### Assess fit ####

#Parameters
nb_params <- list()
nb_params$size <- fit_back_nb$estimate["size"]
nb_params$mu   <- fit_back_nb$estimate["mu"]

#Maximum count used for fitting
obs <- table(factor(data_back, levels = 0:Tmax))

#Expected probabilities under the fitted truncated NB
p <- dnbinom(0:Tmax, size = nb_params$size, mu = nb_params$mu)
p <- p / sum(p)

#Expected counts
exp <- length(data_back) * p
plot_df <- data.frame(
  Count = 0:Tmax,
  Observed = as.numeric(obs),
  Expected = exp
)

ggplot(plot_df, aes(x = Count)) +
  geom_col(
    aes(y = Observed),
    fill = "grey80",
    color = "black",
    width = 0.9
  ) +
  geom_point(
    aes(y = Expected, color = "Fit"),
    size = 3
  ) +
  geom_line(
    aes(y = Expected, color = "Fit", group = 1),
    linewidth = 1
  ) +
  scale_color_manual(
    name = NULL,
    values = c(Fit = "steelblue")
  ) +
  labs(
    x = "ADT count",
    y = "Frequency"
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = c(0.98, 0.98),
    legend.justification = c(1, 1),
    legend.background = element_blank()
  )


################################
#### Learn the BEACON model ####
################################
source(".../R/beacon_cite_seq_suite.R")
data_unfilt <- data
m <- nrow(data_unfilt)
n <- ncol(data_unfilt)

# BEACON parameters
sig_max <- m # Maximum allowable signal barcodes
min_count <- 2 # Minimum allowable counts for a signal barcode
max_iter <- 50 # Maximum number of iterations
nb_params <- nb_params

#Initialize background distribution
thetas_b0 <- rowSums(data_unfilt) / sum(rowSums(data_unfilt))

#Learn BEACON parameters
output <- learn_beacon_params_cite(data_unfilt,thetas_b0,sig_max,min_count,max_iter,nb_params)

#Remove the background
thetas_b <- output$thetas_b
alphas <- output$alphas
thetas_b_mat <- matrix(thetas_b,m,n,byrow=FALSE)
alphas_mat <- matrix(alphas,m,n,byrow=TRUE)
counts <- colSums(data_unfilt)
counts_mat <- matrix(counts,m,n,byrow=TRUE)
thetas_s <- output$thetas_s
data_beacon <- ceiling(counts_mat * (1-alphas_mat) * thetas_s)
rownames(data_beacon) <- rownames(data_unfilt)
colnames(data_beacon) <- colnames(data_unfilt)


