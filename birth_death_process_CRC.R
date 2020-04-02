

# Libraries ---------------------------------------------------------------

library(ggplot2)
library(stringr)
library(plyr)
library(dplyr)
library(reshape2)
library(tidyr)
library(stringi)
library(data.table)
library(xlsx)



# Model 2: Deterministic N1, Probabilistic N2 -------------------------------------------------------------------

#Initial cell
N0 <- 1
N2cases <- 1000

#Growth rates
wAPC <- 0.05
wPOLE <- 0
wKRASsm <- 0.05
wKRASwm <- 0.025

#Time of population births
tMAX <- 14
step_size <- .001
t <- seq(0, tMAX, by = step_size) 
tPOLEm <- t[length(t)/10]
tHist <- 5

# #Mutation rates. The probabiility will be in the sample function
muKRASsm <- 10^-8
muKRASwm <- 10^-8
muKRASsm_POLEm <- 10^-7
muKRASwm_POLEm <- 10^-6

ii <- 1
while (ii<=length(t)) {
  
  if(t[ii] == 0) {
    nAPC <- vector(mode = "numeric", length = length(t))
    nAPC[1] <- N0
    nAPC.POLEmMean <- vector(mode = "numeric", length = length(t))
    variance <- vector(mode = "numeric", length = length(t))
    N2sim <- matrix(0, nrow = N2cases, ncol = length(t))
    M2 <- 1:nrow(N2sim)
  }
  
  if(t[ii] > 0 & t[ii] < tPOLEm) {
    nAPC[ii] <- N0*exp((1+wAPC)*t[ii]) 
  }
  
  if(t[ii] == tPOLEm) {
    nAPC[ii] <- N0*exp((1+wAPC)*t[ii])
    nAPC.POLEmMean[ii] <- 1
    variance[ii] <- N0*exp((1+wAPC)*(t[ii]-tPOLEm))*( exp((1+wAPC)*(t[ii]-tPOLEm)) - 1 )
    last(nAPC.POLEmMean)*(exp(1+wAPC)-1)
    N2sim[,ii] <- 1 
  }
  
  if(t[ii] > tPOLEm) {
    nAPC[ii] <- N0*exp((1+wAPC)*t[ii])
    nAPC.POLEmMean[ii] <- N0*exp((1+wAPC)*(t[ii]-tPOLEm))
    variance[ii] <- N0*exp((1+wAPC)*(t[ii]-tPOLEm))*( exp((1+wAPC)*(t[ii]-tPOLEm)) - 1 )
    
    x <- runif(N2cases,0,1)
    P2 <- (1+wAPC)*step_size*N2sim[,ii-1]
    
    positives <- which(x <= P2)
    negatives <- which(x > P2)
    
    N2sim[positives,ii] <- N2sim[positives, ii-1] + 1
    N2sim[negatives,ii] <- N2sim[negatives, ii-1] 
  }
  ii <- ii+1
}


# Plots -------------------------------------------------------------------


# 1. Mean and deviation ------------------------------------------------------

N2means <- colMeans(N2sim)

z <- data.frame(t, nAPC.POLEmMean,
                  N2means)

z <- reshape2::melt(z, id.vars = "t") %>%
   rename(Time = t, Population = variable, Cell_count = value) %>%
   filter(Time %in% seq(min(t), max(t), 0.1))

ggplot(z,
        aes(x = Time, y = Cell_count)) +
   geom_line(aes(colour = Population), size = 1.5) +
   labs(title = paste(N2cases, "Simulations from Model 2: Deterministic N1, Mean N2")) +
   xlab("Time") + ylab("Number of cells") +
   theme(
     plot.title = element_text(face = "bold"),
     axis.title.x = element_text(color="black", size=12, face="bold"),
     axis.title.y = element_text(color="black", size=12, face="bold"),
     legend.title = element_text(color="black", size=12, face="bold")
   ) +
   scale_colour_manual(values = c("red", "blue"))


# 2. Histogram --------------------------------------------------------------

z <- data.frame(t, t(N2sim)) %>%
  filter(t == tHist)

z <- reshape2::melt(z, id.vars = "t") %>%
  rename(Time = t, Population = variable, Cell_count = value)

x_Neg_bino <- seq(min(z$Cell_count), max(z$Cell_count),1)
times <- tHist
aaa <- exp(-times*(1+wAPC))
bbb <- (1 - exp( -times*(1+wAPC) ) )^(x_Neg_bino-1)
y_Neg_bino <- aaa*bbb
df_Neg_bino <- data.frame(x_Neg_bino, y_Neg_bino)
sum(df_Neg_bino$y_Neg_bino)


ggplot(z) +
  geom_histogram(aes(x = Cell_count, y = ..density..),
                 binwidth = 1, fill = "blue", colour = "black") +
  labs(title = paste0("Time = ", tHist)) +
  ylab("Count") + xlab("Number of cells") +
  theme(
    plot.title = element_text(face = "bold"),
    axis.title.x = element_text(color="black", size=12, face="bold"),
    axis.title.y = element_text(color="black", size=12, face="bold"),
    legend.title = element_text(color="black", size=12, face="bold")
  ) +
  geom_line(data = df_Neg_bino,
            aes(x = x_Neg_bino, y = y_Neg_bino), colour = "red", size = 1.5)



# Saving plots and data ---------------------------------------------------

#ggsave(paste0(N2cases, "Runs", "Hist_time", tHist, ".pdf"))


