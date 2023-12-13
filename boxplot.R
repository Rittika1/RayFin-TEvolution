# Load the necessary library
library(ggplot2)
library(ggplot)
library(viridis)
library(grid)
library(gridExtra)

setwd("/home/rittika/Documents/PolypterusBichir/TEproject/")
# Read data from a file (replace 'your_data_file.csv' with your actual file path)
data <- read.csv('TEdataWithLifeHistory_updatedPCA.csv')

# Create a box plot using ggplot2
ggplot(data, aes(x = teleostOrNot, y = GenomeSize, fill = teleostOrNot)) +
  geom_boxplot() +
  labs(x = "Teleost Status", y = "Genome Size") +
 scale_fill_manual(values = c("teleost" = "lightblue", "NonTeleost" = "violet")) +
ggtitle("Box Plot of Genome Size for Teleost and Non-Teleost Species")

par(mfrow = c(2,2))
plots <- list()
for (i in c("DNA","LTR", "LINE", "SINE")) {
 bplot <- ggplot(data, aes(x = teleostOrNot, y = i, fill = teleostOrNot)) +
    geom_boxplot() +
    labs(x = "Teleost Status", y = paste(i)) +
    scale_fill_manual(values = c("teleost" = "steelblue", "NonTeleost" = "lightgreen")) +
    ggtitle(paste("Box Plot of", i, "for Teleost and Non-Teleost Species"))
 plots[[i]] <- bplot
} 
multiplot(plotlist = plots, cols = 2)

#Boxplot for LTR
ggplot(data, aes(x = teleostOrNot, y = LTR, fill = teleostOrNot)) +
  geom_boxplot() +
  labs(x = "Teleost Status", y = "LTR") +
  scale_fill_manual(values = c("teleost" = "steelblue", "NonTeleost" = "lightgreen")) +
  ggtitle("Box Plot of LTR for Teleost and Non-Teleost Species")

#Boxplot for LINE
ggplot(data, aes(x = teleostOrNot, y = LINE, fill = teleostOrNot)) +
  geom_boxplot() +
  labs(x = "Teleost Status", y = "LINE") +
  scale_fill_manual(values = c("teleost" = "steelblue", "NonTeleost" = "lightgreen")) +
  ggtitle("Box Plot of LINE for Teleost and Non-Teleost Species")

#Boxplot for SINE
ggplot(data, aes(x = teleostOrNot, y = SINE, fill = teleostOrNot)) +
  geom_boxplot() +
  labs(x = "Teleost Status", y = "SINE") +
  scale_fill_manual(values = c("teleost" = "steelblue", "NonTeleost" = "lightgreen")) +
  ggtitle("Box Plot of SINE for Teleost and Non-Teleost Species")

#Boxplot for DNA
ggplot(data, aes(x = teleostOrNot, y = DNA, fill = teleostOrNot)) +
  geom_boxplot() +
  labs(x = "Teleost Status", y = "DNA") +
  scale_fill_manual(values = c("teleost" = "steelblue", "NonTeleost" = "lightgreen")) +
  ggtitle("Box Plot of DNA for Teleost and Non-Teleost Species")

