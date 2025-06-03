
#d = read.table("/Users/subhashmahamkali/Downloads/2.txt")
d = read.table("/Users/subhashmahamkali/Downloads/0.comb_selec_gwas_genes.txt")
library(readxl)
b <- read_excel("/Users/subhashmahamkali/Downloads/1.miscellaneous/sorghum project/7.sorghum.annotations.xlsx", skip = 2)

b <- b[, c(1, 9)]
b <- b[-c(1, 2), ]
b <- b[!duplicated(b$Gene), ]

colnames(b) <- c("Gene", "Description")
colnames(d)[22] <- "Gene"
merged_data <- merge(d, b, by = "Gene", all.x = TRUE)
write.table(merged_data, "/Users/subhashmahamkali/Downloads/1.miscellaneous/sorghum project/166_genes.txt", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

library(VennDiagram)
Bal = 271
pos = 734
selec = 76
GWAS = 1023

Gwas to selec = 225
genes =166


# Define values
Bal <- 271     # Unique to Bal
Pos <- 734     # Unique to Pos
Overlap <- 76  # Shared between Bal and Pos
selection <- 76
GWAS <- 1023
overla <- 25

# Create the Venn diagram
selec <- draw.pairwise.venn(
  area1 = Bal + Overlap,   # Total in Bal (including overlap)
  area2 = Pos + Overlap,   # Total in Pos (including overlap)
  cross.area = Overlap,    # Intersection (overlap)
  category = c("BS", "PS"), # Labels
  fill = c("#00FFFF","#8A2BE2"), # Colors
  alpha = 0.5,             # Transparency
  lty = "solid",           # Line type
  cex = 1,               # Font size
  cat.cex = 1,           # Category font size
  cat.col = c("#00FFFF","#8A2BE2") # Category text color
)
ai_colors2 <- c("#FF1493", "#00FF7F", "#1E90FF")  # Deep Pink, Spring Green, and Dodger Blue

# Save or display the plot
grid.draw(selec)
ai_colors <- c("#FF1493","#32CD32")  # Cyan, Blue-Violet, and Neon Green
s <- draw.pairwise.venn(
  area1 = selection + overla,   # Total in Bal (including overlap)
  area2 = GWAS + overla,   # Total in Pos (including overlap)
  cross.area = overla,    # Intersection (overlap)
  category = c("GWAS", "selection"), # Labels
  fill = c("#FF1493","#32CD32"), # Colors
  alpha = 0.5,             # Transparency
  lty = "solid",           # Line type
  cex = 1,               # Font size
  cat.cex = 1, ext.text = FALSE,           # Category font size
  cat.col = c("#FF1493","#32CD32") # Category text color
)
grid.draw(s)




library(data.table)
dp1=fread("/Users/subhashmahamkali/Downloads/1.miscellaneous/sorghum project/positive_selection/dp1_output.csv", header=T,data.table=F)
dp1[,1] = as.numeric(dp1[,1])  
col = ifelse(dp1[,1] %% 2 == 1, "#00000066", "#BEBEBE99")
dp1 = na.omit(dp1)  # Removes rows with any NA
thr = quantile(dp1[,8], 0.99)  # Compute quantile on cleaned data

png("SM_Fst.png", height = 8, width = 18, res = 600, units = "in")
plot(dp1[,2], dp1[,8], col=col, pch=16, cex=0.4, 
     bty="l", xlim=c(0,691),axes=F,cex.lab=0.6, 
     xlab="", ylab="", font.lab=2)
axis(2,las=2,tck=-.03,cex.axis=1.5,font.axis=2)
mtext("Fst", side=2, line=3.5, font=2, cex=1)
segments(x0=0, x1=720, y0=thr, y1=thr, col="red", lty=2, lwd=2)
axis(1,at=(ch[,3]-ch[,2]/2)/1e6,labels=1:10,tck=-0.03,cex.axis=1.5,font.axis=2)
mtext("Chromosome", side=1, line=2, font=2, cex=1)
dev.off()
