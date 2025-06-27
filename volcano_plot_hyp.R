
#group 1!!!
dev.off()
# xlim=limit of x axis, you can change that based on your data
with(resultsproteomics_hyp, plot(diff_CN_A_Yes_Hyp_No_vs_CN_A_No, -log10(p_CN_A_Yes_Hyp_No_vs_CN_A_No), pch=20, col="grey41", main="Volcano plot", xlim=c(-1.1,1.1), ylim=c(0, 5)))


# Add colored points: red if padj<0.05, orange of log2FC>1, green if both)
with(subset(resultsproteomics_hyp, p_CN_A_Yes_Hyp_No_vs_CN_A_No<.05), points(diff_CN_A_Yes_Hyp_No_vs_CN_A_No, -log10(p_CN_A_Yes_Hyp_No_vs_CN_A_No), pch=20, col="tomato3"))


# Label points with the textxy function from the calibrate plot, you need to look at the p value for the top20 and modify it here in the script to only label the top 20 proteins
library(calibrate)

with(subset(resultsproteomics_hyp, p_CN_A_Yes_Hyp_No_vs_CN_A_No<.00475), points(diff_CN_A_Yes_Hyp_No_vs_CN_A_No, -log10(p_CN_A_Yes_Hyp_No_vs_CN_A_No), pch=20, col="dodgerblue3"))
#Here you need to include the p-value of the 21th protein so that it only names the top 20 proteins (with highest p-value than the 21th)
#the one with the lowest p-value was <0.00001 but was changed to 0.00001 in order to be able to include it in the graph. NOT WORKING
with(subset(resultsproteomics_hyp, p_CN_A_Yes_Hyp_No_vs_CN_A_No<.00475), textxy(diff_CN_A_Yes_Hyp_No_vs_CN_A_No, -log10(p_CN_A_Yes_Hyp_No_vs_CN_A_No), labs=Protein, cex=.4))

#group 2!!!
dev.off()
# xlim=limit of x axis, you can change that based on your data
with(resultsproteomics_hyp, plot(diff_CN_A_Yes_Hyp_Yes_vs_CN_A_No, -log10(p_CN_A_Yes_Hyp_Yes_vs_CN_A_No), pch=20, col="grey41", main="Volcano plot", xlim=c(-1.7,1.3), ylim=c(0, 5)))

# Add colored points: red if padj<0.05, orange of log2FC>1, green if both)
with(subset(resultsproteomics_hyp, p_CN_A_Yes_Hyp_Yes_vs_CN_A_No<.05), points(diff_CN_A_Yes_Hyp_Yes_vs_CN_A_No, -log10(p_CN_A_Yes_Hyp_Yes_vs_CN_A_No), pch=20, col="tomato3"))


# Label points with the textxy function from the calibrate plot, you need to look at the p value for the top20 and modify it here in the script to only label the top 20 proteins
library(calibrate)
with(subset(resultsproteomics_hyp, p_CN_A_Yes_Hyp_Yes_vs_CN_A_No<.01429), points(diff_CN_A_Yes_Hyp_Yes_vs_CN_A_No, -log10(p_CN_A_Yes_Hyp_Yes_vs_CN_A_No), pch=20, col="dodgerblue3"))
#Here you need to include the p-value of the 21th protein so that it only names the top 20 proteins (with highest p-value than the 21th)

with(subset(resultsproteomics_hyp, p_CN_A_Yes_Hyp_Yes_vs_CN_A_No<.01429), textxy(diff_CN_A_Yes_Hyp_Yes_vs_CN_A_No, -log10(p_CN_A_Yes_Hyp_Yes_vs_CN_A_No), labs=Protein, cex=.4))

#group 3!!! adjust p value
dev.off()
# xlim=limit of x axis, you can change that based on your data
with(resultsproteomics_hyp, plot(diff_MCI_A_Yes_Hyp_No_vs_CN_A_No, -log10(p_MCI_A_Yes_Hyp_No_vs_CN_A_No), pch=20, col="grey41", main="Volcano plot", xlim=c(-0.8,-0.6), ylim=c(4.2, 4.6)))

# Add colored points: red if padj<0.05, orange of log2FC>1, green if both)
with(subset(resultsproteomics_hyp, p_MCI_A_Yes_Hyp_No_vs_CN_A_No<.05), points(diff_MCI_A_Yes_Hyp_No_vs_CN_A_No, -log10(p_MCI_A_Yes_Hyp_No_vs_CN_A_No), pch=20, col="tomato3"))


# Label points with the textxy function from the calibrate plot, you need to look at the p value for the top20 and modify it here in the script to only label the top 20 proteins
library(calibrate)
with(subset(resultsproteomics_hyp, p_MCI_A_Yes_Hyp_No_vs_CN_A_No<.00005), points(diff_MCI_A_Yes_Hyp_No_vs_CN_A_No, -log10(p_MCI_A_Yes_Hyp_No_vs_CN_A_No), pch=20, col="dodgerblue3"))
#Here you need to include the p-value of the 21th protein so that it only names the top 20 proteins (with highest p-value than the 21th)

with(subset(resultsproteomics_hyp, p_MCI_A_Yes_Hyp_No_vs_CN_A_No<.00005), textxy(diff_MCI_A_Yes_Hyp_No_vs_CN_A_No, -log10(p_MCI_A_Yes_Hyp_No_vs_CN_A_No), labs=Protein, cex=.4))


#group 4!!!
dev.off()
# xlim=limit of x axis, you can change that based on your data
with(resultsproteomics_hyp, plot(diff_MCI_A_Yes_Hyp_Yes_vs_CN_A_No, -log10(p_MCI_A_Yes_Hyp_Yes_vs_CN_A_No), pch=20, col="grey41", main="Volcano plot", xlim=c(-1.4,1.4), ylim=c(0, 5)))

# Add colored points: red if padj<0.05, orange of log2FC>1, green if both)
with(subset(resultsproteomics_hyp, p_MCI_A_Yes_Hyp_Yes_vs_CN_A_No<.05), points(diff_MCI_A_Yes_Hyp_Yes_vs_CN_A_No, -log10(p_MCI_A_Yes_Hyp_Yes_vs_CN_A_No), pch=20, col="tomato3"))


# Label points with the textxy function from the calibrate plot, you need to look at the p value for the top20 and modify it here in the script to only label the top 20 proteins
library(calibrate)
with(subset(resultsproteomics_hyp, p_MCI_A_Yes_Hyp_Yes_vs_CN_A_No<.00096), points(diff_MCI_A_Yes_Hyp_Yes_vs_CN_A_No, -log10(p_MCI_A_Yes_Hyp_Yes_vs_CN_A_No), pch=20, col="dodgerblue3"))
#Here you need to include the p-value of the 21th protein so that it only names the top 20 proteins (with highest p-value than the 21th)

with(subset(resultsproteomics_hyp, p_MCI_A_Yes_Hyp_Yes_vs_CN_A_No<.00096), textxy(diff_MCI_A_Yes_Hyp_Yes_vs_CN_A_No, -log10(p_MCI_A_Yes_Hyp_Yes_vs_CN_A_No), labs=Protein, cex=.4))

#group 5!!!
dev.off()
# xlim=limit of x axis, you can change that based on your data
with(resultsproteomics_hyp, plot(diff_AD_A_Yes_Hyp_No_vs_CN_A_No, -log10(p_AD_A_Yes_Hyp_No_vs_CN_A_No), pch=20, col="grey41", main="Volcano plot", xlim=c(-3.05,3.05), ylim=c(0, 5)))

# Add colored points: red if padj<0.05, orange of log2FC>1, green if both)
with(subset(resultsproteomics_hyp, p_AD_A_Yes_Hyp_No_vs_CN_A_No<.05), points(diff_AD_A_Yes_Hyp_No_vs_CN_A_No, -log10(p_AD_A_Yes_Hyp_No_vs_CN_A_No), pch=20, col="tomato3"))


# Label points with the textxy function from the calibrate plot, you need to look at the p value for the top20 and modify it here in the script to only label the top 20 proteins
library(calibrate)
with(subset(resultsproteomics_hyp, p_AD_A_Yes_Hyp_No_vs_CN_A_No<.01017), points(diff_AD_A_Yes_Hyp_No_vs_CN_A_No, -log10(p_AD_A_Yes_Hyp_No_vs_CN_A_No), pch=20, col="dodgerblue3"))
#Here you need to include the p-value of the 21th protein so that it only names the top 20 proteins (with highest p-value than the 21th)

with(subset(resultsproteomics_hyp, p_AD_A_Yes_Hyp_No_vs_CN_A_No<.01017), textxy(diff_AD_A_Yes_Hyp_No_vs_CN_A_No, -log10(p_AD_A_Yes_Hyp_No_vs_CN_A_No), labs=Protein, cex=.4))

#group 6!!!
dev.off()
# xlim=limit of x axis, you can change that based on your data
with(resultsproteomics_hyp, plot(diff_AD_A_Yes_Hyp_Yes_vs_CN_A_No, -log10(p_AD_A_Yes_Hyp_Yes_vs_CN_A_No), pch=20, col="grey41", main="Volcano plot", xlim=c(-1.8,1.8), ylim=c(0, 4)))

# Add colored points: red if padj<0.05, orange of log2FC>1, green if both)
with(subset(resultsproteomics_hyp, p_AD_A_Yes_Hyp_Yes_vs_CN_A_No<.05), points(diff_AD_A_Yes_Hyp_Yes_vs_CN_A_No, -log10(p_AD_A_Yes_Hyp_Yes_vs_CN_A_No), pch=20, col="tomato3"))


# Label points with the textxy function from the calibrate plot, you need to look at the p value for the top20 and modify it here in the script to only label the top 20 proteins
library(calibrate)
with(subset(resultsproteomics_hyp, p_AD_A_Yes_Hyp_Yes_vs_CN_A_No<.00758), points(diff_AD_A_Yes_Hyp_Yes_vs_CN_A_No, -log10(p_AD_A_Yes_Hyp_Yes_vs_CN_A_No), pch=20, col="dodgerblue3"))
#Here you need to include the p-value of the 21th protein so that it only names the top 20 proteins (with highest p-value than the 21th)

with(subset(resultsproteomics_hyp, p_AD_A_Yes_Hyp_Yes_vs_CN_A_No<.00758), textxy(diff_AD_A_Yes_Hyp_Yes_vs_CN_A_No, -log10(p_AD_A_Yes_Hyp_Yes_vs_CN_A_No), labs=Protein, cex=.4))

#group 7!!!
dev.off()
# xlim=limit of x axis, you can change that based on your data
with(resultsproteomics_hyp, plot(diff_CN_A_Yes_Hyp_Yes_vs_CN_A_Yes_Hyp_No, -log10(p_CN_A_Yes_Hyp_Yes_vs_CN_A_Yes_Hyp_No), pch=20, col="grey41", main="Volcano plot", xlim=c(-1.8,1.8), ylim=c(0, 4)))

# Add colored points: red if padj<0.05, orange of log2FC>1, green if both)
with(subset(resultsproteomics_hyp, p_CN_A_Yes_Hyp_Yes_vs_CN_A_Yes_Hyp_No<.05), points(diff_CN_A_Yes_Hyp_Yes_vs_CN_A_Yes_Hyp_No, -log10(p_CN_A_Yes_Hyp_Yes_vs_CN_A_Yes_Hyp_No), pch=20, col="tomato3"))


# Label points with the textxy function from the calibrate plot, you need to look at the p value for the top20 and modify it here in the script to only label the top 20 proteins
library(calibrate)
with(subset(resultsproteomics_hyp, p_CN_A_Yes_Hyp_Yes_vs_CN_A_Yes_Hyp_No<.02686), points(diff_CN_A_Yes_Hyp_Yes_vs_CN_A_Yes_Hyp_No, -log10(p_CN_A_Yes_Hyp_Yes_vs_CN_A_Yes_Hyp_No), pch=20, col="dodgerblue3"))
#Here you need to include the p-value of the 21th protein so that it only names the top 20 proteins (with highest p-value than the 21th)

with(subset(resultsproteomics_hyp, p_CN_A_Yes_Hyp_Yes_vs_CN_A_Yes_Hyp_No<.02686), textxy(diff_CN_A_Yes_Hyp_Yes_vs_CN_A_Yes_Hyp_No, -log10(p_CN_A_Yes_Hyp_Yes_vs_CN_A_Yes_Hyp_No), labs=Protein, cex=.4))

#group 8!!!
dev.off()
# xlim=limit of x axis, you can change that based on your data
with(resultsproteomics_hyp, plot(diff_MCI_A_Yes_Hyp_Yes_vs_MCI_A_Yes_Hyp_No, -log10(p_MCI_A_Yes_Hyp_Yes_vs_MCI_A_Yes_Hyp_No), pch=20, col="grey41", main="Volcano plot", xlim=c(-1.6,1.6), ylim=c(0, 4.5)))

# Add colored points: red if padj<0.05, orange of log2FC>1, green if both)
with(subset(resultsproteomics_hyp, p_MCI_A_Yes_Hyp_Yes_vs_MCI_A_Yes_Hyp_No<.05), points(diff_MCI_A_Yes_Hyp_Yes_vs_MCI_A_Yes_Hyp_No, -log10(p_MCI_A_Yes_Hyp_Yes_vs_MCI_A_Yes_Hyp_No), pch=20, col="tomato3"))


# Label points with the textxy function from the calibrate plot, you need to look at the p value for the top20 and modify it here in the script to only label the top 20 proteins
library(calibrate)
with(subset(resultsproteomics_hyp, p_MCI_A_Yes_Hyp_Yes_vs_MCI_A_Yes_Hyp_No<.01545), points(diff_MCI_A_Yes_Hyp_Yes_vs_MCI_A_Yes_Hyp_No, -log10(p_MCI_A_Yes_Hyp_Yes_vs_MCI_A_Yes_Hyp_No), pch=20, col="dodgerblue3"))
#Here you need to include the p-value of the 21th protein so that it only names the top 20 proteins (with highest p-value than the 21th)

with(subset(resultsproteomics_hyp, p_MCI_A_Yes_Hyp_Yes_vs_MCI_A_Yes_Hyp_No<.01545), textxy(diff_MCI_A_Yes_Hyp_Yes_vs_MCI_A_Yes_Hyp_No, -log10(p_MCI_A_Yes_Hyp_Yes_vs_MCI_A_Yes_Hyp_No), labs=Protein, cex=.4))

#group 9!!!
dev.off()
# xlim=limit of x axis, you can change that based on your data
with(resultsproteomics_hyp, plot(diff_AD_A_Yes_Hyp_Yes_vs_AD_A_Yes_Hyp_No, -log10(p_AD_A_Yes_Hyp_Yes_vs_AD_A_Yes_Hyp_No), pch=20, col="grey41", main="Volcano plot", xlim=c(-2.9,2.9), ylim=c(0, 4)))

# Add colored points: red if padj<0.05, orange of log2FC>1, green if both)
with(subset(resultsproteomics_hyp, p_AD_A_Yes_Hyp_Yes_vs_AD_A_Yes_Hyp_No<.05), points(diff_AD_A_Yes_Hyp_Yes_vs_AD_A_Yes_Hyp_No, -log10(p_AD_A_Yes_Hyp_Yes_vs_AD_A_Yes_Hyp_No), pch=20, col="tomato3"))


# Label points with the textxy function from the calibrate plot, you need to look at the p value for the top20 and modify it here in the script to only label the top 20 proteins
library(calibrate)
with(subset(resultsproteomics_hyp, p_AD_A_Yes_Hyp_Yes_vs_AD_A_Yes_Hyp_No<.04141), points(diff_AD_A_Yes_Hyp_Yes_vs_AD_A_Yes_Hyp_No, -log10(p_AD_A_Yes_Hyp_Yes_vs_AD_A_Yes_Hyp_No), pch=20, col="dodgerblue3"))
#Here you need to include the p-value of the 21th protein so that it only names the top 20 proteins (with highest p-value than the 21th)

with(subset(resultsproteomics_hyp, p_AD_A_Yes_Hyp_Yes_vs_AD_A_Yes_Hyp_No<.04141), textxy(diff_AD_A_Yes_Hyp_Yes_vs_AD_A_Yes_Hyp_No, -log10(p_AD_A_Yes_Hyp_Yes_vs_AD_A_Yes_Hyp_No), labs=Protein, cex=.4))
