
##CN A+ CVD- vs CN A-
dev.off()
with(resultsproteomics_CVD, plot(diff_CN_A_Yes_CVD_No_vs_CN_A_No, -log10(p_CN_A_Yes_CVD_No_vs_CN_A_No), pch=20, col="grey41", main="Volcano plot", xlim=c(-1.5, 1.5), ylim=c(0, 4)))

# Add colored points
with(subset(resultsproteomics_CVD, p_CN_A_Yes_CVD_No_vs_CN_A_No<.05 ), points(diff_CN_A_Yes_CVD_No_vs_CN_A_No, -log10(p_CN_A_Yes_CVD_No_vs_CN_A_No), pch=20, col="tomato3"))

# Label and color the top 20 proteins
library(calibrate)
with(subset(resultsproteomics_CVD, p_CN_A_Yes_CVD_No_vs_CN_A_No<.01506), points(diff_CN_A_Yes_CVD_No_vs_CN_A_No, -log10(p_CN_A_Yes_CVD_No_vs_CN_A_No), pch=20, col="dodgerblue3"))
with(subset(resultsproteomics_CVD, p_CN_A_Yes_CVD_No_vs_CN_A_No<.01506), textxy(diff_CN_A_Yes_CVD_No_vs_CN_A_No, -log10(p_CN_A_Yes_CVD_No_vs_CN_A_No), labs=Protein, cex=.4))

##CN A+ CVD+ vs CN A-
dev.off()
with(resultsproteomics_CVD, plot(diff_CN_A_Yes_CVD_Yes_vs_CN_A_No, -log10(p_CN_A_Yes_CVD_Yes_vs_CN_A_No), pch=20, col="grey41", main="Volcano plot", xlim=c(-1.6, 1.6), ylim=c(0, 2.5)))

# Add colored points
with(subset(resultsproteomics_CVD, p_CN_A_Yes_CVD_Yes_vs_CN_A_No<.05), points(diff_CN_A_Yes_CVD_Yes_vs_CN_A_No, -log10(p_CN_A_Yes_CVD_Yes_vs_CN_A_No), pch=20, col="tomato3"))

# Label and color the top 20 proteins
library(calibrate)
with(subset(resultsproteomics_CVD, p_CN_A_Yes_CVD_Yes_vs_CN_A_No<.05), points(diff_CN_A_Yes_CVD_Yes_vs_CN_A_No, -log10(p_CN_A_Yes_CVD_Yes_vs_CN_A_No), pch=20, col="dodgerblue3"))
with(subset(resultsproteomics_CVD, p_CN_A_Yes_CVD_Yes_vs_CN_A_No<.05), textxy(diff_CN_A_Yes_CVD_Yes_vs_CN_A_No, -log10(p_CN_A_Yes_CVD_Yes_vs_CN_A_No), labs=Protein, cex=.4))

##MCI A+ CVD- vs CN A-
dev.off()
with(resultsproteomics_CVD, plot(diff_MCI_A_Yes_CVD_No_vs_CN_A_No, -log10(p_MCI_A_Yes_CVD_No_vs_CN_A_No), pch=20, col="grey41", main="Volcano plot", xlim=c(-1.6, 1.6), ylim=c(0, 5)))

# Add colored points
with(subset(resultsproteomics_CVD, p_MCI_A_Yes_CVD_No_vs_CN_A_No<.05), points(diff_MCI_A_Yes_CVD_No_vs_CN_A_No, -log10(p_MCI_A_Yes_CVD_No_vs_CN_A_No), pch=20, col="tomato3"))

# Label and color the top 20 proteins
library(calibrate)
with(subset(resultsproteomics_CVD, p_MCI_A_Yes_CVD_No_vs_CN_A_No<.00023), points(diff_MCI_A_Yes_CVD_No_vs_CN_A_No, -log10(p_MCI_A_Yes_CVD_No_vs_CN_A_No), pch=20, col="dodgerblue3"))
with(subset(resultsproteomics_CVD, p_MCI_A_Yes_CVD_No_vs_CN_A_No<.00023), textxy(diff_MCI_A_Yes_CVD_No_vs_CN_A_No, -log10(p_MCI_A_Yes_CVD_No_vs_CN_A_No), labs=Protein, cex=.4))

##MCI A+ CVD+ vs CN A-
dev.off()
with(resultsproteomics_CVD, plot(diff_MCI_A_Yes_CVD_Yes_vs_CN_A_No, -log10(p_MCI_A_Yes_CVD_Yes_vs_CN_A_No), pch=20, col="grey41", main="Volcano plot", xlim=c(-2.6, 2.6), ylim=c(0, 5)))

# Add colored points
with(subset(resultsproteomics_CVD, p_MCI_A_Yes_CVD_Yes_vs_CN_A_No<.05), points(diff_MCI_A_Yes_CVD_Yes_vs_CN_A_No, -log10(p_MCI_A_Yes_CVD_Yes_vs_CN_A_No), pch=20, col="tomato3"))

# Label and color the top 20 proteins
library(calibrate)
with(subset(resultsproteomics_CVD, p_MCI_A_Yes_CVD_Yes_vs_CN_A_No<.00048), points(diff_MCI_A_Yes_CVD_Yes_vs_CN_A_No, -log10(p_MCI_A_Yes_CVD_Yes_vs_CN_A_No), pch=20, col="dodgerblue3"))
with(subset(resultsproteomics_CVD, p_MCI_A_Yes_CVD_Yes_vs_CN_A_No<.00048), textxy(diff_MCI_A_Yes_CVD_Yes_vs_CN_A_No, -log10(p_MCI_A_Yes_CVD_Yes_vs_CN_A_No), labs=Protein, cex=.4))

##AD A+ CVD- vs CN A-
dev.off()
with(resultsproteomics_CVD, plot(diff_AD_A_Yes_CVD_No_vs_CN_A_No, -log10(p_AD_A_Yes_CVD_No_vs_CN_A_No), pch=20, col="grey41", main="Volcano plot", xlim=c(-1.5, 1.5), ylim=c(0, 5)))

# Add colored points
with(subset(resultsproteomics_CVD, p_AD_A_Yes_CVD_No_vs_CN_A_No<.05), points(diff_AD_A_Yes_CVD_No_vs_CN_A_No, -log10(p_AD_A_Yes_CVD_No_vs_CN_A_No), pch=20, col="tomato3"))

# Label and color the top 20 proteins
library(calibrate)
with(subset(resultsproteomics_CVD, p_AD_A_Yes_CVD_No_vs_CN_A_No<.01158), points(diff_AD_A_Yes_CVD_No_vs_CN_A_No, -log10(p_AD_A_Yes_CVD_No_vs_CN_A_No), pch=20, col="dodgerblue3"))
with(subset(resultsproteomics_CVD, p_AD_A_Yes_CVD_No_vs_CN_A_No<.01158), textxy(diff_AD_A_Yes_CVD_No_vs_CN_A_No, -log10(p_AD_A_Yes_CVD_No_vs_CN_A_No), labs=Protein, cex=.4))

##AD A+ CVD+ vs CN A-
dev.off()
with(resultsproteomics_CVD, plot(diff_AD_A_Yes_CVD_Yes_vs_CN_A_No, -log10(p_AD_A_Yes_CVD_Yes_vs_CN_A_No), pch=20, col="grey41", main="Volcano plot", xlim=c(-2, 2), ylim=c(0, 3.5)))

# Add colored points
with(subset(resultsproteomics_CVD, p_AD_A_Yes_CVD_Yes_vs_CN_A_No<.05), points(diff_AD_A_Yes_CVD_Yes_vs_CN_A_No, -log10(p_AD_A_Yes_CVD_Yes_vs_CN_A_No), pch=20, col="tomato3"))

# Label and color the top 20 proteins
library(calibrate)
with(subset(resultsproteomics_CVD, p_AD_A_Yes_CVD_Yes_vs_CN_A_No<.02701), points(diff_AD_A_Yes_CVD_Yes_vs_CN_A_No, -log10(p_AD_A_Yes_CVD_Yes_vs_CN_A_No), pch=20, col="dodgerblue3"))
with(subset(resultsproteomics_CVD, p_AD_A_Yes_CVD_Yes_vs_CN_A_No<.02701), textxy(diff_AD_A_Yes_CVD_Yes_vs_CN_A_No, -log10(p_AD_A_Yes_CVD_Yes_vs_CN_A_No), labs=Protein, cex=.4))

##CN A+ CVD+ vs CN A+ CVD-
dev.off()
with(resultsproteomics_CVD, plot(diff_CN_A_Yes_CVD_Yes_vs_CN_A_Yes_CVD_No, -log10(p_CN_A_Yes_CVD_Yes_vs_CN_A_Yes_CVD_No), pch=20, col="grey41", main="Volcano plot", xlim=c(-2.2, 2.2), ylim=c(0, 3.5)))

# Add colored points
with(subset(resultsproteomics_CVD, p_CN_A_Yes_CVD_Yes_vs_CN_A_Yes_CVD_No<.05), points(diff_CN_A_Yes_CVD_Yes_vs_CN_A_Yes_CVD_No, -log10(p_CN_A_Yes_CVD_Yes_vs_CN_A_Yes_CVD_No), pch=20, col="tomato3"))

# Label and color the top 20 proteins
library(calibrate)
with(subset(resultsproteomics_CVD, p_CN_A_Yes_CVD_Yes_vs_CN_A_Yes_CVD_No<.02837), points(diff_CN_A_Yes_CVD_Yes_vs_CN_A_Yes_CVD_No, -log10(p_CN_A_Yes_CVD_Yes_vs_CN_A_Yes_CVD_No), pch=20, col="dodgerblue3"))
with(subset(resultsproteomics_CVD, p_CN_A_Yes_CVD_Yes_vs_CN_A_Yes_CVD_No<.02837), textxy(diff_CN_A_Yes_CVD_Yes_vs_CN_A_Yes_CVD_No, -log10(p_CN_A_Yes_CVD_Yes_vs_CN_A_Yes_CVD_No), labs=Protein, cex=.4))

##MCI A+ CVD+ vs MCI A+ CVD-
dev.off()
with(resultsproteomics_CVD, plot(diff_MCI_A_Yes_CVD_Yes_vs_MCI_A_Yes_CVD_No, -log10(p_MCI_A_Yes_CVD_Yes_vs_MCI_A_Yes_CVD_No), pch=20, col="grey41", main="Volcano plot", xlim=c(-1.5, 1.5), ylim=c(0, 4)))

# Add colored points
with(subset(resultsproteomics_CVD, p_MCI_A_Yes_CVD_Yes_vs_MCI_A_Yes_CVD_No<.05), points(diff_MCI_A_Yes_CVD_Yes_vs_MCI_A_Yes_CVD_No, -log10(p_MCI_A_Yes_CVD_Yes_vs_MCI_A_Yes_CVD_No), pch=20, col="tomato3"))

# Label and color the top 20 proteins
library(calibrate)
with(subset(resultsproteomics_CVD, p_MCI_A_Yes_CVD_Yes_vs_MCI_A_Yes_CVD_No<.01272), points(diff_MCI_A_Yes_CVD_Yes_vs_MCI_A_Yes_CVD_No, -log10(p_MCI_A_Yes_CVD_Yes_vs_MCI_A_Yes_CVD_No), pch=20, col="dodgerblue3"))
with(subset(resultsproteomics_CVD, p_MCI_A_Yes_CVD_Yes_vs_MCI_A_Yes_CVD_No<.01272), textxy(diff_MCI_A_Yes_CVD_Yes_vs_MCI_A_Yes_CVD_No, -log10(p_MCI_A_Yes_CVD_Yes_vs_MCI_A_Yes_CVD_No), labs=Protein, cex=.4))

##AD A+ CVD+ vs AD A+ CVD-
dev.off()
with(resultsproteomics_CVD, plot(diff_AD_A_Yes_CVD_Yes_vs_AD_A_Yes_CVD_No, -log10(p_AD_A_Yes_CVD_Yes_vs_AD_A_Yes_CVD_No), pch=20, col="grey41", main="Volcano plot", xlim=c(-2.2, 2.2), ylim=c(0, 4)))

# Add colored points
with(subset(resultsproteomics_CVD, p_AD_A_Yes_CVD_Yes_vs_AD_A_Yes_CVD_No<.05), points(diff_AD_A_Yes_CVD_Yes_vs_AD_A_Yes_CVD_No, -log10(p_AD_A_Yes_CVD_Yes_vs_AD_A_Yes_CVD_No), pch=20, col="tomato3"))

# Label and color the top 20 proteins
library(calibrate)
with(subset(resultsproteomics_CVD, p_AD_A_Yes_CVD_Yes_vs_AD_A_Yes_CVD_No<.03342), points(diff_AD_A_Yes_CVD_Yes_vs_AD_A_Yes_CVD_No, -log10(p_AD_A_Yes_CVD_Yes_vs_AD_A_Yes_CVD_No), pch=20, col="dodgerblue3"))
with(subset(resultsproteomics_CVD, p_AD_A_Yes_CVD_Yes_vs_AD_A_Yes_CVD_No<.03342), textxy(diff_AD_A_Yes_CVD_Yes_vs_AD_A_Yes_CVD_No, -log10(p_AD_A_Yes_CVD_Yes_vs_AD_A_Yes_CVD_No), labs=Protein, cex=.4))
