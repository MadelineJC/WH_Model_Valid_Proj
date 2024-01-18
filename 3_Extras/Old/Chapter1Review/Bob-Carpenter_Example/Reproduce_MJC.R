# Read in data
lynx_hare_df <- read.csv("3_Extras/Old/Chapter1Review/Bob-Carpenter_Example/hudson-bay-lynx-hare.csv")

# Plotting
plot(lynx_hare_df$Year, lynx_hare_df$Lynx, type = "l", xlab = "Year", ylab = "Population Abundance", ylim = c(0, max(lynx_hare_df$Hare)), col = "firebrick2", lwd = 2)
lines(lynx_hare_df$Year, lynx_hare_df$Hare, col = "orchid", lwd = 2)
legend("topright", legend = c("Lynx", "Hare"), col = c("firebrick2", "orchid"), lwd = 2)

