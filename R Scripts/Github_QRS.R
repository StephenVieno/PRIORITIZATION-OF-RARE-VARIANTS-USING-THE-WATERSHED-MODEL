# Create QRS code 
# install.packages("qrcode", repos = "https://thierryo.r-universe.dev")
library(qrcode)
# Github URL 
code <- qr_code("https://github.com/StephenVieno/PRIORITIZATION-OF-RARE-VARIANTS-USING-THE-WATERSHED-MODEL")
# Plot QRS
plot(code)

# Save QRS 
generate_svg(code, filename = "Plots/Github_QRS.svg")
pdf("Plots/Github_QRS.pdf", width = 7, height = 7)
