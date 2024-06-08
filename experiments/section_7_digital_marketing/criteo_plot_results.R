library(grf)

seed = 20211013

q = seq(0.001, 1, by = 0.003)
load("rates_and_tocs.Rdata")

pdf(file = "toc_comparison_visit.pdf")
par(mar = c(5, 6, 2, 1) + .1)
plot(
  autoc_visit,
  main = NULL,
  sub = NULL,
  ylim = c(-0.05, 0.2),
  lwd = 2,
  ylab = "Effect on Visit Rate",
  xlab = "Fraction Treated",
  cex.lab = 2,
  cex.axis = 1.5
)
legend(
  0.5,
  0.2,
  legend = c("Baseline", "CATE"),
  col = c(1, 2),
  lwd = 2,
  cex = 2
)
dev.off()
pdf(file = "toc_comparison_conversion.pdf")
par(mar = c(5, 6, 2, 1) + .1)
plot(
  autoc_conversion,
  main = NULL,
  sub = NULL,
  ylim = c(-0.01, 0.06),
  lwd = 2,
  ylab = "Effect on Conversion Rate",
  xlab = "Fraction Treated",
  cex.lab = 2,
  cex.axis = 1.5
)
legend(
  0.5,
  0.06,
  legend = c("Baseline", "CATE"),
  col = c(1, 2),
  lwd = 2,
  cex = 2
)
dev.off()

print("AUTOC Visit")
print(autoc_visit)
print("AUTOC Conversion")
print(autoc_conversion)

print("Qini visit")
print(qini_visit)
print("Qini Conversion")
print(qini_conversion)
