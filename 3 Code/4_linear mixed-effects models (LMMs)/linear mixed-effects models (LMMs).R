library(glmmTMB)
library(DHARMa)

df <- read.csv("nematode_data.csv")
df[c("SampleID", "Family")] <- lapply(df[c("SampleID", "Family")], as.factor)
df <- df[complete.cases(df), ]

# ===== 关键：处理MeanMass中的0值（替换而非删除）=====
min_mass <- min(df$MeanMass[df$MeanMass > 0], na.rm = TRUE)
df$MeanMass[df$MeanMass <= 0] <- min_mass / 2
df$Abundance <- pmax(df$Abundance, 0)

eps <- 0.0001
df$DNA_rel <- pmin(pmax(df$DNA_rel, eps), 1 - eps)

model <- glmmTMB(DNA_rel ~ log(Abundance + 1) + log(MeanMass) + (1 | SampleID) + (1 | Family),
                 data = df, family = beta_family())

summary(model)
plot(simulateResiduals(model))

fixef <- summary(model)$coefficients$cond
write.csv(data.frame(Variable = rownames(fixef), fixef, check.names = FALSE), 
          "fixed_effects.csv", row.names = FALSE)

var_df <- data.frame(Group = names(VarCorr(model)$cond), 
                     Variance = as.numeric(VarCorr(model)$cond))
write.csv(transform(var_df, StdDev = sqrt(Variance)), "random_effects.csv", row.names = FALSE)

df$Fitted <- fitted(model)
df$Residuals <- residuals(model, type = "response")
write.csv(df, "data_with_fitted.csv", row.names = FALSE)

mass_p <- summary(model)$coefficients$cond["log(MeanMass)", "Pr(>|z|)"]
cat(sprintf("AIC: %.1f | BIC: %.1f | log(MeanMass) p = %.4f %s\n", 
            AIC(model), BIC(model), mass_p, ifelse(mass_p < 0.05, "*", "ns")))





library(ggplot2)

# 读取数据
fixef <- read.csv("fixed_effects.csv")

# 查看列名（调试用，可删除）
print(colnames(fixef))

# 用列索引提取数据（简单可靠，无需担心列名细节）
fixef <- fixef[fixef[, 1] != "(Intercept)", ]  # 去掉截距（第1列是Variable）
est <- fixef[, 2]                              # 第2列是Estimate
se <- fixef[, 3]                               # 第3列是Std.Error

# 计算置信区间
ci_lower <- est - 1.96 * se
ci_upper <- est + 1.96 * se

# 构建绘图数据
plot_df <- data.frame(
  Variable = c("Abundance", "Mean mass"),
  Estimate = est,
  CI_lower = ci_lower,
  CI_upper = ci_upper
)

# 绘图
p <- ggplot(plot_df, aes(x = Estimate, y = Variable)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.6) +
  geom_errorbarh(aes(xmin = CI_lower, xmax = CI_upper), height = 0.15, color = "black", linewidth = 0.8) +
  geom_point(size = 4, color = "#2b8cbe") +
  labs(x = "Effect size (Estimate ± 95% CI)", 
       y = "",
       title = "Fixed effects on DNA relative abundance") +
  xlim(-0.5, 1.8) +
  theme_minimal(base_size = 13) +
  theme(axis.text.y = element_text(size = 13, face = "bold"),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank())

print(p)
ggsave("forest_plot.png", p, width = 6, height = 3, dpi = 300)
