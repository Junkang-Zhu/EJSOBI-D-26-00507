library(tidyverse)

# 读取数据
data <- read.csv("nematode_data.csv")

# 按方法拆分
methods <- unique(data$Method)
families <- colnames(data)[!colnames(data) %in% c("SampleID", "Method")]

# 定义配对
pairs <- list(
  list(test = "Nematode-extracted DNA", ref = "Morphological identification"),
  list(test = "Nematode-extracted DNA", ref = "Mass proportion"),
  list(test = "eDNA metabarcoding", ref = "Morphological identification"),
  list(test = "eDNA metabarcoding", ref = "Mass proportion")
)

# 计算每个配对的差值和均值
ba_data <- map_dfr(pairs, function(p) {
  test_dat <- data %>% filter(Method == p$test) %>% select(-SampleID, -Method)
  ref_dat  <- data %>% filter(Method == p$ref)  %>% select(-SampleID, -Method)
  
  map_dfr(families, function(fam) {
    t_vec <- test_dat[[fam]]
    r_vec <- ref_dat[[fam]]
    data.frame(
      Family = fam,
      Pair   = paste0(p$test, " vs ", p$ref),
      Test   = p$test,
      Ref    = p$ref,
      Mean   = (t_vec + r_vec) / 2,
      Diff   = t_vec - r_vec
    )
  })
}) %>% filter(!is.na(Mean) & !is.na(Diff))

# ---- 计算每个配对的LOA及其Bootstrap CI ----
set.seed(123)
boot_loa_limits <- function(diffs, B = 2000) {
  n <- length(diffs)
  boot_lower <- replicate(B, {
    s <- sample(diffs, n, replace = TRUE)
    quantile(s, 0.025)
  })
  boot_upper <- replicate(B, {
    s <- sample(diffs, n, replace = TRUE)
    quantile(s, 0.975)
  })
  data.frame(
    lower_loa_ci_lower = quantile(boot_lower, 0.025),
    lower_loa_ci_upper = quantile(boot_lower, 0.975),
    upper_loa_ci_lower = quantile(boot_upper, 0.025),
    upper_loa_ci_upper = quantile(boot_upper, 0.975)
  )
}

loa_summary <- ba_data %>%
  group_by(Pair) %>%
  summarise(
    n = n(),
    median_diff = median(Diff),
    lower_loa = quantile(Diff, 0.025),
    upper_loa = quantile(Diff, 0.975),
    loa_width = upper_loa - lower_loa,
    .groups = "drop"
  )

# 为每个配对添加LOA的Bootstrap CI
loa_ci <- ba_data %>%
  group_by(Pair) %>%
  summarise(boot_ci = list(boot_loa_limits(Diff)), .groups = "drop") %>%
  unnest_wider(boot_ci)

loa_summary <- loa_summary %>% left_join(loa_ci, by = "Pair")

# ---- 输出完整表格 ----
loa_display <- loa_summary %>%
  mutate(
    LOA_95 = paste0("[", round(lower_loa, 4), ", ", round(upper_loa, 4), "]"),
    lower_LOA_CI = paste0("[", round(lower_loa_ci_lower, 4), ", ", round(lower_loa_ci_upper, 4), "]"),
    upper_LOA_CI = paste0("[", round(upper_loa_ci_lower, 4), ", ", round(upper_loa_ci_upper, 4), "]"),
    Width = round(loa_width, 4)
  ) %>%
  select(Pair, n, Median_Diff = median_diff, LOA_95, lower_LOA_CI, upper_LOA_CI, Width)

print(loa_display)
write.csv(loa_display, "ba_loa_summary_with_CI.csv", row.names = FALSE)

# ---- Bootstrap检验LOA宽度差异 ----
boot_loa_width <- function(diffs, B = 2000) {
  n <- length(diffs)
  replicate(B, {
    s <- sample(diffs, n, replace = TRUE)
    quantile(s, 0.975) - quantile(s, 0.025)
  })
}

dna_morpho <- ba_data %>% filter(Pair == "Nematode-extracted DNA vs Morphological identification")
dna_mass   <- ba_data %>% filter(Pair == "Nematode-extracted DNA vs Mass proportion")
boot_dna_morpho <- boot_loa_width(dna_morpho$Diff)
boot_dna_mass   <- boot_loa_width(dna_mass$Diff)
p_dna <- 2 * min(mean(boot_dna_morpho > boot_dna_mass), mean(boot_dna_morpho < boot_dna_mass))

edna_morpho <- ba_data %>% filter(Pair == "eDNA metabarcoding vs Morphological identification")
edna_mass   <- ba_data %>% filter(Pair == "eDNA metabarcoding vs Mass proportion")
boot_edna_morpho <- boot_loa_width(edna_morpho$Diff)
boot_edna_mass   <- boot_loa_width(edna_mass$Diff)
p_edna <- 2 * min(mean(boot_edna_morpho > boot_edna_mass), mean(boot_edna_morpho < boot_edna_mass))

sig_summary <- data.frame(
  Comparison = c("DNA: Morpho vs Mass", "eDNA: Morpho vs Mass"),
  LOA_width_ref1 = c(loa_summary$loa_width[1], loa_summary$loa_width[3]),
  LOA_width_ref2 = c(loa_summary$loa_width[2], loa_summary$loa_width[4]),
  p_value = c(p_dna, p_edna)
)
write.csv(sig_summary, "ba_significance_summary.csv", row.names = FALSE)

# ---- 绝对残差比较 ----
cat("\n=== 绝对残差比较（Mann-Whitney U检验）===\n")
dna_morpho_abs <- abs(ba_data$Diff[ba_data$Pair == "Nematode-extracted DNA vs Morphological identification"])
dna_mass_abs   <- abs(ba_data$Diff[ba_data$Pair == "Nematode-extracted DNA vs Mass proportion"])
wt_dna <- wilcox.test(dna_morpho_abs, dna_mass_abs, paired = FALSE, exact = FALSE)
cat(sprintf("DNA: Morpho中位数=%.4f, Mass中位数=%.4f, p=%.4f\n", 
            median(dna_morpho_abs), median(dna_mass_abs), wt_dna$p.value))

edna_morpho_abs <- abs(ba_data$Diff[ba_data$Pair == "eDNA metabarcoding vs Morphological identification"])
edna_mass_abs   <- abs(ba_data$Diff[ba_data$Pair == "eDNA metabarcoding vs Mass proportion"])
wt_edna <- wilcox.test(edna_morpho_abs, edna_mass_abs, paired = FALSE, exact = FALSE)
cat(sprintf("eDNA: Morpho中位数=%.4f, Mass中位数=%.4f, p=%.4f\n", 
            median(edna_morpho_abs), median(edna_mass_abs), wt_edna$p.value))

# ---- 导出数据用于Origin作图 ----
write.csv(ba_data, "ba_data_for_origin.csv", row.names = FALSE)

# ---- 按科汇总 ----
loa_by_family <- ba_data %>%
  group_by(Pair, Family) %>%
  summarise(
    n = n(),
    median_diff = median(Diff),
    lower_loa = quantile(Diff, 0.025),
    upper_loa = quantile(Diff, 0.975),
    loa_width = upper_loa - lower_loa,
    .groups = "drop"
  )
write.csv(loa_by_family, "ba_loa_by_family.csv", row.names = FALSE)

cat("\n文件已导出: ba_loa_summary_with_CI.csv, ba_significance_summary.csv, ba_data_for_origin.csv, ba_loa_by_family.csv\n")






