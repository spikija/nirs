# analysis of nirs cohort per Collette 2022
# LIBRARIES
# ----
pacman::p_load(
  ggplot2,
  readxl,
  Hmisc,
  gtsummary,
  summarytools,
  kableExtra,
  expss,
  rmarkdown,
  skimr,        # get overview of data
  tidyverse,    # data management + ggplot2 graphics 
  gtsummary,    # summary statistics and tests
  rstatix,      # summary statistics and statistical tests
  janitor,      # adding totals and percents to tables
  scales,       # easily convert proportions to percents  
  flextable,     # converting tables to pretty images
  ggstatsplot,
  palmerpenguins,
  survival,
  DiagrammeR,
  glue,
  igraph,
  ggraph,
  progressr,
  zoo,
  changepoint,
  ggpubr,
  MASS,
  mice,
  nnet,
  DescTools,
  officer,
  dplyr,
  pROC,
  patchwork
)
# ----
# LIBRARIES ENDE


# call usual functions
db.source.file.path <- ifelse(Sys.info()[[4]] == "LAPTOP-NFRMGQDO" | Sys.info()[[4]] == "NEURONODE", "c:/sci/cdk.stat/", "C:/app/r.apps/cdk.stat/")
source(paste0(db.source.file.path, "sp.utilities.r"))
source(paste0(db.source.file.path, "nirs_helper.r"))

# ----

db.language <- ifelse(Sys.info()[[4]] == "LAPTOP-NFRMGQDO" | Sys.info()[[4]] == "NEURONODE", "Meine Ablage", "My Drive")
db.path <- paste0("C:/sci/rcode/nirs/")

# DBDBDBDBDBDBDBDBDBDBDBDBDBDBDBDBDBDBDBDBDBDBDBDBDBDBDBDBDB
# LAST STABLE DB
db_clin <- read.csv2(paste0(db.path, "clinical_170_clean_v1.csv"))

# load main database
db_msr <- read.csv2(paste0(db.path, "nirs_measurements_clean_v2.csv"))

# methodology from Hametner et al. 2015:
# We calculated interhemispheric rSO2 differences before and after the intervention as the difference in the
# averaged 5-minute rSO2 values of the affected hemisphere minus values in
# the nonaffected hemisphere, yielding positive or negative values

# reviewer
# Evaluate and display individual affected-hemisphere identification. 
# Report the number and percentage of patients in whom rSO₂ was lower, equal or higher 
# on the affected side. 
# Show the raw affected–unaffected rSO₂ difference for every patient 
# in a dot plot, with zero indicating no asymmetry. 
# In addition, report how often the hemisphere with the lower rSO₂ value 
# correctly corresponded to the angiographically affected hemisphere, 
# including correct, incorrect and indeterminate classifications with 95% CI. 


# find 5 minute intervals for following points:
# 1. between 10 and 15 minutes after start of recording
# - calculate average for each hemisphere
# - take 4% difference as cut-off (Annus et al.)
# 2. after 10 minutes, and then 5 minutes period after revascularisation

# Take 2% as significant difference between hemispheres (Damian et al 2007, Annus et al 2019)
# Take 4% as significant diffeence for hemisphere

# INPUT VARIABLES

# initial delay in minutes after start of procedure
init_delay <- 10

# duration in minutes upon which rSO2 will be averages
sample_duration <- 10

# duration in minutes after recanalisation end after which the
# average will be calculated
recan_delay <- 5 

# INPUT VARIABLES - END

# MAIN VARIABLES
# avg_start - average of raw rSO2 values 10 minutes after begin,
# for sample_duration duration in seconds

# avg_post - average of raw rSO2 values 5 minutes after recanalisation
# for sample_duration in seconds

avg_start <- db_msr %>% group_by(patient_id, hemisphere_status) %>% 
  filter(elapsed_seconds >= 60*init_delay & elapsed_seconds <= 60*(init_delay+sample_duration)) %>% 
  summarise(avg = mean(nirs_value_raw, na.rm = TRUE)) %>% 
  pivot_wider(names_from = hemisphere_status, values_from = avg) %>% 
  # detect 4% difference
  mutate(
    difference = Affected - Unaffected,
    hemisphere_ratio = (abs(Affected - Unaffected) / Unaffected) * 100,
    ratio_signif = case_when(
      hemisphere_ratio >= 4 & Affected > Unaffected ~ "Affected higher", 
      hemisphere_ratio >= 4 & Affected < Unaffected ~ "Unaffected higher",
      TRUE ~ "no")
  )

avg_start

# ------------------------------------------------------------------
# Affected-hemisphere identification
# ------------------------------------------------------------------

tb1 <- table(avg_start$ratio_signif)

# Counts
n_aff_higher   <- unname(tb1["Affected higher"])
n_equal        <- unname(tb1["no"])
n_unaff_higher <- unname(tb1["Unaffected higher"])
n_total        <- sum(tb1)

# Percentages
pct_aff_higher   <- 100 * n_aff_higher / n_total
pct_equal        <- 100 * n_equal / n_total
pct_unaff_higher <- 100 * n_unaff_higher / n_total

# Overall 3-category distribution
chi_test <- chisq.test(tb1)

# Among patients with directional asymmetry:
# lower rSO2 correctly identifies affected hemisphere
binom_test <- binom.test(
  x = n_unaff_higher,                     # correct
  n = n_unaff_higher + n_aff_higher,      # determinate cases
  p = 0.5
)

# Useful values for RMarkdown
n_determinate <- n_unaff_higher + n_aff_higher

ident_accuracy <- 100 * unname(binom_test$estimate)
ident_ci_low   <- 100 * binom_test$conf.int[1]
ident_ci_high  <- 100 * binom_test$conf.int[2]

chi_p   <- chi_test$p.value
binom_p <- binom_test$p.value
# so with these tests we have answered at least one question:

# the frequency of affected hemisphere having lower rSO2 is significant
# p = 0.041 (95% CI 0.30 - 0.49)

# producing dot plot
# graph will show differences between hemispheres in the 5 minute window 10 minutes
# after procedure start
# ----
# proof of concept
fig1 <- avg_start %>% 
  ggplot(aes(x = 0, y = difference)) + 
  geom_dotplot(binwidth = 0.5, dotsize = 1.5, binaxis = "y", stackdir = "center") 
fig1 

# production graph
# horizontal layout positions
x_dot  <- -0.75
x_text <-  -0.35

fig1 <- avg_start %>%
  ggplot(aes(x = 0, y = difference)) +
  
  # shaded similarity zone
  annotate(
    "rect",
    xmin = -Inf, xmax = Inf,
    ymin = -4, ymax = 4,
    alpha = 0.10,
    fill = "grey60"
  ) +
  
  # patient-level dot distribution
  geom_dotplot(
    binwidth = 0.5,
    dotsize = 1.5,
    binaxis = "y",
    stackdir = "center",
    position = position_nudge(x = x_dot)
  ) +
  
  # ±4 percentage-point thresholds
  geom_hline(
    yintercept = c(-4, 4),
    linetype = "dotted",
    linewidth = 0.45
  ) +
  
  # annotations
  annotate(
    "text",
    x = x_text,
    y = 12,
    label = paste0("Higher on\naffected side\nN = ",n_aff_higher," (", 
                   format(pct_aff_higher, digits=2), "%)"),
    hjust = 0,
    size = 4
  ) +
  
  # 1. Heading
  annotate(
    "text",
    x = x_text,
    y = 3.0,
    label = "Similar between hemispheres",
    hjust = 0,
    size = 4
  ) +
  
  # 2. Delta-rSO2 criterion
  annotate(
    "text",
    x = x_text,
    y = 0.0,
    label = "group('|', Delta*rSO[2], '|') <= 4~'percentage points'",
    parse = TRUE,
    hjust = 0,
    size = 4
  ) +
  
  # 3. N and percentage
  annotate(
    "text",
    x = x_text,
    y = -3.0,
    label = paste0(
      "N = ", n_equal,
      " (", format(pct_equal, digits = 2), "%)"
    ),
    hjust = 0,
    size = 4
  )+
  annotate(
    "text",
    x = x_text,
    y = -12,
    label = paste0("Higher on\nunaffected side\nN = ",n_unaff_higher," (", 
                   format(pct_unaff_higher, digits=2), "%)"),,
    hjust = 0,
    size = 4
  ) +
  
  labs(
    x = NULL,
    y = expression(Delta*rSO[2]~"(percentage points)")
  ) +
  
  # THIS is the important part:
  # explicitly control the meaningless horizontal dimension
  scale_x_continuous(
    limits = c(-1.25, 0.65),
    expand = expansion(mult = 0)
  ) +
  
  scale_y_continuous(
    breaks = c(-20, -10, -4, 0, 4, 10, 20),
    limits=c(-22, 22),
    expand = expansion(mult = c(0.03, 0.05))
  ) +
  
  theme_classic(base_size = 12) +
  
  theme(
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x  = element_blank(),
    panel.grid   = element_blank(),
    
    axis.title.y = element_text(size = 13),
    axis.text.y  = element_text(size = 11)
  )

fig1
# ----
# Block for Q1


# Q2 : 

avg_2 <- db_msr %>% group_by(patient_id, elapsed_seconds) %>% 
  filter(recanalisation_point == TRUE) %>% dplyr::select(elapsed_seconds) %>% 
  mutate(start_value = elapsed_seconds + (60*(recan_delay)), 
         upper_value = elapsed_seconds + (60*(recan_delay + sample_duration)))
avg_2 


# calculate average rSO2 raw values based on cut-off (sample_duration)
calc_avg <- function(db, db_search) {
  print(db_search$patient_id)
  db <- db %>% filter(db$patient_id == db_search$patient_id)
  db_row <- db %>% group_by(patient_id, hemisphere_status) %>% 
    filter(elapsed_seconds >= db_search$elapsed_seconds & elapsed_seconds <= 
             db_search$upper_value) %>% 
     summarise(avg = mean(nirs_value_raw, na.rm = TRUE)) %>% 
     pivot_wider(names_from = hemisphere_status, values_from = avg) %>% 
     # detect 4% difference
     mutate(
       difference = Affected - Unaffected,
       hemisphere_ratio = (abs(Affected - Unaffected) / Unaffected) * 100,
       ratio_signif = case_when(
         hemisphere_ratio >= 4 & Affected > Unaffected ~ "Affected higher", 
         hemisphere_ratio >= 4 & Affected < Unaffected ~ "Unaffected higher",
         TRUE ~ "no")
     )
  return(db_row)
}

avg_post <- tibble()

for (i in 1:nrow(avg_2)) { 
  avg_post <- bind_rows(avg_post, (calc_avg(db_msr, avg_2[i,])))
}


# join avg_pre and avg_post
# change column names
avg_s1 <- avg_start[,c("patient_id", "Affected", "Unaffected")]
colnames(avg_s1) <- c("patient_id", "Affected_start", "Unaffected_start")
avg_s2 <- avg_post[,c("patient_id", "Affected", "Unaffected")]
colnames(avg_s2) <- c("patient_id", "Affected_post", "Unaffected_post")

# inner join (the not recanalised patients will not be shown!)
avg_graph <- inner_join(avg_s2, avg_s1, by="patient_id")

# make line-dot graph
# spread long
avg_graph_pivot <- avg_graph |> 
  pivot_longer(!patient_id, names_to=c("hemisphere", "period"), names_sep = "_", 
               values_to="avg_10min_rSO2")

avg_graph_pivot

str(avg_graph_pivot)

avg_graph_pivot$period <- factor(avg_graph_pivot$period, c("start", "post"))

str(avg_graph_pivot)

# RESULTS

# ------------------------------------------------------------
# Paired start vs post analysis
# ------------------------------------------------------------

# ------------------------------------------------------------
# Prepare paired start/post data
# ------------------------------------------------------------

prepost_wide <- avg_graph_pivot |>
  dplyr::select(patient_id, hemisphere, period, avg_10min_rSO2) |>
  pivot_wider(
    names_from = period,
    values_from = avg_10min_rSO2
  ) |>
  filter(!is.na(start), !is.na(post)) |>
  mutate(
    change = post - start
  )

# ------------------------------------------------------------
# Wilcoxon signed-rank test + Hodges-Lehmann estimate
# ------------------------------------------------------------

wilcox_results <- prepost_wide |>
  group_by(hemisphere) |>
  group_modify(~ {
    
    wt <- wilcox.test(
      .x$post,
      .x$start,
      paired = TRUE,
      conf.int = TRUE,
      conf.level = 0.95,
      exact = FALSE
    )
    
    tibble(
      n = nrow(.x),
      
      start_mean   = mean(.x$start, na.rm = TRUE),
      start_sd     = sd(.x$start, na.rm = TRUE),
      
      post_mean    = mean(.x$post, na.rm = TRUE),
      post_sd      = sd(.x$post, na.rm = TRUE),
      
      start_median = median(.x$start, na.rm = TRUE),
      start_q1     = quantile(.x$start, 0.25, na.rm = TRUE),
      start_q3     = quantile(.x$start, 0.75, na.rm = TRUE),
      
      post_median  = median(.x$post, na.rm = TRUE),
      post_q1      = quantile(.x$post, 0.25, na.rm = TRUE),
      post_q3      = quantile(.x$post, 0.75, na.rm = TRUE),
      
      median_change = median(.x$change, na.rm = TRUE),
      
      HL_estimate = unname(wt$estimate),
      ci_low      = wt$conf.int[1],
      ci_high     = wt$conf.int[2],
      
      p_value     = wt$p.value
    )
  }) |>
  ungroup()

wilcox_results

affected_stats <- wilcox_results |>
  filter(hemisphere == "Affected")

unaffected_stats <- wilcox_results |>
  filter(hemisphere == "Unaffected")

# GRAPHICS

# proof of concept
fig2_affected <- avg_graph_pivot |> filter(hemisphere=="Affected") |>  
  ggplot(aes(x = period, y= avg_10min_rSO2, group=patient_id)) +
                                    geom_point() + geom_line(aes(group=patient_id))
fig2_affected

# figure for paper
# Figure 2, 2 panels, affected and unaffected hemisphere
# ----
# ------------------------------------------------------------------
# Data prep
# ------------------------------------------------------------------
plot_df <- avg_graph_pivot %>%
  mutate(
    period = factor(period, levels = c("start", "post")),
    hemisphere = factor(hemisphere, levels = c("Affected", "Unaffected"))
  ) %>%
  # keep only complete start/post pairs within each patient + hemisphere
  group_by(patient_id, hemisphere) %>%
  filter(all(c("start", "post") %in% period)) %>%
  ungroup()

# ------------------------------------------------------------------
# Paired p values per hemisphere
# Use Wilcoxon signed-rank test (robust choice for paired clinical data).
# If you prefer paired t-test, replace wilcox.test with t.test.
# ------------------------------------------------------------------
pval_df <- plot_df %>%
  dplyr::select(patient_id, hemisphere, period, avg_10min_rSO2) %>%
  pivot_wider(names_from = period, values_from = avg_10min_rSO2) %>%
  group_by(hemisphere) %>%
  summarise(
    p_value = wilcox.test(start, post, paired = TRUE)$p.value,
    y_max_data = max(c(start, post), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    p_label = paste0("p = ", format.pval(p_value, digits = 3, eps = 0.001)),
    y_bracket = y_max_data + 2.0,
    y_text = y_max_data + 3.0
  )

library(dplyr)
library(ggplot2)
library(patchwork)

# make sure period order is correct
plot_df <- plot_df %>%
  mutate(
    period = factor(period, levels = c("start", "post"))
  )

# ------------------------------------------------------------
# Panel A: affected hemisphere
# ------------------------------------------------------------
pA <- plot_df %>%
  filter(hemisphere == "Affected") %>%
  ggplot(
    aes(
      x = period,
      y = avg_10min_rSO2,
      group = patient_id
    )
  ) +
  
  geom_line(
    alpha = 0.10,
    linewidth = 0.35,
    colour = "grey40"
  ) +
  
  geom_point(
    alpha = 0.18,
    size = 1.2,
    colour = "grey30"
  ) +
  
  stat_summary(
    aes(group = 1),
    fun = mean,
    geom = "line",
    linewidth = 1.1
  ) +
  
  stat_summary(
    aes(group = 1),
    fun = mean,
    geom = "point",
    size = 3
  ) +
  
  stat_summary(
    aes(group = 1),
    fun.data = mean_cl_normal,
    geom = "errorbar",
    width = 0.08,
    linewidth = 0.8
  ) +
  
  labs(
    tag = "A",
    title = "Affected hemisphere",
    x = NULL,
    y = expression("Mean rSO"[2]~"(%)")
  ) +
  
  scale_x_discrete(
    labels = c(start = "Start", post = "Post")
  ) +
  
  geom_text(
    data = pval_df[1,],
    aes(x = 1.5, y = y_text, label = p_label),
    inherit.aes = FALSE,
    size = 4
  ) +
  
  theme_classic(base_size = 12) +
  
  theme(
    plot.tag = element_text(
      face = "bold",
      size = 15
    ),
    plot.tag.position = c(0, 1),
    plot.title = element_text(
      hjust = 0.5,
      size = 12
    ),
    panel.grid = element_blank()
  )

pA
# ------------------------------------------------------------
# Panel B: unaffected hemisphere
# ------------------------------------------------------------
pB <- plot_df %>%
  filter(hemisphere == "Unaffected") %>%
  ggplot(
    aes(
      x = period,
      y = avg_10min_rSO2,
      group = patient_id
    )
  ) +
  
  geom_line(
    alpha = 0.10,
    linewidth = 0.35,
    colour = "grey40"
  ) +
  
  geom_point(
    alpha = 0.18,
    size = 1.2,
    colour = "grey30"
  ) +
  
  stat_summary(
    aes(group = 1),
    fun = mean,
    geom = "line",
    linewidth = 1.1
  ) +
  
  stat_summary(
    aes(group = 1),
    fun = mean,
    geom = "point",
    size = 3
  ) +
  
  stat_summary(
    aes(group = 1),
    fun.data = mean_cl_normal,
    geom = "errorbar",
    width = 0.08,
    linewidth = 0.8
  ) +
  
  labs(
    tag = "B",
    title = "Unaffected hemisphere",
    x = NULL,
    y = NULL
  ) +

  geom_text(
    data = pval_df[2,],
    aes(x = 1.5, y = y_text, label = p_label),
    inherit.aes = FALSE,
    size = 4
  ) +
  
  scale_x_discrete(
    labels = c(start = "Start", post = "Post")
  ) +
  
  theme_classic(base_size = 12) +
  
  theme(
    plot.tag = element_text(
      face = "bold",
      size = 15
    ),
    plot.tag.position = c(0, 1),
    plot.title = element_text(
      hjust = 0.5,
      size = 12
    ),
    panel.grid = element_blank()
  )


# ------------------------------------------------------------
# Combine
# ------------------------------------------------------------
fig_prepost <- pA + pB +
  plot_layout(ncol = 2)

fig_prepost
# ----
# END FIGURE pre-post


# DIFERENCE in DIFFERENCE (Reviewer request 3)
change_between_hemispheres <- avg_graph_pivot |>
  dplyr::select(
    patient_id,
    hemisphere,
    period,
    avg_10min_rSO2
  ) |>
  
  # first make start/post columns
  pivot_wider(
    names_from = period,
    values_from = avg_10min_rSO2
  ) |>
  
  # only complete pre/post pairs
  filter(
    !is.na(start),
    !is.na(post)
  ) |>
  
  # change within each hemisphere
  mutate(
    delta_rSO2 = post - start
  ) |>
  
  dplyr::select(
    patient_id,
    hemisphere,
    delta_rSO2
  ) |>
  
  # then make affected/unaffected columns
  pivot_wider(
    names_from = hemisphere,
    values_from = delta_rSO2
  ) |>
  
  # only patients having BOTH hemispheres
  filter(
    !is.na(Affected),
    !is.na(Unaffected)
  ) |>
  
  mutate(
    delta_delta = Affected - Unaffected
  )

change_between_hemispheres
# Between-hemisphere change comparison
n_change <- nrow(change_between_hemispheres)

aff_change_median <- median(
  change_between_hemispheres$Affected,
  na.rm = TRUE
)

aff_change_q1 <- quantile(
  change_between_hemispheres$Affected,
  0.25,
  na.rm = TRUE
)

aff_change_q3 <- quantile(
  change_between_hemispheres$Affected,
  0.75,
  na.rm = TRUE
)
unaff_change_median <- median(
  change_between_hemispheres$Unaffected,
  na.rm = TRUE
)

unaff_change_q1 <- quantile(
  change_between_hemispheres$Unaffected,
  0.25,
  na.rm = TRUE
)

unaff_change_q3 <- quantile(
  change_between_hemispheres$Unaffected,
  0.75,
  na.rm = TRUE
)

# paired test: change in affected vs change in unaffected
between_hemi_test <- wilcox.test(
  change_between_hemispheres$Affected,
  change_between_hemispheres$Unaffected,
  paired = TRUE,
  conf.int = TRUE,
  conf.level = 0.95,
  exact = FALSE
)

HL_between <- unname(between_hemi_test$estimate)

CI_between_low <- between_hemi_test$conf.int[1]
CI_between_high <- between_hemi_test$conf.int[2]

p_between <- between_hemi_test$p.value


p_label <- ifelse(
  p_between < 0.001,
  "p < 0.001",
  paste0("p = ", sprintf("%.3f", p_between))
)

HL_estimate <- unname(between_hemi_test$estimate)
CI_low  <- between_hemi_test$conf.int[1]
CI_high <- between_hemi_test$conf.int[2]

# long format
change_long <- change_between_hemispheres |>
  dplyr::select(patient_id, Affected, Unaffected) |>
  tidyr::pivot_longer(
    cols = c(Affected, Unaffected),
    names_to = "hemisphere",
    values_to = "delta_rSO2"
  ) |>
  dplyr::mutate(
    hemisphere = factor(
      hemisphere,
      levels = c("Affected", "Unaffected")
    )
  )

y_max <- max(change_long$delta_rSO2, na.rm = TRUE)
y_bracket <- y_max + 2
y_text <- y_max + 3.5

fig_between_hemi <- ggplot(
  change_long,
  aes(
    x = hemisphere,
    y = delta_rSO2,
    group = patient_id
  )
) +
  
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.5
  ) +
  
  # individual paired data - subtle
  geom_line(
    alpha = 0.08,
    linewidth = 0.25,
    colour = "grey50"
  ) +
  
  geom_point(
    alpha = 0.30,
    size = 1.5,
    colour = "grey40"
  ) +
  
  # median summary
  stat_summary(
    aes(group = 1),
    fun = median,
    geom = "line",
    linewidth = 1.1
  ) +
  
  stat_summary(
    aes(group = 1),
    fun = median,
    geom = "point",
    size = 3.5
  ) +
  
  # significance bracket
  annotate(
    "segment",
    x = 1, xend = 2,
    y = y_bracket, yend = y_bracket,
    linewidth = 0.6
  ) +
  
  annotate(
    "segment",
    x = 1, xend = 1,
    y = y_bracket - 0.8, yend = y_bracket,
    linewidth = 0.6
  ) +
  
  annotate(
    "segment",
    x = 2, xend = 2,
    y = y_bracket - 0.8, yend = y_bracket,
    linewidth = 0.6
  ) +
  
  annotate(
    "text",
    x = 1.5,
    y = y_text,
    label = p_label,
    size = 4
  ) +
  
  labs(
    x = NULL,
    y = expression(Delta*rSO[2]~"(post - pre, percentage points)")
  ) +
  
  coord_cartesian(
    ylim = c(
      min(change_long$delta_rSO2, na.rm = TRUE) - 2,
      y_text + 2
    )
  ) +
  
  theme_classic(base_size = 12) +
  theme(
    panel.grid = element_blank()
  )

fig_between_hemi


# REVIEWER Request 5
# ----
# ============================================================
# Reviewer analysis:
# Does the patient-level NIRS effect differ by
# (1) collateral status and
# (2) meaningful reperfusion status?
#
# NIRS effect = affected - unaffected mean rSO2 at procedure start
# ============================================================


# ------------------------------------------------------------
# 1. Prepare start-of-procedure NIRS effect
# ------------------------------------------------------------

nirs_start_effect <- avg_start |>
  ungroup() |>
  transmute(
    patient_id,
    delta_rSO2_start = difference
  )


# ============================================================
# A. COLLATERAL STATUS
# ============================================================

collateral_df <- nirs_start_effect |>
  inner_join(
    db_clin |>
      distinct(patient_id, grp_coll_status),
    by = "patient_id"
  ) |>
  filter(
    grp_coll_status %in% c("bad", "good"),
    !is.na(delta_rSO2_start)
  ) |>
  mutate(
    grp_coll_status = factor(
      grp_coll_status,
      levels = c("bad", "good")
    )
  )

# Check sample size
table(collateral_df$grp_coll_status)
nrow(collateral_df)

# Descriptive statistics
collateral_summary <- collateral_df |>
  group_by(grp_coll_status) |>
  summarise(
    n = n(),
    median = median(delta_rSO2_start, na.rm = TRUE),
    q1 = quantile(delta_rSO2_start, 0.25, na.rm = TRUE),
    q3 = quantile(delta_rSO2_start, 0.75, na.rm = TRUE),
    .groups = "drop"
  )


# Direct group comparison
collateral_test <- wilcox.test(
  delta_rSO2_start ~ grp_coll_status,
  data = collateral_df,
  conf.int = TRUE,
  conf.level = 0.95,
  exact = FALSE
)

collateral_HL <- unname(collateral_test$estimate)
collateral_CI_low <- collateral_test$conf.int[1]
collateral_CI_high <- collateral_test$conf.int[2]
collateral_p <- collateral_test$p.value

collateral_bad <- collateral_summary |>
  filter(grp_coll_status == "bad")

collateral_good <- collateral_summary |>
  filter(grp_coll_status == "good")




# ============================================================
# FORMATTED P-VALUES
# ============================================================

collateral_p_txt <- ifelse(
  collateral_p < 0.001,
  "< 0.001",
  paste0("= ", sprintf("%.3f", collateral_p))
)


#----
# request 4 end

# RENDER
# ----
# render our statistics
rmarkdown::render(
  "c:/sci/rcode/nirs/nirs_rev1.Rmd",
  output_file = "affected_hemisphere_results.html",
  envir = globalenv()
)

# ----
# END RENDER







tb2 <- table(avg_post$ratio_signif)
tb2


db_row <- db_msr %>% filter(db_msr$patient_id == avg_2$patient_id[1]) 
db_row

avg_3

  



table(avg_start$ratio_signif)

sum(is.na(avg_start$Unaffected))










# old code 16.09.2026 - cleaning one patient
# ad row id
db_msr <- read.csv2(paste0(db.path, "nirs_measurements_clean_v2.csv"))

db_msr <- db_msr %>% mutate(row_id = row_number())

# all patients, only angiography_start_point == TRUE
length(unique(db_msr$patient_id))

id_true_first_try <- db_msr %>% group_by(patient_id) %>% 
  filter(first_try_point == TRUE) %>% 
  dplyr::select(patient_id)

no_first_try <- (setdiff(as_vector(db_clin$patient_id), as_vector(id_true_first_try)))
# which patients do not have angiography_start_point == TRUE
# 29 patients do not have this flag!
no_first_try

db_msr %>% group_by(patient_id) %>% 
  filter(recanalisation_point == TRUE) %>% summarise(n= n())

# 2 patients do not have flag for recanalisation_point
id_true <- db_msr %>% group_by(patient_id) %>% 
  filter(recanalisation_point == TRUE)

id_norecan <- setdiff(as_vector(db_clin$patient_id), as_vector(id_true))
id_norecan

# do 2 patients have 
db_msr %>% filter(patient_id %in% id_norecan & first_try_point == TRUE) %>% 
  dplyr::select(patient_id, datetime, first_try_point, recanalisation_point)

# the patient 1001590120 has wrong year (2025 instead of 2015)
# correct recan time is: 2015-10-22 12:15:00
row_tochange <- db_msr %>% filter(patient_id == 1001590120 &
                    datetime >= "2015-10-22 12:15:00"
                  & datetime <= "2015-10-22 12:16:00") %>% 
  slice_head(n=1) %>% dplyr::select(row_id)
as.numeric(row_tochange)

db_msr[db_msr$row_id == as.numeric(row_tochange),]$recanalisation_point <- TRUE

# the patient 1001249114 has wrong date order
row_tochange <- db_msr %>% filter(patient_id == 1001249114 &
                                    datetime >= "2013-07-01 09:15:00"
                                  & datetime <= "2013-07-01 09:16:00") %>% 
  slice_head(n=1) %>% dplyr::select(row_id)
as.numeric(row_tochange)

db_msr[db_msr$row_id == as.numeric(row_tochange),]$recanalisation_point <- TRUE

no_open <- db_clin %>% filter(opening_datetime == "") %>% dplyr::select(patient_id)
no_open <- as_vector(no_open$patient_id)
no_open
# change 13 patients to no opening time!
db_msr <- db_msr %>% mutate(
  recanalisation_point = case_when(patient_id %in% no_open ~ FALSE, TRUE ~ recanalisation_point)
)


# save new database
write.csv2(db_msr, paste0(db.path, "nirs_measurements_clean_v2.csv"))
