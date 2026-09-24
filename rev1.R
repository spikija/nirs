# analysis Revision 1
# LIBRARIES ----
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
# LIBRARIES ENDE----


# call usual functions
db.source.file.path <- ifelse(Sys.info()[[4]] == "LAPTOP-NFRMGQDO" | Sys.info()[[4]] == "NEURONODE", "c:/sci/cdk.stat/", "C:/app/r.apps/cdk.stat/")
source(paste0(db.source.file.path, "sp.utilities.r"))
source(paste0(db.source.file.path, "nirs_helper.r"))

db.language <- ifelse(Sys.info()[[4]] == "LAPTOP-NFRMGQDO" | Sys.info()[[4]] == "NEURONODE", "Meine Ablage", "My Drive")
db.path <- paste0("C:/sci/rcode/nirs/")

# DBDBDBDBDBDBDBDBDBDBDBDBDBDBDBDBDBDBDBDBDBDBDBDBDBDBDBDBDB
# LAST STABLE DB
db_clin <- read.csv2(paste0(db.path, "clinical_nirs_24092026_3.csv"))

# load main database
db_msr <- read.csv2(paste0(db.path, "nirs_24092026_1.csv"))

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

# QUESTION 1----
avg_start <- db_msr %>% group_by(patient_id, hemisphere_status) %>% 
  filter(elapsed_seconds >= 60*init_delay & elapsed_seconds <= 60*(init_delay+sample_duration)) %>% 
  #summarise(avg = mean(nirs_value_raw, na.rm = TRUE)) %>% 
  summarise(avg = median(nirs_value_raw, na.rm = TRUE)) %>% 
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

# median across entire brain at the baseline (preparing for multivariable analysis)
avg_overall <- db_msr %>% group_by(patient_id) %>% 
  filter(elapsed_seconds >= 60*init_delay & elapsed_seconds <= 60*(init_delay+sample_duration)) %>% 
  summarise(rso2_med = median(nirs_value_raw, na.rm = TRUE)) 
avg_overall <- inner_join(avg_overall, avg_start |> dplyr::select(patient_id, difference, ratio_signif), 
by="patient_id") 
avg_overall$difference <- abs(avg_overall$difference)

# Affected-hemisphere identification----

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


# Exact 95% CI for each classification in the full cohort

ci_correct <- binom.test(
  x = n_unaff_higher,
  n = n_total
)

ci_incorrect <- binom.test(
  x = n_aff_higher,
  n = n_total
)

ci_indeterminate <- binom.test(
  x = n_equal,
  n = n_total
)

correct_ci_low  <- 100 * ci_correct$conf.int[1]
correct_ci_high <- 100 * ci_correct$conf.int[2]

incorrect_ci_low  <- 100 * ci_incorrect$conf.int[1]
incorrect_ci_high <- 100 * ci_incorrect$conf.int[2]

indeterminate_ci_low  <- 100 * ci_indeterminate$conf.int[1]
indeterminate_ci_high <- 100 * ci_indeterminate$conf.int[2]


# Among patients with directional asymmetry:
# lower rSO2 correctly identifies affected hemisphere

n_determinate <- n_unaff_higher + n_aff_higher

binom_test <- binom.test(
  x = n_unaff_higher,
  n = n_determinate,
  p = 0.5
)

ident_accuracy <- 100 * unname(binom_test$estimate)
ident_ci_low   <- 100 * binom_test$conf.int[1]
ident_ci_high  <- 100 * binom_test$conf.int[2]

binom_p <- binom_test$p.value

# so with these tests we have answered at least one question:

# the frequency of affected hemisphere having lower rSO2 is significant
# p = 0.041 (95% CI 0.30 - 0.49)

# producing dot plot
# graph will show differences between hemispheres in the 5 minute window 10 minutes
# after procedure start

# proof of concept
fig1 <- avg_start %>% 
  ggplot(aes(x = 0, y = difference)) + 
  geom_dotplot(binwidth = 0.5, dotsize = 1.5, binaxis = "y", stackdir = "center") 
fig1 

# production graph ----
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
# Block for Q1


# Q2 : 

# include only 2b-3 recanalisations

avg_2 <- db_msr |>
  semi_join(
    db_clin |>
      filter(tici_grp == "2b-3") |>
      dplyr::select(patient_id),
    by = "patient_id"
  ) |> 
  group_by(patient_id, elapsed_seconds) %>% 
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
     #summarise(avg = mean(nirs_value_raw, na.rm = TRUE)) %>% 
    summarise(avg = median(nirs_value_raw, na.rm = TRUE)) %>% 
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

# RESULTS----

# *************************
# Paired start vs post analysis
# *************************

# *************************
# Prepare paired start/post data
# *************************

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

# *************************
# Wilcoxon signed-rank test + Hodges-Lehmann estimate
# *************************

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
# Figure 2, 2 panels, affected and unaffected hemisphere----
# *******************************
# Data prep
# *******************************
plot_df <- avg_graph_pivot %>%
  mutate(
    period = factor(period, levels = c("start", "post")),
    hemisphere = factor(hemisphere, levels = c("Affected", "Unaffected"))
  ) %>%
  # keep only complete start/post pairs within each patient + hemisphere
  group_by(patient_id, hemisphere) %>%
  filter(all(c("start", "post") %in% period)) %>%
  ungroup()

# *******************************
# Paired p values per hemisphere
# Use Wilcoxon signed-rank test (robust choice for paired clinical data).
# If you prefer paired t-test, replace wilcox.test with t.test.
# *******************************
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

# *************************
# Panel A: affected hemisphere
# *************************
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
# *************************
# Panel B: unaffected hemisphere
# *************************
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


# *************************
# Combine
# *************************
fig_prepost <- pA + pB +
  plot_layout(ncol = 2)

fig_prepost
# END FIGURE pre-post ----


# DIFERENCE in DIFFERENCE (Reviewer request 3)----
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


# REVIEWER Request 5----

# **********************************
# Reviewer analysis:
# Does the patient-level NIRS effect differ by
# (1) collateral status and
# (2) meaningful reperfusion status?
#
# NIRS effect = affected - unaffected mean rSO2 at procedure start
# **********************************


# *************************
# 1. Prepare start-of-procedure NIRS effect
# *************************

nirs_start_effect <- avg_start |>
  ungroup() |>
  transmute(
    patient_id,
    delta_rSO2_start = difference
  )


# **********************************
# A. COLLATERAL STATUS
# **********************************

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




# **********************************
# FORMATTED P-VALUES
# **********************************

collateral_p_txt <- ifelse(
  collateral_p < 0.001,
  "< 0.001",
  paste0("= ", sprintf("%.3f", collateral_p))
)


# request 4 end


# MULTIARIABLE ANALYSIS ----
# join db_clin and db_overall

db_stat <- inner_join(db_clin |> dplyr::select(patient_id, nihss_bei_aufnahme_nur_summe, 
                                               alter,tici_grp, mRS_3mo_grp), avg_overall, by="patient_id")
# this is for change in hemispheres
db_stat_delta_delta <- inner_join(db_stat, change_between_hemispheres[,c("patient_id", "delta_delta")], by="patient_id")

df <- db_stat %>% 
  mutate(
    # Outcome: mRS 0-2 vs 3-6
    outcome = factor(mRS_3mo_grp, levels = c("0–2", "3–6")),
    favorable_outcome = as.integer(outcome=="0–2"),
    
    # Covariates
    nihss_adm = nihss_bei_aufnahme_nur_summe, #case_when(nihss_cat == "0–6" | nihss_cat == "7–15" ~ "0-15", TRUE ~ "16+"),
    tici_grp  = factor(tici_grp),               # categorical
    age = alter
    
  )

df_delta_delta <- 
  db_stat_delta_delta %>% 
  mutate(
    # Outcome: mRS 0-2 vs 3-6
    outcome = factor(mRS_3mo_grp, levels = c("0–2", "3–6")),
    favorable_outcome = as.integer(outcome=="0–2"),
    
    # Covariates
    nihss_adm = nihss_bei_aufnahme_nur_summe, #case_when(nihss_cat == "0–6" | nihss_cat == "7–15" ~ "0-15", TRUE ~ "16+"),
    tici_grp  = factor(tici_grp),               # categorical
    age = alter
    
  )

# Helper: function to fit model + produce OR table
fit_logit_report <- function(data, formula) {
  m <- glm(formula, data = data, family = binomial())
  
  # Odds ratios with 95% CI
  or_tbl <- broom::tidy(m, conf.int = TRUE, exponentiate = TRUE) %>%
    mutate(across(where(is.numeric), ~ round(.x, 3)))
  
  list(
    model = m,
    OR_table = or_tbl,
    vif = tryCatch(car::vif(m), error = function(e) NA),
    pseudoR2 = tryCatch(pscl::pR2(m), error = function(e) NA)
  )
}

# Refit Model A1 mRS 0–2 vs 3–6∼age+NIHSS+mTICI+overall rSO₂ level+∣affected−unaffected rSO₂∣
df_A1 <- df %>% dplyr::select(favorable_outcome, age, nihss_adm, tici_grp, rso2_med, difference) %>% na.omit()
res_A1_simpl <- fit_logit_report(df_A1, favorable_outcome ~ rso2_med + difference + age + nihss_adm + tici_grp)

# Refit Model A2 instead of asymmetry use direction ratio_signif
df_A2 <- df %>% dplyr::select(favorable_outcome, age, nihss_adm, tici_grp, rso2_med, ratio_signif) %>% na.omit()
res_A2_simpl <- fit_logit_report(df_A2, favorable_outcome ~ rso2_med + ratio_signif + age + nihss_adm + tici_grp)

# Refit Model A3 delta-delta on 151 patients (only tici grp 3)
df_A3 <- df_delta_delta %>% dplyr::select(favorable_outcome, age, nihss_adm, rso2_med, delta_delta) %>% na.omit()
res_A3_simpl <- fit_logit_report(df_A3, favorable_outcome ~ age + nihss_adm + rso2_med + delta_delta)


res_A1_simpl$OR_table
res_A1_simpl$vif
res_A1_simpl$pseudoR2

res_A2_simpl$OR_table
res_A2_simpl$vif
res_A2_simpl$pseudoR2

res_A3_simpl$OR_table
res_A3_simpl$vif
res_A3_simpl$pseudoR2

format_or_table <- function(tbl, model_label) {
  tbl %>%
    transmute(
      Model = model_label,
      Variable = term,
      `OR (95% CI)` = sprintf(
        "%.2f (%.2f–%.2f)",
        estimate, conf.low, conf.high
      ),
      `p value` = ifelse(p.value < 0.001, "<0.001", sprintf("%.3f", p.value))
    )
}

tbl_A1 <- format_or_table(res_A1_simpl$OR_table, "Baseline overall oxygenation and abs interhemispheric difference")
tbl_A2 <- format_or_table(res_A2_simpl$OR_table, "Baseline interhemispheric difference direction")
tbl_A3 <- format_or_table(res_A3_simpl$OR_table, "Intrahemispheric difference in recanalised population (start-post)")

tbl_all <- bind_rows(tbl_A1, tbl_A2, tbl_A3)

label_map <- c(
  "(Intercept)" = "Intercept",
  
  # Common covariates
  "age" = "Age, per year",
  "nihss_adm" = "Admission NIHSS, per point",
  "rso2_med" = "Overall rSO₂, per 1%-point",
  "tici_grp2b-3" = "Successful reperfusion (TICI 2b–3)",
  "ratio_signif" = "Direction of interhemispheric difference",
  
  # Model A1
  "difference" = "Absolute interhemispheric rSO₂ difference, per 1%-point",
  
  # Model A2
  "ratio_signifAffected higher" = "Affected hemisphere higher",
  "ratio_signifUnaffected higher" = "Unaffected hemisphere higher",
  
  # Model A3
  "delta_delta" = "Change in interhemispheric rSO₂ difference, per 1%-point"
)

tbl_all <- tbl_all %>%
  mutate(
    Variable = dplyr::recode(Variable, !!!label_map)
  )

ft_multivariable <- flextable(tbl_all) %>%
  set_header_labels(
    Model = "",
    Variable = "Variable",
    `OR (95% CI)` = "Odds ratio (95% CI)",
    `p value` = "p value"
  ) %>%
  merge_v(j = "Model") %>%
  valign(j = "Model", valign = "top") %>%
  autofit() %>%
  align(j = c("OR (95% CI)", "p value"), align = "center") %>%
  bold(part = "header") %>%
  fontsize(size = 9, part = "all") %>%
  set_caption(
    caption = "Supplemental Table X. Multivariable logistic regression models evaluating associations of baseline NIRS-derived rSO₂ measures and peri-procedural changes in rSO₂ with 3-month functional outcome (mRS 0–2 vs 3–6)."
  )

library(officer)

doc <- read_docx() %>%
  body_add_flextable(ft)

print(doc, target = "Supplemental_Table_Multivariable_NIRS.docx")



# Multivariable analysis END----





# RENDER----

# render our statistics
rmarkdown::render(
  "c:/sci/rcode/nirs/nirs_rev1.Rmd",
  output_file = "affected_hemisphere_results.html",
  envir = globalenv()
)

# END RENDER

# FLOW CHART GRAPH----
# SECTION 1
# patient flow chart
# prepare data for rmarkdown
flowChartData <- list()

#orig_db <- read.csv2(paste0(db.path, "db_total_cleaned_11022026.csv2"))


orig_db <- read.delim(
  paste0(db.path, "db_total_cleaned_11022026.csv2"),
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE
) |> 
  as_tibble()

flowChartData$total_patients <- nrow(orig_db) + 30 # 30 due to newly added patients after Zensus with ArchiMed (04.02.2026)
#flowChartData$non_anteriors <- as.numeric(orig_db %>% filter(side == "other") %>% summarise(n = n()))
flowChartData$non_anteriors <- as.numeric(orig_db %>% filter(!vessel_type_angio %in% c("ACI Occlusion", "M1", "M2")) %>% summarise(n = n()))
db_fc <- read_excel(paste0(db.path, "patient_flow_chart_03012026.xlsx"))
# db <- right_join(db, db_fc[db_fc$status == "merge",], by="patient_id")
flowChartData$nirs_data <- nrow(db_fc)


f_rows <- paste0( "<TR>", "<TD ALIGN='LEFT'>", esc(names(other_tab)), "</TD>", "<TD ALIGN='RIGHT'><B>", as.integer(other_tab), "</B></TD>", "</TR>", collapse = "" )

# --- exclusions table label; set a LARGE font inside the HTML TABLE ---
f_label <- paste0(
  "<TABLE BORDER='0' CELLBORDER='0' CELLPADDING='8' BGCOLOR='White' FONTSIZE='28'>",
  "<TR><TD ALIGN='LEFT' COLSPAN='2'><B>Exclusions: ", excl_total, "</B></TD></TR>",
  f_rows,
  "</TABLE>"
)

library(glue)
graph_code <- glue("
digraph flowchart {{

  // 1) Orthogonal edges + bigger base fonts + more spacing
  graph [layout=dot, rankdir=TB, fontname=\"Calibri\", fontsize=28,
         splines=ortho, nodesep=0.6, ranksep=0.7];
  node  [shape=rectangle, style=filled, fillcolor=White,
         fontname=\"Calibri\", fontsize=28, margin=\"0.25,0.15\"];
  edge  [fontname=\"Calibri\", fontsize=28, arrowsize=1.4];

  // existing nodes
  A [label=<
      Total Patients<BR/><B>{flowChartData$total_patients}</B><BR/><FONT POINT-SIZE=\"18\">April 2010 to December 2023</FONT>
  >];

  B [label=<
      Anterior occlusive stroke<BR/><B>{flowChartData$total_patients - flowChartData$non_anteriors}</B>
  >];

  C [label=<
      Other vessels<BR/><B>{flowChartData$non_anteriors}</B>
  >];

  C2 [label=<
      No NIRS measurements<BR/><B>{(flowChartData$total_patients - flowChartData$non_anteriors)-flowChartData$nirs_data}</B>
  >];
  
  D [label=<
      With NIRS data<BR/><B>{flowChartData$nirs_data}</B>
  >];

  // spacer/connector under D
  E [shape=point, label=\"\", width=0];

  // exclusions (HTML table already has FONTSIZE=28 inside)
  F [shape=box, label=<{f_label}>, margin=\"0.18,0.18\"];

  // merged
  G [label=<Final cases to analyze<BR/><B>{merged_n}</B>>];

  // ------- main stem -------
  A -> blank1 [dir=none];
  blank1 [shape=point, label=\"\", width=0];
  blank1 -> C [minlen = 3];
  {{ rank = same; blank1 C }}
  
  blank1 -> B;
  
  B -> blank2 [dir=none];
  
  blank2 [shape=point, label=\"\", width=0];
  blank2 -> C2 [minlen = 3];
  {{ rank = same; blank2 C2 }}
  
  blank2 -> D;
  
  // split to exclusions + merged
  D -> E [dir=none];
  {{ rank = same; E F }}
  E -> F [constraint=false, minlen=3];
  E -> G

}}
")

DiagrammeR::grViz(graph_code, width = 1100, height = 900)

library(DiagrammeR)
library(DiagrammeRsvg)
library(rsvg)

#svg <- export_svg(grViz(DiagrammeR::grViz(graph_code, width = 1100, height = 900)))
svg <- DiagrammeRsvg::export_svg(
  DiagrammeR::grViz(graph_code)
)
# DiagrammeR::grViz(graph_code)
# rsvg_png(
#   charToRaw(svg),
#   file = paste0(db.path,"Figure_1_flowchart.png"),
#   width = 1100,
#   height = 900
# )

library(DiagrammeR)
library(DiagrammeRsvg)
library(rsvg)
library(magick)

# Create graph
g <- DiagrammeR::grViz(graph_code)

# Export to SVG
svg <- DiagrammeRsvg::export_svg(g)

# Desired final size
dpi <- 300
width_in  <- 6
height_in <- 5

width_px  <- width_in  * dpi
height_px <- height_in * dpi

# Save PNG
png_file <- paste0("G:/Meine Ablage/sci/graz/retrospective/NIRS EVT 2025/Manuscript/Submissions/jcbfm/", "Figure_1_flowchart.png")

rsvg::rsvg_png(
  charToRaw(svg),
  file = png_file,
  width = 1800,
  height = 1500
)

# Optional but useful: write 300 dpi metadata into the PNG
magick::image_write(img, path = png_file, format = "png")



# DEMOGRAPHIC TABLE----
# TABLE DEMOGRAPHICS MAIN TEXT----
# Variables to describe
library(dplyr)
library(gtsummary)
library(flextable)

# ============================================================
# Demographic / clinical Table 1
# ============================================================

var_desc <- c(
  "alter", "geschlecht", "pre_mRS", "stroke_typ", "side", "zeitpunkt",
  "hypertonus2", "hypercholesterinaemie", "diabetes",
  "kardiovaskulaere_erkrankungen", "andere_risikofaktoren",
  "periphere_avk", "nikotinabusus", "vorhofflimmern2",
  "vorbehandlung_mit_thrombozytenfh",
  "vorbehandlung_orale_antikoagulanzien",
  "thrombozytenzahl_bei_aufnahme",
  "serumglucose_bei_aufnahme", "crp", "haematokrit",
  "iv_thrombolyse_mit_rtpa", "a_thrombolyse",
  "anzahl_verwendeter_devices", "extrakranieller_stent",
  "vessel_type_preangio",
  "ipsilat_extrakran_hochgr_aci_stenose_verschluss",
  "tici_grp", "nihss_bei_aufnahme_nur_summe",
  "initiale_bildgebung", "infarktfruehzeichen",
  "periprozedurale_dissektion_nicht_flussrel",
  "periprozedurale_dissektion_flussrel",
  "intrazerebrale_blutung", "sab",
  "thrombembolische_komplikationen",
  "reokklusion_des_zielgefaesses",
  "komplikationen_an_punktionsstelle",
  "ecass_code", "hemikraniektomie",
  "final_infarct",
  "stenose_im_therapierten_gefaess",
  "schlaganfall_aetiologie_toast",
  "nihss_entl_max",
  "mt_duration",
  "time_diff_onset_thrombectomy_time",
  "time_diff_onset_lyse",
  "mRS_3mo_grp",
  "grp_coll_status"
)

var_labels <- c(
  alter = "Age (years)",
  geschlecht = "Sex",
  pre_mRS = "Premorbid mRS",
  stroke_typ = "Stroke onset type",
  side = "Affected hemisphere",
  zeitpunkt = "Time of admission",
  hypertonus2 = "Hypertension",
  hypercholesterinaemie = "Hypercholesterolemia",
  diabetes = "Diabetes mellitus",
  kardiovaskulaere_erkrankungen = "Cardiovascular disease",
  andere_risikofaktoren = "Other vascular risk factors",
  periphere_avk = "Peripheral arterial disease",
  nikotinabusus = "Smoking",
  vorhofflimmern2 = "Atrial fibrillation",
  vorbehandlung_mit_thrombozytenfh = "Antiplatelet treatment prior to stroke",
  vorbehandlung_orale_antikoagulanzien = "Oral anticoagulation prior to stroke",
  thrombozytenzahl_bei_aufnahme = "Platelet count (G/L)",
  serumglucose_bei_aufnahme = "Serum glucose (mg/dL)",
  crp = "CRP (mg/dL)",
  haematokrit = "Hematocrit (%)",
  iv_thrombolyse_mit_rtpa = "Intravenous thrombolysis with alteplase",
  a_thrombolyse = "Intra-arterial thrombolysis",
  anzahl_verwendeter_devices = "Number of devices used",
  extrakranieller_stent = "Extracranial stent implantation",
  vessel_type_preangio = "Vessel occlusion",
  ipsilat_extrakran_hochgr_aci_stenose_verschluss =
    "Ipsilateral high-grade ICA stenosis or occlusion",
  tici_grp = "Final mTICI",
  nihss_bei_aufnahme_nur_summe = "NIHSS at admission",
  initiale_bildgebung = "Initial brain imaging",
  infarktfruehzeichen = "Early ischemic changes",
  periprozedurale_dissektion_nicht_flussrel =
    "Dissection, not hemodynamically relevant",
  periprozedurale_dissektion_flussrel =
    "Dissection, hemodynamically relevant",
  intrazerebrale_blutung = "Intracerebral hemorrhage",
  sab = "Subarachnoid hemorrhage",
  thrombembolische_komplikationen = "Thromboembolic complications",
  reokklusion_des_zielgefaesses = "Reocclusion of target vessel",
  komplikationen_an_punktionsstelle = "Groin puncture complications",
  ecass_code = "ECASS hemorrhage classification",
  hemikraniektomie = "Hemicraniectomy",
  final_infarct = "Final infarct size",
  stenose_im_therapierten_gefaess = "Stenosis grade in target vessel",
  schlaganfall_aetiologie_toast = "Stroke etiology (TOAST)",
  nihss_entl_max = "NIHSS at discharge",
  mt_duration = "Intervention duration (min)",
  time_diff_onset_thrombectomy_time = "Onset to groin puncture (min)",
  time_diff_onset_lyse = "Onset to intravenous thrombolysis (min)",
  mRS_3mo_grp = "mRS at 3 months",
  grp_coll_status = "Collateral status"
)

# ------------------------------------------------------------
# Check which requested variables exist in db_clin
# ------------------------------------------------------------

vars_present <- intersect(var_desc, names(db_clin))
vars_missing <- setdiff(var_desc, names(db_clin))

cat(
  "Variables not found in db_clin:\n",
  paste(vars_missing, collapse = "\n"),
  "\n\n"
)

# ------------------------------------------------------------
# Variables to treat explicitly as continuous
# ------------------------------------------------------------

continuous_vars <- c(
  "alter",
  "nihss_bei_aufnahme_nur_summe",
  "nihss_entl_max",
  "thrombozytenzahl_bei_aufnahme",
  "serumglucose_bei_aufnahme",
  "crp",
  "haematokrit",
  "anzahl_verwendeter_devices",
  "mt_duration",
  "time_diff_onset_thrombectomy_time",
  "time_diff_onset_lyse"
)

continuous_vars <- intersect(
  continuous_vars,
  vars_present
)

# ------------------------------------------------------------
# Create Table 1
# ------------------------------------------------------------


tbl_demographic <- db_clin |>
  dplyr::select(
    dplyr::all_of(vars_present)
  ) |>
  gtsummary::tbl_summary(
    type = list(
      dplyr::all_of(continuous_vars) ~ "continuous"
    ),
    label = as.list(
      var_labels[vars_present]
    ),
    statistic = list(
      gtsummary::all_continuous() ~ "{median} ({p25}, {p75})",
      gtsummary::all_categorical() ~ "{n} ({p}%)"
    ),
    digits = list(
      gtsummary::all_continuous() ~ 1,
      gtsummary::all_categorical() ~ c(0, 1)
    ),
    missing = "ifany",
    missing_text = "Missing"
  ) |>
  gtsummary::modify_header(
    label ~ "**Characteristic**",
    stat_0 ~ "**N = {N}**"
  ) |>
  gtsummary::bold_labels()

# ------------------------------------------------------------
# Convert to flextable
# ------------------------------------------------------------

ft_demographic <- tbl_demographic |>
  gtsummary::as_flex_table() |>
  flextable::set_caption(
    caption = paste0(
      "Table 1. Demographic, clinical, imaging, and procedural ",
      "characteristics of the study population."
    )
  ) |>
  flextable::autofit()

ft_demographic




# TABLE 2----
# Patient-level interhemispheric rSO2 asymmetry

# Rename conceptually -- avoids "correct/incorrect" terminology
n_aff_lower <- n_unaff_higher
n_aff_higher_tbl <- n_aff_higher
n_no_asym <- n_equal

n_total <- n_aff_lower + n_aff_higher_tbl + n_no_asym
n_relevant_asym <- n_aff_lower + n_aff_higher_tbl

# Exact binomial 95% CI helper
get_binom_ci <- function(x, n) {
  ci <- stats::binom.test(x, n)$conf.int * 100
  c(low = ci[1], high = ci[2])
}

ci_aff_lower <- get_binom_ci(n_aff_lower, n_total)
ci_no_asym <- get_binom_ci(n_no_asym, n_total)
ci_aff_higher <- get_binom_ci(n_aff_higher_tbl, n_total)
ci_any_asym <- get_binom_ci(n_relevant_asym, n_total)

# Direction among patients with relevant asymmetry
ci_lower_among_asym <- get_binom_ci(
  n_aff_lower,
  n_relevant_asym
)

table2_df <- tibble::tibble(
  `Interhemispheric pattern` = c(
    "Affected hemisphere lower",
    "No relevant asymmetry",
    "Affected hemisphere higher",
    "Any relevant asymmetry"
  ),
  
  Definition = c(
    "Affected − unaffected < −4 pp",
    "|Affected − unaffected| ≤ 4 pp",
    "Affected − unaffected > +4 pp",
    "|Affected − unaffected| > 4 pp"
  ),
  
  `n/N` = c(
    sprintf("%d/%d", n_aff_lower, n_total),
    sprintf("%d/%d", n_no_asym, n_total),
    sprintf("%d/%d", n_aff_higher_tbl, n_total),
    sprintf("%d/%d", n_relevant_asym, n_total)
  ),
  
  `%` = c(
    100 * n_aff_lower / n_total,
    100 * n_no_asym / n_total,
    100 * n_aff_higher_tbl / n_total,
    100 * n_relevant_asym / n_total
  ),
  
  `95% CI` = c(
    sprintf("%.1f–%.1f", ci_aff_lower["low"], ci_aff_lower["high"]),
    sprintf("%.1f–%.1f", ci_no_asym["low"], ci_no_asym["high"]),
    sprintf("%.1f–%.1f", ci_aff_higher["low"], ci_aff_higher["high"]),
    sprintf("%.1f–%.1f", ci_any_asym["low"], ci_any_asym["high"])
  )
)

ft_table2 <- table2_df |>
  flextable::flextable() |>
  flextable::colformat_double(
    j = "%",
    digits = 1
  ) |>
  flextable::set_header_labels(
    `Interhemispheric pattern` = "Interhemispheric pattern",
    Definition = "Definition",
    `n/N` = "n/N",
    `%` = "%",
    `95% CI` = "95% CI"
  ) |>
  flextable::bold(i = 4, bold = TRUE) |>
  flextable::align(
    j = c("n/N", "%", "95% CI"),
    align = "center",
    part = "all"
  ) |>
  flextable::add_footer_lines(
    values = paste0(
      "Among patients with a relevant asymmetry (>4 percentage points), ",
      "rSO₂ was lower in the affected hemisphere in ",
      n_aff_lower, "/", n_relevant_asym,
      " (", sprintf("%.1f", 100 * n_aff_lower / n_relevant_asym), "%; ",
      "95% CI ", sprintf("%.1f", ci_lower_among_asym["low"]), "–",
      sprintf("%.1f", ci_lower_among_asym["high"]), "%). ",
      "Direction of asymmetry is reported descriptively and should not be ",
      "interpreted as correct versus incorrect lateralisation."
    )
  ) |>
  flextable::set_caption(
    caption = paste0(
      "Table 2. Patient-level interhemispheric rSO₂ asymmetry at ",
      "procedural baseline."
    )
  ) |>
  flextable::autofit()

ft_table2


# TABLE 3----
# Patient-level temporal rSO2 response

# ------------------------------------------------------------
# 1. Create one pre/post pair per patient and hemisphere
# ------------------------------------------------------------

prepost_wide <- plot_df |>
  dplyr::select(
    patient_id,
    hemisphere,
    period,
    avg_10min_rSO2
  ) |>
  tidyr::pivot_wider(
    names_from = period,
    values_from = avg_10min_rSO2
  ) |>
  dplyr::filter(
    !is.na(start),
    !is.na(post)
  ) |>
  dplyr::mutate(
    change = post - start
  )

# ------------------------------------------------------------
# 2. Hemisphere-specific descriptive statistics
# ------------------------------------------------------------

hemi_summary <- prepost_wide |>
  dplyr::group_by(hemisphere) |>
  dplyr::summarise(
    n = dplyr::n(),
    
    pre_median = median(start, na.rm = TRUE),
    pre_q1 = quantile(start, 0.25, na.rm = TRUE),
    pre_q3 = quantile(start, 0.75, na.rm = TRUE),
    
    post_median = median(post, na.rm = TRUE),
    post_q1 = quantile(post, 0.25, na.rm = TRUE),
    post_q3 = quantile(post, 0.75, na.rm = TRUE),
    
    change_median = median(change, na.rm = TRUE),
    change_q1 = quantile(change, 0.25, na.rm = TRUE),
    change_q3 = quantile(change, 0.75, na.rm = TRUE),
    
    .groups = "drop"
  )

# ------------------------------------------------------------
# 3. Within-hemisphere Wilcoxon estimates
#
# Testing CHANGE = post - pre ensures:
# positive estimate = increase after reperfusion
# ------------------------------------------------------------

get_wilcox_change <- function(x) {
  
  wt <- stats::wilcox.test(
    x,
    mu = 0,
    conf.int = TRUE,
    exact = FALSE
  )
  
  tibble::tibble(
    HL_estimate = unname(wt$estimate),
    ci_low = wt$conf.int[1],
    ci_high = wt$conf.int[2],
    p_value = wt$p.value
  )
}

hemi_tests <- prepost_wide |>
  dplyr::group_by(hemisphere) |>
  dplyr::group_modify(
    ~ get_wilcox_change(.x$change)
  ) |>
  dplyr::ungroup()

hemi_results <- hemi_summary |>
  dplyr::left_join(
    hemi_tests,
    by = "hemisphere"
  )

# ------------------------------------------------------------
# 4. Direct between-hemisphere comparison of change
#
# delta_delta =
# (post - pre) affected - (post - pre) unaffected
# ------------------------------------------------------------

delta_df <- change_between_hemispheres |>
  dplyr::select(
    patient_id,
    Affected,
    Unaffected,
    delta_delta
  ) |>
  dplyr::filter(
    !is.na(Affected),
    !is.na(Unaffected),
    !is.na(delta_delta)
  )

delta_test <- stats::wilcox.test(
  delta_df$delta_delta,
  mu = 0,
  conf.int = TRUE,
  exact = FALSE
)

delta_summary <- tibble::tibble(
  n = nrow(delta_df),
  change_median = median(
    delta_df$delta_delta,
    na.rm = TRUE
  ),
  change_q1 = quantile(
    delta_df$delta_delta,
    0.25,
    na.rm = TRUE
  ),
  change_q3 = quantile(
    delta_df$delta_delta,
    0.75,
    na.rm = TRUE
  ),
  HL_estimate = unname(delta_test$estimate),
  ci_low = delta_test$conf.int[1],
  ci_high = delta_test$conf.int[2],
  p_value = delta_test$p.value
)

# ------------------------------------------------------------
# 5. Construct publication table
# ------------------------------------------------------------

format_med_iqr <- function(med, q1, q3) {
  sprintf("%.1f [%.1f–%.1f]", med, q1, q3)
}

format_est_ci <- function(est, low, high) {
  sprintf("%.1f (%.1f to %.1f)", est, low, high)
}

format_p <- function(p) {
  ifelse(
    p < 0.001,
    "<0.001",
    sprintf("%.3f", p)
  )
}

table3_hemi <- hemi_results |>
  dplyr::mutate(
    Measure = as.character(hemisphere),
    
    `Pre-reperfusion rSO₂, %` =
      format_med_iqr(
        pre_median,
        pre_q1,
        pre_q3
      ),
    
    `Post-reperfusion rSO₂, %` =
      format_med_iqr(
        post_median,
        post_q1,
        post_q3
      ),
    
    `ΔrSO₂, pp` =
      format_med_iqr(
        change_median,
        change_q1,
        change_q3
      ),
    
    `HL estimate (95% CI), pp` =
      format_est_ci(
        HL_estimate,
        ci_low,
        ci_high
      ),
    
    `P value` =
      format_p(p_value)
  ) |>
  dplyr::select(
    Measure,
    n,
    `Pre-reperfusion rSO₂, %`,
    `Post-reperfusion rSO₂, %`,
    `ΔrSO₂, pp`,
    `HL estimate (95% CI), pp`,
    `P value`
  )

table3_delta <- tibble::tibble(
  Measure = "Differential change (ΔΔrSO₂)",
  n = delta_summary$n,
  `Pre-reperfusion rSO₂, %` = "—",
  `Post-reperfusion rSO₂, %` = "—",
  `ΔrSO₂, pp` =
    format_med_iqr(
      delta_summary$change_median,
      delta_summary$change_q1,
      delta_summary$change_q3
    ),
  `HL estimate (95% CI), pp` =
    format_est_ci(
      delta_summary$HL_estimate,
      delta_summary$ci_low,
      delta_summary$ci_high
    ),
  `P value` =
    format_p(delta_summary$p_value)
)

table3_df <- dplyr::bind_rows(
  table3_hemi,
  table3_delta
)

# ------------------------------------------------------------
# 6. Flextable
# ------------------------------------------------------------

ft_table3 <- table3_df |>
  flextable::flextable() |>
  flextable::set_header_labels(
    Measure = "Measure",
    n = "n",
    `Pre-reperfusion rSO₂, %` = "Pre-reperfusion rSO₂, %",
    `Post-reperfusion rSO₂, %` = "Post-reperfusion rSO₂, %",
    `ΔrSO₂, pp` = "Change ΔrSO₂, pp",
    `HL estimate (95% CI), pp` = "HL estimate (95% CI), pp",
    `P value` = "P value"
  ) |>
  flextable::bold(
    i = 3,
    bold = TRUE
  ) |>
  flextable::align(
    j = 2:7,
    align = "center",
    part = "all"
  ) |>
  flextable::add_footer_lines(
    values = paste0(
      "Values are median [IQR]. ΔrSO₂ = post-reperfusion minus ",
      "pre-reperfusion rSO₂. Differential change (ΔΔrSO₂) = ",
      "(post − pre)affected − (post − pre)unaffected. ",
      "Positive ΔΔrSO₂ values indicate a greater increase, or a smaller ",
      "decrease, in the affected hemisphere relative to the unaffected ",
      "hemisphere. HL = Hodges–Lehmann location-shift estimate; ",
      "pp = percentage points."
    )
  ) |>
  flextable::set_caption(
    caption = paste0(
      "Table 3. Patient-level rSO₂ changes around successful reperfusion ",
      "and direct comparison of temporal changes between hemispheres."
    )
  ) |>
  flextable::autofit()

ft_table3








# 
# tb2 <- table(avg_post$ratio_signif)
# tb2
# 
# 
# db_row <- db_msr %>% filter(db_msr$patient_id == avg_2$patient_id[1]) 
# db_row
# 
# avg_3
# 
#   
# 
# 
# 
# table(avg_start$ratio_signif)
# 
# sum(is.na(avg_start$Unaffected))
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# # old code 16.09.2026 - cleaning one patient
# # ad row id
# db_msr <- read.csv2(paste0(db.path, "nirs_measurements_clean_v2.csv"))
# 
# db_msr <- db_msr %>% mutate(row_id = row_number())
# 
# # all patients, only angiography_start_point == TRUE
# length(unique(db_msr$patient_id))
# 
# id_true_first_try <- db_msr %>% group_by(patient_id) %>% 
#   filter(first_try_point == TRUE) %>% 
#   dplyr::select(patient_id)
# 
# no_first_try <- (setdiff(as_vector(db_clin$patient_id), as_vector(id_true_first_try)))
# # which patients do not have angiography_start_point == TRUE
# # 29 patients do not have this flag!
# no_first_try
# 
# db_msr %>% group_by(patient_id) %>% 
#   filter(recanalisation_point == TRUE) %>% summarise(n= n())
# 
# # 2 patients do not have flag for recanalisation_point
# id_true <- db_msr %>% group_by(patient_id) %>% 
#   filter(recanalisation_point == TRUE)
# 
# id_norecan <- setdiff(as_vector(db_clin$patient_id), as_vector(id_true))
# id_norecan
# 
# # do 2 patients have 
# db_msr %>% filter(patient_id %in% id_norecan & first_try_point == TRUE) %>% 
#   dplyr::select(patient_id, datetime, first_try_point, recanalisation_point)
# 
# # the patient 1001590120 has wrong year (2025 instead of 2015)
# # correct recan time is: 2015-10-22 12:15:00
# row_tochange <- db_msr %>% filter(patient_id == 1001590120 &
#                     datetime >= "2015-10-22 12:15:00"
#                   & datetime <= "2015-10-22 12:16:00") %>% 
#   slice_head(n=1) %>% dplyr::select(row_id)
# as.numeric(row_tochange)
# 
# db_msr[db_msr$row_id == as.numeric(row_tochange),]$recanalisation_point <- TRUE
# 
# # the patient 1001249114 has wrong date order
# row_tochange <- db_msr %>% filter(patient_id == 1001249114 &
#                                     datetime >= "2013-07-01 09:15:00"
#                                   & datetime <= "2013-07-01 09:16:00") %>% 
#   slice_head(n=1) %>% dplyr::select(row_id)
# as.numeric(row_tochange)
# 
# db_msr[db_msr$row_id == as.numeric(row_tochange),]$recanalisation_point <- TRUE
# 
# no_open <- db_clin %>% filter(opening_datetime == "") %>% dplyr::select(patient_id)
# no_open <- as_vector(no_open$patient_id)
# no_open
# # change 13 patients to no opening time!
# db_msr <- db_msr %>% mutate(
#   recanalisation_point = case_when(patient_id %in% no_open ~ FALSE, TRUE ~ recanalisation_point)
# )
# 
# 
# # save new database
# write.csv2(db_msr, paste0(db.path, "nirs_measurements_clean_v2.csv"))
