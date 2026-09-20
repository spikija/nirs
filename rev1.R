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
  pROC
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

# duration
sample_duration <- 10

avg_1 <- db_msr %>% group_by(patient_id, hemisphere_status) %>% 
  filter(elapsed_seconds >= 60*10 & elapsed_seconds <= 60*10*sample_duration) %>% 
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

avg_1

# ------------------------------------------------------------------
# Affected-hemisphere identification
# ------------------------------------------------------------------

tb1 <- table(avg_1$ratio_signif)

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
fig1 <- avg_1 %>% 
  ggplot(aes(x = 0, y = difference)) + 
  geom_dotplot(binwidth = 0.5, dotsize = 1.5, binaxis = "y", stackdir = "center") 
fig1 

library(dplyr)
library(ggplot2)

# horizontal layout positions
x_dot  <- -0.75
x_text <-  -0.35

fig1 <- avg_1 %>%
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
  
  # # zero difference
  # geom_hline(
  #   yintercept = 0,
  #   linetype = "dashed",
  #   linewidth = 0.6
  # ) +
  
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
    y = 9,
    label = "Higher on\naffected side",
    hjust = 0,
    size = 4
  ) +
  
  annotate(
    "text",
    x = x_text,
    y = 0.0,
    label = "Similar between hemispheres\n(|ΔrSO2| ≤ 4 percentage points)",
    hjust = 0,
    size = 3.7
  ) +
  
  annotate(
    "text",
    x = x_text,
    y = -9,
    label = "Higher on\nunaffected side",
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

rmarkdown::render(
  "c:/sci/rcode/nirs/nirs_rev1.Rmd",
  output_file = "affected_hemisphere_results.html",
  envir = globalenv()
)

# Q2 : 

avg_2 <- db_msr %>% group_by(patient_id, elapsed_seconds) %>% 
  filter(recanalisation_point == TRUE) %>% dplyr::select(elapsed_seconds) %>% 
  mutate(start_value = elapsed_seconds + (60*sample_duration), upper_value = elapsed_seconds + (60*10*sample_duration))
avg_2 


calc_avg <- function(db, db_search) {
  print(db_search$patient_id)
  db <- db %>% filter(db$patient_id == db_search$patient_id)
  db_row <- db %>% group_by(patient_id, hemisphere_status) %>% 
    filter(elapsed_seconds >= db_search$elapsed_seconds & elapsed_seconds <= db_search$upper_value) %>% 
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

avg_3 <- tibble()
for (i in 1:nrow(avg_2)) { 
  #db_msr %>% filter(db$patient_id == db_search$patient_id)
  #avg_3 <- add_row(calc_avg(db_msr, avg_2[i,]))
  avg_3 <- bind_rows(avg_3, (calc_avg(db_msr, avg_2[i,])))
}
avg_3

tb2 <- table(avg_3$ratio_signif)
tb2


db_row <- db_msr %>% filter(db_msr$patient_id == avg_2$patient_id[1]) 
db_row

avg_3

  



table(avg_1$ratio_signif)

sum(is.na(avg_1$Unaffected))










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
