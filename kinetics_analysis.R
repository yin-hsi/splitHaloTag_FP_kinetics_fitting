# Install R package "jwtools"
library(devtools)
install_github("Jonas-Wilhelm/jwtools")


## First specify the packages of interest
packages = c("tidyverse", "jwtools",
             "tsoutliers", "furrr", "minpack.lm", "broom", "viridis")
## Now load or install&load all
package.check <- lapply(
   packages,
    FUN = function(x) {
        if (!require(x, character.only = TRUE)) {
        install.packages(x, dependencies = TRUE)
        library(x, character.only = TRUE)
          }
        }
     )

## Needed R packages
library(tidyverse)
library(jwtools)
library(tsoutliers)
library(furrr)
library(minpack.lm)
library(broom)
library(viridis)

# set working directory to path of script

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
#setwd("Desktop/20220727_copy")


# list data files and layout files, remove temporary excel files

data_files   <- list.files("01_raw_data/", full.names = T, pattern = "xlsx$") %>%
  .[!str_detect(., "~\\$")]
layout_files <- list.files("02_plate_layouts/", full.names = T, pattern = "xlsx$") %>%
  .[!str_detect(., "~\\$")]


# read data and layout files
# set max time

data <- read_kinetic_TECAN(data_files, layout_files, 
                           plate_type = "384", 
                           n_cond = 2, 
                           delay = 40)

data <- data %>%
  filter(time < 10*60*60)


# create raw data overview plot

dir.create("03_raw_plots", showWarnings = F)
plot_kinetics_plate(data, plate_type = "384", path = "03_raw_plots", width = 450, height = 300)


# detect and remove strong outliers. sometimes (one in ~10000 points) the TECAN 
# produces strong outliers for some reason

find_outliers <- function(d){
  outliers <- tso(ts(d$value), types = "AO", cval = 10)$times
  d <- d %>%
    add_column(outlier = FALSE)
  d[outliers,]$outlier <- TRUE
  return(d)
}

plan(multisession, workers = availableCores()-1)

data <- data %>%
  group_by(p_well) %>%
  nest() %>%
  ungroup() %>%
  mutate(data = future_map(data, find_outliers, .progress = T)) %>%
  unnest(data)

data$outlier %>% sum()

data <- data %>%
  filter(!outlier) %>%
  select(-outlier) %>%
  select(-cycle, -temp)


# round concentration and convert to factor for plotting

data <- data %>%
  mutate(cond_2_r = round(cond_2, 0) %>% factor() %>% fct_rev())


# plot the data

dir.create("04_plots", showWarnings = F)

p <- data %>%
  ggplot(aes(x = time/60/60, y = value, color = cond_2_r, fill = cond_2_r))+
  stat_summary(fun.data ="mean_sdl", fun.args = list(mult = 1, na.rm = T), geom = "smooth", size = 0.33, se = TRUE, orientation = "x")+
  facet_wrap(~cond_1)+
  scale_color_viridis(discrete = T, direction = -1)+
  scale_fill_viridis(discrete = T, direction = -1, guide = "none")+
  guides(color = guide_legend(reverse = F, override.aes = list(size = 1, fill=NA)))+
  theme_light()+
  labs(x = "time [h]",
       y = "fluorescence polarization [mFP]",
       color = "peptide [nM]")

ggsave("04_plots/plot_01.pdf", p, width = 250, height = 170, units = "mm")


# fit second order reaction model or linear model to data

Xzero <- 0
Yzero <- mean(filter(data, cond_1 == "no_protein", time < 60*60)$value)
A0    <- 5e-10
B0    <- 2.5e-9

newdat <- tibble(time = seq(0, max(data$time), length.out = 500))

# Fit model to EACH replicate (p_well) individually
fit <- data %>%
  group_by(cond_1, cond_2, cond_2_r, p_well) %>% # Added p_well here
  nest() %>%
  mutate(model = map(.x = data, .f = function(x){
    if(max(x$value) > 280){
      nls <- nls(value ~ Plateau + ((Yzero-Plateau)/A0)*(A0*(A0-B0)*exp((A0-B0)*K*(time-Xzero))) / (A0*exp((A0-B0)*K*(time-Xzero))-B0),
                 data = x, 
                 start = list(
                   Plateau = max(x$value),
                   K       = 1e4
                 ), 
                 control=nls.control(minFactor=1e-10, maxiter = 5000, warnOnly = T))
      return(nls)
    }else{
      lm <- lm(value ~ time, data = x)
      return(lm)
    }
  })) %>%
  mutate(class = map_chr(model, ~class(.)[1])) %>% # Robust class check
  mutate(prediction = map(.x = model, ~ tibble(value = predict(., newdat),
                                               time = newdat$time))) %>%
  ungroup()

table(fit$class)


# average data for each time point for plotting fit

data <- data %>%
  group_by(cond_1, cond_2, cond_2_r, time) %>%
  summarise(value = mean(value), .groups = "drop")


# plot data and fit

p <- data %>%
  ggplot(aes(x = time/60/60, y = value, color = cond_2_r, fill = cond_2_r))+
  # The dots are the averaged values
  geom_point(alpha = 0.2, size = 0.2)+ 
  # THE FIX: Include p_well in the data and the grouping aesthetic
  geom_line(data = fit %>% 
              select(cond_1, cond_2_r, p_well, prediction) %>% # Keep p_well
              unnest(prediction),
            aes(group = p_well), # This tells R to draw 2 separate lines
            size = 0.5, 
            alpha = 0.8)+ 
  facet_wrap(~cond_1)+
  scale_color_viridis(discrete = T, direction = -1)+
  scale_fill_viridis(discrete = T, direction = -1, guide = "none")+
  guides(color = guide_legend(reverse = F, override.aes = list(size = 1, fill=NA)))+
  theme_light()+
  labs(x = "time [h]",
       y = "fluorescence polarization [mFP]",
       color = "peptide [nM]",
       caption = "Lines represent individual fits for each replicate well.")

ggsave("04_plots/plot_02.pdf", p, width = 250, height = 170, units = "mm")


# run monte carlo simulation to estimate uncertainties of non linear fit parameters
# estimate uncertainty of linear fit via broom::tidy

plan(multisession, workers = availableCores()-1)

fit <- fit %>% 
  ungroup() %>%
  mutate(param = future_map_if(
    .x = model, 
    .p = class == "nls", 
    .f = ~nls_MC(.x, runs = 1000) %>% summarise_nls_MC(conf.level = 0.95),
    .else = ~(broom::tidy(.x, conf.int = TRUE, conf.level = 0.95)), 
    .progress = T, 
    .options = furrr_options(seed = TRUE))
  )

dir.create("05_fitted_parameters")

fit %>%
  filter(class == "nls") %>%
  select(cond_1, cond_2, cond_2_r, param) %>%
  unnest(param) %>%
  write_csv("05_fitted_parameters/param_nls.csv")

fit %>%
  filter(class == "lm") %>%
  select(cond_1, cond_2, cond_2_r, param) %>%
  unnest(param) %>%
  write_csv("05_fitted_parameters/param_lm.csv")


# calculate initial slopes using derivative of non linear model at t=0 and
# propagate uncertainty 

# 1. Calculate individual slopes per well
slopes_individual <- fit %>%
  mutate(param = map2(param, class, function(par, t){
    if(t == "nls"){
      p  <- par %>% filter(term == "Plateau") %>% pull(estimate)
      k  <- par %>% filter(term == "K") %>% pull(estimate)
      s  <- (p - Yzero)*k*B0
      return(tibble(slope = s))
    }else{
      s  <- par %>% filter(term == "time") %>% pull(estimate)
      return(tibble(slope = s))
    }
  })) %>%
  select(cond_1, cond_2, cond_2_r, p_well, class, param) %>%
  unnest(param)

# 2. Format for final CSV: Calculate Mean, SD, and keep individual replicates
slopes_final <- slopes_individual %>%
  group_by(cond_1, cond_2, cond_2_r, class) %>%
  mutate(rep_id = paste0("Rep_", row_number())) %>% # Create labels Rep_1, Rep_2
  summarise(
    s_mean = mean(slope, na.rm = TRUE),
    s_sd = sd(slope, na.rm = TRUE),
    # This captures the individual values into new columns
    replicates = list(setNames(slope, rep_id)), 
    n = n(),
    .groups = "drop"
  ) %>%
  unnest_wider(replicates)

# Save the detailed file
write_csv(slopes_final, "05_fitted_parameters/slopes.csv")
