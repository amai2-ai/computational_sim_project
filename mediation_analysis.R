# please work github!
# To install MEPS package in R, you need to do a couple of things.

### Step 1: install all needed packages
#install.packages(“devtools”)
#devtools::install_github("e-mitchell/meps_r_pkg/MEPS")
### step 3: first we need to load the proper packages and libraries.
library("MEPS") ## You need to load the library every time you restart R
library(gt)
library(gtsummary)
library(commonmark)
library("dplyr")          # Data wrangling
library("expss")          # Relabel variables
library("bda")            # Perform Sobel text
library("mediation")        # Perform the bootstrap approach


#### There are two ways to load data from AHRQ MEPS website, I will be using the method of defining year and type.
hc2021 = read_MEPS(year = 2021, type = "FYC")

### Step 5: Change column names to lowercase
names(hc2021) <- tolower(names(hc2021))

### Step 6: Select specific variables from MEPS dataset
### 2021
hc2021p = hc2021 %>%
  rename(
    workdays = ddnwrk21,
    diabetes = diabdx_m18,
    health_status = rthlth31) %>%
  dplyr::select(                        ## NOTE: there is a weird issue when you don't involve dplyr::select bc MASS is preferred 
    dupersid, 
    workdays, 
    diabetes, 
    health_status,
    sex)
hc2021p$year <- 2021

### Step 7: We want to clean the data so we don’t include any missing or NA responses)
hc2021p = hc2021p %>%
  filter(workdays >= 0,
         diabetes >= 1,
         health_status >= 1)

# "No diabetes" = value of 0 to simplify.
hc2021p$diabetes[hc2021p$diabetes == 2] = 0

### Step 8: Convert to factor and label
hc2021p$sex <- factor(hc2021p$sex,
                      levels = c(1, 2),
                      labels = c("Male", "Female"))

hc2021p$health_status <- factor(hc2021p$health_status, 
                                levels = c(1, 2, 3, 4, 5),
                                labels = c("Excellent", "Very good", "Good", "Fair", "Poor"))

hc2021p$diabetes <- factor(hc2021p$diabetes,
                           levels = c(0, 1),
                           labels = c("No diabetes", "Diabetes"))

### Step 9: Relabel the variables
hc2021p <- apply_labels(hc2021p,
                        diabetes = "Diabetes status", 
                        workdays = "Days missed from work", 
                        health_status = "Perceived health status", 
                        sex = "Sex")





# Descriptive Analysis
hc2021p %>%
  tbl_summary(include = c(health_status, sex, workdays), 
              by = diabetes, 
              statistic = list(all_continuous() ~ "{mean} ({sd})")) %>%
  modify_header(label = "**Variable***") %>%
  bold_labels()

# Regression model – Total effect
direct.model <- glm(workdays ~ diabetes, data = hc2021p, family = gaussian(link = "identity"))
round(cbind(coef(direct.model), confint(direct.model)), 3)
tbl_regression(direct.model, estimate_fun = ~ style_number(.x, digits = 3))

indirect.model1 <- glm(as.numeric(health_status) ~ diabetes, data = hc2021p, family = gaussian(link = "identity")) ## We need to make sure to change the `health_status` variable to numeric to be consistent of the linear regression model, this goes for other variables as well.
round(cbind(coef(indirect.model1), confint(indirect.model1)), 3)
tbl_regression(indirect.model1, estimate_fun = ~ style_number(.x, digits = 3))

indirect.model2 <- glm(workdays ~ diabetes + as.numeric(health_status), data = hc2021p, family = gaussian(link = "identity"))
round(cbind(coef(indirect.model2), confint(indirect.model2)), 3)
tbl_regression(indirect.model2, estimate_fun = ~ style_number(.x, digits = 3))

t1 <- tbl_regression(direct.model, estimate_fun = ~ style_number(.x, digits = 3)) %>% modify_column_hide(columns = p.value) %>% modify_column_hide(ci)
t2 <- tbl_regression(indirect.model1, estimate_fun = ~ style_number(.x, digits = 3)) %>% modify_column_hide(columns = p.value) %>% modify_column_hide(ci)
t3 <- tbl_regression(indirect.model2, estimate_fun = ~ style_number(.x, digits = 3)) %>% modify_column_hide(columns = p.value) %>% modify_column_hide(ci)
# the second model helps us compare the results between total effect and direct + indirect. The total effect is equal to the direct and indirect effects. 
tbl_merge(
  tbls = list(t1, t2, t3),
  tab_spanner = c("**Total effect**", "**First-part indirect effect**", "**Second-part direct + indirect effects**")
)
# This establishes that perceived health status is explaining some effect between diabetes and workdays but now we need to perform the inferential test to see if the phenomenon is significant. This is done using Sobel and Bootstrapping.
# evaluating the Proportion of total effect mediated by mediator variable. This is where the bda package comes in handy.
## mediation.test(mv,iv,dv), where mv is the mediator variable, iv is the independent varible, and dv is the dependent variable)
mediation <- mediation.test(as.numeric(hc2021p$health_status), hc2021p$diabetes, hc2021p$workdays)
round(mediation, 3)
# the z-value is 13.84, and the p-value is < 0.001 which is statistically significant.

# after using the sobel we also can use the Bootstrap approach. This is where the package mediation comes in.
hc2021p$health_status <- as.numeric(hc2021p$health_status)

## we will be using both First-part indirect effect model and Second-part direct + indirect effect model:
indirect.model1 <- glm(health_status ~ diabetes, data = hc2021p, family = gaussian(link = "identity"))
indirect.model2 <- glm(workdays ~ diabetes + health_status, data = hc2021p, family = gaussian(link = "identity"))
## Mediation analysis with 1000 simulations
mediation.results <- mediate(indirect.model1, indirect.model2, treat = 'diabetes', mediator = 'health_status', boot = TRUE, sims = 1000)
summary(mediation.results)
# we are most interested in the ACME (or average causal mediation effect”) which tells us how much of the total effect is explained by the mediator. Our estimated ACME has a p-value less than 0.001 which is significant. 

# now we install.packages(c("gganimate", "gifski", "ggplot2", "dplyr", "scales")) and also load them.
library(ggplot2)
library(gganimate)
library(gifski)
library(dplyr)

# we animate two different graphs.
# ANIMATION 1 is a Bootstrap distribution of the indirect effect (ACME)
# builds up simulation by simulation, with a running mean
med <- mediation.results

sims_df <- data.frame(draw = seq_along(med$d0.sims),
                      acme = med$d0.sims)

frames <- seq(10, nrow(sims_df), by = 10)   # one frame every 10 simulations

anim_df <- bind_rows(lapply(frames, function(k) {
  sims_df %>%
    filter(draw <= k) %>%
    mutate(frame = k, running_mean = mean(acme))
}))

# Fix the axes so the plot doesn't jump around between frames
bw     <- diff(range(sims_df$acme)) / 30
breaks <- seq(min(sims_df$acme) - bw, max(sims_df$acme) + bw, by = bw)
max_ct <- max(hist(sims_df$acme, breaks = breaks, plot = FALSE)$counts) * 1.15

p_anim1 <- ggplot(anim_df, aes(x = acme)) +
  geom_histogram(breaks = breaks, fill = "#E6550D", color = "white", alpha = 0.85) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  geom_vline(aes(xintercept = running_mean), color = "#08519C", linewidth = 1) +
  coord_cartesian(xlim = range(breaks), ylim = c(0, max_ct)) +
  labs(title = "Bootstrap Distribution of the Indirect Effect (ACME)",
       subtitle = "Simulation {current_frame} of {max(frames)}  |  blue line = running mean",
       x = "Indirect effect: Diabetes -> Health status -> Days missed",
       y = "Count") +
  theme_minimal(base_size = 14) +
  transition_manual(frame)

anim1 <- animate(p_anim1, nframes = length(frames), fps = 10,
                 width = 800, height = 500, renderer = gifski_renderer())
anim1
anim_save("animation1_bootstrap_acme.gif", animation = anim1)
cat("GIF saved in:", getwd(), "\n")


# now we move on to ANIMATION 2: Perceived health status distribution
# morphing between "No diabetes" and "Diabetes"
# ------------------------------------------------------------
plot_data <- hc2021p
plot_data$health_status_f <- factor(plot_data$health_status,
                                    levels = 1:5,
                                    labels = c("Excellent", "Very good", "Good", "Fair", "Poor"))

health_pct <- plot_data %>%
  count(diabetes, health_status_f) %>%
  group_by(diabetes) %>%
  mutate(pct = n / sum(n)) %>%
  ungroup()

p_anim2 <- ggplot(health_pct, aes(x = health_status_f, y = pct, fill = diabetes)) +
  geom_col(width = 0.7, show.legend = FALSE) +
  geom_text(aes(label = scales::percent(pct, accuracy = 0.1)), vjust = -0.5, size = 4.5) +
  scale_y_continuous(labels = scales::percent, limits = c(0, max(health_pct$pct) * 1.15)) +
  scale_fill_manual(values = c("No diabetes" = "#6BAED6", "Diabetes" = "#E6550D")) +
  labs(title = "Perceived Health Status: {closest_state}",
       x = "Perceived health status", y = "Percent within group") +
  theme_minimal(base_size = 14) +
  transition_states(diabetes, transition_length = 2, state_length = 2) +
  ease_aes("cubic-in-out")

anim2 <- animate(p_anim2, fps = 20, duration = 6,
                 width = 800, height = 500, renderer = gifski_renderer())
anim2
anim_save("animation2_health_status.gif", animation = anim2)
cat("GIF saved in:", getwd(), "\n")

