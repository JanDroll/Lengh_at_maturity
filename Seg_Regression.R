library(tidyverse)
library(segmented)

## Read in the data and calculate the growth rate slope ----

gray_mr <- read_csv2("Gray_mrt.csv") %>% 
  mutate(date_m = as.Date(date_m, format = "%d.%m.%y"),
         date_r = as.Date(date_r, format = "%d.%m.%y"),
         dt = as.numeric(difftime(date_r, date_m, units = "days")),
         dL = Lr - Lm,
         slope = dL / dt)

## Estimating length at maturation ----
# Estimate the average trajectory among individuals
m0 <- lm(slope ~ Lm, data = gray_mr)

# Estimate breakpoint with the segmented package
m_seg <- segmented(m0, seg.Z = ~Lm)
summary(m_seg)

bp <- m_seg$psi[,"Est."]

confint(m_seg) # Getting CI intervals for the breakpoint

# Estimating CIs for the model
pred_df <- tibble(Lm = seq(min(gray_mr$Lm), max(gray_mr$Lm), length.out = 400))

pred <- predict(m_seg, newdata = pred_df, se.fit = TRUE)

pred_df %<>% 
  mutate(fit = pred$fit,
         se = pred$se.fit,
         l_ci = pred$fit - 1.96 * pred$se.fit,
         u_ci = pred$fit + 1.96 * pred$se.fit)

immature <- filter(pred_df, Lm <= bp + .3)
mature <- filter(pred_df, Lm >= bp - .3)

## Plot the model results ----
seg_plot <- ggplot() +
  geom_point(data = gray_mr, aes(x = Lm, y = slope, size = dt), shape = 21, color = "black", fill = "grey20", alpha = .4) +
  geom_ribbon(data = immature, aes(x = Lm, ymin = l_ci, ymax = u_ci), alpha = .2) +
  geom_ribbon(data = mature, aes(x = Lm, ymin = l_ci, ymax = u_ci), alpha = .2) +
  geom_line(data = immature, aes(x = Lm, y = fit), size = 1, col = "#03879A") +
  geom_line(data = mature, aes(x = Lm, y = fit), size = 1, col = "#FCA800") +
  geom_point(aes(x = 294.1579, y = 0.2959574401), size = 5, fill = "#FC1D00", col = "black", shape = 23) +
  theme_light(base_size = 13) +
  annotate("curve" , xend = 297, yend = .303, x = 340, y = .38, curvature = .28, arrow = arrow(length = unit(0.3, "cm"), type = "closed"), colour = "black") +
  annotate("text", x = 340, y = .38, label = "Estimated length\n at maturation\n = 294.4 mm", col = "#FC1D00", hjust = 0, size = 4.5) +
  labs(size = expression(Delta*"t"~"(Days)"), x = "Total length (mm) at capture", y = "Slope") +
  theme(legend.position = c(.09, .18),
        legend.background = element_rect(linetype = "solid", color = "grey70")) "solid", color = "grey70"))