## check that uninfected do not differ much in behaviour


data_x = data %>% group_by(fullID) %>% 
  mutate(has_nematodes=mean(Nematode_live_all>0,na.rm=TRUE)==1) %>% #needed to circumvent the NAs on the second behaviour
  ungroup()

mod_mvt_infec <- brm(
  bf(scale_mvt | subset(is_active == TRUE & has_nematodes==FALSE) ~ ## we exclude non-moving snails
       0 + band_number + band_number:population +
       scale_size + (1 | series) + (1 | fullID)),
  data = data_x, chains = 4, iter = 2000, warmup = 1000,
  prior = c(
    set_prior("normal(0,1)", class = "b"),
    set_prior("normal(0,1)", class = "sd"),
    set_prior("normal(0,1)", class = "sigma")
  ),
  backend = "cmdstanr", seed = 42, control = list(adapt_delta = 0.99, max_treedepth = 15)
)

mod_mvt_infec2 <- brm(
  bf(scale_mvt | subset(is_active == TRUE & has_nematodes==FALSE) ~ ## we exclude non-moving snails
       0 + band_number +
       scale_size + (1 | series) + (1 | fullID)),
  data = data_x, chains = 4, iter = 2000, warmup = 1000,
  prior = c(
    set_prior("normal(0,1)", class = "b"),
    set_prior("normal(0,1)", class = "sd"),
    set_prior("normal(0,1)", class = "sigma")
  ),
  backend = "cmdstanr", seed = 42, control = list(adapt_delta = 0.99, max_treedepth = 15)
)

loo(mod_mvt_infec,mod_mvt_infec2)

preds_mvt_x <- data_x %>%
  select(population, band_number, mean_mvt, sd_mvt) %>%
  distinct() %>%
  mutate(scale_size = 0, session = 1, is_active = TRUE, has_nematodes=FALSE) %>%
  add_epred_draws(mod_mvt_infec, re_formula = NA) %>%
  ungroup() %>%
  mutate(population = fct_recode(population,
                                 `shaded habitat (no nematodes)` = "shaded",
                                 `open habitat (nematodes)` = "sun_exposed"
  ))


preds_mvt_x %>%   
  mutate(group = paste(band_number, ", ", population, sep = "")) %>%
  compare_levels(.epred, by = group) %>%
  mean_hdi() %>% 
  print(n=Inf)


mod_food_infec <- brm(
  bf(scale_food | subset( has_nematodes==FALSE)~
       0 + band_number + band_number:population +
       scale_size + (1 | series) + (1 | fullID)),
  data = data_x, chains = 4, iter = 4000, warmup = 2000,
  prior = c(
    set_prior("normal(0,1)", class = "b"),
    set_prior("normal(0,1)", class = "sd")
  ),
  backend = "cmdstanr", seed = 42, control = list(adapt_delta = 0.99, max_treedepth = 15)
)

mod_food_infec2 <- brm(
  bf(scale_food | subset( has_nematodes==FALSE)~
       0 + band_number +
       scale_size + (1 | series) + (1 | fullID)),
  data = data_x, chains = 4, iter = 4000, warmup = 2000,
  prior = c(
    set_prior("normal(0,1)", class = "b"),
    set_prior("normal(0,1)", class = "sd")
  ),
  backend = "cmdstanr", seed = 42, control = list(adapt_delta = 0.99, max_treedepth = 15)
)


loo(mod_food_infec,mod_food_infec2)

preds_food_x <- data_x %>%
  select(population, band_number, mean_mvt, sd_mvt) %>%
  distinct() %>%
  mutate(scale_size = 0, session = 1, is_active = TRUE, has_nematodes=FALSE) %>%
  add_epred_draws(mod_food_infec, re_formula = NA) %>%
  ungroup() %>%
  mutate(population = fct_recode(population,
                                 `shaded habitat (no nematodes)` = "shaded",
                                 `open habitat (nematodes)` = "sun_exposed"
  ))


preds_food_x %>%   
  mutate(group = paste(band_number, ", ", population, sep = "")) %>%
  compare_levels(.epred, by = group) %>%
  mean_hdi() %>% 
  print(n=Inf)

### figure showing the correl food - nemat?
data %>% group_by(fullID) %>% 
  mutate(nematodes=mean(Nematode_live_all,na.rm=TRUE)) %>% #needed to circumvent the NAs on the second behaviour
  ggplot()+
  geom_point(aes(nematodes,scale_food))+facet_wrap(~population)
## hard to group the two food obs of a snail on this plot (versus other figures)

### showing using the average? (for demo only, maybe as info to reviewer)
data %>% group_by(fullID,population) %>% 
  summarise(nematodes=mean(Nematode_live_all,na.rm=TRUE),
            mean_food=mean(scale_food,na.rm=TRUE)) %>% #needed to circumvent the NAs on the second behaviour
  ggplot()+
  geom_point(aes(nematodes,mean_food))+facet_wrap(~population)



