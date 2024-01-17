library(tidyverse)
library(deSolve)
library(FME)

# parameters
pars <- c(alpha = 1, beta = 0.2, delta = 0.5, gamma = 0.2)
# initial state 
init <- c(x = 1, y = 2)
# times
times <- seq(0, 100, by = 1)

deriv <- function(t, state, pars) {
  with(as.list(c(state, pars)), {
    d_x <- alpha * x - beta * x * y
    d_y <- delta * beta * x * y - gamma * y
    return(list(c(x = d_x, y = d_y)))
  })
}
lv_results <- ode(init, times, deriv, pars)

lv_model <- function(pars, times = seq(0, 50, by = 1)) {
  # initial state 
  state <- c(x = 1, y = 2)
  # derivative
  deriv <- function(t, state, pars) {
    with(as.list(c(state, pars)), {
      d_x <- alpha * x - beta * x * y
      d_y <- delta * beta * x * y - gamma * y
      return(list(c(x = d_x, y = d_y)))
    })
  }
  # solve
  ode(y = state, times = times, func = deriv, parms = pars)
}
lv_results <- lv_model(pars = pars, times = seq(0, 50, by = 0.25))

my_palette <- c("springgreen4", "cornflowerblue")
lv_results %>% 
  data.frame() %>% 
  gather(var, pop, -time) %>% 
  mutate(var = if_else(var == "x", "Prey", "Predator")) %>% 
  ggplot(aes(x = time, y = pop)) +
  geom_line(aes(color = var)) +
  scale_color_manual(values = brewer.pal(4, "Paired")[c(2,4)]) +
  labs(title = " ",
       x = "Time", y = "Population Abundance") +
  theme_minimal()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank())