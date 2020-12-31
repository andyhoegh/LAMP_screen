##### Plotting functions

lengthen <- function (output.list) {
  inf <- c(output.list$inf) #takes an outputed matrix and converts to vector
  symptrep <- c(output.list$symptrep) #""
  true_test_positives <- c(output.list$true_test_positives) #""
  case <- c(output.list$case) #""
  isolation <- c(output.list$isolation) #""
  quarantine <- c(output.list$quarantine) #""
  day <- rep(1:100, 500) # adds day for each simulation in vector format
  long <- data.frame(inf, symptrep, true_test_positives, 
                     case, isolation, quarantine, day) 
  return(long)
}

running_tots <- function(output) {
  return(output %>% 
           group_by(group, tests, compliance) %>% 
           mutate(cum_cases.on = cumsum(new.cases.on),
                  cum_reporting.symptoms.on = cumsum(reporting.symptoms.on),
                  cum_all.symptomatics.on = cumsum(all.symptomatics.on),
                  cum_all.asymptomatics.on = cumsum(positive.asympt.on),
                  cum_cases.off = cumsum(new.cases.off),
                  cum_reporting.symptoms.off = cumsum(reporting.symptoms.off),
                  cum_all.symptomatics.off = cumsum(all.symptomatics.off),
                  cum_all.asymptomatics.off = cumsum(positive.asympt.off)
           )
  )
}

plotting.function <- function(output, var, log.10 = T) {
  if (log.10 == T) {
    plot(ggplot(output,aes_string(x = "day", y = var, color = "factor(tests)")) +
           geom_line(aes(group = group), alpha = 0.025) +
           geom_smooth() +
           facet_grid(compliance~.) +
           labs(
             x = "Day of School Year") +
           scale_y_log10(n.breaks = 10)+
           theme_classic())
  }
  if (log.10 == F) {
    plot(ggplot(output,aes_string(x = "day", y = var, color = "factor(tests)")) +
           geom_line(aes(group = group), alpha = 0.025) +
           geom_smooth() +
           facet_grid(compliance~.) +
           labs(
             x = "Day of School Year") +
           theme_classic())
  }
}

plotting.summaries <- function (out) {
  out$lower <- out$mean - (1.96*out$sd)
  out$upper <- out$mean + (1.96*out$sd)
  out %>% 
    ggplot(aes(x=factor(tests), y = mean, fill = factor(tests))) +
    geom_bar(stat = "identity") +
    geom_errorbar(aes(ymin = lower, ymax = upper)) +
    facet_grid(.~compliance) + 
    scale_y_continuous(n.breaks = 10)+
    theme_classic()
}