#### University Simulation Functions


uni_sim <- function(tst = 500, test.timeline = c("Initial", "Sustained", "Both"),
                    compliance = 0.75, init.prev = .03, ppn_sympt = .8, 
                    care.seeking = 0.5, R0.on = 3, R0.off = 1.5, 
                    test.scenario = c("2 Days","1 Day","No Delay"),
                    sens.pcr = .99, spec.pcr = .99, sens.lamp = c(.8, .9, 1), spec.lamp = .99, 
                    lamp.diagnostic = F, community.intro.daily.on = 1, 
                    community.prob.daily.on = 0.1,
                    community.intro.daily.off = 1, 
                    community.prob.daily.off = 0.1,
                    immunity = 0.1, N0 = 16750, on.campus.prop = .25, 
                    contact.tracing.limit = 100, pooling = 4, pooling.multi = 1,
                    days = 100, sims = 200,
                    days.to.isolate = 10, days.to.quarantine = 10){
  
  Tsim <- as.numeric(days)        # time to simulate over, we only care about start
  sims <- as.numeric(sims)    # number of simulations
  lamp.diagnostic <- lamp.diagnostic
  ppn_sympt <- as.numeric(ppn_sympt)
  compliance <- as.numeric(compliance)
  care.seeking <- as.numeric(care.seeking) 
  init.prev <- as.numeric(init.prev)
  tst <- as.numeric(tst)
  test.scenario <- test.scenario
  sens.pcr <- as.numeric(sens.pcr)
  spec.pcr <- as.numeric(spec.pcr)
  sens.lamp <- as.numeric(sens.lamp)
  spec.lamp <- as.numeric(spec.lamp)
  R0.on <- as.numeric(R0.on)   
  R0.off <- as.numeric(R0.off)
  contact.tracing.limit <- as.numeric(contact.tracing.limit)
  pooling <- as.numeric(pooling)
  pooling.multi <- as.numeric(pooling.multi)
  immunity <- as.numeric(immunity)
  N0 <- as.numeric(N0)
  on.campus.prop <- as.numeric(on.campus.prop)
  community.intro.daily.on <- as.numeric(community.intro.daily.on)
  community.prob.daily.on <- as.numeric(community.prob.daily.on)
  community.intro.daily.off <- as.numeric(community.intro.daily.off)
  community.prob.daily.off <- as.numeric(community.prob.daily.off)
  
  days.to.quarantine <- as.numeric(days.to.quarantine)
  days.to.isolate <- as.numeric(days.to.isolate)
  
  # storage vectors
  RE.on <- R0.on*(1-immunity)  # RE assuming some fraction of population is already immune
  RE.off <- R0.off*(1-immunity) # make argument !!!!!!!!!!!
  beta.normal.on <- RE.on * (1/9)  # calculate Beta
  beta.normal.off <- RE.off * (1/9)  # calculate Beta
  # beta.outbreak <- RE * (1/9) * (1-distancing_reduction)   # calculate Beta
  beta_vec.on <- rep(beta.normal.on,sims)
  beta_vec.off <- rep(beta.normal.off,sims)
  theta <- 1/5  # 5 days from infection to infectious
  gamma_I1I2 <- 1/2 # 2 days asymptomatic infectious
  gamma_I2R <- 1/7 # 7 days infectious (this is probably too short)
  
  # storage 
  inf.on <- matrix(NA,Tsim,sims)  # storage for infectious class
  case.on <- matrix(NA,Tsim,sims) # storage for daily cases
  S.on <- matrix(round(N0 * on.campus.prop * (1-immunity)),1,sims)  # start with 0.85 susceptible # change to variable!!!!!!
  E.on <- matrix(floor(N0 * init.prev * on.campus.prop * (1-immunity)),1,sims)
  I1.on <- matrix(0,1,sims)                 # start with 5 asymptomatic infectious
  I2.on <- matrix(0,1,sims)
  R.on <- matrix((N0 * on.campus.prop * immunity) ,1,sims) # change to variable!!!!!!
  sympt1.on <- matrix(0,1,sims)
  sympt2.on <- matrix(0,1,sims)
  symptrep.on <- matrix(0,1,sims)
  symptactual.on <- matrix(0,1,sims)
  new_contacts.on <- matrix(0,1,sims)
  isolation.on <- matrix(0,1,sims)
  quarantine.on <- matrix(0,1,sims)
  new_cases.on <- matrix(0,1,sims)
  true_test_positives.on <- matrix(0,1,sims)
  total_traces.on <- matrix(0,1,sims)
  comply_test_positives.on <- matrix(0,1,sims)
  N.on <- S.on+E.on+I1.on+I2.on+R.on
  
  inf.off <- matrix(NA,Tsim,sims)  # storage for infectious class
  case.off <- matrix(NA,Tsim,sims) # storage for daily cases
  S.off <- matrix(round(N0 * (1 - on.campus.prop) * (1 - immunity)),1,sims)  # start with 0.85 susceptible # change to variable!!!!!!
  E.off <- matrix(floor(N0 * init.prev * (1 - on.campus.prop) * (1-immunity)),1,sims)
  I1.off <- matrix(0,1,sims)                 # start with 5 asymptomatic infectious
  I2.off <- matrix(0,1,sims)
  R.off <- matrix((N0 * (1 - on.campus.prop) * immunity),1,sims)
  sympt1.off <- matrix(0,1,sims)
  sympt2.off <- matrix(0,1,sims)
  symptrep.off <- matrix(0,1,sims)
  symptactual.off <- matrix(0,1,sims)
  new_contacts.off <- matrix(0,1,sims)
  isolation.off <- matrix(0,1,sims)
  quarantine.off <- matrix(0,1,sims)
  new_cases.off <- matrix(0,1,sims)
  true_test_positives.off <- matrix(0,1,sims)
  total_traces.off <- matrix(0,1,sims)
  comply_test_positives.off <- matrix(0,1,sims)
  symp.pcr <- matrix(0,1,sims)
  asymp.pcr <- matrix(0,1,sims)
  cases.caught <- matrix(0,1,sims)
  N.off <- S.off+E.off+I1.off+I2.off+R.off
  
  atest.wait.3 <- array(0,c(1,sims,10))
  atest.wait.2 <- array(0,c(1,sims,10))
  atest.wait.1 <- array(0,c(1,sims,10))
  contact.wait.3 <- array(0,c(1,sims,10))
  contact.wait.2 <- array(0,c(1,sims,10))
  contact.wait.1 <- array(0,c(1,sims,10))
  
  intro.on <- rbinom(1, community.intro.daily.on, community.prob.daily.on)
  intro.off <- rbinom(1, community.intro.daily.off, community.prob.daily.off)
  
  # distancing_reduction <- 0.5 # if trigger is crossed, NPIs are imposed and transmission is reduced by this fraction
  # qi_trigger <- numeric(sims)
  
  for(ts in 2:Tsim){
    if(test.timeline == "Initial" & ts > 20){
      tests = 0
    }
    if(test.timeline == "Initial" & ts <= 20){
      tests = floor(5*tst) 
    }
    
    if(test.timeline == "Sustained"){
      tests = floor(tst)
    }
    
    if(test.timeline == "Both" & ts > 20){
      tests = floor(tst*2.5)
    }
    if(test.timeline == "Both" & ts <= 20){
      tests = floor(tst*0.625)
    }
    out <- sir_lamp(sims, 
                    S.on[ts-1,], E.on[ts-1,], I1.on[ts-1,], I2.on[ts-1,], R.on[ts-1,], 
                    N.on[ts-1,], sympt1.on[ts-1,], sympt2.on[ts-1,], beta_vec.on, 
                    S.off[ts-1,], E.off[ts-1,], I1.off[ts-1,], I2.off[ts-1,], R.off[ts-1,], 
                    N.off[ts-1,], sympt1.off[ts-1,], sympt2.off[ts-1,], beta_vec.off, 
                    theta, gamma_I1I2, gamma_I2R, delta.t=1, tests, contacts.on = 5, contacts.off = 5,
                    ppn_sympt = ppn_sympt, compliance = compliance, care.seeking = care.seeking,
                    atest.wait.3[ts-1,,],atest.wait.2[ts-1,,],atest.wait.1[ts-1,,],
                    contact.wait.3[ts-1,,], contact.wait.2[ts-1,,],contact.wait.1[ts-1,,],
                    test.scenario, sens.pcr = sens.pcr, spec.pcr = spec.pcr, 
                    sens.lamp = sens.lamp, spec.lamp = spec.lamp, lamp.diagnostic = lamp.diagnostic,
                    contact.tracing.limit = contact.tracing.limit, intro.on, intro.off, pooling,
                    pooling.multi) # call to SIR step function above
    S.on <- rbind(S.on,out[,1])  # update state
    E.on <- rbind(E.on,out[,2])  # update state
    I1.on <- rbind(I1.on,out[,3])  # update state
    I2.on <- rbind(I2.on,out[,4])  # update state
    R.on <- rbind(R.on,out[,5])  # update state
    N.on <- rbind(N.on,N.on[ts-1])  # update state
    sympt1.on <- rbind(sympt1.on,out[,23])  # update state
    sympt2.on <- rbind(sympt2.on,out[,25])  # update state
    symptrep.on <- rbind(symptrep.on,out[,27])  # update state
    symptactual.on <- rbind(symptactual.on,out[,70])
    new_contacts.on <- rbind(new_contacts.on,apply(out[,29:33],1,sum))  # update state
    true_test_positives.on <- rbind(true_test_positives.on, apply(out[,15:16],1,sum))  # update state
    comply_test_positives.on <- rbind(comply_test_positives.on, apply(out[,52:53],1,sum))  # update state
    new_cases.on <- rbind(new_cases.on,out[,6])  # update state
    total_traces.on <- rbind(total_traces.on, apply(out[,39:43],1,sum))# update state
    
    S.off <- rbind(S.off,out[,7])  # update state
    E.off <- rbind(E.off,out[,8])  # update state
    I1.off <- rbind(I1.off,out[,9])  # update state
    I2.off <- rbind(I2.off,out[,10])  # update state
    R.off <- rbind(R.off,out[,11])  # update state
    N.off <- rbind(N.off,N.off[ts-1])  # update state
    sympt1.off <- rbind(sympt1.off,out[,24])  # update state
    sympt2.off <- rbind(sympt2.off,out[,26])  # update state
    symptrep.off <- rbind(symptrep.off,out[,28])  # update state
    symptactual.off <- rbind(symptactual.off,out[,71])
    new_contacts.off <- rbind(new_contacts.off,apply(out[,34:38],1,sum))  # update state
    true_test_positives.off <- rbind(true_test_positives.off, apply(out[,20:21],1,sum))  # update state
    comply_test_positives.off <- rbind(comply_test_positives.off, apply(out[,57:58],1,sum))  # update state
    new_cases.off <- rbind(new_cases.off,out[,12])  # update state
    total_traces.off <- rbind(total_traces.off, apply(out[,44:48],1,sum))# update state

    symp.pcr <- rbind(symp.pcr, out[,132])
    asymp.pcr <- rbind(asymp.pcr, out[,133])
    cases.caught <- rbind(cases.caught, out[,134])
    
    atest.wait.3 <- abind(atest.wait.3, array(out[,72:81], c(1,sims,10)), along = 1)
    atest.wait.2 <- abind(atest.wait.2, array(out[,82:91], c(1,sims,10)), along = 1)
    atest.wait.1 <- abind(atest.wait.1, array(out[,92:101], c(1,sims,10)), along = 1)
    contact.wait.3 <- abind(contact.wait.3, array(out[,102:111], c(1,sims,10)), along = 1)
    contact.wait.2 <- abind(contact.wait.2, array(out[,112:121], c(1,sims,10)), along = 1)
    contact.wait.1 <- abind(contact.wait.1, array(out[,122:131], c(1,sims,10)), along = 1)
    ####################################################################################
    # Total in Isolation/Qurantine
    isolation.on <- rbind(isolation.on,apply(comply_test_positives.on[(max(1,ts-days.to.isolate)):ts,],2,sum) + apply(symptrep.on[(max(1,ts-10)):ts,],2,sum)) # isolate for 10 days
    quarantine.on <- rbind(quarantine.on, apply(new_contacts.on[(max(1,ts-days.to.quarantine)):ts,],2,sum) ) # quarantine for 14 days
    
    isolation.off <- rbind(isolation.off,apply(comply_test_positives.off[(max(1,ts-days.to.isolate)):ts,],2,sum) + apply(symptrep.off[(max(1,ts-10)):ts,],2,sum)) # isolate for 10 days
    quarantine.off <- rbind(quarantine.off, apply(new_contacts.off[(max(1,ts-days.to.quarantine)):ts,],2,sum) )
  }
  
  inf.on <- I1.on+I2.on   # total infectious on campus
  inf.off <- I1.off+I2.off   # total infectious on campus
  case.on <- new_cases.on # daily cases
  case.off <- new_cases.off
  return(list("active.inf.on" = inf.on, 
              "reporting.symptoms.on" = symptrep.on, 
              "all.symptomatics.on" = symptactual.on,
              "positive.asympt.on" = true_test_positives.on, 
              "new.cases.on" = case.on, 
              "isolation.complying.on" = isolation.on,
              "quarantine.complying.on" = quarantine.on,
              "total_traces.on" = total_traces.on,
              "active.inf.off" = inf.off, 
              "reporting.symptoms.off" = symptrep.off, 
              "all.symptomatics.off" = symptactual.off,
              "positive.asympt.off" = true_test_positives.off, 
              "new.cases.off" = case.off, 
              "isolation.complying.off" = isolation.off,
              "quarantine.complying.off" = quarantine.off,
              "total_traces.off" = total_traces.off,
              "symp.pcr" = symp.pcr,
              "asymp.pcr" = asymp.pcr,
              "cases.caught" = cases.caught,
              "tests"= matrix(tst,Tsim,sims),
              "compliance" = matrix(compliance,Tsim,sims),
              "init.prev"= matrix(init.prev,Tsim,sims),
              "ppn_sympt"= matrix(ppn_sympt,Tsim,sims),
              'care.seeking'= matrix(care.seeking,Tsim,sims),
              "R0.on" = matrix(R0.on,Tsim,sims),
              "R0.off" = matrix(R0.off,Tsim,sims),
              "test.scenario" = matrix(test.scenario,Tsim,sims),
              "sens.pcr" = matrix(sens.pcr,Tsim,sims),
              "spec.pcr" = matrix(spec.pcr,Tsim,sims),
              "sens.lamp" = matrix(sens.lamp,Tsim,sims),
              "spec.lamp" = matrix(spec.lamp,Tsim,sims),
              "lamp.diagnostic" = matrix(lamp.diagnostic,Tsim,sims),
              "community.intro.daily.on" = matrix(community.intro.daily.on,Tsim,sims), 
              "community.prob.daily.on" = matrix(community.prob.daily.on,Tsim,sims),
              "community.intro.daily.off" = matrix(community.intro.daily.off,Tsim,sims), 
              "community.prob.daily.off" = matrix(community.prob.daily.off,Tsim,sims),
              "immunity" = matrix(immunity,Tsim,sims), 
              "N0" = matrix(N0,Tsim,sims), 
              "on.campus.prop" = matrix(on.campus.prop,Tsim,sims), 
              "contact.tracing.limit" = matrix(contact.tracing.limit,Tsim,sims), 
              "pooling" = matrix(pooling,Tsim,sims), 
              "pooling.multi" = matrix(pooling.multi,Tsim,sims),
              "days" = matrix(days,Tsim,sims), 
              "sims" = matrix(sims,Tsim,sims),
              "test.timeline" = test.timeline,
              "days.to.isolate" = days.to.isolate,
              "days.to.quarantine" = days.to.quarantine
              ))
}
