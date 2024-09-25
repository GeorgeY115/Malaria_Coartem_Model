#################################################################################
################# EXTENDED COST-EFFECTIVENESS ANALYSIS MODEL ####################

## This model establishes the burden of disease and mortality of a disease on a
## population divided into wealth quintiles. It allows for the implementation of
## a vaccine for which there is assumed life-long immunity.

## Read in data, which consists of, for each quintile, 5-year age groups, proportion
## of population in each age group, incidence, case fatality ratios, and the cost per
## case of treatment.
x <- read.csv("Model Test Data.csv")


## Establish size of total population and per-quintile population.
cohort_size <- 50000
q_cohort_size <- cohort_size/5


## Establish vaccine efficacy (assumed constant across Qs) and a vector with the
## coverages for each quintile. Multiply to get impact per quintile.
vaccine_efficacy <- 0.85
vaccine_coverage <- c(Q1=0.1,Q2=0.2,Q3=0.3,Q4=0.4,Q5=0.5)
vaccine_impact <- vaccine_efficacy*vaccine_coverage


## Establish life expectancy.
life_expectancy <- 90


## Establish the time frame (in years) over which the model will run.
timestep <- seq(0,50,by=1)


## Establish the disability weight used in the DALYs calculations.
disweight = 0.14


## Create vector which lists each quintile (simlpy 1-5). This is used for the outer
## loop, so results can be quintile-specific.
qunique <- unique(x$quintile)


## Create empty dataframe in which results will be stored.
outputtab <- data.frame(quintile=numeric(0),cases=numeric(0),deaths=numeric(0),
                        dalys=numeric(0),cost=numeric(0),postcases=numeric(0),
                        postdeaths=numeric(0),postdalys=numeric(0),postcost=numeric(0))


################# Start the model.
## The outer loop iterates through each quintile.
for (q in qunique){
  
  
  ## Subset the main data into one quintile per iteration of outer loop.
  subq <- x[x$quintile == q, ]

  
  ## Set all the outputs to 0 between each quintile. Results stored at the end
  ## of each one. Discount at start should be 1 (100% - results are multiplied by this).
  total_cases <- 0
  total_deaths <- 0
  total_cost <- 0
  total_dalys <- 0
  total_cases_postvax <- 0
  total_deaths_postvax <- 0
  total_cost_postvax <- 0
  total_dalys_postvax <- 0
  discount <- 1
  
  
  ## The previously produced vaccine impact vector has the value which is relevant
  ## to the current quintile iteration selected.
  vax=vaccine_impact[q]

  
  ## The mid loop iterates through the time steps over which the model runs.
  for (t in timestep){

    
    ## The inner loop iterates through each 5-year age group. 
    for (i in seq_along(subq$agegroup)) {
      
      
      ## 'w' here represents whether those vaccinated while <5 are within the age
      ## of the current inner loop iteration, thereby being affected by at least
      ## partial immunity (according to the vax variable). If the programme began
      ## 8 years ago, for example, those 13 years old or younger will be vaccinated,
      ## so the <5 and 5-9 categories should have the full vax impact applied, the
      ## 10-14 category should have 80% of it applied, and any older age group should
      ## not yet experience any artificial immunity. 
      w = i-((t+1)/5)
      
      if (w>1) {
        unaff_prop <- 1
      } else if (w<0) {
        unaff_prop <- 0
      } else {
        unaff_prop <- w
      }

      
      ## Calculate the number of person years within each age category.
      person_years <- q_cohort_size * subq$ageprop[i]
      
     
      ## Calculate the number of cases pre- and post-vaccine, the latter having
      ## the effect of the previous if statement and the vax variable applied.
      cases <- person_years * subq$incidence[i]
      casespostvax <- person_years * subq$incidence[i]*(unaff_prop*(1-vax))

      
      ## Calculate the number of deaths based on cases and case fatality rate.
      deaths <- cases * subq$cfr[i]
      deathspostvax <- casespostvax * subq$cfr[i]
      
      
      ## Calculate costs, encorporating the discount rate of the current time step.
      newcost <- (cases * subq$casecost[i])*discount
      newcostpostvax <- (casespostvax * subq$casecost[i])*discount
      
      
      ## Calculate DALYs, encorporating the discount rate of the current time step.
      dalys <- (cases * disweight + cases * (life_expectancy-subq$avage[i]))*discount
      dalyspostvax <- (casespostvax * disweight + casespostvax * (life_expectancy-subq$avage[i]))*discount
    
      ## Keep a counter for each inner loop iteration for totals.
      total_cases <- total_cases + cases
      total_deaths <- total_deaths + deaths
      total_cost <- total_cost + newcost
      total_dalys <- total_dalys + dalys
      
      
      ## Do the same for the post-vaccine ones.
      total_cases_postvax <- total_cases_postvax + casespostvax
      total_deaths_postvax <- total_deaths_postvax + deathspostvax
      total_cost_postvax <- total_cost_postvax + newcostpostvax
      total_dalys_postvax <- total_dalys_postvax + dalyspostvax

    }
    
    ## At the end of each time iteration, add a 3% discount.
    discount <- discount * 0.97

  }
  
  ## At the end of each quintile iteration, add the totals to a new row, and combine
  ## this with the output table which was constructed before the model started running.
  newrow <- c(q,total_cases,total_deaths,total_dalys,total_cost,total_cases_postvax,
              total_deaths_postvax,total_dalys_postvax,total_cost_postvax)
  outputtab <- rbind(outputtab,newrow)

  
  ## Name the output table columns.
  colnames(outputtab) <- c("Quintile", "(Pre) Cases", "(Pre) Deaths", "(Pre) DALYs", "(Pre) Cost", "(Post) Cases", "(Post) Deaths", "(Post) DALYs", "(Post) Cost")
}


## Calculate percentage decrease in cases for each quintile following vaccine.
for (i in outputtab){
  Impact <- 100-(outputtab$`(Post) Cases`/outputtab$`(Pre) Cases`)*100
}

## Add impact column to final output.
outputtab <- cbind(outputtab,impact)


## Display the results.
outputtab


