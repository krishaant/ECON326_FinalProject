d |> select(Gini, average_AQHI, high_AQHI_days, avg_income, male, vismin, 
            immigrant, married, not_eng_fr, employed, university, 
            college, labor_force) |> st(digits = 6) #generate summary table

d |> select(Gini, average_AQHI, high_AQHI_days, avg_income, male, vismin, 
            immigrant, married, not_eng_fr, employed, university, 
            college, labor_force) |> na.omit() |> cor() #generate correlation matrix



mainreg <- lm(Gini ~ average_AQHI + average_AQHI:income_cat + male + immigrant
              + married + employed + university + college + labor_force, data = d)

summary(mainreg)
coeftest(mainreg, vcov = vcovHC(mainreg, type = "HC1"))



altreg1 <- lm(Gini ~ average_AQHI + I(average_AQHI^2) + average_AQHI:income_cat
              + male + immigrant + married + employed + university + college 
              + labor_force, data = d)

summary(altreg1)
coeftest(altreg1, vcov = vcovHC(altreg1, type = "HC1"))



altreg2 <- lm(Gini ~ average_AQHI + I(average_AQHI^2) + male + immigrant 
              + married + employed + university + college+labor_force, data = d)

summary(altreg2)
coeftest(altreg2, vcov = vcovHC(altreg2, type = "HC1"))



altreg3 <- lm(Gini ~ average_AQHI + average_AQHI:income_cat + male + immigrant
              + married + employed + labor_force, data = d)

summary(altreg3)
coeftest(altreg3, vcov = vcovHC(altreg3, type = "HC1"))



altreg4 <- lm(Gini ~ average_AQHI + average_AQHI:income_cat + male + immigrant
              + married + employed + university + college + labor_force 
              + factor(CMATYPE), data = d)

summary(altreg4)
coeftest(altreg4, vcov = vcovHC(altreg4, type = "HC1"))


altreg5 <- lm(Gini ~ average_AQHI + average_AQHI:income_cat + male + immigrant
              + married + employed + university + college + labor_force 
              + LANDAREA, data = d)

summary(altreg5)
coeftest(altreg5, vcov = vcovHC(altreg5, type = "HC1"))