library(here)
devtools::install_github("jtextor/dagitty/r")
library(dagitty)

g1 <- dagitty('
dag {
"Age/sex/GP practice" [adjusted,pos="0.179,-1.255"]
"Attending GP practice" [adjusted,pos="-0.332,0.834"]
"Calendar period" [adjusted,pos="0.105,-1.100"]
"Genetic risk" [latent,pos="-1.694,-1.148"]
"Harmful alcohol use" [pos="-0.078,0.144"]
"High dose oral glucocorticoids" [pos="0.172,0.599"]
"Physical activity" [latent,pos="-1.394,0.113"]
"Rash position on body" [latent,pos="-1.352,0.298"]
"Sleep problems" [pos="0.126,0.458"]
"anxiety or depression" [outcome,pos="0.525,-0.275"]
Asthma [adjusted,pos="0.035,-0.952"]
BMI [pos="-0.138,0.041"]
Deprivation [adjusted,pos="0.296,-1.416"]
Eczema [exposure,pos="-1.717,-0.283"]
Ethnicity [latent,pos="-1.546,-0.856"]
Multimorbidity [adjusted,pos="-0.029,-0.804"]
Smoking [pos="0.014,0.309"]
"Age/sex/GP practice" -> "anxiety or depression" [pos="0.842,-0.793"]
"Age/sex/GP practice" -> Eczema [pos="-0.396,-1.079"]
"Calendar period" -> "anxiety or depression" [pos="0.864,-0.848"]
"Calendar period" -> Eczema [pos="-0.960,-0.765"]
"Genetic risk" -> "anxiety or depression" [pos="-1.172,-0.514"]
"Genetic risk" -> Eczema [pos="-1.740,-0.713"]
"Harmful alcohol use" -> "anxiety or depression" [pos="0.006,0.085"]
"High dose oral glucocorticoids" -> "anxiety or depression" [pos="0.521,0.584"]
"Physical activity" -> "anxiety or depression" [pos="-1.256,-0.174"]
"Rash position on body" -> "anxiety or depression" [pos="-1.218,-0.172"]
"Sleep problems" -> "anxiety or depression" [pos="0.176,0.150"]
"anxiety or depression" -> "Attending GP practice" [pos="0.920,0.922"]
Asthma -> "anxiety or depression" [pos="0.828,-0.721"]
Asthma -> Eczema [pos="-0.964,-0.752"]
BMI -> "anxiety or depression" [pos="-0.018,-0.068"]
Deprivation -> "anxiety or depression" [pos="0.828,-0.874"]
Deprivation -> Eczema [pos="-0.357,-1.159"]
Eczema -> "Attending GP practice" [pos="-1.810,0.837"]
Eczema -> "Harmful alcohol use" [pos="-0.350,0.078"]
Eczema -> "High dose oral glucocorticoids" [pos="-0.526,0.375"]
Eczema -> "Physical activity" [pos="-1.595,0.026"]
Eczema -> "Rash position on body" [pos="-1.687,0.312"]
Eczema -> "Sleep problems" [pos="-0.533,0.344"]
Eczema -> "anxiety or depression"
Eczema -> BMI [pos="-0.346,0.004"]
Eczema -> Smoking [pos="-0.244,0.181"]
Ethnicity -> "anxiety or depression" [pos="-0.981,-0.516"]
Ethnicity -> Eczema [pos="-1.747,-0.695"]
Multimorbidity -> "anxiety or depression" [pos="0.768,-0.619"]
Multimorbidity -> Eczema [pos="-0.932,-0.774"]
Smoking -> "anxiety or depression" [pos="0.038,0.091"]
}
')
plot(g1)
exposures(g1) <- "Eczema"
outcomes(g1) <- "anxiety or depression"
adjustmentSets(g1)

## marking adjusted/unobserved variables was done by hand on dagitty.net
