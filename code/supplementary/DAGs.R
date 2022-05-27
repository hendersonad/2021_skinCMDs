dag {
	"Age/sex/calendar period" [pos="-0.360,-1.252"]
	"Drug abuse" [pos="-0.640,0.947"]
	"Lifestyle factors" [pos="-0.757,0.125"]
	"Psoriasis TRT" [pos="-0.156,1.282"]
	CMDs [outcome,pos="1.222,-0.172"]
	Comorbidities [pos="1.092,0.668"]
	Deprivation [pos="-0.552,-0.446"]
	Ethnicity [pos="-0.665,-0.802"]
	Genetics [latent,pos="-1.746,-1.433"]
	Inflammation [pos="-1.976,1.090"]
	Psoriasis [exposure,pos="-1.838,-0.133"]
	Stigma [pos="-0.849,1.331"]
	Stress [pos="-0.319,-1.614"]
	"Age/sex/calendar period" -> CMDs
	"Age/sex/calendar period" -> Psoriasis
	"Drug abuse" -> CMDs
	"Lifestyle factors" -> CMDs
	"Psoriasis TRT" -> CMDs
	Comorbidities -> CMDs
	Deprivation -> CMDs
	Deprivation -> Psoriasis
	Ethnicity -> CMDs
	Ethnicity -> Psoriasis
	Genetics -> CMDs
	Genetics -> Psoriasis
	Inflammation -> CMDs
	Psoriasis -> "Drug abuse"
	Psoriasis -> "Lifestyle factors"
	Psoriasis -> "Psoriasis TRT"
	Psoriasis -> CMDs
	Psoriasis -> Comorbidities
	Psoriasis -> Inflammation
	Psoriasis -> Stigma
	Stigma -> CMDs
	Stress -> CMDs
	Stress -> Psoriasis
}


dag {
	"Age/sex/calendar period" [adjusted,pos="-0.488,-1.302"]
	"Psoriasis/eczema" [exposure,pos="-1.324,-0.203"]
	CMDs [outcome,pos="0.951,-0.244"]
	Deprivation [adjusted,pos="-0.503,-1.511"]
	Ethnicity [adjusted,pos="-0.488,-1.077"]
	Genetics [adjusted,pos="-2.091,-1.687"]
	alcohol [pos="-0.296,0.019"]
	bmi [pos="-0.274,0.198"]
	comorbidity [adjusted,pos="-0.451,-0.856"]
	sleep [pos="-0.231,0.585"]
	smoking [pos="-0.248,0.370"]
	steroids [pos="-0.205,0.794"]
	"Age/sex/calendar period" -> "Psoriasis/eczema"
	"Age/sex/calendar period" -> CMDs
	"Psoriasis/eczema" -> CMDs
	"Psoriasis/eczema" -> alcohol
	"Psoriasis/eczema" -> bmi
	"Psoriasis/eczema" -> sleep
	"Psoriasis/eczema" -> smoking
	"Psoriasis/eczema" -> steroids
	Deprivation -> "Psoriasis/eczema"
	Deprivation -> CMDs
	Ethnicity -> "Psoriasis/eczema"
	Ethnicity -> CMDs
	Genetics -> "Psoriasis/eczema"
	Genetics -> CMDs
	alcohol -> CMDs
	bmi -> CMDs
	comorbidity -> "Psoriasis/eczema"
	comorbidity -> CMDs
	sleep -> CMDs
	smoking -> CMDs
	steroids -> CMDs
}

# v2 ----------------------------------------------------------------------


