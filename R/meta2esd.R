# author Rasmus E. Benestad
# Last-Update 2013-07-23 by Abdelkader Mezghani
 
decryptcn <- function(codes,src="ECAD") {
  n <- length(codes)
  country <- rep("",n)
  for (i in 1:n) {
    country[i]= switch(codes[i],'SE'='Sweden','AT'=,'Austria','BE'='Belgium',
            'HR'='Croatia','CY'='Cypros','CZ'='Chzeck Republic',
            'FI'='Finland','FR'='France','DE'='Germany',
            'IS'='Iceland','RU'='Russia','DK'='Denmark',
            'IE'='Ireland','NL'='the Netherlands','IT'='Italy',
            'NO'='Norway','LV'='Latvia','LT'='Lituania',
            'PT'='Portugal','RO'='Romania','SK'='Slovakia',
            'SI'='Slovenia','ES'='Spain','CH'='Switzerland',
            'RS'='Serbia','EE'='Estonia','MK'='Makedonia',
            'GB'='Great Britain','BA'='Bosnia-Hertsogovina',
            'AL'='Albania','DZ'='Algeria','LU'='Luxemburg',
            'AM'='Armenia','GL'='Greenland','AZ'='Azerbaijan',
            'EG'='Egypt','GR'='Greece','PL'='Poland','IL'='Israel',
            'BY'='Belarus','GE'='Georgia','HU'='Hungary',
            'IQ'='Iraq','KZ'='Khazakstan','LY'='Libya',
            'MD'='Moldova','MO'='Morocco','MT'='Malta',
            'SA'='Saudi Arabia','SY'='Syria','TJ'='Tajikistan',
            'TR'="Turkey",'UA'='Ukraina','UZ'='Uzbekistan',
            'B'='Belgia','FIN'='Finland','FRI'='Faroe Islands',
            'G'='Greenland','IRL'='Ireland','IS'='Iceland',
            'N'='Norway','S'='Sweden')
  }
  invisible(country)
}


