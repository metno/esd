ghcn.country.code <- function(code = NULL , name = NULL) {

if (!is.null(code)) icode <- as.numeric(code) else icode <- NULL
if (!is.null(name)) iname <- tolower(as.character(name)) else iname <- NULL

x <- read.csv2("~/SHARED/data/country-codes",head=FALSE)
x.code <- as.numeric(substr(x$V1,1,3))
x.name <- substr(x$V1,5,nchar(as.character(x$V1)))
x.name <- as.character(gsub("  ",x=x.name,replacement=""))
nx <- length(x.name)

# remove extra characters such as empty spaces
id <- (substr(x.name,nchar(x.name),nchar(x.name))==" ")
x.name[id] <- substr(x.name[id],1,nchar(x.name[id])-1)

# create a data frame
df <- data.frame(code = x.code , name = toupper(x.name))

# sub selection based on specified code or name
if (!is.null(icode)) df <- subset(df, code==icode)
if (!is.null(iname)) df <- subset(df, name==iname) 

# Upper case name
df$name <- toupper(as.character(df$name))
return(df)
}

ghcnd.country.abb  <- function(abb = NULL , name = NULL) {

iabb <- abb
iname <- name

x      <- read.csv2("~/SHARED/data/ghcnd/ghcnd-countries.txt",head=FALSE)
x.abb  <- substr(x$V1,1,2)
x.name <- substr(x$V1,4,nchar(as.character(x$V1)))

# remove extra characters such as empty spaces
id <- (substr(x.name,nchar(x.name),nchar(x.name))==" ")
x.name[id] <- substr(x.name[id],1,nchar(x.name[id])-1)

# create a data frame
df <- data.frame(abb = x.abb , name = tolower(x.name))

# sub selection based on specified code or name
if (!is.null(iabb)) {
   id <- is.element(df$abb,iabb)
   df <- subset(df, id)
}
if (!is.null(iname)) {
   id <- is.element(df$name,iname)
   df <- subset(df, name==iname) 
}
df$abb <- as.character(df$abb)
df$name <- toupper(as.character(df$name))
return(df)
}

ecad.country.abb  <- function(abb = NULL , name = NULL) {

iabb <- abb
iname <- name

x      <- read.csv2("~/SHARED/data/ecad/ecad-countries.txt",head=FALSE)
x.abb  <- substr(x$V1,1,2)
x.name <- substr(x$V1,4,nchar(as.character(x$V1)))

# remove extra characters such as empty spaces
id <- (substr(x.name,nchar(x.name),nchar(x.name))==" ")
x.name[id] <- substr(x.name[id],1,nchar(x.name[id])-1)

# create a data frame
df <- data.frame(abb = x.abb , name = tolower(x.name))

# sub selection based on specified code or name
if (!is.null(iabb)) {
   id <- is.element(df$abb,iabb)
   df <- subset(df, id)
}
if (!is.null(iname)) {
   id <- is.element(df$name,iname)
   df <- subset(df, name==iname) 
}
df$abb <- as.character(df$abb)
df$name <- toupper(as.character(df$name))
return(df)
}


