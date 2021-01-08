if(!require(readxl)) install.packages("readxl")

altPaths <- read.csv("C:/Users/Guillem/Documents/PhD/comput/data/mtbcModelingRes/iEK1011_altFluxes4escher.csv")
altPathsComp <- read.csv("C:/Users/Guillem/Documents/PhD/comput/data/mtbcModelingRes/iEK1011_altFluxes.csv")
getLog2FC <- function(x){
        a <- as.numeric(x[2])
        b <- as.numeric(x[3])
        if(a <= 0 | b <= 0){
                m <- abs(min(c(a, b)))
                a <- a + m
                b <- b + m
        }
        r <- (b/a + 1*b)/(a/a + 1*a)
        log2FC <- log2(b) - log2(a)
        return(log2FC)
}

a <- as.numeric(altPaths[13, 2])
b <- as.numeric(altPaths[13, 3])
if(a <= 0 | b <= 0){
        m <- abs(min(c(a, b)))
        a <- a + m
        b <- b + m
}
r <- (b/a + 1*b)/(a/a + 1*a)
r <- b/a + b*a/(a^2)
a
b
log2FC <- log2(r)
log2FC


2 + 0.1
3 + 3*3

(2/3 + 1*2)/(3/3 + 1*3)

2/3 + 3*2/(3^2)

(2 + 3*2)/(3 + 3^2)
altPaths$log2FC <- apply(altPaths, 
                         1, 
                         getLog2FC)

altPaths$diff <- apply(altPaths, 
                       1, 
                       function(x) as.numeric(x[3]) - as.numeric(x[2]))
max(altPaths$diff)
min(altPaths$diff)


mostAlt <- head(altPaths[order(altPaths$diff), ], 20)
altPathsComp[match(mostAlt$ID, altPathsComp$id), c("id", "name", "subsystem")]
