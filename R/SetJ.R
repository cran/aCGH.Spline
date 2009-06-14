`SetJ` <-
function() { 
library(rJava)
Jpath = system.file("java", package="aCGH.Spline")
.jinit(classpath=Jpath, parameters="-Xmx1600m")
#.jinit(classpath="/Users/tf2/Msc_Thesis/R_JPack/my/", parameters="-Xmx512m")
#.jinit(classpath="/Users/tf2/Msc_Thesis/R_JPack/my/", parameters="-Xmx1500m")
#.jinit(classpath="/Users/tf2/Msc_Thesis/R_JPack/final_package/", parameters="-Xmx1600m")
}

