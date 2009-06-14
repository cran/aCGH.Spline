`SetJ` <-
function(JavaHeapSize="-Xmx1000m" ) { 
library(rJava)
Jpath = system.file("java", package="aCGH.Spline")
.jinit(classpath=Jpath, parameters=JavaHeapSize)
}

