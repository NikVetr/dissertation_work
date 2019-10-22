

 par(mfrow = c(1, 1), mar = c(2, 2, 2, 2))
 #
   # simple cylinders
   emptyplot(c(-1.2, 1.2), c(-1, 1), main = "filledcylinder")
 col <- c(rev(greycol(n = 20)), greycol(n = 20))
 col2 <- shadepalette("red", "blue", n = 20)
 col3 <- shadepalette("yellow", "black", n = 20)
 filledcylinder(rx = 0., ry = 0.2, len = 0.25, angle = 0,
                  col = col, mid = c(-1, 0), dr = 0.1)
 filledcylinder(rx = 0.0, ry = 0.2, angle = 90, col = col,
                  mid = c(-0.5, 0), dr = 0.1)
 
 filledcylinder(rx = 0.1, ry = 0.3, angle = 90, col = c(col2, rev(col2)),
                  mid = c(0.45, 0), topcol = col2[10], dr = 0.1)
 
 filledcylinder(rx = 0.05, ry = 0.2, angle = 90, col = c(col3, rev(col3)),
                  mid = c(0.9, 0), topcol = col3[10], dr = 0.1)
 
 filledcylinder(rx = 0.06, ry = 0.1, angle = 90, col = "white",
                  lcol = "black", lcolint = "grey", dr = 0.2)
 #
   # more complex cylinders
   emptyplot(c(-1, 1), c(-1, 1), main = "filledcylinder")
 col <- shadepalette("blue", "black", n = 20)
 col2 <- shadepalette("red", "black", n = 20)
 col3 <- shadepalette("yellow", "black", n = 20)
 filledcylinder(rx = 0.025, ry = 0.2, angle = 90,
                  col = c(col2, rev(col2)), dr = 0.1, mid = c(-0.8, 0),
                  topcol = col2[10], delt = -1., lcol = "black")
 filledcylinder(rx = 0.1, ry = 0.2, angle = 00,
                  col = c(col, rev(col)), dr = 0.1, mid = c(0.0, 0.0),
                  topcol = col, delt = -1.2, lcol = "black")
 filledcylinder(rx = 0.075, ry = 0.2, angle = 90,
                  col = c(col3, rev(col3)), dr = 0.1, mid = c(0.8, 0),
                  topcol = col3[10], delt = 0.0, lcol = "black")
 #
   # rectangles
   color <- shadepalette(grey(0.3), "blue", n = 20)
   emptyplot(c(-1, 1), main = "filledrectangle")
   filledrectangle(wx = 0.5, wy = 0.5, col = color,
                     mid = c(0, 0), angle = 0)
   filledrectangle(wx = 0.5, wy = 0.5, col = color,
                     mid = c(0.5, 0.5), angle = 90)
   filledrectangle(wx = 0.5, wy = 0.5, col = color,
                     mid = c(-0.5, -0.5), angle = -90)
   filledrectangle(wx = 0.5, wy = 0.5, col = color,
                    Karline Soetaert 7
                     mid = c(0.5, -0.5), angle = 180)
   filledrectangle(wx = 0.5, wy = 0.5, col = color,
                     mid = c(-0.5, 0.5), angle = 270)
   #
     # multigonal
     color <- shadepalette(grey(0.3), "blue", n = 20)
     emptyplot(c(-1, 1))
     filledmultigonal(rx = 0.25, ry = 0.25,
                        col = shadepalette(grey(0.3), "blue", n = 20),
                        nr = 3, mid = c(0, 0), angle = 0)
     filledmultigonal(rx = 0.25, ry = 0.25,
                        col = shadepalette(grey(0.3), "darkgreen", n = 20),
                        nr = 4, mid = c(0.5, 0.5), angle = 90)
     filledmultigonal(rx = 0.25, ry = 0.25,
                        col = shadepalette(grey(0.3), "orange", n = 20),
                        nr = 5, mid = c(-0.5, -0.5), angle = -90)
     filledmultigonal(rx = 0.25, ry = 0.25, col = "black",
                        nr = 6, mid = c(0.5, -0.5), angle = 180)
     filledmultigonal(rx = 0.25, ry = 0.3, col = "white", lcol = "black",
                        nr = 7, mid = c(-0.5, 0.5), angle = 270)
     title("filledmultigonal")
    
      
      
     
     
      library(TeachingDemos) 
      
        ms.cylinder <- function(height, width, angle) { 
         theta <- seq(2*pi, 0, length=200) 
         x1 <- cos(theta) * width 
         y1 <- sin(theta) * width * angle 
          
           x <- c(x1, x1[1:100], x1[100]) 
           y <- c(y1 + height/2, y1[1:100]-height/2, y1[100]+height/2) 
            
          return(cbind(x,y)) 
           }
     
      
        x <- 0:8 %/% 3 
      y <- 0:8 %% 3 
      
      
        h <- runif(9) 
      w <- runif(9) 
      a <- runif(9) 
      
        
      my.symbols(x,y,ms.cylinder, height=h, width=w, angle=a)
      

      
      
plot(x,y, xlim=c(-.5,2.5), ylim=c(-.5,2.5), type='n') 
x <- 1
y <- 1
nsteps <- 100
h <- c(4, rnorm(n = nsteps, mean = 0, sd = 0.1))
h <- cumsum(h)
w <- c(4, rnorm(n = nsteps, mean = 0, sd = 0.1))
w <- cumsum(w)
for(i in 1:length(h)){
  plot(x,y, xlim=c(-.5,2.5), ylim=c(-.5,2.5), type='n') 
  my.symbols(x,y,ms.cylinder, height=h[i], width=w[i], angle=.5)
}
