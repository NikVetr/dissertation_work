library(latex2exp)
library(pBrackets)
rootStates <- rbind(c(30,100), c(100,70), c(30,1), c(100,1), c(150,1))
rootStates <- rbind(c(20,5), c(70,5), c(100,5), c(170,5), c(210,5))
plot(0:270, 0:270*.5, col = 0, axes = F, xlab = "", ylab = "")

text2 <- function(label, coordinates, cex, font) {text(label = label, x = coordinates[1], y = coordinates[2], cex = cex, font = font)}
par(mar = c(1,1,1,1))
ang <- 60
angh <- ang / 2
angr <- angh / 180 * pi
liwd <- 5
#tree1
v1 <- 50
v2 <- 30
v3 <- 35
v4 <- 40
E <- rootStates[1,]
A <- E + c(-v1 * sin(angr), v1 * cos(angr))
D <- E + c(v4 * sin(angr), v4 * cos(angr))
B <- D + c(-v2 * sin(angr), v2 * cos(angr))
C <- D + c(v3 * sin(angr), v3 * cos(angr))

lines(rbind(E, A), lwd = liwd)
lines(rbind(E, D), lwd = liwd)
lines(rbind(D, B), lwd = liwd)
lines(rbind(D, C), lwd = liwd)

text2("A", A + c(0, 5), cex = 1.3, font = 4)
text2("B", B + c(0, 5), cex = 1.3, font = 4)
text2("C", C + c(0, 5), cex = 1.3, font = 4)
text2("D", D + c(3, -3), cex = 1.3, font = 4)
text2("E", E + c(0, -5), cex = 1.3, font = 4)

text2(TeX("v_1"), E + c(-v1 * sin(angr), v1 * cos(angr)) / 2 + c(-4, -3), cex = 1.1, font = 3)
text2(TeX("v_2"), D + c(-v2 * sin(angr), v2 * cos(angr)) / 2 + c(-4, -3), cex = 1.1, font = 3)
text2(TeX("v_3"), D + c(v3 * sin(angr), v3 * cos(angr)) / 2 + c(4, -3), cex = 1.1, font = 3)
text2(TeX("v_4"), E + c(v1 * sin(angr), v1 * cos(angr)) / 2 + c(4, -3), cex = 1.1, font = 3)

#tree2, reroot on A

# plot(0:200, 0:200, col = 0, axes = F, xlab = "", ylab = "")
A <- rootStates[2,]
E <- A + c(0, v1)
D <- E + c(0, v4)
B <- D + c(-v2 * sin(angr), v2 * cos(angr))
C <- D + c(v3 * sin(angr), v3 * cos(angr))

lines(rbind(A, E), lwd = liwd)
lines(rbind(E, D), lwd = liwd)
lines(rbind(D, B), lwd = liwd)
lines(rbind(D, C), lwd = liwd)

text2("A", A + c(0, -5), cex = 1.3, font = 4)
text2("B", B + c(0, 5), cex = 1.3, font = 4)
text2("C", C + c(0, 5), cex = 1.3, font = 4)
text2("D", D + c(4, -2), cex = 1.3, font = 4)
text2("E", E + c(4, 0), cex = 1.3, font = 4)

text2(TeX("v_1"), A + c(0, v1) / 2 + c(-4, -3), cex = 1.1, font = 3)
text2(TeX("v_2"), D + c(-v2 * sin(angr), v2 * cos(angr)) / 2 + c(-4, -3), cex = 1.1, font = 3)
text2(TeX("v_3"), D + c(v3 * sin(angr), v3 * cos(angr)) / 2 + c(4, -3), cex = 1.1, font = 3)
text2(TeX("v_4"), E + c(0, v4) / 2 + c(-4, -3), cex = 1.1, font = 3)

#tree2, reroot on B

# plot(0:200, 0:200, col = 0, axes = F, xlab = "", ylab = "")
B <- rootStates[3,]
D <- B + c(0, v2)
E <- D + c(v4 * sin(angr), v4 * cos(angr))
A <- E + c(v1 * sin(angr), v1 * cos(angr))
C <- D + c(-v3 * sin(angr), v3 * cos(angr))

lines(rbind(A, E), lwd = liwd)
lines(rbind(E, D), lwd = liwd)
lines(rbind(D, B), lwd = liwd)
lines(rbind(D, C), lwd = liwd)

text2("A", A + c(0, 5), cex = 1.3, font = 4)
text2("B", B + c(0, -5), cex = 1.3, font = 4)
text2("C", C + c(0, 5), cex = 1.3, font = 4)
text2("D", D + c(4, -2), cex = 1.3, font = 4)
text2("E", E + c(4, 0), cex = 1.3, font = 4)

text2(TeX("v_1"), E + c(v1 * sin(angr), v1 * cos(angr)) / 2 + c(4, -4), cex = 1.1, font = 3)
text2(TeX("v_2"), B + c(0, v2) / 2 + c(-5, 0), cex = 1.1, font = 3)
text2(TeX("v_3"), D + c(-v3 * sin(angr), v3 * cos(angr)) / 2 + c(-4, -4), cex = 1.1, font = 3)
text2(TeX("v_4"), D + c(v4 * sin(angr), v4 * cos(angr)) / 2 + c(4, -4), cex = 1.1, font = 3)

#tree4, reroot along v4

v1 <- 50
v2 <- 30
v3 <- 35
v4 <- 40
subBL <- v4 / 2

R <- rootStates[4,]
E <- R + c(-subBL * sin(angr), subBL * cos(angr))
A <- E + c(-v1 * sin(angr), v1 * cos(angr))
D <- R + c((v4-subBL) * sin(angr), (v4-subBL) * cos(angr))
B <- D + c(-v2 * sin(angr), v2 * cos(angr))
C <- D + c(v3 * sin(angr), v3 * cos(angr))

# plot(0:200, 0:200, col = 0, axes = F, xlab = "", ylab = "")
lines(rbind(R, E), lwd = liwd)
lines(rbind(E, A), lwd = liwd)
lines(rbind(R, D), lwd = liwd)
lines(rbind(D, B), lwd = liwd)
lines(rbind(D, C), lwd = liwd)

text2("A", A + c(0, 5), cex = 1.3, font = 4)
text2("B", B + c(0, 5), cex = 1.3, font = 4)
text2("C", C + c(0, 5), cex = 1.3, font = 4)
text2("D", D + c(3, -3), cex = 1.3, font = 4)
text2("E", E + c(-5, 0), cex = 1.3, font = 4)

text2(TeX("v_1"), E + c(-v1 * sin(angr), v1 * cos(angr)) / 2 + c(-4, -3), cex = 1.1, font = 3)
text2(TeX("v_2"), D + c(-v2 * sin(angr), v2 * cos(angr)) / 2 + c(-4, -3), cex = 1.1, font = 3)
text2(TeX("v_3"), D + c(v3 * sin(angr), v3 * cos(angr)) / 2 + c(4, -3), cex = 1.1, font = 3)
text2(TeX("r"), R + c(-subBL * sin(angr), subBL * cos(angr)) / 2 + c(-4, -1), cex = 1.1, font = 3)
text2(TeX("v_4-r"), R + c((v4-subBL) * sin(angr), (v4-subBL) * cos(angr)) / 2 + c(5, -3), cex = 1.1, font = 3)

#tree5, reroot along v3

v1 <- 50
v2 <- 30
v3 <- 35
v4 <- 40
subBL <- v3 / 2

R <- rootStates[5,]
D <- R + c((v3-subBL) * sin(angr), (v3-subBL) * cos(angr))
C <- R + c(-subBL * sin(angr), subBL * cos(angr))
E <- D + c(v4 * sin(angr), v4 * cos(angr))
B <- D + c(-v2 * sin(angr), v2 * cos(angr))
A <- E + c(v1 * sin(angr), v1 * cos(angr))

# plot(0:200, 0:200, col = 0, axes = F, xlab = "", ylab = "")
lines(rbind(R, C), lwd = liwd)
lines(rbind(R, D), lwd = liwd)
lines(rbind(D, E), lwd = liwd)
lines(rbind(D, B), lwd = liwd)
lines(rbind(E, A), lwd = liwd)

text2("A", A + c(0, 5), cex = 1.3, font = 4)
text2("B", B + c(0, 5), cex = 1.3, font = 4)
text2("C", C + c(0, 5), cex = 1.3, font = 4)
text2("D", D + c(3, -3), cex = 1.3, font = 4)
text2("E", E + c(-5, 0), cex = 1.3, font = 4)

text2(TeX("v_1"), E + c(v1 * sin(angr), v1 * cos(angr)) / 2 + c(-4, 3), cex = 1.1, font = 3)
text2(TeX("v_2"), D + c(-v2 * sin(angr), v2 * cos(angr)) / 2 + c(4, 3), cex = 1.1, font = 3)
text2(TeX("v_4"), D + c(v3 * sin(angr), v3 * cos(angr)) / 2 + c(4, -3), cex = 1.1, font = 3)
text2(TeX("r"), R + c(-subBL * sin(angr), subBL * cos(angr)) / 2 + c(-4, -1), cex = 1.1, font = 3)
text2(TeX("v_3-r"), R + c((v4-subBL) * sin(angr), (v4-subBL) * cos(angr)) / 2 + c(5, -3), cex = 1.1, font = 3)


## let's try to show the pruning algorithm ##
rootStates <- rbind(c(20,5), c(110,5), c(200,5))
plot(0:270, 0:270*.5, col = 0, axes = F, xlab = "", ylab = "")

par(mar = c(1,1,1,1))
ang <- 60
angh <- ang / 2
angr <- angh / 180 * pi
liwd <- 5
#tree1
v1 <- 50
v2 <- 30
v3 <- 35
v4 <- 40
E <- rootStates[1,]
A <- E + c(-v1 * sin(angr), v1 * cos(angr))
D <- E + c(v4 * sin(angr), v4 * cos(angr))
B <- D + c(-v2 * sin(angr), v2 * cos(angr))
C <- D + c(v3 * sin(angr), v3 * cos(angr))

lines(rbind(E, A), lwd = liwd)
lines(rbind(E, D), lwd = liwd)
lines(rbind(D, B), lwd = liwd)
lines(rbind(D, C), lwd = liwd)

text2("A", A + c(0, 5), cex = 1.3, font = 4)
text2("B", B + c(0, 5), cex = 1.3, font = 4)
text2("C", C + c(0, 5), cex = 1.3, font = 4)
text2("D", D + c(3, -3), cex = 1.3, font = 4)
text2("E", E + c(0, -5), cex = 1.3, font = 4)

text2(TeX("v_1"), E + c(-v1 * sin(angr), v1 * cos(angr)) / 2 + c(-4, -3), cex = 1.1, font = 3)
text2(TeX("v_2"), D + c(-v2 * sin(angr), v2 * cos(angr)) / 2 + c(-4, -3), cex = 1.1, font = 3)
text2(TeX("v_3"), D + c(v3 * sin(angr), v3 * cos(angr)) / 2 + c(4, -3), cex = 1.1, font = 3)
text2(TeX("v_4"), E + c(v1 * sin(angr), v1 * cos(angr)) / 2 + c(4, -3), cex = 1.1, font = 3)

arrows(x0 = 60, y0 = 35, x1 = 75, y1 = 35, lwd = 3)

#tree2
v1 <- 50
v2 <- 30
v3 <- 35
v4 <- 40
v4c <- v4 + v2*v3/(v2+v3)
E <- rootStates[2,]
A <- E + c(-v1 * sin(angr), v1 * cos(angr))
D <- E + c(v4 * sin(angr), v4 * cos(angr))
Dp <- E + c(v4c * sin(angr), v4c * cos(angr))

lines(rbind(E, A), lwd = liwd)
lines(rbind(E, D), lwd = liwd)
lines(rbind(D, Dp), lwd = liwd/2.4, lty = 3)

text2("A", A + c(0, 5), cex = 1.3, font = 4)
text2("D'", Dp + c(0, 5), cex = 1.3, font = 4)
text2("E", E + c(0, -5), cex = 1.3, font = 4)

text2(TeX("v_1"), E + c(-v1 * sin(angr), v1 * cos(angr)) / 2 + c(-4, -3), cex = 1.1, font = 3)
text2(TeX("v_4"), E + c(v1 * sin(angr), v1 * cos(angr)) / 2 + c(4, -3), cex = 1.1, font = 3)
text2(TeX("$\\frac{v_2v_3}{v_2+v_3}"), D + c(v2*v3/(v2+v3) * sin(angr), v2*v3/(v2+v3) * cos(angr)) / 2 + c(-10, 3), cex = 1.1, font = 3)

lines(rbind(c(140, 35),c(148,35)), lwd = 2); lines(rbind(c(144, 31),c(144,39)), lwd = 2) 

lines(rbind(c(160, 5),c(160, v2+v3+5)), lwd = liwd)
text2("C", c(160, v2+v3+5) + c(0, 5), cex = 1.3, font = 4)
text2("B", c(160, 5) + c(0, -5), cex = 1.3, font = 4)
text2(TeX("v_2+v_3"), c(160, 5) + c(0, v2+v3) / 2 + c(9, 0), cex = 1.1, font = 3)

arrows(x0 = 180, y0 = 35, x1 = 195, y1 = 35, lwd = 3)

#full decomposition

lines(rbind(c(140, 35),c(148,35)), lwd = 2); lines(rbind(c(144, 31),c(144,39)), lwd = 2) 

lines(rbind(c(225, 5),c(225, 5+v1+v4+v2*v3/(v2+v3))), lwd = liwd)
text2("D'", c(225, 5+v1+v4+v2*v3/(v2+v3)) + c(0, 5), cex = 1.3, font = 4)
text2("A", c(225, 5) + c(0, -5), cex = 1.3, font = 4)
text2(TeX("v_1+v_4+$\\frac{v_2v_3}{v_2+v_3}"), c(225, 5) + c(0, 5+v1+v4+v2*v3/(v2+v3)) / 2 + c(-17, 0), cex = 1.1, font = 3)

lines(rbind(c(250, 5),c(250, v2+v3)), lwd = liwd)
text2("C", c(250, v2+v3) + c(0, 5), cex = 1.3, font = 4)
text2("B", c(250, 5) + c(0, -5), cex = 1.3, font = 4)
text2(TeX("v_2+v_3"), c(250, 5) + c(0, v2+v3) / 2 + c(9, 0), cex = 1.1, font = 3)

lines(rbind(c(234, 35),c(242,35)), lwd = 2); lines(rbind(c(238, 31),c(238,39)), lwd = 3) 

plotMatrix <- function(mobject, size, location, lwd = 2, grid = T, font = 1, cex = 1, rownames = T, colnames = T, title = T, title.label = "Matrix Object"){
  lines(rbind(location, location + c(0,size[2])), lwd = lwd)
  lines(rbind(location, location + c(size[1]/8,0)), lwd = lwd)
  lines(rbind(location + c(0, size[2]), location + c(size[1]/8,size[2])), lwd = lwd)
  lines(rbind(location + c(size[1],0), location + size), lwd = lwd)
  lines(rbind(location + size, location + size - c(size[1]/8,0)), lwd = lwd)
  lines(rbind(location + c(size[1],0), location + c(size[1],0) - c(size[1]/8,0)), lwd = lwd)
  if(grid == T){
    for(i in 1:(dim(mobject)[1]-1)){
      lines(rbind(location + c(0,i*size[2]/dim(mobject)[1]), location + c(size[1], i*size[2]/dim(mobject)[1])))
    }
    for(j in 1:(dim(mobject)[2]-1)){
      lines(rbind(location + c(j*size[1]/dim(mobject)[2],0), location + c(j*size[1]/dim(mobject)[2], size[2])))
    }
  }
  if(class(mobject[1,1]) != "expression" & class(mobject[1,1]) != "character"){mobject <- matrix(as.character(mobject), nrow = dim(mobject)[1], ncol = dim(mobject)[2])}
  for(i in 1:(dim(mobject)[1])){
    for(j in 1:dim(mobject)[2]){
      text(labels = mobject[i,j], x = location[1] + (j-1/2)*size[1]/dim(mobject)[2], y = location[2] + size[2] - (i-1/2)*size[2]/dim(mobject)[1], font = font, cex = cex)
    }
  }
  if(title){
    text(title.label, x = location[1] + size[1]/2, y = location[2] + size[2] + strheight(title.label, font = 2, cex = 1.5)/1.5, cex = 1.5, font = 2)
  }
  if(rownames){
    for(i in 1:dim(mobject)[1]){
      text(rownames(mobject)[i], x = location[1] - strwidth(rownames(mobject)[i])/2 - size[1]/(ncol(mobject)*6), y = location[2] + size[2] - (i-1/2)*size[2]/dim(mobject)[2])
    }
  }
  if(colnames){
    for(i in 1:dim(mobject)[1]){
      text(colnames(mobject)[i], x = location[1] + (i-1/2)*size[1]/dim(mobject)[1], y = location[2] - strheight(colnames(mobject)[i])/2- size[2]/(nrow(mobject)*6))
    }
  }
}

Q <- matrix("", nrow = 2, ncol = 2)
plotMatrix(mobject = Q, size = c(25,20), location = c(1,110), title.label = "Q")
y <- matrix("", nrow = 2, ncol = 1)
plotMatrix(mobject = y, size = c(25,20), location = c(31,110), title.label = "y")
P <- matrix(c(TeX("v_1"), 0, 0, 0, TeX("v_4+v_2"), TeX("v_4"), 0, TeX("v_4"), TeX("v_4+v_3")), nrow = 3, ncol = 3); rownames(P) <- colnames(P) <- LETTERS[1:3]
plotMatrix(mobject = P, size = c(40,20), location = c(1,80), title.label = "P")
x <- matrix(c(TeX("x_A"),TeX("x_B"),TeX("x_C")), nrow = 3, ncol = 1)
plotMatrix(mobject = x, size = c(10,20), location = c(46,80), title.label = "x")

arrows(x0 = 60, y0 = 105, x1 = 75, y1 = 105, lwd = 3)

Q <- matrix(c(TeX("v_2+v_3"), 0, 0, ""), nrow = 2, ncol = 2)
plotMatrix(mobject = Q, size = c(30,20), location = c(80,110), title.label = "Q")
y <- matrix(c(TeX("x_B-x_C"), ""), nrow = 2, ncol = 1)
plotMatrix(mobject = y, size = c(20,20), location = c(116,110), title.label = "y")
P <- matrix(c(TeX("v_1"), 0, 0, 0, TeX("v_4+v_2"), TeX("v_4"), 0, TeX("v_4"), TeX("v_4+v_3")), nrow = 3, ncol = 3); rownames(P) <- colnames(P) <- LETTERS[1:3]
plotMatrix(mobject = P, size = c(40,20), location = c(80,80), title.label = "P")
x <- matrix(c(TeX("x_A"),TeX("x_B"),TeX("x_C")), nrow = 3, ncol = 1)
plotMatrix(mobject = x, size = c(10,20), location = c(126,80), title.label = "x")

arrows(x0 = 140, y0 = 105, x1 = 155, y1 = 105, lwd = 3)

Q <- matrix(c(TeX("v_2+v_3"), 0, 0, TeX("v_1+v_4+$\\frac{v_2v_3}{v_2+v_3}")), nrow = 2, ncol = 2)
plotMatrix(mobject = Q, size = c(60,50), location = c(160,80), title.label = "Q")
y <- matrix(c(TeX("x_B-x_C"), TeX("x_a - $\\frac{v_3x_b+v_2x_c}{v_2+v_3}")), nrow = 2, ncol = 1)
plotMatrix(mobject = y, size = c(30,50), location = c(230,80), title.label = "y")
