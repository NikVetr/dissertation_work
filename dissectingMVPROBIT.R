library(Matrix)
library(matrixcalc)
library(adephylo)
library(phytools)
library(phytools)
library(mvMORPH)
library(mvtnorm)
library(mnormt)
library(ggplot2)
library(phangorn)
library(distory)
library(psych)
library(abind)
library(sirt)
library(irtoys)
library(mvProbit)
library(miscTools)
library(bayesm)

tree <- pbtree(n=5)
plot(tree)
#plot(tree); axisPhylo(); tiplabels(); nodelabels(); edgelabels()
nTraits <- 5
traitsToUse <- sample(linMeasNames, nTraits)
sig <- cov2cor(cov)[traitsToUse, traitsToUse]
sig <- matrix(nrow = 5, ncol = 5, data = c(1, 0.1179964, 0.5630421, 0.1665352, 0.03472295, 0.1179964, 1, 0.07666353, 0.752195, 0.2024577, 0.5630421, 0.07666353, 1, 0.1907083, 0.01605284, 0.1665352, 0.752195, 0.1907083, 1, 0.1626185, 0.03472295, 0.2024577, 0.01605284, 0.1626185, 1))
meanLiabs <- mvSIM(tree = tree, nsim = 1, model = "BM1", param = list(ntraits = nTraits, sigma=sig, mu= rep(0,nTraits)))
populations <- tree$tip.label
numTips <- length(populations)
nObs <- 20

#construct design matrix
xMat <- matrix(data = 0, ncol = numTips, nrow = nObs * numTips)
for (i in 1:numTips){xMat[((nObs*(i-1)+1):(nObs*i)),i] <- rep(1, nObs)}
colnames(xMat) <- paste0("x", 1:numTips)
xMat <- cbind(0, xMat)

# model coefficients
beta <- meanLiabs
beta <- rbind(0, beta)

# covariance matrix of error terms
sigma <- sig

# generate dependent variables
yMatLin <- xMat %*% beta
yMat <- ( yMatLin + rmvnorm( nObs*numTips, sigma = sigma ) ) > 0
colnames( yMat ) <- paste0( "y", 1:nTraits)

# log likelihood values
myData <- as.data.frame( cbind( xMat, yMat ) )
formulaMvProbit <- as.formula(paste0("cbind( ", paste0("y", 1:(nTraits-1), sep = ",", collapse = ""), paste0("y", nTraits), 
                                     ") ~ ", paste0("x", 1:(numTips-1), sep = "+", collapse = ""), paste0("x", numTips)))
formulaMvProbit

logLikValTrue <- mvProbitLogLik( formula = formulaMvProbit, coef = c( beta ), sigma = sigma, data = myData, algorithm = "GHK" )

#find density of yMat[1,] in mvn(mean = yMatLin[1,], cov = sigma)


mvProbitLogLik <- function (formula, coef, sigma = NULL, data, algorithm = "GHK", 
                            nGHK = 1000, returnGrad = oneSidedGrad, oneSidedGrad = FALSE, 
                            eps = 1e-06, random.seed = 123, ...) 
    {
      if (is.list(formula)) {
        stop("using different regressors for the dependent variables", 
             " has not been implemented yet. Sorry!")
      }
      else if (class(formula) != "formula") {
        stop("argument 'formula' must be a formula")
      }
      if (!is.data.frame(data)) {
        stop("argument 'data' must be a data frame")
      }
      mc <- match.call(expand.dots = FALSE)
      m <- match("data", names(mc), 0)
      mf <- mc[c(1, m)]
      mf$formula <- formula
      attributes(mf$formula) <- NULL
      mf$na.action <- na.pass
      mf[[1]] <- as.name("model.frame")
      mf <- eval(mf, parent.frame())
      mt <- attr(mf, "terms")
      xMat <- model.matrix(mt, mf)
      yMat <- model.response(mf)
      if (!is.matrix(yMat)) {
        stop("at least two dependent variables", " must be specified in argument 'formula'", 
             " (e.g. by 'cbind( y1, y2 ) ~ ...')")
      }
      else if (!all(yMat %in% c(0, 1, TRUE, FALSE))) {
        stop("all dependent variables must be either 0, 1, TRUE, or FALSE")
      }
      result <- mvProbitLogLikInternal(yMat = yMat, xMat = xMat, 
                                       coef = coef, sigma = sigma, algorithm = algorithm, nGHK = nGHK, 
                                       returnGrad = returnGrad, oneSidedGrad = oneSidedGrad, 
                                       eps = eps, randomSeed = random.seed, ...)
      return(result)
    }


mvProbitLogLikInternal <- function( yMat, xMat, coef, sigma,
                                    algorithm, nGHK, returnGrad, oneSidedGrad, eps, randomSeed, ... ) {
  
  # number of regressors
  nReg <- ncol( xMat )
  
  # checking and preparing model coefficients and correlation coefficients
  coef <- mvProbitPrepareCoef( yMat = yMat, nReg = nReg, coef = coef, 
                               sigma = sigma )
  
  # checking argument 'returnGrad'
  if( length( returnGrad ) != 1 ) {
    stop( "argument 'returnGrad' must be a single logical value" )
  } else if( !is.logical( returnGrad ) ) {
    stop( "argument 'returnGrad' must be logical" )
  }
  
  # checking argument 'oneSidedGrad'
  if( length( oneSidedGrad ) != 1 ) {
    stop( "argument 'oneSidedGrad' must be a single logical value" )
  } else if( !is.logical( oneSidedGrad ) ) {
    stop( "argument 'oneSidedGrad' must be logical" )
  }
  
  # checking argument 'eps'
  if( oneSidedGrad ) {
    if( length( eps ) != 1 ) {
      stop( "argument 'eps' must be a single numeric value" )
    } else if( !is.numeric( eps ) ) {
      stop( "argument 'eps' must be numeric" )
    }
  }
  
  # number of dependent variables
  nDep <- ncol( coef$sigma )
  
  # number of observations
  nObs <- nrow( xMat )
  
  # calculating linear predictors
  xBeta <- matrix( NA, nrow = nObs, ncol = nDep )
  for( i in 1:nDep ) {
    xBeta[ , i ] <- xMat %*% coef$betaEq[[ i ]]
  }
  
  # calculate log likelihood values (for each observation)
  result <- rep( NA, nObs )
  for( i in 1:nObs ){
    ySign <- 2 * yMat[ i, ] - 1
    xBetaTmp <- xBeta[ i, ] * ySign
    sigmaTmp <- diag( ySign ) %*% coef$sigma %*% diag( ySign )
    result[ i ] <- log( pmvnormWrap( upper = xBetaTmp, sigma = sigmaTmp, 
                                     algorithm = algorithm, nGHK = nGHK, random.seed = randomSeed) )
  }
  
  if( returnGrad ) {
    allCoef <- c( coef$beta, coef$sigma[ lower.tri( coef$sigma ) ] )
    grad <- matrix( NA, nrow = length( result ), 
                    ncol = length( allCoef ) )
    for( i in 1:nDep ) {
      # gradients of intercepts
      coefLower <- coefUpper <- allCoef
      InterceptNo <- ( i - 1 ) * nReg + 1
      if( oneSidedGrad ) {
        coefUpper[ InterceptNo ] <- allCoef[ InterceptNo ] + eps
        llLower <- result
      } else {
        coefLower[ InterceptNo ] <- allCoef[ InterceptNo ] - eps / 2
        coefUpper[ InterceptNo ] <- allCoef[ InterceptNo ] + eps / 2
        llLower <- mvProbitLogLikInternal( yMat = yMat, xMat = xMat, 
                                           coef = coefLower, sigma = NULL, 
                                           algorithm = algorithm, nGHK = nGHK,
                                           returnGrad = FALSE, oneSidedGrad = FALSE, eps = 0, 
                                           randomSeed = randomSeed, ... )
      }
      llUpper <- mvProbitLogLikInternal( yMat = yMat, xMat = xMat, 
                                         coef = coefUpper, sigma = NULL, 
                                         algorithm = algorithm, nGHK = nGHK,
                                         returnGrad = FALSE, oneSidedGrad = FALSE, eps = 0, 
                                         randomSeed = randomSeed, ... )
      grad[ , InterceptNo ] <- ( llUpper - llLower ) / eps
      # gradients of coefficients of other explanatory variables
      if( nReg > 1 ) {
        for( j in 2:nReg ) {
          grad[ , InterceptNo + j - 1 ] <- 
            grad[ , InterceptNo ] * xMat[ , j ] 
        }
      }
    }
    # gradients of correlation coefficients
    for( i in ( nDep * nReg + 1 ):length( allCoef ) ) {
      coefLower <- coefUpper <- allCoef
      if( oneSidedGrad ) {
        coefUpper[ i ] <- allCoef[ i ] + eps
        llLower <- result
      } else {
        coefLower[ i ] <- allCoef[ i ] - eps / 2
        coefUpper[ i ] <- allCoef[ i ] + eps / 2
        llLower <- mvProbitLogLikInternal( yMat = yMat, xMat = xMat, 
                                           coef = coefLower, sigma = NULL, 
                                           algorithm = algorithm, nGHK = nGHK,
                                           returnGrad = FALSE, oneSidedGrad = FALSE, eps = 0, 
                                           randomSeed = randomSeed, ... )
      }
      llUpper <- mvProbitLogLikInternal( yMat = yMat, xMat = xMat, 
                                         coef = coefUpper, sigma = NULL, 
                                         algorithm = algorithm, nGHK = nGHK,
                                         returnGrad = FALSE, oneSidedGrad = FALSE, eps = 0, 
                                         randomSeed = randomSeed, ... )
      grad[ , i ] <- ( llUpper - llLower ) / eps
    }
    colnames( grad ) <- mvProbitCoefNames( nDep = nDep, nReg = nReg )
    attr( result, "gradient" ) <- grad
  }
  
  return( result )
}

mvProbitPrepareCoef <- function( yMat, nReg, coef, sigma ) {
  
  # checking argument 'coef'
  if( !is.vector( coef, mode = "numeric" ) ) {
    stop( "argument 'coef' must be a numeric vector" )
  }
  
  if( !is.null( sigma ) ) {
    # checking argument 'sigma'
    if( !is.matrix( sigma ) ) {
      stop( "argument 'sigma' must be a matrix" )
    } else if( nrow( sigma ) != ncol( sigma ) ) {
      stop( "argument 'sigma' must be a quadratic matrix" )
    } else if( !isSymmetric( sigma ) ) {
      stop( "argument 'sigma' must be a symmetric matrix" )
    } else if( any( abs( diag( sigma ) - 1 ) > 1e-7 ) ) {
      stop( "argument 'sigma' must have ones on its diagonal" )
    } else if( !is.null( yMat ) ) {
      if( ncol( sigma ) != ncol( yMat ) ) {
        stop( "the number of dependent variables specified in argument",
              " 'formula' must be equal to the number of rows and colums",
              " of the matrix specified by argument 'sigma'" )
      }
    }
    # number of dependent variables
    nDep <- ncol( sigma )
    # number of model coefficients
    nCoef <- nDep * nReg
    if( length( coef ) != nCoef ) {
      stop( "given that argument 'sigma' has been specified",
            " argument coef must have ", nCoef, " elements" )
    }
  } else {
    if( !is.null( yMat ) ) {
      # number of dependent variables
      nDep <- ncol( yMat )
      # number of model coefficients
      nCoef <- nDep * nReg
      # number of parameters including sigma
      nCoefSigma <- nCoef + nDep * ( nDep - 1 ) / 2
      if( length( coef ) != nCoefSigma ) {
        stop( "given that argument 'sigma' is 'NULL'",
              " argument coef must have ", nCoefSigma, " elements" )
      }
    } else {
      # number of dependent variables
      nDep <- round( - nReg + 0.5 + 
                       sqrt( ( nReg - 0.5 )^2 + 2 * length( coef ) ) )
      # number of model coefficients
      nCoef <- nDep * nReg
      # number of parameters including sigma
      nCoefSigma <- nCoef + nDep * ( nDep - 1 ) / 2
      if( length( coef ) != nCoefSigma ) {
        stop( "given that argument 'sigma' is 'NULL'",
              " argument coef must have ", nCoefSigma, " elements",
              " if the model has ", nDep, " dependent variables" )
      }
    }
    # extracting correlation coefficients from 'coef' if they are there
    sigma <- diag( nDep )
    sigma[ lower.tri( sigma ) ] <- coef[ -c( 1:nCoef ) ]
    sigma[ upper.tri( sigma ) ] <- t( sigma )[ upper.tri( sigma ) ]
    coef <- coef[ 1:nCoef ]
  }
  
  # separating model coefficients for different equations
  betaEq <- list()
  for( i in 1:nDep ) {
    betaEq[[ i ]] <- coef[ ( ( i - 1 ) * nReg + 1 ):( i * nReg ) ]
  }
  
  
  result <- list()
  result$beta <- coef
  result$sigma <- sigma
  result$betaEq <- betaEq
  
  return( result )
}


pmvnormWrap <- function( lower = -Inf, upper = Inf, sigma, algorithm, 
                         random.seed, nGHK = NULL, ... ) {
  
  # checking argument 'sigma'
  if( !is.matrix( sigma ) ) {
    stop( "argument 'sigma' must be a matrix" )
  } else if( nrow( sigma ) != ncol( sigma ) ) {
    stop( "argument 'sigma' must be a square matrix" )
  } else if( !all( is.numeric( sigma ) ) ) {
    stop( "argument 'sigma must be a numeric matrix" )
  }
  
  # checking argument 'lower'
  if( !all( is.numeric( lower ) ) ) {
    stop( "argument 'lower' must be numeric" )
  } else if( length( lower ) == 1 ) {
    lower <- rep( lower, nrow( sigma ) )
  } else if( length( lower ) != nrow( sigma ) ) {
    stop( "argument 'lower' must either be a single numeric value",
          " or a vector with length equal to the number of rows/columns",
          " of argument 'sigma'" )
  }
  
  # checking argument 'lower'
  if( !all( is.numeric( upper ) ) ) {
    stop( "argument 'upper' must be numeric" )
  } else if( length( upper ) == 1 ) {
    upper <- rep( upper, nrow( sigma ) )
  } else if( length( upper ) != nrow( sigma ) ) {
    stop( "argument 'upper' must either be a single numeric value",
          " or a vector with length equal to the number of rows/columns",
          " of argument 'sigma'" )
  }
  
  # check argument 'algorithm'
  algOkay <- TRUE
  ghk <- FALSE
  if( !is.list( algorithm ) && length( algorithm ) != 1 ) {
    stop( "argument 'algorithm' must be a single function",
          " or a single character string" )
  } else if( is.function( algorithm ) ) {
    algResult <- do.call( algorithm, list() )
    if( ! class( algResult )[ 1 ] %in% c( "GenzBretz", "Miwa", "TVPACK" ) ) {
      algOkay <- FALSE
    }
  } else if( is.character( algorithm ) ) {
    if( tolower( algorithm ) == "ghk" ) {
      ghk <- TRUE
    } else if( ! algorithm %in% c( "GenzBretz", "Miwa", "TVPACK" ) ) {
      algOkay <- FALSE
    }
  } else if( ! class( algorithm ) %in% c( "GenzBretz", "Miwa", "TVPACK" ) ) { 
    algOkay <- FALSE
  }
  if( !algOkay ) {
    stop( "argument 'algorithm' must be either one of the functions",
          " 'GenzBretz()', 'Miwa()', or 'TVPACK()'",
          " or one of the character strings",
          " \"GenzBretz\", \"Miwa\", or \"TVPACK\"" )
  }
  
  
  # checking argument 'random.seed'
  if( length( random.seed ) != 1 ) {
    stop( "argument 'random.seed' must be a single numerical value" )
  } else if( !is.numeric( random.seed ) ) {
    stop( "argument 'random.seed' must be numerical" )
  }
  
  # save seed of the random number generator
  if( exists( ".Random.seed" ) ) {
    savedSeed <- .Random.seed
  }
  
  # set seed for the random number generator (used by pmvnorm)
  set.seed( random.seed )
  
  # restore seed of the random number generator on exit
  # (end of function or error)
  if( exists( "savedSeed" ) ) {
    on.exit( assign( ".Random.seed", savedSeed, envir = sys.frame() ) )
  } else {
    on.exit( rm( .Random.seed, envir = sys.frame() ) )
  }
  
  if( ghk ) {
    if( is.null( nGHK ) ) {
      stop( "if the GHK algorithm is used,",
            " argument 'nGHK' must be specified" )
    } else if( length( nGHK ) != 1 ) {
      stop( "argument 'nGHK' must be a single integer value" )
    } else if( !is.numeric( nGHK ) ) {
      stop( "argument 'nGHK' must be numeric" )
    } else if( nGHK <= 0 ) {
      stop( "argument 'nGHK' must be positive" )
    }
    L <- try( t( chol( sigma ) ), silent = TRUE )
    if( class( L ) == "try-error" ) {
      warning( "the correlation matrix is not positive definite" )
      return( NA )
    }
    trunpt <- rep( NA, length( lower ) )
    above <- rep( NA, length( lower ) )
    for( i in 1:length( lower ) ) {
      if( lower[ i ] == -Inf ) {
        trunpt[ i ] <- upper[ i ]
        above[ i ] <- 1
      } else if( upper[ i ] == Inf ) {
        trunpt[ i ] <- lower[ i ]
        above[ i ] <- 0
      } else {
        stop( "if algorithm 'GHK' is used,",
              " either the lower truncation point must be '-Inf'",
              " or the upper truncation point must be 'Inf'" )
      }
    }
    sink(tempfile())
    on.exit( sink(), add = TRUE )
    result <- ghkvec( L = L, trunpt = trunpt, above = above, r = nGHK )
    result <- drop( result )
  } else {
    result <- pmvnorm( lower = lower, upper = upper, sigma = sigma,
                       algorithm = algorithm, ... )
  }
  
  return( result )
}
