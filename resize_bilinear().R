
resize_bilinear <- function(xin,yin,xout,yout,z)
	{
		library(fields)
		obj <- list(x=1:xin, y=1:yin, z = z)
		tempx <- seq(1,xin,,xout)
		tempy <- seq(1,yin,,yout)
		loc <- make.surface.grid( list( tempx,tempy))

		look <- interp.surface( obj, loc)
		return(as.surface( loc, look)$z)
	}

