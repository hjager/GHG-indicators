##--------------------------------------------------------------
## 
## ReservoirFunctions.R - Functions to load data
##
## Author: Yetta Jager
## Created 8/7/2020
## Last changed 9/19/2023
##
## Functions include
## Data loading:
##  GHG data:  1) get.Hertwich, 2) get.Barros, 3) get.EPA-Ohio
##  Reservoir data:  4) get.ReGeom, 5) get.NLAframe
##
##--------------------------------------------------------------
##--------------------------------------------------------------
##
## get.WBS reads in high-resolution NHD+ waterbodies for each
## of six TVA reservoirs
## get.area uses ReGeom parameters to estimate area from depth
## get.lookup creates a lookup table to estimate area and volume from depth based on ReGeom
## get.littoral.area estiates shallow area - untested
## sim.year simulated reservoir geometry time series based on flow in and out
##
## get.indicators calculates Draw_area and Draw_frac 
##
##--------------------------------------------------------------
get.indicators <- function(ts.in, by.year=TRUE)
{
 	 ICA = ts.in %>%
	    group_by(year) %>%
	    summarize(
		  draw_area = max(elevation, na.rm=TRUE)-min(elevation, na.rm=TRUE),
		  draw_frac = draw_area / max(elevation, na.rm=TRUE),
		  aspect_ratio = draw_area / max(geom.area, na.rm=TRUE)
		)
		
	return(ICA)
}

##-----------------------------------------------------
##
## Calculate runs of days with drawdown (elevation below threshold)
##  get.runs calculates ndrops, draw_dur, low_day_frac, fpool_dur, nrises
##
##-----------------------------------------------------
get.runs <- function(x, Xc) 
{
#   if (sum(is.na(x))>0) browser()
	y  <- x[!is.na(x)]
	R  <- rle(y <= Xc)
	ndays = sum(R$lengths)
#   focus on drawdown period(s)
	ibelow <- which(R$values)
	iabove <- which(!R$values)
	Rbelow <- R$lengths[ibelow]
	Rabove <- R$lengths[iabove]
	
	# initialize returns
	ndrops = 0; draw_dur = 0; low_day_frac = 0; nrises = 0; fpool_dur = 0
	if (length(R) > 0) {
    	# longest drawdown run
    	ndrops  = length(Rbelow)
		if (ndrops > 0) {
			k <- which.max(Rbelow)
			draw_dur  = Rbelow[k]
			nlow.days = sum(Rbelow)
			low_day_frac = nlow.days/ndays
        }
		
		# longest fullpool run
		nrises    = length(Rabove)
		if (nrises > 0) {
			k <- which.max(Rabove)
			fpool_dur = Rabove[k]
		}
	}

   	return (list(ndrops, draw_dur, low_day_frac, fpool_dur, nrises))
}

##---------------------------------------------------------------	
##
## Estimate ramping rates, up and down, from runs of like slope
##
##---------------------------------------------------------------	
get.ramps <- function(x, direction=-1) 
{	
   #direction = -1
   #x = as.numeric(x$elevation)
   ##----------------------------------------------------------
   ## Run of differences needed to estimate seasonal ramping
   ## R = both positive and negative changes in elevation
   ##----------------------------------------------------------
    slope <- c(0, diff(x)); # tack on zero to preserve length and julian dates
	if (direction==1) {
	  right.way = slope>0
	} else {
	  right.way = slope<=0
	}
	
  	R  <- rle(right.way)
	R.dt <- as.data.table(unclass(R))
	R.dt[, run_days := lengths*values]; # identifies only run lengths in desired direction
    all_days = sum(R.dt$run_days, na.rm=TRUE)
    
  	# data frame with endpoints of each run
  	run_end   = cumsum(R$lengths)
    run_start = c(1, stats::lag(run_end)[-1] + 1)
    run_start = c(0, run_end[1:length(run_end)-1]) + 1
   	ramp_dates <- data.table(run_start, run_end)
	R.dt <- cbind(R.dt, ramp_dates)
	
  	# initialize returns
  	jday_start = -1; jday_end = -1; ramp_rate = 0
  	if (length(R)>1) {
      # identify longest directional run
	   # kk refers to position in list of both draw & fill runs, not a julian date
   	   kk = which.max(R.dt$run_days)
	   # make sure julian date is in fall for drawdown, spring for refill
	   jday_start = ramp_dates[kk]$run_start
	   if ((direction==1 & jday_start <= 180) || (direction==-1 & jday_start > 180)) {
		   jday_end   = ramp_dates[kk]$run_end
		   ramp_rate  = (x[jday_end]-x[jday_start])/(jday_end-jday_start+1)
	   } else {
	     print(paste0('Longest run not in appropriate season; jday_start=', jday_start, ' direction=', direction))
	   }
   	}
  	
  	return(list(all_days, jday_start, jday_end, ramp_rate))
}

operate <- function(lookup.dt, ts) 
{
    
	# Generate time series, picking one reservoir and year for now
	# Qin <- wbm.tva.dt %>%
		# filter(year==2010 & GRanD_ID %in% tva.GRanD_IDs[3]) %>%
		# select(RiverDischarge) %>%
		# unlist() %>% as.numeric()

	# Outflows assumed equal with 3 SE random variation
	ndays <- length(Qin)
	Qse  = sd(Qin, na.rm=TRUE)/sqrt(ndays)
	Qout <- max(Qin + rnorm(ndays)*Qse*3, 0)
}

get.area <- function(zD, coeff_a3, shape)
{
    # coeff_a3 = Ltop x Wtop
	area = coeff_a3 * (1-(zD*zD))
	if (shape=='Rectangular Prism') return(area)
	else if (shape=='Rectangular Bowl') return(area * sqrt(1-zD))
	else if (shape=='Rectangular Wedge') return(area * (1-zD))
	else if (shape=='Parabolic Wedge') return(2/3 * area * (1-zD))
	else if (shape=='Elliptical Bowl') return(0.25*pi*area*sqrt(1-zD))
	else my.print('shape not found', shape)
}

# Lookup table of surface area, volume and depth
# Used in Operation-GHG.Rmd
get.lookup <- function(parms, nlayers=30, increasing=TRUE)
{
  # pass parameters
  # z = depth from surface
  shape <- as.character(parms$shape)
  coeff_a3 <- with(parms, coeff_a * Mean_L_km * Mean_W_km)
  coeff_v  <- parms$coeff_v
  Dmax     <- parms$DAM_H_m
  #Vmax    <- parms$V_GRanD_mcm

  layer.depth = Dmax/nlayers
  depth  <- vector()
  area   <- vector()
  volume <- vector()
  #for each layer 
  for (iz in seq(1,nlayers,1))
  {
      # layer depth from surface (m)
      depth[[iz]] = (iz/nlayers)*Dmax
  	  area[[iz]] <- get.area(depth[[iz]]/Dmax, coeff_a3, shape)
  	  if (iz>1) volume[[iz]] <- layer.depth * (area[iz]+area[iz-1])/2 
	  else volume[[iz]] <- layer.depth * area[iz]
  }
  # layer area and volume decrease with depth
  volume    <- coeff_v * rev(cumsum(rev(volume)))
  layer.id  <- seq(1, nlayers, 1)
  lookup.dt <- as.data.frame(cbind(layer.id, depth, area, volume)) 
  # add extra row - bottom of reservoir to aid interpolation
  # Alternatively, could set volume to zero if area is zero
  bottom.row   <- data.frame(layer.id=nlayers+1, depth=Dmax, area=0, volume=0)
  lookup.dt <- rbind(lookup.dt, bottom.row)
  if (increasing==FALSE) lookup.dt <- lookup.dt[rev(seq_len(nrow(lookup.dt))), , drop=FALSE] 
  return(lookup.dt)
}

##--------------------------------------------------------
## Simulate time series of reservoir inflow based on 
##   historical outflow, volume, and evaporation data
##--------------------------------------------------------
sim.inflow.ts <- function(ts.in)
{
  # Create time series 
   start = ts.in[1, 'date']; end = ts.in[nrow(ts.in), 'date']
   Volume      = ts(ts.in[, 'storage'], start=start, end=end) 
   Outflow     = ts(ts.in[, 'outflow'], start=start, end=end) 
   Evaporation = ts(ts.in[, 'evaporation'], start=start, end=end)
   
   # Impute missing evaporation, zero if all missing
   mean_evap   = mean(Evaporation, na.rm=TRUE)
   if (is.na(mean_evap)) mean_evap = 0
   Evaporation[is.na(Evaporation)] = mean_evap   
   
   # Back-out inflow (change in volume = outflow+evaporation-inflow)
   # NOTE: lag function conflict with tidyverse, so specify
   # This does not depend on ReGeom
   # Produces a small number of negative inflow values
   Inflow <- Outflow + Evaporation + (Volume - stats::lag(Volume))
   # Assume it wraps around...
   first.inflow <- Outflow[1] + Evaporation[1] + (Volume[1]-Volume[nrow(ts.in)])
   Inflow <- c(first.inflow, Inflow)
  
   # Augment ts data table 
   ts.in[, evaporation := ifelse(is.na(evaporation), Evaporation, evaporation)]
   ts.in[, inflow := ifelse(is.na(inflow), Inflow, inflow)]
    
   return(ts.in)
}   

##--------------------------------------------------------
# Littoral area 0-15 feet or 4.572 m (max depth of light penetration)
# Used in Operation-GHG.Rmd, called from sim.year and sim.ts
##--------------------------------------------------------
get.littoral.area <- function(depth.ts, area.ts, lookup.dt)
{
    Ldepth = 4.572
    # sort depth ascending 
    #Dlookup.dt   <- data.frame(lookup.dt %>% map_df(rev))
    # add littoral depth to today's depth
    row.id       <- findInterval(depth.ts+Ldepth, lookup.dt$depth)
    shallow.area <- area.ts - lookup.dt[row.id, 'area']
    return(unlist(shallow.area))
}

##--------------------------------------------------------
## Simulate time series of reservoir littoral area based on 
##   historical outflow, volume, and evaporation data
##--------------------------------------------------------
sim.littoral.ts <- function(ts.in, lookup.dt)
{
  # Create time series (using end date creates new records)
   start = ts.in[1, 'date']; # end = ts.in[nrow(ts.in), 'date']
   # Assign  
   Volume    = as.ts(ts.in[, 'storage'],   start=start, frequency=1); #end=end)
   Elevation = as.ts(ts.in[, 'elevation'], start=start, frequency=1); #end=end)
   
   # Impute missing Elevation
   lagElev = stats::lag(Elevation)
   Elevation = ifelse(is.na(Elevation), lagElev, Elevation) 
   
   # Derive time series of daily area and depth from ReGeom lookup table
   row.ids <- findInterval(Volume, vec=rev(lookup.dt$volume), left.open=TRUE)
   row.ids <- nrow(lookup.dt) - row.ids
   Depth   <- unlist(lookup.dt[row.ids, 'depth'])
   Area    <- unlist(lookup.dt[row.ids, 'area'])
   
   # Determine area at depth of littoral zone
   # Estimate shallow area as difference between current surface area and surface area at depth of 5' below current depth
   shallow.area <- get.littoral.area(Depth, Area, lookup.dt)
   
   # Augment ts data table from simulated geometry and outflow data above
   ts.in[, littoral.area := shallow.area][, elevation := Elevation]
   ts.in[, geom.depth := Depth][, geom.area := Area]
   
   return(ts.in)
} 
 
##--------------------------------------------------------
## Simulate one year of reservoir geometry based on 
##    daily inflows and outflows
##--------------------------------------------------------
sim.year <- function(Vmax, lookup.dt, Qin, Qout)
{
   ndays <- length(Qin)
   
   # Initialize vectors
   Area    <- vector(length=ndays)
   dArea   <- vector(length=ndays)
   Depth   <- vector(length=ndays)
   dDepth  <- vector(length=ndays)
   Volume  <- vector(length=ndays)
   dVolume <- vector(length=ndays)
   ALitt   <- vector(length=ndays)
   dALitt  <- vector(length=ndays)

   # Initialize
   # Its winter, so start with half-full reservoir
   Volume[1] = Vmax/2; 
   row.id    <- findInterval(Volume[1], lookup.dt$volume)
   Area[1]   <- unlist(lookup.dt[row.id, 'area'])
   Depth[1]  <- unlist(lookup.dt[row.id, 'depth'])
   # Determine area at depth of littoral zone
   ALitt[1]  <- get.littoral.area(Depth[1], Area[1], lookup.dt)
   for (t in 2:ndays)
   {
       # Update volume (neglecting evaporation)
   	   Volume[t] = Volume[(t-1)] + Qin[t] - Qout[t]
   	   # Not sure this is necessary, but...
   	   if (Volume[t] > Vmax) {
   	      Qout[t] = Qout[t] + (Volume[t]-Vmax)
   	      Volume[t] = Vmax
   	   }
   	   dVolume[t] <- Volume[t] - Volume[(t-1)]
   	   
   	   # Lookup depth, area, volume based on new volume
   	   row.id    <- findInterval(Volume[t], lookup.dt$volume)
   	   # Estimate shallow area as difference between 
	   #  current surface area and surface area at depth of 5' below current depth
   	   Area[t]   <- unlist(lookup.dt[row.id, 'area'])
   	   dArea[t]  <- Area[t] - Area[(t-1)]
   	   Depth[t]  <- unlist(lookup.dt[row.id, 'depth'])
   	   dDepth[t] <- Depth[t] - Depth[(t-1)]
   	   ALitt[t]  <- get.littoral.area(Depth[t], Area[t], lookup.dt)
   	   dALitt[t] <- ALitt[t] - ALitt[(t-1)]
   }

  # Create dataframe from simulated geometry and flow data above
  Res.dt <- data.frame(DOY=1:ndays, Qin=Qin, Qout=Qout, depth=Depth, dDepth=dDepth, area=Area, dArea=dArea, volume=Volume, dVolume=dVolume, ALitt=ALitt, dALitt=dALitt)
  return(Res.dt)
}

get.WBS <- function(res.name, tva6.huc.sf)
{

#res.name <- 'Douglas Lake'

	hrfile <- paste0(ghg.data.path, 'TVA sample reservoirs/',res.name,'/hiresWBS.Rdata')
	if (file.exists(hrfile)) {
		load(file=hrfile, verbose=TRUE)
		return(hrWBS.sfc)
	} else {

	this.sf  <- subset(tva6.huc.sf, RES_NAME==res.name)
	this.res <- droplevels(data.table(st_drop_geometry(this.sf)))
	
	# How many HUC4s associated with this reservoir?
	levels <- unique(na.omit(levels(this.res$HUC4)))
	print(paste0('HUC4s associated with ', res.name, ': ', levels))
	nhuc4 <- length(levels)
	if (nhuc4==0)	{
	   print(paste0('No associated HUC4s, returning NULL'))
	   return(NULL)
	} 
	if (nhuc4==1) {
		# Only one HUC4 associated with reservoir
	    # versus combination of HUC4s
		hr.file <- paste0(ghg.data.path,'NHDplus/NHDPLUS_H_', levels[[1]], '_HU4_GDB.gdb')
		  if (file.exists(hr.file)) {
		  print(paste0('Reading ', hr.file))
			hrWBS.sf <- get_hr_data(gdb=hr.file, layer=c('NHDWaterbody'), proj=latlon.proj)
		  }
	 } else {
		# Read high-res polygon for this reservoir
		# There can be multiple polygons and huc4s per polygon
		# but we don't care about multiple polygons here
		resWBS <- list()
		for (kk in 1:nhuc4) 
		{
			hr.file <- paste0(ghg.data.path,'NHDplus/NHDPLUS_H_', levels[[kk]], '_HU4_GDB.gdb')
			if (file.exists(hr.file)) {
			   print(paste0('Reading ', hr.file))
			   resWBS[[kk]] <- get_hr_data(gdb=hr.file, layer=c('NHDWaterbody'), proj=latlon.proj)
			} else {
			  print(paste0('Missing file ', hr.file))
			  return(0)
			}
		}
		hrWBS.sf <- bind_rows(resWBS)
	 }
	 
	 # End of reading boundary data

	# Combine if reservoir spans more than one HUC4
	# Focus on higher order reaches
	setnames(hrWBS.sf, c('VPUID','REACHCODE'), c('HUC04','ReachCode'), skip_absent=TRUE)
	cols <- c('HUC04','ReachCode','COMID','GNIS_Name')	
	hrWBS.sf <- hrWBS.sf  %>%
  	select(c(cols, 'Shape')) %>%
	mutate(across(cols, as.factor))
	#mutate(across(where(is.character), as.factor))

	# May need to project before using st_join?
	# For some reservoirs (Douglas), need to cast to avoid 
	# Error in CPL_geos_make_valid(x) : Evaluation error: ParseException: Unknown WKB type 12.
	hrWBS.sfc <- hrWBS.sf %>%
	     st_cast('MULTIPOLYGON') %>%
		 st_make_valid() %>%
		 st_transform(tn.proj) 
		 
    this.sfc <- this.sf %>% 
		st_transform(tn.proj) %>% 
		st_make_valid()
	
    hrWBS.sfc <- st_join(this.sfc, hrWBS.sfc, join=st_overlaps, left=TRUE) %>% st_make_valid() 
	  	
	# save data
	print(paste0('Saving data for ', res.name))
	hrfile <- paste0(ghg.data.path, 'TVA sample reservoirs/',res.name,'/hiresWBS.Rdata')
	save(file=hrfile, hrWBS.sfc)
  }

  return(hrWBS.sfc)
  
}
##--------------------------------------------------------------
##
## get.Lines reads in high-resolution NHD+ flowlines for each
## of six TVA reservoirs
##
##--------------------------------------------------------------
get.Lines <- function(res.name, tva6.huc.sf)
{

	hrfile <- paste0(ghg.data.path, 'TVA sample reservoirs/',res.name,'/hiresLines.Rdata')
	if (file.exists(hrfile)) {
		load(file=hrfile, verbose=TRUE)
		return(hrLines.sfc)
	} else {

	this.res <- droplevels(data.table(st_drop_geometry(subset(tva6.huc.sf, RES_NAME==res.name))))
	
	# How many HUC4s associated with this reservoir?
	levels <- unique(na.omit(levels(this.res$HUC4)))
	print(paste0('HUC4s associated with ', res.name, ': ', levels))
	if (length(levels)==0)	{
		  print(paste0('No associated HUC4s, returning NULL'))
		  return(NULL)
	}

	# Only one HUC4 associated with reservoir
	# versus combination of HUC4s
	nhuc4 <- length(levels)
 if (nhuc4==1) {
   	hr.file <- paste0(ghg.data.path,'NHDplus/NHDPLUS_H_', levels[[1]], '_HU4_GDB.gdb')
 	  if (file.exists(hr.file)) {
   	  print(paste0('Reading ', hr.file))
 		   hrLines.sf <- get_hr_data(gdb=hr.file, layer=c('NHDFlowline'), proj=latlon.proj)
 	  }
 } else {
	# Read high-res flowlines for this reservoir
	# There can be multiple polygons and huc4s per polygon
 	# but we don't care about multiple polygons here
	resLines <- list()
	for (kk in 1:nhuc4) 
	{
		hr.file <- paste0(ghg.data.path,'NHDplus/NHDPLUS_H_', levels[[kk]], '_HU4_GDB.gdb')
		if (file.exists(hr.file)) {
	 	   print(paste0('Reading ', hr.file))
		   resLines[[kk]] <- get_hr_data(gdb=hr.file, layer=c('NHDFlowline'), proj=latlon.proj)
 	    } else {
		  print(paste0('Missing file ', hr.file))
		  return(0)
		}
	}
	hrLines.sf <- bind_rows(resLines)
 }

	# Combine if reservoir spans more than one HUC4
	# Focus on higher order reaches
	min.order <- 3
	hrLines.sf  %>%
		filter(StreamOrde > min.order) %>%
		make_standalone() %>%
		mutate(across(where(is.character), as.factor))

	# Plot downloaded high-res data
	hrLines.sfc <- hrLines.sf %>%
		  st_transform(tn.proj)

	# save data
	hrfile <- paste0(ghg.data.path, 'TVA sample reservoirs/',res.name,'/hiresLines.Rdata')
	save(file=hrfile, hrLines.sfc, hrLines.sf)
  }

  return(hrLines.sfc)
  
}

##--------------------------------------------------------------
##
## Load Hertwich global methane emissions data for reservoirs
## 
##--------------------------------------------------------------
get.Hertwich <- function(in.path)
{
	# if needed, retrieve data
	if (!exists('Hertwich.dt'))
	{
	   rfilename <- paste0(in.path,'Barros-Hertwich-Scherer/','Hertwich.Rdata')
	   if (file.exists(rfilename)) {
		 load(file=rfilename, verbose=TRUE); 
		 print(paste0('Loaded ', rfilename))
	   } else {
		 filename <- paste0(in.path,'Barros-Hertwich-Scherer/', 'Hertwich-es401820p_si_002.xlsx')
		 Hertwich.dt <- as.data.table(xlsx::read.xlsx(filename, sheetName='MyDams', startRow=3, stringsAsFactors=TRUE))
		 Hertwich.dt$source <- as.factor(Hertwich.dt$source)
		 old.names <- names(Hertwich.dt)
		 new.names <- c('source','reservoir','lat.dd','long.dd','id','age.y','CO2.mgC.m2.d','CH4.mgC.m2.d','depth.m','area.km2','vol.km3','capacity.MW','elec.GWh.y','NPP.gC.m2.y','residence.time.d','TP.mg.m2.y','DOC.mg.m2.y','catchment.area.km2','emission.code')
		 setnames(Hertwich.dt, old.names, new.names)
		 Hertwich.dt <- Hertwich.dt[!is.na(long.dd),]
		 Hertwich.dt$source <- as.factor(Hertwich.dt$source)

#		 Hertwich.dt$source <- rep('Hertwich (2013)', nrow(Hertwich.dt))
		 save(file=rfilename, Hertwich.dt)
		 print(paste0('Saved Hertwich.dt in ', rfilename))
	   }
	}
	return(Hertwich.dt)
}
##--------------------------------------------------------------
##
## Load National Inventory of Dams dataset
## if use.QA=TRUE, use Carly Hansen's data
## 
##--------------------------------------------------------------
get.NID <- function(in.path, use.QA=TRUE)
{
	if (!exists('NID.dt'))
	{
	   nfilename <- paste0(in.path, 'NID/NID.Rdata')
	   if (file.exists(nfilename)) {
		 load(file=nfilename, verbose=TRUE); 
		 print(paste0('Loaded ', rfilename))
	   } else {
	     # New QA'd dataset from Carly 2/22/2021 does not include
		 # EHA data
		 if (use.QA) filename <- paste0(in.path, 'NID/NID_2019_NPD_QA.csv')
		 else filename <- paste0(in.path, 'NID/NID2019_U.csv')
		 # This takes a fairly long time
		 NID.dt <- data.table::data.table(read.csv(filename, header=TRUE, na.strings=' ', stringsAsFactors=TRUE))
		 #num.list <- c('LATITUDE','LONGITUDE','SURFACE_AREA','NID_STORAGE','VOLUME','NID_HEIGHT','DAM_HEIGHT','YEAR_COMPLETED','WIDTH_OF_LOCKS','SPILLWAY_WIDTH','NUMBER_OF_LOCKS','LENGTH_OF_LOCKS','DAM_LENGTH','STRUCTURAL_HEIGHT','HYDRAULIC_HEIGHT','MAX_DISCHARGE','MAX_STORAGE','NORMAL_STORAGE','DRAINAGE_AREA')
		 NID.dt[, LONGITUDE := -abs(LONGITUDE)]; # assume error in positive values
		 NID.dt <- NID.dt[LATITUDE>0 & LONGITUDE<0,]; # just remove dams without location data
		 save(file=nfilename, NID.dt)
		 print(paste0('Saved NID.dt in ', nfilename))
	   }
	}
	return(NID.dt)
}

##--------------------------------------------------------------
## Same as Barros data?
## Scherer-Pfister GHG data, global, from Illsa Ocko
## Paper: Hydropower’s biogenic carbon footprint 
## Laura Scherer Stephan Pfister Institute of Environmental Engineering, ETH Zurich, 8093 Zurich, Switzerland
## Laura Scherer, e-mail: scherer@ifu.baug.ethz.ch
## 
##--------------------------------------------------------------
get.SP <- function(in.path, sheet)
{
	# if needed, retrieve data
	if (!exists('SchererPfister.dt'))
	{
	   bfilename <- paste0(out.path,'Barros-Hertwich-Scherer/', 'SchererPfister.Rdata')
	   if (file.exists(bfilename)) {
		 load(file=bfilename, verbose=TRUE); 
		 print(paste0('Loaded ', rfilename))
	   } else {
		 xfile <- paste0(in.path,'Barros-Hertwich-Scherer/','Scherer-S1 Table.xlsx')
		 SchererPfister.dt <- data.table(openxlsx::read.xlsx(xfile, sheet=sheet, startRow=2))
		 save(file=bfilename, SchererPfister.dt)
		 print(paste0('Saved SchererPfister.dt in ', bfilename))
	   }
	} else print('SchererPfister.dt already exists')
	return(SchererPfister.dt)
}

##--------------------------------------------------------------
## 
## Deemer et al. data, more data
## 
##--------------------------------------------------------------
get.deemer <- function(in.path)
{
	   dfilename <- paste0(in.path, 'Deemer.Rdata')
	   if (file.exists(dfilename)) {
		 load(file=dfilename, verbose=TRUE); 
		 print(paste0('Loaded ', dfilename))
	   } else {
		 xfile <- paste0(in.path, 'Supplementary_Reservoir_GHG_Data_20200105.xlsx')
		 print(xfile)
		 deemer.dt <- data.table(openxlsx::read.xlsx(xfile, sheet='Supplementary_Reservoir_GHG_Dat', startRow=2))

		 old.names <- c("Latitude","Longitude","System","Alpine.(cutoff.1000m.elevation)","Installed.Hydropower.Capacity?",
		 "Primary.Use.of.Reservoir",
		 "mg.CH4-C.m-2.d-1.Diffusive.+.Ebullitive","mg.CH4-C.m-2.d-1.Diffusive.Only","mg.CH4-C.m-2.d-1.Ebullitive.Only",
		 "Reported.Mean.Annual.Precipitation.(mm)",
		 "Abs.Latitude","Mean.Water.Temperature.(degrees.C)",
		 "Modeled.Mean.Annual.Air.Temperature.(degrees.C)",
		 "DOC.(mg/L)",'Age.(yrs)','Chlorophyll.a.(ug/L)','Modeled.Mean.Annual.Precipitation.(mm)','Modeled.Runoff.(mm.y-1)','Mean.Depth.(m)','Surface.Area.(km2)','Catchment.Area.(km2)','Residence.Time.(days)','Volume.(km3)','Citations.for.GHG.Data')
		 
		 new.names <- c("lat.dd","long.dd","reservoir",'Alpine','Capacity','Purpose','Both.CH4.mgC.m2.d','Diff.CH4.mgC.m2.d','Ebu.CH4.mgC.m2.d','MAP.mm','abs.lat.dd','Water.temp.degC','Air.Temp.degC','DOC.mgL','age.y','Chla.ugL','Mod.MAP.mm','Runoff.mm.y','depth.m','reservoir.area.km2','catchment.area.km2','residence.time.d','Volume.km3','source')
		 setnames(deemer.dt, old.names, new.names)
  
		 droplist <- c('Amazon','X38','Supplementary.Citations','Alpine', 'abs.lat.dd','Catchment.Area.:.Surface.Area')
		 deemer.dt <- subset(deemer.dt, select=!(names(deemer.dt) %in% droplist))
		 
		 # Add variables consistent with EPA-Ohio
 		 deemer.dt$emission.code <- with(deemer.dt, ifelse(is.na(Ebu.CH4.mgC.m2.d), 'D', 'D&B'))

		 # from https://www.usbr.gov/projects/index.php?id=294
		 # P-14116
		 deemer.dt[reservoir=='Keechelus', Capacity := 'yes']
		 deemer.dt[reservoir=='Keechelus', reservoir.area.km2 := 10.36]
		 deemer.dt[reservoir=='Keechelus', Volume.km3 := 0.1946]
		 deemer.dt[reservoir=='Keechelus', catchment.area.km2 := 141]
		 # www.env.nm.gov/wp-content/uploads/sites/25/2019/10/Lakes-2006.pdf
		 deemer.dt[reservoir=='Lower Charrette Lake', reservoir:='Lower Charette Lake']
		 deemer.dt[reservoir=='Lower Charette Lake', Capacity := 'no']
		 deemer.dt[reservoir=='Lower Charette Lake', reservoir.area.km2 := 1.2]
		 deemer.dt[reservoir=='Lower Charette Lake', catchment.area.km2 := 560.004]
		 deemer.dt[reservoir=='Lower Charette Lake', Volume.km3 := 0.000370]
		 deemer.dt[reservoir=='Lower Charette Lake', depth.m := 14.9]
		 deemer.dt[reservoir=='Lower Charette Lake', Chla.ugL := 2.74]
		 deemer.dt[reservoir=='Lower Charette Lake', Trophic.Status := 'Mesotrophic']
		 deemer.dt[reservoir=='Upper Charette Lake', Trophic.Status := 'Eutrophic']

		 # factorize
		 cols <- c('Capacity','Country','Region','Trophic.Status','Purpose','reservoir','emission.code','source')
		 deemer.dt[, (cols) := lapply(.SD, 'as.factor'), .SDcols=cols]
		 
		 save(file=dfilename, deemer.dt)
		 print(paste0('Saved deemer.dt in ', dfilename))
	   }
	return(deemer.dt)
}

##--------------------------------------------------------------
##
## Load Barros global methane emissions data for reservoirs
## 
##--------------------------------------------------------------
get.Barros <- function(in.path)
{
	# if needed, retrieve data
	if (!exists('Barros.dt'))
	{
	   bfilename <- paste0(in.path, 'Barros.training.Rdata')
	   if (file.exists(bfilename)) {
		 load(file=bfilename, verbose=TRUE); 
		 print(paste0('Loaded ', rfilename))
	   } else {
		 # sheet training
		 filename  <- paste0(in.path, 'Barros.xlsx')
		 Btrain.dt <- data.table(xlsx::read.xlsx(filename, sheetName='training', startRow=2, stringsAsFactors=TRUE))
		 setnames(Btrain.dt, c('ID','Name'), c('id','reservoir'))
		 Btrain.dt[, 3:17] <- Btrain.dt[, lapply(.SD, as.numeric), .SDcols=3:17]

		 # Plant data were estimated to get a global estimate
	     # sheet plant; data were predicted, not measured
		 # filename  <- paste0(in.path, 'Barros.xlsx')
		 # Bplant.dt <- data.table(xlsx::read.xlsx(filename, sheetName='plant', startRow=1))
		 # Bplant.dt[, name := as.factor(name)]
		 # Bplant.dt[, purpose := as.factor(purpose)]
#		 Bplant.dt[, lapply(.SD, 'as.factor'), .SDcols=sapply(Bplant.dt, is.character)]
	     #Barros.dt <- merge(Btrain.dt, Bplant.dt, by=c('name'))

		 Barros.dt <- Btrain.dt[, volume.area := NULL]
		 Barros.dt$mean.energy.GWh.y <- Barros.dt$area.km2/(Barros.dt$Area.Electricity)

		 Barros.dt$source  <- rep('Barros et al. 2011', nrow(Barros.dt))
#
	     save(file=bfilename, Barros.dt)
		 print(paste0('Saved Barros.dt in ', bfilename))
	   }
	} else print('Barros.dt already exists')
	return(Barros.dt)
}

##--------------------------------------------------------------
##
## Load GRanD global reservoir dataset
## 
##--------------------------------------------------------------
# GRaND database
get.ReGeom <- function(in.path)
{
	if (!exists('ReGeom.dt'))
	{
	   rfilename <- paste0(in.path, 'USGeom.Rdata')
	   if (file.exists(rfilename)) {
		 load(file=rfilename, verbose=TRUE); 
		 print(paste0('Loaded ', rfilename))
	   } else {
		 filename <- paste0(in.path, 'ReGeomData_WOW_V1.xlsx')
		 ReGeom.dt <- as.data.table(read.xlsx(filename, sheet='ReGeomData', header=TRUE))
		 ReGeom.dt$RES_NAME <- as.factor(ReGeom.dt$RES_NAME)
		 ReGeom.dt$DAM_NAME <- as.factor(ReGeom.dt$DAM_NAME)
		 ReGeom.dt$COUNTRY  <- as.factor(ReGeom.dt$COUNTRY)
		 ReGeom.dt$shape    <- as.factor(ReGeom.dt$shape)
		 USReGeom.dt <- subset(ReGeom.dt, COUNTRY=='United States')
		 save(file=rfilename, USReGeom.dt)
		 print(paste0('Saved USReGeom.dt in ', rfilename))
	   }
	}
	return(USReGeom.dt)
}
##--------------------------------------------------------------
##
## Get projected and clipped state boundaries
## 
##--------------------------------------------------------------
get.states <- function(my.proj, bbox.latlon)
{
	data(us_states); 
	states.sf  <- st_set_crs(us_states, st_crs(bbox.latlon))
	states.sf  <- st_crop(states.sf, bbox.latlon)
	states.sfc <- st_transform(us_states, my.proj)
	return(states.sfc)
}
##--------------------------------------------------------------
##
## Load EPA GHG dataset
## 
##--------------------------------------------------------------
get.EPA.Ohio <- function(ghg.data.path)
{
   if (!exists('EPA-Ohio.dt'))
   {
	   sfilename <- paste0(ghg.data.path,'EPA-GHG-Data/EPA-Ohio.Rdata')
	   if (file.exists(sfilename)) {
		 load(file=sfilename, verbose=TRUE)
		 print(paste0('Loaded ', sfilename))
	   } else {
		 filename <- paste0(ghg.data.path, 'EPA-GHG-Data/EPA2020-siTable.txt')
		 EPA-Ohio.dt <- data.table(read.table(filename, skip=0, sep=' ', header=TRUE, stringsAsFactors=TRUE))
		 old.names <- c('Lake_Name') 
		 new.names <- c('reservoir')
		 setnames(EPA-Ohio.dt, old.names, new.names, skip_absent=TRUE)
		 EPA-Ohio.dt[, reservoir.area.km2 := reservoir.area.m2*1.0e-6][, reservoir.area.m2 := NULL]
		 save(file=sfilename, EPA.Ohio.dt)
		 print(paste0('Saved EPA.Ohio.dt in ', sfilename))
	   }
   } else print('EPA-Ohio.dt already exists')
   return(EPA.Ohio.dt)          
}
##---------------------------------------------------------------
## 
## Project and clean polygon layers
## precision is <1 (lower values are coarser, inverse of tolerance? go figure)
## grid size should be smaller than tolerance
## CHECK RESULTS OF THIS TO GET BEST TOLERANCE
##
##---------------------------------------------------------------
clean <- function(sf.layer, my.tol)
{
	if (my.tol > 0)
	{
	   sf.layer <- sf.layer %>% 
		  st_set_precision(1/my.tol) %>% 
		  lwgeom::st_snap_to_grid(my.tol/2) %>%
		  st_make_valid()
	}
	return(sf.layer)
}
##--------------------------------------------------------------
##
## Read NHDplus-high resolution waterbody boundary data 
## Subset to Region 6 first to save processing time
## NO LONGER USING THIS - webscrape data
## 
##--------------------------------------------------------------
get.nhdwbhr <- function()
{
	name.list <- c('douglas','wattsBar','allatoona','fontana','guntersville','hartwell')
	SF <- list()
	ii <- 0
	for (res in name.list)
	{
		ii <- ii+1
		#res = 'fontana'
		hfilename <- paste0(ghg.data.path, 'TVA sample reservoirs/reservoirPolys_nhdplushr/', res,'_nhdplushr.shp')
		wbhr.sf <- sf::st_read(hfilename, stringsAsFactors=TRUE) %>% mutate(res.name = res) 
		NHDPlusID <- paste0(unique(as.factor(wbhr.sf$NHDPlusID)))
		ReachCode <- paste0(unique(as.factor(wbhr.sf$ReachCode)))

		# Get rid of z-dimension and project 
		wbhr.sf <- wbhr.sf %>%
		   st_zm(drop=TRUE, what='ZM')  %>%
		   mutate(geom_latlon = geometry) %>%
		   st_transform(crs=tn.proj)

		# combine multiple polygons
		if (nrow(wbhr.sf) > 1)
		{
			wbhr.sf <- wbhr.sf %>%
			    st_union() %>%
				clean(my.tol=10) %>%
				st_as_sf() 
		}
		
		#wbhr.sf$NHDPlusID <- NHDPlusID
		SF[[ii]] <- wbhr.sf 
	}
	return(SF)
}

##--------------------------------------------------------------
##
## Load LakeCat data 
## www.epa.gov/national-aquatic-resource-surveys/lakecat-dataset-readme#nhd
## 
##--------------------------------------------------------------
get.LakeCat <- function(in.path, get.atts=FALSE)
{
    require('reticulate')
	if (!exists('NLAframe.sf'))
	{
	   out.file <- paste0(in.path, 'LakeCat-NHD+/LakeCat.Rdata')
	   if (file.exists(out.file)) {
		 load(file=out.file, verbose=TRUE); 
		 print(paste0('Loaded ', out.file))
	   } else {
	   # LakeCat - allBasins.shp
 	     rfilename <- paste0(in.path, 'LakeCat-NHD+/LkCat_Frame_min/shps/allBasins.shp')
		 LakeCat.sf <- sf::st_read(rfilename, stringsAsFactors=TRUE)
		 if (get.atts==TRUE)
		 {
			 np <- import("numpy")
			 rfilename <- paste0(in.path, 'LakeCat-NHD+/LkCat_Frame_min/LakeCat_npy/onNet_LakeCat.npz')
			 npz1 <- np$load(rfilename)
			 npz1$f[["vpus"]]
			 LakeCat.sf$InStreamCat <- as.factor(LakeCat.sf$InStreamCat)
			 LakeCat.sf <- subset(LakeCat.sf, !is.na(COMID))
		 }
		 
		 ## Save results
		 save(file=out.file, LakeCat.sf)
		 print(paste0('Saved LakeCat.sf in ', out.file))
	   }
	} else print('LakeCat.dt already exists')	
	return(LakeCat.sf)
}
##--------------------------------------------------------------
##
## Load HydroLakes data
## www.hydrosheds.org/page/hydrolakes data
## global lakes with surface area >= 10 ha 
## lake pour points
## Messager, M.L., Lehner, B., Grill, G., Nedeva, I., Schmitt, O. (2016): Estimating the volume and age of water stored in global lakes using a geo-statistical approach. Nature Communications: 13603. doi: 10.1038/ncomms13603.
## 
##--------------------------------------------------------------
get.HydroLAKES <- function(data.path)
{
#data.path <- ghg.data.path
	if (!exists('HydroLAKES.sf'))
	{
	   #data.path <- ghg.data.path
	   out.file.h <- paste0(data.path, 'HydroLAKES/HydroLAKES.Rdata')
	   if (file.exists(out.file.h)) {
		 load(file=out.file.h, verbose=TRUE); 
		 print(paste0('Loaded ', out.file.h))
	   } else {
	   # LakeCat - allBasins.shp
 	     hfilename <- paste0(data.path, 'HydroLAKES/HydroLAKES_polys_v10.shp')
		 HydroLAKES.sf <- sf::st_read(hfilename, stringsAsFactors=TRUE)
		 
		 ## Subset to CONUS first to save processing time
		 HydroLAKES.sf <- HydroLAKES.sf %>% filter(Country=='United States of America')
		 st_set_crs(HydroLAKES.sf, 4326)
		 
		 ## Project coordinates
		 my.proj <- 5070
		 data(us_states); #4269
		 st_set_crs(us_states, 4269)
		 conus.sf <- st_transform(us_states, my.proj)
 
         print(paste0('Projecting HydroLAKES.sf to ', my.proj))
		 HydroLAKES.sf <- st_transform(HydroLAKES.sf, crs=my.proj)
		 HydroLAKES.sf <- sf::st_buffer(HydroLAKES.sf, dist=0)
		 HydroLAKES.sf <- st_make_valid(HydroLAKES.sf)
		
		 HydroLAKES.sf <- st_crop(x=HydroLAKES.sf, y=conus.sf)

 		# Replace -9999 with NA, probably a way to do this using replace_na
		 HydroLAKES.sf$Dis_avg   <- ifelse(HydroLAKES.sf$Dis_avg<0, NA, HydroLAKES.sf$Dis_avg)
		 HydroLAKES.sf$Res_time  <- ifelse(HydroLAKES.sf$Res_time<0, NA, HydroLAKES.sf$Res_time)
		 HydroLAKES.sf$Wshd_area <- ifelse(HydroLAKES.sf$Wshd_area<0, NA, HydroLAKES.sf$Wshd_area)
		 HydroLAKES.sf$Elevation <- ifelse(HydroLAKES.sf$Elevation<0, NA, HydroLAKES.sf$Elevation)
  		 HydroLAKES.sf$Lake_type <- as.factor(HydroLAKES.sf$Lake_type)

		## Save results
		 save(file=out.file.h, HydroLAKES.sf)
		 print(paste0('Saved HydroLAKES.sf in ', out.file.h))
	   }
	} else print('HydroLAKES.sf already exists')
	
	return(HydroLAKES.sf)
}
##--------------------------------------------------------------
##
## Load NLA list frame (2012 survey)
## Note: this is a large geodatabase (.gdb) read with sf
## ftype='Reservoir' does not mean what we think though
## huc8 and huc2 have 'H' pre-appended to code
## cntyname, statecty are available
## area_ha and perim_km
## 
##--------------------------------------------------------------
get.NLAframe <- function(data.path, layer='point')
{
    #data.path <- ghg.data.path
	gdb.file <- paste0(data.path, 'EPA-GHG-Data/NLA_Sample_Frame.gdb')
	
	# Dataset has both points (centroids?) and lake polygons
	print(sf::st_layers(gdb.file))
	if (layer=='point') {
		NLAframe.sf <- sf::st_read(gdb.file, layer='NLA_Sample_Frame_Points')
	} else if (layer=='polygon') {
		NLAframe.sf <- sf::st_read(gdb.file, layer='NLA_Sample_Frame_Polys')
	}	
		 
	# Get rid of Z and M dimensions, zero information here
	NLAframe.sf <- sf::st_zm(NLAframe.sf)
		 
	# Set projection to latlong
	NLAframe.sf  <- st_set_crs(NLAframe.sf, "+proj=longlat")
		 
	# Alaska doesn't seem to be in here.
	NLAframe.sf$state <- as.factor(NLAframe.sf$state)
	NLAframe.sf$area.km2 <- NLAframe.sf$area_ha/100
	setnames(NLAframe.sf, c('lat_dd83','lon_dd83','statecty','perim_km','Shape'), c('lat.dd','long.dd','fips','perim.km','geometry'), skip_absent=TRUE)
		 #with(NLAframe.sf, table(state))
		 #with(NLAframe.sf, table(ftype))
	# Tell sf that we renamed geometry column
	st_geometry(NLAframe.sf) <- 'geometry'
	NLAframe.sf <- subset(NLAframe.sf, source %in% c('NLA2012Frame Only') & ftype %in% c('LakePond','Reservoir') & state != 'HI')
	return(NLAframe.sf)
}
##--------------------------------------------------------------
##
## Load CONUS? reservoir data from HydroSource, EHA
## 
##--------------------------------------------------------------
get.EHA <- function(in.path)
{
	# if needed, retrieve data
	if (!exists('EHA.dt'))
	{
	   rfilename <- paste0(out.path, 'EHA.Rdata')
	   if (file.exists(rfilename)) {
		 load(file=rfilename, verbose=TRUE); 
		 print(paste0('Loaded ', rfilename))
	   } else {
		 #shape.file <- paste0(in.path, 'EHA/ORNL_EHAHydroPlant_FY2020_revised/Plant_operational_ExternalGISFY2020revised.shp')
		 attr.file <- paste0(in.path,'EHA/ORNL_EHAHydroPlant_FY2020_revised/ORNL_EHAHydroPlant_FY2020revised.xlsx')
		 EHA.dt <- as.data.table(xlsx::read.xlsx(attr.file, startRow=1, sheetName='Operational', stringsAsFactors=TRUE))
		 names(EHA.dt) <- toupper(names(EHA.dt))
		 old.names <- c('PTNAME','OWTYPE','PT_OWN','HUC','NID_ID')
		 NID.names <- c('DAM_NAME','OWNER_TYPE','OWNER_NAME','HUC10','NIDID')
		 setnames(EHA.dt, old.names, NID.names)
#		 EHA.sf <- st_read(shape.file, stringsAsFactors=TRUE)
#		 print(st_layers(shape.file))
		 save(file=rfilename, EHA.dt)
		 print(paste0('Saved EHA.dt in ', rfilename))
	   }
	}
	return(EHA.dt)
}
##--------------------------------------------------------------
##
## Load linkage map to EHA, NID, HydroLAKES
## 
##--------------------------------------------------------------
get.HILARRI <- function(in.path)
{
	hfile <- paste0(in.path, 'HILARRI/HILARRI_202104_Public.Rdata')
	load(hfile, verbose=TRUE)
	hilarri.dt <- data.table(hilarri)
	return(hilarri.dt)
}

##--------------------------------------------------------------
##
## Map continuous values as bubbles (size and color)
## 
##--------------------------------------------------------------
map.conus <- function(SF, state, fill.var, fill.name, pal, zmax, upper=7)
{
  require('ggsn')
  #SF <- subset(SF, !is.na(get(fill.var)))
  L1 <- ggplot() + geom_sf(data=SF, aes(size=get(fill.var), fill=get(fill.var)), color='black', shape=21, alpha=1) + scale_size(name=fill.name, range=c(1,upper), guide=FALSE) 
  L2 <- L1 + geom_sf(data=state, color='black', fill=NA)
  if (is.na(pal)) {
      L3 <- L2 + scale_fill_viridis(name=fill.name, option="A", na.value=NA, discrete=FALSE, direction=-1, limits=c(0, zmax)) 
      #+ scale_color_viridis(name=fill.name, option="A", na.value='grey70', discrete=FALSE) 
  } else {
      L3 <- L2 + scale_fill_distiller(name=fill.name, type='seq', palette=pal, direction=1, na.value='grey70', limits=c(0, zmax))  
  }
  L4 <- L3 + ggsn::north(data=state, location="bottomleft", scale=0.15)	
  #scalebar(data=SF, location='bottomleft', box.color='black', transform=TRUE, dist=1000, dist_unit='km', st.dist=0.1, st.bottom=TRUE, st.size=2, height=0.05, model = 'WGS84'); # st.size is text 
  L5 <- L4 + theme(axis.title=element_text(size=14), axis.text=element_text(size=14), legend.title=element_text(size=12), legend.text=element_text(size=12), panel.grid.major= element_line(color="grey50", linetype="dashed", size=0.3),  panel.background=element_rect(fill="ghostwhite"))
  return(L5)
}

##--------------------------------------------------------------
##
## Map discrete values as bubbles (size and color)
## 
##--------------------------------------------------------------
map.conus2 <- function(SF, state, fill.var, fill.name, pal, size)
{
  require('ggsn')
  L1 <- ggplot() + geom_sf(data=SF, aes(fill=get(fill.var)), color='black', shape=21, size=size, alpha=1) 
  L2 <- L1 + geom_sf(data=state, color='black', fill=NA)
  if (is.na(pal)) {
     L3 <- L2 + scale_fill_viridis_d(name=fill.name, option="A", na.value=NA, direction=-1) 
  } else {
     L3 <- L2 + scale_fill_discrete(name=fill.name, palette=pal, na.value='grey70')  
  }
  L4 <- L3 + ggsn::north(data=state, location="bottomleft", scale=0.15)	
  #scalebar(data=SF, location='bottomleft', box.color='black', transform=TRUE, dist=1000, dist_unit='km', st.dist=0.1, st.bottom=TRUE, st.size=2, height=0.05, model = 'WGS84'); # st.size is text 
  L5 <- L4 + theme(axis.title=element_text(size=14), axis.text=element_text(size=14), legend.title=element_text(size=12), legend.text=element_text(size=12), panel.grid.major= element_line(color="grey50", linetype="dashed", size=0.3),  panel.background=element_rect(fill="ghostwhite"))
  return(L5)
}

##--------------------------------------------------------------
##
## Map discrete values as pie charts
## 
##--------------------------------------------------------------
map.scatterpie <- function(SF, state, fill.var, fill.name, pal, radius=2)
{
  require('ggsn')
  require('scatterpie')
  L1 <- ggplot() + geom_sf(data=SF, aes(group=get(fill.var), r=radius, color='black', shape=21, alpha=1)) 
  L2 <- L1 + geom_sf(data=state, color='black', fill=NA)
  L3 <- L2 + geom_scatterpie() 
  ## CONUS legend bottom right
  L4 <- L3 + geom_scatterpie_legend(d$radius, x=-75, y=30) 
  if (is.na(pal)) {
      L5 <- L4 + scale_fill_viridis(name=fill.name, option="A", na.value=NA, discrete=FALSE, direction=-1, limits=c(0, zmax)) 
  } else {
      L5 <- L4 + scale_fill_distiller(name=fill.name, type='seq', palette=pal, direction=1, na.value='grey70', limits=c(0, zmax))  
  }
  L6 <- L5 + ggsn::north(data=state, location="bottomleft", scale=0.15)	
  L7 <- L6 + theme(axis.title=element_text(size=14), axis.text=element_text(size=14), legend.title=element_text(size=12), legend.text=element_text(size=12), panel.grid.major= element_line(color="grey50", linetype="dashed", size=0.3),  panel.background=element_rect(fill="ghostwhite"))
  return(L7)
}
