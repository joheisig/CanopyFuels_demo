fast_countover <- function(x, t) {
  .Call(`_lidR_fast_countover`, x, t)
}
fast_table <- function(x, size = 5L) {
  .Call(`_lidR_fast_table`, x, size)
}

lidar_height_metrics = function(x, y, z, dz = 1, th = 2){
  z = z[z < 50 & z >= 0] # exclude outliers
  z = na.omit(z)
  n = length(z)
  zmax = max(z)
  zmean = mean(z)
  zsd = stats::sd(z)
  zcv = zsd / zmean * 100
  ziqr = IQR(z)
  probs = c(0.01, seq(0.05, 0.95, 0.05), 0.99)    # modified for 1 and 99%
  zq = as.list(stats::quantile(z, probs))
  names(zq) = paste0("zq", probs * 100)
  pzabovex = lapply(th, function(x) sum(z > x)/n *100)
  names(pzabovex) = paste0("pzabove", th)
  pzabovemean = sum(z > zmean)/n * 100
  if (zmax <= 0) {
    d = rep(0, 9)
  }
  else {
    breaks = seq(0, zmax, zmax/10)
    d = findInterval(z, breaks)
    #d = fast_table(d, 10)
    d = table(d)
    d = d/sum(d) * 100
    d = cumsum(d)[1:9]
    d = as.list(d)
  }
  names(d) = paste0("zpcum", 1:9)
  
  s = as.numeric(zq)[c(-1,-21)]
  cbh = round(s[which.max(diff(s))+1], 2)  # after Chamberlain et al 2021
  
  zmean_grass = mean(z[z < 0.4])
  zmean_shrub = mean(z[z > 0.4 & z < 4])
  zmean_tree = mean(z[z > 4])
  vert_gap = zmean_tree - zmean_shrub
  
  metrics = list(zmax = zmax, zmean = zmean, cbh = cbh, 
                 zsd = zsd, zcv = zcv, ziqr = ziqr, 
                 zskew = (sum((z - zmean)^3)/n)/(sum((z - zmean)^2)/n)^(3/2), 
                 zkurt = n * sum((z - zmean)^4)/(sum((z - zmean)^2)^2), 
                 zentropy = entropy(z, dz), pzabovezmean = pzabovemean,
                 zmean_grass = zmean_grass, zmean_shrub = zmean_shrub, 
                 zmean_tree = zmean_tree, vert_gap = vert_gap)
  metrics = c(metrics, pzabovex, zq, d, N = n)
  return(metrics)
}

.lidar_height_metrics = ~lidar_height_metrics(X, Y, Z, dz = 1, th = 2)


#============================================================
# CBH

CBH = function(z){
  z = z[z < 50 & z >= 0] # exclude outliers
  z = na.omit(z)
  zmax = max(z)
  probs = c(seq(0.05, 0.95, 0.05)) 
  zq = as.list(stats::quantile(z, probs))
  names(zq) = paste0("zq", probs * 100)
  s = as.numeric(zq)
  cbh = round(s[which.max(diff(s))+1], 2)  # after Chamberlain et al 2021
  return(list(CH = zmax, CBH = cbh))  
}
.CBH = ~CBH(Z)
cbh = ~CBH(Z)
#============================================================

lidar_cover_metrics = function(z, cl){
  n = length(z)
  z = z[cl!=2]
  cover = list(cover = length(z) / n * 100)
  
  bins = c(-Inf,0, 0.4, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 5, 6, 7, 8, 9, 10, 15, 20, 25, 30, Inf)
  zcov = as.list(cumsum(hist(z, breaks=bins, plot=F)$counts / n * 100))
  names(zcov) = c(paste0("cov_", bins[-c(1,length(bins))], "m"), "cov_30m_plus")
  
  # after Novo et al. 2020
  veg_classes = list(
    zcov_grass = zcov$cov_0.4m - zcov$cov_0m,
    zcov_shrub = zcov$cov_4m - zcov$cov_0.4m,
    zcov_tree = zcov$cov_30m_plus - zcov$cov_4m)
  
  return(c(cover, zcov, veg_classes))
}

.lidar_cover_metrics = ~lidar_cover_metrics(Z, Classification)


.rumple_index = ~rumple_index(X,Y,Z)
#============================================================

lidar_density_metrics = function(rn, cl){
  n = length(cl)
  n_ground = sum(cl == 2)
  n_canopy = sum(cl != 2)
  n_canopy_rn1 = sum(cl[rn==1] != 2)
  
  D = n_canopy_rn1 / n_canopy * 100
  pground = n_ground / n * 100
  
  return(list(D = D, pground = pground))
}

.lidar_density_metrics = ~lidar_density_metrics(ReturnNumber, Classification)

#============================================================
# find angle (°N) between 2 coordinates

angle <- function(from,to, roundby = 0){
  f = as.data.frame(sf::st_coordinates(from))
  t = as.data.frame(sf::st_coordinates(to))
  t$X = (t$X - f$X)
  t$Y = (t$Y - f$Y)
  a = atan2(t$X,t$Y) |> 
    units::as_units("radians") |> 
    units::set_units("degrees") 
  a = ((units::drop_units(a) + 360) %% 360) |> round(roundby)
  as.numeric(a)
}
#============================================================

# tree-level metrics for FuelCalc: CH, CBH
# takes a clipped, normalized, segmented PC surrounding a BWI plot
# 1. find tree tops
# 2. constrain to maximum Azimut angle from BWI plot
# 3. exclude trees with negligible share of points
# 4. return CH, CBH, and Azimut angle

tree_level_metrics = function(segPCpath, crs = 25832, radius = 12.5, limit = 0.01){
  segPC = readALSLAS(segPCpath)
  tt = segPC@data |> 
    dplyr::group_by(treeID) |> 
    dplyr::filter(Z == max(Z, na.rm=T)) |> 
    dplyr::ungroup() |> 
    dplyr::select(treeID, X,Y, Zmax=Z) |> 
    dplyr::arrange(treeID) |> 
    na.omit() |> 
    distinct(treeID, .keep_all = T) |> 
    sfheaders::sf_point("X","Y", keep = T) |> 
    sf::st_set_crs(crs)
  
  # exclude trees whose center is outside the BWI plot
  center = sf::st_bbox(segPC) |> sf::st_as_sfc() |> sf::st_centroid()
  tt$dist2c = sf::st_distance(center, tt) |> units::drop_units() |> as.vector()
  keep = which(tt$dist2c < radius)
  filt1 = filter_poi(segPC, treeID %in% keep)
  tt1 = dplyr::filter(tt, treeID %in% keep)
  
  # exclude trees that that only make up a small part of all points
  tt1$share = filt1$treeID |> table() / nrow(filt1)
  keep = tt1$treeID[which(tt1$share > limit)]
  filt2 = filter_poi(filt1, treeID %in% keep)
  tt2 = dplyr::filter(tt1, treeID %in% keep)
  
  if (! nrow(tt2) == 0){
    tm = tree_metrics(filt2, func = .CBH) |> sf::st_as_sf()
    tm$Azimut = angle(center, tm)
    tm$fname = basename(segPCpath)
    tm = dplyr::select(tm, fname, everything(), geometry)
    message(paste(basename(segPCpath), "complete."))
  } else {
    message(paste("No trees detected in", segPCpath))
    tm = NULL
  }
  return(tm)
}

####
####
# optimized for catalog_apply()

tree_level_metrics2 = function(chunk, limit = 0.01, fixed_radius = NA){
  las <- readLAS(chunk)
  if (is.empty(las)) return(NULL)
  radius = ifelse(is.na(fixed_radius), las$radius[1], fixed_radius)
  #radius = 
  if (any(all(is.na(las$treeID)), is.na(radius)) ) {
    message("no trees found :'(")
    return(NULL)
  }
  else {
    crs = epsg(object = las)
    
    segPC = data.frame(X=las$X, Y=las$Y, Z=las$Z, treeID=las$treeID) |> 
      dplyr::filter(Z < 50 & Z >= 0) # exclude outliers
    segPC = na.omit(segPC)
    segPCgroup = segPC |> dplyr::group_by(treeID)  
    tt = segPCgroup |>
      dplyr::filter(Z == max(Z, na.rm=T)) |> 
      dplyr::ungroup() |> 
      dplyr::select(treeID, X,Y, Zmax=Z) |> 
      dplyr::arrange(treeID) |> 
      na.omit() |> 
      dplyr::distinct(treeID, .keep_all = T) |> 
      sfheaders::sf_point("X","Y", keep = T) |> 
      sf::st_set_crs(crs)
    
    # exclude trees whose center is outside the BWI plot
    center = sfheaders::sf_bbox(segPC) |> sf::st_as_sfc() |> sf::st_centroid() |> sf::st_set_crs(crs)
    tt$dist2c = sf::st_distance(center, tt) |> units::drop_units() |> as.vector() |> round(1)
    keep = which(tt$dist2c < radius)
    filt1 = dplyr::filter(segPC, treeID %in% keep)
    tt1 = dplyr::filter(tt, treeID %in% keep)
    
    # exclude trees that that only make up a small part of all points
    share = filt1$treeID |> table() / nrow(filt1) 
    tt1$share = as.numeric(share) |> round(digits = 3)
    keep = tt1$treeID[which(tt1$share > limit)]
    filt2 = dplyr::filter(filt1, treeID %in% keep)
    tt2 = dplyr::filter(tt1, treeID %in% keep)
    
    if (! nrow(tt2) == 0){
      CBH = function(z){
        probs = c(seq(0.05, 0.95, 0.05)) 
        s = as.list(stats::quantile(z, probs)) |> as.numeric()
        cbh = round(s[which.max(diff(s))+1], 2)  # after Chamberlain et al 2021
        return(cbh)  
      }
      
      tt2 = segPCgroup |> 
        summarize(CBH=CBH(Z)) |> 
        inner_join(x=tt2, by="treeID")
      
      tt2$Azimut = angle(center, tt2)
      tt2$CRS = crs
      tm = tt2 |> 
        #sf::st_drop_geometry() |>                            ## remove geometry
        dplyr::select(treeID, Zmax,CBH, everything())      
      
    } else {
      tm = NULL
    }
    return(tm)
  }
}

.tree_level_metrics2 = ~tree_level_metrics2(X,Y,Z,treeID)

#------------------------------------------------------------


tree_level_metrics3 = function(las, limit = 0.01, fixed_radius = NA){
  if (is.empty(las)) return(NULL)
  radius = ifelse(is.na(fixed_radius), las$radius[1], fixed_radius)
  
  if (any(all(is.na(las$treeID)), is.na(radius)) ) {
    message("no trees found :'(")
    return(NULL)
  }
  else {
    crs = epsg(object = las)
    
    segPC = data.frame(X=las$X, Y=las$Y, Z=las$Z, treeID=las$treeID) |> 
      dplyr::filter(Z < 50 & Z >= 0) # exclude outliers
    segPC = na.omit(segPC)
    segPCgroup = segPC |> dplyr::group_by(treeID)  
    tt = segPCgroup |>
      dplyr::filter(Z == max(Z, na.rm=T)) |> 
      dplyr::ungroup() |> 
      dplyr::select(treeID, X,Y, Zmax=Z) |> 
      dplyr::arrange(treeID) |> 
      na.omit() |> 
      dplyr::distinct(treeID, .keep_all = T) |> 
      sfheaders::sf_point("X","Y", keep = T) |> 
      sf::st_set_crs(crs)
    
    # exclude trees whose center is outside the BWI plot
    center = sfheaders::sf_bbox(segPC) |> sf::st_as_sfc() |> sf::st_centroid() |> sf::st_set_crs(crs)
    tt$dist2c = sf::st_distance(center, tt) |> units::drop_units() |> as.vector() |> round(1)
    keep = which(tt$dist2c < radius)
    filt1 = dplyr::filter(segPC, treeID %in% keep)
    tt1 = dplyr::filter(tt, treeID %in% keep)
    
    # exclude trees that only make up a small part of all points
    share = filt1$treeID |> table() / nrow(filt1) 
    tt1$share = as.numeric(share) |> round(digits = 3)
    keep = tt1$treeID[which(tt1$share > limit)]
    filt2 = dplyr::filter(filt1, treeID %in% keep)
    tt2 = dplyr::filter(tt1, treeID %in% keep)
    
    if (! nrow(tt2) == 0){
      CBH = function(z){
        probs = c(seq(0.05, 0.95, 0.05)) 
        s = as.list(stats::quantile(z, probs)) |> as.numeric()
        cbh = round(s[which.max(diff(s))+1], 2)  # after Chamberlain et al 2021
        return(cbh)  
      }
      
      tt2 = segPCgroup |> 
        summarize(CBH=CBH(Z)) |> 
        inner_join(x=tt2, by="treeID")
      
      tt2$Azimut = angle(center, tt2)
      tt2$CRS = crs
      tm = tt2 |> 
        #sf::st_drop_geometry() |>                            ## remove geometry
        dplyr::select(treeID, Zmax,CBH, everything())      
      
    } else {
      tm = NULL
    }
    return(tm)
  }
}


#===========================================================
# tree list metrics from point cloud
process_tlm = function(i, gedi, .ctg = ctg, outdir, rad = 12.5){
  out = file.path(outdir, paste0("met_", gedi$state[i], "_", gedi$uniqueID[i], ".gpkg"))
  if (!file.exists(out)){
    plt = clip_roi(.ctg, st_buffer(gedi[i,], rad))
    if (is.empty(plt)) return(NULL)
    plt = normalize_height(plt, tin()) |> 
      segment_trees(li2012(R = 0, hmin=2)) |> 
      tree_level_metrics3(fixed_radius = rad)
    if (!is.null(plt)) sf::st_write(plt, out, quiet=T, append=F)
  }
}

#===========================================================
# tree list metrics from point cloud with DTM kriging
process_tlm2 = function(i, gedi, .ctg = ctg, mdir, alsdir, rad = 12.5){
  out = paste0("uID_", gedi$uniqueID[i])
  met_out = file.path(mdir, paste0(out, ".gpkg"))
  if (!file.exists(met_out)){
    plt = clip_roi(.ctg, st_buffer(gedi[i,], rad))
    if (is.empty(plt)) return(NULL)
    
    plt = normalize_height(plt, kriging(k=40)) |> 
      segment_trees(li2012(R = 0, hmin=2)) 
    if (!is.null(plt)){
      writeLAS(plt, file.path(alsdir, paste0(out, ".las")))
      
      met = plt |> tree_level_metrics3(fixed_radius = rad)
      if (!is.null(met)) sf::st_write(met, met_out, quiet=T, append=F)
    }
  }
}

process_tlm_by_tile = function(tile_path, gedi, mdir, alsdir, rad = 12.5){
  
  gedi = dplyr::filter(gedi, tile == tile_path)
  message(basename(tile_path),": ", nrow(gedi), " obs")
  out = paste0("id_", gedi$uniqueID,"_tile_", 
               stringr::str_sub(basename(gedi$tile[1]), 1,-5))
  met_out = file.path(mdir, paste0(out, ".gpkg"))
  
  if (!all(file.exists(met_out))){
    las = readLAS(tile_path)
    st_crs(las) = st_crs(gedi)
    
    plt = clip_roi(las, st_buffer(gedi, rad))
    if (! is.list(plt)) plt = list(plt)
    
    plt = purrr::map(plt, .f = function(x){
      if (!is.empty(x)){ 
        normalize_height(x, kriging(k=40)) |> 
        segment_trees(li2012(R = 0, hmin=2))
    }})
    
    p = purrr::map(1:length(plt), .f = function(i){
      if (!is.null(plt[[i]])){ 
        writeLAS(plt[[i]], file.path(alsdir, paste0(out[i], ".las")))
      }})
    
    p = purrr::map(1:length(plt), .f = function(i){
      met = tree_level_metrics3(plt[[i]], fixed_radius = rad)
      if (!is.null(met)) sf::st_write(met, met_out[i], quiet=T, append=F)
    })
  }
}


#============================================================

bwi_metrics = function(las, lasname, Zname = "Z_norm_direct_tin", write=T, writedir = getwd()){
  
  if (is.empty(las)) return(NULL)
  if (all(is.na(las$treeID))) return(NULL)
  else {
    segPC = select(las@data,X,Y,Z = dplyr::all_of(Zname), treeID) |> 
      dplyr::filter(Z < 55 & Z >= 0)  |> 
      na.omit()
    segPCgroup = segPC |> dplyr::group_by(treeID)  
    tt = segPCgroup |>
      dplyr::filter(Z == max(Z, na.rm=T)) |> 
      dplyr::ungroup() |> 
      dplyr::select(treeID, X,Y, Zmax=Z) |> 
      dplyr::arrange(treeID) |> 
      na.omit() |> 
      distinct(treeID, .keep_all = T) |> 
      sfheaders::sf_point("X","Y", keep = T) 
    
    # exclude trees whose center is outside the BWI plot
    center = sfheaders::sf_bbox(segPC) |> 
      sf::st_as_sfc() |> 
      sf::st_centroid() 
    tt$dist2c = sf::st_distance(center, tt) |> 
      as.vector() |> round(1)
    
    # exclude trees that that only make up a small part of all points
    share = segPC$treeID |> table() / nrow(segPC) 
    tt$share = as.numeric(share) |> round(digits = 4)
    
    if (! nrow(tt) == 0){
      CBH = function(z){
        probs = c(seq(0.05, 0.95, 0.05)) 
        s = as.list(stats::quantile(z, probs)) |> as.numeric()
        cbh = round(s[which.max(diff(s))+1], 2)  # after Chamberlain et al 2021
        return(cbh)  
      }
      
      tt2 = segPCgroup |> 
        summarize(CBH = CBH(Z)) |> 
        inner_join(x = tt, by = "treeID")
      tt2$Azimut = angle(center, tt)
      
      name = strsplit(lasname, "_") 
      tt2$Tnr = name[[1]][1]|> as.numeric()
      tt2$Enr = name[[1]][2]|> as.numeric()
      
      tm = tt2 |> dplyr::select(Tnr, Enr, treeID, Zmax,CBH, everything())      
      
    } else {
      tm = NULL
    }
    if (write) write.csv(tm, file.path(writedir, paste0(lasname, ".csv"))) else return(tm)
  }
}


#============================================================

# create sf hole-polygons from two unequal buffers around points
make_rings = function(points, inner_buffer, outer_buffer){
  stopifnot(all(st_is(plots, "POINT")))
  require(dplyr)
  circ.in = sf::st_buffer(sf::st_geometry(points), inner_buffer)
  circ.out = sf::st_buffer(sf::st_geometry(points), outer_buffer)
  
  ring = lapply(1:nrow(points), 
                function(x) sf::st_sf(sf::st_difference(circ.out[x], 
                                                        circ.in[x]))) %>% 
    do.call(rbind, .) %>% 
    setNames('geometry') %>% 
    sf::st_set_geometry('geometry')
  return(ring)
}

#===========================================================
# clip-normalize-segment routine

all_inc_segmentation = function(ctg, points, radius = 12.5){
  
  las = clip_roi(las, points, radius = radius)
  #if (is.empty(las)) las = NULL
  norm = normalize_height(las, tin())
  seg = segment_trees(ctg, li2012(R = 0, hmin=2))
  return(seg)
}



