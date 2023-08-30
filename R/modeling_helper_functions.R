codes = read.csv("data/Raster/tree_species_code.csv")
plot_fractions = function(pl){
  uid = pl$uniqueID[1]
  x = #inner_join(pl, codes, by=c("value"="ID")) |> 
    right_join(pl, codes,  by=c("value"="ID")) |> 
    mutate(uniqueID = uid) |> 
    group_by(name) |> 
    summarise(species = paste0("frac_", unique(name)), 
              coverage = round(sum(coverage_fraction) / 4.908739, 3),
              uniqueID = unique(uniqueID)) |>     # = (pi*(12.5 ** 2)/100)
    select(-name) |> 
    pivot_wider(names_from = "species", values_from = "coverage") 
  return(x)
}

par_plot_frac = function(g, x){
  p = x[g]
  p = purrr::map_dfr(p, plot_fractions)
  replace(p, is.na(p), 0)
}

library(ggplot2)
varImpPlot_ranger = function(x, var, top = 15){
  df = data.frame(variable = names(x$variable.importance), 
                  importance = x$variable.importance) |> 
    arrange(-importance) |> head(top)
  ggplot(df, aes(x=reorder(variable,importance), y=importance,fill=importance))+ 
    geom_bar(stat="identity", position="dodge")+ coord_flip()+
    ylab("Variable Importance")+
    xlab("")+
    ggtitle(var)+
    guides(fill='none')+
    scale_fill_gradient(low="red", high="blue")
}

varImpPlot_df = function(x, var, top = 15){
  df = data.frame(variable = rownames(x), 
                  importance = x[,var]) |> 
    arrange(-importance) |> head(top)
  ggplot(df, aes(x=reorder(variable,importance), y=importance,fill=importance))+ 
    geom_bar(stat="identity", position="dodge", fill="grey30", color = "white")+ coord_flip()+
    ylab("VI")+ xlab("") + ggtitle(var)+ theme_minimal() +
    scale_y_continuous(limits = c(0,100)) +
    theme(axis.text = element_text(size = 10, color=1), plot.margin = rep(unit(0,"cm"),4),
          panel.grid.major.y = element_blank()) +
    guides(fill='none')#+scale_fill_gradient(low="red", high="blue") theme_light()
}


get_measures = function(pred, test){
  m = caret::postResample(pred |> sf::st_drop_geometry(),
                          test |> sf::st_drop_geometry())
  relRMSE = c("rel" = m[1] / sd(test, na.rm=T))
  return(c(m, relRMSE) |> round(2))
}

scatter = function(df, x, y, facets = NULL, ...){
  sp <- ggpubr::ggscatter(df, x = x, y = y, alpha=0.2, add = "reg.line",
                  add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                  conf.int = TRUE, ...)
  if (!is.null(facets)) sp = ggpubr::facet(sp, facet.by = facets)
  sp + ggpubr::stat_cor(p.accuracy = 0.001, r.accuracy = 0.01)
}

make_even_groups = function(n, ngroups){
  s = seq(0, n, length.out = ngroups+1) |> as.integer()
  lapply(1:ngroups, \(x) (s[x]+1):(s[x+1]))
}

merge_aoa = function(aoa_list){
  A = aoa_list[[1]]
  for (i in 2:length(aoa_list)){
    A$DI = c(A$DI, aoa_list[[i]]$DI)
    A$AOA = c(A$AOA, aoa_list[[i]]$AOA)
    A$uniqueID = c(A$uniqueID, aoa_list[[i]]$uniqueID)
  }
  return(A)
}

# mini helper: counts per state
s_table = function(x, bb="bb|als_", nrw="nrw|nw",sn="sn|lsc",th="th", total=T){
  stopifnot(class(x) == 'character')
  t = c(
    "BB" = grep(bb,x, ignore.case = T) |> length(),
    "NRW" = grep(nrw,x, ignore.case = T) |> length(),
    "SN" = grep(sn,x, ignore.case = T) |> length(),
    "TH" = grep(th,x, ignore.case = T) |> length()
  )
  t = c(t, "BY" = length(x) - sum(t))
  if (total) t = c(t, "total" = sum(t))
  return(t)
}
