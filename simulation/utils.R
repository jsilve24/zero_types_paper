plot_t1 <- function(d, x.lim=NULL, scales = "free_y"){
  cn <- colnames(d)
  d <- d %>% 
    gather(f, dist)
  
  
  d.summary <- d %>%
    group_by(f) %>%
    summarize(p2.5 = quantile(dist, p=0.025), 
              p25 = quantile(dist, p=0.25),
              mean = mean(dist),
              p75 = quantile(dist, p=0.75), 
              p97.5 = quantile(dist, p=0.975)) %>% 
    as.data.frame()
  rownames(d.summary) <- d.summary$f
  
  d.density <- d %>% 
    split(.$f) %>% 
    map(~with(density(.x$dist), data.frame(x,y))) %>% 
    map(as.data.frame) %>% 
    bind_rows(.id="f")  %>% 
    mutate(gt.p25 = x >= d.summary[f,"p25"],
           lt.p75 = x <= d.summary[f,"p75"], 
           include.25.75 = gt.p25 & lt.p75, 
           gt.p2.5 = x >= d.summary[f,"p2.5"],
           lt.p97.5 = x <= d.summary[f,"p97.5"], 
           include.2.5.97.5 = gt.p2.5 & lt.p97.5, 
           mean = d.summary[f, "mean"]) %>% 
    mutate(f = factor(f, levels = cn)) 
  
  tmp <- rep(0, nrow(d.summary))
  for (i in 1:nrow(d.summary)){
    incl <- d.density$f == d.summary[i,"f"]
    tmp[i] <- d.density$y[incl][which.min(abs(d.density$x[incl]-d.summary[i,"mean"]))]
  }
  d.summary$yend <- tmp
  
  p <- d.density %>%
    ggplot(aes(x = x, y= y))+
    geom_area(fill="white", alpha=0.3, color="black") +
    geom_area(data = subset(d.density, include.2.5.97.5==TRUE), fill="#619CFF", alpha=0.5) +
    geom_area(data = subset(d.density, include.25.75==TRUE), fill="#619CFF", alpha=0.8) + 
    geom_segment(data = d.summary, aes(x = mean, y=0, xend = mean, yend=yend), color="darkblue") +
    facet_grid(f~., scales=scales) +
    theme_bw()+
    theme(strip.text.y = element_text(angle = 0)) +
    ylab("Density")+
    xlab("x")
  if (!is.null(x.lim)) p <- p + xlim(x.lim)
  p
}


plot_t3 <- function(d,group_factor, x.lim=NULL, scales = "free_y"){
  cn <- colnames(d)
  d <- d %>% 
    gather(f, dist, -"person")
  
  
  d.summary <- d %>%
    group_by(f, person) %>%
    summarize(p2.5 = quantile(dist, p=0.025), 
              p25 = quantile(dist, p=0.25),
              mean = mean(dist),
              p75 = quantile(dist, p=0.75), 
              p97.5 = quantile(dist, p=0.975)) %>% 
    as.data.frame()
  rownames(d.summary) <- paste0(d.summary$f, "_", d.summary$person)
  
  d %>% 
    split(list(.$f, .$person), sep="_") %>% 
    map(~density(.x$dist)) -> foo
  
  d.density <- d %>% 
    split(list(.$f, .$person), sep="_") %>% 
    map(~with(density(.x$dist, kernel="biweight"), data.frame(x,y))) %>% 
    map(as.data.frame) %>% 
    bind_rows(.id="f")  %>% 
    mutate(gt.p25 = x >= d.summary[f,"p25"],
           lt.p75 = x <= d.summary[f,"p75"], 
           include.25.75 = gt.p25 & lt.p75, 
           gt.p2.5 = x >= d.summary[f,"p2.5"],
           lt.p97.5 = x <= d.summary[f,"p97.5"], 
           include.2.5.97.5 = gt.p2.5 & lt.p97.5, 
           mean = d.summary[f, "mean"]) %>% 
    mutate(f_comb = f) %>% 
    separate(f, c("f", "person")) %>% 
    mutate(f = factor(f, levels = cn), 
           person = factor(person, levels = 1:length(unique(person))))
  
  tmp <- rep(0, nrow(d.summary))
  for (i in 1:nrow(d.summary)){
    incl <- d.density$f_comb == rownames(d.summary)[i]
    tmp[i] <- d.density$y[incl][which.min(abs(d.density$x[incl]-d.summary[i,"mean"]))]
  }
  d.summary$yend <- tmp
  
  p <- d.density %>%
    ggplot(aes(x = x, y= y))+
    geom_area(fill="white", alpha=0.3, color="black") +
    geom_area(data = subset(d.density, include.2.5.97.5==TRUE), fill="#619CFF", alpha=0.5) +
    geom_area(data = subset(d.density, include.25.75==TRUE), fill="#619CFF", alpha=0.8) + 
    geom_segment(data = d.summary, aes(x = mean, y=0, xend = mean, yend=yend), color="darkblue") +
    facet_grid(f~person, scales=scales) +
    theme_bw()+
    theme(strip.text.y = element_text(angle = 0)) +
    ylab("Density")+
    xlab("x")
  if (!is.null(x.lim)) p <- p + xlim(x.lim)
  p
}

