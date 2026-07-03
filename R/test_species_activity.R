
jphenos = c(input$sflw_phenos,
                      input$slf_phenos)

n_phenos = length(jphenos)

stopifnot("species activity curve can only be shown with one or two phenophases selected" = n_phenos <= 2)

sphenos = pheno_df |> 
  sbt(full_nm %iin% jphenos) |> 
  getElement('cln_nm') 

pd = d |> 
  sbt(yr >= input$syr_rng[1] & yr <= input$syr_rng[2]) |> 
  sbt(species_2 %==% input$species) |> 
  slt(c("species_2", "date", "yr", "wk", sphenos)) |> 
  frename(sp = "species_2") |> 
  pivot(ids = 1:4) |> 
  mtt(obs_ind = fctr(c("not observed", "observed")[value+1])) |> 
  roworder(obs_ind)

z_df = expand.grid(yr = input$yr_rng[1]:input$yr_rng[2],
                   wk = 1:52) |>
  mtt(prop = 0, n = 0, k = 0) |>
  sbt(!(yr == latest_wk$yr & wk > latest_wk$wk)) |>
  mtt(from_zf = TRUE) |> 
  qDT()

obs_cnts = pd |> 
  gby(yr, wk, variable) |> 
  smr(prop = fmean(value),
      k = fsum(value),
      n = fnobs(value)) |>
  mtt(from_zf = FALSE)

z_to_use = join(z_df, obs_cnts, 
                on = c("yr", "wk"), 
                how = "anti",
                verbose = FALSE)

z_by_pheno = lapply(sphenos, \(x) {
  z_to_use |> mtt(variable = x)
}) |> 
  rbindlist()

zf_cnts = z_by_pheno |> 
  rbind(obs_cnts) |> 
  roworder(yr, wk) |> 
  qDT() |> 
  mtt(nfail = n - k,
      p = k / (k + nfail),
      qinit = qlogis((k+1) / (n+2))) 


fit_curve = function(zf) {
  agg_cnts = zf |> 
    gby(wk) |> 
    smr(ntot = fsum(n),
        nk = fsum(k)) |> 
    mtt(qinit = qlogis((nk + .3) / (ntot + .6)))
  
  init_fun = \() {
    pvec = c(1,
             fsd(agg_cnts$qinit),
             fmean(agg_cnts$qinit), 
             W(agg_cnts$qinit),
             alloc(0, fnrow(zf)))
    
    list(ell = pvec[1],
         wk_sigma = pvec[2],
         intercept = pvec[3],
         wkv = pvec[1:52 + 3],
         yrwkv = alloc(0, fnrow(zf)))
  }
  
  dl = list(N = fnrow(zf),
            k = zf$k,
            n = zf$n,
            wk_i = zf$wk)
  
  fit = optimizing(m, 
                   dl, 
                   init = init_fun)
  
  ell = fit$par[1]
  wk_sigma = fit$par[2]
  intercept = fit$par[3]
  
  eff = tail(fit$par, -3) 
  wkv_fit = head(eff, 52)
  yrwkv_fit = tail(eff, -52) 
  
  zf |> 
    mtt(fitted = plogis(wkv_fit[zf$wk] + yrwkv_fit + intercept)) |> 
    slt(yr, wk, fitted)
}

fit_by_pheno = zf_cnts[,.(fit_params = .(fit_curve(.SD))), by = "variable"][,rbindlist(fit_params), by = variable]

dodge = 1

yr_diff = max(zf_cnts$yr) - min(zf_cnts$yr)

zf_cnts = zf_cnts |> 
  join(fit_by_pheno, on = c("variable", "yr", "wk")) |> 
  mtt(d_hide = as.Date("2026-01-01") + 7*wk - 3.5 + 
      dodge * (yr - 2023) - (yr_diff*dodge/2),
      pheno_yr = paste(variable, yr, sep = "-"))

get_pal = \(pal_f, max_i) {
  col_dt = data.table(yr = 2023:latest_wk$yr) |> 
    mtt(i = floor(seq(1, max_i, length.out = diff(range(yr)) + 1)),
        col = pal_f(100)[i])
  
  col_dt$col
}

pal_l = list(pals::parula, pals::magma)

cmap_df = expand.grid(yr = seq(input$syr_rng[1],
                               input$syr_rng[2]),
                      ph = sphenos) |> 
  mtt(ph_yr = paste(ph, yr, sep = '-'),
      col_val = as.vector(mapply(FUN = get_pal, pal_l, c(88, 90))))

cmap_v = cmap_df$col_val
names(cmap_v) = cmap_df$ph_yr

zf_cnts |> 
  ggplot(aes(d_hide, p)) + 
  geom_line(aes(color = pheno_yr, y = fitted),
            lwd = .7) + 
  geom_point(aes(color = pheno_yr)) + 
  scale_color_manual(values = cmap_v) + 
  theme_dark() + 
  labs(x = NULL, y = 'proportion',
       title = paste0("*", input$species, "*"),
       color = NULL) + 
  scale_x_date(labels = scales::label_date("%b"),
               breaks = as.Date(paste0("2026-", 1:12, "-15"))) + 
  theme(panel.grid.minor.x = element_blank(),
        axis.title.y = element_text(margin = margin(0,0,0,0, 'pt'),
                                    vjust = -13),
        plot.title = element_markdown(),
        text = element_text(size=18)) 
