# This version uses rstan to update the lines

library(shiny)
library(ggplot2)
library(data.table)
library(collapse)
library(ggtext)
library(forcats)

source("load_rstan.R")

to_fct = \(x, genus) {
  x |> 
    fctr() |> 
    fct_relabel(\(x) paste0("*", x, "*")) |> 
    fct_lump_min(5, other_level = paste0("other *", genus, "*")) |> 
    fct_infreq() |> 
    fct_rev() 
}

print("reading...")
# d = fread("data/anecdata_export_EwA_Pheno_Lite_2026-05-28T20-24-36-075Z.csv") |> 
#   janitor::clean_names()

d = fread("data/adj_2026-06-26T23-02-37-068Z.csv")

uniq_lf = fread("data/uniq_lf.tsv") 

uniq_flw = fread('data/uniq_flw.tsv')

latest_wk = d |> 
  slt(wk, yr) |> 
  funique() |> 
  sbt(yr %==% fmax(yr)) |> 
  sbt(wk %==% fmax(wk))

genera = d |> 
  slt(genus) |> 
  fcount() |> 
  roworder(-N) |> 
  na_omit() |> 
  sbt(N >= 30) # can select from genera with >= 30 observations

gs_df = fread("data/Genus_species.tsv")

pheno_df = rowbind(
 list(
  uniq_flw |> slt(full_nm = flw_pheno, cln_nm = cln_flw),
  uniq_lf |> slt(full_nm = lf_pheno, cln_nm = cln_lf)
 ) 
)

input = list(genus = "Quercus", flw_phenos = "Flowers")
input$species = "Quercus rubra"
input$lf_phenos = "Leaves"
input$log_fun = "AND"
input$sflw_phenos = "Flowers"
input$slf_phenos = "Leaves"
input$yr_rng = c(2023, 2026)
input$syr_rng = c(2023, 2026)

# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("Phenolite"),
    
    tabsetPanel(
      tabPanel("Genus view", 
               sidebarPanel(selectInput("genus", 
                                        label = h4("Select Genus:"),
                                        choices = genera$genus),
                            selectInput('flw_phenos',
                                        label = h4("Select flower phenophase:"),
                                        choices = uniq_flw$flw_pheno,
                                        selected = 'Flowers',
                                        multiple = TRUE),
                            selectInput('lf_phenos',
                                        label = h4("Select leaf phenophase:"),
                                        choices = uniq_lf$lf_pheno,
                                        selected = "Leaves",
                                        multiple = TRUE),
                            selectInput('log_fun',
                                        label = h5("Combine phenophase selection with:"),
                                        choices = c("AND", "OR"),
                                        selected = "AND",
                                        selectize = FALSE),
                            sliderInput('yr_rng', 
                                        label = h5("Year range:"),
                                        min = 2023, 
                                        max = latest_wk$yr,
                                        value = c(2023, 2026),
                                        step = 1,
                                        ticks = TRUE,
                                        sep = "")),
               mainPanel(plotOutput('combinedPlot',
                                    height = '800px'))
               
      ),
      tabPanel("Species view",
               sidebarPanel(selectInput("species",
                                        label = h4("Select Species:"),
                                        choices = gs_df$binom,
                                        selected = "Quercus rubra"),
                            selectInput('sflw_phenos',
                                        label = h4("Select flower phenophase:"),
                                        choices = uniq_flw$flw_pheno,
                                        selected = 'Flowers',
                                        multiple = TRUE),
                            selectInput('slf_phenos',
                                        label = h4("Select leaf phenophase:"),
                                        choices = uniq_lf$lf_pheno,
                                        selected = "Leaves",
                                        multiple = TRUE),
                            sliderInput('syr_rng', 
                                        label = h5("Year range:"),
                                        min = 2023, 
                                        max = latest_wk$yr,
                                        value = c(2023, 2026),
                                        step = 1,
                                        ticks = TRUE,
                                        sep = "")),
               mainPanel(plotOutput('species_activity',
                                    height = '300px'),
                         plotOutput('species_plot', 
                                    height = '500px')))
               
    ),
)


# Define server logic required to draw a histogram
server <- function(input, output) {
  
  get_species_act = function() {
    
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
            axis.title.y = element_text(margin = margin(0,0,0,0, 'pt')),
            legend.byrow = TRUE,
            legend.position = 'top',
            plot.title = element_markdown(),
            text = element_text(size=18)) 
  }
  
  output$species_activity <- renderPlot(get_species_act())
  
  output$combinedPlot <- renderPlot({ 
    
    
    sel_obs = d |>
      sbt(yr >= input$yr_rng[1] & yr <= input$yr_rng[2]) |> 
      sbt(genus %like% input$genus)
    
    logical_fun = if (input$log_fun == "AND") {
      matrixStats::rowProds 
    } else {
      matrixStats::rowSums2 
    }
   
    if (!is.null(input$flw_phenos)) {
      
      flw_cols = uniq_flw |> 
        sbt(flw_pheno %iin% input$flw_phenos) |> 
        get_elem("cln_flw")
      
      flw_mat = sel_obs |> 
        slt(flw_cols) |> 
        qM() 
      
    } else {
      if (input$log_fun == "AND") {
        flw_mat = matrix(TRUE, nrow = fnrow(sel_obs))
      } else {
        flw_mat = matrix(FALSE, nrow = fnrow(sel_obs))
      }
    }
    
    if (!is.null(input$lf_phenos)) {
       lf_cols = uniq_lf |> 
        sbt(lf_pheno %iin% input$lf_phenos) |> 
        get_elem("cln_lf")
      
      lf_mat = sel_obs |> 
        slt(lf_cols) |> 
        qM()
      
    } else {
      if (input$log_fun == "AND") {
        lf_mat = matrix(TRUE, nrow = fnrow(sel_obs))
      } else {
        lf_mat = matrix(FALSE, nrow = fnrow(sel_obs))
      }
    }
    
    cnd_vec = cbind(flw_mat, lf_mat) |> 
      logical_fun() 
    
    sel_obs = sel_obs |> 
      mtt(meets_cnd = as.logical(cnd_vec)) 
    
    z_df = expand.grid(yr = input$yr_rng[1]:input$yr_rng[2],
                       wk = 1:52) |>
      mtt(prop = 0, n = 0, k = 0) |>
      sbt(!(yr == latest_wk$yr & wk > latest_wk$wk)) |>
      mtt(from_zf = TRUE) |> 
      qDT()
    
    obs_cnts = sel_obs |> 
      gby(yr, wk) |> 
      smr(prop = fmean(meets_cnd),
          k = fsum(meets_cnd),
          n = fnobs(meets_cnd)) |> 
      mtt(from_zf = FALSE)
    
    zf_cnts = join(z_df, obs_cnts, 
                   on = c("yr", "wk"), 
                   how = "anti",
                   verbose = FALSE) |> 
      rbind(obs_cnts) |> 
      roworder(yr, wk) |> 
      qDT() |> 
      mtt(nfail = n - k,
          p = k / (k + nfail),
          qinit = qlogis((k+1) / (n+2))) 
    
    agg_cnts = zf_cnts |> 
      gby(wk) |> 
      smr(ntot = fsum(n),
          nk = fsum(k)) |> 
      mtt(qinit = qlogis((nk + .3) / (ntot + .6)))
    
    init_fun = \() {
      pvec = c(1,
               fsd(agg_cnts$qinit),
               fmean(agg_cnts$qinit), 
               W(agg_cnts$qinit),
               alloc(0, fnrow(zf_cnts)))
      
      list(ell = pvec[1],
           wk_sigma = pvec[2],
           intercept = pvec[3],
           wkv = pvec[1:52 + 3],
           yrwkv = alloc(0, fnrow(zf_cnts)))
    }

    dl = list(N = fnrow(zf_cnts),
              k = zf_cnts$k,
              n = zf_cnts$n,
              wk_i = zf_cnts$wk)
    
    fit = optimizing(m, 
                     dl, 
                     init = init_fun)
    
    ell = fit$par[1]
    wk_sigma = fit$par[2]
    intercept = fit$par[3]
    
    eff = tail(fit$par, -3) 
    wkv_fit = head(eff, 52)
    yrwkv_fit = tail(eff, -52) 
    
    zf_cnts = zf_cnts |> 
      mtt(fitted = plogis(wkv_fit[zf_cnts$wk] + yrwkv_fit + intercept))
    
    # avg_df = data.table(
    #   d = as.Date(as.IDate("2026-01-01") + 7*1:52 - 3.5),
    #   fitted = plogis(wkv_fit + intercept)
    # )
    
    dodge = 1
    
    yr_diff = max(zf_cnts$yr) - min(zf_cnts$yr)
    
    plot_input = zf_cnts |> 
      mtt(d = as.IDate("2026-01-01") + 7*wk - 3.5 + 
            dodge * (yr - 2023) - (yr_diff*dodge/2)) |> 
      mtt(yr = factor(yr)) 
    
    tstring = paste0('*', input$genus, "* - ", 
                     paste0(unlist(c(input$flw_phenos, input$lf_phenos)), 
                            collapse = paste0(' ', input$log_fun, ' ')))
    
    # complicated method to get the colors so they don't change when your change
    # the year range
    col_dt = data.table(yr = 2023:latest_wk$yr) |> 
      mtt(i = floor(seq(1, 88, length.out = diff(range(yr)) + 1)),
          col = pals::parula(100)[i])
    
    col_vec = col_dt$col |> setNames(col_dt$yr)
    

    # top plot ----------------------------------------------------------------

    
    p1 = plot_input |> 
      ggplot(aes(d, prop)) + 
      geom_line(aes(y = fitted, color = yr)) +
      geom_point(aes(color = yr, group = yr),
                 pch = 15) + 
      # geom_line(data = avg_df,
      #           aes(y = fitted), color = 'red') +
      scale_color_manual(values = col_vec) + 
      scale_x_date(labels = scales::label_date("%b"),
                   breaks = as.Date(paste0("2026-", 1:12, "-15"))) + 
      ylim(c(0,1)) + 
      theme_bw() + 
      theme(panel.grid.minor.x = element_blank(),
            axis.title.y = element_text(margin = margin(0,0,0,0, 'pt'),
                                        vjust = -13),
            plot.title = element_markdown(),
            text = element_text(size=18)) + 
      labs(x = NULL,
           y = 'proportion',
           color = NULL,
           title = tstring)
    
    # Start to differ from above starting here. 
    
    obs_df = sel_obs |> 
      mtt(gs_fct = to_fct(species_2, input$genus), 
          d_hide = as.Date("2026-01-01") + 7*wk - 3.5,
          obs_ind = c("not observed", "observed")[meets_cnd + 1]) |> 
      roworder(meets_cnd)
    
    tstring = paste0('*', input$genus, "* - ", 
                     paste0(unlist(c(input$flw_phenos, input$lf_phenos)), 
                            collapse = paste0(' ', input$log_fun, ' ')))
    
    # bottom plot ---- 
    p2 = ggplot(obs_df, aes(d_hide, gs_fct)) + 
      geom_point(pch = 15, 
                 aes(color = obs_ind),
                 position = position_jitter(width=.5, height = .2),
                 size = .6) + 
      scale_color_manual(values = c("grey", "black")) + 
      labs(color = NULL,
           x = NULL,
           y = NULL) + 
      scale_x_date(labels = scales::label_date("%b"),
                   breaks = as.Date(paste0("2026-", 1:12, "-15"))) + 
      theme_bw() + 
      theme(panel.grid.minor.x = element_blank(),
            plot.margin = unit(c(.3,.3,.3,1), "cm"),
            axis.text.y = element_markdown(),
            plot.title = element_markdown(),
            text = element_text(size=18)) 
    
    g1 = ggplotGrob(p1)
    g2 = ggplotGrob(p2)
    
    grid::grid.newpage()
    grid::grid.draw(rbind(g1, g2))
  })
  
  output$species_plot <- renderPlot({
    
    sphenos = pheno_df |> 
      sbt(full_nm %iin% c(input$sflw_phenos,
                          input$slf_phenos)) |> 
      getElement('cln_nm') 
    
    pd = d |> 
      sbt(yr >= input$syr_rng[1] & yr <= input$syr_rng[2]) |> 
      sbt(species_2 %==% input$species) |> 
      slt(c("species_2", "date", "yr", "wk", sphenos)) |> 
      frename(sp = "species_2") |> 
      pivot(ids = 1:4) |> 
      mtt(obs_ind = fctr(c("not observed", "observed")[value+1])) |> 
      roworder(obs_ind)
    
    lubridate::year(pd$date) <- 2026
    
    pd |> 
      ggplot(aes(date, variable)) + 
      geom_jitter(aes(color = obs_ind),
                  width = .5,
                  height = .2,
                  size = .8,
                  pch = 15) + 
      facet_wrap(ncol = 1, 
                 vars(yr)) + 
      labs(y = NULL,
           title = NULL,
           x = NULL,
           color = NULL) + 
      theme_bw() + 
      theme(plot.title = element_markdown(),
            panel.grid.major.y = element_blank(),
            panel.spacing = unit(1, "pt"),
            legend.position = 'bottom',
            text = element_text(size=18),
            strip.text = element_text(margin = margin(0,0,0,0))) + 
      scale_color_manual(values = c("grey", "black")) + 
      scale_x_date(labels = scales::label_date("%b"),
                   breaks = as.Date(paste0("2026-", 1:12, "-15"))) 
  })
}

# Run the application 
shinyApp(ui = ui, server = server)
