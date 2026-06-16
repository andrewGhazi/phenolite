#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#

library(shiny)
library(ggplot2)
library(data.table)
library(collapse)
library(ggtext)

print("compiling...")
source("binom_gp_ll.R")

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

d = fread("data/adj_2026-05-28T20-24-36-075Z.csv")

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

# input = list(genus = "Quercus", flw_phenos = "Flowers")
# input$lf_phenos = "Leaves"
# input$log_fun = "AND"

# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("Phenolite"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
          selectInput("genus", 
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
                      label = h5("Combine conditions with:"),
                      choices = c("AND", "OR"),
                      selected = "AND")
        ),

        # Show a plot of the generated distribution
        mainPanel(
          plotOutput('combinedPlot',
                     height = '800px')
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  
  output$combinedPlot <- renderPlot({ 
    
    print('subsetting...')
    
    sel_obs = d |>
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
      flw_mat = matrix(TRUE, nrow = fnrow(sel_obs))
    }
    
    if (!is.null(input$lf_phenos)) {
       lf_cols = uniq_lf |> 
        sbt(lf_pheno %iin% input$lf_phenos) |> 
        get_elem("cln_lf")
      
      lf_mat = sel_obs |> 
        slt(lf_cols) |> 
        qM()
      
    } else {
      lf_mat = matrix(TRUE, nrow = fnrow(sel_obs))
    }
    
    cnd_vec = cbind(flw_mat, lf_mat) |> 
      logical_fun() 
    
    sel_obs = sel_obs |> 
      mtt(meets_cnd = as.logical(cnd_vec)) 
    
    z_df = expand.grid(yr = 2023:latest_wk$yr,
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
    
    # m = glmer(family = "binomial",
    #           cbind(zf_cnts$k+1, zf_cnts$nfail+1) ~ (1|wk:yr),
    #           data = zf_cnts)
    
    pvec = c(1,
             fsd(agg_cnts$qinit),
             fmean(agg_cnts$qinit), 
             W(agg_cnts$qinit),
             alloc(0, fnrow(zf_cnts)))
    
    print("optimizing...")
    
    opt_res = optim(pvec, 
                    fn = lpost2, 
                    N = fnrow(zf_cnts),
                    k = zf_cnts$k,
                    n = zf_cnts$n,
                    wk_i = zf_cnts$wk-1,
                     control = list(fnscale = -1,
                                    maxit = 1000),
                     lower = c(.3,.1, alloc(-Inf, 1+52+178)),
                     method = "L-BFGS-B")
    
    eff = opt_res$par |> tail(-3)
    wkv_fit = eff |> head(52)
    yrwkv_fit = eff |> tail(-52)
    
    zf_cnts = zf_cnts |> 
      mtt(fitted = plogis(wkv_fit[zf_cnts$wk] + yrwkv_fit + opt_res$par[3]))
    
    print('plotting')
    dodge = 1
    
    yr_diff = fmax(zf_cnts$yr) - fmin(zf_cnts$yr)
    
    plot_input = zf_cnts |> 
      mtt(d = as.IDate(paste0("2026-01-01")) + 7*wk - 3.5 + 
            dodge * (yr - 2023) - (yr_diff*dodge/2)) |> 
      mtt(yr = factor(yr)) 
    
    tstring = paste0('*', input$genus, "* - ", 
                     paste0(unlist(c(input$flw_phenos, input$lf_phenos)), 
                            collapse = paste0(' ', input$log_fun, ' ')))
    
    p1 = plot_input |> 
      ggplot(aes(d, prop)) + 
      geom_point(aes(color = yr, group = yr)) + 
      geom_line(aes(y = fitted, color = yr)) +
      scale_color_manual(values = pals::parula(8)[c(1,3,5,7)]) + 
      scale_x_date(labels = scales::label_date("%b"),
                   breaks = as.Date(paste0("2026-", 1:12, "-15"))) + 
      ylim(c(0,1)) + 
      theme_bw() + 
      theme(panel.grid.minor.x = element_blank(),
            axis.title.y = element_text(margin = margin(0,0,-10,0, 'pt')),
            plot.title = element_markdown()) + 
      labs(x = NULL,
           y = 'proportion',
           color = NULL,
           title = tstring)
    
    print('donesk...')
  
    sel_obs = d |>
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
      flw_mat = matrix(TRUE, nrow = fnrow(sel_obs))
    }
    
    if (!is.null(input$lf_phenos)) {
      lf_cols = uniq_lf |> 
        sbt(lf_pheno %iin% input$lf_phenos) |> 
        get_elem("cln_lf")
      
      lf_mat = sel_obs |> 
        slt(lf_cols) |> 
        qM()
      
    } else {
      lf_mat = matrix(TRUE, nrow = fnrow(sel_obs))
    }
    
    cnd_vec = cbind(flw_mat, lf_mat) |> 
      logical_fun() 
    
    sel_obs = sel_obs |> 
      mtt(meets_cnd = as.logical(cnd_vec)) 
    
    # Start to differ from above starting here. 
    
    obs_df = sel_obs |> 
      mtt(gs_fct = to_fct(species_2, input$genus), 
          d_hide = as.Date("2026-01-01") + 7*wk - 3.5) |> 
      roworder(meets_cnd)
    
    tstring = paste0('*', input$genus, "* - ", 
                     paste0(unlist(c(input$flw_phenos, input$lf_phenos)), 
                            collapse = paste0(' ', input$log_fun, ' ')))
    
    p2 = ggplot(obs_df, aes(d_hide, gs_fct)) + 
      geom_point(pch = 15, 
                 aes(color = meets_cnd),
                 position = position_jitter(width=.5, height = .2),
                 size = .6) + 
      scale_color_manual(values = c("grey", "black")) + 
      labs(color = "Specified\ncondition",
           x = NULL,
           y = NULL,
           title = tstring) + 
      scale_x_date(labels = scales::label_date("%b"),
                   breaks = as.Date(paste0("2026-", 1:12, "-15"))) + 
      theme_bw() + 
      theme(panel.grid.minor.x = element_blank(),
            plot.margin = unit(c(.3,.3,.3,1), "cm"),
            axis.text.y = element_markdown(),
            plot.title = element_markdown()) 
    
    g1 = ggplotGrob(p1)
    g2 = ggplotGrob(p2)
    
    grid::grid.newpage()
    grid::grid.draw(rbind(g1, g2))
  })
}

# Run the application 
shinyApp(ui = ui, server = server)
