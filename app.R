#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#

library(shiny)
library(shinyjs) #https://stackoverflow.com/questions/24138108/capturing-cat-output-periodically-for-r-shiny-output-renderprint
## app.R ##
library(shinydashboard)
library(readxl)
library(nlme)
library(tidyverse)
library(lemon)
library(gridExtra)

options(shiny.maxRequestSize = 50 * 1024^2) #파일 업로드 최대 용량 50MB

ui <- dashboardPage(
  dashboardHeader(title = "Linear Mixed Effect Modeling for RNA sequencing longitudinal data"),
  ## Sidebar content
  dashboardSidebar(
    sidebarMenu(
      menuItem(text = "Importing file", tabName = "dataset",icon = icon("file")#text = "Importing file", tabName = "dashboard",icon = icon("file")#,
               #menuSubItem(text = "Importing file", tabName = "dataset")
      ), #https://stackoverflow.com/questions/51523783/trying-to-put-an-excel-csv-file-into-r-shinydashboard/51594257#51594257
      menuItem("Dashboard", tabName = "dashboard", icon = icon("dashboard")),
      menuItem("Widgets", tabName = "widgets", icon = icon("th")),
      menuItem("Save pdfs", tabName = "figures", icon = icon("download"))
      #icons: https://fontawesome.com/icons
    )
  ),
  ## Body content
  dashboardBody(
    tabItems(
      tabItem(tabName = "dataset",
              fluidRow(
                box(width = 12,
                    fileInput(inputId = "file",
                              label = "Choose a file",
                              #accept = c(".xlsx",".csv")
                              accept = c(".csv")
                    ),
                    tableOutput(outputId = "Contents"),
                    verbatimTextOutput(outputId = "Data")
                )
              )
      ),
      # First tab content
      tabItem(tabName = "dashboard",
              fluidRow(
                box(plotOutput("plot1", height = 250)),
                
                box(
                  title = "Controls",
                  sliderInput("slider", "Number of observations:", 1, 100, 50)
                )
              )
      ),
      
      # Second tab content
      tabItem(tabName = "widgets",
              h2("Widgets tab content")
      ),
      
      # (SC) Third tab content
      tabItem(tabName = "figures",
              downloadLink('downloadData', 'Download')
      )
    )
  )
)

server <- function(input, output) {
  ##--------------------------------------------------------------------------
  ## modeling 
  ##--------------------------------------------------------------------------
  ## - mainTCGE        : main function
  ##  - data_wide2long : convert wide data to long
  ##  - report_missing : check if there are missing values
  ##  - fit_LMM        : run linear mixed models
  ##--------------------------------------------------------------------------
  fit_LMM <- function(data, 
                      degree = 2, 
                      time_log = F,
                      time_center = F,
                      ctrl = lmeControl()){
    # data = data
    # degree = degree
    # time_log = time_log
    # time_center = time_center
    # ctrl = ctrl
    
    
    ctrl <- ctrl
    
    ## data processing
    data <- 
      data %>% 
      mutate(
        time = as.numeric(time),
        time = ifelse(rep(time_log, length(time)), log1p(time), time),
        time = ifelse(rep(time_center, length(time)), scale(time, T, F), time)
      )
    
    ## regression formulas
    cn_groups <- setdiff(colnames(data), c("name", "id", "time", "value"))
    if(length(cn_groups) == 0){
      uni_group <- TRUE
      fm <<- as.formula(paste("value ~ poly(time, ", degree, ", raw = TRUE)"))
    } else{
      uni_group <- FALSE
      fm_x <- paste("poly(time, ", degree, ", raw = TRUE) *", cn_groups,
                    collapse = " + ")
      
      fm <<- as.formula(paste("value ~ ", fm_x))
    }
    fm_random <<- as.formula("~ 1 | id")
    cor_struct <<- corCAR1(form = ~ 0 + time | id, fixed = F)
    
    cnt <- 0
    lmm_error_flag <- 1
    while(lmm_error_flag != 0){
      cnt <- cnt + 1
      
      res <- try({
        fit <- lme(fixed = fm,
                   random = fm_random, 
                   correlation = cor_struct,
                   data = data,
                   control = ctrl,
                   method = "ML")
        fit_null <- update(fit, value ~ 1, 
                           control = ctrl)
        pvalue <- anova(fit, fit_null)$`p-value`[2]
        ## return
        list(
          fit = fit,
          # fit_null = fit_null
          pvalue = pvalue
        )
      }, silent = T)
      if (inherits(res, "try-error")) {
        if(grepl("convergence error code = 1", res)){
          opt_old <- ctrl$opt
          opt_new <- setdiff(c("nlminb", "optim"), opt_old)
          ctrl$opt <- opt_new
        }
      } else{
        ## good to go
        fit <- res$fit
        pvalue <- res$pvalue
        lmm_error_flag <- 0
      }
      if(cnt == 3){
        ## convergence error
        fit <- NA
        pvalue <- NA
        lmm_error_flag <- res
        break
      }
    }
    
    return(
      list(fit = fit,
           pvalue_overall = pvalue,
           data = data,
           degree = degree,
           uni_group = uni_group,
           error_flag = lmm_error_flag)
    )
  }
  
  
  data_wide2long <- function(data_exp, data_design){
    data_long <- data_exp %>% 
      mutate(name = rownames(data_exp)) %>% 
      pivot_longer(cols = -name, 
                   names_to = c("id", "time"), 
                   names_sep = "_",
                   values_to = "value")
    
    data_design <- data_design %>% 
      mutate(id = as.character(id),
             time = as.character(time))
    
    data_long <- 
      inner_join(data_long, data_design, by = c("id", "time"))
    
    ## delete missing observations
    data_long <- na.omit(data_long)
    
    return(data_long)
  }
  report_missing <- function(data, N_obs){
    n_miss <- N_obs - nrow(data)
    name <- unique(data$name)
    if(n_miss == 0){
      ## no missing
    } else{
      cat(sprintf("%s: there are missing observations!\n", 
                  name))
    }
    return(data.frame(name = name, 
                      n_miss = n_miss))
  }
  
  mainTCGE <- function(data_exp, data_design,
                       timepoint,
                       degree = 2,
                       alpha = 0.05, 
                       mt_type = "BH", 
                       time_log = T,
                       time_center = T,
                       ctrl = NULL,
                       #save_pdf = NULL){
                       save_pdf = "./"){
    
    # data_exp = data_exp[1:10,]
    # data_design = data_design
    # timepoint = tp
    # degree = 2
    # alpha = 0.05
    # mt_type = "BH"
    # time_log = T
    # time_center = T
    # ctrl = lmeControl(opt = "optim",
    #                   msMaxIter = 2000)
    # save_pdf = "./result/250115_GSE3946_1"
    
    
    
    if(is.null(ctrl)){
      ctrl <- lmeControl()
    }
    
    withProgress(message='(2) Fit LMM', value=0, {
      ## (1) Data preprocessing
      shinyjs::html("(1) Start processing data...")
      cat("(1) Start processing data...")
      data_exp <- as.data.frame(data_exp)
      data_design <- as.data.frame(data_design)
      N_obs <- ncol(data_exp)
      
      data_long <- data_wide2long(data_exp, data_design)
      rm("data_exp", "data_design")
      
      res <- list()
      res$data_info <- lapply(data_long %>% 
                                select(-value), 
                              function(x) str_sort(unique(x), numeric = T))
      bool_group <- length(setdiff(names(res$data_info), c("name", "id", "time"))) > 0
      
      raw_time <- names(timepoint)
      user_time <- as.numeric(timepoint)
      
      data_long <- 
        data_long %>% 
        mutate(time = user_time[as.numeric(factor(time, 
                                                  levels = raw_time,
                                                  labels = user_time))])
      cat(" Done!\n")
      
      
      ## (2) Testing significance of the overall model
      cat("(2) Testing significance of the overall model...\n")
      list_ovr <- list_error <- list()
      pvalue_ovr <- c()
      res_summ <- list()
      cnt <- cnt_err <- 0
      #(SC) 프로그레스 바 (변경 필요)
      #https://shiny.posit.co/r/articles/build/progress/
      #pb <- txtProgressBar(min = 0, max = length(res$data_info$name), style = 3)
      
      for(i in seq_along(res$data_info$name)){
        # i=1
        #setTxtProgressBar(pb, i)
        
        # (SC) progress bar 교체 Increment the progress bar, and update the detail text.
        incProgress(1/length(seq_along(res$data_info$name)), detail = paste("Doing part", i))
        
        data <- 
          data_long %>% 
          filter(name == res$data_info$name[i]) 
        
        res_summ[[i]] <- report_missing(data, N_obs)
        res_summ[[i]]$error <- 0
        
        # data_long <- data_arrange(data)
        f_lmm <- fit_LMM(data = data,
                         degree = degree,
                         time_log = time_log,
                         time_center = time_center,
                         ctrl = ctrl)
        
        pvalue_ovr[i] <- f_lmm$pvalue_overall
        
        if(f_lmm$error_flag != 0){
          cnt_err <- cnt_err + 1
          
          res_summ[[i]]$error <- 1
          
          list_error[[cnt_err]] <- 
            list(name = res$data_info$name[i],
                 pvalue_ovr = pvalue_ovr[i],
                 data = data,
                 error = f_lmm$error_flag)
          
          next
        }
        
        ## to minimize computational redundancy
        if(pvalue_ovr[i] < alpha){
          cnt <- cnt + 1
          # name_genes_ovr[cnt] <- name_genes[i]
          list_ovr[[cnt]] <- 
            list(name = res$data_info$name[i],
                 pvalue_ovr = pvalue_ovr[i],
                 data = data,
                 f_lmm = f_lmm)
        }
      }
      
      cat("\n Summarzing results...")
      res$error <- list_error
      res_summ <- do.call(rbind, res_summ)
      names(list_ovr) <- unlist(lapply(list_ovr, function(x) x$name))
      pvalue_adj <- p.adjust(pvalue_ovr, method = mt_type)
      idx_sig <- which(pvalue_adj < alpha)
      
      res_summ$p_ovr <- pvalue_ovr
      res_summ$p_adj <- pvalue_adj
      res_summ$signf <- 1 * (pvalue_adj < alpha)
      res$df_summ <- res_summ
      cat(" Done!\n")
      
      ## (3) Testing significance of individual effects
      cat("(3) Testing significance of individual effects...")
      name_sig <- res$data_info$name[idx_sig]
      
      
      if(!bool_group){
        ## Done if group does not exist
        pvalue_effect <- data.frame(name = name_sig,
                                    overall = pvalue_ovr[idx_sig])
      } else{
        ## results after multiple testing if group exists
        pvalue_effect <- list()
        for(ii in seq_along(name_sig)){
          name_gene <- name_sig[[ii]]
          
          ## test effect (do not extract this job into a separate user-defined function; it causes an error because of "update.lmm")
          lmm <- list_ovr[[name_gene]]$f_lmm
          degree <- lmm$degree
          cn_groups <- setdiff(colnames(lmm$data), c("name", "id", "time", "value"))
          fm_x_t_g <- paste("poly(time, ", degree, ", raw = TRUE) +", cn_groups,
                            collapse = " + ")
          fm_x_t_tg <- paste("poly(time, ", degree, ", raw = TRUE)+", 
                             paste("poly(time, ", degree, ", raw = TRUE):", cn_groups,
                                   collapse = " + "))
          fm_x_g <- paste(cn_groups, collapse = " + ")
          
          fm_t_g <- as.formula(paste("value ~ ", fm_x_t_g))
          fm_t_tg <- as.formula(paste("value ~ ", fm_x_t_tg))
          fm_g <- as.formula(paste("value ~ ", fm_x_g))
          
          
          list_anova <- 
            list(
              effect_tg = anova(lmm$fit, update(lmm$fit, fm_t_g)),
              effect_g = anova(lmm$fit, update(lmm$fit, fm_t_tg)),
              effect_t = anova(lmm$fit, update(lmm$fit, fm_g))
            )
          pvalue <- unlist(lapply(list_anova, function(x) x$`p-value`[2]))
          names(pvalue) <- c("interaction", "group", "time")
          # }
          pvalue_effect[[ii]] <- pvalue
          # pvalue_effect[[ii]] <- test_effect(list_ovr[[name_gene]]$f_lmm)
        }
        pvalue_effect <- 
          do.call(rbind, pvalue_effect) %>% 
          as.data.frame() %>% 
          mutate(name = name_sig,
                 overall = pvalue_ovr[idx_sig]) %>% 
          relocate(name)
        
        # rownames(pvalue_effect) <- name_sig
      }
      res$df_signf <- pvalue_effect
      cat(" Done!\n")
      
      ## (4) Extracting (and plotting if specified) results
      cat("(4) Extracting (and plotting if specified) results...")
      tmp <- pvalue_effect %>% 
        select(-c(name))
      res$gene_signf <- sapply(tmp, function(x){
        pval_adj <- p.adjust(x, method = mt_type)
        id_gene <- which(pval_adj < alpha)
      })
      if(!is.null(save_pdf)){
        # 이 부분 맞게 수정해야
        dir.create(save_pdf)
        for(jj in 1:ncol(tmp)){
          pval_adj <- p.adjust(tmp[,jj], method = mt_type)
          id_gene <- which(pval_adj < alpha)
          name_sig <- pvalue_effect$name[id_gene]
          
          pdf(sprintf("%s/time_course_%s.pdf", 
                      save_pdf,
                      colnames(tmp)[jj]),
              onefile = T,
              width = 12,
              height = 8)
          for(kk in seq_along(name_sig)){
            if(bool_group){
              ## If group does exist
              data_tmp <- list_ovr[[name_sig[kk]]]$f_lmm$data %>% 
                mutate(inter_group = interaction(across(-c(name, id, time, value))))
            } else{
              data_tmp <- list_ovr[[name_sig[kk]]]$f_lmm$data %>% 
                mutate(inter_group = "No group")
            }
            gg_fig  <- list()
            gg_fig[[1]] <- 
              ggplot(data_tmp, 
                     mapping = aes(x = time, 
                                   y = value, 
                                   color = inter_group)) + 
              geom_point(alpha = 0.5, shape = 4) + 
              geom_line(aes(group = id), alpha = 0.2) + 
              
              labs(title = sprintf("%s", name_sig[kk]),
                   color = "Group") + 
              theme_bw() + 
              theme(plot.title = element_text(size = 30),
                    text = element_text(size = 25), 
                    axis.text = element_text(size = 12),
                    axis.title = element_text(size = 20), 
                    strip.text = element_text(size = 12),
                    legend.position = "top"
              )
            if(!bool_group){
              gg_fig[[1]] <- gg_fig[[1]] + 
                guides(color = "none")
            }
            
            
            
            data_tmp <- list_ovr[[name_sig[kk]]]$f_lmm$fit$data %>% 
              select(-c(name, id, value)) %>% 
              unique()
            vn_grid <- c(list(time = seq(min(data_tmp$time),
                                         max(data_tmp$time),
                                         length.out = 100)),
                         data_tmp %>% 
                           select(-time) %>% 
                           lapply(., unique))
            
            df_grid <- expand.grid(vn_grid)
            df_grid <- data.frame(df_grid,
                                  value = predict(list_ovr[[name_sig[kk]]]$f_lmm$fit, 
                                                  newdata = df_grid, 
                                                  level = 0))
            if(bool_group){
              ## If group does exist
              df_grid <- df_grid %>%
                mutate(inter_group = interaction(across(-c(time, value))))
            } else{
              df_grid <- df_grid %>% 
                mutate(inter_group = "No group")
            }
            gg_fig[[2]] <- gg_fig[[1]] + 
              stat_summary(fun = mean, 
                           geom = "line",
                           linewidth = 1,
                           linetype = "dotted",
                           alpha = 1)
            gg_fig[[3]] <- gg_fig[[1]] + 
              geom_line(data = df_grid, 
                        linetype = "solid", 
                        linewidth = 1,
                        alpha = 1)
            if(bool_group){
              invisible(
                grid_arrange_shared_legend(gg_fig[[2]],
                                           gg_fig[[3]],
                                           nrow = 1,
                                           position = "top"))
            } else{
              invisible(
                grid.arrange(gg_fig[[2]],
                             gg_fig[[3]],
                             nrow = 1)
              )
            }
          } # end of for(kk)
          dev.off()
        } # end of for(jj)
      } # end of if
      rm("tmp", "data_tmp")
      cat(" Done!\n")
    })
    
    return(res)
  }
  
  set.seed(122)
  histdata <- rnorm(500)
  
  output$Data <- renderPrint({
    if(is.null(input$file)){
      print("Import Excel data file")
    } else {
      inFile <- input$file
      df <- read.csv(inFile$datapath, header=T)
      #print(df)
      
      data_exp <- df %>%
        rename(id_gene = 1) %>% 
        as.data.frame
      rownames(data_exp) <- data_exp$id_gene
      data_exp$id_gene <- NULL
      
      data_design <- do.call(rbind, stringr::str_split(colnames(data_exp), "_"))
      colnames(data_design) <- c("MLA1", "Bgh", "time", "rep")
      
      data_design <- as.data.frame(data_design) %>%
        group_by_at(vars(-time)) %>% 
        mutate(id = cur_group_id()) %>% 
        ungroup() %>% 
        select(-rep) %>% 
        relocate(id, time)
      
      print(data_design)
      
      coln <- data_design %>% 
        # unite("coln", everything()) %>% 
        unite("coln", id, time) %>%  
        pull(coln)
      
      colnames(data_exp) <- coln
      tp <- c(6, 12, 18, 24)
      names(tp) <- unique(data_design$time)
      
      ##--------------------------------------------------------------------------
      ## Run mainTCGE
      res_TCGE <- mainTCGE(
        # data_exp = data_exp[1:5,,drop = F],
        data_exp = data_exp[1:1000,,drop = F],
        data_design = data_design,
        timepoint = tp,
        time_log = T,
        time_center = T,
        ctrl = lmeControl(opt = "optim",
                          msMaxIter = 2000), 
        #save_pdf = "./result/250115_GSE3946_1")
        save_pdf = NULL)
      
      ## Error report
      res_TCGE$df_signf
      res_TCGE$df_summ
      
      res_TCGE$error %>% length
      msg_err <- factor(unlist(map(res_TCGE$error, "error")))
      levels(msg_err)
    }
  })
  
  output$plot1 <- renderPlot({
    data <- histdata[seq_len(input$slider)]
    hist(data)
  })
  
  output$downloadData <- downloadHandler(
    filename = function() {
      paste('data-', Sys.Date(), '.csv', sep='')
    },
    content = function(con) {
      write.csv(data, con)
    }
  )
}

# Run the application 
shinyApp(ui = ui, server = server)
