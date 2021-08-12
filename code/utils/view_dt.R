library(DT)

# more functionality
view_dt <- function(dat, pageLength = 10, escape = F, rownames = F){
  DT::datatable(dat,
                escape = escape,
                rownames = rownames,
                filter = "top",
                extensions = c('Buttons', 'FixedColumns'), 
                options = list(scrollX = TRUE,
                               #scrollY="100vh",
                               scrollCollapse = FALSE,
                               pageLength = pageLength,
                               dom = '<"top" Bpift>',
                               lengthMenu = list(c(10, 20, 30, -1), c('10', '20', '30', 'All')),
                               searchHighlight = TRUE,
                               buttons = 
                                 list('pageLength', 
                                      list(
                                        extend = 'collection',
                                        buttons = c('csv'),
                                        text = 'Download'
                                      )))
  )
  
}

# minimal
view_dt_minimal <- function(dat){
  DT::datatable(dat,
                escape = F,
                rownames = F,
                extensions = c('Buttons', 'FixedColumns'), 
                class = 'cell-border stripe',
                options = list(pageLength = nrow(dat),
                               dom = '<"top" Bift>',
                               searchHighlight = TRUE,
                               buttons = c('copy', 'csv', 'excel')
                               ))
  
}
