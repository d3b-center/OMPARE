view_dt <- function(dat, pageLength = 10, escape = F, rownames = F){
  DT::datatable(dat,
                escape = escape,
                rownames = rownames,
                extensions = c('Buttons', 'FixedColumns'), 
                options = list(scrollX = TRUE,
                               pageLength = pageLength,
                               dom = '<"top" Bpift>',
                               buttons = 
                                 list(list(
                                   extend = 'collection',
                                   buttons = c('csv'),
                                   text = 'Download'
                                 )))
  )
  
}
