# generate patient report
rmarkdown::render(input = 'OMPARE.Rmd', 
                  params = list(topDir = 'data/PNOC008/',
                                fusion_method = 'arriba',
                                set_title = 'PNOC008 Report'))

