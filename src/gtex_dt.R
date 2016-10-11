gtex.dt <- DT::datatable(gtex.df,rownames=FALSE,class="display",extensions="KeyTable",options=list(keys=TRUE)) %>% DT::formatRound(c("log10.bf","theta.mean","theta.95lb","theta.95ub"), 3)
