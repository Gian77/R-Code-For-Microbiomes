# A color palette example and how to visualize it

palette_CB9 <-c("#88CCEE","#44AA99","#332288","#117733","#999933","#999999","#CC6677","#882255","#AA4499")

par(mar = rep(0, 4))
pie(rep(1, length(palette_CB9)), labels = sprintf("%d (%s)", 
          seq_along(palette_CB9),palette_CB9), col = palette_CB9)



