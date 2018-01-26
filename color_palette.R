# Here some color blind palette I have used for my data.
# Also, a quick a code for a quick plot to visualize them.

palette_CB3 <-c("#44AA99","#332288","#999933")
palette_CB2 <-c("#CC6677", "#999999","#AA4499")
palette_CB4 = c("#88CCEE","#DDCC77","#117733", "#882255")
palette_CB6 <-c("#332288","#88CCEE","#117733","#DDCC77","#FF899d","#AA4499") #"#CC6677"
palette_CB8 <-c("#999933","#DDCC77","#CC6677","#AA4499","#332288","#88CCEE","#44AA99","#117733")
palette_CB9 <-c("#88CCEE","#44AA99","#332288","#117733","#999933","#999999","#CC6677","#882255","#AA4499")

par(mar = rep(0, 4))
pie(rep(1, length(palette_CB9)), labels = sprintf("%d (%s)", 
          seq_along(palette_CB9),palette_CB9), col = palette_CB9)



