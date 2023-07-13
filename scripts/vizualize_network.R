require(igraph)
require(networkD3)
require(optparse)

option_list <- list(
make_option(c("-n", "--network"), type = 'character',
            help = "network file"))
opt <- parse_args(OptionParser(option_list=option_list))


g<-read_graph(
  opt$network,
  format ="graphml"
)


# Convert to object suitable for networkD3
G_d3 <- igraph_to_networkD3(g)
G_d3$nodes$group<-1
G_d3$nodes$gene<-V(g)$id
# Create force directed network plot
fn<-forceNetwork(G_d3$links, Nodes = G_d3$nodes,
             NodeID = 'gene', Group = 'group',opacityNoHover = 1)


saveNetwork(fn, file = 'covid.html')

