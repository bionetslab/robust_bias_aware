require(igraph)
require(networkD3)

g<-read_graph(
  "covid.graphml",
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

