import networkx as nx

class CoVexEdgeWeight:
    """
    An alternative edge weight that penalizes hubs. Note that the corresponding
    steiner tree no longer minimizes the number of vertices. The solutions could
    be more path like because in especially hubs are great steiner points that are now
    no longer attractive to use.
    """
    def __init__(self, graph: nx.Graph, lambda_):
        self.graph = graph
        self.avg_degree = self._calculate_avg_degree()
        self.lambda_ = lambda_

    def _calculate_avg_degree(self):
        return 2 * self.graph.number_of_edges() / self.graph.number_of_nodes()

    def __getitem__(self, e):
        return (1 - self.lambda_) * self.avg_degree + self.lambda_ * 0.5 * (
                self.graph.degree()[e[0]] + self.graph.degree()[e[1]])


class UnitEdgeWeight:
    """
    Simple unit edge weights. Every edge has the same weight. A Steiner Tree, thus,
    minimizes just the number of vertices.
    """
    def __getitem__(self, e):
        return 1.0




class BiasAwareEdgeWeight_Additive:
    """
    TODO
    """
    def __init__(self, graph: nx.Graph, lambda_):
        self.graph = graph
        self.lambda_=lambda_
        self.average_bias=self._calculate_average_bias()
        
    def _calculate_average_bias(self):
        biases = nx.get_node_attributes(self.graph, 'bias_data')
        sum_=0
        for a in self.graph.nodes:
            sum_ = sum_ + biases[a]
        return sum_ / self.graph.number_of_nodes()
    
    def __getitem__(self, e):
        average_edge_bias = (self.graph.nodes[e[0]]['bias_data']+self.graph.nodes[e[1]]['bias_data']) * 0.5
        return (1-self.lambda_) * self.average_bias + self.lambda_ * average_edge_bias
            
class BiasAwareEdgeWeight_AdditiveMax:
    """
    TODO
    """
    def __init__(self, graph: nx.Graph, lambda_):
        self.graph = graph
        self.lambda_=lambda_
        self.average_max_bias=self._calculate_average_max_bias()
        
    def _calculate_average_max_bias(self):
        biases = nx.get_node_attributes(self.graph, 'bias_data')
        sum_=0
        for u, v in self.graph.edges:
            sum_ = sum_ + max(biases[u], biases[v])
        return sum_ / self.graph.number_of_edges()
    
    def __getitem__(self, e):
        max_edge_bias = max(self.graph.nodes[e[0]]['bias_data'], self.graph.nodes[e[1]]['bias_data'])
        return (1-self.lambda_) * self.average_max_bias + self.lambda_ * max_edge_bias


class BiasAwareEdgeWeight_Exponential:
    """
    TODO
    """
    def __init__(self, graph: nx.Graph, lambda_):
        self.graph = graph
        self.lambda_=lambda_
        
    def __getitem__(self, e):
        sum_edge_bias = self.graph.nodes[e[0]]['bias_data']+self.graph.nodes[e[1]]['bias_data']
        return sum_edge_bias ** self.lambda_
            
  
class BiasAwareEdgeWeight_ExponentialMax:
    """
    TODO
    """
    def __init__(self, graph: nx.Graph, lambda_):
        self.graph = graph
        self.lambda_=lambda_
        
    def __getitem__(self, e):
        max_edge_bias = max(self.graph.nodes[e[0]]['bias_data'], self.graph.nodes[e[1]]['bias_data'])
        return max_edge_bias ** self.lambda_



















