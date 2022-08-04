import networkx as nx
import matplotlib.pyplot as plt

def draw_module(filepath):
    module = nx.read_graphml(filepath)
    node_colors = []
    for _, data in module.nodes(data=True):
        if data['isSeed']:
            node_colors.append('red')
        else:
            node_colors.append('orange')
    fig1, ax1 = plt.subplots()
    nx.draw_networkx(module,node_color=node_colors)

draw_module('../output_miniseeds_1tree.graphml')
draw_module('../output_miniseeds_2tree.graphml')
draw_module('../output_miniseeds_3tree.graphml')
draw_module('../output_miniseeds_5tree.graphml')
draw_module('../output_miniseeds_10tree.graphml')
draw_module('../output_miniseeds_20tree.graphml')
draw_module('../output_miniseeds_30tree.graphml')


# draw_module('30_trees.graphml')