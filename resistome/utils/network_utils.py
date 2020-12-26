import igraph
from collections import defaultdict

"""

Utility class for some network analysis functions: 

1. Conversion to igraph from networkx objects
2. Clustering of various sorts
3. Eigenvector centrality
4. Betweenness centrality

"""


def networkx_to_igraph(networkx_graph, directed=True):

    """
    
    Converts networkx graph to igraph equivalent.
    
    :param networkx_graph: networkx graph object
    :param directed: true/false, default is true
    :return: igraph object, dict of (int: node str), edge weights
    """

    nodes = networkx_graph.nodes()
    edges = networkx_graph.edges()

    nodes_to_int = {}
    int_to_nodes = {}

    for i in range(0, len(nodes)):
        nodes_to_int[nodes[i]] = i
        int_to_nodes[i] = nodes[i]

    # convert edges

    converted_edges = []

    converted_weights = []

    for edge in edges:
        node1 = nodes_to_int[edge[0]]
        node2 = nodes_to_int[edge[1]]

        if 'weight' in networkx_graph[edge[0]][edge[1]]:
            converted_weights.append(networkx_graph[edge[0]][edge[1]]['weight'])
        converted_edges.append((node1, node2))

    # build igraph,
    # compute betweenness,
    # convert btw scores back to networkx node names

    if not directed:
        directed = None

    if len(converted_weights) > 0:
        G = igraph.Graph(n=len(nodes), edges=converted_edges, directed=directed,
                         edge_attrs={'weight': converted_weights})
    else:
        G = igraph.Graph(n=len(nodes), edges=converted_edges, directed=directed)

    return G, int_to_nodes, converted_weights


def compute_igraph_clustering(networkx_graph, method=None):

    """
    
    Wrapper for igraph clustering methods.
    
    # clusters graphs into communities using leading eigenvector by MEJ Newmann method
    # citation: Finding community structure in networks using the eigenvectors of matrices
    # I tried other methods, but they just make one super cluster with >1000 others with very few nodes in them.
    # betweenness: too computationally expensive
    # multilevel: few clusters with many nodes
    # spinglass: doesn't work (bug or user error?)
    # infomap: few clusters with many nodes
    # MCL: haven't tried, need to build a wrapper for calling from python, build the software for windows
    
    :param networkx_graph:  networkx graph object
    :param method: clustering method (leading_eigenvector, infomap, edge_betweeness, multilevel, walktrap, fastgreedy
    label_propagation)-see igraph documentation for more details
    :return: dict (node str: cluster id), iterable of iterables for cluster composition
    """

    G, int_to_nodes, converted_weights = networkx_to_igraph(networkx_graph, directed=True)

    if len(converted_weights) == 0:
        converted_weights = None

    if method == 'leading_eigenvector':
        clusters = G.community_leading_eigenvector(weights=converted_weights)  # .as_clustering()
    elif method == 'infomap':
        clusters = G.community_infomap(edge_weights=converted_weights)
    elif method == 'edge_betweenness':
        clusters = G.community_edge_betweenness(weights=converted_weights).as_clustering()
    elif method == 'multilevel':
        clusters = G.as_undirected().community_multilevel(weights=converted_weights)
    elif method == 'walktrap':
        clusters = G.community_walktrap(weights=converted_weights).as_clustering()
    elif method == 'fastgreedy':
        clusters = G.as_undirected().community_fastgreedy(weights=converted_weights).as_clustering()
    elif method == 'label_propagation':
        clusters = G.as_undirected().community_label_propagation(weights=converted_weights)
    else:
        # break if invalid clustering selection selected
        raise AssertionError('Unrecognized clustering method provided: %s' % method)

    membership = clusters.membership

    nodes_to_clusters = {}

    cluster_membership = defaultdict(list)

    for i in range(0, len(membership)):
        node_id = int_to_nodes[i]
        nodes_to_clusters[node_id] = membership[i]

        cluster_membership[membership[i]].append(node_id)

    return nodes_to_clusters, cluster_membership


def compute_igraph_evcentrality(networkx_graph, directed=True):

    """
    Computes eigenvector centrality of networkx_graph.
    
    :param networkx_graph: networkx graph object
    :param directed: true/false
    :return: dict (node str: float (centrality)
    """

    G, int_to_nodes, converted_weights = networkx_to_igraph(networkx_graph, directed=directed)

    if len(converted_weights) == 0:
        converted_weights = None

    ec = G.eigenvector_centrality(weights=converted_weights)

    ec_dict = {}

    for i in range(0, len(ec)):
        ec_dict[int_to_nodes[i]] = ec[i]

    return ec_dict


def compute_igraph_betweenness(networkx_graph, directed=True):

    """
    Computes betweenness centrality of networkx_graph.

    :param networkx_graph: networkx graph object
    :param directed: true/false
    :return: dict (node str: float (centrality)
    """

    G, int_to_nodes, converted_weights = networkx_to_igraph(networkx_graph, directed=directed)

    if len(converted_weights) == 0:
        converted_weights = None

    bc = G.betweenness(weights=converted_weights)

    # normalizing value
    max_value = max(bc)

    bc_dict = {}
    # 0-1.0 scale for all networks
    for i in range(0, len(bc)):
        bc_dict[int_to_nodes[i]] = float(bc[i]) / float(max_value)

    return bc_dict