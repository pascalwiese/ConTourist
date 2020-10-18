class Mesh:
    def __init__(self, nodes, elements):
        self.nodes = nodes
        self.elements = elements
        self.nodes_conn_element_ids = None
        self.element_neigh_ids = None
