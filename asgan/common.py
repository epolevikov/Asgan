

def inv_sign(sign):
    return ["+", "-"][sign == "+"]


def build_block2edge_dict(adjacency_graph):
    block2edge = dict()

    for (node_from, node_to, data) in adjacency_graph.edges(data=True):
        edge_name = data["name"]

        if edge_name.startswith("+") or edge_name.startswith("-"):
            block2edge[edge_name] = (node_from, node_to)

    return block2edge


def build_id2block_dict(synteny_blocks):
    id2block = dict()

    for blocks in synteny_blocks.values():
        for block in blocks:
            id2block[block.signed_id()] = block

    return id2block
