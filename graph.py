from graphviz import Digraph
import textwrap
import math
from collections import defaultdict, deque



# Colour palette
BRANCH_COLORS = [
    "#1565C0",   # deep blue
    "#2E7D32",   # forest green
    "#AD1457",   # crimson
    "#E65100",   # burnt orange
    "#4527A0",   # deep purple
    "#00695C",   # teal
    "#F9A825",   # amber
    "#37474F",   # slate
]

SHARED_COLOR       = "#546E7A"   
START_COLOR        = "#FFF176"   
TARGET_COLOR       = "#EF9A9A"   
SHARED_NODE_COLOR  = "#B3E5FC"  
DEFAULT_NODE_COLOR = "#F5F5F5"   

# Single-target  (Z-snake)
def create_pathway_pdf(data, nodes_per_row=3, filename="output"):
    dot = Digraph(engine='neato')
    dot.attr('node', shape='rect', style='rounded,filled',
             fillcolor=DEFAULT_NODE_COLOR, color='#000000',
             fontname='Helvetica', fontsize='9')
    dot.attr('edge', fontname='Helvetica', fontsize='9', color=SHARED_COLOR)

    compounds = []
    for start, _enzyme, end in data:
        if start not in compounds: compounds.append(start)
        if end   not in compounds: compounds.append(end)

    h_space, v_space = 3.5, 1.8
    for i, compound in enumerate(compounds):
        row = i // nodes_per_row
        col = i %  nodes_per_row
        x = col * h_space if row % 2 == 0 else (nodes_per_row - 1 - col) * h_space
        y = -(row * v_space)
        dot.node(compound, label=compound, pos=f"{x},{y}!")

    for start, enzyme, end in data:
        wrapped = "\n".join(textwrap.wrap(enzyme, width=15))
        dot.edge(start, end, label=wrapped)

    dot.render(filename, format='pdf', cleanup=True)
    print(f"PDF saved: {filename}.pdf")


# Sugiyama-style DAG layout
def _build_dag(pathway_lists):
    nodes      = set()
    edges      = set()
    edge_enzyme = {}
    edge_paths  = defaultdict(set)
    node_paths  = defaultdict(set)

    for path_idx, path in enumerate(pathway_lists):
        if not path:
            continue
        for s, enzyme, e in path:
            nodes.add(s)
            nodes.add(e)
            edges.add((s, e))
            edge_enzyme.setdefault((s, e), enzyme)
            edge_paths[(s, e)].add(path_idx)
            node_paths[s].add(path_idx)
            node_paths[e].add(path_idx)

    return nodes, edges, edge_enzyme, edge_paths, node_paths


def _assign_ranks(nodes, edges):
    in_degree  = {n: 0 for n in nodes}
    successors = defaultdict(set)
    for s, e in edges:
        successors[s].add(e)
        in_degree[e] += 1

    rank   = {n: 0 for n in nodes}
    queue  = deque(n for n in nodes if in_degree[n] == 0)

    while queue:
        node = queue.popleft()
        for child in successors[node]:
            if rank[node] + 1 > rank[child]:
                rank[child] = rank[node] + 1
            in_degree[child] -= 1
            if in_degree[child] == 0:
                queue.append(child)

    return rank


def _order_within_layers(layers, edges):
    predecessors = defaultdict(list)
    successors   = defaultdict(list)
    for s, e in edges:
        successors[s].append(e)
        predecessors[e].append(s)

    rank_of = {}
    for r, ns in layers.items():
        for n in ns:
            rank_of[n] = r

    max_rank = max(layers.keys())

    def barycenter(node, neighbours):
        if not neighbours:
            return float('inf')
        layer = layers[rank_of[neighbours[0]]]
        positions = {n: i for i, n in enumerate(layer)}
        return sum(positions.get(nb, 0) for nb in neighbours) / len(neighbours)

    # Top-down pass
    for r in range(1, max_rank + 1):
        layer = layers[r]
        layer.sort(key=lambda n: barycenter(n, predecessors[n]))

    # Bottom-up pass
    for r in range(max_rank - 1, -1, -1):
        layer = layers[r]
        layer.sort(key=lambda n: barycenter(n, successors[n]))

    return layers


def _compute_coordinates(layers, v_step=2.0, h_step=3.2):
    coord = {}
    for rank, layer_nodes in layers.items():
        n = len(layer_nodes)
        y = -rank * v_step
        # Centre the layer on x=0
        xs = [(i - (n - 1) / 2.0) * h_step for i in range(n)]
        for node, x in zip(layer_nodes, xs):
            coord[node] = (x, y)
    return coord


def _build_layout(pathway_lists, v_step=2.0, h_step=3.2):
    nodes, edges, edge_enzyme, edge_paths, node_paths = _build_dag(pathway_lists)

    if not nodes:
        return {}, {}, {}, {}

    rank   = _assign_ranks(nodes, edges)

    # Group nodes into layers
    layers = defaultdict(list)
    for node in nodes:
        layers[rank[node]].append(node)

    layers = _order_within_layers(dict(layers), edges)
    coord  = _compute_coordinates(layers, v_step=v_step, h_step=h_step)

    return coord, dict(edge_paths), edge_enzyme, dict(node_paths)



# Multi-target renderer
def create_multi_target_pathway(
        pathway_lists,
        target_names=None,
        start_name=None,
        nodes_per_row=4,            # kept for API compatibility, unused
        filename="multi_target_output",
        v_step=2.0,
        h_step=3.2):
    # ── 0. Guard / defaults 
    if target_names is None:
        target_names = [f"Target {i+1}" for i in range(len(pathway_lists))]

    valid_paths = [p for p in pathway_lists if p]
    if not valid_paths:
        print("No valid paths to draw.")
        return

    if start_name is None:
        start_name = valid_paths[0][0][0]

    # ── 1. Compute layout 
    coord, edge_paths, edge_enzyme, node_paths = _build_layout(
        pathway_lists, v_step=v_step, h_step=h_step)

    if not coord:
        print("Layout produced no nodes — check pathway data.")
        return

    # ── 2. Identify special nodes 
    target_node_names = {p[-1][2] for p in pathway_lists if p}
    shared_nodes      = {n for n, idxs in node_paths.items() if len(idxs) > 1}

    # ── 3. Build graphviz graph 
    dot = Digraph(engine='neato')
    dot.attr(overlap='false', splines='line')   # straight arrows
    dot.attr('node', shape='rect', style='rounded,filled',
             fillcolor=DEFAULT_NODE_COLOR, color='#424242',
             fontname='Helvetica', fontsize='9')
    dot.attr('edge', fontname='Helvetica', fontsize='8', arrowsize='0.7')

    # ── 4. Draw nodes 
    for node, (x, y) in coord.items():
        if node == start_name:
            fill, pw = START_COLOR, "2"
        elif node in target_node_names:
            fill, pw = TARGET_COLOR, "2"
        elif node in shared_nodes:
            fill, pw = SHARED_NODE_COLOR, "1"
        else:
            fill, pw = DEFAULT_NODE_COLOR, "1"

        dot.node(node,
                 label=node,
                 pos=f"{x:.3f},{y:.3f}!",
                 fillcolor=fill,
                 penwidth=pw)

    # ── 5. Draw edges 
    for (s, e), path_set in edge_paths.items():
        enzyme  = edge_enzyme.get((s, e), "")
        wrapped = "\n".join(textwrap.wrap(enzyme, width=14))

        if len(path_set) > 1:
            # Shared edge — draw once, thicker, neutral colour
            dot.edge(s, e, label=wrapped,
                     color=SHARED_COLOR, fontcolor=SHARED_COLOR,
                     penwidth="2.5")
        else:
            path_idx = next(iter(path_set))
            color    = BRANCH_COLORS[path_idx % len(BRANCH_COLORS)]
            dot.edge(s, e, label=wrapped,
                     color=color, fontcolor=color,
                     penwidth="1.8")

    # ── 6. Legend (placed below the lowest node) 
    min_y    = min(y for _, y in coord.values())
    legend_y = min_y - v_step * 2.2
    leg_y2   = legend_y - v_step * 1.1

    legend_items = [
        ("__leg_start__",  start_name or "Start",  START_COLOR,       "2", 0.0),
        ("__leg_target__", "Target compound",       TARGET_COLOR,      "2", 3.5),
        ("__leg_shared__", "Shared intermediate",   SHARED_NODE_COLOR, "1", 7.0),
    ]
    for nid, label, fill, pw, lx in legend_items:
        dot.node(nid, label=label, shape="rect", style="rounded,filled",
                 fillcolor=fill, fontname="Helvetica", fontsize="8",
                 pos=f"{lx:.1f},{legend_y:.2f}!", penwidth=pw)

    for i, (tname, color) in enumerate(zip(target_names, BRANCH_COLORS)):
        dot.node(f"__leg_b{i}__", label=f"→ {tname}",
                 shape="rect", style="rounded,filled",
                 fillcolor=color + "33", color=color, fontcolor=color,
                 fontname="Helvetica", fontsize="8",
                 pos=f"{i * 3.5:.1f},{leg_y2:.2f}!", penwidth="2")

    dot.render(filename, format='pdf', cleanup=True)
    print(f"Multi-target PDF saved: {filename}.pdf")

# Smoke-test (glucose→alanine + glucose→ethanol) 
if __name__ == "__main__":

    path_alanine = [
        ['D-Glucose', 'hexokinase', 'D-Glucose 6-phosphate'],
        ['D-Glucose 6-phosphate', 'glucose-6-phosphate isomerase', 'D-Fructose 6-phosphate'],
        ['D-Fructose 6-phosphate', 'transaldolase', 'D-Glyceraldehyde 3-phosphate'],
        ['D-Glyceraldehyde 3-phosphate', '1-deoxy-D-xylulose-5-phosphate synthase', '1-Deoxy-D-xylulose 5-phosphate'],
        ['1-Deoxy-D-xylulose 5-phosphate', '1-deoxy-D-xylulose-5-phosphate synthase', 'Pyruvate'],
        ['Pyruvate', 'alanine transaminase', 'L-Alanine'],
    ]

    path_ethanol = [
        ['D-Glucose', 'hexokinase', 'D-Glucose 6-phosphate'],
        ['D-Glucose 6-phosphate', 'glucose-6-phosphate isomerase', 'D-Fructose 6-phosphate'],
        ['D-Fructose 6-phosphate', 'transaldolase', 'D-Glyceraldehyde 3-phosphate'],
        ['D-Glyceraldehyde 3-phosphate', '2-dehydro-3-deoxy-phosphogluconate aldolase', '2-Dehydro-3-deoxy-6-phospho-D-gluconate'],
        ['2-Dehydro-3-deoxy-6-phospho-D-gluconate', '2-dehydro-3-deoxy-phosphogluconate aldolase', 'Pyruvate'],
        ['Pyruvate', 'formate C-acetyltransferase', 'Acetyl-CoA'],
        ['Acetyl-CoA', 'acetaldehyde dehydrogenase (acetylating)', 'Acetaldehyde'],
        ['Acetaldehyde', 'alcohol dehydrogenase (NADP+)', 'Ethanol'],
    ]

    create_multi_target_pathway(
        [path_alanine, path_ethanol],
        target_names=["L-Alanine", "Ethanol"],
        start_name="D-Glucose",
        filename="test_glucose_alanine_ethanol",
        v_step=2.0,
        h_step=3.5,
    )
