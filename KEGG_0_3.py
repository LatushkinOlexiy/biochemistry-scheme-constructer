import os
import re
import time
import pickle
import requests
import networkx as nx
import xml.etree.ElementTree as ET
from Bio.KEGG import REST
from rapidfuzz import process
import translate as tr

ORGANISM_CACHE = "organism_cache.pkl"


def load_kegg_organisms():
    if os.path.exists(ORGANISM_CACHE):
        with open(ORGANISM_CACHE, "rb") as f:
            return pickle.load(f)

    data = REST.kegg_list("organism").read()
    organisms = {}
    names = []

    for line in data.strip().split("\n"):
        parts = line.split("\t")
        code = parts[1]
        name = parts[2]
        organisms[name] = code
        names.append(name)

    with open(ORGANISM_CACHE, "wb") as f:
        pickle.dump((organisms, names), f)

    return organisms, names


def organism_name_to_code(user_input):
    organisms, names = load_kegg_organisms()
    match = process.extractOne(user_input, names)

    if match is None:
        raise ValueError("No organism match found")

    best_name = match[0]
    code = organisms[best_name]
    print(f"Matched organism: {best_name} → {code}")
    return code


CACHE_DIR = "kgml_cache"
GRAPH_CACHE = "compound_graph_with_reactions.pkl"
ENZYME_CACHE = "enzyme_cache.pkl"

os.makedirs(CACHE_DIR, exist_ok=True)


 
# Name → KEGG compound ID
def name_to_kegg_id(name):
    result = REST.kegg_find("compound", name).read()

    if not result:
        raise ValueError(f"No compound found for: {name}")

    candidates = []

    for line in result.strip().split("\n"):
        entry, names = line.split("\t")
        compound_id = entry.replace("cpd:", "")
        name_list = [n.strip().lower() for n in names.split(";")]

        if name.lower() in name_list:
            return compound_id

        candidates.append((compound_id, name_list))


 
# Compound ID → name
def compound_id_to_name(cid):
    data = REST.kegg_get(f"cpd:{cid}").read()
    for line in data.split("\n"):
        if line.startswith("NAME"):
            return line.replace("NAME", "").strip().split(";")[0]
    return cid


 
# Reaction → EC numbers
def get_reaction_enzymes(rid):
    rid = rid.split()[0].replace("rn:", "").strip()
    if not rid:
        return []
    try:
        data = REST.kegg_get(f"rn:{rid}").read()
    except Exception as e:
        print(f"  Warning: could not fetch reaction {rid}: {e}")
        return []
    for line in data.split("\n"):
        if line.startswith("ENZYME"):
            return line.replace("ENZYME", "").strip().split()
    return []


# Valid EC numbers
_EC_RE = re.compile(r'^\d+\.\d+\.\d+\.(\d+|-)$')

 
# EC → enzyme name
def get_enzyme_name(ec, cache):
    if ec in cache:
        return cache[ec]

    # Skip malformed tokens (e.g. "-", empty strings, residual prefixes)
    if not _EC_RE.match(ec):
        cache[ec] = ec   
        return ec

    try:
        data = REST.kegg_get(f"ec:{ec}").read()
    except Exception as e:
        print(f"  Warning: could not fetch EC {ec}: {e}")
        cache[ec] = ec
        return ec

    for line in data.split("\n"):
        if line.startswith("NAME"):
            name = line.replace("NAME", "").strip().split(";")[0]
            cache[ec] = name
            time.sleep(0.2)
            return name

    cache[ec] = ec
    return ec


 
# Download KGML
def download_kgml(pathway_id):
    filename = os.path.join(CACHE_DIR, f"{pathway_id}.xml")

    if os.path.exists(filename):
        return filename

    url = f"https://rest.kegg.jp/get/{pathway_id}/kgml"
    r = requests.get(url)

    with open(filename, "wb") as f:
        f.write(r.content)

    time.sleep(0.3)
    return filename


 
# Build reaction-aware graph
def clear_graph_cache():
    for path in [GRAPH_CACHE, ENZYME_CACHE]:
        if os.path.exists(path):
            os.remove(path)
            print(f"Deleted cache: {path}")

def build_graph(organism):
    if os.path.exists(GRAPH_CACHE):
        with open(GRAPH_CACHE, "rb") as f:
            return pickle.load(f)

    pathways_data = REST.kegg_list("pathway", organism).read()
    pathway_ids = [line.split("\t")[0]
                   for line in pathways_data.strip().split("\n")]

    G = nx.DiGraph()

    for pid in pathway_ids:
        kgml_file = download_kgml(pid)
        tree = ET.parse(kgml_file)
        root = tree.getroot()

        for reaction in root.findall("reaction"):
            raw_ids = reaction.attrib["name"].split()
            clean_ids = [r.replace("rn:", "").strip() for r in raw_ids if r.strip()]
            if not clean_ids:
                continue
            rid = clean_ids[0]

            substrates = [s.attrib["name"].replace("cpd:", "")
                          for s in reaction.findall("substrate")]
            products = [p.attrib["name"].replace("cpd:", "")
                        for p in reaction.findall("product")]

            reversible = reaction.attrib.get("type") == "reversible"

            for s in substrates:
                for p in products:
                    G.add_edge(s, p, reaction=rid)
                    if reversible:
                        G.add_edge(p, s, reaction=rid)

    with open(GRAPH_CACHE, "wb") as f:
        pickle.dump(G, f)

    return G


 
# Pretty path output
def format_path(G, path):
    if os.path.exists(ENZYME_CACHE):
        with open(ENZYME_CACHE, "rb") as f:
            enzyme_cache = pickle.load(f)
    else:
        enzyme_cache = {}

    output_lines = []

    for i in range(len(path) - 1):
        s = path[i]
        p = path[i + 1]

        rid = G[s][p]["reaction"]
        ecs = get_reaction_enzymes(rid)
        enzyme_names = [get_enzyme_name(ec, enzyme_cache) for ec in ecs]

        s_name = compound_id_to_name(s)
        p_name = compound_id_to_name(p)

        enzyme_str = enzyme_names[0] if enzyme_names and enzyme_names[0] else "unknown enzyme"
        output_lines.append([s_name, enzyme_str, p_name])

    with open(ENZYME_CACHE, "wb") as f:
        pickle.dump(enzyme_cache, f)

    return output_lines


 
# Single-target path (unchanged)
 
def find_path(start, target, organism):
    if not start.startswith("C"):
        start = name_to_kegg_id(start)
    if not target.startswith("C"):
        target = name_to_kegg_id(target)

    G = build_graph(organism_name_to_code(organism))
    path = nx.shortest_path(G, start, target)
    return format_path(G, path)


 
# Multi-target path finder
def find_paths(start, targets, organism):
    start_id = start if start.startswith("C") else name_to_kegg_id(start)

    org_code = organism_name_to_code(organism)

    G = build_graph(org_code)

    all_paths = []

    for target in targets:
        target_id = target if target.startswith("C") else name_to_kegg_id(target)

        try:
            path = nx.shortest_path(G, start_id, target_id)
            formatted = format_path(G, path)
            all_paths.append(formatted)
            print(f"  ✓ Path found to '{target}' ({len(formatted)} steps)")
        except nx.NetworkXNoPath:
            print(f"  ✗ No path found from '{start}' to '{target}' — skipping.")
            all_paths.append([])
        except nx.NodeNotFound as e:
            print(f"  ✗ Node not found: {e} — skipping '{target}'.")
            all_paths.append([])

    return all_paths


 
# CLI
if __name__ == "__main__":
    mode = input("Mode — (1) single target  (2) multiple targets: ").strip()

    if mode == "1":
        start   = input("Start compound: ")
        target  = input("Target compound: ")
        organism = input("Organism (e.g. Homo sapiens): ")
        result = find_path(start, target, organism)
        print("\nPath:\n", result)

    elif mode == "2":
        start    = input("Start compound: ")
        raw      = input("Target compounds (comma-separated): ")
        targets  = [t.strip() for t in raw.split(",")]
        organism = input("Organism (e.g. Homo sapiens): ")
        results  = find_paths(start, targets, organism)
        for i, (t, r) in enumerate(zip(targets, results)):
            print(f"\n--- Path to {t} ---")
            for step in r:
                print(step)
