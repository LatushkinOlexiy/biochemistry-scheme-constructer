import os
import re
import time
import pickle
import requests
import networkx as nx
import xml.etree.ElementTree as ET
import translate as tr


KEGG_BASE = "https://rest.kegg.jp"

def _kegg_get(endpoint, retries=3, delay=1.0):
    url = f"{KEGG_BASE}/{endpoint}"
    headers = {"User-Agent": "KEGG-Pathway-Drawer/1.0 (academic use)"}
    for attempt in range(retries):
        try:
            r = requests.get(url, headers=headers, timeout=30)
            if r.status_code == 200:
                return r.text
            if r.status_code in (400, 404):
                # Don't retry client errors — they won't improve
                raise ValueError(
                    f"KEGG returned HTTP {r.status_code} for: {url}"
                )
            # 5xx or other — wait and retry
            time.sleep(delay * (attempt + 1))
        except requests.RequestException as e:
            if attempt == retries - 1:
                raise
            time.sleep(delay * (attempt + 1))
    raise ValueError(f"KEGG request failed after {retries} attempts: {url}")

# Known short codes that can be passed directly, skipping the API lookup
_ORG_CODE_RE = re.compile(r'^[a-z]{2,6}$')

def organism_name_to_code(user_input):

    user_input = user_input.strip()

    # If it looks like a code already, validate it and return
    if _ORG_CODE_RE.match(user_input):
        try:
            # A valid code returns pathway data; invalid returns 400/404
            _kegg_get(f"list/pathway/{user_input}")
            print(f"Organism code accepted: {user_input}")
            return user_input
        except ValueError:
            pass  # not a valid code — fall through to name search

    # Search by name using the find/organism endpoint
    try:
        data = _kegg_get(f"find/organism/{requests.utils.quote(user_input)}")
    except ValueError as e:
        raise ValueError(f"Organism not found: '{user_input}'") from e

    lines = [l for l in data.strip().split("\n") if l.strip()]
    if not lines:
        raise ValueError(f"No organism matched: '{user_input}'")

    # Each line: T-number <TAB> code <TAB> name <TAB> taxonomy
    # Pick the first (best) match
    parts = lines[0].split("\t")
    if len(parts) < 2:
        raise ValueError(f"Unexpected organism response format: {lines[0]}")

    code = parts[1].strip()
    name = parts[2].strip() if len(parts) > 2 else code
    print(f"Matched organism: {name} → {code}")
    return code


CACHE_DIR = "kgml_cache"
GRAPH_CACHE = "compound_graph_with_reactions.pkl"
ENZYME_CACHE = "enzyme_cache.pkl"

os.makedirs(CACHE_DIR, exist_ok=True)


# -----------------------------
# Name → KEGG compound ID
# -----------------------------
def name_to_kegg_id(name):
    result = _kegg_get(f"find/compound/{requests.utils.quote(name)}")

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


# -----------------------------
# Compound ID → name
# -----------------------------
def compound_id_to_name(cid):
    data = _kegg_get(f"get/cpd:{cid}")
    for line in data.split("\n"):
        if line.startswith("NAME"):
            return line.replace("NAME", "").strip().split(";")[0]
    return cid


# -----------------------------
# Reaction → EC numbers
# -----------------------------
def get_reaction_enzymes(rid):
    # Guard: if rid somehow still contains spaces or rn: prefix
    # (e.g. from an old cache), clean it before forming the URL.
    rid = rid.split()[0].replace("rn:", "").strip()
    if not rid:
        return []
    try:
        data = _kegg_get(f"get/rn:{rid}")
    except Exception as e:
        print(f"  Warning: could not fetch reaction {rid}: {e}")
        return []
    for line in data.split("\n"):
        if line.startswith("ENZYME"):
            return line.replace("ENZYME", "").strip().split()
    return []


# Valid EC numbers: digits and optional "-" in four dot-separated fields
_EC_RE = re.compile(r'^\d+\.\d+\.\d+\.(\d+|-)$')

# -----------------------------
# EC → enzyme name
# -----------------------------
def get_enzyme_name(ec, cache):
    if ec in cache:
        return cache[ec]

    # Skip malformed tokens (e.g. "-", empty strings, residual prefixes)
    if not _EC_RE.match(ec):
        cache[ec] = ec   # store as-is so we don't retry every time
        return ec

    try:
        data = _kegg_get(f"get/ec:{ec}")
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


# -----------------------------
# Download KGML
# -----------------------------
def download_kgml(pathway_id):
    filename = os.path.join(CACHE_DIR, f"{pathway_id}.xml")

    if os.path.exists(filename):
        return filename

    try:
        text = _kegg_get(f"get/{pathway_id}/kgml")
    except ValueError as e:
        print(f"  Warning: could not download KGML for {pathway_id}: {e}")
        return None   # caller must handle None

    with open(filename, "w", encoding="utf-8") as f:
        f.write(text)

    time.sleep(0.3)
    return filename


# -----------------------------
# Build reaction-aware graph
# -----------------------------
def clear_graph_cache():

    for path in [GRAPH_CACHE, ENZYME_CACHE]:
        if os.path.exists(path):
            os.remove(path)
            print(f"Deleted cache: {path}")

def build_graph(organism):
    if os.path.exists(GRAPH_CACHE):
        with open(GRAPH_CACHE, "rb") as f:
            return pickle.load(f)

    pathways_data = _kegg_get(f"list/pathway/{organism}")
    pathway_ids = [line.split("\t")[0]
                   for line in pathways_data.strip().split("\n")]

    G = nx.DiGraph()

    for pid in pathway_ids:
        kgml_file = download_kgml(pid)
        if kgml_file is None:
            continue   # skip pathways whose KGML couldn't be fetched
        try:
            tree = ET.parse(kgml_file)
        except ET.ParseError as e:
            print(f"  Warning: could not parse KGML for {pid}: {e}")
            continue
        root = tree.getroot()

        for reaction in root.findall("reaction"):
            # name attribute can be "rn:R01278 rn:R01279" (space-separated)
            # Strip the rn: prefix from every token, keep the first valid ID
            raw_ids = reaction.attrib["name"].split()
            clean_ids = [r.replace("rn:", "").strip() for r in raw_ids if r.strip()]
            if not clean_ids:
                continue
            rid = clean_ids[0]   # use the first reaction ID for enzyme lookup

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


# -----------------------------
# Pretty path output
# -----------------------------
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


# -----------------------------
# Single-target path (unchanged)
# -----------------------------
def find_path(start, target, organism):
    if not start.startswith("C"):
        start = name_to_kegg_id(start)
    if not target.startswith("C"):
        target = name_to_kegg_id(target)

    G = build_graph(organism_name_to_code(organism))
    path = nx.shortest_path(G, start, target)
    return format_path(G, path)


# -----------------------------
# NEW: Multi-target path finder
# -----------------------------
def find_paths(start, targets, organism):
    # Resolve start compound ID once
    start_id = start if start.startswith("C") else name_to_kegg_id(start)

    # Resolve organism code once
    org_code = organism_name_to_code(organism)

    # Build (or load) graph once — expensive, share across all targets
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
            all_paths.append([])   # empty list keeps index alignment with targets
        except nx.NodeNotFound as e:
            print(f"  ✗ Node not found: {e} — skipping '{target}'.")
            all_paths.append([])

    return all_paths


# -----------------------------
# CLI
# -----------------------------
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