from dfscode import DFSCode, DFSEdge
from graph import Graph
from custom_queue import ProbabilisticMaxQueue, ProbabilisticMinQueue
from collections import defaultdict
import heapq
import math
import random
import bisect
import numpy as np
import time
import psutil
import os

from tqdm import tqdm

class TUSM():
    def __init__(self, input_file) -> None:
        self.min_sup = 4 
        self.input_file = input_file
        self.frequentSubgraphs = []
        self.exact_alg_flag = False

    def read_graph(self, input_file):
        with open(input_file, 'r') as f:
            graph_db = dict()
            graph = None

            for line in f:
                line = line.strip()
                if len(line) == 0:  # skip empty lines
                    continue

                if line[0] == 't':
                    if graph is not None:
                        graph_db[graph.id] = graph
                    graph = Graph(int(line.split()[2]))
                elif line[0] == 'v':
                    _, id, label = line.split()
                    graph.add_vertex(int(id), label)
                elif line[0] == 'e':
                    _, frm, to, label, prob = line.split()
                    graph.add_edge(int(frm), int(to), label, float(prob))
        
            if graph is not None:
                graph_db[graph.id] = graph
        
        return graph_db


    def mine(self, k = 4, eps=0.1, delta=0.1):
        self.k = k
        graph_db = self.read_graph(self.input_file)
        return self.tkug(graph_db, graph_db.keys(), eps, delta)
        

    def find_isomers(self, dfs_code: DFSCode, graph: Graph):
        isomers = []

        if dfs_code.size == 0:
            return isomers
        
        for v in graph.vertices.values():
            if v.label == dfs_code.edges[0].v1_label:
                isomers.append({0: v.id})

        for edge in dfs_code.edges:
            new_isomers = []
            for isomer in isomers:
                curr_v = isomer[edge.v1]
                for e in graph.vertices[curr_v].edges:
                    if e.label == edge.edge_label and graph.vertices[e.to].label == edge.v2_label:
                        if edge.v1 < edge.v2 and e.to not in isomer.values():
                            new_isomers.append({**isomer, edge.v2: e.to})
                        elif edge.v1 > edge.v2 and e.to in isomer.values() and isomer[edge.v2] == e.to:
                            new_isomers.append({**isomer})
            
            isomers = new_isomers
        
        return isomers


        # for e in dfs_code.edges:
        #     for v in graph.vertices.values():
        #         if v.label == e.v1_label:
    
    def dfs_possible_extensions(self, dfs_code:DFSCode, graph:Graph):
        extensions = {}
        if dfs_code.size == 0:
            for vid, vertex in graph.vertices.items():
                for edge in vertex.edges:
                    vertex_label = vertex.label
                    other_vertex_label = graph.vertices[edge.to].label

                    if vertex_label < other_vertex_label:
                        ee = DFSEdge(0, 1, vertex_label, other_vertex_label, edge.label)
                    else:
                        ee = DFSEdge(0, 1, other_vertex_label, vertex_label, edge.label)

                    tmp = extensions.get(ee, set())
                    tmp.add(graph.id)
                    extensions[ee] = tmp
        else:
            isomers = self.find_isomers(dfs_code, graph)
            rm_node = dfs_code.right_most
            
            for isomer in isomers:
                rm_mapped = isomer[rm_node]
                inverted_map = {v: k for k, v in isomer.items()}

                # backward edge
                for edge in graph.vertices[rm_mapped].edges:
                    # check if there exist edge from rm_mapped to edge.to in the dfscode
                    # check if it falls on rightmost path
                    if edge.to not in inverted_map:
                        continue

                    inv_to = inverted_map[edge.to]
                    if not dfs_code.check_edge(inv_to, rm_node) and dfs_code.check_on_right_most_path(inv_to) and not dfs_code.is_pre_rm(inv_to):
                        ee = DFSEdge(rm_node, inv_to, graph.vertices[rm_mapped].label, graph.vertices[edge.to].label, edge.label)
                        tmp = extensions.get(ee, set())
                        tmp.add(graph.id)
                        extensions[ee] = tmp
                
                # forward edge
                for v in dfs_code.right_most_path:
                    v1_mapped = isomer[v]

                    for edge in graph.vertices[v1_mapped].edges:
                        if edge.to in inverted_map:
                            continue

                        ee = DFSEdge(v, rm_node+1, graph.vertices[v1_mapped].label, graph.vertices[edge.to].label, edge.label)
                        tmp = extensions.get(ee, set())
                        tmp.add(graph.id)
                        extensions[ee] = tmp
        return extensions
    

                
    def possible_extensions(self, graph_db, search_gids, dfs_code):
        extensions = {}
        if dfs_code.size == 0:
            for gid in search_gids:
                graph = graph_db[gid]
                for vid, vertex in graph.vertices.items():
                    for edge in vertex.edges:
                        vertex_label = vertex.label
                        other_vertex_label = graph.vertices[edge.to].label

                        if vertex_label < other_vertex_label:
                            ee = DFSEdge(0, 1, vertex_label, other_vertex_label, edge.label)
                        else:
                            ee = DFSEdge(0, 1, other_vertex_label, vertex_label, edge.label)

                        tmp = extensions.get(ee, set())
                        tmp.add(gid)
                        extensions[ee] = tmp
        else:
            for gid in search_gids:
                isomers = self.find_isomers(dfs_code, graph_db[gid])
                rm_node = dfs_code.right_most

                for isomer in isomers:
                    rm_mapped = isomer[rm_node]
                    inverted_map = {v: k for k, v in isomer.items()}

                    # backward edge
                    for edge in graph_db[gid].vertices[rm_mapped].edges:
                        # check if there exist edge from rm_mapped to edge.to in the dfscode
                        # check if it falls on rightmost path
                        if edge.to not in inverted_map:
                            continue

                        inv_to = inverted_map[edge.to]
                        if not dfs_code.check_edge(inv_to, rm_node) and dfs_code.check_on_right_most_path(inv_to) and not dfs_code.is_pre_rm(inv_to):
                            ee = DFSEdge(rm_node, inv_to, graph_db[gid].vertices[rm_mapped].label, graph_db[gid].vertices[edge.to].label, edge.label)
                            tmp = extensions.get(ee, set())
                            tmp.add(gid)
                            extensions[ee] = tmp
                    
                    # forward edge
                    for v in dfs_code.right_most_path:
                        v1_mapped = isomer[v]

                        for edge in graph_db[gid].vertices[v1_mapped].edges:
                            if edge.to in inverted_map:
                                continue

                            ee = DFSEdge(v, rm_node+1, graph_db[gid].vertices[v1_mapped].label, graph_db[gid].vertices[edge.to].label, edge.label)
                            tmp = extensions.get(ee, set())
                            tmp.add(gid)
                            extensions[ee] = tmp
        
        return extensions 
        
    
    def is_min(self, dfs_code: DFSCode):
        dfs_graph = dfs_code.to_graph()
        tmp_dfs_code = DFSCode()

        

        for edge in dfs_code.edges:
            ext = self.dfs_possible_extensions(tmp_dfs_code, dfs_graph)
            min_ee = None

            for ee, gids in ext.items():
                if min_ee == None:
                    min_ee = ee
                elif min_ee > ee:
                    min_ee = ee
            
            if min_ee is None or min_ee != edge:
                return False
            
            tmp_dfs_code.add(min_ee)

        return True
    
    def construct_dnf_formula(self, dfs_code, embeddings):
        formula = []
        for embedding in embeddings:
            clause_edges = []
            for dfs_code_edge in dfs_code.edges:
                v1 = embedding[dfs_code_edge.v1]
                v2 = embedding[dfs_code_edge.v2]

                if v1 > v2:
                    v1, v2 = v2, v1
                
                clause_edges.append((v1, v2))
 
            formula.append(clause_edges)
        return formula
    
    def sample_assignment_satisfying_clause(self, forced_clause, formula, graph):
        """
        Sample a truth assignment for all edges that appear in any clause,
        forcing the edges in 'forced_clause' to be True.
        Returns:
           assignment: a dict mapping each edge to True/False.
           all_edges: the set of all edges (as tuples) that appear in the formula.
        """
        all_edges = set()
        for clause in formula:
            all_edges.update(clause)

        assignment = {}
        for edge in all_edges:
            if edge in forced_clause:
                assignment[edge] = True
            else:
                v1, v2 = edge if edge[0] < edge[1] else (edge[1], edge[0])
                p = graph.edges[(v1, v2)].eprob
                assignment[edge] = (random.random() < p)
        return assignment, all_edges

    def fpras_for_dnf(self, formula, eps, delta, graph):
        m = len(formula)
        if m == 0:
            return 0.0

        # ----- 1.  Pre-compute clause probabilities  -----------------------------
        clause_probs = []
        for clause in formula:
            p_c = 1.0
            for e in clause:
                v1, v2 = e if e[0] < e[1] else (e[1], e[0])
                p_c *= graph.edges[(v1, v2)].eprob
            clause_probs.append(p_c)

        Z = sum(clause_probs)                # normalising constant
        if Z == 0.0:                         # all clauses impossible
            return 0.0

        # Prepare cumulative distribution for O(log m) weighted sampling
        cumw = []
        acc = 0.0
        for p in clause_probs:
            acc += p / Z                     # divide by Z to normalise
            cumw.append(acc)

        # ----- 2.  Main Monte-Carlo loop ----------------------------------------
        N = math.ceil(4.0 * m * math.log(2.0 / delta) / (eps * eps))
        X = Y = 0.0

        for _ in range(N):
            # (a) pick clause i with probability p_i/Z
            r = random.random()
            idx = bisect.bisect_left(cumw, r)
            forced_clause = formula[idx]

            # (b) sample an assignment conditional on forced_clause being true
            assignment, all_edges = self.sample_assignment_satisfying_clause(
                forced_clause, formula, graph)

            # (c) compute Pr[assignment]
            pr_pi = 1.0
            for e in all_edges:
                v1, v2 = e if e[0] < e[1] else (e[1], e[0])
                p = graph.edges[(v1, v2)].eprob
                pr_pi *= p if assignment[e] else (1.0 - p)
            Y += pr_pi

            # (d) accept only if chosen clause is the first satisfied one
            first = True
            for j in range(idx):
                if all(assignment.get(e, False) for e in formula[j]):
                    first = False
                    break
            if first:
                X += pr_pi

        return 0.0 if Y == 0.0 else (X * Z) / Y

    def compute_dnf_probability(self, formula, graph):
        m = len(formula)
        prob = 0.0
        # Iterate over all nonempty subsets of clauses.
        # We use a bitmask from 1 to (2^m)-1.
        for mask in range(1, 1 << m):
            union_edges = set()
            bits = 0
            for j in range(m):
                if mask & (1 << j):
                    bits += 1
                    # Add all edges from clause j.
                    # (Assuming formula[j] is a list of edge tuples.)
                    union_edges.update(formula[j])
            # For the union of edges, compute the product of edge probabilities.
            p = 1.0
            for edge in union_edges:
                # Ensure the edge is in canonical order.
                v1, v2 = edge if edge[0] < edge[1] else (edge[1], edge[0])
                p *= graph.edges[(v1, v2)].eprob
            # Inclusion-Exclusion: add if odd number of clauses; subtract if even.
            if bits % 2 == 1:
                prob += p
            else:
                prob -= p
        return prob



    def exact_occ_prob(self, dfs_code, graph, embeddings):
        self.exact_alg_flag = True
        formula = self.construct_dnf_formula(dfs_code, embeddings)
        return self.compute_dnf_probability(formula, graph)
    
    def approx_occ_prob(self, dfs_code, graph, embeddings, eps, delta):
        formula = self.construct_dnf_formula(dfs_code, embeddings)
        p_est = self.fpras_for_dnf(formula, eps, delta, graph)
        low = max(0.0, p_est - eps)
        high = min(1.0, p_est + eps)
        return (low, high)
    
    def approx_exp_sup(self, dfs_code, graph_db, search_ids, eps, delta):
        lower = 0
        upper = 0
        n = len(search_ids)
        for gid in search_ids:
            # breakpoint()
            graph = graph_db[gid]
            embeddings = self.find_isomers(dfs_code, graph)
            # embeddings = embedding_map[gid]
            if len(embeddings) > 0:
                lhs =  (len(embeddings) * math.log(2)) - math.log(len(embeddings))
                rhs = math.log(4 * math.log(2/delta)/(eps * eps))
                if lhs >= rhs:
                    alpha, beta = self.approx_occ_prob(dfs_code, graph, embeddings, eps, delta)
                else:
                    alpha = beta = self.exact_occ_prob(dfs_code, graph, embeddings)
            else:
                alpha = beta = 0
            lower += alpha
            upper += beta
        
        return (lower/n, upper/n)
    
    def get_mu_sigma_gaussian(self, a, b):
        mu = (a + b) / 2
        sigma = (mu - a) / 3
        return mu, sigma

    def kl_gaussian(self, mu_p, sigma_p, mu_q, sigma_q):
        return np.log(sigma_q / sigma_p) + (sigma_p**2 + (mu_p - mu_q)**2) / (2 * sigma_q**2) - 0.5

    def priority_function(self, item):
        if item[0] == item[1]:
            return 1 - item[0]
        mu_o, sigma_o = self.get_mu_sigma_gaussian(0.99, 1.01)
        mu_p, sigma_p = self.get_mu_sigma_gaussian(item[0], item[1])
        return self.kl_gaussian(mu_p, sigma_p, mu_o, sigma_o)
    
    def compare_sup(self, new_code, existing_code):
        new_code_val = self.priority_function(new_code)
        existing_code_val = self.priority_function(existing_code)
        return new_code_val < existing_code_val
            
        
    def tkug(self, graph_db, search_gids, eps, delta):
        max_queue = ProbabilisticMaxQueue()
        min_queue = ProbabilisticMinQueue()

        # process = psutil.Process(os.getpid())
        peak_memory = 0

        minsup_l = -0.1
        minsup_u = 0.1

        max_queue.push((0, 0, id(DFSCode()),DFSCode()))

        pbar = tqdm(desc="Processing max_queue")

        while max_queue:
            # current_mem = process.memory_info().rss / (1024 * 1024)  # MB
            # peak_memory = max(peak_memory, current_mem)
            _, _, _, dfs_code = max_queue.pop()
            extensions = self.possible_extensions(graph_db, search_gids, dfs_code)
            for extension, gids in tqdm(extensions.items()):
                new_dfs_code = DFSCode()
                new_dfs_code.edges = dfs_code.edges.copy()
                new_dfs_code.size = dfs_code.size
                new_dfs_code.right_most = dfs_code.right_most
                new_dfs_code.right_most_path = dfs_code.right_most_path.copy()
                new_dfs_code.add(extension)

                esup_l, esup_u = self.approx_exp_sup(new_dfs_code, graph_db, search_gids, eps, delta)

                if self.compare_sup((esup_l, esup_u), (minsup_l, minsup_u)) and self.is_min(new_dfs_code):
                    min_queue.push((esup_l, esup_u, id(new_dfs_code), new_dfs_code))

                    if len(min_queue) > self.k:
                        min_queue.pop()
                        minsup_l = min_queue.peek()[0]
                        minsup_u = min_queue.peek()[1]
                    elif len(min_queue) == self.k:
                        minsup_l = min_queue.peek()[0]
                        minsup_u = min_queue.peek()[1]
                    max_queue.push((esup_l, esup_u, id(new_dfs_code), new_dfs_code))
            pbar.update(1)
        pbar.close()
        return min_queue, peak_memory              

if __name__ == "__main__":
    # ----- Dummy dataset directly in the code -----
    dummy_data = """
t # 0
v 0 1
v 1 2
v 2 2
v 3 7
e 0 3 3 0.7
e 1 2 1 0.4
e 1 3 1 0.9
e 2 3 2 0.2

t # 1
v 0 0
v 1 1
v 2 2
v 3 5
v 4 5
e 0 1 3 0.4
e 0 2 2 0.5
e 1 2 1 0.7
e 1 3 1 0.6
e 2 4 1 0.2

t # 2
v 0 0
v 1 1
v 2 2
v 3 2
e 0 1 2 0.3
e 0 2 1 0.4
e 0 3 1 0.2
e 1 2 1 0.8
e 1 3 1 0.9

t # 3
v 0 0
v 1 1
v 2 2
v 3 2
v 4 7
e 0 1 2 0.5
e 0 2 1 0.6
e 1 2 1 0.1
e 2 3 1 0.3
e 2 4 2 0.4
e 3 4 1 0.5

t # 4
v 0 1
v 1 2
v 2 2
v 3 3
v 4 6
v 5 0
e 0 1 1 0.2
e 0 2 4 0.3
e 0 3 1 0.4
e 0 5 2 0.5
e 2 3 2 0.6
e 2 5 1 0.7
e 3 4 3 0.8

t # 5
v 0 1
v 1 2
v 2 3
v 3 6
e 0 2 1 0.9
e 0 3 1 0.8
e 1 2 2 0.3
e 2 3 3 0.4

t # 6
v 0 5
v 1 2
v 2 1
v 3 0
v 4 0
v 5 2
e 0 1 1 0.1
e 1 2 1 0.2
e 2 3 3 0.2
e 2 4 2 0.4
e 4 5 1 0.9
"""

    # Write dummy data to a temporary file
    tmp_file = "dummy_uncertain_graph.txt"
    with open(tmp_file, "w") as f:
        f.write(dummy_data.strip())

    # ----- Initialize TUSM with the dummy file -----
    g = TUSM(tmp_file)

    # ----- Set parameters -----
    k = 3       # top-3 subgraphs
    eps = 0.1   # approximation error
    delta = 0.3 # failure probability

    # ----- Run mining -----
    min_queue, _ = g.mine(k, eps, delta)

    # ----- Print top-k subgraphs -----
    print(f"Top {k} frequent subgraphs:")
    while not min_queue.is_empty():
        _, _, _, dfs_code = min_queue.pop()
        print(dfs_code)


