import csv
import re
from collections import defaultdict, deque

class Graph:
    def __init__(self):
        """Initialize an empty adjacency list graph."""
        self.graph = defaultdict(set)

    def add_edge(self, node1, node2):
        self.graph[node1].add(node2)
        self.graph[node2].add(node1)

    def load_from_csv(self, file_path):
        with open(file_path, "r") as file:
            reader = csv.reader(file)
            next(reader)  

            for row in reader:
                if len(row) < 3:
                    continue  

                locus, allele, disease = row[:3]  
                self.add_edge(locus.strip(), allele.strip()) 
                self.add_edge(allele.strip(), disease.strip())  

def is_disease_node(node):
    # Filter out mutations, gene names, numbers, and meaningless symbols
    if node.startswith("m."):
        return False
    if "MT-" in node or node.strip() == "-" or node.strip() == "":
        return False
    if any(char.isdigit() for char in node):
        return False
    return True

def bfs_locus_to_disease(graph, start_node, max_depth=3):
    if start_node not in graph:
        print(f"❌ {start_node} not found in the graph.")
        return []

    queue = deque([[start_node]]) 
    paths = []  
    visited = set() 

    print(f"🔎 BFS starting from: {start_node}")

    while queue:
        path = queue.popleft()  
        node = path[-1] 
        if len(path) > max_depth:
            continue

        if is_disease_node(node) and node != start_node:
            paths.append(path)
            print(f"✅ Found Disease Path: {' → '.join(path)}")
        else:
            for neighbor in graph.get(node, []):
                if neighbor not in visited:
                    visited.add(neighbor)
                    new_path = path + [neighbor]
                    queue.append(new_path)
                    print(f"🔄 Expanding Path: {' → '.join(new_path)}")

    return paths
