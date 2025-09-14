from smiles_conversion import smiles_to_gtd
from pathlib import Path
import gspan as gsp
from PAMI.extras.stats import graphDatabase as gdb
from contextlib import redirect_stdout
from cmine import CMine
from visualize import visualize_scps
from pathlib import Path

def read_gspan_output(file_path):
    """
    Reads the gSpan output where each line represents a transaction
    and contains subgraph IDs separated by spaces.
    """
    transactions = []
    with open(file_path, 'r') as f:
        for line in f:
            items = set(line.strip().split())
            transactions.append(items)
    return transactions

def run_algorithm(file_path: str, algorithm: str, min_support: float, max_overlap: float, min_coverage: float):
    extracted_string = file_path.split("\\", 1)[-1]
    file_path = Path("uploads") / f"{extracted_string}"
    smiles_to_gtd(file_path, 'gtd.txt')
    
    dpath = 'gtd.txt'
    obj = gdb.graphDatabase(iFile=dpath)

    with open('stats.txt', 'w') as file:
        with redirect_stdout(file):
            obj.printGraphDatabaseStatistics()

    with open('stats.txt', 'r') as file:
        file_content = file.read().split("\n")

    obj = gsp.GSpan(dpath, min_support, False, float("inf"), False)
    obj.mine()

    obj.saveSubgraphsByGraphId("temp.txt")    
    obj.save('frequentSubgraphs.txt')
    
    frequent_subgraphs = obj.frequentSubgraphs

    num = len(obj.frequentSubgraphs)
    file_content.append(f"number of frequent graphs: {num}")
    
    gspan_output_file = "temp.txt"

    transactions = read_gspan_output(gspan_output_file)

    cm = CMine(transactions, min_support, min_coverage, max_overlap)
    patterns = cm.run()
    
    results_path = Path("outputs") / "cmine_results.txt"
    results_path.parent.mkdir(parents=True, exist_ok=True)
    
    with open(results_path, 'w') as out_file:
        out_file.write("Frequent Subgraph Patterns (CMine Output):\n")
        for itemset, coverage in patterns:
            line = f"Subgraph IDs: {set(itemset)} | Coverage: {coverage:.2f}\n"
            out_file.write(line)
            
    visualize_scps('frequentSubgraphs.txt', Path("outputs") / "cmine_results.txt")
    