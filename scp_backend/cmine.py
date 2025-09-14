import os
import time
import psutil
from itertools import combinations

class CMine:
    def __init__(self, database, min_support, min_coverage, max_overlap):
        """
        Parameters:
            database: list of sets; each set is a transaction
            min_support: minimum number of transactions an itemset must appear in
            min_coverage: minimum fraction of database covered by itemset (GTPC)
            max_overlap: maximum allowed overlap with previously accepted patterns
        """
        self.database = database
        self.min_support = min_support
        self.min_coverage = min_coverage
        self.max_overlap = max_overlap
        self.num_transactions = len(database)

        self.L = []              # final patterns
        self.prev_patterns = [] # to check overlap

    def get_frequent_1_itemsets(self):
        item_counts = {}
        for transaction in self.database:
            for item in transaction:
                item_counts[item] = item_counts.get(item, 0) + 1

        frequent_items = []
        for item, count in item_counts.items():
            if count >= self.min_support:
                frequent_items.append((frozenset([item]), count))
        return frequent_items

    def support(self, itemset):
        count = 0
        for transaction in self.database:
            if itemset.issubset(transaction):
                count += 1
        return count

    def coverage(self, itemset):
        covered = set()
        for idx, transaction in enumerate(self.database):
            if itemset.issubset(transaction):
                covered.add(idx)
        return len(covered) / self.num_transactions

    def overlap(self, itemset):
        overlap_count = 0
        for prev, _ in self.prev_patterns:
            if itemset & prev:
                overlap_count += 1
        return overlap_count / len(self.database)

    def join_step(self, prev_level_itemsets):
        candidates = []
        n = len(prev_level_itemsets)
        seen = set()
        for i in range(n):
            for j in range(i + 1, n):
                a, _ = prev_level_itemsets[i]
                b, _ = prev_level_itemsets[j]
                union = a | b
                if len(union) == len(a) + 1 and union not in seen:
                    seen.add(union)
                    count = self.support(union)
                    if count >= self.min_support:
                        candidates.append((union, count))
        return candidates

    def run(self):
        start_time = time.time()

        level = 1
        current_itemsets = self.get_frequent_1_itemsets()

        while current_itemsets:
            next_level_itemsets = []

            for itemset, count in current_itemsets:
                cov = self.coverage(itemset)
                ovl = self.overlap(itemset)

                if cov >= self.min_coverage and ovl <= self.max_overlap:
                    self.L.append((itemset, cov))
                    self.prev_patterns.append((itemset, cov))
                else:
                    next_level_itemsets.append((itemset, count))

            current_itemsets = self.join_step(next_level_itemsets)
            level += 1

        end_time = time.time()
        process = psutil.Process(os.getpid())

        self.execution_time = end_time - start_time
        self.memory_uss = process.memory_full_info().uss
        self.memory_rss = process.memory_info().rss

        # Sort patterns by coverage and then by length descending
        self.L = sorted(self.L, key=lambda x: (x[1], len(x[0])), reverse=True)

        return self.L
