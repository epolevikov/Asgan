

class DisjointSet:
    def __init__(self, size):
        self.parents = [i for i in range(size)]

    def find(self, i):
        if i != self.parents[i]:
            self.parents[i] = self.find(self.parents[i])
        return self.parents[i]

    def union(self, i, j):
        i_id, j_id = self.find(i), self.find(j)
        if i_id != j_id:
            self.parents[j_id] = i_id

    def get_unique_parents(self):
        return list(set(self.parents))
