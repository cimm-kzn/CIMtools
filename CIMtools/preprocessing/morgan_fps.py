from collections import defaultdict
from numpy import zeros
from sklearn.base import BaseEstimator, TransformerMixin


class MorganFPs(BaseEstimator, TransformerMixin):
    def __init__(self, fingerprint_size=1024, bits_active=2):
        self.fingerprint_size = fingerprint_size
        self.bits_active = bits_active

    def cgr_info(self, cgr, adj, hash_set, radius):
        dict_radius = defaultdict(tuple)
        for n, (k, v) in enumerate(adj.items()):
            values = []
            values.append((radius, hash_set[k]))
            for kk, vv in v.items():
                values.append((vv, hash_set[kk]))
                if radius == 2:
                    for i, j in adj[kk].items():
                        values.append((j, hash_set[i]))
            dict_radius[k] = tuple(values)
        for n, a in cgr.atoms():
            if (n in dict_radius.keys()) is False:
                dict_radius[n] = hash_set[n]
            else:
                dict_radius[n] = hash(dict_radius[n])
        if radius == 0:
            dict_radius = hash_set
        info = {v: None for k, v in dict_radius.items()}
        for v in info:
            lst = []
            for n in dict_radius:
                if v == dict_radius[n]:
                    lst.append((n, radius))
            info[v] = tuple(lst)
        return info

    def cgr2sentence(self, cgr, radius):
        radii = list(range(int(radius) + 1))
        adj = defaultdict(dict)
        for n, m, b in cgr.bonds():
            adj[n][m] = adj[m][n] = int(b)

        h = defaultdict(tuple)
        for n, (k, v) in enumerate(adj.items()):
            values = []
            for kk, vv in v.items():
                values.append((vv))
            h[k] = tuple(values)

        dict_rings = {n: 0 for n, a in cgr.atoms()}
        add = []
        for lst in cgr.aromatic_rings:
            for n in lst:
                add.append(n)
            for n in cgr:
                if n in add:
                    dict_rings[n] = 1
                else:
                    dict_rings[n] = 0

        nodes = {n: (int(a), a.atomic_number, a.atomic_mass, a.neighbors, a.p_neighbors,
                     a.charge, a.p_charge, a.hybridization, a.p_hybridization,
                     dict_rings[n]) for n, a in cgr.atoms()}
        hash_set = {n: hash(i) for n, i in nodes.items()}

        info = defaultdict(tuple)
        for r in radii:
            dict_r = self.cgr_info(cgr, adj, hash_set, r)
            info.update(dict_r)

        cgr_atoms = [n for n, a in cgr.atoms()]
        dict_atoms = {x: {r: None for r in radii} for x in cgr_atoms}

        for element in info:
            for atom_idx, radius_at in info[element]:
                dict_atoms[atom_idx][radius_at] = element  # {atom number: {fp radius: identifier}}

        for key, value in dict_atoms.items():
            for k, v in value.items():
                if v is None:
                    dict_atoms[key][k] = dict_atoms[key][radius]

        # iterate over all atoms and radii
        identifier_sentences = []

        for r in radii:  # iterate over radii to get one sentence per radius
            identifiers = []
            for atom in dict_atoms:  # iterate over atoms
                # get one sentence per radius
                identifiers.append(dict_atoms[atom][r])
            identifier_sentences.append(list(map(str, [x for x in identifiers if x])))

        # merge identifiers alternating radius to sentence: atom 0 radius0, atom 0 radius 1, etc.
        identifiers_alt = []
        for atom in dict_atoms:  # iterate over atoms
            for r in radii:  # iterate over radii
                identifiers_alt.append(dict_atoms[atom][r])

        alternating_sentence = map(str, [x for x in identifiers_alt if x])

        return list(identifier_sentences), list(alternating_sentence)

    def morgan(self, cgr, radius):
        active_bits = set()
        fingerprint = zeros(self.fingerprint_size)

        for i in self.cgr2sentence(cgr, radius)[-1]:
            i = int(i)
            active_bits.add(i)

        for b in range(self.bits_active):
            for a in active_bits:
                fingerprint[(a + b) & 1023] = 1
                fingerprint[((a >> 10) + b) & 1023] = 1

        return fingerprint


__all__ = ['MorganFPs']