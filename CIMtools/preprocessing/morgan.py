from ..base import CIMtoolsTransformerMixin


class MorganFingerprint(CIMtoolsTransformerMixin):
    """
    Implementation of Extended-Connectivity Fingerprints
    """
    __slots__ = ('_atoms', '_bonds', '_radius', '_length', '_counts')

    def __init__(self, radius=2, length=1024, counts=False):
        self._atoms = {}
        self._bonds = {}
        self._radius = radius
        self._length = length
        self._counts = counts

    def transform(self, x):
        x = super().transform(x)[0]

        atoms = {k: v for k, v in x.atoms()}
        bonds = {k: {} for k in atoms}
        for k, v, b in x.bonds():
            bonds[k][v] = bonds[v][k] = b.order

        old = {k: hash(tuple([atom.neighbors, atom.total_hydrogens, atom.atomic_number, atom.atomic_mass,
                                atom.charge, atom.hybridization, atom.in_ring, *sorted(bonds[k].values())]))
               for k, atom in atoms.items()}

        all_ = [x for x in old.values()]
        new = {}

        for _ in range(self._radius):
            for n, a in old.items():
                neibs = [(a, b) for a, b in bonds[n].items()]
                new[n] = [a]
                for pair in neibs:
                    new[n].extend(pair)
                new[n] = hash(tuple(new[n]))

            all_.extend(x for x in new.values() if x not in all_)
            if new == old:
                break

            old, new = new, {}

        fingerprint = [0] * self._length
        for a in old.values():
            fingerprint[a & self._length - 1] = 1
            fingerprint[a >> 10 & self._length - 1] = 1

        return fingerprint


__all__ = ['MorganFingerprint']
