"""The 'graphos' NetworkX backend interface (Level 3).

Algorithms implemented natively run on the C++ engine through the view's
construction data; everything else is declined with NotImplementedError,
which lets NetworkX fall back to its own implementation through
convert_to_nx (enable nx.config.fallback_to_nx)."""

from .. import _core
from . import AdjacencyGraph, HasseGraph, IncidenceGraph


def _native_component_labels(G):
    if isinstance(G, AdjacencyGraph):
        if G._exclude is not None:
            return _core.connected_components(G._c, G._k, G._via, G._exclude)
        return _core.connected_components(G._c, G._k, G._via)
    raise NotImplementedError("graphos: native components need an AdjacencyGraph view")


class GraphosBackend:
    @staticmethod
    def convert_from_nx(graph, *args, **kwargs):
        # our views enter as-is; genuine nx graphs are not converted into
        # graphos (construction is the builder API's job)
        return graph

    @staticmethod
    def convert_to_nx(obj, *args, **kwargs):
        if isinstance(obj, (AdjacencyGraph, HasseGraph, IncidenceGraph)):
            return obj.to_networkx()
        return obj

    # ---- native algorithms ---------------------------------------------
    @staticmethod
    def connected_components(G):
        labels = _native_component_labels(G)
        arr = labels.label
        comps = [set() for _ in range(labels.count)]
        for cell, lab in enumerate(arr):
            comps[int(lab)].add(cell)
        return iter(comps)

    @staticmethod
    def number_connected_components(G):
        return int(_native_component_labels(G).count)

    @staticmethod
    def is_connected(G):
        return int(_native_component_labels(G).count) == 1
