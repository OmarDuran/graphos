// graphos Python bindings (graphos._core).
//
// Lifetime follows the builder/frozen split: accessors on the mutable Complex,
// and on operation results, return NumPy copies, safe against reallocation;
// accessors on the immutable FrozenComplex return zero-copy views whose base
// object keeps the complex alive.

#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <algorithm>

#include "graphos/graphos.hpp"

namespace py = pybind11;
using namespace graphos;

namespace {

template <typename T>
py::array_t<T> copy_array(const std::vector<T>& v) {
  py::array_t<T> a(static_cast<py::ssize_t>(v.size()));
  std::copy(v.begin(), v.end(), a.mutable_data());
  return a;
}

py::tuple csr_copy(const std::vector<Index>& offsets, const std::vector<Index>& indices,
                   const std::vector<Sign>& signs) {
  return py::make_tuple(copy_array(offsets), copy_array(indices), copy_array(signs));
}

std::vector<Identification> to_idents(py::handle seq) {
  std::vector<Identification> out;
  for (py::handle item : seq) {
    py::sequence t = py::reinterpret_borrow<py::sequence>(item);
    Identification id;
    id.from = t[0].cast<Index>();
    id.to = t[1].cast<Index>();
    id.rel_sign = t.size() > 2 ? static_cast<Sign>(t[2].cast<int>()) : Sign{1};
    out.push_back(id);
  }
  return out;
}

std::vector<std::vector<Identification>> to_idents_by_dim(py::handle seq) {
  std::vector<std::vector<Identification>> out;
  for (py::handle level : seq) out.push_back(to_idents(level));
  return out;
}

py::list idents_to_py(const std::vector<Identification>& ids) {
  py::list out;
  for (const Identification& id : ids) {
    out.append(py::make_tuple(id.from, id.to, static_cast<int>(id.rel_sign)));
  }
  return out;
}

void check_dim(std::size_t k, std::size_t n, const char* what) {
  if (k >= n) throw py::index_error(std::string(what) + ": dimension out of range");
}

}  // namespace

PYBIND11_MODULE(_core, m) {
  m.doc() = "graphos: a metric-free engine for finite cell complexes: stratified signed incidence, the boundary and coboundary operators, and an operation calculus that returns the induced chain maps";
  m.attr("INVALID_INDEX") = invalid_index;

  // ---- the complex and its operators -----------------------------------
  py::class_<Complex>(m, "Complex")
      .def(py::init<int>(), py::arg("dim"))
      .def("attach_vertices", &Complex::attach_vertices, py::arg("n"))
      .def(
          "attach_cell",
          [](Complex& c, int k, const std::vector<Index>& bnd, const std::vector<int>& sg) {
            std::vector<Sign> s(sg.begin(), sg.end());
            return c.attach_cell(k, bnd, s);
          },
          py::arg("k"), py::arg("boundary"), py::arg("signs"))
      .def_property_readonly("dim", &Complex::dim)
      .def("count", &Complex::count, py::arg("k"))
      .def("counts", &Complex::counts)
      .def("validate", &Complex::validate)
      .def(
          "boundary",
          [](const Complex& c, int k) {
            const BoundaryOperator& b = c.boundary(k);
            return csr_copy(b.offsets, b.indices, b.signs);
          },
          py::arg("k"), "CSR (offsets, indices, signs) of the boundary operator, as copies")
      .def("__repr__", [](const Complex& c) {
        std::string s = "Complex(dim=" + std::to_string(c.dim()) + ", counts=[";
        for (int k = 0; k <= c.dim(); ++k) {
          s += (k ? ", " : "") + std::to_string(c.count(k));
        }
        return s + "])";
      });

  m.def("euler_characteristic", &euler_characteristic);
  m.def("d_squared_is_zero", &d_squared_is_zero);

  py::class_<Marker>(m, "Marker")
      .def(py::init<const Complex&>(), py::arg("complex"))
      .def("mark", &Marker::mark, py::arg("k"), py::arg("cell"),
           py::return_value_policy::reference_internal)
      .def(
          "mark_where",
          [](Marker& mk, int k, const std::function<bool(Index)>& pred) -> Marker& {
            return mk.mark_where(k, pred);
          },
          py::arg("k"), py::arg("predicate"), py::return_value_policy::reference_internal)
      .def("marked", &Marker::marked, py::arg("k"), py::arg("cell"))
      .def("marked_count", &Marker::marked_count, py::arg("k"))
      .def_property_readonly("dim", &Marker::dim)
      .def_static("from_cells", &Marker::from_cells, py::arg("complex"), py::arg("cells"));

  py::class_<ChainMap>(m, "ChainMap")
      .def_property_readonly("dims", [](const ChainMap& cm) { return cm.index.size(); })
      .def(
          "index",
          [](const ChainMap& cm, std::size_t k) {
            check_dim(k, cm.index.size(), "ChainMap.index");
            return copy_array(cm.index[k]);
          },
          py::arg("k"))
      .def(
          "sign",
          [](const ChainMap& cm, std::size_t k) {
            check_dim(k, cm.sign.size(), "ChainMap.sign");
            return copy_array(cm.sign[k]);
          },
          py::arg("k"));
  m.def("compose", [](const ChainMap& a, const ChainMap& b) { return compose(a, b); });

  // ---- constructors -----------------------------------------------------
  m.def("from_edges", &from_edges, py::arg("n_vertices"), py::arg("segments"));
  m.def("from_polygons", &from_polygons, py::arg("n_vertices"), py::arg("polygons"));
  m.def("from_polyhedra", &from_polyhedra, py::arg("n_vertices"), py::arg("cells"));
  m.def("from_simplices", &from_simplices, py::arg("dim"), py::arg("n_vertices"), py::arg("cells"));

  // ---- derived operators ------------------------------------------------
  m.def(
      "coboundary",
      [](const Complex& c, int k) {
        const CoboundaryOperator cob = coboundary(c, k);
        return csr_copy(cob.offsets, cob.indices, cob.signs);
      },
      py::arg("complex"), py::arg("k"));
  m.def(
      "incidence",
      [](const Complex& c, int k, int j) {
        const Adjacency a = incidence(c, k, j);
        return py::make_tuple(copy_array(a.offsets), copy_array(a.indices));
      },
      py::arg("complex"), py::arg("k"), py::arg("j"));
  m.def(
      "adjacency",
      [](const Complex& c, int k, int via) {
        const Adjacency a = adjacency(c, k, via);
        return py::make_tuple(copy_array(a.offsets), copy_array(a.indices));
      },
      py::arg("complex"), py::arg("k"), py::arg("via"));

  // ---- frozen storage, as zero-copy views -------------------------------
  py::class_<FrozenComplex>(m, "FrozenComplex")
      .def_property_readonly("dim", &FrozenComplex::dim)
      .def_property_readonly("halo_depth", &FrozenComplex::halo_depth)
      .def("count", &FrozenComplex::count, py::arg("k"))
      .def("counts", &FrozenComplex::counts)
      .def(
          "boundary",
          [](py::object self_obj, int k) {
            FrozenComplex& f = self_obj.cast<FrozenComplex&>();
#if defined(GRAPHOS_HAVE_CHAI)
            // CHAI-managed storage: fall back to row-wise copies
            std::vector<Index> offsets{0};
            std::vector<Index> indices;
            std::vector<Sign> signs;
            for (Index e = 0; e < f.count(k); ++e) {
              const auto r = f.boundary_row(k, e);
              indices.insert(indices.end(), r.indices, r.indices + r.size);
              signs.insert(signs.end(), r.signs, r.signs + r.size);
              offsets.push_back(static_cast<Index>(indices.size()));
            }
            return csr_copy(offsets, indices, signs);
#else
            const CsrView v = f.boundary_view(k);
            const py::ssize_t n = f.count(k);
            const py::ssize_t nnz = v.offsets[static_cast<std::size_t>(n)];
            return py::make_tuple(py::array_t<Index>({n + 1}, {sizeof(Index)}, v.offsets, self_obj),
                                  py::array_t<Index>({nnz}, {sizeof(Index)}, v.indices, self_obj),
                                  py::array_t<Sign>({nnz}, {sizeof(Sign)}, v.signs, self_obj));
#endif
          },
          py::arg("k"), "CSR (offsets, indices, signs) of ∂_k — zero-copy views")
      .def(
          "coboundary",
          [](py::object self_obj, int k) {
            FrozenComplex& f = self_obj.cast<FrozenComplex&>();
#if defined(GRAPHOS_HAVE_CHAI)
            std::vector<Index> offsets{0};
            std::vector<Index> indices;
            std::vector<Sign> signs;
            for (Index e = 0; e < f.count(k); ++e) {
              const auto r = f.coboundary_row(k, e);
              indices.insert(indices.end(), r.indices, r.indices + r.size);
              signs.insert(signs.end(), r.signs, r.signs + r.size);
              offsets.push_back(static_cast<Index>(indices.size()));
            }
            return csr_copy(offsets, indices, signs);
#else
            const CsrView v = f.coboundary_view(k);
            const py::ssize_t n = f.count(k);
            const py::ssize_t nnz = v.offsets[static_cast<std::size_t>(n)];
            return py::make_tuple(py::array_t<Index>({n + 1}, {sizeof(Index)}, v.offsets, self_obj),
                                  py::array_t<Index>({nnz}, {sizeof(Index)}, v.indices, self_obj),
                                  py::array_t<Sign>({nnz}, {sizeof(Sign)}, v.signs, self_obj));
#endif
          },
          py::arg("k"), "CSR (offsets, indices, signs) of δ_k — zero-copy views")
      .def("star", &FrozenComplex::star, py::arg("k"), py::arg("cell"))
      .def("closure", &FrozenComplex::closure, py::arg("k"), py::arg("cell"))
      .def("link", &FrozenComplex::link, py::arg("k"), py::arg("cell"));
  m.def("freeze", &freeze, py::arg("complex"), py::arg("halo_depth") = 1);

  py::class_<DualView>(m, "DualView")
      .def_property_readonly("dim", &DualView::dim)
      .def("count", &DualView::count, py::arg("k"))
      .def(
          "boundary_row",
          [](const DualView& d, int k, Index cell) {
            const auto r = d.boundary_row(k, cell);
            return py::make_tuple(copy_array(std::vector<Index>(r.indices, r.indices + r.size)),
                                  copy_array(std::vector<Sign>(r.signs, r.signs + r.size)));
          },
          py::arg("k"), py::arg("cell"))
      .def(
          "coboundary_row",
          [](const DualView& d, int k, Index cell) {
            const auto r = d.coboundary_row(k, cell);
            return py::make_tuple(copy_array(std::vector<Index>(r.indices, r.indices + r.size)),
                                  copy_array(std::vector<Sign>(r.signs, r.signs + r.size)));
          },
          py::arg("k"), py::arg("cell"));
  m.def("dual", &dual, py::keep_alive<0, 1>(), py::arg("frozen"));

  // ---- the operation calculus -------------------------------------------
  py::class_<DisjointUnionResult>(m, "DisjointUnionResult")
      .def_readonly("complex", &DisjointUnionResult::complex)
      .def_readonly("a_map", &DisjointUnionResult::a_map)
      .def_readonly("b_map", &DisjointUnionResult::b_map);
  m.def("disjoint_union", &disjoint_union, py::arg("a"), py::arg("b"));

  py::class_<QuotientResult>(m, "QuotientResult")
      .def_readonly("complex", &QuotientResult::complex)
      .def_readonly("map", &QuotientResult::map);
  m.def(
      "quotient",
      [](const Complex& c, py::sequence by_dim) { return quotient(c, to_idents_by_dim(by_dim)); },
      py::arg("complex"), py::arg("identifications"));
  m.def(
      "find_parallel_cells",
      [](const Complex& c, int k) { return idents_to_py(find_parallel_cells(c, k)); },
      py::arg("complex"), py::arg("k"));

  py::class_<PushoutResult>(m, "PushoutResult")
      .def_readonly("complex", &PushoutResult::complex)
      .def_readonly("a_map", &PushoutResult::a_map)
      .def_readonly("b_map", &PushoutResult::b_map);
  m.def(
      "pushout",
      [](const Complex& a, const Complex& b, py::sequence glue, bool dedup) {
        return pushout(a, b, to_idents(glue), dedup);
      },
      py::arg("a"), py::arg("b"), py::arg("vertex_identifications"), py::arg("deduplicate") = true);

  m.def(
      "lift_identifications",
      [](const Complex& c, py::sequence pairs) {
        py::list out;
        for (const auto& level : lift_identifications(c, to_idents(pairs))) {
          out.append(idents_to_py(level));
        }
        return out;
      },
      py::arg("complex"), py::arg("vertex_pairs"));

  py::class_<ProductResult>(m, "ProductResult")
      .def_readonly("complex", &ProductResult::complex)
      .def_property_readonly("blocks", [](const ProductResult& r) {
        py::list out;
        for (const auto& level : r.blocks) {
          py::list row;
          for (const auto& b : level) {
            row.append(py::make_tuple(b.p, b.q, b.offset, b.a_count, b.b_count));
          }
          out.append(row);
        }
        return out;
      });
  m.def("product", &product, py::arg("a"), py::arg("b"));

  py::class_<JoinResult>(m, "JoinResult")
      .def_readonly("complex", &JoinResult::complex)
      .def_readonly("a_map", &JoinResult::a_map)
      .def_readonly("b_map", &JoinResult::b_map);
  m.def("join", &join, py::arg("a"), py::arg("b"));

  py::class_<StarDeletionResult>(m, "StarDeletionResult")
      .def_readonly("complex", &StarDeletionResult::complex)
      .def_readonly("map", &StarDeletionResult::map);
  m.def("star_deletion", &star_deletion, py::arg("complex"), py::arg("cells"));

  py::class_<SubcomplexResult>(m, "SubcomplexResult")
      .def_readonly("complex", &SubcomplexResult::complex)
      .def_readonly("map", &SubcomplexResult::map)
      .def_readonly("embedding", &SubcomplexResult::embedding);
  m.def("subcomplex", &subcomplex, py::arg("complex"), py::arg("cells"));

  py::class_<CutResult>(m, "CutResult")
      .def_readonly("complex", &CutResult::complex)
      .def_readonly("ancestor", &CutResult::ancestor);
  m.def("cut_along", &cut_along, py::arg("complex"), py::arg("interface"));

  py::class_<ReplaceResult>(m, "ReplaceResult")
      .def_readonly("complex", &ReplaceResult::complex)
      .def_readonly("map", &ReplaceResult::map)
      .def_readonly("patch_map", &ReplaceResult::patch_map);
  m.def(
      "replace",
      [](const Complex& c, const Marker& region, const Complex& patch, py::sequence glue) {
        return replace(c, region, patch, to_idents(glue));
      },
      py::arg("complex"), py::arg("region"), py::arg("patch"), py::arg("vertex_glue"));

  py::class_<OrientationResult>(m, "OrientationResult")
      .def_readonly("complex", &OrientationResult::complex)
      .def_readonly("map", &OrientationResult::map)
      .def_readonly("orientable", &OrientationResult::orientable);
  m.def("orient", &orient, py::arg("complex"));

  py::class_<AgglomerationResult>(m, "AgglomerationResult")
      .def_readonly("complex", &AgglomerationResult::complex)
      .def_readonly("map", &AgglomerationResult::map);
  m.def("agglomerate", &agglomerate, py::arg("complex"), py::arg("labels"));

  py::class_<CoarsenResult>(m, "CoarsenResult")
      .def_readonly("complex", &CoarsenResult::complex)
      .def_readonly("map", &CoarsenResult::map);
  m.def(
      "coarsen",
      [](const Complex& c, const std::vector<Index>& labels, py::object protected_cells) {
        if (protected_cells.is_none()) return coarsen(c, labels);
        return coarsen(c, labels, protected_cells.cast<const Marker&>());
      },
      py::arg("complex"), py::arg("labels"), py::arg("protected_cells") = py::none());

  py::class_<CollapseResult>(m, "CollapseResult")
      .def_readonly("complex", &CollapseResult::complex)
      .def_readonly("map", &CollapseResult::map)
      .def_readonly("removed_pairs", &CollapseResult::removed_pairs);
  m.def("collapse", &collapse, py::arg("complex"));
  m.def("free_faces", &free_faces, py::arg("complex"));

  py::class_<SubdivisionResult>(m, "SubdivisionResult")
      .def_readonly("complex", &SubdivisionResult::complex)
      .def_property_readonly("vertex_offset",
                             [](const SubdivisionResult& r) { return copy_array(r.vertex_offset); })
      .def(
          "carrier_dim",
          [](const SubdivisionResult& r, std::size_t k) {
            check_dim(k, r.carrier_dim.size(), "carrier_dim");
            return copy_array(std::vector<Index>(r.carrier_dim[k].begin(), r.carrier_dim[k].end()));
          },
          py::arg("k"))
      .def(
          "carrier_index",
          [](const SubdivisionResult& r, std::size_t k) {
            check_dim(k, r.carrier_index.size(), "carrier_index");
            return copy_array(r.carrier_index[k]);
          },
          py::arg("k"))
      .def(
          "carrier_sign",
          [](const SubdivisionResult& r, std::size_t k) {
            check_dim(k, r.carrier_sign.size(), "carrier_sign");
            return copy_array(r.carrier_sign[k]);
          },
          py::arg("k"));
  m.def("barycentric_subdivision", &barycentric_subdivision, py::arg("complex"));

  // ---- queries ----------------------------------------------------------
  py::class_<FacetClassification>(m, "FacetClassification")
      .def_readonly("maximal", &FacetClassification::maximal)
      .def_readonly("boundary", &FacetClassification::boundary)
      .def_readonly("interior", &FacetClassification::interior)
      .def_readonly("nonmanifold", &FacetClassification::nonmanifold);
  m.def("classify_facets", &classify_facets, py::arg("complex"));
  m.def("is_closed", &is_closed, py::arg("complex"));

  py::class_<ManifoldReport>(m, "ManifoldReport")
      .def_readonly("manifold_like", &ManifoldReport::manifold_like)
      .def_readonly("pure", &ManifoldReport::pure)
      .def_readonly("facet_condition", &ManifoldReport::facet_condition)
      .def_readonly("links_connected", &ManifoldReport::links_connected)
      .def_readonly("offending", &ManifoldReport::offending);
  m.def("check_manifold", &check_manifold, py::arg("complex"));

  py::class_<LinkClassification>(m, "LinkClassification")
      .def_readonly("sphere", &LinkClassification::sphere)
      .def_readonly("ball", &LinkClassification::ball)
      .def_readonly("other", &LinkClassification::other);
  m.def("classify_vertex_links", &classify_vertex_links, py::arg("complex"));

  m.def("betti_numbers_z2", py::overload_cast<const Complex&>(&betti_numbers_z2),
        py::arg("complex"));
  m.def("betti_numbers_z2", py::overload_cast<const Complex&, const Marker&>(&betti_numbers_z2),
        py::arg("complex"), py::arg("relative_to"));

  py::class_<ComponentLabels>(m, "ComponentLabels")
      .def_property_readonly("label", [](const ComponentLabels& r) { return copy_array(r.label); })
      .def_readonly("count", &ComponentLabels::count);
  py::class_<ComplexComponents>(m, "ComplexComponents")
      .def(
          "label",
          [](const ComplexComponents& r, std::size_t k) {
            check_dim(k, r.label.size(), "label");
            return copy_array(r.label[k]);
          },
          py::arg("k"))
      .def_readonly("count", &ComplexComponents::count);
  m.def("connected_components", py::overload_cast<const Complex&>(&connected_components),
        py::arg("complex"));
  m.def("connected_components", py::overload_cast<const Complex&, int, int>(&connected_components),
        py::arg("complex"), py::arg("k"), py::arg("via"));
  m.def("connected_components",
        py::overload_cast<const Complex&, int, int, const Marker&>(&connected_components),
        py::arg("complex"), py::arg("k"), py::arg("via"), py::arg("exclude_via"));

  m.def("star_of", &star_of, py::arg("complex"), py::arg("cells"));
  m.def("closure_of", &closure_of, py::arg("complex"), py::arg("cells"));
  m.def("link_of", &link_of, py::arg("complex"), py::arg("cells"));
  m.def("frontier_of", &frontier_of, py::arg("complex"), py::arg("cells"));

  py::class_<CommonBoundary>(m, "CommonBoundary")
      .def_readonly("facets", &CommonBoundary::facets)
      .def_readonly("n_facets", &CommonBoundary::n_facets)
      .def_readonly("components", &CommonBoundary::components)
      .def_readonly("betti", &CommonBoundary::betti)
      .def_readonly("acyclic", &CommonBoundary::acyclic);
  m.def("common_boundary", &common_boundary, py::arg("complex"), py::arg("k"), py::arg("a"),
        py::arg("b"));
  m.def("excess_intersection", &excess_intersection, py::arg("complex"), py::arg("k"), py::arg("a"),
        py::arg("b"));
  m.def("amalgamates_to_cell", &amalgamates_to_cell, py::arg("complex"), py::arg("k"), py::arg("a"),
        py::arg("b"));
}
