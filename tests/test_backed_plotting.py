"""Run with `python -m unittest discover -s tests` (no Latch runtime needed)."""
import ast
from pathlib import Path
import shutil
import tempfile
import types
import unittest

import anndata
import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp


ROOT = Path(__file__).resolve().parents[1]


def load_functions(path, namespace):
    tree = ast.parse((ROOT / path).read_text())
    names = {"feature_column_indices", "read_matrix_columns", "read_backed_h5ad",
             "_gene_list_zscore_heatmap", "compute_cluster_marker_heatmap_from_degs",
             "choose_heatmap_layer", "cluster_marker_sort_key", "_first_existing_column",
             "cluster_marker_zscore_heatmap", "plotting_candidates"}
    functions = [node for node in ast.walk(tree)
                 if isinstance(node, ast.FunctionDef) and node.name in names]
    exec(compile(ast.Module(body=functions, type_ignores=[]), str(path), "exec"), namespace)


class BackedPlottingTests(unittest.TestCase):
    def setUp(self):
        self.ns = dict(np=np, pd=pd, sp=sp, anndata=anndata, AnnData=anndata.AnnData,
                       Path=Path, List=list, sc=types.SimpleNamespace(read_h5ad=anndata.read_h5ad))
        load_functions("welcome/init.py", self.ns)
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.path = Path(self.temp.name) / "data.h5ad"
        self.values = np.arange(20, dtype=float).reshape(4, 5)
        self.adata = anndata.AnnData(self.values,
                                    obs=pd.DataFrame({"cluster": ["a", "a", "b", "b"]},
                                                     index=["s0", "s1", "s2", "s3"]),
                                    var=pd.DataFrame(index=["g0", "g1", "g2", "g3", "g4"]))

    def write_backed(self, matrix):
        self.adata.X = matrix
        self.adata.write(self.path)
        backed = anndata.read_h5ad(self.path, backed="r+")
        self.addCleanup(backed.file.close)
        return backed

    def test_requested_columns_dense_and_sparse(self):
        for matrix in (self.values, sp.csr_matrix(self.values), sp.csc_matrix(self.values)):
            with self.subTest(matrix=type(matrix).__name__):
                backed = self.write_backed(matrix)
                result = self.ns["read_matrix_columns"](backed.X, [4, 1, 4, 0])
                if sp.issparse(result):
                    result = result.toarray()
                np.testing.assert_array_equal(result, self.values[:, [4, 1, 4, 0]])
                backed.file.close()

    def test_custom_heatmaps_with_duplicate_names(self):
        self.adata.var_names = ["g0", "g1", "g1", "g3", "g4"]
        backed = self.write_backed(self.values)
        for cell in ("welcome/rna_top_heatmap.py", "ge_h5_viewer/atac_top_heatmap.py"):
            load_functions(cell, self.ns)
            fn = self.ns["_gene_list_zscore_heatmap"]
            actual, missing = fn(backed, ["g4", "g1", "absent"], "cluster", "user-selected", "b,a")
            expected, _ = fn(self.adata, ["g4", "g1"], "cluster", "user-selected", "b,a")
            pd.testing.assert_frame_equal(actual, expected)
            self.assertEqual(missing, ["absent"])
            self.assertEqual(actual.index.tolist(), ["b", "a"])

    def test_marker_heatmap_matches_memory(self):
        backed = self.write_backed(self.values)
        degs = pd.DataFrame({"cluster": ["a", "a", "b"],
                             "names": ["g4", "g1", "g0"],
                             "logfoldchanges": [1., 2., 3.], "pvals": [.001] * 3})
        fn = self.ns["compute_cluster_marker_heatmap_from_degs"]
        args = (degs, 2, "pvals", .05, .25, "user-selected", "b,a")
        pd.testing.assert_frame_equal(fn(backed, *args), fn(self.adata, *args))

    def test_optimized_file_preference_and_legacy_fallback(self):
        load_functions("select_data/select_data.py", self.ns)
        fn = self.ns["plotting_candidates"]
        for stem in ("rna_copro", "atac_gs_copro"):
            names = [stem, f"{stem}_sm", f"{stem}_sm_ge"]
            self.ns["children"] = [types.SimpleNamespace(name=lambda n=n: n + ".h5ad")
                                    for n in names]
            self.assertEqual([name for name, _ in fn(stem)], list(reversed(names)))
            self.ns["children"] = self.ns["children"][:2]
            self.assertEqual(fn(stem)[0][0], f"{stem}_sm")
            self.ns["children"] = []
            with self.assertRaises(ValueError):
                fn(stem)

    def test_corrupt_download_retry_and_save_preserves_chunks(self):
        self.adata.write(self.path)
        with h5py.File(self.path, "r+") as handle:
            attrs = dict(handle["X"].attrs)
            del handle["X"]
            handle.create_dataset("X", data=self.values, chunks=(4, 1)).attrs.update(attrs)
        source = self.path
        calls = []

        class Remote:
            def name(self):
                return "download.h5ad"

            def download(self, destination, cache):
                calls.append(cache)
                if cache:
                    destination.write_bytes(b"corrupt cached download")
                else:
                    shutil.copyfile(source, destination)

        # Resolve the notebook's relative download path in our temporary directory.
        self.ns["Path"] = lambda name: source.parent / name
        backed = self.ns["read_backed_h5ad"](Remote())
        self.addCleanup(backed.file.close)
        self.assertEqual(calls, [True, False])
        self.assertTrue(backed.isbacked)
        backed.obs["saved"] = "yes"
        backed.write(backed.filename)
        with h5py.File(backed.filename, "r") as handle:
            self.assertEqual(handle["X"].chunks, (4, 1))
        reread = anndata.read_h5ad(backed.filename)
        self.assertTrue((reread.obs["saved"] == "yes").all())
        np.testing.assert_array_equal(reread.X, self.values)


if __name__ == "__main__":
    unittest.main()
