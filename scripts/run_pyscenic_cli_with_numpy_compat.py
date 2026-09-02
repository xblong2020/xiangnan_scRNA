from __future__ import annotations

import runpy
import sys
import os
from pathlib import Path

import numpy as np


def patch_pyscenic_prune_from_delayed() -> None:
    import dask.dataframe as dd
    import pyscenic.prune as prune

    original = dd.from_delayed

    def from_delayed_compat(dfs, *args, **kwargs):
        if not hasattr(dfs, "__len__"):
            dfs = list(dfs)
        return original(dfs, *args, **kwargs)

    dd.from_delayed = from_delayed_compat
    prune.from_delayed = from_delayed_compat


def patch_arboreto_create_graph() -> None:
    import arboreto.algo as algo
    import arboreto.core as core
    from dask import delayed
    from dask.dataframe import from_delayed

    def create_graph_compat(
        expression_matrix,
        gene_names,
        tf_names,
        regressor_type,
        regressor_kwargs,
        client,
        target_genes="all",
        limit=None,
        include_meta=False,
        early_stop_window_length=core.EARLY_STOP_WINDOW_LENGTH,
        repartition_multiplier=1,
        seed=core.DEMON_SEED,
    ):
        assert expression_matrix.shape[1] == len(gene_names)
        assert client, "client is required"

        tf_matrix, tf_matrix_gene_names = core.to_tf_matrix(expression_matrix, gene_names, tf_names)
        future_tf_matrix = client.scatter(tf_matrix, broadcast=True)
        [future_tf_matrix_gene_names] = client.scatter([tf_matrix_gene_names], broadcast=True)

        delayed_link_dfs = []
        delayed_meta_dfs = []
        for target_gene_index in core.target_gene_indices(gene_names, target_genes):
            target_gene_name = delayed(gene_names[target_gene_index], pure=True)
            target_gene_expression = delayed(expression_matrix[:, target_gene_index], pure=True)
            if include_meta:
                delayed_link_df, delayed_meta_df = delayed(core.infer_partial_network, pure=True, nout=2)(
                    regressor_type,
                    regressor_kwargs,
                    future_tf_matrix,
                    future_tf_matrix_gene_names,
                    target_gene_name,
                    target_gene_expression,
                    include_meta,
                    early_stop_window_length,
                    seed,
                )
                delayed_link_dfs.append(delayed_link_df)
                delayed_meta_dfs.append(delayed_meta_df)
            else:
                delayed_link_df = delayed(core.infer_partial_network, pure=True)(
                    regressor_type,
                    regressor_kwargs,
                    future_tf_matrix,
                    future_tf_matrix_gene_names,
                    target_gene_name,
                    target_gene_expression,
                    include_meta,
                    early_stop_window_length,
                    seed,
                )
                delayed_link_dfs.append(delayed_link_df)

        if not delayed_link_dfs:
            raise ValueError("No target genes available for arboreto GRN inference.")

        all_links_df = from_delayed(delayed_link_dfs, meta=core._GRN_SCHEMA)
        maybe_limited_links_df = all_links_df.nlargest(limit, columns=["importance"]) if limit else all_links_df
        n_parts = len(client.ncores()) * repartition_multiplier

        if include_meta:
            if not delayed_meta_dfs:
                raise ValueError("No target genes available for arboreto GRN metadata inference.")
            all_meta_df = from_delayed(delayed_meta_dfs, meta=core._META_SCHEMA)
            return maybe_limited_links_df.repartition(npartitions=n_parts), all_meta_df.repartition(npartitions=n_parts)
        return maybe_limited_links_df.repartition(npartitions=n_parts)

    core.create_graph = create_graph_compat
    algo.create_graph = create_graph_compat


def main() -> None:
    # pySCENIC 0.12.x still imports removed NumPy aliases during CLI startup.
    # Keep the compatibility shim local to this process.
    for name, value in {
        "object": object,
        "bool": bool,
        "int": int,
        "float": float,
    }.items():
        if name not in np.__dict__:
            setattr(np, name, value)

    compat_path = str(Path(__file__).resolve().parent / "pyscenic_compat_site")
    existing_pythonpath = os.environ.get("PYTHONPATH", "")
    pythonpath_parts = [compat_path] + ([existing_pythonpath] if existing_pythonpath else [])
    os.environ["PYTHONPATH"] = os.pathsep.join(pythonpath_parts)

    patch_arboreto_create_graph()
    patch_pyscenic_prune_from_delayed()
    sys.argv = ["pyscenic", *sys.argv[1:]]
    runpy.run_module("pyscenic.cli.pyscenic", run_name="__main__")


if __name__ == "__main__":
    main()
