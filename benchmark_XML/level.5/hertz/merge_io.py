#!/usr/bin/env python3
"""Merge two Exodus II files (one per Tahoe element group) into one.

When a Tahoe XML has two element groups (e.g. Tet + Hex with different
<tetrahedron/> / <hexahedron/> geometry tags), each group writes its
own *.exo file (io0, io1, …).  ParaView can load them together but the
ergonomics are awkward.  This script merges them into a single Exodus
file for one-click visualisation.

Usage:
    python3 merge_io.py hertz_tet_hex_implicit.io0.exo \\
                       hertz_tet_hex_implicit.io1.exo \\
                       hertz_tet_hex_implicit_merged.exo

Strategy: append node and element data from io1 into a copy of io0,
remapping io1's node numbers by an offset so the two meshes don't
collide.  Reuses the same nodal-variable names if they appear in
both files; pads with zeros where one file is missing a variable
the other has.

Limitations:
  - Both files must have the same time-step list (true by
    construction when both are written by the same Tahoe run).
  - Element-block IDs are renumbered sequentially across both files.
  - Side sets and node sets are kept from io0 only — the merged file
    isn't intended for re-running, only for visualisation.
"""
import os
import sys
import shutil

import numpy as np
from netCDF4 import Dataset


def merge(io0_path: str, io1_path: str, out_path: str):
    print(f"reading {io0_path}")
    print(f"reading {io1_path}")
    print(f"writing {out_path}")

    a = Dataset(io0_path, "r")
    b = Dataset(io1_path, "r")

    n_nod_a = a.dimensions["num_nodes"].size
    n_nod_b = b.dimensions["num_nodes"].size
    n_el_a  = a.dimensions["num_elem"].size
    n_el_b  = b.dimensions["num_elem"].size
    nblk_a  = a.dimensions["num_el_blk"].size
    nblk_b  = b.dimensions["num_el_blk"].size

    n_nod = n_nod_a + n_nod_b
    n_el  = n_el_a  + n_el_b
    nblk  = nblk_a  + nblk_b

    # Use whichever variable list is the longer (most files are equal anyway).
    nv_a = a.dimensions.get("num_nod_var").size if "num_nod_var" in a.dimensions else 0
    nv_b = b.dimensions.get("num_nod_var").size if "num_nod_var" in b.dimensions else 0
    nv   = max(nv_a, nv_b)

    nt_a = a.dimensions["time_step"].size
    nt_b = b.dimensions["time_step"].size
    if nt_a != nt_b:
        print(f"  warning: time_step dim mismatch — a={nt_a}, b={nt_b}, "
              f"using min")
    nt = min(nt_a, nt_b)

    # ---- create output ----
    o = Dataset(out_path, "w", format="NETCDF3_64BIT_OFFSET")
    o.title = (a.title if hasattr(a, "title") else "") + " + " \
            + (b.title if hasattr(b, "title") else "")
    o.api_version = a.api_version if hasattr(a, "api_version") else 4.98
    o.version     = a.version     if hasattr(a, "version")     else 4.98
    o.floating_point_word_size = 8
    o.file_size = 0

    # dimensions
    o.createDimension("len_string", a.dimensions["len_string"].size)
    o.createDimension("len_line",   a.dimensions["len_line"].size)
    o.createDimension("len_name",   a.dimensions.get("len_name",
                                       a.dimensions["len_string"]).size)
    o.createDimension("four", 4)
    o.createDimension("num_dim",   a.dimensions["num_dim"].size)
    o.createDimension("num_nodes", n_nod)
    o.createDimension("num_elem",  n_el)
    o.createDimension("num_el_blk", nblk)
    o.createDimension("time_step", None)
    o.createDimension("num_qa_rec", 1)
    if nv > 0:
        o.createDimension("num_nod_var", nv)

    # time
    t_in  = a.variables["time_whole"][:nt]
    o.createVariable("time_whole", "f8", ("time_step",))[:nt] = t_in

    # coordinates: stack a then b
    cx_a = a.variables["coordx"][:]; cy_a = a.variables["coordy"][:]
    cz_a = a.variables["coordz"][:] if "coordz" in a.variables else None
    cx_b = b.variables["coordx"][:]; cy_b = b.variables["coordy"][:]
    cz_b = b.variables["coordz"][:] if "coordz" in b.variables else None
    cx = np.concatenate([cx_a, cx_b]); cy = np.concatenate([cy_a, cy_b])
    o.createVariable("coordx", "f8", ("num_nodes",))[:] = cx
    o.createVariable("coordy", "f8", ("num_nodes",))[:] = cy
    if cz_a is not None and cz_b is not None:
        o.createVariable("coordz", "f8", ("num_nodes",))[:] = \
            np.concatenate([cz_a, cz_b])

    # qa records (placeholder — minimal)
    nm = o.createVariable("qa_records", "S1",
                          ("num_qa_rec", "four", "len_string"))
    nm[:] = b" "

    # node_num_map: identity (1..n_nod)
    nm = o.createVariable("node_num_map", "i4", ("num_nodes",))
    nm[:] = np.arange(1, n_nod + 1, dtype="i4")

    # element blocks: blocks from a first, then from b (with element offset)
    eb_status = np.ones(nblk, dtype="i4")
    eb_prop1  = np.arange(1, nblk + 1, dtype="i4")
    o.createVariable("eb_status", "i4", ("num_el_blk",))[:] = eb_status
    o.createVariable("eb_prop1",  "i4", ("num_el_blk",))[:] = eb_prop1

    out_blk = 0
    for src, src_label, node_offset in [(a, "a", 0), (b, "b", n_nod_a)]:
        sblk = src.dimensions["num_el_blk"].size
        for ib in range(sblk):
            ne_dim_name  = f"num_el_in_blk{ib+1}"
            nen_dim_name = f"num_nod_per_el{ib+1}"
            ne  = src.dimensions[ne_dim_name].size
            nen = src.dimensions[nen_dim_name].size
            out_ne_dim  = f"num_el_in_blk{out_blk+1}"
            out_nen_dim = f"num_nod_per_el{out_blk+1}"
            o.createDimension(out_ne_dim,  ne)
            o.createDimension(out_nen_dim, nen)
            connect_in  = src.variables[f"connect{ib+1}"][:]
            connect_out = o.createVariable(
                f"connect{out_blk+1}", "i4",
                (out_ne_dim, out_nen_dim))
            connect_out[:] = connect_in + node_offset
            connect_out.elem_type = src.variables[f"connect{ib+1}"].elem_type
            out_blk += 1

    # nodal variables: union of names across the two files
    name_to_idx_a = {}
    name_to_idx_b = {}
    if nv_a > 0:
        for i in range(nv_a):
            n = b"".join(a.variables["name_nod_var"][i]).decode().strip("\x00 ")
            name_to_idx_a[n] = i
    if nv_b > 0:
        for i in range(nv_b):
            n = b"".join(b.variables["name_nod_var"][i]).decode().strip("\x00 ")
            name_to_idx_b[n] = i
    union_names = sorted(set(name_to_idx_a) | set(name_to_idx_b))
    if union_names:
        o.createDimension if "num_nod_var" not in o.dimensions else None
        if "num_nod_var" not in o.dimensions:
            o.createDimension("num_nod_var", len(union_names))
        nv_out = len(union_names)
        nm = o.createVariable("name_nod_var", "S1",
                              ("num_nod_var", "len_string"))
        for i, name in enumerate(union_names):
            for j, ch in enumerate(name):
                nm[i, j] = ch.encode()
        # values: per-variable per-time
        for vi, name in enumerate(union_names):
            arr = np.zeros((nt, n_nod), dtype="f8")
            if name in name_to_idx_a:
                arr[:, :n_nod_a] = a.variables[
                    f"vals_nod_var{name_to_idx_a[name]+1}"][:nt]
            if name in name_to_idx_b:
                arr[:, n_nod_a:] = b.variables[
                    f"vals_nod_var{name_to_idx_b[name]+1}"][:nt]
            v = o.createVariable(f"vals_nod_var{vi+1}", "f8",
                                 ("time_step", "num_nodes"))
            v[:nt] = arr

    a.close(); b.close(); o.close()
    print(f"  merged: {n_nod} nodes, {n_el} elements ({nblk} blocks), "
          f"{nt} frames, {len(union_names) if union_names else 0} nodal vars")


def main():
    if len(sys.argv) != 4:
        print(__doc__)
        sys.exit(1)
    merge(sys.argv[1], sys.argv[2], sys.argv[3])


if __name__ == "__main__":
    main()
