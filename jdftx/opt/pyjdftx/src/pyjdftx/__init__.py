from __future__ import annotations

__all__ = ["JDFTxWrapper", "initialize", "finalize",  "ase"]

# Dispatch between CPU and GPU versions:
import os
if os.environ.get("CUDA_VISIBLE_DEVICES", ""):
    from ._pyjdftx_gpu import JDFTxWrapper, initialize, finalize
else:
    from ._pyjdftx import JDFTxWrapper, initialize, finalize

from . import ase


def comm_divide_evenly(
    comm: mpi4py.MPI.Comm, n_groups: int
) -> tuple[mpi4py.MPI.Comm, int]:
    """Split a communicator evenly into a specified number of groups,
    returning the communicator for the group the present process
    belongs to, and the 0-based index of that group.
    (comm.size should be divisible by n_groups.)"""
    n_procs = comm.size
    i_proc = comm.rank    
    n_each = n_procs // n_groups
    if n_each * n_groups != n_procs:
        raise ValueError(
            f"Number of MPI processes ({n_procs}) not divisible by {n_groups = }"
        )
    i_group = i_proc // n_each
    group_rank = i_proc % n_each
    return comm.Split(i_group, group_rank), i_group
