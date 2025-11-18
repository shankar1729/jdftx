__all__ = ["JDFTxWrapper", "initialize", "finalize",  "ase"]

# Dispatch between CPU and GPU versions:
import os
if os.environ.get("CUDA_VISIBLE_DEVICES", ""):
    from ._pyjdftx_gpu import JDFTxWrapper, initialize, finalize
else:
    from ._pyjdftx import JDFTxWrapper, initialize, finalize

from . import ase
