from slick import to_cpu, to_gpu
import importlib
import numpy as _np

''' This script provides a CPU fallback for modules that are typically used with GPU acceleration.
It allows for seamless conversion of arguments to CPU-compatible formats when calling functions 
from packages that may not support GPU operations, such as 'esda.moran' or 'scipy.interpolate'. This is useful for ensuring
that code can run on systems without GPU support or when GPU libraries are not available. '''

class CPUFallbackModule:
    def __init__(self, module):
        self._module = module

    def __getattr__(self, name):
        attr = getattr(self._module, name)
        if callable(attr):
            def wrapped(*args, **kwargs):
                return attr(*to_cpu(args), **to_cpu(kwargs))
            return wrapped
        return attr
    
def cpu_safe_import(module_path):
    """Import and wrap a CPU-only module (e.g., 'esda.moran', 'scipy.stats' etc)."""
    mod = importlib.import_module(module_path)
    return CPUFallbackModule(mod)