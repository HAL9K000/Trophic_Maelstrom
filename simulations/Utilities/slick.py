"""
CUDA GPU Acceleration Module with Automatic CPU Fallback
Provides drop-in replacements for numpy/scipy with automatic GPU acceleration where possible
"""

# TO-DO: HANDLE PANDAS DF.ILOC[] CONVERSIONS MORE GRACEFULLY WITHOUT USING .GET() ON GPU ARRAYS

import os
import time
import sys
import warnings
import threading
from functools import wraps
import importlib

print(f"GPU_GLOW_UP: PID={os.getpid()}, __name__={__name__}, importing...")
import traceback
import inspect
#print(f"STACK TRACE for PID {os.getpid()}:")
#traceback.print_stack()
#print("=" * 50)

# Determine GPU usage from environment variable
USE_CUDA = os.getenv("USE_GPU", "0") == "1"
GPU_AVAILABLE = False
RAPIDS_AVAILABLE = False # Check if RAPIDs libraries are available
DASK_AVAILABLE = False # Check if Dask distributed is available.

# Standard CPU imports (always available)
import numpy as _cpu_np
import numpy.random as _cpu_nprandom
import scipy as _cpu_scipy
import scipy.fft as _cpu_fft
import scipy.signal as _cpu_signal
import scipy.interpolate as _cpu_interpolate
import scipy.stats as _cpu_stats
import scipy.ndimage as _cpu_ndimage
import pandas as _cpu_pandas

# Ensure warnings don't show full paths
warnings.formatwarning = lambda message, category, filename, lineno, line=None: f"{os.path.basename(filename)}:{lineno}: {category.__name__}: {message}\n"
#def short_formatwarning(message, category, filename, lineno, line=None):
#    return f"{os.path.basename(filename):{lineno}: {category.__name__}: {message}\n"

# Try GPU imports
if USE_CUDA:
    try:
        if sys.platform.startswith("darwin"):
            raise ImportError("CUDA libraries are not supported on macOS, use Linux machine for GPU support. 👺👺"); time.sleep(10)

        import cupy as _cupy
        import cupy.random as _cupy_random
        import cupyx.scipy.fft as _gpu_fft
        import cupyx.scipy.signal as _gpu_signal
        import cupyx.scipy.ndimage as _gpu_ndimage

        print("CuPY online... GPU acceleration enabled...✅"); time.sleep(1)

        GPU_AVAILABLE = True
        # RAPIDs libraries (WORK ONLY ON LINUX WITH NVIDIA GPUS)
        if sys.platform.startswith("linux"):
            
            try:
                #import cudf.pandas
                #cudf.pandas.install()
                #import cuml.accel
                #cuml.accel.install()
                #RAPIDS_AVAILABLE = True
                print("RAPIDs acceleration enabled... ✅✅"); time.sleep(1)
            except ImportError as e:
                warnings.warn(f"RAPIDs libraries not found ❌: {e} . Using CuPy only for GPU acceleration... 👹")
                time.sleep(1)
        elif sys.platform.startswith("win"):
            warnings.warn("RAPIDs libraries are not supported on Windows, using CuPy only for GPU support. 👺")
            time.sleep(1)
        
        #import pandas
        
    except ImportError as e:
        warnings.warn(f"Error importing GPU libraries: {e}. Falling back to CPU usage...❌❌❌"); time.sleep(5)
        USE_CUDA = False
        GPU_AVAILABLE = False

# Check if Dask distributed is available
try:
    import dask.distributed as _dask_distributed
    import dask as _dask
    from dask.distributed import Client as _dask_Client
    from dask import delayed as _dask_delayed
    from dask.distributed import progress as _dask_progress
    from dask.distributed import LocalCluster as _dask_LocalCluster
    from dask.diagnostics import ProgressBar as _dask_ProgressBar

    #DASK_AVAILABLE = True
    #print("Dask distributed is available and enabled for parallel processing.✅✅✅"); time.sleep(1)
except ImportError as e:
    warnings.warn(f"Dask distributed not found: {e}. Switching to joblib... True CPU-GPU interops not available. ❌"); time.sleep(2)

# Utlity functions and wrapper modules
"""Recursively convert GPU arrays to CPU arrays."""
def to_cpu(obj):
    #SOME STUFF
    if not GPU_AVAILABLE:
        return obj
    if isinstance(obj, ( _cpu_np.ndarray, str, int, float)):
        return obj  # Already safe
    if hasattr(_cupy, 'ndarray') and isinstance(obj, _cupy.ndarray):
        return obj.get()
    elif isinstance(obj, dict):
        return {k: to_cpu(v) for k, v in obj.items()}
    elif isinstance(obj, (list, tuple)):
        return type(obj)(to_cpu(v) for v in obj)
    elif isinstance(obj, set):
        return {to_cpu(v) for v in obj}
    return obj

"""Convert CPU arrays (iterators) to GPU arrays where possible."""
def to_gpu(obj):
    #SOME STUFF
    if not GPU_AVAILABLE:
        return obj
    if hasattr(_cupy, 'ndarray') and isinstance(obj, _cupy.ndarray):
        return obj
    if isinstance(obj, _cpu_np.ndarray):
        try:
            return _cupy.asarray(obj)
        except Exception:
            return obj
    elif isinstance(obj, dict):
        return {k: to_gpu(v) for k, v in obj.items()}
    elif isinstance(obj, (list, tuple)):
        return type(obj)(to_gpu(v) for v in obj)
    elif isinstance(obj, set):
        return {to_gpu(v) for v in obj}
    return obj

def cpu_fallback(func_name, cpu_module, gpu_module=None):
    """Decorator to provide CPU fallback for GPU functions."""
    def decorator(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            if GPU_AVAILABLE and gpu_module:
                try:
                    # Try GPU version first
                    gpu_args = to_gpu(args)
                    gpu_kwargs = to_gpu(kwargs)
                    result = func(*gpu_args, **gpu_kwargs)
                    return to_cpu(result)
                except Exception as e:
                    warnings.warn(f"GPU operation failed for {func_name}: {e}. Falling back to CPU.", stacklevel=2)
            
            # Use CPU version
            cpu_args = to_cpu(args)
            cpu_kwargs = to_cpu(kwargs)
            cpu_func = getattr(cpu_module, func_name)
            return cpu_func(*cpu_args, **cpu_kwargs)
        return wrapper
    return decorator



class SmartModule:
    """A module wrapper that tries GPU first, falls back to CPU."""
    def __init__(self, cpu_module, gpu_module=None, module_name=""):
        self._cpu_module = cpu_module
        self._gpu_module = gpu_module
        self._module_name = module_name

        # Define attributes that should always come from CPU module
        # (types, constants, etc.)
        self._cpu_only_attrs = {
            'ndarray', 'dtype', 'generic',
            'float64', 'float32', 'float16', 'int64', 'int32', 'int16', 'int8',
            'uint64', 'uint32', 'uint16', 'uint8', 'bool_', 'complex64', 'complex128',
            # Warnings and exceptions
            #Constants
            'nan', 'inf', 'pi', 'e', 'euler_gamma', 'newaxis', 'NINF', 'NZERO', 'PZERO' }
        

        # Explicit fallbacks for Numpy exceptions (as they are no longer exposed in the public namespace since Numpy 1.2)
        self._fallback_attrs = {
            'ComplexWarning': self.safe_get_exception(cpu_module, 'ComplexWarning', Warning),
            'VisibleDeprecationWarning': self.safe_get_exception(cpu_module, 'VisibleDeprecationWarning', Warning),
            'TooHardError': self.safe_get_exception(cpu_module, 'TooHardError', Exception),
            'AxisError': self.safe_get_exception(cpu_module, 'AxisError', ValueError),
            'RankWarning': self.safe_get_exception(cpu_module, 'RankWarning', Warning)
            # You can add other exception fallbacks here
        }
            # Add other constants and types as needed
        
        
    def __getattr__(self, name):
        # For types and constants, always use CPU module
        if name in self._cpu_only_attrs:
            return getattr(self._cpu_module, name)
        # Check if it's in fallback exception list
        if name in self._fallback_attrs:
            return self._fallback_attrs[name]
        # Check if attribute exists in GPU module first
        if GPU_AVAILABLE and self._gpu_module and hasattr(self._gpu_module, name):
            gpu_attr = getattr(self._gpu_module, name)
            if callable(gpu_attr):
                @wraps(gpu_attr)
                def smart_func(*args, **kwargs):
                    try:
                        # Convert args to GPU
                        gpu_args = to_gpu(args)
                        gpu_kwargs = to_gpu(kwargs)
                        result = gpu_attr(*gpu_args, **gpu_kwargs)
                                # Special handling for boolean indexing operations
                        if name == '__getitem__' and len(args) == 1:
                            # If we have a GPU boolean mask but CPU data, convert mask to CPU
                            mask = args[0]
                            if is_gpu_array(result) and not is_gpu_array(gpu_args[0]) and hasattr(mask, 'dtype') and mask.dtype == bool:
                                warnings.warn(f"Converting GPU boolean mask to CPU for indexing compatibility", stacklevel=2)
                                return gpu_args[0][to_cpu(mask)]
                            # Convert GPU scalars to CPU to avoid interop issues with numpy methods
                            #if is_gpu_array(result) and result.ndim == 0:  # scalar array
                            #    return to_cpu(result)
                        #if name == 'fftconvolve':
                        #    print(f"GPU fftconvolve succeeded, result type: {type(result)}")
                        # ✅ Keep result on GPU! Don't convert to CPU automatically    
                        return result  # or to_cpu(result) if you want to ensure always returning
                        #CPU arrays for compatibility (MASSIVE performance hit)
                    except Exception as e:
                        # Fall back to CPU
                        trace_str= ''.join(traceback.format_exception(type(e), e, tb=e.__traceback__))
                        caller_trace = inspect.stack()[1]; 
                        err_filename = os.path.basename(caller_trace.filename); lineo = caller_trace.lineno
                        warnings.warn(f"GPU operation failed ❌ for {self._module_name}.{name} \n in 📄: {err_filename}, L {lineo}"
                                      f" with Error: {e}", stacklevel=2)
                        # Enhanced fallback with type conversion
                        if "Implicit conversion" in str(e) or "not allowed" in str(e):
                            warnings.warn(f"GPU/CPU type mismatch in {self._module_name}.{name}, converting to compatible types", stacklevel=2)
                        if name == 'fftconvolve':
                            print(f"GPU fftconvolve FAILED: {e}", stacklevel=2)
                        cpu_attr = getattr(self._cpu_module, name)
                        cpu_args = to_cpu(args)
                        cpu_kwargs = to_cpu(kwargs)
                        return cpu_attr(*cpu_args, **cpu_kwargs)
                return smart_func
            else:
                return gpu_attr
        
        # Fall back to CPU module
        cpu_attr = getattr(self._cpu_module, name)
        if callable(cpu_attr):
            @wraps(cpu_attr)
            def cpu_func(*args, **kwargs):
                cpu_args = to_cpu(args)
                cpu_kwargs = to_cpu(kwargs)
                return cpu_attr(*cpu_args, **cpu_kwargs)
            return cpu_func
        return cpu_attr
    
     # Explicit fallbacks for Numpy exceptions (as they are no longer exposed in the public namespace since Numpy 1.2)
    def safe_get_exception(self, module, exc_name, fallback):
        try:
            return getattr(getattr(module, 'exceptions'), exc_name)
        except AttributeError:
            return fallback

class CPUOnlyModule:
    """Wrapper for modules that don't have GPU equivalents."""
    def __init__(self, cpu_module, module_name=""):
        self._cpu_module = cpu_module
        self._module_name = module_name
        
    def __getattr__(self, name):
        attr = getattr(self._cpu_module, name)
        if callable(attr):
            @wraps(attr)
            def cpu_only_func(*args, **kwargs):
                # Convert all inputs to CPU
                cpu_args = to_cpu(args)
                cpu_kwargs = to_cpu(kwargs)
                return attr(*cpu_args, **cpu_kwargs)
            return cpu_only_func
        return attr

""" IMPORTANT: Decorator to ensure function ALWAYS returns CPU arrays - USE SPARINGLY DUE TO PERFORMANCE HIT! """     
def ensure_cpu(func):
    @wraps(func)
    def wrapper(*args, **kwargs):
        result = func(*args, **kwargs)
        return to_cpu(result)
    return wrapper

def cpu_safe_import(module_path):
    """Import and wrap a CPU-only module."""
    try:
        mod = importlib.import_module(module_path)
        return CPUOnlyModule(mod, module_path)
    except ImportError as e:
        raise ImportError(f"Could not import {module_path}: {e}")

# Create smart modules
if GPU_AVAILABLE:
    numpy = SmartModule(_cpu_np, _cupy, "numpy")
    fft = SmartModule(_cpu_fft, _gpu_fft, "scipy.fft")
    signal = SmartModule(_cpu_signal, _gpu_signal, "scipy.signal")
    ndimage = SmartModule(_cpu_ndimage, _gpu_ndimage, "scipy.ndimage")
    # Pre-wrapped common CPU-only modules for convenience
    interpolate = cpu_safe_import("scipy.interpolate")
    stats = cpu_safe_import("scipy.stats")
else:
    numpy = SmartModule(_cpu_np, None, "numpy")
    fft = SmartModule(_cpu_fft, None, "scipy.fft")
    signal = SmartModule(_cpu_signal, None, "scipy.signal")
    ndimage = SmartModule(_cpu_ndimage, None, "scipy.ndimage")
    interpolate = _cpu_scipy.interpolate
    stats = _cpu_scipy.stats

# Create a numpy alias for convenience
np = numpy

'''# Function to set up CUDA context for DASK workers (for true GPU interop)'''
def setup_daskworker_gpu_context():
    """Initialize CUDA context for this worker thread"""
    if GPU_AVAILABLE:
        _cupy.cuda.Device(0).use()  # or assign different devices if multi-GPU
    # Store context info in thread-local storage if needed

def init_daskworker_cuda_context():
    if not GPU_AVAILABLE:
        warnings.warn("GPU acceleration is not available. Dask workers will run on CPU only.")
        return
    
    import threading
    
    device = _cupy.cuda.Device(0)
    device.use()
    context = _cupy.cuda.runtime.getCurrentContext()
    print(f"Worker {threading.current_thread().name}: Context ID {context}")
    # Should see different context IDs for each worker
    
# OTHER USEFUL UTILITIES FOR EXPLICIT TO_CPU() CONVERSIONS BY USER!
def asnumpy(arr):
    """Explicitly convert array to CPU numpy array - use when you need CPU data"""
    return to_cpu(arr)

def asarray(arr):
    """Convert to appropriate array type (GPU if available, CPU otherwise)"""
    if GPU_AVAILABLE:
        return to_gpu(arr)
    return _cpu_np.asarray(arr)

def is_gpu_array(arr):
    """Check if array is on GPU"""
    if not GPU_AVAILABLE:
        return False
    return hasattr(_cupy, 'ndarray') and isinstance(arr, _cupy.ndarray)



# Expose key functions at module level
__all__ = ['np', 'numpy', 'fft', 'signal', 'interpolate', 'stats', 'ndimage', 'init_daskworker_cuda_context',
           'to_cpu', 'to_gpu', 'asnumpy', 'asarray', 'is_gpu_array', 'setup_daskworker_gpu_context',
           'cpu_safe_import', 'ensure_cpu', 'GPU_AVAILABLE', 'USE_GPU', 'RAPIDS_AVAILABLE', 'DASK_AVAILABLE']
