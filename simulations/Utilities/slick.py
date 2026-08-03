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
import itertools
import zlib
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
import numpy.linalg as _cpu_linalg
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
        import cupy.linalg as _cupy_linalg
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

    DASK_AVAILABLE = True
    print("Dask distributed is available and enabled for parallel processing.✅✅✅"); time.sleep(1)
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
    linalg = SmartModule(_cpu_linalg, _cupy_linalg, "numpy.linalg")
    fft = SmartModule(_cpu_fft, _gpu_fft, "scipy.fft")
    signal = SmartModule(_cpu_signal, _gpu_signal, "scipy.signal")
    ndimage = SmartModule(_cpu_ndimage, _gpu_ndimage, "scipy.ndimage")
    # Pre-wrapped common CPU-only modules for convenience
    interpolate = cpu_safe_import("scipy.interpolate")
    stats = cpu_safe_import("scipy.stats")
else:
    numpy = SmartModule(_cpu_np, None, "numpy")
    linalg = SmartModule(_cpu_linalg, None, "numpy.linalg")
    fft = SmartModule(_cpu_fft, None, "scipy.fft")
    signal = SmartModule(_cpu_signal, None, "scipy.signal")
    ndimage = SmartModule(_cpu_ndimage, None, "scipy.ndimage")
    interpolate = _cpu_scipy.interpolate
    stats = _cpu_scipy.stats

# Create a numpy alias for convenience
np = numpy
la = linalg

# =============================================================================
# MULTI-GPU DASK WORKER CONTEXTS 🎛️
# =============================================================================
# Pinning every worker to Device(0) wastes every card but the first. Slick now
# spreads workers across the visible devices, with the policy under your control.

_LOCAL_RANK_COUNTER = itertools.count()   # per-process, covers threaded workers
_WORKER_DEVICE = threading.local()        # remembers this worker's device


def visible_devices():
    """List the CUDA device ids Slick is allowed to use, honouring CUDA_VISIBLE_DEVICES.

    Override with SLICK_GPU_DEVICES="0,2,3" to carve out a subset.
    """
    if not GPU_AVAILABLE:
        return []
    override = os.getenv("SLICK_GPU_DEVICES", "").strip()
    if override:
        try:
            return [int(d) for d in override.replace(" ", "").split(",") if d != ""]
        except ValueError:
            warnings.warn(f"Could not parse SLICK_GPU_DEVICES={override!r} ⚠️ - using all devices.")
    try:
        return list(range(_cupy.cuda.runtime.getDeviceCount()))
    except Exception as e:
        warnings.warn(f"Could not enumerate CUDA devices: {e} ⚠️ - assuming a single device.")
        return [0]


def worker_rank():
    """Best-effort stable integer rank for the calling worker/thread.

    Dask workers are separate *processes*, so a naive module-level counter gives
    every worker rank 0 and round-robin collapses onto one card. We therefore
    look for real identity first, and only then fall back to a local counter.
    Cascade: SLICK_WORKER_RANK env -> Dask worker name -> CRC32 of the worker
    address (stable across runs, unlike PYTHONHASHSEED-salted hash()) -> a
    per-process counter for threads inside one worker.
    """
    env_rank = os.getenv("SLICK_WORKER_RANK")
    if env_rank is not None:
        try:
            return int(env_rank), "SLICK_WORKER_RANK"
        except ValueError:
            pass
    if DASK_AVAILABLE:
        try:
            worker = _dask_distributed.get_worker()
            name = getattr(worker, "name", None)
            if isinstance(name, int):
                return name, "dask worker name"
            if isinstance(name, str) and name.strip().lstrip("-").isdigit():
                return int(name), "dask worker name"
            ident = str(name or getattr(worker, "address", "") or worker.id)
            if ident:
                return zlib.crc32(ident.encode()), f"crc32({ident})"
        except (ValueError, AttributeError, ImportError):
            pass  # Not inside a Dask worker - fine, keep going.
        except Exception:
            pass
    return next(_LOCAL_RANK_COUNTER), "process-local counter"


def _least_busy_device(devices):
    """Pick the visible device with the most free memory right now."""
    best, best_free = devices[0], -1
    for dev in devices:
        try:
            with _cupy.cuda.Device(dev):
                free, _total = _cupy.cuda.runtime.memGetInfo()
        except Exception:
            continue
        if free > best_free:
            best, best_free = dev, free
    return best


def select_cuda_device(policy="round-robin", device=None, devices=None):
    """Resolve which CUDA device this worker should use. Pure - it selects, it does not bind.

    Parameters
    ----------
    policy : {"round-robin", "explicit", "least-busy", "single"} or callable
        - "round-robin" (default): rank % len(devices), spreading workers evenly.
        - "explicit": use `device` (int, callable, or {worker_name/rank: device} mapping).
        - "least-busy": whichever visible device currently has the most free memory.
        - "single": legacy behaviour, everyone on `device` or device 0.
        - callable: your own `policy(rank, devices) -> device_id`.
    device : int | callable | dict, optional
        User-defined selection. A callable receives `(rank, devices)`.
    devices : sequence of int, optional
        Restrict the pool; defaults to :func:`visible_devices`.
    """
    pool = list(devices) if devices is not None else visible_devices()
    if not pool:
        return None
    rank, source = worker_rank()

    if callable(policy):
        chosen = policy(rank, pool)
    elif callable(device):
        chosen = device(rank, pool)
    elif isinstance(device, dict):
        key = rank if rank in device else str(rank)
        chosen = device.get(key, pool[rank % len(pool)])
    elif device is not None and policy in ("explicit", "single"):
        chosen = device
    elif policy == "single":
        chosen = pool[0]
    elif policy == "least-busy":
        chosen = _least_busy_device(pool)
    elif policy == "explicit":
        warnings.warn("policy='explicit' with no `device` given ⚠️ - falling back to round-robin.")
        chosen = pool[rank % len(pool)]
    else:
        if policy != "round-robin":
            warnings.warn(f"Unknown device policy {policy!r} ⚠️ - falling back to round-robin.")
        chosen = pool[rank % len(pool)]

    chosen = int(chosen)
    if chosen not in pool:
        warnings.warn(f"Device {chosen} is not in the visible pool {pool} ⚠️ - clamping to {pool[0]}.")
        chosen = pool[0]
    select_cuda_device.last_rank = (rank, source)   # for diagnostics / logging
    return chosen

'''# Function to set up CUDA context for DASK workers (for true GPU interop)'''
def setup_daskworker_gpu_context(policy="round-robin", device=None, devices=None):
    """Initialize CUDA context for this worker thread (quiet sibling of
    :func:`init_daskworker_cuda_context`). Returns the bound device id, or None on CPU."""
    return init_daskworker_cuda_context(policy=policy, device=device,
                                        devices=devices, verbose=True)

def init_daskworker_cuda_context(policy="round-robin", device=None, devices=None, verbose=True):
    """Bind this Dask worker to a CUDA device and warm up its context.

    Register it across a cluster with either of::

        client.run(slck.init_daskworker_cuda_context)                       # round-robin
        client.run(slck.init_daskworker_cuda_context, policy="least-busy")
        client.run(slck.init_daskworker_cuda_context, policy="explicit",
                   device={0: 0, 1: 1, 2: 0, 3: 1})                         # user-defined
        client.run(slck.init_daskworker_cuda_context,
                   device=lambda rank, pool: pool[rank % 2])                # your own rule

    Returns the bound device id (or None when running CPU-only).
    """
    if not GPU_AVAILABLE:
        warnings.warn("GPU acceleration is not available. Dask workers will run on CPU only.")
        return
    
    chosen = select_cuda_device(policy=policy, device=device, devices=devices)
    if chosen is None:
        warnings.warn("No CUDA devices visible to this worker ⚠️ - staying on CPU.")
        return None

    try:
        _cupy.cuda.Device(chosen).use()
        # Touch the device so the context is actually created, not merely selected.
        _cupy.zeros(1)
    except Exception as e:
        warnings.warn(f"Could not bind worker to device {chosen}: {e} ❌ - falling back to CPU.")
        return None

    _WORKER_DEVICE.device_id = chosen
    rank, source = getattr(select_cuda_device, "last_rank", ("?", "?"))
    if verbose:
        print(f"Worker {threading.current_thread().name} (rank {rank} via {source}) "
              f"→ CUDA device {chosen} of {visible_devices()} ✅ [policy={policy}]")
    return chosen
    
    #device = _cupy.cuda.Device(0)
    #device.use()
    #context = _cupy.cuda.runtime.getCurrentContext()
    #print(f"Worker {threading.current_thread().name}: Context ID {context}")
    # Should see different context IDs for each worker

def get_worker_device():
    """Which device did Slick bind this worker/thread to? None if unbound."""
    return getattr(_WORKER_DEVICE, "device_id", None)

'''# Function to set up CUDA context for DASK workers (for true GPU interop)
def setup_daskworker_gpu_context():
    """Initialize CUDA context for this worker thread"""
    if GPU_AVAILABLE:
        _cupy.cuda.Device(0).use()  # or assign different devices if multi-GPU
    # Store context info in thread-local storage if needed
#'''

'''# 
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
#'''
    
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
__all__ = ['np', 'numpy', 'nprandom', 'la', 'linalg',  'fft', 'signal', 'interpolate', 'stats', 'ndimage', 'init_daskworker_cuda_context',
           'to_cpu', 'to_gpu', 'asnumpy', 'asarray', 'is_gpu_array', 'setup_daskworker_gpu_context',
           'cpu_safe_import', 'ensure_cpu', 'GPU_AVAILABLE', 'USE_GPU', 'RAPIDS_AVAILABLE', 'DASK_AVAILABLE',
           # Multi-GPU worker placement 🎛️
            'select_cuda_device', 'visible_devices', 'worker_rank', 'get_worker_device']
