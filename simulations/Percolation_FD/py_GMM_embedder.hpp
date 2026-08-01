// pyGMM_embedder.hpp  (revised for Rho_t / const_index interface)
//
// Passes a 2D C++ vector (Sp x L*L) directly to gen_clustered()
// as a zero-copy NumPy view, receives scalar stat arrays and optional
// frequency array back.
//
// Link:  $(python3-config --ldflags --embed)
// CFLAGS: $(python3-config --includes)

#pragma once

#define PY_ARRAY_UNIQUE_SYMBOL BRIDGE_ARRAY_API
// Default: all includers get extern-only NumPy API references.
// The one TU that calls PythonBridge::init() must #define BRIDGE_IMPORT_ARRAY
// before including this header to own the API symbols and enable _import_array().
#ifndef BRIDGE_IMPORT_ARRAY
#define NO_IMPORT_ARRAY
#endif
#include <Python.h>
#include <numpy/arrayobject.h>

#include <cstring>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

// ---------------------------------------------------------------------------
// PythonBridge — initialise Python interpreter, import gmm_classifier.py, and get gen_clustered() object
// ---------------------------------------------------------------------------
class PythonBridge 
{
    public:
        // Only compiled in the TU that defines BRIDGE_IMPORT_ARRAY, since
        // _import_array() is only available when NO_IMPORT_ARRAY is absent.
#ifdef BRIDGE_IMPORT_ARRAY
        static void init(const std::string &script_dir)
        {
            Py_Initialize();
            if (_import_array() < 0)  // Initialising NumPy C API — must be done after Py_Initialize
            {
                PyErr_Print();
                throw std::runtime_error("Failed to initialise NumPy C API.");
            }
            PyObject *sys_path = PySys_GetObject("path");
            PyObject *pdir     = PyUnicode_FromString(script_dir.c_str());
            PyList_Append(sys_path, pdir);  // Add script_dir to Python path so we can import gmm_classifier.py
            Py_DECREF(pdir);

            s_module = PyImport_ImportModule("gmm_classifier");
            if (!s_module) { PyErr_Print(); throw std::runtime_error("Cannot import gmm_classifier.py"); }

            s_func = PyObject_GetAttrString(s_module, "gen_clustered");
            if (!s_func || !PyCallable_Check(s_func)) {
                PyErr_Print(); throw std::runtime_error("gen_clustered not Callable.");
            }
            s_save = PyEval_SaveThread(); // Release GIL until we need it in call_gen_clustered
        }
#endif

        static void finalize()
        {
            PyEval_RestoreThread(s_save);
            Py_XDECREF(s_func);
            Py_XDECREF(s_module);
            Py_Finalize();
        }

        static PyObject *func() { return s_func; }

    private:
        static PyObject      *s_module;
        static PyObject      *s_func;
        static PyThreadState *s_save;
};
inline PyObject      *PythonBridge::s_module = nullptr;
inline PyObject      *PythonBridge::s_func   = nullptr;
inline PyThreadState *PythonBridge::s_save   = nullptr;

// ---------------------------------------------------------------------------
// RAII GIL guard
// ---------------------------------------------------------------------------
struct GILGuard {
    PyGILState_STATE state;
    GILGuard()  { state = PyGILState_Ensure(); } // Constructor acquires GIL
    ~GILGuard() { PyGILState_Release(state);   }  
    // Destructor releases GIL when the GILGuard object goes out of scope.
};

// ---------------------------------------------------------------------------
// Result struct
// Scalar stat arrays are Sp-length; NaN at unselected / failed indices.
// freq_block is Sp x 3 x max_sizes (row-major), empty if not requested.
// freq_max_sizes == 0 means evaluate_clusterfrequencies was false or no
// clusters were found.
// ---------------------------------------------------------------------------
struct FrameResult {
    std::size_t Sp = 0;

    // Seven scalar arrays, each of length Sp
    std::vector<double> clus_density;
    std::vector<double> occupied_site_frac;
    std::vector<double> mean_sq_dist_lcenter;
    std::vector<double> mean_sq_dist_ucom;
    std::vector<double> mean_cluster_size;
    std::vector<double> var_cluster_size;
    std::vector<double> num_clusters;

    // Frequency block: flat row-major Sp*3*max_sizes, NaN-padded.
    // Index: [i*3*max_sizes + row*max_sizes + col]
    // row 0 = sizes, row 1 = counts, row 2 = freqs
    std::vector<double> freq_block;
    std::size_t         freq_max_sizes = 0;
    // GMM-classified frames: flat row-major Sp*(L*L), NaN at unselected rows.
    // Index: [s * L*L + site]
    std::vector<double> gmm_frames;   // empty if gen_GMMframe was false
    bool valid = false;
};

// ---------------------------------------------------------------------------
// Helper: extract a contiguous 1-D float64 NumPy array into vector<double>
// Returns vector of NaNs on failure rather than throwing, so partial results
// are still accessible.
// ---------------------------------------------------------------------------
static std::vector<double> numpy1d_to_vec(PyObject *arr, std::size_t expected_len)
{
    std::vector<double> out(expected_len, std::numeric_limits<double>::quiet_NaN());
    if (!arr || arr == Py_None) return out;
    if (!PyArray_Check(arr))    return out;
    // 
    PyArrayObject *a = reinterpret_cast<PyArrayObject *>(arr);
    // Request a C-contiguous float64 view (no copy if already compatible)
    PyArrayObject *a_cont = reinterpret_cast<PyArrayObject *>(
        PyArray_FROM_OTF(arr, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY));
    if (!a_cont) { PyErr_Clear(); return out; }

    npy_intp n = PyArray_SIZE(a_cont);
    std::size_t copy_n = std::min(static_cast<std::size_t>(n), expected_len);
    std::memcpy(out.data(), PyArray_DATA(a_cont), copy_n * sizeof(double));
    Py_DECREF(a_cont);
    return out;
}

/**  ---------------------------------------------------------------------------
// Main callable — pass Rho_t and const_index, get FrameResult back.
//
// rho_t        : Sp x (L*L) matrix, row-major, caller owns memory.
//                Passed as a zero-copy 2-D NumPy view — must remain valid
//                for the duration of this call.
// const_index  : indices into axis-0 of rho_t to classify
// L            : grid side length
// aval/Tval/Rval : forwarded to Python for logging
// n_clusters   : GMM components
// periodic     : periodic BCs for ndimage.label
// eval_freq    : whether to compute and return cluster frequency array (freq_block in FrameResult, size max_sizes*(Sp*3))
// gen_GMMframe  : whether to compute and return GMM-classified frames (Sp x L*L, row-major, NaN at unselected rows).
Returns FrameResult with scalar stat arrays and optional frequency block.
The scalar stat arrays will be NaN at any indices where classification failed, and are of size Sp
freq_block is Sp x 3 x max_sizes (row-major), empty if not requested.
freq_max_sizes == 0 means evaluate_clusterfrequencies was false or no clusters were found.
// --------------------------------------------------------------------------- */
static FrameResult call_gen_clustered(const std::vector<std::vector<double>> &rho_t, const std::vector<int> &const_index,
    int L, const double aval, const double Tval, const int Rval, 
    int  n_clusters= 2, bool periodic = true, bool eval_freq= true, bool gen_GMMframe = false)
{
    FrameResult result;
    const std::size_t Sp = rho_t.size();
    result.Sp = Sp;
    if (Sp == 0) return result;

    GILGuard gil;

    // ---- Build 2-D NumPy array from rho_t without copying ----------------
    // Armadillo / simulation code stores each species as a contiguous row.
    // We need a single contiguous buffer for a 2-D NumPy array, so we copy
    // row pointers into a flat staging buffer here.
    // If your simulation already stores Rho_t as a flat contiguous block
    // (e.g. std::vector<double> of size Sp*(L*L)), you can avoid this copy
    // by passing that pointer directly to PyArray_SimpleNewFromData instead.
    const std::size_t row_len = static_cast<std::size_t>(L) * L;
    // Verify all rows have the expected length before touching Python.
    for (std::size_t s = 0; s < Sp; ++s) 
    {
        if (rho_t[s].size() != row_len) 
        {
            PyErr_SetString(PyExc_ValueError, "rho_t row length != L*L");
            return result;
        }
    }

    // Flat staging buffer — one allocation, row-major.
    std::vector<double> flat(Sp * row_len);
    for (std::size_t s = 0; s < Sp; ++s)
        std::memcpy(flat.data() + s * row_len, rho_t[s].data(), row_len * sizeof(double));

    npy_intp dims_2d[2] = { static_cast<npy_intp>(Sp),
                             static_cast<npy_intp>(row_len) };
    PyObject *py_rho = PyArray_SimpleNewFromData(2, dims_2d, NPY_DOUBLE, flat.data());
    if (!py_rho) { PyErr_Print(); return result; }

    // ---- Build 1-D int (NPY_INTP) array for const_index ------------------
    npy_intp dims_idx[1] = { static_cast<npy_intp>(const_index.size()) };
    // Copy into a local intp buffer so element width matches what NumPy expects
    std::vector<npy_intp> idx_buf(const_index.begin(), const_index.end());
    PyObject *py_idx = PyArray_SimpleNewFromData(1, dims_idx, NPY_INTP, idx_buf.data());
    if (!py_idx) { Py_DECREF(py_rho); PyErr_Print(); return result; }

    // ---- Call gen_clustered(Rho_t, const_index, L,
    //                                aval, Tval, Rval,
    //                                n_clusters, periodic,
    //                                evaluate_clusterfrequencies) ----------
    PyObject *py_ret = PyObject_CallFunction(
        PythonBridge::func(), "OOiddiiOOO",
        py_rho,
        py_idx,
        L,
        aval,
        Tval,
        Rval,
        n_clusters,
        periodic  ? Py_True : Py_False,
        eval_freq ? Py_True : Py_False,
        gen_GMMframe  ? Py_True : Py_False
    );

    Py_DECREF(py_rho);
    Py_DECREF(py_idx);

    if (!py_ret) { PyErr_Print(); return result; }

    // py_ret is a 3-tuple: (scalar_stats_dict, freq_array_or_None, gmm_frames_or_None)
    if (!PyTuple_Check(py_ret) || PyTuple_GET_SIZE(py_ret) != 3) {
        PyErr_SetString(PyExc_TypeError, "gen_clustered must return a 3-tuple");
        Py_DECREF(py_ret);
        return result;
    }

    PyObject *py_stats = PyTuple_GET_ITEM(py_ret, 0);   // dict
    PyObject *py_freq  = PyTuple_GET_ITEM(py_ret, 1);   // ndarray or None
    PyObject *py_frames = PyTuple_GET_ITEM(py_ret, 2); // ndarray or None, only returned if gen_GMMframe is true

    // ---- Extract scalar stat arrays from the dict ------------------------
    // Key order matches STAT_KEYS in gmm_classifier.py
    auto extract = [&](const char *key) -> std::vector<double> {
        PyObject *arr = PyDict_GetItemString(py_stats, key);  // borrowed ref
        return numpy1d_to_vec(arr, Sp);
    };

    result.clus_density          = extract("clus_density");
    result.occupied_site_frac    = extract("occupied_site_frac");
    result.mean_sq_dist_lcenter  = extract("mean_sq_dist_lcenter");
    result.mean_sq_dist_ucom     = extract("mean_sq_dist_ucom");
    result.mean_cluster_size     = extract("mean_cluster_size");
    result.var_cluster_size      = extract("var_cluster_size");
    result.num_clusters          = extract("num_clusters");

    // ---- Extract frequency block (Sp x 3 x max_sizes) -------------------
    if (eval_freq && py_freq != Py_None && PyArray_Check(py_freq)) {
        PyArrayObject *fa = reinterpret_cast<PyArrayObject *>(
            PyArray_FROM_OTF(py_freq, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY));
        if (fa) {
            if (PyArray_NDIM(fa) == 3) {
                npy_intp *sh = PyArray_SHAPE(fa);
                // sh[0] = Sp, sh[1] = 3, sh[2] = max_sizes
                std::size_t total = static_cast<std::size_t>(sh[0])
                                  * static_cast<std::size_t>(sh[1])
                                  * static_cast<std::size_t>(sh[2]);
                result.freq_block.resize(total);
                std::memcpy(result.freq_block.data(),
                            PyArray_DATA(fa),
                            total * sizeof(double));
                result.freq_max_sizes = static_cast<std::size_t>(sh[2]);
            }
            Py_DECREF(fa);
        }
    }

    // GMM frames block: shape (Sp, L*L)
    if (gen_GMMframe && py_frames != Py_None && PyArray_Check(py_frames)) {
        PyArrayObject *ff = reinterpret_cast<PyArrayObject *>(
            PyArray_FROM_OTF(py_frames, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY));
        if (ff) {
            if (PyArray_NDIM(ff) == 2) {
                std::size_t total = static_cast<std::size_t>(Sp) * row_len;
                result.gmm_frames.resize(total);
                std::memcpy(result.gmm_frames.data(), PyArray_DATA(ff), total * sizeof(double));
            }
            Py_DECREF(ff);
        }
    }

    Py_DECREF(py_ret);
    result.valid = true;
    return result;
}