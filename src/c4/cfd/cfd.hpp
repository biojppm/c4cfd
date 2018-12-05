#ifndef _C4_CFD_HPP_
#define _C4_CFD_HPP_

/** @file grid.hpp this is a skeleton sketch for a CFD solver:
 * -custom number of dimensions
 * -custom storage: Structure-Of-Arrays vs Array-Of-Structures
 * -row-major vs col-major tensors
 * -cartesian and unstructured grids
 * */

#include "c4/config.hpp"
#include "c4/memory_resource.hpp"

#include <vector>


namespace c4 {

// utilities to tell the compiler the memory is SIMD-aligned

/** @todo this assumes AVX512. Needs fix. */
constexpr size_t get_simd_size()
{
    return size_t(64);
}

constexpr const size_t simd_alignment = get_simd_size(); ///< align for AVX512.

#define C4_SIMD_ALIGNED(ptr) C4_ASSUME_ALIGNED(ptr, simd_alignment);
#define C4_SIMD_ALIGNED_OFFS(ptr, offs) C4_ASSUME_ALIGNED_OFFS(ptr, simd_alignment. offs);

#define C4_ASSUME_ALIGNED(ptr, align) __builtin_assume_aligned(ptr, align)
#define C4_ASSUME_ALIGNED_OFFS(ptr, align, offs) __builtin_assume_aligned(ptr, align, offs)

} // namespace c4


// make all pointers restricted by default
#define  $         * __restrict__ ///< a restricted pointer
#define $$         & __restrict__ ///< a restricted reference
#define c$   const * __restrict__ ///< a restricted const pointer
#define c$$  const & __restrict__ ///< a restricted const reference



namespace c4 {
namespace cfd {

template<int N, typename T>
struct vec_
{
    T data[N];
};

template<typename T>
struct vec_<1,T>
{
    union {
        T data[1];
        T x;
    };
};

template<typename T>
struct vec_<2,T>
{
    union {
        T data[2];
        T x, y;
    };
};

template<typename T>
struct vec_<3,T>
{
    union {
        T data[3];
        T x, y, z;
    };
};

template<typename T>
struct vec_<4,T>
{
    union {
        T data[4];
        T x, y, z, w;
    };
};

template<int N, typename T>
struct vec : public vec_<N,T>
{
    template <class I> T  $$ operator[] (I i)       { C4_XASSERT(i >= 0 && i < N); return this->data[N]; }
    template <class I> T c$$ operator[] (I i) const { C4_XASSERT(i >= 0 && i < N); return this->data[N]; }
};


//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------


/** @todo assert that T is memcpy-able */
template< typename T, typename I=size_t >
struct memblock
{
    T $  m_val;
    I    m_size;
    bool m_owner;

    inline T  $$ operator[] (I i)       { C4_ASSERT(i >= 0 && i < m_size); return m_val[i]; }
    inline T c$$ operator[] (I i) const { C4_ASSERT(i >= 0 && i < m_size); return m_val[i]; }

    memblock() : m_val(nullptr), m_size(0), m_owner(false) {}
    memblock(I sz) : memblock() { resize(sz); }
    memblock(I sz, T $ arr) : memblock() { borrow(sz, arr); }
    ~memblock() { release(); }

    // TODO rule of 5

    void borrow(I sz, T $ mem)
    {
        if(m_owner)
        {
            release();
        }
        m_size = sz;
        m_val = mem;
        m_owner = false;
    }

    void release()
    {
        if(m_owner)
        {
            c4::afree(m_val);
        }
        m_val = nullptr;
        m_size = 0;
        m_owner = false;
    }

    void resize(I sz)
    {
        if(sz == m_size) return;
        T *mem = C4_ASSUME_ALIGNED(c4::aalloc(sz * sizeof(T), simd_alignment), simd_alignment);
        if(m_val)
        {
            T *aval = C4_ASSUME_ALIGNED(m_val, simd_alignment);
            I min = sz < m_size ? sz : m_size;
            memcpy(mem, aval, min);
        }
        if(m_owner)
        {
            release();
        }
        m_val = mem;
        m_size = sz;
        m_owner = true;
    }

    // etc
};


//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

namespace detail {

template< typename I, size_t N >
struct mat_addr_indices;

template< typename I >
struct mat_addr_indices<I, 1>
{
    using mpos = I;
    static constexpr const I    s_sym2lin[1][1] = {{0}};
    static constexpr const mpos s_lin2sym_rm[1] = {0};
    static constexpr const mpos s_lin2sym_cm[1] = {0};
};

template< typename I >
struct mat_addr_indices<I, 2>
{
    using mpos = vec<2,I>;
    static constexpr const I    s_sym2lin[2][2] = {{0,1},{1,2}};
    static constexpr const mpos s_lin2sym_rm[3] = {{0,0},{0,1},{1,1}};
    static constexpr const mpos s_lin2sym_cm[3] = {{0,0},{1,0},{1,1}};
};

template< typename I >
struct mat_addr_indices<I, 3>
{
    using mpos = vec<2,I>;
    static constexpr const I    s_sym2lin[3][3] = {{0,1,2}, {1,3,4}, {2,4,5}};
    static constexpr const mpos s_lin2sym_rm[6] = {{0,0},{0,1},{0,2},{1,1},{1,2},{2,2}};
    static constexpr const mpos s_lin2sym_cm[6] = {{0,0},{1,0},{2,0},{1,1},{2,1},{2,2}};
};

} // namespace detail


/** dense matrix addressing
 * lin := linear
 * rm := row major
 * cm := col major
 * sym := symmetric
 */

template<typename I, I N>
inline I rm2lin(I i, I j)
{
    return (i * N) + j;
}
template<typename I, I N>
inline vec<2,I> lin2rm(vec<2,I> p)
{
    return {p/N, p%N};
}

template< typename I, I N>
inline I cm2lin(I i, I j)
{
    return (j * N) + i;
}
template<typename I, I N>
inline vec<2,I> lin2cm(vec<2,I> p)
{
    return {p%N, p/N};
}

// https://stackoverflow.com/questions/19143657/linear-indexing-in-symmetric-matrices
template<typename I, I N>
inline I sym2lin(I i, I j)
{
    C4_XASSERT(i < 3);
    C4_XASSERT(j < 3);
    return detail::mat_addr_indices<I,N>::s_sym2lin[i][j];
}

template<typename I, I N>
inline I sym2lin(vec<2,I> i)
{
    return sym2lin<I,N>(i.x, i.y);
}

template<typename I, I N>
inline vec<2,I> lin2sym_rm(I p)
{
    C4_XASSERT(p < 6);
    return detail::mat_addr_indices<I,N>::s_lin2sym_rm[p];
}

template<typename I, I N>
inline vec<2,I> lin2sym_cm(I p)
{
    C4_XASSERT(p < 6);
    return detail::mat_addr_indices<I,N>::s_lin2sym_cm[p];
}


template<typename I, I N>
struct row_major
{
    using mpos = vec<2,I>;
    static inline I mpos2lin(I i, I j) { return rm2lin<I,N>(i, j); }
    static inline I mpos2lin(mpos i) { return rm2lin<I,N>(i.x, i.y); }
    static inline I mpos2lin_sym(I i, I j) { return sym2lin<I,N>(i, j); }
    static inline I mpos2lin_sym(mpos i) { return sym2lin<I,N>(i.x, i.j); }
    static inline mpos lin2mpos(I i) { return lin2rm<I,N>(i); }
    static inline mpos lin2mpos_sym(I i) { return lin2sym_rm<I,N>(i); }
};

template<typename I, I N>
struct col_major
{
    using mpos = vec<2,I>;

    static inline I mpos2lin(I i, I j) { return cm2lin<I,N>(i, j); }
    static inline I mpos2lin(mpos i) { return cm2lin<I,N>(i.x, i.y); }
    static inline I mpos2lin_sym(I i, I j) { return sym2lin<I,N>(i, j); }
    static inline I mpos2lin_sym(mpos i) { return sym2lin<I,N>(i.x, i.j); }
    static inline mpos lin2mpos(I i) { return lin2cm<I,N>(i); }
    static inline mpos lin2mpos_sym(I i) { return lin2sym_cm<I,N>(i); }
};


//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

template< typename T=double, typename I=size_t >
struct _var
{
    using value_type = T;
    using index_type = I;
};



//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

template< int N, typename T=double, typename I=size_t >
struct basic_var;


/** scalar */
template< typename T, typename I >
struct basic_var<1,T,I>
{
    memblock<T,I> m_val;

    using mpos = vec<1,I>;

    inline T  $$ operator() (I elm)       { return this->m_val[elm]; }
    inline T c$$ operator() (I elm) const { return this->m_val[elm]; }

    inline T  $$ operator() (I elm, mpos dim)       { C4_UNUSED(dim); C4_ASSERT(dim.x == 0); return this->m_val.m_val[elm]; }
    inline T c$$ operator() (I elm, mpos dim) const { C4_UNUSED(dim); C4_ASSERT(dim.x == 0); return this->m_val.m_val[elm]; }
};

/** non-scalar */
template< int N, typename T, typename I >
struct basic_var
{
    memblock<T,I> m_val[N];

    using mpos = vec<1,I>;

    inline T  $$ operator() (I elm, mpos dim)       { C4_UNUSED(dim); C4_ASSERT(dim >= 0 && dim.x < N); return this->m_val[dim].m_val[elm]; }
    inline T c$$ operator() (I elm, mpos dim) const { C4_UNUSED(dim); C4_ASSERT(dim >= 0 && dim.x < N); return this->m_val[dim].m_val[elm]; }
};


//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

#define C4_STORAGE_TYPES()                                  \
    using T = typename Storage::value_type;                 \
    using I = typename Storage::value_type;                 \
    static constexpr const I num_dims = Storage::num_dims;  \
    using value_type = typename Storage::value_type;        \
    using index_type = typename Storage::index_type;        \
    using scalar = typename Storage::scalar;                \
    using vector = typename Storage::vector;                \
    using tensor_rm = typename Storage::tensor_rm;          \
    using tensor_cm = typename Storage::tensor_cm;          \
    using symtensor_rm = typename Storage::symtensor_rm;    \
    using symtensor_cm = typename Storage::symtensor_cm


//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

template< int D, typename T=double, typename I=size_t >
struct aos
{
    static constexpr const I num_dims = D;

    using value_type = T;
    using index_type = I;

    using mpos = vec<D,I>;

    struct scalar : public _var<T,I>
    {
        enum : int { N = 1 };
        basic_var<1, T, I> val;

        inline T  $$ operator() (I elm)       { return val(elm); }
        inline T c$$ operator() (I elm) const { return val(elm); }

        inline T  $$ operator() (I elm, I dim)       { C4_XASSERT(dim == 0); return val(elm); }
        inline T c$$ operator() (I elm, I dim) const { C4_XASSERT(dim == 0); return val(elm); }

        inline T  $$ operator() (I elm, mpos dim)       { C4_XASSERT(dim.x == 0); return val(elm); }
        inline T c$$ operator() (I elm, mpos dim) const { C4_XASSERT(dim.x == 0); return val(elm); }
    };

    struct vector : public _var<T, I>
    {
        enum : int { N = D };
        basic_var<1, T, I> val;

        using storage = basic_var<1,T,I>;
        using mpos = typename storage::mpos;

        inline T  $$ operator() (I elm, I dim)       { return val(elm * D + dim); }
        inline T c$$ operator() (I elm, I dim) const { return val(elm * D + dim); }

        inline T c$$ operator() (I elm, mpos dim)       { return val(elm * D + dim); }
        inline T c$$ operator() (I elm, mpos dim) const { return val(elm * D + dim); }
    };

    template<class Addr>
    struct tensor : public _var<T, I>
    {
        enum : int { N = D*D };
        basic_var<1, T, I> val;

        using storage = basic_var<1, T, I>;
        using mpos = typename storage::mpos;
        using addr = Addr;

        inline T  $$ operator() (I elm, I dim1, I dim2)       { return val(elm * N + addr::mpos2lin(dim1, dim2)); }
        inline T c$$ operator() (I elm, I dim1, I dim2) const { return val(elm * N + addr::mpos2lin(dim1, dim2)); }

        inline T  $$ operator() (I elm, mpos dim)       { return val(elm * N + addr::mpos2lin(dim)); }
        inline T c$$ operator() (I elm, mpos dim) const { return val(elm * N + addr::mpos2lin(dim)); }
    };

    template<class Addr>
    struct symtensor : public _var<T, I>
    {
        enum : int { N = D*(D+1)/2 };
        basic_var<1, T, I> val;

        using storage = basic_var<1, T, I>;
        using mpos = typename storage::mpos;
        using addr = Addr;

        inline T  $$ operator() (I elm, I dim1, I dim2)       { return val(elm * N + addr::mpos2lin_sym(dim1, dim2)); }
        inline T c$$ operator() (I elm, I dim1, I dim2) const { return val(elm * N + addr::mpos2lin_sym(dim1, dim2)); }

        inline T  $$ operator() (I elm, mpos dim)       { return val(elm * N + addr::mpos2lin_sym(dim)); }
        inline T c$$ operator() (I elm, mpos dim) const { return val(elm * N + addr::mpos2lin_sym(dim)); }
    };

    using tensor_rm = tensor<row_major<I,D>>;
    using tensor_cm = tensor<col_major<I,D>>;

    using symtensor_rm = symtensor<row_major<I,D>>;
    using symtensor_cm = symtensor<col_major<I,D>>;
};



//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

template< int D, typename T=double, typename I=size_t >
struct soa
{
    static constexpr const I num_dims = D;

    using value_type = T;
    using index_type = I;

    using mpos = vec<D,I>;

    struct scalar : public _var<T, I>
    {
        enum : int { N = 1 };
        basic_var<N, T, I> val;

        inline T  $$ operator() (I elm)       { return val(elm); }
        inline T c$$ operator() (I elm) const { return val(elm); }

        inline T  $$ operator() (I elm, I dim)       { return val(elm); }
        inline T c$$ operator() (I elm, I dim) const { return val(elm); }

        inline T  $$ operator() (I elm, mpos dim)       { return val(elm); }
        inline T c$$ operator() (I elm, mpos dim) const { return v7al(elm); }
    };

    struct vector : public _var<T, I>
    {
        enum : int { N = D };
        basic_var<N, T, I> val;

        using storage = basic_var<D,T,I>;

        inline T  $$ operator() (I elm, I dim)       { return val(elm, dim); }
        inline T c$$ operator() (I elm, I dim) const { return val(elm, dim); }

        inline T c$$ operator() (I elm, mpos dim)       { return val(elm, dim); }
        inline T c$$ operator() (I elm, mpos dim) const { return val(elm, dim); }
    };

    template<class Addr>
    struct tensor : public _var<T, I>
    {
        enum : int { N = D*D };
        basic_var<N, T, I> val;

        using storage = basic_var<N, T, I>;
        using mpos = typename storage::mpos;
        using addr = Addr;

        inline T  $$ operator() (I elm, I dim1, I dim2)       { return val(elm, addr::mpos2lin(dim1, dim2)); }
        inline T c$$ operator() (I elm, I dim1, I dim2) const { return val(elm, addr::mpos2lin(dim1, dim2)); }

        inline T  $$ operator() (I elm, mpos dim)       { return val(elm, addr::mpos2lin(dim)); }
        inline T c$$ operator() (I elm, mpos dim) const { return val(elm, addr::mpos2lin(dim)); }
    };

    template<class Addr>
    struct symtensor : public _var<T, I>
    {
        enum : int { N = D*(D+1)/2 };
        basic_var<N, T, I> val;

        using storage = basic_var<N, T, I>;
        using mpos = typename storage::mpos;
        using addr = Addr;

        inline T  $$ operator() (I elm, I dim1, I dim2)       { return val(elm, addr::mpos2lin_sym(dim1, dim2)); }
        inline T c$$ operator() (I elm, I dim1, I dim2) const { return val(elm, addr::mpos2lin_sym(dim1, dim2)); }

        inline T  $$ operator() (I elm, mpos dim)       { return val(elm, addr::mpos2lin_sym(dim)); }
        inline T c$$ operator() (I elm, mpos dim) const { return val(elm, addr::mpos2lin_sym(dim)); }
    };

    using tensor_rm = tensor<row_major<I,D>>;
    using tensor_cm = tensor<col_major<I,D>>;

    using symtensor_rm = symtensor<row_major<I,D>>;
    using symtensor_cm = symtensor<col_major<I,D>>;

};



//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

template< typename I >
struct adjlist
{

    memblock<I, I> m_pos;
    memblock<I, I> m_adj;

    adjlist() : m_pos(), m_adj() {}
    adjlist(I num_elms, I num_adj_elms) : adjlist() { resize(num_elms, num_adj_elms); }

    I num_elms() const { I sz = m_pos.size(); return sz > 0 ? sz-1 : sz; }

    void resize(I num_elms, I num_adj_elms)
    {
        m_pos.resize(num_elms + 1);
        m_adj.resize(num_adj_elms);
    }

    struct adj_iter
    {
        I c$ b;
        I c$ e;

        adj_iter(I c$ b_, I c$ e_) : b(b_), e(e_) {}

        I c$ begin() { return *b; }
        I c$ end()   { return *e; }
    };

    adj_iter adj(I elm) const
    {
        C4_XASSERT(elm > 0 && elm < num_elms());
        return adj_iter(&m_adj[m_pos[elm]], &m_adj[m_pos[elm+1]]);
    }
};


//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

namespace grid {

template< class I >
struct face_cell
{
    I left, right;
};

template< class Storage >
struct cartesian
{
    C4_STORAGE_TYPES();

    using face_type = face_cell<I>;

    memblock<face_type, I> m_boundaries;

    vector m_vert_coords;
};

template< class Storage >
struct unstructured
{
    C4_STORAGE_TYPES();

    using face_type = I;

    I m_ncells;
    I m_nfaces;
    I m_nverts;

    adjlist<I>   m_cell_faces;
    adjlist<I>   m_vert_cells;
    adjlist<I>   m_cell_verts;

    memblock<face_cell<I>, I> m_face_cells;

    memblock<face_type, I> m_boundaries;

    vector m_vert_coords;

    vector m_cell_center;
    scalar m_cell_vol;

    vector m_face_center;
    vector m_face_nrml;

};

} // namespace grid


//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

namespace bnd {

typedef enum
{
    DIRICHLET,
    NEUMANN
} MathBoundaryType_e;


template<class T, class I>
struct BoundaryValues
{
    T m_single_val;
    memblock<T, I> m_val_per_face;
};


template<class T, class I>
struct Dirichlet : public BoundaryValues<T,I>
{
};

template<class T, class I>
struct Neumann : public BoundaryValues<T,I>
{
};

template<class T, class I>
struct MathBoundary
{
    MathBoundaryType_e m_math_type;
    union {
        Dirichlet<T,I> m_dirichlet;
        Neumann<T,I> m_neumann;
    };
};

typedef enum BoundaryType_e
{
    WALL,
    INFLOW,
    OUTFLOW,
    INTERIOR,
    AMR,
    CUSTOM
} BoundaryType_e;


template< class DependentVars >
struct Boundary
{
    using T = typename DependentVars::T;
    using I = typename DependentVars::I;
    using value_type = typename DependentVars::T;
    using index_type = typename DependentVars::I;

    BoundaryType_e       m_type;
    MathBoundary<T, I>   m_velocity;
    MathBoundary<T, I>   m_pressure;
};

} // namespace bnd


//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

template
<
    template<class> class DependentVars,
    template<class> class Grid,
    class Storage
>
struct problem : public DependentVars< Storage >
{
    C4_STORAGE_TYPES();

    using vars_type = DependentVars< Storage >;
    using grid_type = Grid< Storage >;
    using bnd_type = bnd::Boundary< DependentVars< Storage > >;

    Grid< Storage > m_grid;
    memblock<bnd_type, I> m_boundaries;
};


//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

namespace dvars {

template<class Storage>
struct incompressible
{
    C4_STORAGE_TYPES();

    vector m_velocity;
    scalar m_pressure;
    T      m_mu;
};

template<class Storage>
struct compressible
{
    C4_STORAGE_TYPES();

    scalar m_density;
    vector m_velocity;
    scalar m_pressure;
    scalar m_temperature;
    scalar m_energy;
    T      m_mu;
};

} // namespace dvars


//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

/** @namespace amr adaptive mesh refinement */
namespace amr {

template< class ProblemImpl >
struct amr_problem_node
{
    ProblemImpl *m_parent;
    std::vector< ProblemImpl > m_children;

    amr_problem_node(ProblemImpl *parent_=nullptr) : m_parent(parent_), m_children() {}
};

template< class ProblemImpl >
struct amr_problem
{
    amr_problem_node< ProblemImpl > m_root;
};

} // namespace amr

} // namespace cfd
} // namespace c4

#endif // _C4_CFD_HPP_
