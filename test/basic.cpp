#include <gtest/gtest.h>

#include "c4/cfd/cfd.hpp"

template<int N, class T, class I> using soa = c4::cfd::soa<N, T, I>;
template<int N, class T, class I> using aos = c4::cfd::aos<N, T, I>;

template<class S> using inc   = c4::cfd::dvars::incompressible<S>;
template<class S> using compr = c4::cfd::dvars::compressible<S>;

template<class S> using cart = c4::cfd::grid::cartesian<S>;
template<class S> using unst = c4::cfd::grid::unstructured<S>;

template
<
    int N,
    template<class> class DVars,
    template<class> class Grid,
    template<int,class,class> class Storage,
    class T=double,
    class I=int32_t
>
void test_instantiation()
{
    c4::cfd::problem<DVars, Grid, Storage<N,T,I>> prob;
}

TEST(cart_inc, soa1d)
{
    test_instantiation<1, inc, cart, soa>();
}

TEST(cart_inc, soa2d)
{
    test_instantiation<2, inc, cart, soa>();
}

TEST(cart_inc, soa3d)
{
    test_instantiation<3, inc, cart, soa>();
}

TEST(cart_inc, aos1d)
{
    test_instantiation<1, inc, cart, aos>();
}

TEST(cart_inc, aos2d)
{
    test_instantiation<2, inc, cart, aos>();
}

TEST(cart_inc, aos3d)
{
    test_instantiation<3, inc, cart, aos>();
}

TEST(cart_compr, soa1d)
{
    test_instantiation<1, compr, cart, soa>();
}

TEST(cart_compr, soa2d)
{
    test_instantiation<2, compr, cart, soa>();
}

TEST(cart_compr, soa3d)
{
    test_instantiation<3, compr, cart, soa>();
}

TEST(cart_compr, aos1d)
{
    test_instantiation<1, compr, cart, aos>();
}

TEST(cart_compr, aos2d)
{
    test_instantiation<2, compr, cart, aos>();
}

TEST(cart_compr, aos3d)
{
    test_instantiation<3, compr, cart, aos>();
}




TEST(unst_inc, soa1d)
{
    test_instantiation<1, inc, unst, soa>();
}

TEST(unst_inc, soa2d)
{
    test_instantiation<2, inc, unst, soa>();
}

TEST(unst_inc, soa3d)
{
    test_instantiation<3, inc, unst, soa>();
}

TEST(unst_inc, aos1d)
{
    test_instantiation<1, inc, unst, aos>();
}

TEST(unst_inc, aos2d)
{
    test_instantiation<2, inc, unst, aos>();
}

TEST(unst_inc, aos3d)
{
    test_instantiation<3, inc, unst, aos>();
}

TEST(unst_compr, soa1d)
{
    test_instantiation<1, compr, unst, soa>();
}

TEST(unst_compr, soa2d)
{
    test_instantiation<2, compr, unst, soa>();
}

TEST(unst_compr, soa3d)
{
    test_instantiation<3, compr, unst, soa>();
}

TEST(unst_compr, aos1d)
{
    test_instantiation<1, compr, unst, aos>();
}

TEST(unst_compr, aos2d)
{
    test_instantiation<2, compr, unst, aos>();
}

TEST(unst_compr, aos3d)
{
    test_instantiation<3, compr, unst, aos>();
}
