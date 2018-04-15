#include <gtest/gtest.h>

#include "c4/cfd/cfd.hpp"

template< int N > using soa = c4::cfd::soa< N, double, uint32_t >;
template< int N > using aos = c4::cfd::aos< N, double, uint32_t >;

template< class S > using inc = c4::cfd::dvars::incompressible< S >;

template< class S > using cart = c4::cfd::grid::cartesian< S >;
template< class S > using unst = c4::cfd::grid::unstructured< S >;

template< int N, template<class> class G, template<int> class S > using prob = c4::cfd::problem< inc, G, S<N> >;

template< class Probl >
void test_instantiation()
{
    Probl p;
}

TEST(cart, soa1d)
{
    prob< 1, cart, soa > p;
}

TEST(cart, soa2d)
{
    prob< 2, cart, soa > p;
}

TEST(cart, soa3d)
{
    prob< 2, cart, soa > p;
}

TEST(cart, aos1d)
{
    prob< 1, cart, aos > p;
}

TEST(cart, aos2d)
{
    prob< 2, cart, aos > p;
}

TEST(cart, aos3d)
{
    prob< 2, cart, aos > p;
}



TEST(unst, soa1d)
{
    prob< 1, unst, soa > p;
}

TEST(unst, soa2d)
{
    prob< 2, unst, soa > p;
}

TEST(unst, soa3d)
{
    prob< 2, unst, soa > p;
}

TEST(unst, aos1d)
{
    prob< 1, unst, aos > p;
}

TEST(unst, aos2d)
{
    prob< 2, unst, aos > p;
}

TEST(unst, aos3d)
{
    prob< 2, unst, aos > p;
}
