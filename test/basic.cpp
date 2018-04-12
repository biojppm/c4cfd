#include "cfd.hpp"

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

int main(int argc, char* argv[])
{

    // cartesian mesh, structure of arrays, 1D, 2D, 3D
    {
        prob< 1, cart, soa > p;
    }
    {
        prob< 2, cart, soa > p;
    }
    {
        prob< 3, cart, soa > p;
    }

    // cartesian mesh, array of structures, 1D, 2D, 3D
    {
        prob< 1, cart, aos > p;
    }
    {
        prob< 2, cart, aos > p;
    }
    {
        prob< 3, cart, aos > p;
    }

    // unstructured mesh, structure of arrays, 1D, 2D, 3D
    {
        prob< 1, unst, soa > p;
    }
    {
        prob< 2, unst, soa > p;
    }
    {
        prob< 3, unst, soa > p;
    }

    // unstructured mesh, array of structures, 1D, 2D, 3D
    {
        prob< 1, unst, aos > p;
    }
    {
        prob< 2, unst, aos > p;
    }
    {
        prob< 3, unst, aos > p;
    }

    return 0;
}
