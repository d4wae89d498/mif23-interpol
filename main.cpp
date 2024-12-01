#include <iostream>
#include <cassert>
#include <cmath>

using namespace std;

uint factorial(uint n) {
    if (n <= 1)
        return 1;
    return (n * factorial(n - 1));
}

template <uint dim>
struct Vecteur
{
    double coord[dim];

    Vecteur() {};

    Vecteur(double given_coord[dim]) {
        for (uint i = 0; i < dim; ++i)
            coord[i] = given_coord[i];
    }

    Vecteur(Vecteur<dim> V1, Vecteur<dim> V2) {
        for (uint i = 0; i < dim; ++i)
            coord[i] = V2.coord[i] - V1.coord[i];
    }

    void print() const {
        cout << "{";
        for (uint i = 0; i < dim; ++i) {
            cout << coord[i];
            if (i < dim - 1) cout << ", ";
        }
        cout << "}" << endl;
    }
};

template <uint dim>
Vecteur<dim - 1> tronq(Vecteur<dim> vec, uint k) {
    double coord[dim - 1];
    uint y = 0;
    for (uint i = 0; i < dim; ++i) {
        if (k != i) {
            coord[y] = vec.coord[i];
            ++y;
        }
    }
    return Vecteur<dim - 1>(coord);
}

template <uint dim = 2>
double det(Vecteur<2> vec[2]) {
    return vec[0].coord[0] * vec[1].coord[1] - vec[0].coord[1] * vec[1].coord[0];
}

template <uint dim>
double det(Vecteur<dim> vec[dim]) {
    double output = 0;
    for (uint k = 0; k < dim; ++k) {
        Vecteur<dim - 1> remaining[dim - 1];
        for (uint l = 1; l < dim; ++l)
            remaining[l - 1] = tronq(vec[l], k);
        output += pow(-1., k) * vec[0].coord[k] * det<dim - 1>(remaining);
    }
    return output;
}

template <uint dim>
struct Simplexe
{
    Vecteur<dim> sommets[dim + 1];

    Simplexe(Vecteur<dim> given_sommets[dim + 1]) {
        for (uint i = 0; i < dim + 1; ++i)
            sommets[i] = given_sommets[i];
    }

    double volume() const {
        Vecteur<dim> psk[dim];
        for (uint i = 0; i < dim; ++i)
            psk[i] = Vecteur<dim>(sommets[0], sommets[i + 1]);
        return abs(det<dim>(psk)) / double(factorial(dim));
    }

  	double barycentrique_i(Vecteur<dim> Q, uint i)
	{
		Simplexe<dim> test = *this;
		test.sommets[i] = Q;

		double vol_full = volume();
		double vol_partial = test.volume();

		return vol_partial / vol_full;
	}

	double interpol(Vecteur<dim> Q, double fi[dim + 1])
	{
		double output = 0;

		for (uint i = 0; i < dim + 1; ++i)
		{
			double barycentric_coord = barycentrique_i(Q, i);
			output += barycentric_coord * fi[i];
		}

		return output;
	}

	bool contains(Vecteur<dim> Q)
	{
		double sum_bary = 0;

		for (uint i = 0; i < dim + 1; ++i)
		{
			double bary = barycentrique_i(Q, i);
			if (bary < 0)
				return false;
			sum_bary += bary;
		}
		return fabs(sum_bary - 1.0) < 1e-6;
	}

};

int main() {
    assert(det((Vecteur<2>[2]){
        (double[2]){1., 0.},
        (double[2]){0., 1.}
    }) == 1.);

    assert(det<3>((Vecteur<3>[3]){
        (double[3]){1., 0., 0.},
        (double[3]){0., 1., 0.},
        (double[3]){0., 0., 1.}
    }) == 1.);

    cout << "Determinant tests passed!" << endl;

 	Simplexe<2> s = Simplexe<2>(
        (Vecteur<2>[3]){
            Vecteur<2>((double[2]){0., 0.}),
            Vecteur<2>((double[2]){0., 10.}),
            Vecteur<2>((double[2]){10., 0.})
        }
    );

    assert(s.contains(Vecteur<2>((double[2]){5., 5.})) == true);
    assert(s.contains(Vecteur<2>((double[2]){15., 5.})) == false);
    assert(s.contains(Vecteur<2>((double[2]){2., 3.})) == true);
    assert(s.contains(Vecteur<2>((double[2]){8., 9.})) == false);
    assert(s.contains(Vecteur<2>((double[2]){0., 0.})) == true);

    cout << "Contains tests passed!" << endl;

    double values[3] = {1.0, 2.0, 3.0};
    Vecteur<2> query((double[2]){5., 5.});
	assert(s.interpol(query, values) == 2.5);

    {
        double values[3] = {0.0, 0.0, 0.0};
        Vecteur<2> query((double[2]){0.0, 0.0});
        double result = s.interpol(query, values);
        assert(std::abs(result) < 1e-6);
    }

    std::cout << "Interpol tests passed!" << std::endl;
    return 0;
}
