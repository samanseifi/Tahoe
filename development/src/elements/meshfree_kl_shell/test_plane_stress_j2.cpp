/* Standalone regression tests for PlaneStressJ2.h. The test links Tahoe's
 * SecantMethodT and ExceptionT implementations; use the normal project include
 * paths (or all toolbox/src subdirectories for an ad-hoc build). */
#include <cmath>
#include <cstdio>
#include "PlaneStressJ2.h"

using namespace Tahoe::KLShell;

static int failures = 0;
static void check(const char* name, double got, double expected, double tolerance)
{
	bool pass = std::fabs(got - expected) <= tolerance;
	std::printf("  [%s] %-34s got=% .12e expected=% .12e\n",
		pass ? "PASS" : "FAIL", name, got, expected);
	if (!pass) ++failures;
}

int main()
{
	/* Eq. 79 is exact at zero, reciprocal for reversed increments, and agrees
	 * with the exponential update through second order. */
	check("Pade(0)", ThicknessStretchPade(0.0), 1.0, 0.0);
	check("Pade reciprocal", ThicknessStretchPade(0.2)*ThicknessStretchPade(-0.2), 1.0, 1.0e-14);
	check("Pade second-order accuracy", ThicknessStretchPade(1.0e-3), std::exp(1.0e-3), 1.0e-9);
	check("Pade rejects pole", ThicknessStretchPade(2.0), 0.0, 0.0);

	/* An elastic uniaxial-strain increment must reproduce the condensed
	 * plane-stress law exactly and leave the plastic history untouched. */
	double sig[3] = {0.0, 0.0, 0.0};
	double deps[3] = {1.0e-4, 0.0, 0.0};
	double ep = 0.0;
	const double E = 210000.0, nu = 0.3;
	PlaneStressJ2Return(sig, deps, ep, E, nu, 1.0e9, 0.0);
	double c = E/(1.0 - nu*nu);
	check("elastic sigma11", sig[0], c*deps[0], 1.0e-10);
	check("elastic sigma22", sig[1], c*nu*deps[0], 1.0e-10);
	check("elastic plastic strain", ep, 0.0, 0.0);

	std::printf("\n%s (%d failures)\n", failures ? "FAILURES" : "ALL PASS", failures);
	return failures;
}
