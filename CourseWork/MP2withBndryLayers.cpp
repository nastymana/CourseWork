#include "pch.h"
#include "MP2withBndryLayers.h"


MP2withBndryLayers::MP2withBndryLayers()
{
}

MP2withBndryLayers::MP2withBndryLayers(const vector<double> &inData, const vector<double>& domainParams)
{
	cout << "\n MP2withBndryLayers::Constructor()\n";
	lyambdaParam = domainParams[0];
	gammaParam = domainParams[1];
	velocityParam1 = domainParams[2];
	velocityParam2 = domainParams[3];
	sourceParam = domainParams[4];
}

MP2withBndryLayers::~MP2withBndryLayers()
{
}
double MP2withBndryLayers::realSol(const double & x, const double & y){
	return 	((1 - exp(-(1 - x) / lymbda(x, y))) * (1 - exp(-(1 - y) / lymbda(x, y))) * cos(PI*(x + y)));
}

double MP2withBndryLayers::BC1(const double & x, const double & y) {
	if (((x == 0.) || (x == 1.)) || ((y == 0.) || (y == 1.))) return 0.0;
}

double MP2withBndryLayers::BC2(const double & x, const double & y) {
	return 0.0;
}

double MP2withBndryLayers::BC3(const double & x, const double & y) {
	return 0.0;
}
double MP2withBndryLayers::source(const double & x, const double & y) {
	if (((x >= 0.) && (x <= 1.)) && ((y >= 0.) && (y <= 1.))) {
		vector<double> v = velocity(x, y);

		double expx = exp((x - 1.0) / lymbda(x,y)),
			expy = exp((y - 1) / lymbda(x, y)),
			cosxy = cos(PI*(x + x)),
			sinxy = sin(PI*(x + y));
		double divgradX = (1.0 - expy)*(cosxy*(expx*(PI*PI*lymbda(x, y) - (1.0 / lymbda(x, y)))) + sinxy * 2.0 * PI*expx);
		double divgradY = (1.0 - expx)*(cosxy*(expy*(PI*PI*lymbda(x, y) - (1.0 / lymbda(x, y)))) + sinxy * 2.0 * PI*expy);
		double ugradX = v[0] * (1.0 - expy)*(-(cosxy*expx / lymbda(x, y)) - PI * (1 - expx)*sinxy),
			ugradY = v[1] * (1.0 - expx)*(-(cosxy*expy / lymbda(x, y)) - PI * (1 - expy)*sinxy);;
		return(-divgradX - divgradY + ugradX + ugradY);

	}
}

double MP2withBndryLayers::lymbda(const double & x, const double & y) {

	return lyambdaParam;
}

double MP2withBndryLayers::gamma(const double & x, const double & y) {
	return gammaParam;
}

vector<double> MP2withBndryLayers::velocity(const double & x, const double & y) {
	return { velocityParam1, velocityParam2 };
}

double MP2withBndryLayers::initialCond(const double & x, const double & y)
{
	return 0.0;
}
