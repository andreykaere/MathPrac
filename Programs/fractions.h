#ifndef __fractions__
#define __fractions__

int gcd(int a, int b);

struct Frac
{
	int num, den;
	Frac(int nw_num, int nw_den) : num(nw_num), den(nw_den)
	{}
};

const Frac operator+(const Frac& q1, const Frac& q2

#endif
