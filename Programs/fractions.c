#include "fractions.h"

int gcd(int a, int b)
{
	if(a == 0)
		return b;
	return gcd(b % a, a);
}

struct Frac
{
	int num, den;
	Frac(int nw_num, int nw_den) : num(nw_num), den(nw_den)
	{}
};

const Frac operator+(const Frac& q1, const Frac& q2)
{
	int a = q1.num * q2.den + q2.num * q1.den;
	int b = q1.den * q2.den;
	int c = gcd(a, b);
	return Frac(a / c, b / c);
}

const Frac operator-(const Frac& q1)
{
	return Frac(-q1.num, q1.den);
}

const Frac operator-(const Frac& q1, const Frac& q2)
{
	Frac q = -q2;
	return q1 + q;
}

