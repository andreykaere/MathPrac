#include <iostream>
#include <stdio.h>
#include <math.h>
#include "fraction.h"

    
int gcd(int a, int b)
{
	if(a == 0)
		return abs(b);
	return gcd(abs(b) % abs(a), abs(a));
}

const Frac operator*(const int &l, const Frac &q)
{
	int a = l*q.num;
	int b = q.den;
	int c = gcd(a, b);
	return Frac(a / c, b / c);
}

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

const Frac operator*(const Frac& q1, const Frac& q2)
{
	int a = q1.num * q2.num;
	int b = q1.den * q2.den;
	int c = gcd(a, b);
	return Frac(a / c, b / c);
}

const Frac operator/(const Frac& q1, const Frac& q2)
{
	int a = q2.den*q1.num;
	int b = q2.num*q1.den;
	int c = gcd(a, b);
	if(b == 0)
	{
		printf("Ты чё на 0 делишь, чепух?\n");
		return Frac(0, 1);
	}
	if(b < 0)
	{
		a = -a;
		b = -b;
	}
	return Frac(a / c, b / c);
}	
