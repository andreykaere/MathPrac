#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdbool.h>
#include <iostream>
#include "fraction.h"
#include "Lezhandr.h"

bool is_square(int x, int a)
{
	for(int i = 0; i < a; ++i)
		if((i*i - x) % a == 0)
			return true;

	return false;
}


bool CL(int a, int b, int c)
{
    bool p = ((gcd(a, b) == 1) && (gcd(a, c) == 1) && (gcd(b, c) == 1));
	return p && (is_square(b*c, a) && is_square(a*c, b) && is_square(-a*b, c));
}

trinity find_solution(int a, int b, int c)
{
	trinity ans;
	int xmax = b*c, ymax = a*c, zmax = a*b;
	
    for(int x0 = 0; x0*x0 <= xmax; ++x0)
		for(int y0 = 0; y0*y0 <= ymax; ++y0)
			for(int z0 = 0; z0*z0 <= zmax; ++z0)
				if( (a*x0*x0 + b*y0*y0 - c*z0*z0 == 0) && (x0 || y0 || z0) )
				{
					ans.x = x0;
					ans.y = y0;
					ans.z = z0;
					return ans;
				}
	
	
	return ans;
}


trinity Lezhandr(trinity input)
{
	double start, end;
	int a = input.x, b = input.y, c = input.z;
    trinity ans;
    ans.x = 0;
    ans.y = 0;
    ans.z = 0;
    if((a == 0) || (b == 0) || (c == 0))
        printf("Вырожденный случай, исследуется просто\n");

    for(int i = 2; i*i < a; ++i)
        if(a % (i*i) == 0)
            a = a / (i*i);
    
    for(int i = 2; i*i < a; ++i)
        if(b % (i*i) == 0)
            b = b / (i*i);
    
    for(int i = 2; i*i < a; ++i)
        if(c % (i*i) == 0)
            c = c / (i*i);
    
    
	if(!(CL(a, b, c)))
	{
		printf("Решений нет.\n");
        ans.x = 0;
        ans.y = 0;
        ans.z = 0;
        return ans;
	}
	
	ans = find_solution(a, b, c);
	printf("Найдено решение (%d, %d, %d).\n", ans.x, ans.y, ans.z);

	return ans;
}

