#include <iostream>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "fraction.h"
#include "Lagranzh.h"
#include "Lezhandr.h"

 
 
int main()
{
    char per[3] = {'x', 'y', 'z'};
    int n = 3;
    Frac a[9];
    int val, A, B, C;
    printf("Введите квадратичную форму:\n");
	
	for(int i = 0; i < n; ++i)
		for(int j = 0; j < n; ++j)
		{
			scanf("%d", &val);
			a[i*n + j] = Frac(val, 1);
		}
    
    triple out = Lagranzh(a);

    int DEN = out.x.den * out.y.den * out.z.den;
    trinity coef;
    coef.x = out.x.num * DEN;
    coef.y = out.y.num * DEN;
    coef.z = out.z.num * DEN;
    
    int g = gcd(coef.x, coef.y);
    g = gcd(g, coef.z);
    coef.x = coef.x / g;
    coef.y = coef.y / g;
    coef.z = coef.z / g;
    
    A = coef.x;
    B = coef.y;
    C = coef.z;
    
    
    
    printf("Итак, получили уравнение:\n");
    printf("%dx^2 + %dy^2 + %dz^2 = 0\n", A, B, C);
    
    if(A * B * C > 0)
    {
        if(A > 0)
        {
            if(B > 0)
            {
                printf("Решений нет\n");
                return 0;
            }
            coef.x = A;
            coef.y = C;
            coef.z = -B;
        }
        else
        {
            if(B > 0)
            {
                coef.x = -A;
                coef.y = -C;
                coef.z = B;
            }
            coef.x = -A;
            coef.y = -B;
            coef.z = C;
        }
    }
    else
    {
        A = -A;
        B = -B;
        C = -C;
        if(A > 0)
        {
            if(B > 0)
            {
                printf("Решений нет\n");
                return 0;
            }
            coef.x = A;
            coef.y = C;
            coef.z = -B;
        }
        else
        {
            if(B > 0)
            {
                coef.x = -A;
                coef.y = -C;
                coef.z = B;
            }
            coef.x = -A;
            coef.y = -B;
            coef.z = C;
        }
    }
    
    printf("Переобозначив переменные, получим уравнение:\n");
    printf("%dx^2+%dy^2 = %dz^2\n", coef.x, coef.y, coef.z);
    
    Lezhandr(coef);
    
    return 0;
}
