#ifndef __lezhandr__
#define __lezhandr

struct trinity
{
    int x, y, z;
};

bool is_square(int x, int a);
bool CL(int a, int b, int c);
trinity find_solution(int a, int b, int c);
trinity Lezhandr(trinity input);


#endif