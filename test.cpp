#include <iostream>
#include <vector>
using namespace std;

// void function(int **pointers, int *values, int N)
// {
//     for (int i = 0; i < N; i++)
//         pointers[i] = &values[i % 2];
//     values[1] = 10;
// }

void function_normal(int x1, int x2, int x3, int x4)
{
    cout << x1 << ' ' << x2 << ' ' << x3 << ' ' << x4 << '\n';
}

void function_guai(int *x1, int *x3)
{
    cout << *x1++ << ' ' << *x1 << ' ' << *x3++ << ' ' << *x3 << '\n';
}

int main(int argc, char *argv[])
{

    // int N = atoi(argv[1]);

    int a[] = {00, 10, 20, 30, 40, 50, 60, 70, 80, 90};

    function_normal(a[1], a[2], a[6], a[7]);
    function_guai(&a[1], &a[6]);
}
