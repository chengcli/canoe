#ifndef NDARRAYS_HPP
#define NDARRAYS_HPP

template<typename T>
void NewCArray(T** &a, int n1, int n2)
{
  a = new T* [n1];
  a[0] = new T [n1*n2];

  for (int i = 0; i < n1; ++i)
    a[i] = a[0] + i*n2;
}

template<typename T>
void FreeCArray(T **a)
{
  delete[] a[0];
  delete[] a;
}

template<typename T>
void FreeCArray2(T **a)
{
  delete[] a[0];
  delete[] a;
}

template<typename T>
void NewCArray(T*** &a, int n1, int n2, int n3)
{
  a = new T** [n1];
  a[0] = new T* [n1*n2];
  a[0][0] = new T [n1*n2*n3];

  for (int i = 0; i < n1; ++i) {
    a[i] = a[0] + i*n2;
    for (int j = 0; j < n2; ++j)
      a[i][j] = a[0][0] + i*n2*n3 + j*n3;
  }
}

template<typename T>
void FreeCArray(T ***a)
{
  delete[] a[0][0];
  delete[] a[0];
  delete[] a;
}

template<typename T>
void FreeCArray3(T ***a)
{
  delete[] a[0][0];
  delete[] a[0];
  delete[] a;
}

template<typename T>
void NewCArray(T**** &a, int n1, int n2, int n3, int n4)
{
  a = new T*** [n1];
  a[0] = new T** [n1*n2];
  a[0][0] = new T* [n1*n2*n3];
  a[0][0][0] = new T [n1*n2*n3*n4];

  for (int i = 0; i < n1; ++i) {
    a[i] = a[0] + i*n2;
    for (int j = 0; j < n2; ++j) {
      a[i][j] = a[0][0] + i*n2*n3 + j*n3;
      for (int k = 0; k < n3; ++k)
        a[i][j][k] = a[0][0][0] + i*n2*n3*n4 + j*n3*n4 + k*n4;
    }
  }
}

template<typename T>
void FreeCArray(T ****a)
{
  delete[] a[0][0][0];
  delete[] a[0][0];
  delete[] a[0];
  delete[] a;
}

template<typename T>
void FreeCArray4(T ****a)
{
  delete[] a[0][0][0];
  delete[] a[0][0];
  delete[] a[0];
  delete[] a;
}

#endif
