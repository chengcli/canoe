// Athena++ headers
#include <athena.hpp>
#undef PI
#undef SQRT2

// IceChem headers
#include <transport/stream_transport.hpp>

void set_boundary_condition(AthenaArray<Real> &data)
{
  int rowsh = data.GetDim2();
  int colsh = data.GetDim1();
  int rows = rowsh - 2*GhostZoneSize;
  int cols = colsh - 2*GhostZoneSize;

  int is = GhostZoneSize, ie = cols + GhostZoneSize - 1;
  int js = GhostZoneSize, je = rows + GhostZoneSize - 1;

  for (int j = 0; j < rowsh; ++j) {
    for (int i = 1; i <= GhostZoneSize; ++i)
      data(j,is-i) = data(j,is);
    for (int i = 1; i <= GhostZoneSize; ++i)
      data(j,ie+i) = data(j,ie);
  }
  for (int i = 0; i < colsh; ++i) {
    for (int j = 1; j <= GhostZoneSize; ++j)
      data(js-j,i) = data(js,i);
    for (int j = 1; j <= GhostZoneSize; ++j)
      data(je+j,i) = data(je,i);
  }
  /*for (int i = 0; i < colsh; ++i) {
    data(0,i) = data(1,i);
    data(rowsh-1,i) = data(rows-1,i);
  }*/
  // neumann boundary matrix set in src/transport/stream_transport.cpp
  for (int j = 0; j < rowsh; ++j) {
    data(j,0) = 0.;
    data(j,colsh-1) = 1.;
  }
}

int main(int argc, char **argv) {
  int rows = 5, cols = 5, ntimes = 5;

  // read command line options
  if (argc > 3) {
      rows = std::atoi(argv[1]);
      cols = std::atoi(argv[2]);
      ntimes  = atoi(argv[3]);
  } else if (argc > 2) {
      rows = std::atoi(argv[1]);
      cols = std::atoi(argv[2]);
  }

  int rowsh = rows + 2*GhostZoneSize;
  int colsh = cols + 2*GhostZoneSize;

  Real dx = 1.;
  Real keddy = 0.1;
  //Real keddy = 0.0;
  bool fourth_order = true;

  StreamTransport model(rows, cols, fourth_order);
  AthenaArray<Real> data(rowsh, colsh);
  AthenaArray<Real> streamf(rowsh, colsh);
  data.ZeroClear();
  streamf.ZeroClear();

  int is = GhostZoneSize, ie = cols + GhostZoneSize - 1;
  int js = GhostZoneSize, je = rows + GhostZoneSize - 1;

  // interior
  for (int j = 0; j < rowsh; ++j)
    for (int i = 0; i < colsh; ++i) {
      streamf(j,i) = - (j + 1)*dx;
      //streamf(j,i) = 0.;
      if (j > rowsh/4. && j < 3.*rowsh/4. && i > colsh/4. && i < 3.*colsh/4.)
        data(j,i) = 1.;
    }

  // boundary
  set_boundary_condition(data);
  //for (int i = 0; i < colsh; ++i)
  //  data(rowsh-1,i) = 1.;

  // initial print
  Real sum = 0.;
  printf("data = \n");
  for (int j = 0; j < rowsh; ++j) {
    for (int i = 0; i < colsh; ++i)
      printf("%7.2f", data(j,i));
    printf("\n");
  }
  for (int j = js; j <= je; ++j)
    for (int i = is; i <= ie; ++i)
      sum += data(j,i);
  printf("sum = %7.2f\n", sum);
  printf("streamf = \n");
  for (int j = 0; j < rowsh; ++j) {
    for (int i = 0; i < colsh; ++i)
      printf("%7.2f", streamf(j,i));
    printf("\n");
  }
  printf("\n");

  model.setAdvectionMatrix(streamf, dx);
  model.setDiffusionMatrix(keddy, dx);

  for (int n = 0; n < ntimes; ++n) {
    model.evolve(data, 1.);
    set_boundary_condition(data);

    Real sum = 0.;
    printf("data = \n");
    for (int j = 0; j < rowsh; ++j) {
      for (int i = 0; i < colsh; ++i)
        printf("%7.2f", data(j,i));
      printf("\n");
    }
    for (int j = js; j <= je; ++j)
      for (int i = is; i <= ie; ++i)
        sum += data(j,i);
    printf("sum = %7.2f\n", sum);
    printf("streamf = \n");
    for (int j = 0; j < rowsh; ++j) {
      for (int i = 0; i < colsh; ++i)
        printf("%7.2f", streamf(j,i));
      printf("\n");
    }
    printf("\n");
  }
}
