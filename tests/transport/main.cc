#include <cstdlib>
#include <deal.II/base/logstream.h>
#include <transport/Kinetics.h>

using namespace std;

int main(int argc, char *argv[]) {
    dealii::deallog.depth_console(0);
    int nrows = 5, ncols = 5, ntimes = 5;

    // read command line options
    if (argc > 3) {
        nrows = atoi(argv[1]);
        ncols = atoi(argv[2]);
        ntimes  = atoi(argv[3]);
    } else if (argc > 2) {
        nrows = atoi(argv[1]);
        ncols = atoi(argv[2]);
    }

    Kinetics model(nrows, ncols);
    model.assemble(1.0, 0.01);
    //model.checkout();
    model.run(ntimes);
    model.assemble(1.0, 0.1);
    model.run(ntimes);
    model.assemble(1.0, 1);
    model.run(ntimes);
    /*
    RectGrid<double, 2> grid(5, 5, 1);
    Variable<double, 2> psi(grid, "psi", "streamfunction");
    Advection<double, 2, 4> adv(grid, psi);
    Diffusion<double, 2, 4> diff(grid);
    psi.set_boundary(3, Dirichlet, 5);
    psi.finish();
    psi.setZero();
    psi.val().setRandom();
    cout << psi << endl;
    cout << psi.m_dirichlet << endl;
    adv.update();
    diff.update();
    adv().print_formatted(cout, 2, false);
    diff().print_formatted(cout, 2, false);
    */
}
