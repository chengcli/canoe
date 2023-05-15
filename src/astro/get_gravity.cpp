// C/C++ headers
#include <cmath>
#include <cstring>
#include <sstream>
#include <stdexcept>

double GetGravity(char const *planet, double pclat)
{
  std::stringstream msg;
  double g_new;

  if (strcmp(planet, "Jupiter") == 0) {
    double grav = 23.3;
    double S = sin(pclat * M_PI / 180.);
    double SS = S * S;
    double CS = S*sqrt(1 - SS);
    double GR = - grav + SS*(-4.26594 + SS*(0.47685 + SS*(-0.100513 + SS*(0.0237067 - 0.00305515*SS))));
    double GTH = CS*(-3.42313 +  SS*(0.119119 + SS*(0.00533106 + SS*(-0.00647658 + SS*0.000785945))));
    g_new = sqrt(GR*GR+GTH*GTH);
  } else if (strcmp(planet, "Saturn")) {
    double S = sin(pclat * M_PI / 180.);
    double SS = S * S;
    double CS = S*sqrt(1 - SS);
    double GR = -9.06656 + SS*(-3.59253 + SS*(0.704538 + SS*(-0.260158 + SS*(0.0923098 - SS*0.0166287))));
    double GTH = CS*(-2.25384 + SS*(.152112 + SS*(-.0102391 + SS*(-.00714765 + SS*.000865634))));
    g_new = sqrt(GR*GR+GTH*GTH);
  } else {
    msg << "### FATAL ERROR in planet_gravity" << std::endl
        << "Name of the planet not recognized" << std::endl
        << "Choose from [Jupiter|Saturn]" << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }

  return g_new;
}
