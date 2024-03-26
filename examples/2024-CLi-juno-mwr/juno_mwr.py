from canoe import def_species, load_configure,index_map, air_parcel, VariableType
from canoe.athena import mesh_block
from canoe.snap import init_thermo, thermodynamics
# from canoe.harp import radiation_band
import numpy as np


# def  atm_extrapolate(
#         air: air_parcel, 
#         dz:float,
#         method="dry_adiabat"):
#     pass

def construct_atmosphere(
                        nlyr: int = 100,
                        comp: dict = {},
                        T0: float = 169.0,
                        P0: float = 1.0e5,                           
                        plim: list[float, float] = [0.1e5, 100.0e5],        
                        Tmin: float = 100.0,
                        ) -> mesh_block:
    
    # construct a airparcel obj
    var_type=VariableType(2) ## molefraction
    air = air_parcel(var_type)



    # exit()
    # get IndexMap
    pindex=index_map.get_instance()
    # iH2O=pindex.get_vapor_id("H2O")
    # iNH3=pindex.get_vapor_id("NH3")
    # print(iH2O)
    # print(iNH3)
    # exit()
    # air.set_vapor(iH2O,30)


    # set property for air cube
    print("Load vapors:")
    for vapor in comp:
        ivapor=pindex.get_vapor_id(vapor)
        air.set_vapor(ivapor,comp[vapor]*1E-6) ## molefrac
        # air.hydro()[ivapor]=comp[vapor]*1E-6  ## molefrac

        print(f"Loaded: {vapor}-- vaporID: {ivapor} -- molefrac: {comp[vapor]*1E-6}")

    print(air.hydro())
    print(air.cloud())
    print(air.tracer())
    
    exit()
    # vertical grid spacing
    dlnp = (np.log(plim[1]) - np.log(plim[0])) / nlyr
    print(dlnp)



    # construct a meshblock instance
    mb = mesh_block(nlyr)


    mb.set_layer(0, air)
    for i in range(1, nlyr):
        # atm_extrapolate(air, -dlnp / 2.0, method="dry_adiabat")
        mb.set_layer(i, air)

    return mb


# def modify_atmosphere(
#     mb: mesh_block,
#     comp: dict = {},
#     plim: List[float, float] = [0.1e5, 100.0e5],
#     zlim: List[float, float] = [0.0, 10.0e3],
#     T0: float = 169.0,
#     P0: float = 1.0e5,
# ) -> mesh_block:
#     pass


if __name__ == "__main__":

    # load config
    config = load_configure("jupiter_atmos.yaml")

    # enroll species
    # def_species(config['species'])
    def_species(vapors=["H2O", "NH3"], clouds=["H2O(s)", "H2O(l)", "NH3(s)", "NH3(l)"],
                tracers=["e-", "Na"] )
    # get IndexMap
    # pindex=index_map.get_instance()
    # print(pindex.get_vapor_id("H2O"))
    # print(pindex.get_cloud_id("H2O(s)"))

    # print(pindex.get_species_id("vapor.H2O"))
    # print(pindex.get_species_id("cloud.H2O(s)"))
    # print(pindex.get_species_id("tracer.e-"))
    # iH2O=pindex.get_vapor_id("H2O")

    # init thermodynamics
    init_thermo(config["thermodynamics"])

    # get pthermo
    pthermo=thermodynamics.get_instance()

    # print(pthermo.get_Rd())
    # print(pthermo.get_Gammad())
    # print(pthermo.get_GammadRef())

    nlyr = 100


    
    # # pressure scale height
    pmax = 100.0e5
    pmin = 0.1e5

    P0 = 1.0e5
    T0 = 169.0

    vapors = {}
    vapors["NH3"] = 300# ppmv
    vapors["H2O"] = 3000# ppmv

    mb = construct_atmosphere(100, vapors, T0, P0, plim=[pmin, pmax])

    # mb.add_radiation(config)
    # rad = mb.cal_radiance()

    # mb = modify_atmosphere(mb, T0, P0, plim=[pmin, pmax])
    # rad = mb.cal_radiance()
