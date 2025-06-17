#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include "read_header.hpp"
#include <boost/math/tools/roots.hpp>

namespace fs = std::filesystem;

int main() {
    std::string home = fs::current_path();
    std::string data_path =  home + "/data/LowRes/plt01000"; 
    std::string header = data_path + "/Header";

    std::cout << "Loading data from path: " << data_path<< "\n";
    HeaderInfo hinfo = read_quokka_header(header);
    std::cout << "Header information loaded successfully.\n";
    std::string level_path = data_path  + "/Level_0";
    std::vector<std::string> cell_files = get_all_cell_files(level_path);
    std::cout<< "Cell files loaded successfully.\n\n";

    // i-th entry of all_block_variables will hold all blocks for the i-th variable, j-th entry of the i-th all_block would hold all blocks for the j-th cell file
    // i=0: density, i=1: momentum_x, i=2: momentum_y, i=3: momentum_z, i=4: internal energy
    std::vector<std::vector<BlockData>> all_blocks_variables(5, std::vector<BlockData>());
    if (cell_files.empty()) {
        std::cerr << "No cell files found in the specified directory." << std::endl<< std::endl;
        return 1;
    }

    for(int i=0; i<5; ++i) { // Read 5 variables: density, momentum_x, momentum_y, momentum_z, internal energy
        if(i>=4) {
            for (const auto& cell_file : cell_files) {
            load_data(cell_file, all_blocks_variables[i], i+1); // Assuming 5 is the index for internal energy density
            }
        }else{
            for (const auto& cell_file : cell_files) {
            load_data(cell_file, all_blocks_variables[i], i);
            }
        }
    }
    std::vector<std::vector<BlockData>> all_blocks_variables_velocity(3, std::vector<BlockData>());

    for (int i=1; i<4; ++i) { // Loop over momentum variables
        for (int j=0; j<all_blocks_variables[i].size(); j++) {
            BlockData velocity_block = all_blocks_variables[i][j] / (all_blocks_variables[0][j]); // velocity = momentum / density (in cm/s)
            all_blocks_variables_velocity[i-1].push_back(velocity_block);
            
        }
    }
    
    // Replace the momentum variables in all_blocks_variables with the velocity blocks
    all_blocks_variables[1] = all_blocks_variables_velocity[0]; // Replace momentum_x with velocity_x
    all_blocks_variables[2] = all_blocks_variables_velocity[1]; // Replace momentum_y with velocity_y
    all_blocks_variables[3] = all_blocks_variables_velocity[2]; // Replace momentum_z with velocity_z
    std::cout << "Momentum blocks converted to velocity in cm/s.\n";


    //Initialize the data array for density and velocity variables
    //Assuming global_nx, global_ny, global_nz are the dimensions of the data
    std::vector<std::vector<std::vector<std::vector<double>>>> phys_var(5, std::vector<std::vector<std::vector<double>>>(hinfo.global_nx, std::vector<std::vector<double>>(hinfo.global_ny, std::vector<double>(hinfo.global_nz, 0.0))));


    // Fill the data arrays with the data from all_blocks_variables
    for (int var = 0; var < all_blocks_variables.size(); ++var) {
        for (int blk = 0; blk < all_blocks_variables[var].size(); ++blk) {
            const BlockData& block = all_blocks_variables[var][blk];
            if (block.start_x < 0 || block.end_x >= hinfo.global_nx ||
                block.start_y < 0 || block.end_y >= hinfo.global_ny ||
                block.start_z < 0 || block.end_z >= hinfo.global_nz) {
                std::cerr << "Block indices out of bounds for global dimensions.\n";
                continue;
            }
            for (int i = block.start_x; i <= block.end_x; ++i) {
                for (int j = block.start_y; j <= block.end_y; ++j) {
                    for (int k = block.start_z; k <= block.end_z; ++k) {
                        phys_var[var][i][j][k] = block.data(i - block.start_x, j - block.start_y, k - block.start_z);
                    }
                }
            }
    }
    }

    // Calculate kinetic energy density in ergs/cm^3
    std::vector<std::vector<std::vector<double>>> kinetic_energy = kinetic_energy_density(phys_var[1], phys_var[2], phys_var[3], phys_var[0]); //input is velocity_x, velocity_y, velocity_z, density
    if (kinetic_energy.empty()) {
        std::cerr << "Failed to calculate kinetic energy density.\n";
        return 1;
    }



    //Calculating the Temperature from the eqn ((gamma-1)*energy density= rho*k_B*T/ mu(n_H,T)*m_p)
    std::vector<double> Temperatures=read_vector_csv(home + "/Temperature.csv");
    if (Temperatures.empty()) {
        std::cerr << "Failed to read Temperature.csv.\n";
        return 1;
    }
    std::cout << "Temperatures loaded successfully.\n";
    std::vector<double> Log_nH=read_vector_csv(home + "/log10_nH.csv");
    if (Log_nH.empty()) {
        std::cerr << "Failed to read log10_nH.csv.\n";
        return 1;
    }
    std::vector<std::vector<double>> Mu_grid_slice= read_matrix_csv(home + "/mu_slice_z0.csv");
    if (Mu_grid_slice.empty()) {
        std::cerr << "Failed to read mu_slice_z0.csv.\n";
        return 1;
    }
    std::cout << "mu_slice_z0 loaded successfully.\n";

    // Define constants
    constexpr double cloudy_H_mass_fraction = 1. / (1. + 0.1 * 3.971);
    constexpr double X = cloudy_H_mass_fraction;
    constexpr double Z = 0.02; // metal fraction by mass
    constexpr double Y = 1. - X - Z;
    constexpr double mean_metals_A = 16.; // mean atomic weight of metals
    constexpr double k_B = 1.380649e-23;  // J/K
    constexpr double m_p = 1.6726219e-27; // kg
    constexpr double electron_mass_cgs = 9.10938356e-28; // g
    constexpr double m_e = electron_mass_cgs; // Electron mass in g
    constexpr double gamma = 5.0 / 3.0;

    auto mu = make_mu_interpolator(Log_nH, Temperatures, Mu_grid_slice);

    auto residual = [=](double rho, double T, double e_int) -> double {

        double n_H = rho * X / m_p;
        double mu_val = mu(n_H, T);

        return ((gamma - 1.0) * e_int * mu_val * m_p) / (rho * k_B * T) - T;
    };

    double T_min = 10.0; // Minimum temperature in K
    double T_max = 1e9; // Maximum temperature in K

    std::vector<std::vector<std::vector<double>>> Temperature(std::vector<std::vector<std::vector<double>>>(hinfo.global_nx, std::vector<std::vector<double>>(hinfo.global_ny, std::vector<double>(hinfo.global_nz, 0.0))));
    std::vector<std::vector<std::vector<double>>> electron_density(hinfo.global_nx, std::vector<std::vector<double>>(hinfo.global_ny, std::vector<double>(hinfo.global_nz, 0.0)));
    std::vector<std::vector<std::vector<double>>> hydrogen_number_density(hinfo.global_nx, std::vector<std::vector<double>>(hinfo.global_ny, std::vector<double>(hinfo.global_nz, 0.0)));
    std::vector<std::vector<std::vector<double>>> helium_number_density(hinfo.global_nx, std::vector<std::vector<double>>(hinfo.global_ny, std::vector<double>(hinfo.global_nz, 0.0)));
    // Assuming a constant ratio of He to H for simplicity
    double He_to_H_ratio = 0.1; // Assuming a constant ratio for simplicity
    std::cout << "Calculating temperature and electron density and densities of Hydrogen and Helium\n";
    for (int i = 0; i < hinfo.global_nx; ++i) {
        for (int j = 0; j < hinfo.global_ny; ++j) {
            for (int k = 0; k < hinfo.global_nz; ++k) {
                double rho = phys_var[0][i][j][k]; // Density
                double e_int = phys_var[4][i][j][k]; // Internal energy density
                if (rho > 0 && e_int > 0) { // Avoid division by zero
                    auto temp_pair = boost::math::tools::bisect(
                        [&](double T) { return residual(rho, T, e_int); },
                        T_min, T_max,
                        boost::math::tools::eps_tolerance<double>(30)
                    );
                    Temperature[i][j][k] = (temp_pair.first + temp_pair.second) / 2.0;
                    double nH = rho * X / m_p; // Number density of hydrogen
                    double mu_val = mu(nH, Temperature[i][j][k]);
                                    // compute electron density
	                // N.B. it is absolutely critical to include the metal contribution here!
	                // the approximation for the metals contribution to e- fails at high densities (~1e3 or higher)
                    double n_e = (rho / (m_p + m_e)) * (1.0 - mu_val * (X + Y / 4. + Z / mean_metals_A)) / (mu_val - (electron_mass_cgs / (m_p + m_e)));
                    if (n_e > 1.0e-4 * nH) {
                        electron_density[i][j][k] = 1.0e-4 * nH; // Set a minimum electron density
                    }else if (n_e >=0.0) {
                        electron_density[i][j][k] = n_e; // Set the calculated electron density
                    } else {
                        std::cerr << "Negative electron density calculated at (" << i << ", " << j << ", " << k << "). Setting to zero.\n";
                        electron_density[i][j][k] = 0.0; // Set to zero or handle as needed
                    }
                    
                } else {
                    std::cerr << "Invalid density or internal energy at (" << i << ", " << j << ", " << k << "). Setting temperature to zero.\n";
                    Temperature[i][j][k] = 0.0; // Set to zero or handle as needed
                    electron_density[i][j][k] = 0.0; // Set to zero or handle as needed
                }
                hydrogen_number_density[i][j][k] = rho * X / m_p; // Density * X / m_p
                helium_number_density[i][j][k] = rho * X * He_to_H_ratio / m_p; // n_He/n_H=0.1

            }
        }
    }
    phys_var.push_back(Temperature); // Add the Temperature array to phys_var for further processing if needed
    phys_var.push_back(electron_density); // Add the electron density array to phys_var for further processing if needed
    phys_var.push_back(hydrogen_number_density); // Add the hydrogen number density array to phys_var for further processing if needed
    phys_var.push_back(helium_number_density); // Add the helium number density array to phys_var for further processing if needed
    std::cout << "Temperature calculated successfully and added to phys_var at index " << phys_var.size() - 2 << "\n";
    std::cout << "Electron density calculated successfully and added to phys_var at index " << phys_var.size() - 1 << "\n";
    std::cout << "Hydrogen number density calculated successfully and added to phys_var at index " << phys_var.size() << "\n";
    std::cout << "Helium number density calculated successfully and added to phys_var at index " << phys_var.size() + 1 << "\n";
    

    // Calculate the density of ions: n_H+, n_He+, and n_He++
    std::vector<std::vector<std::vector<double>>> H_plus_density(hinfo.global_nx, std::vector<std::vector<double>>(hinfo.global_ny, std::vector<double>(hinfo.global_nz, 0.0)));
    std::vector<std::vector<std::vector<double>>> He_plus_density(hinfo.global_nx, std::vector<std::vector<double>>(hinfo.global_ny, std::vector<double>(hinfo.global_nz, 0.0)));
    std::vector<std::vector<std::vector<double>>> He_double_plus_density(hinfo.global_nx, std::vector<std::vector<double>>(hinfo.global_ny, std::vector<double>(hinfo.global_nz, 0.0)));\
    std::vector<std::vector<std::vector<double>>> Neutral_hydrogen_density(hinfo.global_nx, std::vector<std::vector<double>>(hinfo.global_ny, std::vector<double>(hinfo.global_nz, 0.0)));
    std::vector<std::vector<std::vector<double>>> Neutral_helium_density(hinfo.global_nx, std::vector<std::vector<double>>(hinfo.global_ny, std::vector<double>(hinfo.global_nz, 0.0)));
    std::cout << "Calculating ion densities...\n";
    for (int i = 0; i < hinfo.global_nx; ++i) {
        for (int j = 0; j < hinfo.global_ny; ++j) {
            for (int k = 0; k < hinfo.global_nz; ++k) {
                double n_e = electron_density[i][j][k];
                double n_H = hydrogen_number_density[i][j][k];
                double n_He = helium_number_density[i][j][k];
                if(n_e < (n_H * (1 + He_to_H_ratio))) {
                    // Assuming the ionization fraction, chi is same for both hydrogen and helium, we have n_e=n_H+ + n_He+  implying n_H+=n_H*n_e/(n_H+n_He)
                    H_plus_density[i][j][k] = n_H * n_e / (n_H + n_He);
                    He_plus_density[i][j][k] = He_to_H_ratio * n_H * n_e / (n_H + n_He);
                    He_double_plus_density[i][j][k] = 0.0; // Assuming no He++ in this case
                    Neutral_hydrogen_density[i][j][k] = n_H - H_plus_density[i][j][k];
                    Neutral_helium_density[i][j][k] = n_He - He_plus_density[i][j][k];
                } else {   //n_e= n_H + n_He + n_He++
                    He_double_plus_density[i][j][k] = n_e - n_H - n_He; // Assuming all excess electrons go to He++
                    Neutral_hydrogen_density[i][j][k] = 0.0; // Assuming all hydrogen is ionized
                    Neutral_helium_density[i][j][k] = 0.0; // Assuming all helium is ionized
                    H_plus_density[i][j][k] = n_H;
                    He_plus_density[i][j][k] = n_He - (n_e - n_H - n_He); // Remaining He+ after accounting for He++
                }
            }
        }
    }

    std::cout << "Ion densities calculated successfully.\n";
    // Add the ion densities to phys_var for further processing if needed
    phys_var.push_back(Neutral_hydrogen_density);
    phys_var.push_back(Neutral_helium_density);
    phys_var.push_back(H_plus_density);
    phys_var.push_back(He_plus_density);
    phys_var.push_back(He_double_plus_density);

    //print the indices of the variables in phys_var
    std::cout << "Physical variables indices in phys_var:\n";
    std::cout << "0: Density (g/cm^3)\n";
    std::cout << "1: Velocity_x (cm/s)\n";
    std::cout << "2: Velocity_y (cm/s)\n";
    std::cout << "3: Velocity_z (cm/s)\n";
    std::cout << "4: Internal Energy Density (ergs/cm^3)\n";
    std::cout << "5: Temperature (K)\n";
    std::cout << "6: Electron Density (cm^-3)\n";
    std::cout << "7: Neutral Hydrogen Density (cm^-3)\n";
    std::cout << "8: Neutral Helium Density (cm^-3)\n";
    std::cout << "9: H+ Density (cm^-3)\n";
    std::cout << "10: He+ Density (cm^-3)\n";
    std::cout << "11: He++ Density (cm^-3)\n";
    
    // Adding the magnetic field variable.
    // Zeroth order approximation of the magnetic field: B= B_0 \hat{z}
    double B_0 = 1e-6; // Example value for the magnetic field strength in Gauss
    std::vector<std::vector<std::vector<double>>> magnetic_field_x(hinfo.global_nx, std::vector<std::vector<double>>(hinfo.global_ny, std::vector<double>(hinfo.global_nz, 0.0)));
    std::vector<std::vector<std::vector<double>>> magnetic_field_y(hinfo.global_nx, std::vector<std::vector<double>>(hinfo.global_ny, std::vector<double>(hinfo.global_nz, 0.0)));
    std::vector<std::vector<std::vector<double>>> magnetic_field_z(hinfo.global_nx, std::vector<std::vector<double>>(hinfo.global_ny, std::vector<double>(hinfo.global_nz, B_0)));
    phys_var.push_back(magnetic_field_x); // Add magnetic field x-component
    phys_var.push_back(magnetic_field_y); // Add magnetic field y-component
    phys_var.push_back(magnetic_field_z); // Add magnetic field z-component
    std::cout << "Magnetic field added to phys_var at indices " << phys_var.size() - 3 << ", " << phys_var.size() - 2 << ", " << phys_var.size() - 1 << "\n";
    // Validate the data
    std::vector<int> indices = {0, 55, 62}; // Example indices to validate
    //last argument is true, which means we want to validate the data along with printing the values at the specified indices. Turn it to false if you only want to validate the data without printing the values.
    std::cout << "Validating data...\n\n";
    validate(phys_var, hinfo, indices, true);
    std::cout<< "Total kinetic energy in the domain is " << total_amount(kinetic_energy, hinfo) << " ergs.\n";
    std::cout << "Data validation completed.\n";


    return 0;
}

