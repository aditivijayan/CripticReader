#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <stdio.h>
#include <filesystem>
#include "read_header.hpp"
#include <boost/math/tools/roots.hpp>

namespace fs = std::filesystem;


int main() {
    std::string home = fs::current_path();
    std::string file_path = "/data/plt80000/";
    std::string data_path =  home + file_path; 
    std::string header = data_path + "/Header";
    HeaderInfo hinfo = read_quokka_header(header);

    std::string level_path = data_path  + "/Level_0";
    std::vector<std::string> cell_files = get_all_cell_files(level_path);


    // i-th entry of all_block_variables will hold all blocks for the i-th variable, j-th entry of the i-th all_block would hold all blocks for the j-th cell file
    // i=0: density, i=1: momentum_x, i=2: momentum_y, i=3: momentum_z, i=4: internal energy
    std::vector<std::vector<BlockData>> all_blocks_variables(5, std::vector<BlockData>());
    if (cell_files.empty()) {
        std::cerr << "No cell files found in the specified directory." << std::endl<< std::endl;
        return {};
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



    //Calculating the Temperature from the eqn ((gamma-1)*energy density= rho*k_B*T/ mu(n_H,T)*m_p)
    std::vector<double> Temperatures=read_vector_csv(home + "/Temperature.csv");
    if (Temperatures.empty()) {
        std::cerr << "Failed to read Temperature.csv.\n";
        return {};
    }
    std::vector<double> Log_nH=read_vector_csv(home + "/log10_nH.csv");
    if (Log_nH.empty()) {
        std::cerr << "Failed to read log10_nH.csv.\n";
        return {};
    }
    std::vector<std::vector<double>> Mu_grid_slice= read_matrix_csv(home + "/mu_slice_z0.csv");
    if (Mu_grid_slice.empty()) {
        std::cerr << "Failed to read mu_slice_z0.csv.\n";
        return {};
    }
    std::cout << "mu_slice_z0 loaded successfully.\n";

    auto mu = make_mu_interpolator(Log_nH, Temperatures, Mu_grid_slice);

    auto residual = [=](double rho, double T, double e_int) -> double {

        double n_H = rho * X / m_p;
        double mu_val = mu(n_H, T);

        return ((2.0/3.0) * e_int * mu_val * m_p) / (rho * k_B * T) - T;
    };

    std::vector<std::vector<std::vector<double>>> Temperature(std::vector<std::vector<std::vector<double>>>(hinfo.global_nx, std::vector<std::vector<double>>(hinfo.global_ny, std::vector<double>(hinfo.global_nz, 0.0))));
    std::vector<std::vector<std::vector<double>>> relative_electron_density(hinfo.global_nx, std::vector<std::vector<double>>(hinfo.global_ny, std::vector<double>(hinfo.global_nz, 0.0)));
    std::vector<std::vector<std::vector<double>>> hydrogen_number_density(hinfo.global_nx, std::vector<std::vector<double>>(hinfo.global_ny, std::vector<double>(hinfo.global_nz, 0.0)));
    std::vector<std::vector<std::vector<double>>> relative_helium_number_density(hinfo.global_nx, std::vector<std::vector<double>>(hinfo.global_ny, std::vector<double>(hinfo.global_nz, He_to_H_ratio)));
        // Calculate the density of ions: n_H+, n_He+, and n_He++
    std::vector<std::vector<std::vector<double>>> relative_H_plus_density(hinfo.global_nx, std::vector<std::vector<double>>(hinfo.global_ny, std::vector<double>(hinfo.global_nz, 0.0)));
    std::vector<std::vector<std::vector<double>>> relative_He_plus_density(hinfo.global_nx, std::vector<std::vector<double>>(hinfo.global_ny, std::vector<double>(hinfo.global_nz, 0.0)));
    std::vector<std::vector<std::vector<double>>> relative_He_double_plus_density(hinfo.global_nx, std::vector<std::vector<double>>(hinfo.global_ny, std::vector<double>(hinfo.global_nz, 0.0)));
    std::vector<std::vector<std::vector<double>>> relative_Neutral_hydrogen_density(hinfo.global_nx, std::vector<std::vector<double>>(hinfo.global_ny, std::vector<double>(hinfo.global_nz, 0.0)));
    std::vector<std::vector<std::vector<double>>> relative_Neutral_helium_density(hinfo.global_nx, std::vector<std::vector<double>>(hinfo.global_ny, std::vector<double>(hinfo.global_nz, 0.0)));
    // Assuming a constant ratio of He to H for simplicity

    for (int i = 0; i < hinfo.global_nx; ++i) {
        for (int j = 0; j < hinfo.global_ny; ++j) {
            for (int k = 0; k < hinfo.global_nz; ++k) {
                double rho = phys_var[0][i][j][k]; // Density
                double e_int = phys_var[4][i][j][k]; // Internal energy density
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
                if (n_e < 1.0e-4 * nH) {
                        relative_electron_density[i][j][k] = 1.0e-4; // Set a minimum electron density
                        n_e = 1.0e-4 * nH; // Set a minimum electron density
                }else {
                        relative_electron_density[i][j][k] = n_e/nH; // Set the calculated electron density
                }

                hydrogen_number_density[i][j][k] = nH; // Density * X / m_p
                double n_He = hydrogen_number_density[i][j][k] * He_to_H_ratio; // Helium number density
                //helium_number_density[i][j][k] = n_He; // n_He/n_H=0.1
                if(n_e < (nH * (1 + He_to_H_ratio))) {
                    // Assuming the ionization fraction, chi is same for both hydrogen and helium, we have n_e=n_H+ + n_He+  implying n_H+=n_H*n_e/(n_H+n_He)
                    relative_H_plus_density[i][j][k] =  n_e / (nH + n_He);
                    relative_He_plus_density[i][j][k] = He_to_H_ratio * n_e / (nH + n_He);
                    relative_He_double_plus_density[i][j][k] = 0.0; // Approximation: no He++ in this case
                    relative_Neutral_hydrogen_density[i][j][k] = 1- relative_H_plus_density[i][j][k];//(nH - H_plus_density[i][j][k])/nH;
                    relative_Neutral_helium_density[i][j][k] = He_to_H_ratio - relative_He_plus_density[i][j][k];
                } else { //All hydrogen and helium atoms are ionized
                    relative_He_double_plus_density[i][j][k] = (n_e - nH - n_He)/nH; // Assuming all excess electrons go to He++
                    relative_Neutral_hydrogen_density[i][j][k] = 0.0; // Assuming all hydrogen is ionized
                    relative_Neutral_helium_density[i][j][k] = 0.0; // Assuming all helium is ionized
                    relative_H_plus_density[i][j][k] = 1.0;
                    relative_He_plus_density[i][j][k] = (n_He - (n_e - nH - n_He))/nH; // Remaining He+ after accounting for He++, assuming that all electrons that are result of ionization of atoms.
                }
                if (relative_He_plus_density[i][j][k] < 0.0) {
                    relative_He_plus_density[i][j][k] = 0.0; // Ensure no negative densities
                }
                if (relative_He_double_plus_density[i][j][k] < 0.0) {
                    relative_He_double_plus_density[i][j][k] = 0.0; // Ensure no negative densities
                }
                if (relative_Neutral_hydrogen_density[i][j][k] < 0.0) {
                    relative_Neutral_hydrogen_density[i][j][k] = 0.0; // Ensure no negative densities
                }
                if (relative_Neutral_helium_density[i][j][k] < 0.0) {
                    relative_Neutral_helium_density[i][j][k] = 0.0; // Ensure no negative densities
                }
                if (relative_H_plus_density[i][j][k] < 0.0) {
                    relative_H_plus_density[i][j][k] = 0.0; // Ensure no negative densities
                }
            }
        }    
    }
    phys_var.push_back(Temperature); // Add the Temperature array to phys_var for further processing if needed
    phys_var.push_back(relative_electron_density); // Add the electron density array to phys_var for further processing if needed
    phys_var.push_back(hydrogen_number_density); // Add the hydrogen number density array to phys_var for further processing if needed
    phys_var.push_back(relative_helium_number_density); // Add the helium number density array to phys_var for further processing if needed

    // Add the ion densities to phys_var for further processing if needed
    phys_var.push_back(relative_Neutral_hydrogen_density);
    phys_var.push_back(relative_Neutral_helium_density);
    phys_var.push_back(relative_H_plus_density);
    phys_var.push_back(relative_He_plus_density);
    phys_var.push_back(relative_He_double_plus_density);

    // Adding the magnetic field variable.
    // Zeroth order approximation of the magnetic field: B= B_0 \hat{z}
    double B_0 = 1e-6; // Example value for the magnetic field strength in Gauss
    std::vector<std::vector<std::vector<double>>> magnetic_field_x(hinfo.global_nx, std::vector<std::vector<double>>(hinfo.global_ny, std::vector<double>(hinfo.global_nz, 0.0)));
    std::vector<std::vector<std::vector<double>>> magnetic_field_y(hinfo.global_nx, std::vector<std::vector<double>>(hinfo.global_ny, std::vector<double>(hinfo.global_nz, 0.0)));
    std::vector<std::vector<std::vector<double>>> magnetic_field_z(hinfo.global_nx, std::vector<std::vector<double>>(hinfo.global_ny, std::vector<double>(hinfo.global_nz, B_0)));
    phys_var.push_back(magnetic_field_x); // Add magnetic field x-component
    phys_var.push_back(magnetic_field_y); // Add magnetic field y-component
    phys_var.push_back(magnetic_field_z); // Add magnetic field z-component



    // Calculate kinetic energy density in ergs/cm^3
    std::vector<std::vector<std::vector<double>>> kinetic_energy = kinetic_energy_density(phys_var[1], phys_var[2], phys_var[3], phys_var[0]); //input is velocity_x, velocity_y, velocity_z, density
    if (kinetic_energy.empty()) {
        std::cerr << "Failed to calculate kinetic energy density.\n";
        return 1;
    }

    //print the indices of the variables in phys_var
    std::cout << "Physical variables indices in phys_var:\n";
    std::cout << "0: Density (g/cm^3)\n";
    std::cout << "1: Velocity_x (cm/s)\n";
    std::cout << "2: Velocity_y (cm/s)\n";
    std::cout << "3: Velocity_z (cm/s)\n";
    std::cout << "4: Internal Energy Density (ergs/cm^3)\n";
    std::cout << "5: Temperature (K)\n";
    std::cout << "6: Relative Electron Density\n";
    std::cout << "7: Hydrogen Number Density (cm^-3)\n";
    std::cout << "8: Relative Helium Number Density\n";
    std::cout << "9: Relative Neutral Hydrogen Density\n";
    std::cout << "10: Relative Neutral Helium Density\n";
    std::cout << "11: Relative H+ Density\n";
    std::cout << "12: Relative He+ Density\n";
    std::cout << "13: Relative He++ Density\n";
    std::cout << "14: Magnetic Field X-component (G)\n";
    std::cout << "15: Magnetic Field Y-component (G)\n";
    std::cout << "16: Magnetic Field Z-component (G)\n";

    std::vector<int> indices = {29, 6, 20}; // Example indices to validate
    //last argument is true, which means we want to validate the data along with printing the values at the specified indices. Turn it to false if you only want to validate the data without printing the values.
    std::cout << "\n\nValidating data...";
    validate(phys_var, hinfo, indices, false);
    std::cout<< "Total kinetic energy in the domain is " << total_amount(kinetic_energy, hinfo) << " ergs.\n";
    std::cout << "Data validation completed.\n";
    
    // Write the phys_var data to a single CSV file
    std::string plt_filename= "plt80000";

    // Check if the plt_filename is part of the file_path
    if (file_path.find(plt_filename) == std::string::npos) {
        std::cerr << "Error: plt_filename '" << plt_filename << "' not found in file_path '" << file_path << "'.\n";
        return 1;
    }
    // check if the files are already present
    if (!fs::exists(plt_filename + ".csv")) {
        write_phys_var_to_single_csv(phys_var, plt_filename);
        std::cout << "Physical variables written to CSV file successfully.\n";
        write_phys_var_to_multi_csv(phys_var, plt_filename);
        std::cout << "Physical variables written to multiple CSV files successfully.\n";
    }
    std::cout << " n_x: " << hinfo.global_nx << ", n_y: " << hinfo.global_ny << ", n_z: " << hinfo.global_nz << "\n";
    // print the domain size, range, cell size and time
    std::cout << "Domain size: (" << hinfo.domain_hi[0] - hinfo.domain_lo[0] << ", "
              << hinfo.domain_hi[1] - hinfo.domain_lo[1] << ", "
              << hinfo.domain_hi[2] - hinfo.domain_lo[2] << ")\n";
    std::cout << "Domain range: (" << hinfo.domain_lo[0] << ", "
              << hinfo.domain_lo[1] << ", "
              << hinfo.domain_lo[2] << ") to ("
              << hinfo.domain_hi[0] << ", "
              << hinfo.domain_hi[1] << ", "
              << hinfo.domain_hi[2] << ")\n";
    std::cout << "Cell size(in cm): (" << (hinfo.domain_hi[0]-hinfo.domain_lo[0])/hinfo.global_nx << ", "
              << (hinfo.domain_hi[1]-hinfo.domain_lo[1])/hinfo.global_ny << ", "
              << (hinfo.domain_hi[2]-hinfo.domain_lo[2])/hinfo.global_nz << ")\n";
    double parsec_to_cm = 3085677581491367000; // 1 parsec in cm
    std::cout << "Cell size(in parsec): (" 
              << (hinfo.domain_hi[0]-hinfo.domain_lo[0])/(hinfo.global_nx * parsec_to_cm) << ", "
              << (hinfo.domain_hi[1]-hinfo.domain_lo[1])/(hinfo.global_ny * parsec_to_cm) << ", "
              << (hinfo.domain_hi[2]-hinfo.domain_lo[2])/(hinfo.global_nz * parsec_to_cm) << ")\n";
    return 0;
}

