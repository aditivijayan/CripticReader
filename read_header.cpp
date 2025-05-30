#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include "read_header.hpp"

namespace fs = std::filesystem;

void print_2D_data(const std::vector<std::vector<double>>& data_2D_vector) {
    for (const auto& row : data_2D_vector) {
        for (const auto& value : row) {
            std::cout << value << " ";
        }
        std::cout << "\n";
    }
}
void print_var_value_at_coordinates(const std::vector<std::vector<std::vector<std::vector<double>>>>& phys_var, int var_index, int x, int y, int z) {
    if (var_index < 0 || var_index >= phys_var.size()) {
        std::cerr << "Variable index out of bounds.\n";
        return;
    }
    if (x < 0 || x >= phys_var[var_index].size() ||
        y < 0 || y >= phys_var[var_index][0].size() ||
        z < 0 || z >= phys_var[var_index][0][0].size()) {
        std::cerr << "Coordinates out of bounds.\n";
        return;
    }
    if (var_index==0){
        std::cout << "Density at (" << x << ", " << y << ", " << z << "): ";
    } else if (var_index==1) {
        std::cout << "Momentum_x at (" << x << ", " << y << ", " << z << "): ";
    } else if (var_index==2) {
        std::cout << "Momentum_y at (" << x << ", " << y << ", " << z << "): ";
    } else if (var_index==3) {
        std::cout << "Momentum_z at (" << x << ", " << y << ", " << z << "): ";
    } else if (var_index==4) {
        std::cout << "Velocity_x at (" << x << ", " << y << ", " << z << "): ";
    } else if (var_index==5) {
        std::cout << "Velocity_y at (" << x << ", " << y << ", " << z << "): ";
    } else if (var_index==6) {
        std::cout << "Velocity_z at (" << x << ", " << y << ", " << z << "): ";
    } else if (var_index==7) {
        std::cout << "Kinetic_energy at (" << x << ", " << y << ", " << z << "): ";
    } else {
        std::cerr << "Unknown variable index.\n";
        return;
    }
    std::cout << phys_var[var_index][x][y][z] << "\n";
}

double total_amount(const std::vector<std::vector<std::vector<double>>>& variable_density, 
                  const HeaderInfo& hinfo) {
    double total = 0.0;
    for (const auto& plane : variable_density) {
        for (const auto& row : plane) {
            for (const auto& value : row) {
                total += value;
            }
        }
    }
    double dV= (hinfo.domain_hi[0] - hinfo.domain_lo[0]) * 
         (hinfo.domain_hi[1] - hinfo.domain_lo[1]) * 
         (hinfo.domain_hi[2] - hinfo.domain_lo[2]) / 
         (hinfo.global_nx * hinfo.global_ny * hinfo.global_nz);
    return total * dV;
}

std::vector<std::vector<std::vector<double>>> kinetic_energy_density(
    const std::vector<std::vector<std::vector<double>>>& velociy_x,
    const std::vector<std::vector<std::vector<double>>>& velociy_y,
    const std::vector<std::vector<std::vector<double>>>& velociy_z,
    const std::vector<std::vector<std::vector<double>>>& density
){
    //check if the dimensions of the input arrays match
    if (velociy_x.size() != velociy_y.size() || 
        velociy_x.size() != velociy_z.size() || 
        velociy_x.size() != density.size()) {
        std::cerr << "Dimension mismatch in kinetic energy density calculation.\n";
        return {};
    }

    std::vector<std::vector<std::vector<double>>> kinetic_energy(density.size(),
        std::vector<std::vector<double>>(density[0].size(),
            std::vector<double>(density[0][0].size(), 0.0)));

    for (size_t i = 0; i < density.size(); ++i) {
        for (size_t j = 0; j < density[i].size(); ++j) {
            for (size_t k = 0; k < density[i][j].size(); ++k) {
                double vel_sq = velociy_x[i][j][k] * velociy_x[i][j][k] +
                                velociy_y[i][j][k] * velociy_y[i][j][k] +
                                velociy_z[i][j][k] * velociy_z[i][j][k];
                kinetic_energy[i][j][k] = 0.5 * density[i][j][k] * vel_sq;
            }
        }
    }

    return kinetic_energy;
}

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
    //std::cout << "Number of cell files found: " << cell_files.size() << "\n";
    std::vector<std::vector<BlockData>> all_blocks_variables(4, std::vector<BlockData>());  // i-th entry of all_block_variables will hold all blocks for the i-th variable, j-th entry of the i-th all_block would hold all blocks for the j-th cell file
    //std::cout << "Number of cell files found: " << cell_files.size() << "\n";
    if (cell_files.empty()) {
        std::cerr << "No cell files found in the specified directory." << std::endl<< std::endl;
        return 1;
    }
    //int i=;
    for(int i=0; i<4; ++i) { // Read 4 variables: density, momentum_x, momentum_y, momentum_z
        for (const auto& cell_file : cell_files) {
         
            //std::cout << "Loading data from file: " << cell_file << " for variable index: " << i << "\n";
            load_data(cell_file, all_blocks_variables[i], i);
        }
        // std::cout << "Loading data from file: " << cell_file << "\n";
        //load_data(cell_file, all_blocks);
        //i+=1;
    }
    std::cout << "Total number of variables loaded: " << all_blocks_variables.size() << "\n";
    // all_blocks[0] would would contain the density data, all_blocks[1] would contain the momentum_x data, etc.
    // double data = all_blocks_variables[1][0].data(12, 1, 12);      // data of the first Cell_D_ file for the second variable (momentum_x)
    // std::cout << "data of the first Cell_D_ file for the second variable (momentum_x) at (12,1,12) " << data << "\n";

    std::vector<std::vector<double>> data_2D_vector;
    // store the data from the all_blocks_variable[i][j], i.e., the data of i-th physical variable in the j-th Cell_D dile as the i,j-th entry of the data_vector

    for (int i = 0; i < 4; ++i) {  // Loop over each variable
        std::vector<double> temp_vector;
        for (const auto& block : all_blocks_variables[i]) {  // Loop over each Cell_D file
            // Assuming we want to extract data at a specific index, e.g., (12, 1, 12)
            temp_vector.push_back(block.data(12, 1, 12)); // Example indices, adjust as needed
        }
        data_2D_vector.push_back(temp_vector);
    }
    // for (const auto& block : all_blocks_variables[1]) {
    //     data_vector.push_back(block.data(12, 1, 12)); // Example indices, adjust as needed
    // }
    //std::cout << "Data vector size: " << data_2D_vector.size() << "\n";
    //std::cout << "Data vector size of the first variable,i.e., number of Cell_D files: " << data_2D_vector[0].size() << "\n";



    // Print the entries of the data_2D_vector for each variable
    // std::cout << "Data for each variable at (12, 1, 12):\n";
    // print_2D_data(data_2D_vector);


    //Initialize the data array for density and momentum_x, momentum_y, momentum_z
    //Assuming global_nx, global_ny, global_nz are the dimensions of the data
    std::vector<std::vector<std::vector<double>>> density(hinfo.global_nx, std::vector<std::vector<double>>(hinfo.global_ny,std::vector<double>(hinfo.global_nz, 0.0)));
    std::vector<std::vector<std::vector<double>>> p_x(hinfo.global_nx, std::vector<std::vector<double>>(hinfo.global_ny,std::vector<double>(hinfo.global_nz, 0.0)));
    std::vector<std::vector<std::vector<double>>> p_y(hinfo.global_nx, std::vector<std::vector<double>>(hinfo.global_ny,std::vector<double>(hinfo.global_nz, 0.0)));
    std::vector<std::vector<std::vector<double>>> p_z(hinfo.global_nx, std::vector<std::vector<double>>(hinfo.global_ny,std::vector<double>(hinfo.global_nz, 0.0)));

    std::vector<std::vector<std::vector<std::vector<double>>>> phys_var;
    phys_var.push_back(density);
    phys_var.push_back(p_x);
    phys_var.push_back(p_y);
    phys_var.push_back(p_z);
    std::cout << "Initialized data arrays for density and momentum variables.\n";

    // Fill the data arrays with the data from all_blocks
    //std::cout<< "size of all_blocks_variables: " << all_blocks_variables.size() << "\n";

    for (int var = 0; var < all_blocks_variables.size(); ++var) {
        for (int blk = 0; blk < all_blocks_variables[var].size(); ++blk) {
            const BlockData& block = all_blocks_variables[var][blk];
            //std::cout << "Processing variable: " << var << ", block: " << blk << " with start indices (" << block.start_x << ", " << block.start_y << ", " << block.start_z << ") and end indices (" 
            //          << block.end_x << ", " << block.end_y << ", " << block.end_z << ")\n";
            if (block.start_x < 0 || block.end_x >= hinfo.global_nx ||
                block.start_y < 0 || block.end_y >= hinfo.global_ny ||
                block.start_z < 0 || block.end_z >= hinfo.global_nz) {
                std::cerr << "Block indices out of bounds for global dimensions.\n";
                continue;
            }
            //std::cout << "Filling data for variable: " << var << ", block: " << blk << "\n";
            for (int i = block.start_x; i <= block.end_x; ++i) {
                for (int j = block.start_y; j <= block.end_y; ++j) {
                    for (int k = block.start_z; k <= block.end_z; ++k) {
                        phys_var[var][i][j][k] = block.data(i - block.start_x, j - block.start_y, k - block.start_z);
                    }
                }
            }
    }
    }

    int i=12;
    int j=1;
    int k=12;
    std::cout << "Values at (" << i << ", " << j << ", " << k << "):\n";
    for (int var = 0; var < phys_var.size(); ++var) {
        // Print the value of the variable at (12, 1, 12)
        print_var_value_at_coordinates(phys_var, var, i, j, k);
    }






    // Velocities from momentum variables//velocity=momentum/density
    // Assuming the momentum variables are in the second, third, and fourth entries of all_blocks_variables
    // std::vector<BlockData> v_x;
    // std::vector<BlockData> v_y;
    // std::vector<BlockData> v_z;
    std::vector<std::vector<BlockData>> all_blocks_variables_velocity(3, std::vector<BlockData>());

    for (int i=1; i<4; ++i) { // Loop over momentum variables
        for (int j=0; j<all_blocks_variables[i].size(); j++) {
            BlockData velocity_block = all_blocks_variables[i][j] / (all_blocks_variables[0][j]); // velocity = momentum / density (in cm/s)
            all_blocks_variables_velocity[i-1].push_back(velocity_block);
            
        }
        //std::cout << "Velocity block for variable " << i << " added.\n";
    }
    
    //To get velocity in km/s, we can divide the velocity by 100000
    for (auto& block : all_blocks_variables_velocity) {
        for (auto& data_block : block) {
            data_block /= 100000.0; // Convert to km/s
        }
    }
    std::cout << "Velocity blocks created and converted to km/s.\n";

    // Convert velocity blocks to 3D arrays and add them to phys_var
    for (const auto& velocity_blocks : all_blocks_variables_velocity) {
        std::vector<std::vector<std::vector<double>>> velocity_array(
            hinfo.global_nx, std::vector<std::vector<double>>(hinfo.global_ny, std::vector<double>(hinfo.global_nz, 0.0)));
        for (const auto& block : velocity_blocks) {
            if (block.start_x < 0 || block.end_x >= hinfo.global_nx ||
                block.start_y < 0 || block.end_y >= hinfo.global_ny ||
                block.start_z < 0 || block.end_z >= hinfo.global_nz) {
                std::cerr << "Velocity block indices out of bounds for global dimensions.\n";
                continue;
            }
            for (int i = block.start_x; i <= block.end_x; ++i) {
                for (int j = block.start_y; j <= block.end_y; ++j) {
                    for (int k = block.start_z; k <= block.end_z; ++k) {
                        velocity_array[i][j][k] = block.data(i - block.start_x, j - block.start_y, k - block.start_z);
                    }
                }
            }
        }
        phys_var.push_back(velocity_array);
    }
    std::cout << "Velocity arrays added to phys_var.\n";
    std::cout << "Total number of physical variables (size of phys_var (7 parameters*3D array of doubles)): " << phys_var.size() << "\n";

    // Print the velocity data at a specific index, e.g., (12, 1, 76)
    // std::cout << "Velocity values at (12, 1, 76):\n";
    // for (int var = 0; var < 3; ++var) { // Loop over velocity variables
    //     std::cout << "Velocity Variable " << var+4 << ": " << phys_var[4 + var][12][1][76] << " km/s\n"; // 4 + var to account for density and momentum variables
    // }

    //join the all_blocks_variables and all_blocks_variables_velocity into a single vector
    all_blocks_variables.insert(all_blocks_variables.end(), all_blocks_variables_velocity.begin(), all_blocks_variables_velocity.end());
    std::cout << "All blocks variables and velocity blocks merged.\n";
    std::cout << "Total number of physical variables(size of all_blocks_variables (2D array of BlockData)): " << all_blocks_variables.size() << "\n";




    double mass = total_amount(phys_var[0], hinfo);
    std::cout << "Total mass in the domain: " << mass << " g\n";

    // Calculate kinetic energy density
    std::vector<std::vector<std::vector<double>>> kinetic_energy = kinetic_energy_density(phys_var[4], phys_var[5], phys_var[6], phys_var[0]);
    if (kinetic_energy.empty()) {
        std::cerr << "Failed to calculate kinetic energy density.\n";
        return 1;
    }

    //append kinetic energy density to phys_var
    phys_var.push_back(kinetic_energy);
    std::cout << "Kinetic energy density calculated and added to phys_var at position " << phys_var.size() - 1 << ".\n";

    // Print the kinetic energy density at a specific index, e.g., (12, 1, 76)
    print_var_value_at_coordinates(phys_var, 7, 12, 1, 12); // 7 is the index for kinetic energy density

    double total_kinetic_energy = total_amount(kinetic_energy, hinfo);
    std::cout << "Total kinetic energy in the domain: " << total_kinetic_energy << " erg\n";

    return 0;
}
