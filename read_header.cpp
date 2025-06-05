#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include "read_header.hpp"

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
    
    for(int i=0; i<4; ++i) { // Read 4 variables: density, momentum_x, momentum_y, momentum_z
        for (const auto& cell_file : cell_files) {
            //std::cout << "Loading data from file: " << cell_file << " for variable index: " << i << "\n";
            load_data(cell_file, all_blocks_variables[i], i);
        }
    }
    for (const auto& cell_file : cell_files) {
        // Load internal energy data (assuming it's the 4th variable loaded into the index 4)
        load_data(cell_file, all_blocks_variables[4], 5); // Assuming 5 is the index for internal energy density 
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

    // Calculate kinetic energy density in ergs/cm^3
    std::vector<std::vector<std::vector<double>>> kinetic_energy = kinetic_energy_density(phys_var[1], phys_var[2], phys_var[3], phys_var[0]); //input is velocity_x, velocity_y, velocity_z, density
    if (kinetic_energy.empty()) {
        std::cerr << "Failed to calculate kinetic energy density.\n";
        return 1;
    }

    // Validate the data
    std::vector<int> indices = {0, 55, 63}; // Example indices to validate
    //last argument is true, which means we want to validate the data along with printing the values at the specified indices. Turn it to false if you only want to validate the data without printing the values.
    std::cout << "Validating data...\n\n";
    validate(phys_var, hinfo, indices, true);
    std::cout << "Data validation completed.\n";

    return 0;
}
