#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include "read_header.hpp"



int main() {
    std::string path =  "/Users/aditivijayan/Projects/SKA-obs/data/plt00001"; 
    std::string header = path + "/Header";

    std::cout << "Loading data from path: " << path << "\n";
    HeaderInfo hinfo = read_quokka_header(header);
    
    // std::cout << "Dimensionality: " << hinfo.dim << "\n";
    // std::cout << "Variables (" << hinfo.num_components << "):\n";
    // std::cout << "Dimensions : nx: " << hinfo.global_nx << " : ny :" << hinfo.global_ny << "\n";
    // for (const auto& name : hinfo.variable_names) {
    //     std::cout << "  " << name << "\n";
    // }

        // std::cout << "Boxes:\n";
        // for (const auto& b : hinfo.boxes) {
        //     std::cout << "  [" << b.xlo << "," << b.ylo << "," << b.zlo << "] to ["
        //               << b.xhi << "," << b.yhi << "," << b.zhi << "]\n";
        // }

    
   

    std::string level_path = "/Users/aditivijayan/Projects/SKA-obs/data/plt00001/Level_0";
    std::vector<std::string> cell_files = get_all_cell_files(level_path);

    std::vector<BlockData> all_blocks;
    int i=0;
    for (const auto& cell_file : cell_files) {
        load_data(cell_file, all_blocks);
        i+=1;
    }
    double data = all_blocks[0].data(12, 1, 12);
    std::cout << "Data from all_blocks: " << data << "\n";

    //Initialize the data array for density and vx
    std::vector<std::vector<std::vector<double>>> density(hinfo.global_nx,
        std::vector<std::vector<double>>(hinfo.global_ny,
          std::vector<double>(hinfo.global_nz, 0.0)));
      
    // std::vector<std::vector<std::vector<double>>> vx = density; // copy structure

    BlockData block  = all_blocks[0];
    std::cout << "Block data: " << block.data(12, 1, 12) << "\n";
   
    return 0;
}
