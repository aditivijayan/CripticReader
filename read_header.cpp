#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include "read_header.hpp"

namespace fs = std::filesystem;

int main() {
    std::string home = fs::current_path();
    std::string data_path =  home + "/data/plt00001"; 
    std::string header = data_path + "/Header";

    std::cout << "Loading data from path: " << data_path << "\n";
    HeaderInfo hinfo = read_quokka_header(header);

    std::string level_path = data_path  + "/Level_0";
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
