#include <iostream>
#include <fstream>

int main() {
    std::string filename = "/Users/aditivijayan/Projects/SKA-obs/data/plt4900000/Level_0/Cell_D_00018";  // Replace with the actual filename
    std::ifstream file(filename, std::ios::binary);  // Open file in binary mode
    std::string header;
    // std::getline(filename, header, "\n"); 
    // std::cout << "Header: " << header << std::endl;

    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return 1;  // Exit with error code
    }
    int block_size_x = 2;  // Replace with actual block size
    int block_size_y = 2;  // Replace with actual block size
    int block_size_z = 2;  // Replace with actual block size

    double data[block_size_x][block_size_y][block_size_z];  // Replace with actual sizes
    for (int i = 0; i < block_size_x; ++i) {
        for (int j = 0; j < block_size_y; ++j) {
            for (int k = 0; k < block_size_z; ++k) {
                file.read(reinterpret_cast<char*>(&data[i][j][k]), sizeof(double));
                std::cout << "Data read: " << data[i][j][k] << std::endl;
            }
        }
    }

    std::cout << "File opened successfully: " << filename << std::endl;

    // You can now proceed with reading the file...

    file.close();  // Close the file when done
    return 0;  // Exit successfully
}
