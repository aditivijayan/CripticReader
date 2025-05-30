#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <filesystem>
#include <algorithm>

// // The purpose of Array3DView to be able to access data by specifying its i,j,k
// struct Array3DView {
//     std::vector<double>& data;
//     int nx, ny, nz;

//     Array3DView(std::vector<double>& data, int nx, int ny, int nz)
//         : data(data), nx(nx), ny(ny), nz(nz) {}

//     double& operator()(int i, int j, int k) {
//         return data[k * nx * ny + j * nx + i];
//     }

//     const double& operator()(int i, int j, int k) const {
//         return data[k * nx * ny + j * nx + i];
//     }
// };

struct Array3DView {
    const double* data_ptr;  // Pointer to external data
    //int nx, ny, nz;
    int nx=64;
    int ny=64;
    int nz=128;
    Array3DView(const std::vector<double>& data, int nx_, int ny_, int nz_)
        : data_ptr(data.data()), nx(nx_), ny(ny_), nz(nz_) {
            
        }

    double operator()(int i, int j, int k) const {
        return data_ptr[k * ny * nx + j * nx + i];
    }

    double& operator()(int i, int j, int k) {
        //std::cout << std::endl << "nx: " << nx << ", ny: " << ny << ", nz: " << nz << std::endl;
        return const_cast<double&>(data_ptr[k * nx * ny + j * nx+ i]);
    }
};


struct Box {
    double xlo, ylo, zlo, xhi, yhi, zhi;
};

struct BlockData {
    int start_x, start_y, start_z;
    int end_x, end_y, end_z;
    std::vector<double> storage;
    Array3DView data;

    // Constructor from rvalue vector
    BlockData(std::vector<double>&& d, int nx, int ny, int nz)
        : storage(std::move(d)), data(storage, nx, ny, nz) {}

    // Copy constructor
    BlockData(const BlockData& other)
        : start_x(other.start_x), start_y(other.start_y), start_z(other.start_z),
          end_x(other.end_x), end_y(other.end_y), end_z(other.end_z),
          storage(other.storage),  // copies vector
          data(storage, other.data.nx, other.data.ny, other.data.nz) {}

    // Move constructor
    BlockData(BlockData&& other) noexcept
        : start_x(other.start_x), start_y(other.start_y), start_z(other.start_z),
          end_x(other.end_x), end_y(other.end_y), end_z(other.end_z),
          storage(std::move(other.storage)),
          data(storage, other.data.nx, other.data.ny, other.data.nz) {}

    // Copy assignment
    BlockData& operator=(const BlockData& other) {
        if (this != &other) {
            start_x = other.start_x;
            start_y = other.start_y;
            start_z = other.start_z;
            end_x = other.end_x;
            end_y = other.end_y;
            end_z = other.end_z;
            storage = other.storage;
            data = Array3DView(storage, other.data.nx, other.data.ny, other.data.nz);
        }
        return *this;
    }

    // Move assignment
    BlockData& operator=(BlockData&& other) noexcept {
        if (this != &other) {
            start_x = other.start_x;
            start_y = other.start_y;
            start_z = other.start_z;
            end_x = other.end_x;
            end_y = other.end_y;
            end_z = other.end_z;
            storage = std::move(other.storage);
            data = Array3DView(storage, other.data.nx, other.data.ny, other.data.nz);
        }
        return *this;
    }

    // Divide operator
    BlockData& operator/=(double divisor) {
        for (auto& value : storage) {
            value /= divisor;
        }
        return *this;
    }
    // Divide operator for another BlockData object
BlockData operator/(const BlockData& other) const {
    if (storage.size() != other.storage.size()) {
        std::cerr << "Sizes: " << storage.size() << " and " << other.storage.size() << std::endl;
        throw std::runtime_error("BlockData sizes do not match for division.");
    }
    BlockData result = *this;
    for (size_t i = 0; i < storage.size(); ++i) {
        result.storage[i] /= other.storage[i];
    }
    return result;
}

};

struct HeaderInfo {
    int dim;
    int num_components;
    int global_nx;
    int global_ny;
    int global_nz;
    float curr_time;
    double domain_lo[3];
    double domain_hi[3];
    std::vector<std::string> variable_names;
    std::vector<Box> boxes;
};

//The following function isn't used in the code
void read_cell_data(const std::string& cell_file, const Box& box, std::vector<float>& data) { // check if you need to modify this for v_x instead of density
    std::ifstream infile(cell_file, std::ios::binary);
    if (!infile) {
        throw std::runtime_error("Unable to open cell file: " + cell_file);
    }

    int nx = box.xhi - box.xlo + 1;
    int ny = box.yhi - box.ylo + 1;
    int nz = box.zhi - box.zlo + 1;
    int num_cells = nx * ny * nz;
    

    //print the location of the pointer
    std::cout << "Current position: " << infile.tellg() << std::endl;

    // Seek to the offset (skip first num_cells floats) // To read the velocity information
    //infile.seekg(num_cells * sizeof(float), std::ios::beg);

    //print the location of the pointer
    std::cout << "Current position after seek: " << infile.tellg() << std::endl;

    // Resize the data array
    data.resize(num_cells);

    // Read binary data into the array
    infile.read(reinterpret_cast<char*>(data.data()), num_cells * sizeof(float));  // Check this line, especially the starting pointer 

    if (!infile) {
        throw std::runtime_error("Error reading data from: " + cell_file);
    }
}


//Go over a single cell file and store all the data in a BlockData object

void load_data(const std::string& filename, std::vector<BlockData>& blocks, int phys_var) { //phys_var is the variable to read, e.g. 0 for density, 1 for momentum_x, etc.
    std::ifstream infile(filename);
    //std::cout << "Loading data from: " << filename << std::endl;
    std::string line;

    if (!infile.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }

    while (std::getline(infile, line)) {
        // Check for lines starting with "FAB"
        if (line.find("FAB") == 0) {
            // Parse the FAB line for start and finish indices
            std::istringstream line_stream(line);
            std::string fab;
            int num_var;
            int block_id, start_x, start_y, start_z;
            int end_x, end_y, end_z;
            std::size_t coords_start = line.find(")))((") + 5;
            std::istringstream coord_stream(line.substr(coords_start));
            char ch;
            coord_stream >>  start_x >> ch >> start_y >> ch >>  start_z >> ch; // (63,63,127)
            coord_stream >> ch >> end_x >> ch >> end_y >> ch >> end_z >> num_var;
            
            int nx = end_x - start_x + 1;
            int ny = end_y - start_y + 1;
            int nz = end_z - start_z + 1;

            // Now read the data for this block (following the "FAB" line)
            std::vector<std::vector<std::vector<double>>> data(end_x - start_x + 1,
            std::vector<std::vector<double>>(end_y - start_y + 1,
            std::vector<double>(end_z - start_z + 1)));
            int num_cells_block = (end_x - start_x + 1) * (end_y - start_y + 1) * (end_z - start_z + 1);
            
            std::vector<double> flat_data(num_cells_block);
            //print the location of the pointer
            //std::cout << "Current position: " << infile.tellg() << std::endl;
            infile.seekg(phys_var * num_cells_block * sizeof(double), std::ios::cur);  // Skip first block, for reading velocity information
            //print the location of the pointer
            //std::cout << "Current position after seek: " << infile.tellg() << std::endl;
            infile.read(reinterpret_cast<char*>(flat_data.data()), num_cells_block * sizeof(double));  //check this line, especially the starting pointer
                // int i = 0;
                // int j = 0;
                // int k = 0;
                // for (const auto& value : flat_data) {
                //     int index = k * nx * ny + j * nx + i; // Assuming row-major order

                //     double value = flat_data[index];
                //     printf("Data: %.10e\n", value);

            // BlockData block(flat_data, nx, ny, nz);
            BlockData block(std::move(flat_data), nx, ny, nz);
            block.start_x = start_x;
            block.start_y = start_y;
            block.start_z = start_z;
            block.end_x = end_x;
            block.end_y = end_y;
            block.end_z = end_z;
            blocks.push_back(block);
            
        }
    }

    infile.close();
}



HeaderInfo read_quokka_header(const std::string& header_path) {
    std::ifstream infile(header_path);
    if (!infile) {
        throw std::runtime_error("Unable to open header file.");
    }

    HeaderInfo info;
    std::string line;

    // Skip format line
    std::getline(infile, line);
    std::getline(infile, line);
    info.num_components = 0;
    //Read variable names
    while(line.find("z-velocity")!=0 ) {    
        std::getline(infile, line);
        info.variable_names.push_back(line);
        info.num_components++;
    }
    
    // Read dimensionality
    std::getline(infile, line);
    info.dim = std::stoi(line);
    
    // Read currrent time
    std::getline(infile, line); 
    info.curr_time = std::stol(line);

    std::getline(infile, line); // Skip the next line 

    // Read domain limits
    
    std::getline(infile, line);
    std::istringstream lo_stream(line);

    // Read lower bounds
    lo_stream >> info.domain_lo[0] >> info.domain_lo[1] >> info.domain_lo[2];

    // Read upper bounds
    std::getline(infile, line);
    std::istringstream hi_stream(line);
    hi_stream >> info.domain_hi[0] >> info.domain_hi[1] >> info.domain_hi[2];

    std::getline(infile, line);
    std::getline(infile, line);

    //Read domain size 
    std::istringstream domsize_stream(line);
    char ch;
    int dummy;
    domsize_stream >> ch >> ch ; //skip "(("
    domsize_stream >> dummy >> ch >> dummy >> ch >> dummy >> ch;
    domsize_stream >> ch;
    domsize_stream >> info.global_nx >> ch >> info.global_ny >> ch >> info.global_nz;
    info.global_nx += 1; // Adjust for 0-based indexing
    info.global_ny += 1;
    info.global_nz += 1;

    //Read extra lines 
    for (int i = 19; i<23; ++i) {
        std::getline(infile, line);
    }
    std::getline(infile, line);
    
    std::istringstream stream(line); 
    int dum, num_boxes;
    stream >> dum >> num_boxes;
    std::getline(infile, line); // file number 
    
    // Read boxes
    for (int i = 0; i < num_boxes; ++i) {
        Box b;
        std::istringstream ss;
        //Read X limits
        std::getline(infile, line);
        ss.clear(); ss.str(line);
        ss  >> b.xlo >> b.xhi; 

        //Read Y limits
        std::getline(infile, line);
        ss.clear(); ss.str(line);
        ss >> b.ylo >> b.yhi;

        //Read Z limits
        std::getline(infile, line);
        ss.clear(); ss.str(line);
        ss >> b.zlo >> b.zhi; 

        //Push back the box
        info.boxes.push_back(b);
    }

    // Read number of components (nvar)
    std::getline(infile, line);
    return info;
}


namespace fs = std::filesystem;

std::vector<std::string> get_all_cell_files(const std::string& level_dir) {
    std::vector<std::string> cell_files;
    for (const auto& entry : std::filesystem::directory_iterator(level_dir)) {
        if (entry.is_regular_file()) {
            std::string name = entry.path().filename().string();
            if (name.rfind("Cell_D_", 0) == 0) {  // name starts with "Cell_D_"
                cell_files.push_back(entry.path().string());
            }
        }
    }
    std::sort(cell_files.begin(), cell_files.end());  // optional but useful
    return cell_files;
}