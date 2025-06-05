#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <filesystem>
#include <algorithm>


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
    
    // Read current time
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

// Function to print the value of a physical variable at specific coordinates
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
        std::cout << "Density (in g/cm^3) at (" << x << ", " << y << ", " << z << "): ";
    } else if (var_index==1) {
        std::cout << "Velocity_x (in cm/s) at (" << x << ", " << y << ", " << z << "): ";
    } else if (var_index==2) {
        std::cout << "Velocity_y (in cm/s) at (" << x << ", " << y << ", " << z << "): ";
    } else if (var_index==3) {
        std::cout << "Velocity_z (in cm/s) at (" << x << ", " << y << ", " << z << "): ";
    } else if (var_index==4) {
        std::cout << "Internal_energy (in ergs) at (" << x << ", " << y << ", " << z << "): ";
    } 
    else {
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
    const std::vector<std::vector<std::vector<double>>>& velocity_x,
    const std::vector<std::vector<std::vector<double>>>& velocity_y,
    const std::vector<std::vector<std::vector<double>>>& velocity_z,
    const std::vector<std::vector<std::vector<double>>>& density
){
    //check if the dimensions of the input arrays match
    if (velocity_x.size() != velocity_y.size() || 
        velocity_x.size() != velocity_z.size() || 
        velocity_x.size() != density.size()) {
        std::cerr << "Dimension mismatch in kinetic energy density calculation.\n";
        return {};
    }

    std::vector<std::vector<std::vector<double>>> kinetic_energy(density.size(),
        std::vector<std::vector<double>>(density[0].size(),
            std::vector<double>(density[0][0].size(), 0.0)));

    for (size_t i = 0; i < density.size(); ++i) {
        for (size_t j = 0; j < density[i].size(); ++j) {
            for (size_t k = 0; k < density[i][j].size(); ++k) {
                double vel_sq = velocity_x[i][j][k] * velocity_x[i][j][k] +
                                velocity_y[i][j][k] * velocity_y[i][j][k] +
                                velocity_z[i][j][k] * velocity_z[i][j][k];
                kinetic_energy[i][j][k] = 0.5 * density[i][j][k] * vel_sq;
            }
        }
    }

    return kinetic_energy;
}


void validate(const std::vector<std::vector<std::vector<std::vector<double>>>>& phys_var, HeaderInfo& hinfo, std::vector<int> indices ={12,1,12}, bool c=true) {
    if (c) {
        // Check if the dimensions of the physical variables are consistent
        for (const auto& var : phys_var) {
            if (var.size() != phys_var[0].size() || var[0].size() != phys_var[0][0].size()) {
                std::cerr << "Inconsistent dimensions found in physical variables.\n";
                return;
            }
        }
        // Check if the indices are within bounds
        for (const auto& index : indices) {
            if (index < 0 || index >= phys_var[0].size()) {
                std::cerr << "Index out of bounds: " << index << "\n";
                return;
            }
        }


        std::cout << "All physical variables have consistent dimensions.\n";
        for (int i=0; i<phys_var.size(); i++){
            print_var_value_at_coordinates(phys_var, i, indices[0],indices[1], indices[2]);
        }
    }

    std::cout<< "Total mass of the domain is "<<total_amount(phys_var[0], hinfo)<< std::endl;
    std::cout<< "Total Internal Energy of the domain is "<<total_amount(phys_var[4], hinfo)<< std::endl;

        //Print the maximum and minimum values and positions of various physical parameters
    for (int var_idx = 0; var_idx < phys_var.size(); ++var_idx) {
            double min_val = std::numeric_limits<double>::max();
            double max_val = std::numeric_limits<double>::lowest();
            std::tuple<int, int, int> min_pos, max_pos;

        for (int i = 0; i < phys_var[var_idx].size(); ++i) {
            for (int j = 0; j < phys_var[var_idx][i].size(); ++j) {
                for (int k = 0; k < phys_var[var_idx][i][j].size(); ++k) {
                double val = phys_var[var_idx][i][j][k];
                if (val < min_val) {
                    min_val = val;
                    min_pos = std::make_tuple(i, j, k);
                }
                if (val > max_val) {
                    max_val = val;
                    max_pos = std::make_tuple(i, j, k);
                }
                }
            }
            }
    std::string var_name = hinfo.variable_names[var_idx];
        if (var_name == "gasDensity") {
            var_name = "Density (g/cm^3)";
        } else if (var_name == "x-GasMomentum") {
            var_name = "Velocity_x (cm/s)";
        } else if (var_name == "y-GasMomentum") {
            var_name = "Velocity_y (cm/s)";
        } else if (var_name == "z-GasMomentum") {
            var_name = "Velocity_z (cm/s)";
        } else if (var_name == "gasEnergy") {
            var_name = "Internal Energy (ergs)";
        }
        std::cout << "Variable: " << var_name << " (index " << var_idx << "): "
                << "min = " << min_val
                << " at (" << std::get<0>(min_pos) << ", " << std::get<1>(min_pos) << ", " << std::get<2>(min_pos) << ")"
                << ", max = " << max_val
                << " at (" << std::get<0>(max_pos) << ", " << std::get<1>(max_pos) << ", " << std::get<2>(max_pos) << ")\n";
    }
}

