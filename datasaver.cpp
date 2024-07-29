//
// Created by konrad_guest on 25/07/2024.
//

#include "Headers/utils/datasaver.h"

namespace mds {
    std::string DataSaver2D::generateDirectoryName(const std::shared_ptr<Simulation>& simulation)
    {
        std::shared_ptr<simulation2dap> sim2D = std::dynamic_pointer_cast<simulation2dap>(simulation);

        auto now = std::chrono::system_clock::now();
        std::time_t now_time = std::chrono::system_clock::to_time_t(now);

        std::stringstream ss;
        ss << std::put_time(std::localtime(&now_time), "%a_%b_%d_%H_%M");

        std::string date_time = ss.str();

        std::string num_atoms,density,iterations;
        try {
            num_atoms = std::to_string(sim2D->get<int>("numberOfAtoms"));
            iterations = std::to_string(sim2D->get<int>("stepLimit"));
            density = std::to_string(sim2D->get<double>("density"));
        }catch(...) {
            num_atoms = "default_num_atoms";
            iterations = "default_iterations";
            density = "default_density";
            throw mds::eeFileIO;
        }

        std::string current_path = std::filesystem::current_path().string();

        std::string dir_name = current_path + "\\Resources\\" + "2D" + "_" + date_time + "_" + num_atoms + "_" + iterations + "_" + density;

        return dir_name;
    }

    void DataSaver2D::createDirectory(const std::string& dir_name)
    {
        std::filesystem::path dir_path(dir_name);
        if (!std::filesystem::exists(dir_path))
        {
            std::filesystem::create_directories(dir_path);
            std::cout << "Directory created: " << dir_path << std::endl;
        }
    }

    void DataSaver2D::openFile(std::ofstream& file, const std::string& filePath, unsigned int mode)
    {
        if(mode == 0) {
            file.open(filePath, std::ios::out | std::ios::app);
        }else if(mode == 1) {
            file.open(filePath, std::ios::out | std::ios::trunc);
        }
        if (!file.is_open()) {
            std::cerr << "Error: Could not open file for writing." << "\n";
        }
    }

    void DataSaver2D::closeFile(std::ofstream& file)
    {
        if(file.is_open()) {file.close();}
    }


    void DataSaver2D::savePropertiesData(
        const std::string&                                                                                          filename,
        const std::vector<std::tuple<int, double, Eigen::VectorXd, double, double, double, double, double, double>>&    data,
        unsigned int                                                                                                    mode,
        const std::shared_ptr<Simulation>&                                                                        simulation,
        std::ofstream&                                                                                  propertiesDataOutput
        )
    {
        std::string dir_name            = generateDirectoryName(simulation);
        createDirectory(dir_name);

        std::string relativePath        = "properties.csv";
        std::filesystem::path filePath  = dir_name + "\\" + relativePath;

        openFile(propertiesDataOutput, filePath.string(), mode);

        propertiesDataOutput << "Step Count,Time Now,Velocity Sum Value X,Velocity Sum Value Y,Potential Energy Sum,"
             << "Potential Energy Sum Squared,Kinetic Energy Sum,Kinetic Energy Sum Squared,Pressure Sum,Pressure Sum Squared\n";

        for (const auto& entry : data)
        {
            propertiesDataOutput << std::setw(5) << std::get<0>(entry) << ","
                                 << std::setw(8) << std::fixed << std::setprecision(4) << std::get<1>(entry) << ","
                                 << std::setw(8) << std::fixed << std::setprecision(4) << std::get<2>(entry).x() << ","
                                 << std::setw(8) << std::fixed << std::setprecision(4) << std::get<2>(entry).y() << ","
                                 << std::setw(8) << std::fixed << std::setprecision(4) << std::get<3>(entry) << ","
                                 << std::setw(8) << std::fixed << std::setprecision(4) << std::get<4>(entry) << ","
                                 << std::setw(8) << std::fixed << std::setprecision(4) << std::get<5>(entry) << ","
                                 << std::setw(8) << std::fixed << std::setprecision(4) << std::get<6>(entry) << ","
                                 << std::setw(8) << std::fixed << std::setprecision(4) << std::get<7>(entry) << ","
                                 << std::setw(8) << std::fixed << std::setprecision(4) << std::get<8>(entry) << "\n";
        }

        closeFile(propertiesDataOutput);
    }

    void DataSaver2D::saveVelocityData(
        const std::string&                                filename,
        const std::vector<std::vector<double>>&               data,
        unsigned int                                          mode,
        const std::shared_ptr<Simulation>&              simulation,
        std::ofstream&              VelocityDistributionDataOutput
        )
    {
        std::string dir_name            = generateDirectoryName(simulation);
        createDirectory(dir_name);

        std::string relativePath        = "velocity_distribution.csv";
        std::filesystem::path filePath  = dir_name + "\\" + relativePath;

        openFile(VelocityDistributionDataOutput, filePath.string(), mode);

        VelocityDistributionDataOutput << "Velocity distribution attitude\n";

        for (const auto& row : data)
        {
            for (const auto& column : row)
            {
                VelocityDistributionDataOutput << std::setw(8) << std::fixed << std::setprecision(4) << column << ",";
            }
            VelocityDistributionDataOutput << "\n---\n";
        }

        closeFile(VelocityDistributionDataOutput);
    }

    void DataSaver2D::saveCoordinateData(
            const std::string&                                    filename,
            const std::vector<std::tuple<int, Eigen::VectorXd>>&      data,
            unsigned int                                              mode,
            const std::shared_ptr<Simulation>&                  simulation,
            std::ofstream&                           coordinatesDataOutput
            )
    {
        std::string dir_name            = generateDirectoryName(simulation);
        createDirectory(dir_name);

        std::string relativePath        = "results.csv";
        std::filesystem::path filePath  = dir_name + "\\" + relativePath;

        openFile(coordinatesDataOutput, filePath.string(), mode);

        for (const auto& entry : data)
        {
            coordinatesDataOutput << std::setw(5) << std::get<0>(entry) << "\t"
                                  << std::setw(8) << std::fixed << std::setprecision(4) << std::get<1>(entry)[0] << "\t"
                                  << std::setw(8) << std::fixed << std::setprecision(4) << std::get<1>(entry)[1] << "\n";
        }
    }
}