//WRITTEN BY: github.com/kgorze

#ifndef DATASAVER_H
#define DATASAVER_H

#include "Headers/simulation/simulation2dap.h"
#include "Headers/core/errorcodes.h"

#include <vector>
#include <string>
#include <tuple>
#include <chrono>
#include <filesystem>
#include <iostream>

#include <Eigen/Dense>


namespace mds {
    class DataSaver
    {
    public:
        virtual ~DataSaver() = default;

        virtual std::string     generateDirectoryName(const std::shared_ptr<Simulation>& simulation)           = 0;
        virtual void            createDirectory(const std::string& dir_name)                                   = 0;
        virtual void            openFile(std::ofstream& file, const std::string& filePath, unsigned int mode)  = 0;
        virtual void            closeFile(std::ofstream& file)                                                 = 0;


        virtual void savePropertiesData(
            const std::string&                                                                                       filename,
            const std::vector<std::tuple<int, double, Eigen::VectorXd, double, double, double, double, double, double>>& data,
            unsigned int                                                                                                 mode,
            const std::shared_ptr<Simulation>&                                                                     simulation,
            std::ofstream&                                                                               propertiesDataOutput
            ) = 0;
        virtual void saveVelocityData(
            const std::string&                            filename,
            const std::vector<std::vector<double>>&           data,
            unsigned int                                      mode,
            const std::shared_ptr<Simulation>&          simulation,
            std::ofstream&          VelocityDistributionDataOutput
            ) = 0;
        virtual void saveCoordinateData(
            const std::string&                                  filename,
            const std::vector<std::tuple<int, Eigen::VectorXd>>&    data,
            unsigned int                                            mode,
            const std::shared_ptr<Simulation>&                simulation,
            std::ofstream&                         coordinatesDataOutput
            ) = 0;
    };

    class DataSaver2D : public DataSaver
    {
    public:

        std::string     generateDirectoryName(const std::shared_ptr<Simulation>& simulation) override;
        void            createDirectory(const std::string& dir_name) override;
        void            openFile(std::ofstream& file, const std::string& filePath, unsigned int mode) override;
        void            closeFile(std::ofstream& file) override;

        void savePropertiesData(
                const std::string&                                                                                       filename,
                const std::vector<std::tuple<int, double, Eigen::VectorXd, double, double, double, double, double, double>>& data,
                unsigned int                                                                                                 mode,
                const std::shared_ptr<Simulation>&                                                                     simulation,
                std::ofstream&                                                                               propertiesDataOutput
                ) override;
        void saveVelocityData(
            const std::string&                                filename,
            const std::vector<std::vector<double>>&               data,
            unsigned int                                          mode,
            const std::shared_ptr<Simulation>&              simulation,
            std::ofstream&              VelocityDistributionDataOutput
            ) override;
        void saveCoordinateData(
            const std::string&                                    filename,
            const std::vector<std::tuple<int, Eigen::VectorXd>>&      data,
            unsigned int                                              mode,
            const std::shared_ptr<Simulation>&                  simulation,
            std::ofstream&                           coordinatesDataOutput
            ) override;
    };

    class DataSaver3D : public DataSaver
    {
    public:
        //void savePropertiesData(const std::string& filename, const std::vector<std::tuple<int, double, Eigen::VectorXd, double, double, double, double, double, double>>& data, unsigned int mode) override;
        //void saveVelocityData(const std::string& filename, const std::vector<std::vector<double>>& data, unsigned int mode) override;
        //void saveCoordinateData(const std::string& filename, const std::vector<std::tuple<int, Eigen::VectorXd>>& data, unsigned int mode) override;
    };
}

#endif //DATASAVER_H