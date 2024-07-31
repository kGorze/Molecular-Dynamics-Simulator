//WRITTEN BY: github.com/kgorze

#ifndef ABSTRACT_FACTORY_H
#define ABSTRACT_FACTORY_H

#include <vector>
#include <memory>

#include <Eigen/Core>

class Atom {
protected:
    Eigen::VectorXd coordinates;
    Eigen::VectorXd velocities;
    Eigen::VectorXd accelerations;

public:
    virtual ~Atom() = default;

    [[nodiscard]] virtual Eigen::VectorXd getCoordinates() const = 0;
    [[nodiscard]] virtual Eigen::VectorXd getVelocities() const = 0;
    [[nodiscard]] virtual Eigen::VectorXd getAccelerations() const = 0;

    virtual void setCoordinates(const Eigen::VectorXd& coords) = 0;
    virtual void setVelocities(const Eigen::VectorXd& vels) = 0;
    virtual void setAccelerations(const Eigen::VectorXd& accels) = 0;
};


class Atom2D : public Atom {
protected:
    Eigen::VectorXd coordinates = Eigen::VectorXd(2);
    Eigen::VectorXd velocities = Eigen::VectorXd(2);
    Eigen::VectorXd accelerations = Eigen::VectorXd(2);
public:
    virtual ~Atom2D() = default;

    [[nodiscard]] Eigen::VectorXd getCoordinates() const override;
    [[nodiscard]] Eigen::VectorXd getVelocities() const override;
    [[nodiscard]] Eigen::VectorXd getAccelerations() const override;

    void setCoordinates(const Eigen::VectorXd& coords) override;
    void setVelocities(const Eigen::VectorXd& vels) override;
    void setAccelerations(const Eigen::VectorXd& accels) override;


};
class Atom3D : public Atom{};

class AtomFactory {
public:
    virtual std::shared_ptr<Atom2D> createAtom() = 0;
};

class ArgonFactory : public AtomFactory {
public:
    std::shared_ptr<Atom2D> createAtom() override;
};






#endif // ABSTRACT_FACTORY_H
