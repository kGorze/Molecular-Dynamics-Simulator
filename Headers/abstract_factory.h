//
// Created by konrad_guest on 16/05/2024.
//

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

    virtual Eigen::VectorXd getCoordinates() const;
    virtual Eigen::VectorXd getVelocities() const;
    virtual Eigen::VectorXd getAccelerations() const;

    virtual void setCoordinates(const Eigen::VectorXd& coords);
    virtual void setVelocities(const Eigen::VectorXd& vels);
    virtual void setAccelerations(const Eigen::VectorXd& accels);
};


class Atom2D : public Atom {};
class Atom3D : public Atom{};

class Argon2D : public Atom {};

class AtomFactory {
public:
    virtual std::shared_ptr<Atom> createAtom() = 0;
};

class ArgonFactory : public AtomFactory {
public:
    std::shared_ptr<Atom> createAtom() override;
};






#endif // ABSTRACT_FACTORY_H
