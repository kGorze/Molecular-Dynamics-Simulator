//
// Created by konrad_guest on 16/05/2024.
//

#ifndef ABSTRACT_FACTORY_H
#define ABSTRACT_FACTORY_H

#include <vector>
#include <memory>
#include <Eigen/Core>

class Atom2D {
protected:
    Eigen::Vector2d coordinates;
    Eigen::Vector2d velocities;
    Eigen::Vector2d accelerations;

public:
    virtual ~Atom2D() = default;

    virtual Eigen::Vector2d getCoordinates() const;
    virtual Eigen::Vector2d getVelocities() const;
    virtual Eigen::Vector2d getAccelerations() const;

    virtual void setCoordinates(double x, double y);
    virtual void setVelocities(double x, double y);
    virtual void setAccelerations(double x, double y);
};

class Argon2D : public Atom2D {};

class AtomFactory {
public:
    virtual std::shared_ptr<Atom2D> createAtom() = 0;
};

class ArgonFactory : public AtomFactory {
public:
    std::shared_ptr<Atom2D> createAtom() override;
};

#endif // ABSTRACT_FACTORY_H
