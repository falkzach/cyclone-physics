/*
 * Implementation file for the particle force generators.
 *
 * Part of the Cyclone physics system.
 *
 * Copyright (c) Icosagon 2003. All Rights Reserved.
 *
 * This software is distributed under licence. Use of this software
 * implies agreement with all terms and conditions of the accompanying
 * software licence.
 */

#include <cyclone/pfgen.h>
#include <math.h>

#include <iostream>

using namespace cyclone;


void ParticleForceRegistry::updateForces(real duration)
{
    Registry::iterator i = registrations.begin();
    for (; i != registrations.end(); i++)
    {
        i->fg->updateForce(i->particle, duration);
    }
}

void ParticleForceRegistry::add(Particle* particle, ParticleForceGenerator *fg)
{
    ParticleForceRegistry::ParticleForceRegistration registration;
    registration.particle = particle;
    registration.fg = fg;
    registrations.push_back(registration);
}

ParticleUplift::ParticleUplift(const Vector3& location, const real radius, const real strength)
: location(location), radius(radius), strength(strength)
{
}

void ParticleUplift::updateForce(Particle* particle, real duration)
{
    real distance = sqrt(pow((location.x - particle->getPosition().x), 2.0) + pow((location.z - particle->getPosition().z), 2.0));
    if (distance <= radius)
    {
        Vector3 up = Vector3(0.0, 1.0, 0.0);

        particle->addForce(up * strength);
    }
}

ParticleGravity::ParticleGravity(const Vector3& gravity)
: gravity(gravity)
{
}

void ParticleGravity::updateForce(Particle* particle, real duration)
{
    // Check that we do not have infinite mass
    if (!particle->hasFiniteMass()) return;

    // Apply the mass-scaled force to the particle
    particle->addForce(gravity * particle->getMass());
}

ParticlePointGravity::ParticlePointGravity(const Vector3& point, const real strength)
: point(point), strength(strength)
{
}

void ParticlePointGravity::updateForce(Particle* particle, real duration)
{
    // Check that we do not have infinite mass
    if (!particle->hasFiniteMass()) return;

    // Calculate direction vector
    Vector3 direction = point - particle->getPosition();

    // Apply the mass-scaled force to the particle
    particle->addForce(direction * strength * particle->getMass());
}

DistanceSquareGravity::DistanceSquareGravity(const Vector3& point, const real mass)
: point(point), mass(mass)
{
}

void DistanceSquareGravity::updateForce(Particle* particle, real duration)
{
    // Check that we do not have infinite mass
    if (!particle->hasFiniteMass()) return;

    // Calculate direction vector
    Vector3 direction = point - particle->getPosition();

    // Calculate distance between points
    real distance = sqrt(pow((point.x - particle->getPosition().x), 2.0) + pow((point.y - particle->getPosition().y), 2.0) + pow((point.z - particle->getPosition().z), 2.0));

    real inverse_distance = 1 / (distance * distance);

    // Normalize direction vector
    // direction = direction * inverse_distance;

    Vector3 force = (direction * c_G * particle->getMass() * mass) * inverse_distance;
    //TODO: scale by 1/distance squared
    //f = G m1 m2 /r^2

    // Apply the mass-scaled force to the particle
    particle->addForce(force);
}

ParticleDrag::ParticleDrag(real k1, real k2)
: k1(k1), k2(k2)
{
}

void ParticleDrag::updateForce(Particle* particle, real duration)
{
    Vector3 force;
    particle->getVelocity(&force);

    // Calculate the total drag coefficient
    real dragCoeff = force.magnitude();
    dragCoeff = k1 * dragCoeff + k2 * dragCoeff * dragCoeff;

    // Calculate the final force and apply it
    force.normalise();
    force *= -dragCoeff;
    particle->addForce(force);
}

ParticleAirBrake::ParticleAirBrake(real k1, real k2, bool enabled)
: k1(k1), k2(k2), enabled(enabled)
{
}

void ParticleAirBrake::toggleEnabled()
{
    enabled = !enabled;
}

void ParticleAirBrake::updateForce(Particle *particle, real duration)
{
    if(enabled) {
        Vector3 force;
        particle->getVelocity(&force);

        // Calculate the total drag coefficient
        real dragCoeff = force.magnitude();
        dragCoeff = k1 * dragCoeff + k2 * dragCoeff * dragCoeff;

        // Calculate the final force and apply it
        force.normalise();
        force *= -dragCoeff;
        particle->addForce(force);
    }
}

ParticleSpring::ParticleSpring(Particle *other, real sc, real rl)
: other(other), springConstant(sc), restLength(rl)
{
}

void ParticleSpring::updateForce(Particle* particle, real duration)
{
    // Calculate the vector of the spring
    Vector3 force;
    particle->getPosition(&force);
    force -= other->getPosition();

    // Calculate the magnitude of the force
    real magnitude = force.magnitude();
    magnitude = real_abs(magnitude - restLength);
    magnitude *= springConstant;

    // Calculate the final force and apply it
    force.normalise();
    force *= -magnitude;
    particle->addForce(force);
}

ParticleBuoyancy::ParticleBuoyancy(real maxDepth,
                                 real volume,
                                 real waterHeight,
                                 real liquidDensity)
:
maxDepth(maxDepth), volume(volume),
waterHeight(waterHeight), liquidDensity(liquidDensity)
{
}

void ParticleBuoyancy::updateForce(Particle* particle, real duration)
{
    // Calculate the submersion depth
    real depth = particle->getPosition().y;

    // Check if we're out of the water
    if (depth >= waterHeight + maxDepth) return;
    Vector3 force(0,0,0);

    // Check if we're at maximum depth
    if (depth <= waterHeight - maxDepth)
    {
        force.y = liquidDensity * volume;
        particle->addForce(force);
        return;
    }

    // Otherwise we are partly submerged
    force.y = liquidDensity * volume *
        (depth - maxDepth - waterHeight) / 2 * maxDepth;
    particle->addForce(force);
}

ParticleBungee::ParticleBungee(Particle *other, real sc, real rl)
: other(other), springConstant(sc), restLength(rl)
{
}

void ParticleBungee::updateForce(Particle* particle, real duration)
{
    // Calculate the vector of the spring
    Vector3 force;
    particle->getPosition(&force);
    force -= other->getPosition();

    // Check if the bungee is compressed
    real magnitude = force.magnitude();
    if (magnitude <= restLength) return;

    // Calculate the magnitude of the force
    magnitude = springConstant * (restLength - magnitude);

    // Calculate the final force and apply it
    force.normalise();
    force *= -magnitude;
    particle->addForce(force);
}

ParticleFakeSpring::ParticleFakeSpring(Vector3 *anchor, real sc, real d)
: anchor(anchor), springConstant(sc), damping(d)
{
}

void ParticleFakeSpring::updateForce(Particle* particle, real duration)
{
    // Check that we do not have infinite mass
    if (!particle->hasFiniteMass()) return;

    // Calculate the relative position of the particle to the anchor
    Vector3 position;
    particle->getPosition(&position);
    position -= *anchor;

    // Calculate the constants and check they are in bounds.
    real gamma = 0.5f * real_sqrt(4 * springConstant - damping*damping);
    if (gamma == 0.0f) return;
    Vector3 c = position * (damping / (2.0f * gamma)) +
        particle->getVelocity() * (1.0f / gamma);

    // Calculate the target position
    Vector3 target = position * real_cos(gamma * duration) +
        c * real_sin(gamma * duration);
    target *= real_exp(-0.5f * duration * damping);

    // Calculate the resulting acceleration and therefore the force
    Vector3 accel = (target - position) * ((real)1.0 / (duration*duration)) -
        particle->getVelocity() * ((real)1.0/duration);
    particle->addForce(accel * particle->getMass());
}

ParticleAnchoredSpring::ParticleAnchoredSpring()
{
}

ParticleAnchoredSpring::ParticleAnchoredSpring(Vector3 *anchor,
                                               real sc, real rl)
: anchor(anchor), springConstant(sc), restLength(rl)
{
}

void ParticleAnchoredSpring::init(Vector3 *anchor, real springConstant,
                                  real restLength)
{
    ParticleAnchoredSpring::anchor = anchor;
    ParticleAnchoredSpring::springConstant = springConstant;
    ParticleAnchoredSpring::restLength = restLength;
}

void ParticleAnchoredBungee::updateForce(Particle* particle, real duration)
{
    // Calculate the vector of the spring
    Vector3 force;
    particle->getPosition(&force);
    force -= *anchor;

    // Calculate the magnitude of the force
    real magnitude = force.magnitude();
    if (magnitude < restLength) return;

    magnitude = magnitude - restLength;
    magnitude *= springConstant;

    // Calculate the final force and apply it
    force.normalise();
    force *= -magnitude;
    particle->addForce(force);
}

void ParticleAnchoredSpring::updateForce(Particle* particle, real duration)
{
    // Calculate the vector of the spring
    Vector3 force;
    particle->getPosition(&force);
    force -= *anchor;

    // Calculate the magnitude of the force
    real magnitude = force.magnitude();
    magnitude = (restLength - magnitude) * springConstant;

    // Calculate the final force and apply it
    force.normalise();
    force *= magnitude;
    particle->addForce(force);
}

ParticleLimitedAnchoredSpring::ParticleLimitedAnchoredSpring(Vector3 *anchor, real sc, real rl, real ll)
: anchor(anchor), springConstant(sc), restLength(rl), limitLength(ll)
{
}

void ParticleLimitedAnchoredSpring::updateForce(Particle* particle, real duration)
{
    // Calculate the vector of the spring
    Vector3 force;
    particle->getPosition(&force);
    force -= *anchor;

    // Calculate the magnitude of the force
    real magnitude = force.magnitude();
    if (magnitude < restLength) return;

    real scaler = 1.0;
    if (magnitude > limitLength)
    {
        scaler = limitLength/magnitude;
    }

    magnitude = magnitude - restLength;
    magnitude *= springConstant * scaler;

    // Calculate the final force and apply it
    force.normalise();
    force *= -magnitude;
    particle->addForce(force);
}


ParticleLighter::ParticleLighter(real volume, real airDensity)
        :
        volume(volume), airDensity(airDensity)
{
}

void ParticleLighter::updateForce(Particle* particle, real duration)
{
    // Calculate the submersion depth
    real height = particle->getPosition().y;

    Vector3 force(0,0,0);
    if (height <= 0) return;
    force.y = (particle->getMass() / volume) / (airDensity * height);
    if (force.y >= 1.0) return;
    particle->addForce(force);
    return;
}