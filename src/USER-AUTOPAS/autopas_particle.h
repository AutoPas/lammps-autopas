#pragma once

#include <vector>

#include "autopas/particles/Particle.h"

namespace LAMMPS_NS {

/**
 * Molecule class for the LJFunctor.
 */
    template <typename floatType = double>
class MoleculeLJLammps final : public autopas::Particle {
    public:
        MoleculeLJLammps() = default;

        /**
         * Constructor of lennard jones molecule with initialization of typeID.
         * @param pos Position of the molecule.
         * @param v Velocitiy of the molecule.
         * @param moleculeId Global Id of the molecule.
         * @param localID Local Id of the molecule.
         * @param typeId TypeId of the molecule.
         */
        explicit MoleculeLJLammps(const std::array<floatType, 3> pos, const std::array<floatType, 3> v, unsigned long moleculeId,
                            int localID, unsigned long typeId = 0)
                : autopas::Particle(pos, v, moleculeId), _typeId(typeId), _localId(localID) {}

        ~MoleculeLJLammps() final = default;

        /**
         * Enums used as ids for accessing and creating a dynamically sized SoA.
         */
        enum AttributeNames : int {
            ptr,
            id,
            posX,
            posY,
            posZ,
            forceX,
            forceY,
            forceZ,
            typeId,
            ownershipState
        };

        /**
         * The type for the SoA storage.
         *
         * @note The attribute owned is of type float but treated as a bool.
         * This means it shall always only take values 0.0 (=false) or 1.0 (=true).
         * The reason for this is the easier use of the value in calculations (See LJFunctor "energyFactor")
         */
        using SoAArraysType = typename autopas::utils::SoAType<
                MoleculeLJLammps<floatType> *, size_t /*id*/, floatType /*x*/, floatType /*y*/, floatType /*z*/,
                floatType /*fx*/, floatType /*fy*/, floatType /*fz*/, size_t /*typeid*/,
                autopas::OwnershipState /*ownershipState*/>::Type;

        /**
         * Non-const getter for the pointer of this object.
         * @tparam attribute Attribute name.
         * @return this.
         */
        template <AttributeNames attribute, std::enable_if_t<attribute == AttributeNames::ptr, bool> = true>
        constexpr typename std::tuple_element<attribute, SoAArraysType>::type::value_type get() {
            return this;
        }
        /**
         * Getter, which allows access to an attribute using the corresponding attribute name (defined in AttributeNames).
         * @tparam attribute Attribute name.
         * @return Value of the requested attribute.
         * @note The value of owned is return as floating point number (true = 1.0, false = 0.0).
         */
        template <AttributeNames attribute, std::enable_if_t<attribute != AttributeNames::ptr, bool> = true>
        constexpr typename std::tuple_element<attribute, SoAArraysType>::type::value_type get() const {
            if constexpr (attribute == AttributeNames::id) {
                return getID();
            } else if constexpr (attribute == AttributeNames::posX) {
                return getR()[0];
            } else if constexpr (attribute == AttributeNames::posY) {
                return getR()[1];
            } else if constexpr (attribute == AttributeNames::posZ) {
                return getR()[2];
            } else if constexpr (attribute == AttributeNames::forceX) {
                return getF()[0];
            } else if constexpr (attribute == AttributeNames::forceY) {
                return getF()[1];
            } else if constexpr (attribute == AttributeNames::forceZ) {
                return getF()[2];
            } else if constexpr (attribute == AttributeNames::typeId) {
                return getTypeId();
            } else if constexpr (attribute == AttributeNames::ownershipState) {
                return this->_ownershipState;
            } else {
                autopas::utils::ExceptionHandler::exception("MoleculeLJLammps::get() unknown attribute {}", attribute);
            }
        }

        /**
         * Setter, which allows set an attribute using the corresponding attribute name (defined in AttributeNames).
         * @tparam attribute Attribute name.
         * @param value New value of the requested attribute.
         * @note The value of owned is extracted from a floating point number (true = 1.0, false = 0.0).
         */
        template <AttributeNames attribute>
        constexpr void set(typename std::tuple_element<attribute, SoAArraysType>::type::value_type value) {
            if constexpr (attribute == AttributeNames::id) {
                setID(value);
            } else if constexpr (attribute == AttributeNames::posX) {
                _r[0] = value;
            } else if constexpr (attribute == AttributeNames::posY) {
                _r[1] = value;
            } else if constexpr (attribute == AttributeNames::posZ) {
                _r[2] = value;
            } else if constexpr (attribute == AttributeNames::forceX) {
                _f[0] = value;
            } else if constexpr (attribute == AttributeNames::forceY) {
                _f[1] = value;
            } else if constexpr (attribute == AttributeNames::forceZ) {
                _f[2] = value;
            } else if constexpr (attribute == AttributeNames::typeId) {
                setTypeId(value);
            } else if constexpr (attribute == AttributeNames::ownershipState) {
                this->_ownershipState = value;
            } else {
                autopas::utils::ExceptionHandler::exception("MoleculeLJLammps::set() unknown attribute {}", attribute);
            }
        }

        /**
         * Get TypeId.
         * @return
         */
        [[nodiscard]] size_t getTypeId() const { return _typeId; }

        /**
         * Set the type id of the Molecule.
         * @param typeId
         */
        void setTypeId(size_t typeId) { _typeId = typeId; }

        [[nodiscard]] int getLocalID() const {
            return _localId;
        }

        void setLocalID(int localId) {
            _localId = localId;
        }

    private:
        /**
         * Particle type id.
         */
        size_t _typeId = 0;

        int _localId = -1;
    };

}  // namespace LAMMPS_NS
