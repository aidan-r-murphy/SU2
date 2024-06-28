/*!
 * \file CTurbWAVariable.hpp
 * \brief Declaration of the variables of the WA turbulence model.
 * \author A. Murphy
 * \version 8.0.1 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2024, SU2 Contributors (cf. AUTHORS.md)
 *
 * SU2 is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * SU2 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with SU2. If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once

#include "CTurbVariable.hpp"

/*!
 * \class CTurbWAVariable
 * \brief Main class for defining the variables of the turbulence model.
 * \ingroup Turbulence_Model
 * \author A. Murphy.
 */

class CTurbWAVariable final : public CTurbVariable {
protected:
  su2double C_1komega;
  su2double C_1kepsilon;
  su2double sigma_komega;
  su2double sigma_kepsilon;
  su2double C_mu;
  VectorType f1;    /*!< \brief WA-2017 switching function. */  
  WA_ParsedOptions waParsedOptions;
public:
  /*!
   * \brief Constructor of the class.
   * \param[in] val_R - Turbulent variable value (initialization value).
   * \param[in] val_muT  - The eddy viscosity
   * \param[in] npoint - Number of points/nodes/vertices in the domain.
   * \param[in] ndim - Number of dimensions of the problem.
   * \param[in] nvar - Number of variables of the problem.
   * \param[in] constants - wa model constants
   * \param[in] config - Definition of the particular problem.
   */
  CTurbWAVariable(su2double val_R, su2double val_muT, unsigned long npoint, unsigned long ndim, 
                  unsigned long nvar, const su2double* constants, CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CTurbWAVariable() override = default;

  /*!
   * \brief Set the switching function for the WA model.
   * \param[in] iPoint - Point index.
   * \param[in] Vorticity - Vorticity vector.
   * \param[in] StrainMag_i - Value of strain rate magnitude.
   * \param[in] val_dist - Value of the distance to the wall.
   * \param[in] val_density - Value of the density.
   * \param[in] val_viscosity - Value of the vicosity.
   */
  void SetSwitchingFunc(unsigned long iPoint, const su2double VorticityMag, const su2double StrainMag_i,
                        su2double val_dist, su2double val_density, su2double val_viscosity) override;

  /*!
   * \brief Get the WA model switching function
   * \param[in] iPoint - Point index.
   * \return Value of the f1 switching function
   */
  inline su2double Getf1Switching(unsigned long iPoint) const override { return f1(iPoint); }

};
