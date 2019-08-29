/*!
 * \file CMeshVariable.hpp
 * \brief Declaration and inlines of the class
 *        to define the variables of the mesh movement.
 * \author Ruben Sanchez
 * \version 6.2.0 "Falcon"
 *
 * The current SU2 release has been coordinated by the
 * SU2 International Developers Society <www.su2devsociety.org>
 * with selected contributions from the open-source community.
 *
 * The main research teams contributing to the current release are:
 *  - Prof. Juan J. Alonso's group at Stanford University.
 *  - Prof. Piero Colonna's group at Delft University of Technology.
 *  - Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *  - Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *  - Prof. Rafael Palacios' group at Imperial College London.
 *  - Prof. Vincent Terrapon's group at the University of Liege.
 *  - Prof. Edwin van der Weide's group at the University of Twente.
 *  - Lab. of New Concepts in Aeronautics at Tech. Institute of Aeronautics.
 *
 * Copyright 2012-2019, Francisco D. Palacios, Thomas D. Economon,
 *                      Tim Albring, and the SU2 contributors.
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

#include "CVariable.hpp"

class CMeshVariable : public CVariable {
protected:

  Vec_t WallDistance;  /*!< \brief Store the wall distance in reference coordinates. */
  Mat_t Mesh_Coord;    /*!< \brief Store the reference coordinates of the mesh. */

  /*!
   * \brief Constructor of the class.
   * \param[in] npoint - Number of points/nodes/vertices in the domain.
   * \param[in] ndim - Number of dimensions of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CMeshVariable(Idx_t npoint, Idx_t ndim, CConfig *config);

public:
  /*!
   * \brief Destructor of the class.
   */
  ~CMeshVariable() = default;

  /*!
   * \brief Get the value of the undeformed coordinates.
   * \param[in] iDim - Index of Mesh_Coord[nDim]
   * \return Value of the original coordinate iDim.
   */
  inline su2double GetMesh_Coord(Idx_t iPoint, Idx_t iDim) const final { return Mesh_Coord(iPoint,iDim); }

  /*!
   * \brief Get the undeformed coordinates.
   * \return Pointer to the reference coordinates.
   */
  inline const su2double *GetMesh_Coord(Idx_t iPoint) const final { return Mesh_Coord[iPoint]; }

  /*!
   * \brief Set the value of the undeformed coordinates.
   * \param[in] iDim - Index of Mesh_Coord[nDim]
   * \param[in] val_coord - Value of Mesh_Coord[nDim]
   */
  inline void SetMesh_Coord(Idx_t iPoint, Idx_t iDim, su2double val_coord) final {
    Mesh_Coord(iPoint,iDim) = val_coord;
  }

  /*!
   * \brief Get the value of the wall distance in reference coordinates.
   * \param[in] iDim - Index of Mesh_Coord[nDim]
   * \return Value of the wall distance in reference coordinates.
   */
  inline su2double GetWallDistance(Idx_t iPoint) const final { return WallDistance(iPoint); }

  /*!
   * \brief Set the value of the wall distance in reference coordinates.
   * \param[in] val_dist - Value of wall distance.
   */
  inline void SetWallDistance(Idx_t iPoint, su2double val_dist) final { WallDistance(iPoint) = val_dist; }

  /*!
   * \brief Register the reference coordinates of the mesh.
   * \param[in] input - Defines whether we are registering the variable as input or as output.
   */
  inline void Register_MeshCoord(Idx_t iPoint, bool input) final {
    if (input) {
      for (Idx_t iDim = 0; iDim < nDim; iDim++)
        AD::RegisterInput(Mesh_Coord(iPoint,iDim));
    }
    else {
      for (Idx_t iDim = 0; iDim < nDim; iDim++)
        AD::RegisterOutput(Mesh_Coord(iPoint,iDim));
    }
  }

  /*!
   * \brief Recover the value of the adjoint of the mesh coordinates.
   */
  inline void GetAdjoint_MeshCoord(Idx_t iPoint, su2double *adj_mesh) const final {
    for (Idx_t iDim = 0; iDim < nDim; iDim++)
      adj_mesh[iDim] = SU2_TYPE::GetDerivative(Mesh_Coord(iPoint,iDim));
  }

};
