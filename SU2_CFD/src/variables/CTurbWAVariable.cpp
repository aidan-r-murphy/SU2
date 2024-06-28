/*!
 * \file CTurbWAVariable.cpp
 * \brief Definition of the solution fields.
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


#include "../../include/variables/CTurbWAVariable.hpp"


CTurbWAVariable::CTurbWAVariable(su2double val_R, su2double val_muT, unsigned long npoint, unsigned long ndim, 
                                 unsigned long nvar, const su2double* constants, CConfig *config) :
                 CTurbVariable(npoint, ndim, nvar, config) {

  waParsedOptions = config->GetWAParsedOptions();

  for(unsigned long iPoint=0; iPoint<nPoint; ++iPoint)
  {
    Solution(iPoint,0) = val_R;
  }

  Solution_Old = Solution;

  C_1komega = constants[0];
  C_1kepsilon = constants[1];
  sigma_komega = constants[2];
  sigma_kepsilon = constants[3];
  C_mu = constants[9];

  f1.resize(nPoint) = su2double(1.0);

  muT.resize(nPoint) = val_muT;

  nAuxVar = 1;
  AuxVar.resize(nPoint,nAuxVar) = su2double(0.0);
  Grad_AuxVar.resize(nPoint,nAuxVar,nDim);
}

void CTurbWAVariable::SetSwitchingFunc(unsigned long iPoint, const su2double VorticityMag, const su2double StrainMag_i,
                                       su2double val_dist, su2double val_density, su2double val_viscosity) {
  su2double arg1;
  
  su2double W, S;
  su2double numerator, denominator;
  su2double k, omega, eta;
  su2double EPS = 1.0e-16; /*!< \brief Error scale for mean strain. */

  AD::StartPreacc();
  AD::SetPreaccIn(val_viscosity);
  AD::SetPreaccIn(val_density);
  AD::SetPreaccIn(val_dist);
  AD::SetPreaccIn(Solution[iPoint], nVar);
  /*--- Eddy viscosity ---*/
  AD::SetPreaccIn(muT(iPoint));

  W = VorticityMag;
  S = max(StrainMag_i, EPS);

  /* --- f1 --- */
  if (waParsedOptions.version == WA_OPTIONS::V2018) { // WA WDF 2018
    k = muT(iPoint)/val_density * S / sqrt(C_mu);
    omega = S / sqrt(C_mu);
    eta = S * max(1.0,abs(W/S));

    arg1 = (val_viscosity/val_density + Solution(iPoint,0))/2.0 * pow(eta,2.0)/(C_mu * k * omega);

    f1(iPoint) = tanh(pow(arg1,4.0));
  } else { // WA 2017 or 2017m
    numerator = 1.0 + val_dist*sqrt(Solution(iPoint,0)*S)*val_density/val_viscosity;
    denominator = 1.0 + pow((max(val_dist*sqrt(Solution(iPoint,0)*S),1.5*Solution(iPoint,0))) / 
                            (20.0*val_viscosity/val_density),2.0);
    arg1 = numerator/denominator;
    f1(iPoint) = min(tanh(pow(arg1,4.0)),0.9);
  }

  AD::SetPreaccOut(f1(iPoint));
  AD::EndPreacc();
}
