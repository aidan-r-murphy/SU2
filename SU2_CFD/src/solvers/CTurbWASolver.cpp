/*!
 * \file CTurbWASolver.cpp
 * \brief Main subroutines of CTurbWASolver class
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

#include "../../include/solvers/CTurbWASolver.hpp"
#include "../../include/variables/CTurbWAVariable.hpp"
#include "../../include/variables/CFlowVariable.hpp"
#include "../../../Common/include/parallelization/omp_structure.hpp"
#include "../../../Common/include/toolboxes/geometry_toolbox.hpp"


CTurbWASolver::CTurbWASolver(CGeometry *geometry, CConfig *config, unsigned short iMesh, CFluidModel* FluidModel)
             : CTurbSolver(geometry, config, false) {
  unsigned long iPoint;

  bool multizone = config->GetMultizone_Problem();
  waParsedOptions = config->GetWAParsedOptions();

  /*--- Dimension of the problem --> dependent of the turbulent model ---*/

  nVar = 1;
  nPrimVar = 1;
  nPoint = geometry->GetnPoint();
  nPointDomain = geometry->GetnPointDomain();

  /*--- Initialize nVarGrad for deallocation ---*/

  nVarGrad = nVar;

  /*--- Define geometry constants in the solver structure ---*/

  nDim = geometry->GetnDim();

  /*--- Single grid simulation ---*/

  if (iMesh == MESH_0 || config->GetMGCycle() == FULLMG_CYCLE) {

    /*--- Define some auxiliar vector related with the residual ---*/

    Residual_RMS.resize(nVar,0.0);
    Residual_Max.resize(nVar,0.0);
    Point_Max.resize(nVar,0);
    Point_Max_Coord.resize(nVar,nDim) = su2double(0.0);

    /*--- Initialization of the structure of the whole Jacobian ---*/

    if (rank == MASTER_NODE) cout << "Initialize Jacobian structure (WA model)." << endl;
    Jacobian.Initialize(nPoint, nPointDomain, nVar, nVar, true, geometry, config, ReducerStrategy);
    LinSysSol.Initialize(nPoint, nPointDomain, nVar, 0.0);
    LinSysRes.Initialize(nPoint, nPointDomain, nVar, 0.0);
    System.SetxIsZero(true);

    if (ReducerStrategy)
      EdgeFluxes.Initialize(geometry->GetnEdge(), geometry->GetnEdge(), nVar, nullptr);

    /*--- Initialize the BGS residuals in multizone problems. ---*/
    if (multizone){
      Residual_BGS.resize(nVar,0.0);
      Residual_Max_BGS.resize(nVar,0.0);
      Point_Max_BGS.resize(nVar,0);
      Point_Max_Coord_BGS.resize(nVar,nDim) = su2double(0.0);
    }

  }

  /* --- Initialize value for WA model constants --- */
  constants[0] = 0.0829;    // C_1komega
  constants[1] = 0.1127;    // C_1kepsilon
  constants[2] = 0.72;      // sigma_komega
  constants[3] = 1.0;       // sigma_kepsilon
  constants[4] = 0.41;      // kappa
  constants[5] = 8.54;      // C_omega

  constants[6] = constants[0]/(pow(constants[4],2.0)) + constants[2];   // C_2komega
  constants[7] = constants[1]/(pow(constants[4],2.0)) + constants[3];   // C_2kepsilon

  constants[8] = 8.0;       // C_m

  /*--- Read farfield conditions from config ---*/
  su2double Density_Inf, LaminarViscosity_Inf, Factor_R_Inf, muT_Inf;

  Density_Inf   = config->GetDensity_FreeStreamND();
  LaminarViscosity_Inf = config->GetViscosity_FreeStreamND();

  /*--- Factor_R_Inf in [3.0, 5.0] ---*/
  Factor_R_Inf = config->GetRFactor_FreeStream();
  su2double R_Inf  = Factor_R_Inf*LaminarViscosity_Inf/Density_Inf;

  Solution_Inf[0] = R_Inf;

  /*--- Eddy viscosity at infinity ---*/
  su2double Xi, Xi_3, fmu, Cw_3 = pow(constants[5],3.0);
  Xi = R_Inf/LaminarViscosity_Inf*Density_Inf;
  Xi_3 = Xi*Xi*Xi;
  fmu = Xi_3/(Xi_3+Cw_3);
  muT_Inf = Density_Inf*fmu*R_Inf;

  /*--- Initialize the solution to the far-field state everywhere. ---*/

  nodes = new CTurbWAVariable(R_Inf, muT_Inf, nPoint, nDim, nVar, constants, config);
  SetBaseClassPointerToNodes();

  /*--- MPI solution ---*/

  InitiateComms(geometry, config, SOLUTION_EDDY);
  CompleteComms(geometry, config, SOLUTION_EDDY);

  /*--- Initializate quantities for SlidingMesh Interface ---*/

  SlidingState.resize(nMarker);
  SlidingStateNodes.resize(nMarker);

  for (unsigned long iMarker = 0; iMarker < nMarker; iMarker++) {
    if (config->GetMarker_All_KindBC(iMarker) == FLUID_INTERFACE) {
      SlidingState[iMarker].resize(nVertex[iMarker], nPrimVar+1) = nullptr;
      SlidingStateNodes[iMarker].resize(nVertex[iMarker],0);
    }
  }

  /*-- Allocation of inlets has to happen in derived classes (not CTurbSolver),
   * due to arbitrary number of turbulence variables ---*/

  Inlet_TurbVars.resize(nMarker);
  for (unsigned long iMarker = 0; iMarker < nMarker; iMarker++)
    Inlet_TurbVars[iMarker].resize(nVertex[iMarker],nVar) = R_Inf;

  /*--- Store the initial CFL number for all grid points. ---*/

  const su2double CFL = config->GetCFL(MGLevel)*config->GetCFLRedCoeff_Turb();
  for (iPoint = 0; iPoint < nPoint; iPoint++) {
    nodes->SetLocalCFL(iPoint, CFL);
  }
  Min_CFL_Local = CFL;
  Max_CFL_Local = CFL;
  Avg_CFL_Local = CFL;

  /*--- Add the solver name. ---*/
  SolverName = "WA";

}

void CTurbWASolver::Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config,
        unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem, bool Output) {
  SU2_OMP_SAFE_GLOBAL_ACCESS(config->SetGlobalParam(config->GetKind_Solver(), RunTime_EqSystem);)

  /*--- Clear Residual and Jacobian. Upwind second order reconstruction and gradients ---*/
  CommonPreprocessing(geometry, config, Output);

}

void CTurbWASolver::Postprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh) {

  const su2double Cw_3 = pow(constants[5],3.0);

  /*--- Compute eddy viscosity ---*/

  AD::StartNoSharedReading();

  auto* flowNodes = su2staticcast_p<CFlowVariable*>(solver_container[FLOW_SOL]->GetNodes());

  SU2_OMP_FOR_STAT(omp_chunk_size)
  for (unsigned long iPoint = 0; iPoint < nPoint; iPoint++) {

    const su2double rho = flowNodes->GetDensity(iPoint);
    const su2double mu = flowNodes->GetLaminarViscosity(iPoint);

    const su2double nu = mu/rho;
    const su2double R = nodes->GetSolution(iPoint,0);
    const su2double dist = geometry->nodes->GetWall_Distance(iPoint);

    const su2double VorticityMag = GeometryToolbox::Norm(3, flowNodes->GetVorticity(iPoint));
    const su2double StrainMag = max(flowNodes->GetStrainMag(iPoint), 1e-16);
    nodes->SetSwitchingFunc(iPoint, VorticityMag, StrainMag, dist, rho, mu);

    const su2double Xi = R/nu;
    const su2double Xi_3 = Xi*Xi*Xi;
    const su2double fmu  = Xi_3/(Xi_3+Cw_3);

    const su2double muT = rho*fmu*R;

    nodes->SetmuT(iPoint,muT);

  }
  END_SU2_OMP_FOR


//  TODO: Implement transition model for WA (WA-AT)
  /*--- Compute turbulence index ---*/
  if (config->GetKind_Trans_Model() != TURB_TRANS_MODEL::NONE) {
    SU2_MPI::Error("Transition model not implemented for WA model", CURRENT_FUNCTION);
  }

  AD::EndNoSharedReading();
}

void CTurbWASolver::Viscous_Residual(unsigned long iEdge, CGeometry* geometry, CSolver** solver_container,
                                     CNumerics* numerics, CConfig* config) {

  /*--- Define an object to set solver specific numerics contribution. ---*/
  auto SolverSpecificNumerics = [&](unsigned long iPoint, unsigned long jPoint) {
    /*--- WA model switching function (only WA) ---*/
    numerics->Setf1Switching(nodes->Getf1Switching(iPoint), nodes->Getf1Switching(jPoint));
    /*--- Roughness heights. ---*/
    numerics->SetRoughness(geometry->nodes->GetRoughnessHeight(iPoint), geometry->nodes->GetRoughnessHeight(jPoint));
  };

  /*--- Now instantiate the generic implementation with the functor above. ---*/

  Viscous_Residual_impl(SolverSpecificNumerics, iEdge, geometry, solver_container, numerics, config);
}

void CTurbWASolver::Source_Residual(CGeometry *geometry, CSolver **solver_container,
                                    CNumerics **numerics_container, CConfig *config, unsigned short iMesh) {

  const bool implicit = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);

  auto* flowNodes = su2staticcast_p<CFlowVariable*>(solver_container[FLOW_SOL]->GetNodes());

  /*--- Pick one numerics object per thread. ---*/
  auto* numerics = numerics_container[SOURCE_FIRST_TERM + omp_get_thread_num()*MAX_TERMS];

  AD::StartNoSharedReading();
  
  /* --- Set strain mag as auxiliary variable --- */
  SU2_OMP_FOR_STAT(omp_chunk_size)
  for (unsigned long iPoint = 0; iPoint < nPoint; iPoint++) {
    const su2double StrainMag = max(flowNodes->GetStrainMag(iPoint), 1e-16);
    nodes->SetAuxVar(iPoint, 0, StrainMag);
  }
  END_SU2_OMP_FOR
  
  /*--- calculate the gradient of the vorticity magnitude (AuxVarGradient) ---*/

  if (config->GetKind_Gradient_Method() == GREEN_GAUSS) SetAuxVar_Gradient_GG(geometry, config);
  if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) SetAuxVar_Gradient_LS(geometry, config);

  /*--- Loop over all points. ---*/

  SU2_OMP_FOR_DYN(omp_chunk_size)
  for (unsigned long iPoint = 0; iPoint < nPointDomain; iPoint++) {

    /*--- Conservative variables w/o reconstruction ---*/

    numerics->SetPrimitive(flowNodes->GetPrimitive(iPoint), nullptr);

    /*--- Gradient of the primitive and conservative variables ---*/

    numerics->SetPrimVarGradient(flowNodes->GetGradient_Primitive(iPoint), nullptr);

    /*--- Set vorticity and strain rate magnitude ---*/

    numerics->SetVorticity(flowNodes->GetVorticity(iPoint), nullptr);

    numerics->SetStrainMag(flowNodes->GetStrainMag(iPoint), 0.0);

    /*--- Turbulent variables w/o reconstruction, and its gradient ---*/

    numerics->SetScalarVar(nodes->GetSolution(iPoint), nullptr);
    numerics->SetScalarVarGradient(nodes->GetGradient(iPoint), nullptr);

    /*--- Set volume ---*/

    numerics->SetVolume(geometry->nodes->GetVolume(iPoint));

    /*--- Set distance to the surface ---*/

    numerics->SetDistance(geometry->nodes->GetWall_Distance(iPoint), 0.0);

    /*--- WA model switching function ---*/

    numerics->Setf1Switching(nodes->Getf1Switching(iPoint), 0.0);

    /*--- calculate the gradient of the vorticity magnitude (AuxVarGradient) ---*/
    
    numerics->SetAuxVarGrad(nodes->GetAuxVarGradient(iPoint), nullptr);

    /*--- Compute the source term ---*/

    auto residual = numerics->ComputeResidual(config);

    /*--- Subtract residual and the Jacobian ---*/

    LinSysRes.SubtractBlock(iPoint, residual);
    if (implicit) Jacobian.SubtractBlock2Diag(iPoint, residual.jacobian_i);

  }
  END_SU2_OMP_FOR

  AD::EndNoSharedReading();

}

void CTurbWASolver::Source_Template(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                                    CConfig *config, unsigned short iMesh) {
}

void CTurbWASolver::BC_HeatFlux_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics,
                                     CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {

  const bool implicit = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);
  
  // TODO: Implement wall roughness for WA model
  bool rough_wall = false;
  string Marker_Tag = config->GetMarker_All_TagBound(val_marker);
  WALL_TYPE WallType; su2double Roughness_Height;
  tie(WallType, Roughness_Height) = config->GetWallRoughnessProperties(Marker_Tag);
  if (WallType == WALL_TYPE::ROUGH) rough_wall = true;

  /*--- Evaluate R at the closest point to the surface using the wall functions. ---*/

  if (config->GetWall_Functions()) {
    SU2_OMP_SAFE_GLOBAL_ACCESS(SetTurbVars_WF(geometry, solver_container, config, val_marker);)
    return;
  }

  SU2_OMP_FOR_STAT(OMP_MIN_SIZE)
  for (auto iVertex = 0u; iVertex < geometry->nVertex[val_marker]; iVertex++) {

    const auto iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/
    if (geometry->nodes->GetDomain(iPoint)) {
      if (!rough_wall) { // smooth wall
        for (auto iVar = 0u; iVar < nVar; iVar++) {
          nodes->SetSolution_Old(iPoint,iVar,0.0);
          nodes->SetSolution(iPoint,iVar,0.0);
        }

        LinSysRes.SetBlock_Zero(iPoint);

        /*--- Includes 1 in the diagonal ---*/

        if (implicit) Jacobian.DeleteValsRowi(iPoint);
       } else {
         // TODO: Implement wall roughness for WA
         SU2_MPI::Error("Wall roughness not implemented for WA model", CURRENT_FUNCTION);
      }
    }
  }
  END_SU2_OMP_FOR
}

void CTurbWASolver::SetTurbVars_WF(CGeometry *geometry, CSolver **solver_container,
                                  const CConfig *config, unsigned short val_marker) {

  const bool implicit = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);

  /*--- We use a very high max nr of iterations, but we only need this the first couple of iterations ---*/
  const unsigned short max_iter = config->GetwallModel_MaxIter();

  /* --- tolerance has LARGE impact on convergence, do not increase this value! --- */
  const su2double tol = 1e-12;

  /*--- Typical constants from boundary layer theory ---*/

  const su2double Cw_3 = 8.54*8.54*8.54;

  CVariable* flow_nodes = solver_container[FLOW_SOL]->GetNodes();

  /*--- Loop over all of the vertices on this boundary marker ---*/

  for (auto iVertex = 0u; iVertex < geometry->nVertex[val_marker]; iVertex++) {

    const auto iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    const auto iPoint_Neighbor = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();

    /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/

    if (geometry->nodes->GetDomain(iPoint_Neighbor)) {

      su2double Y_Plus = solver_container[FLOW_SOL]->GetYPlus(val_marker, iVertex);

      /*--- Do not use wall model at the ipoint when y+ < "limit", use zero flux (Neumann) conditions. ---*/

      if (Y_Plus < config->GetwallModel_MinYPlus()) continue;

      su2double Lam_Visc_Normal = flow_nodes->GetLaminarViscosity(iPoint_Neighbor);
      su2double Density_Normal = flow_nodes->GetDensity(iPoint_Neighbor);
      su2double Kin_Visc_Normal = Lam_Visc_Normal/Density_Normal;

      su2double Eddy_Visc = solver_container[FLOW_SOL]->GetEddyViscWall(val_marker, iVertex);

      /*--- Solve for the new value of R given the eddy viscosity and using a Newton method ---*/

      // start with positive value of R_old
      su2double R = 0.0;
      su2double R_old = nodes->GetSolution(iPoint,0);

      unsigned short counter = 0;
      su2double diff = 1.0;
      su2double relax = config->GetwallModel_RelFac();
      while (diff > tol) {
        // note the error in Nichols and Nelson
        su2double func = pow(R_old,4) - (Eddy_Visc/Density_Normal)*(pow(R_old,3) + pow(Kin_Visc_Normal,3)*Cw_3);
        su2double func_prim = 4.0 * pow(R_old,3) - 3.0*(Eddy_Visc/Density_Normal)*pow(R_old,2);

        // damped Newton method
        R = R_old - relax*(func/func_prim);

        diff = fabs(R-R_old);
        R_old = R;

        // sometimes we get negative values when the solution has not converged yet, we just reset the R in that case.
        if (R_old<tol) {
          relax /= 2.0;
          R_old = nodes->GetSolution(iPoint,0)/relax;
        }

        counter++;
        if (counter > max_iter) break;
      }

      nodes->SetSolution_Old(iPoint_Neighbor, &R);
      LinSysRes.SetBlock_Zero(iPoint_Neighbor);

      /*--- includes 1 in the diagonal ---*/

      if (implicit) Jacobian.DeleteValsRowi(iPoint_Neighbor);
    }
  }
}

void CTurbWASolver::BC_Isothermal_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics,
                                       CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {

  BC_HeatFlux_Wall(geometry, solver_container, conv_numerics, visc_numerics, config, val_marker);

}

void CTurbWASolver::BC_Inlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics,
                             CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {

  const bool implicit = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);

  /*--- Loop over all the vertices on this boundary marker ---*/

  SU2_OMP_FOR_STAT(OMP_MIN_SIZE)
  for (auto iVertex = 0u; iVertex < geometry->nVertex[val_marker]; iVertex++) {

    const auto iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    /*--- Check if the node belongs to the domain (i.e., not a halo node) ---*/

    if (geometry->nodes->GetDomain(iPoint)) {

      /*--- Normal vector for this vertex (negate for outward convention) ---*/

      su2double Normal[MAXNDIM] = {0.0};
      for (auto iDim = 0u; iDim < nDim; iDim++)
        Normal[iDim] = -geometry->vertex[val_marker][iVertex]->GetNormal(iDim);
      conv_numerics->SetNormal(Normal);

      /*--- Allocate the value at the inlet ---*/

      auto V_inlet = solver_container[FLOW_SOL]->GetCharacPrimVar(val_marker, iVertex);

      /*--- Retrieve solution at the farfield boundary node ---*/

      auto V_domain = solver_container[FLOW_SOL]->GetNodes()->GetPrimitive(iPoint);

      /*--- Set various quantities in the solver class ---*/

      conv_numerics->SetPrimitive(V_domain, V_inlet);

      su2double Inlet_Vars[MAXNVAR];
      Inlet_Vars[0] = Inlet_TurbVars[val_marker][iVertex][0];
      if (config->GetInlet_Profile_From_File()) {
         /*--- Non-dimensionalize Inlet_TurbVars if Inlet-Files are used. ---*/
         Inlet_Vars[0] *= config->GetDensity_Ref() / config->GetViscosity_Ref();
      } else {
         /*--- Obtain fluid model for computing the R to impose at the inlet boundary. ---*/
         CFluidModel* FluidModel = solver_container[FLOW_SOL]->GetFluidModel();

         /*--- Obtain density and laminar viscosity at inlet boundary node ---*/

         su2double Density_Inlet;
         if (config->GetKind_Regime() == ENUM_REGIME::COMPRESSIBLE) {
           Density_Inlet = V_inlet[prim_idx.Density()];
           FluidModel->SetTDState_Prho(V_inlet[prim_idx.Pressure()], Density_Inlet);
         } else {
           const su2double* Scalar_Inlet = nullptr;
           if (config->GetKind_Species_Model() != SPECIES_MODEL::NONE) {
            Scalar_Inlet = config->GetInlet_SpeciesVal(config->GetMarker_All_TagBound(val_marker));
           }
           FluidModel->SetTDState_T(V_inlet[prim_idx.Temperature()], Scalar_Inlet);
           Density_Inlet = FluidModel->GetDensity();
         }
         const su2double Laminar_Viscosity_Inlet = FluidModel->GetLaminarViscosity();
         const su2double* Turb_Properties = config->GetInlet_TurbVal(config->GetMarker_All_TagBound(val_marker));
         const su2double R_Factor = Turb_Properties[0];
         Inlet_Vars[0] = R_Factor * Laminar_Viscosity_Inlet / Density_Inlet;
      }

      /*--- Load the inlet turbulence variable (uniform by default). ---*/

      conv_numerics->SetScalarVar(nodes->GetSolution(iPoint), Inlet_Vars);

      /*--- Set various other quantities in the conv_numerics class ---*/

      if (dynamic_grid)
        conv_numerics->SetGridVel(geometry->nodes->GetGridVel(iPoint),
                                  geometry->nodes->GetGridVel(iPoint));

      if (conv_numerics->GetBoundedScalar()) {
        const su2double* velocity = &V_inlet[prim_idx.Velocity()];
        const su2double density = solver_container[FLOW_SOL]->GetNodes()->GetDensity(iPoint);
        conv_numerics->SetMassFlux(BoundedScalarBCFlux(iPoint, implicit, density, velocity, Normal));
      }

      /*--- Compute the residual using an upwind scheme ---*/

      auto residual = conv_numerics->ComputeResidual(config);
      LinSysRes.AddBlock(iPoint, residual);

      /*--- Jacobian contribution for implicit integration ---*/

      if (implicit) Jacobian.AddBlock2Diag(iPoint, residual.jacobian_i);

//      /*--- Viscous contribution, commented out because serious convergence problems ---*/
//
//      su2double Coord_Reflected[MAXNDIM];
//      GeometryToolbox::PointPointReflect(nDim, geometry->nodes->GetCoord(Point_Normal),
//                                               geometry->nodes->GetCoord(iPoint), Coord_Reflected);
//      visc_numerics->SetCoord(geometry->nodes->GetCoord(iPoint), Coord_Reflected);
//      visc_numerics->SetNormal(Normal);
//
//      /*--- Conservative variables w/o reconstruction ---*/
//
//      visc_numerics->SetPrimitive(V_domain, V_inlet);
//
//      /*--- Turbulent variables w/o reconstruction, and its gradients ---*/
//
//      visc_numerics->SetScalarVar(Solution_i, Solution_j);
//      visc_numerics->SetScalarVarGradient(node[iPoint]->GetGradient(), node[iPoint]->GetGradient());
//
//      /*--- WA model switching function ---*/
//
//      visc_numerics->Setf1Switching(node[iPoint]->Getf1Switching(), node[iPoint]->Getf1Switching());
//
//      /*--- Compute residual, and Jacobians ---*/
//
//      auto residual = visc_numerics->ComputeResidual(config);
//
//      /*--- Subtract residual, and update Jacobians ---*/
//
//      LinSysRes.SubtractBlock(iPoint, residual);
//      Jacobian.SubtractBlock2Diag(iPoint, residual.jacobian_i);

    }
  }
  END_SU2_OMP_FOR
}

void CTurbWASolver::BC_Outlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics,
                              CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {

  const bool implicit = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);

  /*--- Loop over all the vertices on this boundary marker ---*/

  SU2_OMP_FOR_STAT(OMP_MIN_SIZE)
  for (auto iVertex = 0u; iVertex < geometry->nVertex[val_marker]; iVertex++) {

    const auto iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    /*--- Check if the node belongs to the domain (i.e., not a halo node) ---*/

    if (geometry->nodes->GetDomain(iPoint)) {

      /*--- Allocate the value at the outlet ---*/

      auto V_outlet = solver_container[FLOW_SOL]->GetCharacPrimVar(val_marker, iVertex);

      /*--- Retrieve solution at the farfield boundary node ---*/

      auto V_domain = solver_container[FLOW_SOL]->GetNodes()->GetPrimitive(iPoint);

      /*--- Set various quantities in the solver class ---*/

      conv_numerics->SetPrimitive(V_domain, V_outlet);

      /*--- Set the turbulent variables. Here we use a Neumann BC such
       that the turbulent variable is copied from the interior of the
       domain to the outlet before computing the residual. ---*/

      conv_numerics->SetScalarVar(nodes->GetSolution(iPoint), nodes->GetSolution(iPoint));

      /*--- Set Normal (negate for outward convention) ---*/

      su2double Normal[MAXNDIM] = {0.0};
      for (auto iDim = 0u; iDim < nDim; iDim++)
        Normal[iDim] = -geometry->vertex[val_marker][iVertex]->GetNormal(iDim);
      conv_numerics->SetNormal(Normal);

      if (dynamic_grid)
        conv_numerics->SetGridVel(geometry->nodes->GetGridVel(iPoint),
                                  geometry->nodes->GetGridVel(iPoint));

      if (conv_numerics->GetBoundedScalar()) {
        const su2double* velocity = &V_outlet[prim_idx.Velocity()];
        const su2double density = solver_container[FLOW_SOL]->GetNodes()->GetDensity(iPoint);
        conv_numerics->SetMassFlux(BoundedScalarBCFlux(iPoint, implicit, density, velocity, Normal));
      }

      /*--- Compute the residual using an upwind scheme ---*/

      auto residual = conv_numerics->ComputeResidual(config);
      LinSysRes.AddBlock(iPoint, residual);

      /*--- Jacobian contribution for implicit integration ---*/

      if (implicit) Jacobian.AddBlock2Diag(iPoint, residual.jacobian_i);

//      /*--- Viscous contribution, commented out because serious convergence problems ---*/
//
//      su2double Coord_Reflected[MAXNDIM];
//      GeometryToolbox::PointPointReflect(nDim, geometry->nodes->GetCoord(Point_Normal),
//                                               geometry->nodes->GetCoord(iPoint), Coord_Reflected);
//      visc_numerics->SetCoord(geometry->nodes->GetCoord(iPoint), Coord_Reflected);
//      visc_numerics->SetNormal(Normal);
//
//      /*--- Conservative variables w/o reconstruction ---*/
//
//      visc_numerics->SetPrimitive(V_domain, V_outlet);
//
//      /*--- Turbulent variables w/o reconstruction, and its gradients ---*/
//
//      visc_numerics->SetScalarVar(Solution_i, Solution_j);
//      visc_numerics->SetScalarVarGradient(node[iPoint]->GetGradient(), node[iPoint]->GetGradient());
//
//      /*--- WA model switching function ---*/
//
//      visc_numerics->Setf1Switching(node[iPoint]->Getf1Switching(), node[iPoint]->Getf1Switching());
//
//      /*--- Compute residual, and Jacobians ---*/
//
//      auto residual = visc_numerics->ComputeResidual(config);
//
//      /*--- Subtract residual, and update Jacobians ---*/
//
//      LinSysRes.SubtractBlock(iPoint, residual);
//      Jacobian.SubtractBlock2Diag(iPoint, residual.jacobian_i);

    }
  }
  END_SU2_OMP_FOR
}

void CTurbWASolver::BC_Inlet_MixingPlane(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics,
                                         CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {

  const bool implicit = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);
  const auto nSpanWiseSections = config->GetnSpanWiseSections();

  /*--- Loop over all the vertices on this boundary marker ---*/
  for (auto iSpan = 0u; iSpan < nSpanWiseSections ; iSpan++){

    su2double extAverageR = solver_container[FLOW_SOL]->GetExtAverageR(val_marker, iSpan);

    /*--- Loop over all the vertices on this boundary marker ---*/

    SU2_OMP_FOR_STAT(OMP_MIN_SIZE)
    for (auto iVertex = 0u; iVertex < geometry->GetnVertexSpan(val_marker,iSpan); iVertex++) {

      /*--- find the node related to the vertex ---*/
      const auto iPoint = geometry->turbovertex[val_marker][iSpan][iVertex]->GetNode();

      /*--- using the other vertex information for retrieving some information ---*/
      const auto oldVertex = geometry->turbovertex[val_marker][iSpan][iVertex]->GetOldVertex();

      /*--- Index of the closest interior node ---*/
      const auto Point_Normal = geometry->vertex[val_marker][oldVertex]->GetNormal_Neighbor();

      /*--- Normal vector for this vertex (negate for outward convention) ---*/

      su2double Normal[MAXNDIM] = {0.0};
      for (auto iDim = 0u; iDim < nDim; iDim++)
        Normal[iDim] = -geometry->vertex[val_marker][oldVertex]->GetNormal(iDim);
      conv_numerics->SetNormal(Normal);

      /*--- Allocate the value at the inlet ---*/
      auto V_inlet = solver_container[FLOW_SOL]->GetCharacPrimVar(val_marker, oldVertex);

      /*--- Retrieve solution at the farfield boundary node ---*/

      auto V_domain = solver_container[FLOW_SOL]->GetNodes()->GetPrimitive(iPoint);

      /*--- Set various quantities in the solver class ---*/

      conv_numerics->SetPrimitive(V_domain, V_inlet);

      /*--- Set the turbulent variable states (prescribed for an inflow) ---*/

      conv_numerics->SetScalarVar(nodes->GetSolution(iPoint), &extAverageR);

      /*--- Set various other quantities in the conv_numerics class ---*/

      if (dynamic_grid)
        conv_numerics->SetGridVel(geometry->nodes->GetGridVel(iPoint),
                                  geometry->nodes->GetGridVel(iPoint));

      /*--- Compute the residual using an upwind scheme ---*/

      auto conv_residual = conv_numerics->ComputeResidual(config);

      /*--- Jacobian contribution for implicit integration ---*/

      LinSysRes.AddBlock(iPoint, conv_residual);
      if (implicit) Jacobian.AddBlock2Diag(iPoint, conv_residual.jacobian_i);

      /*--- Viscous contribution ---*/

      su2double Coord_Reflected[MAXNDIM];
      GeometryToolbox::PointPointReflect(nDim, geometry->nodes->GetCoord(Point_Normal),
                                               geometry->nodes->GetCoord(iPoint), Coord_Reflected);
      visc_numerics->SetCoord(geometry->nodes->GetCoord(iPoint), Coord_Reflected);
      visc_numerics->SetNormal(Normal);

      /*--- Conservative variables w/o reconstruction ---*/

      visc_numerics->SetPrimitive(V_domain, V_inlet);

      /*--- Turbulent variables w/o reconstruction, and its gradients ---*/

      visc_numerics->SetScalarVar(nodes->GetSolution(iPoint), &extAverageR);
      visc_numerics->SetScalarVarGradient(nodes->GetGradient(iPoint), nodes->GetGradient(iPoint));

      /*--- WA model switching function ---*/

      visc_numerics->Setf1Switching(nodes->Getf1Switching(iPoint), nodes->Getf1Switching(iPoint));

      /*--- Compute residual, and Jacobians ---*/

      auto visc_residual = visc_numerics->ComputeResidual(config);

      /*--- Subtract residual, and update Jacobians ---*/

      LinSysRes.SubtractBlock(iPoint, visc_residual);
      if (implicit) Jacobian.SubtractBlock2Diag(iPoint, visc_residual.jacobian_i);

    }
    END_SU2_OMP_FOR
  }

}

void CTurbWASolver::BC_Inlet_Turbo(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics,
                                   CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {

  const bool implicit = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);
  const auto nSpanWiseSections = config->GetnSpanWiseSections();

  CFluidModel *FluidModel = solver_container[FLOW_SOL]->GetFluidModel();

  su2double Factor_R_Inf = config->GetRFactor_FreeStream();

  /*--- Loop over all the spans on this boundary marker ---*/
  for (auto iSpan = 0; iSpan < nSpanWiseSections; iSpan++) {

    su2double rho       = solver_container[FLOW_SOL]->GetAverageDensity(val_marker, iSpan);
    su2double pressure  = solver_container[FLOW_SOL]->GetAveragePressure(val_marker, iSpan);

    FluidModel->SetTDState_Prho(pressure, rho);
    su2double muLam = FluidModel->GetLaminarViscosity();

    su2double R  = Factor_R_Inf*muLam/rho;

    /*--- Loop over all the vertices on this boundary marker ---*/

    SU2_OMP_FOR_STAT(OMP_MIN_SIZE)
    for (auto iVertex = 0u; iVertex < geometry->GetnVertexSpan(val_marker,iSpan); iVertex++) {

      /*--- find the node related to the vertex ---*/
      const auto iPoint = geometry->turbovertex[val_marker][iSpan][iVertex]->GetNode();

      /*--- using the other vertex information for retrieving some information ---*/
      const auto oldVertex = geometry->turbovertex[val_marker][iSpan][iVertex]->GetOldVertex();

      /*--- Index of the closest interior node ---*/
      const auto Point_Normal = geometry->vertex[val_marker][oldVertex]->GetNormal_Neighbor();

      /*--- Normal vector for this vertex (negate for outward convention) ---*/

      su2double Normal[MAXNDIM] = {0.0};
      for (auto iDim = 0u; iDim < nDim; iDim++)
        Normal[iDim] = -geometry->vertex[val_marker][oldVertex]->GetNormal(iDim);
      conv_numerics->SetNormal(Normal);

      /*--- Allocate the value at the inlet ---*/
      auto V_inlet = solver_container[FLOW_SOL]->GetCharacPrimVar(val_marker, oldVertex);

      /*--- Retrieve solution at the farfield boundary node ---*/

      auto V_domain = solver_container[FLOW_SOL]->GetNodes()->GetPrimitive(iPoint);

      /*--- Set various quantities in the solver class ---*/

      conv_numerics->SetPrimitive(V_domain, V_inlet);

      /*--- Set the turbulent variable states (prescribed for an inflow) ---*/

      conv_numerics->SetScalarVar(nodes->GetSolution(iPoint), &R);

      if (dynamic_grid)
        conv_numerics->SetGridVel(geometry->nodes->GetGridVel(iPoint),
                                  geometry->nodes->GetGridVel(iPoint));

      /*--- Compute the residual using an upwind scheme ---*/

      auto conv_residual = conv_numerics->ComputeResidual(config);

      /*--- Jacobian contribution for implicit integration ---*/

      LinSysRes.AddBlock(iPoint, conv_residual);
      if (implicit) Jacobian.AddBlock2Diag(iPoint, conv_residual.jacobian_i);

      /*--- Viscous contribution ---*/

      su2double Coord_Reflected[MAXNDIM];
      GeometryToolbox::PointPointReflect(nDim, geometry->nodes->GetCoord(Point_Normal),
                                               geometry->nodes->GetCoord(iPoint), Coord_Reflected);
      visc_numerics->SetCoord(geometry->nodes->GetCoord(iPoint), Coord_Reflected);
      visc_numerics->SetNormal(Normal);

      /*--- Conservative variables w/o reconstruction ---*/

      visc_numerics->SetPrimitive(V_domain, V_inlet);

      /*--- Turbulent variables w/o reconstruction, and its gradients ---*/

      visc_numerics->SetScalarVar(nodes->GetSolution(iPoint), &R);
      visc_numerics->SetScalarVarGradient(nodes->GetGradient(iPoint), nodes->GetGradient(iPoint));

      /*--- WA model switching function ---*/

      visc_numerics->Setf1Switching(nodes->Getf1Switching(iPoint), nodes->Getf1Switching(iPoint));

      /*--- Compute residual, and Jacobians ---*/

      auto visc_residual = visc_numerics->ComputeResidual(config);

      /*--- Subtract residual, and update Jacobians ---*/

      LinSysRes.SubtractBlock(iPoint, visc_residual);
      if (implicit) Jacobian.SubtractBlock2Diag(iPoint, visc_residual.jacobian_i);

    }
    END_SU2_OMP_FOR
  }

}

void CTurbWASolver::SetInletAtVertex(const su2double *val_inlet,
                                    unsigned short iMarker,
                                    unsigned long iVertex) {

  Inlet_TurbVars[iMarker][iVertex][0] = val_inlet[nDim+2+nDim];

}

su2double CTurbWASolver::GetInletAtVertex(su2double *val_inlet,
                                          unsigned long val_inlet_point,
                                          unsigned short val_kind_marker,
                                          string val_marker,
                                          const CGeometry *geometry,
                                          const CConfig *config) const {
  /*--- Local variables ---*/

  unsigned short iMarker;
  unsigned long iPoint, iVertex;
  su2double Area = 0.0;
  su2double Normal[3] = {0.0,0.0,0.0};

  /*--- Alias positions within inlet file for readability ---*/

  if (val_kind_marker == INLET_FLOW) {

    unsigned short position = nDim+2+nDim;

    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      if ((config->GetMarker_All_KindBC(iMarker) == INLET_FLOW) &&
          (config->GetMarker_All_TagBound(iMarker) == val_marker)) {

        for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++){

          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();

          if (iPoint == val_inlet_point) {

            /*-- Compute boundary face area for this vertex. ---*/

            geometry->vertex[iMarker][iVertex]->GetNormal(Normal);
            Area = GeometryToolbox::Norm(nDim, Normal);

            /*--- Access and store the inlet variables for this vertex. ---*/

            val_inlet[position] = Inlet_TurbVars[iMarker][iVertex][0];

            /*--- Exit once we find the point. ---*/

            return Area;

          }
        }
      }
    }

  }

  /*--- If we don't find a match, then the child point is not on the
   current inlet boundary marker. Return zero area so this point does
   not contribute to the restriction operator and continue. ---*/

  return Area;

}

void CTurbWASolver::SetUniformInlet(const CConfig* config, unsigned short iMarker) {
  if (config->GetMarker_All_KindBC(iMarker) == INLET_FLOW) {
    for (unsigned long iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
      Inlet_TurbVars[iMarker][iVertex][0] = GetR_Inf();
    }
  }
}

void CTurbWASolver::ComputeUnderRelaxationFactor(const CConfig *config) {

  /* Loop over the solution update given by relaxing the linear
   system for this nonlinear iteration. */

  su2double localUnderRelaxation =  1.00;
  const su2double allowableRatio =  0.99;

  SU2_OMP_FOR_STAT(omp_chunk_size)
  for (unsigned long iPoint = 0; iPoint < nPointDomain; iPoint++) {

    localUnderRelaxation = 1.0;
    su2double ratio = fabs(LinSysSol[iPoint]) / (fabs(nodes->GetSolution(iPoint, 0)) + EPS);
    /* We impose a limit on the maximum percentage that the
      turbulence variables can change over a nonlinear iteration. */
    if (ratio > allowableRatio) {
      localUnderRelaxation = min(allowableRatio / ratio, localUnderRelaxation);
    }

    /* Threshold the relaxation factor in the event that there is
     a very small value. This helps avoid catastrophic crashes due
     to non-realizable states by canceling the update. */

    if (localUnderRelaxation < 1e-10) localUnderRelaxation = 0.0;

    /* Store the under-relaxation factor for this point. */

    nodes->SetUnderRelaxation(iPoint, localUnderRelaxation);

  }
  END_SU2_OMP_FOR

}
