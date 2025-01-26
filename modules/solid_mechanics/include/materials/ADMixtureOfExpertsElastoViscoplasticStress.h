//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "Material.h"
#include "ADComputeStressBase.h"
#include "RankTwoTensor.h"
#include <torch/script.h>
#include <torch/torch.h>



class ADMixtureOfExpertsElastoViscoplasticStress : public ADComputeStressBase
{
public:
  static InputParameters validParams();

  ADMixtureOfExpertsElastoViscoplasticStress(const InputParameters & parameters);

  torch::Tensor minmax_transform(const torch::Tensor& x, double lb, double ub, double xmin, double xmax);

  torch::Tensor logtransform_transform(const torch::Tensor& x, double lb, double ub, double logxmin, double logxmax, double eps);

  torch::Tensor logtransform_inverse_transform(const torch::Tensor& z, double lb, double ub, double logxmin, double logxmax, double eps);

  torch::Tensor symlog_inverse_transform(const torch::Tensor& z, double a, double zmin, double zmax, double zbar);

  torch::Tensor minmax_forward_derivative(const torch::Tensor& x, double lb, double ub, double xmin, double xmax);

  torch::Tensor logtransform_forward_derivative(const torch::Tensor& x, double lb, double ub, double logxmin, double logxmax, double eps);

  torch::Tensor logtransform_inverse_derivative(const torch::Tensor& x, double lb, double ub, double logxmin, double logxmax, double eps);

  torch::Tensor symlog_inverse_derivative(const torch::Tensor& z, double a, double zmin, double zmax, double zbar);

  torch::Tensor eval(const torch::Tensor& tensor);

  torch::Tensor eval_jacobian(const torch::Tensor& tensor);

protected:
  torch::jit::script::Module  model;

  virtual void initQpStatefulProperties() override;

  virtual void computeQpStress() override;

  /// Name of the elasticity tensor material property
  const std::string _elasticity_tensor_name;
  /// Elasticity tensor material property
  const ADMaterialProperty<RankFourTensor> & _elasticity_tensor;
  /// Old state of the stress tensor material property
  const MaterialProperty<RankTwoTensor> & _stress_old;
  /// Old state of the mechanical strain material property
  const MaterialProperty<RankTwoTensor> & _mechanical_strain_old;

  const ADMaterialProperty<Real> & _mu;

  const ADMaterialProperty<Real> & _temp;

  ADMaterialProperty<Real> & _rhoc;

  ADMaterialProperty<Real> & _rhow;

  const MaterialProperty<Real> & _rhoc_old;

  const MaterialProperty<Real> & _rhow_old;

  const ADMaterialProperty<Real> & _flux;

  ADMaterialProperty<Real> & _total_effective_plastic_strain;

  const MaterialProperty<Real> & _total_effective_plastic_strain_old;

  const Real &_initial_rhoc;

  const Real &_initial_rhow;

  unsigned int _switch;

  ADReal effective_plastic_strain_increment;

  ADReal rhoc_rate;

  ADReal rhow_rate;

  ADMaterialProperty<RankTwoTensor> & _incremental_plastic_strain_tensor;

  const MaterialProperty<RankTwoTensor> & _incremental_plastic_strain_tensor_old;
};
