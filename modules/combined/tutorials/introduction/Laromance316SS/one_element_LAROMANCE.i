[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]

[Mesh]
  [cube]
    type = GeneratedMeshGenerator
    dim = 3
    xmin = 0
    xmax = 1e-3
    ymin = 0
    ymax = 1e-3
    zmin = 0
    zmax = 1e-3
    nx = 1
    ny = 1
    nz = 1
    elem_type = HEX8
  []
  [pin]
    type = ExtraNodesetGenerator
    input = cube
    new_boundary = pin
    coord = '0 0 0'
  []
[]
[AuxVariables]
  [temperature]
    initial_condition = 900.0
  []
[]

[Physics]
  [SolidMechanics]
    [QuasiStatic]
      [elasto_viscoplastic]
        # new_system = true
        add_variables = true
        strain = FINITE
        # incremental = true
        generate_output = "stress_xx stress_yy stress_zz
                           stress_xy stress_xz stress_yz
                           vonmises_stress
                           mechanical_strain_xx mechanical_strain_yy mechanical_strain_zz
                           mechanical_strain_xy mechanical_strain_xz mechanical_strain_yz"
        # additional_generate_output = 'vonmises_cauchy_stress'
        # use_automatic_differentiation = true
      []
    []
  []
[]

[Materials]
  [elasticity_tensor]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = 330e9
    poissons_ratio = 0.3
  []
  [stress]
    type = ComputeMultipleInelasticStress
    inelastic_models = rom_stress_prediction
  []
  [rom_stress_prediction]
    type = LAROMANCEStressUpdate
    temperature = temperature
    initial_cell_dislocation_density = 5.0e12
    initial_wall_dislocation_density = 6.0e11
    model = solid_mechanics:laromance/test/SS316H.json
    outputs = all
  []
[]

[Postprocessors]
  [effective_strain_avg]
    type = ElementAverageValue
    variable = effective_creep_strain
  []
  [temperature]
    type = ElementAverageValue
    variable = temperature
  []
  [cell_dislocations]
    type = ElementAverageValue
    variable = cell_dislocations
  []
  [wall_disloactions]
    type = ElementAverageValue
    variable = wall_dislocations
  []
  [vonmises_stress]
    type = ElementAverageValue
    variable = vonmises_stress
  []
[]
# [NodalKernels]
#   [force]
#     type = UserForcingFunctionNodalKernel
#     variable = disp_z
#     function = -12.6
#     boundary = 'front'
#   []
# []

[BCs]
  [fix_x]
    type = DirichletBC
    variable = disp_x
    value = 0
    boundary = 'pin'
    # preset = true
  []
  [fix_y]
    type = DirichletBC
    variable = disp_y
    value = 0
    boundary = 'pin'
    # preset = true
  []
  [fix_z]
    type = DirichletBC
    variable = disp_z
    value = 0
    boundary = 'back'
    # preset = true
  []
  [pressure_z]
    type = Pressure
    variable = disp_z
    boundary = 'front'
    # function = t
    factor = -40e6
  []
[]

[Executioner]
  type = Transient
  solve_type = NEWTON
  # petsc_options_iname = '-pc_type'
  # petsc_options_value = 'lu'
  petsc_options_iname = '-pc_type -pc_factor_shift_type'
  petsc_options_value = 'lu       nonzero              '
  snesmf_reuse_base = false
  automatic_scaling = true
  end_time = 500
  dt = 5
  nl_rel_tol = 1e-6
  # nl_abs_tol = 1e-18
  nl_forced_its = 5
  line_search = None
[]

[Outputs]
  exodus = true
  csv = true
  perf_graph =true
[]
