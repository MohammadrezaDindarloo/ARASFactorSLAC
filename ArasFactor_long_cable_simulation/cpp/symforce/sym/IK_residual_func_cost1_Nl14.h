// -----------------------------------------------------------------------------
// This file was autogenerated by symforce from template:
//     function/FUNCTION.h.jinja
// Do NOT modify by hand.
// -----------------------------------------------------------------------------

#pragma once

#include <Eigen/Dense>

#include <sym/rot3.h>

namespace sym {

/**
 * This function was autogenerated from a symbolic function. Do not modify by hand.
 *
 * Symbolic function: IK_residual_func_cost1_Nl14
 *
 * Args:
 *     fh1: Scalar
 *     fv1: Scalar
 *     DeltaRot: Rot3
 *     position_vector: Matrix31
 *     Rot_init: Rot3
 *     epsilon: Scalar
 *
 * Outputs:
 *     res: Matrix41
 */
template <typename Scalar>
Eigen::Matrix<Scalar, 4, 1> IkResidualFuncCost1Nl14(
    const Scalar fh1, const Scalar fv1, const sym::Rot3<Scalar>& DeltaRot,
    const Eigen::Matrix<Scalar, 3, 1>& position_vector, const sym::Rot3<Scalar>& Rot_init,
    const Scalar epsilon) {
  // Total ops: 502

  // Unused inputs
  (void)epsilon;

  // Input arrays
  const Eigen::Matrix<Scalar, 4, 1>& _DeltaRot = DeltaRot.Data();
  const Eigen::Matrix<Scalar, 4, 1>& _Rot_init = Rot_init.Data();

  // Intermediate terms (158)
  const Scalar _tmp0 = Scalar(1.0) / (fh1);
  const Scalar _tmp1 = -_DeltaRot[0] * _Rot_init[1] + _DeltaRot[1] * _Rot_init[0] +
                       _DeltaRot[2] * _Rot_init[3] + _DeltaRot[3] * _Rot_init[2];
  const Scalar _tmp2 = -2 * std::pow(_tmp1, Scalar(2));
  const Scalar _tmp3 = _DeltaRot[0] * _Rot_init[2] + _DeltaRot[1] * _Rot_init[3] -
                       _DeltaRot[2] * _Rot_init[0] + _DeltaRot[3] * _Rot_init[1];
  const Scalar _tmp4 = 1 - 2 * std::pow(_tmp3, Scalar(2));
  const Scalar _tmp5 = Scalar(0.20999999999999999) * _tmp2 + Scalar(0.20999999999999999) * _tmp4;
  const Scalar _tmp6 = _DeltaRot[0] * _Rot_init[3] - _DeltaRot[1] * _Rot_init[2] +
                       _DeltaRot[2] * _Rot_init[1] + _DeltaRot[3] * _Rot_init[0];
  const Scalar _tmp7 = 2 * _tmp6;
  const Scalar _tmp8 = _tmp1 * _tmp7;
  const Scalar _tmp9 = -_DeltaRot[0] * _Rot_init[0] - _DeltaRot[1] * _Rot_init[1] -
                       _DeltaRot[2] * _Rot_init[2] + _DeltaRot[3] * _Rot_init[3];
  const Scalar _tmp10 = 2 * _tmp3;
  const Scalar _tmp11 = _tmp10 * _tmp9;
  const Scalar _tmp12 = _tmp11 + _tmp8;
  const Scalar _tmp13 = -Scalar(0.010999999999999999) * _tmp12;
  const Scalar _tmp14 = _tmp10 * _tmp6;
  const Scalar _tmp15 = 2 * _tmp1 * _tmp9;
  const Scalar _tmp16 = Scalar(0.20999999999999999) * _tmp14 - Scalar(0.20999999999999999) * _tmp15;
  const Scalar _tmp17 = _tmp13 + _tmp16;
  const Scalar _tmp18 = _tmp17 + _tmp5;
  const Scalar _tmp19 = _tmp18 + position_vector(0, 0);
  const Scalar _tmp20 = -2 * std::pow(_tmp6, Scalar(2));
  const Scalar _tmp21 = Scalar(0.20999999999999999) * _tmp2 + Scalar(0.20999999999999999) * _tmp20 +
                        Scalar(0.20999999999999999);
  const Scalar _tmp22 = _tmp1 * _tmp10;
  const Scalar _tmp23 = _tmp7 * _tmp9;
  const Scalar _tmp24 = _tmp22 - _tmp23;
  const Scalar _tmp25 = -Scalar(0.010999999999999999) * _tmp24;
  const Scalar _tmp26 = Scalar(0.20999999999999999) * _tmp14 + Scalar(0.20999999999999999) * _tmp15;
  const Scalar _tmp27 = _tmp25 + _tmp26;
  const Scalar _tmp28 = _tmp21 + _tmp27;
  const Scalar _tmp29 = _tmp28 + position_vector(1, 0);
  const Scalar _tmp30 = std::asinh(_tmp0 * fv1);
  const Scalar _tmp31 = Scalar(1.4083112389913199) * fh1;
  const Scalar _tmp32 = -Scalar(0.20999999999999999) * _tmp11 + Scalar(0.20999999999999999) * _tmp8;
  const Scalar _tmp33 =
      -Scalar(0.010999999999999999) * _tmp20 - Scalar(0.010999999999999999) * _tmp4;
  const Scalar _tmp34 = Scalar(0.20999999999999999) * _tmp22 + Scalar(0.20999999999999999) * _tmp23;
  const Scalar _tmp35 = _tmp33 + _tmp34;
  const Scalar _tmp36 = _tmp32 + _tmp35;
  const Scalar _tmp37 = -_tmp21;
  const Scalar _tmp38 = _tmp25 - _tmp26;
  const Scalar _tmp39 = _tmp37 + _tmp38;
  const Scalar _tmp40 = _tmp39 + position_vector(1, 0);
  const Scalar _tmp41 = _tmp40 + Scalar(110.0);
  const Scalar _tmp42 = -_tmp5;
  const Scalar _tmp43 = _tmp13 - _tmp16;
  const Scalar _tmp44 = _tmp42 + _tmp43;
  const Scalar _tmp45 = _tmp44 + position_vector(0, 0);
  const Scalar _tmp46 = _tmp45 + Scalar(125.0);
  const Scalar _tmp47 = std::pow(Scalar(std::pow(_tmp41, Scalar(2)) + std::pow(_tmp46, Scalar(2))),
                                 Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp48 = _tmp46 * _tmp47;
  const Scalar _tmp49 = _tmp43 + _tmp5;
  const Scalar _tmp50 = _tmp27 + _tmp37;
  const Scalar _tmp51 = _tmp50 + position_vector(1, 0);
  const Scalar _tmp52 = _tmp51 + Scalar(110.0);
  const Scalar _tmp53 = _tmp49 + position_vector(0, 0);
  const Scalar _tmp54 = _tmp53 + Scalar(-125.0);
  const Scalar _tmp55 =
      std::sqrt(Scalar(std::pow(_tmp52, Scalar(2)) + std::pow(_tmp54, Scalar(2))));
  const Scalar _tmp56 = Scalar(1.0) / (_tmp55);
  const Scalar _tmp57 = Scalar(1.0) / (_tmp54);
  const Scalar _tmp58 = _tmp55 * _tmp57;
  const Scalar _tmp59 = _tmp58 * (_tmp49 * _tmp52 * _tmp56 - _tmp50 * _tmp54 * _tmp56);
  const Scalar _tmp60 = _tmp41 * _tmp47;
  const Scalar _tmp61 = _tmp39 * _tmp48 - _tmp44 * _tmp60 + _tmp48 * _tmp59;
  const Scalar _tmp62 = _tmp52 * _tmp57;
  const Scalar _tmp63 = Scalar(1.0) / (_tmp48 * _tmp62 - _tmp60);
  const Scalar _tmp64 = Scalar(1.0) * _tmp63;
  const Scalar _tmp65 = -_tmp32;
  const Scalar _tmp66 = _tmp33 - _tmp34;
  const Scalar _tmp67 = _tmp65 + _tmp66;
  const Scalar _tmp68 = _tmp32 + _tmp66;
  const Scalar _tmp69 = _tmp48 * _tmp68;
  const Scalar _tmp70 = _tmp60 * _tmp67 - _tmp62 * _tmp69;
  const Scalar _tmp71 = _tmp64 * _tmp70;
  const Scalar _tmp72 = Scalar(1.0) * _tmp50;
  const Scalar _tmp73 = -_tmp72;
  const Scalar _tmp74 = Scalar(1.0) / (_tmp39 + _tmp73);
  const Scalar _tmp75 = Scalar(1.0) * _tmp49;
  const Scalar _tmp76 = _tmp74 * (-_tmp44 + _tmp75);
  const Scalar _tmp77 = -_tmp48 * _tmp67 + _tmp69;
  const Scalar _tmp78 = -_tmp64 * _tmp77 + _tmp71 * _tmp76;
  const Scalar _tmp79 = _tmp35 + _tmp65;
  const Scalar _tmp80 = _tmp17 + _tmp42;
  const Scalar _tmp81 = _tmp80 + position_vector(0, 0);
  const Scalar _tmp82 = _tmp81 + Scalar(125.0);
  const Scalar _tmp83 = _tmp21 + _tmp38;
  const Scalar _tmp84 = _tmp83 + position_vector(1, 0);
  const Scalar _tmp85 = _tmp84 + Scalar(-110.0);
  const Scalar _tmp86 = std::pow(Scalar(std::pow(_tmp82, Scalar(2)) + std::pow(_tmp85, Scalar(2))),
                                 Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp87 = _tmp82 * _tmp86;
  const Scalar _tmp88 = _tmp68 * _tmp87;
  const Scalar _tmp89 = _tmp85 * _tmp86;
  const Scalar _tmp90 = _tmp62 * _tmp87 - _tmp89;
  const Scalar _tmp91 = _tmp63 * _tmp90;
  const Scalar _tmp92 = -_tmp62 * _tmp88 - _tmp70 * _tmp91 + _tmp79 * _tmp89;
  const Scalar _tmp93 = -_tmp76 * _tmp92 - _tmp77 * _tmp91 - _tmp79 * _tmp87 + _tmp88;
  const Scalar _tmp94 = Scalar(1.0) / (_tmp93);
  const Scalar _tmp95 = _tmp59 * _tmp87 - _tmp61 * _tmp91 - _tmp80 * _tmp89 + _tmp83 * _tmp87;
  const Scalar _tmp96 = _tmp94 * _tmp95;
  const Scalar _tmp97 = Scalar(1.0) / (_tmp95);
  const Scalar _tmp98 = _tmp93 * _tmp97;
  const Scalar _tmp99 = _tmp98 * (-_tmp61 * _tmp64 - _tmp78 * _tmp96);
  const Scalar _tmp100 = _tmp78 + _tmp99;
  const Scalar _tmp101 = _tmp90 * _tmp94;
  const Scalar _tmp102 = -_tmp100 * _tmp101 + Scalar(1.0);
  const Scalar _tmp103 = _tmp48 * _tmp63;
  const Scalar _tmp104 = _tmp87 * _tmp94;
  const Scalar _tmp105 = _tmp29 + Scalar(-110.0);
  const Scalar _tmp106 = _tmp19 + Scalar(-125.0);
  const Scalar _tmp107 =
      std::pow(Scalar(std::pow(_tmp105, Scalar(2)) + std::pow(_tmp106, Scalar(2))),
               Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp108 = _tmp105 * _tmp107;
  const Scalar _tmp109 = _tmp108 * fh1;
  const Scalar _tmp110 = Scalar(1.0) * _tmp97;
  const Scalar _tmp111 = _tmp64 * _tmp90 * _tmp97;
  const Scalar _tmp112 = _tmp106 * _tmp107;
  const Scalar _tmp113 = fh1 * (_tmp108 * _tmp18 - _tmp112 * _tmp28);
  const Scalar _tmp114 = _tmp62 * _tmp63;
  const Scalar _tmp115 = _tmp114 * _tmp70 + _tmp62 * _tmp68;
  const Scalar _tmp116 = _tmp114 * _tmp77 - _tmp115 * _tmp76 - _tmp68;
  const Scalar _tmp117 = _tmp98 * (_tmp114 * _tmp61 - _tmp116 * _tmp96 - _tmp59);
  const Scalar _tmp118 = _tmp116 + _tmp117;
  const Scalar _tmp119 = -_tmp101 * _tmp118 - _tmp62;
  const Scalar _tmp120 = _tmp112 * fh1;
  const Scalar _tmp121 = Scalar(333.54000000000002) - fv1;
  const Scalar _tmp122 = _tmp72 * _tmp76 + _tmp75;
  const Scalar _tmp123 = 0;
  const Scalar _tmp124 = _tmp123 * _tmp94;
  const Scalar _tmp125 = _tmp101 * _tmp123;
  const Scalar _tmp126 = -_tmp109 * _tmp58 * (_tmp100 * _tmp104 + _tmp102 * _tmp103) -
                         _tmp113 * _tmp58 * (_tmp110 * _tmp87 - _tmp111 * _tmp48) -
                         _tmp120 * _tmp58 * (_tmp103 * _tmp119 + _tmp104 * _tmp118 + Scalar(1.0)) -
                         _tmp121 * _tmp58 * (-_tmp103 * _tmp125 + _tmp124 * _tmp87);
  const Scalar _tmp127 = Scalar(1.0) / (_tmp126);
  const Scalar _tmp128 = _tmp36 * fh1;
  const Scalar _tmp129 = _tmp112 * _tmp128 + Scalar(40.024799999999999) * _tmp12 + _tmp18 * fv1;
  const Scalar _tmp130 = _tmp73 + _tmp83;
  const Scalar _tmp131 = _tmp130 * _tmp76;
  const Scalar _tmp132 = Scalar(1.0) / (-_tmp131 + _tmp75 - _tmp80);
  const Scalar _tmp133 = Scalar(1.0) * _tmp132;
  const Scalar _tmp134 = _tmp130 * _tmp74;
  const Scalar _tmp135 = _tmp92 * _tmp94;
  const Scalar _tmp136 = _tmp130 * _tmp132;
  const Scalar _tmp137 = -_tmp100 * _tmp135 + _tmp136 * _tmp99 - _tmp71;
  const Scalar _tmp138 = Scalar(1.0) * _tmp74;
  const Scalar _tmp139 = _tmp115 + _tmp117 * _tmp136 - _tmp118 * _tmp135;
  const Scalar _tmp140 = _tmp74 * (_tmp131 * _tmp133 + Scalar(1.0));
  const Scalar _tmp141 = _tmp133 * _tmp76;
  const Scalar _tmp142 = -_tmp108 * _tmp128 - Scalar(40.024799999999999) * _tmp24 - _tmp28 * fv1;
  const Scalar _tmp143 = _tmp133 * _tmp98;
  const Scalar _tmp144 = -_tmp110 * _tmp92 + _tmp130 * _tmp143;
  const Scalar _tmp145 = _tmp122 * _tmp132;
  const Scalar _tmp146 = -_tmp124 * _tmp92 - _tmp130 * _tmp145 + _tmp73;
  const Scalar _tmp147 = std::asinh(
      _tmp127 * (Scalar(1.0) * _tmp109 * (_tmp133 * _tmp99 - _tmp137 * _tmp138) +
                 Scalar(1.0) * _tmp113 * (-_tmp138 * _tmp144 + _tmp143) +
                 Scalar(1.0) * _tmp120 * (_tmp117 * _tmp133 - _tmp138 * _tmp139) +
                 Scalar(1.0) * _tmp121 * (-_tmp122 * _tmp133 - _tmp138 * _tmp146 + Scalar(1.0)) +
                 Scalar(1.0) * _tmp129 * (_tmp133 * _tmp134 - _tmp133) +
                 Scalar(1.0) * _tmp142 * (-Scalar(1.0) * _tmp140 + _tmp141)));
  const Scalar _tmp148 = Scalar(1.4083112389913199) * _tmp126;
  const Scalar _tmp149 = _tmp102 * _tmp109 * _tmp63 - _tmp111 * _tmp113 +
                         _tmp119 * _tmp120 * _tmp63 - _tmp121 * _tmp125 * _tmp63;
  const Scalar _tmp150 = Scalar(1.0) / (_tmp149);
  const Scalar _tmp151 = _tmp129 * _tmp133;
  const Scalar _tmp152 =
      std::asinh(_tmp150 * (_tmp109 * _tmp137 * _tmp74 + _tmp113 * _tmp144 * _tmp74 +
                            _tmp120 * _tmp139 * _tmp74 + _tmp121 * _tmp146 * _tmp74 -
                            _tmp134 * _tmp151 + _tmp140 * _tmp142));
  const Scalar _tmp153 = Scalar(1.4083112389913199) * _tmp149;
  const Scalar _tmp154 = _tmp100 * _tmp109 * _tmp94 + _tmp110 * _tmp113 +
                         _tmp118 * _tmp120 * _tmp94 + _tmp121 * _tmp124;
  const Scalar _tmp155 = Scalar(1.0) / (_tmp154);
  const Scalar _tmp156 = std::asinh(_tmp155 * (-_tmp109 * _tmp132 * _tmp99 - _tmp113 * _tmp143 -
                                               _tmp117 * _tmp120 * _tmp132 + _tmp121 * _tmp145 -
                                               _tmp141 * _tmp142 + _tmp151));
  const Scalar _tmp157 = Scalar(1.4083112389913199) * _tmp154;

  // Output terms (1)
  Eigen::Matrix<Scalar, 4, 1> _res;

  _res(0, 0) =
      -_tmp31 *
          (Scalar(34.083374946563197) * _tmp0 + std::cosh(Scalar(1.0) * _tmp30) -
           std::cosh(
               Scalar(0.71007031138673404) * _tmp0 *
               (-_tmp30 * _tmp31 -
                Scalar(125.0) *
                    std::sqrt(Scalar(
                        std::pow(Scalar(1 - Scalar(0.0080000000000000002) * _tmp19), Scalar(2)) +
                        Scalar(0.77439999999999998) *
                            std::pow(Scalar(1 - Scalar(0.0090909090909090905) * _tmp29),
                                     Scalar(2))))))) +
      _tmp36 + position_vector(2, 0);
  _res(1, 0) =
      -_tmp148 *
          (Scalar(34.083374946563197) * _tmp127 + std::cosh(Scalar(1.0) * _tmp147) -
           std::cosh(
               Scalar(0.71007031138673404) * _tmp127 *
               (-_tmp147 * _tmp148 -
                Scalar(125.0) *
                    std::sqrt(Scalar(
                        std::pow(Scalar(1 - Scalar(0.0080000000000000002) * _tmp53), Scalar(2)) +
                        Scalar(0.77439999999999998) *
                            std::pow(Scalar(-Scalar(0.0090909090909090905) * _tmp51 - 1),
                                     Scalar(2))))))) +
      _tmp68 + position_vector(2, 0);
  _res(2, 0) =
      -_tmp153 *
          (Scalar(34.083374946563197) * _tmp150 + std::cosh(Scalar(1.0) * _tmp152) -
           std::cosh(Scalar(0.71007031138673404) * _tmp150 *
                     (-_tmp152 * _tmp153 -
                      Scalar(125.0) *
                          std::sqrt(Scalar(
                              Scalar(0.77439999999999998) *
                                  std::pow(Scalar(-Scalar(0.0090909090909090905) * _tmp40 - 1),
                                           Scalar(2)) +
                              std::pow(Scalar(-Scalar(0.0080000000000000002) * _tmp45 - 1),
                                       Scalar(2))))))) +
      _tmp67 + position_vector(2, 0);
  _res(3, 0) =
      -_tmp157 *
          (Scalar(34.083374946563197) * _tmp155 + std::cosh(Scalar(1.0) * _tmp156) -
           std::cosh(
               Scalar(0.71007031138673404) * _tmp155 *
               (-_tmp156 * _tmp157 -
                Scalar(125.0) * std::sqrt(Scalar(
                                    Scalar(0.77439999999999998) *
                                        std::pow(Scalar(1 - Scalar(0.0090909090909090905) * _tmp84),
                                                 Scalar(2)) +
                                    std::pow(Scalar(-Scalar(0.0080000000000000002) * _tmp81 - 1),
                                             Scalar(2))))))) +
      _tmp79 + position_vector(2, 0);

  return _res;
}  // NOLINT(readability/fn_size)

// NOLINTNEXTLINE(readability/fn_size)
}  // namespace sym
